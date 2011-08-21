###############################################################################
# Copyright (C) 2011 Sang Chul Choi
#
# This file is part of RNASeq Analysis.
# 
# Mauve Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Mauve Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Mauve Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

function job-truncate-reads {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SPECIESS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      HOWMANYREPETITION=95
      #HOWMANYREPETITION=8
      BATCHFILE=batch.sh
      BWA=./bwa
      SAMTOOLS=./samtools

      echo "#!/bin/bash" > $BATCHFILE

      global-variable $SPECIES 1
      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      #for g in $(eval echo {1..$NUMFASTQFILE}); do
        #FASTQNUM=FASTQ$(printf "%02d" $g)
        #GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
#
        #COMMAND1="gunzip $GZIPFASTAQFILE" 
        #echo $COMMAND1 >> $BATCHFILE
      #done

      for REPETITION in $(eval echo {1..$HOWMANYREPETITION}); do
        BATCHJOBFILE=batch$(printf "%02d" $REPETITION)
        echo "nsub $BATCHJOBFILE" >> $BATCHFILE

        ROOTANALYSISDIR=/v4scratch/sc2265/rnaseq
        global-variable $SPECIES $REPETITION
        read-species
        
        CUTSIZE=$REPETITION
        GENOMEFASTA=$(basename $REFGENOMEFASTA)
        REFGENOMELENGTH=$(grep REFGENOMELENGTH $SPECIESFILE | cut -d":" -f2)
        READDEPTH=$(grep READDEPTH $SPECIESFILE | cut -d":" -f2)
cat>$BATCHJOBFILE<<EOF
#!/bin/bash
#PBS -l walltime=59:00:00,nodes=1
#PBS -A acs4_0001
#PBS -j oe
#PBS -N rnaseq
#PBS -q v4
#PBS -m e
#PBS -M schoi@cornell.edu

# Create input and output directory
OUTPUTDIR=\$TMPDIR/output
INPUTDIR=\$TMPDIR/input
mkdir \$INPUTDIR
mkdir \$OUTPUTDIR

DATADIR=/v4scratch/sc2265/rnaseq/data
PLDIR=/v4scratch/sc2265/rnaseq/pl
# copy data
cp \$DATADIR/NC_004350.fna \$INPUTDIR/
cp -r \$PLDIR \$TMPDIR
cp /v4scratch/sc2265/rnaseq/samtools-0.1.16/samtools \$TMPDIR
cp /v4scratch/sc2265/rnaseq/bwa-0.5.9/bwa \$TMPDIR
cp \$PBS_O_WORKDIR/batch* \$TMPDIR

cd \$TMPDIR

# Create a number directory
mkdir -p $BASERUNANALYSIS
mkdir -p $NUMBERDIR
mkdir -p $DATADIR
mkdir -p $BWADIR
mkdir -p $MAQDIR
mkdir -p $BOWTIEDIR
# Copy the genome file
cp $REFGENOMEFASTA $DATADIR
# BWA index the genome
./bwa index -p $BWADIR/$GENOMEFASTA-bwa -a is $DATADIR/$GENOMEFASTA

for g in \$(eval echo {1..$NUMFASTQFILE}); do
  BATCHJOBFILENUM=$BATCHJOBFILE\$(printf "%02d" \$g)
  bash \$BATCHJOBFILENUM&
done

wait

rm -f $DATADIR/$GENOMEFASTA* $BWADIR/$GENOMEFASTA-bwa.*

cp -r output \$PBS_O_WORKDIR

cd
rm -rf \$TMPDIR
EOF

        for g in $(eval echo {1..$NUMFASTQFILE}); do
          ############################################
          # Shorten the reads
          FASTQNUM=FASTQ$(printf "%02d" $g)
          BATCHJOBFILENUM=$BATCHJOBFILE$(printf "%02d" $g)

          GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
          FASTAQFILE=${GZIPFASTAQFILE%.*}
          # GZIPFASTAQDIR=$DATADIR
          OUTFILE=$DATADIR/$FASTQNUM.cut-$CUTSIZE

            #--fastq $FASTAQFILE \
          COMMAND1="perl pl/fastq-sample.pl cut \
            --fastq $GZIPFASTAQFILE \
            --out $OUTFILE \
            --cutsize $CUTSIZE"
          echo $COMMAND1 > $BATCHJOBFILENUM

          echo "gzip $OUTFILE" >> $BATCHJOBFILENUM

          ############################################
          # BWA align
          GZIPFASTAQFILE=$DATADIR/$FASTQNUM.cut-$CUTSIZE.gz
          COMMAND1="$BWA aln -I -t $NUMBERCPU \
                    $BWADIR/$GENOMEFASTA-bwa \
                    $GZIPFASTAQFILE > $BWADIR/$FASTQNUM.sai"
          COMMAND2="$BWA samse -n 1 \
                    -f $BWADIR/$FASTQNUM.sam \
                    $BWADIR/$GENOMEFASTA-bwa \
                    $BWADIR/$FASTQNUM.sai \
                    $GZIPFASTAQFILE"
          COMMAND3="$SAMTOOLS view -bS -o $BWADIR/$FASTQNUM.bam \
                    $BWADIR/$FASTQNUM.sam"

          COMMAND4="$SAMTOOLS sort $BWADIR/$FASTQNUM.bam \
                    $BWADIR/$FASTQNUM.sorted"

          echo $COMMAND1 >> $BATCHJOBFILENUM
          echo $COMMAND2 >> $BATCHJOBFILENUM
          echo $COMMAND3 >> $BATCHJOBFILENUM
          echo $COMMAND4 >> $BATCHJOBFILENUM


          ############################################
          # Pileup
          COMMAND1="$SAMTOOLS mpileup -C50 -d $READDEPTH \
                    -f $DATADIR/$GENOMEFASTA \
                    $BWADIR/$FASTQNUM.sorted.bam \
                    > $BWADIR/$FASTQNUM.pileup"
          COMMAND2="perl pl/samtools-pileup.pl \
                    wiggle \
                    -refgenome $DATADIR/$GENOMEFASTA \
                    -in $BWADIR/$FASTQNUM.pileup \
                    -out $BWADIR/$FASTQNUM.wig" 

          echo $COMMAND1 >> $BATCHJOBFILENUM
          echo $COMMAND2 >> $BATCHJOBFILENUM

          ############################################
          # Cleaning up.
          DELETE1="rm $BWADIR/$FASTQNUM.sai \
                      $BWADIR/$FASTQNUM.sam \
                      $BWADIR/$FASTQNUM.bam \
                      $BWADIR/$FASTQNUM.sorted.bam \
                      $BWADIR/$FASTQNUM.pileup \
                      $DATADIR/$FASTQNUM.cut-$CUTSIZE.gz"
          echo $DELETE1 >> $BATCHJOBFILENUM
        done
                      
      done
      break
    fi
  done
}
