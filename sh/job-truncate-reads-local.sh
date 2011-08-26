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

function job-truncate-reads-local {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SPECIESS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      HOWMANYREPETITION=99
      #HOWMANYREPETITION=8
      BATCHFILE=batch.sh
      # BWA=./bwa
      # SAMTOOLS=./samtools

      echo "#!/bin/bash" > $BATCHFILE

      global-variable $SPECIES 1
      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)

      for REPETITION in $(eval echo {1..$HOWMANYREPETITION}); do
        BATCHJOBFILE=batch$(printf "%02d" $REPETITION)
        echo "bash $BATCHJOBFILE" >> $BATCHFILE

        # For CAC cluster
        global-variable $SPECIES $REPETITION
        read-species
        
        SIZE=$REPETITION
        GENOMEFASTA=$(basename $REFGENOMEFASTA)
        REFGENOMELENGTH=$(grep REFGENOMELENGTH $SPECIESFILE | cut -d":" -f2)
        READDEPTH=$(grep READDEPTH $SPECIESFILE | cut -d":" -f2)

cat>$BATCHJOBFILE<<EOF
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
$BWA index -p $BWADIR/$GENOMEFASTA-bwa -a is $DATADIR/$GENOMEFASTA

g=1
BATCHJOBFILENUM=$BATCHJOBFILE\$(printf "%02d" \$g)
bash \$BATCHJOBFILENUM

rm -f $DATADIR/$GENOMEFASTA* $BWADIR/$GENOMEFASTA-bwa.*

EOF

        for g in $(eval echo {1..$NUMFASTQFILE}); do
          ############################################
          # Shorten the reads
          FASTQNUM=FASTQ$(printf "%02d" $g)
          BATCHJOBFILENUM=$BATCHJOBFILE$(printf "%02d" $g)

          GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
          FASTAQFILE=${GZIPFASTAQFILE%.*}
          # GZIPFASTAQDIR=$DATADIR
          OUTFILE=$DATADIR/$FASTQNUM.sample-$SIZE

          COMMAND1="perl pl/fastq-sample.pl sample \
            --fastq $GZIPFASTAQFILE \
            --out $OUTFILE \
            --samplesize $SIZE"
          echo $COMMAND1 > $BATCHJOBFILENUM
          echo "gzip $OUTFILE" >> $BATCHJOBFILENUM

          ############################################
          # BWA align
          GZIPFASTAQFILE=$OUTFILE.gz
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
          COMMAND1="$SAMTOOLS mpileup -q 20 -d $READDEPTH \
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
                      $OUTFILE.gz"
          echo $DELETE1 >> $BATCHJOBFILENUM
        done
                      
      done
      break
    fi
  done
}
