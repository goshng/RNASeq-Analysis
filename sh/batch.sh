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

SPECIES=smutans12
REPETITION=1
REPLICATE=1
function batch {
  batch-variable
  batch-output
  batch-speciesfile 

  batch-run
  batch-run-bwa-align
  batch-rmessage
  scp $BASEDIR/* cac:run/rnaseq/110111
}

function batch-variable {
  SPECIESFILE=species/$SPECIES
  # Local output directories.
  OUTPUTDIR=$ROOTANALYSISDIR/output
  BASEDIR=$OUTPUTDIR/$SPECIES
  BASERUNANALYSIS=$BASEDIR/run-analysis
  NUMBERDIR=$BASEDIR/$REPETITION
  DATADIR=$NUMBERDIR/data
  BWADIR=$NUMBERDIR/bwa
  RUNANALYSIS=$NUMBERDIR/run-analysis
  # Remote output directories.
  RINPUTDIR=/v4scratch/sc2265/rnaseq/input
  ROUTPUTDIR=/v4scratch/sc2265/rnaseq/output
  RBASEDIR=$ROUTPUTDIR/$SPECIES
  RNUMBERDIR=$ROUTPUTDIR/$SPECIES/$REPETITION
  RANALYSISDIR=$RNUMBERDIR/run-analysis
  RDATADIR=$RNUMBERDIR/data
  RBWADIR=$RNUMBERDIR/bwa
  # Compute node output directories.
  CBASEDIR=output/$SPECIES
  CNUMBERDIR=$CBASEDIR/$REPETITION
  CANALYSISDIR=$CNUMBERDIR/run-analysis
  CDATADIR=$CNUMBERDIR/data
  CBWADIR=$CNUMBERDIR/bwa
}

function batch-output {
  mkdir -p $BASEDIR/data
  mkdir -p $BASERUNANALYSIS
  mkdir -p $DATADIR
  mkdir -p $BWADIR
  mkdir -p $RUNANALYSIS
}

function batch-rmessage {
  echo "Create cac:$RDATADIR"
  echo "Create cac:$RBWADIR"
  echo "Copy $GZIPFASTAQFILE to cac:$RDATADIR"
  echo "Copy $REFGENOMEFASTA to cac:$RDATADIR"
  echo "Copy FASTQxx.fq.gz to cac:$RDATADIR"
  echo "bash run.sh at cac:run/rnaseq/110111"
}

function batch-speciesfile {
  REFGENOMEFASTA=$(grep ^REFGENOMEFASTA\: $SPECIESFILE | cut -d":" -f2)
  BWAHOWMANYNODE=$(grep ^BWAHOWMANYNODE\: $SPECIESFILE | cut -d":" -f2)
  BWAWALLTIME=$(grep ^BWAWALLTIME\: $SPECIESFILE | cut -d":" -f2)
}

function batch-run {
cat>$BASEDIR/run.sh<<EOF
#!/bin/bash
STATUS=bwa-align
bash run-\$STATUS.sh
while [ 1 ]; do
  if [ "\$STATUS" == "bwa-align" ]; then
    NJOBS=\$(qstat | grep $PROJECTNAME-BWA | wc -l) 
    if [ \$NJOBS -eq 0 ]; then
      STATUS=exit
    fi  
  fi  

  if [ "\$STATUS" == "exit" ]; then
    break
  fi

  sleep 60
done
EOF
}

function batch-run-bwa-align {

  BWA=bwa
  SAMTOOLS=samtools
  GENOMEFASTA=$(basename $REFGENOMEFASTA)

cat>$BASEDIR/run-bwa-align.sh<<EOF
#!/bin/bash
sed s/PBSARRAYSIZE/$BWAHOWMANYNODE/g < batch-bwa-align.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF

cat>$BASEDIR/batch-bwa-align.sh<<EOF
#!/bin/bash
#PBS -l walltime=${BWAWALLTIME}:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-BWA
#PBS -q ${QUEUENAME}
#PBS -m e
#PBS -M ${BATCHEMAIL}
#PBS -t 1-PBSARRAYSIZE

# The full path of the clonal origin executable.

BASEDIR=output/$SPECIES
NUMBERDIR=\$BASEDIR/$REPETITION
BWADIR=\$NUMBERDIR/bwa
DATADIR=\$NUMBERDIR/data
ANALYSISDIR=\$NUMBERDIR/run-analysis

function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  # cp \$WARG .
  # cp -r \$PBS_O_WORKDIR/bpp/Mauve-Analysis/pl .

  # All of the batchjob scripts.
  cp \$PBS_O_WORKDIR/batch-bwa-*.sh .

  # Create output directories at the compute node.
  mkdir -p \$NUMBERDIR
  cp -r $RDATADIR \$NUMBERDIR
  mkdir \$NUMBERDIR/bwa

  bwa index -p $CBWADIR/$GENOMEFASTA-bwa -a is \\
    $CDATADIR/$GENOMEFASTA
}

function retrieve-data {
  cp $CBWADIR/* $RBWADIR
}

function process-data {
  cd \$TMPDIR
  CORESPERNODE=8
  for (( i=1; i<=CORESPERNODE; i++))
  do
    g=\$((8 * (PBS_ARRAYID-1) + i + 14))
    BATCHFILENAME=batch-bwa-\$(printf "%02d" \$g)
    bash \$BATCHFILENAME&
  done
}

copy-data
process-data; wait
retrieve-data
EOF

  for g in $(eval echo {15..27}); do
    FASTQNUM=FASTQ$(printf "%02d" $g)
    GZIPFASTAQFILE=$CDATADIR/$FASTQNUM.fq.gz
    BATCHFILE=$BASEDIR/batch-bwa-$(printf "%02d" $g)
    rm -f $BATCHFILE

    # COMMAND1="$BWA aln -I -t $NUMBERCPU \
    COMMAND1="$BWA aln -t $NUMBERCPU \
              $CBWADIR/$GENOMEFASTA-bwa \
              $GZIPFASTAQFILE > $CBWADIR/$FASTQNUM.sai"
    COMMAND2="$BWA samse -n 1 \
              -f $CBWADIR/$FASTQNUM.sam \
              $CBWADIR/$GENOMEFASTA-bwa \
              $CBWADIR/$FASTQNUM.sai \
              $GZIPFASTAQFILE"
    COMMAND3="$SAMTOOLS view -bS -o $CBWADIR/$FASTQNUM.bam \
              $CBWADIR/$FASTQNUM.sam"
    COMMAND4="$SAMTOOLS sort $CBWADIR/$FASTQNUM.bam \
              $CBWADIR/$FASTQNUM.sorted"
    DELETE1="rm $CBWADIR/$FASTQNUM.sai \
                $CBWADIR/$FASTQNUM.sam \
                $CBWADIR/$FASTQNUM.bam"
    echo $COMMAND1 >> $BATCHFILE
    echo $COMMAND2 >> $BATCHFILE
    echo $COMMAND3 >> $BATCHFILE
    echo $COMMAND4 >> $BATCHFILE
    echo $DELETE1 >> $BATCHFILE
  done 
}

