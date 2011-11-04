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

SPECIES=cornell
REPETITION=1
REPLICATE=1
function batch {
  batch-variable
  batch-output
  batch-speciesfile 

  copy-data
  batch-run
  batch-run-fastq-qc
  batch-run-bwa-align
  batch-run-de
  #batch-run-bwa-sum
  scp -q $BASEDIR/* $CAC_USERHOST:$CACWORKDIR
  scp -qr pl $CAC_USERHOST:$CACWORKDIR
  batch-rmessage
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
  echo "bash $BASEDIR/copy-data.sh"
}

function batch-speciesfile {
  REFGENOMEFASTA=$(grep ^REFGENOMEFASTA\: $SPECIESFILE | cut -d":" -f2)
  BWAALIGNNNODE=$(grep ^BWAALIGNNNODE\: $SPECIESFILE | cut -d":" -f2)
  BWAALIGNWALLTIME=$(grep ^BWAALIGNWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  FASTQQCNNODE=$(grep ^FASTQQCNNODE\: $SPECIESFILE | cut -d":" -f2)
  FASTQQCWALLTIME=$(grep ^FASTQQCWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  FASTQFILES=$(grep ^FASTQFILES\: $SPECIESFILE | cut -d":" -f2)
  CACWORKDIR=$(grep ^CACWORKDIR\: $SPECIESFILE | cut -d":" -f2)
  READDEPTH=$(grep ^READDEPTH\: $SPECIESFILE | cut -d":" -f2)
  REFGENOMEPTT=$(grep ^REFGENOMEPTT\: $SPECIESFILE | cut -d":" -f2)
  # Use the following line to add more configurations
  # xxx=$(grep ^xxx\: $SPECIESFILE | cut -d":" -f2)
}

function batch-run {
cat>$BASEDIR/run.sh<<EOF
#!/bin/bash
STATUS=fastq-qc
bash run-\$STATUS.sh
while [ 1 ]; do

  if [ "\$STATUS" == "fastq-qc" ]; then
    NJOBS=\$(qstat | grep $PROJECTNAME-QC | wc -l) 
    if [ \$NJOBS -eq 0 ]; then
      STATUS=bwa-align
    fi  
  fi  

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

function copy-data {
cat>$BASEDIR/copy-data.sh<<EOF
#!/bin/bash
ssh -x $CAC_USERHOST mkdir -p $RBWADIR
ssh -x $CAC_USERHOST mkdir -p $RDATADIR
scp $REFGENOMEFASTA $CAC_USERHOST:$RDATADIR
scp $REFGENOMEPTT $CAC_USERHOST:$RDATADIR
for g in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$g)
  FASTQFILE=\$(grep ^\$FASTQNUM\: $SPECIESFILE | cut -d":" -f2)
  scp \$FASTQFILE $CAC_USERHOST:$RDATADIR/\$FASTQNUM.fq.gz
done
EOF
}
  
function batch-run-fastq-qc {
  
cat>$BASEDIR/run-fastq-qc.sh<<EOF
#!/bin/bash
sed s/PBSARRAYSIZE/$FASTQQCNNODE/g < batch-fastq-qc.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF

cat>$BASEDIR/batch-fastq-qc.sh<<EOF
#!/bin/bash
#PBS -l walltime=${FASTQQCWALLTIME}:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-QC
#PBS -q ${QUEUENAME}
#PBS -m e
#PBS -M ${BATCHEMAIL}
#PBS -t 1-PBSARRAYSIZE

BASEDIR=output/$SPECIES
NUMBERDIR=\$BASEDIR/$REPETITION
BWADIR=\$NUMBERDIR/bwa
DATADIR=\$NUMBERDIR/data
ANALYSISDIR=\$NUMBERDIR/run-analysis

function copy-data {
  cd \$TMPDIR

  # Programs and scripts.

  # All of the batchjob scripts.
  cp \$PBS_O_WORKDIR/job-fastq-qc* .

  # Create output directories at the compute node.
  mkdir -p $CBWADIR
  cp -r $RDATADIR $CNUMBERDIR
}

function retrieve-data {
  cp $CBWADIR/* $RBWADIR
}

function process-data {
  cd \$TMPDIR
  CORESPERNODE=8
  FASTQFILES=( $FASTQFILES )
  for (( i=0; i<CORESPERNODE; i++))
  do
    g=\$((8 * (PBS_ARRAYID-1) + i))
    if [ \$g -lt 5 ]; then
      bash job-fastq-qc \$(printf "%03d" \${FASTQFILES[\$g]})&
    else
      bash job-fastq-qc2 \$(printf "%03d" \${FASTQFILES[\$g]})&
    fi
  done
}

copy-data
process-data; wait
retrieve-data
EOF

# For the first round of RNASeq
grep ^ADAPTER $SPECIESFILE | sed s/:/=/ > $BASEDIR/job-fastq-qc
cat>>$BASEDIR/job-fastq-qc<<EOF
ADAPTERSEQUENCE=ADAPTER\$1
mv $CDATADIR/FASTQ\$1.fq.gz $CBWADIR/temp.FASTQ\$1.fq.gz
/opt/epd/bin/python2.7 \$PBS_O_WORKDIR/cutadapt-1.0/cutadapt --minimum-length=25 \\
  -a \${!ADAPTERSEQUENCE} \\
  -o $CBWADIR/FASTQ\$1.cutadapt.fq.gz \\
  $CBWADIR/temp.FASTQ\$1.fq.gz
EOF

# For the 2nd round
grep ^ADAPTER $SPECIESFILE | sed s/:/=/ > $BASEDIR/job-fastq-qc2
cat>>$BASEDIR/job-fastq-qc2<<EOF
ADAPTERSEQUENCE=ADAPTER\$1
zcat $CDATADIR/FASTQ\$1.fq.gz \\
  | grep -A 3 '^@.* [^:]*:N:[^:]*:' \\
  | sed '/^--$/d' | gzip > $CBWADIR/temp.FASTQ\$1.fq.gz
/opt/epd/bin/python2.7 \$PBS_O_WORKDIR/cutadapt-1.0/cutadapt --minimum-length=25 \\
  -a \${!ADAPTERSEQUENCE} \\
  -o $CBWADIR/FASTQ\$1.cutadapt.fq.gz \\
  $CBWADIR/temp.FASTQ\$1.fq.gz
EOF
}

function batch-run-bwa-align {
  GENOMEFASTA=$(basename $REFGENOMEFASTA)
cat>$BASEDIR/run-bwa-align.sh<<EOF
#!/bin/bash
sed s/PBSARRAYSIZE/$BWAALIGNNNODE/g < batch-bwa-align.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF

cat>$BASEDIR/batch-bwa-align.sh<<EOF
#!/bin/bash
#PBS -l walltime=${BWAALIGNWALLTIME}:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-BWA
#PBS -q ${QUEUENAME}
#PBS -m e
# #PBS -M ${BATCHEMAIL}
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
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$PBS_O_WORKDIR/job-bwa-align* .
  cp \$PBS_O_WORKDIR/samtools .
  cp \$PBS_O_WORKDIR/bwa .

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  cp -r $RBWADIR $CNUMBERDIR
  cp $RDATADIR/$GENOMEFASTA $CDATADIR

  bwa index -p $CBWADIR/$GENOMEFASTA-bwa -a is \\
    $CDATADIR/$GENOMEFASTA
}

function retrieve-data {
  cp $CBWADIR/*.pileup $RBWADIR
  cp $CBWADIR/*.wig $RBWADIR
  cp $CBWADIR/*.sum-pos $RBWADIR
}

function process-data {
  cd \$TMPDIR
  CORESPERNODE=8
  FASTQFILES=( $FASTQFILES )
  for (( i=0; i<CORESPERNODE; i++))
  do
    g=\$((8 * (PBS_ARRAYID-1) + i))
    if [ \$g -lt 5 ]; then
      bash job-bwa-align \$(printf "%03d" \${FASTQFILES[\$g]})&
    else
      bash job-bwa-align2 \$(printf "%03d" \${FASTQFILES[\$g]})&
    fi
  done
}

copy-data
process-data; wait
retrieve-data
EOF

cat>$BASEDIR/job-bwa-align-common<<EOF
./bwa samse -n 1 \\
  -f $CBWADIR/FASTQ\$1.sam \\
  $CBWADIR/$GENOMEFASTA-bwa \\
  $CBWADIR/FASTQ\$1.sai \\
  \$GZIPFASTAQFILE
./samtools view -bS -o $CBWADIR/FASTQ\$1.bam \\
  $CBWADIR/FASTQ\$1.sam
./samtools sort $CBWADIR/FASTQ\$1.bam \\
  $CBWADIR/FASTQ\$1.sorted
rm $CBWADIR/FASTQ\$1.sai \\
  $CBWADIR/FASTQ\$1.sam \\
  $CBWADIR/FASTQ\$1.bam
./samtools mpileup -q 15 -d $READDEPTH \\
  -f $CDATADIR/$GENOMEFASTA \\
  $CBWADIR/FASTQ\$1.sorted.bam \\
  > $CBWADIR/FASTQ\$1.pileup
perl pl/samtools-pileup.pl \\
  wiggle \\
  -refgenome $CDATADIR/$GENOMEFASTA \\
  -in $CBWADIR/FASTQ\$1.pileup \\
  -out $CBWADIR/FASTQ\$1.wig
./samtool view $CBWADIR/FASTQ\$1.sorted.bam \\
  | perl pl/bwa-summary.pl pos > $CBWADIR/FASTQ\$1-sum.pos
EOF

# For the first round of RNASeq
cat>$BASEDIR/job-bwa-align<<EOF
GZIPFASTAQFILE=$CBWADIR/FASTQ\$1.cutadapt.fq.gz
./bwa aln -I -t $NUMBERCPU \\
  $CBWADIR/$GENOMEFASTA-bwa \\
  \$GZIPFASTAQFILE > $CBWADIR/FASTQ\$1.sai
EOF
cat>$BASEDIR/job-bwa-align2<<EOF
GZIPFASTAQFILE=$CBWADIR/FASTQ\$1.cutadapt.fq.gz
./bwa aln -t $NUMBERCPU \\
  $CBWADIR/$GENOMEFASTA-bwa \\
  \$GZIPFASTAQFILE > $CBWADIR/FASTQ\$1.sai
EOF
  cat $BASEDIR/job-bwa-align-common >> $BASEDIR/job-bwa-align
  cat $BASEDIR/job-bwa-align-common >> $BASEDIR/job-bwa-align2
}

function batch-run-de {
  STATUS=de
  GENOMEFASTA=$(basename $REFGENOMEFASTA)
  GENOMEPTT=$(basename $REFGENOMEPTT)
cat>$BASEDIR/run-$STATUS.sh<<EOF
#!/bin/bash
sed s/PBSARRAYSIZE/$DENNODE/g < batch-$STATUS.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF

cat>$BASEDIR/batch-$STATUS.sh<<EOF
#!/bin/bash
#PBS -l walltime=${DEWALLTIME}:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-DE
#PBS -q ${QUEUENAME}
#PBS -m e
# #PBS -M ${BATCHEMAIL}
#PBS -t 1-PBSARRAYSIZE

function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$PBS_O_WORKDIR/job-bwa-align* .
  cp \$PBS_O_WORKDIR/samtools .
  cp \$PBS_O_WORKDIR/bwa .

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  cp -r $RBWADIR $CNUMBERDIR
  cp $RDATADIR/$GENOMEFASTA $CDATADIR
  cp $RDATADIR/$GENOMEPTT $CDATADIR
}

function retrieve-data {
  cp $CBWADIR/feature-genome.* $RBWADIR
  cp $CBWADIR/*.de $RBWADIR
}

function process-data {
  cd \$TMPDIR
  CORESPERNODE=8
  FASTQFILES=( $FASTQFILES )
  for (( i=0; i<CORESPERNODE; i++))
  do
    g=\$((8 * (PBS_ARRAYID-1) + i))
    bash job-$STATUS \$(printf "%03d" \${FASTQFILES[\$g]})&
  done
}

copy-data
process-data; wait
retrieve-data
EOF

cat>$BASEDIR/job-$STATUS<<EOF
perl pl/feature-genome.pl ptt2 \\
  -geneonly \\
  -in $CDATADIR/$GENOMEPTT \\
  -out $CBWADIR/feature-genome.out-geneonly
perl pl/de-count.pl join \\
  -shortread $CBWADIR/FASTQ\$1-sum.pos \\
  -genepos $CBWADIR/feature-genome.out-geneonly \\
  -o $CBWADIR/FASTQ\$1.de
EOF

}

