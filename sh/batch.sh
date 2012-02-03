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
select SPECIES in ${SPECIESS[@]}; do 
if [ "$SPECIES" == "" ];  then
  echo -e "You need to enter something\n"
  continue
else  
  batch-variable
  batch-output
  batch-speciesfile 

  copy-data
  get-data
  ucsc-data
  batch-run
  batch-run-fastq-qc
  batch-run-bwa-align
  batch-run-de
  batch-run-parsernaseq

  batch-run-rnaz
  prepare-data-rnaz
  postprocess-rnaz
  
  #batch-run-bwa-sum
  scp -q $BASEDIR/*.sh $CAC_USERHOST:$CACWORKDIR
  scp -q $BASEDIR/job* $CAC_USERHOST:$CACWORKDIR
  scp -qr pl $CAC_USERHOST:$CACWORKDIR
  batch-rmessage
  break
fi
done
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
  echo "work at cac:$CACWORKDIR"
}

function batch-speciesfile {
  REFGENOMEID=$(grep ^REFGENOMEID\: $SPECIESFILE | cut -d":" -f2)
  REFGENOMELENGTH=$(grep ^REFGENOMELENGTH\: $SPECIESFILE | cut -d":" -f2)
  REFGENOMEFASTA=$(grep ^REFGENOMEFASTA\: $SPECIESFILE | cut -d":" -f2)
  REFGENOMEGBK=$(grep ^REFGENOMEGBK\: $SPECIESFILE | cut -d":" -f2)
  REFGENOMEGFF=$(grep ^REFGENOMEGFF\: $SPECIESFILE | cut -d":" -f2)
  RNAZNNODE=$(grep ^RNAZNNODE\: $SPECIESFILE | cut -d":" -f2)
  RNAZWALLTIME=$(grep ^RNAZWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  DENNODE=$(grep ^DENNODE\: $SPECIESFILE | cut -d":" -f2)
  DEWALLTIME=$(grep ^DEWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  PARSERNASEQNNODE=$(grep ^PARSERNASEQNNODE\: $SPECIESFILE | cut -d":" -f2)
  PARSERNASEQWALLTIME=$(grep ^PARSERNASEQWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  BWAALIGNNNODE=$(grep ^BWAALIGNNNODE\: $SPECIESFILE | cut -d":" -f2)
  BWAALIGNWALLTIME=$(grep ^BWAALIGNWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  FASTQQCNNODE=$(grep ^FASTQQCNNODE\: $SPECIESFILE | cut -d":" -f2)
  FASTQQCWALLTIME=$(grep ^FASTQQCWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  FASTQFILES=$(grep ^FASTQFILES\: $SPECIESFILE | cut -d":" -f2)
  CACWORKDIR=$(grep ^CACWORKDIR\: $SPECIESFILE | cut -d":" -f2)
  READDEPTH=$(grep ^READDEPTH\: $SPECIESFILE | cut -d":" -f2)
  REFGENOMEPTT=$(grep ^REFGENOMEPTT\: $SPECIESFILE | cut -d":" -f2)
  BWAALIGNNCPU=$(grep ^BWAALIGNNCPU\: $SPECIESFILE | cut -d":" -f2)
  # Use the following line to add more configurations
  # xxx=$(grep ^xxx\: $SPECIESFILE | cut -d":" -f2)
}

function batch-run {
cat>$BASEDIR/run.sh<<EOF
#!/bin/bash
STATUS=fastq-qc
bash run-\$STATUS.sh
while [ 1 ]; do
  sleep 60

  if [ "\$STATUS" == "fastq-qc" ]; then
    NJOBS=\$(qstat | grep $PROJECTNAME-QC | wc -l) 
    if [ \$NJOBS -eq 0 ]; then
      STATUS=bwa-align
      bash run-\$STATUS.sh
    fi  
  fi  
  sleep 30

  if [ "\$STATUS" == "bwa-align" ]; then
    NJOBS=\$(qstat | grep $PROJECTNAME-BWA | wc -l) 
    if [ \$NJOBS -eq 0 ]; then
      STATUS=de-all
      bash run-\$STATUS.sh
    fi  
  fi  
  sleep 30

  if [ "\$STATUS" == "de-all" ]; then
    NJOBS=\$(qstat | grep $PROJECTNAME-DE | wc -l) 
    if [ \$NJOBS -eq 0 ]; then
      STATUS=sum-de-all
#      STATUS=parsernaseq
      bash \$STATUS.sh
    fi  
  fi  
  sleep 30

  if [ "\$STATUS" == "parsernaseq" ]; then
    NJOBS=\$(qstat | grep $PROJECTNAME-PARSE | wc -l) 
    if [ \$NJOBS -eq 0 ]; then
      STATUS=rnaz
      bash \$STATUS.sh
    fi  
  fi  
  sleep 30

  if [ "\$STATUS" == "rnaz" ]; then
    NJOBS=\$(qstat | grep $PROJECTNAME-RNAZ | wc -l) 
    if [ \$NJOBS -eq 0 ]; then
      STATUS=sum-de-all
      bash \$STATUS.sh
    fi  
  fi  
  sleep 30

  if [ "\$STATUS" == "sum-de-all" ]; then
    STATUS=exit
  fi  
  sleep 30

  if [ "\$STATUS" == "exit" ]; then
    break
  fi
done
EOF
}

function copy-data {
cat>$BASEDIR/copy-data.sh<<EOF
#!/bin/bash
ssh -x $CAC_USERHOST mkdir -p $RBWADIR
ssh -x $CAC_USERHOST mkdir -p $RDATADIR
scp $REFGENOMEFASTA $CAC_USERHOST:$RDATADIR
REFGENOMEGBK=$REFGENOMEGBK
if [ \$REFGENOMEGBK != "NA" ]; then
  scp $REFGENOMEGFF $CAC_USERHOST:$RDATADIR
  scp $REFGENOMEPTT $CAC_USERHOST:$RDATADIR
fi
scp $BWADIR/feature-genome.out-geneonly $CAC_USERHOST:$RBWADIR
scp output/data/bacteria.fa $CAC_USERHOST:$RDATADIR
for g in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$g)
  FASTQFILE=\$(grep ^\$FASTQNUM\: $SPECIESFILE | cut -d":" -f2)
  scp \$FASTQFILE $CAC_USERHOST:$RDATADIR/\$FASTQNUM.fq.gz
done
EOF
}

function get-data {
cat>$BASEDIR/get-data.sh<<EOF
#!/bin/bash
scp $CAC_USERHOST:$RBWADIR/*.de $BWADIR
scp $CAC_USERHOST:$RBWADIR/*rrna $BWADIR
scp $CAC_USERHOST:$RBWADIR/*-sum.pos $BWADIR
scp $CAC_USERHOST:$RBWADIR/*.wig $BWADIR
scp $CAC_USERHOST:$RBWADIR/intergeniconly.maf $BWADIR
rm -f $BWADIR/rrna.tex
for g in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$g)
  perl pl/bwa-summary.pl rrnaToTex -rrna $BWADIR/\$FASTQNUM-sum.rrna >> $BWADIR/rrna.tex
done
EOF
}

function ucsc-data {
cat>$BASEDIR/send-ucsc-data.sh<<EOF
#!/bin/bash
scp -qr $BWADIR $X11_USERNAME@$X11_LOGIN:public_html/rnaseq/bwa-$SPECIES
EOF

cat>$BASEDIR/ucsc-data.sh<<EOF
#!/bin/bash
DBNAME=SmuUA159
for g in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$g)
  RNASEQNUM=rnaSeq\$(printf "%03d" \$g)
  wigEncode bwa-$SPECIES/\$FASTQNUM.wig bwa-$SPECIES/\$FASTQNUM.temp.wig bwa-$SPECIES/\$FASTQNUM.wib
  hgLoadWiggle \$DBNAME \$RNASEQNUM bwa-$SPECIES/\$FASTQNUM.temp.wig
  rm bwa-$SPECIES/\$FASTQNUM.temp.wig
  mv bwa-$SPECIES/\$FASTQNUM.wib /gbdb_cornell/\$DBNAME/wib/
  hgsql \$DBNAME -e "update \$RNASEQNUM set file='/gbdb_cornell/\$DBNAME/wib/\$FASTQNUM.wib'"
done

TRACKDBRA=trackDb.ra-$SPECIES
rm -f \$TRACKDBRA
echo "track RNAseq" > \$TRACKDBRA
echo "shortLabel RNA-seq" >> \$TRACKDBRA
echo "longLabel RNA-seq coverage" >> \$TRACKDBRA
echo "compositeTrack on" >> \$TRACKDBRA
echo "group x" >> \$TRACKDBRA
echo "visibility full" >> \$TRACKDBRA
echo "autoScale on" >> \$TRACKDBRA
echo "allButtonPair on" >> \$TRACKDBRA
echo "dragAndDrop subTracks" >> \$TRACKDBRA
echo "type wig" >> \$TRACKDBRA
echo -e "\n" >> \$TRACKDBRA
for g in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$g)
  RNASEQNUM=rnaSeq\$(printf "%03d" \$g)

  echo "  track \$RNASEQNUM" >> \$TRACKDBRA
  echo "  parent RNAseq" >> \$TRACKDBRA
  echo "  shortLabel RNA-seq \$g" >> \$TRACKDBRA
  echo "  longLabel RNA-seq sample \$g" >> \$TRACKDBRA
  echo "  priority \$g" >> \$TRACKDBRA
  echo "  type wig " >> \$TRACKDBRA
  echo "  configurable on" >> \$TRACKDBRA
  echo -e "" >> \$TRACKDBRA
done
EOF
  scp -q $BASEDIR/ucsc-data.sh $X11_USERNAME@$X11_LOGIN:public_html/rnaseq
  echo "To send the data to the genome browser:"
  echo "bash $BASEDIR/send-ucsc-data.sh"
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

function copy-data {
  cd \$TMPDIR

  # Programs and scripts.

  # All of the batchjob scripts.
  cp \$PBS_O_WORKDIR/job-fastq-qc* .

  # Create output directories at the compute node.
  mkdir -p $CBWADIR
  mkdir -p $CDATADIR
}

function retrieve-data {
  cp $CBWADIR/*cutadapt.fq.gz $RBWADIR
}

function process-data {
  cd \$TMPDIR
  CORESPERNODE=8
  FASTQFILES=( $FASTQFILES )
  for (( i=0; i<CORESPERNODE; i++))
  do
    g=\$((8 * (PBS_ARRAYID-1) + i))
    if [ \$g -lt \${#FASTQFILES[@]} ]; then
      # if [ \$g -lt 5 ]; then
        # bash job-fastq-qc \$(printf "%03d" \${FASTQFILES[\$g]})&
      # else
        bash job-fastq-qc2 \$(printf "%03d" \${FASTQFILES[\$g]})&
      # fi
    fi
  done
}

copy-data
process-data; wait
retrieve-data
cd
rm -rf \$TMPDIR
EOF

# For the first round of RNASeq
grep ^ADAPTER $SPECIESFILE | sed s/:/=/ > $BASEDIR/job-fastq-qc
cat>>$BASEDIR/job-fastq-qc<<EOF
ADAPTERSEQUENCE=ADAPTER\$1
cp $RDATADIR/FASTQ\$1.fq.gz $CBWADIR/temp.FASTQ\$1.fq.gz
/opt/epd/bin/python2.7 \$PBS_O_WORKDIR/cutadapt-1.0/cutadapt --minimum-length=25 \\
  -a \${!ADAPTERSEQUENCE} \\
  -o $CBWADIR/FASTQ\$1.cutadapt.fq.gz \\
  $CBWADIR/temp.FASTQ\$1.fq.gz
rm $CBWADIR/temp.FASTQ\$1.fq.gz
EOF

# For the 2nd round
grep ^ADAPTER $SPECIESFILE | sed s/:/=/ > $BASEDIR/job-fastq-qc2
cat>>$BASEDIR/job-fastq-qc2<<EOF
ADAPTERSEQUENCE=ADAPTER\$1
cp $RDATADIR/FASTQ\$1.fq.gz $CDATADIR
zcat $CDATADIR/FASTQ\$1.fq.gz \\
  | grep -A 3 '^@.* [^:]*:N:[^:]*:' \\
  | sed '/^--$/d' | gzip > $CBWADIR/temp.FASTQ\$1.fq.gz
rm $CDATADIR/FASTQ\$1.fq.gz
/opt/epd/bin/python2.7 \$PBS_O_WORKDIR/cutadapt-1.0/cutadapt --minimum-length=25 \\
  -a \${!ADAPTERSEQUENCE} \\
  -o $CBWADIR/FASTQ\$1.cutadapt.fq.gz \\
  $CBWADIR/temp.FASTQ\$1.fq.gz
rm $CBWADIR/temp.FASTQ\$1.fq.gz
EOF
}

function batch-run-bwa-align {
  GENOMEFASTA=$(basename $REFGENOMEFASTA)
  GENOMEGFF=$(basename $REFGENOMEGFF)
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

function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$PBS_O_WORKDIR/job-bwa-align* .
  cp \$PBS_O_WORKDIR/samtools .
  cp \$PBS_O_WORKDIR/bwa .

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CBWADIR
  cp $RDATADIR/$GENOMEFASTA $CDATADIR
  cp $RDATADIR/$GENOMEGFF $CDATADIR

  ./bwa index -p $CBWADIR/$GENOMEFASTA-bwa -a is \\
    $CDATADIR/$GENOMEFASTA
}

function retrieve-data {
  cp $CBWADIR/*.bam $RBWADIR
  cp $CBWADIR/*.pileup $RBWADIR
  cp $CBWADIR/*.wig $RBWADIR
  cp $CBWADIR/*-sum.pos $RBWADIR
  cp $CBWADIR/*-sum.rrna $RBWADIR
}

function process-data {
  cd \$TMPDIR
  CORESPERNODE=1
  FASTQFILES=( $FASTQFILES )
  for (( i=0; i<CORESPERNODE; i++))
  do
    g=\$((CORESPERNODE * (PBS_ARRAYID-1) + i))
    if [ \$g -lt \${#FASTQFILES[@]} ]; then
      # if [ \$g -lt 5 ]; then
        # bash job-bwa-align \$(printf "%03d" \${FASTQFILES[\$g]})&
      # else
        bash job-bwa-align2 \$(printf "%03d" \${FASTQFILES[\$g]})&
      # fi
    fi
  done
}

copy-data
process-data; wait
retrieve-data
cd
rm -rf \$TMPDIR
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
./samtools mpileup -q 15 -d $READDEPTH \\
  -f $CDATADIR/$GENOMEFASTA \\
  $CBWADIR/FASTQ\$1.sorted.bam \\
  > $CBWADIR/FASTQ\$1.pileup
perl pl/samtools-pileup.pl \\
  wiggle \\
  -refgenome $CDATADIR/$GENOMEFASTA \\
  -in $CBWADIR/FASTQ\$1.pileup \\
  -out $CBWADIR/FASTQ\$1.wig
./samtools view $CBWADIR/FASTQ\$1.sorted.bam \\
  | perl pl/bwa-summary.pl pos > $CBWADIR/FASTQ\$1-sum.pos
./samtools view $CBWADIR/FASTQ\$1.sorted.bam \\
  | perl pl/bwa-summary.pl rrna \\
  -gff $CDATADIR/$GENOMEGFF > $CBWADIR/FASTQ\$1-sum.rrna
EOF

# For the first round of RNASeq
cat>$BASEDIR/job-bwa-align<<EOF
GZIPFASTAQFILE=$CBWADIR/FASTQ\$1.cutadapt.fq.gz
cp $RBWADIR/FASTQ\$1.cutadapt.fq.gz $CBWADIR
./bwa aln -I -t $BWAALIGNNCPU \\
  $CBWADIR/$GENOMEFASTA-bwa \\
  \$GZIPFASTAQFILE > $CBWADIR/FASTQ\$1.sai
EOF

cat>$BASEDIR/job-bwa-align2<<EOF
GZIPFASTAQFILE=$CBWADIR/FASTQ\$1.cutadapt.fq.gz
cp $RBWADIR/FASTQ\$1.cutadapt.fq.gz $CBWADIR
./bwa aln -t $BWAALIGNNCPU \\
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

  # Find the genome name in the FASTA file.
  STR=$(head -n 1 $REFGENOMEFASTA)
  STR_ARRAY=(`echo $STR | tr " " "\n"`)
  STR_ARRAY=${STR_ARRAY[0]}
  STR_ARRAY=(`echo $STR | tr ">" "\n"`)
  GENOMENAME=${STR_ARRAY[0]}

  scp -q cac/sim/batch_task_gui.sh $CAC_USERHOST:$CACWORKDIR/batchjob.sh

cat>$BASEDIR/sum-$STATUS.sh<<EOF
#!/bin/bash
ID=\$(printf "%03d" \$1)
cut -f 4 $RBWADIR/feature-genome.out-geneonly > x.gene\$ID
paste $RBWADIR/splitdir\$ID/y* | awk '{for(i=t=0;i<NF;) t+=$++i; \$0=t}1' > x.value\$ID
paste x.gene\$ID x.value\$ID > $RBWADIR/FASTQ\$ID.de
rm x.gene\$ID x.value\$ID
EOF

cat>$BASEDIR/run-$STATUS.sh<<EOF
#!/bin/bash
ID=\$(printf "%03d" \$1)
rm -rf splitdir\$ID
mkdir splitdir\$ID
cd splitdir\$ID
split -a 5 -l 100000 -d $RBWADIR/FASTQ\$ID-sum.pos
NUMBERFILE=\$(echo \`ls|wc -l\`)
LASTFILE=\$((NUMBERFILE - 1))
cd ..
JOBIDFILE=de.jobidfile\$ID
rm -f \$JOBIDFILE
h=0
SPLITID=\$(printf "%05d" \$h)
echo "perl pl/de-count.pl join \\
  -first \\
  -shortread $CBWADIR/splitdir\$ID/x\$SPLITID \\
  -genepos $CBWADIR/feature-genome.out-geneonly \\
  -o $CBWADIR/splitdir\$ID/y\$SPLITID" >> \$JOBIDFILE
for h in \$(eval echo {1..\$LASTFILE}); do
  SPLITID=\$(printf "%05d" \$h)
  echo "perl pl/de-count.pl join \\
    -shortread $CBWADIR/splitdir\$ID/x\$SPLITID \\
    -genepos $CBWADIR/feature-genome.out-geneonly \\
    -o $CBWADIR/splitdir\$ID/y\$SPLITID" >> \$JOBIDFILE
done
REFGENOMEGBK=$REFGENOMEGBK
if [ \$REFGENOMEGBK != "NA" ]; then
  perl pl/feature-genome.pl ptt2 \\
    -geneonly \\
    -chromosome $GENOMENAME \\
    -in $RDATADIR/$GENOMEPTT \\
    -out $RBWADIR/feature-genome.out-geneonly
fi
# echo -n "How many computing nodes do you wish to use? (e.g., 3) "
# read HOW_MANY_NODE
sed s/PBSARRAYSIZE/$DENNODE/g < batch-$STATUS.sh > tbatch.sh
sed s/FASTQIDENTIFIER/\$ID/g < tbatch.sh > tbatch2.sh
nsub tbatch2.sh
rm tbatch*.sh
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

ID=FASTQIDENTIFIER
function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$PBS_O_WORKDIR/batchjob.sh .

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CBWADIR
  cp $RBWADIR/feature-genome.out-geneonly $CBWADIR
  cp -r \$PBS_O_WORKDIR/splitdir\$ID $CBWADIR
}

function retrieve-data {
  cp -r $CBWADIR/splitdir\$ID $RBWADIR
}

function process-data {
  cd \$TMPDIR
  CORESPERNODE=8
  for (( i=0; i<CORESPERNODE; i++))
  do
    bash batchjob.sh \\
      \$i \\
      \$PBS_O_WORKDIR/de.jobidfile\$ID \\
      \$PBS_O_WORKDIR/de.lockfile\$ID \\
      \$PBS_O_WORKDIR/status/\$PBS_ARRAYID \\
      PBSARRAYSIZE&
  done
}

copy-data
process-data; wait
retrieve-data
cd
rm -rf \$TMPDIR
EOF

cat>$BASEDIR/sum-de-all.sh<<EOF
for i in $FASTQFILES; do
  bash sum-de.sh \$i    
  ID=\$(printf "%03d" \$i)
  rm -rf splitdir\$ID
  rm -rf $RBWADIR/splitdir\$ID 
done 
EOF

cat>$BASEDIR/run-de-all.sh<<EOF
for i in $FASTQFILES; do
  bash run-de.sh \$i    
done 
EOF

cat>$BASEDIR/job-$STATUS<<EOF
cp $RBWADIR/FASTQ\$1-sum.pos $CBWADIR
perl pl/de-count.pl join \\
  -shortread $CBWADIR/FASTQ\$1-sum.pos \\
  -genepos $CBWADIR/feature-genome.out-geneonly \\
  -o $CBWADIR/FASTQ\$1.de
EOF

}

function batch-run-parsernaseq {
  GENOMEFASTA=$(basename $REFGENOMEFASTA)
  GENOMEGFF=$(basename $REFGENOMEGFF)
  STATUS=parsernaseq
cat>$BASEDIR/run-$STATUS.sh<<EOF
#!/bin/bash
sed s/PBSARRAYSIZE/$PARSERNASEQNNODE/g < batch-$STATUS.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF

cat>$BASEDIR/batch-$STATUS.sh<<EOF
#!/bin/bash
#PBS -l walltime=${PARSERNASEQWALLTIME}:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-PARSE
#PBS -q ${QUEUENAME}
#PBS -m e
# #PBS -M ${BATCHEMAIL}
#PBS -t 1-PBSARRAYSIZE

function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$PBS_O_WORKDIR/job-$STATUS .
  cp \$PBS_O_WORKDIR/ParseRNAseq .

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CBWADIR
  cp $RBWADIR/feature-genome.out-geneonly $CBWADIR
  cp $RDATADIR/$GENOMEFASTA $CDATADIR
  cp $RDATADIR/$GENOMEGFF $CDATADIR
}

function retrieve-data {
  cp $CBWADIR/*parsernaseq* $RBWADIR
  cp $CBWADIR/*.bed* $RBWADIR
  cp $CBWADIR/*.operon* $RBWADIR
}

function process-data {
  cd \$TMPDIR
  CORESPERNODE=1
  FASTQFILES=( $FASTQFILES )
  for (( i=0; i<CORESPERNODE; i++))
  do
    g=\$((CORESPERNODE * (PBS_ARRAYID-1) + i))
    if [ \$g -lt \${#FASTQFILES[@]} ]; then
      bash job-$STATUS \$(printf "%03d" \${FASTQFILES[\$g]})&
    fi
  done
}

copy-data
process-data; wait
retrieve-data
cd
rm -rf \$TMPDIR
EOF

cat>$BASEDIR/job-$STATUS<<EOF
FASTQNUM=FASTQ\$1
cp $RBWADIR/\$FASTQNUM.wig $CBWADIR
cp $RBWADIR/\$FASTQNUM-sum.pos $CBWADIR
# This is okay.
perl pl/transcript-parsernaseq.pl pileup \\
  -wiggle $CBWADIR/\$FASTQNUM.wig \\
  -out $CBWADIR/\$FASTQNUM.parsernaseq.pileup
# This is okay.
perl pl/transcript-parsernaseq.pl gene \\
  -feature $CBWADIR/feature-genome.out-geneonly \\
  -out $CBWADIR/feature-genome.out-parsernaseq
# This is okay.
./ParseRNAseq \\
  $CBWADIR/\$FASTQNUM.parsernaseq.pileup \\
  $CDATADIR/$GENOMEFASTA \\
  $CBWADIR/feature-genome.out-parsernaseq \\
  $CBWADIR/\$FASTQNUM.parsernaseq \\
  -c 10 -b 25 -force_gp -fmt

# This did not work.
perl pl/transcript-parsernaseq.pl bed \\
  -parsernaseq $CBWADIR/\$FASTQNUM.parsernaseq \\
  -out $CBWADIR/\$FASTQNUM.bed

# This is okay.
perl pl/transcript-parsernaseq.pl gff \\
  -parsernaseq $CBWADIR/\$FASTQNUM.parsernaseq \\
  -out $CBWADIR/\$FASTQNUM.bed2
# This is okay.
perl pl/transcript-parsernaseq.pl operon \\
  -feature $CBWADIR/feature-genome.out-geneonly \\
  -parsernaseq $CBWADIR/\$FASTQNUM.parsernaseq \\
  -out $CBWADIR/\$FASTQNUM.operon

# This works.
perl pl/bwa-pos2wig.pl end \\
  -genomeLength $REFGENOMELENGTH \\
  -in $CBWADIR/\$FASTQNUM-sum.pos \\
  -out $CBWADIR/\$FASTQNUM-end.wig

# This works.
perl pl/transcript-parsernaseq.pl adjust \\
  -end $CBWADIR/\$FASTQNUM-end.wig \\
  -operon $CBWADIR/\$FASTQNUM.operon \\
  -out $CBWADIR/\$FASTQNUM.operon2
# This works.
perl pl/transcript-parsernaseq.pl slope \\
  -windowsize 95 \\
  -end $CBWADIR/\$FASTQNUM-end.wig \\
  -operon $CBWADIR/\$FASTQNUM.operon \\
  -out $CBWADIR/\$FASTQNUM.operon.slope
EOF
}

BASEDIR=test/RNAz
NCBIBACTERIADIR=/Volumes/Elements/Documents/Projects/mauve/bacteria
MAKEBLASTDB=build/ncbi-blast-2.2.25+/bin/makeblastdb
BLASTN=build/ncbi-blast-2.2.25+/bin/blastn
BLASTP=build/ncbi-blast-2.2.25+/bin/blastp
MUSCLE=build/muscle3.8.31_i86linux64

RNAPLEX=downloads/build/RNAplex-0.2/Progs/RNAplex

INTERGENICFA=output/smutans12/1/data/intergeniconly.fa
REFGENOMEFASTA=/Volumes/Elements/Documents/Projects/mauve/bacteria/Streptococcus_mutans_UA159_uid57947/NC_004350.fna
GENEONLYFA=output/smutans12/1/data/geneonly.fa
REFPROTEINFASTA=/Volumes/Elements/Documents/Projects/mauve/bacteria/Streptococcus_mutans_UA159_uid57947/NC_004350.faa
REFGENOMEPTT=/Volumes/Elements/Documents/Projects/mauve/bacteria/Streptococcus_mutans_UA159_uid57947/NC_004350.ptt
REFGENOMEGFF=/Volumes/Elements/Documents/Projects/mauve/bacteria/Streptococcus_mutans_UA159_uid57947/NC_004350.gff

function prepare-data-rnaz {
cat>$BASEDIR/prepare-data-rnaz.sh<<EOF
mkdir output/data
cat `find $NCBIBACTERIADIR -name *.fna | grep -v NC_013928` > output/data/bacteria.fa
EOF
}





function batch-run-rnaz {
  GENOMEFASTA=$(basename $REFGENOMEFASTA)
  GENOMEPTT=$(basename $REFGENOMEPTT)
  GENOMEGFF=$(basename $REFGENOMEGFF)
  STATUS=rnaz
cat>$BASEDIR/run-$STATUS.sh<<EOF
#!/bin/bash
sed s/PBSARRAYSIZE/$RNAZNNODE/g < batch-$STATUS.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF

cat>$BASEDIR/batch-$STATUS.sh<<EOF
#!/bin/bash
#PBS -l walltime=${RNAZWALLTIME}:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-RNAZ
#PBS -q ${QUEUENAME}
#PBS -m e
# #PBS -M ${BATCHEMAIL}
#PBS -t 1-PBSARRAYSIZE

function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$PBS_O_WORKDIR/job-$STATUS* .
  cp -r \$PBS_O_WORKDIR/build .

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CBWADIR
  cp $RBWADIR/feature-genome.out-geneonly $CBWADIR
  cp $RDATADIR/$GENOMEFASTA $CDATADIR
  cp $RDATADIR/$GENOMEPTT $CDATADIR
  cp $RDATADIR/$GENOMEGFF $CDATADIR
  cp $RDATADIR/bacteria.fa $CDATADIR
}

function retrieve-data {
  cp $CBWADIR/feature-genome.out-intergeniconly $RBWADIR
  cp $CBWADIR/intergeniconly.* $RBWADIR
  cp smutans.gene2go $RBWADIR
  cp smutans.go2ngene $RBWADIR
}

function process-data {
  cd \$TMPDIR
  bash job-$STATUS &
  bash job-${STATUS}2 &
}

copy-data
process-data; wait
retrieve-data
cd
rm -rf \$TMPDIR
EOF

cat>$BASEDIR/job-$STATUS<<EOF

perl pl/feature-genome.pl ptt \\
  -intergenicregiononly \\
  -in $CDATADIR/$GENOMEPTT \\
  -out $CBWADIR/feature-genome.out-intergeniconly

perl pl/feature-genome.pl extract \\
  -bed $CBWADIR/feature-genome.out-intergeniconly \\
  -in $CDATADIR/$GENOMEFASTA \\
  -out $CBWADIR/intergeniconly.fa

# Create BLAST DB of the FASTA file.
$MAKEBLASTDB -in $CDATADIR/bacteria.fa \\
  -dbtype nucl -title bacteria -input_type fasta \\
  -out $CDATADIR/bacteria

# BLAST intergenic regions against the BLAST DB.
$BLASTN -db $CDATADIR/bacteria -query $CBWADIR/intergeniconly.fa \\
  -task blastn \\
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qseq slen sseq' \\
  -num_threads 8 \\
  -evalue 0.001 \\
  -out $CBWADIR/intergeniconly.blast

# Extract sequence alignments.
rm -rf $CBWADIR/intergeniconly 
mkdir $CBWADIR/intergeniconly
perl  pl/extract-msa.pl clust -in $CBWADIR/intergeniconly.blast \\
  -reference $REFGENOMEID \\
  -muscle $MUSCLE \\
  -outdir $CBWADIR/intergeniconly

# Align sequences using MUSCLE
rm -rf $CBWADIR/aln
mkdir $CBWADIR/aln
for f in \`ls $CBWADIR/intergeniconly\`; do
  $MUSCLE -in $CBWADIR/intergeniconly/\$f -out $CBWADIR/aln/\$f -quiet
done

# Convert FASTA-format alignments to MSA-foramt ones.
perl pl/extract-msa.pl muscle \\
  -genomeLength $REFGENOMELENGTH \\
  -indir $CBWADIR/aln \\
  -out $CBWADIR/intergeniconly.maf

EOF
cat>$BASEDIR/job-${STATUS}2<<EOF
# Enrichment Test
# Find functional categories that enrich for the predicted small RNAs.
# How? I have values for each proteins
# Useful files.
# cp /Users/goshng/Documents/Projects/mauve-analysis/emails/to/mellisa/051811/mannwhitney.R .
# cp /Users/goshng/Documents/Projects/mauve-analysis/emails/to/mellisa/051811/in.gene .
# scp swiftgen:/usr/projects/strep/SpyMGAS315/* .
#
# Download gene_association.goa_uniprot_noiea at Thu Sep 15 17:14:12 EDT 2011
# Download uniref90.fa.gz at Thu Sep 15 16:28:42 EDT 2011
# Download gene_ontology.1_2.obo at Fri Sep 16 10:31:37 EDT 2011
# 1. Prepare protein sequences of S. mutans in FASTA format (NC_004350.faa)
# 2. Prepare protein sequences of uniref90 database in FASTA format, and create BLASTDB file.
# 3. blastp S. mutans sequences against uniref90 database.
# 4. Create a file like SpyMGAS315_go_bacteria.txt using gene_association.goa_uniprot_noiea
# 5. Create another file like SpyMGAS315_go_category_names.txt
# 6. Create a file file like in.gene 
# 7. Use mannwhitney.R or a modified one to test the enrichment.
#############
#cp $REFPROTEINFASTA $BASEDIR
#
# Download at the cac's directory
# wget http://www.geneontology.org/gene-associations/submission/gene_association.goa_uniprot.gz -P .
# wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
# wget http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo -P .
cp \$PBS_O_WORKDIR/gene_association.goa_uniprot.gz .
cp \$PBS_O_WORKDIR/uniref90.fasta.gz .
cp \$PBS_O_WORKDIR/gene_ontology.1_2.obo .
gunzip gene_association.goa_uniprot.gz 
gunzip uniref90.fasta.gz
$MAKEBLASTDB -in uniref90.fasta \\
  -dbtype prot -title uniref90 -input_type fasta \\
  -out uniref90
  
$BLASTP -db uniref90 -query $CDATADIR/$GENOMEFASTA \\
  -task blastp \\
  -outfmt 6 \\
  -num_threads 8 \\
  -evalue 0.001 \\
  -out geneonly.blast
# A machine at the cluster was used.
perl pl/geneontology.pl gene2go \\
  -gff $CDATADIR/$GENOMEGFF \\
  -blast geneonly.blast \\
  -goa gene_association.goa_uniprot \\
  -out smutans.gene2go

perl pl/geneontology.pl go2ngene \\
  -gene2go smutans.gene2go \\
  -obo gene_ontology.1_2.obo \\
  -out smutans.go2ngene
# Rscript R/mannwhitney.R > target_go.txt

EOF
}

function postprocess-rnaz {
  CBWARNAZDIR=$CBWADIR/rnaz
  TARGETRNA=$CDATADIR/startcodon.fa
cat>$BASEDIR/postprocess-rnaz.sh<<EOF
mkdir $CBWARNAZDIR
# Apply RNAz to the mulitple sequence alignments.
# From Stefan:
# I guess the input to rnazCluster.pl was not sorted correctly.
# Please sort the RNAz results first by chrom and start position and try again.
rnazWindow.pl --max-seq=15 $CBWADIR/intergeniconly.maf > $CBWARNAZDIR/intergeniconly.rnaz
# RNAz --show-gaps --both-strands --predict-strand --cutoff=0.5 $CBWARNAZDIR/intergeniconly.rnaz > $CBWARNAZDIR/rnaz.out
RNAz --show-gaps --cutoff=0.5 $CBWARNAZDIR/intergeniconly.rnaz > $CBWARNAZDIR/rnaz.out
rm -rf results
rnazOutputSort.pl $CBWARNAZDIR/rnaz.out > $CBWARNAZDIR/rnaz.out.sorted
rnazCluster.pl --html $CBWARNAZDIR/rnaz.out.sorted > $CBWARNAZDIR/results.dat
rnazAnnotate.pl --bed $CBWARNAZDIR/rfam.bed $CBWARNAZDIR/results.dat > $CBWARNAZDIR/annotated.dat
rnazIndex.pl --bed $CBWARNAZDIR/annotated.dat > $CBWARNAZDIR/results.bed
rnazIndex.pl --html $CBWARNAZDIR/annotated.dat > results/index.html
grep "^locus" $CBWARNAZDIR/annotated.dat | cut -f 9 | sort > $CBWARNAZDIR/results.known
open results/index.html

# Apply RNAplex
cp $REFGENOMEPTT $CDATADIR
cp $REFGENOMEFASTA $CDATADIR
perl pl/feature-genome.pl ptt2 \\
  -startcodon \\
  -in $CDATADIR/$GENOMEPTT \\
  -out $CDATADIR/feature-genome.out-startcodon
perl pl/feature-genome.pl extract \\
  -bed $CDATADIR/feature-genome.out-startcodon \\
  -in $CDATADIR/$GENOMEFASTA \\
  -out $CDATADIR/startcodon.fa
perl pl/feature-genome.pl extract -bed $CBWARNAZDIR/results.bed \\
     -in $CDATADIR/$GENOMEFASTA \\
     -out $CBWARNAZDIR/results.fa
# I need to create $TARGETRNA
# This take time ...
# NOTE: What values should I use?
RNAplfold -W 240 -L 160 -u 30 -O < $TARGETRNA > $CBWARNAZDIR/target.rnaplfold
# RNAplfold -W 40 -L 160 -u 30 -O < $CBWARNAZDIR/results.fa
RNAplfold -W 40 -L 30 -u 30 -O < $CBWARNAZDIR/results.fa
mkdir $CBWARNAZDIR/profiles
mv *.ps $CBWARNAZDIR/profiles
mv *openen $CBWARNAZDIR/profiles

#$RNAPLEX -t $TARGETRNA -q $CBWARNAZDIR/results.fa > $CBWARNAZDIR/target.rnaplex
$RNAPLEX -t $TARGETRNA -q $CBWARNAZDIR/results.fa -l 25 -a $CBWARNAZDIR/profiles > $CBWARNAZDIR/target.rnaplex.acc
perl pl/rnaplex.pl rank -in $CBWARNAZDIR/target.rnaplex.acc -out $CBWARNAZDIR/target.rnaplex.acc.out
# This R script is only for output/cornell/1/bwa/rnaz
Rscript R/mannwhitney.R > $CBWARNAZDIR/target_go.txt
echo "See $CBWARNAZDIR/target_go.txt"

EOF
  echo "bash $BASEDIR/postprocess-rnaz.sh to analyze RNAz data"
}
