###############################################################################
# Copyright (C) 2012 Sang Chul Choi
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

REPETITION=1
function batch2 {
  select SPECIES in ${SPECIESS[@]}; do 
  if [ "$SPECIES" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  else  
    batch2-variable
    batch2-output
    batch2-speciesfile 
    batch2-push-data
    batch2-get-data
    batch2-make-job

    ucsc-data
    batch2-run-fastqc
    batch2-run-cram
    batch2-run-cram2fastq
    batch2-run-bwa-align
    batch2-run-bias 
    batch2-run-coverage
    batch2-run-samtools-pileup

    # This would replace fastqc, bwa-align, and de.
    batch2-run-qcalignde
    # This would create a tux input file
    batch2-run-tux

    create-index

    # batch2-run-blast
    batch-run-parsernaseq
  #
  #  batch-run-rnaz
  #  prepare-data-rnaz
  #  postprocess-rnaz
  #  
  #  #batch-run-bwa-sum
    batch2-copy-scripts
    batch2-rmessage
    break
  fi
  done
}

function batch2-rmessage {
  echo "bash $BASEDIR/push-data.sh"
  echo "work at cac:$CACWORKDIR"
  echo "bash $BASEDIR/get-data.sh"
}

function batch2-copy-scripts {
  ssh -x $CAC_USERHOST mkdir -p $CACWORKDIR
  scp -q $BASEDIR/*.sh $CAC_USERHOST:$CACWORKDIR
  scp -q $BASEDIR/*.R $CAC_USERHOST:$CACWORKDIR
  scp -q $BASEDIR/job* $CAC_USERHOST:$CACWORKDIR
  scp -qr pl $CAC_USERHOST:$CACWORKDIR
}

# Set the shell variables for batch2 SHELL script
function batch2-variable {
  SPECIESFILE=species/$SPECIES

  # ROOTANALYSISDIR is set by run.sh.
  OUTPUTDIR=$ROOTANALYSISDIR/output

  # Local output directories.
  BASEDIR=$OUTPUTDIR/$SPECIES
  NUMBERDIR=$BASEDIR/$REPETITION
  DATADIR=$NUMBERDIR/data
  BWADIR=$NUMBERDIR/bwa
  SUBREADDIR=$NUMBERDIR/subread
  RUNANALYSIS=$NUMBERDIR/run-analysis

  # Remote output directories.
  # ROUTPUTDIR is set by conf.sh.
  RBASEDIR=$ROUTPUTDIR/$SPECIES
  RNUMBERDIR=$ROUTPUTDIR/$SPECIES/$REPETITION
  RDATADIR=$RNUMBERDIR/data
  RBWADIR=$RNUMBERDIR/bwa
  RSUBREADDIR=$RNUMBERDIR/subread
  RANALYSISDIR=$RNUMBERDIR/run-analysis

  # Compute node output directories.
  CBASEDIR=output/$SPECIES
  CNUMBERDIR=$CBASEDIR/$REPETITION
  CDATADIR=$CNUMBERDIR/data
  CBWADIR=$CNUMBERDIR/bwa
  CSUBREADDIR=$CNUMBERDIR/subread
  CANALYSISDIR=$CNUMBERDIR/run-analysis
}

# Create output directories for the species directory
# FIXME: we did not used $BASEDIR/data
function batch2-output {
  mkdir -p $DATADIR
  mkdir -p $BWADIR
  mkdir -p $SUBREADDIR
  mkdir -p $RUNANALYSIS
}

# Set species file variables
function batch2-speciesfile {
  # . fastq files
  # . adapter sequences
  # . quality score scheme
  # . cram data directory in _CRAMDIR_
  # . fastq file indices in _FASTQFILES_
  FASTQFILES=$(grep ^FASTQFILES\: $SPECIESFILE | cut -d":" -f2)
  # . fastq file labels in _FASTQLABEL_
  FASTQLABEL=$(grep ^FASTQLABEL\: $SPECIESFILE | cut -d":" -f2)
  # . cluster working directory in _CACWORKDIR_
  CACWORKDIR=$(grep ^CACWORKDIR\: $SPECIESFILE | cut -d":" -f2)
  # . reference genome ID in _REFGENOMEID_
  REFGENOMEID=$(grep ^REFGENOMEID\: $SPECIESFILE | cut -d":" -f2)
  # . reference genome genbank file in _REFGENOMEGENBANK_
  REFGENOMEGENBANK=$ROOTANALYSISDIR/$(grep ^REFGENOMEGENBANK\: $SPECIESFILE | cut -d":" -f2)
  GENOMEGENBANK=$(basename $REFGENOMEGENBANK)
  # . reference genome fasta file in _REFGENOMEFASTA_
  REFGENOMEFASTA=$ROOTANALYSISDIR/$(grep ^REFGENOMEFASTA\: $SPECIESFILE | cut -d":" -f2)
  # . reference genome annotation file in _REFGENOMEGFF_
  REFGENOMEGFF=$ROOTANALYSISDIR/$(grep ^REFGENOMEGFF\: $SPECIESFILE | cut -d":" -f2)
  # . reference genome annotation TrancriptDb file in _REFGENOMETXDB_
  REFGENOMETXDB=$ROOTANALYSISDIR/$(grep ^REFGENOMETXDB\: $SPECIESFILE | cut -d":" -f2)
  REFGENOMETXDBBASE=$(basename $REFGENOMETXDB)
  # . cram reference genome fasta file in _CRAMGENOMEFASTA_
  CRAMGENOMEFASTA=$ROOTANALYSISDIR/$(grep ^CRAMGENOMEFASTA\: $SPECIESFILE | cut -d":" -f2)
  # . cram wall time in _CRAMWALLTIME_
  CRAMWALLTIME=$(grep ^CRAMWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  # . pipeline wall time in _QCALIGNDEWALLTIME_
  QCALIGNDEWALLTIME=$(grep ^QCALIGNDEWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  # . minimum alignment mapping quality in _MINMAPQ_
  MINMAPQ=$(grep ^MINMAPQ\: $SPECIESFILE | cut -d":" -f2)
  # . minimum base quality to trim 3' end in _MINTRIMQ_
  MINTRIMQ=$(grep ^MINTRIMQ\: $SPECIESFILE | cut -d":" -f2)
  # . set bwa options in _BWAOPTION_
  BWAOPTION=$(grep ^BWAOPTION\: $SPECIESFILE | cut -d":" -f2)
  # . Rscript path in _CACRSCRIPT_
  # CACRSCRIPT=$(grep ^CACRSCRIPT\: $SPECIESFILE | cut -d":" -f2)
  # . set _SINGLECHROMOSOME_ to TRUE or FALSE
  SINGLECHROMOSOME=$(grep ^SINGLECHROMOSOME\: $SPECIESFILE | cut -d":" -f2)
  # . set _TESTFASTQ_ to integers
  TESTFASTQ=$(grep ^TESTFASTQ\: $SPECIESFILE | cut -d":" -f2)

  #################################################################
  # We may use these.
  BWAALIGNWALLTIME=$(grep ^BWAALIGNWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  CRAMDIR=$(grep ^CRAMDIR\: $SPECIESFILE | cut -d":" -f2)

  #################################################################
  # Many of the following variables are not needed. They may have to be moved to
  # somewhere else.
  REFGENOMENAME=$(grep ^REFGENOMENAME\: $SPECIESFILE | cut -d":" -f2)
  REFGENOMELENGTH=$(grep ^REFGENOMELENGTH\: $SPECIESFILE | cut -d":" -f2)
  REFGENOMEGBK=$(grep ^REFGENOMEGBK\: $SPECIESFILE | cut -d":" -f2)
  RNAZNNODE=$(grep ^RNAZNNODE\: $SPECIESFILE | cut -d":" -f2)
  RNAZWALLTIME=$(grep ^RNAZWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  PILEUPWALLTIME=$(grep ^PILEUPWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  PARSERNASEQNNODE=$(grep ^PARSERNASEQNNODE\: $SPECIESFILE | cut -d":" -f2)
  PARSERNASEQWALLTIME=$(grep ^PARSERNASEQWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  READDEPTH=$(grep ^READDEPTH\: $SPECIESFILE | cut -d":" -f2)

  REFGENOMEPTT=$ROOTANALYSISDIR/$(grep ^REFGENOMEPTT\: $SPECIESFILE | cut -d":" -f2)
  DBNAME=$(grep ^DBNAME\: $SPECIESFILE | cut -d":" -f2)
  CLADENAME=$(grep ^CLADENAME\: $SPECIESFILE | cut -d":" -f2)
  GENOMENAME=$(grep ^GENOMENAME\: $SPECIESFILE | cut -d":" -f2)
  ASSEMBLYNAME=$(grep ^ASSEMBLYNAME\: $SPECIESFILE | cut -d":" -f2)
}

function batch2-push-data {
cat>$BASEDIR/push-data.sh<<EOF
#!/bin/bash
ssh -x $CAC_USERHOST mkdir -p $RBWADIR
ssh -x $CAC_USERHOST mkdir -p $RSUBREADDIR
ssh -x $CAC_USERHOST mkdir -p $RDATADIR

scp $REFGENOMEFASTA $CAC_USERHOST:$RDATADIR
scp $REFGENOMEGFF $CAC_USERHOST:$RDATADIR
scp $CRAMGENOMEFASTA $CAC_USERHOST:$RDATADIR
scp $REFGENOMETXDB $CAC_USERHOST:$RDATADIR
scp $REFGENOMEPTT $CAC_USERHOST:$RDATADIR

echo Edit push-data.sh to send cram or fastq files.
for g in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$g)
  FASTQFILE=$ROOTANALYSISDIR/\$(grep ^\$FASTQNUM\: $SPECIESFILE | cut -d":" -f2)
  # scp \$FASTQFILE $CAC_USERHOST:$RDATADIR/\$FASTQNUM.fq.gz
done

EOF
}

function batch2-quality-score {
  QUALITYSCORELIST=""
  for g in $FASTQFILES; do
    QUALITYSCORENUM=QUALITYSCORE$(printf "%03d" $g)
    QUALITYSCORE=$(grep $QUALITYSCORENUM $SPECIESFILE | cut -d":" -f2)
    QUALITYSCORELIST="$QUALITYSCORELIST \"$QUALITYSCORE\""
  done
  QUALITYSCORELISTINR=$(echo $QUALITYSCORELIST | sed -e 's/[ ]/,/g')
}


function batch2-run-fastqc {
  STATUS=fastqc
  GENOMEFASTA=$(basename $REFGENOMEFASTA)
  
cat>$BASEDIR/run-$STATUS.sh<<EOF
#!/bin/bash
FASTQFILES=( $FASTQFILES )
sed s/PBSARRAYSIZE/\${#FASTQFILES[@]}/g < batch-$STATUS.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF

cat>$BASEDIR/batch-$STATUS.sh<<EOF
#!/bin/bash
#PBS -l walltime=${QCALIGNDEWALLTIME}:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-QC
#PBS -q ${QUEUENAME}
#PBS -m e
$EMAILON#PBS -M ${BATCHEMAIL}
#PBS -t 1-PBSARRAYSIZE

function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$HOME/$PRINSEQ pl
  cp \$HOME/$SAMTOOLS samtools
  cp \$HOME/$BWA bwa

  # All of the batchjob scripts.
  cp \$PBS_O_WORKDIR/job-fastqc .
  cp \$PBS_O_WORKDIR/job-stat* .

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CBWADIR

  # Copy common data
  cp $RDATADIR/$GENOMEFASTA $CDATADIR
}

function process-data {
  cd \$TMPDIR
  FASTQFILES=( $FASTQFILES )
  g=\$((PBS_ARRAYID-1))
  NUM=\$(printf "%03d" \${FASTQFILES[\$g]})

  cp $RDATADIR/FASTQ\$NUM.fq.gz $CDATADIR

  bash job-stat \$NUM $CDATADIR/FASTQ\$NUM.fq.gz \\
    $CBWADIR/FASTQ\$NUM.fq.qualPlot.RData
  cp $CBWADIR/FASTQ\$NUM.fq.qualPlot.RData $RBWADIR

  bash job-fastqc \$NUM \\
    $CDATADIR/FASTQ\$NUM.fq.gz \\
    $CBWADIR/FASTQ\$NUM.prinseq.fq.gz 
  # cp $CBWADIR/FASTQ\$NUM.prinseq.fq.gz $RBWADIR
  
  bash job-stat \$NUM $CBWADIR/FASTQ\$NUM.prinseq.fq.gz \\
    $CBWADIR/FASTQ\$NUM.prinseq.fq.qualPlot.RData
  cp $CBWADIR/FASTQ\$NUM.prinseq.fq.qualPlot.RData $RBWADIR
}

copy-data
process-data
cd
rm -rf \$TMPDIR
EOF
}

function batch2-run-cram {
  CRAMGENOMEFASTAFILENAME=$(basename $CRAMGENOMEFASTA)
  CRAMGENOMEFASTABASENAME=${CRAMGENOMEFASTAFILENAME%.fna}
cat>$BASEDIR/run-cram.sh<<EOF
#!/bin/bash
FASTQFILES=( $FASTQFILES )
sed s/PBSARRAYSIZE/\${#FASTQFILES[@]}/g < batch-cram.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF

cat>$BASEDIR/batch-cram.sh<<EOF
#!/bin/bash
#PBS -l walltime=${CRAMWALLTIME}:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-CRAM
#PBS -q ${QUEUENAME}
#PBS -m e
$EMAILON#PBS -M ${BATCHEMAIL}
#PBS -t 1-PBSARRAYSIZE

function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$PBS_O_WORKDIR/job-cram .
  cp \$HOME/$SAMTOOLS samtools
  cp \$HOME/$BWA bwa

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CBWADIR
}

function process-data {
  cd \$TMPDIR
  CORESPERNODE=1
  FASTQFILES=( $FASTQFILES )
  g=\$((PBS_ARRAYID-1))
  NUM=\$(printf "%03d" \${FASTQFILES[\$g]})

  cp $RBWADIR/FASTQ\$NUM.prinseq.fq.gz $CBWADIR/FASTQ\$NUM.fq.gz 
  bash job-cram \$NUM \\
    $CBWADIR/FASTQ\$NUM.fq.gz \\
    $CBWADIR/FASTQ\$NUM.cram
  cp $CBWADIR/FASTQ\$NUM.cram $RBWADIR
}

copy-data
process-data
cd
rm -rf \$TMPDIR
EOF
}

function batch2-run-bwa-align {
  GENOMEFASTA=$(basename $REFGENOMEFASTA)
  GENOMEGFF=$(basename $REFGENOMEGFF)
  status=bwa-align
cat>$BASEDIR/run-$status.sh<<EOF
#!/bin/bash
FASTQFILES=( $FASTQFILES )
sed s/PBSARRAYSIZE/\${#FASTQFILES[@]}/g < batch-$status.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF

cat>$BASEDIR/batch-$status.sh<<EOF
#!/bin/bash
#PBS -l walltime=${BWAALIGNWALLTIME}:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-BWA
#PBS -q ${QUEUENAME}
#PBS -m e
$EMAILON#PBS -M ${BATCHEMAIL}
#PBS -t 1-PBSARRAYSIZE

function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$HOME/$SAMTOOLS samtools
  cp \$HOME/$BWA bwa
  cp \$HOME/$SUBREADBUILDINDEX subread-buildindex
  cp \$HOME/$SUBREADALIGN subread-align

  # All of the batchjob scripts.
  cp \$PBS_O_WORKDIR/job-cram2fastq .
  cp \$PBS_O_WORKDIR/job-bwa-align .

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CBWADIR
  mkdir -p $CSUBREADDIR

  # Copy common data
  cp $RDATADIR/$GENOMEFASTA $CDATADIR
}

function process-data {
  cd \$TMPDIR
  CORESPERNODE=1
  FASTQFILES=( $FASTQFILES )
  g=\$((PBS_ARRAYID-1))
  NUM=\$(printf "%03d" \${FASTQFILES[\$g]})

  cp $RBWADIR/FASTQ\$NUM.recovered.fq.gz $CBWADIR
  bash job-bwa-align \$NUM \\
    $CBWADIR/FASTQ\$NUM.recovered.fq.gz \\
    $CBWADIR/FASTQ\$NUM.sorted \\
    $CSUBREADDIR/FASTQ\$NUM.sorted

  cp $CBWADIR/FASTQ\$NUM.sorted.bam $RBWADIR
}

copy-data
process-data
cd
rm -rf \$TMPDIR
EOF
}

function batch2-run-bias {
  GENOMEFASTA=$(basename $REFGENOMEFASTA)
  GENOMEGFF=$(basename $REFGENOMEGFF)
  status=bias
cat>$BASEDIR/run-$status.sh<<EOF
#!/bin/bash
FASTQFILES=( $FASTQFILES )
sed s/PBSARRAYSIZE/\${#FASTQFILES[@]}/g < batch-$status.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF

cat>$BASEDIR/batch-$status.sh<<EOF
#!/bin/bash
#PBS -l walltime=24:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-BIAS
#PBS -q ${QUEUENAME}
#PBS -m e
$EMAILON#PBS -M ${BATCHEMAIL}
#PBS -t 1-PBSARRAYSIZE

RSCRIPT=$CACRSCRIPT

function copy-data {
  cd \$TMPDIR
  cp \$PBS_O_WORKDIR/job-bias.R .
}

function process-data {
  cd \$TMPDIR
  FASTQFILES=( $FASTQFILES )
  g=\$((PBS_ARRAYID-1))
  NUM=\$(printf "%03d" \${FASTQFILES[\$g]})
  \$RSCRIPT job-bias.R FASTQ\$NUM
}

copy-data
process-data
cd
rm -rf \$TMPDIR
EOF

cat>$BASEDIR/job-$status.R<<EOF
library(seqbias) 
library(Rsamtools) 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{
  cat ("Rscript 1.R FASTQ001\\n")
  quit("no")
}

ref_fn <- "$RDATADIR/$GENOMEFASTA"
ref_f <- FaFile( ref_fn ) 
open.FaFile( ref_f )
bamFile <- sprintf("$RBWADIR/%s.sorted.bam", args[1])
indexBam(bamFile)
ref_seqs <- scanFaIndex( ref_f )

ref_seq <- getSeq(ref_f)
I.all <- GRanges(seqnames=Rle(names(ref_seq),c(2)),
                 ranges=IRanges(c(1,1),width=rep(width(ref_seq),2)),
                 strand=Rle(strand(c("+","-")),c(1,1)))
seqlengths(I.all) <- c(width(ref_seq))
sb <- seqbias.fit( ref_fn, bamFile, L = 5, R = 20 )
ymlFile <- sprintf("$RBWADIR/%s.yml", args[1])
seqbias.save(sb,ymlFile)
EOF
}

function batch2-run-coverage {
  GENOMEFASTA=$(basename $REFGENOMEFASTA)
  status=coverage
cat>$BASEDIR/run-$status.sh<<EOF
#!/bin/bash
FASTQFILES=( $FASTQFILES )
sed s/PBSARRAYSIZE/80/g < batch-$status.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF

cat>$BASEDIR/run-$status-pull.sh<<EOF
#!/bin/bash
FASTQFILES=( $FASTQFILES )
RSCRIPT=$CACRSCRIPT
for k in ${FASTQFILES[@]}; do
  TESTFASTQNUM=FASTQ\$(printf "%03d" \$k)
  \$RSCRIPT job-$status-pull.R \$TESTFASTQNUM
  rm $RBWADIR/\$TESTFASTQNUM.cvg.*
done
EOF

cat>$BASEDIR/batch-$status.sh<<EOF
#!/bin/bash
#PBS -l walltime=24:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-COV
#PBS -q ${QUEUENAME}
#PBS -m e
$EMAILON#PBS -M ${BATCHEMAIL}
#PBS -t 1-PBSARRAYSIZE

NUMSPLIT=\$((PBSARRAYSIZE * 8))
RSCRIPT=$CACRSCRIPT

function copy-data {
  cd \$TMPDIR
  cp \$PBS_O_WORKDIR/job-coverage.R .
  cp \$HOME/$SAMTOOLS samtools
  mkdir -p $CBWADIR
  mkdir -p $CDATADIR
  cp $RDATADIR/$GENOMEFASTA $CDATADIR
  ./samtools faidx $CDATADIR/$GENOMEFASTA 
}

function process-data {
  cd \$TMPDIR
  FASTQFILES=( $FASTQFILES )
  for k in \${FASTQFILES[@]}; do
    FASTQNUM=FASTQ\$(printf "%03d" \$k)
    cp $RBWADIR/\$FASTQNUM.sorted.bam $CBWADIR

    CORESPERNODE=8
    for (( i=1; i<=CORESPERNODE; i++)); do
      g=\$((CORESPERNODE * (PBS_ARRAYID-1) + i))
      cp $CBWADIR/\$FASTQNUM.sorted.bam $CBWADIR/\$FASTQNUM.\$i.sorted.bam
      \$RSCRIPT job-coverage.R \$FASTQNUM $CBWADIR/\$FASTQNUM.\$i.sorted.bam \$NUMSPLIT \$g &
    done
    wait

    rm $CBWADIR/\$FASTQNUM.*sorted.bam 
  done
}

copy-data
process-data
wait
cd
rm -rf \$TMPDIR
EOF


cat>$BASEDIR/job-coverage.R<<EOF
library(seqbias)
library(Rsamtools)
library(ShortRead)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4)
{
  cat ("Rscript 1.R FASTQ001 8 1\\n")
  quit("no")
}

# bamFile <- sprintf("$CBWADIR/%s.sorted.bam", args[2])
bamFile <- args[2]
aln <- readGappedAlignments(bamFile)
aln <- as(aln, "GRanges")

n <- as.integer(args[3]) + 1L
a <- as.integer(seq(1L,length(aln),length.out=n))
a[n] <- a[n] + 1
s.1 <- a[as.integer(args[4])]
s.2 <- a[as.integer(args[4]) + 1] - 1

ref_fn <- "$CDATADIR/$GENOMEFASTA"
ref_f <- FaFile( ref_fn )
open.FaFile( ref_f )
ref_seqs <- scanFaIndex( ref_f )
ref_seq <- getSeq(ref_f)
I.all <- GRanges(seqnames=Rle(c("chr1"),c(2)),
                 ranges=IRanges(c(1,1),width=rep(width(ref_seq),2)),
                 strand=Rle(strand(c("+","-")),c(1,1)))
seqlengths(I.all) <- c(width(ref_seq))

ymlFile <- sprintf("$RBWADIR/%s.yml", args[1])
sb <- seqbias.load( ref_fn, ymlFile)
bias <- seqbias.predict( sb, I.all )

s1 <- start(aln)
s2 <- end(aln)
s3 <- strand(aln)

cvg <- rep(0,width(ref_seq))
for (i in s.1:s.2) {
  if (runValue(s3[i]) == "+") {
    cvg[s1[i]:s2[i]] <- cvg[s1[i]:s2[i]] + 1/bias[[1]][s1[i]]
  } else {
    cvg[s1[i]:s2[i]] <- cvg[s1[i]:s2[i]] + 1/bias[[2]][s2[i]]
  }
}
cvgFile <- sprintf("$RBWADIR/%s.cvg.%s",args[1],args[4])
save(cvg,file=cvgFile)
close.FaFile( ref_f )
EOF

cat>$BASEDIR/job-coverage-pull.R<<EOF
library(seqbias)
library(Rsamtools)
library(ShortRead)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{
  cat ("Rscript 1.R FASTQ001\\n")
  quit("no")
}
ref_fn <- "$RDATADIR/$GENOMEFASTA"
ref_f <- FaFile( ref_fn )
open.FaFile( ref_f )
ref_seqs <- scanFaIndex( ref_f )
ref_seq <- getSeq(ref_f)

a <- list.files(path = "$RBWADIR", pattern=paste(args[1],".cvg.",sep="") , full.names = TRUE)

cvg.org <- rep(0,width(ref_seq))
for (i in a) {
  load(i)
  cvg.org <- cvg.org + cvg
}
cvgFile <- sprintf("$RBWADIR/%s.coverage.RData",args[1])
save(cvg.org,file=cvgFile)

close.FaFile( ref_f )

wigFile <- sprintf("$RBWADIR/%s.wig",args[1])

cat("track type=wiggle_0 name=\\"RNA-seq\\" description=\\"RNA-seq\\" visibility=full autoScale=on color=0,200,100 priority=30\\n",
    file=wigFile)
cat("fixedStep chrom=chr1 start=1 step=1\\n",
    file=wigFile,append=TRUE)
write(cvg.org,file=wigFile,ncolumns=1,append=TRUE)
EOF

cat>$BASEDIR/job-$status-start-end.sh<<EOF
#!/bin/bash
FASTQFILES=( $FASTQFILES )
RSCRIPT=$CACRSCRIPT
for k in \${FASTQFILES[@]}; do
  TESTFASTQNUM=FASTQ\$(printf "%03d" \$k)
  \$RSCRIPT job-$status-start-end.R \$TESTFASTQNUM $RBWADIR/\$TESTFASTQNUM.sorted.bam
done
EOF
cat>$BASEDIR/job-$status-start-end.R<<EOF
library(seqbias)
library(Rsamtools)
library(ShortRead)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2)
{
  cat ("Rscript 1.R FASTQ001 FASTQ001.bam\\n")
  quit("no")
}

bamFile <- args[2]
aln <- readGappedAlignments(bamFile)
aln <- as(aln, "GRanges")

ref_fn <- "$RDATADIR/$GENOMEFASTA"
ref_f <- FaFile( ref_fn )
open.FaFile( ref_f )
ref_seqs <- scanFaIndex( ref_f )
ref_seq <- getSeq(ref_f)
I.all <- GRanges(seqnames=Rle(c("chr1"),c(2)),
                 ranges=IRanges(c(1,1),width=rep(width(ref_seq),2)),
                 strand=Rle(strand(c("+","-")),c(1,1)))
seqlengths(I.all) <- c(width(ref_seq))

ymlFile <- sprintf("$RBWADIR/%s.yml", args[1])
sb <- seqbias.load( ref_fn, ymlFile)
bias <- seqbias.predict( sb, I.all )

s1 <- start(aln)
s2 <- end(aln)
s3 <- strand(aln)

lenseq <- width(ref_seq)
s1 <- start(aln)
s1 <- s1[as.logical(s3 == '+')]
s1 <- IRanges(s1, width = 1)
s1 <- coverage(s1)
s1 <- as.vector(s1) 
s1 <- c(s1,rep(0,lenseq - length(s1)))
s1 <- s1/bias[[1]]
s2 <- start(aln)
s2 <- s2[as.logical(s3 == '-')]
s2 <- IRanges(s2, width = 1)
s2 <- coverage(s2)
s2 <- as.vector(s2)
s2 <- c(s2,rep(0,lenseq - length(s2)))
s2 <- s2/bias[[2]]
cvgStart <- s1 + s2

map.readFile <- sprintf("$RBWADIR/%s.start",args[1])
write(cvgStart,file=map.readFile)
close.FaFile( ref_f )
print(paste("Check",map.readFile))
EOF

}

function batch2-run-samtools-pileup {
  GENOMEFASTA=$(basename $REFGENOMEFASTA)
  GENOMEGFF=$(basename $REFGENOMEGFF)
  status=samtools-pileup
cat>$BASEDIR/run-$status.sh<<EOF
#!/bin/bash
FASTQFILES=( $FASTQFILES )
sed s/PBSARRAYSIZE/\${#FASTQFILES[@]}/g < batch-$status.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF

cat>$BASEDIR/batch-$status.sh<<EOF
#!/bin/bash
#PBS -l walltime=${PILEUPWALLTIME}:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-PILEUP
#PBS -q ${QUEUENAME}
#PBS -m e
$EMAILON#PBS -M ${BATCHEMAIL}
#PBS -t 1-PBSARRAYSIZE

function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$HOME/$SAMTOOLS samtools

  # All of the batchjob scripts.
  cp \$PBS_O_WORKDIR/job-$status .

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CBWADIR
  mkdir -p $CSUBREADDIR

  # Copy common data
  cp $RDATADIR/$GENOMEFASTA $CDATADIR
}

function process-data {
  cd \$TMPDIR
  CORESPERNODE=1
  FASTQFILES=( $FASTQFILES )
  g=\$((PBS_ARRAYID-1))
  NUM=\$(printf "%03d" \${FASTQFILES[\$g]})

  cp $RBWADIR/FASTQ\$NUM.sorted.bam $CBWADIR
  cp $RSUBREADDIR/FASTQ\$NUM.sorted.bam $CSUBREADDIR
  bash job-$status \$NUM \\
    $CBWADIR/FASTQ\$NUM.sorted.bam \\
    $CBWADIR/FASTQ\$NUM.mpileup
  cp $CBWADIR/FASTQ\$NUM.mpileup $RBWADIR
}

copy-data
process-data
cd
rm -rf \$TMPDIR
EOF
}

function create-index {
  echo $FASTQLABEL > $RUNANALYSIS/count.txt.index
}

function batch2-run-cram2fastq {
  STATUS=cram2fastq
  GENOMEFASTA=$(basename $REFGENOMEFASTA)
  GENOMEGFF=$(basename $REFGENOMEGFF)
  CRAMGENOMEFASTAFILENAME=$(basename $CRAMGENOMEFASTA)
  CRAMGENOMEFASTABASENAME=${CRAMGENOMEFASTAFILENAME%.fna}
  
cat>$BASEDIR/run-$STATUS.sh<<EOF
#!/bin/bash
FASTQFILES=( $FASTQFILES )
sed s/PBSARRAYSIZE/\${#FASTQFILES[@]}/g < batch-$STATUS.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF

cat>$BASEDIR/batch-$STATUS.sh<<EOF
#!/bin/bash
#PBS -l walltime=${QCALIGNDEWALLTIME}:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-$STATUS
#PBS -q ${QUEUENAME}
#PBS -m e
$EMAILON#PBS -M ${BATCHEMAIL}
#PBS -t 1-PBSARRAYSIZE


function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$HOME/$PRINSEQ pl
  cp \$HOME/$SAMTOOLS samtools
  cp \$HOME/$BWA bwa

  # All of the batchjob scripts.
  cp \$PBS_O_WORKDIR/job-cram2fastq .

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CBWADIR

  # Copy common data
  cp $RDATADIR/$GENOMEFASTA $CDATADIR
}

function process-data {
  cd \$TMPDIR
  FASTQFILES=( $FASTQFILES )
  g=\$((PBS_ARRAYID-1))
  NUM=\$(printf "%03d" \${FASTQFILES[\$g]})
  # input: cram, and output: bam

  cp $RBWADIR/FASTQ\$NUM.cram $CBWADIR
  bash job-cram2fastq \$NUM \\
    $CBWADIR/FASTQ\$NUM.cram \\
    $CBWADIR/FASTQ\$NUM.recovered.fq
  gzip $CBWADIR/FASTQ\$NUM.recovered.fq

  cp $CBWADIR/FASTQ\$NUM.recovered.fq.gz $RBWADIR
}

copy-data
process-data
cd
rm -rf \$TMPDIR
EOF
}

# This creates a shell script for running a submission of a job.
# For example, run-bwa-align.sh.
function batch2-make-run-dot-sh {
cat>$BASEDIR/run-$1.sh<<EOF
#!/bin/bash
FASTQFILES=( $FASTQFILES )
sed s/PBSARRAYSIZE/\${#FASTQFILES[@]}/g < batch-$1.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF
}

# We modify a copy of batch2-run-qcalignde function.
function batch2-run-tux {
  STATUS=tux
  GENOMEFASTA=$(basename $REFGENOMEFASTA)
  batch2-make-run-dot-sh $STATUS
cat>$BASEDIR/batch-$STATUS.sh<<EOF
#!/bin/bash
#PBS -l walltime=${QCALIGNDEWALLTIME}:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-tux
#PBS -q ${QUEUENAME}
#PBS -m e
$EMAILON#PBS -M ${BATCHEMAIL}
#PBS -t 1-PBSARRAYSIZE

function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$HOME/$PRINSEQ pl
  cp \$HOME/$SAMTOOLS samtools
  cp \$HOME/$BWA bwa

  # All of the batchjob scripts.
  cp \$PBS_O_WORKDIR/job-fastqc .
  cp \$PBS_O_WORKDIR/job-bwa-align .
  cp \$PBS_O_WORKDIR/job-tux* .

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CBWADIR
  mkdir -p $CSUBREADDIR

  # Copy common data
  cp $RDATADIR/$GENOMEFASTA $CDATADIR
}

function process-data {
  cd \$TMPDIR
  FASTQFILES=( $FASTQFILES )
  g=\$((PBS_ARRAYID-1))
  NUM=\$(printf "%03d" \${FASTQFILES[\$g]})

  cp $RDATADIR/FASTQ\$NUM.fq.gz $CDATADIR

  bash job-fastqc \$NUM \\
    $CDATADIR/FASTQ\$NUM.fq.gz \\
    $CBWADIR/FASTQ\$NUM.prinseq.fq.gz 

  NUMBER_READPRIN4=\$(zcat $CBWADIR/FASTQ\$NUM.prinseq.fq.gz | wc -l)
  NUMBER_READPRIN=\$((NUMBER_READPRIN4 / 4))
   
  # 2. Align the sort reads
  bash job-bwa-align \$NUM \\
    $CBWADIR/FASTQ\$NUM.prinseq.fq.gz \\
    $CBWADIR/FASTQ\$NUM.sorted \\
    $CSUBREADDIR/FASTQ\$NUM.sorted

  # 3. Count short reads
  BAMFILE1=$CBWADIR/FASTQ\$NUM.sorted.bam
  NUMBER_READ4=\$(zcat $CDATADIR/FASTQ\$NUM.fq.gz | wc -l)
  NUMBER_READ=\$((NUMBER_READ4 / 4))
  bash job-tux \$BAMFILE1 $CDATADIR/$GENOMEFASTA $CBWADIR/FASTQ\$NUM.bam.tux
  cp $CBWADIR/FASTQ\$NUM.bam.tux $RBWADIR

  echo -en "  - [" > 2.tmp
  cat $CBWADIR/FASTQ\$NUM.bam.tux | tr "\\n" "," >> 2.tmp
  sed s/,$// 2.tmp > $CBWADIR/FASTQ\$NUM.bam.tux.array
  echo -en "]\n" >> $CBWADIR/FASTQ\$NUM.bam.tux.array
  rm 2.tmp
  cp $CBWADIR/FASTQ\$NUM.bam.tux.array $RBWADIR

  echo "============================="
  echo "The content of BWA directory:"
  ls $CBWADIR
  echo "============================="
}

copy-data
process-data
cd
rm -rf \$TMPDIR
EOF
}

function batch2-run-qcalignde {
  STATUS=qcalignde 
  GENOMEFASTA=$(basename $REFGENOMEFASTA)
  GENOMEGFF=$(basename $REFGENOMEGFF)
  CRAMGENOMEFASTAFILENAME=$(basename $CRAMGENOMEFASTA)
  CRAMGENOMEFASTABASENAME=${CRAMGENOMEFASTAFILENAME%.fna}
  
cat>$BASEDIR/run-$STATUS.sh<<EOF
#!/bin/bash
FASTQFILES=( $FASTQFILES )
sed s/PBSARRAYSIZE/\${#FASTQFILES[@]}/g < batch-$STATUS.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF

cat>$BASEDIR/batch-$STATUS.sh<<EOF
#!/bin/bash
#PBS -l walltime=${QCALIGNDEWALLTIME}:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-QAD
#PBS -q ${QUEUENAME}
#PBS -m e
$EMAILON#PBS -M ${BATCHEMAIL}
#PBS -t 1-PBSARRAYSIZE

function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$HOME/$PRINSEQ pl
  cp \$HOME/$SAMTOOLS samtools
  cp \$HOME/$BWA bwa

  # All of the batchjob scripts.
  cp \$PBS_O_WORKDIR/job-cram .
  cp \$PBS_O_WORKDIR/job-cram2fastq .
  cp \$PBS_O_WORKDIR/job-fastqc .
  cp \$PBS_O_WORKDIR/job-bwa-align .
  cp \$PBS_O_WORKDIR/job-de* .
  cp \$PBS_O_WORKDIR/job-de*.R .
  cp \$PBS_O_WORKDIR/job-stat* .

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CBWADIR
  mkdir -p $CSUBREADDIR

  # Copy common data
  cp $RDATADIR/$GENOMEFASTA $CDATADIR
  cp $RDATADIR/$GENOMEGFF $CDATADIR
}

function process-data {
  cd \$TMPDIR
  FASTQFILES=( $FASTQFILES )
  g=\$((PBS_ARRAYID-1))
  NUM=\$(printf "%03d" \${FASTQFILES[\$g]})

  cp $RDATADIR/FASTQ\$NUM.fq.gz $CDATADIR

  bash job-stat \$NUM $CDATADIR/FASTQ\$NUM.fq.gz \\
    $CBWADIR/FASTQ\$NUM.fq.qualPlot.RData
  cp $CBWADIR/FASTQ\$NUM.fq.qualPlot.RData $RBWADIR

  bash job-fastqc \$NUM \\
    $CDATADIR/FASTQ\$NUM.fq.gz \\
    $CBWADIR/FASTQ\$NUM.prinseq.fq.gz 

#  bash job-stat \$NUM $CBWADIR/FASTQ\$NUM.prinseq.fq.gz \\
#    $CBWADIR/FASTQ\$NUM.prinseq.fq.qualPlot.RData
#  cp $CBWADIR/FASTQ\$NUM.prinseq.fq.qualPlot.RData $RBWADIR


  NUMBER_READPRIN4=\$(zcat $CBWADIR/FASTQ\$NUM.prinseq.fq.gz | wc -l)
  NUMBER_READPRIN=\$((NUMBER_READPRIN4 / 4))
   
  # 2. Align the sort reads
  bash job-bwa-align \$NUM \\
    $CBWADIR/FASTQ\$NUM.prinseq.fq.gz \\
    $CBWADIR/FASTQ\$NUM.sorted \\
    $CSUBREADDIR/FASTQ\$NUM.sorted

  # 3. Count short reads
  BAMFILE1=$CBWADIR/FASTQ\$NUM.sorted.bam
  NUMBER_READ4=\$(zcat $CDATADIR/FASTQ\$NUM.fq.gz | wc -l)
  NUMBER_READ=\$((NUMBER_READ4 / 4))
  bash job-de \$BAMFILE1 \$NUMBER_READ \$NUMBER_READPRIN
#  cp \$BAMFILE1 $RBWADIR
  cp \$BAMFILE1.cl $RBWADIR
}

copy-data
process-data
cd
rm -rf \$TMPDIR
EOF

# $1: fastq file in gzip
# $2: RData
# $3: FASTQ number
grep ^QUALITYSCORE $SPECIESFILE | sed s/:/=/ > $BASEDIR/job-stat
cat>>$BASEDIR/job-stat<<EOF
QUALITYSCORESEQUENCE=QUALITYSCORE\$1

# We need more preprocessing of fastq files.
if [ "\${!QUALITYSCORESEQUENCE}" == "illumina" ]; then
  PHRED64=illumina
else
  PHRED64=sanger
fi

RSCRIPT=$CACRSCRIPT
\$RSCRIPT job-stat.R \$2 \$3 \$PHRED64
EOF

cat>$BASEDIR/job-stat.R<<EOF
library(qrqc)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3)
{
  cat ("Rscript job-stat.R 1.fq.gz 1.fq.RData sanger\n")
  quit("yes")
}
fq.name <- args[1]
fq.quality <- args[3]
fq.file <- readSeqFile(fq.name,quality=fq.quality)
toplot <- qualPlot(fq.file)
fq.plot <- args[2]
save(list="toplot", file = fq.plot)
EOF

# bash job-stat-plot
cat>$BASEDIR/job-stat-plot<<EOF
RSCRIPT=$CACRSCRIPT
for i in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$i)
  \$RSCRIPT job-stat-plot.R $RBWADIR/\$FASTQNUM.fq.qualPlot.RData \\
    $RBWADIR/\$FASTQNUM.fq.qualPlot.pdf
  \$RSCRIPT job-stat-plot.R $RBWADIR/\$FASTQNUM.prinseq.fq.qualPlot.RData \\
    $RBWADIR/\$FASTQNUM.prinseq.fq.qualPlot.pdf
done 
EOF

cat>$BASEDIR/job-stat-plot.R<<EOF
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2)
{
  cat ("Rscript job-stat-plot.R 1.RData 2.pdf\n")
  quit("yes")
}
fq.plot <- args[1]
fq.pdf <- args[2]
load(fq.plot)
pdf(fq.pdf)
plot(toplot)
dev.off()
EOF

cat>$BASEDIR/feature-txnc.txt<<EOF
###################################################################
# 4. Number of short reads uniquely mapped on CDS regions
# Map the uniquely mapped reads on CDS regions
# Remove strands information
# 4.1 Reads mapped on CDS tx
# 4.2 Reads mapped on non CDS tx
# 4.3 Reads mapped on non-genic tx
# 4.4 Reads mapped on tx and non-genic tx -> Is this all?
feature.tx <- transcripts(txdb)
feature.cds <- cds(txdb,columns="tx_id")
x <- elementMetadata(feature.tx)\$tx_id %in% unlist(elementMetadata(feature.cds)\$tx_id)
# 4.1 Reads mapped on CDS tx
feature.cds <- feature.tx[x]
elementMetadata(feature.cds) <- data.frame(tx_id=elementMetadata(feature.cds)\$tx_id, tx_name=elementMetadata(feature.cds)\$tx_name,type="CDS")
feature.txnc <- feature.cds

# 4.2 Reads mapped on non CDS tx
feature.nocds <- feature.tx[!x]
if (length(feature.nocds) > 0) {
  elementMetadata(feature.nocds) <- data.frame(tx_id=elementMetadata(feature.nocds)\$tx_id, tx_name=elementMetadata(feature.nocds)\$tx_name,type="NOCDS")
  feature.txnc <- c(feature.txnc,feature.nocds)
}

# 4.3 Reads mapped on non-genic tx
feature.ng <- transcripts(txdb)
strand(feature.ng) <- '*'
feature.ng <- gaps(reduce(feature.ng))
feature.ng <- feature.ng[strand(feature.ng)=='*']
# 4.4 Reads mapped on tx and non-genic tx
if (length(feature.ng) > 0) {
  x <- seq(from=length(feature.tx)+1,to=length(feature.tx)+length(feature.ng))
  y <- paste("NC",x,sep="_")
  elementMetadata(feature.ng) <- data.frame(tx_id=x,tx_name=y,type="NG")
  rm(x,y)
  feature.txnc <- c(feature.txnc,feature.ng)
}
###################################################################
EOF

cat>$BASEDIR/job-check<<EOF
RSCRIPT=$CACRSCRIPT

for i in $FASTQFILES; do
  TESTFASTQNUM=FASTQ\$(printf "%03d" \$i)
  \$RSCRIPT job-check.R \$TESTFASTQNUM
done 

EOF

cat>$BASEDIR/job-check.R<<EOF
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{
  cat ("Rscript job-simulate.R FASTQ001\n")
  quit("yes")
}
load(sprintf("$RBWADIR/%s.sorted.bam.cl", args[1]))
cl.bwa <- cl
rm(cl)

load(sprintf("$RBWADIR/%s.cl", args[1]))
cl.sim <- cl
rm(cl)

print(paste("Rate of correction for count:",sum(cl.bwa==cl.sim)/length(cl.sim)))
q("no")
x <- scan(sprintf("$RBWADIR/%s.tux", args[1]))
y <- scan(sprintf("$RBWADIR/%s.bam.tux", args[1]))
print(sum(x==y)/length(x))
print(paste("Rate of correction for position of reads:",sum(x==y)/length(x)))

EOF

cat>$BASEDIR/job-simulate<<EOF
RSCRIPT=$CACRSCRIPT
for i in $TESTFASTQ; do
  TESTFASTQNUM=FASTQ\$(printf "%03d" \$i)
  rm $RDATADIR/\$TESTFASTQNUM.fq.gz
  \$RSCRIPT job-simulate.R \$TESTFASTQNUM
  gzip $RDATADIR/\$TESTFASTQNUM.fq
  echo Created File: $RDATADIR/\$TESTFASTQNUM.fq.gz
done 
EOF

cat>$BASEDIR/job-simulate.R<<EOF
############################################################################
# Create a test fastq files by sampling 10 100-bp short reads for each gene.
library(DESeq)
library(ShortRead)
library(rtracklayer)
library(GenomicRanges)
library(VariantAnnotation)
library(GenomicFeatures)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{
  cat ("Rscript job-simulate.R FASTQ001\n")
  quit("yes")
}
cl.file <- sprintf("$RBWADIR/%s.cl", args[1])
tux.file <- sprintf("$RBWADIR/%s.tux", args[1])

txdb.file <- "$RDATADIR/$REFGENOMETXDBBASE"


# saveFeatures and loadFeatures
txdb <- loadFeatures(txdb.file)
EOF

cat $BASEDIR/feature-txnc.txt >> $BASEDIR/job-simulate.R 
cat>>$BASEDIR/job-simulate.R<<EOF
strand(feature.txnc) <- '*'

genome.file <- "$RDATADIR/$GENOMEFASTA"
fastq.file <- sprintf("$RDATADIR/%s.fq", args[1])
file.create(fastq.file)
read.GR <- GRanges()
read.GA <- GappedAlignments()

read.loc <- function(x) { 
  s <- sample(c(TRUE,FALSE),size=1) 
  y <- x + 99
  if (s == FALSE) {
    y <- x - 99
  }
  if (y < 1) {
    y <- 1
  }
  if (y > length(s.mutans)) {
    y <- length(s.mutans)
  }
  z <- strand("+")
  if (x > y) {
    z <- x
    x <- y
    y <- z
    z <- strand("-") 
  }
  list(x,y,z)
}

read.extract <- function (x,y,z,w) {
  s <- subseq(s.mutans,x,y)
  if (z == "-") {
    s <- reverseComplement(s)
  }
  s.quality <- paste((rep("B",times=length(s))),collapse="")
  oneRead <- sprintf("@HWI-ST397:%09d 1:N:0:ACAGTG\\n%s\\n+\\n%s\\n", w, s, s.quality)
  cat(oneRead, file=fastq.file, append=TRUE)
}

s.mutans.sequence <- read.DNAStringSet(genome.file)
chrom.list <- names(s.mutans.sequence)

lengthChr1 <- 0
for (i in names(s.mutans.sequence)) { 
  cat (i,"\\n")

  geneIR <- ranges(feature.txnc[seqnames(feature.txnc)==i])

  geneIR <- geneIR[width(geneIR)>1000]
  geneIR <- geneIR - 100
  if (length(geneIR) > 0) {
    read.pos <- c(mapply(function(x,y) sample(x:y,size=10,replace=TRUE),start(geneIR),end(geneIR)))
    s.mutans <- s.mutans.sequence[[i]]
    lengthChr1 <- length(s.mutans)

    read.start.end.strand <- mapply(read.loc, read.pos)
    
    one.GA <- GappedAlignments( seqnames = Rle(i, length(read.pos)),
                                pos = as.integer(unlist(read.start.end.strand[1,])),
                                cigar = rep("100M", length(read.pos)),
                                strand = unlist(read.start.end.strand[3,]) )
    read.GA <- c(read.GA,one.GA)

    read.simulated <- mapply(read.extract, 
                             unlist(read.start.end.strand[1,]),
                             unlist(read.start.end.strand[2,]), 
                             unlist(read.start.end.strand[3,]),
                             seq(length(unlist(read.start.end.strand[3,]))))
  }
}
olap <- summarizeOverlaps(feature.txnc,read.GA,mode="IntersectionStrict")
cl <- assays(olap)\$counts[,1]
save(cl,file=cl.file)
print(paste("See cl file",cl.file))
q("no")

# We assume that a single chromosome genome.
stopifnot( length(names(s.mutans.sequence)) == 1 )
count.read <- rep.int(x=0,times=lengthChr1)
bv <- read.GA
a <- table(c(start(bv[strand(bv)=='+']), end(bv[strand(bv)=='-'])))
count.read[as.integer(names(a))] <- as.vector(a)
print(paste("Creating tux file",tux.file))
write(count.read,file=tux.file,ncolumns=1)

print(paste("See tux file",tux.file))
EOF

# CRAM
# $1: a three-digit number
# $2: a gzipped fastq file
# $3: a cram file
grep ^QUALITYSCORE $SPECIESFILE | sed s/:/=/ > $BASEDIR/job-cram
cat>>$BASEDIR/job-cram<<EOF
#  bash job-cram \$NUM \\
#    $CBWADIR/FASTQ\$NUM.fq.gz \\
#    $RDATADIR/FASTQ\$NUM.cram
QUALITYSCORESEQUENCE=QUALITYSCORE\$1

cp $RDATADIR/$CRAMGENOMEFASTAFILENAME $CDATADIR
cp $CDATADIR/$CRAMGENOMEFASTABASENAME.fna $CDATADIR/$CRAMGENOMEFASTABASENAME.fa

FASTAQFILE=\${GZIPFASTAQFILE%.gz}

./bwa index -p $CBWADIR/$CRAMGENOMEFASTABASENAME-bwa -a is \\
  $CDATADIR/$CRAMGENOMEFASTABASENAME.fa

# We need more preprocessing of fastq files.
if [ "\${!QUALITYSCORESEQUENCE}" == "illumina" ]; then
  cp \$2 $CBWADIR/temp.FASTQ\$1.fq.gz
else
  zcat \$2 \\
    | grep -A 3 '^@.* [^:]*:N:[^:]*:' \\
    | sed '/^--$/d' | gzip > $CBWADIR/temp.FASTQ\$1.fq.gz
fi

if [ "\${!QUALITYSCORESEQUENCE}" == "illumina" ]; then
  ./bwa aln -I -t $NUMBERCPU \\
    $CBWADIR/$CRAMGENOMEFASTABASENAME-bwa \\
    $CBWADIR/temp.FASTQ\$1.fq.gz > $CBWADIR/FASTQ\$1.sai
else
  ./bwa aln -t $NUMBERCPU \\
    $CBWADIR/$CRAMGENOMEFASTABASENAME-bwa \\
    $CBWADIR/temp.FASTQ\$1.fq.gz > $CBWADIR/FASTQ\$1.sai
fi
./bwa samse -n 1 \\
  -r "@RG\\tID:$PROJECTNAME\\tSM:CRAM" \\
  $CBWADIR/$CRAMGENOMEFASTABASENAME-bwa \\
  $CBWADIR/FASTQ\$1.sai \\
  $CBWADIR/temp.FASTQ\$1.fq.gz \\
  | ./samtools view -Sb - > $CBWADIR/FASTQ\$1.bam
rm $CBWADIR/temp.FASTQ\$1.fq.gz

./samtools sort $CBWADIR/FASTQ\$1.bam \\
  $CBWADIR/FASTQ\$1.sorted
rm $CBWADIR/FASTQ\$1.bam

./samtools faidx $CDATADIR/$CRAMGENOMEFASTABASENAME.fa
./samtools index $CBWADIR/FASTQ\$1.sorted.bam

java -jar \$HOME/$CRAMTOOLS cram \\
  --capture-all-quality-scores \\
  --include-unmapped-reads \\
  --capture-unmapped-quality-scores \\
  --capture-all-tags \\
  --input-bam-file $CBWADIR/FASTQ\$1.sorted.bam \\
  --reference-fasta-file $CDATADIR/$CRAMGENOMEFASTABASENAME.fa \\
  --output-cram-file $CDATADIR/FASTQ\$1.cram
rm $CBWADIR/FASTQ\$1.sorted.bam

cp $CDATADIR/FASTQ\$1.cram \$3

EOF

# Unzip fastq files
# $1: a three-digit number
# $2: a cram file
# $3: a fastq file
cat>$BASEDIR/job-cram2fastq<<EOF
cp $RDATADIR/$CRAMGENOMEFASTAFILENAME $CDATADIR
cp $CDATADIR/$CRAMGENOMEFASTABASENAME.fna $CDATADIR/$CRAMGENOMEFASTABASENAME.fa

./bwa index -p $CBWADIR/$CRAMGENOMEFASTABASENAME-bwa -a is \\
  $CDATADIR/$CRAMGENOMEFASTABASENAME.fa &> /dev/null
./samtools faidx $CDATADIR/$CRAMGENOMEFASTABASENAME.fa &> /dev/null

java -jar \$HOME/$CRAMTOOLS bam \\
  --input-cram-file \$2 \\
  --reference-fasta-file $CDATADIR/$CRAMGENOMEFASTABASENAME.fa \\
  --output-bam-file $CDATADIR/FASTQ\$1.bam &> /dev/null

java -Xmx10g -jar \\
  \$HOME/$PICARD/SamToFastq.jar \\
  INPUT=$CDATADIR/FASTQ\$1.bam \\
  FASTQ=\$3 &> /dev/null

rm $CDATADIR/FASTQ\$1.bam
EOF

# QC
# $1: a three-digit number
# $2: a gzipped fastq file
# $3: a QC'ed gzipped fastq file
grep ^ADAPTER $SPECIESFILE | sed s/:/=/ > $BASEDIR/job-fastqc
grep ^QUALITYSCORE $SPECIESFILE | sed s/:/=/ >> $BASEDIR/job-fastqc
cat>>$BASEDIR/job-fastqc<<EOF
ADAPTERSEQUENCE=ADAPTER\$1
QUALITYSCORESEQUENCE=QUALITYSCORE\$1

# We need more preprocessing of fastq files.
if [ "\${!QUALITYSCORESEQUENCE}" == "illumina" ]; then
  PHRED64=-phred64
  cp \$2 $CBWADIR/temp.FASTQ\$1.fq.gz
else
  PHRED64=
  zcat \$2 \\
    | grep -A 3 '^@.* [^:]*:N:[^:]*:' \\
    | sed '/^--$/d' | gzip > $CBWADIR/temp.FASTQ\$1.fq.gz
fi

# Split fastq files to as many files as compute nodes.
# Run them simultaneously and concatenate their resulting files.
l=\$(zcat $CBWADIR/temp.FASTQ\$1.fq.gz | wc -l) 
s1=\$((l / 4)) 
s2=\$((l % 4)) 
if [ \$s2 -ne 0 ]; then
  echo "Length of FASTQ files must be multiples of 4"
  exit
fi
s3=\$((s1 / $NUMBERCPU + $NUMBERCPU)) 
s4=\$((s3 * 4)) 
zcat $CBWADIR/temp.FASTQ\$1.fq.gz | split -d -a 2 -l \$s4 - $CBWADIR/temp.FASTQ\$1.

y=\$(($NUMBERCPU - 1)) 
for i in \$(eval echo {0..\$y}); do
  FASTQSPLITNUM=$CBWADIR/temp.FASTQ\$1.\$(printf "%02d" \$i) 
  mv \$FASTQSPLITNUM \$FASTQSPLITNUM.fq
  gzip \$FASTQSPLITNUM.fq &
done
wait

for i in \$(eval echo {0..\$y}); do
  FASTQSPLITNUM=$CBWADIR/temp.FASTQ\$1.\$(printf "%02d" \$i) 
  FASTQCUTADAPTNUM=$CBWADIR/FASTQ\$1.cutadapt.\$(printf "%02d" \$i).fq.gz

  $PYTHON $CUTADAPT --minimum-length=25 \\
    -a \${!ADAPTERSEQUENCE} \\
    -o \$FASTQCUTADAPTNUM \\
    \$FASTQSPLITNUM.fq.gz &> /dev/null &

done
wait

for i in \$(eval echo {0..\$y}); do
  FASTQSPLITNUM=$CBWADIR/temp.FASTQ\$1.\$(printf "%02d" \$i) 
  FASTQCUTADAPTNUM=$CBWADIR/FASTQ\$1.cutadapt.\$(printf "%02d" \$i).fq.gz
  FASTQPRINSEQNUM=$CBWADIR/FASTQ\$1.prinseq.\$(printf "%02d" \$i).fq.gz

  gzip -dc \$FASTQCUTADAPTNUM | \\
    perl pl/prinseq-lite.pl \\
      -no_qual_header \\
      -ns_max_n 0 \$PHRED64 \\
      -fastq stdin \\
      -trim_qual_right $MINTRIMQ \\
      -out_good stdout | \\
    gzip > \$FASTQPRINSEQNUM &
done
wait

CONFASTQPRINSEQNUM=""
for i in \$(eval echo {0..\$y}); do
  FASTQPRINSEQNUM=$CBWADIR/FASTQ\$1.prinseq.\$(printf "%02d" \$i).fq.gz
  CONFASTQPRINSEQNUM="\$CONFASTQPRINSEQNUM \$FASTQPRINSEQNUM"
done
cat \$CONFASTQPRINSEQNUM > $CBWADIR/outfile

# Single cpu version
#$PYTHON $CUTADAPT --minimum-length=25 \\
#  -a \${!ADAPTERSEQUENCE} \\
#  -o $CBWADIR/FASTQ\$1.cutadapt.fq.gz \\
#  $CBWADIR/temp.FASTQ\$1.fq.gz &> /dev/null
#
#gzip -dc $CBWADIR/FASTQ\$1.cutadapt.fq.gz | \\
#  perl pl/prinseq-lite.pl \\
#    -ns_max_n 0 \$PHRED64 \\
#    -fastq stdin \\
#    -trim_qual_right $MINTRIMQ \\
#    -out_good stdout | \\
#  gzip > $CBWADIR/outfile

#rm $CBWADIR/FASTQ\$1.cutadapt.fq.gz
mv $CBWADIR/outfile \$3
EOF

# Align
# $1: a three-digit number
# $2: a QC'ed gzipped fastq file
# $3: a sorted bam file from bwa
# $4: a sorted bam file from subread
grep ^QUALITYSCORE $SPECIESFILE | sed s/:/=/ > $BASEDIR/job-bwa-align
cat>>$BASEDIR/job-bwa-align<<EOF
QUALITYSCORESEQUENCE=QUALITYSCORE\$1

# We need more preprocessing of fastq files.
if [ "\${!QUALITYSCORESEQUENCE}" == "illumina" ]; then
  PHRED64=-I
else
  PHRED64=
fi

./bwa index -p $CBWADIR/$GENOMEFASTA-bwa -a is \\
  $CDATADIR/$GENOMEFASTA &> /dev/null

# BWA (in CRAM step) converts Illumina 1.5+ quality scores to Sanger scores.
# We no longer need to call bwa with -I option.
./bwa aln -t $NUMBERCPU \\
  $BWAOPTION \\
  $CBWADIR/$GENOMEFASTA-bwa \\
  \$2 > $CBWADIR/FASTQ\$1.sai 2> /dev/null

./bwa samse -n 1 \\
  -r "@RG\\tID:$PROJECTNAME\\tSM:BWA" \\
  $CBWADIR/$GENOMEFASTA-bwa \\
  $CBWADIR/FASTQ\$1.sai \\
  \$2 \\
  | ./samtools view -Sb -q $MINMAPQ - > $CBWADIR/FASTQ\$1.bam 
#  | ./samtools view -Sb -q 0 - > $CBWADIR/FASTQ\$1.bam 

./samtools sort $CBWADIR/FASTQ\$1.bam \$3 &> /dev/null
rm $CBWADIR/FASTQ\$1.sai
rm $CBWADIR/FASTQ\$1.bam
# FIXME:subread
exit

./subread-buildindex -o \\
  $CSUBREADDIR/$GENOMEFASTA-subread \\
  $CDATADIR/$GENOMEFASTA &> /dev/null

# Subread alignment
gunzip \$2
FASTAQFILE=\${2%.gz}
./subread-align \\
  --threads $NUMBERCPU \\
  --phred 3 \\
  --unique \\
  -i $CSUBREADDIR/$GENOMEFASTA-subread \\
  -r \$FASTAQFILE \\
  -o $CSUBREADDIR/FASTQ\$1.sam &> /dev/null

./samtools view -bS -o $CSUBREADDIR/FASTQ\$1.bam \\
  -q 0 \\
  $CSUBREADDIR/FASTQ\$1.sam &> /dev/null

./samtools sort $CSUBREADDIR/FASTQ\$1.bam \$4 &> /dev/null

rm $CSUBREADDIR/FASTQ\$1.sam
rm $CSUBREADDIR/FASTQ\$1.bam
rm \$FASTAQFILE
EOF

# Pileup
# $1: a three-digit number
# $2: a sorted bam file
# $3: a pileup file
cat>$BASEDIR/job-samtools-pileup<<EOF

# ./samtools mpileup -q 15 -d $READDEPTH \\

./samtools faidx $CDATADIR/$GENOMEFASTA &> /dev/null

./samtools mpileup \\
  -B \\
  -f $CDATADIR/$GENOMEFASTA \\
  \$2 \\
  > \$3

EOF


# DE count
cat>$BASEDIR/job-de.R<<EOF
library(DESeq)
library(ShortRead)
library(rtracklayer)
library(GenomicRanges)
library(VariantAnnotation)
library(GenomicFeatures)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3)
{
  cat ("Rscript job-de.R bam.file 1000 1000\n")
  quit("yes")
}
txdb.file <- "$RDATADIR/$REFGENOMETXDBBASE"
bam.file <- paste(args[1])
cl.file <- paste(args[1], "cl", sep=".")

# saveFeatures and loadFeatures
txdb <- loadFeatures(txdb.file)

# 1. Total number of short reads in fastq: zcat FASTQ051.fq.gz | wc -l
num.total.read <- args[2]
# 1.2 Total number of reads that pass filter :N:
num.n.read <- as.numeric(args[3])
# 1.3 Show reads mapped in multiple places, then where?
# ?
# 2. Number of short reads mapped 
indexBam(bam.file)
#bam.what <- c("qname","flag","pos","mapq","cigar")
#bvAll <- scanBam(bam.file, param=ScanBamParam(what=bam.what))
bv <- readBamGappedAlignments(bam.file,use.names=TRUE,param=ScanBamParam(what=c("mapq")))
# > length(bv)
num.mapped.read <- length(bv)
# 3. Number of uniquely mapped short reads
# Select mapped reads with mapq greater than or equal to $MINMAPQ
bv.multiple <- bv[elementMetadata(bv)["mapq"][,1] < $MINMAPQ]
bv <- bv[elementMetadata(bv)["mapq"][,1] >= $MINMAPQ]
num.unique.read <- length(bv)

# This takes much time.
# Find overlapped genes
feature.tx <- transcripts(txdb)
for (i in seq(length(feature.tx))) {
  x <- findOverlaps(feature.tx[i],feature.tx)
  stopifnot(length(x) > 0) 
  if (length(x) > 1) {
    x <- findOverlaps(feature.tx[i],feature.tx,minoverlap=100L)
    if (length(x) > 1) {
      print(i)
      print(x)
    }
  }
}
EOF

cat $BASEDIR/feature-txnc.txt >> $BASEDIR/job-de.R 
cat>>$BASEDIR/job-de.R<<EOF
# Remove strands from both feature and BAM alignments.
strand(feature.txnc) <- '*'

strand(bv) <- '*'
olap <- summarizeOverlaps(feature.txnc,bv,mode="IntersectionStrict")
cl <- assays(olap)\$counts[,1]
cl.nc <- cl[elementMetadata(feature.txnc)\$type=="NG"]

strand(bv.multiple) <- '*'
olap <- summarizeOverlaps(feature.txnc,bv.multiple,mode="IntersectionStrict")
cl.multiple <- assays(olap)\$counts[,1]
cl.multiple.nc <- cl.multiple[elementMetadata(feature.txnc)\$type=="NG"]
cl.multiple.nocds <- cl.multiple[elementMetadata(feature.txnc)\$type=="NOCDS"]

save(num.total.read, num.n.read, num.mapped.read, num.unique.read, cl, cl.nc, cl.multiple, cl.multiple.nc, cl.multiple.nocds, file=cl.file)

EOF

cat>$BASEDIR/job-de<<EOF
RSCRIPT=$CACRSCRIPT
\$RSCRIPT job-de.R \$1 \$2 \$3
EOF

cat>$BASEDIR/job-de-sum<<EOF
#!/bin/bash
RSCRIPT=$CACRSCRIPT
PREFIX1=BWA
PREFIX2=SUBREAD

RALIGNDIR1=$RBWADIR
RALIGNDIR2=$RSUBREADDIR
FILE=/tmp/\$(basename \$0).\$RANDOM.txt
echo $FASTQFILES > \$FILE
for k in {1..1}; do
  RALIGNDIR=RALIGNDIR\$k
  \$RSCRIPT job-de-sum.R \${!RALIGNDIR} \$FILE > \${!RALIGNDIR}/stat\$k.tex
done
rm \$FILE
EOF

cat>$BASEDIR/job-de-sum.R<<EOF
library(rtracklayer)
library(GenomicFeatures)
library(xtable)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2)
{
  cat ("Rscript job-de-sum.R output/ua159/1/bwa fastq_index_file\n")
  quit("yes")
}
txdb.file <- "$RDATADIR/$REFGENOMETXDBBASE"
txdb <- loadFeatures(txdb.file)
EOF
cat $BASEDIR/feature-txnc.txt >> $BASEDIR/job-de-sum.R 
cat>>$BASEDIR/job-de-sum.R<<EOF
count.table <- data.frame(gene=elementMetadata(feature.txnc)\$tx_name)
count.table.multiple <- data.frame(gene=elementMetadata(feature.txnc)\$tx_name)

#feature.tx <- transcripts(txdb)
#feature.cds <- cds(txdb, columns="exon_name")
#count.table <- data.frame(gene=unlist(elementMetadata(feature.cds)["exon_name"][,1]))

table.read.statistics <- c()
fastQIndex <- scan(args[2])
for (i in fastQIndex) {
  cl.file <- sprintf("%s/FASTQ%03d.sorted.bam.cl", args[1], i)
  load(cl.file)
  num.total.read <- as.integer(num.total.read)
  table.read.statistics <- 
    rbind(table.read.statistics, 
    c(i, num.total.read, sprintf("%d (%d)",num.n.read, round(num.n.read/num.total.read*100)), 
                         sprintf("%d (%d)",num.mapped.read, round(num.mapped.read/num.total.read*100)),
                         sprintf("%d (%d)",sum(cl.multiple.nocds), round(sum(cl.multiple.nocds)/num.total.read*100)),
                         sprintf("%d (%d)",num.unique.read, round(num.unique.read/num.total.read*100)),
			 sprintf("%d (%d)",sum(cl), round(sum(cl)/num.total.read*100)), 
			 sprintf("%d (%d)",sum(cl.nc), round(sum(cl.nc)/num.total.read*100))))
  count.table <- data.frame(count.table,cl)
  colnames(count.table)[ncol(count.table)] <- paste("X",i,sep="")
  count.table.multiple <- data.frame(count.table.multiple,cl.multiple)
  colnames(count.table.multiple)[ncol(count.table.multiple)] <- paste("X",i,sep="")
}
colnames(count.table) <- sub("X","",colnames(count.table))
count.table.file <- sprintf("%s/count.txt", args[1])
write.table(count.table,file=count.table.file,quote=FALSE,sep="\\t",row.names=FALSE)

colnames(count.table.multiple) <- sub("X","",colnames(count.table.multiple))
count.table.file <- sprintf("%s/count.multiple.txt", args[1])
write.table(count.table.multiple,file=count.table.file,quote=FALSE,sep="\\t",row.names=FALSE)

# Total reads: Total number of reads in the FASTQ raw file.
# N pass reads: Number of reads that pass the Illumina :N: filter
# Mapped reads: Number of reads that mapped on the reference genome
# Unique reads: Number of reads that are uniquely mapped
# Reads within regions: Number of reads that are uniquely mapped within regions
# Reads within nongenic regions: Number of reads that are uniquely mapped within
#                                nongenic regions
# Reads within non-CDS: Number of reads that are multiply mapped within non-CDS
colnames(table.read.statistics) <- c("Sample ID", "Total reads", "N pass reads", "Mapped reads", "Reads within non-CDS", "Unique reads", "Reads within regions", "Reads within nongenic regions")
x.big <- xtable( as.data.frame(table.read.statistics),
                 display=c("d","d","s","s","s","s","s","s","s"),
                 align=c('l','r','r','r','r','r','r','r','r'),
                 label='count',
                 caption='{\\\\bf Summary statistics of short reads.}'
               )
print( x.big,
       caption.placement="top",
       include.rownames=FALSE )
EOF
}

function batch2-get-data {
  GENOMEGFF=$(basename $REFGENOMEGFF)
cat>$BASEDIR/get-data.sh<<EOF
#!/bin/bash
#scp $CAC_USERHOST:$RBWADIR/*.wig $BWADIR
#scp $CAC_USERHOST:$RBWADIR/*.bed2 $BWADIR
#scp $CAC_USERHOST:$RBWADIR/*.operon $BWADIR
# scp $CAC_USERHOST:$RBWADIR/*.sorted.bam $BWADIR

# scp $CAC_USERHOST:$RBWADIR/*rrna $BWADIR
# scp $CAC_USERHOST:$RBWADIR/*-sum.pos $BWADIR
# scp $CAC_USERHOST:$RBWADIR/*.wig $BWADIR
# scp $CAC_USERHOST:$RBWADIR/intergeniconly.maf $BWADIR
# rm -f $BWADIR/rrna.tex
#for g in $FASTQFILES; do
#  FASTQNUM=FASTQ\$(printf "%03d" \$g)
  # perl pl/bwa-summary.pl rrnaToTex -rrna $BWADIR/\$FASTQNUM-sum.rrna >> $BWADIR/rrna.tex
  #scp $CAC_USERHOST:$RBWADIR/\$FASTQNUM.cutadapt.fq.gz $BWADIR
  #scp $CAC_USERHOST:$RBWADIR/\$FASTQNUM.prinseq.fq.gz $BWADIR
#  scp $CAC_USERHOST:$RSUBREADDIR/\$FASTQNUM.sorted.bam $SUBREADDIR
#done

scp $CAC_USERHOST:$RBWADIR/count.txt $BWADIR/count-$REFGENOMEID.txt
scp $CAC_USERHOST:$RBWADIR/count.multiple.txt $BWADIR/count-$REFGENOMEID.multiple.txt
scp $CAC_USERHOST:$RBWADIR/*.qualPlot.pdf $BWADIR
scp $CAC_USERHOST:$RBWADIR/stat1.tex $BWADIR
scp $CAC_USERHOST:$RBWADIR/tux.in $BWADIR
# scp $CAC_USERHOST:$RSUBREADDIR/count.txt $SUBREADDIR/count-$REFGENOMEID.txt
# scp $CAC_USERHOST:$RDATADIR/*.cram $DATADIR
# scp $CAC_USERHOST:$RDATADIR/$GENOMEGFF.fa $DATADIR

EOF
}

# UCSC
# 1. We need to know directories: local_tracks and /gbdb.
# 2. We might need to run hgsql.
# 3. We might need to place files at /gbdb.
# 4. We need to run make cornell at local_tracks.
# 
# We first create a database for a genome browser: job-ucsc-db. Copy this to a
# linux machine where UCSC genome browser is up and running. Just run it.
# We could have a run directory as we do in CAC. We need a GenBank file.
function ucsc-data {

cat>$BASEDIR/run-ucsc.sh<<EOF
bash job-ucsc-db
bash job-ucsc-bed
hgsql hgcentral_protected < dbDbInsert.sql
echo make cornell DBS=$DBNAME
# samtools faidx NC_004350.fna
# samtools view FASTQ001.sorted.bam | sed s/gi\|347750429\|ref\|NC_004350.2\|/chr1/ | samtools view -bS -T NC_004350.fna - > FASTQ001rnaseqBam.bam
# samtools index FASTQ001rnaseqBam.bam 
# mkdir /gbdb_cornell/SmuUA159v2/bbi
# mv FASTQ001rnaseqBam.bam* /gbdb_cornell/SmuUA159v2/bbi
# hgBbiDbLink SmuUA159v2 FASTQ001rnaseqBam /gbdb_cornell/SmuUA159v2/bbi/FASTQ001rnaseqBam.bam
EOF

cat>$BASEDIR/job-ucsc-db<<EOF
KENT=/usr/local/software/kent
gbToFaRa /dev/null $DBNAME.fna $DBNAME.ra $DBNAME.ta $DBNAME.gbk
toUpper $DBNAME.fna $DBNAME.fna.upper
faSize $DBNAME.fna.upper

rm $DBNAME.fna $DBNAME.ra $DBNAME.ta
echo ">chr1" > $DBNAME.fna
grep -v ">" $DBNAME.fna.upper >> $DBNAME.fna
rm $DBNAME.fna.upper

hgFakeAgp -minContigGap=1 $DBNAME.fna $DBNAME.agp
faToTwoBit $DBNAME.fna $DBNAME.2bit
mkdir -p /gbdb_cornell/$DBNAME/html
cp $DBNAME.2bit /gbdb_cornell/$DBNAME

echo "  creating a database ..."
twoBitInfo $DBNAME.2bit stdout | sort -k2nr > chrom.sizes
rm -rf bed
mkdir -p bed/chromInfo
awk '{printf "%s\\t%d\\t/gbdb_cornell/DBNAME/DBNAME.2bit\\n", \$1, \$2}' \\
  chrom.sizes > bed/chromInfo/chromInfo.tab.tmp
sed s/DBNAME/$DBNAME/g < bed/chromInfo/chromInfo.tab.tmp > bed/chromInfo/chromInfo.tab
hgsql -e "create database $DBNAME;" mysql

echo "  creating grp, chromInfo tables ..."
hgsql $DBNAME < \$KENT/src/hg/lib/grp.sql

cp bed/chromInfo/chromInfo.tab /tmp/
hgLoadSqlTab $DBNAME chromInfo \$KENT/src/hg/lib/chromInfo.sql \\
  /tmp/chromInfo.tab
rm /tmp/chromInfo.tab
hgGoldGapGl $DBNAME $DBNAME.agp

echo "  creating GC5 track ..."
mkdir bed/gc5Base
hgGcPercent -wigOut -doGaps -file=stdout -win=5 -verbose=0 $DBNAME \\
  $DBNAME.2bit | wigEncode stdin bed/gc5Base/gc5Base.{wig,wib}
hgLoadWiggle -pathPrefix=/gbdb_cornell/$DBNAME/wib \\
  $DBNAME gc5Base bed/gc5Base/gc5Base.wig
mkdir -p /gbdb_cornell/$DBNAME/wib/bed/gc5Base
cp bed/gc5Base/gc5Base.wib /gbdb_cornell/$DBNAME/wib/bed/gc5Base
rm -rf bed
done
EOF

cat>$BASEDIR/job-ucsc-bed<<EOF
KENT=/usr/local/software/kent
hgLoadBed $DBNAME knownGenes $DBNAME.knownGenes.bed

FASTQFILES=( $FASTQFILES )
for g in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$g)
  # BED tracks
  hgLoadBed $DBNAME \${FASTQNUM}rnaseqTx \${FASTQNUM}rnaseqTx.bed

  # Bam tracks
  hgBbiDbLink $DBNAME \${FASTQNUM}rnaseqBam \\
    /gbdb_cornell/$DBNAME/bbi/\${FASTQNUM}rnaseqBam .bam

  # Wiggle tracks
  wigEncode \$FASTQNUM.wig \$FASTQNUM.temp.wig \$FASTQNUM.wib
  hgLoadWiggle $DBNAME \${FASTQNUM}rnaseqWiggle \$FASTQNUM.temp.wig
  rm \$FASTQNUM.temp.wig
  mkdir /gbdb_cornell/$DBNAME/wib/
  mv \$FASTQNUM.wib /gbdb_cornell/$DBNAME/wib/
  hgsql $DBNAME -e "update \${FASTQNUM}rnaseqWiggle set file='/gbdb_cornell/$DBNAME/wib/\$FASTQNUM.wib'"

done

EOF

cat>$BASEDIR/dbDbInsert.sql<<EOF
INSERT INTO dbDb
    (name, description, nibPath, organism,
     defaultPos, active, orderKey, genome, scientificName,
     htmlPath, hgNearOk, hgPbOk, sourceName)
VALUES
    ("$DBNAME", "$ASSEMBLYNAME", "/gbdb_cornell/$DBNAME", "$GENOMENAME",
     "chr1:1-100000", 1, 1, "$GENOMENAME", "$GENOMENAME",
     "/gbdb_cornell/$DBNAME/html/description.html", 0, 0, "NCBI");
INSERT INTO defaultDb (genome, name) VALUES ("$GENOMENAME", "$DBNAME");
INSERT INTO genomeClade (genome, clade, priority) VALUES ("$GENOMENAME", "$CLADENAME", 9);
EOF

cat>$BASEDIR/send-ucsc-data.sh<<EOF
#!/bin/bash
FASTQFILES=( $FASTQFILES )
ssh -x $X11_USERNAME@$X11_LOGIN mkdir -p $CACWORKDIR

scp $BASEDIR/ucsc-trackDb.sh $X11_USERNAME@$X11_LOGIN:$CACWORKDIR
scp $BASEDIR/run-ucsc.sh $X11_USERNAME@$X11_LOGIN:$CACWORKDIR
scp $BASEDIR/job-ucsc-db $X11_USERNAME@$X11_LOGIN:$CACWORKDIR
scp $BASEDIR/job-ucsc-bed $X11_USERNAME@$X11_LOGIN:$CACWORKDIR
scp $BASEDIR/dbDbInsert.sql $X11_USERNAME@$X11_LOGIN:$CACWORKDIR
scp $BASEDIR/run-ucsc.sh $X11_USERNAME@$X11_LOGIN:$CACWORKDIR
scp $BASEDIR/1/data/feature-genome.out-geneonly $X11_USERNAME@$X11_LOGIN:$CACWORKDIR/$DBNAME.knownGenes.bed

for g in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$g)
  scp $BASEDIR/1/bwa/\$FASTQNUM.operon \\
    $X11_USERNAME@$X11_LOGIN:$CACWORKDIR/\${FASTQNUM}rnaseqTx.bed
  scp $BASEDIR/1/bwa/\$FASTQNUM.sorted.bam \\
    $X11_USERNAME@$X11_LOGIN:$CACWORKDIR
  scp $BASEDIR/1/bwa/\$FASTQNUM.wig \\
    $X11_USERNAME@$X11_LOGIN:$CACWORKDIR
done

scp $REFGENOMEGENBANK $X11_USERNAME@$X11_LOGIN:$CACWORKDIR/$DBNAME.gbk
# scp -qr $BWADIR $X11_USERNAME@$X11_LOGIN:public_html/rnaseq/bwa-$SPECIES
EOF

cat>$BASEDIR/ucsc-trackDb.sh<<EOF
track knownGenes
shortLabel Known Genes
longLabel Known Genes
group genes
priority 50
visibility pack
type bed 6
colorByStrand 255,0,0 0,0,255
itemRgb on
EOF
# Add more RNA-seq transcripts tracks
for g in $FASTQFILES; do
  FASTQNUM=FASTQ$(printf "%03d" $g)
cat>>$BASEDIR/ucsc-trackDb.sh<<EOF 

track ${FASTQNUM}rnaseqTx
shortLabel RNA-seq Transcripts ${FASTQNUM}
longLabel RNA-seq Transcripts ${FASTQNUM}
group genes
priority $g
visibility pack
type bed 6
colorByStrand 255,0,0 0,0,255
itemRgb on
EOF
done
# Add more BAM tracks
for g in $FASTQFILES; do
  FASTQNUM=FASTQ$(printf "%03d" $g)
cat>>$BASEDIR/ucsc-trackDb.sh<<EOF 

track ${FASTQNUM}rnaseqBam
shortLabel RNA-seq BAM ${FASTQNUM}
longLabel RNA-seq BAM ${FASTQNUM}
group genes
priority $g
visibility pack
type bam
EOF
done
# Add more Wiggle tracks
for g in $FASTQFILES; do
  FASTQNUM=FASTQ$(printf "%03d" $g)
cat>>$BASEDIR/ucsc-trackDb.sh<<EOF 

track ${FASTQNUM}rnaseqWiggle
shortLabel RNA-seq Wiggle ${FASTQNUM}
longLabel RNA-seq Wiggle ${FASTQNUM}
group genes
priority $g
visibility full
type wig
autoScale on
alwaysZero on
color=0,200,100 
EOF
done

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
# FIXME:public_html
  scp -q $BASEDIR/ucsc-data.sh $X11_USERNAME@$X11_LOGIN:public_html/rnaseq
  echo "To send the data to the genome browser:"
  echo "bash $BASEDIR/send-ucsc-data.sh"
}

function batch-run-parsernaseq {
  GENOMEFASTA=$(basename $REFGENOMEFASTA)
  GENOMEGFF=$(basename $REFGENOMEGFF)
  GENOMEPTT=$(basename $REFGENOMEPTT)
  STATUS=parsernaseq
cat>$BASEDIR/run-$STATUS.sh<<EOF
#!/bin/bash
FASTQFILES=( $FASTQFILES )
c=\$((\${#FASTQFILES[@]} / $NUMBERCPU))
b=\$((\${#FASTQFILES[@]} % $NUMBERCPU))
if [ \$b -ne 0 ]; then
  c=\$((c+1))
fi
sed s/PBSARRAYSIZE/\$c/g < batch-$STATUS.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF

cat>$BASEDIR/batch-$STATUS.sh<<EOF
#!/bin/bash
#PBS -l walltime=24:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-PARSE
#PBS -q ${QUEUENAME}
#PBS -m e
$EMAILON#PBS -M ${BATCHEMAIL}
#PBS -t 1-PBSARRAYSIZE

function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$PBS_O_WORKDIR/job-$STATUS .
  cp \$HOME/$PARSERNASEQ ParseRNAseq

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CBWADIR
  cp $RBWADIR/feature-genome.out-geneonly $CBWADIR
  cp $RDATADIR/$GENOMEFASTA $CDATADIR
  cp $RDATADIR/$GENOMEGFF $CDATADIR
  cp $RDATADIR/$GENOMEPTT $CDATADIR
}

#function retrieve-data {
  # cp $CBWADIR/*parsernaseq* $RBWADIR
  # cp $CBWADIR/*.bed* $RBWADIR
  # cp $CBWADIR/*.operon* $RBWADIR
#}

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
cd
rm -rf \$TMPDIR
EOF

cat>$BASEDIR/job-$STATUS<<EOF
FASTQNUM=FASTQ\$1
cp $RBWADIR/\$FASTQNUM.wig $CBWADIR
# cp $RBWADIR/\$FASTQNUM-sum.pos $CBWADIR

# in BED format
perl pl/feature-genome.pl ptt2 \\
  -geneonly \\
  -in $CDATADIR/$GENOMEPTT \\
  -out $CBWADIR/feature-genome.out-geneonly

# This is okay.
perl pl/transcript-parsernaseq.pl pileup \\
  -wiggle $CBWADIR/\$FASTQNUM.wig \\
  -out $CBWADIR/\$FASTQNUM.parsernaseq.pileup

# not a BED format, inclusive.
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
cp $CBWADIR/\$FASTQNUM.parsernaseq $RBWADIR

# This did not work.
#perl pl/transcript-parsernaseq.pl bed \\
#  -parsernaseq $CBWADIR/\$FASTQNUM.parsernaseq \\
#  -out $CBWADIR/\$FASTQNUM.bed

# To create a BED file from the prediction.
# This is okay.
perl pl/transcript-parsernaseq.pl gff \\
  -parsernaseq $CBWADIR/\$FASTQNUM.parsernaseq \\
  -out $CBWADIR/\$FASTQNUM.bed2
cp $CBWADIR/\$FASTQNUM.bed2 $RBWADIR

# To create operon with strands.
# This is okay.
perl pl/transcript-parsernaseq.pl operon \\
  -feature $CBWADIR/feature-genome.out-geneonly \\
  -parsernaseq $CBWADIR/\$FASTQNUM.parsernaseq \\
  -out $CBWADIR/\$FASTQNUM.operon
cp $CBWADIR/\$FASTQNUM.operon $RBWADIR

# This works.
#perl pl/bwa-pos2wig.pl end \\
#  -genomeLength $REFGENOMELENGTH \\
#  -in $CBWADIR/\$FASTQNUM-sum.pos \\
#  -out $CBWADIR/\$FASTQNUM-end.wig
#
## This works.
#perl pl/transcript-parsernaseq.pl adjust \\
#  -end $CBWADIR/\$FASTQNUM-end.wig \\
#  -operon $CBWADIR/\$FASTQNUM.operon \\
#  -out $CBWADIR/\$FASTQNUM.operon2
## This works.
#perl pl/transcript-parsernaseq.pl slope \\
#  -windowsize 95 \\
#  -end $CBWADIR/\$FASTQNUM-end.wig \\
#  -operon $CBWADIR/\$FASTQNUM.operon \\
#  -out $CBWADIR/\$FASTQNUM.operon.slope
EOF
}

####################################################################
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




function batch2-run-blast {
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
$EMAILON#PBS -M ${BATCHEMAIL}
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


###############################################################################
# job
###############################################################################
function batch2-make-job {
# 1: a sorted bam file
# 2: a fasta file for the reference genome
# 3: an output file
cat>$BASEDIR/job-tux-sum<<EOF
#!/bin/bash
rm $RBWADIR/tux.in
for g in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$g)
  echo -en "  - [" > 2.tmp
  cat $RBWADIR/\$FASTQNUM.bam.tux | tr "\\n" "," >> 2.tmp
  sed s/,$// 2.tmp > $RBWADIR/\$FASTQNUM.bam.tux.array
  echo -en "]\n" >> $RBWADIR/\$FASTQNUM.bam.tux.array
  rm 2.tmp
  cat $RBWADIR/\$FASTQNUM.bam.tux.array >> $RBWADIR/tux.in
done
echo Check $RBWADIR/tux.in
EOF
cat>$BASEDIR/job-tux<<EOF
RSCRIPT=$CACRSCRIPT
\$RSCRIPT job-tux.R \$1 \$2 \$3
EOF
cat>$BASEDIR/job-tux.R<<EOF
library(Rsamtools)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3)
{
  cat ("Rscript job-tux.R bam.file a.fasta\n")
  quit("yes")
}
bam.file <- args[1]
tux.file <- args[3]

# Read the length of a sequence, and create a vector of int of the size being
# the length.
ref_fn <-  args[2]
ref_f <- FaFile( ref_fn ) 
open.FaFile( ref_f )
ref_seq <- getSeq(ref_f)
count.read <- rep.int(x=0,times=width(ref_seq))

indexBam(bam.file)

bv <- readBamGappedAlignments(bam.file)
a <- table(c(start(bv[strand(bv)=='+']), end(bv[strand(bv)=='-'])))
# Print this count.read
# Test with a samle fasta and a sorted bam file
count.read[as.integer(names(a))] <- as.vector(a)
write(count.read,file=tux.file,ncolumns=1)
close.FaFile( ref_f )
EOF
}
