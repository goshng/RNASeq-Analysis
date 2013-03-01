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
function batch3 {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SPECIESS[@]}; do 
  if [ "$SPECIES" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  else  
    batch3-detect-os
    batch3-repetition
    batch3-variable
    batch3-output
    batch3-speciesfile 
    if [ "$RUNMODE" == "local" ]; then
      batch3-push-localdata
      batch3-get-localdata
      batch3-clean-localdata
    else
      batch3-push-data
      batch3-get-data
      batch3-clean-data
    fi
    batch3-make-job
    batch3-run-qcalignde
    create-index
    if [ "$RUNMODE" == "local" ]; then
      batch3-copy-localscripts
      batch3-localmessage
    else
      batch3-copy-scripts
      batch3-rmessage
    fi
    break
  fi
  done
}

function batch3-detect-os {
  if [ "$REMOTEMACHINE" == 'Darwin' ]; then
    ZCAT=gzcat
    SPLIT=/opt/local/bin/gsplit
  else
    ZCAT=zcat
    SPLIT=split
  fi 
}

function batch3-repetition {
  SPECIESFILE=$ROOTANALYSISDIR/species/$SPECIES
  REPETITION=$(grep ^REPETITION\: $SPECIESFILE | cut -d":" -f2)
  if [ -z "$REPETITION" ]; then
    REPETITION=1
  fi
}

# Set the shell variables for batch3 script
function batch3-variable {
  # Local output directories.
  BASEDIR=$ROOTANALYSISDIR/output/$SPECIES
  NUMBERDIR=$BASEDIR/$REPETITION
  DATADIR=$NUMBERDIR/data
  BWADIR=$NUMBERDIR/bwa
  RUNANALYSIS=$NUMBERDIR/run-analysis

  # Remote output directories.
  RBASEDIR=$ROUTPUTDIR/$SPECIES
  RNUMBERDIR=$ROUTPUTDIR/$SPECIES/$REPETITION
  RDATADIR=$RNUMBERDIR/data
  RBWADIR=$RNUMBERDIR/bwa
  RRUNANALYSIS=$RNUMBERDIR/run-analysis

  # Compute node output directories.
  CBASEDIR=output/$SPECIES
  CNUMBERDIR=$CBASEDIR/$REPETITION
  CDATADIR=$CNUMBERDIR/data
  CBWADIR=$CNUMBERDIR/bwa
  CRUNANALYSIS=$CNUMBERDIR/run-analysis
}

# Create output directories for the species directory
function batch3-output {
  mkdir -p $DATADIR
  mkdir -p $BWADIR
  mkdir -p $RUNANALYSIS
}

# Set species file variables
function batch3-speciesfile {
  FASTQFILES=$(grep ^FASTQFILES\: $SPECIESFILE | cut -d":" -f2)
  FASTQLABEL=$(grep ^FASTQLABEL\: $SPECIESFILE | cut -d":" -f2)
  CACWORKDIR=$(grep ^CACWORKDIR\: $SPECIESFILE | cut -d":" -f2)
  REFGENOMEFASTA=$(grep ^REFGENOMEFASTA\: $SPECIESFILE | cut -d":" -f2)
  REFGENOMEFASTA=$ROOTANALYSISDIR/$REFGENOMEFASTA
  GENOMEFASTA=$(basename $REFGENOMEFASTA)
  REFGENOMETXDB=$(grep ^REFGENOMETXDB\: $SPECIESFILE | cut -d":" -f2)  
  REFGENOMETXDB=$ROOTANALYSISDIR/$REFGENOMETXDB
  GENOMETXDB=$(basename $REFGENOMETXDB)
  REFGENOMERRNATXDB=$(grep ^REFGENOMERRNATXDB\: $SPECIESFILE | cut -d":" -f2)  
  REFGENOMERRNATXDB=$ROOTANALYSISDIR/$REFGENOMERRNATXDB
  GENOMERRNATXDB=$(basename $REFGENOMERRNATXDB)
  QCALIGNDEWALLTIME=$(grep ^QCALIGNDEWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  MINMAPQ=$(grep ^MINMAPQ\: $SPECIESFILE | cut -d":" -f2)
  MINTRIMQ=$(grep ^MINTRIMQ\: $SPECIESFILE | cut -d":" -f2)
  TRIMLEFT=$(grep ^TRIMLEFT\: $SPECIESFILE | cut -d":" -f2)
  KEEPBAM=$(grep ^KEEPBAM\: $SPECIESFILE | cut -d":" -f2)
  COUNTMODE=$(grep ^COUNTMODE\: $SPECIESFILE | cut -d":" -f2)
  RUNMODE=$(grep ^RUNMODE\: $SPECIESFILE | cut -d":" -f2)

  # Check whether the requried species file vairables exist.
  if [ -z "$FASTQLABEL" ]; then
    echo Add FASTQLABEL to the species file
    exit 1
  fi
  if [ -z "$FASTQFILES" ]; then
    echo Add FASTQFILES to the species file
    exit 1
  fi
  if [ -z "$REFGENOMEFASTA" ]; then
    echo No reference genome sequence in the file species/$SPECIES
    echo Add REFGENOMEFASTA to the species file
    exit 1
  fi
  if [ -z "$REFGENOMETXDB" ]; then
    echo No txdb in the file species/$SPECIES
    echo Add REFGENOMETXDB to the species file
    exit 1
  fi
  if [ -z "$MINTRIMQ" ]; then
    echo No minimum trim quality in the file species/$SPECIES
    echo Add MINTRIMQ to the species file
    exit 1
  fi
  if [ -z "$TRIMLEFT" ]; then
    TRIMLEFT=0
  fi
  if [ -z "$MINMAPQ" ]; then
    echo No minimum mapping quality in the file species/$SPECIES
    echo Add MINMAPQ to the species file
    exit 1
  fi
  if [ -z "$COUNTMODE" ]; then
    echo No count mode in the file species/$SPECIES
    echo Add COUNTMODE to the species file
    exit 1
  fi
  if [ -z "$RUNMODE" ]; then
    RUNMODE=local
  fi
  if [ "$RUNMODE" == "local" ]; then
    CACWORKDIR=""
  else
    if [ -z "$CACWORKDIR" ]; then
      echo No remote work directory in the file species/$SPECIES
      echo Add CACWORKDIR to the species file
      exit 1
    fi
  fi 
}

# Send input files to the remote machine.
function batch3-push-localdata {
cat>$BASEDIR/push-data.sh<<EOF
#!/bin/bash
mkdir -p $RBWADIR
mkdir -p $RDATADIR
mkdir -p $RRUNANALYSIS
ln -s $REFGENOMEFASTA $RDATADIR
ln -s $REFGENOMETXDB $RDATADIR
for g in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$g)
  FASTQFILE=$ROOTANALYSISDIR/\$(grep ^\$FASTQNUM\: $SPECIESFILE | cut -d":" -f2)
  if [ -f \$FASTQFILE ]; then
    ln -s \$FASTQFILE $RDATADIR/\$FASTQNUM.fq.gz
  else
    echo No such file: \$FASTQFILE
  fi
done
echo "Check $RDATADIR"
EOF
}

function batch3-push-data {
cat>$BASEDIR/push-data.sh<<EOF
#!/bin/bash
ssh -x $CAC_USERHOST mkdir -p $RBWADIR
ssh -x $CAC_USERHOST mkdir -p $RDATADIR
ssh -x $CAC_USERHOST mkdir -p $RRUNANALYSIS
scp $REFGENOMEFASTA $CAC_USERHOST:$RDATADIR
scp $REFGENOMETXDB $CAC_USERHOST:$RDATADIR
# scp $REFGENOMERRNATXDB $CAC_USERHOST:$RDATADIR
for g in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$g)
  FASTQFILE=$ROOTANALYSISDIR/\$(grep ^\$FASTQNUM\: $SPECIESFILE | cut -d":" -f2)
  if [ -f \$FASTQFILE ]; then
    scp \$FASTQFILE $CAC_USERHOST:$RDATADIR/\$FASTQNUM.fq.gz
  else
    echo No such file: \$FASTQFILE
  fi
done
echo "Check cac:$RDATADIR"
EOF
}

function batch3-get-localdata {
cat>$BASEDIR/get-data.sh<<EOF
#!/bin/bash
echo $GENOMETXDB > $BWADIR/txdb.txt
cp $RBWADIR/count $BWADIR
cp $RBWADIR/count.cds $BWADIR
cp $RBWADIR/count.nocds $BWADIR
cp $RBWADIR/count.ng $BWADIR
cp $RBWADIR/count.tx $BWADIR
cp $RBWADIR/count.m $BWADIR
cp $RBWADIR/count.m.cds $BWADIR
cp $RBWADIR/count.m.nocds $BWADIR
cp $RBWADIR/count.m.ng $BWADIR
cp $RBWADIR/count.m.tx $BWADIR
cp $RBWADIR/count.m.rrna $BWADIR

for i in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$i)
  cp $RBWADIR/\$FASTQNUM*fq.qualPlot.* $BWADIR
  if [ "$KEEPBAM" == "YES" ]; then
    cp $RBWADIR/\$FASTQNUM.sorted.bam* $BWADIR
  fi
done 
cp $RBWADIR/stat1.tex $BWADIR
echo "Check $BWADIR"
EOF
}

function batch3-get-data {
cat>$BASEDIR/get-data.sh<<EOF
#!/bin/bash
echo $GENOMETXDB > $BWADIR/txdb.txt
scp $CAC_USERHOST:$RBWADIR/count $BWADIR
scp $CAC_USERHOST:$RBWADIR/count.cds $BWADIR
scp $CAC_USERHOST:$RBWADIR/count.nocds $BWADIR
scp $CAC_USERHOST:$RBWADIR/count.ng $BWADIR
scp $CAC_USERHOST:$RBWADIR/count.tx $BWADIR
scp $CAC_USERHOST:$RBWADIR/count.m $BWADIR
scp $CAC_USERHOST:$RBWADIR/count.m.cds $BWADIR
scp $CAC_USERHOST:$RBWADIR/count.m.nocds $BWADIR
scp $CAC_USERHOST:$RBWADIR/count.m.ng $BWADIR
scp $CAC_USERHOST:$RBWADIR/count.m.tx $BWADIR
scp $CAC_USERHOST:$RBWADIR/count.m.rrna $BWADIR
for i in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$i)
  scp $CAC_USERHOST:$RBWADIR/\$FASTQNUM*fq.qualPlot.* $BWADIR
  if [ "$KEEPBAM" == "YES" ]; then
    scp $CAC_USERHOST:$RBWADIR/\$FASTQNUM.sorted.bam* $BWADIR
  fi
done 
scp $CAC_USERHOST:$RBWADIR/stat1.tex $BWADIR
echo "Check $BWADIR"
EOF
}

function batch3-clean-localdata {
cat>$BASEDIR/clean-data.sh<<EOF
#!/bin/bash
echo rm -rf $RBASEDIR
EOF
}

function batch3-clean-data {
cat>$BASEDIR/clean-data.sh<<EOF
#!/bin/bash
echo ssh -x $CAC_USERHOST rm -rf $RBASEDIR
EOF
}

################################################################################
# jobs
# 1. job-stat
# 2. job-fastqc
# 3. job-
# 4. job-
################################################################################
function batch3-make-job {

################################################################################
# Plot quality scores of a FASTQ file.
# $1: a three-digit number
# $2: a gzipped fastq file
# $3: RData
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

################################################################################
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
  $ZCAT \$2 \\
    | grep -A 3 '^@.* [^:]*:N:[^:]*:' \\
    | sed '/^--$/d' | gzip > $CBWADIR/temp.FASTQ\$1.fq.gz
fi

# Split fastq files to as many files as compute nodes.
# Run them simultaneously and concatenate their resulting files.
l=\$($ZCAT $CBWADIR/temp.FASTQ\$1.fq.gz | wc -l) 
s1=\$((l / 4)) 
s2=\$((l % 4)) 
if [ \$s2 -ne 0 ]; then
  echo "Length of FASTQ files must be multiples of 4"
  exit
fi
s3=\$((s1 / $NUMBERCPU + $NUMBERCPU)) 
s4=\$((s3 * 4)) 
$ZCAT $CBWADIR/temp.FASTQ\$1.fq.gz | $SPLIT -d -a 2 -l \$s4 - $CBWADIR/temp.FASTQ\$1.

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
      -trim_left $TRIMLEFT \\
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

################################################################################
# Align
# $1: a three-digit number
# $2: a QC'ed gzipped fastq file
# $3: a sorted bam file from bwa
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

./bwa aln \\
  -t $NUMBERCPU \\
  $PHRED64 \\
  $CBWADIR/$GENOMEFASTA-bwa \\
  \$2 > $CBWADIR/FASTQ\$1.sai 2> /dev/null

./bwa samse -n 1 \\
  -r "@RG\\tID:$PROJECTNAME\\tSM:BWA" \\
  $CBWADIR/$GENOMEFASTA-bwa \\
  $CBWADIR/FASTQ\$1.sai \\
  \$2 \\
  | ./samtools view -Sb -q 0 - > $CBWADIR/FASTQ\$1.bam 

./samtools sort $CBWADIR/FASTQ\$1.bam \$3 &> /dev/null
rm $CBWADIR/FASTQ\$1.sai
rm $CBWADIR/FASTQ\$1.bam

EOF

cat>$BASEDIR/feature-txnc.txt<<EOF
# 1 Reads mapped on CDS tx, feature.cds
# 2 Reads mapped on non CDS tx, feature.nocds
# 3 Reads mapped on non-genic tx, feature.ng
# 4 Reads mapped on tx and non-genic tx, feature.txnc
# KEY
# 1. Transcripts
feature.tx <- transcripts(txdb)
# 2. CDS
feature.cds <- cds(txdb,columns="tx_id")
stopifnot(length(feature.cds) == length(unlist(elementMetadata(feature.cds)\$tx_id)))
# 3. Grab transcripts that contain CDS. The x contains TRUE or FALSE that
# denotes whether a transcript contains a CDS or not.
x <- elementMetadata(feature.tx)\$tx_id %in% unlist(elementMetadata(feature.cds)\$tx_id)
# 4. Subset CDS-containing transcripts of all of the transcripts.
feature.cds <- feature.tx[x]
# 5. Add a metadata column, type="CDS".
elementMetadata(feature.cds) <- data.frame(tx_id=elementMetadata(feature.cds)\$tx_id, tx_name=elementMetadata(feature.cds)\$tx_name,type="CDS")

# 6. feature.txnc contains everything.
feature.txnc <- feature.cds

# 7. Subset transcripts without CDS of all of the transcripts.
feature.nocds <- feature.tx[!x]
# 8. Add a meta data column, type="NOCDS".
if (length(feature.nocds) > 0) {
  elementMetadata(feature.nocds) <- data.frame(tx_id=elementMetadata(feature.nocds)\$tx_id, tx_name=elementMetadata(feature.nocds)\$tx_name,type="NOCDS")
  feature.txnc <- c(feature.txnc,feature.nocds)
}

# 9. Reads mapped on non-genic tx
feature.ng <- transcripts(txdb)
strand(feature.ng) <- '*'
feature.ng <- gaps(reduce(feature.ng))
feature.ng <- feature.ng[strand(feature.ng)=='*']
# 10. Reads mapped on tx and non-genic tx
if (length(feature.ng) > 0) {
  x <- seq(from=length(feature.tx)+1,to=length(feature.tx)+length(feature.ng))
  y <- paste("NC",x,sep="_")
  elementMetadata(feature.ng) <- data.frame(tx_id=x,tx_name=y,type="NG")
  rm(x,y)
  feature.txnc <- c(feature.txnc,feature.ng)
}

# 11. rRNA
feature.rrna <- transcripts(rrna.txdb)
###################################################################
EOF

################################################################################
# DE
# $1: a three-digit number
# $2: the number of reads in FASTQ file
# $3: the number of reads in FASTQ file that is used as an input to the alignment
cat>$BASEDIR/job-de<<EOF
RSCRIPT=$CACRSCRIPT
\$RSCRIPT job-de.R \$1 \$2 \$3
EOF

cat>$BASEDIR/job-de.R<<EOF
library(DESeq)
library(ShortRead)
library(rtracklayer)
library(GenomicRanges)
library(VariantAnnotation)
library(GenomicFeatures)
library(Rsamtools)         # readBamGappedAlignments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3)
{
  cat ("Rscript job-de.R bam.file 1000 1000\n")
  quit("yes")
}
rrna.txdb.file <- "$CDATADIR/$GENOMERRNATXDB"
rrna.txdb <- loadDb(rrna.txdb.file)

txdb.file <- "$CDATADIR/$GENOMETXDB"
bam.file <- paste(args[1])
cl.file <- paste(args[1], "cl", sep=".")
txdb <- loadDb(txdb.file)
num.total.read <- as.numeric(args[2])
num.n.read <- as.numeric(args[3])
indexBam(bam.file)
bv <- readBamGappedAlignments(bam.file,
                              use.names=TRUE,
                              param=ScanBamParam(what=c("mapq")))
num.mapped.read <- length(bv)
bv.multiple <- bv[elementMetadata(bv)["mapq"][,1] < $MINMAPQ]
bv <- bv[elementMetadata(bv)["mapq"][,1] >= $MINMAPQ]
num.unique.read <- length(bv)

# Report overlapping genes.
feature.tx <- transcripts(txdb)
x <- countOverlaps(feature.tx,feature.tx,minoverlap=100L)
overlapping.tx <- elementMetadata(feature.tx)\$tx_name[seq(length(x))[x > 1]]
EOF
cat $BASEDIR/feature-txnc.txt >> $BASEDIR/job-de.R 
cat>>$BASEDIR/job-de.R<<EOF
# Remove strands from both feature and BAM alignments.
strand(feature.txnc) <- '*'
strand(bv) <- '*'
olap <- summarizeOverlaps(feature.txnc,bv,mode="$COUNTMODE")
cl <- assays(olap)\$counts[,1]
cl.cds <- cl[elementMetadata(feature.txnc)\$type=="CDS"]
cl.nocds <- cl[elementMetadata(feature.txnc)\$type=="NOCDS"]
cl.tx <- cl[!elementMetadata(feature.txnc)\$type=="NG"]
cl.ng <- cl[elementMetadata(feature.txnc)\$type=="NG"]

strand(bv.multiple) <- '*'
olap <- summarizeOverlaps(feature.txnc,bv.multiple,mode="$COUNTMODE")
cl.m <- assays(olap)\$counts[,1]
cl.m.cds <- cl.m[elementMetadata(feature.txnc)\$type=="CDS"]
cl.m.nocds <- cl.m[elementMetadata(feature.txnc)\$type=="NOCDS"]
cl.m.tx <- cl.m[!elementMetadata(feature.txnc)\$type=="NG"]
cl.m.ng <- cl.m[elementMetadata(feature.txnc)\$type=="NG"]

strand(feature.rrna) <- '*'
olap.rrna <- summarizeOverlaps(feature.rrna,bv.multiple,mode="$COUNTMODE")
cl.m.rrna <- assays(olap.rrna)\$counts[,1]

save(overlapping.tx, 
     num.total.read, num.n.read, num.mapped.read, num.unique.read, 
     cl, cl.cds, cl.nocds, cl.tx, cl.ng, 
     cl.m, cl.m.cds, cl.m.nocds, cl.m.tx, cl.m.ng, cl.m.rrna,
     file=cl.file)
EOF

cat>$BASEDIR/job-de-sum<<EOF
#!/bin/bash
RSCRIPT=$CACRSCRIPT
PREFIX1=BWA
RALIGNDIR1=$RBWADIR
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
rrna.txdb.file <- "$RDATADIR/$GENOMERRNATXDB"
rrna.txdb <- loadDb(rrna.txdb.file)

txdb.file <- "$RDATADIR/$GENOMETXDB"
txdb <- loadDb(txdb.file)
EOF
cat $BASEDIR/feature-txnc.txt >> $BASEDIR/job-de-sum.R 
cat>>$BASEDIR/job-de-sum.R<<EOF
# Create the first column.
count.table <- data.frame(gene=elementMetadata(feature.txnc)\$tx_name)
x <- elementMetadata(feature.txnc)\$type=="CDS"
count.table.cds <- data.frame(gene=elementMetadata(feature.txnc)\$tx_name[x])
x <- elementMetadata(feature.txnc)\$type=="NOCDS"
count.table.nocds <- data.frame(gene=elementMetadata(feature.txnc)\$tx_name[x])
x <- elementMetadata(feature.txnc)\$type=="NG"
count.table.tx <- data.frame(gene=elementMetadata(feature.txnc)\$tx_name[!x])
count.table.ng <- data.frame(gene=elementMetadata(feature.txnc)\$tx_name[x])

count.m.table <- count.table
count.m.table.cds <- count.table.cds
count.m.table.nocds <- count.table.nocds
count.m.table.tx <- count.table.tx
count.m.table.ng <- count.table.ng

count.m.table.rrna <- data.frame(gene=elementMetadata(feature.rrna)\$tx_name)

table.read.statistics <- c()
fastQIndex <- scan(args[2])
for (i in fastQIndex) {
  cl.file <- sprintf("%s/FASTQ%03d.sorted.bam.cl", args[1], i)
  load(cl.file)
  table.read.statistics <- 
    rbind(table.read.statistics, 
          c(i, num.total.read, 
            sprintf("%d (%d)",num.n.read, 
                    round(num.n.read/num.total.read*100)), 
            sprintf("%d (%d)",num.mapped.read, 
                    round(num.mapped.read/num.total.read*100)),
            sprintf("%d (%d)",sum(cl.m.nocds), 
                    round(sum(cl.m.nocds)/num.total.read*100)),
            sprintf("%d (%d)",sum(cl.m.rrna), 
                    round(sum(cl.m.rrna)/num.total.read*100)),
            sprintf("%d (%d)",num.unique.read, 
                    round(num.unique.read/num.total.read*100)),
            sprintf("%d (%d)",sum(cl), 
                    round(sum(cl)/num.total.read*100)), 
            sprintf("%d (%d)",sum(cl.ng), 
                    round(sum(cl.ng)/num.total.read*100)))
         )
  count.table <- data.frame(count.table,cl)
  colnames(count.table)[ncol(count.table)] <- paste("X",i,sep="")
  count.table.cds <- data.frame(count.table.cds,cl.cds)
  colnames(count.table.cds)[ncol(count.table.cds)] <- paste("X",i,sep="")
  count.table.nocds <- data.frame(count.table.nocds,cl.nocds)
  colnames(count.table.nocds)[ncol(count.table.nocds)] <- paste("X",i,sep="")
  count.table.tx <- data.frame(count.table.tx,cl.tx)
  colnames(count.table.tx)[ncol(count.table.tx)] <- paste("X",i,sep="")
  count.table.ng <- data.frame(count.table.ng,cl.ng)
  colnames(count.table.ng)[ncol(count.table.ng)] <- paste("X",i,sep="")

  count.m.table <- data.frame(count.m.table,cl.m)
  colnames(count.m.table)[ncol(count.m.table)] <- paste("X",i,sep="")
  count.m.table.cds <- data.frame(count.m.table.cds,cl.m.cds)
  colnames(count.m.table.cds)[ncol(count.m.table.cds)] <- paste("X",i,sep="")
  count.m.table.nocds <- data.frame(count.m.table.nocds,cl.m.nocds)
  colnames(count.m.table.nocds)[ncol(count.m.table.nocds)] <- paste("X",i,sep="")
  count.m.table.rrna <- data.frame(count.m.table.rrna,cl.m.rrna)
  colnames(count.m.table.rrna)[ncol(count.m.table.rrna)] <- paste("X",i,sep="")
  count.m.table.tx <- data.frame(count.m.table.tx,cl.m.tx)
  colnames(count.m.table.tx)[ncol(count.m.table.tx)] <- paste("X",i,sep="")
  count.m.table.ng <- data.frame(count.m.table.ng,cl.m.ng)
  colnames(count.m.table.ng)[ncol(count.m.table.ng)] <- paste("X",i,sep="")
}
colnames(count.table) <- sub("X","",colnames(count.table))
count.table.file <- sprintf("%s/count", args[1])
write.table(count.table,file=count.table.file,quote=FALSE,sep="\\t",row.names=FALSE)
colnames(count.table.cds) <- sub("X","",colnames(count.table.cds))
count.table.file <- sprintf("%s/count.cds", args[1])
write.table(count.table.cds,file=count.table.file,quote=FALSE,sep="\\t",row.names=FALSE)
colnames(count.table.nocds) <- sub("X","",colnames(count.table.nocds))
count.table.file <- sprintf("%s/count.nocds", args[1])
write.table(count.table.nocds,file=count.table.file,quote=FALSE,sep="\\t",row.names=FALSE)
colnames(count.table.tx) <- sub("X","",colnames(count.table.tx))
count.table.file <- sprintf("%s/count.tx", args[1])
write.table(count.table.tx,file=count.table.file,quote=FALSE,sep="\\t",row.names=FALSE)
colnames(count.table.ng) <- sub("X","",colnames(count.table.ng))
count.table.file <- sprintf("%s/count.ng", args[1])
write.table(count.table.ng,file=count.table.file,quote=FALSE,sep="\\t",row.names=FALSE)

colnames(count.m.table) <- sub("X","",colnames(count.m.table))
count.m.table.file <- sprintf("%s/count.m", args[1])
write.table(count.m.table,file=count.m.table.file,quote=FALSE,sep="\\t",row.names=FALSE)
colnames(count.m.table.cds) <- sub("X","",colnames(count.m.table.cds))
count.m.table.file <- sprintf("%s/count.m.cds", args[1])
write.table(count.m.table.cds,file=count.m.table.file,quote=FALSE,sep="\\t",row.names=FALSE)
colnames(count.m.table.nocds) <- sub("X","",colnames(count.m.table.nocds))
count.m.table.file <- sprintf("%s/count.m.nocds", args[1])
write.table(count.m.table.nocds,file=count.m.table.file,quote=FALSE,sep="\\t",row.names=FALSE)
colnames(count.m.table.rrna) <- sub("X","",colnames(count.m.table.rrna))
count.m.table.file <- sprintf("%s/count.m.rrna", args[1])
write.table(count.m.table.rrna,file=count.m.table.file,quote=FALSE,sep="\\t",row.names=FALSE)
colnames(count.m.table.tx) <- sub("X","",colnames(count.m.table.tx))
count.m.table.file <- sprintf("%s/count.m.tx", args[1])
write.table(count.m.table.tx,file=count.m.table.file,quote=FALSE,sep="\\t",row.names=FALSE)
colnames(count.m.table.ng) <- sub("X","",colnames(count.m.table.ng))
count.m.table.file <- sprintf("%s/count.m.ng", args[1])
write.table(count.m.table.ng,file=count.m.table.file,quote=FALSE,sep="\\t",row.names=FALSE)

colnames(table.read.statistics) <- c("ID", "a", "b", "c", "d", "e", "f", "g", "h")
x.big <- xtable( as.data.frame(table.read.statistics),
                 display=c("d","d","s","s","s","s","s","s","s","s"),
                 align=c('l','r','r','r','r','r','r','r','r','r'),
                 label='count',
                 caption='{\\\\bf Summary statistics of short reads.}'
               )
print( x.big,
       caption.placement="top",
       include.rownames=FALSE )
print("Overlapping transcript IDs")
if (length(overlapping.tx) > 0) {
  print(overlapping.tx)
}
EOF

################################################################################
# Simulation and check
cat>$BASEDIR/job-simulate<<EOF
RSCRIPT=$CACRSCRIPT
for i in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$i)
  rm -f $RDATADIR/\$FASTQNUM.fq.gz
  \$RSCRIPT job-simulate.R \$FASTQNUM
  gzip $RDATADIR/\$FASTQNUM.fq
  echo Created File: $RDATADIR/\$FASTQNUM.fq.gz
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
txdb.file <- "$RDATADIR/$GENOMETXDB"
txdb <- loadDb(txdb.file)
rrna.txdb.file <- "$RDATADIR/$GENOMERRNATXDB"
rrna.txdb <- loadDb(rrna.txdb.file)
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

s.mutans.sequence <- readDNAStringSet(genome.file)
chrom.list <- names(s.mutans.sequence)

lengthChr1 <- 0
for (i in chrom.list) { 
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
    # Make sure to always combine/compare objects based on the same reference
    # genome (use suppressWarnings() to suppress this warning).
    read.GA <- c(read.GA,one.GA)

    read.simulated <- mapply(read.extract, 
                             unlist(read.start.end.strand[1,]),
                             unlist(read.start.end.strand[2,]), 
                             unlist(read.start.end.strand[3,]),
                             seq(length(unlist(read.start.end.strand[3,]))))
  }
}
olap <- summarizeOverlaps(feature.txnc,read.GA,mode="$COUNTMODE")
cl <- assays(olap)\$counts[,1]
save(cl,file=cl.file)
print(paste("See cl file",cl.file))
EOF

cat>$BASEDIR/job-check<<EOF
RSCRIPT=$CACRSCRIPT
for i in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$i)
  \$RSCRIPT job-check.R \$FASTQNUM
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
EOF
}

function batch3-run-qcalignde {
cat>$BASEDIR/run-qcalignde.sh<<EOF
#!/bin/bash
FASTQFILES=( $FASTQFILES )
sed s/PBSARRAYSIZE/\${#FASTQFILES[@]}/g < batch-qcalignde.sh > tbatch.sh
nsub tbatch.sh 
rm tbatch.sh
EOF

cat>$BASEDIR/batch-qcalignde.sh<<EOF
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
  cp \$PBS_O_WORKDIR/job-fastqc .
  cp \$PBS_O_WORKDIR/job-bwa-align .
  cp \$PBS_O_WORKDIR/job-de* .
  cp \$PBS_O_WORKDIR/job-de*.R .
  cp \$PBS_O_WORKDIR/job-stat* .

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CBWADIR

  # Copy common data
  cp $RDATADIR/$GENOMEFASTA $CDATADIR
  cp $RDATADIR/$GENOMETXDB $CDATADIR
  cp $RDATADIR/$GENOMERRNATXDB $CDATADIR
}

function process-data {
  cd \$TMPDIR
  FASTQFILES=( $FASTQFILES )
  g=\$((PBS_ARRAYID-1))
  NUM=\$(printf "%03d" \${FASTQFILES[\$g]})

  cp $RDATADIR/FASTQ\$NUM.fq.gz $CDATADIR

  echo 0. Plot quality scores per site of the chosen FASTQ file.
  bash job-stat \$NUM $CDATADIR/FASTQ\$NUM.fq.gz \\
    $CBWADIR/FASTQ\$NUM.fq.qualPlot.RData
  cp $CBWADIR/FASTQ\$NUM.fq.qualPlot.RData $RBWADIR

  echo 1. Filter out parts of low quality in the FASTQ file.
  bash job-fastqc \$NUM \\
    $CDATADIR/FASTQ\$NUM.fq.gz \\
    $CBWADIR/FASTQ\$NUM.prinseq.fq.gz 

#  bash job-stat \$NUM $CBWADIR/FASTQ\$NUM.prinseq.fq.gz \\
#    $CBWADIR/FASTQ\$NUM.prinseq.fq.qualPlot.RData
#  cp $CBWADIR/FASTQ\$NUM.prinseq.fq.qualPlot.RData $RBWADIR

  echo Count reads that pass the quality score filter.
  NUMBER_READ4=\$($ZCAT $CDATADIR/FASTQ\$NUM.fq.gz | wc -l)
  NUMBER_READ=\$((NUMBER_READ4 / 4))
  NUMBER_READPRIN4=\$($ZCAT $CBWADIR/FASTQ\$NUM.prinseq.fq.gz | wc -l)
  NUMBER_READPRIN=\$((NUMBER_READPRIN4 / 4))
   
  echo 2. Align the sort reads
  bash job-bwa-align \$NUM \\
    $CBWADIR/FASTQ\$NUM.prinseq.fq.gz \\
    $CBWADIR/FASTQ\$NUM.sorted
  if [ "$KEEPBAM" == "YES" ]; then
    ./samtools index $CBWADIR/FASTQ\$NUM.sorted.bam 
    cp $CBWADIR/FASTQ\$NUM.sorted.bam* $RBWADIR
  fi

  echo 3. Count short reads
  BAMFILE1=$CBWADIR/FASTQ\$NUM.sorted.bam
  bash job-de \$BAMFILE1 \$NUMBER_READ \$NUMBER_READPRIN

  echo 4. Copy the output file of read counts.
  cp \$BAMFILE1.cl $RBWADIR

  # To see the output files in the compute nodes.
  tree
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

function batch3-copy-localscripts {
  scp -qr pl $BASEDIR
}

function batch3-copy-scripts {
  ssh -x $CAC_USERHOST mkdir -p $CACWORKDIR
  scp -q $BASEDIR/*.sh $CAC_USERHOST:$CACWORKDIR
  scp -q $BASEDIR/*.R $CAC_USERHOST:$CACWORKDIR
  scp -q $BASEDIR/job* $CAC_USERHOST:$CACWORKDIR
  scp -qr pl $CAC_USERHOST:$CACWORKDIR
}

function batch3-rmessage {
  echo "bash $BASEDIR/push-data.sh"
  echo "work at cac:$CACWORKDIR"
  echo "bash $BASEDIR/get-data.sh"
  echo "bash $BASEDIR/clean-data.sh"
}

function batch3-localmessage {
  echo "bash $CBASEDIR/push-data.sh"
  echo work at $CBASEDIR
  echo "bash $CBASEDIR/get-data.sh"
  echo "bash $CBASEDIR/clean-data.sh"
}

# END OF BATCH3
################################################################################
