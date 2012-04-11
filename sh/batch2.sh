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
#  ucsc-data
    batch2-run-cram
    batch2-run-bwa-align

    batch2-run-samtools-pileup

    # This would replace fastq-qc, bwa-align, and de.
    # Because this creates similar job scripts to those created by fastq-qc,
    # bwa-align, or de, place qcalignde after them.
    batch2-run-qcalignde

    create-index

    # batch2-run-blast

  #  batch-run-parsernaseq
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
  # . reference genome fasta file in _REFGENOMEFASTA_
  REFGENOMEFASTA=$(grep ^REFGENOMEFASTA\: $SPECIESFILE | cut -d":" -f2)
  # . reference genome annotation file in _REFGENOMEGFF_
  REFGENOMEGFF=$(grep ^REFGENOMEGFF\: $SPECIESFILE | cut -d":" -f2)
  # . reference genome annotation TrancriptDb file in _REFGENOMETXDB_
  REFGENOMETXDB=$(grep ^REFGENOMETXDB\: $SPECIESFILE | cut -d":" -f2)
  REFGENOMETXDBBASE=$(basename $REFGENOMETXDB)
  # . cram reference genome fasta file in _CRAMGENOMEFASTA_
  CRAMGENOMEFASTA=$(grep ^CRAMGENOMEFASTA\: $SPECIESFILE | cut -d":" -f2)
  # . cram wall time in _CRAMWALLTIME_
  CRAMWALLTIME=$(grep ^CRAMWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  # . pipeline wall time in _QCALIGNDEWALLTIME_
  QCALIGNDEWALLTIME=$(grep ^QCALIGNDEWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  # . split file size limit in _MAXLINEDE_
  MAXLINEDE=$(grep ^MAXLINEDE\: $SPECIESFILE | cut -d":" -f2)
  # . minimum alignment mapping quality in _MINMAPQ_
  MINMAPQ=$(grep ^MINMAPQ\: $SPECIESFILE | cut -d":" -f2)
  # . number CPUs in a compute node in _BWAALIGNNCPU_
  BWAALIGNNCPU=$(grep ^BWAALIGNNCPU\: $SPECIESFILE | cut -d":" -f2)
  # . set bwa options in _BWAOPTION_
  BWAOPTION=$(grep ^BWAOPTION\: $SPECIESFILE | cut -d":" -f2)
  # . Rscript path in _CACRSCRIPT_
  CACRSCRIPT=$(grep ^CACRSCRIPT\: $SPECIESFILE | cut -d":" -f2)
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
  REFGENOMEPTT=$(grep ^REFGENOMEPTT\: $SPECIESFILE | cut -d":" -f2)
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

#scp output/data/bacteria.fa $CAC_USERHOST:$RDATADIR

echo Edit push-data.sh to send cram or fastq files.
for g in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$g)
  FASTQFILE=\$(grep ^\$FASTQNUM\: $SPECIESFILE | cut -d":" -f2)
  # scp \$FASTQFILE $CAC_USERHOST:$RDATADIR/\$FASTQNUM.fq.gz
  # scp $CRAMDIR/\$FASTQNUM.cram $CAC_USERHOST:$RDATADIR
done

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
# #PBS -M ${BATCHEMAIL}
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
  bash job-cram \$(printf "%03d" \${FASTQFILES[\$g]})
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
cat>$BASEDIR/run-bwa-align.sh<<EOF
#!/bin/bash
FASTQFILES=( $FASTQFILES )
sed s/PBSARRAYSIZE/\${#FASTQFILES[@]}/g < batch-bwa-align.sh > tbatch.sh
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
  cp \$HOME/$SAMTOOLS samtools
  cp \$HOME/$BWA bwa
  cp \$HOME/$SUBREADBUILDINDEX subread-buildindex
  cp \$HOME/$SUBREADALIGN subread-align

  # All of the batchjob scripts.
  cp \$PBS_O_WORKDIR/job-cram .
  cp \$PBS_O_WORKDIR/job-cram2fastq .
  cp \$PBS_O_WORKDIR/job-fastq-qc .
  cp \$PBS_O_WORKDIR/job-bwa-align .
  cp \$PBS_O_WORKDIR/job-de* .
  cp \$PBS_O_WORKDIR/job-de*.R .

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CBWADIR
  mkdir -p $CSUBREADDIR

  # Copy common data
  cp $RDATADIR/$GENOMEFASTA $CDATADIR
  cp $RDATADIR/$GENOMEGFF $CDATADIR
}

function retrieve-data {
  echo No Copy!
}

function process-data {
  cd \$TMPDIR
  CORESPERNODE=1
  FASTQFILES=( $FASTQFILES )
  g=\$((PBS_ARRAYID-1))
  NUM=\$(printf "%03d" \${FASTQFILES[\$g]})

  bash job-cram2fastq \$NUM \\
    $RDATADIR/FASTQ\$NUM.cram \\
    $CDATADIR/FASTQ\$NUM.recovered.fq
  gzip $CDATADIR/FASTQ\$NUM.recovered.fq

  bash job-bwa-align \$NUM \\
    $CDATADIR/FASTQ\$NUM.recovered.fq.gz \\
    $CBWADIR/FASTQ\$NUM.sorted \\
    $CSUBREADDIR/FASTQ\$NUM.sorted

  cp $CBWADIR/FASTQ\$NUM.sorted.bam $RBWADIR
  cp $CSUBREADDIR/FASTQ\$NUM.sorted.bam $RSUBREADDIR
}

copy-data
process-data
cd
rm -rf \$TMPDIR
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
# #PBS -M ${BATCHEMAIL}
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
# #PBS -M ${BATCHEMAIL}
#PBS -t 1-PBSARRAYSIZE


function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$HOME/$PRINSEQ pl
  cp \$HOME/$SAMTOOLS samtools
  cp \$HOME/$BWA bwa
  # cp \$HOME/$SUBREADBUILDINDEX subread-buildindex
  # cp \$HOME/$SUBREADALIGN subread-align

  # All of the batchjob scripts.
  cp \$PBS_O_WORKDIR/job-cram .
  cp \$PBS_O_WORKDIR/job-cram2fastq .
  cp \$PBS_O_WORKDIR/job-fastq-qc .
  cp \$PBS_O_WORKDIR/job-bwa-align .
  cp \$PBS_O_WORKDIR/job-de* .
  cp \$PBS_O_WORKDIR/job-de*.R .

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CBWADIR
  mkdir -p $CSUBREADDIR

  # Copy common data
  cp $RDATADIR/$GENOMEFASTA $CDATADIR
  cp $RDATADIR/$GENOMEGFF $CDATADIR
}

function retrieve-data {
  # cp $CBWADIR/*.cutadapt.fq.gz $RBWADIR
  # cp $CBWADIR/*.prinseq.fq.gz $RBWADIR
  echo No Copy!
}

function process-data {
  cd \$TMPDIR
  FASTQFILES=( $FASTQFILES )
  g=\$((PBS_ARRAYID-1))
  NUM=\$(printf "%03d" \${FASTQFILES[\$g]})
  # input: cram, and output: bam
  bash job-cram2fastq \$NUM \\
    $RDATADIR/FASTQ\$NUM.cram \\
    $CDATADIR/FASTQ\$NUM.recovered.fq
  gzip $CDATADIR/FASTQ\$NUM.recovered.fq
  bash job-fastq-qc \$NUM \\
    $CDATADIR/FASTQ\$NUM.recovered.fq.gz \\
    $CBWADIR/FASTQ\$NUM.prinseq.fq.gz 

  bash job-bwa-align \$NUM \\
    $CBWADIR/FASTQ\$NUM.prinseq.fq.gz \\
    $CBWADIR/FASTQ\$NUM.sorted \\
    $CSUBREADDIR/FASTQ\$NUM.sorted

# 1. Total number of short reads in fastq: zcat FASTQ051.fq.gz | wc -l
  BAMFILE1=$CBWADIR/FASTQ\$NUM.sorted.bam
  NUMBER_READ4=\$(zcat $CDATADIR/FASTQ\$NUM.recovered.fq.gz | wc -l)
  NUMBER_READ=\$((NUMBER_READ4 / 4))
  echo bash job-de \$BAMFILE1 \$NUMBER_READ
  bash job-de \$BAMFILE1 \$NUMBER_READ
  cp \$BAMFILE1.cl $RBWADIR

#  BAMFILE2=$CSUBREADDIR/FASTQ\$NUM.sorted.bam
#  CLFILE1=$RBWADIR/FASTQ\$NUM.cl
#  CLFILE2=$RSUBREADDIR/FASTQ\$NUM.cl
#  PREFIX1=BWA-
#  PREFIX2=SUBREAD-
#  for k in {1..2}; do
#    BAMFILE=BAMFILE\$k
#    CLFILE=CLFILE\$k
#    PREFIX=PREFIX\$k
#    ./samtools view \${!BAMFILE} | split -d -l $MAXLINEDE - \${!PREFIX}
#    for i in \`ls \${!PREFIX}*\`; do ./samtools view -Sb -T $CDATADIR/$GENOMEFASTA \$i > \$i.bam; done
#    for i in \`ls \${!PREFIX}*.bam\`; do
#      bash job-de \$i
#    done
#    bash job-de2 \${!CLFILE} \${!PREFIX}
#  done 
}

copy-data
process-data
retrieve-data
cd
rm -rf \$TMPDIR
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

print(sum(cl.bwa==cl.sim)/length(cl.sim))
EOF

cat>$BASEDIR/job-simulate<<EOF
RSCRIPT=$CACRSCRIPT
for i in $TESTFASTQ; do
  TESTFASTQNUM=FASTQ\$(printf "%03d" \$i)
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

txdb.file <- "$RDATADIR/$REFGENOMETXDBBASE"
gff.file <- "$RDATADIR/$GENOMEGFF"

# saveFeatures and loadFeatures
txdb <- loadFeatures(txdb.file)
feature.cds <- cds(txdb, columns="exon_name")
strand(feature.cds) <- '*'

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

stopifnot( length(names(s.mutans.sequence)) == 1 )
for (i in names(s.mutans.sequence)) { 
  cat (i,"\\n")

  geneIR <- ranges(feature.cds[seqnames(feature.cds)==i])

  geneIR <- geneIR[width(geneIR)>1000]
  geneIR <- geneIR - 100
  if (length(geneIR) > 0) {
    read.pos <- c(mapply(function(x,y) sample(x:y,size=10,replace=TRUE),start(geneIR),end(geneIR)))
    s.mutans <- s.mutans.sequence[[i]]

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

olap <- summarizeOverlaps(feature.cds,read.GA,mode="IntersectionStrict")

cl <- assays(olap)\$counts[,1]
save(cl,file=cl.file)
EOF

grep ^QUALITYSCORE $SPECIESFILE | sed s/:/=/ > $BASEDIR/job-cram
cat>>$BASEDIR/job-cram<<EOF
QUALITYSCORESEQUENCE=QUALITYSCORE\$1

cp $RDATADIR/$CRAMGENOMEFASTAFILENAME $CDATADIR
cp $CDATADIR/$CRAMGENOMEFASTABASENAME.fna $CDATADIR/$CRAMGENOMEFASTABASENAME.fa

cp $RDATADIR/FASTQ\$1.fq.gz $CBWADIR
GZIPFASTAQFILE=$CBWADIR/FASTQ\$1.fq.gz
FASTAQFILE=\${GZIPFASTAQFILE%.gz}

./bwa index -p $CBWADIR/$CRAMGENOMEFASTABASENAME-bwa -a is \\
  $CDATADIR/$CRAMGENOMEFASTABASENAME.fa

# We need more preprocessing of fastq files.
if [ "\${!QUALITYSCORESEQUENCE}" == "illumina" ]; then
  ADAPTERSEQUENCE=ADAPTER\$1
  mv \$GZIPFASTAQFILE $CBWADIR/temp.FASTQ\$1.fq.gz
else
  zcat \$GZIPFASTAQFILE \\
    | grep -A 3 '^@.* [^:]*:N:[^:]*:' \\
    | sed '/^--$/d' | gzip > $CBWADIR/temp.FASTQ\$1.fq.gz
  rm \$GZIPFASTAQFILE
fi

if [ "\${!QUALITYSCORESEQUENCE}" == "illumina" ]; then
  ./bwa aln -I -t $BWAALIGNNCPU \\
    $CBWADIR/$CRAMGENOMEFASTABASENAME-bwa \\
    $CBWADIR/temp.FASTQ\$1.fq.gz > $CBWADIR/FASTQ\$1.sai
else
  ./bwa aln -t $BWAALIGNNCPU \\
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
  --input-bam-file $CBWADIR/FASTQ\$1.sorted.bam \\
  --reference-fasta-file $CDATADIR/$CRAMGENOMEFASTABASENAME.fa \\
  --output-cram-file $CDATADIR/FASTQ\$1.cram
rm $CBWADIR/FASTQ\$1.sorted.bam

cp $CDATADIR/FASTQ\$1.cram $RDATADIR

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
grep ^ADAPTER $SPECIESFILE | sed s/:/=/ > $BASEDIR/job-fastq-qc
grep ^QUALITYSCORE $SPECIESFILE | sed s/:/=/ >> $BASEDIR/job-fastq-qc
cat>>$BASEDIR/job-fastq-qc<<EOF
ADAPTERSEQUENCE=ADAPTER\$1
QUALITYSCORESEQUENCE=QUALITYSCORE\$1

$PYTHON $CUTADAPT --minimum-length=25 \\
  -a \${!ADAPTERSEQUENCE} \\
  -o $CBWADIR/FASTQ\$1.cutadapt.fq.gz \\
  \$2 &> /dev/null
# rm \$2

gzip -dc $CBWADIR/FASTQ\$1.cutadapt.fq.gz | \\
  perl pl/prinseq-lite.pl \\
    -fastq stdin \\
    -trim_qual_right 20 \\
    -rm_header \\
    -out_good stdout | \\
  gzip > $CBWADIR/outfile
rm $CBWADIR/FASTQ\$1.cutadapt.fq.gz
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

./bwa index -p $CBWADIR/$GENOMEFASTA-bwa -a is \\
  $CDATADIR/$GENOMEFASTA &> /dev/null

./subread-buildindex -o \\
  $CSUBREADDIR/$GENOMEFASTA-subread \\
  $CDATADIR/$GENOMEFASTA &> /dev/null

./bwa aln -t $BWAALIGNNCPU \\
  $BWAOPTION \\
  $CBWADIR/$GENOMEFASTA-bwa \\
  \$2 > $CBWADIR/FASTQ\$1.sai 2> /dev/null

./bwa samse -n 1 \\
  -r "@RG\\tID:$PROJECTNAME\\tSM:BWA" \\
  $CBWADIR/$GENOMEFASTA-bwa \\
  $CBWADIR/FASTQ\$1.sai \\
  \$2 \\
  | ./samtools view -Sb -q $MINMAPQ - > $CBWADIR/FASTQ\$1.bam 

./samtools sort $CBWADIR/FASTQ\$1.bam \$3 &> /dev/null
rm $CBWADIR/FASTQ\$1.sai
rm $CBWADIR/FASTQ\$1.bam
# FIXME:subread
exit

# Subread alignment
gunzip \$2
FASTAQFILE=\${2%.gz}
./subread-align \\
  --threads $BWAALIGNNCPU \\
  --phred 3 \\
  --unique \\
  -i $CSUBREADDIR/$GENOMEFASTA-subread \\
  -r \$FASTAQFILE \\
  -o $CSUBREADDIR/FASTQ\$1.sam &> /dev/null

./samtools view -bS -o $CSUBREADDIR/FASTQ\$1.bam \\
  -q $MINMAPQ \\
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
if (length(args) != 2)
{
  cat ("Rscript job-de.R bam.file 1000\n")
  quit("yes")
}
txdb.file <- "$RDATADIR/$REFGENOMETXDBBASE"
bam.file <- paste(args[1])
num.total.read <- args[2]
cl.file <- paste(args[1], "cl", sep=".")

# saveFeatures and loadFeatures
txdb <- loadFeatures(txdb.file)

# 1. Total number of short reads in fastq: zcat FASTQ051.fq.gz | wc -l
indexBam(bam.file)
bam.what <- c("qname","flag","pos","mapq","cigar")
bvAll <- scanBam(bam.file, param=ScanBamParam(what=bam.what))
bv <- readBamGappedAlignments(bam.file,use.names=TRUE,param=ScanBamParam(what=c("mapq")))
# 2. Number of short reads mapped 
# > length(bv)
num.mapped.read <- length(bv)
# 3. Number of uniquely mapped short reads
# Select mapped reads with mapq greater than or equal to $MINMAPQ
bv <- bv[elementMetadata(bv)["mapq"][,1] >= $MINMAPQ]
num.unique.read <- length(bv)
# 4. Number of short reads uniquely mapped on CDS regions
# Map the uniquely mapped reads on CDS regions
# Remove strands information
feature.cds <- cds(txdb, columns="exon_name")
strand(feature.cds) <- '*'
strand(bv) <- '*'
stopifnot(length(runValue(seqnames(feature.cds)))==1)
seqnames(bv) <- rep(runValue(seqnames(feature.cds)),length(bv))
olap <- summarizeOverlaps(feature.cds, bv,mode="IntersectionStrict")
sum(assays(olap)\$counts)
cl <- assays(olap)\$counts[,1]

# 5. Number of short reads uniquely mapped on non-coding regions
# Flattened tx
feature.tx.flat <- transcripts(txdb)
strand(feature.tx.flat) <- '*'
feature.tx.flat <- reduce(feature.tx.flat)

feature.nc <- transcripts(txdb)
strand(feature.nc) <- '*'
feature.nc <- gaps(reduce(feature.nc))
feature.nc <- feature.nc[-2:-1]
# sum(width(feature.nc))
stopifnot(length(runValue(seqnames(feature.nc)))==1)
seqnames(bv) <- rep(runValue(seqnames(feature.nc)),length(bv))
olap <- summarizeOverlaps(feature.nc,bv,mode="IntersectionStrict")
cl.nc <- assays(olap)\$counts[,1]
# [1] 510054
sum(width(feature.nc)) + sum(width(feature.tx.flat)) 

save(num.total.read, num.mapped.read, num.unique.read, cl, cl.nc,file=cl.file)

EOF

cat>$BASEDIR/job-de<<EOF
RSCRIPT=$CACRSCRIPT
\$RSCRIPT job-de.R \$1 \$2
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
  \$RSCRIPT job-de-sum.R \${!RALIGNDIR} \$FILE
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
feature.cds <- cds(txdb, columns="exon_name")

count.table <- data.frame(gene=unlist(elementMetadata(feature.cds)["exon_name"][,1]))

table.read.statistics <- c()
fastQIndex <- scan(args[2])
for (i in fastQIndex) {
  cl.file <- sprintf("%s/FASTQ%03d.sorted.bam.cl", args[1], i)
  load(cl.file)
  table.read.statistics <- 
    rbind(table.read.statistics, c(i, num.total.read, num.mapped.read, num.unique.read, sum(cl), sum(cl.nc)))
  count.table <- data.frame(count.table,cl)
  colnames(count.table)[ncol(count.table)] <- paste("X",i,sep="")
}
colnames(count.table) <- sub("X","",colnames(count.table))
count.table.file <- sprintf("%s/count.txt", args[1])
write.table(count.table,file=count.table.file,quote=FALSE,sep="\\t",row.names=FALSE)

colnames(table.read.statistics) <- c("Sample ID", "Total number of reads", "Mapped reads", "Unique reads", "Reads on CDS regions", "Reads on noncoding regions")
x.big <- xtable( as.data.frame(table.read.statistics),
                 display=c("d","d","d","d","d","d"),
                 align=c('r','r','r','r','r','r'),
                 label='count',
                 caption='{\\bf Summary statistics of short reads.}'
               )
print( x.big,
       caption.placement="top",
       include.rownames=TRUE )
EOF
}

function batch2-get-data {
  GENOMEGFF=$(basename $REFGENOMEGFF)
cat>$BASEDIR/get-data.sh<<EOF
#!/bin/bash
# scp $CAC_USERHOST:$RBWADIR/*.de $BWADIR
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
# scp $CAC_USERHOST:$RSUBREADDIR/count.txt $SUBREADDIR/count-$REFGENOMEID.txt
# scp $CAC_USERHOST:$RDATADIR/*.cram $DATADIR
# scp $CAC_USERHOST:$RDATADIR/$GENOMEGFF.fa $DATADIR

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
##################################################



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
