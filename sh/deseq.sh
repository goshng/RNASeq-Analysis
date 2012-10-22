###############################################################################
# Copyright (C) 2011-2012 Sang Chul Choi
#
# This file is part of RNAseq Analysis.
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
function deseq {
  select SPECIES in ${SPECIESS[@]}; do 
  if [ "$SPECIES" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  else  
    batch3-repetition
    batch3-variable
    deseq-speciesfile 

    #deseq-push-data
    #deseq-get-data
    deseq-make-job
    deseq-run-deseq
    #create-index
    #deseq-copy-scripts
    deseq-rmessage
    break
  fi
  done
}

# Set species file variables
function deseq-speciesfile {
  DESEQ=$(grep ^DESEQ\: $SPECIESFILE | cut -d":" -f2)
  # Use the following to find pairs of comparisons.
  OIFS=$IFS; IFS=','; DESEQCOMPARISON=($DESEQ); IFS=$OIFS
  # echo ${DESEQCOMPARISON[0]}
  # echo ${DESEQCOMPARISON[1]}

  # Two CDS (this may not be a correct term) FASTA files
  REFCDSFASTA1=$(grep ^REFCDSFASTA1\: $SPECIESFILE | cut -d":" -f2)
  REFCDSFASTA1=$ROOTANALYSISDIR/$REFCDSFASTA1
  CDSFASTA1=$(basename $REFCDSFASTA1)
  REFCDSFASTA2=$(grep ^REFCDSFASTA2\: $SPECIESFILE | cut -d":" -f2)
  REFCDSFASTA2=$ROOTANALYSISDIR/$REFCDSFASTA2
  CDSFASTA2=$(basename $REFCDSFASTA2)
  REF1GFF=$(grep ^REF1GFF\: $SPECIESFILE | cut -d":" -f2)
  REF1GFF=$ROOTANALYSISDIR/$REF1GFF

  # No used ?
  FASTQFILES=$(grep ^FASTQFILES\: $SPECIESFILE | cut -d":" -f2)
  FASTQLABEL=$(grep ^FASTQLABEL\: $SPECIESFILE | cut -d":" -f2)
  CACWORKDIR=$(grep ^CACWORKDIR\: $SPECIESFILE | cut -d":" -f2)
  REFGENOMEFASTA=$(grep ^REFGENOMEFASTA\: $SPECIESFILE | cut -d":" -f2)
  REFGENOMEFASTA=$ROOTANALYSISDIR/$REFGENOMEFASTA
  GENOMEFASTA=$(basename $REFGENOMEFASTA)
  REFGENOMETXDB=$(grep ^REFGENOMETXDB\: $SPECIESFILE | cut -d":" -f2)  
  REFGENOMETXDB=$ROOTANALYSISDIR/$REFGENOMETXDB
  GENOMETXDB=$(basename $REFGENOMETXDB)
  QCALIGNDEWALLTIME=$(grep ^QCALIGNDEWALLTIME\: $SPECIESFILE | cut -d":" -f2)
  MINMAPQ=$(grep ^MINMAPQ\: $SPECIESFILE | cut -d":" -f2)
  MINTRIMQ=$(grep ^MINTRIMQ\: $SPECIESFILE | cut -d":" -f2)
}

# Send input files to the remote machine.
function deseq-push-data {
cat>$BASEDIR/push-data.sh<<EOF
#!/bin/bash
ssh -x $CAC_USERHOST mkdir -p $RBWADIR
ssh -x $CAC_USERHOST mkdir -p $RDATADIR
scp $REFGENOMEFASTA $CAC_USERHOST:$RDATADIR
scp $REFGENOMETXDB $CAC_USERHOST:$RDATADIR
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

function deseq-get-data {
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
for i in $FASTQFILES; do
  FASTQNUM=FASTQ\$(printf "%03d" \$i)
  scp $CAC_USERHOST:$RBWADIR/\$FASTQNUM*fq.qualPlot.* $BWADIR
done 
scp $CAC_USERHOST:$RBWADIR/stat1.tex $BWADIR
echo "Check $BWADIR"
EOF
}

################################################################################
# jobs
# 1. job-deseq
################################################################################
function deseq-make-job {

################################################################################
# 
################################################################################
cat>$BASEDIR/job-deseq<<EOF
RSCRIPT=$CACRSCRIPT
DESEQ=$DESEQ
OIFS=\$IFS; IFS=','; DESEQCOMPARISON=(\$DESEQ); IFS=\$OIFS
for i in "\${DESEQCOMPARISON[@]}"; do
  OIFS=\$IFS; IFS='/'; SAMPLEAB=(\$i); IFS=\$OIFS
  CSVFILE=$RUNANALYSIS/\${SAMPLEAB[0]}-\${SAMPLEAB[1]}
  \$RSCRIPT $BASEDIR/job-deseq.R \${SAMPLEAB[0]} \${SAMPLEAB[1]} \$CSVFILE
done 
echo "Check $RUNANALYSIS"
EOF

cat>$BASEDIR/job-deseq.R<<EOF
library(smutans)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3)
{
  cat ("Rscript job-deseq.R UA159 1SM1 out\n")
  quit("yes")
}
sampleA <- args[1]
sampleB <- args[2]
csvFile <- paste(args[3],"csv",sep=".")
clustFile <- paste(args[3],"pdf",sep=".")
sampleTitle <- paste(sampleA, sampleB)
sGenes <-
  readSmutans (countsFile="$BWADIR/count.cds", 
               indexFile ="$RUNANALYSIS/count.txt.index", 
               condition = c(sampleA, sampleB),
               firstFactorLabel = c("reference"),
               secondFactorLabel = c(sampleA, sampleB))
omz <- newSmutans( sGenes, title=sampleTitle )

pdf(clustFile)
smutans.de2Clust( omz )
dev.off()
omz <- smutans.de2( omz, type="reference", 
                    condA=sampleA, condB=sampleB )
smutans.de2List( omz, csvFile )
EOF

cat>$BASEDIR/job-blast<<EOF
# perl pl/extractcds.pl -f Genbank ../NC_004350.gbk > NC_004350.faa

$HOME/$MAKEBLASTDB -in $REFCDSFASTA1 \\
  -dbtype prot -title chr1 -input_type fasta \\
  -out $RUNANALYSIS/chr1
$HOME/$MAKEBLASTDB -in $REFCDSFASTA2 \\
  -dbtype prot -title chr2 -input_type fasta \\
  -out $RUNANALYSIS/chr2

# BLASTP chr1 and chr2
# Check 'blastp -help' for the output format
$HOME/$BLASTP -query $REFCDSFASTA1 -db $RUNANALYSIS/chr2 \\
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qseq slen sseq' \\
  -num_threads 1 \\
  -evalue 0.001 \\
  -out $RUNANALYSIS/chr1-chr2
$HOME/$BLASTP -query $REFCDSFASTA2 -db $RUNANALYSIS/chr1 \\
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qseq slen sseq' \\
  -num_threads 1 \\
  -evalue 0.001 \\
  -out $RUNANALYSIS/chr2-chr1

# Extract sequence alignments.
perl pl/extract-msa.pl unique \\
  -gene $REFCDSFASTA1 \\
  -coverage 0.5 \\
  -in $RUNANALYSIS/chr1-chr2 \\
  -out $RUNANALYSIS/genes-chr1
perl pl/extract-msa.pl unique \\
  -gene $REFCDSFASTA2 \\
  -coverage 0.5 \\
  -in $RUNANALYSIS/chr2-chr1 \\
  -out $RUNANALYSIS/genes-chr2
rm $RUNANALYSIS/genes-chr1 $RUNANALYSIS/chr1-chr2 
rm $RUNANALYSIS/chr1.phr $RUNANALYSIS/chr1.pin $RUNANALYSIS/chr1.psq 
rm $RUNANALYSIS/chr2-chr1 $RUNANALYSIS/genes-chr2
rm $RUNANALYSIS/chr2.phr $RUNANALYSIS/chr2.pin $RUNANALYSIS/chr2.psq 
# Find 
Rscript $BASEDIR/job-blast.R $RUNANALYSIS/genes-chr2-similar.csv \\
  $RUNANALYSIS/genes-chr2-similar.csv
echo "Check the following files:"
echo "  $RUNANALYSIS/genes-chr1-unique.csv"
echo "  $RUNANALYSIS/genes-chr1-similar.csv"
echo "  $RUNANALYSIS/genes-chr2-unique.csv"
echo "  $RUNANALYSIS/genes-chr2-similar.csv"
echo "END OF FILES"
EOF

cat>$BASEDIR/job-blast.R<<EOF
library(GenomicFeatures)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2)
{
  cat ("Rscript job-blast.R genes-chr2-similar.csv output.csv\n")
  quit("yes")
}

sm.gff <- import.gff3("$REF1GFF")
sm.gene <- sm.gff[sm.gff\$type=="gene",]
sm.CDS <- sm.gff[sm.gff\$type=="CDS",]
gene <- read.csv(args[1])

cat("Please, wait!\n")
gene.products <- c()
if (length(unlist(sm.CDS\$Parent)) == nrow(sm.CDS)) {
  for (i in gene\$chr1) {
    gene.id <- sm.gene[sm.gene\$locus_tag == i,"ID"]\$ID
    gene.product <- sm.CDS[unlist(sm.CDS\$Parent)==gene.id,"product"]\$product
    gene.products <- c(gene.products, gene.product)
  }
} else {
  print("Error: Multple CDS in genes")
}

x <- data.frame(gene, product=gene.products)
colnames(x) <- c("chr2","chr1","Coverage","Product")
write.table(x,file=args[2],sep=",",quote=FALSE,row.names=FALSE)
EOF
}

function deseq-run-deseq {
cat>$BASEDIR/run-deseq.sh<<EOF

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
}

function process-data {
  cd \$TMPDIR
  FASTQFILES=( $FASTQFILES )
  g=\$((PBS_ARRAYID-1))
  NUM=\$(printf "%03d" \${FASTQFILES[\$g]})

  cp $RDATADIR/FASTQ\$NUM.fq.gz $CDATADIR

  # 0. Plot quality scores per site of the chosen FASTQ file.
  bash job-stat \$NUM $CDATADIR/FASTQ\$NUM.fq.gz \\
    $CBWADIR/FASTQ\$NUM.fq.qualPlot.RData
  cp $CBWADIR/FASTQ\$NUM.fq.qualPlot.RData $RBWADIR

  # 1. Filter out parts of low quality in the FASTQ file.
  bash job-fastqc \$NUM \\
    $CDATADIR/FASTQ\$NUM.fq.gz \\
    $CBWADIR/FASTQ\$NUM.prinseq.fq.gz 

#  bash job-stat \$NUM $CBWADIR/FASTQ\$NUM.prinseq.fq.gz \\
#    $CBWADIR/FASTQ\$NUM.prinseq.fq.qualPlot.RData
#  cp $CBWADIR/FASTQ\$NUM.prinseq.fq.qualPlot.RData $RBWADIR

  # Count reads that pass the quality score filter.
  NUMBER_READ4=\$(zcat $CDATADIR/FASTQ\$NUM.fq.gz | wc -l)
  NUMBER_READ=\$((NUMBER_READ4 / 4))
  NUMBER_READPRIN4=\$(zcat $CBWADIR/FASTQ\$NUM.prinseq.fq.gz | wc -l)
  NUMBER_READPRIN=\$((NUMBER_READPRIN4 / 4))
   
  # 2. Align the sort reads
  bash job-bwa-align \$NUM \\
    $CBWADIR/FASTQ\$NUM.prinseq.fq.gz \\
    $CBWADIR/FASTQ\$NUM.sorted

  # 3. Count short reads
  BAMFILE1=$CBWADIR/FASTQ\$NUM.sorted.bam
  bash job-de \$BAMFILE1 \$NUMBER_READ \$NUMBER_READPRIN

  # 4. Copy the output file of read counts.
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

function deseq-copy-scripts {
  ssh -x $CAC_USERHOST mkdir -p $CACWORKDIR
  scp -q $BASEDIR/*.sh $CAC_USERHOST:$CACWORKDIR
  scp -q $BASEDIR/*.R $CAC_USERHOST:$CACWORKDIR
  scp -q $BASEDIR/job* $CAC_USERHOST:$CACWORKDIR
  scp -qr pl $CAC_USERHOST:$CACWORKDIR
}

function deseq-rmessage {
  echo "bash $BASEDIR/job-deseq"
  echo "bash $BASEDIR/job-blast"
}
# END OF DESEQ
################################################################################
