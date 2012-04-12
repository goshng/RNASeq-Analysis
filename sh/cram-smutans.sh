#!/bin/bash -x

# FASTQ=FASTQ001.fq
FASTQ=test.fq
REFFASTA=NC_004350.fa

BWA_BIN=$HOME/bin/bwa
SAMTOOLS_BIN=$HOME/bin/samtools
CRAM_BIN=$HOME/usr/bin/cramtools.jar
PICARD_SAMTOFASTAQ=$HOME/usr/bin/picard-tools-1.62/SamToFastq.jar 
PICARD_VALIDATESAMFILE=$HOME/usr/bin/picard-tools-1.62/ValidateSamFile.jar
PICARD_READGROUP=$HOME/usr/bin/picard-tools-1.62/AddOrReplaceReadGroups.jar

CPU=8

java -Xmx10g -jar \
  $PICARD_SAMTOFASTAQ \
  INPUT=$FASTQ.readgroup.bam \
  FASTQ=$FASTQ.recovered.fq
exit

# $BWA_BIN index $REFFASTA
$BWA_BIN aln -I -t $CPU $REFFASTA $FASTQ > $FASTQ.sai
$BWA_BIN samse $REFFASTA $FASTQ.sai $FASTQ | $SAMTOOLS_BIN view -Sb - > $FASTQ.bam
# $SAMTOOLS_BIN faidx $REFFASTA

$SAMTOOLS_BIN sort $FASTQ.bam $FASTQ.sorted
$SAMTOOLS_BIN index $FASTQ.sorted.bam

java -jar $CRAM_BIN cram \
  --input-bam-file $FASTQ.sorted.bam \
  --reference-fasta-file $REFFASTA \
  --output-cram-file $FASTQ.cram

$SAMTOOLS_BIN view -h $FASTQ.sorted.bam -o $FASTQ.sam

java -Xmx2g -jar \
  $PICARD_READGROUP \
  RGLB=testrglb \
  RGPL=illumina \
  RGPU=testrgpu \
  RGSM=testrgsm \
  INPUT=$FASTQ.sam \
  OUTPUT=$FASTQ.readgroup.bam

java -Xmx2g -jar \
  $PICARD_VALIDATESAMFILE \
  INPUT=$FASTQ.readgroup.bam

java -Xmx10g -jar \
  $PICARD_SAMTOFASTAQ \
  INPUT=$FASTQ.readgroup.bam \
  FASTQ=$FASTQ.recovered.fq

