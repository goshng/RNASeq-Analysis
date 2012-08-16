# I was asked to add genome tracks for RNAseq BAM files.

# Run this to create tables for BAM tracks at a computer where UCSC genome
# browser is setup.
for i in 2 3 4 5 15 16 17 18 19 20 21 22 24 25 26 27 23 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55; do
  FASTQNUM=FASTQ$(printf "%03d" $i)
  hgBbiDbLink SmuUA159v2 ${FASTQNUM}rnaseqBam /gbdb_cornell/SmuUA159v2/bbi/${FASTQNUM}.sorted.bam
done

