#!/bin/bash
# To sample about 5 percent run the following command:
# sample-fastq.sh input.fq.gz output.fq.gz 0.05
zcat $1 \
  | awk 'BEGIN {srand()} !/^$/ {if (rand() <= "'"$3"'" && NR%4==0) for ( o = 1; o <= 4; o++ ) {getline; print}}' \
  | gzip > $2 
