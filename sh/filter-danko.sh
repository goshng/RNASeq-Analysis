#mkdir filtered
for fastq in *.fastq.gz; do
  zcat $fastq | grep -A 3 '^@.* [^:]*:N:[^:]*:' | sed '/^--$/d' > filtered/$fastq
done

