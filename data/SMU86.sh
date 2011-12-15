# Reorder 
awk '{print $2"\t"$3"\t"$4"\t"$5"\t0\t"$7}' data/smu86-u2a.gene > output/smu86/1/bwa/feature-genome.out-geneonly
echo Check output/smu86/1/bwa/feature-genome.out-geneonly
head output/smu86/1/bwa/feature-genome.out-geneonly

