# Reorder 
awk '{print $1"\t"$4"\t"$5"\t"$9"\t0\t"$7}' data/SMU109-OMZ175.gff > output/omz/1/bwa/feature-genome.out-geneonly
