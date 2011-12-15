# (((((((((SmuUA159 SMU109-OMZ175) SMU21-1SM1) SMU44-11VS1) SMU20-15JP3)
# SMU86-U2A) SMU56-N29) SMU69-NLML4)      SMU57-NMT4863) SMU60-U138)

twoBitToFa Smutans.2bit Smutans.fa
for smu in SMU109-OMZ175 SMU21-1SM1 SMU44-11VS1 SMU20-15JP3 SMU86-U2A SMU56-N29 SMU69-NLML4 SMU57-NMT4863 SMU60-U138; do
  grep $smu Smutans.fa | sed "s/>//" > $smu
  twoBitToFa -seqList=$smu Smutans.2bit $smu.fa 
  rm $smu
done

# Download Smutans.2bit from the genome browser.
# or output/smu86
twoBitToFa Smutans.2bit Smutans.fa
grep SMU86 Smutans.fa | sed "s/>//" > SMU86
twoBitToFa -seqList=SMU86 Smutans.2bit Smutans.fa 
grep SMU86 Smutans.fa | sed "s/>//" > SMU86.only
diff SMU86 SMU86.only

# Extract known genes.
featureBits Smutans knownGenes -bed=output.bed -chrom=SMU86-U2A.contig1
# In SQL
hgsql
use Smutans
select * from knownGenes where chrom LIKE 'SMU86%' into outfile '/tmp/smu86-u2a.gene'
