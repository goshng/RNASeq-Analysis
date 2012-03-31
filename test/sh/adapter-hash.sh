declare -A adapter
adapter=( ["ATCACG"]=AR001 ["CGATGT"]=AR002 ["TTAGGC"]=AR003 ["TGACCA"]=AR004 ["ACAGTG"]=AR005 ["GCCAAT"]=AR006 ["CAGATC"]=AR007 ["ACTTGA"]=AR008 ["GATCAG"]=AR009 ["TAGCTT"]=AR010 ["GGCTAC"]=AR011 ["CTTGTA"]=AR012 )

declare -A as
as=( ["AR001"]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG ["AR002"]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG ["AR003"]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG ["AR004"]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG ["AR005"]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG ["AR006"]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG ["AR007"]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG ["AR008"]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG ["AR009"]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG ["AR010"]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG ["AR011"]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG ["AR012"]=GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG )

# TruSeq Adapter, Index 1                 ******
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
# TruSeq Adapter, Index 2
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG
# TruSeq Adapter, Index 3
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG
# TruSeq Adapter, Index 4
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG
# TruSeq Adapter, Index 5
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG
# TruSeq Adapter, Index 6
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG
# TruSeq Adapter, Index 7
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
# TruSeq Adapter, Index 8
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG
# TruSeq Adapter, Index 9
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG
# # TruSeq Adapter, Index 10
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG
# TruSeq Adapter, Index 11
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG
# TruSeq Adapter, Index 12
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG
# TruSeq Adapter, Index 13
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG
# TruSeq Adapter, Index 14
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG
# TruSeq Adapter, Index 15
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG
# TruSeq Adapter, Index 16
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG
# TruSeq Adapter, Index 18
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG 
# TruSeq Adapter, Index 19 
# 5’ GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG





j=55
for i in /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2397_SAG_82.1M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2397_SAG_82.2M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2397_SAG_82.3M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2397_SAG_61.1M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2397_SAG_61.2M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2397_SAG_61.3M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2399_SAG_61.1B.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2399_SAG_61.2B.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2396_SAG_56.1M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2396_SAG_56.2M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2396_SAG_56.3M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2399_SAG_56.2B.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2399_SAG_56.3B.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/8_2395_SAG_2.1M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/8_2395_SAG_2.2M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/8_2395_SAG_2.3M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2398_SAG_2.1B.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2398_SAG_2.2B.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2398_SAG_2.3B.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/8_2395_SAG_13.1M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/8_2395_SAG_13.2M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/8_2395_SAG_13.3M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2398_SAG_13.1B.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2398_SAG_13.2B.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2398_SAG_13.3B.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2396_SAG_14.1M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2396_SAG_14.2M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2396_SAG_14.3M.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2399_SAG_14.1B.fastq.gz /Volumes/Elements/Documents/Projects/rnaseq/data/RNA_seq.4Aug11/2399_SAG_14.2B.fastq.gz ; do
  j=$((j+1))
  sixletters=$(zcat $i | head -n 1 | cut -c49-)
  # echo "$j $sixletters ${adapter[$sixletters]}"
  #echo "$sixletters ${adapter[$sixletters]}"
  NUM=$(printf "%03d" $j)
  arnum=${adapter[$sixletters]}
  echo "$arnum ADAPTER$NUM:${as[$arnum]}"
  # echo $sixletters
done
