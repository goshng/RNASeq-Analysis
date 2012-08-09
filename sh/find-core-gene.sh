# We create a BED file using NC_004350 txdb file.
# We replace the underscore of a gene name with a dot.
R --vanilla <<RSCRIPT
library(rtracklayer)
library(GenomicFeatures)
txdbFile <- "data/analysis/NC_004350.txdb"
txdb <- loadFeatures(txdbFile)
feature.tx <- transcripts(txdb)
feature.cds <- cds(txdb,columns="tx_id")
stopifnot(length(feature.cds) == length(unlist(elementMetadata(feature.cds)\$tx_id)))
x <- elementMetadata(feature.tx)\$tx_id %in% unlist(elementMetadata(feature.cds)\$tx_id)
feature.cds <- feature.tx[x]
elementMetadata(feature.cds) <- data.frame(tx_id=elementMetadata(feature.cds)\$tx_id, tx_name=elementMetadata(feature.cds)\$tx_name,type="CDS")
names(feature.cds) <- elementMetadata(feature.cds)[,"tx_name"]
names(feature.cds) <- sub("_",".",names(feature.cds))
export(feature.cds,"/tmp/cds.bed")
RSCRIPT

# Create core genes of UA159.
perl pl/find-core-gene.pl core \
  -orthomclgroupname \
	-mcl output/email/from/michael-stanhope/080312/gene_content.tab \
	-bed /tmp/cds.bed \
	-coregenome \
"SMU SMU101 SMU102 SMU103 SMU104 SMU105 SMU107 SMU108 SMU109 SMU10 SMU20 SMU21 \
SMU22 SMU26 SMU29 SMU33 SMU36 SMU3 SMU40 SMU41 SMU44 SMU50 SMU52 SMU53 SMU54 \
SMU56 SMU57 SMU58 SMU60 SMU61 SMU62 SMU63 SMU66 SMU68 SMU69 SMU70 SMU72 SMU74 \
SMU75 SMU76 SMU77 SMU78 SMU80 SMU81 SMU82 SMU83 SMU85 SMU86 SMU88 SMU89 SMU92 \
SMU93 SMU94 SMU95 SMU97 SMU98 SMU99 SMU9" \
	-out /tmp/core.gene
mv /tmp/core.gene output/data/UA159.core.gene
echo Check output/data/UA159.core.gene for core genes.

cat>/tmp/t.R<<EOF
library(rtracklayer)
library(GenomicFeatures)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{
  cat ("Rscript file.R SMU21\n")
	quit("yes")
}
txdbFile <- paste("output/denovo/", args[1], ".txdb", sep="")
txdb <- loadFeatures(txdbFile)
feature.tx <- transcripts(txdb)
feature.cds <- cds(txdb,columns="tx_id")
stopifnot(length(feature.cds) == length(unlist(elementMetadata(feature.cds)\$tx_id)))
x <- elementMetadata(feature.tx)\$tx_id %in% unlist(elementMetadata(feature.cds)\$tx_id)
feature.cds <- feature.tx[x]
elementMetadata(feature.cds) <- data.frame(tx_id=elementMetadata(feature.cds)\$tx_id, tx_name=elementMetadata(feature.cds)\$tx_name,type="CDS")
names(feature.cds) <- elementMetadata(feature.cds)[,"tx_name"]
export(feature.cds,"/tmp/cds.bed")
EOF

# Repeat for 7 strains of the Smu.
for i in SMU21 SMU44 SMU20 SMU56 SMU69 SMU57 SMU86; do
	Rscript /tmp/t.R $i
	# Create core genes of UA159.
	perl pl/find-core-gene.pl core \
		-orthomclgroupname \
		-mcl output/email/from/michael-stanhope/080312/gene_content.tab \
		-bed /tmp/cds.bed \
		-coregenome \
	"SMU SMU101 SMU102 SMU103 SMU104 SMU105 SMU107 SMU108 SMU109 SMU10 SMU20 SMU21 \
	SMU22 SMU26 SMU29 SMU33 SMU36 SMU3 SMU40 SMU41 SMU44 SMU50 SMU52 SMU53 SMU54 \
	SMU56 SMU57 SMU58 SMU60 SMU61 SMU62 SMU63 SMU66 SMU68 SMU69 SMU70 SMU72 SMU74 \
	SMU75 SMU76 SMU77 SMU78 SMU80 SMU81 SMU82 SMU83 SMU85 SMU86 SMU88 SMU89 SMU92 \
	SMU93 SMU94 SMU95 SMU97 SMU98 SMU99 SMU9" \
		-out /tmp/core.gene
	mv /tmp/core.gene output/data/$i.core.gene
	echo Check output/data/$i.core.gene for core genes.
done

# Use a core gene file to append a column to a DESeq's output file.
for i in UA159 SMU21 SMU44 SMU20 SMU56 SMU69 SMU57 SMU86; do
	perl -e '
	my $n = 0;
	while (<>)
	{
		if (/^Core genes\s+(\d+)/)
		{
			$n = $1;
			last;
		}
	}
	my $i = 0;
	while (<>)
	{
		if ($i > 0)
		{
			print;
		}
		if ($i == $n)
		{
			last;
		}
		$i++;
	}
	' < output/data/$i.core.gene > output/data/$i.core.gene.only
done

cat>/tmp/t.R<<EOF
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3)
{
  cat ("Rscript file.R 1.csv core.gene.only 2.csv\n")
	quit("yes")
}
x <- read.csv(args[1])
x <- x[-1]
y <- read.table(args[2])
z <- x\$id %in% y\$V1
xo <- data.frame(x,core=z)
write.csv(xo,file=args[3])
EOF

# The 6 strains genome annotations
Rscript /tmp/t.R output/email/to/sara-palmer/071712/UA159-11VS1-SMU44.csv \
  output/data/SMU44.core.gene.only \
	output/email/to/sara-palmer/080612/UA159-11VS1-SMU44.csv

Rscript /tmp/t.R output/email/to/sara-palmer/071712/UA159-15JP3-SMU20.csv \
  output/data/SMU20.core.gene.only \
  output/email/to/sara-palmer/080612/UA159-15JP3-SMU20.csv 

Rscript /tmp/t.R output/email/to/sara-palmer/071712/UA159-1SM1-SMU21.csv \
  output/data/SMU21.core.gene.only \
  output/email/to/sara-palmer/080612/UA159-1SM1-SMU21.csv

Rscript /tmp/t.R output/email/to/sara-palmer/071712/UA159-N29-SMU56.csv \
  output/data/SMU56.core.gene.only \
  output/email/to/sara-palmer/080612/UA159-N29-SMU56.csv

Rscript /tmp/t.R output/email/to/sara-palmer/071712/UA159-NLML4-SMU69.csv \
  output/data/SMU69.core.gene.only \
  output/email/to/sara-palmer/080612/UA159-NLML4-SMU69.csv

Rscript /tmp/t.R output/email/to/sara-palmer/071712/UA159-NMT4863-SMU57.csv \
  output/data/SMU57.core.gene.only \
  output/email/to/sara-palmer/080612/UA159-NMT4863-SMU57.csv

# Revert the gene annotation format back to new format's.
sed 's/\./_/' < output/data/UA159.core.gene.only > output/data/UA159.core.gene.only.tmp 
mv output/data/UA159.core.gene.only.tmp output/data/UA159.core.gene.only

# Using UA159 genome annotations
Rscript /tmp/t.R output/email/to/sara-palmer/071712/UA159-11VS1.csv \
  output/data/UA159.core.gene.only \
	output/email/to/sara-palmer/080612/UA159-11VS1.csv

Rscript /tmp/t.R output/email/to/sara-palmer/071712/UA159-15JP3.csv \
  output/data/UA159.core.gene.only \
  output/email/to/sara-palmer/080612/UA159-15JP3.csv 

Rscript /tmp/t.R output/email/to/sara-palmer/071712/UA159-1SM1.csv \
  output/data/UA159.core.gene.only \
  output/email/to/sara-palmer/080612/UA159-1SM1.csv

Rscript /tmp/t.R output/email/to/sara-palmer/071712/UA159-N29.csv \
  output/data/UA159.core.gene.only \
  output/email/to/sara-palmer/080612/UA159-N29.csv

Rscript /tmp/t.R output/email/to/sara-palmer/071712/UA159-NLML4.csv \
  output/data/UA159.core.gene.only \
  output/email/to/sara-palmer/080612/UA159-NLML4.csv

Rscript /tmp/t.R output/email/to/sara-palmer/071712/UA159-NMT4863.csv \
  output/data/UA159.core.gene.only \
  output/email/to/sara-palmer/080612/UA159-NMT4863.csv

# For Sara's work of CSP and pH
for i in allSample allSampleCsp allSampleph allSampleType \
  allSampleSmu21 allSampleSmu86 \
  csp galactoseOnly \
  glucoseOnly nocsp ph5only ph7only smu21ph smu86csp tw1only ua159csp \
	ua159only ua159ph; do
	Rscript /tmp/t.R output/email/to/sara-palmer/071412/$i.csv \
		output/data/UA159.core.gene.only \
		output/email/to/sara-palmer/080612/$i.csv
done

for i in smu21phSMU21; do
	Rscript /tmp/t.R output/email/to/sara-palmer/071412/$i.csv \
		output/data/SMU21.core.gene.only \
		output/email/to/sara-palmer/080612/$i.csv
done

for i in smu86cspSMU86; do
	Rscript /tmp/t.R output/email/to/sara-palmer/071412/$i.csv \
		output/data/SMU86.core.gene.only \
		output/email/to/sara-palmer/080612/$i.csv
done 

# Have not done for these.
#/Volumes/Elements/Documents/Projects/RNASeq-Analysis/output/email/to/sara-palmer/071412/npAndP.csv
#/Volumes/Elements/Documents/Projects/RNASeq-Analysis/output/email/to/sara-palmer/071412/npAndUA159.csv
#/Volumes/Elements/Documents/Projects/RNASeq-Analysis/output/email/to/sara-palmer/071412/omz.csv
#/Volumes/Elements/Documents/Projects/RNASeq-Analysis/output/email/to/sara-palmer/071412/omzSMU109.csv
#/Volumes/Elements/Documents/Projects/RNASeq-Analysis/output/email/to/sara-palmer/071412/pAndUA159.csv


