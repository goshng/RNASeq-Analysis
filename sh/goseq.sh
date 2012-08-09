###############################################################################
# Copyright (C) 2012 Sang Chul Choi
#
# This file is part of RNAseq Analysis.
# 
# RNAseq Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RNAseq Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RNAseq Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

REPETITION=1
function goseq {
  select SPECIES in ${SPECIESS[@]}; do 
  if [ "$SPECIES" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  else  
    batch3-repetition
    batch3-variable
    batch3-speciesfile 
    goseq-speciesfile 
    goseq-push-data
    goseq-get-data
    goseq-make-job
    goseq-run-goseq
    batch3-copy-scripts
    goseq-rmessage
    break
  fi
  done
}

# Set species file variables
function goseq-speciesfile {
  GOAUNIPROT=$(grep ^GOAUNIPROT\: $SPECIESFILE | cut -d":" -f2)
  GOAUNIPROT=$ROOTANALYSISDIR/$GOAUNIPROT
  GOBO=$(grep ^GOBO\: $SPECIESFILE | cut -d":" -f2)
  GOBO=$ROOTANALYSISDIR/$GOBO
  UNIREFFASTA=$(grep ^UNIREFFASTA\: $SPECIESFILE | cut -d":" -f2)
  UNIREFFASTA=$ROOTANALYSISDIR/$UNIREFFASTA
  GOSEQFASTA=$(grep ^GOSEQFASTA\: $SPECIESFILE | cut -d":" -f2)
  GOSEQFASTA=$ROOTANALYSISDIR/$GOSEQFASTA
  REF1GFF=$(grep ^REF1GFF\: $SPECIESFILE | cut -d":" -f2)
  REF1GFF=$ROOTANALYSISDIR/$REF1GFF
  REFGENOMETXDB=$(grep ^REFGENOMETXDB\: $SPECIESFILE | cut -d":" -f2)  
  REFGENOMETXDB=$ROOTANALYSISDIR/$REFGENOMETXDB
  GENOMETXDB=$(basename $REFGENOMETXDB)

  BLASTEVALUE=$(grep ^BLASTEVALUE\: $SPECIESFILE | cut -d":" -f2)  

  DESEQ=$(grep ^DESEQ\: $SPECIESFILE | cut -d":" -f2)
  OIFS=$IFS; IFS=','; DESEQCOMPARISON=($DESEQ); IFS=$OIFS
}

# Send input files to the remote machine.
function goseq-push-data {
cat>$BASEDIR/goseq-push-data.sh<<EOF
#!/bin/bash
ssh -x $CAC_USERHOST mkdir -p $RDATADIR
ssh -x $CAC_USERHOST mkdir -p $RRUNANALYSIS
scp $REF1GFF $CAC_USERHOST:$RDATADIR
scp $GOAUNIPROT $CAC_USERHOST:$RDATADIR
scp $GOBO $CAC_USERHOST:$RDATADIR
scp $UNIREFFASTA $CAC_USERHOST:$RDATADIR
scp $GOSEQFASTA $CAC_USERHOST:$RDATADIR
scp $REFGENOMETXDB $CAC_USERHOST:$RDATADIR
echo "Check cac:$RDATADIR"
EOF
}

function goseq-get-data {
cat>$BASEDIR/goseq-get-data.sh<<EOF
#!/bin/bash
scp $CAC_USERHOST:$RRUNANALYSIS/goseq.blast $RUNANALYSIS
scp $CAC_USERHOST:$RRUNANALYSIS/smutans.go2ngene $RUNANALYSIS
scp $CAC_USERHOST:$RRUNANALYSIS/smutans.gene2go $RUNANALYSIS
echo "Check $RUNANALYSIS"
EOF
}

function goseq-make-job {
cat>$BASEDIR/job-feature<<EOF
RSCRIPT=Rscript
\$RSCRIPT $BASEDIR/job-feature.R
EOF

cat>$BASEDIR/job-feature.R<<EOF
library(GenomicFeatures)
library(rtracklayer)
txdb <- loadFeatures("$REFGENOMETXDB")
feature.tx <- transcripts(txdb) 
feature.cds <- cds(txdb,columns="tx_id")
x <- elementMetadata(feature.tx)\$tx_id %in% 
     unlist(elementMetadata(feature.cds)\$tx_id)
feature.cds <- feature.tx[x]
names(feature.cds) <- elementMetadata(feature.cds)[,"tx_name"]
export(feature.cds,"$RUNANALYSIS/feature.genes",format="bed")
cat("Check $RUNANALYSIS/feature.genes\n")
EOF

cat>$BASEDIR/job-goseq<<EOF
RSCRIPT=Rscript
DESEQ=$DESEQ
OIFS=\$IFS; IFS=','; DESEQCOMPARISON=(\$DESEQ); IFS=\$OIFS
for i in "\${DESEQCOMPARISON[@]}"; do
  OIFS=\$IFS; IFS='/'; SAMPLEAB=(\$i); IFS=\$OIFS
  CSVFILE=$RUNANALYSIS/\${SAMPLEAB[0]}-\${SAMPLEAB[1]}
  \$RSCRIPT $BASEDIR/job-goseq.R \${SAMPLEAB[0]} \${SAMPLEAB[1]} \\
    \$CSVFILE \\
    $RUNANALYSIS/feature.genes \\
    $RUNANALYSIS/smutans.gene2go \\
    $RUNANALYSIS/smutans.go2ngene > \$CSVFILE.goseq
  echo Check \$CSVFILE.goseq
done 
echo "Check $RUNANALYSIS"
EOF

cat>$BASEDIR/job-goseq.R<<EOF
library(smutans)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6)
{
  cat ("Rscript job-simulate.R UA159 1SM1 out feature.genes gene2go go2ngene\n")
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

f <- args[4]
smutans.feature.genes = read.table(file=f,head=F)
stopifnot(ncol(smutans.feature.genes)==6) 
f <- args[5]
smutans.go.genes <- read.table(file=f,head=F)
stopifnot(ncol(smutans.go.genes)==3) 
f <- args[6]
smutans.cat.desc = read.table(file=f,head=F,sep="\t",quote="")
stopifnot(ncol(smutans.cat.desc)==3) 

cat("===============================\n")
smutans.de2List( omz, qval=0.05 )
cat("===============================\n")
smutans.de2Goseq ( omz, qval=0.05, 
                   feature.gene=smutans.feature.genes,
                   go.genes=smutans.go.genes,
                   cat.desc=smutans.cat.desc )
cat("===============================\n")

EOF

}

function goseq-run-goseq {
cat>$BASEDIR/run-goseq.sh<<EOF
#!/bin/bash
nsub batch-goseq.sh
EOF

cat>$BASEDIR/batch-goseq.sh<<EOF
#!/bin/bash
#PBS -l walltime=${QCALIGNDEWALLTIME}:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-GOSEQ
#PBS -q ${QUEUENAME}
#PBS -m e
$EMAILON#PBS -M ${BATCHEMAIL}

function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$HOME/$MAKEBLASTDB .
  cp \$HOME/$BLASTP .

  # All of the batchjob scripts.
  # No batchjob script.

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CRUNANALYSIS

  # No data files.
}

function process-data {
  cd \$TMPDIR

  GOAUNIPROT=\$(basename $GOAUNIPROT)
  GOBO=\$(basename $GOBO)
  UNIREFFASTA=\$(basename $UNIREFFASTA)
  GOSEQFASTA=\$(basename $GOSEQFASTA)
  REF1GFF=\$(basename $REF1GFF)

  # 1. Get uniref90.fasta
  cp $RDATADIR/\$UNIREFFASTA $CDATADIR
  gunzip $CDATADIR/\$UNIREFFASTA
  UNIREFFASTA=\${UNIREFFASTA%.gz}

  # 2. Get gene_association.goa_uniprot
  cp $RDATADIR/\$GOAUNIPROT $CDATADIR
  gunzip $CDATADIR/\$GOAUNIPROT
  GOAUNIPROT=\${GOAUNIPROT%.gz}

  # 3. Make BLASTDB fo uniref90.fasta
  ./makeblastdb -in $CDATADIR/\$UNIREFFASTA \\
   -dbtype prot -title uniref90 -input_type fasta \\
   -out $CDATADIR/uniref90

  # 4. BLAST protein sequences
  cp $RDATADIR/\$GOSEQFASTA $CDATADIR
  ./blastp -db $CDATADIR/uniref90 -query $CDATADIR/\$GOSEQFASTA \\
    -task blastp \\
    -outfmt 6 \\
    -num_threads $NUMBERCPU \\
    -evalue $BLASTEVALUE \\
    -out $CRUNANALYSIS/goseq.blast
  cp $CRUNANALYSIS/goseq.blast $RRUNANALYSIS

  # 5. Create a file that associates genes with gene ontology terms
  cp $RDATADIR/\$REF1GFF $CDATADIR
  perl pl/geneontology.pl gene2go \\
    -gff $DATADIR/\$REF1GFF \\
    -blast $CRUNANALYSIS/geneseq.blast \\
    -goa $CRUNANALYSIS/gene_association.goa_uniprot \\
    -out $CRUNANALYSIS/smutans.gene2go
  cp $CRUNANALYSIS/smutans.gene2go $RRUNANALYSIS

  # 6. Create a file that associates gene ontology terms with descriptions
  perl pl/geneontology.pl go2ngene \\
    -gene2go $CRUNANALYSIS/smutans.gene2go \\
    -obo $CRUNANALYSIS/gene_ontology.1_2.obo \\
    -out $CRUNANALYSIS/smutans.go2ngene
  cp $CRUNANALYSIS/smutans.go2ngene $RRUNANALYSIS

  # To see the output files in the compute nodes.
  tree
}
copy-data
process-data
cd
rm -rf \$TMPDIR
EOF
}

function goseq-rmessage {
  echo "bash $BASEDIR/goseq-push-data.sh"
  echo "work at cac:$CACWORKDIR"
  echo "bash $BASEDIR/goseq-get-data.sh"
  echo "bash $BASEDIR/job-feature"
  echo "bash $BASEDIR/job-goseq"
}
# END OF GOSEQ
################################################################################
