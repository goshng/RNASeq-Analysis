###############################################################################
# Copyright (C) 2011 Sang Chul Choi
#
# This file is part of RNASeq Analysis.
# 
# Mauve Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Mauve Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Mauve Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

function deseq {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SPECIESS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -n "What repetition do you wish to run? (e.g., 1) "
      read REPETITION
      global-variable $SPECIES $REPETITION
      read-species

      echo -n "What is the count file name? (e.g., count-1-11) "
      read DESEQIN
      echo -n "Enter the first sample? (e.g., UA159) "
      read FRISTSAMPLE
      echo -n "Enter the second sample? (e.g., TW1Glucose) "
      read SECONDSAMPLE
      # DESEQIN=count

      RTEMP=$BWADIR/$RANDOM.R
      OUT=$BWADIR/$DESEQIN.out
      COMMAND="Rscript $RTEMP > $OUT"
cat>$RTEMP<<EOF
library(DESeq)
countsFile <- "$BWADIR/$DESEQIN.txt"
countsTable <- read.delim (countsFile, header=TRUE, stringsAsFactors=TRUE)
rownames(countsTable) <- countsTable\$gene
countsTable <- countsTable[,-1]
countsTable <- countsTable[countsTable[,1]+countsTable[,2]>5,]
conds <- c("$FRISTSAMPLE", "$SECONDSAMPLE")
cds <- newCountDataSet(countsTable, conds)
cds <- estimateSizeFactors(cds)
cds <- estimateVarianceFunctions(cds,method="blind")
res <- nbinomTest(cds, "$FRISTSAMPLE", "$SECONDSAMPLE")
      

plotDE <- function( res )
{
   plot(res\$baseMean, res\$log2FoldChange, log="x", pch=20, cex=.5, 
        col = ifelse( res\$padj < .1, "red", "black" ) )
}

cat ("SCV file: $DESEQIN-de.ps\n")
postscript ("$BWADIR/$DESEQIN-de.ps",  width=10, height=10, horizontal = FALSE, onefile = FALSE,
paper = "special")
plotDE(res)
y <- dev.off()

cat ("SCV file: $DESEQIN-scv.ps\n")
postscript ("$BWADIR/$DESEQIN-scv.ps",  width=10, height=10, horizontal = FALSE, onefile = FALSE,
paper = "special")
scvPlot(cds)
y <- dev.off()

resSig <- res[res\$padj < .1,]
y <- length(countsTable[,1])
cat ("Total number of genes is ", y, "\n",sep="")
y <- length(resSig[,1])
cat ("The number of differentially expressed genes is ", y, "\n",sep="")
y <- length(res[res\$baseMean < 100,1])
cat ("The number of genes with less than 100 of mean counts is ", y, "\n",sep="")
cat ("\nList of differentially expressed genes\n")
options(width = 1000) 
print(resSig[order(resSig\$pval),])
cat ("\nThe same table as above but in order of down-regulated first of the differentially expressed genes\n")
print(resSig[order(resSig\$foldChange,-resSig\$baseMean),])
cat ("\nThe same table as above but in order of up-regulated first of the differentially expressed genes\n")
print(resSig[order(-resSig\$foldChange,-resSig\$baseMean),])

EOF
      if [ "$BATCH" == "YES" ]; then
        echo $COMMAND >> $BATCHFILE
      else
        #echo $COMMAND | bash
        #rm $RTEMP
        #echo $RTEMP
        echo "Edit and run $BWADIR/$RTEMP"
	echo "e.g., Rscript $RTEMP > $DESEQIN.out"
      fi
      break
    fi
  done
}
