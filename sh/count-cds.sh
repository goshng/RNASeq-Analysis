###############################################################################
# Copyright (C) 2012 Sang Chul Choi
#
# This file is part of RNASeq Analysis.
# 
# RNASeq Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RNASeq Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RNASeq Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

function count-cds {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SPECIESS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -n "What repetition do you wish to run? (e.g., 1) "
      read REPETITION
      global-variable

      REFGENOMEID=$(grep ^REFGENOMEID\: $SPECIESFILE | cut -d":" -f2)
      REFGENOMEGFF=$(grep ^REFGENOMEGFF\: $SPECIESFILE | cut -d":" -f2)
      REFGENOMETXDB=$(grep ^REFGENOMETXDB\: $SPECIESFILE | cut -d":" -f2)
      REFGENOMETXDBBASE=$(basename $REFGENOMETXDB)

      RTEMP=$BWADIR/$RANDOM.R
      COMMAND="Rscript $RTEMP"
cat>$RTEMP<<EOF
library(GenomicFeatures)
txdbFile <- "$REFGENOMETXDB"
txdb <- loadFeatures(txdbFile)
count.table.file <- "$BWADIR/count-$REFGENOMEID.txt"
feature.cds <- cds(txdb,columns="tx_name")

count.table <- data.frame(gene=c("gene",unlist(elementMetadata(feature.cds)\$tx_name)))    
x <- read.table(file=count.table.file)                                          
count.table <- c("gene",count.table)                                            
x <- x[x\$V1 %in% count.table\$gene,]                                             
count.cds.file <- "$BWADIR/count-$REFGENOMEID.cds.txt"
write.table(x,file=count.cds.file,quote=FALSE,sep="\\t",row.names=FALSE,col.names=FALSE)
print(paste("Check file",count.cds.file))
EOF
      if [ "$BATCH" == "YES" ]; then
        echo $COMMAND >> $BATCHFILE
      else
        echo $COMMAND | bash
        rm $RTEMP
        echo $RTEMP
        echo "Edit and run $RTEMP"
	echo "e.g., Rscript $RTEMP"
      fi
      break
    fi
  done
}
