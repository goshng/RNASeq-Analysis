function edgeR {
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
      echo "Check file $BASERUNANALYSIS/refFlat.txt"
      echo -n "Enter the first RNA-Seq sample: "
      read FIRST
      FASTQ1ST=FASTQ$(printf "%02d" $FIRST)
      echo -n "Enter the second RNA-Seq sample: "
      read SECOND
      FASTQ2ND=FASTQ$(printf "%02d" $SECOND)
      DEGSEQOUT=$BWADIR/${FIRST}vs${SECOND}

      # Check $DEGSEQOUT directory.
      perl pl/edgeR.pl --degseqout $DEGSEQOUT

      RTEMP=$DEGSEQOUT/$RANDOM.R
      ROUT=$DEGSEQOUT/edgeR.out
      COMMAND="Rscript $RTEMP > $ROUT"
cat>$RTEMP<<EOF
library(edgeR)
setwd("$DEGSEQOUT")
targets <- read.delim(file = "datalist", stringsAsFactors = FALSE)
d <- readDGE(targets)
colnames(d) <- c("$FIRST", "$SECOND")
d <- estimateCommonDisp(d)
de.common <- exactTest(d)
topTags(de.common)
summary(decideTestsDGE(de.common, p.value = 0.05))
detags500.com <- rownames(topTags(de.common, n = 500)$table)
EOF
      if [ "$BATCH" == "YES" ]; then
        echo $COMMAND >> $BATCHFILE
      else
        echo $COMMAND | bash
        # rm $RTEMP
        echo "Check $DEGSEQOUT"
      fi
      break
    fi
  done
}
