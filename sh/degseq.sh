function degseq {
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
      RTEMP=$BWADIR/$RANDOM.R
      COMMAND="Rscript $RTEMP"
cat>$RTEMP<<EOF
library(DEGseq)
rnaseq1sample <- "$BWADIR/$FASTQ1ST.bed.degseq"
rnaseq2sample <- "$BWADIR/$FASTQ2ND.bed.degseq"
refFlat <- "$BASERUNANALYSIS/refFlat.txt"
mapResultBatch1 <- c(rnaseq1sample)
mapResultBatch2 <- c(rnaseq2sample)
DEGseq(mapResultBatch1, mapResultBatch2, fileFormat="bed", refFlat=refFlat, outputDir="$DEGSEQOUT", method="LRT")
EOF
      if [ "$BATCH" == "YES" ]; then
        echo $COMMAND >> $BATCHFILE
      else
        $COMMAND
        rm $RTEMP
        echo "Check $DEGSEQOUT"
      fi
      break
    fi
  done
}
