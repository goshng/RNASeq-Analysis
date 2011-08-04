# Author: Sang Chul Choi

function fastq-summary {
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

      BATCH=YES
      BATCHFILE=batch.sh

      GENOMEFASTA=$(basename $REFGENOMEFASTA)
      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      for g in $(eval echo {1..$NUMFASTQFILE}); do
        FASTQNUM=FASTQ$(printf "%02d" $g)
        GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)

        echo "FASTQ: $FASTQNUM"
        echo "FILE: $GZIPFASTAQFILE"
        TOTALLENGTH=$(zcat $GZIPFASTAQFILE|wc -l)
        TOTALLENGTH=$(($TOTALLENGTH / 4))
        echo "The total number of short reads is $TOTALLENGTH"

        NUMBERMAPPEDREAD=$(trim $(cat $BWADIR/$FASTQNUM.bed|wc -l))
        PERCENTMAPPEDREAD=$(($NUMBERMAPPEDREAD * 100 / $TOTALLENGTH))
        echo "The number of mapped reads using BWA is $NUMBERMAPPEDREAD ($PERCENTMAPPEDREAD%)"
        NUMBERMAPPEDREAD=$(trim $(cat $BOWTIEDIR/$FASTQNUM.bed|wc -l))
        PERCENTMAPPEDREAD=$(($NUMBERMAPPEDREAD * 100 / $TOTALLENGTH))
        echo "The number of mapped reads using Bowtie is $NUMBERMAPPEDREAD ($PERCENTMAPPEDREAD%)"
      done

      break
    fi
  done

}
