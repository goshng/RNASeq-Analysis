# Author: Sang Chul Choi
# Date  : 

function bwa-samtools-view {
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

      # bwa samse [-n maxOcc] <in.db.fasta> <in.sai> <in.fq> > <out.sam>
      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      for g in $(eval echo {1..$NUMFASTQFILE}); do
        FASTQNUM=FASTQ$(printf "%02d" $g)
        GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
        $SAMTOOLS view -b -S $DATADIR/$FASTQNUM.sam > $DATADIR/$FASTQNUM.bam
        echo "Check $DATADIR/$FASTQNUM.bam"
      done

      break
    fi
  done

}
