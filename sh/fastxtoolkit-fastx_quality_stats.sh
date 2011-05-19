# Author: Sang Chul Choi
# Date  : 

function fastxtoolkit-fastx_quality_stats {
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

      in=$DATADIR/SRR03113S.fastq 
      out=$DATADIR/SRR03113S.quality
      fastx_quality_stats -i $in -o $out
      echo "Check $out"

      break
    fi
  done
}
