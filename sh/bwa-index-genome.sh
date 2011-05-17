# Author: Sang Chul Choi
# Date  : Tue May 17 17:48:48 EDT 2011

function bwa-index-genome {
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

      $BWA index $REFGENOMEFASTA

      break
    fi
  done

}
