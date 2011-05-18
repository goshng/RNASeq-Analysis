# Author: Sang Chul Choi
# Date  : 

function segemehl-index-genome {
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

      GENOMEFASTA=$(basename $REFGENOMEFASTA)
      $SEGEMEHL -x $DATADIR/$GENOMEFASTA-segemehl.idx -d $DATADIR/$GENOMEFASTA

      break
    fi
  done

}
