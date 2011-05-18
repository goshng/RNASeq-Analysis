# Author: Sang Chul Choi
# Date  : Tue May 17 22:47:18 EDT 2011

function choose-species {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SPECIESS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -n "What repetition do you wish to run? (e.g., 1) "
      read REPETITION
      global-variable $SPECIES $REPETITION

      mkdir -p $NUMBERDIR
      mkdir -p $DATADIR
      mkdir -p $BWADIR

      break
    fi
  done

}
