# Author: Sang Chul Choi
# Date  : 

function feature-genome {
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
      CHROMOSOME="gi|15644634|ref|NC_000915.1|"

      FEATURES=(gene genestart)
      PS3="What features do you wish to extract? "
      select FEATURE in ${FEATURES[@]}; do 
        if [ "$FEATURES" == "" ];  then
          echo -e "You need to enter something\n"
          continue
        else  
          OUT=$DATADIR/$FUNCNAME.out-$FEATURE
          echo "Coverting $REFGENOMEGFF to $OUT ..."
          perl pl/$FUNCNAME.pl \
            -gff $REFGENOMEGFF \
            -chromosome $CHROMOSOME \
            -feature $FEATURE\
            -out $OUT
          break
        fi
      done 
      break
    fi
  done
}
