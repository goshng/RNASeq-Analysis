GENOMEDATADIR=/Volumes/Elements/Documents/Projects/mauve/bacteria

function build-index {
  Hpylori26695=$GEHOMEDATADIR/Helicobacter_pylori_26695_uid57787/NC_000915.fna
  segemehl.x -x output/NC_000915.idx -d $Hpylori26695
}


#####################################################################
# Main part of the script.
#####################################################################
PS3="Select what you want to do with rnaseq-analysis: "
CHOICES=( build-index )
select CHOICE in ${CHOICES[@]}; do 
 
  if [ "$CHOICE" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  elif [ "$CHOICE" == "build-index" ];  then
    build-index
    break
  elif [ "$CHOICE" == "generate-species" ];  then
    break
  elif [ "$CHOICE" == "preparation" ];  then
    break
  elif [ "$CHOICE" == "receive-run-mauve" ];  then
    break
  elif [ "$CHOICE" == "prepare-run-clonalframe" ];  then
    break
  elif [ "$CHOICE" == "compute-watterson-estimate-for-clonalframe" ];  then
    break
  elif [ "$CHOICE" == "receive-run-clonalframe" ];  then
    break
  elif [ "$CHOICE" == "prepare-run-clonalorigin" ];  then
    break
  elif [ "$CHOICE" == "receive-run-clonalorigin" ];  then
    break
  elif [ "$CHOICE" == "receive-run-2nd-clonalorigin" ];  then
    break
  else
    echo -e "You need to enter something\n"
    continue
  fi
done

