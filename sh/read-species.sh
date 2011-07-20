function read-species {
  SPECIESFILE=species/$SPECIES
  # echo -n "  Reading REFGENOMEFASTA from $SPECIESFILE..."
  REFGENOMEFASTA=$(grep REFGENOMEFASTA $SPECIESFILE | cut -d":" -f2)
  # echo " $REFGENOMEFASTA"
  # echo -n "  Reading REFGENOMEGFF from $SPECIESFILE..."
  REFGENOMEGFF=$(grep REFGENOMEGFF $SPECIESFILE | cut -d":" -f2)
  # echo " $REFGENOMEGFF"
}
