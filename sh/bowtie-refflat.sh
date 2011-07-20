function bowtie-refflat {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SPECIESS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      BASERUNANALYSIS=output/$SPECIES/run-analysis
      echo "You need a refFlat file to convert bowtie mapped reads to"
      echo "expression data that DEGseq can process."
      echo "Create a file with chrom, chromStart, chromEnd, name,"
      echo "and strand using a UCSC genome browser or other method"
      echo "you would like. Go to strep-genome.bscb.cornell.edu."
      echo "Clock Tables. Select S. mutans UA159. Select Genes and Gene"
      echo "Prediction Tracks. Select Known Genes. Then, click get output"
      echo "Choose chrom, chromStart, chromEnd, name, and strand."
      echo "Click get output."
      echo "Note that you have to make sure that"
      echo "genomic positions are 0-based or 1-based."
      echo -n "What is the file name? (e.g., data/hgTables.txt) "
      read HGTABLES
      # HGTABLES=hgTables.txt
      for i in {1..5}; do
        cut -f$i $HGTABLES > $i
      done
      awk '{s=$1+1;print s}' 2 > 22
      mv 22 2
      awk '{print $0","}' 2 > 22
      awk '{print $0","}' 3 > 33
      awk '{print 1}' 1 > 11
      paste 4 4 1 5 2 3 2 3 11 22 33 > x
      sed '1d;$d' x > $BASERUNANALYSIS/refFlat.txt
      rm 1 2 3 4 5 11 22 33 x
      echo "Check file $BASERUNANALYSIS/refFlat.txt"
      break
    fi
  done
}
