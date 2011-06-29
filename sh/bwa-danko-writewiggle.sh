# Author: Sang Chul Choi
# Date  : 

function bwa-danko-writewiggle {
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

      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      for g in $(eval echo {1..$NUMFASTQFILE}); do
        FASTQNUM=FASTQ$(printf "%02d" $g)
        COMMAND1="Rscript $DATADIR/$FASTQNUM.R"
        # COMMAND2="sed s/gi\\\|24378532\\\|ref\\\|NC_004350.1\\\|/chr1/ < $DATADIR/$FASTQNUM.wig > $DATADIR/$FASTQNUM.tmp"
        COMMAND3="mv $DATADIR/$FASTQNUM.tmp $DATADIR/$FASTQNUM.wig"

cat>$DATADIR/$FASTQNUM.R<<EOF
require(GROseq)
load ("$DATADIR/$FASTQNUM.RData")
# WriteWiggle(p=data[,c(1:3,6)], file="$DATADIR/$FASTQNUM", size=5, reverse=FALSE, str="+", debug=FALSE, track.type.line=TRUE)
WriteWiggle(p=data[,c(1:3,6)], file="$DATADIR/$FASTQNUM", size=5, reverse=FALSE, debug=FALSE, track.type.line=TRUE)
EOF

        if [ "$BATCH" == "YES" ]; then
          echo $COMMAND1 >> $BATCHFILE
          echo $COMMAND2 >> $BATCHFILE
          echo $COMMAND3 >> $BATCHFILE
        else
          # fixedStep chrom=gi|24378532|ref|NC_004350.1| start=1 step=5 span=5
          $COMMAND1
          sed s/gi\\\|24378532\\\|ref\\\|NC_004350.1\\\|/chr1/ < $DATADIR/$FASTQNUM.wig > $DATADIR/$FASTQNUM.tmp
          $COMMAND3
          scp $DATADIR/$FASTQNUM.wig $X11_LOGIN:public_html/
          echo "Check $DATADIR/$FASTQNUM.wig"
          echo "Check at http://strep-genome.bscb.cornell.edu/cgi-bin/hgTracks?db=SmuUA159&position=chr1&hgt.customText=http://compgen.bscb.cornell.edu/~choi/$FASTQNUM.wig"
        fi
      done

      break
    fi
  done
}
