###############################################################################
# Copyright (C) 2011 Sang Chul Choi
#
# This file is part of RNASeq Analysis.
# 
# Mauve Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Mauve Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Mauve Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

function de-count {
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

      echo -n "Do you wish to count short reads? (e.g., y/n) "
      read WISH
      if [ "$WISH" == "y" ]; then
        echo -n "Do you wish to run a batch? (e.g., y/n) "
        read WISH
        if [ "$WISH" == "y" ]; then
          BATCH=YES
          BATCHFILE=batch.sh
          echo "#!/bin/bash" > $BATCHFILE
        fi
        NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
        for g in $(eval echo {1..$NUMFASTQFILE}); do
          FASTQNUM=FASTQ$(printf "%02d" $g)
          COMMAND1="perl pl/$FUNCNAME.pl join \
                    -shortread $BWADIR/$FASTQNUM-sum.pos \
                    -genepos $DATADIR/gene.pos \
                    -o $BWADIR/$FASTQNUM.de"

          if [ "$BATCH" == "YES" ]; then
            echo $COMMAND1 >> $BATCHFILE
          else
            echo $COMMAND1 | bash
            echo "Check $BWADIR/$FASTQNUM-sum files"
          fi
        done
      fi

      echo -n "Do you wish to merge de files? (e.g., y/n) "
      read WISH
      if [ "$WISH" == "y" ]; then
        # Create a gene count data file
        COUNTFILE=$BWADIR/count.txt
        COLNAME="0"
        printf "gene" > x
        cut -f1 $BWADIR/FASTQ01.de > 0
        # for g in $(eval echo {1..$NUMFASTQFILE}); do
        for g in 1 3 7 9 11; do
          FASTQNUM=FASTQ$(printf "%02d" $g)
          cut -f2 $BWADIR/$FASTQNUM.de > $g
          COLNAME="$COLNAME $g"
          printf "\t$g" >> x
        done
        printf "\n" >> x
        paste $COLNAME > y
        cat x y > $COUNTFILE
        rm $COLNAME x y
        echo "Check $COUNTFILE"
      fi

      break
    fi
  done

}
