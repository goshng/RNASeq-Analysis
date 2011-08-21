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

function fastq-sample {
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

      echo -n "Do you wish to sample RNA-Seq? (e.g., y/n) "
      read WISH
      if [ "$WISH" == "y" ]; then
        echo -n "In what intervals do you want to sample? (e.g., 100 for 1%) "
        read INTERVAL
        NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
        for g in $(eval echo {1..$NUMFASTQFILE}); do
          FASTQNUM=FASTQ$(printf "%02d" $g)
          GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
          GZIPFASTAQDIR=${GZIPFASTAQFILE%/*}
          OUTFILE=$DATADIR/$FASTQNUM.subsample-$INTERVAL
          perl pl/$FUNCNAME.pl sample \
            --fastq $GZIPFASTAQFILE \
            --out $OUTFILE \
            --interval $INTERVAL
          gzip $OUTFILE
          mv $OUTFILE.gz $GZIPFASTAQDIR
          echo "Check $GZIPFASTAQDIR/$FASTQNUM.subsample-$INTERVAL.gz"
        done
      fi

      echo -n "Do you wish to cut RNA-Seq? (e.g., y/n) "
      read WISH
      if [ "$WISH" == "y" ]; then
        echo -n "How many nucleotides do you wish to cut? (e.g., 65) "
        read CUTSIZE
        NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
        for g in $(eval echo {1..$NUMFASTQFILE}); do
          FASTQNUM=FASTQ$(printf "%02d" $g)
          GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
          GZIPFASTAQDIR=${GZIPFASTAQFILE%/*}
          OUTFILE=$DATADIR/$FASTQNUM.cut-$CUTSIZE
          perl pl/$FUNCNAME.pl cut \
            --fastq $GZIPFASTAQFILE \
            --out $OUTFILE \
            --cutsize $CUTSIZE
          gzip $OUTFILE
          mv $OUTFILE.gz $GZIPFASTAQDIR
          echo "Check $GZIPFASTAQDIR/$FASTQNUM.cut-$CUTSIZE.gz"
        done
      fi

      break
    fi
  done

}
