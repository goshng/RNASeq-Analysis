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

      echo -n "Do you wish to trim RNA-Seq (trimmed)? (e.g., y/n) "
      read WISH
      if [ "$WISH" == "y" ]; then
        NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
        for g in $(eval echo {1..$NUMFASTQFILE}); do
          FASTQNUM=FASTQ$(printf "%02d" $g)
          ADAPTERNUM=ADAPTER$(printf "%02d" $g)
          ADAPTERSEQUENCE=$(grep $ADAPTERNUM $SPECIESFILE | cut -d":" -f2)
          GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
          GZIPFASTAQDIR=${GZIPFASTAQFILE%/*}
          OUTFILE=$DATADIR/$FASTQNUM.trimmed.gz
          echo cutadapt -a $ADAPTERSEQUENCE -o $OUTFILE $GZIPFASTAQFILE
          cutadapt -a $ADAPTERSEQUENCE -o $OUTFILE $GZIPFASTAQFILE
          mv $OUTFILE $GZIPFASTAQDIR
          echo "Check $GZIPFASTAQDIR/$FASTQNUM.trimmed.gz"
        done
        echo "Change raw sequence file names in $GZIPFASTAQDIR"
      fi

      echo -n "Do you wish to sample RNA-Seq (sample)? (e.g., y/n) "
      read WISH
      if [ "$WISH" == "y" ]; then
        echo -n "What proportion do you want to sample? (e.g., 35 for 35%) "
        read PROPORTION
        NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
        for g in $(eval echo {1..$NUMFASTQFILE}); do
          FASTQNUM=FASTQ$(printf "%02d" $g)
          GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
          GZIPFASTAQDIR=${GZIPFASTAQFILE%/*}
          OUTFILE=$DATADIR/$FASTQNUM.sample-$PROPORTION
          perl pl/$FUNCNAME.pl sample \
            --fastq $GZIPFASTAQFILE \
            --out $OUTFILE \
            --proportion $PROPORTION
          echo "gzipping $OUTFILE"
          gzip $OUTFILE
          mv $OUTFILE.gz $GZIPFASTAQDIR
          echo "Check $GZIPFASTAQDIR/$FASTQNUM.sample-$PROPORTION.gz"
        done
      fi

      echo -n "Do you wish to cut RNA-Seq (cut)? (e.g., y/n) "
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
          echo "gzipping $OUTFILE"
          gzip $OUTFILE
          mv $OUTFILE.gz $GZIPFASTAQDIR
          echo "Check $GZIPFASTAQDIR/$FASTQNUM.cut-$CUTSIZE.gz"
        done
      fi

      echo -n "Do you wish to keep longer RNA-Seq (keep)? (e.g., y/n) "
      read WISH
      if [ "$WISH" == "y" ]; then
        echo -n "What is minimum length of short reads? (e.g., 80) "
        read KEEPSIZE
        NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
        for g in $(eval echo {1..$NUMFASTQFILE}); do
          FASTQNUM=FASTQ$(printf "%02d" $g)
          GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
          GZIPFASTAQDIR=${GZIPFASTAQFILE%/*}
          OUTFILE=$DATADIR/$FASTQNUM.keep-$KEEPSIZE
          perl pl/$FUNCNAME.pl keep \
            --fastq $GZIPFASTAQFILE \
            --out $OUTFILE \
            --keepsize $KEEPSIZE
          gzip $OUTFILE
          echo "gzipping $OUTFILE"
          mv $OUTFILE.gz $GZIPFASTAQDIR
          echo "Check $GZIPFASTAQDIR/$FASTQNUM.keep-$KEEPSIZE.gz"
        done
      fi

      echo -n "Do you wish to remove longer RNA-Seq (remove)? (e.g., y/n) "
      read WISH
      if [ "$WISH" == "y" ]; then
        echo -n "What is minimum length of short reads to remove? (e.g., 80) "
        read REMOVESIZE
        NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
        for g in $(eval echo {1..$NUMFASTQFILE}); do
          FASTQNUM=FASTQ$(printf "%02d" $g)
          GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
          GZIPFASTAQDIR=${GZIPFASTAQFILE%/*}
          OUTFILE=$DATADIR/$FASTQNUM.remove-$REMOVESIZE
          perl pl/$FUNCNAME.pl remove \
            --fastq $GZIPFASTAQFILE \
            --out $OUTFILE \
            --removesize $REMOVESIZE
          gzip $OUTFILE
          echo "gzipping $OUTFILE"
          mv $OUTFILE.gz $GZIPFASTAQDIR
          echo "Check $GZIPFASTAQDIR/$FASTQNUM.remove-$REMOVESIZE.gz"
        done
      fi

      echo -n "Do you wish to keep the first some bases of RNA-Seq (max)? (e.g., y/n) "
      read WISH
      if [ "$WISH" == "y" ]; then
        echo -n "What is maximum length of short reads? (e.g., 15) "
        read MAXSIZE
        NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
        for g in $(eval echo {1..$NUMFASTQFILE}); do
          FASTQNUM=FASTQ$(printf "%02d" $g)
          GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
          GZIPFASTAQDIR=${GZIPFASTAQFILE%/*}
          OUTFILE=$DATADIR/$FASTQNUM.max-$MAXSIZE
          perl pl/$FUNCNAME.pl max \
            --fastq $GZIPFASTAQFILE \
            --out $OUTFILE \
            --maxsize $MAXSIZE
          gzip $OUTFILE
          echo "gzipping $OUTFILE"
          mv $OUTFILE.gz $GZIPFASTAQDIR
          echo "Check $GZIPFASTAQDIR/$FASTQNUM.max-$MAXSIZE.gz"
        done
      fi


      break
    fi
  done

}
