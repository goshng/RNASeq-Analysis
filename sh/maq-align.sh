###############################################################################
# Copyright (C) 2011 Sang Chul Choi
#
# This file is part of Mauve Analysis.
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

function maq-align {
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

      echo -n "Do you wish to run a batch? (e.g., y/n) "
      read WISH
      if [ "$WISH" == "y" ]; then
        BATCH=YES
        BATCHFILE=batch.sh
        echo "#!/bin/bash" > $BATCHFILE
      fi

      GENOMEFASTA=$(basename $REFGENOMEFASTA)
      GENOMEFASTABFA=$MAQDIR/$GENOMEFASTA.bfa
      COMMAND1="$MAQ fasta2bfa $REFGENOMEFASTA $GENOMEFASTABFA"
      if [ "$BATCH" == "YES" ]; then
        echo $COMMAND1 >> $BATCHFILE
      else
        echo $COMMAND1 | bash
        echo "Converting $GENOMEFASTA to $GENOMEFASTABFA ..."
      fi
      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      for g in $(eval echo {1..$NUMFASTQFILE}); do
        FASTQNUM=FASTQ$(printf "%02d" $g)
        GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
        FASTQFILE=$MAQDIR/$FASTQNUM.fastq
        MAQFASTQFILE=$MAQDIR/$FASTQNUM.maqfastq
        BFQFILE=$MAQDIR/$FASTQNUM.bfq
        COMMAND1="gunzip -c $GZIPFASTAQFILE > $FASTQFILE"
        COMMAND2="$MAQ sol2sanger $FASTQFILE $MAQFASTQFILE"
        COMMAND3="$MAQ fastq2bfq $MAQFASTQFILE $BFQFILE" 
        MAPFILE=$MAQDIR/$FASTQNUM.map
        COMMAND4="$MAQ match $MAPFILE $GENOMEFASTABFA $BFQFILE" 
        MAPCHECKFILE=$MAQDIR/$FASTQNUM.mapcheck
        COMMAND5="$MAQ mapcheck $GENOMEFASTABFA $MAPFILE > $MAPCHECKFILE"
        CNSFILE=$MAQDIR/$FASTQNUM.cns
        LOGFILE=$MAQDIR/$FASTQNUM.log
        COMMAND6="$MAQ assemble $CNSFILE $GENOMEFASTABFA $MAPFILE 2>$LOGFILE"
        CNSFQFILE=$MAQDIR/$FASTQNUM.cnsfq
        COMMAND7="$MAQ cns2fq $CNSFILE > $CNSFQFILE"

        COMMAND8="rm $FASTQFILE $GENOMEFASTABFA"
        COMMAND9="rm $MAQFASTQFILE"

        if [ "$BATCH" == "YES" ]; then
          echo $COMMAND1 >> $BATCHFILE
          echo $COMMAND2 >> $BATCHFILE
          echo $COMMAND3 >> $BATCHFILE
          echo $COMMAND4 >> $BATCHFILE
          echo $COMMAND5 >> $BATCHFILE
          echo $COMMAND6 >> $BATCHFILE
          echo $COMMAND7 >> $BATCHFILE
          echo $COMMAND8 >> $BATCHFILE
          echo $COMMAND9 >> $BATCHFILE
        else
          echo $COMMAND1 | bash
          echo $COMMAND2 | bash
          echo $COMMAND3 | bash
          echo $COMMAND4 | bash
          echo $COMMAND5 | bash
          echo $COMMAND6 | bash
          echo $COMMAND7 | bash
          echo $COMMAND8 | bash
          echo $COMMAND9 | bash
          echo "Check $MAQDIR/$FASTQNUM"
        fi
      done

      break
    fi
  done
}
