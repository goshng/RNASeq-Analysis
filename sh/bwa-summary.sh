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

function bwa-summary {
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

      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      for g in $(eval echo {1..$NUMFASTQFILE}); do
        FASTQNUM=FASTQ$(printf "%02d" $g)

        COMMAND2="$SAMTOOLS view $BWADIR/$FASTQNUM.sorted.bam \
          | perl pl/bwa-sam.pl parse > $BWADIR/$FASTQNUM-sum.parse"

        COMMAND3="$SAMTOOLS view $BWADIR/$FASTQNUM.sorted.bam \
          | perl pl/bwa-sam.pl pos > $BWADIR/$FASTQNUM-sum.pos"

        if [ "$BATCH" == "YES" ]; then
          echo $COMMAND2 >> $BATCHFILE
          #echo $COMMAND3 >> $BATCHFILE
        else
          echo $COMMAND2 | bash
          #echo $COMMAND3 | bash
          echo "Check $BWADIR/$FASTQNUM-sum files"
        fi
        # Find the total number of reads and mapped reads.
        SUMMARY=$BWADIR/$FASTQNUM.sum
        GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
        echo "FASTQ: $FASTQNUM" > $SUMMARY
        echo "FILE: $GZIPFASTAQFILE" >> $SUMMARY
        TOTALLENGTH=$(zcat $GZIPFASTAQFILE|wc -l)
        TOTALLENGTH=$(($TOTALLENGTH / 4))
        echo "The total number of short reads is $TOTALLENGTH" >> $SUMMARY
        NUMBERBAM=$(trim $($SAMTOOLS view $BWADIR/$FASTQNUM.sorted.bam|wc -l))
        echo "The total number of short reads in sorted.bam is $NUMBERBAM" >> $SUMMARY
        NUMBERBAM=$(trim $(wc -l $BWADIR/$FASTQNUM-sum.pos))
        echo "The total number of short reads in mapped is $NUMBERBAM" >> $SUMMARY

        echo "Check $SUMMARY"
      done

      break
    fi
  done
}