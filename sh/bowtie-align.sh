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

function bowtie-align {
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

      echo -n "Do you wish to run a batch? (e.g., y/n)"
      read WISH
      if [ "$WISH" == "y" ]; then
        BATCH=YES
        BATCHFILE=batch.sh
        echo "#!/bin/bash" > $BATCHFILE
      fi

      GENOMEFASTA=$(basename $REFGENOMEFASTA)
      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      for g in $(eval echo {1..$NUMFASTQFILE}); do
        FASTQNUM=FASTQ$(printf "%02d" $g)
        GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
        COMMAND1="zcat $GZIPFASTAQFILE | \
                 BOWTIE_INDEXES=$BOWTIEDIR \
                 $BOWTIE \
                 -t \
                 --solexa1.3-quals \
                 --threads $NUMBERCPU \
                 -S \
                 $GENOMEFASTA \
                 - $BOWTIEDIR/$FASTQNUM.sam"
        COMMAND2="$SAMTOOLS view -bS -o $BOWTIEDIR/$FASTQNUM.bam \
                  $BOWTIEDIR/$FASTQNUM.sam" 
        COMMAND3="$SAMTOOLS sort $BOWTIEDIR/$FASTQNUM.bam \
                  $BOWTIEDIR/$FASTQNUM.sorted"
        COMMAND4="$SAMTOOLS view $BOWTIEDIR/$FASTQNUM.sorted.bam \
                  | perl pl/sam2bed.pl > $BOWTIEDIR/$FASTQNUM.bed"

        if [ "$BATCH" == "YES" ]; then
          echo $COMMAND1 >> $BATCHFILE
          echo $COMMAND2 >> $BATCHFILE
          echo $COMMAND3 >> $BATCHFILE
          echo $COMMAND4 >> $BATCHFILE
        else
          echo $COMMAND1 | bash
          echo $COMMAND2 | bash
          echo $COMMAND3 | bash
          echo $COMMAND4 | bash
          echo "Check $BOWTIEDIR/$FASTQNUM"
        fi
      done

      break
    fi
  done
}
