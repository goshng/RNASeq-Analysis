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

# Usages:
# perl pl/transcript-summary.pl summary -subcmd unannotated -feature /Users/goshng/Documents/Projects/rnaseq/output/smutans12/1/data/feature.pos -transcript /Users/goshng/Documents/Projects/rnaseq/output/smutans12/1/transcript/FASTQ01.bed > 1.unannotated
# perl pl/transcript-summary.pl getsequence -in 1 -col 2 -size 50 -fasta /Users/goshng/Documents/Projects/rnaseq/output/smutans12/1/data/NC_004350.fna > 2

function transcript-genecoverage {
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

      REFGENOMELENGTH=$(grep REFGENOMELENGTH $SPECIESFILE | cut -d":" -f2)
      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      # for g in $(eval echo {1..$NUMFASTQFILE}); do
      for g in $(eval echo {1..1}); do
        FASTQNUM=FASTQ$(printf "%02d" $g)

        COMMAND1="perl pl/$FUNCNAME.pl coverage \
          -feature $DATADIR/feature-genome.out-geneonly \
          -wiggle $BWADIR/$FASTQNUM.wig"
  
        if [ "$BATCH" == "YES" ]; then
          echo $COMMAND1 >> $BATCHFILE
        else
          # echo $COMMAND1 | bash
          echo $COMMAND1
          echo "Check $BWADIR/$FASTQNUM-sum files"
        fi
      done

      break
    fi
  done
}
