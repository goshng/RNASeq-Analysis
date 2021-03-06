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

function feature-genome {
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

      REFGENOMEPTT=$(grep ^REFGENOMEPTT\: $SPECIESFILE | cut -d":" -f2)
      COMMAND1="perl pl/$FUNCNAME.pl ptt2 \
        -geneonly \
        -in $REFGENOMEPTT \
        -out $DATADIR/$FUNCNAME.out-geneonly"
      COMMAND1a="perl pl/$FUNCNAME.pl ptt2 \
        -withproduct \
        -geneonly \
        -in $REFGENOMEPTT \
        -out $DATADIR/$FUNCNAME.out-geneonly.txt"
      COMMAND2="perl pl/$FUNCNAME.pl ptt \
        -in $REFGENOMEPTT \
        -intergenicregion \
        -out $DATADIR/$FUNCNAME.out-intergenic"
      COMMAND3="perl pl/$FUNCNAME.pl ptt \
        -intergenicregiononly \
        -in $REFGENOMEPTT \
        -out $DATADIR/$FUNCNAME.out-intergeniconly"
      COMMAND4="perl pl/$FUNCNAME.pl ptt2 \
        -startcodon \
        -in $REFGENOMEPTT \
        -out $DATADIR/$FUNCNAME.out-startcodon"

      if [ "$BATCH" == "YES" ]; then
        echo $COMMAND1 >> $BATCHFILE
        echo $COMMAND1a >> $BATCHFILE
        #echo $COMMAND2 >> $BATCHFILE
        echo $COMMAND3 >> $BATCHFILE
        echo $COMMAND4 >> $BATCHFILE
      else
        echo $COMMAND1 | bash
        echo $COMMAND1a | bash
        #echo $COMMAND2 | bash
        echo $COMMAND3 | bash
        echo $COMMAND4 | bash
      fi

      echo -n "Do you wish to extract DNA sequences as well? (e.g., y/n) "
      read WISH
      if [ "$WISH" == "y" ]; then
        REFGENOMEFASTA=$(grep REFGENOMEFASTA $SPECIESFILE | cut -d":" -f2)
        COMMAND1="perl pl/$FUNCNAME.pl extract \
          -bed $DATADIR/$FUNCNAME.out-geneonly \
          -in $REFGENOMEFASTA \
          -out $DATADIR/geneonly.fa"
        COMMAND3="perl pl/$FUNCNAME.pl extract \
          -bed $DATADIR/$FUNCNAME.out-intergeniconly \
          -in $REFGENOMEFASTA \
          -out $DATADIR/intergeniconly.fa"
        COMMAND4="perl pl/$FUNCNAME.pl extract \
          -bed $DATADIR/$FUNCNAME.out-startcodon \
          -in $REFGENOMEFASTA \
          -out $DATADIR/startcodon.fa"
        echo $COMMAND1 | bash
        # echo $COMMAND2 | bash
        echo $COMMAND3 | bash
        echo $COMMAND4 | bash
      fi

      break

      # CHROMOSOME="gi|15644634|ref|NC_000915.1|"
      CHROMOSOME="gi|24378532|ref|NC_004350.1|"

      FEATURES=(gene genestart)
      PS3="What features do you wish to extract? "
      select FEATURE in ${FEATURES[@]}; do 
        if [ "$FEATURES" == "" ];  then
          echo -e "You need to enter something\n"
          continue
        else  
          OUT=$DATADIR/$FUNCNAME.out-$FEATURE
          echo "Coverting $REFGENOMEGFF to $OUT ..."
          perl pl/$FUNCNAME.pl \
            -gff $REFGENOMEGFF \
            -chromosome $CHROMOSOME \
            -feature $FEATURE\
            -out $OUT
          break
        fi
      done 
      break
    fi
  done
}
