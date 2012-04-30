###############################################################################
# Copyright (C) 2011-2012 Sang Chul Choi
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

function transcript-parsernaseq {
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
      REFGENOMEFASTA=$(grep REFGENOMEFASTA $SPECIESFILE | cut -d":" -f2)
      # for g in $(eval echo {1..$NUMFASTQFILE}); do
      # for g in $(eval echo {1..1}); do
# TW1-Glu/Gal Comparison   : 2,24,25     vs. 3,26,27
# UA159-Glu/Gal Comparison : 1,19,20     vs. 21,22,23,28
      for g in 1 2 3 19 20 21 22 23 24 25 26 27 28; do
        FASTQNUM=FASTQ$(printf "%03d" $g)

        COMMAND1="perl pl/$FUNCNAME.pl pileup \
          -wiggle $BWADIR/$FASTQNUM.wig \
          -out $BWADIR/$FASTQNUM.parsernaseq.pileup1"
        COMMAND2="perl pl/$FUNCNAME.pl gene \
          -feature $BWADIR/feature-genome.out-geneonly \
          -out $BWADIR/feature-genome.out-parsernaseq"

        # scp $BWADIR/$FASTQNUM.parsernaseq.pileup swiftgen:Documents/Projects/rnaseq/src/ParseRNAseq-1.1
        # scp $BWADIR/feature-genome.out-parsernaseq swiftgen:Documents/Projects/rnaseq/src/ParseRNAseq-1.1

        COMMAND3="ParseRNAseq \
          $BWADIR/$FASTQNUM.parsernaseq.pileup \
          $REFGENOMEFASTA \
          $BWADIR/feature-genome.out-parsernaseq \
          $BWADIR/$FASTQNUM.parsernaseq \
          -c 10 -b 25 -force_gp -fmt"
        COMMAND3a="perl pl/$FUNCNAME.pl remove \
          -feature $BWADIR/feature-genome.out-parsernaseq \
          -inpileup $BWADIR/$FASTQNUM.parsernaseq.pileup \
          -parsernaseq $BWADIR/$FASTQNUM.parsernaseq \
          -out $BWADIR/$FASTQNUM.parsernaseq.pileup2"
        COMMAND3b="perl pl/$FUNCNAME.pl remove \
          -inpileup $BWADIR/$FASTQNUM.parsernaseq.pileup2 \
          -parsernaseq $BWADIR/$FASTQNUM.parsernaseq2 \
          -out $BWADIR/$FASTQNUM.parsernaseq.pileup3"
        # scp $BWADIR/$FASTQNUM.parsernaseq.pileup3 swiftgen:Documents/Projects/rnaseq/src/ParseRNAseq-1.1
        # scp $BWADIR/$FASTQNUM.parsernaseq.pileup2 swiftgen:Documents/Projects/rnaseq/src/ParseRNAseq-1.1

        # scp swiftgen:Documents/Projects/rnaseq/src/ParseRNAseq-1.1/$FASTQNUM.parsernaseq $BWADIR
        # scp swiftgen:Documents/Projects/rnaseq/src/ParseRNAseq-1.1/$FASTQNUM.parsernaseq2 $BWADIR

        # scp swiftgen:Documents/Projects/rnaseq/src/ParseRNAseq-1.1/$FASTQNUM.bed2-* $BWADIR
        #scp swiftgen:Documents/Projects/rnaseq/src/ParseRNAseq-1.1/$FASTQNUM.parsernaseq? $BWADIR
        #scp swiftgen:Documents/Projects/rnaseq/src/ParseRNAseq-1.1/$FASTQNUM.parsernaseq?? $BWADIR

        PARSERNASEQOUT=$BWADIR/$FASTQNUM.parsernaseq1
        COMMAND4="perl pl/$FUNCNAME.pl bed \
          -parsernaseq $PARSERNASEQOUT \
          -out $BWADIR/$FASTQNUM.bed"
        COMMAND5="perl pl/$FUNCNAME.pl gff \
          -parsernaseq $PARSERNASEQOUT \
          -out $BWADIR/$FASTQNUM.bed2"
          # -toplevel \
        COMMAND5a="perl pl/$FUNCNAME.pl combine \
          -parsernaseq $PARSERNASEQOUT \
          -out $BWADIR/$FASTQNUM.bed2"
        echo $BWADIR/$FASTQNUM.bed2
        COMMAND6="perl pl/$FUNCNAME.pl operon \
          -feature $BWADIR/feature-genome.out-geneonly \
          -parsernaseq $PARSERNASEQOUT \
          -out $BWADIR/$FASTQNUM.operon"

        COMMAND7="perl pl/$FUNCNAME.pl adjust \
          -end $BWADIR/$FASTQNUM-end.wig \
          -operon $BWADIR/$FASTQNUM.operon \
          -out $BWADIR/$FASTQNUM.operon2"
        COMMAND8="perl pl/$FUNCNAME.pl slope \
          -windowsize 95 \
          -end $BWADIR/$FASTQNUM-end.wig \
          -operon $BWADIR/$FASTQNUM.operon \
          -out $BWADIR/$FASTQNUM.operon.slope"
        
        if [ "$BATCH" == "YES" ]; then
          echo $COMMAND1 >> $BATCHFILE
          echo $COMMAND2 >> $BATCHFILE
          echo $COMMAND3 >> $BATCHFILE
          echo $COMMAND3a >> $BATCHFILE
          echo $COMMAND3b >> $BATCHFILE
          echo $COMMAND4 >> $BATCHFILE
          echo $COMMAND5 >> $BATCHFILE
          echo $COMMAND5a >> $BATCHFILE
          echo $COMMAND6 >> $BATCHFILE
          echo $COMMAND7 >> $BATCHFILE
          echo $COMMAND8 >> $BATCHFILE

          # scp $BATCHFILE swiftgen:Documents/Projects/rnaseq/src/ParseRNAseq-1.1
        else
          echo $COMMAND1 | bash
          #echo $COMMAND1
          #echo $COMMAND2
          #echo $COMMAND3
          #echo $COMMAND3a
          #echo $COMMAND4
          #echo $COMMAND5 | bash
          #echo $COMMAND5a | bash
          #echo $COMMAND6
          #echo $COMMAND7
          #echo $COMMAND8
        fi
      done

      break
    fi
  done
}
