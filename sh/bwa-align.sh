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

function bwa-align {
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
      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      for g in $(eval echo {1..$NUMFASTQFILE}); do
        FASTQNUM=FASTQ$(printf "%02d" $g)
        GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
        COMMAND1="$BWA aln -I -t $NUMBERCPU \
                  $BWADIR/$GENOMEFASTA-bwa \
                  $GZIPFASTAQFILE > $BWADIR/$FASTQNUM.sai"
        COMMAND2="$BWA samse -n 1 \
                  -f $BWADIR/$FASTQNUM.sam \
                  $BWADIR/$GENOMEFASTA-bwa \
                  $BWADIR/$FASTQNUM.sai \
                  $GZIPFASTAQFILE"
        COMMAND3="$SAMTOOLS view -bS -o $BWADIR/$FASTQNUM.bam \
                  $BWADIR/$FASTQNUM.sam"
        COMMAND4="$SAMTOOLS sort $BWADIR/$FASTQNUM.bam \
                  $BWADIR/$FASTQNUM.sorted"
        DELETE1="rm $BWADIR/$FASTQNUM.sai \
                    $BWADIR/$FASTQNUM.sam \
                    $BWADIR/$FASTQNUM.bam"
        if [ "$BATCH" == "YES" ]; then
          echo $COMMAND1 >> $BATCHFILE
          echo $COMMAND2 >> $BATCHFILE
          echo $COMMAND3 >> $BATCHFILE
          echo $COMMAND4 >> $BATCHFILE
          echo $DELETE1 >> $BATCHFILE
        else
          echo $COMMAND1 | bash
          echo $COMMAND2 | bash
          echo $COMMAND3 | bash
          echo $COMMAND4 | bash
          echo $DELETE1 | bash
          echo "Check $BWADIR/$FASTQNUM.sorted.bam"
          echo "Use $SAMTOOLS view $BWADIR/$FASTQNUM.sorted.bam"
          echo "to view the alignment."
        fi
      done
      break
    fi
  done
}

function bwa-align-per-fastq {
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

      echo -n "What FASTQ file do you wish to run? (e.g., 1) "
      read g

      GENOMEFASTA=$(basename $REFGENOMEFASTA)
      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      FASTQNUM=FASTQ$(printf "%02d" $g)
      GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
      GZIPFASTAQDIR=${GZIPFASTAQFILE%/*}

      OUTFILE=input/$FASTQNUM.gz
      COMMANDFILTER="zcat $GZIPFASTAQFILE | \
                     grep -A 3 '^@.* [^:]*:N:[^:]*:' | \
                     sed '/^--$/d' | \
                     gzip > $OUTFILE"
      GZIPFASTAQFILE=$OUTFILE

      ADAPTERNUM=ADAPTER$(printf "%02d" $g)
      ADAPTERSEQUENCE=$(grep ^$ADAPTERNUM $SPECIESFILE | cut -d":" -f2)
      OUTFILE=input/$FASTQNUM.cutadapt.gz
      COMMANDCUTADAPT="cutadapt --minimum-length=25 -a $ADAPTERSEQUENCE -o \
                       $OUTFILE $GZIPFASTAQFILE"
      COMMANDCUTADAPT2="mv $OUTFILE $GZIPFASTAQFILE"

      # COMMAND1="$BWA aln -I -t $NUMBERCPU \
      COMMAND1="$BWA aln -t $NUMBERCPU \
                $BWADIR/$GENOMEFASTA-bwa \
                $GZIPFASTAQFILE > $BWADIR/$FASTQNUM.sai"
      COMMAND2="$BWA samse -n 1 \
                -f $BWADIR/$FASTQNUM.sam \
                $BWADIR/$GENOMEFASTA-bwa \
                $BWADIR/$FASTQNUM.sai \
                $GZIPFASTAQFILE"
      COMMAND3="$SAMTOOLS view -bS -o $BWADIR/$FASTQNUM.bam \
                $BWADIR/$FASTQNUM.sam"
      COMMAND4="$SAMTOOLS sort $BWADIR/$FASTQNUM.bam \
                $BWADIR/$FASTQNUM.sorted"
      DELETE1="rm $BWADIR/$FASTQNUM.sai \
                  $BWADIR/$FASTQNUM.sam \
                  $BWADIR/$FASTQNUM.bam"
      BATCH=YES
      if [ "$BATCH" == "YES" ]; then
        BATCHFILE=batch.sh
        echo "#!/bin/bash" > $BATCHFILE
        echo $COMMANDFILTER >> $BATCHFILE
        echo $COMMANDCUTADAPT >> $BATCHFILE
        echo $COMMANDCUTADAPT2 >> $BATCHFILE
        echo $COMMAND1 >> $BATCHFILE
        echo $COMMAND2 >> $BATCHFILE
        echo $COMMAND3 >> $BATCHFILE
        echo $COMMAND4 >> $BATCHFILE
        echo $DELETE1 >> $BATCHFILE
      else
        $COMMANDFILTER | bash
        $COMMANDCUTADAPT | bash
        $COMMANDCUTADAPT2 | bash
        $COMMAND1 | bash
        $COMMAND2 | bash
        $COMMAND3 | bash
        $COMMAND4 | bash
        $DELETE1 | bash
        echo "Check $BWADIR/$FASTQNUM.sorted.bam"
        echo "Use $SAMTOOLS view $BWADIR/$FASTQNUM.sorted.bam"
        echo "to view the alignment."
      fi
      break
    fi
  done
}
