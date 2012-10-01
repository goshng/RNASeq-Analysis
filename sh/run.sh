#!/bin/bash
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

#####################################################################
# External shell scripts
#####################################################################
source sh/copyright.sh
source sh/utility.sh
source sh/conf.sh
source sh/global-variable.sh
source sh/read-species.sh
source sh/init-file-system.sh
source sh/choose-species.sh
source sh/cp-genome.sh
source sh/bwa-index-genome.sh
source sh/bwa-align.sh
source sh/bwa-samse.sh
source sh/bwa-samtools-view.sh
source sh/bwa-samtools-sort.sh
source sh/bwa-samtools-bed.sh
source sh/bwa-degseq-bed.sh
source sh/bwa-batch-align-to-bed.sh
source sh/bwa-samtools-wig.sh
source sh/fastq-summary.sh

source sh/bowtie-index-genome.sh
source sh/bowtie-align.sh
source sh/bowtie-refflat.sh

source sh/bwa-R-saveimage.sh
source sh/bwa-danko-writewiggle.sh
source sh/bwa-danko-countreadsininterval.sh
source sh/bwa-danko-metagene.sh
source sh/bwa-danko-detecttranscriptsem.sh
source sh/segemehl-index-genome.sh
source sh/feature-genome.sh
source sh/fastxtoolkit-fastq_to_fasta.sh
source sh/fastxtoolkit-fastx_quality_stats.sh

source sh/fastq-sample.sh
source sh/degseq.sh
source sh/edgeR.sh
source sh/transcript-cufflinks.sh
source sh/transcript-genecoverage.sh
source sh/transcript-parsernaseq.sh
source sh/maq-align.sh
source sh/bwa-summary.sh
source sh/bwa-pos2wig.sh
source sh/de-count.sh
source sh/bwa-mpileup.sh
source sh/transcript-summary.sh
source sh/job-truncate-reads.sh 
source sh/job-truncate-reads-local.sh 
source sh/convert-gff2txdb.sh
source sh/count-cds.sh
source sh/batch.sh
source sh/batch2.sh
source sh/batch3.sh
source sh/deseq.sh
source sh/goseq.sh
source sh/mcl.sh
source sh/sara-batch.sh

#####################################################################
# Read configuration file
#####################################################################
conf

#####################################################################
# Read directories
#####################################################################
SPECIESS=$(ls $ROOTANALYSISDIR/species|grep -v ^sim)

#####################################################################
# Menus
#####################################################################
PS3="Select the menu : "
CHOICES=( choose-species \
          convert-gff2txdb \
          batch3 \
          deseq \
          goseq \
          mcl \
          fastq-sample \
          sara-batch \
          warranty \
          copyright \
          quit )
if [ $# -eq 1 ]; then
  CHOICEID=$(($1 - 1))
  CHOICE=${CHOICES[$CHOICEID]}
  $CHOICE
else
  select CHOICE in ${CHOICES[@]}; do 
    if [ "$CHOICE" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else
      $CHOICE
      break
    fi
  done
fi
