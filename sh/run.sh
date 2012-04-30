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

ROOTANALYSISDIR=`pwd`

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
source sh/deseq.sh
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
source sh/batch.sh
source sh/batch2.sh

#####################################################################
# Read configuration file
#####################################################################
conf

#####################################################################
# Read directories
#####################################################################
SPECIESS=$(ls species|grep -v ^sim)

#####################################################################
# Menus
#####################################################################
PS3="Select the menu : "
CHOICES=( init-file-system \
          choose-species \
          ---FASTQ---\
          fastq-summary \
          fastq-sample \
          ---BATCH---\
	  batch \
	  batch2 \
          ---BWA-ALIGN---\
          cp-genome \
          bwa-index-genome \
          bwa-align \
          bwa-align-per-fastq \
          bwa-mpileup \
          ---BWA-OTHER---\
          bwa-samtools-wig \
          bwa-summary \
          bwa-pos2wig \
          ---DE---\
          feature-genome \
          ---TRANSCRIPT---\
          transcript-genecoverage \
          transcript-parsernaseq \
          transcript-summary \
          transcript-cufflinks \
          ---SMALL-RNA---\
          rna-rnaz \
          ---READ-QUALITY-CHECK--- \
          warranty \
          copyright \
          quit )
select CHOICE in ${CHOICES[@]}; do 
  if [ "$CHOICE" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  elif [ "$CHOICE" == "init-file-system" ]; then $CHOICE; break
  elif [ "$CHOICE" == "choose-species" ]; then $CHOICE; break
  elif [ "$CHOICE" == "cp-genome" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-index-genome" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bowtie-index-genome" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bowtie-align" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bowtie-refflat" ]; then $CHOICE; break
  elif [ "$CHOICE" == "segemehl-index-genome" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-align" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-align-per-fastq" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-samse" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-mpileup" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-samtools-view" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-samtools-sort" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-samtools-bed" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-batch-align-to-bed" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-degseq-bed" ]; then $CHOICE; break
  elif [ "$CHOICE" == "feature-genome" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-danko-countreadsininterval" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-R-saveimage" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-danko-writewiggle" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-danko-metagene" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-danko-detecttranscriptsem" ]; then $CHOICE; break
  elif [ "$CHOICE" == "fastxtoolkit-fastq_to_fasta" ]; then $CHOICE; break
  elif [ "$CHOICE" == "fastxtoolkit-fastx_quality_stats" ]; then $CHOICE; break
  elif [ "$CHOICE" == "fastq-summary" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-samtools-wig" ]; then $CHOICE; break
  elif [ "$CHOICE" == "fastq-sample" ]; then $CHOICE; break
  elif [ "$CHOICE" == "degseq" ]; then $CHOICE; break
  elif [ "$CHOICE" == "deseq" ]; then $CHOICE; break
  elif [ "$CHOICE" == "edgeR" ]; then $CHOICE; break
  elif [ "$CHOICE" == "warranty" ]; then $CHOICE; break
  elif [ "$CHOICE" == "copyright" ]; then $CHOICE; break
  elif [ "$CHOICE" == "transcript-genecoverage" ]; then $CHOICE; break
  elif [ "$CHOICE" == "transcript-parsernaseq" ]; then $CHOICE; break
  elif [ "$CHOICE" == "transcript-cufflinks" ]; then $CHOICE; break
  elif [ "$CHOICE" == "transcript-summary" ]; then $CHOICE; break
  elif [ "$CHOICE" == "maq-align" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-summary" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-pos2wig" ]; then $CHOICE; break
  elif [ "$CHOICE" == "de-count" ]; then $CHOICE; break
  elif [ "$CHOICE" == "job-truncate-reads" ]; then $CHOICE; break
  elif [ "$CHOICE" == "job-truncate-reads-local" ]; then $CHOICE; break
  elif [ "$CHOICE" == "batch" ]; then $CHOICE; break
  elif [ "$CHOICE" == "batch2" ]; then $CHOICE; break
  elif [ "$CHOICE" == "quit" ]; then break
  else
    echo -e "You need to enter something\n"
    continue
  fi
done

