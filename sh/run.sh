#!/bin/bash

ROOTANALYSISDIR=`pwd`

#####################################################################
# External shell scripts
#####################################################################
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
source sh/segemehl-index-genome.sh


#####################################################################
# Read configuration file
#####################################################################
conf

#####################################################################
# Read directories
#####################################################################
SPECIESS=$(ls species|grep -v ^s)

#####################################################################
# Menus
#####################################################################
PS3="Select what you want to do with : "
CHOICES=( init-file-system \
          choose-species \
          cp-genome \
          bwa-index-genome \
          bwa-align \
          bwa-samse \
          bwa-samtools-view \
          bwa-samtools-sort \
          bwa-samtools-bed \
          segemehl-index-genome 
          )
select CHOICE in ${CHOICES[@]}; do 
  if [ "$CHOICE" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  elif [ "$CHOICE" == "init-file-system" ]; then $CHOICE; break
  elif [ "$CHOICE" == "choose-species" ]; then $CHOICE; break
  elif [ "$CHOICE" == "cp-genome" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-index-genome" ]; then $CHOICE; break
  elif [ "$CHOICE" == "segemehl-index-genome" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-align" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-samse" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-samtools-view" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-samtools-sort" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-samtools-bed" ]; then $CHOICE; break
  elif [ "$CHOICE" == "xxx" ]; then $CHOICE; break
  else
    echo -e "You need to enter something\n"
    continue
  fi
done

