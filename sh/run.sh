#!/bin/bash

ROOTANALYSISDIR=`pwd`

#####################################################################
# External shell scripts
#####################################################################
source sh/conf.sh
source sh/global-variable.sh
source sh/init-file-system.sh
source sh/bwa-index-genome.sh

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
          bwa-index-genome 
          )
select CHOICE in ${CHOICES[@]}; do 
  if [ "$CHOICE" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  elif [ "$CHOICE" == "init-file-system" ]; then $CHOICE; break
  elif [ "$CHOICE" == "bwa-index-genome" ]; then $CHOICE; break
  elif [ "$CHOICE" == "xxx" ]; then $CHOICE; break
  else
    echo -e "You need to enter something\n"
    continue
  fi
done

