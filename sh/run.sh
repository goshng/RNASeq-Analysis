#!/bin/bash

#####################################################################
# 
#####################################################################

source sh/bwa-help.sh

#####################################################################
# Menus
#####################################################################
PS3="Select what you want to do with : "
CHOICES=( bwa-help
          )
select CHOICE in ${CHOICES[@]}; do 
  if [ "$CHOICE" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  elif [ "$CHOICE" == "bwa-help" ]; then $CHOICE; break
  else
    echo -e "You need to enter something\n"
    continue
  fi
done

