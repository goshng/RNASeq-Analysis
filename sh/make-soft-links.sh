#!/bin/bash

echo Create symbolic links at the local base directory.
echo Choose your local and storage base directoires.
echo Delete the line of exit in the script, and execute this shell script again.
LOCALBASEDIR=$HOME/Documents/Projects/RNASeq-Analysis
STORAGEBASEDIR=$HOME/data/RNASeq-Analysis
# exit 

mkdir -p $LOCALBASEDIR

BASEDIR=`pwd`
echo Creating soft links to directories in the source code base directory...
for i in doc pl sh src run; do
rm $LOCALBASEDIR/$i 
ln -s $BASEDIR/$i $LOCALBASEDIR/$i 
done

echo Creating soft links to directories in the storage base directory...
for i in conf data downloads email log output species; do
rm $LOCALBASEDIR/$i 
ln -s $STORAGEBASEDIR/$i $LOCALBASEDIR/$i
done

echo Creating conf directory in the local base directory...
echo Copying a local conf file to the create conf directory...
echo Edit the conf file available at $LOCALBASEDIR/conf/README
# mkdir $LOCALBASEDIR/routput
# mkdir $LOCALBASEDIR/conf
# cp $BASEDIR/conf/README.local $LOCALBASEDIR/conf/README

