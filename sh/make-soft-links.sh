#!/bin/bash

echo Create symbolic links at the local base directory.
echo Choose your local and storage base directoires.
echo Delete the line of exit in the script, and execute this shell script again.
LOCALBASEDIR=/home/sc2265/Documents/Projects/RNASeq-Analysis
STORAGEBASEDIR=/home/sc2265/data/RNASeq-Analysis
# exit 

BASEDIR=`pwd`
echo Creating soft links to directories in the source code base directory...
for i in pl sh src run; do
ln -s $BASEDIR/$i $LOCALBASEDIR/$i 
done

echo Creating soft links to directories in the storage base directory...
for i in data downloads email log output; do
ln -s $STORAGEBASEDIR/$i $LOCALBASEDIR/$i
done

echo Creating conf directory in the local base directory...
echo Copying a local conf file to the create conf directory...
echo Edit the conf file available at $LOCALBASEDIR/conf/README
mkdir $LOCALBASEDIR/conf
cp $BASEDIR/conf/README.local $LOCALBASEDIR/conf/README
