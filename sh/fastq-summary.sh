###############################################################################
# Copyright (C) 2011,2012 Sang Chul Choi
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

function fastq-summary {
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

      GENOMEFASTA=$(basename $REFGENOMEFASTA)
      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      FASTQFILES=$(grep ^FASTQFILES\: $SPECIESFILE | cut -d":" -f2)

      echo -n "Do you wish to copy fastq files to the data directory? (e.g., y/n) "
      read WISH
      if [ "$WISH" == "y" ]; then
        for g in $FASTQFILES; do
          FASTQNUM=FASTQ$(printf "%03d" $g)
          GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
          cp $GZIPFASTAQFILE $DATADIR/$FASTQNUM.fq.gz
          echo See $DATADIR/$FASTQNUM.fq.gz
        done
      fi

      echo -n "Do you wish to make reports for the named fastq files? (e.g., y/n) "
      read WISH
      if [ "$WISH" == "y" ]; then
        QUALITYSCORELIST=""
        for g in $FASTQFILES; do
          QUALITYSCORENUM=QUALITYSCORE$(printf "%03d" $g)
          QUALITYSCORE=$(grep $QUALITYSCORENUM $SPECIESFILE | cut -d":" -f2)
          QUALITYSCORELIST="$QUALITYSCORELIST \"$QUALITYSCORE\""
        done
        QUALITYSCORELISTINR=$(echo $QUALITYSCORELIST | sed -e 's/[ ]/,/g')
        echo $QUALITYSCORELISTINR
        fastq-summary-using-qrqc
      fi

      break
    fi
  done
}

# Try this.
# As an example, using the data that come with the qrqc package:
# 
# First create the plot you want, on the node.
# > s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))
# > toplot <- qualPlot(s.fastq)
# > save(list = "toplot", file = "qualPlot.Rdata")
# 
# You could do a whole bunch of different plots, and then save them all using
# > save(list = ("firstplot","secondplot","thirdplot"), file = "myplots.Rdata")
# 
# Now move the qualPlot.Rdata file to the head node (or just read it in from the compute node).
# 
# > load("qualPlot.Rdata")
# > library(ggplot2)
# > plot(toplot)
# > ggsave("./qualplot.png")

function fastq-summary-using-qrqc {
FASTQFILESINR=$(echo $FASTQFILES | sed -e 's/[ ]/,/g')
cat>$RUNANALYSIS/$FUNCNAME.R<<EOF
library(qrqc)
fq.quality <- c($QUALITYSCORELISTINR)
fq.number <- c($FASTQFILESINR)
datadir <- "$DATADIR"
bwadir <- "$BWADIR"
for (i in 1:length(fq.number)) {
  fq.name <- sprintf("%s/FASTQ%03d.fq.gz", datadir, fq.number[i])
  fq.file <- readSeqFile(fq.name,quality=fq.quality[i])
  # makeReport(fq.file, outputDir=bwadir)
  toplot <- qualPlot(fq.file)
  fq.plot <- sprintf("%s/FASTQ%03d.qualPlot.RData", bwadir, fq.number[i])
  save(list="toplot", file = fq.plot)
}
EOF
cat>$RUNANALYSIS/$FUNCNAME-plot.R<<EOF
library(qrqc)
library(ggplot2)
fq.quality <- c($QUALITYSCORELISTINR)
fq.number <- c($FASTQFILESINR)
datadir <- "$DATADIR"
bwadir <- "$BWADIR"
for (i in 1:length(fq.number)) {
  fq.plot <- sprintf("%s/FASTQ%03d.qualPlot.RData", bwadir, fq.number[i])
  fq.pdf <- sprintf("%s/FASTQ%03d.qualPlot.pdf", bwadir, fq.number[i])
  load(fq.plot)
  pdf(fq.pdf)
  plot(toplot)
  dev.off()
}
EOF
  echo Use R 2.15 and Bioconductor 2.10!
  echo Rscript-2.15 $RUNANALYSIS/$FUNCNAME.R
}
