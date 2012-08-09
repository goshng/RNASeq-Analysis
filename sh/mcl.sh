###############################################################################
# Copyright (C) 2012 Sang Chul Choi
#
# This file is part of RNAseq Analysis.
# 
# RNAseq Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RNAseq Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RNAseq Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

REPETITION=1
function mcl {
  select SPECIES in ${SPECIESS[@]}; do 
  if [ "$SPECIES" == "" ];  then
    echo -e "You need to enter something\n"
    continue
  else  
    batch3-repetition
    batch3-variable
    batch3-speciesfile 
    mcl-speciesfile 
    mcl-push-data
    mcl-get-data
    mcl-make-job
    mcl-run-mcl
    batch3-copy-scripts
    mcl-rmessage
    break
  fi
  done
}

# Set species file variables
function mcl-speciesfile {
  # Menu mcl
  MCLFASTA=$(grep ^MCLFASTA\: $SPECIESFILE | cut -d":" -f2)
  MCLFASTA=$ROOTANALYSISDIR/$MCLFASTA
  REFGENOMETXDB=$(grep ^REFGENOMETXDB\: $SPECIESFILE | cut -d":" -f2)  
  REFGENOMETXDB=$ROOTANALYSISDIR/$REFGENOMETXDB
  GENOMETXDB=$(basename $REFGENOMETXDB)
  BLASTEVALUE=$(grep ^BLASTEVALUE\: $SPECIESFILE | cut -d":" -f2)  
  MCLINFLATION=$(grep ^MCLINFLATION\: $SPECIESFILE | cut -d":" -f2)
}

# Send input files to the remote machine.
function mcl-push-data {
cat>$BASEDIR/mcl-push-data.sh<<EOF
#!/bin/bash
ssh -x $CAC_USERHOST mkdir -p $RDATADIR
ssh -x $CAC_USERHOST mkdir -p $RRUNANALYSIS
scp $MCLFASTA.gz $CAC_USERHOST:$RDATADIR
scp $REFGENOMETXDB $CAC_USERHOST:$RDATADIR
echo "Check cac:$RDATADIR"
EOF
}

function mcl-get-data {
cat>$BASEDIR/mcl-get-data.sh<<EOF
#!/bin/bash
scp $CAC_USERHOST:$RRUNANALYSIS/dump.out.all.blastp.tab.* $RUNANALYSIS
scp $CAC_USERHOST:$RRUNANALYSIS/out.all.blastp.tab.* $RUNANALYSIS
echo "Check $RUNANALYSIS"
EOF
}

function mcl-make-job {
cat>$BASEDIR/job-mcl<<EOF
for i in 12 14 16 18 20 22 24 26 28; do
  perl pl/find-core-gene.pl core \\
    -bed $RUNANALYSIS/feature.genes \\
    -mcl $RUNANALYSIS/dump.out.all.blastp.tab.I\$i \\
    -out $RUNANALYSIS/core.I\$i
done
echo "Check $RUNANALYSIS/core.I*"
EOF
}

function mcl-run-mcl {
cat>$BASEDIR/run-mcl.sh<<EOF
#!/bin/bash
nsub batch-mcl.sh
EOF

cat>$BASEDIR/batch-mcl.sh<<EOF
#!/bin/bash
#PBS -l walltime=${QCALIGNDEWALLTIME}:00:00,nodes=1
#PBS -A ${BATCHACCESS}
#PBS -j oe
#PBS -N $PROJECTNAME-MCL
#PBS -q ${QUEUENAME}
#PBS -m e
$EMAILON#PBS -M ${BATCHEMAIL}

function copy-data {
  cd \$TMPDIR

  # Programs and scripts.
  cp -r \$PBS_O_WORKDIR/pl .
  cp \$HOME/$MAKEBLASTDB .
  cp \$HOME/$BLASTP .
  cp \$HOME/$MCXDEBLAST .
  cp \$HOME/$MCXASSEMBLE .
  cp \$HOME/$MCLBLASTLINE .

  # All of the batchjob scripts.
  # No batchjob script.

  # Create output directories at the compute node.
  mkdir -p $CDATADIR
  mkdir -p $CRUNANALYSIS

  # No data files.
}

function process-data {
  cd \$TMPDIR

  MCLFASTA=\$(basename $MCLFASTA)

  # 1. Get a FASTA file of protein sequences.
  cp $RDATADIR/\$MCLFASTA.gz .
  gunzip \$MCLFASTA.gz

  # 2. Make BLASTDB of the FASTA file.
  ./makeblastdb -dbtype prot \\
    -input_type fasta \\
    -title mcl \\
    -in \$MCLFASTA

  # 3. BLAST protein sequences
  ./blastp -query \$MCLFASTA \\
    -db \$MCLFASTA \\
    -evalue $BLASTEVALUE \\
    -num_threads $NUMBERCPU \\
    -outfmt 7 \\
    -out all.blastp.tab

  perl mcxdeblast --m9 --score=e --sort=a all.blastp.tab    
  ./mcxassemble -b all.blastp.tab -r max --map
  for i in $MCLINFLATION; do
    perl mclblastline --start-mcl --mcl-I=\$i all.blastp.tab
  done    
  cp dump.out.all.blastp.tab.* $RRUNANALYSIS
  cp out.all.blastp.tab.* $RRUNANALYSIS

  # To see the output files in the compute nodes.
  tree
}
copy-data
process-data
cd
rm -rf \$TMPDIR
EOF
}

function mcl-rmessage {
  echo "bash $BASEDIR/mcl-push-data.sh"
  echo "work at cac:$CACWORKDIR"
  echo "bash $BASEDIR/mcl-get-data.sh"
  echo "bash $BASEDIR/job-mcl"
}
# END OF GOSEQ
################################################################################
