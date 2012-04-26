#!/bin/bash
#PBS -l walltime=8:00:00,nodes=1
#PBS -A acs4_0001
#PBS -j oe
#PBS -N sm-cov
#PBS -q v4
#PBS -m e
# #PBS -M schoi@cornell.edu
#PBS -t 1-70

# number of nodes * 8
# 13 * 8 = 104 > 100
# 26 * 8 = 208 > 200
# 39 * 8 = 312
# 50 * 8 = 400
# 60 * 8 = 480
# 70 * 8 = 560
NUMSPLIT=560

function copy-data {
  cd $TMPDIR

  # Programs and scripts.
  cp $PBS_O_WORKDIR/batch-coverage.R .
}

function process-data {
  cd $TMPDIR
  CORESPERNODE=8
  for (( i=1; i<=CORESPERNODE; i++)); do
    g=$((CORESPERNODE * (PBS_ARRAYID-1) + i))
    if [ $g -le $NUMSPLIT ]; then
      /home/fs01/sc2265/Downloads/r-devel/b/bin/Rscript batch-coverage.R $NUMSPLIT $g &
    fi
  done
}

copy-data
process-data
wait
cd
rm -rf $TMPDIR
