#!/bin/bash
#PBS -l walltime=1:00:00,nodes=1
#PBS -A put_access_id
#PBS -j oe
#PBS -N sm-easyRNAseq
#PBS -q v4dev
#PBS -I

cd $PBS_O_WORKDIR

