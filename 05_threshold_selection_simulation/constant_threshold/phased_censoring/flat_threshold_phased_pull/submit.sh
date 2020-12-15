#!/bin/bash
#PBS -N RBatchFlatThresh
#PBS -l ncpus=1
#PBS -o Output/
#PBS -j oe
#PBS -t 1-500%10

echo Moving to project directory
cd ~/flat_threshold_phased
echo Running R script in directory $PWD with command line argument: $PBS_ARRAYID
Rscript main.r $PBS_ARRAYID > Output/console/$PBS_ARRAYID.Rout

## Redirect the console output to the Output directory
