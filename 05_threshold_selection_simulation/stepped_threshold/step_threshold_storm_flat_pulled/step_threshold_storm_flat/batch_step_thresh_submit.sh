#!/bin/bash
#PBS -N RBatchStepThresh
#PBS -l ncpus=1
#PBS -o Output/
#PBS -j oe
#PBS -t 1-100%8

echo Moving to project directory
cd ~/step_threshold_storm_flat
echo Running R script in directory $PWD with command line argument: $PBS_ARRAYID
Rscript batch_step_thresh.r $PBS_ARRAYID > Output/console/$PBS_ARRAYID.Rout

## Redirect the console output to the Output directory
