#!/bin/bash
#PBS -N RChangepointThresholdHard
#PBS -l ncpus=1
#PBS -o changepoint_threshold_hard/Output/bashthings
#PBS -j oe
#PBS -t 1-10%5

echo Moving to project directory
cd ~/changepoint_threshold_hard
echo Running R script in directory $PWD with command line argument: $PBS_ARRAYID
Rscript main.r $PBS_ARRAYID > Output/console/$PBS_ARRAYID.Rout

## Redirect the console output to the Output directory
