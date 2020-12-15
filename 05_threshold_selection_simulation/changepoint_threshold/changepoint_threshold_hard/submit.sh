#!/bin/bash
#PBS -N RChangepointThresholdHard
#PBS -l ncpus=4
#PBS -o Output/
#PBS -j oe
#PBS -t 1-10%5
#PBS -M z.varty@lancaster.ac.uk
#PBS -m abe

## Set num_threads to the same number as you set ncpus (4 by default)

num_threads=4

### Make sure the following is before you run Rscript (or other language)

export MKL_NUM_THREADS=$num_threads ###Limits Intel math kernel

export OPENBLAS_NUM_THREADS=$num_threads ### Limits OPENBLAS, this is the most important one

export MC_CORES=$num_threads ###Limits some packages in R

export OMP_NUM_THREADS=$num_threads ### Limits OpenMP

export NUMEXPR_NUM_THREADS=$num_threads ### Limits NUMEXR in python 

echo Moving to project directory
cd ~/changepoint_threshold_hard
echo Running R script in directory $PWD with command line argument: $PBS_ARRAYID
Rscript main.r $PBS_ARRAYID > Output/console/$PBS_ARRAYID.Rout

## Redirect the console output to the Output directory
