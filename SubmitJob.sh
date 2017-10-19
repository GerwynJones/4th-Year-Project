#!/bin/bash

#PBS -l select=1:ncpus=1:mem=12gb
#PBS -N myjob
#PBS -l walltime=1:00:00
#PBS -q serial
#PBS -P PR329

# Run your command here
echo "Hello from job number $PBS_JOBID"

# Create an output directory
OUT_DIR=/scratch/$USER/$PBS_JOBID
mkdir -p $OUT_DIR
cd $OUT_DIR

echo "The runtime directory for this job is $OUT_DIR"

# Copy the job submission script into the output directory
cp ${PBS_O_WORKDIR}/jobsub3.sh ${OUT_DIR}/