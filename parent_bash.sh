#!/bin/bash
#SBATCH -p scavenge
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -t 280:00:00
#SBATCH -a 1-100
module load GROMACS/2021.5-foss-2020b-CUDA-11.1.1-PLUMED-2.7.3
f="$1"
runs="$2"
leftbound="$3"
rightbound="$4"
filepath="$5"
inputType="$6"
constantTop="$7"
python3 runner.py "$f" ${SLURM_ARRAY_TASK_ID} "$runs" "$leftbound" "$rightbound" "$filepath" "$inputType" "$constantTop"
