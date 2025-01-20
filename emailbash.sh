#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -J REX1basin_1
#SBATCH -t 0:10:00
#SBATCH --mail-user=manuel.lopez@yale.edu
#SBATCH --mail-type=end
#SBATCH -p scavenge

echo "Sending email . . ."
sleep 1
echo "Email sent"