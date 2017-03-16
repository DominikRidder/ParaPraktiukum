#!/bin/bash -x
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=8
#SBATCH --output=cannon_out.%j
#SBATCH --error=cannon_err.%j
#SBATCH --time=00:00:30
#SBATCH --partition=batch

srun cannon
