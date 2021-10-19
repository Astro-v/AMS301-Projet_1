#!/bin/bash

#SBATCH --job-name=MyJob
#SBATCH --output=output_%j.txt
#SBATCH --error=error_%j.txt

##SBATCH --partition=cpu_dist
#SBATCH --partition=cpu_test
##SBATCH --reservation=ams301_csp
#SBATCH --account=ams301

#SBATCH --ntasks=80
#SBATCH --time=00:10:00

## load modules
module load cmake/3.19.7
module load gcc/10.2.0
module load gmsh/4.8.4
module load openmpi/4.1.0

## execution
mpirun -display-map ${SLURM_SUBMIT_DIR}/a.out 1 1 1 0.5 6325 6325