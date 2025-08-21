#!/bin/bash
#SBATCH --job-name="shallow_mpi"
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1          
#SBATCH --time=05:00
#SBATCH --output=mpi_%j.out

# define some variables
PROJECT_DIR="${HOME}/project_info0939-1"
EXEC_DIR="${PROJECT_DIR}/bin"
INPUT_DIR="${PROJECT_DIR}"
PARAM_FILE="param_simple.txt"

# move to the directory from which we submitted the job
cd ${SLURM_SUBMIT_DIR}

# run the MPI simulation
mpirun ${EXEC_DIR}/shallow_mpi ${PARAM_FILE}