#!/bin/bash
#SBATCH --job-name="shallow_openmp"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00
#SBATCH --output=openmp_%j.out

# define some variables
PROJECT_DIR="${HOME}/project_info0939-1"
EXEC_DIR="${PROJECT_DIR}/bin"
INPUT_DIR="${PROJECT_DIR}"
PARAM_FILE="param_simple.txt"


# move to the directory from which we submitted the job
cd ${SLURM_SUBMIT_DIR}

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# run the simulation
${EXEC_DIR}/shallow_openmp ${PARAM_FILE}