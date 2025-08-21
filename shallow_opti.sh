#!/bin/bash
#SBATCH --job-name="shallow opti"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=15:00
#SBATCH --output=opti_%j.out

# define some variables
PROJECT_DIR="${HOME}/project_info0939-1"
EXEC_DIR="${PROJECT_DIR}/bin"
INPUT_DIR="${PROJECT_DIR}"
PARAM_FILE="param_simple.txt"



# move to the directory from which we submitted the job
cd ${SLURM_SUBMIT_DIR}

# run the simulation
${EXEC_DIR}/shallow_opti ${PARAM_FILE}