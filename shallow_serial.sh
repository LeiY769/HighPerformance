#!/bin/bash
#SBATCH --job-name="shallow serial"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00
#SBATCH --output=shallow_serial_%j.out

# define some variables
PROJECT_DIR="${HOME}/project_info0939-1"
EXEC_DIR="${PROJECT_DIR}/bin"
INPUT_DIR="${PROJECT_DIR}"
PARAM_FILE="param_simple.txt"

# copy all files from the input directory to
# the directory from which we submitted the job
cp ${INPUT_DIR}/* ${SLURM_SUBMIT_DIR}

# move to the directory from which we submitted the job
cd ${SLURM_SUBMIT_DIR}

# run the simulation
${EXEC_DIR}/shallow ${PARAM_FILE}