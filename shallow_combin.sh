#!/bin/bash
#SBATCH --job-name="shampi_omp"
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=8     
#SBATCH --time=10:00      
#SBATCH --output=combi_%j.out  


export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}   
export OMP_PROC_BIND=spread                    
export OMP_PLACES=threads

PROJECT_DIR="${HOME}/project_info0939-1"
EXEC_DIR="${PROJECT_DIR}/bin"
INPUT_DIR="${PROJECT_DIR}"
PARAM_FILE="param_simple.txt"


cd ${SLURM_SUBMIT_DIR}

mpirun ${EXEC_DIR}/shallow_combin ${PARAM_FILE}