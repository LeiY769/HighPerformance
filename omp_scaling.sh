#!/bin/bash
#
#SBATCH --exclusive
#SBATCH --partition=hmem
#SBATCH --output=job_output_strong_openmp_%j.out
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --mem=0

module purge
module load GCC

PROJECT_DIR="${HOME}/project_info0939-1"
EXEC_DIR="${PROJECT_DIR}/bin"
INPUT_DIR="${PROJECT_DIR}"
PARAM_FILE="param_simple.txt"

MAX_THREADS="${SLURM_CPUS_ON_NODE}"
EXE="${EXEC_DIR}/shallow_openmp"
EXE_ARGS="${INPUT_DIR}/${PARAM_FILE}"

echo ""
echo " Strong scaling experiment start (OpenMP)"
echo " ----------------------------------------"
echo "        Time: $(date)"
echo "      Job id: ${SLURM_JOBID}"
echo "       Nodes: $(scontrol show hostnames ${SLURM_JOB_NODELIST} | paste -s -d ' ')"
echo " Max threads: ${MAX_THREADS}"
echo " ---------------------------------------"
echo ""

for OMP_BIND in close spread; do
  export OMP_PROC_BIND=$OMP_BIND

  echo   ""
  echo   " +------------------------------------------------------+"
  printf " | Strong scaling with OMP_PROC_BIND = %6s           |\n" $OMP_BIND
  echo   " +--------------+-----------------+---------------------+"
  echo   " | Num. threads | Time (from app) | Time (time command) |"
  echo   " +--------------+-----------------+---------------------+"

  for (( NTHREADS=1; NTHREADS<=${MAX_THREADS}; NTHREADS*=2 )); do
    export OMP_NUM_THREADS=${NTHREADS}

    OUTPUT_DIR="results_strong_openmp_${SLURM_JOBID}/run-${NTHREADS}threads-${OMP_PROC_BIND}"
    OUTPUT_FILE=$(basename ${EXE})-${NTHREADS}threads.out

    mkdir -p ${OUTPUT_DIR}
    pushd ${OUTPUT_DIR} &> /dev/null

    cp ${INPUT_DIR}/h_simple.dat .

    /usr/bin/time -f '%e' -o 'time.out' ${EXE} ${EXE_ARGS} > ${OUTPUT_FILE}

    TIME_FROM_OUTPUT=$(sed -n -E 's/CPU Time: (.*) sec/\1/p' ${OUTPUT_FILE})
    TIME_FROM_TIME_CMD=$(cat time.out)

    printf " | %12s | %15.3f | %19.3f |\n" ${NTHREADS} ${TIME_FROM_OUTPUT} ${TIME_FROM_TIME_CMD}

    popd &> /dev/null
  done

  echo " +--------------+-----------------+---------------------+"
  echo ""
done

echo ""
echo " -----------------------------"
echo " End Time: $(date)"
echo " -----------------------------"
echo " Strong scaling experiment end"
