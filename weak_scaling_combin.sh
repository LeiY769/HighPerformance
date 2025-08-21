#!/bin/bash
#
#SBATCH --exclusive
#SBATCH --partition=hmem
#SBATCH --output=job_output_weak_hybrid_%j.out
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --mem=0

module purge
module load OpenMPI

PROJECT_DIR="${HOME}/project_info0939-1"
EXEC_DIR="${PROJECT_DIR}/bin"
INPUT_DIR="${PROJECT_DIR}"

MIN_THREADS=2
MAX_THREADS=8
MIN_RANKS=${SLURM_NNODES}
EXE="${EXEC_DIR}/shallow_combin"

echo ""
echo " Weak scaling experiment start (hybrid)"
echo " ----------------------------------------"
echo "        Time: $(date)"
echo "      Job id: ${SLURM_JOBID}"
echo "       Nodes: $(scontrol show hostnames ${SLURM_JOB_NODELIST} | paste -s -d ' ')"
echo " Max threads: ${MAX_THREADS}"
echo " ---------------------------------------"
echo ""

for (( NTHREADS=${MIN_THREADS}; NTHREADS<=${MAX_THREADS}; NTHREADS*=2 )); do
  export OMP_NUM_THREADS=${NTHREADS}
  export OMP_PROC_BIND=close

  echo   ""
  echo   " +----------------------------------------------------+"
  printf " | Weak scaling with %2s threads per rank              |\n" ${NTHREADS}
  echo   " +------------+-----------------+---------------------+"
  echo   " | Num. ranks | Time (from app) | Time (time command) |"
  echo   " +------------+-----------------+---------------------+"

  MAX_RANKS=$((SLURM_CPUS_ON_NODE / NTHREADS * SLURM_NNODES))

  for (( NRANKS=${MIN_RANKS}; NRANKS<=${MAX_RANKS}; NRANKS*=2 )); do
    SRUN_OPTS="--nodes=${SLURM_NNODES}"
    SRUN_OPTS="${SRUN_OPTS} --ntasks-per-node=$((NRANKS / SLURM_NNODES))"
    SRUN_OPTS="${SRUN_OPTS} --cpus-per-task=${NTHREADS}"

    # total cores = ranks × threads
    TOTAL_CORES=$((NRANKS * NTHREADS))

    # choisir le fichier param_simple_TOTAL_CORES.txt
    PARAM_FILE="param_simple_${TOTAL_CORES}.txt"
    EXE_ARGS="${INPUT_DIR}/${PARAM_FILE}"

    OUTPUT_DIR="results_weak_hybrid_${SLURM_JOBID}/run-${NRANKS}ranks-${NTHREADS}threads"
    OUTPUT_FILE=$(basename ${EXE})-${NRANKS}ranks-${NTHREADS}threads.out

    mkdir -p ${OUTPUT_DIR}
    pushd ${OUTPUT_DIR} &> /dev/null

    cp ${INPUT_DIR}/h_simple.dat .

    /usr/bin/time -f '%e' -o 'time.out' srun ${SRUN_OPTS} ${EXE} ${EXE_ARGS} > ${OUTPUT_FILE}

    TIME_FROM_OUTPUT=$(sed -n -E 's/CPU Time: (.*) sec/\1/p' ${OUTPUT_FILE})
    TIME_FROM_TIME_CMD=$(cat time.out)

    if [[ -z "${TIME_FROM_OUTPUT}" ]]; then TIME_FROM_OUTPUT="NaN"; fi
    if [[ -z "${TIME_FROM_TIME_CMD}" ]]; then TIME_FROM_TIME_CMD="NaN"; fi

    printf " | %10s | %15s | %19s |\n" ${NRANKS} ${TIME_FROM_OUTPUT} ${TIME_FROM_TIME_CMD}

    popd &> /dev/null
  done

  echo " +------------+-----------------+---------------------+"
  echo ""
done

echo ""
echo " -----------------------------"
echo " End Time: $(date)"
echo " -----------------------------"
echo " Weak scaling experiment end"
