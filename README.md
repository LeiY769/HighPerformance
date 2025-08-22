INFO 0939

```plaintext
.
├── Part1/
│   ├── shallow.c           # Base sequential version
│   ├── shallow_opti.c      # Optimized sequential version
│   ├── shallow_mpi.c       # MPI parallel version
│   ├── shallow_openmp.c    # OpenMP parallel version
│   ├── shallow_combin.c    # Hybrid MPI+OpenMP version
│   └── lucia_gpu.job       # Job script for GPU execution
└── Part2/
    └── New_Source/
        ├── generate.c          # Create a new file
        Coriolis/
        │   shallow_opti.c     # Optimized version for Coriolis
        │   shallow_combin.c   # Hybrid MPI+OpenMP version for Coriolis
        │   shallow_gpu.c      # GPU version for Coriolis
        MPI_GPU/
            shallow_gpu_mpi.c  # GPU version for MPI


gcc shallow.c -o shallow -O3 -lm
gcc shallow_opti.c -o shallow_opti -O3 -lm
mpicc shallow_mpi.c -o shallow_mpi -O3 -lm
gcc shallow_openmp.c -o shallow_openmp -O3 -lm -fopenmp
mpicc shallow_combin.c -o shallow_combin -O3 -lm -fopenmp

gcc -o generate generate.c -lm

./generate "<input_file>"

All the files in the submission has an impact on what I have done during the project.

For the part 2, the principal investigation was focused on the Coriolis the 2 others was to explore for majority of myself.