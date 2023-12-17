#!/bin/bash -x
#PJM -g gt00
#PJM -L rscgrp=lecture-a
#PJM -L gpu=1
#PJM --mpi proc=1
#PJM --omp thread=36
#PJM -L elapse=00:10:00

module purge
module load nvidia nvmpi                                                             

BIN=01_naive

mpiexec -np ${PJM_MPI_PROC} ./${BIN} 


