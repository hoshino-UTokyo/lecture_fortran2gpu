#!/bin/bash -x                                                                                                            
#PJM -g gt00
#PJM -L rscgrp=lecture-a
#PJM -L gpu=1
#PJM --mpi proc=1
#PJM -L elapse=00:10:00                                                                                                   
module purge
module load nvidia nvmpi

for exe in CUDA_U_f_a 
do
    echo ${exe}
    mpiexec -np ${PJM_MPI_PROC} ./wrapper.sh ./${exe}
done

