#!/bin/bash -x                                                                                                            
#PJM -g gt00
#PJM -L rscgrp=lecture-a
#PJM -L gpu=1
#PJM --mpi proc=1
#PJM -L elapse=00:10:00                                                                                                   
module purge
module load nvidia nvmpi

export NVIDIA_ACC_TIME=1

for exe in ACC_U_f_a 
do
    echo ${exe}
    mpiexec -np ${PJM_MPI_PROC} ./wrapper.sh ./${exe}
done

