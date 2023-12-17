#!/bin/bash -x                                                                                                            
#PJM -g gt00
#PJM -L rscgrp=lecture-a
#PJM -L gpu=1
#PJM --mpi proc=1
#PJM -L elapse=00:10:00                                                                                                   
module purge
module load nvidia nvmpi                                                                                   
for exe in Omp_U_f_a Omp_U_f_c Omp_U_s_a Omp_U_s_c
do
    echo ${exe}
    mpiexec -np ${PJM_MPI_PROC} ./wrapper.sh ./${exe}
done

