#!/bin/bash
#PJM -L rscgrp=lecture-a
#PJM -L gpu=4
#PJM --mpi proc=4
#PJM -L elapse=00:10:00
#PJM -g gt00

module purge
module load nvidia nvmpi

mpiexec -machinefile $PJM_O_NODEINF -n 4 ./acc_managed
mpiexec -machinefile $PJM_O_NODEINF -n 4 ./acc
mpiexec -machinefile $PJM_O_NODEINF -n 4 ./acc_async_overlap
