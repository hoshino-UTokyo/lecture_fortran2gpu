#!/bin/bash -x                                                                                                            
#PJM -g gz00
#PJM -L rscgrp=regular-a
#PJM -L node=1
#PJM --mpi proc=1
#PJM -L elapse=01:00:00                                                                                                   

module purge
module load nvidia cuda ompi-cuda                                                                                   

BIN=01_naive

NP=${3:-${PJM_MPI_PROC}}
NPNODE=$(( ${NP} / ${PJM_VNODE} ))
mpiexec -np ${NP} -npernode ${NPNODE} ./wrap-aquarius.sh ./${BIN} 


