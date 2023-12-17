#!/bin/bash
export CUDA_VISIBLE_DEVICES=${OMPI_COMM_WORLD_LOCAL_RANK}
# Launch kernels with "omp nowait" to the same stream
export NV_OMP_AUTO_STREAMS=FALSE
$1
