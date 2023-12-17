#!/bin/bash
#PJM -L rscgrp=lecture-a
#PJM -L gpu=1
#PJM -L elapse=00:10:00
#PJM -g gt00

module load nvidia 

export OMP_NUM_THREADS=36
./run
