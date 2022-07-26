#!/bin/bash

export OMP_NUM_THREADS=$1
ID=$2

./bfs matrices/rgg_n_2_20_s0.mtx bfs_outputs/bfs.rgg.t${OMP_NUM_THREADS}.${ID}.txt 4 4 4 8 8 8 16 16 16 2>&1 | tee bfs_outputs/bfs.rgg.t${OMP_NUM_THREADS}.${ID}.log
