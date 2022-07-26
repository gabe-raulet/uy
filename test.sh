#!/bin/bash
./bfs matrices/rgg_n_2_20_s0.mtx bfs_outputs/bfs.rgg.t1.txt 4 4 4 8 8 8 16 16 16 2>&1 | tee bfs_outputs/bfs.rgg.t1.log
