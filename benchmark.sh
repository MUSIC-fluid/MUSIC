#!/usr/bin/env bash

export OMP_NUM_THREADS=16
export OMP_PROC_BIND=true
export OMP_PLACES=threads


echo "doing benchmark small ..."
time ./mpihydro benchmark/music_input_Gubser_small > benchmark/small.log
#perf stat -e branches,branch-misses,cache-misses,cache-references,cycles,faults,instructions,L1-dcache-load-misses,L1-dcache-loads,L1-dcache-stores,LLC-load-misses,LLC-loads,LLC-store-misses,LLC-stores,migrations,cycles ./mpihydro benchmark/music_input_Gubser_small
echo "doing benchmark middle ..."
time ./mpihydro benchmark/music_input_Gubser_middle > benchmark/middle.log
echo "doing benchmark large ..."
time ./mpihydro benchmark/music_input_Gubser_large > benchmark/large.log
