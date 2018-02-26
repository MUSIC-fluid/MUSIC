#!/usr/bin/env bash

export OMP_NUM_THREADS=16
export OMP_PROC_BIND=true
export OMP_PLACES=threads


echo "doing benchmark small ..."
time ./mpihydro benchmark/music_input_Gubser_small > benchmark/small.log
echo "doing benchmark middle ..."
time ./mpihydro benchmark/music_input_Gubser_middle > benchmark/middle.log
echo "doing benchmark large ..."
time ./mpihydro benchmark/music_input_Gubser_large > benchmark/large.log
