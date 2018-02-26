Profiling
==============================

To get profiling access (on a local machine):

    sudo sh -c 'echo 1 >/proc/sys/kernel/perf_event_paranoid'

For timing:

    perf record ./mpihydro benchmark/music_input_Gubser_small

For cache misses and timing:

    perf record -e branches,cache-misses,cache-references,cycles,faults,instructions,L1-dcache-load-misses,L1-dcache-loads,L1-dcache-stores,l2_rqsts.miss,LLC-load-misses,LLC-loads,LLC-prefetch-misses,LLC-store-misses,LLC-stores,migrations,cycles:uppp ./mpihydro benchmark/music_input_Gubser_small

To see results:

    perf report