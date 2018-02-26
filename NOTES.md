Profiling with perf
==============================

To get profiling access (on a local machine):

    sudo sh -c 'echo 1 >/proc/sys/kernel/perf_event_paranoid'

For timing:

    perf record ./mpihydro benchmark/music_input_Gubser_small

For cache misses and timing:

    perf record -e branches,cache-misses,cache-references,cycles,faults,instructions,L1-dcache-load-misses,L1-dcache-loads,L1-dcache-stores,l2_rqsts.miss,LLC-load-misses,LLC-loads,LLC-prefetch-misses,LLC-store-misses,LLC-stores,migrations,cycles:uppp ./mpihydro benchmark/music_input_Gubser_small

Quick stats:

    perf stat -e branches,cache-misses,cache-references,cycles,faults,instructions,L1-dcache-load-misses,L1-dcache-loads,L1-dcache-stores,l2_rqsts.miss,LLC-load-misses,LLC-loads,LLC-prefetch-misses,LLC-store-misses,LLC-stores,migrations,cycles ./mpihydro benchmark/music_input_Gubser_small

To see results:

    perf report



Running on Cori
==============================

Starting an interactive job:

    salloc -N1 -q debug --perf=vtune -C knl


Profiling with VTune
==============================

Node-level counter based performance profiler

Details [here](http://www.nersc.gov/users/software/performance-and-debugging-tools/vtune/)

Load VTune:

    module load vtune

Do a general exploration:

    amplxe-cl -finalization-mode=deferred -collect advanced-hotspots -r $SCRATCH/201802-knl-profiling ./mpihydro benchmark/music_input_Gubser_small

Things to collect:

 * advanced-hotspots
 * memory-access
 * general-exploration
 * hpc-performance

Other profiling tunes
=====================

aps #Application performing snapshot (included with VTune)
--------------------------------

    aps ./my_app


Intel advisor
-------------

For getting vectorization advice

    advixe-cl -c survey
    advixe-cl -c tripcounts -flop 