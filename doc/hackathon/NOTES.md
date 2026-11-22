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

Generate a report ([details](https://software.intel.com/en-us/vtune-amplifier-help-report)):

    amplxe-cl -report summary -r $SCRATCH/201802-knl-profiling

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


Regular Expressions
===================

Change the phrase `arena[asdf][bob][ted]` into `arena(bob,ted,asdf)` in all files:

    ls *cpp *h | xargs -n 1 sed -i -E 's/arena\[([^]]*)\]\[([^]]*)\]\[([^]]*)\]/arena(\2,\3,\1)/g'


Compiler Flags
===================

 * `-fsanitize=address`       : Check for memory leaks/errors (don't use for production runs)
 * `-O3`                      : Many optimizations
 * `-g`                       : Debugging info (use this always)
 * `-xmic-avx512`             : Use for KNL     (icc)
 * `-xavx2`                   : Use for Haswell (icc)
 * `-march=native`            : Use full capabilities of current CPU (gcc)
 * `-debug inline-debug-info` : For running with VTune (icc) (Can use this always, only makes EXE file larger)
 * `-ipo`                     : Interprocedural optimization (icc) (Same as LTO). Longer compilation time, optimizing between files.
 * `-flto`                    : Link-Time Optimization (gcc) (Same as IPO). Longer compilation time, optimizing between files.
 * `-fma`                     : Fused multiply-add: a faster, more precise way of doing `a=a+b*c`
 * `-align`                   : Ensure memory is allocated to cache-line boundaries, makes accesses faster
 * `-fimf-use-svml`           : Power implementation same as vector mass

Suggested flags INtel compiler:

    -g -std=c11 -O3 -xMIC-AVX512 -fma -align -finline-functions

Suggested flags GCC compiler:

    -std=c11 -march=knl -O3 -mavx512f -mavx512pf -mavx512er -mavx512cd -mfma -malign-data=cacheline -finline-functions



What's vectorizing? What isn't? And why?
========================================

Compiler outputs
----------------

    -qopt-report=5 -qopt-report-phase:vec



Intel advisor
-------------

    /hpcgpfs01/software/Intel/psxe2018.u1/parallel_studio_xe_2018/advisor_2018/bin64/advixe-cl


    /hpcgpfs01/software/Intel/psxe2018.u1/parallel_studio_xe_2018/advisor_2018/bin64/advixe-cl --collect survey -project-dir=profile/  ./a.out 3001 120 out.dem 123

    /hpcgpfs01/software/Intel/psxe2018.u1/parallel_studio_xe_2018/advisor_2018/bin64/advixe-cl --collect tripcounts -flop -project-dir=profile/  ./a.out 3001 120 out.dem 123

    /hpcgpfs01/software/Intel/psxe2018.u1/parallel_studio_xe_2018/advisor_2018/bin64/advixe-cl --collect map -mark-up-list=26 -project-dir=profile/ ./a.out 501 120 out.dem 123



OpenMP Scheduling
===================

 * `schedule(static)`  : Allocate iterations to threads and never changed
 * `schedule(dynamic)` : Threads take new work after each iteration (probably slow for small workloads)
 * `schedule(guided)`  : Threads take new chunks of work after each iteration. Chunksize automagically adjusts.


BNL KNL Cluster
===================

    salloc -A hackathon --reservation=hackathon -p long -t hh:mm:ss

Advisor on Cori via NX and ssh -X
=================================

 * first run command line advisor with your code in an interactive session - this will generate the needed data in the myproj folder:
 
    salloc -N 1 -C knl,quad,cache -t 30:00 -q debug  
    module load advisor
    export OMP_NUM_THREADS=8
    module load impi
    mpiicc -g -openmp -O3 -o mycode.exe mycode.c  #make with debug flags
    export I_MPI_PMI_LIBRARY=/opt/slurm/default/lib/pmi/libpmi.so
    srun -n 1 -c 8 advixe-cl --collect survey --trace-mpi --project-dir ./myproj  -- ./mpihydro input
    srun -n 1 -c 8 advixe-cl --collect tripcounts -flop --trace-mpi --project-dir ./myproj  -- ./mpihydro input

 * open NX session either via browser or app (see http://www.nersc.gov/users/connecting-to-nersc/using-nx/)
 * open konsole or some such and ssh -X to cori
 * Alt-Ctrl-0 opens the NX options where you can resize and maximize the window (on a mac at least)
 * on cori 

    module load advisor

  then 

    advixe-gui ./myproj
