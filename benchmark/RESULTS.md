Benchmarking Results
====================

Date        Time      Machine       Test      Wall-time
----        -----     -------       ----      ---------
2018-02-26  12:00     Cori KNL      Small     01:03.668
2018-02-26  12:00     Cori KNL      Medium    03:26.858
2018-02-26  12:00     Cori KNL      Large     06:31.306

2018-02-27  17:00     Cori KNL      Small     01:36.577
2018-02-27  17:00     Cori KNL      Medium    03:57.133
2018-02-27  17:00     Cori KNL      Large     06:52.356






Perf output:

2018-02-26 12:00 (Cori KNL)
----------------
Performance counter stats for './mpihydro benchmark/music_input_Gubser_small':

      104768124880      branches                                                      (30.00%)
        2458174963      branch-misses             #    2.35% of all branches          (30.00%)
         477794393      cache-misses              #   16.528 % of all cache refs      (20.01%)
        2890880173      cache-references                                              (20.01%)
     1045456695954      cycles                                                        (30.00%)
             20542      faults                                                      
      618850615563      instructions              #    0.59  insns per cycle          (40.00%)
        1193338490      L1-dcache-load-misses     #    0.00% of all L1-dcache hits    (39.99%)
   <not supported>      L1-dcache-loads          
   <not supported>      L1-dcache-stores         
   <not supported>      LLC-load-misses          
        3715962692      LLC-loads                                                     (40.00%)
   <not supported>      LLC-store-misses         
         408176152      LLC-stores                                                    (20.00%)
               277      migrations                                                  
     1045751180536      cycles                                                        (30.00%)

     703.467977811 seconds time elapsed


2018-02-27 17:00 (Cori KNL)
----------------
Performance counter stats for './mpihydro benchmark/music_input_Gubser_small':

      199409159376      branches                                                      (30.01%)
        2398172640      branch-misses             #    1.20% of all branches          (30.00%)
        1390208457      cache-misses              #   28.773 % of all cache refs      (20.00%)
        4831597633      cache-references                                              (20.01%)
     1686078406363      cycles                                                        (30.00%)
             44794      faults                                                      
      972876912335      instructions              #    0.58  insns per cycle          (40.00%)
        1928954862      L1-dcache-load-misses     #    0.00% of all L1-dcache hits    (40.00%)
   <not supported>      L1-dcache-loads          
   <not supported>      L1-dcache-stores         
   <not supported>      LLC-load-misses          
        4593474395      LLC-loads                                                     (40.00%)
   <not supported>      LLC-store-misses         
        2193255311      LLC-stores                                                    (20.01%)
               292      migrations                                                  
     1687274462725      cycles                                                        (30.01%)

      98.379163794 seconds time elapsed


