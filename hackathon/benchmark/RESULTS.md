Benchmarking Results
====================

Date        Time      Machine       Test      Wall-time
----        -----     -------       ----      ---------
2018-02-26  12:00     Cori KNL      Small     01:03.668
2018-02-26  12:00     Cori KNL      Medium    03:26.858
2018-02-26  12:00     Cori KNL      Large     06:31.306

new memory layout:
2018-02-27  17:00     Cori KNL      Small     01:36.577
2018-02-27  17:00     Cori KNL      Medium    03:57.133
2018-02-27  17:00     Cori KNL      Large     06:52.356

new minmod:
2018-02-27  17:40     Cori KNL      Small     02:08.459
2018-02-27  17:40     Cori KNL      Medium    04:22.419
2018-02-27  17:40     Cori KNL      Large     07:08.197

new compiler flags
2018-02-28  12:50     Cori KNL      Small     01:06.188
2018-02-28  12:50     Cori KNL      Medium    02:28.531
2018-02-28  12:50     Cori KNL      Large     04:14.172

-ipo flag in there
2018-02-28  18:00     Cori KNL      Small     01:09.359
2018-02-28  18:00     Cori KNL      Medium    02:23.585
2018-02-28  18:00     Cori KNL      Large     04:00.812

with #pragma in uWHRS
2018-02-28  18:00     Cori KNL      Small     01:17.038
2018-02-28  18:00     Cori KNL      Medium    02:41.463
2018-02-28  18:00     Cori KNL      Large     04:32.225

no more halo copying
2018-02-28  18:35     Cori KNL      Small     00:50.342
2018-02-28  18:35     Cori KNL      Medium    02:02.215
2018-02-28  18:35     Cori KNL      Large     03:40.763

16 threads no -g flag
2018-02-28  18:35     Cori KNL      Small     00:45.130
2018-02-28  18:35     Cori KNL      Medium    02:00.469
2018-02-28  18:35     Cori KNL      Large     03:39.614

32 threads
2018-02-28  22:35     Cori KNL      Small     00:35.816
2018-02-28  22:35     Cori KNL      Medium    01:12.237
2018-02-28  22:35     Cori KNL      Large     02:01.966

64 threads
2018-02-28  23:00     Cori KNL      Small     00:29.711
2018-02-28  23:00     Cori KNL      Medium    00:48.650
2018-02-28  23:00     Cori KNL      Large     01:18.369

16 threads SIMD uWRHS no -ipo
2018-03-01  12:00     Cori KNL      Small      00:46.454
2018-03-01  12:00     Cori KNL      Medium     02:01.209
2018-03-01  12:00     Cori KNL      Large      03:38.343

16 threads SIMD uWRHS with -ipo
2018-03-01  12:10     Cori KNL      Small      0m47.986s
2018-03-01  12:10     Cori KNL      Medium     1m58.253s
2018-03-01  12:10     Cori KNL      Large      3m34.379s

16 threads new code
2018-03-01  22:10     Cori KNL      Small      00:42.803
2018-03-01  22:10     Cori KNL      Medium     01:48.224
2018-03-01  22:10     Cori KNL      Large      03:15.796

16 threads new code - improved reconst
2018-03-02  10:40     Cori KNL      Small      00:44.908
2018-03-02  10:40     Cori KNL      Medium     01:45.722
2018-03-02  10:40     Cori KNL      Large      03:10.951

272 threads
2018-02-28  23:10     Cori KNL      Small      00:29.106
2018-02-28  23:10     Cori KNL      Medium     00:39.570
2018-02-28  23:10     Cori KNL      Large      00:51.547
2018-03-01  08:00     Cori KNL      L(32)      02:13.885
2018-03-01  09:16     Cori KNL      XL (128)   07:42.895
2018-03-02  11:00     Cori KNL      XL (128)   07:23.849

64x200x200                          Threads
----------                          -------
2018-03-02  13:00    Cori KNL       256 th.    03:59.553  
2018-03-02  13:00    Cori KNL       128 th.    04:40.912
2018-03-02  13:00    Cori KNL       64  th.    06:30.982    
2018-03-02  13:45    Cori Haswell   64  th.    03:19.622    

128x200x200
2018-03-02  11:00     Cori KNL      272 th.    07:23.849
2018-03-02  13:50     Cori Haswell  64  th.    06:29.131



Perf output:
============

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

2018-02-27 17:40 (Cori KNL)
----------------
Performance counter stats for './mpihydro benchmark/music_input_Gubser_small':

      195687337866      branches                                                      (30.01%)
        2189188186      branch-misses             #    1.12% of all branches          (29.99%)
        1410245681      cache-misses              #   29.923 % of all cache refs      (20.01%)
        4712898527      cache-references                                              (20.01%)
     1692102563645      cycles                                                        (30.01%)
             44778      faults                                                      
      959960802812      instructions              #    0.57  insns per cycle          (40.01%)
        1948435735      L1-dcache-load-misses     #    0.00% of all L1-dcache hits    (39.99%)
   <not supported>      L1-dcache-loads          
   <not supported>      L1-dcache-stores         
   <not supported>      LLC-load-misses          
        4407360971      LLC-loads                                                     (40.00%)
   <not supported>      LLC-store-misses         
        2235303274      LLC-stores                                                    (20.02%)
               291      migrations                                                  
     1693595999563      cycles                                                        (30.02%)

     102.181073494 seconds time elapsed

2018-02-28 12:50 (Cori KNL)
----------------
 Performance counter stats for './mpihydro benchmark/music_input_Gubser_small':

      124670443637      branches                                                      (29.99%)
         910413939      branch-misses             #    0.73% of all branches          (30.00%)
        1441473878      cache-misses              #   28.550 % of all cache refs      (20.01%)
        5049009695      cache-references                                              (20.01%)
     1079077235170      cycles                                                        (30.01%)
             44749      faults                                                      
      747109355418      instructions              #    0.69  insns per cycle          (40.02%)
        1862820999      L1-dcache-load-misses     #    0.00% of all L1-dcache hits    (40.01%)
   <not supported>      L1-dcache-loads          
   <not supported>      L1-dcache-stores         
   <not supported>      LLC-load-misses          
        4221388742      LLC-loads                                                     (40.01%)
   <not supported>      LLC-store-misses         
        2207814819      LLC-stores                                                    (20.00%)
               293      migrations                                                  
     1079899030978      cycles                                                        (30.00%)

      67.545363014 seconds time elapsed
