Changes
=======

 align gmunu in data for faster access - also changed **gmunu in data to gmunu[4][4]
 InitData DATA __attribute__ ((aligned (64))); in main.cpp (64 is for AVX512, 32 for AVX2...)
 moved two sums into w_rhs at top and bottom into one sum in the end
 do not symmetrize in the end, just go through all nu's
 vectorized innermost loop more efficiently by iterating over 4 indices instead of 3 to avoid masking

ipo problem
===========
 DATA-gmunu[3][nu] and DATA-gmunu[0][nu] in Make_uWRHS are aligned for alignment by (64) w/o ipo but not with ipo



