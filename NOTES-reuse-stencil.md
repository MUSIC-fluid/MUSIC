Changes
=======

 align gmunu in data for faster access - also changed **gmunu in data to gmunu[4][4]
 moved two sums into w_rhs at top and bottom into one sum in the end
 do not symmetrize in the end, just go through all nu's
 vectorized innermost loop more efficiently by iterating over 4 indices instead of 3 to avoid masking

