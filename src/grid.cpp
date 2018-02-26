// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include "./util.h"
#include "./grid.h"

using namespace std;

Grid::Grid() {
   prev_epsilon = 0.0;
   epsilon_t    = 0.0;
   epsilon      = 0.0;

   prev_rhob = 0.0;
   rhob_t    = 0.0;
   rhob      = 0.0;
}

Grid ***Grid::grid_c_malloc(int n1, int n2, int n3) {
   int i,j;
   //int k,inc;
   Grid ***d1_ptr;
   //Grid *tmp_ptr;

   d1_ptr = new Grid **[n1];

   for (i=0; i<n1; i++) {
      d1_ptr[i] = new Grid *[n2];
   } 
   
   for (i=0; i<n1; i++) {
      for (j=0; j<n2; j++) {
         d1_ptr[i][j] = new Grid[n3];
      }
   }
   return d1_ptr;
}/* grid_c_malloc */
