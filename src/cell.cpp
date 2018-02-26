// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include "./util.h"
#include "./cell.h"

using namespace std;

Cell::Cell() {
   prev_epsilon = 0.0;
   epsilon_t    = 0.0;
   epsilon      = 0.0;

   prev_rhob = 0.0;
   rhob_t    = 0.0;
   rhob      = 0.0;
}

Cell ***Cell::grid_c_malloc(int n1, int n2, int n3) {
   int i,j;
   //int k,inc;
   Cell ***d1_ptr;
   //Cell *tmp_ptr;

   d1_ptr = new Cell **[n1];

   for (i=0; i<n1; i++) {
      d1_ptr[i] = new Cell *[n2];
   } 
   
   for (i=0; i<n1; i++) {
      for (j=0; j<n2; j++) {
         d1_ptr[i][j] = new Cell[n3];
      }
   }
   return d1_ptr;
}/* grid_c_malloc */
