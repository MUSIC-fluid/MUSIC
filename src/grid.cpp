#include "util.h"
#include "grid.h"

using namespace std;

Grid::Grid()
{
   prev_epsilon = 0.0;
   epsilon_t = 0.0;
   epsilon = 0.0;

   prev_rhob = 0.0;
   rhob_t = 0.0;
   rhob = 0.0;

   p = 0.0;
   p_t = 0.0;

}

Grid::~Grid()
{

}


Grid *Grid::grid_v_malloc(int n1)
{
   Grid *d1_ptr;
   //int i;
  
   /* pointer to the n1 array */
   d1_ptr = new Grid[n1];
   //(Grid *) malloc (sizeof(Grid )*n1);
   
   return d1_ptr;
}/* grid_v_malloc */


Grid **Grid::grid_m_malloc(int n1, int n2)
{
   int i;
   //int j;
   Grid **d1_ptr, *tmp_ptr;

   tmp_ptr = (Grid *)malloc(sizeof(Grid)*n1*n2);
   d1_ptr = (Grid **) malloc (sizeof(Grid *)*n1);

   for(i=0; i<n1; i++) 
   {
      d1_ptr[i] = &(tmp_ptr[i*n2]);
   }
    
   return d1_ptr;
}/* grid_m_malloc */

Grid ***Grid::grid_c_malloc(int n1, int n2, int n3)
{
   int i,j;
   //int k,inc;
   Grid ***d1_ptr;
   //Grid *tmp_ptr;

   d1_ptr = new Grid **[n1];

   for(i=0; i<n1; i++) 
   {
      d1_ptr[i] = new Grid *[n2];
   } 
   
   for(i=0; i<n1; i++)
   {
      for(j=0; j<n2; j++) 
      {
         d1_ptr[i][j] = new Grid[n3];
      }
   }
   return d1_ptr;
}/* grid_c_malloc */
