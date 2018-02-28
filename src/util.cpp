// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

#include "util.h"
#include <iostream>
#include <string>
#include <execinfo.h>

using namespace std;

namespace Util {

double ***cube_malloc(int n1, int n2, int n3)
{
    int i,j,k;
    double ***d1_ptr; 
    n1+=1;
    n2+=1;
    n3+=1;


    /* pointer to the n1*n2*n3 memory */
    d1_ptr = new double **[n1];

    for(i=0; i<n1; i++) 
     {
       d1_ptr[i] = new double *[n2];
     } 

    for(i=0; i<n1; i++)
    {
     for(j=0; j<n2; j++) 
      {
       d1_ptr[i][j] = new double [n3];
      }
    }

    for(i=0; i<n1; i++)
    {
     for(j=0; j<n2; j++) 
      {
	for(k=0; k<n3; k++) 
	  {
	    d1_ptr[i][j][k] = 0.0;
	  }
      }
    }

    return d1_ptr;
}/* cube_malloc */


void cube_free(double ***cube, int n1, int n2, int n3)
{
  int i,j;
  n1+=1;
  n2+=1;
  n3+=1;
  

  for(i=0; i<n1; i++)
    {
      for(j=0; j<n2; j++) 
	{
	  delete [] cube[i][j];
	}
    }
  for(j=0; j<n1; j++) 
    {
      delete [] cube[j];
    }
  
  delete [] cube;

}/* cube_free */


double **mtx_malloc(int n1, int n2)
{
    int i, j;
    double **d1_ptr; 
    d1_ptr = new double *[n1];
    
    for(i=0; i<n1; i++) 
     {
       d1_ptr[i] = new double [n2];
     }
    
    for(i=0; i<n1; i++) 
     {
      for(j=0; j<n2; j++) 
	{
	  d1_ptr[i][j] = 0.0;
	}
     }

return d1_ptr;
}


void mtx_free(double **m, int n1, int n2)
{
  int j;
  
 for(j=0; j<n1; j++) 
   {
     delete [] m[j];
   }

 delete [] m;

}


double *vector_malloc(int n1)
{
 double *d1_ptr;
 int i;

    /* pointer to the n1 array */
 d1_ptr = new double[n1];
 for(i=0; i<n1; i++) d1_ptr[i] = 0.0; 
    
 return d1_ptr;
}


void vector_free(double *vec)
{
  delete [] vec;
}

int IsFile(string file_name)
{
 FILE *temp;

 if( (temp = fopen(file_name.c_str(),"r")) == NULL) return 0;
 else 
  {
   fclose(temp);
   return 1;
  }
}/* IsFile */

// support comments in the parameters file
// comments need to start with #
string StringFind4(string file_name, string str_in) {
    string inputname = file_name;
    string str = str_in;

    string tmpfilename;
    tmpfilename = "input.default";
    
    // check whether the input parameter file is exist or not
    if (!IsFile(file_name)) {
        if (file_name == "") {
            fprintf(stderr, "No input file name specified.\n");
            fprintf(stderr, "Creating a default file named input.default\n");
        } else {
            cerr << "The file named " << file_name << " is absent." << endl;
            cout << "Creating " << file_name << "..." << endl;
            tmpfilename = file_name;
        }
        ofstream tmp_file(tmpfilename.c_str());
        tmp_file << "EndOfData" << endl;
        tmp_file.close();
        exit(1);
    }/* if isfile */
  
    // pass checking, now read in the parameter file
    string temp_string;
    ifstream input(inputname.c_str());
    getline(input, temp_string);  // read in the first entry

    int ind = 0;
    string para_name;
    string para_val;
    while (temp_string.compare("EndOfData") != 0) {
        // check whether it is the end of the file
        string para_string;
        stringstream temp_ss(temp_string);
        getline(temp_ss, para_string, '#');  // remove the comments
        if (para_string.compare("") != 0
                && para_string.find_first_not_of(' ') != std::string::npos) {
            // check the read in string is not empty
            stringstream para_stream(para_string);
            para_stream >> para_name >> para_val;
            if (para_name.compare(str) == 0) {
                // find the desired parameter
                ind++;
                input.close();
                return(para_val);
            }  /* if right, return */
        }
        getline(input, temp_string);  // read in the next entry
    }/* while */
    input.close(); // finish read in and close the file
    
    // the desired parameter is not in the parameter file, then return "empty"
    if (ind == 0) {
        return("empty");
    }
    // should not cross here !!!
    cout << "Error in StringFind4 !!!\n";
    return("empty");
}/* StringFind4 */


double lin_int(double x1,double x2,double f1,double f2,double x)
{
  double aa, bb;
  
  if (x2 == x1) 
    aa = 0.0;
  else
    aa =(f2-f1)/(x2-x1);
  bb = f1 - aa * x1;
  
  return aa*x + bb;
}

double four_dimension_linear_interpolation(
            double* lattice_spacing, double** fraction, double**** cube) {
    double denorm = 1.0;
    double results = 0.0;
    for (int i = 0; i < 4; i++) {
        denorm *= lattice_spacing[i];
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                for (int l = 0; l < 2; l++) {
                    results += (cube[i][j][k][l]*fraction[i][0]*fraction[j][1]
                                *fraction[k][2]*fraction[l][3]);
                }
            }
        }
    }
    results = results/denorm;
    return (results);
}

double three_dimension_linear_interpolation(
            double* lattice_spacing, double** fraction, double*** cube) {
    double denorm = 1.0;
    double results = 0.0;
    for (int i = 0; i < 3; i++) {
        denorm *= lattice_spacing[i];
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                results += (cube[i][j][k]*fraction[i][0]
                            *fraction[j][1]*fraction[k][2]);
            }
        }
    }
    results = results/denorm;
    return(results);
}

//! this function return the left index of the array where x sits in 
//! between array[idx] and array[idx+1]
//! this function assumes that the input array is monotonic 
int binary_search(double* array, int length, double x) {
    if (length < 3)  // array is too short
        return(0);

    int low_idx, mid_idx, high_idx;
    low_idx = 0;
    high_idx = length - 1;
    // first check the boundaries
    if ((array[low_idx] - x)*(array[high_idx] - x) > 0.) {
        fprintf(stderr, "Util::binary_search: can not find idx!\n");
        fprintf(stderr, "a[0] = %e, a[end] = %e, a = %e \n",
                array[low_idx], array[high_idx], x);
        exit(-1);
    }
    
    // find the index
    while (high_idx - low_idx > 1) {
        mid_idx = (int)((high_idx + low_idx)/2.);
        if ((array[low_idx] - x)*(array[mid_idx] - x) > 0.) {
            low_idx = mid_idx;
        } else {
            high_idx = mid_idx;
        }
    }
    return(low_idx);
}

void print_backtrace_errors() {
    int nptrs;
    void *buffer[BT_BUF_SIZE];
    char **strings;

    nptrs = backtrace(buffer, BT_BUF_SIZE);
    fprintf(stderr, "backtrace() return %d addresses\n", nptrs);
    strings = backtrace_symbols(buffer, nptrs);
    if (strings == NULL) {
        fprintf(stderr, "error in backtrace_symbols!");
        exit(1);
    }
    for (int j = 0; j < nptrs; j++) {
        fprintf(stderr, "%s\n", strings[j]);
    }
    free(strings);
    exit(1);
}

int map_2d_idx_to_1d(int a, int b) {
    // this function maps the 2d indices of a symmetric matrix to the index
    // in a 1-d array, which only stores the 10 independent components
    if (a == 4)
        return(10 + b);
    else if (a > b)  // symmetric matrix
        return(map_2d_idx_to_1d(b, a));
    if (b > 3) {
        cout << "Util::map_2d_idx_to_1d: index exceed dimension! "
             << "a = " << a << ", b = " << b << endl;
        exit(1);
    }
    if (a == 0)
        return(b);
    else if (a == 1)
        return(3 + b);
    else if (a == 2)
        return(5 + b);
    else if (a == 3)
        return(9);
    else {
        cout << "Util::map_2d_idx_to_1d: index exceed dimension! "
             << "a = " << a << ", b = " << b << endl;
        exit(1);
    }
}

}