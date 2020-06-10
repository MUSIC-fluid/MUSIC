// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

#include "util.h"
#include <iostream>
#include <string>
#include <execinfo.h>
#include <algorithm>

using std::string;

namespace Util {


double theta(const double x) {
    if (x < 0.) {
        return(0.0);
    } else {
        return(1.0);
    }
}


double gmn(const int a) {
    if (a == 0) {
        return(-1.0);
    } else {
        return(1.0);
    }
}


double **mtx_malloc(const int n1, const int n2) {
    double **d1_ptr; 
    d1_ptr = new double *[n1];

    for(int i=0; i<n1; i++) 
        d1_ptr[i] = new double [n2];

    for(int i=0; i<n1; i++) 
    for(int j=0; j<n2; j++) 
        d1_ptr[i][j] = 0.0;

    return d1_ptr;
}



void mtx_free(double **m, const int n1, const int n2) {
    for (int j = 0; j < n1; j++) 
        delete [] m[j];
    delete [] m;
}



int IsFile(string file_name) {
    FILE *temp;

    if ((temp = fopen(file_name.c_str(),"r")) == NULL) {
        return(0);
    } else  {
        fclose(temp);
        return(1);
    }
}


// support comments in the parameters file
// comments need to start with #
// case-insensitive
string StringFind4(string file_name, string str_in) {
    string inputname = file_name;
    string str = convert_to_lowercase(str_in);

    string tmpfilename;
    tmpfilename = "input.default";

    // check whether the input parameter file is exist or not
    if (!IsFile(file_name)) {
        if (file_name == "") {
            fprintf(stderr, "No input file name specified.\n");
            fprintf(stderr, "Creating a default file named input.default\n");
        } else {
            std::cerr << "The file named " << file_name << " is absent."
                      << std::endl;
            std::cout << "Creating " << file_name << "..." << std::endl;
            tmpfilename = file_name;
        }
        std::ofstream tmp_file(tmpfilename.c_str());
        tmp_file << "EndOfData" << std::endl;
        tmp_file.close();
        exit(1);
    }/* if isfile */

    // pass checking, now read in the parameter file
    string temp_string;
    std::ifstream input(inputname.c_str());
    getline(input, temp_string);  // read in the first entry

    int ind = 0;
    string para_name;
    string para_val;
    while (convert_to_lowercase(temp_string).compare("endofdata") != 0) {
        // check whether it is the end of the file
        string para_string;
        std::stringstream temp_ss(temp_string);
        getline(temp_ss, para_string, '#');  // remove the comments
        if (para_string.compare("") != 0
                && para_string.find_first_not_of(' ') != std::string::npos) {
            // check the read in string is not empty
            std::stringstream para_stream(para_string);
            para_stream >> para_name >> para_val;
            if (convert_to_lowercase(para_name).compare(str) == 0) {
                // find the desired parameter
                ind++;
                input.close();
                return(para_val);
            }  /* if right, return */
        }
        getline(input, temp_string);  // read in the next entry
    }
    input.close(); // finish read in and close the file

    // the desired parameter is not in the parameter file, then return "empty"
    if (ind == 0) {
        return("empty");
    }
    // should not cross here !!!
    std::cout << "Error in StringFind4 !!!\n";
    return("empty");
}


// this function convert string to lower case
string convert_to_lowercase(string str_in) {
    std::transform(str_in.begin(), str_in.end(), str_in.begin(), ::tolower);
    return(str_in);
}


double lin_int(double x1,double x2,double f1,double f2,double x) {
    double aa, bb;

    if (x2 == x1) {
        aa = 0.0;
    } else {
        aa =(f2-f1)/(x2-x1);
    }
    bb = f1 - aa * x1;

    return aa*x + bb;
}

double four_dimension_linear_interpolation(
            double* lattice_spacing, double fraction[2][4], double**** cube) {
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
            double* lattice_spacing, double fraction[2][3], double*** cube) {
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
    static const int index_map[5][4] = {{ 0,  1,  2,  3},
                                        { 1,  4,  5,  6},
                                        { 2,  5,  7,  8},
                                        { 3,  6,  8,  9},
                                        {10, 11, 12, 13}};
    return index_map[a][b];
}


void map_1d_idx_to_2d(int idx_1d, int &a, int &b) {
    static const int index_1d_a[5] = {1, 1, 1, 2, 2};
    static const int index_1d_b[5] = {1, 2, 3, 2, 3};
    a = index_1d_a[idx_1d - 4];
    b = index_1d_b[idx_1d - 4];
}


Mat4x4 UnpackVecToMatrix(const Arr10 &in_vector) {
    Mat4x4 out_matrix;
    out_matrix[0][0] = in_vector[0];
    out_matrix[0][1] = in_vector[1];
    out_matrix[0][2] = in_vector[2];
    out_matrix[0][3] = in_vector[3];
    out_matrix[1][0] = in_vector[1];
    out_matrix[1][1] = in_vector[4];
    out_matrix[1][2] = in_vector[5];
    out_matrix[1][3] = in_vector[6];
    out_matrix[2][0] = in_vector[2];
    out_matrix[2][1] = in_vector[5];
    out_matrix[2][2] = in_vector[7];
    out_matrix[2][3] = in_vector[8];
    out_matrix[3][0] = in_vector[3];
    out_matrix[3][1] = in_vector[6];
    out_matrix[3][2] = in_vector[8];
    out_matrix[3][3] = in_vector[9];
    return out_matrix;
}


Mat4x4 UnpackVecToMatrix(const ViscousVec &in_vector) {
    Mat4x4 out_matrix;
    out_matrix[0][0] = in_vector[0];
    out_matrix[0][1] = in_vector[1];
    out_matrix[0][2] = in_vector[2];
    out_matrix[0][3] = in_vector[3];
    out_matrix[1][0] = in_vector[1];
    out_matrix[1][1] = in_vector[4];
    out_matrix[1][2] = in_vector[5];
    out_matrix[1][3] = in_vector[6];
    out_matrix[2][0] = in_vector[2];
    out_matrix[2][1] = in_vector[5];
    out_matrix[2][2] = in_vector[7];
    out_matrix[2][3] = in_vector[8];
    out_matrix[3][0] = in_vector[3];
    out_matrix[3][1] = in_vector[6];
    out_matrix[3][2] = in_vector[8];
    out_matrix[3][3] = in_vector[9];
    return out_matrix;
}


Mat4x4 UnpackVecToMatrix(const VorticityVec &in_vector) {
    Mat4x4 out_matrix;
    out_matrix[0][0] = 0.0;
    out_matrix[0][1] = in_vector[0];
    out_matrix[0][2] = in_vector[1];
    out_matrix[0][3] = in_vector[2];
    out_matrix[1][0] = -out_matrix[0][1];
    out_matrix[1][1] = 0.0;
    out_matrix[1][2] = in_vector[3];
    out_matrix[1][3] = in_vector[4];
    out_matrix[2][0] = -out_matrix[0][2];
    out_matrix[2][1] = -out_matrix[1][2];
    out_matrix[2][2] = 0.0;
    out_matrix[2][3] = in_vector[5];
    out_matrix[3][0] = -out_matrix[0][3];
    out_matrix[3][1] = -out_matrix[1][3];
    out_matrix[3][2] = -out_matrix[2][3];
    out_matrix[3][3] = 0.0;
    return out_matrix;
}

}
