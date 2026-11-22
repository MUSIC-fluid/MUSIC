// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion
// collisions Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon,
// Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

#ifndef UTIL_H
#define UTIL_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include <fstream>
#include <memory>
#include <sstream>
#include <string>

#include "data_struct.h"

//! This is a utility class which contains a collection of helper functions.
namespace Util {
const double hbarc = 0.19733;
const double default_tol = 1.0e-8;
const double small_eps = 1e-16;

// the mass of a nucleon (averaged over proton and neutron)
const double m_N = 0.939;  // [GeV]

const int BT_BUF_SIZE = 500;

double theta(const double x);
double gmn(const int a);

double **mtx_malloc(const int n1, const int n2);

void mtx_free(double **m, const int n1, const int n2);

int IsFile(std::string);

int getParameter(
    std::string paramFileName, std::string paramName, int defaultValue);
double getParameter(
    std::string paramFileName, std::string paramName, double defaultValue);
std::string StringFind4(std::string file_name, std::string str_in);
std::string convert_to_lowercase(std::string str_in);

double lin_int(double x1, double x2, double f1, double f2, double x);

double four_dimension_linear_interpolation(
    double *lattice_spacing, double fraction[2][4], double ****cube);
double three_dimension_linear_interpolation(
    double *lattice_spacing, double fraction[2][3], double ***cube);
int binary_search(double *array, int length, double x);
void print_backtrace_errors();

// inline: called per cell and component in the hot loops
inline int map_2d_idx_to_1d(int a, int b) {
    static constexpr int index_map[5][4] = {
        {0, 1, 2, 3},
        {1, 4, 5, 6},
        {2, 5, 7, 8},
        {3, 6, 8, 9},
        {10, 11, 12, 13}};
    return index_map[a][b];
}

inline void map_1d_idx_to_2d(int idx_1d, int &a, int &b) {
    static constexpr int index_1d_a[5] = {1, 1, 1, 2, 2};
    static constexpr int index_1d_b[5] = {1, 2, 3, 2, 3};
    a = index_1d_a[idx_1d - 4];
    b = index_1d_b[idx_1d - 4];
}

Mat4x4 UnpackVecToMatrix(const Arr10 &in_vector);
Mat4x4 UnpackVecToMatrix(const ViscousVec &in_vector);
Mat4x4 UnpackVecToMatrix(const VorticityVec &in_vector);

// check whether a weak pointer is initialized or not
template <typename T>
bool weak_ptr_is_uninitialized(std::weak_ptr<T> const &weak) {
    using wt = std::weak_ptr<T>;
    return !weak.owner_before(wt {}) && !wt {}.owner_before(weak);
}

}  // namespace Util

#endif
