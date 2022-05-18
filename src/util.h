// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <fstream>
#include <string>
#include <memory>
#include <sys/stat.h>
#include "data_struct.h"

//! This is a utility class which contains a collection of helper functions.
namespace Util {
    const double hbarc = 0.19733;
    const double default_tol = 1.0e-8;
    const double small_eps = 1e-16;

    //the mass of a nucleon (averaged over proton and neutron)
    const double m_N = 0.939;   // [GeV]

    const int BT_BUF_SIZE = 500;

    double theta(const double x);
    double gmn(const int a);

    double **mtx_malloc(const int n1, const int n2);

    void mtx_free(double **m, const int n1, const int n2);

    int IsFile(std::string);

    std::string StringFind4(std::string file_name, std::string str_in);
    std::string convert_to_lowercase(std::string str_in);

    double lin_int(double x1,double x2,double f1,double f2,double x);

    double four_dimension_linear_interpolation(
            double* lattice_spacing, double fraction[2][4], double**** cube);
    double three_dimension_linear_interpolation(
            double* lattice_spacing, double fraction[2][3], double*** cube);
    int binary_search(double* array, int length, double x);
    void print_backtrace_errors();

    int map_2d_idx_to_1d(int a, int b);
    void map_1d_idx_to_2d(int idx_1d, int &a, int &b);

    Mat4x4 UnpackVecToMatrix(const Arr10 &in_vector);
    Mat4x4 UnpackVecToMatrix(const ViscousVec &in_vector);
    Mat4x4 UnpackVecToMatrix(const VorticityVec &in_vector);

    // check whether a weak pointer is initialized or not
    template <typename T>
    bool weak_ptr_is_uninitialized(std::weak_ptr<T> const& weak) {
        using wt = std::weak_ptr<T>;
        return !weak.owner_before(wt{}) && !wt{}.owner_before(weak);
    }

}

#endif
