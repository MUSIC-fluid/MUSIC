// Copyright 2021 Chun Shen
#ifndef FIELDS_H_
#define FIELDS_H_

#include <vector>
#include "cell.h"

class Fields {
 private:
    const int NWmunu_ = 14;  // dimension for pi^{\mu\nu}
    const int Nu_ = 4;       // dimension for u^\mu

    // field size in 3D
    int Nx_   = 0;
    int Ny_   = 0;
    int Neta_ = 0;

 public:
    // Scalar fields
    std::vector<double> e_;
    std::vector<double> rhob_;
    std::vector<double> piBulk_;

    // Vector & tensor fields
    std::vector<std::vector<double>> u_;
    std::vector<std::vector<double>> Wmunu_;

    Fields() = default;
    Fields(int Nx, int Ny, int Neta) {
        resizeFields(Nx, Ny, Neta);
    }

    ~Fields();

    void resizeFields(int Nx, int Ny, int Neta);

    int getFieldIdx(int ix, int iy, int ieta) {
        return(ix + Nx_*(iy + Ny_*ieta));
    }

    Cell_small getCell(int ix, int iy, int ieta);
};

#endif
