// Copyright 2021 Chun Shen
#ifndef FIELDS_H_
#define FIELDS_H_

#include <vector>
#include "cell.h"
#include "data_struct.h"

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
    std::vector<double> rhoq_;
    std::vector<double> rhos_;
    
    std::vector<double> piBulk_;

    // Vector & tensor fields
    std::vector<std::vector<double>> u_;
    std::vector<std::vector<double>> Wmunu_;

    Fields() = default;
    Fields(int Nx, int Ny, int Neta) {
        resizeFields(Nx, Ny, Neta);
    }

    ~Fields();

    int nX() const {return(Nx_);}
    int nY() const {return(Ny_);}
    int nEta() const {return(Neta_);}

    void resizeFields(int Nx, int Ny, int Neta);

    void clearFields();

    int getFieldIdx(int ix, int iy, int ieta) {
        return(ix + Nx_*(iy + Ny_*ieta));
    }

    int getFieldIdxHalo(int ix, int iy, int ieta);

    Cell_small getCell(const int idx);
    Cell_small getCell(int ix, int iy, int ieta);

    ReconstCell getCellIdeal(const int idx);
    ReconstCell getCellIdeal(int ix, int iy, int ieta);
};


template<class Func>
void FieldNeighbourLoopIdeal1(Fields &arena, int cx, int cy, int ceta,
                             Func func) {
    const std::array<int, 6> dx   = {-1, 1,  0, 0,  0, 0};
    const std::array<int, 6> dy   = { 0, 0, -1, 1,  0, 0};
    const std::array<int, 6> deta = { 0, 0,  0, 0, -1, 1};
    for(int dir = 0; dir < 3; dir++) {
        int m1nx   = dx  [2*dir];
        int m1ny   = dy  [2*dir];
        int m1neta = deta[2*dir];
        int p1nx   = dx  [2*dir+1];
        int p1ny   = dy  [2*dir+1];
        int p1neta = deta[2*dir+1];
        const auto  Ic = arena.getFieldIdxHalo(cx,      cy,      ceta       );
        const auto Ip1 = arena.getFieldIdxHalo(cx+p1nx, cy+p1ny, ceta+p1neta);
        const auto Im1 = arena.getFieldIdxHalo(cx+m1nx, cy+m1ny, ceta+m1neta);
        func(Ic, Ip1, Im1, dir+1);
    }
}


template<class Func>
void FieldNeighbourLoopIdeal2(Fields &arena, int cx, int cy, int ceta,
                             Func func) {
    const std::array<int, 6> dx   = {-1, 1,  0, 0,  0, 0};
    const std::array<int, 6> dy   = { 0, 0, -1, 1,  0, 0};
    const std::array<int, 6> deta = { 0, 0,  0, 0, -1, 1};
    for(int dir = 0; dir < 3; dir++) {
        int m1nx   = dx  [2*dir];
        int m1ny   = dy  [2*dir];
        int m1neta = deta[2*dir];
        int p1nx   = dx  [2*dir+1];
        int p1ny   = dy  [2*dir+1];
        int p1neta = deta[2*dir+1];
        int m2nx   = 2*m1nx;
        int m2ny   = 2*m1ny;
        int m2neta = 2*m1neta;
        int p2nx   = 2*p1nx;
        int p2ny   = 2*p1ny;
        int p2neta = 2*p1neta;
              auto  c = arena.getCellIdeal(cx,      cy,      ceta       );
        const auto p1 = arena.getCellIdeal(cx+p1nx, cy+p1ny, ceta+p1neta);
        const auto p2 = arena.getCellIdeal(cx+p2nx, cy+p2ny, ceta+p2neta);
        const auto m1 = arena.getCellIdeal(cx+m1nx, cy+m1ny, ceta+m1neta);
        const auto m2 = arena.getCellIdeal(cx+m2nx, cy+m2ny, ceta+m2neta);
        func(c, p1, p2, m1, m2, dir+1);
    }
}


template<class Func>
void FieldNeighbourLoop1(Fields &arena, int cx, int cy, int ceta, Func func) {
    const std::array<int, 6> dx   = {-1, 1,  0, 0,  0, 0};
    const std::array<int, 6> dy   = { 0, 0, -1, 1,  0, 0};
    const std::array<int, 6> deta = { 0, 0,  0, 0, -1, 1};
    for(int dir = 0; dir < 3; dir++) {
        int m1nx   = dx  [2*dir];
        int m1ny   = dy  [2*dir];
        int m1neta = deta[2*dir];
        int p1nx   = dx  [2*dir+1];
        int p1ny   = dy  [2*dir+1];
        int p1neta = deta[2*dir+1];
        const auto  Ic = arena.getFieldIdxHalo(cx,      cy,      ceta       );
        const auto Ip1 = arena.getFieldIdxHalo(cx+p1nx, cy+p1ny, ceta+p1neta);
        const auto Im1 = arena.getFieldIdxHalo(cx+m1nx, cy+m1ny, ceta+m1neta);
        func(Ic, Ip1, Im1, dir+1);
    }
}


template<class Func>
void FieldNeighbourLoop2(Fields &arena, int cx, int cy, int ceta, Func func) {
    const std::array<int, 6> dx   = {-1, 1,  0, 0,  0, 0};
    const std::array<int, 6> dy   = { 0, 0, -1, 1,  0, 0};
    const std::array<int, 6> deta = { 0, 0,  0, 0, -1, 1};
    for(int dir = 0; dir < 3; dir++) {
        int m1nx   = dx  [2*dir];
        int m1ny   = dy  [2*dir];
        int m1neta = deta[2*dir];
        int p1nx   = dx  [2*dir+1];
        int p1ny   = dy  [2*dir+1];
        int p1neta = deta[2*dir+1];
        int m2nx   = 2*m1nx;
        int m2ny   = 2*m1ny;
        int m2neta = 2*m1neta;
        int p2nx   = 2*p1nx;
        int p2ny   = 2*p1ny;
        int p2neta = 2*p1neta;
        const auto  Ic = arena.getFieldIdxHalo(cx,      cy,      ceta       );
        const auto Ip1 = arena.getFieldIdxHalo(cx+p1nx, cy+p1ny, ceta+p1neta);
        const auto Ip2 = arena.getFieldIdxHalo(cx+p2nx, cy+p2ny, ceta+p2neta);
        const auto Im1 = arena.getFieldIdxHalo(cx+m1nx, cy+m1ny, ceta+m1neta);
        const auto Im2 = arena.getFieldIdxHalo(cx+m2nx, cy+m2ny, ceta+m2neta);
        func(Ic, Ip1, Ip2, Im1, Im2, dir+1);
    }
}

#define FNLILAMBDAS1 [&](const int& Ic, const int& Ip1, const int& Im1, const int direction)

#define FNLILAMBDAS2 [&](ReconstCell& c, const ReconstCell& p1, const ReconstCell& p2, const ReconstCell& m1, const ReconstCell& m2, const int direction)

#define FNLLAMBDAS1 [&](const int& Ic, const int& Ip1, const int& Im1, const int direction)

#define FNLLAMBDAS2 [&](const int& Ic, const int& Ip1, const int& Ip2, const int& Im1, const int& Im2, const int direction)

#endif
