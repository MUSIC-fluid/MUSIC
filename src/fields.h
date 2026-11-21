// Copyright 2021 Chun Shen
#ifndef FIELDS_H_
#define FIELDS_H_

#include <array>
#include <vector>

#include "cell.h"
#include "data_struct.h"

class Fields {
  private:
    const int NWmunu_ = 14;  // dimension for pi^{\mu\nu}
    const int Nu_ = 4;       // dimension for u^\mu

    // field size in 3D
    int Nx_ = 0;
    int Ny_ = 0;
    int Neta_ = 0;

  public:
    // Scalar fields
    std::vector<double> e_;

    std::vector<double> rhob_;
    std::vector<double> rhoq_;
    std::vector<double> rhos_;

    std::vector<double> piBulk_;
    std::vector<double> piBulkChem_;

    // Vector & tensor fields
    std::vector<std::vector<double>> u_;
    std::vector<std::vector<double>> Wmunu_;

    Fields() = default;
    Fields(int Nx, int Ny, int Neta) { resizeFields(Nx, Ny, Neta); }

    ~Fields();

    int nX() const { return (Nx_); }
    int nY() const { return (Ny_); }
    int nEta() const { return (Neta_); }

    void clearFields();
    void resizeFields(int Nx, int Ny, int Neta);

    int getFieldIdx(int ix, int iy, int ieta) const {
        return (ix + Nx_ * (iy + Ny_ * ieta));
    }

    int getFieldIdxHalo(int ix, int iy, int ieta) const;

    Cell_small getCell(const int idx) const;
    Cell_small getCell(int ix, int iy, int ieta) const;

    ReconstCell getCellIdeal(const int idx) const;
    ReconstCell getCellIdeal(int ix, int iy, int ieta) const;
};

//! the KT fluxes of the 7 conserved currents through the cell faces normal
//! to x (direction 1), y (2) and eta (3). Face f along a direction lies
//! between the cells f - 1 and f, so there are n + 1 faces along it.
//! Stored per direction and component, faces contiguous.
class FaceFluxes {
  private:
    static constexpr int Ncomp_ = 7;
    int Nx_ = 0;
    int Ny_ = 0;
    int Neta_ = 0;
    std::array<int, 4> nFaces_ = {0, 0, 0, 0};
    std::array<std::vector<double>, 4> flux_;  // flux_[0] unused

  public:
    void resize(int Nx, int Ny, int Neta) {
        Nx_ = Nx;
        Ny_ = Ny;
        Neta_ = Neta;
        nFaces_ = {
            0, (Nx + 1) * Ny * Neta, Nx * (Ny + 1) * Neta,
            Nx * Ny * (Neta + 1)};
        for (int dir = 1; dir < 4; dir++) {
            flux_[dir].resize(static_cast<size_t>(Ncomp_) * nFaces_[dir]);
        }
    }

    //! index of the face (fx, fy, feta) normal to direction
    int getFaceIdx(int direction, int fx, int fy, int feta) const {
        const int nfx = Nx_ + (direction == 1);
        const int nfy = Ny_ + (direction == 2);
        return (fx + nfx * (fy + nfy * feta));
    }

    double& flux(int direction, int alpha, int faceIdx) {
        return flux_[direction]
                    [static_cast<size_t>(alpha) * nFaces_[direction] + faceIdx];
    }
    double flux(int direction, int alpha, int faceIdx) const {
        return flux_[direction]
                    [static_cast<size_t>(alpha) * nFaces_[direction] + faceIdx];
    }
};

template <class Func>
void FieldNeighbourLoopIdeal1(
    Fields& arena, int cx, int cy, int ceta, Func func) {
    const std::array<int, 6> dx = {-1, 1, 0, 0, 0, 0};
    const std::array<int, 6> dy = {0, 0, -1, 1, 0, 0};
    const std::array<int, 6> deta = {0, 0, 0, 0, -1, 1};
    for (int dir = 0; dir < 3; dir++) {
        int m1nx = dx[2 * dir];
        int m1ny = dy[2 * dir];
        int m1neta = deta[2 * dir];
        int p1nx = dx[2 * dir + 1];
        int p1ny = dy[2 * dir + 1];
        int p1neta = deta[2 * dir + 1];
        const auto Ic = arena.getFieldIdxHalo(cx, cy, ceta);
        const auto Ip1 =
            arena.getFieldIdxHalo(cx + p1nx, cy + p1ny, ceta + p1neta);
        const auto Im1 =
            arena.getFieldIdxHalo(cx + m1nx, cy + m1ny, ceta + m1neta);
        func(Ic, Ip1, Im1, dir + 1);
    }
}

template <class Func>
void FieldNeighbourLoop1(Fields& arena, int cx, int cy, int ceta, Func func) {
    const std::array<int, 6> dx = {-1, 1, 0, 0, 0, 0};
    const std::array<int, 6> dy = {0, 0, -1, 1, 0, 0};
    const std::array<int, 6> deta = {0, 0, 0, 0, -1, 1};
    for (int dir = 0; dir < 3; dir++) {
        int m1nx = dx[2 * dir];
        int m1ny = dy[2 * dir];
        int m1neta = deta[2 * dir];
        int p1nx = dx[2 * dir + 1];
        int p1ny = dy[2 * dir + 1];
        int p1neta = deta[2 * dir + 1];
        const auto Ic = arena.getFieldIdxHalo(cx, cy, ceta);
        const auto Ip1 =
            arena.getFieldIdxHalo(cx + p1nx, cy + p1ny, ceta + p1neta);
        const auto Im1 =
            arena.getFieldIdxHalo(cx + m1nx, cy + m1ny, ceta + m1neta);
        func(Ic, Ip1, Im1, dir + 1);
    }
}

template <class Func>
void FieldNeighbourLoop2(Fields& arena, int cx, int cy, int ceta, Func func) {
    const std::array<int, 6> dx = {-1, 1, 0, 0, 0, 0};
    const std::array<int, 6> dy = {0, 0, -1, 1, 0, 0};
    const std::array<int, 6> deta = {0, 0, 0, 0, -1, 1};
    for (int dir = 0; dir < 3; dir++) {
        int m1nx = dx[2 * dir];
        int m1ny = dy[2 * dir];
        int m1neta = deta[2 * dir];
        int p1nx = dx[2 * dir + 1];
        int p1ny = dy[2 * dir + 1];
        int p1neta = deta[2 * dir + 1];
        int m2nx = 2 * m1nx;
        int m2ny = 2 * m1ny;
        int m2neta = 2 * m1neta;
        int p2nx = 2 * p1nx;
        int p2ny = 2 * p1ny;
        int p2neta = 2 * p1neta;
        const auto Ic = arena.getFieldIdxHalo(cx, cy, ceta);
        const auto Ip1 =
            arena.getFieldIdxHalo(cx + p1nx, cy + p1ny, ceta + p1neta);
        const auto Ip2 =
            arena.getFieldIdxHalo(cx + p2nx, cy + p2ny, ceta + p2neta);
        const auto Im1 =
            arena.getFieldIdxHalo(cx + m1nx, cy + m1ny, ceta + m1neta);
        const auto Im2 =
            arena.getFieldIdxHalo(cx + m2nx, cy + m2ny, ceta + m2neta);
        func(Ic, Ip1, Ip2, Im1, Im2, dir + 1);
    }
}

#define FNLILAMBDAS1 \
    [&](const int& Ic, const int& Ip1, const int& Im1, const int direction)

#define FNLLAMBDAS1 \
    [&](const int& Ic, const int& Ip1, const int& Im1, const int direction)

#define FNLLAMBDAS2                                                    \
    [&](const int& Ic, const int& Ip1, const int& Ip2, const int& Im1, \
        const int& Im2, const int direction)

#endif
