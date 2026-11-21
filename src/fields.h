// Copyright 2021 Chun Shen
#ifndef FIELDS_H_
#define FIELDS_H_

#include <algorithm>
#include <array>
#include <vector>

#include "cell.h"
#include "data_struct.h"

//! one component of Fields: a non-owning view of its block in the
//! contiguous Fields buffer, indexed by Fields::getFieldIdx
class FieldComponent {
  private:
    double* p_ = nullptr;
    int n_ = 0;  // cells in the block, including the ghost cells
    friend class Fields;

  public:
    FieldComponent() = default;
    FieldComponent(const FieldComponent&) = delete;
    FieldComponent& operator=(const FieldComponent&) = delete;

    double& operator[](const int idx) { return p_[idx]; }
    const double& operator[](const int idx) const { return p_[idx]; }

    //! copies all cells (ghost cells included) of a component of the same size
    void copyFrom(const FieldComponent& other);

    //! sets all cells (ghost cells included) to value
    void fill(const double value) { std::fill(p_, p_ + n_, value); }
};

//! the hydrodynamic fields on the grid, stored as one contiguous, 64-byte
//! aligned buffer of Ncomp_ components (structure of arrays). Each block
//! holds the nx x ny x neta cells and Nghost_ ghost layers on both sides
//! along every direction with more than one cell, so that the stencils
//! can read up to Nghost_ cells past the grid edge without clamping. The
//! ghost cells hold copies of the nearest edge cell (updateGhosts), which
//! is the zero-gradient boundary condition that clamping used to give. A
//! direction with a single cell has no ghost layers and stride 0, so any
//! offset along it maps back to that cell.
class Fields {
  public:
    static constexpr int Nghost_ = 2;

  private:
    static constexpr int NWmunu_ = 14;  // dimension for pi^{\mu\nu}
    static constexpr int Nu_ = 4;       // dimension for u^\mu
    static constexpr int Ncomp_ = 6 + Nu_ + NWmunu_;
    static constexpr int alignDoubles_ = 8;  // 64 bytes

    // field size in 3D
    int Nx_ = 0;
    int Ny_ = 0;
    int Neta_ = 0;

    // ghost layers, padded sizes, and strides of the blocks
    int gx_ = 0, gy_ = 0, geta_ = 0;
    int Sx_ = 0, Sy_ = 0, Seta_ = 0;
    int strideX_ = 0, strideY_ = 0, strideEta_ = 0;
    int offset0_ = 0;      // index of the cell (0, 0, 0)
    int nStorage_ = 0;     // cells in a block, ghost cells included
    int blockStride_ = 0;  // nStorage_ rounded up to the alignment

    struct AlignedBuffer {
        double* ptr = nullptr;
        size_t size = 0;
        void allocate(size_t n);
        void release();
        ~AlignedBuffer() { release(); }
    };
    AlignedBuffer data_;

    void bindComponents();
    FieldComponent& component(int comp);
    const FieldComponent& component(int comp) const;

  public:
    // Scalar fields
    FieldComponent e_;

    FieldComponent rhob_;
    FieldComponent rhoq_;
    FieldComponent rhos_;

    FieldComponent piBulk_;
    FieldComponent piBulkChem_;

    // Vector & tensor fields
    std::array<FieldComponent, Nu_> u_;
    std::array<FieldComponent, NWmunu_> Wmunu_;

    Fields() = default;
    Fields(int Nx, int Ny, int Neta) { resizeFields(Nx, Ny, Neta); }
    Fields(const Fields& other);
    Fields& operator=(const Fields& other);

    ~Fields();

    int nX() const { return (Nx_); }
    int nY() const { return (Ny_); }
    int nEta() const { return (Neta_); }

    //! number of cells in a component block, ghost cells included: the
    //! range of every index from getFieldIdx and getFieldIdxHalo
    int nStorage() const { return (nStorage_); }

    void clearFields();
    void resizeFields(int Nx, int Ny, int Neta);

    //! copies all components of an arena of the same size
    void copyDataFrom(const Fields& other);

    //! sets the ghost cells to the nearest edge cell; needed after the
    //! cells were changed and before a stencil reads past the grid edge
    void updateGhosts();

    int getFieldIdx(int ix, int iy, int ieta) const {
        return (offset0_ + strideX_ * ix + strideY_ * iy + strideEta_ * ieta);
    }

    //! index of a cell up to Nghost_ cells past the grid edge (a ghost
    //! cell); it reads the right value only if the ghosts are up to date
    int getFieldIdxHalo(int ix, int iy, int ieta) const {
#ifdef MUSIC_CHECK_GHOSTS
        checkGhost(ix, iy, ieta);
#endif
        return (getFieldIdx(ix, iy, ieta));
    }

    //! aborts if the cell (ix, iy, ieta) is a ghost cell that differs from
    //! its edge cell or lies beyond the ghost layers
    void checkGhost(int ix, int iy, int ieta) const;

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
