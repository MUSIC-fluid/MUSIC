// Copyright 2021 Chun Shen

#include "fields.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>

#include "cell.h"

void FieldComponent::copyFrom(const FieldComponent& other) {
    assert(n_ == other.n_);
    std::copy(other.p_, other.p_ + n_, p_);
}

void Fields::AlignedBuffer::allocate(size_t n) {
    release();
    void* mem = nullptr;
    if (n > 0
        && posix_memalign(
               &mem, alignDoubles_ * sizeof(double), n * sizeof(double))
               != 0) {
        throw std::bad_alloc();
    }
    ptr = static_cast<double*>(mem);
    size = n;
}

void Fields::AlignedBuffer::release() {
    free(ptr);
    ptr = nullptr;
    size = 0;
}

Fields::~Fields() { clearFields(); }

Fields::Fields(const Fields& other) {
    resizeFields(other.Nx_, other.Ny_, other.Neta_);
    copyDataFrom(other);
}

Fields& Fields::operator=(const Fields& other) {
    if (this != &other) {
        if (Nx_ != other.Nx_ || Ny_ != other.Ny_ || Neta_ != other.Neta_) {
            resizeFields(other.Nx_, other.Ny_, other.Neta_);
        }
        copyDataFrom(other);
    }
    return (*this);
}

FieldComponent& Fields::component(int comp) {
    if (comp < 6) {
        FieldComponent* scalars[6] = {&e_,    &rhob_,   &rhoq_,
                                      &rhos_, &piBulk_, &piBulkChem_};
        return (*scalars[comp]);
    } else if (comp < 6 + Nu_) {
        return (u_[comp - 6]);
    }
    return (Wmunu_[comp - 6 - Nu_]);
}

const FieldComponent& Fields::component(int comp) const {
    return (const_cast<Fields*>(this)->component(comp));
}

void Fields::bindComponents() {
    for (int comp = 0; comp < Ncomp_; comp++) {
        FieldComponent& c = component(comp);
        c.p_ =
            (data_.ptr == nullptr
                 ? nullptr
                 : data_.ptr + static_cast<size_t>(comp) * blockStride_);
        c.n_ = nStorage_;
    }
}

void Fields::clearFields() {
    data_.release();
    Nx_ = Ny_ = Neta_ = 0;
    gx_ = gy_ = geta_ = 0;
    Sx_ = Sy_ = Seta_ = 0;
    strideX_ = strideY_ = strideEta_ = 0;
    offset0_ = nStorage_ = blockStride_ = 0;
    bindComponents();
}

void Fields::resizeFields(int Nx, int Ny, int Neta) {
    clearFields();
    Nx_ = Nx;
    Ny_ = Ny;
    Neta_ = Neta;
    gx_ = (Nx > 1 ? Nghost_ : 0);
    gy_ = (Ny > 1 ? Nghost_ : 0);
    geta_ = (Neta > 1 ? Nghost_ : 0);
    Sx_ = Nx + 2 * gx_;
    Sy_ = Ny + 2 * gy_;
    Seta_ = Neta + 2 * geta_;
    strideX_ = (Nx > 1 ? 1 : 0);
    strideY_ = (Ny > 1 ? Sx_ : 0);
    strideEta_ = (Neta > 1 ? Sx_ * Sy_ : 0);
    offset0_ = gx_ + Sx_ * (gy_ + Sy_ * geta_);
    nStorage_ = Sx_ * Sy_ * Seta_;
    blockStride_ =
        (nStorage_ + alignDoubles_ - 1) / alignDoubles_ * alignDoubles_;
    data_.allocate(static_cast<size_t>(Ncomp_) * blockStride_);
    std::fill(data_.ptr, data_.ptr + data_.size, 0.);
    bindComponents();
    std::fill(u_[0].p_, u_[0].p_ + nStorage_, 1.);
}

void Fields::copyDataFrom(const Fields& other) {
    assert(Nx_ == other.Nx_ && Ny_ == other.Ny_ && Neta_ == other.Neta_);
    std::copy(other.data_.ptr, other.data_.ptr + other.data_.size, data_.ptr);
}

void Fields::updateGhosts() {
    if (nStorage_ == 0) return;
    // the x ghosts of the rows inside the grid ...
    if (gx_ > 0) {
#pragma omp parallel for collapse(2)
        for (int comp = 0; comp < Ncomp_; comp++) {
            for (int peta = geta_; peta < geta_ + Neta_; peta++) {
                double* p = component(comp).p_;
                for (int py = gy_; py < gy_ + Ny_; py++) {
                    double* row = p + Sx_ * (py + Sy_ * peta);
                    for (int k = 0; k < gx_; k++) {
                        row[k] = row[gx_];
                        row[gx_ + Nx_ + k] = row[gx_ + Nx_ - 1];
                    }
                }
            }
        }
    }
    // ... then the ghost rows: copies of the nearest row inside the grid,
    // with its x ghosts
#pragma omp parallel for collapse(2)
    for (int comp = 0; comp < Ncomp_; comp++) {
        for (int peta = 0; peta < Seta_; peta++) {
            double* p = component(comp).p_;
            const int ceta =
                geta_ + std::max(0, std::min(Neta_ - 1, peta - geta_));
            for (int py = 0; py < Sy_; py++) {
                const int cy = gy_ + std::max(0, std::min(Ny_ - 1, py - gy_));
                if (ceta == peta && cy == py) continue;
                const double* src = p + Sx_ * (cy + Sy_ * ceta);
                std::copy(src, src + Sx_, p + Sx_ * (py + Sy_ * peta));
            }
        }
    }
}

void Fields::checkGhost(int ix, int iy, int ieta) const {
    const int n[3] = {Nx_, Ny_, Neta_};
    const int i[3] = {ix, iy, ieta};
    int c[3];
    bool ghost = false;
    for (int d = 0; d < 3; d++) {
        if (i[d] < -Nghost_ || i[d] >= n[d] + Nghost_) {
            std::cerr << "Fields::checkGhost: cell (" << ix << ", " << iy
                      << ", " << ieta << ") is beyond the ghost layers"
                      << std::endl;
            std::abort();
        }
        c[d] = std::max(0, std::min(n[d] - 1, i[d]));
        if (c[d] != i[d] && n[d] > 1) ghost = true;
    }
    if (!ghost) return;
    const int idx = getFieldIdx(ix, iy, ieta);
    const int cidx = getFieldIdx(c[0], c[1], c[2]);
    for (int comp = 0; comp < Ncomp_; comp++) {
        const double g = component(comp)[idx];
        const double v = component(comp)[cidx];
        if (std::memcmp(&g, &v, sizeof(double)) != 0) {
            std::cerr << "Fields::checkGhost: stale ghost cell (" << ix << ", "
                      << iy << ", " << ieta << "), component " << comp << ": "
                      << g << " vs " << v << std::endl;
            std::abort();
        }
    }
}

Cell_small Fields::getCell(int ix, int iy, int ieta) const {
    ix = std::max(0, std::min(Nx_ - 1, ix));
    iy = std::max(0, std::min(Ny_ - 1, iy));
    ieta = std::max(0, std::min(Neta_ - 1, ieta));
    int fieldIdx = getFieldIdx(ix, iy, ieta);
    return (getCell(fieldIdx));
}

Cell_small Fields::getCell(const int idx) const {
    Cell_small cell;
    cell.epsilon = e_[idx];

    cell.rhob = rhob_[idx];
    cell.rhoq = rhoq_[idx];
    cell.rhos = rhos_[idx];

    for (int i = 0; i < Nu_; i++) {
        cell.u[i] = u_[i][idx];
    }
    for (int i = 0; i < NWmunu_; i++) {
        cell.Wmunu[i] = Wmunu_[i][idx];
    }
    cell.pi_b = piBulk_[idx];
    cell.pi_b_chem = piBulkChem_[idx];
    return (cell);
}

Cell_small Fields::getCellWithoutW(const int idx) const {
    Cell_small cell;
    cell.epsilon = e_[idx];
    cell.rhob = rhob_[idx];
    cell.rhoq = rhoq_[idx];
    cell.rhos = rhos_[idx];
    for (int i = 0; i < Nu_; i++) {
        cell.u[i] = u_[i][idx];
    }
    // kept as is when turn_on_bulk = 1 and turn_on_bulk_chem = 0
    cell.pi_b_chem = piBulkChem_[idx];
    return (cell);
}

ReconstCell Fields::getCellIdeal(const int idx) const {
    ReconstCell cell;
    cell.e = e_[idx];

    cell.rhob = rhob_[idx];
    cell.rhoq = rhoq_[idx];
    cell.rhos = rhos_[idx];

    for (int i = 0; i < Nu_; i++) {
        cell.u[i] = u_[i][idx];
    }
    return (cell);
}

ReconstCell Fields::getCellIdeal(int ix, int iy, int ieta) const {
    ix = std::max(0, std::min(Nx_ - 1, ix));
    iy = std::max(0, std::min(Ny_ - 1, iy));
    ieta = std::max(0, std::min(Neta_ - 1, ieta));
    int fieldIdx = getFieldIdx(ix, iy, ieta);
    return (getCellIdeal(fieldIdx));
}
