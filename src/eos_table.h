// Copyright 2026 Chun Shen

#ifndef SRC_EOS_TABLE_H_
#define SRC_EOS_TABLE_H_

#include <cassert>
#include <cstddef>
#include <vector>

//! one tabulated quantity of an EoS made of several (n_B, e) tables:
//! table itable holds nb_length x e_length values, stored row by row
//! (e fastest), and all tables share one contiguous buffer. table[it][i][j]
//! is the value at the i-th n_B and j-th e of table it, as with the
//! double*** tables this replaces; data() and offset() give the flat
//! layout, e.g. for a copy to a device.
class EosTable {
  private:
    std::vector<double> data_;
    std::vector<size_t> offset_;
    std::vector<int> nbLength_;
    std::vector<int> eLength_;

  public:
    //! the rows of one table: rows[i] points to the i-th n_B row
    template <typename T>
    class Rows {
      private:
        T *base_;
        int eLength_;

      public:
        Rows(T *base, int eLength) : base_(base), eLength_(eLength) {}
        T *operator[](const int inb) const {
            return (base_ + static_cast<size_t>(inb) * eLength_);
        }
    };

    //! removes all tables and makes room for ntables of them
    void resize(const int ntables) {
        data_.clear();
        offset_.assign(ntables, 0);
        nbLength_.assign(ntables, 0);
        eLength_.assign(ntables, 0);
    }

    //! allocates table itable with nbLength x eLength values, set to 0
    void allocate(const int itable, const int nbLength, const int eLength) {
        assert(itable < static_cast<int>(offset_.size()));
        offset_[itable] = data_.size();
        nbLength_[itable] = nbLength;
        eLength_[itable] = eLength;
        data_.resize(
            data_.size() + static_cast<size_t>(nbLength) * eLength, 0.);
    }

    bool empty() const { return (offset_.empty()); }
    int nTables() const { return (static_cast<int>(offset_.size())); }

    Rows<double> operator[](const int itable) {
        return (Rows<double>(data_.data() + offset_[itable], eLength_[itable]));
    }
    Rows<const double> operator[](const int itable) const {
        return (Rows<const double>(
            data_.data() + offset_[itable], eLength_[itable]));
    }

    const double *data() const { return (data_.data()); }
    size_t size() const { return (data_.size()); }
    size_t offset(const int itable) const { return (offset_[itable]); }
};

#endif  // SRC_EOS_TABLE_H_
