// Copyright 2026 Chun Shen
//
// DeviceHydro runs the ideal hydro update of Evolve::AdvanceRK on a device
// (input parameter gpu_backend):
//   0  off (default): the CPU code of Advance
//   1  the device kernels on the CPU (OpenMP), in double they reproduce
//      Advance bit for bit: the reference of the GPU backends
//   2  Metal (Apple GPUs, float)
//   3  CUDA (NVIDIA GPUs, double by default)
// gpu_float = 1 runs backends 1 and 3 in float (Metal always is);
// gpu_store_q (default: 1 in float, 0 in double) carries tau T^{tau mu} on
// the device instead of rebuilding it from the primitive variables in every
// stage, which in float accumulates rounding errors (doc/gpu_porting_*).
//
// The state stays on the device between steps. After each step the new
// slice is copied into the Fields that AdvanceRK would have filled, so
// freeze-out, output and the diagnostics run on the host as before. The
// host Fields are uploaded again only when the device copy is stale: the
// first step, a new slice pointer or grid, and after invalidate(), which
// Evolve calls when it changes the fields itself (beastMode regrid,
// change of delta_tau).
//
// Supported: ideal hydro (viscosity_flag 0) without source terms, EoS 0,
// 9, 91 and 10 - 14. Otherwise the CPU code runs and a warning says why.

#ifndef SRC_DEVICE_DEVICE_HYDRO_H_
#define SRC_DEVICE_DEVICE_HYDRO_H_

#include <memory>
#include <string>
#include <vector>

#include "../data.h"
#include "../pretty_ostream.h"

class EOS;
class Fields;

namespace music_device {

//! one device: see device_api.h, here as a C++ interface
class Backend {
  public:
    virtual ~Backend() {}
    virtual std::string name() const = 0;
    virtual int realSize() const = 0;
    virtual bool resize(int nx, int ny, int neta, std::string &err) = 0;
    virtual bool setEos(
        const void *desc, size_t descBytes, const void *table,
        size_t tableBytes, std::string &err) = 0;
    virtual void *state() = 0;
    virtual bool upload(std::string &err) = 0;
    virtual bool step(
        const void *params, size_t paramBytes, int &flags,
        std::string &err) = 0;
    virtual bool download(std::string &err) = 0;
};

}  // namespace music_device

class DeviceHydro {
  public:
    //! hasHydroSource: Advance adds source terms (not on the device)
    DeviceHydro(
        const InitData &DATA, const EOS &eos, double absErr, double relErr,
        int maxIter, double vCritical, bool hasHydroSource);
    ~DeviceHydro();

    //! the device advances the fluid
    bool active() const { return (backend_ != nullptr); }

    //! one Heun step from *fpCurr at tau; rotates the pointers as
    //! Evolve::AdvanceRK does
    void advanceStep(
        double tau, Fields *&fpPrev, Fields *&fpCurr, Fields *&fpNext);

    //! the host changed the fields: upload before the next step
    void invalidate() { valid_ = false; }

  private:
    const InitData &DATA;
    const EOS &eos_;
    mutable pretty_ostream music_message;
    std::unique_ptr<music_device::Backend> backend_;
    bool float_ = false;
    bool storeQ_ = false;
    double absErr_, relErr_, vCritical_;
    int maxIter_;

    bool valid_ = false;
    const Fields *lastCurr_ = nullptr;
    int nx_ = 0, ny_ = 0, neta_ = 0;
    int nSteps_ = 0;
    double stepSeconds_ = 0., downloadSeconds_ = 0.;

    // the EoS descriptor and pressure table, in the backend's precision
    std::vector<unsigned char> eosDesc_;
    std::vector<unsigned char> eosTable_;

    std::string unsupportedReason(bool hasHydroSource) const;
    bool buildEos(std::string &err);
    void upload(const Fields &arena, double tau);
    void download(Fields &arena) const;
    template <typename Real>
    void uploadT(const Fields &arena, double tau, Real *s) const;
    template <typename Real>
    void downloadT(Fields &arena, const Real *s) const;
    template <typename Real>
    void fillParams(double tau, std::vector<unsigned char> &out) const;
    [[noreturn]] void fail(const std::string &what) const;
};

#endif  // SRC_DEVICE_DEVICE_HYDRO_H_
