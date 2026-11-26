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
// The state stays on the device between steps. After a step the host
// Fields are not filled: Evolve asks for what it reads (syncHost for all
// of the current and the previous slice on output and freeze-out steps,
// syncEtaWindow and syncShell for the diagnostics that read part of the
// grid), and takes the conservation sums and the grid maxima
// (conservationRows) and the cells of the boundary outflow (outflowFaces)
// from the device. Each filled cell holds what AdvanceRK would
// have written. The host Fields are uploaded again only when the device
// copy is stale: the first step, a new slice pointer or grid, and after
// invalidate(), which Evolve calls when it changes the fields itself
// (beastMode regrid, change of delta_tau), after syncHost.
//
// Supported, with or without source terms: ideal hydro (viscosity_flag 0)
// and viscous hydro (viscosity_flag 1, shear, bulk, baryon diffusion,
// Pi_chem, FlagResetCausality) with EoS 0, 9, 91, 10 - 14 and 20, except
// the transport coefficients from files (type 23). The strings source
// (HydroSourceStrings) is evaluated on the device from the strings of the
// tau frame, sent once per step; the other sources on the host for each
// RK stage from the u of the device slice, added on the device.
// Otherwise the CPU code runs and a warning says why.

#ifndef SRC_DEVICE_DEVICE_HYDRO_H_
#define SRC_DEVICE_DEVICE_HYDRO_H_

#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "../data.h"
#include "../pretty_ostream.h"

class EOS;
class Fields;
class HydroSourceBase;
class HydroSourceStrings;
namespace music_device {
template <typename Real>
struct StepParams;
}

namespace music_device {

//! one device: see device_api.h, here as a C++ interface
class Backend {
  public:
    virtual ~Backend() {}
    virtual std::string name() const = 0;
    virtual int realSize() const = 0;
    //! nComp: kNComp, or kNCompVisc for the viscous update
    virtual bool resize(
        int nx, int ny, int neta, int nComp, std::string &err) = 0;
    virtual bool setEos(
        const void *desc, size_t descBytes, const void *table,
        size_t tableBytes, std::string &err) = 0;
    //! which = 0: the current slice; 1: the previous one (tau_{n-1}, read
    //! by the viscous update); 2: the predictor
    virtual void *state(int which) = 0;
    virtual bool upload(std::string &err) = 0;
    //! one RK stage (rk = 0, 1) of the Heun step: see device_api.h
    //! the flags of a stage the backend has not finished come from finish()
    virtual bool stage(
        const void *params, size_t paramBytes, int rk, int &flags,
        std::string &err) = 0;
    //! waits for the second stage of the step (the backends may return
    //! from stage() before it has run); flags: those not returned yet
    virtual bool finish(int &flags, std::string &) {
        flags = 0;
        return (true);
    }
    //! the source terms of the next stage: kNTJ components of the cells
    virtual void *source() = 0;
    virtual bool uploadSource(std::string &err) = 0;
    //! makes components [firstComp, firstComp + nComp) of slice which
    //! current in state(which)
    virtual bool fetch(
        int which, int firstComp, int nComp, std::string &err) = 0;
    //! all components of the listed cells of slice which, gathered into
    //! *out (component-major, see device_api.h). False with an empty err:
    //! not supported, use fetch()
    virtual bool fetchCells(
        int, const std::vector<int> &, const void *&out, std::string &err) {
        out = nullptr;
        err.clear();
        return (false);
    }
    virtual bool download(std::string &err) = 0;
    //! the packed strings (device_strings.h) for the stages with
    //! StepParams::stringSource
    virtual bool setStrings(
        const void *data, size_t bytes, std::string &err) = 0;
    //! kernelConservationRow for every (eta, x) row of the current slice
    //! with the previous one (device_diag.h)
    virtual bool conservation(
        const void *params, size_t paramBytes, void *out, size_t outBytes,
        std::string &err) = 0;
    //! kernelCoarseningRow for every (eta, x) row of the current slice
    //! (device_diag.h)
    virtual bool coarsening(
        double eThreshold, void *out, size_t outBytes, std::string &err) = 0;
    //! kernelEdgeMaxRow for every (eta, x) row of the current slice
    virtual bool edgeMax(
        int bandX, int bandY, void *out, size_t outBytes,
        std::string &err) = 0;
    //! the switches of the run (DeviceHydro::specKey); the CPU kernels
    //! read them at run time
    virtual bool specialize(const std::string &, std::string &) {
        return (true);
    }
    //! fetch() of the slice of tau_n (which = 1) while the step that
    //! stage() started may still run: see MusicDeviceApi::fetchUnwritten
    virtual bool fetchUnwritten(
        int which, int firstComp, int nComp, std::string &err) {
        return (fetch(which, firstComp, nComp, err));
    }
    //! conservation() of the slices after a step, in two halves (see
    //! MusicDeviceApi::conservationBegin): here computed at once and kept
    virtual bool conservationBegin(
        const void *params, size_t paramBytes, size_t outBytes,
        std::string &err) {
        consKept_.resize(outBytes);
        return (conservation(
            params, paramBytes, consKept_.data(), outBytes, err));
    }
    virtual bool conservationEnd(void *out, size_t outBytes, std::string &err) {
        if (outBytes != consKept_.size()) {
            err = "conservationEnd: not the size of conservationBegin";
            return (false);
        }
        std::memcpy(out, consKept_.data(), outBytes);
        return (true);
    }

  private:
    std::vector<unsigned char> consKept_;
};

}  // namespace music_device

class DeviceHydro {
  public:
    //! source: the source terms Advance adds (null if none); the host
    //! evaluates them for each RK stage, the device adds them
    DeviceHydro(
        const InitData &DATA, const EOS &eos, double absErr, double relErr,
        int maxIter, double vCritical,
        std::shared_ptr<HydroSourceBase> source);
    ~DeviceHydro();

    //! the device advances the fluid
    bool active() const { return (backend_ != nullptr); }

    //! one Heun step from *fpCurr at tau; rotates the pointers as
    //! Evolve::AdvanceRK does. Finishes the step beginStep started, if any
    void advanceStep(
        double tau, Fields *&fpPrev, Fields *&fpCurr, Fields *&fpNext);
    //! starts the step from *fpCurr at tau on the device and returns while
    //! it runs: the host may read and write its Fields (not the device)
    //! until advanceStep finishes the step
    void beginStep(double tau, Fields *fpPrev, Fields *fpCurr);
    //! a step is running (beginStep without advanceStep)
    bool stepPending() const { return (pending_); }

    //! the host changed the fields: upload before the next step
    void invalidate() {
        valid_ = false;
        consAheadValid_ = false;
    }

    //! the host Fields hold the state of the last step
    bool hostCurrent() const { return (hostCurrent_); }
    //! fills curr and prev (the slices of the last step, as advanceStep
    //! left them) from the device, ghost cells included
    void syncHost(Fields &curr, Fields &prev);
    //! fills curr, the Fields of tau_n that beginStep started from, while
    //! the step runs (the stages read that slice and write the others;
    //! with Metal the copy overlaps the step), ghost cells included. For
    //! the host work of a step that reads all of fpCurr but not fpPrev
    void syncCurrentDuringStep(Fields &curr);
    //! fills the cells of curr in the eta window (eta_min, eta_max) of the
    //! Cell_info diagnostics (all eta slices if boost-invariant)
    void syncEtaWindow(Fields &curr, double etaMin, double etaMax);
    //! fills the cells of curr (and prev) within bandX (bandY) cells of
    //! the x (y) edges, and with bandEta > 0 the first and last bandEta
    //! eta slices
    void syncShell(Fields &curr, Fields *prev, int bandX, int bandY,
                   int bandEta);
    //! the values of the cells at the outflow faces of slice which (0:
    //! the current one, whose host Fields are curr; 1: the previous one) in
    //! the order of Cell_info::outflowFace, for boundary_outflow_rate:
    //! what it would read from Fields holding the slice. The previous
    //! slice's values are those of the current slice of the step before,
    //! when they were taken then. Valid until the next step
    const std::vector<double> &outflowFaces(const Fields &arena, int which);
    //! the rows of Cell_info::check_conservation_law (kConsVals per (eta,
    //! x) row, see device_diag.h) of the device state at tau: those the
    //! last step computed (conservationAhead) if it computed them for tau,
    //! else computed now
    void conservationRows(double tau, std::vector<double> &rows);
    //! each step computes the conservation rows of its new slice at
    //! tau + delta_tau right after its second stage, on the device, so
    //! that conservationRows of the next iteration does not wait for a
    //! kernel (Evolve: in beastMode 0, which reads them every step)
    void conservationAhead(bool on) { consAhead_ = on; }
    //! how many conservationRows calls took the rows of conservationAhead
    int conservationRowsAhead() const { return (nConsAhead_); }
    //! Evolve::coarseningError of the device state, reduced on the device:
    //! the host receives kCoarsenVals values per (eta, x) row, not the
    //! slice. Sums: e, e de, e du; maxima: de, du; the cell count
    void coarseningSums(double eThreshold, double sums[6]);
    //! the largest e of the current slice of the last step (the Fields
    //! curr) within bandX (bandY) cells of the x (y) edges, reduced on the
    //! device (Evolve::fluidNearGridEdge, which used syncShell)
    double edgeMaxE(const Fields &curr, int bandX, int bandY);

  private:
    const InitData &DATA;
    const EOS &eos_;
    mutable pretty_ostream music_message{std::string("🚀 DeviceHydro")};
    std::unique_ptr<music_device::Backend> backend_;
    std::shared_ptr<HydroSourceBase> source_;
    bool float_ = false;
    bool viscous_ = false;
    bool storeQ_ = false;
    double absErr_, relErr_, vCritical_;
    int maxIter_;

    bool valid_ = false;
    bool hostCurrent_ = true;
    // the step beginStep started
    bool pending_ = false;
    double pendingTau_ = 0.;
    const Fields *pendingCurr_ = nullptr, *pendingPrev_ = nullptr;
    int pendingFlags_ = 0;
    //! fails if a step is running: the device state is being written
    void requireIdle(const char *what) const;
    const HydroSourceStrings *strings_ = nullptr;
    const Fields *lastCurr_ = nullptr;
    const Fields *lastPrev_ = nullptr;
    int nx_ = 0, ny_ = 0, neta_ = 0;
    int nSteps_ = 0;
    double stepSeconds_ = 0., downloadSeconds_ = 0., sourceSeconds_ = 0.;
    double diagSeconds_ = 0.;
    double coarsenSeconds_ = 0.;
    int nCoarsenChecks_ = 0;
    double edgeSeconds_ = 0.;
    int nEdgeChecks_ = 0;
    int nFullSyncs_ = 0, nPartialSyncs_ = 0, nOverlapSyncs_ = 0;
    // conservationAhead: the rows the last step queued, and their tau
    bool consAhead_ = false;
    bool consAheadValid_ = false;
    double consAheadTau_ = 0.;
    int nConsAhead_ = 0;
    //! the backend slice of the Fields of tau_n for a read, and whether it
    //! is read during a step (then slice 1, through fetchUnwritten): curr
    //! must be the current Fields of the last step, or the one the running
    //! step started from
    int currentSlice(const Fields &curr, const char *what, bool &during) const;
    // outflowFaces of the current and the previous slice, and whether they
    // hold the device state (after a step the current one becomes the
    // previous one)
    std::vector<double> outflow_[2];
    bool outflowValid_[2] = {false, false};
    // ideal update: outflow_ holds the Fields' Pi and W, which only the
    // host changes (and then uploads)
    bool outflowStatic_[2] = {false, false};

    // the EoS descriptor and pressure table, in the backend's precision
    std::vector<unsigned char> eosDesc_;
    std::vector<unsigned char> eosTable_;

    std::string unsupportedReason() const;
    //! qi_source of Advance::FirstRKStepT for every cell of the slice s at
    //! tau_rk, into the backend's source buffer; false if all are zero
    template <typename Real>
    bool computeSource(const Real *s, double tauRk);
    bool buildEos(std::string &err);
    void upload(const Fields &curr, const Fields &prev, double tau);
    //! the cells of slice which (0: current, 1: previous) into arena:
    //! keep(iy, ieta, r) gives the x ranges [r[0], r[1]), [r[2], r[3]) of
    //! the row and returns their number (0 - 2)
    template <typename Keep>
    void downloadCells(
        Fields &arena, int which, const Keep &keep, bool during = false);
    //! the strings of the current tau frame, packed for setStrings
    template <typename Real>
    void packStrings(std::vector<unsigned char> &out, int &nActive) const;
    template <typename Real>
    void uploadT(const Fields &arena, double tau, Real *s) const;
    template <typename Real>
    void fillViscParams(music_device::StepParams<Real> &P) const;
    //! the StepParams and EosDesc fields that MD_RUN and MD_EOS_TYPE read,
    //! as "name=value;...": fixed for the run
    std::string specKey() const;
    template <typename Real, typename Keep>
    void downloadT(Fields &arena, const Real *s, const Keep &keep) const;
    //! the cell at Fields index idx from s, component c at s[c stride + at]
    template <typename Real>
    void storeCell(
        Fields &arena, int idx, const Real *s, size_t stride,
        size_t at) const;
    //! the cells gathered by Backend::fetchCells into arena
    template <typename Real>
    void downloadGathered(
        Fields &arena, const Real *g, const std::vector<int> &cells) const;
    std::vector<int> cells_;  // the cells of the last partial download
    //! the cells of the outflow faces (Cell_info::outflowFace order)
    std::vector<int> outflowCells_;
    int outflowCellsNx_ = 0, outflowCellsNy_ = 0, outflowCellsNeta_ = 0;
    template <typename Real>
    //! skipVacuum: StepParams::eVacuum = vacuumSkipEnergy, else 0
    void fillParams(
        double tau, bool hasSource, std::vector<unsigned char> &out,
        bool skipVacuum = false) const;
    [[noreturn]] void fail(const std::string &what) const;
};

#endif  // SRC_DEVICE_DEVICE_HYDRO_H_
