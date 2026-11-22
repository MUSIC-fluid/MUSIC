// Copyright 2026 Chun Shen

#include "device_hydro.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <utility>

#include "../eos.h"
#include "../eos_idealgas.h"
#include "../fields.h"
#include "../util.h"
#include "device_api.h"
#include "device_visc.h"

#ifndef MUSIC_WITH_METAL
const MusicDeviceApi *music_metal_api(void) { return nullptr; }
#endif
#ifndef MUSIC_WITH_CUDA
const MusicDeviceApi *music_cuda_api(void) { return nullptr; }
#endif

namespace music_device {
namespace {

//! gpu_backend = 1: the device kernels in OpenMP loops on the host
template <typename Real>
class CpuBackend : public Backend {
  private:
    int nx_ = 0, ny_ = 0, neta_ = 0, n_ = 0, nComp_ = kNComp;
    // cur: tau_n; pred: the predictor; next: tau_{n-1} at the start of a
    // step (the previous slice of the viscous update), then tau_{n+1}
    std::vector<Real> slice_[3];
    int cur_ = 0, pred_ = 1, next_ = 2;
    std::vector<Real> pCur_, pPred_;
    std::vector<Real> flux_[3];
    std::vector<Real> wFace_[3];
    EosDesc<Real> desc_;
    std::vector<Real> table_;
    int dummyTable_ = 0;

  public:
    std::string name() const {
        return (sizeof(Real) == 8 ? "CPU (device kernels, double)"
                                  : "CPU (device kernels, float)");
    }
    int realSize() const { return (sizeof(Real)); }

    bool resize(int nx, int ny, int neta, int nComp, std::string &) {
        nx_ = nx;
        ny_ = ny;
        neta_ = neta;
        n_ = nx * ny * neta;
        nComp_ = nComp;
        for (auto &s : slice_) s.assign(static_cast<size_t>(nComp) * n_, Real(0));
        pCur_.assign(n_, Real(0));
        pPred_.assign(n_, Real(0));
        for (int d = 0; d < 3; d++) {
            const size_t nf = faceCount(nx, ny, neta, d + 1);
            flux_[d].assign(static_cast<size_t>(kNTJ) * nf, Real(0));
            if (nComp == kNCompVisc) {
                wFace_[d].assign(static_cast<size_t>(kNWFaceVals) * nf, Real(0));
            } else {
                wFace_[d].clear();
            }
        }
        cur_ = 0;
        pred_ = 1;
        next_ = 2;
        return (true);
    }

    bool setEos(
        const void *desc, size_t descBytes, const void *table,
        size_t tableBytes, std::string &err) {
        if (descBytes != sizeof(desc_)) {
            err = "EoS descriptor size mismatch";
            return (false);
        }
        std::memcpy(&desc_, desc, sizeof(desc_));
        table_.resize(tableBytes / sizeof(Real));
        if (tableBytes > 0) std::memcpy(table_.data(), table, tableBytes);
        return (true);
    }

    void *state(int which) {
        return (slice_[which == 0 ? cur_ : next_].data());
    }
    bool upload(std::string &) { return (true); }
    bool download(std::string &) { return (true); }

    bool step(
        const void *params, size_t paramBytes, int &flags, std::string &err) {
        if (paramBytes != sizeof(StepParams<Real>)) {
            err = "StepParams size mismatch";
            return (false);
        }
        StepParams<Real> P;
        std::memcpy(&P, params, sizeof(P));
        if (P.nComp != nComp_) {
            err = "the state slices have the wrong number of components";
            return (false);
        }
        const bool viscous = (nComp_ == kNCompVisc);
        const EosDesc<Real> &D = desc_;
        const Real *tab =
            (table_.empty() ? reinterpret_cast<const Real *>(&dummyTable_)
                            : table_.data());
        int f = 0;
        for (int rk = 0; rk < 2; rk++) {
            P.rkFlag = rk;
            const Real *sIn = slice_[rk == 0 ? cur_ : pred_].data();
            // the second stage averages with tau_n; the viscous update of
            // the first stage takes its time derivatives against tau_{n-1}
            const Real *sPrev = slice_[rk == 0 && viscous ? next_ : cur_].data();
            Real *pIn = (rk == 0 ? pCur_ : pPred_).data();
            Real *sOut = slice_[rk == 0 ? pred_ : next_].data();
            const int n = n_;
#pragma omp parallel
            {
#pragma omp for schedule(static)
                for (int i = 0; i < n; i++)
                    kernelPressure(P, D, tab, sIn, pIn, i);
                for (int d = 1; d <= 3; d++) {
                    const int nf = faceCount(nx_, ny_, neta_, d);
                    Real *fl = flux_[d - 1].data();
#pragma omp for schedule(dynamic, 64)
                    for (int i = 0; i < nf; i++)
                        kernelFace(P, D, tab, sIn, pIn, d, i, fl, &f);
                    if (viscous) {
                        Real *wf = wFace_[d - 1].data();
#pragma omp for schedule(dynamic, 64)
                        for (int i = 0; i < nf; i++)
                            kernelWFace(P, sIn, d, i, wf);
                    }
                }
                if (viscous) {
#pragma omp for schedule(dynamic, 64)
                    for (int i = 0; i < n; i++)
                        kernelCellVisc(
                            P, D, tab, sIn, sPrev, pIn, pCur_.data(),
                            flux_[0].data(), flux_[1].data(), flux_[2].data(),
                            wFace_[0].data(), wFace_[1].data(),
                            wFace_[2].data(), sOut, i, &f);
                } else {
#pragma omp for schedule(dynamic, 64)
                    for (int i = 0; i < n; i++)
                        kernelCell(
                            P, D, tab, sIn, sPrev, pIn, pCur_.data(),
                            flux_[0].data(), flux_[1].data(), flux_[2].data(),
                            sOut, i, &f);
                }
            }
        }
        std::swap(cur_, next_);
        flags = f;
        return (true);
    }
};

//! gpu_backend = 2, 3: a backend behind the C interface of device_api.h
class CApiBackend : public Backend {
  private:
    const MusicDeviceApi *api_;
    void *dev_;
    int realSize_;
    char err_[1024];

    bool check(int status, std::string &err) {
        if (status == 0) return (true);
        err = err_;
        return (false);
    }

  public:
    CApiBackend(const MusicDeviceApi *api, int realSize, std::string &err)
        : api_(api), dev_(nullptr), realSize_(realSize) {
        err_[0] = 0;
        dev_ = api_->create(realSize, err_, sizeof(err_));
        if (!dev_) err = err_;
    }
    ~CApiBackend() {
        if (dev_) api_->destroy(dev_);
    }
    bool ok() const { return (dev_ != nullptr); }

    std::string name() const { return (api_->name(dev_)); }
    int realSize() const { return (realSize_); }
    bool resize(int nx, int ny, int neta, int nComp, std::string &err) {
        return (check(
            api_->resize(dev_, nx, ny, neta, nComp, err_, sizeof(err_)), err));
    }
    bool setEos(
        const void *desc, size_t descBytes, const void *table,
        size_t tableBytes, std::string &err) {
        return (check(
            api_->setEos(
                dev_, desc, descBytes, table, tableBytes, err_, sizeof(err_)),
            err));
    }
    void *state(int which) { return (api_->state(dev_, which)); }
    bool upload(std::string &err) {
        return (check(api_->upload(dev_, err_, sizeof(err_)), err));
    }
    bool step(
        const void *params, size_t paramBytes, int &flags, std::string &err) {
        return (check(
            api_->step(dev_, params, paramBytes, &flags, err_, sizeof(err_)),
            err));
    }
    bool download(std::string &err) {
        return (check(api_->download(dev_, err_, sizeof(err_)), err));
    }
};

}  // namespace
}  // namespace music_device

using music_device::EosDesc;
using music_device::StepParams;

static_assert(
    offsetof(StepParams<float>, rkFlag) == MUSIC_DEVICE_RKFLAG_OFFSET
        && offsetof(StepParams<double>, rkFlag) == MUSIC_DEVICE_RKFLAG_OFFSET,
    "MUSIC_DEVICE_RKFLAG_OFFSET does not match StepParams");

DeviceHydro::DeviceHydro(
    const InitData &DATA_in, const EOS &eos, double absErr, double relErr,
    int maxIter, double vCritical, bool hasHydroSource)
    : DATA(DATA_in),
      eos_(eos),
      absErr_(absErr),
      relErr_(relErr),
      vCritical_(vCritical),
      maxIter_(maxIter) {
    const int which = DATA.gpu_backend;
    if (which == 0) return;
    const std::string reason = unsupportedReason(hasHydroSource);
    if (!reason.empty()) {
        music_message << "gpu_backend = " << which << ": " << reason
                      << "; the CPU code runs instead";
        music_message.flush("warning");
        return;
    }
    viscous_ = (DATA.viscosity_flag == 1);
    float_ = (which == 2 || DATA.gpu_float == 1);
    storeQ_ = (DATA.gpu_store_q < 0 ? float_ : DATA.gpu_store_q == 1);
    if (float_) {
        // MUSIC's beastMode-0 tolerances (1e-16, 1e-14) are below the float
        // resolution: use those of beastMode > 0
        absErr_ = std::max(absErr_, 1e-8);
        relErr_ = std::max(relErr_, 1e-6);
    }
    std::string err;
    if (which == 1) {
        if (float_) {
            backend_.reset(new music_device::CpuBackend<float>());
        } else {
            backend_.reset(new music_device::CpuBackend<double>());
        }
    } else if (which == 2 || which == 3) {
        const MusicDeviceApi *api =
            (which == 2 ? music_metal_api() : music_cuda_api());
        if (api == nullptr) {
            music_message << "gpu_backend = " << which << ": this MUSIC was "
                          << "built without "
                          << (which == 2 ? "Metal" : "CUDA")
                          << "; the CPU code runs instead";
            music_message.flush("warning");
            return;
        }
        music_device::CApiBackend *b = new music_device::CApiBackend(
            api, (float_ ? 4 : 8), err);
        if (!b->ok()) {
            delete b;
            fail("cannot start the GPU backend: " + err);
        }
        backend_.reset(b);
    } else {
        fail("unknown gpu_backend " + std::to_string(which));
    }
    if (!buildEos(err)) fail(err);
    music_message << "hydro on the device: " << backend_->name() << ", "
                  << (viscous_ ? "viscous, " : "ideal, ")
                  << (float_ ? "float" : "double")
                  << (storeQ_ ? ", stored tau T^{tau mu}" : "");
    music_message.flush("info");
}

DeviceHydro::~DeviceHydro() {
    if (backend_ == nullptr || nSteps_ == 0) return;
    music_message << "device: " << nSteps_ << " steps, "
                  << stepSeconds_ << " s in the kernels, "
                  << downloadSeconds_ << " s copying to the host Fields";
    music_message.flush("info");
}

std::string DeviceHydro::unsupportedReason(bool hasHydroSource) const {
    if (hasHydroSource) return ("source terms are not on the device yet");
    const int id = eos_.get_eos_id();
    if (!(id == 0 || id == 9 || id == 91 || (id >= 10 && id <= 14)))
        return ("EoS " + std::to_string(id) + " is not on the device yet");
    if (DATA.viscosity_flag == 0) return ("");
    if (DATA.viscosity_flag != 1)
        return ("viscosity_flag " + std::to_string(DATA.viscosity_flag)
                + " is not on the device");
    if (DATA.turn_on_diff == 1)
        return ("baryon diffusion is not on the device yet");
    if (DATA.turn_on_bulk_chem == 1)
        return ("the chemical bulk pressure is not on the device yet");
    if (!(id == 0 || id == 9 || id == 91))
        return ("viscous hydro with EoS " + std::to_string(id)
                + " is not on the device yet");
    if (DATA.FlagCausalityMethod == 1)
        return ("FlagCausalityMethod = 1 is not on the device yet");
    if (DATA.T_dependent_shear_to_s == 23 || DATA.T_dependent_bulk_to_s == 23)
        return ("transport coefficients from files are not on the device");
    return ("");
}

namespace {
template <typename Real>
void fillEos(
    const EOS &eos, bool withThermal, std::vector<unsigned char> &descOut,
    std::vector<unsigned char> &tableOut, std::string &err) {
    EosDesc<Real> d;
    std::memset(&d, 0, sizeof(d));
    const int id = eos.get_eos_id();
    d.type = (id == 0 ? 0 : (id == 9 || id == 91 ? 1 : 2));
    std::vector<Real> table;
    if (d.type != 0) {
        const EOS_base &b = eos.get_eos_base();
        d.nTables = b.get_number_of_tables();
        if (d.nTables > music_device::kMaxEosTables) {
            err = "the EoS has more tables than the device supports";
            return;
        }
        for (int t = 0; t < d.nTables; t++) {
            d.eLength[t] = b.e_length[t];
            d.nbLength[t] = b.nb_length[t];
            d.offset[t] = static_cast<int>(b.pressure_tb.offset(t));
            d.eBounds[t] = Real(b.e_bounds[t]);
            d.eSpacing[t] = Real(b.e_spacing[t]);
            d.nbBounds[t] = Real(b.nb_bounds[t]);
            d.nbSpacing[t] = Real(b.nb_spacing[t]);
        }
        const double *data = b.pressure_tb.data();
        table.assign(data, data + b.pressure_tb.size());
        if (withThermal) {
            // EOS_hotQCD::get_temperature: T^5 on the e grid of table 0
            d.tOffset = static_cast<int>(
                table.size() + b.temperature_tb.offset(0));
            const double *t = b.temperature_tb.data();
            table.insert(table.end(), t, t + b.temperature_tb.size());
        }
    } else if (withThermal) {
        // EOS_idealgas::get_muB_over_T
        const EOS_idealgas *g =
            dynamic_cast<const EOS_idealgas *>(&eos.get_eos_base());
        if (g == nullptr) {
            err = "EoS 0 is not the ideal gas";
            return;
        }
        const std::vector<double> &x = g->get_x_table();
        d.xOffset = static_cast<int>(table.size());
        d.xLength = static_cast<int>(x.size());
        table.insert(table.end(), x.begin(), x.end());
        d.coeffA = Real(g->get_coeff_a());
        d.coeffB = Real(g->get_coeff_b());
        d.rCap = Real(g->get_r_cap());
        d.xMax = Real(EOS_idealgas::get_muB_over_T_max());
    }
    descOut.resize(sizeof(d));
    std::memcpy(descOut.data(), &d, sizeof(d));
    tableOut.resize(table.size() * sizeof(Real));
    if (!table.empty())
        std::memcpy(tableOut.data(), table.data(), tableOut.size());
}
}  // namespace

bool DeviceHydro::buildEos(std::string &err) {
    if (float_) {
        fillEos<float>(eos_, viscous_, eosDesc_, eosTable_, err);
    } else {
        fillEos<double>(eos_, viscous_, eosDesc_, eosTable_, err);
    }
    if (!err.empty()) return (false);
    return (backend_->setEos(
        eosDesc_.data(), eosDesc_.size(),
        eosTable_.empty() ? nullptr : eosTable_.data(), eosTable_.size(),
        err));
}

template <typename Real>
void DeviceHydro::fillParams(
    double tau, std::vector<unsigned char> &out) const {
    StepParams<Real> P;
    std::memset(&P, 0, sizeof(P));
    P.nx = nx_;
    P.ny = ny_;
    P.neta = neta_;
    P.mode = DATA.face_reconstruction;
    P.rkFlag = 0;
    P.storeQ = (storeQ_ ? 1 : 0);
    P.maxIter = maxIter_;
    P.tau = Real(tau);
    P.dtau = Real(DATA.delta_tau);
    P.dx = Real(DATA.delta_x);
    P.dy = Real(DATA.delta_y);
    P.deta = Real(DATA.delta_eta);
    P.theta = Real(DATA.minmod_theta);
    P.absErr = Real(absErr_);
    P.relErr = Real(relErr_);
    P.vCritical = Real(vCritical_);
    // as in Advance::MakeDeltaQI
    const double delta3 = DATA.delta_eta;
    double coshDeta = cosh(delta3 / 2.) / std::max(delta3, Util::small_eps);
    double sinhDeta = sinh(delta3 / 2.) / std::max(delta3, Util::small_eps);
    sinhDeta = std::max(0.5, sinhDeta);
    if (DATA.boost_invariant) {
        coshDeta = 0.0;
        sinhDeta = 0.5;
    }
    P.coshDeta = Real(coshDeta);
    P.sinhDeta = Real(sinhDeta);
    P.nComp = (viscous_ ? music_device::kNCompVisc : music_device::kNComp);
    if (viscous_) fillViscParams(P);
    out.resize(sizeof(P));
    std::memcpy(out.data(), &P, sizeof(P));
}

template <typename Real>
void DeviceHydro::fillViscParams(StepParams<Real> &P) const {
    P.shear = DATA.turn_on_shear;
    P.bulk = DATA.turn_on_bulk;
    P.secondOrder = DATA.include_second_order_terms;
    P.vorticity = (DATA.include_vorticity_terms ? 1 : 0);
    P.duIdeal = DATA.Du_from_ideal_EoM;
    P.resum = (DATA.FlagResumTransportCoeff ? 1 : 0);
    P.muBShear = DATA.muB_dependent_shear_to_s;
    P.bulkRelaxType = DATA.bulk_relaxation_type;
    P.questRevert = (DATA.Initial_profile != 0 && DATA.Initial_profile != 1);
    P.shearT = DATA.T_dependent_shear_to_s;
    P.bulkT = DATA.T_dependent_bulk_to_s;
    // Diss::nWComponents (no chemical bulk, no diffusion)
    P.nWComp = (DATA.turn_on_bulk == 1 ? 6 : 5);
    for (int a = 0; a < 4; a++)
        for (int b = 0; b < 4; b++) P.g[a][b] = Real(DATA.gmunu[a][b]);
    P.duIdealEDilute = Real(DATA.Du_from_ideal_EoM_e_dilute);
    P.resumAlpha = Real(DATA.resumTransCoeffAlpha);
    P.questRevertStrength = Real(DATA.quest_revert_strength);
    P.shearRelaxFactor = Real(DATA.shear_relax_time_factor);
    P.bulkRelaxFactor = Real(DATA.bulk_relax_time_factor);
    P.shearToS = Real(DATA.shear_to_s);
    P.shear2Min = Real(DATA.shear_2_min);
    P.shear2Slope = Real(DATA.shear_2_slope);
    P.shear2Curv = Real(DATA.shear_2_curv);
    P.shear3TKink = Real(DATA.shear_3_T_kink_in_GeV);
    P.shear3LowSlope = Real(DATA.shear_3_low_T_slope_in_GeV);
    P.shear3HighSlope = Real(DATA.shear_3_high_T_slope_in_GeV);
    P.shear3AtKink = Real(DATA.shear_3_at_kink);
    P.shearMuBf0p2 = Real(DATA.shear_muBf0p2);
    P.shearMuBf0p4 = Real(DATA.shear_muBf0p4);
    P.shearMuBAlpha = Real(DATA.shear_muBDep_alpha);
    P.shearMuBSlope = Real(DATA.shear_muBDep_slope);
    P.shearMuBScale = Real(DATA.shear_muBDep_scale);
    P.bulk2Norm = Real(DATA.bulk_2_normalisation);
    P.bulk2Width = Real(DATA.bulk_2_width_in_GeV);
    P.bulk2Peak = Real(DATA.bulk_2_peak_in_GeV);
    P.bulk3Max = Real(DATA.bulk_3_max);
    P.bulk3Width = Real(DATA.bulk_3_width_in_GeV);
    P.bulk3TPeak = Real(DATA.bulk_3_T_peak_in_GeV);
    P.bulk3Lambda = Real(DATA.bulk_3_lambda_asymm);
    P.bulk10Max = Real(DATA.bulk_10_max);
    P.bulk10MaxMuB0p2 = Real(DATA.bulk_10_max_muB0p2);
    P.bulk10MaxMuB0p4 = Real(DATA.bulk_10_max_muB0p4);
    P.bulk10Tpeak = Real(DATA.bulk_10_Tpeak);
    P.bulk10TpeakMuBCurv = Real(DATA.bulk_10_Tpeak_muBcurv);
    P.bulk10WidthLow = Real(DATA.bulk_10_width_low);
    P.bulk10WidthHigh = Real(DATA.bulk_10_width_high);
}

template <typename Real>
void DeviceHydro::uploadT(const Fields &arena, double tau, Real *s) const {
    const int n = nx_ * ny_ * neta_;
    const bool storeQ = storeQ_;
#pragma omp parallel for collapse(2)
    for (int ieta = 0; ieta < neta_; ieta++)
        for (int iy = 0; iy < ny_; iy++)
            for (int ix = 0; ix < nx_; ix++) {
                const int idx = arena.getFieldIdx(ix, iy, ieta);
                const int i = ix + nx_ * (iy + ny_ * ieta);
                music_device::Cell<double> c;
                c.e = arena.e_[idx];
                c.rhob = arena.rhob_[idx];
                c.rhoq = arena.rhoq_[idx];
                c.rhos = arena.rhos_[idx];
                for (int k = 0; k < 4; k++) c.u[k] = arena.u_[k][idx];
                s[i] = Real(c.e);
                s[n + i] = Real(c.rhob);
                s[2 * n + i] = Real(c.rhoq);
                s[3 * n + i] = Real(c.rhos);
                for (int k = 0; k < 4; k++) s[(4 + k) * n + i] = Real(c.u[k]);
                if (viscous_) {
                    for (int k = 0; k < 14; k++)
                        s[(music_device::kOffW + k) * n + i] =
                            Real(arena.Wmunu_[k][idx]);
                    s[music_device::kOffPi * n + i] = Real(arena.piBulk_[idx]);
                    s[music_device::kOffPiChem * n + i] =
                        Real(arena.piBulkChem_[idx]);
                }
                if (storeQ) {
                    const double p =
                        eos_.get_pressure(c.e, c.rhob, c.rhoq, c.rhos);
                    for (int a = 0; a < music_device::kNTJ; a++)
                        s[(music_device::kNPrim + a) * n + i] =
                            Real(music_device::TJb(c, a, 0, p) * tau);
                }
            }
}

template <typename Real>
void DeviceHydro::downloadT(Fields &arena, const Real *s) const {
    const int n = nx_ * ny_ * neta_;
#pragma omp parallel for collapse(2)
    for (int ieta = 0; ieta < neta_; ieta++)
        for (int iy = 0; iy < ny_; iy++)
            for (int ix = 0; ix < nx_; ix++) {
                const int idx = arena.getFieldIdx(ix, iy, ieta);
                const int i = ix + nx_ * (iy + ny_ * ieta);
                arena.e_[idx] = s[i];
                arena.rhob_[idx] = s[n + i];
                arena.rhoq_[idx] = s[2 * n + i];
                arena.rhos_[idx] = s[3 * n + i];
                for (int k = 0; k < 4; k++) arena.u_[k][idx] = s[(4 + k) * n + i];
                if (viscous_) {
                    for (int k = 0; k < 14; k++)
                        arena.Wmunu_[k][idx] = s[(music_device::kOffW + k) * n + i];
                    arena.piBulk_[idx] = s[music_device::kOffPi * n + i];
                    arena.piBulkChem_[idx] = s[music_device::kOffPiChem * n + i];
                }
            }
}

void DeviceHydro::upload(const Fields &curr, const Fields &prev, double tau) {
    std::string err;
    if (curr.nX() != nx_ || curr.nY() != ny_ || curr.nEta() != neta_) {
        nx_ = curr.nX();
        ny_ = curr.nY();
        neta_ = curr.nEta();
        const int nComp =
            (viscous_ ? music_device::kNCompVisc : music_device::kNComp);
        if (!backend_->resize(nx_, ny_, neta_, nComp, err)) fail(err);
    }
    // the viscous update also reads the previous slice (tau_{n-1}): its tau
    // T^{tau mu} is never read
    const int nSlices = (viscous_ ? 2 : 1);
    for (int k = 0; k < nSlices; k++) {
        const Fields &arena = (k == 0 ? curr : prev);
        if (float_) {
            uploadT(arena, tau, static_cast<float *>(backend_->state(k)));
        } else {
            uploadT(arena, tau, static_cast<double *>(backend_->state(k)));
        }
    }
    if (!backend_->upload(err)) fail(err);
}

void DeviceHydro::download(Fields &arena) const {
    if (float_) {
        downloadT(arena, static_cast<const float *>(backend_->state(0)));
    } else {
        downloadT(arena, static_cast<const double *>(backend_->state(0)));
    }
}

void DeviceHydro::advanceStep(
    double tau, Fields *&fpPrev, Fields *&fpCurr, Fields *&fpNext) {
    if (!valid_ || fpCurr != lastCurr_ || (viscous_ && fpPrev != lastPrev_)
        || fpCurr->nX() != nx_ || fpCurr->nY() != ny_
        || fpCurr->nEta() != neta_) {
        if (viscous_ && (fpPrev->nX() != fpCurr->nX()
                         || fpPrev->nY() != fpCurr->nY()
                         || fpPrev->nEta() != fpCurr->nEta())) {
            fail("the previous and current slices differ in size");
        }
        upload(*fpCurr, *fpPrev, tau);
    }
    std::vector<unsigned char> params;
    if (float_) {
        fillParams<float>(tau, params);
    } else {
        fillParams<double>(tau, params);
    }
    int flags = 0;
    std::string err;
    const auto t0 = std::chrono::steady_clock::now();
    if (!backend_->step(params.data(), params.size(), flags, err)) fail(err);
    const auto t1 = std::chrono::steady_clock::now();
    if (flags & music_device::kFatalFlags) {
        fail("MaxSpeed: the signal speed of a KT face is out of range "
             "(where the CPU code exits)");
    }
    if (DATA.echo_level > 5 && (flags & ~music_device::kFatalFlags)) {
        music_message << "device step flags: " << flags
                      << " (2 revert, 4 regulate, 8 momentum cap, 16 "
                      << "velocity bracket failed)";
        music_message.flush("info");
    }
    if (!backend_->download(err)) fail(err);
    // Evolve::AdvanceRK: after both stages the slice that was the previous
    // one holds the new state, the current one becomes the previous one
    Fields *oldPrev = fpPrev;
    fpPrev = fpCurr;
    fpCurr = oldPrev;
    (void)fpNext;
    download(*fpCurr);
    fpCurr->updateGhosts();
    lastCurr_ = fpCurr;
    lastPrev_ = fpPrev;
    valid_ = true;
    const auto t2 = std::chrono::steady_clock::now();
    stepSeconds_ += std::chrono::duration<double>(t1 - t0).count();
    downloadSeconds_ += std::chrono::duration<double>(t2 - t1).count();
    nSteps_++;
}

void DeviceHydro::fail(const std::string &what) const {
    music_message << "gpu_backend: " << what;
    music_message.flush("error");
    exit(1);
}
