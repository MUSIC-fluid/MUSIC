// Copyright 2026 Chun Shen

#include "device_hydro.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <utility>

#include "../eos.h"
#include "../eos_4D.h"
#include "../eos_idealgas.h"
#include "../fields.h"
#include "../grid_info.h"
#include "../hydro_source_base.h"
#include "../hydro_source_strings.h"
#include "../util.h"
#include "device_api.h"
#include "device_diag.h"
#include "device_strings.h"
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
    std::vector<Real> src_;
    std::vector<Real> strings_;
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
        src_.assign(static_cast<size_t>(kNTJ) * n_, Real(0));
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
        return (slice_[which == 0 ? cur_ : (which == 1 ? next_ : pred_)].data());
    }
    bool upload(std::string &) { return (true); }
    bool download(std::string &) { return (true); }
    void *source() { return (src_.data()); }
    bool uploadSource(std::string &) { return (true); }
    bool fetch(int, int, int, std::string &) { return (true); }
    bool setStrings(const void *data, size_t bytes, std::string &) {
        strings_.resize(bytes / sizeof(Real));
        if (bytes > 0) std::memcpy(strings_.data(), data, bytes);
        return (true);
    }
    bool conservation(
        const void *params, size_t paramBytes, void *out, size_t outBytes,
        std::string &err) {
        if (paramBytes != sizeof(StepParams<Real>)) {
            err = "StepParams size mismatch";
            return (false);
        }
        StepParams<Real> P;
        std::memcpy(&P, params, sizeof(P));
        const int rows = nx_ * neta_;
        if (outBytes < sizeof(Real) * kConsVals * rows) {
            err = "the conservation buffer is too small";
            return (false);
        }
        const Real *tab =
            (table_.empty() ? reinterpret_cast<const Real *>(&dummyTable_)
                            : table_.data());
        Real *o = static_cast<Real *>(out);
        const Real *sCurr = slice_[cur_].data();
        const Real *sPrev = slice_[next_].data();
#pragma omp parallel for schedule(static)
        for (int row = 0; row < rows; row++)
            kernelConservationRow(P, desc_, tab, sCurr, sPrev, o, row);
        return (true);
    }

    bool stage(
        const void *params, size_t paramBytes, int rk, int &flags,
        std::string &err) {
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
        {
            P.rkFlag = rk;
            const Real *sIn = slice_[rk == 0 ? cur_ : pred_].data();
            // the second stage averages with tau_n; the viscous update of
            // the first stage takes its time derivatives against tau_{n-1}
            const Real *sPrev = slice_[rk == 0 && viscous ? next_ : cur_].data();
            Real *pIn = (rk == 0 ? pCur_ : pPred_).data();
            Real *sOut = slice_[rk == 0 ? pred_ : next_].data();
            const int n = n_;
            if (P.stringSource && strings_.empty()) {
                err = "the strings were not set";
                return (false);
            }
#pragma omp parallel
            {
                if (P.stringSource) {
                    // the strings source from the u of the stage's input
#pragma omp for schedule(dynamic, 64)
                    for (int i = 0; i < n; i++)
                        kernelStringSource(
                            P, strings_.data(), sIn, src_.data(), i, &f);
                }
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
                            wFace_[2].data(), src_.data(), sOut, i, &f);
                } else {
#pragma omp for schedule(dynamic, 64)
                    for (int i = 0; i < n; i++)
                        kernelCell(
                            P, D, tab, sIn, sPrev, pIn, pCur_.data(),
                            flux_[0].data(), flux_[1].data(), flux_[2].data(),
                            src_.data(), sOut, i, &f);
                }
            }
        }
        // after the second stage the new slice is the current one
        if (rk == 1) std::swap(cur_, next_);
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
    bool finish(int &flags, std::string &err) {
        flags = 0;
        return (check(api_->finish(dev_, &flags, err_, sizeof(err_)), err));
    }
    bool stage(
        const void *params, size_t paramBytes, int rk, int &flags,
        std::string &err) {
        return (check(
            api_->stage(
                dev_, params, paramBytes, rk, &flags, err_, sizeof(err_)),
            err));
    }
    void *source() { return (api_->source(dev_)); }
    bool uploadSource(std::string &err) {
        return (check(api_->uploadSource(dev_, err_, sizeof(err_)), err));
    }
    bool fetch(int which, int firstComp, int nComp, std::string &err) {
        return (check(
            api_->fetch(
                dev_, which, firstComp, nComp, err_, sizeof(err_)),
            err));
    }
    bool fetchUnwritten(
        int which, int firstComp, int nComp, std::string &err) {
        return (check(
            api_->fetchUnwritten(
                dev_, which, firstComp, nComp, err_, sizeof(err_)),
            err));
    }
    bool conservationBegin(
        const void *params, size_t paramBytes, size_t outBytes,
        std::string &err) {
        if (api_->conservationBegin == nullptr)
            return (Backend::conservationBegin(
                params, paramBytes, outBytes, err));
        return (check(
            api_->conservationBegin(
                dev_, params, paramBytes, outBytes, err_, sizeof(err_)),
            err));
    }
    bool conservationEnd(void *out, size_t outBytes, std::string &err) {
        if (api_->conservationEnd == nullptr)
            return (Backend::conservationEnd(out, outBytes, err));
        return (check(
            api_->conservationEnd(dev_, out, outBytes, err_, sizeof(err_)),
            err));
    }
    bool download(std::string &err) {
        return (check(api_->download(dev_, err_, sizeof(err_)), err));
    }
    bool fetchCells(
        int which, const std::vector<int> &cells, const void *&out,
        std::string &err) {
        out = nullptr;
        err.clear();
        if (api_->fetchCells == nullptr) return (false);
        void *g = nullptr;
        if (!check(api_->fetchCells(dev_, which, cells.data(),
                                    static_cast<int>(cells.size()), &g, err_,
                                    sizeof(err_)),
                   err))
            return (false);
        out = g;
        return (true);
    }
    bool setStrings(const void *data, size_t bytes, std::string &err) {
        return (check(
            api_->setStrings(dev_, data, bytes, err_, sizeof(err_)), err));
    }
    bool specialize(const std::string &key, std::string &err) {
        return (check(
            api_->specialize(dev_, key.c_str(), err_, sizeof(err_)), err));
    }
    bool conservation(
        const void *params, size_t paramBytes, void *out, size_t outBytes,
        std::string &err) {
        return (check(
            api_->conservation(
                dev_, params, paramBytes, out, outBytes, err_, sizeof(err_)),
            err));
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
static_assert(
    offsetof(StepParams<float>, stringSource) == MUSIC_DEVICE_STRINGS_OFFSET
        && offsetof(StepParams<double>, stringSource)
               == MUSIC_DEVICE_STRINGS_OFFSET,
    "MUSIC_DEVICE_STRINGS_OFFSET does not match StepParams");

DeviceHydro::DeviceHydro(
    const InitData &DATA_in, const EOS &eos, double absErr, double relErr,
    int maxIter, double vCritical, std::shared_ptr<HydroSourceBase> source)
    : DATA(DATA_in),
      eos_(eos),
      source_(source),
      absErr_(absErr),
      relErr_(relErr),
      vCritical_(vCritical),
      maxIter_(maxIter) {
    const int which = DATA.gpu_backend;
    if (which == 0) return;
    const std::string reason = unsupportedReason();
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
    if (!backend_->specialize(specKey(), err)) fail(err);
    strings_ = dynamic_cast<const HydroSourceStrings *>(source_.get());
    music_message << "hydro on the device: " << backend_->name() << ", "
                  << (viscous_ ? "viscous, " : "ideal, ")
                  << (float_ ? "float" : "double")
                  << (storeQ_ ? ", stored tau T^{tau mu}" : "")
                  << (strings_ ? ", strings source on the device"
                               : (source_ ? ", source terms from the host"
                                          : ""));
    music_message.flush("info");
}

DeviceHydro::~DeviceHydro() {
    if (backend_ != nullptr && pending_) {
        // a step Evolve did not finish (it returned to re-run): wait
        int f = 0;
        std::string err;
        backend_->finish(f, err);
        pending_ = false;
    }
    if (backend_ == nullptr || nSteps_ == 0) return;
    music_message << "device: " << nSteps_ << " steps, "
                  << stepSeconds_
                  << " s launching and waiting for the device (the steps "
                     "overlap the host work), "
                  << downloadSeconds_ << " s copying to the host Fields ("
                  << nFullSyncs_ << " full, " << nPartialSyncs_
                  << " partial, " << nOverlapSyncs_
                  << " of the current slice during the step), "
                  << diagSeconds_ << " s for the conservation sums ("
                  << nConsAhead_ << " computed right after their step)";
    if (source_)
        music_message << ", " << sourceSeconds_
                      << (strings_ ? " s packing the strings"
                                   : " s for the source terms on the host");
    music_message.flush("info");
}

std::string DeviceHydro::unsupportedReason() const {
    const int id = eos_.get_eos_id();
    if (!(id == 0 || id == 9 || id == 91 || (id >= 10 && id <= 14)
          || id == 20))
        return ("EoS " + std::to_string(id) + " is not on the device yet");
    if (DATA.viscosity_flag == 0) return ("");
    if (DATA.viscosity_flag != 1)
        return ("viscosity_flag " + std::to_string(DATA.viscosity_flag)
                + " is not on the device");
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
    d.type = (id == 0 ? 0 : (id == 9 || id == 91 ? 1 : (id == 20 ? 3 : 2)));
    std::vector<Real> table;
    if (d.type == 3) {
        // EOS_4D: T5diff_ in Real, then the float tables, written straight
        // into tableOut (P alone is 73 MB)
        const EOS_4D *e4 = dynamic_cast<const EOS_4D *>(&eos.get_eos_base());
        if (e4 == nullptr) {
            err = "EoS 20 is not EOS_4D";
            return;
        }
        float x0[4], dx[4], xMax[4];
        int nn[4];
        e4->get_grid(x0, dx, nn, xMax);
        d.T04 = x0[0];
        d.mub04 = x0[1];
        d.muq04 = x0[2];
        d.mus04 = x0[3];
        d.dT4 = dx[0];
        d.dmub4 = dx[1];
        d.dmuq4 = dx[2];
        d.dmus4 = dx[3];
        d.nT4 = nn[0];
        d.nMub4 = nn[1];
        d.nMuq4 = nn[2];
        d.nMus4 = nn[3];
        d.Tmax4 = xMax[0];
        d.mubMax4 = xMax[1];
        d.muqMax4 = xMax[2];
        d.musMax4 = xMax[3];
        d.TrefDilute4 = EOS_4D::get_Ttilde_dilute();
        d.alphaNf4 = Real(e4->get_alphaNf());
        d.oneOverAlphaNf4 = Real(e4->get_one_over_alphaNf());
        const std::vector<double> &t5 = e4->get_T5diff();
        const size_t size = e4->get_table(0).size();
        // ideal hydro: P, and T for the grid maxima of the diagnostics
        const size_t nTab = (withThermal ? 5 : 2);
        if (size != static_cast<size_t>(nn[0]) * nn[1] * nn[2] * nn[3]
            || nTab * size > static_cast<size_t>(INT_MAX)) {
            err = "the 4D EoS tables do not fit the device index";
            return;
        }
        for (size_t k = 1; k < nTab; k++) {
            if (e4->get_table(static_cast<int>(k)).size() != size) {
                err = "the 4D EoS tables differ in size";
                return;
            }
        }
        d.t5Offset = 0;
        d.f4Offset = static_cast<int>(t5.size());
        d.f4Size = static_cast<int>(size);
        d.f4NTables = static_cast<int>(nTab);
        const size_t floatBytes = nTab * size * sizeof(float);
        const size_t floatReals = (floatBytes + sizeof(Real) - 1) / sizeof(Real);
        tableOut.assign((t5.size() + floatReals) * sizeof(Real), 0);
        Real *out = reinterpret_cast<Real *>(tableOut.data());
        for (size_t k = 0; k < t5.size(); k++) out[k] = Real(t5[k]);
        unsigned char *f = tableOut.data() + t5.size() * sizeof(Real);
        for (size_t k = 0; k < nTab; k++) {
            std::memcpy(
                f + k * size * sizeof(float),
                e4->get_table(static_cast<int>(k)).data(),
                size * sizeof(float));
        }
        // the fast division of the Metal kernels (mdDivF, device_defs.h)
        // takes these divisors normal and below 2^126, and the tables
        // finite: the 4D lookup divides by the grid spacings, Tref, T5diff
        // and powers of the clamped T (device_eos.h, eos4DEval)
        {
            auto okDivisor = [](double x) {
                const float f = static_cast<float>(x);
                return (std::isnormal(f) && std::fabs(f) < 8.507059173023462e37f);
            };
            const double pi = 3.14159265358979323846;
            const double TLo = std::max<double>(
                d.T04, std::min<double>(d.Tmax4, d.TrefDilute4));
            const double THi = std::max<double>(d.T04, d.Tmax4);
            bool ok = okDivisor(d.dT4) && okDivisor(d.dmub4)
                      && okDivisor(d.dmuq4) && okDivisor(d.dmus4)
                      && okDivisor(d.TrefDilute4) && TLo > 0.
                      && okDivisor(TLo * TLo)
                      && okDivisor(19. * pi * pi * TLo * TLo * TLo)
                      && okDivisor(19. * pi * pi * THi * THi * THi);
            for (size_t k = 0; ok && k < t5.size(); k++)
                ok = okDivisor(t5[k]);
            const float *tf = reinterpret_cast<const float *>(f);
            for (size_t k = 0; ok && k < nTab * size; k++)
                ok = std::isfinite(tf[k]);
            if (!ok) {
                err = "the 4D EoS grid or tables are outside the range of "
                      "the device division";
                return;
            }
        }
        descOut.resize(sizeof(d));
        std::memcpy(descOut.data(), &d, sizeof(d));
        return;
    }
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
        d.muSOffset = d.muQOffset = d.edgeOffset = -1;
        // the thermal tables also for ideal hydro: the T of the grid maxima
        if (d.type == 1) {
            // EOS_hotQCD::get_temperature: T^5 on the e grid of table 0
            d.tOffset = static_cast<int>(
                table.size() + b.temperature_tb.offset(0));
            const double *t = b.temperature_tb.data();
            table.insert(table.end(), t, t + b.temperature_tb.size());
        } else {
            // EOS_neos: T^5, mu_B (mu_S, mu_Q if tabulated) as the pressure,
            // then the n_B edges
            const auto append = [&](const EosTable &tb) {
                const int off = static_cast<int>(table.size());
                const double *t = tb.data();
                table.insert(table.end(), t, t + tb.size());
                return (off);
            };
            d.t2Offset = append(b.temperature_tb);
            d.muBOffset = append(b.mu_B_tb);
            if (b.get_flag_muS()) d.muSOffset = append(b.mu_S_tb);
            if (b.get_flag_muQ()) d.muQOffset = append(b.mu_Q_tb);
            if (!b.nb_edge.empty()) {
                d.edgeOffset = static_cast<int>(table.size());
                for (int t = 0; t < d.nTables; t++) {
                    if (static_cast<int>(b.nb_edge[t].size()) != d.eLength[t]) {
                        err = "the n_B edges do not match the EoS tables";
                        return;
                    }
                    table.insert(
                        table.end(), b.nb_edge[t].begin(), b.nb_edge[t].end());
                }
            }
        }
    } else {
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
    const bool ok = backend_->setEos(
        eosDesc_.data(), eosDesc_.size(),
        eosTable_.empty() ? nullptr : eosTable_.data(), eosTable_.size(), err);
    // the backend keeps its own copy
    std::vector<unsigned char>().swap(eosTable_);
    return (ok);
}

template <typename Real>
void DeviceHydro::fillParams(
    double tau, bool hasSource, std::vector<unsigned char> &out) const {
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
    P.nComp = (viscous_ ? int(music_device::kNCompVisc)
                        : int(music_device::kNComp));
    if (viscous_) fillViscParams(P);
    P.hasSource = (hasSource ? 1 : 0);
    P.stringSource = (hasSource && strings_ != nullptr ? 1 : 0);
    P.etaMin = Real(DATA.etaMin());
    out.resize(sizeof(P));
    std::memcpy(out.data(), &P, sizeof(P));
}

std::string DeviceHydro::specKey() const {
    std::vector<unsigned char> bytes;
    fillParams<double>(DATA.tau0, false, bytes);
    StepParams<double> P;
    std::memcpy(&P, bytes.data(), sizeof(P));
    const int id = eos_.get_eos_id();
    // as in fillEos
    const int eosType =
        (id == 0 ? 0 : (id == 9 || id == 91 ? 1 : (id == 20 ? 3 : 2)));
    std::ostringstream k;
    k << "mode=" << P.mode << ";storeQ=" << P.storeQ << ";nComp=" << P.nComp
      << ";shear=" << P.shear << ";bulk=" << P.bulk
      << ";secondOrder=" << P.secondOrder << ";vorticity=" << P.vorticity
      << ";duIdeal=" << P.duIdeal << ";resum=" << P.resum
      << ";causality=" << P.causality << ";bulkChem=" << P.bulkChem
      << ";diff=" << P.diff << ";profile1=" << P.profile1
      << ";muBShear=" << P.muBShear << ";bulkRelaxType=" << P.bulkRelaxType
      << ";questRevert=" << P.questRevert << ";nWComp=" << P.nWComp
      << ";eosType=" << eosType;
    return (k.str());
}

template <typename Real>
void DeviceHydro::fillViscParams(StepParams<Real> &P) const {
    P.shear = DATA.turn_on_shear;
    P.bulk = DATA.turn_on_bulk;
    P.secondOrder = DATA.include_second_order_terms;
    P.vorticity = (DATA.include_vorticity_terms ? 1 : 0);
    P.duIdeal = DATA.Du_from_ideal_EoM;
    P.resum = (DATA.FlagResumTransportCoeff ? 1 : 0);
    P.causality = (DATA.FlagCausalityMethod == 1 ? 1 : 0);
    P.bulkChem = DATA.turn_on_bulk_chem;
    P.diff = DATA.turn_on_diff;
    P.profile1 = (DATA.Initial_profile == 1 ? 1 : 0);
    P.kappaCoeff = Real(DATA.kappa_coefficient);
    P.chemRateC = Real(DATA.chem_rate_C);
    P.muBShear = DATA.muB_dependent_shear_to_s;
    P.bulkRelaxType = DATA.bulk_relaxation_type;
    P.questRevert = (DATA.Initial_profile != 0 && DATA.Initial_profile != 1);
    P.shearT = DATA.T_dependent_shear_to_s;
    P.bulkT = DATA.T_dependent_bulk_to_s;
    // Diss::nWComponents (no chemical bulk, no diffusion)
    // Diss::nWComponents
    P.nWComp = 5;
    if (DATA.turn_on_bulk == 1) P.nWComp = 6;
    if (DATA.turn_on_bulk_chem == 1) P.nWComp = 7;
    if (DATA.turn_on_diff == 1) P.nWComp = 10;
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

template <typename Real, typename Keep>
void DeviceHydro::downloadT(
    Fields &arena, const Real *s, const Keep &keep) const {
    const int n = nx_ * ny_ * neta_;
#pragma omp parallel for collapse(2) schedule(dynamic, 16)
    for (int ieta = 0; ieta < neta_; ieta++)
        for (int iy = 0; iy < ny_; iy++) {
            // the x ranges [r[0], r[1]) and [r[2], r[3]) of the row
            int r[4] = {0, 0, 0, 0};
            const int nr = keep(iy, ieta, r);
            for (int k = 0; k < nr; k++)
            for (int ix = r[2 * k]; ix < r[2 * k + 1]; ix++) {
                const int i = ix + nx_ * (iy + ny_ * ieta);
                storeCell(arena, arena.getFieldIdx(ix, iy, ieta), s, n, i);
            }
        }
}

template <typename Real>
void DeviceHydro::storeCell(
    Fields &arena, int idx, const Real *s, size_t stride, size_t at) const {
    const auto v = [s, stride, at](int c) { return (s[c * stride + at]); };
    arena.e_[idx] = v(0);
    arena.rhob_[idx] = v(1);
    arena.rhoq_[idx] = v(2);
    arena.rhos_[idx] = v(3);
    for (int k = 0; k < 4; k++) arena.u_[k][idx] = v(4 + k);
    if (viscous_) {
        for (int k = 0; k < 14; k++)
            arena.Wmunu_[k][idx] = v(music_device::kOffW + k);
        arena.piBulk_[idx] = v(music_device::kOffPi);
        arena.piBulkChem_[idx] = v(music_device::kOffPiChem);
    }
}

template <typename Real>
void DeviceHydro::downloadGathered(
    Fields &arena, const Real *g, const std::vector<int> &cells) const {
    const int nCells = static_cast<int>(cells.size());
#pragma omp parallel for schedule(static)
    for (int k = 0; k < nCells; k++) {
        const int i = cells[k];
        const int ix = i % nx_;
        const int iy = (i / nx_) % ny_;
        const int ieta = i / (nx_ * ny_);
        storeCell(arena, arena.getFieldIdx(ix, iy, ieta), g, nCells, k);
    }
}

template <typename Real>
bool DeviceHydro::computeSource(const Real *s, double tauRk) {
    const int n = nx_ * ny_ * neta_;
    Real *out = static_cast<Real *>(backend_->source());
    const HydroSourceBase &src = *source_;
    bool any = false;
    bool isNan = false;
#pragma omp parallel for collapse(2) reduction(|| : any, isNan)
    for (int ieta = 0; ieta < neta_; ieta++)
        for (int iy = 0; iy < ny_; iy++)
            for (int ix = 0; ix < nx_; ix++) {
                const int i = ix + nx_ * (iy + ny_ * ieta);
                // as in Advance::AdvanceIt and FirstRKStepT
                double eta_s_local = DATA.etaAt(ieta);
                double x_local = -DATA.x_size / 2. + ix * DATA.delta_x;
                double y_local = -DATA.y_size / 2. + iy * DATA.delta_y;
                double qi_source[music_device::kNTJ] = {0.};
                EnergyFlowVec j_mu = {0};
                FlowVec u_local;
                for (int ii = 0; ii < 4; ii++)
                    u_local[ii] = s[(4 + ii) * n + i];
                src.get_hydro_energy_source(
                    tauRk, x_local, y_local, eta_s_local, u_local, j_mu);
                for (int ii = 0; ii < 4; ii++) {
                    qi_source[ii] = tauRk * j_mu[ii];
                    if (std::isnan(qi_source[ii])) isNan = true;
                }
                if (DATA.turn_on_rhob == 1) {
                    qi_source[4] =
                        (tauRk
                         * src.get_hydro_rhob_source(
                             tauRk, x_local, y_local, eta_s_local, u_local));
                }
                if (DATA.turn_on_QS == 1) {
                    qi_source[5] =
                        (tauRk
                         * src.get_hydro_rhoQ_source(
                             tauRk, x_local, y_local, eta_s_local, u_local));
                }
                for (int a = 0; a < music_device::kNTJ; a++) {
                    out[a * n + i] = Real(qi_source[a]);
                    if (qi_source[a] != 0.) any = true;
                }
            }
    // Advance::FirstRKStepT exits here
    if (isNan) fail("qi_source is nan");
    return (any);
}

void DeviceHydro::upload(const Fields &curr, const Fields &prev, double tau) {
    std::string err;
    consAheadValid_ = false;
    if (curr.nX() != nx_ || curr.nY() != ny_ || curr.nEta() != neta_) {
        nx_ = curr.nX();
        ny_ = curr.nY();
        neta_ = curr.nEta();
        const int nComp =
            (viscous_ ? int(music_device::kNCompVisc)
                      : int(music_device::kNComp));
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
    outflowValid_[0] = outflowValid_[1] = false;
    outflowStatic_[0] = outflowStatic_[1] = false;
}

template <typename Keep>
void DeviceHydro::downloadCells(
    Fields &arena, int which, const Keep &keep, bool during) {
    std::string err;
    const int nComp =
        (viscous_ ? int(music_device::kNCompVisc)
                  : int(music_device::kNComp));
    // the listed cells only, where the backend can gather them (CUDA:
    // the diagnostics read about a tenth of the slice). During a step
    // (slice 1 only, as for fetchUnwritten) the backend gathers without
    // waiting for the stages
    cells_.clear();
    for (int ieta = 0; ieta < neta_; ieta++)
        for (int iy = 0; iy < ny_; iy++) {
            int r[4] = {0, 0, 0, 0};
            const int nr = keep(iy, ieta, r);
            for (int k = 0; k < nr; k++)
                for (int ix = r[2 * k]; ix < r[2 * k + 1]; ix++)
                    cells_.push_back(ix + nx_ * (iy + ny_ * ieta));
        }
    if (cells_.size() < static_cast<size_t>(nx_) * ny_ * neta_) {
        const void *g = nullptr;
        if (backend_->fetchCells(which, cells_, g, err)) {
            if (float_) {
                downloadGathered(arena, static_cast<const float *>(g), cells_);
            } else {
                downloadGathered(arena, static_cast<const double *>(g), cells_);
            }
            return;
        }
        if (!err.empty()) fail(err);
    }
    if (during ? !backend_->fetchUnwritten(which, 0, nComp, err)
               : !backend_->fetch(which, 0, nComp, err))
        fail(err);
    if (float_) {
        downloadT(arena, static_cast<const float *>(backend_->state(which)),
                  keep);
    } else {
        downloadT(arena, static_cast<const double *>(backend_->state(which)),
                  keep);
    }
}

void DeviceHydro::syncHost(Fields &curr, Fields &prev) {
    requireIdle("syncHost");
    if (hostCurrent_) return;
    if (&curr != lastCurr_ || &prev != lastPrev_)
        fail("syncHost: not the slices of the last step");
    const auto t0 = std::chrono::steady_clock::now();
    const int nx = nx_;
    const auto all = [nx](int, int, int *r) {
        r[1] = nx;
        return (1);
    };
    downloadCells(curr, 0, all);
    downloadCells(prev, 1, all);
    curr.updateGhosts();
    prev.updateGhosts();
    hostCurrent_ = true;
    nFullSyncs_++;
    downloadSeconds_ += std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t0)
                            .count();
}

int DeviceHydro::currentSlice(
    const Fields &curr, const char *what, bool &during) const {
    // during a step the slice of tau_n is the backend's previous one, which
    // the stages read and do not write (fetchUnwritten)
    during = pending_;
    if (&curr != (pending_ ? pendingCurr_ : lastCurr_))
        fail(std::string(what) + ": not the current slice");
    return (pending_ ? 1 : 0);
}

void DeviceHydro::syncEtaWindow(Fields &curr, double etaMin, double etaMax) {
    if (hostCurrent_) return;
    bool during = false;
    const int which = currentSlice(curr, "syncEtaWindow", during);
    const auto t0 = std::chrono::steady_clock::now();
    // the eta of the Cell_info diagnostics
    const bool boostInvariant = DATA.boost_invariant;
    const int nx = nx_;
    downloadCells(curr, which, [&](int, int ieta, int *r) {
        double eta = 0.0;
        if (!boostInvariant) {
            eta = DATA.etaAt(ieta);
        }
        r[1] = nx;
        return (eta < etaMax && eta > etaMin ? 1 : 0);
    }, during);
    nPartialSyncs_++;
    downloadSeconds_ += std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t0)
                            .count();
}

void DeviceHydro::syncShell(
    Fields &curr, Fields *prev, int bandX, int bandY, int bandEta) {
    requireIdle("syncShell");
    if (hostCurrent_) return;
    if (&curr != lastCurr_ || (prev != nullptr && prev != lastPrev_))
        fail("syncShell: not the slices of the last step");
    const auto t0 = std::chrono::steady_clock::now();
    const int nx = nx_, ny = ny_, neta = neta_;
    const auto keep = [&](int iy, int ieta, int *r) {
        if (iy < bandY || iy >= ny - bandY
            || (bandEta > 0 && (ieta < bandEta || ieta >= neta - bandEta))) {
            r[1] = nx;  // the whole row
            return (1);
        }
        r[1] = std::min(bandX, nx);
        r[2] = std::max(nx - bandX, r[1]);
        r[3] = nx;
        return (2);
    };
    downloadCells(curr, 0, keep);
    if (prev != nullptr) downloadCells(*prev, 1, keep);
    nPartialSyncs_++;
    downloadSeconds_ += std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t0)
                            .count();
}

//! the cells of the outflow faces from the device slice s (see
//! DeviceHydro::outflowFaces). In the ideal update the Fields' Pi and W
//! stay as they are: hostStatic says out already holds them for arena
//! the accessors of gatherOutflow: component c of cell i (the cell of face
//! j) in a whole slice, or of face j in the cells gathered for the faces
template <typename Real>
struct SliceCells {
    const Real *s;
    size_t n;
    double operator()(int c, int, int i) const { return (s[c * n + i]); }
};
template <typename Real>
struct FaceCells {
    const Real *g;
    size_t nFaces;
    double operator()(int c, int j, int) const { return (g[c * nFaces + j]); }
};

template <typename Get>
static void gatherOutflow(
    const Get &s, const Fields &arena, bool viscous, bool hostStatic,
    int nx, int ny, int neta, std::vector<double> &out) {
    const int nFaces = Cell_info::outflowFaceCount(nx, ny, neta);
    const int nv = Cell_info::kOutflowVals;
    out.resize(static_cast<size_t>(nv) * nFaces);
    // scattered reads: dynamic, for the slower efficiency cores
#pragma omp parallel for schedule(dynamic, 512)
    for (int j = 0; j < nFaces; j++) {
        int ix, iy, ieta, mu, side;
        Cell_info::outflowFace(j, nx, ny, neta, ix, iy, ieta, mu, side);
        const int i = ix + nx * (iy + ny * ieta);
        double *v = out.data() + static_cast<size_t>(nv) * j;
        // as downloadT fills the Fields; s(c, j, i): component c of cell i,
        // the cell of face j
        for (int k = 0; k < 4; k++) v[k] = s(k, j, i);
        for (int k = 0; k < 4; k++) v[6 + k] = s(4 + k, j, i);
        if (viscous) {
            v[4] = s(music_device::kOffPi, j, i);
            v[5] = s(music_device::kOffPiChem, j, i);
            for (int k = 0; k < 14; k++)
                v[10 + k] = s(music_device::kOffW + k, j, i);
        } else if (!hostStatic) {
            const int idx = arena.getFieldIdx(ix, iy, ieta);
            v[4] = arena.piBulk_[idx];
            v[5] = arena.piBulkChem_[idx];
            for (int k = 0; k < 14; k++) v[10 + k] = arena.Wmunu_[k][idx];
        }
    }
}

const std::vector<double> &DeviceHydro::outflowFaces(
    const Fields &arena, int which) {
    if (which != 0 && which != 1) fail("outflowFaces: which is 0 or 1");
    // during a step: the current slice only (its backend slice is the
    // previous one, which no stage writes); the previous slice's faces
    // only if they were taken before
    bool during = false;
    int slice = which;
    if (which == 0) {
        slice = currentSlice(arena, "outflowFaces", during);
    } else {
        if (&arena != (pending_ ? pendingPrev_ : lastPrev_))
            fail("outflowFaces: not the previous slice");
        if (pending_ && !outflowValid_[1])
            fail("outflowFaces: the previous slice is being written by the "
                 "step");
    }
    if (outflowValid_[which]) return (outflow_[which]);
    const auto t0 = std::chrono::steady_clock::now();
    std::string err;
    const int nComp = (viscous_ ? int(music_device::kNCompVisc)
                                : int(music_device::kNComp));
    const bool hostStatic = outflowStatic_[which];
    // the cells of the faces only, where the backend can gather them (CUDA:
    // a few percent of the slice instead of all of it over PCIe)
    const int nFaces = Cell_info::outflowFaceCount(nx_, ny_, neta_);
    if (static_cast<int>(outflowCells_.size()) != nFaces
        || outflowCellsNx_ != nx_ || outflowCellsNy_ != ny_
        || outflowCellsNeta_ != neta_) {
        outflowCells_.resize(nFaces);
        for (int j = 0; j < nFaces; j++) {
            int ix, iy, ieta, mu, side;
            Cell_info::outflowFace(j, nx_, ny_, neta_, ix, iy, ieta, mu, side);
            outflowCells_[j] = ix + nx_ * (iy + ny_ * ieta);
        }
        outflowCellsNx_ = nx_;
        outflowCellsNy_ = ny_;
        outflowCellsNeta_ = neta_;
    }
    const void *g = nullptr;
    if (backend_->fetchCells(slice, outflowCells_, g, err)) {
        const size_t nf = static_cast<size_t>(nFaces);
        if (float_) {
            gatherOutflow(FaceCells<float>{static_cast<const float *>(g), nf},
                          arena, viscous_, hostStatic, nx_, ny_, neta_,
                          outflow_[which]);
        } else {
            gatherOutflow(
                FaceCells<double>{static_cast<const double *>(g), nf}, arena,
                viscous_, hostStatic, nx_, ny_, neta_, outflow_[which]);
        }
    } else {
        if (!err.empty()) fail(err);
        if (during ? !backend_->fetchUnwritten(slice, 0, nComp, err)
                   : !backend_->fetch(slice, 0, nComp, err))
            fail(err);
        const size_t n = static_cast<size_t>(nx_) * ny_ * neta_;
        if (float_) {
            gatherOutflow(
                SliceCells<float>{
                    static_cast<const float *>(backend_->state(slice)), n},
                arena, viscous_, hostStatic, nx_, ny_, neta_,
                outflow_[which]);
        } else {
            gatherOutflow(
                SliceCells<double>{
                    static_cast<const double *>(backend_->state(slice)), n},
                arena, viscous_, hostStatic, nx_, ny_, neta_,
                outflow_[which]);
        }
    }
    outflowValid_[which] = true;
    outflowStatic_[which] = true;
    nPartialSyncs_++;
    downloadSeconds_ += std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t0)
                            .count();
    return (outflow_[which]);
}

void DeviceHydro::conservationRows(double tau, std::vector<double> &rows) {
    requireIdle("conservationRows");
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<unsigned char> params;
    const size_t nRows = static_cast<size_t>(nx_) * neta_;
    const size_t nVals = music_device::kConsVals * nRows;
    rows.resize(nVals);
    std::string err;
    // the rows the last step computed, if it computed them for this tau
    // (the same kernel, the same parameters)
    if (consAheadValid_ && tau == consAheadTau_) {
        consAheadValid_ = false;
        if (float_) {
            std::vector<float> out(nVals);
            if (!backend_->conservationEnd(
                    out.data(), out.size() * sizeof(float), err))
                fail(err);
            for (size_t k = 0; k < nVals; k++) rows[k] = out[k];
        } else {
            if (!backend_->conservationEnd(
                    rows.data(), rows.size() * sizeof(double), err))
                fail(err);
        }
        nConsAhead_++;
        diagSeconds_ += std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t0)
                            .count();
        return;
    }
    consAheadValid_ = false;
    if (float_) {
        fillParams<float>(tau, false, params);
        std::vector<float> out(nVals);
        if (!backend_->conservation(
                params.data(), params.size(), out.data(),
                out.size() * sizeof(float), err))
            fail(err);
        for (size_t k = 0; k < nVals; k++) rows[k] = out[k];
    } else {
        fillParams<double>(tau, false, params);
        if (!backend_->conservation(
                params.data(), params.size(), rows.data(),
                rows.size() * sizeof(double), err))
            fail(err);
    }
    diagSeconds_ += std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - t0)
                        .count();
}

template <typename Real>
void DeviceHydro::packStrings(
    std::vector<unsigned char> &out, int &nActive) const {
    using namespace music_device;
    const HydroSourceStrings &S = *strings_;
    const auto &all = S.get_strings();
    const HydroSourceStrings::StringList *lists[4] = {
        &S.get_energy_list(), &S.get_remnant_list(), &S.get_baryon_list(),
        &S.get_electric_list()};
    // the union of the four lists in the order of the full list; each
    // list is a subsequence of it
    size_t cursor[4] = {0, 0, 0, 0};
    std::vector<double> buf(kStrHead, 0.);
    int nStr = 0;
    for (const auto &it : all) {
        bool in[4];
        bool any = false;
        for (int l = 0; l < 4; l++) {
            in[l] = (cursor[l] < lists[l]->size()
                     && (*lists[l])[cursor[l]].get() == it.get());
            if (in[l]) cursor[l]++;
            any = any || in[l];
        }
        if (!any) continue;
        const QCD_string &q = *it;
        double r[kStrRec];
        for (int l = 0; l < 4; l++) r[kRInEnergy + l] = (in[l] ? 1. : 0.);
        r[kRNorm] = q.norm;
        r[kRRapidityShift] = q.rapidity_shift;
        r[kRTilt] = q.string_tilt;
        r[kRERemL] = q.E_remnant_norm_L;
        r[kRERemR] = q.E_remnant_norm_R;
        r[kRTauForm] = q.tau_form;
        r[kRSigmaX] = q.sigma_x;
        r[kRSigmaEta] = q.sigma_eta;
        r[kRTau0] = q.tau_0;
        r[kREta0] = q.eta_s_0;
        r[kRXl] = q.x_pl;
        r[kRYl] = q.y_pl;
        r[kRXr] = q.x_pr;
        r[kRYr] = q.y_pr;
        r[kRTauEndL] = q.tau_end_left;
        r[kRTauEndR] = q.tau_end_right;
        r[kREtaL] = q.eta_s_left;
        r[kREtaR] = q.eta_s_right;
        r[kRYL] = q.y_l;
        r[kRYR] = q.y_r;
        r[kRRemnantYL] = q.y_l + q.remnant_shift_L;
        r[kRRemnantYR] = q.y_r + q.remnant_shift_R;
        r[kRTauBL] = q.tau_baryon_left;
        r[kRTauBR] = q.tau_baryon_right;
        r[kREtaBL] = q.eta_s_baryon_left;
        r[kREtaBR] = q.eta_s_baryon_right;
        r[kRYBL] = q.y_l_baryon;
        r[kRYBR] = q.y_r_baryon;
        r[kRBaryonNormL] = q.B_norm_L;
        r[kRBaryonNormR] = q.B_norm_R;
        r[kRTauQL] = q.tau_Qe_left;
        r[kRTauQR] = q.tau_Qe_right;
        r[kRElectricNormL] = q.Q_norm_L;
        r[kRElectricNormR] = q.Q_norm_R;
        r[kREtaQL] = q.eta_s_Qe_left;
        r[kREtaQR] = q.eta_s_Qe_right;
        r[kRSkipX] = HydroSourceStrings::smearingCutoff * q.sigma_x;
        r[kRSkipEta] = HydroSourceStrings::smearingCutoff * q.sigma_eta;
        buf.insert(buf.end(), r, r + kStrRec);
        nStr++;
    }
    buf[kStrNStrings] = nStr;
    for (int l = 0; l < 4; l++) {
        if (cursor[l] != lists[l]->size())
            fail("a string list is not in the order of the full list");
        buf[kStrNEnergy + l] = static_cast<double>(lists[l]->size());
    }
    buf[kStrPreEqFlow] = S.get_preEqFlowFactor();
    buf[kStrQuench] = S.get_parton_quench_factor();
    buf[kStrUseBQ] = (DATA.use_BQ_ratios ? 1. : 0.);
    buf[kStrRhob] = (DATA.turn_on_rhob == 1 ? 1. : 0.);
    buf[kStrQS] = (DATA.turn_on_QS == 1 ? 1. : 0.);
    buf[kStrDtau] = DATA.delta_tau;
    buf[kStrXSize] = DATA.x_size;
    buf[kStrYSize] = DATA.y_size;
    buf[kStrEtaMin] = DATA.etaMin();
    buf[kStrDx] = DATA.delta_x;
    buf[kStrDy] = DATA.delta_y;
    buf[kStrDeta] = DATA.delta_eta;
    buf[kStrHbarc] = Util::hbarc;
    out.resize(buf.size() * sizeof(Real));
    Real *o = reinterpret_cast<Real *>(out.data());
    for (size_t k = 0; k < buf.size(); k++) o[k] = Real(buf[k]);
    nActive = 0;
    for (int l = 0; l < 4; l++) nActive += static_cast<int>(lists[l]->size());
}


void DeviceHydro::syncCurrentDuringStep(Fields &curr) {
    if (!pending_) fail("syncCurrentDuringStep: no step is running");
    if (&curr != pendingCurr_)
        fail("syncCurrentDuringStep: not the slice the step started from");
    if (hostCurrent_) return;
    const auto t0 = std::chrono::steady_clock::now();
    const int nx = nx_;
    const auto all = [nx](int, int, int *r) {
        r[1] = nx;
        return (1);
    };
    // after beginStep the slice of tau_n is the backend's previous one
    std::string err;
    const int nComp =
        (viscous_ ? int(music_device::kNCompVisc)
                  : int(music_device::kNComp));
    if (!backend_->fetchUnwritten(1, 0, nComp, err)) fail(err);
    if (float_) {
        downloadT(curr, static_cast<const float *>(backend_->state(1)), all);
    } else {
        downloadT(curr, static_cast<const double *>(backend_->state(1)), all);
    }
    curr.updateGhosts();
    nOverlapSyncs_++;
    downloadSeconds_ += std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t0)
                            .count();
}

void DeviceHydro::requireIdle(const char *what) const {
    if (pending_)
        fail(std::string(what)
             + ": the device is running a step (beginStep without "
               "advanceStep)");
}

void DeviceHydro::advanceStep(
    double tau, Fields *&fpPrev, Fields *&fpCurr, Fields *&fpNext) {
    if (!pending_) {
        beginStep(tau, fpPrev, fpCurr);
    } else if (std::abs(tau - pendingTau_) > 1e-12 || fpCurr != pendingCurr_
               || fpPrev != pendingPrev_) {
        fail("advanceStep: not the step beginStep started");
    }
    const auto t0 = std::chrono::steady_clock::now();
    int flags = pendingFlags_, f = 0;
    std::string err;
    if (!backend_->finish(f, err)) fail(err);
    flags |= f;
    pending_ = false;
    stepSeconds_ += std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - t0)
                        .count();
    if (flags & music_device::kFlagSourceNan) fail("qi_source is nan");
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
    // Evolve::AdvanceRK: after both stages the slice that was the previous
    // one holds the new state, the current one becomes the previous one.
    // The host Fields are filled when Evolve asks (syncHost, ...)
    Fields *oldPrev = fpPrev;
    fpPrev = fpCurr;
    fpCurr = oldPrev;
    (void)fpNext;
    lastCurr_ = fpCurr;
    lastPrev_ = fpPrev;
    valid_ = true;
    hostCurrent_ = false;
    // the previous slice is the device buffer that was the current one
    // (and fpPrev the Fields that were fpCurr: both rotate together)
    outflow_[1].swap(outflow_[0]);
    std::swap(outflowStatic_[0], outflowStatic_[1]);
    outflowValid_[1] = outflowValid_[0];
    outflowValid_[0] = false;
    nSteps_++;
}

void DeviceHydro::beginStep(double tau, Fields *fpPrev, Fields *fpCurr) {
    requireIdle("beginStep");
    if (!valid_ || fpCurr != lastCurr_ || (viscous_ && fpPrev != lastPrev_)
        || fpCurr->nX() != nx_ || fpCurr->nY() != ny_
        || fpCurr->nEta() != neta_) {
        if (viscous_ && (fpPrev->nX() != fpCurr->nX()
                         || fpPrev->nY() != fpCurr->nY()
                         || fpPrev->nEta() != fpCurr->nEta())) {
            fail("the previous and current slices differ in size");
        }
        // the host Fields must hold the last step: Evolve fills them
        // (syncHost) before it changes them and calls invalidate()
        if (!hostCurrent_)
            fail("the host Fields are older than the device state");
        upload(*fpCurr, *fpPrev, tau);
    }
    int flags = 0;
    std::string err;
    double sourceSeconds = 0.;
    const auto t0 = std::chrono::steady_clock::now();
    // the strings source: the strings of this tau frame, once for both
    // stages (Evolve called prepare_list_for_current_tau_frame)
    bool stringSource = false;
    if (strings_) {
        const auto s0 = std::chrono::steady_clock::now();
        std::vector<unsigned char> packed;
        int nActive = 0;
        if (float_) {
            packStrings<float>(packed, nActive);
        } else {
            packStrings<double>(packed, nActive);
        }
        stringSource = (nActive > 0);
        if (stringSource
            && !backend_->setStrings(packed.data(), packed.size(), err))
            fail(err);
        sourceSeconds += std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - s0)
                             .count();
    }
    for (int rk = 0; rk < 2; rk++) {
        // the sources at tau_rk from the u of the stage's slice: tau_n, then
        // the predictor
        bool hasSource = stringSource;
        if (source_ && !strings_) {
            const auto s0 = std::chrono::steady_clock::now();
            const int which = (rk == 0 ? 0 : 2);
            // u of the slice, in both stages: with a copy of its own
            // (CUDA) state(0) is only as recent as the last fetch, which
            // the diagnostics do not make every step (on Metal and the CPU
            // state() is the slice and fetch() copies nothing)
            if (!backend_->fetch(which, 4, 4, err)) fail(err);
            const double tauRk = tau + rk * DATA.delta_tau;
            if (float_) {
                hasSource = computeSource(
                    static_cast<const float *>(backend_->state(which)), tauRk);
            } else {
                hasSource = computeSource(
                    static_cast<const double *>(backend_->state(which)),
                    tauRk);
            }
            if (hasSource && !backend_->uploadSource(err)) fail(err);
            sourceSeconds += std::chrono::duration<double>(
                                 std::chrono::steady_clock::now() - s0)
                                 .count();
        }
        std::vector<unsigned char> params;
        if (float_) {
            fillParams<float>(tau, hasSource, params);
        } else {
            fillParams<double>(tau, hasSource, params);
        }
        int f = 0;
        if (!backend_->stage(params.data(), params.size(), rk, f, err))
            fail(err);
        flags |= f;
    }
    const auto t1 = std::chrono::steady_clock::now();
    // the conservation rows of the new slice, for conservationRows of the
    // next iteration (Evolve advances tau by DATA.delta_tau)
    consAheadValid_ = false;
    if (consAhead_) {
        const double tauNext = tau + DATA.delta_tau;
        std::vector<unsigned char> params;
        const size_t nVals =
            music_device::kConsVals * static_cast<size_t>(nx_) * neta_;
        size_t outBytes = 0;
        if (float_) {
            fillParams<float>(tauNext, false, params);
            outBytes = nVals * sizeof(float);
        } else {
            fillParams<double>(tauNext, false, params);
            outBytes = nVals * sizeof(double);
        }
        if (!backend_->conservationBegin(
                params.data(), params.size(), outBytes, err))
            fail(err);
        consAheadValid_ = true;
        consAheadTau_ = tauNext;
    }
    const auto t2 = std::chrono::steady_clock::now();
    diagSeconds_ += std::chrono::duration<double>(t2 - t1).count();
    // the device state of tau_n stays as it is; advanceStep checks the
    // flags and rotates the slices
    valid_ = true;
    pending_ = true;
    pendingTau_ = tau;
    pendingCurr_ = fpCurr;
    pendingPrev_ = fpPrev;
    pendingFlags_ = flags;
    stepSeconds_ +=
        std::chrono::duration<double>(t1 - t0).count() - sourceSeconds;
    sourceSeconds_ += sourceSeconds;
}

void DeviceHydro::fail(const std::string &what) const {
    music_message << "gpu_backend: " << what;
    music_message.flush("error");
    exit(1);
}
