// Copyright 2019 Chun Shen

#include "hydro_source_strings.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include "device/device_strings.h"
#include "util.h"

using music_device::stringBinEnd;
using music_device::stringBinStart;
using std::string;

constexpr double HydroSourceStrings::smearingCutoff;

HydroSourceStrings::HydroSourceStrings(InitData &DATA_in) : DATA(DATA_in) {
    music_message.set_tag("🧵 HydroSourceStrings");
    set_source_tau_min(std::numeric_limits<double>::infinity());
    set_source_tau_max(0.0);
    set_source_tauStart_max(0.0);
    set_sigma_tau(0.1);
    set_sigma_x(DATA.stringSourceSigmaX);
    set_sigma_eta(DATA.stringSourceSigmaEta);
    string_dump_mode = DATA.string_dump_mode;
    string_quench_factor = DATA.string_quench_factor;
    parton_quench_factor = 1.0;  // no diffusion current from the source
    preEqFlowFactor_ = DATA.preEqFlowFactor;
    read_in_QCD_strings_and_partons();
}

HydroSourceStrings::~HydroSourceStrings() { QCD_strings_list.clear(); }

double HydroSourceStrings::getStringEndTau(
    const double tau_0, const double tau_form, const double eta_s_0,
    const double eta_s) const {
    const double delta = eta_s - eta_s_0;
    if (tau_0 == 0.) return tau_form;
    return tau_0 * std::cosh(delta)
           + std::hypot(tau_form, tau_0 * std::sinh(delta));
}

//! This function reads in the spatal information of the strings and partons
//! which are produced from the MC-Glauber-LEXUS model
void HydroSourceStrings::read_in_QCD_strings_and_partons() {
    string QCD_strings_filename = DATA.initName;
    music_message << "read in QCD strings list from " << QCD_strings_filename;
    music_message.flush("info");
    string text_string;

    std::ifstream QCD_strings_file(QCD_strings_filename.c_str());
    if (!QCD_strings_file) {
        music_message << "hydro_source::read_in_QCD_strings_and_partons: "
                      << "can not open QCD strings file: "
                      << QCD_strings_filename;
        music_message.flush("error");
        exit(1);
    }

    // Collision metadata is optional in older string files.
    getline(QCD_strings_file, text_string);
    std::replace(text_string.begin(), text_string.end(), '=', ' ');
    std::stringstream coll_info(text_string);
    std::string key;
    double total_energy = std::numeric_limits<double>::quiet_NaN();
    double netPz = std::numeric_limits<double>::quiet_NaN();
    double nPart = std::numeric_limits<double>::quiet_NaN();
    while (coll_info >> key) {
        if (key != "total_energy" && key != "net_Pz" && key != "Npart")
            continue;
        std::string value;
        double parsed;
        if (coll_info >> value && std::istringstream(value) >> parsed) {
            if (key == "total_energy")
                total_energy = parsed;
            else if (key == "net_Pz")
                netPz = parsed;
            else
                nPart = parsed;
        }
    }

    getline(QCD_strings_file, text_string);  // column names
    int line = 2;
    while (getline(QCD_strings_file, text_string)) {
        ++line;
        std::stringstream text_stream(text_string);
        if ((text_stream >> std::ws).eof()) continue;
        auto new_string = std::make_shared<QCD_string>();
        if (!(text_stream >> new_string->mass >> new_string->m_over_sigma
              >> new_string->tau_form >> new_string->tau_0
              >> new_string->eta_s_0 >> new_string->x_perp >> new_string->y_perp
              >> new_string->x_pl >> new_string->y_pl >> new_string->x_pr
              >> new_string->y_pr >> new_string->eta_s_left
              >> new_string->eta_s_right >> new_string->y_l >> new_string->y_r
              >> new_string->remnant_l >> new_string->remnant_r
              >> new_string->y_l_i >> new_string->y_r_i
              >> new_string->eta_s_baryon_left >> new_string->eta_s_baryon_right
              >> new_string->y_l_baryon >> new_string->y_r_baryon
              >> new_string->baryon_frac_l >> new_string->baryon_frac_r
              >> new_string->Qe_left >> new_string->Qe_right
              >> new_string->eta_s_Qe_left >> new_string->eta_s_Qe_right)
            || !(text_stream >> std::ws).eof()
            || new_string->eta_s_left > new_string->eta_s_right
            || new_string->mass < 0. || new_string->tau_form < 0.
            || new_string->tau_0 < 0. || new_string->remnant_l < 0.
            || new_string->remnant_r < 0.) {
            throw std::runtime_error(
                "Invalid string in " + QCD_strings_filename + " at line "
                + std::to_string(line));
        }

        if (DATA.Initial_profile == 131) {
            new_string->tau_0 = 0.;
            new_string->eta_s_0 = 0.;
        }
        new_string->sigma_x = get_sigma_x();
        new_string->sigma_eta = get_sigma_eta();

        // compute the string end tau
        new_string->tau_end_left = getStringEndTau(
            new_string->tau_0, new_string->tau_form, new_string->eta_s_0,
            new_string->eta_s_left);
        new_string->tau_end_right = getStringEndTau(
            new_string->tau_0, new_string->tau_form, new_string->eta_s_0,
            new_string->eta_s_right);

        // compute the baryon number tau
        new_string->tau_baryon_left = getStringEndTau(
            new_string->tau_0, new_string->tau_form, new_string->eta_s_0,
            new_string->eta_s_baryon_left);
        new_string->tau_baryon_right = getStringEndTau(
            new_string->tau_0, new_string->tau_form, new_string->eta_s_0,
            new_string->eta_s_baryon_right);

        // compute the electric charge tau
        new_string->tau_Qe_left = getStringEndTau(
            new_string->tau_0, new_string->tau_form, new_string->eta_s_0,
            new_string->eta_s_Qe_left);
        new_string->tau_Qe_right = getStringEndTau(
            new_string->tau_0, new_string->tau_form, new_string->eta_s_0,
            new_string->eta_s_Qe_right);

        // determine the tau_start and eta_s_start of the string
        if (new_string->eta_s_left > new_string->eta_s_0) {
            new_string->tau_start = new_string->tau_end_left;
            new_string->eta_s_start = new_string->eta_s_left;
        } else if (new_string->eta_s_right < new_string->eta_s_0) {
            new_string->tau_start = new_string->tau_end_right;
            new_string->eta_s_start = new_string->eta_s_right;
        } else {
            new_string->tau_start = new_string->tau_0 + new_string->tau_form;
            new_string->eta_s_start = new_string->eta_s_0;
        }

        // Sources later than string_source_tau_max are deposited at it.
        const double tauCap = DATA.stringSourceTauMax;
        for (double *t :
             {&new_string->tau_end_left, &new_string->tau_end_right,
              &new_string->tau_baryon_left, &new_string->tau_baryon_right,
              &new_string->tau_Qe_left, &new_string->tau_Qe_right,
              &new_string->tau_start})
            *t = std::min(*t, tauCap);

        // read in one string properly
        QCD_strings_list.push_back(new_string);

        // record the proper time of the first and last string sources
        double source_tau =
            std::max(new_string->tau_end_left, new_string->tau_end_right);
        if (new_string->baryon_frac_l != 0.)
            source_tau = std::max(source_tau, new_string->tau_baryon_left);
        if (new_string->baryon_frac_r != 0.)
            source_tau = std::max(source_tau, new_string->tau_baryon_right);
        if (!DATA.use_rhoQS_to_rhoB_ratios && new_string->Qe_left != 0.)
            source_tau = std::max(source_tau, new_string->tau_Qe_left);
        if (!DATA.use_rhoQS_to_rhoB_ratios && new_string->Qe_right != 0.)
            source_tau = std::max(source_tau, new_string->tau_Qe_right);

        if (new_string->eta_s_left < 1 && new_string->eta_s_right > -1) {
            // the string in the mid-rapidity
            if (new_string->tau_start > get_source_tauStart_max()) {
                set_source_tauStart_max(new_string->tau_start);
            }
        }

        if (!std::isfinite(source_tau))
            throw std::runtime_error(
                "Non-finite formation time in " + QCD_strings_filename);
        set_source_tau_max(std::max(get_source_tau_max(), source_tau));
        double source_start = new_string->tau_start;
        if (new_string->baryon_frac_l != 0.)
            source_start = std::min(source_start, new_string->tau_baryon_left);
        if (new_string->baryon_frac_r != 0.)
            source_start = std::min(source_start, new_string->tau_baryon_right);
        if (!DATA.use_rhoQS_to_rhoB_ratios && new_string->Qe_left != 0.)
            source_start = std::min(source_start, new_string->tau_Qe_left);
        if (!DATA.use_rhoQS_to_rhoB_ratios && new_string->Qe_right != 0.)
            source_start = std::min(source_start, new_string->tau_Qe_right);
        set_source_tau_min(std::min(get_source_tau_min(), source_start));
    }
    if (QCD_strings_list.empty())
        throw std::runtime_error("No strings in " + QCD_strings_filename);

    // The beam rapidity of this event: the energy per participant from the
    // header (Npart, or the baryon number of the strings when Npart is
    // missing), otherwise the one from ecm.
    double beamRapidity = DATA.beam_rapidity;
    double participants = nPart;
    if (!(participants > 0.)) {
        participants = 0.;
        for (auto const &it : QCD_strings_list)
            participants += it->baryon_frac_l + it->baryon_frac_r;
    }
    if (std::isfinite(total_energy) && participants > 0.
        && total_energy > participants * Util::m_N) {
        beamRapidity = acosh(total_energy / (participants * Util::m_N));
        music_message << "beam rapidity of the event = "
                      << beamRapidity << " (from ecm: " << DATA.beam_rapidity
                      << ")";
        music_message.flush("info");
    }
    double minSigmaEta = std::numeric_limits<double>::infinity();
    for (auto &it : QCD_strings_list) {
        const double outerEta =
            std::max(std::abs(it->eta_s_left), std::abs(it->eta_s_right));
        if (outerEta + 3. * it->sigma_eta >= beamRapidity) {
            // Keep the three-sigma edge strictly inside the beam unless the
            // minimum width already extends past it.
            const double insideBeam = std::nextafter(beamRapidity, 0.);
            const double allowedWidth =
                std::nextafter((insideBeam - outerEta) / 3., 0.);
            it->sigma_eta =
                std::max(0.1, std::min(it->sigma_eta, allowedWidth));
        }
        minSigmaEta = std::min(minSigmaEta, it->sigma_eta);
    }
    music_message << "tau_min = " << get_source_tau_min()
                  << " fm/c.";
    music_message.flush("info");
    music_message << "tau_max = " << get_source_tau_max()
                  << " fm/c.";
    music_message.flush("info");

    double total_baryon_number = 0;
    double total_electric_charges = 0;
    for (auto const &it : QCD_strings_list) {
        total_baryon_number += it->baryon_frac_l + it->baryon_frac_r;
        total_electric_charges += it->Qe_left + it->Qe_right;
    }
    music_message << "total baryon number = " << total_baryon_number;
    music_message << ", total electric charges = "
                  << (DATA.use_rhoQS_to_rhoB_ratios
                          ? DATA.ratio_q * total_baryon_number
                          : total_electric_charges);
    music_message.flush("info");
    if (DATA.dynamicGridDetermination) {
        double xMax = 0.;
        double yMax = 0.;
        double etaMin = std::numeric_limits<double>::infinity();
        double etaMax = -etaMin;
        for (auto const &it : QCD_strings_list) {
            xMax = std::max(xMax, std::abs(it->x_perp));
            xMax = std::max(xMax, std::abs(it->x_pl));
            xMax = std::max(xMax, std::abs(it->x_pr));
            yMax = std::max(yMax, std::abs(it->y_perp));
            yMax = std::max(yMax, std::abs(it->y_pl));
            yMax = std::max(yMax, std::abs(it->y_pr));
            const auto includeEta = [&](double eta) {
                const double padding = 5. * it->sigma_eta;
                etaMin = std::min(etaMin, eta - padding);
                etaMax = std::max(etaMax, eta + padding);
            };
            includeEta(it->eta_s_left);
            includeEta(it->eta_s_right);
            if (it->baryon_frac_l != 0.) includeEta(it->eta_s_baryon_left);
            if (it->baryon_frac_r != 0.) includeEta(it->eta_s_baryon_right);
            if (!DATA.use_rhoQS_to_rhoB_ratios && it->Qe_left != 0.)
                includeEta(it->eta_s_Qe_left);
            if (!DATA.use_rhoQS_to_rhoB_ratios && it->Qe_right != 0.)
                includeEta(it->eta_s_Qe_right);

            // Bound the tilted centerline as well as the fixed end sources.
            const double xc = 0.5 * (it->x_pl + it->x_pr);
            const double yc = 0.5 * (it->y_pl + it->y_pr);
            if (DATA.string_tilt_update) {
                if (std::abs(it->mass - MASS_PARTON) < MASS_PARTON_TOLERANCE) {
                    if (it->x_pl != it->x_pr)
                        xMax = std::max(xMax, std::abs(xc) + 2.);
                    if (it->y_pl != it->y_pr)
                        yMax = std::max(yMax, std::abs(yc) + 2.);
                }
            } else {
                for (double frac : {0., 1.}) {
                    xMax = std::max(
                        xMax, std::abs(getStringTransverseCoord(
                                  it->x_pl, it->x_pr, frac,
                                  DATA.stringTransverseShiftFrac)));
                    yMax = std::max(
                        yMax, std::abs(getStringTransverseCoord(
                                  it->y_pl, it->y_pr, frac,
                                  DATA.stringTransverseShiftFrac)));
                }
            }
        }

        // adjust transverse grid size
        double energyAddRadius = sqrt(
            (std::isfinite(total_energy) && total_energy > 0. ? total_energy
                                                              : 0.)
            / std::max(1., total_baryon_number) / 2760.);
        double npartAddRadius = total_baryon_number / 500.;
        double reRunAddRadius = DATA.reRunCount * 2.;
        double gridOffset =
            (DATA.gridPadding + energyAddRadius + npartAddRadius
             + reRunAddRadius);
        gridOffset = std::max(
            gridOffset,
            (5. + std::abs(preEqFlowFactor_)) * DATA.stringSourceSigmaX);
        DATA.x_size = 2. * (xMax + gridOffset);
        DATA.y_size = 2. * (yMax + gridOffset);
        const auto resolve = [](int requested, int inputCount, double size,
                                double maxSpacing) {
            double count =
                std::max<double>(requested, std::ceil(size / maxSpacing) + 1.);
            if (inputCount > 0) count = std::min(count, 2. * inputCount);
            if (!std::isfinite(size) || !std::isfinite(count)
                || count > std::numeric_limits<int>::max())
                throw std::runtime_error(
                    "String-source grid exceeds the supported number of cells");
            return static_cast<int>(count);
        };
        DATA.nx = resolve(
            DATA.nx, DATA.inputGridPoints[0], DATA.x_size, get_sigma_x());
        DATA.ny = resolve(
            DATA.ny, DATA.inputGridPoints[1], DATA.y_size, get_sigma_x());
        DATA.delta_x = DATA.x_size / (DATA.nx - 1);
        DATA.delta_y = DATA.y_size / (DATA.ny - 1);

        if (!DATA.boost_invariant) {
            DATA.eta_min = etaMin;
            DATA.eta_size = etaMax - etaMin;
            DATA.neta = resolve(
                DATA.neta, DATA.inputGridPoints[2], DATA.eta_size,
                minSigmaEta / 2.);
            DATA.delta_eta = DATA.eta_size / (DATA.neta - 1);
        }

        if (DATA.resetDtau) {
            // make sure delta_tau is not too large for delta_x and delta_y
            DATA.delta_tau = std::min(
                DATA.delta_tau, std::min(
                                    DATA.delta_x * DATA.dtaudxRatio,
                                    DATA.delta_y * DATA.dtaudxRatio));
            if (DATA.delta_tau > 0.001)
                DATA.delta_tau =
                    (static_cast<int>(DATA.delta_tau * 1000)) / 1000.;
            DATA.nt = static_cast<int>(DATA.tau_size / DATA.delta_tau + 0.5);
        }
        music_message << "Grid info: x_size = " << DATA.x_size
                      << ", y_size = " << DATA.y_size
                      << ", dx = " << DATA.delta_x
                      << " fm, dy = " << DATA.delta_y
                      << " fm, eta_min = " << DATA.etaMin()
                      << ", eta_max = " << DATA.etaMax()
                      << ", neta = " << DATA.neta
                      << ", deta = " << DATA.delta_eta
                      << ", dtau = " << DATA.delta_tau << " fm/c";
        music_message.flush("info");
    }
    compute_norm_for_strings(total_energy, netPz, minSigmaEta);
}

namespace {
constexpr double sourceCutoff = HydroSourceStrings::smearingCutoff;

double gaussian(double displacement) {
    return std::exp(-0.5 * displacement * displacement);
}

struct SourceAxis {
    double origin, spacing;
    int first, last;
    int originalLast;

    int firstWithin(double low) const {
        return std::max(
            first, static_cast<int>(std::ceil((low - origin) / spacing)));
    }
    int lastWithin(double high) const {
        return std::min(
            last, static_cast<int>(std::floor((high - origin) / spacing)));
    }
    double at(int index) const { return origin + index * spacing; }

    // A fixed grid can be a window into the event. Normalize a source
    // centered outside it on the continuation of the same lattice.
    bool includeOutside(
        double center, double padding, bool fullSupport = false) {
        if (!fullSupport && center >= origin && center <= at(originalLast))
            return false;
        first = std::min(
            first, static_cast<int>(
                       std::floor((center - padding - origin) / spacing)));
        last = std::max(
            last,
            static_cast<int>(std::ceil((center + padding - origin) / spacing)));
        return true;
    }

    double gaussianIntegral(double center, double sigma) const {
        const double radius = sourceCutoff * sigma;
        double sum = 0.;
        for (int i = firstWithin(center - radius);
             i <= lastWithin(center + radius); ++i)
            sum += gaussian((at(i) - center) / sigma);
        return sum * spacing;
    }
};

double transverseIntegral(
    const SourceAxis &x, const SourceAxis &y, double xc, double yc,
    double sigma, double flow) {
    if (flow == 0.)
        return x.gaussianIntegral(xc, sigma) * y.gaussianIntegral(yc, sigma);
    const double radius = sourceCutoff * sigma;
    const int yFirst = y.firstWithin(yc - radius);
    const int yLast = y.lastWithin(yc + radius);
    double sum = 0.;
    for (int ix = x.firstWithin(xc - radius); ix <= x.lastWithin(xc + radius);
         ++ix) {
        const double dx = (x.at(ix) - xc) / sigma;
        for (int iy = yFirst; iy <= yLast; ++iy) {
            const double dy = (y.at(iy) - yc) / sigma;
            const double r2 = dx * dx + dy * dy;
            sum += std::exp(-0.5 * r2) * std::cosh(flow * std::sqrt(r2));
        }
    }
    return sum * x.spacing * y.spacing;
}

double sourceNormalization(double strength, double integral) {
    if (strength == 0.) return 0.;
    if (!(integral > 0.) || !std::isfinite(strength / integral))
        throw std::runtime_error("String source has no usable grid support");
    return strength / integral;
}
}  // namespace

void HydroSourceStrings::compute_norm_for_strings(
    double event_energy, double event_pz, double minSigmaEta) {
    SourceAxis x {-DATA.x_size / 2., DATA.delta_x, 0, DATA.nx - 1, DATA.nx - 1};
    SourceAxis y {-DATA.y_size / 2., DATA.delta_y, 0, DATA.ny - 1, DATA.ny - 1};
    SourceAxis eta {
        DATA.etaMin(), DATA.delta_eta, 0, DATA.neta - 1, DATA.neta - 1};
    bool extended = false;
    if (DATA.boost_invariant) eta.spacing = minSigmaEta / 8.;

    if (!DATA.dynamicGridDetermination || DATA.boost_invariant) {
        for (const auto &it : QCD_strings_list) {
            const auto includeEta = [&](double center) {
                extended |= eta.includeOutside(
                    center, sourceCutoff * it->sigma_eta, DATA.boost_invariant);
            };
            includeEta(it->eta_s_left);
            includeEta(it->eta_s_right);
            if (it->baryon_frac_l != 0.) includeEta(it->eta_s_baryon_left);
            if (it->baryon_frac_r != 0.) includeEta(it->eta_s_baryon_right);
            if (!DATA.use_rhoQS_to_rhoB_ratios && it->Qe_left != 0.)
                includeEta(it->eta_s_Qe_left);
            if (!DATA.use_rhoQS_to_rhoB_ratios && it->Qe_right != 0.)
                includeEta(it->eta_s_Qe_right);
            if (!DATA.dynamicGridDetermination) {
                const double padding = sourceCutoff * it->sigma_x;
                extended |= x.includeOutside(it->x_pl, padding);
                extended |= x.includeOutside(it->x_pr, padding);
                extended |= y.includeOutside(it->y_pl, padding);
                extended |= y.includeOutside(it->y_pr, padding);
            }
        }
    }
    if (extended && !DATA.boost_invariant)
        music_message.warning(
            "The fixed grid excludes some source centers. "
            "Normalization includes the continuation of the grid.");

    struct EnergySource {
        double energy, plus, minus, rapidity;
        double &norm, &shift;
    };
    std::vector<EnergySource> strings, remnants;
    strings.reserve(QCD_strings_list.size());
    double stringEnergy = 0., stringPz = 0.;
    double remnantEnergy = 0., remnantPz = 0.;
    for (const auto &it : QCD_strings_list) {
        const double sy = it->sigma_eta;
        const double sx = it->sigma_x;
        const double width = it->eta_s_right - it->eta_s_left;
        const double radius = sourceCutoff * sy;
        const int first = eta.firstWithin(it->eta_s_left - radius);
        const int last = eta.lastWithin(it->eta_s_right + radius);
        const auto fraction = [&](double at) {
            return music_device::stringEtaFraction(
                at, it->eta_s_left, it->eta_s_right);
        };
        const auto profile = [&](double at) {
            return width == 0. ? gaussian((at - it->eta_s_left) / sy)
                               : music_device::stringEtaProfile(
                                     at, it->eta_s_left, it->eta_s_right, sy);
        };
        const double lossL =
            it->mass * (std::sinh(it->y_l_i) - std::sinh(it->y_l));
        const double lossR =
            it->mass * (std::sinh(it->y_r_i) - std::sinh(it->y_r));
        const double energy = it->mass
                              * (std::cosh(it->y_l_i) - std::cosh(it->y_l)
                                 + std::cosh(it->y_r_i) - std::cosh(it->y_r));
        if (energy < 0. || !std::isfinite(energy))
            throw std::runtime_error("Invalid energy loss in string source");
        stringEnergy += energy;
        stringPz += lossL + lossR;
        it->string_tilt = DATA.stringTransverseShiftFrac;
        if (DATA.string_tilt_update) {
            it->string_tilt = 0.;
            if (energy > 0.
                && std::abs(it->mass - MASS_PARTON) < MASS_PARTON_TOLERANCE) {
                double moment = 0., weight = 0.;
                for (int i = first; i <= last; ++i) {
                    const double frac = fraction(eta.at(i));
                    const double rap = it->y_l + (it->y_r - it->y_l) * frac;
                    const double shape = profile(eta.at(i));
                    moment += (2. * frac - 1.) * shape * std::sinh(rap);
                    weight += shape * std::cosh(rap);
                }
                if (std::abs(moment) > Util::small_eps) {
                    it->string_tilt =
                        (lossR - lossL) * weight / (energy * moment);
                    const double separation = std::max(
                        std::abs(it->x_pl - it->x_pr),
                        std::abs(it->y_pl - it->y_pr));
                    if (it->string_tilt < 0. && separation > 0.)
                        it->string_tilt =
                            std::max(it->string_tilt, -4. / separation);
                }
            }
        }

        double plus = 0., minus = 0.;
        if (energy > 0.) {
            const auto transverse = [&](double frac) {
                return transverseIntegral(
                    x, y,
                    getStringTransverseCoord(
                        it->x_pl, it->x_pr, frac, it->string_tilt),
                    getStringTransverseCoord(
                        it->y_pl, it->y_pr, frac, it->string_tilt),
                    sx, preEqFlowFactor_);
            };
            const double atLeft = transverse(0.);
            const double atRight =
                it->string_tilt == 0. ? atLeft : transverse(1.);
            for (int i = first; i <= last; ++i) {
                const double at = eta.at(i);
                const double frac = fraction(at);
                const double xy =
                    (frac == 0. || it->string_tilt == 0.)
                        ? atLeft
                        : (frac == 1. ? atRight : transverse(frac));
                const double shape = profile(at) * eta.spacing * xy;
                const double rap = it->y_l + (it->y_r - it->y_l) * frac;
                plus += shape * std::exp(rap);
                minus += shape * std::exp(-rap);
            }
            if (!(plus > 0. && minus > 0.) || !std::isfinite(plus + minus))
                throw std::runtime_error(
                    "String energy has no usable grid support");
            strings.push_back(
                {energy, plus, minus, 0.5 * (std::log(plus) - std::log(minus)),
                 it->norm, it->rapidity_shift});
        }

        const auto normalizeEnd =
            [&](double xc, double yc, double stringEnd, double rap,
                double remnant, double baryonEnd, double baryon,
                double chargeEnd, double charge, double &remNorm,
                double &remShift, double &bNorm, double &qNorm) {
                const double remnantMass = remnant * it->mass;
                if (remnantMass != 0.) {
                    const double e = remnantMass * std::cosh(rap);
                    remnantEnergy += e;
                    remnantPz += remnantMass * std::sinh(rap);
                    const double integral =
                        transverseIntegral(x, y, xc, yc, sx, preEqFlowFactor_)
                        * eta.gaussianIntegral(stringEnd, sy);
                    remnants.push_back(
                        {e, integral * std::exp(rap), integral * std::exp(-rap),
                         rap, remNorm, remShift});
                }
                if (baryon != 0. || charge != 0.) {
                    const double xy =
                        x.gaussianIntegral(xc, sx) * y.gaussianIntegral(yc, sx);
                    if (baryon != 0.)
                        bNorm = sourceNormalization(
                            baryon, xy * eta.gaussianIntegral(baryonEnd, sy));
                    if (charge != 0.)
                        qNorm = sourceNormalization(
                            charge, xy * eta.gaussianIntegral(chargeEnd, sy));
                }
            };
        normalizeEnd(
            it->x_pl, it->y_pl, it->eta_s_left, it->y_l, it->remnant_l,
            it->eta_s_baryon_left, it->baryon_frac_l, it->eta_s_Qe_left,
            DATA.use_rhoQS_to_rhoB_ratios ? 0. : it->Qe_left,
            it->E_remnant_norm_L, it->remnant_shift_L, it->B_norm_L,
            it->Q_norm_L);
        normalizeEnd(
            it->x_pr, it->y_pr, it->eta_s_right, it->y_r, it->remnant_r,
            it->eta_s_baryon_right, it->baryon_frac_r, it->eta_s_Qe_right,
            DATA.use_rhoQS_to_rhoB_ratios ? 0. : it->Qe_right,
            it->E_remnant_norm_R, it->remnant_shift_R, it->B_norm_R,
            it->Q_norm_R);
    }

    if (!std::isfinite(event_energy))
        event_energy = stringEnergy + remnantEnergy;
    if (!std::isfinite(event_pz)) event_pz = stringPz + remnantPz;
    const double tolerance = 1e-9 * std::max(1., event_energy);
    const auto fit = [&](const std::vector<EnergySource> &sources,
                         double inputEnergy, double targetEnergy,
                         double targetPz) {
        if (std::abs(targetEnergy) <= tolerance
            && std::abs(targetPz) <= tolerance) {
            // nothing is left for these sources to deposit
            for (const auto &source : sources) {
                source.norm = 0.;
                source.shift = 0.;
            }
            if (inputEnergy > tolerance)
                music_message.warning(
                    "The event four-momentum leaves no energy for the "
                    "string sources: their energy loss is not deposited.");
            return true;
        }
        if (inputEnergy <= 0. || targetEnergy <= std::abs(targetPz))
            return false;

        // The energy shares stay fixed as a common shift in momentum
        // rapidity controls the total Pz monotonically.
        const double targetVelocity = targetPz / targetEnergy;
        const auto velocity = [&](double shift) {
            double sum = 0.;
            for (const auto &source : sources)
                sum += (source.energy / inputEnergy)
                       * std::tanh(source.rapidity + shift);
            return sum;
        };
        double low = -1., high = 1., shift = 0.;
        while (velocity(low) > targetVelocity && low > -2048.) low *= 2.;
        while (velocity(high) < targetVelocity && high < 2048.) high *= 2.;
        for (int i = 0; i < 64; ++i) {
            shift = 0.5 * (low + high);
            const double current = velocity(shift);
            if (std::abs(current - targetVelocity) < 4e-15) break;
            if (current < targetVelocity)
                low = shift;
            else
                high = shift;
        }
        const double boost = std::exp(shift);
        for (const auto &source : sources) {
            source.shift = shift;
            source.norm = sourceNormalization(
                (source.energy / inputEnergy) * targetEnergy,
                0.5 * (source.plus * boost + source.minus / boost));
        }
        return true;
    };

    // Preserve the remnant four-momenta whenever the remaining string
    // four-momentum can be represented by a positive energy source.
    if (fit(strings, stringEnergy, event_energy - remnantEnergy,
            event_pz - remnantPz)) {
        for (const auto &remnant : remnants)
            remnant.norm = sourceNormalization(
                remnant.energy, 0.5 * (remnant.plus + remnant.minus));
    } else {
        for (const auto &remnant : remnants) strings.push_back(remnant);
        if (!fit(strings, stringEnergy + remnantEnergy, event_energy, event_pz))
            throw std::runtime_error(
                "Invalid total energy or net Pz in string event");
        music_message.warning(
            "Matching the event four-momentum also shifts remnant rapidities.");
    }
    music_message << "String source totals: E = " << event_energy
                  << " GeV, net Pz = " << event_pz << " GeV.";
    music_message.flush("info");
}

void HydroSourceStrings::prepare_list_for_current_tau_frame(
    const double tau_local) {
    const double dtau = DATA.delta_tau;
    const double begin = stringBinStart(tau_local, dtau);
    const double end = stringBinEnd(tau_local + dtau, dtau);
    const auto active = [=](double time, double norm) {
        return norm != 0. && time >= begin && time < end;
    };
    QCD_strings_list_current_tau.clear();
    QCD_strings_baryon_list_current_tau.clear();
    QCD_strings_electric_list_current_tau.clear();
    QCD_strings_remnant_list_current_tau.clear();
    for (const auto &it : QCD_strings_list) {
        if (active(it->tau_baryon_left, it->B_norm_L)
            || active(it->tau_baryon_right, it->B_norm_R))
            QCD_strings_baryon_list_current_tau.push_back(it);
        if (active(it->tau_Qe_left, it->Q_norm_L)
            || active(it->tau_Qe_right, it->Q_norm_R))
            QCD_strings_electric_list_current_tau.push_back(it);
        if (active(it->tau_end_left, it->E_remnant_norm_L)
            || active(it->tau_end_right, it->E_remnant_norm_R))
            QCD_strings_remnant_list_current_tau.push_back(it);
        if (it->norm != 0. && it->tau_start <= end
            && std::max(it->tau_end_left, it->tau_end_right) >= begin)
            QCD_strings_list_current_tau.push_back(it);
    }
    // only for a time step in which strings deposit
    if (QCD_strings_list_current_tau.empty()
        && QCD_strings_baryon_list_current_tau.empty()
        && QCD_strings_electric_list_current_tau.empty()
        && QCD_strings_remnant_list_current_tau.empty())
        return;
    music_message << "hydro_source: tau = " << tau_local << " fm."
                  << " number of strings for energy density: "
                  << QCD_strings_list_current_tau.size()
                  << " number of strings for net baryon density: "
                  << QCD_strings_baryon_list_current_tau.size()
                  << " number of strings for net electric charge density: "
                  << QCD_strings_electric_list_current_tau.size();
    music_message.flush("info");
}

void HydroSourceStrings::get_hydro_energy_source(
    const double tau, const double x, const double y, const double eta_s,
    const FlowVec &u_mu, EnergyFlowVec &j_mu) const {
    j_mu = {0};
    if (QCD_strings_list_current_tau.size() == 0
        && QCD_strings_remnant_list_current_tau.size() == 0)
        return;

    const double dtau = DATA.delta_tau;
    const double n_sigma_skip = smearingCutoff;
    const double exp_tau = 1. / tau;
    for (auto const &it : QCD_strings_list_current_tau) {
        const double sigma_x = it->sigma_x;
        if (it->norm == 0.) continue;
        const double sigma_eta = it->sigma_eta;
        const double skip_dis_x = n_sigma_skip * sigma_x;
        const double skip_dis_eta = n_sigma_skip * sigma_eta;
        if (eta_s < it->eta_s_left - skip_dis_eta
            || eta_s > it->eta_s_right + skip_dis_eta)
            continue;
        const double eta_frac = music_device::stringEtaFraction(
            eta_s, it->eta_s_left, it->eta_s_right);
        const double x_perp = getStringTransverseCoord(
            it->x_pl, it->x_pr, eta_frac, it->string_tilt);
        double x_dis = x - x_perp;
        if (std::abs(x_dis) > skip_dis_x) continue;

        const double y_perp = getStringTransverseCoord(
            it->y_pl, it->y_pr, eta_frac, it->string_tilt);
        double y_dis = y - y_perp;
        if (std::abs(y_dis) > skip_dis_x) continue;

        const double exp_eta_s = music_device::stringEtaProfileForStep(
            eta_s, it->eta_s_left, it->eta_s_right, it->eta_s_0, sigma_eta, tau,
            dtau, it->tau_0, it->tau_form, it->tau_end_left,
            DATA.stringSourceTauMax);
        if (exp_eta_s == 0.) continue;

        double exp_xperp =
            exp(-(x_dis * x_dis + y_dis * y_dis) / (2. * sigma_x * sigma_x));
        double cosh_perp = (cosh(
            preEqFlowFactor_ * sqrt(x_dis * x_dis + y_dis * y_dis) / sigma_x));
        double sinh_perp = (sinh(
            preEqFlowFactor_ * sqrt(x_dis * x_dis + y_dis * y_dis) / sigma_x));
        double phi_perp = atan2(y_dis, x_dis);

        double e_local = exp_tau * exp_xperp * exp_eta_s * it->norm;
        const double y_string =
            it->y_l + (it->y_r - it->y_l) * eta_frac + it->rapidity_shift;
        const double cosh_long = cosh(y_string - eta_s);
        const double sinh_long = sinh(y_string - eta_s);
        const double local_eperp = e_local;
        j_mu[0] += local_eperp * cosh_long * cosh_perp;
        j_mu[1] += local_eperp * sinh_perp * cos(phi_perp);
        j_mu[2] += local_eperp * sinh_perp * sin(phi_perp);
        j_mu[3] += local_eperp * sinh_long * cosh_perp;
    }

    for (auto const &it : QCD_strings_remnant_list_current_tau) {
        const double sigma_x = it->sigma_x;
        const double sigma_eta = it->sigma_eta;
        const double skip_dis_x = n_sigma_skip * sigma_x;
        const double skip_dis_eta = n_sigma_skip * sigma_eta;
        for (int side = 0; side < 2; side++) {
            const double tauEnd =
                (side == 0 ? it->tau_end_left : it->tau_end_right);
            const double norm =
                (side == 0 ? it->E_remnant_norm_L : it->E_remnant_norm_R);
            if (norm == 0. || tauEnd < stringBinStart(tau, dtau)
                || tauEnd >= stringBinEnd(tau, dtau))
                continue;
            const double eta_dis =
                eta_s - (side == 0 ? it->eta_s_left : it->eta_s_right);
            const double x_dis = x - (side == 0 ? it->x_pl : it->x_pr);
            const double y_dis = y - (side == 0 ? it->y_pl : it->y_pr);
            if (std::abs(eta_dis) > skip_dis_eta || std::abs(x_dis) > skip_dis_x
                || std::abs(y_dis) > skip_dis_x)
                continue;
            const double exp_eta =
                exp(-eta_dis * eta_dis / (2. * sigma_eta * sigma_eta));
            const double r2 = x_dis * x_dis + y_dis * y_dis;
            const double exp_xperp = exp(-r2 / (2. * sigma_x * sigma_x));
            const double flow = preEqFlowFactor_ * sqrt(r2) / sigma_x;
            const double yEndRap =
                (side == 0 ? it->y_l + it->remnant_shift_L
                           : it->y_r + it->remnant_shift_R);
            const double local_e = exp_tau * exp_eta * norm * exp_xperp;
            const double phi = atan2(y_dis, x_dis);
            j_mu[0] += local_e * cosh(yEndRap - eta_s) * cosh(flow);
            j_mu[1] += local_e * sinh(flow) * cos(phi);
            j_mu[2] += local_e * sinh(flow) * sin(phi);
            j_mu[3] += local_e * sinh(yEndRap - eta_s) * cosh(flow);
        }
    }
    const double prefactor_tau = 1. / dtau;
    const double unit_convert = 1.0 / Util::hbarc;
    const double prefactors = prefactor_tau * unit_convert;
    j_mu[0] *= prefactors;  // 1/fm^4
    j_mu[1] *= prefactors;  // 1/fm^4
    j_mu[2] *= prefactors;  // 1/fm^4
    j_mu[3] *= prefactors;  // 1/fm^4
}

double HydroSourceStrings::get_hydro_rhob_source(
    const double tau, const double x, const double y, const double eta_s,
    const FlowVec &u_mu) const {
    double res = 0.;
    if (QCD_strings_baryon_list_current_tau.size() == 0) return (res);

    // flow velocity
    const double gamma_perp_flow =
        sqrt(1. + u_mu[1] * u_mu[1] + u_mu[2] * u_mu[2]);
    const double y_perp_flow = acosh(gamma_perp_flow);
    const double y_long_flow = asinh(u_mu[3] / gamma_perp_flow) + eta_s;
    // u^x = sinh(y_perp) cos(phi), u^y = sinh(y_perp) sin(phi)
    const double u_perp_flow = sqrt(u_mu[1] * u_mu[1] + u_mu[2] * u_mu[2]);
    const double cos_phi_flow =
        (u_perp_flow > Util::small_eps ? u_mu[1] / u_perp_flow : 1.);
    const double sin_phi_flow =
        (u_perp_flow > Util::small_eps ? u_mu[2] / u_perp_flow : 0.);
    const double dtau = DATA.delta_tau;

    const double exp_tau = 1.0 / tau;
    const double n_sigma_skip = smearingCutoff;
    for (auto &it : QCD_strings_baryon_list_current_tau) {
        const double sigma_x = it->sigma_x;
        const double sigma_eta = it->sigma_eta;
        const double skip_dis_x = n_sigma_skip * sigma_x;
        const double skip_dis_eta = n_sigma_skip * sigma_eta;

        // skip the evaluation if the strings is too far away in the
        // space-time grid
        // dumping energy into the medium from the active strings
        // double tau_dis_left = fabs(tau - it->tau_end_left);
        // double tau_dis_right = fabs(tau - it->tau_end_right);
        int flag_left = 0;
        if (it->tau_baryon_left >= stringBinStart(tau, dtau)
            && it->tau_baryon_left < stringBinEnd(tau, dtau)
            && it->baryon_frac_l != 0.) {
            flag_left = 1;
        }

        int flag_right = 0;
        if (it->tau_baryon_right >= stringBinStart(tau, dtau)
            && it->tau_baryon_right < stringBinEnd(tau, dtau)
            && it->baryon_frac_r != 0.) {
            flag_right = 1;
        }

        if (flag_left == 0 && flag_right == 0) continue;

        if (flag_left == 1) {
            double eta_dis = eta_s - it->eta_s_baryon_left;
            if (std::abs(eta_dis) <= skip_dis_eta) {
                double exp_eta_left =
                    exp(-(eta_dis * eta_dis) / (2. * sigma_eta * sigma_eta));

                double x_dis = x - it->x_pl;
                if (std::abs(x_dis) <= skip_dis_x) {
                    double y_dis = y - it->y_pl;
                    if (std::abs(y_dis) <= skip_dis_x) {
                        double exp_xperp =
                            exp(-(x_dis * x_dis + y_dis * y_dis)
                                / (2. * sigma_x * sigma_x));

                        double fsmear_left =
                            exp_tau
                            * (exp_xperp * exp_eta_left * (it->B_norm_L));

                        if (fsmear_left != 0.) {
                            double rapidity_local = it->y_l_baryon;

                            double y_dump =
                                ((1. - parton_quench_factor) * rapidity_local
                                 + parton_quench_factor * y_long_flow);

                            double y_dump_perp =
                                parton_quench_factor * y_perp_flow;

                            double p_dot_u = 1.;
                            if (parton_quench_factor < 1.) {
                                p_dot_u =
                                    (u_mu[0] * cosh(y_dump - eta_s)
                                         * cosh(y_dump_perp)
                                     - u_mu[1] * sinh(y_dump_perp)
                                           * cos_phi_flow
                                     - u_mu[2] * sinh(y_dump_perp)
                                           * sin_phi_flow
                                     - u_mu[3] * sinh(y_dump - eta_s)
                                           * cosh(y_dump_perp));
                            }

                            res += p_dot_u * fsmear_left;
                        }
                    }
                }
            }
        }

        if (flag_right == 1) {
            double eta_dis = eta_s - it->eta_s_baryon_right;
            if (std::abs(eta_dis) <= skip_dis_eta) {
                double exp_eta_right =
                    exp(-(eta_dis * eta_dis) / (2. * sigma_eta * sigma_eta));

                double x_dis = x - it->x_pr;
                if (std::abs(x_dis) <= skip_dis_x) {
                    double y_dis = y - it->y_pr;
                    if (std::abs(y_dis) <= skip_dis_x) {
                        double exp_xperp =
                            exp(-(x_dis * x_dis + y_dis * y_dis)
                                / (2. * sigma_x * sigma_x));

                        double fsmear_right =
                            exp_tau
                            * (exp_xperp * exp_eta_right * (it->B_norm_R));

                        if (fsmear_right != 0.) {
                            double rapidity_local = it->y_r_baryon;

                            double y_dump =
                                ((1. - parton_quench_factor) * rapidity_local
                                 + parton_quench_factor * y_long_flow);

                            double y_dump_perp =
                                parton_quench_factor * y_perp_flow;

                            double p_dot_u = 1.;
                            if (parton_quench_factor < 1.) {
                                p_dot_u =
                                    (u_mu[0] * cosh(y_dump - eta_s)
                                         * cosh(y_dump_perp)
                                     - u_mu[1] * sinh(y_dump_perp)
                                           * cos_phi_flow
                                     - u_mu[2] * sinh(y_dump_perp)
                                           * sin_phi_flow
                                     - u_mu[3] * sinh(y_dump - eta_s)
                                           * cosh(y_dump_perp));
                            }

                            res += p_dot_u * fsmear_right;
                        }
                    }
                }
            }
        }
    }

    const double prefactor_tau = 1. / dtau;
    res *= prefactor_tau;
    return (res);
}

double HydroSourceStrings::get_hydro_rhoQ_source(
    const double tau, const double x, const double y, const double eta_s,
    const FlowVec &u_mu) const {
    double res = 0.;
    if (DATA.use_rhoQS_to_rhoB_ratios) {
        double nb = get_hydro_rhob_source(tau, x, y, eta_s, u_mu);
        res = DATA.ratio_q * nb;
        return (res);
    }
    if (QCD_strings_electric_list_current_tau.size() == 0) return (res);

    // flow velocity
    const double gamma_perp_flow =
        sqrt(1. + u_mu[1] * u_mu[1] + u_mu[2] * u_mu[2]);
    const double y_perp_flow = acosh(gamma_perp_flow);
    const double y_long_flow = asinh(u_mu[3] / gamma_perp_flow) + eta_s;
    // u^x = sinh(y_perp) cos(phi), u^y = sinh(y_perp) sin(phi)
    const double u_perp_flow = sqrt(u_mu[1] * u_mu[1] + u_mu[2] * u_mu[2]);
    const double cos_phi_flow =
        (u_perp_flow > Util::small_eps ? u_mu[1] / u_perp_flow : 1.);
    const double sin_phi_flow =
        (u_perp_flow > Util::small_eps ? u_mu[2] / u_perp_flow : 0.);
    const double dtau = DATA.delta_tau;

    // the charges are deposited at rest in the fluid (y = y_long_flow);
    // only the transverse flow is quenched
    const double y_dump = y_long_flow;
    const double y_dump_perp = parton_quench_factor * y_perp_flow;
    double p_dot_u = 1.;
    if (parton_quench_factor < 1.) {
        p_dot_u =
            (u_mu[0] * cosh(y_dump - eta_s) * cosh(y_dump_perp)
             - u_mu[1] * sinh(y_dump_perp) * cos_phi_flow
             - u_mu[2] * sinh(y_dump_perp) * sin_phi_flow
             - u_mu[3] * sinh(y_dump - eta_s) * cosh(y_dump_perp));
    }

    const double exp_tau = 1.0 / tau;
    const double n_sigma_skip = smearingCutoff;
    for (auto &it : QCD_strings_electric_list_current_tau) {
        const double sigma_x = it->sigma_x;
        const double sigma_eta = it->sigma_eta;
        const double skip_dis_x = n_sigma_skip * sigma_x;
        const double skip_dis_eta = n_sigma_skip * sigma_eta;

        const bool flag_left =
            (it->tau_Qe_left >= stringBinStart(tau, dtau)
             && it->tau_Qe_left < stringBinEnd(tau, dtau) && it->Qe_left != 0.);
        const bool flag_right =
            (it->tau_Qe_right >= stringBinStart(tau, dtau)
             && it->tau_Qe_right < stringBinEnd(tau, dtau)
             && it->Qe_right != 0.);
        if (!flag_left && !flag_right) continue;

        // each string end deposits its charge, as for the net baryon number
        // (see get_hydro_rhob_source)
        auto addEnd = [&](const double eta_end, const double x_end,
                          const double y_end, const double Qe) {
            const double eta_dis = eta_s - eta_end;
            if (std::abs(eta_dis) > skip_dis_eta) return;
            const double x_dis = x - x_end;
            if (std::abs(x_dis) > skip_dis_x) return;
            const double y_dis = y - y_end;
            if (std::abs(y_dis) > skip_dis_x) return;
            const double exp_eta =
                exp(-(eta_dis * eta_dis) / (2. * sigma_eta * sigma_eta));
            const double exp_xperp = exp(
                -(x_dis * x_dis + y_dis * y_dis) / (2. * sigma_x * sigma_x));
            res += p_dot_u * exp_tau * exp_xperp * exp_eta * Qe;
        };
        if (flag_left)
            addEnd(it->eta_s_Qe_left, it->x_pl, it->y_pl, it->Q_norm_L);
        if (flag_right) {
            addEnd(it->eta_s_Qe_right, it->x_pr, it->y_pr, it->Q_norm_R);
        }
    }
    const double prefactor_tau = 1. / dtau;
    res *= prefactor_tau;
    return (res);
}

double HydroSourceStrings::getStringTransverseCoord(
    const double xl, const double xr, const double etaFrac,
    const double string_tilt) const {
    const double x_c = 0.5 * (xl + xr);

    const double S_transverse_position =
        x_c + string_tilt * (0.5 - etaFrac) * (xl - xr);

    double xT;
    double sub_nucleon_distance = 2.0;  // fm

    xT = S_transverse_position;

    if (string_tilt > 0.0) {
        if ((S_transverse_position - x_c) > sub_nucleon_distance) {
            xT = x_c + sub_nucleon_distance;
        }

        else if ((S_transverse_position - x_c) < -sub_nucleon_distance) {
            xT = x_c - sub_nucleon_distance;
        }
    }
    return xT;
}
