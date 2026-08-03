#include "eos_chem.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "util.h"

using std::string;
using std::stringstream;

namespace {
// grid parameters of the table produced by
// eos_construction/code/build_full_range_eos_table.py
const int N_E_CHEM = 100000;
const int N_YQ_CHEM = 11;
const double YQ_MIN_CHEM = 0.0;
const double YQ_SPACING_CHEM = 0.1;
const double YQ_TOLERANCE_CHEM = 1e-6;
}  // namespace

EOS_chem::EOS_chem(const int eos_id_in) : eos_id(eos_id_in) {
    set_EOS_id(eos_id);
    set_number_of_tables(0);
    set_eps_max(1e5);
    set_flag_muB(false);
    set_flag_muS(false);
    set_flag_muQ(false);
    entropy_tb = nullptr;
}

EOS_chem::~EOS_chem() {
    int ntables = get_number_of_tables();
    for (int itable = 0; itable < ntables; itable++) {
        Util::mtx_free(
            entropy_tb[itable], nb_length[itable], e_length[itable]);
    }
    if (ntables > 0) {
        delete[] entropy_tb;
    }
}

void EOS_chem::initialize_eos() {
    // read the Y_q-dependent chemical equilibration EOS table
    music_message.info("reading EOS chem (Y_q dependent EOS) ...");

    auto envPath = get_hydro_env_path();
    stringstream slocalpath;
    slocalpath << envPath << "/EOS/chem_eq";

    string path = slocalpath.str();
    music_message << "from path " << path;
    music_message.flush("info");

    set_number_of_tables(1);
    resize_table_info_arrays();

    int ntables = get_number_of_tables();

    pressure_tb = new double **[ntables];
    temperature_tb = new double **[ntables];
    entropy_tb = new double **[ntables];
    for (int itable = 0; itable < ntables; itable++) {
        std::ifstream eos_file;
        eos_file.open(path + "/eos_s_T_vs_e_Yq_full_range.dat");

        if (!eos_file) {
            music_message.error("Can not find the chem EoS file.");
            exit(1);
        }

        // skip the comment lines (starting with #) and the one
        // column-name header line that follow them
        string dummy;
        while (std::getline(eos_file, dummy)) {
            if (!dummy.empty() && dummy[0] != '#') break;
        }

        e_length[itable] = N_E_CHEM;
        nb_length[itable] = N_YQ_CHEM;
        nb_bounds[itable] = YQ_MIN_CHEM;
        nb_spacing[itable] = YQ_SPACING_CHEM;

        pressure_tb[itable] =
            Util::mtx_malloc(nb_length[itable], e_length[itable]);
        temperature_tb[itable] =
            Util::mtx_malloc(nb_length[itable], e_length[itable]);
        entropy_tb[itable] =
            Util::mtx_malloc(nb_length[itable], e_length[itable]);

        double yq_val, e_val, p_val, s_val, t_val;
        for (int iyq = 0; iyq < N_YQ_CHEM; iyq++) {
            double yq_expected = YQ_MIN_CHEM + iyq * YQ_SPACING_CHEM;
            for (int ie = 0; ie < N_E_CHEM; ie++) {
                eos_file >> yq_val >> e_val >> p_val >> s_val >> t_val;
                if (std::abs(yq_val - yq_expected) > YQ_TOLERANCE_CHEM) {
                    music_message
                        << "chem EoS file Y_q ordering mismatch: expected "
                        << yq_expected << " got " << yq_val;
                    music_message.flush("error");
                    exit(1);
                }

                e_val /= Util::hbarc;  // [1/fm^4]
                if (iyq == 0 && ie == 0) e_bounds[itable] = e_val;
                if (iyq == 0 && ie == 1) {
                    e_spacing[itable] = e_val - e_bounds[itable];
                }
                if (iyq == 0 && ie == N_E_CHEM - 1) set_eps_max(e_val);

                pressure_tb[itable][iyq][ie] = p_val / Util::hbarc;  // 1/fm^4
                entropy_tb[itable][iyq][ie] = s_val;                 // 1/fm^3
                temperature_tb[itable][iyq][ie] =
                    t_val / Util::hbarc;  // 1/fm
            }
        }

        if (!eos_file.good() && !eos_file.eof()) {
            music_message.error(
                "chem EoS file ended unexpectedly while reading.");
            exit(1);
        }
    }
    music_message.info("Done reading EOS chem.");
}

// the rhob argument in the functions below is accepted for interface
// compatibility and ignored: this EOS has no baryon density dependence
// (see eos_chem.h). Y_q is passed to the generic 2D table machinery in
// the axis slot that stores Y_q for this table.
double EOS_chem::p_e_func(double e, double rhob, double Y_q) const {
    return (get_dpOverde3(e, rhob, Y_q));
}

// returns the local temperature in [1/fm]
// input local energy density e [1/fm^4] and Y_q in [0,1]
double EOS_chem::get_temperature(double e, double rhob, double Y_q) const {
    double T = interpolate2D(e, Y_q, 0, temperature_tb);  // 1/fm
    return (std::max(Util::small_eps, T));
}

// returns the local pressure in [1/fm^4]
// input local energy density e [1/fm^4] and Y_q in [0,1]
double EOS_chem::get_pressure(double e, double rhob, double Y_q) const {
    double f = interpolate2D(e, Y_q, 0, pressure_tb);  // 1/fm^4
    return (std::max(Util::small_eps, f));
}

// returns the local entropy density in [1/fm^3], read directly from the
// tabulated s(e,Y_q) rather than derived from P and T (see the note in
// eos_chem.h about this hiding, not overriding, EOS_base::get_entropy)
double EOS_chem::get_entropy(double e, double rhob, double Y_q) const {
    double s = interpolate2D(e, Y_q, 0, entropy_tb);  // 1/fm^3
    return (std::max(Util::small_eps, s));
}

// returns the local pressure in [1/fm^4] together with dP/de and dP/dY_q
// evaluated from the tabulated P(e,Y_q); dP/dY_q is returned in the
// dpdrhob output slot of the generic interface, since this EOS has no
// rhob dependence of its own
void EOS_chem::get_pressure_with_gradients(
    double e, double rhob, double &p, double &dpde, double &dpdYq,
    double &cs2, double Y_q) const {
    interpolate2D_with_gradients(e, Y_q, 0, pressure_tb, p, dpde, dpdYq);
    p = std::max(Util::small_eps, p);  // [1/fm^4]
    cs2 = std::max(0.01, std::min(1. / 3., dpde));
}

double EOS_chem::get_s2e(double s, double rhob, double Y_q) const {
    double e = get_s2e_finite_rhob(s, rhob, Y_q);
    return (e);
}

double EOS_chem::get_T2e(double T_in_GeV, double rhob, double Y_q) const {
    double e = get_T2e_finite_rhob(T_in_GeV, rhob, Y_q);
    return (e);
}
