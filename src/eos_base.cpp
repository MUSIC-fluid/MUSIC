
#include "eos_base.h"
#include "util.h"

#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>

using std::ostringstream;
using std::setw;
using std::setprecision;
using std::scientific;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using Util::hbarc;

EOS_base::~EOS_base() {
    for (int itable = 0; itable < number_of_tables; itable++) {
        Util::mtx_free(pressure_tb[itable],
                       nb_length[itable], e_length[itable]);
        Util::mtx_free(temperature_tb[itable],
                       nb_length[itable], e_length[itable]);
    }
    if (number_of_tables > 0) {
        delete[] pressure_tb;
        delete[] temperature_tb;
    }
}


double EOS_base::interpolate1D(double e, int table_idx, double ***table) const {
// This is a generic linear interpolation routine for EOS at zero mu_B
// it assumes the class has already read in
//        P(e), T(e), s(e)
// as one-dimensional arrays on an equally spacing lattice grid
// units: e is in 1/fm^4
    //double local_ed = e*hbarc;  // [GeV/fm^3]
    double local_ed = e;

    const double e0       = e_bounds[table_idx];
    const double delta_e  = e_spacing[table_idx];
    const int N_e         = e_length[table_idx];

    // compute the indices
    int idx_e  = static_cast<int>((local_ed - e0)/delta_e);

    // treatment for overflow, use the last two points to do extrapolation
    idx_e  = std::min(N_e - 2, idx_e);

    // check underflow
    idx_e  = std::max(0, idx_e);

    const double frac_e = (local_ed - (idx_e*delta_e + e0))/delta_e;

    double result;
    double temp1 = table[table_idx][0][idx_e];
    double temp2 = table[table_idx][0][idx_e + 1];
    result = temp1*(1. - frac_e) + temp2*frac_e;
    return(result);
}


double EOS_base::interpolate2D(double e, double rhob, int table_idx, double ***table) const {
// This is a generic bilinear interpolation routine for EOS at finite mu_B
// it assumes the class has already read in
//        P(e, rho_b), T(e, rho_b), s(e, rho_b), mu_b(e, rho_b)
// as two-dimensional arrays on an equally spacing lattice grid
// units: e is in 1/fm^4, rhob is in 1/fm^3
    //double local_ed = e*hbarc;  // [GeV/fm^3]
    double local_ed = e;
    double local_nb = rhob;     // [1/fm^3]

    double e0       = e_bounds[table_idx];
    double nb0      = nb_bounds[table_idx];
    double delta_e  = e_spacing[table_idx];
    double delta_nb = nb_spacing[table_idx];

    int N_e  = e_length[table_idx];
    int N_nb = nb_length[table_idx];

    // compute the indices
    int idx_e  = static_cast<int>((local_ed - e0)/delta_e);
    int idx_nb = static_cast<int>((local_nb - nb0)/delta_nb);

    // treatment for overflow, use the last two points to do extrapolation
    idx_e  = std::min(N_e - 2, idx_e);
    idx_nb = std::min(N_nb - 2, idx_nb);

    // check underflow
    idx_e  = std::max(0, idx_e);
    idx_nb = std::max(0, idx_nb);

    double frac_e    = (local_ed - (idx_e*delta_e + e0))/delta_e;
    double frac_rhob = (local_nb - (idx_nb*delta_nb + nb0))/delta_nb;

    double result;
    double temp1 = table[table_idx][idx_nb][idx_e];
    double temp2 = table[table_idx][idx_nb][idx_e + 1];
    double temp3 = table[table_idx][idx_nb + 1][idx_e + 1];
    double temp4 = table[table_idx][idx_nb + 1][idx_e];
    result = ((temp1*(1. - frac_e) + temp2*frac_e)*(1. - frac_rhob)
              + (temp3*frac_e + temp4*(1. - frac_e))*frac_rhob);
    return(result);
}


//! This function returns entropy density in [1/fm^3]
//! The input local energy density e [1/fm^4], rhob[1/fm^3]
double EOS_base::get_entropy(double epsilon, double rhob) const {
    auto P    = get_pressure(epsilon, rhob);
    auto T    = get_temperature(epsilon, rhob);
    auto muB  = get_muB(epsilon, rhob);
    auto muS  = get_muS(epsilon, rhob);
    auto muC  = get_muC(epsilon, rhob);
    auto rhoS = get_rhoS(epsilon, rhob);
    auto rhoC = get_rhoC(epsilon, rhob);
    auto f    = (epsilon + P - muB*rhob - muS*rhoS - muC*rhoC)/(T + 1e-15);
    return(std::max(1e-15, f));
}


double EOS_base::get_cs2(double e, double rhob) const {
    double f = calculate_velocity_of_sound_sq(e, rhob);
    return(f);
}


double EOS_base::calculate_velocity_of_sound_sq(double e, double rhob) const {
    double v_min = 0.01;
    double v_max = 1./3;
    double dpde = p_e_func(e, rhob);
    double dpdrho = p_rho_func(e, rhob);
    double pressure = get_pressure(e, rhob);
    double v_sound = dpde + rhob/(e + pressure + 1e-15)*dpdrho;
    v_sound = std::max(v_min, std::min(v_max, v_sound));
    return(v_sound);
}


double EOS_base::get_dpOverde3(double e, double rhob) const {
   double eLeft = 0.9*e;
   double eRight = 1.1*e;

   double pL = get_pressure(eLeft, rhob);   // 1/fm^4
   double pR = get_pressure(eRight, rhob);  // 1/fm^4
      
   double dpde = (pR - pL)/(eRight - eLeft);
   return dpde;
}


double EOS_base::get_dpOverdrhob2(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double deltaRhob = nb_spacing[table_idx];
    //double rhob_max = nb_bounds[table_idx] + nb_length[table_idx]*deltaRhob;
    
    double rhobLeft  = rhob - deltaRhob*0.5;
    double rhobRight = rhob + deltaRhob*0.5;

    double pL = get_pressure(e, rhobLeft);      // 1/fm^4
    double pR = get_pressure(e, rhobRight);     // 1/fm^4
      
    double dpdrho = (pR - pL)/(rhobRight - rhobLeft);  // 1/fm
    return (dpdrho);   // in 1/fm
}


int EOS_base::get_table_idx(double e) const {
    //double local_ed = e*hbarc;  // [GeV/fm^3]
    double local_ed = e;  // [GeV/fm^3]
    for (int itable = 1; itable < number_of_tables; itable++) {
        if (local_ed < e_bounds[itable]) {
            return(itable - 1);
        }
    }
    return(std::max(0, number_of_tables - 1));
}


//! This function returns local energy density [1/fm^4] from
//! a given temperature T [GeV] and rhob [1/fm^3] using binary search
double EOS_base::get_T2e_finite_rhob(const double T, const double rhob) const {
    double T_goal = T/Util::hbarc;         // convert to 1/fm
    double eps_lower = 1e-15;
    double eps_upper = eps_max;
    double eps_mid   = (eps_upper + eps_lower)/2.;
    double T_lower   = get_temperature(eps_lower, rhob);
    double T_upper   = get_temperature(eps_upper, rhob);
    int ntol         = 1000;
    if (T_goal < 0.0 || T_goal > T_upper) {
        cout << "get_T2e:: T is out of bound, "
             << "T = " << T << ", T_upper = " << T_upper*Util::hbarc
             << ", T_lower = " << T_lower*Util::hbarc << endl;
        exit(1);
    }
    if (T_goal < T_lower) return(eps_lower);

    double rel_accuracy = 1e-8;
    double abs_accuracy = 1e-15;
    double T_mid;
    int iter = 0;
    while (((eps_upper - eps_lower)/eps_mid > rel_accuracy
            && (eps_upper - eps_lower) > abs_accuracy) && iter < ntol) {
        T_mid = get_temperature(eps_mid, rhob);
        if (T_goal < T_mid)
            eps_upper = eps_mid;
        else 
            eps_lower = eps_mid;
        eps_mid = (eps_upper + eps_lower)/2.;
        iter++;
    }
    if (iter == ntol) {
        cout << "get_T2e_finite_rhob:: max iteration reached, "
             << "T = " << T << ", rhob = " << rhob << endl;;
        cout << "T_upper = " << T_upper*Util::hbarc
             << " , T_lower = " << T_lower*Util::hbarc << endl;
        cout << "eps_upper = " << eps_upper
             << " , eps_lower = " << eps_lower
             << ", diff = " << (eps_upper - eps_lower) << endl;
        exit(1);
    }
    return (eps_mid);
}

//! This function returns local energy density [1/fm^4] from
//! a given entropy density [1/fm^3] and rhob [1/fm^3]
//! using binary search
double EOS_base::get_s2e_finite_rhob(double s, double rhob) const {
    double eps_lower = 1e-15;
    double eps_upper = eps_max;
    double eps_mid   = (eps_upper + eps_lower)/2.;
    double s_lower   = get_entropy(eps_lower, rhob);
    double s_upper   = get_entropy(eps_upper, rhob);
    int ntol         = 1000;
    if (s < 0.0 || s > s_upper) {
        cout << "get_s2e_finite_rhob:: s is out of bound, "
             << "s = " << s << ", s_upper = " << s_upper
             << ", s_lower = " << s_lower << endl;
        exit(1);
    }
    if (s < s_lower) return(eps_lower);

    double rel_accuracy = 1e-8;
    double abs_accuracy = 1e-15;
    double s_mid;
    int iter = 0;
    while (((eps_upper - eps_lower)/eps_mid > rel_accuracy
            && (eps_upper - eps_lower) > abs_accuracy) && iter < ntol) {
        s_mid = get_entropy(eps_mid, rhob);
        if (s < s_mid)
            eps_upper = eps_mid;
        else 
            eps_lower = eps_mid;
        eps_mid = (eps_upper + eps_lower)/2.;
        iter++;
    }
    if (iter == ntol) {
        cout << "get_s2e_finite_rhob:: max iteration reached, "
             << "s = " << s << ", rhob = " << rhob << endl;
        cout << "s_upper = " << get_entropy(eps_upper, rhob)
             << " , s_lower = " << get_entropy(eps_lower, rhob) << endl;
        cout << "eps_upper = " << eps_upper
             << " , eps_lower = " << eps_lower
             << ", diff = " << (eps_upper - eps_lower) << endl;
        exit(1);
    }
    return (eps_mid);
}


std::string EOS_base::get_hydro_env_path() const {
    const char *EOSPATH = "HYDROPROGRAMPATH";
    char *pre_envPath = getenv(EOSPATH);
    std::string envPath;
    if (pre_envPath == 0) {
	    envPath=".";
    }
    else {
	    envPath=pre_envPath;
    }
    return(envPath);
}


void EOS_base::resize_table_info_arrays() {
    nb_bounds.resize(number_of_tables, 0.0);
    nb_spacing.resize(number_of_tables, 0.0);
    nb_length.resize(number_of_tables, 0);
    e_bounds.resize(number_of_tables, 0.0);
    e_spacing.resize(number_of_tables, 0.0);
    e_length.resize(number_of_tables, 0);
}


void EOS_base::check_eos_no_muB() const {
    // output EoS as function of e
    ostringstream file_name;
    file_name << "check_EoS_" << whichEOS << "_PST.dat";
    ofstream check_file(file_name.str().c_str());
    check_file << "#e(GeV/fm^3) P(GeV/fm^3) s(1/fm^3) T(GeV) cs^2" << endl;
    double e0 = 1e-3;
    double emax = 100;
    double de = 0.01;
    int ne = (emax - e0)/de + 1;
    for (int i = 0; i < ne; i++) {
        double e_local = (e0 + i*de)/hbarc;
        double p_local = get_pressure(e_local, 0.0);
        double s_local = get_entropy(e_local, 0.0);
        double T_local = get_temperature(e_local, 0.0);
        double cs2_local = get_cs2(e_local, 0.0);
        check_file << scientific << setw(18) << setprecision(8)
                   << e_local*hbarc << "   " << p_local*hbarc << "   " 
                   << s_local << "   " << T_local*hbarc << "   "
                   << cs2_local << endl;
    }
    check_file.close();
}


void EOS_base::check_eos_with_finite_muB() const {
    // output EoS as function of e for several rhob
    double rhob_pick[6] = {0.0, 0.02, 0.05, 0.1, 0.2, 0.5};
    for (int i = 0; i < 6; i++) {
        double rhob_local = rhob_pick[i];
        ostringstream file_name;
        file_name << "check_EoS_PST_rhob_" << rhob_pick[i] << ".dat";
        ofstream check_file(file_name.str().c_str());
        check_file << "#e(GeV/fm^3)  P(GeV/fm^3)  s(1/fm^3)  T(GeV)  cs^2  "
                   << "mu_B(GeV)  mu_S(GeV)  mu_C(GeV)" << endl;
        double e0 = 1e-3;
        double emax = 100;
        double de = 0.01;
        int ne = (emax - e0)/de + 1;
        for (int i = 0; i < ne; i++) {
            double e_local    = (e0 + i*de)/hbarc;
            double p_local    = get_pressure(e_local, rhob_local);
            double s_local    = get_entropy(e_local, rhob_local);
            double T_local    = get_temperature(e_local, rhob_local);
            double cs2_local  = get_cs2(e_local, rhob_local);
            double mu_b_local = get_muB(e_local, rhob_local);
            double mu_s_local = get_muS(e_local, rhob_local);
            double mu_c_local = get_muC(e_local, rhob_local);
            check_file << scientific << setw(18) << setprecision(8)
                       << e_local*hbarc << "   " << p_local*hbarc << "   " 
                       << s_local << "   " << T_local*hbarc << "   "
                       << cs2_local << "   " << mu_b_local*hbarc << "   "
                       << mu_s_local*hbarc << "   "
                       << mu_c_local*hbarc << endl;
        }
        check_file.close();
    }

    // output EoS as a function of rho_b for several energy density
    double e_pick[12] = {0.05, 0.1, 0.2, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7,
                         1.0, 3.0, 5.0};
    for (int i = 0; i < 12; i++) {
        double e_local = e_pick[i]/hbarc;
        ostringstream file_name;
        file_name << "check_EoS_PST_e_" << e_pick[i] << ".dat";
        ofstream check_file(file_name.str().c_str());
        check_file << "#rho_B(1/fm^3)  P(GeV/fm^3)  s(1/fm^3)  T(GeV)  cs^2  "
                   << "mu_B(GeV)  mu_S(GeV)  mu_C(GeV)" << endl;
        double rhob_0 = 0.0;
        double rhob_max = 1.0;
        double drhob = 0.01;
        int nrhob = (rhob_max - rhob_0)/drhob + 1;
        for (int i = 0; i < nrhob; i++) {
            double rhob_local = rhob_0 + i*drhob;
            double p_local    = get_pressure(e_local, rhob_local);
            double s_local    = get_entropy(e_local, rhob_local);
            double T_local    = get_temperature(e_local, rhob_local);
            double cs2_local  = get_cs2(e_local, rhob_local);
            double mu_b_local = get_muB(e_local, rhob_local);
            double mu_s_local = get_muS(e_local, rhob_local);
            double mu_c_local = get_muC(e_local, rhob_local);
            check_file << scientific << setw(18) << setprecision(8)
                       << rhob_local << "   " << p_local*hbarc << "   " 
                       << s_local << "   " << T_local*hbarc << "   "
                       << cs2_local << "   " << mu_b_local*hbarc << "   "
                       << mu_s_local*hbarc << "   "
                       << mu_c_local*hbarc << endl;
        }
        check_file.close();
    }

    // output EoS as a 2D function of e and rho_B
    string file_name1 = "check_EoS_pressure_2D.dat";
    string file_name2 = "check_EoS_cs2_2D.dat";
    ofstream check_file1(file_name1.c_str());
    ofstream check_file2(file_name2.c_str());
    double e_0 = 0.0;           // GeV/fm^3
    double e_max = 100.0;       // GeV/fm^3
    double de = 0.1;            // GeV/fm^3
    int ne = static_cast<int>((e_max - e_0)/de) + 1;
    double rhob_0 = 0.0;        // 1/fm^3
    double rhob_max = 1.0;      // 1/fm^3
    double drhob = 0.01;        // 1/fm^3
    int nrhob = static_cast<int>((rhob_max - rhob_0)/drhob) + 1;
    for (int i = 0; i < ne; i++) {
        double e_local = e_0 + i*de;
        for (int j = 0; j < nrhob; j++) {
            double rhob_local = rhob_0 + j*drhob;
            double p_local = get_pressure(e_local, rhob_local);
            double cs2_local = get_cs2(e_local, rhob_local);
            check_file1 << scientific << setw(18) << setprecision(8)
                        << p_local << "  ";
            check_file2 << scientific << setw(18) << setprecision(8)
                        << cs2_local << "  ";
        }
        check_file1 << endl;
        check_file2 << endl;
    }
    check_file1.close();
    check_file2.close();

    double sovernB[] = {10.0, 20.0, 30.0, 51.0, 70.0, 94.0, 144.0, 420.0};
    int array_length = sizeof(sovernB)/sizeof(double);
    double s_0 = 0.00;         // 1/fm^3
    double s_max = 100.0;      // 1/fm^3
    double ds = 0.005;         // 1/fm^3
    int ns = static_cast<int>((s_max - s_0)/ds) + 1;
    for (int i = 0; i < array_length; i++) {
        ostringstream file_name;
        file_name << "check_EoS_cs2_vs_e_sovernB_" << sovernB[i] << ".dat";
        ofstream check_file9(file_name.str().c_str());
        check_file9 << "# e(GeV/fm^3)  T(GeV)  cs^2  mu_B(GeV)  "
                    << "s(1/fm^3)  rho_B(1/fm^3)  dP/de  dP/drho  "
                    << "mu_S(GeV)  mu_C(GeV)" << endl;
        for (int j = 0; j < ns; j++) {
            double s_local     = s_0 + j*ds;
            double nB_local    = s_local/sovernB[i];
            double e_local     = get_s2e(s_local, nB_local);
            double s_check     = get_entropy(e_local, nB_local);
            double cs2_local   = get_cs2(e_local, nB_local);
            double dpde        = p_e_func(e_local, nB_local);
            double dpdrho      = p_rho_func(e_local, nB_local);
            double temperature = get_temperature(e_local, nB_local)*hbarc;
            double mu_B        = get_muB(e_local, nB_local)*hbarc;
            double mu_S        = get_muS(e_local, nB_local)*hbarc;
            double mu_C        = get_muC(e_local, nB_local)*hbarc;
            check_file9 << scientific << setw(18) << setprecision(8)
                        << e_local*hbarc << "  " << temperature << "  "
                        << cs2_local << "  " << mu_B << "  " 
                        << s_check << "  " << nB_local << "  "
                        << dpde << "  " << dpdrho << "  "
                        << mu_S << "  " << mu_C << endl;
        }
        check_file9.close();
    }
}
