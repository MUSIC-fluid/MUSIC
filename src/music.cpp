// Original copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
// Massively cleaned up and improved by Chun Shen 2015-2016
#include <stdio.h>
#include <sys/stat.h>
#include <vector>
#include <memory>
#include <iomanip>
#include <string>
#include "music.h"
#include "init.h"
#include "evolve.h"
#include "dissipative.h"
#include "data_struct.h"
#include "hydro_source_strings.h"
#include "hydro_source_ampt.h"
#include "hydro_source_smash.h"
#include "hydro_source_TATB.h"

#ifdef GSL
    #include "freeze.h"
#endif

using std::vector;
using std::setw;
using std::setprecision;
using std::scientific;


MUSIC::MUSIC(std::string input_file) :
    DATA(ReadInParameters::read_in_parameters(input_file)),
    eos(DATA.whichEOS) {

    DATA.reRunHydro = false;
    DATA.reRunCount = 0;
    mode                   = DATA.mode;
    flag_hydro_run         = 0;
    flag_hydro_initialized = 0;

    // setup hydro evolution information
    hydro_info_ptr         = nullptr;
    if (DATA.store_hydro_info_in_memory == 1) {
        hydro_info_ptr = std::make_shared<HydroinfoMUSIC> ();
    }

    // setup source terms
    hydro_source_terms_ptr = nullptr;
}


MUSIC::~MUSIC() {
}


//! This function adds hydro source terms pointer
void MUSIC::add_hydro_source_terms(
            std::shared_ptr<HydroSourceBase> hydro_source_ptr_in) {
    hydro_source_terms_ptr = hydro_source_ptr_in;
}


//! This function setup source terms from dynamical initialization
void MUSIC::generate_hydro_source_terms() {
    if (DATA.Initial_profile == 13 || DATA.Initial_profile == 131) {
        // MC-Glauber-LEXUS
        auto hydro_source_ptr = std::shared_ptr<HydroSourceStrings> (
                                            new HydroSourceStrings (DATA));
        add_hydro_source_terms(hydro_source_ptr);
    } else if (DATA.Initial_profile == 30) {  // AMPT
        auto hydro_source_ptr = std::shared_ptr<HydroSourceAMPT> (
                                            new HydroSourceAMPT (DATA));
        add_hydro_source_terms(hydro_source_ptr);
    } else if (DATA.Initial_profile == 31){ // SMASH
        auto hydro_source_ptr = std::shared_ptr<HydroSourceSMASH> (
                                            new HydroSourceSMASH (DATA));
        add_hydro_source_terms(hydro_source_ptr);
    } else if (DATA.Initial_profile == 112 || DATA.Initial_profile == 113) {
        // source from TA and TB
        auto hydro_source_ptr = std::shared_ptr<HydroSourceTATB> (
                                            new HydroSourceTATB (DATA));
        add_hydro_source_terms(hydro_source_ptr);
    } 
}


void MUSIC::clean_all_the_surface_files() {
    system_status_ = system("rm surface*.dat 2> /dev/null");
}


//! This function change the parameter value in DATA
void MUSIC::set_parameter(std::string parameter_name, double value) {
    ReadInParameters::set_parameter(DATA, parameter_name, value);
}


//! This function initialize hydro
void MUSIC::initialize_hydro() {
    clean_all_the_surface_files();

    generate_hydro_source_terms();

    Init initialization(eos, DATA, hydro_source_terms_ptr);
    initialization.InitArena(arenaFieldsPrev, arenaFieldsCurr, arenaFieldsNext);
    flag_hydro_initialized = 1;
}


//! this is a shell function to run hydro
int MUSIC::run_hydro() {
    Evolve evolve_local(eos, DATA, hydro_source_terms_ptr);

    if (hydro_info_ptr == nullptr && DATA.store_hydro_info_in_memory == 1) {
        hydro_info_ptr = std::make_shared<HydroinfoMUSIC> ();
    }
    evolve_local.EvolveIt(arenaFieldsPrev, arenaFieldsCurr, arenaFieldsNext,
                          (*hydro_info_ptr));
    flag_hydro_run = 1;
    return(0);
}


//! this is a shell function to run Cooper-Frye
int MUSIC::run_Cooper_Frye() {
#ifdef GSL
    Freeze cooper_frye(&DATA);
    cooper_frye.CooperFrye_pseudo(DATA.particleSpectrumNumber, mode,
                                  &DATA, &eos);
#endif
    return(0);
}


void MUSIC::check_eos() {
    music_message << "check eos ...";
    music_message.flush("info");
    eos.check_eos();
}

void MUSIC::check_source() {
    generate_hydro_source_terms();

    double tau0 = hydro_source_terms_ptr->get_source_tau_min();
    //double tau0 = 0.5;
    //double tau0 = hydro_source_terms_ptr->get_source_tau_max();
    double nx = DATA.nx;
    double ny = DATA.ny;
    double neta = DATA.neta;

    FlowVec u                   = {0};
    u[0]                        = 1.0;
    EnergyFlowVec j_mu          = {0};
    EnergyFlowVec j_mu_one_step = {0};

    std::ostringstream file_name;
    file_name << "check_hydro_source_tau0.dat";
    std::ofstream check_file(file_name.str().c_str());
    check_file << "# tau(fm) x(fm) y(fm) eta ";
    check_file << "e(GeV/fm^3) j_1 j_2 j_3 ";
    check_file << "rhoB(1/fm^3)  rhoQ(1/fm^3)  rhoS(1/fm^3)  " << std::endl;

    hydro_source_terms_ptr->prepare_list_for_current_tau_frame(tau0);

    double tau_i = 0.0;
    double dtau = 0.005;

    // start with tau0
    int n_tau_steps = static_cast<int>((tau0 - tau_i)/dtau);

    for (int ieta = 0; ieta < neta; ieta++){
        double eta = (DATA.delta_eta)*ieta - (DATA.eta_size)/2.0;
        //double eta = 0;
        for (int ix = 0; ix < nx; ix++) {
            double x_local = - DATA.x_size/2. + ix*DATA.delta_x;
            for (int iy = 0; iy < ny; iy++) {
                double y_local = - DATA.y_size/2. + iy*DATA.delta_y;
                double epsilon = 0.0;double j_1 = 0.0;
                double j_2 = 0.0;double j_3 = 0.0;
                double rhob = 0.0;double rhoq = 0.0;double rhos = 0.0;
                double rhob_local = 0.0;double rhoq_local = 0.0;double rhos_local = 0.0;
                //hydro_source_terms_ptr->get_hydro_energy_source(
                //                   tau0, x_local, y_local, eta, u, j_mu);
                //epsilon = j_mu[0];
                //j_1 = j_mu[1];
                //j_2 = j_mu[2];
                //j_3 = j_mu[3];

                j_mu = {0.0};

                for (int i = 0; i < n_tau_steps; i++) {
                    double tau_local = tau_i + (i + 0.5)*dtau;

                    j_mu_one_step = {0};
                    hydro_source_terms_ptr->get_hydro_energy_source(tau_local, x_local, y_local, eta, u, j_mu_one_step);
                    for (int j = 0; j < 4; j++) {
                        j_mu[j] += tau_local*j_mu_one_step[j]*dtau;
                    }

                    if (DATA.turn_on_rhob == 1) {
                        rhob_local = hydro_source_terms_ptr->get_hydro_rhob_source(
                                tau_local, x_local, y_local, eta, u);
                    } else {
                        rhob_local = 0.0;
                    }
                    if (DATA.turn_on_QS == 1) {
                        rhoq_local = hydro_source_terms_ptr->get_hydro_rhoq_source(
                                tau_local, x_local, y_local, eta, u);
                        rhos_local = hydro_source_terms_ptr->get_hydro_rhos_source(
                                tau_local, x_local, y_local, eta, u);
                    } else {
                        rhoq_local = 0.0;
                        rhos_local = 0.0;
                    }
                    rhob += tau_local*rhob_local*dtau;
                    rhoq += tau_local*rhoq_local*dtau;
                    rhos += tau_local*rhos_local*dtau;
                }
                for (int j = 0; j < 4; j++) {
                    j_mu[j] /= tau0;
                    j_mu[j] *= DATA.delta_tau;
                }
                epsilon = j_mu[0];
                j_1 = j_mu[1];
                j_2 = j_mu[2];
                j_3 = j_mu[3];
                rhob /= tau0;
                rhoq /= tau0;
                rhos /= tau0;

                rhob *= DATA.delta_tau;
                rhoq *= DATA.delta_tau;
                rhos *= DATA.delta_tau;

                check_file << scientific << setw(18) << std::setprecision(8)
                    << tau0 << "   " << x_local << "   " << y_local << " " << eta << " "
                    <<  epsilon << "   " << j_1 << "   " << j_2 << "  " << j_3
                    << "  " << rhob << "  " << rhoq << "  " << rhos << std::endl;
            }

        }

    }
    check_file.close();
}








//! this is a test function to output the transport coefficients as
//! function of T and mu_B
void MUSIC::output_transport_coefficients() {
    music_message << "output transport coefficients as functions of T and muB";
    music_message.flush("info");
    Diss temp_dissipative_ptr(eos, DATA);
    temp_dissipative_ptr.output_eta_over_s_T_and_muB_dependence();
    temp_dissipative_ptr.output_eta_over_s_along_const_sovernB();
    temp_dissipative_ptr.output_kappa_T_and_muB_dependence();
    temp_dissipative_ptr.output_kappa_along_const_sovernB();
    temp_dissipative_ptr.output_zeta_over_s_T_and_muB_dependence();
    temp_dissipative_ptr.output_zeta_over_s_along_const_sovernB();
}


void MUSIC::initialize_hydro_from_jetscape_preequilibrium_vectors(
        const double dx, const double dz, const double z_max, const int nz,
        vector<double> e_in, vector<double> P_in,
        vector<double> u_tau_in, vector<double> u_x_in,
        vector<double> u_y_in,   vector<double> u_eta_in,
        vector<double> pi_00_in, vector<double> pi_01_in,
        vector<double> pi_02_in, vector<double> pi_03_in,
        vector<double> pi_11_in, vector<double> pi_12_in,
        vector<double> pi_13_in, vector<double> pi_22_in,
        vector<double> pi_23_in, vector<double> pi_33_in,
        vector<double> Bulk_pi_in) {

    DATA.Initial_profile = 42;
    clean_all_the_surface_files();

    DATA.neta = nz;
    if (nz > 1) {
        DATA.boost_invariant = false;
        DATA.delta_eta = dz;
        DATA.eta_size = nz*dz;
    } else {
        DATA.boost_invariant = true;
        DATA.delta_eta = 0.1;
        DATA.eta_size = 0.;
    }
    DATA.delta_x = dx;
    DATA.delta_y = dx;

    Init initialization(eos, DATA, hydro_source_terms_ptr);
    initialization.get_jetscape_preequilibrium_vectors(
        e_in, P_in, u_tau_in, u_x_in, u_y_in, u_eta_in,
        pi_00_in, pi_01_in, pi_02_in, pi_03_in, pi_11_in, pi_12_in, pi_13_in,
        pi_22_in, pi_23_in, pi_33_in, Bulk_pi_in);
    initialization.InitArena(arenaFieldsPrev, arenaFieldsCurr, arenaFieldsNext);
    flag_hydro_initialized = 1;
}

void MUSIC::get_hydro_info(
        const double x, const double y, const double z, const double t,
        fluidCell* fluid_cell_info) {
    if (DATA.store_hydro_info_in_memory == 0 || hydro_info_ptr == nullptr) {
        music_message << "hydro evolution informaiton is not stored "
                      << "in the momeory! Please set the parameter "
                      << "store_hydro_info_in_memory to 1~";
        music_message.flush("error");
        exit(1);
    }
    hydro_info_ptr->getHydroValues(x, y, z, t, fluid_cell_info);
}

void MUSIC::clear_hydro_info_from_memory() {
    if (DATA.store_hydro_info_in_memory == 0 || hydro_info_ptr == nullptr) {
        music_message << "The parameter store_hydro_info_in_memory is 0. "
                      << "No need to clean memory~";
        music_message.flush("warning");
    } else {
        hydro_info_ptr->clean_hydro_event();
    }
}
