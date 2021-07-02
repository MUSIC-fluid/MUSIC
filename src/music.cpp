// Original copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
// Massively cleaned up and improved by Chun Shen 2015-2016
#include <stdio.h>
#include <sys/stat.h>
#include <vector>
#include <memory>

#include <string>
#include "music.h"
#include "init.h"
#include "evolve.h"
#include "dissipative.h"
#include "data_struct.h"
#include "hydro_source_strings.h"
#include "hydro_source_ampt.h"
#include "hydro_source_TATB.h"

#ifdef GSL
    #include "freeze.h"
#endif

using std::vector;

MUSIC::MUSIC(std::string input_file) :
    DATA(ReadInParameters::read_in_parameters(input_file)),
    eos(DATA.whichEOS) {

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
    generate_hydro_source_terms();
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
    } else if (DATA.Initial_profile == 112) {  // source from TA and TB
        auto hydro_source_ptr = std::shared_ptr<HydroSourceTATB> (
                                            new HydroSourceTATB (DATA));
        add_hydro_source_terms(hydro_source_ptr);
    }
}


void MUSIC::clean_all_the_surface_files() {
    system_status_ = system(
            "rm surface.dat surface?.dat surface??.dat 2> /dev/null");
}


//! This function change the parameter value in DATA
void MUSIC::set_parameter(std::string parameter_name, double value) {
    ReadInParameters::set_parameter(DATA, parameter_name, value);
}


//! This function initialize hydro
void MUSIC::initialize_hydro() {
    clean_all_the_surface_files();

    Init initialization(eos, DATA, hydro_source_terms_ptr);
    initialization.InitArena(arena_prev, arena_current, arena_future);
    flag_hydro_initialized = 1;
}


//! this is a shell function to run hydro
int MUSIC::run_hydro() {
    Evolve evolve_local(eos, DATA, hydro_source_terms_ptr);

    if (hydro_info_ptr == nullptr && DATA.store_hydro_info_in_memory == 1) {
        hydro_info_ptr = std::make_shared<HydroinfoMUSIC> ();
    }
    evolve_local.EvolveIt(arena_prev, arena_current, arena_future,
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

    if (nz > 1) {
        DATA.boost_invariant = false;
        DATA.delta_eta       = dz;
        DATA.neta            = nz;
        DATA.eta_size        = nz*dz;
    }
    DATA.delta_x = dx;
    DATA.delta_y = dx;

    Init initialization(eos, DATA, hydro_source_terms_ptr);
    initialization.get_jetscape_preequilibrium_vectors(
        e_in, P_in, u_tau_in, u_x_in, u_y_in, u_eta_in,
        pi_00_in, pi_01_in, pi_02_in, pi_03_in, pi_11_in, pi_12_in, pi_13_in,
        pi_22_in, pi_23_in, pi_33_in, Bulk_pi_in);
    initialization.InitArena(arena_prev, arena_current, arena_future);
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
