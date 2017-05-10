// Original copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
// Massively cleaned up and improved by Chun Shen 2015-2016
#include <stdio.h>
#include <sys/stat.h>

#include <cstring>
#include "./music.h"
#include "./dissipative.h"

using namespace std;

MUSIC::MUSIC(InitData *DATA_in, string input_file) {
    DATA = DATA_in;
    reader.read_in_parameters(DATA, input_file);
    mode = DATA->mode;
    eos = new EOS(DATA);
    util = new Util();
    if (DATA->Initial_profile == 12 || DATA->Initial_profile == 13) {
        hydro_source_ptr = new hydro_source(DATA_in);
    } else if (DATA->Initial_profile == 30) {
        hydro_source_ptr = new hydro_source(DATA_in);
    }
    flag_hydro_run = 0;
    flag_hydro_initialized = 0;
}


MUSIC::~MUSIC() {
    if (flag_hydro_initialized == 1) {
        delete init;
    }
    if (flag_hydro_run == 1) {
        delete evolve;
    }
    delete eos;
    delete util;
    if (DATA->Initial_profile == 12 || DATA->Initial_profile == 13) {
        delete hydro_source_ptr;
    } else if (DATA->Initial_profile == 30) {
        delete hydro_source_ptr;
    }
}

int MUSIC::initialize_hydro() {
    // clean all the surface files
    system("rm surface.dat surface?.dat surface??.dat 2> /dev/null");

    init = new Init(eos, DATA, hydro_source_ptr);
    init->InitArena(DATA, &arena);
    flag_hydro_initialized = 1;
    return(0);
}


int MUSIC::run_hydro() {
    // this is a shell function to run hydro
    evolve = new Evolve(eos, DATA, hydro_source_ptr);
    evolve->EvolveIt(DATA, arena);
    flag_hydro_run = 1;
    return(0);
}


void MUSIC::output_transport_coefficients() {
    // this is a test function to output the transport coefficients as
    // function of T and mu_B
    cout << "output transport coefficients as functions of T and mu_B" << endl;
    Diss *temp_dissipative_ptr = new Diss(eos, DATA);
    temp_dissipative_ptr->output_eta_over_s_T_and_muB_dependence(DATA);
    temp_dissipative_ptr->output_eta_over_s_along_const_sovernB(DATA);
    temp_dissipative_ptr->output_kappa_T_and_muB_dependence(DATA);
    temp_dissipative_ptr->output_kappa_along_const_sovernB(DATA);
    delete temp_dissipative_ptr;
}
