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
}


MUSIC::~MUSIC() {
    if (mode == 1 && mode == 2) {
        delete init;
        delete evolve;
    }
    delete eos;
    delete util;
}

int MUSIC::initialize_hydro() {
    // clean all the surface files
    system("rm surface.dat surface?.dat surface??.dat 2> /dev/null");

    init = new Init(eos, DATA);
    init->InitArena(DATA, &arena);
    return(0);
}


int MUSIC::run_hydro() {
    // this is a shell function to run hydro
    evolve = new Evolve(eos, DATA);
    evolve->EvolveIt(DATA, arena);
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
