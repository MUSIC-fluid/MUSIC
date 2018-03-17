// Original copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
// Massively cleaned up and improved by Chun Shen 2015-2016
#include <stdio.h>
#include <sys/stat.h>

#include <string>
#include "music.h"
#include "dissipative.h"

MUSIC::MUSIC(string input_file) : 
    DATA(ReadInParameters::read_in_parameters(input_file)),
    eos(DATA),
    hydro_source_terms(DATA) {
    mode = DATA.mode;
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
}


//! This function initialize hydro
int MUSIC::initialize_hydro() {
    // clean all the surface files
    int status = system(
                    "rm surface.dat surface?.dat surface??.dat 2> /dev/null");

    init = new Init(eos, DATA, hydro_source_terms);
    init->InitArena(arena_prev, arena_current, arena_future);
    flag_hydro_initialized = 1;
    return(status);
}


//! this is a shell function to run hydro
int MUSIC::run_hydro() {
    if (evolve != nullptr) {
        delete evolve;
    }

    evolve = new Evolve(eos, DATA, hydro_source_terms);

    evolve->EvolveIt(arena_prev, arena_current, arena_future);
        
    flag_hydro_run = 1;
    return(0);
}


//! this is a shell function to run Cooper-Frye
int MUSIC::run_Cooper_Frye() {
    if (freeze != nullptr) {
        delete freeze;
    }
    freeze = new Freeze(&DATA);
    freeze->CooperFrye_pseudo(DATA.particleSpectrumNumber, mode, &DATA, &eos);
    return(0);
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

