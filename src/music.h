// Copyright @ Bjoern Schenke, Sangyong Jeon, Charles Gale, and Chun Shen
#ifndef SRC_MUSIC_H_
#define SRC_MUSIC_H_

#include "util.h"
#include "cell.h"
#include "grid.h"
#include "data.h"
#include "freeze.h"
#include "init.h"
#include "eos.h"
#include "evolve.h"
#include "hydro_source.h"
#include "read_in_parameters.h"
#include "pretty_ostream.h"

//! This is a wrapper class for the MUSIC hydro
class MUSIC {
 private:
    //! records running mode
    int mode;

    //! flag to tell whether hydro is initialized
    int flag_hydro_initialized;
    
    //! flag to tell whether hydro is run
    int flag_hydro_run;

    InitData *DATA = nullptr;

    EOS *eos = nullptr;

    Grid arena;

    Init *init     = nullptr;
    Evolve *evolve = nullptr;
    Freeze *freeze = nullptr;

    hydro_source *hydro_source_ptr;

    pretty_ostream music_message;

 public:
    MUSIC(InitData *DATA_in, std::string input_file);
    ~MUSIC();

    //! this function returns the running mode
    int get_running_mode() {return(mode);}

    //! This function initialize hydro
    int initialize_hydro();

    //! this is a shell function to run hydro
    int run_hydro();

    //! this is a shell function to run Cooper-Frye
    int run_Cooper_Frye();

    //! this is a test function to output the transport coefficients as
    //! function of T and mu_B
    void output_transport_coefficients();

    //! This function prints out the welcome message
    void welcome_message();

    //! This function prints out code desciprtion and copyright information
    void display_code_description_and_copyright();

    //! This function prints out the program logo
    void display_logo(int selector);
};

#endif  // SRC_MUSIC_H_
