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

    InitData DATA;

    EOS eos;

    SCGrid arena_prev;
    SCGrid arena_current;
    SCGrid arena_future;

    Init *init     = nullptr;
    Evolve *evolve = nullptr;
    Freeze *freeze = nullptr;

    hydro_source hydro_source_terms;

    pretty_ostream music_message;

 public:
    MUSIC(std::string input_file);
    ~MUSIC();

    //! this function returns the running mode
    int get_running_mode() {return(mode);}

    //! This function initialize hydro
    int initialize_hydro();

    //! this is a shell function to run hydro
    int run_hydro();

    //! this is a shell function to run Cooper-Frye
    int run_Cooper_Frye();

    void check_eos();
    //! this is a test function to output the transport coefficients as
    //! function of T and mu_B
    void output_transport_coefficients();
};

#endif  // SRC_MUSIC_H_
