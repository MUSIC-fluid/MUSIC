// Copyright @ Bjoern Schenke, Sangyong Jeon, Charles Gale, and Chun Shen
#ifndef SRC_MUSIC_H_
#define SRC_MUSIC_H_

#include <memory>

#include "util.h"
#include "cell.h"
#include "grid.h"
#include "data.h"
#include "eos.h"
#include "hydro_source_base.h"
#include "read_in_parameters.h"
#include "pretty_ostream.h"
#include "HydroinfoMUSIC.h"

//! This is a wrapper class for the MUSIC hydro
class MUSIC {
 private:
    //! records running mode
    int mode;

    //! flag to tell whether hydro is initialized
    int flag_hydro_initialized;

    //! flag to tell whether hydro is run
    int flag_hydro_run;

    int system_status_;

    InitData DATA;

    EOS eos;

    SCGrid arena_prev;
    SCGrid arena_current;
    SCGrid arena_future;

    std::shared_ptr<HydroSourceBase> hydro_source_terms_ptr;

    std::shared_ptr<HydroinfoMUSIC> hydro_info_ptr;

    pretty_ostream music_message;

 public:
    MUSIC(std::string input_file);
    ~MUSIC();

    //! this function returns the running mode
    int get_running_mode() {return(mode);}

    //! This function initialize hydro
    void initialize_hydro();

    //! This function change the parameter value in DATA
    void set_parameter(std::string parameter_name, double value);

    //! this is a shell function to run hydro
    int run_hydro();

    //! this is a shell function to run Cooper-Frye
    int run_Cooper_Frye();

    //! this function adds hydro source terms pointer
    void add_hydro_source_terms(
            std::shared_ptr<HydroSourceBase> hydro_source_ptr_in);

    //! This function setup source terms from dynamical initialization
    void generate_hydro_source_terms();

    //! This function calls routine to check EoS
    void check_eos();

    //! this is a test function to output the transport coefficients as
    //! function of T and mu_B
    void output_transport_coefficients();

    void clean_all_the_surface_files();

    void initialize_hydro_from_jetscape_preequilibrium_vectors(
        const double dx, const double dz, const double z_max, const int nz,
        std::vector<double> e_in, std::vector<double> P_in,
        std::vector<double> u_tau_in, std::vector<double> u_x_in,
        std::vector<double> u_y_in,   std::vector<double> u_eta_in,
        std::vector<double> pi_00_in, std::vector<double> pi_01_in,
        std::vector<double> pi_02_in, std::vector<double> pi_03_in,
        std::vector<double> pi_11_in, std::vector<double> pi_12_in,
        std::vector<double> pi_13_in, std::vector<double> pi_22_in,
        std::vector<double> pi_23_in, std::vector<double> pi_33_in,
        std::vector<double> Bulk_pi_in);

    void get_hydro_info(
        const double x, const double y, const double z, const double t,
        fluidCell* fluid_cell_info);
    int get_number_of_fluid_cells() const {
        return(hydro_info_ptr->get_number_of_fluid_cells());
    }
    void get_fluid_cell_with_index(const int idx, fluidCell *info) const {
        return(hydro_info_ptr->get_fluid_cell_with_index(idx, info));
    }
    void clear_hydro_info_from_memory();

    double get_hydro_tau0() const {return(hydro_info_ptr->get_hydro_tau0());}
    double get_hydro_dtau() const {return(hydro_info_ptr->get_hydro_dtau());}
    double get_hydro_tau_max() const {
        return(hydro_info_ptr->get_hydro_tau_max());
    }
    double get_hydro_dx() const {return(hydro_info_ptr->get_hydro_dx());}
    double get_hydro_x_max() const {return(hydro_info_ptr->get_hydro_x_max());}
    double get_hydro_deta() const {return(hydro_info_ptr->get_hydro_deta());}
    double get_hydro_eta_max() const {
        return(hydro_info_ptr->get_hydro_eta_max());
    }
    int get_ntau() const {return(hydro_info_ptr->get_ntau());}
    int get_neta() const {return(hydro_info_ptr->get_neta());}
    int get_nx()   const {return(hydro_info_ptr->get_nx()  );}
    bool is_boost_invariant() const {
        return(hydro_info_ptr->is_boost_invariant());
    }
};

#endif  // SRC_MUSIC_H_
