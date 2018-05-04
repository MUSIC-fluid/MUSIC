// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

#include "util.h"
#include "eos.h"
#include "data.h"
using namespace std;

#define ideal_cs2 (1.0/3.0)

EOS::EOS(const InitData &para_in) : parameters_ptr(para_in)  {
    whichEOS = parameters_ptr.whichEOS;
    number_of_tables = 0;
    eps_max = 1e5;  // [1/fm^4]
    initialize_eos();
}

// destructor
EOS::~EOS() {
    for (int itable = 0; itable < number_of_tables; itable++) {
        Util::mtx_free(pressure_tb[itable],
                       nb_length[itable], e_length[itable]);
        Util::mtx_free(temperature_tb[itable],
                       nb_length[itable], e_length[itable]);
        if (parameters_ptr.whichEOS == 1 || parameters_ptr.whichEOS > 9) {
            Util::mtx_free(mu_B_tb[itable],
                           nb_length[itable], e_length[itable]);
        }
        if (parameters_ptr.whichEOS == 11) {
            Util::mtx_free(mu_S_tb[itable],
                           nb_length[itable], e_length[itable]);
        }
    }
    delete [] pressure_tb;
    delete [] temperature_tb;
    if (parameters_ptr.whichEOS == 1 || parameters_ptr.whichEOS > 9) {
        delete [] mu_B_tb;
    }
    if (parameters_ptr.whichEOS == 11) {
        delete [] mu_S_tb;
    }
}


void EOS::initialize_eos() {
    if (parameters_ptr.Initial_profile == 0) {
        music_message.info("Using the ideal gas EOS");
        number_of_tables = 0;
    } else if (parameters_ptr.whichEOS == 1) {
        music_message.info("Using EOS-Q from AZHYDRO");
        init_eos();
    } else if (parameters_ptr.whichEOS == 2) {
        music_message.info("Using lattice EOS from Huovinen/Petreczky");
        init_eos_s95p(0);
    } else if (parameters_ptr.whichEOS == 3) {
        music_message << "Using lattice EOS from Huovinen/Petreczky with "
                      << "partial chemical equilibrium (PCE) "
                      << "chem. f.o. at 150 MeV";
        music_message.flush("info");
        init_eos_s95p(1);
    } else if (parameters_ptr.whichEOS == 4) {
        music_message << "Using lattice EOS from Huovinen/Petreczky with "
                      << "partial chemical equilibrium (PCE) "
                      << "chem. f.o. at 155 MeV";
        music_message.flush("info");
        init_eos_s95p(2);
    } else if (parameters_ptr.whichEOS == 5) {
        music_message << "Using lattice EOS from Huovinen/Petreczky with "
                      << "partial chemical equilibrium (PCE) "
                      << "chem. f.o. at 160 MeV";
        music_message.flush("info");
        init_eos_s95p(3);
    } else if (parameters_ptr.whichEOS == 6) {
        music_message << "Using lattice EOS from Huovinen/Petreczky with "
                      << "partial chemical equilibrium (PCE) chem. f.o. "
                       << "at 165 MeV";
        music_message.flush("info");
        init_eos_s95p(4);
    } else if (parameters_ptr.whichEOS == 7) {
        music_message.info(
            "Using lattice EOS from Huovinen/Petreczky s95p-v1.2 (for UrQMD)");
        init_eos_s95p(5);
    } else if (parameters_ptr.whichEOS == 8) {
        music_message.info("Using lattice EOS parameterization from WB");
    } else if (parameters_ptr.whichEOS == 10) {
        music_message.info("Using lattice EOS from A. Monnai");
        init_eos10();
    } else if (parameters_ptr.whichEOS == 11) {
        music_message.info("Using lattice EOS from Pasi");
        init_eos11();
    } else if (parameters_ptr.whichEOS == 12) {
        music_message.info("Using lattice EOS from A. Monnai (up to mu_B^6)");
        init_eos12();
    } else {
        music_message << "No EOS for whichEOS = " << parameters_ptr.whichEOS
             << ". Use EOS_to_use = 0 (ideal gas) 1 (AZHYDRO EOS-Q), "
             << "2 (s95p-v1), 3 (s95p-PCE150-v1), 4 (s95p-PCE155-v1), "
             << "5 (s95p-PCE160-v1), 6 (s95p-PCE165-v1),"
             << "7 (s95p-v1.2), "
             << "8 (WB), "
             << "10(lattice EOS at finite muB), "
             << "11(lattice EoS at finite muB from Pasi), "
             << "12(lattice EOS at finite muB from A. Monnai up to mu_B^6)";
        music_message.flush("error");
        exit(1);
    }

    if (whichEOS >= 2 && whichEOS < 8) {
        eps_max = (e_bounds[6] + e_spacing[6]*e_length[6])/hbarc;  // [1/fm^4]
    } else if (whichEOS == 10) {
        eps_max = (e_bounds[6] + e_spacing[6]*e_length[6])/hbarc;  // [1/fm^4]
    } else if (whichEOS == 11) {
        eps_max = (e_bounds[3] + e_spacing[3]*e_length[3])/hbarc;  // [1/fm^4]
    } else if (whichEOS == 12) {
        eps_max = (e_bounds[6] + e_spacing[6]*e_length[6])/hbarc;  // [1/fm^4]
    }
}


void EOS::init_eos() {
// read the azhydro pressure, temperature, and 
// baryon chemical potential from file
    whichEOS = 1;
    music_message.info("reading EOS...");
    auto envPath = get_hydro_env_path();
    music_message << "from path " << envPath.c_str() << "/EOS";
    music_message.flush("info");
    
    number_of_tables = 2;
    resize_table_info_arrays();

    string eos_file_string_array[2] = {"1", "2"};
    pressure_tb    = new double** [number_of_tables];
    temperature_tb = new double** [number_of_tables];
    mu_B_tb        = new double** [number_of_tables];
    for (int itable = 0; itable < number_of_tables; itable++) {
        std::ifstream eos_p(envPath + "/EOS/EOS-Q/aa"
                            + eos_file_string_array[itable] + "_p.dat");
        std::ifstream eos_T(envPath + "/EOS/EOS-Q/aa"
                            + eos_file_string_array[itable] + "_t.dat");
        std::ifstream eos_mub(envPath + "/EOS/EOS-Q/aa"
                              + eos_file_string_array[itable] + "_mb.dat");
  
        // read the first two lines:
        // first value of rhob, first value of epsilon
        // deltaRhob, deltaEpsilon, number of rhob steps, number of epsilon steps
        int N_e, N_rhob;
        eos_p >> nb_bounds[itable] >> e_bounds[itable];
        eos_p >> nb_spacing[itable] >> e_spacing[itable]
              >> N_rhob >> N_e;
        nb_length[itable] = N_rhob + 1;
        e_length[itable]  = N_e + 1;

        // skip the header in T and mu_B files
        string dummy;
        std::getline(eos_T, dummy);
        std::getline(eos_T, dummy);
        std::getline(eos_mub, dummy);
        std::getline(eos_mub, dummy);
 
        // allocate memory for pressure arrays
        pressure_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        temperature_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);
        mu_B_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                           e_length[itable]);

        // read pressure, temperature and chemical potential values
        for (int j = 0; j < e_length[itable]; j++) {
            for (int i = 0; i < nb_length[itable]; i++) {
                eos_p >> pressure_tb[itable][i][j];
                eos_T >> temperature_tb[itable][i][j];
                eos_mub >> mu_B_tb[itable][i][j];
            }
        }
    }
    music_message.info("Done reading EOS.");
}


void EOS::init_eos_s95p(int selector) {
    // read the lattice EOS pressure, temperature, and 
    // baryon chemical potential from file
    music_message.info("reading EOS s95p ...");

    auto envPath = get_hydro_env_path();
    stringstream spath;
    spath << envPath;
    if (selector == 0) {
        spath << "/EOS/s95p-v1/";
    } else if (selector == 1) {
        spath << "/EOS/s95p-PCE-v1/";
    } else if (selector == 2) {
        spath << "/EOS/s95p-PCE155/";
    } else if (selector == 3) {
        spath << "/EOS/s95p-PCE160/";
    } else if (selector == 4) {
        spath << "/EOS/s95p-PCE165-v0/";
    } else if (selector == 5) {
        spath << "/EOS/s95p-v1.2/";
    }
    
    music_message << "from path " << spath.str();
    music_message.flush("info");

    if (selector == 0) {
        spath << "s95p-v1_";
    } else if (selector == 5) {
        spath << "s95p-v1.2_";
    }
    
    number_of_tables = 7;
    resize_table_info_arrays();
    
    string eos_file_string_array[7] = {"1", "2", "3", "4", "5", "6", "7"};
    pressure_tb    = new double** [number_of_tables];
    temperature_tb = new double** [number_of_tables];
    for (int itable = 0; itable < number_of_tables; itable++) {
        std::ifstream eos_d(spath.str() + "dens"
                            + eos_file_string_array[itable] + ".dat");
        std::ifstream eos_T(spath.str() + "par"
                            + eos_file_string_array[itable] + ".dat");
        // read the first two lines with general info:
        // lowest value of epsilon
        // deltaEpsilon, number of epsilon steps (i.e. # of lines)
        eos_d >> e_bounds[itable];
        eos_d >> e_spacing[itable] >> e_length[itable];
        
        // skip the header in T file
        string dummy;
        std::getline(eos_T, dummy);
        std::getline(eos_T, dummy);
        
        // no rho_b dependence at the moment
        nb_length[itable] = 1;
        
        // allocate memory for pressure arrays
        pressure_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        temperature_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);

        // read pressure, temperature and chemical potential values
        // files have it backwards, so I start with maximum j and count down
        int i = 0;
        double d_dummy;
        for (int j = e_length[itable] - 1; j >= 0; j--) {
            eos_d >> d_dummy;
            eos_d >> pressure_tb[itable][i][j];
            eos_d >> d_dummy >> dummy >> dummy;;
            eos_T >> temperature_tb[itable][i][j] >> dummy >> dummy;
        }
    }
    music_message.info("Done reading EOS.");
}


void EOS::init_eos10() {
    // read the lattice EOS at finite muB
    // pressure, temperature, and baryon chemical potential from file
    music_message.info("reading EOS...");

    stringstream slocalpath;
    slocalpath << "./EOS/neos_2/";

    string path = slocalpath.str();
    music_message << "from path " << path;
    music_message.flush("info");
    
    number_of_tables = 7;
    resize_table_info_arrays();

    string eos_file_string_array[7] = {"0a", "0b", "0c", "1a", "2", "3", "4"};
    pressure_tb    = new double** [number_of_tables];
    temperature_tb = new double** [number_of_tables];
    mu_B_tb        = new double** [number_of_tables];

    for (int itable = 0; itable < number_of_tables; itable++) {
        std::ifstream eos_p(path + "neos" + eos_file_string_array[itable]
                            + "_p.dat");
        std::ifstream eos_T(path + "neos" + eos_file_string_array[itable]
                            + "_t.dat");
        std::ifstream eos_mub(path + "neos" + eos_file_string_array[itable]
                              + "_mb.dat");
        // read the first two lines with general info:
        // first value of rhob, first value of epsilon
        // deltaRhob, deltaE, number of rhob points, number of epsilon points
        // the table size is
        // (number of rhob points + 1, number of epsilon points + 1)
        int N_e, N_rhob;
        eos_p >> nb_bounds[itable] >> e_bounds[itable];
        eos_p >> nb_spacing[itable] >> e_spacing[itable]
              >> N_rhob >> N_e;
        nb_length[itable] = N_rhob + 1;
        e_length[itable]  = N_e + 1;

        // skip the header in T and mu_B files
        string dummy;
        std::getline(eos_T, dummy);
        std::getline(eos_T, dummy);
        std::getline(eos_mub, dummy);
        std::getline(eos_mub, dummy);

        // allocate memory for EOS arrays
        pressure_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        temperature_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);
        mu_B_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                           e_length[itable]);

        // read pressure, temperature and chemical potential values
        for (int j = 0; j < e_length[itable]; j++) {
            for (int i = 0; i < nb_length[itable]; i++) {
                eos_p >> pressure_tb[itable][i][j];
                eos_T >> temperature_tb[itable][i][j];
                eos_mub >> mu_B_tb[itable][i][j];
            }
        }
    }
    music_message.info("Done reading EOS.");
}


void EOS::init_eos11() {
    // read Pasi's lattice EOS at finite muB
    // pressure, temperature, and baryon chemical potential from file
    music_message.info("reading EOS (Pasi) at finite mu_B ...");

    stringstream slocalpath;
    slocalpath << "./EOS/s95p-finite_muB/";
    string path = slocalpath.str();

    music_message << "from path " << path;
    music_message.flush("info");
    
    number_of_tables = 4;
    resize_table_info_arrays();
    
    string eos_file_string_array[4] = {"1", "2", "3", "4"};
    pressure_tb    = new double** [number_of_tables];
    temperature_tb = new double** [number_of_tables];
    mu_B_tb        = new double** [number_of_tables];
    mu_S_tb        = new double** [number_of_tables];
    
    for (int itable = 0; itable < number_of_tables; itable++) {
        ifstream eos_p(path + "p" + eos_file_string_array[itable] + ".dat");
        ifstream eos_T(path + "t" + eos_file_string_array[itable] + ".dat");
        ifstream eos_mub(path + "mb" + eos_file_string_array[itable] + ".dat");
        ifstream eos_mus(path + "ms" + eos_file_string_array[itable] + ".dat");

        // read the first two lines with general info:
        // first value of rhob, first value of epsilon
        // deltaRhob, deltaE, number of rhob points, number of epsilon points
        // the table size is
        // (number of rhob points + 1, number of epsilon points + 1)
        int N_e, N_rhob;
        eos_p >> nb_bounds[itable] >> e_bounds[itable];
        eos_p >> nb_spacing[itable] >> e_spacing[itable]
              >> N_rhob >> N_e;

        nb_length[itable] = N_rhob;
        e_length[itable]  = N_e + 1;

        // skip the header in T and mu_B files
        string dummy;
        std::getline(eos_T, dummy);
        std::getline(eos_T, dummy);
        std::getline(eos_mub, dummy);
        std::getline(eos_mub, dummy);
        std::getline(eos_mus, dummy);
        std::getline(eos_mus, dummy);

        // allocate memory for pressure arrays
        pressure_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        temperature_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);
        mu_B_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                           e_length[itable]);
        mu_S_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                           e_length[itable]);

        // read pressure, temperature and chemical potential values
        for (int j = N_e; j >= 0; j--) {
            for (int i = 0; i < N_rhob; i++) {
                eos_p >> pressure_tb[itable][i][j];
                eos_T >> temperature_tb[itable][i][j];
                eos_mub >> mu_B_tb[itable][i][j];
                eos_mus >> mu_S_tb[itable][i][j];
            }
        }
    }
    music_message.info("Done reading EOS.");
}

void EOS::init_eos12() {
    // read the lattice EOS at finite muB
    // pressure, temperature, and baryon chemical potential from file
    music_message.info("reading EOS ...");
    
    stringstream slocalpath;
    slocalpath << "./EOS/neos_3/";

    string path = slocalpath.str();
    music_message << "from path " << path;
    music_message.flush("info");
    
    number_of_tables = 7;
    resize_table_info_arrays();

    string eos_file_string_array[7] = {"1", "2", "3", "4", "5", "6", "7"};
    pressure_tb    = new double** [number_of_tables];
    temperature_tb = new double** [number_of_tables];
    mu_B_tb        = new double** [number_of_tables];

    for (int itable = 0; itable < number_of_tables; itable++) {
        std::ifstream eos_p(path + "neos"
                            + eos_file_string_array[itable] + "_p.dat");
        std::ifstream eos_T(path + "neos"
                            + eos_file_string_array[itable] + "_t.dat");
        std::ifstream eos_muB(path + "neos"
                              + eos_file_string_array[itable] + "_mub.dat");

        // read the first two lines with general info:
        // first value of rhob, first value of epsilon
        // deltaRhob, deltaE, number of rhob points, number of epsilon points
        // the table size is
        // (number of rhob points + 1, number of epsilon points + 1)
        int N_e, N_rhob;
        eos_p >> nb_bounds[itable] >> e_bounds[itable];
        eos_p >> nb_spacing[itable] >> e_spacing[itable]
              >> N_rhob >> N_e;
        nb_length[itable] = N_rhob + 1;
        e_length[itable]  = N_e + 1;
        
        // skip the header in T and mu_B files
        string dummy;
        std::getline(eos_T, dummy);
        std::getline(eos_T, dummy);
        std::getline(eos_muB, dummy);
        std::getline(eos_muB, dummy);
        
        // allocate memory for EOS arrays
        pressure_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        temperature_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);
        mu_B_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                           e_length[itable]);

        // read pressure, temperature and chemical potential values
        for (int j = 0; j < e_length[itable]; j++) {
            for (int i = 0; i < nb_length[itable]; i++) {
                eos_p >> pressure_tb[itable][i][j];
                eos_T >> temperature_tb[itable][i][j];
                eos_muB >> mu_B_tb[itable][i][j];
            }
        }
    }
    music_message.info("Done reading EOS.");
}


double EOS::get_dpOverde3(double e, double rhob) const {
   double eLeft = 0.9*e;
   double eRight = 1.1*e;

   double pL = get_pressure(eLeft, rhob);   // 1/fm^4
   double pR = get_pressure(eRight, rhob);  // 1/fm^4
      
   double dpde = (pR - pL)/(eRight - eLeft);
   return dpde;
}


double EOS::get_dpOverdrhob2(double e, double rhob) const {
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


double EOS::get_cs2(double e, double rhob) const {
    double f = calculate_velocity_of_sound_sq(e, rhob);
    return(f);
}


//! this function output the EoS matrix on the grid to files for checking
//! purpose
void EOS::output_eos_matrix(int ne, int nrhob, double** matrix_ptr,
                            string filename) const {
    ofstream output_file(filename.c_str());
    for (int i = 0; i < nrhob; i++) {
        for (int j = 0; j < ne; j++) {
            output_file << scientific << setw(18) << setprecision(8)
                        << matrix_ptr[i][j] << "  ";
        }
        output_file << endl;
    }
    output_file.close();
}


double EOS::calculate_velocity_of_sound_sq(double e, double rhob) const {
    double v_min = 0.01;
    double v_max = 1./3;
    double dpde = p_e_func(e, rhob);
    double dpdrho = p_rho_func(e, rhob);
    double pressure = get_pressure(e, rhob);
    double v_sound = dpde + rhob/(e + pressure + 1e-15)*dpdrho;
    v_sound = std::max(v_min, std::min(v_max, v_sound));
    return(v_sound);
}

    
//! This function returns the local pressure in [1/fm^4]
//! the input local energy density [1/fm^4], rhob [1/fm^3]
double EOS::get_pressure(double e, double rhob) const {
    double f = 0.0;
    if (whichEOS == 0) {
        f = ideal_cs2*e;
    } else if (whichEOS == 1) {
        int table_idx = get_table_idx(e);
        f = interpolate2D(e, std::abs(rhob), table_idx, pressure_tb);
        f = f/hbarc;  // 1/fm^4
    } else if (whichEOS >= 2 && whichEOS < 8) {
        int table_idx = get_table_idx(e);
        f = interpolate1D(e, table_idx, pressure_tb);
        f = f/hbarc;  // 1/fm^4
    } else if (whichEOS == 8) {
        f = get_pressure_WB(e);
    } else if (whichEOS >= 10) {
        // EOS is symmetric in rho_b for pressure
        int table_idx = get_table_idx(e);
        f = interpolate2D(e, std::abs(rhob), table_idx, pressure_tb);
        f = f/hbarc;  // 1/fm^4
    }
    return(f);
}


double EOS::p_rho_func(double e, double rhob) const {
    // return dP/drho_b (in 1/fm)
    double f = 0.0;
    if (whichEOS == 1 || whichEOS >= 10) {
        f = get_dpOverdrhob2(e, rhob);
    }
    return(f);
}


double EOS::p_e_func(double e, double rhob) const {
    // return dP/de
    double f;
    if (whichEOS == 8) {
        f = get_dpOverde_WB(e);
    } else {
        f = get_dpOverde3(e, rhob);
    }
    return(f);
}


double EOS::get_dpOverde_WB(double e_local) const {
    double cs2_local;
    double e1 = e_local;
	double e2 = e1*e1;
	double e3 = e2*e1;
	double e4 = e3*e1;
	double e5 = e4*e1;
	double e6 = e5*e1;
	double e7 = e6*e1;
	double e8 = e7*e1;
	double e9 = e8*e1;
	double e10 = e9*e1;
	double e11 = e10*e1;
	double e12 = e11*e1;
	double e13 = e12*e1;
	cs2_local = ((5.191934309650155e-32 + 4.123605749683891e-23*e1
                 + 3.1955868410879504e-16*e2 + 1.4170364808063119e-10*e3
                 + 6.087136671592452e-6*e4 + 0.02969737949090831*e5
                 + 15.382615282179595*e6 + 460.6487249985994*e7
                 + 1612.4245252438795*e8 + 275.0492627924299*e9
                 + 58.60283714484669*e10 + 6.504847576502024*e11
                 + 0.03009027913262399*e12 + 8.189430244031285e-6*e13)
		        /(1.4637868900982493e-30 + 6.716598285341542e-22*e1
                  + 3.5477700458515908e-15*e2 + 1.1225580509306008e-9*e3
                  + 0.00003551782901018317*e4 + 0.13653226327408863*e5
                  + 60.85769171450653*e6 + 1800.5461219450308*e7
                  + 15190.225535036281*e8 + 590.2572000057821*e9
                  + 293.99144775704605*e10 + 21.461303090563028*e11
                  + 0.09301685073435291*e12 + 0.000024810902623582917*e13));
    return(cs2_local);
}


double EOS::get_pressure_WB(double e_local) const {
    double p;
    double e1 = e_local;
    double e2 = e1*e_local;
    double e3 = e2*e_local;
    double e4 = e3*e_local;
    double e5 = e4*e_local;
    double e6 = e5*e_local;
    double e7 = e6*e_local;
    double e8 = e7*e_local;
    double e9 = e8*e_local;
    double e10 = e9*e_local;
    double e11 = e10*e_local;
    double e12 = e11*e_local;
	
	p = ((  1.9531729608963267e-11*e12 + 3.1188455176941583e-7*e11
          + 0.0009417586777847889*e10 + 0.7158279081255019*e9
          + 141.5073484468774*e8 + 6340.448389300905*e7
          + 41913.439282708554*e6 + 334334.4309240126*e5
          + 1.6357487344679043e6*e4 + 3.1729694865420084e6*e3
          + 1.077580993288114e6*e2 + 9737.845799644809*e1
          - 0.25181736420168666)
         /(  3.2581066229887368e-18*e12 + 5.928138360995685e-11*e11
           + 9.601103399348206e-7*e10 + 0.002962497695527404*e9
           + 2.3405487982094204*e8 + 499.04919730607065*e7
           + 26452.34905933697*e6 + 278581.2989342773*e5
           + 1.7851642641834426e6*e4 + 1.3512402226067686e7*e3
           + 2.0931169138134286e7*e2 + 4.0574329080826794e6*e1
           + 45829.44617893836));
    p = std::max(1e-16, p);
    return(p);
}


double EOS::get_temperature_WB(double e_local) const {
    double temperature;
    double e1 = e_local;
	double e2 = e1*e1;
	double e3 = e2*e1;
	double e4 = e3*e1;
	double e5 = e4*e1;
	double e6 = e5*e1;
	double e7 = e6*e1;
	double e8 = e7*e1;
	double e9 = e8*e1;
	double e10 = e9*e1;
	double e11 = e10*e1;
	temperature = ((1.510073201405604e-29 + 8.014062800678687e-18*e1
                    + 2.4954778310451065e-10*e2 + 0.000063810382643387*e3
                    + 0.4873490574161924*e4 + 207.48582344326206*e5
                    + 6686.07424325115*e6 + 14109.766109389702*e7
                    + 1471.6180520527757*e8 + 14.055788949565482*e9
                    + 0.015421252394182246*e10 + 1.5780479034557783e-6*e11)
                   /(7.558667139355393e-28 + 1.3686372302041508e-16*e1
                     + 2.998130743142826e-9*e2 + 0.0005036835870305458*e3
                     + 2.316902328874072*e4 + 578.0778724946719*e5
                     + 11179.193315394154*e6 + 17965.67607192861*e7
                     + 1051.0730543534657*e8 + 5.916312075925817*e9
                     + 0.003778342768228011*e10 + 1.8472801679382593e-7*e11));
    return(temperature);
}


int EOS::get_table_idx(double e) const {
    double local_ed = e*hbarc;  // [GeV/fm^3]
    for (int itable = 1; itable < number_of_tables; itable++) {
        if (local_ed < e_bounds[itable]) {
            return(itable - 1);
        }
    }
    return(number_of_tables - 1);
}

double EOS::interpolate2D(double e, double rhob, int table_idx, double ***table) const {
// This is a generic bilinear interpolation routine for EOS at finite mu_B
// it assumes the class has already read in
//        P(e, rho_b), T(e, rho_b), s(e, rho_b), mu_b(e, rho_b)
// as two-dimensional arrays on an equally spacing lattice grid
// units: e is in 1/fm^4, rhob is in 1/fm^3
    double local_ed = e*hbarc;  // [GeV/fm^3]
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
    double temp1 = std::max(table[table_idx][idx_nb][idx_e], 0.0);
    double temp2 = std::max(table[table_idx][idx_nb][idx_e + 1], 0.0);
    double temp3 = std::max(table[table_idx][idx_nb + 1][idx_e + 1], 0.0);
    double temp4 = std::max(table[table_idx][idx_nb + 1][idx_e], 0.0);
    result = ((temp1*(1. - frac_e) + temp2*frac_e)*(1. - frac_rhob)
              + (temp3*frac_e + temp4*(1. - frac_e))*frac_rhob);
    result = std::max(result, 1e-15);
    return(result);
}


double EOS::interpolate1D(double e, int table_idx, double ***table) const {
// This is a generic linear interpolation routine for EOS at zero mu_B
// it assumes the class has already read in
//        P(e), T(e), s(e)
// as one-dimensional arrays on an equally spacing lattice grid
// units: e is in 1/fm^4
    double local_ed = e*hbarc;  // [GeV/fm^3]

    double e0       = e_bounds[table_idx];
    double delta_e  = e_spacing[table_idx];
    int N_e  = e_length[table_idx];

    // compute the indices
    int idx_e  = static_cast<int>((local_ed - e0)/delta_e);

    // treatment for overflow, use the last two points to do extrapolation
    idx_e  = std::min(N_e - 2, idx_e);

    // check underflow
    idx_e  = std::max(0, idx_e);

    double frac_e    = (local_ed - (idx_e*delta_e + e0))/delta_e;

    double result;
    double temp1 = std::max(table[table_idx][0][idx_e], 0.0);
    double temp2 = std::max(table[table_idx][0][idx_e + 1], 0.0);
    result = temp1*(1. - frac_e) + temp2*frac_e;
    result = std::max(1e-15, result);
    return(result);
}


double EOS::T_from_eps_ideal_gas(double eps) const {
    // Define number of colours and of flavours
    const double Nc = 3;
    const double Nf = 2.5;
    return pow(90.0/M_PI/M_PI*(eps/3.0)/(2*(Nc*Nc-1)+7./2*Nc*Nf), .25);
}


double EOS::s2e_ideal_gas(double s) const {
    // Define number of colours and of flavours
    double Nc = 3;
    double Nf = 2.5;

    //e=T*T*T*T*(M_PI*M_PI*3.0*(2*(Nc*Nc-1)+7./2*Nc*Nf)/90.0);
    //s = 4 e / (3 T)
    //s =4/3 T*T*T*(M_PI*M_PI*3.0*(2*(Nc*Nc-1)+7./2*Nc*Nf)/90.0);
    //T = pow(3. * s / 4. / (M_PI*M_PI*3.0*(2*(Nc*Nc-1)+7./2*Nc*Nf)/90.0), 1./3.);
    return 3./4.*s*pow(3.*s/4./(M_PI*M_PI*3.0*(2*(Nc*Nc-1)+7./2*Nc*Nf)/90.0), 1./3.); //in 1/fm^4
}


//! This function returns entropy density in [1/fm^3]
//! The input local energy density e [1/fm^4], rhob[1/fm^3]
double EOS::get_entropy(double epsilon, double rhob) const {
    double P = get_pressure(epsilon, rhob);
    double T = get_temperature(epsilon, rhob);
    double mu = get_mu(epsilon, rhob);
    double f = (epsilon + P - mu*rhob)/(T + 1e-15);
    return(std::max(1e-16, f));
}/* get_entropy */

//! This function returns the local temperature in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS::get_temperature(double eps, double rhob) const {
    double T = 0.0;
    if (whichEOS == 0) {
        T = T_from_eps_ideal_gas(eps);
    } else if (whichEOS == 1) {
        int table_idx = get_table_idx(eps);
        T = interpolate2D(eps, std::abs(rhob), table_idx,
                          temperature_tb)/hbarc;  // 1/fm
    } else if (whichEOS < 8) {
        int table_idx = get_table_idx(eps);
        T = interpolate1D(eps, table_idx, temperature_tb)/hbarc;  // 1/fm
    } else if (whichEOS == 8) {
        T = get_temperature_WB(eps);
    } else if (whichEOS >= 10) {
        int table_idx = get_table_idx(eps);
        T = interpolate2D(eps, std::abs(rhob), table_idx,
                          temperature_tb)/hbarc;  // 1/fm
    }
    return(std::max(1e-15, T));
}


//! This function returns the local baryon chemical potential  mu_B in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS::get_mu(double eps, double rhob) const {
    double mu = 0.0;
    if (whichEOS >= 10) {
        int table_idx = get_table_idx(eps);
        double sign = rhob/(std::abs(rhob) + 1e-15);
        mu = sign*interpolate2D(eps, std::abs(rhob), table_idx,
                                    mu_B_tb)/hbarc;  // 1/fm
    }
    return(mu);
}


double EOS::get_muS(double eps, double rhob) const {
    // return mu_S in [1/fm]
    double mu = 0.0;
    if (whichEOS == 11) {
        int table_idx = get_table_idx(eps);
        double sign = rhob/(std::abs(rhob) + 1e-15);
        mu = sign*interpolate2D(eps, std::abs(rhob), table_idx,
                                    mu_S_tb)/hbarc;  // 1/fm
    }
    return(mu);
}


double EOS::get_s2e(double s, double rhob) const {
    // s - entropy density in 1/fm^3
    double e;  // epsilon - energy density
    if (whichEOS == 0) {
        e = s2e_ideal_gas(s);
    } else {
        e = get_s2e_finite_rhob(s, rhob);
    }
    return(e);  // in 1/fm^4
}


//! This function returns local energy density [1/fm^4] from
//! a given entropy density [1/fm^3] and rhob [1/fm^3]
//! using binary search
double EOS::get_s2e_finite_rhob(double s, double rhob) const {
    double eps_lower = 1e-15;
    double eps_upper = eps_max;
    double eps_mid   = (eps_upper + eps_lower)/2.;
    double s_lower   = get_entropy(eps_lower, rhob);
    double s_upper   = get_entropy(eps_upper, rhob);
    int ntol         = 1000;
    if (s < 0.0 || s > s_upper) {
        fprintf(stderr, "get_s2e_finite_rhob:: s is out of bound, "
                        "s = %.5e, s_upper = %.5e, s_lower = %.5e \n",
                        s, s_upper, s_lower);
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
        fprintf(stderr, "get_s2e_finite_rhob:: max iteration reached, "
                        "s = %.5e, rhob = %.5e \n", s, rhob);
        fprintf(stderr, "s_upper = %.5e, s_lower = %.5e \n",
                get_entropy(eps_upper, rhob), get_entropy(eps_lower, rhob));
        fprintf(stderr, "eps_upper = %.5e, eps_lower = %.5e, diff = %.10e \n",
                eps_upper, eps_lower, (eps_upper - eps_lower));
        exit(1);
    }
    return (eps_mid);
}

//! This function returns local net baryon density rhob [1/fm^3]
//! from given local energy density e [1/fm^4] and mu_b [1/fm]
double EOS::get_rhob_from_mub(double e, double mub) const {
    double local_ed = e*hbarc;      // GeV/fm^3
    double local_mub = mub*hbarc;   // GeV

    int table_idx    = get_table_idx(e);
    int NEps         = e_length[table_idx];
    int Nrhob        = nb_length[table_idx];
    double eps0      = e_bounds[table_idx];
    double deltaEps  = e_spacing[table_idx];
    double rhob0     = nb_bounds[table_idx];
    double deltaRhob = nb_spacing[table_idx];

    // compute the indices
    int idx_e = static_cast<int>((local_ed - eps0)/deltaEps);
    double frac_e = (local_ed - (idx_e*deltaEps + eps0))/deltaEps; 
    
    // check overflow and underflow
    idx_e = std::max(0, std::min(NEps - 2, idx_e));

    double *array_left  = new double [Nrhob];
    double *array_right = new double [Nrhob];

    for (int i = 0; i < Nrhob; i++) {
       array_left[i]  = mu_B_tb[table_idx][i][idx_e];
       array_right[i] = mu_B_tb[table_idx][i][idx_e+1];
    }

    int idx_rhob_left    = Util::binary_search(array_left, Nrhob, local_mub);
    int idx_rhob_right   = Util::binary_search(array_right, Nrhob, local_mub);
    double rhob_left_1   = rhob0 + idx_rhob_left*deltaRhob;
    double rhob_left_2   = rhob0 + (idx_rhob_left+1)*deltaRhob;
    double mub_left_1    = array_left[idx_rhob_left];
    double mub_left_2    = array_left[idx_rhob_left+1];
    double frac_mub_left = (local_mub - mub_left_1)/(mub_left_2 - mub_left_1);
    double rhob_left     = (rhob_left_1*(1. - frac_mub_left)
                            + rhob_left_2*frac_mub_left);

    double rhob_right_1   = rhob0 + idx_rhob_right*deltaRhob;
    double rhob_right_2   = rhob0 + (idx_rhob_right+1)*deltaRhob;
    double mub_right_1    = array_right[idx_rhob_right];
    double mub_right_2    = array_right[idx_rhob_right+1];
    double frac_mub_right = ((local_mub - mub_right_1)
                             /(mub_right_2 - mub_right_1));
    double rhob_right = (rhob_right_1*(1. - frac_mub_right)
                         + rhob_right_2*frac_mub_right);

    double rhob = rhob_left*(1. - frac_e) + rhob_right*frac_e;   // 1/fm^3
    return(rhob);
}


//! This is a shell function to check EoS
void EOS::check_eos() const {
    if (whichEOS >= 2 && whichEOS < 10) {
        check_eos_no_muB();
    } else if (whichEOS == 10) {
        check_eos_with_finite_muB();
    } else if (whichEOS == 12) {
        check_eos_with_finite_muB();
    }
}


void EOS::check_eos_no_muB() const {
    // output EoS as function of e
    ostringstream file_name;
    file_name << "check_EoS_PST.dat";
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

void EOS::check_eos_with_finite_muB() const {
    // output EoS as function of e for several rhob
    double rhob_pick[6] = {0.0, 0.02, 0.05, 0.1, 0.2, 0.5};
    for (int i = 0; i < 6; i++) {
        double rhob_local = rhob_pick[i];
        ostringstream file_name;
        file_name << "check_EoS_PST_rhob_" << rhob_pick[i] << ".dat";
        ofstream check_file(file_name.str().c_str());
        check_file << "#e(GeV/fm^3)  P(GeV/fm^3)  s(1/fm^3)  T(GeV)  cs^2  "
                   << "mu_B(GeV)  mu_S(GeV)" << endl;
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
            double mu_b_local = get_mu(e_local, rhob_local);
            double mu_s_local = get_muS(e_local, rhob_local);
            check_file << scientific << setw(18) << setprecision(8)
                       << e_local*hbarc << "   " << p_local*hbarc << "   " 
                       << s_local << "   " << T_local*hbarc << "   "
                       << cs2_local << "   " << mu_b_local*hbarc << "   "
                       << mu_s_local*hbarc << endl;
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
                   << "mu_B(GeV)  mu_S(GeV)" << endl;
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
            double mu_b_local = get_mu(e_local, rhob_local);
            double mu_s_local = get_muS(e_local, rhob_local);
            check_file << scientific << setw(18) << setprecision(8)
                       << rhob_local << "   " << p_local*hbarc << "   " 
                       << s_local << "   " << T_local*hbarc << "   "
                       << cs2_local << "   " << mu_b_local*hbarc << "   "
                       << mu_s_local*hbarc << endl;
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

    output_eos_matrix(e_length[0], nb_length[0], temperature_tb[0],
                      "check_EoS_T_table1.dat");
    output_eos_matrix(e_length[0], nb_length[0], mu_B_tb[0],
                      "check_EoS_muB_table1.dat");
    
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
                    << "s(1/fm^3)  rho_B(1/fm^3)  dP/de  dP/drho" << endl;
        for (int j = 0; j < ns; j++) {
            double s_local     = s_0 + j*ds;
            double nB_local    = s_local/sovernB[i];
            double e_local     = get_s2e_finite_rhob(s_local, nB_local);
            double s_check     = get_entropy(e_local, nB_local);
            double cs2_local   = get_cs2(e_local, nB_local);
            double dpde        = p_e_func(e_local, nB_local);
            double dpdrho      = p_rho_func(e_local, nB_local);
            double temperature = get_temperature(e_local, nB_local)*hbarc;
            double mu_B        = get_mu(e_local, nB_local)*hbarc;
            check_file9 << scientific << setw(18) << setprecision(8)
                        << e_local*hbarc << "  " << temperature << "  "
                        << cs2_local << "  " << mu_B << "  " 
                        << s_check << "  " << nB_local << "  "
                        << dpde << "  " << dpdrho << endl;
        }
        check_file9.close();
    }

    return;
}

std::string EOS::get_hydro_env_path() const {
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

void EOS::resize_table_info_arrays() {
    nb_bounds.resize(number_of_tables, 0.0);
    nb_spacing.resize(number_of_tables, 0.0);
    nb_length.resize(number_of_tables, 0);
    e_bounds.resize(number_of_tables, 0.0);
    e_spacing.resize(number_of_tables, 0.0);
    e_length.resize(number_of_tables, 0);
}
