// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
// Copyright 2014-2016 Chun Shen
#ifndef SRC_GRID_INFO_H_
#define SRC_GRID_INFO_H_

#include <iostream>
#include <iomanip>
#include "./data.h"
#include "./eos.h"
#include "./grid.h"
#include "./pretty_ostream.h"

class Cell_info {
 private:
    InitData* DATA_ptr;
    EOS* eos_ptr;
    pretty_ostream music_message;
    
    int deltaf_qmu_coeff_table_length_T;
    int deltaf_qmu_coeff_table_length_mu;
    double delta_qmu_coeff_table_T0, delta_qmu_coeff_table_mu0;
    double delta_qmu_coeff_table_dT, delta_qmu_coeff_table_dmu;
    double **deltaf_qmu_coeff_tb;
    int deltaf_coeff_table_14mom_length_T;
    int deltaf_coeff_table_14mom_length_mu;
    double delta_coeff_table_14mom_T0, delta_coeff_table_14mom_mu0;
    double delta_coeff_table_14mom_dT, delta_coeff_table_14mom_dmu;
    double **deltaf_coeff_tb_14mom_DPi, **deltaf_coeff_tb_14mom_BPi;
    double **deltaf_coeff_tb_14mom_BPitilde;
    double **deltaf_coeff_tb_14mom_BV, **deltaf_coeff_tb_14mom_DV;
    double **deltaf_coeff_tb_14mom_Bpi_shear;

 public:
    Cell_info(InitData* DATA_in, EOS *eos_ptr_in);
    ~Cell_info();

    //! This function outputs a header files for JF and Gojko's EM programs
    void Output_hydro_information_header(InitData *DATA);

    //! This function outputs hydro evolution file in binary format
    void OutputEvolutionDataXYEta(Cell ***arena, InitData *DATA,
                                  double tau);

    //! This function outputs hydro evolution file in binary format
    void OutputEvolutionDataXYEta_chun(Cell ***arena, InitData *DATA,
                                       double tau);

    void load_deltaf_qmu_coeff_table(string filename);
    void load_deltaf_qmu_coeff_table_14mom(string filename);
    double get_deltaf_qmu_coeff(double T, double muB);
    double get_deltaf_coeff_14moments(double T, double muB, double type);

    //! This function putputs files to check with Gubser flow solution
    void Gubser_flow_check_file(Cell ***arena, double tau);

    //! This function outputs files to cross check with 1+1D simulation
    void output_1p1D_check_file(Cell ***arena, double tau);

    //! This function prints to the screen the maximum local energy density,
    //! the maximum temperature in the current grid
    void get_maximum_energy_density(Cell ***arena);

    //! This function outputs energy density and n_b for making movies
    void output_evolution_for_movie(Cell ***arena, double tau);

    //! This function outputs average T and mu_B as a function of proper tau
    //! within a given space-time rapidity range
    void output_average_phase_diagram_trajectory(
                double tau, double eta_min, double eta_max, Cell ***arena);

    //! This function dumps the energy density and net baryon density
    void output_energy_density_and_rhob_disitrubtion(Cell ***arena,
                                                     string filename);
    
    //! This function checks the total energy and total net baryon number
    //! at a give proper time
    void check_conservation_law(Cell ***arena, InitData *DATA, double tau);

    //! This function outputs the evolution of hydrodynamic variables at a
    //! give fluid cell
    void monitor_fluid_cell(Cell ***arena, int ix, int iy, int ieta,
                            double tau);
};

#endif  // SRC_GRID_INFO_H_
