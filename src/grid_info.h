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

class Grid_info {
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
    Grid_info(InitData* DATA_in, EOS *eos_ptr_in);
    ~Grid_info();

    void Output_hydro_information_header(InitData *DATA);
    void OutputEvolutionDataXYEta(Grid ***arena, InitData *DATA,
                                  double tau);
    void OutputEvolutionDataXYEta_chun(Grid ***arena, InitData *DATA,
                                       double tau);
    void Gubser_flow_check_file(Grid ***arena, double tau);
    void output_1p1D_check_file(Grid ***arena, double tau);
    void load_deltaf_qmu_coeff_table(string filename);
    void load_deltaf_qmu_coeff_table_14mom(string filename);
    double get_deltaf_qmu_coeff(double T, double muB);
    double get_deltaf_coeff_14moments(double T, double muB, double type);

    //! This function prints to the screen the maximum local energy density,
    //! the maximum temperature in the current grid
    void get_maximum_energy_density(Grid ***arena);

    void check_conservation_law(Grid ***arena, InitData *DATA, double tau);
    void check_velocity_shear_tensor(Grid ***arena, double tau);
    void output_evolution_for_movie(Grid ***arena, double tau);
    void output_average_phase_diagram_trajectory(
                double tau, double eta_min, double eta_max, Grid ***arena);
};

#endif  // SRC_GRID_INFO_H_
