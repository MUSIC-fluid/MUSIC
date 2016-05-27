// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_GRID_INFO_H_
#define SRC_GRID_INFO_H_
#include <iostream>
#include <iomanip>
#include "./data.h"
#include "./eos.h"
#include "./grid.h"

class Grid_info
{
    private:
        InitData* DATA_ptr;
        
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
        Grid_info(InitData* DATA_in);
        ~Grid_info();

        void OutputEvolutionDataXYEta(Grid ***arena, InitData *DATA,
                                      EOS *eos, double tau);
        void Gubser_flow_check_file(Grid ***arena, EOS *eos, double tau);
        void load_deltaf_qmu_coeff_table(string filename);
        void load_deltaf_qmu_coeff_table_14mom(string filename);
        double get_deltaf_qmu_coeff(double T, double muB);
        double get_deltaf_coeff_14moments(double T, double muB, double type);
        void check_conservation_law(Grid ***arena, InitData *DATA, double tau);
        void check_velocity_shear_tensor(Grid ***arena, double tau);
};

#endif  // SRC_GRID_INFO_H_
