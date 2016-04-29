#ifndef GRID_INFO_H
#define GRID_INFO_H
#include "data.h"
#include "eos.h"
#include "grid.h"
#include <mpi.h>
#include <iostream>
#include <iomanip>

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

        void ComputeV2(InitData *DATA, Grid ***arena, double tau); //added
        void ComputeEccentricity(InitData *DATA, Grid ***arena, double tau); //added
        void print_rhob_evolution(InitData *DATA, Grid ***arena, double tau, EOS* eos, int rank);
        void print_rhob_evolution_3d(InitData *DATA, Grid ***arena, double tau, EOS* eos, int rank);
        void print_qmu_evolution(InitData *DATA, Grid ***arena, double tau, EOS* eos, int rank);
        void print_fireball_evolution_on_phasediagram(InitData *DATA, Grid ***arena, double tau, EOS* eos, int rank);
        void ComputeAnisotropy(InitData *DATA, Grid ***arena, double tau); //added
        void PrintGrid(Grid *grid_p, int rk_order);
        void LinkNeighbors(InitData *DATA, Grid ****arena);
        void InitTJb(InitData *DATA, Grid ****arena);
        void PrintAxy(InitData *DATA, Grid ***arena, double tau);
        void PrintAxy2(InitData *DATA, Grid ***arena, double tau);
        void PrintdEdEta(InitData *DATA, Grid ***arena);
        void OutputEvolutionDataXYZ(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank);
        void OutputEvolutionDataXYEta(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank);
        void OutputEvolutionDataXYEta_finite_muB(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank);
        void OutputPlotDataXYZ(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank);
        void OutputEvolutionOSCAR(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank);
        void OutputXY(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank);
        void PrintArena(Grid ***arena, InitData *DATA, double tau);
        void PrintEtaEpsilon(Grid ***arena, InitData *DATA, double tau, int size, int rank);
        void PrintxEpsilon(Grid ***arena, InitData *DATA, double tau, int size, int rank);
        void ComputeEnergyConservation(InitData *DATA, Grid ***arena, double tau);
        void getAverageTandPlasmaEvolution(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank);
        void Output_hydro_information_header(InitData *DATA, EOS *eos);
        void Tmax_profile(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank);

        void load_deltaf_qmu_coeff_table(string filename);
        void load_deltaf_qmu_coeff_table_14mom(string filename);
        double get_deltaf_qmu_coeff(double T, double muB);
        double get_deltaf_coeff_14moments(double T, double muB, double type);
        void check_conservation_law(Grid ***arena, InitData *DATA,
                                    double tau, int size, int rank);
};
#endif
