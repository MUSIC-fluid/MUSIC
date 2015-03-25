#ifndef GRID_INFO_H
#define GRID_INFO_H
#include "data.h"
#include "eos.h"
#include "Grid.h"
#include <mpi.h>
#include <iostream>
#include <iomanip>

class Grid_info
{
    private:

    public:
        Grid_info();
        ~Grid_info();

        void ComputeV2(InitData *DATA, Grid ***arena, double tau); //added
        void ComputeEccentricity(InitData *DATA, Grid ***arena, double tau); //added
        void print_rhob_evolution(InitData *DATA, Grid ***arena, double tau, EOS* eos, int rank);
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
};
#endif
