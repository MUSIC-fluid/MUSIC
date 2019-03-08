// Copyright @ Chun Shen

#include "doctest.h"
#include "read_in_parameters.h"
#include "critical_modes.h"
#include "grid.h"
#include <iostream>
#include <iomanip>


TEST_CASE("Check CriticalSlowModes initialization") {
    EOS eos_ideal(0);
    InitData DATA(ReadInParameters::read_in_parameters(
                            "tests/unittest_files/music_input_criticalmodes"));
    CriticalSlowModes test(eos_ideal, DATA);
    SCGrid arena_current(10, 10, 5);

    test.InitializeFields(10, arena_current);
    CHECK(test.get_Qvec_size() == 10);
    test.InitializeFields(20, arena_current);
    CHECK(test.get_Qvec_size() == 20);
}

TEST_CASE("TEST phiQ evolution static medium") {
    EOS eos_ideal(0);
    InitData DATA(ReadInParameters::read_in_parameters(
                            "tests/unittest_files/music_input_criticalmodes"));
    CriticalSlowModes test(eos_ideal, DATA);
    const int grid_neta = 3;
    const int grid_nx   = 3;
    const int grid_ny   = 3;

    SCGrid arena_current(grid_nx, grid_ny, grid_neta);
    SCGrid arena_prev   (grid_nx, grid_ny, grid_neta);
    SCGrid arena_future (grid_nx, grid_ny, grid_neta);
    SCGrid *ap_current = &arena_current;
    SCGrid *ap_prev    = &arena_prev;
    SCGrid *ap_future  = &arena_future;
  

    // initialize the background
    for (int ieta = 0; ieta < grid_neta; ieta++) {
        for (int ix = 0; ix < grid_nx; ix++) {
            for (int iy = 0; iy < grid_ny; iy++) {
                arena_current(ix, iy, ieta).epsilon = 1.0;
                arena_current(ix, iy, ieta).rhob    = 0.0;
                arena_current(ix, iy, ieta).u[0] = 1.0;
                arena_current(ix, iy, ieta).u[1] = 0.0;
                arena_current(ix, iy, ieta).u[2] = 0.0;
                arena_current(ix, iy, ieta).u[3] = 0.0;
            }
        }
    }
    arena_prev   = arena_current;
    arena_future = arena_current;

    const int nQ = 50;
    test.InitializeFields(nQ, arena_current);
    test.InitializeFields(nQ, arena_prev);
    test.InitializeFields(nQ, arena_future);
    
    const int rk_order = 2;
    for (int it = 0; it <= 1000; it++) {
        double tau_local = DATA.tau0 + it*DATA.delta_tau;
        // rk evolution
        for (int rk_flag = 0; rk_flag < rk_order; rk_flag++) {
            for (int ieta = 0; ieta < grid_neta; ieta++) {
                for (int ix = 0; ix < grid_nx; ix++) {
                    for (int iy = 0; iy < grid_ny; iy++) {
                        double theta_local = 1./(tau_local 
                                                 + rk_flag*DATA.delta_tau);
                        test.evolve_phiQfields(tau_local,
                                               *ap_prev, *ap_current,
                                               *ap_future, theta_local,
                                               ix, iy, ieta, rk_flag);
                    }
                }
            }
            if (rk_flag == 0) {
                SCGrid *temp = ap_prev;
                ap_prev      = ap_current;
                ap_current   = ap_future;
                ap_future    = temp;
            } else {
                SCGrid *temp = ap_future;
                ap_future    = ap_current;
                ap_current   = temp;
            }
        }
        if (it % 100 == 0) {
            for (int iQ = 0; iQ < test.get_Qvec_size(); iQ++) {
                const double Q_local = test.get_Qi(iQ);
                const double phiQ_eq = test.compute_phiQ_equilibrium(
                                Q_local, (*ap_current)(0, 0, 0).epsilon,
                                (*ap_current)(0, 0, 0).rhob);
                std::cout << std::scientific
                          << it << "  " << iQ << "  " << phiQ_eq << "  "
                          << (*ap_current)(0, 0, 0).phi_Q[iQ] << "  "
                          << (*ap_current)(0, 0, 0).phi_Q[iQ]
                          << std::endl;
            }
        }
    }
}

