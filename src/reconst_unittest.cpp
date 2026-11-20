
#include "reconst.h"

#include <iostream>

#include "cell.h"
#include "data_struct.h"
#include "doctest.h"
#include "eos.h"

TEST_CASE("Check Reconst constructor") {
    EOS eos_ideal(0);
    Reconst reconst_test(eos_ideal, 9, 0);  // beastMode 0: tight tolerances
    CHECK(reconst_test.get_max_iter() == 100);
    CHECK(reconst_test.get_echo_level() == 9);
    CHECK(reconst_test.get_abs_err() == doctest::Approx(1e-16).epsilon(1e-16));
    CHECK(
        reconst_test.get_v_critical()
        == doctest::Approx(0.563624).epsilon(0.00001));
}

TEST_CASE("Test Newton solver v") {
    EOS eos_ideal(0);
    Reconst reconst_test(eos_ideal, 9, 0);  // beastMode 0: tight tolerances
    // the solvers stop at rel_err = 1e-14, so compare to that
    double v_sol = 10.0;
    // v = 0 is the root for M = 0 (the finite-difference derivative used
    // to be 0/0 there)
    CHECK(reconst_test.solve_velocity_Newton(0.0, 1.0, 0.0, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(0.0).epsilon(1e-16));

    double M0 = 1.0;
    double M = 0.5;
    double v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));

    M0 = 1.0;
    M = 0.99;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));

    M0 = 1.0;
    M = 1.0 - 1e-6;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));

    M0 = 1.0;
    M = 1e-6;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));

    M0 = 1.0;
    M = 1.0 - 1e-10;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));

    M0 = 1.0;
    M = 1e-10;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));

    M0 = 1e-8;
    M = 1e-12;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));

    M0 = 1e-12;
    M = 3.23e-15;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));

    M0 = 1e-14;
    M = M0 - 2.43e-18;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));
}

TEST_CASE("Test hybrid solver v") {
    EOS eos_ideal(0);
    Reconst reconst_test(eos_ideal, 9, 0);  // beastMode 0: tight tolerances
    // the solvers stop at rel_err = 1e-14, so compare to that
    double v_sol = 10.0;
    CHECK(reconst_test.solve_v_Hybrid(0.0, 1.0, 0.0, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(0.0).epsilon(1e-16));

    double M0 = 1.0;
    double M = 0.5;
    double v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_v_Hybrid(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));

    M0 = 1.0;
    M = 0.99;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_v_Hybrid(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));

    M0 = 1.0;
    M = 1.0 - 1e-6;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_v_Hybrid(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));

    M0 = 1.0;
    M = 1e-6;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_v_Hybrid(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));
    M0 = 1.0;
    M = 1.0 - 1e-10;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_v_Hybrid(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));

    M0 = 1.0;
    M = 1e-10;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_v_Hybrid(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));

    M0 = 1e-8;
    M = 1e-12;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_v_Hybrid(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));

    M0 = 1e-12;
    M = 3.23e-15;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_v_Hybrid(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));

    M0 = 1e-14;
    M = M0 - 2.43e-18;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    std::cout << "check v = " << v_correct << std::endl;
    CHECK(reconst_test.solve_v_Hybrid(0.0, M0, M, 0.0, v_sol) == 1);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-14));
}

TEST_CASE("Test Newton solver u0") {
    EOS eos_ideal(0);
    Reconst reconst_test(eos_ideal, 9, 0);  // beastMode 0: tight tolerances

    double u0_sol = 0.0;
    double M0 = 1.0;
    double M = 0.0;
    double v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    double u0_correct = 1. / sqrt(1. - v_correct * v_correct);
    double utol = u0_correct * 1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    reconst_test.solve_u0_Newton(100.0, M0, M * M, M, 0.0, u0_sol);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));

    M0 = 1.0;
    M = 0.5;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    u0_correct = 1. / sqrt(1. - v_correct * v_correct);
    utol = u0_correct * 1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    reconst_test.solve_u0_Newton(100.0, M0, M * M, M, 0.0, u0_sol);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));

    M0 = 1.0;
    M = 1e-4;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    u0_correct = 1. / sqrt(1. - v_correct * v_correct);
    utol = u0_correct * 1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    reconst_test.solve_u0_Newton(100.0, M0, M * M, M, 0.0, u0_sol);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));

    M0 = 1.0;
    M = 0.999;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    u0_correct = 1. / sqrt(1. - v_correct * v_correct);
    utol = u0_correct * 1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    reconst_test.solve_u0_Newton(100.0, M0, M * M, M, 0.0, u0_sol);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));

    M0 = 1.0;
    M = 1.0 - 1e-4;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    u0_correct = 1. / sqrt(1. - v_correct * v_correct);
    utol = u0_correct * 1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    reconst_test.solve_u0_Newton(100.0, M0, M * M, M, 0.0, u0_sol);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));

    M0 = 1e-15;
    M = M0 - 3.463e-18;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    u0_correct = 1. / sqrt(1. - v_correct * v_correct);
    utol = u0_correct * 1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    reconst_test.solve_u0_Newton(100.0, M0, M * M, M, 0.0, u0_sol);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));

    M0 = 1e-15;
    M = 3.463e-18;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    u0_correct = 1. / sqrt(1. - v_correct * v_correct);
    utol = u0_correct * 1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    reconst_test.solve_u0_Newton(100.0, M0, M * M, M, 0.0, u0_sol);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));
}

TEST_CASE("Test Newton solver u0 hybrid") {
    EOS eos_ideal(0);
    Reconst reconst_test(eos_ideal, 9, 0);  // beastMode 0: tight tolerances

    double u0_sol = 0.0;
    double M0 = 1.0;
    double M = 0.0;
    double v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    double u0_correct = 1. / sqrt(1. - v_correct * v_correct);
    double utol = u0_correct * 1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    CHECK(reconst_test.solve_u0_Hybrid(0.0, M0, M * M, M, 0.0, u0_sol) == 1);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));

    M0 = 1.0;
    M = 0.5;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    u0_correct = 1. / sqrt(1. - v_correct * v_correct);
    utol = u0_correct * 1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    // a guess >= 1 only brackets [0.99, 1.01] x guess (in MUSIC the guess
    // comes from the velocity solution); a far-off guess must report failure
    CHECK(reconst_test.solve_u0_Hybrid(2.0, M0, M * M, M, 0.0, u0_sol) == 0);
    CHECK(
        reconst_test.solve_u0_Hybrid(
            1.01 * u0_correct, M0, M * M, M, 0.0, u0_sol)
        == 1);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));

    M0 = 1.0;
    M = 1e-4;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    u0_correct = 1. / sqrt(1. - v_correct * v_correct);
    utol = u0_correct * 1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    CHECK(reconst_test.solve_u0_Hybrid(1.0, M0, M * M, M, 0.0, u0_sol) == 1);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));

    M0 = 1.0;
    M = 0.999;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    u0_correct = 1. / sqrt(1. - v_correct * v_correct);
    utol = u0_correct * 1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    CHECK(reconst_test.solve_u0_Hybrid(0.0, M0, M * M, M, 0.0, u0_sol) == 1);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));

    M0 = 1.0;
    M = 1.0 - 1e-4;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    u0_correct = 1. / sqrt(1. - v_correct * v_correct);
    utol = u0_correct * 1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    CHECK(reconst_test.solve_u0_Hybrid(0.0, M0, M * M, M, 0.0, u0_sol) == 1);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));

    M0 = 1e-15;
    M = M0 - 3.463e-18;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    u0_correct = 1. / sqrt(1. - v_correct * v_correct);
    utol = u0_correct * 1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    CHECK(reconst_test.solve_u0_Hybrid(0.0, M0, M * M, M, 0.0, u0_sol) == 1);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));

    M0 = 1e-15;
    M = 3.463e-18;
    v_correct = 3. * M / (2. * M0 + sqrt(4. * M0 * M0 - 3. * M * M));
    u0_correct = 1. / sqrt(1. - v_correct * v_correct);
    utol = u0_correct * 1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    CHECK(reconst_test.solve_u0_Hybrid(0.0, M0, M * M, M, 0.0, u0_sol) == 1);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));
}

TEST_CASE("Test reconst main function") {
    EOS eos_ideal(0);
    Reconst reconst_test(eos_ideal, 9, 0);  // beastMode 0: tight tolerances

    const double rel_tol = 1e-12;
    const double abs_tol = 1e-16;

    double tau = 2.0;
    double e_local = 1.0;
    double rhob_local = 0.0;
    double p_local = eos_ideal.get_pressure(e_local, rhob_local);
    double ux = 0.0;
    double uy = 0.0;
    double ueta = 0.0;
    double utau = sqrt(1. + ux * ux + uy * uy + ueta * ueta);
    TJbVec tauq_vec = {
        tau * ((e_local + p_local) * utau * utau - p_local),
        tau * (e_local + p_local) * utau * ux,
        tau * (e_local + p_local) * utau * uy,
        tau * (e_local + p_local) * utau * ueta, tau * rhob_local * utau};
    ReconstCell temp_grid;
    temp_grid.e = 1e-5;
    temp_grid.rhob = 0.0;
    temp_grid.u = {1.0, 0.0, 0.0, 0.0};

    auto cell_sol = reconst_test.ReconstIt_shell(tau, tauq_vec, temp_grid);
    CHECK(
        cell_sol.e
        == doctest::Approx(e_local).epsilon(
            std::max(abs_tol, e_local * rel_tol)));
    CHECK(
        cell_sol.rhob
        == doctest::Approx(rhob_local)
               .epsilon(std::max(abs_tol, rhob_local * rel_tol)));
    CHECK(
        cell_sol.u[0]
        == doctest::Approx(utau).epsilon(std::max(abs_tol, utau * rel_tol)));
    CHECK(
        cell_sol.u[1]
        == doctest::Approx(ux).epsilon(std::max(abs_tol, ux * rel_tol)));
    CHECK(
        cell_sol.u[2]
        == doctest::Approx(uy).epsilon(std::max(abs_tol, uy * rel_tol)));
    CHECK(
        cell_sol.u[3]
        == doctest::Approx(ueta).epsilon(std::max(abs_tol, ueta * rel_tol)));

    tau = 2.4;
    e_local = 12.3;
    rhob_local = 0.2;
    p_local = eos_ideal.get_pressure(e_local, rhob_local);
    ux = 0.1;
    uy = 0.4;
    ueta = 0.2;
    utau = sqrt(1. + ux * ux + uy * uy + ueta * ueta);
    tauq_vec = {
        tau * ((e_local + p_local) * utau * utau - p_local),
        tau * (e_local + p_local) * utau * ux,
        tau * (e_local + p_local) * utau * uy,
        tau * (e_local + p_local) * utau * ueta, tau * rhob_local * utau};
    cell_sol = reconst_test.ReconstIt_shell(tau, tauq_vec, temp_grid);
    CHECK(
        cell_sol.e
        == doctest::Approx(e_local).epsilon(
            std::max(abs_tol, e_local * rel_tol)));
    CHECK(
        cell_sol.rhob
        == doctest::Approx(rhob_local)
               .epsilon(std::max(abs_tol, rhob_local * rel_tol)));
    CHECK(
        cell_sol.u[0]
        == doctest::Approx(utau).epsilon(std::max(abs_tol, utau * rel_tol)));
    CHECK(
        cell_sol.u[1]
        == doctest::Approx(ux).epsilon(std::max(abs_tol, ux * rel_tol)));
    CHECK(
        cell_sol.u[2]
        == doctest::Approx(uy).epsilon(std::max(abs_tol, uy * rel_tol)));
    CHECK(
        cell_sol.u[3]
        == doctest::Approx(ueta).epsilon(std::max(abs_tol, ueta * rel_tol)));

    tau = 5.2;
    e_local = 2.6;
    rhob_local = 23.;
    p_local = eos_ideal.get_pressure(e_local, rhob_local);
    ux = 1.3;
    uy = 2.4;
    ueta = 2.2;
    utau = sqrt(1. + ux * ux + uy * uy + ueta * ueta);
    tauq_vec = {
        tau * ((e_local + p_local) * utau * utau - p_local),
        tau * (e_local + p_local) * utau * ux,
        tau * (e_local + p_local) * utau * uy,
        tau * (e_local + p_local) * utau * ueta, tau * rhob_local * utau};
    cell_sol = reconst_test.ReconstIt_shell(tau, tauq_vec, temp_grid);
    CHECK(
        cell_sol.e
        == doctest::Approx(e_local).epsilon(
            std::max(abs_tol, e_local * rel_tol)));
    CHECK(
        cell_sol.rhob
        == doctest::Approx(rhob_local)
               .epsilon(std::max(abs_tol, rhob_local * rel_tol)));
    CHECK(
        cell_sol.u[0]
        == doctest::Approx(utau).epsilon(std::max(abs_tol, utau * rel_tol)));
    CHECK(
        cell_sol.u[1]
        == doctest::Approx(ux).epsilon(std::max(abs_tol, ux * rel_tol)));
    CHECK(
        cell_sol.u[2]
        == doctest::Approx(uy).epsilon(std::max(abs_tol, uy * rel_tol)));
    CHECK(
        cell_sol.u[3]
        == doctest::Approx(ueta).epsilon(std::max(abs_tol, ueta * rel_tol)));

    tau = 5.2;
    e_local = 1e-4;
    rhob_local = 0.0;
    p_local = eos_ideal.get_pressure(e_local, rhob_local);
    ux = 10.3;
    uy = 2.4;
    ueta = 2.2;
    utau = sqrt(1. + ux * ux + uy * uy + ueta * ueta);
    tauq_vec = {
        tau * ((e_local + p_local) * utau * utau - p_local),
        tau * (e_local + p_local) * utau * ux,
        tau * (e_local + p_local) * utau * uy,
        tau * (e_local + p_local) * utau * ueta, tau * rhob_local * utau};
    cell_sol = reconst_test.ReconstIt_shell(tau, tauq_vec, temp_grid);
    CHECK(
        cell_sol.e
        == doctest::Approx(e_local).epsilon(
            std::max(abs_tol, e_local * rel_tol)));
    CHECK(
        cell_sol.rhob
        == doctest::Approx(rhob_local)
               .epsilon(std::max(abs_tol, rhob_local * rel_tol)));
    CHECK(
        cell_sol.u[0]
        == doctest::Approx(utau).epsilon(std::max(abs_tol, utau * rel_tol)));
    CHECK(
        cell_sol.u[1]
        == doctest::Approx(ux).epsilon(std::max(abs_tol, ux * rel_tol)));
    CHECK(
        cell_sol.u[2]
        == doctest::Approx(uy).epsilon(std::max(abs_tol, uy * rel_tol)));
    CHECK(
        cell_sol.u[3]
        == doctest::Approx(ueta).epsilon(std::max(abs_tol, ueta * rel_tol)));

    tau = 15.2;
    e_local = 2.45e-8;
    rhob_local = 0.0;
    p_local = eos_ideal.get_pressure(e_local, rhob_local);
    ux = 1.3;
    uy = 2.4;
    ueta = 20.2;
    utau = sqrt(1. + ux * ux + uy * uy + ueta * ueta);
    tauq_vec = {
        tau * ((e_local + p_local) * utau * utau - p_local),
        tau * (e_local + p_local) * utau * ux,
        tau * (e_local + p_local) * utau * uy,
        tau * (e_local + p_local) * utau * ueta, tau * rhob_local * utau};
    cell_sol = reconst_test.ReconstIt_shell(tau, tauq_vec, temp_grid);
    CHECK(
        cell_sol.e
        == doctest::Approx(e_local).epsilon(
            std::max(abs_tol, e_local * rel_tol)));
    CHECK(
        cell_sol.rhob
        == doctest::Approx(rhob_local)
               .epsilon(std::max(abs_tol, rhob_local * rel_tol)));
    CHECK(
        cell_sol.u[0]
        == doctest::Approx(utau).epsilon(std::max(abs_tol, utau * rel_tol)));
    CHECK(
        cell_sol.u[1]
        == doctest::Approx(ux).epsilon(std::max(abs_tol, ux * rel_tol)));
    CHECK(
        cell_sol.u[2]
        == doctest::Approx(uy).epsilon(std::max(abs_tol, uy * rel_tol)));
    CHECK(
        cell_sol.u[3]
        == doctest::Approx(ueta).epsilon(std::max(abs_tol, ueta * rel_tol)));

    tau = 7.43;
    e_local = 3.21e-13;
    rhob_local = 0.0;
    p_local = eos_ideal.get_pressure(e_local, rhob_local);
    ux = 1.3;
    uy = 65.4;
    ueta = 0.2;
    utau = sqrt(1. + ux * ux + uy * uy + ueta * ueta);
    tauq_vec = {
        tau * ((e_local + p_local) * utau * utau - p_local),
        tau * (e_local + p_local) * utau * ux,
        tau * (e_local + p_local) * utau * uy,
        tau * (e_local + p_local) * utau * ueta, tau * rhob_local * utau};
    cell_sol = reconst_test.ReconstIt_shell(tau, tauq_vec, temp_grid);
    CHECK(
        cell_sol.e
        == doctest::Approx(e_local).epsilon(
            std::max(abs_tol, e_local * rel_tol)));
    CHECK(
        cell_sol.rhob
        == doctest::Approx(rhob_local)
               .epsilon(std::max(abs_tol, rhob_local * rel_tol)));
    CHECK(
        cell_sol.u[0]
        == doctest::Approx(utau).epsilon(std::max(abs_tol, utau * rel_tol)));
    CHECK(
        cell_sol.u[1]
        == doctest::Approx(ux).epsilon(std::max(abs_tol, ux * rel_tol)));
    CHECK(
        cell_sol.u[2]
        == doctest::Approx(uy).epsilon(std::max(abs_tol, uy * rel_tol)));
    CHECK(
        cell_sol.u[3]
        == doctest::Approx(ueta).epsilon(std::max(abs_tol, ueta * rel_tol)));

    temp_grid.u = {100.0, 0.0, 0.0, 0.0};
    tau = 1.0;
    e_local = 3.21e-13;
    rhob_local = 0.0;
    p_local = eos_ideal.get_pressure(e_local, rhob_local);
    ux = 5642.;
    uy = 65.4;
    ueta = 0.2;
    utau = sqrt(1. + ux * ux + uy * uy + ueta * ueta);
    tauq_vec = {
        tau * ((e_local + p_local) * utau * utau - p_local),
        tau * (e_local + p_local) * utau * ux,
        tau * (e_local + p_local) * utau * uy,
        tau * (e_local + p_local) * utau * ueta, tau * rhob_local * utau};
    cell_sol = reconst_test.ReconstIt_shell(tau, tauq_vec, temp_grid);
    CHECK(cell_sol.e == doctest::Approx(e_local).epsilon(1e-20));
    CHECK(cell_sol.rhob == doctest::Approx(rhob_local).epsilon(1e-20));
    CHECK(cell_sol.u[0] == doctest::Approx(utau).epsilon(utau * 1e-11));
    CHECK(cell_sol.u[1] == doctest::Approx(ux).epsilon(ux * 1e-11));
    CHECK(cell_sol.u[2] == doctest::Approx(uy).epsilon(uy * 1e-9));
    CHECK(cell_sol.u[3] == doctest::Approx(ueta).epsilon(ueta * 2e-8));
}

TEST_CASE("A regulated cell keeps its net baryon density") {
    EOS eos_ideal(0);
    Reconst reconst_test(eos_ideal, 9, 0);
    // T^{tau tau} < 0: the cell is reset to rest at e = abs_err; its
    // J^tau stays, so N_B is conserved
    const double tau = 1.8;
    const double nB = 0.3;
    TJbVec tauq_vec = {-tau * 1e-3, tau * 2e-3, 0., 0., tau * nB};
    ReconstCell prev;
    prev.e = 1e-6;
    prev.rhob = 0.05;
    prev.u = {1.0, 0.0, 0.0, 0.0};
    auto cell = reconst_test.ReconstIt_shell(tau, tauq_vec, prev);
    CHECK(cell.e <= 1e-10);
    CHECK(cell.u[0] == 1.);
    CHECK(cell.u[1] == 0.);
    CHECK(cell.rhob == doctest::Approx(nB));
}
