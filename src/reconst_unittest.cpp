
#include <iostream>
#include "eos.h"
#include "reconst.h"
#include "doctest.h"

TEST_CASE("Check Reconst constructor") {
    EOS eos_ideal(0);
    Reconst reconst_test(eos_ideal, 9);
    CHECK(reconst_test.get_max_iter() == 100);
    CHECK(reconst_test.get_echo_level() == 9);
    CHECK(reconst_test.get_abs_err() == doctest::Approx(1e-16).epsilon(1e-16));
    CHECK(reconst_test.get_v_critical() == doctest::Approx(0.563624).epsilon(0.00001));
}

TEST_CASE("Test Newton solver v") {
    EOS eos_ideal(0);
    Reconst reconst_test(eos_ideal, 9);
    double v_sol = 10.0;
    reconst_test.solve_velocity_Newton(0.0, 1.0, 0.0, 0.0, v_sol);
    CHECK(v_sol == doctest::Approx(0.0).epsilon(1e-16));

    double M0 = 1.0;
    double M = 0.5;
    double v_correct = 3.*M/(2.*M0 + sqrt(4.*M0*M0 - 3.*M*M));
    std::cout << "check v = " << v_correct << std::endl;
    reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-16));
    
    M0 = 1.0;
    M = 0.99;
    v_correct = 3.*M/(2.*M0 + sqrt(4.*M0*M0 - 3.*M*M));
    std::cout << "check v = " << v_correct << std::endl;
    reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-16));
    
    M0 = 1.0;
    M = 1.0 - 1e-6;
    v_correct = 3.*M/(2.*M0 + sqrt(4.*M0*M0 - 3.*M*M));
    std::cout << "check v = " << v_correct << std::endl;
    reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-16));
    
    M0 = 1.0;
    M = 1e-6;
    v_correct = 3.*M/(2.*M0 + sqrt(4.*M0*M0 - 3.*M*M));
    std::cout << "check v = " << v_correct << std::endl;
    reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-16));
    
    M0 = 1.0;
    M = 1.0 - 1e-10;
    v_correct = 3.*M/(2.*M0 + sqrt(4.*M0*M0 - 3.*M*M));
    std::cout << "check v = " << v_correct << std::endl;
    reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-16));
    
    M0 = 1.0;
    M = 1e-10;
    v_correct = 3.*M/(2.*M0 + sqrt(4.*M0*M0 - 3.*M*M));
    std::cout << "check v = " << v_correct << std::endl;
    reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-16));
    
    M0 = 1e-8;
    M = 1e-12;
    v_correct = 3.*M/(2.*M0 + sqrt(4.*M0*M0 - 3.*M*M));
    std::cout << "check v = " << v_correct << std::endl;
    reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-16));
    
    M0 = 1e-12;
    M = 3.23e-15;
    v_correct = 3.*M/(2.*M0 + sqrt(4.*M0*M0 - 3.*M*M));
    std::cout << "check v = " << v_correct << std::endl;
    reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-16));
    
    M0 = 1e-14;
    M = M0 - 2.43e-18;
    v_correct = 3.*M/(2.*M0 + sqrt(4.*M0*M0 - 3.*M*M));
    std::cout << "check v = " << v_correct << std::endl;
    reconst_test.solve_velocity_Newton(0.0, M0, M, 0.0, v_sol);
    CHECK(v_sol == doctest::Approx(v_correct).epsilon(1e-16));
}

TEST_CASE("Test Newton solver u0") {
    EOS eos_ideal(0);
    Reconst reconst_test(eos_ideal, 9);
    
    double u0_sol = 0.0;
    double M0 = 1.0;
    double M  = 0.0;
    double v_correct = 3.*M/(2.*M0 + sqrt(4.*M0*M0 - 3.*M*M));
    double u0_correct = 1./sqrt(1. - v_correct*v_correct);
    double utol = u0_correct*1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    reconst_test.solve_u0_Newton(100.0, M0, M*M, M, 0.0, u0_sol);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));
    
    M0 = 1.0;
    M  = 0.5;
    v_correct = 3.*M/(2.*M0 + sqrt(4.*M0*M0 - 3.*M*M));
    u0_correct = 1./sqrt(1. - v_correct*v_correct);
    utol = u0_correct*1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    reconst_test.solve_u0_Newton(100.0, M0, M*M, M, 0.0, u0_sol);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));
    
    M0 = 1.0;
    M  = 1e-4;
    v_correct = 3.*M/(2.*M0 + sqrt(4.*M0*M0 - 3.*M*M));
    u0_correct = 1./sqrt(1. - v_correct*v_correct);
    utol = u0_correct*1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    reconst_test.solve_u0_Newton(100.0, M0, M*M, M, 0.0, u0_sol);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));
    
    M0 = 1.0;
    M  = 0.999;
    v_correct = 3.*M/(2.*M0 + sqrt(4.*M0*M0 - 3.*M*M));
    u0_correct = 1./sqrt(1. - v_correct*v_correct);
    utol = u0_correct*1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    reconst_test.solve_u0_Newton(100.0, M0, M*M, M, 0.0, u0_sol);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));
    
    M0 = 1.0;
    M  = 1.0 - 1e-4;
    v_correct = 3.*M/(2.*M0 + sqrt(4.*M0*M0 - 3.*M*M));
    u0_correct = 1./sqrt(1. - v_correct*v_correct);
    utol = u0_correct*1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    reconst_test.solve_u0_Newton(100.0, M0, M*M, M, 0.0, u0_sol);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));
    
    M0 = 1e-15;
    M  = M0 - 3.463e-18;
    v_correct = 3.*M/(2.*M0 + sqrt(4.*M0*M0 - 3.*M*M));
    u0_correct = 1./sqrt(1. - v_correct*v_correct);
    utol = u0_correct*1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    reconst_test.solve_u0_Newton(100.0, M0, M*M, M, 0.0, u0_sol);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));
    
    M0 = 1e-15;
    M  = 3.463e-18;
    v_correct = 3.*M/(2.*M0 + sqrt(4.*M0*M0 - 3.*M*M));
    u0_correct = 1./sqrt(1. - v_correct*v_correct);
    utol = u0_correct*1e-14;
    std::cout << "check u0 = " << u0_correct << std::endl;
    reconst_test.solve_u0_Newton(100.0, M0, M*M, M, 0.0, u0_sol);
    CHECK(u0_sol == doctest::Approx(u0_correct).epsilon(utol));
}

