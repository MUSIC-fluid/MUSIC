// Copyright 2018 @ Chun Shen
#include "eos_idealgas.h"
#include "doctest.h"

#include <cassert>
#include <iostream>

TEST_CASE("test constructor") {
    EOS_idealgas test;
    CHECK(test.get_EOS_id() == 0);
    CHECK(test.get_number_of_tables() == 0);
}


TEST_CASE("test speed of sound squared") {
    // speed of sound squared is a constant 1/3
    EOS_idealgas test;
    CHECK(test.get_cs2(1.0, 1.0)  == 1./3.);
    CHECK(test.get_cs2(2.0, 1.0)  == 1./3.);
    CHECK(test.get_cs2(1.0, -1.0) == 1./3.);
}


TEST_CASE("test dPdrho") {
    // dPdrho = 0.
    EOS_idealgas test;
    CHECK(test.p_rho_func(1.0, 1.0)  == 0.);
    CHECK(test.p_rho_func(2.0, 1.0)  == 0.);
    CHECK(test.p_rho_func(1.0, -1.0) == 0.);
}


TEST_CASE("test dPde") {
    // dPdrho = 1/3.
    EOS_idealgas test;
    CHECK(test.p_e_func(1.0, 1.0)  == 1./3.);
    CHECK(test.p_e_func(2.0, 1.0)  == 1./3.);
    CHECK(test.p_e_func(1.0, -1.0) == 1./3.);
}


TEST_CASE("test muB") {
    // muB = 0.0
    EOS_idealgas test;
    CHECK(test.get_mu(1.0, 1.0)  == 0.);
    CHECK(test.get_mu(2.0, 1.0)  == 0.);
    CHECK(test.get_mu(1.0, -1.0) == 0.);
}


TEST_CASE("test muS") {
    // muS = 0.0
    EOS_idealgas test;
    CHECK(test.get_muS(1.0, 1.0)  == 0.);
    CHECK(test.get_muS(2.0, 1.0)  == 0.);
    CHECK(test.get_muS(1.0, -1.0) == 0.);
}


TEST_CASE("test pressure") {
    // P = e/3.
    EOS_idealgas test;
    CHECK(test.get_pressure(1.0, 1.0)  == 1./3.);
    CHECK(test.get_pressure(2.0, 1.0)  == 2./3.);
    CHECK(test.get_pressure(1.0, -1.0) == 1./3.);
}
