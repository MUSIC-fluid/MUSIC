// Copyright 2018 @ Chun Shen

#include "eos_s95p.h"
#include "doctest.h"

#include <cassert>
#include <iostream>

TEST_CASE("test constructor") {
    EOS_s95p test;
    CHECK(test.get_EOS_id() == 2);
    CHECK(test.get_number_of_tables() == 0);
}


TEST_CASE("test initialize eos") {
    EOS_s95p test;
    test.initialize_eos(3);
    CHECK(test.get_EOS_id() == 3);
    CHECK(test.get_number_of_tables() == 7);
    CHECK(test.get_eps_max() == doctest::Approx(769.82367).epsilon(0.0001));
}


TEST_CASE("test speed of sound squared") {
    // speed of sound squared is a constant 1/3
    EOS_s95p test;
    test.initialize_eos(7);
    CHECK(test.get_EOS_id() == 7);
    CHECK(test.get_cs2(1.0, 0.0) == doctest::Approx(0.15963).epsilon(0.0001));
    CHECK(test.get_cs2(1.0, 2.0) == doctest::Approx(0.15963).epsilon(0.0001));
    CHECK(test.get_cs2(1.0, -3.0) == doctest::Approx(0.15963).epsilon(0.0001));

}


TEST_CASE("test dPdrho") {
    // dPdrho = 0.
    EOS_s95p test;
    test.initialize_eos(2);
    CHECK(test.get_EOS_id() == 2);
    CHECK(test.p_rho_func(1.0, 1.0)  == 0.);
    CHECK(test.p_rho_func(2.0, 1.0)  == 0.);
    CHECK(test.p_rho_func(1.0, -1.0) == 0.);
}


TEST_CASE("test dPde") {
    EOS_s95p test;
}


TEST_CASE("test muB") {
    // muB = 0.0
    EOS_s95p test;
    test.initialize_eos(4);
    CHECK(test.get_EOS_id() == 4);
    CHECK(test.get_mu(1.0, 1.0)  == 0.);
    CHECK(test.get_mu(2.0, 1.0)  == 0.);
    CHECK(test.get_mu(1.0, -1.0) == 0.);
}


TEST_CASE("test muS") {
    // muS = 0.0
    EOS_s95p test;
    test.initialize_eos(5);
    CHECK(test.get_EOS_id() == 5);
    CHECK(test.get_muS(1.0, 1.0)  == 0.);
    CHECK(test.get_muS(2.0, 1.0)  == 0.);
    CHECK(test.get_muS(1.0, -1.0) == 0.);
}


TEST_CASE("test pressure") {
    EOS_s95p test;
}


TEST_CASE("check eos") {
    EOS_s95p test;
    test.initialize_eos(2);
    test.check_eos();
}
