// Copyright 2018 @ Chun Shen

#include "eos.h"
#include "doctest.h"

#include <cassert>
#include <iostream>

TEST_CASE("test constructor") {
    EOS test(0);
    CHECK(test.get_pressure(1.0, 0.0) == 1./3.);
    CHECK(test.get_dpde(1.0, 0.0)     == 1./3.);
    CHECK(test.get_dpdrhob(1.0, 0.0)  == 0.0);
    CHECK(test.get_muB(1.0, 1.0)      == 0.0);
    CHECK(test.get_muS(1.0, -1.0)     == 0.0);
}
