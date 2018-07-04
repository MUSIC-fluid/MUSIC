// Copyright 2018 @ Chun Shen

#include "eos.h"
#include "doctest.h"

#include <cassert>
#include <iostream>

TEST_CASE("test constructor") {
    EOS test;
    CHECK(test.get_pressure(1.0, 0.0) == 1./3.);
}
