// Copyright 2018 @ Chun Shen

#include "eos_wrapper.h"
#include "doctest.h"

#include <cassert>
#include <iostream>

TEST_CASE("test constructor") {
    EOS_wrapper test;
    CHECK(test.get_pressure(1.0, 0.0) == 1./3.);
}
