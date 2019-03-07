// Copyright @ Chun Shen
#include "doctest.h"
#include "read_in_parameters.h"
#include "critical_modes.h"
#include "grid.h"

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
