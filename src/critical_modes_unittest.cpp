// Copyright @ Chun Shen
#include "doctest.h"
#include "read_in_parameters.h"
#include "critical_modes.h"

TEST_CASE("Check CriticalSlowModes initialization") {
    EOS eos_ideal(0);
    InitData DATA(ReadInParameters::read_in_parameters(
                            "tests/unittest_files/music_input_criticalmodes"));
    CriticalSlowModes test(eos_ideal, DATA);
    test.InitializeFields(10, 3, 4, 5);
    CHECK(test.get_Qvec_size() == 10);
    CHECK(test.get_nphiQfields() == 10);
    test.InitializeFields(20, 5, 10, 3);
    CHECK(test.get_Qvec_size() == 20);
    CHECK(test.get_nphiQfields() == 20);
}
