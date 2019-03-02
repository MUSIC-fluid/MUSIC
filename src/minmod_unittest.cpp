#include "doctest.h"
#include "minmod.h"

TEST_CASE("Check Minmod theta") {
    Minmod test(1.5);
    auto theta_out = test.get_theta();
    CHECK(theta_out == 1.5);
}

TEST_CASE("Does Minmod work") {
    Minmod test(1.8);
    auto test_dx = test.minmod_dx(0.0, 1.0, 0.0);
    CHECK(test_dx == 0.0);
    test_dx = test.minmod_dx(2.0, 1.0, 0.0);
    CHECK(test_dx == 1.0);
    test_dx = test.minmod_dx(-2.0, -1.0, 0.0);
    CHECK(test_dx == -1.0);
    test_dx = test.minmod_dx(0.0, 0.1, 2.0);
    CHECK(test_dx == doctest::Approx(-0.18).epsilon(0.0001));
    test_dx = test.minmod_dx(0.0, 1.9, 2.0);
    CHECK(test_dx == doctest::Approx(-0.18).epsilon(0.0001));
}
