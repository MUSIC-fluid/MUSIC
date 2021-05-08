#include "cell.h"
#include "doctest.h"
#include <cassert>
#include <iostream>


TEST_CASE("Does cell = work") {
    Cell_small cell1;
    cell1.epsilon = 1.0;
    cell1.rhob = 2.0;
    cell1.u[1] = 0.2;
    cell1.Wmunu[3] = 0.2;

    Cell_small cell2 = cell1;
    CHECK(cell2.epsilon == cell1.epsilon);
    CHECK(cell2.rhob == cell1.rhob);
    CHECK(cell2.u == cell1.u);
    CHECK(cell2.Wmunu == cell1.Wmunu);
    CHECK(cell2.pi_b == cell1.pi_b);
}


TEST_CASE("Does cell + work") {
    Cell_small cell1;
    Cell_small cell2;
    cell1.epsilon = 1.0;
    cell2.epsilon = 1.0;
    Cell_small cell3 = cell1 + cell2;
    CHECK(cell3.epsilon == (cell1.epsilon + cell2.epsilon));
    CHECK(cell3.u[0] == 1.0);
}


TEST_CASE("Does cell * work") {
    Cell_small cell1;
    cell1.epsilon = 1.0;
    cell1.rhob = 2.0;
    cell1.u[1] = 0.2;

    double factor = 2.;
    Cell_small cell2 = cell1*factor;
    CHECK(cell2.epsilon == cell1.epsilon*factor);
    CHECK(cell2.rhob == cell1.rhob*factor);
    CHECK(cell2.u[1] == cell1.u[1]*factor);
}
