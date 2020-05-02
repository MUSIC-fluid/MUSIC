#include "cell.h"
#include "doctest.h"
#include <cassert>
#include <iostream>


TEST_CASE("Does cell + work") {
    Cell_small cell1;
    Cell_small cell2;
    cell1.epsilon = 1.0;
    cell2.epsilon = 1.0;
    Cell_small cell3 = cell1 + cell2;
    CHECK(cell3.epsilon == (cell1.epsilon + cell2.epsilon));
}
