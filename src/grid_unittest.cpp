#include "grid.h"
#include "doctest.h"
#include <cassert>
#include <iostream>

TEST_CASE("Does grid copy work"){
    SCGrid grid(3, 3, 3);
    grid(0, 0, 0).epsilon = 3;
    auto grid2 = grid;
    CHECK(grid2(0, 0, 0).epsilon == grid(0, 0, 0).epsilon);
}

TEST_CASE("check neighbourloop1") {
    SCGrid grid(1, 1, 1);
    grid(0,0,0).epsilon = 3;
    Neighbourloop(grid, 0, 0, 0, NLAMBDAS {
        CHECK(c.epsilon == p1.epsilon);
        CHECK(c.epsilon == p2.epsilon);
        CHECK(c.epsilon == m1.epsilon);
        CHECK(c.epsilon == m2.epsilon);
    });
}
 
TEST_CASE("check neighbourloop2") {
    SCGrid grid1(5, 1, 1);
    for (int i = 0; i < 5; i++) {
        grid1(i, 0, 0).epsilon = i + 1;
    }
    Neighbourloop(grid1, 2, 0, 0, NLAMBDAS {
        if (direction == 1) {
            CHECK(p1.epsilon == 4);
            CHECK(p2.epsilon == 5);
            CHECK(m1.epsilon == 2);
            CHECK(m2.epsilon == 1);
        } else {
            CHECK(p1.epsilon == c.epsilon);
            CHECK(p2.epsilon == c.epsilon);
            CHECK(m1.epsilon == c.epsilon);
            CHECK(m2.epsilon == c.epsilon);
        }
    });
}
 
TEST_CASE("check neighbourloop3") {
    SCGrid grid2(3, 3, 3);
    for (int i = 0; i < 3; i++) 
    for (int j = 0; j < 3; j++) 
    for (int k = 0; k < 3; k++) {
        grid2(i, j, k).epsilon = 1;
    }

    Neighbourloop(grid2, 1, 1, 1, NLAMBDAS {
        int sum = 0;
        sum += c.epsilon;
        sum += p1.epsilon;
        sum += p2.epsilon;
        sum += m1.epsilon;
        sum += m2.epsilon;
        CHECK(sum == 5);
    });
}

TEST_CASE("check neighbourloop4"){
    SCGrid grid(1, 1, 1);
    grid(0, 0, 0).epsilon = 1;

    int sum = 0;
    Neighbourloop(grid, 0, 0, 0, NLAMBDAS {
      sum += c.epsilon;
      sum += p1.epsilon;
      sum += p2.epsilon;
      sum += m1.epsilon;
      sum += m2.epsilon;
    });
    CHECK(sum == 15);
}

TEST_CASE("check dimension"){
    SCGrid grid(1, 2,3);

    CHECK(grid.nX()   == 1);
    CHECK(grid.nY()   == 2);
    CHECK(grid.nEta() == 3);
}

