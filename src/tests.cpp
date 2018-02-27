#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include <iostream>
#include "cell.h"
#include "./grid.h"

TEST_CASE("Test constructor and copy"){
  Cell cell_a;

  cell_a.epsilon = 3;
  cell_a.u[1][3] = 1.0;

  auto cell_b = cell_a;

  CHECK(cell_a.epsilon==cell_b.epsilon);
}

TEST_CASE("Test the dimension of the 1D arrays"){
  Cell cell_a;

  CHECK(cell_a.prev_pi_b.size()==1);
}

TEST_CASE("Test the dimension of the 2D arrays"){
  Cell cell_a;

  CHECK(cell_a.u.size()==2);
  CHECK(cell_a.u[0].size()==4);
  CHECK(cell_a.prev_u.size()==1);
  CHECK(cell_a.prev_u[0].size()==4);
  CHECK(cell_a.dUsup.size()==5);
  CHECK(cell_a.dUsup[0].size()==4);
  CHECK(cell_a.Wmunu.size()==2);
  CHECK(cell_a.Wmunu[0].size()==14);
  CHECK(cell_a.prevWmunu.size()==1);
  CHECK(cell_a.prevWmunu[0].size()==14);
}
