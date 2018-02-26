#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "cell.h"

TEST_CASE("Test constructor and copy"){
  Cell cell_a;

  cell_a.epsilon = 3;

  auto cell_b = cell_a;

  CHECK(cell_a.epsilon==cell_b.epsilon);
}