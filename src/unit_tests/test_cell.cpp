
#include <iostream>
#include "../cell.h"


int main() {
    Cell test_cell_1;
    Cell test_cell_2;

    // test copy
    test_cell_1.epsilon = 0.1;
    test_cell_1.u_rk0[0] = 1.;
    test_cell_1.dUsup[1][3] = 2.;
    
    test_cell_2 = test_cell_1;
    std::cout << "Cell1: " << test_cell_1.epsilon
              << ", Cell2: " << test_cell_2.epsilon << std::endl;
    std::cout << "Cell1: " << test_cell_1.u_rk0[0]
              << ", Cell2: " << test_cell_2.u_rk0[0] << std::endl;
    std::cout << "Cell1: " << test_cell_1.dUsup[1][3]
              << ", Cell2: " << test_cell_2.dUsup[1][3] << std::endl;
}

