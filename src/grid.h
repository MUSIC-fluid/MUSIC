#ifndef _SRC_GRID_H_
#define _SRC_GRID_H_

#include <vector>
#include "cell.h"

class Grid {
 private:
  std::vector<Cell> grid;
  int Nx;
  int Ny;
  int Neta;
  
 public:
  Grid(int Nx0, int Ny0, int Neta0);

  Cell& operator()(int x, int y, int eta);
  const Cell& operator()(int x, int y, int eta) const;

  void updateHalo();

};

#endif
