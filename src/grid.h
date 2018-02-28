#ifndef _SRC_GRID_H_
#define _SRC_GRID_H_

#include <cassert>
#include <vector>
#include "cell.h"
#include "./grid.h"

class Grid {
 private:
  std::vector<Cell> grid;
  int Nx   = 0;
  int Ny   = 0;
  int Neta = 0;
  
  Cell& get(int x, int y, int eta);
  
 public:
  Grid() = default;
  Grid(int Nx0, int Ny0, int Neta0);

  int nX()   const;
  int nY()   const;
  int nEta() const;

  const Cell& operator()(int x, int y, int eta) const {
    return const_cast<Cell&>(static_cast<const Cell &>((*this)(x,y,eta)));
  }

  Cell& operator()(int x, int y, int eta) {
    x   += 2;  assert(x>=0);   assert(x<Nx);
    y   += 2;  assert(y>=0);   assert(y<Ny);
    eta += 2;  assert(eta>=0); assert(eta<Neta);
    return grid[Nx*(Ny*eta+y)+x];
  }

  void updateHalo();

};

template<class Func>
void Neighbourloop(Grid &arena, int cx, int cy, int ceta, Func func){
  const std::array<int,6> dx   = {-1,1, 0,0, 0,0};
  const std::array<int,6> dy   = { 0,0,-1,1, 0,0};
  const std::array<int,6> deta = { 0,0, 0,0,-1,1};

  for(int dir=0;dir<3;dir++){
    const int m1nx   = dx  [2*dir];
    const int m1ny   = dy  [2*dir];
    const int m1neta = deta[2*dir];
    const int p1nx   = dx  [2*dir+1];
    const int p1ny   = dy  [2*dir+1];
    const int p1neta = deta[2*dir+1];
    const int m2nx   = 2*m1nx;  
    const int m2ny   = 2*m1ny;  
    const int m2neta = 2*m1neta;
    const int p2nx   = 2*p1nx;  
    const int p2ny   = 2*p1ny;  
    const int p2neta = 2*p1neta;
          auto&  c   = arena(cx,      cy,      ceta       );
    const auto& p1   = arena(cx+p1nx, cy+p1ny, ceta+p1neta);
    const auto& p2   = arena(cx+p2nx, cy+p2ny, ceta+p2neta);
    const auto& m1   = arena(cx+m1nx, cy+m1ny, ceta+m1neta);
    const auto& m2   = arena(cx+m2nx, cy+m2ny, ceta+m2neta);
    func(c,p1,p2,m1,m2,dir+1);
  }
}

#define NLAMBDA [&](Cell& c, const Cell& p1, const Cell& p2, const Cell& m1, const Cell& m2, const int direction) 

#endif
