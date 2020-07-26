#ifndef _SRC_GRID_H_
#define _SRC_GRID_H_

#include <cassert>
#include <vector>
#include "cell.h"
#include "grid.h"

template<class T>
class GridT {
 private:
    std::vector<T> grid;

    int Nx   = 0;
    int Ny   = 0;
    int Neta = 0;

    T& get(int x, int y, int eta) {
        return grid[Nx*(Ny*eta+y)+x];
    }

 public:
    GridT() = default;
    GridT(int Nx0, int Ny0, int Neta0) {
        Nx   = Nx0  ;
        Ny   = Ny0  ;
        Neta = Neta0;
        grid.resize(Nx*Ny*Neta);
    }

    int nX()   const {return(Nx );  }
    int nY()   const {return(Ny );  }
    int nEta() const {return(Neta );}
    int size() const {return Nx*Ny*Neta;}

    T& getHalo(int x, int y, int eta){
        assert(-2<=x  ); assert(x  <Nx  +2);
        assert(-2<=y  ); assert(y  <Ny  +2);
        assert(-2<=eta); assert(eta<Neta+2);
        if(x  <0)   x  =0;  else if(x  >=Nx)   x  = Nx   - 1;
        if(y  <0)   y  =0;  else if(y  >=Ny)   y  = Ny   - 1;
        if(eta<0)   eta=0;  else if(eta>=Neta) eta= Neta - 1;
        return get(x,y,eta);
    }

    const T& getHalo(int x, int y, int eta) const {
        return getHalo(x, y, eta);
    }

    T& operator()(const int x, const int y, const int eta) {
        assert(0<=x  ); assert(x  <Nx);
        assert(0<=y  ); assert(y  <Ny);
        assert(0<=eta); assert(eta<Neta);
        return get(x, y, eta);
    }

    const T& operator()(int x, int y, int eta) const {
        assert(0<=x  ); assert(x  <Nx);
        assert(0<=y  ); assert(y  <Ny);
        assert(0<=eta); assert(eta<Neta);
        return get(x, y, eta);
    }

    T& operator()(const int i) {
        assert(0<=i  ); assert(i<Nx*Ny*Neta);
        return grid[i];
    }

    const T& operator()(const int i) const {
        assert(0<=i  ); assert(i<Nx*Ny*Neta);
        return grid[i];
    }

    void clear() {
        grid.clear();
        grid.shrink_to_fit();
    }
};

typedef GridT<Cell_small> SCGrid;

template<class T, class Func>
void Neighbourloop(GridT<T> &arena, int cx, int cy, int ceta, Func func) {
    const std::array<int, 6> dx   = {-1, 1,  0, 0,  0, 0};
    const std::array<int, 6> dy   = { 0, 0, -1, 1,  0, 0};
    const std::array<int, 6> deta = { 0, 0,  0, 0, -1, 1};
    for(int dir = 0; dir < 3; dir++) {
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
              auto&  c   = arena        (cx,      cy,      ceta       );
        const auto& p1   = arena.getHalo(cx+p1nx, cy+p1ny, ceta+p1neta);
        const auto& p2   = arena.getHalo(cx+p2nx, cy+p2ny, ceta+p2neta);
        const auto& m1   = arena.getHalo(cx+m1nx, cy+m1ny, ceta+m1neta);
        const auto& m2   = arena.getHalo(cx+m2nx, cy+m2ny, ceta+m2neta);
        func(c,p1,p2,m1,m2,dir+1);
    }
}

#define NLAMBDAS [&](Cell_small& c, const Cell_small& p1, const Cell_small& p2, const Cell_small& m1, const Cell_small& m2, const int direction) 

#endif
