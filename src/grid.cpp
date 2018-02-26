#include "grid.h"
#include "doctest.h"

Grid::Grid(int Nx0, int Ny0, int Neta0){
  Nx   = Nx0  +4;
  Ny   = Ny0  +4;
  Neta = Neta0+4;
  grid.resize(Nx*Ny*Neta);
}

Cell& Grid::operator()(int x, int y, int eta)
{
  x   += 2;
  y   += 2;
  eta += 2;
  return grid[Nx*(Ny*eta+y)+x];
}

const Cell& Grid::operator()(int x, int y, int eta) const
{
  x   += 2;
  y   += 2;
  eta += 2;
  return grid[Nx*(Ny*eta+y)+x];
}

void Grid::updateHalo()
{
  auto self = *this;

  //x-y planes
  auto xycopy = [&](int from,int to){
    for(int y=0;y<Ny;y++)
    for(int x=0;x<Nx;x++)
      self(x,y,to)=self(x,y,from);
  };
  
  xycopy(2,1);
  xycopy(2,0);
  xycopy(Neta-3,Neta-2);
  xycopy(Neta-3,Neta-1);
  
  //x-eta planes
  auto xetacopy = [&](int from,int to){
    for(int eta=0;eta<Neta;eta++)
    for(int x  =0;x  <Nx;  x++  )
      self(x,to,eta)=self(x,from,eta);
  };
  
  xetacopy(2,1);
  xetacopy(2,0);
  xetacopy(Ny-3,Ny-2);
  xetacopy(Ny-3,Ny-1);
  
  //y-eta planes
  auto yetacopy = [&](int from,int to){
    for(int eta=0;eta<Neta;eta++)
    for(int y  =0;y  <Ny;  y++  )
      self(to,y,eta)=self(from,y,eta);
  };
  
  yetacopy(2,1);
  yetacopy(2,0);
  yetacopy(Nx-3,Nx-2);
  yetacopy(Nx-3,Nx-1); 
}



TEST_CASE("Does grid copy work"){
  Grid grid(3,3,3);

  grid(0,0,0).epsilon = 3;

  auto grid2 = grid;

  CHECK(grid2(0,0,0).epsilon == grid(0,0,0).epsilon);
}