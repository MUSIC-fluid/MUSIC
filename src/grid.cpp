#include "Grid.h"

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
  //x-y planes
  auto xycopy = [&](int from,int to){
    for(int y=0;y<Ny;y++)
      for(int x=0;x<Nx;x++)
        grid(x,y,to)=grid(x,y,from);
  };
  
  xycopy(2,1);
  xycopy(2,0);
  xycopy(Neta-3,Neta-2);
  xycopy(Neta-3,Neta-1);
  
  //x-eta planes
  auto xetacopy = [&](int from,int to){
    for(int eta=0;eta<Neta;eta++)
      for(int x=0;x<Nx;x++)
        grid(x,to,eta)=grid(x,from,eta);
  };
  
  xetacopy(2,1);
  xetacopy(2,0);
  xetacopy(Ny-3,Ny-2);
  xetacopy(Ny-3,Ny-1);
  
  //y-eta planes
  auto xetacopy = [&](int from,int to){
    for(int eta=0;eta<Neta;eta++)
      for(int y=0;y<Ny;y++)
        grid(to,y,eta)=grid(from,y,eta);
  };
  
  xetacopy(2,1);
  xetacopy(2,0);
  xetacopy(Nx-3,Nx-2);
  xetacopy(Nx-3,Nx-1); 
}
