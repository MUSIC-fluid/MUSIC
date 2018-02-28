#include "grid.h"
#include "doctest.h"
#include <cassert>
#include <iostream>

Grid::Grid(int Nx0, int Ny0, int Neta0){
  Nx   = Nx0  +4;
  Ny   = Ny0  +4;
  Neta = Neta0+4;
  grid.resize(Nx*Ny*Neta);
}

const Cell& Grid::operator()(int x, int y, int eta) const
{
  return const_cast<Cell&>(static_cast<const Cell &>((*this)(x,y,eta)));
}

Cell& Grid::operator()(int x, int y, int eta)
{
  x   += 2;  assert(x>=0);   assert(x<Nx);
  y   += 2;  assert(y>=0);   assert(y<Ny);
  eta += 2;  assert(eta>=0); assert(eta<Neta);
  return grid[Nx*(Ny*eta+y)+x];
}

Cell& Grid::get(int x, int y, int eta){
  return grid[Nx*(Ny*eta+y)+x];
}

int Grid::nX() const {
  return(Nx - 4);
}

int Grid::nY() const {
  return(Ny - 4);
}

int Grid::nEta() const {
  return(Neta - 4);
}

void Grid::updateHalo()
{
  //x-y planes
  auto xycopy = [&](int from,int to){
    for(int y=0;y<Ny;y++)
    for(int x=0;x<Nx;x++)
      get(x,y,to)=get(x,y,from);
  };
  
  xycopy(2,1);
  xycopy(2,0);
  xycopy(Neta-3,Neta-2);
  xycopy(Neta-3,Neta-1);
  
  //x-eta planes
  auto xetacopy = [&](int from,int to){
    for(int eta=0;eta<Neta;eta++)
    for(int x  =0;x  <Nx;  x++  )
      get(x,to,eta)=get(x,from,eta);
  };
  
  xetacopy(2,1);
  xetacopy(2,0);
  xetacopy(Ny-3,Ny-2);
  xetacopy(Ny-3,Ny-1);
  
  //y-eta planes
  auto yetacopy = [&](int from,int to){
    for(int eta=0;eta<Neta;eta++)
    for(int y  =0;y  <Ny;  y++  )
      get(to,y,eta)=get(from,y,eta);
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

TEST_CASE("Check halo"){
 Grid grid(1,1,1);

 grid(0,0,0).epsilon = 3;
 grid.updateHalo();

 for(int i=-2; i<3; i++)
   for(int j=-2; j<3; j++)
     for(int k=-2; k<3; k++)
	{
	  CHECK(grid(i,j,k).epsilon == grid(0,0,0).epsilon);
	}
}

TEST_CASE("check neighbourloop1"){
 Grid grid(1,1,1);

 grid(0,0,0).epsilon = 3;
 grid.updateHalo();

 Neighbourloop(grid, 0, 0, 0, NLAMBDA{
    CHECK(c.epsilon == p1.epsilon);
    CHECK(c.epsilon == p2.epsilon);
    CHECK(c.epsilon == m1.epsilon);
    CHECK(c.epsilon == m2.epsilon);
  });
}
 
TEST_CASE("check neighbourloop2"){
 Grid grid1(5,1,1);
 for (int i = 0; i < 5; i++) {
    grid1(i,0,0).epsilon = i + 1;
 }
 grid1.updateHalo();
 
 Neighbourloop(grid1, 2, 0, 0, NLAMBDA{
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
 
TEST_CASE("check neighbourloop3"){
  Grid grid2(3,3,3);
  for (int i = 0; i < 3; i++) 
  for (int j = 0; j < 3; j++) 
  for (int k = 0; k < 3; k++) {
    grid2(i,j,k).epsilon = 1;
  }
  grid2.updateHalo();

  Neighbourloop(grid2, 1, 1, 1, NLAMBDA{
    int sum = 0;
    sum += c.epsilon;
    sum += p1.epsilon;
    sum += p2.epsilon;
    sum += m1.epsilon;
    sum += m2.epsilon;
    CHECK(sum == 5);
  });

  int a = 0;
}

TEST_CASE("check neighbourloop4"){
  Grid grid(1,1,1);
  grid(0,0,0).epsilon = 1;
  grid.updateHalo();

  int sum = 0;
  Neighbourloop(grid, 0, 0, 0, NLAMBDA{
    sum += c.epsilon;
    sum += p1.epsilon;
    sum += p2.epsilon;
    sum += m1.epsilon;
    sum += m2.epsilon;
  });
  CHECK(sum == 15);
}

TEST_CASE("check dimension"){
  Grid grid(1,2,3);

  CHECK(grid.nX()   == 1);
  CHECK(grid.nY()   == 2);
  CHECK(grid.nEta() == 3);
}

