#ifndef INIT_H
#define INIT_H

#include <stdio.h>
#include <iostream>
#include "data.h"
#include "grid.h"
#include "glauber.h"
#include <vector>
#include <time.h>
#include <cmath>

class Init{
 private:
  Random *random;
  EOS *eos;
  Util * util;
  Glauber *glauber;

  vector<ReturnValue> nucleusA;  // list of x and y coordinates of nucleons in nucleus A      
  vector<ReturnValue> nucleusB;  // list of x and y coordinates of nucleons in nucleus B 
  
 public:
  Init(EOS *eos, Glauber* glauber);//constructor
  ~Init();//destructor

  void sampleTA();
  void InitArena(InitData *DATA, Grid ****arena, Grid ****Lneighbor, Grid ****Rneighbor, int size, int rank);
  void LinkNeighbors(InitData *DATA, Grid ****arena, int size, int rank);
  int InitTJb(InitData *DATA, Grid ****arena, Grid ****Lneighbor, Grid ****Rneighbor, int size, int rank);
  // The following two functions return T_A and T_B. 
  // Normalization: \int r T_A(r) dr dphi = A
  double TATarget(InitData *DATA, double r);
  double TAProjectile(InitData *DATA, double r);
  double eta_profile_normalisation(InitData *DATA, double eta);
};
#endif
