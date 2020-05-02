// Copyright 2017 Chun Shen, Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef _SRC_CELL_H_
#define _SRC_CELL_H_

#include "data_struct.h"
#include <array>
#include <cmath>

class Cell_small {
 public:
    double epsilon = 0;
    double rhob = 0;
    FlowVec u = {1., 0., 0., 0.};

    ViscousVec Wmunu = {0.};
    double pi_b = 0.;


    Cell_small operator + (Cell_small const &obj) {
        Cell_small res;
        res.epsilon = epsilon + obj.epsilon;
        res.rhob = rhob + obj.rhob;
        res.u[1] = u[1] + obj.u[1];
        res.u[2] = u[2] + obj.u[2];
        res.u[3] = u[3] + obj.u[3];
        res.u[0] = sqrt(1. + res.u[1]*res.u[1]
                        + res.u[2]*res.u[2] + res.u[3]*res.u[3]);
        for (unsigned int i = 0; i < Wmunu.size(); i++) {
            res.Wmunu[i] = Wmunu[i] + obj.Wmunu[i];
        }
        res.pi_b = pi_b + obj.pi_b;
        return(res);
    }


    Cell_small operator * (const double a) {
        Cell_small res;
        res.epsilon = epsilon*a;
        res.rhob = rhob*a;
        res.u[1] = u[1]*a;
        res.u[2] = u[2]*a;
        res.u[3] = u[3]*a;
        res.u[0] = sqrt(1. + res.u[1]*res.u[1]
                        + res.u[2]*res.u[2] + res.u[3]*res.u[3]);
        for (unsigned int i = 0; i < Wmunu.size(); i++) {
            res.Wmunu[i] = Wmunu[i]*a;
        }
        res.pi_b = pi_b*a;
        return(res);
    }
};

#endif  // SRC_CELL_H_
