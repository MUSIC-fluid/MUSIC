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


class Cell_aux {
 public:
    VorticityVec omega_kSP = {0.};
    VorticityVec omega_k = {0.};
    VorticityVec omega_th = {0.};
    VorticityVec omega_T = {0.};
    VelocityShearVec sigma_th = {0.};
    DmuMuBoverTVec DbetaMu = {0.};


    Cell_aux operator + (Cell_aux const &obj) {
        Cell_aux res;
        for (unsigned int i = 0; i < omega_k.size(); i++) {
            res.omega_kSP[i] = omega_kSP[i] + obj.omega_kSP[i];
            res.omega_k[i] = omega_k[i] + obj.omega_k[i];
            res.omega_th[i] = omega_th[i] + obj.omega_th[i];
            res.omega_T[i] = omega_T[i] + obj.omega_T[i];
        }
        for (unsigned int i = 0; i < sigma_th.size(); i++)
            res.sigma_th[i] = sigma_th[i] + obj.sigma_th[i];
        for (unsigned int i = 0; i < DbetaMu.size(); i++)
            res.DbetaMu[i] = DbetaMu[i] + obj.DbetaMu[i];
        return(res);
    }


    Cell_aux operator * (const double a) {
        Cell_aux res;
        for (unsigned int i = 0; i < omega_k.size(); i++) {
            res.omega_kSP[i] = omega_kSP[i]*a;
            res.omega_k[i] = omega_k[i]*a;
            res.omega_th[i] = omega_th[i]*a;
            res.omega_T[i] = omega_T[i]*a;
        }
        for (unsigned int i = 0; i < sigma_th.size(); i++)
            res.sigma_th[i] = sigma_th[i]*a;
        for (unsigned int i = 0; i < DbetaMu.size(); i++)
            res.DbetaMu[i] = DbetaMu[i]*a;
        return(res);
    }
};

#endif  // SRC_CELL_H_
