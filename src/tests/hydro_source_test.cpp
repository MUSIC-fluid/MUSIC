// Copyright 2017 Chun Shen
// Test hydro_source class

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "./hydro_source.h"
#include "./read_in_parameters.h"
#include "./util.h"
#include "./data.h"

using namespace std;

int main(int argc, char *argv[]) {
    string input_file;
    InitData DATA;
    
    if (argc > 1) {
        input_file = *(argv+1);
    } else {
        cout << "usage: hydro_source_test inputfile" << endl;
        exit(0);
    }

    ReadInParameters reader;
    reader.read_in_parameters(&DATA, input_file);
    hydro_source *hydro_source_ptr = new hydro_source(&DATA);

    double ***ed = new double** [64];
    double ***rhob = new double** [64];
    for (int ieta = 0; ieta < 64; ieta++) {
        ed[ieta] = new double* [201];
        rhob[ieta] = new double* [201];
        for (int ix = 0; ix < 201; ix++) {
            ed[ieta][ix] = new double[201];
            rhob[ieta][ix] = new double[201];
            for (int iy = 0; iy < 201; iy++) {
                ed[ieta][ix][iy] = 0.0;
                rhob[ieta][ix][iy] = 0.0;
            }
        }
    }

    double dtau = 0.01;
    double dx = 0.1;
    double dy = 0.1;
    double deta = 0.2;
    double *u = new double[4];
    u[0] = 1.0;
    u[1] = 0.0;
    u[2] = 0.0;
    u[3] = 0.0;
    double *j_mu = new double[4];
    for (int i = 0; i < 3; i++) {
        j_mu[i] = 0.0;
    }
    double tau0 = 0.4;
    for (int itau = 0; itau < 401; itau++) {
        double tau_local = tau0 + itau*dtau;
        double total_e = 0.0;
        double total_rhob = 0.0;
        for (int ieta = 0; ieta < 64; ieta++) {
            double eta_local = -6.4 + ieta*deta;
            for (int ix = 0; ix < 201; ix++) {
                double x_local = -10. + ix*dx;
                for (int iy = 0; iy < 201; iy++) {
                    double y_local = -10. + iy*dy;
                    hydro_source_ptr->get_hydro_energy_source(
                            tau_local, x_local, y_local, eta_local, u, j_mu);
                    ed[ieta][ix][iy] += j_mu[0]*tau_local;
                    rhob[ieta][ix][iy] += (
                        hydro_source_ptr->get_hydro_rhob_source(
                            tau_local, x_local, y_local, eta_local, u)
                        *tau_local);
                    total_e += ed[ieta][ix][iy]*cosh(eta_local);
                    total_rhob += rhob[ieta][ix][iy];
                }
            }
        }
        total_e *= dx*dy*dtau*deta;
        total_rhob *= dx*dy*dtau*deta;
        cout << "tau = " << tau_local << ", total_e = " << total_e
             << ", total_rhob = " << total_rhob << endl;
    }

    ofstream of("check_rapidity.dat");
    for (int ieta = 0; ieta < 64; ieta++) {
        double eta_local = -6.4 + ieta*deta;
        double total_e = 0.0;
        double total_rhob = 0.0;
        for (int ix = 0; ix < 201; ix++) {
            for (int iy = 0; iy < 201; iy++) {
                total_e += ed[ieta][ix][iy]*cosh(eta_local);
                total_rhob += rhob[ieta][ix][iy];
            }
        }
        of << scientific << setw(18) << setprecision(8)
           << eta_local << "  " << total_e*dx*dy*dtau*deta << "  "
           << total_rhob*dx*dy*dtau*deta << endl;
    }
    of.close();

    delete[] u;
    delete[] j_mu;
    for (int ieta = 0; ieta < 64; ieta++) {
        for (int ix = 0; ix < 201; ix++) {
            delete[] ed[ieta][ix];
            delete[] rhob[ieta][ix];
        }
        delete[] ed[ieta];
        delete[] rhob[ieta];
    }
    delete[] ed;
    delete[] rhob;
    delete hydro_source_ptr;
    return(0);
}

