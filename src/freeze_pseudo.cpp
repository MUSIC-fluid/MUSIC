// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

// Copyright @ Bjoern Schenke, Sangyong Jeon, and Charles Gale

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

#include<sys/stat.h>
#include<iomanip>
#include<cstring>
#include "freeze.h"

using Util::hbarc;
using std::string;
using std::stringstream;
using std::ofstream;
using std::endl;

// read in thermal spectra from file to then perform resonance decays with them
// Must set verbose to 1 if you want particleMax to be set by this routine.
void Freeze::ReadSpectra_pseudo(InitData* DATA, int full, int verbose) {
    // read in thermal spectra from file:
    int number, iptmax, iphimax;
    double deltaeta, etamax, ptmin, ptmax;
    int ietamax;
    int pseudo_steps = 0;
    int ip = 0;
    if (verbose) {
        music_message.info("reading spectra");
    }
    // open particle information file:
    FILE *p_file;
    const char* p_name = "particleInformation.dat";
    const char* pf_name = "FparticleInformation.dat";
    if (full)
        p_file = fopen(pf_name, "r");
    else
        p_file = fopen(p_name, "r");
    int count;
    count = 0;
    if (verbose) {
        music_message << "NumberOfParticlesToInclude "
                      << DATA->NumberOfParticlesToInclude;
        music_message.flush("info");
    }
    // read particle information:
    while (fscanf(p_file,"%d %lf %d %lf %lf %d %d ",
                  &number, &etamax, &ietamax, &ptmin, &ptmax,
                  &iptmax, &iphimax) == 7) {
        count++;
        if (count >= DATA->NumberOfParticlesToInclude) {
            break;
        }
        pseudo_steps = ietamax - 1;
        ip = partid[MHALF+number];
        particleList[ip].ny = ietamax;
        particleList[ip].npt = iptmax;
        particleList[ip].nphi = iphimax;
        particleList[ip].phimin = 0;
        particleList[ip].phimax = 2*M_PI;
        particleList[ip].slope = 1;
        particleList[ip].ymax = etamax;
        deltaeta = 0.;
        if (pseudo_steps > 0) {
            deltaeta = 2*etamax/pseudo_steps;
        }
        particleList[ip].deltaY = deltaeta;
        for (int ipt = 0; ipt < iptmax; ipt++) {
            particleList[ip].pt[ipt] = (
                ptmin + (ptmax - ptmin)*pow(static_cast<double>(ipt), 2.)
                        /pow(static_cast<double>(iptmax - 1), 2.));
        }
        for (int ieta = 0; ieta <= pseudo_steps; ieta++) {
            particleList[ip].y[ieta] = ieta*deltaeta - etamax;
        }
    
    }
    particleMax = ip + 1;

    fclose(p_file);

    FILE *s_file;
    string s_name = "yptphiSpectra.dat";
    string sf_name = "FyptphiSpectra.dat";
    if (full)
        s_file = fopen(sf_name.c_str(), "r");
    else
        s_file = fopen(s_name.c_str(), "r");
  
    if (verbose) {
        music_message << "ietamax = " << pseudo_steps+1;
        music_message.flush("info");
        music_message << "cells = " << (pseudo_steps+1)*(iptmax)*iphimax;
        music_message.flush("info");
        music_message << "particleMax = " << particleMax;
        music_message.flush("info");
    }
    for (ip = 1; ip < particleMax; ip++) {
        if (verbose) {
            music_message << "reading particle " << ip << ": "
                          << particleList[ip].number << "  "
                          << particleList[ip].name << " ...";
            music_message.flush("info");
        }
        for(int ieta = 0; ieta <= pseudo_steps; ieta++) {
            for (int ipt = 0; ipt < iptmax; ipt++) {
                for (int iphi = 0; iphi < iphimax; iphi++) {
                    system_status_ = fscanf(
                           s_file, "%lf",
                           &particleList[ip].dNdydptdphi[ieta][ipt][iphi]);
                    if (particleList[ip].dNdydptdphi[ieta][ipt][iphi] < 0.) {
                        particleList[ip].dNdydptdphi[ieta][ipt][iphi] = 0;
                    }
                }
            }
        }
    }

    // phiArray for use in Edndp3 interpolation function
    // during resonance decay calculation
    if (phiArray == NULL) {
        phiArray = new double[iphimax];
        for (int iphi = 0; iphi < iphimax; iphi++) {
            phiArray[iphi] = iphi*2*M_PI/iphimax;
        }
    }
    fclose(s_file);
}

// Modified spectra calculation by ML 05/2013
// Calculates on fixed grid in pseudorapidity, pt, and phi
// adapted from ML and improved on performance (C. Shen 2015)
void Freeze::ComputeParticleSpectrum_pseudo_improved(InitData *DATA,
                                                     int number) {
    double y_minus_eta_cut = 4.0;
    int j = partid[MHALF+number];
    // set some parameters
    double etamax = DATA->max_pseudorapidity;
    int ietamax = DATA->pseudo_steps + 1;  // pseudo_steps is number of steps.
                                           // Including edges
                                           // number of points is steps + 1
    double deltaeta = 0;
    if (ietamax > 1) {
        deltaeta = 2.*etamax/DATA->pseudo_steps;
    }
    double ptmax = DATA->max_pt;
    double ptmin = DATA->min_pt;
    int iptmax = DATA->pt_steps+1;  // Number of points is iptmax + 1
                                    // (ipt goes from 0 to iptmax)
    int iphimax = DATA->phi_steps;  // number of points
                                    // (phi=2pi equal to phi=0)
    double deltaphi = 2*M_PI/iphimax;
  
    // Reuse rapidity variables (Need to reuse variable y
    // so resonance decay routine can be used as is.
    // Might as well misuse ymax and deltaY too)
    particleList[j].ymax = etamax; 
    particleList[j].deltaY = deltaeta;

    music_message << "Doing " << j << ": "
                  << particleList[j].name << "("
                  << particleList[j].number << ") ... ";
    music_message.flush("info");
 
    particleList[j].ny = DATA->pseudo_steps + 1;
    particleList[j].npt = iptmax;
    particleList[j].nphi = iphimax;
  
    // set particle properties
    double m = particleList[j].mass;
    int deg = particleList[j].degeneracy;
    int baryon = particleList[j].baryon;
    // int s = particleList[j].strange;
    // int c = particleList[j].charge;
    double mu_PCE = particleList[j].muAtFreezeOut;
    int sign;
    if (baryon == 0) {
        sign = -1.;
    } else {
        sign = 1.;
    }

    // open files to write
    FILE *d_file;
    char *buf = new char[10];
    sprintf (buf, "%d", 0);
    
    char *specString=new char[30];
    strcpy(specString, "yptphiSpectra");
    strcat(specString, buf);
    strcat(specString, ".dat");
    char* s_name = specString;
    FILE *s_file;
    s_file = fopen(s_name, "w");

    delete[] specString;
    delete[] buf;

    // -----------------------------------------------------------------------
    // write information in the particle information file
    const char* d_name = "particleInformation.dat";
    d_file = fopen(d_name, "a");
    fprintf(d_file, "%d %e %d %e %e %d %d \n",
            number,  etamax, DATA->pseudo_steps+1, ptmin, ptmax,
            iptmax, iphimax);

    // caching 
    double* cos_phi = new double [iphimax];
    double* sin_phi = new double [iphimax];
    for (int iphi = 0; iphi < iphimax; iphi++) {
        double phi_local = deltaphi*iphi;
        cos_phi[iphi] = cos(phi_local);
        sin_phi[iphi] = sin(phi_local);
    }
    double* pt_array = new double [iptmax];
    for (int ipt = 0; ipt < iptmax; ipt++) {
        double pt =  (ptmin + (ptmax - ptmin)
                               *pow(static_cast<double>(ipt), 2.)
                               /pow(static_cast<double>(iptmax - 1), 2.));
        pt_array[ipt] = pt;
        particleList[j].pt[ipt] = pt;
    }
    
    double *bulk_deltaf_coeffs = new double[3];
    if (DATA_ptr->turn_on_bulk == 1 && DATA_ptr->include_deltaf_bulk == 1) {
        if (bulk_deltaf_kind == 0) {
            bulk_deltaf_coeffs = new double[3];
        } else {
            bulk_deltaf_coeffs = new double[2];
        }
    }
      
    double alpha = 0.0;
    // main loop begins ...
    // store E dN/d^3p as function of phi,
    // pt and eta (pseudorapidity) in sumPtPhi:
    int ieta;
    for (ieta = 0; ieta < ietamax; ieta++) {
        double eta = -etamax + ieta*deltaeta;

        // Use this variable to store pseudorapidity instead of rapidity
        // May cause confusion in the future,
        // but easier to to share code for both options:
        // calculating on a fixed grid in rapidity or pseudorapidity
        particleList[j].y[ieta] = eta;  // store particle pseudo-rapidity

        double* rapidity = new double [iptmax];
        double* cosh_y = new double [iptmax];
        double* sinh_y = new double [iptmax];
        for (int ipt = 0; ipt < iptmax; ipt++) {
            double pt = pt_array[ipt];
            double y_local;
            // rapidity as a function of pseudorapidity:
            if (DATA->pseudofreeze == 1) {
                y_local = Rap(eta, pt, m);
            } else {
                y_local = eta;
            }
            rapidity[ipt] = y_local;
            cosh_y[ipt] = cosh(y_local);
            sinh_y[ipt] = sinh(y_local);
        }

        double** temp_sum = new double* [iptmax];
        for (int ii = 0; ii < iptmax; ii++) {
            temp_sum[ii] = new double [iphimax];
            for (int jj = 0; jj < iphimax; jj++) {
                temp_sum[ii][jj] = 0.0;
            }
        }

        #pragma omp parallel
        {
            // need to perform reduction on a 2d array
            // first create a private copy for every omp_thread
            // use critical clause at the end to add the private 2d arrays
            // together from individual omp_thread
            double temp_sum_private[iptmax][iphimax];
            for (int ii = 0; ii < iptmax; ii++) {
                for (int jj = 0; jj < iphimax; jj++) {
                    temp_sum_private[ii][jj] = 0.0;
                }
            }

            #pragma omp for
            for (int icell = 0; icell < NCells; icell++) {
                double tau        = surface[icell].x[0];
                double eta_s      = surface[icell].x[3];
                double cosh_eta_s = surface[icell].cosh_eta_s;
                double sinh_eta_s = surface[icell].sinh_eta_s;

                double T   = surface[icell].T_f*hbarc;  // GeV
                double muB = surface[icell].mu_B*hbarc;  // GeV
                double mu  = baryon*muB;  // GeV
                if (DATA->whichEOS>=3 && DATA->whichEOS < 10) {
                    // for PCE use the previously computed mu
                    // at the freeze-out energy density
                    mu += mu_PCE;  // GeV
                }

                double sigma_mu[4];
                double u_flow[4];
                for (int ii = 0; ii < 4; ii++) {
                    sigma_mu[ii] = surface[icell].s[ii];
                    u_flow[ii] = surface[icell].u[ii];
                }

                double W00 = 0.0;
                double W01 = 0.0;
                double W02 = 0.0;
                double W03 = 0.0;
                double W11 = 0.0;
                double W12 = 0.0;
                double W13 = 0.0;
                double W22 = 0.0;
                double W23 = 0.0;
                double W33 = 0.0;
                int flag_shear_deltaf = 0;
                if (DATA->turn_on_shear == 1 && DATA->include_deltaf == 1) {
                    flag_shear_deltaf = 1;
                    W00 = surface[icell].W[0][0];
                    W01 = surface[icell].W[0][1];
                    W02 = surface[icell].W[0][2];
                    W03 = surface[icell].W[0][3];
                    W11 = surface[icell].W[1][1];
                    W12 = surface[icell].W[1][2];
                    W13 = surface[icell].W[1][3];
                    W22 = surface[icell].W[2][2];
                    W23 = surface[icell].W[2][3];
                    W33 = surface[icell].W[3][3];
                }

                double Pi_bulk = 0.0;
                int flag_bulk_deltaf = 0;
                if (DATA->turn_on_bulk == 1 && DATA->include_deltaf_bulk == 1) {
                    flag_bulk_deltaf = 1;
                    Pi_bulk = surface[icell].pi_b;
                    getbulkvisCoefficients(T, bulk_deltaf_coeffs);
                }
                
                double qmu_0 = 0.0;
                double qmu_1 = 0.0;
                double qmu_2 = 0.0;
                double qmu_3 = 0.0;
                double deltaf_qmu_coeff = 1.0;
                double deltaf_qmu_coeff_14mom_DV = 0.0;
                double deltaf_qmu_coeff_14mom_BV = 0.0;
                int flag_qmu_deltaf = 0;
                if (DATA->turn_on_diff == 1 && DATA->include_deltaf_qmu == 1) {
                    flag_qmu_deltaf = 1;
                    qmu_0 = surface[icell].q[0];
                    qmu_1 = surface[icell].q[1];
                    qmu_2 = surface[icell].q[2];
                    qmu_3 = surface[icell].q[3];
                    if (DATA->deltaf_14moments == 0) {
                        deltaf_qmu_coeff = get_deltaf_qmu_coeff(T, muB);
                    } else {
                        deltaf_qmu_coeff_14mom_DV =
                                        get_deltaf_coeff_14moments(T, muB, 3);
                        deltaf_qmu_coeff_14mom_BV =
                                        get_deltaf_coeff_14moments(T, muB, 4);
                    }
                }
                double rhoB = 0.0;
                if (DATA->turn_on_rhob == 1) {
                    rhoB = surface[icell].rho_B;
                }

                double eps_plus_P_over_T = surface[icell].eps_plus_p_over_T_FO;
                double prefactor_shear = 1./(2.*eps_plus_P_over_T*T*T*T)*hbarc; 
                                                                // fm^4/GeV^2
                double prefactor_qmu = rhoB/(eps_plus_P_over_T*T);   // 1/GeV

                for (int ipt = 0; ipt < iptmax; ipt++) {
                    double y = rapidity[ipt];
                    if (fabs(y - eta_s) < y_minus_eta_cut) {
                        double pt = pt_array[ipt];
                        double cosh_y_local = cosh_y[ipt];
                        double sinh_y_local = sinh_y[ipt];
                        double mt = sqrt(m*m + pt*pt);     // all in GeV
                        double ptau = mt*(cosh_y_local*cosh_eta_s
                                          - sinh_y_local*sinh_eta_s); 
                        double peta = mt*(sinh_y_local*cosh_eta_s
                                          - cosh_y_local*sinh_eta_s); 
                        for (int iphi = 0; iphi < iphimax; iphi++) {
                            double px = pt*cos_phi[iphi];
                            double py = pt*sin_phi[iphi];
                  
                            double sum;
                            // compute p^mu*dSigma_mu [fm^3*GeV]
                            double pdSigma = tau*(ptau*sigma_mu[0]
                                                  + px*sigma_mu[1]
                                                  + py*sigma_mu[2]
                                                  + peta/tau*sigma_mu[3]);
                            double E = (ptau*u_flow[0] - px*u_flow[1]
                                        - py*u_flow[2] - peta*u_flow[3]);
                            // this is the equilibrium f, f_0:
                            double f = 1./(exp(1./T*(E - mu)) + sign);
                     
                            // now comes the delta_f: check if still correct
                            // at finite mu_b 
                            // we assume here the same C=eta/s for
                            // all particle species because
                            // it is the simplest way to do it.
                            // also we assume Xi(p)=p^2, the quadratic Ansatz
                            double Wfactor = 0.0;
                            double delta_f_shear = 0.0;
                            if (flag_shear_deltaf == 1) {
                                Wfactor = (ptau*W00*ptau - 2.*ptau*W01*px
                                           - 2.*ptau*W02*py
                                           - 2.*ptau*W03*peta
                                           + px*W11*px + 2.*px*W12*py
                                           + 2.*px*W13*peta
                                           + py*W22*py + 2.*py*W23*peta
                                           + peta*W33*peta);
                                delta_f_shear = (f*(1. - sign*f)
                                                 *prefactor_shear*Wfactor);
                                if (DATA->include_deltaf==2) {
                                    // if delta f is proportional
                                    // to p^(2-alpha):
                                    delta_f_shear = (delta_f_shear
                                                     *pow((T/E), 1.*alpha)
                                                     *120./(tgamma(6.-alpha))); 
                                }
                            }
                            double delta_f_bulk = 0.0;
                            if (flag_bulk_deltaf == 1) {
                                if (bulk_deltaf_kind == 0) {
                                    delta_f_bulk = (
                                        - f*(1. - sign*f)*Pi_bulk
                                        *(bulk_deltaf_coeffs[0]*m*m
                                          + bulk_deltaf_coeffs[1]*E
                                          + bulk_deltaf_coeffs[2]*E*E));
                                } else if (bulk_deltaf_kind == 1) {
                                    double E_over_T = E/T;
                                    double mass_over_T = m/T;
                                    delta_f_bulk = (
                                        - f*(1. - sign*f)/E_over_T
                                        *bulk_deltaf_coeffs[0]
                                        *(mass_over_T*mass_over_T/3.
                                          - bulk_deltaf_coeffs[1]
                                            *E_over_T*E_over_T)
                                        *Pi_bulk);
                                } else if (bulk_deltaf_kind == 2) {
                                    double E_over_T = E/T;
                                    delta_f_bulk = (
                                        - f*(1. - sign*f)
                                        *(-bulk_deltaf_coeffs[0]
                                          + bulk_deltaf_coeffs[1]*E_over_T)
                                        *Pi_bulk);
                                } else if (bulk_deltaf_kind == 3) {
                                    double E_over_T = E/T;
                                    delta_f_bulk = (
                                        - f*(1.-sign*f)/sqrt(E_over_T)
                                            *(- bulk_deltaf_coeffs[0]
                                              + bulk_deltaf_coeffs[1]*E_over_T)
                                            *Pi_bulk);
                                } else if (bulk_deltaf_kind == 4) {
                                    double E_over_T = E/T;
                                    delta_f_bulk = (
                                        - f*(1.-sign*f)
                                            *(bulk_deltaf_coeffs[0]
                                              - bulk_deltaf_coeffs[1]/E_over_T)
                                            *Pi_bulk);
                                }
                            }

                            // delta f for qmu
                            double qmufactor = 0.0;
                            double delta_f_qmu = 0.0;
                            if (flag_qmu_deltaf == 1) {
                                // p^\mu q_\mu
                                qmufactor = (ptau*qmu_0 - px*qmu_1 - py*qmu_2
                                             - peta*qmu_3);
                                if (DATA->deltaf_14moments == 0) {
                                    delta_f_qmu = (
                                        f*(1. - sign*f)
                                         *(prefactor_qmu - baryon/E)*qmufactor
                                         /deltaf_qmu_coeff);
                                } else {
                                    delta_f_qmu = (
                                        f*(1. - sign*f)
                                         *(baryon*deltaf_qmu_coeff_14mom_DV
                                           + 2.*deltaf_qmu_coeff_14mom_BV*E)
                                         *qmufactor);
                                }
                            }

                            double max_ratio = 1.0;
                            double total_deltaf = (delta_f_shear
                                                   + delta_f_bulk
                                                   + delta_f_qmu);

                            if (fabs(total_deltaf)/f > max_ratio) {
                                total_deltaf *= f/fabs(total_deltaf);
                            }
                            sum = (f + total_deltaf)*pdSigma;
                     
                            if (sum > 10000) {
                                music_message << "sum>10000 in summation. sum = "
                                     << sum 
                                     << ", f=" << f << ", deltaf="
                                     << delta_f_shear
                                     << ", pdSigma=" << pdSigma << ", T=" << T 
                                     << ", E=" << E << ", mu=" << mu;
                                music_message.flush("warning");
                            }
                            if (f < 0.) {
                                music_message << " f_eq < 0.! f_eq = " << f
                                              << ", T = " << T << " GeV, mu = "
                                              << mu << " GeV, E = "
                                              << E << " GeV";
                                music_message.flush("error");
                            }
                            temp_sum_private[ipt][iphi] += sum;
                        }
                    }
                }
            }
            #pragma omp critical
            {
                for (int ipt = 0; ipt < iptmax; ipt++) {
                    for (int iphi = 0; iphi < iphimax; iphi++) {
                        temp_sum[ipt][iphi] += temp_sum_private[ipt][iphi];
                    }
                }
            }
        }
        double prefactor = deg/(pow(2.*M_PI,3.)*pow(hbarc,3.));
        // store the final results
        for (int ipt = 0; ipt < iptmax; ipt++) {
            for (int iphi = 0; iphi < iphimax; iphi++) {
                double sum = temp_sum[ipt][iphi]*prefactor;   // in GeV^(-2)
                particleList[j].dNdydptdphi[ieta][ipt][iphi] = sum;
                fprintf(s_file,"%e ", sum);
            }
            fprintf(s_file,"\n");
        }
        // clean up
        delete[] rapidity;
        delete[] cosh_y;
        delete[] sinh_y;
        for (int ipt = 0; ipt < iptmax; ipt++)
            delete[] temp_sum[ipt];
        delete[] temp_sum;
    }
    
    if (DATA_ptr->turn_on_bulk == 1 && DATA_ptr->include_deltaf_bulk == 1) {
        delete[] bulk_deltaf_coeffs;
    }

    // clean up
    delete[] cos_phi;
    delete[] sin_phi;
    delete[] pt_array;
    fclose(s_file);
    fclose(d_file);
}


void Freeze::ComputeParticleSpectrum_pseudo_boost_invariant(
                                        InitData *DATA, int number) {
    // this function compute thermal paritcle spectra assuming a
    // boost-invarianat hyper-surface from hydro simulations
    int j = partid[MHALF+number];
    // set some parameters
    double etamax = DATA->max_pseudorapidity;
    int ietamax = DATA->pseudo_steps + 1;  // pseudo_steps is number of steps.
                                           // Including edges
                                           // number of points is steps + 1
    double deltaeta = 1.0;
    if (ietamax > 1)
        deltaeta = 2.*etamax/DATA->pseudo_steps;
    double ptmax = DATA->max_pt;
    double ptmin = DATA->min_pt;
    int iptmax = DATA->pt_steps+1;  // Number of points is iptmax + 1
                                    // (ipt goes from 0 to iptmax)
    int iphimax = DATA->phi_steps;  // number of points
                                    // (phi=2pi equal to phi=0)
    double deltaphi = 2*M_PI/iphimax;
  
    // Reuse rapidity variables (Need to reuse variable y
    // so resonance decay routine can be used as is.
    // Might as well misuse ymax and deltaY too)
    particleList[j].ymax = etamax; 
    particleList[j].deltaY = deltaeta;

    music_message << "Doing " << j << ": "
                  << particleList[j].name << "("
                  << particleList[j].number << ") ...";
    music_message.flush("info");
 
    particleList[j].ny = DATA->pseudo_steps + 1;
    particleList[j].npt = iptmax;
    particleList[j].nphi = iphimax;
  
    // set particle properties
    double m = particleList[j].mass;
    int deg = particleList[j].degeneracy;
    int baryon = particleList[j].baryon;
    // int s = particleList[j].strange;
    // int c = particleList[j].charge;
    double mu_PCE = particleList[j].muAtFreezeOut;
    int sign;
    if (baryon == 0) {
        sign = -1.;
    } else {
        sign = 1.;
    }

    // open files to write
    FILE *d_file;
    string d_name = "particleInformation.dat";
    d_file = fopen(d_name.c_str(), "a");

    string s_name = "yptphiSpectra0.dat";
    FILE *s_file;
    s_file = fopen(s_name.c_str(), "w");

    // -----------------------------------------------------------------------
    // write information in the particle information file
    fprintf(d_file, "%d %e %d %e %e %d %d \n",
            number,  etamax, DATA->pseudo_steps+1, ptmin, ptmax,
            iptmax, iphimax);
    // caching 
    double *cos_phi = new double[iphimax];
    double *sin_phi = new double[iphimax];
    for (int iphi = 0; iphi < iphimax; iphi++) {
        double phi_local = deltaphi*iphi;
        cos_phi[iphi] = cos(phi_local);
        sin_phi[iphi] = sin(phi_local);
    }
    double* pt_array = new double[iptmax];
    for (int ipt = 0; ipt < iptmax; ipt++) {
        double pt =  (ptmin + (ptmax - ptmin)
                              *pow(static_cast<double>(ipt), 2.)
                              /pow(static_cast<double>(iptmax - 1), 2.));
        pt_array[ipt] = pt;
        particleList[j].pt[ipt] = pt;
    }
    
    double *bulk_deltaf_coeffs = new double[3];
    if (DATA_ptr->turn_on_bulk == 1 && DATA_ptr->include_deltaf_bulk == 1) {
        if (bulk_deltaf_kind == 0) {
            bulk_deltaf_coeffs = new double[3];
        } else {
            bulk_deltaf_coeffs = new double[2];
        }
    }
      
    double alpha = 0.0;
    double** temp_sum = new double* [iptmax];
    for (int ii = 0; ii < iptmax; ii++) {
        temp_sum[ii] = new double [iphimax];
        for (int jj = 0; jj < iphimax; jj++) {
            temp_sum[ii][jj] = 0.0;
        }
    }
    
    // main loop begins ...
    // store E dN/d^3p as function of phi,
    // pt and eta (pseudorapidity) in sumPtPhi:

    #pragma omp parallel
    {
        // need to perform reduction on a 2d array
        // first create a private copy for every omp_thread
        // use critical clause at the end to add the private 2d arrays
        // together from individual omp_thread
        double temp_sum_private[iptmax][iphimax];
        for (int ii = 0; ii < iptmax;  ii++) 
        for (int jj = 0; jj < iphimax; jj++) {
            temp_sum_private[ii][jj] = 0.0;
        }

        #pragma omp for
        for (int icell = 0; icell < NCells; icell++) {
            double tau = surface[icell].x[0];

            double T = surface[icell].T_f*hbarc;  // GeV
            double muB = 0.0;
            double mu = baryon*muB;  // GeV
            if (DATA->whichEOS>=3 && DATA->whichEOS < 10) {
                // for PCE use the previously computed mu
                // at the freeze-out energy density
                mu += mu_PCE;  // GeV
            }

            double sigma_mu[4];
            double u_flow[4];
            for (int ii = 0; ii < 4; ii++) {
                sigma_mu[ii] = surface[icell].s[ii];
                u_flow[ii] = surface[icell].u[ii];
            }

            double W00 = 0.0;
            double W01 = 0.0;
            double W02 = 0.0;
            double W03 = 0.0;
            double W11 = 0.0;
            double W12 = 0.0;
            double W13 = 0.0;
            double W22 = 0.0;
            double W23 = 0.0;
            double W33 = 0.0;
            int flag_shear_deltaf = 0;
            if (DATA->turn_on_shear == 1 && DATA->include_deltaf == 1) {
                flag_shear_deltaf = 1;
                W00 = surface[icell].W[0][0];
                W01 = surface[icell].W[0][1];
                W02 = surface[icell].W[0][2];
                W03 = surface[icell].W[0][3];
                W11 = surface[icell].W[1][1];
                W12 = surface[icell].W[1][2];
                W13 = surface[icell].W[1][3];
                W22 = surface[icell].W[2][2];
                W23 = surface[icell].W[2][3];
                W33 = surface[icell].W[3][3];
            }

            double Pi_bulk = 0.0;
            int flag_bulk_deltaf = 0;
            if (DATA->turn_on_bulk == 1 && DATA->include_deltaf_bulk == 1) {
                flag_bulk_deltaf = 1;
                Pi_bulk = surface[icell].pi_b;
                getbulkvisCoefficients(T, bulk_deltaf_coeffs);
            }

            double eps_plus_P_over_T = surface[icell].eps_plus_p_over_T_FO;
            double prefactor_shear = 1./(2.*eps_plus_P_over_T*T*T*T)*hbarc; 
                                                                // fm^4/GeV^2

            for (int ieta_s = 0; ieta_s < n_eta_s_integral; ieta_s++) {
                double cosh_eta_s = cosh_eta_s_inte[ieta_s];
                double sinh_eta_s = sinh_eta_s_inte[ieta_s];
                for (int ipt=0; ipt<iptmax; ipt++) {
                    double pt = pt_array[ipt];
                    double mt = sqrt(m*m + pt*pt);     // all in GeV
                    double ptau = mt*cosh_eta_s;
                    // sinh(y - eta_s) = - sinh(eta_s)
                    double peta = - mt*sinh_eta_s;
                    for (int iphi = 0; iphi < iphimax; iphi++) {
                        double px = pt*cos_phi[iphi];
                        double py = pt*sin_phi[iphi];
                  
                        double sum;
                        // compute p^mu*dSigma_mu [fm^3*GeV]
                        double pdSigma = tau*(ptau*sigma_mu[0]
                                              + px*sigma_mu[1]
                                              + py*sigma_mu[2]
                                              + peta*sigma_mu[3]/tau);
                        double E = (ptau*u_flow[0] - px*u_flow[1]
                                    - py*u_flow[2]- peta*u_flow[3]);
                        // this is the equilibrium f, f_0:
                        double f = 1./(exp(1./T*(E - mu)) + sign);
                    
                        // now comes the delta_f: check if still correct
                        // at finite mu_b 
                        // we assume here the same C=eta/s for
                        // all particle species because
                        // it is the simplest way to do it.
                        // also we assume Xi(p)=p^2, the quadratic Ansatz
                        double Wfactor = 0.0;
                        double delta_f_shear = 0.0;
                        if (flag_shear_deltaf == 1) {
                            Wfactor = (ptau*W00*ptau - 2.*ptau*W01*px
                                       - 2.*ptau*W02*py
                                       - 2.*ptau*W03*peta
                                       + px*W11*px + 2.*px*W12*py
                                       + 2.*px*W13*peta
                                       + py*W22*py + 2.*py*W23*peta
                                       + peta*W33*peta);
                            delta_f_shear = (f*(1. - sign*f)
                                             *prefactor_shear*Wfactor);
                            if (DATA->include_deltaf==2) {
                                // if delta f is supposed to be proportional
                                // to p^(2-alpha):
                                delta_f_shear = (delta_f_shear
                                                 *pow((T/E), 1.*alpha)
                                                 *120./(tgamma(6.-alpha))); 
                            }
                        }
                        double delta_f_bulk = 0.0;
                        if (flag_bulk_deltaf == 1) {
                            if (bulk_deltaf_kind == 0) {
                                delta_f_bulk = (
                                    - f*(1. - sign*f)*Pi_bulk
                                    *(bulk_deltaf_coeffs[0]*m*m
                                      + bulk_deltaf_coeffs[1]*E
                                      + bulk_deltaf_coeffs[2]*E*E));
                            } else if (bulk_deltaf_kind == 1) {
                                double E_over_T = E/T;
                                double mass_over_T = m/T;
                                delta_f_bulk = (
                                    - f*(1. - sign*f)/E_over_T
                                    *bulk_deltaf_coeffs[0]
                                    *(mass_over_T*mass_over_T/3.
                                      - bulk_deltaf_coeffs[1]
                                        *E_over_T*E_over_T)
                                    *Pi_bulk);
                            } else if (bulk_deltaf_kind == 2) {
                                double E_over_T = E/T;
                                delta_f_bulk = (
                                    - f*(1. - sign*f)
                                    *(-bulk_deltaf_coeffs[0]
                                      + bulk_deltaf_coeffs[1]*E_over_T)
                                    *Pi_bulk);
                            } else if (bulk_deltaf_kind == 3) {
                                double E_over_T = E/T;
                                delta_f_bulk = (
                                    - f*(1.-sign*f)/sqrt(E_over_T)
                                        *(- bulk_deltaf_coeffs[0]
                                          + bulk_deltaf_coeffs[1]*E_over_T)
                                        *Pi_bulk);
                            } else if (bulk_deltaf_kind == 4) {
                                double E_over_T = E/T;
                                delta_f_bulk = (
                                    - f*(1.-sign*f)
                                        *(bulk_deltaf_coeffs[0]
                                          - bulk_deltaf_coeffs[1]/E_over_T)
                                        *Pi_bulk);
                            }
                        }
                        double total_deltaf = delta_f_shear + delta_f_bulk;
                        double max_ratio = 1.0;
                        if (fabs(total_deltaf)/f > max_ratio) {
                            total_deltaf *= f/fabs(total_deltaf);
                        }
                        sum = (f + total_deltaf)*pdSigma;
                    
                        if (sum > 10000) {
                            music_message << "sum>10000 in summation. sum = "
                                 << sum 
                                 << ", f=" << f << ", deltaf=" << delta_f_shear
                                 << ", pdSigma=" << pdSigma << ", T=" << T 
                                 << ", E=" << E << ", mu=" << mu;
                            music_message.flush("warning");
                        }
                        temp_sum_private[ipt][iphi] += (
                                                sum*eta_s_inte_weight[ieta_s]);
                    }
                }
            }
        }
        #pragma omp critical
        {
            for (int ipt = 0; ipt < iptmax; ipt++) {
                for (int iphi = 0; iphi < iphimax; iphi++) {
                    temp_sum[ipt][iphi] += temp_sum_private[ipt][iphi];
                }
            }
        }
    }
    double prefactor = deg/(pow(2.*M_PI, 3.)*pow(hbarc, 3.));

    for (int ieta = 0; ieta < ietamax; ieta++) {
        double eta = -etamax + ieta*deltaeta;
        // Use this variable to store pseudorapidity instead of rapidity
        // May cause confusion in the future,
        // but easier to to share code for both options:
        // calculating on a fixed grid in rapidity or pseudorapidity
        particleList[j].y[ieta] = eta;  // store particle pseudo-rapidity

        // store the final results
        for (int ipt = 0; ipt < iptmax; ipt++) {
            for(int iphi = 0; iphi < iphimax; iphi++) {
                double sum = temp_sum[ipt][iphi]*prefactor;
                particleList[j].dNdydptdphi[ieta][ipt][iphi] = sum;
                fprintf(s_file, "%e ", sum);
            }
            fprintf(s_file, "\n");
        }
    }
    // clean up
    for (int ipt = 0; ipt < iptmax; ipt++) {
        delete[] temp_sum[ipt];
    }
    delete[] temp_sum;
    
    if (DATA_ptr->turn_on_bulk == 1 && DATA_ptr->include_deltaf_bulk == 1) {
        delete [] bulk_deltaf_coeffs;
    }

    // clean up
    delete[] cos_phi;
    delete[] sin_phi;
    delete[] pt_array;
    fclose(s_file);
    fclose(d_file);
}

void Freeze::OutputFullParticleSpectrum_pseudo(InitData *DATA, int number,
                                               int anti, int full) {
    FILE *d_file;
    const char* d_name = "FparticleInformation.dat";
    d_file = fopen(d_name, "a");
  
    fprintf(d_file, "%d %e %d %e %e %d %d \n",
            number,  DATA->max_pseudorapidity, DATA->pseudo_steps+1,
            DATA->min_pt, DATA->max_pt, DATA->pt_steps+1, DATA->phi_steps);

    fclose(d_file);
  
    FILE *s_file;
    const char* s_name = "FyptphiSpectra.dat";
    s_file = fopen(s_name, "a");

    int j = partid[MHALF+number];
    for (int iy = 0; iy <= DATA->pseudo_steps; iy++) {
        for (int ipt = 0; ipt <= DATA->pt_steps; ipt++) {
            for (int iphi = 0; iphi < DATA->phi_steps; iphi++) {
                fprintf(s_file, "%e ",
                        particleList[j].dNdydptdphi[iy][ipt][iphi]);
            }
            fprintf(s_file,"\n");
        }
    }
    fclose(s_file);
}


//! this function computes particle thermal spectra
void Freeze::compute_thermal_spectra(int particleSpectrumNumber,
                                     InitData* DATA) {
    double mass_tol = 1e-3;
    double mu_tol = 1e-3;

    // clean up
    system_status_ = system("rm yptphiSpectra.dat yptphiSpectra?.dat "
           "yptphiSpectra??.dat particleInformation.dat 2> /dev/null");

    ReadFreezeOutSurface(DATA);  // read freeze out surface
    if (particleSpectrumNumber == 0) {
        // do all particles up to particleMax
        music_message.info("Doing all particles. May take a while ...");

        int particleMax_copy = particleMax;
        for (int i = 1; i < particleMax_copy; i++) {
            int number = particleList[i].number;
            int computespectrum = 1;

            // Only calculate particles with unique mass
            for(int part = 1; part < i; part++) {
                double mass_diff = fabs(particleList[i].mass
                                        - particleList[part].mass);
                double mu_diff = fabs(particleList[i].muAtFreezeOut
                                      - particleList[part].muAtFreezeOut);
                if (mass_diff < mass_tol && mu_diff < mu_tol
                    && (DATA->turn_on_rhob == 0
                        || particleList[i].baryon == particleList[part].baryon)
                   ) {
                    // here we assume zero mu_B
                    computespectrum = 0;

                    // If there is more than one processor,
                    // this processor doesn't have all pseudorapidity
                    // values in memory
                    music_message << "Copying " << i << ":"
                                  << particleList[i].name << " ("
                                  << particleList[i].number << ") from "
                                  << particleList[part].name;
                    music_message.flush("info");

                    int iphimax = DATA->phi_steps;
                    int iptmax = DATA->pt_steps + 1;
                    int ietamax = DATA->pseudo_steps + 1;
                    double ptmax = DATA->max_pt;
                    double ptmin = DATA->min_pt;
                    double etamax = DATA->max_pseudorapidity;
                    // If the particles have a different degeneracy,
                    // we have to multiply by the ratio when copying.
                    double degen_ratio = (
                        static_cast<double>(particleList[i].degeneracy)
                        /static_cast<double>(particleList[part].degeneracy)
                    );

                    // open files to write
                    FILE *d_file;
                    const char* d_name = "particleInformation.dat";
                    d_file = fopen(d_name, "a");
                    FILE *s_file;
                    const char* s_name = "yptphiSpectra.dat";
                    s_file = fopen(s_name, "a");
                    particleList[i].ymax = particleList[part].ymax; 
                    particleList[i].deltaY = particleList[part].deltaY;
                    particleList[i].ny = particleList[part].ny;
                    particleList[i].npt = particleList[part].npt;
                    particleList[i].nphi = particleList[part].nphi;
                    fprintf(d_file, "%d %e %d %e %e %d %d \n",
                            number, etamax, ietamax, ptmin, ptmax,
                            iptmax, iphimax);
                    for (int ieta = 0; ieta < ietamax; ieta++) {
                        for (int ipt = 0; ipt < iptmax; ipt++) {
                            particleList[i].pt[ipt] =
                                            particleList[part].pt[ipt];
                            particleList[i].y[ieta] =
                                            particleList[part].y[ieta];  
                            for (int iphi = 0; iphi < iphimax; iphi++) {
                                particleList[i].dNdydptdphi[ieta][ipt][iphi] =
                                    (degen_ratio
                                     *particleList[part].dNdydptdphi[ieta][ipt][iphi]);
                                fprintf(s_file, "%e ",
                                        particleList[i].dNdydptdphi[ieta][ipt][iphi]);
                            }
                            fprintf(s_file, "\n");
                        }
                    }
                    fclose(s_file);
                    fclose(d_file);
                    part = particleMax_copy; // break out of particle loop
                }  // if particles have same mass
            }  // loop over particles that have already been calculated

            if (computespectrum) {
                if (boost_invariant) {
                    ComputeParticleSpectrum_pseudo_boost_invariant(
                                                                DATA, number);
                } else {
                    ComputeParticleSpectrum_pseudo_improved(DATA, number);
                }

                system_status_ = system(
                        "cat yptphiSpectra?.dat >> yptphiSpectra.dat");
                system_status_ = system(
                        "cat yptphiSpectra??.dat >> yptphiSpectra.dat "
                        "2> /dev/null");
            }
        }
    } else {
        // compute single one particle with pid = particleSpectrumNumber
        if (particleSpectrumNumber >= particleMax) {
            music_message << "No particle has the number "
                          << particleSpectrumNumber << " Exiting.";
            music_message.flush("error");
            exit(1);
        }
        int number = particleList[particleSpectrumNumber].number;
        music_message.info("COMPUTE");
        if (boost_invariant) {
            ComputeParticleSpectrum_pseudo_boost_invariant(DATA, number);
        } else {
            ComputeParticleSpectrum_pseudo_improved(DATA, number);
        }
        system_status_ = system("cat yptphiSpectra?.dat >> yptphiSpectra.dat");
        system_status_ = system("cat yptphiSpectra??.dat >> yptphiSpectra.dat "
                                "2> /dev/null");
    }
}

void Freeze::perform_resonance_decays(InitData *DATA) {
    ReadSpectra_pseudo(DATA, 0, 1);
    int bound = 211; //number of lightest particle to calculate. 
    music_message << "particleMax = " << particleMax;
    music_message.flush("info");
    music_message << "doing all from " << particleMax << ": "
                  << particleList[particleMax].name << "to "
                  << partid[MHALF+bound] << ": "
                  << particleList[partid[MHALF+bound]].name;
    music_message.flush("info");

    cal_reso_decays(particleMax, decayMax, bound);

    system_status_ = system(
            "rm FyptphiSpectra.dat FparticleInformation.dat 2> /dev/null");
    for (int i = 1; i < particleMax; i++) {
        int number = particleList[i].number;
        int b = particleList[i].baryon;
        OutputFullParticleSpectrum_pseudo(DATA, number, b, 1);
    }
}

void Freeze::compute_thermal_particle_spectra_and_vn(InitData* DATA) {
    // take tabulated spectra
    // and compute various observables and integrated quantities
    int status = mkdir("./outputs", 0755);
    if (status != 0) {
        music_message << "Can not create folder ./outputs ...";
        music_message.flush("warning");
    }
    ReadSpectra_pseudo(DATA, 0, 1);

    //Make sure the pion, kaon and proton spectra are available
    if (particleMax < 21) {
        music_message << __func__ 
                      << ": cannot compute pion, kaon and proton observables, "
                      << "some of the spectra are not available.";
        music_message.flush("warning");
    } else {
        int pid_list [] = {211, -211, 321, -321, 2212, -2212};
        int pid_list_length = sizeof(pid_list)/sizeof(int);
        for (int k = 0; k < pid_list_length; k++) {
            int number = pid_list[k];
            int id = partid[MHALF+number];
            if (id < particleMax) {
                OutputDifferentialFlowAtMidrapidity(DATA, number, 0);
                OutputDifferentialFlowNearMidrapidity(DATA, number, 0);
                OutputIntegratedFlow_vs_y(DATA, number, 0, 0.01, 3.0);
            }
        }
    }
}

void Freeze::compute_final_particle_spectra_and_vn(InitData* DATA) {
    // take tabulated post-decay spectra and
    // compute various observables and integrated quantities
    int status = mkdir("./outputs", 0755);   // output results folder
    if (status != 0) {
        music_message << "Can not create folder ./outputs ...";
        music_message.flush("warning");
    }
    ReadSpectra_pseudo(DATA, 1, 1);  // read in particle spectra information
    
    if (DATA->NumberOfParticlesToInclude > 200) {
        if ((DATA-> whichEOS != 7 && particleMax < 320)
                || ((DATA->whichEOS == 7) && particleMax < 327)) {
            music_message << "Not all hadronic species are computed "
                          << "particleMax = " << particleMax;
            music_message.flush("warning");
        }
        // calculate spectra and vn for charged hadrons
        Output_charged_hadrons_eta_differential_spectra(DATA, 1, 0.01, 3.0);
        Output_charged_hadrons_pT_differential_spectra(DATA, 1, -0.5, 0.5);
        Output_charged_IntegratedFlow(DATA, 0.01, 3.0, -0.5, 0.5);
        // for PHENIX pT cut
        Output_charged_IntegratedFlow(DATA, 0.15, 3.0, -0.5, 0.5);
        // for ALICE or STAR pT cut
        Output_charged_IntegratedFlow(DATA, 0.2, 3.0, -0.5, 0.5);
        // for CMS pT cut
        Output_charged_IntegratedFlow(DATA, 0.3, 3.0, -0.5, 0.5);

        Output_charged_IntegratedFlow(DATA, 0.4, 3.0, -0.5, 0.5);
        // for ATLAS pT cut
        Output_charged_IntegratedFlow(DATA, 0.5, 3.0, -0.5, 0.5);
    }

    if (particleMax < 21) {
        music_message << __func__ 
                      << "cannot compute pion, kaon and proton observables, "
                      << "some of the spectra are not available.";
        music_message.flush("warning");
    } else {
        // calculate spectra and vn for identified particles
        int pid_list [] = {211, -211, 321, -321, 2212, -2212};
        int pid_list_length = sizeof(pid_list)/sizeof(int);
        for (int k = 0; k < pid_list_length; k++) {
            int number = pid_list[k];
            int id = partid[MHALF+number];
            if (id < particleMax) {
                OutputDifferentialFlowAtMidrapidity(DATA, number, 1);
                OutputDifferentialFlowNearMidrapidity(DATA, number, 1);
                OutputIntegratedFlow_vs_y(DATA, number, 1, 0.01, 3.0);
            }
        }
    }
}


void Freeze::CooperFrye_pseudo(int particleSpectrumNumber, int mode,
                               InitData *DATA, EOS *eos) {
    // this is a shell function for Cooper-Frye routine
    // -- spectra calculated on an equally-spaced grid
    // in phi and (pseudo)rapidity
    ReadParticleData(DATA, eos); // read in data for Cooper-Frye
    if (mode == 3 || mode == 1) {  // compute thermal spectra
        compute_thermal_spectra(particleSpectrumNumber, DATA);
    }
    if (mode==4 || mode==1) { //  do resonance decays
        perform_resonance_decays(DATA);
    } 
    if (mode==13 || mode == 3 || mode == 1) {
        compute_thermal_particle_spectra_and_vn(DATA);
    } 
    if (mode==14 || mode == 4 || mode == 1) {
        compute_final_particle_spectra_and_vn(DATA);
    } 

    // clean up
    delete[] partid;
    free(particleList);
}


// returns pseudorapidity for a given rapidity, transverse momentum, and mass
double Freeze::PseudoRap(double y, double pt, double m) {
    double eta = acosh(2*m*m/pt/pt*sinh(y)*sinh(y) + cosh(2*y))/2.;
    if (y < 0)
        eta *= -1.;
    return eta;
}

// returns rapidity for a given pseudorapidity, transverse momentum, and mass
double Freeze::Rap(double eta, double pt, double m) {
    double y = log((sqrt(m*m + pt*pt*cosh(eta)*cosh(eta)) + pt*sinh(eta))
                   /sqrt(m*m+pt*pt));
    return y;
}

// jacobian to convert dN/dy to dN/deta:  dN/deta = dydeta*dN/dy
double Freeze::dydeta(double eta, double pt, double m) {
    return pt*cosh(eta)/sqrt(m*m + pt*pt*cosh(eta)*cosh(eta));
}


// calculates eta= or y-integrated flow versus rapidity.
// format is vn[n][i(real=0 or imaginary part=1)][pt] .
// yflag = 0 takes minrap and maxrap as a range in psuedorapidity,
// while yflag = 1 uses rapidity
// Yield (n=0) is dN/dpt/dy or dN/dpt/deta for minrap==maxrap,
// or dN/dpt otherwise 
void Freeze::rapidity_integrated_flow(
        InitData *DATA, int number, int yflag, double minrap, double maxrap,
        double vn[nharmonics][2][etasize]) {
    int j = partid[MHALF+number];
    int nphi = particleList[j].nphi;
    int npt = particleList[j].npt;
    int neta = particleList[j].ny;
    double intvn[ptsize][nharmonics][2] = {};
    double m = particleList[j].mass;
   
    double testmin;
    if (yflag) {
        if (DATA->pseudofreeze == 1)
            testmin = Rap(particleList[j].y[0], particleList[j].pt[0], m);
        else
            testmin = particleList[j].y[0];
    } else {
        if (DATA->pseudofreeze == 1)
            testmin = particleList[j].y[0];
        else
            testmin = PseudoRap(particleList[j].y[0], particleList[j].pt[0],
                                m);
    }

    if (minrap < testmin) {
        music_message << "Error: called out of range rapidity in rap_integrated_flow, "
             << minrap << " < minimum " << testmin << "\n";
        music_message << particleList[j].y[0] << "  " << particleList[j].pt[0] << "   "
             << m;
        music_message.flush("error");
        exit(1);
    }
    double testmax;
    if (yflag) {
        if (DATA->pseudofreeze == 1)
            testmax = Rap(particleList[j].y[neta-1], particleList[j].pt[0], m);
        else
            testmax = particleList[j].y[neta-1];
    } else {
        if (DATA->pseudofreeze == 1)
            testmax = particleList[j].y[neta-1];
        else
            testmax = Rap(particleList[j].y[neta-1], particleList[j].pt[0], m);
    }
    if (maxrap > testmax) {
        music_message << "Error: called out of range rapidity in rap_integrated_flow, "
             << maxrap << " > maximum " << testmax;
        music_message.flush("error");
        exit(1);
    }
    if (minrap > maxrap) {
        music_message << "Error in rap_integrated_flow: "
              << "minrap must be less than or equal to maxrap";
        music_message.flush("error");
        exit(1);
    }
   
    // loop over pt
    for (int ipt = 0; ipt < npt; ipt++) {
        double pt = particleList[j].pt[ipt];
     
        // Integrate over phi using trapezoid rule
        for (int iphi = 0; iphi < nphi; iphi++) {
            double ylist[etasize] = {0};
            double etalist[etasize] = {0};
           
            // Integrate over pseudorapidity using gsl
            double dndpt[etasize] = {0};
            for (int ieta = 0; ieta < neta; ieta++) {
                double rap_local = particleList[j].y[ieta];

                double eta_local, y_local;
                if (DATA->pseudofreeze == 1) {
                     eta_local = rap_local;
                     etalist[ieta] = eta_local;
                     y_local = Rap(eta_local, pt, m);   // get rapidity
                     ylist[ieta] = y_local;   // get rapidity
                } else {
                     y_local = rap_local;
                     ylist[ieta] = y_local;
                     // get Pseudo-rapidity
                     eta_local = PseudoRap(y_local, pt, m);
                     etalist[ieta] = eta_local;   // get Pseudo-rapidity
                }
                if (yflag) {
                     // dN/(dydpTdphi)
                     dndpt[ieta] = (
                             pt*particleList[j].dNdydptdphi[ieta][ipt][iphi]);
                } else {
                     // dN/(detadpTdphi)
                     dndpt[ieta] = (
                             pt*dydeta(eta_local, pt, m)
                             *particleList[j].dNdydptdphi[ieta][ipt][iphi]);
                }
            }
            gsl_interp_accel *acc = gsl_interp_accel_alloc ();
            gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, neta);
            if (yflag) {
                // rapidity
                gsl_spline_init (spline, ylist, dndpt , neta);
            } else {
                // pseudo-rapidity
                gsl_spline_init (spline, etalist, dndpt , neta);
            }
   
            double dNdp;
            if (minrap < maxrap) {
                // integrated from minrap to maxrap
                dNdp = gsl_spline_eval_integ(spline, minrap, maxrap, acc);
            } else {
                // value at minrap = maxrap
                dNdp = gsl_spline_eval(spline, maxrap, acc);
            }
           
            double phi = iphi*2*M_PI/nphi;
            for (int i = 0; i < nharmonics; i++) {
                intvn[ipt][i][0] += cos(i*phi)*dNdp*2*M_PI/nphi;
                intvn[ipt][i][1] += sin(i*phi)*dNdp*2*M_PI/nphi;
            }
            gsl_spline_free (spline);
            gsl_interp_accel_free (acc);
        }  // phi loop
   
        for (int k = 0; k < 2; k++) {
             vn[0][k][ipt] = intvn[ipt][0][k];
        }
        
        for (int i = 1; i < nharmonics; i++) {
             for(int k = 0; k < 2; k++) {
                 vn[i][k][ipt] = intvn[ipt][i][k]/intvn[ipt][0][0];
             }
        }
    }  // pt loop
}


// calculates pt- and (pseudo)rapidity-integrated flow for a given range
// in pt and eta or y
// format is vn[n][i(real=0 or imaginary part=1)]
// the rapidity integral is nested inside the phi integral
// inside the pt integral
void Freeze::pt_and_rapidity_integrated_flow(InitData *DATA, int number,
        double minpt, double maxpt, int yflag,  double minrap, double maxrap,
        double vn[nharmonics][2]) {
    int j = partid[MHALF+number];
    int npt = particleList[j].npt;  
    
    if (minpt < particleList[j].pt[0]) {
        music_message << "Error: called out of range pt in "
             << "pt_and_rapidity_integrated_flow, " 
             << minpt << " < minimum " << particleList[j].pt[0];
        music_message.flush("error");
        exit(1);
    }
    if (maxpt > particleList[j].pt[npt-1]) {
        music_message << "Error: called out of range pt in "
             << "pt_and_rapidity_integrated_flow, " 
             << maxpt << " > maximum " << particleList[j].pt[npt-1];
        music_message.flush("error");
        exit(1);
    }
    if (minpt > maxpt) {
        music_message << "Error in pt_and_rapidity_integrated_flow: "
             << "minpt must be less than or equal to maxpt";
        music_message.flush("error");
        exit(1);
    }
  
    // do rapidity-integral first
    double vnpt[nharmonics][2][etasize] = {};
    rapidity_integrated_flow(DATA, number, yflag,  minrap, maxrap, vnpt);
    double dndpt[ptsize];
    for (int ipt = 0; ipt < npt; ipt++) {
        dndpt[ipt] = vnpt[0][0][ipt];
    }
    for (int n = 0; n < nharmonics; n++) {
        double vncos[ptsize] = {0};
        double vnsin[ptsize] = {0};
        for(int ipt = 0; ipt < npt; ipt++) {
            double weight;
            if (n != 0)
                weight = dndpt[ipt];
            else
                weight = 1;
            vncos[ipt] = vnpt[n][0][ipt]*weight;
            vnsin[ipt] = vnpt[n][1][ipt]*weight;
        }
        gsl_interp_accel *cosacc = gsl_interp_accel_alloc ();
        gsl_spline *cosspline = gsl_spline_alloc (gsl_interp_linear, npt);
        gsl_spline_init (cosspline, particleList[j].pt ,vncos , npt);

        gsl_interp_accel *sinacc = gsl_interp_accel_alloc ();
        gsl_spline *sinspline = gsl_spline_alloc (gsl_interp_linear, npt);
        gsl_spline_init (sinspline, particleList[j].pt ,vnsin , npt);

        if (fabs(minpt - maxpt) > 1e-6) {
            vn[n][0] = gsl_spline_eval_integ(cosspline, minpt, maxpt, cosacc);
            vn[n][1] = gsl_spline_eval_integ(sinspline, minpt, maxpt, sinacc);
        } else {
            vn[n][0] = gsl_spline_eval(cosspline, minpt, cosacc);
            vn[n][1] = gsl_spline_eval(sinspline, minpt, sinacc);
        }
    
        if (n != 0) {
            for (int i = 0; i < 2; i++) {
                vn[n][i] /= vn[0][0];
            }
        }
    
        // clean up
        gsl_spline_free (cosspline);
        gsl_interp_accel_free (cosacc);
        gsl_spline_free (sinspline);
        gsl_interp_accel_free (sinacc);
    }
}


// Output yield dN/ptdydpt and v_n at eta=0 as a function of pT
void Freeze::OutputDifferentialFlowAtMidrapidity(
        InitData *DATA, int number, int full) {
    // Define index j used in particleList[j]
    int j = partid[MHALF+number];
    int npt = particleList[j].npt;

    music_message << "Calculating flow at midrapidity for particle " << number;
    music_message.flush("info");

    //Set output file name
    string fname;
    stringstream tmpStr;
    fname = "./outputs/";
    if (full)
        fname += "F";
    fname += "vnpt_y-";
    tmpStr << number;
    fname += tmpStr.str();
    fname += ".dat";  
    
    string fname2;
    stringstream tmpStr2;
    fname2 = "./outputs/";
    if (full)
        fname2 += "F";
    fname2 += "vnpt_eta-";
    tmpStr2 << number;
    fname2 += tmpStr.str();
    fname2 += ".dat"; 
    
    // Open output file for vn
    ofstream outfilevn;
    outfilevn.open(fname.c_str());
    
    ofstream outfilevn2;
    outfilevn2.open(fname2.c_str());
    
    outfilevn << "#pt  dN/ptdYdptdphi  v1cos  v1sin  v2cos  v2sin  "
              << "v3cos  v3sin  v4cos  v4sin  v5cos  v5sin  v6cos  v6sin  "
              << "v7cos  v7sin" << endl;
    
    outfilevn2 << "#pt  dN/ptdYdptdphi  v1cos  v1sin  v2cos  v2sin  "
               << "v3cos  v3sin  v4cos  v4sin  v5cos  v5sin  v6cos  v6sin  "
               << "v7cos  v7sin" << endl;
    
    double vn[nharmonics][2];

    double eta = 0.0;
    // Loop over pT
    for (int ipt = 0; ipt < npt; ipt++) {
        double pt=particleList[j].pt[ipt];
        pt_and_rapidity_integrated_flow(DATA, number, pt, pt, 1,
                                        eta, eta, vn);  // in rapidity

        //Output result
        outfilevn << pt;
        outfilevn << "  " << vn[0][0]/pt/(2*M_PI);
        for (int i = 1;i<nharmonics;i++) 
            for (int k =0;k<2;k++) 
                outfilevn << "  " << vn[i][k];
        outfilevn << endl;
 
        pt_and_rapidity_integrated_flow(DATA, number, pt, pt, 0,
                                        eta, eta, vn);  // pseudo-rapidity 
        outfilevn2 << pt;
        outfilevn2 << "  " << vn[0][0]/pt/(2*M_PI);
        for (int i = 1;i<nharmonics;i++) 
            for (int k =0;k<2;k++) 
                outfilevn2 << "  " << vn[i][k];
        outfilevn2 << endl;
    }
    
    // Close file
    outfilevn.close();
    outfilevn2.close();
}

// Output yield dN/ptdydpt and v_n for |y|<1 as a function of pT
void Freeze::OutputDifferentialFlowNearMidrapidity(
        InitData *DATA, int number, int full) {
    // Define index j used in particleList[j]
    int j = partid[MHALF+number];
    int npt = particleList[j].npt;
    
    double y_min = DATA->dNdyptdpt_y_min;
    double y_max = DATA->dNdyptdpt_y_max;
    double eta_min = DATA->dNdyptdpt_eta_min;
    double eta_max = DATA->dNdyptdpt_eta_max;
    
    music_message << "Calculating flow near midrapidity for particle "
                  << number;
    music_message.flush("info");
    
    //Set output file name
    stringstream fname;
    if (full) {
        fname << "./outputs/Fvnpt-" << number << "_y_"
              << y_min << "_" << y_max << ".dat";
    } else {
        fname << "./outputs/vnpt-" << number << "_y_"
              << y_min << "_" << y_max << ".dat";
    }
    
    stringstream fname2;
    if (full) {
        fname2 << "./outputs/Fvnpt-" << number << "_eta_"
               << eta_min << "_" << eta_max << ".dat";
    } else {
        fname2 << "./outputs/vnpt-" << number << "_eta_"
               << eta_min << "_" << eta_max << ".dat";
    }
    
    //Open output file for vn
    ofstream outfilevn;
    outfilevn.open(fname.str().c_str());
    
    ofstream outfilevn2;
    outfilevn2.open(fname2.str().c_str());

    //Set the format of the output
    outfilevn << "#pt  dN/ptdYdptdphi  v1cos  v1sin  v2cos  v2sin  "
              << "v3cos  v3sin  v4cos  v4sin  v5cos  v5sin  v6cos  v6sin  "
              << "v7cos  v7sin" << endl;   // in rapidity
    
    outfilevn2 << "#pt  dN/ptdYdptdphi  v1cos  v1sin  v2cos  v2sin  "
               << "v3cos  v3sin  v4cos  v4sin  v5cos  v5sin  v6cos  v6sin  "
               << "v7cos  v7sin" << endl;  // in pseudo-rapidity

    double vn[nharmonics][2];
    
    // Loop over pT
    for (int ipt = 0;ipt < npt; ipt++) {
        double pt = particleList[j].pt[ipt];
        
        pt_and_rapidity_integrated_flow(DATA, number, pt, pt, 1,
                                        y_min, y_max, vn);  // rapidity

        // Output result
        outfilevn << pt;
        outfilevn << "  " << vn[0][0]/pt/(y_max - y_min)/(2*M_PI);
        for (int i = 1;i<nharmonics;i++)
            for (int k =0;k<2;k++) 
                outfilevn << "  " << vn[i][k];
        outfilevn << endl;
        
        // pseudo-rapidity
        pt_and_rapidity_integrated_flow(DATA, number, pt, pt, 0,
                                        eta_min, eta_max, vn);
        // Output result
        outfilevn2 << pt;
        outfilevn2 << "  " << vn[0][0]/pt/(eta_max - eta_min)/(2*M_PI);
        for (int i = 1;i<nharmonics;i++)
            for (int k =0;k<2;k++) 
                outfilevn2 << "  " << vn[i][k];
        outfilevn2 << endl;
    }
    
    // Close file
    outfilevn.close();
    outfilevn2.close();
}


// Output yield dN/deta and v_n integrated over pT,
// as a function of rapidity and pseudorapidity
void Freeze::OutputIntegratedFlow_vs_y(
        InitData *DATA, int number, int full, double pT_min, double pT_max) {
    music_message << "Calculating flow integrated between " << pT_min
         << " < p_T < "
         << pT_max << " vs. pseudorapidity for particle " << number;
    music_message.flush("info");
   
    // Set output file name
    stringstream fname;
    if (full) {
       fname << "./outputs/Fvn_y-"
             << number << "_pT_" << pT_min << "_" << pT_max << ".dat";
    } else {
       fname << "./outputs/vn_y-"
             << number << "_pT_" << pT_min << "_" << pT_max << ".dat";
    }
    
    stringstream fname2;
    if (full) {
       fname2 << "./outputs/Fvn_eta-"
              << number << "_pT_" << pT_min << "_" << pT_max << ".dat";
    } else {
       fname2 << "./outputs/vn_eta-"
              << number << "_pT_" << pT_min << "_" << pT_max << ".dat";
    }
    
    // Open output file for vn
    ofstream outfilevn;
    outfilevn.open(fname.str().c_str());
    
    ofstream outfilevn2;
    outfilevn2.open(fname2.str().c_str());
    
    // header of the files
    outfilevn << "#y  dN/dy dummy  v1cos  v1sin  v2cos  v2sin  v3cos  v3sin  "
              << "v4cos  v4sin  v5cos  v5sin  v6cos  v6sin  v7cos  v7sin\n";
    outfilevn2 << "#eta  dN/deta dummy  v1cos  v1sin  v2cos  v2sin  v3cos  v3sin  "
               << "v4cos  v4sin  v5cos  v5sin  v6cos  v6sin  v7cos  v7sin\n";
   
    double vn[nharmonics][2];
    
    double y_min = DATA->dNdy_y_min;
    double y_max = DATA->dNdy_y_max;
    double eta_min = DATA->dNdy_eta_min;
    double eta_max = DATA->dNdy_eta_max;
    int nrap = DATA->dNdy_nrap;
    double dy = (y_max - y_min)/(nrap - 1);
    double deta = (eta_max - eta_min)/(nrap - 1);
    // Loop over eta or y
    for (int irap = 0; irap < nrap; irap++) {
        double y_local = y_min + irap*dy;
        double eta_local = eta_min + irap*deta;

        // rapidity
        pt_and_rapidity_integrated_flow(DATA, number, pT_min, pT_max, 1,
                                        y_local, y_local, vn);
   
        // Output result
        outfilevn << std::scientific << std::setprecision(8) << y_local;
        for (int i = 0; i < nharmonics; i++) {
            for (int k = 0; k < 2; k++) {
                outfilevn << "  " << vn[i][k];
            }
        }
        outfilevn << endl;
    
        // pseudo-rapidity
        pt_and_rapidity_integrated_flow(DATA, number, pT_min, pT_max, 0,
                                        eta_local, eta_local, vn);
    
        outfilevn2 << std::scientific << std::setprecision(8) << eta_local;
        for (int i = 0; i < nharmonics; i++) {
            for (int k = 0; k < 2; k++) {
                outfilevn2 << "  " << vn[i][k];
            }
        }
        outfilevn2 << endl;
    }
   
    // Close file
    outfilevn.close();
    outfilevn2.close();
}

// Output yield dN/deta and v_n integrated over pT and rapidity
void Freeze::Output_charged_IntegratedFlow(
        InitData *DATA, double pT_min, double pT_max,
        double eta_min, double eta_max) {
    music_message << "Calculating charged hadron flow integrated between "
         << pT_min << " < p_T < " << pT_max << " and "
         << eta_min << " < eta < " << eta_max << "...";
    music_message.flush("info");
    
    // Set output file name
    stringstream fname;
    fname << "./outputs/vnch_pT_" << pT_min << "_" << pT_max 
          << "_eta_" << eta_min << "_" << eta_max << ".dat";
    
    // Open output file for vn
    ofstream outfile;
    outfile.open(fname.str().c_str());
    
    // header of the files
    outfile << "#dN/deta  v1cos  v1sin  v2cos  v2sin  v3cos  v3sin  "
            << "v4cos  v4sin  v5cos  v5sin  v6cos  v6sin  v7cos  v7sin\n";
    
    double dNch = get_Nch(DATA, pT_min, pT_max, 0, eta_min, eta_max);
    outfile << std::scientific << std::setw(18) << std::setprecision(8)
            << dNch << "  ";
    int max_flow_order = 7;
    for (int iorder = 1; iorder < max_flow_order; iorder++) {
        double *vn_temp = new double [2];
        get_vn_ch(DATA, pT_min, pT_max, 0, eta_min, eta_max, iorder, vn_temp);
        outfile << vn_temp[0] << "  " << vn_temp[1] << "  ";
        delete [] vn_temp;
    }
    outfile << endl;
    outfile.close();
}

// Output pT and phi integrated charged hadron spectra as a function of eta
void Freeze::Output_charged_hadrons_eta_differential_spectra(
        InitData *DATA, int full, double pT_min, double pT_max) {
    music_message << "Collecting charged hadrons dN/deta and vn vs eta...";
    music_message.flush("info");
    stringstream tmpStr;
    tmpStr << "./outputs/vnchdeta_pT_" << pT_min << "_" << pT_max << ".dat";
    
    ofstream outfile;
    outfile.open(tmpStr.str().c_str(), std::ios::trunc);
    
    outfile << "#eta  dNch/deta  v1cos  v1sin  v2cos  v2sin  v3cos  v3sin  "
            << "v4cos  v4sin  v5cos  v5sin  v6cos  v6sin  v7cos  v7sin\n";

    double eta_min = DATA->dNdy_eta_min;
    double eta_max = DATA->dNdy_eta_max;
    int nrap = DATA->dNdy_nrap;
    double deta = (eta_max - eta_min)/(nrap - 1);
    int max_flow_order = 7;
    for (int ieta = 0; ieta < nrap; ieta++) {
        double eta_local = eta_min + ieta*deta;
        double dNch = get_Nch(DATA, pT_min, pT_max, 0, eta_local, eta_local);
        outfile << eta_local << "  "  << dNch << "  ";
        for (int iorder = 1; iorder < max_flow_order; iorder++) {
            double *vn_temp = new double [2];
            get_vn_ch(DATA, pT_min, pT_max, 0, eta_local, eta_local,
                      iorder, vn_temp);
            outfile << vn_temp[0] << "  " << vn_temp[1] << "  ";
            delete [] vn_temp;
        }
        outfile << endl;
    }
    outfile.close();
}

// output eta and phi integerated charged hadron spectra and vn as
// functions of pT
void Freeze::Output_charged_hadrons_pT_differential_spectra(
        InitData *DATA, int full, double eta_min, double eta_max) {
    double pT_min = particleList[1].pt[0];
    double pT_max = particleList[1].pt[particleList[1].npt-1];
    int npT = 30;
    double dpT = (pT_max - pT_min)/(npT - 1);
    music_message << "Collecting charged hadrons dN/detapTdpT and vn(pT) vs pT...";
    music_message.flush("info");

    stringstream tmpStr;
    tmpStr << "./outputs/vnchpT_eta_" << eta_min << "_" << eta_max << ".dat";
    ofstream outfile;
    outfile.open(tmpStr.str().c_str(), std::ios::trunc);
    // header
    outfile << "#pT  dNch/detapTdpTdphi  v1cos  v1sin  v2cos  v2sin  "
            << "v3cos  v3sin  v4cos  v4sin  v5cos  v5sin  v6cos  v6sin  "
            << "v7cos  v7sin" << endl;

    int max_flow_order = 7;
    for (int ipT = 0; ipT < npT; ipT++) {
        double pT_local = pT_min + ipT*dpT;
        double dNch = get_Nch(DATA, pT_local, pT_local, 0, eta_min, eta_max);
        outfile << pT_local << "  "  << dNch/pT_local/(2*M_PI) << "  ";
        for (int iorder = 1; iorder < max_flow_order; iorder++) {
            double *vn_temp = new double [2];
            get_vn_ch(DATA, pT_local, pT_local, 0, eta_min, eta_max,
                      iorder, vn_temp);
            outfile << vn_temp[0] << "  " << vn_temp[1] << "  ";
            delete [] vn_temp;
        }
        outfile << endl;
    }
    outfile.close();
}

// Return yield in specified range of phase space
double Freeze::get_N(InitData *DATA, int number, double minpt, double maxpt,
                     int yflag, double minrap, double maxrap) {
    double vn[nharmonics][2];
    pt_and_rapidity_integrated_flow(DATA, number, minpt, maxpt, yflag,
                                    minrap, maxrap, vn);
    return(vn[0][0]);
}

// Return charged hadron yield in specified range of phase space
// You'll almost certainly want to set yflat=0 to calculate 
// over a fixed range in eta
// if minrap==maxrap, returns dN/eta.  If minpt==maxpt, dN/dpt.
// If both, dN/dpt/deta.  Otherwise, total yield N
double Freeze::get_Nch(InitData *DATA, double minpt, double maxpt, int yflag,
                       double minrap, double maxrap) {
    double N = 0;
    for (int k = 0; k < charged_hadron_list_length; k++) {
        int number = charged_hadron_list[k];
        N+= get_N(DATA, number, minpt, maxpt, yflag, minrap, maxrap);
    }
    return N;
}


// Return charged hadron vn for specified range of phase space
double Freeze::get_vn_ch(InitData *DATA, double minpt, double maxpt, int yflag,
                         double minrap, double maxrap, int n,
                         double* vn_results) {
    if (n > nharmonics - 1) {
        music_message << "harmonic (" << n << ") too large.  "
             << "Must increase nharmonics in freeze_pseudo.cpp";
        music_message.flush("error");
        exit(1);
    }
    double numr = 0.;  // real part (x projection, \sum N*v_n*cos(Psi_n))
    double numi = 0.;  // imaginary part (y projection, \sum N*v_n*sin(Psi_n))
    double den = 0.;   // denominator (\sum N)
    for (int k=0; k< charged_hadron_list_length; k++) {
        int number = charged_hadron_list[k];
        double vn[nharmonics][2];
        pt_and_rapidity_integrated_flow(DATA, number, minpt, maxpt, yflag,
                                        minrap, maxrap, vn);
        numr+= vn[0][0]*vn[n][0];
        numi+= vn[0][0]*vn[n][1];
        den+= vn[0][0];
    }
    vn_results[0] = numr/den;
    vn_results[1] = numi/den;
    return sqrt(numr*numr+numi*numi)/den;
}


// Return <pt> in specified range of phase space
double Freeze::get_meanpt(InitData *DATA, int number,
                          double minpt, double maxpt,
                          int yflag, double minrap, double maxrap) {
    int j = partid[MHALF+number];
    int npt = particleList[j].npt;
    if (minpt < particleList[j].pt[0]) {
        music_message << "Error: called out of range pt in get_meanpt, "
             << minpt << " < minimum " << particleList[j].pt[0];
        music_message.flush("error");
        exit(1);
    }
    if (maxpt > particleList[j].pt[npt-1]) {
        music_message << "Error: called out of range pt in get_meanpt, "
             << maxpt << " > maximum " << particleList[j].pt[npt-1];
        music_message.flush("error");
        exit(1);
    }
    if (minpt > maxpt) {
        music_message << "Error in get_meanpt:  "
             << "minpt must be less than or equal to maxpt";
        music_message.flush("error");
        exit(1);
    }
  
    double ptdndpt[ptsize];
    double dndpt[ptsize];
    for (int ipt = 0; ipt < npt; ipt++) {
        double pt = particleList[j].pt[ipt];
        dndpt[ipt] = get_N(DATA, number, pt, pt, yflag, minrap, maxrap);
        ptdndpt[ipt] = pt*dndpt[ipt];
    }
    gsl_interp_accel *numacc = gsl_interp_accel_alloc ();
    gsl_spline *numspline = gsl_spline_alloc (gsl_interp_linear, npt);
    gsl_spline_init (numspline, particleList[j].pt ,ptdndpt , npt);
    
    gsl_interp_accel *denacc = gsl_interp_accel_alloc ();
    gsl_spline *denspline = gsl_spline_alloc (gsl_interp_linear, npt);
    gsl_spline_init (denspline, particleList[j].pt ,dndpt , npt);
    
    double num = gsl_spline_eval_integ(numspline, minpt, maxpt, numacc);
    double den = gsl_spline_eval_integ(denspline, minpt, maxpt, denacc);
    
    gsl_spline_free (numspline);
    gsl_interp_accel_free (numacc);
    gsl_spline_free (denspline);
    gsl_interp_accel_free (denacc);
    
    return num/den;
}


// Return charged hadron Psi_n for specified range of phase space
double Freeze::get_psi_n_ch(InitData *DATA, double minpt, double maxpt,
                            int yflag, double minrap, double maxrap, int n) {
    if (n > nharmonics-1) {
        music_message << "harmonic (" << n << ") too large.  "
             << "Must increase nharmonics in freeze_pseudo.cpp";
        music_message.flush("error");
        exit(1);
    }
    double numr = 0.;  // real part (x projection, \sum N*v_n*cos(Psi_n))
    double numi = 0.;  // imaginary part (y projection, \sum N*v_n*sin(Psi_n))
    for (int k = 0; k < charged_hadron_list_length; k++) {
        int number = charged_hadron_list[k];
        double vn[nharmonics][2];
        pt_and_rapidity_integrated_flow(DATA, number, minpt, maxpt,
                                        yflag, minrap, maxrap, vn);
        numr += vn[0][0]*vn[n][0];
        numi += vn[0][0]*vn[n][1];
    }
    return atan2(numi, numr)/n;
}

// Return v_n for specified range of phase space
double Freeze::get_vn(InitData *DATA, int number, double minpt, double maxpt,
                      int yflag, double minrap, double maxrap, int n) {
    if (n > nharmonics - 1) {
        music_message << "harmonic (" << n << ") too large.  "
             << "Must increase nharmonics in freeze_pseudo.cpp";
        music_message.flush("error");
        exit(1);
    }
    double vn[nharmonics][2];
    pt_and_rapidity_integrated_flow(DATA, number, minpt, maxpt,
                                    yflag, minrap, maxrap, vn);
    return sqrt(vn[n][0]*vn[n][0] + vn[n][1]*vn[n][1]);
}

// Return Psi_n for specified range of phase space
double Freeze::get_psi_n(InitData *DATA, int number,
                         double minpt, double maxpt,
                         int yflag, double minrap, double maxrap, int n) {
    if (n > nharmonics - 1) {
        music_message << "harmonic (" << n << ") too large.  "
             << "Must increase nharmonics in freeze_pseudo.cpp";
        music_message.flush("error");
        exit(1);
    }
    double vn[nharmonics][2];
    pt_and_rapidity_integrated_flow(DATA, number, minpt, maxpt,
                                    yflag, minrap, maxrap, vn);
    return atan2(vn[n][1],vn[n][0])/n;
}


// Calculates v1 and psi_1 with pt-dependent weight (pt - <pt^2>/<pt>)
// in specified range of phase space
// ch == 0 calculates for particle "number",
// while ch > 0 calculates for all charged hadrons
void Freeze::weighted_v1(InitData *DATA, int number,
                         double minpt, double maxpt,
                         int yflag, double minrap, double maxrap,
                         double vn[2], int ch) {
    int j = partid[MHALF+number];
    int npt = particleList[j].npt;
    if (minpt < particleList[j].pt[0]) {
        music_message << "Error: called out of range pt in get_meanpt, "
             << minpt << " < minimum " << particleList[j].pt[0];
        music_message.flush("error");
        exit(1);
    }
    if (maxpt > particleList[j].pt[npt-1]) {
        music_message << "Error: called out of range pt in get_meanpt, "
             << maxpt << " > maximum " << particleList[j].pt[npt-1];
        music_message.flush("error");
        exit(1);
    }
    if (minpt > maxpt) {
        music_message << "Error in get_meanpt:  "
             << "minpt must be less than or equal to maxpt";
        music_message.flush("error");
        exit(1);
    }
  
    // start by finding <pt> and <pt^2>
    double ptptdndpt[ptsize];
    double ptdndpt[ptsize];
    double dndpt[ptsize];
    for (int ipt = 0; ipt < npt; ipt++) {
        double pt = particleList[j].pt[ipt];
        if (ch) {
            dndpt[ipt] = get_Nch(DATA, pt, pt, yflag, minrap, maxrap);
        } else {
            dndpt[ipt] = get_N(DATA, number, pt, pt, yflag, minrap, maxrap);
        }
        ptdndpt[ipt] = pt*dndpt[ipt];
        ptptdndpt[ipt] = pt*pt*dndpt[ipt];
    }
  
    gsl_interp_accel *numacc = gsl_interp_accel_alloc ();
    gsl_spline *numspline = gsl_spline_alloc (gsl_interp_linear, npt);
    gsl_spline_init (numspline, particleList[j].pt ,ptdndpt , npt);
    
    gsl_interp_accel *num2acc = gsl_interp_accel_alloc ();
    gsl_spline *num2spline = gsl_spline_alloc (gsl_interp_linear, npt);
    gsl_spline_init (num2spline, particleList[j].pt, ptptdndpt , npt);
    
    gsl_interp_accel *denacc = gsl_interp_accel_alloc ();
    gsl_spline *denspline = gsl_spline_alloc (gsl_interp_linear, npt);
    gsl_spline_init (denspline, particleList[j].pt ,dndpt , npt);
    
    double num2 = gsl_spline_eval_integ(num2spline, minpt, maxpt, num2acc);
    double num = gsl_spline_eval_integ(numspline, minpt, maxpt, numacc);
    double den = gsl_spline_eval_integ(denspline, minpt, maxpt, denacc);
    
    gsl_spline_free (numspline);
    gsl_interp_accel_free (numacc);
    gsl_spline_free (num2spline);
    gsl_interp_accel_free (num2acc);
    gsl_spline_free (denspline);
    gsl_interp_accel_free (denacc);
    
    double meanpt = num/den;
    double meanpt2 = num2/den;
    
    // next calculate integrated v1 and psi1 with weight pt - <pt^2>/<pt>
    double v1cos[ptsize];
    double v1sin[ptsize];
    double v1den[ptsize];
    for (int ipt = 0; ipt < npt; ipt++) {
        double pt = particleList[j].pt[ipt];
        double weight = (pt - meanpt2/meanpt)*dndpt[ipt];
        
        double v1;
        double psi1;
        double * vntemp = new double [2];
        if (ch) {
            v1 = get_vn_ch(DATA, pt, pt, yflag, minrap, maxrap, 1, vntemp);
        } else {
            v1 = get_vn(DATA, number, pt, pt, yflag, minrap, maxrap, 1);
        }
        if (ch) {
            psi1 = get_psi_n_ch(DATA, pt, pt, yflag, minrap, maxrap, 1);
        } else {
            psi1 = get_psi_n(DATA, number, pt, pt, yflag, minrap, maxrap, 1);
        }
        v1cos[ipt] = weight*v1*cos(psi1);
        v1sin[ipt] = weight*v1*sin(psi1);
        v1den[ipt] = weight;
        delete [] vntemp;
    }
    
    gsl_spline_init (numspline, particleList[j].pt ,v1cos , npt);
    gsl_spline_init (num2spline, particleList[j].pt, v1sin , npt);
    gsl_spline_init (denspline, particleList[j].pt ,v1den , npt);
    
    num2 = gsl_spline_eval_integ(num2spline, minpt, maxpt, numacc);
    num = gsl_spline_eval_integ(numspline, minpt, maxpt, numacc);
    den = gsl_spline_eval_integ(denspline, minpt, maxpt, denacc);
    
    gsl_spline_free (numspline);
    gsl_interp_accel_free (numacc);
    gsl_spline_free (num2spline);
    gsl_interp_accel_free (num2acc);
    gsl_spline_free (denspline);
    gsl_interp_accel_free (denacc);
  
    vn[0] = sqrt(num*num+num2*num2)/den;
    vn[1] = atan2(num2,num);
}


double Freeze::get_weighted_v1(InitData *DATA, int number,
                               double minpt, double maxpt,
                               int yflag, double minrap, double maxrap,
                               int ch) {
    double vn[2];
    weighted_v1(DATA, number, minpt, maxpt, yflag, minrap, maxrap, vn, ch);
    return vn[0];
}

double Freeze::get_weighted_psi1(InitData *DATA, int number,
                                 double minpt, double maxpt,
                                 int yflag, double minrap, double maxrap,
                                 int ch) {
    double vn[2];
    weighted_v1(DATA, number, minpt, maxpt, yflag, minrap, maxrap, vn, ch);
    return vn[1];
}

void Freeze::load_deltaf_qmu_coeff_table(string filename) {
    std::ifstream table(filename.c_str());
    deltaf_qmu_coeff_table_length_T = 150;
    deltaf_qmu_coeff_table_length_mu = 100;
    delta_qmu_coeff_table_T0 = 0.05;
    delta_qmu_coeff_table_mu0 = 0.0;
    delta_qmu_coeff_table_dT = 0.001;
    delta_qmu_coeff_table_dmu = 0.007892;
    deltaf_qmu_coeff_tb = new double* [deltaf_qmu_coeff_table_length_T];
    for (int i = 0; i < deltaf_qmu_coeff_table_length_T; i++) {
        deltaf_qmu_coeff_tb[i] = new double [deltaf_qmu_coeff_table_length_mu];
    }

    double dummy;
    for (int j = 0; j < deltaf_qmu_coeff_table_length_mu; j++) {
        for (int i = 0; i < deltaf_qmu_coeff_table_length_T; i++) {
            table >> dummy >> dummy >> deltaf_qmu_coeff_tb[i][j];
        }
    }
    table.close();
}

void Freeze::load_deltaf_qmu_coeff_table_14mom(string filename) {
    std::ifstream table(filename.c_str());
    deltaf_coeff_table_14mom_length_T = 190;
    deltaf_coeff_table_14mom_length_mu = 160;
    delta_coeff_table_14mom_T0 = 0.01;
    delta_coeff_table_14mom_mu0 = 0.0;
    delta_coeff_table_14mom_dT = 0.001;
    delta_coeff_table_14mom_dmu = 0.005;

    deltaf_coeff_tb_14mom_DPi = new double* [deltaf_coeff_table_14mom_length_T];
    deltaf_coeff_tb_14mom_BPi = new double* [deltaf_coeff_table_14mom_length_T];
    deltaf_coeff_tb_14mom_BPitilde = 
                            new double* [deltaf_coeff_table_14mom_length_T];
    deltaf_coeff_tb_14mom_DV = new double* [deltaf_coeff_table_14mom_length_T];
    deltaf_coeff_tb_14mom_BV = new double* [deltaf_coeff_table_14mom_length_T];
    deltaf_coeff_tb_14mom_Bpi_shear =
                            new double* [deltaf_coeff_table_14mom_length_T];
    for (int i = 0; i < deltaf_coeff_table_14mom_length_T; i++) {
        deltaf_coeff_tb_14mom_DPi[i] =
                            new double[deltaf_coeff_table_14mom_length_mu];
        deltaf_coeff_tb_14mom_BPi[i] =
                            new double[deltaf_coeff_table_14mom_length_mu];
        deltaf_coeff_tb_14mom_BPitilde[i] =
                            new double[deltaf_coeff_table_14mom_length_mu];
        deltaf_coeff_tb_14mom_DV[i] =
                            new double[deltaf_coeff_table_14mom_length_mu];
        deltaf_coeff_tb_14mom_BV[i] =
                            new double[deltaf_coeff_table_14mom_length_mu];
        deltaf_coeff_tb_14mom_Bpi_shear[i] =
                            new double[deltaf_coeff_table_14mom_length_mu];
    }

    double dummy;
    for (int i = 0; i < deltaf_coeff_table_14mom_length_T; i++) {
        for (int j = 0; j < deltaf_coeff_table_14mom_length_mu; j++) {
            table >> dummy >> dummy >> deltaf_coeff_tb_14mom_DPi[i][j]
                  >> deltaf_coeff_tb_14mom_BPi[i][j] 
                  >> deltaf_coeff_tb_14mom_BPitilde[i][j]
                  >> deltaf_coeff_tb_14mom_DV[i][j] 
                  >> deltaf_coeff_tb_14mom_BV[i][j]
                  >> deltaf_coeff_tb_14mom_Bpi_shear[i][j];
        }
    }
    table.close();

    // convert units
    double hbarc3 = hbarc*hbarc*hbarc;
    double hbarc4 = hbarc3*hbarc;
    for (int i = 0; i < deltaf_coeff_table_14mom_length_T; i++) {
        for (int j = 0; j < deltaf_coeff_table_14mom_length_mu; j++) {
            deltaf_coeff_tb_14mom_DPi[i][j] *= hbarc4;          // fm^4/GeV
            deltaf_coeff_tb_14mom_BPi[i][j] *= hbarc4;          // fm^4/(GeV^2)
            deltaf_coeff_tb_14mom_BPitilde[i][j] *= hbarc4;     // fm^4/(GeV^2)
            deltaf_coeff_tb_14mom_DV[i][j] *= hbarc3;           // fm^3/GeV
            deltaf_coeff_tb_14mom_BV[i][j] *= hbarc3;           // fm^3/(GeV^2)
            deltaf_coeff_tb_14mom_Bpi_shear[i][j] *= hbarc4;    // fm^4/(GeV^2)
        }
    }
}

double Freeze::get_deltaf_qmu_coeff(double T, double muB) {
    int idx_T = static_cast<int>((T - delta_qmu_coeff_table_T0)
                                 /delta_qmu_coeff_table_dT);
    int idx_mu = static_cast<int>((muB - delta_qmu_coeff_table_mu0)
                                  /delta_qmu_coeff_table_dmu);
    double x_fraction = ((T - delta_qmu_coeff_table_T0)
                         /delta_qmu_coeff_table_dT - idx_T);
    double y_fraction = ((muB - delta_qmu_coeff_table_mu0)
                         /delta_qmu_coeff_table_dmu - idx_mu);

    // avoid overflow: return a large number so that delta f = 0
    if (idx_mu > deltaf_qmu_coeff_table_length_mu - 2) {
        return(1e30);
    }
    if (idx_T > deltaf_qmu_coeff_table_length_T - 2) {
        return(1e30);
    }
    
    // avoid underflow: return a large number so that delta f = 0
    if (idx_mu < 0) {
        return(1e30);
    }
    if (idx_T < 0) {
        return(1e30);
    }
    
    double f1 = deltaf_qmu_coeff_tb[idx_T][idx_mu];
    double f2 = deltaf_qmu_coeff_tb[idx_T][idx_mu+1];
    double f3 = deltaf_qmu_coeff_tb[idx_T+1][idx_mu+1];
    double f4 = deltaf_qmu_coeff_tb[idx_T+1][idx_mu];

    double coeff = f1*(1. - x_fraction)*(1. - y_fraction) 
                   + f2*(1. - x_fraction)*y_fraction
                   + f3*x_fraction*y_fraction
                   + f4*x_fraction*(1. - y_fraction);
    return(coeff);
}

double Freeze::get_deltaf_coeff_14moments(double T, double muB, double type) {
    int idx_T = static_cast<int>((T - delta_coeff_table_14mom_T0)
                                 /delta_coeff_table_14mom_dT);
    int idx_mu = static_cast<int>((muB - delta_coeff_table_14mom_mu0)
                                  /delta_coeff_table_14mom_dmu);
    double x_fraction = ((T - delta_coeff_table_14mom_T0)
                         /delta_coeff_table_14mom_dT - idx_T);
    double y_fraction = ((muB - delta_coeff_table_14mom_mu0)
                         /delta_coeff_table_14mom_dmu - idx_mu);

    double **deltaf_table = NULL;
    if (type == 0) {
       deltaf_table = deltaf_coeff_tb_14mom_DPi;
    } else if (type == 1) {
       deltaf_table = deltaf_coeff_tb_14mom_BPi;
    } else if (type == 2) {
       deltaf_table = deltaf_coeff_tb_14mom_BPitilde;
    } else if (type == 3) {
       deltaf_table = deltaf_coeff_tb_14mom_DV;
    } else if (type == 4) {
       deltaf_table = deltaf_coeff_tb_14mom_BV;
    } else if (type == 5) {
       deltaf_table = deltaf_coeff_tb_14mom_Bpi_shear;
    } else {
       music_message << "Freeze::get_deltaf_coeff_14moments: unknown type: " 
                     << type;
       music_message.flush("error");
       exit(-1);
    }
    
    double f1 = deltaf_table[idx_T][idx_mu];
    double f2 = deltaf_table[idx_T][idx_mu+1];
    double f3 = deltaf_table[idx_T+1][idx_mu+1];
    double f4 = deltaf_table[idx_T+1][idx_mu];

    double coeff = f1*(1. - x_fraction)*(1. - y_fraction) 
                   + f2*(1. - x_fraction)*y_fraction
                   + f3*x_fraction*y_fraction
                   + f4*x_fraction*(1. - y_fraction);
    return(coeff);
}

void Freeze::getbulkvisCoefficients(double Tdec, double* bulkvisCoefficients) {
    double Tdec_fm = Tdec/hbarc;  // [1/fm]
    double Tdec_fm_power[11];     // cache the polynomial power of Tdec_fm
    Tdec_fm_power[1] = Tdec_fm;
    for (int ipower = 2; ipower < 11; ipower++)
        Tdec_fm_power[ipower] = Tdec_fm_power[ipower-1]*Tdec_fm;
    if (bulk_deltaf_kind == 0) {       // 14 moment expansion
        // parameterization for mu = 0
        // B0[fm^3/GeV^3]
        bulkvisCoefficients[0] = (
                exp(-15.04512474*Tdec_fm + 11.76194266)/pow(hbarc, 3)); 
        // D0 [fm^3/GeV^2]
        bulkvisCoefficients[1] = (
                exp( -12.45699277*Tdec_fm + 11.4949293)/hbarc/hbarc);  
        // E0 [fm^3/GeV^3]
        bulkvisCoefficients[2] = (
                -exp(-14.45087586*Tdec_fm + 11.62716548)/pow(hbarc, 3));  
    } else if(bulk_deltaf_kind == 1) {     // relaxation type 1
        // parameterization from JF
        // A Polynomial fit to each coefficient -- temperature in fm^-1
        // Both fits are reliable between T=100 -- 180 MeV
        // do not trust it beyond
        bulkvisCoefficients[0] = (  642096.624265727
                                  - 8163329.49562861*Tdec_fm_power[1]
                                  + 47162768.4292073*Tdec_fm_power[2]
                                  - 162590040.002683*Tdec_fm_power[3]
                                  + 369637951.096896*Tdec_fm_power[4]
                                  - 578181331.809836*Tdec_fm_power[5]
                                  + 629434830.225675*Tdec_fm_power[6]
                                  - 470493661.096657*Tdec_fm_power[7]
                                  + 230936465.421*Tdec_fm_power[8]
                                  - 67175218.4629078*Tdec_fm_power[9]
                                  + 8789472.32652964*Tdec_fm_power[10]);

        bulkvisCoefficients[1] = (  1.18171174036192
                                  - 17.6740645873717*Tdec_fm_power[1]
                                  + 136.298469057177*Tdec_fm_power[2]
                                  - 635.999435106846*Tdec_fm_power[3]
                                  + 1918.77100633321*Tdec_fm_power[4]
                                  - 3836.32258307711*Tdec_fm_power[5]
                                  + 5136.35746882372*Tdec_fm_power[6]
                                  - 4566.22991441914*Tdec_fm_power[7]
                                  + 2593.45375240886*Tdec_fm_power[8]
                                  - 853.908199724349*Tdec_fm_power[9]
                                  + 124.260460450113*Tdec_fm_power[10]);
    } else if (bulk_deltaf_kind == 2) {     // relaxation type 2
        // A Polynomial fit to each coefficient -- temperature in fm^-1
        // Both fits are reliable between T=100 -- 180 MeV
        // do not trust it beyond
        bulkvisCoefficients[0] = (
                21091365.1182649 - 290482229.281782*Tdec_fm_power[1]
              + 1800423055.01882*Tdec_fm_power[2]
              - 6608608560.99887*Tdec_fm_power[3]
              + 15900800422.7138*Tdec_fm_power[4]
              - 26194517161.8205*Tdec_fm_power[5]
              + 29912485360.2916*Tdec_fm_power[6]
              - 23375101221.2855*Tdec_fm_power[7]
              + 11960898238.0134*Tdec_fm_power[8]
              - 3618358144.18576*Tdec_fm_power[9]
              + 491369134.205902*Tdec_fm_power[10]);

        bulkvisCoefficients[1] = (
                4007863.29316896 - 55199395.3534188*Tdec_fm_power[1]
              + 342115196.396492*Tdec_fm_power[2]
              - 1255681487.77798*Tdec_fm_power[3]
              + 3021026280.08401*Tdec_fm_power[4]
              - 4976331606.85766*Tdec_fm_power[5]
              + 5682163732.74188*Tdec_fm_power[6]
              - 4439937810.57449*Tdec_fm_power[7]
              + 2271692965.05568*Tdec_fm_power[8]
              - 687164038.128814*Tdec_fm_power[9]
              + 93308348.3137008*Tdec_fm_power[10]);
    } else if (bulk_deltaf_kind == 3) {    // relaxation type 3
        bulkvisCoefficients[0] = (
                160421664.93603 - 2212807124.97991*Tdec_fm_power[1]
              + 13707913981.1425*Tdec_fm_power[2]
              - 50204536518.1767*Tdec_fm_power[3]
              + 120354649094.362*Tdec_fm_power[4]
              - 197298426823.223*Tdec_fm_power[5]
              + 223953760788.288*Tdec_fm_power[6]
              - 173790947240.829*Tdec_fm_power[7]
              + 88231322888.0423*Tdec_fm_power[8]
              - 26461154892.6963*Tdec_fm_power[9]
              + 3559805050.19592*Tdec_fm_power[10]);
        bulkvisCoefficients[1] = (
                33369186.2536556 - 460293490.420478*Tdec_fm_power[1]
              + 2851449676.09981*Tdec_fm_power[2]
              - 10443297927.601*Tdec_fm_power[3]
              + 25035517099.7809*Tdec_fm_power[4]
              - 41040777943.4963*Tdec_fm_power[5]
              + 46585225878.8723*Tdec_fm_power[6]
              - 36150531001.3718*Tdec_fm_power[7]
              + 18353035766.9323*Tdec_fm_power[8]
              - 5504165325.05431*Tdec_fm_power[9]
              + 740468257.784873*Tdec_fm_power[10]);
    } else if (bulk_deltaf_kind == 4)  {    // relaxation type 4
        bulkvisCoefficients[0] = (
                1167272041.90731 - 16378866444.6842*Tdec_fm_power[1]
              + 103037615761.617*Tdec_fm_power[2]
              - 382670727905.111*Tdec_fm_power[3]
              + 929111866739.436*Tdec_fm_power[4]
              - 1540948583116.54*Tdec_fm_power[5]
              + 1767975890298.1*Tdec_fm_power[6]
              - 1385606389545*Tdec_fm_power[7]
              + 709922576963.213*Tdec_fm_power[8]
              - 214726945096.326*Tdec_fm_power[9]
              + 29116298091.9219*Tdec_fm_power[10]);
        bulkvisCoefficients[1] = (
                5103633637.7213 - 71612903872.8163*Tdec_fm_power[1]
              + 450509014334.964*Tdec_fm_power[2]
              - 1673143669281.46*Tdec_fm_power[3]
              + 4062340452589.89*Tdec_fm_power[4]
              - 6737468792456.4*Tdec_fm_power[5]
              + 7730102407679.65*Tdec_fm_power[6]
              - 6058276038129.83*Tdec_fm_power[7]
              + 3103990764357.81*Tdec_fm_power[8]
              - 938850005883.612*Tdec_fm_power[9]
              + 127305171097.249*Tdec_fm_power[10]);
    }
    return;
}
