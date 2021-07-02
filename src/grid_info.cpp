// Copyright Chun Shen @ 2014-2016
#include <string>
#include <iomanip>
#include <cmath>
#include <vector>

#include "util.h"
#include "grid_info.h"

using Util::hbarc;
using Util::small_eps;
using std::string;
using std::scientific;
using std::setw;
using std::setprecision;
using std::endl;
using std::ofstream;
using std::ostringstream;

Cell_info::Cell_info(const InitData &DATA_in, const EOS &eos_in) :
    DATA(DATA_in),
    eos(eos_in),
    u_derivative_helper(DATA_in, eos_in) {

    // read in tables for delta f coefficients
    if (DATA.turn_on_diff == 1) {
        if (DATA.deltaf_14moments == 1) {
            load_deltaf_qmu_coeff_table_14mom(
                    "tables/deltaf_coefficients_14moments.dat");
        } else {
            if (DATA.include_deltaf_qmu == 1)
                load_deltaf_qmu_coeff_table(
                        "tables/Coefficients_RTA_diffusion.dat");
        }
    }
    Pmu_edge_prev = {0.};
    outflow_flux = {0.};
}

Cell_info::~Cell_info() {
    if (DATA.turn_on_diff == 1) {
        if (DATA.deltaf_14moments == 1) {
            for (int i = 0; i < deltaf_coeff_table_14mom_length_T; i++) {
                delete [] deltaf_coeff_tb_14mom_DPi[i];
                delete [] deltaf_coeff_tb_14mom_BPi[i];
                delete [] deltaf_coeff_tb_14mom_BPitilde[i];
                delete [] deltaf_coeff_tb_14mom_DV[i];
                delete [] deltaf_coeff_tb_14mom_BV[i];
                delete [] deltaf_coeff_tb_14mom_Bpi_shear[i];
            }
            delete [] deltaf_coeff_tb_14mom_DPi;
            delete [] deltaf_coeff_tb_14mom_BPi;
            delete [] deltaf_coeff_tb_14mom_BPitilde;
            delete [] deltaf_coeff_tb_14mom_DV;
            delete [] deltaf_coeff_tb_14mom_BV;
            delete [] deltaf_coeff_tb_14mom_Bpi_shear;
        } else {
            if (DATA.include_deltaf_qmu == 1) {
                for (int i = 0; i < deltaf_qmu_coeff_table_length_T; i++)
                    delete [] deltaf_qmu_coeff_tb[i];
                delete [] deltaf_qmu_coeff_tb;
            }
        }
    }
}


//! This function outputs a header files for JF and Gojko's EM program
void Cell_info::Output_hydro_information_header() {
    string fname = "hydro_info_header_h";

    // Open output file
    ofstream outfile;
    outfile.open(fname.c_str());

    int grid_nx = ceil(
        (static_cast<double>(DATA.nx + 1))/DATA.output_evolution_every_N_x);
    int grid_ny = ceil(
        (static_cast<double>(DATA.ny + 1))/DATA.output_evolution_every_N_y);
    int grid_neta = ceil(
        (static_cast<double>(DATA.neta))/DATA.output_evolution_every_N_eta);

    outfile << "const int MUSIC_real_nx = " << grid_nx << ";" << endl;
    outfile << "const int MUSIC_real_ny = " << grid_ny << ";" << endl;
    outfile << "const int MUSIC_real_neta = " << grid_neta << ";" << endl;

    outfile << "const double MUSIC_tau0 = " << DATA.tau0 << ";" << endl;

    outfile << "const double MUSIC_dx = "
            << DATA.delta_x*DATA.output_evolution_every_N_x << ";" << endl;
    outfile << "const double MUSIC_dy = "
            << DATA.delta_y*DATA.output_evolution_every_N_y << ";" << endl;
    outfile << "const double MUSIC_deta = "
            << DATA.delta_eta*DATA.output_evolution_every_N_eta << ";"
            << endl;
    outfile << "const double MUSIC_dtau = "
            << DATA.output_evolution_every_N_timesteps*DATA.delta_tau << ";"
            << endl;

    outfile << "const bool MUSIC_with_shear_viscosity = "
            << ((DATA.viscosity_flag) && (DATA.turn_on_shear)) << ";\n";
    outfile << "const bool MUSIC_with_bulk_viscosity = "
            << ((DATA.viscosity_flag) && (DATA.turn_on_bulk)) << ";\n";
    outfile << "const bool MUSIC_with_rhob = " << (DATA.turn_on_rhob) << ";\n";
    outfile << "const bool MUSIC_with_diffusion = "
            << ((DATA.viscosity_flag) && (DATA.turn_on_diff)) << ";\n";

    outfile << "const bool MUSIC_outputBinaryEvolution="
            << DATA.outputBinaryEvolution << ";" << endl;

    outfile.close();
}


//! This function outputs hydro evolution file in binary format
void Cell_info::OutputEvolutionDataXYEta(SCGrid &arena, double tau) {
    const string out_name_xyeta = "evolution_xyeta.dat";
    const string out_name_W_xyeta =
                        "evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
    const string out_name_bulkpi_xyeta = "evolution_bulk_pressure_xyeta.dat";
    const string out_name_q_xyeta = "evolution_qmu_xyeta.dat";
    string out_open_mode;
    FILE *out_file_xyeta        = NULL;
    FILE *out_file_W_xyeta      = NULL;
    FILE *out_file_bulkpi_xyeta = NULL;
    FILE *out_file_q_xyeta      = NULL;

    // If it's the first timestep, overwrite the previous file
    if (tau == DATA.tau0) {
        out_open_mode = "w";
    } else {
        out_open_mode = "a";
    }
    // If we output in binary, set the mode accordingly
    if (0 == DATA.outputBinaryEvolution) {
        out_open_mode += "b";
    }
    out_file_xyeta = fopen(out_name_xyeta.c_str(), out_open_mode.c_str());
    if (DATA.turn_on_shear == 1) {
        out_file_W_xyeta = fopen(out_name_W_xyeta.c_str(),
                                 out_open_mode.c_str());
    }
    if (DATA.turn_on_bulk == 1) {
        out_file_bulkpi_xyeta = fopen(out_name_bulkpi_xyeta.c_str(),
                                      out_open_mode.c_str());
    }
    if (DATA.turn_on_diff == 1) {
        out_file_q_xyeta = fopen(out_name_q_xyeta.c_str(),
                                 out_open_mode.c_str());
    }
    const int n_skip_x   = DATA.output_evolution_every_N_x;
    const int n_skip_y   = DATA.output_evolution_every_N_y;
    const int n_skip_eta = DATA.output_evolution_every_N_eta;
    for (int ieta = 0; ieta < arena.nEta(); ieta += n_skip_eta) {
        double eta = 0.0;
        if (!DATA.boost_invariant) {
            eta = ((static_cast<double>(ieta))*(DATA.delta_eta)
                    - (DATA.eta_size)/2.0);
        }
        double cosh_eta = cosh(eta);
        double sinh_eta = sinh(eta);
        for (int iy = 0; iy < arena.nY(); iy += n_skip_y) {
            for (int ix = 0; ix < arena.nX(); ix += n_skip_x) {
                double e_local    = arena(ix, iy, ieta).epsilon;  // 1/fm^4
                double rhob_local = arena(ix, iy, ieta).rhob;     // 1/fm^3
                double p_local = eos.get_pressure(e_local, rhob_local);
                double utau = arena(ix, iy, ieta).u[0];
                double ux   = arena(ix, iy, ieta).u[1];
                double uy   = arena(ix, iy, ieta).u[2];
                double ueta = arena(ix, iy, ieta).u[3];
                double ut = utau*cosh_eta + ueta*sinh_eta;  // gamma factor
                double vx = ux/ut;
                double vy = uy/ut;
                double uz = ueta*cosh_eta + utau*sinh_eta;
                double vz = uz/ut;

                double T_local   = eos.get_temperature(e_local, rhob_local);
                double cs2_local = eos.get_cs2(e_local, rhob_local);
                double muB_local = eos.get_muB(e_local, rhob_local);
                double enthropy  = e_local + p_local;  // [1/fm^4]

                double Wtautau = 0.0;
                double Wtaux   = 0.0;
                double Wtauy   = 0.0;
                double Wtaueta = 0.0;
                double Wxx     = 0.0;
                double Wxy     = 0.0;
                double Wxeta   = 0.0;
                double Wyy     = 0.0;
                double Wyeta   = 0.0;
                double Wetaeta = 0.0;
                if (DATA.turn_on_shear == 1) {
                    Wtautau = arena(ix, iy, ieta).Wmunu[0]/enthropy;
                    Wtaux   = arena(ix, iy, ieta).Wmunu[1]/enthropy;
                    Wtauy   = arena(ix, iy, ieta).Wmunu[2]/enthropy;
                    Wtaueta = arena(ix, iy, ieta).Wmunu[3]/enthropy;
                    Wxx     = arena(ix, iy, ieta).Wmunu[4]/enthropy;
                    Wxy     = arena(ix, iy, ieta).Wmunu[5]/enthropy;
                    Wxeta   = arena(ix, iy, ieta).Wmunu[6]/enthropy;
                    Wyy     = arena(ix, iy, ieta).Wmunu[7]/enthropy;
                    Wyeta   = arena(ix, iy, ieta).Wmunu[8]/enthropy;
                    Wetaeta = arena(ix, iy, ieta).Wmunu[9]/enthropy;
                }

                double bulk_Pi = 0.0;
                if (DATA.turn_on_bulk == 1) {
                    bulk_Pi = arena(ix, iy, ieta).pi_b;  // [1/fm^4]
                }

                // outputs for baryon diffusion part
                double common_term_q = 0.0;
                double qtau          = 0.0;
                double qx            = 0.0;
                double qy            = 0.0;
                double qeta          = 0.0;
                if (DATA.turn_on_diff == 1) {
                    common_term_q = rhob_local*T_local/enthropy;
                    double kappa_hat = get_deltaf_qmu_coeff(T_local,
                                                            muB_local);
                    qtau = arena(ix, iy, ieta).Wmunu[10]/kappa_hat;
                    qx   = arena(ix, iy, ieta).Wmunu[11]/kappa_hat;
                    qy   = arena(ix, iy, ieta).Wmunu[12]/kappa_hat;
                    qeta = arena(ix, iy, ieta).Wmunu[13]/kappa_hat;
                }

                // exclude the actual coordinates from the output to save space:
                if (DATA.outputBinaryEvolution == 0) {
                    fprintf(out_file_xyeta, "%e %e %e %e %e\n",
                            T_local*hbarc, muB_local*hbarc, vx, vy, vz);
                    if (DATA.viscosity_flag == 1) {
                        if (DATA.turn_on_shear) {
                            fprintf(out_file_W_xyeta,
                                    "%e %e %e %e %e %e %e %e %e %e\n",
                                    Wtautau, Wtaux, Wtauy, Wtaueta, Wxx, Wxy,
                                    Wxeta, Wyy, Wyeta, Wetaeta);
                        }
                        if (DATA.turn_on_bulk) {
                            fprintf(out_file_bulkpi_xyeta,"%e %e %e\n", bulk_Pi, enthropy, cs2_local);
                        }
                    }
                } else {
                    float array[] = {static_cast<float>(T_local*hbarc),
                                     static_cast<float>(muB_local*hbarc),
                                     static_cast<float>(vx),
                                     static_cast<float>(vy),
                                     static_cast<float>(vz)};
                    fwrite(array, sizeof(float), 5, out_file_xyeta);
                    if (DATA.viscosity_flag == 1) {
                        if (DATA.turn_on_shear == 1) {
                            float array2[] = {static_cast<float>(Wtautau),
                                              static_cast<float>(Wtaux),
                                              static_cast<float>(Wtauy),
                                              static_cast<float>(Wtaueta),
                                              static_cast<float>(Wxx),
                                              static_cast<float>(Wxy),
                                              static_cast<float>(Wxeta),
                                              static_cast<float>(Wyy),
                                              static_cast<float>(Wyeta),
                                              static_cast<float>(Wetaeta)};
                            fwrite(array2, sizeof(float), 10,
                                   out_file_W_xyeta);
                        }
                        if (DATA.turn_on_bulk == 1) {
                            float array1[] = {static_cast<float>(bulk_Pi),
                                              static_cast<float>(enthropy),
                                              static_cast<float>(cs2_local)};
                            fwrite(array1, sizeof(float), 3,
                                   out_file_bulkpi_xyeta);
                        }
                        if (DATA.turn_on_diff == 1) {
                            float array3[] = {static_cast<float>(common_term_q),
                                              static_cast<float>(qtau),
                                              static_cast<float>(qx),
                                              static_cast<float>(qy),
                                              static_cast<float>(qeta)};
                            fwrite(array3, sizeof(float), 5,
                                   out_file_q_xyeta);
                        }
                    }
                }
            }
        }
    }
    fclose(out_file_xyeta);
    if (DATA.turn_on_shear == 1) {
        fclose(out_file_W_xyeta);
    }
    if (DATA.turn_on_bulk == 1) {
        fclose(out_file_bulkpi_xyeta);
    }
    if (DATA.turn_on_diff == 1) {
        fclose(out_file_q_xyeta);
    }
}


void Cell_info::OutputEvolution_Knudsen_Reynoldsnumbers(
        SCGrid &arena, const double tau) const {
    const string out_name_xyeta = "evolution_KRnumbers.dat";
    FILE *out_file_xyeta        = NULL;

    // If it's the first timestep, overwrite the previous file
    string out_open_mode;
    if (tau == DATA.tau0) {
        out_open_mode = "w";
    } else {
        out_open_mode = "a";
    }

    // If we output in binary, set the mode accordingly
    if (DATA.outputBinaryEvolution == 0) {
        out_open_mode += "b";
    }
    out_file_xyeta = fopen(out_name_xyeta.c_str(), out_open_mode.c_str());

    const int n_skip_x   = DATA.output_evolution_every_N_x;
    const int n_skip_y   = DATA.output_evolution_every_N_y;
    const int n_skip_eta = DATA.output_evolution_every_N_eta;
    for (int ieta = 0; ieta < arena.nEta(); ieta += n_skip_eta) {
        for (int iy = 0; iy < arena.nY(); iy += n_skip_y) {
            for (int ix = 0; ix < arena.nX(); ix += n_skip_x) {
                double R_pi = 0.0;
                double R_Pi = 0.0;
                calculate_inverse_Reynolds_numbers(arena, ieta, ix, iy,
                                                   R_pi, R_Pi);

                if (DATA.outputBinaryEvolution == 0) {
                    fprintf(out_file_xyeta, "%e %e\n", R_pi, R_Pi);
                } else {
                    float array[] = {static_cast<float>(R_pi),
                                     static_cast<float>(R_Pi)};
                    fwrite(array, sizeof(float), 2, out_file_xyeta);
                } 
            }
        }
    }
    fclose(out_file_xyeta);
}


void Cell_info::calculate_inverse_Reynolds_numbers(
                                SCGrid &arena_current,
                                const int ieta, const int ix, const int iy,
                                double &R_pi, double &R_Pi) const {
    const auto grid_pt = arena_current(ix, iy, ieta);

    const double e_local  = grid_pt.epsilon;
    const double rhob     = grid_pt.rhob;
    const double pressure = eos.get_pressure(e_local, rhob);

    const double pi_00 = grid_pt.Wmunu[0];
    const double pi_01 = grid_pt.Wmunu[1];
    const double pi_02 = grid_pt.Wmunu[2];
    const double pi_03 = grid_pt.Wmunu[3];
    const double pi_11 = grid_pt.Wmunu[4];
    const double pi_12 = grid_pt.Wmunu[5];
    const double pi_13 = grid_pt.Wmunu[6];
    const double pi_22 = grid_pt.Wmunu[7];
    const double pi_23 = grid_pt.Wmunu[8];
    const double pi_33 = grid_pt.Wmunu[9];

    const double pisize = (
           pi_00*pi_00 + pi_11*pi_11 + pi_22*pi_22 + pi_33*pi_33
         - 2.*(pi_01*pi_01 + pi_02*pi_02 + pi_03*pi_03)
         + 2.*(pi_12*pi_12 + pi_13*pi_13 + pi_23*pi_23));

    const double pi_local = grid_pt.pi_b;

    R_pi = sqrt(pisize)/pressure;
    R_Pi = pi_local/pressure;
}


//! This function outputs hydro evolution file into memory for JETSCAPE
void Cell_info::OutputEvolutionDataXYEta_memory(
            SCGrid &arena, const double tau, HydroinfoMUSIC &hydro_info_ptr) {
    const int n_skip_x   = DATA.output_evolution_every_N_x;
    const int n_skip_y   = DATA.output_evolution_every_N_y;
    const int n_skip_eta = DATA.output_evolution_every_N_eta;
    for (int ix = 0; ix < arena.nX(); ix += n_skip_x) {
        for (int iy = 0; iy < arena.nY(); iy += n_skip_y) {
            for (int ieta = 0; ieta < arena.nEta(); ieta += n_skip_eta) {
                double eta = 0.0;
                if (!DATA.boost_invariant) {
                    eta = ((static_cast<double>(ieta))*(DATA.delta_eta)
                            - (DATA.eta_size)/2.0);
                }
                double cosh_eta = cosh(eta);
                double sinh_eta = sinh(eta);

                double e_local    = arena(ix, iy, ieta).epsilon;  // 1/fm^4
                double rhob_local = arena(ix, iy, ieta).rhob;     // 1/fm^3
                double p_local = eos.get_pressure(e_local, rhob_local);
                double utau = arena(ix, iy, ieta).u[0];
                double ux   = arena(ix, iy, ieta).u[1];
                double uy   = arena(ix, iy, ieta).u[2];
                double ueta = arena(ix, iy, ieta).u[3];
                double ut = utau*cosh_eta + ueta*sinh_eta;  // gamma factor
                double vx = ux/ut;
                double vy = uy/ut;
                double uz = ueta*cosh_eta + utau*sinh_eta;
                double vz = uz/ut;

                double T_local   = eos.get_temperature(e_local, rhob_local);
                double s_local   = eos.get_entropy(e_local, rhob_local);

                hydro_info_ptr.dump_ideal_info_to_memory(
                    tau, eta, e_local, p_local, s_local, T_local, vx, vy, vz);
            }
        }
    }
}


//! This function outputs hydro evolution file in binary format
void Cell_info::OutputEvolutionDataXYEta_chun(SCGrid &arena, double tau) {
    // the format of the file is as follows,
    //    itau ix iy ieta e P T cs^2 ux uy ueta
    // if turn_on_shear == 1:
    //    itau ix iy ieta e P T cs^2 ux uy ueta Wxx Wxy Wxeta Wyy Wyeta
    // if turn_on_shear == 1 and turn_on_bulk == 1:
    //    itau ix iy ieta e P T cs^2 ux uy ueta Wxx Wxy Wxeta Wyy Wyeta pi_b
    // if turn_on_rhob == 1:
    //    itau ix iy ieta e P T cs^2 ux uy ueta rho_B mu_B
    // if turn_on_rhob == 1 and turn_on_shear == 1:
    //    itau ix iy ieta e P T cs^2 ux uy ueta rho_B mu_B Wxx Wxy Wxeta Wyy Wyeta
    // if turn_on_rhob == 1 and turn_on_shear == 1 and turn_on_diff == 1:
    //    itau ix iy ieta e P T cs^2 ux uy ueta rho_B mu_B Wxx Wxy Wxeta Wyy Wyeta qx qy qeta
    // if turn_on_rhob == 1 and turn_on_shear == 1 and turn_on_bulk == 1 and turn_on_diff == 1:
    //    itau ix iy ieta e P T cs^2 ux uy ueta rho_B mu_B Wxx Wxy Wxeta Wyy Wyeta pi_b qx qy qeta
    // Here ueta = tau*ueta, Wieta = tau*Wieta, Wetaeta = tau^2*Wetaeta, qeta = tau*qeta
    // Here Wij is reduced variables Wij/(e+P) in the fluid rest frame
    // and qi is reduced variables qi/kappa_hat in the fluid rest frame
    const string out_name_xyeta = "evolution_all_xyeta.dat";
    string out_open_mode;
    FILE *out_file_xyeta;
    // If it's the first timestep, overwrite the previous file
    if (tau == DATA.tau0) {
        out_open_mode = "wb";
    } else {
        out_open_mode = "ab";
    }
    out_file_xyeta = fopen(out_name_xyeta.c_str(), out_open_mode.c_str());

    int n_skip_tau     = DATA.output_evolution_every_N_timesteps;
    double output_dtau = DATA.delta_tau*n_skip_tau;
    int itau = static_cast<int>((tau - DATA.tau0)/(output_dtau) + 0.1);

    int n_skip_x       = DATA.output_evolution_every_N_x;
    int n_skip_y       = DATA.output_evolution_every_N_y;
    int n_skip_eta     = DATA.output_evolution_every_N_eta;

    // write out header
    const int output_nx        = static_cast<int>(arena.nX()/n_skip_x);
    const int output_ny        = static_cast<int>(arena.nY()/n_skip_y);
    const int output_neta      = static_cast<int>(arena.nEta()/n_skip_eta);
    const double output_dx     = DATA.delta_x*n_skip_x;
    const double output_dy     = DATA.delta_y*n_skip_y;
    const double output_deta   = DATA.delta_eta*n_skip_eta;
    const double output_xmin   = - DATA.x_size/2.;
    const double output_ymin   = - DATA.y_size/2.;
    const double output_etamin = - DATA.eta_size/2.;

    const int nVar_per_cell = (11 + DATA.turn_on_rhob*2 + DATA.turn_on_shear*5
                                  + DATA.turn_on_bulk*1 + DATA.turn_on_diff*3);
    if (tau == DATA.tau0) {
        float header[] = {
            static_cast<float>(DATA.tau0), static_cast<float>(output_dtau),
            static_cast<float>(output_nx), static_cast<float>(output_dx),
            static_cast<float>(output_xmin),
            static_cast<float>(output_ny), static_cast<float>(output_dy),
            static_cast<float>(output_ymin),
            static_cast<float>(output_neta), static_cast<float>(output_deta),
            static_cast<float>(output_etamin),
            static_cast<float>(DATA.turn_on_rhob),
            static_cast<float>(DATA.turn_on_shear),
            static_cast<float>(DATA.turn_on_bulk),
            static_cast<float>(DATA.turn_on_diff),
            static_cast<float>(nVar_per_cell)};
        fwrite(header, sizeof(float), 16, out_file_xyeta);
    }
    for (int ieta = 0; ieta < arena.nEta(); ieta += n_skip_eta) {
        double eta_local = - DATA.eta_size/2. + ieta*DATA.delta_eta;
        double cosh_eta = cosh(eta_local);
        double sinh_eta = sinh(eta_local);
        for (int iy = 0; iy < arena.nY(); iy += n_skip_y) {
            for (int ix = 0; ix < arena.nX(); ix += n_skip_x) {
                double e_local    = arena(ix, iy, ieta).epsilon;  // 1/fm^4
                double rhob_local = arena(ix, iy, ieta).rhob;     // 1/fm^3

                if (e_local*hbarc < DATA.output_evolution_e_cut) continue;
                // only ouput fluid cells that are above cut-off temperature

                double p_local    = eos.get_pressure(e_local, rhob_local);
                double cs2        = eos.get_cs2(e_local, rhob_local);

                double ux = arena(ix, iy, ieta).u[1];
                double uy = arena(ix, iy, ieta).u[2];
                double uz = (  arena(ix, iy, ieta).u[3]*cosh_eta
                             + arena(ix, iy, ieta).u[0]*sinh_eta);


                // T_local is in 1/fm
                double T_local = eos.get_temperature(e_local, rhob_local);

                double muB_local = 0.0;
                if (DATA.turn_on_rhob == 1)
                    muB_local = eos.get_muB(e_local, rhob_local);

                ShearVisVecLRF piLRF;
                get_LRF_shear_stress_tensor(arena(ix, iy, ieta), eta_local,
                                            piLRF);
                double div_factor = e_local + p_local;  // 1/fm^4
                double Wxx = 0.0;
                double Wxy = 0.0;
                double Wxz = 0.0;
                double Wyy = 0.0;
                double Wyz = 0.0;
                if (DATA.turn_on_shear == 1) {
                    Wxx = piLRF[0]/div_factor;
                    Wxy = piLRF[1]/div_factor;
                    Wxz = piLRF[2]/div_factor;
                    Wyy = piLRF[3]/div_factor;
                    Wyz = piLRF[4]/div_factor;
                }

                double pi_b = 0.0;
                if (DATA.turn_on_bulk == 1) {
                    pi_b = arena(ix, iy, ieta).pi_b/div_factor;
                }

                // outputs for baryon diffusion part
                //double common_term_q = 0.0;
                double qx = 0.0;
                double qy = 0.0;
                double qz = 0.0;
                if (DATA.turn_on_diff == 1) {
                    //common_term_q = rhob_local*T_local/div_factor;
                    double kappa_hat = get_deltaf_qmu_coeff(T_local,
                                                            muB_local);
                    qx = piLRF[5]/kappa_hat;
                    qy = piLRF[6]/kappa_hat;
                    qz = piLRF[7]/kappa_hat;
                }

                float ideal[] = {static_cast<float>(itau),
                                 static_cast<float>(ix/n_skip_x),
                                 static_cast<float>(iy/n_skip_y),
                                 static_cast<float>(ieta/n_skip_eta),
                                 static_cast<float>(e_local*hbarc),
                                 static_cast<float>(p_local*hbarc),
                                 static_cast<float>(T_local*hbarc),
                                 static_cast<float>(cs2),
                                 static_cast<float>(ux),
                                 static_cast<float>(uy),
                                 static_cast<float>(uz)};

                fwrite(ideal, sizeof(float), 11, out_file_xyeta);

                if (DATA.turn_on_rhob == 1) {
                    float mu[] = {static_cast<float>(rhob_local),
                                  static_cast<float>(muB_local*hbarc)};
                    fwrite(mu, sizeof(float), 2, out_file_xyeta);
                }

                if (DATA.turn_on_shear == 1) {
                    float shear_pi[] = {static_cast<float>(Wxx),
                                        static_cast<float>(Wxy),
                                        static_cast<float>(Wxz),
                                        static_cast<float>(Wyy),
                                        static_cast<float>(Wyz)};
                    fwrite(shear_pi, sizeof(float), 5, out_file_xyeta);
                }

                if (DATA.turn_on_bulk == 1) {
                    float bulk_pi[] = {static_cast<float>(pi_b)};
                    fwrite(bulk_pi, sizeof(float), 1, out_file_xyeta);
                }

                if (DATA.turn_on_diff == 1) {
                    float diffusion[] = {static_cast<float>(qx),
                                         static_cast<float>(qy),
                                         static_cast<float>(qz)};
                    fwrite(diffusion, sizeof(float), 3, out_file_xyeta);
                }
            }
        }
    }
    fclose(out_file_xyeta);
}


//! This function outputs hydro evolution file in binary format for photon production
void Cell_info::OutputEvolutionDataXYEta_photon(SCGrid &arena, double tau) {
    // volume = tau*dtau*dx*dy*deta
    // the format of the file is as follows,
    //    volume T ux uy ueta
    // if turn_on_shear == 1:
    //    volume T ux uy ueta Wxx Wxy Wxeta Wyy Wyeta
    // if turn_on_shear == 1 and turn_on_bulk == 1:
    //    volume T ux uy ueta Wxx Wxy Wxeta Wyy Wyeta pi_b
    // if turn_on_rhob == 1:
    //    volume T ux uy ueta mu_B
    // if turn_on_rhob == 1 and turn_on_shear == 1:
    //    volume T ux uy ueta mu_B Wxx Wxy Wxeta Wyy Wyeta
    // if turn_on_rhob == 1 and turn_on_shear == 1 and turn_on_diff == 1:
    //    volume T ux uy ueta mu_B Wxx Wxy Wxeta Wyy Wyeta qx qy qeta
    // Here ueta = tau*ueta, Wieta = tau*Wieta, qeta = tau*qeta
    // Here Wij is reduced variables Wij/(e+P) used in delta f
    // and qi is reduced variables qi/kappa_hat
    const string out_name_xyeta = "evolution_for_photon_xyeta.dat";
    string out_open_mode;
    FILE *out_file_xyeta;
    // If it's the first timestep, overwrite the previous file
    if (tau == DATA.tau0) {
        out_open_mode = "wb";
    } else {
        out_open_mode = "ab";
    }
    out_file_xyeta = fopen(out_name_xyeta.c_str(), out_open_mode.c_str());

    int n_skip_tau = DATA.output_evolution_every_N_timesteps;
    int n_skip_x = DATA.output_evolution_every_N_x;
    int n_skip_y = DATA.output_evolution_every_N_y;
    int n_skip_eta = DATA.output_evolution_every_N_eta;
    double dtau = DATA.delta_tau;
    double dx = DATA.delta_x;
    double dy = DATA.delta_y;
    double deta = DATA.delta_eta;
    double volume = tau*n_skip_tau*dtau*n_skip_x*dx*n_skip_y*dy*n_skip_eta*deta;

    for (int ieta = 0; ieta < arena.nEta(); ieta += n_skip_eta) {
        double eta_local = - DATA.eta_size/2. + ieta*deta;
        for (int iy = 0; iy < arena.nY(); iy += n_skip_y) {
            for (int ix = 0; ix < arena.nX(); ix += n_skip_x) {
                double e_local = arena(ix, iy, ieta).epsilon;  // 1/fm^4

                if (e_local*hbarc < DATA.output_evolution_e_cut) continue;
                // only ouput fluid cells that are above cut-off temperature

                double rhob_local = arena(ix, iy, ieta).rhob;  // 1/fm^3

                double ux   = arena(ix, iy, ieta).u[1];
                double uy   = arena(ix, iy, ieta).u[2];
                double ueta = arena(ix, iy, ieta).u[3];

                // T_local is in 1/fm
                double T_local = eos.get_temperature(e_local, rhob_local);
                double muB_local = 0.0;
                if (DATA.turn_on_rhob == 1) {
                    muB_local = eos.get_muB(e_local, rhob_local);
                }

                double p_local = eos.get_pressure(e_local, rhob_local);
                double div_factor = e_local + p_local;  // 1/fm^4
                double Wxx = 0.0;
                double Wxy = 0.0;
                double Wxeta = 0.0;
                double Wyy = 0.0;
                double Wyeta = 0.0;
                if (DATA.turn_on_shear == 1) {
                    Wxx   = arena(ix, iy, ieta).Wmunu[4]/div_factor;
                    Wxy   = arena(ix, iy, ieta).Wmunu[5]/div_factor;
                    Wxeta = arena(ix, iy, ieta).Wmunu[6]/div_factor;
                    Wyy   = arena(ix, iy, ieta).Wmunu[7]/div_factor;
                    Wyeta = arena(ix, iy, ieta).Wmunu[8]/div_factor;
                }

                double pi_b = 0.0;
                if (DATA.turn_on_bulk == 1) {
                    pi_b = arena(ix, iy, ieta).pi_b;   // 1/fm^4
                }

                // outputs for baryon diffusion part
                double common_term_q = 0.0;
                double qx = 0.0;
                double qy = 0.0;
                double qeta = 0.0;
                if (DATA.turn_on_diff == 1) {
                    common_term_q = rhob_local*T_local/div_factor;
                    double kappa_hat = get_deltaf_qmu_coeff(T_local,
                                                            muB_local);
                    qx   = arena(ix, iy, ieta).Wmunu[11]/kappa_hat;
                    qy   = arena(ix, iy, ieta).Wmunu[12]/kappa_hat;
                    qeta = arena(ix, iy, ieta).Wmunu[13]/kappa_hat;
                }

                float ideal[] = {static_cast<float>(volume),
                                 static_cast<float>(eta_local),
                                 static_cast<float>(T_local*hbarc),
                                 static_cast<float>(ux),
                                 static_cast<float>(uy),
                                 static_cast<float>(ueta)};

                fwrite(ideal, sizeof(float), 6, out_file_xyeta);

                if (DATA.turn_on_rhob == 1) {
                    float mu[] = {static_cast<float>(muB_local*hbarc)};
                    fwrite(mu, sizeof(float), 1, out_file_xyeta);
                }

                if (DATA.turn_on_shear == 1) {
                    float shear_pi[] = {static_cast<float>(Wxx),
                                        static_cast<float>(Wxy),
                                        static_cast<float>(Wxeta),
                                        static_cast<float>(Wyy),
                                        static_cast<float>(Wyeta)};
                    fwrite(shear_pi, sizeof(float), 5, out_file_xyeta);
                }

                if (DATA.turn_on_bulk == 1) {
                    float bulk_pi[] = {static_cast<float>(pi_b)};
                    fwrite(bulk_pi, sizeof(float), 1, out_file_xyeta);
                }

                if (DATA.turn_on_diff == 1) {
                    float diffusion[] = {static_cast<float>(common_term_q),
                                         static_cast<float>(qx),
                                         static_cast<float>(qy),
                                         static_cast<float>(qeta)};
                    fwrite(diffusion, sizeof(float), 4, out_file_xyeta);
                }
            }
        }
    }
    fclose(out_file_xyeta);
}


//! This function outputs hydro evolution file in binary format
void Cell_info::OutputEvolutionDataXYEta_vorticity(
        SCGrid &arena_curr, SCGrid &arena_prev, double tau) {
    // the format of the file is as follows,
    //    itau ix iy ieta e P T ux uy ueta mu_B
    //    omega^tx omega^ty omega^tz omega^xy omega^xz omega^yz
    const string out_name_xyeta = "evolution_all_xyeta.dat";
    string out_open_mode;
    FILE *out_file_xyeta;
    // If it's the first timestep, overwrite the previous file
    if (tau == DATA.tau0) {
        out_open_mode = "wb";
    } else {
        out_open_mode = "ab";
    }
    out_file_xyeta = fopen(out_name_xyeta.c_str(), out_open_mode.c_str());

    int n_skip_tau     = DATA.output_evolution_every_N_timesteps;
    double output_dtau = DATA.delta_tau*n_skip_tau;
    int itau = static_cast<int>((tau - DATA.tau0)/(output_dtau) + 0.1);

    int n_skip_x       = DATA.output_evolution_every_N_x;
    int n_skip_y       = DATA.output_evolution_every_N_y;
    int n_skip_eta     = DATA.output_evolution_every_N_eta;

    // write out header
    const int output_nx        = static_cast<int>(arena_curr.nX()/n_skip_x);
    const int output_ny        = static_cast<int>(arena_curr.nY()/n_skip_y);
    const int output_neta      = static_cast<int>(arena_curr.nEta()/n_skip_eta);
    const double output_dx     = DATA.delta_x*n_skip_x;
    const double output_dy     = DATA.delta_y*n_skip_y;
    const double output_deta   = DATA.delta_eta*n_skip_eta;
    const double output_xmin   = - DATA.x_size/2.;
    const double output_ymin   = - DATA.y_size/2.;
    const double output_etamin = - DATA.eta_size/2.;

    const int nVar_per_cell = 35;

    if (tau == DATA.tau0) {
        float header[] = {
            static_cast<float>(DATA.tau0), static_cast<float>(output_dtau),
            static_cast<float>(output_nx), static_cast<float>(output_dx),
            static_cast<float>(output_xmin),
            static_cast<float>(output_ny), static_cast<float>(output_dy),
            static_cast<float>(output_ymin),
            static_cast<float>(output_neta), static_cast<float>(output_deta),
            static_cast<float>(output_etamin),
            static_cast<float>(nVar_per_cell)};
        fwrite(header, sizeof(float), 12, out_file_xyeta);
    }
    for (int ieta = 0; ieta < arena_curr.nEta(); ieta += n_skip_eta) {
        double eta_local = - DATA.eta_size/2. + ieta*DATA.delta_eta;
        for (int iy = 0; iy < arena_curr.nY(); iy += n_skip_y) {
            for (int ix = 0; ix < arena_curr.nX(); ix += n_skip_x) {
                double e_local    = arena_curr(ix, iy, ieta).epsilon;  // 1/fm^4
                double rhob_local = arena_curr(ix, iy, ieta).rhob;     // 1/fm^3
                double p_local    = eos.get_pressure(e_local, rhob_local);

                double ux   = arena_curr(ix, iy, ieta).u[1];
                double uy   = arena_curr(ix, iy, ieta).u[2];
                double ueta = arena_curr(ix, iy, ieta).u[3];

                // T_local is in GeV
                double T_local = eos.get_temperature(e_local, rhob_local)*hbarc;

                if (T_local < DATA.output_evolution_T_cut) continue;
                // only ouput fluid cells that are above cut-off temperature

                double muB_local = eos.get_muB(e_local, rhob_local);

                VorticityVec omega_kSP = {0.0};
                VorticityVec omega_k   = {0.0};
                VorticityVec omega_th  = {0.0};
                VorticityVec omega_T   = {0.0};
                u_derivative_helper.compute_vorticity_shell(
                    tau, arena_prev, arena_curr, ieta, ix, iy, eta_local,
                    omega_kSP, omega_k, omega_th, omega_T);

                float ideal[] = {static_cast<float>(itau),
                                 static_cast<float>(ix/n_skip_x),
                                 static_cast<float>(iy/n_skip_y),
                                 static_cast<float>(ieta/n_skip_eta),
                                 static_cast<float>(e_local*hbarc),
                                 static_cast<float>(p_local*hbarc),
                                 static_cast<float>(T_local),
                                 static_cast<float>(ux),
                                 static_cast<float>(uy),
                                 static_cast<float>(ueta),
                                 static_cast<float>(muB_local*hbarc)};

                fwrite(ideal, sizeof(float), 11, out_file_xyeta);

                float vor_vec[6];
                for (int i = 0; i < 6; i++) {
                    // no minus sign because it has an opposite sign compared
                    // to kinetic vorticity
                    vor_vec[i] = static_cast<float>(omega_kSP[i]*hbarc);  // GeV
                }
                fwrite(vor_vec, sizeof(float), 6, out_file_xyeta);
                for (int i = 0; i < 6; i++) {
                    // the minus sign is from metric
                    // output quantities for g = (1, -1, -1 , -1)
                    vor_vec[i] = static_cast<float>(-omega_k[i]*hbarc);  // GeV
                }
                fwrite(vor_vec, sizeof(float), 6, out_file_xyeta);
                for (int i = 0; i < 6; i++) {
                    // the minus sign is from metric
                    // output quantities for g = (1, -1, -1 , -1)
                    vor_vec[i] = static_cast<float>(-omega_th[i]);  // 1
                }
                fwrite(vor_vec, sizeof(float), 6, out_file_xyeta);
                for (int i = 0; i < 6; i++) {
                    // the minus sign is from metric
                    // output quantities for g = (1, -1, -1 , -1)
                    vor_vec[i] = static_cast<float>(-omega_T[i]*hbarc*hbarc);  // GeV^2
                }
                fwrite(vor_vec, sizeof(float), 6, out_file_xyeta);
            }
        }
    }
    fclose(out_file_xyeta);
}


//! This function prints to the screen the maximum local energy density,
//! the maximum temperature in the current grid
void Cell_info::get_maximum_energy_density(
        SCGrid &arena, double &e_max, double &nB_max, double &Tmax) {
    double eps_max  = 0.0;
    double rhob_max = 0.0;
    double T_max    = 0.0;

    // get the grid information
    const int neta = arena.nEta();
    const int nx   = arena.nX();
    const int ny   = arena.nY();

    #pragma omp parallel for collapse(3) reduction(max:eps_max, rhob_max, T_max)
    for (int ieta = 0; ieta < neta; ieta++)
    for (int ix = 0; ix < nx; ix++) 
    for (int iy = 0; iy < ny; iy++) {
        const auto eps_local  = arena(ix, iy, ieta).epsilon;
        const auto rhob_local = arena(ix, iy, ieta).rhob;
        eps_max  = std::max(eps_max,  eps_local );
        rhob_max = std::max(rhob_max, rhob_local);
        T_max    = std::max(T_max,    eos.get_temperature(eps_local, rhob_local));
    }
    eps_max *= Util::hbarc;   // GeV/fm^3
    T_max   *= Util::hbarc;   // GeV

    if (eps_max > 1e5) {
        music_message << "The maximum e = " << eps_max << " < 1e5 GeV/fm^3";
        music_message.flush("error");
        music_message.error("This normally should not happen!");
        music_message.error("Exiting ...");
        exit(1);
    }
    music_message << "eps_max = " << eps_max << " GeV/fm^3, "
                  << "rhob_max = " << rhob_max << " 1/fm^3, "
                  << "T_max = " << T_max << " GeV.";
    music_message.flush("info");
    e_max = eps_max;
    nB_max = rhob_max;
    Tmax = T_max;
}


//! This function computes global angular momentum at a give proper time
void Cell_info::compute_angular_momentum(
        SCGrid &arena, SCGrid &arena_prev, const double tau,
        const double eta_min, const double eta_max) {
    ostringstream filename;
    filename << "global_angular_momentum_eta_"
             << eta_min << "_" << eta_max << ".dat";
    ofstream output_file;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        // create new files at the first time step
        output_file.open(filename.str().c_str(), std::ofstream::out);
        output_file << "# tau[fm]  Lx[hbarc]  Ly[hbarc]  Lz[hbarc]  "
                    << "L^{tx}[hbarc]  L^{ty}[hbarc]  L^{tz}[hbarc]"
                    << std::endl;
    } else {
        output_file.open(filename.str().c_str(),
                         std::fstream::out | std::fstream::app);
    }
    double Lx  = 0.0;
    double Ly  = 0.0;
    double Lz  = 0.0;
    double Ltx = 0.0;
    double Lty = 0.0;
    double Ltz = 0.0;
    const double deta = DATA.delta_eta;
    const double dx   = DATA.delta_x;
    const double dy   = DATA.delta_y;
    const int neta    = arena.nEta();
    const int nx      = arena.nX();
    const int ny      = arena.nY();
    #pragma omp parallel for collapse(3) reduction(+:Lx, Ly, Lz, Ltx, Lty, Ltz)
    for (int ieta = 0; ieta < neta; ieta++)
    for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++) {
        const auto& c      = arena     (ix, iy, ieta);
        const auto& c_prev = arena_prev(ix, iy, ieta);

        double eta_s = deta*ieta - (DATA.eta_size)/2.0;
        if (DATA.boost_invariant) {
            eta_s = 0.0;
        }

        const double cosh_eta = cosh(eta_s);
        const double sinh_eta = sinh(eta_s);
        const double t_local = tau*cosh_eta;
        const double x_local = DATA.x_size/2. + ix*dx;
        const double y_local = DATA.x_size/2. + iy*dy;
        const double z_local = tau*sinh_eta;

        const double e_local   = c.epsilon;
        const double rhob      = c.rhob;
        const double pressure  = eos.get_pressure(e_local, rhob);
        const double u0        = c.u[0];
        const double u1        = c.u[1];
        const double u2        = c.u[2];
        const double u3        = c.u[3];

        const double T00_local = (e_local + pressure)*u0*u0 - pressure;
        const double Pi00_rk_0 = (c_prev.pi_b
                                  *(-1.0 + c_prev.u[0]*c_prev.u[0]));

        const double T_tau_tau = (T00_local + c_prev.Wmunu[0] + Pi00_rk_0);
        const double T_tau_x   = ((e_local + pressure)*u0*u1 + c_prev.Wmunu[1]
                                  + c_prev.pi_b*c_prev.u[0]*c_prev.u[1]);
        const double T_tau_y   = ((e_local + pressure)*u0*u2 + c_prev.Wmunu[2]
                                  + c_prev.pi_b*c_prev.u[0]*c_prev.u[2]);
        const double T_tau_eta = ((e_local + pressure)*u0*u3 + c_prev.Wmunu[3]
                                  + c_prev.pi_b*c_prev.u[0]*c_prev.u[3]);
        const double T_tau_t = T_tau_tau*cosh_eta + T_tau_eta*sinh_eta;
        const double T_tau_z = T_tau_tau*sinh_eta + T_tau_eta*cosh_eta;

        if (eta_s < eta_max && eta_s > eta_min) {
            Lx  += (y_local*T_tau_z - z_local*T_tau_y);
            Ly  += (z_local*T_tau_x - x_local*T_tau_z);
            Lz  += (x_local*T_tau_y - y_local*T_tau_x);
            Ltx += (t_local*T_tau_x - x_local*T_tau_t);
            Lty += (t_local*T_tau_y - y_local*T_tau_t);
            Ltz += (t_local*T_tau_z - z_local*T_tau_t);
        }
    }
    // add units
    double factor = tau*dx*dy*deta;
    Lx  *= factor;
    Ly  *= factor;
    Lz  *= factor;
    Ltx *= factor;
    Lty *= factor;
    Ltz *= factor;

    // output results
    output_file << scientific << setprecision(6)
                << tau << "  " << Lx << "  " << Ly << "  " << Lz << "  "
                << Ltx << "  " << Lty << "  " << Ltz << std::endl;
    output_file.close();
}


//! This function checks the total energy and total net baryon number
//! at a give proper time
void Cell_info::check_conservation_law(SCGrid &arena, SCGrid &arena_prev,
                                       const double tau) {
    std::string filename = "global_conservation_laws.dat";
    ofstream output_file;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        output_file.open(filename.c_str(), std::ofstream::out);
        output_file << "# tau(fm)  E(GeV)  Px(GeV)  Py(GeV)  Pz(GeV)  N_B "
                    << std::endl;
    } else {
        output_file.open(filename.c_str(),
                         std::fstream::out | std::fstream::app);
    }
    double N_B     = 0.0;
    double T_tau_t = 0.0;
    double T_tau_x = 0.0;
    double T_tau_y = 0.0;
    double T_tau_z = 0.0;
    double N_B_edge     = 0.0;
    double T_tau_t_edge = 0.0;
    double T_tau_x_edge = 0.0;
    double T_tau_y_edge = 0.0;
    double T_tau_z_edge = 0.0;
    double deta    = DATA.delta_eta;
    double dx      = DATA.delta_x;
    double dy      = DATA.delta_y;
    const int neta = arena.nEta();
    const int nx   = arena.nX();
    const int ny   = arena.nY();

    #pragma omp parallel for collapse(3) reduction(+:N_B, T_tau_t, T_tau_x, T_tau_y, T_tau_z, N_B_edge, T_tau_t_edge, T_tau_x_edge, T_tau_y_edge, T_tau_z_edge)
    for (int ieta = 0; ieta < neta; ieta++)
    for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++) {
        const auto& c      = arena     (ix, iy, ieta);
        const auto& c_prev = arena_prev(ix, iy, ieta);

        const double eta_s = deta*ieta - (DATA.eta_size)/2.0;
        const double cosh_eta = cosh(eta_s);
        const double sinh_eta = sinh(eta_s);
        N_B += (c.rhob*c.u[0] + c_prev.Wmunu[10]);
        const double e_local   = c.epsilon;
        const double rhob      = c.rhob;
        const double pressure  = eos.get_pressure(e_local, rhob);
        const double u0        = c.u[0];
        const double u1        = c.u[1];
        const double u2        = c.u[2];
        const double u3        = c.u[3];
        const double T00_local = (e_local + pressure)*u0*u0 - pressure;
        const double Pi00_rk_0 = (c_prev.pi_b
                                  *(-1.0 + c_prev.u[0]*c_prev.u[0]));

        const double T_tau_tau = (T00_local + c_prev.Wmunu[0] + Pi00_rk_0);
        const double T01_local = ((e_local + pressure)*u0*u1 + c_prev.Wmunu[1]
                                  + c_prev.pi_b*c_prev.u[0]*c_prev.u[1]);
        const double T02_local = ((e_local + pressure)*u0*u2 + c_prev.Wmunu[2]
                                  + c_prev.pi_b*c_prev.u[0]*c_prev.u[2]);
        const double T_tau_eta = ((e_local + pressure)*u0*u3 + c_prev.Wmunu[3]
                                  + c_prev.pi_b*c_prev.u[0]*c_prev.u[3]);
        T_tau_t += T_tau_tau*cosh_eta + T_tau_eta*sinh_eta;
        T_tau_x += T01_local;
        T_tau_y += T02_local;
        T_tau_z += T_tau_tau*sinh_eta + T_tau_eta*cosh_eta;

        // compute the energy-momentum vector on the edge
        if (ieta == 0 || ieta == neta - 1 || ix == 0 || ix == nx - 1
            || iy == 0 || iy == ny - 1) {
            N_B_edge     += c.rhob*c.u[0] + c_prev.Wmunu[10];
            T_tau_t_edge += T_tau_tau*cosh_eta + T_tau_eta*sinh_eta;
            T_tau_x_edge += T01_local;
            T_tau_y_edge += T02_local;
            T_tau_z_edge += T_tau_tau*sinh_eta + T_tau_eta*cosh_eta;
        }
    }
    // add units
    double factor = tau*dx*dy*deta;
    N_B *= factor;
    T_tau_t *= factor*Util::hbarc;  // GeV
    T_tau_x *= factor*Util::hbarc;  // GeV
    T_tau_y *= factor*Util::hbarc;  // GeV
    T_tau_z *= factor*Util::hbarc;  // GeV
    N_B_edge *= factor;
    T_tau_t_edge *= factor*Util::hbarc;  // GeV
    T_tau_x_edge *= factor*Util::hbarc;  // GeV
    T_tau_y_edge *= factor*Util::hbarc;  // GeV
    T_tau_z_edge *= factor*Util::hbarc;  // GeV

    // compute the outflow flux
    if (tau > DATA.tau0) {
        outflow_flux[0] += T_tau_t_edge - Pmu_edge_prev[0];
        outflow_flux[1] += T_tau_x_edge - Pmu_edge_prev[1];
        outflow_flux[2] += T_tau_y_edge - Pmu_edge_prev[2];
        outflow_flux[3] += T_tau_z_edge - Pmu_edge_prev[3];
        outflow_flux[4] += N_B_edge - Pmu_edge_prev[4];

        N_B     += outflow_flux[4];
        T_tau_t += outflow_flux[0];  // GeV
        T_tau_x += outflow_flux[1];  // GeV
        T_tau_y += outflow_flux[2];  // GeV
        T_tau_z += outflow_flux[3];  // GeV
    }
    Pmu_edge_prev[0] = T_tau_t_edge;
    Pmu_edge_prev[1] = T_tau_x_edge;
    Pmu_edge_prev[2] = T_tau_y_edge;
    Pmu_edge_prev[3] = T_tau_z_edge;
    Pmu_edge_prev[4] = N_B_edge;

    // output results
    music_message << "total energy T^{taut} = " << T_tau_t << " GeV";
    music_message.flush("info");
    music_message << "net longitudinal momentum Pz = " << T_tau_z << " GeV";
    music_message.flush("info");
    music_message << "net baryon number N_B = " << N_B;
    if (N_B > 0. || N_B < 500.) {
        music_message.flush("info");
    } else {
        music_message.flush("error");
        exit(1);
    }
    output_file << scientific << setprecision(6)
                << tau << "  " << T_tau_t << "  " << T_tau_x << "  "
                << T_tau_y << "  " << T_tau_z << "  " << N_B << std::endl;
    output_file.close();
}


//! This function putputs files to check with Gubser flow solution
void Cell_info::Gubser_flow_check_file(SCGrid &arena, const double tau) {
    if (tau > 1.) {
        ostringstream filename_analytic;
        filename_analytic << "tests/Gubser_flow/y=x_tau="
                          << tau << "_SemiAnalytic.dat";

        double T_analytic[201], ux_analytic[201], uy_analytic[201];
        double pixx_analytic[201], pixy_analytic[201];
        double piyy_analytic[201], pizz_analytic[201];
        double dummy;
        std::ifstream input_file(filename_analytic.str().c_str());
        for (int i = 0; i < 201; i++) {
            input_file >> dummy >> dummy >> T_analytic[i] >> ux_analytic[i]
                       >> uy_analytic[i] >> pixx_analytic[i]
                       >> piyy_analytic[i] >> pixy_analytic[i]
                       >> pizz_analytic[i];
        }
        input_file.close();

        double T_diff    = 0.0;
        double ux_diff   = 0.0;
        double uy_diff   = 0.0;
        double pixx_diff = 0.0;
        double pixy_diff = 0.0;
        double piyy_diff = 0.0;
        double pizz_diff = 0.0;
        double T_sum     = 0.0;
        double ux_sum    = 0.0;
        double uy_sum    = 0.0;
        double pixx_sum  = 0.0;
        double pixy_sum  = 0.0;
        double piyy_sum  = 0.0;
        double pizz_sum  = 0.0;
        for (int i = 0; i < arena.nX(); i++) {
            double e_local = arena(i,i,0).epsilon;
            double T_local = (
                    eos.get_temperature(e_local, 0.0)*Util::hbarc);
            T_diff += fabs(T_analytic[i] - T_local);
            T_sum += fabs(T_analytic[i]);
            ux_diff += fabs(ux_analytic[i] - arena(i,i,0).u[1]);
            ux_sum += fabs(ux_analytic[i]);
            uy_diff += fabs(uy_analytic[i] - arena(i,i,0).u[2]);
            uy_sum += fabs(uy_analytic[i]);
            pixx_diff += (fabs(pixx_analytic[i]
                               - arena(i,i,0).Wmunu[4]*Util::hbarc));
            pixx_sum += fabs(pixx_analytic[i]);
            pixy_diff += (fabs(pixx_analytic[i]
                               - arena(i,i,0).Wmunu[5]*Util::hbarc));
            pixy_sum += fabs(pixx_analytic[i]);
            piyy_diff += (fabs(piyy_analytic[i]
                               - arena(i,i,0).Wmunu[7]*Util::hbarc));
            piyy_sum += fabs(piyy_analytic[i]);
            pizz_diff += (fabs(pizz_analytic[i]
                               - arena(i,i,0).Wmunu[9]*Util::hbarc));
            pizz_sum += fabs(pizz_analytic[i]);
        }
        music_message << "Autocheck: T_diff = " << T_diff/T_sum
                      << ", ux_diff = " << ux_diff/ux_sum
                      << ", uy_diff = " << uy_diff/uy_sum
                      << ", pixx_diff = " << pixx_diff/pixx_sum
                      << ", pixy_diff = " << pixy_diff/pixy_sum
                      << ", piyy_diff = " << piyy_diff/piyy_sum
                      << ", pizz_diff = " << pizz_diff/pizz_sum;
        music_message.flush("info");
    }

    ostringstream filename;
    filename << "Gubser_flow_check_tau_" << tau << ".dat";
    ofstream output_file(filename.str().c_str());

    double dx = DATA.delta_x;
    double x_min = -DATA.x_size/2.;
    double dy = DATA.delta_y;
    double y_min = -DATA.y_size/2.;
    for (int ix = 0; ix < arena.nX(); ix++)
    for (int iy = 0; iy < arena.nY(); iy++) {
        double x_local = x_min + ix*dx;
        double y_local = y_min + iy*dy;
        double e_local = arena(ix,iy,0).epsilon;
        double rhob_local = arena(ix,iy,0).rhob;
        double T_local = eos.get_temperature(e_local, 0.0);
        output_file << scientific << setprecision(8) << setw(18)
                    << x_local << "  " << y_local << "  "
                    << e_local*Util::hbarc << "  " << rhob_local << "  "
                    << T_local*Util::hbarc << "  "
                    << arena(ix,iy,0).u[1] << "  "
                    << arena(ix,iy,0).u[2] << "  "
                    << arena(ix,iy,0).Wmunu[4]*Util::hbarc << "  "
                    << arena(ix,iy,0).Wmunu[7]*Util::hbarc << "  "
                    << arena(ix,iy,0).Wmunu[5]*Util::hbarc << "  "
                    << arena(ix,iy,0).Wmunu[9]*Util::hbarc << "  "
                    << endl;
    }
    output_file.close();
}


//! This function outputs files to cross check with 1+1D simulation
void Cell_info::output_1p1D_check_file(SCGrid &arena, const double tau) {
    ostringstream filename;
    filename << "1+1D_check_tau_" << tau << ".dat";
    ofstream output_file(filename.str().c_str());

    double deta = DATA.delta_eta;
    double eta_min = -6.94;
    for (int ieta = 0; ieta < arena.nEta(); ieta++) {
        double eta_local = eta_min + ieta*deta;
        double e_local = arena(1, 1, ieta).epsilon;
        double rhob_local = arena(1, 1, ieta).rhob;
        output_file << scientific << setprecision(8) << setw(18)
                    << eta_local << "  "
                    << e_local*Util::hbarc << "  " << rhob_local
                    << endl;
    }
    output_file.close();
}


//! This function outputs energy density and n_b for making movies
void Cell_info::output_evolution_for_movie(SCGrid &arena, const double tau) {
    const string out_name_xyeta = "evolution_for_movie_xyeta.dat";
    string out_open_mode;
    FILE *out_file_xyeta;
    // If it's the first timestep, overwrite the previous file
    if (tau == DATA.tau0) {
        out_open_mode = "wb";
    } else {
        out_open_mode = "ab";
    }
    out_file_xyeta = fopen(out_name_xyeta.c_str(), out_open_mode.c_str());

    int n_skip_tau = DATA.output_evolution_every_N_timesteps;
    double output_dtau = DATA.delta_tau*n_skip_tau;
    int itau = static_cast<int>((tau - DATA.tau0)/(output_dtau) + 0.1);

    int n_skip_x   = DATA.output_evolution_every_N_x;
    int n_skip_y   = DATA.output_evolution_every_N_y;
    int n_skip_eta = DATA.output_evolution_every_N_eta;
    double dx      = DATA.delta_x;
    double dy      = DATA.delta_y;
    double deta    = DATA.delta_eta;
    double volume  = tau*output_dtau*n_skip_x*dx*n_skip_y*dy*n_skip_eta*deta;

    const int nVar_per_cell = 14;
    if (tau == DATA.tau0) {
        // write out header
        const int output_nx   = static_cast<int>(arena.nX()/n_skip_x);
        const int output_ny   = static_cast<int>(arena.nY()/n_skip_y);
        const int output_neta = (
                std::min(1, static_cast<int>(arena.nEta()/n_skip_eta)));
        const double output_dx     = DATA.delta_x*n_skip_x;
        const double output_dy     = DATA.delta_y*n_skip_y;
        const double output_deta   = DATA.delta_eta*n_skip_eta;
        const double output_xmin   = - DATA.x_size/2.;
        const double output_ymin   = - DATA.y_size/2.;
        const double output_etamin = - DATA.eta_size/2.;
        float header[] = {
            static_cast<float>(DATA.tau0), static_cast<float>(output_dtau),
            static_cast<float>(output_nx), static_cast<float>(output_dx),
            static_cast<float>(output_xmin),
            static_cast<float>(output_ny), static_cast<float>(output_dy),
            static_cast<float>(output_ymin),
            static_cast<float>(output_neta), static_cast<float>(output_deta),
            static_cast<float>(output_etamin),
            static_cast<float>(nVar_per_cell)};
        fwrite(header, sizeof(float), 12, out_file_xyeta);
    }
    for (int ieta = 0; ieta < arena.nEta(); ieta += n_skip_eta) {
        double eta_local = - DATA.eta_size/2. + ieta*deta;
        if (DATA.boost_invariant) eta_local = 0.;
        for (int iy = 0; iy < arena.nY(); iy += n_skip_y) {
            for (int ix = 0; ix < arena.nX(); ix += n_skip_x) {
                double e_local = arena(ix, iy, ieta).epsilon;  // 1/fm^4
                double rhob_local = arena(ix, iy, ieta).rhob;  // 1/fm^3

                // T_local is in 1/fm
                double T_local   = eos.get_temperature(e_local, rhob_local);
                if (T_local*hbarc < DATA.output_evolution_T_cut) continue;

                double muB_local = eos.get_muB(e_local, rhob_local);  // 1/fm

                double pressure  = eos.get_pressure(e_local, rhob_local);
                double u0        = arena(ix, iy, ieta).u[0];
                double u1        = arena(ix, iy, ieta).u[1];
                double u2        = arena(ix, iy, ieta).u[2];
                double u3        = arena(ix, iy, ieta).u[3];
                double T00_ideal = (e_local + pressure)*u0*u0 - pressure;
                double T03_ideal = (e_local + pressure)*u0*u3;
                double Pi00      = arena(ix, iy, ieta).pi_b*(-1.0 + u0*u0);
                double Pi03      = arena(ix, iy, ieta).pi_b*u0*u3;
                double T00_full  = (  T00_ideal
                                    + arena(ix, iy, ieta).Wmunu[0] + Pi00);
                double T03_full  = (  T03_ideal
                                    + arena(ix, iy, ieta).Wmunu[3] + Pi03);
                double Ttaut     = (  T00_full*cosh(eta_local)
                                    + T03_full*sinh(eta_local));
                double JBtau     = (  rhob_local*u0
                                    + arena(ix, iy, ieta).Wmunu[10]);
                float array[] = {static_cast<float>(itau),
                                 static_cast<float>(ix/n_skip_x),
                                 static_cast<float>(iy/n_skip_y),
                                 static_cast<float>(ieta/n_skip_eta),
                                 static_cast<float>(volume),
                                 static_cast<float>(e_local*hbarc),
                                 static_cast<float>(rhob_local),
                                 static_cast<float>(T_local*hbarc),
                                 static_cast<float>(muB_local*hbarc),
                                 static_cast<float>(u1),
                                 static_cast<float>(u2),
                                 static_cast<float>(u3),
                                 static_cast<float>(Ttaut*hbarc),
                                 static_cast<float>(JBtau)};
                fwrite(array, sizeof(float), 14, out_file_xyeta);
            }
        }
    }
    fclose(out_file_xyeta);
}


//! This function dumps the energy density and net baryon density
void Cell_info::output_energy_density_and_rhob_disitrubtion(SCGrid &arena,
                                                            string filename) {
    ofstream output_file(filename.c_str());
    const int n_skip_x   = DATA.output_evolution_every_N_x;
    const int n_skip_y   = DATA.output_evolution_every_N_y;
    const int n_skip_eta = DATA.output_evolution_every_N_eta;
    for (int ieta = 0; ieta < arena.nEta(); ieta += n_skip_eta)
    for (int ix   = 0; ix   < arena.nX();   ix += n_skip_x)
    for (int iy   = 0; iy   < arena.nY();   iy += n_skip_y) {
        double e_local = arena(ix, iy, ieta).epsilon*Util::hbarc;
        double rhob_local = arena(ix, iy, ieta).rhob;
        output_file << scientific << setprecision(5) << setw(18)
                    << e_local << "  " << rhob_local << endl;
    }
    output_file.close();
}


//! This function outputs the evolution of hydrodynamic variables at a
//! give fluid cell
void Cell_info::monitor_a_fluid_cell(SCGrid &arena_curr, SCGrid &arena_prev,
                                     const int ix, const int iy,
                                     const int ieta, const double tau) {
    ostringstream filename;
    filename << "monitor_fluid_cell_ix_" << ix << "_iy_" << iy
             << "_ieta_" << ieta << ".dat";
    ofstream output_file;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        output_file.open(filename.str().c_str(), std::ofstream::out);
        output_file << "# tau(fm)  e(GeV/fm^3)  rhob(1/fm^3)  "
                    << "u^tau  u^x  u^y  tau*u^eta  "
                    << "Bulk_Pi(GeV/fm^3)  theta(1/fm)  "
                    << "pi^{\\mu\\nu}(GeV/fm^3)  sigma^{\\mu\\nu}(1/fm)  "
                    << "omega^{\\mu\\nu}(1/fm)" << std::endl;
    } else {
        output_file.open(filename.str().c_str(),
                         std::fstream::out | std::fstream::app);
    }
    output_file << scientific << setprecision(8)
                << tau << "  " << arena_curr(ix,iy,ieta).epsilon*Util::hbarc
                << "  " << arena_curr(ix,iy,ieta).rhob << "  ";
    for (int i = 0; i < 4; i++) {
        output_file << scientific << setprecision(8)
                    << arena_curr(ix, iy, ieta).u[i] << "  ";
    }

    u_derivative_helper.MakedU(tau, arena_prev, arena_curr, ix, iy, ieta);
    auto theta_local = u_derivative_helper.calculate_expansion_rate(
                                            tau, arena_curr, ieta, ix, iy);
    output_file << scientific << setprecision(8)
                << arena_curr(ix, iy, ieta).pi_b*Util::hbarc << "  "
                << theta_local << "  ";
    DumuVec a_local;
    u_derivative_helper.calculate_Du_supmu(tau, arena_curr,
                                           ieta, ix, iy, a_local);
    VelocityShearVec sigma_local;
    u_derivative_helper.calculate_velocity_shear_tensor(
                    tau, arena_curr, ieta, ix, iy, a_local, sigma_local);
    VorticityVec omega_local;
    u_derivative_helper.calculate_kinetic_vorticity_with_spatial_projector(
                    tau, arena_curr, ieta, ix, iy, a_local, omega_local);
    for (int i = 0; i < 10; i++) {
        output_file << scientific << setprecision(8)
                    << arena_curr(ix, iy, ieta).Wmunu[i]*Util::hbarc << "  "
                    << sigma_local[i] << "  ";
    }
    for (int i = 0; i < 6; i++) {
        output_file << scientific << setprecision(8)
                    << omega_local[i] << "  ";
    }
    output_file << endl;
    output_file.close();
}

void Cell_info::output_vorticity_distribution(
                SCGrid &arena_curr, SCGrid &arena_prev, const double tau,
                const double eta_min, const double eta_max) {
    // This function outputs the vorticity tensor at a given tau
    ostringstream filename1;
    filename1 << "vorticity_dis_kinetic_wSP_eta_" << eta_min
              << "_" << eta_max << "_tau_" << tau << ".dat";
    std::fstream of1;
    of1.open(filename1.str().c_str(), std::fstream::out);
    // write the header
    of1 << "# x[fm]  y[fm]  T[GeV]  muB[GeV]  "
        << "omega^{tx}/T  omega^{ty}/T  "
        << "omega^{tz}/T  omega^{xy}/T  omega^{xz}/T  "
        << "omega^{yz}/T" << std::endl;
    ostringstream filename2;
    filename2 << "vorticity_dis_kinetic_eta_" << eta_min
              << "_" << eta_max << "_tau_" << tau << ".dat";
    std::fstream of2;
    of2.open(filename2.str().c_str(), std::fstream::out);
    // write the header
    of2 << "# x[fm]  y[fm]  T[GeV]  muB[GeV]  "
        << "omega^{tx}/T  omega^{ty}/T  "
        << "omega^{tz}/T  omega^{xy}/T  omega^{xz}/T  "
        << "omega^{yz}/T" << std::endl;
    ostringstream filename3;
    filename3 << "vorticity_dis_thermal_eta_" << eta_min
              << "_" << eta_max << "_tau_" << tau << ".dat";
    std::fstream of3;
    of3.open(filename3.str().c_str(), std::fstream::out);
    // write the header
    of3 << "# x[fm]  y[fm]  T[GeV]  muB[GeV]  "
        << "omega^{tx}  omega^{ty}  "
        << "omega^{tz}  omega^{xy}  omega^{xz}  "
        << "omega^{yz}" << std::endl;
    ostringstream filename4;
    filename4 << "vorticity_dis_T_eta_" << eta_min
              << "_" << eta_max << "_tau_" << tau << ".dat";
    std::fstream of4;
    of4.open(filename4.str().c_str(), std::fstream::out);
    // write the header
    of4 << "# x[fm]  y[fm]  T[GeV]  muB[GeV]  "
        << "omega^{tx}/T^2  omega^{ty}/T^2  "
        << "omega^{tz}/T^2  omega^{xy}/T^2  omega^{xz}/T^2  "
        << "omega^{yz}/T^2" << std::endl;

    for (int ix = 0; ix < arena_curr.nX(); ix++) {
        for (int iy = 0; iy < arena_curr.nY(); iy++) {
            const double x_local = -DATA.x_size/2. + ix*DATA.delta_x;
            const double y_local = -DATA.y_size/2. + iy*DATA.delta_y;

            VorticityVec omega_kSP = {0.0};
            VorticityVec omega_k   = {0.0};
            VorticityVec omega_th  = {0.0};
            VorticityVec omega_T   = {0.0};
            double T_avg = 0.0;
            double muB_avg = 0.0;
            double weight   = 0.0;
            for (int ieta = 0; ieta < arena_curr.nEta(); ieta++) {
                double eta_local = - DATA.eta_size/2. + ieta*DATA.delta_eta;
                if (DATA.boost_invariant)
                    eta_local = 0.0;
                if (eta_local < eta_max && eta_local > eta_min) {
                    const double e_local = arena_curr(ix, iy, ieta).epsilon;
                    if (e_local < 0.1) continue;
                    const double rhob_local = arena_curr(ix, iy, ieta).rhob;
                    const double T_local = (
                            eos.get_temperature(e_local, rhob_local));
                    const double muB_local = eos.get_muB(e_local, rhob_local);
                    T_avg += e_local*T_local*hbarc;
                    muB_avg += e_local*muB_local*hbarc;

                    VorticityVec omega_local_1, omega_local_2;
                    VorticityVec omega_local_3, omega_local_4;
                    u_derivative_helper.compute_vorticity_shell(
                        tau, arena_prev, arena_curr, ieta, ix, iy, eta_local,
                        omega_local_1, omega_local_2,
                        omega_local_3, omega_local_4);
                    for (unsigned int ii = 0; ii < omega_k.size(); ii++) {
                        omega_kSP[ii] += e_local*omega_local_1[ii]/T_local;
                        omega_k[ii]   += e_local*omega_local_2[ii]/T_local;
                        omega_th[ii]  += e_local*omega_local_3[ii];
                        omega_T[ii]   += (e_local*omega_local_4[ii]
                                          /T_local/T_local);
                    }
                    weight += e_local;
                }
            }
            weight = std::max(weight, small_eps);
            of1 << scientific << setprecision(8) << setw(18)
                << x_local << "  " << y_local << "  "
                << T_avg/weight << "  " << muB_avg/weight << "  ";
            of2 << scientific << setprecision(8) << setw(18)
                << x_local << "  " << y_local << "  "
                << T_avg/weight << "  " << muB_avg/weight << "  ";
            of3 << scientific << setprecision(8) << setw(18)
                << x_local << "  " << y_local << "  "
                << T_avg/weight << "  " << muB_avg/weight << "  ";
            of4 << scientific << setprecision(8) << setw(18)
                << x_local << "  " << y_local << "  "
                << T_avg/weight << "  " << muB_avg/weight << "  ";
            for (unsigned int ii = 0; ii < omega_k.size(); ii++) {
                // no minus sign because it has an opposite sign to
                // the kinetic vorticity
                of1 << scientific << setprecision(8) << setw(18)
                    << omega_kSP[ii]/weight << "  ";
                // minus sign from the metric
                // output quantities in g = (1, -1, -1 , -1)
                of2 << scientific << setprecision(8) << setw(18)
                    << -omega_k[ii]/weight << "  ";
                of3 << scientific << setprecision(8) << setw(18)
                    << -omega_th[ii]/weight << "  ";
                of4 << scientific << setprecision(8) << setw(18)
                    << -omega_T[ii]/weight << "  ";
            }
            of1 << std::endl;
            of2 << std::endl;
            of3 << std::endl;
            of4 << std::endl;
        }
    }
    of1.close();
    of2.close();
    of3.close();
    of4.close();
}

void Cell_info::output_vorticity_time_evolution(
                SCGrid &arena_curr, SCGrid &arena_prev, const double tau,
                const double eta_min, const double eta_max) {
    // This function outputs the time evolution of the vorticity tensor
    ostringstream filename1;
    filename1 << "vorticity_evo_kinetic_wSP_eta_" << eta_min
              << "_" << eta_max << ".dat";
    std::fstream of1;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        of1.open(filename1.str().c_str(), std::fstream::out);
        of1 << "# tau[fm]  omega^{tx}/T  omega^{ty}/T  "
            << "omega^{tz}/T  omega^{xy}/T  omega^{xz}/T  "
            << "omega^{yz}/T" << std::endl;
    } else {
        of1.open(filename1.str().c_str(),
                 std::fstream::out | std::fstream::app);
    }
    ostringstream filename2;
    filename2 << "vorticity_evo_kinetic_eta_" << eta_min
              << "_" << eta_max << ".dat";
    std::fstream of2;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        of2.open(filename2.str().c_str(), std::fstream::out);
        of2 << "# tau[fm]  omega^{tx}/T  omega^{ty}/T  "
            << "omega^{tz}/T  omega^{xy}/T  omega^{xz}/T  "
            << "omega^{yz}/T" << std::endl;
    } else {
        of2.open(filename2.str().c_str(),
                 std::fstream::out | std::fstream::app);
    }
    ostringstream filename3;
    filename3 << "vorticity_evo_thermal_eta_" << eta_min
              << "_" << eta_max << ".dat";
    std::fstream of3;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        of3.open(filename3.str().c_str(), std::fstream::out);
        of3 << "# tau[fm]  omega^{tx}  omega^{ty}  "
            << "omega^{tz}  omega^{xy}  omega^{xz}  "
            << "omega^{yz}" << std::endl;
    } else {
        of3.open(filename3.str().c_str(),
                 std::fstream::out | std::fstream::app);
    }
    ostringstream filename4;
    filename4 << "vorticity_evo_T_eta_" << eta_min
              << "_" << eta_max << ".dat";
    std::fstream of4;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        of4.open(filename4.str().c_str(), std::fstream::out);
        of4 << "# tau[fm]  omega^{tx}/T^2  omega^{ty}/T^2  "
            << "omega^{tz}/T^2  omega^{xy}/T^2  "
            << "omega^{xz}/T^2  omega^{yz}/T^2" << std::endl;
    } else {
        of4.open(filename4.str().c_str(),
                 std::fstream::out | std::fstream::app);
    }

    VorticityVec omega_kSP = {0.0};
    VorticityVec omega_k   = {0.0};
    VorticityVec omega_th  = {0.0};
    VorticityVec omega_T   = {0.0};
    double weight = 0.0;
    for (int ieta = 0; ieta < arena_curr.nEta(); ieta++) {
        double eta = 0.0;
        if (!DATA.boost_invariant) {
            eta = ((static_cast<double>(ieta))*(DATA.delta_eta)
                    - (DATA.eta_size)/2.0);
        }
        if (eta < eta_max && eta > eta_min) {
            for (int iy = 0; iy < arena_curr.nY(); iy++)
            for (int ix = 0; ix < arena_curr.nX(); ix++) {
                const double e_local = arena_curr(ix, iy, ieta).epsilon;
                if (e_local < 0.1) continue;
                const double rhob_local = arena_curr(ix, iy, ieta).rhob;
                const double T_local = (
                            eos.get_temperature(e_local, rhob_local));
                VorticityVec omega_local_1, omega_local_2;
                VorticityVec omega_local_3, omega_local_4;
                u_derivative_helper.compute_vorticity_shell(
                    tau, arena_prev, arena_curr, ieta, ix, iy, eta,
                    omega_local_1, omega_local_2,
                    omega_local_3, omega_local_4);
                for (unsigned int ii = 0; ii < omega_k.size(); ii++) {
                    omega_kSP[ii] += e_local*omega_local_1[ii]/T_local;
                    omega_k[ii]   += e_local*omega_local_2[ii]/T_local;
                    omega_th[ii]  += e_local*omega_local_3[ii];
                    omega_T[ii]   += e_local*omega_local_4[ii]/T_local/T_local;
                }
                weight += e_local;
            }
        }
    }
    weight = std::max(weight, small_eps);

    of1 << scientific << setw(18) << setprecision(8) << tau << "  ";
    of2 << scientific << setw(18) << setprecision(8) << tau << "  ";
    of3 << scientific << setw(18) << setprecision(8) << tau << "  ";
    of4 << scientific << setw(18) << setprecision(8) << tau << "  ";
    for (unsigned int ii = 0; ii < omega_k.size(); ii++) {
        // no minus sign because it has opposite sign to the kinetic vorcitity
        of1 << scientific << setprecision(8) << setw(18)
            << omega_kSP[ii]/weight << "  ";
        // minus sign from the metric
        // output quantities in g = (1, -1, -1 , -1)
        of2 << scientific << setprecision(8) << setw(18)
            << -omega_k[ii]/weight << "  ";
        of3 << scientific << setprecision(8) << setw(18)
            << -omega_th[ii]/weight << "  ";
        of4 << scientific << setprecision(8) << setw(18)
            << -omega_T[ii]/weight << "  ";
    }
    of1 << std::endl;
    of2 << std::endl;
    of3 << std::endl;
    of4 << std::endl;

    of1.close();
    of2.close();
    of3.close();
    of4.close();
}


void Cell_info::load_deltaf_qmu_coeff_table(string filename) {
    std::ifstream table(filename.c_str());
    deltaf_qmu_coeff_table_length_T = 150;
    deltaf_qmu_coeff_table_length_mu = 100;
    delta_qmu_coeff_table_T0 = 0.05;
    delta_qmu_coeff_table_mu0 = 0.0;
    delta_qmu_coeff_table_dT = 0.001;
    delta_qmu_coeff_table_dmu = 0.007892;
    deltaf_qmu_coeff_tb = new double* [deltaf_qmu_coeff_table_length_T];
    for (int i = 0; i < deltaf_qmu_coeff_table_length_T; i++)
       deltaf_qmu_coeff_tb[i] = new double[deltaf_qmu_coeff_table_length_mu];

    double dummy;
    for (int j = 0; j < deltaf_qmu_coeff_table_length_mu; j++)
       for (int i = 0; i < deltaf_qmu_coeff_table_length_T; i++)
          table >> dummy >> dummy >> deltaf_qmu_coeff_tb[i][j];
    table.close();
}


void Cell_info::load_deltaf_qmu_coeff_table_14mom(string filename) {
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
    for (int i = 0; i < deltaf_coeff_table_14mom_length_T; i++)
        for (int j = 0; j < deltaf_coeff_table_14mom_length_mu; j++)
            table >> dummy >> dummy >> deltaf_coeff_tb_14mom_DPi[i][j]
                  >> deltaf_coeff_tb_14mom_BPi[i][j]
                  >> deltaf_coeff_tb_14mom_BPitilde[i][j]
                  >> deltaf_coeff_tb_14mom_DV[i][j]
                  >> deltaf_coeff_tb_14mom_BV[i][j]
                  >> deltaf_coeff_tb_14mom_Bpi_shear[i][j];
    table.close();

    // convert units
    double hbarc3 = hbarc*hbarc*hbarc;
    double hbarc4 = hbarc3*hbarc;
    for (int i = 0; i < deltaf_coeff_table_14mom_length_T; i++) {
        for (int j = 0; j < deltaf_coeff_table_14mom_length_mu; j++) {
            deltaf_coeff_tb_14mom_DPi[i][j] =
                    deltaf_coeff_tb_14mom_DPi[i][j]*hbarc4;   // fm^4/GeV
            deltaf_coeff_tb_14mom_BPi[i][j] =
                    deltaf_coeff_tb_14mom_BPi[i][j]*hbarc4;   // fm^4/(GeV^2)
            deltaf_coeff_tb_14mom_BPitilde[i][j] =
                    deltaf_coeff_tb_14mom_BPitilde[i][j]*hbarc4;  // fm^4/(GeV^2)
            deltaf_coeff_tb_14mom_DV[i][j] =
                    deltaf_coeff_tb_14mom_DV[i][j]*hbarc3;   // fm^3/GeV
            deltaf_coeff_tb_14mom_BV[i][j] =
                    deltaf_coeff_tb_14mom_BV[i][j]*hbarc3;   // fm^3/(GeV^2)
            deltaf_coeff_tb_14mom_Bpi_shear[i][j] =
                    deltaf_coeff_tb_14mom_Bpi_shear[i][j]*hbarc4;  // fm^4/(GeV^2)
        }
    }
}


double Cell_info::get_deltaf_qmu_coeff(double T, double muB) {
    if (muB < 0) {
       muB = -muB;
    }
    int idx_T = static_cast<int>(
                    (T - delta_qmu_coeff_table_T0)/delta_qmu_coeff_table_dT);
    int idx_mu = static_cast<int>(
                (muB - delta_qmu_coeff_table_mu0)/delta_qmu_coeff_table_dmu);

    if (idx_T > deltaf_qmu_coeff_table_length_T - 2 || idx_T < 0)
        return(1.0);
    if (idx_mu > deltaf_qmu_coeff_table_length_mu - 2)
        return(1.0);

    double x_fraction = ((T - delta_qmu_coeff_table_T0)
                         /delta_qmu_coeff_table_dT - idx_T);
    double y_fraction = ((muB - delta_qmu_coeff_table_mu0)
                         /delta_qmu_coeff_table_dmu - idx_mu);

    double f1 = deltaf_qmu_coeff_tb[idx_T][idx_mu];
    double f2 = deltaf_qmu_coeff_tb[idx_T][idx_mu+1];
    double f3 = deltaf_qmu_coeff_tb[idx_T+1][idx_mu+1];
    double f4 = deltaf_qmu_coeff_tb[idx_T+1][idx_mu];

    double coeff = (f1*(1. - x_fraction)*(1. - y_fraction)
                    + f2*(1. - x_fraction)*y_fraction
                    + f3*x_fraction*y_fraction
                    + f4*x_fraction*(1. - y_fraction));
    return(coeff);
}


double Cell_info::get_deltaf_coeff_14moments(double T, double muB,
                                             double type) {
    int idx_T = static_cast<int>(
                (T - delta_coeff_table_14mom_T0)/delta_coeff_table_14mom_dT);
    int idx_mu = static_cast<int>(
            (muB - delta_coeff_table_14mom_mu0)/delta_coeff_table_14mom_dmu);
    double x_fraction = (
        (T - delta_coeff_table_14mom_T0)/delta_coeff_table_14mom_dT - idx_T);
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
       music_message << "Cell_info::get_deltaf_coeff_14moments: unknown type: "
                     << type;
       music_message.flush("error");
       exit(-1);
    }

    double f1 = deltaf_table[idx_T][idx_mu];
    double f2 = deltaf_table[idx_T][idx_mu+1];
    double f3 = deltaf_table[idx_T+1][idx_mu+1];
    double f4 = deltaf_table[idx_T+1][idx_mu];

    double coeff = (f1*(1. - x_fraction)*(1. - y_fraction)
                    + f2*(1. - x_fraction)*y_fraction
                    + f3*x_fraction*y_fraction
                    + f4*x_fraction*(1. - y_fraction));
    return(coeff);
}


//! This function outputs average T and mu_B as a function of proper tau
//! within a given space-time rapidity range
void Cell_info::output_average_phase_diagram_trajectory(
        const double tau, const double eta_min, const double eta_max,
        SCGrid &arena) {
    ostringstream filename;
    filename << "averaged_phase_diagram_trajectory_eta_" << eta_min
             << "_" << eta_max << ".dat";
    std::fstream of;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        of.open(filename.str().c_str(), std::fstream::out);
        of << "# tau(fm)  <T>(GeV)  std(T)(GeV)  <mu_B>(GeV)  std(mu_B)(GeV)  "
           << "V4 (fm^4)" << endl;
    } else {
        of.open(filename.str().c_str(),
                std::fstream::out | std::fstream::app);
    }
    double avg_T  = 0.0;
    double avg_mu = 0.0;
    double std_T  = 0.0;
    double std_mu = 0.0;
    double weight = 0.0;
    double V4     = 0.0;
    const double unit_volume = tau*DATA.delta_x*DATA.delta_y*DATA.delta_eta;
    for (int ieta = 0; ieta < arena.nEta(); ieta++) {
        double eta = 0.0;
        if (!DATA.boost_invariant) {
            eta = ((static_cast<double>(ieta))*(DATA.delta_eta)
                    - (DATA.eta_size)/2.0);
        }
        if (eta < eta_max && eta > eta_min) {
            double cosh_eta = cosh(eta);
            double sinh_eta = sinh(eta);
            for (int iy = 0; iy < arena.nY(); iy++)
            for (int ix = 0; ix < arena.nX(); ix++) {
                double e_local      = arena(ix, iy, ieta).epsilon;  // 1/fm^4
                if (e_local > 0.16/hbarc)
                    V4 += unit_volume;
                double rhob_local   = arena(ix, iy, ieta).rhob;     // 1/fm^3
                double utau         = arena(ix, iy, ieta).u[0];
                double ueta         = arena(ix, iy, ieta).u[3];
                double ut           = utau*cosh_eta + ueta*sinh_eta;  // gamma factor
                double T_local      = eos.get_temperature(e_local, rhob_local);
                double muB_local    = eos.get_muB(e_local, rhob_local);
                double weight_local = e_local*ut;
                avg_T  += T_local*weight_local;
                avg_mu += muB_local*weight_local;
                std_T  += T_local*T_local*weight_local;
                std_mu += muB_local*muB_local*weight_local;
                weight += weight_local;
            }
        }
    }
    avg_T  = avg_T/std::max(weight, small_eps)*hbarc;
    avg_mu = avg_mu/std::max(weight, small_eps)*hbarc;
    std_T  = sqrt(std_T/std::max(weight, small_eps)*hbarc*hbarc - avg_T*avg_T);
    std_mu = sqrt(std_mu/std::max(weight, small_eps)*hbarc*hbarc - avg_mu*avg_mu);
    of << scientific << setw(18) << setprecision(8)
       << tau << "  " << avg_T << "  " << std_T << "  "
       << avg_mu << "  " << std_mu << "  " << V4 << endl;
    of.close();
}


//! This function outputs system's eccentricity and momentum anisotropy
//! as functions of eta_s
void Cell_info::output_momentum_anisotropy_vs_etas(
                const double tau, SCGrid &arena) const {
    ostringstream filename;
    filename << "momentum_anisotropy_tau_" << tau << ".dat";
    std::fstream of;
    of.open(filename.str().c_str(), std::fstream::out);
    of << "# tau(fm)  epsilon_p(ideal)(cos)  epsilon_p(ideal)(sin)  "
       << "epsilon_p(shear)(cos)  epsilon_p(shear)(sin)  "
       << "epsilon_p(full)(cos)  epsilon_p(full)(sin)  " << endl;

    ostringstream filename1;
    filename1 << "eccentricities_evo_ed_tau_" << tau << ".dat";
    std::fstream of1;
    of1.open(filename1.str().c_str(), std::fstream::out);
    of1 << "# tau(fm)  ed(GeV/fm^3)  ecc_n(cos)  ecc_n(sin) (n=1-6)"<< endl;

    ostringstream filename2;
    filename2 << "eccentricities_evo_nB_tau_" << tau << ".dat";
    std::fstream of2;
    of2.open(filename2.str().c_str(), std::fstream::out);
    of2 << "# tau(fm)  nB(1/fm^3)  ecc_n(cos)  ecc_n(sin) (n=1-6)"<< endl;

    const int norder = 6;
    for (int ieta = 0; ieta < arena.nEta(); ieta++) {
        double eta = 0.0;
        if (!DATA.boost_invariant) {
            eta = ((static_cast<double>(ieta))*(DATA.delta_eta)
                    - (DATA.eta_size)/2.0);
        }

        // compute the central of mass position
        double x_ed_o = 0.0, x_nB_o = 0.0;
        double y_ed_o = 0.0, y_nB_o = 0.0;
        double w_ed_sum = 0.0, w_nB_sum = 0.0;
        for (int iy = 0; iy < arena.nY(); iy++)
        for (int ix = 0; ix < arena.nX(); ix++) {
            double x_local    = - DATA.x_size/2. + ix*DATA.delta_x;
            double y_local    = - DATA.y_size/2. + iy*DATA.delta_y;
            double e_local    = arena(ix, iy, ieta).epsilon;  // 1/fm^4
            double nB_local   = arena(ix, iy, ieta).rhob;     // 1/fm^3
            double gamma_perp = arena(ix, iy, ieta).u[0];
            x_ed_o   += x_local*e_local*gamma_perp;
            y_ed_o   += y_local*e_local*gamma_perp;
            w_ed_sum += e_local*gamma_perp;
            x_nB_o   += x_local*nB_local*gamma_perp;
            y_nB_o   += y_local*nB_local*gamma_perp;
            w_nB_sum += nB_local*gamma_perp;
        }
        x_ed_o /= std::max(small_eps, w_ed_sum);
        y_ed_o /= std::max(small_eps, w_ed_sum);
        x_nB_o /= std::max(small_eps, w_nB_sum);
        y_nB_o /= std::max(small_eps, w_nB_sum);

        // compute epsilon_p with ideal, ideal + shear, and full T^{\munu}
        std::vector<double> ep_num1(3, 0.0);
        std::vector<double> ep_num2(3, 0.0);
        std::vector<double> ep_den (3, 0.0);

        // spatial eccentricity arrays
        std::vector<double> eccn_ed_num1(norder, 0.0);
        std::vector<double> eccn_ed_num2(norder, 0.0);
        std::vector<double> eccn_ed_den (norder, 0.0);
        std::vector<double> eccn_nB_num1(norder, 0.0);
        std::vector<double> eccn_nB_num2(norder, 0.0);
        std::vector<double> eccn_nB_den (norder, 0.0);
        for (int iy = 0; iy < arena.nY(); iy++)
        for (int ix = 0; ix < arena.nX(); ix++) {
            double x_ed = - DATA.x_size/2. + ix*DATA.delta_x - x_ed_o;
            double y_ed = - DATA.y_size/2. + iy*DATA.delta_y - y_ed_o;
            double x_nB = - DATA.x_size/2. + ix*DATA.delta_x - x_nB_o;
            double y_nB = - DATA.y_size/2. + iy*DATA.delta_y - y_nB_o;
            double r_ed = sqrt(x_ed*x_ed + y_ed*y_ed);
            double r_nB = sqrt(x_nB*x_nB + y_nB*y_nB);
            double phi_ed = atan2(y_ed, x_ed);
            double phi_nB = atan2(y_nB, x_nB);

            double e_local    = arena(ix, iy, ieta).epsilon;  // 1/fm^4
            double rhob_local = arena(ix, iy, ieta).rhob;     // 1/fm^3
            double P_local    = eos.get_pressure(e_local, rhob_local);
            double enthopy    = e_local + P_local;
            double u0         = arena(ix, iy, ieta).u[0];
            double ux         = arena(ix, iy, ieta).u[1];
            double uy         = arena(ix, iy, ieta).u[2];
            double pi_xx      = arena(ix, iy, ieta).Wmunu[4];
            double pi_xy      = arena(ix, iy, ieta).Wmunu[5];
            double pi_yy      = arena(ix, iy, ieta).Wmunu[7];
            double bulk_Pi    = arena(ix, iy, ieta).pi_b;

            double T_xx_ideal   = enthopy*ux*ux + P_local;
            double T_xy_ideal   = enthopy*ux*uy;
            double T_yy_ideal   = enthopy*uy*uy + P_local;

            double T_xx_shear   = T_xx_ideal + pi_xx;
            double T_xy_shear   = T_xy_ideal + pi_xy;
            double T_yy_shear   = T_yy_ideal + pi_yy;

            double T_xx_full    = T_xx_shear - bulk_Pi*(-1 - ux*ux);
            double T_xy_full    = T_xy_shear + bulk_Pi*ux*uy;
            double T_yy_full    = T_yy_shear - bulk_Pi*(-1 - uy*uy);

            ep_num1[0] += T_xx_ideal - T_yy_ideal;
            ep_num2[0] += 2.*T_xy_ideal;
            ep_den[0]  += T_xx_ideal + T_yy_ideal;
            ep_num1[1] += T_xx_shear - T_yy_shear;
            ep_num2[1] += 2.*T_xy_shear;
            ep_den[1]  += T_xx_shear + T_yy_shear;
            ep_num1[2] += T_xx_full - T_yy_full;
            ep_num2[2] += 2.*T_xy_full;
            ep_den[2]  += T_xx_full + T_yy_full;

            double w_ed = 0.0, w_nB = 0.0;
            for (int i = 0; i < norder; i++) {
                if (i == 1) {
                    w_ed = u0*e_local*pow(r_ed, 3);
                    w_nB = u0*rhob_local*pow(r_nB, 3);
                } else {
                    w_ed = u0*e_local*pow(r_ed, i);
                    w_nB = u0*rhob_local*pow(r_nB, i);
                }
                if (i == 0) {
                    eccn_ed_num1[i] += e_local;
                    eccn_nB_num1[i] += rhob_local;
                } else {
                    eccn_ed_num1[i] += w_ed*cos(i*phi_ed);
                    eccn_ed_num2[i] += w_ed*sin(i*phi_ed);
                    eccn_ed_den [i] += w_ed;
                    eccn_nB_num1[i] += w_nB*cos(i*phi_nB);
                    eccn_nB_num2[i] += w_nB*sin(i*phi_nB);
                    eccn_nB_den [i] += w_nB;
                }
            }
        }

        // output results
        of << scientific << setw(18) << setprecision(8)
           << eta << "  ";
        for (int i = 0; i < 3; i++) {
            of << ep_num1[i]/std::max(ep_den[i], small_eps) << "  "
               << ep_num2[i]/std::max(ep_den[i], small_eps)<< "  ";
        }
        of << endl;
        of1 << scientific << setw(18) << setprecision(8)
            << eta << "  "
            << eccn_ed_num1[0]*DATA.delta_x*DATA.delta_y*hbarc << "  ";
        for (int i = 1; i < norder; i++) {
            // the minus sign ensure the vector points to the short axis
            of1 << -eccn_ed_num1[i]/std::max(eccn_ed_den[i], small_eps) << "  "
                << -eccn_ed_num2[i]/std::max(eccn_ed_den[i], small_eps) << "  ";
        }
        of1 << endl;
        of2 << scientific << setw(18) << setprecision(8)
            << eta << "  "
            << eccn_nB_num1[0]*DATA.delta_x*DATA.delta_y << "  ";
        for (int i = 1; i < norder; i++) {
            // the minus sign ensure the vector points to the short axis
            of2 << -eccn_nB_num1[i]/std::max(eccn_nB_den[i], small_eps) << "  "
                << -eccn_nB_num2[i]/std::max(eccn_nB_den[i], small_eps) << "  ";
        }
        of2 << endl;
    }
    of.close();
    of1.close();
    of2.close();
}


//! This function outputs system's momentum anisotropy as a function of tau
void Cell_info::output_momentum_anisotropy_vs_tau(
                const double tau, const double eta_min, const double eta_max,
                SCGrid &arena) const {
    ostringstream filename;
    filename << "momentum_anisotropy_eta_" << eta_min
             << "_" << eta_max << ".dat";
    std::fstream of;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        of.open(filename.str().c_str(), std::fstream::out);
        of << "# tau(fm)  epsilon_p(ideal)(cos)  epsilon_p(ideal)(sin)  "
           << "epsilon_p(shear)(cos)  epsilon_p(shear)(sin)  "
           << "epsilon_p(full)(cos)  epsilon_p(full)(sin)  "
           << "epsilon_2p(ideal)(cos)  epsilon_2p(ideal)(sin)  "
           << "epsilon_2p(shear)(cos)  epsilon_2p(shear)(sin)  "
           << "epsilon_2p(full)(cos)  epsilon_2p(full)(sin)  "
           << "epsilon_3p(ideal)(cos)  epsilon_3p(ideal)(sin)  "
           << "epsilon_3p(shear)(cos)  epsilon_3p(shear)(sin)  "
           << "epsilon_3p(full)(cos)  epsilon_3p(full)(sin)  "
           << endl;
    } else {
        of.open(filename.str().c_str(),
                std::fstream::out | std::fstream::app);
    }

    ostringstream filename1;
    filename1 << "eccentricities_evo_eta_" << eta_min
              << "_" << eta_max << ".dat";
    std::fstream of1;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        of1.open(filename1.str().c_str(), std::fstream::out);
        of1 << "# tau(fm)  ecc_n(cos)  ecc_n(sin) (n=1-6)"<< endl;
    } else {
        of1.open(filename1.str().c_str(),
                 std::fstream::out | std::fstream::app);
    }

    ostringstream filename2;
    filename2 << "inverse_Reynolds_number_eta_" << eta_min
             << "_" << eta_max << ".dat";
    std::fstream of2;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        of2.open(filename2.str().c_str(), std::fstream::out);
        of2 << "# tau(fm)  R_shearpi  R_Pi  gamma  T[GeV]" << endl;
    } else {
        of2.open(filename2.str().c_str(),
                 std::fstream::out | std::fstream::app);
    }

    double ideal_num1 = 0.0;
    double ideal_num2 = 0.0;
    double ideal_den  = 0.0;
    double shear_num1 = 0.0;
    double shear_num2 = 0.0;
    double shear_den  = 0.0;
    double full_num1  = 0.0;
    double full_num2  = 0.0;
    double full_den   = 0.0;
    double u_perp_num = 0.0;
    double u_perp_den = 0.0;
    double T_avg_num  = 0.0;
    double T_avg_den  = 0.0;
    double R_Pi_num   = 0.0;
    double R_Pi_den   = 0.0;
    double R_shearpi_num   = 0.0;
    double R_shearpi_den   = 0.0;

    // compute epsilon_{p2} and epsilon_{p3} using T^{0\mu} vector
    // epsilon_{pn} = (\int (T^{0r} exp(i n \phi_u))/(\int T^{0r})
    // for every n, we compute T^{0\mu} for ideal, ideal + shear, and full
    std::vector<double> ep_num1(6, 0.0);
    std::vector<double> ep_num2(6, 0.0);
    std::vector<double> ep_den (6, 0.0);

    const int norder = 6;
    std::vector<double> eccn_num1(norder, 0.0);
    std::vector<double> eccn_num2(norder, 0.0);
    std::vector<double> eccn_den (norder, 0.0);
    for (int ieta = 0; ieta < arena.nEta(); ieta++) {
        double eta = 0.0;
        if (!DATA.boost_invariant) {
            eta = ((static_cast<double>(ieta))*(DATA.delta_eta)
                    - (DATA.eta_size)/2.0);
        }
        if (eta < eta_max && eta > eta_min) {
            double x_o   = 0.0;
            double y_o   = 0.0;
            double w_sum = 0.0;
            for (int iy = 0; iy < arena.nY(); iy++)
            for (int ix = 0; ix < arena.nX(); ix++) {
                double x_local    = - DATA.x_size/2. + ix*DATA.delta_x;
                double y_local    = - DATA.y_size/2. + iy*DATA.delta_y;
                double e_local    = arena(ix, iy, ieta).epsilon;  // 1/fm^4
                double gamma_perp = arena(ix, iy, ieta).u[0];
                x_o   += x_local*e_local*gamma_perp;
                y_o   += y_local*e_local*gamma_perp;
                w_sum += e_local*gamma_perp;
            }
            x_o /= w_sum;
            y_o /= w_sum;
            for (int iy = 0; iy < arena.nY(); iy++)
            for (int ix = 0; ix < arena.nX(); ix++) {
                double x_local   = (- DATA.x_size/2. + ix*DATA.delta_x - x_o);
                double y_local   = (- DATA.y_size/2. + iy*DATA.delta_y - y_o);
                double r_local   = sqrt(x_local*x_local + y_local*y_local);
                double phi_local = atan2(y_local, x_local);

                double e_local    = arena(ix, iy, ieta).epsilon;  // 1/fm^4
                double rhob_local = arena(ix, iy, ieta).rhob;     // 1/fm^3
                double P_local    = eos.get_pressure(e_local, rhob_local);
                double enthopy    = e_local + P_local;
                double T_local    = eos.get_temperature(e_local, rhob_local);
                double u0         = arena(ix, iy, ieta).u[0];
                double ux         = arena(ix, iy, ieta).u[1];
                double uy         = arena(ix, iy, ieta).u[2];
                double pi_0x      = arena(ix, iy, ieta).Wmunu[1];
                double pi_0y      = arena(ix, iy, ieta).Wmunu[2];
                double pi_xx      = arena(ix, iy, ieta).Wmunu[4];
                double pi_xy      = arena(ix, iy, ieta).Wmunu[5];
                double pi_yy      = arena(ix, iy, ieta).Wmunu[7];
                double bulk_Pi    = arena(ix, iy, ieta).pi_b;

                double T_0x_ideal   = enthopy*u0*ux;
                double T_0y_ideal   = enthopy*u0*uy;
                double T_0r_ideal   = sqrt(  T_0x_ideal*T_0x_ideal
                                           + T_0y_ideal*T_0y_ideal);
                double phi_u_ideal  = atan2(T_0y_ideal, T_0x_ideal);
                double T_xx_ideal   = enthopy*ux*ux + P_local;
                double T_xy_ideal   = enthopy*ux*uy;
                double T_yy_ideal   = enthopy*uy*uy + P_local;

                double T_0x_shear   = T_0x_ideal + pi_0x;
                double T_0y_shear   = T_0y_ideal + pi_0y;
                double T_0r_shear   = sqrt(  T_0x_shear*T_0x_shear
                                           + T_0y_shear*T_0y_shear);
                double phi_u_shear  = atan2(T_0y_shear, T_0x_shear);
                double T_xx_shear   = T_xx_ideal + pi_xx;
                double T_xy_shear   = T_xy_ideal + pi_xy;
                double T_yy_shear   = T_yy_ideal + pi_yy;

                double T_0x_full    = T_0x_shear + bulk_Pi*u0*ux;
                double T_0y_full    = T_0y_shear + bulk_Pi*u0*uy;
                double T_0r_full    = sqrt(  T_0x_full*T_0x_full
                                           + T_0y_full*T_0y_full);
                double phi_u_full   = atan2(T_0y_full, T_0x_full);
                double T_xx_full    = T_xx_shear - bulk_Pi*(-1 - ux*ux);
                double T_xy_full    = T_xy_shear + bulk_Pi*ux*uy;
                double T_yy_full    = T_yy_shear - bulk_Pi*(-1 - uy*uy);

                ideal_num1 += T_xx_ideal - T_yy_ideal;
                ideal_num2 += 2.*T_xy_ideal;
                ideal_den  += T_xx_ideal + T_yy_ideal;
                shear_num1 += T_xx_shear - T_yy_shear;
                shear_num2 += 2.*T_xy_shear;
                shear_den  += T_xx_shear + T_yy_shear;
                full_num1  += T_xx_full - T_yy_full;
                full_num2  += 2.*T_xy_full;
                full_den   += T_xx_full + T_yy_full;

                double weight_local = e_local;
                u_perp_num += weight_local*u0;
                u_perp_den += weight_local;
                T_avg_num  += weight_local*T_local;
                T_avg_den  += weight_local;

                if (e_local > 1e-3) {
                    double r_shearpi_tmp, r_bulkPi_tmp;
                    calculate_inverse_Reynolds_numbers(arena, ieta, ix, iy,
                                                       r_shearpi_tmp,
                                                       r_bulkPi_tmp);
                    R_shearpi_num += weight_local*r_shearpi_tmp;
                    R_shearpi_den += weight_local;
                    R_Pi_num      += weight_local*r_bulkPi_tmp;
                    R_Pi_den      += weight_local;
                }

                for (int i = 0; i < 2; i++) {
                    int idx = 3*i;
                    int iorder = 2+i;
                    ep_num1[idx]   += T_0r_ideal*cos(iorder*phi_u_ideal);
                    ep_num2[idx]   += T_0r_ideal*sin(iorder*phi_u_ideal);
                    ep_den [idx]   += T_0r_ideal;
                    ep_num1[idx+1] += T_0r_shear*cos(iorder*phi_u_shear);
                    ep_num2[idx+1] += T_0r_shear*sin(iorder*phi_u_shear);
                    ep_den [idx+1] += T_0r_shear;
                    ep_num1[idx+2] += T_0r_full*cos(iorder*phi_u_full);
                    ep_num2[idx+2] += T_0r_full*sin(iorder*phi_u_full);
                    ep_den [idx+2] += T_0r_full;
                }
                for (int i = 1; i <= norder; i++) {
                    if (i == 1) {
                        weight_local = u0*e_local*pow(r_local, 3);
                    } else {
                        weight_local = u0*e_local*pow(r_local, i);
                    }
                    eccn_num1[i-1] += weight_local*cos(i*phi_local);
                    eccn_num2[i-1] += weight_local*sin(i*phi_local);
                    eccn_den [i-1] += weight_local;
                }
            }
        }
    }
    double R_shearpi = R_shearpi_num/std::max(R_shearpi_den, small_eps);
    double R_Pi      = R_Pi_num/std::max(R_Pi_den, small_eps);
    double u_avg     = u_perp_num/std::max(u_perp_den, small_eps);
    double T_avg     = T_avg_num/std::max(T_avg_den, small_eps)*hbarc;

    of << scientific << setw(18) << setprecision(8)
       << tau << "  "
       << ideal_num1/std::max(ideal_den, small_eps) << "  "
       << ideal_num2/std::max(ideal_den, small_eps) << "  "
       << shear_num1/std::max(shear_den, small_eps) << "  "
       << shear_num2/std::max(shear_den, small_eps) << "  "
       << full_num1/std::max(full_den, small_eps) << "  "
       << full_num2/std::max(full_den, small_eps) << "  ";
    for (int i = 0; i < 6; i++) {
        of << ep_num1[i]/std::max(ep_den[i], small_eps) << "  "
           << ep_num2[i]/std::max(ep_den[i], small_eps)<< "  ";
    }
    of << endl;
    of.close();

    of1 << scientific << setw(18) << setprecision(8)
        << tau << "  ";
    for (int i = 0; i < norder; i++) {
        // the minus sign ensure the vector points to the short axis
        of1 << -eccn_num1[i]/std::max(eccn_den[i], small_eps) << "  "
            << -eccn_num2[i]/std::max(eccn_den[i], small_eps)<< "  ";
    }
    of1 << endl;
    of1.close();

    of2 << scientific << setw(18) << setprecision(8)
        << tau << "  " << R_shearpi << "  " << R_Pi << "  "
        << u_avg << "  " << T_avg << endl;
    of2.close();
}


void Cell_info::get_LRF_shear_stress_tensor(const Cell_small &cell,
                                            const double eta_s,
                                            ShearVisVecLRF &res) {
    const double cosh_eta = cosh(eta_s);
    const double sinh_eta = sinh(eta_s);

    double u0 = cell.u[0];
    double ux = cell.u[1];
    double uy = cell.u[2];
    double u3 = cell.u[3];
    double ut = u0*cosh_eta + u3*sinh_eta;
    double uz = u3*cosh_eta + u0*sinh_eta;
    double LorentzBoost[4][4] = {
        {ut, -ux, -uy, -uz},
        {-ux, 1. + ux*ux/(ut+1.), ux*uy/(ut+1.), ux*uz/(ut+1.)},
        {-uy, ux*uy/(ut+1.), 1. + uy*uy/(ut+1.), uy*uz/(ut+1.)},
        {-uz, ux*uz/(ut+1.), uy*uz/(ut+1.), 1. + uz*uz/(ut+1.)}
    };

    auto ShearVisVec = cell.Wmunu;
    double pi_tz[4][4];
    pi_tz[0][0] = (  ShearVisVec[0]*cosh_eta*cosh_eta
                   + 2.*ShearVisVec[3]*cosh_eta*sinh_eta
                   + ShearVisVec[9]*sinh_eta*sinh_eta);
    pi_tz[0][1] = ShearVisVec[1]*cosh_eta + ShearVisVec[6]*sinh_eta;
    pi_tz[0][2] = ShearVisVec[2]*cosh_eta + ShearVisVec[8]*sinh_eta;
    pi_tz[0][3] = (  ShearVisVec[0]*cosh_eta*sinh_eta
                   + ShearVisVec[3]*(cosh_eta*cosh_eta
                                       + sinh_eta*sinh_eta)
                     + ShearVisVec[9]*sinh_eta*cosh_eta);
    pi_tz[1][0] = pi_tz[0][1];
    pi_tz[1][1] = ShearVisVec[4];
    pi_tz[1][2] = ShearVisVec[5];
    pi_tz[1][3] = ShearVisVec[1]*sinh_eta + ShearVisVec[6]*cosh_eta;
    pi_tz[2][0] = pi_tz[0][2];
    pi_tz[2][1] = pi_tz[1][2];
    pi_tz[2][2] = ShearVisVec[7];
    pi_tz[2][3] = ShearVisVec[2]*sinh_eta + ShearVisVec[8]*cosh_eta;
    pi_tz[3][0] = pi_tz[0][3];
    pi_tz[3][1] = pi_tz[1][3];
    pi_tz[3][2] = pi_tz[2][3];
    pi_tz[3][3] = pi_tz[0][0] - pi_tz[1][1] - pi_tz[2][2];
    double pi_LRF[4][4];
    for (int i = 1; i < 3; i++) {
        for (int j = i; j < 4; j++) {
            pi_LRF[i][j] = 0.;
            for (int a = 0; a < 4; a++) {
                for (int b = 0; b < 4; b++) {
                    pi_LRF[i][j] += (LorentzBoost[i][a]*pi_tz[a][b]
                                     *LorentzBoost[b][j]);
                }
            }
        }
    }
    res[0] = pi_LRF[1][1];
    res[1] = pi_LRF[1][2];
    res[2] = pi_LRF[1][3];
    res[3] = pi_LRF[2][2];
    res[4] = pi_LRF[2][3];

    double q_tz[4];
    q_tz[0] = ShearVisVec[10]*cosh_eta + ShearVisVec[13]*sinh_eta;
    q_tz[1] = ShearVisVec[11];
    q_tz[2] = ShearVisVec[12];
    q_tz[3] = ShearVisVec[10]*sinh_eta + ShearVisVec[13]*cosh_eta;
    double q_LRF[4];
    for (int i = 1; i < 4; i++) {
        q_LRF[i] = 0.;
        for (int a = 0; a < 4; a++)
            q_LRF[i] += LorentzBoost[i][a]*q_tz[a];
    }
    res[5] = q_LRF[1];
    res[6] = q_LRF[2];
    res[7] = q_LRF[3];
}
