// Copyright Chun Shen @ 2014-2016
#include "grid_info.h"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "util.h"

using std::endl;
using std::ofstream;
using std::ostringstream;
using std::scientific;
using std::setprecision;
using std::setw;
using std::string;
using Util::hbarc;
using Util::small_eps;

Cell_info::Cell_info(const InitData &DATA_in, const EOS &eos_in)
    : DATA(DATA_in), eos(eos_in), u_derivative_helper(DATA_in, eos_in) {
    // read in tables for delta f coefficients
    if (DATA.turn_on_diff) {
        if (DATA.deltaf_14moments) {
            load_deltaf_qmu_coeff_table_14mom(
                "tables/deltaf_coefficients_14moments.dat");
        } else {
            if (DATA.include_deltaf_qmu)
                load_deltaf_qmu_coeff_table(
                    "tables/Coefficients_RTA_diffusion.dat");
        }
    }
    outflow_flux = {0.};
}

Cell_info::~Cell_info() {
    if (DATA.turn_on_diff) {
        if (DATA.deltaf_14moments) {
            for (int i = 0; i < deltaf_coeff_table_14mom_length_T; i++) {
                delete[] deltaf_coeff_tb_14mom_DPi[i];
                delete[] deltaf_coeff_tb_14mom_BPi[i];
                delete[] deltaf_coeff_tb_14mom_BPitilde[i];
                delete[] deltaf_coeff_tb_14mom_DV[i];
                delete[] deltaf_coeff_tb_14mom_BV[i];
                delete[] deltaf_coeff_tb_14mom_Bpi_shear[i];
            }
            delete[] deltaf_coeff_tb_14mom_DPi;
            delete[] deltaf_coeff_tb_14mom_BPi;
            delete[] deltaf_coeff_tb_14mom_BPitilde;
            delete[] deltaf_coeff_tb_14mom_DV;
            delete[] deltaf_coeff_tb_14mom_BV;
            delete[] deltaf_coeff_tb_14mom_Bpi_shear;
        } else {
            if (DATA.include_deltaf_qmu) {
                for (int i = 0; i < deltaf_qmu_coeff_table_length_T; i++)
                    delete[] deltaf_qmu_coeff_tb[i];
                delete[] deltaf_qmu_coeff_tb;
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
        (static_cast<double>(DATA.nx + 1)) / DATA.output_evolution_every_N_x);
    int grid_ny = ceil(
        (static_cast<double>(DATA.ny + 1)) / DATA.output_evolution_every_N_y);
    int grid_neta = ceil(
        (static_cast<double>(DATA.neta)) / DATA.output_evolution_every_N_eta);

    outfile << "const int MUSIC_real_nx = " << grid_nx << ";" << endl;
    outfile << "const int MUSIC_real_ny = " << grid_ny << ";" << endl;
    outfile << "const int MUSIC_real_neta = " << grid_neta << ";" << endl;

    outfile << "const double MUSIC_tau0 = " << DATA.tau0 << ";" << endl;

    outfile << "const double MUSIC_dx = "
            << DATA.delta_x * DATA.output_evolution_every_N_x << ";" << endl;
    outfile << "const double MUSIC_dy = "
            << DATA.delta_y * DATA.output_evolution_every_N_y << ";" << endl;
    outfile << "const double MUSIC_deta = "
            << DATA.delta_eta * DATA.output_evolution_every_N_eta << ";"
            << endl;
    outfile << "const double MUSIC_dtau = "
            << DATA.output_evolution_every_N_timesteps * DATA.delta_tau << ";"
            << endl;

    outfile << "const bool MUSIC_with_shear_viscosity = " << DATA.turn_on_shear
            << ";\n";
    outfile << "const bool MUSIC_with_bulk_viscosity = " << DATA.turn_on_bulk
            << ";\n";
    outfile << "const bool MUSIC_with_rhob = " << (DATA.turn_on_rhob) << ";\n";
    outfile << "const bool MUSIC_with_diffusion = " << DATA.turn_on_diff
            << ";\n";
    outfile << "const bool MUSIC_outputBinaryEvolution="
            << DATA.outputBinaryEvolution << ";" << endl;

    outfile.close();
}

//! This function outputs hydro evolution file in binary format
void Cell_info::OutputEvolutionDataXYEta(Fields &arena, double tau) {
    const string out_name_xyeta = "evolution_xyeta.dat";
    const string out_name_W_xyeta =
        "evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
    const string out_name_bulkpi_xyeta = "evolution_bulk_pressure_xyeta.dat";
    const string out_name_q_xyeta = "evolution_qmu_xyeta.dat";
    string out_open_mode;
    FILE *out_file_xyeta = NULL;
    FILE *out_file_W_xyeta = NULL;
    FILE *out_file_bulkpi_xyeta = NULL;
    FILE *out_file_q_xyeta = NULL;

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
    if (DATA.turn_on_shear) {
        out_file_W_xyeta =
            fopen(out_name_W_xyeta.c_str(), out_open_mode.c_str());
    }
    if (DATA.turn_on_bulk) {
        out_file_bulkpi_xyeta =
            fopen(out_name_bulkpi_xyeta.c_str(), out_open_mode.c_str());
    }
    if (DATA.turn_on_diff) {
        out_file_q_xyeta =
            fopen(out_name_q_xyeta.c_str(), out_open_mode.c_str());
    }
    const int n_skip_x = DATA.output_evolution_every_N_x;
    const int n_skip_y = DATA.output_evolution_every_N_y;
    const int n_skip_eta = DATA.output_evolution_every_N_eta;
    ThermalVec thermalVec;
    for (int ieta = 0; ieta < arena.nEta(); ieta += n_skip_eta) {
        double eta = 0.0;
        if (!DATA.boost_invariant) {
            eta = DATA.etaAt(ieta);
        }
        double cosh_eta = cosh(eta);
        double sinh_eta = sinh(eta);
        for (int iy = 0; iy < arena.nY(); iy += n_skip_y) {
            for (int ix = 0; ix < arena.nX(); ix += n_skip_x) {
                int fieldIdx = arena.getFieldIdx(ix, iy, ieta);
                eos.getThermalVariables(
                    arena.e_[fieldIdx], arena.rhob_[fieldIdx],
                    arena.rhoq_[fieldIdx], arena.rhos_[fieldIdx], thermalVec);
                double utau = arena.u_[0][fieldIdx];
                double ux = arena.u_[1][fieldIdx];
                double uy = arena.u_[2][fieldIdx];
                double ueta = arena.u_[3][fieldIdx];
                double ut = utau * cosh_eta + ueta * sinh_eta;  // gamma factor
                double vx = ux / ut;
                double vy = uy / ut;
                double uz = ueta * cosh_eta + utau * sinh_eta;
                double vz = uz / ut;

                // Entropy with other conserved charges ??
                double enthropy = thermalVec[0] + thermalVec[2];  // [1/fm^4]
                // e + P >= 0, = 0 only in empty cells
                const double enthalpyDiv = std::max(enthropy, Util::tiny_normal);
                double trace_anomaly = thermalVec[0] - 3.0 * thermalVec[2];

                double Wtautau = 0.0;
                double Wtaux = 0.0;
                double Wtauy = 0.0;
                double Wtaueta = 0.0;
                double Wxx = 0.0;
                double Wxy = 0.0;
                double Wxeta = 0.0;
                double Wyy = 0.0;
                double Wyeta = 0.0;
                double Wetaeta = 0.0;
                if (DATA.turn_on_shear) {
                    Wtautau = arena.Wmunu_[0][fieldIdx] / enthalpyDiv;
                    Wtaux = arena.Wmunu_[1][fieldIdx] / enthalpyDiv;
                    Wtauy = arena.Wmunu_[2][fieldIdx] / enthalpyDiv;
                    Wtaueta = arena.Wmunu_[3][fieldIdx] / enthalpyDiv;
                    Wxx = arena.Wmunu_[4][fieldIdx] / enthalpyDiv;
                    Wxy = arena.Wmunu_[5][fieldIdx] / enthalpyDiv;
                    Wxeta = arena.Wmunu_[6][fieldIdx] / enthalpyDiv;
                    Wyy = arena.Wmunu_[7][fieldIdx] / enthalpyDiv;
                    Wyeta = arena.Wmunu_[8][fieldIdx] / enthalpyDiv;
                    Wetaeta = arena.Wmunu_[9][fieldIdx] / enthalpyDiv;
                }

                double bulk_Pi = 0.0;
                double bulk_Pi_kinetic = 0.0;
                double bulk_Pi_chem = 0.0;
                if (DATA.turn_on_bulk) {
                    // Output separate components for analysis
                    bulk_Pi_kinetic = arena.piBulk_[fieldIdx];   // [1/fm^4]
                    bulk_Pi_chem = arena.piBulkChem_[fieldIdx];  // [1/fm^4]
                    bulk_Pi = bulk_Pi_kinetic + bulk_Pi_chem;  // Total [1/fm^4]
                }

                // outputs for baryon diffusion part
                double common_term_q = 0.0;
                double qtau = 0.0;
                double qx = 0.0;
                double qy = 0.0;
                double qeta = 0.0;
                if (DATA.turn_on_diff) {
                    common_term_q = thermalVec[1] * thermalVec[6] / enthalpyDiv;
                    double kappa_hat =
                        get_deltaf_qmu_coeff(thermalVec[6], thermalVec[7]);
                    qtau = arena.Wmunu_[10][fieldIdx] / kappa_hat;
                    qx = arena.Wmunu_[11][fieldIdx] / kappa_hat;
                    qy = arena.Wmunu_[12][fieldIdx] / kappa_hat;
                    qeta = arena.Wmunu_[13][fieldIdx] / kappa_hat;
                }

                // exclude the actual coordinates from the output to save space:
                if (!DATA.outputBinaryEvolution) {
                    fprintf(
                        out_file_xyeta, "%e %e %e %e %e\n",
                        thermalVec[6] * hbarc, thermalVec[7] * hbarc, vx, vy,
                        vz);
                    if (DATA.viscosity_flag) {
                        if (DATA.turn_on_shear) {
                            fprintf(
                                out_file_W_xyeta,
                                "%e %e %e %e %e %e %e %e %e %e\n", Wtautau,
                                Wtaux, Wtauy, Wtaueta, Wxx, Wxy, Wxeta, Wyy,
                                Wyeta, Wetaeta);
                        }
                        if (DATA.turn_on_bulk) {
                            fprintf(
                                out_file_bulkpi_xyeta, "%e %e %e %e %e %e %e\n",
                                tau, bulk_Pi_kinetic, bulk_Pi_chem, bulk_Pi,
                                enthropy, thermalVec[5], trace_anomaly);
                        }
                    }
                } else {
                    float array[] = {
                        static_cast<float>(thermalVec[6] * hbarc),
                        static_cast<float>(thermalVec[7] * hbarc),
                        static_cast<float>(vx), static_cast<float>(vy),
                        static_cast<float>(vz)};
                    fwrite(array, sizeof(float), 5, out_file_xyeta);
                    if (DATA.viscosity_flag) {
                        if (DATA.turn_on_shear) {
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
                            fwrite(array2, sizeof(float), 10, out_file_W_xyeta);
                        }
                        if (DATA.turn_on_bulk) {
                            float array1[] = {
                                static_cast<float>(tau),
                                static_cast<float>(bulk_Pi_kinetic),
                                static_cast<float>(bulk_Pi_chem),
                                static_cast<float>(bulk_Pi),
                                static_cast<float>(enthropy),
                                static_cast<float>(thermalVec[5]),
                                static_cast<float>(trace_anomaly)};
                            fwrite(
                                array1, sizeof(float), 7,
                                out_file_bulkpi_xyeta);
                        }
                        if (DATA.turn_on_diff) {
                            float array3[] = {
                                static_cast<float>(common_term_q),
                                static_cast<float>(qtau),
                                static_cast<float>(qx), static_cast<float>(qy),
                                static_cast<float>(qeta)};
                            fwrite(array3, sizeof(float), 5, out_file_q_xyeta);
                        }
                    }
                }
            }
        }
    }
    fclose(out_file_xyeta);
    if (DATA.turn_on_shear) {
        fclose(out_file_W_xyeta);
    }
    if (DATA.turn_on_bulk) {
        fclose(out_file_bulkpi_xyeta);
    }
    if (DATA.turn_on_diff) {
        fclose(out_file_q_xyeta);
    }
}

void Cell_info::OutputEvolution_Knudsen_Reynoldsnumbers(
    Fields &arena, const double tau) const {
    const string out_name_xyeta = "evolution_KRnumbers.dat";
    FILE *out_file_xyeta = NULL;

    // If it's the first timestep, overwrite the previous file
    string out_open_mode;
    if (tau == DATA.tau0) {
        out_open_mode = "w";
    } else {
        out_open_mode = "a";
    }

    // If we output in binary, set the mode accordingly
    if (!DATA.outputBinaryEvolution) {
        out_open_mode += "b";
    }
    out_file_xyeta = fopen(out_name_xyeta.c_str(), out_open_mode.c_str());

    const int n_skip_x = DATA.output_evolution_every_N_x;
    const int n_skip_y = DATA.output_evolution_every_N_y;
    const int n_skip_eta = DATA.output_evolution_every_N_eta;
    for (int ieta = 0; ieta < arena.nEta(); ieta += n_skip_eta) {
        for (int iy = 0; iy < arena.nY(); iy += n_skip_y) {
            for (int ix = 0; ix < arena.nX(); ix += n_skip_x) {
                int fieldIdx = arena.getFieldIdx(ix, iy, ieta);

                double P_local = eos.get_pressure(
                    arena.e_[fieldIdx], arena.rhob_[fieldIdx],
                    arena.rhoq_[fieldIdx], arena.rhos_[fieldIdx]);
                double R_pi = 0.0;
                double R_Pi = 0.0;
                calculate_inverse_Reynolds_numbers(
                    arena, fieldIdx, arena.e_[fieldIdx] + P_local, R_pi, R_Pi);

                if (!DATA.outputBinaryEvolution) {
                    fprintf(
                        out_file_xyeta, "%e %e %e\n", arena.e_[fieldIdx], R_pi,
                        R_Pi);
                } else {
                    float array[] = {
                        static_cast<float>(arena.e_[fieldIdx]),
                        static_cast<float>(R_pi), static_cast<float>(R_Pi)};
                    fwrite(array, sizeof(float), 3, out_file_xyeta);
                }
            }
        }
    }
    fclose(out_file_xyeta);
}

void Cell_info::calculate_inverse_Reynolds_numbers(
    Fields &arena, const int Idx, const double P_local, double &R_pi,
    double &R_Pi) const {
    const double pi_00 = arena.Wmunu_[0][Idx];
    const double pi_01 = arena.Wmunu_[1][Idx];
    const double pi_02 = arena.Wmunu_[2][Idx];
    const double pi_03 = arena.Wmunu_[3][Idx];
    const double pi_11 = arena.Wmunu_[4][Idx];
    const double pi_12 = arena.Wmunu_[5][Idx];
    const double pi_13 = arena.Wmunu_[6][Idx];
    const double pi_22 = arena.Wmunu_[7][Idx];
    const double pi_23 = arena.Wmunu_[8][Idx];
    const double pi_33 = arena.Wmunu_[9][Idx];

    const double pisize =
        (pi_00 * pi_00 + pi_11 * pi_11 + pi_22 * pi_22 + pi_33 * pi_33
         - 2. * (pi_01 * pi_01 + pi_02 * pi_02 + pi_03 * pi_03)
         + 2. * (pi_12 * pi_12 + pi_13 * pi_13 + pi_23 * pi_23));

    // P_local is e + P >= 0 here, = 0 only in empty cells
    const double enthalpy = std::max(P_local, Util::tiny_normal);
    R_pi = sqrt(pisize) / enthalpy;
    R_Pi = arena.piBulk_[Idx] / enthalpy;
}

//! This function outputs hydro evolution file into memory for JETSCAPE
void Cell_info::OutputEvolutionDataXYEta_memory(
    Fields &arena, const double tau, HydroinfoMUSIC &hydro_info_ptr) {
    const int n_skip_x = DATA.output_evolution_every_N_x;
    const int n_skip_y = DATA.output_evolution_every_N_y;
    const int n_skip_eta = DATA.output_evolution_every_N_eta;
    for (int ix = 0; ix < arena.nX(); ix += n_skip_x) {
        for (int iy = 0; iy < arena.nY(); iy += n_skip_y) {
            for (int ieta = 0; ieta < arena.nEta(); ieta += n_skip_eta) {
                int fieldIdx = arena.getFieldIdx(ix, iy, ieta);
                double eta = 0.0;
                if (!DATA.boost_invariant) {
                    eta = DATA.etaAt(ieta);
                }
                double e_local = arena.e_[fieldIdx];        // 1/fm^4
                double rhob_local = arena.rhob_[fieldIdx];  // 1/fm^3
                double rhoq_local = arena.rhoq_[fieldIdx];  // 1/fm^3
                double rhos_local = arena.rhos_[fieldIdx];  // 1/fm^3
                double p_local = eos.get_pressure(
                    e_local, rhob_local, rhoq_local, rhos_local);
                double T_local = eos.get_temperature(
                    e_local, rhob_local, rhoq_local, rhos_local);
                double s_local = eos.get_entropy(
                    e_local, rhob_local, rhoq_local, rhos_local);
                // store the Milne components u^x, u^y, tau*u^eta;
                // HydroinfoMUSIC converts them to Cartesian velocities
                double ux = arena.u_[1][fieldIdx];
                double uy = arena.u_[2][fieldIdx];
                double ueta = arena.u_[3][fieldIdx];

                hydro_info_ptr.dump_ideal_info_to_memory(
                    tau, eta, e_local, p_local, s_local, T_local, ux, uy, ueta);
            }
        }
    }
}

//! This function outputs hydro evolution file in binary format
int Cell_info::OutputEvolutionDataXYEta_chun(Fields &arena, double tau) {
    // the format of the file is as follows,
    //    itau ix iy ieta e P T cs^2 ux uy ueta
    // if turn_on_shear == 1:
    //    itau ix iy ieta e P T cs^2 ux uy ueta Wxx Wxy Wxeta Wyy Wyeta
    // if turn_on_shear == 1 and turn_on_bulk == 1:
    //    itau ix iy ieta e P T cs^2 ux uy ueta Wxx Wxy Wxeta Wyy Wyeta pi_b
    // if turn_on_rhob == 1:
    //    itau ix iy ieta e P T cs^2 ux uy ueta rho_B mu_B
    // if turn_on_rhob == 1 and turn_on_shear == 1:
    //    itau ix iy ieta e P T cs^2 ux uy ueta rho_B mu_B Wxx Wxy Wxeta Wyy
    //    Wyeta
    // if turn_on_rhob == 1 and turn_on_shear == 1 and turn_on_diff == 1:
    //    itau ix iy ieta e P T cs^2 ux uy ueta rho_B mu_B Wxx Wxy Wxeta Wyy
    //    Wyeta qx qy qeta
    // if turn_on_rhob == 1 and turn_on_shear == 1 and turn_on_bulk == 1 and
    // turn_on_diff == 1:
    //    itau ix iy ieta e P T cs^2 ux uy ueta rho_B mu_B Wxx Wxy Wxeta Wyy
    //    Wyeta pi_b qx qy qeta
    // Here ueta = tau*ueta, Wieta = tau*Wieta, Wetaeta = tau^2*Wetaeta, qeta =
    // tau*qeta Here Wij is reduced variables Wij/(e+P) in the fluid rest frame
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

    int n_skip_tau = DATA.output_evolution_every_N_timesteps;
    double output_dtau = DATA.delta_tau * n_skip_tau;
    int itau = static_cast<int>((tau - DATA.tau0) / (output_dtau) + 0.1);

    int n_skip_x = DATA.output_evolution_every_N_x;
    int n_skip_y = DATA.output_evolution_every_N_y;
    int n_skip_eta = DATA.output_evolution_every_N_eta;

    // write out header
    const int output_nx = static_cast<int>(arena.nX() / n_skip_x);
    const int output_ny = static_cast<int>(arena.nY() / n_skip_y);
    const int output_neta = static_cast<int>(arena.nEta() / n_skip_eta);
    const double output_dx = DATA.delta_x * n_skip_x;
    const double output_dy = DATA.delta_y * n_skip_y;
    const double output_deta = DATA.delta_eta * n_skip_eta;
    const double output_xmin = -DATA.x_size / 2.;
    const double output_ymin = -DATA.y_size / 2.;
    const double output_etamin = DATA.etaMin();

    int nVar_per_cell =
        (11 + DATA.turn_on_rhob * 2 + DATA.turn_on_QS * 4
         + DATA.turn_on_shear * 5 + DATA.turn_on_bulk * 1
         + DATA.turn_on_diff * 3);
    if (DATA.output_evolution_ideal_only) {
        nVar_per_cell = 11 + DATA.turn_on_rhob * 2 + DATA.turn_on_QS * 4;
    }

    if (tau == DATA.tau0) {
        float header[] = {static_cast<float>(DATA.tau0),
                          static_cast<float>(output_dtau),
                          static_cast<float>(output_nx),
                          static_cast<float>(output_dx),
                          static_cast<float>(output_xmin),
                          static_cast<float>(output_ny),
                          static_cast<float>(output_dy),
                          static_cast<float>(output_ymin),
                          static_cast<float>(output_neta),
                          static_cast<float>(output_deta),
                          static_cast<float>(output_etamin),
                          static_cast<float>(DATA.turn_on_rhob),
                          static_cast<float>(DATA.turn_on_shear),
                          static_cast<float>(DATA.turn_on_bulk),
                          static_cast<float>(DATA.turn_on_diff),
                          static_cast<float>(nVar_per_cell),
                          static_cast<float>(DATA.turn_on_QS)};
        if (DATA.output_evolution_ideal_only) {
            header[12] = 0.0;  // turn on shear
            header[13] = 0.0;  // turn on bulk
            header[14] = 0.0;  // turn on diff
        }
        fwrite(header, sizeof(float), 17, out_file_xyeta);
    }
    ThermalVec thermalVec;
    for (int ieta = 0; ieta < arena.nEta(); ieta += n_skip_eta) {
        double eta_local = DATA.etaAt(ieta);
        double cosh_eta = cosh(eta_local);
        double sinh_eta = sinh(eta_local);
        for (int iy = 0; iy < arena.nY(); iy += n_skip_y) {
            for (int ix = 0; ix < arena.nX(); ix += n_skip_x) {
                int fieldIdx = arena.getFieldIdx(ix, iy, ieta);
                double e_local = arena.e_[fieldIdx];        // 1/fm^4
                double rhob_local = arena.rhob_[fieldIdx];  // 1/fm^3
                double rhoq_local = arena.rhoq_[fieldIdx];  // 1/fm^3
                double rhos_local = arena.rhos_[fieldIdx];  // 1/fm^3

                if (e_local * hbarc < DATA.output_evolution_e_cut) continue;
                // only ouput fluid cells that are above cut-off temperature
                // check whether the cells are reaching the grid edges
                if (ix == 0 || ix == arena.nX() - 1 || iy == 0
                    || iy == arena.nY() - 1) {
                    std::cout << "Physical cells with e > "
                              << DATA.output_evolution_e_cut
                              << " GeV/fm^3, are reaching the grid edge! "
                              << "ix = " << ix << ", iy = " << iy << std::endl;
                    fclose(out_file_xyeta);
                    return (-1);
                }

                // Adapted for 4D EoS. ThermalVec shape has not been changed.
                eos.getThermalVariables(
                    e_local, rhob_local, rhoq_local, rhos_local, thermalVec);
                double muB_local = thermalVec[7];
                double muQ_local = thermalVec[10];
                double muS_local = thermalVec[8];
                double p_local = thermalVec[2];
                double cs2 = thermalVec[5];

                double ux = arena.u_[1][fieldIdx];
                double uy = arena.u_[2][fieldIdx];
                double uz =
                    (arena.u_[3][fieldIdx] * cosh_eta
                     + arena.u_[0][fieldIdx] * sinh_eta);

                // T_local is in 1/fm
                double T_local = thermalVec[6];

                ShearVisVecLRF piLRF;
                get_LRF_shear_stress_tensor(arena, fieldIdx, eta_local, piLRF);
                // e + P >= 0, = 0 only in empty cells
                double div_factor =
                    std::max(e_local + p_local, Util::tiny_normal);  // 1/fm^4
                double Wxx = 0.0;
                double Wxy = 0.0;
                double Wxz = 0.0;
                double Wyy = 0.0;
                double Wyz = 0.0;
                if (DATA.turn_on_shear) {
                    Wxx = piLRF[0] / div_factor;
                    Wxy = piLRF[1] / div_factor;
                    Wxz = piLRF[2] / div_factor;
                    Wyy = piLRF[3] / div_factor;
                    Wyz = piLRF[4] / div_factor;
                }

                double pi_b = 0.0;
                if (DATA.turn_on_bulk) {
                    pi_b = arena.piBulk_[fieldIdx] / div_factor;
                }

                // outputs for baryon diffusion part
                // double common_term_q = 0.0;
                double qx = 0.0;
                double qy = 0.0;
                double qz = 0.0;
                if (DATA.turn_on_diff) {
                    // common_term_q = rhob_local*T_local/div_factor;
                    double kappa_hat = get_deltaf_qmu_coeff(T_local, muB_local);
                    qx = piLRF[5] / kappa_hat;
                    qy = piLRF[6] / kappa_hat;
                    qz = piLRF[7] / kappa_hat;
                }

                float ideal[] = {
                    static_cast<float>(itau),
                    static_cast<float>(ix / n_skip_x),
                    static_cast<float>(iy / n_skip_y),
                    static_cast<float>(ieta / n_skip_eta),
                    static_cast<float>(e_local * hbarc),
                    static_cast<float>(p_local * hbarc),
                    static_cast<float>(T_local * hbarc),
                    static_cast<float>(cs2),
                    static_cast<float>(ux),
                    static_cast<float>(uy),
                    static_cast<float>(uz)};

                fwrite(ideal, sizeof(float), 11, out_file_xyeta);

                if (DATA.turn_on_rhob) {
                    float mu_B[] = {
                        static_cast<float>(rhob_local),
                        static_cast<float>(muB_local * hbarc)};
                    fwrite(mu_B, sizeof(float), 2, out_file_xyeta);
                }

                if (!DATA.output_evolution_ideal_only) {
                    if (DATA.turn_on_shear) {
                        float shear_pi[] = {
                            static_cast<float>(Wxx), static_cast<float>(Wxy),
                            static_cast<float>(Wxz), static_cast<float>(Wyy),
                            static_cast<float>(Wyz)};
                        fwrite(shear_pi, sizeof(float), 5, out_file_xyeta);
                    }

                    if (DATA.turn_on_bulk) {
                        float bulk_pi[] = {static_cast<float>(pi_b)};
                        fwrite(bulk_pi, sizeof(float), 1, out_file_xyeta);
                    }

                    if (DATA.turn_on_diff) {
                        float diffusion[] = {
                            static_cast<float>(qx), static_cast<float>(qy),
                            static_cast<float>(qz)};
                        fwrite(diffusion, sizeof(float), 3, out_file_xyeta);
                    }
                }
                // rho_Q, mu_Q, rho_S, mu_S (4 floats, as counted in
                // nVar_per_cell; mu_Q = mu_S = 0 for EoS without them)
                if (DATA.turn_on_QS) {
                    float mu_QS[] = {
                        static_cast<float>(rhoq_local),
                        static_cast<float>(muQ_local * hbarc),
                        static_cast<float>(rhos_local),
                        static_cast<float>(muS_local * hbarc)};
                    fwrite(mu_QS, sizeof(float), 4, out_file_xyeta);
                }
            }
        }
    }
    fclose(out_file_xyeta);
    return 0;
}

//! This function outputs hydro evolution file in binary format
void Cell_info::OutputEvolutionDataXYEta_MLtraining(Fields &arena, double tau) {
    // the format of the file is as follows,
    //    e P ux uy ueta Wxx Wxy Wxeta Wyy Wyeta bulk_pi rhob
    const string out_name_xyeta = "evolution_all_xyeta_ML.dat";
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
    double output_dtau = DATA.delta_tau * n_skip_tau;

    int n_skip_x = DATA.output_evolution_every_N_x;
    int n_skip_y = DATA.output_evolution_every_N_y;
    int n_skip_eta = DATA.output_evolution_every_N_eta;

    // write out header
    const int output_nx = static_cast<int>(arena.nX() / n_skip_x);
    const int output_ny = static_cast<int>(arena.nY() / n_skip_y);
    const int output_neta = static_cast<int>(arena.nEta() / n_skip_eta);
    const double output_dx = DATA.delta_x * n_skip_x;
    const double output_dy = DATA.delta_y * n_skip_y;
    const double output_deta = DATA.delta_eta * n_skip_eta;
    const double output_xmin = -DATA.x_size / 2.;
    const double output_ymin = -DATA.y_size / 2.;
    const double output_etamin = DATA.etaMin();

    int nVar_per_cell = 13;

    if (tau == DATA.tau0) {
        float header[] = {static_cast<float>(DATA.tau0),
                          static_cast<float>(output_dtau),
                          static_cast<float>(output_nx),
                          static_cast<float>(output_dx),
                          static_cast<float>(output_xmin),
                          static_cast<float>(output_ny),
                          static_cast<float>(output_dy),
                          static_cast<float>(output_ymin),
                          static_cast<float>(output_neta),
                          static_cast<float>(output_deta),
                          static_cast<float>(output_etamin),
                          static_cast<float>(nVar_per_cell)};
        fwrite(header, sizeof(float), 12, out_file_xyeta);
    }
    ThermalVec thermalVec;
    for (int ieta = 0; ieta < arena.nEta(); ieta += n_skip_eta) {
        for (int iy = 0; iy < arena.nY(); iy += n_skip_y) {
            for (int ix = 0; ix < arena.nX(); ix += n_skip_x) {
                int fieldIdx = arena.getFieldIdx(ix, iy, ieta);
                double e_local = arena.e_[fieldIdx];        // 1/fm^4
                double rhob_local = arena.rhob_[fieldIdx];  // 1/fm^3

                eos.getThermalVariables(
                    e_local, rhob_local, arena.rhoq_[fieldIdx],
                    arena.rhos_[fieldIdx], thermalVec);
                double p_local = thermalVec[2];
                double T_local = thermalVec[6] * hbarc;  // GeV
                // double muB_local = thermalVec[7] * hbarc;  // GeV

                double ux = arena.u_[1][fieldIdx];
                double uy = arena.u_[2][fieldIdx];
                double ueta = arena.u_[3][fieldIdx];

                // e + P >= 0, = 0 only in empty cells
                double div_factor =
                    std::max(e_local + p_local, Util::tiny_normal);  // 1/fm^4
                double Wxx = arena.Wmunu_[4][fieldIdx] / div_factor;
                double Wxy = arena.Wmunu_[5][fieldIdx] / div_factor;
                double Wxeta = arena.Wmunu_[6][fieldIdx] / div_factor;
                double Wyy = arena.Wmunu_[7][fieldIdx] / div_factor;
                double Wyeta = arena.Wmunu_[8][fieldIdx] / div_factor;
                double pi_b = arena.piBulk_[fieldIdx] / div_factor;

                float fluidCell[] = {
                    static_cast<float>(e_local * hbarc),
                    static_cast<float>(p_local * hbarc),
                    static_cast<float>(ux),
                    static_cast<float>(uy),
                    static_cast<float>(ueta),
                    static_cast<float>(Wxx),
                    static_cast<float>(Wxy),
                    static_cast<float>(Wxeta),
                    static_cast<float>(Wyy),
                    static_cast<float>(Wyeta),
                    static_cast<float>(pi_b),
                    static_cast<float>(rhob_local),
                    static_cast<float>(T_local),
                };

                fwrite(fluidCell, sizeof(float), nVar_per_cell, out_file_xyeta);
            }
        }
    }
    fclose(out_file_xyeta);
}

//! This function outputs hydro evolution file in binary format for photon
//! production
void Cell_info::OutputEvolutionDataXYEta_photon(Fields &arena, double tau) {
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
    double volume = tau * n_skip_tau * dtau * n_skip_x * dx * n_skip_y * dy
                    * n_skip_eta * deta;

    for (int ieta = 0; ieta < arena.nEta(); ieta += n_skip_eta) {
        double eta_local = DATA.etaAt(ieta);
        for (int iy = 0; iy < arena.nY(); iy += n_skip_y) {
            for (int ix = 0; ix < arena.nX(); ix += n_skip_x) {
                int fieldIdx = arena.getFieldIdx(ix, iy, ieta);
                double e_local = arena.e_[fieldIdx];  // 1/fm^4

                if (e_local * hbarc < DATA.output_evolution_e_cut) continue;
                // only ouput fluid cells that are above cut-off energy density

                double rhob_local = arena.rhob_[fieldIdx];  // 1/fm^3
                double rhoq_local = arena.rhoq_[fieldIdx];  // 1/fm^3
                double rhos_local = arena.rhos_[fieldIdx];  // 1/fm^3

                double ux = arena.u_[1][fieldIdx];
                double uy = arena.u_[2][fieldIdx];
                double ueta = arena.u_[3][fieldIdx];

                // T_local is in 1/fm
                double T_local = eos.get_temperature(
                    e_local, rhob_local, rhoq_local, rhos_local);
                double p_local = eos.get_pressure(
                    e_local, rhob_local, rhoq_local, rhos_local);
                double muB_local = 0.0;
                double muQ_local = 0.0;
                double muS_local = 0.0;
                if (DATA.turn_on_rhob) {
                    muB_local = eos.get_muB(
                        e_local, rhob_local, rhoq_local, rhos_local);
                    if (DATA.turn_on_QS && DATA.whichEOS == 20) {
                        muQ_local = eos.get_muQ(
                            e_local, rhob_local, rhoq_local, rhos_local);
                        muS_local = eos.get_muS(
                            e_local, rhob_local, rhoq_local, rhos_local);
                    }
                }

                // e + P >= 0, = 0 only in empty cells
                double div_factor =
                    std::max(e_local + p_local, Util::tiny_normal);  // 1/fm^4
                double Wxx = 0.0;
                double Wxy = 0.0;
                double Wxeta = 0.0;
                double Wyy = 0.0;
                double Wyeta = 0.0;
                if (DATA.turn_on_shear) {
                    Wxx = arena.Wmunu_[4][fieldIdx] / div_factor;
                    Wxy = arena.Wmunu_[5][fieldIdx] / div_factor;
                    Wxeta = arena.Wmunu_[6][fieldIdx] / div_factor;
                    Wyy = arena.Wmunu_[7][fieldIdx] / div_factor;
                    Wyeta = arena.Wmunu_[8][fieldIdx] / div_factor;
                }

                double pi_b = 0.0;
                if (DATA.turn_on_bulk) {
                    pi_b = arena.piBulk_[fieldIdx];  // 1/fm^4
                }

                // outputs for baryon diffusion part
                double common_term_q = 0.0;
                double qx = 0.0;
                double qy = 0.0;
                double qeta = 0.0;
                if (DATA.turn_on_diff) {
                    common_term_q = rhob_local * T_local / div_factor;
                    double kappa_hat = get_deltaf_qmu_coeff(T_local, muB_local);
                    qx = arena.Wmunu_[11][fieldIdx] / kappa_hat;
                    qy = arena.Wmunu_[12][fieldIdx] / kappa_hat;
                    qeta = arena.Wmunu_[13][fieldIdx] / kappa_hat;
                }

                float ideal[] = {
                    static_cast<float>(volume),
                    static_cast<float>(eta_local),
                    static_cast<float>(T_local * hbarc),
                    static_cast<float>(ux),
                    static_cast<float>(uy),
                    static_cast<float>(ueta)};

                fwrite(ideal, sizeof(float), 6, out_file_xyeta);

                if (DATA.turn_on_rhob) {
                    float mu[] = {static_cast<float>(muB_local * hbarc)};
                    fwrite(mu, sizeof(float), 1, out_file_xyeta);
                }

                if (DATA.turn_on_shear) {
                    float shear_pi[] = {
                        static_cast<float>(Wxx), static_cast<float>(Wxy),
                        static_cast<float>(Wxeta), static_cast<float>(Wyy),
                        static_cast<float>(Wyeta)};
                    fwrite(shear_pi, sizeof(float), 5, out_file_xyeta);
                }

                if (DATA.turn_on_bulk) {
                    float bulk_pi[] = {static_cast<float>(pi_b)};
                    fwrite(bulk_pi, sizeof(float), 1, out_file_xyeta);
                }

                if (DATA.turn_on_diff) {
                    float diffusion[] = {
                        static_cast<float>(common_term_q),
                        static_cast<float>(qx), static_cast<float>(qy),
                        static_cast<float>(qeta)};
                    fwrite(diffusion, sizeof(float), 4, out_file_xyeta);
                }

                // rho_Q, mu_Q, rho_S, mu_S (mu_Q = mu_S = 0 for EoS without
                // them)
                if (DATA.turn_on_QS) {
                    float mu_QS[] = {
                        static_cast<float>(rhoq_local),
                        static_cast<float>(muQ_local * hbarc),
                        static_cast<float>(rhos_local),
                        static_cast<float>(muS_local * hbarc)};
                    fwrite(mu_QS, sizeof(float), 4, out_file_xyeta);
                }
            }
        }
    }
    fclose(out_file_xyeta);
}

//! This function outputs hydro evolution file in binary format
void Cell_info::OutputEvolutionDataXYEta_vorticity(
    Fields &arena_curr, Fields &arena_prev, double tau) {
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

    int n_skip_tau = DATA.output_evolution_every_N_timesteps;
    double output_dtau = DATA.delta_tau * n_skip_tau;
    int itau = static_cast<int>((tau - DATA.tau0) / (output_dtau) + 0.1);

    int n_skip_x = DATA.output_evolution_every_N_x;
    int n_skip_y = DATA.output_evolution_every_N_y;
    int n_skip_eta = DATA.output_evolution_every_N_eta;

    // write out header
    const int output_nx = static_cast<int>(arena_curr.nX() / n_skip_x);
    const int output_ny = static_cast<int>(arena_curr.nY() / n_skip_y);
    const int output_neta = static_cast<int>(arena_curr.nEta() / n_skip_eta);
    const double output_dx = DATA.delta_x * n_skip_x;
    const double output_dy = DATA.delta_y * n_skip_y;
    const double output_deta = DATA.delta_eta * n_skip_eta;
    const double output_xmin = -DATA.x_size / 2.;
    const double output_ymin = -DATA.y_size / 2.;
    const double output_etamin = DATA.etaMin();

    const int nVar_per_cell = 35;

    if (tau == DATA.tau0) {
        float header[] = {static_cast<float>(DATA.tau0),
                          static_cast<float>(output_dtau),
                          static_cast<float>(output_nx),
                          static_cast<float>(output_dx),
                          static_cast<float>(output_xmin),
                          static_cast<float>(output_ny),
                          static_cast<float>(output_dy),
                          static_cast<float>(output_ymin),
                          static_cast<float>(output_neta),
                          static_cast<float>(output_deta),
                          static_cast<float>(output_etamin),
                          static_cast<float>(nVar_per_cell)};
        fwrite(header, sizeof(float), 12, out_file_xyeta);
    }
    for (int ieta = 0; ieta < arena_curr.nEta(); ieta += n_skip_eta) {
        double eta_local = DATA.etaAt(ieta);
        for (int iy = 0; iy < arena_curr.nY(); iy += n_skip_y) {
            for (int ix = 0; ix < arena_curr.nX(); ix += n_skip_x) {
                int fieldIdx = arena_curr.getFieldIdx(ix, iy, ieta);
                double e_local = arena_curr.e_[fieldIdx];        // 1/fm^4
                double rhob_local = arena_curr.rhob_[fieldIdx];  // 1/fm^3
                double rhoq_local = arena_curr.rhoq_[fieldIdx];  // 1/fm^3
                double rhos_local = arena_curr.rhos_[fieldIdx];  // 1/fm^3

                double p_local = eos.get_pressure(
                    e_local, rhob_local, rhoq_local, rhos_local);
                // T_local, muB_local are in 1/fm
                double T_local = eos.get_temperature(
                    e_local, rhob_local, rhoq_local, rhos_local);
                double muB_local = 0.0;
                double muQ_local = 0.0;
                double muS_local = 0.0;
                muB_local =
                    eos.get_muB(e_local, rhob_local, rhoq_local, rhos_local);
                muQ_local =
                    eos.get_muQ(e_local, rhob_local, rhoq_local, rhos_local);
                muS_local =
                    eos.get_muS(e_local, rhob_local, rhoq_local, rhos_local);

                double ux = arena_curr.u_[1][fieldIdx];
                double uy = arena_curr.u_[2][fieldIdx];
                double ueta = arena_curr.u_[3][fieldIdx];

                if (e_local * hbarc < DATA.output_evolution_e_cut) continue;
                // only ouput fluid cells that are above cut-off energy density

                VorticityVec omega_kSP = {0.0};
                VorticityVec omega_k = {0.0};
                VorticityVec omega_th = {0.0};
                VorticityVec omega_T = {0.0};
                VelocityShearVec sigma_munu = {0.0};
                DmuMuBoverTVec DbetaMu = {0.0};
                u_derivative_helper.compute_vorticity_shell(
                    tau, arena_prev, arena_curr, ieta, ix, iy, fieldIdx,
                    eta_local, omega_kSP, omega_k, omega_th, omega_T,
                    sigma_munu, DbetaMu);

                float ideal[] = {
                    static_cast<float>(itau),
                    static_cast<float>(ix / n_skip_x),
                    static_cast<float>(iy / n_skip_y),
                    static_cast<float>(ieta / n_skip_eta),
                    static_cast<float>(e_local * hbarc),
                    static_cast<float>(p_local * hbarc),
                    static_cast<float>(T_local * hbarc),
                    static_cast<float>(ux),
                    static_cast<float>(uy),
                    static_cast<float>(ueta),
                    static_cast<float>(muB_local * hbarc),
                    static_cast<float>(muQ_local * hbarc),
                    static_cast<float>(muS_local * hbarc)};
                fwrite(ideal, sizeof(float), 13, out_file_xyeta);

                float vor_vec[6];
                for (int i = 0; i < 6; i++) {
                    // no minus sign because it has an opposite sign compared
                    // to kinetic vorticity
                    vor_vec[i] =
                        static_cast<float>(omega_kSP[i] * hbarc);  // GeV
                }
                fwrite(vor_vec, sizeof(float), 6, out_file_xyeta);
                for (int i = 0; i < 6; i++) {
                    // the minus sign is from metric
                    // output quantities for g = (1, -1, -1 , -1)
                    vor_vec[i] =
                        static_cast<float>(-omega_k[i] * hbarc);  // GeV
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
                    vor_vec[i] = static_cast<float>(
                        -omega_T[i] * hbarc * hbarc);  // GeV^2
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
    Fields &arena, double &e_max, double &nB_max, double &Tmax,
    const GridMaxima *known) {
    double eps_max = 0.0;
    double rhob_max = 0.0;
    double rhoq_max = -10.0;
    double rhos_max = -10.0;
    double T_max = 0.0;
    const bool sweep = (known == nullptr || !known->valid);
    if (!sweep) {
        eps_max = known->e;
        rhob_max = known->rhob;
        rhoq_max = known->rhoq;
        rhos_max = known->rhos;
        T_max = known->T;
    }

    // get the grid information
    const int neta = arena.nEta();
    const int nx = arena.nX();
    const int ny = arena.nY();

// every maximum is a reduction: shared updates of rhoq_max and rhos_max
// raced between the threads
#pragma omp parallel for collapse(3) if (sweep) \
    reduction(max : eps_max, rhob_max, rhoq_max, rhos_max, T_max)
    for (int ieta = 0; ieta < (sweep ? neta : 0); ieta++)
        for (int iy = 0; iy < ny; iy++)
            for (int ix = 0; ix < nx; ix++) {
                int fieldIdx = arena.getFieldIdx(ix, iy, ieta);
                const auto eps_local = arena.e_[fieldIdx];
                const auto rhob_local = arena.rhob_[fieldIdx];
                const auto rhoq_local = arena.rhoq_[fieldIdx];
                const auto rhos_local = arena.rhos_[fieldIdx];

                eps_max = std::max(eps_max, eps_local);
                rhob_max = std::max(rhob_max, rhob_local);
                T_max = std::max(
                    T_max, eos.get_temperature(
                               eps_local, rhob_local, rhoq_local, rhos_local));
                rhoq_max = std::max(rhoq_max, rhoq_local);
                rhos_max = std::max(rhos_max, rhos_local);
            }
    eps_max *= Util::hbarc;  // GeV/fm^3
    T_max *= Util::hbarc;    // GeV

    if (eps_max > 1e5) {
        music_message << "The maximum e = " << eps_max << " < 1e5 GeV/fm^3";
        music_message.flush("error");
        music_message.error("This normally should not happen!");
        music_message.error("Exiting ...");
        exit(1);
    }
    music_message << "eps_max = " << eps_max << " GeV/fm^3, "
                  << "rhob_max = " << rhob_max << " 1/fm^3, ";
    if (DATA.turn_on_QS) {
        music_message << "rhoq_max = " << rhoq_max << " 1/fm^3, "
                      << "rhos_max = " << rhos_max << " 1/fm^3, ";
    }
    music_message << "T_max = " << T_max << " GeV.";

    music_message.flush("info");
    e_max = eps_max;
    nB_max = rhob_max;
    Tmax = T_max;
}

//! This function computes OAM area density at a give proper time
void Cell_info::compute_Lmunu(
    Fields &arena, Fields &arena_prev, const double tau) {
    const double deta = DATA.delta_eta;
    const double dx = DATA.delta_x;
    const double dy = DATA.delta_y;
    const int neta = arena.nEta();
    const int nx = arena.nX();
    const int ny = arena.nY();

    ostringstream tau_suffix;
    tau_suffix << std::fixed << std::setfill('0') << std::setw(7)
               << std::setprecision(3) << tau;
    ostringstream filename;
    filename << "Lmunu_angular_momentum_tau_" << tau_suffix.str() << ".dat";

    ofstream output_file;
    output_file.open(filename.str().c_str(), std::ofstream::out);
    output_file << "# Tau = " << tau << "\n";
    output_file << "# deta = " << deta << "\n";
    output_file << "# dx = " << dx << "\n";
    output_file << "# dy = " << dy << "\n";
    output_file << "# neta = " << neta << "\n";
    output_file << "# nx = " << nx << "\n";
    output_file << "# ny = " << ny << "\n";

    output_file
        << "Lx_xy[hbar fm^-2]   Ly_xy[hbar fm^-2]   Lz_xy[hbar fm^-2]   "
        << "L^{tx}_xy[hbar fm^-2]   L^{ty}_xy[hbar fm^-2]   L^{tz}_xy[hbar "
           "fm^-2]"
        << std::endl;

    std::vector<std::vector<std::vector<double>>> Lx(
        neta,
        std::vector<std::vector<double>>(nx, std::vector<double>(ny, 0.0)));
    std::vector<std::vector<std::vector<double>>> Ly = Lx;
    std::vector<std::vector<std::vector<double>>> Lz = Lx;
    std::vector<std::vector<std::vector<double>>> Ltx = Lx;
    std::vector<std::vector<std::vector<double>>> Lty = Lx;
    std::vector<std::vector<std::vector<double>>> Ltz = Lx;

    std::vector<std::vector<double>> Lx_xy(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double>> Ly_xy(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double>> Lz_xy(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double>> Ltx_xy(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double>> Lty_xy(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double>> Ltz_xy(nx, std::vector<double>(ny, 0.0));

    for (int ieta = 0; ieta < neta; ieta++) {
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                const int Idx = arena.getFieldIdx(ix, iy, ieta);

                double eta_s = DATA.etaAt(ieta);
                const double cosh_eta = cosh(eta_s);
                const double sinh_eta = sinh(eta_s);
                const double t_local = tau * cosh_eta;
                const double x_local = -DATA.x_size / 2.0 + ix * dx;
                const double y_local = -DATA.y_size / 2.0 + iy * dy;
                const double z_local = tau * sinh_eta;

                const double e_local = arena.e_[Idx];
                const double pressure = eos.get_pressure(
                    e_local, arena.rhob_[Idx], arena.rhoq_[Idx],
                    arena.rhos_[Idx]);
                const double u0 = arena.u_[0][Idx];
                const double u1 = arena.u_[1][Idx];
                const double u2 = arena.u_[2][Idx];
                const double u3 = arena.u_[3][Idx];
                const double uPrev0 = arena_prev.u_[0][Idx];
                const double uPrev1 = arena_prev.u_[1][Idx];
                const double uPrev2 = arena_prev.u_[2][Idx];
                const double uPrev3 = arena_prev.u_[3][Idx];

                const double T00_local =
                    (e_local + pressure) * u0 * u0 - pressure;
                const double Pi00_rk_0 =
                    arena_prev.piBulk_[Idx] * (-1. + uPrev0 * uPrev0);

                const double T_tau_tau =
                    (T00_local + arena_prev.Wmunu_[0][Idx] + Pi00_rk_0);
                const double T_tau_x =
                    ((e_local + pressure) * u0 * u1 + arena_prev.Wmunu_[1][Idx]
                     + arena_prev.piBulk_[Idx] * uPrev0 * uPrev1);
                const double T_tau_y =
                    ((e_local + pressure) * u0 * u2 + arena_prev.Wmunu_[2][Idx]
                     + arena_prev.piBulk_[Idx] * uPrev0 * uPrev2);
                const double T_tau_eta =
                    ((e_local + pressure) * u0 * u3 + arena_prev.Wmunu_[3][Idx]
                     + arena_prev.piBulk_[Idx] * uPrev0 * uPrev3);
                const double T_tau_t =
                    T_tau_tau * cosh_eta + T_tau_eta * sinh_eta;
                const double T_tau_z =
                    T_tau_tau * sinh_eta + T_tau_eta * cosh_eta;

                Lx[ieta][ix][iy] = y_local * T_tau_z - z_local * T_tau_y;
                Ly[ieta][ix][iy] = z_local * T_tau_x - x_local * T_tau_z;
                Lz[ieta][ix][iy] = x_local * T_tau_y - y_local * T_tau_x;

                Ltx[ieta][ix][iy] = t_local * T_tau_x - x_local * T_tau_t;
                Lty[ieta][ix][iy] = t_local * T_tau_y - y_local * T_tau_t;
                Ltz[ieta][ix][iy] = t_local * T_tau_z - z_local * T_tau_t;
            }
        }
    }

    // Integration over eta
    // Lmunu area densities are written

    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int ieta = 0; ieta < neta; ++ieta) {
                const double dz = tau * deta;
                Lx_xy[ix][iy] += Lx[ieta][ix][iy] * dz;
                Ly_xy[ix][iy] += Ly[ieta][ix][iy] * dz;
                Lz_xy[ix][iy] += Lz[ieta][ix][iy] * dz;
                Ltx_xy[ix][iy] += Ltx[ieta][ix][iy] * dz;
                Lty_xy[ix][iy] += Lty[ieta][ix][iy] * dz;
                Ltz_xy[ix][iy] += Ltz[ieta][ix][iy] * dz;
            }

            output_file << std::scientific << std::setprecision(6)
                        << Lx_xy[ix][iy] << " " << Ly_xy[ix][iy] << " "
                        << Lz_xy[ix][iy] << " " << Ltx_xy[ix][iy] << " "
                        << Lty_xy[ix][iy] << " " << Ltz_xy[ix][iy] << std::endl;
        }
    }
    output_file.close();
}

//! This function computes global angular momentum at a give proper time
void Cell_info::compute_angular_momentum(
    Fields &arena, Fields &arena_prev, const double tau, const double eta_min,
    const double eta_max) {
    ostringstream filename;
    filename << "global_angular_momentum_eta_" << eta_min << "_" << eta_max
             << ".dat";
    ofstream output_file;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        // create new files at the first time step
        output_file.open(filename.str().c_str(), std::ofstream::out);
        output_file << "# tau[fm]  Lx[hbarc]  Ly[hbarc]  Lz[hbarc]  "
                    << "L^{tx}[hbarc]  L^{ty}[hbarc]  L^{tz}[hbarc]"
                    << std::endl;
    } else {
        output_file.open(
            filename.str().c_str(), std::fstream::out | std::fstream::app);
    }
    double Lx = 0.0;
    double Ly = 0.0;
    double Lz = 0.0;
    double Ltx = 0.0;
    double Lty = 0.0;
    double Ltz = 0.0;
    const double deta = DATA.delta_eta;
    const double dx = DATA.delta_x;
    const double dy = DATA.delta_y;
    const int neta = arena.nEta();
    const int nx = arena.nX();
    const int ny = arena.nY();
#pragma omp parallel for collapse(3) reduction(+ : Lx, Ly, Lz, Ltx, Lty, Ltz)
    for (int ieta = 0; ieta < neta; ieta++)
        for (int ix = 0; ix < nx; ix++)
            for (int iy = 0; iy < ny; iy++) {
                const int Idx = arena.getFieldIdx(ix, iy, ieta);

                double eta_s = DATA.etaAt(ieta);
                if (DATA.boost_invariant) {
                    eta_s = 0.0;
                }

                const double cosh_eta = cosh(eta_s);
                const double sinh_eta = sinh(eta_s);
                const double t_local = tau * cosh_eta;
                const double x_local = -DATA.x_size / 2. + ix * dx;
                const double y_local = -DATA.y_size / 2. + iy * dy;
                const double z_local = tau * sinh_eta;

                const double e_local = arena.e_[Idx];
                const double pressure = eos.get_pressure(
                    e_local, arena.rhob_[Idx], arena.rhoq_[Idx],
                    arena.rhos_[Idx]);

                const double u0 = arena.u_[0][Idx];
                const double u1 = arena.u_[1][Idx];
                const double u2 = arena.u_[2][Idx];
                const double u3 = arena.u_[3][Idx];
                const double uPrev0 = arena_prev.u_[0][Idx];
                const double uPrev1 = arena_prev.u_[1][Idx];
                const double uPrev2 = arena_prev.u_[2][Idx];
                const double uPrev3 = arena_prev.u_[3][Idx];

                // Adapt these formula to 4D EoS ?
                const double T00_local =
                    (e_local + pressure) * u0 * u0 - pressure;
                const double Pi00_rk_0 =
                    arena_prev.piBulk_[Idx] * (-1. + uPrev0 * uPrev0);

                const double T_tau_tau =
                    (T00_local + arena_prev.Wmunu_[0][Idx] + Pi00_rk_0);
                const double T_tau_x =
                    ((e_local + pressure) * u0 * u1 + arena_prev.Wmunu_[1][Idx]
                     + arena_prev.piBulk_[Idx] * uPrev0 * uPrev1);
                const double T_tau_y =
                    ((e_local + pressure) * u0 * u2 + arena_prev.Wmunu_[2][Idx]
                     + arena_prev.piBulk_[Idx] * uPrev0 * uPrev2);
                const double T_tau_eta =
                    ((e_local + pressure) * u0 * u3 + arena_prev.Wmunu_[3][Idx]
                     + arena_prev.piBulk_[Idx] * uPrev0 * uPrev3);
                const double T_tau_t =
                    T_tau_tau * cosh_eta + T_tau_eta * sinh_eta;
                const double T_tau_z =
                    T_tau_tau * sinh_eta + T_tau_eta * cosh_eta;

                if (eta_s < eta_max && eta_s > eta_min) {
                    Lx += (y_local * T_tau_z - z_local * T_tau_y);
                    Ly += (z_local * T_tau_x - x_local * T_tau_z);
                    Lz += (x_local * T_tau_y - y_local * T_tau_x);
                    Ltx += (t_local * T_tau_x - x_local * T_tau_t);
                    Lty += (t_local * T_tau_y - y_local * T_tau_t);
                    Ltz += (t_local * T_tau_z - z_local * T_tau_t);
                }
            }
    // add units
    double factor = tau * dx * dy * deta;
    Lx *= factor;
    Ly *= factor;
    Lz *= factor;
    Ltx *= factor;
    Lty *= factor;
    Ltz *= factor;

    // output results
    output_file << scientific << setprecision(6) << tau << "  " << Lx << "  "
                << Ly << "  " << Lz << "  " << Ltx << "  " << Lty << "  " << Ltz
                << std::endl;
    output_file.close();
}

//! This function checks the total energy and total net baryon, electric charge
//! and strangeness number at a give proper time
void Cell_info::check_conservation_law(
    Fields &arena, Fields &arena_prev, const double tau,
    const std::vector<double> *pressureIn, GridMaxima *maxima,
    const std::vector<double> *deviceRows,
    const std::vector<double> *outflowCurr,
    const std::vector<double> *outflowPrev) {
    std::string filename = "global_conservation_laws.dat";
    ofstream output_file;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        output_file.open(filename.c_str(), std::ofstream::out);
        output_file << "# tau(fm)  E(GeV)  Px(GeV)  Py(GeV)  Pz(GeV)  N_B  "
                    << "N_Q  N_S" << std::endl;
    } else {
        output_file.open(
            filename.c_str(), std::fstream::out | std::fstream::app);
    }
    double N_B = 0.0;
    double N_Q = 0.0;
    double N_S = 0.0;
    double T_tau_t = 0.0;
    double T_tau_x = 0.0;
    double T_tau_y = 0.0;
    double T_tau_z = 0.0;
    double deta = DATA.delta_eta;
    double dx = DATA.delta_x;
    double dy = DATA.delta_y;
    const int neta = arena.nEta();
    const int nx = arena.nX();
    const int ny = arena.nY();

    // cosh and sinh of eta_s, once per eta slice
    std::vector<double> coshEta(neta), sinhEta(neta);
    for (int ieta = 0; ieta < neta; ieta++) {
        const double eta_s = DATA.etaAt(ieta);
        coshEta[ieta] = cosh(eta_s);
        sinhEta[ieta] = sinh(eta_s);
    }
    // the sums of every (eta, x) row, added up in a fixed order below
    const int nSums = (deviceRows != nullptr ? 12 : 7);
    std::vector<double> rowSums;
    double eps_max = 0.0;
    double rhob_max = 0.0;
    double rhoq_max = -10.0;
    double rhos_max = -10.0;
    double T_max = 0.0;
    const bool withMaxima = (maxima != nullptr);
    const size_t nRows = static_cast<size_t>(neta) * nx;
    if (deviceRows != nullptr) {
        if (deviceRows->size() != nSums * nRows) {
            music_message.error(
                "check_conservation_law: the device rows do not match the "
                "grid");
            exit(1);
        }
        for (size_t row = 0; row < nRows; row++) {
            const double *r = deviceRows->data() + nSums * row;
            eps_max = std::max(eps_max, r[7]);
            rhob_max = std::max(rhob_max, r[8]);
            rhoq_max = std::max(rhoq_max, r[9]);
            rhos_max = std::max(rhos_max, r[10]);
            T_max = std::max(T_max, r[11]);
        }
    } else {
        rowSums.assign(static_cast<size_t>(nSums) * nRows, 0.);
    }
    const std::vector<double> &rows =
        (deviceRows != nullptr ? *deviceRows : rowSums);

#pragma omp parallel for collapse(2) schedule(static) \
    reduction(max : eps_max, rhob_max, rhoq_max, rhos_max, T_max)
    for (int ieta = 0; ieta < (deviceRows != nullptr ? 0 : neta); ieta++)
        for (int ix = 0; ix < nx; ix++) {
            const double cosh_eta = coshEta[ieta];
            const double sinh_eta = sinhEta[ieta];
            double rN_B = 0.0, rN_Q = 0.0, rN_S = 0.0;
            double rT_t = 0.0, rT_x = 0.0, rT_y = 0.0, rT_z = 0.0;
            for (int iy = 0; iy < ny; iy++) {
                const int Idx = arena.getFieldIdx(ix, iy, ieta);

                // The KT/RK2 update conserves tau_n T_ideal^{tau mu}(n) +
                // 1/2 [tau_{n+1} V(n) + tau_n V(n-1)], with V^mu = W^{tau mu}
                // + Pi (g^{tau mu} + u^tau u^mu) and Pi = Pi_bulk + Pi_chem
                // at each time level (see Advance::AdvanceIt). Divided by
                // tau here; the factor below multiplies by tau.
                const double wCurr = 0.5 * (1. + DATA.delta_tau / tau);
                const double wPrev = 0.5;
                const double e_local = arena.e_[Idx];
                const double pressure =
                    (pressureIn != nullptr
                         ? (*pressureIn)[Idx]
                         : eos.get_pressure(
                               e_local, arena.rhob_[Idx], arena.rhoq_[Idx],
                               arena.rhos_[Idx]));
                const double PiCurr =
                    arena.piBulk_[Idx] + arena.piBulkChem_[Idx];
                const double PiPrev =
                    arena_prev.piBulk_[Idx] + arena_prev.piBulkChem_[Idx];
                double T[4];
                for (int a = 0; a < 4; a++) {
                    const double g = (a == 0 ? -1. : 0.);
                    const double VCurr =
                        arena.Wmunu_[a][Idx]
                        + PiCurr * (g + arena.u_[0][Idx] * arena.u_[a][Idx]);
                    const double VPrev =
                        arena_prev.Wmunu_[a][Idx]
                        + PiPrev
                              * (g
                                 + arena_prev.u_[0][Idx]
                                       * arena_prev.u_[a][Idx]);
                    T[a] = (e_local + pressure) * arena.u_[0][Idx]
                               * arena.u_[a][Idx]
                           + (a == 0 ? -pressure : 0.) + wCurr * VCurr
                           + wPrev * VPrev;
                }
                rN_B += arena.rhob_[Idx] * arena.u_[0][Idx]
                        + wCurr * arena.Wmunu_[10][Idx]
                        + wPrev * arena_prev.Wmunu_[10][Idx];
                // no diffusion currents for Q and S
                rN_Q += arena.rhoq_[Idx] * arena.u_[0][Idx];
                rN_S += arena.rhos_[Idx] * arena.u_[0][Idx];
                const double T_tau_tau = T[0];
                const double T_tau_eta = T[3];
                rT_t += T_tau_tau * cosh_eta + T_tau_eta * sinh_eta;
                rT_x += T[1];
                rT_y += T[2];
                rT_z += T_tau_tau * sinh_eta + T_tau_eta * cosh_eta;

                if (withMaxima) {
                    const double rhob_local = arena.rhob_[Idx];
                    const double rhoq_local = arena.rhoq_[Idx];
                    const double rhos_local = arena.rhos_[Idx];
                    eps_max = std::max(eps_max, e_local);
                    rhob_max = std::max(rhob_max, rhob_local);
                    T_max = std::max(
                        T_max,
                        eos.get_temperature(
                            e_local, rhob_local, rhoq_local, rhos_local));
                    rhoq_max = std::max(rhoq_max, rhoq_local);
                    rhos_max = std::max(rhos_max, rhos_local);
                }
            }
            double *r =
                rowSums.data() + nSums * (ix + static_cast<size_t>(nx) * ieta);
            r[0] = rN_B;
            r[1] = rN_Q;
            r[2] = rN_S;
            r[3] = rT_t;
            r[4] = rT_x;
            r[5] = rT_y;
            r[6] = rT_z;
        }
    for (size_t row = 0; row < nRows; row++) {
        const double *r = rows.data() + nSums * row;
        N_B += r[0];
        N_Q += r[1];
        N_S += r[2];
        T_tau_t += r[3];
        T_tau_x += r[4];
        T_tau_y += r[5];
        T_tau_z += r[6];
    }
    if (withMaxima) {
        maxima->valid = true;
        maxima->e = eps_max;
        maxima->rhob = rhob_max;
        maxima->rhoq = rhoq_max;
        maxima->rhos = rhos_max;
        maxima->T = T_max;
    }
    // add units
    double factor = tau * dx * dy * deta;
    N_B *= factor;
    N_Q *= factor;
    N_S *= factor;
    T_tau_t *= factor * Util::hbarc;  // GeV
    T_tau_x *= factor * Util::hbarc;  // GeV
    T_tau_y *= factor * Util::hbarc;  // GeV
    T_tau_z *= factor * Util::hbarc;  // GeV
    // add what has left the grid: the time integral (trapezoidal rule) of
    // the flux through the grid faces since tau0
    // the outflow rate of arena (or arena_prev) at t
    const auto rateOf = [&](Fields &a, const std::vector<double> *faces,
                                 double t, TJbVec &rate) {
        if (faces != nullptr) {
            boundary_outflow_rate(*faces, a.nX(), a.nY(), a.nEta(), t, rate);
        } else {
            boundary_outflow_rate(a, t, rate);
        }
    };
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        outflow_flux = {0.};
        rateOf(arena, outflowCurr, tau, outflowRateLast_);
        outflowTauLast_ = tau;
    } else {
        TJbVec fluxCurr, fluxPrev;
        rateOf(arena, outflowCurr, tau, fluxCurr);
        // the grid was remapped at tau: the fluid left through the faces of
        // the old grid, whose rate at tau was recorded before the remap
        const bool regridded =
            (regridTau_ > 0. && std::abs(tau - regridTau_) < 1e-10);
        // the interval since the previous record: delta_tau, except in the
        // step where beastMode changed it (arena_prev is then extrapolated
        // to tau - delta_tau_new) or remapped the grid; use the rate stored
        // at that record
        double dt = DATA.delta_tau;
        if (outflowTauLast_ > 0.
            && (regridded
                || std::abs(tau - outflowTauLast_ - DATA.delta_tau) > 1e-8)) {
            dt = tau - outflowTauLast_;
            fluxPrev = outflowRateLast_;
        } else {
            rateOf(arena_prev, outflowPrev, tau - DATA.delta_tau, fluxPrev);
        }
        const TJbVec &fluxEnd = regridded ? outflowRateBeforeRegrid_ : fluxCurr;
        for (unsigned int i = 0; i < outflow_flux.size(); i++) {
            outflow_flux[i] += 0.5 * dt * (fluxEnd[i] + fluxPrev[i]);
        }
        outflowRateLast_ = fluxCurr;
        outflowTauLast_ = tau;
    }
    T_tau_t += outflow_flux[0];  // GeV
    T_tau_x += outflow_flux[1];  // GeV
    T_tau_y += outflow_flux[2];  // GeV
    T_tau_z += outflow_flux[3];  // GeV
    N_B += outflow_flux[4];
    N_Q += outflow_flux[5];
    N_S += outflow_flux[6];

    // output results
    music_message << "total energy T^{taut} = " << T_tau_t << " GeV";
    music_message.flush("info");
    music_message << "net longitudinal momentum Pz = " << T_tau_z << " GeV";
    music_message.flush("info");
    music_message << "net baryon number N_B = " << N_B;
    // N_B is exactly 0 without baryon transport, so 0 must pass
    if (N_B >= 0. && N_B < 500.) {
        music_message.flush("info");
    } else {
        music_message.flush("error");
        exit(1);
    }
    music_message << "net electric charges N_Q = " << N_Q
                  << ", net strangeness N_S = " << N_S;
    music_message.flush("info");
    output_file << scientific << setprecision(10) << tau << "  " << T_tau_t
                << "  " << T_tau_x << "  " << T_tau_y << "  " << T_tau_z << "  "
                << N_B << "  " << N_Q << "  " << N_S << std::endl;
    output_file.close();
}

void Cell_info::noteRegrid(Fields &arena, const double tau) {
    if (regridTau_ > 0. && std::abs(tau - regridTau_) < 1e-10) return;
    boundary_outflow_rate(arena, tau, outflowRateBeforeRegrid_);
    regridTau_ = tau;
}

//! Rate (per unit tau) at which E, Px, Py, Pz [GeV], N_B, N_Q and N_S leave
//! the grid
//! through its faces. With the components tau u^eta, tau pi^{mu eta}, ...
//! (hat), the Cartesian components nu = t, x, y, z obey
//!   d_tau(tau T^{tau nu}) + d_i(tau T^{i nu}) + d_eta(hat T^{eta nu}) = 0,
//! hat T^{eta t} = cosh(eta) hat T^{eta tau} + sinh(eta) hat T^{eta eta}.
//! The ghost cells copy the edge cells, so the flux through a face is the
//! flux of the edge cell.
void Cell_info::outflowFace(
    int j, int nx, int ny, int neta, int &ix, int &iy, int &ieta, int &mu,
    int &side) {
    // per eta slice: the x faces (iy, then -/+), then the y faces (ix,
    // then -/+); then, with neta > 1, the eta faces ((iy, ix), then -/+)
    const int perSlice = 2 * ny + 2 * nx;
    if (j < neta * perSlice) {
        ieta = j / perSlice;
        int r = j % perSlice;
        side = r % 2;
        if (r < 2 * ny) {
            mu = 1;
            iy = r / 2;
            ix = (side == 0 ? 0 : nx - 1);
        } else {
            r -= 2 * ny;
            mu = 2;
            ix = r / 2;
            iy = (side == 0 ? 0 : ny - 1);
        }
    } else {
        const int r = j - neta * perSlice;
        side = r % 2;
        mu = 3;
        iy = (r / 2) / nx;
        ix = (r / 2) % nx;
        ieta = (side == 0 ? 0 : neta - 1);
    }
}

int Cell_info::outflowFaceCount(int nx, int ny, int neta) {
    return (neta * (2 * ny + 2 * nx) + (neta > 1 ? 2 * nx * ny : 0));
}

template <typename Values>
void Cell_info::outflowRate(
    int nx, int ny, int neta, const double tau, TJbVec &rate,
    const Values &values) const {
    rate = {0.};
    const double dx = DATA.delta_x;
    const double dy = DATA.delta_y;
    const double deta = DATA.delta_eta;
    // W^{mu nu} index for (mu, nu) in tau, x, y, eta
    const int widx[4][4] = {
        {0, 1, 2, 3}, {1, 4, 5, 6}, {2, 5, 7, 8}, {3, 6, 8, 9}};
    const double ghat[4] = {1., -1., -1., -1.};

    // adds the outward flux of the cell with the values v (kOutflowVals:
    // e, rhob, rhoq, rhos, Pi, Pi_chem, u^mu, W^0 ... W^13) through the
    // face normal to direction mu (1 = x, 2 = y, 3 = eta) with sign = +1 / -1
    auto addFace = [&](const double *v, double eta, int mu, double sign,
                       double area) {
        const double e = v[0];
        const double rhob = v[1];
        const double rhoq = v[2];
        const double rhos = v[3];
        const double P = eos.get_pressure(e, rhob, rhoq, rhos);
        // total bulk pressure, as in MakeWSource and the cell content
        const double Pi = v[4] + v[5];
        double u[4];
        for (int i = 0; i < 4; i++) u[i] = v[6 + i];
        double T[4];  // hat T^{mu nu}, nu = tau, x, y, eta
        for (int nu = 0; nu < 4; nu++) {
            T[nu] = (e + P + Pi) * u[mu] * u[nu] + v[10 + widx[mu][nu]];
            if (nu == mu) T[nu] -= (P + Pi) * ghat[mu];
        }
        const double J = rhob * u[mu] + v[10 + 10 + mu];
        const double ch = cosh(eta);
        const double sh = sinh(eta);
        const double f = sign * area;
        rate[0] += f * (ch * T[0] + sh * T[3]) * Util::hbarc;
        rate[1] += f * T[1] * Util::hbarc;
        rate[2] += f * T[2] * Util::hbarc;
        rate[3] += f * (sh * T[0] + ch * T[3]) * Util::hbarc;
        rate[4] += f * J;
        rate[5] += f * rhoq * u[mu];
        rate[6] += f * rhos * u[mu];
    };

    // The eta faces sit half a cell beyond the outermost cell centers.
    const double eta_m = DATA.etaMin() - deta / 2.0;
    const double eta_p = DATA.etaMax() + deta / 2.0;
    const int nFaces = outflowFaceCount(nx, ny, neta);
    for (int j = 0; j < nFaces; j++) {
        int ix, iy, ieta, mu, side;
        outflowFace(j, nx, ny, neta, ix, iy, ieta, mu, side);
        const double sign = (side == 0 ? -1. : +1.);
        const double *v = values(j, ix, iy, ieta);
        if (mu == 3) {
            addFace(v, side == 0 ? eta_m : eta_p, 3, sign, dx * dy);
        } else {
            const double eta = DATA.etaAt(ieta);
            addFace(v, eta, mu, sign,
                    (mu == 1 ? tau * dy * deta : tau * dx * deta));
        }
    }
}

void Cell_info::boundary_outflow_rate(
    Fields &arena, const double tau, TJbVec &rate) const {
    double v[kOutflowVals];
    outflowRate(
        arena.nX(), arena.nY(), arena.nEta(), tau, rate,
        [&](int, int ix, int iy, int ieta) {
            const int Idx = arena.getFieldIdx(ix, iy, ieta);
            v[0] = arena.e_[Idx];
            v[1] = arena.rhob_[Idx];
            v[2] = arena.rhoq_[Idx];
            v[3] = arena.rhos_[Idx];
            v[4] = arena.piBulk_[Idx];
            v[5] = arena.piBulkChem_[Idx];
            for (int i = 0; i < 4; i++) v[6 + i] = arena.u_[i][Idx];
            for (int i = 0; i < 14; i++) v[10 + i] = arena.Wmunu_[i][Idx];
            return (static_cast<const double *>(v));
        });
}

void Cell_info::boundary_outflow_rate(
    const std::vector<double> &faces, int nx, int ny, int neta,
    const double tau, TJbVec &rate) const {
    if (faces.size()
        != static_cast<size_t>(kOutflowVals) * outflowFaceCount(nx, ny, neta)) {
        std::cerr << "boundary_outflow_rate: the face values do not match "
                  << "the grid" << std::endl;
        exit(1);
    }
    outflowRate(nx, ny, neta, tau, rate, [&](int j, int, int, int) {
        return (faces.data() + static_cast<size_t>(kOutflowVals) * j);
    });
}

//! This function putputs files to check with Gubser flow solution
void Cell_info::Gubser_flow_check_file(Fields &arena, const double tau) {
    if (tau > 1.) {
        ostringstream filename_analytic;
        filename_analytic << "tests/Gubser_flow/y=x_tau=" << tau
                          << "_SemiAnalytic.dat";

        double T_analytic[201], ux_analytic[201], uy_analytic[201];
        double pixx_analytic[201], pixy_analytic[201];
        double piyy_analytic[201], pizz_analytic[201];
        double dummy;
        std::ifstream input_file(filename_analytic.str().c_str());
        for (int i = 0; i < 201; i++) {
            input_file >> dummy >> dummy >> T_analytic[i] >> ux_analytic[i]
                >> uy_analytic[i] >> pixx_analytic[i] >> piyy_analytic[i]
                >> pixy_analytic[i] >> pizz_analytic[i];
        }
        input_file.close();

        double T_diff = 0.0;
        double ux_diff = 0.0;
        double uy_diff = 0.0;
        double pixx_diff = 0.0;
        double pixy_diff = 0.0;
        double piyy_diff = 0.0;
        double pizz_diff = 0.0;
        double T_sum = 0.0;
        double ux_sum = 0.0;
        double uy_sum = 0.0;
        double pixx_sum = 0.0;
        double pixy_sum = 0.0;
        double piyy_sum = 0.0;
        double pizz_sum = 0.0;
        for (int i = 0; i < arena.nX(); i++) {
            int fieldIdx = arena.getFieldIdx(i, i, 0);
            double e_local = arena.e_[fieldIdx];

            double T_local =
                (eos.get_temperature(e_local, 0.0, 0.0, 0.0) * Util::hbarc);

            T_diff += fabs(T_analytic[i] - T_local);
            T_sum += fabs(T_analytic[i]);
            ux_diff += fabs(ux_analytic[i] - arena.u_[1][fieldIdx]);
            ux_sum += fabs(ux_analytic[i]);
            uy_diff += fabs(uy_analytic[i] - arena.u_[2][fieldIdx]);
            uy_sum += fabs(uy_analytic[i]);
            pixx_diff += (fabs(
                pixx_analytic[i] - arena.Wmunu_[4][fieldIdx] * Util::hbarc));
            pixx_sum += fabs(pixx_analytic[i]);
            pixy_diff += (fabs(
                pixy_analytic[i] - arena.Wmunu_[5][fieldIdx] * Util::hbarc));
            pixy_sum += fabs(pixy_analytic[i]);
            piyy_diff += (fabs(
                piyy_analytic[i] - arena.Wmunu_[7][fieldIdx] * Util::hbarc));
            piyy_sum += fabs(piyy_analytic[i]);
            pizz_diff += (fabs(
                pizz_analytic[i] - arena.Wmunu_[9][fieldIdx] * Util::hbarc));
            pizz_sum += fabs(pizz_analytic[i]);
        }
        music_message << "Autocheck: T_diff = " << T_diff / std::max(T_sum, small_eps)
                      << ", ux_diff = " << ux_diff / std::max(ux_sum, small_eps)
                      << ", uy_diff = " << uy_diff / std::max(uy_sum, small_eps)
                      << ", pixx_diff = " << pixx_diff / std::max(pixx_sum, small_eps)
                      << ", pixy_diff = " << pixy_diff / std::max(pixy_sum, small_eps)
                      << ", piyy_diff = " << piyy_diff / std::max(piyy_sum, small_eps)
                      << ", pizz_diff = " << pizz_diff / std::max(pizz_sum, small_eps);
        music_message.flush("info");
    }

    ostringstream filename;
    filename << "Gubser_flow_check_tau_" << tau << ".dat";
    ofstream output_file(filename.str().c_str());

    double dx = DATA.delta_x;
    double x_min = -DATA.x_size / 2.;
    double dy = DATA.delta_y;
    double y_min = -DATA.y_size / 2.;
    for (int ix = 0; ix < arena.nX(); ix++)
        for (int iy = 0; iy < arena.nY(); iy++) {
            double x_local = x_min + ix * dx;
            double y_local = y_min + iy * dy;
            int fieldIdx = arena.getFieldIdx(ix, iy, 0);
            double e_local = arena.e_[fieldIdx];

            double rhob_local = arena.rhob_[fieldIdx];
            double rhoq_local = arena.rhoq_[fieldIdx];
            double rhos_local = arena.rhos_[fieldIdx];

            double T_local = eos.get_temperature(e_local, 0.0);
            output_file << scientific << setprecision(8) << setw(18) << x_local
                        << "  " << y_local << "  " << e_local * Util::hbarc
                        << "  " << rhob_local << "  " << T_local * Util::hbarc
                        << "  " << arena.u_[1][fieldIdx] << "  "
                        << arena.u_[2][fieldIdx] << "  "
                        << arena.Wmunu_[4][fieldIdx] * Util::hbarc << "  "
                        << arena.Wmunu_[7][fieldIdx] * Util::hbarc << "  "
                        << arena.Wmunu_[5][fieldIdx] * Util::hbarc << "  "
                        << arena.Wmunu_[9][fieldIdx] * Util::hbarc << "  "
                        << rhoq_local << "  " << rhos_local << "  " << endl;
        }
    output_file.close();
}

//! This function outputs files to cross check with 1+1D simulation
void Cell_info::output_1p1D_check_file(Fields &arena, const double tau) {
    ostringstream filename;
    filename << "1+1D_check_tau_" << tau << ".dat";
    ofstream output_file(filename.str().c_str());

    for (int ieta = 0; ieta < arena.nEta(); ieta++) {
        const double eta_local = DATA.etaAt(ieta);
        int fieldIdx = arena.getFieldIdx(1, 1, ieta);
        double e_local = arena.e_[fieldIdx];

        double rhob_local = arena.rhob_[fieldIdx];
        double rhoq_local = arena.rhoq_[fieldIdx];
        double rhos_local = arena.rhos_[fieldIdx];

        output_file << scientific << setprecision(8) << setw(18) << eta_local
                    << "  " << e_local * Util::hbarc << "  " << rhob_local
                    << "  " << rhoq_local << "  " << rhos_local << "  " << endl;
    }
    output_file.close();
}

//! This function outputs energy density and n_b, n_q and n_s for making movies
void Cell_info::output_evolution_for_movie(Fields &arena, const double tau) {
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
    double output_dtau = DATA.delta_tau * n_skip_tau;
    int itau = static_cast<int>((tau - DATA.tau0) / (output_dtau) + 0.1);

    int n_skip_x = DATA.output_evolution_every_N_x;
    int n_skip_y = DATA.output_evolution_every_N_y;
    int n_skip_eta = DATA.output_evolution_every_N_eta;
    double dx = DATA.delta_x;
    double dy = DATA.delta_y;
    double deta = DATA.delta_eta;
    double volume =
        tau * output_dtau * n_skip_x * dx * n_skip_y * dy * n_skip_eta * deta;

    const int nVar_per_cell = 14;
    if (tau == DATA.tau0) {
        // write out header
        const int output_nx = static_cast<int>(arena.nX() / n_skip_x);
        const int output_ny = static_cast<int>(arena.nY() / n_skip_y);
        const int output_neta =
            (std::max(1, static_cast<int>(arena.nEta() / n_skip_eta)));
        const double output_dx = DATA.delta_x * n_skip_x;
        const double output_dy = DATA.delta_y * n_skip_y;
        const double output_deta = DATA.delta_eta * n_skip_eta;
        const double output_xmin = -DATA.x_size / 2.;
        const double output_ymin = -DATA.y_size / 2.;
        const double output_etamin = DATA.etaMin();
        float header[] = {static_cast<float>(DATA.tau0),
                          static_cast<float>(output_dtau),
                          static_cast<float>(output_nx),
                          static_cast<float>(output_dx),
                          static_cast<float>(output_xmin),
                          static_cast<float>(output_ny),
                          static_cast<float>(output_dy),
                          static_cast<float>(output_ymin),
                          static_cast<float>(output_neta),
                          static_cast<float>(output_deta),
                          static_cast<float>(output_etamin),
                          static_cast<float>(nVar_per_cell)};
        fwrite(header, sizeof(float), 12, out_file_xyeta);
    }
    for (int ieta = 0; ieta < arena.nEta(); ieta += n_skip_eta) {
        double eta_local = DATA.etaAt(ieta);
        if (DATA.boost_invariant) eta_local = 0.;
        for (int iy = 0; iy < arena.nY(); iy += n_skip_y) {
            for (int ix = 0; ix < arena.nX(); ix += n_skip_x) {
                int fieldIdx = arena.getFieldIdx(ix, iy, ieta);
                double e_local = arena.e_[fieldIdx];        // 1/fm^4
                double rhob_local = arena.rhob_[fieldIdx];  // 1/fm^3
                double rhoq_local = arena.rhoq_[fieldIdx];  // 1/fm^3
                double rhos_local = arena.rhos_[fieldIdx];  // 1/fm^3

                double T_local = eos.get_temperature(
                    e_local, rhob_local, rhoq_local, rhos_local);
                if (T_local * hbarc < DATA.output_evolution_T_cut) continue;
                double pressure = eos.get_pressure(
                    e_local, rhob_local, rhoq_local, rhos_local);
                double muB_local = 0.0;
                double muQ_local = 0.0;
                double muS_local = 0.0;
                muB_local = eos.get_muB(
                    e_local, rhob_local,  // 1/fm
                    rhoq_local, rhos_local);
                muQ_local = eos.get_muQ(
                    e_local, rhob_local,  // 1/fm
                    rhoq_local, rhos_local);
                muS_local = eos.get_muS(
                    e_local, rhob_local,  // 1/fm
                    rhoq_local, rhos_local);

                double u0 = arena.u_[0][fieldIdx];
                double u1 = arena.u_[1][fieldIdx];
                double u2 = arena.u_[2][fieldIdx];
                double u3 = arena.u_[3][fieldIdx];
                double T00_ideal = (e_local + pressure) * u0 * u0 - pressure;
                double T03_ideal = (e_local + pressure) * u0 * u3;
                double Pi00 = arena.piBulk_[fieldIdx] * (-1.0 + u0 * u0);
                double Pi03 = arena.piBulk_[fieldIdx] * u0 * u3;
                double T00_full =
                    (T00_ideal + arena.Wmunu_[0][fieldIdx] + Pi00);
                double T03_full =
                    (T03_ideal + arena.Wmunu_[3][fieldIdx] + Pi03);
                double Ttaut =
                    (T00_full * cosh(eta_local) + T03_full * sinh(eta_local));
                double JBtau = (rhob_local * u0 + arena.Wmunu_[10][fieldIdx]);
                float array[] = {
                    static_cast<float>(itau),
                    static_cast<float>(ix / n_skip_x),
                    static_cast<float>(iy / n_skip_y),
                    static_cast<float>(ieta / n_skip_eta),
                    static_cast<float>(volume),
                    static_cast<float>(e_local * hbarc),
                    static_cast<float>(rhob_local),
                    static_cast<float>(T_local * hbarc),
                    static_cast<float>(muB_local * hbarc),
                    static_cast<float>(u1),
                    static_cast<float>(u2),
                    static_cast<float>(u3),
                    static_cast<float>(Ttaut * hbarc),
                    static_cast<float>(JBtau),
                    static_cast<float>(rhoq_local),
                    static_cast<float>(rhos_local),
                    static_cast<float>(muQ_local * hbarc),
                    static_cast<float>(muS_local * hbarc)};
                fwrite(array, sizeof(float), 18, out_file_xyeta);
            }
        }
    }
    fclose(out_file_xyeta);
}

//! This function dumps the energy density and net baryon density
void Cell_info::output_energy_density_and_rhob_disitrubtion(
    Fields &arena, string filename) {
    ofstream output_file(filename.c_str());
    const int n_skip_x = DATA.output_evolution_every_N_x;
    const int n_skip_y = DATA.output_evolution_every_N_y;
    const int n_skip_eta = DATA.output_evolution_every_N_eta;
    for (int ieta = 0; ieta < arena.nEta(); ieta += n_skip_eta)
        for (int ix = 0; ix < arena.nX(); ix += n_skip_x)
            for (int iy = 0; iy < arena.nY(); iy += n_skip_y) {
                int fieldIdx = arena.getFieldIdx(ix, iy, ieta);
                double e_local = arena.e_[fieldIdx] * Util::hbarc;
                double rhob_local = arena.rhob_[fieldIdx];
                output_file << scientific << setprecision(5) << setw(18)
                            << e_local << "  " << rhob_local << endl;
            }
    output_file.close();
}

//! This function dumps the energy density and net densities
void Cell_info::output_energy_density_and_netdens_disitrubtion(
    Fields &arena, string filename) {
    ofstream output_file(filename.c_str());
    const int n_skip_x = DATA.output_evolution_every_N_x;
    const int n_skip_y = DATA.output_evolution_every_N_y;
    const int n_skip_eta = DATA.output_evolution_every_N_eta;
    for (int ieta = 0; ieta < arena.nEta(); ieta += n_skip_eta)
        for (int ix = 0; ix < arena.nX(); ix += n_skip_x)
            for (int iy = 0; iy < arena.nY(); iy += n_skip_y) {
                int fieldIdx = arena.getFieldIdx(ix, iy, ieta);
                double e_local = arena.e_[fieldIdx] * Util::hbarc;
                double rhob_local = arena.rhob_[fieldIdx];
                double rhoq_local = arena.rhoq_[fieldIdx];
                double rhos_local = arena.rhos_[fieldIdx];
                output_file << scientific << setprecision(5) << setw(18)
                            << e_local << "  " << rhob_local << "  "
                            << rhoq_local << "  " << rhos_local << endl;
            }
    output_file.close();
}

//! This function outputs the evolution of hydrodynamic variables at a
//! give fluid cell
void Cell_info::monitor_a_fluid_cell(
    Fields &arena_curr, Fields &arena_prev, const int ix, const int iy,
    const int ieta, const double tau) {
    int fieldIdx = arena_curr.getFieldIdx(ix, iy, ieta);
    ostringstream filename;
    filename << "monitor_fluid_cell_ix_" << ix << "_iy_" << iy << "_ieta_"
             << ieta << ".dat";
    ofstream output_file;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        output_file.open(filename.str().c_str(), std::ofstream::out);
        output_file << "# tau(fm)  e(GeV/fm^3)  rhob(1/fm^3)  "
                    << "u^tau  u^x  u^y  tau*u^eta  "
                    << "Bulk_Pi(GeV/fm^3)  theta(1/fm)  "
                    << "pi^{\\mu\\nu}(GeV/fm^3)  sigma^{\\mu\\nu}(1/fm)  "
                    << "omega^{\\mu\\nu}(1/fm)  "
                    << "rhoq(1/fm^3) rhos(1/fm^3)" << std::endl;
    } else {
        output_file.open(
            filename.str().c_str(), std::fstream::out | std::fstream::app);
    }
    output_file << scientific << setprecision(8) << tau << "  "
                << arena_curr.e_[fieldIdx] * Util::hbarc << "  "
                << arena_curr.rhob_[fieldIdx] << "  ";
    for (int i = 0; i < 4; i++) {
        output_file << scientific << setprecision(8)
                    << arena_curr.u_[i][fieldIdx] << "  ";
    }

    u_derivative_helper.MakedU(
        tau, arena_prev, arena_curr, fieldIdx, ix, iy, ieta);
    auto theta_local =
        u_derivative_helper.calculate_expansion_rate(tau, arena_curr, fieldIdx);
    output_file << scientific << setprecision(8)
                << arena_curr.piBulk_[fieldIdx] * Util::hbarc << "  "
                << theta_local << "  ";
    DumuVec a_local;
    u_derivative_helper.calculate_Du_supmu(tau, arena_curr, fieldIdx, a_local);
    VelocityShearVec sigma_local;
    u_derivative_helper.calculate_velocity_shear_tensor(
        tau, arena_curr, fieldIdx, theta_local, a_local, sigma_local);
    VorticityVec omega_local;
    u_derivative_helper.calculate_kinetic_vorticity_with_spatial_projector(
        tau, arena_curr, fieldIdx, a_local, omega_local);
    for (int i = 0; i < 10; i++) {
        output_file << scientific << setprecision(8)
                    << arena_curr.Wmunu_[i][fieldIdx] * Util::hbarc << "  "
                    << sigma_local[i] << "  ";
    }
    for (int i = 0; i < 6; i++) {
        output_file << scientific << setprecision(8) << omega_local[i] << "  ";
    }
    output_file << arena_curr.rhoq_[fieldIdx] << "  "
                << arena_curr.rhos_[fieldIdx] << "  ";
    output_file << endl;
    output_file.close();
}

void Cell_info::output_vorticity_distribution(
    Fields &arena_curr, Fields &arena_prev, const double tau,
    const double eta_min, const double eta_max) {
    // This function outputs the vorticity tensor at a given tau
    ostringstream filename1;
    filename1 << "vorticity_dis_kinetic_wSP_eta_" << eta_min << "_" << eta_max
              << "_tau_" << tau << ".dat";
    std::fstream of1;
    of1.open(filename1.str().c_str(), std::fstream::out);
    // write the header
    of1 << "# x[fm]  y[fm]  T[GeV]  muB[GeV]  "
        << "omega^{tx}/T  omega^{ty}/T  "
        << "omega^{tz}/T  omega^{xy}/T  omega^{xz}/T  "
        << "omega^{yz}/T" << std::endl;
    ostringstream filename2;
    filename2 << "vorticity_dis_kinetic_eta_" << eta_min << "_" << eta_max
              << "_tau_" << tau << ".dat";
    std::fstream of2;
    of2.open(filename2.str().c_str(), std::fstream::out);
    // write the header
    of2 << "# x[fm]  y[fm]  T[GeV]  muB[GeV]  "
        << "omega^{tx}/T  omega^{ty}/T  "
        << "omega^{tz}/T  omega^{xy}/T  omega^{xz}/T  "
        << "omega^{yz}/T" << std::endl;
    ostringstream filename3;
    filename3 << "vorticity_dis_thermal_eta_" << eta_min << "_" << eta_max
              << "_tau_" << tau << ".dat";
    std::fstream of3;
    of3.open(filename3.str().c_str(), std::fstream::out);
    // write the header
    of3 << "# x[fm]  y[fm]  T[GeV]  muB[GeV]  "
        << "omega^{tx}  omega^{ty}  "
        << "omega^{tz}  omega^{xy}  omega^{xz}  "
        << "omega^{yz}" << std::endl;
    ostringstream filename4;
    filename4 << "vorticity_dis_T_eta_" << eta_min << "_" << eta_max << "_tau_"
              << tau << ".dat";
    std::fstream of4;
    of4.open(filename4.str().c_str(), std::fstream::out);
    // write the header
    of4 << "# x[fm]  y[fm]  T[GeV]  muB[GeV]  "
        << "omega^{tx}/T^2  omega^{ty}/T^2  "
        << "omega^{tz}/T^2  omega^{xy}/T^2  omega^{xz}/T^2  "
        << "omega^{yz}/T^2" << std::endl;

    for (int ix = 0; ix < arena_curr.nX(); ix++) {
        for (int iy = 0; iy < arena_curr.nY(); iy++) {
            const double x_local = -DATA.x_size / 2. + ix * DATA.delta_x;
            const double y_local = -DATA.y_size / 2. + iy * DATA.delta_y;

            VorticityVec omega_kSP = {0.0};
            VorticityVec omega_k = {0.0};
            VorticityVec omega_th = {0.0};
            VorticityVec omega_T = {0.0};
            double T_avg = 0.0;
            double muB_avg = 0.0;
            double weight = 0.0;
            for (int ieta = 0; ieta < arena_curr.nEta(); ieta++) {
                int fieldIdx = arena_curr.getFieldIdx(ix, iy, ieta);
                double eta_local = DATA.etaAt(ieta);
                if (DATA.boost_invariant) eta_local = 0.0;
                if (eta_local < eta_max && eta_local > eta_min) {
                    const double e_local = arena_curr.e_[fieldIdx];
                    if (e_local < 0.1) continue;

                    const double rhob_local = arena_curr.rhob_[fieldIdx];
                    const double rhoq_local = arena_curr.rhoq_[fieldIdx];
                    const double rhos_local = arena_curr.rhos_[fieldIdx];
                    const double T_local = eos.get_temperature(
                        e_local, rhob_local, rhoq_local, rhos_local);
                    const double muB_local = eos.get_muB(
                        e_local, rhob_local, rhoq_local, rhos_local);

                    T_avg += e_local * T_local * hbarc;
                    muB_avg += e_local * muB_local * hbarc;

                    VorticityVec omega_local_1, omega_local_2;
                    VorticityVec omega_local_3, omega_local_4;
                    VelocityShearVec sigma_local = {0.0};
                    DmuMuBoverTVec DbetaMu = {0.0};
                    u_derivative_helper.compute_vorticity_shell(
                        tau, arena_prev, arena_curr, ieta, ix, iy, fieldIdx,
                        eta_local, omega_local_1, omega_local_2, omega_local_3,
                        omega_local_4, sigma_local, DbetaMu);
                    for (unsigned int ii = 0; ii < omega_k.size(); ii++) {
                        omega_kSP[ii] += e_local * omega_local_1[ii] / T_local;
                        omega_k[ii] += e_local * omega_local_2[ii] / T_local;
                        omega_th[ii] += e_local * omega_local_3[ii];
                        omega_T[ii] +=
                            (e_local * omega_local_4[ii] / T_local / T_local);
                    }
                    weight += e_local;
                }
            }
            weight = std::max(weight, small_eps);
            of1 << scientific << setprecision(8) << setw(18) << x_local << "  "
                << y_local << "  " << T_avg / weight << "  " << muB_avg / weight
                << "  ";
            of2 << scientific << setprecision(8) << setw(18) << x_local << "  "
                << y_local << "  " << T_avg / weight << "  " << muB_avg / weight
                << "  ";
            of3 << scientific << setprecision(8) << setw(18) << x_local << "  "
                << y_local << "  " << T_avg / weight << "  " << muB_avg / weight
                << "  ";
            of4 << scientific << setprecision(8) << setw(18) << x_local << "  "
                << y_local << "  " << T_avg / weight << "  " << muB_avg / weight
                << "  ";
            for (unsigned int ii = 0; ii < omega_k.size(); ii++) {
                // no minus sign because it has an opposite sign to
                // the kinetic vorticity
                of1 << scientific << setprecision(8) << setw(18)
                    << omega_kSP[ii] / weight << "  ";
                // minus sign from the metric
                // output quantities in g = (1, -1, -1 , -1)
                of2 << scientific << setprecision(8) << setw(18)
                    << -omega_k[ii] / weight << "  ";
                of3 << scientific << setprecision(8) << setw(18)
                    << -omega_th[ii] / weight << "  ";
                of4 << scientific << setprecision(8) << setw(18)
                    << -omega_T[ii] / weight << "  ";
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
    Fields &arena_curr, Fields &arena_prev, const double tau,
    const double eta_min, const double eta_max) {
    // This function outputs the time evolution of the vorticity tensor
    ostringstream filename1;
    filename1 << "vorticity_evo_kinetic_wSP_eta_" << eta_min << "_" << eta_max
              << ".dat";
    std::fstream of1;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        of1.open(filename1.str().c_str(), std::fstream::out);
        of1 << "# tau[fm]  omega^{tx}/T  omega^{ty}/T  "
            << "omega^{tz}/T  omega^{xy}/T  omega^{xz}/T  "
            << "omega^{yz}/T" << std::endl;
    } else {
        of1.open(
            filename1.str().c_str(), std::fstream::out | std::fstream::app);
    }
    ostringstream filename2;
    filename2 << "vorticity_evo_kinetic_eta_" << eta_min << "_" << eta_max
              << ".dat";
    std::fstream of2;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        of2.open(filename2.str().c_str(), std::fstream::out);
        of2 << "# tau[fm]  omega^{tx}/T  omega^{ty}/T  "
            << "omega^{tz}/T  omega^{xy}/T  omega^{xz}/T  "
            << "omega^{yz}/T" << std::endl;
    } else {
        of2.open(
            filename2.str().c_str(), std::fstream::out | std::fstream::app);
    }
    ostringstream filename3;
    filename3 << "vorticity_evo_thermal_eta_" << eta_min << "_" << eta_max
              << ".dat";
    std::fstream of3;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        of3.open(filename3.str().c_str(), std::fstream::out);
        of3 << "# tau[fm]  omega^{tx}  omega^{ty}  "
            << "omega^{tz}  omega^{xy}  omega^{xz}  "
            << "omega^{yz}" << std::endl;
    } else {
        of3.open(
            filename3.str().c_str(), std::fstream::out | std::fstream::app);
    }
    ostringstream filename4;
    filename4 << "vorticity_evo_T_eta_" << eta_min << "_" << eta_max << ".dat";
    std::fstream of4;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        of4.open(filename4.str().c_str(), std::fstream::out);
        of4 << "# tau[fm]  omega^{tx}/T^2  omega^{ty}/T^2  "
            << "omega^{tz}/T^2  omega^{xy}/T^2  "
            << "omega^{xz}/T^2  omega^{yz}/T^2" << std::endl;
    } else {
        of4.open(
            filename4.str().c_str(), std::fstream::out | std::fstream::app);
    }

    VorticityVec omega_kSP = {0.0};
    VorticityVec omega_k = {0.0};
    VorticityVec omega_th = {0.0};
    VorticityVec omega_T = {0.0};
    double weight = 0.0;
    for (int ieta = 0; ieta < arena_curr.nEta(); ieta++) {
        double eta = 0.0;
        if (!DATA.boost_invariant) {
            eta = DATA.etaAt(ieta);
        }
        if (eta < eta_max && eta > eta_min) {
            for (int iy = 0; iy < arena_curr.nY(); iy++)
                for (int ix = 0; ix < arena_curr.nX(); ix++) {
                    int fieldIdx = arena_curr.getFieldIdx(ix, iy, ieta);
                    const double e_local = arena_curr.e_[fieldIdx];  // 1/fm^4
                    if (e_local < 0.1) continue;

                    const double rhob_local =
                        arena_curr.rhob_[fieldIdx];  // 1/fm^3
                    const double rhoq_local =
                        arena_curr.rhoq_[fieldIdx];  // 1/fm^3
                    const double rhos_local =
                        arena_curr.rhos_[fieldIdx];  // 1/fm^3
                    const double T_local = eos.get_temperature(
                        e_local, rhob_local, rhoq_local, rhos_local);
                    VorticityVec omega_local_1, omega_local_2;
                    VorticityVec omega_local_3, omega_local_4;
                    VelocityShearVec sigma_local = {0.0};
                    DmuMuBoverTVec DbetaMu = {0.0};
                    u_derivative_helper.compute_vorticity_shell(
                        tau, arena_prev, arena_curr, ieta, ix, iy, fieldIdx,
                        eta, omega_local_1, omega_local_2, omega_local_3,
                        omega_local_4, sigma_local, DbetaMu);
                    for (unsigned int ii = 0; ii < omega_k.size(); ii++) {
                        omega_kSP[ii] += e_local * omega_local_1[ii] / T_local;
                        omega_k[ii] += e_local * omega_local_2[ii] / T_local;
                        omega_th[ii] += e_local * omega_local_3[ii];
                        omega_T[ii] +=
                            e_local * omega_local_4[ii] / T_local / T_local;
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
            << omega_kSP[ii] / weight << "  ";
        // minus sign from the metric
        // output quantities in g = (1, -1, -1 , -1)
        of2 << scientific << setprecision(8) << setw(18)
            << -omega_k[ii] / weight << "  ";
        of3 << scientific << setprecision(8) << setw(18)
            << -omega_th[ii] / weight << "  ";
        of4 << scientific << setprecision(8) << setw(18)
            << -omega_T[ii] / weight << "  ";
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
    deltaf_qmu_coeff_tb = new double *[deltaf_qmu_coeff_table_length_T];
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

    deltaf_coeff_tb_14mom_DPi = new double *[deltaf_coeff_table_14mom_length_T];
    deltaf_coeff_tb_14mom_BPi = new double *[deltaf_coeff_table_14mom_length_T];
    deltaf_coeff_tb_14mom_BPitilde =
        new double *[deltaf_coeff_table_14mom_length_T];
    deltaf_coeff_tb_14mom_DV = new double *[deltaf_coeff_table_14mom_length_T];
    deltaf_coeff_tb_14mom_BV = new double *[deltaf_coeff_table_14mom_length_T];
    deltaf_coeff_tb_14mom_Bpi_shear =
        new double *[deltaf_coeff_table_14mom_length_T];
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
    double hbarc3 = hbarc * hbarc * hbarc;
    double hbarc4 = hbarc3 * hbarc;
    for (int i = 0; i < deltaf_coeff_table_14mom_length_T; i++) {
        for (int j = 0; j < deltaf_coeff_table_14mom_length_mu; j++) {
            deltaf_coeff_tb_14mom_DPi[i][j] =
                deltaf_coeff_tb_14mom_DPi[i][j] * hbarc4;  // fm^4/GeV
            deltaf_coeff_tb_14mom_BPi[i][j] =
                deltaf_coeff_tb_14mom_BPi[i][j] * hbarc4;  // fm^4/(GeV^2)
            deltaf_coeff_tb_14mom_BPitilde[i][j] =
                deltaf_coeff_tb_14mom_BPitilde[i][j] * hbarc4;  // fm^4/(GeV^2)
            deltaf_coeff_tb_14mom_DV[i][j] =
                deltaf_coeff_tb_14mom_DV[i][j] * hbarc3;  // fm^3/GeV
            deltaf_coeff_tb_14mom_BV[i][j] =
                deltaf_coeff_tb_14mom_BV[i][j] * hbarc3;  // fm^3/(GeV^2)
            deltaf_coeff_tb_14mom_Bpi_shear[i][j] =
                deltaf_coeff_tb_14mom_Bpi_shear[i][j] * hbarc4;  // fm^4/(GeV^2)
        }
    }
}

double Cell_info::get_deltaf_qmu_coeff(double T, double muB) {
    if (muB < 0) {
        muB = -muB;
    }
    int idx_T = static_cast<int>(
        (T - delta_qmu_coeff_table_T0) / delta_qmu_coeff_table_dT);
    int idx_mu = static_cast<int>(
        (muB - delta_qmu_coeff_table_mu0) / delta_qmu_coeff_table_dmu);

    if (idx_T > deltaf_qmu_coeff_table_length_T - 2 || idx_T < 0) return (1.0);
    if (idx_mu > deltaf_qmu_coeff_table_length_mu - 2) return (1.0);

    double x_fraction =
        ((T - delta_qmu_coeff_table_T0) / delta_qmu_coeff_table_dT - idx_T);
    double y_fraction =
        ((muB - delta_qmu_coeff_table_mu0) / delta_qmu_coeff_table_dmu
         - idx_mu);

    double f1 = deltaf_qmu_coeff_tb[idx_T][idx_mu];
    double f2 = deltaf_qmu_coeff_tb[idx_T][idx_mu + 1];
    double f3 = deltaf_qmu_coeff_tb[idx_T + 1][idx_mu + 1];
    double f4 = deltaf_qmu_coeff_tb[idx_T + 1][idx_mu];

    double coeff =
        (f1 * (1. - x_fraction) * (1. - y_fraction)
         + f2 * (1. - x_fraction) * y_fraction + f3 * x_fraction * y_fraction
         + f4 * x_fraction * (1. - y_fraction));
    return (coeff);
}

double Cell_info::get_deltaf_coeff_14moments(
    double T, double muB, double type) {
    int idx_T = static_cast<int>(
        (T - delta_coeff_table_14mom_T0) / delta_coeff_table_14mom_dT);
    int idx_mu = static_cast<int>(
        (muB - delta_coeff_table_14mom_mu0) / delta_coeff_table_14mom_dmu);
    double x_fraction =
        ((T - delta_coeff_table_14mom_T0) / delta_coeff_table_14mom_dT - idx_T);
    double y_fraction =
        ((muB - delta_coeff_table_14mom_mu0) / delta_coeff_table_14mom_dmu
         - idx_mu);

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
    double f2 = deltaf_table[idx_T][idx_mu + 1];
    double f3 = deltaf_table[idx_T + 1][idx_mu + 1];
    double f4 = deltaf_table[idx_T + 1][idx_mu];

    double coeff =
        (f1 * (1. - x_fraction) * (1. - y_fraction)
         + f2 * (1. - x_fraction) * y_fraction + f3 * x_fraction * y_fraction
         + f4 * x_fraction * (1. - y_fraction));
    return (coeff);
}

//! This function outputs average T and mu_B as a function of proper tau
//! within a given space-time rapidity range
void Cell_info::output_average_phase_diagram_trajectory(
    const double tau, const double eta_min, const double eta_max,
    Fields &arena) {
    ostringstream filename;
    filename << "averaged_phase_diagram_trajectory_eta_" << eta_min << "_"
             << eta_max << ".dat";
    std::fstream of;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        of.open(filename.str().c_str(), std::fstream::out);
        of << "# tau(fm)  <T>(GeV)  std(T)(GeV)  <mu_B>(GeV)  std(mu_B)(GeV)  "
           << "V4 (fm^4)" << endl;
    } else {
        of.open(filename.str().c_str(), std::fstream::out | std::fstream::app);
    }
    double avg_T = 0.0;
    double avg_muB = 0.0;
    double avg_muQ = 0.0;
    double avg_muS = 0.0;
    double std_T = 0.0;
    double std_muB = 0.0;
    double std_muQ = 0.0;
    double std_muS = 0.0;
    double weight = 0.0;
    double V4 = 0.0;
    const double unit_volume =
        tau * DATA.delta_x * DATA.delta_y * DATA.delta_eta;
    for (int ieta = 0; ieta < arena.nEta(); ieta++) {
        double eta = 0.0;
        if (!DATA.boost_invariant) {
            eta = DATA.etaAt(ieta);
        }
        if (eta < eta_max && eta > eta_min) {
            double cosh_eta = cosh(eta);
            double sinh_eta = sinh(eta);
            for (int iy = 0; iy < arena.nY(); iy++)
                for (int ix = 0; ix < arena.nX(); ix++) {
                    int fieldIdx = arena.getFieldIdx(ix, iy, ieta);

                    double e_local = arena.e_[fieldIdx];  // 1/fm^4
                    if (e_local * hbarc > DATA.output_evolution_e_cut)
                        V4 += unit_volume;
                    double rhob_local = arena.rhob_[fieldIdx];  // 1/fm^3
                    double rhoq_local = arena.rhoq_[fieldIdx];  // 1/fm^3
                    double rhos_local = arena.rhos_[fieldIdx];  // 1/fm^3
                                                                //
                    double utau = arena.u_[0][fieldIdx];
                    double ueta = arena.u_[3][fieldIdx];
                    double ut =
                        utau * cosh_eta + ueta * sinh_eta;  // gamma factor
                                                            //
                    double T_local = eos.get_temperature(
                        e_local, rhob_local, rhoq_local, rhos_local);
                    double muB_local = eos.get_muB(
                        e_local, rhob_local, rhoq_local, rhos_local);
                    double muQ_local = eos.get_muQ(
                        e_local, rhob_local, rhoq_local, rhos_local);
                    double muS_local = eos.get_muS(
                        e_local, rhob_local, rhoq_local, rhos_local);

                    double weight_local = e_local * ut;
                    avg_T += T_local * weight_local;
                    avg_muB += muB_local * weight_local;
                    avg_muQ += muQ_local * weight_local;
                    avg_muS += muS_local * weight_local;
                    std_T += T_local * T_local * weight_local;
                    std_muB += muB_local * muB_local * weight_local;
                    std_muQ += muQ_local * muQ_local * weight_local;
                    std_muS += muS_local * muS_local * weight_local;
                    weight += weight_local;
                }
        }
    }
    avg_T = avg_T / std::max(weight, small_eps) * hbarc;
    avg_muB = avg_muB / std::max(weight, small_eps) * hbarc;
    avg_muQ = avg_muQ / std::max(weight, small_eps) * hbarc;
    avg_muS = avg_muS / std::max(weight, small_eps) * hbarc;
    std_T = sqrt(
        std_T / std::max(weight, small_eps) * hbarc * hbarc - avg_T * avg_T);
    std_muB = sqrt(
        std_muB / std::max(weight, small_eps) * hbarc * hbarc
        - avg_muB * avg_muB);
    std_muQ = sqrt(
        std_muQ / std::max(weight, small_eps) * hbarc * hbarc
        - avg_muQ * avg_muQ);
    std_muS = sqrt(
        std_muS / std::max(weight, small_eps) * hbarc * hbarc
        - avg_muS * avg_muS);
    of << scientific << setw(18) << setprecision(8) << tau << "  " << avg_T
       << "  " << std_T << "  " << avg_muB << "  " << std_muB << "  " << V4
       << " " << avg_muQ << "  " << std_muQ << "  " << "  " << avg_muS << "  "
       << std_muS << "  " << endl;
    of.close();
}

//! This function outputs system's eccentricity and momentum anisotropy
//! as functions of eta_s
void Cell_info::output_momentum_anisotropy_vs_etas(
    const int itau, const double tau, Fields &arena) const {
    ostringstream filename;
    if (itau == 0) {
        filename << "momentum_anisotropy_tau_INITIAL.dat";
    } else {
        filename << "momentum_anisotropy_tau_" << tau << ".dat";
    }
    std::fstream of;
    of.open(filename.str().c_str(), std::fstream::out);
    of << "# eta_s  epsilon_p(ideal)(cos)  epsilon_p(ideal)(sin)  "
       << "epsilon_p(shear)(cos)  epsilon_p(shear)(sin)  "
       << "epsilon_p(full)(cos)  epsilon_p(full)(sin)  " << endl;

    ostringstream filename1;
    if (itau == 0) {
        filename1 << "eccentricities_evo_ed_tau_INITIAL.dat";
    } else {
        filename1 << "eccentricities_evo_ed_tau_" << tau << ".dat";
    }
    std::fstream of1;
    of1.open(filename1.str().c_str(), std::fstream::out);
    of1 << "# eta_s  ed(GeV/fm^3)  ecc_n(cos)  ecc_n(sin) (n=1-6)  "
        << "x_o(ed)(fm)  y_o(ed)(fm)" << endl;

    ostringstream filename2;
    if (itau == 0) {
        filename2 << "eccentricities_evo_nB_tau_INITIAL.dat";
    } else {
        filename2 << "eccentricities_evo_nB_tau_" << tau << ".dat";
    }
    std::fstream of2;
    of2.open(filename2.str().c_str(), std::fstream::out);
    of2 << "# eta_s  nB(1/fm^3)  ecc_n(cos)  ecc_n(sin) (n=1-6)  "
        << "x_o(nB)(fm)  y_o(nB)(fm)" << endl;

    ostringstream filename2Q;
    filename2Q << "eccentricities_evo_nQ_tau_" << tau << ".dat";
    std::fstream of2Q;
    of2Q.open(filename2Q.str().c_str(), std::fstream::out);
    of2Q << "# eta_s  nQ(1/fm^3)  ecc_n(cos)  ecc_n(sin) (n=1-6)" << endl;

    ostringstream filename3;
    if (itau == 0) {
        filename3 << "meanpT_estimators_tau_INITIAL.dat";
    } else {
        filename3 << "meanpT_estimators_tau_" << tau << ".dat";
    }
    std::fstream of3;
    of3.open(filename3.str().c_str(), std::fstream::out);
    of3 << "# eta_s  dS/deta_s  dE/deta_s (GeV) [s] (1/fm^-3)  [r^2] (fm^2)"
        << endl;

    const int norder = 6;
    // each eta slice is summed on its own, so the slices run in parallel;
    // their lines are collected here and written in order below
    const int neta = arena.nEta();
    std::vector<std::array<std::ostringstream, 5>> lines(neta);
#pragma omp parallel for schedule(dynamic)
    for (int ieta = 0; ieta < neta; ieta++) {
        ThermalVec thermalVec;
        double eta = 0.0;
        if (!DATA.boost_invariant) {
            eta = DATA.etaAt(ieta);
        }

        // compute the central of mass position
        double x_ed_o = 0.0, x_nB_o = 0.0, x_nQ_o = 0.0;
        double y_ed_o = 0.0, y_nB_o = 0.0, y_nQ_o = 0.0;
        double w_ed_sum = 0.0, w_nB_sum = 0.0, w_nQ_sum = 0.0;
        for (int iy = 0; iy < arena.nY(); iy++)
            for (int ix = 0; ix < arena.nX(); ix++) {
                int fieldIdx = arena.getFieldIdx(ix, iy, ieta);
                double x_local = -DATA.x_size / 2. + ix * DATA.delta_x;
                double y_local = -DATA.y_size / 2. + iy * DATA.delta_y;
                double e_local = arena.e_[fieldIdx];      // 1/fm^4
                double nB_local = arena.rhob_[fieldIdx];  // 1/fm^3
                double nQ_local = arena.rhoq_[fieldIdx];  // 1/fm^3
                double gamma_perp = arena.u_[0][fieldIdx];
                x_ed_o += x_local * e_local * gamma_perp;
                y_ed_o += y_local * e_local * gamma_perp;
                w_ed_sum += e_local * gamma_perp;
                x_nB_o += x_local * nB_local * gamma_perp;
                y_nB_o += y_local * nB_local * gamma_perp;
                w_nB_sum += nB_local * gamma_perp;
                x_nQ_o += x_local * nQ_local * gamma_perp;
                y_nQ_o += y_local * nQ_local * gamma_perp;
                w_nQ_sum += nQ_local * gamma_perp;
            }
        x_ed_o /= std::max(small_eps, w_ed_sum);
        y_ed_o /= std::max(small_eps, w_ed_sum);
        x_nB_o /= std::max(small_eps, w_nB_sum);
        y_nB_o /= std::max(small_eps, w_nB_sum);
        x_nQ_o /= std::max(small_eps, w_nQ_sum);
        y_nQ_o /= std::max(small_eps, w_nQ_sum);

        // compute epsilon_p with ideal, ideal + shear, and full T^{\munu}
        std::vector<double> ep_num1(3, 0.0);
        std::vector<double> ep_num2(3, 0.0);
        std::vector<double> ep_den(3, 0.0);

        // spatial eccentricity arrays
        std::vector<double> eccn_ed_num1(norder, 0.0);
        std::vector<double> eccn_ed_num2(norder, 0.0);
        std::vector<double> eccn_ed_den(norder, 0.0);
        std::vector<double> eccn_nB_num1(norder, 0.0);
        std::vector<double> eccn_nB_num2(norder, 0.0);
        std::vector<double> eccn_nB_den(norder, 0.0);
        std::vector<double> eccn_nQ_num1(norder, 0.0);
        std::vector<double> eccn_nQ_num2(norder, 0.0);
        std::vector<double> eccn_nQ_den(norder, 0.0);

        std::vector<double> meanpT_est_num(4, 0.0);
        std::vector<double> meanpT_est_den(1, 0.0);
        for (int iy = 0; iy < arena.nY(); iy++)
            for (int ix = 0; ix < arena.nX(); ix++) {
                int fieldIdx = arena.getFieldIdx(ix, iy, ieta);
                double x_ed = -DATA.x_size / 2. + ix * DATA.delta_x - x_ed_o;
                double y_ed = -DATA.y_size / 2. + iy * DATA.delta_y - y_ed_o;
                double x_nB = -DATA.x_size / 2. + ix * DATA.delta_x - x_nB_o;
                double y_nB = -DATA.y_size / 2. + iy * DATA.delta_y - y_nB_o;
                double x_nQ = -DATA.x_size / 2. + ix * DATA.delta_x - x_nQ_o;
                double y_nQ = -DATA.y_size / 2. + iy * DATA.delta_y - y_nQ_o;
                double r_ed = sqrt(x_ed * x_ed + y_ed * y_ed);
                double r_nB = sqrt(x_nB * x_nB + y_nB * y_nB);
                double r_nQ = sqrt(x_nQ * x_nQ + y_nQ * y_nQ);
                double phi_ed = atan2(y_ed, x_ed);
                double phi_nB = atan2(y_nB, x_nB);
                double phi_nQ = atan2(y_nQ, x_nQ);

                double e_local = arena.e_[fieldIdx];        // 1/fm^4
                double rhob_local = arena.rhob_[fieldIdx];  // 1/fm^3
                double rhoq_local = arena.rhoq_[fieldIdx];  // 1/fm^3
                double rhos_local = arena.rhos_[fieldIdx];  // 1/fm^3

                eos.getThermalVariables(
                    e_local, rhob_local, rhoq_local, rhos_local, thermalVec);

                double P_local = thermalVec[2];
                double enthopy = e_local + P_local;
                double s_local = thermalVec[12];
                double u0 = arena.u_[0][fieldIdx];
                double ux = arena.u_[1][fieldIdx];
                double uy = arena.u_[2][fieldIdx];
                double pi_xx = arena.Wmunu_[4][fieldIdx];
                double pi_xy = arena.Wmunu_[5][fieldIdx];
                double pi_yy = arena.Wmunu_[7][fieldIdx];
                double bulk_Pi = arena.piBulk_[fieldIdx];

                double T_00_ideal = enthopy * u0 * u0 - P_local;
                double T_xx_ideal = enthopy * ux * ux + P_local;
                double T_xy_ideal = enthopy * ux * uy;
                double T_yy_ideal = enthopy * uy * uy + P_local;

                double T_xx_shear = T_xx_ideal + pi_xx;
                double T_xy_shear = T_xy_ideal + pi_xy;
                double T_yy_shear = T_yy_ideal + pi_yy;

                double T_xx_full = T_xx_shear - bulk_Pi * (-1 - ux * ux);
                double T_xy_full = T_xy_shear + bulk_Pi * ux * uy;
                double T_yy_full = T_yy_shear - bulk_Pi * (-1 - uy * uy);

                ep_num1[0] += T_xx_ideal - T_yy_ideal;
                ep_num2[0] += 2. * T_xy_ideal;
                ep_den[0] += T_xx_ideal + T_yy_ideal;
                ep_num1[1] += T_xx_shear - T_yy_shear;
                ep_num2[1] += 2. * T_xy_shear;
                ep_den[1] += T_xx_shear + T_yy_shear;
                ep_num1[2] += T_xx_full - T_yy_full;
                ep_num2[2] += 2. * T_xy_full;
                ep_den[2] += T_xx_full + T_yy_full;

                double w_ed = 0.0, w_nB = 0.0, w_nQ = 0.0;
                for (int i = 0; i < norder; i++) {
                    if (i == 1) {
                        w_ed = u0 * e_local * pow(r_ed, 3);
                        w_nB = u0 * rhob_local * pow(r_nB, 3);
                        w_nQ = u0 * rhoq_local * pow(r_nQ, 3);
                    } else {
                        w_ed = u0 * e_local * pow(r_ed, i);
                        w_nB = u0 * rhob_local * pow(r_nB, i);
                        w_nQ = u0 * rhoq_local * pow(r_nQ, i);
                    }
                    if (i == 0) {
                        eccn_ed_num1[i] += e_local;
                        eccn_nB_num1[i] += rhob_local;
                        eccn_nQ_num1[i] += rhoq_local;
                    } else {
                        eccn_ed_num1[i] += w_ed * cos(i * phi_ed);
                        eccn_ed_num2[i] += w_ed * sin(i * phi_ed);
                        eccn_ed_den[i] += w_ed;
                        eccn_nB_num1[i] += w_nB * cos(i * phi_nB);
                        eccn_nB_num2[i] += w_nB * sin(i * phi_nB);
                        eccn_nB_den[i] += w_nB;
                        eccn_nQ_num1[i] += w_nQ * cos(i * phi_nQ);
                        eccn_nQ_num2[i] += w_nQ * sin(i * phi_nQ);
                        eccn_nQ_den[i] += w_nQ;
                    }
                }

                meanpT_est_num[0] += tau * s_local * u0;          // dS/deta_s
                meanpT_est_num[1] += tau * T_00_ideal;            // dE/deta_s
                meanpT_est_num[2] += s_local * u0 * e_local;      // [s]
                meanpT_est_num[3] += r_ed * r_ed * u0 * e_local;  // [r^2]
                meanpT_est_den[0] += u0 * e_local;
            }

        // output results
        lines[ieta][0] << scientific << setw(18) << setprecision(8) << eta
                       << "  ";
        for (int i = 0; i < 3; i++) {
            lines[ieta][0] << ep_num1[i] / std::max(ep_den[i], small_eps)
                           << "  "
                           << ep_num2[i] / std::max(ep_den[i], small_eps)
                           << "  ";
        }
        lines[ieta][0] << endl;
        lines[ieta][1] << scientific << setw(18) << setprecision(8) << eta
                       << "  "
                       << eccn_ed_num1[0] * DATA.delta_x * DATA.delta_y * hbarc
                       << "  ";
        for (int i = 1; i < norder; i++) {
            // the minus sign ensure the vector points to the short axis
            lines[ieta][1] << -eccn_ed_num1[i]
                                  / std::max(eccn_ed_den[i], small_eps)
                           << "  "
                           << -eccn_ed_num2[i]
                                  / std::max(eccn_ed_den[i], small_eps)
                           << "  ";
        }
        lines[ieta][1] << x_ed_o << "  " << y_ed_o;
        lines[ieta][1] << endl;
        lines[ieta][2] << scientific << setw(18) << setprecision(8) << eta
                       << "  " << eccn_nB_num1[0] * DATA.delta_x * DATA.delta_y
                       << "  ";
        for (int i = 1; i < norder; i++) {
            // the minus sign ensure the vector points to the short axis
            lines[ieta][2] << -eccn_nB_num1[i]
                                  / std::max(eccn_nB_den[i], small_eps)
                           << "  "
                           << -eccn_nB_num2[i]
                                  / std::max(eccn_nB_den[i], small_eps)
                           << "  ";
        }
        lines[ieta][2] << x_nB_o << "  " << y_nB_o;
        lines[ieta][2] << endl;

        lines[ieta][3] << scientific << setw(18) << setprecision(8) << eta
                       << "  " << eccn_nQ_num1[0] * DATA.delta_x * DATA.delta_y
                       << "  ";
        for (int i = 1; i < norder; i++) {
            // the minus sign ensure the vector points to the short axis
            lines[ieta][3] << -eccn_nQ_num1[i]
                                  / std::max(eccn_nQ_den[i], small_eps)
                           << "  "
                           << -eccn_nQ_num2[i]
                                  / std::max(eccn_nQ_den[i], small_eps)
                           << "  ";
        }
        lines[ieta][3] << endl;

        lines[ieta][4] << scientific << setw(18) << setprecision(8) << eta
                       << "  "
                       << meanpT_est_num[0] * DATA.delta_x * DATA.delta_y
                       << "  "
                       << meanpT_est_num[1] * hbarc * DATA.delta_x
                              * DATA.delta_y
                       << "  " << meanpT_est_num[2] / std::max(meanpT_est_den[0], small_eps) << "  "
                       << meanpT_est_num[3] / std::max(meanpT_est_den[0], small_eps) << endl;
    }
    for (int ieta = 0; ieta < neta; ieta++) {
        of << lines[ieta][0].str();
        of1 << lines[ieta][1].str();
        of2 << lines[ieta][2].str();
        of2Q << lines[ieta][3].str();
        of3 << lines[ieta][4].str();
    }
    of.close();
    of1.close();
    of2.close();
    of2Q.close();
    of3.close();
}

//! This function outputs system's momentum anisotropy as a function of tau
void Cell_info::output_momentum_anisotropy_vs_tau(
    const double tau, const double eta_min, const double eta_max,
    Fields &arena) const {
    ostringstream filename;
    filename << "momentum_anisotropy_eta_" << eta_min << "_" << eta_max
             << ".dat";
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
           << "epsilon_3p(full)(cos)  epsilon_3p(full)(sin)  " << endl;
    } else {
        of.open(filename.str().c_str(), std::fstream::out | std::fstream::app);
    }

    ostringstream filename1;
    filename1 << "eccentricities_evo_eta_" << eta_min << "_" << eta_max
              << ".dat";
    std::fstream of1;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        of1.open(filename1.str().c_str(), std::fstream::out);
        of1 << "# tau(fm)  ecc_n(cos)  ecc_n(sin) (n=1-6) x_o (fm)  y_o (fm)"
            << endl;
    } else {
        of1.open(
            filename1.str().c_str(), std::fstream::out | std::fstream::app);
    }

    ostringstream filename2;
    filename2 << "inverse_Reynolds_number_eta_" << eta_min << "_" << eta_max
              << ".dat";
    std::fstream of2;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        of2.open(filename2.str().c_str(), std::fstream::out);
        of2 << "# tau(fm)  R_shearpi  R_Pi  gamma  T[GeV]" << endl;
    } else {
        of2.open(
            filename2.str().c_str(), std::fstream::out | std::fstream::app);
    }

    ostringstream filename3;
    filename3 << "meanpT_estimators_eta_" << eta_min << "_" << eta_max
              << ".dat";
    std::fstream of3;
    if (std::abs(tau - DATA.tau0) < 1e-10) {
        of3.open(filename3.str().c_str(), std::fstream::out);
        of3 << "# tau (fm)  dS/deta_s  dE/deta_s (GeV) [s] (1/fm^-3)  "
            << "[r^2] (fm^2)" << endl;
    } else {
        of3.open(
            filename3.str().c_str(), std::fstream::out | std::fstream::app);
    }

    // The sums run over the rows (eta, y) of the eta window in parallel:
    // each row is summed over x on its own, and the row sums are added in
    // row order, so the result does not depend on the number of threads.
    const int norder = 6;
    enum {
        kIdealNum1,
        kIdealNum2,
        kIdealDen,
        kShearNum1,
        kShearNum2,
        kShearDen,
        kFullNum1,
        kFullNum2,
        kFullDen,
        kUPerpNum,
        kUPerpDen,
        kTAvgNum,
        kTAvgDen,
        kRPiNum,
        kRPiDen,
        kRShearNum,
        kRShearDen,
        // epsilon_{p2} and epsilon_{p3} from the T^{0 mu} vector,
        // epsilon_{pn} = (int T^{0r} exp(i n phi_u)) / (int T^{0r}), for
        // ideal, ideal + shear, and full T^{mu nu}
        kEpNum1,
        kEpNum2 = kEpNum1 + 6,
        kEpDen = kEpNum2 + 6,
        kEccNum1 = kEpDen + 6,
        kEccNum2 = kEccNum1 + norder,
        kEccDen = kEccNum2 + norder,
        kMeanNum = kEccDen + norder,
        kMeanDen = kMeanNum + 4,
        kNSums = kMeanDen + 1
    };
    std::vector<int> rowEta;
    for (int ieta = 0; ieta < arena.nEta(); ieta++) {
        double eta = 0.0;
        if (!DATA.boost_invariant) {
            eta = DATA.etaAt(ieta);
        }
        if (eta < eta_max && eta > eta_min) rowEta.push_back(ieta);
    }
    const int ny = arena.nY();
    const int nRows = static_cast<int>(rowEta.size()) * ny;

    // the centre of mass
    std::vector<std::array<double, 3>> rowCm(nRows);
#pragma omp parallel for schedule(dynamic)
    for (int row = 0; row < nRows; row++) {
        const int ieta = rowEta[row / ny];
        const int iy = row % ny;
        std::array<double, 3> cm = {0.0, 0.0, 0.0};
        for (int ix = 0; ix < arena.nX(); ix++) {
            int fieldIdx = arena.getFieldIdx(ix, iy, ieta);
            double x_local = -DATA.x_size / 2. + ix * DATA.delta_x;
            double y_local = -DATA.y_size / 2. + iy * DATA.delta_y;
            double e_local = arena.e_[fieldIdx];  // 1/fm^4
            double gamma_perp = arena.u_[0][fieldIdx];
            cm[0] += x_local * e_local * gamma_perp;
            cm[1] += y_local * e_local * gamma_perp;
            cm[2] += e_local * gamma_perp;
        }
        rowCm[row] = cm;
    }
    double x_o = 0.0;
    double y_o = 0.0;
    double w_sum = 0.0;
    for (int row = 0; row < nRows; row++) {
        x_o += rowCm[row][0];
        y_o += rowCm[row][1];
        w_sum += rowCm[row][2];
    }
    x_o /= w_sum;
    y_o /= w_sum;

    std::vector<std::array<double, kNSums>> rowSums(nRows);
#pragma omp parallel for schedule(dynamic)
    for (int row = 0; row < nRows; row++) {
        const int ieta = rowEta[row / ny];
        const int iy = row % ny;
        std::array<double, kNSums> sums;
        sums.fill(0.0);
        double &ideal_num1 = sums[kIdealNum1];
        double &ideal_num2 = sums[kIdealNum2];
        double &ideal_den = sums[kIdealDen];
        double &shear_num1 = sums[kShearNum1];
        double &shear_num2 = sums[kShearNum2];
        double &shear_den = sums[kShearDen];
        double &full_num1 = sums[kFullNum1];
        double &full_num2 = sums[kFullNum2];
        double &full_den = sums[kFullDen];
        double &u_perp_num = sums[kUPerpNum];
        double &u_perp_den = sums[kUPerpDen];
        double &T_avg_num = sums[kTAvgNum];
        double &T_avg_den = sums[kTAvgDen];
        double &R_Pi_num = sums[kRPiNum];
        double &R_Pi_den = sums[kRPiDen];
        double &R_shearpi_num = sums[kRShearNum];
        double &R_shearpi_den = sums[kRShearDen];
        double *ep_num1 = &sums[kEpNum1];
        double *ep_num2 = &sums[kEpNum2];
        double *ep_den = &sums[kEpDen];
        double *eccn_num1 = &sums[kEccNum1];
        double *eccn_num2 = &sums[kEccNum2];
        double *eccn_den = &sums[kEccDen];
        double *meanpT_est_num = &sums[kMeanNum];
        double *meanpT_est_den = &sums[kMeanDen];
        ThermalVec thermalVec;
        for (int ix = 0; ix < arena.nX(); ix++) {
            int fieldIdx = arena.getFieldIdx(ix, iy, ieta);
            double x_local = (-DATA.x_size / 2. + ix * DATA.delta_x - x_o);
            double y_local = (-DATA.y_size / 2. + iy * DATA.delta_y - y_o);
            double r_local = sqrt(x_local * x_local + y_local * y_local);
            double phi_local = atan2(y_local, x_local);

            double e_local = arena.e_[fieldIdx];  // 1/fm^4
            eos.getThermalVariables(
                e_local, arena.rhob_[fieldIdx], arena.rhoq_[fieldIdx],
                arena.rhos_[fieldIdx], thermalVec);

            double P_local = thermalVec[2];
            double enthopy = e_local + P_local;
            double u0 = arena.u_[0][fieldIdx];
            double ux = arena.u_[1][fieldIdx];
            double uy = arena.u_[2][fieldIdx];
            double pi_0x = arena.Wmunu_[1][fieldIdx];
            double pi_0y = arena.Wmunu_[2][fieldIdx];
            double pi_xx = arena.Wmunu_[4][fieldIdx];
            double pi_xy = arena.Wmunu_[5][fieldIdx];
            double pi_yy = arena.Wmunu_[7][fieldIdx];
            double bulk_Pi = arena.piBulk_[fieldIdx];

            double T_00_ideal = enthopy * u0 * u0 - P_local;
            double T_0x_ideal = enthopy * u0 * ux;
            double T_0y_ideal = enthopy * u0 * uy;
            double T_0r_ideal =
                sqrt(T_0x_ideal * T_0x_ideal + T_0y_ideal * T_0y_ideal);
            double phi_u_ideal = atan2(T_0y_ideal, T_0x_ideal);
            double T_xx_ideal = enthopy * ux * ux + P_local;
            double T_xy_ideal = enthopy * ux * uy;
            double T_yy_ideal = enthopy * uy * uy + P_local;

            double T_0x_shear = T_0x_ideal + pi_0x;
            double T_0y_shear = T_0y_ideal + pi_0y;
            double T_0r_shear =
                sqrt(T_0x_shear * T_0x_shear + T_0y_shear * T_0y_shear);
            double phi_u_shear = atan2(T_0y_shear, T_0x_shear);
            double T_xx_shear = T_xx_ideal + pi_xx;
            double T_xy_shear = T_xy_ideal + pi_xy;
            double T_yy_shear = T_yy_ideal + pi_yy;

            double T_0x_full = T_0x_shear + bulk_Pi * u0 * ux;
            double T_0y_full = T_0y_shear + bulk_Pi * u0 * uy;
            double T_0r_full =
                sqrt(T_0x_full * T_0x_full + T_0y_full * T_0y_full);
            double phi_u_full = atan2(T_0y_full, T_0x_full);
            double T_xx_full = T_xx_shear - bulk_Pi * (-1 - ux * ux);
            double T_xy_full = T_xy_shear + bulk_Pi * ux * uy;
            double T_yy_full = T_yy_shear - bulk_Pi * (-1 - uy * uy);

            ideal_num1 += T_xx_ideal - T_yy_ideal;
            ideal_num2 += 2. * T_xy_ideal;
            ideal_den += T_xx_ideal + T_yy_ideal;
            shear_num1 += T_xx_shear - T_yy_shear;
            shear_num2 += 2. * T_xy_shear;
            shear_den += T_xx_shear + T_yy_shear;
            full_num1 += T_xx_full - T_yy_full;
            full_num2 += 2. * T_xy_full;
            full_den += T_xx_full + T_yy_full;

            double weight_local = e_local;
            u_perp_num += weight_local * u0;
            u_perp_den += weight_local;
            T_avg_num += weight_local * thermalVec[6];
            T_avg_den += weight_local;

            if (e_local > 1e-3) {
                double r_shearpi, r_bulkPi;
                calculate_inverse_Reynolds_numbers(
                    arena, fieldIdx, e_local + P_local, r_shearpi, r_bulkPi);
                R_shearpi_num += weight_local * r_shearpi;
                R_shearpi_den += weight_local;
                R_Pi_num += weight_local * r_bulkPi;
                R_Pi_den += weight_local;
            }

            for (int i = 0; i < 2; i++) {
                int idx = 3 * i;
                int iorder = 2 + i;
                ep_num1[idx] += T_0r_ideal * cos(iorder * phi_u_ideal);
                ep_num2[idx] += T_0r_ideal * sin(iorder * phi_u_ideal);
                ep_den[idx] += T_0r_ideal;
                ep_num1[idx + 1] += T_0r_shear * cos(iorder * phi_u_shear);
                ep_num2[idx + 1] += T_0r_shear * sin(iorder * phi_u_shear);
                ep_den[idx + 1] += T_0r_shear;
                ep_num1[idx + 2] += T_0r_full * cos(iorder * phi_u_full);
                ep_num2[idx + 2] += T_0r_full * sin(iorder * phi_u_full);
                ep_den[idx + 2] += T_0r_full;
            }
            for (int i = 1; i <= norder; i++) {
                if (i == 1) {
                    weight_local = u0 * e_local * pow(r_local, 3);
                } else {
                    weight_local = u0 * e_local * pow(r_local, i);
                }
                eccn_num1[i - 1] += weight_local * cos(i * phi_local);
                eccn_num2[i - 1] += weight_local * sin(i * phi_local);
                eccn_den[i - 1] += weight_local;
            }

            meanpT_est_num[0] += tau * thermalVec[12] * u0;         // dS/deta_s
            meanpT_est_num[1] += tau * T_00_ideal;                  // dE/deta_s
            meanpT_est_num[2] += thermalVec[12] * u0 * e_local;     // [s]
            meanpT_est_num[3] += r_local * r_local * u0 * e_local;  // [r^2]
            meanpT_est_den[0] += u0 * e_local;
        }
        rowSums[row] = sums;
    }
    std::array<double, kNSums> total;
    total.fill(0.0);
    for (int row = 0; row < nRows; row++) {
        for (int k = 0; k < kNSums; k++) total[k] += rowSums[row][k];
    }
    const double ideal_num1 = total[kIdealNum1];
    const double ideal_num2 = total[kIdealNum2];
    const double ideal_den = total[kIdealDen];
    const double shear_num1 = total[kShearNum1];
    const double shear_num2 = total[kShearNum2];
    const double shear_den = total[kShearDen];
    const double full_num1 = total[kFullNum1];
    const double full_num2 = total[kFullNum2];
    const double full_den = total[kFullDen];
    const double u_perp_num = total[kUPerpNum];
    const double u_perp_den = total[kUPerpDen];
    const double T_avg_num = total[kTAvgNum];
    const double T_avg_den = total[kTAvgDen];
    const double R_Pi_num = total[kRPiNum];
    const double R_Pi_den = total[kRPiDen];
    const double R_shearpi_num = total[kRShearNum];
    const double R_shearpi_den = total[kRShearDen];
    const double *ep_num1 = &total[kEpNum1];
    const double *ep_num2 = &total[kEpNum2];
    const double *ep_den = &total[kEpDen];
    const double *eccn_num1 = &total[kEccNum1];
    const double *eccn_num2 = &total[kEccNum2];
    const double *eccn_den = &total[kEccDen];
    const double *meanpT_est_num = &total[kMeanNum];
    const double *meanpT_est_den = &total[kMeanDen];
    double R_shearpi = R_shearpi_num / std::max(R_shearpi_den, small_eps);
    double R_Pi = R_Pi_num / std::max(R_Pi_den, small_eps);
    double u_avg = u_perp_num / std::max(u_perp_den, small_eps);
    double T_avg = T_avg_num / std::max(T_avg_den, small_eps) * hbarc;

    of << scientific << setw(18) << setprecision(8) << tau << "  "
       << ideal_num1 / std::max(ideal_den, small_eps) << "  "
       << ideal_num2 / std::max(ideal_den, small_eps) << "  "
       << shear_num1 / std::max(shear_den, small_eps) << "  "
       << shear_num2 / std::max(shear_den, small_eps) << "  "
       << full_num1 / std::max(full_den, small_eps) << "  "
       << full_num2 / std::max(full_den, small_eps) << "  ";
    for (int i = 0; i < 6; i++) {
        of << ep_num1[i] / std::max(ep_den[i], small_eps) << "  "
           << ep_num2[i] / std::max(ep_den[i], small_eps) << "  ";
    }
    of << endl;
    of.close();

    of1 << scientific << setw(18) << setprecision(8) << tau << "  ";
    for (int i = 0; i < norder; i++) {
        // the minus sign ensure the vector points to the short axis
        of1 << -eccn_num1[i] / std::max(eccn_den[i], small_eps) << "  "
            << -eccn_num2[i] / std::max(eccn_den[i], small_eps) << "  ";
    }
    of1 << x_o << "  " << y_o;
    of1 << endl;
    of1.close();

    of2 << scientific << setw(18) << setprecision(8) << tau << "  " << R_shearpi
        << "  " << R_Pi << "  " << u_avg << "  " << T_avg << endl;
    of2.close();

    of3 << scientific << setw(18) << setprecision(8) << tau << "  "
        << meanpT_est_num[0] * DATA.delta_x * DATA.delta_y << "  "
        << meanpT_est_num[1] * hbarc * DATA.delta_x * DATA.delta_y << "  "
        << meanpT_est_num[2] / std::max(meanpT_est_den[0], small_eps) << "  "
        << meanpT_est_num[3] / std::max(meanpT_est_den[0], small_eps) << endl;
    of3.close();
}

void Cell_info::get_LRF_shear_stress_tensor(
    const Fields &arena, const int Idx, const double eta_s,
    ShearVisVecLRF &res) {
    const double cosh_eta = cosh(eta_s);
    const double sinh_eta = sinh(eta_s);

    double u0 = arena.u_[0][Idx];
    double ux = arena.u_[1][Idx];
    double uy = arena.u_[2][Idx];
    double u3 = arena.u_[3][Idx];
    double ut = u0 * cosh_eta + u3 * sinh_eta;
    double uz = u3 * cosh_eta + u0 * sinh_eta;
    double LorentzBoost[4][4] = {
        {ut, -ux, -uy, -uz},
        {-ux, 1. + ux * ux / (ut + 1.), ux * uy / (ut + 1.),
         ux * uz / (ut + 1.)},
        {-uy, ux * uy / (ut + 1.), 1. + uy * uy / (ut + 1.),
         uy * uz / (ut + 1.)},
        {-uz, ux * uz / (ut + 1.), uy * uz / (ut + 1.),
         1. + uz * uz / (ut + 1.)}};

    ViscousVec ShearVisVec;
    for (unsigned int i = 0; i < ShearVisVec.size(); i++) {
        ShearVisVec[i] = arena.Wmunu_[i][Idx];
    }
    double pi_tz[4][4];
    pi_tz[0][0] =
        (ShearVisVec[0] * cosh_eta * cosh_eta
         + 2. * ShearVisVec[3] * cosh_eta * sinh_eta
         + ShearVisVec[9] * sinh_eta * sinh_eta);
    pi_tz[0][1] = ShearVisVec[1] * cosh_eta + ShearVisVec[6] * sinh_eta;
    pi_tz[0][2] = ShearVisVec[2] * cosh_eta + ShearVisVec[8] * sinh_eta;
    pi_tz[0][3] =
        (ShearVisVec[0] * cosh_eta * sinh_eta
         + ShearVisVec[3] * (cosh_eta * cosh_eta + sinh_eta * sinh_eta)
         + ShearVisVec[9] * sinh_eta * cosh_eta);
    pi_tz[1][0] = pi_tz[0][1];
    pi_tz[1][1] = ShearVisVec[4];
    pi_tz[1][2] = ShearVisVec[5];
    pi_tz[1][3] = ShearVisVec[1] * sinh_eta + ShearVisVec[6] * cosh_eta;
    pi_tz[2][0] = pi_tz[0][2];
    pi_tz[2][1] = pi_tz[1][2];
    pi_tz[2][2] = ShearVisVec[7];
    pi_tz[2][3] = ShearVisVec[2] * sinh_eta + ShearVisVec[8] * cosh_eta;
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
                    pi_LRF[i][j] +=
                        (LorentzBoost[i][a] * pi_tz[a][b] * LorentzBoost[b][j]);
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
    q_tz[0] = ShearVisVec[10] * cosh_eta + ShearVisVec[13] * sinh_eta;
    q_tz[1] = ShearVisVec[11];
    q_tz[2] = ShearVisVec[12];
    q_tz[3] = ShearVisVec[10] * sinh_eta + ShearVisVec[13] * cosh_eta;
    double q_LRF[4];
    for (int i = 1; i < 4; i++) {
        q_LRF[i] = 0.;
        for (int a = 0; a < 4; a++) q_LRF[i] += LorentzBoost[i][a] * q_tz[a];
    }
    res[5] = q_LRF[1];
    res[6] = q_LRF[2];
    res[7] = q_LRF[3];
}
