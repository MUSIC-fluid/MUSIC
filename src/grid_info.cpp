// Copyright Chun Shen @ 2014-2016
#include <string>
#include <iomanip>

#include "./util.h"
#include "./grid_info.h"

using namespace std;

Cell_info::Cell_info(InitData* DATA_in, EOS *eos_ptr_in) {
    DATA_ptr = DATA_in;
    eos_ptr = eos_ptr_in;

    // read in tables for delta f coefficients
    if (DATA_ptr->turn_on_diff == 1) {
        if (DATA_ptr->deltaf_14moments == 1) {
            load_deltaf_qmu_coeff_table_14mom(
                    "tables/deltaf_coefficients_14moments.dat");
        } else {
            if (DATA_ptr->include_deltaf_qmu == 1)
                load_deltaf_qmu_coeff_table(
                        "tables/Coefficients_RTA_diffusion.dat");
        }
    }
}

Cell_info::~Cell_info() {
    if (DATA_ptr->turn_on_diff == 1) {
        if (DATA_ptr->deltaf_14moments == 1) {
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
            if (DATA_ptr->include_deltaf_qmu == 1) {
                for (int i = 0; i < deltaf_qmu_coeff_table_length_T; i++)
                    delete [] deltaf_qmu_coeff_tb[i];
                delete [] deltaf_qmu_coeff_tb;
            }
        }
    }
}


//! This function outputs a header files for JF and Gojko's EM program
void Cell_info::Output_hydro_information_header(InitData *DATA) {
    string fname = "hydro_info_header_h";

    // Open output file
    ofstream outfile;
    outfile.open(fname.c_str());

    int grid_nx = ceil(
        (static_cast<double>(DATA->nx + 1))/DATA->output_evolution_every_N_x);
    int grid_ny = ceil(
        (static_cast<double>(DATA->ny + 1))/DATA->output_evolution_every_N_y);
    int grid_neta = ceil(
        (static_cast<double>(DATA->neta))/DATA->output_evolution_every_N_eta);

    outfile << "const int MUSIC_real_nx = " << grid_nx << ";" << endl;
    outfile << "const int MUSIC_real_ny = " << grid_ny << ";" << endl;
    outfile << "const int MUSIC_real_neta = " << grid_neta << ";" << endl;

    outfile << "const double MUSIC_tau0 = " << DATA->tau0 << ";" << endl;

    outfile << "const double MUSIC_dx = "
            << DATA->delta_x*DATA->output_evolution_every_N_x << ";" << endl;
    outfile << "const double MUSIC_dy = "
            << DATA->delta_y*DATA->output_evolution_every_N_y << ";" << endl;
    outfile << "const double MUSIC_deta = "
            << DATA->delta_eta*DATA->output_evolution_every_N_eta << ";"
            << endl;
    outfile << "const double MUSIC_dtau = "
            << DATA->output_evolution_every_N_timesteps*DATA->delta_tau << ";"
            << endl;

    outfile << "const bool MUSIC_with_shear_viscosity = "
            << ((DATA->viscosity_flag) && (DATA->turn_on_shear)) << ";\n";
    outfile << "const bool MUSIC_with_bulk_viscosity = "
            << ((DATA->viscosity_flag) && (DATA->turn_on_bulk)) << ";\n";
    outfile << "const bool MUSIC_with_rhob = " << (DATA->turn_on_rhob) << ";\n";
    outfile << "const bool MUSIC_with_diffusion = "
            << ((DATA->viscosity_flag) && (DATA->turn_on_diff)) << ";\n";

    outfile << "const bool MUSIC_outputBinaryEvolution="
            << DATA->outputBinaryEvolution << ";" << endl;

    outfile.close();
}


//! This function outputs hydro evolution file in binary format
void Cell_info::OutputEvolutionDataXYEta(Grid &arena, InitData *DATA,
                                         double tau) {
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
    if (tau == DATA->tau0) {
        out_open_mode = "w";
    } else {
        out_open_mode = "a";
    }
    // If we output in binary, set the mode accordingly
    if (0 == DATA->outputBinaryEvolution) {
        out_open_mode += "b";
    }
    out_file_xyeta = fopen(out_name_xyeta.c_str(), out_open_mode.c_str());
    if (DATA->turn_on_shear == 1) {
        out_file_W_xyeta = fopen(out_name_W_xyeta.c_str(),
                                 out_open_mode.c_str());
    }
    if (DATA->turn_on_bulk == 1) {
        out_file_bulkpi_xyeta = fopen(out_name_bulkpi_xyeta.c_str(),
                                      out_open_mode.c_str());
    }
    if (DATA->turn_on_diff == 1) {
        out_file_q_xyeta = fopen(out_name_q_xyeta.c_str(),
                                 out_open_mode.c_str());
    }
    int ix, iy, ieta;
    int n_skip_x = DATA->output_evolution_every_N_x;
    int n_skip_y = DATA->output_evolution_every_N_y;
    int n_skip_eta = DATA->output_evolution_every_N_eta;
    for (ieta = 0; ieta < DATA->neta; ieta += n_skip_eta) {
        double eta;
        if (DATA->boost_invariant == 1) {
            eta = 0.0;
        } else {
            eta = ((static_cast<double>(ieta))*(DATA->delta_eta)
                    - (DATA->eta_size)/2.0);
        }
        double cosh_eta = cosh(eta);
        double sinh_eta = sinh(eta);
        for (iy = 0; iy <= DATA->ny; iy += n_skip_y) {
            for (ix = 0; ix <= DATA->nx; ix += n_skip_x) {
                double e_local = arena(ix,iy,ieta).epsilon;  // 1/fm^4
                double rhob_local = arena(ix,iy,ieta).rhob;  // 1/fm^3
                double p_local = eos_ptr->get_pressure(e_local, rhob_local);
                double utau = arena(ix,iy,ieta).u[0][0];
                double ux = arena(ix,iy,ieta).u[0][1];
                double uy = arena(ix,iy,ieta).u[0][2];
                double ueta = arena(ix,iy,ieta).u[0][3];
                double ut = utau*cosh_eta + ueta*sinh_eta;  // gamma factor
                double vx = ux/ut;
                double vy = uy/ut;
                double uz = ueta*cosh_eta + utau*sinh_eta;
                double vz = uz/ut;

                double T_local = eos_ptr->get_temperature(e_local, rhob_local);
                double cs2_local = eos_ptr->get_cs2(e_local, rhob_local);
                double muB_local = eos_ptr->get_mu(e_local, rhob_local);
                double entropy = e_local + p_local;  // [1/fm^4]

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
                if (DATA->turn_on_shear == 1) {
                    Wtautau = arena(ix,iy,ieta).Wmunu[0][0]/entropy;
                    Wtaux   = arena(ix,iy,ieta).Wmunu[0][1]/entropy;
                    Wtauy   = arena(ix,iy,ieta).Wmunu[0][2]/entropy;
                    Wtaueta = arena(ix,iy,ieta).Wmunu[0][3]/entropy;
                    Wxx     = arena(ix,iy,ieta).Wmunu[0][4]/entropy;
                    Wxy     = arena(ix,iy,ieta).Wmunu[0][5]/entropy;
                    Wxeta   = arena(ix,iy,ieta).Wmunu[0][6]/entropy;
                    Wyy     = arena(ix,iy,ieta).Wmunu[0][7]/entropy;
                    Wyeta   = arena(ix,iy,ieta).Wmunu[0][8]/entropy;
                    Wetaeta = arena(ix,iy,ieta).Wmunu[0][9]/entropy;
                }

                double bulk_Pi = 0.0;
                if (DATA->turn_on_bulk == 1) {
                    bulk_Pi = arena(ix,iy,ieta).pi_b[0];  // [1/fm^4]
                }

                // outputs for baryon diffusion part
                double common_term_q = 0.0;
                double qtau          = 0.0;
                double qx            = 0.0;
                double qy            = 0.0;
                double qeta          = 0.0;
                if (DATA->turn_on_diff == 1) {
                    common_term_q = rhob_local*T_local/entropy;
                    double kappa_hat = get_deltaf_qmu_coeff(T_local,
                                                            muB_local);
                    qtau = arena(ix,iy,ieta).Wmunu[0][10]/kappa_hat;
                    qx   = arena(ix,iy,ieta).Wmunu[0][11]/kappa_hat;
                    qy   = arena(ix,iy,ieta).Wmunu[0][12]/kappa_hat;
                    qeta = arena(ix,iy,ieta).Wmunu[0][13]/kappa_hat;
                }

                // exclude the actual coordinates from the output to save space:
                if (DATA->outputBinaryEvolution == 0) {
                    fprintf(out_file_xyeta, "%e %e %e %e %e\n",
                            T_local*hbarc, muB_local*hbarc, vx, vy, vz);
                    if (DATA->viscosity_flag == 1) {
                        fprintf(out_file_W_xyeta,
                                "%e %e %e %e %e %e %e %e %e %e\n",
                                Wtautau, Wtaux, Wtauy, Wtaueta, Wxx, Wxy,
                                Wxeta, Wyy, Wyeta, Wetaeta);
                    }
                } else {
                    double array[] = {T_local*hbarc, muB_local*hbarc,
                                      vx, vy, vz};
                    fwrite(array, sizeof(double), 5, out_file_xyeta);
                    if (DATA->viscosity_flag == 1) {
                        if (DATA->turn_on_shear == 1) {
                            double array2[] = {Wtautau, Wtaux, Wtauy, Wtaueta,
                                               Wxx, Wxy, Wxeta, Wyy, Wyeta,
                                               Wetaeta};
                            fwrite(array2, sizeof(double), 10,
                                   out_file_W_xyeta);
                        }
                        if (DATA->turn_on_bulk == 1) {
                            double array1[] = {bulk_Pi, entropy, cs2_local};
                            fwrite(array1, sizeof(double), 3,
                                   out_file_bulkpi_xyeta);
                        }
                        if (DATA->turn_on_diff == 1) {
                            double array3[] = {common_term_q,
                                               qtau, qx, qy, qeta};
                            fwrite(array3, sizeof(double), 5,
                                   out_file_q_xyeta);
                        }
                    }
                }
            }/* ix */
        }/* iy */
    }/* ieta */
    fclose(out_file_xyeta);
    if (DATA->turn_on_shear == 1) {
        fclose(out_file_W_xyeta);
    }
    if (DATA->turn_on_bulk == 1) {
        fclose(out_file_bulkpi_xyeta);
    }
    if (DATA->turn_on_diff == 1) {
        fclose(out_file_q_xyeta);
    }
}/* OutputEvolutionDataXYEta */

    
//! This function outputs hydro evolution file in binary format
void Cell_info::OutputEvolutionDataXYEta_chun(Grid &arena, InitData *DATA,
                                              double tau) {
    // the format of the file is as follows,
    //    itau ix iy ieta T ux uy ueta
    // if turn_on_shear == 1:
    //    itau ix iy ieta T ux uy ueta Wxx Wxy Wxeta Wyy Wyeta
    // if turn_on_shear == 1 and turn_on_bulk == 1:
    //    itau ix iy ieta T ux uy ueta Wxx Wxy Wxeta Wyy Wyeta pi_b
    // if turn_on_rhob == 1:
    //    itau ix iy ieta T ux uy ueta mu_B
    // if turn_on_rhob == 1 and turn_on_shear == 1:
    //    itau ix iy ieta T ux uy ueta mu_B Wxx Wxy Wxeta Wyy Wyeta
    // if turn_on_rhob == 1 and turn_on_shear == 1 and turn_on_diff == 1:
    //    itau ix iy ieta T ux uy ueta mu_B Wxx Wxy Wxeta Wyy Wyeta qx qy qeta
    // Here ueta = tau*ueta, Wieta = tau*Wieta, qeta = tau*qeta
    // Here Wij is reduced variables Wij/(e+P) used in delta f
    // and qi is reduced variables qi/kappa_hat
    const string out_name_xyeta = "evolution_all_xyeta.dat";
    string out_open_mode;
    FILE *out_file_xyeta;
    // If it's the first timestep, overwrite the previous file
    if (tau == DATA->tau0) {
        out_open_mode = "wb";
    } else {
        out_open_mode = "ab";
    }
    out_file_xyeta = fopen(out_name_xyeta.c_str(), out_open_mode.c_str());

    int itau = static_cast<int>((tau - DATA->tau0)/DATA->delta_tau);

    int n_skip_x   = DATA->output_evolution_every_N_x;
    int n_skip_y   = DATA->output_evolution_every_N_y;
    int n_skip_eta = DATA->output_evolution_every_N_eta;
    for (int ieta = 0; ieta < DATA->neta; ieta += n_skip_eta) {
        for (int iy = 0; iy <= DATA->ny; iy += n_skip_y) {
            for (int ix = 0; ix <= DATA->nx; ix += n_skip_x) {
                double e_local = arena(ix,iy,ieta).epsilon;  // 1/fm^4
                double rhob_local = arena(ix,iy,ieta).rhob;  // 1/fm^3
                double p_local = eos_ptr->get_pressure(e_local, rhob_local);

                double ux   = arena(ix,iy,ieta).u[0][1];
                double uy   = arena(ix,iy,ieta).u[0][2];
                double ueta = arena(ix,iy,ieta).u[0][3];

                // T_local is in 1/fm
                double T_local = eos_ptr->get_temperature(e_local, rhob_local);


                if (T_local*hbarc < DATA->output_evolution_T_cut) continue;
                // only ouput fluid cells that are above cut-off temperature

                double muB_local = 0.0;
                if (DATA->turn_on_rhob == 1)
                    muB_local = eos_ptr->get_mu(e_local, rhob_local);

                double div_factor = e_local + p_local;  // 1/fm^4
                double Wxx = 0.0;
                double Wxy = 0.0;
                double Wxeta = 0.0;
                double Wyy = 0.0;
                double Wyeta = 0.0;
                if (DATA->turn_on_shear == 1) {
                    Wxx   = arena(ix,iy,ieta).Wmunu[0][4]/div_factor;
                    Wxy   = arena(ix,iy,ieta).Wmunu[0][5]/div_factor;
                    Wxeta = arena(ix,iy,ieta).Wmunu[0][6]/div_factor;
                    Wyy   = arena(ix,iy,ieta).Wmunu[0][7]/div_factor;
                    Wyeta = arena(ix,iy,ieta).Wmunu[0][8]/div_factor;
                }

                double pi_b = 0.0;
                if (DATA->turn_on_bulk == 1) {
                    pi_b = arena(ix,iy,ieta).pi_b[0];   // 1/fm^4
                }

                // outputs for baryon diffusion part
                //double common_term_q = 0.0;
                double qx = 0.0;
                double qy = 0.0;
                double qeta = 0.0;
                if (DATA->turn_on_diff == 1) {
                    //common_term_q = rhob_local*T_local/div_factor;
                    double kappa_hat = get_deltaf_qmu_coeff(T_local,
                                                            muB_local);
                    qx = arena(ix,iy,ieta).Wmunu[0][11]/kappa_hat;
                    qy = arena(ix,iy,ieta).Wmunu[0][12]/kappa_hat;
                    qeta = arena(ix,iy,ieta).Wmunu[0][13]/kappa_hat;
                }

                int pos[] = {itau, ix, iy, ieta};
                double ideal[] = {T_local*hbarc, ux, uy, ueta};

                fwrite(pos, sizeof(int), 4, out_file_xyeta);
                fwrite(ideal, sizeof(double), 4, out_file_xyeta);

                if (DATA->turn_on_rhob == 1) {
                    double mu[] = {muB_local*hbarc};
                    fwrite(mu, sizeof(double), 1, out_file_xyeta);
                }

                if (DATA->turn_on_shear == 1) {
                    double shear_pi[] = {Wxx, Wxy, Wxeta, Wyy, Wyeta};
                    fwrite(shear_pi, sizeof(double), 5, out_file_xyeta);
                }

                if (DATA->turn_on_bulk == 1) {
                    double bulk_pi[] = {pi_b};
                    fwrite(bulk_pi, sizeof(double), 1, out_file_xyeta);
                }

                if (DATA->turn_on_diff == 1) {
                    double diffusion[] = {qx, qy, qeta};
                    fwrite(diffusion, sizeof(double), 3, out_file_xyeta);
                }
            }  /* ix */
        }  /* iy */
    }  /* ieta */
    fclose(out_file_xyeta);
}/* OutputEvolutionDataXYEta */


//! This function prints to the screen the maximum local energy density,
//! the maximum temperature in the current grid
void Cell_info::get_maximum_energy_density(Grid &arena) {
    double eps_max  = 0.0;
    double rhob_max = 0.0;
    double T_max    = 0.0;

    // get the grid information
    int neta = DATA_ptr->neta;
    int nx   = DATA_ptr->nx;
    int ny   = DATA_ptr->ny;

    #pragma omp parallel for collapse(3) reduction(max:eps_max, rhob_max, T_max)
    for (int ieta = 0; ieta < neta; ieta++)
    for (int ix = 0; ix <= nx; ix++) 
    for (int iy = 0; iy <= ny; iy++) {
        const auto eps_local  = arena(ix,iy,ieta).epsilon;
        const auto rhob_local = arena(ix,iy,ieta).rhob;
        eps_max  = std::max(eps_max,  eps_local );
        rhob_max = std::max(rhob_max, rhob_local);
        T_max    = std::max(T_max,    eos_ptr->get_temperature(eps_local, rhob_local) );
    }
    eps_max *= 0.19733;   // GeV/fm^3
    T_max *= 0.19733;     // GeV
    music_message << "eps_max = " << eps_max << " GeV/fm^3, "
                  << "rhob_max = " << rhob_max << " 1/fm^3, "
                  << "T_max = " << T_max << " GeV.";
    music_message.flush("info");
}


//! This function checks the total energy and total net baryon number
//! at a give proper time
void Cell_info::check_conservation_law(Grid &arena, InitData *DATA,
                                      double tau) {
    double N_B     = 0.0;
    double T_tau_t = 0.0;
    int neta       = DATA->neta;
    double deta    = DATA->delta_eta;
    int    nx      = DATA->nx;
    double dx      = DATA->delta_x;
    int    ny      = DATA->ny;
    double dy      = DATA->delta_y;

    #pragma omp parallel for collapse(3) reduction(+:N_B, T_tau_t)
    for (int ieta = 0; ieta < neta; ieta++)     
    for (int ix = 0; ix <= nx; ix++) 
    for (int iy = 0; iy <= ny; iy++) {
        const auto& c = arena(ix,iy,ieta);

        const double eta_s = deta*ieta - (DATA->eta_size)/2.0;
        const double cosh_eta = cosh(eta_s);
        const double sinh_eta = sinh(eta_s);
        N_B += (c.rhob
                *c.u[0][0]
                + c.prevWmunu[0][10]);
        const double Pi00_rk_0 = (
            c.prev_pi_b[0]
            *(-1.0 + c.prev_u[0][0]
                     *(c.prev_u[0][0])));
        const double e_local   = c.epsilon;
        const double rhob      = c.rhob;
        const double pressure  = eos_ptr->get_pressure(e_local, rhob);
        const double u0        = c.u[0][0];
        const double u3        = c.u[0][3];
        const double T00_local = (e_local + pressure)*u0*u0 - pressure;
        const double T03_local = (e_local + pressure)*u0*u3;
        const double T_tau_tau = (T00_local
                + c.prevWmunu[0][0] + Pi00_rk_0);

        const double Pi03_rk_0 = (
            c.prev_pi_b[0]
            *(c.prev_u[0][0]
              *(c.prev_u[0][3])));
        const double T_tau_eta = (T03_local
                + c.prevWmunu[0][3] + Pi03_rk_0);
        T_tau_t += T_tau_tau*cosh_eta + T_tau_eta*sinh_eta;
    }
    double factor = tau*dx*dy*deta;
    N_B *= factor;
    T_tau_t *= factor*0.19733;  // GeV
    music_message << "net baryon number N_B = " << N_B;
    music_message.flush("info");
    music_message << "total energy T^{taut} = " << T_tau_t << " GeV";
    music_message.flush("info");
}



//! This function putputs files to check with Gubser flow solution
void Cell_info::Gubser_flow_check_file(Grid &arena, double tau) {
    double unit_convert = 0.19733;  // hbarC
    if (tau > 1.) {
        ostringstream filename_analytic;
        filename_analytic << "tests/Gubser_flow/y=x_tau="
                          << tau << "_SemiAnalytic.dat";

        double T_analytic[201], ux_analytic[201], uy_analytic[201];
        double pixx_analytic[201], pixy_analytic[201];
        double piyy_analytic[201], pizz_analytic[201];
        double dummy;
        ifstream input_file(filename_analytic.str().c_str());
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
        for (int i = 0; i < 201; i++) {
            double e_local = arena(i,i,0).epsilon;
            double T_local = (
                    eos_ptr->get_temperature(e_local, 0.0)*unit_convert);
            T_diff += fabs(T_analytic[i] - T_local);
            T_sum += fabs(T_analytic[i]);
            ux_diff += fabs(ux_analytic[i] - arena(i,i,0).u[0][1]);
            ux_sum += fabs(ux_analytic[i]);
            uy_diff += fabs(uy_analytic[i] - arena(i,i,0).u[0][2]);
            uy_sum += fabs(uy_analytic[i]);
            pixx_diff += (fabs(pixx_analytic[i]
                               - arena(i,i,0).Wmunu[0][4]*unit_convert));
            pixx_sum += fabs(pixx_analytic[i]);
            pixy_diff += (fabs(pixx_analytic[i]
                               - arena(i,i,0).Wmunu[0][5]*unit_convert));
            pixy_sum += fabs(pixx_analytic[i]);
            piyy_diff += (fabs(piyy_analytic[i]
                               - arena(i,i,0).Wmunu[0][7]*unit_convert));
            piyy_sum += fabs(piyy_analytic[i]);
            pizz_diff += (fabs(pizz_analytic[i]
                               - arena(i,i,0).Wmunu[0][9]*unit_convert));
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

    double dx = DATA_ptr->delta_x;
    double x_min = -DATA_ptr->x_size/2.;
    double dy = DATA_ptr->delta_y;
    double y_min = -DATA_ptr->y_size/2.;
    for (int ix = 0; ix <= DATA_ptr->nx; ix++) {
        double x_local = x_min + ix*dx;
        for (int iy = 0; iy <= DATA_ptr->ny; iy++) {
            double y_local = y_min + iy*dy;
            double e_local = arena(ix,iy,0).epsilon;
            double rhob_local = arena(ix,iy,0).rhob;
            double T_local = eos_ptr->get_temperature(e_local, 0.0);
            output_file << scientific << setprecision(8) << setw(18)
                        << x_local << "  " << y_local << "  "
                        << e_local*unit_convert << "  " << rhob_local << "  "
                        << T_local*unit_convert << "  "
                        << arena(ix,iy,0).u[0][1] << "  "
                        << arena(ix,iy,0).u[0][2] << "  "
                        << arena(ix,iy,0).Wmunu[0][4]*unit_convert << "  "
                        << arena(ix,iy,0).Wmunu[0][7]*unit_convert << "  "
                        << arena(ix,iy,0).Wmunu[0][5]*unit_convert << "  "
                        << arena(ix,iy,0).Wmunu[0][9]*unit_convert << "  "
                        << endl;
        }
    }
    output_file.close();
}


//! This function outputs files to cross check with 1+1D simulation
void Cell_info::output_1p1D_check_file(Grid &arena, double tau) {
    ostringstream filename;
    filename << "1+1D_check_tau_" << tau << ".dat";
    ofstream output_file(filename.str().c_str());

    double unit_convert = 0.19733;  // hbarC
    double deta = DATA_ptr->delta_eta;
    double eta_min = -6.94;
    for (int ieta = 0; ieta < DATA_ptr->neta; ieta++) {
        double eta_local = eta_min + ieta*deta;
        double e_local = arena(1,1,ieta).epsilon;
        double rhob_local = arena(1,1,ieta).rhob;
        output_file << scientific << setprecision(8) << setw(18)
                    << eta_local << "  "
                    << e_local*unit_convert << "  " << rhob_local
                    << endl;
    }
    output_file.close();
}


//! This function outputs energy density and n_b for making movies
void Cell_info::output_evolution_for_movie(Grid &arena, double tau) {
    ostringstream filename;
    filename << "movie_tau_"
             << setprecision(2) << fixed << tau << ".dat";
    output_energy_density_and_rhob_disitrubtion(arena, filename.str());
}


//! This function dumps the energy density and net baryon density
void Cell_info::output_energy_density_and_rhob_disitrubtion(Grid &arena,
                                                            string filename) {
    ofstream output_file(filename.c_str());
    double unit_convert = 0.19733;  // hbarC [GeV*fm]
    int n_skip_x = DATA_ptr->output_evolution_every_N_x;
    int n_skip_y = DATA_ptr->output_evolution_every_N_y;
    int n_skip_eta = DATA_ptr->output_evolution_every_N_eta;
    for (int ieta = 0; ieta < DATA_ptr->neta; ieta += n_skip_eta) {
        for (int ix = 0; ix <= DATA_ptr->nx; ix += n_skip_x) {
            for (int iy = 0; iy <= DATA_ptr->ny; iy += n_skip_y) {
                double e_local = arena(ix,iy,ieta).epsilon*unit_convert;
                double rhob_local = arena(ix,iy,ieta).rhob;
                output_file << scientific << setprecision(5) << setw(18)
                            << e_local << "  " << rhob_local << endl;
            }
        }
    }
    output_file.close();
}


//! This function outputs the evolution of hydrodynamic variables at a
//! give fluid cell
void Cell_info::monitor_fluid_cell(Grid &arena, int ix, int iy, int ieta,
                                   double tau) {
    ostringstream filename;
    filename << "monitor_fluid_cell_ix_" << ix << "_iy_" << iy
             << "_ieta_" << ieta << ".dat";
    ofstream output_file(filename.str().c_str(),
                         std::ofstream::out | std::ofstream::app);
    output_file << scientific << setprecision(8)
                << tau << "  " << arena(ix,iy,ieta).epsilon << "  "
                << arena(ix,iy,ieta).rhob << "  ";
    for (int i = 0; i < 4; i++) {
        output_file << scientific << setprecision(8)
                    << arena(ix,iy,ieta).u[0][i] << "  ";
    }
    for (int i = 0; i < 10; i++) {
        output_file << scientific << setprecision(8)
                    << arena(ix,iy,ieta).Wmunu[0][i] << "  ";
    }
    output_file << scientific << setprecision(8)
                << arena(ix,iy,ieta).pi_b[0] << "  ";
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            output_file << scientific << setprecision(8)
                        << arena(ix,iy,ieta).dUsup[i][j] << "  ";
        }
    }
    output_file << endl;
    output_file.close();
}

void Cell_info::load_deltaf_qmu_coeff_table(string filename) {
    ifstream table(filename.c_str());
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
    ifstream table(filename.c_str());
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
                double tau, double eta_min, double eta_max, Grid &arena) {
    ostringstream filename;
    filename << "averaged_phase_diagram_trajectory_eta_" << eta_min
             << "_" << eta_max << ".dat";
    fstream of(filename.str().c_str(), std::fstream::app | std::fstream::out);
    if (fabs(tau - DATA_ptr->tau0) < 1e-10) {
        of << "# tau(fm)  <T>(GeV)  std(T)(GeV)  <mu_B>(GeV)  std(mu_B)(GeV)"
           << endl;
    }
    double avg_T = 0.0;
    double avg_mu = 0.0;
    double std_T = 0.0;
    double std_mu = 0.0;
    double weight = 0.0;
    for (int ieta = 0; ieta < DATA_ptr->neta; ieta++) {
        double eta;
        if (DATA_ptr->boost_invariant == 1) {
            eta = 0.0;
        } else {
            eta = ((static_cast<double>(ieta))*(DATA_ptr->delta_eta)
                    - (DATA_ptr->eta_size)/2.0);
        }
        if (eta < eta_max && eta > eta_min) {
            double cosh_eta = cosh(eta);
            double sinh_eta = sinh(eta);
            for (int iy = 0; iy <= DATA_ptr->ny; iy++) {
                for (int ix = 0; ix <= DATA_ptr->nx; ix++) {
                    double e_local    = arena(ix,iy,ieta).epsilon;  // 1/fm^4
                    double rhob_local = arena(ix,iy,ieta).rhob;  // 1/fm^3
                    double utau       = arena(ix,iy,ieta).u[0][0];
                    double ueta       = arena(ix,iy,ieta).u[0][3];
                    double ut         = utau*cosh_eta + ueta*sinh_eta;  // gamma factor
                    double T_local    = eos_ptr->get_temperature(e_local,
                                                              rhob_local);
                    double muB_local = eos_ptr->get_mu(e_local, rhob_local);
                    double weight_local = e_local*ut;
                    avg_T += T_local*weight_local;
                    avg_mu += muB_local*weight_local;
                    std_T += T_local*T_local*weight_local;
                    std_mu += muB_local*muB_local*weight_local;
                    weight += weight_local;
                }
            }
        }
    }
    avg_T = avg_T/(weight + 1e-15)*hbarc;
    avg_mu = avg_mu/(weight + 1e-15)*hbarc;
    std_T = sqrt(std_T/(weight + 1e-15)*hbarc*hbarc - avg_T*avg_T);
    std_mu = sqrt(std_mu/(weight + 1e-15)*hbarc*hbarc - avg_mu*avg_mu);
    of << scientific << setw(18) << setprecision(8)
       << tau << "  " << avg_T << "  " << std_T << "  "
       << avg_mu << "  " << std_mu << endl;
    of.close();
}

