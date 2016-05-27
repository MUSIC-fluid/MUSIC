#include "util.h"
#include "grid_info.h"

using namespace std;

Grid_info::Grid_info(InitData* DATA_in) {
    DATA_ptr = DATA_in;

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

Grid_info::~Grid_info() {
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

// Output for MARTINI when the option Coordinates is set to "taueta".
// Because fluctuations in initial conditions are taken into account here,
// the range is now all x, y, and eta, not restricted to positive values.
// However, to save space, x, y, and eta are not written into the output, and
// the interpolation in MARTINI must match properly with the output generated
// here. -CFY 11/12/2010
void Grid_info::OutputEvolutionDataXYEta(Grid ***arena, InitData *DATA,
                                         EOS *eos, double tau) {
    const string out_name_xyeta = "evolution_xyeta.dat";
    const string out_name_W_xyeta = 
                        "evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
    const string out_name_q_xyeta = "evolution_qmu_xyeta.dat";
    string out_open_mode;
    FILE *out_file_xyeta;
    FILE *out_file_W_xyeta;
    FILE *out_file_q_xyeta;

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
    if (DATA->turn_on_diff == 1) {
        out_file_q_xyeta = fopen(out_name_q_xyeta.c_str(),
                                 out_open_mode.c_str());
    }
      
    int ix, iy, ieta;
    int n_skip_x = DATA->output_evolution_every_N_x;
    int n_skip_y = DATA->output_evolution_every_N_y;
    int n_skip_eta = DATA->output_evolution_every_N_eta;
    for (ieta = 0; ieta < DATA->neta; ieta+=n_skip_eta) {
	    double eta = ((double)ieta)*(DATA->delta_eta)-(DATA->eta_size)/2.0;
        double cosh_eta = cosh(eta);
        double sinh_eta = sinh(eta);
	    for (iy = 0; iy <= DATA->ny; iy += n_skip_y) {
	        for (ix = 0; ix <= DATA->nx; ix += n_skip_x) {
                double e_local = arena[ieta][ix][iy].epsilon;  // 1/fm^4
                double p_local = arena[ieta][ix][iy].p;        // 1/fm^4
                double rhob_local = arena[ieta][ix][iy].rhob;  // 1/fm^3
                double utau = arena[ieta][ix][iy].u[0][0];
                double ux = arena[ieta][ix][iy].u[0][1];
                double uy = arena[ieta][ix][iy].u[0][2];
                double ueta = arena[ieta][ix][iy].u[0][3];
                double ut = utau*cosh_eta + ueta*sinh_eta;  // gamma factor
                double vx = ux/ut;
                double vy = uy/ut;
                double uz = ueta*cosh_eta + utau*sinh_eta;
                double vz = uz/ut;
                
                double T_local = arena[ieta][ix][iy].T;  // 1/fm
                double muB_local = arena[ieta][ix][iy].mu;  // 1/fm
                double div_factor = e_local + p_local;  // 1/fm^4

                double Wtautau, Wtaux, Wtauy, Wtaueta;
                double Wxx, Wxy, Wxeta;
                double Wyy, Wyeta;
                double Wetaeta;
                if (DATA->turn_on_shear == 1) {
                    Wtautau = arena[ieta][ix][iy].Wmunu[0][0][0]/div_factor;
                    Wtaux = arena[ieta][ix][iy].Wmunu[0][0][1]/div_factor;
                    Wtauy = arena[ieta][ix][iy].Wmunu[0][0][2]/div_factor;
                    Wtaueta = arena[ieta][ix][iy].Wmunu[0][0][3]/div_factor;
                    Wxx = arena[ieta][ix][iy].Wmunu[0][1][1]/div_factor;
                    Wxy = arena[ieta][ix][iy].Wmunu[0][1][2]/div_factor;
                    Wxeta = arena[ieta][ix][iy].Wmunu[0][1][3]/div_factor;
                    Wyy = arena[ieta][ix][iy].Wmunu[0][2][2]/div_factor;
                    Wyeta = arena[ieta][ix][iy].Wmunu[0][2][3]/div_factor;
                    Wetaeta = arena[ieta][ix][iy].Wmunu[0][3][3]/div_factor;
                }
                
                // outputs for baryon diffusion part
                double common_term_q;
                double qtau, qx, qy, qeta;
                if (DATA->turn_on_diff == 1) {
                    common_term_q = rhob_local*T_local/div_factor;
                    double kappa_hat = get_deltaf_qmu_coeff(T_local,
                                                            muB_local);
                    qtau = arena[ieta][ix][iy].Wmunu[0][4][0]/kappa_hat;
                    qx = arena[ieta][ix][iy].Wmunu[0][4][1]/kappa_hat;
                    qy = arena[ieta][ix][iy].Wmunu[0][4][2]/kappa_hat;
                    qeta = arena[ieta][ix][iy].Wmunu[0][4][3]/kappa_hat;
                }

                // exclude the actual coordinates from the output to save space:
                if (0 == DATA->outputBinaryEvolution) {
                    fprintf(out_file_xyeta, "%e %e %e %e %e\n",
                            T_local*hbarc, muB_local*hbarc, vx, vy, vz);
                    if (1 == DATA->viscosity_flag) {
		                fprintf(out_file_W_xyeta,
                                "%e %e %e %e %e %e %e %e %e %e\n",
                                Wtautau, Wtaux, Wtauy, Wtaueta, Wxx, Wxy,
                                Wxeta, Wyy, Wyeta, Wetaeta); 
                    }
                } else {
		            double array[] = {T_local*hbarc, muB_local*hbarc,
                                      vx, vy, vz};
                    fwrite(array, sizeof(double), 5, out_file_xyeta);
		            if (1 == DATA->viscosity_flag) {
		  	            double array2[] = {
                            Wtautau, Wtaux, Wtauy, Wtaueta, Wxx, Wxy, Wxeta,
                            Wyy, Wyeta, Wetaeta};
		  	            fwrite(array2, sizeof(double), 10, out_file_W_xyeta);
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
    if (DATA->turn_on_diff == 1) {
        fclose(out_file_q_xyeta);
    }
}/* OutputEvolutionDataXYEta */

void Grid_info::check_conservation_law(Grid ***arena, InitData *DATA,
                                       double tau) {
    int nx = DATA->nx;
    int ny = DATA->ny;
    int neta = DATA->neta;
    double dx = DATA->delta_x;
    double dy = DATA->delta_y;
    double deta = DATA->delta_eta;

    double N_B = 0.0;
    double T_tau_t = 0.0;
    for (int ieta = 0; ieta < neta; ieta++) {
        double eta_s = deta*ieta - (DATA->eta_size)/2.0;
        for (int ix = 0; ix <= nx; ix++) {
            for (int iy = 0; iy <= ny; iy++) {
                double cosh_eta = cosh(eta_s);
                double sinh_eta = sinh(eta_s);
                N_B += (arena[ix][iy][ieta].rhob
                        *arena[ix][iy][ieta].u[0][0]
                        + (arena[ix][iy][ieta].prevWmunu[0][4][0]
                           + arena[ix][iy][ieta].prevWmunu[1][4][0])*0.5);
                double T_tau_tau = (
                    arena[ix][iy][ieta].TJb[0][0][0]
                    + (arena[ix][iy][ieta].prevWmunu[0][0][0]
                       + arena[ix][iy][ieta].prevWmunu[1][0][0])*0.5);
                double T_tau_eta = (
                    arena[ix][iy][ieta].TJb[0][0][3]
                    + (arena[ix][iy][ieta].prevWmunu[0][0][3]
                       + arena[ix][iy][ieta].prevWmunu[1][0][3])*0.5);
                T_tau_t += T_tau_tau*cosh_eta + T_tau_eta*sinh_eta;
            }
        }
    }
    N_B *= tau*dx*dy*deta;
    T_tau_t *= tau*dx*dy*deta*0.19733;
    cout << "check: net baryon number N_B = " << N_B << endl;
    cout << "check: total energy T^{taut} = " << T_tau_t << " GeV" << endl;
}

void Grid_info::check_velocity_shear_tensor(Grid ***arena, double tau) {
    ostringstream filename;
    filename << "Check_velocity_shear_tensor_tau_" << tau << ".dat";
    ofstream output_file(filename.str().c_str());

    double unit_convert = 0.19733;  // hbarC
    double dx = DATA_ptr->delta_x;
    double x_min = - DATA_ptr->nx/2.*dx;
    for (int ix = 0; ix <= DATA_ptr->nx; ix++) {
        double x_local = x_min + ix*dx;
        double e_local = arena[0][ix][ix].epsilon;
        output_file << scientific << setprecision(8) << setw(18)
                    << x_local << "  "
                    << e_local*unit_convert << "  "
                    << arena[0][ix][ix].u[0][0] << "  "
                    << arena[0][ix][ix].u[0][1] << "  "
                    << arena[0][ix][ix].u[0][2] << "  "
                    << arena[0][ix][ix].u[0][3] << "  "
                    << arena[0][ix][ix].sigma[0][0][0] << "  "
                    << arena[0][ix][ix].sigma[0][0][1] << "  "
                    << arena[0][ix][ix].sigma[0][0][2] << "  "
                    << arena[0][ix][ix].sigma[0][0][3] << "  "
                    << arena[0][ix][ix].sigma[0][1][1] << "  "
                    << arena[0][ix][ix].sigma[0][1][2] << "  "
                    << arena[0][ix][ix].sigma[0][1][3] << "  "
                    << arena[0][ix][ix].sigma[0][2][2] << "  "
                    << arena[0][ix][ix].sigma[0][2][3] << "  "
                    << arena[0][ix][ix].sigma[0][3][3] << "  "
                    << arena[0][ix][ix].prev_u[0][0] << "  "
                    << arena[0][ix][ix].prev_u[0][1] << "  "
                    << arena[0][ix][ix].prev_u[0][2] << "  "
                    << arena[0][ix][ix].prev_u[0][3]
                    << endl;
    }
    output_file.close();
}

void Grid_info::Gubser_flow_check_file(Grid ***arena, EOS *eos, double tau) {
    ostringstream filename;
    filename << "Gubser_flow_check_tau_" << tau << ".dat";
    ofstream output_file(filename.str().c_str());

    double unit_convert = 0.19733;  // hbarC
    double dx = 0.05;
    double x_min = -5.0;
    double dy = 0.05;
    double y_min = -5.0;
    for (int ix = 0; ix <= DATA_ptr->nx; ix++) {
        double x_local = x_min + ix*dx;
        for (int iy = 0; iy <= DATA_ptr->ny; iy++) {
            double y_local = y_min + iy*dy;
            double e_local = arena[0][ix][iy].epsilon;
            double T_local = eos->get_temperature(e_local, 0.0);
            output_file << scientific << setprecision(8) << setw(18)
                        << x_local << "  " << y_local << "  "
                        << e_local*unit_convert << "  "
                        << T_local*unit_convert << "  "
                        << arena[0][ix][iy].u[0][1] << "  "
                        << arena[0][ix][iy].u[0][2] << "  "
                        << arena[0][ix][iy].Wmunu[0][1][1]*unit_convert << "  "
                        << arena[0][ix][iy].Wmunu[0][2][2]*unit_convert << "  "
                        << arena[0][ix][iy].Wmunu[0][1][2]*unit_convert << "  "
                        << arena[0][ix][iy].Wmunu[0][3][3]*unit_convert << "  "
                        << endl;
        }
    }
    output_file.close();
}

void Grid_info::load_deltaf_qmu_coeff_table(string filename)
{
    ifstream table(filename.c_str());
    deltaf_qmu_coeff_table_length_T = 150;
    deltaf_qmu_coeff_table_length_mu = 100;
    delta_qmu_coeff_table_T0 = 0.05;
    delta_qmu_coeff_table_mu0 = 0.0;
    delta_qmu_coeff_table_dT = 0.001;
    delta_qmu_coeff_table_dmu = 0.007892;
    deltaf_qmu_coeff_tb = new double* [deltaf_qmu_coeff_table_length_T];
    for(int i = 0; i < deltaf_qmu_coeff_table_length_T; i++)
       deltaf_qmu_coeff_tb[i] = new double [deltaf_qmu_coeff_table_length_mu];

    double dummy;
    for(int j = 0; j < deltaf_qmu_coeff_table_length_mu; j++)
       for(int i = 0; i < deltaf_qmu_coeff_table_length_T; i++)
          table >> dummy >> dummy >> deltaf_qmu_coeff_tb[i][j];
    table.close();
}

void Grid_info::load_deltaf_qmu_coeff_table_14mom(string filename)
{
    ifstream table(filename.c_str());
    deltaf_coeff_table_14mom_length_T = 190;
    deltaf_coeff_table_14mom_length_mu = 160;
    delta_coeff_table_14mom_T0 = 0.01;
    delta_coeff_table_14mom_mu0 = 0.0;
    delta_coeff_table_14mom_dT = 0.001;
    delta_coeff_table_14mom_dmu = 0.005;

    deltaf_coeff_tb_14mom_DPi = new double* [deltaf_coeff_table_14mom_length_T];
    deltaf_coeff_tb_14mom_BPi = new double* [deltaf_coeff_table_14mom_length_T];
    deltaf_coeff_tb_14mom_BPitilde = new double* [deltaf_coeff_table_14mom_length_T];
    deltaf_coeff_tb_14mom_DV = new double* [deltaf_coeff_table_14mom_length_T];
    deltaf_coeff_tb_14mom_BV = new double* [deltaf_coeff_table_14mom_length_T];
    deltaf_coeff_tb_14mom_Bpi_shear = new double* [deltaf_coeff_table_14mom_length_T];
    for(int i = 0; i < deltaf_coeff_table_14mom_length_T; i++)
    {
       deltaf_coeff_tb_14mom_DPi[i] = new double [deltaf_coeff_table_14mom_length_mu];
       deltaf_coeff_tb_14mom_BPi[i] = new double [deltaf_coeff_table_14mom_length_mu];
       deltaf_coeff_tb_14mom_BPitilde[i] = new double [deltaf_coeff_table_14mom_length_mu];
       deltaf_coeff_tb_14mom_DV[i] = new double [deltaf_coeff_table_14mom_length_mu];
       deltaf_coeff_tb_14mom_BV[i] = new double [deltaf_coeff_table_14mom_length_mu];
       deltaf_coeff_tb_14mom_Bpi_shear[i] = new double [deltaf_coeff_table_14mom_length_mu];
    }

    double dummy;
    for(int i = 0; i < deltaf_coeff_table_14mom_length_T; i++)
       for(int j = 0; j < deltaf_coeff_table_14mom_length_mu; j++)
          table >> dummy >> dummy >> deltaf_coeff_tb_14mom_DPi[i][j]
                >> deltaf_coeff_tb_14mom_BPi[i][j] >> deltaf_coeff_tb_14mom_BPitilde[i][j]
                >> deltaf_coeff_tb_14mom_DV[i][j] >> deltaf_coeff_tb_14mom_BV[i][j]
                >> deltaf_coeff_tb_14mom_Bpi_shear[i][j];
    table.close();

    // convert units
    double hbarc3 = hbarc*hbarc*hbarc;
    double hbarc4 = hbarc3*hbarc;
    for(int i = 0; i < deltaf_coeff_table_14mom_length_T; i++)
    {
       for(int j = 0; j < deltaf_coeff_table_14mom_length_mu; j++)
       {
          deltaf_coeff_tb_14mom_DPi[i][j] = deltaf_coeff_tb_14mom_DPi[i][j]*hbarc4;   // fm^4/GeV
          deltaf_coeff_tb_14mom_BPi[i][j] = deltaf_coeff_tb_14mom_BPi[i][j]*hbarc4;   // fm^4/(GeV^2)
          deltaf_coeff_tb_14mom_BPitilde[i][j] = deltaf_coeff_tb_14mom_BPitilde[i][j]*hbarc4;   // fm^4/(GeV^2)
          deltaf_coeff_tb_14mom_DV[i][j] = deltaf_coeff_tb_14mom_DV[i][j]*hbarc3;   // fm^3/GeV
          deltaf_coeff_tb_14mom_BV[i][j] = deltaf_coeff_tb_14mom_BV[i][j]*hbarc3;   // fm^3/(GeV^2)
          deltaf_coeff_tb_14mom_Bpi_shear[i][j] = deltaf_coeff_tb_14mom_Bpi_shear[i][j]*hbarc4;   // fm^4/(GeV^2)
       }
    }
}

double Grid_info::get_deltaf_qmu_coeff(double T, double muB) {
    if (muB < 0) {
       muB = -muB;
    }
    int idx_T = (int)((T - delta_qmu_coeff_table_T0)/delta_qmu_coeff_table_dT);
    int idx_mu = (int)((muB - delta_qmu_coeff_table_mu0)/delta_qmu_coeff_table_dmu);

    if (idx_T > deltaf_qmu_coeff_table_length_T - 2 || idx_T < 0)
        return(1.0);
    if (idx_mu > deltaf_qmu_coeff_table_length_mu - 2)
        return(1.0);

    double x_fraction = (T - delta_qmu_coeff_table_T0)/delta_qmu_coeff_table_dT - idx_T;
    double y_fraction = (muB - delta_qmu_coeff_table_mu0)/delta_qmu_coeff_table_dmu - idx_mu;

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

double Grid_info::get_deltaf_coeff_14moments(double T, double muB,
                                             double type) {
    int idx_T = (int)((T - delta_coeff_table_14mom_T0)/delta_coeff_table_14mom_dT);
    int idx_mu = (int)((muB - delta_coeff_table_14mom_mu0)/delta_coeff_table_14mom_dmu);
    double x_fraction = (T - delta_coeff_table_14mom_T0)/delta_coeff_table_14mom_dT - idx_T;
    double y_fraction = (muB - delta_coeff_table_14mom_mu0)/delta_coeff_table_14mom_dmu - idx_mu;

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
       cout << "Grid_info::get_deltaf_coeff_14moments: unknown type: "
            << type << endl;
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
