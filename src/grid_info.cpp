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
    int nx = DATA->nx;
    int ny = DATA->ny;
    int neta = DATA->neta;
    
    const string out_name_xyeta = "evolution_xyeta.dat";
    const string out_name_W_xyeta = 
                        "evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
    const string out_name_q_xyeta = "evolution_qmu_xyeta.dat";
    string out_open_mode;
    FILE *out_file_xyeta;
    FILE *out_file_W_xyeta;

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
                double muB = arena[ieta][ix][iy].mu;  // 1/fm
                double QGPfrac1 = 0.0;
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
                                Wtautau, Wtaux, Wtauy, Wtaueta, Wxx, Wxy, Wxz,
                                Wyy, Wyeta, Wetaeta); 
                    }
                } else {
		            double array[] = {T_local*hbarc, muB_local*hbarc,
                                      vx, vy, vz};
                    fwrite(array, sizeof(double), 5, out_file_xyeta);
		            if (1 == DATA->viscosity_flag) {
		  	            double array2[] = {
                            Wtautau, Wtaux, Wtauy, Wtaueta, Wxx, Wxy, Wxz,
                            Wyy, Wyeta, Wetaeta};
		  	            fwrite(array2, sizeof(double), 10, out_file_W_xyeta);
                        if (DATA->turn_on_diff == 1) {
                            float array3[] = {common_term_q,
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
      fclose(out_file_W_xyeta);
}/* OutputEvolutionDataXYEta */

void Grid_info::check_conservation_law(Grid ***arena, InitData *DATA,
                                       double tau, int size, int rank) {
    int position;
    int nx = DATA->nx;
    int ny = DATA->ny;
    int neta = DATA->neta;
    double dx = DATA->delta_x;
    double dy = DATA->delta_y;
    double deta = DATA->delta_eta;

    if (rank > 0) {
        int to = 0;
        double N_B = 0.0;
        double T_tau_t = 0.0;
        for (int ix = 0; ix <= nx; ix++) {
            for (int iy = 0; iy <= ny; iy++) {
                for (int ieta = 0; ieta < neta; ieta++) {
                    double eta_s = (deta*(ieta + neta*rank)
                                    - (DATA->eta_size)/2.0);
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
        MPI::COMM_WORLD.Send(&N_B, 1 , MPI::DOUBLE, to, 1);
        MPI::COMM_WORLD.Send(&T_tau_t, 1 , MPI::DOUBLE, to, 2);
    }
  
    if (rank == 0) {
        double N_B = 0.0;
        double T_tau_t = 0.0;
        for (int ix = 0; ix<=nx; ix++) {
            for(int iy = 0; iy <= ny; iy++) {
                for(int ieta = 0; ieta < neta; ieta++) {
                    double eta_s = (deta*(ieta + neta*rank)
                                    - (DATA->eta_size)/2.0);
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
        for (int irank = 1; irank < size; irank++) {
            int from = irank;
            double N_B_rank;
            double T_tau_t_rank;
            MPI::COMM_WORLD.Recv(&N_B_rank, 1, MPI::DOUBLE, from, 1);
            MPI::COMM_WORLD.Recv(&T_tau_t_rank, 1, MPI::DOUBLE, from, 2);
            N_B += N_B_rank;
            T_tau_t += T_tau_t_rank;
        }
        N_B *= tau*dx*dy*deta;
        T_tau_t *= tau*dx*dy*deta*0.19733;
        cout << "check: net baryon number N_B = " << N_B << endl;
        cout << "check: total energy T^{taut} = " << T_tau_t << " GeV" << endl;
	}
}


void Grid_info::Tmax_profile(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank)
{
  Util *util;
  util = new Util();

  int position;
  int ix, iy, ieta, nx, ny, neta;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;

  int sizeOfData = (nx+1)*(ny+1)*(neta);
  int to, from;
  
  double *eps;
  double *rhob;

  eps = new double[sizeOfData];
  rhob = new double[sizeOfData];

   if (rank>0)
    {
      to = 0;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  position = ieta+(neta*(ix + ((nx+1)*iy)));
		  eps[position] = arena[ix][iy][ieta].epsilon;
		  rhob[position] = arena[ix][iy][ieta].rhob;
		}
	    }
	}
    
      MPI::COMM_WORLD.Send(eps,sizeOfData,MPI::DOUBLE,to,1);
      MPI::COMM_WORLD.Send(rhob,sizeOfData,MPI::DOUBLE,to,2);
    }
  
  if (rank==0) 
    {
      double ***epsFrom;
      double ***rhobFrom;
      
      epsFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      rhobFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  epsFrom[ix][iy][ieta] = arena[ix][iy][ieta].epsilon;
		  rhobFrom[ix][iy][ieta] = arena[ix][iy][ieta].rhob;
		}
	    }
	}
	

      for (int irank=1; irank<size; irank++)
	{
	  from = irank;
	  MPI::COMM_WORLD.Recv(eps,sizeOfData,MPI::DOUBLE,from,1);
	  MPI::COMM_WORLD.Recv(rhob,sizeOfData,MPI::DOUBLE,from,2);
	  
	  for(ix=0; ix<=nx; ix++)
	    {
	      for(iy=0; iy<=ny; iy++)
		{
		  for(ieta=0; ieta<neta; ieta++)
		    {
		      position = ieta+(neta*(ix + ((nx+1)*iy)));
		      epsFrom[ix][iy][ieta+irank*neta] = eps[position];
		      rhobFrom[ix][iy][ieta+irank*neta] = rhob[position];
		    }
		}
	    }
	}

      double e;

      double * e_max_eta = util->vector_malloc(size*neta);
      for(ieta=0; ieta<size*neta; ieta++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ix=0; ix<=nx; ix++)
		{
		  e = epsFrom[ix][iy][ieta];
		  if (e > e_max_eta[ieta]) e_max_eta[ieta] = e;
		}
	    }
	}
      double * e_max_x = util->vector_malloc(nx+1);
      for(ix=0; ix<=nx; ix++)
      {
            for(ieta=0; ieta<size*neta; ieta++)
	    {
	        for(iy=0; iy<=ny; iy++)
	        {
		  e = epsFrom[ix][iy][ieta];
		  if (e > e_max_x[ix]) e_max_x[ix] = e;
		}
	    }
      }
      double * e_max_y = util->vector_malloc(ny+1);
      for(iy=0; iy<=ny; iy++)
      {
            for(ieta=0; ieta<size*neta; ieta++)
	    {
	        for(ix=0; ix<=nx; ix++)
	        {
		  e = epsFrom[ix][iy][ieta];
		  if (e > e_max_y[iy]) e_max_y[iy] = e;
		}
	    }
      }

	int Nfiles = 3;
	string strTmaxAll[] = {"T_max_eta.txt", "T_max_x.txt", "T_max_y.txt"};
	int NRec[] = {neta*size, nx+1, ny+1};
	double * arrRec;
	for (int ifile = 0; ifile < Nfiles; ifile++)
	{
		if (! ( util->fileExists( strTmaxAll[ifile] ) ) ) 
		{
			std::ofstream Tmax (strTmaxAll[ifile].c_str(), std::ofstream::out | std::fstream::app);
			Tmax
				<< "#tau0= " <<  DATA->tau0
				<< " nx= " << DATA->nx 
				<< " ny= " << DATA->ny 
				<< " neta= " << DATA->neta
				<< " x_size= " << DATA->x_size
				<< " y_size= " << DATA->y_size
				<< " eta_size= " << DATA->eta_size
			<< endl;
			Tmax.close();
		}
			
		std::ofstream Tmax (strTmaxAll[ifile].c_str(), std::ofstream::out | std::fstream::app);
		Tmax << std::setprecision (4) << setw(5);

		switch (ifile)
		{
			case 0: arrRec = e_max_eta; break;
			case 1: arrRec = e_max_x; break;
			case 2: arrRec = e_max_y; break;
		}
		for (int iRec = 0; iRec < NRec[ifile]; iRec++) 
		Tmax << eos->get_temperature(arrRec[iRec], 0.)*hbarc << "\t";
		Tmax << endl;
		Tmax.close();
	}
      
      delete [] e_max_x;
      delete [] e_max_y;
      delete [] e_max_eta;
      util->cube_free(epsFrom,nx+1,ny+1,size*neta);
      util->cube_free(rhobFrom,nx+1,ny+1,size*neta);
    }

  delete [] eps;
  delete [] rhob;
  delete(util);
  
}


void Grid_info::getAverageTandPlasmaEvolution(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank)
{
  Util *util;
  util = new Util();

  int position;
  int ix, iy, ieta, nx, ny, neta;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;

  /// MPI send and receive

  int sizeOfData = (nx+1)*(ny+1)*(neta);
  int to, from;
  
  double *eps;
  double *rhob;

  eps = new double[sizeOfData];
  rhob = new double[sizeOfData];

   if (rank>0) //send all to rank 0
    {
      to = 0;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  position = ieta+(neta*(ix + ((nx+1)*iy)));
		  eps[position] = arena[ix][iy][ieta].epsilon;
		  rhob[position] = arena[ix][iy][ieta].rhob;
		}
	    }
	}
      //      cout << "eps[5]=" <<eps[5]<< endl;
    
      MPI::COMM_WORLD.Send(eps,sizeOfData,MPI::DOUBLE,to,1);
      MPI::COMM_WORLD.Send(rhob,sizeOfData,MPI::DOUBLE,to,2);
    }
  
  if (rank==0) 
    {
      double ***epsFrom;
      double ***rhobFrom;
      
      epsFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      rhobFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  epsFrom[ix][iy][ieta] = arena[ix][iy][ieta].epsilon;
		  rhobFrom[ix][iy][ieta] = arena[ix][iy][ieta].rhob;
		}
	    }
	}
	

      for (int irank=1; irank<size; irank++)
	{
	  from = irank;
	  MPI::COMM_WORLD.Recv(eps,sizeOfData,MPI::DOUBLE,from,1);
	  MPI::COMM_WORLD.Recv(rhob,sizeOfData,MPI::DOUBLE,from,2);
	  
	  for(ix=0; ix<=nx; ix++)
	    {
	      for(iy=0; iy<=ny; iy++)
		{
		  for(ieta=0; ieta<neta; ieta++)
		    {
		      position = ieta+(neta*(ix + ((nx+1)*iy)));
		      epsFrom[ix][iy][ieta+irank*neta] = eps[position];
		      rhobFrom[ix][iy][ieta+irank*neta] = rhob[position];
		    }
		}
	    }
	  //	  cout << "epsFrom[0][0][irank*neta+5]=" <<epsFrom[0][0][irank*neta+5]<< endl;
	}
      
      //cout << " RECEIVED ALL DATA FROM OTHER ranks" << endl;

      neta*=size;
      ofstream out_file;
      string out_name = "avgT.dat";
      out_file.open(out_name.c_str(), ios::out | ios::app );
      ofstream out_file2;
      string out_name2 = "plasmaEvolutionTime.dat";
      out_file2.open(out_name2.c_str(), ios::out | ios::app );

//       int iz;
      int nCells=0, nCells2=0;
//       double x, y, z, eta, delta_z, z_size, etafrac;
      double T, eps;
      double totalT=0, totalT2=0;
//       double ux, uy, ueta, utau, epsilon;
      double rhob, QGPfrac;
//       double eta_lower, eps_lower, eps_higher, rhob_lower, rhob_higher;
      for(ieta=0; ieta<neta; ieta++)
	{
	  for(iy=0; iy<DATA->ny; iy++) // do all y
	    {
	      for(ix=0; ix<DATA->nx; ix++) // do all x
		{
		  eps = epsFrom[ix][iy][ieta];
		  rhob = rhobFrom[ix][iy][ieta];
			  
		  T=eos->get_temperature(eps,rhob);
		  //QGPfrac=eos->get_qgp_frac(eps,rhob);
              QGPfrac = 0.0;
		  
		  if(T>0.16/hbarc && QGPfrac==1)
		    {
		      totalT+=T;
		      nCells+=1;
		    }
	
		  if(T>0.16/hbarc)
		    {
		      totalT2+=T;
		      nCells2+=1;
		    }
		  
		}/* ix */
	    }/* iy */
	}/* iz */
	
      if(nCells==0 && DATA->plasmaEvolutionTime==0)
	{
	  DATA->plasmaEvolutionTime=tau-DATA->tau0;
	  out_file2 << setprecision(10) << tau-DATA->tau0;
     	}

      if(nCells2==0 && DATA->plasmaEvolutionTime2==0)
	{
	  DATA->plasmaEvolutionTime2=tau-DATA->tau0;
	  out_file2 << setprecision(10) << tau-DATA->tau0 << endl;
	}
      
      if (nCells>0)
	{
	  totalT/=nCells;
	  DATA->avgT+=totalT;
	  DATA->nSteps+=1;
	}

      if (nCells2>0)
	{
	  totalT2/=nCells2;
	  DATA->avgT2+=totalT2;
	  DATA->nSteps2+=1;
	}
      
      out_file << setprecision(10) << setw(10) << tau << " " << setw(10) <<  totalT*hbarc << " " << setw(10) << totalT2*hbarc << endl;
      
      out_file.close();
      out_file2.close();
      util->cube_free(epsFrom,nx+1,ny+1,size*neta);
      util->cube_free(rhobFrom,nx+1,ny+1,size*neta);
    }

//   delete [] eps;
//   delete [] rhob;
  delete(util);
  
}/* OutputEvolutionDataXYZ */



void Grid_info::OutputXY(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank)
{
  Util *util;
  util = new Util();

  int position;
  int ix, iy, ieta, nx, ny, neta;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;

  /// MPI send and receive

  int sizeOfData = (nx+1)*(ny+1)*(neta);
  int to, from;
  double s;
  
  double *trouble;
  double *eps;
  double *Txx;
  double *Tyy;
  double *Txy;
  double *rhob;
  double *utau;
  double *ux;
  double *uy;
  double *ueta;
  //  double *corr;
  
  trouble = (double *)malloc(sizeof(double)*sizeOfData);
  eps = (double *)malloc(sizeof(double)*sizeOfData);
  Txx = (double *)malloc(sizeof(double)*sizeOfData);
  Tyy = (double *)malloc(sizeof(double)*sizeOfData);
  Txy = (double *)malloc(sizeof(double)*sizeOfData);
  rhob = (double *)malloc(sizeof(double)*sizeOfData);
  utau = (double *)malloc(sizeof(double)*sizeOfData);
  ux = (double *)malloc(sizeof(double)*sizeOfData);
  uy = (double *)malloc(sizeof(double)*sizeOfData);
  ueta = (double *)malloc(sizeof(double)*sizeOfData);
  //  corr = (double *)malloc(sizeof(double)*sizeOfData);

  if (rank>0) //send all to rank 0
    {
      to = 0;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  position = ieta+(neta*(ix + ((nx+1)*iy)));
		  eps[position] = arena[ix][iy][ieta].epsilon;
		  trouble[position] = static_cast<double>(arena[ix][iy][ieta].trouble);
		  Txx[position] = arena[ix][iy][ieta].TJb[0][1][1];
		  Tyy[position] = arena[ix][iy][ieta].TJb[0][2][2];
		  Txy[position] = arena[ix][iy][ieta].TJb[0][1][2];
		  rhob[position] = arena[ix][iy][ieta].rhob;
		  utau[position] = arena[ix][iy][ieta].u[0][0];
		  ux[position] = arena[ix][iy][ieta].u[0][1];
		  uy[position] = arena[ix][iy][ieta].u[0][2];
		  ueta[position] = arena[ix][iy][ieta].u[0][3];
		  //corr[position] = arena[ix][iy][ieta].rhob;
		}
	    }
	}
      //      cout << "eps[5]=" <<eps[5]<< endl;
    
      MPI::COMM_WORLD.Send(trouble,sizeOfData,MPI::DOUBLE,to,0);
      MPI::COMM_WORLD.Send(eps,sizeOfData,MPI::DOUBLE,to,1);
      MPI::COMM_WORLD.Send(rhob,sizeOfData,MPI::DOUBLE,to,2);
      MPI::COMM_WORLD.Send(utau,sizeOfData,MPI::DOUBLE,to,3);
      MPI::COMM_WORLD.Send(ux,sizeOfData,MPI::DOUBLE,to,4);
      MPI::COMM_WORLD.Send(uy,sizeOfData,MPI::DOUBLE,to,5);
      MPI::COMM_WORLD.Send(ueta,sizeOfData,MPI::DOUBLE,to,6);
      MPI::COMM_WORLD.Send(Txx,sizeOfData,MPI::DOUBLE,to,7);
      MPI::COMM_WORLD.Send(Tyy,sizeOfData,MPI::DOUBLE,to,8);
      MPI::COMM_WORLD.Send(Txy,sizeOfData,MPI::DOUBLE,to,9);
      //    MPI::COMM_WORLD.Send(corr,sizeOfData,MPI::DOUBLE,to,10);
    }
  
  if (rank==0) 
    {
      double ***troubleFrom;
      double ***epsFrom;
      double ***TxxFrom;
      double ***TyyFrom;
      double ***TxyFrom;
      double ***rhobFrom;
      double ***utauFrom;
      double ***uxFrom;
      double ***uyFrom;
      double ***uetaFrom;
      //      double ***corrFrom;
      
      troubleFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      //corrFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      epsFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      TxxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      TyyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      TxyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      rhobFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      utauFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  troubleFrom[ix][iy][ieta] = static_cast<double>(arena[ix][iy][ieta].trouble);
		  //corrFrom[ix][iy][ieta] = arena[ix][iy][ieta].rhob;
		  epsFrom[ix][iy][ieta] = arena[ix][iy][ieta].epsilon;
		  TxxFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][1][1];
		  TyyFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][2][2];
		  TxyFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][1][2];
		  rhobFrom[ix][iy][ieta] = arena[ix][iy][ieta].rhob;
		  utauFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][0];
		  uxFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][1];
		  uyFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][2];
		  uetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][3];
		}
	    }
	}
	

      for (int irank=1; irank<size; irank++)
	{
	  from = irank;
	  MPI::COMM_WORLD.Recv(trouble,sizeOfData,MPI::DOUBLE,from,0);
	  MPI::COMM_WORLD.Recv(eps,sizeOfData,MPI::DOUBLE,from,1);
	  MPI::COMM_WORLD.Recv(rhob,sizeOfData,MPI::DOUBLE,from,2);
	  MPI::COMM_WORLD.Recv(utau,sizeOfData,MPI::DOUBLE,from,3);
	  MPI::COMM_WORLD.Recv(ux,sizeOfData,MPI::DOUBLE,from,4);
	  MPI::COMM_WORLD.Recv(uy,sizeOfData,MPI::DOUBLE,from,5);
	  MPI::COMM_WORLD.Recv(ueta,sizeOfData,MPI::DOUBLE,from,6);
	  MPI::COMM_WORLD.Recv(Txx,sizeOfData,MPI::DOUBLE,from,7);
	  MPI::COMM_WORLD.Recv(Tyy,sizeOfData,MPI::DOUBLE,from,8);
	  MPI::COMM_WORLD.Recv(Txy,sizeOfData,MPI::DOUBLE,from,9);
	  //	  MPI::COMM_WORLD.Recv(corr,sizeOfData,MPI::DOUBLE,from,10);
	  
	  for(ix=0; ix<=nx; ix++)
	    {
	      for(iy=0; iy<=ny; iy++)
		{
		  for(ieta=0; ieta<neta; ieta++)
		    {
		      position = ieta+(neta*(ix + ((nx+1)*iy)));
		      //	      corrFrom[ix][iy][ieta+irank*neta] = corr[position];
		      troubleFrom[ix][iy][ieta+irank*neta] = trouble[position];
		      epsFrom[ix][iy][ieta+irank*neta] = eps[position];
		      rhobFrom[ix][iy][ieta+irank*neta] = rhob[position];
		      utauFrom[ix][iy][ieta+irank*neta] = utau[position];
		      uxFrom[ix][iy][ieta+irank*neta] = ux[position];
		      uyFrom[ix][iy][ieta+irank*neta] = uy[position];
		      uetaFrom[ix][iy][ieta+irank*neta] = ueta[position];
		      TxxFrom[ix][iy][ieta+irank*neta] = Txx[position];
		      TyyFrom[ix][iy][ieta+irank*neta] = Tyy[position];
		      TxyFrom[ix][iy][ieta+irank*neta] = Txy[position];
		    }
		}
	    }
	  //	  cout << "epsFrom[0][0][irank*neta+5]=" <<epsFrom[0][0][irank*neta+5]<< endl;
	}
      
      //cout << " RECEIVED ALL DATA FROM OTHER ranks" << endl;
      FILE *ecc_file;
      const char* ecc_name = "eccentricity.dat";
      ecc_file = fopen(ecc_name, "a");
      
//       double value, epsp;
      double rA;
//       double numerator, denominator;
      double avcos, avsin, avcos3, avsin3;
      double Psi2, phiA, Psi3, avrSq, avr3;
      double eccentricity2, eccentricity3;

      FILE *out_file;
      const char* out_name = "e_x_y_profile.dat";
      out_file = fopen(out_name, "a");
//       fprintf(out_file,"");
      FILE *s_file;
      const char* s_name = "entropy-eta.dat";
      s_file = fopen(s_name, "a");
//       fprintf(s_file,"");
      FILE *out_file_2;
      const char* out_name_2 = "e_x_y_profile_05.dat";
      out_file_2 = fopen(out_name_2, "a");
//       fprintf(out_file_2,"");
      int iz, nz;
//       double trouble_lower;
      double T, x, y, z, eta, delta_z, z_size, eps, Txx, Tyy, Txy,etafrac;
//       double trouble, trouble_higher;
      //double corr_lower, corr_higher;
      double ux, uy, ueta, utau, rhob, QGPfrac;
//       double epsilon;
      double eta_lower, eps_lower, eps_higher, rhob_lower, rhob_higher, Txx_lower, Txx_higher, Tyy_lower, Tyy_higher, Txy_lower, Txy_higher;
      double utau_lower, utau_higher, ux_lower, ux_higher, uy_lower, uy_higher, ueta_lower, ueta_higher, u0, uz;
      nz = 80;
      delta_z = 0.5;
      z_size = 40.; // z_size is fixed to 40 fm: corresponds to Hydro:zmax=20 in MARTINI!
      //iz=nz/2+1;
      iz=nz/2;//nz-1;
      // get temperature value and others by interpolation
      z = iz*(delta_z) - (z_size)/2.0;
      eta = asinh(z/tau);
      //eta = 1.;
      //      eta+=5;
      //      cerr << "eta =" << eta << endl;
      //cout << "eta_size=" << (DATA->eta_size) << endl;
      ieta = floor((eta + (DATA->eta_size)/2.0)/(DATA->delta_eta));
      eta_lower = ieta*(DATA->delta_eta) - (DATA->eta_size)/2.0;
      etafrac = (eta-eta_lower)/DATA->delta_eta;
      
      for(iy=0; iy<DATA->ny; iy++) 
	{
	  y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
	  for(ix=0; ix<DATA->nx; ix++)
	    {
	      x = ix*(DATA->delta_x) - (DATA->x_size)/2.0;
// 	      trouble = troubleFrom[ix][iy][ieta];
	      //corr_lower = corrFrom[ix][iy][ieta];
	      eps_lower = epsFrom[ix][iy][ieta];
	      Txx_lower = TxxFrom[ix][iy][ieta];
	      Tyy_lower = TyyFrom[ix][iy][ieta];
	      Txy_lower = TxyFrom[ix][iy][ieta];
	      rhob_lower = rhobFrom[ix][iy][ieta];
	      utau_lower = utauFrom[ix][iy][ieta];
	      ux_lower = uxFrom[ix][iy][ieta];
	      uy_lower = uyFrom[ix][iy][ieta];
	      ueta_lower = uetaFrom[ix][iy][ieta];
	      
	      if (ieta+1 < DATA->neta*size)
		{
		  eps_higher = epsFrom[ix][iy][ieta+1];
// 		  trouble_higher = troubleFrom[ix][iy][ieta+1];
		  //corr_higher = corrFrom[ix][iy][ieta+1];
		  Txx_higher = TxxFrom[ix][iy][ieta+1];
		  Tyy_higher = TyyFrom[ix][iy][ieta+1];
		  Txy_higher = TxyFrom[ix][iy][ieta+1];
		  rhob_higher = rhobFrom[ix][iy][ieta+1];
		  utau_higher = utauFrom[ix][iy][ieta+1];
		  ux_higher = uxFrom[ix][iy][ieta+1];
		  uy_higher = uyFrom[ix][iy][ieta+1];
		  ueta_higher = uetaFrom[ix][iy][ieta+1];
		}
	      else
		{
		  eps_higher = eps_lower;
// 		  trouble_higher = trouble_lower;
		  //corr_higher = corr_lower;
		  Txx_higher = Txx_lower;
		  Tyy_higher = Tyy_lower;
		  Txy_higher = Txy_lower;
		  rhob_higher = rhob_lower;
		  utau_higher = utau_lower;
		  ux_higher = ux_lower;
		  uy_higher = uy_lower;
		  ueta_higher = ueta_lower;
		}
	      
	      eps = eps_lower * (1.-etafrac) + (etafrac)*eps_higher;
// 	      trouble = trouble_lower * (1.-etafrac) + (etafrac)*trouble_higher;
	      //corr = corr_lower * (1.-etafrac) + (etafrac)*corr_higher;
	      Txx = Txx_lower * (1.-etafrac) + (etafrac)*Txx_higher;
	      Tyy = Tyy_lower * (1.-etafrac) + (etafrac)*Tyy_higher;
	      Txy = Txy_lower * (1.-etafrac) + (etafrac)*Txy_higher;
	      rhob = rhob_lower * (1.-etafrac) + (etafrac)*rhob_higher;
	      utau = utau_lower * (1.-etafrac) + (etafrac)*utau_higher;
	      ux = ux_lower * (1.-etafrac) + (etafrac)*ux_higher;
	      uy = uy_lower * (1.-etafrac) + (etafrac)*uy_higher;
	      ueta = ueta_lower * (1.-etafrac) + (etafrac)*ueta_higher;
	      
	      T = eos->get_temperature(eps, rhob);
	      //QGPfrac = eos->get_qgp_frac(eps, rhob);
            QGPfrac = 0.0;
	      
	      //Now these are the flow velocities as e,g, MARTINI needs them
	      u0 = ueta*sinh(eta)+utau*cosh(eta); // = gamma factor
	      ux = ux/u0;
	      uy = uy/u0;
	      uz = ueta*cosh(eta)+utau*sinh(eta);
	      uz /= u0;
	      
	      /* 	     epsilon = arena[ix][iy][ieta].epsilon; */
	      /* 	     rhob = arena[ix][iy][ieta].rhob; */
	      /* 	     utau = arena[ix][iy][ieta].u[0][0]; */
	      /* 	     ux = arena[ix][iy][ieta].u[0][1]; */
	      /* 	     uy = arena[ix][iy][ieta].u[0][2]; */
	      /* 	     ueta = arena[ix][iy][ieta].u[0][3]; */
	      
	      
	      // output with coordinates for testing:
	      //fprintf(out_file,"%e %e %e %e %e %e %e %e %e %e\n", 
	      //  tau, x, y, z, eta, T*hbarc, QGPfrac, ux, uy, uz);
	      //cout << "(tau,x,y,z)=(" << tau << ", " << x << ", " << y << ", " << z << ")" << endl;
	      

	      //	      if(ix%5 == 0 && iy%5 == 0)
	      fprintf(out_file,"%e %e %e %e %e %e %e %e %e %e %e %e \n",tau, x, y, 
		      T*hbarc, QGPfrac, eps, Txx, Tyy, Txy, ux, uy, ueta);
	      
	      if (ix==DATA->nx-1)
	      //if (ix==DATA->nx-1 && iy%5==0)
		fprintf(out_file,"\n");
	    }/*  */
	}/* iy */
      fclose(out_file);

      // compute and print eccentricity and angles:
//       numerator = 0.;
//       denominator = 0.;
      avcos = 0.;
      avsin = 0.;
      avcos3 = 0.;
      avsin3 = 0.;
      avrSq=0.;
      avr3=0.;
      ieta = DATA->neta/2;
      for(ix=0; ix<=DATA->nx; ix++)
	{
	  x = ix*(DATA->delta_x) - (DATA->x_size/2.0);  
	  for(iy=0; iy<=DATA->ny; iy++)
	    {
	      y = iy*(DATA->delta_x) - (DATA->y_size/2.0);
	      if (x>=0)
		{
		  phiA = atan(y/x);
		  if (x==0) 
		    {
		      if (y>=0) phiA=PI/2.;
		      else if (y<0) phiA=3.*PI/2.;
		    }
		}
	      else
		{
		  phiA = atan(y/x)+PI;
		}
	      rA = sqrt( x*x + y*y );
	      avrSq += rA*rA*(arena[ix][iy][ieta].epsilon); // compute average r^2
	      avr3 += rA*rA*rA*(arena[ix][iy][ieta].epsilon);
	      
	      avcos  += rA*rA*cos(2.*phiA)*(arena[ix][iy][ieta].epsilon);
	      avsin  += rA*rA*sin(2.*phiA)*(arena[ix][iy][ieta].epsilon);
	      avcos3 += rA*rA*rA*cos(3.*phiA)*(arena[ix][iy][ieta].epsilon);
	      avsin3 += rA*rA*rA*sin(3.*phiA)*(arena[ix][iy][ieta].epsilon);
	    }
	}
      Psi2 = (atan(avsin/avcos))/2.+PI/2.;
      Psi3 = (atan(avsin3/avcos3)+PI)/3.;
      eccentricity2 = sqrt(avcos*avcos+avsin*avsin)/avrSq;
      eccentricity3 = sqrt(avcos3*avcos3+avsin3*avsin3)/avr3;
      
      cout << "ecc2=" << eccentricity2 << endl;
      cout << "Psi2=" << Psi2 << endl;
      cout << "ecc3=" << eccentricity3 << endl;
      cout << "Psi3=" << Psi3 << endl;
      fprintf(ecc_file, "%f %f %f %f %f \n", tau, eccentricity2, Psi2, eccentricity3, Psi3);
      
      fclose(ecc_file);


      // compute entropy as a function of rapidity
      for(ieta=0; ieta<DATA->neta*size; ieta++) 
	{
	  s = 0.;
	  eta = ieta*(DATA->delta_eta) - (DATA->eta_size)/2.0;
	  for(iy=0; iy<DATA->ny; iy++) 
	    {
	      y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
	      for(ix=0; ix<DATA->nx; ix++)
		{
		  eps = epsFrom[ix][iy][ieta];
		  rhob=0.;
		  if (DATA->whichEOS==1)
		    {
		      cout << " not implemented for EOS=1" << endl;
		    }
		  else if (DATA->whichEOS>1)
		    {
		      s += eos->interpolate2(eps, rhob, 2)*utauFrom[ix][iy][ieta]*tau;
		    }
		}/*  */
	    }/* iy */
	  fprintf(s_file,"%e %e %e \n", tau, eta, s);
	}//ieta
      fprintf(s_file,"\n");
      fclose(s_file);

//       eta = 0.5;
//       //      eta+=5;
//       cerr << "eta =" << eta << endl;
//       //cout << "eta_size=" << (DATA->eta_size) << endl;
//       ieta = floor((eta + (DATA->eta_size)/2.0)/(DATA->delta_eta));
//       eta_lower = ieta*(DATA->delta_eta) - (DATA->eta_size)/2.0;
//       etafrac = (eta-eta_lower)/DATA->delta_eta;
      
//       for(iy=0; iy<DATA->ny; iy++) 
// 	{
// 	  y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
// 	  for(ix=0; ix<DATA->nx; ix++)
// 	    {
// 	      x = ix*(DATA->delta_x) - (DATA->x_size)/2.0;
// 	      trouble = troubleFrom[ix][iy][ieta];
// 	      eps_lower = epsFrom[ix][iy][ieta];
// 	      Txx_lower = TxxFrom[ix][iy][ieta];
// 	      Tyy_lower = TyyFrom[ix][iy][ieta];
// 	      Txy_lower = TxyFrom[ix][iy][ieta];
// 	      rhob_lower = rhobFrom[ix][iy][ieta];
// 	      utau_lower = utauFrom[ix][iy][ieta];
// 	      ux_lower = uxFrom[ix][iy][ieta];
// 	      uy_lower = uyFrom[ix][iy][ieta];
// 	      ueta_lower = uetaFrom[ix][iy][ieta];
	      
// 	      if (ieta+1 < DATA->neta*size)
// 		{
// 		  eps_higher = epsFrom[ix][iy][ieta+1];
// 		  trouble_higher = troubleFrom[ix][iy][ieta+1];
// 		  Txx_higher = TxxFrom[ix][iy][ieta+1];
// 		  Tyy_higher = TyyFrom[ix][iy][ieta+1];
// 		  Txy_higher = TxyFrom[ix][iy][ieta+1];
// 		  rhob_higher = rhobFrom[ix][iy][ieta+1];
// 		  utau_higher = utauFrom[ix][iy][ieta+1];
// 		  ux_higher = uxFrom[ix][iy][ieta+1];
// 		  uy_higher = uyFrom[ix][iy][ieta+1];
// 		  ueta_higher = uetaFrom[ix][iy][ieta+1];
// 		}
// 	      else
// 		{
// 		  eps_higher = eps_lower;
// 		  trouble_higher = trouble_lower;
// 		  Txx_higher = Txx_lower;
// 		  Tyy_higher = Tyy_lower;
// 		  Txy_higher = Txy_lower;
// 		  rhob_higher = rhob_lower;
// 		  utau_higher = utau_lower;
// 		  ux_higher = ux_lower;
// 		  uy_higher = uy_lower;
// 		  ueta_higher = ueta_lower;
// 		}
	      
// 	      eps = eps_lower * (1.-etafrac) + (etafrac)*eps_higher;
// 	      trouble = trouble_lower * (1.-etafrac) + (etafrac)*trouble_higher;
// 	      Txx = Txx_lower * (1.-etafrac) + (etafrac)*Txx_higher;
// 	      Tyy = Tyy_lower * (1.-etafrac) + (etafrac)*Tyy_higher;
// 	      Txy = Txy_lower * (1.-etafrac) + (etafrac)*Txy_higher;
// 	      rhob = rhob_lower * (1.-etafrac) + (etafrac)*rhob_higher;
// 	      utau = utau_lower * (1.-etafrac) + (etafrac)*utau_higher;
// 	      ux = ux_lower * (1.-etafrac) + (etafrac)*ux_higher;
// 	      uy = uy_lower * (1.-etafrac) + (etafrac)*uy_higher;
// 	      ueta = ueta_lower * (1.-etafrac) + (etafrac)*ueta_higher;
	      
// 	      if (DATA->whichEOS==1)
// 		{
// 		  T = eos->interpolate(eps, rhob, 0);
// 		  QGPfrac=(eps*hbarc-0.45)/(1.6-0.45); // e=(1-QGPfrac)*e_H + QGPfrac*e_QGP in the mixed phase
// 		  //cout << T << " " << QGPfrac << endl;
// 		  if (QGPfrac>1.) QGPfrac = 1;
// 		  else if (QGPfrac<0.) QGPfrac=0.;
// 		}
// 	      else if (DATA->whichEOS>1)
// 		{
// 		  T = eos->interpolate2(eps, rhob, 1);
// 		  QGPfrac = eos->interpolate2(eps, rhob, 3);
// 		}
	      
// 	      //Now these are the flow velocities as e,g, MARTINI needs them
// 	      u0 = ueta*sinh(eta)+utau*cosh(eta); // = gamma factor
// 	      ux = ux/u0;
// 	      uy = uy/u0;
// 	      uz = ueta*cosh(eta)+utau*sinh(eta);
// 	      uz /= u0;
	      
// 	      /* 	     epsilon = arena[ix][iy][ieta].epsilon; */
// 	      /* 	     rhob = arena[ix][iy][ieta].rhob; */
// 	      /* 	     utau = arena[ix][iy][ieta].u[0][0]; */
// 	      /* 	     ux = arena[ix][iy][ieta].u[0][1]; */
// 	      /* 	     uy = arena[ix][iy][ieta].u[0][2]; */
// 	      /* 	     ueta = arena[ix][iy][ieta].u[0][3]; */
	      
	      
// 	      // output with coordinates for testing:
// 	      //fprintf(out_file,"%e %e %e %e %e %e %e %e %e %e\n", 
// 	      //  tau, x, y, z, eta, T*hbarc, QGPfrac, ux, uy, uz);
// 	      //cout << "(tau,x,y,z)=(" << tau << ", " << x << ", " << y << ", " << z << ")" << endl;
	      

// 	      if(ix%5 == 0 && iy%5 == 0)
// 		fprintf(out_file_2,"%e %e %e %e %e %e %e %e %e \n",tau, x, y, 
// 			T*hbarc, QGPfrac, eps, Txx, Tyy, Txy);
	      
// 	      if (ix==DATA->nx-1 && iy%5==0)
// 		fprintf(out_file_2,"\n");
// 	    }/* ix */
// 	}/* iy */
      fclose(out_file_2);

      util->cube_free(troubleFrom,nx+1,ny+1,size*neta);
      //     util->cube_free(corrFrom,nx+1,ny+1,size*neta);
      util->cube_free(epsFrom,nx+1,ny+1,size*neta);
      util->cube_free(TxxFrom,nx+1,ny+1,size*neta);
      util->cube_free(TyyFrom,nx+1,ny+1,size*neta);
      util->cube_free(TxyFrom,nx+1,ny+1,size*neta);
      util->cube_free(rhobFrom,nx+1,ny+1,size*neta);
      util->cube_free(utauFrom,nx+1,ny+1,size*neta);
      util->cube_free(uxFrom,nx+1,ny+1,size*neta);
      util->cube_free(uyFrom,nx+1,ny+1,size*neta);
      util->cube_free(uetaFrom,nx+1,ny+1,size*neta);
    }
  free(trouble);
  //  free(corr);
  free(eps);
  free(Txx);
  free(Tyy);
  free(Txy);
  free(rhob);
  free(utau);
  free(ux);
  free(uy);
  free(ueta);
  delete(util);
}/* OutputXY */


void Grid_info::PrintArena(Grid ***arena, InitData *DATA, double tau)
{
 int ix, iy, ieta;
 double x, y, eta;
 double ux, uy, ueta, utau, epsilon, rhob;
 for(ix=0; ix<=DATA->nx; ix++)
  {
   x = ix*(DATA->delta_x) - (DATA->x_size)/2.0;
   for(iy=0; iy<=DATA->ny; iy++)
    {
     y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
     for(ieta=0; ieta<DATA->neta; ieta++)
      {
       eta = ieta*(DATA->delta_eta) - (DATA->eta_size)/2.0;
       epsilon = arena[ix][iy][ieta].epsilon;
       rhob = arena[ix][iy][ieta].rhob;
       utau = arena[ix][iy][ieta].u[0][0];
       ux = arena[ix][iy][ieta].u[0][1];
       uy = arena[ix][iy][ieta].u[0][2];
       ueta = arena[ix][iy][ieta].u[0][3];

       printf("%e  %e  %e  %e  %e  %e  %e  %e  %e\n", 
               tau, x, y, eta, epsilon, rhob, ux/utau, uy/utau, ueta/utau);
      }/* ieta */
    }/* iy */
  }/* ix */
}/* PrintArena */

void Grid_info::PrintEtaEpsilon(Grid ***arena, InitData *DATA, double tau, int size, int rank)
{
 int ix, iy, ieta;
//  double x, y;
 double eta;
//  double ux, uy, ueta, utau, rhob;
 double epsilon;
 string d_name;
 ofstream d_file;

 if(rank==0) d_name = "e_profile.dat";
 else d_name = "e_profile2.dat";
 d_file.open(d_name.c_str(), ios::out | ios::app );

//  x=0.;
//  y=0.;
 ix = floor(DATA->x_size/2/DATA->delta_x);
 iy = floor(DATA->y_size/2/DATA->delta_y);

 for(ieta=0; ieta<DATA->neta; ieta++)
   {
     eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
     //fprintf(stderr,"%d  %d  %d \n", ix, iy, ieta);
     epsilon = arena[ix][iy][ieta].epsilon;
     d_file << tau << " " << eta << " " << setprecision(8) <<epsilon*hbarc << endl;
   }

 d_file << endl;
 d_file.close();
}/* PrintEtaEpsilon */

void Grid_info::PrintxEpsilon(Grid ***arena, InitData *DATA, double tau, int size, int rank)
{
 int ix, iy, ieta;
 double x, y;
//  double ux, uy, ueta, utau, rhob, eta;
 double epsilon;
 if(size==1 || (size==2 && rank==1))
   {
     string d_name;
     ofstream d_file;
     string d2_name;
     ofstream d2_file;

     d_name = "e_x_profile.dat";
     d_file.open(d_name.c_str(), ios::out | ios::app );
     d2_name = "e_y_profile.dat";
     d2_file.open(d2_name.c_str(), ios::out | ios::app );

     if (size==1) ieta = floor(DATA->eta_size/2/DATA->delta_eta+0.00001);
     else if (size==2) ieta = 0;
     iy = floor(DATA->y_size/2/DATA->delta_y);
     
     for(ix=0; ix<=DATA->nx; ix++)
       {
	 x = ix*(DATA->delta_x) - (DATA->x_size)/2.0;
	 //fprintf(stderr,"%d  %d  %d \n", ix, iy, ieta);
	 epsilon = arena[ix][iy][ieta].epsilon;
	 d_file << setw(10) << tau << " " << setw(10) << x << " " << setw(10) << setprecision(8) <<epsilon*hbarc << endl;
       }
     d_file << endl; 
     d_file.close();
     
     ix = floor(DATA->x_size/2/DATA->delta_x);
     
     for(iy=0; iy<=DATA->ny; iy++)
       {
	 y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
	 //fprintf(stderr,"%d  %d  %d \n", ix, iy, ieta);
	 epsilon = arena[ix][iy][ieta].epsilon;
	 d_file << setw(10) << tau << " " << setw(10) << y << " " << setw(10) << setprecision(8) <<epsilon*hbarc << endl;
       }
     d2_file << endl; 
     d2_file.close();
   }
}/* PrintxEpsilon */


void Grid_info::PrintGrid(Grid *grid_p, int rk_order)
{
 int r, i, a;
//  int m, n, j;

 fprintf(stderr, "\nStarting PrintGrid...\n");
 fprintf(stderr, "TJb\n");
 for(r=0; r<rk_order; r++)
  {
   for(a=0; a<5; a++)
    {
     for(i=0; i<4; i++)
      {
       fprintf(stderr, "TJb[%d][%d][%d] = %e\n", r, a, i, grid_p->TJb[r][a][i]);
      }/* i */
    }/* a */
  }/* r */ 
 fprintf(stderr, "\n");
 
 fprintf(stderr, "u_mu\n");
 for(r=0; r<rk_order; r++)
  {
     for(i=0; i<4; i++)
      {
       fprintf(stderr, "u[%d][%d] = %e\n", r, i, grid_p->u[r][i]);
      }/* i */
  }/* r */ 
 fprintf(stderr, "\n");

 fprintf(stderr, "p = %e\n", grid_p->p);
 fprintf(stderr, "epsilon = %e\n", grid_p->epsilon);
 fprintf(stderr, "rhob = %e\n", grid_p->rhob);
 fprintf(stderr, "Done PrintGrid\n\n");
}/* PrintGrid */

void Grid_info::PrintAxy2(InitData *DATA, Grid ***arena, double tau)
{
 int ix, iy, ieta;
 double x, y, fxx, fyy, gxx, gyy, f00, g00, hxx, hyy, h00, kxx, kyy, k00;

 ieta = DATA->neta-1; //if i want midrapidity when using 2 processors
 for(ix=0; ix<=DATA->nx; ix++)
  {
   x = ix*(DATA->delta_x) - (DATA->x_size/2.0);  
   for(iy=0; iy<=DATA->ny; iy++)
    {
     y = iy*(DATA->delta_x) - (DATA->y_size/2.0);  
/* epsilon is in fm^-4 */
     gxx = arena[ix][iy][ieta].TJb[0][1][1];
     fxx = arena[ix][iy][ieta].Wmunu[0][1][1];
     hxx = arena[ix][iy][ieta].Pimunu[0][1][1];
     kxx = gxx + fxx + hxx;
     
     gyy = arena[ix][iy][ieta].TJb[0][2][2];
     fyy = arena[ix][iy][ieta].Wmunu[0][2][2];
     hyy = arena[ix][iy][ieta].Pimunu[0][2][2];
     kyy = gyy + fyy + hyy;
    
     g00 = arena[ix][iy][ieta].TJb[0][0][0];
     f00 = arena[ix][iy][ieta].Wmunu[0][0][0];
     h00 = arena[ix][iy][ieta].Pimunu[0][0][0];
     k00 = g00 + f00 + h00;
     
     printf(
     "%e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n", 
     x, y, 
     gxx, fxx, hxx, kxx,
     gyy, fyy, hyy, kyy,
     g00, f00, h00, k00,
     arena[ix][iy][ieta].epsilon);
    }
    printf("\n");
  }
  return;
}/* PrintAxy2 */


void Grid_info::PrintAxy(InitData *DATA, Grid ***arena, double tau)
{
 FILE *d_file;
 const char* d_name = "e_x_y_profile.dat";
 d_file = fopen(d_name, "a");
 int ix, iy, ieta;
 double x, y;

 ieta = DATA->neta/2;
 for(ix=0; ix<=DATA->nx; ix++)
  {
   x = ix*(DATA->delta_x) - (DATA->x_size/2.0);  
   for(iy=0; iy<=DATA->ny; iy++)
    {
     y = iy*(DATA->delta_y) - (DATA->y_size/2.0);  
/* epsilon is in fm^-4 */
     fprintf(d_file, "%e  %e  %e  %e  %e %e %e\n", 
               x, 
	       y, 
	       arena[ix][iy][ieta].TJb[0][1][1], 
	       arena[ix][iy][ieta].TJb[0][2][2],
	       arena[ix][iy][ieta].epsilon,
	       arena[ix][iy][ieta].T,
	       arena[ix][iy][ieta].mu
	    );
     if(iy==DATA->ny)
       fprintf(d_file,"\n");
    }
  }
 fclose(d_file);
  return;
}/* PrintAxy */

void Grid_info::print_qmu_evolution(InitData *DATA, Grid ***arena, double tau, EOS *eos, int rank)
{
  string d_name;
  ofstream d_file;
  d_name = "qmu_evo.dat";
  d_file.open(d_name.c_str(), ios::out | ios::app );
  
  int ieta = DATA->neta/2;
  double eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
  int iy = (int)(DATA->ny/2);
  for(int ix=0; ix<=DATA->nx; ix++)
  {
    double x = ix*(DATA->delta_x) - (DATA->x_size/2.0);
    double y = iy*(DATA->delta_x) - (DATA->y_size/2.0);
    double rhob = arena[ix][iy][ieta].rhob;
    double epsilon = arena[ix][iy][ieta].epsilon;
    double muB = eos->get_mu(epsilon, rhob);
    double kappa = 0.2*rhob/muB;
    d_file << scientific << setw(18) << setprecision(8)
           << tau << "  " << x << "  " << y << "  " << eta << "  " 
           << rhob << "  " ;
    double D_mu = arena[ix][iy][ieta].a[0][4];
    for(int nu = 0; nu < 4; nu++)
    {
        double q_nu = arena[ix][iy][ieta].Wmunu[0][4][nu];
        double partial_mu = arena[ix][iy][ieta].dUsup[0][4][nu];
        double flow_u = arena[ix][iy][ieta].u[0][nu];
        double NS_term = -kappa*(partial_mu + flow_u*D_mu);
        d_file << scientific << setw(18) << setprecision(8)
               << q_nu << "  " << NS_term << "  " ;
               //<< partial_mu << "  " << flow_u*D_mu << "  " << NS_term << "  " ;
    }
    d_file << endl;
  }
  int ix = (int)(DATA->ny/2);
  for(int iy=0; iy<=DATA->ny; iy++)
  {
    double x = ix*(DATA->delta_x) - (DATA->x_size/2.0);
    double y = iy*(DATA->delta_x) - (DATA->y_size/2.0);
    double rhob = arena[ix][iy][ieta].rhob;
    double epsilon = arena[ix][iy][ieta].epsilon;
    double muB = eos->get_mu(epsilon, rhob);
    double kappa = 0.2*rhob/muB;
    d_file << scientific << setw(18) << setprecision(8)
           << tau << "  " << x << "  " << y << "  " << eta << "  "
           << rhob << "  " ;
    double D_mu = arena[ix][iy][ieta].a[0][4];
    for(int nu = 0; nu < 4; nu++)
    {
        double q_nu = arena[ix][iy][ieta].Wmunu[0][4][nu];
        double partial_mu = arena[ix][iy][ieta].dUsup[0][4][nu];
        double flow_u = arena[ix][iy][ieta].u[0][nu];
        double NS_term = -kappa*(partial_mu + flow_u*D_mu);
        d_file << scientific << setw(18) << setprecision(8)
               << q_nu << "  " << NS_term << "  " ;
               //<< partial_mu << "  " << flow_u*D_mu << "   " ;
    }
    d_file << endl;
  }
  d_file.close();
  
  return;
}/* print qmu evolution */

void Grid_info::print_rhob_evolution(InitData *DATA, Grid ***arena, double tau, EOS* eos, int rank)
{
  string d_name;
  ofstream d_file;
  d_name = "rhoB_evo.dat";
  d_file.open(d_name.c_str(), ios::out | ios::app );
  
  int ieta = DATA->neta/2;
  double eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
  for(int ix=0; ix<=DATA->nx; ix++)
  {
    double x = ix*(DATA->delta_x) - (DATA->x_size/2.0);
    int iy = (int)(DATA->ny/2);
    double y = iy*(DATA->delta_x) - (DATA->y_size/2.0);
    double rhob = arena[ix][iy][ieta].rhob;
    double ed = arena[ix][iy][ieta].epsilon;
    double mub = eos->get_mu(ed, rhob);
    double temperature = eos->get_temperature(ed, rhob);
    d_file << scientific << setw(18) << setprecision(8)
           << tau << "  " << x << "  " << y << "  " << eta << "  "
           << ed << "  " << rhob << "  " 
           << mub << "  " << temperature << "  "
           << mub/temperature << endl;
  }
  for(int iy=0; iy<=DATA->ny; iy++)
  {
    double y = iy*(DATA->delta_x) - (DATA->y_size/2.0);
    int ix = (int)(DATA->nx/2);
    double x = ix*(DATA->delta_x) - (DATA->x_size/2.0);
    double rhob = arena[ix][iy][ieta].rhob;
    double ed = arena[ix][iy][ieta].epsilon;
    double mub = eos->get_mu(ed, rhob);
    double temperature = eos->get_temperature(ed, rhob);
    d_file << scientific << setw(18) << setprecision(8)
           << tau << "  " << x << "  " << y << "  "  << eta << "  "
           << ed << "  " << rhob << "  "
           << mub << "  " << temperature << "  "
           << mub/temperature << endl;
  }
  d_file.close();
  
  return;
}/* print rhob evolution */

void Grid_info::print_rhob_evolution_3d(InitData *DATA, Grid ***arena, double tau, EOS* eos, int rank)
{
  int n_skip = 5;

  string d_name;
  ofstream d_file;
  d_name = "rhoB_evo_3d_xy.dat";
  d_file.open(d_name.c_str(), ios::out | ios::app );

  // first in x-y plane at eta = 0
  double eta_goal = 0.0;
  if((int)((tau - DATA->tau0)/DATA->delta_tau) % n_skip == 0)
  {
     for(int ieta = 0; ieta < DATA->neta; ieta++)
     {
        double eta = DATA->delta_eta*(ieta + DATA->neta*rank) - DATA->eta_size/2.0;
        if(fabs(eta - eta_goal) < 1e-5)  // mid rapidity
        {
           for(int ix=0; ix<=DATA->nx; ix++)
           {
              double x = ix*(DATA->delta_x) - (DATA->x_size/2.0);
              for(int iy = 0; iy <= DATA->ny; iy++)
              {
                 double y = iy*(DATA->delta_y) - (DATA->y_size/2.0);
                 double rhob = arena[ix][iy][ieta].rhob;
                 double ed = arena[ix][iy][ieta].epsilon;
                 double mub = eos->get_mu(ed, rhob);
                 double temperature = eos->get_temperature(ed, rhob);
                 d_file << scientific << setw(18) << setprecision(8)
                        << tau << "  " << x << "  " << y << "  "
                        << ed << "  " << rhob << "  " 
                        << mub << "  " << temperature 
                        << endl;
              }
           }
        }
     }
  }
  d_file.close();

  ofstream d_file_eta;
  d_name = "rhoB_evo_3d_xeta.dat";
  d_file_eta.open(d_name.c_str(), ios::out | ios::app );
  // then in x-eta plane at y = 0
  if((int)((tau - DATA->tau0)/DATA->delta_tau) % n_skip == 0)
  {
     for(int ieta = 0; ieta < DATA->neta; ieta++)
     {
        double eta = DATA->delta_eta*(ieta + DATA->neta*rank) - DATA->eta_size/2.0;
        for(int ix=0; ix<=DATA->nx; ix++)
        {
           double x = ix*(DATA->delta_x) - (DATA->x_size/2.0);
           int iy = DATA->ny/2;
           double rhob = arena[ix][iy][ieta].rhob;
           double ed = arena[ix][iy][ieta].epsilon;
           double mub = eos->get_mu(ed, rhob);
           double temperature = eos->get_temperature(ed, rhob);
           d_file_eta << scientific << setw(18) << setprecision(8)
                      << tau << "  " << x << "  " << eta << "  "
                      << ed << "  " << rhob << "  " 
                      << mub << "  " << temperature 
                      << endl;
        }
     }
  }
  d_file_eta.close();
  
  return;
}/* print rhob evolution */

void Grid_info::print_fireball_evolution_on_phasediagram(InitData *DATA, Grid ***arena, double tau, EOS* eos, int rank)
{
  int n_skip = 25;

  int nT = 401;
  int nmu = 801;
  double T0 = 0.0;
  double T_max = 0.4;
  double dT = (T_max - T0)/(nT - 1);
  double mu0 = 0.0;
  double mu_max = 0.8;
  double dmu = (mu_max - mu0)/(nmu - 1);

  double eta_goal = 0.0;
  double e_dec = 0.1;   // GeV/fm^3
  double volume_element = DATA->delta_x * DATA->delta_y;

  if((int)((tau - DATA->tau0)/DATA->delta_tau) % n_skip == 0)
  {
     for(int ieta = 0; ieta < DATA->neta; ieta++)
     {
        double eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
        if(fabs(eta - eta_goal) < 1e-5)  // mid rapidity
        {
            double** dVdeta = new double* [nT];
            for(int iT = 0; iT < nT; iT++)
            {
               dVdeta[iT] = new double [nmu];
               for(int imu = 0; imu < nmu; imu++)
                  dVdeta[iT][imu] = 0.0;
            }

            ostringstream filename;
            filename << "fireball_evo_phasediagram_tau_" << tau << ".dat";
            ofstream d_file;
            d_file.open(filename.str().c_str());

            for(int ix=0; ix<=DATA->nx; ix++)
            {
                for(int iy = 0; iy <= DATA->ny; iy++)
                {
                    double rhob = arena[ix][iy][ieta].rhob;
                    double ed = arena[ix][iy][ieta].epsilon;
                    if(ed*hbarc > e_dec)
                    {
                       double mub = eos->get_mu(ed, rhob)*hbarc;   // GeV
                       double temperature = eos->get_temperature(ed, rhob)*hbarc;  // GeV

                       int T_idx = (int)((temperature - T0)/dT);
                       int mu_idx = (int)((mub - mu0)/dmu);
                       
                       if(T_idx >= 0 && T_idx < nT)
                          if(mu_idx >=0 && mu_idx < nmu)
                             dVdeta[T_idx][mu_idx] += volume_element;
                    }
               }
            }

            for(int iT = 0; iT < nT; iT++)
            {
               for(int imu = 0; imu < nmu; imu++)
               {
                   d_file << scientific << setw(18) << setprecision(8)
                          << dVdeta[iT][imu] << "   ";
               }
               d_file << endl;
            }
            d_file.close();

            for(int iT = 0; iT < nT; iT++)
                delete [] dVdeta[iT];
            delete [] dVdeta;
        }
     }
  }
  return;
}

void Grid_info::ComputeAnisotropy(InitData *DATA, Grid ***arena, double tau)
{
  
  string d_name;
  ofstream d_file;
  d_name = "aniso.dat";
  d_file.open(d_name.c_str(), ios::out | ios::app );
  
  int ix, iy, ieta;
  double v2;
//   double x, y, epsp;
  double numerator, denominator;
  numerator = 0.;
  denominator = 0.;
  
  ieta = DATA->neta/2;
  for(ix=0; ix<=DATA->nx; ix++)
    {
//       x = ix*(DATA->delta_x) - (DATA->x_size/2.0);  
      for(iy=0; iy<=DATA->ny; iy++)
	{
// 	  y = iy*(DATA->delta_x) - (DATA->y_size/2.0);
	  numerator += (arena[ix][iy][ieta].TJb[0][1][1]-arena[ix][iy][ieta].TJb[0][2][2]);
	  denominator += (arena[ix][iy][ieta].TJb[0][1][1]+arena[ix][iy][ieta].TJb[0][2][2]);
	}
    }
  v2 = numerator/denominator;
  d_file << setw(10) << tau << " " << setw(10) << setprecision(10) << v2 << endl;
  d_file.close();
  return;
}/* ComputeV2 */

void Grid_info::ComputeEccentricity(InitData *DATA, Grid ***arena, double tau)
{
  FILE *ecc_file;
  const char* ecc_name = "eccentricity.dat";
  ecc_file = fopen(ecc_name, "a");
  
  int ix, iy, ieta;
  double value, x, y;
//   double epsp;
  double numerator, denominator;
  numerator = 0.;
  denominator = 0.;
  
  ieta = DATA->neta/2;
  for(ix=0; ix<=DATA->nx; ix++)
    {
      x = ix*(DATA->delta_x) - (DATA->x_size/2.0);  
      for(iy=0; iy<=DATA->ny; iy++)
	{
	  y = iy*(DATA->delta_x) - (DATA->y_size/2.0);
	  numerator += (arena[ix][iy][ieta].epsilon)*(y*y-x*x);
	  denominator += (arena[ix][iy][ieta].epsilon)*(y*y+x*x);
	}
    }
  value = numerator/denominator;
  fprintf(ecc_file, "%f %f \n", tau, value);
  fclose(ecc_file);
  return;
}/* ComputeEccentricity */


void Grid_info::ComputeEnergyConservation(InitData *DATA, Grid ***arena, double tau)
{
  FILE *ecc_file;
  const char* ecc_name = "energyConservation.dat";
  ecc_file = fopen(ecc_name, "a");
  
  int ix, iy, ieta;
  double value;
//   double x, y, epsp, eta;
  double TtautauPart, TetaetaPart;//, TetaetaPrev;
  TtautauPart = 0.;
  TetaetaPart = 0.;
//   TetaetaPrev = 0.;
  
  for(ieta=0; ieta<DATA->neta; ieta++)
    {
      for(ix=0; ix<=DATA->nx; ix++)
	{
	  for(iy=0; iy<=DATA->ny; iy++)
	    {
	      if (ieta==0 || ieta==DATA->neta-1)
		{
		  if (ix==0 || ix==DATA->nx)
		    {
		      if (iy==0 || iy==DATA->ny)
			{
			  TtautauPart += 0.125*(arena[ix][iy][ieta].TJb[0][0][0]*(tau+DATA->delta_tau)-arena[ix][iy][ieta].prev_T00*(tau));
			  TetaetaPart += 0.125*(DATA->delta_tau*0.5*(arena[ix][iy][ieta].TJb[0][3][3]+arena[ix][iy][ieta].prev_T33));
			}
		      else
			{
			  TtautauPart += 0.25*(arena[ix][iy][ieta].TJb[0][0][0]*(tau+DATA->delta_tau)-arena[ix][iy][ieta].prev_T00*(tau));
			  TetaetaPart += 0.25*(DATA->delta_tau*0.5*(arena[ix][iy][ieta].TJb[0][3][3]+arena[ix][iy][ieta].prev_T33));
			}
		    }
		  else
		    {
		      if (iy==0 || iy==DATA->ny)
			{
			  TtautauPart += 0.25*(arena[ix][iy][ieta].TJb[0][0][0]*(tau+DATA->delta_tau)-arena[ix][iy][ieta].prev_T00*(tau));
			  TetaetaPart += 0.25*(DATA->delta_tau*0.5*(arena[ix][iy][ieta].TJb[0][3][3]+arena[ix][iy][ieta].prev_T33));
			}
		      else
			{
			  TtautauPart += 0.5*(arena[ix][iy][ieta].TJb[0][0][0]*(tau+DATA->delta_tau)-arena[ix][iy][ieta].prev_T00*(tau));
			  TetaetaPart += 0.5*(DATA->delta_tau*0.5*(arena[ix][iy][ieta].TJb[0][3][3]+arena[ix][iy][ieta].prev_T33));
			}
		    }
		}
	      else
		{  if (ix==0 || ix==DATA->nx)
		    {
		      if (iy==0 || iy==DATA->ny)
			{
			  TtautauPart += 0.25*(arena[ix][iy][ieta].TJb[0][0][0]*(tau+DATA->delta_tau)-arena[ix][iy][ieta].prev_T00*(tau));
			  TetaetaPart += 0.25*(DATA->delta_tau*0.5*(arena[ix][iy][ieta].TJb[0][3][3]+arena[ix][iy][ieta].prev_T33));
			}
		      else
			{
			  TtautauPart += 0.5*(arena[ix][iy][ieta].TJb[0][0][0]*(tau+DATA->delta_tau)-arena[ix][iy][ieta].prev_T00*(tau));
			  TetaetaPart += 0.5*(DATA->delta_tau*0.5*(arena[ix][iy][ieta].TJb[0][3][3]+arena[ix][iy][ieta].prev_T33));
			}
		    }
		  else
		    {
		      if (iy==0 || iy==DATA->ny)
			{
			  TtautauPart += 0.5*(arena[ix][iy][ieta].TJb[0][0][0]*(tau+DATA->delta_tau)-arena[ix][iy][ieta].prev_T00*(tau));
			  TetaetaPart += 0.5*(DATA->delta_tau*0.5*(arena[ix][iy][ieta].TJb[0][3][3]+arena[ix][iy][ieta].prev_T33));
			}
		      else
			{
			  TtautauPart += arena[ix][iy][ieta].TJb[0][0][0]*(tau+DATA->delta_tau)-arena[ix][iy][ieta].prev_T00*(tau);
			  TetaetaPart += DATA->delta_tau*0.5*(arena[ix][iy][ieta].TJb[0][3][3]+arena[ix][iy][ieta].prev_T33);
			}
		    }
		}
	    }	  
	}
    }
  value = (TtautauPart+TetaetaPart)/TtautauPart*100.;
  cout << "energy changes by " << value << "% per time step" << endl;
  cout << "Ttautau=" << TtautauPart << endl;
  cout << "Tetaeta=" << TetaetaPart << endl;
  cout << "sum=" << TtautauPart+TetaetaPart << endl;
  fprintf(ecc_file, "%f %f \n", tau, value);
  fclose(ecc_file);
  return;
}/* ComputeEnergyConservation */


void Grid_info::PrintdEdEta(InitData *DATA, Grid ***arena)
{
 int ix, iy, ieta;
//  int i;
 double eta, f, g, h, l, k;
//  double x, y;

 ix = DATA->nx/2; 
 iy = DATA->ny/2; 
 
 for(ieta=0; ieta<DATA->neta; ieta++)
  {
//   x = ix*(DATA->delta_x) - (DATA->x_size/2.0);  
   eta = ieta*(DATA->delta_eta) - (DATA->eta_size/2.0);  
   f = arena[ix][iy][ieta].u[0][1];
   g = arena[ix][iy][ieta].u[0][2];
   h = arena[ix][iy][ieta].u[0][3];
   l = arena[ix][iy][ieta].u[0][0];
   k = arena[ix][iy][ieta].epsilon;
   /* everything is in fm. epsilon is in fm^-4 */
   printf("%e  %e  %e  %e  %e\n", eta, f/l, g/l, h/l, k*hbarc);
  }
}

void Grid_info::Output_hydro_information_header(InitData *DATA, EOS *eos)
{
	string fname = "hydro_info_header_h";

	//Open output file
	ofstream outfile;
	outfile.open(fname.c_str());

      int grid_nx = ceil(((double)(DATA->nx+1))/DATA->output_evolution_every_N_x);
      int grid_ny = ceil(((double)(DATA->ny+1))/DATA->output_evolution_every_N_y);
      int grid_neta = ceil(((double)(DATA->neta*MPI::COMM_WORLD.Get_size()))/DATA->output_evolution_every_N_eta);

	outfile << "const int MUSIC_real_nx=" << grid_nx << ";" << endl;
	outfile << "const int MUSIC_real_ny=" << grid_ny << ";" << endl;
	//DATA->neta is _not_ the actual number of cells in eta, it is the number of cells in eta _per_ processor
	outfile << "const int MUSIC_real_neta=" << grid_neta << ";" << endl;

	//double x_size; /* in fermi -x_size/2 < x < x_size/2 */
	//double y_size; /* in fermi, ditto */
	//double eta_size; /* ditto */

	outfile << "const double MUSIC_tau0=" << DATA->tau0 << ";" << endl;

	outfile << "const double MUSIC_dx=" << DATA->delta_x*DATA->output_evolution_every_N_x << ";" << endl;
	outfile << "const double MUSIC_dy=" << DATA->delta_y*DATA->output_evolution_every_N_y << ";" << endl;
	outfile << "const double MUSIC_deta=" << DATA->delta_eta*DATA->output_evolution_every_N_eta << ";" << endl;
	outfile << "const double MUSIC_effective_dtau=" << DATA->output_evolution_every_N_timesteps*DATA->delta_tau << ";" << endl;

	outfile << "const bool MUSIC_with_shear_viscosity=" << ((DATA->viscosity_flag)&&(DATA->turn_on_shear)) << ";\n";
	outfile << "const bool MUSIC_with_bulk_viscosity=" << ((DATA->viscosity_flag)&&(DATA->turn_on_bulk)) << ";\n";

      if(DATA->turn_on_rhob == 0)
      {
	    outfile << "const double MUSIC_kinetic_FO_temperature_in_GeV=";
	    if (DATA->useEpsFO) {
	    	outfile << eos->get_temperature(DATA->epsilonFreeze/hbarc,0.0)*hbarc;
	    }
	    else {
	    	outfile << DATA->TFO;
	    }
	    outfile << ";" << endl;
	    outfile << "const int MUSIC_use_temperature_FO=1;" << endl;
      }
      else
      {
	    outfile << "const double MUSIC_kinetic_FO_energy_density_in_GeV_over_fm3=" << DATA->eps_freeze_min << ";" << endl;
	    outfile << "const int MUSIC_use_temperature_FO=0;" << endl;
      }

	outfile << "const bool MUSIC_outputBinaryEvolution=" << DATA->outputBinaryEvolution << ";" << endl;
	outfile << "const bool MUSIC_with_baryon_diffusion=" << ((DATA->viscosity_flag)&&(DATA->turn_on_diff)) << ";" << endl;

	outfile.close();
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

double Grid_info::get_deltaf_qmu_coeff(double T, double muB)
{
    if(muB < 0)
    {
       muB = -muB;
    }
    int idx_T = (int)((T - delta_qmu_coeff_table_T0)/delta_qmu_coeff_table_dT);
    int idx_mu = (int)((muB - delta_qmu_coeff_table_mu0)/delta_qmu_coeff_table_dmu);

    if(idx_T > deltaf_qmu_coeff_table_length_T - 2 || idx_T < 0)
        return(1.0);
    if(idx_mu > deltaf_qmu_coeff_table_length_mu - 2)
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

double Grid_info::get_deltaf_coeff_14moments(double T, double muB, double type)
{
    int idx_T = (int)((T - delta_coeff_table_14mom_T0)/delta_coeff_table_14mom_dT);
    int idx_mu = (int)((muB - delta_coeff_table_14mom_mu0)/delta_coeff_table_14mom_dmu);
    double x_fraction = (T - delta_coeff_table_14mom_T0)/delta_coeff_table_14mom_dT - idx_T;
    double y_fraction = (muB - delta_coeff_table_14mom_mu0)/delta_coeff_table_14mom_dmu - idx_mu;

    double **deltaf_table = NULL;
    if (type == 0)
       deltaf_table = deltaf_coeff_tb_14mom_DPi;
    else if (type == 1)
       deltaf_table = deltaf_coeff_tb_14mom_BPi;
    else if (type == 2)
       deltaf_table = deltaf_coeff_tb_14mom_BPitilde;
    else if (type == 3)
       deltaf_table = deltaf_coeff_tb_14mom_DV;
    else if (type == 4)
       deltaf_table = deltaf_coeff_tb_14mom_BV;
    else if (type == 5)
       deltaf_table = deltaf_coeff_tb_14mom_Bpi_shear;
    else
    {
       cout << "Freeze::get_deltaf_coeff_14moments: unknown type: " << type << endl;
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
