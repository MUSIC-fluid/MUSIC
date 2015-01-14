#include "freeze.h"

// read in thermal spectra from file to then perform resonance decays with them
// Must set verbose to 1 if you want particleMax to be set by this routine.
void Freeze::ReadSpectra_pseudo(InitData* DATA, int full, int verbose)
{
  // read in thermal spectra from file:
  int number, iptmax, iphimax;
  double deltaeta, etamax, ptmin, ptmax;
  int ietamax;
  int pseudo_steps;
  int ip;
//   int bytes_read;
  fprintf(stderr,"reading spectra\n");
  // open particle information file:
  FILE *p_file;
  const char* p_name = "particleInformation.dat";
  const char* pf_name = "FparticleInformation.dat";
//   if(full) strcpy(p_name, "FparticleInformation.dat");
  if(full) p_file = fopen(pf_name, "r");
  else p_file = fopen(p_name, "r");
//   checkForReadError(p_file,p_name);
  int count;
  count = 0;
  if(verbose) cout << "NumberOfParticlesToInclude " << DATA->NumberOfParticlesToInclude << endl;
  // read particle information:
  while( fscanf(p_file,"%d %lf %d %lf %lf %d %d ",&number, &etamax, &ietamax, &ptmin, &ptmax, &iptmax, &iphimax) == 7)
    {
      count ++;
      if (count>DATA->NumberOfParticlesToInclude) break;
// 	  fprintf(stderr,"%d %e %d %e %e %d %d \n", number, etamax, pseudo_steps, ptmin, ptmax, iptmax, iphimax);
      pseudo_steps = ietamax-1;
      ip = partid[MHALF+number];
      particleList[ip].ny = ietamax;
      particleList[ip].npt = iptmax;
      particleList[ip].nphi = iphimax;
      particleList[ip].phimin = 0;
      particleList[ip].phimax = 2*PI;
      particleList[ip].slope = 1;
      particleList[ip].ymax = etamax;
      deltaeta = 0.;
      if(pseudo_steps>0) deltaeta = 2*etamax/pseudo_steps;
      particleList[ip].deltaY = deltaeta;
// 	  cout << "ptmin = " << ptmin << endl;
      for (int ipt=0; ipt<iptmax; ipt++ )
	{
// 	      particleList[ip].pt[ipt] =  ptmin + (ptmax - ptmin)*(exp(ipt)-1)/(exp(iptmax)-1); // log distributed values
	  particleList[ip].pt[ipt] =  ptmin + (ptmax - ptmin)*pow((double)ipt,(double)2)/pow((double)iptmax-1,(double)2); // power law
// 	  particleList[ip].pt[ipt] =  0.01 + ipt*ptmax/(iptmax-1); // compare to UVH2+1
// 	  particleList[ip].pt[ipt] = gala15x[ipt]/12; // gauss laguerre abissas
	}
      for (int ieta=0; ieta<=pseudo_steps; ieta++ )
	{
	  particleList[ip].y[ieta] =  ieta*deltaeta-etamax; // store pseudorapidity here
	  //	  if (ip==1) cout << "read particleList[ip].y[" << ieta << "] = " <<  particleList[ip].y[ieta] << endl;
	}
	
	
      // phiArray for use in Edndp3 interpolation function during resonance decay calculation
      phiArray = util->vector_malloc(iphimax);
      for(int iphi=0; iphi<iphimax; iphi++) phiArray[iphi] = iphi*2*PI/iphimax;
    }
  double pMax = ip;
  if(verbose)
  {
    cout << "particleMax = " << particleMax << endl;
    particleMax = ip;
  }

  fclose(p_file);

  FILE *s_file;
  const char* s_name = "yptphiSpectra.dat";
  const char* sf_name = "FyptphiSpectra.dat";
  if(full) s_file = fopen(sf_name, "r");
  else s_file = fopen(s_name, "r");
//   checkForReadError(s_file,s_name);
  
  if(verbose)
  {
    cout << "ietamax=" << pseudo_steps+1 << endl;
    cout << "cells=" << (pseudo_steps+1)*(iptmax)*iphimax << endl;
    cout << "particleMax = " << particleMax << endl;
  }
  for ( ip=1; ip<=pMax; ip++ )
    {
      //cout << ip << endl;
      if(verbose) fprintf(stderr,"reading particle %d: %d %s\n", ip, particleList[ip].number, particleList[ip].name);
      for (int ieta=0; ieta<=pseudo_steps; ieta++)
	{
	  for (int ipt=0; ipt<iptmax; ipt++)
	    {
	      for (int iphi=0; iphi<iphimax; iphi++)
		{
// 		  bytes_read=fscanf(s_file, "%lf", &particleList[ip].dNdydptdphi[ieta][ipt][iphi]);
		  fscanf(s_file, "%lf", &particleList[ip].dNdydptdphi[ieta][ipt][iphi]);
		  if(particleList[ip].dNdydptdphi[ieta][ipt][iphi]<0.)
		    particleList[ip].dNdydptdphi[ieta][ipt][iphi]=0;
			  //		  cout << particleList[ip].y[ieta] << " " << particleList[ip].pt[ipt] << " " << phiArray[iphi] << " " << particleList[ip].dNdydptdphi[ieta][ipt][iphi] << endl; 
		  //	  printf("%f %f %f \n",particleList[ip].y[ieta],particleList[ip].pt[ipt],phiArray[iphi]);
		}
	    }
	}
    }
  //particleMax=2;
  fclose(s_file);
}

// Improved accuracy when freeze out method 3 is used.  
// Contribution from each surface segment is integrated 
// over its extent in eta using Riemann sum,
// instead of counting entire contribution from center of segment.
double Freeze::summation3(double px, double py, double y, double m, int deg, int baryon, double muAtFreezeOut, InitData *DATA)
{
  
  double sum = 0.;
  double dSigma[4] = {0};
  int i;
  double ptau, peta;
  double f, T, mu, tau, eta, E, sign, delta_f;
  double pdSigma, Wfactor;
  double mt = sqrt(m*m+px*px+py*py); // all in GeV
  double alpha = 0.; // make a parameter
  int subsections;
  double maxDETA = DATA->max_delta_eta2;
  // if max_delta_eta was set smaller than delta_eta, the freezeout surface
  // segment was already subdivided in FindFreezeOutSurface3.  The extent
  // in eta (for cells not oriented in the eta direction) is actually:
  double DETA=DATA->delta_eta/(floor(DATA->delta_eta/DATA->max_delta_eta) + 1);
    // Bose or Fermi statistics.
  if (baryon==0)
    sign = -1.;
  else
    sign = 1.;

  //fprintf(stderr,"sign=%f\n",sign);
  //fprintf(stderr,"baryon=%d\n",baryon);
  for (i=0; i<NCells; i++)
    {
      for (int j=0; j<4; j++) dSigma[j] = surface[i].s[j];
      if(dSigma[3]==0) 
      {
	subsections = floor(DETA/maxDETA) + 1;
	for (int j=0; j<3; j++) dSigma[j]/=subsections;
//       if (subsections > 1) cout << "Splitting surface element into " << subsections << " segments in eta\n";
      }
      else subsections = 1;
      
      tau = surface[i].x[0];

      T = surface[i].T_f*hbarc; // GeV
      mu = baryon*surface[i].mu_B*hbarc; //GeV
      if(DATA->whichEOS>=3) // for PCE use the previously computed mu at the freeze-out energy density
	mu=muAtFreezeOut; //GeV
      
      for (int k=0; k < subsections; k++)
      {
	  eta = surface[i].x[3] - DETA/2 + DETA/2/subsections + (k*DETA)/(subsections);
	  
	  ptau = mt*cosh(y-eta); // GeV    this is p^tau
	  peta = mt/tau*sinh(y-eta); // GeV/fm     this is p^eta

	  // compute p^mu*dSigma_mu
	  pdSigma = tau*(ptau*dSigma[0]+px*dSigma[1]+py*dSigma[2]+peta*dSigma[3]); //fm^3*GeV
	  E = (ptau*surface[i].u[0]-px*surface[i].u[1]-py*surface[i].u[2]-tau*tau*peta*surface[i].u[3]/tau);
	  // this is the equilibrium f, f_0:
	  f = 1./(exp(1./T*(E-mu))+sign);
	  // now comes the delta_f: check if still correct at finite mu_b 
	  // we assume here the same C=eta/s for all particle species because it is the simplest way to do it.
	  // also we assume Xi(p)=p^2, the quadratic Ansatz

	  if(f<0)
	  {
	    cerr << "Mistake in thermal spectrum calculation, f<0\n";
	    f=0;
	  }

	  if(E<=0)
	    f=0;
	  
	  if(T<0)
	    f=0;
	
	  if (DATA->include_deltaf>=1 && DATA->viscosity_flag==1)
	    {
	      Wfactor=(ptau*surface[i].W[0][0]*ptau
		      -2.*ptau*surface[i].W[0][1]*px
		      -2.*ptau*surface[i].W[0][2]*py
		      -2.*tau*tau*ptau*surface[i].W[0][3]/tau*peta
		      +px*surface[i].W[1][1]*px
		      +2.*px*surface[i].W[1][2]*py
		      +2.*tau*tau*px*surface[i].W[1][3]/tau*peta
		      +py*surface[i].W[2][2]*py
		      +2.*tau*tau*py*surface[i].W[2][3]/tau*peta
		      +tau*tau*tau*tau*peta*surface[i].W[3][3]/tau/tau*peta)
		*pow(hbarc,4.); // W is like energy density
	      
	      
	      delta_f = f*(1.-sign*f)/(2.*surface[i].eps_plus_p_over_T_FO*pow(hbarc,3.)*pow(T,3.))*Wfactor;
	    
	      if (DATA->include_deltaf==2) // if delta f is supposed to be proportional to p^(2-alpha):
		{
		  delta_f = delta_f * pow((T/E),1.*alpha)*120./(tgamma(6.-alpha)); 
		}
	      
	    }
	  else
	    {
	      delta_f=0.;
	    }

	    sum += 1/pow(2.*PI,3.) * (f+delta_f) * pdSigma;
	    //	cout << "sum=" << sum << ", pdSigma=" << pdSigma << endl;
	    if (sum>10000)
	      cout << "WARNING: sum>10000 in summation. sum=" << sum << ", f=" << f << ", deltaf=" << delta_f << ", pdSigma=" << pdSigma 
		  << ", T=" << T << ", E=" << E << ", mu=" << mu << endl;
      }// loop over subsections
    }

  sum *= deg/pow(hbarc,3.); // in GeV^(-2)

  return sum;
}


//Modified spectra calculation by ML 05/2013
//Calculates on fixed grid in pseudorapidity, pt, and phi
void Freeze::ComputeParticleSpectrum_pseudo(InitData *DATA, int number, int anti, int size, int rank)
{
//   char *specString;
//   specString = util->char_malloc(30);
  int j;
  double pt, phi, px, py;
  double eta;
 

  j = partid[MHALF+number];
  // set some parameters
  double etamax = DATA->max_pseudorapidity;
  int ietamax = DATA->pseudo_steps + 1;// pseudo_steps is number of steps.  Including edges, number of points is steps + 1
  double deltaeta = 0;
  if(ietamax>1) deltaeta = 2*etamax/DATA->pseudo_steps;
  double ptmax = DATA->max_pt;
  double ptmin = DATA->min_pt;
  int iptmax = DATA->pt_steps+1; // Number of points is iptmax + 1 (ipt goes from 0 to iptmax)
//   double deltapt = ptmax/iptmax;
  int iphimax = DATA->phi_steps; // Number of steps equal to number of points (phi=2pi is understood as equal to phi=0)
  double deltaphi = 2*PI/iphimax;
  
  //  Parallellize in pseudorapidity.
  //  Allow for unequal load division if number of points does
  //  not divide equally among processors.  Not efficient, but
  //  might be convenient to have the option (e.g. for mode 1).
  int rem = ietamax%size;
  ietamax/=size;
  if(rank<rem) ietamax++;
//   cout << "rank = " << rank << ", ietamax = " << ietamax << endl;

//   if (size>1)
//     {
//       if (ietamax%size!=1)
// 	{
// 	  fprintf(stderr,"number of steps in pseudorapidity (pseudo_steps) is not a multiple of the number of processors. Exiting.\n");
// 	  MPI::Finalize();
// 	  exit(1);
// 	}
//       ietamax/=size;
// //       cout << "r" << rank << " ietamax=" << ietamax << endl;
//     }
	  
	  
// Reuse rapidity variables (Need to reuse variable y so resonance decay routine can be used as is.
// 	  Might as well misuse ymax and deltaY too)
  particleList[j].ymax = etamax; 
  particleList[j].deltaY = deltaeta;

  fprintf(stderr,"Doing %d: %s (%d)\n", j, particleList[j].name, particleList[j].number);
 
//   particleList[j].ny = ietamax*size;
  particleList[j].ny = DATA->pseudo_steps + 1;
  particleList[j].npt = iptmax;
  particleList[j].nphi = iphimax;
  
  // set particle properties
  double m = particleList[j].mass;
  int d = particleList[j].degeneracy;
  int b = particleList[j].baryon;
//   int s = particleList[j].strange;
//   double c = particleList[j].charge;
  double mu = particleList[j].muAtFreezeOut;

  char *buf = new char[10];
  sprintf (buf, "%d", number);  
 
  // open files to write
  FILE *d_file;
  const char* d_name = "particleInformation.dat";
  d_file = fopen(d_name, "a");

  char *specString=new char[30];
  
  FILE *s_file;
  sprintf (buf, "%d", rank);
  
  strcpy(specString, "yptphiSpectra");
  strcat(specString, buf);
  strcat(specString, ".dat");
  char* s_name = specString;
  s_file = fopen(s_name, "w");
  delete[] specString;
  delete[] buf;

  // --------------------------------------------------------------------------
  
  // write information in the particle information file
  if (rank==0) fprintf(d_file,"%d %e %d %e %e %d %d \n", number,  etamax, DATA->pseudo_steps+1, ptmin, ptmax, iptmax, iphimax);
//   if (rank==0) fprintf(d_file,"%d %e %e %d %e %d %d \n", number, deltaeta, etamax,  DATA->pseudo_steps, ptmax, iptmax, iphimax);

  
  
  // store E dN/d^3p as function of phi, pt and eta (pseudorapidity) in sumPtPhi:
  
  for (int ieta=0; ieta<ietamax; ieta++)
    {
//       eta = -etamax + deltaeta/2. + ieta*deltaeta + rank*(etamax/size*2.);
      double offset;
      if(rank<=rem) offset = deltaeta*rank*floor((DATA->pseudo_steps+1)/size + 1);
      else offset = deltaeta*(rank*floor((DATA->pseudo_steps+1)/size) + rem);
      eta = -etamax + ieta*deltaeta + offset;
      //     if (j==1) cout << " do particleList[ip].y[" << iy << "] = " <<  particleList[j].y[iy] << " y = " << y << endl;

      for (int ipt=0; ipt<iptmax; ipt++)
	{
// 	  pt = deltapt/2. + ipt*deltapt;
// 	  pt = ipt*deltapt;
// 	  pt = ptmin + (ptmax - ptmin)*(exp(ipt)-1)/(exp(iptmax)-1);
// 	      pt =  ptmin + (ptmax - ptmin)*(exp(ipt)-1)/(exp(iptmax)-1); // log distributed values
	      pt =  ptmin + (ptmax - ptmin)*pow((double)ipt,(double)2)/pow((double)iptmax-1,(double)2); // power law
// 	      pt =  0.01 + ipt*ptmax/(iptmax-1); // compare to UVH2+1
// 	      pt = gala15x[ipt]/12.; // gauss laguerre absissas
	  particleList[j].pt[ipt] = pt;
	  
	  
	  //rapidity as a function of pseudorapidity:
	  double y = Rap(eta,pt,m);
	  
	  // Use this variable to store pseudorapidity instead of rapidity
	  // May cause confusion in the future, but easier to to share code for both options:
	  // calculating on a fixed grid in rapidity or pseudorapidity
	  particleList[j].y[ieta] = eta;
	  
	  
	  for (int iphi=0; iphi<iphimax; iphi++)
	    {
	      phi = deltaphi*iphi;
	      px = pt*cos(phi);
	      py = pt*sin(phi);
	      double sum;
// 	      if(DATA->freezeOutMethod==3) cout << "Calling summation3\n";
	      if(DATA->freezeOutMethod==3) sum = summation3(px, py, y, m, d, b, mu, DATA);
	      else sum = summation(px, py, y, m, d, b, mu, DATA);
	      particleList[j].dNdydptdphi[ieta][ipt][iphi] = sum;
	      fprintf(s_file,"%e ", sum);
	    }
	  fprintf(s_file,"\n");
	}
      
    }
  fclose(s_file);
  fclose(d_file);
  
}

// adapted from ML
void Freeze::ComputeParticleSpectrum_pseudo_improved(InitData *DATA, int number, int anti, int size, int rank)
{
  double y_minus_eta_cut = 10.0;
  //char *specString;
  //specString = util->char_malloc(30);
  int j = partid[MHALF+number];
  // set some parameters
  double etamax = DATA->max_pseudorapidity;
  int ietamax = DATA->pseudo_steps + 1;// pseudo_steps is number of steps.  Including edges, number of points is steps + 1
  double deltaeta = 0;
  if(ietamax>1)
     deltaeta = 2*etamax/DATA->pseudo_steps;
  double ptmax = DATA->max_pt;
  double ptmin = DATA->min_pt;
  int iptmax = DATA->pt_steps+1; // Number of points is iptmax + 1 (ipt goes from 0 to iptmax)
  //double deltapt = ptmax/iptmax;
  int iphimax = DATA->phi_steps; // Number of steps equal to number of points (phi=2pi is understood as equal to phi=0)
  double deltaphi = 2*PI/iphimax;
  
  //  Parallellize in pseudorapidity.
  //  Allow for unequal load division if number of points does
  //  not divide equally among processors.  Not efficient, but
  //  might be convenient to have the option (e.g. for mode 1).
  int rem = ietamax%size;
  ietamax/=size;
  if(rank<rem) ietamax++;
//   cout << "rank = " << rank << ", ietamax = " << ietamax << endl;

//   if (size>1)
//     {
//       if (ietamax%size!=1)
// 	{
// 	  fprintf(stderr,"number of steps in pseudorapidity (pseudo_steps) is not a multiple of the number of processors. Exiting.\n");
// 	  MPI::Finalize();
// 	  exit(1);
// 	}
//       ietamax/=size;
// //       cout << "r" << rank << " ietamax=" << ietamax << endl;
//     }
	  
	  
  // Reuse rapidity variables (Need to reuse variable y so resonance decay routine can be used as is.
  // Might as well misuse ymax and deltaY too)
  particleList[j].ymax = etamax; 
  particleList[j].deltaY = deltaeta;

  fprintf(stderr,"Doing %d: %s (%d)\n", j, particleList[j].name, particleList[j].number);
 
  //particleList[j].ny = ietamax*size;
  particleList[j].ny = DATA->pseudo_steps + 1;
  particleList[j].npt = iptmax;
  particleList[j].nphi = iphimax;
  
  // set particle properties
  double m = particleList[j].mass;
  int deg = particleList[j].degeneracy;
  int baryon = particleList[j].baryon;
  //int s = particleList[j].strange;
  //double c = particleList[j].charge;
  double mu_PCE = particleList[j].muAtFreezeOut;
  int sign;
  if (baryon==0)
    sign = -1.;
  else
    sign = 1.;

  // open files to write
  FILE *d_file;
  const char* d_name = "particleInformation.dat";
  d_file = fopen(d_name, "a");

  char *buf = new char[10];
  sprintf (buf, "%d", number);  
  sprintf (buf, "%d", rank);
  
  char *specString=new char[30];
  strcpy(specString, "yptphiSpectra");
  strcat(specString, buf);
  strcat(specString, ".dat");
  char* s_name = specString;
  FILE *s_file;
  s_file = fopen(s_name, "w");

  delete[] specString;
  delete[] buf;

  // --------------------------------------------------------------------------
  // write information in the particle information file
  if (rank==0)
     fprintf(d_file,"%d %e %d %e %e %d %d \n", number,  etamax, DATA->pseudo_steps+1, ptmin, ptmax, iptmax, iphimax);
  //if (rank==0) fprintf(d_file,"%d %e %e %d %e %d %d \n", number, deltaeta, etamax,  DATA->pseudo_steps, ptmax, iptmax, iphimax);

  // caching 
  double* cos_phi = new double [iphimax];
  double* sin_phi = new double [iphimax];
  for(int iphi = 0; iphi < iphimax; iphi++)
  {
     double phi_local = deltaphi*iphi;
     cos_phi[iphi] = cos(phi_local);
     sin_phi[iphi] = sin(phi_local);
  }
  double* pt_array = new double [iptmax];
  for(int ipt = 0; ipt < iptmax; ipt++)
  {
     //pt = deltapt/2. + ipt*deltapt;
     //pt = ipt*deltapt;
     //pt = ptmin + (ptmax - ptmin)*(exp(ipt)-1)/(exp(iptmax)-1);
     //pt =  ptmin + (ptmax - ptmin)*(exp(ipt)-1)/(exp(iptmax)-1); // log distributed values
     //pt =  0.01 + ipt*ptmax/(iptmax-1); // compare to UVH2+1
     //pt = gala15x[ipt]/12.; // gauss laguerre absissas
     double pt =  ptmin + (ptmax - ptmin)*pow((double)ipt,(double)2)/pow((double)iptmax-1,(double)2); // power law
     pt_array[ipt] = pt;
     particleList[j].pt[ipt] = pt;
  }
      
  double alpha = 0.0;
  // main loop begins ...
  // store E dN/d^3p as function of phi, pt and eta (pseudorapidity) in sumPtPhi:
  for (int ieta=0; ieta<ietamax; ieta++)
  {
     //eta = -etamax + deltaeta/2. + ieta*deltaeta + rank*(etamax/size*2.);
     double offset;
     if(rank <= rem)
        offset = deltaeta*rank*floor((DATA->pseudo_steps+1)/size + 1);
     else 
        offset = deltaeta*(rank*floor((DATA->pseudo_steps+1)/size) + rem);
     double eta = -etamax + ieta*deltaeta + offset;

     // Use this variable to store pseudorapidity instead of rapidity
     // May cause confusion in the future, but easier to to share code for both options:
     // calculating on a fixed grid in rapidity or pseudorapidity
     particleList[j].y[ieta] = eta; // store particle pseudo-rapidity

     double* rapidity = new double [iptmax];
     double* cosh_y = new double [iptmax];
     double* sinh_y = new double [iptmax];
     for (int ipt=0; ipt<iptmax; ipt++)
     {
        double pt = pt_array[ipt];
        //rapidity as a function of pseudorapidity:
        double y_local = Rap(eta, pt, m);
        rapidity[ipt] = y_local;
        cosh_y[ipt] = cosh(y_local);
        sinh_y[ipt] = sinh(y_local);
     }

     double** temp_sum = new double* [iptmax];
     for(int ii = 0; ii < iptmax; ii++)
     {
        temp_sum[ii] = new double [iphimax];
        for(int jj = 0; jj < iphimax; jj++)
            temp_sum[ii][jj] = 0.0;
     }
       
     for (int icell = 0; icell < NCells; icell++)
     {
        double tau = surface[icell].x[0];
        double eta_s = surface[icell].x[3];
        double cosh_eta_s = surface[icell].cosh_eta_s;
        double sinh_eta_s = surface[icell].sinh_eta_s;

        double T = surface[icell].T_f*hbarc; // GeV
        double mu = baryon*surface[icell].mu_B*hbarc; //GeV
        if(DATA->whichEOS>=3 && DATA->whichEOS < 10) // for PCE use the previously computed mu at the freeze-out energy density
	     mu += mu_PCE; //GeV

        double sigma_mu[4];
        double u_flow[4];
        for(int ii = 0; ii < 4; ii++)
        {
           sigma_mu[ii] = surface[icell].s[ii];
           u_flow[ii] = surface[icell].u[ii];
        }

        double W00 = surface[icell].W[0][0];
        double W01 = surface[icell].W[0][1];
        double W02 = surface[icell].W[0][2];
        double W03 = surface[icell].W[0][3];
        double W11 = surface[icell].W[1][1];
        double W12 = surface[icell].W[1][2];
        double W13 = surface[icell].W[1][3];
        double W22 = surface[icell].W[2][2];
        double W23 = surface[icell].W[2][3];
        double W33 = surface[icell].W[3][3];

        double eps_plus_P_over_T = surface[icell].eps_plus_p_over_T_FO;
        double prefactor_shear = 1./(2.*eps_plus_P_over_T*T*T*T)*hbarc; // fm^4/GeV^2

        for (int ipt=0; ipt<iptmax; ipt++)
        {
           double pt = pt_array[ipt];
           double cosh_y_local = cosh_y[ipt];
           double sinh_y_local = sinh_y[ipt];
           double mt = sqrt(m*m + pt*pt); // all in GeV
           double y = rapidity[ipt];
           if( fabs(y - eta_s) < y_minus_eta_cut)
           {
              // ptau = mt*cosh(y - eta_s) // GeV
              double ptau = mt*(cosh_y_local*cosh_eta_s - sinh_y_local*sinh_eta_s); 
              // peta = mt/tau*sinh(y - eta_s); // GeV/fm
              double peta = mt/tau*(sinh_y_local*cosh_eta_s - cosh_y_local*sinh_eta_s); 

              for (int iphi=0; iphi<iphimax; iphi++)
              {
                 double px = pt*cos_phi[iphi];
                 double py = pt*sin_phi[iphi];
              
                 double sum;

                 // compute p^mu*dSigma_mu
                 double pdSigma = tau*(ptau*sigma_mu[0]+px*sigma_mu[1]+py*sigma_mu[2]+peta*sigma_mu[3]); //fm^3*GeV
                 //double E = (ptau*u_flow[0] - px*u_flow[1] - py*u_flow[2]- tau*tau*peta*u_flow[3]/tau);
                 double E = (ptau*u_flow[0] - px*u_flow[1] - py*u_flow[2]- tau*peta*u_flow[3]);
                 // this is the equilibrium f, f_0:
                 double f = 1./(exp(1./T*(E - mu)) + sign);

                 // now comes the delta_f: check if still correct at finite mu_b 
                 // we assume here the same C=eta/s for all particle species because it is the simplest way to do it.
                 // also we assume Xi(p)=p^2, the quadratic Ansatz
                 double Wfactor = 0.0;
                 double delta_f_shear = 0.0;
                 if (DATA->include_deltaf>=1 && DATA->viscosity_flag==1)
	           {
                    Wfactor = (  ptau*W00*ptau - 2.*ptau*W01*px - 2.*ptau*W02*py - 2.*tau*tau*ptau*W03/tau*peta
	           	          + px*W11*px + 2.*px*W12*py + 2.*tau*tau*px*W13/tau*peta
	           	          + py*W22*py + 2.*tau*tau*py*W23/tau*peta
	           	          +tau*tau*tau*tau*peta*W33/tau/tau*peta); // W is like energy density
	             
	             delta_f_shear = f*(1.-sign*f)*prefactor_shear*Wfactor;
	           
	             if (DATA->include_deltaf==2) // if delta f is supposed to be proportional to p^(2-alpha):
	             {
	                delta_f_shear = delta_f_shear * pow((T/E),1.*alpha)*120./(tgamma(6.-alpha)); 
	             }
	             
	           }
                 else
	           {
	              delta_f_shear = 0.;
	           }
	           sum = (f + delta_f_shear) * pdSigma;
	           if (sum>10000)
	              cout << "WARNING: sum>10000 in summation. sum=" << sum 
                         << ", f=" << f << ", deltaf=" << delta_f_shear
                         << ", pdSigma=" << pdSigma << ", T=" << T 
                         << ", E=" << E << ", mu=" << mu << endl;
                 temp_sum[ipt][iphi] += sum;
              }
           }
        }
     }

     double prefactor = deg/(pow(2.*PI,3.)*pow(hbarc,3.));
     // store the final results
     for(int ipt = 0; ipt < iptmax; ipt++)
     {
        for(int iphi = 0; iphi < iphimax; iphi++)
        {
            double sum = temp_sum[ipt][iphi]*prefactor;   // in GeV^(-2)
            particleList[j].dNdydptdphi[ieta][ipt][iphi] = sum;
            fprintf(s_file,"%e ", sum);
        }
        fprintf(s_file,"\n");
     }

     //clean up
     delete [] rapidity;
     delete [] cosh_y;
     delete [] sinh_y;
     for(int ipt = 0; ipt < iptmax; ipt++)
        delete [] temp_sum[ipt];
     delete [] temp_sum;
  }
  //clean up
  delete [] cos_phi;
  delete [] sin_phi;
  delete [] pt_array;
  fclose(s_file);
  fclose(d_file);
}


// part of new cooper frye calculation from ML 5/2013
void Freeze::OutputFullParticleSpectrum_pseudo(InitData *DATA, int number, int anti, int full)
{

  FILE *d_file;
  
//   string d_name;
//   if(full) d_name = "F";
//   else d_name = "";
//   string s_name;
//   if(full) s_name = "F";
//   else s_name = "";
//       // open files to write
//   d_name += "particleInformation.dat";
  const char* d_name = "FparticleInformation.dat";
  d_file = fopen(d_name, "a");
  
  fprintf(d_file,"%d %e %d %e %e %d %d \n", number,  DATA->max_pseudorapidity, DATA->pseudo_steps+1, DATA->min_pt, DATA->max_pt, DATA->pt_steps+1, DATA->phi_steps);

//   fprintf(d_file,"%d %e %e %e %e %e %d %d %d \n", number, DATA->, ymax, slope, phimin, phimax, iymax, iptmax, iphimax);
   fclose(d_file);

  
  FILE *s_file;
      
//   s_name += "yptphiSpectra.dat";
  const char* s_name = "FyptphiSpectra.dat";
  s_file = fopen(s_name, "a");

  int j = partid[MHALF+number];
  for (int iy=0; iy<=DATA->pseudo_steps; iy++)
    {
      for (int ipt=0; ipt<=DATA->pt_steps; ipt++)
	{
	  for (int iphi=0; iphi<DATA->phi_steps; iphi++)
	    {

	      fprintf(s_file,"%e ", particleList[j].dNdydptdphi[iy][ipt][iphi]);
	    }
	  fprintf(s_file,"\n");
	}
    }     
  fclose(s_file);
}


//Cooper-Frye routine adapted by ML 05/2013
//-- spectra calculated on an equally-spaced grid in phi and pseudorapidity,
//for ease in comparing to experimental data (and improved accuracy in azimuthal integrals)

void Freeze::CooperFrye_pseudo(int particleSpectrumNumber, int mode, InitData *DATA, EOS *eos, int size, int rank)
{
  if(DATA->pseudofreeze) pseudofreeze = 1;
  else pseudofreeze =0;
  int alreadyread = 0;
  ReadParticleData(DATA, eos); // read in data for Cooper-Frye
  int i, b, number;
  if (mode == 3 || mode == 1) // compute thermal spectra
    {
      if (size > (DATA->pseudo_steps + 1)) 
      {
	cout << "Cannot run Cooper-Frye with " << DATA->pseudo_steps << " steps in pseudorapidity on " << size << " processors.  Exiting\n";
	exit(1);
      }

      if (rank == 0) system("rm yptphiSpectra.dat yptphiSpectra?.dat yptphiSpectra??.dat particleInformation.dat 2> /dev/null");
//       MPI_Barrier(MPI_COMM_WORLD);
      ReadFreezeOutSurface(DATA); // read freeze out surface (has to be done after the evolution of course)
      if (particleSpectrumNumber==0) // do all particles
	{
	  fprintf(stderr,"Doing all particles. May take a while ... \n");
	  for ( i=1; i<particleMax; i++ )
	    {
	      number = particleList[i].number;
	      b = particleList[i].baryon;
	      int computespectrum = 1;
	      
	      // Only calculate particles with unique mass (and/or baryon number)
	      for(int part = 1; part < i; part++)
	      {
		if(
		  (particleList[i].mass == particleList[part].mass) &&
		  (DATA->turn_on_rhob == 0 || particleList[i].baryon==particleList[part].baryon)
		)
		{
// 		  cout << "rank " << rank << ", same mass " << particleList[i].number << ", " << particleList[part].number << endl;
		  computespectrum = 0;
		  if(rank==0)
		  {
		    // If there is more than one processor, this processor doesn't have all pseudorapidity values in memory
		    if(size > 1 && part > alreadyread)  
		    {
		      alreadyread = part;
		      ReadSpectra_pseudo(DATA, 0, 0);
		    }
		    fprintf(stderr,"Copying %d: %s (%d) from %s\n", i, particleList[i].name, particleList[i].number, particleList[part].name);
		    int iphimax = DATA->phi_steps;
		    int iptmax = DATA->pt_steps + 1;
		    int ietamax = DATA->pseudo_steps + 1;
		    double ptmax = DATA->max_pt;
		    double ptmin = DATA->min_pt;
		    double etamax = DATA->max_pseudorapidity;
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
		    fprintf(d_file,"%d %e %d %e %e %d %d \n", number,  etamax, ietamax, ptmin, ptmax, iptmax, iphimax);
		    for (int ieta=0; ieta<ietamax; ieta++)
		    for (int ipt=0; ipt<iptmax; ipt++)
		      {
			particleList[i].pt[ipt] = particleList[part].pt[ipt];
			particleList[i].y[ieta] = particleList[part].y[ieta];  
			for (int iphi=0; iphi<iphimax; iphi++)
			  {
			    particleList[i].dNdydptdphi[ieta][ipt][iphi] = particleList[part].dNdydptdphi[ieta][ipt][iphi];
			    fprintf(s_file,"%e ", particleList[i].dNdydptdphi[ieta][ipt][iphi]);
			  }
			fprintf(s_file,"\n");
		      }
		    fclose(s_file);
		    fclose(d_file);
		  }// if rank==0
		  part=particleMax; // break out of particle loop
		}// if particles have same mass
	      }// loop over particles that have already been calculated
	      

	      if(computespectrum) 
	      {
		//ComputeParticleSpectrum_pseudo(DATA, number, b, size, rank);
		ComputeParticleSpectrum_pseudo_improved(DATA, number, b, size, rank);

  
		// Wait until all processors are finished and concatonate results
		// (Don't want any processors to start writing the next particle to 
		// file until the concatonation is done)
		MPI_Barrier(MPI_COMM_WORLD);
		
		if(rank==0) 
		{
		  system("cat yptphiSpectra?.dat >> yptphiSpectra.dat");
		  system("cat yptphiSpectra??.dat >> yptphiSpectra.dat 2> /dev/null");
		}
		
		MPI_Barrier(MPI_COMM_WORLD); // occasionally another processor will open file before concatonation is done if this is removed
	      }
	      
	    }
	}
      else
	{
	  if (particleSpectrumNumber>=particleMax)
	    {
	      fprintf(stderr,"No particle has the number %d. Exiting.\n",particleSpectrumNumber); exit(1);
	    }  
	  number = particleList[particleSpectrumNumber].number;
	  b = particleList[particleSpectrumNumber].baryon;
	  cout << "COMPUTE" << endl;
	  //ComputeParticleSpectrum_pseudo(DATA, number,  b,  size, rank);
	  ComputeParticleSpectrum_pseudo_improved(DATA, number,  b,  size, rank);
	  
	  // send something just to make sure that rank 0 waits for all others to be done:
	  int check;
	  check=rank;
	  
	  if (rank > 0)
	    MPI::COMM_WORLD.Send(&check,1,MPI::INT,0,1);
	    
	  if (rank == 0)
	    {
	      for (int from=1; from < size; from ++)
		MPI::COMM_WORLD.Recv(&check,1,MPI::INT,from,1);
	      
	      system("cat yptphiSpectra?.dat >> yptphiSpectra.dat");
	      system("cat yptphiSpectra??.dat >> yptphiSpectra.dat 2> /dev/null");
	      
	    }
	}
	free(surface);
    }
  if (mode==4 || mode==1)  //  do resonance decays
    {
      ReadSpectra_pseudo(DATA, 0, 1);
      int bound = 211; //number of lightest particle to calculate. 
	  cout << "particleaMax = " << particleMax << endl;
	  fprintf(stderr,"doing all from %i: %s to %i: %s.\n",particleMax,particleList[particleMax].name,
		  partid[MHALF+bound],particleList[partid[MHALF+bound]].name);
//       if(rank==0)  cal_reso_decays_pseudo(particleMax,decayMax,bound,mode, DATA->pseudofreeze);
      if(rank==0)  cal_reso_decays(particleMax,decayMax,bound,mode);

      if (rank == 0) system("rm FyptphiSpectra.dat FparticleInformation.dat 2> /dev/null");
	  for ( i=1; i<particleMax; i++ )
	    {
	      number = particleList[i].number;
	      b = particleList[i].baryon;
	      if(rank==0) OutputFullParticleSpectrum_pseudo(DATA, number, b, 1);
	    }
    }
  else if (mode==13) // take tabulated spectra and compute various observables and integrated quantities
    {
    
      ReadSpectra_pseudo(DATA, 0, 1);
// //       for ( i=1; i<particleMax; i++ )
 	const int stable_charged_hadron_list[] = {211,-211,321,-321,2212,-2212};
 	Output_charged_hadrons_eta_differential_spectra(DATA, 0, stable_charged_hadron_list,sizeof(stable_charged_hadron_list)/sizeof(int));
 	int other_hadron = 211;
//	Output_charged_hadrons_eta_differential_spectra(DATA, 0,  &other_hadron,1);
        Output_midrapidity_hadrons_spectra(DATA, 0,  &other_hadron,1);
// 	int other_hadron = 211;
// 	Output_charged_hadrons_eta_differential_spectra(DATA, 0, &other_hadron,1);
      
      int neta = particleList[1].ny;

//       int chargedhd[6] = {1,3,4,5,17,18};
      int positivech[3] = {1,4,17};

      for ( int k=0; k<3; k++ )
	{
	  i = positivech[k];
	  number = particleList[i].number;
	  if(particleMax>=i)
	  {
	    OutputDifferentialFlowAtMidrapidity(DATA, number, 0);
	    OutputDifferentialFlowNearMidrapidity(DATA, number, 0);
	    OutputIntegratedFlow(DATA, number, 0);
	  }
	}
	
      
      //Set output file name for eta differential yield
      string fname;
      stringstream tmpStr;
      fname="./outputs/dNchdeta.dat";	
      
      ofstream outfile;

      outfile.open(fname.c_str(),ios::trunc);

       outfile << "#eta\tdNch/deta\n";
	for(int ieta=0; ieta<neta; ieta++) 
	  outfile 
// 	  << "eta, dNch/deta = " 
	  << particleList[1].y[ieta] << "\t"  
	  << get_Nch(DATA, particleList[1].pt[0], particleList[1].pt[particleList[i].npt-1], 0, particleList[1].y[ieta], particleList[1].y[ieta]) 
	  << endl;


	outfile.close();
	
	cout << "pion <pt> for 0.25 < pt < 1.8 GeV, |y| < 1 = " 
	<< get_meanpt(DATA, 211, 0.25, 1.8, 1, -1, 1) << endl;
	
    }
  else if (mode==14) // take tabulated post-decay spectra and compute various observables and integrated quantities
    {
      ReadSpectra_pseudo(DATA, 1, 1);
      
	const int stable_charged_hadron_list[] = {211,-211,321,-321,2212,-2212};
	Output_charged_hadrons_eta_differential_spectra(DATA, 1, stable_charged_hadron_list,sizeof(stable_charged_hadron_list)/sizeof(int));

 	int other_hadron = 211;
//	Output_charged_hadrons_eta_differential_spectra(DATA, 0,  &other_hadron,1);
        Output_midrapidity_hadrons_spectra(DATA, 1,  &other_hadron,1);
	//Output_charged_hadrons_eta_differential_spectra(DATA, 1);
	
      int positivech[3] = {1,4,17};
      for ( int k=0; k<3; k++ )
	{
	  i = positivech[k];
	  number = particleList[i].number;
	  if(particleMax>=i)
	  {
	    OutputDifferentialFlowAtMidrapidity(DATA, number, 1);
	    OutputDifferentialFlowNearMidrapidity(DATA, number, 1);
	    OutputIntegratedFlow(DATA, number, 1);
	  }
	}
	
      //Set output file name for eta differential yield
      string fname;
      stringstream tmpStr;
      fname="./outputs/FdNchdeta.dat";	
      
      ofstream outfile;

      outfile.open(fname.c_str(),ios::trunc);
      
       int neta = particleList[1].ny;
       outfile << "#eta\tdNch/deta\n";
	for(int ieta=0; ieta<neta; ieta++) 
	  outfile 
// 	  << "eta, dNch/deta = " 
	  << particleList[1].y[ieta] << "\t" 
	  << get_Nch(DATA, particleList[i].pt[0], particleList[1].pt[particleList[i].npt-1], 0, particleList[1].y[ieta], particleList[1].y[ieta]) 
	  << endl;
      
	outfile.close();
	
	cout << "pion <pt> for 0.25 < pt < 1.8 GeV, |y| < 1 = " 
	<< get_meanpt(DATA, 211, 0.25, 1.8, 1, -1, 1) << endl;
	
    }
  else if (mode!=3)
    {
      cerr << "Mode " << mode << " not supported for pseudofreeze = 1.  Exiting.\n";
      exit(1);
    }
  delete[] partid;
  for (int i=0;i<=DATA->NumberOfParticlesToInclude+1;i++) util->char_free(particleList[i].name);
  free(particleList);
}


// returns pseudorapidity for a given rapidity, transverse momentum, and mass
double Freeze::PseudoRap(double y, double pt, double m)
{
  double eta = acosh(
		  2*m*m/pt/pt*sinh(y)*sinh(y)
		  + cosh(2*y)
	       )/2.;
  if (y<0) eta*=-1.;
  return eta;
}

// returns rapidity for a given pseudorapidity, transverse momentum, and mass
double Freeze::Rap(double eta, double pt, double m)
{
  double y = log
  (
    (
      sqrt(m*m + pt*pt*cosh(eta)*cosh(eta))
      + pt*sinh(eta)
    )
    /
    sqrt(m*m+pt*pt)
  );
  return y;
}

// jacobian to convert dN/dy to dN/deta:  dN/deta = dydeta*dN/dy
double Freeze::dydeta(double eta, double pt, double m)
{
  return pt*cosh(eta)/sqrt(m*m + pt*pt*cosh(eta)*cosh(eta));
}



// calculates pt-integrated flow versus pseudorapidity with dn/deta weights.
// format is vn[n][i(real=0 or imaginary part=1)][eta]
// E.g.,
// n=0 i=1,2 is the phi-integrated spectrum, dN/dpt/deta for minpt==maxpt, or dN/deta otherwise
// n=1 real part(i=0) is v_1 cos(Psi_1)
// n=1 imaginary part(i=1) is v_1 sin(Psi_1)
// n=2 i=0 is v_2 cos(Psi_2)
// n=2 i=1 is v_2 sin(Psi_2)
// etc.
// Setting ptmin=ptmax evaluates flow at a fixed transverse momentum, in which case the spectrum is dN/deta/dpt
void Freeze::pt_integrated_flow(InitData *DATA, int number, double minpt, double maxpt, double vn[nharmonics][2][etasize])
{
// 	cout << "Calculating integrated flow for " << minpt << " < p_T > " << maxpt << " for particle " << number << endl;

	//Define index j used in particleList[j]
	int j = partid[MHALF+number];
	double  pt;
  //       double intvn[8][2] = {0};
	int nphi = particleList[j].nphi;
	int npt = particleList[j].npt;
	int neta = particleList[j].ny;
	double intvneta[etasize][nharmonics][2] = {};
	double intvny[etasize][nharmonics][2] = {};
	double m = particleList[j].mass;
	
	if(minpt < particleList[j].pt[0]) 
	{
	  cerr << "Error: called out of range pt in pt_integrated_flow, " 
	  << minpt << " < minimum " << particleList[j].pt[0] << endl;
	  exit(1);
	}
	if(maxpt > particleList[j].pt[npt-1]) 
	{
	  cerr << "Error: called out of range pt in pt_integrated_flow, " 
	  << maxpt << " > maximum " << particleList[j].pt[npt-1] << endl;
	  exit(1);
	}
	if (minpt > maxpt)
	{
	  cerr << "Error in pt_integrated_flow:  minpt must be less than or equal to maxpt\n";
	  exit(1);
	}
	
	
	//loop over pseudorapidity
// 	cout << "ietamax = " << particleList[j].ny << endl;
	for(int ieta=0;ieta<neta;ieta++)
	{
// 	  for(int i = 0;i<8;i++) for(int k =0;k<2;k++) intvn[ieta][i][k]=0;
	  double eta = particleList[j].y[ieta];
	  
	    //Integrate over phi using trapezoid rule
	    for(int iphi=0;iphi<nphi;iphi++) 
	    {
	      
	      // Integrate over pt using gsl
	      double dnetadpt[etasize] = {0};
	      double dnydpt[etasize] = {0};
	      for(int ipt=0;ipt<npt;ipt++) 
	      {
		pt = particleList[j].pt[ipt];
		// jacobian -- dN/deta = jac*dN/dy
		double jac = pt*cosh(eta)/sqrt(m*m + pt*pt*cosh(eta)*cosh(eta));
		dnetadpt[ipt] = pt*jac*particleList[j].dNdydptdphi[ieta][ipt][iphi];
		dnydpt[ipt] = pt*particleList[j].dNdydptdphi[ieta][ipt][iphi];
// 		cout << "spectra = " << dndpt[ipt] << endl;
	      }
	      gsl_interp_accel *ptacc = gsl_interp_accel_alloc ();
	      gsl_interp_accel *yptacc = gsl_interp_accel_alloc ();
	      gsl_spline *ptspline = gsl_spline_alloc (gsl_interp_cspline, npt);
	      gsl_spline *yptspline = gsl_spline_alloc (gsl_interp_cspline, npt);
	      gsl_spline_init (ptspline, particleList[j].pt ,dnetadpt , npt);
	      gsl_spline_init (yptspline, particleList[j].pt ,dnydpt , npt);
	      
	      
	      double dNdeta;
	      if (minpt!=maxpt) dNdeta = gsl_spline_eval_integ(ptspline, minpt, maxpt, ptacc);
	      else dNdeta = gsl_spline_eval(ptspline, minpt, ptacc);
	      double dNdy;
	      if (minpt!=maxpt) dNdy = gsl_spline_eval_integ(yptspline, minpt, maxpt, yptacc);
	      else dNdy = gsl_spline_eval(yptspline, minpt, yptacc);
// 	      cout << "dN = " << dN << endl;
	      
	      double phi = iphi*2*PI/nphi;
	      for(int n = 0;n<nharmonics;n++)
	      {
		intvneta[ieta][n][0] += cos(n*phi)*dNdeta*2*PI/nphi;
		intvneta[ieta][n][1] += sin(n*phi)*dNdeta*2*PI/nphi;
		intvny[ieta][n][0] += cos(n*phi)*dNdy*2*PI/nphi;
		intvny[ieta][n][1] += sin(n*phi)*dNdy*2*PI/nphi;
	      }
	      gsl_spline_free (ptspline);
	      gsl_interp_accel_free (ptacc);
	      gsl_spline_free (yptspline);
	      gsl_interp_accel_free (yptacc);
	    }// phi loop
// 	  for(int i = 1;i<8;i++) for(int k =0;k<2;k++) intvn[ieta][i][k]/=intvn[ieta][0][0];


	  for(int k =0;k<2;k++) 
	  {
	    vn[0][k][ieta] = intvneta[ieta][0][k];
// 	    vn[1][0][k][ieta] = intvny[ieta][0][k];
	  }
	  
	  for(int n = 1;n<nharmonics;n++) for(int k =0;k<2;k++) 
	  {
	    vn[n][k][ieta] = intvneta[ieta][n][k]/intvneta[ieta][0][0];
// 	    vn[1][n][k][ieta] = intvny[ieta][n][k]/intvny[ieta][0][0];
	  }
	  
	  
	}// eta loop

}

// calculates eta= or y-integrated flow versus rapidity.
// format is vn[n][i(real=0 or imaginary part=1)][pt] .
// yflag = 0 takes minrap and maxrap as a range in psuedorapidity,
// while yflag = 1 uses rapidity
// Yield (n=0) is dN/dpt/dy or dN/dpt/deta for minrap==maxrap, or dN/dpt otherwise 
void Freeze::rapidity_integrated_flow(InitData *DATA, int number, int yflag, double minrap, double maxrap, double vn[nharmonics][2][etasize])
{
// 	cout << "Calculating integrated flow for " << miney << " < y < " << maxy << " for particle " << number << endl;

    
	//Define index j used in particleList[j]
	int j = partid[MHALF+number];
// 	double fac, pt;
  //       double intvn[8][2] = {0};
	int nphi = particleList[j].nphi;
	int npt = particleList[j].npt;
	int neta = particleList[j].ny;
	double intvn[ptsize][nharmonics][2] = {};
	double m = particleList[j].mass;
	
	double testmin;
	if(yflag) testmin = Rap(particleList[j].y[0],particleList[j].pt[0],m);
	else testmin = particleList[j].y[0];
	if(minrap < testmin) 
	{
	  cerr << "Error: called out of range rapidity in rap_integrated_flow, " 
	  << minrap << " < minimum " << testmin << endl;
	  exit(1);
	}
	double testmax;
	if(yflag) testmax = Rap(particleList[j].y[neta-1],particleList[j].pt[0],m);
	else testmax = particleList[j].y[neta-1];
	if(maxrap > testmax) 
	{
	  cerr << "Error: called out of range rapidity in rap_integrated_flow, " 
	  << maxrap << " > maximum " << testmax << endl;
	  exit(1);
	}
	if (minrap > maxrap)
	{
	  cerr << "Error in rap_integrated_flow:  minrap must be less than or equal to maxrap\n";
	  exit(1);
	}
	
	
	//loop over pt
// 	cout << "ietamax = " << particleList[j].ny << endl;
	for(int ipt=0;ipt<npt;ipt++)
	{
// 	  for(int i = 0;i<8;i++) for(int k =0;k<2;k++) intvn[ieta][i][k]=0;
	  double pt = particleList[j].pt[ipt];
	  
	    //Integrate over phi using trapezoid rule
	    for(int iphi=0;iphi<nphi;iphi++) 
	    {
	      double ylist[etasize] = {0};
	      
	      // Integrate over pseudorapidity using gsl
	      double dndpt[etasize] = {0};
	      for(int ieta=0;ieta<neta;ieta++) 
	      {
		double eta = particleList[j].y[ieta];
		ylist[ieta] = Rap(eta, pt, m);
		
		if(yflag)dndpt[ieta] = pt*particleList[j].dNdydptdphi[ieta][ipt][iphi];
		else dndpt[ieta] = pt*dydeta(eta, pt, m)*particleList[j].dNdydptdphi[ieta][ipt][iphi];
		
	      }
	      gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	      gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, neta);
	      if(yflag) gsl_spline_init (spline, ylist ,dndpt , neta);
	      else gsl_spline_init (spline, particleList[j].y ,dndpt , neta);

	      
	      double dNdp;
	      if (minrap!=maxrap) dNdp = gsl_spline_eval_integ(spline, minrap, maxrap, acc);
	      else dNdp = gsl_spline_eval(spline, maxrap, acc);
	      
	      
	      double phi = iphi*2*PI/nphi;
	      for(int i = 0;i<nharmonics;i++)
	      {
		intvn[ipt][i][0] += cos(i*phi)*dNdp*2*PI/nphi;
		intvn[ipt][i][1] += sin(i*phi)*dNdp*2*PI/nphi;
	      }
	      gsl_spline_free (spline);
	      gsl_interp_accel_free (acc);
	    }// phi loop


	  for(int k =0;k<2;k++) 
	  {
	    vn[0][k][ipt] = intvn[ipt][0][k];
	  }
	  
	  for(int i = 1;i<nharmonics;i++) for(int k =0;k<2;k++) 
	  {
	    vn[i][k][ipt] = intvn[ipt][i][k]/intvn[ipt][0][0];
	  }
	  
	  
	}// pt loop

}


// calculates pt- and eta-integrated flow for a given range in pt and eta
// format is vn[n][i(real=0 or imaginary part=1)]
// this one has the pt integral nested inside the phi integral inside the eta integral
void Freeze::pt_and_eta_integrated_flow(InitData *DATA, int number, double minpt, double maxpt, double mineta, double maxeta, double vn[nharmonics][2])
{
  int j = partid[MHALF+number];
//   int npt = particleList[j].npt;
  int neta = particleList[j].ny;
  // do pt integration first
//   double vneta[2][8][2][100] = {0};
  
  if(mineta < particleList[j].y[0]) 
  {
    cerr << "Error: called out of range eta in pt_and_eta_integrated_flow, " 
    << mineta << " < minimum " << particleList[j].y[0] << endl;
    exit(1);
  }
  if(maxeta > particleList[j].y[neta-1]) 
  {
    cerr << "Error: called out of range eta in pt_and_eta_integrated_flow, " 
    << maxeta << " > maximum " << particleList[j].y[neta-1] << endl;
    exit(1);
  }
  if (mineta > maxeta)
  {
    cerr << "Error in pt_and_eta_integrated_flow:  mineta must be less than or equal to maxeta\n";
    exit(1);
  }
  
  double vneta[nharmonics][2][etasize] = {};
  pt_integrated_flow(DATA, number, minpt, maxpt, vneta);
  double dndeta[etasize];
  for(int ieta=0;ieta<neta;ieta++) 
  {
    dndeta[ieta] = vneta[0][0][ieta];
  }
  for(int n = 0; n < nharmonics; n++)
  {
    double vncos[etasize] = {0};
    double vnsin[etasize] = {0};
    for(int ieta=0;ieta<neta;ieta++) 
    {
      double weight;
      if (n!=0) weight = dndeta[ieta];
      else weight = 1;
      vncos[ieta] = vneta[n][0][ieta]*weight;
      vnsin[ieta] = vneta[n][1][ieta]*weight;
    }
    gsl_interp_accel *cosacc = gsl_interp_accel_alloc ();
    gsl_spline *cosspline = gsl_spline_alloc (gsl_interp_cspline, neta);
    gsl_spline_init (cosspline, particleList[j].y ,vncos , neta);
    
    gsl_interp_accel *sinacc = gsl_interp_accel_alloc ();
    gsl_spline *sinspline = gsl_spline_alloc (gsl_interp_cspline, neta);
    gsl_spline_init (sinspline, particleList[j].y ,vnsin , neta);
    
    if (mineta!=maxeta)
    {
      vn[n][0] = gsl_spline_eval_integ(cosspline, mineta, maxeta, cosacc);
      vn[n][1] = gsl_spline_eval_integ(sinspline, mineta, maxeta, sinacc);
    }
    else 
    {
      vn[n][0] = gsl_spline_eval(cosspline, mineta, cosacc);
      vn[n][1] = gsl_spline_eval(sinspline, mineta, sinacc);
    }
    
    if(n!=0) for(int i = 0; i<2; i++) vn[n][i]/=vn[0][0];
    
    gsl_spline_free (cosspline);
    gsl_interp_accel_free (cosacc);
    gsl_spline_free (sinspline);
    gsl_interp_accel_free (sinacc);
  }
}

// calculates pt- and (pseudo)rapidity-integrated flow for a given range in pt and eta or y
// format is vn[n][i(real=0 or imaginary part=1)]
//  the rapidity integral is nested inside the phi integral inside the pt integral
void Freeze::pt_and_rapidity_integrated_flow(InitData *DATA, int number, double minpt, double maxpt, int yflag,  double minrap, double maxrap, double vn[nharmonics][2])
{
  int j = partid[MHALF+number];
  int npt = particleList[j].npt;  
  
  if(minpt < particleList[j].pt[0]) 
  {
    cerr << "Error: called out of range pt in pt_and_rapidity_integrated_flow, " 
    << minpt << " < minimum " << particleList[j].pt[0] << endl;
    exit(1);
  }
  if(maxpt > particleList[j].pt[npt-1]) 
  {
    cerr << "Error: called out of range pt in pt_and_rapidity_integrated_flow, " 
    << maxpt << " > maximum " << particleList[j].pt[npt-1] << endl;
    exit(1);
  }
  if (minpt > maxpt)
  {
    cerr << "Error in pt_and_rapidity_integrated_flow:  minpt must be less than or equal to maxpt\n";
    exit(1);
  }
  
  // do  rapidity-integral first
  double vnpt[nharmonics][2][ptsize] = {};
  rapidity_integrated_flow(DATA, number, yflag,  minrap, maxrap, vnpt);
  double dndpt[ptsize];
  for(int ipt=0;ipt<npt;ipt++) 
  {
    dndpt[ipt] = vnpt[0][0][ipt];
  }
  for(int n = 0; n < nharmonics; n++)
  {
    double vncos[ptsize] = {0};
    double vnsin[ptsize] = {0};
    for(int ipt=0;ipt<npt;ipt++) 
    {
      double weight;
      if (n!=0) weight = dndpt[ipt];
      else weight = 1;
      vncos[ipt] = vnpt[n][0][ipt]*weight;
      vnsin[ipt] = vnpt[n][1][ipt]*weight;
    }
    gsl_interp_accel *cosacc = gsl_interp_accel_alloc ();
    gsl_spline *cosspline = gsl_spline_alloc (gsl_interp_cspline, npt);
    gsl_spline_init (cosspline, particleList[j].pt ,vncos , npt);
    
    gsl_interp_accel *sinacc = gsl_interp_accel_alloc ();
    gsl_spline *sinspline = gsl_spline_alloc (gsl_interp_cspline, npt);
    gsl_spline_init (sinspline, particleList[j].pt ,vnsin , npt);
    
    if (minpt!=maxpt)
    {
      vn[n][0] = gsl_spline_eval_integ(cosspline, minpt, maxpt, cosacc);
      vn[n][1] = gsl_spline_eval_integ(sinspline, minpt, maxpt, sinacc);
    }
    else 
    {
      vn[n][0] = gsl_spline_eval(cosspline, minpt, cosacc);
      vn[n][1] = gsl_spline_eval(sinspline, minpt, sinacc);
    }
    
    if(n!=0) for(int i = 0; i<2; i++) vn[n][i]/=vn[0][0];
    
    gsl_spline_free (cosspline);
    gsl_interp_accel_free (cosacc);
    gsl_spline_free (sinspline);
    gsl_interp_accel_free (sinacc);
  }
}

// Return yield in specified range of phase space
double Freeze::get_N(InitData *DATA, int number, double minpt, double maxpt, int yflag, double minrap, double maxrap)
{
  double vn[nharmonics][2];
  pt_and_rapidity_integrated_flow(DATA, number, minpt, maxpt, yflag, minrap, maxrap, vn);
  return vn[0][0];
}

// Return yield in specified range of phase space
// if mineta==maxeta, returns dN/eta.  If minpt==maxpt, dN/dpt.  If both, dN/dpt/deta.  Otherwise, total yield N
double Freeze::get_yield(InitData *DATA, int number, double minpt, double maxpt, double mineta, double maxeta)
{
  double vn[nharmonics][2];
  pt_and_eta_integrated_flow(DATA, number, minpt, maxpt, mineta, maxeta, vn);
  return vn[0][0];
}


// Return <pt> in specified range of phase space
double Freeze::get_meanpt(InitData *DATA, int number, double minpt, double maxpt, int yflag, double minrap, double maxrap)
{
  int j = partid[MHALF+number];
  int npt = particleList[j].npt;
//   int neta = particleList[j].ny;
  
  
  if(minpt < particleList[j].pt[0]) 
  {
    cerr << "Error: called out of range pt in get_meanpt, " 
    << minpt << " < minimum " << particleList[j].pt[0] << endl;
    exit(1);
  }
  if(maxpt > particleList[j].pt[npt-1]) 
  {
    cerr << "Error: called out of range pt in get_meanpt, " 
    << maxpt << " > maximum " << particleList[j].pt[npt-1] << endl;
    exit(1);
  }
  if (minpt > maxpt)
  {
    cerr << "Error in get_meanpt:  minpt must be less than or equal to maxpt\n";
    exit(1);
  }
  

  double ptdndpt[ptsize];
  double dndpt[ptsize];
  for(int ipt=0;ipt<npt;ipt++) 
  {
    double pt = particleList[j].pt[ipt];
    dndpt[ipt] = get_N(DATA, number, pt, pt, yflag, minrap, maxrap);
    ptdndpt[ipt] = pt*dndpt[ipt];
  }
    gsl_interp_accel *numacc = gsl_interp_accel_alloc ();
    gsl_spline *numspline = gsl_spline_alloc (gsl_interp_cspline, npt);
    gsl_spline_init (numspline, particleList[j].pt ,ptdndpt , npt);
    
    gsl_interp_accel *denacc = gsl_interp_accel_alloc ();
    gsl_spline *denspline = gsl_spline_alloc (gsl_interp_cspline, npt);
    gsl_spline_init (denspline, particleList[j].pt ,dndpt , npt);
    
    double num = gsl_spline_eval_integ(numspline, minpt, maxpt, numacc);
    double den = gsl_spline_eval_integ(denspline, minpt, maxpt, denacc);
    
    gsl_spline_free (numspline);
    gsl_interp_accel_free (numacc);
    gsl_spline_free (denspline);
    gsl_interp_accel_free (denacc);
    
    return num/den;
  
}


// Return v_n for specified range of phase space
double Freeze::get_vn(InitData *DATA, int number, double minpt, double maxpt, int yflag, double minrap, double maxrap, int n)
{
  if( n > nharmonics-1 ) 
  {
    cout << "harmonic (" << n << ") too large.  Must increase nharmonics in freeze_pseudo.cpp\n";
    exit(1);
  }
  double vn[nharmonics][2];
  pt_and_rapidity_integrated_flow(DATA, number, minpt, maxpt, yflag, minrap, maxrap, vn);
  return sqrt(vn[n][0]*vn[n][0] + vn[n][1]*vn[n][1]);
}

// Return Psi_n for specified range of phase space
double Freeze::get_psi_n(InitData *DATA, int number, double minpt, double maxpt, int yflag, double minrap, double maxrap, int n)
{
  if( n > nharmonics-1 ) 
  {
    cout << "harmonic (" << n << ") too large.  Must increase nharmonics in freeze_pseudo.cpp\n";
    exit(1);
  }
  double vn[nharmonics][2];
  pt_and_rapidity_integrated_flow(DATA, number, minpt, maxpt, yflag, minrap, maxrap, vn);
  return atan2(vn[n][1],vn[n][0])/n;
}

// Return charged hadron yield in specified range of phase space
// You'll almost certainly want to set yflat=0 to calculate over a fixed range in eta
// if minrap==maxrap, returns dN/eta.  If minpt==maxpt, dN/dpt.  If both, dN/dpt/deta.  Otherwise, total yield N
double Freeze::get_Nch(InitData *DATA, double minpt, double maxpt, int yflag, double minrap, double maxrap)
{
  if(particleMax<18)
  {
    cout << "Cannot compute charged hadron yield.  \
	  Spectra for all charged hadrons have not been computed. \
	  particleMax = " << particleMax << endl;
    exit(1);
  }
  double N=0;
  int chargedhd[6] = {1,3,4,5,17,18};
  for ( int k=0; k<6; k++ )
    {
      int i = chargedhd[k];
      int number = particleList[i].number;
//       double vn[nharmonics][2];
//       pt_and_eta_integrated_flow(DATA, number, minpt, maxpt, mineta, maxeta, vn);
//       N+= vn[0][0];
      N+= get_N(DATA, number, minpt, maxpt, yflag, minrap, maxrap);
    }
    return N;
}

// Return charged hadron vn for specified range of phase space
double Freeze::get_vn_ch(InitData *DATA, double minpt, double maxpt, int yflag, double minrap, double maxrap, int n)
{
  if(particleMax<18)
  {
    cout << "Cannot compute charged hadron vn.  \
	  Spectra for all charged hadrons have not been computed. \
	  particleMax = " << particleMax << endl;
    exit(1);
  }
  if( n > nharmonics-1 ) 
  {
    cout << "harmonic (" << n << ") too large.  Must increase nharmonics in freeze_pseudo.cpp\n";
    exit(1);
  }
  double numr=0.;//real part (x projection, \sum N*v_n*cos(Psi_n))
  double numi=0.;//imaginary part (y projection, \sum N*v_n*sin(Psi_n))
  double den=0.;// denominator (\sum N)
  int chargedhd[6] = {1,3,4,5,17,18};
  for ( int k=0; k<6; k++ )
    {
      int i = chargedhd[k];
      int number = particleList[i].number;
      double vn[nharmonics][2];
      pt_and_rapidity_integrated_flow(DATA, number, minpt, maxpt, yflag,  minrap, maxrap, vn);
      numr+= vn[0][0]*vn[n][0];
      numi+= vn[0][0]*vn[n][1];
      den+= vn[0][0];
    }
  return sqrt(numr*numr+numi*numi)/den;
}

// Return charged hadron Psi_n for specified range of phase space
double Freeze::get_psi_n_ch(InitData *DATA, double minpt, double maxpt, int yflag, double minrap, double maxrap, int n)
{
  if(particleMax<18)
  {
    cout << "Cannot compute charged hadron event plane.  \
	  Spectra for all charged hadrons have not been computed. \
	  particleMax = " << particleMax << endl;
    exit(1);
  }
  if( n > nharmonics-1 ) 
  {
    cout << "harmonic (" << n << ") too large.  Must increase nharmonics in freeze_pseudo.cpp\n";
    exit(1);
  }
  double numr=0.;//real part (x projection, \sum N*v_n*cos(Psi_n))
  double numi=0.;//imaginary part (y projection, \sum N*v_n*sin(Psi_n))
  int chargedhd[6] = {1,3,4,5,17,18};
  for ( int k=0; k<6; k++ )
    {
      int i = chargedhd[k];
      int number = particleList[i].number;
      double vn[nharmonics][2];
      pt_and_rapidity_integrated_flow(DATA, number, minpt, maxpt, yflag, minrap, maxrap, vn);
      numr+= vn[0][0]*vn[n][0];
      numi+= vn[0][0]*vn[n][1];
    }
  return atan2(numi, numr)/n;
}



// Calculates v1 and psi_1 with pt-dependent weight (pt - <pt^2>/<pt>) in specified range of phase space
// ch == 0 calculates for particle "number", while ch > 0 calculates for all charged hadrons
void Freeze::weighted_v1(InitData *DATA, int number, double minpt, double maxpt, int yflag,  double minrap, double maxrap, double vn[2], int ch)
{
  int j = partid[MHALF+number];
  int npt = particleList[j].npt;
//   int neta = particleList[j].ny;
  
  
  if(minpt < particleList[j].pt[0]) 
  {
    cerr << "Error: called out of range pt in weighted_v1, " 
    << minpt << " < minimum " << particleList[j].pt[0] << endl;
    exit(1);
  }
  if(maxpt > particleList[j].pt[npt-1]) 
  {
    cerr << "Error: called out of range pt in weighted_v1, " 
    << maxpt << " > maximum " << particleList[j].pt[npt-1] << endl;
    exit(1);
  }
  if (minpt > maxpt)
  {
    cerr << "Error in weighted_v1:  minpt must be less than or equal to maxpt\n";
    exit(1);
  }
  
  
  // start by finding <pt> and <pt^2>
  double ptptdndpt[ptsize];
  double ptdndpt[ptsize];
  double dndpt[ptsize];
  for(int ipt=0;ipt<npt;ipt++) 
  {
    double pt = particleList[j].pt[ipt];
    if(ch) dndpt[ipt] = get_Nch(DATA, pt, pt, yflag, minrap, maxrap);
    else dndpt[ipt] = get_N(DATA, number, pt, pt, yflag, minrap, maxrap);
    ptdndpt[ipt] = pt*dndpt[ipt];
    ptptdndpt[ipt] = pt*pt*dndpt[ipt];
  }
  
  
    gsl_interp_accel *numacc = gsl_interp_accel_alloc ();
    gsl_spline *numspline = gsl_spline_alloc (gsl_interp_cspline, npt);
    gsl_spline_init (numspline, particleList[j].pt ,ptdndpt , npt);
    
    gsl_interp_accel *num2acc = gsl_interp_accel_alloc ();
    gsl_spline *num2spline = gsl_spline_alloc (gsl_interp_cspline, npt);
    gsl_spline_init (num2spline, particleList[j].pt, ptptdndpt , npt);
    
    gsl_interp_accel *denacc = gsl_interp_accel_alloc ();
    gsl_spline *denspline = gsl_spline_alloc (gsl_interp_cspline, npt);
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
    for(int ipt=0;ipt<npt;ipt++) 
    {
      double pt = particleList[j].pt[ipt];
      double weight = (pt - meanpt2/meanpt)*dndpt[ipt];
      
      double v1;
      double psi1;
      if(ch) v1 = get_vn_ch(DATA, pt, pt, yflag, minrap, maxrap, 1);
      else v1 = get_vn(DATA, number, pt, pt, yflag, minrap, maxrap, 1);
      if(ch) psi1 = get_psi_n_ch(DATA, pt, pt, yflag, minrap, maxrap, 1);
      else psi1 = get_psi_n(DATA, number, pt, pt, yflag, minrap, maxrap, 1);
      v1cos[ipt] = weight*v1*cos(psi1);
      v1sin[ipt] = weight*v1*sin(psi1);
      v1den[ipt] = weight;
    }
    
  
//     gsl_interp_accel *numacc = gsl_interp_accel_alloc ();
//     gsl_spline *numspline = gsl_spline_alloc (gsl_interp_cspline, npt);
    gsl_spline_init (numspline, particleList[j].pt ,v1cos , npt);
    
//     gsl_interp_accel *num2acc = gsl_interp_accel_alloc ();
//     gsl_spline *num2spline = gsl_spline_alloc (gsl_interp_cspline, npt);
    gsl_spline_init (num2spline, particleList[j].pt, v1sin , npt);
    
//     gsl_interp_accel *denacc = gsl_interp_accel_alloc ();
//     gsl_spline *denspline = gsl_spline_alloc (gsl_interp_cspline, npt);
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

double Freeze::get_weighted_v1(InitData *DATA, int number, double minpt, double maxpt, int yflag, double minrap, double maxrap, int ch)
{
  double vn[2];
  weighted_v1(DATA, number, minpt, maxpt, yflag, minrap, maxrap, vn, ch);
  return vn[0];
}

double Freeze::get_weighted_psi1(InitData *DATA, int number, double minpt, double maxpt, int yflag, double minrap, double maxrap, int ch)
{
  double vn[2];
  weighted_v1(DATA, number, minpt, maxpt, yflag, minrap, maxrap, vn, ch);
  return vn[1];
}



// Output yield dN/ptdydpt and v_n at eta=0 as a function of pT
void Freeze::OutputDifferentialFlowAtMidrapidity(InitData *DATA, int number, int full) 
{
  
	//Define index j used in particleList[j]
	int j = partid[MHALF+number];
	double m = particleList[j].mass;
// 	int nphi = particleList[j].nphi;
	int npt = particleList[j].npt;

// 	double minpt = particleList[j].pt[0];
// 	double maxpt = particleList[j].pt[npt-1];
	double eta = 0.0;
	
	cout << "Calculating flow at midrapidity for particle " << number << endl;

// 	for (int iphi=0;iphi<nphi;iphi++) phipbuff[iphi] = iphi*2*PI/nphi;

	//Set output file name
	string fname;
	stringstream tmpStr;
	fname="./outputs/";
	if (full) {
		fname+="F";
	}
	fname+="vnpteta02-";
	tmpStr << number;
	fname+=tmpStr.str();
	fname+=".dat";	
	
	string fname2;
	stringstream tmpStr2;
	fname2="./outputs/";
	if (full) {
		fname2+="F";
	}
	fname2+="vnpteta03-";
	tmpStr2 << number;
	fname2+=tmpStr.str();
	fname2+=".dat";	
	
		//Open output file for vn
	ofstream outfilevn;
	outfilevn.open(fname.c_str());
	
	ofstream outfilevn2;
	outfilevn2.open(fname2.c_str());

	//Set the format of the output
	//outfile.width (10);
// 	outfile.precision(6);
// 	outfile.setf(ios::scientific);
// 	outfilevn.precision(6);
// 	outfilevn.setf(ios::scientific);
	outfilevn << "#pt\tdN/ptdYdptdphi\tv1cos\tv1sin\tv2cos\tv2sin\tv3cos\tv3sin\tv4cos\tv4sin\tv5cos\tv5sin\tv6cos\tv6sin\tv7cos\tv7sin\n";
	
// 	outfilevn.precision(6);
// 	outfilevn.setf(ios::scientific);
	outfilevn2 << "#pt\tdN/ptdYdptdphi\tv1cos\tv1sin\tv2cos\tv2sin\tv3cos\tv3sin\tv4cos\tv4sin\tv5cos\tv5sin\tv6cos\tv6sin\tv7cos\tv7sin\n";

	double vn[nharmonics][2];

	//Loop over pT
// 	cout << "npt = " << npt << endl;
	for(int ipt=0;ipt<npt;ipt++) {
		double pt=particleList[j].pt[ipt];
// 		cout << "pt = " << pt << endl;
		
// 		pt_and_eta_integrated_flow(DATA, number, pt, pt, eta, eta, vn);
// 		pt_and_y_integrated_flow(DATA, number, pt, pt, 0., 0., vn);
		pt_and_rapidity_integrated_flow(DATA, number, pt, pt, 1, 0., 0., vn);

		//Output result
		outfilevn << pt;
		outfilevn << "\t" << vn[0][0]/pt;
		for(int i = 1;i<nharmonics;i++) for(int k =0;k<2;k++) outfilevn << "\t" << vn[i][k];
		outfilevn << endl;
		
// 		pt_and_eta_integrated_flow2(DATA, number, pt, pt, eta, eta, vn);
		pt_and_rapidity_integrated_flow(DATA, number, pt, pt, 0, eta, eta, vn);
		
		outfilevn2 << pt;
		outfilevn2 << "\t" << vn[0][0]/dydeta(eta, pt, m)/pt;
		for(int i = 1;i<nharmonics;i++) for(int k =0;k<2;k++) outfilevn2 << "\t" << vn[i][k];
		outfilevn2 << endl;
	}
	

	//Close file
 	outfilevn.close();
	outfilevn2.close();
// 	}

}

// Output yield dN/ptdydpt and v_n for |y|<1 as a function of pT
void Freeze::OutputDifferentialFlowNearMidrapidity(InitData *DATA, int number, int full) 
{
  
	//Define index j used in particleList[j]
	int j = partid[MHALF+number];
	double m = particleList[j].mass;
	int npt = particleList[j].npt;

	double mineta = -1.;
	double maxeta = 1.;
	
	cout << "Calculating flow near midrapidity for particle " << number << endl;

// 	for (int iphi=0;iphi<nphi;iphi++) phipbuff[iphi] = iphi*2*PI/nphi;

	//Set output file name
	string fname;
	stringstream tmpStr;
	fname="./outputs/";
	if (full) {
		fname+="F";
	}
	fname+="vnptTPC2-";
	tmpStr << number;
	fname+=tmpStr.str();
	fname+=".dat";	
	
	string fname2;
	stringstream tmpStr2;
	fname2="./outputs/";
	if (full) {
		fname2+="F";
	}
	fname2+="vnptTPC3-";
	tmpStr2 << number;
	fname2+=tmpStr.str();
	fname2+=".dat";	
	
		//Open output file for vn
	ofstream outfilevn;
	outfilevn.open(fname.c_str());
	
	ofstream outfilevn2;
// 	outfilevn2.open(fname2.c_str());

	//Set the format of the output
	//outfile.width (10);
// 	outfile.precision(6);
// 	outfile.setf(ios::scientific);
// 	outfilevn.precision(6);
// 	outfilevn.setf(ios::scientific);
	outfilevn << "#pt\tdN/ptdYdptdphi\tv1cos\tv1sin\tv2cos\tv2sin\tv3cos\tv3sin\tv4cos\tv4sin\tv5cos\tv5sin\tv6cos\tv6sin\tv7cos\tv7sin\n";
	
// 	outfilevn.precision(6);
// 	outfilevn.setf(ios::scientific);
// 	outfilevn2 << "#pt\tdN/ptdYdptdphi\tv1cos\tv1sin\tv2cos\tv2sin\tv3cos\tv3sin\tv4cos\tv4sin\tv5cos\tv5sin\tv6cos\tv6sin\tv7cos\tv7sin\n";

	double vn[nharmonics][2];

	//Loop over pT
// 	cout << "npt = " << npt << endl;
	for(int ipt=0;ipt<npt;ipt++) {
		double pt=particleList[j].pt[ipt];
// 		cout << "pt = " << pt << endl;
		
// 		pt_and_eta_integrated_flow(DATA, number, pt, pt, mineta, maxeta, vn);
		double miny = Rap(mineta,pt,m);
		double maxy = Rap(maxeta,pt,m);
// 		pt_and_y_integrated_flow(DATA, number, pt, pt, miny, maxy, vn);
		pt_and_rapidity_integrated_flow(DATA, number, pt, pt, 1, miny, maxy, vn);

		//Output result
		outfilevn << pt;
		outfilevn << "\t" << vn[0][0]/pt/(maxy-miny);
		for(int i = 1;i<nharmonics;i++) for(int k =0;k<2;k++) outfilevn << "\t" << vn[i][k];
		outfilevn << endl;
/*		
		pt_and_eta_integrated_flow2(DATA, number, pt, pt, mineta, maxeta, vn);
		
		outfilevn2 << pt;
		outfilevn2 << "\t" << vn[0][0]/dydeta(eta, pt, m)/pt/(maxeta-mineta);
		for(int i = 1;i<nharmonics;i++) for(int k =0;k<2;k++) outfilevn2 << "\t" << vn[i][k];
		outfilevn2 << endl;*/
	}
	

	//Close file
 	outfilevn.close();
// 	outfilevn2.close();
// 	}

}


// Output yield dN/deta and v_n integrated over pT, as a function of pseudorapidity
void Freeze::OutputIntegratedFlow(InitData *DATA, int number, int full) 
{
  
	//Define index j used in particleList[j]
	int j = partid[MHALF+number];
	int neta = particleList[j].ny;
	int npt = particleList[j].npt;

	double minpt = particleList[j].pt[0];
	double maxpt = particleList[j].pt[npt-1];
// 	double eta = 0.0;
	
	cout << "Calculating flow integrated between " << minpt << " < p_T < " << maxpt << " vs. pseudorapidity for particle " << number << endl;

// 	for (int iphi=0;iphi<nphi;iphi++) phipbuff[iphi] = iphi*2*PI/nphi;

	//Set output file name
	string fname;
	stringstream tmpStr;
	fname="./outputs/";
	if (full) {
		fname+="F";
	}
	fname+="vneta2-";
	tmpStr << number;
	fname+=tmpStr.str();
	fname+=".dat";	
	
	string fname2;
	stringstream tmpStr2;
	fname2="./outputs/";
	if (full) {
		fname2+="F";
	}
	fname2+="vneta3-";
	tmpStr2 << number;
	fname2+=tmpStr.str();
	fname2+=".dat";	
	
		//Open output file for vn
	ofstream outfilevn;
	outfilevn.open(fname.c_str());
	
	ofstream outfilevn2;
	outfilevn2.open(fname2.c_str());

	//Set the format of the output
	//outfile.width (10);
// 	outfile.precision(6);
// 	outfile.setf(ios::scientific);
// 	outfilevn.precision(6);
// 	outfilevn.setf(ios::scientific);
	outfilevn << "#eta\tdN/ptdYdptdphi\tv1cos\tv1sin\tv2cos\tv2sin\tv3cos\tv3sin\tv4cos\tv4sin\tv5cos\tv5sin\tv6cos\tv6sin\tv7cos\tv7sin\n";
	
// 	outfilevn.precision(6);
// 	outfilevn.setf(ios::scientific);
	outfilevn2 << "#eta\tdN/ptdYdptdphi\tv1cos\tv1sin\tv2cos\tv2sin\tv3cos\tv3sin\tv4cos\tv4sin\tv5cos\tv5sin\tv6cos\tv6sin\tv7cos\tv7sin\n";

	double vn[nharmonics][2];

	//Loop over pT
// 	cout << "npt = " << npt << endl;
	for(int ieta=0;ieta<neta;ieta++) {
		double eta=particleList[j].y[ieta];
// 		cout << "pt = " << pt << endl;
		
// 		pt_and_eta_integrated_flow(DATA, number, minpt, maxpt, eta, eta, vn);
		pt_and_rapidity_integrated_flow(DATA, number, minpt, maxpt, 0, eta, eta, vn);

		//Output result
		outfilevn << eta;
		for(int i = 0;i<nharmonics;i++) for(int k =0;k<2;k++) outfilevn << "\t" << vn[i][k];
		outfilevn << endl;
		
// 		pt_and_eta_integrated_flow2(DATA, number, minpt, maxpt, eta, eta, vn);
		pt_and_rapidity_integrated_flow(DATA, number, minpt, maxpt, 1, eta, eta, vn);
		
		outfilevn2 << eta;
		for(int i = 0;i<nharmonics;i++) for(int k =0;k<2;k++) outfilevn2 << "\t" << vn[i][k];
		outfilevn2 << endl;
	}
	

	//Close file
 	outfilevn.close();
	outfilevn2.close();
// 	}

}

//Output pT and phi integrated charged hadron spectra as a function of eta
// Function by J-F 
void Freeze::Output_charged_hadrons_eta_differential_spectra(InitData *DATA, int full, const int * hadron_list, int nb_hadrons) {
	//
	int nb_hadron = nb_hadrons;
	
	int j, number, nphi, npt, neta;
	double eta;
	double m;
	double pt;
	double tmp_dNdeta, * tmp_dNdetadpT, * tmp_ptList;

	//Make sure all the necessary spectra are available
	for(int ihadron=0;ihadron<nb_hadron;ihadron++) {

		j=partid[MHALF+hadron_list[ihadron]];

		if (j >= particleMax) {
			cout << "Can't compute charged hadron spectra, some spectra are not available\n";
			exit(1);
		}
	}

	//Set output file name
	stringstream fname;
	//stringstream tmpStr;
	fname << "./outputs/";
	if (full) {
		fname <<"F";
	}
	fname <<"dNdeta";
	for(int ihadron=0;ihadron<nb_hadron;ihadron++) fname << "_" << hadron_list[ihadron]; 
	fname << ".dat";
	
	//Open output file for vn
	ofstream outfile;
	outfile.open(fname.str().c_str());

	//Set the format of the output
	outfile.precision(6);
	outfile.setf(ios::scientific);
	outfile << "#pt\tdN/deta\n";

	//Assume all particles have the same discretization in pT, eta and phi
	j=partid[MHALF+hadron_list[0]];
	nphi = particleList[j].nphi;
	npt = particleList[j].npt;
	neta = particleList[j].ny;

	tmp_dNdetadpT = (double *) malloc(npt*sizeof(double));
	tmp_ptList = (double *) malloc(npt*sizeof(double));

	gsl_interp * interp_pT = gsl_interp_alloc(gsl_interp_cspline, npt);
	gsl_interp_accel * accel_pT = gsl_interp_accel_alloc();

	//cout << "Calculating dN/dy for charged hadrons" << endl;
	
	for(int ieta=0;ieta<neta;ieta++) {

		eta = particleList[j].y[ieta];

		tmp_dNdeta=0.0;

		//Negative and zero indices are mesons and antimesons, even positive indices are baryons, odd positive indices are antibaryons
		for(int ihadron=0;ihadron<nb_hadron;ihadron++) {

			//Define index j used in particleList[j]
			number=hadron_list[ihadron];
			j = partid[MHALF+number];
			m = particleList[j].mass;

			//cout << "Adding hadron " << j << " with mass=" << m << "..." << endl;
			
			//Loop over pT
			for(int ipt=0;ipt<npt;ipt++) {
				pt=particleList[j].pt[ipt];
			
				tmp_dNdetadpT[ipt]=0.0;
				tmp_ptList[ipt]=pt;
				//jacobian to switch from dN/dY to dN/deta
				double jac = pt*cosh(eta)/sqrt(m*m + pt*pt*cosh(eta)*cosh(eta));
				//double jac = 1.;
				
				//Integrate over phi using trapezoid rule
				for(int iphi=0;iphi<nphi;iphi++) {
					tmp_dNdetadpT[ipt]+=pt*jac*particleList[j].dNdydptdphi[ieta][ipt][iphi]*(2*M_PI)/(nphi);
					//tmp_dNdetadpT[ipt]+=exp(-1*particleList[j].pt[ipt])*pow(cos(iphi*2*M_PI/(nphi)),2)*(2*M_PI)/(nphi);
				}
			}

			//Integrate in pT with GSL
			gsl_interp_init(interp_pT, tmp_ptList, tmp_dNdetadpT, npt);
			tmp_dNdeta+=gsl_interp_eval_integ(interp_pT, tmp_ptList, tmp_dNdetadpT, DATA->min_pt, DATA->max_pt, accel_pT);

		}

		//Output result
		outfile << eta << "\t" << tmp_dNdeta << "\n";

	}

	free(tmp_dNdetadpT);
	free(tmp_ptList);
	gsl_interp_accel_reset(accel_pT);
	gsl_interp_free(interp_pT);

	//Close file
	outfile.close();

}

//Output midrapidity, phi integrated hadron spectra as a function of pT
void Freeze::Output_midrapidity_hadrons_spectra(InitData *DATA, int full, const int * hadron_list, int nb_hadrons) {
	//
	int nb_hadron = nb_hadrons;
	
	int j, number, nphi, npt, neta;
	double eta;
	double m;
	double pt;
	double tmp_dNdy, tmp_dNdeta, tmp;

	//Make sure all the necessary spectra are available
	for(int ihadron=0;ihadron<nb_hadron;ihadron++) {

		j=partid[MHALF+hadron_list[ihadron]];

		if (j >= particleMax) {
			cout << "Can't compute charged hadron spectra, some spectra are not available\n";
			exit(1);
		}
	}

	//Set output file name
	stringstream fname;
	//stringstream tmpStr;
	fname << "./outputs/";
	if (full) {
		fname <<"F";
	}
	fname <<"spectra_eta0";
	for(int ihadron=0;ihadron<nb_hadron;ihadron++) fname << "_" << hadron_list[ihadron]; 
	fname << ".dat";
	
	//Open output file for vn
	ofstream outfile;
	outfile.open(fname.str().c_str());

	//Set the format of the output
	outfile.precision(6);
	outfile.setf(ios::scientific);
	outfile << "#pt\t1/pT dN/dpt dy\t1/pT dN/dpt deta\n";

	//Assume all particles have the same discretization in pT, eta and phi
	j=partid[MHALF+hadron_list[0]];
	nphi = particleList[j].nphi;
	npt = particleList[j].npt;
	neta = particleList[j].ny;

	//cout << "Calculating dN/dy for charged hadrons" << endl;
	int ieta=neta/2;
	eta = particleList[j].y[ieta];
	
	//Loop over pT
	for(int ipt=0;ipt<npt;ipt++) {

		pt=particleList[j].pt[ipt];

		tmp_dNdy=0.0;
		tmp_dNdeta=0.0;

		for(int ihadron=0;ihadron<nb_hadron;ihadron++) {

			//Define index j used in particleList[j]
			number=hadron_list[ihadron];
			j = partid[MHALF+number];
			m = particleList[j].mass;

			//cout << "Adding hadron " << j << " with mass=" << m << "..." << endl;
			
			//jacobian to switch from dN/dY to dN/deta
			double jac = pt*cosh(eta)/sqrt(m*m + pt*pt*cosh(eta)*cosh(eta));
			//double jac = 1.;
			
			//Integrate over phi using trapezoid rule
			for(int iphi=0;iphi<nphi;iphi++) {
				tmp=particleList[j].dNdydptdphi[ieta][ipt][iphi]*(2*M_PI)/(nphi);
				tmp_dNdy+=tmp;
				tmp_dNdeta+=jac*tmp;
				//tmp_dNdetadpT[ipt]+=exp(-1*particleList[j].pt[ipt])*pow(cos(iphi*2*M_PI/(nphi)),2)*(2*M_PI)/(nphi);
			}
		}

		//Output result
		outfile << pt << "\t" << tmp_dNdy << "\t" << tmp_dNdeta << "\n";

	}

	//Close file
	outfile.close();

}
