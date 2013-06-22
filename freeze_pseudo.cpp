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
  int bytes_read;
  double  neta;
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
	  particleList[ip].pt[ipt] =  ptmin + (ptmax - ptmin)*pow(ipt,2)/pow(iptmax-1,2); // power law
// 	  particleList[ip].pt[ipt] =  0.01 + ipt*ptmax/(iptmax-1); // compare to UVH2+1
// 	  particleList[ip].pt[ipt] = gala15x[ipt]/12; // gauss laguerre abissas
	}
      for (int ieta=0; ieta<=pseudo_steps; ieta++ )
	{
	  particleList[ip].y[ieta] =  ieta*deltaeta-etamax; // store pseudorapidity here
	  //	  if (ip==1) cout << "read particleList[ip].y[" << ieta << "] = " <<  particleList[ip].y[ieta] << endl;
	}
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
		  bytes_read=fscanf(s_file, "%lf", &particleList[ip].dNdydptdphi[ieta][ipt][iphi]);
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
  double f, T, mu, tau, x, eta, E, sign, delta_f;
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
  int s = particleList[j].strange;
  double c = particleList[j].charge;
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
	      pt =  ptmin + (ptmax - ptmin)*pow(ipt,2)/pow(iptmax-1,2); // power law
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

  

// --------------------- resonance decays. adapted from azhydro ------------------------------------------------


/*************************************************
*
*	Edndp3
*
* 
**************************************************/
// This function interpolates the needed spectra for a given y, pt and phi.
// Uses trilinear interpolation

double Freeze::Edndp3_pseudo(double yr, double ptr, double phirin, int res_num)
/* 				/\* supersedes during test the right one *\/ */
/* 	double	yr;		/\* y  of resonance *\/ */
/* 	double	ptr;		/\* pt of resonance *\/ */
/* 	double	phirin;		/\* phi angle  of resonance *\/ */
/* 	int	res_num;	/\* Montecarlo number of resonance 	*\/ */
{
  double	phir, val, val1, val2;
  double        f1, f2, f1s, f2s;
  int     	pn, ny, npt, nphi;
  
  if(phirin < 0.0){
    printf("ERROR: phir %15.8le < 0 !!! \n", phirin);exit(0);}
  if(phirin > 2.0*PI){
    printf("ERROR: phir %15.8le > 2PI !!! \n", phirin);exit(0);}
  phir= phirin;
//  if(phirin < 0.5*PI) phir = phirin;
//  else{
//    if(phirin < PI) phir = PI - phirin;
//    else{
//      if(phirin < 1.5*PI) phir = phirin - PI;
//      else phir = 2.0*PI - phirin;
//    }
//  }
  
  pn = partid[MHALF + res_num];

  double yrtemp = yr;
  // If pseudofreeze flag is set,  dNdydptdphi is on a fixed grid in pseudorapidity. 
  // Set yr to the *pseudorapidity* of the resonance, and then interpolate the yield
  // at that value.
  yr = PseudoRap(yrtemp, ptr, particleList[pn].mass);
  
  if (yr < -particleList[pn].ymax || yr > particleList[pn].ymax)
    {
      //      fprintf(stderr,"yr=%f out of range ymax=%f\n", yr,particleList[pn].ymax);

      return 0.;
    }

  nphi = 1; 
  while((phir > phiArray[nphi])&&(nphi<(particleList[pn].nphi-1))) nphi++; 
  npt = 1; 
  while((ptr > particleList[pn].pt[npt]) && npt<(particleList[pn].npt - 1)) npt++; 
  ny = 1; 
  while((yr > particleList[pn].y[ny]) && ny<(particleList[pn].ny - 1)) ny++; 

  /* phi interpolation */
  f1 = util->lin_int(phiArray[nphi-1], phiArray[nphi], 
	       particleList[pn].dNdydptdphi[ny-1][npt-1][nphi-1], 
	       particleList[pn].dNdydptdphi[ny-1][npt-1][nphi], phir);
  f2 = util->lin_int(phiArray[nphi-1], phiArray[nphi], 
	       particleList[pn].dNdydptdphi[ny-1][npt][nphi-1], 
	       particleList[pn].dNdydptdphi[ny-1][npt][nphi], phir);

  if (f1<0.) f1=0.; // security: if for some reason we got a negative number of particles (happened in the viscous code at large eta sometimes)
  if (f2<0.) f2=0.;

  f1s=f1;
  f2s=f2;

  if(ptr > PTCHANGE && f1!=0 && f2!=0){
    f1 = log(f1); 
    f2 = log(f2);
  }
  val1 = util->lin_int(particleList[pn].pt[npt-1],particleList[pn].pt[npt], 
		f1, f2, ptr);

  if(ptr > PTCHANGE && f1!=0 && f2!=0)
    val1 = exp(val1);
  
  

  if (isnan(val1))
    {
      fprintf(stderr,"\n number=%d\n\n",res_num);
      //      fprintf(stderr,"val=%f\n",val);
      fprintf(stderr,"val1=%f\n",val1);
      fprintf(stderr,"val2=%f\n",val2);
      fprintf(stderr,"f1=%f\n",f1);
      fprintf(stderr,"f2=%f\n",f2);
      fprintf(stderr,"f1s=%f\n",f1s);
      fprintf(stderr,"f2s=%f\n",f2s);
      
      fprintf(stderr,"pn=%d\n",pn);
      fprintf(stderr,"ny=%d\n",ny);
      fprintf(stderr,"npt=%d\n",npt);
      fprintf(stderr,"nphi=%d\n",nphi);
      fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny-1][npt-1][nphi-1]);
      fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny][npt-1][nphi-1]);
      fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny-1][npt][nphi-1]);
      fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny-1][npt-1][nphi]);
      fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny][npt][nphi-1]);
      fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny-1][npt][nphi]);
      fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny][npt-1][nphi]);
      fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny][npt][nphi]);
      fprintf(stderr,"phi1=%f\n",phiArray[nphi-1]);
      fprintf(stderr,"phi2=%f\n",phiArray[nphi]);
      fprintf(stderr,"pt1=%f\n",particleList[pn].pt[npt-1]);
      fprintf(stderr,"pt2=%f\n",particleList[pn].pt[npt]);
      fprintf(stderr,"y1=%f\n",particleList[pn].y[ny-1]);
      fprintf(stderr,"y2=%f\n",particleList[pn].y[ny]);
  }

  f1 = util->lin_int(phiArray[nphi-1], phiArray[nphi], 
	       particleList[pn].dNdydptdphi[ny][npt-1][nphi-1], 
	       particleList[pn].dNdydptdphi[ny][npt-1][nphi], phir);
  f2 = util->lin_int(phiArray[nphi-1], phiArray[nphi], 
	       particleList[pn].dNdydptdphi[ny][npt][nphi-1], 
	       particleList[pn].dNdydptdphi[ny][npt][nphi], phir);
  
  if (f1<0.) f1=0.; // security: if for some reason we got a negative number of particles (happened in the viscous code at large eta sometimes)
  if (f2<0.) f2=0.;

  if(ptr > PTCHANGE && f1!=0 && f2!=0){
    f1 = log(f1); 
    f2 = log(f2);
  }
  val2 = util->lin_int(particleList[pn].pt[npt-1],particleList[pn].pt[npt], 
		f1, f2, ptr);
  if(ptr > PTCHANGE && f1!=0 && f2!=0)
    val2 = exp(val2);
  
  val = util->lin_int(particleList[pn].y[ny-1],particleList[pn].y[ny],val1,val2,yr);

  /*
    printf(" nphi  %i npt %i \n", nphi,npt);
    printf(" f1  %15.8le %15.8le  \n", f1, f2);
    printf(" phi  %15.8lf %15.8lf  \n", phiArray[nphi-1], phiArray[nphi]); 
    printf(" pt   %15.8lf %15.8lf  \n", particleList[pn].pt[npt-1],particleList[pn].pt[npt]);
    printf(" phi  %15.8lf pt %15.8lf    val %15.8lf \n", phir, ptr,val); 
    printf(" phi %15.8le %15.8le \n",particleList[pn].dNdydptdphi[npt][nphi-1],
    particleList[pn].dNdydptdphi[npt][nphi]);
    printf(" pt  %15.8le %15.8le \n",particleList[pn].dNdydptdphi[npt-1][nphi-1],
    particleList[pn].dNdydptdphi[npt-1][nphi]);
    
    exit(0);
  */
  if (isnan(val))
    {
      fprintf(stderr,"val=%f\n",val);
      fprintf(stderr,"val1=%f\n",val1);
      fprintf(stderr,"val2=%f\n",val2);
      fprintf(stderr,"f1=%f\n",f1);
      fprintf(stderr,"f2=%f\n",f2);
      fprintf(stderr,"f1s=%f\n",f1s);
      fprintf(stderr,"f2s=%f\n",f2s);
      fprintf(stderr,"ny=%d\n",ny);
      fprintf(stderr,"npt=%d\n",npt);
      fprintf(stderr,"nphi=%d\n",nphi);
      fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny-1][npt-1][nphi-1]);
      fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny][npt-1][nphi-1]);
      fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny-1][npt][nphi-1]);
      fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny-1][npt-1][nphi]);
      fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny][npt][nphi-1]);
      fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny-1][npt][nphi]);
      fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny][npt-1][nphi]);
      fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny][npt][nphi]);
      fprintf(stderr,"phi1=%f\n",phiArray[nphi-1]);
      fprintf(stderr,"phi2=%f\n",phiArray[nphi]);
      fprintf(stderr,"pt1=%f\n",particleList[pn].pt[npt-1]);
      fprintf(stderr,"pt2=%f\n",particleList[pn].pt[npt]);
      fprintf(stderr,"y1=%f\n",particleList[pn].y[ny-1]);
      fprintf(stderr,"y2=%f\n",particleList[pn].y[ny]);
      fprintf(stderr,"yR=%f\n",yr);
    }
  
/*   if (val2>10*val1) */
/*     { */
/*       fprintf(stderr,"y1=%f\n",particleList[pn].y[ny-1]); */
/*       fprintf(stderr,"y2=%f\n",particleList[pn].y[ny]); */
/*       fprintf(stderr,"yR=%f\n",yr); */
/*       fprintf(stderr,"val1=%f\n",val1); */
/*       fprintf(stderr,"val2=%f\n",val2); */
/*       fprintf(stderr,"val=%f\n",val); */
/*       fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny][npt-1][nphi-1]); */
/*       fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny][npt][nphi-1]); */
/*       fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny][npt-1][nphi]); */
/*       fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny][npt][nphi]); */
/*     } */

//  fprintf(stderr,"yr=%f, ptr=%f, phir=%f, val=%e\n",yr,ptr,phirin,val);
 /*  if(ptr<0.5 && fabs(yr)<0.5) */
/*     { */
/*       fprintf(stderr,"\n number=%d\n\n",res_num); */
/*       fprintf(stderr,"val=%f\n",val); */
/*       fprintf(stderr,"val1=%f\n",val1); */
/*       fprintf(stderr,"val2=%f\n",val2); */
/*       fprintf(stderr,"f1=%f\n",f1); */
/*       fprintf(stderr,"f2=%f\n",f2); */
/*       fprintf(stderr,"f1s=%f\n",f1s); */
/*       fprintf(stderr,"f2s=%f\n",f2s); */
      
/*       fprintf(stderr,"pn=%d\n",pn); */
/*       fprintf(stderr,"ny=%d\n",ny); */
/*       fprintf(stderr,"npt=%d\n",npt); */
/*       fprintf(stderr,"nphi=%d\n",nphi); */
/*       fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny-1][npt-1][nphi-1]); */
/*       fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny][npt-1][nphi-1]); */
/*       fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny-1][npt][nphi-1]); */
/*       fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny-1][npt-1][nphi]); */
/*       fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny][npt][nphi-1]); */
/*       fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny-1][npt][nphi]); */
/*       fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny][npt-1][nphi]); */
/*       fprintf(stderr,"dN..=%e\n",particleList[pn].dNdydptdphi[ny][npt][nphi]); */
/*       fprintf(stderr,"phi1=%f\n",phiArray[nphi-1]); */
/*       fprintf(stderr,"phi2=%f\n",phiArray[nphi]); */
/*       fprintf(stderr,"pt1=%f\n",particleList[pn].pt[npt-1]); */
/*       fprintf(stderr,"pt2=%f\n",particleList[pn].pt[npt]); */
/*       fprintf(stderr,"yr=%f\n",yr); */
/*       fprintf(stderr,"y1=%f\n",particleList[pn].y[ny-1]); */
/*       fprintf(stderr,"y2=%f\n",particleList[pn].y[ny]); */
/*     }  */
  return val;
}


double Freeze::dnpir2N_pseudo (double phi, void *para1)    
{
  pblockN *para = (pblockN *) para1;
  double D;
  double eR, plR, ptR, yR, phiR, sume, jac;
  double cphiR, sphiR;
  double dnr;			/* dn/mtdmt of resonance */
  
  sume = para->e + para->e0;

  D = para->e * para->e0 + para->pl * para->p0 * para->costh +
    para->pt * para->p0 * para->sinth * cos (phi) + para->m1 * para->m1;

  eR = para->mr * (sume * sume / D - 1.0);
  jac = para->mr + eR;
  plR = para->mr * sume * (para->pl - para->p0 * para->costh) / D;
  ptR = (eR * eR - plR * plR - para->mr * para->mr);

  if (ptR < 0.0)
    ptR = 0.0;

  else
    ptR = sqrt (ptR);

  yR = 0.5 * log ((eR + plR) / (eR - plR));
  cphiR = -jac * (para->p0 * para->sinth * cos (phi + para->phi)
		  - para->pt * cos (para->phi)) / (sume * ptR);
  sphiR = -jac * (para->p0 * para->sinth * sin (phi + para->phi)
		  - para->pt * sin (para->phi)) / (sume * ptR);

  if ((fabs (cphiR) > 1.000) || (fabs (sphiR) > 1.000))
    {
      if ((fabs (cphiR) > 1.01) || (fabs (sphiR) > 1.01))
	{
	  //  printf ("  |phir| = %15.8lf  > 1 ! \n", phiR);
	  printf (" phi %15.8le D %15.8le \n", phi, D);
	  printf (" eR %15.8le plR %15.8le \n", eR, plR);
	  printf (" ptR %15.8le jac %15.8le \n", ptR, jac);
	  printf (" sume %15.8le costh %15.8le \n", sume, para->costh);

	  printf (" pt %15.8le \n", para->pt);
	  printf (" mt  %15.8le \n", para->mt);
	  printf (" y %15.8le \n", para->y);
	  printf (" e %15.8le \n", para->e);
	  printf (" e0 %15.8le \n", para->e0);
	  printf (" p0 %15.8le \n", para->p0);
	  printf (" pl %15.8le \n", para->pl);
	  printf (" phi %15.8le \n", para->phi);

	  printf (" m1 %15.8le \n", para->m1);
	  printf (" m2 %15.8le \n", para->m2);
	  printf (" m3 %15.8le \n", para->m3);
	  printf (" mr %15.8le \n", para->mr);
	  if (cphiR > 1.0)
	    cphiR = 1.0;
	  if (cphiR < -1.0)
	    cphiR = -1.0;
	  //exit (0);
	}
      else
	{
	  if (cphiR > 1.0)
	    cphiR = 1.0;
	  if (cphiR < -1.0)
	    cphiR = -1.0;
	}
    }

  phiR = acos (cphiR);
  if (sphiR < 0.0)
    phiR = 2.0 * PI - phiR;

  dnr = Edndp3_pseudo (yR, ptR, phiR, para->res_num);

  /*printf(" phir = %15.8lf  ! ", phiR);
     printf(" ptR %15.8le jac %15.8le ", ptR, jac );
     printf(" dnr %15.8le \n", dnr); */

  return dnr * jac * jac / (2.0 * sume * sume);
}

double Freeze::dnpir1N_pseudo (double costh, void* para1)	       
{
  pblockN *para = (pblockN *) para1;
  double r;
  para->costh = costh;
  para->sinth = sqrt (1.0 - para->costh * para->costh);
  r = gauss (PTN2, &Freeze::dnpir2N_pseudo, 0.0, 2.0 * PI, para); //Integrates the "dnpir2N" kernel over phi using gaussian integration
  return r;
}

double Freeze::dn2ptN_pseudo (double w2, void* para1)
{
  pblockN *para = (pblockN *) para1;
  para->e0 = (para->mr * para->mr + para->m1 * para->m1 - w2) / (2 * para->mr); //particle one energy in resonance rest frame
  para->p0 = sqrt (para->e0 * para->e0 - para->m1 * para->m1); // particle one absolute value of three momentum on resonance rest frame
  return gauss (PTN1, &Freeze::dnpir1N_pseudo, -1.0, 1.0, para); //Integrate the "dnpir1N" kernel over cos(theta) using gaussian integration
}

double Freeze::dn3ptN_pseudo (double x, void* para1)  //The integration kernel for "W" in 3-body decays. x=invariant mass of other particles squared
{
  pblockN *para = (pblockN *) para1;
  double e0 =(para->mr * para->mr + para->m1 * para->m1 - x) / (2 * para->mr);
  double p0 = sqrt (e0 * e0 - para->m1 * para->m1);
  double a = (para->m2 + para->m3) * (para->m2 + para->m3);
  double b = (para->m2 - para->m3) * (para->m2 - para->m3);
  double re = p0 * sqrt ((x - a) * (x - b)) / x * dn2ptN_pseudo (x, para);
  return re;
}


/********************************************************************
*
*	Edndp3_2bodyN()
*
* transverse momentum spectrum in GeV^-2 from pions out of resonances
*********************************************************************/
double Freeze::Edndp3_2bodyN_pseudo (double y, double pt, double phi, double m1, double m2, double mr, int res_num)
/* 		/\* in units of GeV^-2,includes phasespace and volume, */
/* 		   does not include degeneracy factors  *\/ */
/*      double y;			/\* rapidity of particle 1       *\/ */
/*      double pt;			/\* transverse momentum of particle 1    *\/ */
/*      double phi;		/\* phi angle of particle 1      *\/ */
/*      double m1, m2;		/\* restmasses of decay particles in MeV *\/ */
/*      double mr;			/\* restmass of resonance MeV            *\/ */
/*      int res_num;		/* Montecarlo number of the Resonance   */ 

{
  double mt = sqrt (pt * pt + m1 * m1);
  double norm2;			/* 2-body normalization         */
  pblockN para;
  double res2;

  para.pt = pt;
  para.mt = mt;
  para.e = mt * cosh (y);
  para.pl = mt * sinh (y);
  para.y = y;
  para.phi = phi;
  para.m1 = m1;
  para.m2 = m2;
  para.mr = mr;

  para.res_num = res_num;

  norm2 = 1.0 / (2.0 * PI);
  res2 = norm2 * dn2ptN_pseudo (m2 * m2, &para); //Calls the integration routines for 2-body
  if (res2<0.) res2=0.;
  return res2;			/* like Ed3ndp3_2body() */
}


/********************************************************************
*
*	Edndp3_3bodyN()
*
* transverse momentum spectrum in GeV^-2 from pions out of resonances
*********************************************************************/
double Freeze::Edndp3_3bodyN_pseudo (double y, double pt, double phi, double m1, double m2,
		      double m3, double mr, double norm3, int res_num)
		/* in units of GeV^-2,includes phasespace and volume,
		   does not include degeneracy factors  */
{
  double mt = sqrt (pt * pt + m1 * m1);
  pblockN para;
  double wmin, wmax;
  double res3;
  double slope;			/* slope of resonance for high mt */
  int pn;

  para.pt = pt;
  para.mt = mt;
  para.y = y;
  para.e = mt * cosh (y);
  para.pl = mt * sinh (y);
  para.phi = phi;

  para.m1 = m1;
  para.m2 = m2;
  para.m3 = m3;
  para.mr = mr;

  pn = partid[MHALF + res_num];

  para.res_num = res_num;

  wmin = (m2 + m3) * (m2 + m3); 
  wmax = (mr - m1) * (mr - m1);
  res3 = 2.0 * norm3 * gauss (PTS4, &Freeze::dn3ptN_pseudo, wmin, wmax, &para) / mr;  //Integrates "W" using gaussian 
  if (res3<0.) res3=0.; 
  return res3;
}


/**************************************************************************
*									  *
*	add_reso							  *
*									  *
* computes the pt, mt distribution including resonance decays		  *
***************************************************************************/


void Freeze::add_reso_pseudo(int pn, int pnR, int k, int j, int pseudofreeze)
{
  nblock paranorm;		/* for 3body normalization integral */
  double y, eta;
  double m1, m2, m3, mr;
  double norm3;			/* normalisation of 3-body integral */
  int pn2, pn3, pn4;		/* internal numbers for resonances */
  int part;
  int l, i, n;
  int neta, npt, nphi;

  neta = particleList[pn].ny; // for pseudofreeze, this variable stores number of points in pseudorapidity
  npt = particleList[pn].npt;
  nphi = particleList[pn].nphi;
  double deltaphi = 2*PI/nphi;


  // Determine the number of particles involved in the decay with the switch
  switch (abs (decay[j].numpart))
    {
    case 1: //Only 1 particle, if it gets here, by accident, this prevents any integration for 1 particle chains
      break;

    case 2: // 2-body decay 
      {
	if (k == 0)
	  pn2 = partid[MHALF + decay[j].part[1]];

	else
	  pn2 = partid[MHALF + decay[j].part[0]];

	//printf ("case 2:  i %3i j %3i k %3i \n", pn, j, k);
	m1 = particleList[pn].mass;
	m2 = particleList[pn2].mass;
	mr = particleList[pnR].mass;

	while ((m1 + m2) > mr)
	  {
	    mr += 0.25 * particleList[pnR].width;
	    m1 -= 0.5 * particleList[pn].width;
	    m2 -= 0.5 * particleList[pn2].width;
	  }
// 	fprintf(stderr,"mr=%f\n",mr);
// 	fprintf(stderr,"m1=%f\n",m1);
// 	fprintf(stderr,"m2=%f\n",m2);
		
	for (n = 0; n < neta; n++)
	  {
	    eta = particleList[pn].y[n];// for pseudofreeze, this variable stores pseudorapidity
	    for (l = 0; l < npt; l++)
	      {
		y = Rap(eta,particleList[pn].pt[l],m1);
		for (i = 0; i < nphi; i++)
		  {
		    double phi = i*deltaphi;
		    double spectrum = Edndp3_2bodyN_pseudo (y, particleList[pn].pt[l], phi, m1, m2, mr, particleList[pnR].number);
		    if (isnan(spectrum)
			)
		      //	||Edndp3_2bodyN (y, particleList[pn].pt[l], phiArray[i], m1, m2, mr, particleList[pnR].number)<0
		      {
			fprintf(stderr,"2 pt=%f\n",particleList[pn].pt[l]);
			fprintf(stderr,"2 number=%d\n",particleList[pnR].number);
			fprintf(stderr,"2 Edn..=%f\n", spectrum);
		      }
		    else
		      // Call the 2-body decay integral and add its contribution to the daughter particle of interest
		      {   
			particleList[pn].dNdydptdphi[n][l][i] += decay[j].branch *
			  spectrum; 
			if(n==neta/2 && i==0) {
			  // fprintf(stderr,"m1=%f, m2=%f, mr=%f, pnR=%d\n",m1,m2,mr,pnR);

			  particleList[pn].resCont[n][l][i]+= decay[j].branch *
			  spectrum; 
// 			  fprintf(stderr," %d %f %e %e %e %e\n", n, eta, particleList[pn].pt[l], decay[j].branch *
// 				  spectrum,
// 				  particleList[pn].dNdydptdphi[n][l][i]-decay[j].branch *
// 				  spectrum,particleList[pn].resCont[n][l][i]); 
			}
		      }
		  }
	      }
	  }
	break;
      }

    case 3: //3-body decay
      {
	if (k == 0)
	  {
	    pn2 = partid[MHALF + decay[j].part[1]];
	    pn3 = partid[MHALF + decay[j].part[2]];
	  }
	else
	  {
	    if (k == 1)
	      {
		pn2 = partid[MHALF + decay[j].part[0]];
		pn3 = partid[MHALF + decay[j].part[2]];
	      }
	    else
	      {
		pn2 = partid[MHALF + decay[j].part[0]];
		pn3 = partid[MHALF + decay[j].part[1]];
	      }
	  }
	
	m1 = particleList[pn].mass;
	m2 = particleList[pn2].mass;
	m3 = particleList[pn3].mass;
	mr = particleList[pnR].mass;
	paranorm.a = (mr + m1) * (mr + m1);
	paranorm.b = (mr - m1) * (mr - m1);
	paranorm.c = (m2 + m3) * (m2 + m3);
	paranorm.d = (m2 - m3) * (m2 - m3);
	norm3 = mr * mr / (2 * PI * gauss (PTS3, &Freeze::norm3int, paranorm.c,
					   paranorm.b, &paranorm));
	
	// printf("case 3:  i %3i j %3i k %3i \n",pn,j,k); 
	
	for (n = 0; n < neta; n++)
	  {
	    eta = particleList[pn].y[n];// for pseudofreeze, this variable stores pseudorapidity
	    for (l = 0; l < npt; l++)
	      {
// 		if(pseudofreeze)  y = Rap(particleList[pn].y[n],particleList[pn].pt[l],m1);
		y = Rap(eta,particleList[pn].pt[l],m1);
		for (i = 0; i < nphi; i++)
		  {
		    double phi = i*deltaphi;
		    double spectrum = Edndp3_3bodyN(y, particleList[pn].pt[l], phi,
					    m1, m2, m3, mr, norm3, particleList[pnR].number);
		    if (isnan(spectrum))
		      {
			fprintf(stderr,"3 number=%d\n",particleList[pnR].number);
			// Call the 3-body decay integral and add its contribution to the daughter particle of interest 
			fprintf(stderr,"3 Edn..=%f\n",   spectrum);
		      }
		    else
		      {
			particleList[pn].dNdydptdphi[n][l][i] += decay[j].branch *
		      	spectrum;
		      }
		    
		    if(n==neta/2 && i==0)
		      {
			//	fprintf(stderr,"m1=%f, m2=%f, m3=%f, mr=%f, pnR=%d\n",m1,m2,m3,mr,pnR);
			particleList[pn].resCont[n][l][i]+= decay[j].branch *
			  spectrum;
		       
/*			fprintf(stderr,"%d %f %e %e %e %e\n",n,eta, particleList[pn].pt[l], decay[j].branch *
				spectrum,
				particleList[pn].dNdydptdphi[n][l][i]-decay[j].branch *
				spectrum,particleList[pn].resCont[n][l][i]);*/ 
		      }
		  }
	      }
	  }
	break;
      }

    case 4: //4-body decay (rare and low contribution)
      {
	if (k == 0)
	  {
	    pn2 = partid[MHALF + decay[j].part[1]];
	    pn3 = partid[MHALF + decay[j].part[2]];
	    pn4 = partid[MHALF + decay[j].part[3]];
	  }
	else
	  {
	    if (k == 1)
	      {
		pn2 = partid[MHALF + decay[j].part[0]];
		pn3 = partid[MHALF + decay[j].part[2]];
		pn4 = partid[MHALF + decay[j].part[3]];
	      }
	    else
	      {
		if (k == 2)
		  {
		    pn2 = partid[MHALF + decay[j].part[0]];
		    pn3 = partid[MHALF + decay[j].part[1]];
		    pn4 = partid[MHALF + decay[j].part[3]];
		  }
		else
		  {
		    pn2 = partid[MHALF + decay[j].part[0]];
		    pn3 = partid[MHALF + decay[j].part[1]];
		    pn4 = partid[MHALF + decay[j].part[2]];
		  }
	      }
	  }
	//approximate the 4-body with a 3-body decay with the 4th particle being the center of mass of 2 particles.
	m1 = particleList[pn].mass;
	m2 = particleList[pn2].mass;
	mr = particleList[pnR].mass;
	m3 = 0.5 * (particleList[pn3].mass + particleList[pn4].mass + mr - m1 - m2);
	paranorm.a = (mr + m1) * (mr + m1);
	paranorm.b = (mr - m1) * (mr - m1);
	paranorm.c = (m2 + m3) * (m2 + m3);
	paranorm.d = (m2 - m3) * (m2 - m3);
	norm3 = mr * mr / (2 * PI * gauss (PTS3, &Freeze::norm3int, paranorm.c,
					   paranorm.b, &paranorm));
	// printf("case 3:  i %3i j %3i k %3i \n",pn,j,k); 
	
	
	for (n = 0; n < neta; n++)
	  {
	    eta = particleList[pn].y[n];// for pseudofreeze, this variable stores pseudorapidity
	    for (i = 0; i < nphi; i++)
	      {
		double phi = i*deltaphi;
		for (l = 0; l < npt; l++)
		  {
// 		    if(pseudofreeze)  y = Rap(particleList[pn].y[n],particleList[pn].pt[l],m1);
		    y = Rap(eta,particleList[pn].pt[l],m1);
		    double spectrum = Edndp3_3bodyN(y, particleList[pn].pt[l], phi,
					    m1, m2, m3, mr, norm3, particleList[pnR].number);
		     if (isnan(spectrum))
		      {
			fprintf(stderr,"3 number=%d\n",particleList[pnR].number);
			// Call the 3-body decay integral and add its contribution to the daughter particle of interest 
			particleList[pn].dNdydptdphi[n][l][i] += decay[j].branch *
			  spectrum;
		      }
		     else
		      particleList[pn].dNdydptdphi[n][l][i] += decay[j].branch *
			spectrum;
		     //fprintf(stderr,"4 Edn..=%f\n", Edndp3_2bodyN (y, particleList[pn].pt[l], phiArray[i],
		    //					m1, m2, mr, particleList[pnR].number));
		    // the 4-body decay approximated by the 3-body decay routine
		  }
	      }
	  }
	break;
      }
      
    default:
      printf ("ERROR in add_reso! \n");
      printf ("%i decay not implemented ! \n", abs (decay[j].numpart));
      exit (0);
    }
}

void Freeze::cal_reso_decays_pseudo(int maxpart, int maxdecay, int bound, int mode, int pseudofreeze)
{

  int i, j, k, l, ll;
  int pn, pnR, pnaR;
  int part,n1,n2,n3,ny,npt,nphi;
  
  fprintf (stderr," CALCULATE RESONANCE DECAYS (as fast as I can) \n");
  pn = partid[MHALF + bound];
  if (mode==5) pn = maxpart-1;

  ny = particleList[pn].ny;
  npt = particleList[pn].npt;
  nphi = particleList[pn].nphi;
  
  for(i=maxpart-1;i > pn-1;i--)  //Cycle the particles known from the particle.dat input
    {

      for (n1 = 0; n1 < ny; n1++)
	{
	  for (n2 = 0; n2 < npt; n2++)
	    {
	      for (n3 = 0; n3 < nphi; n3++)
		{
		  particleList[pn].resCont[n1][n2][n3] =0.;
		}
	    }
	}

      part = particleList[i].number;
      fprintf (stderr,"Calculating the decays with ");
      fprintf (stderr,"%s \n", particleList[i].name);
      //fprintf (stderr,"%i, %i, b=%d\n", part, maxdecay,particleList[i].baryon);
      switch (particleList[i].baryon)  // Check to see whether or not the particle is baryon, anti-baryon or meson
	{
	case 1: //Baryon
	  {
	    //fprintf(stderr,"Is a baryon. \n");
	   
	    for (j = 0; j < maxdecay; j++) // Cycle through every decay channel known (as given in resoweak.dat)
	      {                            // to see if the particle was a daughter particle in a decay channel
		pnR = partid[MHALF + decay[j].reso];
		//fprintf(stderr,"Partid is %i.\n",pnR);
		for (k = 0; k < abs (decay[j].numpart); k++)
		  {
		    if ((part == decay[j].part[k]) && (decay[j].numpart != 1))// Make sure that the decay channel isn't trivial
		      {                                                       // and contains the daughter particle
			fprintf(stderr,"Partid is %i. %s into %s \n",pnR,particleList[pnR].name,particleList[i].name);
			//fprintf(stderr,"Calculating a decay \n");
			add_reso_pseudo (i, pnR, k, j, pseudofreeze);
		      }
		  }
	      }
	    break;
	  }

	case -1: //Anti-Baryon
	  {
	    //fprintf(stderr,"Is an anti-baryon.\n");
	    for (j = 0; j < maxdecay; j++)// Cycle through every decay channel known (as given in resoweak.dat)
	      {                            // to see if the particle was a daughter particle in a decay channel
		pnaR = partid[MHALF - decay[j].reso];
		//if (pnaR==-1) continue;
		for (k = 0; k < abs (decay[j].numpart); k++)
		  {
		    if ((-part == decay[j].part[k]) && (decay[j].numpart != 1))// Make sure that the decay channel isn't trivial
		      {                                                        // and contains the daughter particle
			fprintf(stderr,"Partid is %i. %s into %s \n",pnaR,particleList[pnaR].name,particleList[i].name);
			//fprintf(stderr,"Calculating a decay \n");
			add_reso_pseudo (i, pnaR, k, j, pseudofreeze);
		      }
		  }
	      }
	    break;
	  }

	case 0:// Meson
	  {
	    //fprintf(stderr,"Is a meson. \n");
	    
	    for (j = 0; j < maxdecay; j++)
	      {
		pnR = partid[MHALF + decay[j].reso];
		//fprintf(stderr,"Partid is %i.\n",pnR);
		for (k = 0; k < abs (decay[j].numpart); k++)
		  {
		    if (particleList[pnR].baryon == 1)
		      {
			pnaR = partid[MHALF - decay[j].reso];
			if ((particleList[i].charge == 0)
			    && (particleList[i].strange == 0))
			  {
			    if ((part == decay[j].part[k])
				&& (decay[j].numpart != 1))
			      {
				fprintf(stderr,"Partid is %i, %s into %s \n",pnR,particleList[pnR].name,particleList[i].name);
				fprintf(stderr,"and %i, %s into %s \n",pnaR,particleList[pnaR].name,particleList[i].name);
				//fprintf(stderr,"Calculating a decay \n");
				add_reso_pseudo (i, pnR, k, j, pseudofreeze);
				add_reso_pseudo (i, pnaR, k, j, pseudofreeze);
			      }
			  }
			else
			  {
			    if ((part == decay[j].part[k])
				&& (decay[j].numpart != 1))
			      {
				fprintf(stderr,"Partid is %i, %s into %s \n",pnR,particleList[pnR].name,particleList[i].name);
				//fprintf(stderr,"Calculating a decay \n");
				add_reso_pseudo (i, pnR, k, j, pseudofreeze);
			      }
			    if ((-part == decay[j].part[k])
				&& (decay[j].numpart != 1))
			      {
				fprintf(stderr,"Partid is %i, %s into %s \n",pnaR,particleList[pnaR].name,particleList[i].name);
				//fprintf(stderr,"Calculating a decay \n");
				add_reso_pseudo (i, pnaR, k, j, pseudofreeze);
			      }
			  }
		      }
		    else
		      {
			if ((part == decay[j].part[k])
			    && (decay[j].numpart != 1))
			  {
			    fprintf(stderr,"Partid is %i, %s into %s \n",pnR,particleList[pnR].name,particleList[i].name);
			    //fprintf(stderr,"Calculating a decay \n");
			    add_reso_pseudo (i, pnR, k, j, pseudofreeze);
			  }
		      }
		  }
	      }
	    break;
	  }

	  
	default:
	  fprintf (stderr,"Error in switch in func partden_wdecay \n");
	  exit (0);

	}
    }
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
      int ret;
      if (rank == 0) ret = system("rm yptphiSpectra.dat yptphiSpectra?.dat yptphiSpectra??.dat particleInformation.dat 2> /dev/null");
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
		ComputeParticleSpectrum_pseudo(DATA, number, b, size, rank);

  
		// Wait until all processors are finished and concatonate results
		// (Don't want any processors to start writing the next particle to 
		// file until the concatonation is done)
		MPI_Barrier(MPI_COMM_WORLD);
		
		if(rank==0) 
		{
		  ret = system("cat yptphiSpectra?.dat >> yptphiSpectra.dat");
		  ret = system("cat yptphiSpectra??.dat >> yptphiSpectra.dat 2> /dev/null");
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
	  ComputeParticleSpectrum_pseudo(DATA, number,  b,  size, rank);
	  
	  // send something just to make sure that rank 0 waits for all others to be done:
	  int check;
	  check=rank;
	  
	  if (rank > 0)
	    MPI::COMM_WORLD.Send(&check,1,MPI::INT,0,1);
	    
	  if (rank == 0)
	    {
	      for (int from=1; from < size; from ++)
		MPI::COMM_WORLD.Recv(&check,1,MPI::INT,from,1);
	      
	      ret = system("cat yptphiSpectra?.dat >> yptphiSpectra.dat");
	      ret = system("cat yptphiSpectra??.dat >> yptphiSpectra.dat 2> /dev/null");
/*	      
	      ReadSingleSpectrum(DATA);
	      cout << "output results..." << endl;
	      OutputFullParticleSpectrum2(DATA, number, 4., b, 0);*/
	    }
	}
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
      if (rank == 0) int ret = system("rm FyptphiSpectra.dat FparticleInformation.dat 2> /dev/null");
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


      //Set output file name for total multiplicity
      string fname;
      stringstream tmpStr;
      fname="./outputs/NCMS.dat";	
      
      ofstream outfile;

      outfile.open(fname.c_str(),ios::trunc);
      
      int neta = particleList[i].ny;
      
//       double dndeta[etasize] = {0};
      double dNchdeta[100] = {0};
      int chargedhd[6] = {1,3,4,5,17,18};
      for ( int k=0; k<6; k++ )
	{
	  i = chargedhd[k];
	  number = particleList[i].number;
	  if(particleMax>=i)
	  {
	    int npt = particleList[i].npt;
	    double vn[8][2];
	    for(int ieta=0; ieta<neta; ieta++)
	    {
	      double eta = particleList[i].y[ieta];
	      double minpt = particleList[i].pt[0];
	      double maxpt = particleList[i].pt[npt-1];
	      pt_and_eta_integrated_flow(DATA,number,minpt,maxpt,eta,eta,vn);
	      dNchdeta[ieta]+=vn[0][0];
	    }
	    pt_and_eta_integrated_flow(DATA,number,0.5,particleList[i].pt[npt-1],-2.4,2.4,vn);
	    double N = vn[0][0];
	    outfile << number << "\t" << N << endl;
	  }
	}
	for(int ieta=0; ieta<neta; ieta++) cout << "eta, yield = " << particleList[i].y[ieta] << ", " << dNchdeta[ieta] << endl;


	outfile.close();
	
    }
  else if (mode==14) // take tabulated post-decay spectra and compute various observables and integrated quantities
    {
      ReadSpectra_pseudo(DATA, 1, 1);
	const int stable_charged_hadron_list[] = {211,-211,321,-321,2212,-2212};
	Output_charged_hadrons_eta_differential_spectra(DATA, 1, stable_charged_hadron_list,sizeof(stable_charged_hadron_list)/sizeof(int));
	//Output_charged_hadrons_eta_differential_spectra(DATA, 1);

//       for ( i=1; i<particleMax; i++ )
      int chargedhd[6] = {1,3,4,5,17,18};
      

      //Set output file name for total multiplicity
      string fname;
      stringstream tmpStr;
      fname="./outputs/FNchCMS.dat";	
      
      ofstream outfile;

      outfile.open(fname.c_str(),ios::trunc);

      //Set the format of the output
//       outfile.precision(4);
//       outfile.setf(ios::scientific);
      
      double N = 0;
      for ( int k=0; k<6; k++ )
	{
	  i = chargedhd[k];
	  number = particleList[i].number;
	  if(particleMax>=i)
	  {
// 	    OutputDifferentialFlowAtMidrapidity(DATA, number,1);
  //  	  OutputIntegratedFlowForCMS(DATA, number,1);
	    N += OutputYieldForCMS(DATA, number,1);
	  }
	}
      cout << "Nch = " << N << endl;
      outfile << N << endl;
      outfile.close();
      
      int npt = particleList[1].npt;
      cout << "Nch from get_Nch = " << get_Nch(DATA, 0.4, particleList[1].pt[npt-1], -8,8) << endl;
      
      cout << "Nch in CMS acceptance (e.g., for p-Pb paper) = " << get_Nch(DATA, 0.4, particleList[1].pt[npt-1], -2.4,2.4) << endl;
      
      int neta = particleList[1].ny;

     double dNchdeta[100] = {0};
//       int chargedhd[6] = {1,3,4,5,17,18};
      for ( int k=0; k<6; k++ )
	{
	  i = chargedhd[k];
	  number = particleList[i].number;
	  if(particleMax>=i)
	  {

	    int npt = particleList[i].npt;
	    double vn[8][2];
	    for(int ieta=0; ieta<neta; ieta++)
	    {
	      double eta = particleList[i].y[ieta];
	      double minpt = particleList[i].pt[0];
	      double maxpt = particleList[i].pt[npt-1];
	      pt_and_eta_integrated_flow(DATA,number,minpt,maxpt,eta,eta,vn);
	      dNchdeta[ieta]+=vn[0][0];
	    }
	    pt_and_eta_integrated_flow(DATA,number,0.5,particleList[i].pt[npt-1],-2.4,2.4,vn);
	    double N = vn[0][0];
	    outfile << number << "\t" << N << endl;
	  }
	}
// 	for(int ieta=0; ieta<neta; ieta++) cout << "eta, yield = " << particleList[i].y[ieta] << ", " << dNchdeta[ieta] << endl;
	for(int ieta=0; ieta<neta; ieta++) cout << "eta, yield = " << particleList[i].y[ieta] << ", " << get_Nch(DATA,particleList[i].pt[0],particleList[i].pt[npt-1],particleList[1].y[ieta],particleList[1].y[ieta]) << endl;

      
      
    }
  else if (mode!=3)
    {
      cerr << "Mode " << mode << " not supported for pseudofreeze = 1.  Exiting.\n";
      exit(1);
    }
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

const int nharmonics = 8; // calculate up to maximum harmonic (n-1) -- for nharmonics = 8, calculate from v_0 o v_7
const int etasize = 100; // max number of points in eta for array
const int ptsize = 100; // max number of points in pt for array


// calculates pt-integrated flow versus pseudorapidity with dn/deta weights.
// format is vn[n][i(real=0 or imaginary part=1)][eta] .
// n=0 i=1,2 is the phi-integrated spectrum, dN/deta.
// n=1 real part(i=0) is v_1 cos(Psi_1)
// n=1 imaginary part(i=1) is v_1 sin(Psi_1)
// n=2 i=0 is v_2 cos(Psi_2)
// n=2 i=1 is v_2 sin(Psi_2)
// etc.
// Setting ptmin=ptmax evaluates flow at a fixed transverse momentum, in which case the spectrum is dN/deta/dpt
// void Freeze::pt_integrated_flow(InitData *DATA, int number, double minpt, double maxpt, double ****vn)
void Freeze::pt_integrated_flow(InitData *DATA, int number, double minpt, double maxpt, double vn[nharmonics][2][etasize])
{
// 	cout << "Calculating integrated flow for " << minpt << " < p_T > " << maxpt << " for particle " << number << endl;

	//Define index j used in particleList[j]
	int j = partid[MHALF+number];
	double fac, pt;
  //       double intvn[8][2] = {0};
	int nphi = particleList[j].nphi;
	int npt = particleList[j].npt;
	int neta = particleList[j].ny;
	double intvneta[etasize][nharmonics][2] = {0};
	double intvny[etasize][nharmonics][2] = {0};
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
	      for(int i = 0;i<nharmonics;i++)
	      {
		intvneta[ieta][i][0] += cos(i*phi)*dNdeta*2*PI/nphi;
		intvneta[ieta][i][1] += sin(i*phi)*dNdeta*2*PI/nphi;
		intvny[ieta][i][0] += cos(i*phi)*dNdy*2*PI/nphi;
		intvny[ieta][i][1] += sin(i*phi)*dNdy*2*PI/nphi;
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
	  
	  for(int i = 1;i<nharmonics;i++) for(int k =0;k<2;k++) 
	  {
	    vn[i][k][ieta] = intvneta[ieta][i][k]/intvneta[ieta][0][0];
// 	    vn[1][i][k][ieta] = intvny[ieta][i][k]/intvny[ieta][0][0];
	  }
	  
	  
	}// eta loop

}

// calculates eta-integrated yeild dN/dpt and eta-averaged flow versus pt.
// format is vn[n][i(real=0 or imaginary part=1)][pt] .
// Yield (n=0) is dN/dpt. 
// void Freeze::eta_integrated_flow(InitData *DATA, int number, double mineta, double maxeta, double ***vn)
// Setting etamin=etamax evaluates flow at a fixed pseudorapidity, with yield dN/dpt/deta
void Freeze::eta_integrated_flow(InitData *DATA, int number, double mineta, double maxeta, double vn[nharmonics][2][ptsize])
{
// 	cout << "Calculating integrated flow for " << mineta << " < eta < " << maxeta << " for particle " << number << endl;

    
	//Define index j used in particleList[j]
	int j = partid[MHALF+number];
	double fac, pt;
  //       double intvn[8][2] = {0};
	int nphi = particleList[j].nphi;
	int npt = particleList[j].npt;
	int neta = particleList[j].ny;
	double intvn[ptsize][nharmonics][2] = {0};
	double m = particleList[j].mass;
	
	if(mineta < particleList[j].y[0]) 
	{
	  cerr << "Error: called out of range eta in eta_integrated_flow, " 
	  << mineta << " < minimum " << particleList[j].y[0] << endl;
	  exit(1);
	}
	if(maxeta > particleList[j].y[neta-1]) 
	{
	  cerr << "Error: called out of range eta in eta_integrated_flow, " 
	  << maxeta << " > maximum " << particleList[j].y[neta-1] << endl;
	  exit(1);
	}
	if (mineta > maxeta)
	{
	  cerr << "Error in eta_integrated_flow:  mineta must be less than or equal to maxeta\n";
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
	      
	      // Integrate over pseudorapidity using gsl
	      double dndpt[ptsize] = {0};
	      for(int ieta=0;ieta<neta;ieta++) 
	      {
		double eta = particleList[j].y[ieta];
		// jacobian -- dN/deta = jac*dN/dy
		double jac = pt*cosh(eta)/sqrt(m*m + pt*pt*cosh(eta)*cosh(eta));
		dndpt[ieta] = pt*jac*particleList[j].dNdydptdphi[ieta][ipt][iphi];
// 		cout << "spectra = " << dndpt[ipt] << endl;
	      }
	      gsl_interp_accel *etaacc = gsl_interp_accel_alloc ();
	      gsl_spline *etaspline = gsl_spline_alloc (gsl_interp_cspline, neta);
	      gsl_spline_init (etaspline, particleList[j].y ,dndpt , neta);
	      
	      double dNdp;
	      if (mineta!=maxeta) dNdp = gsl_spline_eval_integ(etaspline, mineta, maxeta, etaacc);
	      else dNdp = gsl_spline_eval(etaspline, maxeta, etaacc);
	      
	      double phi = iphi*2*PI/nphi;
	      for(int i = 0;i<nharmonics;i++)
	      {
		intvn[ipt][i][0] += cos(i*phi)*dNdp*2*PI/nphi;
		intvn[ipt][i][1] += sin(i*phi)*dNdp*2*PI/nphi;
	      }
	      gsl_spline_free (etaspline);
	      gsl_interp_accel_free (etaacc);
	    }// phi loop


	  for(int k =0;k<2;k++) 
	  {
	    vn[0][k][ipt] = intvn[ipt][0][k];
	  }
	  
	  for(int i = 1;i<nharmonics;i++) for(int k =0;k<2;k++) 
	  {
	    vn[i][k][ipt] = intvn[ipt][i][k]/intvn[ipt][0][0];
	  }
	  
	  
	}// eta loop

}


// calculates y-integrated flow versus rapidity.
// format is vn[n][i(real=0 or imaginary part=1)][pt] .
// Yield (n=0) is dN/dpt. 
// void Freeze::eta_integrated_flow(InitData *DATA, int number, double mineta, double maxeta, double ***vn)
void Freeze::y_integrated_flow(InitData *DATA, int number, double miny, double maxy, double vn[nharmonics][2][etasize])
{
// 	cout << "Calculating integrated flow for " << miney << " < y < " << maxy << " for particle " << number << endl;

    
	//Define index j used in particleList[j]
	int j = partid[MHALF+number];
	double fac, pt;
  //       double intvn[8][2] = {0};
	int nphi = particleList[j].nphi;
	int npt = particleList[j].npt;
	int neta = particleList[j].ny;
	double intvn[etasize][nharmonics][2] = {0};
	double m = particleList[j].mass;
	
	
	if(miny < Rap(particleList[j].y[0],particleList[j].pt[npt-1],m)) 
	{
	  cerr << "Error: called out of range y in y_integrated_flow, " 
	  << miny << " < minimum " << Rap(particleList[j].y[0],particleList[j].pt[npt-1],m) << endl;
	  exit(1);
	}
	if(maxy > Rap(particleList[j].y[neta-1],particleList[j].pt[npt-1],m)) 
	{
	  cerr << "Error: called out of range y in y_integrated_flow, " 
	  << maxy << " > maximum " << Rap(particleList[j].y[neta-1],particleList[j].pt[npt-1],m) << endl;
	  exit(1);
	}
	if (miny > maxy)
	{
	  cerr << "Error in y_integrated_flow:  miny must be less than or equal to maxy\n";
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
		dndpt[ieta] = pt*particleList[j].dNdydptdphi[ieta][ipt][iphi];
	      }
	      gsl_interp_accel *yacc = gsl_interp_accel_alloc ();
	      gsl_spline *yspline = gsl_spline_alloc (gsl_interp_cspline, npt);
	      gsl_spline_init (yspline, ylist ,dndpt , npt);
	      
	      
	      double dNdp;
	      if (miny!=maxy) dNdp = gsl_spline_eval_integ(yspline, miny, maxy, yacc);
	      else dNdp = gsl_spline_eval(yspline, maxy, yacc);
	      
	      double phi = iphi*2*PI/nphi;
	      for(int i = 0;i<nharmonics;i++)
	      {
		intvn[ipt][i][0] += cos(i*phi)*dNdp*2*PI/nphi;
		intvn[ipt][i][1] += sin(i*phi)*dNdp*2*PI/nphi;
	      }
	      gsl_spline_free (yspline);
	      gsl_interp_accel_free (yacc);
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
  
  double vneta[nharmonics][2][etasize] = {0};
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
    
    gsl_spline_free (cosspline);
    gsl_interp_accel_free (cosacc);
    gsl_spline_free (sinspline);
    gsl_interp_accel_free (sinacc);
  }
}

  
// calculates pt- and eta-integrated flow for a given range in pt and eta
// format is vn[n][i(real=0 or imaginary part=1)]
// this one has the eta integral nested inside the phi integral inside the pt integral
void Freeze::pt_and_eta_integrated_flow2(InitData *DATA, int number, double minpt, double maxpt, double mineta, double maxeta, double vn[nharmonics][2])
{
  int j = partid[MHALF+number];
  int npt = particleList[j].npt;
//   int neta = particleList[j].ny;
  
  
  if(minpt < particleList[j].pt[0]) 
  {
    cerr << "Error: called out of range pt in pt_and_eta_integrated_flow2, " 
    << minpt << " < minimum " << particleList[j].pt[0] << endl;
    exit(1);
  }
  if(maxpt > particleList[j].pt[npt-1]) 
  {
    cerr << "Error: called out of range pt in pt_and_eta_integrated_flow2, " 
    << maxpt << " > maximum " << particleList[j].pt[npt-1] << endl;
    exit(1);
  }
  if (minpt > maxpt)
  {
    cerr << "Error in pt_and_eta_integrated_flow2:  minpt must be less than or equal to maxpt\n";
    exit(1);
  }
  
  
  
  // do  eta-integral first
  double vnpt[nharmonics][2][etasize] = {0};
  eta_integrated_flow(DATA, number, mineta, maxeta, vnpt);
  double dndpt[etasize];
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
    
    gsl_spline_free (cosspline);
    gsl_interp_accel_free (cosacc);
    gsl_spline_free (sinspline);
    gsl_interp_accel_free (sinacc);
  }
}

// calculates pt- and y-integrated flow for a given range in pt and y
// format is vn[n][i(real=0 or imaginary part=1)]
//  the y integral is nested inside the phi integral inside the pt integral
void Freeze::pt_and_y_integrated_flow(InitData *DATA, int number, double minpt, double maxpt, double miny, double maxy, double vn[nharmonics][2])
{
  int j = partid[MHALF+number];
  int npt = particleList[j].npt;
  int neta = particleList[j].ny;
  
  
  if(minpt < particleList[j].pt[0]) 
  {
    cerr << "Error: called out of range pt in pt_and_y_integrated_flow, " 
    << minpt << " < minimum " << particleList[j].pt[0] << endl;
    exit(1);
  }
  if(maxpt > particleList[j].pt[npt-1]) 
  {
    cerr << "Error: called out of range pt in pt_and_y_integrated_flow, " 
    << maxpt << " > maximum " << particleList[j].pt[npt-1] << endl;
    exit(1);
  }
  if (minpt > maxpt)
  {
    cerr << "Error in pt_and_y_integrated_flow:  minpt must be less than or equal to maxpt\n";
    exit(1);
  }
  
  // do  eta-integral first
  double vnpt[nharmonics][2][ptsize] = {0};
  y_integrated_flow(DATA, number, miny, maxy, vnpt);
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
      if (n==0) weight = dndpt[ipt];
      else weight = 1;
      vncos[ipt] = vnpt[n][0][ipt]*weight;
      vnsin[ipt] = vnpt[n][1][ipt]*weight;
    }
    gsl_interp_accel *cosacc = gsl_interp_accel_alloc ();
    gsl_spline *cosspline = gsl_spline_alloc (gsl_interp_cspline, npt);
    gsl_spline_init (cosspline, particleList[j].y ,vncos , npt);
    
    gsl_interp_accel *sinacc = gsl_interp_accel_alloc ();
    gsl_spline *sinspline = gsl_spline_alloc (gsl_interp_cspline, npt);
    gsl_spline_init (sinspline, particleList[j].y ,vnsin , npt);
    
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
    
    gsl_spline_free (cosspline);
    gsl_interp_accel_free (cosacc);
    gsl_spline_free (sinspline);
    gsl_interp_accel_free (sinacc);
  }
}

// Return yield in specified range of phase space
// if mineta==maxeta, returns dN/eta.  If minpt==maxpt, dN/dpt.  If both, dN/dpt/deta.  Otherwise, total yield N
double Freeze::get_yield(InitData *DATA, int number, double minpt, double maxpt, double mineta, double maxeta)
{
  double vn[nharmonics][2];
  pt_and_eta_integrated_flow(DATA, number, minpt, maxpt, mineta, maxeta, vn);
  return vn[0][0];
}


// Return v_n for specified range of phase space
double Freeze::get_vn(InitData *DATA, int number, double minpt, double maxpt, double mineta, double maxeta, int n)
{
  if( n > nharmonics-1 ) 
  {
    cout << "harmonic (" << n << ") too large.  Must increase nharmonics in freeze_pseudo.cpp\n";
    exit(1);
  }
  double vn[nharmonics][2];
  pt_and_eta_integrated_flow(DATA, number, minpt, maxpt, mineta, maxeta, vn);
  return sqrt(vn[n][0]*vn[n][0] + vn[n][1]*vn[n][1]);
}

// Return Psi_n for specified range of phase space
double Freeze::get_psi_n(InitData *DATA, int number, double minpt, double maxpt, double mineta, double maxeta, int n)
{
  if( n > nharmonics-1 ) 
  {
    cout << "harmonic (" << n << ") too large.  Must increase nharmonics in freeze_pseudo.cpp\n";
    exit(1);
  }
  double vn[nharmonics][2];
  pt_and_eta_integrated_flow(DATA, number, minpt, maxpt, mineta, maxeta, vn);
  return atan2(vn[n][1],vn[n][0]);
}

// Return charged hadron yield in specified range of phase space
// if mineta==maxeta, returns dN/eta.  If minpt==maxpt, dN/dpt.  If both, dN/dpt/deta.  Otherwise, total yield N
double Freeze::get_Nch(InitData *DATA, double minpt, double maxpt, double mineta, double maxeta)
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
      N+= get_yield(DATA, number, minpt, maxpt, mineta, maxeta);
    }
    return N;
}

// Return charged hadron vn for specified range of phase space
double Freeze::get_vn_ch(InitData *DATA, double minpt, double maxpt, double mineta, double maxeta, int n)
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
      pt_and_eta_integrated_flow(DATA, number, minpt, maxpt, mineta, maxeta, vn);
      numr+= vn[0][0]*vn[n][0];
      numi+= vn[0][0]*vn[n][1];
      den+= vn[0][0];
    }
  return sqrt(numr*numr+numi*numi)/den;
}

// Return charged hadron Psi_n for specified range of phase space
double Freeze::get_psi_n_ch(InitData *DATA, double minpt, double maxpt, double mineta, double maxeta, int n)
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
      pt_and_eta_integrated_flow(DATA, number, minpt, maxpt, mineta, maxeta, vn);
      numr+= vn[0][0]*vn[n][0];
      numi+= vn[0][0]*vn[n][1];
    }
  return atan2(numi, numr);
}

// Return yield and v_n at eta=0 as a function of pT
void Freeze::OutputDifferentialFlowAtMidrapidity2(InitData *DATA, int number, int full) 
{
  
	//Define index j used in particleList[j]
	int j = partid[MHALF+number];
	double fac, pt;
	double intvn[nharmonics][2] = {0};
	double m = particleList[j].mass;
	int nphi = particleList[j].nphi;
	int npt = particleList[j].npt;

	double minpt = particleList[j].pt[0];
	double maxpt = particleList[j].pt[npt-1];
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
		pt=particleList[j].pt[ipt];
// 		cout << "pt = " << pt << endl;
		
		pt_and_eta_integrated_flow(DATA, number, pt, pt, eta, eta, vn);

		//Output result
		outfilevn << pt;
		for(int i = 0;i<nharmonics;i++) for(int k =0;k<2;k++) outfilevn << "\t" << vn[i][k];
		outfilevn << endl;
		
		pt_and_eta_integrated_flow2(DATA, number, pt, pt, eta, eta, vn);
		
		outfilevn2 << pt;
		for(int i = 0;i<nharmonics;i++) for(int k =0;k<2;k++) outfilevn2 << "\t" << vn[i][k];
		outfilevn2 << endl;
	}
	

	//Close file
 	outfilevn.close();
	outfilevn2.close();
// 	}

}


// Return yield and v_n at eta~0 as a function of pT
void Freeze::OutputDifferentialFlowAtMidrapidity(InitData *DATA, int number, int full) 
{

  
	//Define index j used in particleList[j]
	int j = partid[MHALF+number];
	double fac, pt;
	double intvn[nharmonics][2] = {0};
	double m = particleList[j].mass;
	int nphi = particleList[j].nphi;
	int npt = particleList[j].npt;
//         double phipbuff[nphi];
// 	double resbuff[npt+1][nphi]; 
/*	
	gsl_fft_real_wavetable * real;
	gsl_fft_halfcomplex_wavetable * hc;
	gsl_fft_real_workspace * work;*/
	
//       cout << "here?\n";
//       ReadSpectra(DATA);
//       for (int j=1; j<particleMax; j++ )
// 	{
// 	number = particleList[j].number;
	cout << "Calculating flow at midrapidity for particle " << number << endl;

// 	for (int iphi=0;iphi<nphi;iphi++) phipbuff[iphi] = iphi*2*PI/nphi;

	//Set output file name
	string fname;
	stringstream tmpStr;
	fname="./outputs/";
	if (full) {
		fname+="F";
	}
	fname+="vnpteta0-";
	tmpStr << number;
	fname+=tmpStr.str();
	fname+=".dat";	
	
		//Open output file for vn
	ofstream outfilevn;
	outfilevn.open(fname.c_str());

	//Set the format of the output
	//outfile.width (10);
// 	outfile.precision(6);
// 	outfile.setf(ios::scientific);
	outfilevn.precision(6);
	outfilevn.setf(ios::scientific);
	outfilevn << "#pt\tdN/ptdYdptdphi\tv1cos\tv1sin\tv2cos\tv2sin\tv3cos\tv3sin\tv4cos\tv4sin\tv5cos\tv5sin\tv6cos\tv6sin\tv7cos\tv7sin\n";

	int ieta = particleList[j].ny/2;
	double eta = particleList[j].y[ieta];

	//Loop over pT
// 	cout << "npt = " << npt << endl;
	for(int ipt=0;ipt<npt;ipt++) {
		pt=particleList[j].pt[ipt];
// 		cout << "pt = " << pt << endl;
		
		//jacobian to switch from dN/dY to dN/deta
		double jac = pt*cosh(eta)/sqrt(m*m + pt*pt*cosh(eta)*cosh(eta))/pt/cosh(eta);
// 		double jac = 1.;
		
		for(int i = 0;i<nharmonics;i++) for(int k =0;k<2;k++) intvn[i][k]=0;

// 		cout << "ipt = " << ipt << endl;
		//Integrate over phi using trapezoid rule at closest point to midrapidity
		for(int iphi=0;iphi<nphi;iphi++) {
// 			int ieta = particleList[j].ny/2;
		        double dN = jac*particleList[j].dNdydptdphi[ieta][ipt][iphi];
			if (iphi==0) cout << "pt, dndydpt = " << pt << ", " << 2*PI*dN << endl;
			fac = 1.;
			double phi = iphi*2*PI/nphi;
			for(int i = 0;i<nharmonics;i++)
			{
			  intvn[i][0] += cos(i*phi)*fac*dN;
			  intvn[i][1] += sin(i*phi)*fac*dN;
			}
		}


		//Output result
		outfilevn << pt << "\t" << intvn[0][0]*2*PI/nphi;
		for(int i = 1;i<nharmonics;i++) for(int k =0;k<2;k++) outfilevn << "\t" << intvn[i][k]/intvn[0][0];
		outfilevn << endl;
	}

	//Close file
// 	outfile.close();
	outfilevn.close();
// 	}

}


// Return yield and v_n integrated over pT, as a function of pseudorapidity
void Freeze::OutputIntegratedFlow(InitData *DATA, int number, int full) 
{
  
	//Define index j used in particleList[j]
	int j = partid[MHALF+number];
	double intvn[nharmonics][2] = {0};
	double m = particleList[j].mass;
	int neta = particleList[j].ny;
	int npt = particleList[j].npt;

	double minpt = particleList[j].pt[0];
	double maxpt = particleList[j].pt[npt-1];
// 	double eta = 0.0;
	
	cout << "Calculating integrated flow vs. midrapidity for particle " << number << endl;

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
		
		pt_and_eta_integrated_flow(DATA, number, minpt, maxpt, eta, eta, vn);

		//Output result
		outfilevn << eta;
		for(int i = 0;i<nharmonics;i++) for(int k =0;k<2;k++) outfilevn << "\t" << vn[i][k];
		outfilevn << endl;
		
		pt_and_eta_integrated_flow2(DATA, number, minpt, maxpt, eta, eta, vn);
		
		outfilevn2 << eta;
		for(int i = 0;i<nharmonics;i++) for(int k =0;k<2;k++) outfilevn2 << "\t" << vn[i][k];
		outfilevn2 << endl;
	}
	

	//Close file
 	outfilevn.close();
	outfilevn2.close();
// 	}

}



//Output yield and v_n integrated over pT>500MeV as a function of pseudorapidity using gsl interpolation and integration
void Freeze::OutputIntegratedFlowForCMS(InitData *DATA, int number, int full) 
{

  
      //Define index j used in particleList[j]
      int j = partid[MHALF+number];
      double fac, pt;
//       double intvn[8][2] = {0};
      int nphi = particleList[j].nphi;
      int npt = particleList[j].npt;
      double m = particleList[j].mass;
      double minpt=0.5;
      double maxpt = particleList[j].pt[npt-1];
      double maxeta = 2.4;

	cout << "Calculating integrated flow for |p_T| > " << minpt << " for particle " << number << endl;

// 	for (int iphi=0;iphi<nphi;iphi++) phipbuff[iphi] = iphi*2*PI/nphi;

	//Set output file name
	string fname;
	stringstream tmpStr;
	fname="./outputs/";
	if (full) {
		fname+="F";
	}
	fname+="vnetaCMS-";
	tmpStr << number;
	fname+=tmpStr.str();
	fname+=".dat";	
	
		//Open output file for vn
	ofstream outfilevn;
	outfilevn.open(fname.c_str());

	//Set the format of the output
	outfilevn.precision(6);
	outfilevn.setf(ios::scientific);
	outfilevn << "#eta\tdN/ptdetadptdphi\tv1cos\tv1sin\tv2cos\tv2sin\tv3cos\tv3sin\tv4cos\tv4sin\tv5cos\tv5sin\tv6cos\tv6sin\tv7cos\tv7sin\n";

	double intvn[etasize][nharmonics][2] = {0};
	double etalist[etasize];
	
	
	//loop over pseudorapidity
// 	cout << "ietamax = " << particleList[j].ny << endl;
	for(int ieta=0;ieta<particleList[j].ny;ieta++)
	{
// 	  for(int i = 0;i<8;i++) for(int k =0;k<2;k++) intvn[ieta][i][k]=0;
	  double eta = particleList[j].y[ieta];
	  etalist[ieta] = eta;
	  
	  
	    //Integrate over phi using trapezoid rule
	    for(int iphi=0;iphi<nphi;iphi++) 
	    {
	      
	      // Integrate over pt using gsl
	      double dndpt[ptsize] = {0};
	      for(int ipt=0;ipt<npt;ipt++) 
	      {
		pt = particleList[j].pt[ipt];
		double jac = pt*cosh(eta)/sqrt(m*m + pt*pt*cosh(eta)*cosh(eta));
		dndpt[ipt] = jac*particleList[j].dNdydptdphi[ieta][ipt][iphi];
// 		cout << "spectra = " << dndpt[ipt] << endl;
	      }
	      gsl_interp_accel *ptacc = gsl_interp_accel_alloc ();
	      gsl_spline *ptspline = gsl_spline_alloc (gsl_interp_cspline, npt);
	      gsl_spline_init (ptspline, particleList[j].pt ,dndpt , npt);
	      
	      
	      double dN = gsl_spline_eval_integ(ptspline, minpt, maxpt,ptacc);
// 	      cout << "dN = " << dN << endl;
	      
	      double phi = iphi*2*PI/nphi;
	      for(int i = 0;i<nharmonics;i++)
	      {
		intvn[ieta][i][0] += cos(i*phi)*dN*2*PI/nphi;
		intvn[ieta][i][1] += sin(i*phi)*dN*2*PI/nphi;
	      }
	      gsl_spline_free (ptspline);
	      gsl_interp_accel_free (ptacc);
	    }// phi loop
// 	  for(int i = 1;i<8;i++) for(int k =0;k<2;k++) intvn[ieta][i][k]/=intvn[ieta][0][0];


	  
	  //Output result
	  outfilevn << eta << "\t" << intvn[ieta][0][0];
	  for(int i = 1;i<nharmonics;i++) for(int k =0;k<2;k++) outfilevn << "\t" << intvn[ieta][i][k]/intvn[ieta][0][0];
	  outfilevn << endl;
	  
	  
	}// eta loop
	
	//Close file
// 	outfile.close();
	outfilevn.close();
	
	
	double N = 0;
	double neta = particleList[j].ny;
	//integrate over pseudorapidity |eta| < 2.4
	if (particleList[j].ymax >=2.4)
	for (int i = 0;i<nharmonics;i++)
	{
	  double vncos[etasize] = {0};
	  double vnsin[etasize] = {0};
	  for(int ieta=0;ieta<neta;ieta++)
	  {
	    vncos[ieta] = intvn[ieta][i][0];
// 	    cout << "vncos at eta " << etalist[ieta] << " = " << vncos[ieta] << endl;
	    vnsin[ieta] = intvn[ieta][i][1];
	  }
	  gsl_interp_accel *etacacc = gsl_interp_accel_alloc ();
	  gsl_interp_accel *etasacc = gsl_interp_accel_alloc ();
	  gsl_spline *etacspline = gsl_spline_alloc (gsl_interp_cspline, neta);
	  gsl_spline *etasspline = gsl_spline_alloc (gsl_interp_cspline, neta);
	  gsl_spline_init (etacspline, etalist ,vncos , neta);
	  gsl_spline_init (etasspline, etalist ,vnsin , neta);
	  
	  double vnc = gsl_spline_eval_integ(etacspline, -maxeta, maxeta, etacacc)/2/maxeta;
	  double vns = gsl_spline_eval_integ(etasspline, -maxeta, maxeta, etacacc)/2/maxeta;
	  
	  if (i==0) N = vnc;
	  if (i==0) cout << "pt and eta integrated dN/deta = " << N << endl;
	  else 
	    cout << "pt and eta integrated v_" << i 
	    << ", Psi_" << i 
	    << " = " << sqrt(vnc*vnc+vns*vns) 
	    << ", " << atan2(vns, vnc)/i << endl;
	  
	  gsl_spline_free (etacspline);
	  gsl_spline_free (etasspline);
	  gsl_interp_accel_free (etacacc);
	  gsl_interp_accel_free (etasacc);
	  
	  
	}

// 	}

}


//Output total yield integrated over pT>400MeV and |eta|>max_pseudorapidity
double Freeze::OutputYieldForCMS(InitData *DATA, int number, int full) 
{

    
	//Define index j used in particleList[j]
	int j = partid[MHALF+number];
	double fac, pt;
	int nphi = particleList[j].nphi;
	int npt = particleList[j].npt;
	int neta = particleList[j].ny;
	double m = particleList[j].mass;
//	double *dndpt = new double[npt](); // This initializes array to zero
//	double *dndpt;
//       	dndpt = new double[npt]; // Not initialized
	double dndpt[50];
	double minpt=0.4; // in GeV
	double etamax = particleList[j].ymax;

	cout << "Calculating yield for " << minpt << " < p_T < " << particleList[j].pt[npt-1] 
	      << " and |eta| < " << etamax << " for particle " << number << ", " 
	      << particleList[j].name << endl;

	for(int ipt=0;ipt<npt;ipt++) 
	{
		pt = particleList[j].pt[ipt];
		dndpt[ipt] = 0.;
		double norm = 0.;
		
		//Integrate over pseudorapidity using trapezoid rule
		for(int ieta=0;ieta<neta;ieta++)
		{
			double dndpteta = 0;
			double eta = particleList[j].y[ieta];
			double jac = pt*cosh(eta)/sqrt(m*m + pt*pt*cosh(eta)*cosh(eta));
// 			cout << "jac = " << jac << endl;
			//Integrate over phi using trapezoid rule
			for(int iphi=0;iphi<nphi;iphi++) 
			{
				double dNdY = particleList[j].dNdydptdphi[ieta][ipt][iphi];
				double dNdeta = dNdY*jac;
				dndpteta+=dNdeta;
			}
			
			if(ieta==0 || ieta==neta-1) fac = 0.5;
			else fac = 1.0;
			dndpt[ipt]+= fac*dndpteta;
// 			cout << "dndpt = " << dndpt[ipt] << endl;

			
		}// eta loop
		if(0!=neta) dndpt[ipt]*=pt*2*PI/nphi*2*etamax/neta;
		else dndpt[ipt]*=pt*2*PI/nphi*2;
// 		cout << particleList[j].pt[ipt] << ", " << dndpt[ipt] << endl;	
	}// pt loop
	
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
        gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, npt);
//         gsl_spline *spline = gsl_spline_alloc (gsl_interp_akima, npt); //other types of splines to test
// 	gsl_spline *spline = gsl_spline_alloc (gsl_interp_polynomial, npt);
	gsl_spline_init (spline, particleList[j].pt ,dndpt , npt);
	
	
// 	for(int ipt=0;ipt<npt;ipt++) //for testing
// 	{
// 		double pt = particleList[j].pt[ipt];
// 		cout << particleList[j].pt[ipt] << ", " << 
// 		gsl_spline_eval(spline, pt, acc) << endl;
// 	}
// 	
// 	cout << "minpt = " << minpt << endl;
// 	cout << "maxpt = " << particleList[j].pt[npt] << endl;
// 	cout << gsl_spline_eval(spline, 0.250464, acc) << endl;
	double N = 
	  gsl_spline_eval_integ(spline, minpt, particleList[j].pt[npt-1],acc);
// 	  gsl_spline_eval_integ(spline, 1.0, 2.0, acc);
// 	  gsl_spline_eval(spline, 0.250464, acc);

// 	cout << "N = " << N << endl;
	
	gsl_spline_free (spline);
        gsl_interp_accel_free (acc);

//	delete [] dndpt;
	  
	return N;

}


//Output pT and phi integrated charged hadron spectra as a function of eta
void Freeze::Output_charged_hadrons_eta_differential_spectra(InitData *DATA, int full, const int * hadron_list, int nb_hadrons) {
	//
	int nb_hadron = nb_hadrons;
	
	int j, number, nphi, npt, neta;
	double eta;
	double m;
	double fac, pt;
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
	fname <<"dNdy";
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
