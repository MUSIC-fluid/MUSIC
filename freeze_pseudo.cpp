#include "freeze.h"

// read in thermal spectra from file to then perform resonance decays with them
void Freeze::ReadSpectra_pseudo(InitData* DATA, int full)
{
  // read in thermal spectra from file:
  int number, iptmax, iphimax;
  double deltaeta, etamax, ptmin, ptmax;
  int pseudo_steps;
  int ip;
  int bytes_read;
  double  neta;
  fprintf(stderr,"reading spectra\n");
  // open particle information file:
  FILE *p_file;
  char* p_name = "particleInformation.dat";
  char* pf_name = "FparticleInformation.dat";
//   if(full) strcpy(p_name, "FparticleInformation.dat");
  if(full) p_file = fopen(pf_name, "r");
  else p_file = fopen(p_name, "r");
  checkForReadError(p_file,p_name);
  int count;
  count = 0;
  cout << "NumberOfParticlesToInclude " << DATA->NumberOfParticlesToInclude << endl;
  // read particle information:
  while( fscanf(p_file,"%d %lf %d %lf %lf %d %d ",&number, &etamax, &pseudo_steps, &ptmin, &ptmax, &iptmax, &iphimax) == 7)
    {
      count ++;
      if (count>DATA->NumberOfParticlesToInclude) break;
// 	  fprintf(stderr,"%d %e %d %e %e %d %d \n", number, etamax, pseudo_steps, ptmin, ptmax, iptmax, iphimax);
      ip = partid[MHALF+number];
      particleList[ip].ny = pseudo_steps;
      particleList[ip].npt = iptmax;
      particleList[ip].nphi = iphimax;
      particleList[ip].phimin = 0;
      particleList[ip].phimax = 2*PI;
      particleList[ip].slope = 1;
      particleList[ip].ymax = etamax;
      deltaeta = 0.;
      if(pseudo_steps>1) deltaeta = 2*etamax/pseudo_steps;
      particleList[ip].deltaY = deltaeta;
// 	  cout << "ptmin = " << ptmin << endl;
      for (int ipt=0; ipt<=iptmax; ipt++ )
	{
// 	      particleList[ip].pt[ipt] =  ptmin + (ptmax - ptmin)*(exp(ipt)-1)/(exp(iptmax)-1); // log distributed values
	  particleList[ip].pt[ipt] =  ptmin + (ptmax - ptmin)*pow(ipt,2)/pow(iptmax,2); // power law
	}
      for (int ieta=0; ieta<=pseudo_steps; ieta++ )
	{
	  particleList[ip].y[ieta] =  ieta*deltaeta-etamax; // store pseudorapidity here
	  //	  if (ip==1) cout << "read particleList[ip].y[" << ieta << "] = " <<  particleList[ip].y[ieta] << endl;
	}
      phiArray = util->vector_malloc(iphimax);
      for(int iphi=0; iphi<iphimax; iphi++) phiArray[iphi] = iphi*2*PI/iphimax;
    }
  
  particleMax = ip;

  fclose(p_file);

  FILE *s_file;
  char* s_name = "yptphiSpectra.dat";
  char* sf_name = "FyptphiSpectra.dat";
  if(full) s_file = fopen(sf_name, "r");
  else s_file = fopen(s_name, "r");
  checkForReadError(s_file,s_name);
  
  cout << "ietamax=" << pseudo_steps+1 << endl;
  cout << "cells=" << (pseudo_steps+1)*(iptmax+1)*iphimax << endl;
  cout << "particleMax = " << particleMax << endl;
  for ( ip=1; ip<=particleMax; ip++ )
    {
      //cout << ip << endl;
      fprintf(stderr,"reading particle %d: %d %s\n", ip, particleList[ip].number, particleList[ip].name);
      for (int ieta=0; ieta<=pseudo_steps; ieta++)
	{
	  for (int ipt=0; ipt<=iptmax; ipt++)
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
  int iptmax = DATA->pt_steps; // Number of points is iptmax + 1 (ipt goes from 0 to iptmax)
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
 
  particleList[j].ny = ietamax*size;
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
  char* d_name = "particleInformation.dat";
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
  if (rank==0) fprintf(d_file,"%d %e %d %e %e %d %d \n", number,  etamax, DATA->pseudo_steps, ptmin, ptmax, iptmax, iphimax);
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

      for (int ipt=0; ipt<=iptmax; ipt++)
	{
// 	  pt = deltapt/2. + ipt*deltapt;
// 	  pt = ipt*deltapt;
// 	  pt = ptmin + (ptmax - ptmin)*(exp(ipt)-1)/(exp(iptmax)-1);
// 	      pt =  ptmin + (ptmax - ptmin)*(exp(ipt)-1)/(exp(iptmax)-1); // log distributed values
	      pt =  ptmin + (ptmax - ptmin)*pow(ipt,2)/pow(iptmax,2); // power law
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
	      double sum = summation(px, py, y, m, d, b, mu, DATA);
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
  char* d_name = "FparticleInformation.dat";
  d_file = fopen(d_name, "a");
  
  fprintf(d_file,"%d %e %d %e %e %d %d \n", number,  DATA->max_pseudorapidity, DATA->pseudo_steps, DATA->min_pt, DATA->max_pt, DATA->pt_steps, DATA->phi_steps);

//   fprintf(d_file,"%d %e %e %e %e %e %d %d %d \n", number, DATA->, ymax, slope, phimin, phimax, iymax, iptmax, iphimax);
   fclose(d_file);

  
  FILE *s_file;
      
//   s_name += "yptphiSpectra.dat";
  char* s_name = "FyptphiSpectra.dat";
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
	fprintf(stderr,"mr=%f\n",mr);
	fprintf(stderr,"m1=%f\n",m1);
	fprintf(stderr,"m2=%f\n",m2);
		
	for (n = 0; n <= neta; n++)
	  {
	    eta = particleList[pn].y[n];// for pseudofreeze, this variable stores pseudorapidity
	    for (l = 0; l <= npt; l++)
	      {
// 		if(pseudofreeze)  eta = Rap(y,particleList[pn].pt[l],m1);
		for (i = 0; i < nphi; i++)
		  {
		    double phi = i*deltaphi;
		    double spectrum = Edndp3_2bodyN (eta, particleList[pn].pt[l], phi, m1, m2, mr, particleList[pnR].number);
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
			  fprintf(stderr," %d %f %e %e %e %e\n", n, eta, particleList[pn].pt[l], decay[j].branch *
				  spectrum,
				  particleList[pn].dNdydptdphi[n][l][i]-decay[j].branch *
				  spectrum,particleList[pn].resCont[n][l][i]); 
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
	
	for (n = 0; n <= neta; n++)
	  {
	    eta = particleList[pn].y[n];// for pseudofreeze, this variable stores pseudorapidity
	    for (l = 0; l <= npt; l++)
	      {
// 		if(pseudofreeze)  y = Rap(particleList[pn].y[n],particleList[pn].pt[l],m1);
		for (i = 0; i < nphi; i++)
		  {
		    double phi = i*deltaphi;
		    double spectrum = Edndp3_3bodyN(eta, particleList[pn].pt[l], phi,
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
		       
			fprintf(stderr,"%d %f %e %e %e %e\n",n,eta, particleList[pn].pt[l], decay[j].branch *
				spectrum,
				particleList[pn].dNdydptdphi[n][l][i]-decay[j].branch *
				spectrum,particleList[pn].resCont[n][l][i]); 
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
	
	
	for (n = 0; n <= neta; n++)
	  {
	    eta = particleList[pn].y[n];// for pseudofreeze, this variable stores pseudorapidity
	    for (i = 0; i < nphi; i++)
	      {
		double phi = i*deltaphi;
		for (l = 0; l <= npt; l++)
		  {
// 		    if(pseudofreeze)  y = Rap(particleList[pn].y[n],particleList[pn].pt[l],m1);
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
  ReadParticleData(DATA, eos); // read in data for Cooper-Frye
  int i, b, number;
  if (mode == 3) // compute thermal spectra
    {
      int ret;
      if (rank == 0) ret = system("rm yptphiSpectra.dat yptphiSpectra?.dat yptphiSpectra??.dat particleInformation.dat 2> /dev/null");
      ReadFreezeOutSurface(DATA); // read freeze out surface (has to be done after the evolution of course)
      if (particleSpectrumNumber==0) // do all particles
	{
	  fprintf(stderr,"Doing all particles on this processor. May take a while ... \n");
	  for ( i=1; i<particleMax; i++ )
	    {
	      number = particleList[i].number;
	      b = particleList[i].baryon;
// 	      	  cout << "number2 = " << number << endl;
	      ComputeParticleSpectrum_pseudo(DATA, number, b, size, rank);
// 	      	  cout << "number3 = " << number << endl;
  
	      // Wait until all processors are finished and concatonate results
	      // (Don't want any processors to start writing the next particle to 
	      // file until the concatonation is done)
	      MPI_Barrier(MPI_COMM_WORLD);
	      
	      if(rank==0) ret = system("cat yptphiSpectra?.dat yptphiSpectra??.dat >> yptphiSpectra.dat 2> /dev/null");
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
	      
	      ret = system("cat yptphiSpectra?.dat yptphiSpectra??.dat >> yptphiSpectra.dat 2> /dev/null");
/*	      
	      ReadSingleSpectrum(DATA);
	      cout << "output results..." << endl;
	      OutputFullParticleSpectrum2(DATA, number, 4., b, 0);*/
	    }
	}
    }
  else if (mode==4) //  do resonance decays
    {
      ReadSpectra_pseudo(DATA, 0);
      int bound = 211; //number of lightest particle to calculate. 
	  cout << "particleaMax = " << particleMax << endl;
	  fprintf(stderr,"doing all from %i: %s to %i: %s.\n",particleMax,particleList[particleMax].name,
		  partid[MHALF+bound],particleList[partid[MHALF+bound]].name);
	  cal_reso_decays_pseudo(particleMax,decayMax,bound,mode, DATA->pseudofreeze);
      if (rank == 0) int ret = system("rm FyptphiSpectra.dat FparticleInformation.dat 2> /dev/null");
	  for ( i=1; i<particleMax; i++ )
	    {
	      number = particleList[i].number;
	      b = particleList[i].baryon;
	      OutputFullParticleSpectrum_pseudo(DATA, number, b, 1);
// 	      ComputeAveragePT(number,4.);
	    }
	  if(particleMax>=20)
	    {
// 	      ReadFullSpectra(DATA);
// 	      ComputeChargedHadrons(DATA,4.);
	    }
    }
  else if (mode==13) // take tabulated spectra and compute various observables and integrated quantities
    {
      ReadSpectra_pseudo(DATA, 0);
//       for ( i=1; i<particleMax; i++ )
      

      //Set output file name for total multiplicity
      string fname;
      stringstream tmpStr;
      fname="./outputs/NCMS.dat";	
      
      ofstream outfile;

      outfile.open(fname.c_str(),ios::trunc);

      //Set the format of the output
//       outfile.precision(4);
//       outfile.setf(ios::scientific);
      
      
      int chargedhd[5] = {1,3,4,5,17};
      for ( int k=0; k<5; k++ )
	{
	  i = chargedhd[k];
	  number = particleList[i].number;
	  OutputDifferentialFlowAtMidrapidity(DATA, number,0);
// 	  OutputIntegratedFlowForCMS(DATA, number,0);
	  
	  double N = OutputYieldForCMS(DATA, number,0);
	  
	  outfile << number << "\t" << N << endl;
	}

	outfile.close();
    }
  else if (mode==14) // take tabulated post-decay spectra and compute various observables and integrated quantities
    {
      ReadSpectra_pseudo(DATA, 1);
//       for ( i=1; i<particleMax; i++ )
      int chargedhd[5] = {1,3,4,5,17};
      

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
      for ( int k=0; k<5; k++ )
	{
	  i = chargedhd[k];
	  number = particleList[i].number;
 	  OutputDifferentialFlowAtMidrapidity(DATA, number,1);
//  	  OutputIntegratedFlowForCMS(DATA, number,1);
 	  N += OutputYieldForCMS(DATA, number,1);
	}
      cout << "Nch = " << N << endl;
      outfile << N << endl;
      outfile.close();
    }
    else 
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

//Output yield and v_n at eta~0 as a function of pT
void Freeze::OutputDifferentialFlowAtMidrapidity(InitData *DATA, int number, int full) 
{

  
	//Define index j used in particleList[j]
	int j = partid[MHALF+number];
	double fac, pt;
	double intvn[8][2] = {0};
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
	for(int ipt=0;ipt<=npt;ipt++) {
		pt=particleList[j].pt[ipt];
// 		cout << "pt = " << pt << endl;
		
		//jacobian to switch from dN/dY to dN/deta
		double jac = sqrt(m*m + pt*pt*cosh(eta)*cosh(eta))/pt/cosh(eta);

		for(int i = 0;i<8;i++) for(int k =0;k<2;k++) intvn[i][k]=0;

// 		cout << "ipt = " << ipt << endl;
		//Integrate over phi using trapezoid rule at closest point to midrapidity
		for(int iphi=0;iphi<nphi;iphi++) {
// 			int ieta = particleList[j].ny/2;
		        double dN = jac*particleList[j].dNdydptdphi[ieta][ipt][iphi];
			fac = 1.;
			double phi = iphi*2*PI/nphi;
			for(int i = 0;i<8;i++)
			{
			  intvn[i][0] += cos(i*phi)*fac*dN;
			  intvn[i][1] += sin(i*phi)*fac*dN;
			}
		}


		//Output result
		outfilevn << pt << "\t" << intvn[0][0]/2/PI;
		for(int i = 1;i<8;i++) for(int k =0;k<2;k++) outfilevn << "\t" << intvn[i][k]/intvn[0][0];
		outfilevn << endl;
	}

	//Close file
// 	outfile.close();
	outfilevn.close();
// 	}

}

//Output yield and v_n integrated over pT>500MeV (or the closest point in pT that has been calculated) as a function of pseudorapidity
void Freeze::OutputIntegratedFlowForCMS(InitData *DATA, int number, int full) 
{

  
      //Define index j used in particleList[j]
      int j = partid[MHALF+number];
      double fac, pt;
      double intvn[8][2] = {0};
      int nphi = particleList[j].nphi;
      int npt = particleList[j].npt;
      int minipt=0;
      double minpt=0.;

      // Simple trapezoid rule needs the value at the boundary.  
      // Find nearest lower bound to desired value and integrate to highest pT available
      for(int ipt=0;ipt<=npt;ipt++) 
      {
	pt=particleList[j].pt[ipt];
	if(pt<=0.5)
	{
	  minipt = ipt;
	  minpt = pt;
	}
	if(pt>0.5) break;
      }
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
	outfilevn << "#eta\tdN/ptdYdptdphi\tv1cos\tv1sin\tv2cos\tv2sin\tv3cos\tv3sin\tv4cos\tv4sin\tv5cos\tv5sin\tv6cos\tv6sin\tv7cos\tv7sin\n";


	//loop over pseudorapidity
	for(int ieta=0;ieta<=particleList[j].ny;ieta++)
	{
	  for(int i = 0;i<8;i++) for(int k =0;k<2;k++) intvn[i][k]=0;
	  double eta = particleList[j].y[ieta];
// 	  cout << "eta = " << eta << endl;
	  //Loop over pT
  // 	cout << "npt = " << npt << endl;
	  for(int ipt=minipt;ipt<=npt;ipt++) 
	  {
		  pt=particleList[j].pt[ipt];
		  double intvnphi[8][2] = {0};
		  //Integrate over phi using trapezoid rule
		  for(int iphi=0;iphi<nphi;iphi++) {
  // 			int ieta = particleList[j].ny/2;
			  double dN = particleList[j].dNdydptdphi[ieta][ipt][iphi];
			  fac = 1.;
			  double phi = iphi*2*PI/nphi;
			  for(int i = 0;i<8;i++)
			  {
			    intvnphi[i][0] += cos(i*phi)*fac*dN;
			    intvnphi[i][1] += sin(i*phi)*fac*dN;
			  }
		  }
		  for(int i = 1;i<8;i++) for(int k =0;k<2;k++) intvnphi[i][k]/=intvnphi[0][0];
		  if(ipt==minipt || ipt == npt) fac = 0.5;
		  else fac = 1.0;
		  for(int i = 0;i<8;i++) for(int k =0;k<2;k++) intvn[i][k]+=fac*intvnphi[i][k];
	  }// pt loop
	  
	  //Output result
	  outfilevn << eta << "\t" << intvn[0][0]/2/PI;
	  for(int i = 1;i<8;i++) for(int k =0;k<2;k++) outfilevn << "\t" << intvn[i][k]/intvn[0][0];
	  outfilevn << endl;
	}// eta loop

	//Close file
// 	outfile.close();
	outfilevn.close();
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

	cout << "Calculating yield for " << minpt << " < p_T < " << particleList[j].pt[npt] 
	      << " and |eta| < " << etamax << " for particle " << number << ", " 
	      << particleList[j].name << endl;

	for(int ipt=0;ipt<=npt;ipt++) 
	{
		pt = particleList[j].pt[ipt];
		dndpt[ipt] = 0.;
		double norm = 0.;
		
		//Integrate over pseudorapidity using trapezoid rule
		for(int ieta=0;ieta<=neta;ieta++)
		{
			double dndpteta = 0;
			double eta = particleList[j].y[ieta];
			double jac = sqrt(m*m + pt*pt*cosh(eta)*cosh(eta))/pt/cosh(eta);
			//Integrate over phi using trapezoid rule
			for(int iphi=0;iphi<nphi;iphi++) 
			{
				double dNdY = particleList[j].dNdydptdphi[ieta][ipt][iphi];
				double dNdeta = dNdY*jac;
				dndpteta+=dNdeta;
			}
			
			if(ieta==0 || ieta==neta) fac = 0.5;
			else fac = 1.0;
			dndpt[ipt]+= fac*dndpteta;

			
		}// eta loop
		dndpt[ipt]*=2*PI/nphi*2*etamax/neta;
// 		cout << particleList[j].pt[ipt] << ", " << dndpt[ipt] << endl;	
	}// pt loop
	
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
        gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, npt+1);
//         gsl_spline *spline = gsl_spline_alloc (gsl_interp_akima, npt+1); //other types of splines to test
// 	gsl_spline *spline = gsl_spline_alloc (gsl_interp_polynomial, npt+1);
	gsl_spline_init (spline, particleList[j].pt ,dndpt , npt+1);
	
	
// 	for(int ipt=0;ipt<=npt;ipt++) //for testing
// 	{
// 		double pt = particleList[j].pt[ipt];
// 		cout << particleList[j].pt[ipt] << ", " << 
// 		gsl_spline_eval(spline, pt, acc) << endl;
// 	}
// 	
// 	cout << "minpt = " << minpt << endl;
// 	cout << "maxpt = " << particleList[j].pt[npt] << endl;
	double N = 
	  gsl_spline_eval_integ(spline, minpt, particleList[j].pt[npt],acc);
// 	  gsl_spline_eval_integ(spline, 1.0, 2.0, acc);
// 	  gsl_spline_eval(spline, 0.250464, acc);

// 	cout << "N = " << N << endl;
	
	gsl_spline_free (spline);
        gsl_interp_accel_free (acc);

//	delete [] dndpt;
	  
	return N;

}
