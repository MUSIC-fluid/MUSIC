#include "util.h"
#include "eos.h"
using namespace std;

#define cs2 (1.0/3.0)

EOS::EOS()
{
  util = new Util;
}

// destructor
EOS::~EOS()
{
  delete util;
}

void EOS::checkForReadError(FILE *file, char* name)
{
  if(!(file))
    {
      fprintf(stderr, "file %s not found.\n", name);
      fprintf(stderr, "Exiting...\n");
      exit(0);
    }
}

void EOS::init_eos()
{
  // read the azhydro pressure, temperature, and 
  // baryon chemical potential from file
  whichEOS = 1;
  fprintf(stderr,"reading EOS... \n");
  int bytes_read;
  int i, j;
  FILE *eos_p1, *eos_p2;
  FILE *eos_T1, *eos_T2;
  FILE *eos_mu1, *eos_mu2;
  const char* EOSPATH = "HYDROPROGRAMPATH";
  char* envPath = getenv(EOSPATH);
  char* eos_p1_name;
  char* eos_p2_name;
  char* eos_T1_name;
  char* eos_T2_name;
  char* eos_mu1_name;
  char* eos_mu2_name;
      
  if (envPath != 0 && *envPath != '\0') 
    {
      fprintf(stderr,"from path %s/EOS \n", envPath);
      eos_p1_name = util->char_malloc(100);
      strcat(eos_p1_name,envPath);
      strcat(eos_p1_name,"/EOS/aa1_p.dat");

      eos_p2_name = util->char_malloc(100);
      strcat(eos_p2_name,envPath);
      strcat(eos_p2_name,"/EOS/aa2_p.dat");

      eos_T1_name = util->char_malloc(100);
      strcat(eos_T1_name,envPath);
      strcat(eos_T1_name,"/EOS/aa1_t.dat");

      eos_T2_name = util->char_malloc(100);
      strcat(eos_T2_name,envPath);
      strcat(eos_T2_name,"/EOS/aa2_t.dat");

      eos_mu1_name = util->char_malloc(100);
      strcat(eos_mu1_name,envPath);
      strcat(eos_mu1_name,"/EOS/aa1_mb.dat");
 
      eos_mu2_name = util->char_malloc(100);
      strcat(eos_mu2_name,envPath);
      strcat(eos_mu2_name,"/EOS/aa2_mb.dat");
    }
  else
    {
      eos_p1_name = util->char_malloc(100);
      strcat(eos_p1_name,".");
      strcat(eos_p1_name,"/EOS/aa1_p.dat");

      eos_p2_name = util->char_malloc(100);
      strcat(eos_p2_name,".");
      strcat(eos_p2_name,"/EOS/aa2_p.dat");

      eos_T1_name = util->char_malloc(100);
      strcat(eos_T1_name,".");
      strcat(eos_T1_name,"/EOS/aa1_t.dat");

      eos_T2_name = util->char_malloc(100);
      strcat(eos_T2_name,".");
      strcat(eos_T2_name,"/EOS/aa2_t.dat");

      eos_mu1_name = util->char_malloc(100);
      strcat(eos_mu1_name,".");
      strcat(eos_mu1_name,"/EOS/aa1_mb.dat");
 
      eos_mu2_name = util->char_malloc(100);
      strcat(eos_mu2_name,".");
      strcat(eos_mu2_name,"/EOS/aa2_mb.dat");
    }
  
  eos_p1 = fopen(eos_p1_name, "r");
  eos_p2 = fopen(eos_p2_name, "r");
  eos_T1 = fopen(eos_T1_name, "r");
  eos_T2 = fopen(eos_T2_name, "r");
  eos_mu1 = fopen(eos_mu1_name, "r");
  eos_mu2 = fopen(eos_mu2_name, "r");
  
  checkForReadError(eos_p1,eos_p1_name);
  checkForReadError(eos_p2,eos_p2_name);
  checkForReadError(eos_T1,eos_T1_name);
  checkForReadError(eos_T2,eos_T2_name);
  checkForReadError(eos_mu1,eos_mu1_name);
  checkForReadError(eos_mu2,eos_mu2_name);
 
  //read the first two lines:
  // first value of rhob, first value of epsilon
  // deltaRhob, deltaEpsilon, number of rhob steps, number of epsilon steps
  bytes_read=fscanf(eos_p1,"%lf %lf",&BNP1,&EPP1);
  bytes_read=fscanf(eos_p1,"%lf %lf %d %d",&deltaBNP1,&deltaEPP1,&NBNP1,&NEPP1);
  bytes_read=fscanf(eos_p2,"%lf %lf",&BNP2,&EPP2);
  bytes_read=fscanf(eos_p2,"%lf %lf %d %d",&deltaBNP2,&deltaEPP2,&NBNP2,&NEPP2);
  bytes_read=fscanf(eos_T1,"%lf %lf",&BNP1,&EPP1);
  bytes_read=fscanf(eos_T1,"%lf %lf %d %d",&deltaBNP1,&deltaEPP1,&NBNP1,&NEPP1);
  bytes_read=fscanf(eos_T2,"%lf %lf",&BNP2,&EPP2);
  bytes_read=fscanf(eos_T2,"%lf %lf %d %d",&deltaBNP2,&deltaEPP2,&NBNP2,&NEPP2);
  bytes_read=fscanf(eos_mu1,"%lf %lf",&BNP1,&EPP1);
  bytes_read=fscanf(eos_mu1,"%lf %lf %d %d",&deltaBNP1,&deltaEPP1,&NBNP1,&NEPP1);
  bytes_read=fscanf(eos_mu2,"%lf %lf",&BNP2,&EPP2);
  bytes_read=fscanf(eos_mu2,"%lf %lf %d %d",&deltaBNP2,&deltaEPP2,&NBNP2,&NEPP2);
 
  // allocate memory for pressure arrays
  pressure1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  pressure2=util->mtx_malloc(NBNP2+1,NEPP2+1);

  // allocate memory for temperature arrays
  // we assume here that all files have the same structure (i.e., same EPP1, deltaEPP1, etc....)
  temperature1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  temperature2=util->mtx_malloc(NBNP2+1,NEPP2+1);

  // allocate memory for baryon chemical potential arrays
  mu1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  mu2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  
  // read pressure, temperature and chemical potential values
  for(j=0;j<=NEPP1;j++)
    for(i=0;i<=NBNP1;i++)
      {
	bytes_read=fscanf(eos_p1,"%lf",&pressure1[i][j]);
	bytes_read=fscanf(eos_T1,"%lf",&temperature1[i][j]);
	bytes_read=fscanf(eos_mu1,"%lf",&mu1[i][j]);
      }

  for(j=0;j<=NEPP2;j++)
    for(i=0;i<=NBNP2;i++)
      {
	bytes_read=fscanf(eos_p2,"%lf",&pressure2[i][j]);
	bytes_read=fscanf(eos_T2,"%lf",&temperature2[i][j]);
	bytes_read=fscanf(eos_mu2,"%lf",&mu2[i][j]);
      }

  fclose(eos_p1);
  fclose(eos_p2);
  fclose(eos_T1);
  fclose(eos_T2);
  fclose(eos_mu1);
  fclose(eos_mu2);
  fprintf(stderr,"done reading\n");
}

void EOS::init_eos2()
{
  // read the lattice EOS pressure, temperature, and 
  // baryon chemical potential from file
  fprintf(stderr,"reading EOS... \n");
  whichEOS = 2; 
   int bytes_read;
   int i, j;
  FILE *eos_d1, *eos_d2, *eos_d3, *eos_d4;
  FILE *eos_T1, *eos_T2, *eos_T3, *eos_T4;
  const char* EOSPATH = "HYDROPROGRAMPATH";
  char* envPath = getenv(EOSPATH);
  char* eos_d1_name;
  char* eos_d2_name;
  char* eos_d3_name;
  char* eos_d4_name;
  char* eos_T1_name;
  char* eos_T2_name;
  char* eos_T3_name;
  char* eos_T4_name;
  double eps, baryonDensity; //dummies for now
      
  if (envPath != 0 && *envPath != '\0') // if path is set in the environment
    {
      fprintf(stderr,"from path %s/EOS/s95p-v1 \n", envPath);
      eos_d1_name = util->char_malloc(100);
      strcat(eos_d1_name,envPath);
      strcat(eos_d1_name,"/EOS/s95p-v1/s95p-v1_dens1.dat");
      eos_d2_name = util->char_malloc(100);
      strcat(eos_d2_name,envPath);
      strcat(eos_d2_name,"/EOS/s95p-v1/s95p-v1_dens2.dat");
      eos_d3_name = util->char_malloc(100);
      strcat(eos_d3_name,envPath);
      strcat(eos_d3_name,"/EOS/s95p-v1/s95p-v1_dens3.dat");
      eos_d4_name = util->char_malloc(100);
      strcat(eos_d4_name,envPath);
      strcat(eos_d4_name,"/EOS/s95p-v1/s95p-v1_dens4.dat");
      eos_T1_name = util->char_malloc(100);
      strcat(eos_T1_name,envPath);
      strcat(eos_T1_name,"/EOS/s95p-v1/s95p-v1_par1.dat");
      eos_T2_name = util->char_malloc(100);
      strcat(eos_T2_name,envPath);
      strcat(eos_T2_name,"/EOS/s95p-v1/s95p-v1_par2.dat");
      eos_T3_name = util->char_malloc(100);
      strcat(eos_T3_name,envPath);
      strcat(eos_T3_name,"/EOS/s95p-v1/s95p-v1_par3.dat");
      eos_T4_name = util->char_malloc(100);
      strcat(eos_T4_name,envPath);
      strcat(eos_T4_name,"/EOS/s95p-v1/s95p-v1_par4.dat");
    }
  else //if path is not set in the environment use current folder and then /EOS/s95p-v1 subfolder
    {
      eos_d1_name = util->char_malloc(100);
      strcat(eos_d1_name,".");
      strcat(eos_d1_name,"/EOS/s95p-v1/s95p-v1_dens1.dat");
      eos_d2_name = util->char_malloc(100);
      strcat(eos_d2_name,".");
      strcat(eos_d2_name,"/EOS/s95p-v1/s95p-v1_dens2.dat");
      eos_d3_name = util->char_malloc(100);
      strcat(eos_d3_name,".");
      strcat(eos_d3_name,"/EOS/s95p-v1/s95p-v1_dens3.dat");
      eos_d4_name = util->char_malloc(100);
      strcat(eos_d4_name,".");
      strcat(eos_d4_name,"/EOS/s95p-v1/s95p-v1_dens4.dat");
      eos_T1_name = util->char_malloc(100);
      strcat(eos_T1_name,".");
      strcat(eos_T1_name,"/EOS/s95p-v1/s95p-v1_par1.dat");
      eos_T2_name = util->char_malloc(100);
      strcat(eos_T2_name,".");
      strcat(eos_T2_name,"/EOS/s95p-v1/s95p-v1_par2.dat");
      eos_T3_name = util->char_malloc(100);
      strcat(eos_T3_name,".");
      strcat(eos_T3_name,"/EOS/s95p-v1/s95p-v1_par3.dat");
      eos_T4_name = util->char_malloc(100);
      strcat(eos_T4_name,".");
      strcat(eos_T4_name,"/EOS/s95p-v1/s95p-v1_par4.dat");
    }
  
  eos_d1 = fopen(eos_d1_name, "r");
  eos_d2 = fopen(eos_d2_name, "r");
  eos_d3 = fopen(eos_d3_name, "r");
  eos_d4 = fopen(eos_d4_name, "r");
  eos_T1 = fopen(eos_T1_name, "r");
  eos_T2 = fopen(eos_T2_name, "r");
  eos_T3 = fopen(eos_T3_name, "r");
  eos_T4 = fopen(eos_T4_name, "r");
  
  checkForReadError(eos_d1,eos_d1_name);
  checkForReadError(eos_d2,eos_d2_name);
  checkForReadError(eos_d3,eos_d3_name);
  checkForReadError(eos_d4,eos_d4_name);
  checkForReadError(eos_T1,eos_T1_name);
  checkForReadError(eos_T2,eos_T2_name);
  checkForReadError(eos_T3,eos_T3_name);
  checkForReadError(eos_T4,eos_T4_name);
 
  // read the first two lines with general info:
  // lowest value of epsilon
  // deltaEpsilon, number of epsilon steps (i.e. # of lines)
  bytes_read=fscanf(eos_T1,"%lf",&EPP1);
  bytes_read=fscanf(eos_T1,"%lf %d",&deltaEPP1,&NEPP1);
  bytes_read=fscanf(eos_T2,"%lf",&EPP2);
  bytes_read=fscanf(eos_T2,"%lf %d",&deltaEPP2,&NEPP2);
  bytes_read=fscanf(eos_T3,"%lf",&EPP3);
  bytes_read=fscanf(eos_T3,"%lf %d",&deltaEPP3,&NEPP3);
  bytes_read=fscanf(eos_T4,"%lf",&EPP4);
  bytes_read=fscanf(eos_T4,"%lf %d",&deltaEPP4,&NEPP4);
  bytes_read=fscanf(eos_d1,"%lf",&EPP1);
  bytes_read=fscanf(eos_d1,"%lf %d",&deltaEPP1,&NEPP1);
  bytes_read=fscanf(eos_d2,"%lf",&EPP2);
  bytes_read=fscanf(eos_d2,"%lf %d",&deltaEPP2,&NEPP2);
  bytes_read=fscanf(eos_d3,"%lf",&EPP3);
  bytes_read=fscanf(eos_d3,"%lf %d",&deltaEPP3,&NEPP3);
  bytes_read=fscanf(eos_d4,"%lf",&EPP4);
  bytes_read=fscanf(eos_d4,"%lf %d",&deltaEPP4,&NEPP4);
 
  // no rho_b dependence at the moment
  NBNP1=0; 
  NBNP2=0;
  NBNP3=0;
  NBNP4=0;

  // allocate memory for pressure arrays
  pressure1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  pressure2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  pressure3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  pressure4=util->mtx_malloc(NBNP4+1,NEPP4+1);

  // allocate memory for entropy density arrays
  entropyDensity1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  entropyDensity2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  entropyDensity3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  entropyDensity4=util->mtx_malloc(NBNP4+1,NEPP4+1);

  // allocate memory for QGP fraction arrays
  QGPfraction1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  QGPfraction2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  QGPfraction3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  QGPfraction4=util->mtx_malloc(NBNP4+1,NEPP4+1);

  // allocate memory for temperature arrays
  temperature1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  temperature2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  temperature3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  temperature4=util->mtx_malloc(NBNP4+1,NEPP4+1);

  // allocate memory for baryon chemical potential arrays
  // currently always zero
  /*   mu1=mtx_malloc(NBNP1+1,NEPP1+1); */
  /*   mu2=mtx_malloc(NBNP2+1,NEPP2+1); */
  /*   mu3=mtx_malloc(NBNP3+1,NEPP3+1); */
  /*   mu4=mtx_malloc(NBNP4+1,NEPP4+1); */
  
  // read pressure, temperature and chemical potential values
  // files have it backwards, so I start with maximum j and count down
  i=0;
  for(j=NEPP1-1; j>=0; j--)
    {
      bytes_read=fscanf(eos_d1,"%lf",&eps);
      bytes_read=fscanf(eos_d1,"%lf",&pressure1[i][j]);
      bytes_read=fscanf(eos_d1,"%lf",&entropyDensity1[i][j]);
      bytes_read=fscanf(eos_d1,"%lf",&baryonDensity);
      bytes_read=fscanf(eos_d1,"%lf",&QGPfraction1[i][j]);
      bytes_read=fscanf(eos_T1,"%lf",&temperature1[i][j]);
      bytes_read=fscanf(eos_T1,"%lf",&eps); //dummy
      bytes_read=fscanf(eos_T1,"%lf",&eps); //dummy
    }

  for(j=NEPP2-1; j>=0; j--)
    {
      bytes_read=fscanf(eos_d2,"%lf",&eps);
      bytes_read=fscanf(eos_d2,"%lf",&pressure2[i][j]);
      bytes_read=fscanf(eos_d2,"%lf",&entropyDensity2[i][j]);
      bytes_read=fscanf(eos_d2,"%lf",&baryonDensity);
      bytes_read=fscanf(eos_d2,"%lf",&QGPfraction2[i][j]);
      bytes_read=fscanf(eos_T2,"%lf",&temperature2[i][j]);
      bytes_read=fscanf(eos_T2,"%lf",&eps); //dummy
      bytes_read=fscanf(eos_T2,"%lf",&eps); //dummy
    }

  for(j=NEPP3-1; j>=0; j--)
    {
      bytes_read=fscanf(eos_d3,"%lf",&eps);
      bytes_read=fscanf(eos_d3,"%lf",&pressure3[i][j]);
      bytes_read=fscanf(eos_d3,"%lf",&entropyDensity3[i][j]);
      bytes_read=fscanf(eos_d3,"%lf",&baryonDensity);
      bytes_read=fscanf(eos_d3,"%lf",&QGPfraction3[i][j]);
      bytes_read=fscanf(eos_T3,"%lf",&temperature3[i][j]);
      bytes_read=fscanf(eos_T3,"%lf",&eps); //dummy
      bytes_read=fscanf(eos_T3,"%lf",&eps); //dummy
    }

  for(j=NEPP4-1; j>=0; j--)
    {
      bytes_read=fscanf(eos_d4,"%lf",&eps);
      bytes_read=fscanf(eos_d4,"%lf",&pressure4[i][j]);
      bytes_read=fscanf(eos_d4,"%lf",&entropyDensity4[i][j]);
      bytes_read=fscanf(eos_d4,"%lf",&baryonDensity);
      bytes_read=fscanf(eos_d4,"%lf",&QGPfraction4[i][j]);
      bytes_read=fscanf(eos_T4,"%lf",&temperature4[i][j]);
      bytes_read=fscanf(eos_T4,"%lf",&eps); //dummy
      bytes_read=fscanf(eos_T4,"%lf",&eps); //dummy
    }

//test if the reading worked:
/*   for(j=0; j<NEPP1; j++) */
/*     { */
/*       fprintf(stderr,"%lf %lf %lf %lf\n",pressure1[i][j],entropyDensity1[i][j],QGPfraction1[i][j],temperature1[i][j]); */
/*     } */
/*   for(j=0; j<NEPP2; j++) */
/*     { */
/*       fprintf(stderr,"%lf %lf %lf %lf\n",pressure2[i][j],entropyDensity2[i][j],QGPfraction2[i][j],temperature2[i][j]); */
/*     } */
/*   for(j=0; j<NEPP3; j++) */
/*     { */
/*       fprintf(stderr,"%lf %lf %lf %lf\n",pressure3[i][j],entropyDensity3[i][j],QGPfraction3[i][j],temperature3[i][j]); */
/*     } */
/*   for(j=0; j<NEPP4; j++) */
/*     { */
/*       fprintf(stderr,"%lf %lf %lf %lf\n",pressure4[i][j],entropyDensity4[i][j],QGPfraction4[i][j],temperature4[i][j]); */
/*     } */

  fclose(eos_d1);
  fclose(eos_d2);
  fclose(eos_d3);
  fclose(eos_d4);
  fclose(eos_T1);
  fclose(eos_T2);
  fclose(eos_T3);
  fclose(eos_T4);
}




double EOS::interpolate_pressure(double e, double rhob)
{
  double p, pa, pb, pa1, pa2, pb1, pb2;
  //use linear interpolation
  int ie1, ie2, inb1, inb2;
  double frace, fracnb;
  
  //if(rhob>0.00001) fprintf(stderr,"rhob=%lf\n", rhob);
  
  e*=hbarc; // in the files epsilon is in GeV/fm^3

  if(e>24.) return (e-4.*.3642)/3./hbarc;

  if(e<EPP2) //use first file for small epsilon values
    {
      if(e<EPP1) 
	{
	  ie1 = 0;
	  ie2 = 1;
	  inb1 = floor((rhob-BNP1)/deltaBNP1);
	  inb2 = floor((rhob-BNP1)/deltaBNP1+1);
	  frace = e/(EPP1);
	}
      else
	{
	  ie1 = floor((e-EPP1)/deltaEPP1);
	  ie2 = floor((e-EPP1)/deltaEPP1+1);
	  inb1 = floor((rhob-BNP1)/deltaBNP1);
	  inb2 = floor((rhob-BNP1)/deltaBNP1+1);
	  frace = (e-(ie1*deltaEPP1+EPP1))/deltaEPP1; 
	}
      
      if(ie1>NEPP1 || inb1>NBNP1 || ie1<0 || inb1<0 )
	{
	  fprintf(stderr,"ERROR in inperpolate_pressure. out of range.\n");
	  fprintf(stderr,"ie1=%d,NEPP1=%d; inb1=%d, NBNP1=%d, e=%f, rhob=%f \n", ie1, NEPP1, inb1, NBNP1,e,rhob);
	  exit(0);
	}
      if(ie2>NEPP1 || inb2>NBNP1 || ie2<0 || inb2<0 )
	{
	  fprintf(stderr,"ERROR in inperpolate_pressure. out of range.\n");
	  fprintf(stderr,"ie2=%d,NEPP1=%d; inb2=%d, NBNP1=%d, e=%f, rhob=%f \n", ie2, NEPP1, inb2, NBNP1,e,rhob);
	  exit(0);
	}

      pa1 = pressure1[inb1][ie1];
      pa2 = pressure1[inb2][ie1];

      fracnb = (rhob-(inb1*deltaBNP1+BNP1))/deltaBNP1; 
      
      pa = pa1*(1-fracnb) + pa2*fracnb;
      
      pb1 = pressure1[inb1][ie2];
      pb2 = pressure1[inb2][ie2];
      
      pb = pb1*(1-fracnb) + pb2*fracnb;
      
     
      if(e<EPP1) 
	{
	  p = pa*(frace);
	  //if (p<0) fprintf(stderr,"pa=%lf\n", pa);
	  //if (p<0) fprintf(stderr,"p=%lf\n", p);
	}
      else
	{
	  p = pa*(1-frace) + pb*frace;
	}
    }
  else // use second file for larger epsilon values
    {
      ie1 = floor((e-EPP2)/deltaEPP2);
      ie2 = floor((e-EPP2)/deltaEPP2+1);
      
      inb1 = floor((rhob-BNP2)/deltaBNP2);
      inb2 = floor((rhob-BNP2)/deltaBNP2+1);
      
      if(ie1>NEPP2 || inb1>NBNP2 || ie1<0 || inb1<0 )
	{
	  fprintf(stderr,"ERROR in inperpolate_pressure. out of range.\n");
	  fprintf(stderr,"ie1=%d,NEPP2=%d; inb1=%d, NBNP2=%d\n", ie1, NEPP2, inb1, NBNP2);
	  fprintf(stderr,"e=%f,EPP2=%f; maxe=%f, deltaEPP2=%f\n", e, EPP2, NEPP2*deltaEPP2+EPP2, deltaEPP2);
	  exit(0);
	}
      if(ie2>NEPP2 || inb2>NBNP2 || ie2<0 || inb2<0 )
	{
	  fprintf(stderr,"ERROR in inperpolate_pressure. out of range.\n");
	  fprintf(stderr,"ie2=%d,NEPP2=%d; inb2=%d, NBNP2=%d\n", ie2, NEPP2, inb2, NBNP2);
	  exit(0);
	}
/*       fprintf(stderr,"ie1=%d\n", ie1); */
/*       fprintf(stderr,"ie2=%d\n", ie2); */
/*       fprintf(stderr,"e=%lf\n", e); */
      
/*       fprintf(stderr,"inb1=%d\n", inb1); */
/*       fprintf(stderr,"inb2=%d\n", inb2); */
/*       fprintf(stderr,"nb=%lf\n", rhob); */
      
      pa1 = pressure2[inb1][ie1];
      pa2 = pressure2[inb2][ie1];
      
      //     fprintf(stderr,"e1=%lf\n", ie1*deltaEPP2+EPP2);
      //fprintf(stderr,"e2=%lf\n", (ie1+1)*deltaEPP2+EPP2);
      
      
      frace = (e-(ie1*deltaEPP2+EPP2))/deltaEPP2; 
      //fprintf(stderr,"frace=%lf\n", frace);
      fracnb = (rhob-(inb1*deltaBNP2+BNP2))/deltaBNP2; 
      
      pa = pa1*(1-fracnb) + pa2*fracnb;
      
      pb1 = pressure2[inb1][ie2];
      pb2 = pressure2[inb2][ie2];
      
      pb = pb1*(1-fracnb) + pb2*fracnb;
      
   /*    fprintf(stderr,"fracnb=%lf\n", fracnb); */
/*       fprintf(stderr,"inb1=%d\n", inb1); */
/*       fprintf(stderr,"inb2=%d\n", inb2); */
/*       fprintf(stderr,"nb1=%lf\n", (inb1*deltaBNP2+BNP2)); */
/*       fprintf(stderr,"nb2=%lf\n", (inb2*deltaBNP2+BNP2)); */
/*       fprintf(stderr,"nb=%lf\n", rhob); */
      
      p = pa*(1-frace) + pb*frace;
      
      //fprintf(stderr,"pa(ie1)=%lf\n", pa1);
      //fprintf(stderr,"pa(ie2)=%lf\n", pa2);
      
      //fprintf(stderr,"pa(e)=%lf\n", pa);  
      //fprintf(stderr,"pb(e)=%lf\n", pb);
      //fprintf(stderr,"p(e)=%lf\n", p);
  
    }
  return p/hbarc;
}


double EOS::interpolate2(double e, double rhob, int selector)
{
  double p, pa, pb, pa1, pa2, pb1, pb2;
  //use linear interpolation
  int ie1, ie2, inb1, inb2, NEps;
  double frace, fracnb;
  double eps0, deltaEps;
  double **array; 
  //selector = 0 : pressure
  //selector = 1 : temperature
  //selector = 2 : entropy density
  //selector = 3 : QGP fraction
  
  //if(rhob>0.00001) fprintf(stderr,"rhob=%lf\n", rhob);
  //fprintf(stderr,"e=%f\n",e);
  
  e*=hbarc; // in the files epsilon is in GeV/fm^3
  
  //fprintf(stderr,"e=%f\n",e);

  if(e<EPP2) //use first file for small epsilon values
    {
      if(e<EPP1) 
	{
	  ie1 = 0;
	  ie2 = 1;
	  frace = e/(EPP1);
	}
      else
	{
	  ie1 = floor((e-EPP1)/deltaEPP1);
	  ie2 = floor((e-EPP1)/deltaEPP1+1);
	  frace = (e-(ie1*deltaEPP1+EPP1))/deltaEPP1; 
	}
      
      if(ie1>NEPP1)
	{
	  fprintf(stderr,"ERROR in interpolate2. out of range.\n");
	  fprintf(stderr,"ie1=%d,NEPP1=%d\n", ie1, NEPP1);
	  exit(0);
	}
      if(ie2>NEPP1)
	{
	  fprintf(stderr,"ERROR in interpolate2. out of range.\n");
	  fprintf(stderr,"ie2=%d,NEPP1=%d\n", ie2, NEPP1);
	  exit(0);
	}

      switch (selector) 
	{
	case 0: array = pressure1; break;
	case 1: array = temperature1; break;
	case 2: array = entropyDensity1; break;
	case 3: array = QGPfraction1; break;
	default: fprintf(stderr,"ERROR in interpolate2 - selector must be 0,1,2, or 3\n"); exit(1);
	}

      pa = array[0][ie1];
      pb = array[0][ie2];
      
      if(e<EPP1) 
	{
	  p = pa*(frace);
	  //if (p<0) fprintf(stderr,"pa=%lf\n", pa);
	  //if (p<0) fprintf(stderr,"p=%lf\n", p);
	}
      else
	{
	  p = pa*(1-frace) + pb*frace;
	}
    }
  
  else // use other files for larger epsilon values
    {
      if (e<EPP3)
	{
	  eps0 = EPP2;
	  NEps = NEPP2;
	  deltaEps = deltaEPP2;
	  switch (selector) 
	    {
	    case 0: array = pressure2; break;
	    case 1: array = temperature2; break;
	    case 2: array = entropyDensity2; break;
	    case 3: array = QGPfraction2; break;
	    default: fprintf(stderr,"ERROR in interpolate2 - selector must be 0,1,2, or 3\n"); exit(1);
	    }
	}
      else if (e<EPP4)
	{
	  eps0 = EPP3;
	  NEps = NEPP3;
	  deltaEps = deltaEPP3;
	  switch (selector) 
	    {
	    case 0: array = pressure3; break;
	    case 1: array = temperature3; break;
	    case 2: array = entropyDensity3; break;
	    case 3: array = QGPfraction3; break;
	    default: fprintf(stderr,"ERROR in interpolate2 - selector must be 0,1,2, or 3\n"); exit(1);
	    }
	}
      else
	{
	  eps0 = EPP4;
	  NEps = NEPP4;
	  deltaEps = deltaEPP4;
	  switch (selector) 
	    {
	    case 0: array = pressure4; break;
	    case 1: array = temperature4; break;
	    case 2: array = entropyDensity4; break;
	    case 3: array = QGPfraction4; break;
	    default: fprintf(stderr,"ERROR in interpolate2 - selector must be 0,1,2, or 3\n"); exit(1);
	    }
	}

      ie1 = floor((e-eps0)/deltaEps);
      ie2 = floor((e-eps0)/deltaEps+1);
      
      if(ie1>NEps)
	{
	  fprintf(stderr,"ERROR in inperpolate2. out of range.\n");
	  fprintf(stderr,"ie1=%d,NEPP2=%d\n", ie1, NEps);
	  fprintf(stderr,"e=%f,eps0=%f; maxe=%f, deltaEps=%f\n", e, eps0, NEps*deltaEps+eps0, deltaEps);
	  exit(0);
	}
      if(ie2>NEps)
	{
	  fprintf(stderr,"ERROR in inperpolate2. out of range.\n");
	  fprintf(stderr,"ie2=%d,NEPP2=%d\n", ie2, NEps);
	  exit(0);
	}

      pa = array[0][ie1];
           
      frace = (e-(ie1*deltaEps+eps0))/deltaEps; 
      
      pb = array[0][ie2];
      
      p = pa*(1-frace) + pb*frace;
      //fprintf(stderr,"p=%f\n",p);
    }
  switch (selector) 
    {
    case 0: p/=hbarc; break;
    case 1: p/=hbarc; break;
    case 2: break;
    case 3: break;
    default: fprintf(stderr,"ERROR in interpolate2 - selector must be 0,1,2, or 3\n"); exit(1);
    }

  return p;
}


double EOS::get_dpOverde(double e, double rhob)
{
  double dp, pL, pR;
  //use linear interpolation
  double eLeft, eRight;
  
  e*=hbarc;

  //fprintf(stderr,"EPP2=%lf\n", EPP2);
  if(e<EPP2) //use first file for small epsilon values
    {
      eLeft=e-deltaEPP1/2.;
      eRight=e+deltaEPP1/2.;
      if(eLeft<EPP1)
	{
	  eLeft=EPP1;
	  eRight=EPP1+deltaEPP1;
	}
      
      if(eRight>EPP1+NEPP1*deltaEPP1)
	{
	  //eLeft=EPP1+(NEPP1-1)*deltaEPP1;
	  eRight=EPP1+NEPP1*deltaEPP1;
	  //fprintf(stderr,"e=%lf, EPP2=%lf\n", e, EPP2); 
	  //fprintf(stderr,"eLeft=%lf\n", eLeft); 
	  //fprintf(stderr,"eRight=%lf\n", eRight); 
	  //if(eRight-eLeft<deltaEPP1) fprintf(stderr,"delta=%lf, EPP1=%lf\n", eRight-eLeft, deltaEPP1); 
	}

      pL = p_func(eLeft/hbarc, rhob);
      pR = p_func(eRight/hbarc,  rhob);
      
      dp = pR-pL;
      
       if(dp<0) 
       { 
 	  fprintf(stderr,"1dp/de=%lf\n", dp/((eRight-eLeft)/hbarc)); 
 	  fprintf(stderr,"1pL=%lf\n", pL); 
 	  fprintf(stderr,"1pR=%lf\n", pR); 
 	  fprintf(stderr,"1e=%lf\n", e);  
 	  fprintf(stderr,"1eLeft=%lf\n", eLeft/hbarc);  
 	  fprintf(stderr,"1eRight=%lf\n", eRight/hbarc);  
       } 
      //if(dp<0) dp=0.;
	
      return dp/((eRight-eLeft)/hbarc);
    } 
  else
    {
      eLeft=e-deltaEPP2;
      eRight=e+deltaEPP2;
      if(eLeft<EPP2)
	{
	  eLeft=EPP2;
	  eRight=EPP2+deltaEPP2;
	}
      
      if(eRight>EPP2+NEPP2*deltaEPP2)
	{
	  eRight=EPP2+NEPP2*deltaEPP2;
	}
      
      pL = p_func(eLeft/hbarc, rhob);
      pR = p_func(eRight/hbarc, rhob);
      
      dp = pR-pL;
      
      if(dp<0 && pL*hbarc>0.1 && eLeft<24.)
	{
 	  fprintf(stderr,"*2dp=%lf\n", dp/hbarc);
 	  fprintf(stderr,"2dp/de=%lf\n", dp/((eRight-eLeft)/hbarc));
 	  fprintf(stderr,"2pL=%lf\n", pL*hbarc);
 	  fprintf(stderr,"2pR=%lf\n", pR*hbarc);
 	  fprintf(stderr,"2rhob=%lf\n", rhob);
 	  fprintf(stderr,"2e=%lf\n", e);
 	  fprintf(stderr,"EPP2=%lf\n", EPP2);
 	  fprintf(stderr,"2eLeft=%lf\n", eLeft);
 	  fprintf(stderr,"2eRight=%lf\n", eRight);
 	}  
      //fprintf(stderr,"dp/de=%lf\n", dp/((eRight-eLeft)/hbarc));  
      if(dp<0) dp=0.;

      return dp/((eRight-eLeft)/hbarc);
    }
}

double EOS::get_dpOverdrhob(double e, double rhob)
{
  double dp, pL, pR;
  //use linear interpolation
  double rhoLeft, rhoRight;
  
  e*=hbarc;

  //fprintf(stderr,"EPP2=%lf\n", EPP2);
  if(e<EPP2) //use first file for small epsilon values
    {
      rhoLeft=rhob-deltaBNP1/2.;
      rhoRight=rhob+deltaBNP1/2.;
      
      if(rhoLeft<BNP1)
	{
	  rhoLeft=BNP1;
	  rhoRight=BNP1+deltaBNP1;
	}
      
      if(rhoRight>BNP1+NBNP1*deltaBNP1)
	{
	  rhoRight=BNP1+NBNP1*deltaBNP1;
	}
/*       if(rhob>0) */
/* 	{ */
/* 	  fprintf(stderr,"rhoLeft=%lf\n", rhoLeft); */
/* 	  fprintf(stderr,"rhoRight=%lf\n", rhoRight); */
/* 	  fprintf(stderr,"rho=%lf\n", rhob); */
/* 	  fprintf(stderr,"BNP1=%lf\n", BNP1); */
/* 	} */
      pL = p_func(e/hbarc, rhoLeft);
      pR = p_func(e/hbarc, rhoRight);

      dp = pR-pL;
      return dp/((rhoRight-rhoLeft));
    } 
  else
    {
      rhoLeft=rhob-deltaBNP2/2.;
      rhoRight=rhob+deltaBNP2/2.;
      if(rhoLeft<BNP2)
	{
	  rhoLeft=BNP2;
	}
      
      if(rhoRight>BNP2+NBNP2*deltaBNP2)
	{
	  rhoRight=BNP2+NBNP2*deltaBNP2;
	}

      pL = p_func(e/hbarc, rhoLeft); 
      pR = p_func(e/hbarc, rhoRight);

/*       fprintf(stderr,"rhoLeft=%lf\n", rhoLeft); */
/*       fprintf(stderr,"rhoRight=%lf\n", rhoRight); */
/*       fprintf(stderr,"rho=%lf\n", rhob); */
/*       fprintf(stderr,"BNP2=%lf\n", BNP2); */
/*       fprintf(stderr,"pL=%lf\n", pL*hbarc); */
/*       fprintf(stderr,"pR=%lf\n", pR*hbarc); */
      
      dp = pR-pL;
      //      fprintf(stderr,"dp=%lf\n", dp);
      //if(dp/(rhoRight-rhoLeft)<1.) fprintf(stderr,"dp/drho2=%lf\n", dp/((rhoRight-rhoLeft)));
      
      return dp/((rhoRight-rhoLeft));
    }
}

double EOS::get_dpOverde2(double e, double rhob)
{
  double dp, pL, pR, deltaEps, eLeft, eRight;
  //use linear interpolation
  //CFY trying something: are the units right?:
  e *= hbarc;

  if (e<EPP2)
    {
      deltaEps = deltaEPP1;
    }
  else if (e<EPP3)
    {
      deltaEps = deltaEPP2;
    }
  else if (e<EPP4)
    {
      deltaEps = deltaEPP3;
    }
  else
    {
      deltaEps = deltaEPP4;
    }

  //CFY again:
  e /= hbarc;
  deltaEps /= hbarc;

  eLeft=e-deltaEps/2.;
  eRight=e+deltaEps/2.;
   
  //if(eLeft<EPP1)
  if(eLeft*hbarc<EPP1)
    {
      //eLeft=EPP1;
      //eRight=EPP1+deltaEPP1;
      //CFY again:
      eLeft=EPP1/hbarc;
      eRight=(EPP1+deltaEPP1)/hbarc;
      deltaEps = deltaEPP1/hbarc;
      //cout << "eLeft = " << eLeft << ", eRight = " << eRight << endl;
    }
   
  pL = p_func(eLeft, rhob);
  pR = p_func(eRight,  rhob);
  //cout << "pL = " << pL << ", pR = " << pR << endl;
  dp = pR-pL;
      
  if(dp<0) 
    { 
      fprintf(stderr,"1dp/de=%lf\n", dp/((eRight-eLeft)/hbarc)); 
      fprintf(stderr,"1pL=%lf\n", pL); 
      fprintf(stderr,"1pR=%lf\n", pR); 
      fprintf(stderr,"1e=%lf\n", e);  
      fprintf(stderr,"1eLeft=%lf\n", eLeft/hbarc);  
      fprintf(stderr,"1eRight=%lf\n", eRight/hbarc);  
    } 
   
  //return dp/((deltaEps)/hbarc);
  //CFY again:
  return dp/((deltaEps));
}

double EOS::get_dsOverde2(double e, double rhob, InitData *DATA)
{
  double dp, pL, pR, ds, sL, sR, deltaEps, eLeft, eRight;
  //use linear interpolation
  //CFY trying something: are the units right?:
  e *= hbarc;

  if (e<EPP2)
    {
      deltaEps = deltaEPP1;
    }
  else if (e<EPP3)
    {
      deltaEps = deltaEPP2;
    }
  else if (e<EPP4)
    {
      deltaEps = deltaEPP3;
    }
  else
    {
      deltaEps = deltaEPP4;
    }

  //CFY again:
  e /= hbarc;
  deltaEps /= hbarc;

  eLeft=e-deltaEps/2.;
  eRight=e+deltaEps/2.;

  //if(eLeft<EPP1)
  if(eLeft*hbarc<EPP1)
    {
      //eLeft=EPP1;
      //eRight=EPP1+deltaEPP1;
      //CFY again:
      eLeft=EPP1/hbarc;
      eRight=(EPP1+deltaEPP1)/hbarc;
      deltaEps = deltaEPP1/hbarc;
      //cout << "eLeft = " << eLeft << ", eRight = " << eRight << endl;
    }

  pL = p_func(eLeft, rhob);
  pR = p_func(eRight, rhob);
  sL = s_func(eLeft, pL, rhob, DATA);
  sR = s_func(eRight, pR,  rhob, DATA);
  //cout << "pL = " << pL << ", pR = " << pR << endl;
  ds = sR-sL;
      
  if(ds<0) 
    { 
      fprintf(stderr,"1ds/de=%lf\n", ds/((eRight-eLeft)/hbarc)); 
      fprintf(stderr,"1sL=%lf\n", sL); 
      fprintf(stderr,"1sR=%lf\n", sR); 
      fprintf(stderr,"1e=%lf\n", e);  
      fprintf(stderr,"1eLeft=%lf\n", eLeft/hbarc);  
      fprintf(stderr,"1eRight=%lf\n", eRight/hbarc);  
    } 
   
  //return dp/((deltaEps)/hbarc);
  //CFY again:
  return ds/((deltaEps));
}

double EOS::p_func(double e, double rhob)
{
 double f;
 if (whichEOS==0)
   f = cs2*e;
 else if (whichEOS==1)
   f = interpolate_pressure(e,rhob);
 else if (whichEOS==2)
   f = interpolate2(e,rhob,0); //selector 0 means get pressure 
 return f;
}/* p_func */

double EOS::p_rho_func(double e, double rhob)
{
 double f;
 if (whichEOS==0) 
   f = 0.0;
 else if (whichEOS==1)
   f = get_dpOverdrhob(e, rhob);
 else if (whichEOS==2)
   f = 0.0;
 return f;
}/* p_rho_func */

double EOS::p_e_func(double e, double rhob)
{
 double f;
 
 if (whichEOS==0)
   f = cs2;
 else if (whichEOS==1)
   f = get_dpOverde(e, rhob);
 else if (whichEOS==2)
   f = get_dpOverde2(e, rhob);
 return f;
}/* p_e_func */


double EOS::interpolate(double e, double rhob, int selector) // for eos 1, get temperature or mu_B
{
  // selector=0 : return temperature
  // selector=1 : return baryon chemical potential
  double T, pa, pb, pa1, pa2, pb1, pb2;
  //use linear interpolation
  int ie1, ie2, inb1, inb2;
  double frace, fracnb;
  
  e*=hbarc; // in the files epsilon is in GeV/fm^3

  if(e>24.)
    {
      if (selector==0)
	{
	  return pow(((e-.3642)*120./(pow(PI,2.))/169.*(pow(hbarc,3.))),0.25)/hbarc;
	}    
      else if (selector==1)
	{
	  return 0.;
	}
    }
  if(e<EPP2) //use first file for small epsilon values
    {
      if(e<EPP1) 
	{
	  ie1 = 0;
	  ie2 = 1;
	  inb1 = floor((rhob-BNP1)/deltaBNP1);
	  inb2 = floor((rhob-BNP1)/deltaBNP1+1);
	  frace = e/(EPP1);
  	}
      else
	{
	  ie1 = floor((e-EPP1)/deltaEPP1);
	  ie2 = floor((e-EPP1)/deltaEPP1+1);
	  inb1 = floor((rhob-BNP1)/deltaBNP1);
	  inb2 = floor((rhob-BNP1)/deltaBNP1+1);
	  frace = (e-(ie1*deltaEPP1+EPP1))/deltaEPP1; 
  	}
      
      if(ie1>NEPP1 || inb1>NBNP1)
	{
	  fprintf(stderr,"ERROR in inperpolate_temperature. out of range.\n");
	  fprintf(stderr,"ie1=%d,NEPP1=%d; inb1=%d, NBNP1=%d\n", ie1, NEPP1, inb1, NBNP1);
	  exit(0);
	}
      if(ie2>NEPP1 || inb2>NBNP1)
	{
	  fprintf(stderr,"ERROR in inperpolate_temperature. out of range.\n");
	  fprintf(stderr,"ie2=%d,NEPP1=%d; inb2=%d, NBNP1=%d\n", ie2, NEPP1, inb2, NBNP1);
	  exit(0);
	}
      
      if (selector==0)
	{
	  pa1 = temperature1[inb1][ie1];
	  pa2 = temperature1[inb2][ie1];
	}
      else if (selector==1)
	{
	  pa1 = mu1[inb1][ie1];
	  pa2 = mu1[inb2][ie1];
	}

      if (pa1<0) pa1=0.;
      if (pa2<0) pa2=0.;
    
      fracnb = (rhob-(inb1*deltaBNP1+BNP1))/deltaBNP1; 
      
      pa = pa1*(1-fracnb) + pa2*fracnb;
      
      if (selector==0)
	{
	  pb1 = temperature1[inb1][ie2];
	  pb2 = temperature1[inb2][ie2];
	}
      else if (selector==1)
	{
	  pb1 = mu1[inb1][ie2];
	  pb2 = mu1[inb2][ie2];
	}

      if (pb1<0) pb1=0.;
      if (pb2<0) pb2=0.;
      
      pb = pb1*(1-fracnb) + pb2*fracnb;
     
      //fprintf(stderr,"fracnb=%lf\n", fracnb);

      if(e<EPP1) 
	{
	  T = pa*(frace);
	}
      else 
	{
	  T = pa*(1-frace) + pb*frace;
	  
	}
    }
  else // use second file for larger epsilon values
    {
      ie1 = floor((e-EPP2)/deltaEPP2);
      ie2 = floor((e-EPP2)/deltaEPP2+1);
      
      inb1 = floor((rhob-BNP2)/deltaBNP2);
      inb2 = floor((rhob-BNP2)/deltaBNP2+1);
      
      if(ie1>NEPP2 || inb1>NBNP2)
	{
	  fprintf(stderr,"ERROR in inperpolate_temperature. out of range.\n");
	  fprintf(stderr,"ie1=%d,NEPP2=%d; inb1=%d, NBNP2=%d\n", ie1, NEPP2, inb1, NBNP2);
	  fprintf(stderr,"e=%f,EPP2=%f; maxe=%f, deltaEPP2=%f\n", e, EPP2, NEPP2*deltaEPP2+EPP2, deltaEPP2);
	  exit(0);
	}
      if(ie2>NEPP2 || inb2>NBNP2)
	{
	  fprintf(stderr,"ERROR in inperpolate_temperature. out of range.\n");
	  fprintf(stderr,"ie2=%d,NEPP2=%d; inb2=%d, NBNP2=%d\n", ie2, NEPP2, inb2, NBNP2);
	  exit(0);
	}

      if (selector==0)
	{
	  pa1 = temperature2[inb1][ie1];
	  pa2 = temperature2[inb2][ie1];
	}
      else if (selector==1)
	{
	  pa1 = mu2[inb1][ie1];
	  pa2 = mu2[inb2][ie1];
	}
      
      if (pa1<0) pa1=0.;
      if (pa2<0) pa2=0.;
      
      frace = (e-(ie1*deltaEPP2+EPP2))/deltaEPP2; 
      fracnb = (rhob-(inb1*deltaBNP2+BNP2))/deltaBNP2; 
      
      pa = pa1*(1-fracnb) + pa2*fracnb;
      
      if (selector==0)
	{
	  pb1 = temperature2[inb1][ie2];
	  pb2 = temperature2[inb2][ie2];
	}     
      else if (selector==1)
	{
	  pb1 = mu2[inb1][ie2];
	  pb2 = mu2[inb2][ie2];
	}

      if (pb1<0) pb1=0.;
      if (pb2<0) pb2=0.;
      
      pb = pb1*(1-fracnb) + pb2*fracnb;

      T = pa*(1-frace) + pb*frace;
    }
  return T/hbarc;
}

double EOS::s_func(double epsilon, double p, double rhob, InitData *DATA)
{
 double f;
 double P, T, mu;
  
 P = p_func(epsilon, rhob);

 if (whichEOS==0)
   {
    f = pow(epsilon + p, 0.75)/2.1365;
    return f;
   }
 else if (whichEOS==1)
   {
     T = interpolate(epsilon, rhob, 0);
     mu = interpolate(epsilon, rhob, 1);
   }
 else if (whichEOS==2)
   {
     T = interpolate2(epsilon, rhob, 1);
     mu = 0.0;
   }
 
 if (T!=0)
   f = (epsilon + P - mu*rhob)/T;
 else
   f = 0.;

 return f;
 
}/* s_func */

double EOS::ssolve(double e, double rhob, double s)
{
  // takes e in GeV/fm^3 and passes it on in 1/fm^4 ...
  double P, T, mu;
  P = p_func(e/hbarc, rhob);
  if (whichEOS==1)
    {
      T = interpolate(e/hbarc, rhob, 0);
      mu = interpolate(e/hbarc, rhob, 1);
    }
  else if (whichEOS==2)
    {
      T = interpolate2(e/hbarc, rhob, 1);
      mu = 0.0;
    }
      
  //fprintf(stderr,"T=%f\n",T);
  //fprintf(stderr,"P=%f\n",P);
  //fprintf(stderr,"mu=%f\n",mu);
  //fprintf(stderr,"ssolve=%f\n",T*s-e-P+mu*rhob);

  return T*s-e-P+mu*rhob;
}

double EOS::findRoot(double (EOS::*func)(double, double, double), double rhob, double s, double e1, double e2, double eacc)
{
  int j, jmax;
  jmax=40;
  double emid, de, value, f, fmid;
  fmid = (this->*func)(e2, rhob, s);
  f = (this->*func)(e1, rhob, s);
 
  //  fprintf(stderr,"fmid=%f\n",fmid);
  //fprintf(stderr,"fabs(f)=%f\n",fabs(f));
  //fprintf(stderr,"eacc=%f\n",eacc);
  
  if(f*fmid>=0)
    {
      if( fabs(f) < eacc )
	{
	 return 0.;
	}
      fprintf(stderr,"root must be bracketed in findRoot\n");
      fprintf(stderr,"f=%f\n",f);
      fprintf(stderr,"fmid=%f\n",fmid);
      fprintf(stderr,"fabs(f)=%f\n",fabs(f));
      fprintf(stderr,"eacc=%f\n",eacc);
    }
     
  if (f<0)
    {
      value=e1;
      de=e2-e1;
    }
  else
    {
      value=e2;
      de=e1-e2;
    }
  for(j=1; j<=jmax; j++)
    {
      de*=0.5;
      emid = value+de;
      fmid = (this->*func)(emid, rhob, s);
      //fprintf(stderr,"fmid(emid)=%f\n",fmid);
      //fprintf(stderr,"emid=%f\n",emid);
      //fprintf(stderr,"value=%f\n",value);
      if (fmid<=0.) value = emid;
      if (fabs(de)<eacc || fmid==0.) return value;
    }
  fprintf(stderr,"too many bisections in findRoot\n");
}

