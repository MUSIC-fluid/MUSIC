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
  gsl_interp_free(interp_s2e);
  gsl_interp_accel_free(accel_s2e);
}

void EOS::checkForReadError(FILE *file, const char* name)
{
  if(!(file))
    {
      fprintf(stderr, "file %s not found.\n", name);
      fprintf(stderr, "Exiting...\n");
      exit(0);
    }
}

void EOS::init_eos0() {
  whichEOS = 0;
}

void EOS::init_eos()
{
  // read the azhydro pressure, temperature, and 
  // baryon chemical potential from file
  whichEOS = 1;
  fprintf(stderr,"reading EOS... \n");
//   int bytes_read;
  int i, j;
  FILE *eos_p1, *eos_p2;
  FILE *eos_T1, *eos_T2;
  FILE *eos_mu1, *eos_mu2;
  const char* EOSPATH = "HYDROPROGRAMPATH";
  char* envPath;
  envPath = util->char_malloc(100);
  envPath = getenv(EOSPATH);
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
      strcat(eos_p1_name,"/EOS/EOS-Q/aa1_p.dat");

      eos_p2_name = util->char_malloc(100);
      strcat(eos_p2_name,envPath);
      strcat(eos_p2_name,"/EOS/EOS-Q/aa2_p.dat");

      eos_T1_name = util->char_malloc(100);
      strcat(eos_T1_name,envPath);
      strcat(eos_T1_name,"/EOS/EOS-Q/aa1_t.dat");

      eos_T2_name = util->char_malloc(100);
      strcat(eos_T2_name,envPath);
      strcat(eos_T2_name,"/EOS/EOS-Q/aa2_t.dat");

      eos_mu1_name = util->char_malloc(100);
      strcat(eos_mu1_name,envPath);
      strcat(eos_mu1_name,"/EOS/EOS-Q/aa1_mb.dat");
 
      eos_mu2_name = util->char_malloc(100);
      strcat(eos_mu2_name,envPath);
      strcat(eos_mu2_name,"/EOS/EOS-Q/aa2_mb.dat");
    }
  else
    {
      eos_p1_name = util->char_malloc(100);
      strcat(eos_p1_name,".");
      strcat(eos_p1_name,"/EOS/EOS-Q/aa1_p.dat");

      eos_p2_name = util->char_malloc(100);
      strcat(eos_p2_name,".");
      strcat(eos_p2_name,"/EOS/EOS-Q/aa2_p.dat");

      eos_T1_name = util->char_malloc(100);
      strcat(eos_T1_name,".");
      strcat(eos_T1_name,"/EOS/EOS-Q/aa1_t.dat");

      eos_T2_name = util->char_malloc(100);
      strcat(eos_T2_name,".");
      strcat(eos_T2_name,"/EOS/EOS-Q/aa2_t.dat");

      eos_mu1_name = util->char_malloc(100);
      strcat(eos_mu1_name,".");
      strcat(eos_mu1_name,"/EOS/EOS-Q/aa1_mb.dat");
 
      eos_mu2_name = util->char_malloc(100);
      strcat(eos_mu2_name,".");
      strcat(eos_mu2_name,"/EOS/EOS-Q/aa2_mb.dat");
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
  fscanf(eos_p1,"%lf %lf",&BNP1,&EPP1);
  fscanf(eos_p1,"%lf %lf %d %d",&deltaBNP1,&deltaEPP1,&NBNP1,&NEPP1);
  fscanf(eos_p2,"%lf %lf",&BNP2,&EPP2);
  fscanf(eos_p2,"%lf %lf %d %d",&deltaBNP2,&deltaEPP2,&NBNP2,&NEPP2);
  fscanf(eos_T1,"%lf %lf",&BNP1,&EPP1);
  fscanf(eos_T1,"%lf %lf %d %d",&deltaBNP1,&deltaEPP1,&NBNP1,&NEPP1);
  fscanf(eos_T2,"%lf %lf",&BNP2,&EPP2);
  fscanf(eos_T2,"%lf %lf %d %d",&deltaBNP2,&deltaEPP2,&NBNP2,&NEPP2);
  fscanf(eos_mu1,"%lf %lf",&BNP1,&EPP1);
  fscanf(eos_mu1,"%lf %lf %d %d",&deltaBNP1,&deltaEPP1,&NBNP1,&NEPP1);
  fscanf(eos_mu2,"%lf %lf",&BNP2,&EPP2);
  fscanf(eos_mu2,"%lf %lf %d %d",&deltaBNP2,&deltaEPP2,&NBNP2,&NEPP2);
 
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
  
  // allocate memory for baryon chemical potential arrays
  cs2_1 = util->mtx_malloc(NBNP1+1, NEPP1+1);
  cs2_2 = util->mtx_malloc(NBNP2+1, NEPP2+1);

  // read pressure, temperature and chemical potential values
  for(j=0;j<=NEPP1;j++)
    for(i=0;i<=NBNP1;i++)
      {
	fscanf(eos_p1,"%lf",&pressure1[i][j]);
	fscanf(eos_T1,"%lf",&temperature1[i][j]);
	fscanf(eos_mu1,"%lf",&mu1[i][j]);
      }

  for(j=0;j<=NEPP2;j++)
    for(i=0;i<=NBNP2;i++)
      {
	fscanf(eos_p2,"%lf",&pressure2[i][j]);
	fscanf(eos_T2,"%lf",&temperature2[i][j]);
	fscanf(eos_mu2,"%lf",&mu2[i][j]);
      }

  build_velocity_of_sound_sq_matrix();
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
//   int bytes_read;
  int i, j;
  FILE *eos_d1, *eos_d2, *eos_d3, *eos_d4, *eos_d5, *eos_d6, *eos_d7;
  FILE *eos_T1, *eos_T2, *eos_T3, *eos_T4, *eos_T5, *eos_T6, *eos_T7;
  const char* EOSPATH = "HYDROPROGRAMPATH";
  char* envPath;
  envPath = util->char_malloc(100);
  envPath = getenv(EOSPATH);
  char* eos_d1_name;
  char* eos_d2_name;
  char* eos_d3_name;
  char* eos_d4_name;
  char* eos_d5_name;
  char* eos_d6_name;
  char* eos_d7_name;
  char* eos_T1_name;
  char* eos_T2_name;
  char* eos_T3_name;
  char* eos_T4_name;
  char* eos_T5_name;
  char* eos_T6_name;
  char* eos_T7_name;
  double eps, baryonDensity; //dummies for now
  char* temp;
  temp = new char[80];//util->char_malloc(100);

    
 if (envPath != 0 && *envPath != '\0') // if path is set in the environment
    {
      envPath[strlen(envPath)]='\0';
      fprintf(stderr,"from path %s/EOS/s95p-v1 \n", envPath);
      
      eos_d1_name = new char[strlen(envPath)+strlen("/EOS/s95p-v1/s95p-v1_dens1.dat")+1];//util->char_malloc(100);
      eos_d1_name[0] = '\0';
      strcat(eos_d1_name,envPath);
      strcpy(temp,"/EOS/s95p-v1/s95p-v1_dens1.dat");
      strcat(eos_d1_name,temp);
 
      
      strcpy(temp,"/EOS/s95p-v1/s95p-v1_dens2.dat");
      eos_d2_name = new char[strlen(envPath)+strlen("/EOS/s95p-v1/s95p-v1_dens2.dat")+1];//util->char_malloc(100);
      eos_d2_name[0] = '\0';
      strcat(eos_d2_name,envPath);
      strcat(eos_d2_name,temp);
  
      strcpy(temp,"/EOS/s95p-v1/s95p-v1_dens3.dat");
      eos_d3_name = new char[strlen(envPath)+strlen("/EOS/s95p-v1/s95p-v1_dens3.dat")+1];//util->char_malloc(100);
      eos_d3_name[0] = '\0';
      strcat(eos_d3_name,envPath);
      strcat(eos_d3_name,temp);
      
      strcpy(temp,"/EOS/s95p-v1/s95p-v1_dens4.dat");
      eos_d4_name = new char[strlen(envPath)+strlen("/EOS/s95p-v1/s95p-v1_dens4.dat")+1];//util->char_malloc(100);
      eos_d4_name[0] = '\0';
      strcat(eos_d4_name,envPath);
      strcat(eos_d4_name,temp);
      
      strcpy(temp,"/EOS/s95p-v1/s95p-v1_dens5.dat");
      eos_d5_name = new char[strlen(envPath)+strlen("/EOS/s95p-v1/s95p-v1_dens5.dat")+1];//util->char_malloc(100);
      eos_d5_name[0] = '\0';
      strcat(eos_d5_name,envPath);
      strcat(eos_d5_name,temp);
      
      strcpy(temp,"/EOS/s95p-v1/s95p-v1_dens6.dat");
      eos_d6_name = new char[strlen(envPath)+strlen("/EOS/s95p-v1/s95p-v1_dens6.dat")+1];//util->char_malloc(100);
      eos_d6_name[0] = '\0';
      strcat(eos_d6_name,envPath);
      strcat(eos_d6_name,temp);
      
      strcpy(temp,"/EOS/s95p-v1/s95p-v1_dens7.dat");
      eos_d7_name = new char[strlen(envPath)+strlen("/EOS/s95p-v1/s95p-v1_dens7.dat")+1];//util->char_malloc(100);
      eos_d7_name[0] = '\0';
      strcat(eos_d7_name,envPath);
      strcat(eos_d7_name,temp);
      
      strcpy(temp,"/EOS/s95p-v1/s95p-v1_par1.dat");
      eos_T1_name = new char[strlen(envPath)+strlen("/EOS/s95p-v1/s95p-v1_par1.dat")+1];//util->char_malloc(100);
      eos_T1_name[0] = '\0';
      strcat(eos_T1_name,envPath);
      strcat(eos_T1_name,temp);
      
      strcpy(temp,"/EOS/s95p-v1/s95p-v1_par2.dat");
      eos_T2_name = new char[strlen(envPath)+strlen("/EOS/s95p-v1/s95p-v1_par2.dat")+1];//util->char_malloc(100);
      eos_T2_name[0] = '\0';
      strcat(eos_T2_name,envPath);
      strcat(eos_T2_name,temp);
      
      strcpy(temp,"/EOS/s95p-v1/s95p-v1_par3.dat");
      eos_T3_name = new char[strlen(envPath)+strlen("/EOS/s95p-v1/s95p-v1_par3.dat")+1];//util->char_malloc(100);
      eos_T3_name[0] = '\0';
      strcat(eos_T3_name,envPath);
      strcat(eos_T3_name,temp);
      
      strcpy(temp,"/EOS/s95p-v1/s95p-v1_par4.dat");
      eos_T4_name = new char[strlen(envPath)+strlen("/EOS/s95p-v1/s95p-v1_par4.dat")+1];//util->char_malloc(100);
      eos_T4_name[0] = '\0';
      strcat(eos_T4_name,envPath);
      strcat(eos_T4_name,temp);
      
      strcpy(temp,"/EOS/s95p-v1/s95p-v1_par5.dat");
      eos_T5_name = new char[strlen(envPath)+strlen("/EOS/s95p-v1/s95p-v1_par5.dat")+1];//util->char_malloc(100);
      eos_T5_name[0] = '\0';
      strcat(eos_T5_name,envPath);
      strcat(eos_T5_name,temp);
      
      strcpy(temp,"/EOS/s95p-v1/s95p-v1_par6.dat");
      eos_T6_name = new char[strlen(envPath)+strlen("/EOS/s95p-v1/s95p-v1_par6.dat")+1];//util->char_malloc(100);
      eos_T6_name[0] = '\0';
      strcat(eos_T6_name,envPath);
      strcat(eos_T6_name,temp);
      
      strcpy(temp,"/EOS/s95p-v1/s95p-v1_par7.dat");
      eos_T7_name = new char[strlen(envPath)+strlen("/EOS/s95p-v1/s95p-v1_par7.dat")+1];//util->char_malloc(100);
      eos_T7_name[0] = '\0';
      strcat(eos_T7_name,envPath);
      strcat(eos_T7_name,temp);
    }
 else //if path is not set in the environment use current folder and then /EOS/s95p-v1 subfolder
    {
      eos_d1_name = util->char_malloc(300);
      strcat(eos_d1_name,".");
      strcat(eos_d1_name,"/EOS/s95p-v1/s95p-v1_dens1.dat");
      eos_d2_name = util->char_malloc(300);
      strcat(eos_d2_name,".");
      strcat(eos_d2_name,"/EOS/s95p-v1/s95p-v1_dens2.dat");
      eos_d3_name = util->char_malloc(300);
      strcat(eos_d3_name,".");
      strcat(eos_d3_name,"/EOS/s95p-v1/s95p-v1_dens3.dat");
      eos_d4_name = util->char_malloc(300);
      strcat(eos_d4_name,".");
      strcat(eos_d4_name,"/EOS/s95p-v1/s95p-v1_dens4.dat");
      eos_d5_name = util->char_malloc(300);
      strcat(eos_d5_name,".");
      strcat(eos_d5_name,"/EOS/s95p-v1/s95p-v1_dens5.dat");
      eos_d6_name = util->char_malloc(300);
      strcat(eos_d6_name,".");
      strcat(eos_d6_name,"/EOS/s95p-v1/s95p-v1_dens6.dat");
      eos_d7_name = util->char_malloc(300);
      strcat(eos_d7_name,".");
      strcat(eos_d7_name,"/EOS/s95p-v1/s95p-v1_dens7.dat");
      eos_T1_name = util->char_malloc(300);
      strcat(eos_T1_name,".");
      strcat(eos_T1_name,"/EOS/s95p-v1/s95p-v1_par1.dat");
      eos_T2_name = util->char_malloc(300);
      strcat(eos_T2_name,".");
      strcat(eos_T2_name,"/EOS/s95p-v1/s95p-v1_par2.dat");
      eos_T3_name = util->char_malloc(300);
      strcat(eos_T3_name,".");
      strcat(eos_T3_name,"/EOS/s95p-v1/s95p-v1_par3.dat");
      eos_T4_name = util->char_malloc(300);
      strcat(eos_T4_name,".");
      strcat(eos_T4_name,"/EOS/s95p-v1/s95p-v1_par4.dat");
      eos_T5_name = util->char_malloc(300);
      strcat(eos_T5_name,".");
      strcat(eos_T5_name,"/EOS/s95p-v1/s95p-v1_par5.dat");
      eos_T6_name = util->char_malloc(300);
      strcat(eos_T6_name,".");
      strcat(eos_T6_name,"/EOS/s95p-v1/s95p-v1_par6.dat");
      eos_T7_name = util->char_malloc(300);
      strcat(eos_T7_name,".");
      strcat(eos_T7_name,"/EOS/s95p-v1/s95p-v1_par7.dat");
    }
  
  eos_d1 = fopen(eos_d1_name, "r");
  eos_d2 = fopen(eos_d2_name, "r");
  eos_d3 = fopen(eos_d3_name, "r");
  eos_d4 = fopen(eos_d4_name, "r");
  eos_d5 = fopen(eos_d5_name, "r");
  eos_d6 = fopen(eos_d6_name, "r");
  eos_d7 = fopen(eos_d7_name, "r");
  eos_T1 = fopen(eos_T1_name, "r");
  eos_T2 = fopen(eos_T2_name, "r");
  eos_T3 = fopen(eos_T3_name, "r");
  eos_T4 = fopen(eos_T4_name, "r");
  eos_T5 = fopen(eos_T5_name, "r");
  eos_T6 = fopen(eos_T6_name, "r");
  eos_T7 = fopen(eos_T7_name, "r");
  
  checkForReadError(eos_d1,eos_d1_name);
  checkForReadError(eos_d2,eos_d2_name);
  checkForReadError(eos_d3,eos_d3_name);
  checkForReadError(eos_d4,eos_d4_name);
  checkForReadError(eos_d5,eos_d5_name);
  checkForReadError(eos_d6,eos_d6_name);
  checkForReadError(eos_d7,eos_d7_name);
  checkForReadError(eos_T1,eos_T1_name);
  checkForReadError(eos_T2,eos_T2_name);
  checkForReadError(eos_T3,eos_T3_name);
  checkForReadError(eos_T4,eos_T4_name);
  checkForReadError(eos_T5,eos_T5_name);
  checkForReadError(eos_T6,eos_T6_name);
  checkForReadError(eos_T7,eos_T7_name);
 
  //read the first two lines with general info:
  // lowest value of epsilon
  // deltaEpsilon, number of epsilon steps (i.e. # of lines)
  fscanf(eos_T1,"%lf",&EPP1);
  fscanf(eos_T1,"%lf %d",&deltaEPP1,&NEPP1);
  fscanf(eos_T2,"%lf",&EPP2);
  fscanf(eos_T2,"%lf %d",&deltaEPP2,&NEPP2);
  fscanf(eos_T3,"%lf",&EPP3);
  fscanf(eos_T3,"%lf %d",&deltaEPP3,&NEPP3);
  fscanf(eos_T4,"%lf",&EPP4);
  fscanf(eos_T4,"%lf %d",&deltaEPP4,&NEPP4);
  fscanf(eos_T5,"%lf",&EPP5);
  fscanf(eos_T5,"%lf %d",&deltaEPP5,&NEPP5);
  fscanf(eos_T6,"%lf",&EPP6);
  fscanf(eos_T6,"%lf %d",&deltaEPP6,&NEPP6);
  fscanf(eos_T7,"%lf",&EPP7);
  fscanf(eos_T7,"%lf %d",&deltaEPP7,&NEPP7);
  fscanf(eos_d1,"%lf",&EPP1);
  fscanf(eos_d1,"%lf %d",&deltaEPP1,&NEPP1);
  fscanf(eos_d2,"%lf",&EPP2);
  fscanf(eos_d2,"%lf %d",&deltaEPP2,&NEPP2);
  fscanf(eos_d3,"%lf",&EPP3);
  fscanf(eos_d3,"%lf %d",&deltaEPP3,&NEPP3);
  fscanf(eos_d4,"%lf",&EPP4);
  fscanf(eos_d4,"%lf %d",&deltaEPP4,&NEPP4);
  fscanf(eos_d5,"%lf",&EPP5);
  fscanf(eos_d5,"%lf %d",&deltaEPP5,&NEPP5);
  fscanf(eos_d6,"%lf",&EPP6);
  fscanf(eos_d6,"%lf %d",&deltaEPP6,&NEPP6);
  fscanf(eos_d7,"%lf",&EPP7);
  fscanf(eos_d7,"%lf %d",&deltaEPP7,&NEPP7);
  
 
  // no rho_b dependence at the moment
  NBNP1=0; 
  NBNP2=0;
  NBNP3=0;
  NBNP4=0;
  NBNP5=0;
  NBNP6=0;
  NBNP7=0;

  // allocate memory for pressure arrays
  pressure1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  pressure2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  pressure3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  pressure4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  pressure5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  pressure6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  pressure7=util->mtx_malloc(NBNP7+1,NEPP7+1);

  // allocate memory for entropy density arrays
  entropyDensity1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  entropyDensity2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  entropyDensity3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  entropyDensity4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  entropyDensity5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  entropyDensity6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  entropyDensity7=util->mtx_malloc(NBNP7+1,NEPP7+1);

  // allocate memory for QGP fraction arrays
  QGPfraction1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  QGPfraction2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  QGPfraction3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  QGPfraction4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  QGPfraction5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  QGPfraction6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  QGPfraction7=util->mtx_malloc(NBNP7+1,NEPP7+1);

  // allocate memory for temperature arrays
  temperature1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  temperature2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  temperature3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  temperature4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  temperature5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  temperature6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  temperature7=util->mtx_malloc(NBNP7+1,NEPP7+1);
  
  // allocate memory for velocity of sound squared arrays
  cs2_1 = util->mtx_malloc(NBNP1+1, NEPP1+1);
  cs2_2 = util->mtx_malloc(NBNP2+1, NEPP2+1);
  cs2_3 = util->mtx_malloc(NBNP3+1, NEPP3+1);
  cs2_4 = util->mtx_malloc(NBNP4+1, NEPP4+1);
  cs2_5 = util->mtx_malloc(NBNP5+1, NEPP5+1);
  cs2_6 = util->mtx_malloc(NBNP6+1, NEPP6+1);
  cs2_7 = util->mtx_malloc(NBNP7+1, NEPP7+1);

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
      fscanf(eos_d1,"%lf",&eps);
      fscanf(eos_d1,"%lf",&pressure1[i][j]);
      fscanf(eos_d1,"%lf",&entropyDensity1[i][j]);
      fscanf(eos_d1,"%lf",&baryonDensity);
      fscanf(eos_d1,"%lf",&QGPfraction1[i][j]);
      fscanf(eos_T1,"%lf",&temperature1[i][j]);
      fscanf(eos_T1,"%lf",&eps); //dummy
      fscanf(eos_T1,"%lf",&eps); //dummy
    }

  for(j=NEPP2-1; j>=0; j--)
    {
      fscanf(eos_d2,"%lf",&eps);
      fscanf(eos_d2,"%lf",&pressure2[i][j]);
      fscanf(eos_d2,"%lf",&entropyDensity2[i][j]);
      fscanf(eos_d2,"%lf",&baryonDensity);
      fscanf(eos_d2,"%lf",&QGPfraction2[i][j]);
      fscanf(eos_T2,"%lf",&temperature2[i][j]);
      fscanf(eos_T2,"%lf",&eps); //dummy
      fscanf(eos_T2,"%lf",&eps); //dummy
    }

  for(j=NEPP3-1; j>=0; j--)
    {
      fscanf(eos_d3,"%lf",&eps);
      fscanf(eos_d3,"%lf",&pressure3[i][j]);
      fscanf(eos_d3,"%lf",&entropyDensity3[i][j]);
      fscanf(eos_d3,"%lf",&baryonDensity);
      fscanf(eos_d3,"%lf",&QGPfraction3[i][j]);
      fscanf(eos_T3,"%lf",&temperature3[i][j]);
      fscanf(eos_T3,"%lf",&eps); //dummy
      fscanf(eos_T3,"%lf",&eps); //dummy
    }

  for(j=NEPP4-1; j>=0; j--)
    {
      fscanf(eos_d4,"%lf",&eps);
      fscanf(eos_d4,"%lf",&pressure4[i][j]);
      fscanf(eos_d4,"%lf",&entropyDensity4[i][j]);
      fscanf(eos_d4,"%lf",&baryonDensity);
      fscanf(eos_d4,"%lf",&QGPfraction4[i][j]);
      fscanf(eos_T4,"%lf",&temperature4[i][j]);
      fscanf(eos_T4,"%lf",&eps); //dummy
      fscanf(eos_T4,"%lf",&eps); //dummy
    }

  for(j=NEPP5-1; j>=0; j--)
    {
      fscanf(eos_d5,"%lf",&eps);
      fscanf(eos_d5,"%lf",&pressure5[i][j]);
      fscanf(eos_d5,"%lf",&entropyDensity5[i][j]);
      fscanf(eos_d5,"%lf",&baryonDensity);
      fscanf(eos_d5,"%lf",&QGPfraction5[i][j]);
      fscanf(eos_T5,"%lf",&temperature5[i][j]);
      fscanf(eos_T5,"%lf",&eps); //dummy
      fscanf(eos_T5,"%lf",&eps); //dummy
    } 

  for(j=NEPP6-1; j>=0; j--)
    {
      fscanf(eos_d6,"%lf",&eps);
      fscanf(eos_d6,"%lf",&pressure6[i][j]);
      fscanf(eos_d6,"%lf",&entropyDensity6[i][j]);
      fscanf(eos_d6,"%lf",&baryonDensity);
      fscanf(eos_d6,"%lf",&QGPfraction6[i][j]);
      fscanf(eos_T6,"%lf",&temperature6[i][j]);
      fscanf(eos_T6,"%lf",&eps); //dummy
      fscanf(eos_T6,"%lf",&eps); //dummy
    } 

  for(j=NEPP7-1; j>=0; j--)
    {
      fscanf(eos_d7,"%lf",&eps);
      fscanf(eos_d7,"%lf",&pressure7[i][j]);
      fscanf(eos_d7,"%lf",&entropyDensity7[i][j]);
      fscanf(eos_d7,"%lf",&baryonDensity);
      fscanf(eos_d7,"%lf",&QGPfraction7[i][j]);
      fscanf(eos_T7,"%lf",&temperature7[i][j]);
      fscanf(eos_T7,"%lf",&eps); //dummy
      fscanf(eos_T7,"%lf",&eps); //dummy
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
  fclose(eos_d5);
  fclose(eos_d6);
  fclose(eos_d7);
  fclose(eos_T1);
  fclose(eos_T2);
  fclose(eos_T3);
  fclose(eos_T4);
  fclose(eos_T5);
  fclose(eos_T6);
  fclose(eos_T7);

  cout << "Done reading EOS." << endl;

  build_velocity_of_sound_sq_matrix();

  //Allocate and fill arrays for get_s2e()
  //First, create new arrays that will permit to extract information we need from each EOS table in a unified way
  const int number_of_EOS_tables = 7;
  double lowest_eps_list[number_of_EOS_tables] = {EPP1,EPP2,EPP3,EPP4,EPP5,EPP6,EPP7};
  double delta_eps_list[number_of_EOS_tables] = {deltaEPP1,deltaEPP2,deltaEPP3,deltaEPP4,deltaEPP5,deltaEPP6,deltaEPP7};
  int nb_elements_list[number_of_EOS_tables] = {NEPP1,NEPP2,NEPP3,NEPP4,NEPP5,NEPP6,NEPP7};
  double ** s_list[number_of_EOS_tables] = {entropyDensity1,entropyDensity2,entropyDensity3,entropyDensity4,entropyDensity5,entropyDensity6,entropyDensity7};
  //Compute the maximum number of entropyDensity elements (the actual number of entropyDensity elements we will use 
  //will be smaller since we skip duplicate entropyDensities)
  int s_list_rho0_maximal_length=0;
  for(int kk=0;kk<number_of_EOS_tables;kk++) s_list_rho0_maximal_length+=nb_elements_list[kk];
  double last=0.0, tmp_s;
  s_list_rho0_length=0;
  //Allocate memory
  eps_list_rho0= new double [s_list_rho0_maximal_length];
  s_list_rho0= new double [s_list_rho0_maximal_length];

  for(int table_id=0; table_id<number_of_EOS_tables;table_id++) {
	  for(i=0;i<nb_elements_list[table_id];i++) {
		tmp_s=s_list[table_id][0][i];
		if (tmp_s > last) {
			eps_list_rho0[s_list_rho0_length]=(lowest_eps_list[table_id]+i*delta_eps_list[table_id])/hbarc;
			s_list_rho0[s_list_rho0_length]=tmp_s;
			s_list_rho0_length++;
		}
		last=tmp_s;
	  }
  }
  //Initialise the gsl interpolator
  interp_s2e = gsl_interp_alloc(gsl_interp_cspline, s_list_rho0_length);
  gsl_interp_init(interp_s2e, s_list_rho0, eps_list_rho0, s_list_rho0_length);
  accel_s2e = gsl_interp_accel_alloc();
}

void EOS::init_eos3(int selector)
{
  // read the lattice EOS pressure, temperature, and 
  // baryon chemical potential from file
  cout << "reading EOS..." << endl;
  whichEOS = 3; 
  int i, j;
  string envPath;
  if (getenv("HYDROPROGRAMPATH") != 0) {
    envPath=getenv("HYDROPROGRAMPATH");
  }
  else {
	envPath="";
  }
  
  stringstream spath;
  stringstream slocalpath;
  stringstream streos_d1_name;
  stringstream streos_d2_name;
  stringstream streos_d3_name;
  stringstream streos_d4_name;
  stringstream streos_d5_name;
  stringstream streos_d6_name;
  stringstream streos_d7_name;
  stringstream streos_T1_name;
  stringstream streos_T2_name;
  stringstream streos_T3_name;
  stringstream streos_T4_name;
  stringstream streos_T5_name;
  stringstream streos_T6_name;
  stringstream streos_T7_name;
  string eos_d1_name;
  string eos_d2_name;
  string eos_d3_name;
  string eos_d4_name;
  string eos_d5_name;
  string eos_d6_name;
  string eos_d7_name;
  string eos_T1_name;
  string eos_T2_name;
  string eos_T3_name;
  string eos_T4_name;
  string eos_T5_name;
  string eos_T6_name;
  string eos_T7_name;
  double eps, baryonDensity; //dummies for now

  spath << envPath;
  if(selector == 1)
    {
      spath << "/EOS/s95p-PCE-v1/";
      slocalpath << "./EOS/s95p-PCE-v1/";
    }
  else if(selector ==2)
    {
      spath << "/EOS/s95p-PCE155/";
      slocalpath << "./EOS/s95p-PCE155/";
     }  
  else if(selector ==3)
    {
      spath << "/EOS/s95p-PCE160/";
      slocalpath << "./EOS/s95p-PCE160/";
    }
  else if(selector ==4)
    {
      spath << "/EOS/s95p-PCE165-v0/";
      slocalpath << "./EOS/s95p-PCE165-v0/";
    }

  string path = spath.str();
  string localpath = slocalpath.str();


  if (envPath != "") // if path is set in the environment
    {
      streos_d1_name << path << "dens1.dat";
      eos_d1_name = streos_d1_name.str();
      streos_d2_name << path << "dens2.dat";
      eos_d2_name = streos_d2_name.str();
      streos_d3_name << path << "dens3.dat";
      eos_d3_name = streos_d3_name.str();
      streos_d4_name << path << "dens4.dat";
      eos_d4_name = streos_d4_name.str();
      streos_d5_name << path << "dens5.dat";
      eos_d5_name = streos_d5_name.str();
      streos_d6_name << path << "dens6.dat";
      eos_d6_name = streos_d6_name.str();
      streos_d7_name << path << "dens7.dat";
      eos_d7_name = streos_d7_name.str();

      streos_T1_name << path << "par1.dat";
      eos_T1_name = streos_T1_name.str();
      streos_T2_name << path << "par2.dat";
      eos_T2_name = streos_T2_name.str();
      streos_T3_name << path << "par3.dat";
      eos_T3_name = streos_T3_name.str();
      streos_T4_name << path << "par4.dat";
      eos_T4_name = streos_T4_name.str();
      streos_T5_name << path << "par5.dat";
      eos_T5_name = streos_T5_name.str();
      streos_T6_name << path << "par6.dat";
      eos_T6_name = streos_T6_name.str();
      streos_T7_name << path << "par7.dat";
      eos_T7_name = streos_T7_name.str();
      
      cout << "from path " << path << endl;
    }
  else //if path is not set in the environment use current folder and then /EOS/s95p-PCE160 subfolder
    {
      streos_d1_name << localpath << "dens1.dat";
      eos_d1_name = streos_d1_name.str();
      streos_d2_name << localpath << "dens2.dat";
      eos_d2_name = streos_d2_name.str();
      streos_d3_name << localpath << "dens3.dat";
      eos_d3_name = streos_d3_name.str();
      streos_d4_name << localpath << "dens4.dat";
      eos_d4_name = streos_d4_name.str();
      streos_d5_name << localpath << "dens5.dat";
      eos_d5_name = streos_d5_name.str();
      streos_d6_name << localpath << "dens6.dat";
      eos_d6_name = streos_d6_name.str();
      streos_d7_name << localpath << "dens7.dat";
      eos_d7_name = streos_d7_name.str();

      streos_T1_name << localpath << "par1.dat";
      eos_T1_name = streos_T1_name.str();
      streos_T2_name << localpath << "par2.dat";
      eos_T2_name = streos_T2_name.str();
      streos_T3_name << localpath << "par3.dat";
      eos_T3_name = streos_T3_name.str();
      streos_T4_name << localpath << "par4.dat";
      eos_T4_name = streos_T4_name.str();
      streos_T5_name << localpath << "par5.dat";
      eos_T5_name = streos_T5_name.str();
      streos_T6_name << localpath << "par6.dat";
      eos_T6_name = streos_T6_name.str();
      streos_T7_name << localpath << "par7.dat";
      eos_T7_name = streos_T7_name.str();
       
      cout << "from path " << localpath << endl;
    
    }
 

  ifstream eos_d1(eos_d1_name.c_str());
  ifstream eos_d2(eos_d2_name.c_str());
  ifstream eos_d3(eos_d3_name.c_str());
  ifstream eos_d4(eos_d4_name.c_str());
  ifstream eos_d5(eos_d5_name.c_str());
  ifstream eos_d6(eos_d6_name.c_str());
  ifstream eos_d7(eos_d7_name.c_str());
  ifstream eos_T1(eos_T1_name.c_str());
  ifstream eos_T2(eos_T2_name.c_str());
  ifstream eos_T3(eos_T3_name.c_str());
  ifstream eos_T4(eos_T4_name.c_str());
  ifstream eos_T5(eos_T5_name.c_str());
  ifstream eos_T6(eos_T6_name.c_str());
  ifstream eos_T7(eos_T7_name.c_str());

  //checkForReadError(eos_d1,eos_d1_name);
  //checkForReadError(eos_d2,eos_d2_name);
  //checkForReadError(eos_d3,eos_d3_name);
  //checkForReadError(eos_d4,eos_d4_name);
  //checkForReadError(eos_T1,eos_T1_name);
  //checkForReadError(eos_T2,eos_T2_name);
  //checkForReadError(eos_T3,eos_T3_name);
  //checkForReadError(eos_T4,eos_T4_name);
  
  //read the first two lines with general info:
  // lowest value of epsilon
  // deltaEpsilon, number of epsilon steps (i.e. # of lines)
  eos_T1 >> EPP1;
  eos_T1 >> deltaEPP1;
  eos_T1 >> NEPP1;
  eos_T2 >> EPP2;
  eos_T2 >> deltaEPP2;
  eos_T2 >> NEPP2;
  eos_T3 >> EPP3;
  eos_T3 >> deltaEPP3;
  eos_T3 >> NEPP3;
  eos_T4 >> EPP4;
  eos_T4 >> deltaEPP4;
  eos_T4 >> NEPP4;
  eos_T5 >> EPP5;
  eos_T5 >> deltaEPP5;
  eos_T5 >> NEPP5;
  eos_T6 >> EPP6;
  eos_T6 >> deltaEPP6;
  eos_T6 >> NEPP6;
  eos_T7 >> EPP7;
  eos_T7 >> deltaEPP7;
  eos_T7 >> NEPP7;
  eos_d1 >> EPP1;
  eos_d1 >> deltaEPP1;
  eos_d1 >> NEPP1;
  eos_d2 >> EPP2;
  eos_d2 >> deltaEPP2;
  eos_d2 >> NEPP2;
  eos_d3 >> EPP3;
  eos_d3 >> deltaEPP3;
  eos_d3 >> NEPP3;
  eos_d4 >> EPP4;
  eos_d4 >> deltaEPP4;
  eos_d4 >> NEPP4;
  eos_d5 >> EPP5;
  eos_d5 >> deltaEPP5;
  eos_d5 >> NEPP5;
  eos_d6 >> EPP6;
  eos_d6 >> deltaEPP6;
  eos_d6 >> NEPP6;
  eos_d7 >> EPP7;
  eos_d7 >> deltaEPP7;
  eos_d7 >> NEPP7;
  
  // fscanf(eos_d1,"%lf",&EPP1);
  // fscanf(eos_d1,"%lf %d",&deltaEPP1,&NEPP1);
  // cout << "check 4" << endl;

  // no rho_b dependence at the moment
  NBNP1=0; 
  NBNP2=0;
  NBNP3=0;
  NBNP4=0;
  NBNP5=0;
  NBNP6=0;
  NBNP7=0;

  // allocate memory for pressure arrays
  pressure1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  pressure2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  pressure3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  pressure4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  pressure5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  pressure6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  pressure7=util->mtx_malloc(NBNP7+1,NEPP7+1);
 
  // allocate memory for entropy density arrays
  entropyDensity1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  entropyDensity2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  entropyDensity3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  entropyDensity4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  entropyDensity5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  entropyDensity6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  entropyDensity7=util->mtx_malloc(NBNP7+1,NEPP7+1);

  // allocate memory for QGP fraction arrays
  QGPfraction1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  QGPfraction2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  QGPfraction3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  QGPfraction4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  QGPfraction5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  QGPfraction6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  QGPfraction7=util->mtx_malloc(NBNP7+1,NEPP7+1);

  // allocate memory for temperature arrays
  temperature1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  temperature2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  temperature3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  temperature4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  temperature5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  temperature6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  temperature7=util->mtx_malloc(NBNP7+1,NEPP7+1);
  
  // allocate memory for velocity of sound squared arrays
  cs2_1 = util->mtx_malloc(NBNP1+1, NEPP1+1);
  cs2_2 = util->mtx_malloc(NBNP2+1, NEPP2+1);
  cs2_3 = util->mtx_malloc(NBNP3+1, NEPP3+1);
  cs2_4 = util->mtx_malloc(NBNP4+1, NEPP4+1);
  cs2_5 = util->mtx_malloc(NBNP5+1, NEPP5+1);
  cs2_6 = util->mtx_malloc(NBNP6+1, NEPP6+1);
  cs2_7 = util->mtx_malloc(NBNP7+1, NEPP7+1);

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
      eos_d1 >> eps;
      eos_d1 >> pressure1[i][j];
      eos_d1 >> entropyDensity1[i][j];
      eos_d1 >> baryonDensity;
      eos_d1 >> QGPfraction1[i][j];
      //      fscanf(eos_d1,"%lf",&eps);
      //fscanf(eos_d1,"%lf",&pressure1[i][j]);
      //fscanf(eos_d1,"%lf",&entropyDensity1[i][j]);
      //fscanf(eos_d1,"%lf",&baryonDensity);
      //fscanf(eos_d1,"%lf",&QGPfraction1[i][j]);
      eos_T1 >> temperature1[i][j];
      eos_T1 >> eps;
      eos_T1 >> eps;
      //fscanf(eos_T1,"%lf",&temperature1[i][j]);
      //fscanf(eos_T1,"%lf",&eps); //dummy
      //fscanf(eos_T1,"%lf",&eps); //dummy
    }

  for(j=NEPP2-1; j>=0; j--)
    {
      eos_d2 >> eps;
      eos_d2 >> pressure2[i][j];
      eos_d2 >> entropyDensity2[i][j];
      eos_d2 >> baryonDensity;
      eos_d2 >> QGPfraction2[i][j];
      eos_T2 >> temperature2[i][j];
      eos_T2 >> eps;
      eos_T2 >> eps;
    }
  for(j=NEPP3-1; j>=0; j--)
    {
      eos_d3 >> eps;
      eos_d3 >> pressure3[i][j];
      eos_d3 >> entropyDensity3[i][j];
      eos_d3 >> baryonDensity;
      eos_d3 >> QGPfraction3[i][j];
      eos_T3 >> temperature3[i][j];
      eos_T3 >> eps;
      eos_T3 >> eps;
    }
  for(j=NEPP4-1; j>=0; j--)
    {
      eos_d4 >> eps;
      eos_d4 >> pressure4[i][j];
      eos_d4 >> entropyDensity4[i][j];
      eos_d4 >> baryonDensity;
      eos_d4 >> QGPfraction4[i][j];
      eos_T4 >> temperature4[i][j];
      eos_T4 >> eps;
      eos_T4 >> eps;
    }
  for(j=NEPP5-1; j>=0; j--)
    {
      eos_d5 >> eps;
      eos_d5 >> pressure5[i][j];
      eos_d5 >> entropyDensity5[i][j];
      eos_d5 >> baryonDensity;
      eos_d5 >> QGPfraction5[i][j];
      eos_T5 >> temperature5[i][j];
      eos_T5 >> eps;
      eos_T5 >> eps;
    }
  for(j=NEPP6-1; j>=0; j--)
    {
      eos_d6 >> eps;
      eos_d6 >> pressure6[i][j];
      eos_d6 >> entropyDensity6[i][j];
      eos_d6 >> baryonDensity;
      eos_d6 >> QGPfraction6[i][j];
      eos_T6 >> temperature6[i][j];
      eos_T6 >> eps;
      eos_T6 >> eps;
    }
  for(j=NEPP7-1; j>=0; j--)
    {
      eos_d7 >> eps;
      eos_d7 >> pressure7[i][j];
      eos_d7 >> entropyDensity7[i][j];
      eos_d7 >> baryonDensity;
      eos_d7 >> QGPfraction7[i][j];
      eos_T7 >> temperature7[i][j];
      eos_T7 >> eps;
      eos_T7 >> eps;
    }

 //  //test if the reading worked:
//   for(j=0; j<NEPP1; j++) 
//     { 
//       fprintf(stderr,"%lf %lf %lf %lf\n",pressure1[i][j],entropyDensity1[i][j],QGPfraction1[i][j],temperature1[i][j]); 
//     } 
//   for(j=0; j<NEPP2; j++) 
//     { 
//       fprintf(stderr,"%lf %lf %lf %lf\n",pressure2[i][j],entropyDensity2[i][j],QGPfraction2[i][j],temperature2[i][j]); 
//     } 
//   for(j=0; j<NEPP3; j++) 
//     { 
//       fprintf(stderr,"%lf %lf %lf %lf\n",pressure3[i][j],entropyDensity3[i][j],QGPfraction3[i][j],temperature3[i][j]); 
//     } 
//   cout << " check before 4" << endl;
//   for(j=0; j<NEPP4; j++) 
//     { 
//       fprintf(stderr,"%lf %lf %lf %lf\n",pressure4[i][j],entropyDensity4[i][j],QGPfraction4[i][j],temperature4[i][j]); 
//     } 
//   cout << " check before 5" << endl;
//   for(j=0; j<NEPP5; j++) 
//     { 
//       fprintf(stderr,"%lf %lf %lf %lf\n",pressure5[i][j],entropyDensity5[i][j],QGPfraction5[i][j],temperature5[i][j]); 
//     } 
//   cout << " check before 6" << endl;
//   for(j=0; j<NEPP6; j++) 
//     { 
//       fprintf(stderr,"%lf %lf %lf %lf\n",pressure6[i][j],entropyDensity6[i][j],QGPfraction6[i][j],temperature6[i][j]); 
//     } 
//   cout << " check before 7" << endl;
//   for(j=0; j<NEPP7; j++) 
//     { 
//       fprintf(stderr,"%lf %lf %lf %lf\n",pressure7[i][j],entropyDensity7[i][j],QGPfraction7[i][j],temperature7[i][j]); 
//     } 

  
  eos_d1.close();
  eos_d2.close();
  eos_d3.close();
  eos_d4.close();
  eos_d5.close();
  eos_d6.close();
  eos_d7.close();
  eos_T1.close();
  eos_T2.close();
  eos_T3.close();
  eos_T4.close();
  eos_T5.close();
  eos_T6.close();
  eos_T7.close();

  cout << "Done reading EOS." << endl;

  build_velocity_of_sound_sq_matrix();

  //Allocate and fill arrays for get_s2e()
  //First, create new arrays that will permit to extract information we need from each EOS table in a unified way
  const int number_of_EOS_tables = 7;
  double lowest_eps_list[number_of_EOS_tables] = {EPP1,EPP2,EPP3,EPP4,EPP5,EPP6,EPP7};
  double delta_eps_list[number_of_EOS_tables] = {deltaEPP1,deltaEPP2,deltaEPP3,deltaEPP4,deltaEPP5,deltaEPP6,deltaEPP7};
  int nb_elements_list[number_of_EOS_tables] = {NEPP1,NEPP2,NEPP3,NEPP4,NEPP5,NEPP6,NEPP7};
  double ** s_list[number_of_EOS_tables] = {entropyDensity1,entropyDensity2,entropyDensity3,entropyDensity4,entropyDensity5,entropyDensity6,entropyDensity7};
  //Compute the maximum number of entropyDensity elements (the actual number of entropyDensity elements we will use 
  //will be smaller since we skip duplicate entropyDensities)
  int s_list_rho0_maximal_length=0;
  for(int kk=0;kk<number_of_EOS_tables;kk++) s_list_rho0_maximal_length+=nb_elements_list[kk];
  double last=0.0, tmp_s;
  s_list_rho0_length=0;
  //Allocate memory
  eps_list_rho0= new double [s_list_rho0_maximal_length];
  s_list_rho0= new double [s_list_rho0_maximal_length];

  for(int table_id=0; table_id<number_of_EOS_tables;table_id++) {
	  for(i=0;i<nb_elements_list[table_id];i++) {
		tmp_s=s_list[table_id][0][i];
		if (tmp_s > last) {
			eps_list_rho0[s_list_rho0_length]=(lowest_eps_list[table_id]+i*delta_eps_list[table_id])/hbarc;
			s_list_rho0[s_list_rho0_length]=tmp_s;
			s_list_rho0_length++;
		}
		last=tmp_s;
	  }
  }
  //Initialise the gsl interpolator
  interp_s2e = gsl_interp_alloc(gsl_interp_cspline, s_list_rho0_length);
  gsl_interp_init(interp_s2e, s_list_rho0, eps_list_rho0, s_list_rho0_length);
  accel_s2e = gsl_interp_accel_alloc();
}

void EOS::init_eos10(int selector)
// read the lattice EOS at finite muB
// pressure, temperature, and baryon chemical potential from file
{
  cout << "reading EOS..." << endl;
  whichEOS = 10; 

  // get environment path
  string envPath;
  if (getenv("HYDROPROGRAMPATH") != 0)
    envPath=getenv("HYDROPROGRAMPATH");
  else 
    envPath="";
  
  stringstream spath;
  stringstream slocalpath;

  if(selector == 0)  // lattice EOS from A. Monnai
  {
    spath << envPath << "/EOS/neos_2/";
    slocalpath << "./EOS/neos_2/";
  }
  else
  {
    fprintf(stderr, "EOS::init_eos10: unrecognized selector = %d \n", selector);
    exit(-1);
  }

  string path;
  if (envPath != "") // if path is set in the environment
      path = spath.str();
  else
      path = slocalpath.str();
  
  stringstream streos_p1_name;
  stringstream streos_p2_name;
  stringstream streos_p3_name;
  stringstream streos_p4_name;
  stringstream streos_p5_name;
  stringstream streos_T1_name;
  stringstream streos_T2_name;
  stringstream streos_T3_name;
  stringstream streos_T4_name;
  stringstream streos_T5_name;
  stringstream streos_mub1_name;
  stringstream streos_mub2_name;
  stringstream streos_mub3_name;
  stringstream streos_mub4_name;
  stringstream streos_mub5_name;

  streos_p1_name << path << "neos0_p.dat";
  streos_p2_name << path << "neos1a_p.dat";
  streos_p3_name << path << "neos2_p.dat";
  streos_p4_name << path << "neos3_p.dat";
  streos_p5_name << path << "neos4_p.dat";

  streos_T1_name << path << "neos0_t.dat";
  streos_T2_name << path << "neos1a_t.dat";
  streos_T3_name << path << "neos2_t.dat";
  streos_T4_name << path << "neos3_t.dat";
  streos_T5_name << path << "neos4_t.dat";
  
  streos_mub1_name << path << "neos0_mb.dat";
  streos_mub2_name << path << "neos1a_mb.dat";
  streos_mub3_name << path << "neos2_mb.dat";
  streos_mub4_name << path << "neos3_mb.dat";
  streos_mub5_name << path << "neos4_mb.dat";
  
  cout << "from path " << path << endl;

  ifstream eos_p1(streos_p1_name.str().c_str());
  ifstream eos_p2(streos_p2_name.str().c_str());
  ifstream eos_p3(streos_p3_name.str().c_str());
  ifstream eos_p4(streos_p4_name.str().c_str());
  ifstream eos_p5(streos_p5_name.str().c_str());
  ifstream eos_T1(streos_T1_name.str().c_str());
  ifstream eos_T2(streos_T2_name.str().c_str());
  ifstream eos_T3(streos_T3_name.str().c_str());
  ifstream eos_T4(streos_T4_name.str().c_str());
  ifstream eos_T5(streos_T5_name.str().c_str());
  ifstream eos_mub1(streos_mub1_name.str().c_str());
  ifstream eos_mub2(streos_mub2_name.str().c_str());
  ifstream eos_mub3(streos_mub3_name.str().c_str());
  ifstream eos_mub4(streos_mub4_name.str().c_str());
  ifstream eos_mub5(streos_mub5_name.str().c_str());
  
  // read the first two lines with general info:
  // first value of rhob, first value of epsilon
  // deltaRhob, deltaEpsilon, number of rhob steps, number of epsilon steps
  eos_p1 >> BNP1 >> EPP1;
  eos_p1 >> deltaBNP1 >> deltaEPP1 >> NBNP1 >> NEPP1;
  eos_p2 >> BNP2 >> EPP2;
  eos_p2 >> deltaBNP2 >> deltaEPP2 >> NBNP2 >> NEPP2;
  eos_p3 >> BNP3 >> EPP3;
  eos_p3 >> deltaBNP3 >> deltaEPP3 >> NBNP3 >> NEPP3;
  eos_p4 >> BNP4 >> EPP4;
  eos_p4 >> deltaBNP4 >> deltaEPP4 >> NBNP4 >> NEPP4;
  eos_p5 >> BNP5 >> EPP5;
  eos_p5 >> deltaBNP5 >> deltaEPP5 >> NBNP5 >> NEPP5;
  eos_T1 >> BNP1 >> EPP1;
  eos_T1 >> deltaBNP1 >> deltaEPP1 >> NBNP1 >> NEPP1;
  eos_T2 >> BNP2 >> EPP2;
  eos_T2 >> deltaBNP2 >> deltaEPP2 >> NBNP2 >> NEPP2;
  eos_T3 >> BNP3 >> EPP3;
  eos_T3 >> deltaBNP3 >> deltaEPP3 >> NBNP3 >> NEPP3;
  eos_T4 >> BNP4 >> EPP4;
  eos_T4 >> deltaBNP4 >> deltaEPP4 >> NBNP4 >> NEPP4;
  eos_T5 >> BNP5 >> EPP5;
  eos_T5 >> deltaBNP5 >> deltaEPP5 >> NBNP5 >> NEPP5;
  eos_mub1 >> BNP1 >> EPP1;
  eos_mub1 >> deltaBNP1 >> deltaEPP1 >> NBNP1 >> NEPP1;
  eos_mub2 >> BNP2 >> EPP2;
  eos_mub2 >> deltaBNP2 >> deltaEPP2 >> NBNP2 >> NEPP2;
  eos_mub3 >> BNP3 >> EPP3;
  eos_mub3 >> deltaBNP3 >> deltaEPP3 >> NBNP3 >> NEPP3;
  eos_mub4 >> BNP4 >> EPP4;
  eos_mub4 >> deltaBNP4 >> deltaEPP4 >> NBNP4 >> NEPP4;
  eos_mub5 >> BNP5 >> EPP5;
  eos_mub5 >> deltaBNP5 >> deltaEPP5 >> NBNP5 >> NEPP5;

  EPP6 = 1e4;  // take a large enough value to make sure only the first 5 tables will be used

  // allocate memory for pressure arrays
  pressure1=util->mtx_malloc(NBNP1+1, NEPP1+1);
  pressure2=util->mtx_malloc(NBNP2+1, NEPP2+1);
  pressure3=util->mtx_malloc(NBNP3+1, NEPP3+1);
  pressure4=util->mtx_malloc(NBNP4+1, NEPP4+1);
  pressure5=util->mtx_malloc(NBNP5+1, NEPP5+1);
 
  // allocate memory for entropy density arrays
  entropyDensity1=util->mtx_malloc(NBNP1+1, NEPP1+1);
  entropyDensity2=util->mtx_malloc(NBNP2+1, NEPP2+1);
  entropyDensity3=util->mtx_malloc(NBNP3+1, NEPP3+1);
  entropyDensity4=util->mtx_malloc(NBNP4+1, NEPP4+1);
  entropyDensity5=util->mtx_malloc(NBNP5+1, NEPP5+1);

  // allocate memory for temperature arrays
  temperature1=util->mtx_malloc(NBNP1+1, NEPP1+1);
  temperature2=util->mtx_malloc(NBNP2+1, NEPP2+1);
  temperature3=util->mtx_malloc(NBNP3+1, NEPP3+1);
  temperature4=util->mtx_malloc(NBNP4+1, NEPP4+1);
  temperature5=util->mtx_malloc(NBNP5+1, NEPP5+1);
  
  // allocate memory for mu_B arrays
  mu1=util->mtx_malloc(NBNP1+1, NEPP1+1);
  mu2=util->mtx_malloc(NBNP2+1, NEPP2+1);
  mu3=util->mtx_malloc(NBNP3+1, NEPP3+1);
  mu4=util->mtx_malloc(NBNP4+1, NEPP4+1);
  mu5=util->mtx_malloc(NBNP5+1, NEPP5+1);
  
  // allocate memory for velocity of sound squared arrays
  cs2_1 = util->mtx_malloc(NBNP1+1, NEPP1+1);
  cs2_2 = util->mtx_malloc(NBNP2+1, NEPP2+1);
  cs2_3 = util->mtx_malloc(NBNP3+1, NEPP3+1);
  cs2_4 = util->mtx_malloc(NBNP4+1, NEPP4+1);
  cs2_5 = util->mtx_malloc(NBNP5+1, NEPP5+1);

  // read pressure, temperature and chemical potential values
  for(int j = 0; j < NEPP1 + 1; j++)
  {
    double ed = EPP1 + j*deltaEPP1;
    for(int i = 0; i < NBNP1 + 1; i++)
    {
      double rhob = BNP1 + i*deltaBNP1;
      eos_p1 >> pressure1[i][j];
      eos_T1 >> temperature1[i][j];
      eos_mub1 >> mu1[i][j];
      double sd = (ed + pressure1[i][j] - mu1[i][j]*rhob)/temperature1[i][j];
      entropyDensity1[i][j] = sd;
    }
  }
  
  for(int j = 0; j < NEPP2 + 1; j++)
  {
    double ed = EPP2 + j*deltaEPP2;
    for(int i = 0; i < NBNP2 + 1; i++)
    {
      double rhob = BNP2 + i*deltaBNP2;
      eos_p2 >> pressure2[i][j];
      eos_T2 >> temperature2[i][j];
      eos_mub2 >> mu2[i][j];
      double sd = (ed + pressure2[i][j] - mu2[i][j]*rhob)/temperature2[i][j];
      entropyDensity2[i][j] = sd;
    }
  }
  
  for(int j = 0; j < NEPP3 + 1; j++)
  {
    double ed = EPP3 + j*deltaEPP3;
    for(int i = 0; i < NBNP3 + 1; i++)
    {
      double rhob = BNP3 + i*deltaBNP3;
      eos_p3 >> pressure3[i][j];
      eos_T3 >> temperature3[i][j];
      eos_mub3 >> mu3[i][j];
      double sd = (ed + pressure3[i][j] - mu3[i][j]*rhob)/temperature3[i][j];
      entropyDensity3[i][j] = sd;
    }
  }
  
  for(int j = 0; j < NEPP4 + 1; j++)
  {
    double ed = EPP4 + j*deltaEPP4;
    for(int i = 0; i < NBNP4 + 1; i++)
    {
      double rhob = BNP4 + i*deltaBNP4;
      eos_p4 >> pressure4[i][j];
      eos_T4 >> temperature4[i][j];
      eos_mub4 >> mu4[i][j];
      double sd = (ed + pressure4[i][j] - mu4[i][j]*rhob)/temperature4[i][j];
      entropyDensity4[i][j] = sd;
    }
  }
  
  for(int j = 0; j < NEPP5 + 1; j++)
  {
    double ed = EPP5 + j*deltaEPP5;
    for(int i = 0; i < NBNP5 + 1; i++)
    {
      double rhob = BNP5 + i*deltaBNP5;
      eos_p5 >> pressure5[i][j];
      eos_T5 >> temperature5[i][j];
      eos_mub5 >> mu5[i][j];
      double sd = (ed + pressure5[i][j] - mu5[i][j]*rhob)/temperature5[i][j];
      entropyDensity5[i][j] = sd;
    }
  }

  eos_p1.close();
  eos_p2.close();
  eos_p3.close();
  eos_p4.close();
  eos_p5.close();
  eos_T1.close();
  eos_T2.close();
  eos_T3.close();
  eos_T4.close();
  eos_T5.close();
  eos_mub1.close();
  eos_mub2.close();
  eos_mub3.close();
  eos_mub4.close();
  eos_mub5.close();

  cout << "Done reading EOS." << endl;

  build_velocity_of_sound_sq_matrix();

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
	  fprintf(stderr,"ie1=%d,BNP1=%lf; inb1=%d, NBNP1=%d, e=%f, rhob=%f \n", ie1, BNP1, inb1, NBNP1,e,rhob);
	 
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
	  fprintf(stderr,"ie1=%d,BNP1=%lf; inb1=%d, NBNP1=%d, e=%f, rhob=%f \n", ie1, BNP1, inb1, NBNP1,e,rhob);
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
  double p, pa, pb;

  //use linear interpolation
  int ie1, ie2, NEps;
  double frace;
  double eps0, deltaEps;
  double **array; 
  //selector = 0 : pressure
  //selector = 1 : temperature
  //selector = 2 : entropy density
  //selector = 3 : QGP fraction
  //selector = 4 : velocity of sound squared
  
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
	  fprintf(stderr,"ERROR interpolate2: someting's wrong.\n");
	  fprintf(stderr,"ie1=%d,NEPP1=%d\n", ie1, NEPP1);
	  exit(0);
	}
      if(ie2>NEPP1)
	{
	  fprintf(stderr,"ERROR interpolate2: someting's wrong.\n");
	  fprintf(stderr,"ie2=%d,NEPP1=%d\n", ie2, NEPP1);
	  exit(0);
	}

      switch (selector) 
	{
	case 0: array = pressure1; break;
	case 1: array = temperature1; break;
	case 2: array = entropyDensity1; break;
	case 3: array = QGPfraction1; break;
	case 4: array = cs2_1; break;
	default: fprintf(stderr,"ERROR in interpolate2 - selector must be 0,1,2,3, or 4\n"); exit(1);
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
	    case 4: array = cs2_2; break;
	    default: fprintf(stderr,"ERROR in interpolate2 - selector must be 0,1,2,3, or 4\n"); exit(1);
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
	    case 4: array = cs2_3; break;
	    default: fprintf(stderr,"ERROR in interpolate2 - selector must be 0,1,2,3, or 4\n"); exit(1);
	    }
	}
      else if (e<EPP5)
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
	    case 4: array = cs2_4; break;
	    default: fprintf(stderr,"ERROR in interpolate2 - selector must be 0,1,2,3, or 4\n"); exit(1);
	    }
	}
      else if (e<EPP6)
	{
	  eps0 = EPP5;
	  NEps = NEPP5;
	  deltaEps = deltaEPP5;
	  switch (selector) 
	    {
	    case 0: array = pressure5; break;
	    case 1: array = temperature5; break;
	    case 2: array = entropyDensity5; break;
	    case 3: array = QGPfraction5; break;
	    case 4: array = cs2_5; break;
	    default: fprintf(stderr,"ERROR in interpolate2 - selector must be 0,1,2,3, or 4\n"); exit(1);
	    }
	}
      else if (e<EPP7)
	{
	  eps0 = EPP6;
	  NEps = NEPP6;
	  deltaEps = deltaEPP6;
	  switch (selector) 
	    {
	    case 0: array = pressure6; break;
	    case 1: array = temperature6; break;
	    case 2: array = entropyDensity6; break;
	    case 3: array = QGPfraction6; break;
	    case 4: array = cs2_6; break;
	    default: fprintf(stderr,"ERROR in interpolate2 - selector must be 0,1,2,3, or 4\n"); exit(1);
	    }
	}
      else
	{
	  eps0 = EPP7;
	  NEps = NEPP7;
	  deltaEps = deltaEPP7;
	  switch (selector) 
	    {
	    case 0: array = pressure7; break;
	    case 1: array = temperature7; break;
	    case 2: array = entropyDensity7; break;
	    case 3: array = QGPfraction7; break;
	    case 4: array = cs2_7; break;
	    default: fprintf(stderr,"ERROR in interpolate2 - selector must be 0,1,2,3, or 4\n"); exit(1);
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
      if(ie2>=NEps)
	{
	  fprintf(stderr,"ERROR in inperpolate2. out of range.\n");
	  fprintf(stderr,"ie2=%d,NEPP2=%d\n", ie2, NEps);
	  exit(0);
	}

      //if (ie1 ==NEps)
	//cout << " last e=" << e << " T=" << T1*pow(e,T2)/hbarc << endl;

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
    case 4: break;
    default: fprintf(stderr,"ERROR in interpolate2 - selector must be 0,1,2,3, or 4\n"); exit(1);
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

      pL = get_pressure(eLeft/hbarc, rhob);
      pR = get_pressure(eRight/hbarc,  rhob);
      
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
      
      pL = get_pressure(eLeft/hbarc, rhob);
      pR = get_pressure(eRight/hbarc, rhob);
      
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
      pL = get_pressure(e/hbarc, rhoLeft);
      pR = get_pressure(e/hbarc, rhoRight);

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

      pL = get_pressure(e/hbarc, rhoLeft); 
      pR = get_pressure(e/hbarc, rhoRight);

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

   //The energy has to be in GeV in what follows
   e=e*hbarc;

   double eps0, deltaEps;
   int Neps;

   if(e < EPP1)
   {
       eps0 = 0.0;
       deltaEps = EPP1;
       Neps = 2;
   }
   if(e<EPP2)
   {
       eps0 = EPP1;
       deltaEps = deltaEPP1;
       Neps = NEPP1;
   }
   else if (e<EPP3)
   {
       eps0 = EPP2;
       deltaEps = deltaEPP2;
       Neps = NEPP2;
   }
   else if (e<EPP4)
   {
       eps0 = EPP3;
       deltaEps = deltaEPP3;
       Neps = NEPP3;
   }
   else if (e<EPP5)
   {
       eps0 = EPP4;
       deltaEps = deltaEPP4;
       Neps = NEPP4;
   }
   else if (e<EPP6)
   {
       eps0 = EPP5;
       deltaEps = deltaEPP5;
       Neps = NEPP5;
   }
   else if (e<EPP7)
   {
       eps0 = EPP6;
       deltaEps = deltaEPP6;
       Neps = NEPP6;
   }
   else 
   {
       eps0 = EPP7;
       deltaEps = deltaEPP7;
       Neps = NEPP7;
   }
   double eps_end = eps0 + (Neps - 1)*deltaEps;

   double eLeft = e - deltaEps*0.1;    // GeV/fm^3
   double eRight = e + deltaEps*0.1;   // GeV/fm^3
  
   // deal with boundary, avoid to exceed the table
   if(eLeft < eps0)
       eLeft = eps0;
   if(eRight > eps_end)
       eRight = eps_end;

   double pL = get_pressure(eLeft/hbarc, rhob);  // 1/fm^4
   double pR = get_pressure(eRight/hbarc, rhob); // 1/fm^4
      
   double dpde = (pR - pL)*hbarc/(eRight - eLeft);
      
   if(dpde < 0) 
   { 
       fprintf(stderr, "dp/de=%lf\n", dpde); 
       fprintf(stderr, "pL=%lf\n", pL); 
       fprintf(stderr, "pR=%lf\n", pR); 
       fprintf(stderr, "e=%lf\n", e);  
       fprintf(stderr, "eLeft=%lf\n", eLeft/hbarc);  
       fprintf(stderr, "eRight=%lf\n", eRight/hbarc);  
   } 
   return dpde;
}

double EOS::get_dpOverdrhob2(double e, double rhob)
{
    double local_ed = e*hbarc;      // GeV/fm^3
    double local_rhob = rhob;       // 1/fm^3
    
    double deltaRhob;               // in 1/fm^3
    double rhob_end;
    if (local_ed < EPP2)
    {
        deltaRhob = deltaBNP1;
        rhob_end = BNP1 + (NBNP1 - 1)*deltaRhob;
    }
    else if (local_ed < EPP3)
    {
        deltaRhob = deltaBNP2;
        rhob_end = BNP2 + (NBNP2 - 1)*deltaRhob;
    }
    else if (local_ed < EPP4)
    {
        deltaRhob = deltaBNP3;
        rhob_end = BNP3 + (NBNP3 - 1)*deltaRhob;
    }
    else if (local_ed < EPP5)
    {
        deltaRhob = deltaBNP4;
        rhob_end = BNP4 + (NBNP4 - 1)*deltaRhob;
    }
    else if (local_ed < EPP6)
    {
        deltaRhob = deltaBNP5;
        rhob_end = BNP5 + (NBNP5 - 1)*deltaRhob;
    }
    else if (local_ed < EPP7)
    {
        deltaRhob = deltaBNP6;
        rhob_end = BNP6 + (NBNP6 - 1)*deltaRhob;
    }
    else 
    {
        deltaRhob = deltaBNP7;
        rhob_end = BNP7 + (NBNP7 - 1)*deltaRhob;
    }
    
    double rhobLeft = local_rhob - deltaRhob*0.1;
    double rhobRight = local_rhob + deltaRhob*0.1;

    // avoid exceeding the table
    if(rhobRight > rhob_end)
        rhobRight = rhob_end;
   
    double pL = get_pressure(e, rhobLeft);  // 1/fm^4
    double pR = get_pressure(e, rhobRight); // 1/fm^4
      
    double dpdrho = (pR - pL)/(rhobRight - rhobLeft);  // 1/fm
   
    return dpdrho;   // in 1/fm
}

double EOS::get_cs2(double e, double rhob)
{
    double f;
    if (whichEOS==0)
        f = cs2;
    else if (whichEOS==1)
        f = interpolate(e, rhob, 2);
    else if (whichEOS>=2 && whichEOS < 10)
        f = interpolate2(e, rhob, 4);
    else if (whichEOS >= 10)
        f = interpolate2D(e, fabs(rhob), 4);
    else
    {
        fprintf(stderr,"EOS::get_cs2: whichEOS = %d is out of range!\n", whichEOS);
        exit(0);
    }
    return f;
}

void EOS::build_velocity_of_sound_sq_matrix()
{
    // whiichEOS == 0 : ideal gas EoS, cs^2 = 1/3 always
    if (whichEOS==1)  // 2 tables
    {
        fill_cs2_matrix(EPP1, deltaEPP1, NEPP1, BNP1, deltaBNP1, NBNP1, cs2_1);
        fill_cs2_matrix(EPP2, deltaEPP2, NEPP2, BNP2, deltaBNP2, NBNP2, cs2_2);
    }
    else if (whichEOS>=2 && whichEOS < 10)  // 7 tables
    {
        fill_cs2_matrix(EPP1, deltaEPP1, NEPP1, BNP1, deltaBNP1, NBNP1, cs2_1);
        fill_cs2_matrix(EPP2, deltaEPP2, NEPP2, BNP2, deltaBNP2, NBNP2, cs2_2);
        fill_cs2_matrix(EPP3, deltaEPP3, NEPP3, BNP3, deltaBNP3, NBNP3, cs2_3);
        fill_cs2_matrix(EPP4, deltaEPP4, NEPP4, BNP4, deltaBNP4, NBNP4, cs2_4);
        fill_cs2_matrix(EPP5, deltaEPP5, NEPP5, BNP5, deltaBNP5, NBNP5, cs2_5);
        fill_cs2_matrix(EPP6, deltaEPP6, NEPP6, BNP6, deltaBNP6, NBNP6, cs2_6);
        fill_cs2_matrix(EPP7, deltaEPP7, NEPP7, BNP7, deltaBNP7, NBNP7, cs2_7);
    }
    else if (whichEOS >= 10)  // 5 tables
    {
        fill_cs2_matrix(EPP1, deltaEPP1, NEPP1, BNP1, deltaBNP1, NBNP1, cs2_1);
        fill_cs2_matrix(EPP2, deltaEPP2, NEPP2, BNP2, deltaBNP2, NBNP2, cs2_2);
        fill_cs2_matrix(EPP3, deltaEPP3, NEPP3, BNP3, deltaBNP3, NBNP3, cs2_3);
        fill_cs2_matrix(EPP4, deltaEPP4, NEPP4, BNP4, deltaBNP4, NBNP4, cs2_4);
        fill_cs2_matrix(EPP5, deltaEPP5, NEPP5, BNP5, deltaBNP5, NBNP5, cs2_5);
    }
    else
    {
        fprintf(stderr,"EOS::build_velocity_of_sound_sq_matrix: whichEOS = %d is out of range!\n", whichEOS);
        exit(0);
    }
}

void EOS::fill_cs2_matrix(double e0, double de, int ne, double rhob0, double drhob, int nrhob, double** cs2_ptr)
{
    for(int i = 0; i < nrhob; i++)
    {
        double rhob_local = rhob0 + i*drhob;
        for(int j = 0; j < ne; j++)
        {
            double epsilon_local = (e0 + j*de)/hbarc;
            double cs2_local = calculate_velocity_of_sound_sq(epsilon_local, rhob_local);
            cs2_ptr[i][j] = cs2_local;
        }
    }
}

double EOS::calculate_velocity_of_sound_sq(double e, double rhob)
{
    double v_min = 0.0;
    double dpde = p_e_func(e, rhob);
    double dpdrho = p_rho_func(e, rhob);
    double pressure = get_pressure(e, rhob);
    double v_sound = dpde + rhob/(e + pressure + 1e-15)*dpdrho;

    if(v_sound < v_min)
        v_sound = v_min;
    /*
    if(v_sound < 0.)
    {
        fprintf(stderr, "EOS::get_velocity_of_sound:velocity of sound is negative!\n");
        fprintf(stderr, "v_sound = %lf \n", v_sound);
        fprintf(stderr, "e = %lf \n", e*hbarc);
        fprintf(stderr, "rhob = %lf \n", rhob);
        fprintf(stderr, "pressure = %lf \n", pressure);
        fprintf(stderr, "dP/de = %lf \n", dpde);
        fprintf(stderr, "dP/drho = %lf \n", dpdrho);
    }
    */
    return(v_sound);
}

double EOS::get_pressure(double e, double rhob)
// return pressure in [1/fm^4]
{
    double f;
    if (whichEOS==0)
        f = cs2*e;
    else if (whichEOS==1)
        f = interpolate_pressure(e,rhob);
    else if (whichEOS>=2 && whichEOS < 10)
        f = interpolate2(e,rhob,0);    //selector 0 means get pressure 
    else if (whichEOS >= 10)
        f = interpolate2D(e, fabs(rhob), 0);   // EOS is symmetric in rho_b for pressure
    else
    {
        fprintf(stderr,"EOS::get_pressure: whichEOS = %d is out of range!\n", whichEOS);
        exit(0);
    }
    return f;
}/* get_pressure */


double EOS::p_rho_func(double e, double rhob)
// return dP/drho_b (in 1/fm)
{
    double f;
    if (whichEOS==0) 
        f = 0.0;
    else if (whichEOS==1)
        f = get_dpOverdrhob(e, rhob);
    else if (whichEOS>=2 && whichEOS < 10)
        f = 0.0;
    else if (whichEOS >= 10)
        f = get_dpOverdrhob2(e, rhob);
    else
    {
        fprintf(stderr,"EOS::p_rho_func: whichEOS = %d is out of range!\n", whichEOS);
        exit(0);
    }
    return f;
}/* p_rho_func */

double EOS::p_e_func(double e, double rhob)
// return dP/de
{
    double f;
    if (whichEOS==0)
        f = cs2;
    else if (whichEOS==1)
        f = get_dpOverde(e, rhob);
    else if (whichEOS>=2 && whichEOS < 10)
        f = get_dpOverde2(e, rhob);
    else if (whichEOS >= 10)
        f = get_dpOverde2(e, rhob);
    else
    {
        fprintf(stderr,"EOS::p_e_func: whichEOS = %d is out of range!\n", whichEOS);
        exit(0);
    }
    return f;
}/* p_e_func */


double EOS::interpolate(double e, double rhob, int selector)
{
  // selector=0 : return temperature
  // selector=1 : return baryon chemical potential
  // selector=2 : return cs^2
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
      else if (selector == 2)
      {
	  pa1 = cs2_1[inb1][ie1];
	  pa2 = cs2_1[inb2][ie1];
      }
   	else
      {
        fprintf(stderr,"selector = %d is out of range.\n", selector); 
        exit(0);
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
      else if (selector == 2)
      {
	  pb1 = cs2_1[inb1][ie2];
	  pb2 = cs2_1[inb2][ie2];
      }
   	else
      {
        fprintf(stderr,"selector = %d is out of range.\n", selector); 
        exit(0);
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
      else if (selector == 2)
      {
	  pa1 = cs2_1[inb1][ie1];
	  pa2 = cs2_1[inb2][ie1];
      }
   	else
      {
        fprintf(stderr,"selector = %d is out of range.\n", selector); 
        exit(0);
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
      else if (selector == 2)
      {
	  pb1 = cs2_1[inb1][ie2];
	  pb2 = cs2_1[inb2][ie2];
      }
   	else
      {
        fprintf(stderr,"selector = %d is out of range.\n", selector); 
        exit(0);
      }

      if (pb1<0) pb1=0.;
      if (pb2<0) pb2=0.;
      
      pb = pb1*(1-fracnb) + pb2*fracnb;

      T = pa*(1-frace) + pb*frace;
    }
  return T/hbarc;
}

double EOS::interpolate2D(double e, double rhob, int selector)
// This is a generic bilinear interpolation routine for EOS at finite mu_B
// it assumes the class has already read in 
//        P(e, rho_b), T(e, rho_b), s(e, rho_b), mu_b(e, rho_b) 
// as two-dimensional arrays on an equally spacing lattice grid
// units: e is in 1/fm^4, rhob is in 1/fm^3
// selector: return specific type of thermodynamical quantities
//           0: pressure [1/fm^4]
//           1: temperature [1/fm]
//           2: entropy density [1/fm^3]
//           3: mu_B [1/fm]
//           4: velocity of sound squared cs^2
{
    double **array; 
  
    double local_ed = e*hbarc; // convert energy density from 1/fm^4 to GeV/fm^3
    double local_rhob = rhob;  // [1/fm^3]

    // first choosing the right table
    double eps0, rhob0, deltaEps, deltaRhob;
    int NEps, Nrhob;
    if(local_ed < EPP1)  // energy density is smaller than the smallest value in the table use linear extrapolation
    {
	  eps0 = EPP1;
	  NEps = NEPP1;
	  deltaEps = deltaEPP1;
        rhob0 = BNP1;
	  Nrhob = NBNP1;
	  deltaRhob = deltaBNP1;
        switch (selector) 
	  {
	      case 0: array = pressure1; break;
	      case 1: array = temperature1; break;
	      case 2: array = entropyDensity1; break;
	      case 3: array = mu1; break;
	      case 4: array = cs2_1; break;
	      default: fprintf(stderr,"ERROR in interpolate2D - selector must be 0,1,2,3, or 4, selector = %d \n", selector); exit(1);
	  }
    }
    else if(local_ed < EPP2)
    {
	  eps0 = EPP1;
	  NEps = NEPP1;
	  deltaEps = deltaEPP1;
        rhob0 = BNP1;
	  Nrhob = NBNP1;
	  deltaRhob = deltaBNP1;
        switch (selector) 
	  {
	      case 0: array = pressure1; break;
	      case 1: array = temperature1; break;
	      case 2: array = entropyDensity1; break;
	      case 3: array = mu1; break;
	      case 4: array = cs2_1; break;
	      default: fprintf(stderr,"ERROR in interpolate2D - selector must be 0,1,2,3, or 4, selector = %d \n", selector); exit(1);
	  }
    }
    else if(local_ed < EPP3)
    {
	  eps0 = EPP2;
	  NEps = NEPP2;
	  deltaEps = deltaEPP2;
        rhob0 = BNP2;
	  Nrhob = NBNP2;
	  deltaRhob = deltaBNP2;
	  switch (selector) 
	  {
	      case 0: array = pressure2; break;
	      case 1: array = temperature2; break;
	      case 2: array = entropyDensity2; break;
	      case 3: array = mu2; break;
	      case 4: array = cs2_2; break;
	      default: fprintf(stderr,"ERROR in interpolate2D - selector must be 0,1,2,3, or 4, selector = %d \n", selector); exit(1);
	  }
    }
    else if(local_ed < EPP4)
    {
	  eps0 = EPP3;
	  NEps = NEPP3;
	  deltaEps = deltaEPP3;
        rhob0 = BNP3;
	  Nrhob = NBNP3;
	  deltaRhob = deltaBNP3;
	  switch (selector) 
	  {
	      case 0: array = pressure3; break;
	      case 1: array = temperature3; break;
	      case 2: array = entropyDensity3; break;
	      case 3: array = mu3; break;
	      case 4: array = cs2_3; break;
	      default: fprintf(stderr,"ERROR in interpolate2D - selector must be 0,1,2,3, or 4, selector = %d \n", selector); exit(1);
	  }
    }
    else if(local_ed < EPP5)
    {
	  eps0 = EPP4;
	  NEps = NEPP4;
	  deltaEps = deltaEPP4;
        rhob0 = BNP4;
	  Nrhob = NBNP4;
	  deltaRhob = deltaBNP4;
	  switch (selector) 
	  {
	      case 0: array = pressure4; break;
	      case 1: array = temperature4; break;
	      case 2: array = entropyDensity4; break;
	      case 3: array = mu4; break;
	      case 4: array = cs2_4; break;
	      default: fprintf(stderr,"ERROR in interpolate2D - selector must be 0,1,2,3, or 4, selector = %d \n", selector); exit(1);
	  }
    }
    else if(local_ed < EPP6)
    {
	  eps0 = EPP5;
	  NEps = NEPP5;
	  deltaEps = deltaEPP5;
        rhob0 = BNP5;
	  Nrhob = NBNP5;
	  deltaRhob = deltaBNP5;
	  switch (selector) 
	  {
	      case 0: array = pressure5; break;
	      case 1: array = temperature5; break;
	      case 2: array = entropyDensity5; break;
	      case 3: array = mu5; break;
	      case 4: array = cs2_5; break;
	      default: fprintf(stderr,"ERROR in interpolate2D - selector must be 0,1,2,3, or 4, selector = %d \n", selector); exit(1);
	  }
    }
    else if(local_ed < EPP7)
    {
	  eps0 = EPP6;
	  NEps = NEPP6;
	  deltaEps = deltaEPP6;
        rhob0 = BNP6;
	  Nrhob = NBNP6;
	  deltaRhob = deltaBNP6;
	  switch (selector) 
	  {
	      case 0: array = pressure6; break;
	      case 1: array = temperature6; break;
	      case 2: array = entropyDensity6; break;
	      case 3: array = mu6; break;
	      case 4: array = cs2_6; break;
	      default: fprintf(stderr,"ERROR in interpolate2D - selector must be 0,1,2,3, or 4, selector = %d \n", selector); exit(1);
	  }
    }
    else
    {
	  eps0 = EPP7;
	  NEps = NEPP7;
	  deltaEps = deltaEPP7;
        rhob0 = BNP7;
	  Nrhob = NBNP7;
	  deltaRhob = deltaBNP7;
	  switch (selector) 
	  {
	      case 0: array = pressure7; break;
	      case 1: array = temperature7; break;
	      case 2: array = entropyDensity7; break;
	      case 3: array = mu7; break;
	      case 4: array = cs2_7; break;
	      default: fprintf(stderr,"ERROR in interpolate2D - selector must be 0,1,2,3, or 4, selector = %d \n", selector); exit(1);
	  }
    }

    // compute the indices
    int idx_e = (int)((local_ed - eps0)/deltaEps);
    double frac_e = (local_ed - (idx_e*deltaEps + eps0))/deltaEps; 
    int idx_nb = (int)((local_rhob - rhob0)/deltaRhob);
    double frac_rhob = (local_rhob - (idx_nb*deltaRhob + rhob0))/deltaRhob; 

    // check overflow
    if(idx_e > (NEps-1) || idx_nb > (Nrhob-1))
    {
        fprintf(stderr, "ERROR in inperpolate2D: out of range of the table! \n");
        fprintf(stderr, "e = %e, rhob = %e, eps0 = %e, rhob0 = %e\n", 
                local_ed, local_rhob, eps0, rhob0);
        fprintf(stderr, "idx_e=%d, NEPP1=%d; idx_nb=%d, NBNP1=%d\n", 
                idx_e, NEps, idx_nb, Nrhob);
        fprintf(stderr, "deps = %f, drho = %f \n", deltaEps, deltaRhob);
        exit(0);
    }
    
    // check underflow
    if(idx_nb < 0)
    {
        fprintf(stderr, "ERROR in inperpolate2D: out of range of the table! \n");
        fprintf(stderr, "e = %e, rhob = %e, eps0 = %e, rhob0 = %e\n", 
                local_ed, local_rhob, eps0, rhob0);
        fprintf(stderr, "idx_e=%d, NEPP1=%d; idx_nb=%d, NBNP1=%d\n", 
                idx_e, NEps, idx_nb, Nrhob);
        fprintf(stderr, "deps = %f, drho = %f \n", deltaEps, deltaRhob);
        exit(0);
    }

    if(idx_e < 0)
    {
        if(local_ed > EPP1 + 1e-15 || local_ed < 0.)
        {
            fprintf(stderr, "ERROR in inperpolate2D: out of range of the table! \n");
            fprintf(stderr, "e = %e, rhob = %e, eps0 = %e, rhob0 = %e\n", 
                    local_ed, local_rhob, eps0, rhob0);
            fprintf(stderr, "idx_e=%d, NEPP1=%d; idx_nb=%d, NBNP1=%d\n", 
                    idx_e, NEps, idx_nb, Nrhob);
            fprintf(stderr, "deps = %f, drho = %f \n", deltaEps, deltaRhob);
            exit(0);
        }
        else
            idx_e = 0;  // do linear extrapolation for small energy density
    }

    double result;
    if(local_ed < EPP1)
    {
        double temp11 = max(array[idx_nb][idx_e], 0.0);
        double temp22 = max(array[idx_nb+1][idx_e], 0.0);
        result = (temp11*(1. - frac_rhob) + temp22*frac_rhob)*local_ed/EPP1;
    }
    else
    {
        double temp1 = max(array[idx_nb][idx_e], 0.0);
        double temp2 = max(array[idx_nb][idx_e+1], 0.0);
        double temp3 = max(array[idx_nb+1][idx_e+1], 0.0);
        double temp4 = max(array[idx_nb+1][idx_e], 0.0);

        result = (  (temp1*(1. - frac_e) + temp2*frac_e)*(1. - frac_rhob)
                  + (temp3*frac_e + temp4*(1. - frac_e))*frac_rhob);
    }
    
    // convert back to fm unit
    switch (selector) 
    {
        case 0: result /= hbarc; break;   // pressure in [1/fm^4]
        case 1: result /= hbarc; break;   // temperature in [1/fm]
        case 2: break;                    // entropy density in [1/fm^3]
        case 3: result /= hbarc; break;   // mu_B in [1/fm]
        case 4: break;                    // cs^2
	  default: fprintf(stderr,"ERROR in interpolate2D - selector must be 0,1,2,3, or 4, selector = %d \n", selector); exit(1);
    }

    return result;
}

double EOS::T_from_eps_ideal_gas(double eps)
{

	//Define number of colours and of flavours
 	const double Nc=3, Nf=2.5;
	
	return pow(90.0/M_PI/M_PI*(eps/3.0)/(2*(Nc*Nc-1)+7./2*Nc*Nf),.25);

}
double EOS::s2e_ideal_gas(double s)
{

	//Define number of colours and of flavours
 	const double Nc=3, Nf=2.5;

	//e=T*T*T*T*(M_PI*M_PI*3.0*(2*(Nc*Nc-1)+7./2*Nc*Nf)/90.0);
	//s = 4 e / (3 T)
	//s =4/3 T*T*T*(M_PI*M_PI*3.0*(2*(Nc*Nc-1)+7./2*Nc*Nf)/90.0);
	//T = pow(3. * s / 4. / (M_PI*M_PI*3.0*(2*(Nc*Nc-1)+7./2*Nc*Nf)/90.0), 1./3.);
	return 3. / 4. * s * pow(3. * s / 4. / (M_PI*M_PI*3.0*(2*(Nc*Nc-1)+7./2*Nc*Nf)/90.0), 1./3.); //in 1/fm^4

}

double EOS::get_entropy(double epsilon, double rhob)
// return entropy density in [1/fm^3]
{
    double f;
    double P, T, mu;
    
    if (whichEOS >=0 && whichEOS <=1)
    {
        P = get_pressure(epsilon, rhob);
        T = get_temperature(epsilon,rhob);
        mu = get_mu(epsilon, rhob);
   
        if (T!=0)
            f = (epsilon + P - mu*rhob)/(T + 1e-15);
        else
            f = 0.;
    }
    else if (whichEOS >= 2 && whichEOS < 10)
        f = interpolate2(epsilon,rhob,2);
    else if (whichEOS >= 10)
        f = interpolate2D(epsilon, fabs(rhob), 2);  // EOS is symmetric in rho_b
    else
    {
        fprintf(stderr,"EOS::get_entropy: whichEOS = %d is out of range!\n", whichEOS);
        exit(0);
    }
    return f;
}/* get_entropy */

double EOS::ssolve(double e, double rhob, double s)
{
  // takes e in GeV/fm^3 and passes it on in 1/fm^4 ...
  double P, T, mu;
  P = get_pressure(e/hbarc, rhob);
  T = get_temperature(e/hbarc, rhob);
  mu = get_mu(e/hbarc, rhob);
      
//   fprintf(stderr,"T=%f\n",T*hbarc);
//   fprintf(stderr,"P=%f\n",P*hbarc);
//   fprintf(stderr,"e=%f\n",e);
//   fprintf(stderr,"s=%f\n",s);
//   fprintf(stderr,"ssolve=%f\n",T*s-e/hbarc-P+mu*rhob);

  return T*s-e/hbarc-P+mu*rhob;
}

double EOS::Tsolve(double e, double rhob, double T)
{
  // takes e in GeV/fm^3 and passes it on in 1/fm^4 ...
  return T-get_temperature(e/hbarc,rhob);
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
      if (fabs(de)<eacc || fmid==0.) return value/hbarc;
    }
  fprintf(stderr,"too many bisections in findRoot\n");
  return 0.;
}


double EOS::get_temperature(double eps, double rhob)
// return temperature in [1/fm]
{
    double T;
    if (whichEOS==0)
        T = T_from_eps_ideal_gas(eps);
    else if (whichEOS==1)
        T = interpolate(eps, rhob, 0);
    else if (whichEOS>=2 && whichEOS < 10)
        T = interpolate2(eps, rhob, 1);
    else if (whichEOS >= 10)
        T = interpolate2D(eps, fabs(rhob), 1);  // EOS is symmetric in rho_b
    else
    {
        fprintf(stderr,"EOS::get_temperature: whichEOS = %d is out of range!\n", whichEOS);
        exit(0);
    }
    return T;
}


double EOS::get_mu(double eps, double rhob)
// return mu_B in [1/fm]
{
    double mu;
    if (whichEOS==0)
        mu=0.0;
    else if (whichEOS==1)
        mu = 0.0;
    else if (whichEOS>=2 && whichEOS < 10)
        mu = 0.0;
    else if (whichEOS >= 10)
    {
        if(rhob < 0.0)    // EOS is anti-symmetric in rho_b for mu_B
            mu = -interpolate2D(eps, -rhob, 3);
        else
            mu = interpolate2D(eps, rhob, 3);
    }
    else
    {
        fprintf(stderr,"EOS::get_mu: whichEOS = %d is out of range!\n", whichEOS);
        exit(0);
    }
    return mu;
}


double EOS::get_qgp_frac(double eps, double rhob) {

 double frac;

 if (whichEOS==0)
   {
     frac=-1.0;
   }
 else if (whichEOS==1)
   {
     frac=(eps*hbarc-0.45)/(1.6-0.45); // e=(1-QGPfrac)*e_H + QGPfrac*e_QGP in the mixed phase
     if (frac>1.) frac = 1;
     else if (frac<0.) frac=0.;
   }
 else if (whichEOS>=2)
   {
     frac = interpolate2(eps, rhob, 3);
   }
   	else {fprintf(stderr,"whichEOS out of range.\n");exit(0);}

  return frac;

}


double EOS::get_s2e(double s, double rhob) { //s - entropy density in 1/fm^3

	double e; //epsilon - energy density
	int status;

	//This function does not work for non-zero baryon density
	if (rhob > 0.0) {
		cerr << "Function get_s2e() is only implemented for rhob=0. Aborting...\n";
		exit(1);
	}

	if (whichEOS==0) {
		e=s2e_ideal_gas(s);
	}
	else if (whichEOS>=2 && whichEOS<=6) {
		status=gsl_interp_eval_e(interp_s2e, s_list_rho0, eps_list_rho0, s, accel_s2e, &e);
		if (status == GSL_EDOM) {
			cerr << "Error: can't get energy from entropy, entropy s="<<s<<" fm^-3 is outside the current tabulation of the EOS...\n";
			exit(1);
		}
	}
   	else {fprintf(stderr,"whichEOS out of range.\n");exit(0);}

	return e; //in 1/fm^4

}
