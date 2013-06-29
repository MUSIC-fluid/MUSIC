#include "freeze.h"

using namespace std;

Freeze::Freeze()
{
  integral = new Int;
  util = new Util;
}

// destructors
Freeze::~Freeze()
{
  delete integral;
  delete util;
}

void Freeze::checkForReadError(FILE *file, const char* name)
{
  if(!(file))
    {
      fprintf(stderr, "file %s not found.\n", name);
      fprintf(stderr, "Exiting...\n");
      exit(0);
    }
}


void Freeze::ReadParticleData(InitData *DATA, EOS *eos)
{
  // read in particle and decay information from file:
  // partid = (int *)malloc(MAXINTV * sizeof(int)); 
  partid = new int[MAXINTV * sizeof(int)]; 
  int bytes_read;
  static char *s;
  s = new char[120];
  int i, j, k, d1, d2, d3, decays, h;
  double b, npi, nK, neta, dummy;
  fprintf(stderr,"reading particle data\n");
  char *anti;
  // open particle data file:
  const char* EOSPATH = "HYDROPROGRAMPATH";
  char* envPath = getenv(EOSPATH);
  char* p_name;
  if (envPath != 0 && *envPath != '\0') 
    {
      p_name = util->char_malloc(100);
      strcat(p_name, envPath);
      strcat(p_name,"/EOS/pdg05.dat");
    }  
  else
    {
      p_name = util->char_malloc(100);
      strcat(p_name, ".");
      strcat(p_name,"/EOS/pdg05.dat");
    }
  
  cout << "from " << p_name << endl; 

  FILE *p_file;
  p_file = fopen(p_name, "r");
  checkForReadError(p_file,p_name);
  
  for(k=0;k<MAXINTV;k++) 
    partid[k] = -1; 

  i=0;
  j=0;
  cout << "before" << endl; cout << "sizeofParticle=" << sizeof(Particle)/1000000 << endl;
  particleList = (Particle *)malloc((DATA->NumberOfParticlesToInclude+1) * sizeof(Particle));
  
//   particleList = new Particle[(DATA->NumberOfParticlesToInclude)];
  cout <<"after first (check if there is enough memory... seg fault may be due to lack of memory)" << endl; 
  
  // read particle data:
  while(i<DATA->NumberOfParticlesToInclude+1)
    {
      particleList[i].name = util->char_malloc(50);
      //particleList[i].name = new char[50];
      bytes_read=fscanf(p_file, "%d", &particleList[i].number);
      bytes_read=fscanf(p_file, "%s", particleList[i].name);
      bytes_read=fscanf(p_file, "%lf", &particleList[i].mass);
      bytes_read=fscanf(p_file, "%lf", &particleList[i].width);
      bytes_read=fscanf(p_file, "%d", &particleList[i].degeneracy);
      bytes_read=fscanf(p_file, "%d", &particleList[i].baryon);
      bytes_read=fscanf(p_file, "%d", &particleList[i].strange);
      bytes_read=fscanf(p_file, "%d", &particleList[i].charm);
      bytes_read=fscanf(p_file, "%d", &particleList[i].bottom);
      bytes_read=fscanf(p_file, "%d", &particleList[i].isospin);
      bytes_read=fscanf(p_file, "%lf", &particleList[i].charge);
      bytes_read=fscanf(p_file, "%d", &particleList[i].decays); // number of decays
      
      //Bytes_read=fscanf(p_file, "%lf", &npi);
      //bytes_read=fscanf(p_file, "%lf", &nK);
      //bytes_read=fscanf(p_file, "%lf", &neta);
   
   /*    fprintf(stderr,"%s %i \n",particleList[i].name,particleList[i].number); */
/*       fprintf(stderr,"%lf %lf %d %d %d %d %d %d %lf %d \n",particleList[i].mass,particleList[i].width, */
/* 	      particleList[i].degeneracy,particleList[i].baryon,particleList[i].strange, */
/* 	      particleList[i].charm,particleList[i].bottom,particleList[i].isospin,particleList[i].charge,particleList[i].decays); */
      
      partid[MHALF + particleList[i].number] = i;
      
      particleList[i].stable = 0;
      
      for(k=0;k<particleList[i].decays;k++) 
	{
	  h=fscanf(p_file,"%i%i%lf%i%i%i%i%i",
		   &decay[j].reso, &decay[j].numpart, &decay[j].branch, 
		   &decay[j].part[0], &decay[j].part[1], &decay[j].part[2],
		   &decay[j].part[3], &decay[j].part[4]);
/* 	  fprintf(stderr,"%i %i %lf %i %i %i %i %i\n", */
/* 		   decay[j].reso, decay[j].numpart, decay[j].branch, */
/* 		   decay[j].part[0], decay[j].part[1], decay[j].part[2], */
/* 		   decay[j].part[3], decay[j].part[4]); */
	 
	  if (h != 8) 
	    {
	      printf("Error in scanf decay \n");
	      exit(0);
	    }
	  if (decay[j].numpart == 1) particleList[i].stable = 1; //"decays" into one particle, i.e. is stable 
	  j++; // increase the decay counting variable "j" by 1
	}
	
      // include anti-baryons (there are none in the file)
      //fprintf(stderr,"[%i] has pid %i \n",MHALF + particleList[i].number,partid[MHALF + particleList[i].number]);
      //if (particleList[i].number==960225) sleep(1);
      if ( particleList[i].baryon!=0 )
	{
	  //fprintf(stderr,"b=%d\n",particleList[i].baryon);
	  i++;
	  particleList[i].name = util->char_malloc(20);
	  anti = util->char_malloc(30);
	  strcat(anti,"Anti-");
	  strcat(anti,particleList[i-1].name);
	  particleList[i].width    =  particleList[i-1].width;
	  particleList[i].charm    = -particleList[i-1].charm;
	  particleList[i].bottom   = -particleList[i-1].bottom;
	  particleList[i].isospin =  particleList[i-1].isospin;
	  particleList[i].charge   = -particleList[i-1].charge;
	  particleList[i].decays   = particleList[i-1].decays;
	  particleList[i].stable =  particleList[i-1].stable;
	  particleList[i].number = -particleList[i-1].number;
	  particleList[i].name = anti; 
	  particleList[i].mass = particleList[i-1].mass;
	  particleList[i].degeneracy = particleList[i-1].degeneracy;
	  particleList[i].baryon = -particleList[i-1].baryon;
	  particleList[i].strange = -particleList[i-1].strange;
	  particleList[i].charge = -particleList[i-1].charge;
	  partid[MHALF + particleList[i].number] = i;
	  	  //  fprintf(stderr,"%s %i \n",particleList[i].name,particleList[i].number);
	  //fprintf(stderr,"[%i] has pid %i \n",MHALF + particleList[i].number,partid[MHALF + particleList[i].number]);
	}
      i++;
    } 
  decayMax = j;
  particleMax = i-1;
  fclose(p_file);
  // here read the stable particles' chemical potential at freeze-out
  if (DATA->whichEOS>=3)
    {
      double ef;
      if (1==DATA->useEpsFO) {
        ef=DATA->epsilonFreeze;
      }
      else {
        cout << "determining epsFO from TFO=" << DATA->TFO << endl; 
        ef= eos->findRoot(&EOS::Tsolve, 0., DATA->TFO/hbarc, 0.001, 300.,0.001)*hbarc;
        cout << "freeze out energy density (assuming rhob=0) is " << ef << endl;
      }

      cout << "Determining chemical potentials at freeze out energy density " << ef << " GeV/fm^3." << endl;
      
      const char* EOSPATH = "HYDROPROGRAMPATH";
      char* envPath = getenv(EOSPATH);
      char* mu_name;
      if (envPath != 0 && *envPath != '\0') 
	{
	  mu_name = util->char_malloc(100);
	  strcat(mu_name, envPath);
	  if (DATA->whichEOS==3)
	    strcat(mu_name,"/EOS/s95p-PCE-v1/s95p-PCE-v1_pichem1.dat");
	  else if (DATA->whichEOS==4)
	    strcat(mu_name,"/EOS/s95p-PCE155/pichem1.dat");
	  else if (DATA->whichEOS==5)
	    strcat(mu_name,"/EOS/s95p-PCE160/pichem1.dat");
	  else if (DATA->whichEOS==6)
	    strcat(mu_name,"/EOS/s95p-PCE165-v0/s95p-PCE165-v0_pichem1.dat");
	}  
      else
	{
	  mu_name = util->char_malloc(100);
	  strcat(mu_name, ".");
	  if (DATA->whichEOS==3)
	    strcat(mu_name,"/EOS/s95p-PCE-v1/s95p-PCE-v1_pichem1.dat");
	  else if (DATA->whichEOS==4)
	    strcat(mu_name,"/EOS/s95p-PCE155/pichem1.dat");
	  else if (DATA->whichEOS==5)
	    strcat(mu_name,"/EOS/s95p-PCE160/pichem1.dat");
	  else if (DATA->whichEOS==6)
	    strcat(mu_name,"/EOS/s95p-PCE165-v0/s95p-PCE165-v0_pichem1.dat");
	}
      
      cout << "Reading chemical potentials from file\n " << mu_name << endl; 
      
      FILE *mu_file;
      mu_file = fopen(mu_name, "r");
      checkForReadError(mu_file,mu_name);
      double BNP1, EPP1;            // start value for \mu_B and epsilon
      double deltaBNP1, deltaEPP1;  // step size for \mu_B and epsilon
      int NBNP1, NEPP1;             // number of entries for \mu_B and epsilon
      int numStable;                // number of stable particles (number of columns in the file)
      double **chemPot;

      bytes_read=fscanf(mu_file,"%lf",&EPP1);
      bytes_read=fscanf(mu_file,"%lf %d",&deltaEPP1,&NEPP1);
      bytes_read=fscanf(mu_file,"%d",&numStable);
      
      //      cout << "EPP1=" << EPP1 << ", deltaEPP1=" << deltaEPP1 << ", NEPP1=" << NEPP1 << ", numStable=" << numStable << endl;
      
      chemPot=util->mtx_malloc(numStable+1,NEPP1+1); // chemical potential for every stable particle 

      for(j=NEPP1-1; j>=0; j--)
	{
	  for(i=0; i<numStable; i++)
	    {
	      bytes_read=fscanf(mu_file,"%lf",&chemPot[i][j]);
	      //      cout << chemPot[i][j] << " ";
	    }
	  //	  cout << endl;
	}

      double frace;
      int ie1, ie2;
      if(ef<EPP1) 
	{
	  ie1 = 0;
	  ie2 = 1;
	  frace = ef/(EPP1);
	}
      else
	{
	  ie1 = floor((ef-EPP1)/deltaEPP1);
	  ie2 = floor((ef-EPP1)/deltaEPP1+1);
	  frace = (ef-(ie1*deltaEPP1+EPP1))/deltaEPP1; 
	}
      
      if(ie1>NEPP1)
	{
	  fprintf(stderr,"ERROR in ReadParticleData. out of range.\n");
	  fprintf(stderr,"ie1=%d,NEPP1=%d\n", ie1, NEPP1);
	  exit(0);
	}
      if(ie2>NEPP1)
	{
	  fprintf(stderr,"ERROR in ReadParticleData. out of range.\n");
	  fprintf(stderr,"ie2=%d,NEPP1=%d\n", ie2, NEPP1);
	  exit(0);
	}

      double pa, pb;
      double mu[numStable+1];
      cout << "numStable=" << numStable << endl;
      
      for(i=1; i<=numStable; i++)
	{
	  pa = chemPot[i-1][ie1];
	  pb = chemPot[i-1][ie2];
	  
	  if(ef<EPP1) 
	    {
	      mu[i] = pa*(frace);
	      //if (p<0) fprintf(stderr,"pa=%lf\n", pa);
	      //if (p<0) fprintf(stderr,"p=%lf\n", p);
	    }
	  else
	    {
	      mu[i] = pa*(1-frace) + pb*frace;
	    }
	  //	  cout << "mu of stable particle " << i << "=" << mu[i] << endl;
	}
     
      for(i=0; i<DATA->NumberOfParticlesToInclude; i++)
        {
          particleList[i].muAtFreezeOut = 0.;      
        }
      
      if (DATA->NumberOfParticlesToInclude>=17)
        {
          for(i=1; i<=9; i++)
            {
              particleList[i].muAtFreezeOut = mu[i];
            }
        }
      else
        {
          cout << "Need at least 9 particles. Increase number of particles to include. Exiting." << endl;
          exit(1);
        }
      if (DATA->NumberOfParticlesToInclude>=17)
      particleList[17].muAtFreezeOut = mu[10];
      if (DATA->NumberOfParticlesToInclude>=18)
      particleList[18].muAtFreezeOut = mu[11];
      if (DATA->NumberOfParticlesToInclude>=19)
      particleList[19].muAtFreezeOut = mu[12];
      if (DATA->NumberOfParticlesToInclude>=20)
      particleList[20].muAtFreezeOut = mu[13];
      if (DATA->NumberOfParticlesToInclude>=21)
      particleList[21].muAtFreezeOut = mu[14];

      if (DATA->NumberOfParticlesToInclude>=26)
      particleList[26].muAtFreezeOut = mu[15];
      if (DATA->NumberOfParticlesToInclude>=27)
      particleList[27].muAtFreezeOut = mu[16];
      if (DATA->NumberOfParticlesToInclude>=28)
      particleList[28].muAtFreezeOut = mu[17];

      if (DATA->NumberOfParticlesToInclude>=30)
      particleList[30].muAtFreezeOut = mu[18];
      if (DATA->NumberOfParticlesToInclude>=31)
      particleList[31].muAtFreezeOut = mu[19];
      if (DATA->NumberOfParticlesToInclude>=32)
      particleList[32].muAtFreezeOut = mu[20];
      if (DATA->NumberOfParticlesToInclude>=33)
      particleList[33].muAtFreezeOut = mu[21];
      if (DATA->NumberOfParticlesToInclude>=34)
      particleList[34].muAtFreezeOut = mu[22];
      if (DATA->NumberOfParticlesToInclude>=35)
      particleList[35].muAtFreezeOut = mu[23];

      if (DATA->NumberOfParticlesToInclude>=60)
      particleList[60].muAtFreezeOut = mu[24];
      if (DATA->NumberOfParticlesToInclude>=61)
      particleList[61].muAtFreezeOut = mu[25];
      if (DATA->NumberOfParticlesToInclude>=62)
      particleList[62].muAtFreezeOut = mu[26];
      if (DATA->NumberOfParticlesToInclude>=63)
      particleList[63].muAtFreezeOut = mu[27];

      if (DATA->NumberOfParticlesToInclude>=110)
      particleList[110].muAtFreezeOut = mu[28];
      if (DATA->NumberOfParticlesToInclude>=111)
      particleList[111].muAtFreezeOut = mu[29];

      if (DATA->NumberOfParticlesToInclude>=117)
      particleList[117].muAtFreezeOut = mu[30];
      if (DATA->NumberOfParticlesToInclude>=118)
      particleList[118].muAtFreezeOut = mu[31];
      if (DATA->NumberOfParticlesToInclude>=119)
      particleList[119].muAtFreezeOut = mu[32];
      if (DATA->NumberOfParticlesToInclude>=120)
      particleList[120].muAtFreezeOut = mu[33];

      if (DATA->NumberOfParticlesToInclude>=170)
      particleList[170].muAtFreezeOut = mu[34];
      if (DATA->NumberOfParticlesToInclude>=171)
      particleList[171].muAtFreezeOut = mu[35];

      
      cout << "Got the chemical potentials at freeze out for stable particles." << endl;
 
      k=0;
      for(i=1; i<DATA->NumberOfParticlesToInclude; i++) //skip the photon (i=0)
	{
	  if (particleList[i].muAtFreezeOut==0)
	    {
	      for(j=1; j<=particleList[i].decays; j++)
		{
		  if (particleList[i].baryon > -1)
		    k++;

		  for (int m=0; m<abs(decay[k].numpart); m++)
		    {
		      
		      particleList[i].muAtFreezeOut += decay[k].branch*particleList[partid[MHALF+decay[k].part[m]]].muAtFreezeOut;
		      
		      //        cout <<  particleList[partid[MHALF+decay[k].part[m]]].name << endl;
		    }
		}
	    }
	  else
	    {
	      for(j=1; j<=particleList[i].decays; j++)
		{
		  if (particleList[i].baryon > -1)
		    k++;
		  //cout << "k=" << k << " " << particleList[partid[MHALF+decay[k].reso]].name << "=" << decay[k].reso << endl;
		  //cout << particleList[i].name << endl;		  
		}
	    }
	}

//       for(i=1; i<DATA->NumberOfParticlesToInclude; i++)
// 	{
// 	  cout << particleList[i].name << " " << particleList[i].muAtFreezeOut << endl;
// 	}
      
      cout << "Got the chemical potentials at freeze-out for the resonances." << endl;

    }

  

  cout << "Done reading particle data." << endl;
  // fprintf(stderr,"decayMax=%d\n",decayMax);
  // fprintf(stderr,"particleMax=%d\n",particleMax);
}


int Freeze::countLines (std::istream& in)
{
    return std::count(std::istreambuf_iterator<char>(in),
                      std::istreambuf_iterator<char>(),
                      '\n');
}

void Freeze::ReadFreezeOutSurface(InitData *DATA)
{
  size_t nbytes=100;
  int i=0;
  int bytes_read;
  char * dummy;
  double test;
  fprintf(stderr,"reading freeze-out surface\n");
  // open particle data file:
  FILE *s_file;
  char* line;
  const char* s_name = "./surface.dat";
  s_file = fopen(s_name, "r");
  checkForReadError(s_file,s_name);
//   // open geometry data file:
//   FILE *g_file;
//   char* g_name = "./geometry.dat"; 
//   g_file = fopen(g_name, "r");
//   double tecc2, tecc3, tecc3r3, tPsi2, tPsi3, tPsi3r3;
//   
//   if(!(g_file))
//      {
//        cout << "file " << g_name << " not found." << endl;
//        tecc2 = 1.; 
//        tecc3 = 1.; 
//        tecc3r3 = 1.; 
//        tPsi2 = 0.;
//        tPsi3 = 0.; 
//        tPsi3r3 = 0.; 
//      }
//    else
//      {
//        bytes_read=fscanf(g_file, "%lf", &tecc2);
//        bytes_read=fscanf(g_file, "%lf", &tPsi2);
//        bytes_read=fscanf(g_file, "%lf", &tecc3);
//        bytes_read=fscanf(g_file, "%lf", &tPsi3);
//        //bytes_read=fscanf(g_file, "%lf", &tecc3r3);
//        //bytes_read=fscanf(g_file, "%lf", &tPsi3r3);
//        fclose(g_file);
//    }
//   
//    DATA->ecc2 = tecc2;
//    DATA->ecc3 = tecc3;
//    DATA->ecc3r3 = tecc3r3;
//    DATA->Psi2 = tPsi2;
//    //DATA->Psi3 = tPsi3;
//    //DATA->Psi3r3 = tPsi3r3;
//    
//    cout << "e2=" << DATA->ecc2 << ", Psi2=" << DATA->Psi2 << ", e3=" << DATA->ecc3 << ", Psi3=" << DATA->Psi3 << endl;
//    //	<< ", e3r3=" << DATA->ecc3 << ", Psi3r3=" << DATA->Psi3r3 << endl;
//    
  NCells=0;
  /*   int bytes_read; */
  //  fprintf(s_file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",
  //  tauf, xf, yf, etaf, FULLSU[0],FULLSU[1],FULLSU[2],FULLSU[3],
  //  utau, ux, uy, ueta, epsFO, TFO, muB);
  //  line = (char *) malloc(nbytes + 1);
  //   while( getline(&line, &nbytes, s_file)!=-1 )
  //     {
  //       NCells++;
  //     }

  // new counting, mac compatible ...
  ifstream in;
  in.open(s_name);
  NCells = 0;
  NCells += countLines(in);
  fprintf(stderr,"NCells=%d\n",NCells);
  fclose(s_file);

  s_file = fopen(s_name, "r");
  checkForReadError(s_file,s_name);

  // Now allocate memory: array of surfaceElements with length NCells
  surface = (SurfaceElement *) malloc((NCells)*sizeof(SurfaceElement));
//   bytes_read=fscanf(s_file, "%s", dummy);
//   bytes_read=fscanf(s_file, "%s", dummy);
//   bytes_read=fscanf(s_file, "%s", dummy);
//   bytes_read=fscanf(s_file, "%s", dummy);
//   bytes_read=fscanf(s_file, "%s", dummy);
//   bytes_read=fscanf(s_file, "%s", dummy);
//   bytes_read=fscanf(s_file, "%s", dummy);
//   bytes_read=fscanf(s_file, "%s", dummy);
//   bytes_read=fscanf(s_file, "%s", dummy);
//   bytes_read=fscanf(s_file, "%s", dummy);
//   bytes_read=fscanf(s_file, "%s", dummy);
//   bytes_read=fscanf(s_file, "%s", dummy);
//   bytes_read=fscanf(s_file, "%s", dummy);
//   bytes_read=fscanf(s_file, "%s", dummy);
//   bytes_read=fscanf(s_file, "%s", dummy);
 
  while(i<NCells)
    {
      // position in (tau, x, y, eta)
      bytes_read=fscanf(s_file, "%lf", &surface[i].x[0]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].x[1]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].x[2]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].x[3]);
      // hypersurface vector in (tau, x, y, eta)
      bytes_read=fscanf(s_file, "%lf", &surface[i].s[0]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].s[1]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].s[2]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].s[3]);
      // flow velocity in (tau, x, y, eta)
      bytes_read=fscanf(s_file, "%lf", &surface[i].u[0]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].u[1]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].u[2]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].u[3]);
      // freeze-out energy density
      bytes_read=fscanf(s_file, "%lf", &surface[i].epsilon_f);
      if(surface[i].epsilon_f<0) 
	cout << "WARNING: epsilon-f<0." << endl;
      // freeze-out temperature
      bytes_read=fscanf(s_file, "%lf", &surface[i].T_f);
      if(surface[i].T_f<0) 
	cout << "WARNING: T_f<0." << endl;
      // freeze-out baryon chemical potential
      bytes_read=fscanf(s_file, "%lf", &surface[i].mu_B);
      // freeze-out entropy density s
      bytes_read=fscanf(s_file, "%lf", &surface[i].eps_plus_p_over_T_FO);
      // freeze-out Wmunu
      bytes_read=fscanf(s_file, "%lf", &surface[i].W[0][0]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].W[0][1]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].W[0][2]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].W[0][3]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].W[1][1]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].W[1][2]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].W[1][3]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].W[2][2]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].W[2][3]);
      bytes_read=fscanf(s_file, "%lf", &surface[i].W[3][3]);
      i++;
    }				
  fclose(s_file);
 
//  for(i=0; i<NCells; i++) 
//     { 
//       fprintf(stderr,"%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", 
// 	      surface[i].x[0],surface[i].x[1],surface[i].x[2],surface[i].x[3],surface[i].s[0],surface[i].s[1],surface[i].s[2],surface[i].s[3], 
//  	      surface[i].u[0],surface[i].u[1],surface[i].u[2],surface[i].u[3],surface[i].epsilon_f,surface[i].T_f,surface[i].mu_B,surface[i].sFO,
// 	      surface[i].W[0][0],surface[i].W[0][1],surface[i].W[0][2],surface[i].W[0][3],surface[i].W[1][1],
// 	      surface[i].W[1][2],surface[i].W[1][3],surface[i].W[2][2],surface[i].W[2][3],surface[i].W[3][3]); 
//       sleep(2);
//     }
}

// read in thermal spectra from file to then perform resonance decays with them
void Freeze::ReadSpectra(InitData* DATA)
{
  // read in thermal spectra from file:
  int number, iymax, iptmax, iphimax;
  double deltaY, ymax, slope, phimax, phimin;
  double deltaeta, etamax, ptmin, ptmax;
  int pseudo_steps;
  int ip, iphi, ipt, i;
  double *p, *w;		        // pointing to data for Gaussian integration in phi 
  int bytes_read;
  int iy, j, k, d1, d2, d3, decays, h;
  double b, npi, nK, neta, dummy;
  fprintf(stderr,"reading spectra\n");
  char *anti;
  // open particle information file:
  FILE *p_file;
  const char* p_name = "particleInformation.dat";
  p_file = fopen(p_name, "r");
  checkForReadError(p_file,p_name);
  int count;
  count = 0;
  cout << "NumberOfParticlesToInclude " << DATA->NumberOfParticlesToInclude << endl;
  // read particle information:
  while( fscanf(p_file,"%d %lf %lf %lf %lf %lf %d %d %d ",&number, &deltaY, &ymax, &slope, &phimin, &phimax, &iymax, &iptmax, &iphimax) == 9)
    {
      count ++;
      if (count>DATA->NumberOfParticlesToInclude) break;
      fprintf(stderr,"%d %e %e %e %e %d %d %d \n", number, deltaY, ymax, slope, phimax, iymax, iptmax, iphimax);
      ip = partid[MHALF+number];
      particleList[ip].ny = iymax;
      particleList[ip].npt = iptmax;
      particleList[ip].nphi = iphimax;
      particleList[ip].phimax = phimax;
      particleList[ip].phimin = phimin;
      particleList[ip].slope = slope;
      particleList[ip].ymax = ymax;
      particleList[ip].deltaY = deltaY;
      for ( i=0; i<iptmax; i++ )
	{
	  particleList[ip].pt[i] =  gala15x[i]/slope;
	}
      for ( i=0; i<iymax; i++ )
	{
	  particleList[ip].y[i] =  i*deltaY-ymax+deltaY/2.;
	  //	  if (ip==1) cout << "read particleList[ip].y[" << i << "] = " <<  particleList[ip].y[i] << endl;
	}
 
	  switch (iphimax) 
	    {
	    case 4: p= gaulep4; w= gaulew4; break;
	    case 8: p= gaulep8; w= gaulew8; break;
	    case 10: p= gaulep10; w= gaulew10; break;
	    case 12: p= gaulep12; w= gaulew12; break;
	    case 16: p= gaulep16; w= gaulew16; break;
	    case 20: p= gaulep20; w= gaulew20; break;
	    case 48: p= gaulep48; w= gaulew48; break;
	    default: fprintf(stderr,"specified number of phi-points not available\n"); exit(1);
	    }
	  
	  phiArray = util->vector_malloc(iphimax);
	  for(iphi=0; iphi<iphimax; iphi++)
	    {
	      if ( iphi < iphimax/2 )
		{
		  phiArray[iphi] = (phimax-phimin)/2.*(1.-p[iphi])+phimin;
		}
	      else
		{
		  phiArray[iphi] = (phimax-phimin)/2.*(1.+p[iphimax-iphi-1])+phimin;
		}
	      //      cout << "phi[" << iphi << "]=" << phiArray[iphi] << ", " <<  phimax/2.*(1.-p[iphi]) << ", " 
	      //   << (phimax)/2.*(1.+p[iphimax-iphi-1]) << endl;
	    }
      }
  
  particleMax = ip;

  fclose(p_file);

  FILE *s_file;
  const char* s_name = "yptphiSpectra.dat";
  s_file = fopen(s_name, "r");
  checkForReadError(s_file,s_name);
  
  cout << "iymax=" << iymax << endl;
  cout << "cells=" << iymax*iptmax*iphimax << endl;
  for ( ip=1; ip<=particleMax; ip++ )
    {
      //cout << ip << endl;
      fprintf(stderr,"reading particle %d: %d %s\n", ip, particleList[ip].number, particleList[ip].name);
      for (iy=0; iy<iymax; iy++)
	{
	  for (ipt=0; ipt<iptmax; ipt++)
	    {
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  bytes_read=fscanf(s_file, "%lf", &particleList[ip].dNdydptdphi[iy][ipt][iphi]);
		  if(particleList[ip].dNdydptdphi[iy][ipt][iphi]<0.)
		    particleList[ip].dNdydptdphi[iy][ipt][iphi]=0;
			  //		  cout << particleList[ip].y[iy] << " " << particleList[ip].pt[ipt] << " " << phiArray[iphi] << " " << particleList[ip].dNdydptdphi[iy][ipt][iphi] << endl; 
		  //	  printf("%f %f %f \n",particleList[ip].y[iy],particleList[ip].pt[ipt],phiArray[iphi]);
		}
	    }
	}
    }
  //particleMax=2;
  fclose(s_file);
}


// read in thermal spectra from file to then perform resonance decays or something else with them
void Freeze::Read3Spectra(InitData* DATA) // read pion, kaon, proton
{
  // read in thermal spectra from file:
  int number, iymax, iptmax, iphimax;
  double deltaY, ymax, slope, phimax, phimin;
  int ip, iphi, ipt, i;
  double *p, *w;		        // pointing to data for Gaussian integration in phi 
  int bytes_read;
  static char *s;
  s = util->char_malloc(120);
  int iy, j, k, d1, d2, d3, decays, h;
  double b, npi, nK, neta, dummy;
  fprintf(stderr,"reading spectra\n");
  char *anti;
  // open particle information file:
  FILE *p_file;
  const char* p_name = "particleInformation.dat";
  p_file = fopen(p_name, "r");
  checkForReadError(p_file,p_name);
  int count;
  count = 0;
  // read particle information:
  while( fscanf(p_file,"%d %lf %lf %lf %lf %lf %d %d %d ",&number, &deltaY, &ymax, &slope, &phimin, &phimax, &iymax, &iptmax, &iphimax) == 9)
    {
      count ++;
      if (count>DATA->NumberOfParticlesToInclude) break;
      fprintf(stderr,"%d %e %e %e %e %d %d %d \n", number, deltaY, ymax, slope, phimax, iymax, iptmax, iphimax);
      ip = count;
      particleList[ip].number = number;
      particleList[ip].ny = iymax;
      particleList[ip].npt = iptmax;
      particleList[ip].nphi = iphimax;
      particleList[ip].phimax = phimax;
      particleList[ip].phimin = phimin;
      particleList[ip].slope = slope;
      particleList[ip].ymax = ymax;
      particleList[ip].deltaY = deltaY;
      for ( i=0; i<iptmax; i++ )
	{
	  particleList[ip].pt[i] =  gala15x[i]/slope;
	}
      for ( i=0; i<iymax; i++ )
	{
	  particleList[ip].y[i] =  i*deltaY-ymax+deltaY/2.;
	  //	  if (ip==1) cout << "read particleList[ip].y[" << i << "] = " <<  particleList[ip].y[i] << endl;
	}
      
      switch (iphimax) 
	{
	case 4: p= gaulep4; w= gaulew4; break;
	case 8: p= gaulep8; w= gaulew8; break;
	case 10: p= gaulep10; w= gaulew10; break;
	case 12: p= gaulep12; w= gaulew12; break;
	case 16: p= gaulep16; w= gaulew16; break;
	case 20: p= gaulep20; w= gaulew20; break;
	case 48: p= gaulep48; w= gaulew48; break;
	default: fprintf(stderr,"specified number of phi-points not available\n"); exit(1);
	}
      
      phiArray = util->vector_malloc(iphimax);
      for(iphi=0; iphi<iphimax; iphi++)
	{
	  if ( iphi < iphimax/2 )
	    {
	      phiArray[iphi] = (phimax-phimin)/2.*(1.-p[iphi])+phimin;
	    }
	  else
	    {
	      phiArray[iphi] = (phimax-phimin)/2.*(1.+p[iphimax-iphi-1])+phimin;
	    }
	  //      cout << "phi[" << iphi << "]=" << phiArray[iphi] << ", " <<  phimax/2.*(1.-p[iphi]) << ", " 
	  //   << (phimax)/2.*(1.+p[iphimax-iphi-1]) << endl;
	}
    }
  
  particleMax = ip;

  fclose(p_file);

  FILE *s_file;
  const char* s_name = "yptphiSpectra.dat";
  s_file = fopen(s_name, "r");
  checkForReadError(s_file,s_name);
  
  cout << "iymax=" << iymax << endl;
  cout << "cells=" << iymax*iptmax*iphimax << endl;
  for ( ip=1; ip<=particleMax; ip++ )
    {
      //cout << ip << endl;
      fprintf(stderr,"reading particle %d: %d \n", ip, particleList[ip].number);
      for (iy=0; iy<iymax; iy++)
	{
	  for (ipt=0; ipt<iptmax; ipt++)
	    {
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  bytes_read=fscanf(s_file, "%lf", &particleList[ip].dNdydptdphi[iy][ipt][iphi]);
		}
	    }
	}
    }
  //particleMax=2;
  fclose(s_file);
}


// read in thermal spectra from file to then perform resonance decays with them
void Freeze::ReadSingleSpectrum(InitData* DATA)
{
  // read in thermal spectra from file:
  int number, iymax, iptmax, iphimax;
  double deltaY, ymax, slope, phimax, phimin;
  int ip, iphi, ipt, i;
  double *p, *w;		        // pointing to data for Gaussian integration in phi 
  int bytes_read;
  static char *s;
  s = util->char_malloc(120);
  int iy, j, k, d1, d2, d3, decays, h;
  double b, npi, nK, neta, dummy;
  fprintf(stderr,"reading spectra\n");
  char *anti;
  // open particle information file:
  FILE *p_file;
  const char* p_name = "particleInformation.dat";
  p_file = fopen(p_name, "r");
  checkForReadError(p_file,p_name);
  int count;
  count = 0;
  // read particle information:
  while( fscanf(p_file,"%d %lf %lf %lf %lf %lf %d %d %d ",&number, &deltaY, &ymax, &slope, &phimin, &phimax, &iymax, &iptmax, &iphimax) == 9)
    {
      count ++;
      if (count>DATA->NumberOfParticlesToInclude) break;
      fprintf(stderr,"%d %e %e %e %e %d %d %d \n", number, deltaY, ymax, slope, phimax, iymax, iptmax, iphimax);
      ip = partid[MHALF+number];
      particleList[ip].ny = iymax;
      particleList[ip].npt = iptmax;
      particleList[ip].nphi = iphimax;
      particleList[ip].phimax = phimax;
      particleList[ip].phimin = phimin;
      particleList[ip].slope = slope;
      particleList[ip].ymax = ymax;
      particleList[ip].deltaY = deltaY;
      for ( i=0; i<iptmax; i++ )
	{
	  particleList[ip].pt[i] =  gala15x[i]/slope;
	}
      for ( i=0; i<iymax; i++ )
	{
	  particleList[ip].y[i] =  i*deltaY-ymax+deltaY/2.;
	  //	  if (ip==1) cout << "read particleList[ip].y[" << i << "] = " <<  particleList[ip].y[i] << endl;
	}

      switch (iphimax) 
	{
	case 4: p= gaulep4; w= gaulew4; break;
	case 8: p= gaulep8; w= gaulew8; break;
	case 10: p= gaulep10; w= gaulew10; break;
	case 12: p= gaulep12; w= gaulew12; break;
	case 16: p= gaulep16; w= gaulew16; break;
	case 20: p= gaulep20; w= gaulew20; break;
	case 48: p= gaulep48; w= gaulew48; break;
	default: fprintf(stderr,"specified number of phi-points not available\n"); exit(1);
	}
      
      phiArray = util->vector_malloc(iphimax);
      for(iphi=0; iphi<iphimax; iphi++)
	{
	  if ( iphi < iphimax/2 )
	    {
	      phiArray[iphi] = (phimax-phimin)/2.*(1.-p[iphi])+phimin;
	    }
	  else
	    {
	      phiArray[iphi] = (phimax-phimin)/2.*(1.+p[iphimax-iphi-1])+phimin;
	    }
	  //      cout << "phi[" << iphi << "]=" << phiArray[iphi] << ", " <<  phimax/2.*(1.-p[iphi]) << ", " 
	  //   << (phimax)/2.*(1.+p[iphimax-iphi-1]) << endl;
	}
    }
  
  particleMax = ip;

  fclose(p_file);

  FILE *s_file;
  const char* s_name = "yptphiSpectra.dat";
  s_file = fopen(s_name, "r");
  checkForReadError(s_file,s_name);
  
  cout << "iymax=" << iymax << endl;
  cout << "cells=" << iymax*iptmax*iphimax << endl;

  //fprintf(stderr,"reading particle %d: %d %s\n", ip, particleList[ip].number, particleList[ip].name);
  for (iy=0; iy<iymax; iy++)
    {
      for (ipt=0; ipt<iptmax; ipt++)
	{
	  for (iphi=0; iphi<iphimax; iphi++)
	    {
	      bytes_read=fscanf(s_file, "%lf", &particleList[ip].dNdydptdphi[iy][ipt][iphi]);
	      //	      cout << particleList[ip].y[iy] << " " << particleList[ip].pt[ipt] << " " << phiArray[iphi] << " " << particleList[ip].dNdydptdphi[iy][ipt][iphi] << endl; 
	      //	  printf("%f %f %f \n",particleList[ip].y[iy],particleList[ip].pt[ipt],phiArray[iphi]);
	    }
	}
    }
  
  //particleMax=2;
  fclose(s_file);
}

double Freeze::summation(double px, double py, double y, double m, int deg, int baryon, double muAtFreezeOut, InitData *DATA)
{
  double sum = 0.;
  int i;
  double ptau, peta;
  double f, T, mu, tau, x, eta, E, sign, delta_f;
  double pdSigma, Wfactor;
  double mt = sqrt(m*m+px*px+py*py); // all in GeV
  double alpha = 0.; // make a parameter
    // Bose or Fermi statistics.
  if (baryon==0)
    sign = -1.;
  else
    sign = 1.;

  //fprintf(stderr,"sign=%f\n",sign);
  //fprintf(stderr,"baryon=%d\n",baryon);
  for (i=0; i<NCells; i++)
    {
      tau = surface[i].x[0];
      eta = surface[i].x[3];
      
      ptau = mt*cosh(y-eta); // GeV    this is p^tau
      peta = mt/tau*sinh(y-eta); // GeV/fm     this is p^eta

      // compute p^mu*dSigma_mu
      pdSigma = tau*(ptau*surface[i].s[0]+px*surface[i].s[1]+py*surface[i].s[2]+peta*surface[i].s[3]); //fm^3*GeV
      //pdSigma = tau*(ptau*surface[i].s[0]-px*surface[i].s[1]-py*surface[i].s[2]-peta*surface[i].s[3]); //fm^3*GeV
     
      // compute f
      T = surface[i].T_f*hbarc; // GeV
      mu = baryon*surface[i].mu_B*hbarc; //GeV
      if(DATA->whichEOS>=3) // for PCE use the previously computed mu at the freeze-out energy density
	mu=muAtFreezeOut; //GeV
      // cout << "mu=" << mu << endl;
      // exit(1);
      //     cout << "mu=" << mu << " GeV " << endl;
      // E = u^mu*p_mu
      E = (ptau*surface[i].u[0]-px*surface[i].u[1]-py*surface[i].u[2]-tau*tau*peta*surface[i].u[3]/tau);
      //E = (ptau*surface[i].u[0]+px*surface[i].u[1]+py*surface[i].u[2]+tau*tau*peta*surface[i].u[3]/tau);
      // this is the equilibrium f, f_0:
      f = 1./(exp(1./T*(E-mu))+sign);
      // now comes the delta_f: check if still correct at finite mu_b 
      // we assume here the same C=eta/s for all particle species because it is the simplest way to do it.
      // also we assume Xi(p)=p^2, the quadratic Ansatz

      if(f<=0)
	f=0;

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
      
 //      if( fabs(surface[i].W[3][3]/surface[i].sFO/(T/hbarc))>1.5) 
// 	{cout << surface[i].W[0][0]/surface[i].sFO/(T/hbarc) << " " 
// 	      << surface[i].W[0][1]/surface[i].sFO/(T/hbarc) << " " 
// 	      << surface[i].W[0][2]/surface[i].sFO/(T/hbarc) << " " 
// 	      << surface[i].W[0][3]/surface[i].sFO/(T/hbarc) << " " 
// 	      << surface[i].W[1][1]/surface[i].sFO/(T/hbarc) << " " 
// 	      << surface[i].W[1][2]/surface[i].sFO/(T/hbarc) << " " 
// 	      << surface[i].W[1][3]/surface[i].sFO/(T/hbarc) << " " 
// 	      << surface[i].W[2][2]/surface[i].sFO/(T/hbarc) << " " 
// 	      << surface[i].W[2][3]/surface[i].sFO/(T/hbarc) << " " 
// 	      << surface[i].W[3][3]/surface[i].sFO/(T/hbarc) << " " 
// 	      << endl;
// 	}
//       cout << ptau*surface[i].W[0][0]*ptau << " " 
// 	   << -2.*ptau*surface[i].W[0][1]*px << " " 
// 	   << -2.*ptau*surface[i].W[0][2]*py << " "
// 	   << -2.*tau*tau*ptau*surface[i].W[0][3]/tau*peta << " " 
// 	   << px*surface[i].W[1][1]*px << " " 
// 	   << 2.*px*surface[i].W[1][2]*py << " " 
// 	   << 2.*tau*tau*px*surface[i].W[1][3]/tau*peta << " " 
// 	   << py*surface[i].W[2][2]*py << " " 
// 	   << 2.*tau*tau*py*surface[i].W[2][3]/tau*peta << " " 
// 	   << tau*tau*tau*tau*peta*surface[i].W[3][3]/tau/tau*peta << " " << endl;
//       cout << Wfactor/pow(hbarc,4.) << endl;
	


    //   if (delta_f>f && f>0.001 && fabs(y)>3.)
// 	{
// 	  cout << "f=" << f << ", delta_f=" << delta_f << ", delta_f/f=" << delta_f/f <<endl;
// 	  cout << "y=" << y << ", px=" << px << ", py=" << py << ", Wfactor=" << Wfactor << endl;
// 	  cout << "f*(1.-sign*f)=" << f*(1.-sign*f) << endl;
// 	  cout << "1/(2.*surface[i].sFO*pow(hbarc,3.)*pow(T,3.))=" <<1/(2.*surface[i].sFO*pow(hbarc,3.)*pow(T,3.)) << endl;
// 	  cout << "surface[i].sFO*pow(hbarc,3.)=" << surface[i].sFO*pow(hbarc,3.) << endl;
// 	  cout << "surface[i].sFO=" << surface[i].sFO << endl;
// 	  cout << "pow(T,3.)=" << pow(T,3.) << endl;
// 	  cout << "T=" << T << endl;

// 	  cout << surface[i].W[0][0] << " " << surface[i].W[0][1] << " " << surface[i].W[0][2] << " " << surface[i].W[0][3] << " " << 
// 	    surface[i].W[1][1] << " " << surface[i].W[1][2] << " " << surface[i].W[1][3] << " " << surface[i].W[2][2] << " " << 
// 	    surface[i].W[2][3] << " " << surface[i].W[3][3] << endl;
	
// 	  cout << surface[i].W[1][1] << " " << surface[i].W[2][2] << " " << surface[i].W[3][3] << endl;
// 	  cout << ptau << " " << px << " " << py << " " << peta << endl;
// 	  cout << "y=" << y << ", eta=" << eta << endl;
// 	}

	sum += 1/pow(2.*PI,3.) * (f+delta_f) * pdSigma;
	//	cout << "sum=" << sum << ", pdSigma=" << pdSigma << endl;
	if (sum>10000)
	  cout << "WARNING: sum>10000 in summation. sum=" << sum << ", f=" << f << ", deltaf=" << delta_f << ", pdSigma=" << pdSigma 
	       << ", T=" << T << ", E=" << E << ", mu=" << mu << endl;
    }

  sum *= deg/pow(hbarc,3.); // in GeV^(-2)

  return sum;
}

void Freeze::ComputeParticleSpectrum(InitData *DATA, int number, double ptmax, int anti, int iptmax, int iphimax, int size, int rank)
{
  char *specString;
  specString = util->char_malloc(30);
  char *numberStringy;
  char *numberString;
  char *numberStringpty;
  char *numberStringPhi;
  char *numberStringeta;
  char *numberStringv2;
  char *numberStringv2eta;
  char *numberStringv4;
  char *numberStringv4eta;
  char buf[10];
  double *p, *w;		        // pointing to data for Gaussian integration in phi 
  double slope, slope1, slope2, fleft1, fleft2, fright1, fright2;
  double ymaxIntegral;
  int ipt, iphi, iymax, iy, iymaxIntegral, ieta, ietamax;
  int i,j;
  double pt, phi, px, py, y, deltaPT, deltaPhi, deltaY, ymax, phimin, phimax, sum, sumpt, sumpt2, sumv[5];
  double phiOffs, phiDiff, deltaEta, etamax;
  int returnValue;
  double eta, mt, etaMaxIntegral;
  double half;
 
  if (iptmax != 15) 
    {
      fprintf(stderr,"only iptmax=15 possible. you picked %d. exiting.\n", iptmax);
      exit(1);
    }

  j = partid[MHALF+number];
  // set some parameters
  deltaY = DATA->deltaY; //use 0.05
  ymax = DATA->ymax;
  etaMaxIntegral = 1.3;
  ymaxIntegral = 0.5-deltaY/2.;
  iymax = floor(2.*ymax/deltaY+0.0001);
  fprintf(stderr,"iymax=%d\n",iymax);
  if (size>0)
    {
      if (iymax%size!=0)
	{
	  fprintf(stderr,"number of steps in y (iymax) is not a multiple of the number of processors. Exiting.\n");
	  MPI::Finalize();
	  exit(1);
	}
      iymax=iymax/size;
      cout << "r" << rank << " iymax=" << iymax << endl;
    }
  iymaxIntegral = floor(2.*ymaxIntegral/deltaY);
  particleList[j].ymax = ymax;
  particleList[j].deltaY = deltaY;
  
  if (DATA->rotateBy45degrees==0)
    {
      phimin = 0.;
      phimax = 2.*PI;
    }
  else
    {
      phimin = PI/4.;
      phimax = 3.*PI/4.; // can select larger phi.
    }
  
  phiArray = util->vector_malloc(iphimax);
  phiOffs = 0.5 * ( phimin + phimax );
  phiDiff = 0.5 * ( phimax - phimin );

  fprintf(stderr,"Doing %d: %s (%d)\n", j, particleList[j].name, particleList[j].number);
 
  particleList[j].ny = iymax*size;
  particleList[j].npt = iptmax;
  particleList[j].nphi = iphimax;
  
  // set particle properties
  double m = particleList[j].mass;
  int d = particleList[j].degeneracy;
  int b = particleList[j].baryon;
  int s = particleList[j].strange;
  double c = particleList[j].charge;
  double mu = particleList[j].muAtFreezeOut;
  //fprintf(stderr,"m=%f \n", m);
  //fprintf(stderr,"d=%d \n", d);
  //fprintf(stderr,"b=%d \n", b);
  sumYPtPhi = util->cube_malloc(iymax, iptmax, iphimax);

  // set up file names
  numberStringy = util->char_malloc(30);
  numberString = util->char_malloc(30);
  numberStringpty = util->char_malloc(30);
  numberStringPhi = util->char_malloc(30);
  numberStringeta = util->char_malloc(30);
  numberStringv2 = util->char_malloc(30);
  numberStringv2eta = util->char_malloc(30);
  numberStringv4 = util->char_malloc(30);
  numberStringv4eta = util->char_malloc(30);
  sprintf (buf, "%d", number);
  
  //itoa(number, buf);
  
  if( chdir("./outputs") != 0 )
    {
      fprintf(stderr,"directory \"outputs\" does not exist. Exiting.\n");
      exit(1);
    }
  else returnValue=chdir("..");

  strcat(numberString, "./outputs/Npt-");
  strcat(numberString, buf);
  strcat(numberString,".dat");

  strcat(numberStringpty, "./outputs/Npteta-");
  strcat(numberStringpty, buf);
  strcat(numberStringpty,".dat");

  strcat(numberStringPhi, "./outputs/NphiPT2-3-");
  strcat(numberStringPhi, buf);
  strcat(numberStringPhi,".dat");

  strcat(numberStringeta, "./outputs/Neta-");
  strcat(numberStringeta, buf);
  strcat(numberStringeta,".dat");

  strcat(numberStringv2eta, "./outputs/v2eta-");
  strcat(numberStringv2eta, buf);
  strcat(numberStringv2eta,".dat");
  
  strcat(numberStringv4eta, "./outputs/v4eta-");
  strcat(numberStringv4eta, buf);
  strcat(numberStringv4eta,".dat");
  
  strcat(numberStringy, "./outputs/Ny-");
  strcat(numberStringy, buf);
  strcat(numberStringy,".dat");
 
  strcat(numberStringv2, "./outputs/v2pt-");
  strcat(numberStringv2, buf);
  strcat(numberStringv2,".dat");

  strcat(numberStringv4, "./outputs/v4pt-");
  strcat(numberStringv4, buf);
  strcat(numberStringv4,".dat");

  // open files to write
  FILE *d_file;
  const char* d_name = "particleInformation.dat";
  d_file = fopen(d_name, "a");

  FILE *s_file;
  sprintf (buf, "%d", rank);
  
  strcat(specString, "yptphiSpectra");
  strcat(specString, buf);
  strcat(specString, ".dat");
  char* s_name = specString;
  s_file = fopen(s_name, "a");

  FILE *y_file;
  char* y_name = numberStringy;
  y_file = fopen(y_name, "w");

  FILE *pty_file;
  char* pty_name = numberStringpty;
  pty_file = fopen(pty_name, "w");

  FILE *eta_file;
  char* eta_name = numberStringeta;
  eta_file = fopen(eta_name, "w");

  FILE *v2_eta_file;
  char* v2_eta_name = numberStringv2eta;
  v2_eta_file = fopen(v2_eta_name, "w");

  FILE *v4_eta_file;
  char* v4_eta_name = numberStringv4eta;
  v4_eta_file = fopen(v4_eta_name, "w");

  FILE *p_file;
  char* p_name = numberString;
  p_file = fopen(p_name, "w");

  FILE *phi_file;
  char* phi_name = numberStringPhi;
  phi_file = fopen(phi_name, "w");

  FILE *v2_file;
  char* v2_name = numberStringv2;
  v2_file = fopen(v2_name, "w");

  FILE *v4_file;
  char* v4_name = numberStringv4;
  v4_file = fopen(v4_name, "w");

  FILE *a_file;
  const char* a_name = "angle";
  a_file = fopen(a_name, "w");

  // set spacing in pt and phi:

  // --------------------------------------------------------------------------
  // slope for the pt-spacing
//   fright1 = summation(ptmax, 0., y, m, d, b, mu, DATA);
//   fleft1  = summation(0.9*ptmax, 0., y, m, d, b, mu, DATA);
//   slope1  = -(log(fright1)-log(fleft1))/(ptmax-ptmax*0.9);
//   fright2 = summation(0., ptmax, y, m, d, b, mu, DATA);
//   fleft2  = summation(0., 0.9*ptmax, y, m, d, b, mu, DATA);
//   slope2  = -(log(fright2)-log(fleft2))/(ptmax-ptmax*0.9);
//   slope = 0.5 * (slope1+slope2); // slope at the high-pt end of the spectrum
//   // --------------------------------------------------------------------------
//   if (slope!=slope)
//     {
//       cout << "[Freeze:ComputePartcileSpectrum]: WARNING: slope not a number." << endl;
//       slope=3.;
//     }

  // fixed slope. makes more sense right now
  slope=3;


  slope = slope*4.;
  //slope = 15.943641;
  particleList[j].slope = slope;
  
  fprintf(stderr,"slope=%f\n",slope);
  // --------------------------------------------------------------------------
  // define phi spacings
  switch (iphimax) 
    {
    case 4: p= gaulep4; w= gaulew4; break;
    case 8: p= gaulep8; w= gaulew8; break;
    case 10: p= gaulep10; w= gaulew10; break;
    case 12: p= gaulep12; w= gaulew12; break;
    case 16: p= gaulep16; w= gaulew16; break;
    case 20: p= gaulep20; w= gaulew20; break;
    case 48: p= gaulep48; w= gaulew48; break;
    default: fprintf(stderr,"specified number of phi-points not available\n"); exit(1);
    }

  //  deltaPhi=(phimax-phimin)/iphimax;
  for(iphi=0; iphi<iphimax; iphi++)
    {
      if ( iphi < iphimax/2 )
	{
	  phiArray[iphi] = (phimax-phimin)/2.*(1.-p[iphi])+phimin;
	}
      else
	{
	  phiArray[iphi] = (phimax-phimin)/2.*(1.+p[iphimax-iphi-1])+phimin;
	}
      //            cout << "phi[" << iphi << "]=" << phiArray[iphi] << endl;
    }
  // --------------------------------------------------------------------------
  
  // write information in the particle information file
  if (rank==0) fprintf(d_file,"%d %e %e %e %e %e %d %d %d \n", number, deltaY, ymax, slope, phimin, phimax, iymax*size, iptmax, iphimax);

  // store value as function of phi, pt and y in sumPtPhi:
  
  for (iy=0; iy<iymax; iy++)
    {
      y = iy*deltaY-ymax+rank*(ymax/size*2.)+deltaY/2.;
      particleList[j].y[iy] = y;
      //     if (j==1) cout << " do particleList[ip].y[" << iy << "] = " <<  particleList[j].y[iy] << " y = " << y << endl;
      sumpt=0.;
      for (ipt=0; ipt<iptmax; ipt++)
	{
	  pt = gala15x[ipt]/slope; // legendre abscissa / slope
	  particleList[j].pt[ipt] = pt;
	  for (iphi=0; iphi<iphimax; iphi++)
	    {
	      phi = phiArray[iphi];
	      px = pt*cos(phi);
	      py = pt*sin(phi);
	      //cout << "phi=" << phi << endl;
	      sum = summation3(px, py, y, m, d, b, mu, DATA);
	      sumYPtPhi[iy][ipt][iphi] = sum;
	      particleList[j].dNdydptdphi[iy][ipt][iphi] = sum;
	      fprintf(s_file,"%e ", sum);
	      //	      sumpt += deltaPhi * ( sum ) * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
	      if (iphi<iphimax/2)
		{
		  sumpt += w[iphi] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		}
	      else
		{
		  sumpt += w[iphimax-iphi-1] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		}
	    }
	  fprintf(s_file,"\n");
	}
      
      fprintf(stderr,"r=%d,       y=%f, sumpt=%f\n",rank, y, sumpt/slope);
      fprintf(y_file,"%f %f\n", y, sumpt/slope);
    }



  // uncomment this for correlations

//   //as function of phi for 'correlations'
//   deltaEta=deltaY;
//   double sumPhi;
//   etamax=ymax;
//   ietamax=floor(2.*etamax/deltaEta);
//   double assocLower = 2.;
//   double assocHigher = 3.;
//   double ptRange = assocHigher-assocLower;
//   double MydeltaPT;
//   double MydeltaPhi;
//   int steps = 40;
//   MydeltaPT = ptRange/steps;
//   MydeltaPhi = 2*PI/iphimax;
//   for (iphi=0; iphi<iphimax; iphi++)
//     {
//       sumPhi=0.;
//       phi = MydeltaPhi*iphi;
//       for (ipt=0; ipt<steps; ipt++)
// 	{
// 	  pt = assocLower+ipt*MydeltaPT;
// 	  y = 0.;
// 	  cout << "phi=" << phi << ", pt=" << pt << endl;
// 	  iy = floor((y+ymax)/deltaY+0.0001);

// 	  //cout << "y=" << iy*deltaY-ymax << endl;
	  
// 	  px = pt*cos(phi);
// 	  py = pt*sin(phi);
// 	  mt = sqrt(pt*pt+m*m);
// 	  sum = summation(px, py, y, m, d, b)*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
// 	  if(ipt==0 || ipt==steps)
// 	    sumPhi += 0.5 * sum * pt * MydeltaPT *0.7 ; // 0.7 for eta range -0.35 < eta < 0.35
// 	  else
// 	    sumPhi += sum * pt * MydeltaPT *0.7 ; // 0.7 for eta range -0.35 < eta < 0.35
// 	}
      
//       fprintf(phi_file,"%e %e\n", phi, sumPhi);
//       cout << phi << " " << sumPhi << endl;
//     }

// ---- correlations end



 

//    for (ipt=0; ipt<iptmax; ipt++)
//     {
//       int countY = 0;
//       sumpt = 0.;
//       pt = gala15x[ipt]/slope; // legendre abscissa / slope
//       for (iy=0; iy<iymax; iy++)
// 	{
// 	  y = iy*deltaY-ymax+deltaY/2.;
// 	  if (fabs(y)>ymaxIntegral) continue;
// 	  countY++;
// 	  fprintf(stderr,"y=%f\n",y);
// 	  // fprintf(stderr,"countY=%d\n",countY);
// 	  // fprintf(stderr,"iymaxIntegral=%d\n",iymaxIntegral);
// 	  for (iphi=0; iphi<iphimax; iphi++)
// 	    {
// 	      phi = phiArray[iphi];
// 	      sum = sumYPtPhi[iy][ipt][iphi];
// 	      if (iphi<iphimax/2)
// 		{
// 		  if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
// 		    {
// 		      sumpt += 0.5 * w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/phimax;
// 		    }
// 		  else
// 		    {
// 		      sumpt += w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/phimax;
// 		    }
// 		  //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphi], iphi, p[iphi]);
// 		}
// 	      else
// 		{
// 		  if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
// 		    {
// 		      sumpt += 0.5 * w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/phimax ;
// 		    }
// 		  else
// 		    {
// 		      sumpt += w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/phimax ;
// 		    }
// 		  //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphimax-iphi-1], iphi, -p[iphimax-iphi-1]);
// 		}
// 	      //if (iphi==iphimax-1) fprintf(stderr,"value for y=%f: %f\n",y, sumpt);
// 	    }
// 	}
//       fprintf(stderr,"pt=%f, sumpt=%f\n", pt, sumpt / (2.*ymaxIntegral));
//       fprintf(p_file,"%e %e \n", pt, sumpt / (2.*ymaxIntegral));
//     }

//   // v2 and v4
//   for (ipt=0; ipt<iptmax; ipt++)
//     {
//       sumpt = 0.;
//       sumv[2] = 0.;
//       sumv[4] = 0.;
//       int countY = 0;
//       pt = particleList[j].pt[ipt];
//       ymaxIntegral = 0.5*log((sqrt(pt*pt*cosh(etaMaxIntegral)*cosh(etaMaxIntegral)+m*m)
// 			      +pt*sinh(etaMaxIntegral))/(sqrt(pt*pt*cosh(etaMaxIntegral)*cosh(etaMaxIntegral)+m*m)-pt*sinh(etaMaxIntegral)));
//       iymaxIntegral = floor(2.*ymaxIntegral/deltaY);
//       for (iy=0; iy<iymax; iy++)
// 	{
// 	  y = iy*deltaY-ymax;
// 	  if (fabs(y)>ymaxIntegral) continue;
// 	  countY++;
// 	  for (iphi=0; iphi<iphimax; iphi++)
// 	    {
// 	      phi = phiArray[iphi];
// 	      sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
// 	      if (iphi<iphimax/2)
// 		{
// 		  if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
// 		    {
// 		      sumpt += 0.5 * w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/phimax;
// 		      sumv[2] += 0.5 * w[iphi] * ( sum*cos(2.*phi) ) * deltaY * phiDiff * 2.*PI/phimax;
// 		      sumv[4] += 0.5 * w[iphi] * ( sum*cos(4.*phi) ) * deltaY * phiDiff * 2.*PI/phimax;
// 		    }
// 		  else
// 		    {
// 		      sumpt += w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/phimax;
// 		      sumv[2] += w[iphi] * ( sum*cos(2.*phi) ) * deltaY * phiDiff * 2.*PI/phimax;
// 		      sumv[4] += w[iphi] * ( sum*cos(4.*phi) ) * deltaY * phiDiff * 2.*PI/phimax;
// 		    }
// 		  //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphi], iphi, p[iphi]);
// 		}
// 	      else
// 		{
// 		  if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
// 		    {
// 		      sumpt += 0.5 * w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/phimax ;
// 		      sumv[2] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(2.*phi) ) * deltaY * phiDiff * 2.*PI/phimax;
// 		      sumv[4] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(4.*phi) ) * deltaY * phiDiff * 2.*PI/phimax;
// 		    }
// 		  else
// 		    {
// 		      sumpt += w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/phimax ;
// 		      sumv[2] += w[iphimax-iphi-1] * ( sum*cos(2.*phi) ) * deltaY * phiDiff * 2.*PI/phimax;
// 		      sumv[4] += w[iphimax-iphi-1] * ( sum*cos(4.*phi) ) * deltaY * phiDiff * 2.*PI/phimax;
// 		    }
// 		  //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphimax-iphi-1], iphi, -p[iphimax-iphi-1]);
// 		}
// 	    }
// 	}
//       //fprintf(stderr,"pt=%f, sumpt=%f\n", pt, sumpt / (2.*ymaxIntegral));
//       //fprintf(stderr,"pt=%f, sumv%d=%f\n", pt, 2, sumv[2]/sumpt);
//       fprintf(v2_file,"%e %e \n", pt, sumv[2]/sumpt);
//       fprintf(v4_file,"%e %e \n", pt, sumv[4]/sumpt);
//     }



//   //pseudorapidity2
//   deltaEta=deltaY;
//   etamax=ymax;
//   ietamax=floor(2.*etamax/deltaEta);
//   for (ieta=0; ieta<ietamax; ieta++)
//     {
//       eta = ieta*deltaEta-etamax+deltaEta/2.;
//       sumpt=0.;
//       sumv[2]=0.;
//       for (ipt=0; ipt<iptmax; ipt++)
// 	{
// 	  pt = particleList[j].pt[ipt];
// 	  y = 0.5*log((sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)+pt*sinh(eta))/(sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)-pt*sinh(eta)));
	 
// 	  for(i=0; i<iymax-1; i++) // find closest iy
// 	    {
// 	      if ( particleList[j].y[i]< y &&  particleList[j].y[i+1] > y )
// 		{
// 		  if ( fabs(particleList[j].y[i]-y) <  fabs(particleList[j].y[i+1]-y)) 
// 		    iy = i;
// 		  else
// 		    iy = i+1;
// 		}
// 	    }

// 	  mt = sqrt(pt*pt+m*m);
// 	  for (iphi=0; iphi<iphimax; iphi++)
// 	    {
// 	      phi = phiArray[iphi];
// 	      px = pt*cos(phi);
// 	      py = pt*sin(phi);
// 	      sum = particleList[j].dNdydptdphi[iy][ipt][iphi]*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
// 	      if (iphi<iphimax/2)
// 		{
// 		  sumpt += w[iphi] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/phimax;
// 		  sumv[2] += w[iphi] * ( sum*cos(2.*phi) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/phimax;
// 		  sumv[4] += w[iphi] * ( sum*cos(4.*phi) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/phimax;
// 		}
// 	      else
// 		{
// 		  sumpt += w[iphimax-iphi-1] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/phimax;
// 		  sumv[2] += w[iphimax-iphi-1] * ( sum*cos(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/phimax;
// 		  sumv[4] += w[iphimax-iphi-1] * ( sum*cos(4.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/phimax;
// 		}
// 	    }
// 	}
//       //   fprintf(stderr,"eta=%f, sumpt=%f\n", eta, sumpt/slope);
//       fprintf(eta_file,"%f %f\n", eta, sumpt/slope);
//       fprintf(v2_eta_file,"%f %f\n", eta, sumv[2]/sumpt);
//       fprintf(v4_eta_file,"%f %f\n", eta, sumv[4]/sumpt);
//     }
  

//   //angle
//   ipt = floor((double)(iptmax)/2.);
//   pt = particleList[j].pt[ipt];
//   iy = floor((double)(iymax)/2.);
//   y = iy*deltaY-ymax;
//   for (iphi=0; iphi<iphimax; iphi++)
//     {
//       phi = phiArray[iphi];
//       sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
//       fprintf(a_file,"%e %e %e %e %e %e %e %e \n", pt, y, phi, particleList[j].dNdydptdphi[iy-4][ipt][iphi], particleList[j].dNdydptdphi[iy-2][ipt][iphi], sum,particleList[j].dNdydptdphi[iy+2][ipt][iphi],particleList[j].dNdydptdphi[iy+4][ipt][iphi] );
//     }



/*   //trapezoid */
/*   fprintf(y_file,"\n"); */
/*   deltaPhi = 0.1; */
/*   iphimax = floor(phimax/deltaPhi); */
/*   deltaPT = 0.1; */
/*   iptmax = floor(5./deltaPT); */
  
/*   for (iy=0; iy<iymax; iy++) */
/*     { */
/*       y = iy*deltaY-ymax+deltaY/2.; */
/*       sumpt=0.; */
/*       for (ipt=0; ipt<iptmax; ipt++) */
/* 	{ */
/* 	  pt = ipt*deltaPT; */
/* 	  for (iphi=0; iphi<iphimax; iphi++) */
/* 	    { */
/* 	      phi = iphi*deltaPhi; */
/* 	      px = pt*cos(phi); */
/* 	      py = pt*sin(phi); */
/* 	      sum = summation(px, py, y, m, d, b); */
/* 	      if (ipt==0 || ipt==iptmax-1) // trapezoid rule: only half the edges */
/* 		{ */
/* 		  if (iphi==0 || iphi==iphimax-1) // trapezoid rule: only half the edges */
/* 		    { */
/* 		      half = 0.25; */
/* 		    } */
/* 		  else */
/* 		    { */
/* 		      half = 0.5; */
/* 		    } */
/* 		} */
/* 	      else */
/* 		{ */
/* 		  if (iphi==0 || iphi==iphimax-1) // trapezoid rule: only half the edges */
/* 		    { */
/* 		      half = 0.5; */
/* 		    } */
/* 		  else */
/* 		    { */
/* 		      half = 1.; */
/* 		    } */
/* 		} */
/* 	      sumpt += half * ( sum ) * deltaPhi * deltaPT * pt * 2.*PI/phimax; */
/* 	    } */
/* 	} */
/*       fprintf(stderr,"y=%f, sumpt=%f\n", y, sumpt); */
/*       fprintf(y_file,"%f %f\n", y, sumpt); */
/*     } */
/*   exit(1); */

   
  // simple integration for check
/*   for (ipt=0; ipt<2.*iptmax; ipt++) */
/*     {    */
/*       deltaPT = 0.05; */
/*       sumpt = 0.; */
/*       sumv[2] = 0.; */
/*       for (iy=0; iy<iymax; iy++) */
/* 	{ */
/* 	  y = iy*deltaY-ymax; */
/* 	  pt = ipt*deltaPT; // legendre abscissa / slope */
/* 	  deltaPhi=0.01; */
/* 	  iphimax = floor(phimax/deltaPhi); */
/* 	  for (iphi=0; iphi<iphimax; iphi++) */
/* 	    { */
/* 	      phi = iphi*deltaPhi; */
/* 	      px = pt*cos(phi); */
/* 	      py = pt*sin(phi); */
/* 	      sum = summation(px, py, y, m, d, b); */
/* 	      sumpt += ( sum ) * deltaY * deltaPhi * 2.*PI/phimax; */
/* 	      sumv[2] += ( sum*cos(2.*phi) ) * deltaY * deltaPhi * 2.*PI/phimax; */
/* 	    } */
/* 	} */
/*       fprintf(stderr,"pt=%f, sumpt=%f\n", pt, sumpt / (2.*ymax)); */
/*       fprintf(stderr,"pt=%f, sumv%d=%f\n", pt, 2, sumv[2]/sumpt); */
/*     } */

  fclose(a_file);
  fclose(d_file);
  fclose(s_file);
  fclose(p_file);
  fclose(phi_file);
  fclose(y_file);
  fclose(eta_file);
  fclose(v2_file);
  fclose(v2_eta_file);
  fclose(v4_file);
  fclose(v4_eta_file);
  //  util->vector_free(phiArray);
  //util->cube_free(sumYPtPhi,iptmax,iphimax,iymax);
}


  void Freeze::OutputFullParticleSpectrum(InitData *DATA, int number, double ptmax, int anti, int full)
{
  char *numberStringy;
  char *numberStringpteta;
  char *numberString;
  char *numberStringy0;
  char *numberStringPhi;
  char *numberStringv2;
  char *numberStringv2y0;
  char *numberStringv2eta;
  char *numberStringv2pteta;
  char *numberStringv2tot;
  char *numberStringv3;
  char *numberStringv3r3;
  char *numberStringv3eta;
  char *numberStringv3pteta;
  char *numberStringv3tot;
  char *numberStringv3r3eta;
  char *numberStringv3r3tot;
  char *numberStringv4;
  char *numberStringv4eta;
  char *numberStringv4pteta;
  char *numberStringv4tot;
  char *numberStringeta;

  char *numberStringpt;
  char *numberStringv1pt;
  char *numberStringv2pt;
  char *numberStringv3pt;
  char *numberStringv4pt;
  char *numberStringv5pt;
  char *numberStringv6pt;

  FILE *s_file;
  FILE *d_file;
  char buf[10];
  double *p, *w;		        // pointing to data for Gaussian integration in phi 
  double slope, slope1, slope2, fleft1, fleft2, fright1, fright2;
  int ipt, iphi, iymax, iy, iptmax, iphimax, ieta, ietamax;
  int i,j, iymaxIntegral;
  double pt, phi, px, py, y, deltaPT, deltaPhi, deltaY, ymax, phimin, phimax, sum, sumpt, sumpt2, sumv[7], sumv3r3, sumPhi, sumpteta, sumvpteta[7];
  double phiOffs, phiDiff, ymaxIntegral, eta, deltaEta, etamax, mt, etaMaxIntegral;
  int returnValue;

  // set some parameters
  
  j = partid[MHALF+number];
  fprintf(stderr,"Doing %s (%d)\n", particleList[j].name, particleList[j].number);
  phimax = particleList[j].phimax;
  phimin = particleList[j].phimin;
  phiOffs = 0.5 * ( phimin + phimax );
  phiDiff = 0.5 * ( phimax - phimin );

  ymax = particleList[j].ymax;
  deltaY = particleList[j].deltaY;
  ymaxIntegral = 0.5;
  etaMaxIntegral = 1.3; // for v2 integration
  iymax = particleList[j].ny;
  iymaxIntegral = floor(2.*ymaxIntegral/deltaY+0.00000000000001);
  iptmax = particleList[j].npt;
  iphimax = particleList[j].nphi;

  if (iptmax != 15) 
    {
      fprintf(stderr,"only iptmax=15 possible. you picked %d. fixing...\n", iptmax);
      fprintf(stderr,"iphimax[%d]=%d, phimax[%d]=%f, phimin[%d]=%f,  \n", j, iphimax, j, phimax, j, phimin);
      fprintf(stderr,"ymax[%d]=%f, deltaY[%d]=%f, iymax[%d]=%d \n", j, ymax, j, deltaY, j, iymax);
      //      exit(1);
      iptmax = 15;
    }

  
  // set particle properties
  double m = particleList[j].mass;
  int d = particleList[j].degeneracy;
  int b = particleList[j].baryon;
  int s = particleList[j].strange;
  double c = particleList[j].charge;
  double mu = particleList[j].muAtFreezeOut;
  //fprintf(stderr,"m=%f \n", m);
  //fprintf(stderr,"d=%d \n", d);
  //fprintf(stderr,"b=%d \n", b);
  //fprintf(stderr,"iptmax=%d \n",iptmax );
  //fprintf(stderr,"phimax=%f \n", phimax);
  //fprintf(stderr,"iymax=%d \n", iymax);

  // set up file names
  numberStringy = util->char_malloc(30);
  numberString = util->char_malloc(30);
  numberStringy0 = util->char_malloc(30);
  numberStringpteta = util->char_malloc(30);
  numberStringPhi = util->char_malloc(30);
  numberStringeta = util->char_malloc(30);
  numberStringv2 = util->char_malloc(30);
  numberStringv2y0 = util->char_malloc(30);
  numberStringv3 = util->char_malloc(30);
  numberStringv3r3 = util->char_malloc(30);
  numberStringv4 = util->char_malloc(30);
  numberStringv2eta = util->char_malloc(30);
  numberStringv3eta = util->char_malloc(30);
  numberStringv3r3eta = util->char_malloc(30);
  numberStringv4eta = util->char_malloc(30);
  numberStringv2pteta = util->char_malloc(30);
  numberStringv3pteta = util->char_malloc(30);
  numberStringv4pteta = util->char_malloc(30);
  numberStringv2tot = util->char_malloc(30);
  numberStringv3tot = util->char_malloc(30);
  numberStringv3r3tot = util->char_malloc(30);
  numberStringv4tot = util->char_malloc(30);
  
  numberStringpt = util->char_malloc(30);
  numberStringv1pt = util->char_malloc(30);
  numberStringv2pt = util->char_malloc(30);
  numberStringv3pt = util->char_malloc(30);
  numberStringv4pt = util->char_malloc(30);
  numberStringv5pt = util->char_malloc(30);
  numberStringv6pt = util->char_malloc(30);
  
  if( chdir("./outputs") != 0 )
    {
      fprintf(stderr,"directory \"outputs\" does not exist. Exiting.\n");
      exit(1);
    }
  else returnValue=chdir("..");

  if (full==1)
    {
      sprintf (buf, "%d", number);
      //itoa(number, buf);
      strcat(numberString, "./outputs/FNpt-");
      strcat(numberString, buf);
      strcat(numberString,".dat");

      strcat(numberStringy0, "./outputs/FNpty0-");
      strcat(numberStringy0, buf);
      strcat(numberStringy0,".dat");

      strcat(numberStringpteta, "./outputs/FNetapt-");
      strcat(numberStringpteta, buf);
      strcat(numberStringpteta,".dat");
      
      strcat(numberStringPhi, "./outputs/FNphi-");
      strcat(numberStringPhi, buf);
      strcat(numberStringPhi,".dat");
      
      strcat(numberStringy, "./outputs/FNy-");
      strcat(numberStringy, buf);
      strcat(numberStringy,".dat");
      
      strcat(numberStringv2, "./outputs/Fv2pt-");
      strcat(numberStringv2, buf);
      strcat(numberStringv2,".dat");
      
      strcat(numberStringv2y0, "./outputs/Fv2pty0-");
      strcat(numberStringv2y0, buf);
      strcat(numberStringv2y0,".dat");
      
      strcat(numberStringv3, "./outputs/Fv3pt-");
      strcat(numberStringv3, buf);
      strcat(numberStringv3,".dat");
      
      strcat(numberStringv3r3, "./outputs/Fv3r3pt-");
      strcat(numberStringv3r3, buf);
      strcat(numberStringv3r3,".dat");
      
      strcat(numberStringv4, "./outputs/Fv4pt-");
      strcat(numberStringv4, buf);
      strcat(numberStringv4,".dat");
      
      strcat(numberStringeta, "./outputs/FNeta-");
      strcat(numberStringeta, buf);
      strcat(numberStringeta,".dat");
      
      strcat(numberStringv2eta, "./outputs/Fv2eta-");
      strcat(numberStringv2eta, buf);
      strcat(numberStringv2eta,".dat");
      
      strcat(numberStringv3eta, "./outputs/Fv3eta-");
      strcat(numberStringv3eta, buf);
      strcat(numberStringv3eta,".dat");
      
      strcat(numberStringv3r3eta, "./outputs/Fv3r3eta-");
      strcat(numberStringv3r3eta, buf);
      strcat(numberStringv3r3eta,".dat");
      
      strcat(numberStringv4eta, "./outputs/Fv4eta-");
      strcat(numberStringv4eta, buf);
      strcat(numberStringv4eta,".dat");
      
      strcat(numberStringv2pteta, "./outputs/Fv2etapt-");
      strcat(numberStringv2pteta, buf);
      strcat(numberStringv2pteta,".dat");
      
      strcat(numberStringv3pteta, "./outputs/Fv3etapt-");
      strcat(numberStringv3pteta, buf);
      strcat(numberStringv3pteta,".dat");
      
      strcat(numberStringv4pteta, "./outputs/Fv4etapt-");
      strcat(numberStringv4pteta, buf);
      strcat(numberStringv4pteta,".dat");
      
      strcat(numberStringv2tot, "./outputs/Fv2tot-");
      strcat(numberStringv2tot, buf);
      strcat(numberStringv2tot,".dat");
      
      strcat(numberStringv3tot, "./outputs/Fv3tot-");
      strcat(numberStringv3tot, buf);
      strcat(numberStringv3tot,".dat");
      
      strcat(numberStringv3r3tot, "./outputs/Fv3r3tot-");
      strcat(numberStringv3r3tot, buf);
      strcat(numberStringv3r3tot,".dat");
      
      strcat(numberStringv4tot, "./outputs/Fv4tot-");
      strcat(numberStringv4tot, buf);
      strcat(numberStringv4tot,".dat");


      strcat(numberStringpt, "./outputs/FNptExp-");
      strcat(numberStringpt, buf);
      strcat(numberStringpt,".dat");
      
      strcat(numberStringv1pt, "./outputs/Fv1ptExp-");
      strcat(numberStringv1pt, buf);
      strcat(numberStringv1pt,".dat");
      
      strcat(numberStringv2pt, "./outputs/Fv2ptExp-");
      strcat(numberStringv2pt, buf);
      strcat(numberStringv2pt,".dat");
      
      strcat(numberStringv3pt, "./outputs/Fv3ptExp-");
      strcat(numberStringv3pt, buf);
      strcat(numberStringv3pt,".dat");
      
      strcat(numberStringv4pt, "./outputs/Fv4ptExp-");
      strcat(numberStringv4pt, buf);
      strcat(numberStringv4pt,".dat");
      
      strcat(numberStringv5pt, "./outputs/Fv5ptExp-");
      strcat(numberStringv5pt, buf);
      strcat(numberStringv5pt,".dat");
      
      strcat(numberStringv6pt, "./outputs/Fv6ptExp-");
      strcat(numberStringv6pt, buf);
      strcat(numberStringv6pt,".dat");
      
      // open files to write
      const char* d_name = "FparticleInformation.dat";
      d_file = fopen(d_name, "a");
      
      const char* s_name = "FyptphiSpectra.dat";
      s_file = fopen(s_name, "a");
    }
  else if (full==0)
    {
      sprintf (buf, "%d", number);
      //itoa(number, buf);
      strcat(numberString, "./outputs/Npt-");
      strcat(numberString, buf);
      strcat(numberString,".dat");
      
      strcat(numberStringy0, "./outputs/Npty0-");
      strcat(numberStringy0, buf);
      strcat(numberStringy0,".dat");
      
      strcat(numberStringpteta, "./outputs/Netapt-");
      strcat(numberStringpteta, buf);
      strcat(numberStringpteta,".dat");
      
      strcat(numberStringPhi, "./outputs/Nphi-");
      strcat(numberStringPhi, buf);
      strcat(numberStringPhi,".dat");
      
      strcat(numberStringy, "./outputs/Ny-");
      strcat(numberStringy, buf);
      strcat(numberStringy,".dat");
      
      strcat(numberStringv2, "./outputs/v2pt-");
      strcat(numberStringv2, buf);
      strcat(numberStringv2,".dat");
      
      strcat(numberStringv2y0, "./outputs/v2pty0-");
      strcat(numberStringv2y0, buf);
      strcat(numberStringv2y0,".dat");
      
      strcat(numberStringv3, "./outputs/v3pt-");
      strcat(numberStringv3, buf);
      strcat(numberStringv3,".dat");
      
      strcat(numberStringv3r3, "./outputs/v3r3pt-");
      strcat(numberStringv3r3, buf);
      strcat(numberStringv3r3,".dat");
      
      strcat(numberStringv4, "./outputs/v4pt-");
      strcat(numberStringv4, buf);
      strcat(numberStringv4,".dat");
      
      strcat(numberStringv2pteta, "./outputs/v2etapt-");
      strcat(numberStringv2pteta, buf);
      strcat(numberStringv2pteta,".dat");
      
      strcat(numberStringv3pteta, "./outputs/v3etapt-");
      strcat(numberStringv3pteta, buf);
      strcat(numberStringv3pteta,".dat");
      
      strcat(numberStringv4pteta, "./outputs/v4etapt-");
      strcat(numberStringv4pteta, buf);
      strcat(numberStringv4pteta,".dat");
      
      strcat(numberStringeta, "./outputs/Neta-");
      strcat(numberStringeta, buf);
      strcat(numberStringeta,".dat");
      
      strcat(numberStringv2eta, "./outputs/v2eta-");
      strcat(numberStringv2eta, buf);
      strcat(numberStringv2eta,".dat");
      
      strcat(numberStringv3eta, "./outputs/v3eta-");
      strcat(numberStringv3eta, buf);
      strcat(numberStringv3eta,".dat");
      
      strcat(numberStringv3r3eta, "./outputs/v3r3eta-");
      strcat(numberStringv3r3eta, buf);
      strcat(numberStringv3r3eta,".dat");
      
      strcat(numberStringv4eta, "./outputs/v4eta-");
      strcat(numberStringv4eta, buf);
      strcat(numberStringv4eta,".dat");

      strcat(numberStringv2tot, "./outputs/v2tot-");
      strcat(numberStringv2tot, buf);
      strcat(numberStringv2tot,".dat");
      
      strcat(numberStringv3tot, "./outputs/v3tot-");
      strcat(numberStringv3tot, buf);
      strcat(numberStringv3tot,".dat");
      
      strcat(numberStringv3r3tot, "./outputs/v3r3tot-");
      strcat(numberStringv3r3tot, buf);
      strcat(numberStringv3r3tot,".dat");
      
      strcat(numberStringv4tot, "./outputs/v4tot-");
      strcat(numberStringv4tot, buf);
      strcat(numberStringv4tot,".dat");

      strcat(numberStringpt, "./outputs/NptExp-");
      strcat(numberStringpt, buf);
      strcat(numberStringpt,".dat");
      
      strcat(numberStringv1pt, "./outputs/v1ptExp-");
      strcat(numberStringv1pt, buf);
      strcat(numberStringv1pt,".dat");
      
      strcat(numberStringv2pt, "./outputs/v2ptExp-");
      strcat(numberStringv2pt, buf);
      strcat(numberStringv2pt,".dat");
      
      strcat(numberStringv3pt, "./outputs/v3ptExp-");
      strcat(numberStringv3pt, buf);
      strcat(numberStringv3pt,".dat");
      
      strcat(numberStringv4pt, "./outputs/v4ptExp-");
      strcat(numberStringv4pt, buf);
      strcat(numberStringv4pt,".dat");
      
      strcat(numberStringv5pt, "./outputs/v5ptExp-");
      strcat(numberStringv5pt, buf);
      strcat(numberStringv5pt,".dat");
      
      strcat(numberStringv6pt, "./outputs/v6ptExp-");
      strcat(numberStringv6pt, buf);
      strcat(numberStringv6pt,".dat");
    }
  else
    fprintf(stderr,"[ERROR in Freeze]: Wrong option 'full=%d' in OutputFullParticleSpectrum\n",full);
  
  FILE *y_file;
  char* y_name = numberStringy;
  y_file = fopen(y_name, "w");
  
  FILE *eta_file;
  char* eta_name = numberStringeta;
  eta_file = fopen(eta_name, "w");
  
  FILE *pteta_file;
  char* pteta_name = numberStringpteta;
  pteta_file = fopen(pteta_name, "w");
  
  FILE *v2_eta_file;
  char* v2_eta_name = numberStringv2eta;
  v2_eta_file = fopen(v2_eta_name, "w");

  FILE *v3_eta_file;
  char* v3_eta_name = numberStringv3eta;
  v3_eta_file = fopen(v3_eta_name, "w");

  FILE *v3r3_eta_file;
  char* v3r3_eta_name = numberStringv3r3eta;
  v3r3_eta_file = fopen(v3r3_eta_name, "w");

  FILE *v4_eta_file;
  char* v4_eta_name = numberStringv4eta;
  v4_eta_file = fopen(v4_eta_name, "w");

  FILE *v2_pteta_file;
  char* v2_pteta_name = numberStringv2pteta;
  v2_pteta_file = fopen(v2_pteta_name, "w");

  FILE *v3_pteta_file;
  char* v3_pteta_name = numberStringv3pteta;
  v3_pteta_file = fopen(v3_pteta_name, "w");

  FILE *v4_pteta_file;
  char* v4_pteta_name = numberStringv4pteta;
  v4_pteta_file = fopen(v4_pteta_name, "w");

  FILE *v2_tot_file;
  char* v2_tot_name = numberStringv2tot;
  v2_tot_file = fopen(v2_tot_name, "w");

  FILE *v3_tot_file;
  char* v3_tot_name = numberStringv3tot;
  v3_tot_file = fopen(v3_tot_name, "w");

  FILE *v3r3_tot_file;
  char* v3r3_tot_name = numberStringv3r3tot;
  v3r3_tot_file = fopen(v3r3_tot_name, "w");

  FILE *v4_tot_file;
  char* v4_tot_name = numberStringv4tot;
  v4_tot_file = fopen(v4_tot_name, "w");

  FILE *p_file;
  char* p_name = numberString;
  p_file = fopen(p_name, "w");

  FILE *py0_file;
  char* py0_name = numberStringy0;
  py0_file = fopen(py0_name, "w");

  FILE *phi_file;
  char* phi_name = numberStringPhi;
  phi_file = fopen(phi_name, "w");

  FILE *v2_file;
  char* v2_name = numberStringv2;
  v2_file = fopen(v2_name, "w");

  FILE *v2y0_file;
  char* v2y0_name = numberStringv2y0;
  v2y0_file = fopen(v2y0_name, "w");

  FILE *v3_file;
  char* v3_name = numberStringv3;
  v3_file = fopen(v3_name, "w");

  FILE *v3r3_file;
  char* v3r3_name = numberStringv3r3;
  v3r3_file = fopen(v3r3_name, "w");

  FILE *v4_file;
  char* v4_name = numberStringv4;
  v4_file = fopen(v4_name, "w");

  FILE *pt_file;
  char* pt_name = numberStringpt;
  pt_file = fopen(pt_name, "w");
  
  FILE *v1_pt_file;
  char* v1_pt_name = numberStringv1pt;
  v1_pt_file = fopen(v1_pt_name, "w");
  
  FILE *v2_pt_file;
  char* v2_pt_name = numberStringv2pt;
  v2_pt_file = fopen(v2_pt_name, "w");
  
  FILE *v3_pt_file;
  char* v3_pt_name = numberStringv3pt;
  v3_pt_file = fopen(v3_pt_name, "w");
  
  FILE *v4_pt_file;
  char* v4_pt_name = numberStringv4pt;
  v4_pt_file = fopen(v4_pt_name, "w");
  
  FILE *v5_pt_file;
  char* v5_pt_name = numberStringv5pt;
  v5_pt_file = fopen(v5_pt_name, "w");
  
  FILE *v6_pt_file;
  char* v6_pt_name = numberStringv6pt;
  v6_pt_file = fopen(v6_pt_name, "w");


  slope = particleList[j].slope;
  switch (iphimax) 
    {
    case 4: p= gaulep4; w= gaulew4; break;
    case 8: p= gaulep8; w= gaulew8; break;
    case 10: p= gaulep10; w= gaulew10; break;
    case 12: p= gaulep12; w= gaulew12; break;
    case 16: p= gaulep16; w= gaulew16; break;
    case 20: p= gaulep20; w= gaulew20; break;
    case 48: p= gaulep48; w= gaulew48; break;
    default: fprintf(stderr,"specified number of phi-points not available\n"); exit(1);
    }

  cout << "phimin=" << phimin << endl;
  cout << "phimax=" << phimax << endl;



//   deltaPhi=(phimax-phimin)/iphimax;
//   for(iphi=0; iphi<iphimax; iphi++)
//     {
//       phiArray[iphi] = iphi*(phimax-phimin)/iphimax;
//       cout << "phi[" << iphi << "]=" << phiArray[iphi] << endl;
//     }


  for(iphi=0; iphi<iphimax; iphi++)
    {
      if ( iphi < iphimax/2 )
	{
	  phiArray[iphi] = (phimax-phimin)/2.*(1.-p[iphi])+phimin;
	}
      else
	{
	  phiArray[iphi] = (phimax-phimin)/2.*(1.+p[iphimax-iphi-1])+phimin;
	}
      //      cout << "phi[" << iphi << "]=" << phiArray[iphi] << ", " <<  phimax/2.*(1.-p[iphi]) << ", " 
      //   << (phimax)/2.*(1.+p[iphimax-iphi-1]) << endl;
    }

  // write information in the particle information file
  if (full==1) fprintf(d_file,"%d %e %e %e %e %e %d %d %d \n", number, deltaY, ymax, slope, phimin, phimax, iymax, iptmax, iphimax);

  // retrieve value as function of phi, pt and y in sumPtPhi:
  double Psi6pSin = 0.;
  double Psi6pCos = 0.;
  double Psi6p;
  double Psi5pSin = 0.;
  double Psi5pCos = 0.;
  double Psi5p;
  double Psi4pSin = 0.;
  double Psi4pCos = 0.;
  double Psi4p;
  double Psi3pSin = 0.;
  double Psi3pCos = 0.;
  double Psi3p;
  double Psi2pSin = 0.;
  double Psi2pCos = 0.;
  double Psi2p;
  double Psi1pSin = 0.;
  double Psi1pCos = 0.;
  double Psi1p;
  for (iy=0; iy<iymax; iy++)
    {
      y =  particleList[j].y[iy];
      //      fprintf(stderr,"y=%f \n", y);
      sumpt=0.;
      for (ipt=0; ipt<iptmax; ipt++)
	{
	  pt = particleList[j].pt[ipt];
	  //fprintf(stderr,"pt=%f \n", pt);
	  for (iphi=0; iphi<iphimax; iphi++)
	    {
	      phi = phiArray[iphi];
	      sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
	      
	      if (sum>1000) 
		cout << "WARNING: sum=" << sum << ">1000 for j=" << j << ", iy=" << iy << ", ipt=" << ipt << ", iphi=" << iphi << endl; 

	      if (full==1) fprintf(s_file,"%e ", sum);
	      // integrate over pt and phi
	      //     sumpt += deltaPhi * ( sum ) * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
	      // Psi3pSin += deltaPhi * ( sum*pt*sin(3.*phi) ) * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
	      //Psi3pCos += deltaPhi * ( sum*pt*cos(3.*phi) ) * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
	  if (iphi<iphimax/2)
		    {
		      sumpt += w[iphi] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi6pSin += w[iphi] * ( sum*sin(6.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi6pCos += w[iphi] * ( sum*cos(6.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi5pSin += w[iphi] * ( sum*sin(5.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi5pCos += w[iphi] * ( sum*cos(5.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi4pSin += w[iphi] * ( sum*sin(4.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi4pCos += w[iphi] * ( sum*cos(4.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi3pSin += w[iphi] * ( sum*sin(3.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi3pCos += w[iphi] * ( sum*cos(3.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi2pSin += w[iphi] * ( sum*sin(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi2pCos += w[iphi] * ( sum*cos(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi1pSin += w[iphi] * ( sum*sin(1.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi1pCos += w[iphi] * ( sum*cos(1.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		    }
		  else
		    {
		      sumpt += w[iphimax-iphi-1] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi6pSin += w[iphimax-iphi-1] * ( sum*sin(6.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi6pCos += w[iphimax-iphi-1] * ( sum*cos(6.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi5pSin += w[iphimax-iphi-1] * ( sum*sin(5.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi5pCos += w[iphimax-iphi-1] * ( sum*cos(5.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi4pSin += w[iphimax-iphi-1] * ( sum*sin(4.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi4pCos += w[iphimax-iphi-1] * ( sum*cos(4.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi3pSin += w[iphimax-iphi-1] * ( sum*sin(3.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi3pCos += w[iphimax-iphi-1] * ( sum*cos(3.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi2pSin += w[iphimax-iphi-1] * ( sum*sin(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi2pCos += w[iphimax-iphi-1] * ( sum*cos(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi1pSin += w[iphimax-iphi-1] * ( sum*sin(1.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi1pCos += w[iphimax-iphi-1] * ( sum*cos(1.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		    }
	    }
	  if (full==1) fprintf(s_file,"\n");
	}
      fprintf(stderr,"y=%f, sumpt=%f\n", y, sumpt/slope);
      fprintf(y_file,"%f %f\n", y, sumpt/slope);
    }

  if (Psi1pCos<0.)
    Psi1p = (atan(Psi1pSin/Psi1pCos)+PI);
  else
    Psi1p = (atan(Psi1pSin/Psi1pCos));

  if (Psi2pCos<0.)
    Psi2p = 1./2.*(atan(Psi2pSin/Psi2pCos)+PI);
  else
    Psi2p = 1./2.*(atan(Psi2pSin/Psi2pCos));
    
  if (Psi3pCos<0.) 
    Psi3p = 1./3.*(atan(Psi3pSin/Psi3pCos)+PI);
  else
    Psi3p = 1./3.*(atan(Psi3pSin/Psi3pCos));
  
  if (Psi4pCos<0.) 
    Psi4p = 1./4.*(atan(Psi4pSin/Psi4pCos)+PI);
  else
    Psi4p = 1./4.*(atan(Psi4pSin/Psi4pCos));
  
  if (Psi5pCos<0.) 
    Psi5p = 1./5.*(atan(Psi5pSin/Psi5pCos)+PI);
  else
    Psi5p = 1./5.*(atan(Psi5pSin/Psi5pCos));

  if (Psi6pCos<0.) 
    Psi6p = 1./6.*(atan(Psi6pSin/Psi6pCos)+PI);
  else
    Psi6p = 1./6.*(atan(Psi6pSin/Psi6pCos));
  
  cout << "Psi2p=" << Psi2p << endl; 
  cout << "Psi2pSin=" << Psi2pSin << endl; 
  cout << "Psi2pCos=" << Psi2pCos << endl; 
  cout << "Psi3p=" << Psi3p << endl; 
  cout << "Psi3pSin=" << Psi3pSin << endl; 
  cout << "Psi3pCos=" << Psi3pCos << endl; 

  ofstream fout1("Psin.dat",ios::out);
  fout1 << Psi1p << " " << Psi2p << " " << Psi3p << " " << Psi4p << " " << Psi5p << " " << Psi6p << endl;
  fout1.close();

  ofstream fout2("PsinCorrelatorsH+-.dat",ios::out);
  fout2 << cos(4.*(Psi2p-Psi4p)) << " "
	<< cos(8.*(Psi2p-Psi4p)) << " " 
	<< cos(12.*(Psi2p-Psi4p)) << " " 
	<< cos(6.*(Psi2p-Psi3p)) << " "
	<< cos(6.*(Psi2p-Psi6p)) << " " 
	<< cos(6.*(Psi3p-Psi6p)) << " "
	<< cos(12.*(Psi3p-Psi4p)) << " "
	<< cos(10.*(Psi2p-Psi5p)) << " "
	<< cos(2.*Psi2p+3.*Psi3p-5.*Psi5p) << " "
	<< cos(2.*Psi2p+4.*Psi4p-6.*Psi6p) << " "
	<< cos(2.*Psi2p-6.*Psi3p+4.*Psi4p) << " " 
	<< cos(-8.*Psi2p+3.*Psi3p+5.*Psi5p) << " "
	<< cos(-10.*Psi2p+4.*Psi4p+6.*Psi6p) << " "
	<< cos(-10.*Psi2p+6.*Psi3p+4.*Psi4p) << " "
	<< sin(3.*(Psi2p-Psi3p))*sin(5.*(Psi2p-Psi5p)) << " "
	<< cos(3.*(Psi2p-Psi3p))*cos(5.*(Psi2p-Psi5p)) << " "
	<< sin(4.*(Psi2p-Psi4p))*sin(6.*(Psi2p-Psi6p)) << " "
	<< cos(4.*(Psi2p-Psi4p))*cos(6.*(Psi2p-Psi6p)) << " "
	<< sin(6.*(Psi2p-Psi3p))*sin(4.*(Psi2p-Psi4p)) << " "
	<< cos(6.*(Psi2p-Psi3p))*cos(4.*(Psi2p-Psi4p)) 
	<<endl;
  fout2.close();

  cout << "doing pt spectrum at y=0" << endl;
  //pt spectra at y=0
  for (ipt=0; ipt<iptmax; ipt++)
    {
      sumpt = 0.;
      sumv[2] = 0.;
      int countY = 1;//= 0;
      pt = particleList[j].pt[ipt];
      for (iy=iymax/2; iy<=iymax/2; iy++) // modify
	{
	  y = iy*deltaY-ymax;
	  if (fabs(y)>ymaxIntegral) continue;
	  countY++;
	  for (iphi=0; iphi<iphimax; iphi++)
	    {
	      phi = phiArray[iphi];
	      sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
	      
	      if (iphi<iphimax/2)
		{
		  if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
		    {
		      sumpt += 0.5 * w[iphi] * ( sum ) * phiDiff  * 2.*PI/(phimax-phimin);// * deltaY;
		    }
		  else
		    {
		      sumpt += w[iphi] * ( sum ) * phiDiff  * 2.*PI/(phimax-phimin);// * deltaY
		      sumv[2] += w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * phiDiff * 2.*PI/(phimax-phimin);
		    }
		}
	      else
		{
		  if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
		    {
		      sumpt += 0.5 * w[iphimax-iphi-1] * ( sum ) * phiDiff * 2.*PI/(phimax-phimin);// * deltaY;
		    }
		  else
		    {
		      sumpt += w[iphimax-iphi-1] * ( sum ) * phiDiff * 2.*PI/(phimax-phimin);// * deltaY;
		      sumv[2] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * phiDiff * 2.*PI/(phimax-phimin);
		    }
		}
 	    }
	}
      fprintf(py0_file,"%e %e %f\n", pt, sumpt, y);
      fprintf(v2y0_file,"%e %e \n", pt, sumv[2]/sumpt);
    }

  cout << "doing pt spectrum integrated over y" << endl;
  //pt spectra integrated over y.
  for (ipt=0; ipt<iptmax; ipt++)
    {
      sumpt = 0.;
      int countY = 0;
      pt = particleList[j].pt[ipt];
      for (iy=0; iy<=iymax; iy++) // modify
	{
	  y = iy*deltaY-ymax;
	  if (fabs(y)>ymaxIntegral) continue;
	  countY++;
	  for (iphi=0; iphi<iphimax; iphi++)
	    {
	      phi = phiArray[iphi];
	      sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
	      
	      if (iphi<iphimax/2)
		{
		  if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
		    {
		      sumpt += 0.5 * w[iphi] * ( sum ) * phiDiff  * 2.*PI/(phimax-phimin)* deltaY;
		    }
		  else
		    {
		      sumpt += w[iphi] * ( sum ) * phiDiff  * 2.*PI/(phimax-phimin)* deltaY;
		    }
		  //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphi], iphi, p[iphi]);
		}
	      else
		{
		  if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
		    {
		      sumpt += 0.5 * w[iphimax-iphi-1] * ( sum ) * phiDiff * 2.*PI/(phimax-phimin)* deltaY;
		    }
		  else
		    {
		      sumpt += w[iphimax-iphi-1] * ( sum ) * phiDiff * 2.*PI/(phimax-phimin)* deltaY;
		    }
		  //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphimax-iphi-1], iphi, -p[iphimax-iphi-1]);
		}
 	    }
	}
      //fprintf(stderr,"pt=%f, sumpt=%f\n", pt, sumpt / (2.*ymaxIntegral));
      //      fprintf(p_file,"%e %e \n", pt, sumpt / (2.*ymaxIntegral));
      fprintf(p_file,"%e %e\n", pt, sumpt/(2.*ymaxIntegral));
      //fprintf(stderr,"pt=%f, sumv%d=%f\n", pt, 2, sumv[2]/sumpt);
    }


  cout << "doing v2,v3,v4" << endl;
  // v2, v3 and v4
  for (ipt=0; ipt<iptmax; ipt++)
    {
      sumpt = 0.;
      sumv[2] = 0.;
      sumv[3] = 0.;
      sumv3r3 = 0.;
      sumv[4] = 0.;
      int countY = 0;
      pt = particleList[j].pt[ipt];
      ymaxIntegral = 0.5*log((sqrt(pt*pt*cosh(etaMaxIntegral)*cosh(etaMaxIntegral)+m*m)
			      +pt*sinh(etaMaxIntegral))/(sqrt(pt*pt*cosh(etaMaxIntegral)*cosh(etaMaxIntegral)+m*m)-pt*sinh(etaMaxIntegral)));
      iymaxIntegral = floor(2.*ymaxIntegral/deltaY);
      for (iy=0; iy<iymax; iy++)
	{
	  y = iy*deltaY-ymax;
	  if (fabs(y)>ymaxIntegral) continue;
	  countY++;
	  for (iphi=0; iphi<iphimax; iphi++)
	    {
	      phi = phiArray[iphi];
	      sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
	      
	//       if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
// 		    {
// 		      sumpt += 0.5 * deltaPhi * ( sum ) * deltaY  * 2.*PI/(phimax-phimin);
// 		      sumv[2] += 0.5 * deltaPhi * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * 2.*PI/(phimax-phimin);
// 		      sumv[3] += 0.5 * deltaPhi * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * 2.*PI/(phimax-phimin);
// 		      sumv3r3 += 0.5 * deltaPhi * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * 2.*PI/(phimax-phimin);
// 		      sumv[4] += 0.5 * deltaPhi * ( sum*cos(4.*(phi-phimin)) ) * deltaY * 2.*PI/(phimax-phimin);
// 		    }
// 		  else
// 		    {
// 		      sumpt += deltaPhi * ( sum ) * deltaY * 2.*PI/(phimax-phimin);
// 		      sumv[2] += deltaPhi * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * 2.*PI/(phimax-phimin);
// 		      sumv[3] += deltaPhi * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * 2.*PI/(phimax-phimin);
// 		      sumv3r3 += deltaPhi * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * 2.*PI/(phimax-phimin);
// 		      sumv[4] += deltaPhi * ( sum*cos(4.*(phi-phimin)) ) * deltaY * 2.*PI/(phimax-phimin);
// 		    }
		
	      if (iphi<iphimax/2)
		{
		  if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
		    {
		      sumpt += 0.5 * w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
		      sumv[2] += 0.5 * w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[3] += 0.5 * w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      //sumv3r3 += 0.5 * w[iphi] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[4] += 0.5 * w[iphi] * ( sum*cos(4.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		    }
		  else
		    {
		      sumpt += w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
		      sumv[2] += w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[3] += w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      //sumv3r3 += w[iphi] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[4] += w[iphi] * ( sum*cos(4.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		    }
		  //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphi], iphi, p[iphi]);
		}
	      else
		{
		  if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
		    {
		      sumpt += 0.5 * w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[2] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[3] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      //sumv3r3 += 0.5 * w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[4] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		    }
		  else
		    {
		      sumpt += w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[2] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[3] += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      //sumv3r3 += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[4] += w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		    }
		  //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphimax-iphi-1], iphi, -p[iphimax-iphi-1]);
		}
 	    }
	}
      //fprintf(stderr,"pt=%f, sumpt=%f\n", pt, sumpt / (2.*ymaxIntegral));
      //fprintf(stderr,"pt=%f, sumv%d=%f\n", pt, 2, sumv[2]/sumpt);
      fprintf(v2_file,"%e %e \n", pt, sumv[2]/sumpt);
      fprintf(v3_file,"%e %e \n", pt, sumv[3]/sumpt);
      //fprintf(v3r3_file,"%e %e \n", pt, sumv3r3/sumpt);
      fprintf(v4_file,"%e %e \n", pt, sumv[4]/sumpt);
    }

  cout << "check" << endl;

  //pseudorapidity2
  deltaEta=deltaY;
  etamax=ymax;
  ietamax=floor(2.*etamax/deltaEta);
  for (ieta=0; ieta<ietamax; ieta++)
    {
      eta = ieta*deltaEta-etamax+deltaEta/2.;
      sumpt=0.;
      sumv[2]=0.;
      sumv[3]=0.;
      sumv3r3=0.;
      sumv[4]=0.;
      for (ipt=0; ipt<iptmax; ipt++)
	{
	  sumpteta=0.;
	  sumvpteta[2]=0.;
	  sumvpteta[3]=0.;
	  sumvpteta[4]=0.;
	  pt = particleList[j].pt[ipt];
	  y = 0.5*log((sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)+pt*sinh(eta))/(sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)-pt*sinh(eta)));
	 
	  for(i=0; i<iymax-1; i++) // find closest iy
	    {
	      if ( particleList[j].y[i]< y &&  particleList[j].y[i+1] > y )
		{
		  if ( fabs(particleList[j].y[i]-y) <  fabs(particleList[j].y[i+1]-y)) 
		    iy = i;
		  else
		    iy = i+1;
		}
	    }
	  mt = sqrt(pt*pt+m*m);
	  for (iphi=0; iphi<iphimax; iphi++)
	    {
	      phi = phiArray[iphi];
	      px = pt*cos(phi);
	      py = pt*sin(phi);
	      sum = particleList[j].dNdydptdphi[iy][ipt][iphi]*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
	      if (iphi<iphimax/2)
		{
	       	  sumpt += w[iphi] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		  sumv[2] += w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		  sumv[3] += w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		  //sumv3r3 += w[iphi] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		  sumv[4] += w[iphi] * ( sum*cos(4.*(phi-phimin)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		  sumpteta += w[iphi] * ( sum ) * phiDiff * 2.*PI/(phimax-phimin);
		  sumvpteta[2] += w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * phiDiff * 2.*PI/(phimax-phimin);
		  sumvpteta[3] += w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * phiDiff * 2.*PI/(phimax-phimin);
		  sumvpteta[4] += w[iphi] * ( sum*cos(4.*(phi-phimin)) ) * phiDiff * 2.*PI/(phimax-phimin);		
		}
	      else
		{
		  sumpt += w[iphimax-iphi-1] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		  sumv[2] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		  sumv[3] += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		  //sumv3r3 += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		  sumv[4] += w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		  sumpteta += w[iphimax-iphi-1] * ( sum ) * phiDiff * 2.*PI/(phimax-phimin);
		  sumvpteta[2] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * phiDiff * 2.*PI/(phimax-phimin);
		  sumvpteta[3] += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * phiDiff * 2.*PI/(phimax-phimin);
		  sumvpteta[4] += w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin)) ) * phiDiff * 2.*PI/(phimax-phimin);
		}
	    }
	  fprintf(pteta_file,"%f %f %e\n", eta, pt, sumpteta);
	  fprintf(v2_pteta_file,"%f %f %e\n", eta, pt, sumvpteta[2]/sumpteta);
	  fprintf(v3_pteta_file,"%f %f %e\n", eta, pt, sumvpteta[3]/sumpteta);
	  fprintf(v4_pteta_file,"%f %f %e\n", eta, pt, sumvpteta[4]/sumpteta);
	}
      //   fprintf(stderr,"eta=%f, sumpt=%f\n", eta, sumpt/slope);
      fprintf(eta_file,"%f %f\n", eta, sumpt/slope);
      fprintf(v2_eta_file,"%f %f\n", eta, sumv[2]/sumpt);
      fprintf(v3_eta_file,"%f %f\n", eta, sumv[3]/sumpt);
      // fprintf(v3r3_eta_file,"%f %f\n", eta, sumv3r3/sumpt);
      fprintf(v4_eta_file,"%f %f\n", eta, sumv[4]/sumpt);
      fprintf(pteta_file,"\n");
      fprintf(v2_pteta_file,"\n");
      fprintf(v3_pteta_file,"\n");
      fprintf(v4_pteta_file,"\n");
    
    }


  //total v2, v3 and v4
  deltaEta=deltaY;
  etamax=ymax;
  ietamax=floor(2.*etamax/deltaEta);
  double sumtot =0.;
  double sumtotv2 =0.;
  double sumtotv3 =0.;
  double sumtotv3r3 =0.;
  double sumtotv4 =0.;
  int eta1=0;
  for (ieta=0; ieta<ietamax; ieta++)
    {
      eta = ieta*deltaEta-etamax+deltaEta/2.;
      if (fabs(eta)>1.5) continue; // experimental cut
      eta1+=1;
      sumpt=0.;
      sumv[2]=0.;
      sumv[3]=0.;
      sumv3r3=0.;
      sumv[4]=0.;
      for (ipt=0; ipt<iptmax; ipt++)
	{
	  pt = particleList[j].pt[ipt];
	  if (pt < 0.8 || pt > 4.) continue; // experimental cut 
	  y = 0.5*log((sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)+pt*sinh(eta))/(sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)-pt*sinh(eta)));
	 
	  for(i=0; i<iymax-1; i++) // find closest iy
	    {
	      if ( particleList[j].y[i]< y &&  particleList[j].y[i+1] > y )
		{
		  if ( fabs(particleList[j].y[i]-y) <  fabs(particleList[j].y[i+1]-y)) 
		    iy = i;
		  else
		    iy = i+1;
		}
	    }

	  mt = sqrt(pt*pt+m*m);
	  for (iphi=0; iphi<iphimax; iphi++)
	    {
	      phi = phiArray[iphi];
	      px = pt*cos(phi);
	      py = pt*sin(phi);
	      sum = particleList[j].dNdydptdphi[iy][ipt][iphi]*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
	      if (iphi<iphimax/2)
		{
		  sumpt += w[iphi] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		  sumv[2] += w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		  sumv[3] += w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		  //sumv3r3 += w[iphi] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		  sumv[4] += w[iphi] * ( sum*cos(4.*(phi-phimin)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		}
	      else
		{
		  sumpt += w[iphimax-iphi-1] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		  sumv[2] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		  sumv[3] += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		  //sumv3r3 += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		  sumv[4] += w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		}
	    }
	}

      sumtot += sumpt * deltaEta;
      if (eta1 == 1)
	sumtot -= 0.5 * sumpt * deltaEta;
      sumtotv2 += sumv[2] * deltaEta;
      if (eta1 == 1)
	sumtotv2 -= 0.5 * sumv[2] * deltaEta;
      sumtotv3 += sumv[3] * deltaEta;
      if (eta1 == 1)
	sumtotv3 -= 0.5 * sumv[3] * deltaEta;
      // sumtotv3r3 += sumv3r3 * deltaEta;
      //if (eta1 == 1)
      //	sumtotv3r3 -= 0.5 * sumv3r3 * deltaEta;
      sumtotv4 += sumv[4] * deltaEta;
      if (eta1 == 1)
	sumtotv4 -= 0.5 * sumv[4] * deltaEta;
    }
  sumtot -= 0.5 * sumpt * deltaEta;
  sumtotv2 -= 0.5 * sumv[2] * deltaEta;
  sumtotv3 -= 0.5 * sumv[3] * deltaEta;
  //sumtotv3r3 -= 0.5 * sumv3r3 * deltaEta;
  sumtotv4 -= 0.5 * sumv[4] * deltaEta;
  
  sumtotv2/=sumtot;
  sumtotv3/=sumtot;
  //sumtotv3r3/=sumtot;
  sumtotv4/=sumtot;
  
  fprintf(v2_tot_file,"%f\n", sumtotv2);
  fprintf(v3_tot_file,"%f\n", sumtotv3);
  //fprintf(v3r3_tot_file,"%f\n", sumtotv3r3);
  fprintf(v4_tot_file,"%f\n", sumtotv4);
  
  
  //as function of phi for 'correlations'
  deltaEta=deltaY;
  etamax=ymax;
  ietamax=floor(2.*etamax/deltaEta);
  for (iphi=0; iphi<iphimax; iphi++)
    {
      sumPhi=0.;
      phi = phiArray[iphi];
      for (ipt=iptmax-3; ipt<iptmax-2; ipt++)
	{
	  pt = particleList[j].pt[ipt];
	  y = 0.;
	  cout << "pt=" << pt << endl;
	  iy = floor((y+ymax)/deltaY+0.0001);

	  cout << "y=" << iy*deltaY-ymax << endl;
	  
	  px = pt*cos(phi);
	  py = pt*sin(phi);
	  mt = sqrt(pt*pt+m*m);
	  sum = particleList[j].dNdydptdphi[iy][ipt][iphi]*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
	  
	  sumPhi = sum * 0.7; // 0.7 for eta range -0.35 < eta < 0.35
	}
      
      fprintf(phi_file,"%e %e\n", phi, sumPhi);
    }
  

  //pt integrated over eta
  int countY = 0;
  deltaEta=deltaY;
  etamax=1.; //ATLAS
  ietamax=floor(2.*etamax/deltaEta);
  
  for (ipt=0; ipt<iptmax; ipt++)
    {
      sumpt=0.;
      sumv[1]=0.;
      sumv[2]=0.;
      sumv[3]=0.;
      sumv[4]=0.;
      sumv[5]=0.;
      sumv[6]=0.;
      for (ieta=0; ieta<ietamax; ieta++)
	{
	  eta = ieta*deltaEta-etamax+deltaEta/2.;
	  // j = partid[MHALF+setOfNumbers[ip]];
	  m = particleList[j].mass;
	  pt = particleList[j].pt[ipt];
	  y = 0.5*log((sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)+pt*sinh(eta))/(sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)-pt*sinh(eta)));
	  
	  for(i=0; i<iymax-1; i++) // find closest iy
	    {
	      if ( particleList[j].y[i]< y && particleList[j].y[i+1] > y )
		{
		  // if ( fabs(particleList[j].y[i]-y) < fabs(particleList[j].y[i+1]-y)) 
		  iy = i;
		  //	      else
		  //	iy = i+1;
		}
	    }
	  
	  mt = sqrt(pt*pt+m*m);
	  for (iphi=0; iphi<iphimax; iphi++)
	    {
	      phi = phiArray[iphi];
	      px = pt*cos(phi);
	      py = pt*sin(phi);
	      sum = (
		     (particleList[j].dNdydptdphi[iy+1][ipt][iphi])*(particleList[j].y[iy]-y)/(particleList[j].y[iy+1]-particleList[j].y[iy])
		     +(particleList[j].dNdydptdphi[iy][ipt][iphi])*(1.-(particleList[j].y[iy]-y)/(particleList[j].y[iy+1]-particleList[j].y[iy]))
		     )
		*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
	      //interpolation in eta
	      
	      if (iphi<iphimax/2)
		{
		  if (ieta==0 || ieta==ietamax) // trapezoid rule: only half the edges
		    {
		      sumpt += 0.5 * w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
		      sumv[1] += 0.5 * w[iphi] * ( sum*cos(1.*(phi-phimin-Psi1p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[2] += 0.5 * w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[3] += 0.5 * w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[4] += 0.5 * w[iphi] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[5] += 0.5 * w[iphi] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[6] += 0.5 * w[iphi] * ( sum*cos(6.*(phi-phimin-Psi6p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		    }
		  else
		    {
		      sumpt += w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
		      sumv[1] += w[iphi] * ( sum*cos(1.*(phi-phimin-Psi1p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[2] += w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[3] += w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[4] += w[iphi] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[5] += w[iphi] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[6] += w[iphi] * ( sum*cos(6.*(phi-phimin-Psi6p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		    }
		}
	      else
		{
		  if (ieta==0 || ieta==ietamax) // trapezoid rule: only half the edges
		    {
		      sumpt += 0.5 * w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[1] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(1.*(phi-phimin-Psi1p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[2] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[3] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[4] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[5] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[6] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(6.*(phi-phimin-Psi6p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		    }
		  else
		    {
		      sumpt += w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[1] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi1p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[2] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[3] += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[4] += w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[5] += w[iphimax-iphi-1] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[6] += w[iphimax-iphi-1] * ( sum*cos(6.*(phi-phimin-Psi6p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
		    }
		}
	    }
	}
      fprintf(pt_file,"%e %e\n", pt, sumpt/(2*etamax));
      fprintf(v1_pt_file,"%e %e\n", pt, sumv[1]/sumpt);
      fprintf(v2_pt_file,"%e %e\n", pt, sumv[2]/sumpt);
      fprintf(v3_pt_file,"%e %e\n", pt, sumv[3]/sumpt);
      fprintf(v4_pt_file,"%e %e\n", pt, sumv[4]/sumpt);
      fprintf(v5_pt_file,"%e %e\n", pt, sumv[5]/sumpt);
      fprintf(v6_pt_file,"%e %e\n", pt, sumv[6]/sumpt);
    }

  
  // simple integration for check
/*   for (ipt=0; ipt<2.*iptmax; ipt++) */
/*     {    */
/*       deltaPT = 0.05; */
/*       sumpt = 0.; */
/*       sumv[2] = 0.; */
/*       for (iy=0; iy<iymax; iy++) */
/* 	{ */
/* 	  y = iy*deltaY-ymax; */
/* 	  pt = ipt*deltaPT; // legendre abscissa / slope */
/* 	  deltaPhi=0.01; */
/* 	  iphimax = floor(phimax/deltaPhi); */
/* 	  for (iphi=0; iphi<iphimax; iphi++) */
/* 	    { */
/* 	      phi = iphi*deltaPhi; */
/* 	      px = pt*cos(phi); */
/* 	      py = pt*sin(phi); */
/* 	      sum = summation(px, py, y, m, d, b); */
/* 	      sumpt += ( sum ) * deltaY * deltaPhi * 2.*PI/phimax; */
/* 	      sumv[2] += ( sum*cos(2.*phi) ) * deltaY * deltaPhi * 2.*PI/phimax; */
/* 	    } */
/* 	} */
/*       fprintf(stderr,"pt=%f, sumpt=%f\n", pt, sumpt / (2.*ymax)); */
/*       fprintf(stderr,"pt=%f, sumv%d=%f\n", pt, 2, sumv[2]/sumpt); */
/*     } */

  if (full==1) fclose(s_file);
  if (full==1) fclose(d_file);
  fclose(p_file);
  fclose(phi_file);
  fclose(v2_eta_file);
  fclose(v4_eta_file);
  fclose(v2_tot_file);
  fclose(v3_tot_file);
  fclose(v4_tot_file);
  fclose(eta_file);
  fclose(y_file);
  fclose(v2_file);
  fclose(v3_file);
  fclose(v3r3_file);
  fclose(v4_file);
  fprintf(stderr,"Done with %s (%d)\n", particleList[j].name, particleList[j].number);
  fclose(pteta_file);
  fclose(v3_eta_file);
  fclose(v3r3_eta_file);
  fclose(v2_pteta_file);
  fclose(v3_pteta_file);
  fclose(v4_pteta_file);
  fclose(v3r3_tot_file);
  fclose(py0_file);
  fclose(v2y0_file);
  fclose(pt_file);
  fclose(v1_pt_file);
  fclose(v2_pt_file);
  fclose(v3_pt_file);
  fclose(v4_pt_file);
  fclose(v5_pt_file);
  fclose(v6_pt_file);
}


// read in full spectra from file to then perform anything with them
void Freeze::ReadFullSpectra(InitData* DATA)
{
  // read in thermal spectra from file:
  int number, iymax, iptmax, iphimax;
  double deltaY, ymax, slope, phimax, phimin;
  int ip, iphi, ipt, i;
  double *p, *w;		        // pointing to data for Gaussian integration in phi 
  int bytes_read;
  static char *s;
  s = util->char_malloc(120);
  int iy, j, k, d1, d2, d3, decays, h;
  double b, npi, nK, neta, dummy;
  fprintf(stderr,"reading spectra\n");
  char *anti;
  // open particle information file:
  FILE *p_file;
  const char* p_name = "FparticleInformation.dat";
  p_file = fopen(p_name, "r");
  checkForReadError(p_file,p_name);
  int count;
  count = 0;
  // read particle information:
  while( fscanf(p_file,"%d %lf %lf %lf %lf %lf %d %d %d ",&number, &deltaY, &ymax, &slope, &phimin, &phimax, &iymax, &iptmax, &iphimax) == 9)
    //while( fscanf(p_file,"%d %lf %lf %lf %lf %d %d %d ",&number, &deltaY, &ymax, &slope, &phimax, &iymax, &iptmax, &iphimax) == 8)
    {
      count ++;
     
      //phimin=0.;
     
      if (count>DATA->NumberOfParticlesToInclude) break;
      fprintf(stderr,"%d %e %e %e %e %e %d %d %d \n", number, deltaY, ymax, slope, phimin, phimax, iymax, iptmax, iphimax);
      ip = partid[MHALF+number];
      particleList[ip].ny = iymax;
      particleList[ip].npt = iptmax;
      particleList[ip].nphi = iphimax;
      particleList[ip].phimin = phimin;
      particleList[ip].phimax = phimax;
      particleList[ip].slope = slope;
      particleList[ip].ymax = ymax;
      particleList[ip].deltaY = deltaY;
      for ( i=0; i<iptmax; i++ )
	{
	  particleList[ip].pt[i] =  gala15x[i]/slope;
	}
      for ( i=0; i<iymax; i++ )
	{
	  particleList[ip].y[i] =  i*deltaY-ymax+deltaY/2.;
	}

      switch (iphimax) 
	{
	case 4: p= gaulep4; w= gaulew4; break;
	case 8: p= gaulep8; w= gaulew8; break;
	case 10: p= gaulep10; w= gaulew10; break;
	case 12: p= gaulep12; w= gaulew12; break;
	case 16: p= gaulep16; w= gaulew16; break;
	case 20: p= gaulep20; w= gaulew20; break;
	case 48: p= gaulep48; w= gaulew48; break;
	default: fprintf(stderr,"specified number of phi-points %d not available\n",iphimax); exit(1);
	}
      
      phiArray = util->vector_malloc(iphimax);

      for(iphi=0; iphi<iphimax; iphi++)
	{
	  if ( iphi < iphimax/2 )
	    {
	      phiArray[iphi] = (phimax-phimin)/2.*(1.-p[iphi])+phimin;
	    }
	  else
	    {
	      phiArray[iphi] = (phimax-phimin)/2.*(1.+p[iphimax-iphi-1])+phimin;
	    }
	  //      cout << "phi[" << iphi << "]=" << phiArray[iphi] << ", " <<  phimax/2.*(1.-p[iphi]) << ", " 
	  //   << (phimax)/2.*(1.+p[iphimax-iphi-1]) << endl;
	}
    }

  particleMax = ip;

  fclose(p_file);

  FILE *s_file;
  const char* s_name = "FyptphiSpectra.dat";
  s_file = fopen(s_name, "r");
  checkForReadError(s_file,s_name);
  
  for ( ip=1; ip<=particleMax; ip++ )
    {
      //fprintf(stderr,"reading particle %d: %d %s\n", ip, particleList[ip].number, particleList[ip].name);
      for (iy=0; iy<iymax; iy++)
	{
	  for (ipt=0; ipt<iptmax; ipt++)
	    {
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  bytes_read=fscanf(s_file, "%lf", &particleList[ip].dNdydptdphi[iy][ipt][iphi]);
		  //printf("%f %f %f \n",particleList[ip].y[iy],particleList[ip].pt[ipt],phiArray[iphi]);
		}
	    }
	}
    }
  fclose(s_file);
}

// read in full spectra from file to then perform anything with them
void Freeze::ReadFullSpectra2(InitData* DATA)
{
  // read in thermal spectra from file:
  int number, iymax, iptmax, iphimax;
  double deltaY, ymax, slope, phimax, phimin;
  int ip, iphi, ipt, i;
  double *p, *w;		        // pointing to data for Gaussian integration in phi 
  int bytes_read;
  static char *s;
  s = util->char_malloc(120);
  int iy, j, k, d1, d2, d3, decays, h;
  double b, npi, nK, neta, dummy;
  fprintf(stderr,"reading spectra\n");
  char *anti;
  // open particle information file:
  FILE *p_file;
  const char* p_name = "FparticleInformation.dat";
  p_file = fopen(p_name, "r");
  checkForReadError(p_file,p_name);
  int count;
  count = 0;
  // read particle information:
  while( fscanf(p_file,"%d %lf %lf %lf %lf %lf %d %d %d ",&number, &deltaY, &ymax, &slope, &phimin, &phimax, &iymax, &iptmax, &iphimax) == 9)
    //while( fscanf(p_file,"%d %lf %lf %lf %lf %d %d %d ",&number, &deltaY, &ymax, &slope, &phimax, &iymax, &iptmax, &iphimax) == 8)
    {
      count ++;
     
      //phimin=0.;
     
      if (count>DATA->NumberOfParticlesToInclude) break;
      fprintf(stderr,"%d %e %e %e %e %e %d %d %d \n", number, deltaY, ymax, slope, phimin, phimax, iymax, iptmax, iphimax);
      ip = partid[MHALF+number];
      particleList[ip].ny = iymax;
      particleList[ip].npt = iptmax;
      particleList[ip].nphi = iphimax;
      particleList[ip].phimin = phimin;
      particleList[ip].phimax = phimax;
      particleList[ip].slope = slope;
      particleList[ip].ymax = ymax;
      particleList[ip].deltaY = deltaY;
      for ( i=0; i<iptmax; i++ )
	{
	  particleList[ip].pt[i] =  gala15x[i]/slope;
	}
      for ( i=0; i<iymax; i++ )
	{
	  particleList[ip].y[i] =  i*deltaY-ymax+deltaY/2.;
	}

      switch (iphimax) 
	{
	case 4: p= gaulep4; w= gaulew4; break;
	case 8: p= gaulep8; w= gaulew8; break;
	case 10: p= gaulep10; w= gaulew10; break;
	case 12: p= gaulep12; w= gaulew12; break;
	case 16: p= gaulep16; w= gaulew16; break;
	case 20: p= gaulep20; w= gaulew20; break;
	case 48: p= gaulep48; w= gaulew48; break;
	default: fprintf(stderr,"specified number of phi-points %d not available\n",iphimax); exit(1);
	}
      
      phiArray = util->vector_malloc(iphimax);

      for(iphi=0; iphi<iphimax; iphi++)
	{
	  if ( iphi < iphimax/2 )
	    {
	      phiArray[iphi] = (phimax-phimin)/2.*(1.-p[iphi])+phimin;
	    }
	  else
	    {
	      phiArray[iphi] = (phimax-phimin)/2.*(1.+p[iphimax-iphi-1])+phimin;
	    }
	  //      cout << "phi[" << iphi << "]=" << phiArray[iphi] << ", " <<  phimax/2.*(1.-p[iphi]) << ", " 
	  //   << (phimax)/2.*(1.+p[iphimax-iphi-1]) << endl;
	}
    }

  particleMax = ip;

  fclose(p_file);

  FILE *s_file;
  const char* s_name = "FyptphiSpectra2.dat";
  s_file = fopen(s_name, "r");
  checkForReadError(s_file,s_name);
  
  for ( ip=1; ip<=particleMax; ip++ )
    {
      //fprintf(stderr,"reading particle %d: %d %s\n", ip, particleList[ip].number, particleList[ip].name);
      for (iy=0; iy<iymax; iy++)
	{
	  for (ipt=0; ipt<iptmax; ipt++)
	    {
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  bytes_read=fscanf(s_file, "%lf", &particleList[ip].dNdydptdphi2[iy][ipt][iphi]);
		  //printf("%f %f %f \n",particleList[ip].y[iy],particleList[ip].pt[ipt],phiArray[iphi]);
		}
	    }
	}
    }
  fclose(s_file);
}


void Freeze::ComputeAveragePT(int number, double ptmax)
{
  char *numberString;
  char buf[10];
  double *p, *w;		        // pointing to data for Gaussian integration in phi 
  double slope;
  int ipt, iphi, iymax, iy, iptmax, iphimax, ieta, ietamax;
  int i,j, iymaxIntegral;
  double pt, phi, px, py, y, deltaPT, deltaPhi, deltaY, ymax, phimin, phimax, sum, sumpt, sumnorm, sumv[5];
  double phiOffs, phiDiff, ymaxIntegral, eta, deltaEta, etamax, mt, etaMaxIntegral;
  int returnValue;
  int ip;

  // set some parameters
  
  j = partid[MHALF+number];
  fprintf(stderr,"Doing %s (%d)\n", particleList[j].name, particleList[j].number);
  phimax = particleList[j].phimax;
  phimin = particleList[j].phimin;
  phiOffs = 0.5 * ( phimin + phimax );
  phiDiff = 0.5 * ( phimax - phimin );
  iphimax = particleList[j].nphi;
  iptmax = particleList[j].npt;

  if (iptmax != 15) 
    {
      fprintf(stderr,"only iptmax=15 possible. you picked %d. exiting.\n", iptmax);
      exit(1);
    }

  ymax = particleList[j].ymax;
  deltaY = particleList[j].deltaY;
  ymaxIntegral = 0.5-deltaY/2.;
  etaMaxIntegral = 1.3; // for v2 integration
  iymax = particleList[j].ny;
  iymaxIntegral = floor(2.*ymaxIntegral/deltaY);
  iptmax = particleList[j].npt;
  iphimax = particleList[j].nphi;
  
  // set particle properties
  double m = particleList[j].mass;
  int d = particleList[j].degeneracy;
  int b = particleList[j].baryon;
  int s = particleList[j].strange;
  double c = particleList[j].charge;
  double mu = particleList[j].muAtFreezeOut;
 
  // set up file names
  numberString = util->char_malloc(30);
  
  if( chdir("./outputs") != 0 )
    {
      fprintf(stderr,"directory \"outputs\" does not exist. Exiting.\n");
      exit(1);
    }
  else returnValue=chdir("..");

  sprintf (buf, "%d", number);
  strcat(numberString, "./outputs/FavPT-");
  strcat(numberString, buf);
  strcat(numberString,".dat");
    
  FILE *p_file;
  char* p_name = numberString;
  p_file = fopen(p_name, "w");

  slope = particleList[j].slope;
  switch (iphimax) 
    {
    case 4: p= gaulep4; w= gaulew4; break;
    case 8: p= gaulep8; w= gaulew8; break;
    case 10: p= gaulep10; w= gaulew10; break;
    case 12: p= gaulep12; w= gaulew12; break;
    case 16: p= gaulep16; w= gaulew16; break;
    case 20: p= gaulep20; w= gaulew20; break;
    case 48: p= gaulep48; w= gaulew48; break;
    default: fprintf(stderr,"specified number of phi-points not available\n"); exit(1);
    }

  for(iphi=0; iphi<iphimax; iphi++)
    {
      if ( iphi < iphimax/2 )
	{
	  phiArray[iphi] = (phimax-phimin)/2.*(1.-p[iphi])+phimin;
	}
      else
	{
	  phiArray[iphi] = (phimax-phimin)/2.*(1.+p[iphimax-iphi-1])+phimin;
	}
    }

  cout << "starting..." << endl;

  //compute average p_T as function of y:
  for (iy=0; iy<iymax; iy++)
    {
      y =  particleList[j].y[iy];
      //      fprintf(stderr,"y=%f \n", y);
      sumpt=0.;
      sumnorm=0.;
      for (ipt=0; ipt<iptmax; ipt++)
	{
	  pt = particleList[j].pt[ipt];
	  //fprintf(stderr,"pt=%f \n", pt);
	  for (iphi=0; iphi<iphimax; iphi++)
	    {
	      phi = phiArray[iphi];
	      sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
	      // integrate over pt and phi
	      if (iphi<iphimax/2)
		{
		  sumpt += w[iphi] * ( sum ) * phiDiff * gala15w[ipt] * pt * pt * 2.*PI/(phimax-phimin);
		  sumnorm += w[iphi] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		  //fprintf(stderr,"sumpt=%f \n", sumpt);
		}
	      else
		{
		  sumpt += w[iphimax-iphi-1] * ( sum ) * phiDiff * gala15w[ipt] * pt * pt * 2.*PI/(phimax-phimin);
		  sumnorm += w[iphimax-iphi-1] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		  //fprintf(stderr,"sumpt=%f \n", sumpt);
		}
	    }
	}
      fprintf(stderr,"y=%f, sumpt=%f\n", y, sumpt/sumnorm);
      fprintf(p_file,"%f %f\n", y, sumpt/sumnorm);
    }
  fclose(p_file);
  util->char_free(numberString);
}

void Freeze::ComputeChargedHadrons(InitData* DATA, double ptmax)
{
  char *numberStringy;
  char *numberString;
  char *numberStringPhi;
  char *numberStringv2;
  char *numberStringv2reac;
  char *numberStringv2eta;
  char *numberStringv2tot;
  char *numberStringv3;
  char *numberStringv3reac;
  char *numberStringv3r3;
  char *numberStringv3eta;
  char *numberStringv3tot;
  char *numberStringv3r3eta;
  char *numberStringv3r3tot;
  char *numberStringv4;
  char *numberStringv4reac;
  char *numberStringv4eta;
  char *numberStringv4tot;
  char *numberStringv5;
  char *numberStringv5reac;
  char *numberStringv5eta;
  char *numberStringv5tot;
  char *numberStringv6tot;
  char *numberStringv6eta;
  char *numberStringv1tot;
  char *numberStringv1eta;
  char *numberStringeta;
  char *numberStringpteta;
  char *numberStringv2pteta;
  char *numberStringv3pteta;
  char *numberStringv4pteta;
  char *numberStringv5pteta;

  char *numberStringpt;
  char *numberStringv1pt;
  char *numberStringv2pt;
  char *numberStringv3pt;
  char *numberStringv4pt;
  char *numberStringv5pt;
  char *numberStringv6pt;

  char buf[10];
  double *p, *w;		        // pointing to data for Gaussian integration in phi 
  double slope, slope1, slope2, fleft1, fleft2, fright1, fright2;
  int ipt, iphi, iymax, iy, iptmax, iphimax, ieta, ietamax;
  int i,j, iymaxIntegral;
  double pt, phi, px, py, y, deltaPT, deltaPhi, deltaY, ymax, phimin, phimax, sum, sumpt, sumpt2, sumv[7], m, sumv3r3, sumpteta, sumvpteta[7];
  double phiOffs, phiDiff, ymaxIntegral, eta, deltaEta, etamax, mt, etaMaxIntegral;
  int returnValue;
  int number;
  int setOfNumbers[6] = {211,-211,2212,-2212,321,-321};
  int ip;
  
  // set some parameters

  number = setOfNumbers[0]; //use pion to get settings

  j = 1;
  fprintf(stderr,"Doing charged hadrons\n");
  phimax = particleList[j].phimax;
  phimin = particleList[j].phimin;
  
  phiOffs = 0.5 * ( phimin + phimax );
  phiDiff = 0.5 * ( phimax - phimin );
  iphimax = particleList[j].nphi;
  iptmax = particleList[j].npt;

  if (iptmax != 15) 
    {
      fprintf(stderr,"only iptmax=15 possible. you picked %d. exiting.\n", iptmax);
      exit(1);
    }

  ymax = particleList[j].ymax;
  deltaY = particleList[j].deltaY;
  ymaxIntegral = 0.5-deltaY/2.;
  etaMaxIntegral = 0.35; // for v2 integration
  iymax = particleList[j].ny;
  iymaxIntegral = floor(2.*ymaxIntegral/deltaY);
  iptmax = particleList[j].npt;
  iphimax = particleList[j].nphi;

  cout << "setting up file names..."  << endl;
  // set up file names
  numberStringy = util->char_malloc(30);
  numberString = util->char_malloc(30);
  numberStringPhi = util->char_malloc(30);
  numberStringeta = util->char_malloc(30);
  numberStringv2 = util->char_malloc(30);
  numberStringv3 = util->char_malloc(30);
  numberStringv2reac = util->char_malloc(30);
  numberStringv3reac = util->char_malloc(30);
  numberStringv4reac = util->char_malloc(30);
  numberStringv5reac = util->char_malloc(30);
  numberStringv3r3 = util->char_malloc(30);
  numberStringv4 = util->char_malloc(30);
  numberStringv5 = util->char_malloc(30);
  numberStringv2eta = util->char_malloc(30);
  numberStringv3eta = util->char_malloc(30);
  numberStringv3r3eta = util->char_malloc(30);
  numberStringv4eta = util->char_malloc(30);
  numberStringv5eta = util->char_malloc(30);
  numberStringv6eta = util->char_malloc(30);
  numberStringv2tot = util->char_malloc(30);
  numberStringv3tot = util->char_malloc(30);
  numberStringv3r3tot = util->char_malloc(30);
  numberStringv4tot = util->char_malloc(30);
  numberStringv5tot = util->char_malloc(30);
  numberStringv6tot = util->char_malloc(30);
  numberStringv1tot = util->char_malloc(30);
  numberStringv1eta = util->char_malloc(30);
  numberStringpteta = util->char_malloc(30);
  numberStringv2pteta = util->char_malloc(30);
  numberStringv3pteta = util->char_malloc(30);
  numberStringv4pteta = util->char_malloc(30);
  numberStringv5pteta = util->char_malloc(30);

  numberStringpt = util->char_malloc(30);
  numberStringv1pt = util->char_malloc(30);
  numberStringv2pt = util->char_malloc(30);
  numberStringv3pt = util->char_malloc(30);
  numberStringv4pt = util->char_malloc(30);
  numberStringv5pt = util->char_malloc(30);
  numberStringv6pt = util->char_malloc(30);
  
  if( chdir("./outputs") != 0 )
    {
      fprintf(stderr,"directory \"outputs\" does not exist. Exiting.\n");
      exit(1);
    }
  else returnValue=chdir("..");
  
  
  strcat(numberString, "./outputs/FNpt-H+-");
  strcat(numberString,".dat");
  
  strcat(numberStringpteta, "./outputs/FNetapt-H+-");
  strcat(numberStringpteta,".dat");

  strcat(numberStringv2pteta, "./outputs/Fv2etapt-H+-");
  strcat(numberStringv2pteta,".dat");

  strcat(numberStringv3pteta, "./outputs/Fv3etapt-H+-");
  strcat(numberStringv3pteta,".dat");

  strcat(numberStringv4pteta, "./outputs/Fv4etapt-H+-");
  strcat(numberStringv4pteta,".dat");
      
  strcat(numberStringv5pteta, "./outputs/Fv5etapt-H+-");
  strcat(numberStringv5pteta,".dat");
      
  strcat(numberStringPhi, "./outputs/FNphi-H+-");
  strcat(numberStringPhi,".dat");
  
  strcat(numberStringy, "./outputs/FNy-H+-");
  strcat(numberStringy,".dat");
  
  strcat(numberStringv2, "./outputs/Fv2pt-H+-");
  strcat(numberStringv2,".dat");
  
  strcat(numberStringv3, "./outputs/Fv3pt-H+-");
  strcat(numberStringv3,".dat");
  
  strcat(numberStringv2reac, "./outputs/Fv2pt-react-H+-");
  strcat(numberStringv2reac,".dat");
  
  strcat(numberStringv3reac, "./outputs/Fv3pt-react-H+-");
  strcat(numberStringv3reac,".dat");
  
  strcat(numberStringv4reac, "./outputs/Fv4pt-react-H+-");
  strcat(numberStringv4reac,".dat");
  
  strcat(numberStringv5reac, "./outputs/Fv5pt-react-H+-");
  strcat(numberStringv5reac,".dat");
  
  strcat(numberStringv3r3, "./outputs/Fv3r3pt-H+-");
  strcat(numberStringv3r3,".dat");
  
  strcat(numberStringv4, "./outputs/Fv4pt-H+-");
  strcat(numberStringv4,".dat");
  
  strcat(numberStringv5, "./outputs/Fv5pt-H+-");
  strcat(numberStringv5,".dat");
  
  strcat(numberStringeta, "./outputs/FNeta-H+-");
  strcat(numberStringeta,".dat");
  
  strcat(numberStringv2eta, "./outputs/Fv2eta-H+-");
  strcat(numberStringv2eta,".dat");
  
  strcat(numberStringv3eta, "./outputs/Fv3eta-H+-");
  strcat(numberStringv3eta,".dat");
  
  strcat(numberStringv3r3eta, "./outputs/Fv3r3eta-H+-");
  strcat(numberStringv3r3eta,".dat");
  
  strcat(numberStringv4eta, "./outputs/Fv4eta-H+-");
  strcat(numberStringv4eta,".dat");
  
  strcat(numberStringv5eta, "./outputs/Fv5eta-H+-");
  strcat(numberStringv5eta,".dat");
  
  strcat(numberStringv6eta, "./outputs/Fv6eta-H+-");
  strcat(numberStringv6eta,".dat");
  
  strcat(numberStringv1eta, "./outputs/Fv1eta-H+-");
  strcat(numberStringv1eta,".dat");
  
  strcat(numberStringv2tot, "./outputs/Fv2tot-H+-");
  strcat(numberStringv2tot,".dat");
  
  strcat(numberStringv3tot, "./outputs/Fv3tot-H+-");
  strcat(numberStringv3tot,".dat");
  
  strcat(numberStringv3r3tot, "./outputs/Fv3r3tot-H+-");
  strcat(numberStringv3r3tot,".dat");
  
  strcat(numberStringv4tot, "./outputs/Fv4tot-H+-");
  strcat(numberStringv4tot,".dat");
  
  strcat(numberStringv5tot, "./outputs/Fv5tot-H+-");
  strcat(numberStringv5tot,".dat");
  
  strcat(numberStringv6tot, "./outputs/Fv6tot-H+-");
  strcat(numberStringv6tot,".dat");
  
  strcat(numberStringv1tot, "./outputs/Fv1tot-H+-");
  strcat(numberStringv1tot,".dat");
  

  strcat(numberStringpt, "./outputs/FNptExp-H+-");
  strcat(numberStringpt,".dat");

  strcat(numberStringv1pt, "./outputs/Fv1ptExp-H+-");
  strcat(numberStringv1pt,".dat");

  strcat(numberStringv2pt, "./outputs/Fv2ptExp-H+-");
  strcat(numberStringv2pt,".dat");

  strcat(numberStringv3pt, "./outputs/Fv3ptExp-H+-");
  strcat(numberStringv3pt,".dat");

  strcat(numberStringv4pt, "./outputs/Fv4ptExp-H+-");
  strcat(numberStringv4pt,".dat");
      
  strcat(numberStringv5pt, "./outputs/Fv5ptExp-H+-");
  strcat(numberStringv5pt,".dat");

  strcat(numberStringv6pt, "./outputs/Fv6ptExp-H+-");
  strcat(numberStringv6pt,".dat");

  FILE *y_file;
  char* y_name = numberStringy;
  y_file = fopen(y_name, "w");
  
  FILE *eta_file;
  char* eta_name = numberStringeta;
  eta_file = fopen(eta_name, "w");

  FILE *pteta_file;
  char* pteta_name = numberStringpteta;
  pteta_file = fopen(pteta_name, "w");
  
  FILE *v2_pteta_file;
  char* v2_pteta_name = numberStringv2pteta;
  v2_pteta_file = fopen(v2_pteta_name, "w");
  
  FILE *v3_pteta_file;
  char* v3_pteta_name = numberStringv3pteta;
  v3_pteta_file = fopen(v3_pteta_name, "w");
  
  FILE *v4_pteta_file;
  char* v4_pteta_name = numberStringv4pteta;
  v4_pteta_file = fopen(v4_pteta_name, "w");
  
  FILE *v5_pteta_file;
  char* v5_pteta_name = numberStringv5pteta;
  v5_pteta_file = fopen(v5_pteta_name, "w");
  
  FILE *v2_eta_file;
  char* v2_eta_name = numberStringv2eta;
  v2_eta_file = fopen(v2_eta_name, "w");

  FILE *v3_eta_file;
  char* v3_eta_name = numberStringv3eta;
  v3_eta_file = fopen(v3_eta_name, "w");

  FILE *v3r3_eta_file;
  char* v3r3_eta_name = numberStringv3r3eta;
  v3r3_eta_file = fopen(v3r3_eta_name, "w");

  FILE *v4_eta_file;
  char* v4_eta_name = numberStringv4eta;
  v4_eta_file = fopen(v4_eta_name, "w");

  FILE *v5_eta_file;
  char* v5_eta_name = numberStringv5eta;
  v5_eta_file = fopen(v5_eta_name, "w");

  FILE *v6_eta_file;
  char* v6_eta_name = numberStringv6eta;
  v6_eta_file = fopen(v6_eta_name, "w");

  FILE *v1_eta_file;
  char* v1_eta_name = numberStringv1eta;
  v1_eta_file = fopen(v1_eta_name, "w");

  FILE *v2_tot_file;
  char* v2_tot_name = numberStringv2tot;
  v2_tot_file = fopen(v2_tot_name, "w");

  FILE *v3_tot_file;
  char* v3_tot_name = numberStringv3tot;
  v3_tot_file = fopen(v3_tot_name, "w");

  FILE *v3r3_tot_file;
  char* v3r3_tot_name = numberStringv3r3tot;
  v3r3_tot_file = fopen(v3r3_tot_name, "w");

  FILE *v4_tot_file;
  char* v4_tot_name = numberStringv4tot;
  v4_tot_file = fopen(v4_tot_name, "w");

  FILE *v5_tot_file;
  char* v5_tot_name = numberStringv5tot;
  v5_tot_file = fopen(v5_tot_name, "w");

  FILE *v6_tot_file;
  char* v6_tot_name = numberStringv6tot;
  v6_tot_file = fopen(v6_tot_name, "w");

  FILE *v1_tot_file;
  char* v1_tot_name = numberStringv1tot;
  v1_tot_file = fopen(v1_tot_name, "w");

  FILE *p_file;
  char* p_name = numberString;
  p_file = fopen(p_name, "w");

  FILE *phi_file;
  char* phi_name = numberStringPhi;
  phi_file = fopen(phi_name, "w");

  FILE *v2_file;
  char* v2_name = numberStringv2;
  v2_file = fopen(v2_name, "w");

  FILE *v3_file;
  char* v3_name = numberStringv3;
  v3_file = fopen(v3_name, "w");

  FILE *v2_reac_file;
  char* v2_reac_name = numberStringv2reac;
  v2_reac_file = fopen(v2_reac_name, "w");

  FILE *v3_reac_file;
  char* v3_reac_name = numberStringv3reac;
  v3_reac_file = fopen(v3_reac_name, "w");

  FILE *v4_reac_file;
  char* v4_reac_name = numberStringv4reac;
  v4_reac_file = fopen(v4_reac_name, "w");

  FILE *v5_reac_file;
  char* v5_reac_name = numberStringv5reac;
  v5_reac_file = fopen(v5_reac_name, "w");

  FILE *v3r3_file;
  char* v3r3_name = numberStringv3r3;
  v3r3_file = fopen(v3r3_name, "w");

  FILE *v4_file;
  char* v4_name = numberStringv4;
  v4_file = fopen(v4_name, "w");

  FILE *v5_file;
  char* v5_name = numberStringv5;
  v5_file = fopen(v5_name, "w");


  FILE *pt_file;
  char* pt_name = numberStringpt;
  pt_file = fopen(pt_name, "w");
  
  FILE *v1_pt_file;
  char* v1_pt_name = numberStringv1pt;
  v1_pt_file = fopen(v1_pt_name, "w");
  
  FILE *v2_pt_file;
  char* v2_pt_name = numberStringv2pt;
  v2_pt_file = fopen(v2_pt_name, "w");
  
  FILE *v3_pt_file;
  char* v3_pt_name = numberStringv3pt;
  v3_pt_file = fopen(v3_pt_name, "w");
  
  FILE *v4_pt_file;
  char* v4_pt_name = numberStringv4pt;
  v4_pt_file = fopen(v4_pt_name, "w");
  
  FILE *v5_pt_file;
  char* v5_pt_name = numberStringv5pt;
  v5_pt_file = fopen(v5_pt_name, "w");
  
  FILE *v6_pt_file;
  char* v6_pt_name = numberStringv6pt;
  v6_pt_file = fopen(v6_pt_name, "w");
  
  cout << " done." << endl;

  slope = particleList[j].slope;
  switch (iphimax) 
    {
    case 4: p= gaulep4; w= gaulew4; break;
    case 8: p= gaulep8; w= gaulew8; break;
    case 10: p= gaulep10; w= gaulew10; break;
    case 12: p= gaulep12; w= gaulew12; break;
    case 16: p= gaulep16; w= gaulew16; break;
    case 20: p= gaulep20; w= gaulew20; break;
    case 48: p= gaulep48; w= gaulew48; break;
    default: fprintf(stderr,"specified number of phi-points not available\n"); exit(1);
    }

  cout << " filling phi array." << endl;
  for(iphi=0; iphi<iphimax; iphi++)
    {
      if ( iphi < iphimax/2 )
	{
	  phiArray[iphi] = (phimax-phimin)/2.*(1.-p[iphi])+phimin;
	}
      else
	{
	  phiArray[iphi] = (phimax-phimin)/2.*(1.+p[iphimax-iphi-1])+phimin;
	}
      //      cout << "phi[" << iphi << "]=" << phiArray[iphi] << ", " <<  phimax/2.*(1.-p[iphi]) << ", " 
      //   << (phimax)/2.*(1.+p[iphimax-iphi-1]) << endl;
    }

  cout << "starting..." << endl;

  double Psi6pSin = 0.;
  double Psi6pCos = 0.;
  double Psi6p;
  double Psi5pSin = 0.;
  double Psi5pCos = 0.;
  double Psi5p;
  double Psi4pSin = 0.;
  double Psi4pCos = 0.;
  double Psi4p;
  double Psi3pSin = 0.;
  double Psi3pCos = 0.;
  double Psi3p;
  double Psi2pSin = 0.;
  double Psi2pCos = 0.;
  double Psi2p;
  double Psi1pSin = 0.;
  double Psi1pCos = 0.;
  double Psi1p;
  for (iy=0; iy<iymax; iy++)
    {
      y =  particleList[j].y[iy];
      //      fprintf(stderr,"y=%f \n", y);
      sumpt=0.;
      for (ip = 0; ip<6; ip++)
	{
	  j = partid[MHALF+setOfNumbers[ip]];
	  //cout << j << " " << setOfNumbers[ip] << endl;
	  for (ipt=0; ipt<iptmax; ipt++)
	    {
	      pt = particleList[j].pt[ipt];
	      //fprintf(stderr,"pt=%f \n", pt);
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  phi = phiArray[iphi];
		  sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
		  // integrate over pt and phi
		  if (iphi<iphimax/2)
		    {
		      sumpt += w[iphi] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi6pSin += w[iphi] * ( sum*sin(6.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi6pCos += w[iphi] * ( sum*cos(6.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi5pSin += w[iphi] * ( sum*sin(5.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi5pCos += w[iphi] * ( sum*cos(5.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi4pSin += w[iphi] * ( sum*sin(4.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi4pCos += w[iphi] * ( sum*cos(4.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi3pSin += w[iphi] * ( sum*sin(3.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi3pCos += w[iphi] * ( sum*cos(3.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi2pSin += w[iphi] * ( sum*sin(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi2pCos += w[iphi] * ( sum*cos(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi1pSin += w[iphi] * ( sum*sin(1.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi1pCos += w[iphi] * ( sum*cos(1.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		    }
		  else
		    {
		      sumpt += w[iphimax-iphi-1] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi6pSin += w[iphimax-iphi-1] * ( sum*sin(6.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi6pCos += w[iphimax-iphi-1] * ( sum*cos(6.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi5pSin += w[iphimax-iphi-1] * ( sum*sin(5.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi5pCos += w[iphimax-iphi-1] * ( sum*cos(5.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi4pSin += w[iphimax-iphi-1] * ( sum*sin(4.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi4pCos += w[iphimax-iphi-1] * ( sum*cos(4.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi3pSin += w[iphimax-iphi-1] * ( sum*sin(3.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi3pCos += w[iphimax-iphi-1] * ( sum*cos(3.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi2pSin += w[iphimax-iphi-1] * ( sum*sin(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi2pCos += w[iphimax-iphi-1] * ( sum*cos(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi1pSin += w[iphimax-iphi-1] * ( sum*sin(1.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi1pCos += w[iphimax-iphi-1] * ( sum*cos(1.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		    }
		}
	    }
	}
      fprintf(y_file,"%f %f\n", y, sumpt/slope);
    }
  if (Psi1pCos<0.)
    Psi1p = (atan(Psi1pSin/Psi1pCos)+PI);
  else
    Psi1p = (atan(Psi1pSin/Psi1pCos));

  if (Psi2pCos<0.)
    Psi2p = 1./2.*(atan(Psi2pSin/Psi2pCos)+PI);
  else
    Psi2p = 1./2.*(atan(Psi2pSin/Psi2pCos));
    
  if (Psi3pCos<0.) 
    Psi3p = 1./3.*(atan(Psi3pSin/Psi3pCos)+PI);
  else
    Psi3p = 1./3.*(atan(Psi3pSin/Psi3pCos));
  
  if (Psi4pCos<0.) 
    Psi4p = 1./4.*(atan(Psi4pSin/Psi4pCos)+PI);
  else
    Psi4p = 1./4.*(atan(Psi4pSin/Psi4pCos));
  
  if (Psi5pCos<0.) 
    Psi5p = 1./5.*(atan(Psi5pSin/Psi5pCos)+PI);
  else
    Psi5p = 1./5.*(atan(Psi5pSin/Psi5pCos));

  if (Psi6pCos<0.) 
    Psi6p = 1./6.*(atan(Psi6pSin/Psi6pCos)+PI);
  else
    Psi6p = 1./6.*(atan(Psi6pSin/Psi6pCos));
  
  cout << "Psi2p=" << Psi2p << endl; 
  cout << "Psi2pSin=" << Psi2pSin << endl; 
  cout << "Psi2pCos=" << Psi2pCos << endl; 
  cout << "Psi3p=" << Psi3p << endl; 
  cout << "Psi3pSin=" << Psi3pSin << endl; 
  cout << "Psi3pCos=" << Psi3pCos << endl; 

  ofstream fout1("PsinH+-.dat",ios::out);
  fout1 << Psi1p << " " << Psi2p << " " << Psi3p << " " << Psi4p << " " << Psi5p << " " << Psi6p << endl;
  fout1.close();

  ofstream fout2("PsinCorrelatorsH+-.dat",ios::out);
  fout2 << cos(4.*(Psi2p-Psi4p)) << " "
	<< cos(8.*(Psi2p-Psi4p)) << " " 
	<< cos(12.*(Psi2p-Psi4p)) << " " 
	<< cos(6.*(Psi2p-Psi3p)) << " "
	<< cos(6.*(Psi2p-Psi6p)) << " " 
	<< cos(6.*(Psi3p-Psi6p)) << " "
	<< cos(12.*(Psi3p-Psi4p)) << " "
	<< cos(10.*(Psi2p-Psi5p)) << " "
	<< cos(2.*Psi2p+3.*Psi3p-5.*Psi5p) << " "
	<< cos(2.*Psi2p+4.*Psi4p-6.*Psi6p) << " "
	<< cos(2.*Psi2p-6.*Psi3p+4.*Psi4p) << " " 
	<< cos(-8.*Psi2p+3.*Psi3p+5.*Psi5p) << " "
	<< cos(-10.*Psi2p+4.*Psi4p+6.*Psi6p) << " "
	<< cos(-10.*Psi2p+6.*Psi3p+4.*Psi4p) << " "
	<< sin(3.*(Psi2p-Psi3p))*sin(5.*(Psi2p-Psi5p)) << " "
	<< cos(3.*(Psi2p-Psi3p))*cos(5.*(Psi2p-Psi5p)) << " "
	<< sin(4.*(Psi2p-Psi4p))*sin(6.*(Psi2p-Psi6p)) << " "
	<< cos(4.*(Psi2p-Psi4p))*cos(6.*(Psi2p-Psi6p)) << " "
	<< sin(6.*(Psi2p-Psi3p))*sin(4.*(Psi2p-Psi4p)) << " "
	<< cos(6.*(Psi2p-Psi3p))*cos(4.*(Psi2p-Psi4p)) 
	<<endl;
  fout2.close();

  //pt spectra
  for (ipt=0; ipt<iptmax; ipt++)
    {
      sumpt = 0.;
      int countY = 0;
      pt = particleList[j].pt[ipt];
      for (ip = 0; ip<6; ip++)
	{
	  j = partid[MHALF+setOfNumbers[ip]];
	  for (iy=0; iy<iymax; iy++)
	    {
	      y = iy*deltaY-ymax;
	      if (fabs(y)>ymaxIntegral) continue;
	      countY++;
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  phi = phiArray[iphi];
		  sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
		  if (iphi<iphimax/2)
		    {
		      if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
			{
			  sumpt += 0.5 * w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
			}
		      else
			{
			  sumpt += w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
			}
		      //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphi], iphi, p[iphi]);
		    }
		  else
		    {
		      if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
			{
			  sumpt += 0.5 * w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      else
			{
			  sumpt += w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphimax-iphi-1], iphi, -p[iphimax-iphi-1]);
		    }
		}
	    }
	}
      //fprintf(stderr,"pt=%f, sumpt=%f\n", pt, sumpt / (2.*ymaxIntegral));
      fprintf(p_file,"%e %e \n", pt, sumpt / (2.*ymaxIntegral));
      //fprintf(stderr,"pt=%f, sumv%d=%f\n", pt, 2, sumv[2]/sumpt);
    }
  
  cout << "pT done" << endl;


//   // v2, v3 and v4
//   for (ipt=0; ipt<iptmax; ipt++)
//     {
//       sumpt = 0.;
//       sumv[2] = 0.;
//       sumv[3] = 0.;
//       sumv3r3 = 0.;
//       sumv[4] = 0.;
//       int countY = 0;
//       for (j=1; j<=3; j++)
// 	{
// 	  pt = particleList[j].pt[ipt];
// 	  ymaxIntegral = 0.5*log((sqrt(pt*pt*cosh(etaMaxIntegral)*cosh(etaMaxIntegral)+m*m)
// 				  +pt*sinh(etaMaxIntegral))/(sqrt(pt*pt*cosh(etaMaxIntegral)*cosh(etaMaxIntegral)+m*m)-pt*sinh(etaMaxIntegral)));
// 	  iymaxIntegral = floor(2.*ymaxIntegral/deltaY);
// 	  for (iy=0; iy<iymax; iy++)
// 	    {
// 	      y = iy*deltaY-ymax;
// 	      if (fabs(y)>ymaxIntegral) continue;
// 	      countY++;
// 	      for (iphi=0; iphi<iphimax; iphi++)
// 		{
// 		  phi = phiArray[iphi];
// 		  sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
		  
// 		  //       if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
// // 		    {
// // 		      sumpt += 0.5 * deltaPhi * ( sum ) * deltaY  * 2.*PI/(phimax-phimin);
// // 		      sumv[2] += 0.5 * deltaPhi * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * 2.*PI/(phimax-phimin);
// // 		      sumv[3] += 0.5 * deltaPhi * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * 2.*PI/(phimax-phimin);
// // 		      sumv3r3 += 0.5 * deltaPhi * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * 2.*PI/(phimax-phimin);
// // 		      sumv[4] += 0.5 * deltaPhi * ( sum*cos(4.*(phi-phimin)) ) * deltaY * 2.*PI/(phimax-phimin);
// // 		    }
// // 		  else
// // 		    {
// // 		      sumpt += deltaPhi * ( sum ) * deltaY * 2.*PI/(phimax-phimin);
// // 		      sumv[2] += deltaPhi * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * 2.*PI/(phimax-phimin);
// // 		      sumv[3] += deltaPhi * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * 2.*PI/(phimax-phimin);
// // 		      sumv3r3 += deltaPhi * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * 2.*PI/(phimax-phimin);
// // 		      sumv[4] += deltaPhi * ( sum*cos(4.*(phi-phimin)) ) * deltaY * 2.*PI/(phimax-phimin);
// // 		    }
		
// 		  if (iphi<iphimax/2)
// 		    {
// 		      if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
// 			{
// 			  sumpt += 0.5 * w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
// 			  sumv[2] += 0.5 * w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			  sumv[3] += 0.5 * w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			  sumv3r3 += 0.5 * w[iphi] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			  sumv[4] += 0.5 * w[iphi] * ( sum*cos(4.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			}
// 		      else
// 			{
// 			  sumpt += w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
// 			  sumv[2] += w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			  sumv[3] += w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			  sumv3r3 += w[iphi] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			  sumv[4] += w[iphi] * ( sum*cos(4.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			}
// 		      //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphi], iphi, p[iphi]);
// 		    }
// 		  else
// 		    {
// 		      if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
// 			{
// 			  sumpt += 0.5 * w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			  sumv[2] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			  sumv[3] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			  sumv3r3 += 0.5 * w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			  sumv[4] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			}
// 		      else
// 			{
// 			  sumpt += w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			  sumv[2] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			  sumv[3] += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			  sumv3r3 += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			  sumv[4] += w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
// 			}
// 		      //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphimax-iphi-1], iphi, -p[iphimax-iphi-1]);
// 		    }
// 		}
// 	    }
// 	}
//       //fprintf(stderr,"pt=%f, sumpt=%f\n", pt, sumpt / (2.*ymaxIntegral));
//       //fprintf(stderr,"pt=%f, sumv%d=%f\n", pt, 2, sumv[2]/sumpt);
//       fprintf(v2_file,"%e %e \n", pt, sumv[2]/sumpt);
//       fprintf(v3_file,"%e %e \n", pt, sumv[3]/sumpt);
//       fprintf(v3r3_file,"%e %e \n", pt, sumv3r3/sumpt);
//       fprintf(v4_file,"%e %e \n", pt, sumv[4]/sumpt);
//     }

  cout << "v2, v3, v4, v5(pT)" << endl;
  // v2, v3 and v4, v5
  for (ipt=0; ipt<iptmax; ipt++)
    {
      //cout << "pt=" << pt << endl;
      sumpt = 0.;
      sumv[2] = 0.;
      sumv[3] = 0.;
      sumv[4] = 0.;
      sumv[5] = 0.;
      sumv3r3 = 0.;
      int countY = 0;
      pt = particleList[j].pt[ipt];
      //cout << "doing v2 v3 v4 at pT=" << pt << endl;
      for (ip = 0; ip<6; ip++)
	{
	  j = partid[MHALF+setOfNumbers[ip]];
	  m = particleList[j].mass;
	  //cout << "particle=" << ip << endl;
	  for (iy=0; iy<iymax; iy++)
	    {
	      y = iy*deltaY-ymax;
	      ymaxIntegral = 0.5*log((sqrt(pt*pt*cosh(etaMaxIntegral)*cosh(etaMaxIntegral)+m*m)
				      +pt*sinh(etaMaxIntegral))/(sqrt(pt*pt*cosh(etaMaxIntegral)*cosh(etaMaxIntegral)+m*m)-pt*sinh(etaMaxIntegral)));
	      if (fabs(y)>ymaxIntegral) continue;
	      countY++;
	      //cout << "iy=" << iy << endl;
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  //cout << "iphi=" << iphi << endl;
		  phi = phiArray[iphi];
		  sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
		  if (iphi<iphimax/2)
		    {
		      if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
			{
			  sumpt += 0.5 * w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
			  sumv[2] += 0.5 * w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[3] += 0.5 * w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  //sumv3r3 += 0.5 * w[iphi] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[4] += 0.5 * w[iphi] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[5] += 0.5 * w[iphi] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      else
			{
			  sumpt += w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
			  sumv[2] += w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[3] += w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  //sumv3r3 += w[iphi] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[4] += w[iphi] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[5] += w[iphi] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphi], iphi, p[iphi]);
		    }
		  else
		    {
		      if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
			{
			  sumpt += 0.5 * w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[2] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[3] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  //sumv3r3 += 0.5 * w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[4] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[5] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      else
			{
			  sumpt += w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[2] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[3] += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  //sumv3r3 += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[4] += w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[5] += w[iphimax-iphi-1] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphimax-iphi-1], iphi, -p[iphimax-iphi-1]);
		    }
		}
	    }
	}
      //cout << "ok" << endl;
      //fprintf(stderr,"pt=%f, sumpt=%f\n", pt, sumpt / (2.*ymaxIntegral));
      //fprintf(p_file,"%e %e \n", pt, sumpt / (2.*ymaxIntegral));
      //fprintf(stderr,"pt=%f, sumv%d=%f\n", pt, 2, sumv[2]/sumpt);
      fprintf(v2_file,"%e %e \n", pt, sumv[2]/sumpt);
      //cout << "v2file ok" << endl;
      fprintf(v4_file,"%e %e \n", pt, sumv[4]/sumpt);
      //cout << "v4file ok" << endl;
      fprintf(v3_file,"%e %e \n", pt, sumv[3]/sumpt);
      //cout << "v3file ok" << endl;
      fprintf(v5_file,"%e %e \n", pt, sumv[5]/sumpt);
      fprintf(v3r3_file,"%e %e \n", pt, sumv3r3/sumpt);
      //cout << "v3r3file ok" << endl;
    }

  cout << "v2v3v4v5(pT) done" << endl;

  // v2, v3 and v4 in reaction plane
  for (ipt=0; ipt<iptmax; ipt++)
    {
      //cout << "pt=" << pt << endl;
      sumpt = 0.;
      sumv[2] = 0.;
      sumv[4] = 0.;
      sumv[5] = 0.;
      sumv[3] = 0.;
      sumv3r3 = 0.;
      int countY = 0;
      pt = particleList[j].pt[ipt];
      //cout << "doing v2 v3 v4 at pT=" << pt << endl;
      for (ip = 0; ip<6; ip++)
	{
	  j = partid[MHALF+setOfNumbers[ip]];
	  m = particleList[j].mass;
	  //cout << "particle=" << ip << endl;
	  for (iy=0; iy<iymax; iy++)
	    {
	      y = iy*deltaY-ymax;
	      ymaxIntegral = 0.5*log((sqrt(pt*pt*cosh(etaMaxIntegral)*cosh(etaMaxIntegral)+m*m)
				      +pt*sinh(etaMaxIntegral))/(sqrt(pt*pt*cosh(etaMaxIntegral)*cosh(etaMaxIntegral)+m*m)-pt*sinh(etaMaxIntegral)));
	      if (fabs(y)>ymaxIntegral) continue;
	      countY++;
	      //cout << "iy=" << iy << endl;
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  //cout << "iphi=" << iphi << endl;
		  phi = phiArray[iphi];
		  sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
		  if (iphi<iphimax/2)
		    {
		      if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
			{
			  sumpt += 0.5 * w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
			  sumv[2] += 0.5 * w[iphi] * ( sum*cos(2.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[3] += 0.5 * w[iphi] * ( sum*cos(3.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[4] += 0.5 * w[iphi] * ( sum*cos(4.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[5] += 0.5 * w[iphi] * ( sum*cos(5.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      else
			{
			  sumpt += w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
			  sumv[2] += w[iphi] * ( sum*cos(2.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[3] += w[iphi] * ( sum*cos(3.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[4] += w[iphi] * ( sum*cos(4.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[5] += w[iphi] * ( sum*cos(5.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphi], iphi, p[iphi]);
		    }
		  else
		    {
		      if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
			{
			  sumpt += 0.5 * w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[2] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[3] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[4] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[5] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(5.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      else
			{
			  sumpt += w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[2] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[3] += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[4] += w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[5] += w[iphimax-iphi-1] * ( sum*cos(5.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphimax-iphi-1], iphi, -p[iphimax-iphi-1]);
		    }
		}
	    }
	}
      fprintf(v2_reac_file,"%e %e \n", pt, sumv[2]/sumpt);
      fprintf(v3_reac_file,"%e %e \n", pt, sumv[3]/sumpt);
      fprintf(v4_reac_file,"%e %e \n", pt, sumv[4]/sumpt);
      fprintf(v5_reac_file,"%e %e \n", pt, sumv[5]/sumpt);
    }

  cout << "v2v3v4(pT) in reaction plane done" << endl;

  //pseudorapidity-pteta
  deltaEta=deltaY;
  etamax=ymax;
  ietamax=floor(2.*etamax/deltaEta);
  for (ieta=0; ieta<ietamax; ieta++)
    {
      eta = ieta*deltaEta-etamax+deltaEta/2.;
      sumpt=0.;
      sumv[2]=0.;
      sumv[3]=0.;
      sumv3r3=0.;
      sumv[4]=0.;
      sumv[5]=0.;
      for (ipt=0; ipt<iptmax; ipt++)
	{
	  sumpteta=0.;
	  sumvpteta[2]=0.;
	  sumvpteta[3]=0.;
	  sumvpteta[4]=0.;
	  sumvpteta[5]=0.;
	  for (ip = 0; ip<6; ip++)
	    {
	      j = partid[MHALF+setOfNumbers[ip]];
	      m = particleList[j].mass;
	      pt = particleList[j].pt[ipt];
	      y = 0.5*log((sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)+pt*sinh(eta))/(sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)-pt*sinh(eta)));
	      
	      for(i=0; i<iymax-1; i++) // find closest iy
		{
		  if ( particleList[j].y[i]< y &&  particleList[j].y[i+1] > y )
		    {
		      if ( fabs(particleList[j].y[i]-y) <  fabs(particleList[j].y[i+1]-y)) 
			iy = i;
		      else
			iy = i+1;
		    }
		}
	      
	      mt = sqrt(pt*pt+m*m);
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  phi = phiArray[iphi];
		  px = pt*cos(phi);
		  py = pt*sin(phi);
		  sum = particleList[j].dNdydptdphi[iy][ipt][iphi]*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
		  if (iphi<iphimax/2)
		    {
		      sumpteta += w[iphi] * ( sum ) * phiDiff * 2.*PI/(phimax-phimin);
		      sumvpteta[2] += w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * phiDiff * 2.*PI/(phimax-phimin);
		      sumvpteta[3] += w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * phiDiff * 2.*PI/(phimax-phimin);
		      sumvpteta[4] += w[iphi] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * phiDiff * 2.*PI/(phimax-phimin);		
		      sumvpteta[5] += w[iphi] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * phiDiff * 2.*PI/(phimax-phimin);		
		    }
		  else
		    {
		      sumpteta += w[iphimax-iphi-1] * ( sum ) * phiDiff * 2.*PI/(phimax-phimin);
		      sumvpteta[2] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * phiDiff * 2.*PI/(phimax-phimin);
		      sumvpteta[3] += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * phiDiff * 2.*PI/(phimax-phimin);
		      sumvpteta[4] += w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * phiDiff * 2.*PI/(phimax-phimin);
		      sumvpteta[5] += w[iphimax-iphi-1] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * phiDiff * 2.*PI/(phimax-phimin);
		    }
		}
	    }
	  fprintf(pteta_file,"%f %f %e\n", eta, pt, sumpteta);
	  fprintf(v2_pteta_file,"%f %f %e\n", eta, pt, sumvpteta[2]/sumpteta);
	  fprintf(v3_pteta_file,"%f %f %e\n", eta, pt, sumvpteta[3]/sumpteta);
	  fprintf(v4_pteta_file,"%f %f %e\n", eta, pt, sumvpteta[4]/sumpteta);
	  fprintf(v5_pteta_file,"%f %f %e\n", eta, pt, sumvpteta[5]/sumpteta);
	}
      fprintf(pteta_file,"\n");
      fprintf(v2_pteta_file,"\n");
      fprintf(v3_pteta_file,"\n");
      fprintf(v4_pteta_file,"\n");
      fprintf(v5_pteta_file,"\n");
    }
  
  //pseudorapidity
  deltaEta=deltaY;
  etamax=ymax;
  ietamax=floor(2.*etamax/deltaEta);
  for (ieta=0; ieta<ietamax; ieta++)
    {
      eta = ieta*deltaEta-etamax+deltaEta/2.;
      sumpt=0.;
      sumv[2]=0.;
      sumv[3]=0.;
      sumv3r3=0.;
      sumv[4]=0.;
      sumv[5]=0.;
      sumv[6]=0.;
      for (ip = 0; ip<6; ip++)
	{
	  j = partid[MHALF+setOfNumbers[ip]];
	  m = particleList[j].mass;
	  for (ipt=0; ipt<iptmax; ipt++)
	    {
	      pt = particleList[j].pt[ipt];
	      //ATLAS cut
	      //	      if(pt>3. || pt<2.) continue;
	      y = 0.5*log((sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)+pt*sinh(eta))/(sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)-pt*sinh(eta)));
	      
	      for(i=0; i<iymax-1; i++) // find closest iy
		{
		  if ( particleList[j].y[i]< y &&  particleList[j].y[i+1] > y )
		    {
		      //		      if ( fabs(particleList[j].y[i]-y) <  fabs(particleList[j].y[i+1]-y)) 
			iy = i;
			// else
			//iy = i+1;
		    }
		}
	      
	      mt = sqrt(pt*pt+m*m);
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  phi = phiArray[iphi];
		  px = pt*cos(phi);
		  py = pt*sin(phi);
	
 		  sum = (
 			 (particleList[j].dNdydptdphi[iy+1][ipt][iphi])*(y-particleList[j].y[iy])/(particleList[j].y[iy+1]-particleList[j].y[iy])
 			 +(particleList[j].dNdydptdphi[iy][ipt][iphi])*(1.-(y-particleList[j].y[iy])/(particleList[j].y[iy+1]-particleList[j].y[iy]))
			 )
 		    *sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
		  //interpolation in eta
		  //		  sum = (particleList[j].dNdydptdphi[iy][ipt][iphi])*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
		   
		  
		  if (iphi<iphimax/2)
		    {
		      sumpt += w[iphi] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[1] += w[iphi] * ( sum*cos(1.*(phi-phimin-Psi1p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[2] += w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[3] += w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[4] += w[iphi] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[5] += w[iphi] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[6] += w[iphi] * ( sum*cos(6.*(phi-phimin-Psi6p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		    }
		  else
		    {
		      sumpt += w[iphimax-iphi-1] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[1] += w[iphimax-iphi-1] * ( sum*cos(1.*(phi-phimin-Psi1p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[2] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[3] += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[4] += w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[5] += w[iphimax-iphi-1] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[6] += w[iphimax-iphi-1] * ( sum*cos(6.*(phi-phimin-Psi6p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		    }
		}
	    }
	}	
      //   fprintf(stderr,"eta=%f, sumpt=%f\n", eta, sumpt/slope);
      fprintf(eta_file,"%f %f\n", eta, sumpt/slope);
      fprintf(v1_eta_file,"%f %f\n", eta, sumv[1]/sumpt);
      fprintf(v2_eta_file,"%f %f\n", eta, sumv[2]/sumpt);
      fprintf(v3_eta_file,"%f %f\n", eta, sumv[3]/sumpt);
      fprintf(v4_eta_file,"%f %f\n", eta, sumv[4]/sumpt);
      fprintf(v5_eta_file,"%f %f\n", eta, sumv[5]/sumpt);
      fprintf(v6_eta_file,"%f %f\n", eta, sumv[6]/sumpt);
    }

  cout << "pseudo-rapidity done" << endl;

  //pt integrated over eta
  int countY = 0;
  deltaEta=deltaY;
  etamax=1.; //ATLAS
  ietamax=floor(2.*etamax/deltaEta);
  
  for (ipt=0; ipt<iptmax; ipt++)
    {
      sumpt=0.;
      sumv[1]=0.;
      sumv[2]=0.;
      sumv[3]=0.;
      sumv[4]=0.;
      sumv[5]=0.;
      sumv[6]=0.;
      for (ip = 0; ip<6; ip++)
	{
	  for (ieta=0; ieta<ietamax; ieta++)
	    {
	      eta = ieta*deltaEta-etamax+deltaEta/2.;
	      j = partid[MHALF+setOfNumbers[ip]];
	      m = particleList[j].mass;
	      pt = particleList[j].pt[ipt];
	      y = 0.5*log((sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)+pt*sinh(eta))/(sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)-pt*sinh(eta)));
	      
	      for(i=0; i<iymax-1; i++) // find closest iy
		{
		  if ( particleList[j].y[i]< y && particleList[j].y[i+1] > y )
		    {
		      // if ( fabs(particleList[j].y[i]-y) < fabs(particleList[j].y[i+1]-y)) 
			iy = i;
			//	      else
			//	iy = i+1;
		    }
		}
	      
	      mt = sqrt(pt*pt+m*m);
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  phi = phiArray[iphi];
		  px = pt*cos(phi);
		  py = pt*sin(phi);
		  sum = (
			 (particleList[j].dNdydptdphi[iy+1][ipt][iphi])*(y-particleList[j].y[iy])/(particleList[j].y[iy+1]-particleList[j].y[iy])
			 +(particleList[j].dNdydptdphi[iy][ipt][iphi])*(1.-(y-particleList[j].y[iy])/(particleList[j].y[iy+1]-particleList[j].y[iy]))
			 )
		    *sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
			 //interpolation in eta

		  if (iphi<iphimax/2)
		    {
		      if (ieta==0 || ieta==ietamax) // trapezoid rule: only half the edges
			{
			  sumpt += 0.5 * w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
			  sumv[1] += 0.5 * w[iphi] * ( sum*cos(1.*(phi-phimin-Psi1p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[2] += 0.5 * w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[3] += 0.5 * w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[4] += 0.5 * w[iphi] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[5] += 0.5 * w[iphi] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[6] += 0.5 * w[iphi] * ( sum*cos(6.*(phi-phimin-Psi6p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      else
			{
			  sumpt += w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
			  sumv[1] += w[iphi] * ( sum*cos(1.*(phi-phimin-Psi1p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[2] += w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[3] += w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[4] += w[iphi] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[5] += w[iphi] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[6] += w[iphi] * ( sum*cos(6.*(phi-phimin-Psi6p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		    }
		  else
		    {
		      if (ieta==0 || ieta==ietamax) // trapezoid rule: only half the edges
			{
			  sumpt += 0.5 * w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[1] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(1.*(phi-phimin-Psi1p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[2] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[3] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[4] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[5] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[6] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(6.*(phi-phimin-Psi6p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      else
			{
			  sumpt += w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[1] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi1p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[2] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[3] += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[4] += w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[5] += w[iphimax-iphi-1] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[6] += w[iphimax-iphi-1] * ( sum*cos(6.*(phi-phimin-Psi6p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		    }
		}
	    }
	}
      fprintf(pt_file,"%e %e\n", pt, sumpt/(2*etamax));
      fprintf(v1_pt_file,"%e %e\n", pt, sumv[1]/sumpt);
      fprintf(v2_pt_file,"%e %e\n", pt, sumv[2]/sumpt);
      fprintf(v3_pt_file,"%e %e\n", pt, sumv[3]/sumpt);
      fprintf(v4_pt_file,"%e %e\n", pt, sumv[4]/sumpt);
      fprintf(v5_pt_file,"%e %e\n", pt, sumv[5]/sumpt);
      fprintf(v6_pt_file,"%e %e\n", pt, sumv[6]/sumpt);
    }


  //total v2, v3 and v4
  deltaEta=deltaY;
  etamax=ymax;
  ietamax=floor(2.*etamax/deltaEta);
  double sumtot =0.;
  double sumtotv1 =0.;
  double sumtotv2 =0.;
  double sumtotv4 =0.;
  double sumtotv5 =0.;
  double sumtotv6 =0.;
  double sumtotv3 =0.;
  double sumtotv3r3 =0.;
  int eta1=0;
  for (ip = 0; ip<6; ip++)
    {
      j = partid[MHALF+setOfNumbers[ip]];
      m = particleList[j].mass;
      for (ieta=0; ieta<ietamax; ieta++)
	{
	  eta = ieta*deltaEta-etamax+deltaEta/2.;
	  if (fabs(eta)>1.) continue; // experimental cut
	  eta1+=1;
	  sumpt=0.;
	  sumv[1]=0.;
	  sumv[2]=0.;
	  sumv[3]=0.;
	  sumv[4]=0.;
	  sumv[5]=0.;
	  sumv[6]=0.;
	  for (ipt=0; ipt<iptmax; ipt++)
	    {
	      pt = particleList[j].pt[ipt];
	      if (pt < 1.5 || pt > 4.) continue; // experimental cut 
	      y = 0.5*log((sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)+pt*sinh(eta))/(sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)-pt*sinh(eta)));
	      
	      for(i=0; i<iymax-1; i++) // find closest iy
		{
		  if ( particleList[j].y[i]< y &&  particleList[j].y[i+1] > y )
		    {
		      if ( fabs(particleList[j].y[i]-y) <  fabs(particleList[j].y[i+1]-y)) 
			iy = i;
		      else
			iy = i+1;
		    }
		}
	      
	      mt = sqrt(pt*pt+m*m);
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  phi = phiArray[iphi];
		  px = pt*cos(phi);
		  py = pt*sin(phi);
		  sum = particleList[j].dNdydptdphi[iy][ipt][iphi]*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
		  if (iphi<iphimax/2)
		    {
		      sumpt += w[iphi] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[1] += w[iphi] * ( sum*cos(1.*(phi-phimin-Psi1p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[2] += w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[3] += w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[4] += w[iphi] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[5] += w[iphi] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[6] += w[iphi] * ( sum*cos(6.*(phi-phimin-Psi6p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		    }
		  else
		    {
		      sumpt += w[iphimax-iphi-1] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[1] += w[iphimax-iphi-1] * ( sum*cos(1.*(phi-phimin-Psi1p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[2] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[3] += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[4] += w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin-Psi4p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[5] += w[iphimax-iphi-1] * ( sum*cos(5.*(phi-phimin-Psi5p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[6] += w[iphimax-iphi-1] * ( sum*cos(6.*(phi-phimin-Psi6p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		    }
		}
	    }
		  
	  sumtot += sumpt * deltaEta;
	  if (eta1 == 1)
	    sumtot -= 0.5 * sumpt * deltaEta;
	  sumtotv1 += sumv[1] * deltaEta;
	  if (eta1 == 1)
	    sumtotv1 -= 0.5 * sumv[1] * deltaEta;
	  sumtotv2 += sumv[2] * deltaEta;
	  if (eta1 == 1)
	    sumtotv2 -= 0.5 * sumv[2] * deltaEta;
	  sumtotv3 += sumv[3] * deltaEta;
	  if (eta1 == 1)
	    sumtotv3 -= 0.5 * sumv[3] * deltaEta;
	  sumtotv4 += sumv[4] * deltaEta;
	  if (eta1 == 1)
	    sumtotv4 -= 0.5 * sumv[4] * deltaEta;
	  sumtotv5 += sumv[5] * deltaEta;
	  if (eta1 == 1)
	    sumtotv5 -= 0.5 * sumv[5] * deltaEta;
	  sumtotv6 += sumv[6] * deltaEta;
	  if (eta1 == 1)
	    sumtotv6 -= 0.5 * sumv[6] * deltaEta;
	}
      sumtot -= 0.5 * sumpt * deltaEta;
      sumtotv1 -= 0.5 * sumv[1] * deltaEta;
      sumtotv2 -= 0.5 * sumv[2] * deltaEta;
      sumtotv3 -= 0.5 * sumv[3] * deltaEta;
      sumtotv4 -= 0.5 * sumv[4] * deltaEta;
      sumtotv5 -= 0.5 * sumv[5] * deltaEta;
      sumtotv6 -= 0.5 * sumv[6] * deltaEta;
    }
  
  sumtotv1/=sumtot;
  sumtotv2/=sumtot;
  sumtotv3/=sumtot;
  sumtotv4/=sumtot;
  sumtotv5/=sumtot;
  sumtotv6/=sumtot;

  fprintf(v1_tot_file,"%f\n", sumtotv1);
  fprintf(v2_tot_file,"%f\n", sumtotv2);
  fprintf(v3_tot_file,"%f\n", sumtotv3);
  fprintf(v4_tot_file,"%f\n", sumtotv4);
  fprintf(v5_tot_file,"%f\n", sumtotv5);
  fprintf(v6_tot_file,"%f\n", sumtotv6);
  
  fprintf(stderr,"done.\n");
  cout << "done" << endl;
  
  fclose(p_file);
  fclose(v1_eta_file);
  fclose(v2_eta_file);
  fclose(v3_eta_file);
  fclose(v4_eta_file);
  fclose(v5_eta_file);
  fclose(v6_eta_file);
  fclose(v2_pteta_file);
  fclose(v3_pteta_file);
  fclose(v4_pteta_file);
  fclose(v5_pteta_file);
  fclose(v1_pt_file);
  fclose(v2_pt_file);
  fclose(v3_pt_file);
  fclose(v4_pt_file);
  fclose(v5_pt_file);
  fclose(v6_pt_file);
  fclose(v1_tot_file);
  fclose(v2_tot_file);
  fclose(v3_tot_file);
  fclose(v4_tot_file);
  fclose(v5_tot_file);
  fclose(v6_tot_file);
  fclose(v2_reac_file);
  fclose(v3_reac_file);
  fclose(v4_reac_file);
  fclose(v5_reac_file);
  fclose(eta_file);
  fclose(y_file);
  fclose(v2_file);
  fclose(v3_file);
  fclose(v3r3_file);
  fclose(v4_file);
  fclose(v5_file);
}

//compute pair yield and yield as functions of delta eta and delta phi
void Freeze::ComputeCorrelations(InitData* DATA, double ptmax)
{
  char *numberString;
  char buf[10];
  double *p, *w;		        // pointing to data for Gaussian integration in phi 
  double slope, slope1, slope2, fleft1, fleft2, fright1, fright2;
  int ipt, iphi, iymax, iy, iptmax, iphimax, ieta, ietamax;
  int i,j, iymaxIntegral;
  double pt, phi, px, py, y, deltaPT, deltaPhi, deltaY, ymax, phimin, phimax, sum, sumpt, sumpt2, sumv[5], m, sumv3r3;
  double phiOffs, phiDiff, ymaxIntegral, eta, deltaEta, etamax, mt, etaMaxIntegral;
  int returnValue;
  int number;
  int setOfNumbers[6] = {211,-211,2212,-2212,321,-321}; //these are the charged hadrons to include 
  double v2Delta, v3Delta;
  int ip;
  
  // set some parameters

  number = setOfNumbers[0]; //use pion to get settings

  j = 1;
  fprintf(stderr,"Doing correlations\n");
  phimax = particleList[j].phimax;
  phimin = particleList[j].phimin;
  
  phiOffs = 0.5 * ( phimin + phimax );
  phiDiff = 0.5 * ( phimax - phimin );
  iphimax = particleList[j].nphi;
  iptmax = particleList[j].npt;

  if (iptmax != 15) 
    {
      fprintf(stderr,"only iptmax=15 possible. you picked %d. exiting.\n", iptmax);
      exit(1);
    }

  ymax = particleList[j].ymax;
  deltaY = particleList[j].deltaY;
  ymaxIntegral = 0.5-deltaY/2.;
  etaMaxIntegral = 1.3; // for v2 integration
  iymax = particleList[j].ny;
  iymaxIntegral = floor(2.*ymaxIntegral/deltaY);
  iptmax = particleList[j].npt;
  iphimax = particleList[j].nphi;


  // set up file names
  numberString = util->char_malloc(30);
  
  if( chdir("./outputs") != 0 )
    {
      fprintf(stderr,"directory \"outputs\" does not exist. Exiting.\n");
      exit(1);
    }
  else returnValue=chdir("..");
  
  
  strcat(numberString, "./outputs/correlation");
  strcat(numberString,".dat");
  
  
  FILE *p_file;
  char* p_name = numberString;
  p_file = fopen(p_name, "w");

  slope = particleList[j].slope;
  switch (iphimax) 
    {
    case 4: p= gaulep4; w= gaulew4; break;
    case 8: p= gaulep8; w= gaulew8; break;
    case 10: p= gaulep10; w= gaulew10; break;
    case 12: p= gaulep12; w= gaulew12; break;
    case 16: p= gaulep16; w= gaulew16; break;
    case 20: p= gaulep20; w= gaulew20; break;
    case 48: p= gaulep48; w= gaulew48; break;
    default: fprintf(stderr,"specified number of phi-points not available\n"); exit(1);
    }

  phiArray = util->vector_malloc(iphimax+1);
  for(iphi=0; iphi<iphimax; iphi++)
    {
      if ( iphi < iphimax/2 )
	{
	  phiArray[iphi] = (phimax-phimin)/2.*(1.-p[iphi])+phimin;
	}
      else
	{
	  phiArray[iphi] = (phimax-phimin)/2.*(1.+p[iphimax-iphi-1])+phimin;
	}
      //      cout << "phi[" << iphi << "]=" << phiArray[iphi] << ", " <<  phimax/2.*(1.-p[iphi]) << ", " 
      //   << (phimax)/2.*(1.+p[iphimax-iphi-1]) << endl;
    }

  cout << "starting..." << endl;

  double Psi3pSin = 0.;
  double Psi3pCos = 0.;
  double Psi3p;
  double Psi2pSin = 0.;
  double Psi2pCos = 0.;
  double Psi2p;

//   //compute event plane to subtract v_2
//   for (iy=0; iy<iymax; iy++)
//     {
//       y =  particleList[j].y[iy];
//       //      fprintf(stderr,"y=%f \n", y);
//       sumpt=0.;
//       for (ip = 0; ip<6; ip++)
// 	{
// 	  j = partid[MHALF+setOfNumbers[ip]];
// 	  //cout << j << " " << setOfNumbers[ip] << endl;
// 	  for (ipt=0; ipt<iptmax; ipt++)
// 	    {
// 	      pt = particleList[j].pt[ipt];
// 	      //fprintf(stderr,"pt=%f \n", pt);
// 	      for (iphi=0; iphi<iphimax; iphi++)
// 		{
// 		  phi = phiArray[iphi];
// 		  sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
// 		  // integrate over pt and phi
// 		  if (iphi<iphimax/2)
// 		    {
// 		      sumpt += w[iphi] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
// 		      Psi2pSin += w[iphi] * ( sum*pt*sin(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
// 		      Psi2pCos += w[iphi] * ( sum*pt*cos(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
// 		      //fprintf(stderr,"sumpt=%f \n", sumpt);
// 		    }
// 		  else
// 		    {
// 		      sumpt += w[iphimax-iphi-1] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
// 		      Psi2pSin += w[iphimax-iphi-1] * ( sum*pt*sin(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
// 		      Psi2pCos += w[iphimax-iphi-1] * ( sum*pt*cos(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
// 		      //fprintf(stderr,"sumpt=%f \n", sumpt);
// 		    }
// 		}
// 	    }
// 	}
//       //      fprintf(stderr,"y=%f, sumpt=%f\n", y, sumpt/slope);
//     }
//   if (Psi2pCos<0.)
//     Psi2p = 1./2.*(atan(Psi2pSin/Psi2pCos)+PI);
//   else
//     Psi2p = 1./2.*(atan(Psi2pSin/Psi2pCos));
  
  double N;
  
  // next put in interpolation in phi and y to get same distances in phi and eta
  
  // outer loop deltaEta and deltaPhi (eta-eta2 and phi-phi2)
  // average over phi, eta
  // 1/(N*(N-1))*dN/dPhi dEta*dN/d(Phi-DeltaPhi)d(Eta-DeltaEta) = a
  // and 1/N*dN/dPhi dEta = b and 1/N*dN/d(phi-DeltaPhi)d(Eta-DeltaEta) = c
  // in the end average the three quantities over events
  // finally C=<a>-<b><c>

  double sum1, sum1Up, sum1Down, sum1PhiUp, sum1PhiDown;
  double sum2, sum2Up, sum2Down, sum2PhiUp, sum2PhiDown;
  double sum3, sum3Up, sum3Down, sum3PhiUp, sum3PhiDown;
  double sum4, sum4Up, sum4Down, sum4PhiUp, sum4PhiDown;
  double sumv21, sumv24;
  
  etamax=1.5;
  double trueetamax = 4.8;
  ietamax=14;
  double dEta = 2.*static_cast<double>(etamax)/static_cast<double>(ietamax);
  cout << "ietamax=" << ietamax << endl;
  int ideltaeta, ideltaphi, iy2, iphi1, iphi2, pointsForDeltaEta, pointsForDeltaPhi;
  //double deltaEta, deltaPhi;
  double suma, sumb, sumc, y2, phi1, phi2, sum1Average, sum4Average, sumv21Average, sumv24Average, sum1sum4Average;

  v2Delta=0.;
  v3Delta=0.;
  double norm=0.;

  for (ideltaeta=0; ideltaeta<29; ideltaeta++)
    {
      deltaEta = -2.8 + 0.2*ideltaeta; 
      for (ideltaphi=0; ideltaphi<46; ideltaphi++)
	{
	  pointsForDeltaEta=0;
	  pointsForDeltaPhi=0;
	  deltaPhi = -1.6 + 0.2*ideltaphi;
	  cout << "deltaeta=" << deltaEta << endl;
	  cout << "deltaphi=" << deltaPhi << endl;
	  sum1sum4Average=0.;
	  sumv21=0.;
	  sumv24=0.;
	  //for (ieta=ietamax/2; ieta<ietamax/2+1; ieta++)//0
	  for (ieta=0; ieta<ietamax; ieta++)// this is average eta (eta1+eta2)/2
	    {
	      eta = -etamax+dEta*ieta;
	      //cout << " eta-deltaEta=" << eta-deltaEta << endl;
	      if (eta-deltaEta/2.<-trueetamax)
		{
		  //  cout << "skipping for eta=" << eta << " and deltaEta=" << deltaEta << endl;
		  continue;
		}
	      if (eta+deltaEta/2.>trueetamax) 
		{
		  //cout << "skipping for eta=" << eta << " and deltaEta=" << deltaEta << endl;
		  continue;
		}
	      pointsForDeltaEta++;
	      for (ip = 0; ip<6; ip++)
		{
		  j = partid[MHALF+setOfNumbers[ip]];
		  m = particleList[j].mass;
		  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		  for (ipt=7; ipt<iptmax; ipt++) // put cuts right
		    {
		      pt = particleList[j].pt[ipt];
		      cout << "pt=" << pt << endl;
		      y = 0.5*log((sqrt(pt*pt*cosh(eta+deltaEta/2.)
					*cosh(eta+deltaEta/2.)+m*m)+pt*sinh(eta+deltaEta/2.))
				  /(sqrt(pt*pt*cosh(eta+deltaEta/2.)
					 *cosh(eta+deltaEta/2.)+m*m)-pt*sinh(eta+deltaEta/2.)));
		      
		      y2 = 0.5*log((sqrt(pt*pt*cosh(eta-deltaEta/2.)
					 *cosh(eta-deltaEta/2.)+m*m)+pt*sinh(eta-deltaEta/2.))
				  /(sqrt(pt*pt*cosh(eta-deltaEta/2.)
					 *cosh(eta-deltaEta/2.)+m*m)-pt*sinh(eta-deltaEta/2.)));
		      
		      // cout << "y=" << y << ", eta=" << eta << ",pt=" << pt << endl;
		      // cout << "y2=" << y2 << ", eta-deltaEta=" 
		      //    << eta-deltaEta << ", pt=" << pt << endl;
		      
		      
		      int ok=0;
		      for(i=0; i<iymax-1; i++) // for interpolation in y
			{
			  //			  cout << particleList[j].y[i] << endl;
			  if ( particleList[j].y[i]<= y &&  particleList[j].y[i+1] > y )
			    {
			      iy = i;
			      ok+=1;
			      //cout << "iy=" << iy << endl;
			    }
			}
		      if (ok==0)
			{
			  cout << "no iy found for y=" << y << endl; 
			  exit(1);
			}
		      ok=0;
		      for(i=0; i<iymax-1; i++) // for interpolation in y
			{
			  if ( particleList[j].y[i]<= y2 &&  particleList[j].y[i+1] > y2 )
			    {
			      iy2 = i;
			      ok+=1;
			      //cout << "iy2=" << iy2 << endl;
			    }
			}
		      if (ok==0) 
			{
			  cout << "no iy2 found for y2=" << y2 << endl; 
			  exit(1);
			}

		      mt = sqrt(pt*pt+m*m);
		      for (iphi=0; iphi<61; iphi++)
			//for (iphi=30; iphi<31; iphi++)
			{
			  phi = iphi*0.1047+deltaPhi/2.;
			  phi2 = iphi*0.1047-deltaPhi/2.;
			  
			  if (phi2<0) // go in a circle
			    phi2+=2.*M_PI;
			  if (phi>2*M_PI) // go in a circle
			    phi-=2.*M_PI;
			  if (phi<0) // go in a circle
			    phi+=2.*M_PI;
			  if (phi2>2*M_PI) // go in a circle
			    phi2-=2.*M_PI;
			  pointsForDeltaPhi++;
  
			  //if (phi<phiArray[0]) cout << "phi=" << phi << "<phiArray[0]=" << phiArray[0] << endl;
			  //if (phi2<phiArray[0]) cout << "phi2=" << phi2 << "<phiArray[0]=" << phiArray[0] << endl;

			  // cout << "phi=" << phi << endl;
			  //cout << "phi2=" << phi2 << endl;
                          //cout << "phimax=" << phiArray[iphimax-1] << endl;
		          //cout << "phimin=" << phiArray[0] << endl;
			  
			  ok=0;
			  for(i=0; i<iphimax-1; i++) // for interpolation in phi
			    {
			      if ( phiArray[i] <= phi 
				   &&  (phiArray[i+1] > phi || i+1 == iphimax-1))
				{
				  iphi1 = i;
				  //  cout << "phi=" << phi << endl;
				  //cout << "phiArray[i]=" << phiArray[i] << endl;
				  //cout << "phiArray[i+1]=" << phiArray[i+1] << endl;
				  ok+=1;
				}
			    }
			  if (ok==0 && (phi>phiArray[0]&&phi<phiArray[iphimax-1]))
			    {
			      cout << "no iphi found for phi=" << phi << endl; 
			      //cout << "phimax=" << phimax << endl; 
			      exit(1);
			    }
			  else if(phi<phiArray[0])
			    {
			      iphi1=iphimax-1;
			      phiArray[iphimax]=phiArray[0];
			      particleList[j].dNdydptdphi[iy][ipt][iphimax]=particleList[j].dNdydptdphi[iy][ipt][0];
			      particleList[j].dNdydptdphi[iy2][ipt][iphimax]=particleList[j].dNdydptdphi[iy2][ipt][0];
			      particleList[j].dNdydptdphi[iy+1][ipt][iphimax]=particleList[j].dNdydptdphi[iy][ipt][0];
			      particleList[j].dNdydptdphi[iy2+1][ipt][iphimax]=particleList[j].dNdydptdphi[iy2][ipt][0];
			    }
			  else if(phi>phiArray[iphimax-1])
			    {
			      iphi1=iphimax-1;
			      phiArray[iphimax]=phiArray[0];
			      particleList[j].dNdydptdphi[iy][ipt][iphimax]=particleList[j].dNdydptdphi[iy][ipt][0];
			      particleList[j].dNdydptdphi[iy2][ipt][iphimax]=particleList[j].dNdydptdphi[iy2][ipt][0];
			      particleList[j].dNdydptdphi[iy+1][ipt][iphimax]=particleList[j].dNdydptdphi[iy][ipt][0];
			      particleList[j].dNdydptdphi[iy2+1][ipt][iphimax]=particleList[j].dNdydptdphi[iy2][ipt][0];
			    }
			  ok=0;
			  for(i=0; i<iphimax-1; i++) // for interpolation in phi
			    {
			      if ( phiArray[i]<= phi2 
				   && (phiArray[i+1] > phi2 || i+1 == iphimax-1))
				{
				  iphi2 = i;
				  //cout << "phi2=" << phi2 << endl;
				  //cout << "phiArray[i]=" << phiArray[i] << endl;
				  //cout << "phiArray[i+1]=" << phiArray[i+1] << endl;
				  ok++;
				}
			    }
			  if (ok==0 && (phi2>phiArray[0]&&phi2<phiArray[iphimax-1]))
			    {
			      cout << "no iphi2 found for phi2=" << phi2 << endl; 
			      exit(1);
			    }
			  else if(phi2<phiArray[0])
			    {
			      iphi2=iphimax-1;
			      phiArray[iphimax]=phiArray[0];
			      particleList[j].dNdydptdphi2[iy][ipt][iphimax]=particleList[j].dNdydptdphi2[iy][ipt][0];
			      particleList[j].dNdydptdphi2[iy2][ipt][iphimax]=particleList[j].dNdydptdphi2[iy2][ipt][0];
			      particleList[j].dNdydptdphi2[iy+1][ipt][iphimax]=particleList[j].dNdydptdphi2[iy][ipt][0];
			      particleList[j].dNdydptdphi2[iy2+1][ipt][iphimax]=particleList[j].dNdydptdphi2[iy2][ipt][0];
			    }
			      else if(phi2>phiArray[iphimax-1])
			    {
			      iphi2=iphimax-1;
			      phiArray[iphimax]=phiArray[0];
			      particleList[j].dNdydptdphi2[iy][ipt][iphimax]=particleList[j].dNdydptdphi2[iy][ipt][0];
			      particleList[j].dNdydptdphi2[iy2][ipt][iphimax]=particleList[j].dNdydptdphi2[iy2][ipt][0];
			      particleList[j].dNdydptdphi2[iy+1][ipt][iphimax]=particleList[j].dNdydptdphi2[iy][ipt][0];
			      particleList[j].dNdydptdphi2[iy2+1][ipt][iphimax]=particleList[j].dNdydptdphi2[iy2][ipt][0];
			    }
		
			  //cout << "iphimax=" << iphimax << endl;
			  //cout << "iphi=" << iphi1 << endl;
			  //cout << "iphi2=" << iphi2 << endl;
			  //phiArray[iphi];
			  
			  //px = pt*cos(phi);
			  //py = pt*sin(phi);
	

			  // check interpolation for phi<phiArray[0]
			  
			  //if (phi>=phiArray[0])
			  //{
			      //sum1: this is at eta1, phi1
			      sum1Down = particleList[j].dNdydptdphi[iy][ipt][iphi1]
				*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
			      sum1Up = particleList[j].dNdydptdphi[iy+1][ipt][iphi1]
				*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
			      sum1PhiDown = sum1Down*(1.-(y-particleList[j].y[iy])/(particleList[j].y[iy+1]-particleList[j].y[iy]))
				+sum1Up*(y-particleList[j].y[iy])/(particleList[j].y[iy+1]-particleList[j].y[iy]);
			      
			      sum1Down = particleList[j].dNdydptdphi[iy][ipt][iphi1+1]
				*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
			      sum1Up = particleList[j].dNdydptdphi[iy+1][ipt][iphi1+1]
				*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
			      sum1PhiUp = sum1Down*(1-(y-particleList[j].y[iy])/(particleList[j].y[iy+1]-particleList[j].y[iy]))
				+sum1Up*(y-particleList[j].y[iy])/(particleList[j].y[iy+1]-particleList[j].y[iy]);
			      
			      sum1 = sum1PhiDown*(1.-(phi-phiArray[iphi1])/(phiArray[iphi1+1]-phiArray[iphi1]))
				+sum1PhiUp*(phi-phiArray[iphi1])/(phiArray[iphi1+1]-phiArray[iphi1]);
			      
// 			      sumv21=sum1*cos(2.*(phi-phimin-Psi2p));
			      
			      //cout << sum1 << endl;
			      //}
// 			  else
//  			    {
//  			      sum1Down = particleList[j].dNdydptdphi[iy][ipt][0]
//  				*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
// 			      sum1Up = particleList[j].dNdydptdphi[iy+1][ipt][0]
//  				*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
//  			      sum1 = sum1Down*(1-(y-particleList[j].y[iy])/(particleList[j].y[iy+1]-particleList[j].y[iy]))
//  				+sum1Up*(y-particleList[j].y[iy])/(particleList[j].y[iy+1]-particleList[j].y[iy]);
//  			    }

			  //Sum4: this is at eta1-deltaEta, phi1-deltaEta
			  //cout << particleList[j].dNdydptdphi[iy2][ipt][iphi2] << endl;
			  //cout << iy2 << " " << ipt << " " << iphi2 << " " << phi2 << endl;
			      //			  if(phi2>=phiArray[0])
			      //{
			      sum4Down = particleList[j].dNdydptdphi2[iy2][ipt][iphi2]
				*sqrt(1.-(m*m)/(mt*mt*cosh(y2)*cosh(y2)));
			      sum4Up = particleList[j].dNdydptdphi2[iy2+1][ipt][iphi2]
				*sqrt(1.-(m*m)/(mt*mt*cosh(y2)*cosh(y2)));
			      sum4PhiDown = sum4Down*(1-(y2-particleList[j].y[iy2])/(particleList[j].y[iy2+1]-particleList[j].y[iy2]))
				+sum4Up*(y2-particleList[j].y[iy2])/(particleList[j].y[iy2+1]-particleList[j].y[iy2]);
			      
			      sum4Down = particleList[j].dNdydptdphi2[iy2][ipt][iphi2+1]
				*sqrt(1.-(m*m)/(mt*mt*cosh(y2)*cosh(y2)));
			      sum4Up = particleList[j].dNdydptdphi2[iy2+1][ipt][iphi2+1]
				*sqrt(1.-(m*m)/(mt*mt*cosh(y2)*cosh(y2)));
			      sum4PhiUp = sum4Down*(1.-(y2-particleList[j].y[iy2])/(particleList[j].y[iy2+1]-particleList[j].y[iy2]))
				+sum4Up*(y2-particleList[j].y[iy2])/(particleList[j].y[iy2+1]-particleList[j].y[iy2]);
			      
			      sum4 = sum4PhiDown*(1.-(phi2-phiArray[iphi2])/(phiArray[iphi2+1]-phiArray[iphi2]))
				+sum4PhiUp*(phi2-phiArray[iphi2])/(phiArray[iphi2+1]-phiArray[iphi2]);
			      
// 			      sumv24=sum1*cos(2.*(phi-phimin-Psi2p));
			      
			      //cout << sum4 << endl;
			      //}
		// 	  else
//  			    {
//  			      sum4Down = particleList[j].dNdydptdphi[iy2][ipt][0]
//  				*sqrt(1.-(m*m)/(mt*mt*cosh(y2)*cosh(y2)));
//  			      sum4Up = particleList[j].dNdydptdphi[iy2+1][ipt][0]
//  				*sqrt(1.-(m*m)/(mt*mt*cosh(y2)*cosh(y2)));
//  			      sum4 = sum4Down*(1-(y2-particleList[j].y[iy2])/(particleList[j].y[iy2+1]-particleList[j].y[iy2]))
//  				+sum4Up*(y2-particleList[j].y[iy2])/(particleList[j].y[iy2+1]-particleList[j].y[iy2]);
//  			    }

//      cout << "sum=" << sum << ", sum1=" << sum1 << ", sumUp=" << sumUp << endl;
			  
			  // here integrate over pT and average over phi1 and eta1   
			  //0.1047 is step in phi
			  //dEta is step in eta
			  
			  // cout << "sum1=" << sum1 << endl;
			  //cout << "sum4=" << sum4 << endl;

			  if ( iphi == 0 || iphi == 60 )
			    {
			      if ( ieta == 0 || ieta == ietamax-1 )
				{
				  sum1sum4Average+=0.5*0.5*(0.1047 * ( sum1*sum4 ) * gala15w[ipt] * pt * dEta);
				  //sum1Average+=0.5*0.5*(0.1047 * ( sum1 ) * gala15w[ipt] * pt * dEta);
				  //sumv21Average+=0.5*0.5*(0.1047 * ( sumv21 ) * gala15w[ipt] * pt * dEta);
				  //sum4Average+=0.5*0.5*(0.1047 * ( sum4 ) * gala15w[ipt] * pt * dEta);
				  //sumv24Average+=0.5*0.5*(0.1047 * ( sumv24 ) * gala15w[ipt] * pt * dEta);
				}
			      else
				{
				  sum1sum4Average+=0.5*(0.1047 * ( sum1*sum4 ) * gala15w[ipt] * pt * dEta);
				  //sum1Average+=0.5*(0.1047 * ( sum1 ) * gala15w[ipt] * pt * dEta);
				  //sumv21Average+=0.5*(0.1047 * ( sumv21 ) * gala15w[ipt] * pt * dEta);
				  //sum4Average+=0.5*(0.1047 * ( sum4 ) * gala15w[ipt] * pt * dEta);
				  //sumv24Average+=0.5*(0.1047 * ( sumv24 ) * gala15w[ipt] * pt * dEta);
				}
			    }
			  else if ( ieta == 0 || ieta == ietamax-1 )
				 {
				   sum1sum4Average+=0.5*(0.1047 * ( sum1*sum4 ) * gala15w[ipt] * pt * dEta);
				   //sum1Average+=0.5*(0.1047 * ( sum1 ) * gala15w[ipt] * pt * dEta);
				   //sumv21Average+=0.5*(0.1047 * ( sumv21 ) * gala15w[ipt] * pt * dEta);
				   //sum4Average+=0.5*(0.1047 * ( sum4 ) * gala15w[ipt] * pt * dEta);
				   //sumv24Average+=0.5*(0.1047 * ( sumv24 ) * gala15w[ipt] * pt * dEta);
				 }
			  else
			    {
			      sum1sum4Average+=(0.1047 * ( sum1*sum4 ) * gala15w[ipt] * pt * dEta);
			      //sum1Average+=(0.1047 * ( sum1 ) * gala15w[ipt] * pt * dEta);
			      //sumv21Average+=(0.1047 * ( sumv21 ) * gala15w[ipt] * pt * dEta);
			      //sum4Average+=(0.1047 * ( sum4 ) * gala15w[ipt] * pt * dEta);
			      //sumv24Average+=(0.1047 * ( sumv24 ) * gala15w[ipt] * pt * dEta);
			    }
			}
		    }
		}
	    }

	  //compute background for ZYAM

	  // double B = (1.+2.*sumv21Average/sum1Average*sumv24Average/sum4Average*cos(2.*deltaPhi));

	  N=2; //fix this - has to be integrated over all deltaphi deltaeta...
	  
	  //sum1Average/=pointsForDeltaPhi; //average
	  //sum4Average/=pointsForDeltaPhi; //average
	  //sum1sum4Average/=(pointsForDeltaPhi); //average
	
	  //sum1Average/=pointsForDeltaEta; //average
	  //sum4Average/=pointsForDeltaEta; //average
	  //sum1sum4Average/=(pointsForDeltaEta); //average
	  
	  //sum1Average/=N;
	  //sum4Average/=N;
	  //sum1sum4Average/=(N*(N-1));

	  //cout << pointsForDeltaEta << endl;
	  //cout << (sum1sum4Average - sum1Average*sum4Average) << endl;
	  cout << sum1sum4Average << endl;
	  fprintf(p_file, "%f %f %e \n", deltaEta, deltaPhi, sum1sum4Average);

	  // integrate now and compute v2Delta and v3Delta
	  
	  //improve this integration:
	  if(deltaPhi>=0 && deltaPhi<=2.*M_PI)
	    if(deltaEta>=1.2 && deltaEta<=2.)
	      {
		v2Delta+=0.2*0.2*cos(2.*deltaPhi)*sum1sum4Average;
		v3Delta+=0.2*0.2*cos(3.*deltaPhi)*sum1sum4Average;
		norm+=0.2*0.2*sum1sum4Average;
	      }
	}
      fprintf(p_file,"\n");
    }
  
  cout<<"v_2Delta=" << v2Delta/norm << ", v_3Delta=" << v2Delta/norm << ", v3D/v2D=" << v3Delta/v2Delta << endl;
  
  fclose(p_file);
}



void Freeze::Compute3ChargedHadrons(InitData* DATA,double ptmax)
{

  char *numberStringy;
  char *numberString;
  char *numberStringPhi;
  char *numberStringv2;
  char *numberStringv2eta;
  char *numberStringv2tot;
  char *numberStringv3;
  char *numberStringv3r3;
  char *numberStringv3eta;
  char *numberStringv3tot;
  char *numberStringv3r3eta;
  char *numberStringv3r3tot;
  char *numberStringv4;
  char *numberStringv4eta;
  char *numberStringv4tot;
  char *numberStringeta;
  char buf[10];
  double *p, *w;		        // pointing to data for Gaussian integration in phi 
  double slope, slope1, slope2, fleft1, fleft2, fright1, fright2;
  int ipt, iphi, iymax, iy, iptmax, iphimax, ieta, ietamax;
  int i,j, iymaxIntegral;
  double pt, phi, px, py, y, deltaPT, deltaPhi, deltaY, ymax, phimin, phimax, sum, sumpt, sumpt2, sumv[5], m, sumv3r3;
  double phiOffs, phiDiff, ymaxIntegral, eta, deltaEta, etamax, mt, etaMaxIntegral;
  int returnValue;
  int number;
  int setOfNumbers[6] = {211,-211,2212,-2212,321,-321};
  int ip;
  

  // set some parameters

  number = setOfNumbers[0]; //use pion to get settings

  j = 1;
  fprintf(stderr,"Doing charged hadrons\n");
  phimax = particleList[j].phimax;
  phimin = particleList[j].phimin;
  
  phiOffs = 0.5 * ( phimin + phimax );
  phiDiff = 0.5 * ( phimax - phimin );
  iphimax = particleList[j].nphi;
  iptmax = particleList[j].npt;
  
  m = particleList[j].mass;

  if (iptmax != 15) 
    {
      fprintf(stderr,"only iptmax=15 possible. you picked %d. exiting.\n", iptmax);
      exit(1);
    }

  ymax = particleList[j].ymax;
  deltaY = particleList[j].deltaY;
  ymaxIntegral = 0.5-deltaY/2.;
  etaMaxIntegral = 1.3; // for v2 integration
  iymax = particleList[j].ny;
  iymaxIntegral = floor(2.*ymaxIntegral/deltaY);
  iptmax = particleList[j].npt;
  iphimax = particleList[j].nphi;


  // set up file names
  numberStringy = util->char_malloc(30);
  numberString = util->char_malloc(30);
  numberStringPhi = util->char_malloc(30);
  numberStringeta = util->char_malloc(30);
  numberStringv2 = util->char_malloc(30);
  numberStringv3 = util->char_malloc(30);
  numberStringv3r3 = util->char_malloc(30);
  numberStringv4 = util->char_malloc(30);
  numberStringv2eta = util->char_malloc(30);
  numberStringv3eta = util->char_malloc(30);
  numberStringv3r3eta = util->char_malloc(30);
  numberStringv4eta = util->char_malloc(30);
  numberStringv2tot = util->char_malloc(30);
  numberStringv3tot = util->char_malloc(30);
  numberStringv3r3tot = util->char_malloc(30);
  numberStringv4tot = util->char_malloc(30);
  
  if( chdir("./outputs") != 0 )
    {
      fprintf(stderr,"directory \"outputs\" does not exist. Exiting.\n");
      exit(1);
    }
  else returnValue=chdir("..");
  
  strcat(numberString, "./outputs/Npt-H+-");
  strcat(numberString,".dat");
  
  strcat(numberStringPhi, "./outputs/Nphi-H+-");
  strcat(numberStringPhi,".dat");
  
  strcat(numberStringy, "./outputs/Ny-H+-");
  strcat(numberStringy,".dat");
  
  strcat(numberStringv2, "./outputs/v2pt-H+-");
  strcat(numberStringv2,".dat");
  
  strcat(numberStringv3, "./outputs/v3pt-H+-");
  strcat(numberStringv3,".dat");
  
  strcat(numberStringv3r3, "./outputs/v3r3pt-H+-");
  strcat(numberStringv3r3,".dat");
  
  strcat(numberStringv4, "./outputs/v4pt-H+-");
  strcat(numberStringv4,".dat");
  
  strcat(numberStringeta, "./outputs/Neta-H+-");
  strcat(numberStringeta,".dat");
  
  strcat(numberStringv2eta, "./outputs/v2eta-H+-");
  strcat(numberStringv2eta,".dat");
  
  strcat(numberStringv3eta, "./outputs/v3eta-H+-");
  strcat(numberStringv3eta,".dat");
  
  strcat(numberStringv3r3eta, "./outputs/v3r3eta-H+-");
  strcat(numberStringv3r3eta,".dat");
  
  strcat(numberStringv4eta, "./outputs/v4eta-H+-");
  strcat(numberStringv4eta,".dat");
  
  strcat(numberStringv2tot, "./outputs/v2tot-H+-");
  strcat(numberStringv2tot,".dat");
  
  strcat(numberStringv3tot, "./outputs/v3tot-H+-");
  strcat(numberStringv3tot,".dat");
  
  strcat(numberStringv3r3tot, "./outputs/v3r3tot-H+-");
  strcat(numberStringv3r3tot,".dat");
  
  strcat(numberStringv4tot, "./outputs/v4tot-H+-");
  strcat(numberStringv4tot,".dat");
  
  FILE *y_file;
  char* y_name = numberStringy;
  y_file = fopen(y_name, "w");
  
  FILE *eta_file;
  char* eta_name = numberStringeta;
  eta_file = fopen(eta_name, "w");
  
  FILE *v2_eta_file;
  char* v2_eta_name = numberStringv2eta;
  v2_eta_file = fopen(v2_eta_name, "w");

  FILE *v3_eta_file;
  char* v3_eta_name = numberStringv3eta;
  v3_eta_file = fopen(v3_eta_name, "w");

  FILE *v3r3_eta_file;
  char* v3r3_eta_name = numberStringv3r3eta;
  v3r3_eta_file = fopen(v3r3_eta_name, "w");

  FILE *v4_eta_file;
  char* v4_eta_name = numberStringv4eta;
  v4_eta_file = fopen(v4_eta_name, "w");

  FILE *v2_tot_file;
  char* v2_tot_name = numberStringv2tot;
  v2_tot_file = fopen(v2_tot_name, "w");

  FILE *v3_tot_file;
  char* v3_tot_name = numberStringv3tot;
  v3_tot_file = fopen(v3_tot_name, "w");

  FILE *v3r3_tot_file;
  char* v3r3_tot_name = numberStringv3r3tot;
  v3r3_tot_file = fopen(v3r3_tot_name, "w");

  FILE *v4_tot_file;
  char* v4_tot_name = numberStringv4tot;
  v4_tot_file = fopen(v4_tot_name, "w");

  FILE *p_file;
  char* p_name = numberString;
  p_file = fopen(p_name, "w");

  FILE *phi_file;
  char* phi_name = numberStringPhi;
  phi_file = fopen(phi_name, "w");

  FILE *v2_file;
  char* v2_name = numberStringv2;
  v2_file = fopen(v2_name, "w");

  FILE *v3_file;
  char* v3_name = numberStringv3;
  v3_file = fopen(v3_name, "w");

  FILE *v3r3_file;
  char* v3r3_name = numberStringv3r3;
  v3r3_file = fopen(v3r3_name, "w");

  FILE *v4_file;
  char* v4_name = numberStringv4;
  v4_file = fopen(v4_name, "w");

  slope = particleList[j].slope;
  switch (iphimax) 
    {
    case 4: p= gaulep4; w= gaulew4; break;
    case 8: p= gaulep8; w= gaulew8; break;
    case 10: p= gaulep10; w= gaulew10; break;
    case 12: p= gaulep12; w= gaulew12; break;
    case 16: p= gaulep16; w= gaulew16; break;
    case 20: p= gaulep20; w= gaulew20; break;
    case 48: p= gaulep48; w= gaulew48; break;
    default: fprintf(stderr,"specified number of phi-points not available\n"); exit(1);
    }

  for(iphi=0; iphi<iphimax; iphi++)
    {
      if ( iphi < iphimax/2 )
	{
	  phiArray[iphi] = (phimax-phimin)/2.*(1.-p[iphi])+phimin;
	}
      else
	{
	  phiArray[iphi] = (phimax-phimin)/2.*(1.+p[iphimax-iphi-1])+phimin;
	}
      //      cout << "phi[" << iphi << "]=" << phiArray[iphi] << ", " <<  phimax/2.*(1.-p[iphi]) << ", " 
      //   << (phimax)/2.*(1.+p[iphimax-iphi-1]) << endl;
    }

  cout << "starting..." << endl;


  // retrieve value as function of phi, pt and y in sumPtPhi:
  double Psi3pSin = 0.;
  double Psi3pCos = 0.;
  double Psi3p;
  double Psi2pSin = 0.;
  double Psi2pCos = 0.;
  double Psi2p;
  for (iy=0; iy<iymax; iy++)
    {
      sumpt=0.;
      for (j=1; j<=3; j++)
	{
	  y =  particleList[j].y[iy];
	  //      fprintf(stderr,"y=%f \n", y);
	  for (ipt=0; ipt<iptmax; ipt++)
	    {
	      pt = particleList[j].pt[ipt];
	      //fprintf(stderr,"pt=%f \n", pt);
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  phi = phiArray[iphi];
		  sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
		  // integrate over pt and phi
		  //     sumpt += deltaPhi * ( sum ) * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		  // Psi3pSin += deltaPhi * ( sum*pt*sin(3.*phi) ) * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		  //Psi3pCos += deltaPhi * ( sum*pt*cos(3.*phi) ) * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		  if (iphi<iphimax/2)
		    {
		      sumpt += w[iphi] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi3pSin += w[iphi] * ( sum*pt*sin(3.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi3pCos += w[iphi] * ( sum*pt*cos(3.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi2pSin += w[iphi] * ( sum*pt*sin(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi2pCos += w[iphi] * ( sum*pt*cos(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      //fprintf(stderr,"sumpt=%f \n", sumpt);
		    }
		  else
		    {
		      sumpt += w[iphimax-iphi-1] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi3pSin += w[iphimax-iphi-1] * ( sum*pt*sin(3.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi3pCos += w[iphimax-iphi-1] * ( sum*pt*cos(3.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi2pSin += w[iphimax-iphi-1] * ( sum*pt*sin(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      Psi2pCos += w[iphimax-iphi-1] * ( sum*pt*cos(2.*phi) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      //fprintf(stderr,"sumpt=%f \n", sumpt);
		    }
		}
	    }
	}
      fprintf(y_file,"%f %f\n", y, 2.*sumpt/slope);
    }
  if (Psi2pCos<0.)
    Psi2p = 1./2.*(atan(Psi2pSin/Psi2pCos)+PI);
  else
    Psi2p = 1./2.*(atan(Psi2pSin/Psi2pCos));
  
  if (Psi3pCos<0.) 
    Psi3p = 1./3.*(atan(Psi3pSin/Psi3pCos)+PI);
  else
    Psi3p = 1./3.*(atan(Psi3pSin/Psi3pCos));
  
  cout << "Psi2p=" << Psi2p << endl; 
  cout << "Psi2pSin=" << Psi2pSin << endl; 
  cout << "Psi2pCos=" << Psi2pCos << endl; 
  cout << "Psi3p=" << Psi3p << endl; 
  cout << "Psi3pSin=" << Psi3pSin << endl; 
  cout << "Psi3pCos=" << Psi3pCos << endl; 

//   for (iy=0; iy<iymax; iy++)
//     {
//       y =  particleList[j].y[iy];
//       //      fprintf(stderr,"y=%f \n", y);
//       sumpt=0.;
//       for (ip = 1; ip<4; ip++)
// 	{
// 	  j = ip;
// 	  //	  cout << j << " " << particleList[j].name << endl;
// 	  for (ipt=0; ipt<iptmax; ipt++)
// 	    {
// 	      pt = particleList[j].pt[ipt];
// 	      //fprintf(stderr,"pt=%f \n", pt);
// 	      for (iphi=0; iphi<iphimax; iphi++)
// 		{
// 		  phi = phiArray[iphi];
// 		  sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
// 		  // integrate over pt and phi
// 		  if (iphi<iphimax/2)
// 		    {
// 		      sumpt += w[iphi] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
// 		      //fprintf(stderr,"sumpt=%f \n", sumpt);
// 		    }
// 		  else
// 		    {
// 		      sumpt += w[iphimax-iphi-1] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
// 		      //fprintf(stderr,"sumpt=%f \n", sumpt);
// 		    }
// 		}
// 	    }
// 	}
//       //      fprintf(stderr,"y=%f, sumpt=%f\n", y, sumpt/slope);
//       fprintf(y_file,"%f %f\n", y, 2.*sumpt/slope); //times two for anti-particles
//     }

  //pt spectra
  for (ipt=0; ipt<iptmax; ipt++)
    {
      sumpt = 0.;
      int countY = 0;
      pt = particleList[j].pt[ipt];
      for (ip = 1; ip<4; ip++)
	{
	  j = ip;
	  for (iy=0; iy<iymax; iy++)
	    {
	      y = iy*deltaY-ymax;
	      if (fabs(y)>ymaxIntegral) continue;
	      countY++;
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  phi = phiArray[iphi];
		  sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
		  if (iphi<iphimax/2)
		    {
		      if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
			{
			  sumpt += 0.5 * w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
			}
		      else
			{
			  sumpt += w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
			}
		      //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphi], iphi, p[iphi]);
		    }
		  else
		    {
		      if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
			{
			  sumpt += 0.5 * w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      else
			{
			  sumpt += w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphimax-iphi-1], iphi, -p[iphimax-iphi-1]);
		    }
		}
	    }
	}
      //fprintf(stderr,"pt=%f, sumpt=%f\n", pt, sumpt / (2.*ymaxIntegral));
      fprintf(p_file,"%e %e \n", pt, 2.*sumpt / (2.*ymaxIntegral));
      //fprintf(stderr,"pt=%f, sumv%d=%f\n", pt, 2, sumv[2]/sumpt);
    }


  // v2, v3 and v4
  for (ipt=0; ipt<iptmax; ipt++)
    {
      sumpt = 0.;
      sumv[2] = 0.;
      sumv[3] = 0.;
      sumv3r3 = 0.;
      sumv[4] = 0.;
      int countY = 0;
      for (j=1; j<=3; j++)
	{
	  pt = particleList[j].pt[ipt];
	  ymaxIntegral = 0.5*log((sqrt(pt*pt*cosh(etaMaxIntegral)*cosh(etaMaxIntegral)+m*m)
				  +pt*sinh(etaMaxIntegral))/(sqrt(pt*pt*cosh(etaMaxIntegral)*cosh(etaMaxIntegral)+m*m)-pt*sinh(etaMaxIntegral)));
	  iymaxIntegral = floor(2.*ymaxIntegral/deltaY);
	  for (iy=0; iy<iymax; iy++)
	    {
	      y = iy*deltaY-ymax;
	      if (fabs(y)>ymaxIntegral) continue;
	      countY++;
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  phi = phiArray[iphi];
		  sum = particleList[j].dNdydptdphi[iy][ipt][iphi];
		  
		  //       if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
// 		    {
// 		      sumpt += 0.5 * deltaPhi * ( sum ) * deltaY  * 2.*PI/(phimax-phimin);
// 		      sumv[2] += 0.5 * deltaPhi * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * 2.*PI/(phimax-phimin);
// 		      sumv[3] += 0.5 * deltaPhi * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * 2.*PI/(phimax-phimin);
// 		      sumv3r3 += 0.5 * deltaPhi * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * 2.*PI/(phimax-phimin);
// 		      sumv[4] += 0.5 * deltaPhi * ( sum*cos(4.*(phi-phimin)) ) * deltaY * 2.*PI/(phimax-phimin);
// 		    }
// 		  else
// 		    {
// 		      sumpt += deltaPhi * ( sum ) * deltaY * 2.*PI/(phimax-phimin);
// 		      sumv[2] += deltaPhi * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * 2.*PI/(phimax-phimin);
// 		      sumv[3] += deltaPhi * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * 2.*PI/(phimax-phimin);
// 		      sumv3r3 += deltaPhi * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * 2.*PI/(phimax-phimin);
// 		      sumv[4] += deltaPhi * ( sum*cos(4.*(phi-phimin)) ) * deltaY * 2.*PI/(phimax-phimin);
// 		    }
		
		  if (iphi<iphimax/2)
		    {
		      if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
			{
			  sumpt += 0.5 * w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
			  sumv[2] += 0.5 * w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[3] += 0.5 * w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  //sumv3r3 += 0.5 * w[iphi] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[4] += 0.5 * w[iphi] * ( sum*cos(4.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      else
			{
			  sumpt += w[iphi] * ( sum ) * deltaY * phiDiff  * 2.*PI/(phimax-phimin);
			  sumv[2] += w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[3] += w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  //sumv3r3 += w[iphi] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[4] += w[iphi] * ( sum*cos(4.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphi], iphi, p[iphi]);
		    }
		  else
		    {
		      if (countY==1 || countY==iymaxIntegral) // trapezoid rule: only half the edges
			{
			  sumpt += 0.5 * w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[2] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[3] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  //sumv3r3 += 0.5 * w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[4] += 0.5 * w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      else
			{
			  sumpt += w[iphimax-iphi-1] * ( sum ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[2] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[3] += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  //sumv3r3 += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			  sumv[4] += w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin)) ) * deltaY * phiDiff * 2.*PI/(phimax-phimin);
			}
		      //fprintf(stderr,"iphi=%d, w[%d]=%f, p[%d]=%f\n", iphi, iphi, w[iphimax-iphi-1], iphi, -p[iphimax-iphi-1]);
		    }
		}
	    }
	}
      //fprintf(stderr,"pt=%f, sumpt=%f\n", pt, sumpt / (2.*ymaxIntegral));
      //fprintf(stderr,"pt=%f, sumv%d=%f\n", pt, 2, sumv[2]/sumpt);
      fprintf(v2_file,"%e %e \n", pt, sumv[2]/sumpt);
      fprintf(v3_file,"%e %e \n", pt, sumv[3]/sumpt);
      fprintf(v3r3_file,"%e %e \n", pt, sumv3r3/sumpt);
      fprintf(v4_file,"%e %e \n", pt, sumv[4]/sumpt);
    }


  //pseudorapidity2
  deltaEta=deltaY;
  etamax=ymax;
  ietamax=floor(2.*etamax/deltaEta);
  for (ieta=0; ieta<ietamax; ieta++)
    {
      eta = ieta*deltaEta-etamax+deltaEta/2.;
      sumpt=0.;
      sumv[2]=0.;
      sumv[3]=0.;
      sumv3r3=0.;
      sumv[4]=0.;
      for (j=1; j<=3; j++)
	{
	  for (ipt=0; ipt<iptmax; ipt++)
	    {
	      pt = particleList[j].pt[ipt];
	      y = 0.5*log((sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)+pt*sinh(eta))/(sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)-pt*sinh(eta)));
	      
	      for(i=0; i<iymax-1; i++) // find closest iy
		{
		  if ( particleList[j].y[i]< y &&  particleList[j].y[i+1] > y )
		    {
		      if ( fabs(particleList[j].y[i]-y) <  fabs(particleList[j].y[i+1]-y)) 
			iy = i;
		      else
			iy = i+1;
		    }
		}
	      
	      mt = sqrt(pt*pt+m*m);
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  phi = phiArray[iphi];
		  px = pt*cos(phi);
		  py = pt*sin(phi);
		  sum = particleList[j].dNdydptdphi[iy][ipt][iphi]*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
		  if (iphi<iphimax/2)
		    {
		      sumpt += w[iphi] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[2] += w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[3] += w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		      //sumv3r3 += w[iphi] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[4] += w[iphi] * ( sum*cos(4.*(phi-phimin)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		    }
		  else
		    {
		      sumpt += w[iphimax-iphi-1] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[2] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[3] += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      //sumv3r3 += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[4] += w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		    }
		}
	    }
	}
      //   fprintf(stderr,"eta=%f, sumpt=%f\n", eta, sumpt/slope);
      fprintf(eta_file,"%f %f\n", eta, 2.*sumpt/slope);
      fprintf(v2_eta_file,"%f %f\n", eta, sumv[2]/sumpt);
      fprintf(v3_eta_file,"%f %f\n", eta, sumv[3]/sumpt);
      //fprintf(v3r3_eta_file,"%f %f\n", eta, sumv3r3/sumpt);
      fprintf(v4_eta_file,"%f %f\n", eta, sumv[4]/sumpt);
      
    }
  
  
  //total v2, v3 and v4
  deltaEta=deltaY;
  etamax=ymax;
  ietamax=floor(2.*etamax/deltaEta);
  double sumtot =0.;
  double sumtotv2 =0.;
  double sumtotv3 =0.;
  double sumtotv3r3 =0.;
  double sumtotv4 =0.;
  int eta1=0;
  for (ieta=0; ieta<ietamax; ieta++)
    {
      eta = ieta*deltaEta-etamax+deltaEta/2.;
      if (fabs(eta)>1.5) continue; // experimental cut
      eta1+=1;
      sumpt=0.;
      sumv[2]=0.;
      sumv[3]=0.;
      sumv3r3=0.;
      sumv[4]=0.;
      for (j=1; j<=3; j++)
	{
	  for (ipt=0; ipt<iptmax; ipt++)
	    {
	      pt = particleList[j].pt[ipt];
	      if (pt < 0.2 || pt > 4.) continue; // experimental cut 
	      y = 0.5*log((sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)+pt*sinh(eta))/(sqrt(pt*pt*cosh(eta)*cosh(eta)+m*m)-pt*sinh(eta)));
	      
	      for(i=0; i<iymax-1; i++) // find closest iy
		{
		  if ( particleList[j].y[i]< y &&  particleList[j].y[i+1] > y )
		    {
		      if ( fabs(particleList[j].y[i]-y) <  fabs(particleList[j].y[i+1]-y)) 
			iy = i;
		      else
			iy = i+1;
		    }
		}
	      
	      mt = sqrt(pt*pt+m*m);
	      for (iphi=0; iphi<iphimax; iphi++)
		{
		  phi = phiArray[iphi];
		  px = pt*cos(phi);
		  py = pt*sin(phi);
		  sum = particleList[j].dNdydptdphi[iy][ipt][iphi]*sqrt(1.-(m*m)/(mt*mt*cosh(y)*cosh(y)));
		  if (iphi<iphimax/2)
		    {
		      sumpt += w[iphi] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[2] += w[iphi] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[3] += w[iphi] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		      //sumv3r3 += w[iphi] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		      sumv[4] += w[iphi] * ( sum*cos(4.*(phi-phimin)) ) * gala15w[ipt] * pt * phiDiff * 2.*PI/(phimax-phimin);
		    }
		  else
		    {
		      sumpt += w[iphimax-iphi-1] * ( sum ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[2] += w[iphimax-iphi-1] * ( sum*cos(2.*(phi-phimin-Psi2p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[3] += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-Psi3p)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      //sumv3r3 += w[iphimax-iphi-1] * ( sum*cos(3.*(phi-phimin-DATA->Psi3r3)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		      sumv[4] += w[iphimax-iphi-1] * ( sum*cos(4.*(phi-phimin)) ) * phiDiff * gala15w[ipt] * pt * 2.*PI/(phimax-phimin);
		    }
		}
	    }
	}

      sumtot += sumpt * deltaEta;
      if (eta1 == 1)
	sumtot -= 0.5 * sumpt * deltaEta;
      sumtotv2 += sumv[2] * deltaEta;
      if (eta1 == 1)
	sumtotv2 -= 0.5 * sumv[2] * deltaEta;
      sumtotv3 += sumv[3] * deltaEta;
      if (eta1 == 1)
	sumtotv3 -= 0.5 * sumv[3] * deltaEta;
      sumtotv3r3 += sumv3r3 * deltaEta;
      if (eta1 == 1)
	sumtotv3r3 -= 0.5 * sumv3r3 * deltaEta;
      sumtotv4 += sumv[4] * deltaEta;
      if (eta1 == 1)
	sumtotv4 -= 0.5 * sumv[4] * deltaEta;
    }
  sumtot -= 0.5 * sumpt * deltaEta;
  sumtotv2 -= 0.5 * sumv[2] * deltaEta;
  sumtotv3 -= 0.5 * sumv[3] * deltaEta;
  sumtotv3r3 -= 0.5 * sumv3r3 * deltaEta;
  sumtotv4 -= 0.5 * sumv[4] * deltaEta;
  
  sumtotv2/=sumtot;
  sumtotv3/=sumtot;
  sumtotv3r3/=sumtot;
  sumtotv4/=sumtot;
  
  fprintf(v2_tot_file,"%f\n", sumtotv2);
  fprintf(v3_tot_file,"%f\n", sumtotv3);
  fprintf(v3r3_tot_file,"%f\n", sumtotv3r3);
  fprintf(v4_tot_file,"%f\n", sumtotv4);
  
  fclose(p_file);
  fclose(phi_file);
  fclose(v2_eta_file);
  fclose(v4_eta_file);
  fclose(v2_tot_file);
  fclose(v3_tot_file);
  fclose(v4_tot_file);
  fclose(eta_file);
  fclose(y_file);
  fclose(v2_file);
  fclose(v3_file);
  fclose(v3r3_file);
  fclose(v4_file);
}





// --------------------- resonance decays. adapted from azhydro ------------------------------------------------


/*************************************************
*
*	Edndp3
*
* 
**************************************************/
// This function interpolates the needed spectra for a given y, pt and phi.

double Freeze::Edndp3(double yr, double ptr, double phirin, int res_num)
/* 				/\* supersedes during test the right one *\/ */
/* 	double	yr;		/\* y  of resonance *\/ */
// if pseudofreeze flag is set, yr is the *pseudorapidity* of the resonance

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

  
  // If pseudofreeze flag is set,  dNdydptdphi is on a fixed grid in pseudorapidity. 
  // Set yr to the *pseudorapidity* of the resonance, and then interpolate the yield
  // at that value.
  if(pseudofreeze)
  {
    double yrtemp = yr;
    yr = PseudoRap(yrtemp, ptr, particleList[pn].mass);
  }
  
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
//       fprintf(stderr,"val2=%f\n",val2);
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


double Freeze::dnpir2N (double phi, void *para1)    
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

  dnr = Edndp3 (yR, ptR, phiR, para->res_num);

  /*printf(" phir = %15.8lf  ! ", phiR);
     printf(" ptR %15.8le jac %15.8le ", ptR, jac );
     printf(" dnr %15.8le \n", dnr); */

  return dnr * jac * jac / (2.0 * sume * sume);
}

double Freeze::norm3int (double x, void *paranorm) // this computes "Q(m_R,m_1,m_2,m_3)"
{
  nblock *tmp = (nblock *) paranorm;
  double res = sqrt ((tmp->a - x) * (tmp->b - x)
		     * (x - tmp->c) * (x - tmp->d)) / x;
  return res;
}

double Freeze::dnpir1N (double costh, void* para1)	       
{
  pblockN *para = (pblockN *) para1;
  double r;
  para->costh = costh;
  para->sinth = sqrt (1.0 - para->costh * para->costh);
  r = gauss (PTN2, &Freeze::dnpir2N, 0.0, 2.0 * PI, para); //Integrates the "dnpir2N" kernel over phi using gaussian integration
//   r = riemannsum (PTN2, &Freeze::dnpir2N, 0.0, 2.0 * PI, para); //Integrates the "dnpir2N" kernel over phi using trapezoid rule (same as riemann sum for periodic functions)
  return r;
}

double Freeze::dn2ptN (double w2, void* para1)
{
  pblockN *para = (pblockN *) para1;
  para->e0 = (para->mr * para->mr + para->m1 * para->m1 - w2) / (2 * para->mr); //particle one energy in resonance rest frame
  para->p0 = sqrt (para->e0 * para->e0 - para->m1 * para->m1); // particle one absolute value of three momentum on resonance rest frame
  return gauss (PTN1, &Freeze::dnpir1N, -1.0, 1.0, para); //Integrate the "dnpir1N" kernel over cos(theta) using gaussian integration
}

double Freeze::dn3ptN (double x, void* para1)  //The integration kernel for "W" in 3-body decays. x=invariant mass of other particles squared
{
  pblockN *para = (pblockN *) para1;
  double e0 =(para->mr * para->mr + para->m1 * para->m1 - x) / (2 * para->mr);
  double p0 = sqrt (e0 * e0 - para->m1 * para->m1);
  double a = (para->m2 + para->m3) * (para->m2 + para->m3);
  double b = (para->m2 - para->m3) * (para->m2 - para->m3);
  double re = p0 * sqrt ((x - a) * (x - b)) / x * dn2ptN (x, para);
  return re;
}

double Freeze::gauss(int n, double (Freeze::*f)(double, void *), double xlo, double xhi, void *optvec )
	{
	double	xoffs, xdiff; 
	int	ix;
	double	s;		/* summing up */
	double	*p, *w;		/* pointing to active list */

	switch ( n ) {
		case 4:		p= gaulep4; w= gaulew4; break;
		case 8:		p= gaulep8; w= gaulew8; break;
		case 10:	p=gaulep10; w=gaulew10; break;
		case 12:	p=gaulep12; w=gaulew12; break;
		case 16:	p=gaulep16; w=gaulew16; break;
		case 20:	p=gaulep20; w=gaulew20; break;
		case 48:	p=gaulep48; w=gaulew48; break;
		default:	printf("\ngauss():%d points not in list\n",n);
				exit(0);
		}
	xoffs = 0.5 * ( xlo + xhi );
	xdiff = 0.5 * ( xhi - xlo );
	s = 0;
	for( ix=0; ix<n/2; ix++ ) 	/* n is even */
	  s += w[ix] * ( (this->*f)(xoffs+xdiff*p[ix],optvec)
			     + (this->*f)(xoffs-xdiff*p[ix],optvec) );
	return( s * xdiff );
	}


// Left Riemann Sum.  For integration of periodic functions (e.g., over phi)
double Freeze::riemannsum(int n, double (Freeze::*f)(double, void *), double xlo, double xhi, void *optvec )
	{
	double	s=0;		/* summing up */
	double xdiff = xhi - xlo ;
	for(int ix=0; ix<n; ix++ )
	{
	  s +=  ( (this->*f)(ix*xdiff/n,optvec));
	}
	return( s * xdiff/n );
	}
	
	
/********************************************************************
*
*	Edndp3_2bodyN()
*
* transverse momentum spectrum in GeV^-2 from pions out of resonances
*********************************************************************/
double Freeze::Edndp3_2bodyN (double y, double pt, double phi, double m1, double m2, double mr, int res_num)
/* 		/\* in units of GeV^-2,includes phasespace and volume, */
/* 		   does not include degeneracy factors  *\/ */
/*      double y;			/\* rapidity of particle 1       *\/ */
/*      double pt;			/\* transverse momentum of particle 1    *\/ */
/*      double phi;		/\* phi angle of particle 1      *\/ */
/*      double m1, m2;		/\* restmasses of decay particles in MeV *\/ */
/*      double mr;			/\* restmass of resonance MeV            *\/ */
/*      int res_num;		/\* Montecarlo number of the Resonance   */ 

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
  res2 = norm2 * dn2ptN (m2 * m2, &para); //Calls the integration routines for 2-body
  if (res2<0.) res2=0.;
  return res2;			/* like Ed3ndp3_2body() */
}


/********************************************************************
*
*	Edndp3_3bodyN()
*
* transverse momentum spectrum in GeV^-2 from pions out of resonances
*********************************************************************/
double Freeze::Edndp3_3bodyN (double y, double pt, double phi, double m1, double m2,
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
  res3 = 2.0 * norm3 * gauss (PTS4, &Freeze::dn3ptN, wmin, wmax, &para) / mr;  //Integrates "W" using gaussian 
  if (res3<0.) res3=0.; 
  return res3;
}


/**************************************************************************
*									  *
*	add_reso							  *
*									  *
* computes the pt, mt distribution including resonance decays		  *
***************************************************************************/

void Freeze::add_reso (int pn, int pnR, int k, int j)
{
  nblock paranorm;		/* for 3body normalization integral */
  double y;
  double m1, m2, m3, mr;
  double norm3;			/* normalisation of 3-body integral */
  int pn2, pn3, pn4;		/* internal numbers for resonances */
  int part;
  int l, i, n;
  int ny, npt, nphi;

  ny = particleList[pn].ny; // for pseudofreeze, this variable stores number of points in pseudorapidity
  npt = particleList[pn].npt;
  nphi = particleList[pn].nphi;
  double deltaphi = 2*PI/nphi; // for pseudofreeze
  

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
		
	for (n = 0; n < ny; n++)
	  {
	    y = particleList[pn].y[n];
	    for (l = 0; l < npt; l++)
	      {
		if(pseudofreeze) y = Rap(particleList[pn].y[n],particleList[pn].pt[l],m1);
		for (i = 0; i < nphi; i++)
		  {
		    double phi;
		    if (pseudofreeze) phi = i*deltaphi;
		    else phi = phiArray[i];
		    double spectrum = Edndp3_2bodyN (y, particleList[pn].pt[l], phi, m1, m2, mr, particleList[pnR].number);
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
			if(n==ny/2 && i==0) {
			  // fprintf(stderr,"m1=%f, m2=%f, mr=%f, pnR=%d\n",m1,m2,mr,pnR);

			  particleList[pn].resCont[n][l][i]+= decay[j].branch *
			  spectrum; 
// 			  fprintf(stderr," %d %f %e %e %e %e\n", n, y, particleList[pn].pt[l], decay[j].branch *
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
	
	for (n = 0; n < ny; n++)
	  {
	    y = particleList[pn].y[n];
	    for (l = 0; l < npt; l++)
	      {
		if(pseudofreeze) y = Rap(particleList[pn].y[n],particleList[pn].pt[l],m1);
		for (i = 0; i < nphi; i++)
		  {
		    double phi;
		    if (pseudofreeze) phi = i*deltaphi;
		    else phi = phiArray[i];
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
		    
		    if(n==ny/2 && i==0)
		      {
			//	fprintf(stderr,"m1=%f, m2=%f, m3=%f, mr=%f, pnR=%d\n",m1,m2,m3,mr,pnR);
			particleList[pn].resCont[n][l][i]+= decay[j].branch *
			  spectrum;
		       
// 			fprintf(stderr,"%d %f %e %e %e %e\n",n,y, particleList[pn].pt[l], decay[j].branch *
// 				spectrum,
// 				particleList[pn].dNdydptdphi[n][l][i]-decay[j].branch *
// 				spectrum,particleList[pn].resCont[n][l][i]); 
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
	
	
	for (n = 0; n < ny; n++)
	  {
	    y = particleList[pn].y[n];
	    for (i = 0; i < nphi; i++)
	      {
		double phi;
		if (pseudofreeze) phi = i*deltaphi;
		else phi = phiArray[i];
		for (l = 0; l < npt; l++)
		  {
		      if(pseudofreeze) y = Rap(particleList[pn].y[n],particleList[pn].pt[l],m1);
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

void Freeze::cal_reso_decays (int maxpart, int maxdecay, int bound, int mode)
{
  // mode=4: do all
  // mode=5: do one
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
			add_reso (i, pnR, k, j);
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
			add_reso (i, pnaR, k, j);
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
				add_reso (i, pnR, k, j);
				add_reso (i, pnaR, k, j);
			      }
			  }
			else
			  {
			    if ((part == decay[j].part[k])
				&& (decay[j].numpart != 1))
			      {
				fprintf(stderr,"Partid is %i, %s into %s \n",pnR,particleList[pnR].name,particleList[i].name);
				//fprintf(stderr,"Calculating a decay \n");
				add_reso (i, pnR, k, j);
			      }
			    if ((-part == decay[j].part[k])
				&& (decay[j].numpart != 1))
			      {
				fprintf(stderr,"Partid is %i, %s into %s \n",pnaR,particleList[pnaR].name,particleList[i].name);
				//fprintf(stderr,"Calculating a decay \n");
				add_reso (i, pnaR, k, j);
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
			    add_reso (i, pnR, k, j);
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

void Freeze::CooperFrye(int particleSpectrumNumber, int mode, InitData *DATA, EOS *eos, int size, int rank)
{
  if(DATA->pseudofreeze) pseudofreeze = 1;
  else pseudofreeze =0;
  ReadParticleData(DATA, eos); // read in data for Cooper-Frye
  int i, b, number;
  if (mode == 3 || mode == 1) // compute thermal spectra
    {
      char *specString;
      specString = util->char_malloc(30);
      FILE *d_file;
      const char* d_name = "particleInformation.dat";
      d_file = fopen(d_name, "w");
//       fprintf(d_file,"");
      fclose(d_file);
      char buf[10];
      FILE *s_file;
      sprintf (buf, "%d", rank);
      strcat(specString, "yptphiSpectra");
      strcat(specString, buf);
      strcat(specString, ".dat");
      char* s_name = specString;
      s_file = fopen(s_name, "w");
      fclose(s_file);
      ReadFreezeOutSurface(DATA); // read freeze out surface (has to be done after the evolution of course)
      // (particle number, maximum p_T, [baryon=1, anti-baryon=-1], # of pts for pt integration, # of pts for phi integration)
      int ret;
      if (particleSpectrumNumber==0) // do all particles
	{
	  fprintf(stderr,"Doing all particles on this processor. May take a while ... \n");
	  for ( i=1; i<particleMax; i++ )
	    {
	      number = particleList[i].number;
	      b = particleList[i].baryon;
	      ComputeParticleSpectrum(DATA, number, 4., b, 15, 16, size, rank);
	    }
	  ReadSpectra(DATA);
	  for ( i=1; i<particleMax; i++ )
	    {
	      number = particleList[i].number;
	      b = particleList[i].baryon;
	      OutputFullParticleSpectrum(DATA, number, 4., b, 0);
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
	  ComputeParticleSpectrum(DATA, number, 4., b, 15, 16, size, rank);
	  // send something just to make sure that rank 0 waits for all others to be done:
	  int check[1];
	  check[0]=rank;
	  
	  if (rank > 0)
	    MPI::COMM_WORLD.Send(check,1,MPI::INT,0,1);
	    
	  if (rank == 0)
	    {
	      for (int from=1; from < size; from ++)
		MPI::COMM_WORLD.Recv(check,1,MPI::INT,from,1);

	      remove("yptphiSpectra.dat");
	      
	      switch (size) 
		{
		case 1: ret = system("cat yptphiSpectra0.dat > yptphiSpectra.dat"); break;
		case 2: ret = system("cat yptphiSpectra0.dat yptphiSpectra1.dat > yptphiSpectra.dat"); break;
		case 3: ret = system("cat yptphiSpectra0.dat yptphiSpectra1.dat yptphiSpectra2.dat > yptphiSpectra.dat"); break;
		case 4: ret = system("cat yptphiSpectra0.dat yptphiSpectra1.dat yptphiSpectra2.dat yptphiSpectra3.dat > yptphiSpectra.dat"); break;
		case 5: ret = system("cat yptphiSpectra0.dat yptphiSpectra1.dat yptphiSpectra2.dat yptphiSpectra3.dat yptphiSpectra4.dat > yptphiSpectra.dat"); break;
		case 6: ret = system("cat yptphiSpectra0.dat yptphiSpectra1.dat yptphiSpectra2.dat yptphiSpectra3.dat yptphiSpectra4.dat yptphiSpectra5.dat > yptphiSpectra.dat"); break;
		case 7: ret = system("cat yptphiSpectra0.dat yptphiSpectra1.dat yptphiSpectra2.dat yptphiSpectra3.dat yptphiSpectra4.dat yptphiSpectra5.dat yptphiSpectra6.dat> yptphiSpectra.dat"); break;
		case 8: ret = system("cat yptphiSpectra0.dat yptphiSpectra1.dat yptphiSpectra2.dat yptphiSpectra3.dat yptphiSpectra4.dat yptphiSpectra5.dat yptphiSpectra6.dat yptphiSpectra7.dat > yptphiSpectra.dat"); break;
		case 9: ret = system("cat yptphiSpectra0.dat yptphiSpectra1.dat yptphiSpectra2.dat yptphiSpectra3.dat yptphiSpectra4.dat yptphiSpectra5.dat yptphiSpectra6.dat yptphiSpectra7.dat yptphiSpectra8.dat > yptphiSpectra.dat"); break;
		case 10: ret = system("cat yptphiSpectra0.dat yptphiSpectra1.dat yptphiSpectra2.dat yptphiSpectra3.dat yptphiSpectra4.dat yptphiSpectra5.dat yptphiSpectra6.dat yptphiSpectra7.dat yptphiSpectra8.dat yptphiSpectra9.dat > yptphiSpectra.dat"); break;
		case 11: ret = system("cat yptphiSpectra0.dat yptphiSpectra1.dat yptphiSpectra2.dat yptphiSpectra3.dat yptphiSpectra4.dat yptphiSpectra5.dat yptphiSpectra6.dat yptphiSpectra7.dat yptphiSpectra8.dat yptphiSpectra9.dat yptphiSpectra10.dat > yptphiSpectra.dat"); break;
		case 12: ret = system("cat yptphiSpectra0.dat yptphiSpectra1.dat yptphiSpectra2.dat yptphiSpectra3.dat yptphiSpectra4.dat yptphiSpectra5.dat yptphiSpectra6.dat yptphiSpectra7.dat yptphiSpectra8.dat yptphiSpectra9.dat yptphiSpectra10.dat yptphiSpectra11.dat > yptphiSpectra.dat"); break;
		case 13: ret = system("cat yptphiSpectra0.dat yptphiSpectra1.dat yptphiSpectra2.dat yptphiSpectra3.dat yptphiSpectra4.dat yptphiSpectra5.dat yptphiSpectra6.dat yptphiSpectra7.dat yptphiSpectra8.dat yptphiSpectra9.dat yptphiSpectra10.dat yptphiSpectra11.dat yptphiSpectra12.dat > yptphiSpectra.dat"); break;
		case 14: ret = system("cat yptphiSpectra0.dat yptphiSpectra1.dat yptphiSpectra2.dat yptphiSpectra3.dat yptphiSpectra4.dat yptphiSpectra5.dat yptphiSpectra6.dat yptphiSpectra7.dat yptphiSpectra8.dat yptphiSpectra9.dat yptphiSpectra10.dat yptphiSpectra11.dat yptphiSpectra12.dat yptphiSpectra13.dat > yptphiSpectra.dat"); break;
		case 15: ret = system("cat yptphiSpectra0.dat yptphiSpectra1.dat yptphiSpectra2.dat yptphiSpectra3.dat yptphiSpectra4.dat yptphiSpectra5.dat yptphiSpectra6.dat yptphiSpectra7.dat yptphiSpectra8.dat yptphiSpectra9.dat yptphiSpectra10.dat yptphiSpectra11.dat yptphiSpectra12.dat yptphiSpectra13.dat yptphiSpectra14.dat > yptphiSpectra.dat"); break;
		case 16: ret = system("cat yptphiSpectra0.dat yptphiSpectra1.dat yptphiSpectra2.dat yptphiSpectra3.dat yptphiSpectra4.dat yptphiSpectra5.dat yptphiSpectra6.dat yptphiSpectra7.dat yptphiSpectra8.dat yptphiSpectra9.dat yptphiSpectra10.dat yptphiSpectra11.dat yptphiSpectra12.dat yptphiSpectra13.dat yptphiSpectra14.dat yptphiSpectra15.dat > yptphiSpectra.dat"); break;
		case 32: ret = system("cat yptphiSpectra0.dat yptphiSpectra1.dat yptphiSpectra2.dat yptphiSpectra3.dat yptphiSpectra4.dat yptphiSpectra5.dat yptphiSpectra6.dat yptphiSpectra7.dat yptphiSpectra8.dat yptphiSpectra9.dat yptphiSpectra10.dat yptphiSpectra11.dat yptphiSpectra12.dat yptphiSpectra13.dat yptphiSpectra14.dat yptphiSpectra15.dat yptphiSpectra16.dat yptphiSpectra17.dat yptphiSpectra18.dat yptphiSpectra19.dat yptphiSpectra20.dat yptphiSpectra21.dat yptphiSpectra22.dat yptphiSpectra23.dat yptphiSpectra24.dat yptphiSpectra25.dat yptphiSpectra26.dat yptphiSpectra27.dat yptphiSpectra28.dat yptphiSpectra29.dat yptphiSpectra30.dat yptphiSpectra31.dat > yptphiSpectra.dat"); break;
		default: fprintf(stderr,"maximum number of processors is 32 at the moment - intermediate ones are missing!\n");exit(1);
		}
	      
	      ReadSingleSpectrum(DATA);
	      cout << "output results..." << endl;
	      OutputFullParticleSpectrum(DATA, number, 4., b, 0);
	      //Read3Spectra(DATA);
	      //Compute3ChargedHadrons(DATA, 4.);
	    }
	}
    }
  else if (mode==4 || mode==5) //  do resonance decays
    {
      ReadSpectra(DATA);
      FILE *d_file;
      const char* d_name = "FparticleInformation.dat";
      d_file = fopen(d_name, "w");
//       fprintf(d_file,"");
      fclose(d_file);
      FILE *s_file;
      const char* s_name = "FyptphiSpectra.dat";
      s_file = fopen(s_name, "w");
//       fprintf(s_file,"");
      fclose(s_file);
      int bound = 211; //number of lightest particle to calculate. 
      if (mode==4) // do resonance decays
	{
	  fprintf(stderr,"doing all from %i: %s to %i: %s.\n",particleMax,particleList[particleMax].name,
		  partid[MHALF+bound],particleList[partid[MHALF+bound]].name);
	  cal_reso_decays(particleMax,decayMax,bound,mode);
	  for ( i=1; i<particleMax; i++ )
	    {
	      number = particleList[i].number;
	      b = particleList[i].baryon;
	      OutputFullParticleSpectrum(DATA, number, 4., b, 1);
	      ComputeAveragePT(number,4.);
	    }
	  if(particleMax>=20)
	    {
	      ReadFullSpectra(DATA);
	      ComputeChargedHadrons(DATA,4.);
	    }
	}
      else if (mode==5) // only for testing - this will miss the complete chain of decays
	{
	  cal_reso_decays(particleSpectrumNumber+1,decayMax,bound,mode);
          number = particleList[particleSpectrumNumber].number;
	  b = particleList[particleSpectrumNumber].baryon;
	  OutputFullParticleSpectrum(DATA, number, 4., b, 1);
	}
    }
  else if (mode==6) //  do additional manipulation
    {
      ReadFullSpectra(DATA);
      ComputeChargedHadrons(DATA,4.);
      number = particleList[particleSpectrumNumber].number;
      b = particleList[particleSpectrumNumber].baryon;
      OutputFullParticleSpectrum(DATA, number, 4., b, 1);
      for ( i=1; i<particleMax; i++ )
	{
	  number = particleList[i].number;
	  ComputeAveragePT(number,4.);
	}
      //  util->vector_free(phiArray);
      //ComputeChargedHadrons(DATA,4.);
    }
  else if (mode==7) //  do additional manipulation
    {
      Read3Spectra(DATA);
      Compute3ChargedHadrons(DATA, 4.);
      util->vector_free(phiArray);
      //ComputeChargedHadrons(DATA,4.);
    }
  else if (mode==8) //  compute correlations
    {
      ReadFullSpectra(DATA);
      ReadFullSpectra2(DATA);
      ComputeCorrelations(DATA, 4.);
    }
  else if (mode==9) // output full spectra
    {
      ReadFullSpectra(DATA);
      for ( i=1; i<particleMax; i++ )
	{
	  number = particleList[i].number;
	  b = particleList[i].baryon;
	  OutputFullParticleSpectrum(DATA, number, 4., b, 1);
	}
      ComputeChargedHadrons(DATA,4.);
    }
}
