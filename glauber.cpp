#include "glauber.h"

using namespace std;

Glauber::Glauber()
{
  util = new Util;
}

// destructor
Glauber::~Glauber()
{
  delete util;
  remove("tmp.dat");
}


void Glauber::FindNucleusData(Nucleus *nucleus, string target, string file_name, int rank)
{
  string func_name;
 char *tmp_name; 
//  char *buf;
//  double x;
 static int ind;
//  int bytes_read;
 string inputname;
 
 stringstream tmpfilename;
 inputname = file_name;
 
 // s = util->char_malloc(120);
 string s;
 // name = util->char_malloc(120);
 string name;
 ifstream input (inputname.c_str());

 char c;
 tmpfilename << "tmp";
 tmpfilename << rank;
 tmpfilename << ".dat"; 
 
 string fn;
 fn = tmpfilename.str();
 
 ofstream tmp_file(fn.c_str());
 // bytes_read=fscanf(input, "%s", s);

//  target;

 input >> s;
 ind = 0;
 while(s.compare("EndPfData") != 0)
  {
    input >> name;
    if(name.compare(target) == 0)
    {
      ind++;
      tmp_file << "Name " << name << endl;
      c = input.get();
     while(c != 'N')
      {
	tmp_file << c;
	c = input.get();
      }/* while */
     break;
    }/* if target is found */
    input >> s;
  }/* while end of the file is not encountered */
 tmp_file << "EndOfData" << endl;
 
 tmp_file.close();
 input.close();

 nucleus->rho_WS = 0.15; /* default rho.  this WILL change to the right 
			 value that gives integral(rho) = A */
 nucleus->name = util->StringFind(fn, "Name");
 nucleus->A = util->DFind(fn, "A");
 nucleus->Z = util->DFind(fn, "Z");
 
 if(ind != 0)
  {
    nucleus->w_WS = util->DFind(fn, "w_WS");
    nucleus->a_WS = util->DFind(fn, "a_WS");
    nucleus->R_WS = util->DFind(fn, "R_WS");
    func_name = util->StringFind(fn, "density_func");
  }
 else
  {
   nucleus->w_WS = 0.0; 
   nucleus->a_WS = 0.53;
   nucleus->R_WS = 1.15*pow(nucleus->A, 1.0/3.0);
   func_name = "3Fermi";
  }

 string funcname;
 stringstream strfuncname;
 strfuncname << func_name;
 strfuncname >> funcname;
 
 if(funcname.compare("2HO")==0) 
  {
    nucleus->AnumFunc = 1; //Anum2HO;
    nucleus->AnumFuncIntegrand = 1; //Anum2HOInt;
    nucleus->DensityFunc = 1; //NuInt2HO;
  }
 else if(funcname.compare("3Gauss")==0) 
  {
    nucleus->AnumFunc = 2; //Anum3Gauss;
    nucleus->AnumFuncIntegrand = 2; //Anum3GaussInt;
    nucleus->DensityFunc = 2; //NuInt3Gauss;
  }
 else if(funcname.compare("3Fermi")==0) 
  {
    nucleus->AnumFunc = 3; //Anum3Fermi;
    nucleus->AnumFuncIntegrand = 3; //Anum3FermiInt;
    nucleus->DensityFunc = 3; //NuInt3Fermi;
  }


 remove(tmp_name);
}/* FindNucleusData */



void Glauber::PrintLexusData()
{
 fprintf(stderr, "LexusData.SigmaNN = %e\n",  LexusData.SigmaNN);
 fprintf(stderr, "LexusData.InterMax = %d\n", LexusData.InterMax);
 fprintf(stderr, "LexusData.SCutOff = %f\n", LexusData.SCutOff);
}/* PrintLexusData */


void Glauber::PrintNucleusData(Nucleus *nucleus)
{
  cout << "Nucleus Name: " << nucleus->name << endl;
  cout << " Nucleus.A = " << nucleus->A << endl;
  cout << " Nucleus.Z = " << nucleus->Z << endl;
  cout << " Nucleus.w_WS = " << nucleus->w_WS << endl;
  cout << " Nucleus.a_WS = " << nucleus->a_WS << endl;
  cout << " Nucleus.R_WS = " << nucleus->R_WS << endl;

}/* FindNucleusData */


int Glauber::LinearFindXorg(double x, double *Vx, int ymax)
{
/* finds the first of the 4 points, x is between the second and the third */
 
 int x_org;
 double nx;

 nx = ymax*(x - Vx[0])/(Vx[ymax] - Vx[0]);
 
 x_org = (int) nx;
 x_org -= 1;

 if( x_org <= 0 ) return 0;
 else if(x_org >= ymax - 3) return ymax - 3;
 else return x_org;

}/* Linear Find Xorg */


double Glauber::FourPtInterpolate(double x, double *Vx, double *Vy, double h, int x_org, int ymax)
{
 /* interpolating points are x_org, x_org+1, x_org+2, x_org+3 */
 /* cubic polynomial approximation */

 double a, b, c, d, f;

 MakeCoeff(&a, &b, &c, &d,  Vy, Vx, h, x_org);
 
 f = a*pow(x - Vx[x_org], 3.);
 f += b*pow(x - Vx[x_org], 2.);
 f += c*(x - Vx[x_org]);
 f += d;
 
 return f;
}/* FourPtInterpolate */


void Glauber::MakeCoeff(double *a, double *b, double *c, double *d, 
			double *Vy, double *Vx, double h, int x_org)
{
 double f0, f1, f2, f3;
//  double x1, x2, x3;
//  double x1sqr, x2sqr, x3sqr;
//  double x1cube, x2cube, x3cube;

 f0 = Vy[x_org];
 f1 = Vy[x_org+1];
 f2 = Vy[x_org+2];
 f3 = Vy[x_org+3];

 *a =  (-f0 + 3.0*f1 - 3.0*f2 + f3)/(6.0*h*h*h);
 
 *b =  (2.0*f0 - 5.0*f1 + 4.0*f2 - f3)/(2.0*h*h);

 *c =  (-11.0*f0 + 18.0*f1 - 9.0*f2 + 2.0*f3)/(6.0*h);

 *d = f0;

}/* MakeCoeff */

int Glauber::FindXorg(double x, double *Vx, int ymax)
{
 int i,  x_org;
//  int ymid;

 i = 0;
 while(Vx[i] < x) i++;
 
 x_org = i - 2;
 
 if( x_org <= 1 ) return 1;
 else if(x_org >= ymax - 3) return ymax - 3;
 else return x_org;

}/* Find Xorg */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

double Glauber::VInterpolate(double x, double *Vx, double *Vy, int ymax)
{
 int x_org;
 double h;

 if( (x < Vx[0])||(x > Vx[ymax]) )
  {
   fprintf(stderr, 
           "VInterpolate: x = %le is outside the range (%le, %le).\n", 
	    x,Vx[0],Vx[ymax]);
   fprintf(stderr, "This can't happen.  Exiting...\n");
   exit(0);
  }

/* we only deal with evenly spaced Vx */
/* x_org is the first of the 4 points */

 x_org = LinearFindXorg(x, Vx, ymax);

 h = (Vx[ymax] - Vx[0])/ymax;

 return FourPtInterpolate(x, Vx, Vy, h, x_org, ymax);

}/* VInterpolate */


double *Glauber::MakeVx(double down, double up, int maxi_num)
{
 static double dx, *vx;
 int i;

 vx = util->vector_malloc(maxi_num + 1);
 dx = (up - down)/maxi_num;
 
 for(i=0; i<=maxi_num; i++)
  {
   vx[i] = dx*i;
  }

 return vx;

}/* MakeVx */


double *Glauber::MakeVy(string st, double *vx, int maxi_num)
{
 int i, di;
 static double *vy;
//  static char *dst;

 if(maxi_num > 200) di = 100;
 if(maxi_num <= 200) di = 20;
 
 vy = util->vector_malloc(maxi_num + 1);

 ofstream data_file(st.c_str());
 
 data_file << "EndOfData" << endl;
 
 for(i=0; i<=maxi_num; i++)
  {
   vy[i] = NuInS(vx[i]);
   if(i % di == 0)
    {
      cerr << st << "[" << i << "] = " << vy[i] << endl; 
    }
   data_file << vx[i] << " " << vy[i] << endl;
  }

 data_file.close();

 return vy;
}/* MakeVy */


double *Glauber::ReadInVx(char *file_name, int maxi_num, int quiet)
{
 static double x, *vx;
 int i;
 FILE *input;
 static char *s, *sx;
 int bytes_read;
 s = util->char_malloc(120);
 sx = util->char_malloc(120);
 
 vx = util->vector_malloc(maxi_num + 1);

 if(quiet == 1)
  {
   fprintf(stderr, "Reading in Vx from %s ...\n", file_name);
   }
 
 input = fopen(file_name, "r");
 bytes_read=fscanf(input, "%s", s);
 while(strcmp(s, "EndOfData") != 0)
  {
   bytes_read=fscanf(input, "%s", sx);
   bytes_read=fscanf(input, "%s", s);
  }

 for(i=0; i<=maxi_num; i++)
  {
   bytes_read=fscanf(input, "%lf", &x);
   vx[i] = x;
   bytes_read=fscanf(input, "%lf", &x);
  }
 fclose(input);

 util->char_free(sx);
 util->char_free(s);
 return vx;

}/* ReadInVx */


double *Glauber::ReadInVy(char *file_name, int maxi_num, int quiet)
{
 static double y, *vy; 
//  static double x;
 int i;
 FILE *input;
 static char *s, *sy;
 int bytes_read;
 s = util->char_malloc(120);
 sy = util->char_malloc(120);
 
 vy = util->vector_malloc(maxi_num + 1);

 if(quiet == 1)
 {
  fprintf(stderr, "Reading in Vy from %s ...\n", file_name);
 }
 
 input = fopen(file_name, "r");
 bytes_read=fscanf(input, "%s", s);
 while(strcmp(s, "EndOfData") != 0)
  {
   bytes_read=fscanf(input, "%s", sy);
   bytes_read=fscanf(input, "%s", s);
  }

 for(i=0; i<=maxi_num; i++)
  {
   bytes_read=fscanf(input, "%lf", &y);
   bytes_read=fscanf(input, "%lf", &y);
   vy[i] = y;
  }
 fclose(input);
 
 util->char_free(s);
 util->char_free(sy);
 return vy;

}/* ReadInVy */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%% */

double Glauber::InterNuPInSP(double s)
{
 double y;
 static int ind = 0;
 static double up, down; 
 static int maxi_num; 
 static double *vx, *vy;
//  int i;
//  double dx;
//  FILE *output;
 ind++;

 string st;

 if(LexusData.Projectile.A == 1.0) return 0.0;

 CalcRho(&(LexusData.Projectile));

 up = 2.0*LexusData.SCutOff;
 down = 0.0; 
 maxi_num = LexusData.InterMax; 

 st = "./NuPInSP.dat";

 if(ind == 1)
  {
    vx = MakeVx(down, up, maxi_num);
    vy = MakeVy(st, vx, maxi_num);
  }/* if ind */
 
 if(s > up) return 0.0;
 else{
   y = VInterpolate(s, vx, vy, maxi_num);
   if( y < 0.0 ) return 0.0; 
   else 
     {
       return y;
     }
  }
}/* InterNuPInSP */

double Glauber::InterNuTInST(double s)
{
 double y;
 static int ind = 0;
 static double up, down; 
 static int maxi_num; 
 static double *vx, *vy;
//  int i;
//  double dx;
 //static char *st, *dst;
//  FILE *output;

 ind++;
 if(LexusData.Target.A == 1.0) return 0.0;

 
 const char *paf = "./NuTInST.dat";
//  paf = "./NuTInST";
// strcpy(paf, "./NuTInST");

 //st = util->char_malloc(120);
//  dst = ".dat";
 //strcpy(dst, "./NuTInST");

 //st = strcpy(st, paf);
 //st = strcat(st, dst);
 
 if(ind == 1) 
  {
   CalcRho(&(LexusData.Target));
   
   up = 2.0*LexusData.SCutOff;
   down = 0.0; 
   maxi_num = LexusData.InterMax; 
     
//    if(util->IsFile(st))
//     {
//      vx = ReadInVx(st, maxi_num, 1);
//      vy = ReadInVy(st, maxi_num, 1);
//     }
//    else
//     {
   // output = fopen(st, "w");
   // fprintf(stderr, "Name %s\n", LexusData.Target.name);
   //fprintf(stderr, "SigmaNN %e\n", LexusData.SigmaNN);
   //fclose(output);
  
     vx = MakeVx(down, up, maxi_num);
     vy = MakeVy(paf, vx, maxi_num);
//    }
  }/* if ind */
 
 //util->char_free(st);

 if(s > up) return 0.0;
 else{
   y = VInterpolate(s, vx, vy, maxi_num);
   if( y < 0.0 ) return 0.0; 
   else return y;
  }
}/* InterNuTInST */

void Glauber::CalcRho(Nucleus *nucleus)
{
 double f, R_WS;
//  double rho;
 double down;
 int count;
//  static int ind = 0;
/* to pass to AnumIntegrand */ 
 
 Nuc_WS = nucleus;
 
 down = 0.0;
 count = 0;

 R_WS = nucleus->R_WS;   
 
 if (nucleus->AnumFunc==1)
   f = Anum2HO(R_WS)/(nucleus->rho_WS);
 else if (nucleus->AnumFunc==2)
   f = Anum3Gauss(R_WS)/(nucleus->rho_WS);
 else if (nucleus->AnumFunc==3)
   f = Anum3Fermi(R_WS)/(nucleus->rho_WS);
 nucleus->rho_WS = (nucleus->A)/f;
}/* CalcRho */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

double Glauber::NuInS(double s)
{
 double y;
 int count;
//  static int ind = 0;
 int id;

/* to pass to the DensityFunc's */
 NuInS_S = s;

 id = Nuc_WS->DensityFunc;

 count = 0;
 y = integral(id, 0.0, 1.0, TOL, &count); 
 
 return y;
}

double Glauber::Anum3Fermi(double R_WS)
{
 int count=0;
 double up, down, a_WS, w_WS, rho, f;

 a_WS = Nuc_WS->a_WS;
 w_WS = Nuc_WS->w_WS;
 rho = Nuc_WS->rho_WS;

/* to pass to Anumintegrand */
 AnumR = R_WS/a_WS;
 
 down = 0.0;
 up = 1.0;
 
 f = integral(4, down, up, TOL, &count);
 f *= 4.0*M_PI*rho*pow(a_WS, 3.);

 return f;
}/* Anum3Fermi */


double Glauber::Anum3FermiInt(double xi) 
{
 double f;
 double r;
 double R_WS, w_WS;

 if(xi==0.0) xi = tiny;
 if(xi==1.0) xi = 1.0-tiny;

 w_WS = Nuc_WS->w_WS;

/* already divided by a */
 R_WS = AnumR;
 r = -log(xi);

 f = r*r;
 f *= ( 1.0 + w_WS*pow(r/R_WS, 2.) );
 f /= (xi + exp(-R_WS));
  
 return f;
}/* Anum3FermiInt */


double Glauber::NuInt3Fermi(double xi)
{
 double f;
//  double argexp;
 double c;
 double z, r, s;
 double w_WS, R_WS, a_WS, rho;
//  static int ind = 0;

 a_WS = Nuc_WS->a_WS;
 w_WS = Nuc_WS->w_WS;
 R_WS = Nuc_WS->R_WS;
 rho = Nuc_WS->rho_WS;

/* xi = exp(-z/a), r = sqrt(z^2 + s^2) */

   if ( xi == 0.0 ) xi = tiny; 
   if ( xi == 1.0 ) xi = 1.0-tiny; 
   
/* devide by a_WS, make life simpler */
   
   s = NuInS_S/a_WS;
   z = -log(xi);
   r = sqrt(s*s + z*z);
   R_WS /= a_WS;

   c = exp(-R_WS);
   
   f = 2.0*a_WS*rho*(LexusData.SigmaNN);
   f *= 1.0 + w_WS*pow(r/R_WS, 2.);
   f /= xi + c*exp(s*s/(r + z)); 
 
 return f;
}/* NuInt3Fermi */


/* %%%%%%% 3 Parameter Gauss %%%%%%%%%%%% */


double Glauber::Anum3Gauss(double R_WS)
{
 int count=0;
 double up, down, a_WS, w_WS, rho, f;

 a_WS = Nuc_WS->a_WS;
 w_WS = Nuc_WS->w_WS;
 rho = Nuc_WS->rho_WS;

/* to pass to Anumintegrand */
 AnumR = R_WS/a_WS;
 
 down = 0.0;
 up = 1.0;
 
 f = integral(5, down, up, TOL, &count);
 f *= 4.0*M_PI*rho*pow(a_WS, 3.);

 return f;
}/* Anum3Gauss */


double Glauber::Anum3GaussInt(double xi) 
{
 double y;
 double r_sqr;
 double R_WS, w_WS;

 if(xi==0.0) xi = tiny;
 if(xi==1.0) xi = 1.0 - tiny;

 w_WS = Nuc_WS->w_WS;

/* already divided by a */
 
 R_WS = AnumR;
 
 r_sqr = -log(xi);

 y = sqrt(r_sqr);
 
 y *= 1.0 + w_WS*r_sqr/pow(R_WS, 2.);

/* 2 comes from dr^2 = 2 rdr */

 y /= 2.0*(xi + exp(-R_WS*R_WS));
 
 return y;
}/* Anum3GaussInt */


double Glauber::NuInt3Gauss(double xi)
{
 double f;
//  double argexp;
 double c;
 double z_sqr, r_sqr, s;
 double w_WS, R_WS, a_WS, rho;
//  static int ind = 0;

 a_WS = Nuc_WS->a_WS;
 w_WS = Nuc_WS->w_WS;
 R_WS = Nuc_WS->R_WS;
 rho = Nuc_WS->rho_WS;

/* xi = exp(-z*z/a/a), r = sqrt(z^2 + s^2) */

   if ( xi == 0.0 ) xi = tiny; 
   if ( xi == 1.0 ) xi = 1.0-tiny; 

/* devide by a_WS, make life simpler */

   s = NuInS_S/a_WS;
   z_sqr = -log(xi);
   r_sqr = s*s + z_sqr;
   R_WS /= a_WS;

   c = exp(-R_WS*R_WS);

   f = a_WS*rho*(LexusData.SigmaNN);
   f *= 1.0 + w_WS*r_sqr/pow(R_WS, 2.);
   f /= sqrt(z_sqr)*(xi + c*exp(s*s)); 
 
 return f;
}/* NuInt3Gauss */



/* %%%%%%% 2 Parameter HO %%%%%%%%%%%% */


double Glauber::Anum2HO(double R_WS)
{
 int count=0;
 double up, down, a_WS, w_WS, rho, f;

 a_WS = Nuc_WS->a_WS;
 w_WS = Nuc_WS->w_WS; /* take this to be alpha */
 rho = Nuc_WS->rho_WS;

 down = 0.0;
 up = 1.0;
 
 f = integral(6, down, up, TOL, &count);
 f *= 4.0*M_PI*rho*pow(a_WS, 3.);

 return f;
}/* Anum2HO */


double Glauber::Anum2HOInt(double xi) 
{
 double y;
 double r_sqr, r;
 double w_WS;

 if(xi==0.0) xi = tiny;
 if(xi==1.0) xi = 1.0 - tiny;

 w_WS = Nuc_WS->w_WS;

/* already divided by a */
 
 r_sqr = -log(xi);
 
 r = sqrt(r_sqr);
 
/* 2 comes from dr^2 = 2 rdr */
 
 y = r + w_WS*r*r_sqr;
 y /= 2.0;
 
 return y;
}/* Anum2HOInt */


double Glauber::NuInt2HO(double xi)
{
 double f;
//  double argexp;
 double z_sqr, r_sqr, s;
 double w_WS, R_WS, a_WS, rho;
//  static int ind = 0;
 
 a_WS = Nuc_WS->a_WS;
 w_WS = Nuc_WS->w_WS;
 R_WS = Nuc_WS->R_WS;
 rho = Nuc_WS->rho_WS;

/* xi = exp(-z*z/a/a), r = sqrt(z^2 + s^2) */

   if ( xi == 0.0 ) xi = tiny; 
   if ( xi == 1.0 ) xi = 1.0-tiny; 

/* devide by a_WS, make life simpler */
   
   s = NuInS_S/a_WS;
   z_sqr = -log(xi);
   r_sqr = s*s + z_sqr;

/* no need to divide by 2 here because -infty < z < infty and 
   we integrate only over positive z */
   
   if(z_sqr < 0.0) z_sqr = tiny;
   f = a_WS*rho*(LexusData.SigmaNN);
   f *= (1.0 + w_WS*r_sqr)*exp(-s*s)/sqrt(z_sqr);
 
 return f;
}/* NuInt2HO */


double Glauber::integral (int id, double down, double up, double tol, int *count)
{
  double dx, y, g1[7];
  int i;
//   int j;
  
  if (down == up) y = 0.0;
  else
    {
      dx = (up-down)/6.0;
      for( i=0; i<7; i++) 
	{
	  if (id==1) g1[i] = NuInt2HO(down + i*dx);
	  else if (id==2) g1[i] = NuInt3Gauss(down + i*dx);
	  else if (id==3) g1[i] = NuInt3Fermi(down + i*dx);
	  else if (id==4) g1[i] = Anum3FermiInt(down + i*dx);
	  else if (id==5) g1[i] = Anum3GaussInt(down + i*dx);
	  else if (id==6) g1[i] = Anum2HOInt(down + i*dx);
	  else if (id==7) g1[i] = OLSIntegrand(down + i*dx);
	}
      *count = 7;
      y = qnc7(id, tol, down, dx, g1, 0.0, 0.0, count);
    }
  return y;
} /* end of integral */

double Glauber::qnc7(int id, double tol, double down, double dx, double *f_of, 
		     double pre_sum, double area, int *count)
{
  int i;
//   double x;
  double left_sum, right_sum, ans;
  static double w[] = 
    {41.0/140.0, 54.0/35.0, 27.0/140.0, 68.0/35.0, 27.0/140, 54.0/35.0,
       41.0/140.0};
  double fl[7];
  double fr[7];
  /*
    qnc7 calculates integral over left and right half of the given interval
    and branches
    to do so, first halve dx
    */
  
  dx /= 2.0;

  /*
    first calculate the left estimate
    f_of[] contains the evaluated values at down+i, 0< i <7
    store half distanced values for the left sum in fl[]
    */

  
  if (id==1)
    {
      fl[1] = NuInt2HO(down + dx);
      fl[3] = NuInt2HO(down + 3.0*dx);
      fl[5] = NuInt2HO(down + 5.0*dx);
    }
  else if (id==2)
    {
      fl[1] = NuInt3Gauss(down + dx);
      fl[3] = NuInt3Gauss(down + 3.0*dx);
      fl[5] = NuInt3Gauss(down + 5.0*dx);
    }
  else if (id==3)
    {
      fl[1] = NuInt3Fermi(down + dx);
      fl[3] = NuInt3Fermi(down + 3.0*dx);
      fl[5] = NuInt3Fermi(down + 5.0*dx);
    }
  else if (id==4)
    {
      fl[1] = Anum3FermiInt(down + dx);
      fl[3] = Anum3FermiInt(down + 3.0*dx);
      fl[5] = Anum3FermiInt(down + 5.0*dx);
    }  
  else if (id==5)
    {
      fl[1] = Anum3GaussInt(down + dx);
      fl[3] = Anum3GaussInt(down + 3.0*dx);
      fl[5] = Anum3GaussInt(down + 5.0*dx);
    }  
  else if (id==6)
    {
      fl[1] = Anum2HOInt(down + dx);
      fl[3] = Anum2HOInt(down + 3.0*dx);
      fl[5] = Anum2HOInt(down + 5.0*dx);
    }  
  else if (id==7)
    {
      fl[1] = OLSIntegrand(down + dx);
      fl[3] = OLSIntegrand(down + 3.0*dx);
      fl[5] = OLSIntegrand(down + 5.0*dx);
    }  
  
  fl[0] = f_of[0];
  fl[2] = f_of[1];
  fl[4] = f_of[2];
  fl[6] = f_of[3];
  
  *count += 3;

  left_sum = 0.0;
  for(i=0; i<7; i++) left_sum += w[i]*fl[i];
  left_sum *= dx;

/*printf("leftsum is %le\n", left_sum);*/

  /*
    like wise, the right sum is in fr[]
    */

  if (id==1)
    {
      fr[1] = NuInt2HO(down + 7.0*dx);
      fr[3] = NuInt2HO(down + 9.0*dx);
      fr[5] = NuInt2HO(down + 11.0*dx);
    }
  else if (id==2)
    {
      fr[1] = NuInt3Gauss(down + 7.0*dx);
      fr[3] = NuInt3Gauss(down + 9.0*dx);
      fr[5] = NuInt3Gauss(down + 11.0*dx);
    }
  else if (id==3)
    {
      fr[1] = NuInt3Fermi(down + 7.0*dx);
      fr[3] = NuInt3Fermi(down + 9.0*dx);
      fr[5] = NuInt3Fermi(down + 11.0*dx);
    }
  else if (id==4)
    {
      fr[1] = Anum3FermiInt(down + 7.0*dx);
      fr[3] = Anum3FermiInt(down + 9.0*dx);
      fr[5] = Anum3FermiInt(down + 11.0*dx);
    }
  else if (id==5)
    {
      fr[1] = Anum3GaussInt(down + 7.0*dx);
      fr[3] = Anum3GaussInt(down + 9.0*dx);
      fr[5] = Anum3GaussInt(down + 11.0*dx);
    }
  else if (id==6)
    {
      fr[1] = Anum2HOInt(down + 7.0*dx);
      fr[3] = Anum2HOInt(down + 9.0*dx);
      fr[5] = Anum2HOInt(down + 11.0*dx);
    }
  else if (id==7)
    {
      fr[1] = OLSIntegrand(down + 7.0*dx);
      fr[3] = OLSIntegrand(down + 9.0*dx);
      fr[5] = OLSIntegrand(down + 11.0*dx);
    }
  
  fr[0] = f_of[3];
  fr[2] = f_of[4];
  fr[4] = f_of[5];
  fr[6] = f_of[6];

  *count += 3;

  right_sum = 0.0;
  for(i=0; i<7; i++) right_sum += w[i]*fr[i];
  right_sum *= dx;

/*printf("rightsum is %le\n", right_sum);*/

  ans = left_sum + right_sum;

/*printf("ans is %le\n", ans);*/


  /* 
    up date total area subtract previously assigned area for this interval
    and add newly calculated area
    */
/*printf("pre_area is %le\n", area);*/

  area += -fabs(pre_sum) + fabs(left_sum) + fabs(right_sum);

  /* 
    printf("presum is %le\n", pre_sum);
    printf("area is %le\n", area);
    */
  /*
    branch if the refined sum is finer than the previous estimate
    */

  if( fabs(ans - pre_sum) > tol*fabs(area) && (*count < limit))
    {
      /*
	branch by calling the function itself
	by calling the qnc7 twice, we are branching
	since left hand side is being calculated first, until the condition
	is satisfied, the left branch keeps branching
	when finally the condition is met by one left-most interval, 
	qnc7 returns the right hand side of one up level,
	and the same process resumes 
	until the criterion is met by all the branched
	intervals,
	then qnc7 returns to the original right branch and resumes halving 
	until the condition is met by all intervals
	(funk, down, dx, f_of[7], pre_ans, ans)
	*/

      tol /= 1.414;
      left_sum = qnc7(id, tol, down, dx, fl, left_sum, area, count);
      right_sum = qnc7(id, tol, down+dx*6, dx, fr, right_sum, area, count);

      ans = left_sum + right_sum;

      /* printf("ans is %le\n", ans);*/

    }/* belongs to if*/

  return ans;
} /* end of qnc */


double Glauber::OLSIntegrand(double s)
{
  double sum, arg, x, r;
//   double sb;
  int k, m;
  m = 20;
  sum = 0.0;
  for(k=1; k<=m; k++)
    {
      arg = M_PI*(2.0*k - 1.0)/(2.0*m);
      x = cos(arg);
      r = sqrt(s*s + b*b + 2.0*s*b*x);
      sum += InterNuTInST(r);
    }/* k */
  
  return s*sum*M_PI/(1.0*m)*InterNuPInSP(s);
  
}/* OLSIntegrand */


double Glauber::TAB()
{
  double f;
  int count = 0;
  f = integral(7, 0.0, LexusData.SCutOff, TOL, &count); // integrate OLSIntegrand(s)
  f *= 2.0/(LexusData.SigmaNN); //here TAB is the number of binary collisions, dimensionless 
                                //(1/fm^4 integrated over dr_T^2 (gets rid of 1/fm^2), divided by sigma (gets rid of the other))
  return f;
}/* TAB */


double Glauber::PAB(double x, double y)
{
  double s1=sqrt(pow(x+b/2.,2.)+y*y);
  double s2=sqrt(pow(x-b/2.,2.)+y*y);
  return InterNuPInSP(s1)*InterNuTInST(s2)/(currentTAB*LexusData.SigmaNN);
}/* PAB */


void Glauber::initGlauber(double SigmaNN, string Target, string Projectile, double inb, int imax, int size, int rank)
{
  
  //remove("NuTInST.dat");
  //remove("NuPInSP.dat");

  string Target_Name;
  Target_Name = Target;

  string Projectile_Name;
  Projectile_Name = Projectile;

  //  char* p_name;
 
  string p_name;
  stringstream sp_name;

  const char* EOSPATH = "HYDROPROGRAMPATH";
  char* envPath = getenv(EOSPATH);
  if (envPath != 0 && *envPath != '\0') 
    {
      sp_name << envPath;
      sp_name << "/known_nuclei.dat";
      p_name = sp_name.str();
    }  
  else
    {
      sp_name << "./known_nuclei.dat";
      p_name = sp_name.str();
    }

  string paf;
  //char *paf;
  paf = p_name;

//   if(!util->IsFile(paf))
//     {
//       fprintf(stderr, "No known_nuclei.dat.  Please provide one.\n");
//       fprintf(stderr, "Exiting...\n");
//       exit(0);
//     }

  /* pp total cross-section : 40 mb = 4 fm**2 */
  /* LexusData.SigmaNN = 4.0; */
  
  /* energy unit is always GeV and length unit is fm */

  FindNucleusData(&(LexusData.Target), Target_Name, paf, rank);
  //PrintNucleusData(&(LexusData.Target));
  
  FindNucleusData(&(LexusData.Projectile), Projectile_Name, paf, rank);
  PrintNucleusData(&(LexusData.Projectile));

  LexusData.SigmaNN = 0.1*SigmaNN; // sigma in fm^2 
  currentA = LexusData.Projectile.A;
  
  /* Run Specific */
  
  LexusData.InterMax = imax;
  LexusData.SCutOff = 12.;
  
  b=inb; 
  PrintLexusData();
  //currentTAB=TAB();
}/* init */

double Glauber::areaTA(double x, double A)
{
  double f;
  f = A*220.*(1.-exp(-0.025*x*x)); 
  return f;
}

ReturnValue Glauber::SampleTARejection(Random *random)
{
  ReturnValue returnVec;

  double r, x, y, tmp;
//   double ratio;
  double phi;
  double A=1.2*LexusData.SigmaNN/4.21325504715; // increase the envelope for larger sigma_inel (larger root_s) (was originally written
  // for root(s)=200 GeV, hence the cross section of 4.21325504715 fm^2 (=42.13 mb)
  cout.precision(10);
  do
    {
      phi = 2.*M_PI*random->genrand64_real1();
      r = 6.32456*sqrt(-log((-0.00454545*(-220.*A+areaTA(15.,A)*random->genrand64_real1()))/A));
      // here random->genrand64_real1()*areaTA(LexusData.SCutOff) 
      // is a uniform random number on [0, area under f(x)]
      tmp = random->genrand64_real1();
      // x is uniform on [0,1]
      if ( r*InterNuPInSP(r) > A*r*11.*exp(-r*r/40.) ) 
	cout << "WARNING: TA>envelope: " << "TA=" << r*InterNuPInSP(r) 
	     << ", f=" << A*r*11.*exp(-r*r/40.) << endl;
    } while( tmp > r*InterNuPInSP(r)/(A*r*11.*exp(-r*r/40.))); 
  // reject if tmp is larger than the ratio p(y)/f(y), f(y)=A*r*11.*exp(-r*r/40.))
  x=r*cos(phi);
  y=r*sin(phi);
  returnVec.x=x;
  returnVec.y=y;
  returnVec.collided=0;
  return returnVec; 
}

