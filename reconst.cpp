#include "util.h"
#include "data.h"
#include "grid.h"
#include "eos.h"
#include "reconst.h"

using namespace std;

Reconst::Reconst(EOS *eosIn, Grid *gridIn)
{
  eos = new EOS;
  eos = eosIn;
  grid = new Grid;
  grid = gridIn;
  util = new Util;
}

// destructor
Reconst::~Reconst()
{
  delete eos;
  delete grid;
  delete util;
}

int Reconst::ReconstIt(Grid *grid_p, int direc, double tau, double **uq, Grid *grid_pt,
		       double eps_init, double rhob_init, InitData *DATA, int rk_flag)
{
 /* reconstruct TJb from q[0] - q[4] */
 double K00, T00, J0, u[4], epsilon, p, h, rhob;
//  double tempk; 
 double epsilon_prev, rhob_prev, p_prev, p_guess, temperr;
 double epsilon_next, rhob_next, p_next, err, tempf, temph, cs2;
 double eps_guess, scalef;
 int iter, mu, nu, alpha;
 double q[5];
 const double RECONST_PRECISION = 1e-8;


/* prepare for the iteration */
/* uq = qiphL, qiphR, etc 
   qiphL[alpha][direc] means, for instance, TJ[alpha][0] 
   in the cell at x+dx/2 calculated from the left 
   */

/* uq are the conserved charges. That is, the ones appearing in
   d_tau (Ttautau/tau) + d_eta(Ttaueta/tau) + d_perp(Tperptau) = -Tetaeta
   d_tau (Ttaueta) + d_eta(Tetaeta) + d_v(tau Tveta) = -Ttaueta/tau 
   d_tau (Ttauv) + d_eta(Tetav) + d_w(tau Twv) = 0
   d_tau (Jtau) + d_eta Jeta + d_perp(tau Jperp) = 0
   */

/* q[0] = Ttautau/tau, q[1] = Ttaux, q[2] = Ttauy, q[3] = Ttaueta
   q[4] = Jtau */
/* uq = qiphL, qiphR, qimhL, qimhR, qirk */


 for(alpha=0; alpha<5; alpha++)
  {
   q[alpha] = uq[alpha][direc];
  }

 K00 = q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
 K00 /= (tau*tau);
 
 T00 = q[0]/tau;
 J0 = q[4]/tau;
 
 //if(q[4]<0)
 //  cout << q[4] << " direc=" << direc << endl;
//  tempk = T00-K00/T00;

 //remove if for speed
 // if(!finite(T00) || !finite(K00) || !finite(J0))
//   {
//    fprintf(stderr, "T00 = %e\n", T00);
//    fprintf(stderr, "K00 = %e\n", K00);
//    fprintf(stderr, "J0 = %e\n", J0);
//    fprintf(stderr, "q[0] = %e\n", q[0]);
//    fprintf(stderr, "q[1] = %e\n", q[1]);
//    fprintf(stderr, "q[2] = %e\n", q[2]);
//    fprintf(stderr, "q[3] = %e\n", q[3]);
//    fprintf(stderr, "q[4] = %e\n", q[4]);
//    exit(0);
//   }

 if( (T00 < 0.0) || ((T00 - K00/T00) < 0.0) || (T00 < (SMALL)) )
  {
  /* can't make Tmunu with this. restore the previous value */
  /* remember that uq are eigher halfway cells or the final q_next */
  /* at this point, the original values in grid_pt->TJb are not touched. */
    grid_p->epsilon = grid_pt->epsilon;
    grid_p->rhob = grid_pt->rhob;
    //rhob =grid_p->rhob ;
    //  if (rhob>1) cout << "rhob=" << rhob << endl;
     
    grid_p->p = grid_pt->p;
    for(mu=0; mu<4; mu++)
     {
      grid_p->TJb[0][4][mu] = grid_pt->TJb[rk_flag][4][mu];
      //cout << "grid_p->TJb[0][4][" << mu << "]=" << grid_p->TJb[0][4][mu] << endl;
      grid_p->u[0][mu] = grid_pt->u[rk_flag][mu];
 
      for(nu=0; nu<4; nu++)
       {
        grid_p->TJb[0][nu][mu] = grid_pt->TJb[rk_flag][nu][mu];
       }/* nu */
     }/* mu */
   return -1;
  }/* if t00-k00/t00 < 0.0 */

/* Iteration scheme */
 
 
 cs2 = eos->p_e_func(eps_init, rhob_init);
 eps_guess = GuessEps(T00, K00, cs2);
 epsilon_next = eps_guess;
 //cout << "epsilon_next=" << epsilon_next << endl;
 //cout << "rhob_init=" << rhob_init << endl;
     
 if(isnan(epsilon_next)) cout << "problem " << eps_guess << " T00=" << T00 << " K00=" << K00 << " cs2=" 
			      << cs2 << " q[0]=" << q[0] << " uq[0][" << direc << "]=" << uq[0][direc] 
			      << " q[1]=" << q[1] << " q[2]=" << q[2] << endl;
 p_guess = eos->get_pressure(epsilon_next, rhob_init);
 p_next = p_guess;

/* rhob = J0*sqrt( (eps + p)/(T00 + p) ) */

 if(J0 == 0.0)
  {rhob_next = 0.0;}
 else
  { rhob_next = J0*sqrt((eps_guess + p_guess)/(T00 + p_guess));
    if(!finite(rhob_next) || rhob_next<0) 
    {
     rhob_next = 0.0;
    }
  }

 err = 1.0;
 for(iter=0; iter<100; iter++)
  {
   if(err < (RECONST_PRECISION)*0.01)
    {
      if(isnan(epsilon_next)) cout << "problem2" << endl;
      p_next = eos->get_pressure(epsilon_next, rhob_next);
     break;
    }
   else
    {
     epsilon_prev = epsilon_next;
     rhob_prev = rhob_next;
    }
   
   //   cout << "epsilon_next=" << epsilon_prev << endl;
   //cout << "rhob_next=" << rhob_prev << endl;
     
   if(isnan(epsilon_prev)) cout << "problem3" << endl;

   p_prev = eos->get_pressure(epsilon_prev, rhob_prev);
   epsilon_next = T00 - K00/(T00 + p_prev);
  
   err = 0.0;
   
   if(DATA->turn_on_rhob == 1)
    {
     rhob_next = J0*sqrt((epsilon_prev + p_prev)/(T00 + p_prev));
     temperr = fabs((rhob_next-rhob_prev)/(rhob_prev+SMALL));
     if(finite(temperr)) err += temperr;
     else if(isinf(temperr)) err += 1000.0; /* big enough */
     else if(isnan(temperr)) err += 1000.0;
    }
   else
    {
     rhob_next = 0.0;
    }
   
   temperr = fabs((epsilon_next-epsilon_prev)/(epsilon_prev+SMALL));
   if(finite(temperr)) err += temperr;
   else if(isinf(temperr)) err += 1000.0; /* big enough */
   else if(isnan(temperr)) err += 1000.0;

  }/* iter */ 

  if(iter == 100)
    {
      fprintf(stderr, "Reconst didn't converge.\n");
      cout << grid_p->epsilon << endl;
      fprintf(stderr, "Reverting to the previous TJb...\n"); 
      for(mu=0; mu<4; mu++)
       {
        grid_p->TJb[0][4][mu] = grid_pt->TJb[rk_flag][4][mu];
        grid_p->u[0][mu] = grid_pt->u[rk_flag][mu];
        
        for(nu=0; nu<4; nu++)
         {
          grid_p->TJb[0][nu][mu] = grid_pt->TJb[rk_flag][nu][mu];
         }/* nu */
       }/* mu */
     
     grid_p->epsilon = grid_pt->epsilon;
     grid_p->p = grid_pt->p;
     grid_p->rhob = grid_pt->rhob;
     return -2;
    }/* if iteration is unsuccessful, revert */

/* update */
 
 epsilon = grid_p->epsilon = epsilon_next;
 p = grid_p->p = p_next;
 rhob = grid_p->rhob = rhob_next;
 h = p+epsilon;
 //remove if for speed
 // if( (epsilon < 0.0) || !finite(epsilon) )
//   {
//    fprintf(stderr, "In Reconst, reconstructed epsilon = %e.\n", epsilon);
//    fprintf(stderr, "Can't happen.\n"); 
//    exit(0);
//   }

/* q[0] = Ttautau/tau, q[1] = Ttaux, q[2] = Ttauy, q[3] = Ttaueta,
   q[4] = Jtau */

 u[0] = sqrt((q[0]/tau + p)/h);
 //remove if for speed
 if(!finite(u[0]))
  {
   u[0] = 1.0;
   u[1] = 0.0;
   u[2] = 0.0;
   u[3] = 0.0;
  }
 else
  {
   u[1] = q[1]/tau/h/u[0]; 
   u[2] = q[2]/tau/h/u[0]; 
   u[3] = q[3]/tau/h/u[0]; 
  }

   if(u[0] > cosh(DATA->local_y_max))
    {
      fprintf(stderr, "Reconst: u[0] = %e is too large.\n", u[0]);
      if(grid_pt->epsilon > 0.3)
	{
	  fprintf(stderr, "Reconst: u[0] = %e is too large.\n", u[0]);
	  fprintf(stderr, "epsilon = %e\n", grid_pt->epsilon);
	  fprintf(stderr, "Reverting to the previous TJb...\n"); 
	}
      for(mu=0; mu<4; mu++)
	{
	  grid_p->TJb[0][4][mu] = grid_pt->TJb[rk_flag][4][mu];
	  grid_p->u[0][mu] = grid_pt->u[rk_flag][mu];
	  
	  for(nu=0; nu<4; nu++)
	    {
	      grid_p->TJb[0][nu][mu] = grid_pt->TJb[rk_flag][nu][mu];
	    }/* nu */
	}/* mu */
      grid_p->epsilon = grid_pt->epsilon;
      grid_p->p = grid_pt->p;
      grid_p->rhob = grid_pt->rhob;
      
      return -2;
    }/* if u[0] is too large, revert */



 /* Correcting normalization of 4-velocity */
 temph = u[0]*u[0] - u[1]*u[1] - u[2]*u[2] - u[3]*u[3];
 //Correct velocity when unitarity is not satisfied to numerical accuracy (constant "SMALL")
 if(fabs(temph - 1.0) > SMALL)
 {
  //If the deviation is too large, exit MUSIC
  if(fabs(temph - 1.0) > 0.1)
   {
    fprintf(stderr, "In Reconst, reconstructed : u2 = %e\n", temph);
    fprintf(stderr, "Can't happen.\n");
    exit(0);
   }
  //Warn only when the deviation from 1 is relatively large
  else if(fabs(temph - 1.0) > sqrt(SMALL))
  {
      fprintf(stderr, "In Reconst, reconstructed : u2 = %e\n", temph);
      fprintf(stderr, "with u[0] = %e\n", u[0]);
      fprintf(stderr, "Correcting it...\n");
  }   

  //Rescaling spatial components of velocity so that unitarity is exactly satisfied
  //(u[0] is not modified)
  scalef = (u[0]-1.0);
  scalef *= (u[0]+1.0);
  scalef /= (u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
  scalef = sqrt(scalef);
  u[1] *= scalef;
  u[2] *= scalef;
  u[3] *= scalef;
   
 }/* if u^mu u_\mu != 1 */
 /* End: Correcting normalization of 4-velocity */

 for(mu=0; mu<4; mu++)
  {
   
   tempf = grid_p->TJb[0][4][mu] = rhob*u[mu];
   //remove if for speed
 //   if(!finite(tempf)) 
//      {
//        fprintf(stderr, "Update: Jb[%d] is %e.\n", mu, rhob*u[mu]);
//        fprintf(stderr, "Update: rhob is %e.\n", rhob);
//        fprintf(stderr, "Update: u[%d] is %e.\n", mu, u[mu]);
//        exit(0);
//      }
   
   tempf = grid_p->u[0][mu] = u[mu];
   
   //remove if for speed
 //   if(!finite(tempf)) {
//     fprintf(stderr, "Update: u[%d] is %e.\n", mu, u[mu]);
//     exit(0);
//     }
   
   for(nu=0; nu<4; nu++)
    {
     tempf = grid_p->TJb[0][nu][mu] 
     = ((epsilon + p)*u[nu]*u[mu] + p*(DATA->gmunu)[nu][mu]);
     if(!finite(tempf)) {
       fprintf(stderr, "Update: TJb[0][%d][%d] is %e.\n",
                        nu, mu, grid_p->TJb[0][nu][mu]);
       fprintf(stderr, "Update: epsilon is %e.\n", epsilon);
       exit(0);
       }
    }/* nu */
  }/* mu */

 return 1; /* on successful execution */
}/* Reconst */

/* reconstruct TJb from q[0] - q[4] */
/* reconstruct velocity first for finite mu_B case (add by C. Shen Nov. 2014) */
int Reconst::ReconstIt_velocity(Grid *grid_p, int direc, double tau, double **uq, Grid *grid_pt,
		                    double eps_init, double rhob_init, InitData *DATA, int rk_flag)
{
   double K00, T00, J0, u[4], epsilon, p, enthalpy, rhob;
   int iter, mu, nu, alpha;
   double q[5];
   const double RECONST_PRECISION = 1e-8;

/* prepare for the iteration */
/* uq = qiphL, qiphR, etc 
   qiphL[alpha][direc] means, for instance, TJ[alpha][0] 
   in the cell at x+dx/2 calculated from the left */

/* uq are the conserved charges. That is, the ones appearing in
   d_tau (Ttautau/tau) + d_eta(Ttaueta/tau) + d_perp(Tperptau) = -Tetaeta
   d_tau (Ttaueta) + d_eta(Tetaeta) + d_v(tau Tveta) = -Ttaueta/tau 
   d_tau (Ttauv) + d_eta(Tetav) + d_w(tau Twv) = 0
   d_tau (Jtau) + d_eta Jeta + d_perp(tau Jperp) = 0 */

/* q[0] = Ttautau/tau, q[1] = Ttaux, q[2] = Ttauy, q[3] = Ttaueta
   q[4] = Jtau */
/* uq = qiphL, qiphR, qimhL, qimhR, qirk */


 for(alpha=0; alpha<5; alpha++)
  {
   q[alpha] = uq[alpha][direc];
  }

 K00 = q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
 K00 /= (tau*tau);
 
 T00 = q[0]/tau;
 J0 = q[4]/tau;
 
 //if(q[4]<0)
 //  cout << q[4] << " direc=" << direc << endl;
//  tempk = T00-K00/T00;

 //remove if for speed
 // if(!finite(T00) || !finite(K00) || !finite(J0))
//   {
//    fprintf(stderr, "T00 = %e\n", T00);
//    fprintf(stderr, "K00 = %e\n", K00);
//    fprintf(stderr, "J0 = %e\n", J0);
//    fprintf(stderr, "q[0] = %e\n", q[0]);
//    fprintf(stderr, "q[1] = %e\n", q[1]);
//    fprintf(stderr, "q[2] = %e\n", q[2]);
//    fprintf(stderr, "q[3] = %e\n", q[3]);
//    fprintf(stderr, "q[4] = %e\n", q[4]);
//    exit(0);
//   }

 if( (T00 < 0.0) || ((T00 - K00/T00) < 0.0) || (T00 < (SMALL)) )
  {
  /* can't make Tmunu with this. restore the previous value */
  /* remember that uq are eigher halfway cells or the final q_next */
  /* at this point, the original values in grid_pt->TJb are not touched. */
    grid_p->epsilon = grid_pt->epsilon;
    grid_p->rhob = grid_pt->rhob;
    //rhob =grid_p->rhob ;
    //  if (rhob>1) cout << "rhob=" << rhob << endl;
     
    grid_p->p = grid_pt->p;
    for(mu=0; mu<4; mu++)
     {
      grid_p->TJb[0][4][mu] = grid_pt->TJb[rk_flag][4][mu];
      //cout << "grid_p->TJb[0][4][" << mu << "]=" << grid_p->TJb[0][4][mu] << endl;
      grid_p->u[0][mu] = grid_pt->u[rk_flag][mu];
 
      for(nu=0; nu<4; nu++)
       {
        grid_p->TJb[0][nu][mu] = grid_pt->TJb[rk_flag][nu][mu];
       }/* nu */
     }/* mu */
   return -1;
  }/* if t00-k00/t00 < 0.0 */

/* Iteration scheme */
 
 
 cs2 = eos->p_e_func(eps_init, rhob_init);
 eps_guess = GuessEps(T00, K00, cs2);
 epsilon_next = eps_guess;
 //cout << "epsilon_next=" << epsilon_next << endl;
 //cout << "rhob_init=" << rhob_init << endl;
     
 if(isnan(epsilon_next)) cout << "problem " << eps_guess << " T00=" << T00 << " K00=" << K00 << " cs2=" 
			      << cs2 << " q[0]=" << q[0] << " uq[0][" << direc << "]=" << uq[0][direc] 
			      << " q[1]=" << q[1] << " q[2]=" << q[2] << endl;
 p_guess = eos->get_pressure(epsilon_next, rhob_init);
 p_next = p_guess;

/* rhob = J0*sqrt( (eps + p)/(T00 + p) ) */

 if(J0 == 0.0)
  {rhob_next = 0.0;}
 else
  { rhob_next = J0*sqrt((eps_guess + p_guess)/(T00 + p_guess));
    if(!finite(rhob_next) || rhob_next<0) 
    {
     rhob_next = 0.0;
    }
  }

 err = 1.0;
 for(iter=0; iter<100; iter++)
  {
   if(err < (RECONST_PRECISION)*0.01)
    {
      if(isnan(epsilon_next)) cout << "problem2" << endl;
      p_next = eos->get_pressure(epsilon_next, rhob_next);
     break;
    }
   else
    {
     epsilon_prev = epsilon_next;
     rhob_prev = rhob_next;
    }
   
   //   cout << "epsilon_next=" << epsilon_prev << endl;
   //cout << "rhob_next=" << rhob_prev << endl;
     
   if(isnan(epsilon_prev)) cout << "problem3" << endl;

   p_prev = eos->get_pressure(epsilon_prev, rhob_prev);
   epsilon_next = T00 - K00/(T00 + p_prev);
  
   err = 0.0;
   
   if(DATA->turn_on_rhob == 1)
    {
     rhob_next = J0*sqrt((epsilon_prev + p_prev)/(T00 + p_prev));
     temperr = fabs((rhob_next-rhob_prev)/(rhob_prev+SMALL));
     if(finite(temperr)) err += temperr;
     else if(isinf(temperr)) err += 1000.0; /* big enough */
     else if(isnan(temperr)) err += 1000.0;
    }
   else
    {
     rhob_next = 0.0;
    }
   
   temperr = fabs((epsilon_next-epsilon_prev)/(epsilon_prev+SMALL));
   if(finite(temperr)) err += temperr;
   else if(isinf(temperr)) err += 1000.0; /* big enough */
   else if(isnan(temperr)) err += 1000.0;

  }/* iter */ 

  if(iter == 100)
    {
      fprintf(stderr, "Reconst didn't converge.\n");
      cout << grid_p->epsilon << endl;
      fprintf(stderr, "Reverting to the previous TJb...\n"); 
      for(mu=0; mu<4; mu++)
       {
        grid_p->TJb[0][4][mu] = grid_pt->TJb[rk_flag][4][mu];
        grid_p->u[0][mu] = grid_pt->u[rk_flag][mu];
        
        for(nu=0; nu<4; nu++)
         {
          grid_p->TJb[0][nu][mu] = grid_pt->TJb[rk_flag][nu][mu];
         }/* nu */
       }/* mu */
     
     grid_p->epsilon = grid_pt->epsilon;
     grid_p->p = grid_pt->p;
     grid_p->rhob = grid_pt->rhob;
     return -2;
    }/* if iteration is unsuccessful, revert */

/* update */
 
 epsilon = grid_p->epsilon = epsilon_next;
 p = grid_p->p = p_next;
 rhob = grid_p->rhob = rhob_next;
 h = p+epsilon;
 //remove if for speed
 // if( (epsilon < 0.0) || !finite(epsilon) )
//   {
//    fprintf(stderr, "In Reconst, reconstructed epsilon = %e.\n", epsilon);
//    fprintf(stderr, "Can't happen.\n"); 
//    exit(0);
//   }

/* q[0] = Ttautau/tau, q[1] = Ttaux, q[2] = Ttauy, q[3] = Ttaueta,
   q[4] = Jtau */

 u[0] = sqrt((q[0]/tau + p)/h);
 //remove if for speed
 if(!finite(u[0]))
  {
   u[0] = 1.0;
   u[1] = 0.0;
   u[2] = 0.0;
   u[3] = 0.0;
  }
 else
  {
   u[1] = q[1]/tau/h/u[0]; 
   u[2] = q[2]/tau/h/u[0]; 
   u[3] = q[3]/tau/h/u[0]; 
  }

   if(u[0] > cosh(DATA->local_y_max))
    {
      fprintf(stderr, "Reconst: u[0] = %e is too large.\n", u[0]);
      if(grid_pt->epsilon > 0.3)
	{
	  fprintf(stderr, "Reconst: u[0] = %e is too large.\n", u[0]);
	  fprintf(stderr, "epsilon = %e\n", grid_pt->epsilon);
	  fprintf(stderr, "Reverting to the previous TJb...\n"); 
	}
      for(mu=0; mu<4; mu++)
	{
	  grid_p->TJb[0][4][mu] = grid_pt->TJb[rk_flag][4][mu];
	  grid_p->u[0][mu] = grid_pt->u[rk_flag][mu];
	  
	  for(nu=0; nu<4; nu++)
	    {
	      grid_p->TJb[0][nu][mu] = grid_pt->TJb[rk_flag][nu][mu];
	    }/* nu */
	}/* mu */
      grid_p->epsilon = grid_pt->epsilon;
      grid_p->p = grid_pt->p;
      grid_p->rhob = grid_pt->rhob;
      
      return -2;
    }/* if u[0] is too large, revert */



 /* Correcting normalization of 4-velocity */
 temph = u[0]*u[0] - u[1]*u[1] - u[2]*u[2] - u[3]*u[3];
 //Correct velocity when unitarity is not satisfied to numerical accuracy (constant "SMALL")
 if(fabs(temph - 1.0) > SMALL)
 {
  //If the deviation is too large, exit MUSIC
  if(fabs(temph - 1.0) > 0.1)
   {
    fprintf(stderr, "In Reconst, reconstructed : u2 = %e\n", temph);
    fprintf(stderr, "Can't happen.\n");
    exit(0);
   }
  //Warn only when the deviation from 1 is relatively large
  else if(fabs(temph - 1.0) > sqrt(SMALL))
  {
      fprintf(stderr, "In Reconst, reconstructed : u2 = %e\n", temph);
      fprintf(stderr, "with u[0] = %e\n", u[0]);
      fprintf(stderr, "Correcting it...\n");
  }   

  //Rescaling spatial components of velocity so that unitarity is exactly satisfied
  //(u[0] is not modified)
  scalef = (u[0]-1.0);
  scalef *= (u[0]+1.0);
  scalef /= (u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
  scalef = sqrt(scalef);
  u[1] *= scalef;
  u[2] *= scalef;
  u[3] *= scalef;
   
 }/* if u^mu u_\mu != 1 */
 /* End: Correcting normalization of 4-velocity */

 for(mu=0; mu<4; mu++)
  {
   
   tempf = grid_p->TJb[0][4][mu] = rhob*u[mu];
   //remove if for speed
 //   if(!finite(tempf)) 
//      {
//        fprintf(stderr, "Update: Jb[%d] is %e.\n", mu, rhob*u[mu]);
//        fprintf(stderr, "Update: rhob is %e.\n", rhob);
//        fprintf(stderr, "Update: u[%d] is %e.\n", mu, u[mu]);
//        exit(0);
//      }
   
   tempf = grid_p->u[0][mu] = u[mu];
   
   //remove if for speed
 //   if(!finite(tempf)) {
//     fprintf(stderr, "Update: u[%d] is %e.\n", mu, u[mu]);
//     exit(0);
//     }
   
   for(nu=0; nu<4; nu++)
    {
     tempf = grid_p->TJb[0][nu][mu] 
     = ((epsilon + p)*u[nu]*u[mu] + p*(DATA->gmunu)[nu][mu]);
     if(!finite(tempf)) {
       fprintf(stderr, "Update: TJb[0][%d][%d] is %e.\n",
                        nu, mu, grid_p->TJb[0][nu][mu]);
       fprintf(stderr, "Update: epsilon is %e.\n", epsilon);
       exit(0);
       }
    }/* nu */
  }/* mu */

 return 1; /* on successful execution */
}/* Reconst */



void Reconst::ReconstError(const char *str, int i, int rk_flag, double *qi, double **qi2, Grid *grid_pt)
{
 int alpha;

 fprintf(stderr, "Reconst %s in the direction = %d reports an error.\n", 
         str, i); 
 fprintf(stderr, "grid_pt position = (%d, %d, %d).\n", 
         grid_pt->position[1], grid_pt->position[2], grid_pt->position[3]);
 fprintf(stderr, "rk_flag = %d\n", rk_flag); 
 
 for(alpha=0; alpha<5; alpha++)
  { fprintf(stderr, "qi[%d] = %e\n", alpha, qi[alpha]); }
 for(alpha=0; alpha<5; alpha++)
  { fprintf(stderr, "qi2[%d][%d] = %e\n", alpha, i, qi2[alpha][i]); }
 return;
}/* ReconstErr */


double Reconst::GuessEps(double T00, double K00, double cs2)
{
 double f;
 
 if(cs2 < SMALL)
  {
   f = ((-K00 + util->Power(T00,2))*(cs2*K00*util->Power(T00,2) + util->Power(T00,4) 
      + util->Power(cs2,2)*K00*(2*K00 - util->Power(T00,2))))/util->Power(T00,5);
  }
 else
 {
  f = ((-1.0 + cs2)*T00 
     + sqrt(-4.0*cs2*K00 + util->Power(T00,2) + 2.0*cs2*util->Power(T00,2)
     + util->Power(cs2,2)*util->Power(T00,2)))/(2.0*cs2);
 }
 return f;
}/*  GuessEps */


