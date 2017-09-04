/*
	
The program solves the following system of two first order differential equations.
	
These describe the radial mass distribution of a white dwarf star. 
 
m = Mo * m'
rho = rho_o * rho'
r = Ro * r' 


dm'
-  = r' ^ 2 * rho'
dr' 

drho'    - m'rho'
-    =  ----------
dr'        r'^2 gamma(rho')


                        rho'^2/3
gamma(rho') =      -----------------
			       3root(1+rho'^2/3)
			       
note;

m'(0) = 0
rho'(0) = rho_c
rho'(R) = 0 	

also, m' = 1/3 r'^3 rho_c


The step-size of the integrator is automatically adjusted by the controller.

*/

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

int func (double t, const double y[], double f[], void *params);

int main (void) 
{
  size_t neqs = 2;          /* number of equations */
  double eps_abs = 1.e-8, 
    eps_rel = 0.;           /* desired precision */
  double stepsize = 1e-6;   /* initial integration step */
  double rho_c;            /* central density */
  double r, r1 = 12000000000;    /* radius interval */
  int status;
  int i, np = 100;
  double step, rho_min = .1, rho_max = 1000000;
  
  step = (rho_max  - rho_min) / (np-1);
  
  for (i = 0; i<np;i++)
  
  {
  	rho_c = rho_min + i*step;
  	r = 0;
  /* 
   * Initial conditions 
   */
   
  double y[2] = { 0., rho_c };

  /*
   * Explicit embedded Runge-Kutta-Fehlberg (4,5) method. 
   * This method is a good general-purpose integrator.
   */
  gsl_odeiv2_step    *s = gsl_odeiv2_step_alloc 
                            (gsl_odeiv2_step_rkf45, neqs);
  gsl_odeiv2_control *c = gsl_odeiv2_control_y_new (eps_abs, 
						    eps_rel);
  gsl_odeiv2_evolve  *e = gsl_odeiv2_evolve_alloc (neqs);
    
  gsl_odeiv2_system sys = {func, NULL, neqs, &rho_c};
    
  /*
   * Evolution loop 
   */
   
  while ( (r < r1) && (y[1] > 0) )
  {
    status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &r,
                                      r1, &stepsize, y);

    if (status != GSL_SUCCESS) {
      printf ("Troubles: % .5e  % .5e  % .5e  % .5e\n", 
              rho_c, r, y[0], y[1]);
      break;
    }
    
  }

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
  
  printf ("% .5e  % .5e  % .5e  % .5e\n", rho_c, r, y[0], y[1]);
  
  }

  return 0;
}
