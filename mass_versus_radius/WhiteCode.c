/*
 * The function provides the RHS of the following system of two
 * first order dimensionlesss differential equations, which describe 
 * the radial distribution of density, rho(r) and mass, m(r), inside the star.
 *
 *   m' = rho' r^2 
 *   rho' = -m rho' / (gamma(rho') r^2)
 *
 *
 * Here gamma(rho') = rho'^(2/3)/(3 \sqrt(1 + rho'^(2/3))
 * 
 */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

int func (double r, const double y[], double f[], void *params)
 
{

  double rho_c = *(double *) params;

  double m = y[0];
  double rho = y[1];
  
  if (rho>0)
  {
  
  	double p = pow(rho, 2./3.);
  	double gamma = p / ( 3. * sqrt(1.+p)) ;
  	f[0] = rho*r*r;
  	if (r>1.e-6)
  	{
  		f[1] = -m * rho / (gamma * r * r);
  	}
  	else
  	{
  		f[1] = -rho_c/3. * r * rho/gamma;
  	}
  }
  else 
  {
  	f[0] = 0.;
  	f[1] = 0.;
  }
  
  return GSL_SUCCESS;
}
