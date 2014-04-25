/* Esta basado en el antiguo que hiciste de prueba. No usaremos armadillo, sino
   gsl y vector... */

#include <vector>
#include <gsl/gsl_random>
#include <simplectic.h>

void PointIn3Sphere(simplectic &x,simplectic &y, double radio, gsl_rng *r ){
  // Da un punto aleatorio en la superficie de S3. 
  // El punto es un par simplectico.
  // Radio es el radio. 
  double xq,xp, yq,yp;
  double magnitud;
  
  xq= gsl_ran_gaussian(gsl_rng * r, 1.00);
  xp= gsl_ran_gaussian(gsl_rng * r, 1.00);
  yq= gsl_ran_gaussian(gsl_rng * r, 1.00);  
  yp= gsl_ran_gaussian(gsl_rng * r, 1.00);
  
  magnitud=sqrt(xq*xq+xp*xp+yq*yq+yq*yq);
  
  xq/=magnitud;
  


};


