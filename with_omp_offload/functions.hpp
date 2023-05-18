#ifndef __FUNCTIONS_HPP
#define __FUNCTIONS_HPP

#include "classes.hpp"
#include <fstream>
#include <cmath>
#include <omp.h>
//#include <gsl/gsl_rng.h>

// //COMPLEX NUMBERS ----------------------------------------------------
// //#pragma omp declare target
// struct Complex {
//   double real;
//   double imag;

//   Complex() {}
//   Complex(double r, double i){ real = r; imag = i; }

//   Complex& operator=(const Complex& c){
//     real = c.real;
//     imag = c.imag;
//     return *this;
//   }
//   Complex& operator+=(const Complex& c){
//     real += c.real;
//     imag += c.imag;
//     return *this;
//   }
// };
// //#pragma omp declare target

// //#pragma omp declare target
// Complex polar(double rho, double theta) {
//   Complex c;
//   c.real = rho * cos(theta);
//   c.imag = rho * sin(theta);
//   return c;
// }
// //#pragma omp declare target

// //#pragma omp declare target
// double norm(Complex a) {
//   return sqrt(a.real*a.real + a.imag*a.imag);
// }
// //#pragma omp declare target




//Read configuration
bool chrg_conf(Particles& part, double L[3]);

//Coulomb potential using Ewald summation--------------------------------------
//double real_coulomb_kernel(const Particles &part, double L, int i, int j,
//			   double rcut);
void real_coulomb_kernel(double *rad, double *chrg, double *pos, int N, double L,
			 double alpha, double *Ur, double rcut, int *count);
double real_potential(const Particles &part, double L, double alpha,
		      double rcut);

void k_vector(Kvector &Kvec, double L, double alpha, int kmax);
double recip_coulomb(double *chrg, double *pos, int N, double kk2,
		     double K[][3]);
void recip_coulomb_kernel(double *chrg, double *pos, int N, int kmax,
			  const double *Kvec, double L, double alpha,
			  double *Uk/*, int *cnt*/);
double recip_potential(const Particles &part, const Kvector &Kvec,
		       double L, double alpha, int kmax);
//----------------------------------------------------------------------

#endif
