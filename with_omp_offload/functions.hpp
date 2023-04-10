#ifndef __FUNCTIONS_HPP
#define __FUNCTIONS_HPP

#include "classes.hpp"
#include <fstream>
#include <cmath>
#include <omp.h>
//#include <gsl/gsl_rng.h>

//Read configuration
bool chrg_conf(Particles& part, double L[3]);

// //Complex numbers
// #pragma omp declare target
// struct Complex {
//   double real;
//   double imag;

//   Complex() {}
//   Complex(double r, double i);

//   Complex& operator=(const Complex& c);
//   Complex& operator+=(const Complex& c);
// };
// #pragma omp declare target

// #pragma omp declare target
// Complex polar(double rho, double theta);
// #pragma omp declare target

// #pragma omp declare target
// double norm(Complex a);
// #pragma omp declare target


//Coulomb potential using Ewald summation--------------------------------------
double real_coulomb_kernel(const Particles &part, double L, int i, int j,
			   double rcut);
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
