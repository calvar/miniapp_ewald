#ifndef __FUNCTIONS_HPP
#define __FUNCTIONS_HPP

#include "classes.hpp"
#include <sycl/sycl.hpp>
#include <fstream>
#include <cmath>
//#include <gsl/gsl_rng.h>

//Read configuration
bool chrg_conf(Particles& part, double L[3]);

//Coulomb potential using Ewald summation--------------------------------------
double real_coulomb(double *Q, double *X, double L, int i, int j,
		    double alpha, double rcut, int &count);
double real_potential(Particles &part, double L, double alpha,
		      double rcut);

void k_vector(Kvector &Kvec, double L, double alpha, int kmax);
double recip_coulomb(double *Q, double *X, int N, double kk2,
		     const double K[][3]);
double recip_potential(Particles &part, Kvector &Kvec,
		       double L, double alpha, int kmax);
//----------------------------------------------------------------------

#endif