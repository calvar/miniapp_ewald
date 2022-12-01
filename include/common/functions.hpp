#ifndef __FUNCTIONS_HPP
#define __FUNCTIONS_HPP

#include "classes.hpp"
#include <fstream>
#include <cmath>
// #include <gsl/gsl_rng.h>

//Read configuration
bool chrg_conf(Particles& part, double L[3]);

//Coulomb potential using Ewald summation--------------------------------------
double real_coulomb(const Particles &part, double L, int i, int j);
double real_potential(const Particles &part, double L, double alpha);

void k_vector(Kvector &Kvec, double L, double alpha, int kmax);
double recip_coulomb(const Particles &part, int N, double kk2,
		     const double K[][3]);
double recip_potential(const Particles &part, const Kvector &Kvec,
		       double L, double alpha, int kmax);
//----------------------------------------------------------------------

#endif
