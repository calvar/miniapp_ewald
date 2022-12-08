#ifndef __FUNCTIONS_HPP
#define __FUNCTIONS_HPP

#include "classes.hpp"
#include <fstream>
#include <cmath>
//#include <gsl/gsl_rng.h>

//Read configuration
bool chrg_conf(Particles& part, double L[3]);

//Coulomb potential using Ewald summation--------------------------------------

void k_vector(Kvector &Kvec, double L, double alpha, int kmax);

//----------------------------------------------------------------------

#endif
