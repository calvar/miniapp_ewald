#ifndef DEV_FUNCTIONS_HPP
#define DEV_FUNCTIONS_HPP

#include "classes.hpp"

double real_potential(const Particles &part, double L, double alpha);

double recip_potential(const Particles &part, const Kvector &Kvec, double L,
		       double alpha, int kmax);

#endif
