#include "functions.hpp"
#include <complex>
#include <iostream>
#include <cstdio>

//Read configuration
bool chrg_conf(Particles& part, double L[3]){
  std::string line;
  int NC = part.get_Nkinds();
  int Ntot = part.get_Ntot();
  char name[30];
  sprintf(name, "conf.dat");
  std::ifstream Con(name);
  if(! Con){
    Con.close();
    return false;
  }
  
  for(int n = 0; n < NC; ++n){ //charges
    double q;
    Con >> q;
  }
  Con >> L[0] >> L[1] >> L[2];
      
  for(int i = 0; i < Ntot; ++i){
    for(int a = 0; a < 3; ++a){
      double p;
      Con >> p;
      part.set_pos(i, a, p); 
    }
    for(int a = 0; a < 3; ++a){
      double m;
      Con >> m;
      part.set_mom(i, a, m);
    }
  }
  Con.close();
  
  return true;
}


//-------------------------------------------------------------------------

void k_vector(Kvector& Kvec, double L, double alpha, int kmax) {
  double c = 4 * alpha * alpha;
  double P2 = 2 * M_PI / L;
  int sz = kmax+1;
  int sz2 = sz*sz;
  for(int nz = 0; nz <= kmax; ++nz){
    double kz = P2 * nz;
    for(int ny = 0; ny <= kmax; ++ny){
      double ky = P2 * ny;
      for(int nx = 0; nx <= kmax; ++nx){
	double kx = P2 * nx;
	double ksq = kx*kx + ky*ky + kz*kz;
	if(ksq != 0.){
	  Kvec.set(nz*sz2 + ny*sz + nx, 4*M_PI*exp(-ksq/c) / ksq); 
	}else{
	  Kvec.set(nz*sz2 + ny*sz + nx, 0);
	}
      }
    }
  }
}
