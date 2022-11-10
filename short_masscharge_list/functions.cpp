#include "functions.hpp"
#include <complex>
#include <iostream>

//Read configuration
bool chrg_conf(Particles& part, double L[3]){
  std::string line;
  int NC = part.get_Nkinds();
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
    part.set_charge(n, q);
  }
  Con >> L[0] >> L[1] >> L[2];
      
  for(int n = 0; n < NC; ++n){
    int NN = part.get_N(n);
    for(int i = 0; i < NN; i++){
      for(int a = 0; a < 3; ++a){
	double p;
	Con >> p;
	part.set_pos(n, i, a, p); 
      }
      for(int a = 0; a < 3; ++a){
	double m;
	Con >> m;
	part.set_mom(n, i, a, m);
      }
    }
  }
  Con.close();
  
  return true;
}



//Coulomb potential using Ewald summation---------------------------------------
double real_coulomb(const Particles &part, double L, int i, int j,
		    double alpha) {
  double U = 0;
  int m, mi, n, nj;
  part.get_kind(i, m, mi);
  part.get_kind(j, n, nj);
    
  double ri[3], rj[3];
  for(int a = 0; a < 3; ++a){
    ri[a] = part.get_pos(m, mi, a);
    rj[a] = part.get_pos(n, nj, a);
  }
  double rij[3];
  for(int a = 0; a < 3; ++a){
    rij[a] = ri[a] - rj[a];
    rij[a] -= L * floor(rij[a] / L + 0.5);
  }
  double RIJSQ = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
  
  double qi = part.get_charge(m); 
  double qj = part.get_charge(n);
  double r = sqrt(RIJSQ);
  
  U = qi*qj * erfc(alpha*r) / r;
  
  return U;
}

double real_potential(const Particles &part, double L, double alpha) {
  double Ur = 0;

  int N = 0;
  for(int i = 0; i < part.get_Nkinds(); ++i)
    N += part.get_N(i);
  
  for(int i = 0; i < N; i++){
    int mx = static_cast<int>(ceil(static_cast<float>(N-1)/2));
    if(fmod(static_cast<float>(N),2) == 0. && i >= N/2)
      mx = static_cast<int>(floor(static_cast<float>(N-1)/2));

    int j = i+1 - N*static_cast<int>(floor((i+1)/N + 0.5));
    int cnt = 0;
    while(cnt < mx){
      Ur += real_coulomb(part, L, i, j, alpha);
      
      j += 1 - N*static_cast<int>(floor((j+1)/N + 0.5));
      cnt++;
    }
  }
  return Ur;
}

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

double recip_coulomb(const Particles &part, int N, double kk2,
		     const double K[][3]) {
  double U = 0;
  
  //To store the summations over the particles
  std::complex<double> cc(0, 0);
  std::complex<double> sum[4];
  for(int b = 0; b < 4; b++)
    sum[b] = cc;

  //Summation over particles
  for(int i = 0; i < N; i++){
    int n, ni;
    part.get_kind(i, n, ni);

    double ri[3];
    for(int a = 0; a < 3; ++a)
      ri[a] = part.get_pos(n, ni, a);
    double qi = part.get_charge(n);

    for(int b = 0; b < 4; b++){
      double rik = ri[0]*K[b][0] + ri[1]*K[b][1] + ri[2]*K[b][2]; 
      sum[b] += std::polar(qi, rik);
    }
  }

  //Compute energy
  for(int b = 0; b < 4; b++){
    U += kk2 * 0.5 * std::norm(sum[b]);
  }
  
  return U;
}

double recip_potential(const Particles &part, const Kvector &Kvec,
		       double L, double alpha, int kmax) {
  double Uk = 0;

  int N = 0;
  for(int i = 0; i < part.get_Nkinds(); ++i)
    N += part.get_N(i);

  double V = L*L*L;
  double P2 = 2*M_PI/L;
  kmax += 1;
  int kmax2 = kmax*kmax;
  int kmax3 = kmax2*kmax;

  for(int kn = 0; kn < kmax3; ++kn){
    float knf = static_cast<float>(kn);
    int nz = static_cast<int>(floor(knf/kmax2));
    float l = knf - kmax2 * nz;
    int ny = static_cast<int>(floor(l/kmax));
    int nx = static_cast<int>(l) - kmax * ny;
    int nsq = nx*nx + ny*ny + nz*nz;

    double kx = P2*nx;
    double ky = P2*ny;
    double kz = P2*nz;

    if(nsq <= kmax2){ //if image is within a spherical shell...
      double kk2 = 2. * Kvec.get(kn) / V;  //mult by 2 for symmetry
      
      double K[4][3]; //Store kvectors 
      K[0][0] = kx; K[0][1] = ky; K[0][2] = kz;
      K[1][0] = -kx; K[1][1] = ky; K[1][2] = kz;
      K[2][0] = kx; K[2][1] = -ky; K[2][2] = kz;
      K[3][0] = -kx; K[3][1] = -ky; K[3][2] = kz;

      Uk += recip_coulomb(part, N, kk2, K);

      //Correct for the symmetries used
      if((nx == 0 && ny == 0) || (nx == 0 && nz == 0) || (ny == 0 && nz == 0))
	Uk /= 4;
      else if((nz == 0 && nx != 0 && ny != 0) || (ny == 0 && nz != 0 && nx != 0)
	      || (nx == 0 && ny != 0 && nz != 0))
	Uk /= 2;
    }
  }

  //Self energy of the main cell
  double self = 0;
  for(int i = 0; i < N; i++){
    int n, ni;
    part.get_kind(i, n, ni);
    double qi = part.get_charge(n);
    self += qi*qi;
  }
  self *= 0.5 * alpha / sqrt(M_PI);

  Uk -= self;
  
  return Uk;
}

//-------------------------------------------------------------------------
