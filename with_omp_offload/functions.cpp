#include "functions.hpp"
#include <complex>
#include <iostream>
#include <cstdio>

const int TEAM_SIZE = 256; //More than 512 in GTX1050ti gives erroneous results

  

//Read configuration
bool chrg_conf(Particles& part, double L[3]){
  std::string line;
  int NC = part.get_Nkinds();
  int Ntot = part.get_Ntot();
  char name[30];
  sprintf(name, "config.dat");
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


//Coulomb potential using Ewald summation---------------------------------------


//REAL PART-----------------------------------------------------------------
void real_coulomb_kernel(double *rad, double *chrg, double *pos, int N, double L,
		  double alpha, double *Ur, double rcut, int *count) {
  int thrid = omp_get_thread_num();
  int nthr = omp_get_num_threads();
  int tmid = omp_get_team_num();
  int ntm = omp_get_num_teams();
  int i = thrid + nthr * tmid;

  int lcount = 0;
  double U = 0.;
  double ri, qi, rj, qj;
  double posij[3];
  
  if(i < N){
    ri = rad[i];
    qi = chrg[i];

    //even out the load among threads
    int mx = static_cast<int>(ceil(static_cast<float>(N-1)/2));
    if(fmod(static_cast<float>(N),static_cast<float>(2)) == 0. && i >= N/2)
      mx = static_cast<int>(floor(static_cast<float>(N-1)/2));

    int j = i+1 - N*(static_cast<int>(floor((i+1)/N + 0.5)));
    int cnt = 0;
    while(cnt < mx){
      rj = rad[j];
      qj = chrg[j];
      for(int a = 0; a < 3; ++a){
	posij[a] = pos[i*3+a] - pos[j*3+a];
	posij[a] -= L * floor(posij[a] / L + 0.5);
      }
      
      double RIJSQ = posij[0]*posij[0] + posij[1]*posij[1] + posij[2]*posij[2];
      if(RIJSQ < rcut*rcut){
	double RIJ = sqrt(RIJSQ);
	
	U += qi*qj * erfc(alpha*RIJ) / RIJ;
	lcount++;
      }
      j += 1 - N*static_cast<int>(floor((j+1)/N + 0.5));
      cnt++;
    }

    Ur[i] = U;
    count[i] = lcount; 
    //printf("Particle %d taken by thread %d/%d of team %d/%d.\n",i,thrid,nthr,tmid,ntm);
  }
}

double real_potential(const Particles &part, double L, double alpha,
		      double rcut) {
  int N = part.get_Ntot();

  //Get particle data arrays
  double* R = part.get_R();
  double* Q = part.get_Q();
  double* X = part.get_X();

  //Define No. of teams
  int Nteams = N / TEAM_SIZE;
  if(N % TEAM_SIZE) Nteams++;

  //Define output array
  double *Ur = new double[N];
  int *count = new int[N];

  //Prepare arrays in the device
#pragma omp target data map(to: R[:N], Q[:N], X[:3*N]) map(from: Ur[:N], count[:N])
  //Offload kernel
  #pragma omp target teams num_teams(Nteams)
  #pragma omp parallel num_threads(TEAM_SIZE)
  {    
    real_coulomb_kernel(R, Q, X, N, L, alpha, Ur, rcut, count); 
  }

  double tot_Ur = 0;
  int t_count = 0;
  for(int i = 0; i < N; i++){
    tot_Ur += Ur[i];
    t_count += count[i];
  }
  
  delete[] Ur;
  delete[] count;

  printf("No. interactions: %d\n", t_count);
  
  return tot_Ur;
}


//RECIPROCAL PART --------------------------------------------------------

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

double recip_coulomb(double *chrg, double *pos, int N, double kk2,
		     double K[][3]) {
  double U = 0;
  
  //To store the summations over the particles
  std::complex<double> cc(0, 0);
  std::complex<double> sum[4];
  for(int b = 0; b < 4; b++)
    sum[b] = cc;

  //Summation over particles
  for(int i = 0; i < N; i++){
    double ri[3];
    for(int a = 0; a < 3; ++a)
      ri[a] = pos[i*3+a];
    double qi = abs(chrg[i]);
    double shift = chrg[i] < 0 ? M_PI : 0;
    
    for(int b = 0; b < 4; b++){
      double rik = ri[0]*K[b][0] + ri[1]*K[b][1] + ri[2]*K[b][2]; 
      sum[b] += std::polar(qi, rik+shift);
    }
  }

  //Compute energy
  for(int b = 0; b < 4; b++){
    U += kk2 * 0.5 * std::norm(sum[b]);
  }
  
  return U;
}

void recip_coulomb_kernel(double *chrg, double *pos, int N, int kmax,
			  const double *Kvec, double L, double alpha,
			  double *Uk/*, int *cnt*/){
  int thrid = omp_get_thread_num();
  int nthr = omp_get_num_threads();
  int tmid = omp_get_team_num();
  int ntm = omp_get_num_teams();
  int kn = thrid + nthr * tmid;

  //int pcnt = 0;
  double U = 0.;
  
  double V = L*L*L;
  double P2 = 2*M_PI/L;
  int kmax2 = kmax*kmax;
  int kmax3 = kmax2*kmax;
  
  if(kn < kmax3){
    float knf = static_cast<float>(kn);
    int nz = static_cast<int>(floor(knf/kmax2));
    float l = knf - kmax2 * nz;
    int ny = static_cast<int>(floor(l/kmax));
    int nx = static_cast<int>(l) - kmax * ny;
    int nsq = nx*nx + ny*ny + nz*nz;

    if(nsq <= kmax2){ //if image is within a spherical shell...
      double kx = P2*nx;
      double ky = P2*ny;
      double kz = P2*nz;

      double kk2 = 2. * Kvec[kn] / V;  //mult by 2 for symmetry
    
      double K[4][3]; //Store kvectors 
      K[0][0] = kx; K[0][1] = ky; K[0][2] = kz;
      K[1][0] = -kx; K[1][1] = ky; K[1][2] = kz;
      K[2][0] = kx; K[2][1] = -ky; K[2][2] = kz;
      K[3][0] = -kx; K[3][1] = -ky; K[3][2] = kz;
      
      
      U += recip_coulomb(chrg, pos, N, kk2, K);
      //printf("U: %f\n", U);
      
      //Correct for the symmetries used
      if((nx == 0 && ny == 0) || (nx == 0 && nz == 0) || (ny == 0 && nz == 0))
	U /= 4;
      else if((nz == 0 && nx != 0 && ny != 0) || (ny == 0 && nz != 0 && nx != 0)
	      || (nx == 0 && ny != 0 && nz != 0))
	U /= 2;

      //pcnt += 1;
    }

    Uk[kn] = U;
  }
}

double recip_potential(const Particles &part, const Kvector &Kvec,
		       double L, double alpha, int kmax) {
  //int count = 0;
  int N = part.get_Ntot();

  //Get particle data arrays
  double* Q = part.get_Q();
  double* X = part.get_X();
  
  //Get Kvec array
  int ksz = Kvec.size();
  double *kvec = Kvec.get_all();
  kmax += 1;
  int kmax2 = kmax*kmax;
  int kmax3 = kmax2*kmax;

  //Define No. of teams
  int Nteams = N / TEAM_SIZE;
  if(N % TEAM_SIZE) Nteams++;

  //Define output array
  double *Uk = new double[kmax3]; 
  
  //Prepare arrays in the device
  #pragma omp target data map(to: Q[:N], X[:3*N], kvec[:ksz]) map(from: Uk[:kmax3])
  //Offload kernel
  #pragma omp target teams num_teams(Nteams)
  #pragma omp parallel num_threads(TEAM_SIZE)
  {    
     recip_coulomb_kernel(Q, X, N, kmax, kvec, L, alpha, Uk/*,cnt*/);
  }

  double tot_Uk = 0;
  for(int i = 0; i < kmax3; i++)
    tot_Uk += Uk[i];
  
  delete[] Uk;

  // int count = 0;
  // for(int i = 0; i < blocks; ++i)
  //   count += cnt[i];
  // delete[] cnt;
  // printf("count: %d\n",count);

  //Self energy of the main cell
  double self = 0;
  for(int i = 0; i < N; i++){
    double qi = part.get_charge(i);
    self += qi*qi;
  }
  self *= 0.5 * alpha / sqrt(M_PI);

  tot_Uk -= self;
  
  return tot_Uk;
}
