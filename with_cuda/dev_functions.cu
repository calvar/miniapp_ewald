#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cuda_profiler_api.h>
#include "classes.hpp"
#include "dev_functions.hpp"

#define BLOCK_WIDTH 4
#define K_SIM_SIZE 4


//Implementation from the CUDA page: https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#atomic-functions
__device__ double atomicAdd(double* address, double val){
  unsigned long long int* address_as_ull =
    (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
		    __double_as_longlong(val +
					 __longlong_as_double(assumed)));
    
    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
  } while (assumed != old);
  
  return __longlong_as_double(old);
}

//Atomic add for cuDoubleComplex
__device__ void atAddComplex(cuDoubleComplex* a, cuDoubleComplex b){
  //transform the addresses of real and imag. parts to double pointers
  double *x = (double*)a;
  double *y = x+1; //add 1*sizeof(double) bits

  //use atomicAdd for double variables
  atomicAdd(x, cuCreal(b));
  atomicAdd(y, cuCimag(b));
}



//Coulomb potential using Ewald summation---------------------------------------
__global__ void real_coulomb_kernel(double *rad, double *chrg, double *pos,
				    int N, double L, double alpha, double *Ur) {

  int i = threadIdx.x + blockDim.x * blockIdx.x;
  //shared variables for accumulated values
  __shared__ double sum_ur[BLOCK_WIDTH];
  __syncthreads();
  
  double U = 0.;
  double ri, qi, rj, qj;
  double posij[3];

  if(i < N-1){
    ri = rad[i];
    qi = chrg[i];

    //even out the load among threads
    int mx = static_cast<int>(ceil(static_cast<float>(N-1)/2));
    if(fmod(static_cast<float>(N),2) == 0. && i >= N/2)
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
      double RIJ = sqrt(RIJ);

      U += qi*qj * erfc(alpha*RIJ) / RIJ;

      j += 1 - N*static_cast<int>(floor((j+1)/N + 0.5));
      cnt++;
    }
  }

  //Save local contribution
  sum_ur[threadIdx.x] = U;
  __syncthreads();
  //Add contributions
  if(threadIdx.x == 0){
    Ur[blockIdx.x] = 0.;
    for(int k = 0; k < BLOCK_WIDTH; i++)
      Ur += sum_ur[k];
  }
}

double real_potential(const Particles &part, double L, double alpha) {
  int N = part.get_Ntot();

  double *d_rad, *d_chrg, *d_pos, *d_Ur;

  //Get particle data arrays
  double* R = part.get_R();
  double* Q = part.get_Q();
  double* X = part.get_X();

  //Set service
  cudaService(0);

  //Allocate and copy memory to device
  cudaMalloc((void**)&d_rad, sizeof(double)*N);
  cudaMemcpy(d_rad, R, sizeof(double)*N, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_chrg, sizeof(double)*N);
  cudaMemcpy(d_chrg, Q, sizeof(double)*N, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_pos, sizeof(double)*(3*N));
  cudaMemcpy(d_pos, X, sizeof(double)*(3*N), cudaMemcpyHostToDevice);

  //Define grid and block size
  int blocks = N / BLOCK_WIDTH;
  if(N % BLOCK_WIDTH) blocks++;

  //Allocate memory for output
  cudaMalloc((void**)&d_Ur, sizeof(double)*blocks);

  //Kernel invocation
  dim3 dimGrid(blocks, 1, 1);
  dim3 dimBlock(BLOCK_WIDTH, 1, 1);
  real_coulomb_kernel<<<dimGrid,dimBlock>>>(d_rad, d_chrg, d_pos, N, L, alpha,
					    d_Ur);

  //copy output to host
  double *Ur;
  Ur = new double[blocks];
  cudaMemcpy(Ur, d_Ur, sizeof(double)*blocks, cudaMemcpyDeviceToHost);

  double tot_Ur = 0.;
  for(int i = 0; i < blocks; ++i)
    tot_Ur += Ur[i];
  return tot_Ur;
}


/*
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
    double ri[3];
    for(int a = 0; a < 3; ++a)
      ri[a] = part.get_pos(i, a);
    double qi = part.get_charge(i);
    
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

  int N = part.get_Ntot();
  
  double V = L*L*L;
  double P2 = 2*M_PI/L;
  kmax += 1;
  int kmax2 = kmax*kmax;
  int kmax3 = kmax2*kmax;  
  
  #pragma omp parallel
  #pragma omp single
  {
    for(int kn = 0; kn < kmax3; kn++){
      
      #pragma omp task firstprivate(kn) shared(Uk)
      {
	double partUk = 0;
	
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
	  
	  partUk += recip_coulomb(part, N, kk2, K);
	
	  //Correct for the symmetries used
	  if((nx == 0 && ny == 0) || (nx == 0 && nz == 0) || (ny == 0 && nz == 0))
	    partUk /= 4;
	  else if((nz == 0 && nx != 0 && ny != 0) || (ny == 0 && nz != 0 && nx != 0)
		  || (nx == 0 && ny != 0 && nz != 0))
	    partUk /= 2;
	}
	
	#pragma omp atomic
	Uk += partUk;
      }
      
    }
  }
  
  //Self energy of the main cell
  double self = 0;
  for(int i = 0; i < N; i++){
    double qi = part.get_charge(i);
    self += qi*qi;
  }
  self *= 0.5 * alpha / sqrt(M_PI);

  Uk -= self;
  
  return Uk;
}
*/
