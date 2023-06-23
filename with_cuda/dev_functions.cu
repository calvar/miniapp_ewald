#include <cuda_runtime.h>
#include <cuComplex.h>
//#include <cuda_profiler_api.h>
#include <cstdio>
#include <iostream>
#include "dev_functions.hpp"

const int BLOCK_WIDTH_R = 512; //Max threads per block
const int BLOCK_WIDTH_K = 128;  //Max threads per block
const int K_SIM_SIZE = 4;


//Implementation from the CUDA page: https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#atomic-functions
#if defined(__CUDA_ARCH__) && __CUDA_ARCH < 600
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
#endif

//Atomic add for cuDoubleComplex
__device__ void atAddComplex(cuDoubleComplex* a, cuDoubleComplex b){
  //transform the addresses of real and imag. parts to double pointers
  double *x = (double*)a;
  double *y = x+1; //add 1*sizeof(double) bits

  //use atomicAdd for double variables
  atomicAdd(x, cuCreal(b));
  atomicAdd(y, cuCimag(b));
}



//Coulomb potential using Ewald summation*********************************************


//REAL PART-----------------------------------------------------------------
__global__ void real_coulomb_kernel(double *rad, double *chrg, double *pos,
				    int N, double L, double alpha, double *Ur,
				    double rcut/*,int *count*/) {

  int i = threadIdx.x + blockDim.x * blockIdx.x;

  //int lcount = 0;
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
      }
      j += 1 - N*static_cast<int>(floor((j+1)/N + 0.5));
      cnt++;
    }
    //lcount += cnt;
    Ur[i] = U;
  }
  
}

double real_potential(const Particles &part, double L, double alpha,
		      double rcut) {
  int N = part.get_Ntot();

  double *d_rad, *d_chrg, *d_pos, *d_Ur;
  //int *d_count;

  //Get particle data arrays
  double* R = part.get_R();
  double* Q = part.get_Q();
  double* X = part.get_X();

  //Set service
  //cudaService(0);

  //Allocate and copy memory to device
  cudaMalloc((void**)&d_rad, sizeof(double)*N);
  cudaMemcpy(d_rad, R, sizeof(double)*N, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_chrg, sizeof(double)*N);
  cudaMemcpy(d_chrg, Q, sizeof(double)*N, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_pos, sizeof(double)*(3*N));
  cudaMemcpy(d_pos, X, sizeof(double)*(3*N), cudaMemcpyHostToDevice);

  //Define grid and block size
  int blocks = N / BLOCK_WIDTH_R;
  if(N % BLOCK_WIDTH_R) blocks++;

  //Allocate memory for output
  //cudaMalloc((void**)&d_count, sizeof(int)*N);
  cudaMalloc((void**)&d_Ur, sizeof(double)*N);
  
  //Kernel invocation
  dim3 dimGrid(blocks, 1, 1);
  dim3 dimBlock(BLOCK_WIDTH_R, 1, 1);
  real_coulomb_kernel<<<dimGrid,dimBlock>>>(d_rad, d_chrg, d_pos, N, L, alpha,
					    d_Ur, rcut/*, d_count*/);

  // //copy output to host
  // int *count;
  // count = new int[N];
  // cudaMemcpy(count, d_count, sizeof(int)*N, cudaMemcpyDeviceToHost);

  // int tot_count = 0;
  // for(int i = 0; i < N; ++i)
  //   tot_count += count[i];
  // printf("No. interactions: %d\n", tot_count);
  
  //copy output to host
  double *Ur;
  Ur = new double[N];
  cudaMemcpy(Ur, d_Ur, sizeof(double)*N, cudaMemcpyDeviceToHost);

  cudaFree(d_rad);
  cudaFree(d_chrg);
  cudaFree(d_pos);
  cudaFree(d_Ur);
  
  double tot_Ur = 0.;
  for(int i = 0; i < N; ++i)
    tot_Ur += Ur[i];
  
  delete[] Ur;
  
  return tot_Ur;
}


//RECIPROCAL PART --------------------------------------------------------

__device__ double recip_coulomb(double *chrg, double *pos, int N, double kk2, double K[][3]) {
  double U = 0;
  
  //To store the summations over the particles
  cuDoubleComplex cc = make_cuDoubleComplex(0., 0.);
  cuDoubleComplex sum[4];
  for(int b = 0; b < 4; b++)
    sum[b] = cc;

  //Summation over particles
  for(int i = 0; i < N; i++){
    double posi[3];
    for(int a = 0; a < 3; ++a)
      posi[a] = pos[i*3+a];
    double qi = chrg[i];
    
    for(int b = 0; b < 4; b++){
      double posik = posi[0]*K[b][0] + posi[1]*K[b][1] + posi[2]*K[b][2]; 
      sum[b] = cuCadd( sum[b], make_cuDoubleComplex(qi*cos(posik), qi*sin(posik)) );
    }
  }
   
  //Compute energy
  for(int b = 0; b < 4; b++){
    U += kk2 * 0.5 * cuCreal( cuCmul(sum[b], cuConj(sum[b])) );
  }
  
  return U;
}

__global__ void recip_coulomb_kernel(double *chrg, double *pos, int N, int kmax,
				     const double *Kvec, double L, double alpha,
				     double *Uk/*, int *cnt*/) {
  int kn = threadIdx.x + blockDim.x * blockIdx.x;

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

double recip_potential(const Particles &part, const Kvector &Kvec, double L,
		       double alpha, int kmax) {
  int N = part.get_Ntot();

  //int *d_cnt;
  double *d_chrg, *d_pos, *d_kvec, *d_Uk;

  //Get particle data arrays
  double* Q = part.get_Q();
  double* X = part.get_X();
  
  //Get Kvec array
  int ksz = Kvec.size();
  double *kvec = Kvec.get_all();
  
  //Set service
  //cudaService(0);

  //Allocate and copy memory to device
  cudaMalloc((void**)&d_chrg, sizeof(double)*N);
  cudaMemcpy(d_chrg, Q, sizeof(double)*N, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_pos, sizeof(double)*(3*N));
  cudaMemcpy(d_pos, X, sizeof(double)*(3*N), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_kvec, sizeof(double)*(ksz));
  cudaMemcpy(d_kvec, kvec, sizeof(double)*(ksz), cudaMemcpyHostToDevice);

  //Define grid and block size
  kmax += 1;
  int kmax2 = kmax*kmax;
  int kmax3 = kmax2*kmax;
  int blocks = kmax3 / BLOCK_WIDTH_K;
  if(kmax3 % BLOCK_WIDTH_K) blocks++;

  //Allocate memory for output
  //cudaMalloc((void**)&d_cnt, sizeof(int)*kmax3);
  cudaMalloc((void**)&d_Uk, sizeof(double)*kmax3);
  
  //Kernel invocation
  dim3 dimGrid(blocks, 1, 1);
  dim3 dimBlock(BLOCK_WIDTH_K, 1, 1);
  recip_coulomb_kernel<<<dimGrid,dimBlock>>>(d_chrg, d_pos, N, kmax, d_kvec, L,
					     alpha, d_Uk/*, d_cnt*/);
  cudaDeviceSynchronize();
  cudaError_t err{cudaGetLastError()};
  if (err != cudaSuccess){
    std::cerr << cudaGetErrorString(err) << std::endl;
  }

  // //copy output to host
  // int *cnt;
  // cnt = new int[blocks];
  // cudaMemcpy(cnt, d_cnt, sizeof(int)*kmax, cudaMemcpyDeviceToHost);
  // cudaFree(d_cnt);
  
  //copy output to host
  double *Uk;
  Uk = new double[kmax3];
  cudaMemcpy(Uk, d_Uk, sizeof(double)*kmax3, cudaMemcpyDeviceToHost);
  
  cudaFree(d_chrg);
  cudaFree(d_pos);
  cudaFree(d_kvec);
  cudaFree(d_Uk);
  
  double tot_Uk = 0.;
  for(int i = 0; i < kmax3; ++i)
    tot_Uk += Uk[i];
  
  delete[] Uk;

  // int count = 0;
  // for(int i = 0; i < kmax3; ++i)
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

