#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cstdio>

__global__ void comp_kernel() {
  cuDoubleComplex z1 = make_cuDoubleComplex(3, -5);
  cuDoubleComplex z2 = make_cuDoubleComplex(1, 1);

  cuDoubleComplex p = make_cuDoubleComplex(3*cos(2.0), 3*sin(2.0));
  printf("%f,%f\n",cuCreal(p),cuCimag(p));
  
  double n = cuCreal( cuCmul(z1, cuConj(z1)) );
  printf("%f\n",n);

  z1 = cuCadd(z1, z2);
  printf("%f,%f\n",cuCreal(z1),cuCimag(z1));
}

int main() {
  comp_kernel<<<1,1>>>();
  cudaDeviceSynchronize();
  
  return 0;
}
