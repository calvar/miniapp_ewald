#include <cstdio>
#include <complex>
#include <omp.h>
#include <cmath>

int main() {
  //#pragma omp target
  //{
  double q = abs(-1.0);
  double shift = q < 0 ? M_PI : 0.;
  std::complex<double> c = std::polar(q,1.+shift);//(q, 1.+shift);
  
  
  printf("%f %f\n",c.real(),c.imag());
  //}
  return 0;
}
