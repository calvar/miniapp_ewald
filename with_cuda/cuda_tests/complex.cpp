#include <complex>
#include <cstdio>

int main() {
  std::complex<double> z1(3,-5);
  std::complex<double> z2(1,1);

  std::complex<double> p = std::polar(3.0, 2.0);
  printf("%f,%f\n",p.real(),p.imag());
  
  double n = std::norm(z1);
  printf("%f\n",n);

  z1 += z2;
  printf("%f,%f\n",z1.real(),z1.imag());
  
  return 0;
}
