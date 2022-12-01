#include <iostream>
#include <chrono>
#include "functions.hpp"
#include "classes.hpp"


int main() {
  int Nkinds = 2;
  std::vector<int> N = {40, 4800};
  std::vector<double> m = {1.0, 0.022};
  std::vector<double> r = {1.0, 0.06};
  std::vector<double> q = {41.52, -0.346};
  
  Particles part(Nkinds, N, m, r, q);
  double l[3];
  
  //charge configuration
  if(! chrg_conf(part, l)){
    std::cout << "problem opening conf.dat for reading\n";
    return 1;
  }
  double L = l[0];
  double alpha = 0.28;
  double kmax = 15; //has to be >= 1
  
  std::cout << "L: " << L << "\n";
  std::cout << "P0: q= " << part.get_charge(0)
	    << " pos= " << part.get_pos(0, 0, 0) << ","
	    << part.get_pos(0, 0, 1) << ","
	    << part.get_pos(0, 0, 2) << "\n";
  std::cout << "P39: q= " << part.get_charge(0)
	    << " pos= " << part.get_pos(0, 39, 0) << ","
	    << part.get_pos(0, 39, 1) << ","
	    << part.get_pos(0, 39, 2) << "\n";
  std::cout << "P40: q= " << part.get_charge(1)
	    << " pos= " << part.get_pos(1, 0, 0) << ","
	    << part.get_pos(1, 0, 1) << ","
	    << part.get_pos(1, 0, 2) << "\n";
  std::cout << "P4839: q= " << part.get_charge(1)
	    << " pos= " << part.get_pos(1, 4799, 0) << ","
	    << part.get_pos(1, 4799, 1) << ","
	    << part.get_pos(1, 4799, 2) << "\n";


  unsigned ksize = pow((kmax+1),3);
  Kvector Kvec(ksize);
  k_vector(Kvec, L, alpha, kmax);
  
  
  auto start = std::chrono::high_resolution_clock::now();
  double Ur = real_potential(part, L, alpha);
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
  std::cout << "Ur: " << Ur << " time: " << duration.count()*1.e-3 << " ms.\n";

  start = std::chrono::high_resolution_clock::now();
  double Uk = recip_potential(part, Kvec, L, alpha, kmax);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
  std::cout << "Uk: " << Uk << " time: " << duration.count()*1.e-3 << " ms.\n";
  
  return 0;
}
