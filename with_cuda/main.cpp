#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include "functions.hpp"
#include "dev_functions.hpp"
#include "classes.hpp"


int main(int argc, char *argv[]) {
  int Nreps = atoi(argv[1]);

  int Nkinds = 2;
  std::vector<int> N = {50000, 50000};//{40, 4800};
  std::vector<double> m = {1.0, 1.0};//{1.0, 0.022};
  std::vector<double> r = {0.5, 0.5};//{1.0, 0.06};
  std::vector<double> q = {1, -1};//{41.52, -0.346};
  
  Particles part(Nkinds, N, m, r, q);
  double l[3];
  
  //charge configuration
  if(! chrg_conf(part, l)){
    std::cout << "problem opening config.dat for reading\n";
    return 1;
  }
  double L = l[0];
  double alpha = 1.0;
  double kmax = 20; //has to be >= 1
  double rcut_perc = 10;
  
  if(rcut_perc > 100) rcut_perc = 100;
  double rcut = rcut_perc * (L/2) / 100;
  
  std::cout << "L: " << L << "\n";
  std::cout << "alpha: " << alpha << "\n";
  std::cout << "rcut: " << rcut << "\n";
  std::cout << "kmax: " << kmax << "\n";
  // std::cout << "P0: q= " << part.get_charge(0)
  // 	    << " pos= " << part.get_pos(0, 0) << ","
  // 	    << part.get_pos(0, 1) << ","
  // 	    << part.get_pos(0, 2) << "\n";
  // std::cout << "P39: q= " << part.get_charge(39)
  // 	    << " pos= " << part.get_pos(39, 0) << ","
  // 	    << part.get_pos(39, 1) << ","
  // 	    << part.get_pos(39, 2) << "\n";
  // std::cout << "P40: q= " << part.get_charge(40)
  // 	    << " pos= " << part.get_pos(40, 0) << ","
  // 	    << part.get_pos(40, 1) << ","
  // 	    << part.get_pos(40, 2) << "\n";
  // std::cout << "P4839: q= " << part.get_charge(4839)
  // 	    << " pos= " << part.get_pos(4839, 0) << ","
  // 	    << part.get_pos(4839, 1) << ","
  // 	    << part.get_pos(4839, 2) << "\n";


  unsigned ksize = pow((kmax+1),3);
  Kvector Kvec(ksize);
  k_vector(Kvec, L, alpha, kmax);
  
  std::cout.precision(12);
  
  double tac_r = 0;
  double tsqac_r = 0;
  double tac_k = 0;
  double tsqac_k = 0;
  double Ur, Uk;
  for(int n = 0; n < Nreps; n++){
     auto start = std::chrono::high_resolution_clock::now();
     Ur = real_potential(part, L, alpha, rcut);
     auto stop = std::chrono::high_resolution_clock::now();
     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
     double t = duration.count()*1.e-3;
     tac_r += t;
     tsqac_r += t*t;

     start = std::chrono::high_resolution_clock::now();
     Uk = recip_potential(part, Kvec, L, alpha, kmax);
     stop = std::chrono::high_resolution_clock::now();
     duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
     t = duration.count()*1.e-3;
     tac_k += t;
     tsqac_k += t*t;
  
     std::cout << "Repetitions: " << n+1 << "\n";
  }

  double tav_r = tac_r/Nreps;
  double tsd_r = sqrt(tsqac_r/Nreps - tav_r*tav_r);
  std::cout << "Ur/N: " << Ur/(N[0]+N[1]) 
	    << " time: " << tav_r << " +/- " 
	    << tsd_r << " ms.\n";
  
  double tav_k = tac_k/Nreps;
  double tsd_k = sqrt(tsqac_k/Nreps - tav_k*tav_k);
  std::cout << "Uk/N: " << Uk/(N[0]+N[1]) 
	    << " time: " << tav_k << " +/- " 
	    << tsd_k << " ms.\n";
    
  std::cout << "Total energy per particle: " << (Ur+Uk)/(N[0]+N[1]) << "\n";
  
  return 0;
}
