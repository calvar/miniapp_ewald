#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>

int main() {
  double dens = 0.3;
  int Nkinds = 2;
  std::vector<int> N = {50000, 50000};
  std::vector<double> m = {1.0, 1.0};
  std::vector<double> r = {0.5, 0.5};
  std::vector<double> q = {1.0, -1.0};

  double part_vol = 0;
  for(int n = 0; n < N.size(); n++)
    part_vol += N[n]*r[n]*r[n];
  part_vol *= 4*M_PI/3;
  
  double L = pow(part_vol/dens, 1./3);

  
  std::default_random_engine rng;
  std::uniform_real_distribution<double> dist(-L/2,L/2);

  int Ntot = 0;
  for(int i = 0; i < N.size(); i++){
    Ntot += N[i];
  }

  std::vector< std::vector<double> > params_vec;
  for(int n = 0; n < N.size(); n++){
    for(int i = 0; i < N[n]; i++){
      std::vector<double> v = {m[n], r[n], q[n]};
      params_vec.push_back(v);
    }
  }
  
  std::vector< std::vector<double> > part_vec;
  for(int i = 0; i < Ntot; i++){
    //std::cout << i << "\n";
    
    double pos[3];
    bool test = true;
    while(test){
      for(int a = 0; a < 3; a++){
	pos[a] = dist(rng);
      }
      double radi = params_vec[i][1]; 
      
      test = false;
      for(int j = 0; j < part_vec.size(); j++){
	double radj = part_vec[j][1]; 
	
	double rijsq = 0;
	for(int a = 0; a < 3; a++){
	  double x = part_vec[j][3+a] - pos[a];
	  x -= L * floor(x/L + 0.5);
	  rijsq += x*x;
	}
	if(rijsq < (radi+radj)*(radi+radj)){
	  test = true;
	  break;
	}
      }
    }
    
    std::vector<double> v = params_vec[i];
    for(int a = 0; a < 3; a++){
      v.push_back(pos[a]);
    }
    part_vec.push_back(v);
  }

  // for(int i = 0; i < part_vec.size(); i++){
  //   for(int b = 0; b < part_vec[0].size(); b++)
  //     std::cout << part_vec[i][b] << " ";
  //   std::cout << "\n";
  // }
  
  std::vector< std::vector<double> > new_part_vec;
  std::ofstream Out("config.dat");
  for(int n = 0; n < q.size(); n++)
    Out << q[n] << " ";
  Out << "\n";
  Out << L << " " << L << " " << L << "\n"; 
  for(int i = 0; i < part_vec.size(); i++){
    for(int a = 0; a < 3; a++)
      Out << part_vec[i][3+a] << " ";
    Out << 0 << " " << 0 << " " << 0 << "\n";
  }
  Out.close();
  
  return 0;
}
