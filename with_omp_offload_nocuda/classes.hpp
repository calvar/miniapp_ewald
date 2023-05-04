#ifndef __CLASSES_HPP
#define __CLASSES_HPP

#include <vector>
#include <stdexcept>
#include <cmath>

//NeighborCells. The real space is divided into cells. Particles interact with
// others in their cell and neighboring cells.
class NeighborCells {
private:
  int **array;
  unsigned side;
  unsigned surf;
  unsigned n_cells;
  unsigned n_neighbors;
  double L; //simulation box size
  double len; //lenght of the cell

  unsigned addone(int num);
  unsigned subone(int num);
  
public:
  NeighborCells(unsigned sz, double l);
  ~NeighborCells();

  double get_len() const;
  bool find(unsigned i, unsigned j) const;
  unsigned which_cell(double x, double y, double z) const;
  void display() const;
};


//Class that stores particle positions and properties in 
// 1D dynamicaly alocated arays, so that copying them
// to the device is easy.************************************************
class Particles {
private:
  int Nkinds; //How many kinds of particles are there
  int* Npart; //How many particles of each kind are there
  int Ntot;   //Total number of particles
  double* m; //Mass of each particle (linearized 3*N array)
  double* r; //Radius of each particle (linearized 3*N array)
  double* q; //Charge of each particle (linearized 3*N array)

  double* pos; //Position of each particle (linearized 3*N array)
  double* mom; //Momentum of each particle (linearized 3*N array)

  unsigned* cells; //Cells in which each particle is located 
  
public:
  Particles();
  Particles(int, const std::vector<int>&, const std::vector<double>&,
	    const std::vector<double>&, const std::vector<double>&);
  ~Particles();

  //Value getters
  int get_Nkinds() const;
  int get_N(int) const;
  int get_Ntot() const;
  double get_mass(int) const;
  double get_rad(int) const;
  double get_charge(int) const;
  double get_pos(int, int) const;
  double get_mom(int, int) const;
  unsigned get_cell(int) const;
   
  //Array getters (these give access to private members! Think about friend functions)
  int* get_Npart();
  double* get_M();
  double* get_R();
  double* get_Q();
  double* get_X();
  double* get_P();
  unsigned* get_C();
  
  //Setters
  void set_mass(int, double);
  void set_rad(int, double);
  void set_charge(int, double);
  void set_pos(int, int, double);
  void set_mom(int, int, double);
  void set_cells(const NeighborCells&);
};


//Kvector
class Kvector {
private:
  double *kvec;
  unsigned _size;
  
public:
  Kvector(unsigned s);
  ~Kvector();

  unsigned size() const;
  double * get_all() {return kvec ;}
  double get(unsigned i) const;
  void set(unsigned i, double val);
};

#endif
