#include "classes.hpp"
#include <iostream>



//Particles*****************************************************************************
//**************************************************************************************
Particles::Particles(): Nkinds(1) {
  Npart = new int[1];
  m = new double[1];
  r = new double[1];
  q = new double[1];

  pos = new double[3];
  mom = new double[3];  
}

Particles::Particles(int kinds, const std::vector<int>& N, const std::vector<double>& mass,
		     const std::vector<double>& rad, const std::vector<double>& cha): Nkinds(kinds) {
  Npart = new int[Nkinds];
  
  Ntot = 0;
  for(int n = 0; n < Nkinds; ++n){
    Npart[n] = N[n];
    Ntot += N[n];
  }

  m = new double[Ntot];
  r = new double[Ntot];
  q = new double[Ntot];
  int i = 0;
  int sp = 0;
  int Nac = Npart[0];
  while(i < Ntot){
    if(i == Nac){
      sp++;
      Nac += Npart[sp];
    }
    m[i] = mass[sp];
    r[i] = rad[sp];
    q[i] = cha[sp];
    i++;
  }
  pos = new double[3*Ntot];
  mom = new double[3*Ntot];
}

Particles::~Particles() {
  delete[] Npart;
  delete[] m;
  delete[] r;
  delete[] q;
  
  delete[] pos;
  delete[] mom;
}


//Value getters---------------------------------------------------------------
int Particles::get_Nkinds() const{
  return Nkinds;
}

int Particles::get_N(int i) const{
  if(i > Nkinds-1)
    throw std::out_of_range("get_N: Species index out of bounds.");
  return Npart[i];
}

int Particles::get_Ntot() const{
  return Ntot;
}

double Particles::get_mass(int i) const{
  if(i > Ntot-1)
    throw std::out_of_range("get_mass: Particle index out of bounds.");
  return m[i];
}

double Particles::get_rad(int i) const{
  if(i > Ntot-1)
    throw std::out_of_range("get_rad: Particle index out of bounds.");
  return r[i];
}

double Particles::get_charge(int i) const{
  if(i > Ntot-1)
    throw std::out_of_range("get_charge: Particle index out of bounds.");
  return q[i];
}

//Get coordinate "a" from particle "i"; 
double Particles::get_pos(int i, int a) const{
  //check boundaries
  if(i > Ntot-1)
    throw std::out_of_range("get_pos: Particle index out of bounds.");
  if(a > 2)
    throw std::out_of_range("get_pos: Component out of bounds.");
  
  //get component "a" of the position
  return pos[i*3+a];
}

//Get momentum component "a" from particle "i"; 
double Particles::get_mom(int i, int a) const{
  //check boundaries
  if(i > Ntot-1)
    throw std::out_of_range("get_mom: Particle index out of bounds.");
  if(a > 2)
    throw std::out_of_range("get_mom: Component out of bounds.");
  
  //get component "a" of the position
  return mom[i*3+a];
}



//Array getters---------------------------------------------------------------
int* Particles::get_Npart() {
  return Npart;
}

double* Particles::get_M() {
  return m;
}

double* Particles::get_R() {
  return r;
}

double* Particles::get_Q() {
  return q;
}

double* Particles::get_X() {
  return pos;
}

double* Particles::get_P() {
  return mom;
}


//Setters-------------------------------------------------------------------------
//Set mass to particle "i"
void Particles::set_mass(int i, double mass) {
  //check boundaries
  if(i > Ntot-1)
    throw std::out_of_range("set_mass: Particle index out of bounds.");
  m[i] = mass;
}

//Set radius of particle "i"
void Particles::set_rad(int i, double rad) {
  //check boundaries
  if(i > Ntot-1)
    throw std::out_of_range("set_rad: Particle index out of bounds.");
  r[i] = rad;
}

//Set charge to particle "i"
void Particles::set_charge(int i, double cha) {
  //check boundaries
  if(i > Ntot-1)
    throw std::out_of_range("set_pos: Particle index out of bounds.");
  q[i] = cha;
}

//Set coordinate "a" from particle "i"
void Particles::set_pos(int i, int a, double val){
  //check boundaries
  if(i > Ntot-1)
    throw std::out_of_range("set_pos: Particle index out of bounds.");
  if(a > 2)
    throw std::out_of_range("set_pos: Component out of bounds.");
  
  pos[i*3+a] = val;
}

//Set momentum component "a" from particle "i" of species "n"
void Particles::set_mom(int i, int a, double val){
  //check boundaries
  if(i > Ntot-1)
    throw std::out_of_range("set_mom: Particle index out of bounds.");
  if(a > 2)
    throw std::out_of_range("set_mom: Component out of bounds.");
  
  mom[i*3+a] = val;
}




//Kvector*************************************************************************************
//********************************************************************************************

Kvector::Kvector(unsigned s) {
  _size = s;
  kvec = new double[_size];
}

Kvector::~Kvector() {
  delete[] kvec;
}

unsigned Kvector::size() const{
  return _size;
}

double Kvector::get(unsigned i) const{
  if(i >= _size)
    throw std::out_of_range("Kvector get: out of range.");
  return kvec[i];
}

void Kvector::set(unsigned i, double val) {
  if(i >= _size)
    throw std::out_of_range("Kvector set: out of range.");
  kvec[i] = val;
}
