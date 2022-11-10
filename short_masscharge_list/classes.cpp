#include "classes.hpp"

//Particles*****************************************************************************
//**************************************************************************************
Particles::Particles(): Nkinds(1) {
  Npart = new int[1];
  M = new double[1];
  R = new double[1];
  Q = new double[1];

  pos = new double[3];
  mom = new double[3];
}

Particles::Particles(int kinds, const std::vector<int>& N, const std::vector<double>& m,
		     const std::vector<double>& r, const std::vector<double>& q): Nkinds(kinds) {
  Npart = new int[Nkinds];
  M = new double[Nkinds];
  R = new double[Nkinds];  
  Q = new double[Nkinds];
  
  Ntot = 0;
  for(int i = 0; i < Nkinds; ++i){
    Npart[i] = N[i];
    M[i] = m[i];  
    R[i] = r[i];
    Q[i] = q[i];
    Ntot += N[i];
  }

  pos = new double[3*Ntot];
  mom = new double[3*Ntot];
}

Particles::~Particles() {
  delete[] Npart;
  delete[] M;
  delete[] R;
  delete[] Q;
  
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

double Particles::get_mass(int i) const{
  if(i > Nkinds-1)
    throw std::out_of_range("get_mass: Species index out of bounds.");
  return M[i];
}

double Particles::get_rad(int i) const{
  if(i > Nkinds-1)
    throw std::out_of_range("get_rad: Species index out of bounds.");
  return R[i];
}

double Particles::get_charge(int i) const{
  if(i > Nkinds-1)
    throw std::out_of_range("get_charge: Species index out of bounds.");
  return Q[i];
}

//Get coordinate "a" from particle "i" of species "n"; 
double Particles::get_pos(int n, int i, int a) const{
  //check boundaries
  if(i > Npart[n]-1)
    throw std::out_of_range("get_pos: Particle index out of bounds.");
  if(a > 2)
    throw std::out_of_range("get_pos: Component out of bounds.");

  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  //get component "a" of the position
  return pos[part*3+a];
}

//Get momentum component "a" from particle "i" of species "n"; 
double Particles::get_mom(int n, int i, int a) const{
  //check boundaries
  if(i > Npart[n]-1)
    throw std::out_of_range("get_mom: Particle index out of bounds.");
  if(a > 2)
    throw std::out_of_range("get_mom: Component out of bounds.");
    
  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  //get component "a" of the position
  return mom[part*3+a];
}

//Get the kind and number of the particle
void Particles::get_kind(int id, int& n, int& i) const{
  int Nto=0, Ntn=0;
  for(int nn = 0; nn < Nkinds; ++nn){
    Nto = Ntn;
    Ntn += Npart[nn];
    if(id < Ntn){
      n = nn;
      i = id - Nto;
      break;
    }
  }
}


//Array getters---------------------------------------------------------------
int* Particles::get_Npart() {
  return Npart;
}

double* Particles::get_M() {
  return M;
}

double* Particles::get_R() {
  return R;
}

double* Particles::get_Q() {
  return Q;
}

double* Particles::get_X() {
  return pos;
}

double* Particles::get_P() {
  return mom;
}


//Setters-------------------------------------------------------------------------
void Particles::set_charge(int i, double q) {
  if(i > Nkinds-1)
    throw std::out_of_range("set_charge: Species index out of bounds.");
  Q[i] = q;
}

//Set coordinate "a" from particle "i" of species "n"
void Particles::set_pos(int n, int i, int a, double val){
  //check boundaries
  if(i > Npart[n]-1)
    throw std::out_of_range("set_pos: Particle index out of bounds.");
  if(a > 2)
    throw std::out_of_range("set_pos: Component out of bounds.");
  
  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  pos[part*3+a] = val;
}

//Set momentum component "a" from particle "i" of species "n"
void Particles::set_mom(int n, int i, int a, double val){
  //check boundaries
  if(i > Npart[n]-1)
    throw std::out_of_range("set_mom: Particle index out of bounds.");
  if(a > 2)
    throw std::out_of_range("set_mom: Component out of bounds.");
    
  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  mom[part*3+a] = val;
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
