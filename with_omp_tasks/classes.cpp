#include "classes.hpp"
#include <iostream>

//NeighborCells*******************************************************************************
//********************************************************************************************
NeighborCells::NeighborCells(unsigned sz, double l) {
  side = sz;
  surf = side*side;
  n_cells = surf*side;
  n_neighbors = 26;
  L = l;
  len = L/side;
  
  array = new int*[n_cells];
  for(int i = 0; i < n_cells; i++)
    array[i] = new int[n_neighbors];

  for(int i = 0; i < n_cells; i++){
    unsigned lay = floor(i/surf);
    unsigned row = floor((i-lay*surf)/side);
    unsigned col = i - lay*surf - row*side;
    
    array[i][0] = subone(lay)*surf + subone(row)*side + subone(col);
    array[i][1] = subone(lay)*surf + subone(row)*side + col;
    array[i][2] = subone(lay)*surf + subone(row)*side + addone(col);
    array[i][3] = subone(lay)*surf + row*side + subone(col);
    array[i][4] = subone(lay)*surf + row*side + col;
    array[i][5] = subone(lay)*surf + row*side + addone(col);
    array[i][6] = subone(lay)*surf + addone(row)*side + subone(col);
    array[i][7] = subone(lay)*surf + addone(row)*side + col;
    array[i][8] = subone(lay)*surf + addone(row)*side + addone(col);

    array[i][9] = lay*surf + subone(row)*side + subone(col);
    array[i][10] = lay*surf + subone(row)*side + col;
    array[i][11] = lay*surf + subone(row)*side + addone(col);
    array[i][12] = lay*surf + row*side + subone(col);
    array[i][13] = lay*surf + row*side + addone(col);
    array[i][14] = lay*surf + addone(row)*side + subone(col);
    array[i][15] = lay*surf + addone(row)*side + col;
    array[i][16] = lay*surf + addone(row)*side + addone(col);
    
    array[i][17] = addone(lay)*surf + subone(row)*side + subone(col);
    array[i][18] = addone(lay)*surf + subone(row)*side + col;
    array[i][19] = addone(lay)*surf + subone(row)*side + addone(col);
    array[i][20] = addone(lay)*surf + row*side + subone(col);
    array[i][21] = addone(lay)*surf + row*side + col;
    array[i][22] = addone(lay)*surf + row*side + addone(col);
    array[i][23] = addone(lay)*surf + addone(row)*side + subone(col);
    array[i][24] = addone(lay)*surf + addone(row)*side + col;
    array[i][25] = addone(lay)*surf + addone(row)*side + addone(col);
  }
}

NeighborCells::~NeighborCells() {
  for(int i = 0; i < n_cells; i++)
    delete[] array[i];
  delete[] array;
}

unsigned NeighborCells::addone(int num) {
  return static_cast<unsigned>( (num+1)%side );
}

unsigned NeighborCells::subone(int num) {
  int res = num-1;
  if(res < 0)
    return static_cast<unsigned>( side-1 );
  else
    return static_cast<unsigned>( res%side );
}

double NeighborCells::get_len() const{
  return len;
}

bool NeighborCells::find(unsigned i, unsigned val) const{
  for(int j = 0; j < n_neighbors; j++){
    //std::cout << array[i][j] << "\n";
    if(array[i][j] == val)
      return true;
  }
  return false;
}

unsigned NeighborCells::which_cell(double x, double y, double z) const{
  double L2 = L/2;
  // z -= L*floor(z/L+0.5);
  // y -= L*floor(y/L+0.5);
  // x -= L*floor(x/L+0.5);
  unsigned lay = floor((z+L2)/len);
  unsigned row = floor((y+L2)/len);
  unsigned col = floor((x+L2)/len);
  return lay*surf + row*side + col; 
}

void NeighborCells::display() const{
  for(int i = 0;  i < n_cells; i++){
    std::cout << i << ": ";
    for(int j = 0; j < n_neighbors; j++)
      std::cout << array[i][j] << " ";
    std::cout << "\n";
  }
}




//Particles*****************************************************************************
//**************************************************************************************
Particles::Particles(): Nkinds(1) {
  Npart = new int[1];
  m = new double[1];
  r = new double[1];
  q = new double[1];

  pos = new double[3];
  mom = new double[3];
  
  cells = new unsigned[1];
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
  cells = new unsigned[Ntot];
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
    cells[i] = 0;
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

  delete[] cells;
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

//Get cell in which the particle is located
unsigned Particles::get_cell(int i) const{
  return cells[i];
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

unsigned* Particles::get_C() {
  return cells;
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

//Set cell in which the particle is located
void Particles::set_cells(const NeighborCells &ncells) {
  for(int i = 0; i < Ntot; i++)
    cells[i] = ncells.which_cell(pos[i*3+0], pos[i*3+1], pos[i*3+2]);
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
