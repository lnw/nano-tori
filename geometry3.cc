// Copyright (c) 2019, Lukas Wirz
// All rights reserved.

// This file is part of 'nano-tori' which is released under the BSD-2-clause license.
// See file LICENSE in this project.


#include <iomanip>
#include <limits>
#include <cmath>

#include "geometry2.hh"
#include "geometry3.hh"


void coord3d::scale(const double f){
  for(int i=0; i<3; i++){
    x[i] *= f;
  }
}

string structure3d::to_turbomole() const {
  const double aa2bohr = 1.889716164632;
  ostringstream s;
  s << setprecision(8) << fixed;
  s << "$coord" << endl;
  for(size_t i=0; i<size(); ++i){
    s << setw(12) << (*this)[i][0]*aa2bohr << "  "<< setw(12) << (*this)[i][1]*aa2bohr << "  " << setw(12) << (*this)[i][2]*aa2bohr << "  " << atomtypes[i] << endl;
  }
  s << "$end" << endl;

  return s.str();
}

string structure3d::to_xyz() const {
  ostringstream s;
  s << setprecision(6) << fixed;
  s << size() << endl;
  s << "i could write something here" << endl;
  for(size_t i=0; i<size(); ++i){
    s << atomtypes[i] << setw(12) << (*this)[i][0] << "  " << setw(10) << (*this)[i][1] << "  " << setw(10) << (*this)[i][2] << endl;
  }
  return s.str();
}

structure3d structure3d::from_xyz(const string& filename)
{
  ifstream file(filename.c_str());
  unsigned int N;
  string Nstring, comment, element,line;
  structure3d s;

  getline(file, Nstring);
  getline(file, comment);

  N = stoi(Nstring);

  //  cout << "N = " << Nstring << "; comment = " << comment << endl;

  for(size_t i=0; i < N && getline(file,line); i++){
    stringstream l(line);
    coord3d x;

    l >> element;
    for(int j=0;j<3 && l.good(); j++){
      l >> x[j];
    }

    s.push_back(x);
    s.atomtypes.push_back(element);
    //    cout << i << ": " << x << endl;
  }
  file.close();

  assert(s.size() == N);

  return s;
}

void structure3d::scale(const double f){
  for (coord3d &c3d: dat) c3d.scale(f);
}

void structure3d::invert(){
  scale(-1);
}

void structure3d::move(const coord3d m){
  for (coord3d &c3d: dat) c3d+=m;
}

void structure3d::centre_at_origin(){
  coord3d barycentre({0,0,0});
  for (coord3d c3d: dat) barycentre += c3d;
  barycentre /= size();
  move(-barycentre);
}

void structure3d::clear(){
  dat.clear();
  atomtypes.clear();
}

void structure3d::resize(const size_t n){
  dat.resize(n);
  atomtypes.resize(n);
}

void structure3d::erase(const size_t n){
  dat.erase(dat.begin() + n);
  atomtypes.erase(atomtypes.begin() + n);
}

matrix3d coord3d::outer(const coord3d& y) const {
  return matrix3d(x[0]*y[0],x[0]*y[1],x[0]*y[2],  x[1]*y[0],x[1]*y[1],x[1]*y[2],  x[2]*y[0],x[2]*y[1],x[2]*y[2]);
}

void structure3d::push_back(const structure3d& s3d){
  for(coord3d c3d: s3d.dat) push_back(c3d);
  for(string a: s3d.atomtypes) push_back_atom(a);
}

