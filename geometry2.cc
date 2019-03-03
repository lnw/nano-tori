// Copyright (c) 2019, Lukas Wirz
// All rights reserved.

// This file is part of 'nano-tori' which is released under the BSD-2-clause license.
// See file LICENSE in this project.


#include <cmath>
#include <iomanip>

#include "geometry2.hh"


coord2d operator*(const double &d, const coord2d &c2d) {return coord2d(d*c2d.x[0], d*c2d.x[1]);}

// signed area of triangle
// cross product of two edges of the triangle, in 3D, z component of result is the area
double signed_area( double x1, double y1, double x2, double y2, double x3, double y3 ){
  return 0.5 * ((x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3));
}

// vec*mat
coord2d coord2d::operator*(const matrix2d& m) const {
  return coord2d(x[0]*m(0,0)+x[1]*m(1,0),  x[0]*m(0,1)+x[1]*m(1,1));
}

structure2d structure2d::operator*(const matrix2d& m) const{
  structure2d s2d(*this);
  for(coord2d &c2d: s2d){
    c2d = c2d*m;
  }
  return s2d;
}

void structure2d::crop(double x1, double y1, double x2, double y2){
  const double eps = 1.e-8;
  assert(x1 < x2+eps);
  assert(y1 < y2+eps);
  vector<coord2d> tmp;
  for(size_t k=0; k<dat.size(); k++){
    coord2d c2d = dat[k];
    if(c2d[0] >= x1 && c2d[0] <= x2 && c2d[1] >= y1 && c2d[1] <= y2){
      tmp.push_back(c2d);
    }
  }
  dat = tmp;
}

bool structure2d::contains(coord2d p) const {
  const double eps = 1.e-5;
  for(const coord2d &c2d: dat){
    if(abs(c2d[0] - p[0]) < eps && abs(c2d[1] - p[1]) < eps) return true;
  }
  return false;
}


double structure2d::area() const {
  const coord2d p(0,0); // reference point, position is irrelevant
  double A=0;
  for(size_t i=0; i<dat.size(); i++){
    const int ii = (i+1)%dat.size(); // the next element
    A += signed_area(p[0], p[1], dat[i][0], dat[i][1], dat[ii][0], dat[ii][1]);
  }
  return abs(A);
}

// mat*vec
coord2d matrix2d::operator*(const coord2d& x) const {
  const matrix2d &m = *this;
  return coord2d( m(0,0)*x[0]+m(0,1)*x[1],  m(1,0)*x[0]+m(1,1)*x[1] );
}


structure2d matrix2d::operator*(const structure2d& s2d) const{
  const matrix2d &m = *this;
  structure2d s2d_new(s2d);
  for(coord2d &c2d: s2d_new){
    c2d = m*c2d;
  }
  return s2d_new;
}


