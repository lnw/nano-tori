// Copyright (c) 2019, Lukas Wirz
// All rights reserved.

// This file is part of 'nano-tori' which is released under the BSD-2-clause license.
// See file LICENSE in this project.

#ifndef GEOMETRY2_HH
#define GEOMETRY2_HH

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "auxiliary.hh"

using namespace std;

struct matrix2d;

struct coord2d {
  double x[2];

  coord2d(const double y[2]) { x[0] = y[0]; x[1] = y[1]; }
  explicit coord2d(const double x_=0, const double y=0) { x[0] = x_; x[1] = y; }

  coord2d operator/(const double s)   const { return coord2d(*this) /= s; }
  coord2d operator*(const double s)   const { return coord2d(*this) *= s; }
  coord2d operator+(const coord2d& y) const { return coord2d(*this) += y; }
  coord2d operator-(const coord2d& y) const { return coord2d(*this) -= y; }
  coord2d& operator+=(const coord2d& y){ x[0] += y[0]; x[1] += y[1]; return *this; }
  coord2d& operator-=(const coord2d& y){ x[0] -= y[0]; x[1] -= y[1]; return *this; }
  coord2d& operator*=(const double& y){ x[0] *= y; x[1] *= y; return *this; }
  coord2d& operator/=(const double& y){ x[0] /= y; x[1] /= y; return *this; }
  coord2d operator-() const {coord2d y(-x[0],-x[1]); return y;}
  coord2d operator*(const matrix2d& m) const;


  double& operator[](unsigned int i){ return x[i]; }
  double  operator[](unsigned int i) const { return x[i]; }

  double norm() const {return sqrt(x[0]*x[0] + x[1]*x[1]);}


  friend ostream& operator<<(ostream& S, const coord2d &c2) {
    S << "{" << c2[0] << ", " << c2[1] << "}";
    return S;
  }
};

coord2d operator*(const double &d, const coord2d &c2d);


struct structure2d {
  vector<coord2d> dat;

  vector<coord2d>::iterator begin(){return dat.begin();}
  vector<coord2d>::iterator end(){return dat.end();}

  structure2d operator*(const matrix2d& m) const;
  coord2d& operator[](unsigned int i){ return dat[i]; }
  coord2d  operator[](unsigned int i) const { return dat[i]; }
  void push_back(const coord2d &c2d){dat.push_back(c2d);}
  size_t size() const {return dat.size();}

  bool contains(coord2d p) const;
  void crop(double x1, double y1, double x2, double y2);

  // assume that all points are consecutive, and both ends are connected
  // if the curve is self-intersecting, every intersection changes the sign, but the total sign is positive
  double area() const;

  friend ostream& operator<<(ostream& S, const structure2d &s2) {
    S << s2.dat;
    return S;
  }
};


struct matrix2d {
  double values[4];

  explicit matrix2d(const double w=0, const double x=0, const double y=0, const double z=0) { // (0,0), (0,1), (1,0), (1,1)
    values[0]=w; values[1]=x; values[2]=y; values[3]=z;
  }
  double& operator()(int i, int j)         { return values[i*2+j]; } // i: row, j: column
  double  operator()(int i, int j) const   { return values[i*2+j]; }
  double& operator[](unsigned int i)       { return values[i]; }
  double  operator[](unsigned int i) const { return values[i]; }

  coord2d operator*(const coord2d& c2d) const;
  structure2d operator*(const structure2d& s2d) const;

  friend ostream& operator<<(ostream& S, const matrix2d &M) {
    S << "{"; for(int i=0;i<2;i++) S << vector<double>(&M.values[i*2],&M.values[(i+1)*2]) << (i+1<2?",":"}");
    return S;
  }
};


#endif

