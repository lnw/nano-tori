#ifndef SHAPEGEN_HH
#define SHAPEGEN_HH

#include <cassert>
#include <cmath>
#include <fstream>

#include "auxiliary.hh"
#include "geometry3.hh"

using namespace std;


auto tube = +[] (const double t) -> coord3d { // convert lambda to function pointer
  const double x = t;
  const double y = 0;
  const double z = 0;
  return coord3d(x,y,z);
};

auto tube_tangent = +[](const double t) -> coord3d {
  const double x = 1.0;
  const double y = 0;
  const double z = 0;
  return coord3d(x,y,z).normalised();
};

auto circle = +[] (const double t) -> coord3d {
  const double x = cos(t);
  const double y = sin(t);
  const double z = 0;
  return coord3d(x,y,z);
};

auto circle_tangent = +[](const double t) -> coord3d {
  const double x = -sin(t);
  const double y = cos(t);
  const double z = 0;
  return coord3d(x,y,z).normalised();
};

auto trefoil = +[](const double t) -> coord3d {
  const double x = sin(t) + 2*sin(2*t);
  const double y = cos(t) - 2*cos(2*t);
  const double z = -sin(3*t);
  return coord3d(x,y,z);
};

auto trefoil_tangent = +[](const double t){
  const double x = cos(t) + 4*cos(2*t);
  const double y = -sin(t) + 4*sin(2*t);
  const double z = -3*cos(3*t);
  return coord3d(x,y,z).normalised();
};

auto lemniscate = +[](const double t) -> coord3d {
  const double a = 2.0;
  const double b = 0.3;
  const double x = a*cos(t)/(1+sin(t)*sin(t));
  const double y = a*sin(t)*cos(t)/(1+sin(t)*sin(t));
  const double z = b * sin(t);
  return coord3d(x,y,z);
};

auto lemniscate_tangent = +[](const double t){
  const double a = 2.0;
  const double b = 0.3;
  const double x = -(2*a*(5+cos(2*t))*sin(t))/pow(-3+cos(2*t),2);
  const double y = (2*a*(-1+3*cos(2*t)))/pow(-3+cos(2*t),2);
  const double z = b * cos(t);
  return coord3d(x,y,z);
};

auto curve35 = +[](const double t) -> coord3d {
  const double x = cos(3*t) * (1.0 + 0.4*cos(5*t));
  const double y = sin(3*t) * (1.0 + 0.4*cos(5*t));
  const double z = 0.4*sin(5*t);
  return coord3d(x,y,z);
};

auto curve35_tangent = +[](const double t){
  const double x = -3 *(1 + 0.4 *cos(5*t)) * sin(3*t) - 2.*cos(3*t)* sin(5*t);
  const double y = 3*cos(3*t)*(1 + 0.4*cos(5*t)) - 2.*sin(3*t)*sin(5*t);
  const double z = 1.5*cos(5*t);
  return coord3d(x,y,z).normalised();
};



class ShapeGenerator {

public:

  auto get_shape(int c) const {
    switch(c){
    case 1: return tube;
    case 2: return circle;
    case 3: return trefoil;
    case 4: return lemniscate;
    case 5: return curve35;
    default:
      cerr << "invalid curve chosen, aborting ..." << endl;
      break;
    }
  }

  auto get_shape_tan(int c) const {
    switch(c){
    case 1: return tube_tangent;
    case 2: return circle_tangent;
    case 3: return trefoil_tangent;
    case 4: return lemniscate_tangent;
    case 5: return curve35_tangent;
    default:
      cerr << "invalid curve chosen, aborting ..." << endl;
      break;
    }
  }

  bool is_closed_curve(int c) const {
    switch(c){
    case 1: return false;
    case 2:
    case 3:
    case 4:
    case 5: return true;
    default:
      cerr << "invalid curve chosen, aborting ..." << endl;
      break;
    }
  }

  vector<double> get_inverse_mapping(int c, int n_steps, double& length_tot) const;

};



// put points on a line along the curve, equidistant in t.  Measure distances
// along neighbouring points.  Use total length to scale, and individual
// lengths to space equidistantly.
  vector<double> ShapeGenerator::get_inverse_mapping(int c, int n_steps, double& length_tot) const {

  auto shape = get_shape(c);

  const int n_points = n_steps + 1;

  structure3d S_aux; // auxiliary structure following the line shape
  for (int i = 0; i<n_points; i++){
    const double ii = 2*M_PI*i / n_steps;
    const coord3d c3d = shape(ii);
    S_aux.push_back(c3d);
    S_aux.push_back_atom(string("C"));
  }
  // cout << "saux: " << S_aux << endl;

  // distances
  vector<double> dD_dTheta(n_points);
  for(int i = 0; i<n_points; i++){
    //const coord3d previous = S_aux[i-1];
    const coord3d here = S_aux[i];
    // const coord3d next = S_aux[i+1];
    if(i==0)
       dD_dTheta[i] = (here-S_aux[i+1]).norm();
    else if(i==n_points-1)
       dD_dTheta[i] = (here-S_aux[i-1]).norm();
    else
       dD_dTheta[i] = ((here-S_aux[i+1]).norm() + (here-S_aux[i-1]).norm())/2;
  }
  // cout << "dD_dTheta : " << dD_dTheta << endl;

  // distances inverse, total length
  length_tot = 0;
  for(const double d: dD_dTheta){
    length_tot += d;
  }
  length_tot -= dD_dTheta[0]/2; // because there are more points than intervals
  length_tot -= dD_dTheta[n_points-1]/2;
  // cout << "tot: " << length_tot << endl;

  double d_avg = length_tot / n_steps; // average derivative
  // cout << "avg: " << d_avg << endl;

  // distances inverse, total length
  vector<double> dists_inv;
  for(const double d: dD_dTheta){
    const double d_inv = 2*d_avg - d; // a - (x - a)
    dists_inv.push_back(d_inv); // 'inverse' derivative
  }
  // cout << "dists inv: " << dists_inv << endl;

  // accumulated inverse distances
  vector<double> dists_inv_acc(n_points);
  double acc=0;
  const double fac = 2*M_PI/length_tot;
  // cout << "fac: " << fac << endl;
  for(int i=0; i< n_points; i++){
    dists_inv_acc[i] = acc;
    acc += dists_inv[i]*fac;
  }
  // cout << "dists inv acc: " << dists_inv_acc << endl;

  return dists_inv_acc;
}



#endif

