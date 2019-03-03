// Copyright (c) 2019, Lukas Wirz
// All rights reserved.

// This file is part of 'nano-tori' which is released under the BSD-2-clause license.
// See file LICENSE in this project.


#include <cmath>
#include <string>
#include <tuple>
#include <vector>

#include "auxiliary.hh"
#include "geometry2.hh"
#include "geometry3.hh"

using namespace std;

//
// input: m, n, p, q
//        index of hardcoded curve, which is given parametrically
//
// -- sample points on curve
// -- tabulate density function of points on the curve by inverting list of sampled points
// -- use space curve, modulated with the density function


coord3d circle(const double t){
  const double x = cos(t);
  const double y = sin(t);
  const double z = 0;
  return coord3d(x,y,z);
}
coord3d circle_tangent(const double t){
  const double x = -sin(t);
  const double y = cos(t);
  const double z = 0;
  return coord3d(x,y,z).normalised();
}

coord3d trefoil(const double t){
  const double x = sin(t) + 2*sin(2*t);
  const double y = cos(t) - 2*cos(2*t);
  const double z = -sin(3*t);
  return coord3d(x,y,z);
}
coord3d trefoil_tangent(const double t){
  const double x = cos(t) + 4*cos(2*t);
  const double y = -sin(t) + 4*sin(2*t);
  const double z = -3*cos(3*t);
  return coord3d(x,y,z).normalised();
}

coord3d lemniscate(const double t){
  const double a = 2.0;
  const double b = 0.3;
  const double x = a*cos(t)/(1+sin(t)*sin(t));
  const double y = a*sin(t)*cos(t)/(1+sin(t)*sin(t));
  const double z = b * sin(t);
  return coord3d(x,y,z);
}
coord3d lemniscate_tangent(const double t){
  const double a = 2.0;
  const double b = 0.3;
  const double x = -(2*a*(5+cos(2*t))*sin(t))/pow(-3+cos(2*t),2);
  const double y = (2*a*(-1+3*cos(2*t)))/pow(-3+cos(2*t),2);
  const double z = b * cos(t);
  return coord3d(x,y,z);
}

coord3d curve(int c, const double t)
{ 
  switch(c){
  case 1: return circle(t);
  case 2: return trefoil(t);
  case 3: return lemniscate(t);
  default:
    break;
  }
  cerr << "invalid curve chosen, aborting ..." << endl;
  abort();
}

coord3d curve_tan(int c, const double t)
{ 
  switch(c){
  case 1: return circle_tangent(t);
  case 2: return trefoil_tangent(t);
  case 3: return lemniscate_tangent(t);
  default:
    break;
  }
  cerr << "invalid curve chosen, aborting ..." << endl;
  abort();
}



int main(int ac, char **av) {

// cerr << ac << endl;
  if (ac != 7) {
    cout << "like torus-full, without ellipticity, but producing silly knots" << endl;
    cout << "usage: " << av[0] << " <m> <n> <p> <q> <bl> <structure>" << endl;
    cout << "  <m>: first component of the chiral vector" << endl;
    cout << "  <n>: second component of the chiral vector" << endl;
    cout << "  <p>: first component of the second side" << endl;
    cout << "  <q>: second component of the second side" << endl;
    cout << "  <target bond length>: 1.43, maybe?" << endl;
    cout << "  <structure>: select by index" << endl;
    abort();
  }
  int m = stol(av[1], 0, 0);
  int n = stol(av[2], 0, 0);
  int p = stol(av[3], 0, 0);
  int q = stol(av[4], 0, 0);
  const double CC_bondlength = stod(av[5], 0);
  const int n_curve = stol(av[6], 0, 0);

  int m_print = m, n_print = n;

  bool invert;
  if(n>m){
    tie(n,m) = make_tuple(m,n); // swap
    m_print = n, n_print = m;
    invert = true;
  }else invert = false;

  cout << "This program will generate tori which may be locally sheared. The number of atoms can be 4*n or 4*n+2." << endl;
  cout << "m, n, p, q: " << m << ", " << n << ", " << p << ", " << q << endl;
  cout << "target bond length: " << CC_bondlength << endl;
  cout << "chosen curve: " << n_curve << endl;

  const double eps = 1.e-8;

  // the basis vectors in the hexagonal latice
  const coord2d a(sqrt(3), 0);
  const coord2d b(sqrt(3)/2, 1.5);

  // the chiral vector
  const coord2d ch_v(m*a + n*b);
  const double ch_v_length(ch_v.norm());
  const double ch_v_phi(atan2(ch_v[1], ch_v[0]));

  const matrix2d rot_mat(cos(ch_v_phi), sin(ch_v_phi), -sin(ch_v_phi), cos(ch_v_phi));
// cout << "angle: " << angle << endl;
  cout << "rotation matrix: " << rot_mat << endl;

  // the other vector
  const coord2d scnd_v(p*a + q*b);
  const double scnd_v_length(scnd_v.norm());
  const double scnd_v_phi(atan2(scnd_v[1], scnd_v[0]));
// cout << "vec, length, angle: " << scnd_v << ", " << scnd_v_length << ", " << scnd_v_phi << endl;

  const double angle_diff = scnd_v_phi - ch_v_phi;
// cout << "angle1, angle2, diff: " << ch_v_phi << ", " <<  scnd_v_phi << ", " << angle_diff << endl;
  if(angle_diff > M_PI){
    cerr << "the chosen p, q imply an angle larger than pi" << endl;
  }
  if(angle_diff < 0){
    cerr << "the chosen p, q imply an angle smaller than 0" << endl;
  }

  // get r and R
  const double r = ch_v_length/(2*M_PI);
  const double R = (scnd_v_length * sin(angle_diff))/(2*M_PI); // cos (x-pi/2) = sin (x)
  cout << "R (major radius): " << R << endl;
  cout << "r (minor radius): " << r << endl;

  structure2d S2;

  // populate grid
  for (int x=2*min(0,p); x<2*max(m,p); x++){
    for (int y=2*min(0,q); y<2*max(n,q); y++){
// cout << "x,y:" << x << ", " << y << endl;

      coord2d cc1 = coord2d(x*a + y*b) + coord2d(0,1.0); // upper
      coord2d cc2 = coord2d(x*a + y*b); // lower
// cout << cc1 << endl;
// cout << cc2 << endl;

      S2.push_back(cc1);
      S2.push_back(cc2);
    }
  }

  // rotate grid
  S2 = rot_mat * S2;

  // shear grid
  const double lambda = tan(angle_diff - M_PI/2);
  // if(abs(lambda) < eps){ cout << "you could have used torus-simple for this structure." << endl;}
  matrix2d shear_mat(1, lambda, 0, 1); // shear parallel to x axis (which now coincides with the chiral vector)
  cout << "shear matrix: " << shear_mat << endl;
  S2 = shear_mat * S2;

// cout << coord2d(0,0) << ", " << ch_v << ", " << scnd_v << ", " << ch_v + scnd_v << endl;
// cout << rot_mat * coord2d(0,0) << ", " << rot_mat * ch_v << ", " << rot_mat * scnd_v << ", " << rot_mat * (ch_v + scnd_v) << endl;
// cout << shear_mat * (rot_mat * coord2d(0,0)) << ", "
//      << shear_mat * (rot_mat * ch_v) << ", "
//      << shear_mat * (rot_mat * scnd_v) << ", "
//      << shear_mat * (rot_mat * (ch_v + scnd_v)) << endl;
//
// cout << S2.contains(coord2d(0,0))                << endl;
// cout << S2.contains(coord2d(0,2*M_PI*R))         << endl;
// cout << S2.contains(coord2d(2*M_PI*r,0))         << endl;
// cout << S2.contains(coord2d(2*M_PI*r, 2*M_PI*R)) << endl;


  // crop
  S2.crop(0-eps, 0-eps, 2*M_PI*r-eps, 2*M_PI*R-eps);
// cout << "crop2 " << S2.size() << endl;

//  const double integral_over_2pi = 1; //approximate_that_bloody_integral(r, dr, 2*M_PI);
//  const double integral_over_2pi_normalised = integral_over_2pi / (2*M_PI);
  // cout << integral_over_2pi << endl;
  // cout << integral_over_2pi_normalised << endl;

  structure3d S3;

  const int n_points(10000);
  for (int i = 0; i<n_points; i++){
    const double ii = 2*M_PI*i / n_points;
    const coord3d c3d = curve(n_curve, ii);
    S3.push_back(c3d);
    S3.push_back_atom(string("C"));
  }




  // distances
  vector<double> dists;
  for(size_t i = 0; i< S3.size(); i++){
    coord3d previous = S3[(i+S3.size()-1)%S3.size()];
    coord3d here = S3[i];
    coord3d next = S3[(i+1)%S3.size()];
    double dist = ((here-next).norm() + (here-previous).norm())/2;
    //cout << (here-next).norm() << ", " <<  (here-previous).norm() << ", " << ((here-next).norm() + (here-previous).norm())/2 << endl;
    dists.push_back(dist);
  }
  //cout << "dists : " << dists << endl;
  // distances inverse, total length
  double length_tot = 0;
  for(const double dist: dists){
    length_tot += dist;
  }
  cout << "tot: " << length_tot << endl;
  double dist_avg = length_tot / n_points;
  cout << "avg: " << dist_avg << endl;
  // distances inverse, total length
  vector<double> dists_inv;
  for(const double dist: dists){
    double dist_inv = 2*dist_avg - dist; // a - (x - a)
    dists_inv.push_back(dist_inv);
  }
  //cout << "dists inv: " << dists_inv << endl;
  // accumulated inverse distances
  vector<double> dists_inv_acc;
  double tmp=0;
  const double fac = 2*M_PI/length_tot;
  for(const double dist_inv: dists_inv){
    tmp += dist_inv*fac;
    dists_inv_acc.push_back(tmp);
  }
  //cout << "dists inv acc: " << dists_inv_acc << endl;


  const double factor = 2* M_PI * R / length_tot;


  structure3d S4;

  for(size_t i=0; i<S2.size(); i++){
    const double phi  ( (S2[i][0])/r );
    const double theta( (S2[i][1])/R );
    // cout << "phi, theta: " << theta << ", " << phi << endl;
    // modulate phi such that in a torus with ellipsoidal cross section the atoms are less dense at top/bottom and less sparse on the sides
    double theta_prime;
    if (theta*n_points/(2*M_PI) - floor( theta*n_points/(2*M_PI) ) < 1.e-3){
      int index = int(theta*n_points/(2*M_PI) + 0.5);
      theta_prime = dists_inv_acc[index];
    }
    else {
      int index1 = int(theta*n_points/(2*M_PI) + 0.5);
      double weight1 = theta*n_points/(2*M_PI) - index1;
      int index2 = int(theta*n_points/(2*M_PI) + 1.5);
      double weight2 = index2 - theta*n_points/(2*M_PI);
      // cout << index1 << ", " << index2 << endl;
      // cout << weight1 << ", " << weight2 << endl;
      theta_prime = weight1 * dists_inv_acc[index1] + weight2 * dists_inv_acc[index2];
    }

    coord3d nn(curve(n_curve, theta_prime) * factor);
    coord3d tangent(curve_tan(n_curve, theta_prime));
    coord3d ref(0,0,1);
    coord3d va( (tangent.cross(ref)).normalised() );
    coord3d vb( (tangent.cross(va)).normalised() );
    coord3d c3d( nn + va*r*cos(phi) + vb*r*sin(phi) );
    // cout << c3d << endl;
    S4.push_back(c3d);
    S4.push_back_atom(string("C"));
  }


  // S4.atomtypes[0] = "O";
  S4.scale(CC_bondlength);
  if (invert) S4.invert();
// cout << "S: " << S << endl;


  const int atom_count( 2*((m+n)*q-n*(p+q)) ); // Kirby 1993 faraday
  cout << atom_count << " atoms" << endl;

  string basename("knot-" + to_string(n_curve) + "-" + to_string(m_print) + "-" + to_string(n_print) + "-" + to_string(p) + "-" + to_string(q) + "-"
                           + to_string_with_precision(CC_bondlength) + "-" + to_string(atom_count));
  ofstream xyz((basename + ".xyz").c_str());
  ofstream turbo((basename + ".coord").c_str());
  xyz << S4.to_xyz();
  turbo << S4.to_turbomole();
  xyz.close();
  turbo.close();

  return 0;
}

