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
#include "shape-gen.hh"

using namespace std;

//
// input: m, n, p, q
//        index of hardcoded curve, which is given parametrically
//
// -- sample points on curve
// -- tabulate density function of points on the curve by inverting list of sampled points
// -- use space curve, modulated with the density function



int main(int ac, char **av) {

// cerr << ac << endl;
  if (ac != 7) {
    cout << "like torus-full, without ellipticity, but producing silly knots" << endl;
    cout << "usage: " << av[0] << " <m> <n> <p> <q>" << endl;
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

  ShapeGenerator shape_gen;
  auto shape = shape_gen.get_shape(n_curve);
  auto shape_tan = shape_gen.get_shape_tan(n_curve);

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


  const int n_steps(1000);
  double length_tot;
  const vector<double> dists_inv_acc = shape_gen.get_inverse_mapping(n_curve, n_steps, length_tot);
  // cout << dists_inv_acc << endl;

  const double factor = 2* M_PI * R / length_tot;
  // cout << "fac: " << factor << endl;


  structure3d S3;

  for(size_t i=0; i<S2.size(); i++){
    const double phi  ( (S2[i][0])/r );
    const double theta( (S2[i][1])/R );
    // cout << "phi, theta: " << theta << ", " << phi << endl;
    // modulate phi such that in a torus with ellipsoidal cross section the atoms are less dense at top/bottom and less sparse on the sides
    double theta_prime;
    if ( abs(theta*n_steps/(2*M_PI) - round( theta*n_steps/(2*M_PI) )) < 1.e-3){
      const int index = int(theta*n_steps/(2*M_PI) + 0.5);
      theta_prime = dists_inv_acc[index];
    }
    else {
      const int index1 = int(theta*n_steps/(2*M_PI));
      const double weight2 = theta*n_steps/(2*M_PI) - index1;
      const int index2 = index1 + 1;
      const double weight1 = index2 - theta*n_steps/(2*M_PI);
      // cout << index1 << ", " << index2 << endl;
      // cout << weight1 << ", " << weight2 << endl;
      theta_prime = weight1 * dists_inv_acc[index1] + weight2 * dists_inv_acc[index2];
    }

    coord3d nn(shape(theta_prime) * factor);
    coord3d tangent(shape_tan(theta_prime));
    coord3d ref(0,0,1);
    coord3d va( (tangent.cross(ref)).normalised() );
    coord3d vb( (tangent.cross(va)).normalised() );
    coord3d c3d( nn + va*r*cos(phi) + vb*r*sin(phi) );
    // cout << c3d << endl;
    S3.push_back(c3d);
    S3.push_back_atom(string("C"));
  }


  S3.scale(CC_bondlength);
  if (invert) S3.invert();
// cout << "S: " << S << endl;


  const int atom_count( 2*((m+n)*q-n*(p+q)) ); // Kirby 1993 faraday
  cout << atom_count << " atoms" << endl;

  string basename("knot-" + to_string(n_curve) + "-" + to_string(m_print) + "-" + to_string(n_print) + "-" + to_string(p) + "-" + to_string(q) + "-"
                           + to_string_with_precision(CC_bondlength) + "-" + to_string(atom_count));
  ofstream xyz((basename + ".xyz").c_str());
  ofstream turbo((basename + ".coord").c_str());
  xyz << S3.to_xyz();
  turbo << S3.to_turbomole();
  xyz.close();
  turbo.close();

  return 0;
}

