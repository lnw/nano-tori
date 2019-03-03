// Copyright (c) 2019, Lukas Wirz
// All rights reserved.

// This file is part of 'nano-tori' which is released under the BSD-2-clause license.
// See file LICENSE in this project.


#include <cmath>
#include <string>
#include <vector>
#include <tuple>

#include "auxiliary.hh"
#include "rational.hh"
#include "geometry2.hh"
#include "geometry3.hh"

using namespace std;

#define TORUS 1
#define TUBE 0


int main(int ac, char **av) {

// cerr << ac << endl;
  if (ac != 4) {
    cout << "usage: " << av[0] << " <m> <n> <l>" << endl;
    cout << "  <m>: first component of the chiral vector" << endl;
    cout << "  <n>: second component of the chiral vector" << endl;
    cout << "  <l>: length (not all lengths are possible!)" << endl;
    abort();
  }
  int m = stol(av[1], 0, 0);
  int n = stol(av[2], 0, 0);
  const int l = stol(av[3], 0, 0);

  int m_print = m, n_print = n;
  bool invert;
  if(n>m){
    tie(n,m) = make_tuple(m,n); // swap
    m_print = n, n_print = m;
    invert = true;
  }else invert = false;

  cout << "This program will generate tori with locally correct angles and 4*n atoms." << endl;
  cout << "m, n, l: " << m << ", " << n << ", " << l << endl;

  const double eps = 1.e-8;

  // the basis vectors in the hexagonal latice
  const coord2d a(sqrt(3), 0);
  const coord2d b(sqrt(3)/2, 1.5);

  // the chiral vector
  const coord2d ch_v(m*a + n*b);
  const double ch_v_length(ch_v.norm());
  const double ch_v_phi(atan2(ch_v[1], ch_v[0]));
// cout << "angle: " << angle << endl;

  const matrix2d rot_mat(cos(ch_v_phi), sin(ch_v_phi), -sin(ch_v_phi), cos(ch_v_phi));
  cout << "rotation matrix: " << rot_mat << endl;

  const double r = ch_v_length/(2*M_PI) ;
  const double R = (1.5*l / cos(ch_v_phi))/(2*M_PI); // 1/cos = sec
  cout << "R, r: " << R << ", " << r << " (major radius, minor radius)" << endl;

#if TORUS
  rational rat(m + 2*n, 2*m + n);  // after writing a page ...
  const int den(rat.get_den());

  if ( l % den != 0 ){
    cout << "Unfortunately the computer is out of " << l << "s" << endl;
    cout << "For m=" << m << " and n=" << n << ", l must be a multiple of " << rat.get_den() << "." << endl;
    cout << "May we interest you in l=" << floor(double(l)/den)*den << " or " << ceil(double(l)/den)*den << "?" << endl;
    return 0;
  }
  else cout << "all good" << endl;
#endif

  structure2d S2;

  for (int x=-2*l; x<m+n; x++){
    for (int y=0; y<2*l; y++){
// cout << "x,y:" << x << ", " << y << endl;

      coord2d cc1 = x*a + y*b + coord2d(0,1.0); // upper
      coord2d cc2 = x*a + y*b; // lower
// cout << cc1 << endl;
// cout << cc2 << endl;

      S2.push_back(cc1);
      S2.push_back(cc2);
    }
  }

  // rotate grid
  S2 = rot_mat * S2;

  // crop
  S2.crop(0-eps, 0-eps, 2*M_PI*r-eps, 2*M_PI*R-eps);
// cout << "crop2 " << S2.size() << endl;

  structure3d S3;
  for(size_t i=0; i<S2.size(); i++){
    const double phi  ( (S2[i][0])/r );
    const double theta( (S2[i][1])/R );
// cout << "phi, theta: " << theta << ", " << phi << endl;
#if TORUS
    coord3d c3d( (R+r*cos(phi))*cos(theta), (R+r*cos(phi))*sin(theta), r*sin(phi) );
#elif TUBE
    coord3d c3d( r*cos(phi), r*sin(phi), R*theta );
#else
    coord3d c3d( r*phi, R*theta, 0 );
#endif
// cout << c3d << endl;
    S3.push_back(c3d);
    S3.push_back_atom(string("C"));
  }

  const double CC_bondlength = 1.43;
  S3.scale(CC_bondlength);
  if (invert) S3.invert();
// cout << "S: " << S << endl;
  const int atom_count(S3.size());
  cout << atom_count << " atoms" << endl;

  string basename("torus-" + to_string(m_print) + "-" + to_string(n_print) + "-" + to_string(l) + "-" + to_string(atom_count));
  ofstream xyz((basename + ".xyz").c_str());
  ofstream turbo((basename + ".coord").c_str());
  xyz << S3.to_xyz();
  turbo << S3.to_turbomole();
  xyz.close();
  turbo.close();

  return 0;
}

