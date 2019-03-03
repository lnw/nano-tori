// Copyright (c) 2019, Lukas Wirz
// All rights reserved.

// This file is part of 'nano-tori' which is released under the BSD-2-clause license.
// See file LICENSE in this project.


#include <cmath>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>

#include "auxiliary.hh"
#include "rational.hh"
#include "geometry2.hh"
#include "geometry3.hh"


using namespace std;

#define TORUS 1
#define TUBE 0


int factorial(const int i){
  if(i==0) return 1;
  return i * factorial(i-1);
}

int double_factorial(const int i){
  if(i==0 || i==1) return 1;
  return i * double_factorial(i-2);
}

// after Ivory/Bessel
// minor/major half axis, approximation order
double ellipse_circumference(const double a, const double b, const int n){
  const double h = ((a-b)*(a-b))/((a+b)*(a+b));
  double C = 1 + 0.25*h;
  for(int i=2; i<=n; i++){
    C += pow( double_factorial(2*i-3) / (double_factorial(2*i) ), 2) * pow(h, i);
  }
  C *= M_PI*(a+b);
  return C;
}


// what we're doing is: we need the integral of the product of two functions f
// and g.  f is a function of the angle between the radius and the tangent
// vector.  g is a function of the length of the radius.  For reasons of
// simplicity, we get the Fourier-Cos series of f and g, and expand each
// coefficient in powers of alpha/r.  That works fairly well.

double approximate_that_bloody_integral(const double r0, const double alpha, const double phi){

// change sign
// f0 = r + alpha^2 + alpha^4
// f4 = -alpha^2 - alpha^4
// f8 = -0.25 alpha^4
//
// g0 = 1 + 0.25 alpha^2
// g2 = -alpha + 0.25 alpha^4
// g4 = -0.25 alpha^2
// g6 = -0.25 alpha^4

  const double dr = alpha/r0;

  const double f0 = 1 + pow(dr, 2) + pow(dr, 4);
  const double f4 = pow(dr, 2) + pow(dr, 4);
  const double f8 = 0.25 * pow(dr, 4);
// cout << "f048: " << f0 << ", " << f4 << ", " << f8 << endl;

  const double g0 = 1 + 0.25*pow(dr, 2);
  const double g2 = -dr + 0.25*pow(dr, 4);
  const double g4 = -0.25*pow(dr, 2);
  const double g6 = -0.25*pow(dr, 4);
// cout << "g0246: " << g0 << ", " << g2 << ", " << g4 << ", " << g6 << endl;

  double integral = 0;
  // f0 g0
  integral += f0 * g0 * phi;
  // f4 g0
  integral += f4 * g0 * (sin(4*phi)/4);
  // f8 g0
  integral += f8 * g0 * (sin(8*phi)/8);

  // f0 g2
  integral += f0 * g2 * (sin(2*phi)/2);
  // f4 g2
  integral += f4 * g2 * (sin(2*phi)/4 + sin(6*phi)/12);
  // f8 g2
  integral += f8 * g2 * (sin(6*phi)/12 + sin(10*phi)/20);

  // f0 g4
  integral += f0 * g4 * (sin(4*phi)/4);
  // f4 g4
  integral += f4 * g4 * (phi/2 + sin(8*phi)/16);
  // f8 g4
  integral += f8 * g4 * (sin(4*phi)/8 + sin(12*phi)/24);

  // f0 g6
  integral += f0 * g6 * (sin(6*phi)/6);
  // f4 g6
  integral += f4 * g6 * (sin(2*phi)/4 + sin(10*phi)/20);
  // f8 g6
  integral += f8 * g6 * (sin(2*phi)/4 + sin(14*phi)/28);

  return integral;
}



int main(int ac, char **av) {

// cerr << ac << endl;
  if (ac != 8) {
    cout << "usage: " << av[0] << " <m> <n> <p> <q> <phase> <bl> <ell>" << endl;
    cout << "  <m>: first component of the chiral vector" << endl;
    cout << "  <n>: second component of the chiral vector" << endl;
    cout << "  <p>: first component of the second side" << endl;
    cout << "  <q>: second component of the second side" << endl;
    cout << "  <phase>: where is the first atom [in rad], eg `0'?" << endl;
    cout << "  <target bond length>: 1.43, maybe?" << endl;
    cout << "  <elliptic parameter>: choose `1' for none" << endl;
    abort();
  }
  int m = stol(av[1], 0, 0);
  int n = stol(av[2], 0, 0);
  int p = stol(av[3], 0, 0);
  int q = stol(av[4], 0, 0);
  const double phase = stod(av[5], 0);
  const double CC_bondlength = stod(av[6], 0);
  const double alpha = stod(av[7], 0);

  int m_print = m, n_print = n;

  bool invert;
  if(n>m){
    tie(n,m) = make_tuple(m,n); // swap
    m_print = n, n_print = m;
    invert = true;
  }else invert = false;

  cout << "This program will generate tori which may be locally sheared. The number of atoms can be 4*n or 4*n+2." << endl;
  cout << "m, n, p, q: " << m << ", " << n << ", " << p << ", " << q << endl;
  cout << "phase: " << phase << endl;
  cout << "target bond length: " << CC_bondlength << endl;
  cout << "elliptic parameter: " << alpha << endl;

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
  const double rz = r*alpha;
  const double dr = rz - r;
  const double rxy = r - dr;
  const double R = (scnd_v_length * sin(angle_diff))/(2*M_PI); // cos (x-pi/2) = sin (x)
  cout << "R (major radius): " << R << endl;
  cout << "r (minor radius): " << r << endl;
  cout << "rxy (minor radius in horizontal direction): " << rxy << endl;
  cout << "rz (minor radius in vertical direction): " << rz << endl;
  cout << "radius of the hole (R - rxy): " << R-rxy << endl;

// here one would only print some hints, but every choice is possible
#if 0
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
  if(abs(lambda) < eps){ cout << "you could have used torus-simple for this structure." << endl;}
  matrix2d shear_mat(1, lambda, 0, 1); // shear parallel to x axis (which now coincides with the chiral vector)
  cout << "shear matrix: " << shear_mat << endl;
  S2 = shear_mat * S2;

  // crop
  S2.crop(0-eps, 0-eps, 2*M_PI*r-eps, 2*M_PI*R-eps);
// cout << "crop2 " << S2.size() << endl;

  const double integral_over_2pi = approximate_that_bloody_integral(r, dr, 2*M_PI);
  const double integral_over_2pi_normalised = integral_over_2pi / (2*M_PI);
  // cout << integral_over_2pi << endl;
  // cout << integral_over_2pi_normalised << endl;

  structure3d S3;
  for(size_t i=0; i<S2.size(); i++){
    const double phi  ( (S2[i][0])/r + phase );
    const double theta( (S2[i][1])/R );
// cout << "phi, theta: " << theta << ", " << phi << endl;
#if TORUS
    // modulate phi such that in a torus with ellipsoidal cross section the atoms are less dense at top/bottom and less sparse on the sides
    const double phi_prime = approximate_that_bloody_integral(r, dr, phi) / integral_over_2pi_normalised;

    // cout << phi << " --> " << phi_prime << " (" << 2*M_PI << ")" << endl;
    coord3d c3d( (R+rxy*cos(phi_prime))*cos(theta), (R+rxy*cos(phi_prime))*sin(theta), rz*sin(phi_prime) );
#elif TUBE
    coord3d c3d( r*cos(phi), r*sin(phi), R*theta );
#else
    coord3d c3d( r*phi, R*theta, 0 );
#endif
// cout << c3d << endl;
    S3.push_back(c3d);
    S3.push_back_atom(string("C"));
  }

  // S3.atomtypes[0] = "O";

  S3.scale(CC_bondlength);
  if (invert) S3.invert();
// cout << "S: " << S << endl;

  const int atom_count( 2*((m+n)*q-n*(p+q)) ); // Kirby 1993 faraday (but they define m, n, p, q differently)
  cout << atom_count << " atoms" << endl;

  string basename("torus-" + to_string(m_print) + "-" + to_string(n_print) + "-" + to_string(p) + "-" + to_string(q) + "-"
                           + to_string_with_precision(phase) + "-"
                           + to_string_with_precision(CC_bondlength) + "-" + to_string_with_precision(alpha) + "-" + to_string(atom_count));
  ofstream xyz((basename + ".xyz").c_str());
  ofstream turbo((basename + ".coord").c_str());
  xyz << S3.to_xyz();
  turbo << S3.to_turbomole();
  xyz.close();
  turbo.close();

  return 0;
}

