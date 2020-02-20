// Copyright (c) 2019, Lukas Wirz
// All rights reserved.

// This file is part of 'nano-tori' which is released under the BSD-2-clause license.
// See file LICENSE in this project.

#include <array>
#include <cstdio>
#include <cerrno>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include "geometry3.hh"

using namespace std;

std::string exec(const char* cmd) {
  std::array<char, 128> buffer;
  std::string result;
  std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
  if (!pipe) throw std::runtime_error("popen() failed!");
  while (!feof(pipe.get())) {
    if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
      result += buffer.data();
  }
  return result;
}

vector<double> read_energy(const string& filename)
{
  ifstream file(filename.c_str());
  string whatever, line;
  vector<double> energies;
  getline(file,whatever);
  while( getline(file,line) ){
    string l(line);
    try {
      double e(stod(l, 0));
      energies.push_back(e);
    } catch (invalid_argument&){
    }
  }
  cout << "read energies: " << energies << endl;
  file.close();
  return energies;
}


struct params_t {
  int m, n, p, q;
  double phase;
  size_t n_eval;
  vector<double> energies;
  params_t(const int xm, const int xn, const int xp, const int xq, const double phasex, const size_t xne, const vector<double> xe): m(xm), n(xn), p(xp), q(xq), phase(phasex), n_eval(xne), energies(xe) {}
};


double run_torus(const int m, const int n, const int p , const int q, const double phase, const double CC, const double alpha) {

  // make torus
  string make_torus_command ( "./torus-full " + to_string(m) + " " + to_string(n) + " " +
                                              to_string(p) + " " + to_string(q) + " " + to_string_with_precision(phase) + " " +
                                              to_string_with_precision(CC) + " " + to_string_with_precision(alpha) );
  cout << make_torus_command << endl;
  exec(make_torus_command.c_str());
  string cp_command ( "cp torus*.coord coord; mv torus*.coord structures; mv torus*.xyz structures" );
  cout << cp_command << endl;
  exec(cp_command.c_str());
  
  // calculate energy using turbomole
  string tm_command ( "ridft > ridft.out" );
  cout << tm_command << endl;
  exec(tm_command.c_str());

  // get energy from 'energy'
  string get_energy_command ( "tail -n 2 energy | head -n 1 | awk {'print $2'}" );
  cout << get_energy_command << endl;
  const string e_out = exec(get_energy_command.c_str());
  const double E = stod(e_out.c_str(), 0);
  return E;
}


double energy(const gsl_vector *torus_specs, void *parameters) {
  //cout << "entering polyhedron pot" << endl;

  params_t &params = *static_cast<params_t *>(parameters);
  const int m = params.m;
  const int n = params.n;
  const int p = params.p;
  const int q = params.q;
  const double phase = params.phase;
  size_t &n_e = params.n_eval;
  vector<double> energies = params.energies;

  const double cc = gsl_vector_get( torus_specs, 0 );
  const double alpha  = gsl_vector_get( torus_specs, 1 );

  // if counter < number of stored energies, just return energy and increment counter
  n_e++;
  cout << "evaluation no: " << n_e << endl;
  if(n_e < energies.size()){
    cout << "we already know energy[" << n_e << "]: " << energies[n_e] << endl;
    return energies[n_e];
  }
  
  const double E = run_torus(m, n, p, q, phase, cc, alpha);
  cout << "new energy: " << setprecision(14) << E << endl;
  return E;
}


bool optimise(const int m, const int n, const int p, const int q, const double phase, const double cc, const double e,
              const double tol, const double step_r, const double step_e, const size_t max_iterations) {

  vector<double> energies = read_energy("energy");
  int n_eval = -1;

  params_t params(m, n, p, q, phase, n_eval, energies);

  gsl_multimin_function my_func;
  my_func.n = 2;
  my_func.f = energy;
  my_func.params = static_cast<void *>(&params);

  gsl_vector *torus_specs = gsl_vector_alloc(my_func.n);
  gsl_vector_set(torus_specs, 0, cc);
  gsl_vector_set(torus_specs, 1, e);

  gsl_vector *init_stepsize = gsl_vector_alloc(my_func.n);
  gsl_vector_set(init_stepsize, 0, step_r);
  gsl_vector_set(init_stepsize, 1, step_e);

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, my_func.n);
  gsl_multimin_fminimizer_set(s, &my_func, torus_specs, init_stepsize);

  size_t iter = 0;
  int status;
  cout << "setup complete" << endl;
  do {
    iter++;
    if (iter % 10 == 0) cout << iter << endl;

    status = gsl_multimin_fminimizer_iterate(s);
    if (status) break;

	// The minimizer-specific characteristic size is calculated as the average
	// distance from the geometrical center of the simplex to all its vertices. 
    const double size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, tol);
    cout << "size: " << size << endl;
    if (status == GSL_SUCCESS) {
      cerr << "Final energy: " << -s->fval << endl;
    }
    cout << "--- --- --- ---" << endl;
  } while (status == GSL_CONTINUE && iter < max_iterations);

  const double final_cc = gsl_vector_get(s->x, 0);
  const double final_alpha = gsl_vector_get(s->x, 1);
  cout << "Final parameters: " << final_cc << ", " << final_alpha << endl;
  cout << "run calculation with final parameters ..." << endl;
  const double final_E = run_torus(m, n, p, q, phase, final_cc, final_alpha);
  cout << "final energy: " << final_E << endl;

  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(torus_specs);

  cerr << "status: " << status << endl; // 0: success; -2: not converged; 27: 'not making progress towards solution'
  return status == 0 ? false : true;
}


int main(int ac, char **av) {

  if (ac != 11) {
    cout << "very non-robust two-parameter optimisation of nano tori.  Use at your own risk." << endl;
    cout << "usage: " << av[0] << " <m> <n> <p> <q> <ph> <bl> <ell> <tol> <st1> <st2>" << endl;
    cout << "  <m>: first component of the chiral vector" << endl;
    cout << "  <n>: second component of the chiral vector" << endl;
    cout << "  <p>: first component of the second side" << endl;
    cout << "  <q>: second component of the second side" << endl;
    cout << "  <phase>: in multiples of pi" << endl;
    cout << "  <target bond length>: 1.43, maybe?" << endl;
    cout << "  <ellipticity>: choose '1' for none" << endl;
    cout << "  <tolerance>: the convergence criterion.  Try 0.01" << endl;
    cout << "  <step 1>: initial step of the scaling factor.  Try 0.05" << endl;
    cout << "  <step 2>: initial step of the ellipticity.  Try 0.05" << endl;
    abort();
  }
  int m = stol(av[1], 0, 0);
  int n = stol(av[2], 0, 0);
  int p = stol(av[3], 0, 0);
  int q = stol(av[4], 0, 0);
  double phase = stod(av[5], 0);
  double CC_bondlength = stod(av[6], 0);
  double alpha = stod(av[7], 0);
  double tol = stod(av[8], 0);
  double step_r = stod(av[9], 0);
  double step_e = stod(av[10], 0);

  string mkdir_command ( "mkdir -p structures; mv torus*coord structures; mv torus*xyz structures" );
  cout << mkdir_command << endl;
  exec(mkdir_command.c_str());

  const size_t max_it = 50;
  if( optimise(m, n, p , q, phase, CC_bondlength, alpha, tol, step_r, step_e, max_it) )
    cerr << "something went wrong" << endl;
  else
    cerr << "looks good" << endl;

  return 0;
}

