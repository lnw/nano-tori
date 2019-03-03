
# Generation of nano-torus structures

## LICENSING / CITING

This software is published under a very permissive license to make it as useful
as possible for a variety of usecases.  If it is useful to you, especially in a
scientific context, you are encouraged to cite XXX. 

Some parts of geometry{2,3}.{hh,cc} and auxilary.hh have been borrowed from
http://ctcp.massey.ac.nz/index.php?group=&page=fullerenes .


## COMPILING

There are two ways to build these little programs.  Should work on any Linux
and related systems such as macOS.

### using the ordinary Makefile

make


### using CMake (probably more robust)

mkdir mybuild
cd mybuild
[CXX=clang++] cmake ..
make


## Usage

Each of the programs prints the list of requred parameters when it's called
with no or the wrong number of parameters.

### torus-simple

Generate a subset of the possible tori which has undistorted hexagons.

usage: ./torus-simple <m> <n> <l>
  <m>: first component of the chiral vector
  <n>: second component of the chiral vector
  <l>: length (not all lengths are possible!)

The resulting tori will have hexagons with ideal angles, and the atom count
will always be 4\*n.  If you chose an invalid `l', you will be suggested valid
choices for `l'

### torus-full

Generate any circular torus, and choose the ellipticity of the cross section.

usage: ./torus-full <m> <n> <p> <q> <phase> <bl> <ell>
  <m>: first component of the chiral vector
  <n>: second component of the chiral vector
  <p>: first component of the second side
  <q>: second component of the second side
  <phase>: where is the first atom [in rad], eg `0'?
  <target bond length>: 1.43, maybe?
  <elliptic parameter>: choose `1' for none

If the printed shear matrix if far from {{1,0},{0,1}} the obtained structure is
probably not sensible and other m,n,p,q should be chosen Tuning the phase is
mostly useful to get the maximal point group symmetry.  Note that the bond
length is only a target and there will always be a spectrum of shorter and
longer bonds (especially for large r and small R).  The elliptic parameter can
be left at 0.0 for very large R, but should be increased to 1.1 ... 1.5 for
smaller R or large r.

### torus-extra

Generate closed nano-tube rings which need not be tori.  Hardcoded shapes are a
torus (0), a lemniscate (2), and a trefoil (1).  Other shapes can easily be
defined.

like torus-full, without ellipticity, but producing silly knots
usage: ./torus-extra <m> <n> <p> <q> <bl> <structure>
  <m>: first component of the chiral vector
  <n>: second component of the chiral vector
  <p>: first component of the second side
  <q>: second component of the second side
  <target bond length>: 1.43, maybe?
  <structure>: select by index


### torus-grid

Generate grids (as xyz file) below, on, and above the surface of a torus.

This program generates grids to be used with the tori generated by torus-full.
usage: ./torus-grid <m> <n> <p> <q> <bl> <ell>
  <m>: first component of the chiral vector
  <n>: second component of the chiral vector
  <p>: first component of the second side
  <q>: second component of the second side
  <target bond length>: 1.43, maybe?
  <elliptic parameter>: choose '1' for none

(And before someone asks: the bond length is required for scaling, even though
there are no bonds.)

### botch_optimise-torus-full

If one is not sure about the best target bond length or the best elliptical
parameter but needs structures that are optimised at soe expensive qc-level of
theory, it can be efficient to preoptimise the structure in this two-parameter space.

The current version structly aims at turbomole, but it should be possible to
adapt it.  Basically, we calculate single points until simplex downhill is
converged.

very non-robust two-parameter optimisation of nano tori.  Use at your own risk.
usage: ./botch_optimise-torus-full <m> <n> <p> <q> <ph> <bl> <ell> <tol> <st1> <st2>
  <m>: first component of the chiral vector
  <n>: second component of the chiral vector
  <p>: first component of the second side
  <q>: second component of the second side
  <phase>: in multiples of pi
  <target bond length>: 1.43, maybe?
  <ellipticity>: choose '1' for none
  <tolerance>: the convergence criterion.  Try 0.01
  <step 1>: initial step of the scaling factor.  Try 0.05
  <step 2>: initial step of the ellipticity.  Try 0.05

