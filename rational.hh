// Copyright (c) 2019, Lukas Wirz
// All rights reserved.

// This file is part of 'nano-tori' which is released under the BSD-2-clause license.
// See file LICENSE in this project.

#ifndef RATIONAL_HH
#define RATIONAL_HH

#include <cassert>
#include <iostream>

using namespace std;


int gcd(const int a, const int b){
  if (b == 0) return a;
  else return gcd(b, a % b);
}


class rational{
  int num, den;

  void simplify(){
    const int x(gcd(num, den));
    den /= x;
    num /= x;
  }

public:
  rational(const int x=0, const int y=1): num(x), den(y) {simplify();}
  rational(const rational &r): num(r.get_num()), den(r.get_den()) {}

  rational operator*(const int &i) const { return rational(*this) *= i; }
  rational operator/(const int &i) const { return rational(*this) /= i; }
  rational& operator*=(const int& i){ num *= i; simplify(); return *this; }
  rational& operator/=(const int& i){ den *= i; simplify(); return *this; }

//   double operator*(const double &d) const { return (num*d)/den; }
//   double operator/(const double &d) const { return num/(d*den); }

  rational operator+(const rational &r) const { return rational(*this) += r; }
  rational operator-(const rational &r) const { return rational(*this) -= r; }
  rational& operator+=(const rational& r){ num = num*r.den + r.num*den; den *= r.den; simplify(); return *this; }
  rational& operator-=(const rational& r){ num = num*r.den - r.num*den; den *= r.den; simplify(); return *this; }

  int get_num() const {return num;}
  int get_den() const {return den;}

  friend ostream& operator<<(ostream& S, const rational &r) {
    S << "{" << r.num << ", " << r.den << "}";
    return S;
  }

};

#endif

