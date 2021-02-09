//
// Created by Paul Ilardi on 1/30/2021.
//

#include <iostream>
#include "eigen_units.h"
#include "../common/typestring.h"

using namespace std;
using namespace Units;

typedef Eigen::Matrix<Meters, -1, 1> MVec;
typedef Eigen::Matrix<Seconds, -1, 1> SVec;
typedef Eigen::Matrix<MetersPerSecond, -1, 1> MPSVec;

typedef Eigen::Matrix<Miles, -1, 1> MiVec;
typedef Eigen::Matrix<Hours, -1, 1> HVec;
typedef Eigen::Matrix<MPH , -1, 1> MPHVec;

int main(){
  size_t sz = 3;
  MVec m(sz);
  m << Meters(1000), Meters(2000), Meters(3000);
  SVec s(sz);
  s << Seconds(400), Seconds(500), Seconds(600);
  MiVec mi = m.cast<Miles>();
  HVec h = s.cast<Hours>();

  // basic math
  MVec m2 = m+m;
  MVec m3 = m*2;
  cout << "m2: " << m2.transpose() << "\n";
  cout << "m3: " << m3.transpose() << "\n";
  MVec m4 = m + mi.cast<Meters>();
  cout << "m4: " << m4.transpose() << "\n";

  // speed = distance / time
  MPSVec mps = (m.array()/s.array()).cast<MetersPerSecond>();  // converts units properly Meters/Seconds -> MetersPerSecond
  cout << "Meters: " << m.transpose() << "\n";
  cout << "Seconds: " << s.transpose() << "\n";
  cout << "MetersPerSecond: " << mps.transpose() << "\n";

  // calculate speed and convert
  MPHVec mph = (m.array()/s.array()).cast<MPH>();
  cout << "Miles: " << mi.transpose() << "\n";
  cout << "Hours: " << h.transpose() << "\n";
  cout << "MPH: " << mph.transpose() << "\n";

  // set values using IgnoreUnits
  // Last resort function to force eigen to ignore units and treat as double (use carefully)
  IgnoreUnits(m3).setRandom();
  cout << "m3 rand: " << m3.transpose() << "\n";

  // since eigen can't handle Area = Distance * Distance, must IgnoreUnits
  double norm = IgnoreUnits(m).norm();
  cout << "norm: " << norm << "\n";

  // test IgnoreUnits
  const MVec& mr = m;
  cout << "IgnoreUnits readwrite version: " << typestring(IgnoreUnits(m))   << "\n";  // readwrite version
  cout << "IgnoreUnits readonly version:  " << typestring(IgnoreUnits(mr))  << "\n";  // readonly version
  cout << "IgnoreUnits rvalue version:    " << typestring(IgnoreUnits(m+m)) << "\n";  // rvalue version

  // double = Unit/Unit
  // Should work like this, but Eigen doesn't like it.  See OpTraits note above about "ambiguous template"
  //Eigen::Vector3d v = m.array()/m2.array();
  // use IgnoreUnits instead
  Eigen::Vector3d v = IgnoreUnits(m).array()/IgnoreUnits(m2).array();
  cout << "v: " << v.transpose() << "\n";

  // sum works on units
  Meters d = (m+m).sum();
  cout << "(m+m).sum(): " << d << "\n";

  // apply to all
  MVec m5(100);
  m5.setOnes();
  m5 = m5.unaryExpr([](Meters s){return s*2;});
  assert(m5[0].get()==2);
}

