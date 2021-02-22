// Copyright (C) 2021 Paul Ilardi (http://github.com/CodePi)
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, unconditionally.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#include "../units_coords.h"

int main(){
  LLA_ddm lla_ddm = {42.346735, -71.097223, 300};
  ECEF_m ecef_m = lla_ddm;
  LLA_ddm lla_ddm2 = ecef_m;
  ECEF_m ecef_m2 = lla_ddm2;

  std::cout.precision(10);
  std::cout << lla_ddm.lat() << " " << lla_ddm.lng() << " " << lla_ddm.alt() << "\n";
  std::cout << ecef_m.x() << " " << ecef_m.y() << " " << ecef_m.z() << "\n";
  std::cout << lla_ddm2.lat() << " " << lla_ddm2.lng() << " " << lla_ddm2.alt() << "\n";
  std::cout << ecef_m2.x() << " " << ecef_m2.y() << " " << ecef_m2.z() << "\n";
  std::cout << "\n";

  // enu
  LLA_ddm lla_ddm_other = {42.446735, -71.197223, 100};
  ENU_m enu_m(lla_ddm_other, lla_ddm);
  LLA_ddm lla_ddm_other2 = enu_m.get<LLA_ddm>(lla_ddm);
  std::cout << lla_ddm_other.lat() << " " << lla_ddm_other.lng() << " " << lla_ddm_other.alt() << "\n";
  std::cout << enu_m.x() << " " << enu_m.y() << " " << enu_m.z() << "\n";
  std::cout << lla_ddm_other2.lat() << " " << lla_ddm_other2.lng() << " " << lla_ddm_other2.alt() << "\n";
  std::cout << "\n";

  // distance
  std::cout << lla_ddm.distance(lla_ddm_other) << "\n";
  std::cout << "\n";

  // compare
  assert(lla_ddm.distance(ecef_m) < Units::Meters(1));

  // wtf
  LLA<Units::Radians, Units::Furlongs> wtf = ecef_m;
  std::cout << wtf.lat() << " " << wtf.lng() << " " << wtf.alt() << "\n";
}
