//
// Created by paul.ilardi on 2/12/2021.
//

#include "../units.h"

template <typename U=Units::Meters>
struct ECEF{
  U x;
  U y;
  U z;

  using FloatType = typename U::FloatType;

  ECEF(FloatType x_, FloatType y_, FloatType z_) { set(x_,y_,z_); }

  template <typename U2>
  ECEF(const ECEF<U2>& ecef2) { set(ecef2); }

  template <typename U2>
  ECEF operator=(const ECEF<U2>& ecef2) { set(ecef2); return *this; }

  void set(FloatType x_, FloatType y_, FloatType z_){ x.set(x_); y.set(y_); z.set(z_); }

  template <typename U2>
  void set(const ECEF<U2>& ecef2) { x=ecef2.x; y=ecef2.y; z=ecef2.z; };

  FloatType* data(){ return (FloatType*)this; };

  operator U*(){ return (U*)this; }
};

template <typename Ull, typename Ua>
struct LLA{
  Ull lat;
  Ull lng;
  Ua  alt;

  static_assert(std::is_same<typename Ull::FloatType, typename Ua::FloatType>::value,"");

  using FloatType = typename Ull::FloatType;

  LLA(FloatType lat_, FloatType lng_, FloatType alt_) { set(lat_, lng_, alt_); }

  template <typename Ull2, typename Ua2>
  LLA(const LLA<Ull2, Ua2>& lla2) { set(lla2); }

  template <typename Ull2, typename Ua2>
  LLA operator=(const LLA<Ull2, Ua2>& lla2) { set(lla2); return *this; }

  void set(FloatType lat_, FloatType lng_, FloatType alt_){ lat.set(lat_); lng.set(lng_); alt.set(alt_); }

  template <typename Ull2, typename Ua2>
  void set(const LLA<Ull2, Ua2>& lla2) { lat=lla2.lat; lng=lla2.lng; alt=lla2.alt; };

  FloatType* data(){ return (FloatType*)this; };

};
// 42.346735, -71.097223

using ECEF_m = ECEF<Units::Meters>;
using LLA_ddm = LLA<Units::Degrees, Units::Meters>;
using LLA_rrm = LLA<Units::Radians, Units::Meters>;

int main(){
  using namespace Units::literals;
  ECEF_m ecef_m = {1529476, -4466529, 4274147};
  LLA_ddm lla_ddm = {42.346735, -71.097223, 0};

  ECEF<Units::Feet> ecef_ft = ecef_m;
  LLA_rrm llh_rrm = lla_ddm;
}