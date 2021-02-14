//
// Created by paul.ilardi on 2/12/2021.
//

#include "../units.h"
#include <iostream>
#include <iomanip>

template <typename U=Units::Meters>
struct ECEF{
  U x;
  U y;
  U z;

  using FloatType = typename U::FloatType;

  ECEF() = default;
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

using ECEF_m = ECEF<Units::Meters>;
using LLA_ddm = LLA<Units::Degrees, Units::Meters>;
using LLA_rrm = LLA<Units::Radians, Units::Meters>;

inline double squared(double v) noexcept {return v*v;};
inline double cubed(double v) noexcept {return v*v*v;};
static const double WGS84_A = 6378137;
static const double WGS84_F = 1/298.257223563;
static const double WGS84_B = WGS84_A*(1-WGS84_F);
static const double WGS84_E2 = 1 - squared(1-WGS84_F);

ECEF_m lla2ecef(const LLA_rrm& lla){
  double clat = cos(lla.lat);
  double slat = sin(lla.lat);
  double clon = cos(lla.lng);
  double slon = sin(lla.lng);

  double N = WGS84_A / sqrt(1.0 - WGS84_E2 * slat * slat);

  ECEF_m ecef;
  ecef.x.set((N + lla.alt.get()) * clat * clon);
  ecef.y.set((N + lla.alt.get()) * clat * slon);
  ecef.z.set((N * (1.0 - WGS84_E2) + lla.alt.get()) * slat);
  return ecef;
}

LLA_rrm ecef2lla(const ECEF_m& ecef){
  double x = ecef.x.get();
  double y = ecef.y.get();
  double z = ecef.z.get();

  double a = WGS84_A;
  double e2 = WGS84_E2;
  double b = WGS84_B;
  double ep = sqrt((a*a - b*b)/(b*b));
  double p = sqrt(x*x+y*y);
  double th = atan2(a*z, b*p);

  double lon = atan2(y,x);
  double lat = atan2(z + ep*ep*b*cubed(sin(th)), p - e2*a*cubed(cos(th)));
  double N = WGS84_A / sqrt(1.0 - WGS84_E2 * squared(sin(lat)));
  double alt = p / cos(lat) - N;

  return {lat, lon, alt};
}

int main(){
  using namespace Units::literals;
  LLA_ddm lla_ddm = {42.346735, -71.097223, 0};
  ECEF_m ecef_m = lla2ecef(lla_ddm);
  LLA_ddm lla_ddm2 = ecef2lla(ecef_m);
  ECEF_m ecef_m2 = lla2ecef(lla_ddm2);

  std::cout.precision(8);
  std::cout << lla_ddm.lat << " " << lla_ddm.lng << " " << lla_ddm.alt << "\n";
  std::cout << ecef_m.x << " " << ecef_m.y << " " << ecef_m.z << "\n";
  std::cout << lla_ddm2.lat << " " << lla_ddm2.lng << " " << lla_ddm2.alt << "\n";
  std::cout << ecef_m2.x << " " << ecef_m2.y << " " << ecef_m2.z << "\n";
}