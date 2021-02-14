//
// Created by paul.ilardi on 2/12/2021.
//

#include "../units.h"
#include <iostream>
#include <iomanip>

struct ECEFCat{};
struct LLACat{};
struct ENUCat{};

template <typename C_, typename U1, typename U2, typename U3>
struct Pos3d{
  using FloatType = typename U1::FloatType;
  using C = C_;

  Pos3d()=default;

  template <typename Uo>
  Pos3d(const Uo& other) { set(other); }

  Pos3d(FloatType a,FloatType b,FloatType c){ set(a,b,c); }

  template <typename Uo>
  Pos3d operator=(const Uo& ecef2) { set(ecef2); return *this; }

  void set(FloatType a_, FloatType b_, FloatType c_){ a.set(a_); b.set(b_); c.set(c_); }
  void set(U1 a_, U2 b_, U3 c_){ a.set(a_); b.set(b_); c.set(c_); }

  template <typename Uo>
  void set(const Uo& other) {
    convert_from(other);
  }

  FloatType* data(){ return (FloatType*)&a; };

  static_assert(std::is_same<typename U1::FloatType,typename U2::FloatType>::value,"");
  static_assert(std::is_same<typename U1::FloatType,typename U3::FloatType>::value,"");

  U1 a;
  U2 b;
  U3 c;

private:
  static constexpr double squared(double v) noexcept {return v*v;};
  static constexpr double cubed(double v) noexcept {return v*v*v;};
  static constexpr double WGS84_A = 6378137;
  static constexpr double WGS84_F = 1/298.257223563;
  static constexpr double WGS84_B = WGS84_A*(1-WGS84_F);
  static constexpr double WGS84_E2 = 1 - squared(1-WGS84_F);

  template <typename Uo, typename std::enable_if<std::is_same<C, typename Uo::C>::value, bool>::type = true>
  void convert_from(const Uo& other){ set(other.a,other.b,other.c); }

  template <typename Uo, typename std::enable_if<std::is_same<typename Uo::C, LLACat>::value && std::is_same<C, ECEFCat>::value, bool>::type = true>
  void convert_from(const Uo& other){
    Pos3d<ECEFCat,Units::Meters,Units::Meters,Units::Meters> ecef;
    lla2ecef(other, ecef);
    *this = ecef;
  }

  template <typename Uo, typename std::enable_if<std::is_same<typename Uo::C, ECEFCat>::value && std::is_same<C, LLACat>::value, bool>::type = true>
  void convert_from(const Uo& other){
    Pos3d<LLACat,Units::Radians,Units::Radians,Units::Meters> lla;
    ecef2lla(other, lla);
    *this = lla;
  }

  void ecef2lla(const Pos3d<ECEFCat,Units::Meters,Units::Meters,Units::Meters>& ecef,
                Pos3d<LLACat,Units::Radians,Units::Radians,Units::Meters>& lla)
  {
    double x = ecef.a.get();
    double y = ecef.b.get();
    double z = ecef.c.get();

    double a = WGS84_A;
    double e2 = WGS84_E2;
    double b = WGS84_B;
    double ep = sqrt((a*a - b*b)/(b*b));
    double p = sqrt(x*x+y*y);
    double th = atan2(a*z, b*p);

    double lng = atan2(y,x);
    double lat = atan2(z + ep*ep*b*cubed(sin(th)), p - e2*a*cubed(cos(th)));
    double N = WGS84_A / sqrt(1.0 - WGS84_E2 * squared(sin(lat)));
    double alt = p / cos(lat) - N;

    lla.set(lat, lng, alt);
  }

  void lla2ecef(const Pos3d<LLACat,Units::Radians,Units::Radians,Units::Meters>& lla,
                Pos3d<ECEFCat,Units::Meters,Units::Meters,Units::Meters>& ecef)
  {
    double alt = lla.c.get();
    double clat = cos(lla.a);
    double slat = sin(lla.a);
    double clon = cos(lla.b);
    double slon = sin(lla.b);

    double N = WGS84_A / sqrt(1.0 - WGS84_E2 * slat * slat);

    ecef.a.set((N + alt) * clat * clon);
    ecef.b.set((N + alt) * clat * slon);
    ecef.c.set((N * (1.0 - WGS84_E2) + alt) * slat);
  }
};

template <typename U>
struct ECEF : public Pos3d<ECEFCat,U,U,U> {
  using Pos3d<ECEFCat,U,U,U>::Pos3d;
  U& x(){return this->a;}
  U& y(){return this->b;}
  U& z(){return this->c;}
  const U& x() const{return this->a;}
  const U& y() const{return this->b;}
  const U& z() const{return this->c;}
};

template <typename Ull, typename Ua>
struct LLA : public Pos3d<LLACat,Ull,Ull,Ua>{
  using Pos3d<LLACat,Ull,Ull,Ua>::Pos3d;
  Ull& lat(){return this->a;}
  Ull& lng(){return this->b;}
  Ua& alt(){return this->c;}
  const Ull& lat() const{return this->a;}
  const Ull& lng() const{return this->b;}
  const Ua& alt() const{return this->c;}
};

using ECEF_m = ECEF<Units::Meters>;
using LLA_ddm = LLA<Units::Degrees, Units::Meters>;
using LLA_rrm = LLA<Units::Radians, Units::Meters>;

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
}
