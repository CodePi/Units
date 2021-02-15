//
// Created by paul.ilardi on 2/12/2021.
//

#include "../units.h"
#include <iostream>
#include <iomanip>

enum Category{ECEFCat, LLACat, ENUCat};

template <Category C_, typename U1, typename U2, typename U3>
struct Pos3d{
  using FloatType = typename U1::FloatType;
  static constexpr Category C = C_;

  Pos3d()=default;

  template <typename Po>
  Pos3d(const Po& other) { set(other); }

  Pos3d(FloatType a,FloatType b,FloatType c){ set(a,b,c); }

  template <typename Po>
  Pos3d operator=(const Po& other) { set(other); return *this; }

  void set(FloatType a_, FloatType b_, FloatType c_){ v1.set(a_); v2.set(b_); v3.set(c_); }
  void set(U1 a_, U2 b_, U3 c_){ v1.set(a_); v2.set(b_); v3.set(c_); }

  template <typename Po>
  void set(const Po& other) {
    convert_from(other);
  }

  // direct access to values
  FloatType* data(){ return (FloatType*)&v1; };

  // arithmetic operators
  template <typename Po> Pos3d  operator+(const Po& other)   const{ auto out = Pos3d(other); out.v1 += v1; out.v2 += v2; out.v3 += v3;      return out; }
  template <typename Po> Pos3d  operator-(const Po& other)   const{ auto out = Pos3d(other); out.v1 -= v1; out.v2 -= v2; out.v3 -= v3;      return out; }
  template <typename Po> Pos3d& operator+=(const Po& other)  { *this = *this + other;  return *this; }
  template <typename Po> Pos3d& operator-=(const Po& other)  { *this = *this - other;  return *this; }
  Pos3d  operator*(FloatType scalar)  const{ auto out = *this; out.v1 *= scalar; out.v2 *= scalar; out.v3 *= scalar; return out; }
  Pos3d  operator/(FloatType scalar)  const{ auto out = *this; out.v1 /= scalar; out.v2 /= scalar; out.v3 /= scalar; return out; }
  Pos3d& operator*=(FloatType scalar) { *this = *this * scalar; return *this; }
  Pos3d& operator/=(FloatType scalar) { *this = *this / scalar; return *this; }

  // distance if same type
  template <typename Po, typename std::enable_if<C==Po::C && C!=LLACat, bool>::type = true>
  Units::Meters distance(const Po& other){
    auto d = *this - other;
    double x = d.v1.get();
    double y = d.v2.get();
    double z = d.v3.get();
    double dist = sqrt(x*x+y*y+z*z);
    return Units::Meters(dist);
  }

  // distance if one or both are lla
  template <typename Po, typename std::enable_if<C==LLACat||Po::C==LLACat, bool>::type = true>
  Units::Meters distance(const Po& other) {
    Pos3d<ECEFCat, Units::Meters, Units::Meters, Units::Meters> a = *this;
    Pos3d<ECEFCat, Units::Meters, Units::Meters, Units::Meters> b = other;
    a.template distance(b);
  }

  // values
  U1 v1;
  U2 v2;
  U3 v3;

private:
  static constexpr double squared(double v) noexcept {return v*v;};
  static constexpr double cubed(double v) noexcept {return v*v*v;};
  static constexpr double WGS84_A = 6378137;
  static constexpr double WGS84_F = 1/298.257223563;
  static constexpr double WGS84_B = WGS84_A*(1-WGS84_F);
  static constexpr double WGS84_E2 = 1 - squared(1-WGS84_F);

  // if same category, just call set
  template <typename Po, typename std::enable_if<C==Po::C, bool>::type = true>
  void convert_from(const Po& other){ set(other.v1,other.v2,other.v3); }

  // if other is LLA and this is ECEF, call lla2ecef
  template <typename Po, typename std::enable_if<Po::C==LLACat && C==ECEFCat, bool>::type = true>
  void convert_from(const Po& other){
    Pos3d<ECEFCat,Units::Meters,Units::Meters,Units::Meters> ecef;
    lla2ecef(other, ecef);
    *this = ecef;
  }

  // if other is ECEF and this is LLA, call ecef2lla
  template <typename Po, typename std::enable_if<Po::C==ECEFCat && C==LLACat, bool>::type = true>
  void convert_from(const Po& other){
    Pos3d<LLACat,Units::Radians,Units::Radians,Units::Meters> lla;
    ecef2lla(other, lla);
    *this = lla;
  }

  static void ecef2lla(const Pos3d<ECEFCat,Units::Meters,Units::Meters,Units::Meters>& ecef,
                Pos3d<LLACat,Units::Radians,Units::Radians,Units::Meters>& lla)
  {
    double x = ecef.v1.get();
    double y = ecef.v2.get();
    double z = ecef.v3.get();

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

  static void lla2ecef(const Pos3d<LLACat,Units::Radians,Units::Radians,Units::Meters>& lla,
                Pos3d<ECEFCat,Units::Meters,Units::Meters,Units::Meters>& ecef)
  {
    double alt = lla.v3.get();
    double clat = cos(lla.v1);
    double slat = sin(lla.v1);
    double clng = cos(lla.v2);
    double slng = sin(lla.v2);

    double N = WGS84_A / sqrt(1.0 - WGS84_E2 * slat * slat);

    ecef.v1.set((N + alt) * clat * clng);
    ecef.v2.set((N + alt) * clat * slng);
    ecef.v3.set((N * (1.0 - WGS84_E2) + alt) * slat);
  }

  // sanity checking...
  // make sure all values have same FloatType
  static_assert(std::is_same<typename U1::FloatType,typename U2::FloatType>::value,"");
  static_assert(std::is_same<typename U1::FloatType,typename U3::FloatType>::value,"");
};

template <typename U>
struct ECEF : public Pos3d<ECEFCat,U,U,U> {
  using Pos3d<ECEFCat,U,U,U>::Pos3d;
  U& x(){return this->v1;}
  U& y(){return this->v2;}
  U& z(){return this->v3;}
  const U& x() const{return this->v1;}
  const U& y() const{return this->v2;}
  const U& z() const{return this->v3;}
  static_assert(std::is_same<typename U::C, Units::Distance>::value,"");
};
using ECEF_m = ECEF<Units::Meters>;

template <typename Ull, typename Ua>
struct LLA : public Pos3d<LLACat,Ull,Ull,Ua>{
  using Pos3d<LLACat,Ull,Ull,Ua>::Pos3d;
  Ull& lat(){return this->v1;}
  Ull& lng(){return this->v2;}
  Ua& alt(){return this->v3;}
  const Ull& lat() const{return this->v1;}
  const Ull& lng() const{return this->v2;}
  const Ua& alt() const{return this->v3;}
  static_assert(std::is_same<typename Ull::C, Units::Angle>::value,"");
  static_assert(std::is_same<typename Ua::C, Units::Distance>::value,"");
};
using LLA_ddm = LLA<Units::Degrees, Units::Meters>;
using LLA_rrm = LLA<Units::Radians, Units::Meters>;

template <typename U>
struct ENU : public Pos3d<ENUCat,U,U,U> {
  using Pos3d<ENUCat,U,U,U>::Pos3d;
  template <typename P1, typename P2>
  ENU(const P1& pos, const P2& origin){
    // remove offset
    ECEF_m ecef = pos;
    ECEF_m ecefo = origin;
    Units::Meters x = ecef.x()-ecefo.x();
    Units::Meters y = ecef.y()-ecefo.y();
    Units::Meters z = ecef.z()-ecefo.z();
    // rotate
    LLA<Units::Radians,Units::Meters> lla = origin;
    double clat = cos(lla.lat());
    double slat = sin(lla.lat());
    double clng = cos(lla.lng());
    double slng = sin(lla.lng());
    this->x() = -x*slng      + y*clng;
    this->y() = -x*slat*clng - y*slat*slng + z*clat;
    this->z() =  x*clat*clng + y*clat*slng + z*slat;
  }
  template <typename P1, typename P2>
  P1 get(const P2& origin){
    // rotate
    LLA_rrm lla = origin;
    double clat = cos(lla.lat());
    double slat = sin(lla.lat());
    double clng = cos(lla.lng());
    double slng = sin(lla.lng());
    ECEF_m pos;
    pos.x() = -x()*slng - y()*slat*clng + z()*clat*clng;
    pos.y() =  x()*clng - y()*slat*slng + z()*clat*slng;
    pos.z() =           + y()*clat      + z()*slat;
    // add offset
    return pos + origin;
  }
  U& x(){return this->v1;}
  U& y(){return this->v2;}
  U& z(){return this->v3;}
  const U& x() const{return this->v1;}
  const U& y() const{return this->v2;}
  const U& z() const{return this->v3;}
  static_assert(std::is_same<typename U::C, Units::Distance>::value,"");
};
using ENU_m = ENU<Units::Meters>;

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

  // wtf
  LLA<Units::Radians, Units::Furlongs> wtf = ecef_m;
  std::cout << wtf.lat() << " " << wtf.lng() << " " << wtf.alt() << "\n";
}
