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

#pragma once

#include <cmath>
#include <iostream>
#include <cinttypes>

namespace Units {

///////////////////////////////////////////////////////////////////////////////
// Unit Class

template<typename Category, typename Ratio, typename FloatType_=double>
class Unit{
public:

  //////////////////////
  // template parameters
  typedef Ratio R;
  typedef Category C;
  typedef FloatType_ FloatType;
  typedef FloatType_ FT;
  static constexpr long double M = Ratio::M;
  static constexpr long double O = Ratio::O;

  //////////////////////////
  // Constructors and assign

  Unit() = default;

  // Construct from FloatTyoe
  explicit Unit(FloatType val) : m_val(val) {}

  // Construct from other Unit (and convert)
  template<typename R2,typename FT2>
  Unit(Unit<C,R2,FT2> other){
    set(other);
  }

  // Assign from other Unit (and convert)
  template<typename R2,typename FT2>
  Unit& operator=(Unit<C,R2,FT2> other){
    set(other);
    return *this;
  }

  //////////////
  // get methods

  // get value as new unit
  template<typename U>
  FloatType get() const {
    return U(*this).get();
  }

  // get value as current unit
  FloatType get() const { return m_val; }

  //////////////
  // set methods

  // set by FloatTyoe value
  template<typename U>
  void set(FloatType val){
    set(U(val));
  }

  // set by other Unit (and convert)
  template<typename R2, typename FT2>
  void set(Unit<C,R2,FT2> other){
    // Compute multiplier: mult =  R2 / R
    constexpr long double mult = ((long double)R2::M / R::M);
    // Compute offset: offset = R_offset - mult * R2_offset
    constexpr FloatType offset = ((long double)R::O) - mult*((long double)R2::O);
    // val = mult * other + offset
    m_val = (FloatType)mult * other.get();
    if(offset!=0) m_val += offset;
  }

  void set(FloatType val){
    m_val = val;
  }

  ///////////////////////
  // arithmetic operators

  // Unit = Unit +- Unit
  template<typename R2, typename FT2> Unit operator+(Unit<C,R2,FT2> other) const{ return Unit(m_val+Unit(other).get()); }
  template<typename R2, typename FT2> Unit operator-(Unit<C,R2,FT2> other) const{ return Unit(m_val-Unit(other).get()); }

  // Unit = Unit */ FloatType
  Unit operator*(FloatType scalar) const{ return Unit(m_val*scalar); }
  Unit operator/(FloatType scalar) const{ return Unit(m_val/scalar); }

  // FloatType = Unit / Unit
  template<typename R2, typename FT2> FloatType operator/(Unit<C,R2,FT2> other) const{ return m_val / Unit(other).get(); }

  // Unary operators +-
  Unit operator+() const{ return Unit(m_val); }
  Unit operator-() const{ return Unit(-m_val); }

  // operators += -= *= /=
  template<typename R2, typename FT2> Unit& operator+=(Unit<C,R2,FT2> other){ m_val += Unit(other).get(); return *this; }
  template<typename R2, typename FT2> Unit& operator-=(Unit<C,R2,FT2> other){ m_val -= Unit(other).get(); return *this; }
  Unit& operator*=(FloatType other)     { m_val *= other; return *this; }
  Unit& operator/=(FloatType other)     { m_val /= other; return *this; }

  ////////////////////
  // compare operators
  template<typename R2, typename FT2> bool operator==(Unit<C,R2,FT2> other) const{ return m_val == Unit(other).get(); }
  template<typename R2, typename FT2> bool operator!=(Unit<C,R2,FT2> other) const{ return m_val != Unit(other).get(); }
  template<typename R2, typename FT2> bool operator<=(Unit<C,R2,FT2> other) const{ return m_val <= Unit(other).get(); }
  template<typename R2, typename FT2> bool operator>=(Unit<C,R2,FT2> other) const{ return m_val >= Unit(other).get(); }
  template<typename R2, typename FT2> bool operator< (Unit<C,R2,FT2> other) const{ return m_val <  Unit(other).get(); }
  template<typename R2, typename FT2> bool operator> (Unit<C,R2,FT2> other) const{ return m_val >  Unit(other).get(); }

  template<typename R2, typename FT2>
  bool approx(Unit<C,R2,FT2> other, double epsilon=1e-6) const{
    return fabs(m_val-Unit(other).get()) <= epsilon;
  }

  // allow explicit casts to FloatType (but not implicit)
  explicit operator FloatType() const {return m_val;};

private:
  // storage
  FloatType m_val = 0;
};

///////////////////////////////////////////////////////////////////////////////
// Override default FloatType.  e.g. Float32<Meters> is float instead of double
template<typename U> using Float32 = Unit<typename U::C, typename U::R, float>;
template<typename U> using Float64 = Unit<typename U::C, typename U::R, double>;

///////////////////////////////////////////////////////////////////////////////
// Ratio template

// Ratio (and offset, if needed)
template<int64_t N, int64_t D=1, int64_t ON=0, int64_t OD=1>
struct Ratio {
  static constexpr long double M = (long double)N/D;
  static constexpr long double O = (long double)ON/OD;
};

// RatioProd base case
template<typename ... T>
struct RatioProd {
  static constexpr long double M = 1;
  static constexpr long double O = 0;
};

// RatioProd recursive
template<typename R, typename ... Rest>
struct RatioProd<R, Rest ...>{
  static_assert(R::O==0,"RatioProd doesn't support offset");
  static constexpr long double M = R::M*RatioProd<Rest ...>::M;
  static constexpr long double O = 0;
};

// RatioInv
template<typename R>
struct RatioInv{
  static_assert(R::O==0,"RatioInv doesn't support offset");
  static constexpr long double M = 1/R::M;
  static constexpr long double O = 0;
};

///////////////////////////////////////////////////////////////////////////////
// Macros

// helper macro for generating user-defined literals for units
#define GEN_LITERAL(lit, unit)                                              \
  namespace literals {                                                      \
  inline unit operator "" lit(long double v){return unit(v);};              \
  inline unit operator "" lit(unsigned long long v){return unit(v);};       \
  }

// Macro for generating Uo = Ua * Ub
#define GEN_MULT(Uo,Ua,Ub,FT)                                                \
  template <typename Ra, typename Rb,     /* Ra*Rb*Uo/Ua/Ub */               \
            typename Ro = RatioProd<Ra,Rb,Uo,RatioInv<Ua>,RatioInv<Ub>>>     \
  inline Units::Unit<Uo::C, Ro, FT>                                          \
  operator*(Units::Unit<Ua::C,Ra,FT> a, Units::Unit<Ub::C,Rb,FT> b) {        \
    static_assert(Ra::O==0||Rb::O==0, "GEN_MULT doesn't support offset");    \
    Uo::FloatType val = a.get()*b.get();                                     \
    return Units::Unit<Uo::C, Ro, FT>(val);                                  \
  }

// Macro for generating Uo = Ua / Ub
#define GEN_DIV(Uo,Ua,Ub,FT)                                                 \
  template <typename Ra, typename Rb,     /* Ra/Rb*Uo/Ua*Ub */               \
            typename Ro = RatioProd<Ra,RatioInv<Rb>,Uo,RatioInv<Ua>,Ub>>     \
  inline Units::Unit<Uo::C, Ro, FT>                                          \
  operator/(Units::Unit<Ua::C,Ra,FT> a, Units::Unit<Ub::C,Rb,FT> b) {        \
    static_assert(Ra::O==0||Rb::O==0, "GEN_DIV doesn't support offset");     \
    Uo::FloatType val = a.get()/b.get();                                     \
    return Units::Unit<Uo::C, Ro, FT>(val);                                  \
  }

// Macro for generating Uo = sqrt(Ui) (Note: potentially unnecessary conversions because no constexpr sqrt)
#define GEN_SQRT(Uo,Ui)                                              \
  Uo sqrt(Ui in){ return Uo(std::sqrt(in.get())); }                  \

// generate functions U1=U2*U2, U1=U2/U1, and U2=sqrt(U1)
#define GEN_MULT_DIV_SQ(U1,U2)                                       \
  GEN_MULT(U1,U2,U2,double)                                          \
  GEN_MULT(U1,U2,U2,float)                                           \
  GEN_DIV(U2,U1,U2,double)                                           \
  GEN_DIV(U2,U1,U2,float)                                            \
  GEN_SQRT(U2,U1)

// generate functions: U1=U2*U3, U1=U3*U2, U2=U1/U3, amd U3=U1/U2
#define GEN_MULT_DIV(U1,U2,U3)                                       \
  GEN_MULT(U1,U2,U3,double)                                          \
  GEN_MULT(U1,U2,U3,float)                                           \
  GEN_MULT(U1,U3,U2,double)                                          \
  GEN_MULT(U1,U3,U2,float)                                           \
  GEN_DIV(U2,U1,U3,double)                                           \
  GEN_DIV(U2,U1,U3,float)                                            \
  GEN_DIV(U3,U1,U2,double)                                           \
  GEN_DIV(U3,U1,U2,float)

// Macro helper GEN_INVERSE
#define GEN_INVERSE_HELPER(Uo,Ui,FT)                                                 \
  template <typename Ri,                                                             \
            typename Ro = RatioProd<RatioInv<Ri>,Uo,Ui>>                             \
  inline Unit<Uo::C, Ro, FT> inverse(Unit<Ui::C, Ri, FT> s){                         \
    static_assert(Uo::O==0||Ui::O==0, "GEN_INVERSE doesn't support offset");         \
    return Unit<Uo::C, Ro, FT>(1/s.get());                                           \
  }

// Macro helper for generating Uo = inverse(Ui)
#define GEN_INVERSE(U1, U2)           \
  GEN_INVERSE_HELPER(U1, U2, double)  \
  GEN_INVERSE_HELPER(U1, U2, float)   \
  GEN_INVERSE_HELPER(U2, U1, double)  \
  GEN_INVERSE_HELPER(U2, U1, float)

///////////////////////////////////////////////////////////////////////////////
// categories
struct Distance{};
struct Area{};
struct Volume{};
struct Time{};
struct Speed{};
struct Acceleration{};
struct Mass{};
struct Force{};
struct Momentum{};
struct Frequency{};
struct Energy{};
struct Power{};
struct Voltage{};
struct Amperage{};
struct Resistance{};
struct Angle{};
struct AngularVelocity{};
struct Temperature{};
struct FuelEfficiency{};

// SI helper
using Tera  = Ratio<1000000000000>;
using Giga  = Ratio<1000000000>;
using Mega  = Ratio<1000000>;
using Kilo  = Ratio<1000>;
using Base  = Ratio<1,1>;
using Centi = Ratio<1,100>;
using Milli = Ratio<1,1000>;
using Micro = Ratio<1,1000000>;
using Nano  = Ratio<1,1000000000>;

///////////
// Distance
using Kilometers    = Unit<Distance,Kilo>;
using Meters        = Unit<Distance,Base>;
using Centimeters   = Unit<Distance,Centi>;
using Millimeters   = Unit<Distance,Milli>;
using Microns       = Unit<Distance,Micro>;
using Nanometers    = Unit<Distance,Nano>;
using Angstroms     = Unit<Distance,RatioProd<Nanometers,Ratio<1,10>>>;
using Inches        = Unit<Distance,Ratio<254,10000>>;
using Feet          = Unit<Distance,Ratio<254*12,10000>>;
using Yards         = Unit<Distance,Ratio<254*12*3,10000>>;
using Furlongs      = Unit<Distance,Ratio<254*12*3*220,10000>>;
using Miles         = Unit<Distance,Ratio<254*12*3*1760,10000>>;
using NauticalMiles = Unit<Distance,Ratio<1852>>;
using Rods          = Unit<Distance,RatioProd<Feet,Ratio<165,10>>>;
using Bananas       = Unit<Distance,Ratio<178,1000>>;

GEN_LITERAL(_nm, Nanometers)
GEN_LITERAL(_um, Microns)
GEN_LITERAL(_mm, Millimeters)
GEN_LITERAL(_cm, Centimeters)
GEN_LITERAL( _m, Meters)
GEN_LITERAL(_km, Kilometers)
GEN_LITERAL(_in, Inches)
GEN_LITERAL(_ft, Feet)
GEN_LITERAL(_yd, Yards)
GEN_LITERAL(_mi, Miles)

///////
// Area
using SquareMeters     = Unit<Area,Base>;
using SquareInches     = Unit<Area,RatioProd<Inches,Inches>>;
using SquareFeet       = Unit<Area,RatioProd<Feet,Feet>>;
using SquareMiles      = Unit<Area,RatioProd<Miles,Miles>>;
using SquareKilometers = Unit<Area,RatioProd<Kilometers,Kilometers>>;
using Acres            = Unit<Area,RatioProd<SquareFeet,Ratio<43560>>>;

/////////
// Volume
using Liters           = Unit<Volume,Base>;
using Milliliters      = Unit<Volume,Milli>;
using CubicCentimeters = Unit<Volume,Milli>;
using Microliters      = Unit<Volume,Micro>;
using CubicMeters      = Unit<Volume,Ratio<1000>>;
using CubicInches      = Unit<Volume,RatioProd<Inches,Inches,Inches,CubicMeters>>;
using Gallons          = Unit<Volume,RatioProd<CubicInches,Ratio<231>>>;
using Quarts           = Unit<Volume,RatioProd<Gallons,Ratio<1,4>>>;
using Pints            = Unit<Volume,RatioProd<Gallons,Ratio<1,8>>>;
using Cups             = Unit<Volume,RatioProd<Gallons,Ratio<1,16>>>;
using FluidOunces      = Unit<Volume,RatioProd<Gallons,Ratio<1,128>>>;
using HogsHead         = Unit<Volume,RatioProd<Gallons,Ratio<63>>>;

GEN_LITERAL( _l, Liters)
GEN_LITERAL(_ml, Milliliters)
GEN_LITERAL(_ul, Microliters)
GEN_LITERAL(_cc, CubicCentimeters)

///////
// Time
using Seconds      = Unit<Time,Base>;
using Minutes      = Unit<Time,Ratio<60>>;
using Hours        = Unit<Time,Ratio<60*60>>;
using Days         = Unit<Time,Ratio<60*60*24>>;
using Weeks        = Unit<Time,Ratio<60*60*24*7>>;
using Fortnights   = Unit<Time,Ratio<60*60*24*14>>;
using Years        = Unit<Time,RatioProd<Days,Ratio<3652422,10000>>>;
using Milliseconds = Unit<Time,Milli>;
using Microseconds = Unit<Time,Micro>;
using Nanoseconds  = Unit<Time,Nano>;

GEN_LITERAL(_wk, Weeks)
GEN_LITERAL(_dy, Days)
GEN_LITERAL(_hr, Hours)
GEN_LITERAL(_mn, Minutes)
GEN_LITERAL( _s, Seconds)
GEN_LITERAL(_ms, Milliseconds)
GEN_LITERAL(_us, Microseconds)
GEN_LITERAL(_ns, Nanoseconds)

////////
// Speed
using MetersPerSecond   = Unit<Speed,Base>;
using KilometersPerHour = Unit<Speed,Ratio<1000,3600>>;
using MilesPerHour      = Unit<Speed,RatioProd<Miles,RatioInv<Hours>>>;
using MPH               = Unit<Speed,MetersPerSecond>;
using Knots             = Unit<Speed,Ratio<1852,3600>>;
using Mach              = Unit<Speed,Ratio<343>>;
using SpeedOfLight      = Unit<Speed,Ratio<299792458>>;

GEN_LITERAL(_mps, MetersPerSecond)
GEN_LITERAL(_mph, MilesPerHour)
GEN_LITERAL(_kph, KilometersPerHour)

/////////////////////////
// Astronomical Distances
using LightYears        = Unit<Distance,RatioProd<SpeedOfLight,Years>>;
using AstronomicalUnits = Unit<Distance,RatioProd<LightYears,Ratio<10,632411>>>;
using Parsecs           = Unit<Distance,RatioProd<AstronomicalUnits,Ratio<206265>>>;

///////////////
// Acceleration
using MetersPerSecondSquared = Unit<Acceleration,Base>;
using FeetPerSecondSquared   = Unit<Acceleration,Feet>;
using Gravity                = Unit<Acceleration,Ratio<980665,100000>>;

GEN_LITERAL(_G, Gravity)

///////
// Mass
using MetricTons = Unit<Mass,Mega>;
using Kilograms  = Unit<Mass,Kilo>;
using Grams      = Unit<Mass,Base>;
using Milligrams = Unit<Mass,Milli>;
using Micrograms = Unit<Mass,Micro>;
using Pounds     = Unit<Mass,Ratio<45359237,100000>>;
using Ounces     = Unit<Mass,RatioProd<Pounds,Ratio<1,16>>>;
using Tons       = Unit<Mass,RatioProd<Pounds,Ratio<2000>>>;
using Stones     = Unit<Mass,RatioProd<Pounds,Ratio<14>>>;

GEN_LITERAL(_kg, Kilograms)
GEN_LITERAL( _g, Grams)
GEN_LITERAL(_mg, Milligrams)
GEN_LITERAL(_ug, Micrograms)
GEN_LITERAL(_lbs, Pounds)
GEN_LITERAL(_oz, Ounces)

/////////////////////
// Force and momentum
using Newtons                  = Unit<Force,Base>;
using PoundsForce              = Unit<Force, RatioProd<Pounds,Gravity,Ratio<1,1000>>>;
using KilogramsMetersPerSecond = Unit<Momentum,Base>;
GEN_LITERAL(_N, Newtons);

/////////
// Energy
using Gigawatts  = Unit<Power,Giga>;
using Megawatts  = Unit<Power,Mega>;
using Kilowatts  = Unit<Power,Kilo>;
using Watts      = Unit<Power,Base>;
using Milliwatts = Unit<Power,Milli>;

GEN_LITERAL(_GW, Gigawatts);
GEN_LITERAL(_MW, Megawatts);
GEN_LITERAL(_KW, Kilowatts);
GEN_LITERAL( _W, Watts);
GEN_LITERAL(_mW, Milliwatts);

////////
// Power
using Megajoules    = Unit<Energy,Mega>;
using Kilojoules    = Unit<Energy,Kilo>;
using Joules        = Unit<Energy,Base>;
using WattHours     = Unit<Energy,Hours>;
using KilowattHours = Unit<Energy,RatioProd<Kilo,WattHours>>;
using MegawattHours = Unit<Energy,RatioProd<Mega,WattHours>>;

GEN_LITERAL(_MWh, MegawattHours)
GEN_LITERAL(_KWh, KilowattHours)
GEN_LITERAL( _Wh, WattHours)
GEN_LITERAL(_MJ,  Megajoules)
GEN_LITERAL(_KJ,  Kilojoules)
GEN_LITERAL( _J,  Joules)

//////////////
// Electricity
using Kilovolts  = Unit<Voltage,Kilo>;
using Volts      = Unit<Voltage,Base>;
using Millivolts = Unit<Voltage,Milli>;
using Amps       = Unit<Amperage,Base>;
using Milliamps  = Unit<Amperage,Milli>;
using Kiloohms   = Unit<Resistance,Kilo>;
using Ohms       = Unit<Resistance,Base>;
using Milliohms  = Unit<Resistance,Milli>;

GEN_LITERAL(_KV, Kilovolts)
GEN_LITERAL(_V, Volts)
GEN_LITERAL(_mV, Millivolts)
GEN_LITERAL(_A, Amps)
GEN_LITERAL(_mA, Milliamps)
GEN_LITERAL(_ohms, Ohms)

//////////////////
// Fuel Efficiency
using KilometersPerLiter = Unit<FuelEfficiency, Base>;
using MilesPerGallon     = Unit<FuelEfficiency, RatioProd<Miles, RatioInv<Gallons>, Ratio<1,1000>>>;
using RodsToTheHogHead   = Unit<FuelEfficiency, RatioProd<Rods, RatioInv<HogsHead>, Ratio<1,1000>>>;

GEN_LITERAL(_mpg, MilesPerGallon);

////////////
// Frequency
using  Hz = Unit<Frequency,Base>;
using KHz = Unit<Frequency,Kilo>;
using MHz = Unit<Frequency,Mega>;
using GHz = Unit<Frequency,Giga>;
using THz = Unit<Frequency,Tera>;

GEN_LITERAL( _Hz,  Hz)
GEN_LITERAL(_KHz, KHz)
GEN_LITERAL(_MHz, MHz)
GEN_LITERAL(_GHz, GHz)
GEN_LITERAL(_THz, THz)

////////
// Angle
struct RatioPi180 {
  static constexpr double M = M_PIl/180;
  static constexpr double O = 0;
};
using Radians = Unit<Angle,Base>;
using Degrees = Unit<Angle,RatioPi180>;
using RadiansPerSecond = Unit<AngularVelocity,Base>;
using DegreesPerSecond = Unit<AngularVelocity,RatioPi180>;
GEN_LITERAL(_rad, Radians)
GEN_LITERAL(_deg, Degrees)

// trig functions
double sin(Radians r) { return std::sin(r.get()); }
double cos(Radians r) { return std::cos(r.get()); }
double tan(Radians r) { return std::tan(r.get()); }
// Note: using _units here to avoid collision with std namespace since they both take double as arg
Radians asin_units(double v) { return Radians(std::asin(v)); }
Radians acos_units(double v) { return Radians(std::acos(v)); }
Radians atan_units(double v) { return Radians(std::atan(v)); }
Radians atan2_units(double v1, double v2) { return Radians(std::atan2(v1,v2)); }

///////////////
// Temperature
using Kelvin     = Unit<Temperature,Ratio<1,1,      0,  1>>;
using Celsius    = Unit<Temperature,Ratio<1,1, -27315,100>>;
using Fahrenheit = Unit<Temperature,Ratio<5,9, -45967,100>>;
using Rankine    = Unit<Temperature,Ratio<5,9,      0,  1>>;

GEN_LITERAL(_C, Celsius)
GEN_LITERAL(_K, Kelvin)
GEN_LITERAL(_F, Fahrenheit)

///////////////////////////////////////////////////////////////////////////////
// Cross Unit Relationships

// Frequency = inverse(Time) and vice versa
GEN_INVERSE(Hz, Seconds)

// Generate operator* and operator/ functions for interunit math
GEN_MULT_DIV(Meters, MetersPerSecond, Seconds)                      // Distance = Speed*Time
GEN_MULT_DIV(MetersPerSecond, MetersPerSecondSquared, Seconds)      // Speed = Acceleration*Time
GEN_MULT_DIV(Newtons, Kilograms, MetersPerSecondSquared)            // Force = Mass*Acceleration
GEN_MULT_DIV(KilogramsMetersPerSecond, Kilograms, MetersPerSecond)  // Momentum = Mass*Velocity
GEN_MULT_DIV(KilogramsMetersPerSecond, Newtons, Seconds)            // Momentum = Force*Time
GEN_MULT_DIV(CubicMeters, SquareMeters, Meters)                     // Volume = Area*Distance
GEN_MULT_DIV(WattHours, Watts, Hours)                               // Energy = Power*Time
GEN_MULT_DIV(Volts, Amps, Ohms)                                     // Ohm's Law: V=I*R
GEN_MULT_DIV(Watts, Volts, Amps)                                    // Watts = Volts * Amps
GEN_MULT_DIV(Miles, MilesPerGallon, Gallons)                        // Fuel Efficiency
GEN_MULT_DIV(Radians, RadiansPerSecond, Seconds)                    // Angular Velocity

// Generate operator* and operator/ and sqrt()
GEN_MULT_DIV_SQ(SquareMeters, Meters)                               // Area = Dist*Dist

///////////////////////////////////////////////////////////////////////////////
// Stream operator
template <typename R,typename C,typename FloatType>
inline std::ostream& operator<<(std::ostream&os, const Unit<R,C,FloatType>& u){
  os << u.get();
  return os;
}

}  // end namespace Units
