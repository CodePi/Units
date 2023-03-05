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

#if __cplusplus < 202000
#error "C++20 compatible compiler is required for units20.h, use units.h instead"
#endif

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <cinttypes>

namespace Units {

///////////////////////////////////////////////////////////////////////////////
// Unit Class

template<typename Category, typename Conversion, typename FloatType_=double>
class Unit{
public:

  //////////////////////
  // template parameters
  typedef Conversion R;
  typedef Category C;
  typedef FloatType_ FloatType;
  typedef FloatType_ FT;
  static constexpr long double M = Conversion::M;
  static constexpr long double O = Conversion::O;

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
// Conversion template

// Conversion (and offset, if needed)
template<long double Multiplier, long double Offset=0.0L>
struct Conversion {
  static constexpr long double M = Multiplier;
  static constexpr long double O = Offset;
};

// constexpr sqrt
inline constexpr long double sqrtc(long double x, long double curr, long double prev){
  return curr == prev ? curr : sqrtc(x,0.5*(curr+x/curr),curr);
}
inline constexpr long double sqrtc(long double x){ return sqrtc(x,x,0); }
inline constexpr long double squarec(long double x){ return x*x; }
inline constexpr long double cubec(long double x){ return x*x*x; }


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
            typename Ro = Conversion<Ra::M*Rb::M*Uo::M/Ua::M/Ub::M>>              \
  inline Units::Unit<Uo::C, Ro, FT>                                          \
  operator*(Units::Unit<Ua::C,Ra,FT> a, Units::Unit<Ub::C,Rb,FT> b) {        \
    static_assert(Ra::O==0||Rb::O==0, "GEN_MULT doesn't support offset");    \
    Uo::FloatType val = a.get()*b.get();                                     \
    return Units::Unit<Uo::C, Ro, FT>(val);                                  \
  }

// Macro for generating Uo = Ua / Ub
#define GEN_DIV(Uo,Ua,Ub,FT)                                                 \
  template <typename Ra, typename Rb,     /* Ra/Rb*Uo/Ua*Ub */               \
            typename Ro = Conversion<Ra::M/Rb::M*Uo::M/Ua::M*Ub::M>>     \
  inline Units::Unit<Uo::C, Ro, FT>                                          \
  operator/(Units::Unit<Ua::C,Ra,FT> a, Units::Unit<Ub::C,Rb,FT> b) {        \
    static_assert(Ra::O==0||Rb::O==0, "GEN_DIV doesn't support offset");     \
    Uo::FloatType val = a.get()/b.get();                                     \
    return Units::Unit<Uo::C, Ro, FT>(val);                                  \
  }

// Macro for generating Uo = sqrt(Ui)
#define GEN_SQRT(Uo,Ui,FT)                                                     \
  template <typename Ri,          /* Uo*sqrt(Ri/Ui) */                         \
            typename Ro = Conversion<Uo::M*sqrt(Ri::M/Ui::M)>> \
  inline Units::Unit<Uo::C, Ro, FT>                                            \
  sqrt(Units::Unit<Ui::C,Ri,FT> in) {                                          \
    static_assert(Ri::O==0, "GEN_SQRT doesn't support offset");                \
    Uo::FloatType val = std::sqrt(in.get());                                   \
    return Units::Unit<Uo::C, Ro, FT>(val);                                    \
  }

// generate functions U1=U2*U2, U1=U2/U1, and U2=sqrt(U1)
#define GEN_MULT_DIV_SQ(U1,U2)                                       \
  GEN_MULT(U1,U2,U2,double)                                          \
  GEN_MULT(U1,U2,U2,float)                                           \
  GEN_DIV(U2,U1,U2,double)                                           \
  GEN_DIV(U2,U1,U2,float)                                            \
  GEN_SQRT(U2,U1,double)                                             \
  GEN_SQRT(U2,U1,float)

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
            typename Ro = Conversion<Ui::M*Uo::M/Ri::M>>                             \
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
struct AcceleConversionn{};
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
using Tera  = Conversion<1e12L>;
using Giga  = Conversion<1e9L>;
using Mega  = Conversion<1e6L>;
using Kilo  = Conversion<1e3L>;
using Base  = Conversion<1e0L>;
using Centi = Conversion<1e-2L>;
using Milli = Conversion<1e-3L>;
using Micro = Conversion<1e-6L>;
using Nano  = Conversion<1e-9L>;

///////////
// Distance
using Kilometers    = Unit<Distance,Kilo>;
using Meters        = Unit<Distance,Base>;
using Centimeters   = Unit<Distance,Centi>;
using Millimeters   = Unit<Distance,Milli>;
using Microns       = Unit<Distance,Micro>;
using Nanometers    = Unit<Distance,Nano>;
using Angstroms     = Unit<Distance,Conversion<Nanometers::M/10>>;
using Inches        = Unit<Distance,Conversion<2.54L/100>>;
using Feet          = Unit<Distance,Conversion<Inches::M*12>>;
using Yards         = Unit<Distance,Conversion<Inches::M*12*3>>;
using Furlongs      = Unit<Distance,Conversion<Inches::M*12*3*220>>;
using Miles         = Unit<Distance,Conversion<Inches::M*12*3*1760>>;
using NauticalMiles = Unit<Distance,Conversion<1852.0L>>;
using Rods          = Unit<Distance,Conversion<Feet::M*16.5>>;
using Bananas       = Unit<Distance,Conversion<0.178L>>;

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
using SquareInches     = Unit<Area,Conversion<Inches::M*Inches::M>>;
using SquareFeet       = Unit<Area,Conversion<Feet::M*Feet::M>>;
using SquareMiles      = Unit<Area,Conversion<Miles::M*Miles::M>>;
using SquareKilometers = Unit<Area,Conversion<Kilometers::M*Kilometers::M>>;
using Acres            = Unit<Area,Conversion<SquareFeet::M*43560>>;

/////////
// Volume
using Liters           = Unit<Volume,Base>;
using Milliliters      = Unit<Volume,Milli>;
using CubicCentimeters = Unit<Volume,Milli>;
using Microliters      = Unit<Volume,Micro>;
using CubicMeters      = Unit<Volume,Conversion<1000.0L>>;
using CubicInches      = Unit<Volume,Conversion<Inches::M*Inches::M*Inches::M*CubicMeters::M>>;
using Gallons          = Unit<Volume,Conversion<CubicInches::M*231>>;
using Quarts           = Unit<Volume,Conversion<Gallons::M/4>>;
using Pints            = Unit<Volume,Conversion<Gallons::M/8>>;
using Cups             = Unit<Volume,Conversion<Gallons::M/16>>;
using FluidOunces      = Unit<Volume,Conversion<Gallons::M/128>>;
using HogsHead         = Unit<Volume,Conversion<Gallons::M*63>>;

GEN_LITERAL( _l, Liters)
GEN_LITERAL(_ml, Milliliters)
GEN_LITERAL(_ul, Microliters)
GEN_LITERAL(_cc, CubicCentimeters)

///////
// Time
using Seconds      = Unit<Time,Base>;
using Minutes      = Unit<Time,Conversion<60.0L>>;
using Hours        = Unit<Time,Conversion<60.0L*60>>;
using Days         = Unit<Time,Conversion<60.0L*60*24>>;
using Weeks        = Unit<Time,Conversion<60.0L*60*24*7>>;
using Fortnights   = Unit<Time,Conversion<60.0L*60*24*14>>;
using Years        = Unit<Time,Conversion<Days::M*365.2422>>;
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
using KilometersPerHour = Unit<Speed,Conversion<1000.0L/3600>>;
using MilesPerHour      = Unit<Speed,Conversion<Miles::M/Hours::M>>;
using MPH               = Unit<Speed,Conversion<MetersPerSecond::M>>;
using Knots             = Unit<Speed,Conversion<1852.0L/3600>>;
using Mach              = Unit<Speed,Conversion<343.0L>>;
using SpeedOfLight      = Unit<Speed,Conversion<299792458.0L>>;

GEN_LITERAL(_mps, MetersPerSecond)
GEN_LITERAL(_mph, MilesPerHour)
GEN_LITERAL(_kph, KilometersPerHour)

/////////////////////////
// Astronomical Distances
using LightYears        = Unit<Distance,Conversion<SpeedOfLight::M*Years::M>>;
using AstronomicalUnits = Unit<Distance,Conversion<LightYears::M/63241.1>>;
using Parsecs           = Unit<Distance,Conversion<AstronomicalUnits::M*206265>>;

///////////////
// AcceleConversionn
using MetersPerSecondSquared = Unit<AcceleConversionn,Base>;
using FeetPerSecondSquared   = Unit<AcceleConversionn,Feet>;
using Gravity                = Unit<AcceleConversionn,Conversion<9.80665L>>;

GEN_LITERAL(_G, Gravity)

///////
// Mass
using MetricTons = Unit<Mass,Mega>;
using Kilograms  = Unit<Mass,Kilo>;
using Grams      = Unit<Mass,Base>;
using Milligrams = Unit<Mass,Milli>;
using Micrograms = Unit<Mass,Micro>;
using Pounds     = Unit<Mass,Conversion<453.59237L>>;
using Ounces     = Unit<Mass,Conversion<Pounds::M/16>>;
using Tons       = Unit<Mass,Conversion<Pounds::M*2000>>;
using Stones     = Unit<Mass,Conversion<Pounds::M*14>>;

GEN_LITERAL(_kg, Kilograms)
GEN_LITERAL( _g, Grams)
GEN_LITERAL(_mg, Milligrams)
GEN_LITERAL(_ug, Micrograms)
GEN_LITERAL(_lbs, Pounds)
GEN_LITERAL(_oz, Ounces)

/////////////////////
// Force and momentum
using Newtons                  = Unit<Force,Base>;
using PoundsForce              = Unit<Force, Conversion<Pounds::M*Gravity::M/1000>>;
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
using KilowattHours = Unit<Energy,Conversion<Kilo::M,WattHours::M>>;
using MegawattHours = Unit<Energy,Conversion<Mega::M,WattHours::M>>;

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
using MilesPerGallon     = Unit<FuelEfficiency, Conversion<Miles::M/Gallons::M/1000>>;
using RodsToTheHogHead   = Unit<FuelEfficiency, Conversion<Rods::M/HogsHead::M/1000>>;

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
#ifndef M_PIl  // not defined on some systems
#define M_PIl M_PI
#endif
struct ConversionPi180 {
  static constexpr long double M = M_PIl/180;
  static constexpr long double O = 0;
};
using Radians = Unit<Angle,Base>;
using Degrees = Unit<Angle,ConversionPi180>;
using RadiansPerSecond = Unit<AngularVelocity,Base>;
using DegreesPerSecond = Unit<AngularVelocity,ConversionPi180>;
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
using Kelvin     = Unit<Temperature,Conversion<1.0L, 0.0L>>;
using Celsius    = Unit<Temperature,Conversion<1.0L, -273.15L>>;
using Fahrenheit = Unit<Temperature,Conversion<5.0L/9, -459.67L>>;
using Rankine    = Unit<Temperature,Conversion<5.0L/9, 0.0L>>;

GEN_LITERAL(_C, Celsius)
GEN_LITERAL(_K, Kelvin)
GEN_LITERAL(_F, Fahrenheit)

///////////////////////////////////////////////////////////////////////////////
// Cross Unit Relationships

// Frequency = inverse(Time) and vice versa
GEN_INVERSE(Hz, Seconds)

// Generate operator* and operator/ functions for interunit math
GEN_MULT_DIV(Meters, MetersPerSecond, Seconds)                      // Distance = Speed*Time
GEN_MULT_DIV(MetersPerSecond, MetersPerSecondSquared, Seconds)      // Speed = AcceleConversionn*Time
GEN_MULT_DIV(Newtons, Kilograms, MetersPerSecondSquared)            // Force = Mass*AcceleConversionn
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
