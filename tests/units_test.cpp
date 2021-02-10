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

#include "../units.h"

#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <map>

using namespace std;
using namespace Units;
using namespace Units::literals;

// if you don't know the unit at compile time, you can use something like this to convert it
Meters ToMeters(double val, const string& unit){
  if(unit=="m") return Meters(val);
  if(unit=="ft") return Feet(val);
  if(unit=="mi") return Miles(val);
  throw std::runtime_error("ToMeters: unit not recognized");
}

// parses string of format number followed by unit, spaced ignored, anything past units ignored
// e.g. "5ft", "2.7m", "  27 ft  ", "3.14ft this is ignored"
Meters ToMeters(const string& val){
  stringstream ss(val);
  double dval;
  string unit;
  ss >> dval >> unit;
  if(!ss) throw std::runtime_error("ToMeters: error parsing string: "+val);
  return ToMeters(dval, unit);
}

int main(){
  Kilometers km1(5);
  assert(km1.get<Meters>()==5000);
  assert(km1.get<Centimeters>()==500000);

  Centimeters cm1 = km1;
  assert(cm1.get<Kilometers>()==5);
  assert(cm1.get<Centimeters>()==500000);

  Kilometers km2;
  km2.set<Centimeters>(500000);
  assert(km2.get<Kilometers>()==5);
  assert(km2.get<Centimeters>()==500000);

  Centimeters cm2;
  cm2.set<Kilometers>(5);
  assert(cm2.get<Meters>()==5000);
  assert(cm2.get<Centimeters>()==500000);

  Kilometers km3;
  km3 = cm2;
  assert(km3.get<Kilometers>()==5);
  assert(km3.get<Centimeters>()==500000);

  Kilometers km4 = km3+cm2;
  assert(km4.get<Kilometers>()==10);
  assert(km4.get<Centimeters>()==1000000);

  Centimeters cm4 = km3+cm2;
  assert(cm4.get<Kilometers>()==10);
  assert(cm4.get<Centimeters>()==1000000);

  Centimeters cm5 = cm4*2;
  assert(cm5.get<Meters>()==20000);

  Centimeters cm6 = cm5/2+km4*2;
  assert(cm6 == Kilometers(30));

  double ratio = km4/cm1;
  assert(ratio==2);

  assert(cm1==km1);
  assert(cm1!=km4);
  assert(cm1<km4);
  assert(cm1<=cm1);
  assert(cm1<=km4);

  Milliseconds ms(100);
  assert(ms.get<Microseconds>()==100*1000);

  // print out conversions
  cout << "Inches in a meter: " << Meters(1).get<Inches>() << "\n";
  cout << "Inches in a foot: " << Feet(1).get<Inches>() << "\n";
  cout << "Feet in a meter: " << Meters(1).get<Feet>() << "\n";
  cout << "Yards in a meter: " << Meters(1).get<Yards>() << "\n";
  cout << "Feet in a yard: " << Yards(1).get<Feet>() << "\n";
  cout << "Meters in a mile: " << Miles(1).get<Meters>() << "\n";
  cout << "Feet in a mile: " << Miles(1).get<Feet>() << "\n";
  cout << "Yards in a mile: " << Miles(1).get<Yards>() << "\n";
  cout << "Centimeters in an inch: " << Inches(1).get<Centimeters>() << "\n";
  cout << "Meters in a nautical mile: " << NauticalMiles(1).get<Meters>() << "\n";
  cout << "120 KPH to MPH: " << KilometersPerHour(120).get<MPH>() << "\n";
  cout << "100 knots to mps: " << Knots(100).get<MetersPerSecond>() << "\n";
  cout << "Speed of sound in MPH: " << Mach(1).get<MPH>() << "\n";
  cout << "Square Feet in an Square Meter: " << SquareMeters(1).get<SquareFeet>() << "\n";
  cout << "Square Feet in an Acre: " << Acres(1).get<SquareFeet>() << "\n";
  cout << "CubicCentimeters in a Liter: " << Liters(1).get<CubicCentimeters>() << "\n";
  cout << "Liters in a Gallon: " << Gallons(1).get<Liters>() << "\n";
  cout << "Quarts in a Gallon: " << Gallons(1).get<Quarts>() << "\n";
  cout << "Pints in a Gallon: " << Gallons(1).get<Pints>() << "\n";
  cout << "Cups in a Gallon: " << Gallons(1).get<Cups>() << "\n";
  cout << "FluidOunces in a Gallon: " << Gallons(1).get<FluidOunces>() << "\n";
  cout << "CubicInches in a Gallon: " << Gallons(1).get<CubicInches>() << "\n";
  cout << "CubicInches in a Liter: " << Liters(1).get<CubicInches>() << "\n";
  cout << "SquareInches in a SquareMeter: " << SquareMeters(1).get<SquareInches>() << "\n";
  cout << "Pounds in a Kilogram: " << Kilograms(1).get<Pounds>() << "\n";
  cout << "40 C to F: " << Celsius(40).get<Fahrenheit>() << "\n";
  cout << "0 C to F: " << Celsius(0).get<Fahrenheit>() << "\n";
  cout << "-40 C to F: " << Celsius(-40).get<Fahrenheit>() << "\n";
  cout << "0 C to K: " << Celsius(0).get<Kelvin>() << "\n";
  cout << "180 deg in radians: " << Degrees(180).get<Radians>() << "\n";
  cout << "40 rods to the hogshead to miles per gallon: " << RodsToTheHogHead(40).get<MilesPerGallon>() << "\n";
  cout << "1 lightyear in meters: " << LightYears(1).get<Meters>() << "\n";

  // user-defined literals
  Meters m = 10_m + 15_cm;
  assert(m == 10.15_m);

  // Test unit conversions
  assert(1_in == 2.54_cm);
  assert(12_in == 1_ft);
  assert(1_yd == 3_ft);
  assert(381_m == 1250_ft);
  assert(1_mi == 5280_ft);
  assert(1_nm == Angstroms(10));
  assert(NauticalMiles(1)==1852_m);
  assert(Knots(900)==MetersPerSecond(463));
  assert(Mach(1)==MetersPerSecond(343));
  assert(Knots(1)==1.852_kph);
  assert(Gallons(1)==FluidOunces(128));
  assert(Gallons(1).approx(Liters(3.78541)));
  assert(Acres(1)==SquareFeet(43560));
  assert(SquareMeters(145161)==SquareFeet(1562500));
  assert(Milliliters(1)==CubicCentimeters(1));
  assert(Gallons(1)==CubicInches(231));
  assert(Liters(2048383)==CubicInches(125000000));
  assert(CubicMeters(1)==Liters(1000));
  assert(MetersPerSecondSquared(381)==FeetPerSecondSquared(1250));
  assert(45359.237_kg == 100000_lbs);
  assert(16_oz == 1_lbs);
  assert(PoundsForce(1).approx(4.44822_N));
  assert(40_C == 104_F);
  assert(0_C == 32_F);
  assert(0_C == 273.15_K);
  assert(Celsius(-40) == Fahrenheit(-40));
  assert(Rankine(0)==0_K);
  assert(Rankine(459.67)==0_F);
  assert((1_mpg).approx(KilometersPerLiter(0.425143707)));
  assert(Rods(1)==16.5_ft);
  assert(Gallons(63)==HogsHead(1));
  assert(RodsToTheHogHead(40).approx(MilesPerGallon(0.001984127)));
  assert(1000_kg == MetricTons(1));
  assert(2000_lbs == Tons(1));
  assert(MetricTons(1).approx(Tons(1.1023113)));
  assert(Bananas(1)==178_mm);
  assert(SpeedOfLight(1) == MetersPerSecond(299792458));
  assert(LightYears(1).approx(Meters(9.46053e+15)));
  assert(LightYears(10) == AstronomicalUnits(632411));
  assert(Parsecs(1) == AstronomicalUnits(206265));
  assert(86400_s == 1_dy);
  assert(220_yd == Furlong(1));
  assert(Fortnight(1)==14_dy);

  // speed = distance / time
  MPH mph = 100_mi / 8_hr;
  assert(mph == 12.5_mph);
  Knots knots = 1234.5_ft / 10_wk;

  // Energy and power
  assert(WattHours(1) == Joules(3600));
  assert(Watts(6000)*Minutes(30) == KilowattHours(3));
  assert((1.21_GW).approx(1985_MJ/(2015_in/88_mph)/1.260942));
  assert(10_V*120_A == 1.2_KW);
  assert(10_mA * 20_ohms == 200_mV);

  // frequency
  assert(1_KHz==1000_Hz);
  assert(1_MHz==1000_KHz);
  assert(1_GHz==1000_MHz);
  assert(1_THz==1000_GHz);

  // frequency to interval conversion
  assert(inverse(1_Hz)==1_s);
  assert(inverse(1_KHz)==1_ms);
  assert(inverse(1_us)==1_MHz);
  assert(inverse(1_ns)==1_GHz);

  // distance = speed * time
  Feet ft = 50_mph * 20_mn;
  cout << "50 MPH for 20 minutes is " << ft << " ft\n";
  assert(ft == 88000_ft);

  // ostream
  stringstream ss;
  ss << ft;

  // vector conversion
  vector<Feet> vf(10);
  for(int i=0;i<vf.size();i++) vf[i] = Feet(i);
  vector<Inches> vi(vf.begin(), vf.end());  // convert vector<Feet> to vector<Inches>
  for(int i=0;i<vi.size();i++) assert(vi[i].get() == vf[i].get<Inches>());

  // algorithm
  Meters mx = *std::max_element(vf.begin(), vf.end());
  Meters sum = std::accumulate(vf.begin(), vf.end(), Feet{});
  cout << "vf max and sum: " << mx << " " << sum << "\n";

  // Interunit relationships
  assert(SquareFeet(1) == 1_ft*1_ft);
  assert(SquareFeet(1) == 12_in*12_in);
  assert(Liters(1) == 10_cm*10_cm*10_cm);
  assert(CubicInches(30)/Inches(10) == SquareInches(3));
  assert(Acres(1) == 55_yd*88_yd);
  assert(CubicMeters(1000) == 10_m*10_m*10_m);
  assert(CubicInches(1000) == 10_in*10_in*10_in);
  assert((9.80665_m/1_s/1_s) == 1_G);
  assert((1_G).approx(32.174049_ft/1_s/1_s));
  assert(1_N == 1_kg*(1_m/1_s/1_s));
  assert(20_mi/10_hr == 2_mph);
  assert(20_mi == 2_mph*10_hr);
  assert(SquareInches(30)/Inches(5) == Inches(6));
  assert(200_mi/Gallons(10) == MilesPerGallon(20));
  assert(Rods(40)/HogsHead(1) == RodsToTheHogHead(40));
  assert(Radians(20)/Seconds(10) == RadiansPerSecond(2));
  assert(6050_J / 5_us == 1.21_GW);
  assert(11_mi / 450_s == 88_mph);

  // trig
  cout << "sin(180_deg): " << sin(180_deg) << "\n";
  cout << "cos(180_deg): " << cos(180_deg) << "\n";
  assert(Degrees(180)==Radians(M_PI));
  assert(sin(0_deg)==0);
  assert(sin(180_deg)==sin(M_PI));
  assert(sin(90_deg)==1);
  assert(cos(0_deg)==1);
  assert(cos(90_deg)==cos(M_PI/2));
  assert(cos(180_deg)==-1);
  assert(asin_units(0).approx(0_deg));
  assert(asin_units(1).approx(90_deg));
  assert(acos_units(0).approx(90_deg));
  assert(acos_units(1).approx(0_deg));
  assert(atan2_units(2,2).approx(45_deg));

  // override default FloatType type
  Meters mdef(5);
  assert(Float32<Inches>(24).get<Float32<Feet>>()==2);
  Float64<Meters> m64(5);
  Float32<Meters> m32(5);
  assert(sizeof(mdef)==8);
  assert(sizeof(m32)==4);
  assert(sizeof(m64)==8);
  assert(m32 == mdef);
  assert(m64 == mdef);
  Float32<Inches>in = Feet(4);
  assert(in.get()==48);
  auto mphf = Float32<Miles>(10)/Float32<Hours>(5);
  assert(sizeof(mphf)==4);
  assert(mphf.get()==2);

  // string to unit
  assert(ToMeters(5,"ft") == 5_ft);
  assert(ToMeters("5ft") == 5_ft);
  assert(ToMeters("2.7m") == 2.7_m);
  assert(ToMeters("  3.14 ft  foo ") == 3.14_ft);

  // extreme example
  assert(Angstroms(1)*Angstroms(1)*Angstroms(1)==CubicMeters(1e-30));

  // these shouldn't compile
  //Feet f = 10;            // implicit conversion not allowed
  //double d = Feet(10);    // implicit conversion not allowed
  //Miles(1).get<Knots>();
  //Miles m2;
  //m2.set<Knots>(1);
  //Miles m = Seconds(1);
  //auto s = Meters(1)/Hz(1);
}
