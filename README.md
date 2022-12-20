# Units
Header only, C++11, minimal overhead, type enforced conversion and arithmetic library for measurement units.

## Features
* Unit safety enforced at compile time. 
* Most conversion operations compile down to single assembly instruction.
* Zero size overhead.  The unit storage is exactly the size of the floating point type that represents it (float or double).
* Supports more than a dozen measurement categories.  E.g. distance, speed, temperature, etc..
* Supports more than 100 varied [units](units.md).
* Easy to add new units and categories.  Typically one line of code per new unit.
* Automatically converts accurately and quickly between units in the same category.  E.g. meters to miles.
* Understands [relationships](units.md#relationships) between categories.  E.g. `speed = distance/time` and `energy = power*time`.
* Most units have user-defined literals.  E.g. `88_mph` and `1.21_GW`.
* Support for common 3D positional coordinate systems in units_coords.h (LLA, ECEF, ENU).
* Optional (experimental) support for Eigen library.  E.g. can create an eigen vector of type Meters.  And supports bulk conversion and arithmetic on eigen containers of units.

## Example Usage
### Basic usage and math with automatic conversion:
```c++
Meters d1(5.5);         // d1 is 5.1 meters
Centimeters d2(27);     // d2 is 27 centimeters
Millimeters d3 = d1+d2; // d3 is 5370 millimeters
auto d4 = Meters(5)*7;  // d4 is 35 meters (unit*double is unit)
Parsecs kessel_run(12); // kessel_run in 12 parsecs  
assert(Fahrenheit(212)==Celsius(100));  // Boiling point
```

### Built-in user-defined literals
```c++
Meters d = 25.3_cm + 1_km;   // d is 1000.253 meters
Celsius t = 98.6_F;          // t is 37 celsius
```

### Interunit relationships
```c++
MilesPerHour speed = 11_mi / 450_s;  // 11 miles in 450 seconds is 88 mph
Gigawatts power = 6050_J / 5_us;     // 6050 joules in 5 microseconds is 1.21 gigawatts
Millivolts mv = 10_mA * 20_ohms;     // 200 millivolts (Ohm's Law)
```

### Get/Set methods
```c++
Meter u1(5);
double d1 = m.get();                  // get raw value 5
double d2 = m.get<Inches>();          // get value converted to inches (196.85)
Inches u2;        
u2.set(u1);                           // set unit from other.  u2 becomes 196.85 inches
u2.set<Feet>(3);                      // set unit from double with unit template param.  u2 becomes 36 inches.
double d3 = Feet(100).get<Meters>();  // convert double from one unit to another (d3 becomes 30.48)
```

### Adding new units
New units are easy to add.  For existing categories, such as `Distance`, say you want to add the unit `Smoots`.  A Smoot
is exactly 1.702m.  Meters are the base unit for Distance.  So to add Smoots, you need a ratio that when multiplied
by a value in Smoots, you get the equivalent value in Meters.  Use the following format: 
```c++
using Smoots = Unit<Distance, Ratio<1702,1000>>;
```
Alternatively, if the most similar unit is not the base unit, you can combine multiple ratios.  Say you wanted to add
a new unit called `FootballFields`.  A football field is exactly 100 yards.  So you can combine the ratio for yards and
100 using `RatioProd`.  I.e.:
```c++
using FootballFields = Unit<Distance, RatioProd<Yards, Ratio<100,1> >>;
```

## More info
* See [units.md](units.md) for all of the pre-defined units, categories, and relationships.
* See [tests/units_test.cpp](tests/units_test.cpp) for many more examples.
* See [tests/units_coords_test.cpp](tests/units_coords_test.cpp) for LLA, ECEF, and ENU usage examples.
