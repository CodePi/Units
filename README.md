# Units
Header only, C++11, minimal overhead, type enforced conversion library for measurement units.

## Features
* Units safety enforced at compile time. 
* Most conversion operations compile down to single assembly instruction.
* Zero size overhead.  The unit storage is exactly the size of the floating point type that represents it (float or double).
* Supports more than a dozen measurement categories.  E.g. distance, speed, temperature, etc..
* Supports more than 100 varied units.
* Easy to add new units and categories.  Typically one line of code per new unit.
* Automatically converts accurately and quickly between units in the same category.  E.g. meters to miles.
* Understands relationships between categories.  E.g. `speed = distance/time` and `energy = power*time`.
* Most units have user-defined literals.  E.g. `88_mph` and `1.21_GW`. 
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
Gigawatts power = 6050_J / 5_us;     // 60500 joules in 50 microseconds is 1.21 gigawatts
Millivolts mv = 10_mA * 20_ohms;     // 200 millivolts (Ohm's Law)
```

## More info
* See [units.md](units.md) for all of the pre-defined units, categories, and relationships.
* See [tests/units_test.cpp](tests/units_test.cpp) for many more examples.
* TODO: add instructions for new units.
* TODO: add examples of get and set.
