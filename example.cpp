#include <iostream>
#include "studis.hpp"

using namespace studis::literals;
using namespace studis::constants;

int main () {
  auto l1 = 1.5_m, l2 = 2_cm;
  auto t = 3_s;
  auto l3 = l1 + l2;                            // fine
  std::cout << l3 << std::endl;                 // prints '1.52 m' (yes the unit is printed as well)
  std::cout << l1 + l2 << std::endl;            // same
  std::cout << (l1 < l2) << std::endl;          // works
  // std::cout << l1 + t << std::endl;          // error
  // std::cout << (l1 < t) << std::endl;        // error
  auto speed = l1 / t;
  std::cout << speed << std::endl;              // prints '0.5 m/s'
  
  // All math functions that make sense for quantities with dimension are overloaded
  std::cout << abs (-1_A) << std::endl;         // prints '1 A'
  std::cout << atan2 (7_m, 1_km) << std::endl;  // prints some number
  // std::cout << atan2 (1_m, 1_s) << std::endl;// error
  std::cout << hypot (3_m, 4_m) << std::endl;   // prints '5 m'
  // std::cout << hypot (1_m, 1_s) << std::endl;// error
  
  // pow is the only function that has different signature than its std counterpart
  // This can't be avoided since the dimension of output depends on the power
  auto energy = 0.5 * 1_kg * pow<2> (speed);
  std::cout << energy << std::endl;             // prints '0.125 J (m2.kg/s2)'
  
  // sqrt, cbrt are overloaded for quantities whose result is not a fractional dimension
  // pi, standard_gravity and many other constants are defined in the constants namespace
  auto pendulum_frequency = sqrt (standard_gravity / 1_m) / (2*pi);
  std::cout << pendulum_frequency << std::endl; // prints '0.498403 Hz (1/s)'
  std::cout << cbrt (1_litre) << std::endl;     // prints '0.1 m'
  // fractional power dimensions are not supported
  // std::cout << sqrt (1_s) << std::endl;      // error
  // std::cout << cbrt (1_m2) << std::endl;     // error
  
  // Dimensionless quantities can be used with any math function
  // Since they implicitly convert to a floating-point
  auto pos = 1_cm * cos (2*pi*1_s*pendulum_frequency);
  std::cout << pos << std::endl;                // prints '-0.0099995 m'
  // std::cout << cos (1_s) << std::endl;       // error
  
  // There are so many units and prefixes here
  auto resistance = 1.7_kOhm;                   // we don't have greek letters, so that's kiloohm
  auto inductance = 1_uH;                       // same, this is microhenry
  auto capacitance = 1_pF;
  if (resistance > 2*sqrt (inductance/capacitance))
    std::cout << "overdamped" << std::endl;
  else if (resistance == 2*sqrt (inductance/capacitance))
    std::cout << "cricitally damped" << std::endl;
  else std::cout << "underdamped" << std::endl;
  // You can use STUDIS simply as a unit convertor (to SI units)
  std::cout << 10_ly << std::endl;              // prints '9.46073e+16 m'
  std::cout << 30_knot << std::endl;            // prints '15.4333 m/s'
  std::cout << 1_MeV << std::endl;              // prints '1.60218e-13 J (m2.kg/s2)'
  std::cout << 2000_kcal << std::endl;          // prints '8.368e+09 J (m2.kg/s2)'
  std::cout << 120_mmHg << std::endl;           // prints '15998.7 Pa (kg/m.s-2)'
  
  // And so many constants
  auto radiative_power = Stefan_Boltzmann_constant * pow<4>(300_K) * 1_m2;
  std::cout << radiative_power << std::endl;    // prints '459.3 W (m2.kg/s3)'
  std::cout << electron_mass << std::endl;      // prints '9.10938e-31 kg'
  
  // Value of a (non-const) variable can change but its dimension can't
  auto mass = 1_kg;
  mass = 300_g;                                 // fine
  // mass = 1_m3;                               // error
  std::cin >> mass;                             // you can also read its value
  std::cout << mass << std::endl;
  
  // If you don't want to specify an initial value (not recommended)
  // You have to specify the dimension of the quantity
  studis::Density d;
  std::cin >> d;
  std::cout << d << std::endl;
  
  // Many common dimensions are there, but in the case you can't find it there
  // You can specify the power for all 7 base dimensions of the SI yourself.
  studis::Quantity<studis::Dimension<1,0,-3,0,0,0,0>> jerk;
  std::cin >> jerk;
  std::cout << jerk << std::endl;
  
  return 0;
}
