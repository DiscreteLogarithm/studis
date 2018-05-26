STUDIS Strongly Typed Units & Dimensions In SI
===============================================

Copyright 2018 Morteza Jalalvand
Licensed under the NDPL please see [Licence](#licence) for details.

Scientifically valid equations must be dimensionally homogeneous. It means that you can't compare quantities with different dimensions or add or subtract them. The argument of sine and many other mathematical functions must be a dimensionless quantity. Moreover, quantities of the same dimension but differing units should be converted to the same unit before comparing, adding or subtracting them. Breaking these rules in a program results in logical errors that can easily go undetected. STUDIS enforces the concept of dimensional homogeneity as syntax rules so that you get a compile error for violating it. It also internally converts all units to SI units so that quantities with differing units can be easily compared, added or subtracted.

Table of contents
-----------------

- [Usage](#usage)
- [Performance](#performance)
- [How many dimensions are there](#how-many-dimensions-are-there)
- [List of units](#list-of-units)
- [List of constants](#list-of-constants)
- [Acknowledgement](#acknowledgement)
- [Dedication](#dedication)
- [Licence](#licence)

Usage
-----

What you see in this section is basically the content of example.cpp.

You should begin by

```C++
#include "unit.hpp"

using namespace studis::literals;
using namespace studis::constants;
```

Then you can define and use quantities easily

```C++
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
```

All math functions that make sense for quantities with dimension are overloaded

```C++
std::cout << abs (-1_A) << std::endl;         // prints '1 A'
std::cout << atan2 (7_m, 1_km) << std::endl;  // prints some number
// std::cout << atan2 (1_m, 1_s) << std::endl;// error
std::cout << hypot (3_m, 4_m) << std::endl;   // prints '5 m'
// std::cout << hypot (1_m, 1_s) << std::endl;// error
```

`pow` is the only function that has different signature than its `std` counterpart, this can't be avoided since the dimension of the output depends on the power

```C++
auto energy = 0.5 * 1_kg * pow<2> (speed);
std::cout << energy << std::endl;             // prints '0.125 J (m2.kg/s2)'
```

`sqrt`, `cbrt` are overloaded for quantities whose result does not have a non-integer dimensional exponent

```C++
// pi, standard_gravity and many other constants are defined in the constants namespace
auto pendulum_frequency = sqrt (standard_gravity / 1_m) / (2*pi);
std::cout << pendulum_frequency << std::endl; // prints '0.498403 Hz (1/s)'
std::cout << cbrt (1_litre) << std::endl;     // prints '0.1 m'
```

Fractional power dimensions are not supported

```C++
// std::cout << sqrt (1_s) << std::endl;      // error
// std::cout << cbrt (1_m2) << std::endl;     // error
```

Dimensionless quantities can be used with any math function since they implicitly convert to a floating-point

```C++
auto pos = 1_cm * cos (2*pi*1_s*pendulum_frequency);
std::cout << pos << std::endl;                // prints '-0.0099995 m'
// std::cout << cos (1_s) << std::endl;       // error
```

There are so many units and prefixes in STUDIS

```C++
auto resistance = 1.7_kOhm;                   // we don't have greek letters, so that's kiloohm
auto inductance = 1_uH;                       // same, this is microhenry
auto capacitance = 1_pF;
if (resistance > 2*sqrt (inductance/capacitance))
  std::cout << "overdamped" << std::endl;
else if (resistance == 2*sqrt (inductance/capacitance))
  std::cout << "cricitally damped" << std::endl;
else std::cout << "underdamped" << std::endl;
```

You can use STUDIS simply as a unit convertor (to SI units)

```C++
std::cout << 10_ly << std::endl;              // prints '9.46073e+16 m'
std::cout << 30_knot << std::endl;            // prints '15.4333 m/s'
std::cout << 1_MeV << std::endl;              // prints '1.60218e-13 J (m2.kg/s2)'
std::cout << 2000_kcal << std::endl;          // prints '8.368e+09 J (m2.kg/s2)'
std::cout << 120_mmHg << std::endl;           // prints '15998.7 Pa (kg/m.s-2)'
```

And so many constants

```C++
auto radiative_power = Stefan_Boltzmann_constant * pow<4>(300_K) * 1_m2;
std::cout << radiative_power << std::endl;    // prints '459.3 W (m2.kg/s3)'
std::cout << electron_mass << std::endl;      // prints '9.10938e-31 kg'
```

Value of a (non-const) variable can change but its dimension can't

```C++
auto mass = 1_kg;
mass = 300_g;                                 // fine
// mass = 1_m3;                               // error
std::cin >> mass;                             // you can also read its value
std::cout << mass << std::endl;
```

If you don't want to specify an initial value (not recommended), you have to specify the dimension of the quantity

```C++
studis::Density d;
std::cin >> d;
std::cout << d << std::endl;
```

Many common dimensions are there, but in the case you can't find it there, you can specify the power for all 7 base dimensions of the SI yourself

```C++
studis::Quantity<studis::Dimension<1,0,-3,0,0,0,0>> jerk;
std::cin >> jerk;
std::cout << jerk << std::endl;
```

Performance
-----------

STUDIS should not incur any noticeable overhead at runtime. Information about dimension of quantities are encoded in the type system so they are not stored and only the value itself consumes memory. All dimension checks are of course performed during compilation and incur no cost at runtime.

How many dimensions are there?
------------------------------

Really a lot. Much more than any reasonable use case scenario. The dimensional exponents of quantities can always range from -127 to 127 (it could actually be more), so at least about 2<sup>56</sup>. In other words a quantity Q with dimension

dim Q = L<sup>α</sup> M<sup>β</sup> T<sup>γ</sup> I<sup>δ</sup> Θ<sup>ε</sup> N<sup>ζ</sup> J<sup>η</sup>

is guaranteed to be in STUDIS as long as all of α, β, γ, δ, ε, ζ, and η are integers in interval -127 to 127.

Common dimensions have type-aliases for easy access

| type-alias                  | dimension                                     |
|-----------------------------|-----------------------------------------------|
| `Dimmensionless`            | 1                                             |
| `Length`                    | L                                             |
| `Mass`                      | M                                             |
| `Time`, `Duration`          | T                                             |
| `ElectricCurrent`           | I                                             |
| `Temperature`               | Θ                                             |
| `AmountOfSubstance`         | N                                             |
| `LuminousIntensity`         | J                                             |
| `LuminousFlux`              | J                                             |
| `Wavenumber`                | L<sup>-1</sup>                                |
| `Area`                      | L<sup>2</sup>                                 |
| `Volume`                    | L<sup>3</sup>                                 |
| `CurrentDensity`            | L<sup>-2</sup> I                              |
| `Density`                   | L<sup>-3</sup> M                              |
| `Concentration`             | L<sup>-3</sup> N                              |
| `Velocity`, `Speed`         | L T<sup>-1</sup>                              |
| `Acceleration`              | L T<sup>-2</sup>                              |
| `Momentum`                  | L M T<sup>-1</sup>                            |
| `Action`                    | L<sup>2</sup> M T<sup>-1</sup>                |
| `Frequency`                 | T<sup>-1</sup>                                |
| `Radioactivity`             | T<sup>-1</sup>                                |
| `Force`                     | L M T<sup>-2</sup>                            |
| `Pressure`, `Stress`        | L<sup>-1</sup> M T<sup>-2</sup>               |
| `DynamicViscosity`          | L<sup>-1</sup> M T<sup>-1</sup>               |
| `KinematicViscosity`        | L<sup>2</sup> T<sup>-2</sup>                  |
| `Torque`                    | L<sup>2</sup> M T<sup>-2</sup>                |
| `Energy`, `Work`, `Heat`    | L<sup>2</sup> M T<sup>-2</sup>                |
| `Power`, `RadiantFlux`      | L<sup>2</sup> M T<sup>-3</sup>                |
| `HeatCapacity`              | L<sup>2</sup> M T<sup>-2</sup> Θ<sup>-1</sup> |
| `Entropy`                   | L<sup>2</sup> M T<sup>-2</sup> Θ<sup>-1</sup> |
| `ElectricCharge`            | T I                                           |
| `ElectricPotential`, `ElectromotiveForce`, `Voltage` | L<sup>2</sup> M T<sup>-3</sup> I<sup>-1</sup> |
| `Capacitance`               | L<sup>-2</sup> M<sup>-1</sup> T<sup>4</sup> I<sup>2</sup> |
| `Resistance`, `Impedance`   | L<sup>2</sup> M T<sup>-3</sup> I<sup>-2</sup> |
| `Conductance`, `Admittance` | L<sup>-2</sup> M<sup>-1</sup> T<sup>3</sup> I<sup>2</sup> |
| `MagneticFlux`              | L<sup>2</sup> M T<sup>-2</sup> I<sup>-1</sup> |
| `MagneticFluxDensity`       | M T<sup>-3</sup> I<sup>-1</sup>               |
| `Inductance`                | L<sup>2</sup> M T<sup>-2</sup> I<sup>-2</sup> |
| `Illuminance`               | L<sup>-2</sup> J                              |
| `CatalyticActivity`         | T<sup>-1</sup> N                              |

List of units
-------------

| Quantity              | Unit                      | Symbols                                       |
|-----------------------|---------------------------|-----------------------------------------------|
| `Length`              | metre                     | `fm`, `pm`, `nm`, `um`, `mm`, `cm`, `m`, `km`, `micron`|
| `Length`              | angstrom                  | `angstrom`                                    |
| `Length`              | inch                      | `in`                                          |
| `Length`              | foot                      | `ft`                                          |
| `Length`              | yard                      | `yd`                                          |
| `Length`              | mile                      | `mile`                                        |
| `Length`              | nautical mile             | `nautical_mile`                               |
| `Length`              | astronomical unit         | `au`                                          |
| `Length`              | light year                | `ly`, `kly`, `Mly`, `Gly`                     |
| `Length`              | parsec                    | `pc`, `kpc`, `Mpc`, `Gpc`                     |
| `Mass`                | gram                      | `fg`, `pg`, `ng`, `ug`, `mg`, `g`, `gr`, `kg` |
| `Mass`                | dalton                    | `Da`, `kDa`, `MDa`                            |
| `Mass`                | pound                     | `lb`                                          |
| `Mass`                | ounce                     | `oz`                                          |
| `Mass`                | tonne                     | `t`                                           |
| `Time`                | second                    | `fs`, `ps`, `ns`, `us`, `ms`, `s`, `sec`      |
| `Time`                | svedberg                  | `Svedberg`                                    |
| `Time`                | minute                    | `min`                                         |
| `Time`                | hour                      | `h`, `hour`                                   |
| `Time`                | day                       | `d`, `day`                                    |
| `Time`                | Julian year               | `julian_year`                                 |
| `ElectricCurrent`     | ampere                    | `nA`, `uA`, `mA`, `A`, `kA`                   |
| `Temperature`         | kelvin                    | `K`                                           |
| `Temperature`         | degree Celsius            | `deg_C`, `degree_Celsius`                     |
| `Temperature`         | degree Fahrenheit         | `deg_F`, `degree_Fahrenheit`                  |
| `AmountOfSubstance`   | mole                      | `nmol`, `umol`, `mmol`, `mol`, `kmol`         |
| `LuminousIntensity`   | candela                   | `cd`                                          |
| `Area`                |                           | `mm2`, `cm2`, `m2`, `km2`                     |
| `Area`                |                           | `in2`, `ft2`, `yd2`, `mile2`                  |
| `Area`                | barn                      | `barn`                                        |
| `Area`                | hectare                   | `ha`, `hectare`                               |
| `Volume`              |                           | `cm3`, `m3`                                   |
| `Volume`              | litre                     | `ul`, `uL`, `ml`, `mL`, `l`, `L`, `litre`     |
| `Density`             | gram per cubic centimetre | `gr_per_cm3`, `gr_per_ml`, `gr_per_mL`        |
| `Density`             | kilogram per litre        | `kg_per_l`, `kg_per_L`                        |
| `Density`             | kilogram per cubic metre  | `kg_per_m3`                                   |
| `Concentration`       | molar                     | `pM`, `nM`, `uM`, `mM`, `M`                   |
| `Velocity`            | metre per second          | `m_per_s`                                     |
| `Velocity`            | foot per second           | `ft_per_s`, `ft_per_sec`                      |
| `Velocity`            | kilometre per hour        | `km_per_hour`                                 |
| `Velocity`            | mile per hour             | `mile_per_hour`                               |
| `Velocity`            | knot                      | `knot`                                        |
| `Acceleration`        | metre per square second   | `m_per_s2`                                    |
| `Acceleration`        | foot per square second    | `ft_per_s2`                                   |
| `Acceleration`        | gal                       | `Gal`                                         |
| `Momentum`            | metre kilogram per second | `m_kg_per_s`                                  |
| `Action`              | joule second              | `J_s`                                         |
| `Frequency`           | hertz                     | `Hz`, `kHz`, `MHz`, `GHz`, `THz`              |
| `Frequency`           | Baud                      | `Bd`, `kBd`, `MBd`, `GBd`                     |
| `Frequency`           | FLOPS                     | `FLOPS`, `kFLOPS`, `MFLOPS`, `GFLOPS`, `TFLOPS` |
| `Frequency`           | revolutions per minute    | `rpm`                                         |
| `Frequency`           | frames per second         | `fps`                                         |
| `Radioactivity`       | becquerel                 | `Bq`                                          |
| `Force`               | newton                    | `pN`, `nN`, `uN`, `mN`, `N`, `kN`             |
| `Force`               | dyne                      | `dyn`, `dyne`                                 |
| `Force`               | pound force               | `lbf`                                         |
| `Pressure`            | pascal                    | `Pa`, `kPa`, `MPa`, `GPa`                     |
| `Pressure`            | torr                      | `mTorr`, `Torr`                               |
| `Pressure`            | millimetre of mercury     | `mmHg`, `cmHg`                                |
| `Pressure`            | psi                       | `psi`                                         |
| `Pressure`            | bar                       | `mbar`, `bar`                                 |
| `Pressure`            | standard atmosphere       | `atm`                                         |
| `DynamicViscosity`    | pascal second             | `Pa_s`                                        |
| `DynamicViscosity`    | poise                     | `cP`, `P`                                     |
| `KinematicViscosity`  | square metre per second   | `m2_per_s`                                    |
| `KinematicViscosity`  | stokes                    | `cSt`, `St`                                   |
| `Torque`              | newton metre              | `N_m`                                         |
| `Energy`              | joule                     | `J`, `kJ`, `MJ`, `GJ`                         |
| `Energy`              | electronvolt              | `eV`, `keV`, `MeV`, `GeV`                     |
| `Energy`              | erg                       | `erg`                                         |
| `Energy`              | watt hour                 | `Wh`, `kWh`                                   |
| `Energy`              | british thermal unit      | `BTU`                                         |
| `Energy`              | calorie                   | `cal`, `kcal`                                 |
| `Power`               | watt                      | `nW`, `uW`, `mW`, `W`, `kW`, `MW`, `GW`       |
| `ElectricCharge`      | coulomb                   | `pC`, `nC`, `uC`, `mC`, `C`                   |
| `ElectricCharge`      | ampere hour               | `mAh`, `Ah`                                   |
| `ElectricPotential`   | volt                      | `uV`, `mV`, `V`, `kV`, `MV`                   |
| `Capacitance`         | farad                     | `pF`, `nF`, `uF`, `mF`, `F`                   |
| `Resistance`          | ohm                       | `uOhm`, `mOhm`, `Ohm`, `kOhm`, `MOhm`, `GOhm` |
| `Conductance`         | siemens                   | `S`                                           |
| `MagneticFlux`        | weber                     | `nWb`, `uWb`, `mWb`, `Wb`                     |
| `MagneticFlux`        | maxwell                   | `Mx`                                          |
| `MagneticFluxDensity` | tesla                     | `uT`, `mT`, `T`                               |
| `MagneticFluxDensity` | gauss                     | `mG`, `G`                                     |
| `Inductance`          | henry                     | `uH`, `mH`, `H`                               |
| `Illuminance`         | lux                       | `lx`                                          |
| `CatalyticActivity`   | katal                     | `kat`                                         |

List of constants
-----------------

### Fundamental constants

These values are based on the 2014 self-consistent set of values of the constants and conversion factors of physics and chemistry recommended by the Committee on Data for Science and Technology (CODATA). These values are based on a least-squares adjustment that takes into account all data available up to 31 December 2014.

| Constants | Value | Unit |
|-----------|-------|------|
| `speed_of_light` | _c_ = 299792458 | m/s |
| `magnetic_constant`, `vacuum_permeability` | _μ_<sub>0</sub> = 4π * 10<sup>-7</sup> | N/A<sup>2</sup> |
| `electric_constant`, `vacuum_permittivity` | _ε_<sub>0</sub> = 1 / (_μ_<sub>0</sub>_c_<sup>2</sup>) | F/m |
| `Newtonian_constant_of_gravitation`, <br>`universal_gravitational_constant`, <br>`gravitational_constant` | _G_ = 6.67408 * 10<sup>-11</sup> | N/(m<sup>2</sup> kg<sup>2</sup>) |
| `Planck_constant` | ℎ = 6.626070040 * 10<sup>-34</sup> | J s |
| `reduced_Planck_constant` | ℏ = 1.054571800 * 10<sup>-34</sup> | J s |
| `elementary_charge` | _e_ = 1.6021766208 * 10<sup>-34</sup> | C |
| `magnetic_flux_quantum` | 𝛷<sub>0</sub> = ℎ / (2 _e_) = 2.067833831 * 10<sup>-34</sup> | Wb |
| `conductance_quantum` | _G_<sub>0</sub> = 2 _e_<sup>2</sup> / ℎ = 7.748091730 * 10<sup>-5</sup> | S |
| `electron_mass` | _m_<sub>e</sub> = 9.10938356 * 10<sup>-31</sup> | kg |
| `proton_mass` | _m_<sub>p</sub> = 1.672621898 * 10<sup>-27</sup> | kg |
| `proton_electron_mass_ratio` | _m_<sub>p</sub> / _m_<sub>e</sub> = 1836.15267389
| `fine_structure_constant` | _α_ = _e_<sup>2</sup> / (4 π _ε_<sub>0</sub> ℏ _c_) = 0.0072973525664
| `inverse_fine_structure_constant` | _α_<sup>-1</sup> = 137.035999139
| `Rydberg_constant` | _R_<sub>∞</sub> = _α_<sup>2</sup> _m_<sub>e</sub> _c_ / (2 ℎ) = 10973731.568508 | 1/m
| `Avogadro_constant` | _N_<sub>A</sub> = 6.022140857 * 10<sup>23</sup> | 1/mol |
| `Faraday_constant` | _F_ = _e_ _N_<sub>A</sub> = 96485.33289 | C/mol |
| `molar_gas_constant`, `universal_gas_constant`, <br>`gas_constant` | _R_ = 8.3144598 | J/(mol K) |
| `Boltzmann_constant` | _k_ = _R_ / _N_<sub>A</sub> = 1.38064852 * 10<sup>23</sup> | J/K
| `Stefan_Boltzmann_constant` | σ = (π<sup>2</sup> / 60) _k_<sup>4</sup> / (ℏ<sup>3</sup> c<sup>2</sup>) = 5.670367 | W/(m<sup>2</sup> K<sup>4</sup>) |
| `Wien_displacement_law_constant`, `Wien_constant` | b = 2.8977729 * 10<sup>-3</sup> | K m |

### Constants holding the value of non-SI units accepted for use with the International System of Units

| Constant    | Value                                     |
|-------------|-------------------------------------------|
| `minute`    | 1 min = 60 s                              |
| `hour`      | 1 h = 60 min = 3600 s                     |
| `day`       | 1 d = 24 h = 86400 s                      |
| `degree`    | 1° = (π/180) rad                          |
| `arcminute` | 1′ = (1/60)° = (π/10800) rad              |
| `arcsecond` | 1″ = (1/60)′ = (π/648000) rad             |
| `hectare`   | 1 ha = 10<sup>4</sup> m<sup>2</sup>       |
| `litre`     | 1 L = 1 l = 10<sup>-3</sup> m<sup>3</sup> |
| `tonne`     | 1 t = 10<sup>3</sup> kg                   |

### Constants holding the value of non-SI units associated with the CGS and the CGS-Gaussian system of units

| Constant  | Value                                     |
|-----------|-------------------------------------------|
| `erg`     | 1 erg = 10<sup>-7</sup> J
| `dyne`    | 1 dyn = 10<sup>-5</sup> N
| `poise`   | 1 P = 1 dyn s cm<sup>-2</sup> = 0.1 Pa s
| `stokes`  | 1 St = 1 cm<sup>2</sup>/s = 10<sup>-4</sup> m<sup>2</sup>/s
| `gauss`   | 1 G = 1 Mx/cm<sup>2</sup> = 10<sup>-4</sup> T
| `maxwell` | 1 Mx = 1 G cm<sup>2</sup> = 10<sup>-8</sup> Wb

### Constants holding the value of non-SI units defined by the International Astronomical Union (IAU)

| Constant            | Value                                     |
|---------------------|-------------------------------------------|
| `julian_year`       | 365.25 day
| `astronomical_unit` | 149597870700 m
| `light_year`        | Product of Julian year and speed of light
| `parsec`            | (648000/π) astronomical units

### Adopted values

| Constant              | Value                     | Unit             | Remarks |
|-----------------------|---------------------------|------------------|---------|
| `standard_gravity`    | _g_<sub>n</sub> = 9.80665 | m/s<sup>2</sup>  |
| `standard_atmosphere` | atm = 101325              | Pa               |
| `mercury_density`     | ρ<sub>Hg</sub> = 13595.1  | kg/m<sup>3</sup> | Density used in the definition of mmHg

### Constants holding the value of UK and US custmary units

| Constant | Value |
|----------|-------|
| `inch`   | 1 in = 2.54 cm
| `foot`   | 1 ft = 12 in
| `yard`   | 1 yd = 3 ft
| `mile`   | 1 mile = 1760 yd
| `nautical_mile` | 1 nautical mile = 1852 m
| `knot`   | 1 knot = 1 nautical mile per hour
| `pound`  | 1 lb = 0.45359237 kg
| `ounce`  | 1 oz = (1/16) lb
| `pound_force` | 1 lbf = 1 lb * _g_<sub>n</sub>
| `pound_force_per_squared_inch` | 1 psi = 1 lbf/in<sup>2</sup>
| `british_thermal_unit` | 1 BTU = 788169 ft lbf
| `thermochemical_calorie` | 1 cal = 4184 J

### Constants holding the value of other non-SI units

| Constant                | Value |
|-------------------------|-------|
| `angstrom`              | 1 Å = 10<sup>-10</sup> m
| `svedberg`              | 1 S = 10<sup>-13</sup> s
| `torr`                  | 1 Torr = (1/760) atm
| `millimeter_of_mercury` | 1 mmHg = ρ<sub>Hg</sub> * _g_<sub>n</sub> * 1 mm
| `watt_hour`             | 1 Wh = 1 W * 1 h
| `ampere_hour`           | 1 Ah = 1 A * 1 h

Acknowledgement
---------------

This is inspired by the idea of a strongly typed template MKS unit system discussed in the book _The C++ Programming Language_ by _Bjarne Stroustrup_.

Dedication
----------

This library is dedicated to all my mentors particularly Seyed Mehdi Vaez Allaei and Mohammad A. Charsooghi to whom I am grateful for both their teachings and friendship.

Licence
-------

This library is distributed under the terms of Non-Discriminatory Public Licence. You can read the exact licence terms in the 'LICENSE' file, but here is a summary:

- You can use and modify the software
- You can distribute the original or the modified version of the software under the same terms in a non-discriminatory manner if you also provide the source code

If you have to comply with laws that compels you to restrict access of certain groups of people (such as export control laws), you can only use and modify this software for your own purposes, but you can no longer distribute it.
