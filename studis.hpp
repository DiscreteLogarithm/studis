#ifndef UNIT_HPP
#define UNIT_HPP

#include <cmath>
#include <cstdint>
#include <sstream>
#include <type_traits>

#if __cplusplus < 201402L
#error C++14 compatible compiler is required
#endif

namespace studis {
  
  template <int_fast8_t _L, int_fast8_t _M, int_fast8_t _T, int_fast8_t _I, int_fast8_t _Temp, int_fast8_t _N, int_fast8_t _J>
  struct Dimension {
    enum {L=_L, M=_M, T=_T, I=_I, Temp=_Temp, N=_N, J=_J};
  };
  
  template <typename T>
  constexpr bool is_dimension_v = false;
  
  template <int_fast8_t L, int_fast8_t M, int_fast8_t T, int_fast8_t I, int_fast8_t Temp, int_fast8_t N, int_fast8_t J>
  constexpr bool is_dimension_v<Dimension<L,M,T,I,Temp,N,J>> = true;
  
  template <typename D>
  std::string str () {
    static_assert (is_dimension_v<D>, "invalid template parameter");
    constexpr auto L = D::L;
    constexpr auto M = D::M;
    constexpr auto T = D::T;
    constexpr auto I = D::I;
    constexpr auto Temp = D::Temp;
    constexpr auto N = D::N;
    constexpr auto J = D::J;
    int_fast8_t num_pos_dims = (L>0) + (M>0) + (T>0) + (I>0) + (Temp>0) + (N>0) + (J>0);
    int_fast8_t num_neg_dims = (L<0) + (M<0) + (T<0) + (I<0) + (Temp<0) + (N<0) + (J<0);
    std::ostringstream oss;
    if (num_pos_dims==0 && num_neg_dims==0) return "";
    else if (num_pos_dims) {
      if (L>0) {
        oss << "m";
        if (L!=1) oss << L;
        if (--num_pos_dims) oss << '.';
      }
      if (M>0) {
        oss << "kg";
        if (M!=1) oss << M;
        if (--num_pos_dims) oss << '.';
      }
      if (T>0) {
        oss << "s";
        if (T!=1) oss << T;
        if (--num_pos_dims) oss << '.';
      }
      if (I>0) {
        oss << "A";
        if (I!=1) oss << I;
        if (--num_pos_dims) oss << '.';
      }
      if (Temp>0) {
        oss << "K";
        if (Temp!=1) oss << Temp;
        if (--num_pos_dims) oss << '.';
      }
      if (N>0) {
        oss << "mol";
        if (N!=1) oss << N;
        if (--num_pos_dims) oss << '.';
      }
      if (J>0) {
        oss << "cd";
        if (J!=1) oss << J;
      }
    } else oss << '1';
    if (num_neg_dims) {
      oss << '/';
      if (L<0) {
        oss << "m";
        if (L!=-1) oss << -L;
        if (--num_neg_dims) oss << '.';
      }
      if (M<0) {
        oss << "kg";
        if (M!=-1) oss << -M;
        if (--num_neg_dims) oss << '.';
      }
      if (T<0) {
        oss << "s";
        if (T!=-1) oss << -T;
        if (--num_neg_dims) oss << '.';
      }
      if (I<0) {
        oss << "A";
        if (I!=-1) oss << -I;
        if (--num_neg_dims) oss << '.';
      }
      if (Temp<0) {
        oss << "K";
        if (Temp!=-1) oss << -Temp;
        if (--num_neg_dims) oss << '.';
      }
      if (N<0) {
        oss << "mol";
        if (N!=-1) oss << -N;
        if (--num_neg_dims) oss << '.';
      }
      if (J<0) {
        oss << "cd";
        if (J!=-1) oss << -J;
      }
    }
    return oss.str();
  }
  
  template <typename U1, typename U2>
  using dimensions_multiply_t = Dimension<U1::L+U2::L, U1::M+U2::M, U1::T+U2::T, U1::I+U2::I, U1::Temp+U2::Temp, U1::N+U2::N, U1::J+U2::J>;
  
  template <typename U1, typename U2>
  using dimensions_divide_t = Dimension<U1::L-U2::L, U1::M-U2::M, U1::T-U2::T, U1::I-U2::I, U1::Temp-U2::Temp, U1::N-U2::N, U1::J-U2::J>;
  
  template <typename D, int_fast8_t p>
  using dimension_power_t = Dimension<D::L*p, D::M*p, D::T*p, D::I*p, D::Temp*p, D::N*p, D::J*p>;
  
  template <typename T>
  constexpr bool is_dimmensionless_v = false;
  
  template <int_fast8_t L, int_fast8_t M, int_fast8_t T, int_fast8_t I, int_fast8_t Temp, int_fast8_t N, int_fast8_t J>
  constexpr bool is_dimmensionless_v<Dimension<L,M,T,I,Temp,N,J>> = L==0 && M==0 && T==0 && I==0 && Temp==0 && N==0 && J==0;
  
  template <typename D>
  constexpr bool is_squared_v = std::abs(D::L)%2==0 && std::abs(D::M)%2==0 && std::abs(D::T)%2==0 && std::abs(D::I)%2==0 && std::abs(D::Temp)%2==0 && std::abs(D::N)%2==0 && std::abs(D::J)%2==0;
  
  template <typename D>
  constexpr bool is_cubed_v = std::abs(D::L)%3==0 && std::abs(D::M)%3==0 && std::abs(D::T)%3==0 && std::abs(D::I)%3==0 && std::abs(D::Temp)%3==0 && std::abs(D::N)%3==0 && std::abs(D::J)%3==0;
  
  template <typename D>
  using dimension_sqrt_t = Dimension<D::L/2, D::M/2, D::T/2, D::I/2, D::Temp/2, D::N/2, D::J/2>;
  
  template <typename D>
  using dimension_cbrt_t = Dimension<D::L/3, D::M/3, D::T/3, D::I/3, D::Temp/3, D::N/3, D::J/3>;
  
  using No_Dim                  = Dimension<0,0,0,0,0,0,0>;
  using Length_Dim              = Dimension<1,0,0,0,0,0,0>;
  using Mass_Dim                = Dimension<0,1,0,0,0,0,0>;
  using Time_Dim                = Dimension<0,0,1,0,0,0,0>;
  using ElectricCurrent_Dim     = Dimension<0,0,0,1,0,0,0>;
  using Temperature_Dim         = Dimension<0,0,0,0,1,0,0>;
  using AmountOfSubstance_Dim     = Dimension<0,0,0,0,0,1,0>;
  using LuminousIntensity_Dim   = Dimension<0,0,0,0,0,0,1>;
  using LuminousFlux_Dim        = Dimension<0,0,0,0,0,0,1>;
  
  using Wavenumber_Dim          = dimension_power_t<Length_Dim,-1>;
  using Area_Dim                = dimension_power_t<Length_Dim,2>;
  using Volume_Dim              = dimension_power_t<Length_Dim,3>;
  using CurrentDensity_Dim      = dimensions_divide_t<ElectricCurrent_Dim,Area_Dim>;
  using Density_Dim             = dimensions_divide_t<Mass_Dim,Volume_Dim>;
  using Concentration_Dim       = dimensions_divide_t<AmountOfSubstance_Dim,Volume_Dim>;
  using Velocity_Dim            = dimensions_divide_t<Length_Dim,Time_Dim>;
  using Accelaration_Dim        = dimensions_divide_t<Velocity_Dim,Time_Dim>;
  using Momentum_Dim            = dimensions_multiply_t<Mass_Dim,Velocity_Dim>;
  using Action_Dim              = dimensions_multiply_t<Length_Dim,Momentum_Dim>;
  
  using Frequency_Dim           = dimension_power_t<Time_Dim,-1>;
  using Radioactivity_Dim       = dimension_power_t<Time_Dim,-1>;
  using Force_Dim               = dimensions_multiply_t<Mass_Dim,Accelaration_Dim>;
  using Pressure_Dim            = dimensions_divide_t<Force_Dim,Area_Dim>;
  using DynamicViscosity_Dim    = dimensions_divide_t<Pressure_Dim,dimensions_divide_t<Velocity_Dim,Length_Dim>>;
  using KinematicViscosity_Dim  = dimensions_divide_t<DynamicViscosity_Dim,Density_Dim>;
  using Torque_Dim              = dimensions_multiply_t<Force_Dim,Length_Dim>;
  using Energy_Dim              = dimensions_multiply_t<Force_Dim,Length_Dim>;
  using Power_Dim               = dimensions_divide_t<Energy_Dim,Time_Dim>;
  using ElectricCharge_Dim      = dimensions_multiply_t<ElectricCurrent_Dim,Time_Dim>;
  using Voltage_Dim             = dimensions_divide_t<Power_Dim,ElectricCurrent_Dim>;
  using Capacitance_Dim         = dimensions_divide_t<ElectricCharge_Dim,Voltage_Dim>;
  using Impedance_Dim           = dimensions_divide_t<Voltage_Dim,ElectricCurrent_Dim>;
  using Conductance_Dim         = dimension_power_t<Impedance_Dim,-1>;
  using MagneticFlux_Dim        = dimensions_multiply_t<Voltage_Dim,Time_Dim>;
  using MagneticFluxDensity_Dim = dimensions_divide_t<MagneticFlux_Dim,Area_Dim>;
  using Inductance_Dim          = dimensions_divide_t<MagneticFlux_Dim,ElectricCurrent_Dim>;
  using Illuminance_Dim         = dimensions_divide_t<LuminousIntensity_Dim,Area_Dim>;
  using Dose_Dim                = dimensions_divide_t<Energy_Dim,Mass_Dim>;
  using CatalyticActivity_Dim   = dimensions_divide_t<AmountOfSubstance_Dim,Time_Dim>;
  
  using HeatCapacity_Dim        = dimensions_divide_t<Energy_Dim,Temperature_Dim>;
  using Entropy_Dim             = dimensions_divide_t<Energy_Dim,Temperature_Dim>;
  
  template <>
  std::string str<Frequency_Dim> () {return "Hz (1/s)";}
  template <>
  std::string str<Force_Dim> () {return "N (m.kg/s2)";}
  template <>
  std::string str<Pressure_Dim> () {return "Pa (kg/m.s2)";}
  template <>
  std::string str<Energy_Dim> () {return "J (m2.kg/s2)";}
  template <>
  std::string str<Power_Dim> () {return "W (m2.kg/s3)";}
  template <>
  std::string str<ElectricCharge_Dim> () {return "C (s.A)";}
  template <>
  std::string str<Voltage_Dim> () {return "V (m2.kg/s3.A)";}
  template <>
  std::string str<Capacitance_Dim> () {return "F (s4.A2/m2.kg)";}
  template <>
  std::string str<Impedance_Dim> () {return "Ohm (m2.kg/s3.A2)";}
  template <>
  std::string str<Conductance_Dim> () {return "S (s3.A2/m2.kg)";}
  template <>
  std::string str<MagneticFlux_Dim> () {return "Wb (m2.kg/s2.A)";}
  template <>
  std::string str<MagneticFluxDensity_Dim> () {return "T (kg/s2.A)";}
  template <>
  std::string str<Inductance_Dim> () {return "H (m2.kg/s2.A2)";}
//   template <>
//   std::string str<Illuminance_Dim> () {return "lx";}
//   template <>
//   std::string str<Dose_Dim> () {return "Gy|Sv";}
  template <>
  std::string str<CatalyticActivity_Dim> () {return "kat (mol/s)";}
  
  template <typename D, typename float_t = double>
  class Quantity;
  
  template <typename D, typename float_t>
  using quantity_sqrt_t = std::enable_if_t<is_squared_v<D>,Quantity<dimension_sqrt_t<D>,float_t>>;
  
  template <typename D, typename float_t>
  using quantity_cbrt_t = std::enable_if_t<is_cubed_v<D>,Quantity<dimension_cbrt_t<D>,float_t>>;
  
  using Dimmensionless      = Quantity<No_Dim>;
  using Length              = Quantity<Length_Dim>;
  using Mass                = Quantity<Mass_Dim>;
  using Time                = Quantity<Time_Dim>;
  using Duration            = Time;
  using ElectricCurrent     = Quantity<ElectricCurrent_Dim>;
  using Temperature         = Quantity<Temperature_Dim>;
  using AmountOfSubstance   = Quantity<AmountOfSubstance_Dim>;
  using LuminousIntensity   = Quantity<LuminousIntensity_Dim>;
  using LuminousFlux        = Quantity<LuminousFlux_Dim>;
  
  using Wavenumber          = Quantity<Wavenumber_Dim>;
  using Area                = Quantity<Area_Dim>;
  using Volume              = Quantity<Volume_Dim>;
  using CurrentDensity      = Quantity<CurrentDensity_Dim>;
  using Density             = Quantity<Density_Dim>;
  using Concentration       = Quantity<Concentration_Dim>;
  using Velocity            = Quantity<Velocity_Dim>;
  using Accelaration        = Quantity<Accelaration_Dim>;
  using Momentum            = Quantity<Momentum_Dim>;
  using Action              = Quantity<Action_Dim>;
  
  using Frequency           = Quantity<Frequency_Dim>;
  using Radioactivity       = Quantity<Radioactivity_Dim>;
  using Force               = Quantity<Force_Dim>;
  using Pressure            = Quantity<Pressure_Dim>;
  using DynamicViscosity    = Quantity<DynamicViscosity_Dim>;
  using KinematicViscosity  = Quantity<KinematicViscosity_Dim>;
  using Torque              = Quantity<Torque_Dim>;
  using Energy              = Quantity<Energy_Dim>;
  using Work                = Energy;
  using Heat                = Energy;
  using Power               = Quantity<Power_Dim>;
  using RadiantFlux         = Power;
  using ElectricCharge      = Quantity<ElectricCharge_Dim>;
  using Voltage             = Quantity<Voltage_Dim>;
  using ElectricPotential   = Voltage;
  using ElectromotiveForce  = Voltage;
  using Capacitance         = Quantity<Capacitance_Dim>;
  using Impedance           = Quantity<Impedance_Dim>;
  using Conductance         = Quantity<Conductance_Dim>;
  using MagneticFlux        = Quantity<MagneticFlux_Dim>;
  using MagneticFluxDensity = Quantity<MagneticFluxDensity_Dim>;
  using Inductance          = Quantity<Inductance_Dim>;
  using Illuminance         = Quantity<Illuminance_Dim>;
  using Dose                = Quantity<Dose_Dim>;
  using CatalyticActivity   = Quantity<CatalyticActivity_Dim>;
  
  using HeatCapacity        = Quantity<HeatCapacity_Dim>;
  using Entropy             = Quantity<Entropy_Dim>;
  
  template <typename D, typename float_t>
  class Quantity final {
    static_assert (is_dimension_v<D>, "invalid template parameter");
    static_assert (std::is_floating_point<float_t>::value, "invalid template parameter");
    float_t val;
  public:
    explicit constexpr Quantity () noexcept = default;
    explicit constexpr Quantity (float_t d) noexcept : val{d} {}
    
    template <typename arithmetic_t, typename = std::enable_if_t<std::is_arithmetic<arithmetic_t>::value>>
    explicit constexpr Quantity (arithmetic_t d) noexcept : val{static_cast<float_t>(d)} {}
    
    constexpr Quantity (const Quantity &) noexcept = default;
    constexpr Quantity (Quantity &&) noexcept = default;
    
    constexpr Quantity & operator= (const Quantity &) & noexcept = default;
    constexpr Quantity & operator= (const Quantity &) && noexcept = delete;
    constexpr Quantity & operator= (Quantity &&) & noexcept = default;
    constexpr Quantity & operator= (Quantity &&) && noexcept = delete;
    
    using Dimension = D;
    
    constexpr float_t value () const {return val;}
    
    template <typename rhs_float_t>
    constexpr Quantity operator+= (const Quantity<D,rhs_float_t> &rhs) & noexcept {
      val += rhs.val;
      return *this;
    }
    template <typename rhs_float_t>
    constexpr Quantity operator-= (const Quantity<D,rhs_float_t> &rhs) & noexcept {
      val -= rhs.val;
      return *this;
    }
    
    constexpr Quantity operator*= (float_t rhs) & noexcept {
      val *= rhs;
      return *this;
    }
    constexpr Quantity operator/= (float_t rhs) & noexcept {
      val /= rhs;
      return *this;
    }
    
    constexpr Quantity operator+ () const noexcept {
      return *this;
    }
    constexpr Quantity operator- () const noexcept {
      return Quantity{-val};
    }
    
    constexpr auto operator+ (Quantity rhs) const noexcept {
      return Quantity{val+rhs.val};
    }
    constexpr auto operator- (Quantity rhs) const noexcept {
      return Quantity{val-rhs.val};
    }
    
    constexpr auto operator* (float_t rhs) const noexcept {
      return Quantity{val*rhs};
    }
    constexpr auto operator/ (float_t rhs) const noexcept {
      return Quantity{val/rhs};
    }
    
    friend constexpr auto operator* (float_t lhs, Quantity rhs) noexcept {
      return Quantity{lhs*rhs.val};
    }
    friend constexpr auto operator/ (float_t lhs, Quantity rhs) noexcept {
      return Quantity<dimension_power_t<D,-1>>{lhs/rhs.val};
    }
    
    constexpr bool operator== (Quantity rhs) const noexcept {
      return val == rhs.val;
    }
    constexpr bool operator!= (Quantity rhs) const noexcept {
      return val != rhs.val;
    }
    
    constexpr bool operator> (Quantity rhs) const noexcept {
      return val > rhs.val;
    }
    constexpr bool operator< (Quantity rhs) const noexcept {
      return val < rhs.val;
    }
    
    constexpr bool operator>= (Quantity rhs) const noexcept {
      return val >= rhs.val;
    }
    constexpr bool operator<= (Quantity rhs) const noexcept {
      return val <= rhs.val;
    }
    
    friend constexpr std::istream & operator>> (std::istream &lhs, Quantity &rhs) noexcept {
      return lhs >> rhs.val;
    }
  };
  
  template <typename dim, typename lhs_float_t, typename rhs_float_t>
  constexpr auto operator+ (Quantity<dim, lhs_float_t> lhs, Quantity<dim, rhs_float_t> rhs) noexcept {
    return Quantity<dim, std::common_type_t<lhs_float_t,rhs_float_t>> {lhs.value()+rhs.value()};
  }
  template <typename dim, typename lhs_float_t, typename rhs_float_t>
  constexpr auto operator- (Quantity<dim, lhs_float_t> lhs, Quantity<dim, rhs_float_t> rhs) noexcept {
    return Quantity<dim, std::common_type_t<lhs_float_t,rhs_float_t>> {lhs.value()-rhs.value()};
  }
  
  template <typename lhs_dim, typename rhs_dim, typename float_t>
  constexpr auto operator* (Quantity<lhs_dim, float_t> lhs, Quantity<rhs_dim, float_t> rhs) noexcept {
    return Quantity<dimensions_multiply_t<lhs_dim,rhs_dim>, float_t> {lhs.value()*rhs.value()};
  }
  template <typename lhs_dim, typename rhs_dim, typename float_t>
  constexpr auto operator/ (Quantity<lhs_dim, float_t> lhs, Quantity<rhs_dim, float_t> rhs) noexcept {
    return Quantity<dimensions_divide_t<lhs_dim,rhs_dim>, float_t> {lhs.value()/rhs.value()};
  }
  
  template <typename lhs_dim, typename lhs_float_t, typename rhs_dim, typename rhs_float_t>
  constexpr auto operator* (Quantity<lhs_dim, lhs_float_t> lhs, Quantity<rhs_dim, rhs_float_t> rhs) noexcept {
    return Quantity<dimensions_multiply_t<lhs_dim,rhs_dim>, std::common_type_t<lhs_float_t,rhs_float_t>> {lhs.value()*rhs.value()};
  }
  template <typename lhs_dim, typename lhs_float_t, typename rhs_dim, typename rhs_float_t>
  constexpr auto operator/ (Quantity<lhs_dim, lhs_float_t> lhs, Quantity<rhs_dim, rhs_float_t> rhs) noexcept {
    return Quantity<dimensions_divide_t<lhs_dim,rhs_dim>, std::common_type_t<lhs_float_t,rhs_float_t>> {lhs.value()/rhs.value()};
  }
  
  template <typename D, typename float_t>
  constexpr auto abs (Quantity<D, float_t> x) noexcept {
    return Quantity<D, float_t>{std::abs (x.value())};
  }
  
  template <typename D, typename y_float_t, typename x_float_t>
  constexpr auto atan2 (Quantity<D, y_float_t> y, Quantity<D, x_float_t> x) {
    return std::atan2 (y.value(), x.value());
  }
  
  template <typename D, typename float_t>
  constexpr auto hypot (Quantity<D, float_t> x, Quantity<D, float_t> y) {
    return Quantity<D, float_t>{std::hypot (x.value(), y.value())};
  }
  
  template <typename D, typename x_float_t, typename y_float_t>
  constexpr auto hypot (Quantity<D, x_float_t> x, Quantity<D, y_float_t> y) {
    return Quantity<D, std::common_type_t<x_float_t,y_float_t>>{std::hypot (x.value(), y.value())};
  }
  
  template <typename float_t>
  constexpr float_t pos_int_power (float_t x, uint_fast8_t p) {
    return p ? x*pos_int_power (x, p-1) : 1;
  }
  
  template <typename float_t>
  constexpr float_t int_power (float_t x, int_fast8_t p) {
    return pos_int_power (p >= 0 ? x : 1/x, std::abs(p));
  }
  
  template <int_fast8_t P, typename D, typename float_t>
  constexpr auto pow (Quantity<D, float_t> x) noexcept {
    return Quantity<dimension_power_t<D,P>, float_t>{int_power (x.value(), P)};
  }
  
  template <typename D, typename float_t, typename = std::enable_if_t<is_squared_v<D>>>
  constexpr auto sqrt (Quantity<D, float_t> x) noexcept {
    return quantity_sqrt_t<D, float_t>{std::sqrt(x.value())};
  }
  
  template <typename D, typename float_t, typename = std::enable_if_t<is_cubed_v<D>>>
  constexpr auto cbrt (Quantity<D, float_t> x) noexcept {
    return quantity_cbrt_t<D, float_t>{std::cbrt(x.value())};
  }
  
  template <typename D, typename float_t>
  constexpr std::ostream & operator<< (std::ostream &lhs, Quantity<D, float_t> rhs) {
    return lhs << rhs.value() << ' ' << str<D>();
  }
  
  template <typename float_t>
  class Quantity<No_Dim, float_t> final {
    float_t val;
  public:
    explicit constexpr Quantity () noexcept = default;
    explicit constexpr Quantity (float_t d) : val{d} {}
    
    constexpr Quantity (const Quantity &) noexcept = default;
    constexpr Quantity (Quantity &&) noexcept = default;
    
    constexpr Quantity & operator= (const Quantity &) & noexcept = default;
    constexpr Quantity & operator= (const Quantity &) && noexcept = delete;
    constexpr Quantity & operator= (Quantity &&) & noexcept = default;
    constexpr Quantity & operator= (Quantity &&) && noexcept = delete;
    
    constexpr float_t value () const {return val;}
    
    constexpr Quantity operator+= (float_t rhs) & noexcept {
      val += rhs;
    }
    constexpr Quantity operator-= (float_t rhs) & noexcept {
      val -= rhs;
    }
    
    constexpr Quantity operator*= (float_t rhs) & noexcept {
      val *= rhs;
    }
    constexpr Quantity operator/= (float_t rhs) & noexcept {
      val /= rhs;
    }
    
    friend constexpr std::istream & operator>> (std::istream &lhs, Quantity &rhs) noexcept {
      return lhs >> rhs.val;
    }
    
    constexpr operator float_t () const noexcept {return val;}
  };
  
  namespace constants {
    constexpr double pi = 3.141592653589793238L;
  }
  
  namespace literals {
    
    constexpr Length operator"" _m (long double d) {return Length{d};}
    constexpr Length operator"" _m (unsigned long long d) {return Length{d};}
    constexpr Length operator"" _cm (long double d) {return d*1e-2_m;}
    constexpr Length operator"" _cm (unsigned long long d) {return d*1.0_cm;}
    
    constexpr Mass operator"" _kg (long double d) {return Mass{d};}
    constexpr Mass operator"" _kg (unsigned long long d) {return Mass{d};}
    constexpr Mass operator"" _g  (long double d) {return d*1e-3_kg;}
    constexpr Mass operator"" _g  (unsigned long long d) {return d*1.0_g;}
    
    constexpr Time operator"" _s (long double d) {return Time{d};}
    constexpr Time operator"" _s (unsigned long long d) {return Time{d};}
    
    constexpr ElectricCurrent operator"" _A (long double d) {return ElectricCurrent{d};}
    constexpr ElectricCurrent operator"" _A (unsigned long long d) {return ElectricCurrent{d};}
    
    constexpr Temperature operator"" _K (long double d) {return Temperature{d};}
    constexpr Temperature operator"" _K (unsigned long long d) {return Temperature{d};}
    
    constexpr AmountOfSubstance operator"" _mol (long double d) {return AmountOfSubstance{d};}
    constexpr AmountOfSubstance operator"" _mol (unsigned long long d) {return AmountOfSubstance{d};}
    
    constexpr LuminousIntensity operator"" _cd (long double d) {return LuminousIntensity{d};}
    constexpr LuminousIntensity operator"" _cd (unsigned long long d) {return LuminousIntensity{d};}
    
    constexpr Area operator"" _m2 (long double d) {return Area{d};}
    constexpr Area operator"" _m2 (unsigned long long d) {return Area{d};}
    
    constexpr Volume operator"" _m3 (long double d) {return Volume{d};}
    constexpr Volume operator"" _m3 (unsigned long long d) {return Volume{d};}
    
    constexpr Density operator"" _kg_per_m3 (long double d) {return Density{d};}
    constexpr Density operator"" _kg_per_m3 (unsigned long long d) {return Density{d};}
    
    constexpr Concentration operator"" _mol_per_m3 (long double d) {return Concentration{d};}
    constexpr Concentration operator"" _mol_per_m3 (unsigned long long d) {return Concentration{d};}
    
    constexpr Velocity operator"" _m_per_s (long double d) {return Velocity{d};}
    constexpr Velocity operator"" _m_per_s (unsigned long long d) {return Velocity{d};}
    
    constexpr Accelaration operator"" _m_per_s2 (long double d) {return Accelaration{d};}
    constexpr Accelaration operator"" _m_per_s2 (unsigned long long d) {return Accelaration{d};}
    
    constexpr Momentum operator"" _kg_m_per_s (long double d) {return Momentum{d};}
    constexpr Momentum operator"" _kg_m_per_s (unsigned long long d) {return Momentum{d};}
    
    constexpr Action operator"" _J_s (long double d) {return Action{d};}
    constexpr Action operator"" _J_s (unsigned long long d) {return Action{d};}
    
    constexpr Frequency operator"" _Hz (long double d) {return Frequency{d};}
    constexpr Frequency operator"" _Hz (unsigned long long d) {return Frequency{d};}
    
    constexpr Force operator"" _N (long double d) {return Force{d};}
    constexpr Force operator"" _N (unsigned long long d) {return Force{d};}
    
    constexpr Pressure operator"" _Pa (long double d) {return Pressure{d};}
    constexpr Pressure operator"" _Pa (unsigned long long d) {return Pressure{d};}
    
    constexpr DynamicViscosity operator"" _Pa_s (long double d) {return DynamicViscosity{d};}
    constexpr DynamicViscosity operator"" _Pa_s (unsigned long long d) {return DynamicViscosity{d};}
    
    constexpr KinematicViscosity operator"" _m2_per_s (long double d) {return KinematicViscosity{d};}
    constexpr KinematicViscosity operator"" _m2_per_s (unsigned long long d) {return KinematicViscosity{d};}
    
    constexpr Torque operator"" _N_m (long double d) {return Torque{d};}
    constexpr Torque operator"" _N_m (unsigned long long d) {return Torque{d};}
    
    constexpr Energy operator"" _J (long double d) {return Energy{d};}
    constexpr Energy operator"" _J (unsigned long long d) {return Energy{d};}
    
    constexpr Power operator"" _W (long double d) {return Power{d};}
    constexpr Power operator"" _W (unsigned long long d) {return Power{d};}
    
    constexpr ElectricCharge operator"" _C (long double d) {return ElectricCharge{d};}
    constexpr ElectricCharge operator"" _C (unsigned long long d) {return ElectricCharge{d};}
    
    constexpr Voltage operator"" _V (long double d) {return Voltage{d};}
    constexpr Voltage operator"" _V (unsigned long long d) {return Voltage{d};}
    
    constexpr Capacitance operator"" _F (long double d) {return Capacitance{d};}
    constexpr Capacitance operator"" _F (unsigned long long d) {return Capacitance{d};}
    
    constexpr Impedance operator"" _Ohm (long double d) {return Impedance{d};}
    constexpr Impedance operator"" _Ohm (unsigned long long d) {return Impedance{d};}
    
    constexpr Conductance operator"" _S (long double d) {return Conductance{d};}
    constexpr Conductance operator"" _S (unsigned long long d) {return Conductance{d};}
    
    constexpr MagneticFlux operator"" _Wb (long double d) {return MagneticFlux{d};}
    constexpr MagneticFlux operator"" _Wb (unsigned long long d) {return MagneticFlux{d};}
    
    constexpr MagneticFluxDensity operator"" _T (long double d) {return MagneticFluxDensity{d};}
    constexpr MagneticFluxDensity operator"" _T (unsigned long long d) {return MagneticFluxDensity{d};}
    
    constexpr Inductance operator"" _H (long double d) {return Inductance{d};}
    constexpr Inductance operator"" _H (unsigned long long d) {return Inductance{d};}
    
    constexpr LuminousFlux operator"" _lm (long double d) {return LuminousFlux{d/constants::pi/4};}
    constexpr LuminousFlux operator"" _lm (unsigned long long d) {return LuminousFlux{d/constants::pi/4};}
    
    constexpr Illuminance operator"" _lx (long double d) {return Illuminance{d/constants::pi/4};}
    constexpr Illuminance operator"" _lx (unsigned long long d) {return Illuminance{d/constants::pi/4};}
    
    constexpr Radioactivity operator"" _Bq (long double d) {return Radioactivity{d};}
    constexpr Radioactivity operator"" _Bq (unsigned long long d) {return Radioactivity{d};}
    
    constexpr Dose operator"" _Gy (long double d) {return Dose{d};}
    constexpr Dose operator"" _Gy (unsigned long long d) {return Dose{d};}
    
    constexpr Dose operator"" _Sv (long double d) {return Dose{d};}
    constexpr Dose operator"" _Sv (unsigned long long d) {return Dose{d};}
    
    constexpr CatalyticActivity operator"" _kat (long double d) {return CatalyticActivity{d};}
    constexpr CatalyticActivity operator"" _kat (unsigned long long d) {return CatalyticActivity{d};}
    
  }
  
  namespace constants {
    using studis::pow;
    using namespace literals;
    
    constexpr auto speed_of_light = 299792458_m / 1_s;
    constexpr auto Planck_constant = 6.62607015e-34_J * 1_s;
    constexpr auto reduced_Planck_constant = Planck_constant / (2 * pi);
    constexpr auto elementary_charge = 1.602176634e-19_C;
    constexpr auto Avogadro_constant = 6.02214076e23 / 1_mol;
    constexpr auto Boltzmann_constant = 1.380649e-23_J / 1_K;
    constexpr auto hyperfine_transition_frequency_of_Cs_133 = 9192631770_Hz;
    constexpr auto luminous_efficacy = 873_lm / 1_W;
    
    constexpr auto magnetic_flux_quantum = Planck_constant / elementary_charge / 2;
    constexpr auto conductance_quantum = 2*pow<2>(elementary_charge) / Planck_constant;
    constexpr auto Josephson_constant = 1 / magnetic_flux_quantum;
    constexpr auto von_Klitzing_constant = Planck_constant * pow<-2>(elementary_charge);
    constexpr auto Faraday_constant = elementary_charge * Avogadro_constant;
    constexpr auto gas_constant = Avogadro_constant * Boltzmann_constant;
    constexpr auto molar_gas_constant = gas_constant;
    constexpr auto universal_gas_constant = gas_constant;
    constexpr auto Stefan_Boltzmann_constant = pow<2>(pi/speed_of_light) * pow<4>(Boltzmann_constant) * pow<-3>(reduced_Planck_constant) / 60;
    constexpr auto first_radiation_constant = 2 * pi * Planck_constant * pow<2>(speed_of_light);
    constexpr auto second_radiation_constant = Planck_constant * speed_of_light / Boltzmann_constant;
    constexpr auto Wien_displacement_law_constant = 2.897771955185172e-3_K * 1_m;
    constexpr auto Wien_constant = Wien_displacement_law_constant;
    
    constexpr auto magnetic_constant = 1.25663706212e-6_H / 1_m;
    constexpr auto vacuum_permeability = magnetic_constant;
    constexpr auto electric_constant = pow<-2>(speed_of_light) / magnetic_constant;
    constexpr auto vacuum_permittivity = electric_constant;
    constexpr auto characteristic_impedance_of_vacuum = magnetic_constant * speed_of_light;
    constexpr auto gravitational_constant = 6.67430e-11_N * 1_m2 * pow<-2>(1_kg);
    constexpr auto Newtonian_constant_of_gravitation = gravitational_constant;
    constexpr auto universal_gravitational_constant = gravitational_constant;
    constexpr auto electron_mass = 9.1093837015e-31_kg;
    constexpr auto proton_mass = 1.67262192369e-27_kg;
    constexpr auto proton_electron_mass_ratio = 1836.15267343;
    constexpr auto fine_structure_constant = magnetic_constant * speed_of_light*pow<2>(elementary_charge) / Planck_constant / 2;
    constexpr auto inverse_fine_structure_constant = 1 / fine_structure_constant;
    constexpr auto Rydberg_constant = electron_mass * pow<4>(elementary_charge) / (8 * pow<2>(vacuum_permittivity) * pow<3>(Planck_constant) * speed_of_light); //decltype(fine_structure_constant*fine_structure_constant*electron_mass*speed_of_light/Planck_constant/2) {10973731.568508};
    constexpr auto atomic_mass_constant = 1.66053906660e-27_kg;
    constexpr auto atomic_mass_unit = atomic_mass_constant;
    constexpr auto Dalton = atomic_mass_constant;
    constexpr auto electron_volt = elementary_charge*1_V;
    constexpr auto inverse_of_conductance_quantum = 1 / conductance_quantum;
    constexpr auto Bohr_magneton = elementary_charge*reduced_Planck_constant/electron_mass/2;
    constexpr auto nuclear_magneton = elementary_charge*reduced_Planck_constant/proton_mass/2;
    constexpr auto Bohr_radius = reduced_Planck_constant / (electron_mass * speed_of_light * fine_structure_constant);
//     constexpr auto Planck_length = sqrt (reduced_Planck_constant*gravitational_constant*pow<-3>(speed_of_light));
//     constexpr auto Planck_mass = sqrt (reduced_Planck_constant*speed_of_light/gravitational_constant);
//     constexpr auto Planck_time = sqrt (reduced_Planck_constant*gravitational_constant*pow<-5>(speed_of_light));
//     constexpr auto Planck_charge = sqrt (2*electric_constant*Planck_constant*speed_of_light);
//     constexpr auto Planck_temperature = Planck_mass*pow<2>(speed_of_light)/Boltzmann_constant;
    
    constexpr auto minute = 60_s;
    constexpr auto hour = 60*minute;
    constexpr auto day = 24*hour;
    constexpr auto degree = pi/180;
    constexpr auto arcminute = degree/60;
    constexpr auto arcsecond = arcminute/60;
    constexpr auto hectare = 1e4_m2;
    constexpr auto litre = 1e-3_m3;
    constexpr auto tonne = 1e3_kg;
    
    constexpr auto dyne = 1_cm*1_g/pow<2>(1_s);
    constexpr auto erg = dyne*1_cm;
    constexpr auto poise = 1_g/1_cm/1_s;
    constexpr auto stokes = pow<2>(1_cm)/1_s;
    constexpr auto gal = 1_cm/pow<2>(1_s);
    constexpr auto gauss = 1e-4_T;
    constexpr auto maxwell = gauss*pow<2>(1_cm);
    
    constexpr auto standard_gravity = 9.80665_m_per_s2;
    constexpr auto standard_atmosphere = 101325_Pa;
    constexpr auto standard_state_pressure = 100000_Pa;
    constexpr auto mercury_density = 13595.1_kg_per_m3;
    
    constexpr auto inch = 2.54_cm;
    constexpr auto foot = inch*12;
    constexpr auto yard = foot*3;
    constexpr auto mile = yard*1760;
    constexpr auto nautical_mile = 1852_m;
    constexpr auto knot = nautical_mile / hour;
    constexpr auto pound = 0.45359237_kg;
    constexpr auto ounce = pound/16;
    constexpr auto pound_force = pound*standard_gravity;
    constexpr auto pound_force_per_squared_inch = pound_force/inch/inch;
    constexpr auto british_thermal_unit = 788169*foot*pound_force;
    constexpr auto thermochemical_calorie = 4184_J;
    
    constexpr auto julian_year = 365.25*day;
    constexpr auto astronomical_unit = 149597870700_m;
    constexpr auto light_year = speed_of_light*julian_year;
    constexpr auto parsec = astronomical_unit*648000/constants::pi;
    
    constexpr auto angstrom = 1e-10_m;
    constexpr auto svedberg = 1e-13_s;
    constexpr auto barn = 1e-28_m2;
    constexpr auto bar = 1e5_Pa;
    constexpr auto torr = standard_atmosphere/760;
    constexpr auto millimeter_of_mercury = mercury_density*standard_gravity*1e-3_m;
    constexpr auto watt_hour = 1_W*hour;
    constexpr auto ampere_hour = 1_A*hour;
//     constexpr auto Hartree_energy = 2*Rydberg_constant*Planck_constant*speed_of_light;
  }
  
  namespace literals {
    using studis::pow;
    
    constexpr Length operator"" _fm (long double d) {return d*1e-15_m;}
    constexpr Length operator"" _pm (long double d) {return d*1e-12_m;}
    constexpr Length operator"" _nm (long double d) {return d*1e-9_m;}
    constexpr Length operator"" _um (long double d) {return d*1e-6_m;}
    constexpr Length operator"" _mm (long double d) {return d*1e-3_m;}
    constexpr Length operator"" _km (long double d) {return d*1e+3_m;}
    
    constexpr Length operator"" _fm (unsigned long long d) {return d*1e-15_m;}
    constexpr Length operator"" _pm (unsigned long long d) {return d*1e-12_m;}
    constexpr Length operator"" _nm (unsigned long long d) {return d*1e-9_m;}
    constexpr Length operator"" _um (unsigned long long d) {return d*1e-6_m;}
    constexpr Length operator"" _mm (unsigned long long d) {return d*1e-3_m;}
    constexpr Length operator"" _km (unsigned long long d) {return d*1e+3_m;}
    
    constexpr Length operator"" _angstrom (long double d) {return d*constants::angstrom;}
    constexpr Length operator"" _micron (long double d) {return d*1.0_um;}
    constexpr Length operator"" _in (long double d) {return d*constants::inch;}
    constexpr Length operator"" _ft (long double d) {return d*constants::foot;}
    constexpr Length operator"" _yd (long double d) {return d*constants::yard;}
    constexpr Length operator"" _mile (long double d) {return d*constants::mile;}
    constexpr Length operator"" _nautical_mile (long double d) {return d*constants::nautical_mile;}
    constexpr Length operator"" _au (long double d) {return d*constants::astronomical_unit;}
    constexpr Length operator"" _ly (long double d) {return d*constants::light_year;}
    constexpr Length operator"" _kly (long double d) {return d*constants::light_year*1e3;}
    constexpr Length operator"" _Mly (long double d) {return d*constants::light_year*1e6;}
    constexpr Length operator"" _Gly (long double d) {return d*constants::light_year*1e9;}
    constexpr Length operator"" _pc (long double d) {return d*constants::parsec;}
    constexpr Length operator"" _kpc (long double d) {return d*constants::parsec*1e3;}
    constexpr Length operator"" _Mpc (long double d) {return d*constants::parsec*1e6;}
    constexpr Length operator"" _Gpc (long double d) {return d*constants::parsec*1e9;}
    
    constexpr Length operator"" _angstrom (unsigned long long d) {return d*constants::angstrom;}
    constexpr Length operator"" _micron (unsigned long long d) {return d*1.0_um;}
    constexpr Length operator"" _in (unsigned long long d) {return d*constants::inch;}
    constexpr Length operator"" _ft (unsigned long long d) {return d*constants::foot;}
    constexpr Length operator"" _yd (unsigned long long d) {return d*constants::yard;}
    constexpr Length operator"" _mile (unsigned long long d) {return d*constants::mile;}
    constexpr Length operator"" _nautical_mile (unsigned long long d) {return d*constants::nautical_mile;}
    constexpr Length operator"" _au (unsigned long long d) {return d*constants::astronomical_unit;}
    constexpr Length operator"" _ly (unsigned long long d) {return d*constants::light_year;}
    constexpr Length operator"" _kly (unsigned long long d) {return d*constants::light_year*1e3;}
    constexpr Length operator"" _Mly (unsigned long long d) {return d*constants::light_year*1e6;}
    constexpr Length operator"" _Gly (unsigned long long d) {return d*constants::light_year*1e9;}
    constexpr Length operator"" _pc (unsigned long long d) {return d*constants::parsec;}
    constexpr Length operator"" _kpc (unsigned long long d) {return d*constants::parsec*1e3;}
    constexpr Length operator"" _Mpc (unsigned long long d) {return d*constants::parsec*1e6;}
    constexpr Length operator"" _Gpc (unsigned long long d) {return d*constants::parsec*1e9;}
    
    constexpr Mass operator"" _fg (long double d) {return d*1e-18_kg;}
    constexpr Mass operator"" _pg (long double d) {return d*1e-15_kg;}
    constexpr Mass operator"" _ng (long double d) {return d*1e-12_kg;}
    constexpr Mass operator"" _ug (long double d) {return d*1e-9_kg;}
    constexpr Mass operator"" _mg (long double d) {return d*1e-6_kg;}
    
    constexpr Mass operator"" _fg (unsigned long long d) {return d*1e-18_kg;}
    constexpr Mass operator"" _pg (unsigned long long d) {return d*1e-15_kg;}
    constexpr Mass operator"" _ng (unsigned long long d) {return d*1e-12_kg;}
    constexpr Mass operator"" _ug (unsigned long long d) {return d*1e-9_kg;}
    constexpr Mass operator"" _mg (unsigned long long d) {return d*1e-6_kg;}
    
    constexpr Mass operator"" _Da (long double d) {return d*constants::Dalton;}
    constexpr Mass operator"" _kDa (long double d) {return d*constants::Dalton*1e3;}
    constexpr Mass operator"" _MDa (long double d) {return d*constants::Dalton*1e6;}
    constexpr Mass operator"" _gr (long double d) {return d*1.0_g;}
    constexpr Mass operator"" _lb (long double d) {return d*constants::pound;}
    constexpr Mass operator"" _oz (long double d) {return d*constants::ounce;}
    constexpr Mass operator"" _t (long double d) {return d*constants::tonne;}
    
    constexpr Mass operator"" _Da (unsigned long long d) {return d*constants::Dalton;}
    constexpr Mass operator"" _kDa (unsigned long long d) {return d*constants::Dalton*1e3;}
    constexpr Mass operator"" _MDa (unsigned long long d) {return d*constants::Dalton*1e6;}
    constexpr Mass operator"" _gr (unsigned long long d) {return d*1.0_g;}
    constexpr Mass operator"" _lb (unsigned long long d) {return d*constants::pound;}
    constexpr Mass operator"" _oz (unsigned long long d) {return d*constants::ounce;}
    constexpr Mass operator"" _t (unsigned long long d) {return d*constants::tonne;}
    
    constexpr Time operator"" _fs (long double d) {return d*1e-15_s;}
    constexpr Time operator"" _ps (long double d) {return d*1e-12_s;}
    constexpr Time operator"" _ns (long double d) {return d*1e-9_s;}
    constexpr Time operator"" _us (long double d) {return d*1e-6_s;}
    constexpr Time operator"" _ms (long double d) {return d*1e-3_s;}
    
    constexpr Time operator"" _fs (unsigned long long d) {return d*1e-15_s;}
    constexpr Time operator"" _ps (unsigned long long d) {return d*1e-12_s;}
    constexpr Time operator"" _ns (unsigned long long d) {return d*1e-9_s;}
    constexpr Time operator"" _us (unsigned long long d) {return d*1e-6_s;}
    constexpr Time operator"" _ms (unsigned long long d) {return d*1e-3_s;}
    
    constexpr Time operator"" _Svedberg (long double d) {return d*constants::svedberg;}
    constexpr Time operator"" _sec (long double d) {return d*1.0_s;}
    constexpr Time operator"" _min (long double d) {return d*constants::minute;}
    constexpr Time operator"" _h (long double d) {return d*constants::hour;}
    constexpr Time operator"" _hour (long double d) {return d*constants::hour;}
    constexpr Time operator"" _d (long double d) {return d*constants::day;}
    constexpr Time operator"" _day (long double d) {return d*constants::day;}
//     constexpr Time operator"" _week (long double d) {return d*constants::day*7;}
//     constexpr Time operator"" _month (long double d) {return d*constants::day*30;}
    constexpr Time operator"" _year (long double d) {return d*constants::julian_year;}
    
    constexpr Time operator"" _Svedberg (unsigned long long d) {return d*constants::svedberg;}
    constexpr Time operator"" _sec (unsigned long long d) {return d*1.0_s;}
    constexpr Time operator"" _min (unsigned long long d) {return d*constants::minute;}
    constexpr Time operator"" _h (unsigned long long d) {return d*constants::hour;}
    constexpr Time operator"" _hour (unsigned long long d) {return d*constants::hour;}
    constexpr Time operator"" _d (unsigned long long d) {return d*constants::day;}
    constexpr Time operator"" _day (unsigned long long d) {return d*constants::day;}
//     constexpr Time operator"" _week (unsigned long long d) {return d*constants::day*7;}
//     constexpr Time operator"" _month (unsigned long long d) {return d*constants::day*30;}
    constexpr Time operator"" _year (unsigned long long d) {return d*constants::julian_year;}
    
    constexpr ElectricCurrent operator"" _nA (long double d) {return d*1e-9_A;}
    constexpr ElectricCurrent operator"" _uA (long double d) {return d*1e-6_A;}
    constexpr ElectricCurrent operator"" _mA (long double d) {return d*1e-3_A;}
    constexpr ElectricCurrent operator"" _kA (long double d) {return d*1e+3_A;}
    
    constexpr ElectricCurrent operator"" _nA (unsigned long long d) {return d*1e-9_A;}
    constexpr ElectricCurrent operator"" _uA (unsigned long long d) {return d*1e-6_A;}
    constexpr ElectricCurrent operator"" _mA (unsigned long long d) {return d*1e-3_A;}
    constexpr ElectricCurrent operator"" _kA (unsigned long long d) {return d*1e+3_A;}
    
    constexpr Temperature operator"" _deg_C (long double d) {return Temperature{d+273.15};}
    constexpr Temperature operator"" _degree_Celsius (long double d) {return Temperature{d+273.15};}
    constexpr Temperature operator"" _deg_F (long double d) {return Temperature{(d+32.)*5/9.+273.15};}
    constexpr Temperature operator"" _degree_Fahrenheit (long double d) {return Temperature{(d+32.)*5/9.+273.15};}
    
    constexpr Temperature operator"" _deg_C (unsigned long long d) {return Temperature{d+273.15};}
    constexpr Temperature operator"" _degree_Celsius (unsigned long long d) {return Temperature{d+273.15};}
    constexpr Temperature operator"" _deg_F (unsigned long long d) {return Temperature{(d+32.)*5/9.+273.15};}
    constexpr Temperature operator"" _degree_Fahrenheit (unsigned long long d) {return Temperature{(d+32.)*5/9.+273.15};}
    
    constexpr AmountOfSubstance operator"" _nmol (long double d) {return d*1e-9_mol;}
    constexpr AmountOfSubstance operator"" _umol (long double d) {return d*1e-6_mol;}
    constexpr AmountOfSubstance operator"" _mmol (long double d) {return d*1e-3_mol;}
    constexpr AmountOfSubstance operator"" _kmol (long double d) {return d*1e+3_mol;}
    
    constexpr AmountOfSubstance operator"" _nmol (unsigned long long d) {return d*1e-9_mol;}
    constexpr AmountOfSubstance operator"" _umol (unsigned long long d) {return d*1e-6_mol;}
    constexpr AmountOfSubstance operator"" _mmol (unsigned long long d) {return d*1e-3_mol;}
    constexpr AmountOfSubstance operator"" _kmol (unsigned long long d) {return d*1e+3_mol;}
    
    
    constexpr Area operator"" _barn (long double d) {return d*constants::barn;}
    constexpr Area operator"" _mm2 (long double d) {return d*1_mm*1_mm;}
    constexpr Area operator"" _cm2 (long double d) {return d*1_cm*1_cm;}
    constexpr Area operator"" _km2 (long double d) {return d*1_km*1_km;}
    constexpr Area operator"" _in2 (long double d) {return d*1_in*1_in;}
    constexpr Area operator"" _ft2 (long double d) {return d*1_ft*1_ft;}
    constexpr Area operator"" _yd2 (long double d) {return d*1_yd*1_yd;}
    constexpr Area operator"" _mile2 (long double d) {return d*1_mile*1_mile;}
    constexpr Area operator"" _ha (long double d) {return d*constants::hectare;}
    constexpr Area operator"" _hectare (long double d) {return d*constants::hectare;}
    
    constexpr Area operator"" _barn (unsigned long long d) {return d*constants::barn;}
    constexpr Area operator"" _mm2 (unsigned long long d) {return d*1_mm*1_mm;}
    constexpr Area operator"" _cm2 (unsigned long long d) {return d*1_cm*1_cm;}
    constexpr Area operator"" _km2 (unsigned long long d) {return d*1_km*1_km;}
    constexpr Area operator"" _in2 (unsigned long long d) {return d*1_in*1_in;}
    constexpr Area operator"" _ft2 (unsigned long long d) {return d*1_ft*1_ft;}
    constexpr Area operator"" _yd2 (unsigned long long d) {return d*1_yd*1_yd;}
    constexpr Area operator"" _mile2 (unsigned long long d) {return d*1_mile*1_mile;}
    constexpr Area operator"" _ha (unsigned long long d) {return d*constants::hectare;}
    constexpr Area operator"" _hectare (unsigned long long d) {return d*constants::hectare;}
    
    constexpr Volume operator"" _cm3 (long double d) {return d*1_cm*1_cm*1_cm;}
    constexpr Volume operator"" _litre (long double d) {return d*10_cm*10_cm*10_cm;}
    constexpr Volume operator"" _l (long double d) {return d*1.0_litre;}
    constexpr Volume operator"" _L (long double d) {return d*1.0_litre;}
    constexpr Volume operator"" _ml (long double d) {return d*1e-3_litre;}
    constexpr Volume operator"" _mL (long double d) {return d*1e-3_litre;}
    constexpr Volume operator"" _ul (long double d) {return d*1e-6_litre;}
    constexpr Volume operator"" _uL (long double d) {return d*1e-6_litre;}
    
    constexpr Volume operator"" _cm3 (unsigned long long d) {return d*1.0_cm3;}
    constexpr Volume operator"" _litre (unsigned long long d) {return d*1.0_litre;}
    constexpr Volume operator"" _l (unsigned long long d) {return d*1.0_litre;}
    constexpr Volume operator"" _L (unsigned long long d) {return d*1.0_litre;}
    constexpr Volume operator"" _ml (unsigned long long d) {return d*1.0_ml;}
    constexpr Volume operator"" _mL (unsigned long long d) {return d*1.0_ml;}
    constexpr Volume operator"" _ul (unsigned long long d) {return d*1.0_ul;}
    constexpr Volume operator"" _uL (unsigned long long d) {return d*1.0_ul;}
    
    constexpr Density operator"" _gr_per_cm3 (long double d) {return d*1_gr/1_cm3;}
    constexpr Density operator"" _gr_per_ml (long double d) {return d*1_gr/1_ml;}
    constexpr Density operator"" _gr_per_mL (long double d) {return d*1_gr/1_mL;}
    constexpr Density operator"" _kg_per_l (long double d) {return d*1_kg/1_l;}
    constexpr Density operator"" _kg_per_L (long double d) {return d*1_kg/1_L;}
    
    constexpr Density operator"" _gr_per_cm3 (unsigned long long d) {return d*1_gr/1_cm3;}
    constexpr Density operator"" _gr_per_ml (unsigned long long d) {return d*1_gr/1_ml;}
    constexpr Density operator"" _gr_per_mL (unsigned long long d) {return d*1_gr/1_mL;}
    constexpr Density operator"" _kg_per_l (unsigned long long d) {return d*1_kg/1_l;}
    constexpr Density operator"" _kg_per_L (unsigned long long d) {return d*1_kg/1_L;}
    
    constexpr Concentration operator"" _pM (long double d) {return Concentration{d*1e-9};}
    constexpr Concentration operator"" _nM (long double d) {return Concentration{d*1e-6};}
    constexpr Concentration operator"" _uM (long double d) {return Concentration{d*1e-3};}
    constexpr Concentration operator"" _mM (long double d) {return Concentration{d};}
    constexpr Concentration operator"" _M (long double d) {return Concentration{d*1e3};}
    
    constexpr Concentration operator"" _pM (unsigned long long d) {return Concentration{d*1e-9};}
    constexpr Concentration operator"" _nM (unsigned long long d) {return Concentration{d*1e-6};}
    constexpr Concentration operator"" _uM (unsigned long long d) {return Concentration{d*1e-3};}
    constexpr Concentration operator"" _mM (unsigned long long d) {return Concentration{d};}
    constexpr Concentration operator"" _M (unsigned long long d) {return Concentration{d*1e3};}
    
    constexpr Velocity operator"" _ft_per_s (long double d) {return d*1_ft/1_s;}
    constexpr Velocity operator"" _ft_per_sec (long double d) {return d*1_ft/1_sec;}
    constexpr Velocity operator"" _km_per_hour (long double d) {return d*1_km/1_hour;}
    constexpr Velocity operator"" _mile_per_hour (long double d) {return d*1_km/1_hour;}
    constexpr Velocity operator"" _knot (long double d) {return d*constants::knot;}
    
    constexpr Velocity operator"" _ft_per_s (unsigned long long d) {return d*1_ft/1_s;}
    constexpr Velocity operator"" _ft_per_sec (unsigned long long d) {return d*1_ft/1_sec;}
    constexpr Velocity operator"" _km_per_hour (unsigned long long d) {return d*1_km/1_hour;}
    constexpr Velocity operator"" _mile_per_hour (unsigned long long d) {return d*1_km/1_hour;}
    constexpr Velocity operator"" _knot (unsigned long long d) {return d*constants::knot;}
    
//     constexpr Accelaration operator"" _cm_per_s2 (long double d) {return d*1_cm/1_s/1_s;}
    constexpr Accelaration operator"" _ft_per_s2 (long double d) {return d*1_ft/1_s/1_s;}
    constexpr Accelaration operator"" _Gal (long double d) {return d*constants::gal;}
    
//     constexpr Accelaration operator"" _cm_per_s2 (unsigned long long d) {return d*1_cm/1_s/1_s;}
    constexpr Accelaration operator"" _ft_per_s2 (unsigned long long d) {return d*1_ft/1_s/1_s;}
    constexpr Accelaration operator"" _Gal (unsigned long long d) {return d*constants::gal;}
    
    constexpr Frequency operator"" _kHz (long double d) {return d*1e3_Hz;}
    constexpr Frequency operator"" _MHz (long double d) {return d*1e6_Hz;}
    constexpr Frequency operator"" _GHz (long double d) {return d*1e9_Hz;}
    constexpr Frequency operator"" _THz (long double d) {return d*1e12_Hz;}
    
    constexpr Frequency operator"" _kHz (unsigned long long d) {return d*1.0_kHz;}
    constexpr Frequency operator"" _MHz (unsigned long long d) {return d*1.0_MHz;}
    constexpr Frequency operator"" _GHz (unsigned long long d) {return d*1.0_GHz;}
    constexpr Frequency operator"" _THz (unsigned long long d) {return d*1.0_THz;}
    
    constexpr Frequency operator"" _Bd (long double d) {return Frequency{d};}
    constexpr Frequency operator"" _kBd (long double d) {return d*1e3_Bd;}
    constexpr Frequency operator"" _MBd (long double d) {return d*1e6_Bd;}
    constexpr Frequency operator"" _GBd (long double d) {return d*1e9_Bd;}
    
    constexpr Frequency operator"" _Bd (unsigned long long d) {return Frequency{d};}
    constexpr Frequency operator"" _kBd (unsigned long long d) {return d*1e3_Bd;}
    constexpr Frequency operator"" _MBd (unsigned long long d) {return d*1e6_Bd;}
    constexpr Frequency operator"" _GBd (unsigned long long d) {return d*1e9_Bd;}
    
    constexpr Frequency operator"" _FLOPS (long double d) {return Frequency{d};}
    constexpr Frequency operator"" _kFLOPS (long double d) {return d*1e3_FLOPS;}
    constexpr Frequency operator"" _MFLOPS (long double d) {return d*1e6_FLOPS;}
    constexpr Frequency operator"" _GFLOPS (long double d) {return d*1e9_FLOPS;}
    constexpr Frequency operator"" _TFLOPS (long double d) {return d*1e12_FLOPS;}
    
    constexpr Frequency operator"" _FLOPS (unsigned long long d) {return Frequency{d};}
    constexpr Frequency operator"" _kFLOPS (unsigned long long d) {return d*1e3_FLOPS;}
    constexpr Frequency operator"" _MFLOPS (unsigned long long d) {return d*1e6_FLOPS;}
    constexpr Frequency operator"" _GFLOPS (unsigned long long d) {return d*1e9_FLOPS;}
    constexpr Frequency operator"" _TFLOPS (unsigned long long d) {return d*1e12_FLOPS;}
    
    constexpr Frequency operator"" _rpm (long double d) {return Frequency{d/60.0};}
    constexpr Frequency operator"" _fps (long double d) {return Frequency{d};}
    
    constexpr Frequency operator"" _rpm (unsigned long long d) {return Frequency{d/60.0};}
    constexpr Frequency operator"" _fps (unsigned long long d) {return Frequency{d};}
    
    constexpr Force operator"" _pN (long double d) {return d*1e-12_N;}
    constexpr Force operator"" _nN (long double d) {return d*1e-9_N;}
    constexpr Force operator"" _uN (long double d) {return d*1e-6_N;}
    constexpr Force operator"" _mN (long double d) {return d*1e-3_N;}
    constexpr Force operator"" _kN (long double d) {return d*1e+3_N;}
    
    constexpr Force operator"" _pN (unsigned long long d) {return d*1e-12_N;}
    constexpr Force operator"" _nN (unsigned long long d) {return d*1e-9_N;}
    constexpr Force operator"" _uN (unsigned long long d) {return d*1e-6_N;}
    constexpr Force operator"" _mN (unsigned long long d) {return d*1e-3_N;}
    constexpr Force operator"" _kN (unsigned long long d) {return d*1e+3_N;}
    
    constexpr Force operator"" _dyn (long double d) {return d*constants::dyne;}
    constexpr Force operator"" _dyne (long double d) {return d*constants::dyne;}
    constexpr Force operator"" _lbf (long double d) {return d*constants::pound_force;}
    
    constexpr Force operator"" _dyn (unsigned long long d) {return d*constants::dyne;}
    constexpr Force operator"" _dyne (unsigned long long d) {return d*constants::dyne;}
    constexpr Force operator"" _lbf (unsigned long long d) {return d*constants::pound_force;}
    
    constexpr Pressure operator"" _kPa (long double d) {return d*1e3_Pa;}
    constexpr Pressure operator"" _MPa (long double d) {return d*1e6_Pa;}
    constexpr Pressure operator"" _GPa (long double d) {return d*1e9_Pa;}
    
    constexpr Pressure operator"" _kPa (unsigned long long d) {return d*1e3_Pa;}
    constexpr Pressure operator"" _MPa (unsigned long long d) {return d*1e6_Pa;}
    constexpr Pressure operator"" _GPa (unsigned long long d) {return d*1e9_Pa;}
    
    constexpr Pressure operator"" _Torr (long double d) {return d*constants::torr;}
    constexpr Pressure operator"" _mTorr (long double d) {return d*constants::torr*1e-3;}
    constexpr Pressure operator"" _mmHg (long double d) {return d*constants::millimeter_of_mercury;}
    constexpr Pressure operator"" _cmHg (long double d) {return d*constants::millimeter_of_mercury/10;}
    constexpr Pressure operator"" _psi (long double d) {return d*constants::pound_force_per_squared_inch;}
    constexpr Pressure operator"" _bar (long double d) {return d*1e5_Pa;}
    constexpr Pressure operator"" _mbar (long double d) {return d*1e2_Pa;}
    constexpr Pressure operator"" _atm (long double d) {return d*constants::standard_atmosphere;}
    
    constexpr Pressure operator"" _Torr (unsigned long long d) {return d*constants::torr;}
    constexpr Pressure operator"" _mTorr (unsigned long long d) {return d*constants::torr*1e-3;}
    constexpr Pressure operator"" _mmHg (unsigned long long d) {return d*constants::millimeter_of_mercury;}
    constexpr Pressure operator"" _cmHg (unsigned long long d) {return d*constants::millimeter_of_mercury/10;}
    constexpr Pressure operator"" _psi (unsigned long long d) {return d*constants::pound_force_per_squared_inch;}
    constexpr Pressure operator"" _bar (unsigned long long d) {return d*1e5_Pa;}
    constexpr Pressure operator"" _mbar (unsigned long long d) {return d*1e2_Pa;}
    constexpr Pressure operator"" _atm (unsigned long long d) {return d*constants::standard_atmosphere;}
    
    constexpr DynamicViscosity operator"" _cP (long double d) {return d*constants::poise*1e-2;}
    constexpr DynamicViscosity operator"" _P (long double d) {return d*constants::poise;}
    
    constexpr DynamicViscosity operator"" _cP (unsigned long long d) {return d*constants::poise*1e-2;}
    constexpr DynamicViscosity operator"" _P (unsigned long long d) {return d*constants::poise;}
    
    constexpr KinematicViscosity operator"" _cSt (long double d) {return d*constants::stokes*1e-2;}
    constexpr KinematicViscosity operator"" _St (long double d) {return d*constants::stokes;}
    
    constexpr KinematicViscosity operator"" _cSt (unsigned long long d) {return d*constants::stokes*1e-2;}
    constexpr KinematicViscosity operator"" _St (unsigned long long d) {return d*constants::stokes;}
    
    constexpr Energy operator"" _kJ (long double d) {return d*1e3_J;}
    constexpr Energy operator"" _MJ (long double d) {return d*1e6_J;}
    constexpr Energy operator"" _GJ (long double d) {return d*1e9_J;}
    
    constexpr Energy operator"" _kJ (unsigned long long d) {return d*1e3_J;}
    constexpr Energy operator"" _MJ (unsigned long long d) {return d*1e6_J;}
    constexpr Energy operator"" _GJ (unsigned long long d) {return d*1e9_J;}
    
    constexpr Energy operator"" _eV (long double d) {return d*constants::electron_volt;}
    constexpr Energy operator"" _keV (long double d) {return d*constants::electron_volt*1e3;}
    constexpr Energy operator"" _MeV (long double d) {return d*constants::electron_volt*1e6;}
    constexpr Energy operator"" _GeV (long double d) {return d*constants::electron_volt*1e9;}
    constexpr Energy operator"" _erg (long double d) {return d*constants::erg;}
    constexpr Energy operator"" _Wh (long double d) {return d*constants::watt_hour;}
    constexpr Energy operator"" _kWh (long double d) {return d*constants::watt_hour*1e3;}
    constexpr Energy operator"" _MWh (long double d) {return d*constants::watt_hour*1e6;}
    constexpr Energy operator"" _BTU (long double d) {return d*constants::british_thermal_unit;}
    constexpr Energy operator"" _cal (long double d) {return d*constants::thermochemical_calorie;}
    constexpr Energy operator"" _kcal (long double d) {return d*constants::thermochemical_calorie*1e3;}
    
    constexpr Energy operator"" _eV (unsigned long long d) {return d*constants::electron_volt;}
    constexpr Energy operator"" _keV (unsigned long long d) {return d*constants::electron_volt*1e3;}
    constexpr Energy operator"" _MeV (unsigned long long d) {return d*constants::electron_volt*1e6;}
    constexpr Energy operator"" _GeV (unsigned long long d) {return d*constants::electron_volt*1e9;}
    constexpr Energy operator"" _erg (unsigned long long d) {return d*constants::erg;}
    constexpr Energy operator"" _Wh (unsigned long long d) {return d*constants::watt_hour;}
    constexpr Energy operator"" _kWh (unsigned long long d) {return d*constants::watt_hour*1e3;}
    constexpr Energy operator"" _MWh (unsigned long long d) {return d*constants::watt_hour*1e6;}
    constexpr Energy operator"" _BTU (unsigned long long d) {return d*constants::british_thermal_unit;}
    constexpr Energy operator"" _cal (unsigned long long d) {return d*constants::thermochemical_calorie;}
    constexpr Energy operator"" _kcal (unsigned long long d) {return d*constants::thermochemical_calorie*1e3;}
    
    constexpr Power operator"" _nW (long double d) {return d*1e-9_W;}
    constexpr Power operator"" _uW (long double d) {return d*1e-6_W;}
    constexpr Power operator"" _mW (long double d) {return d*1e-3_W;}
    constexpr Power operator"" _kW (long double d) {return d*1e3_W;}
    constexpr Power operator"" _MW (long double d) {return d*1e6_W;}
    constexpr Power operator"" _GW (long double d) {return d*1e9_W;}
    
    constexpr Power operator"" _nW (unsigned long long d) {return d*1e-9_W;}
    constexpr Power operator"" _uW (unsigned long long d) {return d*1e-6_W;}
    constexpr Power operator"" _mW (unsigned long long d) {return d*1e-3_W;}
    constexpr Power operator"" _kW (unsigned long long d) {return d*1e3_W;}
    constexpr Power operator"" _MW (unsigned long long d) {return d*1e6_W;}
    constexpr Power operator"" _GW (unsigned long long d) {return d*1e9_W;}
    
    constexpr ElectricCharge operator"" _pC (long double d) {return d*1e-12_C;}
    constexpr ElectricCharge operator"" _nC (long double d) {return d*1e-9_C;}
    constexpr ElectricCharge operator"" _uC (long double d) {return d*1e-6_C;}
    constexpr ElectricCharge operator"" _mC (long double d) {return d*1e-3_C;}
    
    constexpr ElectricCharge operator"" _pC (unsigned long long d) {return d*1e-12_C;}
    constexpr ElectricCharge operator"" _nC (unsigned long long d) {return d*1e-9_C;}
    constexpr ElectricCharge operator"" _uC (unsigned long long d) {return d*1e-6_C;}
    constexpr ElectricCharge operator"" _mC (unsigned long long d) {return d*1e-3_C;}
    
    constexpr ElectricCharge operator"" _mAh (long double d) {return d*constants::ampere_hour*1e-3;}
    constexpr ElectricCharge operator"" _Ah (long double d) {return d*constants::ampere_hour;}
    
    constexpr ElectricCharge operator"" _mAh (unsigned long long d) {return d*constants::ampere_hour*1e-3;}
    constexpr ElectricCharge operator"" _Ah (unsigned long long d) {return d*constants::ampere_hour;}
    
    constexpr Voltage operator"" _uV (long double d) {return d*1e-6_V;}
    constexpr Voltage operator"" _mV (long double d) {return d*1e-3_V;}
    constexpr Voltage operator"" _kV (long double d) {return d*1e+3_V;}
    constexpr Voltage operator"" _MV (long double d) {return d*1e+6_V;}
    
    constexpr Voltage operator"" _uV (unsigned long long d) {return d*1e-6_V;}
    constexpr Voltage operator"" _mV (unsigned long long d) {return d*1e-3_V;}
    constexpr Voltage operator"" _kV (unsigned long long d) {return d*1e+3_V;}
    constexpr Voltage operator"" _MV (unsigned long long d) {return d*1e+6_V;}
    
    constexpr Capacitance operator"" _pF (long double d) {return d*1e-12_F;}
    constexpr Capacitance operator"" _nF (long double d) {return d*1e-9_F;}
    constexpr Capacitance operator"" _uF (long double d) {return d*1e-6_F;}
    constexpr Capacitance operator"" _mF (long double d) {return d*1e-3_F;}
    
    constexpr Capacitance operator"" _pF (unsigned long long d) {return d*1e-12_F;}
    constexpr Capacitance operator"" _nF (unsigned long long d) {return d*1e-9_F;}
    constexpr Capacitance operator"" _uF (unsigned long long d) {return d*1e-6_F;}
    constexpr Capacitance operator"" _mF (unsigned long long d) {return d*1e-3_F;}
    
    constexpr Impedance operator"" _uOhm (long double d) {return d*1e-6_Ohm;}
    constexpr Impedance operator"" _mOhm (long double d) {return d*1e-3_Ohm;}
    constexpr Impedance operator"" _kOhm (long double d) {return d*1e+3_Ohm;}
    constexpr Impedance operator"" _MOhm (long double d) {return d*1e+6_Ohm;}
    constexpr Impedance operator"" _GOhm (long double d) {return d*1e+9_Ohm;}
    
    constexpr Impedance operator"" _uOhm (unsigned long long d) {return d*1e-6_Ohm;}
    constexpr Impedance operator"" _mOhm (unsigned long long d) {return d*1e-3_Ohm;}
    constexpr Impedance operator"" _kOhm (unsigned long long d) {return d*1e+3_Ohm;}
    constexpr Impedance operator"" _MOhm (unsigned long long d) {return d*1e+6_Ohm;}
    constexpr Impedance operator"" _GOhm (unsigned long long d) {return d*1e+9_Ohm;}
    
    constexpr MagneticFlux operator"" _nWb (long double d) {return d*1e-9_Wb;}
    constexpr MagneticFlux operator"" _uWb (long double d) {return d*1e-6_Wb;}
    constexpr MagneticFlux operator"" _mWb (long double d) {return d*1e-3_Wb;}
    
    constexpr MagneticFlux operator"" _nWb (unsigned long long d) {return d*1e-9_Wb;}
    constexpr MagneticFlux operator"" _uWb (unsigned long long d) {return d*1e-6_Wb;}
    constexpr MagneticFlux operator"" _mWb (unsigned long long d) {return d*1e-3_Wb;}
    
    constexpr MagneticFlux operator"" _Mx (long double d) {return d*constants::maxwell;}
    
    constexpr MagneticFlux operator"" _Mx (unsigned long long d) {return d*constants::maxwell;}
    
    constexpr MagneticFluxDensity operator"" _uT (long double d) {return d*1e-6_T;}
    constexpr MagneticFluxDensity operator"" _mT (long double d) {return d*1e-3_T;}
    
    constexpr MagneticFluxDensity operator"" _uT (unsigned long long d) {return d*1e-6_T;}
    constexpr MagneticFluxDensity operator"" _mT (unsigned long long d) {return d*1e-3_T;}
    
    constexpr MagneticFluxDensity operator"" _mG (long double d) {return d*constants::gauss*1e-3;}
    constexpr MagneticFluxDensity operator"" _G (long double d) {return d*constants::gauss;}
    
    constexpr MagneticFluxDensity operator"" _mG (unsigned long long d) {return d*constants::gauss*1e-3;}
    constexpr MagneticFluxDensity operator"" _G (unsigned long long d) {return d*constants::gauss;}
    
    constexpr Inductance operator"" _uH (long double d) {return d*1e-6_H;}
    constexpr Inductance operator"" _mH (long double d) {return d*1e-3_H;}
    
    constexpr Inductance operator"" _uH (unsigned long long d) {return d*1e-6_H;}
    constexpr Inductance operator"" _mH (unsigned long long d) {return d*1e-3_H;}
    
  }
  
}

#endif
