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

#include <Eigen/Core>
#include "units.h"

namespace Eigen {

template<typename R, typename C>
struct NumTraits<Units::Unit<R,C>>
        : NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
{
  typedef Units::Unit<R,C> Real;
  typedef Units::Unit<R,C> NonInteger;
  typedef Units::Unit<R,C> Nested;

  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 1,
    ReadCost = 1,
    AddCost = 3,
    MulCost = 3
  };
};

// Unit = Unit * double
template<typename R, typename C>
struct ScalarBinaryOpTraits<Units::Unit<R,C>,double,internal::scalar_product_op<Units::Unit<R,C>,double>> {
  typedef Units::Unit<R,C> ReturnType;
};

// don't use this, use IgnoreUnits instead
#if 0 // not working  error: ambiguous template.  Colliding with built in trait of T = T / T.
// double = Unit / Unit
template<typename R1, typename R2, typename C>
struct ScalarBinaryOpTraits<Units::Unit<R1,C>,Units::Unit<R2,C>,internal::scalar_quotient_op<Units::Unit<R1,C>,Units::Unit<R2,C>>> {
  typedef double ReturnType;
};
#elif 0
template<>  // need to be very specific
struct ScalarBinaryOpTraits<Units::Meters,Units::Meters,internal::scalar_quotient_op<Units::Meters,Units::Meters>> {
  typedef double ReturnType;
};
#endif

// Speed = Distance / Time
template<typename R1, typename R2>
struct ScalarBinaryOpTraits<
        Units::Unit<Units::Distance, R1>,
        Units::Unit<Units::Time, R2>,
        internal::scalar_quotient_op<Units::Unit<Units::Distance, R1>,Units::Unit<Units::Time,R2>>>
{
  typedef decltype(Units::Unit<Units::Distance,R1>{}/Units::Unit<Units::Time,R2>{}) ReturnType;
};

}

#if 0
// function overrides
// implementation of functions to allow certain eigen functions (like norm)
// TODO: these break some Units philosophies.  Consider removing and using IgnoreUnits instead
namespace Units{
template<typename R, typename C>
Units::Unit<R,C> sqrt(Unit<R,C> u){ return Unit<R,C>(std::sqrt(u.get()));}
template<typename R, typename C>
Units::Unit<R,C> operator*(Unit<R,C> a, Unit<R,C> b){ return Unit<R,C>(a.get()*b.get());}
}
#endif

// IgnoreUnits - maps or casts to double without units
// Last resort function to force eigen to ignore units and treat as double (use carefully)
// readwrite version
template<int H, int W, typename R, typename T>
Eigen::Map<Eigen::Matrix<double,H,W>> IgnoreUnits(Eigen::Matrix<Units::Unit<R,T>,H,W>& m){
  return Eigen::Map<Eigen::Matrix<double,H,W>>((double*)m.data(),m.rows(),m.cols());
}
// readonly version
template<int H, int W, typename R, typename T>
Eigen::Map<const Eigen::Matrix<double,H,W>> IgnoreUnits(const Eigen::Matrix<Units::Unit<R,T>,H,W>& m){
  return Eigen::Map<const Eigen::Matrix<double,H,W>>((double*)m.data(),m.rows(),m.cols());
}
// rvalue version
template<typename D>
auto IgnoreUnits(const Eigen::MatrixBase<D>& m) -> decltype(m.template cast<double>()) {
  return m.template cast<double>();
}

