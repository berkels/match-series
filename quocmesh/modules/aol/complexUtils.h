#ifndef __COMPLEXUTILS_H
#define __COMPLEXUTILS_H

#include <aol.h>

namespace aol {

//! Mathematical constant i
const std::complex<double> I ( 0, 1 );

//! Allows int * complex<double>
template <class TypeA, class TypeB>
std::complex<TypeB> operator * ( TypeA a, std::complex<TypeB> b ) {
  return std::complex<TypeB> ( static_cast<TypeB> ( a ) ) * b;
}

//! Allows complex<double> * int
template <class TypeA, class TypeB>
std::complex<TypeA> operator * ( std::complex<TypeA> a, TypeB b ) {
  return a * std::complex<TypeA> ( static_cast<TypeA> ( b ) );
}

//! Allows int + complex<double>
template <class TypeA, class TypeB>
std::complex<TypeB> operator + ( TypeA a, std::complex<TypeB> b ) {
  return std::complex<TypeB> ( static_cast<TypeB> ( a ) ) + b;
}

//! Allows complex<double> + int
template <class TypeA, class TypeB>
std::complex<TypeA> operator + ( std::complex<TypeA> a, TypeB b ) {
  return a + std::complex<TypeA> ( static_cast<TypeA> ( b ) );
}

//! Allows int - complex<double>
template <class TypeA, class TypeB>
std::complex<TypeB> operator - ( TypeA a, std::complex<TypeB> b ) {
  return std::complex<TypeB> ( static_cast<TypeB> ( a ) ) - b;
}

//! Allows complex<double> - int
template <class TypeA, class TypeB>
std::complex<TypeA> operator - ( std::complex<TypeA> a, TypeB b ) {
  return a - std::complex<TypeA> ( static_cast<TypeA> ( b ) );
}

//! Allows int / complex<double>
template <class TypeA, class TypeB>
std::complex<TypeB> operator / ( TypeA a, std::complex<TypeB> b ) {
  return std::complex<TypeB> ( static_cast<TypeB> ( a ) ) / b;
}

//! Allows complex<double> / int
template <class TypeA, class TypeB>
std::complex<TypeA> operator / ( std::complex<TypeA> a, TypeB b ) {
  return a / std::complex<TypeA> ( static_cast<TypeA> ( b ) );
}

}

// Implicit use in operation of complex with int
using aol::operator*;
using aol::operator+;
using aol::operator-;
using aol::operator/;

#endif
