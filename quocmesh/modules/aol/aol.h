#ifndef __AOL_H
#define __AOL_H

#ifdef USE_SSE
#include <xmmintrin.h>
#include <emmintrin.h>
#endif

// C standard libraries
#include <cmath>
#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <memory>

// STL
#include <complex>
#include <limits>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <algorithm>
#include <deque>      // necessary?
#include <functional> // necessary?
#include <iterator>   // necessary?
#include <list>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <vector>

#include <typeinfo>

#ifdef USE_DUMA
#include <malloc.h>
#include <duma.h>
#include <new>
#include <dumapp.h>
#endif

// Something that defines the namespace std has to be included
// before you can use the namespace in VC++.
using namespace std;

#include <qmException.h>
#include <platformDependent.h>

namespace aol {

//! needed as base class if dynamic_cast is used.
class Obj {
public:

  virtual ~Obj() {}

  virtual void dummy() {}
};

//! Provide appropriate zero for different data types, not only standard ones
template <class T> struct ZTrait {
  static const T zero;
};

//! Provide appropriate one for different data types, not only standard ones
template <class T> struct ZOTrait : public ZTrait<T> {
  static const T one;
};

/** Provides NaN and Inf as appropriate float types (where applicable) in addition to zero and one
 *  \attention The NumberTrait contains static members that must not be used to initialize other static members (order of initialization is unclear). Use get* methods instead.
 *  \attention Do not compare values to NaN because NaN == NaN is false. Use aol::isNaN to check whether variable is NaN.
 *  \author Schwen
 */
template <class T> struct NumberTrait : public ZOTrait<T> {
  static const T NaN;
  static const T Inf;
  static const T pi;

  // static T getZero ( ); // not implemented yet. necessary?
  // static T getOne ( );  // not implemented yet. necessary?
  static T getNaN ( );
  static T getInf ( );
  static T getPi ( );
};

/**
 * Provides an appropriate initial value for a variable that is used to determine the maximum
 * of a set of numbers, see aol::Vector::getMaxValueMasked for an example.
 *
 * \author Berkels
 */
template <class T> struct MaxInitializerTrait {
  static const T MaxInitializer;
};

//! defines aliasses for certain types, which can be used in algebraic computations
enum RealTypeAlias {
  FLOAT,
  DOUBLE,
  LONG_DOUBLE
};

/** Trait to obtain a proper floating point type for abritrary types.
 *  It maps all types other than double and long double to float.
 **/
template <class T> struct RealTrait {
  typedef float RealType;
  static const RealTypeAlias ALIAS = FLOAT;
};

template <> struct RealTrait<double> {
  typedef double RealType;
  static const RealTypeAlias ALIAS = DOUBLE;
};

template <> struct RealTrait<long double> {
  typedef long double RealType;
  static const RealTypeAlias ALIAS = LONG_DOUBLE;
};

/**
 * Converts a typename into a string with the name of the typename.
 *
 * \author Berkels
 */
template <typename T>
struct TemplateToStringTrait {
  static const char* GetString();
};

/** Possible handling types for overflows. @see ScalarArray<QC_2D>::save()
 */
enum OverflowHandlingType {
  SCALE,
  CLIP,
  REFLECT,
  CLIP_THEN_SCALE
};

enum ExitStatus {
  SUCCESS,
  FAILURE
};

enum OperatorType {
  ONTHEFLY,
  ASSEMBLED
};

//! CopyFlag is for use with aol::Vector, qc::Array etc. and other structures in future ...
enum CopyFlag {
  FLAT_COPY,      // flat or shallow copy: same data pointer
  DEEP_COPY,      // deep copy: copy data
  STRUCT_COPY,    // copy structure of (Multi)Vector or other block structures
  STRUCT_COPY_UNINIT // STRUCT_COPY, but do not initialize vector (e.g. because it will be overwritten immediately)
};

//! IncludeWriteMode is for setting the correct masked matrix apply mode, when using Dirichlet bc.
enum IncludeWriteMode {
  INCLUDE_INT_WRITE_INT,    // include only interior nodes and write only to rows of inner nodes
  INCLUDE_BD_WRITE_INT,     // include only boundary nodes and write only to rows of inner nodes
  INCLUDE_ALL_WRITE_INT,    // include all nodes and write only to rows of inner nodes
  INCLUDE_INT_WRITE_ALL,    // include only interior nodes and write to all rows
  INCLUDE_ALL_WRITE_ALL,    // include all nodes and write to all rows (default case, see below)
  INCLUDE_WRITE_DEFAULT = INCLUDE_ALL_WRITE_ALL
};

enum Verbosity {
  AOL_QUIET,
  AOL_VERBOSE    // else conflict with -DVERBOSE
};

enum GridGlobalIndexMode {
  QUOC_GRID_INDEX_MODE,
  DT_GRID_INDEX_MODE,
  QUOC_ADAPTIVEGRID_INDEX_MODE
};

void assertAlert(string assertion, string file, int line, string message = "");

#ifdef DEBUG
#define QUOC_ASSERT(Expression) if ( !(Expression) ) { aol::assertAlert ( #Expression, __FILE__, __LINE__, "" ); }
#else
#define QUOC_ASSERT(Expression)
#endif

/** This struct, only specialized for true, will produce a compiler error for expressions that can be evaluated to false at compile time.
 *  Example: aol::AssertAtCompileTime < !(std::numeric_limits<RealType>::is_integer) >();
 */
template<bool> struct AssertAtCompileTime;
template<> struct AssertAtCompileTime<true> {};

/** Similar to AssertAtCompileTime, this struct will produce a compiler error if the two types given are not identical.
 *  Example: aol::AssertAtCompileTimeThatTypesAreEqual<DataType,RealType>();
 */
template<typename,typename> struct AssertAtCompileTimeThatTypesAreEqual;
template<typename T> struct AssertAtCompileTimeThatTypesAreEqual<T,T> {};

// \brief These functors are used in the matrix/vector multiplication to mask out entries in the input vector
struct BitMaskFunctorTrue { inline bool operator() ( bool /*expr*/ ) const { return true; } };
struct BitMaskFunctorFalse { inline bool operator() ( bool /*expr*/ ) const { return false; } };
struct BitMaskFunctorIdentity { inline bool operator() ( bool expr ) const { return expr; } };
struct BitMaskFunctorNegate { inline bool operator() ( bool expr ) const { return !expr; } };

/**
 * \brief Just returns the argument.
 *
 * Since the implementation is inside aol.cpp, the compiler can't inline this function
 * (at least if no link time optimization is used). This hopefully ensures that the
 * argument is forced out of the CPU register and only saved with the accuracy of the
 * template T.
 *
 * \author Berkels
 */
template<class T> T returnArgument ( const T a );

/** Check for NaN. Returns true if NaN. Note that NaN != NaN. */
template<class T> inline bool isNaN ( const T a ) { return (a != a); }

/** Check for infinite values. Returns true if +Inf or -Inf */
template<class T> inline bool isInf ( const T a ) { return ( a != 0 && a == 2*a ); }

/** Check for Inf and NaN. Returns true if value is neither +- Inf nor NaN.
 *  Note that Inf == Inf and -Inf == -Inf but NaN != NaN.
 *  Note that a method isfinite or finite is available on some platforms but this implementation is platform independent. */
 template<class T> inline bool isFinite ( const T a ) { T aa = returnArgument ( a ); return ( !isNaN ( aa ) && !isInf ( aa ) ); }

namespace { // invisible outside this file

// helper class for partial specialization
template <typename Type, bool is_signed> class OverflowChecker {};
template <typename Type> class OverflowChecker<Type,true> {
 public:
  static inline bool sum (Type a, Type b) {
    aol::AssertAtCompileTime<std::numeric_limits<Type>::is_integer>();
    aol::AssertAtCompileTime<std::numeric_limits<Type>::is_signed>();
    const Type max = std::numeric_limits<Type>::max();
    const Type min = std::numeric_limits<Type>::min();
    // If b = 0, the result always fits
    // Given the preconditions, the subtraction cannot overflow
    if ( b > 0 ) { if ( a > max - b ) return false; }
    if ( b < 0 ) { if ( a < min - b ) return false; }
    return true;
  };
  static inline bool product (Type a, Type b) {
    aol::AssertAtCompileTime<std::numeric_limits<Type>::is_integer>();
    aol::AssertAtCompileTime<std::numeric_limits<Type>::is_signed>();
    const Type max = (std::numeric_limits<Type>::max)();
    const Type min = (std::numeric_limits<Type>::min)();
    if (a > 0) {
      if (b > 0) {
        if (a > max/b) return false; // a, b > 0
      } else {
        if (b < min/a) return false; // a > 0 >= b
      }
    } else {
      if (b > 0) {
        if (a < min/b) return false; // b > 0 >= a

      } else {
        if ( (a != 0) && (b < max/a) ) return false; // 0 >= a, b
      }
    }

    return true;
  };
};

// helper class for partial specialization for unsigned types to avoid warnings
template <typename Type> class OverflowChecker<Type,false> {
 public:
  static inline bool sum (Type a, Type b) {
    aol::AssertAtCompileTime<std::numeric_limits<Type>::is_integer>();
    aol::AssertAtCompileTime<!std::numeric_limits<Type>::is_signed>();
    const Type max = std::numeric_limits<Type>::max();
    // If b = 0, the result always fits
    // Given the preconditions, the subtraction cannot overflow
    if ( b > 0 ) { if ( a > max - b ) return false; }
    return true;
  };
  static inline bool product (Type a, Type b) {
    aol::AssertAtCompileTime<std::numeric_limits<Type>::is_integer>();
    aol::AssertAtCompileTime<!std::numeric_limits<Type>::is_signed>();
    const Type max = std::numeric_limits<Type>::max();
    if ((b != 0) && a > max/b) return false;
    return true;
  };
};

}

//! checks whether a+b fits into a variable of type Type, Type has to be an integral type
//! (for floating point types, +/-inf are considered acceptable and can reasonalby be checked via isFinite)
template <typename Type> bool inline sumWillFit ( const Type a, const Type b ) {
  return OverflowChecker<Type,std::numeric_limits<Type>::is_signed>::sum (a, b);
}

//! checks whether a*b fits into a variable of type Type, Type has to be an integral type
//! (for floating point types, +/-inf are considered acceptable and can reasonalby be checked via isFinite)
template <typename Type> bool inline productWillFit ( const Type a, const Type b ) {
  return OverflowChecker<Type,std::numeric_limits<Type>::is_signed>::product (a, b);
}

namespace { // invisible outside this file

// helper class for partial specialization
template <typename Type, bool is_integer> class PowerOverflowChecker {};
template <typename Type> class PowerOverflowChecker<Type,true> {
 public:
  static inline Type sqr (Type a) {
    aol::AssertAtCompileTime<std::numeric_limits<Type>::is_integer>();
#ifdef DEBUG
    if ( ! ( aol::productWillFit<Type> ( a, a ) ) ) {
      cerr << a << "^2 != " << static_cast<Type>(a * a) << endl;
      throw aol::Exception ( "aol::Sqr: Overflow detected.", __FILE__, __LINE__ );
    }
#endif
    return ( a * a );
  };
  static inline Type cub (Type a) {
    aol::AssertAtCompileTime<std::numeric_limits<Type>::is_integer>();
#ifdef DEBUG
    if ( ( ! ( aol::productWillFit<Type> ( a, a ) ) || ! ( aol::productWillFit<Type> ( static_cast<Type>( a*a ), a ) ) ) ) { // cast is necessary due to short * short = int
      cerr << a << "^3 != " << static_cast<Type>(a * a * a) << endl;
      throw aol::Exception ( "aol::Cub: Overflow detected.", __FILE__, __LINE__ );
    }
#endif
    return ( a * a * a );
  };
  static inline Type qrt (Type a) {
    aol::AssertAtCompileTime<std::numeric_limits<Type>::is_integer>();
#ifdef DEBUG
    if ( ( ! ( aol::productWillFit<Type> ( a, a ) ) || ! ( aol::productWillFit<Type> ( static_cast<Type>( a*a ), static_cast<Type>( a*a ) ) ) ) ) { // cast is necessary due to short * short = int
      cerr << a << "^4 != " << static_cast<Type>(a * a * a * a) << endl;
      throw aol::Exception ( "aol::Qrt: Overflow detected.", __FILE__, __LINE__ );
    }
#endif
    return ( a * a * a * a);
  };
};
// helper class for partial specialization
template <typename Type> class PowerOverflowChecker<Type,false> {
 public:
  static inline Type sqr (Type a) {
    aol::AssertAtCompileTime<!std::numeric_limits<Type>::is_integer>();
    return ( a * a );
  };
  static inline Type cub (Type a) {
    aol::AssertAtCompileTime<!std::numeric_limits<Type>::is_integer>();
    return ( a * a * a );
  };
  static inline Type qrt (Type a) {
    aol::AssertAtCompileTime<!std::numeric_limits<Type>::is_integer>();
    return ( a * a * a * a);
  };
};
}

/** Returns minimum of two values. */
template<class T> inline T Min ( const T a, const T b ) { return ( ( a < b ) ? a : b ); }
/** Returns minimum of three values -- for more than that you may want to keep them in a Vector and use its minimum method */
template<class T> inline T Min ( const T a, const T b, const T c ) { return ( aol::Min ( aol::Min ( a, b ), c ) ); }

/** Returns maximum of two values. */
template<class T> inline T Max ( const T a, const T b ) { return ( ( a < b ) ? b : a ); }
/** Returns maximum of three values */
template<class T> inline T Max ( const T a, const T b, const T c ) { return ( aol::Max ( aol::Max ( a, b ), c ) ); }

/** Returns Value clamped into [Min,Max]. */
template<class T> inline T Clamp ( const T Value, const T Min, const T Max ) { return ( aol::Max ( aol::Min ( Value, Max ), Min ) ); }

/** Returns square of a. */
template<class T> inline T Sqr ( const T a ) {
  return PowerOverflowChecker<T,std::numeric_limits<T>::is_integer>::sqr (a);
}

/** Returns cubic of a.  */
template<class T> inline T Cub ( const T a ) {
  return PowerOverflowChecker<T,std::numeric_limits<T>::is_integer>::cub (a);
}

/** Returns quartic of a.  */
template<class T> inline T Qrt ( const T a ) {
  return PowerOverflowChecker<T,std::numeric_limits<T>::is_integer>::qrt (a);
}

/** Returns fourth power of Arg computed by squaring twice */
template<typename T> inline T FourthPower ( const T Arg ) {
  return ( aol::Sqr ( aol::Sqr ( Arg ) ) );
}

/** Evaluates a smooth approximation of the max(a,b) function **/
template<class RealType> inline RealType SmoothMax ( const RealType a, const RealType b, const RealType eps = 1e-6 )  { return b + .5 * ( sqrt( aol::Sqr(a-b) + eps ) + (a-b) ); }

/** Evaluates derivative w.r.t. to first argument of the \a SmoothMax function **/
template<class RealType> inline RealType SmoothMaxDiff ( const RealType a, const RealType b, const RealType eps = 1e-6 )  { return .5 * ( (a-b)/sqrt( aol::Sqr(a-b) + eps ) + aol::ZOTrait<RealType>::one ); }

/** Evaluates a smooth approximation of the min(a,b) function **/
template<class RealType> inline RealType SmoothMin ( const RealType a, const RealType b, const RealType eps = 1e-6 )  { return a - .5 * ( sqrt( aol::Sqr(a-b) + eps ) + (a-b) ); }

/** Evaluates derivative w.r.t. to first argument of the \a SmoothMin function **/
template<class RealType> inline RealType SmoothMinDiff ( const RealType a, const RealType b, const RealType eps = 1e-6 )  { return -.5 * ( (a-b)/sqrt( aol::Sqr(a-b) + eps ) - aol::ZOTrait<RealType>::one ); }

/** Returns absolute value of a */
template<class T> inline T Abs ( const T a ) { return ( ( a < aol::ZTrait<T>::zero ) ? -a : a ); }
// Use native abs where available
template<> inline int          Abs ( const int a )          { return std::abs(a); }
template<> inline long         Abs ( const long a )         { return std::abs(a); }
template<> inline float        Abs ( const float a )        { return std::abs(a); }
template<> inline double       Abs ( const double a )       { return std::abs(a); }
template<> inline long double  Abs ( const long double a )  { return std::abs(a); }
// Don't apply the unary minus operator to unsigned types.
template<> inline unsigned int Abs ( const unsigned int a ) { return a; }
template<> inline uint64_t     Abs ( const uint64_t a )     { return a ; }

/** symmetrical round to nearest int */
template <typename T> int Rint ( const T A ) { return ( A >= 0 ) ? static_cast<int> ( A + 0.5 ) : static_cast<int> ( A - 0.5 ); }

//! Signum function template, attention: signum1at0 ( 0 ) = 1
template<class T> inline T signum1at0 ( const T x ) {
  if (x >= aol::ZTrait<T>::zero) return 1;
  else if (x < aol::ZTrait<T>::zero) return -1;
  return x; // NaN
}

//! Signum function template, signum ( 0 ) = 0
template<class T> inline T signum ( const T x ) {
  if (x > aol::ZTrait<T>::zero) return 1;
  else if (x == aol::ZTrait<T>::zero) return 0;
  else if (x < aol::ZTrait<T>::zero) return -1;
  return x; // NaN
}

//! std::modf replacement, necessary since std::modf crashes under some MinGW versions
template<class RealType> RealType Modf ( RealType x, RealType* iptr ) {
  const int sign = aol::signum1at0 ( x );
  if ( aol::isInf ( x ) ) {
    *iptr = x;
    return sign * aol::ZOTrait<RealType>::zero;
  }
  // Handle signed zero.
  else if ( x == 0 ) {
    *iptr = x;
    return x;
  }
  *iptr = sign * std::floor ( sign * x );
  return x - *iptr;
}

/** checks if pointer exists before delete **/
template<typename T> inline void intelligentDelete ( T& pointer ) {
  if ( pointer ) {
    delete pointer;
    pointer = NULL;
  }
}

/** checks if pointer exists before delete (array version) **/
template<typename T> inline void intelligentDeleteArray ( T& pointer ) {
  if ( pointer ) {
    delete [] pointer;
    pointer = NULL;
  }
}

/** Conversion from degrees to radians */
template<class T> inline T DegreesToRadians ( const T a ) { return a / 180. * aol::NumberTrait<T>::pi; }

/** Conversion from radians to degrees */
template<class T> inline T RadiansToDegrees ( const T a ) { return a * 180. / aol::NumberTrait<T>::pi; }

/** Returns approximation of the Gauss error function (which is not provided as erf prior to C99).
 *  (Implementation of Abramowitz and Stegun: Handbook of Mathematical Functions, formula 7.1.26; absolute error is below 1.5e-7)
 */
template<class RealType>
RealType Erf ( const RealType x ) {
  const RealType
    a[5] = { 0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429 },
    p = 0.3275911,
    xa = fabs(x),
    t = 1.0 / ( 1.0 + p*xa ),
    y = 1.0 - ( ( ( ( ( a[4] * t + a[3]) * t ) + a[2]) * t + a[1] ) * t + a[0] ) * t * exp( - aol::Sqr ( xa ) );

  return ( aol::signum ( x ) * y );
}

#ifdef DEBUG
#define HIDE_PARAMETER_IN_OPT_MODE(parameter) parameter
#else
#define HIDE_PARAMETER_IN_OPT_MODE(parameter)
#endif

#ifndef _OPENMP
#define NON_PARALLEL_STATIC static
#else
#define NON_PARALLEL_STATIC
#endif

/**
 * \brief Calculates (Base)^(Exponent), where Base is an integer and Exponent a non negative integer.
 *
 * In contrast to the standard pow function, this one takes integers as arguments and returns an integer.
 * This way you can calculate i^j without any casts or warnings.
 *
 * \author Berkels
 */
inline int Pow( const int Base, const int Exponent ) {
  int power = 1;
  for ( int i = 0; i < Exponent; ++i )
    power *= Base;
  return power;
}

/**
 * \brief Checks whether the specified non-negative integer is a power of two
 *
 * \author Mevenkamp
 */
inline bool isPowerOfTwo ( const int I ) {
  return ((I != 0) && !(I & (I - 1)));
}

/**
 * \brief Calculates the largest Exponent such that (Factor)^(Exponent) is still a divisor of Number.
 *
 * \author Berkels
 */
inline int getFactorExponent ( const int Number, const int Factor ) {
  int exponent = 0;
  int number = Number;

  while ( ( number % Factor ) == 0 ) {
    number = number / Factor;
    ++exponent;
  };

  return exponent;
}

//! Functor for formatting numbers
//! @ingroup Output
class Format {
public:
  //! Initialize format parameters
  Format () {}
  virtual ~Format () {}

  //! Format floating point numbers
  virtual string operator () ( const long double val ) const {
    ostringstream oss;
    oss << val;
    return oss.str ();
  }

  string operator () ( const float  val ) const { return operator () ( static_cast<long double> ( val ) ); }
  string operator () ( const double val ) const { return operator () ( static_cast<long double> ( val ) ); }

  //! Format signed integer numbers
  virtual string operator () ( const int64_t val ) const {
    ostringstream oss;
    oss << val;
    return oss.str ();
  }

  string operator () ( const short val ) const { return operator () ( static_cast<int64_t> ( val ) ); }
  string operator () ( const int   val ) const { return operator () ( static_cast<int64_t> ( val ) ); }

  //! Format unsigned integer numbers
  virtual string operator () ( const uint64_t val ) const {
    ostringstream oss;
    oss << val;
    return oss.str ();
  }

  string operator () ( const unsigned short val ) const { return operator () ( static_cast<uint64_t> ( val ) ); }
  string operator () ( const unsigned int   val ) const { return operator () ( static_cast<uint64_t> ( val ) ); }

  //! Format single character
  virtual string operator () ( const char val ) const {
    ostringstream oss;
    oss << val;
    return oss.str ();
  }

  //! Format complex numbers appropriately
  template <typename RealType> string operator () ( const std::complex<RealType> number ) const {
    return string ( "(" ) + operator () ( real ( number ) ) +
           + "," + operator () ( imag ( number ) ) + ")";
  }

  //! Format anything else
  template <typename AnyType> string operator () ( const AnyType thing ) const {
    ostringstream oss;
    oss << thing;
    return oss.str ();
  }
};

//! Format functor that sets some ostream parameters
//! @ingroup Output
class SimpleFormat : public Format {
public:
  //! Sets field width, precision and flags
  SimpleFormat ( const int wid, const int prec, const ios::fmtflags fmt = ios::fixed | ios::right, const char fill = ' ' ) : width ( wid ), precision ( prec ), format ( fmt ), fillchar ( fill ) {}

  string operator () ( const long double val ) const;

  string operator () ( const int64_t val ) const;

  string operator () ( const uint64_t val ) const;

  string operator () ( const char val ) const {
    return ( string ( width, val ) );
  }


  // for other types
  using Format::operator();

private:
  int width;
  int precision;
  ios::fmtflags format;
  char fillchar;
};

//! Format functor for nice mixed formatting
//! @ingroup Output
class MixedFormat : public Format {
public:
  //! Set space before and after decimal point, at least (2, 3) necessary for scientific notation
  MixedFormat ( const int bef, const int aft ) : before ( max ( bef, 2 ) ), after ( max ( aft, 3 ) ) {}

  string operator () ( const long double val ) const;

  string operator () ( const int64_t val ) const;

  string operator () ( const uint64_t val ) const;

  string operator () ( const char val ) const {
    return ( string ( before + 1 + after, val ) );
  }

  // for other types
  using Format::operator();

private:
  int before;
  int after;
};


//! Format functor for hexadecimal formatting
//! @ingroup Output
class HexFormat {
public:
  template <typename AnyType>
  string operator () ( const AnyType & Arg ) const {
    ostringstream oss;
    oss.fill ('0');
    oss << std::hex;

    // for integer argument, we could simply (and equivalently, up to endianness) write
    //   oss << setw ( 2 * sizeof ( AnyType ) ) << std::hex << Arg;
    // but then floating point numbers (or other stuff) would not be printed in hex format
    const unsigned char* ptr = reinterpret_cast<const unsigned char*> ( &Arg );
    ptr += sizeof ( AnyType ) - 1;
    for ( unsigned short i = 0; i < sizeof ( AnyType ); ++i, --ptr ) {
      unsigned short val = *ptr;
      oss << setw ( 2 ); // this needs to be done each time
      oss << val;
    }
    return oss.str ();
  }
};



extern const SimpleFormat scientificFormat;
extern const SimpleFormat longScientificFormat;
extern const SimpleFormat shortFormat;
extern const SimpleFormat intFormat;
extern const SimpleFormat longIntFormat;
extern const SimpleFormat fillFormat;
extern const MixedFormat mixedFormat;
extern const MixedFormat detailedFormat;
extern const HexFormat hexFormat;

//! Binary input from stream, nothing is done about byte order
template <class DataType>
void readbinary ( istream& is, DataType& data ) {
  is.read ( reinterpret_cast<char*> ( &data ), sizeof ( data ) );
}

//! Binary output to stream, nothing is done about byte order
template <class DataType>
void writebinary ( ostream& os, const DataType& data ) {
  os.write ( reinterpret_cast<const char*> ( &data ), sizeof ( data ) );
}

/**
 * \author Berkels
 */
template <typename InDataType, typename OutDataType>
void readBinaryData ( istream &In, OutDataType* Dest, const int Length ) {
  if ( typeid ( InDataType ) == typeid ( OutDataType ) ) {
    In.read ( reinterpret_cast< char* > ( Dest ), Length * sizeof ( InDataType ) );
  } else {
    InDataType *dummy = new InDataType[ Length];
    In.read ( reinterpret_cast< char* > ( dummy ),  Length * sizeof ( InDataType ) );
    for ( int i = 0; i < Length; ++i ) Dest[i] = static_cast< OutDataType > ( dummy[i] );
    delete[] dummy;
  }
}

/**
 * \author Berkels
 */
template <typename InDataType, typename OutDataType>
OutDataType readBinaryData ( istream &In ) {
  OutDataType temp;
  readBinaryData<InDataType, OutDataType> ( In, &temp, 1 );
  return temp;
}

/**
 * \author Berkels
 */
template <typename InDataType, typename OutDataType>
void writeBinaryData ( const InDataType* InputData, const int Length, ostream &Out ) {
  if ( typeid ( OutDataType ) == typeid ( InDataType ) ) {
    Out.write ( reinterpret_cast<const char *> ( InputData ), Length * sizeof ( OutDataType ) );
  } else {
    OutDataType * buffer = new OutDataType[Length];
    for ( int i = 0; i < Length; ++i ) buffer[i] =  static_cast<OutDataType> ( InputData[i] );
    Out.write ( reinterpret_cast<const char *> ( buffer ), Length * sizeof ( OutDataType ) );
    delete[] buffer;
  }
}

/**
 * \author Berkels
 */
template <typename InDataType, typename OutDataType>
void writeBinaryData ( const InDataType &InputData, ostream &Out ) {
  writeBinaryData<InDataType, OutDataType> ( &InputData, 1, Out );
}

//! Converts big endian to little endian and vice versa
//! \todo code overlaps with aol::Vector<DataType>::swapByteOrder()
template <typename DataType>
DataType swapByteOrder ( const DataType Input ) {
  DataType dest;
  const char* op = reinterpret_cast<const char*> ( &Input );
  char* dp = reinterpret_cast<char*> ( &dest );
  int n = sizeof ( DataType );
  for ( int j = 0; j < n; ++j )
    dp [j] = op [n-1-j];
  return dest;
}

//! STL vector output
template <class DataType>
ostream& operator << ( ostream& os, const std::vector<DataType> vec ) {
  for ( int i = 0; i < static_cast<int> ( vec.size () ); ++i )
    os << ( i ? "; " : "" ) << vec [i];
  return os << endl;
}

namespace {
  //! A (relative) difference of APPEQFACTOR * numeric_limits::epsilon () counts as equal in the relative comparison methods
  const int APPEQFACTOR = 16;
  //! A (absolute) difference of APPEQVALUE counts as equal in the absolute comparison methods
  const double APPEQVALUE  = 1E-15;
}

//! Approximate equality for floating point numbers, relative comparison.
//! Attention: Not reasonable if a and b may have different sign.
template <class DataType>
inline bool appeqRelative ( const DataType a, const DataType b ) {
  aol::AssertAtCompileTime < !(std::numeric_limits<DataType>::is_integer) >();
  return aol::Abs ( a - b ) < APPEQFACTOR * std::numeric_limits<DataType>::epsilon () * Max ( abs(a), abs(b) );
}
//! Approximate equality for floating point numbers, absolute comparison
template <class DataType>
inline bool appeqAbsolute ( const DataType a, const DataType b ) {
  aol::AssertAtCompileTime < !(std::numeric_limits<DataType>::is_integer) >();
  return aol::Abs ( a - b ) < APPEQVALUE;
}

//! Approximate greater-equal (even a tiny bit smaller is OK) for floating point numbers, relative comparison
//! Attention: Not reasonable if a and b may be zero or have different sign.
template <class DataType>
inline bool appgeqRelative ( const DataType a, const DataType b ) {
  aol::AssertAtCompileTime < !(std::numeric_limits<DataType>::is_integer) >();
  return a * ( 1 + aol::signum(a) * APPEQFACTOR * std::numeric_limits<DataType>::epsilon () ) > b;
}

//! Approximate lesser-equal (even a tiny bit larger is OK) for floating point numbers, relative comparison
//! Attention: Not reasonable if a and b may be zero or have different sign.
template <class DataType>
inline bool appleqRelative ( const DataType a, const DataType b ) {
  aol::AssertAtCompileTime < !(std::numeric_limits<DataType>::is_integer) >();
  return a * ( 1 - aol::signum(a) * APPEQFACTOR * std::numeric_limits<DataType>::epsilon () ) < b;
}

//! Approximate greater-equal (even a tiny bit smaller is OK) for floating point numbers, absolute comparison
template <class DataType>
inline bool appgeqAbsolute ( const DataType a, const DataType b ) {
  aol::AssertAtCompileTime < !(std::numeric_limits<DataType>::is_integer) >();
  return a + APPEQVALUE > b;
}

//! Approximate lesser-equal (even a tiny bit larger is OK) for floating point numbers, absolute comparison
template <class DataType>
inline bool appleqAbsolute ( const DataType a, const DataType b ) {
  aol::AssertAtCompileTime < !(std::numeric_limits<DataType>::is_integer) >();
  return a - APPEQVALUE < b;
}

//! Returns 1 if x is approximately equal to 0 and 1/x otherwise.
template <typename RealType>
inline RealType Reciprocal ( const RealType X ) {
  return ( aol::appeqAbsolute ( X, aol::ZOTrait<RealType>::zero ) ? aol::ZOTrait<RealType>::one : aol::ZOTrait<RealType>::one / X );
}

//! Compute 32 bit checksum for size bytes of any kind of data
unsigned int crc32 ( const void* const ptr, const unsigned int size );

//! Calculate time between call of @c start and @c stop
/** StopWatch measures cputime.
 *  @ingroup Profiling Utilities
 */
class StopWatch {
public:
  StopWatch() : t1 ( 0 ), t2 ( 0 ), delta ( 0 ), started ( 0 ), finished ( 0 ), _startedMiliseconds( 0 ), _finishedMiliseconds( 0 ) {}

  //! Reset timer
  void start();

  //! Accumulate times
  void cont();

  //! Measure time from start
  void stop();

  //! Returns time in seconds
  //! \warning This method doesn't seem to work correctly if openmp is used!
  double elapsedCpuTime() {
#ifdef _OPENMP
    if ( !_suppressOpenMPWarning ) cerr << "Warning: elapsedCpuTime doesn't seem to work correctly if openmp is used! \n";
#endif
    return ( delta );
  }

  //! Returns elapsed wall clock time in seconds
  double elapsedWallClockTime() {
    return static_cast<double> ( finished - started + (_finishedMiliseconds - _startedMiliseconds)/1000. );
  }

  //! Returns total run time of program
  double total_runtime();

  //! Returns converts the double value Time to a temporary string containing the elapsed time nicely formatted
  //! Assumes the Time to be given in seconds.
  const char* timeToString( double Time ){
    int minutes = static_cast<int> ( Time ) / 60;
    const double seconds = Time - static_cast<double> ( minutes * 60 );
    if ( minutes / 60 > 0 ) {
      const int hours = minutes / 60;
      minutes = minutes % 60;
      sprintf ( _elapsedString, "%dh %02dm %02.3fs", hours, minutes, seconds );
    } else {
      sprintf ( _elapsedString, "%dm %02.3fs", minutes, seconds );
    }
    return ( _elapsedString );
  }

  //! Returns temporary string containing the elapsed time nicely formatted
  const char* elapsedCpuTimeString() {
    return timeToString( elapsedCpuTime() );
  }

  //! Returns temporary string containing the date and time of starting watch
  const char* startedAt();

  //! Returns temporary string containing the date and time of stopping watch
  const char* stoppedAt();

  //! Outputs CPU and wall clock time nicely formatted to the ostream Out.
  void printReport( ostream &Out ) {
    Out << "CPU time        = " << elapsedCpuTimeString() << endl;
    Out << "Wall clock time = " << timeToString( elapsedWallClockTime() ) << endl;
  }

protected:
  double t1, t2, delta;
  time_t started, finished;
  unsigned short _startedMiliseconds, _finishedMiliseconds;
  char _elapsedString[512], _startedString[512], _stoppedString[512];

public:
  static bool _suppressOpenMPWarning;
};

/**
 * Generates a nicely formatted string containing the current time and date.
 *
 * \author Berkels
 */
string generateCurrentTimeAndDateString ( );

/**
 * \brief Given a filename that possibly contains a path this functions
 *        returns the filename without the path.
 *
 * \author Berkels
 */
const char* getFileName( const char* FileNameWithPath );

/**
 * \brief Given a filename that possibly contains a path this functions
 *        returns the filename without the path and without the suffix,
 *        i.e. without the last dot and anything after it.
 *
 * \author Berkels
 */
string getBaseFileName ( const string &FileNameWithPath );

/**
 * \brief Given a filename that contains a path this functions returns
 *        just the path including the trailing directory separator.
 *
 * \author Berkels
 */
string getPath ( const string &FileNameWithPath );

/**
 * \brief Given a string starting with either "~" or with "${ENV_VAR}",
 *        the tilde or the environment variable is replaced with the
 *        home directory or the value of the environment variable with
 *        name ENV_VAR and the altered string is returned. If the input
 *        string does not start with either of the two, the input string
 *        is returned without changes.
 *
 * \author Berkels
 */
string expandTildeOrEnvVar ( const string &InputString );
  
/**
 * \brief Returns true if FileName ends with Ending, false otherwise.
 *        For example useful to check a filename for a certain suffix,
 *        e.g. fileNameEndsWith ( filename, ".bz2" ). The comparison
 *        takes into account trailing whitespace, i.e. a filename
 *        ending with ".bz2 " will still count as ending with ".bz2".
 *        It also ignores the case of the characters, i.e. a filename
 *        ending with ".PNG" will still count as ending with ".png".∫
 *
 * \author Berkels
 */
bool fileNameEndsWith ( const char* FileName, const char* Ending );

/**
 * \brief Returns true if FileName ends with a suffix that indicates bzip compression, i.e. ".bz2".
 *
 * \author Berkels
 */
bool hasBzipSuffix ( const char *FileName );

/**
 * \brief Returns true if FileName ends with a suffix that indicates zlib compression, i.e. ".zz".
 *
 * \author Geihe
 */
bool hasZlibSuffix ( const char *FileName );

/**
 * \brief Returns true if FileName ends with a suffix that indicates gzip compression, i.e. ".gz".
 *
 * \author Geihe
 */
bool hasGzipSuffix ( const char *FileName );

//! to_string(s1, 10.5); // double to string
//! to_string(s2, 123);  // int to string
//! to_string(s3, true); // bool to string
template <class T>
void to_string ( string & result, const T & t ) {
  ostringstream oss; // create a stream
  oss << t;          // insert value to stream
  result = oss.str();  // extract value and write it to result
}

template <class T>
string to_string ( const T & t ) {
  string result;
  ostringstream oss; // create a stream
  oss << t;          // insert value to stream
  result = oss.str();  // extract value and write it to result
  return result;
}

//! Use like this:
//! \code
//! double d;
//! string salary;
//! string s="12.56";
//! d=convert<double> (s); //d equals 12.56
//! salary=convert<string> (9000.0);//salary equals "9000"
//! \endcode
template <class out_type, class in_value>
out_type convert ( const in_value & t ) {
  stringstream stream;
  stream << t;      // insert value to stream
  out_type result;  // store conversion result here
  stream >> result; // write value to result
  return result;
}

//! Generate next of a sequence of output filenames,
//! adding numbers between the last (or if !fromEnd first) two dots
//! Sequence is like file.ext, file.0000.ext, file.0001.ext, ...
string getOutFileName ( const string &inFileName, bool fromEnd = true );

//! Give back formatted string, analogously to
//! sprintf, but save the long way 'round with
//! char arrays.
string strprintf(const char * format, ...);

/**
 * Reads the next line or string (controlled by the Argument "Line") from the
 * istream In and compares it with ExpectedString. If the comparison fails,
 * an expection is thrown.
 *
 * \author Berkels
 */
void checkNextLineOrString ( istream &In, const char *ExpectedString, const bool Line = true );

/**
 * \author Berkels
 */
int countDigitsOfNumber ( const int N );

template <int length, typename T>
class auto_container {
private:
  struct data_entity {
    data_entity() : _ptr ( NULL ), _isReference ( false ) { }
    data_entity ( T* ptr, bool isReference ) : _ptr ( ptr ), _isReference ( isReference ) { }
    ~data_entity() {
      if ( _ptr && !_isReference ) delete _ptr;
    }
    T* _ptr;
    bool _isReference;
  };

  data_entity data[length];

public:
  auto_container() {}

  ~auto_container() {}

  //! adds a copy to the container, delete the copy on overwrite/destruction
  void set_copy ( int i, const T& obj ) {
    if ( i < 0 || i >= length ) throw aol::Exception ( "hey, your index is out of bounds!", __FILE__, __LINE__ );
    if ( data[i]._ptr != NULL && !(data[i]._isReference) ) delete data[i]._ptr;
    data[i]._ptr = new T ( obj );
    data[i]._isReference = false;
  }

  //! adds the given object to the container, do not take ownership of it
  void set_reference ( int i, const T& obj ) {
    if ( i < 0 || i >= length ) throw aol::Exception ( "hey, your index is out of bounds!", __FILE__, __LINE__ );
    if ( data[i]._ptr != NULL && !(data[i]._isReference) ) delete data[i]._ptr;
    data[i]._ptr =  &obj;
    data[i]._isReference = true;
  }

  //! adds the given object to the container, take ownership
  //! this methods gets a pointer to make obvious that the object will be
  //! destructed via delete().
  void take_over ( int i, T * obj ) {
    if ( i < 0 || i >= length ) throw aol::Exception ( "hey, your index is out of bounds!", __FILE__, __LINE__ );
    if ( data[i]._ptr != NULL && !(data[i]._isReference) ) delete data[i]._ptr;
    data[i]._ptr =  obj;
    data[i]._isReference = false;
  }

  T &operator[] ( int i ) {
    if ( i<0 || i >= length ) throw aol::Exception ( "hey, your index is out of bounds!", __FILE__, __LINE__ );
    return * ( data[i]._ptr );
  }

  const T &operator[] ( int i ) const {
    if ( i<0 || i >= length ) throw aol::Exception ( "hey, your index is out of bounds!", __FILE__, __LINE__ );
    return * ( data[i]._ptr );
  }

};

//! for consistent output of successful selfTests (message should be 80 characters wide)
void printSelfTestSuccessMessage ( const string &message );

//! for consistent output of failed selfTests (message should be 80 characters wide)
void printSelfTestFailureMessage ( const string &message );

//! Check command line arguments for the tools/benchmark syntax.
//! \author Berkels
bool checkForBenchmarkArguments ( const int argn, char **argv, string &ResultFilename );

//! Consistenly log a benchmark rating to the console and a file.
//! \author Berkels
template <typename RealType>
void logBenchmarkResult ( const string &BenchmarkName, const RealType Nupsi, const RealType Wupsi, const string &ResultFilename ) {
  cerr << "Congratulations: your computer has achieved " << Nupsi << " nupsi and " << Wupsi << " wupsi!" << endl;
  if ( ResultFilename.size() > 0 ) {
    std::ofstream resultfile ( ResultFilename.c_str (), std::ios_base::app );
    resultfile << BenchmarkName;
    for ( int i = BenchmarkName.size(); i < 58; ++i )
      resultfile << " ";
    resultfile << "| " << aol::longIntFormat (static_cast<int> (Nupsi)) << " | " << aol::longIntFormat (static_cast<int> (Wupsi)) << endl;
  }
}

template < typename anything >
void doNothingWithArgumentToPreventUnusedParameterWarning ( const anything & ) {
}

#if defined ( USE_CPP11 )
template<typename... ParamTypes>
void doNothingWithArgumentsToPreventUnusedParameterWarning ( const ParamTypes&... ) {
}
#endif

void READ_COMMENTS ( istream &in );

/** magic characters in saved aol::Vector, MultiVector, qc::ScalarArray, indicating which of these data structures is saved in a file
 *  \todo implement BitArray as LO, LP, LQ
 *  \todo use in qc::ScalarArray
 *  \author Schwen (MEVIS)
 */
struct VectorFileMagicChar {
  static const char
    BitVector = 'L',
    Vector = 'M',
    MultiVector = 'N',
    ScalarArray_1D = 'O',
    ScalarArray_2D = 'P',
    ScalarArray_3D = 'Q',
    MultiArray = 'U';
};

/** magic numbers in saved aol data structures, indicating which data type is saved in binary format.
 *  The order is somewhat mixed for historical reasons (pgm magic numbers were first extended for some data types, later for more)
 *  \author Schwen (MEVIS)
 */
enum FileFormatType {
  FF_BOOL            = 38,
  FF_SIGNED_CHAR     = 32,
  FF_UNSIGNED_CHAR   = 5,
  FF_SIGNED_SHORT    = 11,
  FF_UNSIGNED_SHORT  = 10,
  FF_SIGNED_INT      = 33,
  FF_UNSIGNED_INT    = 34,
  FF_SIGNED_64INT    = 35,
  FF_UNSIGNED_64INT  = 36,
  FF_FLOAT           = 8,
  FF_DOUBLE          = 9,
  FF_LONG_DOUBLE     = 37
};

/** magic numbers for data types in file formats
 *  \author Schwen (MEVIS)
 */
template < typename DataType >
class FileFormatMagicNumber {
public:
  static const aol::FileFormatType FFType;
  static const string FFContainedName;
};

//! \brief Copies a file.
//! \param[in] source Source file name.
//! \param[in] dest Destination file name.
//! \author Toelkes
void copyFile ( const char* source, const char* dest );

void printHGChangesetID ( ostream& os );

} // namespace aol

// Use defined output operator without qualification
using aol::operator<<;

// the following will be used to detect whether a module header uses namespaces it shouldn't
namespace aol {
  typedef int   avoid_namespace_collision;
}
namespace qc {
  typedef float avoid_namespace_collision;
}

#endif

