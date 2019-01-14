#ifndef __SMALLVEC_H
#define __SMALLVEC_H

#include <aol.h>

namespace aol {

//! \brief A Vector of any dimension (template-parameter),
//!        class is relatively simple, for beeing able to
//!        build a vector of vectors  (ON)
//! @ingroup Vector, Point, Matrix
template <int dimension, typename _DataType> class Vec {
protected:
  // Allow dummy vectors with dimension == 0, but avoid arrays of size 0.
  _DataType coords[dimension > 0 ? dimension : 1];

  static const Format* _pFormat;
  static bool prettyFormat;

public:
  typedef _DataType DataType;
  static const int Depth = 1;

  //! RealType is used for all operations that require floating point numbers, like calculating the norm of a vector of integers.
  //!TODO remove compile errors when the following line is uncommented and uncomment it afterwards!
  //typedef typename RealTrait<DataType>::RealType RealType;


  //! Constructor for creating Vec from c array
  explicit Vec ( const _DataType rhs[dimension] ) {
    for ( int i = 0; i < dimension; ++i )
      coords[i] = rhs[i];
  }

  //! Constructor for Vec with constant entries
  explicit Vec ( const _DataType Scalar ) {
    for ( int i = 0; i < dimension; ++i )
      coords[i] = Scalar;
  }

  //! Standard constructor
  Vec( ) {
    for ( int i = 0; i < dimension; ++i )
      coords[i] = ZOTrait<_DataType>::zero;
  }

  //! Copy-constructor
  Vec ( const Vec<dimension, _DataType> &rhs ) {
    for ( int i = 0; i < dimension; ++i )
      coords[i] = rhs.coords[i];
  }

  //! Copy-constructor for structure or data copying
  Vec ( const Vec<dimension, _DataType> &rhs, CopyFlag copyFlag ) {
    switch ( copyFlag ) {
      case DEEP_COPY:
        for ( int i = 0; i < dimension; ++i )
          coords[i] = rhs.coords[i];
        break;
      case STRUCT_COPY:
        for ( int i = 0; i < dimension; ++i )
          coords[i] = ZOTrait<_DataType>::zero;
        break;
      default:
        string errorMessage = strprintf( "Copying a Vec is not possible with copyFlag=%d", copyFlag );
        throw Exception( errorMessage, __FILE__, __LINE__);
    }
  }

  //! Conversion-Copy-constructor
  template <class AnotherType>
  explicit Vec ( const Vec<dimension, AnotherType> &rhs ) {
    for ( int i = 0; i < dimension; ++i )
      coords[i] = static_cast<_DataType> ( rhs[i] );
  }

  //! operator=
  Vec<dimension, _DataType>& operator= ( const Vec<dimension, _DataType> &rhs ) {
    for ( int i = 0; i < dimension; ++i )
      coords[i] = rhs.coords[i];
    return *this;
  }

  enum { dim = dimension };

private:
  //! operator=
  // private, to avoid implementation, danger of errors (UD)
  Vec<dimension, _DataType>& operator= ( _DataType /* scalar */ ) {
    throw aol::Exception ( "aol::Vec<dimension, _DataType>::operator= ( _DataType scalar ) should not be used", __FILE__, __LINE__ );
  }

public:
  //! Set output format
  static void setFormat ( const Format& NewFormat ) {
    _pFormat = &NewFormat;
  }

  //! Get output format
  static const Format& getFormat ( ) {
    return *(_pFormat);
  }

  //! Set pretty output formatting on/off
  static void setPrettyFormat ( bool newPrettyFormat ) {
    prettyFormat = newPrettyFormat;
  }

  //! Get pretty output formatting status
  static bool getPrettyFormat ( ) {
    return prettyFormat;
  }

  //! Set Vec to Co.
  void set ( const Vec &Co ) {
    setMultiple ( Co, ZOTrait<DataType>::one );
  }

  //! Set Vec to Fator * Co.
  void setMultiple ( const Vec &Co, DataType Factor ) {
    for ( int i = 0; i < dimension; ++i )
      coords[i] = Co.coords[i] * Factor;
  }

  //! Equality operator.
  bool operator== ( const Vec &Co ) const {
    for ( int i = 0; i < dimension; ++i ) {
      if ( Co.coords[i] != coords[i] )
        return false;
    }
    return true;
  }

  //! Inequality operator.  DOES THIS WORK?
  bool operator!= ( const Vec &Co ) const {
    return !operator== ( Co );
  }

  //! Returns Ith component.
  _DataType get ( int I ) const {
#ifdef BOUNDS_CHECK
    if ( I < 0 || I >= dimension ) {
      char error[1024];
      sprintf ( error, "aol::VecN<_DataType>::get: index %d out of bounds. N = %d\n", I, dimension );
      throw OutOfBoundsException ( error, __FILE__, __LINE__ );
    }
#endif
    return coords[I];
  }

  //! Sets all entries to Val.
  void setAll ( _DataType Val ) {
    for ( int i = 0; i < dimension; ++i )
      coords [i] = Val;
  }

  //! Set the Vec to zero.
  //!
  void setZero( ) {
    for ( int i = 0; i < dimension; ++i )
      coords [i] = ZOTrait<_DataType>::zero;
  }

  // for compatibility with aol::Vector
  // changing size does not make sense, so just set to zero
  void reallocate( const Vec & /*other*/ ) {
    setZero();
  }

  template <typename InDataType>
  void readFromBuffer ( InDataType* buffer ) {
    for ( int i = 0; i < dimension; ++i )
      coords [i] = static_cast<InDataType> ( buffer[i] );
  }

  inline bool checkForNANsAndINFs() const {
    for ( int i = 0; i < dimension; ++i ) {
      if ( !aol::isFinite ( coords[i] ) ) return true;
    }
    return false;
  }

  //! Unary minus
  Vec <dimension, _DataType> operator- () const {
    Vec<dimension, _DataType> res;
    for ( int i = 0; i < dimension; ++i )
      res [i] = - coords [i];
    return res;
  }

  //! Bracket operator.
  inline _DataType& operator[] ( int I ) {
#ifdef BOUNDS_CHECK
    if ( I < 0 || I >= dimension ) {
      char error[1024];
      sprintf ( error, "aol::VecN<_DataType>::operator[]: index %d out of bounds. N = %d\n", I, dimension );
      throw OutOfBoundsException ( error, __FILE__, __LINE__ );
    }
#endif
    return coords[I];
  }

  //! square bracket operator
  inline const _DataType& operator[] ( int I ) const {
#ifdef BOUNDS_CHECK
    if ( I < 0 || I >= dimension ) {
      char error[1024];
      sprintf ( error, "aol::VecN<_DataType>::operator[]: index %d out of bounds. N = %d\n", I, dimension );
      throw OutOfBoundsException ( error, __FILE__, __LINE__ );
    }
#endif
    return coords[I];
  }

  //! Scalar product operator
  _DataType operator* ( const Vec<dimension, _DataType> &Other ) const {
    _DataType scalarProduct = 0;
    for ( int i = 0; i < dimension; ++i )
      scalarProduct += coords[i] * Other.coords[i];
    return scalarProduct;
  }

  //! Dot product
  _DataType dotProduct ( const Vec<dimension, _DataType> &Other ) const {
    return this->operator*( Other );
  }

  //! multiplication by scalar
  Vec<dimension, _DataType> operator * ( _DataType alpha ) const {
    Vec<dimension, _DataType> res;
    for ( int i = 0; i < dimension; ++i )
      res [i] = alpha * coords [i];
    return res;
  }

  //! division by scalar
  Vec<dimension, _DataType> operator / ( _DataType alpha ) const {
    Vec<dimension, _DataType> res;
    for ( int i = 0; i < dimension; ++i )
      res [i] = coords [i] / alpha;
    return res;
  }

  //! add two vectors and return the result
  Vec<dimension, _DataType> operator + ( const Vec<dimension, _DataType> &Other ) const {
    Vec<dimension, _DataType> res;
    for ( int i = 0; i < dimension; ++i )
      res [i] = coords [i] + Other[i];
    return res;
  }

  //! subtract two vectors and return the result
  Vec<dimension, _DataType> operator - ( const Vec<dimension, _DataType> &Other ) const {
    Vec<dimension, _DataType> res;
    for ( int i = 0; i < dimension; ++i )
      res [i] = coords [i] - Other[i];
    return res;
  }

  //! add other vector
  Vec<dimension, _DataType>& operator+= ( const Vec<dimension, _DataType> &Other ) {
    for ( int i = 0; i < dimension; ++i )
      coords[i] += Other.coords[i];
    return *this;
  }

  //! subtract other vector
  Vec<dimension, _DataType>& operator-= ( const Vec<dimension, _DataType> &Other ) {
    for ( int i = 0; i < dimension; ++i )
      coords[i] -= Other.coords[i];
    return *this;
  }

  //! multiply by scalar
  Vec<dimension, _DataType>& operator*= ( _DataType Alpha ) {
    for ( int i = 0; i < dimension; ++i )
      coords[i] *= Alpha;
    return *this;
  }

  //! divide by scalar
  Vec<dimension, _DataType>& operator/= ( _DataType Alpha ) {
    for ( int i = 0; i < dimension; ++i )
      coords[i] /= Alpha;
    return *this;
  }

  //! add multiple of other vector
  void addMultiple ( const Vec<dimension, _DataType>& AddedVec, _DataType Factor ) {
    for ( int i = 0; i < dimension; ++i )
      coords[i] += Factor * AddedVec.coords[i];
  }

  //! add a scalar to all components
  void addToAll ( _DataType Scalar ) {
    for ( int i = 0; i < dimension; ++i )
      coords[i] += Scalar;
  }

  _DataType sum( ) const {
    _DataType res = 0;
    for ( int i = 0; i < dimension; ++i )
      res += coords[i];
    return res;
  }

  _DataType prod( ) const {
    _DataType res = 1;
    for ( int i = 0; i < dimension; ++i )
      res *= coords[i];
    return res;
  }

  _DataType normSqr( ) const {
    _DataType res = 0;
    for ( int i = 0; i < dimension; ++i )
      res += ( coords[i] * coords[i] );
    return res;
  }

  typename RealTrait<_DataType>::RealType norm( ) const {
    return sqrt ( static_cast<typename RealTrait<_DataType>::RealType>(normSqr()) );
  }

  _DataType norm ( _DataType p ) const {
    _DataType res = 0;
    for ( int i = 0; i < dimension; ++i )
      res += pow ( Abs ( coords[i] ), p );
    return pow ( res, 1.0 / p );
  }

  _DataType infinityNorm( ) const {
    _DataType res = 0;
    for ( int i = 0; i < dimension; ++i )
      if ( Abs ( coords[i] ) > res )
        res = Abs ( coords[i] );
    return res;
  }

  void normalize( ) {
    typedef typename RealTrait<_DataType>::RealType T;

    T ns = sqrt ( static_cast<T> ( normSqr() ) );
    if ( !ns ) throw Exception ( "Vec<dimension,_DataType>::normalize: Length of vector is zero", __FILE__, __LINE__ );

    for ( int i = 0; i < dimension; ++i )
      coords[i] = static_cast<_DataType> ( static_cast<T> ( coords[i] ) / ns );
  }

  //! Returns the minimal value in this vector.
  _DataType getMinValue() const {
    _DataType min = coords[0];
    for ( int i = 1; i < dimension; ++i ) {
      if ( coords[ i ] < min ) min = coords[ i ];
    }
    return min;
  }

  //! Returns the maximum value in this vector.
  _DataType getMaxValue() const {
    _DataType max = coords[0];
    for ( int i = 1; i < dimension; ++i ) {
      if ( coords[ i ] > max ) max = coords[ i ];
    }
    return max;
  }

  //! return std::pair of index and value of maximal entry
  std::pair< int, _DataType > getMaxIndexAndValue() const {
    std::pair< int, _DataType > IndVal ( 0,  coords[ 0 ] );
    for ( int i = 1; i < dimension; ++i ) {
      if ( coords[ i ] > IndVal.second ) {
        IndVal.first = i;
        IndVal.second = coords[ i ];
      }
    }
    return ( IndVal );
  }

  //! Returns the minimum absolute value in this vector (not a norm)
  _DataType getMinAbsValue() const {
    _DataType min = Abs ( coords[ 0 ] );
    for ( int i = 1; i < dimension; ++i ) {
      if ( Abs ( coords[ i ] ) < min ) min = Abs ( coords[ i ] );
    }
    return min;
  }

  //! Returns the arithmetic mean value of the values of the vector.
  typename RealTrait<_DataType>::RealType getMeanValue() const {
    typedef typename RealTrait<_DataType>::RealType RealType;
    return ( static_cast<RealType> ( sum ( ) ) / static_cast<RealType> ( dimension ) );
  }

  //! Returns a standard deviation of the values: sqrt ( ( (v[0] - v)^2 + (v[1] - v)^2 + ... + (v[n-1] - v)^2 ) / ( n - 1 ) ) [Sic: n-1]
  typename RealTrait<_DataType>::RealType getStdDev() const {
    typedef typename RealTrait<_DataType>::RealType RealType;
    if ( dimension == 0 ) {
      return ( aol::NumberTrait<RealType>::NaN );
    } else if ( dimension == 1 ) {
      return ( NumberTrait<RealType>::zero );
    } else {
      const RealType mean = getMeanValue();
      RealType sumsqr = 0.0;
      for ( int i = 0; i < dimension; ++i ) {
        sumsqr += aol::Sqr ( coords[ i ] - mean );
      }
      return ( sqrt ( sumsqr / ( dimension - 1 ) ) );
    }
  }

  int numNonZeroes( ) const {
    int nNonZeroes = 0;
    for ( int i = 0; i < dimension; ++i ) {
      if ( coords[i] != 0 )
        ++nNonZeroes;
    }
    return nNonZeroes;
  }

  int getTotalSize() const {
    return dimension;
  }

  ostream& print ( ostream &os ) const {
    os << ( prettyFormat ? "( " : "" );
    for ( int i = 0; i < dimension - 1; ++i )
      os << ( *_pFormat ) ( coords[i] ) << ( prettyFormat ? " ; " : " " );
    os << ( *_pFormat ) ( coords[ dimension - 1 ] )
       << ( prettyFormat ? " )" : "" );
    return os;
  }

  //! Saves values into ASCII-file.
  void saveASCII ( const char *fileName ) const {
    ofstream file ( fileName );
    for ( int i = 0; i < dimension; ++i )
      file << coords[i] << endl;
    file.close();
  }

  istream& read ( istream &is ) {
    char c;
    is >> c;
    if ( c != '(' ) is.putback ( c );
    is >> coords [0];
    for ( int i = 1; i < dimension; ++i ) {
      is >> c;
      if ( c != ';' ) is.putback ( c );
      is >> coords [i];
    }
    is >> c;
    if ( c != ')' ) is.putback ( c );
    return is;
  }

  //! inverse(!)-lexicographic comparison
  bool operator < ( const Vec<dimension, _DataType>& Co ) const {
    for ( int i = dimension-1; i >= 0; --i )
      if ( coords[i] != Co.coords[i] ) return coords[i] < Co.coords[i];
    return false; // this == Co
  }

  //! set value at position I ( needed by aol::DerivativeValidatorBase<> )
  Vec<1, int> setIthComponent ( const int I, const DataType Value ) {
    coords[I] = Value;
    Vec<1, int> res;
    res[0] = I;
    return res;
  }

  const _DataType * getData () const {  return coords;  }
        _DataType * getData ()       {  return coords;  }
};

//! left-multiplication scalar * Vec
template<int dimension, typename _DataType>
Vec<dimension, _DataType> operator* ( const _DataType alpha, const Vec<dimension, _DataType> &vec ) {
  Vec<dimension, _DataType> res;
  for ( int i = 0; i < dimension; ++i )
    res[i] = alpha * vec[i];

  return ( res );
}



//! Euclidian distance of two points, i. e. norm of the difference of two Vecs
template<int dimension, typename _DataType >
_DataType euclidianDist ( const Vec<dimension, _DataType> &vec1,
                         const Vec<dimension, _DataType> &vec2 ) {
  _DataType distSqr = 0.;
  for ( int i = 0; i < dimension; ++i ) distSqr += Sqr ( vec1[i] - vec2[i] );
  return sqrt ( distSqr );
}



//////////////////////////////////////////////////////////////
// THE CLASS VEC2 DERIVED FROM VEC
// / /////////////////////////////////////////////////////////

//! A simple 2D-Vector.
//! Including all methods which aren't in class Vec
//! @ingroup Vector, Point
template <typename _DataType>
class Vec2 : public Vec<2, _DataType> {

public:

  //! Constructor for two values
  explicit Vec2 ( const _DataType X, const _DataType Y ) : Vec<2, _DataType>() {
    this->coords[0] = X;
    this->coords[1] = Y;
  }

  //! Constructor for c-array
  explicit Vec2 ( const _DataType rhs[2] ) : Vec<2, _DataType> ( rhs ) { }

  //! Constructor for constant vector
  explicit Vec2 ( const _DataType Scalar ) : Vec<2, _DataType> ( Scalar ) { }

  //! Standard constructor
  Vec2() : Vec<2, _DataType>() { }

  //! Copy-constructor
  Vec2 ( const Vec2<_DataType> &rhs ) : Vec<2, _DataType> ( rhs ) {}

  //! Copy-constructor
  Vec2 ( const Vec2<_DataType> &rhs, CopyFlag copyFlag ) : Vec<2, _DataType> ( rhs, copyFlag ) { }

  //! \brief Conversion-constructor
  //! \todo Make this explicit. Currently, it is used to implicitly convert Vec<2, DataType> references
  //!       to Vec2<DataType> references involving a temporary copy of the vector.
  Vec2 ( const Vec<2, _DataType> &rhs ) : Vec<2, _DataType> ( rhs ) {}

  //! Conversion-Copy-constructor
  template <class AnotherType>
  explicit Vec2 ( const Vec<2, AnotherType> &rhs ) : Vec<2, _DataType> ( rhs ) { }


  //! Returns x component.
  inline _DataType x() const {
    return this->coords[0];
  }

  //! Returns reference to x component
  inline _DataType &xref() {
    return this->coords[0];
  }

  //! Returns const reference to x component
  inline const _DataType &xrefc() const {
    return this->coords[0];
  }

  //! Returns y component.
  inline _DataType y() const {
    return this->coords[1];
  }

  //! Returns reference to y component
  inline _DataType &yref() {
    return this->coords[1];
  }

  //! Returns const reference to y component
  inline const _DataType &yrefc() const {
    return this->coords[1];
  }

  //! Set Vec2 to Co.
  inline void set ( const Vec2 &Co ) {
    this->coords[0] = Co.coords[0];
    this->coords[1] = Co.coords[1];
  }

  //! Set Vec2 to ( X, Y )
  inline void set ( _DataType X, _DataType Y ) {
    this->coords[0] = X;
    this->coords[1] = Y;
  }

  //! multiplication by scalar (from the right)
  Vec2<_DataType> operator* ( _DataType alpha ) const {
    return ( Vec2<_DataType> ( alpha * this->coords[0], alpha * this->coords[1] ) );
  }

  //! division by scalar
  Vec2<_DataType> operator/ ( _DataType alpha ) const {
    return ( Vec2<_DataType> ( this->coords[0] / alpha, this->coords[1] / alpha ) );
  }

  //! addition of Vec2s. Attention: slow because temporary object must be created.
  Vec2<_DataType> operator+ ( const Vec2<_DataType> &other ) const {
    return ( Vec2<_DataType> ( this->coords[0] + other.coords[0], this->coords[1] + other.coords[1] ) );
  }

  //! subtraction of Vec2s. Attention: slow because temporary object must be created.
  Vec2<_DataType> operator- ( const Vec2<_DataType> &other ) const {
    return ( Vec2<_DataType> ( this->coords[0] - other.coords[0], this->coords[1] - other.coords[1] ) );
  }


  //! add other vector
  Vec2<_DataType>& operator+= ( const Vec2<_DataType> &Other ) {
    this->coords[0] += Other.coords[0];
    this->coords[1] += Other.coords[1];
    return ( *this );
  }

  //! subtract other vector
  Vec2<_DataType>& operator-= ( const Vec2<_DataType> &Other ) {
    this->coords[0] -= Other.coords[0];
    this->coords[1] -= Other.coords[1];
    return *this;
  }

  //! multiply by scalar
  Vec2<_DataType>& operator*= ( _DataType Alpha ) {
    this->coords[0] *= Alpha;
    this->coords[1] *= Alpha;
    return *this;
  }

  //! divide by scalar
  Vec2<_DataType>& operator/= ( _DataType Alpha ) {
    this->coords[0] /= Alpha;
    this->coords[1] /= Alpha;
    return *this;
  }

  //! Scalar product of this Vec2 with another Vec2
  _DataType operator* ( const Vec2<_DataType> &Other ) const {
    return ( this->coords[0] * Other.coords[0] + this->coords[1] * Other.coords[1] );
  }

  //! unary minus
  Vec2<_DataType> operator- () const {
    return ( Vec2<_DataType>( - this->coords[0], - this->coords[1]) );
  }

  //! Rotate clockwise by pi/2
  inline void rotateRight () {
    _DataType temp = this->coords [0];
    this->coords [0] = this->coords [1];
    this->coords [1] = - temp;
  }

  //! Rotate counter-clockwise pi/2
  inline void rotateLeft () {
    _DataType temp = this->coords [1];
    this->coords [1] = this->coords [0];
    this->coords [0] = - temp;
  }

  istream& read ( istream &is ) {
    char c;
    is >> c;
    if ( c != '(' ) is.putback ( c );
    is >> xref ();
    is >> c;
    if ( c != ';' ) is.putback ( c );
    is >> yref ();
    is >> c;
    if ( c != ')' ) is.putback ( c );
    return is;
  }

  //! inverse-lexicographic comparison
  bool operator < ( const Vec2& Co ) const {
    if ( y() != Co.y() ) return y() < Co.y();
    if ( x() != Co.x() ) return x() < Co.x();
    return false; // this == Co
  }

};

//! left-multiplication scalar * Vec2
template<typename _DataType>
Vec2<_DataType> operator* ( const _DataType alpha, const Vec2<_DataType> &vec ) {
  return ( Vec2<_DataType> ( alpha * vec[0], alpha * vec[1] ) );
}


//////////////////////////////////////////////////////////////
// THE CLASS VEC3 DERIVED FROM VEC
// / /////////////////////////////////////////////////////////

//! A simple 3D-Vector.
//! Including all methods which aren't in class Vec
//! @ingroup Vector, Point
template <typename _DataType>
class Vec3 : public Vec<3, _DataType> {


public:
  //! Constructor for three values
  explicit Vec3 ( const _DataType X, const _DataType Y, const _DataType Z ) {
    this->coords[0] = X;
    this->coords[1] = Y;
    this->coords[2] = Z;
  }

  //! Constructor for c array
  explicit Vec3 ( const _DataType rhs[3] ) : Vec<3, _DataType> ( rhs ) { }


  //! Constructor for two values, setting the last value to zero
  explicit Vec3 ( const _DataType X, const _DataType Y ) : Vec<3, _DataType>() {
    this->coords[0] = X;
    this->coords[1] = Y;
    this->coords[2] = ZTrait<_DataType>::zero ;
  }

  //! Constructor for constant vector
  explicit Vec3 ( const _DataType Scalar ) : Vec<3, _DataType> ( Scalar ) { }

  //! Standard constructor
  explicit Vec3() : Vec<3, _DataType>() {  }

  //! Copy-constructor
  Vec3 ( const Vec3<_DataType> &rhs ) : Vec<3, _DataType> ( rhs ) { }

  //! Copy-constructor
  Vec3 ( const Vec3<_DataType> &rhs, CopyFlag copyFlag ) : Vec<3, _DataType> ( rhs, copyFlag ) { }

  //! Conversion-Copy-constructor
  template <class AnotherType>
  explicit Vec3 ( const Vec<3, AnotherType> &rhs ) : Vec<3, _DataType> ( rhs ) { }

  //! Copy-constructor (for om::TriMesh<>::Point)
  template <typename Point>
  explicit Vec3 ( const Point &p ) {
    this->coords[0] = p[0];
    this->coords[1] = p[1];
    this->coords[2] = p[2];
  }

  //! operator= (for om::TriMesh<>::Point )
  template <typename Point>
  Vec3<_DataType>& operator= ( const Point &rhs ) {
    for ( int i = 0; i < 3; ++i )
      this->coords[i] = rhs[i];
    return *this;
  }

  //! Returns x component.
  inline _DataType x() const {
    return this->coords[0];
  }

  //! Returns reference to x component
  inline _DataType &xref() {
    return this->coords[0];
  }

  //! Returns const reference to x component
  inline const _DataType &xrefc() const {
    return this->coords[0];
  }

  //! Returns y component.
  inline _DataType y() const {
    return this->coords[1];
  }

  //! Returns reference to y component
  inline _DataType &yref() {
    return this->coords[1];
  }

  //! Returns const reference to y component
  inline const _DataType &yrefc() const {
    return this->coords[1];
  }

  //! Returns z component.
  inline _DataType z() const {
    return this->coords[2];
  }

  //! Returns reference to z component
  inline _DataType &zref() {
    return this->coords[2];
  }

  //! Returns const reference to y component
  inline const _DataType &zrefc() const {
    return this->coords[2];
  }

  //! Set Vec3 to Co.
  void set ( const Vec3 &Co ) {
    this->coords[0] = Co.coords[0];
    this->coords[1] = Co.coords[1];
    this->coords[2] = Co.coords[2];
  }

  //! Set Vec3 to (X, Y, Z).
  void set ( _DataType X, _DataType Y, _DataType Z ) {
    this->coords[0] = X;
    this->coords[1] = Y;
    this->coords[2] = Z;
  }

  //! multiplication by scalar
  Vec3<_DataType> operator* ( _DataType alpha ) const {
    return ( Vec3<_DataType> ( alpha * this->coords[0], alpha * this->coords[1], alpha * this->coords[2] ) );
  }

  //! division by scalar
  Vec3<_DataType> operator/ ( _DataType alpha ) const {
    return ( Vec3<_DataType> ( this->coords[0] / alpha, this->coords[1] / alpha, this->coords[2] / alpha ) );
  }

  //! addition of Vec3s. Attention: slow because temporary object must be created.
  Vec3<_DataType> operator+ ( const Vec3<_DataType> &other ) const {
    return ( Vec3<_DataType> ( this->coords[0] + other.coords[0], this->coords[1] + other.coords[1], this->coords[2] + other.coords[2] ) );
  }

  //! subtraction of Vec3s. Attention: slow because temporary object must be created.
  Vec3<_DataType> operator- ( const Vec3<_DataType> &other ) const {
    return ( Vec3<_DataType> ( this->coords[0] - other.coords[0], this->coords[1] - other.coords[1], this->coords[2] - other.coords[2] ) );
  }

  //! add other vector
  Vec3<_DataType>& operator+= ( const Vec3<_DataType> &Other ) {
    this->coords[0] += Other.coords[0];
    this->coords[1] += Other.coords[1];
    this->coords[2] += Other.coords[2];
    return ( *this );
  }

  //! subtract other vector
  Vec3<_DataType>& operator-= ( const Vec3<_DataType> &Other ) {
    this->coords[0] -= Other.coords[0];
    this->coords[1] -= Other.coords[1];
    this->coords[2] -= Other.coords[2];
    return *this;
  }

  //! multiply by scalar
  Vec3<_DataType>& operator*= ( _DataType Alpha ) {
    this->coords[0] *= Alpha;
    this->coords[1] *= Alpha;
    this->coords[2] *= Alpha;
    return *this;
  }

  //! divide by scalar
  Vec3<_DataType>& operator/= ( _DataType Alpha ) {
    this->coords[0] /= Alpha;
    this->coords[1] /= Alpha;
    this->coords[2] /= Alpha;
    return *this;
  }

  //! add multiple of other vector
  template< typename Point >
  void addMultiple ( const Point& AddedVec, _DataType Factor ) {
    for ( int i = 0; i < 3; ++i )
      this->coords[i] += Factor * AddedVec[i];
  }

  //! Scalar product of this Vec3 with another Vec3
  _DataType operator* ( const Vec3<_DataType> &Other ) const {
    return ( this->coords[0] * Other.coords[0] + this->coords[1] * Other.coords[1] + this->coords[2] * Other.coords[2] );
  }

  //! inverse-lexicographic comparison
  bool operator < ( const Vec3& Co ) const {
    if ( z() != Co.z() ) return z() < Co.z();
    if ( y() != Co.y() ) return y() < Co.y();
    if ( x() != Co.x() ) return x() < Co.x();
    return false; // this == Co
  }

  //! set this Vec3 to the cross product (vector product) of two other vectors
  void makeCrossProduct ( const Vec3<_DataType> &a, const Vec3<_DataType> &b ) {
    this->coords[0] = a[1] * b[2] - a[2] * b[1];
    this->coords[1] = a[2] * b[0] - a[0] * b[2];
    this->coords[2] = a[0] * b[1] - a[1] * b[0];
  }

  //! compute cross product (vector product) of this Vec3 with another Vec3
  Vec3<_DataType> crossProduct ( const Vec3<_DataType> &b ) const {
    Vec3<_DataType> res;
    res[0] = this->coords[1] * b[2] - this->coords[2] * b[1];
    res[1] = this->coords[2] * b[0] - this->coords[0] * b[2];
    res[2] = this->coords[0] * b[1] - this->coords[1] * b[0];
    return res;
  }

  //! unary minus
  Vec3<_DataType> operator- () const {
    return ( Vec3<_DataType>( - this->coords[0], - this->coords[1], - this->coords[2] ) );
  }

  Vec3<_DataType> normalized ( ) const {
    Vec3<_DataType> res ( *this );
    res.normalize();
    return ( res );
  }

  istream& read ( istream &is ) {
    char c;
    is >> c;
    if ( c != '(' ) is.putback ( c );
    is >> xref ();
    is >> c;
    if ( c != ';' ) is.putback ( c );
    is >> yref ();
    is >> c;
    if ( c != ';' ) is.putback ( c );
    is >> zref ();
    is >> c;
    if ( c != ')' ) is.putback ( c );
    return is;
  }
};

/** component-wise minimum of two Vec3 returned as Vec3 */
template<typename DataType>
inline aol::Vec3<DataType> CompWiseMin ( const aol::Vec3<DataType>& vecA, const aol::Vec3<DataType>& vecB ) {
  aol::Vec3<DataType> ret;
  for ( unsigned short d = 0; d < 3; ++d ) {
    ret[d] = aol::Min ( vecA[d], vecB[d] );
  }
  return ( ret );
}

/** component-wise maximum of two Vec3 returned as Vec3 */
template<typename DataType>
inline aol::Vec3<DataType> CompWiseMax ( const aol::Vec3<DataType>& vecA, const aol::Vec3<DataType>& vecB ) {
  aol::Vec3<DataType> ret;
  for ( unsigned short d = 0; d < 3; ++d ) {
    ret[d] = aol::Max ( vecA[d], vecB[d] );
  }
  return ( ret );
}

/** component-wise multiplication of two Vec3 returned as Vec3 */
template< typename RealType >
inline aol::Vec3<RealType> CompWiseMultiply ( const aol::Vec3<RealType> &Arg0, const aol::Vec3<RealType> &Arg1 ) {
  aol::Vec3<RealType> ret;
  for ( short d = 0; d < 3; ++d ) {
    ret[d] = Arg0[d] * Arg1[d];
  }
  return ( ret );
}

/** component-wise division (regular or integer) of two Vec3 returned as Vec3 */
template< typename RealType >
inline aol::Vec3<RealType> CompWiseDivide ( const aol::Vec3<RealType> &Arg0, const aol::Vec3<RealType> &Arg1 ) {
  aol::Vec3<RealType> ret;
  for ( short d = 0; d < 3; ++d ) {
    ret[d] = Arg0[d] / Arg1[d];
  }
  return ( ret );
}


/////////////////////////////////////////////////////////
// Barycentric coordinates
/////////////////////////////////////////////////////////

template <int Dim, typename RealType>
class BarCoord {};

template <typename RealType>
class BarCoord<2, RealType> : public Vec<3, RealType> {
public:
  BarCoord<2, RealType> () {
    (*this)[0] = ZOTrait<RealType>::one;
    (*this)[1] = ZOTrait<RealType>::zero;
    (*this)[2] = ZOTrait<RealType>::zero;
  }
  BarCoord<2, RealType> ( RealType x0, RealType x1, RealType x2 ) {
    (*this)[0] = x0;
    (*this)[1] = x1;
    (*this)[2] = x2;
  }

  explicit BarCoord<2, RealType> ( int numOneCoord ) {
    this->setAll ( ZOTrait<RealType>::zero );
    (*this)[numOneCoord] = ZOTrait<RealType>::one;
  }

  bool addsToOne () const {
    return ( fabs ( this->sum() - ZOTrait<RealType>::one ) < 1E-14 );
  }

  // checks whether sum is one and all coordinates are within [0,1]
  bool checkForConsistency () const {
    if ( ! this->addsToOne() )
      return false;
    for ( unsigned short int i = 0; i < 3; ++i ) {
      if ( (*this)[i] < ZOTrait<RealType>::zero || (*this)[i] > ZOTrait<RealType>::one )
        return false;
    }
    return true;
  }
 };

template <typename RealType>
class BarCoord<3, RealType> : public Vec<4, RealType> {
public:
  BarCoord<3, RealType> () {
    (*this)[0] = ZOTrait<RealType>::one;
    (*this)[1] = ZOTrait<RealType>::zero;
    (*this)[2] = ZOTrait<RealType>::zero;
    (*this)[3] = ZOTrait<RealType>::zero;
  }
  BarCoord<3, RealType> ( RealType x0, RealType x1,
                          RealType x2, RealType x3 ) {
    (*this)[0] = x0;
    (*this)[1] = x1;
    (*this)[2] = x2;
    (*this)[3] = x3;
  }

  explicit BarCoord<3, RealType> ( int numOneCoord ) {
    this->setAll ( ZOTrait<RealType>::zero );
    (*this)[numOneCoord] = ZOTrait<RealType>::one;
  }

  bool addsToOne () const {
    return ( fabs ( this->sum() - ZOTrait<RealType>::one ) < 1E-14 );
  }

  // checks whether sum is one and all coordinates are within [0,1]
  bool checkForConsistency () const {
    if ( ! this->addsToOne() )
      return false;
    for ( unsigned short int i = 0; i < 4; ++i ) {
      if ( (*this)[i] < ZOTrait<RealType>::zero || (*this)[i] > ZOTrait<RealType>::one )
        return false;
    }
    return true;
  }
};

/////////////////////////////////////////////////////////
// VecDimTrait
/////////////////////////////////////////////////////////

template <typename DataType, int Dim> struct VecDimTrait {
  typedef Vec< Dim, DataType> VecType;
  typedef BarCoord<Dim, DataType> BarCoordType;
};

template <typename DataType> struct VecDimTrait<DataType, 2> {
  typedef Vec2<DataType> VecType;
  typedef BarCoord<2, DataType> BarCoordType;
};

template <typename DataType> struct VecDimTrait<DataType, 3> {
  typedef Vec3<DataType> VecType;
  typedef BarCoord<3, DataType> BarCoordType;
};

//! left-multiplication scalar * Vec3
template<typename _DataType>
Vec3<_DataType> operator* ( const _DataType alpha, const Vec3<_DataType> &vec ) {
  return ( Vec3<_DataType> ( alpha * vec[0], alpha * vec[1], alpha * vec[2] ) );
}

/////////////////////////////////////////////////////////
// Input/Output operators
/////////////////////////////////////////////////////////


template <int dim, class T>
inline ostream &operator<< ( ostream &os, const Vec<dim, T> &v ) {
  return v.print ( os );
}


template <class T>
inline ostream &operator<< ( ostream &os, const Vec2<T> &v ) {
  return v.print ( os );
}

template <int dim, class T>
inline istream &operator>> ( istream &is, Vec<dim, T> &v ) {
  return v.read ( is );
}

template <class T>
inline istream &operator>> ( istream &is, Vec2<T> &v ) {
  return v.read ( is );
}

template <class T>
inline ostream &operator<< ( ostream &os, const Vec3<T> &v ) {
  return v.print ( os );
}

template <class T>
inline istream &operator>> ( istream &is, Vec3<T> &v ) {
  return v.read ( is );
}

//! return component-wise minimum of two aol::Vec<dimension,DataType>
template <int dimension, typename DataType>
aol::Vec<dimension, DataType> CompWiseMin ( const aol::Vec<dimension, DataType>& vecA, const aol::Vec<dimension, DataType>& vecB ) {
  aol::Vec<dimension, DataType> ret;
  for ( unsigned int d = 0; d < dimension; ++d ) {
    ret[d] = aol::Min ( vecA[d], vecB[d] );
  }
  return ( ret );
}


//! return component-wise maximum of two aol::Vec<dimension,DataType>
template <int dimension, typename DataType>
aol::Vec<dimension, DataType> CompWiseMax ( const aol::Vec<dimension, DataType>& vecA, const aol::Vec<dimension, DataType>& vecB ) {
  aol::Vec<dimension, DataType> ret;
  for ( unsigned int d = 0; d < dimension; ++d ) {
    ret[d] = aol::Max ( vecA[d], vecB[d] );
  }
  return ( ret );
}

} // end of namespace aol.

#endif
