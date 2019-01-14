#ifndef __VEC_H
#define __VEC_H

#include <aol.h>
#include <quoc.h>
#include <qmException.h>
#include <smallVec.h>
#include <bitVector.h>
#include <memoryManager.h>

// forward declaration
namespace qc {
class GridStructure;
namespace simplex {
template <typename CubicGridType, qc::Dimension Dim>
class GridStructure;
}
}

namespace aol {

// This class is needed to properly initialize overflowMax in Vector
// in case DataType = signed char, because a signed char can't exceed 127
template <class T> struct OverflowTrait {
  static const T max;
};

template <typename T> class Solver;
template <typename T> class MultiVector;

//! this class is needed when you use operators with DestType = a scalar value
template <typename _DataType>
class Scalar {
public:
  typedef _DataType DataType;

  Scalar ( DataType s = 0 ) : v ( s ) {}

  //! Sets value to zero
  void setZero () {
    v = 0;
  }

  operator DataType () {
    return v;
  }

  Scalar& operator += ( Scalar s ) {
    v += s.v;
    return *this;
  }
  Scalar& operator += ( DataType s ) {
    v += s;
    return *this;
  }

  Scalar& operator -= ( Scalar s ) {
    v -= s.v;
    return *this;
  }
  Scalar& operator -= ( DataType s ) {
    v -= s;
    return *this;
  }

  Scalar& operator *= ( Scalar s ) {
    v *= s.v;
    return *this;
  }
  Scalar& operator /= ( Scalar s ) {
    v /= s.v;
    return *this;
  }

  void add ( Scalar s ) {
    *this += s;
  }

  void addMultiple ( Scalar s, DataType factor ) {
    *this += s * factor;
  }

  inline DataType& operator[] ( int ) {
    return v;
  }
  inline const DataType& operator[] ( int ) const {
    return v;
  }

  DataType v;
};

//! A vector class.
template< typename _DataType >
class Vector : public Obj {
public:
  typedef _DataType DataType;
  static const int Depth = 1;

  //! RealType is used for all operations that require floating point numbers, like calculating the norm of a vector of integers.
  typedef typename RealTrait<DataType>::RealType RealType;

protected:
  int _size, _sizeReserved;
  //! Should the data pointer be deleted during destruction of the class?
  bool _deleteFlag;

  //! contains the entries of the vector
  DataType *_pData;

  aol::OverflowHandlingType  overflowHandling; //!< Type of overflow handling
  DataType                   overflowMin;      //!< Lower bound for overflow handling
  DataType                   overflowMax;      //!< Upper bound for overflow handling

public:
  static bool quietMode;  //< Be verbose or be quiet?

private:
  static bool PrettyFormat;
  static const Format* _pFormat;

public:
  //! Allocate vector of size Length. By convention, vectors (and thus derived classes) are initialized with zero.
  explicit Vector ( int Length ) : _size ( Length ), _sizeReserved ( Length ), _deleteFlag ( true ),
      overflowHandling ( aol::CLIP ), overflowMin ( 0 ), overflowMax ( OverflowTrait<DataType>::max ) {
    _pData = static_cast<DataType*> ( aol::MemoryManager::allocateAtLeast ( _sizeReserved, sizeof ( DataType ) ) );
    // A NULL pointer is ok if we only requested a size of 0. We still need to make sure that at least as much
    // memory as requested was allocated though.
    if ( ( !_pData && _sizeReserved ) || ( _sizeReserved < _size ) )
      throw OutOfMemoryException ( aol::strprintf( "Vector<DataType,Realtype>::Vector: Could not allocate memory for Vector of size. We asked for %d, but got %d", Length, _sizeReserved ) , __FILE__, __LINE__ );

    setZero ();
  }

  //! Allocate vector of size Length and initialize with ConstValue.
  Vector ( int Length, const DataType& ConstValue ) : _size ( Length ), _sizeReserved ( Length ), _deleteFlag ( true ),
      overflowHandling ( aol::CLIP ), overflowMin ( 0 ), overflowMax ( OverflowTrait<DataType>::max ) {
    _pData = static_cast<DataType*> ( aol::MemoryManager::allocateAtLeast ( _sizeReserved, sizeof ( DataType ) ) );
    // A NULL pointer is ok if we only requested a size of 0. We still need to make sure that at least as much
    // memory as requested was allocated though.
    if ( ( !_pData && _sizeReserved ) || ( _sizeReserved < _size ) )
      throw OutOfMemoryException ( "Vector<DataType,Realtype>::Vector: Could not allocate memory for Vector.", __FILE__, __LINE__ );

    setAll ( ConstValue );
  }

  //! Instanciate a vector wrapper for the field Data of N DataTypes.
  Vector ( DataType *Data, int Size, CopyFlag copyFlag = FLAT_COPY ); // default: flat copy - does deep copy make sense? Maybe ...

  //! Copy from multi-vector
  explicit Vector ( const MultiVector<DataType>& mv ); // default: deep copy - flat does not make sense.

  //! constructor from STL vector
  explicit Vector ( const std::vector<DataType> & std_vec );

  //! conversion from aol::Vec<N, DataType>
  template < int N >
  explicit Vector ( const aol::Vec< N, DataType > &SmallVec ) : _size ( N ), _sizeReserved ( N ), _deleteFlag ( true ), _pData ( NULL ),
    overflowHandling ( aol::CLIP ), overflowMin ( 0 ), overflowMax ( aol::OverflowTrait<DataType>::max ) { // must be implemented here in header due to additional template
    _pData = static_cast< DataType* > ( aol::MemoryManager::allocateAtLeast ( _sizeReserved, sizeof ( DataType ) ) );
    for ( int i = 0; i < _size; ++i )
      _pData[i] = SmallVec[i];
  }

  /** \brief Copy constuctor
   *  \attention This constructor will not be used on call-by-value (and other implicit copies), since it is explicit.
   *  Luckily, it also prevents the usage of the automatically generated copy constructor, so that call-by-value
   *  cannot be used for this class.
   */
  // Unfortunately older versions of GCC don't allow aol::Vector to be used inside an std::vector if the copy constructor is explicit.
  // Furthermore, using aol::Vector as "firstprivate" with OpenMP and VC++ doesn't work in this case.
#if ( !defined ( __GNUC__ ) || ( GCC_VERSION >= 40200 ) ) && ( !defined ( _MSC_VER ) || !defined ( _OPENMP ) )
  explicit
#endif
  Vector ( const Vector<DataType> &Vec, CopyFlag copyFlag = DEEP_COPY );

  //! Standard constructor
  Vector() : _size ( 0 ), _sizeReserved ( 0 ), _deleteFlag ( false ), _pData ( NULL ),
      overflowHandling ( aol::CLIP ), overflowMin ( 0 ), overflowMax ( OverflowTrait<DataType>::max ) {}


  //! Constructor creating Vector of correct size for a qc::GridDefinition
  explicit Vector ( const qc::GridStructure &Grid );

  //! Constructor creating Vector of correct size for a simplex grid
  template <typename CubicGridType, qc::Dimension Dim>
  explicit Vector ( const qc::simplex::GridStructure<CubicGridType, Dim> & Grid )
      : _size ( Grid.getNumberOfNodes() ), _sizeReserved ( _size ), _deleteFlag ( true ),
      overflowHandling ( aol::CLIP ), overflowMin ( 0 ), overflowMax ( OverflowTrait<DataType>::max ) {
    _pData = static_cast< DataType* > ( aol::MemoryManager::allocateAtLeast ( _sizeReserved, sizeof ( DataType ) ) );
    if ( !_pData )
      throw OutOfMemoryException ( "Vector<DataType,Realtype>::Vector: Could not allocate memory for Vector.", __FILE__, __LINE__ );

    setZero ();
  }

  virtual ~Vector() {
    if ( _deleteFlag && _pData )
      aol::MemoryManager::deallocate ( _pData, _sizeReserved, sizeof ( DataType ) );
  }


  /** reserve the vector and do not delete contents.
   * \author Notthoff
   */
  void reserve ( int Length );

  void reserveUninit ( int Length );

  /** Change size of the vector, keep old contents (at least those that still fit when shrinking) and set new contents to zero (when growing).
   *  Attention: This function may or may not make sense on derived classes.
   *  \author Notthoff
   */
  virtual void resize ( const int Length );

  /** Adds Data to the end of the vector and changes the size of the vector if necessary
   * \author Teusner
   */
  virtual void pushBack ( const _DataType Data );

  /** Appends contents of otherVector to this Vector
   * \author Schwen
   */
  virtual void pushBackValues ( const aol::Vector<_DataType> &otherVector );

  /** Change size of the vector, setting all entries to zero.
   *  \author Schwen
   */
  virtual void reallocate ( const int Length );

  void reallocate ( const Vector<DataType> &other ) {
    this->reallocate ( other.size() );
  }

  void reallocate ( const qc::GridStructure &grid );

  //! change size of the vector, setting all entries to zero and freeing old / allocating new memory (in order to reduce memory consumption)
  virtual void reallocateClear ( const int Length );

  /** Increase (or decrease) the size of this vector via resize */
  virtual void growBy ( const int addlLength ) {
    resize ( _size + addlLength );
  }

  //! Access the Ith element
  inline DataType& operator[] ( int I ) {
#ifdef BOUNDS_CHECK
    if ( I < 0 || I >= _size ) {
      char error[1024];
      sprintf ( error, "Vector<DataType,Realtype>::operator[]: index %d out of bounds. _size = %d\n", I, _size );
      throw OutOfBoundsException ( error, __FILE__, __LINE__ );
    }
#endif
    return _pData[I];
  }

  //! Access the Ith element
  inline const DataType& operator[] ( int I ) const {
#ifdef BOUNDS_CHECK
    if ( I < 0 || I >= _size ) {
      char error[1024];
      sprintf ( error, "Vector<DataType,Realtype>::operator[]: index %d out of bounds. _size = %d\n", I, _size );
      throw OutOfBoundsException ( error, __FILE__, __LINE__ );
    }
#endif
    return _pData[I];
  }


  /** Very rudimentary write without storing size
   */
  void saveAsVector ( const char *fileName ) const {
    FILE *out;

    if ( ! ( out = fopen ( fileName, "wb" ) ) ) {
      throw Exception ( "Cannot open file for writing", __FILE__, __LINE__ );
    }

    fwrite ( _pData, sizeof ( DataType ), size(), out );
    fclose ( out );
  }

  /** Very rudimentary read, size must be set in advance
   */
  void loadAsVector ( const char *fileName ) {
    FILE *in;

    if ( ! ( in = fopen ( fileName, "rb" ) ) ) {
      throw Exception ( "Cannot open file for reading", __FILE__, __LINE__ );
    }

    const int result = fread ( _pData, sizeof ( DataType ), size(), in );
    if ( result != size() )
      throw IOException ( "Error reading from file", __FILE__, __LINE__ );

    fclose ( in );
  }

  //! Saves values into ASCII-file.
  void saveASCII ( const char *fileName, int prec = 3 ) const {
    ofstream file ( fileName );
    file.precision( prec );

    for ( int i = 0; i < this->size(); ++i )
      file << this->operator[] ( i ) << endl;

    file.close();
  }

  /** Sets the handling of data which cannot be stored in the standard PGM format.
   *  This settings take effect if you call save with type 5 subsequently. The default value is CLIP
   *  @param  type The type of overflow handling.
   *               CLIP clips values to given range [min, max],
   *               SCALE scales values from [minimal value, maximal value] to [0, 255] (value range of PGM),
   *               CLIP_THEN_SCALE does both. (Note that after clipping, minimal value = min, maximal value = max.)
   *  @param  min  For CLIP, CLIP_THEN_SCALE and REFLECT gives the value of the lower bound of the clipping or reflexion interval
   *  @param  max  For CLIP, CLIP_THEN_SCALE and REFLECT gives the value of the upper bound of the clipping or reflexion interval
   *  \author Preusser
   */

  void setOverflowHandling ( const aol::OverflowHandlingType type,
                             const DataType min,
                             const DataType max ) {
    if ( min >= max )
      throw aol::Exception ( "Vector<DataType>::setOverflowHandling(): min >= max",
                             __FILE__, __LINE__ );

    overflowHandling = type;

    overflowMin = min;
    overflowMax = max;
  }

  //! Sets the overflow handling to aol::CLIP_THEN_SCALE, the min value to getMinValue() and the max value to getMaxValue().
  //! In case both values are the same, the max value is set to getMinValue() + 1, since max is required to be bigger than min.
  void setOverflowHandlingToCurrentValueRange ( ) {
    const DataType minValue = getMinValue();
    const DataType maxValue = getMaxValue();
    setOverflowHandling ( aol::CLIP_THEN_SCALE, minValue, ( maxValue > minValue ) ? maxValue : minValue + 1 );
  }

  //! Sets the overflow handling to aol::CLIP_THEN_SCALE, and the min value and max values such that the specified
  //! percentage of entries gets saturated to enhance the contrast.
  void setOverflowHandlingToEnhanceContrast ( const RealType SaturatedEntryPercentage ) {
    const aol::Vec2<DataType> minMax = getSaturatedMinMaxValue ( SaturatedEntryPercentage );
    setOverflowHandling ( aol::CLIP_THEN_SCALE, minMax[0], ( minMax[1] > minMax[0] ) ? minMax[1] : minMax[0] + 1 );
  }

  bool createOverflowHandledData ( aol::Vector<unsigned char> &buffer,
                                   const DataType min,
                                   const DataType max ) const;

  bool createOverflowHandledData ( aol::Vector<unsigned char> &buffer ) const {
    const DataType min = this->getMinValue();
    const DataType max = this->getMaxValue();

    return createOverflowHandledData ( buffer, min, max );
  }

  DataType getOverflowMin ( ) const {
    return overflowMin;
  }

  DataType getOverflowMax ( ) const {
    return overflowMax;
  }

  static void setQuietMode ( bool qmode ) {
    quietMode = qmode;
  }
  
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
    PrettyFormat = newPrettyFormat;
  }

  //! Get pretty output formatting status
  static bool getPrettyFormat ( ) {
    return PrettyFormat;
  }

  /**
   * Uses createOverflowHandledData() for USIGNED_CHAR SaveTypes (which respects selected overflowHandling).
   * For UNSIGNED_SHORT, takes the absolute value, then scales [0,Maximum] to [0,65535].
   * @param Minimum in the old version of this function Minimum equaled zero and was no function parameter
   * @param Maximum it was the third and last parameter in the old version of this function and it's standard behaviour was zero
   */
  void saveRaw ( const char *filename, const qc::SaveType type, const DataType Minimum, const DataType Maximum ) const;

  /**
   * Uses createOverflowHandledData() for USIGNED_CHAR SaveTypes (which respects selected overflowHandling).
   * For UNSIGNED_SHORT, takes the absolute value, then scales [0,Maximum] to [0,65535].
   * @param Minimum in the old version of this function Minimum equaled zero and was no function parameter
   * @param Maximum it was the third and last parameter in the old version of this function and it's standard behaviour was zero
   */
  void saveRaw ( ostream &out, const qc::SaveType type, const DataType Minimum, const DataType Maximum ) const;

  void loadRaw ( istream &in, const int Type );

  //! lossless save in binary quocmesh format
  void saveToFile ( const char *filename ) const;

  //! load from file in quocmesh format, no conversion of data types
  void loadFromFile ( const char *filename );

  //! load from file in quocmesh format with different DataType
  template< typename LoadDataType >
  void loadConvertFromFile ( const char *filename ) {
    aol::Vector<LoadDataType> loadVec;
    loadVec.loadFromFile ( filename );
    this->reallocate ( loadVec.size() );
    this->convertFrom ( loadVec );
  }

  DataType get ( const int I ) const {
#ifdef BOUNDS_CHECK
    if ( I < 0 || I >= _size ) {
      char error[1024];
      sprintf ( error, "Vector<DataType>::get: index %d out of bounds. _size = %d\n", I, _size );
      throw OutOfBoundsException ( error, __FILE__, __LINE__ );
    }
#endif
    return _pData[I];
  }

  void set ( int I, DataType val ) {
#ifdef BOUNDS_CHECK
    if ( I < 0 || I >= _size ) {
      char error[1024];
      sprintf ( error, "Vector<DataType>::set: index %d out of bounds. _size = %d\n", I, _size );
      throw OutOfBoundsException ( error, __FILE__, __LINE__ );
    }
#endif
    _pData[I] = val;
  }

  void apply ( DataType f ( const DataType ) ) {
    const int s = size();
    for ( int i = 0; i < s; ++i ) {
      const DataType v = ( *this ) [i];
      ( *this ) [i] = f ( v );
    }
  }

  void add ( int I, DataType val ) {
#ifdef BOUNDS_CHECK
    if ( I < 0 || I >= _size ) {
      char error[1024];
      sprintf ( error, "Vector<DataType,Realtype>::add: index %d out of bounds. _size = %d\n", I, _size );
      throw OutOfBoundsException ( error, __FILE__, __LINE__ );
    }
#endif
    _pData[I] += val;
  }

  RealType interpolate ( const RealType pos ) const {
    if ( pos < 0 || pos > _size - 1 )
      throw OutOfBoundsException ( "aol::Vector<DataType>::interpolate: position out of range", __FILE__, __LINE__ );

    const int p = static_cast<int> ( pos );
    const RealType pf = pos - p;

    if ( p == _size - 1 )
      return ( static_cast<RealType> ( this->get ( p ) ) );
    else
      return ( ( 1 - pf ) * this->get ( p ) + pf * this->get ( p + 1 ) );
  }

  RealType interpolateInRange ( const RealType Pos, const RealType Rmin, const RealType Rmax ) const {
    if ( Pos < Rmin || Pos > Rmax )
      throw OutOfBoundsException ( "aol::Vector<DataType>::interpolateInRange: position out of range", __FILE__, __LINE__ );

    const RealType bary = ( _size - 1 ) * ( Pos - Rmin ) / ( Rmax - Rmin );
    return interpolate ( bary );
  }

  //! Returns the size of the vector.
  // (NB: Calling this method size() follows the STL-conventions.)
  int size() const {
    return _size;
  }

  //! Returns the size of the vector. For compatibility with getTotalSize() from MultiVector.
  int getTotalSize() const {
    return size();
  }

  //! Returns the capacity of the vector.
  // (NB: Calling this method capacity() follows the STL-conventions.)
  int capacity() const {
    return _sizeReserved;
  }

  bool getDeleteFlag() const {
    return _deleteFlag;
  }

  //! Dot product
  virtual DataType operator* ( const Vector<DataType> &c ) const;

  //! Dot product
  DataType dotProduct ( const Vector<DataType> &Vec ) const;

  //! Print the vector to stdout.
  //! deprecated? use operator << and/or print instead? Note that those methods cannot print interactively!
  void dump ( bool = false ) const ;

  //! Print the vector to os
  ostream& print ( ostream& os ) const {
    os << ( PrettyFormat ? "[ " : "" );
    for ( int i = 0; i < size (); ++i )
      if ( i + 1 < size () ) os << (*_pFormat) ( ( *this ) [i] ) << ( PrettyFormat ? " ; " : " " );
      else os << (*_pFormat) ( ( *this ) [i] );
    return os << ( PrettyFormat ? " ]" : "" );
  }

  istream& read ( istream &is ) {
    resize ( 1 );
    char c = 0;
    DataType value;
    is >> c;
    if ( c != '(' && c != '[' ) is.putback ( c );
    is >> ( *this ) [0];
    while ( !is.eof() ) {
      is >> c;
      if ( c == ')' || c == ']' || is.eof() )
        return is;
      if ( c != ';' && c != ',' )
        is.putback ( c );
      is >> value;
      pushBack ( value );
    }
    return is;
  }

  //! Reads in ASCII files containing lists of numbers
  void read ( const char* FileName ) {
    ifstream fileStream ( FileName );
    read ( fileStream );
    fileStream.close();
  }

  //! Sets all values to zero
  void setZero () {
    memset ( _pData, 0, sizeof ( DataType ) * _size );
  }

  template < int N >
  void copyTo ( aol::Vec<N, DataType> &Dest ) const {
    QUOC_ASSERT ( size() == N );
    for ( int i = 0; i < N; ++i )
      Dest [i] = (*this)[i];
  }

  void copyToBuffer ( DataType* buffer ) const;

  void readFromBuffer ( const DataType* const buffer );
  
  //! \brief Creates a resampled version using linear interpolation
  void resampleFrom ( const aol::Vector<DataType> &Other ) {
    const int n = size ( ), nOther = Other.size ( );
    const RealType s = ( nOther - 1 ) / static_cast<RealType> ( n - 1 );
    for ( int i=0; i<n ; ++i )
      (*this)[i] = Other.interpolate ( i * s );
  }
  
  //! \brief Creates a resampled version assuming other vector is a piece-wise constant function on [0,1] and a constant point spread function
  void resampleFromPiecewiseConstant ( const aol::Vector<DataType> &Other ) {
    const int n = size ( ), nOther = Other.size ( );
    const RealType s = nOther / static_cast<RealType> ( n );
    for ( int i=0; i<n ; ++i ) {
      const RealType xl = i * s, xr = ( i + 1 ) * s;
      const int xlFloor = floor ( xl ), xrFloor = floor ( xr );
      if ( xlFloor == xrFloor ) (*this)[i] = ( xr - xl ) * Other[xlFloor];
      else {
        (*this)[i] = ( 1 - ( xl - xlFloor ) ) * Other[xlFloor];
        for ( int j=xlFloor+1; j<xrFloor ; ++j ) (*this)[i] += Other[j];
        (*this)[i] += ( xrFloor < nOther ) ? ( xr - xrFloor ) * Other[xrFloor] : 0;
      }
    }
    (*this) /= s;
  }

  //! \brief Assigning a vector to another.
  //!
  //! Sizes must match (if necessary, resize before assigning).
  //! Assignment is not virtual.
  //! Cf. S. Meyers: Effective C++, Item 33.
  //! If it would be virtual, you would need to write a operator= on e.g. ScalarArray<QC_2D> that accepts a Vector as right hand side,
  //! because only this version would overwrite Vector's operator=. This operator= would then have to check (e.g. by dynamic_cast)
  //! whether the RHS is in fact a ScalarArray<QC_2D>. Only this would allow to assign a ScalarArray<QC_2D> to a ScalarArray<QC_2D> correctly,
  //! if both are only available as Vector&, but this solution is deemed too complicated.
  Vector<DataType> & operator= ( const Vector<DataType> &Vec );

  //! Assigning a vector to another.
  Vector<DataType>& assignFrom ( const Vector<DataType> &Vec ) {
    operator= ( Vec );
    return *this;
  }

  Vector<DataType> & assignFrom ( const BitVector & bitField );
  Vector<DataType> & copyUnblockedFrom ( const MultiVector<DataType> & multiVector );

  //! Assign only values that are marked "true" (if invertMask == false)
  Vector<DataType> & assignMaskedFrom ( const Vector<DataType> & Vec, const BitVector & mask, bool invertMask = false ) {
    DataType *ptrThis = _pData;
    DataType *ptrRhs  = Vec._pData;
    // Note: ^ is the C notation for XOR. With a short look at the logic table,
    // you can see that (a XOR b) equals a, if b == true, and not(a) if b == false.
    for ( int i = 0; i < _size; ++i ) {
      if ( mask[i] ^ invertMask )
        *ptrThis = *ptrRhs;
      ptrThis++;
      ptrRhs++;
    }
    return *this;
  }

  //! Set all components of the vector to the absolute value of the component.
  void setEntriesToAbs ( ) {
    for ( int i = 0; i < _size; ++i ) _pData[ i ] = Abs ( _pData[i] );
  }

  //! Subtracting vectors.
  //virtual Vector<DataType> operator-( const Vector<DataType> &Vec );

  //! Adding a scalar to all components of the vector
  virtual Vector<DataType>& addToAll ( DataType Scalar ) {
    const int n = _size;
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
    for ( int i = 0; i < n; ++i ) _pData[ i ] += Scalar;
    return *this;
  }

  //! Multiplies each components of this vector with the corresponding component of the argument vector Vec
  void entrywiseMultiplyFrom ( const Vector<DataType> &Vec ) {
    if ( Vec.size() == _size ) {
      for ( int i = 0; i < _size; ++i ) {
        _pData[i] *= Vec.get ( i );
      }
    } else {
      throw aol::Exception ( "aol::Vector::entrywiseMultiply dimensions don't match", __FILE__, __LINE__ );
    }
  }

  //! Subtracting a vector from this vector
  virtual Vector<DataType>& operator-= ( const Vector<DataType> &Vec );

  // operator+ has to create a temporary object (same as operator-),
  // thus it is slow. We agreed not to implement it as it will
  // certainly be used (by accident) at time-critical locations once
  // it is implemented.
  // virtual Vector<DataType> operator+( const Vector<DataType> &Vec );

  //! Adding a vector to this vector
  virtual Vector<DataType>& operator+= ( const Vector<DataType> &Vec );

  //! Multiplication with a scalar
  virtual Vector<DataType>& operator*= ( const DataType Value );

  //! Division by a scalar
  virtual Vector<DataType>& operator/= ( const DataType Value );

  //! Multiplication with a scalar only on nodes where
  //! a given mask is set "true" (in case invertMask == false)
  Vector<DataType> & multMasked ( const DataType Value, const BitVector & mask, bool invertMask = false );

  //! Multiply every component of this vector by the corresponding component of Vec
  virtual Vector<DataType>& operator*= ( const Vector<DataType> &Vec ) {
    const int n = _size;
    for ( int i = 0; i < n; ++i ) {
      _pData[i] *= Vec._pData[i];
    }
    return *this;
  }

  //! Divide every component of this vector by the corresponding component of Vec
  virtual Vector<DataType>& operator/= ( const Vector<DataType> &Vec ) {
    const int n = _size;
    for ( int i = 0; i < n; ++i ) _pData[i] /= Vec._pData[i];
    return *this;
  }

  virtual bool operator== ( const Vector<DataType> &Vec ) const {
    for ( int i = 0; i < size(); ++i ) {
      if ( Vec[ i ] != _pData[ i ] ) {
        return false;
      }
    }
    return true;
  }

  virtual bool operator!= ( const Vector<DataType> &Vec ) const {
    return ! ( *this == Vec );
  }

  //! lexicographic comparison
  bool operator < ( const Vector<DataType>& other ) const {
    for ( int i = 0; i < aol::Min( this->size(), other.size() ); i++ )
      if ( _pData[i] != other[i] ) return _pData[i] < other[i];
    return ( this->size() < other.size() );
  }

  //! add only components that are marked as "true" (if invertMask == false)
  Vector<DataType> & addMasked ( const Vector<DataType> & Vec, const BitVector & mask, bool invertMask = false );

  //! Add Factor*Vec to this vector.
  Vector<DataType> & addMultiple ( const Vector<DataType>& Vec, DataType Factor );
  Vector<DataType> & addMultipleMasked ( const Vector<DataType> & Vec, DataType Factor, const BitVector & mask, bool invertMask = false );

  Vector<DataType> & scaleAndAdd ( const DataType Factor, const aol::Vector<DataType>& Vec ) {
    if ( Vec.size() != _size )
      throw aol::Exception ( "aol::Vector::scaleAndAdd(): dimensions don't match", __FILE__, __LINE__ );

    for ( int i = 0; i < _size; ++i )
      _pData[ i ] = Factor * _pData[ i ] + Vec[ i ];

    return *this;
  }

  Vector<DataType> & scaleAndAddMultiple ( const DataType ScaleFactor, const aol::Vector<DataType>& Vec, const DataType VecFactor ) {
    if ( Vec.size() != _size )
      throw aol::Exception ( "aol::Vector::scaleAndAdd(): dimensions don't match", __FILE__, __LINE__ );

    for ( int i = 0; i < _size; ++i )
      _pData[ i ] = ScaleFactor * _pData[ i ] + VecFactor * Vec[ i ];

    return *this;
  }

  void setSum ( const Vector<DataType>& Vec1, const Vector<DataType>& Vec2, DataType Factor );

  //! Makes this vector the absolute difference of Vec1 and Vec2.
  void absDiff ( Vector<DataType> const& Vec1,
                 Vector<DataType> const& Vec2 );

  //! the euclidian norm
  RealType norm() const {
    return sqrt ( static_cast<RealType> ( normSqr() ) );
  }


  //! the euclidian norm squared
  //  RealType normSqr( ) const { return lpNormPowP(2); }
  DataType normSqr() const {
    return ( ( *this ) * ( *this ) );
  }

  /** Returns \f$ \|\cdot\|_p^p \f$ which is the p-th power of the \f$l^p\f$ norm of this vector
   * \author Preusser
  */
  RealType lpNormPowP ( RealType p = 2.0 ) const;

  /** Returns \f$ \|\cdot\|_p \f$ which is the \f$l^p\f$ norm of this vector
    * \author Preusser
    */
  RealType lpNorm ( RealType p = 2.0 ) const;

  /** Returns \f$ \sum_{i=1}^N x_i^p \f$.
   * \author Preusser
   */
  DataType sum ( int p = 1 ) const;

  /** Masked equivalent of sum for p=1.
   * \author Berkels
   */
  DataType sumMasked ( const BitVector &Mask, bool InvertMask = false ) const {
    DataType *ptrThis = _pData;
    DataType sum = 0;
    for ( int i = 0; i < _size; ++i ) {
      if ( Mask[i] ^ InvertMask )
        sum += *ptrThis;
      ++ptrThis;
    }
    return sum;
  }

  /** Returns \f$ \sum_{i=1}^N x_i y_i \f$.
   * \author Geihe
   */
  DataType sumWeighted ( Vector<DataType> const& weight ) const;

  //! Returns the minimal value in this vector.
  DataType getMinValue() const {
    DataType min = get ( 0 ); // this implementation is intended: if size is 0, bounds checking throws exception.
    for ( int i = 1; i < _size; ++i ) {
      if ( _pData[ i ] < min ) min = _pData[ i ];
    }
    return min;
  }

  //! return std::pair of index and value of minimal entry
  std::pair< int, DataType > getMinIndexAndValue ( ) const {
    std::pair< int, DataType > IndVal ( 0, get ( 0 ) );
    for ( int i = 1; i < _size; ++i ) {
      if ( _pData[ i ] < IndVal.second ) {
        IndVal.first = i;
        IndVal.second = _pData[ i ];
      }
    }
    return ( IndVal );
  }

  //! Returns the maximum value in this vector.
  DataType getMaxValue() const {
    DataType max = get ( 0 );
    for ( int i = 1; i < _size; ++i ) {
      if ( _pData[ i ] > max ) max = _pData[ i ];
    }
    return max;
  }

  //! return std::pair of index and value of maximal entry
  std::pair< int, DataType > getMaxIndexAndValue() const {
    std::pair< int, DataType > IndVal ( 0, get ( 0 ) );
    for ( int i = 1; i < _size; ++i ) {
      if ( _pData[ i ] > IndVal.second ) {
        IndVal.first = i;
        IndVal.second = _pData[ i ];
      }
    }
    return ( IndVal );
  }

  //! Returns a Vec2 with minimal value in the first, maximal value
  //! in the second component
  Vec2<DataType> getMinMaxValue() const {
    DataType min = get ( 0 ), max = get ( 0 );
    for ( int i = 1; i < _size; ++i ) {
      if ( _pData[ i ] > max ) max = _pData[ i ];
      if ( _pData[ i ] < min ) min = _pData[ i ];
    }
    return Vec2<DataType> ( min, max );
  }

  //! Returns the minimum absolute value in this vector (not a norm)
  DataType getMinAbsValue() const {
    DataType min = Abs ( get ( 0 ) );
    for ( int i = 1; i < _size; ++i ) {
      if ( Abs ( _pData[i] ) < min ) min = Abs ( _pData[i] );
    }
    return min;
  }

  //! Returns the maximum absolute value in this vector (l-infinity norm)
  DataType getMaxAbsValue() const {
    DataType max = Abs ( get ( 0 ) );
    for ( int i = 1; i < _size; ++i ) {
      if ( Abs ( _pData[i] ) > max ) max = Abs ( _pData[i] );
    }
    return max;
  }

  //! Returns min and max values such that half of the specified percentage of the entries
  //! are lower than min and half bigger than max. For instance useful to enhance the contrast.
  //! As percentage the expected range for the argument is 0 to 100.
  Vec2<DataType> getSaturatedMinMaxValue( const RealType SaturatedEntryPercentage ) const {
    aol::Vector<DataType> temp ( *this );
    temp.sortValues();
    const int sizeM1 = temp.size()-1;
    const DataType minValue = temp [ aol::Clamp ( static_cast<int> ( SaturatedEntryPercentage / 200 * sizeM1 ), 0, sizeM1 ) ];
    const DataType maxValue = temp[ aol::Clamp ( static_cast<int> ( ( 1 - SaturatedEntryPercentage / 200 ) * sizeM1 ), 0, sizeM1 ) ];
    return Vec2<DataType> ( minValue, maxValue );
  }
  
  //! Returns the minimum value among all finite values in this vector.
  DataType getFiniteMinValue() const {
    int i0 = 0;
    while ( !aol::isFinite<DataType> ( _pData [ i0 ] ) && i0 < _size ) ++i0;
    if ( i0 == _size ) throw aol::Exception ( "Vector does not contain any finite entries!", __FILE__, __LINE__ );
    else {
      DataType min = _pData [ i0 ];
      for ( int i = i0+1; i < _size; ++i ) if ( aol::isFinite<DataType> ( _pData[ i ] ) && _pData[ i ] < min ) min = _pData[ i ];
      return min;
    }
  }
  
  //! Returns the maximum value among all finite values in this vector.
  DataType getFiniteMaxValue() const {
    int i0 = 0;
    while ( !aol::isFinite<DataType> ( _pData [ i0 ] ) && i0 < _size ) ++i0;
    if ( i0 == _size ) throw aol::Exception ( "Vector does not contain any finite entries!", __FILE__, __LINE__ );
    else {
      DataType max = _pData [ i0 ];
      for ( int i = i0+1; i < _size; ++i ) if ( aol::isFinite<DataType> ( _pData[ i ] ) && _pData[ i ] > max ) max = _pData[ i ];
      return max;
    }
  }
  
  //! return std::pair of index and value of maximal of finite entries
  std::pair< int, DataType > getFiniteMinIndexAndValue() const {
    int i0 = 0;
    while ( !aol::isFinite<DataType> ( _pData [ i0 ] ) && i0 < _size ) ++i0;
    if ( i0 == _size ) throw aol::Exception ( "Vector does not contain any finite entries!", __FILE__, __LINE__ );
    std::pair< int, DataType > IndVal ( i0, get ( i0 ) );
    for ( int i = i0+1; i < _size; ++i ) {
      if ( _pData[ i ] < IndVal.second ) {
        IndVal.first = i;
        IndVal.second = _pData[ i ];
      }
    }
    return ( IndVal );
  }
  
  //! return std::pair of index and value of maximal of finite entries
  std::pair< int, DataType > getFiniteMaxIndexAndValue() const {
    int i0 = 0;
    while ( !aol::isFinite<DataType> ( _pData [ i0 ] ) && i0 < _size ) ++i0;
    if ( i0 == _size ) throw aol::Exception ( "Vector does not contain any finite entries!", __FILE__, __LINE__ );
    std::pair< int, DataType > IndVal ( i0, get ( i0 ) );
    for ( int i = i0+1; i < _size; ++i ) {
      if ( _pData[ i ] > IndVal.second ) {
        IndVal.first = i;
        IndVal.second = _pData[ i ];
      }
    }
    return ( IndVal );
  }

  //! Returns true if all values of the vector are in [A,B] or false otherwise. Note: A <= B is assumed.
  bool checkRange ( const DataType A = 0, const DataType B = 1 ) const {
    aol::Vec2<DataType> minMax = getMinMaxValue();
    return ( ( minMax[0] >= A ) && ( minMax[1] <= B ) );
  }

  //! return std::pair of index and value of the maximal absolute entry in this vector
  std::pair< int, DataType > getMaxAbsIndexAndValue() const {
    std::pair< int, DataType > IndVal ( 0, Abs ( get ( 0 ) ) );
    for ( int i = 1; i < _size; ++i ) {
      if ( Abs ( _pData[ i ] ) > IndVal.second ) {
        IndVal.first = i;
        IndVal.second = Abs ( _pData[ i ] );
      }
    }
    return ( IndVal );
  }

  //! Returns the minimal value in this vector.
  DataType getMinValueMasked ( const BitVector & Mask, bool invertMask = false ) const {
#ifdef BOUNDS_CHECK
    if ( Mask.numOccurence ( !invertMask ) == 0 )
      throw aol::Exception ( "Vector::getMinValueMasked(): given mask is empty.", __FILE__, __LINE__ );
#endif
    DataType min = numeric_limits<DataType>::max();
    for ( int i = 0; i < _size; ++i ) {
      if ( Mask[i] ^ invertMask )
        if ( _pData[ i ] < min ) min = _pData[ i ];
    }
    return min;
  }

  //! Returns the maximum value in this vector.
  DataType getMaxValueMasked ( const BitVector & Mask, bool invertMask = false ) const {
#ifdef BOUNDS_CHECK
    if ( Mask.numOccurence ( !invertMask ) == 0 )
      throw aol::Exception ( "Vector::getMinValueMasked(): given mask is empty.", __FILE__, __LINE__ );
#endif
    DataType max = MaxInitializerTrait<DataType>::MaxInitializer;
    for ( int i = 0; i < _size; ++i ) {
      if ( Mask[i] ^ invertMask )
        if ( _pData[ i ] > max ) max = _pData[ i ];
    }
    return max;
  }

  //! Returns the maximum absolute value in this vector (l-infinity norm)
  DataType getMaxAbsValueMasked ( const BitVector & Mask, bool invertMask = false ) const {
#ifdef BOUNDS_CHECK
    if ( Mask.numOccurence ( !invertMask ) == 0 )
      throw aol::Exception ( "Vector::getMinValueMasked(): given mask is empty.", __FILE__, __LINE__ );
#endif
    DataType max = 0;
    for ( int i = 0; i < _size; ++i ) {
      if ( Mask[i] ^ invertMask )
        if ( Abs ( _pData[i] ) > max ) max = Abs ( _pData[i] );
    }
    return max;
  }

  //! Returns the arithmetic mean value of the values of the vector.
  RealType getMeanValue() const {
    return ( static_cast<RealType> ( sum ( 1 ) ) / static_cast<RealType> ( size() ) );
  }

  //! Returns the weighted arithmetic mean value of the values of the vector.
  RealType getWeightedMeanValue( const aol::Vector<RealType> &Weights ) const;

  //! Returns the median of the values of the vector.
  //! Note that this method requires sorting and is slow.
  RealType getMedianValue() const;

  //! Returns the weighted median of the values of the vector.
  //! Note that this method requires sorting and is slow.
  RealType getWeightedMedianValue( const aol::Vector<RealType> &Weights ) const;

  //! Returns a standard deviation of the values: sqrt ( ( (v[0] - v)^2 + (v[1] - v)^2 + ... + (v[n-1] - v)^2 ) / ( n - 1 ) ) [Sic: n-1]
  RealType getStdDev() const {
    return sqrt ( getVariance() );
  }
  
  //! Returns unbiased estimator of the standard deviation of the values: 1 / kappa * sqrt ( ( (v[0] - v)^2 + (v[1] - v)^2 + ... + (v[n-1] - v)^2 ) / ( n - 1 ) ) [Sic: n-1]
  //! kappa = 1 - 1/4 n^-1 + 7/32 n^-2 + O(n^-3) \approx sqrt(2/(n-1)) * Gamma(n/2) / Gamma((n-1)/2)
  RealType getStdDevUnbiased() const {
    return ( 1.0 - 1.0 / ( 4.0 * static_cast<RealType> ( size() ) ) + 7.0 / ( 32.0 * aol::Sqr<RealType> ( static_cast<RealType> ( size() ) ) ) ) * sqrt ( getVariance() );
  }
  
  //! Returns a variance of the values: ( (v[0] - v)^2 + (v[1] - v)^2 + ... + (v[n-1] - v)^2 ) / ( n - 1 ) [Sic: n-1]
  RealType getVariance() const {
    if ( size() == 0 ) {
      return ( aol::NumberTrait<RealType>::NaN );
    } else if ( size() == 1 ) {
      return ( NumberTrait<RealType>::zero );
    } else {
      const RealType mean = getMeanValue();
      RealType sumsqr = 0.0;
      for ( int i = 0; i < size(); ++i ) {
        sumsqr += aol::Sqr ( _pData[i] - mean );
      }
      return sumsqr / ( size() - 1 );
    }
  }
  
  //! Sets all entries to Val.
  void setAll ( DataType Val ) {
    DataType *ptr = _pData;
    for ( int i = 0; i < _size; ++i ) {
      * ( ptr++ ) = Val;
    }
  }

  //! set all components that are marked as "true" (if invertMask == false)
  void setAllMasked ( DataType Val, const BitVector & mask, bool invertMask = false ) {
    DataType *ptr = _pData;
    // Note: ^ is the C notation for XOR. With a short look at the logic table,
    // you can see that (a XOR b) equals a, if b == true, and not(a) if b == false.
    for ( int i = 0; i < _size; ++i ) {
      if ( mask[i] ^ invertMask )
        * ( ptr ) = Val;
      ptr++;
    }
  }

  //! set value at position I ( needed by aol::DerivativeValidatorBase<> )
  Vec<1, int> setIthComponent ( const int I, const DataType Value ) {
    set(I, Value);
    Vec<1, int> res;
    res[0] = I;
    return res;
  }

  //! Clamps all entries into [Min,Max].
  void clamp ( const DataType Min, const DataType Max ) {
    for ( int i = 0; i < _size; ++i )
      _pData[ i ] = aol::Clamp ( _pData[ i ], Min, Max );
  }

  //! All entries bigger than ThresholdValue are set to B, all others are set to A.
  void threshold ( const DataType ThresholdValue, const DataType A, const DataType B ) {
    for ( int i = 0; i < _size; ++i ) {
      if ( _pData[ i ] > ThresholdValue )
        _pData[ i ] = B;
      else
        _pData[ i ] = A;
    }
  }
  
  //! All entries smaller than ThresholdValue are set to A, all others remain the same.
  void thresholdFromBelow ( const DataType ThresholdValue, const DataType A ) {
    for ( int i = 0; i < _size; ++i ) {
      if ( _pData[ i ] < ThresholdValue )
        _pData[ i ] = A;
    }
  }

  //! Scales all values into [0,1]. The minimal value is mapped to 0, the maximal value is mapped to 1,
  //! intermediate values are handled linearly.
  void scaleValuesTo01 ( ) {
    addToAll ( -1*getMinValue() );
    const DataType maxVal = getMaxValue();
    if ( maxVal != 0 )
      ( *this ) /= maxVal;
  }

  //! Like scaleValuesTo01, but scales all values into [A,B].
  void scaleValuesToAB ( const DataType A, const DataType B ) {
    const DataType min = getMinValue();
    const DataType max = getMaxValue();
    ( *this ) *= ( A - B ) / ( min - max );
    addToAll ( ( B*min - A*max ) / ( min - max ) );
  }

  //! Get data array.
  //! This routine provides direct access to protected data.
  //! Only for compatibility to use in interfaces.
  DataType* getData() const {
    return _pData;
  }

  //! performs \f$ x_i \mapsto \max_j x_j - x_i \f$
  void revert() {
    DataType max = getMaxValue();
    for ( int i = 0; i < _size; ++i ) {
      _pData[ i ] = max - _pData[ i ];
    }
  }

  //! rearranges values in backwards order (in-place)
  void revertOrder ( ) {
    int i0 = 0, i1 = _size - 1;
    DataType val;
    while ( i0 < i1 ) {
      val = this->get ( i0 );
      this->set ( i0, this->get ( i1 ) );
      this->set ( i1, val );
      ++i0;
      --i1;
    }
  }
  
  //! imports values from Other in backwards order, size must be set previously
  void revertOrderFrom( const aol::Vector<DataType> &Other ) {
    for ( int i = 0; i < Other.size(); ++i ) {
      this->set ( i, Other.get( _size - 1 - i ) );
    }
  }

  //! Reads a part from the full vector
  void getBlock ( int start, Vector<DataType>& block ) const;

  //! Write a part into the full vector
  void setBlock ( int start, const Vector<DataType>& block );

  //! Add a part into the full vector
  void addBlock ( int start, const Vector<DataType>& block );

  //! checks whether one of the numbers in this vector is not-a-number or infinite
  //! ATTENTION: the name is misleading!
  bool checkForNANsAndINFs() const {
    for ( int i = 0; i < _size; ++i ) {
      if ( !aol::isFinite ( _pData[i] ) ) return true;
    }
    return false;
  }

  bool compareDim ( const Vector<DataType>& v ) const {
    return ( size() == v.size() );
  }

  //! changes the byte order of the data from little to big endian and vice versa
  //! useful if reading and writing were done on PCs with different CPU types
  void swapByteOrder () {
    for ( int i = 0; i < this->size (); ++i ) {
      DataType org = this->get ( i );
      DataType dest;
      char* op = reinterpret_cast<char*> ( &org );
      char* dp = reinterpret_cast<char*> ( &dest );
      int n = sizeof ( DataType );
      for ( int j = 0; j < n; ++j ) {
        dp [j] = op [n-1-j];
      }
      this->set ( i, dest );
    }
  }

  //! returns the number of entries in this vector that are exactly equal to value
  // todo: write similar method for DataType = float that uses tolerance
  int numOccurence ( const DataType value ) const {
    int count = 0;
    for ( int i = 0; i < this->size(); ++i ) {
      if ( _pData[i] == value ) {
        ++count;
      }
    }
    return ( count );
  }
  
  //! returns the number of entries in this vector that are non-zero
  int numNonZeroes( ) const {
    int nNonZeroes = 0;
    for ( int i = 0; i < this->size(); ++i ) {
      if ( _pData[i] != 0 )
        ++nNonZeroes;
    }
    return nNonZeroes;
  }

  //! elementwise casting
  //! \attention must have same size as argument (call resize/reallocate before convertFrom)
  template <typename InputDataType>
  void convertFrom ( const Vector<InputDataType>& v ) {
    if ( _size != v.size() )
      throw aol::Exception ( "Vector::convertFrom(): Input Vector has wrong size.", __FILE__, __LINE__ );
    for ( int i = 0; i < _size; ++i )
      _pData[i] = static_cast<DataType> ( v[i] );
  }

  //! sorts the entries of this vector (in-place sorting using \c std::make_heap/std::sort_heap )
  void sortValues ( );

  //! remove entry at position Index
  void erase ( const int Index );

  //! insert entry at position Index
  void insert ( const int Index, DataType value );

  //! return index of first occurence of value, -1 if not found
  int indexOfFirstOccurence ( const DataType Value ) const;

  //! erase the first occurence of a given value from the vector; do nothing if value is not found
  void eraseFirstOccurence ( const DataType Value );

  //! return 32-bit checksum of the data contained
  unsigned int crc32OfData ( ) const;

  //! Creates a histogram of the data values with linearly distributed bin boundaries ranging from min entry to max entry
  void createHistogramOfValues ( aol::Vector<int> &Histo, const int NumBins ) const {
    createHistogramOfValues ( Histo, NumBins, getMinValue ( ), getMaxValue ( ) );
  }
  
  //! Creates a histogram of the data values with linearly distributed bin boundaries over a specified range
  //! \author berkels
  void createHistogramOfValues ( aol::Vector<int> &Histo, const int NumBins, const RealType MinVal, const RealType MaxVal ) const {
    Histo.reallocate ( NumBins );

    for ( int i = 0; i < size(); ++i ) {
      const RealType roundedDown = floor ( ( get(i) - MinVal ) * NumBins / ( MaxVal - MinVal ) );
      const int roundedDownInt = static_cast<int> ( roundedDown );
      if ( roundedDownInt >= 0 && roundedDownInt < NumBins )
        ++ ( Histo[ roundedDownInt ] );
      // The values on the right boundary belong to the rightmost bin.
      else if ( ( roundedDownInt == NumBins ) && aol::appeqAbsolute ( roundedDown, static_cast<RealType> ( roundedDownInt ) ) )
        ++ ( Histo[ NumBins-1 ] );
    }
  }
  
  //! Creates a histogram of the data values with linearly distributed bin boundaries ranging from min entry to max entry.
  //! This function stores both the mean value and the number of counts for each bin.
  //! \author mevenkamp
  void createHistogramOfValues ( std::vector<std::pair<RealType, int> > &Histo, const int NumBins ) {
    createHistogramOfValues ( Histo, NumBins, getMinValue ( ), getMaxValue ( ) );
  }
  
  //! Creates a histogram of the data values with linearly distributed bin boundaries over a specified range.
  //! This function stores both the mean value and the number of counts for each bin.
  //! \author mevenkamp
  void createHistogramOfValues ( std::vector<std::pair<RealType, int> > &Histo, const int NumBins, const RealType MinVal, const RealType MaxVal ) const {
    Histo.resize ( NumBins );
    
    aol::Vector<int> histo;
    createHistogramOfValues ( histo, NumBins, MinVal, MaxVal );
    
    for ( int i = 0; i < NumBins ; ++i ) {
      Histo [ i ].first = MinVal + ( i + 0.5 ) / NumBins * ( MaxVal - MinVal );
      Histo [ i ].second = histo[ i ];
    }
  }
  
  void createHistogramOfValues ( aol::Vector<RealType> &BinValues, aol::Vector<int> &BinCounts, const int NumBins ) const {
    createHistogramOfValues ( BinValues, BinCounts, NumBins, getMinValue ( ), getMaxValue ( ) );
  }
  
  void createHistogramOfValues ( aol::Vector<RealType> &BinValues, aol::Vector<int> &BinCounts, const int NumBins, const RealType MinVal, const RealType MaxVal ) const {
    BinValues.reallocate ( NumBins );
    BinCounts.reallocate ( NumBins );
    
    std::vector<std::pair<RealType, int> > histo;
    createHistogramOfValues ( histo, NumBins, MinVal, MaxVal );
    
    for ( int i = 0; i < NumBins ; ++i ) {
      BinValues[ i ] = histo[ i ].first;
      BinCounts[ i ] = histo[ i ].second;
    }
  }

  //! Creates a histogram of the data values (without binning), but is only implemented if DataType is unsigned char.
  //! \author berkels
  void createHistogram ( aol::Vector<int> &Histo) const;
}; // end of Vector


//! operator: vector == Vector
template< typename T >
bool operator== ( const vector<T> &arg1, const Vector<T> &arg2 ) {
  if ( static_cast<int> ( arg1.size() ) != arg2.size() ) {
    return false;
  }
  bool equal = true;
  for ( unsigned int i = 0; i < arg1.size(); ++i ) {
    if ( arg1[i] != arg2[i] ) {
      equal = false;
    }
  }
  return equal;
}


//! operator: Vector == vector
template< typename T >
bool operator== ( const Vector<T> &arg1, const vector<T> &arg2 ) {
  return ( arg2 == arg1 );
}


//! operator: vector != Vector
template< typename T >
bool operator!= ( const vector<T> &arg1, const Vector<T> &arg2 ) {
  return ( ! ( arg1 == arg2 ) );
}


//! operator: Vector != vector
template< typename T >
bool operator!= ( const Vector<T> &arg1, const vector<T> &arg2 ) {
  return ( ! ( arg1 == arg2 ) );
}


//! Read Vector from istream
template <typename DataType>
istream &operator>> ( istream &is, Vector<DataType> &vec ) {
  return vec.read ( is );
}


//! Write Vector to ostream
template <typename DataType>
ostream &operator<< ( ostream &os, const Vector<DataType> &vec ) {
  return vec.print ( os );
}


}

#endif

