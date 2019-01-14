#ifndef __VECTOREXTENSIONS_H
#define __VECTOREXTENSIONS_H

#include <aol.h>
#include <vec.h>

namespace aol {

/** A container class similar to an std::vector that works for classes with standard and
 *  explicit copy constructor (e. g. aol::Vector, aol::MultiVector).
 *  RandomAccessContainer does its own memory management and owns the instances contained.
 *  \attention Similar to std::vector, aol::RandomAccessContainer is not thread safe.
 *  \author Schwen, Wirth
 */
template< typename _DataType >
class RandomAccessContainer {
public:
  typedef _DataType DataType;

protected:
  std::vector< DataType* >  _data;

public:
  RandomAccessContainer ( ) {
  }

  explicit RandomAccessContainer ( const int Size ) {
    _data.reserve( Size );
    for ( int i = 0; i < Size; ++i )
      _data.push_back( new DataType );
  }

  RandomAccessContainer ( const int Size, const DataType &Datum ) {
    _data.reserve( Size );
    for ( int i = 0; i < Size; ++i )
      _data.push_back( new DataType( Datum ) );
  }

#if defined ( USE_CPP11 )
  template<typename... ParamTypes>
  RandomAccessContainer ( const int Size, const ParamTypes&... Parameters ) {
    _data.reserve( Size );
    for ( int i = 0; i < Size; ++i )
      _data.push_back( new DataType ( Parameters... ) );
  }

#else

  //! for initialization by something the contained objects can be conversion-constructed from
  template < typename InitializerType >
  RandomAccessContainer ( const int Size, const InitializerType &Init ) {
    _data.reserve( Size );
    for ( int i = 0; i < Size; ++i )
      _data.push_back( new DataType ( Init ) );
  }

  template<typename ParamType1, typename ParamType2>
  RandomAccessContainer ( const int Size, const ParamType1 &Param1, const ParamType2 &Param2 ) {
    _data.reserve( Size );
    for ( int i = 0; i < Size; ++i )
      _data.push_back( new DataType ( Param1, Param2 ) );
  }

#endif

  explicit RandomAccessContainer ( const RandomAccessContainer<DataType>& other ){
    _data.reserve ( other.size() );
    for ( int i = 0; i < other.size(); ++i )
      _data.push_back( new DataType( other[i] ) );
  }

  RandomAccessContainer& operator= ( const RandomAccessContainer<DataType>& other ) {
    // beware of self-assignment
    if ( this != &other ) {
      for ( unsigned int i = 0; i < _data.size ( ); ++i ) {
        delete _data[i];
      }

      _data.resize( other.size() );

      for ( int i = 0; i < other.size(); ++i ) {
        _data[i] = new DataType( other[i] );
      }
    }
    return ( *this );
  }

  ~RandomAccessContainer ( ) {
    for ( unsigned int i = 0; i < _data.size ( ); ++i ) {
      if ( _data[i] !=  NULL )
        delete _data[i];
    }
  }

  int size ( ) const {
    return ( static_cast<int>( _data.size ( ) ) );
  }

  DataType& operator[] ( const int I ) {
#ifdef DEBUG
    if ( I < 0 || I >= size() ) {
      cerr << "I = " << I << ", size = " << size() << endl;
      throw aol::OutOfBoundsException ("aol::RandomAccessContainer::operator[]: index out of bounds", __FILE__, __LINE__ );
    }
#endif
    return ( *(_data[I]) );
  }

  const DataType& operator[] ( const int I ) const {
#ifdef DEBUG
    if ( I < 0 || I >= size() ) {
      cerr << "I = " << I << ", size = " << size() << endl;
      throw aol::OutOfBoundsException ("aol::RandomAccessContainer::operator[]: index out of bounds", __FILE__, __LINE__ );
    }
#endif
    return ( *(_data[I]) );
  }

  void clear ( ) {
    for ( unsigned int i = 0; i < _data.size ( ); ++i )
      delete _data[i];
    _data.clear();
  }

  void reallocate ( const int NewSize ) {
    this->clear();
    _data.reserve ( NewSize );
    for ( int i = 0; i < NewSize; ++i )
      _data.push_back( new DataType );
  }
  
  void reallocate ( const int NewSize, const DataType &Datum ) {
    this->clear();
    _data.reserve ( NewSize );
    for ( int i = 0; i < NewSize; ++i )
      _data.push_back( new DataType( Datum ) );
  }

  void reallocate ( const RandomAccessContainer<DataType>& other ) {
    this->reallocate ( other.size() );
  }

  void pushBack ( const DataType &Datum ) {
    _data.push_back( new DataType ( Datum ) );
  }

  void pushBack ( const aol::RandomAccessContainer < DataType > &other ) {
    for ( int i = 0; i < other.size (); ++i )
      pushBack ( other[i] );
  }

#if defined ( USE_CPP11 )
  template<typename... ParamTypes>
  void constructDatumAndPushBack ( const ParamTypes&... Parameters ) {
    _data.push_back ( new DataType ( Parameters... ) );
  }

#else

  template<typename ParamType1>
  void constructDatumAndPushBack ( const ParamType1 &Param1) {
    _data.push_back( new DataType ( Param1 ) );
  }
  
  template<typename ParamType1, typename ParamType2>
  void constructDatumAndPushBack ( const ParamType1 &Param1, const ParamType2 &Param2 ) {
    _data.push_back( new DataType ( Param1, Param2 ) );
  }

  template<typename ParamType1, typename ParamType2, typename ParamType3>
  void constructDatumAndPushBack ( const ParamType1 &Param1, const ParamType2 &Param2, const ParamType3 &Param3 ) {
    _data.push_back( new DataType ( Param1, Param2, Param3 ) );
  }

  template<typename ParamType1, typename ParamType2, typename ParamType3, typename ParamType4>
  void constructDatumAndPushBack ( const ParamType1 &Param1, const ParamType2 &Param2, const ParamType3 &Param3, const ParamType4 &Param4 ) {
    _data.push_back( new DataType ( Param1, Param2, Param3, Param4 ) );
  }

  template<typename ParamType1, typename ParamType2, typename ParamType3, typename ParamType4, typename ParamType5>
  void constructDatumAndPushBack ( const ParamType1 &Param1, const ParamType2 &Param2, const ParamType3 &Param3, const ParamType4 &Param4,
      const ParamType5 &Param5 ) {
    _data.push_back( new DataType ( Param1, Param2, Param3, Param4, Param5 ) );
  }

#endif

};

/** A RandomAccessContainer that additionally provides functionalities similar to Vector
 *  (and hence can be used just as a Vector).
 *  \author Wirth
 */
template< typename DataType >
class VectorContainer :
  public RandomAccessContainer<DataType> {

public:
  typedef typename DataType::RealType RealType;
  static const int Depth = DataType::Depth + 1;

public:
  VectorContainer() :
    RandomAccessContainer<DataType>() {}

  explicit VectorContainer( const int Size ) :
    RandomAccessContainer<DataType>( Size ) {}

  VectorContainer( const int Size, const DataType &Datum  ) :
    RandomAccessContainer<DataType>( Size, Datum ) {}

  template < typename InitializerType >
  VectorContainer ( const int Size, const InitializerType &Init ) :
    RandomAccessContainer<DataType>( Size, Init ) {}

  template<typename ParamType1, typename ParamType2>
  VectorContainer ( const int Size, const ParamType1 &Param1, const ParamType2 &Param2 ) :
    RandomAccessContainer<DataType>( Size, Param1, Param2 ) {}

  VectorContainer ( const VectorContainer<DataType>& other , CopyFlag copyFlag = DEEP_COPY) :
    RandomAccessContainer<DataType>() {

    // assume that DataType has a copy-constructor wiht copyFlag argument
    if (copyFlag == DEEP_COPY || copyFlag == STRUCT_COPY) {
      this->_data.reserve ( other.size() );
      for ( int i = 0; i < other.size(); ++i )
        this->_data.push_back ( new DataType ( other[i], copyFlag ) );
    }
    else
      throw aol::Exception ( "aol::VectorContainer<DataType>::VectorContainer( VectorContainer& Data, CopyFlag copyFlag ): illegal copy flag", __FILE__, __LINE__ );
  }

  void getSizes (Vector<int> & sizes) const {
    sizes.resize( this->size() );
    for ( int i = 0; i < this->size(); ++i )
      sizes[i] = ( *this )[i].size ();
  }

  int getTotalSize() const {
    int totalSize = 0;
    for ( int i = 0; i < this->size(); ++i )
      totalSize += (*this)[i].getTotalSize();
    return totalSize;
  }

  //! needed by aol::DerivativeValidatorBase<>
  Vec<Depth, int> setIthComponent ( const int I, const RealType Value ) {
    // bounds checking
    if ( !( I < getTotalSize() ) )
      throw aol::Exception ( "aol::VectorContainer<DataType>::setIthComponent: index out of bounds", __FILE__, __LINE__ );
    // assume that all components of *this have the same total size
    int sizeOfData = (*this)[0].getTotalSize();
    if ( getTotalSize() != ( this->size() * sizeOfData ) )
      throw aol::Exception ( "aol::VectorContainer<DataType>::setIthComponent: all components must have equal sizes!", __FILE__, __LINE__ );
    int comp;
    for( comp = 0; (comp+1)*sizeOfData  <= I; comp++) {}
    int J = I - comp * sizeOfData;
    Vec<DataType::Depth,int> aux = (*this)[comp].setIthComponent ( J, Value );
    Vec<Depth, int> result ( comp );
    for( int c = 0; c < DataType::Depth; c++ )
      result[ c+1 ] = aux[c];
    return result;
  }

  bool checkForNANsAndINFs ( ) const {
    for ( int i = 0; i < this->size(); ++i )
      if( (*this)[i].checkForNANsAndINFs() )
        return true;
    return false;
  }

  void setZero() {
    for ( int i = 0; i < this->size(); ++i )
      (*this)[i].setZero();
  }
    
  void setAll ( const RealType value ) {
    for ( int i = 0; i < this->size(); ++i )
      (*this)[i].setAll ( value );
  }

  VectorContainer<DataType> &  operator+=(  const VectorContainer<DataType>& other ){
    return addMultiple( other, aol::ZOTrait<RealType>::one );
  }

  VectorContainer<DataType> &  operator-=(  const VectorContainer<DataType>& other ){
    return addMultiple( other, -aol::ZOTrait<RealType>::one );
  }

  VectorContainer<DataType> &  operator*=(  const RealType value ){
    for ( int i = 0; i < this->size(); ++i )
      (*this)[i] *= value;
    return *this;
  }

  VectorContainer<DataType> &  addMultiple( const VectorContainer<DataType>& other, const RealType value ){
    if ( this->size() != other.size() )
      throw aol::Exception ( "aol::VectorContainer::addMultiple: sizes don't match", __FILE__, __LINE__ );

    for ( int i = 0; i < this->size(); ++i )
      (*this)[i].addMultiple( other[i], value );
    return *this;
  }
    
  VectorContainer<DataType> & scaleAndAdd ( const RealType factor, const VectorContainer<DataType>& other ){
    if ( this->size() != other.size() )
      throw aol::Exception ( "aol::VectorContainer::addMultiple: sizes don't match", __FILE__, __LINE__ );
   
    for ( int i = 0; i < this->size(); ++i )
      (*this)[i].scaleAndAdd( factor, other[i] );
    return *this;
  }
  
  VectorContainer<DataType> & scaleAndAddMultiple ( const RealType scaleFactor, const VectorContainer<DataType>& other, const RealType vecFactor ){
    if ( this->size() != other.size() )
      throw aol::Exception ( "aol::VectorContainer::addMultiple: sizes don't match", __FILE__, __LINE__ );
    
    for ( int i = 0; i < this->size(); ++i )
      (*this)[i].scaleAndAddMultiple( scaleFactor, other[i], vecFactor );
    return *this;
  }

  RealType normSqr() const {
    RealType s = aol::ZOTrait<RealType>::zero;
    for ( int i = 0; i < this->size(); ++i )
      s += static_cast<RealType>( (*this)[i].normSqr() );
    return s;
  }

  RealType norm() const {
    return sqrt ( normSqr() );
  }


  RealType dotProduct ( const VectorContainer<DataType>& other )  const {
    return this->operator*( other );
  }


  RealType operator* ( const VectorContainer< DataType > &  other ) const {
    if( this->size() != other.size())
      throw Exception ( "aol::VectorContainer<DataType>::operator* : dimensions don't match.", __FILE__, __LINE__ );

    RealType dot = aol::ZOTrait<RealType>::zero;
    for ( int i = 0; i < this->size(); ++i )
      dot += (*this)[i].dotProduct( other[i] ) ;
    return dot;
  }

  RealType getMinValue() const {
    RealType min =  (*this)[0].getMinValue();
    for ( int i = 1; i < this->size(); ++i )
      min = aol::Min ( (*this)[i].getMinValue(), min );
    return min;
  }

  RealType getMaxValue() const {
    RealType max =  (*this)[0].getMaxValue();
    for ( int i = 1; i < this->size(); ++i )
      max = aol::Max ( (*this)[i].getMaxValue(), max );
    return max;
  }
};

/** Abstract multi vector class that models a vector with components of ComponentType
 *  In contrast to aol::MultiVector this also allows for AMultiVectors which contain
 *  AMultiVectors
 *  \author Preusser
 */
template <typename ComponentType>
  class AMultiVector {
 public:
  typedef typename ComponentType::DataType DataType;

  explicit AMultiVector( ) {
  }

  /**
   * Note: This does not allocate any memory for the entries of vecs.
   * The only way to make use of the reserved length of vecs is by calling
   * appendReference.
   */
  void reserve(const int size) {
    vecs.reserve(size);
  }

  void info() const {
    cerr << numComponents() << endl;
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it ) {
      cerr << it->ptr << endl;
      cerr << it->ptr->getMinValue() << "/" << it->ptr->getMaxValue() << endl;
    }
    cerr << "=====================================";
    getchar();
  }

  ~AMultiVector() {
    while ( vecs.size() ) {
      if ( vecs.back().deleteFlag ) {
        delete vecs.back().ptr;
      }
      vecs.pop_back();
    }
  }

  /** Copy constuctor
   */
  explicit AMultiVector ( const AMultiVector &rhs, CopyFlag copyFlag = DEEP_COPY ) {
    vecs.reserve ( rhs.numComponents() );
    for ( int i = 0; i < rhs.numComponents(); ++i ) {
      // call copy constructor of vector with same copyFlag, true is deleteFlag for vec_entries,
      // i. e. the vectors will be deleted and it is up to the vectors to free their memory or not
      vecs.push_back ( vec_entry ( new ComponentType ( rhs[i], copyFlag ), true ) );
    }
  }

  void appendReference(ComponentType &v) {
    vec_entry entry ( &v, false );
    vecs.push_back ( entry );
  }

  ComponentType &operator[] ( int Index ) {
    return * ( vecs[ Index ].ptr );
  }

  AMultiVector<ComponentType> &operator-=( const AMultiVector<ComponentType> &Vec ) {
    int i = 0;
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it, ++i ) {
      * ( it->ptr ) -= Vec[ i ];
    }
    return *this;
  }

  AMultiVector<ComponentType> &operator+=( const AMultiVector<ComponentType> &Vec ) {
    int i = 0;
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it, ++i ) {
      * ( it->ptr ) += Vec[ i ];
    }
    return *this;
  }

  DataType normSqr() const {
    DataType s = 0;
    int i = 0;
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it, ++i ) {
      s += it->ptr->normSqr();
    }
    return s;
  }

  DataType norm() const {
    return sqrt ( normSqr() );
  }

  DataType operator*( const AMultiVector<ComponentType> &Vec ) {
    DataType dot = 0;
    int i = 0;
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it, ++i ) {
      dot += * ( it->ptr ) * Vec[ i ];
    }
    return dot;
  }

  AMultiVector<ComponentType> &operator*= ( DataType Scalar ) {
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it ) {
      * ( it->ptr ) *= Scalar;
    }
    return *this;
  }

  //! Adds multiple of Vec to this MultiVector
  AMultiVector<ComponentType> &addMultiple ( const AMultiVector<ComponentType> &Vec, DataType Factor ) {
    int i = 0;
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it, ++i ) {
      it->ptr->addMultiple ( Vec[ i ], Factor );
    }
    return *this;
  }

  int numComponents() const {
    return vecs.size();
  }

  const ComponentType &operator[] ( int Index ) const {
    return * ( vecs[ Index ].ptr );
  }

  void setZero () {
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it ) {
      if ( !it->ptr ) {
        throw Exception ( "AMultiVector::setZero (): it->ptr == NULL !", __FILE__, __LINE__ );
      }
      it->ptr->setZero ();
    }
  }

  void setAll ( const DataType value ) {
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it ) {
      if ( !it->ptr ) {
        throw Exception ( "AMultiVector::setAll( const RealType value ): it->ptr == NULL !", __FILE__, __LINE__ );
      }
      it->ptr->setAll ( value );
    }
  }

  AMultiVector<ComponentType> &operator= ( const AMultiVector<ComponentType> &Vec ) {
    int i = 0;
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it, ++i ) {
      * ( it->ptr ) = Vec[ i ];
    }

    return *this;
  }

protected:
  struct vec_entry {
    vec_entry ( ComponentType *Ptr, bool DeleteFlag = true ) :
        ptr ( Ptr ), deleteFlag ( DeleteFlag ) {}

    vec_entry() : ptr ( NULL ), deleteFlag ( false ) { }

    ComponentType *ptr;
    bool deleteFlag;
  };

  mutable vector<vec_entry> vecs;
};

/**
 * \brief Collection of helper functions that work on different vector classes.
 *
 * \author Berkels
 */
template <typename VectorType>
class VectorHelper {
  typedef typename VectorType::DataType DataType;
public:
  static DataType getMinValueNoNaNOrInf ( const VectorType &Vector ) {
    DataType min = numeric_limits<DataType>::max();
    
    for ( int i = 0; i < Vector.size(); ++i ) {
      if ( aol::isFinite (  Vector[ i ] ) == false )
        continue;
      
      if ( Vector[ i ] < min ) min = Vector[ i ];
    }
    return min;
  }
  
  static DataType getMaxValueNoNaNOrInf ( const VectorType &Vector ) {
    DataType max = aol::MaxInitializerTrait<DataType>::MaxInitializer;
    
    for ( int i = 0; i < Vector.size(); ++i ) {
      if ( aol::isFinite (  Vector[ i ] ) == false )
        continue;
      
      if ( Vector[ i ] > max ) max = Vector[ i ];
    }
    return max;
  }

  static int numNonNaNs ( const VectorType &Vector ) {
    int num = 0;
    for ( int i = 0; i < Vector.size(); ++i ) {
      if ( aol::isNaN (  Vector[ i ] ) == false )
        ++num;
    }
    return num;
  }

  template <typename InputVectorType>
  static void setBlock ( VectorType &TargetVector, const int StartIdx, const InputVectorType& InputBlock ) {
    const int size = InputBlock.size ();

    if ( StartIdx + size > TargetVector.size() )
      throw aol::Exception ( "aol::VectorHelper::setBlock: Parameter mismatch: Block does not fit", __FILE__, __LINE__ );

    for ( int i = 0; i < size; ++i )
      TargetVector[StartIdx + i] = InputBlock[i];
  }

  static void copyUnblockedFrom ( VectorType &TargetVector, const MultiVector<DataType> &MVector ) {
    if ( TargetVector.size() != MVector.getTotalSize() ) {
      cerr << TargetVector.size() << " " << MVector.getTotalSize() << endl;
      throw aol::Exception ( "aol::VectorHelper::copyUnblockedFrom : dimensions do not match.", __FILE__, __LINE__ );
    }

    TargetVector.setZero();
    const int components = MVector.numComponents();
    int n = 0;
    for ( int i = 0; i < components; ++i ) {
      VectorHelper<VectorType>::setBlock ( TargetVector, n, MVector[i] );
      n += MVector[i].size();
    }
  }

  static void copySplitFrom ( MultiVector<DataType> &TargetMVector, const VectorType &InputVector ) {
    QUOC_ASSERT ( InputVector.size() == TargetMVector.getTotalSize() );
    Vector<int> sizes;
    TargetMVector.getSizes(sizes);
    int index = 0;
    for ( int i = 0; i < sizes.size(); ++i )
      for ( int j = 0; j < sizes[i]; ++j, ++index )
        TargetMVector[i][j] = InputVector[index];
  }
};


}
#endif
