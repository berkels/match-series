#ifndef __MULTIVECTOR_H
#define __MULTIVECTOR_H

#include <aol.h>
#include <vec.h>

namespace aol {

template <typename _DataType>
class MultiVector {
public:
  typedef _DataType DataType;
  static const int Depth = 2;

  typedef typename aol::RealTrait<_DataType>::RealType RealType;

protected:
  struct vec_entry {
    vec_entry ( Vector<DataType> *Ptr, bool DeleteFlag = true ) :
        ptr ( Ptr ), deleteFlag ( DeleteFlag ) {}

    vec_entry() : ptr ( NULL ), deleteFlag ( false ) { }

    Vector<DataType> *ptr;
    bool deleteFlag;
  };

  mutable vector<vec_entry> vecs;

public:
  //! Standard constructor creating an empty MultiVector
  MultiVector ( ) {
  }

  //! Constructor creating MultiVector with Num components of size Size.
  //! Do NOT implement default behavior or a constructor with one int only -- it is unclear
  //! whether the user wants to create one component of size N or N components of size 0.
  MultiVector ( const int Num, const int Size ) {
    vecs.reserve ( Num );
    for ( int i = 0; i < Num; ++i ) {
      vecs.push_back ( vec_entry ( new Vector<DataType> ( Size ) ) );
    }
  }

  //! constructor creating MultiVector from ( Num, Size ) with Num components of size Size
  explicit MultiVector ( const std::pair< int, int > &NSPair ) {
    vecs.reserve ( NSPair.first );
    for ( int i = 0; i < NSPair.first; ++i ) {
      vecs.push_back ( vec_entry ( new Vector<DataType> ( NSPair.second ) ) );
    }
  }

  //! This constructor creates a MultiVector with several size.size() components, component i has size Size[i]
  explicit MultiVector ( const Vector<int>& Size ) {
    vecs.reserve ( Size.size () );
    for ( int i = 0; i < Size.size (); ++i )
      vecs.push_back ( vec_entry ( new Vector<DataType> ( Size [i] ) ) );
  }

  //! This constructor creates a MultiVector with N components, component i has size Size[i]
  template < int N >
  explicit MultiVector ( const Vec<N, int>& Size ) {
    vecs.reserve ( N );
    for ( int i = 0; i < N; ++i )
      vecs.push_back ( vec_entry ( new Vector<DataType> ( Size [i] ) ) );
  }

  explicit MultiVector ( const Vector<unsigned int>& Size ) {
    vecs.reserve ( Size.size () );
    for ( int i = 0; i < Size.size (); ++i ) {
      vecs.push_back ( vec_entry ( new Vector<DataType> ( static_cast<int>( Size [i] ) ) ) );
    }
  }

/*  //! Conversion constructor from AMultiVector non-explicit!!
  MultiVector ( const AMultiVector<Vector<DataType> > &avec) {
    const int nc = avec.numComponents();
    vecs.reserve( nc );
    for(int i=0; i<nc; ++i) {
      vec_entry entry ( &const_cast<Vector<DataType>& > ( avec[i] ), false );
      vecs.push_back ( entry );
    }
  }
*/

  //! Split one long vector to MultiVector with shorter components of specified size, making deep copy.
  MultiVector ( const Vector<DataType>& data, const Vector<int>& size ) {
    int k = 0;
    vecs.reserve ( size.size () );
    for ( int i = 0; i < size.size (); ++i ) {
      vecs.push_back ( vec_entry ( new Vector<DataType> ( size [i] ) ) );
      for ( int j = 0; j < size [i]; ++j, ++k )
        ( *vecs [i].ptr ) [j] = data [k];
    }
    int rest = data.size () - k;
    if ( rest > 0 ) {
      Vector<DataType>* temp = new Vector<DataType> ( rest );
      for ( int i = 0; i < rest; ++i, ++k )
        ( *temp ) [i] = data [k];
      appendReference ( *temp );
    }
  }

  //! Copy constuctor
  //!
  //! for copyFlag DEEP_COPY, create deep copy of every vector contained in
  //! \a rhs, for FLAT_COPY only append references to the vectors that
  //! are contained in \a rhs.
  explicit MultiVector ( const MultiVector &rhs, CopyFlag copyFlag = DEEP_COPY, const int ComponentOffset = 0 ) {
    switch ( copyFlag ) {
    case FLAT_COPY:
      for ( int i = ComponentOffset; i < rhs.numComponents(); ++i ) {
        appendReference ( rhs[i] );
      }
      break;

    case DEEP_COPY:
    case STRUCT_COPY:
    case STRUCT_COPY_UNINIT:
      // jump across labels intended
      vecs.reserve ( rhs.numComponents() );
      for ( int i = ComponentOffset; i < rhs.numComponents(); ++i ) {
        vecs.push_back ( vec_entry ( new Vector<DataType> ( rhs[i], copyFlag ), /* delete this = */ true ) );
      }
      break;

    default:
      throw aol::Exception ( "aol::MultiVector( MultiVector&, copyFlag ): invalid copyFlag specified", __FILE__, __LINE__ );
    }
  }

  //! Constructor creating a MultiVector from qc::GridDefinition, using Dimension components of size Number of Nodes
  template<typename GridType>
  explicit MultiVector ( const GridType &Grid ){
    vecs.reserve ( Grid.getDimOfWorld() );
    for ( int i = 0; i < Grid.getDimOfWorld(); ++i )
      vecs.push_back ( vec_entry ( new Vector<DataType> ( Grid.getNumberOfNodes() ) ) );
  }

  //! Instanciate a MultiVector with one component that is vector wrapper for the field Data of ComponentSize DataTypes.
  MultiVector ( DataType *Data, int ComponentSize, CopyFlag copyFlag = FLAT_COPY ) {
    vecs.push_back ( vec_entry ( new Vector<DataType> ( Data, ComponentSize, copyFlag ), true ) );
  }

  //! Split one long data vector to MultiVector with shorter components of specified size.
  MultiVector ( DataType *Data, const Vector<int>& size, CopyFlag copyFlag = FLAT_COPY ) {
    int k = 0;
    vecs.reserve ( size.size () );
    for ( int i = 0; i < size.size (); ++i ) {
      vecs.push_back ( vec_entry ( new Vector<DataType> ( Data+k, size [i], copyFlag ) ) );
      k += size [i];
    }
  }

  ~MultiVector() {
    while ( vecs.size() ) {
      if ( vecs.back().deleteFlag ) {
        delete vecs.back().ptr;
      }
      vecs.pop_back();
    }
  }

  void copySplitFrom ( const Vector<DataType> & singleVector );

  void copySplitFrom ( const Vector<DataType> & singleVector, const Vector<int> & sizes );

  //! only works if all component vectors have same length!
  void copySplitTransposeFrom ( const Vector<DataType> & singleVector );

  template < int N >
  void copyTo ( aol::Vec<N, DataType> &Dest ) const {
    int k = 0;
    for ( int i = 0; i < numComponents (); ++i ) {
      for ( int j = 0; j < (*this)[i].size (); ++j, ++k )
        Dest [k] = (*this)[i][j];
    }
  }

  ostream& print ( ostream& out ) const {
    for ( int i = 0; i < numComponents(); ++i ) {
      ( *this ) [i].print ( out );
      if ( i < numComponents() - 1 ) out << endl;
    }
    return out;
  }

  //! Saves values into ASCII-file (with precision prec)
  void saveASCII ( const char *fileName, int prec = 10 ) const {
    ofstream file ( fileName );
    file.precision( prec );

    Vector<int> sizes;
    getSizes (sizes);
    for (int i = 0; i < sizes.size(); ++i){
      for (int j = 0; j < sizes[i]; ++j)
        file << (*this)[i][j] <<" ";
      file << endl;
    }

    file.close();
  }

  void load ( const char *FileName );
  void save ( const char *FileName ) const;

  //! if deleteFlag==true, it will be deleted when the multivector is destroyed
  void appendReference ( const Vector<DataType> &Vec, bool deleteFlag = false ) {
    vec_entry entry ( &const_cast<Vector<DataType>& > ( Vec ), deleteFlag );
    vecs.push_back ( entry );
  }

  //! calls appendReference on each component of MultiVec
  void appendReference ( const MultiVector<DataType> &MultiVec, bool deleteFlag = false ) {
    for( int i = 0; i < MultiVec.numComponents(); ++i )
      appendReference( MultiVec[i], deleteFlag );
  }

  void removeReference ( int i ) {
    if ( vecs[i].deleteFlag )
      delete vecs[i].ptr;

    vecs.erase ( vecs.begin () + i );
  }

  //! check if *this owns data
  bool isDataOwned ( ) const {
    for ( int i = 0; i < static_cast<int>( vecs.size() ); ++i )
      if( ( !vecs[i].deleteFlag ) && ( (*this)[i].getData() != NULL ) )
        return false;
    return true;
  }

  //! change size, deleting old contents
  void reallocate ( const int NumComponents, const int SizeOfComponents ) {
    resize ( 0, 0 );
    vecs.reserve ( NumComponents );
    for ( int i = 0; i < NumComponents; ++i ) {
      vecs.push_back ( vec_entry ( new Vector<DataType> ( SizeOfComponents ) ) ) ;
    }
  }

  //! change size, preserving old contents to the extent possible
  void resize ( const int NumComponents, const int SizeOfComponents ) {
    if ( !isDataOwned() )
      throw Exception ( "aol::MultiVector<DataType>::resize may not be called on FLAT_COPY vectors!\n", __FILE__, __LINE__ );

    if ( vecs.size() > static_cast<unsigned int>( NumComponents ) ) {
      for ( int i = static_cast<int>( vecs.size() ); i > NumComponents; --i ) {
        if ( vecs.back().deleteFlag ) {
          delete vecs.back().ptr;
        }
        vecs.pop_back();
      }
    }

    vecs.reserve ( NumComponents );
    if ( static_cast<int>( vecs.size() ) < NumComponents ) {
      for ( int i = static_cast<int>( vecs.size() ); i < NumComponents; ++i) {
        vecs.push_back ( vec_entry ( new Vector<DataType> ( SizeOfComponents ) ) ) ;
      }
    }

    for ( int i = 0; i < NumComponents; ++i ) {
      (vecs[i].ptr)->resize ( SizeOfComponents );
    }
  }

  //! change size, preserving old contents to the extent possible
  void resize ( const MultiVector<_DataType> &Other ) {
    if ( !isDataOwned() )
      throw Exception ( "aol::MultiVector<DataType>::resize may not be called on FLAT_COPY vectors!\n", __FILE__, __LINE__ );

    if ( vecs.size() > static_cast<unsigned int>( Other.numComponents() ) ) {
      for ( int i = static_cast<int>( vecs.size() ); i > Other.numComponents(); --i ) {
        if ( vecs.back().deleteFlag ) {
          delete vecs.back().ptr;
        }
        vecs.pop_back();
      }
    }

    vecs.reserve ( Other.numComponents() );
    if ( static_cast<int>( vecs.size() ) < Other.numComponents() ) {
      for ( int i = static_cast<int>( vecs.size() ); i < Other.numComponents(); ++i) {
        vecs.push_back ( vec_entry ( new Vector<DataType> ( Other[i].size() ) ) ) ;
      }
    }

    for ( int i = 0; i < Other.numComponents(); ++i ) {
      (vecs[i].ptr)->resize ( Other[i].size() );
    }
  }

  //! clear the vector
  void clear () {
    for ( unsigned int i = 0; i < vecs.size (); ++i ) {
      if ( vecs[i].deleteFlag ) {
        delete vecs[i].ptr;
      }
    }

    vecs.clear ();
  }

  void reserve ( const int NumComponents, const int SizeOfComponents ) {
    if ( !isDataOwned() )
      throw Exception ( "aol::MultiVector<DataType>::reserve may not be called on FLAT_COPY vectors!\n", __FILE__, __LINE__ );

    if ( static_cast<int>( vecs.size() ) != NumComponents )
      throw aol::Exception ( "aol::MultiVector::reserve cannot change the number of components.", __FILE__, __LINE__ );

    for ( int i = 0; i < static_cast<int>( vecs.size() ); ++i ) {
      (vecs[i].ptr)->reserve ( SizeOfComponents );
    }

  }


  //! change size of MultiVector and its components to values given by Vector<int> size. Old contents are deleted.
  void reallocate ( const Vector<int> &size ) {
    if ( !isDataOwned() )
      throw Exception ( "aol::MultiVector<DataType>::reallocate may not be called on FLAT_COPY vectors!\n", __FILE__, __LINE__ );

    while ( vecs.size() ) {
      if ( vecs.back().deleteFlag ) {
        delete vecs.back().ptr;
      }
      vecs.pop_back();
    }
    // A simple vecs.clear(); leads to a memory leak.

    for ( int i = 0; i < size.size (); ++i ) {
      vecs.push_back ( vec_entry ( new Vector<DataType> ( size [i] ) ) ) ;
    }
  }

  void reallocate ( const MultiVector<DataType> &other ) {
    Vector<int> sizes ( other.numComponents() );
    for ( int i = 0; i < sizes.size(); ++i )
      sizes[i] = other[i].size();
    this->reallocate ( sizes );
  }

  void reallocate ( const qc::GridStructure &grid );

  Vector<DataType> &operator[] ( int Index ) {
#ifdef BOUNDS_CHECK
    if ( Index < 0 || Index >= static_cast<int>(vecs.size()) ) {
      char error[1024];
      sprintf ( error, "MultiVector<Realtype>::operator[]: index %d out of bounds. N = %d\n", Index, static_cast<int>(vecs.size()) );
      throw OutOfBoundsException ( error, __FILE__, __LINE__ );
    }
#endif
    return * ( vecs[ Index ].ptr );
  }

  const Vector<DataType> &operator[] ( int Index ) const {
#ifdef BOUNDS_CHECK
    if ( Index < 0 || Index >= static_cast<int>(vecs.size()) ) {
      char error[1024];
      sprintf ( error, "MultiVector<Realtype>::operator[]: index %d out of bounds. N = %d\n", Index, static_cast<int>(vecs.size()) );
      throw OutOfBoundsException ( error, __FILE__, __LINE__ );
    }
#endif
    return * ( vecs[ Index ].ptr );
  }

  bool compareDim ( const MultiVector<DataType>& mv ) const {
    if ( mv.numComponents () != numComponents () ) return false;
    for ( int i = 0; i < numComponents (); ++i )
      if ( mv [i].size () != ( *this ) [i].size () ) return false;
    return true;
  }

  bool allDimsEqual () const {
    int dim = ( *this ) [0].size ();
    for ( int i = 1; i < numComponents (); ++i )
      if ( ( *this ) [i].size () != dim ) return false;
    return true;
  }

  bool allDimsEqual ( const int dim ) const {
    for ( int i = 0; i < numComponents (); ++i )
      if ( ( *this ) [i].size () != dim ) return false;
    return true;
  }

  //! return size of one component if all components have the same length
  int getEqualComponentSize ( ) const {
    if ( allDimsEqual() == false ) {
      Vector<int> sizes;
      getSizes ( sizes );
      cerr << sizes << endl;
      throw ( aol::Exception ( "Sizes of components differ.", __FILE__, __LINE__ ) );
      return ( -1 );
    }
    return ( ( *this ) [0].size () );
  }

  void getSizes (Vector<int> & sizes) const {
    sizes.resize( numComponents() );
    for ( int i = 0; i < numComponents (); ++i )
      sizes[i] = ( *this ) [i].size ();
  }

  int getTotalSize () const {
    int size = 0;
    for ( int i = 0; i < numComponents (); ++i ) {
      size += ( *this ) [i].size ();
    }
    return size;
  }

  int numComponents() const {
    return static_cast<int>(vecs.size());
  }

  void eraseComponent ( const int Index ) {
    if ( vecs[Index].deleteFlag ) {
      delete vecs[Index].ptr;
    }
    vecs.erase (vecs.begin()+Index);
  }

  MultiVector<DataType> &operator+= ( const MultiVector<DataType> &Vec ) {
    if ( !Vec.compareDim ( *this ) ) {
      throw Exception ( "MultiVector<DataType>::operator+= : dimensions don't match.", __FILE__, __LINE__ );
    }
    int i = 0;
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it, ++i ) {
      * ( it->ptr ) += Vec[ i ];
    }

    return *this;
  }

  MultiVector<DataType> &operator*= ( const MultiVector<DataType> &Vec ) {
    if ( !Vec.compareDim ( *this ) ) {
      throw Exception ( "MultiVector<DataType>::operator*= : dimensions don't match.", __FILE__, __LINE__ );
    }
    int i = 0;
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it, ++i ) {
      * ( it->ptr ) *= Vec[ i ];
    }

    return *this;
  }

  MultiVector<DataType> &operator/= ( const MultiVector<DataType> &Vec ) {
    if ( !Vec.compareDim ( *this ) ) {
      throw Exception ( "MultiVector<DataType>::operator/= : dimensions don't match.", __FILE__, __LINE__ );
    }
    int i = 0;
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it, ++i ) {
      * ( it->ptr ) /= Vec[ i ];
    }

    return *this;
  }

  bool operator== ( const MultiVector<DataType> &Vec ) const {
    if ( !Vec.compareDim ( *this ) ) {
      throw Exception ( "MultiVector<DataType>::operator== : dimensions don't match.", __FILE__, __LINE__ );
    }
    int i = 0;
    for ( typename vector<vec_entry>::const_iterator it = vecs.begin(); it != vecs.end(); ++it, ++i ) {
      if ( * ( it->ptr ) != Vec[ i ] ) {
        return false;
      }
    }
    return true;
  }

  bool operator!= ( const MultiVector<DataType> &Vec ) const {
    if ( !Vec.compareDim ( *this ) ) {
      throw Exception ( "MultiVector<DataType>::operator!= : dimensions don't match.", __FILE__, __LINE__ );
    }
    int i = 0;
    for ( typename vector<vec_entry>::const_iterator it = vecs.begin(); it != vecs.end(); ++it, ++i ) {
      if ( * ( it->ptr ) == Vec[ i ] ) {
        return false;
      }
    }
    return true;
  }

  MultiVector<DataType> &operator-= ( const MultiVector<DataType> &Vec ) {
    if ( !Vec.compareDim ( *this ) ) {
      throw Exception ( "MultiVector<DataType>::operator-= : dimensions don't match.", __FILE__, __LINE__ );
    }
    int i = 0;
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it, ++i ) {
      * ( it->ptr ) -= Vec[ i ];
    }
    return *this;
  }

  //! Adds multiple of Vec to this MultiVector
  MultiVector<DataType> &addMultiple ( const MultiVector<DataType> &Vec, DataType Factor ) {
    if ( !Vec.compareDim ( *this ) ) {
      throw Exception ( "MultiVector<DataType>::addMultiple dimensions don't match.", __FILE__, __LINE__ );
    }
    int i = 0;
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it, ++i ) {
      it->ptr->addMultiple ( Vec[ i ], Factor );
    }
    return *this;
  }

  MultiVector<DataType> &scaleAndAdd ( const DataType Factor, const aol::MultiVector<DataType>& Vec ) {
    if ( !Vec.compareDim ( *this ) )
      throw Exception ( "MultiVector<DataType>::scaleAndAdd dimensions don't match.", __FILE__, __LINE__ );

    for ( int i = 0; i < numComponents (); ++i )
      ( *this ) [i].scaleAndAdd ( Factor, Vec[i] );

    return *this;
  }

  MultiVector<DataType> & scaleAndAddMultiple ( const DataType ScaleFactor, const aol::MultiVector<DataType>& Vec, const DataType VecFactor ) {
    if ( !Vec.compareDim ( *this ) )
      throw Exception ( "MultiVector<DataType>::scaleAndAddMultiple dimensions don't match.", __FILE__, __LINE__ );

    for ( int i = 0; i < numComponents (); ++i )
      ( *this ) [i].scaleAndAddMultiple ( ScaleFactor, Vec[i], VecFactor );

    return *this;
  }

  void setSum ( const MultiVector<DataType>& Vec1, const MultiVector<DataType>& Vec2, DataType Factor ) {
    if ( !Vec1.compareDim ( *this ) || !Vec2.compareDim ( *this ) ) {
      throw Exception ( "MultiVector<DataType>::setSum dimensions don't match.", __FILE__, __LINE__ );
    }
    for ( int i = 0; i < numComponents (); ++i )
      ( *this ) [i].setSum ( Vec1[ i ], Vec2[ i ], Factor );
  }

  //! Adds multiple of Vec to this MultiVector
  MultiVector<DataType> &addMultipleParallel ( const MultiVector<DataType> &Vec, DataType Factor ) {
    if ( !Vec.compareDim ( *this ) ) {
      throw Exception ( "MultiVector<DataType>::addMultipleParallel dimensions don't match.", __FILE__, __LINE__ );
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < this->numComponents(); ++i ) {
      this->operator[](i).addMultiple ( Vec[ i ], Factor );
    }
    return *this;
  }

  //! Adds scalar value to all entries
  MultiVector<DataType> &addToAll ( const DataType Scalar ) {
    for ( int i = 0; i < this->numComponents(); ++i )
      this->operator[]( i ).addToAll( Scalar );
    return *this;
  }

  //! Dot product
  DataType operator* ( const MultiVector<DataType> &Vec ) const {
    if ( !Vec.compareDim ( *this ) ) {
      cerr << this->numComponents() << " " << Vec.numComponents() << endl;
      if ( this->numComponents() == Vec.numComponents() )
        for ( int i = 0; i < this->numComponents(); ++i )
          cerr << (*this)[i].size() << " " << Vec[i].size() << "    ";
      cerr << endl;

      throw Exception ( "MultiVector<DataType>::operator* : dimensions don't match.", __FILE__, __LINE__ );
    }
    DataType dot = 0;
    int i = 0;
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it, ++i ) {
      dot += * ( it->ptr ) * Vec[ i ];
    }
    return dot;
  }

  //! Dot product
  DataType dotProduct( const MultiVector<DataType> &Vec ) const {
    return this->operator*( Vec );
  }

  DataType normSqr() const {
    DataType s = 0;
    int i = 0;
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it, ++i ) {
      s += it->ptr->normSqr();
    }
    return s;
  }

  RealType norm() const {
    return sqrt ( static_cast<RealType> ( normSqr() ) );
  }

  /** Returns lp norm of MultiVector interpreted as Vector
    * \author Preusser
    */
  RealType lpNorm ( RealType p ) const;


  DataType getMaxValue() const {
    typename vector<vec_entry>::iterator it = vecs.begin();
    DataType max = it++->ptr->getMaxValue();
    for ( ; it != vecs.end(); ++it ) {
      max = aol::Max ( it->ptr->getMaxValue(), max );
    }
    return max;
  }

  DataType getMaxAbsValue() const {
    typename vector<vec_entry>::iterator it = vecs.begin();
    DataType maxabs = it++->ptr->getMaxAbsValue();
    for ( ; it != vecs.end(); ++it ) {
      maxabs = aol::Max ( it->ptr->getMaxAbsValue(), maxabs );
    }
    return maxabs;
  }

  DataType getMinValue() const {
    typename vector<vec_entry>::iterator it = vecs.begin();
    DataType min = it++->ptr->getMinValue();
    for ( ; it != vecs.end(); ++it ) {
      min = aol::Min ( it->ptr->getMinValue(), min );
    }
    return min;
  }
  
  void getMeanComponents ( aol::Vector<DataType> &Dest ) const {
    if ( Dest.size ( ) != this->numComponents ( ) ) throw aol::Exception ( "Dimension of destination vector does not match number of components!", __FILE__, __LINE__ );
    for ( int i = 0; i < numComponents() ; ++i ) Dest[i] = (*this)[i].getMeanValue ( );
  }
  
  //! Returns the mininum of the min values for each component and the maximum of the max values for each component,
  //! such that half of the specified percentage of the entries of the corresponding component are lower than min and half bigger than max.
  //! For instance useful to enhance the contrast. As percentage the expected range for the argument is 0 to 100.
  Vec2<DataType> getSaturatedMinMaxValue( const RealType SaturatedEntryPercentage ) const {
    Vector<DataType> minVals ( this->numComponents ( ) ), maxVals ( this->numComponents ( ) );
    Vec2<DataType> minMax;
    for ( int i = 0; i < numComponents(); ++i ) {
      minMax = (*this)[i].getSaturatedMinMaxValue ( SaturatedEntryPercentage );
      minVals[i] = minMax[0];
      maxVals[i] = minMax[1];
    }
    minMax[0] = minVals.getMinValue ( );
    minMax[1] = maxVals.getMaxValue ( );
    return minMax;
  }

  bool checkForNANsAndINFs() const {
    bool nans = false;
    typename vector<vec_entry>::iterator it = vecs.begin();
    for ( ; it != vecs.end(); ++it ) {
      if ( it->ptr->checkForNANsAndINFs() )
        nans = true;
    }
    return nans;
  }

  MultiVector<DataType> &operator= ( const MultiVector<DataType> &Vec ) {
    if ( !Vec.compareDim ( *this ) ) {
      throw Exception ( "MultiVector<DataType>::operator= : dimensions don't match.", __FILE__, __LINE__ );
    }
    int i = 0;
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it, ++i ) {
      * ( it->ptr ) = Vec[ i ];
    }

    return *this;
  }

  MultiVector<DataType> &operator*= ( const DataType Scalar ) {
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it ) {
      * ( it->ptr ) *= Scalar;
    }
    return *this;
  }

  MultiVector<DataType> &operator/= ( const DataType Scalar ) {
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it ) {
      * ( it->ptr ) /= Scalar;
    }
    return *this;
  }

  void setZero () {
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it ) {
      if ( !it->ptr ) {
        throw Exception ( "MultiVector::setZero (): it->ptr == NULL !", __FILE__, __LINE__ );
      }
      it->ptr->setZero ();
    }
  }

  //! Set all components of all vectors in MultiVector to same value
  void setAll ( const DataType value ) {
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it ) {
      if ( !it->ptr ) {
        throw Exception ( "MultiVector::setAll( const DataType value ): it->ptr == NULL !", __FILE__, __LINE__ );
      }
      it->ptr->setAll ( value );
    }
  }

  //! For each component: set all components that are marked as "true" (if invertMask == false)
  void setAllMasked ( const DataType val, const aol::BitVector& mask, bool invertMask = false ) {
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it ) {
      if ( !it->ptr )
        throw Exception ( "MultiVector::setAllMasked: it->ptr == NULL !", __FILE__, __LINE__ );
      if ( mask.size() !=  it->ptr->size() )
        throw Exception ( "MultiVector::setAllMasked: sizes don`t match !", __FILE__, __LINE__ );
      it->ptr->setAllMasked( val, mask, invertMask );
    }
  }

  //! get values at position I in all components, avoiding temporary object
  template < int N >
  void getTo ( const int I, aol::Vec<N, DataType> &vals ) const {
    if ( N != numComponents() )
      throw aol::Exception ( "aol::MultiVector<DataType>::getTo ( int, aol::Vec ): size of Vec does not match", __FILE__, __LINE__ );
    for ( int j = 0; j < N; ++j )
      vals[j] = (*this)[j][I];
  }

  //! get values at position I in all components, avoiding temporary object
  void getTo ( const int I, aol::Vector<DataType> &vals ) const {
    const int vectorSize = vals.size();
    if ( vectorSize != numComponents() )
      throw aol::Exception ( "aol::MultiVector<DataType>::getTo ( int, aol::Vector ): size of vector does not match", __FILE__, __LINE__ );
    for ( int j = 0; j < vectorSize; ++j )
      vals[j] = (*this)[j][I];
  }

  //! set values at position I in all components
  template < int N >
  void set ( const int I, const aol::Vec<N, DataType> &vals ) {
    if ( N != numComponents() )
      throw aol::Exception ( "aol::MultiVector<DataType>::set ( int, aol::Vec ): size of Vec does not match", __FILE__, __LINE__ );
    for ( int j = 0; j < N; ++j )
      (*this)[j][I] = vals[j];
  }

  //! set values at position I in all components
  void set ( const int I, const aol::Vector<DataType> &vals ) {
    const int vectorSize = vals.size ();
    if ( vectorSize != numComponents() )
      throw aol::Exception ( "aol::MultiVector<DataType>::set ( int, aol::Vector ): size of vector does not match", __FILE__, __LINE__ );
    for ( int j = 0; j < vectorSize; ++j )
      (*this)[j][I] = vals[j];
  }
  
  //! add values at position I in all components
  template < int N >
  void addMultiple ( const int I, const aol::Vec<N, DataType> &vals, DataType Factor = 1. ) {
    if ( N != numComponents() )
      throw aol::Exception ( "aol::MultiVector<DataType>::addMultiple ( int, aol::Vec, DataType ): size of Vec does not match", __FILE__, __LINE__ );
    for ( int j = 0; j < N; ++j )
      (*this)[j][I] += Factor * vals[j];
  }

  //! set value at position I, 0 <= I <= getTotalSize()-1  ( needed by aol::DerivativeValidatorBase<> )
  Vec2<int> setIthComponent ( const int I, const DataType Value ) {
    int i, j;
    for (i = 0, j = I; j >= (*this)[i].size(); j -= (*this)[i++].size()) ;
    (*this)[i][j] = Value;
    Vec2<int> res(i,j);
    return res;
  }

  //! Clamps all entries of all vectors in MultiVector into [Min,Max].
  void clamp ( const DataType Min, const DataType Max ) {
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it ) {
      if ( !it->ptr ) {
        throw Exception ( "MultiVector::clamp ( const DataType Min, const DataType Max ): it->ptr == NULL !", __FILE__, __LINE__ );
      }
      it->ptr->clamp( Min, Max );
    }
  }

  //! All entries of all vectors bigger then ThresholdValue are set to B, all others are set to A.
  void threshold ( const DataType ThresholdValue, const DataType A, const DataType B ) {
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it ) {
      if ( !it->ptr ) {
        throw Exception ( "MultiVector::threshold ( const DataType ThresholdValue, const DataType A, const DataType B ): it->ptr == NULL !", __FILE__, __LINE__ );
      }
      it->ptr->threshold( ThresholdValue, A, B );
    }
  }
  
   /**
   * \brief All entries of all vectors smaller then ThresholdValue are set to A, all others remain the same.
   *
   * \author Simon
   */
  void thresholdFromBelow ( const DataType ThresholdValue, const DataType A ) {
    for ( typename vector<vec_entry>::iterator it = vecs.begin(); it != vecs.end(); ++it ) {
      if ( !it->ptr ) {
        throw Exception ( "MultiVector::thresholdFromBelow ( const DataType ThresholdValue, const DataType A ): it->ptr == NULL !", __FILE__, __LINE__ );
      }
      it->ptr->thresholdFromBelow( ThresholdValue, A );
    }
  }

  //! Scales all values into [0,1]. The minimal value is mapped to 0, the maximal value is mapped to 1,
  //! intermediate values are handled linearly.
  void scaleValuesTo01 ( ) {
    addToAll ( -1*getMinValue() );
    ( *this ) /= getMaxValue();
  }

  /**
   * \brief Calculates the postion-wise median over the components (assuming all components have the same size), ignoring inf values.
   *
   * \author Berkels
   */
  void getMedianVecOverComponents ( aol::Vector<RealType> &MedianVec, const RealType FillValue = 0 );

  void setAll ( const aol::Vector<DataType> &vecValue ) {
    for ( int i = 0; i < numComponents(); ++i )
      (*this)[i] = vecValue;
      //*( vecs[ i ].ptr ) = vecValue;
  }
  void setVectorRef ( int Index, Vector<DataType> &Vec, bool deleteFlag = false ) {
    if ( Index >= static_cast<int> ( vecs.size() ) ) {
      throw Exception ( "MultiVector::setVectorRef: index out of bounds.\n", __FILE__, __LINE__ );
    }
    if ( vecs[ Index ].ptr && vecs[ Index ].deleteFlag ) {
      delete vecs[ Index ].ptr;
    }
    vecs[ Index ].ptr = &Vec;
    vecs[ Index ].deleteFlag = deleteFlag;
  }

  //! lossless save in binary quocmesh format
  void saveToFile ( const char *filename ) const;

  //! load from file in quocmesh format, no conversion of data types
  void loadFromFile ( const char *filename );

  //! load from file in quocmesh format with different DataType
  template< typename LoadDataType >
  void loadConvertFromFile ( const char *filename ) {
    aol::MultiVector<LoadDataType> loadVec;
    loadVec.loadFromFile ( filename );
    Vector<int> sizes;
    loadVec.getSizes (sizes);
    this->reallocate ( sizes );
    for ( int i = 0; i < this->numComponents(); ++i ) {
      (*this)[i].convertFrom ( loadVec[i] );
    }
  }
  
  //! load from ASCII file (expecting spaces as separators)
  void loadASCII ( const char* filename ) {
    std::vector<std::vector<std::string> > cells;
    std::ifstream asciiFile ( filename );
    std::string line, cell;
    if ( asciiFile.is_open ( ) ) {
      while ( std::getline ( asciiFile, line ) ) {
        if ( line.length ( ) > 0 && line[0] != '#' ) {
          cells.push_back ( std::vector<string> ( ) );
          std::stringstream lineStream ( line );
            while ( std::getline ( lineStream, cell, ' ' ) )
              cells[cells.size ( ) -1].push_back ( cell );
        }
          
      }
      asciiFile.close ( );
      
      Vector<int> sizes ( static_cast<int> ( cells.size ( ) ) );
      for ( unsigned int i=0; i<cells.size ( ) ; ++i )
        sizes[i] = static_cast<int> ( cells[i].size ( ) );
      this->reallocate ( sizes );
      
      for ( unsigned int i=0; i<cells.size ( ) ; ++i ) {
        for ( unsigned int j=0; j<cells[i].size ( ) ; ++j ) {
          std::istringstream istr ( cells[i][j] );
          DataType number;
          istr >> number;
          (*this)[i][j] = number;
        }
      }
    } else
      throw aol::Exception ( "Unable to open ASCII file!", __FILE__, __LINE__ );
  }
};

//! Read MultiVector from istream
template <typename DataType>
istream &operator>> ( istream &is, MultiVector<DataType> &MVec ) {
  return MVec.read ( is );
}

//! Write MultiVector to ostream
template <typename DataType>
ostream &operator<< ( ostream &os, const MultiVector<DataType> &MVec ) {
  return MVec.print ( os );
}

/**
 * \author Berkels
 */
template <typename VecType>
class VectorInitTrait {};

template <typename RealType>
class VectorInitTrait<aol::Vector<RealType> > {
public:
  typedef int InitType;
};

template <typename RealType>
class VectorInitTrait<aol::MultiVector<RealType> > {
public:
  typedef aol::Vector<int> InitType;
};

// ___________________________________________________________________________________
//
// ***********************************************************************************
// ===================================================================================
//

/** Data array that performs exactly as MultiVector, only that
 *  operator[] returns a (const) reference to a subclass of Vector.
 */
template <typename DerivedType>
class MultiDerivedVector : public MultiVector<typename DerivedType::DataType> {

public:
  // *** constructors ***
  MultiDerivedVector ( const int Num, const int Size )
    : MultiVector<typename DerivedType::DataType> ( Num, Size )
  {}

  explicit MultiDerivedVector ( Vector<int>& Size )
    : MultiVector<typename DerivedType::DataType> ( Size )
  {}

  MultiDerivedVector ( const MultiDerivedVector<DerivedType> &rhs, CopyFlag copyFlag = DEEP_COPY )
    : MultiVector<typename DerivedType::DataType> ( rhs, copyFlag )
  {}

  explicit MultiDerivedVector ( const qc::GridStructure &Grid )
    : MultiVector<typename DerivedType::DataType> ( Grid )
  {}


  DerivedType & operator[] ( int Index ) {
    return dynamic_cast<DerivedType &>
              ( MultiVector<typename DerivedType::DataType>::operator[] ( Index ) );
  }

  const DerivedType & operator[] ( int Index ) const {
    return dynamic_cast<const DerivedType &>
              ( MultiVector<typename DerivedType::DataType>::operator[] ( Index ) );
  }
};

} // end of namespace aol.

#endif
