#ifndef __MULTIARRAY_H
#define __MULTIARRAY_H

#include <multiVector.h>
#include <multilevelArray.h>
#include <scalarArray.h>
#include <smallMat.h>
#include <rgbColorMap.h>

namespace qc {

//! Very crude version of an vector-valued function on a cartesian grid, abstract in both dimensions
//! Still much work and documentation to do
template <class DataType, int rangedim, int imagedim = rangedim>
class MultiArray : public aol::MultiVector<DataType> {
public:

  typedef typename ScalarArrayTrait<DataType, rangedim>::ArrayType ScalarArrayType;
  typedef typename aol::VecDimTrait<DataType, imagedim>::VecType DataVecType;
  typedef typename aol::VecDimTrait<short, rangedim>::VecType IndexVecType;

  enum create_arrays { CREATE_ALL_ARRAYS,
                       CREATE_NO_ARRAYS   };

  //! standard constructor calling standard constructor for all ScalarArrays contained
  MultiArray ( ) : aol::MultiVector<DataType> ( 0, 0 ) {
    for ( int i = 0; i < imagedim; ++i ) {
      this->appendReference ( * ( new ScalarArrayType ), true );
    }
  }

  //! by default, new ScalarArrays are created. If you want to add references, you must make sure to add the correct number of references before using this class.
  MultiArray ( int sizeX, int sizeY, create_arrays ca = CREATE_ALL_ARRAYS ) : aol::MultiVector<DataType> ( 0, 0 ) {
    if ( ca == CREATE_ALL_ARRAYS ) {
      for ( int i = 0; i < imagedim; ++i ) {
        this->appendReference ( * ( new ScalarArrayType ( sizeX, sizeY ) ), true );
      }
    }
  }

  explicit MultiArray ( const qc::GridStructure &grid, create_arrays ca = CREATE_ALL_ARRAYS ) : aol::MultiVector<DataType> ( 0, 0 ) {
    if ( ca == CREATE_ALL_ARRAYS ) {
      for ( int i = 0; i < imagedim; ++i ) {
        this->appendReference ( * ( new ScalarArrayType ( grid ) ), true );
      }
    }
  }

  template <qc::Dimension Dim>
  explicit MultiArray ( const qc::GridSize<Dim> &size, create_arrays ca = CREATE_ALL_ARRAYS ) : aol::MultiVector<DataType> ( 0, 0 ) {
    if ( ca == CREATE_ALL_ARRAYS ) {
      for ( int i = 0; i < imagedim; ++i ) {
        this->appendReference ( * ( new ScalarArrayType ( size ) ), true );
      }
    }
  }

  explicit MultiArray ( const qc::MultiArray<DataType, rangedim, imagedim> &other, aol::CopyFlag copytype = aol::DEEP_COPY ) : aol::MultiVector<DataType> ( 0, 0 ) {
    // Calling the copy constructor of MultiVector with other and copytype as arguments does NOT work
    // since this creates new pointers to aol::Vector, which can't be casted into ScalarArray pointers!
    this->vecs.reserve ( imagedim );
    for ( int i = 0; i < imagedim; i++ ) {
      this->vecs.push_back ( typename aol::MultiVector<DataType>::vec_entry ( new ScalarArrayType ( other[i], copytype ), true ) ); // call copy constructor of ScalarArray with same copytype, true is deleteFlag for vec_entries, i. e. the ScalarArrays will be deleted and it is up to the ScalarArrays to free their memory or not
    }
  }

  //! \todo Get rid of this constructor, the order of the first three arguments if different from the corresponding ScalarArray constructor.
  template <typename GridType>
  MultiArray ( const GridType &grid, const aol::MultiVector<DataType> &multiVector, aol::CopyFlag CopyType = aol::FLAT_COPY, const int ComponentOffset = 0 ) : aol::MultiVector<DataType> ( 0, 0 ) {
    if( ( multiVector.numComponents() + ComponentOffset ) < imagedim )
      throw aol::Exception( "multiVector.numComponents() < imagedim !", __FILE__, __LINE__);
    for ( int i = 0; i < imagedim; ++i ) {
      this->appendReference ( * ( new ScalarArrayType ( multiVector[i+ComponentOffset], grid, CopyType ) ), true );
    }
  }

  template <typename GridType>
  MultiArray ( const aol::MultiVector<DataType> &multiVector, const GridType &grid, aol::CopyFlag CopyType = aol::FLAT_COPY, const int ComponentOffset = 0 ) : aol::MultiVector<DataType> ( 0, 0 ) {
    if( ( multiVector.numComponents() + ComponentOffset ) < imagedim )
      throw aol::Exception( "multiVector.numComponents() < imagedim !", __FILE__, __LINE__);
    for ( int i = 0; i < imagedim; ++i ) {
      this->appendReference ( * ( new ScalarArrayType ( multiVector[i+ComponentOffset], grid, CopyType ) ), true );
    }
  }

  // Note: RealType is an extra template argument here, because we don't want this constructor to be compiled automatically for "unsigned char"
  template <typename RealType, typename ArrayType, typename ProlongOpType, typename RestrictOpType, typename GridType>
  MultiArray ( const qc::MultiDimMultilevelArray<RealType, ArrayType, ProlongOpType, RestrictOpType, GridType> &MDMLArray, aol::CopyFlag CopyType = aol::FLAT_COPY ) : aol::MultiVector<DataType> ( 0, 0 ) {
    if( MDMLArray.numComponents() < imagedim )
      throw aol::Exception( "MDMLArray.numComponents() < imagedim !", __FILE__, __LINE__);
    for ( int i = 0; i < imagedim; ++i ) {
      this->appendReference ( * ( new ScalarArrayType ( MDMLArray.getArray ( i ), CopyType ) ), true );
    }
  }

  explicit MultiArray ( const ScalarArrayType &Array, aol::CopyFlag CopyType ) : aol::MultiVector<DataType> ( 0, 0 ) {
    if( imagedim != 1 )
      throw aol::Exception( "Converting a ScalarArray to a MultiArray is only possible in case imagedim == 1!", __FILE__, __LINE__);

    this->appendReference ( * ( new ScalarArrayType ( Array, CopyType ) ), true );
  }

  template <qc::Dimension Dim>
  explicit MultiArray ( const qc::GridSize<Dim> &Size, DataType **Data , aol::CopyFlag CopyType = aol::FLAT_COPY ) {

    for ( int i = 0; i < imagedim; ++i ) {
      this->appendReference ( * ( new ScalarArrayType ( Data[i], Size, CopyType ) ), true );
    }
  }

  explicit MultiArray ( const char *FileName ) : aol::MultiVector<DataType> ( 0, 0 ) {
    for ( int i = 0; i < imagedim; ++i ) {
      this->appendReference ( * ( new ScalarArrayType ), true );
    }

    if ( imagedim == 1 )
      comp ( 0 ).load ( FileName );
    else if ( aol::fileNameEndsWith ( FileName, ".png" ) )
      loadPNG ( FileName );
    else
      throw aol::Exception( "Constructor from file name only implemented for \"imagedim == 1\" or PNGs", __FILE__, __LINE__);
  }

  explicit MultiArray ( const string &ComponentOneFileName, const string &ComponentTwoFileName ) : aol::MultiVector<DataType> ( 0, 0 ) {
    if( imagedim != 2 )
      throw aol::Exception( "Creating a MultiArray from two files is only possible in case imagedim == 2!", __FILE__, __LINE__);

    this->appendReference ( * ( new ScalarArrayType ( ComponentOneFileName ) ), true );
    this->appendReference ( * ( new ScalarArrayType ( ComponentTwoFileName ) ), true );

    if ( comp ( 0 ).getSize() != comp ( 1 ).getSize() )
      throw aol::DimensionMismatchException( "Both array files need to have the same size!", __FILE__, __LINE__);
  }

  void addToAll ( const DataVecType& Val ) {
    for ( int i = 0; i < imagedim; ++i )
      comp ( i ).addToAll ( Val[i] );
  }

  using aol::MultiVector<DataType>::addToAll;

  template <typename IndexType>
  void set ( const IndexType& pos, const aol::Vec<imagedim, DataType>& val ) {
    for ( int i = 0; i < imagedim; ++i )
      comp ( i ).set ( pos, val [i] );
  }

  // explicitly overload one specific set from MultiVector
  void set ( int pos, const aol::Vec<imagedim, DataType>& val ) {
    this->template set<int> (pos, val);
  }

  using aol::MultiVector<DataType>::set;

  void setAll ( const DataVecType& val ) {
    for ( int i = 0; i < imagedim; ++i )
      comp ( i ).setAll ( val[i] );
  }

  using aol::MultiVector<DataType>::setAll;

  using aol::MultiVector<DataType>::getTo;

  template <typename IndexType>
  DataVecType get ( const IndexType& pos ) const {
    DataVecType res;
    for ( int i = 0; i < imagedim; ++i )
      res [i] = comp ( i ).get ( pos );
    return res;
  }

  DataVecType get ( const int I ) const {
    // this is a method that would belong to MultiVector, but cannot be implemented there (method cannot have guessed template return value)
    DataVecType res;
    this->getTo ( I, res );
    return ( res );
  }

  const ScalarArrayType& comp ( int i ) const {
    return static_cast<const ScalarArrayType&> ( * ( this->vecs[ i ].ptr ) );
  }

  ScalarArrayType& comp ( int i ) {
    return static_cast<ScalarArrayType&> ( * ( this->vecs[ i ].ptr ) );
  }

  const ScalarArrayType& operator[] ( int i ) const {
    return static_cast< const ScalarArrayType& > ( * ( this->vecs[ i ].ptr ) );
  }

  ScalarArrayType& operator[] ( int i ) {
    return static_cast< ScalarArrayType& > ( * ( this->vecs[ i ].ptr ) );
  }

  MultiArray <DataType, rangedim, imagedim >& operator= ( const ScalarArrayType &Other ) {
    if( imagedim != 1 )
      throw aol::Exception( "Converting a ScalarArray to a MultiArray is only possible in case imagedim == 1!", __FILE__, __LINE__);

    comp ( 0 ) = Other;

    return *this;
  }

  int getNumX ( ) const {
    // MultiArray assumes that all components have the same size.
    return ( comp ( 0 ).getNumX() );
  }

  int getNumY ( ) const {
    // MultiArray assumes that all components have the same size.
    return ( comp ( 0 ).getNumY() );
  }

  aol::Vec3<int> getSize() const {
    // MultiArray assumes that all components have the same size.
    return ( comp ( 0 ).getSize() );
  }

  void resampleFrom ( const MultiArray& array ) {
    for ( int i = 0; i < imagedim; ++i )
      comp ( i ).resampleFrom ( array.comp ( i ) );
  }

  //! Works like flipFrom on qc::ScalarArray.
  void flipFrom ( const qc::MultiArray<DataType, rangedim, imagedim> &InputArray, const qc::Comp Component ) {
    for ( int i = 0; i < imagedim; ++i )
      comp ( i ).flipFrom ( InputArray.comp ( i ), Component );
  }

  //! Works like copyBlockTo on qc::ScalarArray.
  void copyBlockTo ( const aol::Vec<rangedim, int> Start,  qc::MultiArray<DataType, rangedim, imagedim> &block ) const {
    for ( int i = 0; i < imagedim; ++i )
      comp ( i ).copyBlockTo ( Start, block.comp ( i ) );
  }

  //! @todo Special case <int,2,3>, generalize!
  void save ( ostream& os ) const {
    os << "P3" << endl
    << comp ( 0 ).getNumX () << " " << comp ( 0 ).getNumY () << endl
    << 256 << endl;

    for ( int y = 0; y < comp ( 0 ).getNumY (); ++y )
      for ( int x = 0; x < comp ( 0 ).getNumX (); ++x ) {
        for ( int i = 0; i < imagedim; ++i )
          os << static_cast<unsigned int> ( comp ( i ).get ( x, y ) ) << " ";
        os << endl;
      }
  }

  //! @todo Special case of binary writing using FILE* (much faster for gcc3 than using ofstream) for <unsigned char,2,3>, generalize!
  void save ( FILE *file ) const {
    const int numX = comp ( 0 ).getNumX();
    const int numY = comp ( 0 ).getNumY();
    fprintf ( file, "P6 %d %d 255\n", numX, numY );
    for ( int y = 0; y < numY; y++ ) {
      for ( int x = 0; x < numX; x++ ) {
        for ( int i = 0;i < 3;i++ ) {
          fwrite ( & ( comp ( i ) [y*numX+x] ), 1, 1, file );
        }
      }
    }
  }


  //! Call load for each component, generating file names from suitable mask
  void load ( const char *fileNameMask ) {
    string fileNameMaskString = fileNameMask;
    const int numPercentSymbols = count(fileNameMaskString.begin(), fileNameMaskString.end(), '%' );

    if ( numPercentSymbols > 1 )
      throw aol::Exception( "MultiArray::load() : Invalid filename mask!", __FILE__, __LINE__);

    if ( ( numPercentSymbols == 0 ) && ( imagedim != 1 ) )
      throw aol::Exception( string ( "MultiArray::load() : Masks without percentage symbol are only valid for imagedim == 1! (specified mask is " ) + fileNameMask + ")", __FILE__, __LINE__);

    // In this special case, MultiArray::load should behave like ScalarArray::load.
    if ( ( numPercentSymbols == 0 ) && ( imagedim == 1 ) ) {
      comp ( 0 ).load ( fileNameMask );
      return;
    }

    for ( int i = 0; i < imagedim; ++i ) {
      char fileName[1024];
      sprintf ( fileName, fileNameMask, i );
      comp ( i ).load ( fileName );
    }
  }

  void loadPNG ( const char *FileName );

  void savePNG ( const char *FileName ) const;

  void loadPPM ( const char *FileName );

  //! Call save for each component, generating file names from suitable mask
  void save ( const char *fileNameMask, const SaveType type, const char* comment = NULL ) const {
    for ( int i = 0; i < imagedim; ++i ) {
      char fileName[1024];
      sprintf ( fileName, fileNameMask, i );
      comp ( i ).save ( fileName, type, comment );
    }
  }

  //! Call saveRaw for each component, generating file names from suitable mask using one clipping value
  void saveRaw ( const char *fileNameMask, const SaveType type, const DataType Minimum, const DataType Maximum ) {
    for ( int i = 0; i < imagedim; ++i ) {
      char fileName[1024];
      sprintf ( fileName, fileNameMask, i );
      comp ( i ).saveRaw ( fileName , type, Minimum, Maximum );
    }
  }

  //! Call saveRaw for each component, generating file names from suitable mask using multiple clipping values
  void saveRaw ( const char *fileNameMask, const SaveType type, const DataVecType &Minima, const DataVecType &Maxima ) {
    for ( int i = 0; i < imagedim; ++i ) {
      char fileName[1024];
      sprintf ( fileName, fileNameMask, i );
      comp ( i ).saveRaw ( fileName , type, Minima[i], Maxima[i] );
    }
  }

  //! Calls setQuietMode on each component.
  void setQuietMode ( bool QuietMode ) {
    for ( int i = 0; i < imagedim; ++i )
      comp ( i ).setQuietMode(QuietMode);
  }

  //! Calls setOverflowHandling on each component.
  void setOverflowHandling ( aol::OverflowHandlingType Type, DataType Min, DataType Max ){
    for ( int i = 0; i < imagedim; ++i )
      comp ( i ).setOverflowHandling(Type, Min, Max);
  }

  void setOverflowHandlingToCurrentValueRange ( ) {
    const DataType minValue = this->getMinValue();
    const DataType maxValue = this->getMaxValue();
    setOverflowHandling ( aol::CLIP_THEN_SCALE, minValue, ( maxValue > minValue ) ? maxValue : minValue + 1 );
  }

  //! Works like pasteFrom on ScalarArray<QC_2D>.
  void pasteFrom ( const qc::MultiArray<DataType, rangedim, imagedim> &Array, const int XPosition, const int YPosition ) {
    for ( int i = 0; i < imagedim; ++i )
      comp ( i ).pasteFrom(Array[i], XPosition, YPosition);
  }

  //! Works like shiftByOffsetFrom on ScalarArray.
  void shiftByOffsetFrom ( const qc::CoordType &Offset, const MultiArray &Original ) {
    for ( int i = 0; i < imagedim; ++i )
      comp ( i ).shiftByOffsetFrom( Offset, Original[i] );
  }

  //! reallocate all ScalarArrays contained, but do not change number of components
  template< typename Structure >
  void reallocate ( const Structure &struc ) {
    for ( int i = 0; i < imagedim; ++i )
      comp(i).reallocate ( struc );
  }

  //! reallocate all ScalarArrays contained, but do not change number of components
  void reallocate ( const int NumX, const int NumY, const int NumZ=1 ) {
    aol::Vec3<int> gridSize( NumX, NumY, NumZ );
    qc::GridStructure tempGrid ( gridSize, (NumZ == 1) ? qc::QC_2D : qc::QC_3D );
    for ( int i = 0; i < imagedim; ++i )
      comp(i).reallocate ( tempGrid );
  }

  void getPointWiseNorm ( typename ScalarArrayTrait<typename aol::RealTrait<DataType>::RealType, rangedim>::ArrayType &dest ) const {
    if ( rangedim != imagedim )
      throw aol::Exception ( "qc::MultiArray::getPointWiseNorm does not make sense for rangedim != imagedim", __FILE__, __LINE__ );

    for ( int posi = 0; posi < comp(0).size(); ++posi ) {
      DataVecType tmp;
      for ( int i = 0; i < imagedim; ++i )
        tmp [i] = comp ( i )[posi];

      dest[posi] = tmp.norm();
    }
  }

  typename aol::RealTrait<DataType>::RealType getMaxOfPointWiseNorm ( ) const {
    typename aol::RealTrait<DataType>::RealType maxNorm = 0;

    for ( int posi = 0; posi < comp(0).size(); ++posi ) {
      DataVecType tmp;
      for ( int i = 0; i < imagedim; ++i )
        tmp [i] = comp ( i )[posi];

      maxNorm = aol::Max ( maxNorm, tmp.norm() );
    }
    return maxNorm;
  }

  aol::Vec<imagedim, typename aol::RealTrait<DataType>::RealType> getMeanValue ( ) const {
    aol::Vec<imagedim, typename aol::RealTrait<DataType>::RealType> meanValue;
    for ( int i = 0; i < imagedim; ++i )
      meanValue[i] = comp ( i ).getMeanValue();
    return meanValue;
  }

  //! Visualize matrix as pixel graphic: color range centered at zero (green), smaller values down to blue, bigger values up to red. Entries that are exactly zero are shown as white.
  template < typename MatrixType >
  void visualizeMatrix ( const MatrixType& mat ) {
    if ( rangedim != 2 || imagedim != 3 )
      throw aol::Exception ( "qc::MultiArray::visualizeMatrix only works for rangedim = 2, imagedim = 3", __FILE__, __LINE__ );

    ScalarArray< typename MatrixType::DataType, qc::QC_2D > arr;
    arr.visualizeMatrix( mat );
    this->reallocate ( arr );

    const typename MatrixType::DataType mx = arr.getMaxAbsValue();
    aol::RGBColorMap< typename MatrixType::DataType > cMap ( -mx, mx, aol::RGBColorMap< typename MatrixType::DataType >::HSV_BLUE_TO_RED );
    for ( int i = 0; i < arr.size(); ++i ) {
      aol::Vec3< typename MatrixType::DataType > col ( 1, 1, 1 ); // white if entry is zero
      if ( arr[i] != 0 ) {
        cMap.scalarToColor ( arr[i], col );
      }
      aol::Vec3<unsigned char> color ( static_cast<unsigned char> ( 255 * col[0] ), static_cast<unsigned char> ( 255 * col[1] ), static_cast<unsigned char> ( 255 * col[2] ) );

      this->set ( i, color );
    }
  }

  //! lossless save in binary quocmesh format
  void saveToFile ( const char *FileName ) const {
    aol::Bzipofstream writer ( FileName );
    writer << aol::VectorFileMagicChar::MultiArray << rangedim << imagedim << aol::FileFormatMagicNumber<double>::FFType << endl;
    writer << "# This is a QuOcMesh file of type " << aol::VectorFileMagicChar::MultiArray << rangedim << imagedim << aol::FileFormatMagicNumber<double>::FFType
    << " storing a qc::MultiArray<" << aol::FileFormatMagicNumber<DataType>::FFContainedName << ", " << rangedim << ", " << imagedim << ">" << endl;

    for ( int i = 0; i < imagedim; ++i ) {
      for ( int j = 0; j < rangedim; ++j ) {
        writer << ( (*this)[i].getSize() )[j] << " ";
      }
    }
    writer << endl;

    for ( int i = 0; i < imagedim; ++i ) {
      const char* buffer = reinterpret_cast<char*> ( (*this)[i].getData() );
      writer.write ( buffer, (*this)[i].size() * sizeof ( DataType ) );
    }
  }

  //! load from file in quocmesh format, no conversion of data types
  void loadFromFile ( const char *FileName ) {
    aol::Bzipifstream reader ( FileName );

    char M = 0, R = 0, I = 0;
    int ident;
    reader >> M;
    reader >> R;
    reader >> I;
    reader >> ident;
    if ( ( M != aol::VectorFileMagicChar::MultiArray ) || ( R != rangedim + '0' ) || ( I != imagedim + '0' ) ||  ( ident != aol::FileFormatMagicNumber<double>::FFType ) ) {
      cerr << M << R << I << ident << ", should be " << aol::VectorFileMagicChar::MultiArray << rangedim << imagedim << aol::FileFormatMagicNumber<double>::FFType << endl;
      throw aol::Exception ( "Illegal magic number for loadMultiArrayFromFile", __FILE__, __LINE__ );
    }
    aol::READ_COMMENTS ( reader );

    for ( int i = 0; i < imagedim; ++i ) {
      qc::CoordType size;
      for ( int j = 0; j < rangedim; ++j ) {
        reader >> size[j];
      }
      qc::Array<DataType> tmp ( size );
      (*this)[i].reallocate ( tmp );
    }

    {
      char buffer[1024];
      reader.getline ( buffer, 2 ); // space and endl
    }

    for ( int i = 0; i < imagedim; ++i ) {
      reader.read ( reinterpret_cast<char*>( (*this)[i].getData() ), (*this)[i].size() * sizeof( DataType ) );
    }
  }
};

} // end namespace

#endif
