#ifndef __SCALARARRAY_H
#define __SCALARARRAY_H

#include <array.h>
#include <quoc.h>
#include <vec.h>
#include <smallMat.h>
#include <auxiliary.h>
#include <kernel1d.h>
#include <kernel2d.h>
#include <bzipiostream.h>
#include <simplexGrid.h>

#include <kernel3d.h>
#include <rectangularGrid.h>
#include <iterators.h>

#ifdef USE_LIB_HDF5
#include <hdf5.h>
#include <H5Cpp.h>
#include <stack>
#endif

namespace qc {

void netCDFgetDepthFirstDataAndGroups ( const char *fileName, std::vector<std::string> &groupNames, std::string &dataName, unsigned int dimension );

template <typename _DataType, qc::Dimension Dim> class ScalarArray;

template <typename _DataType>
class ScalarArray<_DataType, qc::QC_1D> : public Array<_DataType> {
public:
  typedef _DataType DataType;
  typedef typename aol::Vector<DataType>::RealType RealType;
  static const qc::Dimension Dim = qc::QC_1D;
  // constructor
  explicit ScalarArray ( int NumX ) :
    Array<DataType> ( NumX, 1, 1 ) {}

  ScalarArray() :
    Array<DataType> ( 0, 0 ) {}

  explicit ScalarArray ( const aol::Vector<DataType> &Vector, aol::CopyFlag copyFlag = aol::FLAT_COPY ) :
    Array<DataType> ( Vector, Vector.size(), 1, 1, copyFlag ) {}

  ScalarArray ( const aol::Vector<DataType> &Vector, int NumX, aol::CopyFlag copyFlag = aol::FLAT_COPY ) :
    Array<DataType> ( Vector, NumX, 1, 1, copyFlag ) {}

  template <typename GridType>
  ScalarArray ( const aol::Vector<DataType> &Vector, const GridType &Grid, aol::CopyFlag copyFlag = aol::FLAT_COPY ) :
    Array<DataType> ( Vector, Grid.getNumX(), 1, 1, copyFlag ) {}

  ScalarArray ( const aol::Vector<DataType> &Vector, const GridSize<QC_1D> &GridSize, aol::CopyFlag copyFlag = aol::FLAT_COPY ) :
    Array<DataType> ( Vector, GridSize.getNumX(), 1, 1, copyFlag ) {}

  ScalarArray ( const Array<DataType> &A, aol::CopyFlag copyFlag ) :
    Array<DataType> ( A, A.getNumX(), 1, 1, copyFlag ) {}

  ScalarArray ( int NumX, DataType *Data, aol::CopyFlag copyFlag = aol::FLAT_COPY ) :
    Array<DataType> ( NumX, 1, 1, Data, copyFlag ) {}

  explicit ScalarArray ( const GridStructure &Grid ) :
    Array<DataType> ( Grid.getNumX(), 1, 1 ) {
    if ( Grid.getDimOfWorld() != 1 )
      throw aol::Exception ( "ScalarArray<QC_1D>::initialized with a grid whose dimension does not equal 1!\n", __FILE__, __LINE__ );
  }

  explicit ScalarArray ( const GridSize<QC_1D> &GridSize ) :
    Array<DataType> ( GridSize.getNumX(), 1, 1 ) {}

  /** Copy constructor
   */
  explicit ScalarArray ( const qc::ScalarArray<DataType, qc::QC_1D> &org, aol::CopyFlag copyFlag = aol::DEEP_COPY ) :
    Array<DataType> ( org, copyFlag ) {}

  //! Read data from file
  //! and set size of array correctly automatically
  explicit ScalarArray ( const string &filename );

  void reallocate ( const int NumX ) {
    Array<DataType>::reallocate ( NumX, 1 );
  }

  template< typename Structure >
  void reallocate ( const Structure &struc ) {
    reallocate ( struc.getNumX() );
  }

  // We do not want to be able to use all set and add methods from the parent class,
  // so we explicitly implement those that should be usable.
  void set ( const CoordType &Coord, DataType V ) {
    Array<DataType>::set ( Coord, V );
  }

  using aol::Vector<DataType>::set;

  void add ( const CoordType &Pt, DataType V ) {
    Array<DataType>::add ( Pt, V );
  }

  using aol::Vector<DataType>::add;

  DataType get ( const CoordType &Coord ) const {
    return Array<DataType>::get ( Coord );
  }

  DataType get ( const aol::Vec<1, int> &Coord ) const {
    return aol::Vector<DataType>::get ( Coord[0] );
  }

  using aol::Vector<DataType>::get;

  DataType &getReference ( int X ) {
    return Array<DataType>::getReference ( X, 0, 0 );
  }

  DataType &getReference ( const CoordType &Coords ) {
    return Array<DataType>::getReference ( Coords );
  }

  ScalarArray<DataType, qc::QC_1D>& operator= ( const ScalarArray<DataType, qc::QC_1D> &Other ) {
    if ( this == &Other ) // correct self-assignment
      return ( *this );
    if ( this->numX != Other.numX ||  this->numY != Other.numY || this->numZ != Other.numZ )
      throw aol::Exception ( "qc::ScalarArray<QC_1D>::operator= trying to assign incompatible array", __FILE__, __LINE__ );
    aol::Vector<DataType>::operator= ( Other );
    return *this;
  }

  ScalarArray<DataType, qc::QC_1D>& operator= ( const aol::Vector<DataType> &Other ) {
    aol::Vector<DataType>::operator= ( Other );
    return *this;
  }

  RealType interpolate ( const aol::Vec<1, RealType> &Pt ) const {
    return aol::Vector<DataType>::interpolate ( Pt[0] );
  }

  DataType getPeriodic ( int X ) const {
    if ( X < 0 ) X += this->numX;
    else if ( X > this->numX - 1 ) X = X % this->numX;
    return this->get ( X );
  }

  void copyBlockTo ( const aol::Vec<1, int> Start, qc::ScalarArray<DataType, qc::QC_1D> &Block ) const {
    for ( int x = 0; x < Block.getNumX(); ++x ) {
      Block[x] = this->get ( x + Start[0] );
    }
  }

  void resampleFrom ( const ScalarArray<DataType, qc::QC_1D> &Input ) {
    aol::Vector<DataType>::resampleFrom ( static_cast<aol::Vector<DataType> > ( Input ) );
  }

  /** Load array in 1d pgm style format from the given stream. Note that if you would like
   *  to load compressed files directly you must use the method load(char *fileName),
   *  since only this method can open a pipe and determine the type of compression
   *  @throws aol::TypeException if the magic number of the file is unknown, or the dimensions of the array contained in that file do not match this one
   *  @throws aol::IOException if an io error occured
   *  \author Preusser, Droske, Teusner
   */
  void load ( istream &in );

  /** Load array in 1d pgm style format from the file named fileName. If the file contains
   *  compressed data in the bz2, gz or Z format, it will automatically be decompressed
   *  by a pipe stream that is passed to load(istream &in).
   *  @throws aol::FileException if the file could not be opened
   *  \author Preusser, Droske, Teusner
   */
  void load ( const char *fileName );

  void loadRaw ( istream &in, const int Type, const int InWidth );

  void loadMetaImageFile ( const char *fileName );

  void saveASCII ( const char *FileName, const aol::Format &OutFormat = aol::scientificFormat ) const {
    std::ofstream file ( FileName );
    for ( int i = 0; i < this->getNumX(); ++i )
      file << OutFormat ( this->get ( i ) ) << "\t";
    file << std::endl;
    file.close();
  }


  /** Save aray 1d pgm style format of type type
   *  @param out The stream to write the data to
   *  @param type The file type (see qc::SaveType)
   *  @param comment The comment to be written into the header of the file. If comment is NULL the date and the data type will be written as comment
   *  @throws aol::IOException if an io error occured
   *  \author Preusser, Droske, Teusner
   */
  void save ( ostream &out, SaveType type, const char *comment = NULL ) const;

  /** Save aray 1d pgm style format of type type.
   *  This method is capable of directly compressing the data using bzip2, gzip or compress.
   *  The compressed files can be read again directly through a call of load(char *fileName).
   *  See qc::SaveType for details on file types
   *  @param fileName The name of the file to be written. If compressed data shall be saved the suffix should correspond to the compression type (i.e. '.bz2', '.gz' or '.Z'). Other wise a warning will be displayed
   *  @param comment The comment to be written into the header of the file. If comment is NULL the date and the data type will be written as comment
   *  @throws aol::TypeException if the given type is incompatible
   *  @throws aol::IOException if an io error occured
   *  @throws aol::FileException if the file could not be opened
   *  @throws aol::PipeException if the compression pipe could not be opened
   *  \author Preusser, Droske, Teusner
   */
  void save ( const char *fileName, SaveType type, const char *comment = NULL ) const;

  //! \author Berkels
  void saveMetaImageFile ( const char *BaseFileName, SaveType Type ) const;

  //! lossless save in binary quocmesh format
  void saveToFile ( const char *filename ) const;

  //! load from file in quocmesh format, no conversion of data types
  void loadFromFile ( const char *filename );

  DataType getMedianFilterValue ( int X ) const;

  DataType getConvolveValue ( int X, const Kernel1d<RealType> &Kernel ) const;

  //! Apply the linear filter kernel to this array and store the results in FilteredOutput.
  void applyLinearFilterTo ( const Kernel1d<RealType> &Kernel, ScalarArray<RealType, qc::QC_1D> &FilteredOutput ) const {
    for ( int x = 0; x < this->numX; ++x )
      FilteredOutput.set ( x, getConvolveValue( x, Kernel ) );
  }

  void applyMedianFilter();
};

/** Implementation of a 2D array used to contain images mapped on a vector
 * @todo Make clear the problem of mapping indices via IndexMapper
 */
template <typename _DataType>
class ScalarArray<_DataType, qc::QC_2D> : public Array<_DataType> {
public:
  typedef _DataType DataType;
  typedef typename aol::Vector<DataType>::RealType RealType;
  static const qc::Dimension Dim = qc::QC_2D;

  // --------Functor classes used by fillValuesFrom---------------------

  // \brief Standart copying class
  // \author toelkes
  class StandardCopier {
  public:
    inline bool sourceRangeOk ( const ScalarArray<DataType, qc::QC_2D> &source, const int sRanX[2], const int sRanY[2], const int tRanX[2], const int tRanY[2] ) {
      // Check if only valid values of source are requested in getValuesFrom().
      return ( sRanX[0] + ( tRanX[1] - tRanX[0] ) <= source.getNumX() && sRanY[0] + ( tRanY[1] - tRanY[0] ) <= source.getNumY() );
    }

    inline DataType operator() ( int x, int y, const ScalarArray<DataType, qc::QC_2D> &source ) {
      return source.get ( x, y );
    }
  };

  // \brief Standart copying class for adding values instead of overwriting
  // \author toelkes
  class AddCopier {
  public:
    AddCopier ( const ScalarArray<DataType, qc::QC_2D> &origDestValues )
      : _origDestValues ( origDestValues )
    {}

    inline bool sourceRangeOk ( const ScalarArray<DataType, qc::QC_2D> &source, const int sRanX[2], const int sRanY[2], const int tRanX[2], const int tRanY[2] ) {
      // Check if only valid values of source are requested in getValuesFrom().
      return ( sRanX[0] + ( tRanX[1] - tRanX[0] ) <= source.getNumX() && sRanY[0] + ( tRanY[1] - tRanY[0] ) <= source.getNumY() );
    }

    inline DataType operator() ( int x, int y, const ScalarArray<DataType, qc::QC_2D> &source ) {
      return ( _origDestValues.get ( x, y ) + source.get ( x, y ) );
    }
  private:
    // saved as reference, so be careful when using this class
    const ScalarArray<DataType, qc::QC_2D> _origDestValues;
  };

  // \brief Class used in shiftByOffset. Copies and shifts by a given offset.
  // \author toelkes
  class OffsetCopier {
  public:
    OffsetCopier ( const int xOffset, const int yOffset, const int numX, const int numY )
      : _xOffset ( xOffset ), _yOffset ( yOffset ), _numX ( numX ), _numY ( numY ) {
      if ( abs ( _xOffset ) > _numX || abs ( _yOffset ) > _numY ) {
        throw aol::Exception ( "OffsetCopier: Offset larger than array!\n", __FILE__, __LINE__ );
      }
    }

    inline bool sourceRangeOk ( const ScalarArray<DataType, qc::QC_2D> &source, const int sRanX[2], const int sRanY[2], const int tRanX[2], const int tRanY[2] ) {
      // Check if the offsets are valid. Is checked before, but if the exception is catched this could still be called.
      if ( ( abs ( _xOffset ) > _numX ) || ( abs ( _yOffset ) > _numY ) ) {
        throw aol::Exception ( "OffsetCopier: Cannot use OffsetCopier with offsets bigger than the grid size!\n", __FILE__, __LINE__ );
        return false;
      }

      // Check if only valid values of source are requested in getValuesFrom().
      return ( ( sRanX[0] + ( tRanX[1] - tRanX[0] ) <= source.getNumX() ) && ( sRanY[0] + ( tRanY[1] - tRanY[0] ) <= source.getNumY() ) );
    }

    inline DataType operator() ( int x, int y, const ScalarArray<DataType, qc::QC_2D> &source ) {
      return source.get ( ( ( x - _xOffset ) >= 0 ) ? ( x - _xOffset ) % _numX : ( _numX + ( x - _xOffset ) ),
                          ( ( y - _yOffset ) >= 0 ) ? ( y - _yOffset ) % _numY : ( _numY + ( y - _yOffset ) ) );
    }

  private:
    const int _xOffset, _yOffset, _numX, _numY;
  };

  // \brief Class used in resampleFrom. Copies and scales/interpolates.
  // \author toelkes
  class ResampleCopier {
  public:
    ResampleCopier ( const ScalarArray<DataType, qc::QC_2D> &source, const ScalarArray<DataType, qc::QC_2D> &target ) {
      _hx = static_cast<RealType> ( source.getNumX() - 1 ) / ( target.getNumX() - 1 );
      _hy = static_cast<RealType> ( source.getNumY() - 1 ) / ( target.getNumY() - 1 );
    }

    inline bool sourceRangeOk ( const ScalarArray<DataType, qc::QC_2D> &source, const int sRanX[2], const int sRanY[2], const int tRanX[2], const int tRanY[2] ) {
      // Only values inside of source can be interpolated.
      return ( ( static_cast<int> ( ( sRanX[0] + ( tRanX[1] - tRanX[0] ) - 1 ) * _hx ) <= ( source.getNumX() - 1 ) )
            && ( static_cast<int> ( ( sRanY[0] + ( tRanY[1] - tRanY[0] ) - 1 ) * _hy ) <= ( source.getNumY() - 1 ) ) );
    }

    inline DataType operator() ( int x, int y, const ScalarArray<DataType, qc::QC_2D> &source ) {
      const RealType scaledX = x * _hx, scaledY = y * _hy;
      return static_cast<DataType> ( source.interpolate ( scaledX, scaledY ) );
    }

  private:
    // long double precision is crucial here.
    long double _hx, _hy;
  };

  // \brief Class used in flipFrom. Copies and flips.
  // \author toelkes, berkels
  class flipCopier {
  public:
    flipCopier ( const int numX, const int numY, const qc::Comp component )
      : _numX ( numX ), _numY ( numY ), _component ( component )
    {}

    inline bool sourceRangeOk ( const ScalarArray<DataType, qc::QC_2D> &source, const int sRanX[2], const int sRanY[2], const int tRanX[2], const int tRanY[2] ) {
      // Check if only valid values of source are requested in getValuesFrom().
      return ( sRanX[0] + ( tRanX[1] - tRanX[0] ) <= source.getNumX()
               && sRanY[0] + ( tRanY[1] - tRanY[0] ) <= source.getNumY() );
    }

    inline DataType operator() ( int x, int y, const ScalarArray<DataType, qc::QC_2D> &source ) {
      if ( _component == qc::QC_X )
        return source.get ( _numX - 1 - x, y );
      else if ( _component == qc::QC_Y )
        return source.get ( x, _numY - 1 - y );
      else {
        throw aol::Exception ( "flipCopier: Component to flip is invalid!\n", __FILE__, __LINE__ );
        return source.get ( x, y );
      }
    }

  private:
    const int _numX, _numY;
    const qc::Comp _component;
  };

  //--------------------------------------------------------------

  // constructor
  ScalarArray ( int NumX, int NumY ) :
    Array<DataType> ( NumX, NumY ) {
    init();
  }

  ScalarArray ()
    : Array<DataType> ( 0, 0 ) {
    init();
  }

  explicit ScalarArray ( int Width )
    : Array<DataType> ( Width, Width ) {
    init();
  }

  //! this constructor can be used to treat a vector like a ScalarArray
  //! (one modifies the value of the vector by modifying the value of the ScalarArray and does not need to care for the indices)
  ScalarArray ( const aol::Vector<DataType> &Vector, int NumX, int NumY, aol::CopyFlag copyFlag = aol::FLAT_COPY ) :
    Array<DataType> ( Vector, NumX, NumY, 1, copyFlag ) { // default: flat copy
    init();
  }

  ScalarArray ( int NumX, int NumY, DataType *Data, aol::CopyFlag copyFlag = aol::FLAT_COPY ) :
    Array<DataType> ( NumX, NumY, 1, Data, copyFlag ) {
    init();
  }

  ScalarArray ( DataType *Data, const GridSize<QC_2D> &Size, aol::CopyFlag copyFlag = aol::FLAT_COPY ) :
    Array<DataType> ( Size.getNumX(), Size.getNumY(), 1, Data, copyFlag ) {
    init();
  }

  template <typename GridType>
  ScalarArray ( const aol::Vector<DataType> &Vector, const GridType &Grid, aol::CopyFlag copyFlag = aol::FLAT_COPY ) :
    Array<DataType> ( Vector, Grid.getNumX(), Grid.getNumY(), 1, copyFlag ) { // default: flat copy
    init();
  }

  ScalarArray ( const aol::Vector<DataType> &Vector, const GridSize<QC_2D> &GridSize, aol::CopyFlag copyFlag = aol::FLAT_COPY ) :
    Array<DataType> ( Vector, GridSize.getNumX(), GridSize.getNumY(), 1, copyFlag ) { // default: flat copy
    init();
  }

  ScalarArray ( const Array<DataType> &A, aol::CopyFlag copyFlag ) :
    Array<DataType> ( A, A.getNumX(), A.getNumY(), 1, copyFlag ) { // default: flat copy
    init();
  }

  explicit ScalarArray ( const GridStructure &Grid ) :
    Array<DataType> ( Grid.getNumX(), Grid.getNumY() ) {
    if ( Grid.getDimOfWorld() != 2 ) {
      throw aol::Exception ( "ScalarArray<QC_2D>::initialized with a grid whose dimension does not equal 2!\n", __FILE__, __LINE__ );
    }
    init();
  }

  explicit ScalarArray ( const GridSize<QC_2D> &GridSize ) :
    Array<DataType> ( GridSize.getNumX(), GridSize.getNumY() ) {
    init();
  }

  template <typename CubicGrid, Dimension Dim>
  explicit ScalarArray ( const simplex::GridStructure<CubicGrid, Dim> & SimplexGrid ) :
    Array<DataType> ( SimplexGrid.getNumX(), SimplexGrid.getNumY() ) {
    init();
  }

  /** Copy constructor
   */
  explicit ScalarArray ( const qc::ScalarArray<DataType, qc::QC_2D> &org, aol::CopyFlag copyFlag = aol::DEEP_COPY ) :
    Array<DataType> ( org, copyFlag ),
    diffKernelX ( org.diffKernelX  ? new GaussDiffKernel2d<RealType> ( *org.diffKernelX ) : NULL ),
    diffKernelY ( org.diffKernelY  ? new GaussDiffKernel2d<RealType> ( *org.diffKernelY ) : NULL ),
    diffKernelXX ( org.diffKernelXX ? new GaussDiffKernel2d<RealType> ( *org.diffKernelXX ) : NULL ),
    diffKernelXY ( org.diffKernelXY ? new GaussDiffKernel2d<RealType> ( *org.diffKernelXY ) : NULL ),
    diffKernelYY ( org.diffKernelYY ? new GaussDiffKernel2d<RealType> ( *org.diffKernelYY ) : NULL ) {
    // do not call init()
  }


  //! Read data from file
  //! and set size of array correctly automatically
  explicit ScalarArray ( const string &filename );

  ~ScalarArray() {
    if ( diffKernelX ) delete diffKernelX;
    if ( diffKernelXX ) delete diffKernelXX;
    if ( diffKernelY ) delete diffKernelY;
    if ( diffKernelYY ) delete diffKernelYY;
    if ( diffKernelXY ) delete diffKernelXY;
  }

  void reallocate ( const int NumX, const int NumY ) {
    if ( NumY == -1 )
      Array<DataType>::reallocate ( NumX, NumX );
    else
      Array<DataType>::reallocate ( NumX, NumY );
  }
  
//  void reallocate ( const aol::Vec<2, short> &Size ) {
//    reallocate ( Size[0], Size[1] );
//  }

  template< typename Structure >
  void reallocate ( const Structure &struc ) {
    reallocate ( struc.getNumX(), struc.getNumY() );
  }



  //! Initialize some variables, do not touch data!
  void init() {
    // setOverflowHandling( CLIP, 0, static_cast<DataType>( 255 ) );

    diffKernelX = new GaussDiffKernel2d<RealType> ( this->DIFF_STENCIL, this->DIFF_STD_SIGMA, DIFF_X );
    diffKernelY = new GaussDiffKernel2d<RealType> ( this->DIFF_STENCIL, this->DIFF_STD_SIGMA, DIFF_Y );
    diffKernelXX = new GaussDiffKernel2d<RealType> ( this->DIFF_STENCIL, this->DIFF_STD_SIGMA, DIFF_XX );
    diffKernelXY = new GaussDiffKernel2d<RealType> ( this->DIFF_STENCIL, this->DIFF_STD_SIGMA, DIFF_XY );
    diffKernelYY = new GaussDiffKernel2d<RealType> ( this->DIFF_STENCIL, this->DIFF_STD_SIGMA, DIFF_YY );
  }



  /** Returns a value which is interpolated at the center of the given element
   */
  RealType interpolate ( const Element &El ) const {
    return static_cast<RealType> ( .25 * ( this->get ( El.x(), El.y() ) + this->get ( El.x() + 1, El.y() ) + this->get ( El.x(), El.y() + 1 ) + this->get ( El.x() + 1, El.y() + 1 ) ) );
  }

  /** Returns an interpolated value of the image at position X, Y (in pixels)
   * @param X x-coordinate of the interpolated point
   * @param Y y-coordinate of the interpolated point
   */
  RealType interpolate ( RealType X, RealType Y ) const {
    int xL = static_cast<int> ( X ), yL = static_cast<int> ( Y );

    if ( X >= static_cast<RealType> ( this->numX - 1 ) ) xL = this->numX - 2;
    if ( Y >= static_cast<RealType> ( this->numY - 1 ) ) yL = this->numY - 2;
    if ( X < 0. ) xL = 0;
    if ( Y < 0. ) yL = 0;

    X -= xL;
    Y -= yL;

    const int xU = xL + 1, yU = yL + 1;

    const DataType v[4] = {
      this->get ( xL, yL ), this->get ( xU, yL ),
      this->get ( xL, yU ), this->get ( xU, yU )
    };

    return static_cast<RealType> ( ( 1 - X ) * ( 1 - Y ) * v[0] + X * ( 1 - Y ) * v[1] +
                                   ( 1 - X ) *   Y * v[2] + X *   Y * v[3] );
  }

  /** Returns an interpolated value at position x, y (world coordinates) where the array is assumed to discretize [0,1]^2
   * Throws if the scalar array is not quadratic
   * @param x x-coordinate of interpolated point (in [0,1])
   * @param y y-coordinate of interpolated point (in [0,1])
   */
  RealType interpolate_on01 ( RealType x, RealType y ) const {
#ifdef BOUNDS_CHECK
    // assert that dimensions are equal
    GridSize<QC_2D> sizeChecker ( *this );
    sizeChecker.quadraticOrDie ();
#endif
    return ( interpolate ( x * ( this->getNumX() - 1 ), y * ( this->getNumY() - 1 ) ) );
  }

  RealType interpolate_on01_periodic ( RealType x, RealType y ) const {
    return interpolate_on01 ( fabs ( x - floor ( x ) ), fabs ( y - floor ( y ) ) );
  }

  /** Returns an interpolated value of the image at position X, Y
   * @param pt coordinates of the interpolated point
   */
  RealType interpolate ( const aol::Vec<2, RealType> &pt ) const {
    return interpolate ( pt[0], pt[1] );
  }

  RealType interpolate_on01 ( const aol::Vec<2, RealType> & coord ) const {
    return interpolate_on01 ( coord[0], coord[1] );
  }

  RealType interpolate_on01_periodic ( const aol::Vec<2, RealType> & coord ) const {
    return interpolate_on01_periodic ( coord[0], coord[1] );
  }

  template <class T>
  void hessian ( const CoordType &/*Coord*/, aol::Matrix22<T> &/*Hessian*/ ) {
    cerr << "Not implemented" << endl;
  }

  template <class T>
  void cellGradient ( const Element &el, aol::Vec<2, T> &Grad ) const {
    Grad[ 0 ] = dx ( el );
    Grad[ 1 ] = dy ( el );
  }

  template <class T>
  void cellGradient ( const Element &el, RealType LocalX, RealType LocalY, aol::Vec<2, T> &Grad ) const {
    Grad[ 0 ] = dx ( el, LocalX, LocalY );
    Grad[ 1 ] = dy ( el, LocalX, LocalY );
  }

  template <typename T>
  void gradient ( const aol::Vec<2, T> &pos, aol::Vec<2, T> &grad ) const {
    grad[0] = dx ( pos[0], pos[1] );
    grad[1] = dy ( pos[0], pos[1] );
  }

  //! computes the gradient of the image (with pixel size 1) using central finite differences (one-sided at boundary)
  template <typename T>
  void gradientFD ( const aol::Vec<2, int> &pos, aol::Vec<2, T> &grad ) const {
    grad[0] = dxFD ( pos[0], pos[1] );
    grad[1] = dyFD ( pos[0], pos[1] );
  }

  //! computes the x derivative of the image (with pixel size 1) using central finite differences (one-sided at boundary)
  RealType dxFD ( const int X, const int Y ) const {
    if ( X == 0 ) {
      return static_cast<RealType> ( this->get ( X + 1, Y ) - this->get ( X, Y ) );
    } else if ( X == this->numX - 1 ) {
      return static_cast<RealType> ( this->get ( X, Y ) - this->get ( X - 1, Y ) );
    } else {
      return static_cast<RealType> ( 0.5 * ( this->get ( X + 1, Y ) - this->get ( X - 1, Y ) ) );
    }
  }

  RealType dxFD ( const CoordType &Coord ) const {
    return dxFD ( Coord[0], Coord[1] );
  }

  //! computes the y derivative of the image (with pixel size 1) using central finite differences (one-sided at boundary)
  RealType dyFD ( const int X, const int Y ) const {
    if ( Y == 0 ) {
      return static_cast<RealType> ( this->get ( X , Y + 1 ) - this->get ( X, Y ) );
    } else if ( Y == this->numY - 1 ) {
      return static_cast<RealType> ( this->get ( X, Y ) - this->get ( X, Y - 1 ) );
    } else {
      return static_cast<RealType> ( 0.5 * ( this->get ( X, Y + 1 ) - this->get ( X, Y - 1 ) ) );
    }
  }

  RealType dyFD ( const CoordType &Coord ) const {
    return dyFD ( Coord[0], Coord[1] );
  }

  /**
   * @param   coordinates
   * @returns x-derivative at arbitrary coordinates (with pixel size 1)
   */
  RealType dx ( const RealType &X, const RealType &Y ) const {
    int x = static_cast<int> ( X ), y = static_cast<int> ( Y );

    RealType a2 = Y - y;

    if ( x >= this->numX - 1 || y >= this->numY - 1 ) {
      cerr << "ERROR in ScalarArray<QC_2D>::dx() x = " << x << " y = " << y << endl;
    }

    return static_cast<RealType> ( ( 1. - a2 ) * ( this->get ( x + 1, y ) - this->get ( x, y ) ) + a2 * ( this->get ( x + 1, y + 1 ) - this->get ( x, y + 1 ) ) );
  }

  /**
   * @param   coordinates
   * @returns y-derivative at arbitrary coordinates (with pixel size 1)
   */
  RealType dy ( const RealType &X, const RealType &Y ) const {
    int x = static_cast<int> ( X ), y = static_cast<int> ( Y );

    if ( x >= this->numX - 1 && y >= this->numY - 1 ) {
      cerr << "ERROR in ScalarArray<QC_2D>::dy() x = " << x << " y = " << y << endl;
    }

    RealType a1 = X - x;

    return static_cast<RealType> ( ( 1. - a1 ) * ( this->get ( x, y + 1 ) - this->get ( x, y ) ) + a1 * ( this->get ( x + 1, y + 1 ) - this->get ( x + 1, y ) ) );

  }

  /** Returns the element wise central x - derivative of the image in
   * element e (with pixel size 1)
   * @param e The element on which computation is to be performed
   */
  RealType dx ( const Element &e ) const {
    int x = e.x(), y = e.y();

    if ( x == this->numX - 1 ) --x;
    if ( y == this->numY - 1 ) --y;

    return static_cast<RealType> ( ( this->get ( x + 1, y ) - this->get ( x, y ) + this->get ( x + 1, y + 1 ) - this->get ( x, y + 1 ) ) / 2. );
  }

  /** Returns the element wise central y - derivative of the image in
   * element e  (with pixel size 1)
   * @param e The element on which computation is to be performed
   */
  RealType dy ( const Element &e ) const {
    int x = e.x(), y = e.y();

    if ( x == this->numX - 1 ) --x;
    if ( y == this->numY - 1 ) --y;

    return static_cast<RealType> ( ( this->get ( x, y + 1 ) - this->get ( x, y ) + this->get ( x + 1, y + 1 ) - this->get ( x + 1, y ) ) / 2. );
  }

  inline RealType dx ( const Element &El, RealType /*LocalX*/, RealType LocalY ) const {
    const int x = El.x();
    const int y = El.y();
    return static_cast<RealType> ( ( this->get ( x + 1, y ) - this->get ( x, y ) ) * ( 1. - LocalY ) + ( this->get ( x + 1, y + 1 ) - this->get ( x, y + 1 ) ) * LocalY );
  }

  inline RealType dy ( const Element &El, RealType LocalX, RealType /*LocalY*/ ) const {
    const int x = El.x();
    const int y = El.y();
    return static_cast<RealType> ( ( this->get ( x, y + 1 ) - this->get ( x, y ) ) * ( 1. - LocalX ) + ( this->get ( x + 1, y + 1 ) - this->get ( x + 1, y ) ) * LocalX );
  }

  inline RealType dx ( const aol::Vec<2, RealType> &Pt ) const {
    return dx ( Pt[0], Pt[1] );
  }

  inline RealType dy ( const aol::Vec<2, RealType> &Pt ) const {
    return dy ( Pt[0], Pt[1] );
  }

  //! errr. torus!
  DataType getPeriodic ( int X, int Y ) const {
    if ( X < 0 ) X += this->numX;
    else if ( X > this->numX - 1 ) X = X % this->numX;
    if ( Y < 0 ) Y += this->numY;
    else if ( Y > this->numY - 1 ) Y = Y % this->numY;
    return this->get ( X, Y );
  }
  
  //! reflection clipping
  DataType getReflection ( int X, int Y ) const {
    if ( X < 0 ) X = -X;
    else if ( X > this->numX - 1 ) X = 2 * ( this->numX - 1 ) - X;
    if ( Y < 0 ) Y = -Y;
    else if ( Y > this->numY - 1 ) Y = 2 * ( this->numY - 1 ) - Y;
    return this->get ( X, Y );
  }

  DataType getClip ( int X, int Y ) const {
    if ( X < 0 ) X = 0;
    else if ( X > this->numX - 1 ) X = this->numX - 1;
    if ( Y < 0 ) Y = 0;
    else if ( Y > this->numY - 1 ) Y = this->numY - 1;
    return this->get ( X, Y );
  }


  void setMax ( int X, int Y, DataType value ) {
    int in = Array<DataType>::index ( X, Y );
    if ( this->_pData[ in ] < value ) this->_pData[ in ] = value;
  }

  void setMin ( int X, int Y, DataType value ) {
    int in = Array<DataType>::index ( X, Y );
    if ( this->_pData[ in ] > value ) this->_pData[ in ] = value;
  }

  // We do not want to be able to use all set and add methods from the parent class, so we explicitly implement those that should be usable.
  // Unfortunately, this cannot be done via "using".
  void set ( const CoordType &Coord, DataType V ) {
    Array<DataType>::set ( Coord, V );
  }

  void set ( const aol::Vec<2, short> &Coord, DataType V ) {
    Array<DataType>::set ( Coord, V );
  }
  
  void set ( const aol::Vec<2, int> &Coord, DataType V ) {
    Array<DataType>::set ( Coord, V );
  }

  // sets a value of the array
  void set ( int X, int Y, DataType V ) {
    Array<DataType>::set ( X, Y, V );
  }

  using aol::Vector<DataType>::set;

  DataType get ( const CoordType &Coord ) const {
    return Array<DataType>::get ( Coord );
  }

  DataType get ( const aol::Vec<2, short> &Coord ) const {
    return Array<DataType>::get ( Coord );
  }
  
  DataType get ( const aol::Vec<2, int> &Coord ) const {
    return Array<DataType>::get ( Coord );
  }

  DataType get ( int X, int Y ) const {
    return Array<DataType>::get ( X, Y );
  }

  using aol::Vector<DataType>::get;

  DataType &getReference ( int X, int Y ) {
    return Array<DataType>::getReference ( X, Y );
  }

  DataType &getReference ( const CoordType &Coords ) {
    return Array<DataType>::getReference ( Coords );
  }

  DataType &getReference ( const aol::Vec<2, short> &Coords ) {
    return Array<DataType>::getReference ( Coords );
  }

  void add ( int X, int Y, DataType V ) {
    Array<DataType>::add ( X, Y, V );
  }

  void add ( const CoordType &Pt, DataType V ) {
    Array<DataType>::add ( Pt, V );
  }

  void add ( const aol::Vec<2, short> &Pt, DataType V ) {
    Array<DataType>::add ( Pt, V );
  }

  DataType* getRowDataPointer ( const int Row ) const {
    return ( this->getData() + this->index ( 0, Row ) );
  }

  using aol::Vector<DataType>::add;

  /**
   * Pastes Image into this array, using (XPosition, YPosition) as top left
   * coordinate of the pasted image in this array.
   */
  void pasteFrom ( const ScalarArray<DataType, qc::QC_2D> &Image, const int XPosition, const int YPosition );

  /**
   * Like pasteFrom but adds the new values to the existing one instead of overwriting them.
   */
  void pasteAddFrom ( const ScalarArray<DataType, qc::QC_2D> &Image, const int XPosition, const int YPosition );

  /**
   * Put the line into the specified row (Comp == QC_X) or column (Comp == QC_Y).
   *
   * \author Berkels
   */
  void putLine ( const qc::Comp Comp, const int Index, const qc::ScalarArray<DataType, qc::QC_1D> &Line ) {
    switch ( Comp ) {
      case qc::QC_X:
        if ( this->getNumY() != Line.getNumX() ) {
          cerr << "getNumY() = " << this->getNumY()
          << "Line.getNumX() = " << Line.getNumX() << endl;
          throw aol::Exception ( "Dimensions don't match", __FILE__, __LINE__ );
        }
        for ( int y = 0; y < this->getNumY(); ++y )
          this->set ( Index, y, Line.get ( y ) );
        break;
      case qc::QC_Y:
        if ( this->getNumX() != Line.getNumX() ) {
          cerr << "getNumX() = " << this->getNumX()
          << "Line.getNumX() = " << Line.getNumX() << endl;
          throw aol::Exception ( "Dimensions don't match", __FILE__, __LINE__ );
        }
        for ( int x = 0; x < this->getNumX(); ++x )
          this->set ( x, Index, Line.get ( x ) );
        break;
      default:
        throw aol::Exception ( "Unexpected qc::Comp", __FILE__, __LINE__ );
    }
  }

  /**
   * putLine counterpart.
   *
   * \author Berkels
   */
  void getLine ( const qc::Comp Comp, const int Index, qc::ScalarArray<DataType, qc::QC_1D> &Line ) const {
    switch ( Comp ) {
      case qc::QC_X:
        if ( this->getNumY() != Line.getNumX() ) {
          cerr << "getNumY() = " << this->getNumY()
          << "Line.getNumX() = " << Line.getNumX() << endl;
          throw aol::Exception ( "Dimensions don't match", __FILE__, __LINE__ );
        }
        for ( int y = 0; y < this->getNumY(); ++y )
          Line.set ( y, this->get ( Index, y ) );
        break;
      case qc::QC_Y:
        if ( this->getNumX() != Line.getNumX() ) {
          cerr << "getNumX() = " << this->getNumX()
          << "Line.getNumX() = " << Line.getNumX() << endl;
          throw aol::Exception ( "Dimensions don't match", __FILE__, __LINE__ );
        }
        for ( int x = 0; x < this->getNumX(); ++x )
          Line.set ( x, this->get ( x, Index ) );
        break;
      default:
        throw aol::Exception ( "Unexpected qc::Comp", __FILE__, __LINE__ );
    }
  }

  /** Load array in 2d pgm style format from the given stream. Note that if you would like
   *  to load compressed files directly you must use the method load(char *fileName),
   *  since only this method can open a pipe and determine the type of compression
   *  @throws aol::TypeException if the magic number of the file is unknown, or the dimensions of the array contained in that file do not match this one
   *  @throws aol::IOException if an io error occured
   *  \author Preusser, Droske
   */
  void load ( istream &in );

  /** Load array in 2d pgm style format from the file named fileName. If the file contains
   *  compressed data in the bz2, gz or Z format, it will automatically be decompressed
   *  by a pipe stream that is passed to load(istream &in).
   *  @throws aol::FileException if the file could not be opened
   *  \author Preusser, Droske
   */
  void load ( const char *fileName );

  void loadRaw ( istream &in, const int Type, const int InWidth, const int InHeight );

  void loadPNG ( const char *fileName );

  void loadTIFF ( const char *fileName );

  /**
   * Loads a 2D array in the meta image data format. Currently only supports a small
   * subset of this format.
   *
   * \author Berkels
   */
  void loadMetaImageFile ( const char *fileName );

  void saveTIFF ( const char *fileName ) const;

  void savePNG ( const char *fileName ) const;

  //! Load array in the Digital Micrograph 3 file format.
  //! \author Berkels
  void loadDM3 ( const char *FileName );

  //! Load array in the NetCDF file format
  //! Uses depth first search to find a data set within the container
  //! \author Mevenkamp, Berkels
  void loadNetCDF ( const char *fileName ) {
    std::vector<std::string> groupNames;
    std::string dataName;
    netCDFgetDepthFirstDataAndGroups ( fileName, groupNames, dataName, 2 );
    loadNetCDF ( fileName, groupNames, dataName );
  }

  //! Load array in the NetCDF file format
  //! Name of groups and contained dataset is manually specified
  //! \author Mevenkamp, Berkels
  void loadNetCDF ( const char *fileName, const std::vector<std::string> &groupNames, const std::string &dataName );

  //! Save array in NetCDF file format (data will be stored in a variable "data" inside the root group)
  //! \author Mevenkamp, Berkels
  void saveNetCDF ( const char *fileName, const char *comment = NULL ) const;

  void saveVTK ( const char *fileName ) const;

  /** Save aray 2d pgm style format of type type
   *  @param out The stream to write the data to
   *  @param type The file type (see qc::SaveType)
   *  @param comment The comment to be written into the header of the file. If comment is NULL the date and the data type will be written as comment
   *  @throws aol::IOException if an io error occured
   *  \author Preusser, Droske
   */
  void save ( ostream &out, SaveType type, const char *comment = NULL ) const;


  void saveASCII ( const char *FileName, const aol::Format &OutFormat = aol::scientificFormat ) const {
    std::ofstream file ( FileName );
    for ( int j = 0; j < this->getNumY(); ++j ) {
      for ( int i = 0; i < this->getNumX(); ++i )
        file << OutFormat ( this->get ( i, j ) ) << "\t";
      file << std::endl;
    }
    file.close();
  }

  /** Save aray 2d pgm style format of type type.
   *  This method is capable of directly compressing the data using bzip2, gzip or compress.
   *  The compressed files can be read again directly through a call of load(char *fileName).
   *  See qc::SaveType for details on file types
   *  @param fileName The name of the file to be written. If compressed data shall be saved the suffix should correspond to the compression type (i.e. '.bz2', '.gz' or '.Z'). Other wise a warning will be displayed
   *  @param comment The comment to be written into the header of the file. If comment is NULL the date and the data type will be written as comment
   *  @throws aol::TypeException if the given type is incompatible
   *  @throws aol::IOException if an io error occured
   *  @throws aol::FileException if the file could not be opened
   *  @throws aol::PipeException if the compression pipe could not be opened
   *  \author Preusser, Droske
   */
  void save ( const char *fileName, SaveType type, const char *comment = NULL ) const;

  //! \author Berkels
  void saveMetaImageFile ( const char *BaseFileName, SaveType Type ) const;

  //! lossless save in binary quocmesh format
  void saveToFile ( const char *filename ) const;

  //! load from file in quocmesh format, no conversion of data types
  void loadFromFile ( const char *filename );

  // needed for element saturation
  DataType getElementSaturationMax ( int X, int Y );

  DataType getElementSaturationMin ( int X, int Y );

  DataType getMedianFilterValue ( int X, int Y ) const;

  DataType getConvolveValue ( int X, int Y, const Kernel2d<RealType> &Kernel ) const;

  //! Apply the linear filter kernel to this array and store the results in FilteredOutput.
  void applyLinearFilterTo ( const Kernel2d<RealType> &Kernel, ScalarArray<RealType, qc::QC_2D> &FilteredOutput ) const {
    for ( int y = 0; y < this->numY; ++y )
      for ( int x = 0; x < this->numX; ++x )
        FilteredOutput.set ( x, y, getConvolveValue( x, y, Kernel ) );
  }

  /**
   * \brief Computes \f$\frac{\int_\Omega f(x)w(x)K(x-\bar x)dx}{\int_\Omega w(x)K(x-\bar x)dx}\f$
   * for kernel \f$K\f$, weight \f$w\f$, at point \f$\bar x\f$.
   */
  DataType getWeightedConvolveValue ( int X, int Y, const Kernel2d<RealType> &Kernel, const ScalarArray<RealType, qc::QC_2D> &Weight ) const;

  DataType getWeightedConvolveValue ( int X, int Y, int Z, const Kernel2d<RealType> &Kernel, const ScalarArray<RealType, qc::QC_2D> &Weight ) const;

  DataType getCannyEdgeValue ( int X, int Y ) const;

  RealType getFilterGradientMag ( int X, int Y ) const {
    return sqrt ( static_cast<RealType> ( aol::Sqr ( getConvolveValue ( X, Y, *diffKernelX ) )
                                          + aol::Sqr ( getConvolveValue ( X, Y, *diffKernelY ) ) ) );
  }

  void     saltAndPepperNoise ( double Frac );

  void     applyMedianFilter();

  void setDiffSigma ( RealType DiffSigma );

  //! currently only copies data
  ScalarArray<DataType, qc::QC_2D>& operator= ( const ScalarArray<DataType, qc::QC_2D> &Other ) {

    if ( this == &Other ) { // correct self-assignment
      return ( *this );
    }

    if ( this->numX != Other.numX ||  this->numY != Other.numY || this->numZ != Other.numZ ) {
      throw aol::Exception ( "qc::ScalarArray<QC_2D>::operator= trying to assign incompatible array", __FILE__, __LINE__ );
    }

    aol::Vector<DataType>::operator= ( Other );

    if ( diffKernelX ) delete diffKernelX;
    if ( diffKernelXX ) delete diffKernelXX;
    if ( diffKernelY ) delete diffKernelY;
    if ( diffKernelYY ) delete diffKernelYY;
    if ( diffKernelXY ) delete diffKernelXY;

    diffKernelX = ( Other.diffKernelX  ? new GaussDiffKernel2d<RealType> ( *Other.diffKernelX ) : NULL );
    diffKernelY = ( Other.diffKernelY  ? new GaussDiffKernel2d<RealType> ( *Other.diffKernelY ) : NULL );
    diffKernelXX = ( Other.diffKernelXX  ? new GaussDiffKernel2d<RealType> ( *Other.diffKernelXX ) : NULL );
    diffKernelXY = ( Other.diffKernelXY  ? new GaussDiffKernel2d<RealType> ( *Other.diffKernelXY ) : NULL );
    diffKernelYY = ( Other.diffKernelYY  ? new GaussDiffKernel2d<RealType> ( *Other.diffKernelYY ) : NULL );

    return *this;
  }

  ScalarArray<DataType, qc::QC_2D>& operator= ( const aol::Vector<DataType> &Other ) {
    aol::Vector<DataType>::operator= ( Other );
    return *this;
  }

  ostream &print ( ostream &os ) const;

  // \brief Makes this array the scaled version of image
  // \author toelkes
  void resampleFrom ( const ScalarArray<DataType, qc::QC_2D> &image ) {
    const int sRanX[2] = { 0, image.getNumX() };
    const int sRanY[2] = { 0, image.getNumY() };
    const int tRanX[2] = { 0, this->getNumX() };
    const int tRanY[2] = { 0, this->getNumY() };

    ResampleCopier copier ( image, *this );
    fillValuesFrom ( image, sRanX, sRanY, tRanX, tRanY, copier, true );
  }

  //! \brief Resamples the input image into this using bicubic resampling.
  //! \author Berkels
  void bicubicResampleFrom ( const ScalarArray<DataType, qc::QC_2D> &Image );

  //! \brief Universal function for copying, padding, pasting...
  //! \author toelkes
  template <typename CopyFunctor>
  void fillValuesFrom ( const ScalarArray<DataType, qc::QC_2D> &source, const int sRanX[2], const int sRanY[2],
                        const int tRanX[2], const int tRanY[2], CopyFunctor getValue, bool differentRangeSizesOk = false ) {
    // Check if the ranges are valid
    if ( sRanX[0] > sRanX[1] || sRanY[0] > sRanY[1]
         || tRanX[0] > tRanX[1] || tRanY[0] > tRanY[1] ) {
      throw aol::Exception ( "ScalarArray<QC_2D>::fillValuesFrom: One or more of the ranges are incorrect: Ran[1] must not be bigger than Ran[0]!\n", __FILE__, __LINE__ );
    }

    // Check if ranges are >= 0
    if ( sRanX[0] < 0 || sRanY[0] < 0 || tRanX[0] < 0 || tRanY[0] < 0 ) {
      throw aol::Exception ( "ScalarArray<QC_2D>::fillValuesFrom: One or more of the ranges are incorrect: They must not have elements below 0!\n", __FILE__, __LINE__ );
    }

    // Check if source range is valid. This has to be done by the copier because different copiers have different prerequisites.
    if ( !getValue.sourceRangeOk ( source, sRanX, sRanY, tRanX, tRanY ) ) {
      throw aol::Exception ( "ScalarArray<QC_2D>::fillValuesFrom(): Source range invalid!\n", __FILE__, __LINE__ );
    }

    // Check if the target range is valid, i.e. not outside of this.
    if ( tRanX[1] > this->getNumX() ) {
      throw aol::Exception ( "ScalarArray<QC_2D>::fillValuesFrom(): X range does not fit into the current array!\n", __FILE__, __LINE__ );
    }
    if ( tRanY[1] > this->getNumY() ) {
      throw aol::Exception ( "ScalarArray<QC_2D>::fillValuesFrom(): Y range does not fit into the current array!\n", __FILE__, __LINE__ );
    }

    // Optionally check if the target range is as big as the source range.
    if ( ( ( tRanX[1] - tRanX[0] ) != ( sRanX[1] - sRanX[0] ) || ( tRanY[1] - tRanY[0] ) != ( sRanY[1] - sRanY[0] ) ) && !differentRangeSizesOk ) {
      throw aol::Exception ( "ScalarArray<QC_2D>::fillValuesFrom(): Source range and target range have different sizes!\n", __FILE__, __LINE__ );
    }


    // copy
    for ( int x = 0; x < ( tRanX[1] - tRanX[0] ); ++x )
      for ( int y = 0; y < ( tRanY[1] - tRanY[0] ); ++y ) {
        this->set ( tRanX[0] + x, tRanY[0] + y, getValue ( sRanX[0] + x, sRanY[0] + y, source ) );
      }
  }

  //! \brief Pastes image into current array. If the image is smaller in any direction, it is centered and padded with fillValue. If it is bigger, only its central values are copied.
  //! \author toelkes
  void padFrom ( const ScalarArray<DataType, qc::QC_2D> &image, const DataType FillValue = 0 ) {
    const int diffX = this->getNumX() - image.getNumX();
    const int diffY = this->getNumY() - image.getNumY();
    const int padX = diffX / 2;
    const int padY = diffY / 2;

    this->setAll ( FillValue );

    // If this->getNumX() - image.getNumX() (or getNumY respectively) is odd, ( this->getNumX() - image.getNumX() ) / 2 will be rounded down.
    // Therefore tRanX or sRanX (or tRanY or sRanY respectively) will be 1 too large.
    // To fix this the remainder ( this->getNumX() - image.getNumX() ) % 2 is added to/subtracted from padX at one side of the target or source range
    // resulting in a shift of the image (padY analog).
    const int tRanX[2] = { aol::Max ( aol::NumberTrait<int>::zero, padX ),
        this->getNumX() - aol::Max ( aol::NumberTrait<int>::zero, padX + diffX % 2 ) };
    const int tRanY[2] = { aol::Max ( aol::NumberTrait<int>::zero, padY ),
        this->getNumY() - aol::Max ( aol::NumberTrait<int>::zero, padY + diffY % 2 ) };
    const int sRanX[2] = { -aol::Min ( aol::NumberTrait<int>::zero, padX + diffX % 2 ),
        image.getNumX() + aol::Min ( aol::NumberTrait<int>::zero, padX ) };
    const int sRanY[2] = { -aol::Min ( aol::NumberTrait<int>::zero, padY + diffY % 2 ),
        image.getNumY() + aol::Min ( aol::NumberTrait<int>::zero, padY ) };

    StandardCopier copier;
    fillValuesFrom ( image, sRanX, sRanY, tRanX, tRanY, copier );
  }

  void flipFrom ( const ScalarArray<DataType, qc::QC_2D> &image, const qc::Comp component ) {
    // make sure there is enough memory
    if ( image.getSize() != this->getSize() )
      this->reallocate ( image );

    const int sRanX[2] = { 0, image.getNumX() };
    const int sRanY[2] = { 0, image.getNumY() };
    const int tRanX[2] = { 0, this->getNumX() };
    const int tRanY[2] = { 0, this->getNumY() };

    flipCopier copier ( this->getNumX(), this->getNumY(), component );
    fillValuesFrom ( image, sRanX, sRanY, tRanX, tRanY, copier );
  }

  void rotate90From ( const ScalarArray<DataType, qc::QC_2D> &image, const bool counterClockWise = true ) {
    // make sure the dimensions of this array match the dimensions of a 90 degree rotated version of the input image
    if ( this->getNumX ( ) != image.getNumY ( ) || this->getNumY ( ) != image.getNumX ( ) )
      this->reallocate ( image.getNumY ( ), image.getNumX ( ) );

    if ( counterClockWise ) {
      for ( int x=0; x<image.getNumX ( ) ; ++x )
        for ( int y=0; y<image.getNumY ( ) ; ++y )
          this->set ( y, image.getNumX ( ) - x - 1, image.get ( x, y ) );
    } else {
      for ( int x=0; x<image.getNumX ( ) ; ++x )
        for ( int y=0; y<image.getNumY ( ) ; ++y )
          this->set ( image.getNumY ( ) - y - 1, x, image.get ( x, y ) );
    }
  }

  void shiftByOffset ( const int xOffset, const int yOffset, const qc::ScalarArray<DataType, qc::QC_2D> &image ) {
    if ( ( this->getNumX() != image.getNumX() ) || ( this->getNumY() != image.getNumY() ) ) {
      throw aol::Exception ( "ScalarArray<QC_2D>::shiftByOffset(): Dimensions do not match!", __FILE__, __LINE__ );
    }

    const int tRanX[2] = { 0, this->getNumX() };
    const int tRanY[2] = { 0, this->getNumY() };
    const int sRanX[2] = { 0, image.getNumX() };
    const int sRanY[2] = { 0, image.getNumY() };

    OffsetCopier copier ( xOffset, yOffset, this->getNumX(), this->getNumY() );
    fillValuesFrom ( image, sRanX, sRanY, tRanX, tRanY, copier );
  }

  /** Shifts the contents of Image by the given offset
   *  \author wirth
   */
  void shiftByOffsetFrom ( const qc::CoordType Offset,
                           const qc::ScalarArray<DataType, qc::QC_2D> &Image ) {
    shiftByOffset ( Offset[0], Offset[1], Image );
  }

  /** Point-mirrors the image at the center ( numX/2, numY/2 ), i.e.:
   *  this( x , y ) = Image( numX - x - 1 , numY - y - 1 )
   *  \author Berkels
   */
  void pointMirrorAtCenter ( const qc::ScalarArray<DataType, qc::QC_2D> &Image ) {
    for ( int x = 0; x < this->numX; ++x ) {
      for ( int y = 0; y < this->numY; ++y ) {
        this->set ( x, y, Image.get ( this->numX - x - 1 , this->numY - y - 1 ) );
      }
    }
  }

  /** Point-mirrors the image at the origin ( 0, 0 ), compare pointMirrorAtCenter
   *  \author Berkels
   */
  void pointMirrorAtOrigin ( const qc::ScalarArray<DataType, qc::QC_2D> &Image ) {
    // basic idea: shift origin to center, mirror at center, shift center back to origin
    // The temporary array is not necessary, but make the code easier to understand and gives less code
    // duplication since pointMirrorAtCenter is used
    qc::ScalarArray<DataType, qc::QC_2D> temp ( Image, aol::STRUCT_COPY );

    // shift origin to center
    temp.shiftByOffset ( this->numX / 2, this->numY / 2, Image );
    // mirror at center
    this->pointMirrorAtCenter ( temp );
    temp = *this;
    // shift center back to origin: instead of shifting by -numX/2, -numY/2,
    // shift by numX - num/X2, instead, since shiftByOffset does not accept negative offsets
    // and shifting by numX, numY doesn't change the image
    this->shiftByOffset ( this->numX - this->numX / 2, this->numY - this->numY / 2, temp );
  }

  void copyBlockTo ( const int startX, const int startY, qc::ScalarArray<DataType, qc::QC_2D> &block ) const {
    const int sRanX[2] = { startX, startX + aol::Min( this->getNumX(), block.getNumX() ) };
    const int sRanY[2] = { startY, startY + aol::Min( this->getNumY(), block.getNumY() ) };
    const int tRanX[2] = { 0, aol::Min( this->getNumX(), block.getNumX() ) };
    const int tRanY[2] = { 0, aol::Min( this->getNumY(), block.getNumY() ) };

    StandardCopier copier;
    block.fillValuesFrom ( *this, sRanX, sRanY, tRanX, tRanY, copier );
  }

  void copyBlockTo ( const aol::Vec<2, short> Start, qc::ScalarArray<DataType, qc::QC_2D> &block ) const {
    copyBlockTo ( Start[0], Start[1], block );
  }

  void copyBlockTo ( const aol::Vec<2, int> Start, qc::ScalarArray<DataType, qc::QC_2D> &block ) const {
    copyBlockTo ( Start[0], Start[1], block );
  }

  void copyBlockTo ( const CoordType Start, qc::ScalarArray<DataType, qc::QC_2D> &block ) const {
    copyBlockTo ( Start[0], Start[1], block );
  }

  void setBlock ( const aol::Vec<2, int> &BlockStart, const aol::Vec<2, int> &BlockStop, const DataType Value )  {
    for ( int j = BlockStart[1]; j < BlockStop[1]; ++j )
      for ( int i = BlockStart[0]; i < BlockStop[0]; ++i )
        this->set ( i, j, Value );
  }

  void crop ( const aol::Vec<2, int> &CropStart, const aol::Vec<2, int> &CropSize )  {
    qc::ScalarArray<DataType, qc::QC_2D> copy ( *this, aol::DEEP_COPY );
    this->reallocate ( CropSize[0], CropSize[1] );
    copy.copyBlockTo ( CropStart, *this );
  }

  /**
   * \author Berkels
   */
  void cropBorders ( const DataType BorderColor = 0 )  {
    const int numX = this->getNumX();
    const int numY = this->getNumY();
    aol::Vec2<int> cropStart;
    aol::Vec2<int> cropStop ( numX, numY );
    for ( int i = 0; i < numX; ++i ){
      int j = 0;
      while ( j < numY && ( this->get( i, j ) == BorderColor ) )
        ++j;

      if ( j < numY ) {
        cropStart[0] = i;
        cropStop[0] -= cropStart[0];
        break;
      }
    }

    for ( int i = numX-1; i >= 0; --i ){
      int j = 0;
      while ( j < numY && ( this->get( i, j ) == BorderColor ) )
        ++j;

      if ( j < numY ) {
        cropStop[0] = i+1 - cropStart[0];
        break;
      }
    }

    for ( int j = 0; j < numY; ++j ){
      int i = 0;
      while ( i < numX && ( this->get( i, j ) == BorderColor ) )
        ++i;

      if ( i < numX ) {
        cropStart[1] = j;
        cropStop[1] -= cropStart[1];
        break;
      }
    }

    for ( int j = numY-1; j >= 0; --j ){
      int i = 0;
      while ( i < numX && ( this->get( i, j ) == BorderColor ) )
        ++i;

      if ( i < numX ) {
        cropStop[1] = j+1 - cropStart[1];
        break;
      }
    }
    crop ( cropStart, cropStop );
  }

  /**
   * Copies the data of the array into the buffer while flipping the x and y indexing.
   *
   * \author Berkels
   */
  void copyToBufferFlipped ( DataType *buffer ) const {
    for ( int i = 0; i < this->numX; ++i )
      for ( int j = 0; j < this->numY; ++j )
        buffer[i * this->numY + j] = this->get ( i, j );
  }

  /**
   * Reads the data from the buffer into the array while flipping the x and y indexing.
   *
   * \author Berkels
   */
  void readFromBufferFlipped ( DataType *buffer ) {
    for ( int i = 0; i < this->numX; ++i )
      for ( int j = 0; j < this->numY; ++j )
        this->set ( i, j , buffer[i * this->numY + j] );
  }

  /** interpret matrix entries as grey values
   */
  template< typename MatrixType >
  void visualizeMatrix ( const MatrixType &mat ) {
    reallocate ( mat.getNumCols (), mat.getNumRows () );
    this->setOverflowHandling ( aol::SCALE, 0, aol::OverflowTrait<DataType>::max );

    for ( int i = 0; i < mat.getNumRows (); ++i ) {
      for ( int j = 0; j < mat.getNumCols (); ++j ) {
        set ( j, i, mat.get ( i, j ) );
      }
    }
  }

  //! return number of points in all space dimensions, if equal, throw exception otherwise
  //! \todo rename
  int getNumXYZ ( ) const {
    if ( this->getNumX() != this->getNumY() )
      throw aol::Exception ( "qc::ScalarArray<DataType, QC_2D>::getNumXYZ() may not be called for non-square arrays", __FILE__, __LINE__ );

    return ( this->getNumX() );
  }


  //! truncates bottom row and right column and writes results to image
  void truncateAndWriteTo ( aol::Vector<DataType> &image, int gridlevel ) const {
    int TwoToLevPlusOne = aol::Pow ( 2, gridlevel ) + 1;
    if ( image.size() != aol::Pow ( TwoToLevPlusOne - 1, 2 ) )
      throw aol::Exception ( "ScalarArray<_DataType, qc::QC_2D> ::truncateAndWriteTo: size of image does not fit to gridlevel.", __FILE__, __LINE__ );
    if ( this->size() != aol::Pow ( TwoToLevPlusOne,   2 ) )
      throw aol::Exception ( "ScalarArray<_DataType, qc::QC_2D> ::truncateAndWriteTo: *this has wrong size.", __FILE__, __LINE__ );

    int j = 0;
    for ( int i = 1 + TwoToLevPlusOne; i < this->size() + 1; i++ )
      if ( i % TwoToLevPlusOne != 0 )
        image[j++] = this->get ( i - 1 );
  }

  //! loads from truncated image and adds zero to bottom row and right column
  void loadFromTruncatedImage (  const aol::Vector<DataType> &image, int gridlevel ) {
    int TwoToLevPlusOne = aol::Pow ( 2, gridlevel ) + 1;
    if ( image.size() != aol::Pow ( TwoToLevPlusOne - 1, 2 ) )
      throw aol::Exception ( "ScalarArray<_DataType, qc::QC_2D> ::loadFromTruncatedImage: size of image does not fit to gridlevel.", __FILE__, __LINE__ );
    if ( this->size() != aol::Pow ( TwoToLevPlusOne,   2 ) )
      reallocate ( TwoToLevPlusOne, TwoToLevPlusOne );

    // add zero row first
    for ( int i = 0; i < TwoToLevPlusOne; i++ )
      this->set ( i, 0 );

    int j = 0;
    for ( int i = TwoToLevPlusOne + 1; i < this->size() + 1; i++ )
      if ( i % TwoToLevPlusOne != 0  )
        this->set ( i - 1, image[j++] );
      else
        this->set ( i - 1, 0 );

  }

  /**
   * \brief \f$ \text{arg}[j] =\sum_i \text{get}(i,j) \f$
   *
   * \author Berkels
   */
  void sumInXDirectionTo ( aol::Vector<DataType> &DirectionSums ) const {
    DirectionSums.setZero();
    for ( int j = 0; j < this->numY; ++j ) {
      for ( int i = 0; i < this->numX; ++i ) {
        DirectionSums.add ( j, this->get ( i, j ) );
      }
    }
  }

  /**
   * \brief \f$ \text{this}(i,j) = \text{arg}[j] \f$
   *
   * \author Berkels
   */
  void setInXDirectionFrom ( const aol::Vector<DataType> &ArgVec ) {
    for ( int j = 0; j < this->numY; ++j ) {
      const DataType value = ArgVec.get ( j );
      for ( int i = 0; i < this->numX; ++i ) {
        this->set ( i, j, value );
      }
    }
  }
  
  /**
   * \brief \f$ \text{arg}[i] =\sum_j \text{get}(i,j) \f$
   *
   * \author Tatano
   */
  void sumInYDirectionTo ( aol::Vector<DataType> &DirectionSums ) const {
    DirectionSums.setZero();
    for ( int i = 0; i < this->numX; ++i ) {
      for ( int j = 0; j < this->numY; ++j ) {
        DirectionSums.add ( i, this->get ( i, j ) );
      }
    }
  }

  /**
   * Checks whether the values of all neighboring pixels are smaller or equal to the value of the specified pixel.
   *
   * \author Berkels
   */
  bool isPositionLocalMaximum ( const int X, const int Y ) const {
    const int minX = aol::Max ( X - 1, 0 );
    const int minY = aol::Max ( Y - 1, 0 );
    const int maxX = aol::Min ( X + 1, this->getNumX() );
    const int maxY = aol::Min ( Y + 1, this->getNumY() );

    const DataType value = this->get ( X, Y );
    for ( int y = minY; y < maxY; ++y ) {
      for ( int x = minX; x < maxX; ++x ) {
        if ( this->get ( x, y ) > value )
          return false;
      }
    }
    return true;
  }

  const GaussDiffKernel2d<RealType> &getDiffKernelXXReference ( ) const {
    return *diffKernelXX;
  }

  const GaussDiffKernel2d<RealType> &getDiffKernelYYReference ( ) const {
    return *diffKernelYY;
  }

private:
  void resize ( const int ) {
    throw aol::Exception ( "qc::Array::resize( int ) preserving contents does not make sense.", __FILE__, __LINE__ );
  }

  void reallocate ( const int ) {
    throw aol::Exception ( "qc::Array::reallocate( int ) preserving contents does not make sense.", __FILE__, __LINE__ );
  }


protected:

  GaussDiffKernel2d<RealType> *diffKernelX, *diffKernelY,
                    *diffKernelXX, *diffKernelXY, *diffKernelYY;

  static const aol::Format &format;              //! for printing
  static bool prettyFormat;

};

template <typename _DataType>
ostream &operator<< ( ostream &os, const ScalarArray<_DataType, qc::QC_2D> &s ) {
  return s.print ( os );
}

/**
 * A full 3d array implementation for sizes \f$(2^n+k)^3\quad k=0,1\f$.
 *
 * DataType specifies the type of the data to be stored. A corresponding
 * RealType is selected by the base class aol::Vector and used in computations
 * such as L2-projections, curvature, etc. if DataType is an int for instance.
 **/
template <typename _DataType>
class ScalarArray<_DataType, qc::QC_3D> : public qc::Array<_DataType> {
public:
  typedef _DataType DataType;
  typedef typename aol::Vector<DataType>::RealType RealType;
  static const qc::Dimension Dim = qc::QC_3D;
protected:
  typedef RealType ( *FUNC ) ( RealType, RealType, RealType );
  FUNC polyBases[ 10 ];

  static RealType _xx ( RealType X, RealType  , RealType ) {
    return X * X;
  }
  static RealType _yy ( RealType  , RealType Y, RealType ) {
    return Y * Y;
  }
  static RealType _zz ( RealType  , RealType  , RealType Z ) {
    return Z * Z;
  }
  static RealType _yz ( RealType  , RealType Y, RealType Z ) {
    return Y * Z;
  }
  static RealType _xz ( RealType X, RealType  , RealType Z ) {
    return X * Z;
  }
  static RealType _xy ( RealType X, RealType Y, RealType ) {
    return X * Y;
  }
  static RealType _x ( RealType X, RealType  , RealType ) {
    return X;
  }
  static RealType _y ( RealType  , RealType Y, RealType ) {
    return Y;
  }
  static RealType _z ( RealType  , RealType  , RealType Z ) {
    return Z;
  }
  static RealType _1 ( RealType  , RealType  , RealType ) {
    return 1.0;
  }

  void initPolyBases ( void ) {
    polyBases[ 0 ] = _xx;
    polyBases[ 1 ] = _yy;
    polyBases[ 2 ] = _zz;
    polyBases[ 3 ] = _yz;
    polyBases[ 4 ] = _xz;
    polyBases[ 5 ] = _xy;
    polyBases[ 6 ] = _x;
    polyBases[ 7 ] = _y;
    polyBases[ 8 ] = _z;
    polyBases[ 9 ] = _1;
  }


public:
  class Hexahedron {
  public:
    Hexahedron ( const Array<DataType> &Array, int X, int Y, int Z, int Step )
      : array ( Array ) {
      v000 = static_cast<RealType> ( array.get ( X       , Y       , Z ) );
      v100 = static_cast<RealType> ( array.get ( X + Step, Y       , Z ) );
      v010 = static_cast<RealType> ( array.get ( X       , Y + Step, Z ) );
      v110 = static_cast<RealType> ( array.get ( X + Step, Y + Step, Z ) );
      v001 = static_cast<RealType> ( array.get ( X       , Y       , Z + Step ) );
      v101 = static_cast<RealType> ( array.get ( X + Step, Y       , Z + Step ) );
      v011 = static_cast<RealType> ( array.get ( X       , Y + Step, Z + Step ) );
      v111 = static_cast<RealType> ( array.get ( X + Step, Y + Step, Z + Step ) );
    }

    // values should be in [0,1]
    // barycentric coordinates
    RealType localEvaluate ( RealType a1, RealType a2, RealType a3 ) const {
      return ( 1.0f - a3 ) * ( ( 1.0f - a2 ) * ( ( 1.0f - a1 ) * v000 + a1 * v100 )
                               + a2        * ( ( 1.0f - a1 ) * v010 + a1 * v110 ) )
             +  a3 * ( ( 1.0f - a2 ) * ( ( 1.0f - a1 ) * v001 + a1 * v101 )
                       + a2        * ( ( 1.0f - a1 ) * v011 + a1 * v111 ) );
    }

    RealType boundaryLocalEvaluate ( RealType a1, RealType a2, RealType a3 ) const {
      if ( a1 == 0.0f )
        return ( 1.0f - a3 ) * ( ( 1.0f - a2 ) * v000 + a2 * v010 )
               +  a3 * ( ( 1.0f - a2 ) * v001  + a2 * v011 );
      else if ( a1 == 1.0f )
        return ( 1.0f - a3 ) * ( ( 1.0f - a2 ) * v100 + a2        * v110 )
               +  a3 * ( ( 1.0f - a2 ) * v101 + a2 * v111 );

      if ( a2 == 0.0f )
        return ( 1.0f - a3 ) * ( ( 1.0f - a1 ) * v000 + a1 * v100 )
               +  a3 * ( ( 1.0f - a1 ) * v001 + a1 * v101 );
      else if ( a2 == 1.0f )
        return ( 1.0f - a3 ) * ( ( 1.0f - a1 ) * v010 + a1 * v110 )
               +  a3 * ( ( 1.0f - a1 ) * v011 + a1 * v111 );

      if ( a3 == 0.0f )
        return ( ( 1.0f - a2 ) * ( ( 1.0f - a1 ) * v000 + a1 * v100 )
                 + a2        * ( ( 1.0f - a1 ) * v010 + a1 * v110 ) );
      else if ( a3 == 1.0f )
        return ( ( 1.0f - a2 ) * ( ( 1.0f - a1 ) * v001 + a1 * v101 )
                 + a2        * ( ( 1.0f - a1 ) * v011 + a1 * v111 ) );

      return localEvaluate ( a1, a2, a3 );

    }

    Hexahedron &operator-= ( RealType Val ) {
      v000 -= Val;
      v100 -= Val;
      v010 -= Val;
      v110 -= Val;
      v001 -= Val;
      v101 -= Val;
      v011 -= Val;
      v111 -= Val;
      return *this;
    }

  protected:
    const Array<DataType> &array;
    RealType v000, v100, v010, v110,
             v001, v101, v011, v111;
  };

  // --------Functor classes used by fillValuesFrom---------------------

  // \brief Standart copying class
  // \author toelkes
  class StandardCopier {
  public:
    inline bool sourceRangeOk ( const ScalarArray<DataType, qc::QC_3D> &source, const int sRanX[2], const int sRanY[2], const int sRanZ[2],
                                const int tRanX[2], const int tRanY[2], const int tRanZ[2] ) {
      // Check if only valid values of source are requested in getValuesFrom().
      return ( sRanX[0] + ( tRanX[1] - tRanX[0] ) <= source.getNumX()
               && sRanY[0] + ( tRanY[1] - tRanY[0] ) <= source.getNumY()
               && sRanZ[0] + ( tRanZ[1] - tRanZ[0] ) <= source.getNumZ() );
    }

    inline DataType operator() ( int x, int y, int z, const ScalarArray<DataType, qc::QC_3D> &source ) {
      return source.get ( x, y, z );
    }
  };

  // \brief Class used in shiftByOffset
  // \author toelkes
  class OffsetCopier {
  public:
    OffsetCopier ( const int xOffset, const int yOffset, const int zOffset, const int numX, const int numY, const int numZ )
      : _xOffset ( xOffset ), _yOffset ( yOffset ), _zOffset ( zOffset ),  _numX ( numX ), _numY ( numY ), _numZ ( numZ ) {
      if ( abs ( _xOffset ) > _numX || abs ( _yOffset ) > _numY || abs ( _zOffset ) > _numZ ) {
        throw aol::Exception ( "OffsetCopier: Offset larger than array!\n", __FILE__, __LINE__ );
      }
    }

    inline bool sourceRangeOk ( const ScalarArray<DataType, qc::QC_3D> &source, const int sRanX[2], const int sRanY[2], const int sRanZ[2],
                                const int tRanX[2], const int tRanY[2], const int tRanZ[2] ) {
      // Check if the offsets are valid. Is checked before, but if the exception is catched this could still be called.
      if ( abs ( _xOffset ) > _numX || abs ( _yOffset ) > _numY || abs ( _zOffset ) > _numZ ) {
        throw aol::Exception ( "OffsetCopier: Cannot use OffsetCopier with offsets bigger than the grid size!\n", __FILE__, __LINE__ );
      }

      // Check if only valid values of source are requested in getValuesFrom().
      return ( sRanX[0] + ( tRanX[1] - tRanX[0] ) <= source.getNumX()
               && sRanY[0] + ( tRanY[1] - tRanY[0] ) <= source.getNumY()
               && sRanZ[0] + ( tRanZ[1] - tRanZ[0] ) <= source.getNumZ() );
    }

    inline DataType operator() ( int x, int y, int z, const ScalarArray<DataType, qc::QC_3D> &source ) {
      return source.get ( ( ( x - _xOffset ) >= 0 ) ? ( x - _xOffset ) % _numX : ( _numX + ( x - _xOffset ) ),
                          ( ( y - _yOffset ) >= 0 ) ? ( y - _yOffset ) % _numY : ( _numY + ( y - _yOffset ) ),
                          ( ( z - _zOffset ) >= 0 ) ? ( z - _zOffset ) % _numZ : ( _numZ + ( z - _zOffset ) ) );
    }

  private:
    const int _xOffset, _yOffset, _zOffset, _numX, _numY, _numZ;
  };

  // \brief Class used in flipFrom. Copies and flips.
  // \author toelkes
  class flipCopier {
  public:
    flipCopier ( const int numX, const int numY, const int numZ, const qc::Comp component )
      : _numX ( numX ), _numY ( numY ), _numZ ( numZ ), _component ( component )
    {}

    inline bool sourceRangeOk ( const ScalarArray<DataType, qc::QC_3D> &source, const int sRanX[2], const int sRanY[2], const int sRanZ[2],
                                const int tRanX[2], const int tRanY[2], const int tRanZ[2] ) {
      // Check if only valid values of source are requested in getValuesFrom().
      return ( sRanX[0] + ( tRanX[1] - tRanX[0] ) <= source.getNumX()
               && sRanY[0] + ( tRanY[1] - tRanY[0] ) <= source.getNumY()
               && sRanZ[0] + ( tRanZ[1] - tRanZ[0] ) <= source.getNumZ() );
    }

    inline DataType operator() ( int x, int y, int z, const ScalarArray<DataType, qc::QC_3D> &source ) {
      if ( _component == qc::QC_X )
        return source.get ( _numX - 1 - x, y, z );
      else if ( _component == qc::QC_Y )
        return source.get ( x, _numY - 1 - y, z );
      else if ( _component == qc::QC_Z )
        return source.get ( x, y, _numZ - 1 - z );
      else {
        throw aol::Exception ( "flipCopier: Component to flip is invalid!\n", __FILE__, __LINE__ );
        return source.get ( x, y, z );
      }
    }

  private:
    const int _numX, _numY, _numZ;
    const qc::Comp _component;
  };

  // \brief Class used in resampleFrom. Copies and scales/interpolates.
  // \author toelkes
  class ResampleCopier {
  public:
    ResampleCopier ( const ScalarArray<DataType, qc::QC_3D> &source, const ScalarArray<DataType, qc::QC_3D> &target ) {
      _hx = static_cast<RealType> ( source.getNumX() - 1 ) / ( target.getNumX() - 1 );
      _hy = static_cast<RealType> ( source.getNumY() - 1 ) / ( target.getNumY() - 1 );
      _hz = static_cast<RealType> ( source.getNumZ() - 1 ) / ( target.getNumZ() - 1 );
    }

    inline bool sourceRangeOk ( const ScalarArray<DataType, qc::QC_3D> &source, const int sRanX[2], const int sRanY[2], const int sRanZ[2],
                                const int tRanX[2], const int tRanY[2], const int tRanZ[2] ) {
      // Only values inside of source can be interpolated.
      return ( static_cast<int>(( sRanX[0] + ( tRanX[1] - tRanX[0] ) - 1 ) * _hx) <= ( source.getNumX() - 1 )
               &&  static_cast<int> (( sRanY[0] + ( tRanY[1] - tRanY[0] ) - 1 ) * _hy) <= ( source.getNumY() - 1 )
               && static_cast<int>(( sRanZ[0] + ( tRanZ[1] - tRanZ[0] ) - 1 ) * _hz) <= ( source.getNumZ() - 1 ) );
    }

    inline DataType operator() ( int x, int y, int z, const ScalarArray<DataType, qc::QC_3D> &source ) {
      const RealType scaledX = x * _hx, scaledY = y * _hy, scaledZ = z * _hz;
      return static_cast<DataType> ( source.interpolate ( scaledX, scaledY, scaledZ ) );
    }

  private:
    // long double precision is crucial here.
    long double _hx, _hy, _hz;
  };

  //--------------------------------------------------------------

  ScalarArray () :
    Array<DataType> ( 0, 0, 0 ) {
    initPolyBases();
    init();
  }

  /*! the array can be given a pointer with inverse-lexicographically ordered data
      which will not be deleted in the destructor of the array
    */
  ScalarArray ( int NumX, int NumY, int NumZ, DataType *Data, aol::CopyFlag copyFlag = aol::FLAT_COPY ) :
    Array<DataType> ( NumX, NumY, NumZ, Data, copyFlag ) { // default: flat copy!
    initPolyBases();
    init();
  }

  ScalarArray ( DataType *Data, const GridSize<QC_3D> &Size, aol::CopyFlag copyFlag = aol::FLAT_COPY ) :
    Array<DataType> ( Size.getNumX(), Size.getNumY(), Size.getNumZ(), Data, copyFlag ) {
    initPolyBases();
    init();
  }

  /*! Generating an array of dimension \f$NumX\cdot NumY\cdot NumZ\f$.
   */
  ScalarArray ( int NumX, int NumY, int NumZ ) :
    Array<DataType> ( NumX, NumY, NumZ ) {
    initPolyBases();
    init();
  }

  explicit ScalarArray ( const aol::Vec3<int> &Size ) :
    Array<DataType> ( Size[0], Size[1], Size[2] ) {
    initPolyBases();
    init();
  }

  explicit ScalarArray ( int Width ) :
    Array<DataType> ( Width, Width, Width ) {
    initPolyBases();
    init();
  }

  explicit ScalarArray ( const GridStructure &Grid ) :
    Array<DataType> ( Grid.getNumX(), Grid.getNumY(), Grid.getNumZ() ) {
    initPolyBases();
    if ( Grid.getDimOfWorld() != 3 ) {
      throw aol::Exception ( "ScalarArray<QC_3D> initialized with a grid whose dimension does not equal 3!\n", __FILE__, __LINE__ );
    }
    init();
  }

  explicit ScalarArray ( const GridSize<QC_3D> &GridSize ) :
    Array<DataType> ( GridSize.getNumX(), GridSize.getNumY(), GridSize.getNumZ() ) {
    initPolyBases();
    init();
  }

  ScalarArray ( const aol::Vector<DataType> &Vector, int NumX, int NumY, int NumZ, aol::CopyFlag copyFlag = aol::FLAT_COPY ) :
    Array<DataType> ( Vector, NumX, NumY, NumZ, copyFlag ) { // default: flat copy

    if ( Vector.size() != this->getNumX() *this->getNumY() *this->getNumZ() ) {
      char error[1024];
      sprintf ( error, "qc::ScalarArray<QC_3D>::ScalarArray( const aol::Vector<DataType> &, int, int, int ): Vectorlength = %d should be equal to size of array, which is %d", Vector.size(), this->getNumX() *this->getNumY() *this->getNumZ() );
      throw ( aol::Exception ( error, __FILE__, __LINE__ ) );
    }
    init();
  }

  template <typename GridType>
  ScalarArray ( const aol::Vector<DataType> &Vector, const GridType &Grid, aol::CopyFlag copyFlag = aol::FLAT_COPY ) :
    Array<DataType> ( Vector, Grid.getNumX(), Grid.getNumY(), Grid.getNumZ(), copyFlag ) { // default: flat copy

    if ( Vector.size() != this->getNumX() *this->getNumY() *this->getNumZ() ) {
      char error[1024];
      sprintf ( error, "qc::ScalarArray<QC_3D>::ScalarArray( const aol::Vector<DataType> &, const GridType ): Vectorlength = %d should be equal to size of array, which is %d", Vector.size(), this->getNumX() *this->getNumY() *this->getNumZ() );
      throw ( aol::Exception ( error, __FILE__, __LINE__ ) );
    }
    init();
  }

  /** Copy constructor
   */
  explicit ScalarArray ( const ScalarArray<DataType, qc::QC_3D> &org, aol::CopyFlag copyFlag = aol::DEEP_COPY ) :
    Array<DataType> ( org, copyFlag ),
    diffKernelX ( org.diffKernelX ? new GaussDiffKernel3d<RealType> ( *org.diffKernelX ) : NULL ),
    diffKernelY ( org.diffKernelY ? new GaussDiffKernel3d<RealType> ( *org.diffKernelY ) : NULL ),
    diffKernelZ ( org.diffKernelZ ? new GaussDiffKernel3d<RealType> ( *org.diffKernelZ ) : NULL ) {
    initPolyBases();
    // median sort vec has already correct size
  }

  ScalarArray ( const Array<DataType> &A, aol::CopyFlag copyFlag ) :
    Array<DataType> ( A, A.getNumX(), A.getNumY(), A.getNumZ(), copyFlag ) {
    initPolyBases();
    init();
  }

  //! Read data from file
  //! and set size of array correctly automatically
  explicit ScalarArray ( const string &filename );

  virtual ~ScalarArray();

  void reallocate ( const int NumX, const int NumY, const int NumZ ) {
    Array<DataType>::reallocate ( NumX, NumY, NumZ );
    this->numX = NumX;
    this->numY = NumY;
    this->numZ = NumZ;
  }

  void reallocate ( const qc::CoordType &size ) {
    reallocate ( size[0], size[1], size[2] );
  }

  template< typename Structure >
  void reallocate ( const Structure &other ) {
    reallocate ( other.getNumX(), other.getNumY(), other.getNumZ() );
  }

  //! Initialize some variables, do not touch data!
  void init() {
    diffKernelX = new GaussDiffKernel3d<RealType> ( this->DIFF_STENCIL, this->DIFF_STD_SIGMA, DIFF_X );
    diffKernelY = new GaussDiffKernel3d<RealType> ( this->DIFF_STENCIL, this->DIFF_STD_SIGMA, DIFF_Y );
    diffKernelZ = new GaussDiffKernel3d<RealType> ( this->DIFF_STENCIL, this->DIFF_STD_SIGMA, DIFF_Z );
  }

  /** Assignment of an vector as data, retain array dimensions etc.
   */
  ScalarArray<DataType, qc::QC_3D> & operator= ( const aol::Vector<DataType> &V ) {
    aol::Vector<DataType>::operator= ( V );
    return *this;
  }

  /** Assignment
   */
  ScalarArray<DataType, qc::QC_3D> &operator= ( const ScalarArray<DataType, qc::QC_3D>& V ) {

    if ( this == &V ) { // correct self-assignment
      return ( *this );
    }

    if ( this->numX != V.numX ||  this->numY != V.numY || this->numZ != V.numZ ) {
      throw aol::Exception ( "qc::ScalarArray<QC_3D>::operator= trying to assign incompatible array", __FILE__, __LINE__ );
    }

    operator= ( static_cast<const aol::Vector<DataType>& > ( V ) );

    if ( diffKernelX ) delete diffKernelX;
    if ( diffKernelY ) delete diffKernelY;
    if ( diffKernelZ ) delete diffKernelZ;

    this->_offset = V._offset;

    diffKernelX = ( V.diffKernelX ? new GaussDiffKernel3d<RealType> ( *V.diffKernelX ) : NULL );
    diffKernelY = ( V.diffKernelY ? new GaussDiffKernel3d<RealType> ( *V.diffKernelY ) : NULL );
    diffKernelZ = ( V.diffKernelZ ? new GaussDiffKernel3d<RealType> ( *V.diffKernelZ ) : NULL );
    initPolyBases();

    return *this;
  }

  /** Returns an interpolated value of the image at position X, Y, Z
     * @param X x-coordinate of the interpolated point
     * @param Y y-coordinate of the interpolated point
     * @param Z z-coordinate of the interpolated point
     */
  RealType interpolate ( RealType X, RealType Y, RealType Z ) const {

    int xL = static_cast<int> ( X ), yL = static_cast<int> ( Y ), zL = static_cast<int> ( Z );

    if ( X == static_cast<RealType> ( this->getNumX() - 1 ) ) xL -= 1;
    if ( Y == static_cast<RealType> ( this->getNumY() - 1 ) ) yL -= 1;
    if ( Z == static_cast<RealType> ( this->getNumZ() - 1 ) ) zL -= 1;

    X -= xL;
    Y -= yL;
    Z -= zL;

    const int xU = xL + 1, yU = yL + 1, zU = zL + 1;

    RealType v[8] = { static_cast<RealType> ( this->get ( xL, yL, zL ) ),
                      static_cast<RealType> ( this->get ( xU, yL, zL ) ),
                      static_cast<RealType> ( this->get ( xL, yL, zU ) ),
                      static_cast<RealType> ( this->get ( xU, yL, zU ) ),
                      static_cast<RealType> ( this->get ( xL, yU, zL ) ),
                      static_cast<RealType> ( this->get ( xU, yU, zL ) ),
                      static_cast<RealType> ( this->get ( xL, yU, zU ) ),
                      static_cast<RealType> ( this->get ( xU, yU, zU ) )
                    };

    const RealType oneMinusY = 1 - Y;
    const RealType oneMinusZ = 1 - Z;
    return (  ( 1 - X ) * ( oneMinusY * ( oneMinusZ * v[0] + Z * v[2] )
                            + Y * ( oneMinusZ * v[4] + Z * v[6] ) )
              + X * ( oneMinusY * ( oneMinusZ * v[1] + Z * v[3] )
                      + Y * ( oneMinusZ * v[5] + Z * v[7] ) ) );

    //return ( (1-X)*(1-Y)*(1-Z)*v[0] + X*(1-Y)*(1-Z)*v[1] + (1-X)*(1-Y)*Z*v[2] + X*(1-Y)*Z*v[3] + (1-X)*Y*(1-Z)*v[4] + X*Y*(1-Z)*v[5] + (1-X)*Y*Z*v[6] + X*Y*Z*v[7] );
  }

  /** Returns an interpolated value at position x, y, z (world coordinates) where the array is assumed to discretize [0,1]^3
   * Throws if the scalar array is not cubic
   * @param x x-coordinate of interpolated point (in [0,1])
   * @param y y-coordinate of interpolated point (in [0,1])
   * @param z z-coordinate of interpolated point (in [0,1])
   */
  RealType interpolate_on01 ( RealType X, RealType Y, RealType Z ) const {
#ifdef BOUNDS_CHECK
    // assert that dimensions are equal
    GridSize<QC_3D> sizeChecker ( *this );
    sizeChecker.quadraticOrDie ();
#endif
    return interpolate ( X * ( this->getNumX() - 1 ), Y * ( this->getNumY() - 1 ), Z * ( this->getNumZ() - 1 ) );
  }

  RealType interpolate_on01_periodic ( RealType X, RealType Y, RealType Z ) const {
    return interpolate_on01 ( fabs ( X - floor ( X ) ), fabs ( Y - floor ( Y ) ), fabs ( Z - floor ( Z ) ) );
  }

  /** @returns an interpolated value of the image at position Coords
   *  @param Coord coordinates of the interpolated point
   */
  RealType interpolate ( const aol::Vec<3, RealType> &Coord ) const {
    return interpolate ( Coord[0], Coord[1], Coord[2] );
  }

  RealType interpolate_on01 ( const aol::Vec3<RealType> & coord ) const {
    return interpolate_on01 ( coord[0], coord[1], coord[2] );
  }

  RealType interpolate_on01_periodic ( const aol::Vec3<RealType> & coord ) const {
    return interpolate_on01_periodic ( coord[0], coord[1], coord[2] );
  }

  /** Interpolation on a fine grid cell.
   * @returns interpolated value on the barycenter of the element
   * @param El Element
   */
  RealType interpolate ( const Element &El ) const {
    const short &x = El.xrefc();
    const short &y = El.yrefc();
    const short &z = El.zrefc();
    return static_cast<RealType> ( 0.125 * ( this->get ( x, y, z ) + this->get ( x + 1, y, z ) + this->get ( x, y + 1, z ) + this->get ( x + 1, y + 1, z ) +
                                             this->get ( x, y, z + 1 ) + this->get ( x + 1, y, z + 1 ) + this->get ( x, y + 1, z + 1 ) + this->get ( x + 1, y + 1, z + 1 ) ) );
  }

  void gradient ( const aol::Vec3<RealType> &coord, aol::Vec3<RealType> &Grad ) const {
    Grad[0] = dx ( coord );
    Grad[1] = dy ( coord );
    Grad[2] = dz ( coord );
  }

  /** Copmute the whole gradient in the cell containing coord with finite
   *  differences and store the result in Grad (assuming voxel size 1)
   */
  void gradientFD ( const aol::Vec3<DataType> &coord, aol::Vec3<DataType> &Grad ) const {
    Grad[0] = static_cast<DataType> ( dxFD ( static_cast<int> ( coord[0] ), static_cast<int> ( coord[1] ), static_cast<int> ( coord[2] ) ) );
    Grad[1] = static_cast<DataType> ( dyFD ( static_cast<int> ( coord[0] ), static_cast<int> ( coord[1] ), static_cast<int> ( coord[2] ) ) );
    Grad[2] = static_cast<DataType> ( dzFD ( static_cast<int> ( coord[0] ), static_cast<int> ( coord[1] ), static_cast<int> ( coord[2] ) ) );
  }

  /** Copmute the whole gradient in the cell containing coord with finite
   *  differences and store the result in Grad (assuming voxel size 1)
   */
  void gradientFD ( const int x, const int y, const int z, aol::Vec3<DataType> &Grad ) const {
    Grad[0] = static_cast<DataType> ( dxFD ( x, y, z ) );
    Grad[1] = static_cast<DataType> ( dyFD ( x, y, z ) );
    Grad[2] = static_cast<DataType> ( dzFD ( x, y, z ) );
  }

  void cellGradient ( const Element &el, aol::Vec3<RealType> &Grad ) const {
    Grad[ 0 ] = dx ( el );
    Grad[ 1 ] = dy ( el );
    Grad[ 2 ] = dz ( el );
  }

  /** Returns the point wise finite difference in x direction (assuming voxel size 1)
   */
  RealType dxFD ( const int X, const int Y, const int Z ) const {
    if ( X == 0 ) {
      return static_cast<RealType> ( this->get ( X + 1, Y, Z ) - this->get ( X, Y, Z ) );
    } else if ( X == this->numX - 1 ) {
      return static_cast<RealType> ( this->get ( X, Y, Z ) - this->get ( X - 1, Y, Z ) );
    } else {
      return static_cast<RealType> ( 0.5 * ( this->get ( X + 1, Y, Z ) - this->get ( X - 1, Y, Z ) ) );
    }
  }

  RealType dxFD ( const CoordType &Coord ) const {
    return dxFD ( Coord[0], Coord[1], Coord[2] );
  }

  /** Returns the point wise finite difference in y direction (assuming voxel size 1)
  */
  RealType dyFD ( const int X, const int Y, const int Z ) const {
    if ( Y == 0 ) {
      return static_cast<RealType> ( this->get ( X, Y + 1, Z ) - this->get ( X, Y, Z ) );
    } else if ( Y == this->numY - 1 ) {
      return static_cast<RealType> ( this->get ( X, Y, Z ) - this->get ( X, Y - 1, Z ) );
    } else {
      return static_cast<RealType> ( 0.5 * ( this->get ( X, Y + 1, Z ) - this->get ( X, Y - 1, Z ) ) );
    }
  }

  /** Returns the point wise finite difference in z direction (assuming voxel size 1)
   */
  RealType dzFD ( const int X, const int Y, const int Z ) const {
    if ( Z == 0 ) {
      return static_cast<RealType> ( this->get ( X, Y, Z + 1 ) - this->get ( X, Y, Z ) );
    } else if ( Z == this->numZ - 1 ) {
      return static_cast<RealType> ( this->get ( X, Y, Z ) - this->get ( X, Y, Z - 1 ) );
    } else {
      return static_cast<RealType> ( 0.5 * ( this->get ( X, Y, Z + 1 ) - this->get ( X, Y, Z - 1 ) ) );
    }
  }


  /** Returns the element wise central x - derivative of the image in
   *  element e
   *  @param e The element on which computation is to be performed
   */
  RealType dx ( const Element &e ) const {
    int x = e.x(), y = e.y(), z = e.z();

    if ( x == this->getNumX() - 1 ) --x;
    if ( y == this->getNumY() - 1 ) --y;
    if ( z == this->getNumZ() - 1 ) --z;

    return static_cast<RealType> ( ( this->get ( x + 1, y  , z ) - this->get ( x, y  , z ) +
                                     this->get ( x + 1, y  , z + 1 ) - this->get ( x, y  , z + 1 ) +
                                     this->get ( x + 1, y + 1, z ) - this->get ( x, y + 1, z ) +
                                     this->get ( x + 1, y + 1, z + 1 ) - this->get ( x, y + 1, z + 1 ) ) / 4. );
  }

  /** Returns the element wise central y - derivative of the image in
   *  element e
   *  @param e The element on which computation is to be performed
   */
  RealType dy ( const Element &e ) const {
    int x = e.x(), y = e.y(), z = e.z();

    if ( x == this->getNumX() - 1 ) --x;
    if ( y == this->getNumY() - 1 ) --y;
    if ( z == this->getNumZ() - 1 ) --z;

    return static_cast<RealType> ( ( this->get ( x  , y + 1, z ) - this->get ( x  , y, z ) +
                                     this->get ( x  , y + 1, z + 1 ) - this->get ( x  , y, z + 1 ) +
                                     this->get ( x + 1, y + 1, z ) - this->get ( x + 1, y, z ) +
                                     this->get ( x + 1, y + 1, z + 1 ) - this->get ( x + 1, y, z + 1 ) ) / 4. );
  }

  /** Returns the element wise central z - derivative of the image in
   *  element e
   *  @param e The element on which computation is to be performed
   */
  RealType dz ( const Element &e ) const {
    int x = e.x(), y = e.y(), z = e.z();

    if ( x == this->getNumX() - 1 ) --x;
    if ( y == this->getNumY() - 1 ) --y;
    if ( z == this->getNumZ() - 1 ) --z;

    return static_cast<RealType> ( ( this->get ( x  , y  , z + 1 ) - this->get ( x  , y  , z ) +
                                     this->get ( x + 1, y  , z + 1 ) - this->get ( x + 1, y  , z ) +
                                     this->get ( x  , y + 1, z + 1 ) - this->get ( x  , y + 1, z ) +
                                     this->get ( x + 1, y + 1, z + 1 ) - this->get ( x + 1, y + 1, z ) ) / 4. );
  }

  RealType dx ( const aol::Vec3<RealType> &coord ) const {
    return dx ( coord[0], coord[1], coord[2] );
  }

  RealType dy ( const aol::Vec3<RealType> &coord ) const {
    return dy ( coord[0], coord[1], coord[2] );
  }

  RealType dz ( const aol::Vec3<RealType> &coord ) const {
    return dz ( coord[0], coord[1], coord[2] );
  }

  /**
   * @param   coordinates
   * @returns x-derivative at arbitrary coordinates
   */
  RealType dx ( const RealType &X, const RealType &Y, const RealType &Z ) const {
    int x = static_cast<int> ( X ), y = static_cast<int> ( Y ), z = static_cast<int> ( Z );

    RealType a2 = Y - static_cast<RealType> ( y );
    RealType a3 = Z - static_cast<RealType> ( z );


    if ( x >= this->numX - 1 || y >= this->numY - 1 || z >= this->numZ - 1 ) {
      cerr << "ERROR in ScalarArray<QC_3D>::dx() x = " << x << " y = " << y << endl;
    }

    return static_cast<RealType> ( a2 * ( a3 * this->get ( x + 1, y + 1, z + 1 ) + ( 1.0 - a3 ) * this->get ( x + 1, y + 1, z ) ) +
                                   ( 1.0 - a2 ) * ( a3 * this->get ( x + 1, y, z + 1 ) + ( 1.0 - a3 ) * this->get ( x + 1, y, z ) ) -
                                   ( a2 * ( a3 * this->get ( x, y + 1, z + 1 ) + ( 1.0 - a3 ) * this->get ( x, y + 1, z ) ) +
                                     ( 1.0 - a2 ) * ( a3 * this->get ( x, y, z + 1 ) + ( 1.0 - a3 ) * this->get ( x, y, z ) ) ) );
  }

  /**
   * @param   coordinates
   * @returns y-derivative at arbitrary coordinates
   */
  RealType dy ( const RealType &X, const RealType &Y, const RealType &Z ) const {
    int x = static_cast<int> ( X ), y = static_cast<int> ( Y ), z = static_cast<int> ( Z );

    RealType a1 = X - static_cast<RealType> ( x );
    RealType a3 = Z - static_cast<RealType> ( z );

    if ( x >= this->numX - 1 || y >= this->numY - 1 || z >= this->numZ - 1 ) {
      cerr << "ERROR in ScalarArray<QC_3D>::dx() x = " << x << " y = " << y << endl;
    }

    return static_cast<RealType> ( a1 * ( a3 * this->get ( x + 1, y + 1, z + 1 ) + ( 1.0 - a3 ) * this->get ( x + 1, y + 1, z ) ) +
                                   ( 1.0 - a1 ) * ( a3 * this->get ( x, y + 1, z + 1 ) + ( 1.0 - a3 ) * this->get ( x, y + 1, z ) ) -
                                   ( a1 * ( a3 * this->get ( x + 1, y, z + 1 ) + ( 1.0 - a3 ) * this->get ( x + 1, y, z ) ) +
                                     ( 1.0 - a1 ) * ( a3 * this->get ( x, y, z + 1 ) + ( 1.0 - a3 ) * this->get ( x, y, z ) ) ) );
  }

  /**
   * @param   coordinates
   * @returns z-derivative at arbitrary coordinates
   */
  RealType dz ( const RealType &X, const RealType &Y, const RealType &Z ) const {

    int x = static_cast<int> ( X ), y = static_cast<int> ( Y ), z = static_cast<int> ( Z );

    RealType a1 = X - static_cast<RealType> ( x );
    RealType a2 = Y - static_cast<RealType> ( y );


    if ( x >= this->numX - 1 || y >= this->numY - 1 || z >= this->numZ - 1 ) {
      cerr << "ERROR in ScalarArray<QC_3D>::dx() x = " << x << " y = " << y << endl;
    }

    return static_cast<RealType> ( a1 * ( a2 * this->get ( x + 1, y + 1, z + 1 ) + ( 1.0 - a2 ) * this->get ( x + 1, y, z + 1 ) ) +
                                   ( 1.0 - a1 ) * ( a2 * this->get ( x, y + 1, z + 1 ) + ( 1.0 - a2 ) * this->get ( x, y, z + 1 ) ) -
                                   ( a1 * ( a2 * this->get ( x + 1, y + 1, z ) + ( 1.0 - a2 ) * this->get ( x + 1, y, z ) ) +
                                     ( 1.0 - a1 ) * ( a2 * this->get ( x, y + 1, z ) + ( 1.0 - a2 ) * this->get ( x, y, z ) ) ) );
  }

  // \brief Makes this array the scaled version of image
  // \author toelkes
  void resampleFrom ( const ScalarArray<DataType, qc::QC_3D> &image ) {
    const int sRanX[2] = { 0, image.getNumX() };
    const int sRanY[2] = { 0, image.getNumY() };
    const int sRanZ[2] = { 0, image.getNumZ() };
    const int tRanX[2] = { 0, this->getNumX() };
    const int tRanY[2] = { 0, this->getNumY() };
    const int tRanZ[2] = { 0, this->getNumZ() };

    ResampleCopier copier ( image, *this );
    fillValuesFrom ( image, sRanX, sRanY, sRanZ, tRanX, tRanY, tRanZ, copier, true );
  }

  //! \brief Universal function for copying, padding, pasting...
  //! \author toelkes
  template <typename CopyFunctor>
  void fillValuesFrom ( const ScalarArray<DataType, qc::QC_3D> &source, const int sRanX[2], const int sRanY[2],  const int sRanZ[2],
                        const int tRanX[2], const int tRanY[2], const int tRanZ[2], CopyFunctor getValue, bool differentRangeSizesOk = false ) {
    // Check if the ranges are valid
    if ( sRanX[0] > sRanX[1] || sRanY[0] > sRanY[1] || sRanZ[0] > sRanZ[1]
         || tRanX[0] > tRanX[1] || tRanY[0] > tRanY[1] || tRanZ[0] > tRanZ[1] ) {
      throw aol::Exception ( "ScalarArray<QC_3D>::fillValuesFrom: One or more of the ranges are incorrect: Ran[1] must not be bigger than Ran[0]!\n", __FILE__, __LINE__ );
    }

    // Check if ranges are >= 0
    if ( sRanX[0] < 0 || sRanY[0] < 0 || sRanZ[0] < 0 || tRanX[0] < 0 || tRanY[0] < 0 || tRanZ[0] < 0 ) {
      throw aol::Exception ( "ScalarArray<QC_3D>::fillValuesFrom: One or more of the ranges are incorrect: They must not have elements below 0!\n", __FILE__, __LINE__ );
    }

    // Check if source range is valid. This has to be done by the copier because different copiers have different prerequisites.
    if ( !getValue.sourceRangeOk ( source, sRanX, sRanY, sRanZ, tRanX, tRanY, tRanZ ) ) {
      throw aol::Exception ( "ScalarArray<QC_3D>::fillValuesFrom(): Source range invalid!\n", __FILE__, __LINE__ );
    }

    // Check if target range is valied, i.e. not outside of this.
    if ( tRanX[1] > this->getNumX() ) {
      throw aol::Exception ( "ScalarArray<QC_3D>::fillValuesFrom(): X range does not fit into the current array!\n", __FILE__, __LINE__ );
    }
    if ( tRanY[1] > this->getNumY() ) {
      throw aol::Exception ( "ScalarArray<QC_3D>::fillValuesFrom(): Y range does not fit into the current array!\n", __FILE__, __LINE__ );
    }
    if ( tRanZ[1] > this->getNumZ() ) {
      throw aol::Exception ( "ScalarArray<QC_3D>::fillValuesFrom(): Z range does not fit into the current array!\n", __FILE__, __LINE__ );
    }

    // Optionally check if the target range is as big as the source range
    if ( ( ( tRanX[1] - tRanX[0] ) != ( sRanX[1] - sRanX[0] )
           || ( tRanY[1] - tRanY[0] ) != ( sRanY[1] - sRanY[0] )
           || ( tRanZ[1] - tRanZ[0] ) != ( sRanZ[1] - sRanZ[0] ) ) && !differentRangeSizesOk ) {
      throw aol::Exception ( "ScalarArray<QC_3D>::fillValuesFrom(): Source range and target range have different sizes!\n", __FILE__, __LINE__ );
    }


    // copy
    for ( int x = 0; x < ( tRanX[1] - tRanX[0] ); ++x )
      for ( int y = 0; y < ( tRanY[1] - tRanY[0] ); ++y )
        for ( int z = 0; z < ( tRanZ[1] - tRanZ[0] ); ++z ) {
          // this->set ( tRanX[0] + x, tRanY[0] + y, tRanZ[0] + z, source.get ( sRanX[0] + x, sRanY[0] + y, sRanZ[0] + z ) );
          this->set ( tRanX[0] + x, tRanY[0] + y, tRanZ[0] + z, getValue ( sRanX[0] + x, sRanY[0] + y, sRanZ[0] + z, source ) );
        }
  }

  //! \brief Pastes image into current array. If the image is smaller in any direction, it is centered and padded with fillValue. If it is bigger, only its central values are copied.
  //! \author toelkes
  void padFrom ( const ScalarArray<DataType, qc::QC_3D> &image, const DataType fillValue = 0 ) {
    // integer division and conversion on purpose!
    qc::CoordType offsetb ( ( this->getSize() - image.getSize() ) / 2 );
    qc::CoordType offsett ( ( this->getSize() - image.getSize() ) / 2 );

    // If this->getNumX() - image.getNumX() (or getNumY or getNumZ respectively) is odd, ( this->getNumX() - image.getNumX() ) / 2 will be rounded down.
    // Therefore tRanX or sRanX (or the Y or Z ranges) will be 1 too large.
    // To fix this the remainder ( this->getNumX() - image.getNumX() ) % 2 is added to/subtracted from padX at one side of the target or source range
    // resulting in a shift of the image (padY analog).
    // Dimension is 3.
    for ( int i = 0; i < 3; ++i ) {
      if ( offsetb[i] < 0 )
	offsetb[i] += ( this->getSize()[i] - image.getSize()[i] ) % 2;
      else
	offsett[i] += ( this->getSize()[i] - image.getSize()[i] ) % 2;
    }

    this->setAll ( fillValue );

    const int tRanX[2] = { aol::Max ( aol::NumberTrait<short>::zero, offsetb[0] ), this->getNumX() - aol::Max ( aol::NumberTrait<short>::zero, offsett[0] ) };
    const int tRanY[2] = { aol::Max ( aol::NumberTrait<short>::zero, offsetb[1] ), this->getNumY() - aol::Max ( aol::NumberTrait<short>::zero, offsett[1] ) };
    const int tRanZ[2] = { aol::Max ( aol::NumberTrait<short>::zero, offsetb[2] ), this->getNumZ() - aol::Max ( aol::NumberTrait<short>::zero, offsett[2] ) };
    const int sRanX[2] = { -aol::Min ( aol::NumberTrait<short>::zero, offsetb[0] ), image.getNumX() + aol::Min ( aol::NumberTrait<short>::zero, offsett[0] ) };
    const int sRanY[2] = { -aol::Min ( aol::NumberTrait<short>::zero, offsetb[1] ), image.getNumY() + aol::Min ( aol::NumberTrait<short>::zero, offsett[1] ) };
    const int sRanZ[2] = { -aol::Min ( aol::NumberTrait<short>::zero, offsetb[2] ), image.getNumZ() + aol::Min ( aol::NumberTrait<short>::zero, offsett[2] ) };

    StandardCopier copier;
    fillValuesFrom ( image, sRanX, sRanY, sRanZ, tRanX, tRanY, tRanZ, copier );
  }

  void flipFrom ( const ScalarArray<DataType, qc::QC_3D> &image, const qc::Comp component ) {
    // make sure there is enough memory
    this->reallocate ( image );

    const int sRanX[2] = { 0, image.getNumX() };
    const int sRanY[2] = { 0, image.getNumY() };
    const int sRanZ[2] = { 0, image.getNumZ() };
    const int tRanX[2] = { 0, this->getNumX() };
    const int tRanY[2] = { 0, this->getNumY() };
    const int tRanZ[2] = { 0, this->getNumZ() };

    flipCopier copier ( this->getNumX(), this->getNumY(), this->getNumZ(), component );
    fillValuesFrom ( image, sRanX, sRanY, sRanZ, tRanX, tRanY, tRanZ, copier );
  }

  void copyBlockTo ( const int startX, const int startY, const int startZ, qc::ScalarArray<DataType, qc::QC_3D> &Block ) const {
    const int sRanX[2] = { startX, startX + Block.getNumX() };
    const int sRanY[2] = { startY, startY + Block.getNumY() };
    const int sRanZ[2] = { startZ, startZ + Block.getNumZ() };
    const int tRanX[2] = { 0, Block.getNumX() };
    const int tRanY[2] = { 0, Block.getNumY() };
    const int tRanZ[2] = { 0, Block.getNumZ() };

    StandardCopier copier;
    Block.fillValuesFrom ( *this, sRanX, sRanY, sRanZ, tRanX, tRanY, tRanZ, copier );
  }

  void copyBlockTo ( const aol::Vec<3, int> Start, qc::ScalarArray<DataType, qc::QC_3D> &Block ) const {
    copyBlockTo ( Start[0], Start[1], Start[2], Block );
  }
  
  void crop ( const aol::Vec<3, int> &CropStart, const aol::Vec<3, int> &CropSize )  {
    qc::ScalarArray<DataType, qc::QC_3D> copy ( *this, aol::DEEP_COPY );
    this->reallocate ( CropSize[0], CropSize[1], CropSize[2] );
    copy.copyBlockTo ( CropStart, *this );
  }

  //! Returns the value of the array at the given coordinates after
  //! median filtering.
  DataType getMedianFilterValue ( int X, int Y, int Z ) const {
#ifndef INTEL
    vector<DataType> median_sort_vec;

    int XMin, XMax, YMin, YMax, ZMin, ZMax;
    int Off = ( this->MED_FILTER_WIDTH - 1 ) >> 1;

    XMin = aol::Max ( 0, X - Off );
    YMin = aol::Max ( 0, Y - Off );
    ZMin = aol::Max ( 0, Z - Off );
    XMax = aol::Min ( this->getNumX() - 1, X + Off );
    YMax = aol::Min ( this->getNumY() - 1, Y + Off );
    ZMax = aol::Min ( this->getNumZ() - 1, Z + Off );

    for ( int x = XMin; x <= XMax; ++x ) {
      for ( int y = YMin; y <= YMax; ++y ) {
        for ( int z = ZMin; z <= ZMax; ++z ) {
          median_sort_vec.push_back ( this->get ( x, y, z ) );
        }
      }
    }

    sort ( median_sort_vec.begin(), median_sort_vec.end() );
    return median_sort_vec[ ( median_sort_vec.size() >> 1 ) ];
#else
    throw aol::Exception ( "qc::ScalarArray<QC_3D>::getMedianFilterValue() does not do the right thing if INTEL is defined", __FILE__, __LINE__ );
    return this->get ( X, Y, Z );
#endif
  }

  //! Returns the value at the given coordinates after convolution with the
  //! given filter kernel.
  RealType getConvolveValue ( int X, int Y, int Z, const Kernel3d<RealType> &Kernel ) const;

  /**
   * \brief Computes \f$\frac{\int_\Omega f(x)w(x)K(x-\bar x)dx}{\int_\Omega w(x)K(x-\bar x)dx}\f$
   * for kernel \f$K\f$, weight \f$w\f$, at point \f$\bar x\f$.
   */
  RealType getWeightedConvolveValue ( int X, int Y, int Z, const Kernel3d<RealType> &Kernel, const ScalarArray<RealType, qc::QC_3D> &Weight ) const;

  void gradientConv ( const aol::Vec3<short> &Coord, aol::Vec3<RealType> &Grad ) const {
    gradientConv ( Coord.x(), Coord.y(), Coord.z(), Grad );
  }

  void gradientConv ( int X, int Y, int Z, aol::Vec3<RealType> &Grad ) const {
    Grad[ 0 ] = getConvolveValue ( X, Y, Z, *diffKernelX );
    Grad[ 1 ] = getConvolveValue ( X, Y, Z, *diffKernelY );
    Grad[ 2 ] = getConvolveValue ( X, Y, Z, *diffKernelZ );
  }

  void setMax ( int X, int Y, int Z, DataType value ) {
    int in = Array<DataType>::index ( X, Y, Z );
    if ( this->_pData[ in ] < value ) this->_pData[ in ] = value;
  }

  void setMin ( int X, int Y, int Z, DataType value ) {
    int in = Array<DataType>::index ( X, Y, Z );
    if ( this->_pData[ in ] > value ) this->_pData[ in ] = value;
  }

  /** Load array in 2d pgm style format from the given stream. Note that if you would like
     *  to load compressed files directly you must use the method load(char *fileName),
     *  since only this method can open a pipe and determine the type of compression
     *  @throws aol::TypeException if the magic number of the file is unknown, or the dimensions of the array contained in that file do not match this one
     *  @throws aol::IOException if an io error occured
     *  \author Preusser, Droske
     */
  void load ( istream &in );

  void loadRaw ( istream &in, const int Type, const int InWidth, const int InHeight, const int InDepth );

  /** Load array in the standard 2d pgm style format from the file named fileName. If the file contains
   *  compressed data in the bz2, gz or Z format, it will automatically be decompressed
   *  by a pipe stream that is passed to load(istream &in).
   *  @throws aol::FileException if the file could not be opened
   *  \author Preusser, Droske
   */
  void load ( const char *fileName );


  /** Load array in 2d pgm style format from the files indicated by fileNameMask
   *  by calling the ScalarArray<QC_2D>::load function for each slice. The slices will be
   *  centered in the 3d array.
   *  @param fileNameMask The mask of the filenames of the slices in printf syntax, e.g., "data_%05d.pgm"
   *  @param dir the direction of how the slices should be stored in the array. @see putSlice
   *  @param begin the number of the first slice
   *  @param end the number of the last slice
   *  @attention end is the last existing slice
   *  @throws aol::FileException if one of the files could not be opened
   *  \author Droske
   */
  void loadSlices ( const char *fileNameMask, qc::Comp dir, int begin, int end );
  
  /** Load array in 2d pgm style format from the files indicated by fileNameMask
   *  by calling the ScalarArray<QC_2D>::load function for each slice. The slices will be
   *  centered in the 3d array.
   *  @param fileNameList A list of paths to the 2D image slices
   *  @param dir the direction of how the slices should be stored in the array. @see putSlice
   *  @throws aol::FileException if one of the files could not be opened
   *  \author Mevenkamp
   */
  void loadSlices ( const std::vector<std::string> &fileNameList, qc::Comp dir, bool resize = false );

  //! Load array in the MRC (Medical Research Council) file format.
  //! \author Berkels
  void loadMRC ( const char *FileName ) {
    std::ifstream in ( FileName );
    const int numX = aol::readBinaryData<int32_t, int> ( in );
    const int numY = aol::readBinaryData<int32_t, int> ( in );
    const int numZ = aol::readBinaryData<int32_t, int> ( in );
    const int type = aol::readBinaryData<int32_t, int> ( in );
    // The size of the MRC header is 1024 bytes. Since the values above are the only ones
    // we need to import the raw data from the MRC file, we skip directly to the data now.
    in.seekg ( 1024 );
    switch ( type ) {
      case 1:
        cerr << "Warning: Reading of int16 MRC files untested! If you see this warning, test if the data was loaded properly and remove the warning, if it was loaded properly.\n";
        loadRaw ( in, qc::PGM_SHORT_BINARY, numX, numY, numZ );
        break;
      case 2:
        loadRaw ( in, qc::PGM_FLOAT_BINARY, numX, numY, numZ );
        break;
      default:
        throw aol::TypeException ( "qc::ScalarArray<DataType, qc::QC_3D>::loadMRC: illegal type", __FILE__, __LINE__ );
    }
  }
  
  //! Load array in the Digital Micrograph 3 file format.
  //! \author Berkels
  void loadDM3 ( const char *FileName );
  
  //! Load array in the NetCDF file format
  //! Uses depth first search to find a data set within the container
  //! \author Mevenkamp
  void loadNetCDF ( const char *fileName ) {
    std::vector<std::string> groupNames;
    std::string dataName;
    netCDFgetDepthFirstDataAndGroups ( fileName, groupNames, dataName, 3 );
    loadNetCDF ( fileName, groupNames, dataName );
  }
  
  //! Load array in the NetCDF file format
  //! Name of groups and contained dataset is manually specified
  //! \author Mevenkamp
  void loadNetCDF ( const char *fileName, const std::vector<std::string> &groupNames, const std::string &dataName );
  
  //! Save array in NetCDF file format (data will be stored in a variable "data" inside the root group)
  //! \author Mevenkamp
  void saveNetCDF ( const char *fileName );
  
  //! Save array in NetCDF file format while trying to copy all meta data from the specified source file
  //! In case the format of the specified source file is invalid, saveNetCDF ( fileName ) will be called
  //! Else, the data from this ScalarArray will be written to the first suitable variable found by depth first search in the source file
  //! \author Mevenkamp
  void saveNetCDF ( const char *fileName, const char *sourceFileName,
                    const char *comment = "Data was modified by the QuocMesh library" ) {
    std::vector<std::string> groupNames;
    std::string dataName;
    netCDFgetDepthFirstDataAndGroups ( sourceFileName, groupNames, dataName, 3 );
    saveNetCDF ( fileName, sourceFileName, groupNames, dataName, comment );
  }
  
  //! Save array in NetCDF file format while trying to copy all meta data from the specified source file
  //! In case the format of the specified source file is invalid, saveNetCDF ( fileName ) will be called
  //! Else, the data will be written to the specified variable within the specified group
  //! \author Mevenkamp
  void saveNetCDF ( const char *fileName, const char *sourceFileName,
                    const std::vector<std::string> &groupNames, const std::string &dataName,
                    const char *comment = "Data was modified by the QuocMesh library" );
  
  //! Load array in the HDF5 file format
  //! Uses depth first search to find a data set within the container
  //! \author Mevenkamp
  void loadHDF5 ( const char *fileName ) {
    std::string dataName;
    hdf5GetDepthFirstData ( fileName, dataName );
    loadHDF5 ( fileName, dataName );
  }
  
  //! Load array in the HDF5 file format
  //! Name of groups and contained dataset is manually specified
  //! \author Mevenkamp
  void loadHDF5 ( const char *fileName, const std::string &dataName );
  
  //! Save array in HDF5 file format (data will be stored in a variable "data" inside the root group)
  //! \author Mevenkamp
  void saveHDF5 ( const char *fileName );
  
  //! Save array in NetCDF file format while trying to copy all meta data from the specified source file
  //! In case the format of the specified source file is invalid, saveHDF5 ( fileName ) will be called
  //! Else, the data from this ScalarArray will be written to the first suitable variable found by depth first search in the source file
  //! \author Mevenkamp
  void saveHDF5 ( const char *fileName, const char *sourceFileName,
                  const char *comment = "Data was modified by the QuocMesh library" ) {
    std::string dataName;
    hdf5GetDepthFirstData ( sourceFileName, dataName );
    saveHDF5 ( fileName, sourceFileName, dataName, comment );
  }
  
  //! Save array in HDF5 file format while trying to copy all meta data from the specified source file
  //! In case the format of the specified source file is invalid, saveHDF5 ( fileName ) will be called
  //! Else, the data will be written to the specified variable within the specified group
  //! \author Mevenkamp
  void saveHDF5 ( const char *fileName, const char *sourceFileName,
                  const std::string &dataName,
                  const char *comment = "Data was modified by the QuocMesh library" );

  /** Save array in 3d pgm style format of type type.
   *  See Array2d::save(ostream &out, int type, const char *comment) for details on the different file types
   *  @param out The stream to write the data to
   *  @param type The file type (see qc::SaveType)
   *  @param comment The comment to be written into the header of the file. If comment is NULL the date and the data type will be written as comment
   *  @throws aol::IOException if an io error occured
   *  \author Preusser, Droske
   */
  void save ( ostream &out, SaveType type, const char *comment = NULL ) const ;

  /** Save array in 3d pgm style format of type type.
   *  See qc::SaveType for details on the different file types
   *  This method is capable of directly compressing the data using bzip2, gzip or compress.
   *  The compressed files can be read again directly through a call of load(char *fileName).
   *  @param fileName The name of the file to be written. If compressed data shall be saved the suffix should correspond to the compression type (i.e. '.bz2', '.gz' or '.Z'). Other wise a warning will be displayed
   *  @param comment The comment to be written into the header of the file. If comment is NULL the date and the data type will be written as comment
   *  @throws aol::TypeException if the given type is incompatible
   *  @throws aol::IOException if an io error occured
   *  @throws aol::FileException if the file could not be opened
   *  @throws aol::PipeException if the compression pipe could not be opened
   *  \author Preusser, Droske
   */
  void save ( const char *fileName, SaveType type, const char *comment = NULL ) const ;

  //! lossless save in binary quocmesh format
  void saveToFile ( const char *filename ) const;

  //! load from file in quocmesh format, no conversion of data types
  void loadFromFile ( const char *filename );

  /** Save array as slices in 2d pgm style format of type type.
   *  See qc::SaveType for details on file types
   *  The number of slice is printed into the fileNameMask
   *  @param direction gives the slicing direction
   *  \author Preusser, Droske
   */
  void saveSlices ( const char *fileNameMask,
                    Comp direction,
                    SaveType type, const char *comment = NULL,
                    aol::OverflowHandlingType oh = aol::CLIP,
                    DataType overflowMin = aol::ZOTrait<RealType>::zero,
                    DataType overflowMax = aol::OverflowTrait<DataType>::max ) const;

  void saveASCII ( const char *fileName, const aol::Format &format = aol::scientificFormat ) const {
    std::ofstream file ( fileName );
    for ( int k = 0; k < this->getNumZ(); ++k ) {
      for ( int j = 0; j < this->getNumY(); ++j ) {
        for ( int i = 0; i < this->getNumX(); ++i )
          file << format ( this->get ( i, j, k ) ) << "\t";
        file << std::endl;
      }
      file << "====================================================" << std::endl;
    }
    file.close();
  }


  /**
   * Writes the array in the meta image data format (can be read by ParaView for example).
   * This format consists of two files, a header file (BaseFileName.mhd) and a raw data
   * file (BaseFileName.raw).
   *
   * \author Berkels
   */
  void saveMetaImageFile ( const char *BaseFileName, SaveType Type ) const ;

  /** Shifts the contents of Image by the given offset and stores the result in this, i.e.:
   *  this( (x + XOffset) % numX , (y + YOffset) % numY, (z + ZOffset) % numZ ) = Image( x, y, z )
   *  \author Berkels
   */
  void shiftByOffsetFrom ( const int XOffset,
                           const int YOffset,
                           const int ZOffset,
                           const qc::ScalarArray<DataType, qc::QC_3D> &Image ) {
    if ( ( this->numX != Image.getNumX() ) || ( this->numY != Image.getNumY() ) || ( this->numZ != Image.getNumZ() ) ) {
      throw aol::Exception ( "Dimiension don't match", __FILE__, __LINE__ );
    }
    for ( int z = 0; z < this->numZ; ++z ) {
      for ( int y = 0; y < this->numY; ++y ) {
        for ( int x = 0; x < this->numX; ++x ) {
          this->set ( ( x + XOffset ) % this->numX , ( y + YOffset ) % this->numY, ( z + ZOffset ) % this->numZ, Image.get ( x, y, z ) );
        }
      }
    }
  }

  /** Shifts the contents of Image by the given offset
   *  \author wirth
   */
  void shiftByOffsetFrom ( const qc::CoordType Offset,
                           const qc::ScalarArray<DataType, qc::QC_3D> &Image ) {
    shiftByOffsetFrom ( Offset[0], Offset[1], Offset[2], Image );
  }

  //! needed for element saturation
  DataType getElementSaturationMax ( int X, int Y, int Z ) const;

  //! needed for element saturation
  DataType getElementSaturationMin ( int X, int Y, int Z ) const;


  template <class T>
  void convertFrom ( const T &OtherImage ) {
    const int size = this->getNumX() * this->getNumY() * this->getNumZ();
    for ( int i = 0; i < size; ++i ) this->set ( i, static_cast<DataType> ( OtherImage.get ( i ) ) );
  }

  //! return number of points in all space dimensions, if equal, throw exception otherwise
  int getNumXYZ ( ) const {
    if ( ( this->getNumX() != this->getNumY() ) || ( this->getNumY() != this->getNumZ() ) )
      throw aol::Exception ( "qc::ScalarArray<DataType,QC_3D>::getNumXYZ() may not be called for non-cubic arrays", __FILE__, __LINE__ );

    return ( this->getNumX() );
  }

  // We do not want to be able to use all set and add methods from the parent class, so we explicitly implement those that should be usable.
  // Unfortunately, this cannot be done via "using".
  void set ( const CoordType &Coord, DataType V ) {
    Array<DataType>::set ( Coord, V );
  }

  // sets a value of the array
  void set ( int X, int Y, int Z, DataType V ) {
    Array<DataType>::set ( X, Y, Z, V );
  }

  using aol::Vector<DataType>::set;

  DataType get ( const CoordType &Coord ) const {
    return Array<DataType>::get ( Coord );
  }

  // sets a value of the array
  DataType get ( int X, int Y, int Z ) const {
    return Array<DataType>::get ( X, Y, Z );
  }

  using aol::Vector<DataType>::get;

  DataType &getReference ( int X, int Y, int Z ) {
    return Array<DataType>::getReference ( X, Y, Z );
  }

  DataType &getReference ( const CoordType &Coords ) {
    return Array<DataType>::getReference ( Coords );
  }

  void add ( int X, int Y, int Z, DataType V ) {
    Array<DataType>::add ( X, Y, Z, V );
  }

  void add ( const CoordType &Pt, DataType V ) {
    Array<DataType>::add ( Pt, V );
  }

  using aol::Vector<DataType>::add;

  void resize ( const int newNumX, const int newNumY, const int newNumZ );

  /** save as 8-bit integer df3 file that can be read by povray, only tested for floating point data types so far */
  void saveAsPovrayVolume8Bit ( const char* const Fname, const DataType MinValue, const DataType MaxValue ) const {
    this->template saveAsPovrayVolume<unsigned char> ( Fname, MinValue, MaxValue );
  }

  /** save as 16-bit integer df3 file that can be read by povray, only tested for floating point data types so far */
  void saveAsPovrayVolume16Bit ( const char* const Fname, const DataType MinValue, const DataType MaxValue ) const {
    this->template saveAsPovrayVolume<unsigned short> ( Fname, MinValue, MaxValue );
  }

  /** save as 32-bit integer df3 file that can be read by povray, only tested for floating point data types so far */
  void saveAsPovrayVolume32Bit ( const char* const Fname, const DataType MinValue, const DataType MaxValue ) const {
    this->template saveAsPovrayVolume<unsigned int> ( Fname, MinValue, MaxValue );
  }

  void writePovrayVolumeHeader ( std::ostream &Out ) const;

  // only call with unsigned char, short, int (method cannot be made protected as it needs to be called for other template type
  template<typename IntType>
  void saveAsPovrayVolume ( const char* const Fname, const DataType MinValue, const DataType MaxValue ) const {
    std::ofstream out ( Fname );
    this->writePovrayVolumeHeader ( out );

    const IntType maxc = std::numeric_limits<IntType>::max();
    const DataType maxr = static_cast<DataType> ( maxc );

    qc::ScalarArray< RealType, qc::QC_3D > scaledData ( this->getSize() );
    scaledData.convertFrom ( *this );
    scaledData.addToAll ( -MinValue );
    scaledData *= maxr / ( MaxValue - MinValue );

    qc::ScalarArray< IntType, qc::QC_3D > intData ( this->getNumX(), this->getNumY(), this->getNumZ() );
    for ( int i = 0; i < scaledData.size(); ++i ) {
      const DataType value = scaledData[i];
      if ( value < 0.0 ) {
        intData[i] = 0;
      } else if ( value > maxr ) {
        intData[i] = maxc;
      } else {
        intData[i] = static_cast<IntType> ( floor ( scaledData[i] ) );
      }
    }

    intData.swapByteOrder(); // does no harm for one-byte case
    const char* buffer = reinterpret_cast<char*> ( intData.getData() );
    out.write ( buffer, intData.size() * sizeof ( IntType ) );
  }

protected:
  GaussDiffKernel3d<RealType> *diffKernelX, *diffKernelY, *diffKernelZ;

private:
  void hdf5GetDepthFirstData ( const char *fileName, std::string &dataName );

  inline RealType integrateOnSquare ( RealType OffX, RealType OffY, RealType OffZ,
                                      int BasisFuncNum, const Hexahedron &Hex ) const;

  void resize ( const int ) {
    throw aol::Exception ( "qc::Array::resize( int ) preserving contents does not make sense.", __FILE__, __LINE__ );
  }

  void reallocate ( const int ) {
    throw aol::Exception ( "qc::Array::reallocate( int ) preserving contents does not make sense.", __FILE__, __LINE__ );
  }

};


template <class _DataType, int dim>
class ScalarArrayTrait {};

template <class _DataType>
class ScalarArrayTrait<_DataType, 1> {
public:
  typedef qc::ScalarArray<_DataType, qc::QC_1D> ArrayType;
};

template <class _DataType>
class ScalarArrayTrait<_DataType, 2> {
public:
  typedef qc::ScalarArray<_DataType, qc::QC_2D> ArrayType;
};

template <class _DataType>
class ScalarArrayTrait<_DataType, 3> {
public:
  typedef qc::ScalarArray<_DataType, qc::QC_3D> ArrayType;
};

  
  
template <typename DataType>
class NetCDFTrait {
public:
  static int getVar ( int /*ncid*/, int /*varid*/, DataType * /*p*/ ) {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  static int defVar ( int /*ncid*/, const char* /*varname*/, int /*ndims*/, int * /*dimids*/, int * /*varid*/ ) {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  static int putVar ( int /*ncid*/, int /*varid*/, DataType * /*p*/ ) {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};

#ifdef USE_LIB_HDF5
static inline H5::DataSpace hdf5GetDataSet ( const aol::Vec3<int> &Size ) {
  hsize_t dims[3];
  dims[0] = Size[2];
  dims[1] = Size[1];
  dims[2] = Size[0];
  
  return H5::DataSpace ( 3, dims );
}

static inline H5::DSetCreatPropList hdf5GetCreatePList ( const aol::Vec3<int> &Size ) {
  aol::Vec3<int> size;
  size[0] = Size[2];
  size[1] = Size[1];
  size[2] = Size[0];
  
  H5::DSetCreatPropList plist;
  hsize_t chunk_dims[3];
  for ( int i=0; i<3 ; ++i )
    chunk_dims[i] = aol::Min<int> ( size[i], aol::Max<int> ( 10, static_cast<int> ( log2 ( size[i] ) ) ) );
  plist.setChunk ( 3, chunk_dims );
  plist.setDeflate ( 5 );
  plist.setShuffle ( );
  
  return plist;
}

static inline herr_t hdf5GetDataSets ( hid_t /*object_id*/, const char *name, const H5O_info_t *info, void *op_data ) {
  if ( info->type == H5O_TYPE_DATASET )
    static_cast<std::vector<std::string>*> ( op_data )->push_back ( aol::strprintf ( "%s", name ) );
  
  return 0;
}
  
template <typename DataType>
class HDF5Trait {
public:
  static void read ( const H5::DataSet &/*DataSet*/, void */*buf*/ ) {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }

  static void createDataSet ( const H5::H5File &/*File*/, const char* /*name*/, const aol::Vec3<int> &/*Size*/ ) {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  static void write ( const H5::DataSet &/*DataSet*/, void */*buf*/ ) {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};
  
template <>
class HDF5Trait<float> {
public:
  static void read ( const H5::DataSet &DataSet, void *buf ) {
    DataSet.read ( buf, H5::PredType::NATIVE_FLOAT );
  }
  
  static void createDataSet ( const H5::H5File &File, const char* name, const aol::Vec3<int> &Size ) {
    File.createDataSet ( name, H5::PredType::NATIVE_FLOAT, hdf5GetDataSet ( Size ), hdf5GetCreatePList ( Size ) );
  }
  
  static void write ( const H5::DataSet &DataSet, void *buf ) {
    DataSet.write ( buf, H5::PredType::NATIVE_FLOAT );
  }
};
#endif
  
} // end namespace qc

#endif
