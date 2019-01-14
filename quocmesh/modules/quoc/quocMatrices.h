#ifndef __QUOCMATRICES_H
#define __QUOCMATRICES_H

#include <CSRMatrix.h>
#include <UGBMatrix.h>
#include <bandMatrix.h>
#include <sparseMatrices.h>
#include <fastUniformGridMatrix.h>

namespace qc {

/** Special sparse rows in case of uniform grids */
template <typename DataType>
class UniformGridSparseRow : public aol::Row<DataType> {
protected:
  const GridDefinition &_grid;
  DataType *_thisRow;
  int _lastIndex, _middleIndex;

public:

  explicit UniformGridSparseRow ( const GridDefinition &Grid ) : _grid ( Grid ), _thisRow ( NULL ) {
    // line does not need to know which line it is
    // now we can use the same instance for different lines
    if ( _grid.getDimOfWorld() == QC_2D ) {
      _lastIndex = 9;
      _middleIndex = 4;
    } else if ( _grid.getDimOfWorld() == QC_3D ) {
      _lastIndex = 27;
      _middleIndex = 13;
    } else {
      throw aol::Exception ( "UniformGridSparseRow not implemented yet for dimension different from  2, 3", __FILE__, __LINE__ );
    }
    _thisRow = new DataType[ _lastIndex ];
    // initialize this array!
    for ( int i = 0; i < _lastIndex; ++i )
      _thisRow[i] = 0;
  }

  ~UniformGridSparseRow() {
    if ( _thisRow ) delete[] _thisRow;
  }

  explicit UniformGridSparseRow ( const UniformGridSparseRow<DataType> &other );

  int numNonZeroes() const {
    throw aol::Exception ( "qc::UniformGridSparseRow<DataType>::numNonZeroes() does not return what the method name suggests", __FILE__, __LINE__ );
    return _lastIndex;
  }

  DataType get ( int I, int J ) const {
    if ( I == J )
      return _thisRow[_middleIndex];
    for ( int k = 0; k < _middleIndex; ++k )
      if ( J == ( I + _grid.getUniformGridIndexOffset ( k ) ) )
        return _thisRow[k];
    for ( int k = _middleIndex + 1; k < _lastIndex; ++k )
      if ( J == ( I + _grid.getUniformGridIndexOffset ( k ) ) )
        return _thisRow[k];
    // else:
    return aol::NumberTrait<DataType>::zero;
  }

  void set ( int I, int J, DataType Value ) {
#ifdef BOUNDS_CHECK
    bool okay = ( Value == static_cast<DataType> ( 0.0 ) );
#endif
    for ( int k = 0; k < _lastIndex; ++k )
      if ( J == I + _grid.getUniformGridIndexOffset ( k ) ) {
        _thisRow[k] = Value;
#ifdef BOUNDS_CHECK
        okay = true;
#endif
      }
#ifdef BOUNDS_CHECK
    if ( !okay ) {
      cerr << I << ", " << J << endl;
      throw aol::Exception ( "UniformGridSparseRow::set: This position may not be used.", __FILE__, __LINE__ );
    }
#endif
  }

  /** row-vector scalar mulitplication
   */
  DataType mult ( const aol::Vector<DataType> &src, const int rowIndex ) const {
    DataType dst = static_cast<DataType> ( 0.0 );
    for ( int k = 0; k < _lastIndex; ++k )
      dst += _thisRow[k] * src[ rowIndex + _grid.getUniformGridIndexOffset ( k ) ];
    return ( dst );
  }


  virtual DataType multMaskedFunctorTrue ( const aol::Vector<DataType> &Src, int Row, const aol::BitVector & Mask ) {
    return multMasked<aol::BitMaskFunctorTrue> ( Src, Row, Mask );
  }
  virtual DataType multMaskedFunctorFalse ( const aol::Vector<DataType> &Src, int Row, const aol::BitVector & Mask ) {
    return multMasked<aol::BitMaskFunctorFalse> ( Src, Row, Mask );
  }
  virtual DataType multMaskedFunctorIdentity ( const aol::Vector<DataType> &Src, int Row, const aol::BitVector & Mask ) {
    return multMasked<aol::BitMaskFunctorIdentity> ( Src, Row, Mask );
  }
  virtual DataType multMaskedFunctorNegate ( const aol::Vector<DataType> &Src, int Row, const aol::BitVector & Mask ) {
    return multMasked<aol::BitMaskFunctorNegate> ( Src, Row, Mask );
  }

  /** masked row-vector scalar mulitplication (if desired only include nodes that are (not) masked
    */
  template <typename BitMaskFunctorType>
  DataType multMasked ( const aol::Vector<DataType> &src, int rowIndex,
                        const aol::BitVector & Mask) {
    BitMaskFunctorType maskFunctor;
    DataType dst = static_cast<DataType> ( 0.0 );
    for ( int k = 0; k < _lastIndex; ++k )
      if (maskFunctor(Mask[rowIndex + _grid.getUniformGridIndexOffset ( k )]))
        dst += _thisRow[k] * src[ rowIndex + _grid.getUniformGridIndexOffset ( k ) ];
    return ( dst );
  }


  /** row is mapped to scalar * row
   */
  void scale ( DataType factor ) {
    for ( int k = 0; k < _lastIndex; ++k )
      _thisRow[k] *= factor;
  }

  /** compute row sum
   */
  DataType sum ( int ) {
    // I is ignored
    DataType result = 0;
    for ( int k = 0; k < _lastIndex; ++k ) {
      result += _thisRow[k];
    }
    return ( result );
  }

  void scale ( int, DataType factor ) {  // 1st parameter is ignored
    scale ( factor );
  }

  void add ( int I, int J, DataType Value ) {
#ifdef BOUNDS_CHECK
    bool okay = ( Value == static_cast<DataType> ( 0.0 ) );
#endif
    for ( int k = 0; k < _lastIndex; ++k ) {
      if ( J == ( I + _grid.getUniformGridIndexOffset ( k ) ) ) {
        _thisRow[k] += Value;
#ifdef BOUNDS_CHECK
        okay = true;
#endif
      }
    }
#ifdef BOUNDS_CHECK
    if ( !okay ) {
      cerr << I << ", " << J << endl;
      throw aol::Exception ( "UniformGridSparseRow::add: This position may not be used.", __FILE__, __LINE__ );
    }
#endif
  }

  void setZero() {
    memset ( _thisRow, 0, sizeof ( DataType ) *_lastIndex );
  }


  bool checkForNANsAndINFs() const {
    for ( int i = 0; i < _lastIndex; ++i ) {
      if ( !aol::isFinite ( _thisRow[ i ] ) ) {
        return true;
      }
    }
    return false;
  }

  UniformGridSparseRow<DataType>& operator= ( const UniformGridSparseRow<DataType> &from );
    /*
      // This cannot work because there is no member row. Apparently, this method was not used so I won't think about how to fix it.
      if ( this->_lastIndex != from._lastIndex && this->row ) {
      delete this->row;
      this->_lastIndex = from._lastIndex;
      this->row = new DataType[this->_lastIndex];
      }
      for ( int i = 0; i < _lastIndex; ++i ) {
      this->row[i] = from.row[i];
      }
      this->_middleIndex = from._middleIndex;
    */

protected:

  inline int getIndex ( int offset ) {
    if ( _grid.getDimOfWorld() == QC_2D ) {
      if ( offset < -1 ) {
        return offset + _grid.getWidth() + 1;
      } else if ( offset > 1 ) {
        return offset - _grid.getWidth() + 7;
      }
      return offset + 4;
    } else if ( _grid.getDimOfWorld() == QC_3D ) {
      throw aol::UnimplementedCodeException ( "UniformGridSparseRow::getIndex not implemented for 3D", __FILE__, __LINE__ );
      // and it never was ...
    }
  }

  virtual void makeRowEntries ( vector<typename aol::Row<DataType>::RowEntry > &vec, const int RowNum ) const {
    vec.resize ( _lastIndex );
    for ( int i = 0; i < _lastIndex; ++i ) {
      vec[i].col   = RowNum + _grid.getUniformGridIndexOffset ( i ) ;
      vec[i].value = _thisRow[i];
#ifdef BOUNDS_CHECK
      if ( vec[i].col < 0 || vec[i].col > _grid.getNumberOfNodes() ) {
        cerr << vec[i].col << endl;
        throw aol::Exception ( "UniformGridSparseRow::makeRowEntries: Illegal entry in makeRowEntries", __FILE__, __LINE__ );
      }
#endif
    }
  }

  virtual void makeRowSortedEntries ( vector<typename aol::Row<DataType>::RowEntry > &vec, const int RowNum ) const {
    makeRowEntries ( vec, RowNum ); // are already sorted
  }

};

//! Special sparse matrix type for the sparsity structure resulting from quoc grids
/**
 * \ingroup Matrix
 */
template <typename DataType>
class UniformGridSparseMatrix : public aol::GenSparseMatrix<DataType> {
  const GridDefinition _grid; // apparently must have instance?

  UniformGridSparseMatrix ( int Rows, int Columns );
  // : aol::GenSparseMatrix<DataType> ( Rows, Columns ), _grid ( ) {}

public:

  explicit UniformGridSparseMatrix ( const GridDefinition &grid )
      : aol::GenSparseMatrix<DataType> ( grid.getNumberOfNodes(), grid.getNumberOfNodes() ), _grid ( grid ) {
    // initialize rows:
    init ( grid );
  }

  template <Dimension Dim>
  explicit UniformGridSparseMatrix ( const GridSize<Dim> &gridSize )
      : aol::GenSparseMatrix<DataType> ( gridSize.getNumberOfNodes(), gridSize.getNumberOfNodes() ), _grid ( gridSize ) {
    // check if given grid size is equally dimensioned.
    gridSize.quadraticOrDie ();
    // initialize rows:
    init ( _grid );
  }

  virtual UniformGridSparseRow<DataType>* newDefaultRow() const {
    return new UniformGridSparseRow<DataType> ( _grid );
  }

  void init ( const GridDefinition &grid ) {
    GridDefinition::OldFullBoundaryNodeIterator fbnit;
    GridDefinition::OldFullNodeIterator fnit;

    for ( fbnit = grid.begin(); fbnit != grid.end(); ++fbnit ) {
      int i = fbnit->x() +
              grid.getWidth() * fbnit->y() +
              grid.getWidth() * grid.getWidth() * fbnit->z();
      this->rows[ i ] = new aol::SparseRow<DataType>();
    }

    for ( fnit = grid.begin(); fnit != grid.end(); ++fnit ) {
      int i = fnit->x() +
              grid.getWidth() * fnit->y() +
              grid.getWidth() * grid.getWidth() * fnit->z();
      if ( ! ( this->rows[i] ) )
        this->rows[ i ] = new qc::UniformGridSparseRow<DataType> ( grid );
    }

  }

  virtual ~UniformGridSparseMatrix() {
    destroy();
  }

private:

  void destroy() {
    for ( int i = 0; i < static_cast<int> ( this->rows.size() ); ++i ) {
      delete this->rows[ i ];
    }
  }

public:

  void reallocate ( const int NRows, const int NCols ) {
    if ( NRows != NCols ) {throw aol::Exception ( "UniformGridSparseMatrix::reallocate: NRows != NCols", __FILE__, __LINE__ );
    }

    // ATTN: Small errors in pow may cause round-down to wrong number if you do not add an epsilon before casting to int!
    int width = static_cast<int> ( pow ( static_cast<double> ( NRows ), 1. / 3. ) + 0.5 );
    int depth = qc::logBaseTwo ( width );
    if ( NRows == static_cast<int> ( pow ( static_cast<double> ( ( 1 << depth ) + 1 ), 3. ) + 0.5 ) ) {
      GridDefinition def ( depth, QC_3D );
      // ATTN: Do not use reference to nameless temporary!
      reallocate ( def );
    } else {
      width = static_cast<int> ( pow ( static_cast<double> ( NRows ), 1. / 2. ) + 0.5 );
      depth = qc::logBaseTwo ( width );
      if ( NRows == static_cast<int> ( pow ( static_cast<double> ( ( 1 << depth ) + 1 ), 2. ) + 0.5 ) ) {
        GridDefinition def ( depth, QC_2D );
        reallocate ( def );
      } else {
        throw aol::Exception ( "NRows must be (2^n+1)^d for some n and some d=2,3.", __FILE__, __LINE__ );
      }
    }
  }

  //!resize Matrix. ATTENTION: this does not update the grid reference!
  virtual void reallocate ( const GridDefinition &grid ) {
    destroy();
    this->_numRows = this->_numCols = grid.getNumberOfNodes();
    init ( grid );
  }

  void setIdentity() {
    throw aol::Exception ( "UniformGridSparseMatrix<T>::setIdentity(): Standard setIdentity() incompatible with this type of matrix.", __FILE__, __LINE__ );
    // could be done via init and explicitely setting entries. but is it necessary at all?
  }
};


//! Basis class for band matrices for multilinear FE on quoc grids
/**
 * \ingroup Matrix
 */
template <class DataTye, qc::Dimension dim>
class MultilinFEBandMatrix {};

/** \brief A class for band matrices for use with 9-point stencils (e. g. in 2D bilinear finite element computations)
 *  \author Schwen
 *  \ingroup Matrix
 */
template<typename DataType>
class MultilinFEBandMatrix< DataType, qc::QC_2D > : public aol::GenBandMatrix<DataType> {
protected:
  int _wi, _he;

public:
  typedef BitArray<qc::QC_2D> MaskType;

  //! Constructor adapted to this matrix type,
  //@param _wi width of the grid, not! number of columns
  //@param _he height of the grid, not! number of rows
  MultilinFEBandMatrix ( const int wi, const int he ) : aol::GenBandMatrix<DataType> ( wi*he, wi*he ) {
    init ( wi, he );
  }

  explicit MultilinFEBandMatrix ( const GridStructure &Struc ) : aol::GenBandMatrix<DataType> ( Struc.getNumberOfNodes(), Struc.getNumberOfNodes() ) {
    init ( Struc.getNumX(), Struc.getNumY() );
  }

  explicit MultilinFEBandMatrix ( const GridSize<qc::QC_2D> &gridSize ) : aol::GenBandMatrix<DataType> ( gridSize.getNumberOfNodes(), gridSize.getNumberOfNodes() ) {
    init ( gridSize.getNumX(), gridSize.getNumY() );
  }

  //! optimized apply for a matrix of this structure
  virtual void apply ( const aol::Vector<DataType> &src, aol::Vector<DataType> &dst ) const {
    // might want to compare sizes
    dst.setZero();
    applyAdd ( src, dst );
  }

  //! optimized applyAdd for a matrix of this structure
  virtual void applyAdd ( const aol::Vector<DataType> &src, aol::Vector<DataType> &dst ) const {
    // might want to compare sizes
    {      // bottom left corner [0,0] -> 0
      const int k = 0;
      dst[k] +=
        src[k]              * this->_pData[this->map_index ( 4,k ) ] +
        src[k+1]            * this->_pData[this->map_index ( 5,k ) ] +
        src[k+_wi]          * this->_pData[this->map_index ( 7,k ) ] +
        src[k+_wi+1]        * this->_pData[this->map_index ( 8,k ) ];
    }

    {     // bottom right corner [0, wi-1] -> wi-1
      const int k = _wi - 1;
      dst[k] +=
        src[k-1]           * this->_pData[this->map_index ( 3,k ) ] +
        src[k]             * this->_pData[this->map_index ( 4,k ) ] +
        src[k+_wi-1]       * this->_pData[this->map_index ( 6,k ) ] +
        src[k+_wi]         * this->_pData[this->map_index ( 7,k ) ];
    }

    {      // top left corner[_he-1, 0] -> (_he-1) * wi
      const int k = ( _he - 1 ) * _wi;
      dst[k] +=
        src[k-_wi]         * this->_pData[this->map_index ( 1,k ) ] +
        src[k-_wi+1]       * this->_pData[this->map_index ( 2,k ) ] +
        src[k]             * this->_pData[this->map_index ( 4,k ) ] +
        src[k+1]           * this->_pData[this->map_index ( 5,k ) ];
    }

    {      // top right corner[_he-1, wi-1] -> _he * wi - 1
      const int k = _he * _wi - 1;
      dst[k] +=
        src[k-_wi-1]       * this->_pData[this->map_index ( 0,k ) ] +
        src[k-_wi]         * this->_pData[this->map_index ( 1,k ) ] +
        src[k-1]           * this->_pData[this->map_index ( 3,k ) ] +
        src[k]             * this->_pData[this->map_index ( 4,k ) ];
    }

    // left boundary
    for ( int i = 1; i < _he - 1; ++i ) {
      const int k = i * _wi;
      dst[k] +=
        src[k-_wi]         * this->_pData[this->map_index ( 1,k ) ] +
        src[k-_wi+1]       * this->_pData[this->map_index ( 2,k ) ] +
        src[k]             * this->_pData[this->map_index ( 4,k ) ] +
        src[k+1]           * this->_pData[this->map_index ( 5,k ) ] +
        src[k+_wi]         * this->_pData[this->map_index ( 7,k ) ] +
        src[k+_wi+1]       * this->_pData[this->map_index ( 8,k ) ];
    }

    // right boundary
    for ( int i = 1; i < _he - 1; ++i ) {
      const int k = i * _wi + _wi - 1;
      dst[k] +=
        src[k-_wi-1]       * this->_pData[this->map_index ( 0,k ) ] +
        src[k-_wi]         * this->_pData[this->map_index ( 1,k ) ] +
        src[k-1]           * this->_pData[this->map_index ( 3,k ) ] +
        src[k]             * this->_pData[this->map_index ( 4,k ) ] +
        src[k+_wi-1]       * this->_pData[this->map_index ( 6,k ) ] +
        src[k+_wi]         * this->_pData[this->map_index ( 7,k ) ];
    }

    // bottom boundary
    for ( int i = 1; i < _wi - 1; ++i ) {
      const int k = i;
      dst[k] +=
        src[k-1]           * this->_pData[this->map_index ( 3,k ) ] +
        src[k]             * this->_pData[this->map_index ( 4,k ) ] +
        src[k+1]           * this->_pData[this->map_index ( 5,k ) ] +
        src[k+_wi-1]       * this->_pData[this->map_index ( 6,k ) ] +
        src[k+_wi]         * this->_pData[this->map_index ( 7,k ) ] +
        src[k+_wi+1]       * this->_pData[this->map_index ( 8,k ) ];
    }

    // top boundary
    for ( int i = 1; i < _wi - 1; ++i ) {
      const int k = ( _he - 1 ) * _wi + i;
      dst[k] +=
        src[k-_wi-1]       * this->_pData[this->map_index ( 0,k ) ] +
        src[k-_wi]         * this->_pData[this->map_index ( 1,k ) ] +
        src[k-_wi+1]       * this->_pData[this->map_index ( 2,k ) ] +
        src[k-1]           * this->_pData[this->map_index ( 3,k ) ] +
        src[k]             * this->_pData[this->map_index ( 4,k ) ] +
        src[k+1]           * this->_pData[this->map_index ( 5,k ) ];
    }

    // interior points
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 1; i < _he - 1; ++i ) {
      for ( int j = 1; j < _wi - 1; ++j ) {
        int k = i * _wi + j;
        dst[k] +=
          src[k-_wi-1]       * this->_pData[this->map_index ( 0,k ) ] +
          src[k-_wi]         * this->_pData[this->map_index ( 1,k ) ] +
          src[k-_wi+1]       * this->_pData[this->map_index ( 2,k ) ] +
          src[k-1]           * this->_pData[this->map_index ( 3,k ) ] +
          src[k]             * this->_pData[this->map_index ( 4,k ) ] +
          src[k+1]           * this->_pData[this->map_index ( 5,k ) ] +
          src[k+_wi-1]       * this->_pData[this->map_index ( 6,k ) ] +
          src[k+_wi]         * this->_pData[this->map_index ( 7,k ) ] +
          src[k+_wi+1]       * this->_pData[this->map_index ( 8,k ) ];
      }
    }

  }

protected:
  void init ( const int wi, const int he ) {
    _wi = wi;
    _he = he;

    aol::Vector< int > offsets ( 9 );
    offsets[0] = -_wi - 1;
    offsets[1] = -_wi    ;
    offsets[2] = -_wi + 1;
    offsets[3] =       -1;
    offsets[4] =        0;
    offsets[5] =        1;
    offsets[6] =  _wi - 1;
    offsets[7] =  _wi    ;
    offsets[8] =  _wi + 1;

    if ( ( _wi < 4 ) || ( _he < 4 ) ) {
      throw aol::Exception ( "MultilinFEBandMatrix<DataType, qc::QC_2D>: size too small. Will not work and would not be useful.", __FILE__, __LINE__ );
    } else {
      this->reallocate ( _wi*_he, _wi*_he, offsets );
    }
  }

};

/** \brief A class for band matrices for use with 27-point stencils (e. g. in 3D trilinear finite element computations)
 *  \author Schwen
 *  \ingroup Matrix
 */
template<typename DataType>
class MultilinFEBandMatrix< DataType, qc::QC_3D > : public aol::GenBandMatrix<DataType> {

public:
  typedef BitArray<qc::QC_3D> MaskType;


  explicit MultilinFEBandMatrix ( const GridStructure &Struc ) : aol::GenBandMatrix<DataType> ( Struc.getNumberOfNodes(), Struc.getNumberOfNodes() ) {
    init ( Struc.getNumX(), Struc.getNumY(), Struc.getNumZ() );
  }

  explicit MultilinFEBandMatrix ( const GridSize<QC_3D> &gridSize ) : aol::GenBandMatrix<DataType> ( gridSize.getNumberOfNodes(), gridSize.getNumberOfNodes() ) {
    init ( gridSize.getNumX(), gridSize.getNumY(), gridSize.getNumZ() );
  }

protected:
  void init ( const int nx, const int ny, const int nz ) {
    aol::Vector<int> offsets ( 27 );
    offsets[ 0] = qc::ILexCombine3 ( -1, -1, -1, nx, ny );
    offsets[ 1] = qc::ILexCombine3 (  0, -1, -1, nx, ny );
    offsets[ 2] = qc::ILexCombine3 ( +1, -1, -1, nx, ny );
    offsets[ 3] = qc::ILexCombine3 ( -1,  0, -1, nx, ny );
    offsets[ 4] = qc::ILexCombine3 (  0,  0, -1, nx, ny );
    offsets[ 5] = qc::ILexCombine3 ( +1,  0, -1, nx, ny );
    offsets[ 6] = qc::ILexCombine3 ( -1, +1, -1, nx, ny );
    offsets[ 7] = qc::ILexCombine3 (  0, +1, -1, nx, ny );
    offsets[ 8] = qc::ILexCombine3 ( +1, +1, -1, nx, ny );
    offsets[ 9] = qc::ILexCombine3 ( -1, -1,  0, nx, ny );
    offsets[10] = qc::ILexCombine3 (  0, -1,  0, nx, ny );
    offsets[11] = qc::ILexCombine3 ( +1, -1,  0, nx, ny );
    offsets[12] = qc::ILexCombine3 ( -1,  0,  0, nx, ny );
    offsets[13] = qc::ILexCombine3 (  0,  0,  0, nx, ny );
    offsets[14] = qc::ILexCombine3 ( +1,  0,  0, nx, ny );
    offsets[15] = qc::ILexCombine3 ( -1, +1,  0, nx, ny );
    offsets[16] = qc::ILexCombine3 (  0, +1,  0, nx, ny );
    offsets[17] = qc::ILexCombine3 ( +1, +1,  0, nx, ny );
    offsets[18] = qc::ILexCombine3 ( -1, -1, +1, nx, ny );
    offsets[19] = qc::ILexCombine3 (  0, -1, +1, nx, ny );
    offsets[20] = qc::ILexCombine3 ( +1, -1, +1, nx, ny );
    offsets[21] = qc::ILexCombine3 ( -1,  0, +1, nx, ny );
    offsets[22] = qc::ILexCombine3 (  0,  0, +1, nx, ny );
    offsets[23] = qc::ILexCombine3 ( +1,  0, +1, nx, ny );
    offsets[24] = qc::ILexCombine3 ( -1, +1, +1, nx, ny );
    offsets[25] = qc::ILexCombine3 (  0, +1, +1, nx, ny );
    offsets[26] = qc::ILexCombine3 ( +1, +1, +1, nx, ny );

    if ( ny < 4 || nx < 4 ) { // this may not be a sharp bound
      throw aol::Exception ( "MultilinFEBandMatrix<qc::QC_3D>: size too small. Will not work and would not be useful.", __FILE__, __LINE__ );
    } else {
      this->reallocate ( nx*ny*nz, nx*ny*nz, offsets );
    }
  }
};

//! Basis class for band matrices for affine FE on quoc grids
/**
 * \ingroup Matrix
 */
template <class DataTye, qc::Dimension dim>
class AffineFEBandMatrix {};


/** \brief A class for band matrices for use with 15-point stencils (e. g. in 3D affine finite element or CFE computations).
 *  Note that there is no canonical 15-point stencil.
 *  \author Schwen
 *  \ingroup Matrix
 */
template<typename DataType >
class AffineFEBandMatrix< DataType, qc::QC_3D > : public aol::GenBandMatrix<DataType> {

public:
  explicit AffineFEBandMatrix ( const GridStructure &grid ) : aol::GenBandMatrix<DataType> ( grid.getNumberOfNodes(), grid.getNumberOfNodes() ) {

    const int nx = grid.getNumX(), ny = grid.getNumY(), nz = grid.getNumZ();

    aol::Vector<int> offsets ( 15 );
    offsets[ 0] = qc::ILexCombine3 (  0, -1, -1, nx, ny );
    offsets[ 1] = qc::ILexCombine3 ( +1, -1, -1, nx, ny );
    offsets[ 2] = qc::ILexCombine3 ( -1,  0, -1, nx, ny );
    offsets[ 3] = qc::ILexCombine3 (  0,  0, -1, nx, ny );
    offsets[ 4] = qc::ILexCombine3 (  0, -1,  0, nx, ny );
    offsets[ 5] = qc::ILexCombine3 ( +1, -1,  0, nx, ny );
    offsets[ 6] = qc::ILexCombine3 ( -1,  0,  0, nx, ny );
    offsets[ 7] = qc::ILexCombine3 (  0,  0,  0, nx, ny );
    offsets[ 8] = qc::ILexCombine3 ( +1,  0,  0, nx, ny );
    offsets[ 9] = qc::ILexCombine3 ( -1, +1,  0, nx, ny );
    offsets[10] = qc::ILexCombine3 (  0, +1,  0, nx, ny );
    offsets[11] = qc::ILexCombine3 (  0,  0, +1, nx, ny );
    offsets[12] = qc::ILexCombine3 ( +1,  0, +1, nx, ny );
    offsets[13] = qc::ILexCombine3 ( -1, +1, +1, nx, ny );
    offsets[14] = qc::ILexCombine3 (  0, +1, +1, nx, ny );

    if ( ny < 4 || nx < 4 ) { // this may not be a sharp bound
      throw aol::Exception ( "AffineFEBandMatrix: size too small. Will not work and would not be useful.", __FILE__, __LINE__ );
    } else {
      this->reallocate ( nx*ny*nz, nx*ny*nz, offsets );
    }
  }

};



//! A special CSR matrix type for the sparsity structure resulting from quoc grids
/**
 * \ingroup Matrix
 */
template <Dimension Dim>
class UniGridCSR_Matrix { };

template <>
class UniGridCSR_Matrix<QC_2D> : public aol::CSR_Matrix<> {
public:

  explicit UniGridCSR_Matrix ( const GridStructure &Grid ) {
    init( Grid.getNumX(), Grid.getNumY() );
    this->setZero();
  }

  explicit UniGridCSR_Matrix ( const GridSize<QC_2D> &GridSize ) {
    init( GridSize.getNumX(), GridSize.getNumY() );
    this->setZero();
  }

protected:
  void init( const int NumX, const int NumY ) {
    const int size = NumX * NumY;

    const int non_zeroes = ( NumX - 2 ) * ( NumY - 2 ) * 9 + 2 * ( NumX - 2 ) * 6 + 2 * ( NumY - 2 ) * 6 + 4 * 4;
    int count = 0;

    this->_numRows = this->_numCols = size;
    this->_indPointer.reallocate ( this->_numRows + 1 );
    this->_value.reallocate ( non_zeroes );
    this->_index.reallocate ( non_zeroes );
    this->_indPointer[0] = 0;
    int counter = 0;
    for ( int RowNum = 0; RowNum < size; ++RowNum ) {
      int startIndex = RowNum - NumX - 1;
      int y = RowNum / NumX;
      int x = RowNum % NumX;
      int numInRow = 0;
      for ( int j = -1; j <= 1; ++j ) {
        int Y = y + j;
        for ( int i = -1; i <= 1; ++i ) {
          int X = x + i;
          if ( X >= 0 && X < NumX && Y >= 0 && Y < NumY ) {
            this->_index[counter++] = startIndex;
            ++numInRow;
            // std::cerr <<  startIndex << " ";
            ++count;
          }
          ++startIndex;
        }
        startIndex += NumX - 3;
      }
      this->_indPointer[RowNum+1] =  this->_indPointer[RowNum] + numInRow;
      //cerr << endl;
    }
    // cerr << "count = " << count << " non_zeroes = " << non_zeroes << endl;
  }
};


template <>
class UniGridCSR_Matrix<QC_3D> : public aol::CSR_Matrix<> {
public:

  explicit UniGridCSR_Matrix ( const GridStructure &Grid ) {
    init3d ( Grid.getNumX(), Grid.getNumY(), Grid.getNumZ() );
    this->setZero();
  }

  explicit UniGridCSR_Matrix ( const GridSize<QC_3D> &GridSize ) {
    init3d ( GridSize.getNumX(), GridSize.getNumY(), GridSize.getNumZ() );
    this->setZero();
  }

protected:
  //! note this has not been tested yet with numX!=numY!=numZ
  void init3d ( int numX, int numY, int numZ ) {
#ifdef VERBOSE
    cerr << "numZ = " << numZ << endl;
#endif
    const int size = numX * numY * numZ;

    const int non_zeroes = ( numX - 2 ) * ( numY - 2 ) * ( numZ - 2 ) * 27
                           + 2 * ( numX - 2 ) * ( numY - 2 ) * 18 +
                           + 2 * ( numX - 2 ) * ( numZ - 2 ) * 18 +
                           + 2 * ( numY - 2 ) * ( numZ - 2 ) * 18 +
                           + 4 * ( numX - 2 ) * 12 +
                           + 4 * ( numY - 2 ) * 12 +
                           + 4 * ( numZ - 2 ) * 12 +
                           + 8 * 8;
    int count = 0;

    this->_numRows = this->_numCols = size;
    this->_indPointer.resize ( this->_numRows + 1 );
    this->_value.resize ( non_zeroes );
    this->_index.resize ( non_zeroes );
    this->_indPointer[0] = 0;
    int counter = 0;

    int RowNum = 0;
    for ( int z = 0; z < numZ; ++z ) {
      for ( int y = 0; y < numY; ++y ) {
        for ( int x = 0; x < numX; ++x ) {
          int startIndex = RowNum - numX * numY - numX - 1;
          int numInRow = 0;
          for ( int k = -1; k <= 1; ++k ) {
            int Z = z + k;
            for ( int j = -1; j <= 1; ++j ) {
              int Y = y + j;
              for ( int i = -1; i <= 1; ++i ) {
                int X = x + i;
                if ( X >= 0 && X < numX && Y >= 0 && Y < numY && Z >= 0 && Z < numZ ) {
                  this->_index[counter++] = startIndex;
                  ++numInRow;
                  // std::cerr <<  startIndex << " ";
                  ++count;
                }
                ++startIndex;
              }
              startIndex += numX - 3;
            }
            startIndex += numX * numY - 2 * numX - 2 - ( numX - 3 + 1 ) ;
          }
          this->_indPointer[RowNum+1] = this->_indPointer[RowNum] + numInRow;
          ++RowNum;
          // cerr << endl;
        }
      }
    }
#ifdef VERBOSE
    cerr << "count = " << count << " non_zeroes = " << non_zeroes << endl;
#endif
  }

};

}

#endif
