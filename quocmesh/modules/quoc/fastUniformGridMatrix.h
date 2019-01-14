#ifndef __FASTUNIFORMGRIDMATRIX_H
#define __FASTUNIFORMGRIDMATRIX_H

#include <quoc.h>
#include <vec.h>
#include <gridBase.h>
#include <simplexGrid.h>
#include <rows.h>
#include <sparseMatrices.h>
#include <bitArray.h>

namespace qc {

// forward declaration of CellCenteredCubicGrid
template <qc::Dimension Dim> class CellCenteredCubicGrid;

template <class _DataType, Dimension Dim, typename BaseClass = aol::GenSparseOp<_DataType> >
class FastUniformGridMatrix {
public:
  virtual ~FastUniformGridMatrix() {}
};

enum FastAssembleUniformGridMatrixApplyMode { ROWWISE,
                                              BANDWISE,
                                              BLOCK3WISE_SEGMENTED,
                                              BLOCK3WISE,
                                              CAREFUL_EVERYWHERE };

template <class _DataType, FastAssembleUniformGridMatrixApplyMode apply_mode = ROWWISE >
class FastAssembleUniformGridMatrix : public aol::Op<aol::Vector<_DataType> > {
  const int _w;
  const int _size;
public:
  typedef _DataType DataType;
  explicit FastAssembleUniformGridMatrix ( const qc::GridDefinition &Grid )
      : _w ( Grid.getWidth() ), _size ( Grid.getNumberOfNodes() ) {
    if ( Grid.getDimOfWorld() != QC_2D ) {
      throw aol::Exception ( "hey, you try to give a 3d grid to a 2d-version of fast uniform grid matrix.\n", __FILE__, __LINE__ );
    }
    if ( Grid.getGridDepth() == 0 ) {
      throw aol::Exception ( "FastUniformGridMatrix doesn't work for level 0.\n", __FILE__, __LINE__ );
    }
    _rows = new aol::Vec<9, DataType> [ _size ];
    setZero();  // be tidy when born
  }

  explicit FastAssembleUniformGridMatrix ( const GridSize<QC_2D> &GridSize )
      : _w ( GridSize.getNumX() ), _size ( GridSize.getNumberOfNodes() ) {
    // check if given grid size is equally dimensioned.
    GridSize.quadraticOrDie ();
    if ( GridSize.getNumX() < 3 ) {
      throw aol::Exception ( "FastUniformGridMatrix doesn't work for level 0.\n", __FILE__, __LINE__ );
    }
    _rows = new aol::Vec<9, DataType> [ _size ];
    setZero();  // be tidy when born
  }

  FastAssembleUniformGridMatrix ( const int NumRows, const int NumCols )
      : _w ( static_cast<int>( sqrt( NumRows ) ) ), _size ( NumRows ) {
    // check if given grid size is equally dimensioned.
    if ( NumRows != NumCols  ) {
      throw aol::Exception ( "Number of rows must be equal to number of columns.\n", __FILE__, __LINE__ );
    }
    if ( _w*_w != NumRows ) {
      throw aol::Exception ( "Numbers of rows must be square.\n", __FILE__, __LINE__ );
    }
    _rows = new aol::Vec<9, DataType> [ _size ];
    setZero();  // be tidy when born
  }

  virtual ~FastAssembleUniformGridMatrix() {
    delete [] _rows;
  }

  int getNumCols() const {
    return _size;
  }

  int getNumRows() const {
    return _size;
  }

  void setZero() {
    // use memset? may be dangerous..
    for ( int i = 0; i < _size; ++i )  {
      _rows[i].setZero();
    }
  }

  inline DataType get ( int I, int J ) const {
    return _rows[I][ getIndex ( J-I ) ];
  }

  inline void set ( int I, int J, DataType Value ) {
    _rows[I][ getIndex ( J-I ) ] = Value;
  }

  inline void add ( int I, int J, DataType Value ) {
    _rows[I][ getIndex ( J-I ) ] += Value;
  }

  inline DataType getDiag( int I ) const {
    return _rows[I][4];
  }

  FastAssembleUniformGridMatrix<DataType,apply_mode>& addMultiple ( const FastAssembleUniformGridMatrix<DataType,apply_mode> &other, const DataType factor ) {
    if ( _size != other._size ) {
      throw aol::Exception ( "Sizes of matrices do not match", __FILE__, __LINE__ );
    }
    for ( int i = 0; i < _size; ++i ) {
      _rows[i].addMultiple( other._rows[i], factor );
    }
    return *this;
  }

  FastAssembleUniformGridMatrix<DataType,apply_mode>& operator+= ( const FastAssembleUniformGridMatrix<DataType,apply_mode> &other ) {
    addMultiple ( other, aol::NumberTrait<DataType>::one );
    return *this;
  }


  void applyAdd ( const aol::Vector<DataType> &ArgVec, aol::Vector<DataType> &DestVec ) const {
    DataType * Arg  = ArgVec.getData();
    DataType * Dest = DestVec.getData();

    for ( int i = 0; i <= _w + 1; ++i ) {
      multCarefullyAtBoundary ( i, Arg, Dest );
    }

    // try something "cache-optimal" here
    switch ( apply_mode ) {
      case CAREFUL_EVERYWHERE: {
        for ( int i = _w + 2; i < _size - _w - 2; ++i ) {
          multCarefullyAtBoundary ( i, Arg, Dest );
        }
      } break;

      case BANDWISE: {
        int l = 0;
        int globos = -_w - 1;
        for ( int k = 0; k < 3; ++k ) {
          for ( int j = 0; j < 3; ++j ) {
            int g = globos + _w + 2;
            for ( int i = _w + 2; i < _size - _w - 2; ++i ) {
              // cerr << "inner row = " << i << endl;
              Dest[i] += Arg[g++] * _rows[i][l];
            }
            ++l;
            ++globos;
          }
          globos += _w - 3;
        }
      } break;

      case BLOCK3WISE_SEGMENTED: {
        const int seg_len = 100;
        for ( int seg = _w + 2; seg < _size - _w - 2; seg += seg_len ) {
          const int seg_end = aol::Min ( seg + seg_len, _size - _w - 2 );
          for ( int k = 0; k < 3; ++k ) {
            int locos = 3 * k;
            int globos = _w * ( k - 1 ) - 1;
            for ( int i = seg; i < seg_end; ++i ) {
              // cerr << "inner row = " << i << endl;
              int l = locos, g = globos + i;
              for ( int j = 0; j < 3; ++j ) {
                Dest[i] += Arg[g++] * _rows[i][l++];
              }
            }
          }
        }
      } break;

      case BLOCK3WISE: {
        for ( int k = 0; k < 3; ++k ) {
          int locos = 3 * k;
          int globos = _w * ( k - 1 ) - 1;
          for ( int i = _w + 2; i < _size - _w - 2; ++i ) {
            // cerr << "inner row = " << i << endl;
            int l = locos, g = globos + i;
            for ( int j = 0; j < 3; ++j ) {
              Dest[i] += Arg[g++] * _rows[i][l++];
            }
          }
        }
      } break;

      case ROWWISE: {
        for ( int i = _w + 2; i < _size - _w - 2; ++i ) {
          int l = 0, g = i - _w - 1;
          for ( int k = 0; k < 3; ++k ) {
            for ( int j = 0; j < 3; ++j ) {
              Dest[i] += Arg[g++] * _rows[i][l++];
            }
            g += _w - 3;
          }
        }
      } break;

      default:
      throw aol::Exception ( "FastAssembleUniformGridMatrix::applyAdd: illegal multiplication mode", __FILE__, __LINE__ );
    }

    for ( int i = aol::Max ( _w + 2, _size - _w - 2 ); i < _size; ++i ) {
      multCarefullyAtBoundary ( i, Arg, Dest );
    }

  }

  //! Return vector of row entries. Entries need not be sorted with respect to column index and zeros may be contained.
  void makeRowEntries ( vector<typename aol::Row<DataType>::RowEntry > &vec, const int RowNum ) const {
    int startIndex = RowNum - _w - 1;

    unsigned int n = 0;
    const unsigned int vec_size = vec.size();
    int k = 0;

    const int y = RowNum / _w;
    const int x = RowNum % _w;

    for ( int Y = y - 1; Y <= y + 1; ++Y ) {
      if ( Y >= 0 && Y < _w ) {
        for ( int X = x - 1; X <= x + 1; ++X ) {
          if ( X >= 0 && X < _w ) {
            // cerr << "row = " << row << " index = " << startIndex << " matrixentry = " <<  _rows[row][k] << endl;
            if ( n < vec_size ) {
              vec[n].col = startIndex;
              vec[n].value = _rows[RowNum][k];
            } else {
              vec.push_back ( typename aol::Row<DataType>::RowEntry ( startIndex, _rows[RowNum][k] ) );
            }
            ++n;
          }
          ++startIndex;
          ++k;
        }
      } else {
        k += 3;
        startIndex += 3;
      }
      startIndex += _w - 3;
    }
    if ( n < vec.size() ) vec.resize ( n );
  }

  //! this method is broken and will not compile
  FastAssembleUniformGridMatrix<DataType, apply_mode>& operator*= ( const DataType factor ) {
    for ( int i = 0; i < _size; ++i )  {
      _rows[i] *= ( factor );
    }
    return ( *this );
  }

  //! this method is broken and will not compile
  void scaleRow ( const int row, const DataType factor ) {
    _rows[ row ].scale ( factor );
  }

protected:
  void multCarefullyAtBoundary ( int row, const DataType * Arg, DataType * Dest ) const {
    int startIndex = row - _w - 1;
    int k = 0;

    int y = row / _w;
    int x = row % _w;

    for ( int j = -1; j <= 1; ++j ) {
      int Y = y + j;
      for ( int i = -1; i <= 1; ++i ) {
        int X = x + i;
        if ( X >= 0 && X < _w && Y >= 0 && Y < _w ) {
          // cerr << "row = " << row << " index = " << startIndex << " matrixentry = " <<  _rows[row][k] << endl;
          Dest[row] += _rows[row][k] * Arg[ startIndex ];
        }
        ++startIndex;
        ++k;
      }
      startIndex += _w - 3;
    }
  }

  void testGetIndex() const {
    int offsets[9] = { -_w - 1, -_w, -_w + 1, -1, 0, 1, _w - 1, _w, _w + 1 };
    cerr << "test get index: ";
    for ( int i = 0; i < 9; ++i ) {
      cerr << getIndex ( offsets[i] ) << " ";
    }
    cerr << endl;
  }

  inline int getIndex ( int offset ) const {
    if ( offset < -1 ) {
      return offset + _w + 1;
    } else if ( offset > 1 ) {
      return offset - _w + 7;
    }
    return offset + 4;
  }

  aol::Vec<9, DataType> * _rows;

public:
  const aol::Vec<9, DataType> * getInternalDataReference ( ) const {
    return _rows;
  }
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * \ingroup Matrix
 */
template <class _DataType, typename BaseClass>
class FastUniformGridMatrix<_DataType, qc::QC_2D, BaseClass> : public BaseClass {
  const int _w;
  const int _size;
public:
  typedef _DataType DataType;

  explicit FastUniformGridMatrix ( const qc::GridDefinition &Grid )
      : BaseClass(Grid.getNumberOfNodes(), Grid.getNumberOfNodes()), _w ( Grid.getNumX() ), _size ( Grid.getNumberOfNodes() ) {
    if ( Grid.getDimOfWorld() != QC_2D ) {
      throw aol::Exception ( "hey, you try to give a 3d grid to a 2d-version of fast uniform grid matrix.\n", __FILE__, __LINE__ );
    }
    if ( Grid.getGridDepth() == 0 ) {
      throw aol::Exception ( "FastUniformGridMatrix doesn't work for level 0.\n", __FILE__, __LINE__ );
    }
    for ( int i = 0; i < 3; ++i )
      _rows[i] = new aol::Vec<3, DataType> [ _size ];
    setZero();  // be tidy when born
  }

  template <typename CubicGridType>
  explicit FastUniformGridMatrix ( const qc::simplex::GridStructure<CubicGridType, QC_2D> & Grid )
      : BaseClass(Grid.getNumberOfNodes(), Grid.getNumberOfNodes()), _w ( Grid.getNumX() ), _size ( Grid.getNumberOfNodes() ) {
    if ( Grid.getGridDepth() == 0 ) {
      throw aol::Exception ( "FastUniformGridMatrix doesn't work for level 0.\n", __FILE__, __LINE__ );
    }
    for ( int i = 0; i < 3; ++i )
      _rows[i] = new aol::Vec<3, DataType> [ _size ];
    setZero();  // be tidy when born
  }

  explicit FastUniformGridMatrix ( const qc::CellCenteredCubicGrid<qc::QC_2D> &Grid );

  explicit FastUniformGridMatrix ( const qc::GridSize<QC_2D> & GridSize )
      : BaseClass(GridSize.getNumberOfNodes(), GridSize.getNumberOfNodes()), _w ( GridSize.getNumX() ), _size ( GridSize.getNumberOfNodes() ) {
    // check if given grid size is equally dimensioned.
    GridSize.quadraticOrDie ();
    for ( int i = 0; i < 3; ++i )
      _rows[i] = new aol::Vec<3, DataType> [ _size ];
    setZero();  // be tidy when born
  }

  FastUniformGridMatrix ( int NumRows, int NumCols )
      : BaseClass(NumRows, NumCols), _w ( static_cast<int>( sqrt( static_cast<double> ( NumRows ) ) ) ), _size ( NumRows ) {
    if ( NumRows != NumCols  ) {
      throw aol::Exception ( "Number of rows must be equal to number of columns.\n", __FILE__, __LINE__ );
    }
    if ( _w*_w != NumRows ) {
      throw aol::Exception ( "Numbers of rows must be square.\n", __FILE__, __LINE__ );
    }
    for ( int i = 0; i < 3; ++i )
      _rows[i] = new aol::Vec<3, DataType> [ _size ];
    setZero();  // be tidy when born
  }

  virtual ~FastUniformGridMatrix() {
    for ( int i = 0; i < 3; ++i )
      delete [] _rows[i];
  }

  // WARNING: this does not make any compatibility checks.
  // the argument must have a compatible matrix structure, otherwise this function won't work.
  FastUniformGridMatrix<DataType, qc::QC_2D, BaseClass>& operator+= ( const aol::Matrix<DataType> &Matrix ) {
    vector<typename aol::Row<DataType>::RowEntry > vec;
    for ( int i = 0; i < _size; ++i ) {
      Matrix.makeRowEntries ( vec, i );
      for ( typename vector<typename aol::Row<DataType>::RowEntry >::iterator it = vec.begin();
            it != vec.end(); ++it ) {
        add ( i, it->col, it->value );
      }
    }
    return *this;
  }

  void setZero() {
    // use memset? may be dangerous..
    for ( int k = 0; k < 3; ++k ) {
      for ( int i = 0; i < _size; ++i )  {
        _rows[k][i].setZero();
      }
    }
  }

  DataType getMaxValue() {
    DataType max = 0;
    for ( int k = 0; k < 3; ++k )
      for ( int i = 0; i < _size; ++i )
        max = aol::Max( max, _rows[k][i].getMaxValue() );
    return max;
  }

  DataType getMinValue() {
    DataType min = 0;
    for ( int k = 0; k < 3; ++k )
      for ( int i = 0; i < _size; ++i )
        min = aol::Min( min, _rows[k][i].getMinValue() );
    return min;
  }

  inline DataType get ( int I, int J ) const {
    NON_PARALLEL_STATIC int block, index;
    getBlockAndIndex ( J - I, block, index );
    if ( block == -1 || index == -1 ) return 0.;
    return _rows[block][I][ index ];
  }

  inline void set ( int I, int J, DataType Value ) {
    NON_PARALLEL_STATIC int block, index;
    getBlockAndIndex ( J - I, block, index );
#ifdef BOUNDS_CHECK
    // This bounds check is sufficient although after the assembly not all entries of the diagonals are filled.
    // But they are stored!
    if ( block == -1 || index == -1 ) throw aol::OutOfBoundsException( "FastUniformGridMatrix2d: set: Indices out of bound!", __FILE__, __LINE__ );;
#endif
    _rows[block][I][ index ] = Value;
  }

  inline void add ( int I, int J, DataType Value ) {
    NON_PARALLEL_STATIC int block, index;
    getBlockAndIndex ( J - I, block, index );
#ifdef BOUNDS_CHECK
    if ( block == -1 || index == -1 ) throw aol::OutOfBoundsException( "FastUniformGridMatrix2d: add: Indices out of bound!", __FILE__, __LINE__ );;
#endif
    _rows[block][I][ index ] += Value;
  }


  // the diagonal element is in the 2nd Vec3 and there it's the entry with index 1
  inline DataType getDiag( int I ) const {
    return _rows[1][I][1];
  }

  // the diagonal element is in the 2nd Vec3 and there it's the entry with index 1
  inline void setDiag( int I, DataType Value ) {
    _rows[1][I][1] = Value;
  }

  void applyAdd ( const aol::Vector<DataType> &Arg, aol::Vector<DataType> &Dest ) const;

  void setRowColToZero ( const int index ) {
    int i;

    if ( index >= 0 && index < _size ) {
      // Delete row
      for ( i = 0; i < 3; ++i )
        _rows[i][index].setZero();

      // Delete column
      i = index - 1;
      if ( i >= 0 ) _rows[1][i][2] = 0.0;
      i = index - _w + 1;
      if ( i >= 0 ) _rows[2][i][0] = 0.0;
      i = index - _w;
      if ( i >= 0 ) _rows[2][i][1] = 0.0;
      i = index - _w - 1;
      if ( i >= 0 ) _rows[2][i][2] = 0.0;

      i = index + 1;
      if ( i < _size ) _rows[1][i][0] = 0.0;
      i = index + _w - 1;
      if ( i < _size ) _rows[0][i][2] = 0.0;
      i = index + _w;
      if ( i < _size ) _rows[0][i][1] = 0.0;
      i = index + _w + 1;
      if ( i < _size ) _rows[0][i][0] = 0.0;
    }
  }

  void setRowToZero ( const int index ) {
    int i;

    if ( index >= 0 && index < _size ) {
      // Delete row
      for ( i = 0; i < 3; i++ )
        _rows[i][index].setZero();
    }
  }

  void setColToZero ( const int index ) {
    int i;

    if ( index >= 0 && index < _size ) {
      // Delete column
      i = index - 1;
      if ( i >= 0 ) _rows[1][i][2] = 0.0;
      i = index - _w + 1;
      if ( i >= 0 ) _rows[2][i][0] = 0.0;
      i = index - _w;
      if ( i >= 0 ) _rows[2][i][1] = 0.0;
      i = index - _w - 1;
      if ( i >= 0 ) _rows[2][i][2] = 0.0;

      i = index + 1;
      if ( i < _size ) _rows[1][i][0] = 0.0;
      i = index + _w - 1;
      if ( i < _size ) _rows[0][i][2] = 0.0;
      i = index + _w;
      if ( i < _size ) _rows[0][i][1] = 0.0;
      i = index + _w + 1;
      if ( i < _size ) _rows[0][i][0] = 0.0;
    }
  }

  virtual ostream& print ( ostream &out ) const {
    out << "Size is " << _size << std::endl;
    for ( int i = 0; i < _size; ++i ) {
      for ( int k = 0; k < 3; ++k ) {
        aol::Vec<3, DataType> r = _rows[k][i];
        for ( int l = 0; l < 3; ++l ) out << r[l] << " ";
      }
      out << std::endl;
    }
    return out;
  }

  FastUniformGridMatrix<DataType, qc::QC_2D, BaseClass>& operator*= ( const DataType factor ) {
    for ( int k = 0; k < 3; ++k ) {
      for ( int i = 0; i < _size; ++i ) {
        _rows[k][i] *= factor;
      }
    }
    return ( *this );
  }

  void scaleRow ( const int row, const DataType factor ) {
    for ( int k = 0; k < 3; ++k ) {
      _rows[k][row] *= factor;
    }
  }

  //! Return vector of row entries. Entries need not be sorted with respect to column index and zeros may be contained.
  void makeRowEntries ( vector<typename aol::Row<DataType>::RowEntry > &vec, const int RowNum ) const {
    int startIndex = RowNum - _w - 1;

    if ( RowNum <= _w + 1 || RowNum >= _size - _w - 2 ) {
      int n = 0;
      int k = 0, b = 0;
      int y = RowNum / _w;
      int x = RowNum % _w;

      for ( int j = -1; j <= 1; ++j ) {
        int Y = y + j;
        for ( int i = -1; i <= 1; ++i ) {
          int X = x + i;
          if ( X >= 0 && X < _w && Y >= 0 && Y < _w ) {
            // cerr << "row = " << row << " index = " << startIndex << " matrixentry = " <<  _rows[row][k] << endl;
            if ( n < static_cast< int > ( vec.size() ) ) {
              vec[n].col = startIndex;
              vec[n].value = _rows[b][RowNum][k];
            } else {
              vec.push_back ( typename aol::Row<DataType>::RowEntry ( startIndex, _rows[b][RowNum][k] ) );
            }
            ++n;
          }
          ++startIndex;
          ++k;
        }
        startIndex += _w - 3;
        ++b;
        k = 0;
      }
      vec.resize ( n );
    } else {
      int k = 0, b = 0;
      int n = 0;
      for ( int i = 0; i < 3; ++i ) {
        for ( int j = 0; j < 3; ++j ) {
          if ( startIndex >= 0 ) {
            if ( startIndex >= _size ) {
              vec.resize ( n );
              return;
            }
            if ( n < static_cast< int >(  vec.size() ) ) {
              vec[n].col = startIndex;
              vec[n].value = _rows[b][RowNum][k];
            } else {
              vec.push_back ( typename aol::Row<DataType>::RowEntry ( startIndex, _rows[b][RowNum][k] ) );
            }
            ++n;
          }
          ++startIndex;
          ++k;
        }
        startIndex += _w - 3;
        ++b;
        k = 0;
      }
      vec.resize ( n );
    }
  }

  template <typename OtherBaseClass>
  FastUniformGridMatrix<DataType, qc::QC_2D, BaseClass>& addMultiple ( const FastUniformGridMatrix<DataType, qc::QC_2D, OtherBaseClass> &other, const DataType factor ) {
    if ( _size != other.getNumRows() ) {
      throw aol::Exception ( "Sizes of matrices do not match", __FILE__, __LINE__ );
    }
    for ( int b = 0; b < 3; ++b ) {
      for ( int i = 0; i < _size; ++i ) {
        _rows[b][i] += factor * other.getInternalDataReference ( b) [i];
      }
    }
    return *this;
  }

  template <typename OtherBaseClass>
  FastUniformGridMatrix<DataType, qc::QC_2D, BaseClass>& operator+= ( const FastUniformGridMatrix<DataType, qc::QC_2D, OtherBaseClass> &other ) {
    addMultiple ( other, aol::NumberTrait<DataType>::one );
    return *this;
  }

  FastUniformGridMatrix<DataType, qc::QC_2D, BaseClass>& operator= ( const FastUniformGridMatrix<DataType, qc::QC_2D, BaseClass> &other ) {
    if ( _size != other._size ) {
      throw aol::Exception ( "Sizes of matrices do not match", __FILE__, __LINE__ );
    }
    for ( int b = 0; b < 3; ++b ) {
      for ( int i = 0; i < _size; ++i ) {
        _rows[b][i] = other._rows[b][i];
      }
    }
    return *this;
  }

  bool isApproxEqual ( const FastUniformGridMatrix<DataType, qc::QC_2D, BaseClass> &Other, DataType Epsilon ) {
    if ( _size != Other._size ) {
      throw aol::Exception ( "Sizes of matrices do not match", __FILE__, __LINE__ );
    }

    for ( int b = 0; b < 3; ++b ) {
      for ( int i = 0; i < _size; ++i ) {
        for ( unsigned int k = 0; k < 3; ++k ) {
          if ( fabs ( _rows[b][i][k] - Other._rows[b][i][k] ) > Epsilon ) {
            return false;
          }
        }
      }
    }
    return true;
  }

  //! Uses the general function aol::transposeAToB, therefore not speed optimized
  //! for the structure of this matrix class.
  void transposeTo ( FastUniformGridMatrix<DataType, qc::QC_2D, BaseClass> &Dest ) const {
    aol::transposeAToB<FastUniformGridMatrix<DataType, qc::QC_2D, BaseClass>, DataType> ( *this, Dest );
  }

  int numNonZeroes ( int I ) const {
    int nNonZeroes = 0;
    for ( int i = 0; i < 3; i++ ) {
      nNonZeroes += _rows[i][I].numNonZeroes();
    }
    return nNonZeroes;
  }

  int maxNumNonZeroesPerRow () const {
    return aol::maxNumNonZeroesPerRow<FastUniformGridMatrix<DataType, qc::QC_2D, BaseClass>, DataType> ( *this );
  }

protected:
  void multCarefullyAtBoundary ( int row, const aol::Vector<DataType> &ArgVec, aol::Vector<DataType> &DestVec ) const {
    DataType * Arg  = ArgVec.getData();
    DataType * Dest = DestVec.getData();

    int startIndex = row - _w - 1;
    int k = 0, b = 0;

    int y = row / _w;
    int x = row % _w;

    for ( int j = -1; j <= 1; ++j ) {
      int Y = y + j;
      for ( int i = -1; i <= 1; ++i ) {
        int X = x + i;

        if ( X >= 0 && X < _w && Y >= 0 && Y < _w ) {
          // cerr << "row = " << row << " index = " << startIndex << " matrixentry = " <<  _rows[row][k] << endl;
          Dest[row] += _rows[b][row][k] * Arg[ startIndex ];
        }
        ++startIndex;
        ++k;
      }
      k = 0;
      ++b;
      startIndex += _w - 3;
    }
  }

  void testGetIndex() const {
    int offsets[9] = { -_w - 1, -_w, -_w + 1, -1, 0, 1, _w - 1, _w, _w + 1 };
    cerr << "test get index: ";
    for ( int i = 0; i < 9; ++i ) {
      int block, index;
      getBlockAndIndex ( offsets[i], block, index );
      cerr << "offset = " <<  offsets[i] << " " << block << " "  << index << endl;
    }
    cerr << endl;
  }

  inline void getBlockAndIndex ( int offset, int &block, int &index ) const {
    if ( offset >= -_w-1 && offset <= -_w+1 ) {
      block = 0;
      index = offset + _w + 1;
    } else if ( offset >= _w-1 && offset <= _w+1 ) {
      block = 2;
      index = offset - _w + 1;
    } else if ( offset >= -1 && offset <= 1 ) {
      block = 1;
      index = offset + 1;
    } else {
      block = index = -1;
    }
  }

  aol::Vec<3, DataType> * _rows[3];

public:
  const aol::Vec<3, DataType> * getInternalDataReference ( const int Index ) const {
    return _rows[Index];
  }
};






/* ************************************************************************************
 *
 *              FASTUNIFORMGRIDMATRIX in 3D
 *
 * ************************************************************************************ */

//////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * \ingroup Matrix
 */
template <class _DataType, typename BaseClass>
class FastUniformGridMatrix<_DataType, qc::QC_3D, BaseClass> : public BaseClass {
  const int _w;       // width of the mesh in one direction
  const int _wsqr;    // width^2
  const int _size;    // the whole size

public:
  typedef _DataType DataType;

  explicit FastUniformGridMatrix ( const qc::GridDefinition &Grid )
      : BaseClass(Grid.getNumberOfNodes(), Grid.getNumberOfNodes()), _w ( Grid.getNumX() ), _wsqr ( _w*_w ), _size ( Grid.getNumberOfNodes() ) {
    if ( Grid.getDimOfWorld() != QC_3D ) {
      throw aol::Exception ( "hey, you try to give a n-dim. grid, n != 3,  to a 3d-version of fast uniform grid matrix.\n", __FILE__, __LINE__ );
    }
    if ( Grid.getGridDepth() == 0 ) {
      throw aol::Exception ( "FastUniformGridMatrix doesn't work for level 0.\n", __FILE__, __LINE__ );
    }
    // the array rows contains the values of the adjacent nodes
    // in 3D there are 9 of such rows
#ifdef VERBOSE
    cerr << "resizing vectors...";
#endif
    for ( int i = 0; i < 9; ++i )
      _rows[i] = new aol::Vec<3, DataType> [ _size ];
#ifdef VERBOSE
    cerr << "ready resizing vectors...\n";
#endif
    setZero();  // be tidy when born
  }

  template <typename CubicGridType>
  explicit FastUniformGridMatrix ( const qc::simplex::GridStructure<CubicGridType, QC_3D> & Grid )
      : BaseClass(Grid.getNumberOfNodes(), Grid.getNumberOfNodes()), _w ( Grid.getNumX() ), _wsqr ( _w*_w ), _size ( Grid.getNumberOfNodes() ) {
    if ( Grid.getGridDepth() == 0 ) {
      throw aol::Exception ( "FastUniformGridMatrix doesn't work for level 0.\n", __FILE__, __LINE__ );
    }
    // the array rows contains the values of the adjacent nodes
    // in 3D there are 9 of such rows
#ifdef VERBOSE
    cerr << "resizing vectors...";
#endif
    for ( int i = 0; i < 9; ++i )
      _rows[i] = new aol::Vec<3, DataType> [ _size ];
#ifdef VERBOSE
    cerr << "ready resizing vectors...\n";
#endif
    setZero();  // be tidy when born
  }

  explicit FastUniformGridMatrix ( const qc::GridSize<QC_3D> & GridSize )
      : BaseClass(GridSize.getNumberOfNodes(), GridSize.getNumberOfNodes()), _w ( GridSize.getNumX() ), _wsqr ( _w*_w ), _size ( GridSize.getNumberOfNodes() ) {
    // check if given grid size is equally dimensioned.
    GridSize.quadraticOrDie ();

    // the array rows contains the values of the adjacent nodes
    // in 3D there are 9 of such rows
#ifdef VERBOSE
    cerr << "resizing vectors...";
#endif
    for ( int i = 0; i < 9; ++i )
      _rows[i] = new aol::Vec<3, DataType> [ _size ];
#ifdef VERBOSE
    cerr << "ready resizing vectors...\n";
#endif
    setZero();  // be tidy when born
  }

  FastUniformGridMatrix ( int NumRows, int NumCols )
  : BaseClass(NumRows, NumCols), _w ( static_cast<int>( cbrt( static_cast<double> ( NumRows ) ) ) ), _wsqr ( _w*_w ), _size ( NumRows ) {
    if ( NumRows != NumCols  ) {
      throw aol::Exception ( "Number of rows must be equal to number of columns.\n", __FILE__, __LINE__ );
    }
    if ( _w*_w*_w != NumRows ) {
      throw aol::Exception ( "Numbers of rows must be cube.\n", __FILE__, __LINE__ );
    }
    for ( int i = 0; i < 9; ++i )
      _rows[i] = new aol::Vec<3, DataType> [ _size ];
    setZero();  // be tidy when born
  }

  virtual ~FastUniformGridMatrix() {
    for ( int i = 0; i < 9; ++i )
      delete [] _rows[i];
  }

  void setZero() {
    // use memset? may be dangerous..
    for ( int k = 0; k < 9; ++k ) {
      for ( int i = 0; i < _size; ++i )  {
        _rows[k][i].setZero();
      }
    }
  }

  DataType getMaxValue() {
    DataType max = 0;
    for ( int k = 0; k < 9; ++k )
      for ( int i = 0; i < _size; ++i )
        max = aol::Max( max, _rows[k][i].getMaxValue() );
    return max;
  }

  DataType getMinValue() {
    DataType min = 0;
    for ( int k = 0; k < 9; ++k )
      for ( int i = 0; i < _size; ++i )
        min = aol::Min( min, _rows[k][i].getMinValue() );
    return min;
  }

  void setRowColToZero ( const int index ) {
    int i, row, block, plane, blockInPlane, indexInBlock;

    if ( index >= 0 && index < _size ) {
      // Delete row
      for ( i = 0; i < 9; ++i )
        _rows[i][index].setZero();

      // Delete column
      for ( plane = -1; plane <= 1; ++plane )
        for ( blockInPlane = -1; blockInPlane <= 1; ++blockInPlane )
          for ( indexInBlock = -1; indexInBlock <= 1; ++ indexInBlock ) {
            block = 3 * plane + blockInPlane + 4;
            row = index - plane * _wsqr - blockInPlane * _w - indexInBlock;
            if ( ( row >= 0 ) && ( row < _size ) )
              _rows[block][row][indexInBlock+1] = 0.0;
          }
    }
  }

  // get and set entry of the matrix
  inline DataType get ( int I, int J ) const {
      NON_PARALLEL_STATIC int block, index;
      // if position (I,J) is illegal, then block and index will be set to -1
      getBlockAndIndex ( J - I, block, index );
      if ( block == -1 || index == -1 ) return 0.;
      else return _rows[block][I][ index ];
    }

  inline void set ( int I, int J, DataType Value ) {
    NON_PARALLEL_STATIC int block, index;
    // if position (I,J) is illegal, then block and index will be set to -1
    getBlockAndIndex ( J - I, block, index );
#ifdef BOUNDS_CHECK
    // This bounds check is sufficient although after the assembly not all entries of the diagonals are filled.
    // But they are stored!
    if ( block == -1 || index == -1 ) throw aol::OutOfBoundsException( "FastUniformGridMatrix3d: Set: Indices out of bound!", __FILE__, __LINE__ );;
#endif
    _rows[block][I][ index ] = Value;
  }

  inline void add ( int I, int J, DataType Value ) {
    NON_PARALLEL_STATIC int block, index;
    // if position (I,J) is illegal, then block and index will be set to -1
    getBlockAndIndex ( J - I, block, index );
#ifdef BOUNDS_CHECK
    if ( block == -1 || index == -1 ) throw aol::OutOfBoundsException( "FastUniformGridMatrix3d: Add: Indices out of bound!", __FILE__, __LINE__ );;
#endif
    _rows[block][I][ index ] += Value;
  }

  // the diagonal element is in the 4th Vec3 and there it's the entry with index 1
  inline DataType getDiag( int I ) const {
    return _rows[4][I][1];
  }

  // the diagonal element is in the 4th Vec3 and there it's the entry with index 1
  inline void setDiag( int I, DataType Value ) {
    _rows[4][I][1] = Value;
  }

  // -------------------------- The apply-functions ---------------------------------------

  // auxiliary function for subApplyAdd
  // applies a little section (_w elements) on the triple diagonal
  // (analogue to a 1D-apply)
  inline void tripleDiagApplyAdd ( const int startRow, const int startCol,
                                   const int tripDiag, const aol::Vector<DataType> &ArgVec, aol::Vector<DataType> &DestVec ) const {
    DataType * Arg  = ArgVec.getData();
    DataType * Dest = DestVec.getData();
    // the first row only consists of two mutliplications:
    Dest[startRow] += Arg[startCol] * _rows[tripDiag][startRow][1];
    Dest[startRow] += Arg[startCol + 1] * _rows[tripDiag][startRow][2];
    // now the sum of _w-1 products (inner nodes)
    for ( int i = 1; i < _w - 1; ++i )
      for ( int j = 0; j < 3; ++j )
        Dest[startRow + i] += Arg[startCol + i + j - 1] * _rows[tripDiag][startRow + i][j];
    // the last row only consists of two mutliplications too:
    Dest[startRow + _w-1] += Arg[startCol + _w-2 ] * _rows[tripDiag][startRow + _w-1][0];
    Dest[startRow + _w-1] += Arg[startCol + _w-1 ] * _rows[tripDiag][startRow + _w-1][1];
  }

  // auxiliary function for applyAdd
  // applys a smaller matrix to a piece of the arg-vector
  // args are the row of the dest-vec and the column of the
  // arg-vec, where the multiplication should start.
  // the bundle contains three triple diagonals.
  // This scheme is analogue to a 2D-apply!
  inline void subApplyAdd ( const int startRow, const int startCol, const int bundle,
                            const aol::Vector<DataType> &Arg, aol::Vector<DataType> &Dest ) const {
    // the first _w rows are special
    tripleDiagApplyAdd ( startRow, startCol, bundle*3 + 1, Arg, Dest );
    tripleDiagApplyAdd ( startRow, startCol + _w, bundle*3 + 2, Arg, Dest );
    // apply _w-2 times the triple diagonal form
    for ( int i = 1; i < _w - 1; ++i )
      for ( int j = 0; j < 3; ++j )
        tripleDiagApplyAdd ( startRow + i*_w, startCol + ( i + j - 1 ) *_w, bundle*3 + j, Arg, Dest );
    // the last _w rows are special again
    tripleDiagApplyAdd ( startRow + ( _w - 1 ) *_w, startCol + ( _w - 2 ) *_w, bundle*3, Arg, Dest );
    tripleDiagApplyAdd ( startRow + ( _w - 1 ) *_w, startCol + ( _w - 1 ) *_w, bundle*3 + 1, Arg, Dest );
  }

  void applyAdd ( const aol::Vector<DataType> &Arg, aol::Vector<DataType> &Dest ) const;

  // -------------------------- END APPLY --------------------------------------------


  void scaleRow ( const int row, const DataType factor ) {
    for ( int k = 0; k < 9; ++k ) {
      _rows[k][row] *= factor;
    }
  }

  int numNonZeroes ( int I ) const {
    int nNonZeroes = 0;
    for ( int i = 0; i < 9; i++ ) {
      nNonZeroes += _rows[i][I].numNonZeroes();
    }
    return nNonZeroes;
  }


  // -------------------------- START MAKEROWENTRIES ----------------------------------------

  // auxiliary function for makeRowEntries, fills the vector vec with the entries, that
  // are in referred _wsqr * _wsqr - block, localRowNum is the local index of the row in
  // such a block ( 0 <= subRowNum < _wsqr ), colOffset is the "x-position" of the block,
  // and rowOffset the "y-position" in the matrix.
  void subMakeRowEntries ( vector<typename aol::Row<DataType>::RowEntry > &vec, const int localRowNum,
                           const int colOffset, const int rowOffset, const int bundleNr,
                           int &elementNum ) const {
    int startTripDiag   = 0;      // start index of the trip-diagonal ( -1 or 0 )
    int endTripDiag     = 2;      // end index (0 or 1)
    int startLocalIndex = 0;      // start index in one trip-diagonal (-1 or 0)
    int endLocalIndex   = 2;      // end index in one trip-diagonal (0 or 1)


    if ( localRowNum < _w )         startTripDiag = 1;    // only 2 trip-diagonals
    if ( localRowNum >= ( _w - 1 ) *_w ) endTripDiag   = 1;    // only 2 trip-diagonals

    if ( localRowNum % _w == 0 )    startLocalIndex = 1;  // only two elements != 0 in the trip-diag.
    if ( localRowNum % _w == _w - 1 ) endLocalIndex = 1;  // only two elements != 0 in the trip-diag.

    for ( int tripDiag = startTripDiag; tripDiag <= endTripDiag; ++tripDiag ) {
      for ( int localIndex = startLocalIndex; localIndex <= endLocalIndex; ++localIndex ) {
        vec[elementNum].col = colOffset + localRowNum - ( _w + 1 ) + ( tripDiag * _w ) + localIndex;
        vec[elementNum].value =
          _rows[3 * bundleNr + tripDiag][rowOffset + localRowNum][localIndex];

        ++elementNum;
      }
    }
  }

  // method gets a vector and the number of a row and fills this vector with
  // row-entries, that means with the elements which aren't 0 and their column number.
  void makeRowEntries ( vector<typename aol::Row<DataType>::RowEntry > &vec, const int RowNum ) const {
    int blockNum    = RowNum / _wsqr;
    int localRowNum = RowNum % _wsqr;
    int rowOffset   = blockNum * _wsqr;

    int elementNum = 0;           // counter for the elements != 0
    vec.resize ( 27 );            // maximal elements != 0


    // three cases are to differentiate:
    if ( blockNum == 0 ) {
      subMakeRowEntries ( vec, localRowNum, 0, 0, 1, elementNum );
      subMakeRowEntries ( vec, localRowNum, _wsqr, 0, 2, elementNum );
    }
    else if ( blockNum > 0 && blockNum < _w - 1 ) {    // the inner case => three blocks
      subMakeRowEntries ( vec, localRowNum, ( blockNum - 1 ) * _wsqr, rowOffset , 0, elementNum );
      subMakeRowEntries ( vec, localRowNum, blockNum * _wsqr, rowOffset , 1, elementNum );
      subMakeRowEntries ( vec, localRowNum, ( blockNum + 1 ) * _wsqr, rowOffset , 2, elementNum );
    }
    else if ( blockNum == _w - 1 ) {
      subMakeRowEntries ( vec, localRowNum, ( _w - 2 ) * _wsqr, rowOffset, 0, elementNum );
      subMakeRowEntries ( vec, localRowNum, ( _w - 1 ) * _wsqr, rowOffset, 1, elementNum );
    }

    // now resize the vector to the number of elements which were REALLY used
    vec.resize ( elementNum );
  }

  // -------------------------- END MAKEROWENTRIES ----------------------------------------


  FastUniformGridMatrix<DataType, qc::QC_3D, BaseClass>& addMultiple ( const FastUniformGridMatrix<DataType, qc::QC_3D, BaseClass> &other, DataType factor ) {
    if ( _size != other._size ) {
      throw aol::Exception ( "Sizes of matrices do not match", __FILE__, __LINE__ );
    }
    for ( int b = 0; b < 9; ++b ) {
      for ( int i = 0; i < _size; ++i ) {
        _rows[b][i] += factor * other._rows[b][i];
      }
    }
    return *this;
  }


  // WARNING: this does not make any compatibility checks.
  // the argument must have a compatible matrix structure, otherwise this function won't work.
  FastUniformGridMatrix<DataType, qc::QC_3D, BaseClass>& operator+= ( const aol::Matrix<DataType> &Matrix ) {
    vector<typename aol::Row<DataType>::RowEntry > vec;
    for ( int i = 0; i < _size; ++i ) {
      Matrix.makeRowEntries ( vec, i );
      for ( typename vector<typename aol::Row<DataType>::RowEntry >::iterator it = vec.begin();
            it != vec.end(); ++it ) {
        add ( i, it->col, it->value );
      }
    }
    return *this;
  }

  FastUniformGridMatrix<DataType, qc::QC_3D, BaseClass>& operator+= ( const FastUniformGridMatrix<DataType, qc::QC_3D, BaseClass> &other ) {
    if ( _size != other._size ) {
      throw aol::Exception ( "Sizes of matrices do not match", __FILE__, __LINE__ );
    }
    for ( int b = 0; b < 9; ++b ) {
      for ( int i = 0; i < _size; ++i ) {
        _rows[b][i] += other._rows[b][i];
      }
    }
    return *this;
  }

  FastUniformGridMatrix<DataType, qc::QC_3D, BaseClass>& operator= ( const FastUniformGridMatrix<DataType, qc::QC_3D, BaseClass> &other ) {
    if ( _size != other._size ) {
      throw aol::Exception ( "Sizes of matrices do not match", __FILE__, __LINE__ );
    }
    for ( int b = 0; b < 9; ++b ) {
      for ( int i = 0; i < _size; ++i ) {
        _rows[b][i] = other._rows[b][i];
      }
    }
    return *this;
  }

  FastUniformGridMatrix<DataType, qc::QC_3D, BaseClass>& operator*= ( const DataType factor ) {
    for ( int k = 0; k < 9; ++k ) {
      for ( int i = 0; i < _size; ++i )  {
        _rows[k][i] *= factor;
      }
    }
    return *this;
  }

protected:

  // function gets an offset from the diagonal element and returns
  // the number of the row-block (trip-diag) where the belonging element is inside, and
  // returns the local index in this block!
  inline void getBlockAndIndex ( int offset, int &block, int &index ) const {
    // 9 possibilities for the block
    // a) the plane in which our element lives (that are the inner 3 diagonals)
    if ( offset >= -_w - 1 && offset <= _w + 1 ) {
      if ( offset > - _w - 2 && offset < - _w + 2 ) {
        block = 3;
        index = offset + _w + 1;
      }
      else if ( offset > -2 && offset < 2 ) {
        block = 4;
        index = offset + 1;
      }
      else if ( offset > _w - 2 && offset < _w + 2 ) {
        block = 5;
        index = offset - _w + 1;
      }
      else {
        block = index = -1;
      }
    } else {
      // b) the plane in front of  the 3d-element (the 3 diags on the left)
      if ( offset < -_w - 1 ) {
        if ( offset > - _wsqr - _w - 2 && offset < - _wsqr - _w + 2 ) {
          block = 0;
          index = offset + _wsqr + _w + 1;
        }
        else if ( offset > - _wsqr - 2 && offset < - _wsqr + 2 ) {
          block = 1;
          index = offset + _wsqr + 1;
        }
        else if ( offset > - _wsqr + _w - 2 && offset < - _wsqr + _w + 2 ) {
          block = 2;
          index = offset + _wsqr - _w + 1;
        }
        else {
          block = index = -1;
        }
      } else {
        // c) the plane behind the 3d-element (the 3 diags on the right)
        if ( offset > _w + 1 ) {
          if ( offset > _wsqr - _w - 2 && offset < _wsqr - _w + 2 ) {
            block = 6;
            index = offset - _wsqr + _w + 1;
          }
          else if ( offset > _wsqr - 2 && offset < _wsqr + 2 ) {
            block = 7;
            index = offset - _wsqr + 1;
          }
          else if ( offset > _wsqr + _w - 2 && offset < _wsqr + _w + 2 ) {
            block = 8;
            index = offset - _wsqr - _w + 1;
          }
          else {
            block = index = -1;
          }
        }
      }
    }

  }

  // the vector containing the non-zero elements
  // used to be vector<aol::Vec<3, DataType> > _rows[9];
  aol::Vec<3, DataType> * _rows[9];
};

template <typename T, qc::Dimension d>
ostream &operator<< ( ostream &os, const FastUniformGridMatrix<T,d> &mat ) {
  int n = mat.getNumRows();
  vector< typename aol::Row<T>::RowEntry > rowEntries;
  os << setprecision(2);
  for (int i = 0; i < n; ++i)
  {
    int col = 0;
    mat.makeRowEntries(rowEntries, i);
    for ( typename vector< typename aol::Row<T>::RowEntry >::iterator it = rowEntries.begin(); it != rowEntries.end(); ++it ) {
      for(; col < it->col; ++col)
        os << "0" << "\t";
      os << it->value << "\t";
      ++col;
    }
    for(; col < mat.getNumCols(); ++col)
      os << "0" << "\t";
    os << "\n";
  }
  return os;
}


} // end namespace

#endif
