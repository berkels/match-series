#include <fastUniformGridMatrix.h>
#include <cellCenteredGrid.h>

#ifdef _OPENMP
WARNING_OFF ( uninitialized )
#endif

template <class _DataType, typename BaseClass>
qc::FastUniformGridMatrix<_DataType, qc::QC_2D, BaseClass>::FastUniformGridMatrix ( const qc::CellCenteredCubicGrid<qc::QC_2D> &Grid )
: BaseClass(Grid.getNumberOfNodes(), Grid.getNumberOfNodes()), _w ( Grid.getNumXYZ() ), _size ( Grid.getNumberOfNodes() ) {
  if ( Grid.getGridDepth() == 0 ) {
    throw aol::Exception ( "FastUniformGridMatrix doesn't work for level 0.\n", __FILE__, __LINE__ );
  }
  for ( int i = 0; i < 3; ++i )
    _rows[i] = new aol::Vec<3, DataType> [ _size ];
  setZero();  // be tidy when born
}

template <class _DataType, typename BaseClass>
void qc::FastUniformGridMatrix<_DataType, qc::QC_2D, BaseClass>::applyAdd ( const aol::Vector<_DataType> &Arg, aol::Vector<_DataType> &Dest ) const {
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for ( int i = 0; i <= _w + 1; ++i ) {
    multCarefullyAtBoundary ( i, Arg, Dest );
  }
#if 1
  const DataType * ArgPtr  = Arg.getData();
  DataType * DestPtr = Dest.getData();
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for ( int i = _w + 2; i < _size - _w - 2; ++i ) {
    for ( int k = 0; k < 3; ++k ) {
      int globos = _w * ( k - 1 ) - 1;
      int g = globos + i;
      for ( int j = 0; j < 3; ++j ) {
        DestPtr[i] += ArgPtr[g++] * _rows[k][i][j];
      }
    }
  }
#endif
#if 0
  const int seg_len = 10;
  for ( int seg = _w + 2; seg < _size - _w - 2; seg += seg_len ) {
    const int seg_end = aol::Min ( seg + seg_len, _size - _w - 2 );
    for ( int k = 0; k < 3; ++k ) {
      int globos = _w * ( k - 1 ) - 1;
      for ( int i = seg; i < seg_end; ++i ) {
        // cerr << "inner row = " << i << endl;
        int g = globos + i;
        for ( int j = 0; j < 3; ++j ) {
          Dest[i] += Arg[g++] * _rows[k][i][j];
        }
      }
    }
  }
#endif
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for ( int i = aol::Max ( _w + 2, _size - _w - 2 ); i < _size; ++i ) {
    multCarefullyAtBoundary ( i, Arg, Dest );
  }
}


template <class _DataType, typename BaseClass>
void qc::FastUniformGridMatrix<_DataType, qc::QC_3D, BaseClass>::applyAdd ( const aol::Vector<_DataType> &Arg, aol::Vector<_DataType> &Dest ) const {
  // the first _wsqr elements are a special case, with only two sub-applies
  subApplyAdd ( 0, 0, 1, Arg, Dest );
  subApplyAdd ( 0, _wsqr, 2, Arg, Dest );

  // now the inner elements (_w-2 applies with 3 sub-applies)
  // subApplyAdd ( i*_wsqr, ... ) writes to the Dest components (i*_wsqr) till (i*_wsqr +w_sqr-1).
  // So the memory write access of calls to subApplyAdd with different i values doesn't overlap,
  // i.e. we may parallelize.
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for ( int i = 1; i < _w - 1; ++i )
    for ( int j = 0; j < 3; ++j )
      subApplyAdd ( i*_wsqr, ( i + j - 1 ) *_wsqr, j, Arg, Dest );

  // the last _wsqr elements are a special case too, with only two sub-applies
  subApplyAdd ( ( _w - 1 ) *_wsqr, ( _w - 2 ) *_wsqr, 0, Arg, Dest );
  subApplyAdd ( ( _w - 1 ) *_wsqr, ( _w - 1 ) *_wsqr, 1, Arg, Dest );
}

template class qc::FastUniformGridMatrix<float, qc::QC_2D, aol::GenSparseOp<float> >;
template class qc::FastUniformGridMatrix<double, qc::QC_2D, aol::GenSparseOp<double> >;
template class qc::FastUniformGridMatrix<long double, qc::QC_2D, aol::GenSparseOp<long double> >;

template class qc::FastUniformGridMatrix<float, qc::QC_2D, aol::MatrixAbstractBase<float> >;
template class qc::FastUniformGridMatrix<double, qc::QC_2D, aol::MatrixAbstractBase<double> >;
template class qc::FastUniformGridMatrix<long double, qc::QC_2D, aol::MatrixAbstractBase<long double> >;

template class qc::FastUniformGridMatrix<float, qc::QC_3D, aol::GenSparseOp<float> >;
template class qc::FastUniformGridMatrix<double, qc::QC_3D, aol::GenSparseOp<double> >;
template class qc::FastUniformGridMatrix<long double, qc::QC_3D, aol::GenSparseOp<long double> >;

template class qc::FastUniformGridMatrix<float, qc::QC_3D, aol::MatrixAbstractBase<float> >;
template class qc::FastUniformGridMatrix<double, qc::QC_3D, aol::MatrixAbstractBase<double> >;
template class qc::FastUniformGridMatrix<long double, qc::QC_3D, aol::MatrixAbstractBase<long double> >;
