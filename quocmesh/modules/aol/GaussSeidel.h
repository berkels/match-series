#ifndef __GAUSSSEIDEL_H
#define __GAUSSSEIDEL_H

#include <op.h>
#include <vec.h>
#include <multiVector.h>
#include <smallMat.h>
#include <rows.h>
#include <progressBar.h>

namespace aol {
//! Modes for Gauss-Seidel row iterations: symmetric alternates between forward and backward, red-black is even-only and odd-only, zebra2 is first quarter, third quarter, second quarter, fourth quarter (for parallelization)
enum GaussSeidelSweepingMode { GAUSS_SEIDEL_FORWARD,
                               GAUSS_SEIDEL_SYMMETRIC,
                               GAUSS_SEIDEL_RED_BLACK,
                               GAUSS_SEIDEL_BACKWARD,
                               GAUSS_SEIDEL_EVEN_ONLY,
                               GAUSS_SEIDEL_ODD_ONLY,
                               GAUSS_SEIDEL_ZEBRA2,
                               GAUSS_SEIDEL_ZEBRA2_SYMMETRIC,
                               GAUSS_SEIDEL_ZEBRA2_FORWARD,
                               GAUSS_SEIDEL_ZEBRA2_BACKWARD };

}

namespace {
  // nameless namespace that is only visible in this header file

  template < typename RealType, class VectorType, class OpType >
  void performGaussSeidelSweep ( const VectorType & , VectorType & , const aol::GaussSeidelSweepingMode, const OpType & , const RealType );

  template < typename RealType, /*aol::Vector<RealType>,*/ class OpType >
  void performGaussSeidelSweep ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest, const aol::GaussSeidelSweepingMode gss, const OpType &op, const RealType relax ) {
    std::vector<typename aol::Row<RealType>::RowEntry > vec;
    for ( int row = 0; row < static_cast<int> ( Arg.size() ); ++row ) {
      int row_index = row;
      if ( gss == aol::GAUSS_SEIDEL_BACKWARD ){
        row_index =  static_cast<int> ( Arg.size() ) - 1 - row ;
      }

      if ( ( gss == aol::GAUSS_SEIDEL_EVEN_ONLY ) && ( row_index % 2 == 1 ) ){
        continue;
      }

      if ( ( gss == aol::GAUSS_SEIDEL_ODD_ONLY  ) && ( row_index % 2 == 0 ) ){
        continue;
      }

      op.makeRowEntries ( vec, row_index );

      RealType
        diag = 0.0,
        v = 0.0;

      for ( typename std::vector<typename aol::Row<RealType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
        if ( it->col == row_index )
          diag = it->value;

        v += it->value * Dest[ it->col ];

      }

      Dest[ row_index ] += relax * ( Arg[ row_index ] - v ) / diag;

    }
  }

  template < typename RealType, /*aol::MultiVector<RealType>,*/ class OpType >
  void performGaussSeidelSweep ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest, const aol::GaussSeidelSweepingMode gss, const OpType &op, const RealType relax ) {
    switch ( gss ) {
    case aol::GAUSS_SEIDEL_FORWARD:

      for ( int comp_row = 0; comp_row < op.getNumRows(); ++comp_row ) {
        for ( int row = 0; row < Arg[comp_row].size(); ++row ) {
          int row_index = row;

          RealType
            diag = 0.0,
            v = 0.0;

          for ( int comp_col = 0; comp_col < op.getNumCols(); ++comp_col ) {
            std::vector<typename aol::Row<RealType>::RowEntry > vec;
            op.getReference(comp_row, comp_col).makeRowEntries ( vec, row_index );

            for ( typename vector<typename aol::Row<RealType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
              if ( ( comp_row == comp_col ) && ( it->col == row_index ) ) {
                diag = it->value;
              }

              v += it->value * Dest[comp_col][ it->col ];
            }
          }

          Dest[comp_row][row_index] += relax * ( Arg[comp_row][row_index] - v ) / diag;

        }
      }

      break;
    case aol::GAUSS_SEIDEL_BACKWARD:

      for ( int comp_row = op.getNumRows() - 1; comp_row >= 0 ; --comp_row ) {
        for ( int row = Arg[comp_row].size() - 1; row >= 0 ; --row ) {
          int row_index = row;

          RealType
            diag = 0.0,
            v = 0.0;

          for ( int comp_col = 0; comp_col < op.getNumCols(); ++comp_col ) {
            std::vector<typename aol::Row<RealType>::RowEntry > vec;
            op.getReference(comp_row, comp_col).makeRowEntries ( vec, row_index );

            for ( typename vector<typename aol::Row<RealType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
              if ( ( comp_row == comp_col ) && ( it->col == row_index ) ) {
                diag = it->value;
              }

              v += it->value * Dest[comp_col][ it->col ];
            }
          }

          Dest[comp_row][row_index] += relax * ( Arg[comp_row][row_index] - v ) / diag;

        }
      }

      break;
    default:
      throw aol::UnimplementedCodeException ( "aol::GaussSeidelSweeper, MultiVector case: GaussSeidel sweeping mode is not supported", __FILE__, __LINE__ );
    }
  }


  template < typename RealType, class VectorType, class OpType >
  void performSelectiveGaussSeidelSweep ( const VectorType & , VectorType & , const aol::GaussSeidelSweepingMode, const aol::BitVector &, const OpType & , const RealType );

  template < typename RealType, /*aol::Vector<RealType>,*/ class OpType >
  void performSelectiveGaussSeidelSweep ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest, const aol::GaussSeidelSweepingMode gss, const aol::BitVector &sweepMask, const OpType &op, const RealType relax ) {
    std::vector<typename aol::Row<RealType>::RowEntry > vec;
    for ( int row = 0; row < static_cast<int> ( Arg.size() ); ++row ) {
      int row_index = row;
      if ( gss == aol::GAUSS_SEIDEL_BACKWARD ){
        row_index =  static_cast<int> ( Arg.size() ) - 1 - row ;
      }

      if ( ( gss == aol::GAUSS_SEIDEL_EVEN_ONLY ) && ( row_index % 2 == 1 ) ){
        continue;
      }

      if ( ( gss == aol::GAUSS_SEIDEL_ODD_ONLY  ) && ( row_index % 2 == 0 ) ){
        continue;
      }

      if ( sweepMask.get( row_index ) ) {

        op.makeRowEntries ( vec, row_index );

        RealType
          diag = 0.0,
          v = 0.0;

        for ( typename std::vector<typename aol::Row<RealType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
          if ( it->col == row_index )
            diag = it->value;

          v += it->value * Dest[ it->col ];

        }

        Dest[ row_index ] += relax * ( Arg[ row_index ] - v ) / diag;

      }
    }
  }

  template < typename RealType, /*aol::MultiVector<RealType>,*/ class OpType >
  void performSelectiveGaussSeidelSweep ( const aol::MultiVector<RealType>&, aol::MultiVector<RealType>&, const aol::GaussSeidelSweepingMode, const OpType&, const RealType ) {
    throw aol::Exception( "performSelectiveGaussSeidelSweeo< RealType, aol::MultiVector<RealType>, OpType > not implemented", __FILE__, __LINE__ );
  }

}

namespace aol {

/** GaussSeidelSweeper's apply performs one Gauss-Seidel smoothing step.
 *  Such sweepers are used in smoothers (performing a fixed number of iterations) and solvers (performing iterations up to some stopping criterion)
 *  This parent class is useless but necessary as a common base for scalar and vector-valued case.
 *  \author Schwen (based on older code)
 */
template< typename RealType, typename VectorType, typename OpType >
class GaussSeidelSweeper {
protected:
  const OpType&  _op;
  RealType       _relax;

public:
  GaussSeidelSweeper ( const OpType &Op, const RealType Relax ) : _op ( Op ), _relax ( Relax ) {
  }

  void setRelax ( const RealType Relax ) {
    _relax = Relax;
  }

  void apply ( const VectorType &Arg, VectorType &Dest, const aol::GaussSeidelSweepingMode gss ) const {
    performGaussSeidelSweep( Arg, Dest, gss, _op, _relax );
  }

private:
  GaussSeidelSweeper ( ) {
    throw aol::UnimplementedCodeException ( "aol::GaussSeidelSweeper standard constructor not implemented", __FILE__, __LINE__ );
  }

  GaussSeidelSweeper ( const GaussSeidelSweeper<RealType, VectorType, OpType>& ){
    throw aol::UnimplementedCodeException ( "aol::GaussSeidelSweeper copy constructor not implemented", __FILE__, __LINE__ );
  }

  GaussSeidelSweeper<RealType, VectorType, OpType>& operator= ( const GaussSeidelSweeper<RealType, VectorType, OpType>& ) {
    throw aol::UnimplementedCodeException ( "aol::GaussSeidelSweeper::operator= not implemented", __FILE__, __LINE__ );
  }

}; // end class GaussSeidelSweeper


/** Gauss-Seidel sweeper that only operates on subset of vector entries specified by a bitField.
 *  \author Schwen
 */
template< typename RealType, typename VectorType, typename OpType >
class GaussSeidelSelectiveSweeper {
protected:
  const OpType&        _op;
  const aol::BitVector& _sweepMask;
  RealType             _relax;

public:
  GaussSeidelSelectiveSweeper ( const OpType &Op, const aol::BitVector &SweepMask, const RealType Relax ) : _op ( Op ), _sweepMask ( SweepMask ), _relax ( Relax ) {
  }

  void setRelax ( const RealType Relax ) {
    _relax = Relax;
  }

  void apply ( const VectorType &Arg, VectorType &Dest, const aol::GaussSeidelSweepingMode gss ) const {
    performSelectiveGaussSeidelSweep( Arg, Dest, gss, _sweepMask, _op, _relax );
  }

  GaussSeidelSelectiveSweeper ( ) : _op ( *reinterpret_cast< const OpType* >( NULL ) ), _sweepMask ( *reinterpret_cast< const aol::BitVector* > ( NULL ) ) {
    throw aol::UnimplementedCodeException ( "aol::GaussSeidelSelectiveSweeper standard constructor not implemented", __FILE__, __LINE__ );
  }

private:
  GaussSeidelSelectiveSweeper ( const GaussSeidelSelectiveSweeper<RealType, VectorType, OpType>& );

  GaussSeidelSelectiveSweeper<RealType, VectorType, OpType>& operator= ( const GaussSeidelSelectiveSweeper<RealType, VectorType, OpType>& ) {
    throw aol::UnimplementedCodeException ( "aol::GaussSeidelSelectiveSweeper::operator= not implemented", __FILE__, __LINE__ );
  }

}; // end class GaussSeidelSelectiveSweeper



/** The BlockGaussSeidelSweeper performs 3x3-block-wise Gauss-Seidel iterations.
 *  It is useful for vector-valued computations where 3 components are
 *  computed for each geometric node and such iterations can be viewed
 *  as an implicit reordering of the unknowns so that all components
 *  for one node are neighbors in the value vector. Then three
 *  components are treated simultaneously in the Gauss-Seidel
 *  iterations
 *  \author Schwen (based on older code)
 */

template< typename RealType, typename BlockOpType >
class BlockGaussSeidelSweeper {

protected:
  typedef aol::MultiVector<RealType>   VectorType;

  const BlockOpType&  _blockOp;
  RealType            _relax;
  std::vector< aol::Matrix33<RealType> > _inv_diags; // for cacheing 3x3 matrices and their inverses

public:
  BlockGaussSeidelSweeper ( const BlockOpType &BlockOp,
                            const RealType Relax )
    : _blockOp ( BlockOp ),
      _relax ( Relax ),
      _inv_diags ( 0 ) {

    if ( this->_blockOp.getNumCols() != 3 || this->_blockOp.getNumRows() != 3 ) {
      throw aol::Exception ( "BlockGaussSeidelSweeper not implemented for non-3D.", __FILE__, __LINE__ );
    }

    cacheDiagInverses ( this->_blockOp.getReference(0,0).getNumRows() );
  }

  virtual ~BlockGaussSeidelSweeper ( ) {
  }

protected:
  void cacheDiagInverses ( const int n ) {

    _inv_diags.resize( n );

#ifdef VERBOSE
    aol::ProgressBar<> pb( "Cacheing 3x3 inverses" );
    pb.start(n);
#endif
    for ( int row = 0; row < n; ++row ) {
#ifdef VERBOSE
      pb++;
#endif

      aol::Matrix33<RealType> diag;
      std::vector< typename aol::Row<RealType>::RowEntry > vec;

      for ( short a = 0; a < 3; ++a ) {
        for ( short b = 0; b < 3; ++b ) {
          if ( this->_blockOp.getPointer ( a, b ) != NULL ) {
            this->_blockOp.getReference ( a, b ).makeRowEntries ( vec, row );

            for ( typename std::vector<typename aol::Row<RealType>::RowEntry >::iterator it = vec.begin();  it != vec.end(); ++it ) {
              if ( it->col == row )
                diag[a][b] = it->value;
            }
          }
        }
      }

#ifdef VERBOSE
      if ( ( ( fabs ( diag.det() ) > 1.0e10 ) || fabs ( diag.det() < 1.0e-5 ) ) && ( diag.det() != 0.0 ) ) {
        cout << "DIAGDET " << diag.det() << endl << diag << endl;
        aol::Matrix33<RealType> test;
        test = diag.inverse();
        cout << test << endl;
        test *= diag;
        cout << test << endl;
      }
      // todo: think about what to do in case of poorly inverted matrices or how to avoid them
#endif

      if ( diag.det() != 0.0 )
        _inv_diags[row] = diag.inverse();
      else
        _inv_diags[row].setIdentity();

    }
#ifdef VERBOSE
    cerr << endl;
#endif
  }

public:
  void setRelax ( const RealType Relax ) {
    _relax = Relax;
  }

  void apply ( const VectorType &Arg, VectorType &Dest, const aol::GaussSeidelSweepingMode gss ) const {

#ifdef BOUNDS_CHECK
    if ( Dest[0].size() != Dest[1].size() || Dest[0].size() != Dest[2].size() || Arg[0].size() != Dest[0].size() ||  Arg[1].size() != Dest[1].size() ||  Arg[2].size() != Dest[2].size() ) {
      throw aol::Exception ( "Block GaussSeidelSweeper: incompatible MultiVectors.", __FILE__, __LINE__ );
    }
#endif

    const int n = static_cast<int> ( Dest[0].size() );

    switch ( gss ) {
    case GAUSS_SEIDEL_FORWARD:
      for ( int row = 0; row < n; ++row )
        treatRow ( Arg, Dest, row );
      break;

    case GAUSS_SEIDEL_BACKWARD:
      for ( int row = n - 1; row >= 0; --row )
        treatRow ( Arg, Dest, row );
      break;

    case GAUSS_SEIDEL_EVEN_ONLY:
      for ( int row = 0; row < n; row += 2 )
        treatRow ( Arg, Dest, row );
      break;

    case GAUSS_SEIDEL_ODD_ONLY:
      for ( int row = 1; row < n; row += 2 )
        treatRow ( Arg, Dest, row );
      break;

    case GAUSS_SEIDEL_ZEBRA2_FORWARD:
      {

        aol::Vector<int> ubnds(5);
        for ( int k = 0; k < 5; ++k ) {
          ubnds[k] = ( k * n ) / 4;
        }

#ifdef _OPENMP

#pragma omp parallel sections
        {
#pragma omp section
          {
            for ( int row = ubnds[0]; row < ubnds[1]; ++row ) {
              treatRow ( Arg, Dest, row );
            }
          }
#pragma omp section
          {
            for ( int row = ubnds[2]; row < ubnds[3]; ++row ) {
              treatRow ( Arg, Dest, row );
            }
          }
        }

#pragma omp parallel sections
        {
#pragma omp section
          {
            for ( int row = ubnds[1]; row < ubnds[2]; ++row ) {
              treatRow ( Arg, Dest, row );
            }
          }
#pragma omp section
          {
            for ( int row = ubnds[3]; row < ubnds[4]; ++row ) {
              treatRow ( Arg, Dest, row );
            }
          }
        }

#else

        // zebra2 mode does not really make sense if we do not parallelize, but anyway ...
        for ( int row = ubnds[0]; row < ubnds[1]; ++row )
          treatRow ( Arg, Dest, row );
        for ( int row = ubnds[2]; row < ubnds[3]; ++row )
          treatRow ( Arg, Dest, row );
        for ( int row = ubnds[1]; row < ubnds[2]; ++row )
          treatRow ( Arg, Dest, row );
        for ( int row = ubnds[3]; row < ubnds[4]; ++row )
          treatRow ( Arg, Dest, row );

#endif

      }
      break;

    case GAUSS_SEIDEL_ZEBRA2_BACKWARD:
      {

        aol::Vector<int> ubnds(5);
        for ( int k = 0; k < 5; ++k ) {
          ubnds[k] = ( k * n ) / 4; // integer division
        }

#ifdef _OPENMP

#pragma omp parallel sections
        {
#pragma omp section
          {
            for ( int row = ubnds[4] - 1; row >= ubnds[3]; --row ) {
              treatRow ( Arg, Dest, row );
            }
          }
#pragma omp section
          {
            for ( int row = ubnds[2] - 1; row >= ubnds[1]; --row ) {
              treatRow ( Arg, Dest, row );
            }
          }
        }

#pragma omp parallel sections
        {
#pragma omp section
          {
            for ( int row = ubnds[3] - 1; row >= ubnds[2]; --row ) {
              treatRow ( Arg, Dest, row );
            }
          }
#pragma omp section
          {
            for ( int row = ubnds[1] - 1; row >= ubnds[0]; --row ) {
              treatRow ( Arg, Dest, row );
            }
          }
        }

#else

        // zebra2 mode does not really make sense if we do not parallelize, but anyway ...
        for ( int row = ubnds[4] - 1; row >= ubnds[3]; --row )
          treatRow ( Arg, Dest, row );
        for ( int row = ubnds[2] - 1; row >= ubnds[1]; --row )
          treatRow ( Arg, Dest, row );
        for ( int row = ubnds[3] - 1; row >= ubnds[2]; --row )
          treatRow ( Arg, Dest, row );
        for ( int row = ubnds[1] - 1; row >= ubnds[0]; --row )
          treatRow ( Arg, Dest, row );

#endif

      }
      break;

    default:
      throw aol::Exception ( "BlockGaussSeidelSweeper: illegal GaussSeidelSweepingMode", __FILE__, __LINE__ );
    }

  }

  virtual void treatRow ( const VectorType &Arg, VectorType &Dest, const int row_index ) const {
    doTreatRow ( Arg, Dest, row_index );
  }

  void doTreatRow ( const VectorType &Arg, VectorType &Dest, const int row_index ) const {
    aol::Vec3<RealType> v;

    for ( int a = 0; a < 3; ++a ) {
      for ( int b = 0; b < 3; ++b ) {
        if ( this->_blockOp.getPointer ( a, b ) != NULL ) {
          std::vector< typename aol::Row<RealType>::RowEntry > vec;
          this->_blockOp.getReference ( a, b ).makeRowEntries ( vec, row_index );

          for ( typename std::vector<typename aol::Row<RealType>::RowEntry >::iterator it = vec.begin();  it != vec.end(); ++it ) {
            v[ a ] += it->value * Dest[ b ][ it->col ];
          }  // columns done

        }
      }
    }

    aol::Vec3<RealType> dummy, tmp;

    for ( int a = 0; a < 3; ++a ) {
      tmp[a] = Arg[a][row_index];
    }

    tmp -= v;
    tmp *= this->_relax;

    _inv_diags[row_index].mult ( tmp, dummy );

    for ( int a = 0; a < 3; ++a ) {
      Dest[a][row_index] += dummy[a];
    }
  }

private:
  BlockGaussSeidelSweeper ( ) {
    throw aol::UnimplementedCodeException ( "aol::BlockGaussSeidelSweeper standard constructor not implemented", __FILE__, __LINE__ );
  }

  BlockGaussSeidelSweeper ( const BlockGaussSeidelSweeper< RealType, BlockOpType >& );

  BlockGaussSeidelSweeper< RealType, BlockOpType >& operator= ( const BlockGaussSeidelSweeper< RealType, BlockOpType >& ) {
    throw aol::UnimplementedCodeException ( "aol::BlockGaussSeidelSweeper::operator= not implemented", __FILE__, __LINE__ );
  }

}; // end class BlockGaussSeidelSweeper


/** Block Gauss-Seidel sweeper that only operates on subset of vector entries specified by a bitField.
 *  \author Schwen
 */
template< typename RealType, typename BlockOpType >
class BlockGaussSeidelSelectiveSweeper : public aol::BlockGaussSeidelSweeper< RealType, BlockOpType > {
protected:
  typedef aol::MultiVector<RealType> VectorType;

  const aol::BitVector& _sweepMask;

public:
  BlockGaussSeidelSelectiveSweeper ( const BlockOpType &BlockOp, const RealType Relax  )
    : BlockGaussSeidelSweeper< RealType, BlockOpType > ( BlockOp, Relax ),
      _sweepMask ( *reinterpret_cast< const aol::BitVector* >( NULL ) ) {
    throw aol::Exception ( "aol::BlockGaussSeidelSelectiveSweeper: constructor without SweepMask exists only for compatibility and must not be called.", __FILE__, __LINE__ );
  }

  BlockGaussSeidelSelectiveSweeper ( const BlockOpType &BlockOp,
                                     const RealType Relax,
                                     const aol::BitVector &SweepMask )
    : BlockGaussSeidelSweeper< RealType, BlockOpType > ( BlockOp, Relax ),
    _sweepMask ( SweepMask ){
  }

  virtual void treatRow ( const VectorType &Arg, VectorType &Dest, const int row_index ) const {
    if ( _sweepMask.get ( row_index ) ) {
      this->doTreatRow ( Arg, Dest, row_index );
    }
  }
};

}

#endif
