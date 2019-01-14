#ifndef __PRECONDITIONER_H
#define __PRECONDITIONER_H

#include <solver.h>

namespace aol {

/** Basis class for template specialization */
template< typename VectorType >
class DiagonalPreconditioner;


/**
 * \brief A class for diagonally preconditioning matrices (and matrices only)
 * \author Schwen
 * \ingroup preconditioner
 */
template< typename RealType >
class DiagonalPreconditioner < aol::Vector<RealType> > : public Op< aol::Vector<RealType> > {
protected:
  aol::DiagonalMatrix< RealType > _diagInv;

public:
  template< typename MatrixType >
  explicit DiagonalPreconditioner ( const MatrixType &Matrix ) : _diagInv ( Matrix.getNumRows() ) {
    _diagInv.setToAbsInverseDiagonalOf ( Matrix );
  }

  virtual void applyAdd ( const Vector<RealType> &Arg, Vector<RealType> &Dest ) const {
    _diagInv.applyAdd ( Arg, Dest );
  }

  virtual void apply ( const Vector<RealType> &Arg, Vector<RealType> &Dest ) const {
    _diagInv.apply ( Arg, Dest );
  }

  void getPreconditionerMatrix ( aol::DiagonalMatrix<RealType> &mat ) const {
    mat = _diagInv;
  }

  const aol::DiagonalMatrix<RealType> &getPreconditionerMatrixReference ( ) const {
    return _diagInv;
  }

  template< typename MatrixType >
  void setForMatrix ( const MatrixType &Matrix ) {
    _diagInv.reallocate ( Matrix.getNumRows() );
    _diagInv.setToAbsInverseDiagonalOf ( Matrix );
  }
};


/**
 * \brief A class for diagonally preconditioning BlockMatrices
 * \author Teusner, Berkels
 * \ingroup preconditioner
 *
 * The diagonal entries of the BlockOp may contain NULL pointers.
 * In this case the corresponding block of the preconditioner is
 * set to the IdentityOp and getPreconditionerBlockMatrix can't be
 * called.
 */
template< typename RealType >
class DiagonalPreconditioner < aol::MultiVector<RealType> > : public Op< aol::MultiVector<RealType> > {
  aol::IdentityOp<aol::Vector<RealType> > _identity;
  std::vector< aol::DiagonalPreconditioner< aol::Vector<RealType> >* > _diag;

public:
  template< typename BlockMatrixType >
  explicit DiagonalPreconditioner ( const BlockMatrixType &BlockMat ) : _diag ( BlockMat.getNumRows() ) {
    if ( BlockMat.getNumRows() != BlockMat.getNumCols() ) {
      throw aol::Exception ( "DiagonalPreconditioner needs square block structure", __FILE__, __LINE__ );
    }
    for ( int i = 0; i < BlockMat.getNumRows(); ++i ) {
      if ( BlockMat.getPointer ( i, i ) != NULL ) {
        _diag[i] = new DiagonalPreconditioner< aol::Vector<RealType> > ( BlockMat.getReference ( i, i ) );
      } else {
        _diag[i] = NULL;
      }
    }
  }

  ~DiagonalPreconditioner () {
    for ( int i = 0; i < static_cast<int> ( _diag.size() ); ++i ) {
      delete ( _diag[i] );
      _diag[i] = NULL;
    }
  }

  virtual void applyAdd ( const MultiVector<RealType> &Arg, MultiVector<RealType> &Dest ) const {
    for ( int i = 0; i < static_cast<int> ( _diag.size() ); ++i ) {
      if ( _diag[i] ) {
        _diag[i]->applyAdd ( Arg[i], Dest[i] );
      } else {
        Dest[i] += Arg[i];
      }
    }
  }

  virtual void apply ( const MultiVector<RealType> &Arg, MultiVector<RealType> &Dest ) const {
    for ( int i = 0; i < static_cast<int> ( _diag.size() ); ++i ) {
      if ( _diag[i] ) {
        _diag[i]->apply ( Arg[i], Dest[i] );
      } else {
        Dest[i] = Arg[i];
      }
    }
  }

  void getPreconditionerBlockMatrix ( aol::BlockMatrix < aol::DiagonalMatrix< RealType > > &blockMat ) const {
    for ( int i = 0; i < static_cast<int> ( _diag.size() ); ++i ) {
      if ( _diag[i] ) {
        _diag[i]->getPreconditionerMatrix( blockMat.getReference(i, i) );
      } else {
        throw Exception ( "Entry that is not of type DiagonalPreconditioner<RealType> found (most likely an IdenityOp), can't get the preconditioner matrix for this.\n", __FILE__, __LINE__ );
      }
    }
  }

  template< typename BlockMatrixType >
  void setForBlockMatrix ( const BlockMatrixType &BlockMat ) {
    for ( int i = 0; i < static_cast<int> ( _diag.size() ); ++i ) {
      if ( _diag[i] ) {
        delete ( _diag[i] );
        _diag[i] = NULL;
      }
    }
    for ( int i = 0; i < BlockMat.getNumRows(); ++i ) {
      if ( BlockMat.getPointer ( i, i ) != NULL ) {
        _diag[i] = new DiagonalPreconditioner< aol::Vector<RealType> > ( BlockMat.getReference ( i, i ) );
      }
    }
  }

};


/**
 * \brief A class for diagonally preconditioning ABlockOps
 * \author Torben Paetz (JUB)
 * \ingroup preconditioner
 */
template< typename RealType >
class DiagonalPreconditioner < aol::AMultiVector<aol::MultiVector<RealType> > > : public Op< AMultiVector< MultiVector< RealType > > > {
  std::vector< aol::DiagonalPreconditioner< aol::MultiVector<RealType> >* > _diag;

public:
  template< typename BlockMatrixType >
  explicit DiagonalPreconditioner ( const BlockMatrixType &BlockMat ) : Op< AMultiVector< MultiVector< RealType > > > (){
    for ( int i = 0; i < BlockMat.getNumRows(); i++ ) {
       if ( BlockMat.getPointer ( i, i ) != NULL ) {
        _diag[i] = new DiagonalPreconditioner< aol::MultiVector<RealType> > ( BlockMat.getReference ( i, i ) );
      } else {
        _diag[i] = NULL;
      }
    }
  }

public:
  ~DiagonalPreconditioner () {
    for ( int i = 0; i < static_cast<int> ( _diag.size() ); ++i ) {
      delete ( _diag[i] );
      _diag[i] = NULL;
    }
  }

  void getPreconditionerBlockMatrix ( aol::BlockMatrix < aol::DiagonalMatrix< RealType > > &blockMat ) const {

    for ( int i = 0; i < static_cast<int> ( _diag.size() ); ++i ) {
      if ( _diag[i] ) {
        _diag[i]->getPreconditionerMatrix( blockMat.getReference(i, i) );
      } else {
        throw Exception ( "Entry that is not of type DiagonalPreconditioner<RealType> found (most likely an IdenityOp), can't get the preconditioner matrix for this.\n", __FILE__, __LINE__ );
      }
    }
  }

  virtual void applyAdd ( const AMultiVector<MultiVector<RealType> > &Arg, AMultiVector< MultiVector<RealType> > &Dest ) const {
    for ( int i = 0; i < static_cast<int> ( _diag.size() ); ++i ) {
      if ( _diag[i] ) {
        _diag[i]->applyAdd ( Arg[i], Dest[i] );
      } else {
        Dest[i] += Arg[i];
      }
    }
  }

  virtual void apply ( const AMultiVector<MultiVector<RealType> > &Arg, AMultiVector<MultiVector<RealType> > &Dest ) const {
    for ( int i = 0; i < static_cast<int> ( _diag.size() ); ++i ) {
      if ( _diag[i] ) {
        _diag[i]->apply ( Arg[i], Dest[i] );
      } else {
        Dest[i] = Arg[i];
      }
    }
  }


};



/** Basis class for template specialization */
template< typename VectorType >
class GeometricScalingPreconditioner;


/**
 * \brief A class for preconditioning matrices using geometric scaling (division by l2 norm of rows, cf. Gordon & Gordon, arXiv:0812.2769v2 [cs.MS])
 * \author Schwen
 * \ingroup preconditioner
 */
template< typename RealType >
class GeometricScalingPreconditioner < aol::Vector<RealType> > : public Op< aol::Vector<RealType> > {
protected:
  aol::DiagonalMatrix< RealType > _diagInv;

public:
  template< typename MatrixType >
  explicit GeometricScalingPreconditioner ( const MatrixType &Matrix ) : _diagInv ( Matrix.getNumRows() ) {
    for ( int i = 0; i < Matrix.getNumRows(); ++i ) {
      RealType l2sum = aol::NumberTrait<RealType>::zero;
      vector<typename Row<RealType>::RowEntry > vec;
      Matrix.makeRowEntries ( vec, i );
      for ( typename vector<typename Row<RealType>::RowEntry >::const_iterator it = vec.begin(); it != vec.end(); ++it ) {
        l2sum += aol::Sqr ( it->value );
      }

      // make sure we do not divide by zero
      if ( l2sum == aol::NumberTrait<RealType>::zero )
        l2sum = aol::NumberTrait<RealType>::one;

      _diagInv.set ( i, i, aol::NumberTrait<RealType>::one / sqrt ( l2sum ) );
    }
  }

  virtual void applyAdd ( const Vector<RealType> &Arg, Vector<RealType> &Dest ) const {
    _diagInv.applyAdd ( Arg, Dest );
  }

  virtual void apply ( const Vector<RealType> &Arg, Vector<RealType> &Dest ) const {
    _diagInv.apply ( Arg, Dest );
  }

  void getPreconditionerMatrix ( aol::DiagonalMatrix<RealType> &mat ) const {
    mat = _diagInv;
  }
};


/**
 * \brief A class for preconditioning block matrices using geometric scaling (division by l2 norm of unblocked rows)
 * \author Schwen
 * \ingroup preconditioner
 */
template< typename RealType >
class GeometricScalingPreconditioner < aol::MultiVector<RealType> > : public BlockOp< RealType, aol::Op< aol::Vector<RealType> > > {
public:
  template< typename BlockMatrixType >
  explicit GeometricScalingPreconditioner ( const BlockMatrixType &BlockMat ) : aol::BlockOp< RealType, aol::Op< aol::Vector<RealType> > > ( BlockMat.getNumRows(), BlockMat.getNumCols() ) {
    for ( int i = 0; i < BlockMat.getNumRows(); i++ ) {
      DiagonalMatrix<RealType> *pCurrentDiag = new DiagonalMatrix< RealType > ( BlockMat.getReference ( i, i ).getNumRows() );
      this->setPointer ( i, i, pCurrentDiag );

      for ( int ii = 0; ii < pCurrentDiag->getNumRows(); ++ii ) {
        RealType l2sum = aol::NumberTrait<RealType>::zero;
        vector<typename Row<RealType>::RowEntry > vec;
        BlockMat.makeUnblockedRowEntries ( vec, i, ii );

        for ( typename vector<typename Row<RealType>::RowEntry >::const_iterator it = vec.begin(); it != vec.end(); ++it ) {
          l2sum += aol::Sqr ( it->value );
        }

        // make sure we do not divide by zero
        if ( l2sum == aol::NumberTrait<RealType>::zero )
          l2sum = aol::NumberTrait<RealType>::one;

        pCurrentDiag->set ( ii, ii, aol::NumberTrait<RealType>::one / sqrt ( l2sum ) );
      }
    }
  }

public:
  ~GeometricScalingPreconditioner () {
    for ( int i = 0; i < this->getNumRows(); i++ ) {
      delete ( this->getPointer ( i, i ) );
    }
  }

  void getPreconditionerBlockMatrix ( aol::BlockMatrix < aol::DiagonalMatrix< RealType > > &blockMat ) const {
    for ( int i = 0; i < this->getNumRows(); i++ ) {
      blockMat.getReference ( i, i ) = dynamic_cast< const aol::DiagonalMatrix<RealType>& > ( this->getReference ( i, i ) );
    }
  }

};


/** Basis class for template specialization */
template< typename VectorType, typename MatrixType >
class SSORPreconditioner;


/** \brief A class for preconditioning sparse matrices by SSOR.
 *  Represents for a given matrix \f$A\f$ and its decomposition
 *  \f$A=L+D+U\f$ into the lower triangular, diagonal, and upper triangular part
 *  the operator \f$\left(\frac{D}\omega+L\right)\frac\omega{2-\omega}D^{-1}\left(\frac{D}\omega+U\right)\f$,
 *  which is equivalent to computing one single SSOR iterate starting with the initial vector 0.
 *  Note, the SOR-iteration to solve the system \f$Ax=b\f$ is given by
 *  \f$\left(\frac{D}\omega+L\right)x^{m+1}=\frac{2-\omega}\omega Dx^m-\left(\frac{D}\omega+U\right)x^m+b\f$,
 *  the backward SOR-iteration by
 *  \f$\left(\frac{D}\omega+U\right)x^{m+1}=\frac{2-\omega}\omega Dx^m-\left(\frac{D}\omega+L\right)x^m+b\f$,
 *  and the SSOR-iteration by alternatingly applying the SOR-iteration first and then the backward SOR-iteration.
 *  \author Droske
 *  \ingroup preconditioner
 */
template <typename RealType, typename MatrixType>
class SSORPreconditioner< aol::Vector<RealType>, MatrixType > : public Op<aol::Vector<RealType> > {
protected:
  const MatrixType &_sparse;
  RealType _omega;
public:
  SSORPreconditioner ( const MatrixType &SparseMatrix, RealType Omega = 1.2 )
      : _sparse ( SparseMatrix ), _omega ( Omega ) {}

  virtual void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp ( Dest.size() );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  virtual void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> &x = Dest;
    const aol::Vector<RealType> &b = Arg;

    vector<typename Row<RealType>::RowEntry > vec;
    for ( int i = 0; i < x.size(); i++ ) {
      _sparse.makeRowEntries ( vec, i );

#ifdef CHECK_ZERO_DIAG
      RealType v = 0., diag = 1.;
#else
      RealType v = 0., diag = 0.;
#endif
      for ( typename vector<typename Row<RealType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
#ifdef BOUNDS_CHECK
        if ( it->col >= x.size() || it->col < 0 ) {
          cerr << "error: it-> col = " <<  it-> col << endl;
        }
#endif
        if ( it->col == i ) {
          diag = it->value;
#ifdef CHECK_ZERO_DIAG
          if ( Abs ( diag ) < 1e-20 ) {
            diag = 1.;
          }
#endif
        } else if ( it->col < i ) {
          v += it->value * x[it->col];
        }
      }
      x[i] = ( b[i] - _omega * v ) / diag;
    }

    for ( int i = x.size() - 1; i >= 0; --i ) {
      _sparse.makeRowEntries ( vec, i );

#ifdef CHECK_ZERO_DIAG
      RealType v = 0., diag = 1.;
#else
      RealType v = 0., diag = 0.;
#endif
      for ( typename vector<typename Row<RealType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
        if ( it->col == i ) {
          diag = it->value;
#ifdef CHECK_ZERO_DIAG
          if ( Abs ( diag ) < 1e-20 ) {
            diag = 1.;
          }
#endif
          x[i] *= diag;
        } else if ( it->col > i ) {
          v += it->value * x[it->col];
        }
      }
      x[i] = ( x[i] - _omega * v ) / diag;
    }
  }
};

/** \brief SSOR preconditioner for block matrices.
 *  Represents for a given matrix \f$A\f$ and its decomposition
 *  \f$A=L+D+U\f$ into the lower triangular, diagonal, and upper triangular part
 *  the operator \f$\left(\frac{D}\omega+L\right)\frac\omega{2-\omega}D^{-1}\left(\frac{D}\omega+U\right)\f$,
 *  which is equivalent to computing one single SSOR iterate starting with the initial vector 0.
 *  Note, the SOR-iteration to solve the system \f$Ax=b\f$ is given by
 *  \f$\left(\frac{D}\omega+L\right)x^{m+1}=\frac{2-\omega}\omega Dx^m-\left(\frac{D}\omega+U\right)x^m+b\f$,
 *  the backward SOR-iteration by
 *  \f$\left(\frac{D}\omega+U\right)x^{m+1}=\frac{2-\omega}\omega Dx^m-\left(\frac{D}\omega+L\right)x^m+b\f$,
 *  and the SSOR-iteration by alternatingly applying the SOR-iteration first and then the backward SOR-iteration.
 *  \author Wirth
 *  \ingroup preconditioner
 */
template <typename RealType, typename MatrixType>
class SSORPreconditioner< aol::MultiVector<RealType>, MatrixType > :
public aol::Op<MultiVector<RealType> > {
protected:
  const aol::BlockOpBase<RealType, MatrixType> &_blockMat;
  RealType _omega;
public:
  SSORPreconditioner ( const aol::BlockOpBase<RealType, MatrixType> &BlockMat, RealType Omega = 1.2 ) :
    _blockMat( BlockMat ),
    _omega( Omega ) {}

  virtual void applyAdd ( const MultiVector<RealType> &Arg, MultiVector<RealType> &Dest ) const {
    MultiVector<RealType> tmp ( Dest, aol::STRUCT_COPY );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  virtual void apply ( const MultiVector<RealType> &Arg, MultiVector<RealType> &Dest ) const {
    // the following notation is used: L,D,U - lower triangular, diagonal, upper triangular part of _blockMat;
    // b - Arg; A = L+D+U

    vector<typename Row<RealType>::RowEntry > matRow;

    // solve (\omega L+D)y=b for y by forward substitution and write y into Dest
    for ( int k = 0; k < _blockMat.getNumRows(); k++ ) {
      // compute \sum_{j=1}^{k-1}A_{kj}y_j
      Vector<RealType> sum( Dest[k], aol::STRUCT_COPY );
      for ( int j = 0; j < k; j++ )
        _blockMat.getReference( k, j ).applyAdd( Dest[j], sum );
      // solve (\omega L_{kk}+D_{kk})y_k=b_k-\omega\sum_{j=1}^{k-1}A_{kj}y_j by forward substitution
      for ( int i = 0; i < Dest[k].size(); i++ ) {
        _blockMat.getReference( k, k ).makeRowEntries( matRow, i );
        RealType v = 0., diag = 1.;
        for ( typename vector<typename Row<RealType>::RowEntry >::iterator it = matRow.begin(); it != matRow.end(); ++it )
          if ( it->col < i )
            v += it->value * Dest[k][it->col];
          else if ( it->col == i )
            diag = it->value;
        Dest[k][i] = ( Arg[k][i] - _omega * ( sum[i] + v ) ) / diag;
      }
    }

    // solve (D+\omega U)x=Dy for x by backward substitution and write x into Dest
    for ( int k = _blockMat.getNumRows() - 1; k >= 0; k-- ) {
      // compute \sum_{j=k+1}^{n}A_{kj}x_j
      Vector<RealType> sum( Dest[k], aol::STRUCT_COPY );
      for ( int j = k+1; j < _blockMat.getNumCols(); j++ )
        _blockMat.getReference( k, j ).applyAdd( Dest[j], sum );
      // solve (D_{kk}+\omega U_{kk})x_k=D_{kk}y_k-\omega\sum_{j=k+1}^{n}A_{kj}x_j by backward substitution
      for ( int i = Dest[k].size() - 1; i >= 0; --i ) {
        _blockMat.getReference( k, k ).makeRowEntries( matRow, i );
        RealType v = 0., diag = 1.;
        for ( typename vector<typename Row<RealType>::RowEntry >::iterator it = matRow.begin(); it != matRow.end(); ++it )
          if ( it->col > i )
            v += it->value * Dest[k][it->col];
          else if ( it->col == i )
            diag = it->value;
        Dest[k][i] = ( diag * Dest[k][i] - _omega * ( sum[i] + v ) ) / diag;
      }
    }
  }
};


/** \brief SSOR preconditioning which is naively implemented by performing
 *              one forward and one backward Block-Gauss-Seidel sweep. This does
 *              not have optimal complexity.
 *  \todo think about initialization of Dest (OS,BW)
 *  \author Schwen
 *  \ingroup preconditioner
 */
template< typename VectorType, typename MatrixType >
class GaussSeidelPreconditioner : public Op< VectorType > {
protected:
  typedef typename VectorType::DataType DataType;
  aol::GaussSeidelSweeper< DataType, VectorType, MatrixType > _sweeper;
  aol::GaussSeidelSweepingMode _firstSweepMode, _secondSweepMode;

public:
  explicit GaussSeidelPreconditioner ( const MatrixType &Matrix, const DataType Relax = 1.2 ) : _sweeper ( Matrix, Relax ), _firstSweepMode ( aol::GAUSS_SEIDEL_FORWARD ), _secondSweepMode ( aol::GAUSS_SEIDEL_BACKWARD ) {
  }

  void setSweepingModes ( const aol::GaussSeidelSweepingMode firstMode, const aol::GaussSeidelSweepingMode secondMode ) {
    _firstSweepMode = firstMode;
    _secondSweepMode = secondMode;
  }

  virtual void applyAdd ( const VectorType &Arg, VectorType &Dest ) const {
    VectorType tmp ( Dest, aol::DEEP_COPY );
    // ? VectorType tmp ( Dest, aol::STRUCT_COPY );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    // ? Dest.setZero();
    _sweeper.apply ( Arg, Dest, _firstSweepMode );
    _sweeper.apply ( Arg, Dest, _secondSweepMode );
  }
};


/** \brief A class for preconditioning sparse matrices by ILU(0).
 *  \attention The preconditioner must compute a decomposition (in the constructor), so it uses a copy of the matrix
 *             The matrix must be filled at time of the precoditioners construction
 *             In the contructor get (,) is used for the needed column values, values are written via set/add (,)
 *  \author unknown
 *  \ingroup preconditioner
 */
template <typename RealType, typename MatrixType>
class ILU0Preconditioner : public Op<aol::Vector<RealType> > {

protected:
  MatrixType _decomp;

public:
  ILU0Preconditioner ( const MatrixType& matrix )
      : _decomp ( matrix ) {
    vector<typename Row<RealType>::RowEntry > rowi;
    vector<typename Row<RealType>::RowEntry > rowk;

    for ( int i = 1; i < _decomp.getNumRows (); ++i ) { // Nothing is done in first row
      _decomp.makeSortedRowEntries ( rowi, i ); // Rowwise elimination

      for ( typename vector<typename Row<RealType>::RowEntry >::iterator ikentry = rowi.begin ();
            ikentry != rowi.end () && ikentry->col < i; ++ikentry ) { // Only left from diagonal, RowEntries are sorted
        int k = ikentry->col;

        RealType ikval = ikentry->value;
        _decomp.makeSortedRowEntries ( rowk, k );

        typename vector<typename Row<RealType>::RowEntry >::iterator kjentry = rowk.begin ();

        // Find index kk
        while ( kjentry != rowk.end () && kjentry->col < k ) ++kjentry;
        if ( kjentry == rowk.end () || kjentry->col != k ) throw Exception ( "Pivot zero in ILU", __FILE__, __LINE__ );

        // Eliminate index ik, store factor
        ikval /= kjentry->value;
        _decomp.set ( i, k, ikval );

        typename vector<typename Row<RealType>::RowEntry >::iterator ikentrynext = ikentry; ++ikentrynext; // Only right from current column
        for ( typename vector<typename Row<RealType>::RowEntry >::iterator ijentry = ikentrynext; ijentry != rowi.end (); ++ijentry ) {
          int j = ijentry->col;

          // Find index kj
          while ( kjentry != rowk.end () && kjentry->col < j ) ++kjentry;
          if ( kjentry == rowk.end () || kjentry->col != j ) continue; // Zero entry
          RealType kjval = kjentry->value;
          RealType newijval = - ikval * kjval;
          _decomp.add ( i, j, newijval );
          ijentry->value += newijval;
        }
      }
    }
  }

  virtual ~ILU0Preconditioner () {}

  virtual void applyAdd ( const Vector<RealType>& arg, Vector<RealType>& dest ) const {
    Vector<RealType> tmp ( dest.size () );
    apply ( arg, tmp );
    dest += tmp;
  }

  virtual void apply ( const Vector<RealType> &arg, Vector<RealType> &dest ) const {
    vector<typename Row<RealType>::RowEntry > row;

    // Forwards-substitution
    for ( int i = 0; i < _decomp.getNumRows (); ++i ) {
      _decomp.makeSortedRowEntries ( row, i ); // Sorted to allow early break
      typename vector<typename Row<RealType>::RowEntry >::iterator end = row.end ();
      RealType sum = 0;

      for ( typename vector<typename Row<RealType>::RowEntry >::iterator kentry = row.begin (); kentry != end && kentry->col < i; ++kentry ) {
        int k = kentry->col;

        sum += dest [k] * kentry->value;
      }
      dest [i] = arg [i] - sum;
    }

    // Backwards-substitution
    for ( int i = _decomp.getNumRows () - 1; i >= 0; --i ) {
      _decomp.makeRowEntries ( row, i ); // Need not be sorted
      RealType sum = 0;
      RealType diag = 0;

      // Strangely this crashes with the intel-CC if using a normal iterator loop (like above)
      for ( int kentry = 0, end = row.size (); kentry < end; kentry++ ) {
        const int k = row[kentry].col;

        if ( k < i ) continue; // U part
        if ( k == i ) {
          diag = row [kentry].value;
          continue;
        }
        sum += dest [k] * row[kentry].value;
      }
      dest [i] -= sum; dest [i] /= diag;
    }
  }
};

template <typename RealType, typename BlockMatrixType>
class ILU0BlockPreconditioner : public Op<aol::MultiVector<RealType> > {

protected:
  BlockMatrixType _decomp;

public:
  ILU0BlockPreconditioner ( const BlockMatrixType& matrix )
      : _decomp ( matrix ) {
    vector<typename Row<RealType>::RowEntry > rowi;
    vector<typename Row<RealType>::RowEntry > rowk;

    for ( int i = 0; i < _decomp.getNumRows (); ++i ) { // Nothing is done in first row
      for ( int locRow = 0; locRow < _decomp.getReference(0, 0).getNumRows(); locRow++ ) {
        if ( !(i==0 && locRow==0) ) {
          _decomp.makeUnblockedSortedRowEntries( rowi, i, locRow );

          for ( typename vector<typename Row<RealType>::RowEntry >::iterator ikentry = rowi.begin ();
                ikentry != rowi.end () && ikentry->col < i*_decomp.getReference(0, 0).getNumRows()+locRow; ++ikentry ) { // Only left from diagonal, RowEntries are sorted
            int k = ikentry->col;

            RealType ikval = ikentry->value;
            _decomp.makeUnblockedSortedRowEntries( rowk, floor(k/_decomp.getReference(0, 0).getNumRows()), k%_decomp.getReference(0, 0).getNumRows() );

            typename vector<typename Row<RealType>::RowEntry >::iterator kjentry = rowk.begin ();

            // Find index kk
            while ( kjentry != rowk.end () && kjentry->col < k ) ++kjentry;
            if ( kjentry == rowk.end () || kjentry->col != k ) throw Exception ( "Pivot zero in ILU", __FILE__, __LINE__ );

            // Eliminate index ik, store factor
            ikval /= kjentry->value;

            _decomp.getReference(i, floor(k/_decomp.getReference(0, 0).getNumCols())).set ( locRow, k%_decomp.getReference(0, 0).getNumCols(), ikval );

            typename vector<typename Row<RealType>::RowEntry >::iterator ikentrynext = ikentry; ++ikentrynext; // Only right from current column
            for ( typename vector<typename Row<RealType>::RowEntry >::iterator ijentry = ikentrynext; ijentry != rowi.end (); ++ijentry ) {
              int j = ijentry->col;

              // Find index kj
              while ( kjentry != rowk.end () && kjentry->col < j ) ++kjentry;
              if ( kjentry == rowk.end () || kjentry->col != j ) continue; // Zero entry
              RealType kjval = kjentry->value;
              RealType newijval = - ikval * kjval;
              _decomp.getReference(i, floor(j/_decomp.getReference(0, 0).getNumCols())).add ( locRow, j%_decomp.getReference(0, 0).getNumCols(), newijval );
              ijentry->value += newijval;
            }
          }
        }
      }
    }
//   _decomp.getPointer(0,0)->print(cerr);
// _decomp.getPointer(0,1)->print(cerr);
// _decomp.getPointer(1,0)->print(cerr);
// _decomp.getPointer(1,1)->print(cerr);
  }

  virtual ~ILU0BlockPreconditioner () {}

  virtual void applyAdd ( const MultiVector<RealType>& arg, MultiVector<RealType>& dest ) const {
    MultiVector<RealType> tmp ( dest, STRUCT_COPY );
    apply ( arg, tmp );
    dest += tmp;
  }

  virtual void apply ( const MultiVector<RealType> &arg, MultiVector<RealType> &dest ) const {
    vector<typename Row<RealType>::RowEntry > row;

    // Forwards-substitution
    for ( int i = 0; i < _decomp.getNumRows (); ++i ) {
      for ( int locRow = 0; locRow < _decomp.getReference(0, 0).getNumRows(); locRow++ ) {
        _decomp.makeUnblockedSortedRowEntries( row, i, locRow );

        typename vector<typename Row<RealType>::RowEntry >::iterator end = row.end ();
        RealType sum = 0;

        for ( typename vector<typename Row<RealType>::RowEntry >::iterator kentry = row.begin (); kentry != end && kentry->col < i*_decomp.getReference(0, 0).getNumRows()+locRow; ++kentry ) {
          int k = kentry->col;

          sum += dest [floor(k/_decomp.getReference(0, 0).getNumRows())][k%_decomp.getReference(0, 0).getNumRows()] * kentry->value;
        }
        dest [i][locRow] = arg [i][locRow] - sum;
      }
    }

    // Backwards-substitution
    for ( int i = _decomp.getNumRows () - 1; i >= 0; --i ) {
      for ( int locRow = _decomp.getReference(0, 0).getNumRows() - 1; locRow >= 0; --locRow ) {
        _decomp.makeUnblockedRowEntries ( row, i, locRow ); // Need not be sorted
        RealType sum = 0;
        RealType diag = 0;

        // Strangely this crashes with the intel-CC if using a normal iterator loop (like above)
        for ( int kentry = 0, end = row.size (); kentry < end; kentry++ ) {
          const int k = row[kentry].col;

          if ( k < i*_decomp.getReference(0, 0).getNumRows()+locRow ) continue; // U part
          if ( k == i*_decomp.getReference(0, 0).getNumRows()+locRow ) {
            diag = row [kentry].value;
            continue;
          }
          sum += dest [floor(k/_decomp.getReference(0, 0).getNumRows())][k%_decomp.getReference(0, 0).getNumRows()] * row[kentry].value;
        }
        dest [i][locRow] -= sum;
        dest [i][locRow] /= diag;
      }
    }
  }
};


/** \brief A class for block-diagonally preconditioning BlockMatrices with diagonal blocks all the same size and all set.
 *  \author Schwen
 *  \ingroup preconditioner
 *  This preconditioner emulates rearranging the unknowns such that x, y and z component of function at one grid point are
 *  subsequent unknowns in the system and then inverts d x d (d = 2, 3) blocks along the diagonal, then the rearranging is undone.
 */

template< typename RealType, qc::Dimension Dim >
class BlockDiagonalPreconditioner : public Op< aol::MultiVector< RealType > > {

protected:
  aol::BlockMatrix< aol::DiagonalMatrix<RealType> > _blockDiagBlockInv;

public:
  template< typename BlockMatrixType >
  explicit BlockDiagonalPreconditioner ( const BlockMatrixType &BlockMat ) :
    _blockDiagBlockInv ( Dim, Dim, BlockMat.getReference ( 0, 0 ).getNumRows(), BlockMat.getReference ( 0, 0 ).getNumCols() ) {

    setBlockDiagBlockInv ( BlockMat, _blockDiagBlockInv );
  }


  virtual void applyAdd ( const MultiVector<RealType> &Arg, MultiVector<RealType> &Dest ) const {
    _blockDiagBlockInv.applyAdd ( Arg, Dest );
  }

  virtual void apply ( const MultiVector<RealType> &Arg, MultiVector<RealType> &Dest ) const {
    _blockDiagBlockInv.apply ( Arg, Dest );
  }

  void getPreconditionerBlockMatrix ( aol::BlockMatrix < aol::DiagonalMatrix< RealType > > &BlockMat ) const {
    if ( BlockMat.getNumCols() != Dim || BlockMat.getNumRows() != Dim )
      throw aol::Exception ( "aol::BlockDiagonalPreconditioner: illegal number of blocks", __FILE__, __LINE__ );

    for ( int i = 0; i < Dim; ++i ) {
      for ( int j = 0; j < Dim; ++j ) {
        BlockMat.getReference(i, j) = _blockDiagBlockInv.getReference(i, j);
      }
    }
  }

  const aol::BlockMatrix< aol::DiagonalMatrix<RealType> >& getPreconditionerBlockMatrixRef ( ) const {
    return ( _blockDiagBlockInv );
  }

  template< typename BlockMatrixType >
  static void setBlockDiagBlockInv ( const BlockMatrixType &BlockMat, aol::BlockMatrix< aol::DiagonalMatrix<RealType> > &bdbInv ) {
    if ( BlockMat.getNumCols() != Dim || BlockMat.getNumRows() != Dim )
      throw aol::Exception ( "aol::BlockDiagonalPreconditioner::setBlockDiagBlockInv: illegal number of blocks", __FILE__, __LINE__ );

    const int n = BlockMat.getReference ( 0, 0 ).getNumRows(); // we assume squared blocks of all the same size and rely on bounds checking

    typename aol::MatDimTrait<RealType, Dim>::MatType diagBlock, diagBlockInv;

    for ( int i = 0; i < n; ++i ) {
      for ( short j = 0; j < Dim; ++j ) {
        for ( short k = 0; k < Dim; ++k ) {
          diagBlock.set ( j, k, BlockMat.getReference ( j, k ).getDiag ( i ) );
        }
      }

      if ( diagBlock.det() != aol::NumberTrait<RealType>::zero ) {
        diagBlockInv = diagBlock.inverse();
      } else {
        diagBlockInv.setIdentity();
      }

      for ( short j = 0; j < Dim; ++j ) {
        for ( short k = 0; k < Dim; ++k ) {
          bdbInv.getReference ( j, k ).set( i, i, diagBlockInv.get ( j, k ) );
        }
      }
    }

  }

};

/** \brief Block-SSOR preconditioning which is naively implemented by
 *         performing one forward and one backward relaxed Block-Gauss-Seidel
 *         sweep. This does not have optimal complexity.
 *
 *  Currently, it only works for 3D (due to the BlockGaussSeidelSweeper).
 *  Can be parallelized by using parallelized Gauss-Seidel iterations
 *  (which is no longer SSOR!)
 *  \author Schwen
 *  \ingroup preconditioner
 */
template< typename VectorType, typename BlockOpType >
class BlockGaussSeidelPreconditioner : public Op< VectorType > {
protected:
  typedef typename VectorType::DataType DataType;
  aol::BlockGaussSeidelSweeper< DataType, BlockOpType > _sweeper;
  aol::GaussSeidelSweepingMode _firstSweepMode, _secondSweepMode;

public:
  explicit BlockGaussSeidelPreconditioner ( const BlockOpType &BlockMatrix, const DataType Relax = 1.2 ) : _sweeper ( BlockMatrix, Relax ), _firstSweepMode ( aol::GAUSS_SEIDEL_FORWARD ), _secondSweepMode ( aol::GAUSS_SEIDEL_BACKWARD ) {
  }

  void setSweepingModes ( const aol::GaussSeidelSweepingMode firstMode, const aol::GaussSeidelSweepingMode secondMode ) {
    _firstSweepMode = firstMode;
    _secondSweepMode = secondMode;
  }

  virtual void applyAdd ( const VectorType &Arg, VectorType &Dest ) const {
    VectorType tmp ( Dest, aol::DEEP_COPY );
    // ? VectorType tmp ( Dest, aol::STRUCT_COPY );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    // ? Dest.setZero()
    _sweeper.apply ( Arg, Dest, _firstSweepMode );
    _sweeper.apply ( Arg, Dest, _secondSweepMode );
  }
};


/** \brief Block variant of the SSOR preconditioner.
 *  \author Schwen (based on code by Droske and Wirth)
 *  \ingroup preconditioner
 */
template< typename RealType, typename BlockMatrixType, qc::Dimension Dim >
class BlockSSORPreconditioner : public Op< aol::MultiVector< RealType > > {
protected:
  std::vector< typename aol::MatDimTrait<RealType,Dim>::MatType > _invDiags; // for cacheing inverses of the diagonal blocks
  const BlockMatrixType &_bMat;
  RealType _omega;

public:
  BlockSSORPreconditioner ( const BlockMatrixType &BlockMatrix, const RealType Omega = 1.2f )
    : _invDiags ( ), _bMat ( BlockMatrix ), _omega ( Omega ) {

    cacheDiagInverses ( _bMat, _invDiags );

  }

  virtual void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    aol::MultiVector<RealType> tmp ( Dest, aol::STRUCT_COPY );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  virtual void apply ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {

    for ( int i = 0; i < Dest[0].size(); ++i ) { // assuming all have the same size
      aol::Vec<Dim,RealType> v;
      for ( short a = 0; a < Dim; ++a ) {
        for ( short b = 0; b < Dim; ++b ) {
          std::vector< typename aol::Row<RealType>::RowEntry > vec;
          this->_bMat.getReference ( a, b ).makeSortedRowEntries ( vec, i );

          for ( typename std::vector<typename aol::Row<RealType>::RowEntry >::iterator it = vec.begin(); ( it != vec.end() ) && ( it->col < i ); ++it ) {
            v[ a ] += it->value * Dest[ b ][ it->col ];
          }
        }
      }
      aol::Vec<Dim,RealType> tmp, tmpprod;
      for ( short a = 0; a < Dim; ++a )
        tmp[a] = Arg[a][i];
      tmp -= _omega * v;
      _invDiags[i].mult ( tmp, tmpprod );
      for ( short a = 0; a < Dim; ++a )
        Dest[a][i] = tmpprod[a];
    }

    for ( int i = Dest[0].size() - 1; i >= 0; --i ) {
      aol::Vec<Dim,RealType> v;
      for ( short a = 0; a < Dim; ++a ) {
        for ( short b = 0; b < Dim; ++b ) {
          std::vector< typename aol::Row<RealType>::RowEntry > vec;
          this->_bMat.getReference ( a, b ).makeRowEntries ( vec, i ); // sortedness not necessary here

          for ( typename std::vector<typename aol::Row<RealType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
            if ( it->col > i ) {
              v[ a ] += it->value * Dest[ b ][ it->col ];
            }
          }
        }
      }
      v *= _omega;
      aol::Vec<Dim,RealType> tmpmlt;
      _invDiags[i].mult ( v, tmpmlt );
      for ( short a = 0; a < Dim; ++a )
        Dest[a][i] -= tmpmlt[a];
    }
  }

  static void cacheDiagInverses ( const BlockMatrixType &BlockMat, std::vector< typename aol::MatDimTrait<RealType,Dim>::MatType > &invDiags ) {
    const int n = BlockMat.getReference(0,0).getNumRows(); // assuming this is correct ...
    invDiags.resize ( n );

    for ( int row = 0; row < n; ++row ) {
      typename aol::MatDimTrait<RealType,Dim>::MatType diag;

      for ( short a = 0; a < Dim; ++a ) {
        for ( short b = 0; b < Dim; ++b ) {
          std::vector< typename aol::Row<RealType>::RowEntry > vec;
          BlockMat.getReference ( a, b ).makeRowEntries ( vec, row );
          for ( typename std::vector<typename aol::Row<RealType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
            if ( it->col == row ) {
              diag[a][b] = it->value;
              break;
            }
          }
        }
      }

      if ( diag.det() != 0.0 )
        invDiags[row] = diag.inverse();
      else
        invDiags[row].setIdentity();
    }
  }
};


} // end namespace

#endif
