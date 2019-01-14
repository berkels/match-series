#ifndef __OP_H
#define __OP_H

#include <aol.h>
#include <vec.h>
#include <multiVector.h>
#include <bitVector.h>
#include <randomGenerator.h>
#include <progressBar.h>

namespace aol {

template <class DomainType, class RangeType> class LinCombOp;

/** An abstract operator on vectors.
 * \author Droske
 */
template <typename _DomainType, typename _RangeType = _DomainType>
class Op {
  friend class LinCombOp<_DomainType, _RangeType>;

public:
  typedef _DomainType DomainType;
  typedef _RangeType  RangeType;

  Op() { }

  // Destroy polymorphic Ops correctly, important!
  virtual ~Op () {}

  virtual ostream& print ( ostream& /*out*/ ) const {
    throw ( Exception ( "Print not possible on abstract class Op.", __FILE__, __LINE__ ) );
  }

  virtual void operator() ( const DomainType &Arg, RangeType &Dest ) const {
    this->apply ( Arg, Dest );
  }

  virtual void apply ( const DomainType &Arg, RangeType &Dest ) const {
    Dest.setZero ();
    applyAdd ( Arg, Dest );
  }

  //! when using applyAddMasked, implement applyAdd as
  //! applyAddMasked(Arg, Dest, INCLUDE_WRITE_DEFAULT).
  virtual void applyAdd ( const DomainType &Arg, RangeType &Dest ) const = 0;

  virtual void applySingle ( DomainType& ) const {
    throw aol::UnimplementedCodeException ( "aol::Op<>::applySingle ( Arg ) may be overloaded on classes that have a useful behavior for one-argument apply.", __FILE__, __LINE__ );
  }

  virtual void applyAddSingle ( DomainType& ) const {
    throw aol::UnimplementedCodeException ( "aol::Op<>::applyAddSingle ( Arg ) may be overloaded on classes that have a useful behavior for one-argument applyAdd.", __FILE__, __LINE__ );
  }

};

//! \brief Extends the abstract operator aol::Op, requiring a RealType in _RangeType to introduce a applyAddMultiple method
//! \author Toelkes
template <typename _DomainType, typename _RangeType = _DomainType, typename DataType = typename _RangeType::DataType >
class OpWithDataType : public Op<_DomainType, _RangeType> {
public:
  typedef _DomainType DomainType;
  typedef _RangeType  RangeType;

  virtual ~OpWithDataType () {}

  virtual void applyAddMultiple ( const DomainType &Arg, RangeType &Dest, DataType factor = 1.0 ) const = 0;
};

template <class DataType>
class DiagonalBlockOp : public Op<MultiVector<DataType> > {
public:
  DiagonalBlockOp ( const Op<Vector<DataType> > &op )
      : _op ( op ) {  }

  void applyAdd ( const MultiVector<DataType> &arg, MultiVector<DataType> &dest ) const {
    if ( !arg.compareDim ( dest ) ) {
      throw ( aol::Exception ( "dimensions of arg and dest differ", __FILE__, __LINE__ ) );
    }

    if ( &arg == &dest ) {
      throw ( aol::Exception ( "arg and dest should not be the same objects!\n", __FILE__, __LINE__ ) );
    }

    for ( int i = 0; i < static_cast<int> ( arg.numComponents() ); ++i ) {
      _op.applyAdd ( arg[i], dest[i] );
    }
  }

private:
  const Op<Vector<DataType> > &_op;
};

//! \brief Interface for an operator with adjoint
//! \author Effland
template <typename _DomainType, typename _RangeType = _DomainType>
class OpWithAdjoint : public Op< _DomainType, _RangeType > {
public:
  OpWithAdjoint ( ) { }

  virtual ~OpWithAdjoint () {}

  virtual void applyAdjoint ( const _RangeType &Arg, _DomainType &Dest ) const {
    Dest.setZero ();
    this->applyAddAdjoint ( Arg, Dest );
  }

  virtual void applyAddAdjoint ( const _RangeType &Arg, _DomainType &Dest ) const = 0;
};
  
//! \brief Interface for an operator with trasposed matrix assembled
//! \author Tatano
template <typename _DomainType, typename _RangeType = _DomainType>
class OpWithTransposedAssembled : public Op< _DomainType, _RangeType > {
public:
  OpWithTransposedAssembled ( ) { }
    
  virtual ~OpWithTransposedAssembled () {}
    
  virtual void applyAndAssembleTransposed ( const _DomainType &Arg, _RangeType &Dest, _RangeType &DestTrans ) const = 0;
};

/** \brief Useless default template for specialization.
 *  Using the *Base class is a workaround to be able to control parallelization via a template parameter with default behavior (which is not possible when using template specialization)
 */
template < typename VectorType, bool Parallelize = true  >
class BiOpBase : public Op<VectorType> { };

/*! \brief Derive from this class to incorporate standard behaviour for multivectors.
 */
template < typename RealType, bool Parallelize >
class BiOpBase< Vector<RealType>, Parallelize > : public Op<Vector<RealType> > {
public:

  BiOpBase ( ) {
  }

  virtual ~BiOpBase ( ) {
  }

  //! Vector-applyAdd called on all components of MultiVector
  void applyAdd ( const MultiVector<RealType> &arg, MultiVector<RealType> &dest ) const {
#ifdef _OPENMP
#pragma omp parallel for if ( Parallelize )
#endif
    for ( int i = 0; i < arg.numComponents(); ++i ) {
      applyAdd ( arg[i], dest[i] );
    }
  }

  //! Vector-apply called on all components of MultiVector
  void apply ( const MultiVector<RealType> &arg, MultiVector<RealType> &dest ) const {
#ifdef _OPENMP
#pragma omp parallel for if ( Parallelize )
#endif
    for ( int i = 0; i < arg.numComponents(); ++i ) {
      apply ( arg[i], dest[i] );
    }
  }

  //! Vector-applySingle called on all components of MultiVector
  void applySingle ( MultiVector<RealType> &dest ) const {
#ifdef _OPENMP
#pragma omp parallel for if ( Parallelize )
#endif
    for ( int i = 0; i < dest.numComponents(); ++i ) {
      applySingle ( dest[i] );
    }
  }

  //! Vector-applySingle called on all components of MultiVector
  void applyAddSingle ( MultiVector<RealType> &dest ) const {
#ifdef _OPENMP
#pragma omp parallel for if ( Parallelize )
#endif
    for ( int i = 0; i < dest.numComponents(); ++i ) {
      applyAddSingle ( dest[i] );
    }
  }

  using aol::Op< Vector<RealType> >::applySingle;
  using aol::Op< Vector<RealType> >::applyAddSingle;

  //! Vector-applyAdd. Implement this method on derived classes.
  virtual void applyAdd ( const Vector<RealType> &arg, Vector<RealType> &dest ) const  = 0;

  //! Vector-apply. This method is virtual on Op and may be overloaded on derived classes if apply is the standard method for the operator (rather than applyAdd)
  void apply ( const Vector<RealType> &arg, Vector<RealType> &dest ) const {
    Op<Vector<RealType> >::apply ( arg, dest );
  }

};

/*! \brief Derive from this class to incorporate standard behaviour for vectors.
 *  Use this BiOp only if you know what you are doing and compile with --expert :-)
 */
template < typename RealType, bool ParallelizeUnused >
class BiOpBase< MultiVector<RealType>, ParallelizeUnused > : public Op< MultiVector<RealType> > {
public:

  BiOpBase ( ) {
  }

  virtual ~BiOpBase ( ) {
  }

  //! MultiVector-applyAdd. Implement this method on derived classes.
  virtual void applyAdd ( const MultiVector<RealType> &arg, MultiVector<RealType> &dest ) const = 0;

  //! MultiVector-apply. This method is virtual on Op and may be overloaded on derived classes if apply is the standard method for the operator (rather than applyAdd)
  void apply ( const MultiVector<RealType> &arg, MultiVector<RealType> &dest ) const {
    Op<MultiVector<RealType> >::apply ( arg, dest );
  }

  //! Vector-applyAdd calling applyAdd on wrapper-MultiVectors for the arguments passed.
  void applyAdd ( const Vector<RealType> &arg, Vector<RealType> &dest ) const {
    MultiVector<RealType> Arg ( 0, 0 ), Dest ( 0, 0 );
    Arg.appendReference ( arg );
    Dest.appendReference ( dest );
    applyAdd ( Arg, Dest );
  }

  //! Vector-apply calling applyAdd on wrapper-MultiVectors for the arguments passed.
  void apply ( const Vector<RealType> &arg, Vector<RealType> &dest ) const {
    MultiVector<RealType> Arg ( 0, 0 ), Dest ( 0, 0 );
    Arg.appendReference ( arg );
    Dest.appendReference ( dest );
    apply ( Arg, Dest );
  }

};


template < typename VectorType, bool Parallelize = true>
class BiOp : public BiOpBase< VectorType, Parallelize >{

};


/** A class for the combination of different Operators.
 * \author Droske
 */
template <typename DomainType, typename RangeType = DomainType>
class LinCombOp : public Op<DomainType, RangeType> {
public:
  typedef Op<DomainType, RangeType>         OpDR;
  typedef typename RangeType::DataType      DataType;
  //! Constructor
  LinCombOp() : Op<DomainType, RangeType>() {}

  virtual ~LinCombOp() { }

  //! apply the linear combination to Arg (multicomponent-version).
  void applyAdd ( const DomainType &Arg, RangeType &Dest ) const {
    typename list<const aol::Op<DomainType, RangeType>* >::const_iterator         opIt;
    typename list<DataType>::const_iterator                                    coeffIt;

    RangeType tmp ( Dest );

    // qc::Compute all other operators
    for ( opIt = opList.begin(), coeffIt = opCoeff.begin();
          opIt != opList.end(); ++opIt, ++coeffIt ) {
      tmp.setZero();
      ( *opIt )->apply ( Arg, tmp );

      if ( *coeffIt != 1.0 ) {
        // tmp *= (*coeffIt);
        Dest.addMultiple ( tmp, *coeffIt );
      } else {
        Dest += tmp;
      }
    }
  }

  //! Append reference to an arbitrary operator.
  void appendReference ( const OpDR &Op, DataType coeff = aol::NumberTrait<DataType>::one ) {
    opList.push_back ( &Op );
    opCoeff.push_back ( coeff );
  }

  //! Remove all operators.
  void clearAll() {
    //clearCompWiseList( );
    clearOpList();
  }

  //! Remove all other operators from list
  void clearOpList() {
    opList.erase ( opList.begin(), opList.end() );
    opCoeff.erase ( opCoeff.begin(), opCoeff.end() );
  }

  void setZero(){
    clearOpList();
  }

  const list<const OpDR*> &getOpList() const {
    return opList;
  }

  //!
  list<const OpDR*> &getOpList() {
    return opList;
  }

  const list<DataType> &getCoeffList() const {
    return opCoeff;
  }

  //! dummy function to make second derivative validator work, cf. aol::SecondDerivativeValidator<>
  //! \todo re-consider this!
  void transpose() const {
    throw aol::Exception("aol::LinCombOp::transpose: not implemented", __FILE__, __LINE__ );
  }

protected:
  list<const OpDR*>    opList;
  list<DataType>       opCoeff;

};


/**
 *  BlockOpBase is a non-abstract basis class for BlockOp, BlockMatrix and possibly others.
 *  It contains pointers to the operators and read access to those pointers.
 *  It can be applied to instances of ArgType, which has to have operator[] implemented
 *  (the single blocks of type OpType then operate on the single components of ArgType, which can e.g. be a MultiVector or VectorContainer).
 *  Write access to the pointers may only be implemented in derived classes that do not do memory management.
 *  Indexing is done via qc::ILexCombine2, note the order of arguments (for historic reasons): j, i,_n.
 *  \author Droske, Schwen
 */
template <typename RealType, typename OpType = const Op< Vector<RealType> >, typename ArgType = MultiVector< RealType > >
class BlockOpBase : public Op<ArgType> {
public:
  typedef OpType CompOpType;

protected:
  int _m, _n;
  vector< OpType* >  _blockEntryPointers;
  bool _parallelizeBlockOp;

public:
  BlockOpBase ( const int M, const int N ) : _m ( M ), _n ( N ), _blockEntryPointers ( M * N ), _parallelizeBlockOp ( true ) {
    setAllPointersToNULL();
  }

  virtual ~BlockOpBase() {}

  int getNumRows() const {
    return _m;
  }
  int getNumCols() const {
    return _n;
  }

#ifdef BOUNDS_CHECK
  inline bool boundsCheck ( const int i, const int j, const char* msg, const char* fi, const int li ) const {
    const bool isIn = ( i >= 0 && i < _m &&
                        j >= 0 && j < _n );
    if ( !isIn ) {
      char errmsg[1024];
      sprintf( errmsg, "%s %d %d (upper bounds: %d %d)", msg, i, j, _m, _n );
      throw aol::OutOfBoundsException ( errmsg, fi, li );
    }
    return ( isIn );
  }
#endif

  const OpType& getReference ( const int I, const int J ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( I, J, "aol::BlockOpBase::getReference: Index out of bounds", __FILE__, __LINE__ );
#endif
    if ( !_blockEntryPointers[ qc::ILexCombine2( J, I, _n )] ) {
      cerr<<"\n _blockEntryPointers["<<I<<"*"<<_n<<"+"<<J<<"] is not set \n";
    }
    return *_blockEntryPointers[ qc::ILexCombine2( J, I, _n ) ];
  }

  const OpType* getPointer ( const int I, const int J ) const {
#ifdef BOUNDS_CHECK
    boundsCheck ( I, J, "aol::BlockOpBase::getPointer: Index out of bounds", __FILE__, __LINE__ );
#endif
    return _blockEntryPointers[ qc::ILexCombine2( J, I, _n ) ];
  }

  OpType& getReference ( const int I, const int J ) {
#ifdef BOUNDS_CHECK
    boundsCheck ( I, J, "aol::BlockOpBase::getReference: Index out of bounds", __FILE__, __LINE__ );
#endif
    if ( !_blockEntryPointers[ qc::ILexCombine2( J, I, _n ) ] ) {
      cerr<<"\n _blockEntryPointers["<<I<<"*"<<_n<<"+"<<J<<"] is not set \n";
    }
    return *_blockEntryPointers[ qc::ILexCombine2( J, I, _n ) ];
  }

  OpType* getPointer ( const int I, const int J ) {
#ifdef BOUNDS_CHECK
    boundsCheck ( I, J, "aol::BlockOpBase::getPointer: Index out of bounds", __FILE__, __LINE__ );
#endif
    return _blockEntryPointers[ qc::ILexCombine2( J, I, _n ) ];
  }

  void setZero() {
    for ( int i = 0; i < _m; ++i )
      for ( int j = 0; j < _n; ++j )
        if ( getPointer ( i, j ) )
          getReference ( i, j ).setZero();
  }

  void applyAdd ( const ArgType &MArg, ArgType &MDest ) const {
#ifdef BOUNDS_CHECK
    if ( MArg.numComponents() != _n || MDest.numComponents() != _m ) {
      throw aol::Exception("aol::BlockOp::applyAdd: MultiVectors do not have correct numbers of components", __FILE__, __LINE__ );
    }
#endif

#ifdef _OPENMP
#pragma omp parallel for if ( _parallelizeBlockOp )
#endif
    for ( int i = 0; i < _m; ++i ) {
      for ( int j = 0; j < _n; ++j ) {
        if ( getPointer ( i, j ) ) {
          getReference ( i, j ).applyAdd ( MArg[j], MDest[i] );
        }
      }
    }
  }

  //! Change size of the BlockOp, setting all pointers to NULL
  virtual void reallocate ( const int M, const int N ) {
    _m = M;
    _n = N;
    _blockEntryPointers.resize ( _m * _n );
    setAllPointersToNULL();
  }

  void setParallelizeBlockOp ( const bool Parallelize ) {
    _parallelizeBlockOp = Parallelize;
  }

protected:
  void setAllPointersToNULL ( ) {
    for ( typename vector< OpType* >::iterator it = this->_blockEntryPointers.begin(); it != this->_blockEntryPointers.end(); ++it ) {
      ( *it ) = NULL ;
    }
  }

  void setPointer ( const int I, const int J, OpType *op ) {
    _blockEntryPointers[ qc::ILexCombine2( J, I, this->_n ) ] = op;
  }

};


/**
 *  A Block Operator consists of several operators in block structure.
 *  Apply works on MultiVectors in the corresponding block-wise manner.
 *
 *  \author Droske, Schwen
 */
template <typename RealType, typename OpType = const Op< Vector<RealType> > >
class BlockOp : public BlockOpBase< RealType, OpType > {

public:
  typedef typename BlockOpBase< RealType, OpType >::CompOpType CompOpType;

public:
  BlockOp ( const int M, const int N ) : BlockOpBase< RealType, OpType >( M, N ) {
  }

  virtual ~BlockOp() {}

  void setReference ( const int I, const int J, OpType &op ) {
    this->_blockEntryPointers[ qc::ILexCombine2( J, I, this->_n ) ] = &op;
  }

  void unset ( const int I, const int J ) {
    this->_blockEntryPointers[ qc::ILexCombine2( J, I, this->_n ) ] = NULL;
  }

  //! This method should be public here
  using BlockOpBase<RealType, OpType>::setAllPointersToNULL;
  using BlockOpBase<RealType, OpType>::setPointer;

};

/** Abstract block operator class. In contrast to BlockOp this class also enables
 *  of arbitrary block- and argument-types. In particular it is possible to create
 *  block ops whose blocks are block ops
 *
 *  \author Preusser
 */
template <typename VectorType, typename OpType>
class ABlockOp  : public Op<VectorType>{
public:

  ABlockOp( const int M, const int N ) : m ( M ), n ( N ), blockMat ( M * N ) {
    typename vector< OpType* >::iterator it;
    for ( it = blockMat.begin(); it != blockMat.end(); ++it ) {
      ( *it ) = NULL ;
    }
  }

  void applyAdd(const VectorType &MArg, VectorType &MDest) const {
    for ( int i = 0; i < m; ++i ) {
      for ( int j = 0; j < n; ++j ) {
        if ( getPointer ( i, j ) ) {
    getReference ( i, j ).applyAdd ( MArg[j], MDest[i] );
  }
      }
    }
  }

  virtual ~ABlockOp() {}

  int getNumRows() const {
    return m;
  }

  int getNumCols() const {
    return n;
  }

  void setReference ( const int I, const int J, OpType &op ) {
    blockMat[ I * n + J] = &op;
  }

  OpType* getPointer ( const int I, const int J ) const {
    return blockMat[ I * n + J ];
  }

  OpType& getReference ( const int I, const int J ) const {
    if ( !blockMat[ I * n + J ] ) {
      cerr<<"\n blockMat["<<I<<"*"<<n<<"+"<<J<<"] is not set \n";
    }
    return *blockMat[ I * n + J ];
  }

protected:
  const int m, n;
  vector< OpType* >  blockMat;
};

/** Obviously this operator represents a composition of several operators
 *  to be successively applied.
 *
 * \author Droske
 */
template <class DomRanType>
class CompositeOp : public Op<DomRanType> {
public:
  CompositeOp() {}

  virtual ~CompositeOp() {}

  void appendReference ( const Op<DomRanType> &Op ) {
    ops.push_back ( &Op );
  }

  void reset() {
    ops.erase ( ops.begin(), ops.end() );
  }

  virtual void applyAdd ( const DomRanType &Arg, DomRanType &Dest ) const {
    DomRanType tmp ( Dest );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  virtual void apply ( const DomRanType &Arg, DomRanType &Dest ) const {

    // \todo it is assumed, that sizes of all vectors are constant.

    DomRanType tmp ( Dest );

    Dest = Arg;
    for ( typename vector<const Op<DomRanType>* >::const_iterator it = ops.begin(); it != ops.end(); ++it ) {
      tmp = Dest;
      ( *it )->apply ( tmp, Dest );
    }
  }
protected:
  vector<const Op<DomRanType>* > ops;
};


/** What do think this is?
 * \author Droske
 */
template <typename DomRanType>
class IdentityOp : public Op<DomRanType> {
public:
  IdentityOp() {}
  
  template<typename DummyType>
  IdentityOp( const DummyType& /*dummy*/ ) {}
  
  void applyAdd ( const DomRanType& Arg, DomRanType &Dest ) const {
    Dest += Arg;
  }

  void apply ( const DomRanType& Arg, DomRanType &Dest ) const {
    Dest = Arg;
  }

  void applySingle ( DomRanType& /*Arg*/ ) const {
  }

  void setZero() {
    throw aol::Exception( "Trying to call setZero() on the IdentityOp, which makes no sense:", __FILE__, __LINE__);
  }

};


/** Maps everything onto the zero vector.
 * \author Droske
 */
template <typename DomType, typename RangeType = DomType>
class NullOp : public Op<DomType, RangeType> {
public:
  void applyAdd ( const DomType& , RangeType& ) const {
  }

  void applySingle ( DomType& arg ) const {
    arg.setZero();
  }
};



template <typename DomType, typename RangeType = DomType>
class ConstOp : public Op<DomType, RangeType> {
public:

  ConstOp ( const RangeType &Vec ) : _vec ( Vec ) { }

  void applyAdd ( const DomType&, RangeType &Dest ) const {
    Dest += _vec;
  }

protected:
  const RangeType &_vec;
};

/**
 * Creates an Op from MultiVector to Scalar, by applying a given Op from Vector to Scalar
 * on every component of the argument MultiVector, adding up the resulting Scalars and
 * storing this in the destination Scalar.
 *
 * \author Berkels
 */
template <typename RealType, typename OpType = Op<Vector<RealType>, Scalar<RealType> > >
class ScalarVecToScalarMVecOp : public Op<MultiVector<RealType>, Scalar<RealType> > {
  const OpType &_vecOp;
public:
  ScalarVecToScalarMVecOp ( const OpType &VecOp )
    : _vecOp ( VecOp ) {}

  void applyAdd ( const MultiVector<RealType> &MArg, Scalar<RealType> &Dest ) const {
    for ( int i = 0; i < MArg.numComponents(); ++i )
      _vecOp.applyAdd ( MArg[i], Dest );
  }

  void applyDerivative ( const MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    for ( int i = 0; i < MArg.numComponents(); ++i )
      _vecOp.applyDerivative ( MArg[i], MDest[i] );
  }

  template<typename BlockMatrixType>
  void applyAddMultipleSecondDerivative ( const MultiVector<RealType> &MArg, BlockMatrixType &MatDest, const RealType Factor ) const {
    for ( int i = 0; i < MArg.numComponents(); ++i )
      _vecOp.applyAddMultipleSecondDerivative ( MArg[i], MatDest.getReference( i, i ), Factor );
  }
};

/**
 * Creates an Op from MultiVector to MultiVector, by applying a given Op from Vector to Vector
 * on every component of the argument MultiVector and storing the resulting Vectors in the
 * corresponding components of the destination MultiVector.
 *
 * \author Berkels
 */
template <typename RealType>
class VecToMVecOp : public Op<MultiVector<RealType> > {
  const Op<Vector<RealType> > &_vecOp;
public:
  VecToMVecOp ( const Op<Vector<RealType> > &VecOp )
    : _vecOp ( VecOp ) {}

  void applyAdd ( const MultiVector<RealType> &MArg, MultiVector<RealType> &MDest ) const {
    for ( int i = 0; i < MArg.numComponents(); ++i )
      _vecOp.applyAdd ( MArg[i], MDest[i] );
  }

  void apply ( const MultiVector<RealType> &MArg, MultiVector<RealType> &MDest ) const {
    for ( int i = 0; i < MArg.numComponents(); ++i )
      _vecOp.apply ( MArg[i], MDest[i] );
  }
};

/**
 * Like ScalarVecToScalarMVecOp, but the application of the underlying Op can be scaled
 * differently for each component of the argument MultiVector.
 *
 * \author Berkels
 */
template <typename RealType, int NumComponents>
class ScaledScalarVecToScalarMVecOp : public Op<MultiVector<RealType>, Scalar<RealType> > {
  const Op<Vector<RealType>, Scalar<RealType> > &_vecOp;
  // Intentionally make a copy of the scaling vec.
  const aol::Vec<NumComponents, RealType> _scaling;
public:
  ScaledScalarVecToScalarMVecOp ( const Op<Vector<RealType>, Scalar<RealType> > &VecOp, const aol::Vec<NumComponents, RealType> &Scaling )
    : _vecOp ( VecOp ),
      _scaling ( Scaling ) {}

  void applyAdd ( const MultiVector<RealType> &MArg, Scalar<RealType> &Dest ) const {
    for ( int i = 0; i < NumComponents; ++i ) {
      Dest /= _scaling[i];
      _vecOp.applyAdd ( MArg[i], Dest );
      Dest *= _scaling[i];
    }
  }
};

/**
 * Like VecToMVecOp, but the application of the underlying Op can be scaled
 * differently for each component of the argument MultiVector.
 *
 * \author Berkels
 */
template <typename RealType, int NumComponents>
class ScaledVecToMVecOp : public Op<MultiVector<RealType> > {
  const Op<Vector<RealType> > &_vecOp;
  // Intentionally make a copy of the scaling vec.
  const aol::Vec<NumComponents, RealType> _scaling;
public:
  ScaledVecToMVecOp ( const Op<Vector<RealType> > &VecOp, const aol::Vec<NumComponents, RealType> &Scaling )
    : _vecOp ( VecOp ),
      _scaling ( Scaling ) {}

  void applyAdd ( const MultiVector<RealType> &MArg, MultiVector<RealType> &MDest ) const {
    for ( int i = 0; i < NumComponents; ++i ) {
      if ( _scaling[i] != aol::ZOTrait<RealType>::one )
        MDest[i] /= _scaling[i];
      _vecOp.applyAdd ( MArg[i], MDest[i] );
      if ( _scaling[i] != aol::ZOTrait<RealType>::one )
        MDest[i] *= _scaling[i];
    }
  }

  void apply ( const MultiVector<RealType> &MArg, MultiVector<RealType> &MDest ) const {
    for ( int i = 0; i < NumComponents; ++i ) {
      _vecOp.apply ( MArg[i], MDest[i] );
      if ( _scaling[i] != aol::ZOTrait<RealType>::one )
        MDest[i] *= _scaling[i];
    }
  }
};

/**
 * Creates an Op from Vector to Scalar by applying a given Op from MultiVector to Scalar
 * on the argument Vector interpreted as MultiVector with one component.
 *
 * \author Berkels
 */
template <typename RealType>
class ScalarMVecToScalarVecOp : public Op<Vector<RealType>, Scalar<RealType> > {
  const Op<MultiVector<RealType>, Scalar<RealType> > &_scalarMVecOp;
public:
  ScalarMVecToScalarVecOp ( const Op<MultiVector<RealType>, Scalar<RealType> > &ScalarMVecOp )
    : _scalarMVecOp ( ScalarMVecOp ) {}

  void applyAdd ( const Vector<RealType> &Arg, Scalar<RealType> &Dest ) const {
    MultiVector<RealType> mArg;
    mArg.appendReference ( Arg );
    _scalarMVecOp.applyAdd ( mArg, Dest );
  }

  void apply ( const Vector<RealType> &Arg, Scalar<RealType> &Dest ) const {
    MultiVector<RealType> mArg;
    mArg.appendReference ( Arg );
    _scalarMVecOp.apply ( mArg, Dest );
  }
};

/**
 * Creates an Op from Vector to Vector by applying a given Op from MultiVector to MultiVector
 * on the argument and destination Vector both interpreted as MultiVector with one component.
 *
 * \author Berkels
 */
template <typename RealType>
class MVecToVecOp : public Op<Vector<RealType> > {
  const Op<MultiVector<RealType> > &_mVecOp;
public:
  MVecToVecOp ( const Op<MultiVector<RealType> > &MVecOp )
    : _mVecOp ( MVecOp ) {}

  void applyAdd ( const Vector<RealType> &Arg, Vector<RealType> &Dest ) const {
    MultiVector<RealType> mArg, mDest;
    mArg.appendReference ( Arg );
    mDest.appendReference ( Dest );
    _mVecOp.applyAdd ( mArg, mDest );
  }

  void apply ( const Vector<RealType> &Arg, Vector<RealType> &Dest ) const {
    MultiVector<RealType> mArg, mDest;
    mArg.appendReference ( Arg );
    mDest.appendReference ( Dest );
    _mVecOp.apply ( mArg, mDest );
  }
};

/**
 * Creates an Op from MultiVector to Scalar, by applying a given Op from Vector to Scalar
 * on the selected component of the argument MultiVector.
 *
 * \author Berkels
 */
template <typename RealType>
class ScalarVecToScalarComponentMVecOp : public Op<MultiVector<RealType>, Scalar<RealType> > {
  const Op<Vector<RealType>, Scalar<RealType> > &_vecOp;
  const int _componentNum;
public:
  ScalarVecToScalarComponentMVecOp ( const Op<Vector<RealType>, Scalar<RealType> > &VecOp, const int ComponentNum )
  : _vecOp ( VecOp ), _componentNum ( ComponentNum ) {}

  void applyAdd ( const MultiVector<RealType> &MArg, Scalar<RealType> &Dest ) const {
    _vecOp.applyAdd ( MArg[_componentNum], Dest );
  }
};

/**
 * Creates an Op from MultiVector to MultiVector, by applying a given Op from Vector to Vector
 * on the selected component of the argument MultiVector and storing the resulting Vectors in the
 * corresponding component of the destination MultiVector.
 *
 * \author Berkels
 */
template <typename RealType>
class VecToComponentMVecOp : public Op<MultiVector<RealType> > {
  const Op<Vector<RealType> > &_vecOp;
  const int _componentNum;
public:
  VecToComponentMVecOp ( const Op<Vector<RealType> > &VecOp, const int ComponentNum )
  : _vecOp ( VecOp ), _componentNum ( ComponentNum ) {}

  void applyAdd ( const MultiVector<RealType> &MArg, MultiVector<RealType> &MDest ) const {
    _vecOp.applyAdd ( MArg[_componentNum], MDest[_componentNum] );
  }

  void apply ( const MultiVector<RealType> &MArg, MultiVector<RealType> &MDest ) const {
    _vecOp.apply ( MArg[_componentNum], MDest[_componentNum] );
  }
};

/**
 * Creates an Op from a class that has the member function applyDerivative, by calling applyDerivative in apply.
 *
 * \author Berkels
 */
template <typename RealType, typename WrappedOpType, typename VecType = aol::Vector<RealType>, typename DerivativeType = VecType>
class DerivativeWrapper : public aol::Op<VecType, DerivativeType> {
  const WrappedOpType &_op;
public:
  DerivativeWrapper ( const WrappedOpType &Op ) : _op ( Op ) {}

  void applyAdd ( const VecType &Arg, DerivativeType &Dest ) const {
    DerivativeType result ( Dest, aol::STRUCT_COPY );
    apply ( Arg, result );
    Dest += result;
  }

  void apply ( const VecType &Arg, DerivativeType &Dest ) const {
    _op.applyDerivative ( Arg, Dest );
  }
};

/**
 * Creates an Op from a class that has the member function applySecondDerivative, by calling applySecondDerivative in apply.
 *
 * \author Wirth
 */
template <typename RealType, typename WrappedOpType, typename VecType, typename DestType>
class SecondDerivativeWrapper : public aol::Op<VecType,DestType> {
  const WrappedOpType &_op;
public:
  SecondDerivativeWrapper ( const WrappedOpType &Op ) : _op ( Op ) {}

  void applyAdd ( const VecType &Arg, DestType &Dest ) const {
    DestType result ( Dest, aol::STRUCT_COPY );
    apply ( Arg, result );
    Dest += result;
  }

  void apply ( const VecType &Arg, DestType &Dest ) const {
    _op.applySecondDerivative ( Arg, Dest );
  }
};

/**
 * Scalar valued function that maps everything to zero.
 *
 * \author Berkels
 */
template <typename VectorType>
class NullEnergy : public aol::NullOp<VectorType, aol::Scalar<typename VectorType::RealType> > {
  typedef typename VectorType::RealType RealType;
public:
  template <typename InitType>
  explicit NullEnergy ( const InitType & ) { }

  void applyAddDerivative ( const VectorType &/*Arg*/, VectorType &/*Dest*/ ) const {
  }

  template <typename MatrixType>
  void applyAddSecondDerivative ( const VectorType &/*Arg*/, MatrixType &/*Dest*/ ) const {
  }
};

/** Generates some noise. By default, use platform-independent reproducible sequence of random numbers. Call randomize() to randomize.
 * \author Droske
 *
 * \warning NoiseOperator::apply/applySingle ( MultiVector ) is unpredictable with OpenMP because every thread
 *          needs to use the member randGen. Hence, these functions should not be used when OpenMP is used!
 */
template <typename RealType>
class NoiseOperator : public BiOp<Vector<RealType> > {

public:
  enum NOISE_TYPE {
    SALT_AND_PEPPER,
    GAUSSIAN,
    EQUALLY_DISTRIBUTED
  };

protected:
  NOISE_TYPE noiseType;
  RealType min, max;
  double ratio;
  mutable RandomGenerator randGen;

public:

  /*  @param  min  For NoiseType = GAUSSIAN gives the expectation value of the normal distribution generating the noise
   *  @param  max  For NoiseType = GAUSSIAN gives the standard deviation of the normal distribution generating the noise
   */
  NoiseOperator ( NOISE_TYPE NoiseType = SALT_AND_PEPPER,
                  RealType Min = ZOTrait<RealType>::zero,
                  RealType Max = ZOTrait<RealType>::one ) :
    noiseType ( NoiseType ), min ( Min ), max ( Max ),  ratio ( 0.25 ), randGen () {
  }

  void setNoiseRatio ( double Ratio ) {
    ratio = Ratio;
  }

  //! randomize random generator
  void randomize ( ) {
    randGen.randomize();
  }

  RandomGenerator& getRandomGeneratorRef ( ) {
    return randGen;
  }

  void applyAdd ( const Vector<RealType> &Arg, Vector<RealType> &Dest ) const {
    Vector<RealType> tmp ( Dest.size() );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  void apply ( const Vector<RealType>& Arg, Vector<RealType> &Dest ) const {
    if ( Arg.size() != Dest.size() )
      throw aol::Exception ("aol::NoiseOperator::apply: incompatible vector lengths", __FILE__, __LINE__ );

    const int size = Dest.size();

    switch ( noiseType ) {
    case SALT_AND_PEPPER:
      for ( int i = 0; i < size; ++i ) {
        if ( randGen.rReal<RealType> ( 0, 1 ) <= ratio ) {
          Dest[i] = ( randGen.rBool() ? max : min );
        } else {
          Dest[ i ] = Arg.get ( i );
        }
      }
      break;

    case GAUSSIAN:
      for ( int i = 0; i < size; ++i ) {
        Dest[i] = randGen.normalrReal<RealType> ( min, max ); // mean, stdDev
      }
      break;

    case EQUALLY_DISTRIBUTED:
      for ( int i = 0; i < size; ++i ) {
        Dest[i] = randGen.rReal<RealType> ( min, max );
      }
      break;
    
    default:
        throw aol::UnimplementedCodeException ( "NoiseOperator::apply: unknown NoiseType", __FILE__, __LINE__ );
    }
  }

  //! Fill Vector with noise
  void applySingle ( Vector<RealType> &Dest ) const {
    Vector<RealType> dummy ( Dest, STRUCT_COPY );
    this->apply ( dummy, Dest );
  }

  //! Add noise to Vector
  void applyAddSingle ( Vector<RealType> &Dest ) const {
    Vector<RealType> dummy ( Dest, STRUCT_COPY );
    this->applySingle ( dummy );
    Dest += dummy;
  }

  using aol::BiOp< Vector<RealType> >::apply;
  using aol::BiOp< Vector<RealType> >::applySingle;
  using aol::BiOp< Vector<RealType> >::applyAdd;
  using aol::BiOp< Vector<RealType> >::applyAddSingle;
};

/** A Nonlinear Operator, which depends on vector Vec
 * \author diepenbruck,olischlaeger
 */
template <typename DomainType, typename RangeType = DomainType>
class DependOp : public Op< DomainType, RangeType> {
public:

  DependOp() : Vec ( NULL ) { }

  DependOp ( DomainType &vec ) : Vec ( &vec ) { }

  virtual ~DependOp() { }

  virtual void setV ( DomainType &V ) {
    Vec = &V;
  }

  const DomainType &getV() const {
    if ( !Vec ) {
      cerr << " in aol::DependOp: setV first !\n";
      throw aol::Exception ( "in aol::DependOp: setV first !", __FILE__ , __LINE__ );
    }
    return *Vec;
  }

  void apply ( const DomainType &Arg, RangeType &Dest ) const {
    if ( !Vec ) {
      cerr << " in aol::DependOp: setV first !\n";
      throw aol::Exception ( "in aol::DependOp: setV first !", __FILE__ , __LINE__ );
    }
    Dest.clear();
    applyAdd ( Arg, Dest );
  }
protected:
  DomainType *Vec;
};


//! Compares two operators and checks whether they produce the same result up to a given tolerance
template <class RealType>
bool compareOps( const aol::Op<aol::Vector<RealType> > &Op1, const aol::Op<aol::Vector<RealType> > &Op2, int ArgSize, int DestSize, RealType tolerance = 0.0 ) {
  aol::Vector<RealType> Arg ( ArgSize );
  aol::Vector<RealType> Dest1 ( DestSize );
  aol::Vector<RealType> Dest2 ( DestSize );
  aol::Vector<RealType> Diff  ( DestSize );

  RealType col1Norm = 0.;     // 1-Norm of the column-vectors
  RealType max1Norm = 0.;     // 1-operator-Norm of the Difference of the Matrices

  bool success = true;

  Arg.setZero();

#ifdef VERBOSE
  ProgressBar<> pb ( "Comparing Ops" );
  pb.start ( ArgSize );
#endif

  for ( int j = 0; j < ArgSize; ++j ) {
#ifdef VERBOSE
    pb++;
#endif
    Arg[ j ] = 1.;
    if ( j > 0 ) Arg[ j - 1 ] = 0.;

    Op1.apply ( Arg, Dest1 );
    Op2.apply ( Arg, Dest2 );

    Diff = Dest1;
    Diff -= Dest2;

    if ( Diff.norm() > tolerance ) {
      col1Norm = Diff.lpNorm(1);
      if ( col1Norm > max1Norm ) max1Norm = col1Norm;
      success = false;
      cerr << "inconsistency in column " << j << " found, column-l1-norm: " << col1Norm << "\n";
    }
  }
#ifdef VERBOSE
  pb.finish();
#endif

  if ( !success ) cerr << endl << color::red << "1-Operator-Norm of Difference: " << max1Norm << endl << color::reset;

  return success;
}

//! 1-Operatornorm of difference
template <class RealType>
RealType distanceOfOps ( const aol::Op<aol::Vector<RealType> > &Op1, const aol::Op<aol::Vector<RealType> > &Op2, int ArgSize, int DestSize ) {
  aol::Vector<RealType> Arg ( ArgSize );
  aol::Vector<RealType> Dest1 ( DestSize );
  aol::Vector<RealType> Dest2 ( DestSize );

  RealType result = 0;

  Arg.setZero();
  for ( int j = 0; j < ArgSize; ++j ) {
    Arg[ j ] = 1.;
    if ( j > 0 ) Arg[ j - 1 ] = 0.;

    Op1.apply ( Arg, Dest1 );
    Op2.apply ( Arg, Dest2 );

    Dest1 -= Dest2;
    result = aol::Max<RealType> ( result, Dest1.lpNorm ( 1 ) );
  }

  return result;
}

/** This operator is a composition of several operators where
 *  dimensions of range and image of the individual ops may be
 *  different (composition of non-endomorphisms).
 *  CompositeOp cannot provide this because temporary vectors of
 *  appropriate dimension need to be stored.  These must be given
 *  explicitely because some operators cannot provide dimension of
 *  image (e. g. IdentityOp).
 * \author Schwen
 */
template <class DomRanType>
class Composite_differentDim_Op : public Op<DomRanType> {
public:
  Composite_differentDim_Op() {}

  virtual ~Composite_differentDim_Op() {}

  //! @param dst_dim dimension of the range of this operator. If first operator appended with -1, all dimensions are assumed to be those of the destination vector.
  void appendReference ( Op<DomRanType> &Op, int dst_dim = -1 ) {
    ops.push_back ( &Op );
    dims.push_back ( dst_dim );
  }

  void reset() {
    ops.erase ( ops.begin(), ops.end() );
    dims.erase ( dims.begin(), dims.end() );
  }

  virtual void apply ( const DomRanType &Arg, DomRanType &Dest ) const {
    if ( dims[0] == -1 ) { // all dimensions are the same
      DomRanType tmp ( Dest.size() );

      Dest = Arg;
      for ( typename vector<const Op<DomRanType>* >::const_iterator it = ops.begin(); it != ops.end(); ++it ) {
        tmp = Dest;
        ( *it )->apply ( tmp, Dest );
      };

    } else {             // jump between dimensions
      DomRanType tmp_new ( Arg.size() ), tmp_old ( 0 );
      tmp_new = Arg;

      typename vector<int>::const_iterator dit = dims.begin();
      for ( typename vector<const Op<DomRanType>* >::const_iterator it = ops.begin(); it != ops.end(); ++it, ++dit ) {
        tmp_old.resize ( tmp_new.size() );
        tmp_old = tmp_new;
        tmp_new.resize ( ( *dit ) );

        ( *it )->apply ( tmp_old, tmp_new );
      };

      Dest = tmp_new;
    };
  }

  virtual void applyAdd ( const DomRanType &Arg, DomRanType &Dest ) const {
    DomRanType tmp ( Dest );
    apply ( Arg, tmp );
    Dest += tmp;
  }

protected:
  vector<const Op<DomRanType>* > ops;
  vector<int> dims;
};

template <typename DomainType, typename DataType = typename DomainType::DataType>
class GenEnergyOp : public aol::Op<DomainType, Scalar<DataType> > {
public:
  GenEnergyOp() {}
  virtual ~GenEnergyOp() {}
  virtual void getLastEnergy ( Vector<DataType> &Energy ) const = 0;
};

/**
 * \author Berkels
 */
template <typename DomainType, typename DataType = typename DomainType::DataType>
class StandardGenEnergyOp : public aol::GenEnergyOp<DomainType, DataType> {
protected:
  mutable DataType _lastEnergy;
public:
  StandardGenEnergyOp()
    : _lastEnergy ( aol::MaxInitializerTrait<DataType>::MaxInitializer ) {}
  virtual ~StandardGenEnergyOp() {}

  void getLastEnergy( aol::Vector<DataType> &Energy ) const {
    Energy.resize(1);
    Energy[0] = _lastEnergy;
  }
};

/** The MultiVector variant of RestrOp (subset of DOFs, not multigrid)
 */
template <class DataType>
class MultiRestrOp : public Op<MultiVector<DataType> > {
public:
  //! @param op             operator to be restricted, NB: read-only, hence not changed!
  //! @param selectionMask  bitfield, 1 in selected region, 0 outside
  MultiRestrOp ( const Op<MultiVector<DataType> > &op,
     const BitVector &selectionMask )
      : _op ( op ), _selectionMask ( selectionMask ) { }

  //! deprecated, but extensively used by the surfMesh module. Better use BitVector.
  MultiRestrOp ( const Op<MultiVector<DataType> > &op, const vector<bool> &selectionMask )
      : _op ( op ), _selectionMask ( selectionMask.size() ) {
    for( unsigned int i = 0; i < selectionMask.size(); ++i ){
      _selectionMask.set( i, selectionMask[i] );
    };
  }

  virtual ~MultiRestrOp() { }

  virtual ostream& print ( ostream& out ) const {
    cerr << "Instance of RestrOp:\n";
    cerr << "selectionMask: \n";
    for ( int i = 0; i < _selectionMask.size(); ++i ) {
      if ( _selectionMask[i] ) {
        cerr << i << ": true" << endl;
      } else {
        cerr << i << ":  --" << endl;
      }
    }
    _op.print ( out );
    return out;
  }

  virtual void apply ( const MultiVector<DataType> &Arg, MultiVector<DataType> &Dest ) const {
  int comp, i;

  if ( ! ( Arg.compareDim ( Dest ) ) ||
       _selectionMask.size() !=  Arg[0].size()  ||
       _selectionMask.size() !=  Dest[0].size() ||
       Arg.numComponents()   != Dest.numComponents() ) {
    throw aol::Exception ( "Size of arg, dest and selection incompatible.", __FILE__, __LINE__ );
  }

  aol::MultiVector< DataType > tmp ( Arg.numComponents(), _selectionMask.size() );

  // Simulate deletion of column elements of _op outside selection:
  for ( comp = 0; comp < Arg.numComponents(); ++comp ) {
    for ( i = 0; i < Arg[0].size(); ++i ) {
      if ( _selectionMask[i] ) {
        tmp[comp][i] = Arg[comp][i];
      }  // Inside selection.
      else                    {
        tmp[comp][i] = 0.0;
      }           // Outside selection.
    }
  }

  // Apply the simulated operator:
  _op.apply ( tmp, Dest );

  // Simulate the replacement of the diagonal elements of _op by 1
  // AND the deletion of row elements outside selection:
  for ( comp = 0; comp < Arg.numComponents(); ++comp ) {
    for ( i = 0; i < Arg[0].size(); ++i ) {
      if ( ! ( _selectionMask[i] ) ) {
        Dest[comp][i] = Arg[comp][i];   // Outside selection.
      }
    }
  }
}

  virtual void applyAdd ( const MultiVector<DataType> &Arg, MultiVector<DataType> &Dest ) const {
    MultiVector<DataType> tmp ( Dest );
    apply ( Arg, tmp );
    Dest += tmp;
  }

private:
  const Op<MultiVector<DataType> >  &_op;
  aol::BitVector                _selectionMask;
};



/** Class for simulating the restriction of an operator to be only
  * applied on a subset of the DOFs. This is not restriction in the
  * multigrid sense.
  *
  * E.g. let op be a matrix. Then the method apply works as if the
  * rows and columns of op corresponding to the outside of the
  * selected region have been deleted and the corresponding diagonal
  * element of op has been replaced by 1.
  *
  * NB: A read-only reference to op is being saved, so op itself is
  * not changed!
*/
template <class DataType>
class RestrOp : public Op<Vector<DataType> > {
public:
  //! @param op             operator to be restricted, NB: read-only, hence not changed!
  //! @param selectionMask  bitfield, 1 in selected region, 0 outside
  RestrOp ( const Op<Vector<DataType> > &op, const aol::BitVector &selectionMask )
      : _op ( op ),
      _blockOp ( op ),
      _multiOp ( _blockOp, selectionMask ) {}

  //! deprecated, but extensively used by the surfMesh module; better use BitVector.
  RestrOp ( const Op<Vector<DataType> > &op, const vector<bool> &selectionMask )
      : _op ( op ),
      _blockOp ( op ),
      _multiOp ( _blockOp, selectionMask ) {
  }

  virtual ~RestrOp() { }

  virtual void applyAdd ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const {
    MultiVector<DataType> MArg ( 0, Arg.size() );
    MultiVector<DataType> MDest ( 0, Dest.size() );

    MArg.appendReference ( Arg );
    MDest.appendReference ( Dest );

    _multiOp.applyAdd ( MArg, MDest );
  }

private:
  const Op<Vector<DataType> >  &_op;
  DiagonalBlockOp<DataType> _blockOp;
  MultiRestrOp<DataType> _multiOp;
};

template< class T >
ostream& operator<< ( ostream& out, const RestrOp<T>& op ) {
  return op.print ( out );
}

/**
 * This operator gets two vectors Vec1 and Vec2 and computes
 * \f[ Dest = (Vec1 Vec2^T ) Arg = (Vec2 \cdot Arg) Vec1 \f]
 *
 * \author Teusner
 */
template < typename RealType, typename VectorType >
class TensorOp
  : public aol::Op<VectorType, VectorType> {
protected:
  const VectorType &_vec1;
  const VectorType &_vec2;
public:
  TensorOp ( const VectorType &Vec1,
             const VectorType &Vec2 ) :
    _vec1 ( Vec1 ),
    _vec2 ( Vec2 ) {}

  TensorOp ( const VectorType &Vec ) :
    _vec1 ( Vec ),
    _vec2 ( Vec ) {}

  virtual ~TensorOp () {}

  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    Dest = _vec1;
    Dest *= (_vec2 * Arg);
  }

  virtual void applyAdd ( const VectorType &Arg, VectorType &Dest ) const {
    VectorType tmp (Arg, aol::STRUCT_COPY);
    this->apply(Arg, tmp);
    Dest += tmp;
  }
};

/**
 * Given an Op  \f$ A \f$ from VecType to VecType, creates an Op from VecType to aol::Scalar with the mapping  \f$ x \mapsto \frac{1}{2}A[x]\cdot x \f$.
 *
 * \author Berkels
 */
template <typename VecType>
class QuadraticFormOp : public aol::StandardGenEnergyOp<VecType> {
protected:
  const aol::Op<VecType, VecType> &_vectorValuedOp;
public:
  QuadraticFormOp ( const aol::Op<VecType, VecType> &VectorValuedOp )
    : _vectorValuedOp ( VectorValuedOp ) {}
  virtual ~QuadraticFormOp () {}
  virtual void apply ( const VecType &Arg, aol::Scalar<typename VecType::DataType> &Dest ) const {
    VecType temp ( Arg, aol::STRUCT_COPY );
    _vectorValuedOp.apply ( Arg, temp );
    this->_lastEnergy = 0.5 * ( temp * Arg );
    Dest[0] = this->_lastEnergy;
  }
  virtual void applyAdd ( const VecType &Arg, aol::Scalar<typename VecType::DataType> &Dest ) const {
    VecType temp ( Arg, aol::STRUCT_COPY );
    _vectorValuedOp.apply ( Arg, temp );
    this->_lastEnergy = 0.5 * ( temp * Arg );
    Dest[0] += this->_lastEnergy;
  }
};


/**
 * Given an Op  \f$ A \f$ from aol::Vector to aol::Vector, creates an Op from aol::Vector to aol::Scalar with the mapping  \f$ x \mapsto \frac{1}{2}(A[x]\cdot A[x])\f$.
 *
 * \author Berkels
 */
template <typename VecType>
class SquaredNormOp : public aol::StandardGenEnergyOp<VecType> {
  typedef typename VectorInitTrait<VecType>::InitType RangeSizeType;
  typedef typename VecType::RealType RealType;
protected:
  const aol::Op<VecType, VecType> &_vectorValuedOp;
  const RangeSizeType _opDestVecLength;
public:
  SquaredNormOp ( const aol::Op<VecType, VecType> &VectorValuedOp, const RangeSizeType &OpDestVecLength )
    : _vectorValuedOp ( VectorValuedOp ), _opDestVecLength ( OpDestVecLength ) {}

  virtual ~SquaredNormOp () {}

  RealType computeEnergy ( const VecType &Arg ) const {
    VecType temp ( _opDestVecLength );
    _vectorValuedOp.apply ( Arg, temp );
    this->_lastEnergy = 0.5 * temp.normSqr();
    return this->_lastEnergy;
  }

  virtual void apply ( const VecType &Arg, aol::Scalar<RealType> &Dest ) const {
    Dest[0] = computeEnergy ( Arg );
  }

  virtual void applyAdd ( const VecType &Arg, aol::Scalar<RealType> &Dest ) const {
    Dest[0] += computeEnergy ( Arg );;
  }
};


/**
 * Gradient of SquaredNormOp. Useful for derivative tests of vector valued operators used in Gauss-Newton approaches.
 *
 * \author Berkels
 */
template <typename VecType, typename JacobianMatrixType>
class SquaredNormOpGradient : public aol::Op<VecType> {
  typedef typename VectorInitTrait<VecType>::InitType RangeSizeType;
protected:
  const aol::Op<VecType> &_F;
  const aol::Op<VecType, JacobianMatrixType > &_DF;
  const RangeSizeType _dimRangeF;
  mutable VecType _FAtPos;
  mutable JacobianMatrixType _DFAtPos;
public:
  SquaredNormOpGradient ( const aol::Op<VecType> &F,
                          const aol::Op<VecType, JacobianMatrixType > &DF,
                          const RangeSizeType &DimRangeF,
                          const int DFNumRows,
                          const int DFNumCols )
    : _F ( F ), _DF ( DF ), _dimRangeF ( DimRangeF ), _FAtPos ( DimRangeF ), _DFAtPos ( DFNumRows, DFNumCols ) {}

  virtual ~SquaredNormOpGradient () {}

  void applyAdd ( const VecType &Arg, VecType &Dest ) const {
    _F.apply ( Arg, _FAtPos );
    _DF.apply ( Arg, _DFAtPos );
    _DFAtPos.applyAddTransposed ( _FAtPos, Dest );
  }
};

}   // end namespace aol

#endif

