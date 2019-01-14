#ifndef __SOLVER_H
#define __SOLVER_H

#include <op.h>
#include <vec.h>
#include <matrix.h>
#include <QRDecomposition.h>
#include <GaussSeidel.h>
#include <solverInfo.h>
#include <pointerClasses.h>
#include <vectorExtensions.h>

namespace aol {

/** General basis class for inverse operators ("solvers")
 */
template <typename VectorType >
class InverseOp : public aol::Op< VectorType > {};

/** General abstract basis class for iterative solvers
 */
template < typename VectorType, typename OpType = Op<VectorType> >
class IterativeInverseOp : public InverseOp<VectorType> {
public:
  typedef typename VectorType::DataType DataType;

protected:

  const OpType&               _op;        //!< reference to the operator to be inverted
  mutable DeleteFlagPointer<SolverInfo<DataType> > _infoPtr;

public:
  IterativeInverseOp ( const OpType &Op,
                       const DataType Epsilon,
                       const int MaxIter,
                       const StoppingMode stop,
                       const bool Quiet,
                       ostream& Out )
      : _op ( Op )
      , _infoPtr ( new SolverInfo<DataType> ( Epsilon, MaxIter, stop, Quiet, Out ), true ) {}

  IterativeInverseOp ( const OpType &Op,
                       SolverInfo<DataType> & info )
      : _op ( Op )
      , _infoPtr ( &info, false ) {}

  virtual ~IterativeInverseOp () {}

  // ------------------------------------------------------------------------
  //! Set solver accuracy
  void setAccuracy ( const DataType accuracy ) {
    _infoPtr->setAccuracy ( accuracy );
  }
  // ------------------------------------------------------------------------
  DataType getAccuracy ( ) const {
    return _infoPtr->getAccuracy();
  }
  // ------------------------------------------------------------------------
  void setMaxIterations ( int MaxIter ) {
    _infoPtr->setMaxIterations ( MaxIter );
  }
  // ------------------------------------------------------------------------
  //! Select the type of stopping criterion
  void setStopping ( const StoppingMode mode ) {
    _infoPtr->setStoppingMode ( mode );
  }
  // ------------------------------------------------------------------------
  StoppingMode getStopping ( ) const {
    return _infoPtr->getStoppingMode();
  }
  // ------------------------------------------------------------------------
  //! Set whether the solver should print intermediate results
  void setQuietMode ( const bool Quiet = true ) {
    _infoPtr->setQuietMode ( Quiet );
  }
  // ------------------------------------------------------------------------
  //! Nothing is printed!
  void setMegaQuietMode ( const bool Quiet = true ) {
    _infoPtr->setMegaQuietMode ( Quiet );
  }
  // ------------------------------------------------------------------------
  bool getQuietMode() const {
    return _infoPtr->getQuietMode();
  }
  // ------------------------------------------------------------------------
  bool getMegaQuietMode() const {
    return _infoPtr->getMegaQuietMode();
  }
  // ------------------------------------------------------------------------
  //! Return reference to operator being inverted
  const Op<VectorType>& getOpReference ( ) const {
    return ( _op );
  }
  // ------------------------------------------------------------------------
  SolverInfo<DataType> & getSolverInfoReference() {
    return *_infoPtr;
  }
  const SolverInfo<DataType> & getSolverInfoReference() const {
    return *_infoPtr;
  }
  // ------------------------------------------------------------------------
  void setSolverInfoPointer ( SolverInfo<DataType> * info ) {
    _infoPtr = DeleteFlagPointer<SolverInfo<DataType> > ( info, false );
  }
  // ------------------------------------------------------------------------
  void setSolverInfoReference ( SolverInfo<DataType> & info ) {
    _infoPtr = DeleteFlagPointer<SolverInfo<DataType> > ( *info, false );
  }
  // ------------------------------------------------------------------------

  //! Return the number of iterations used for the last run of this solver
  int getCount() const {
    return _infoPtr->getIterationCount();
  }

  //! Apply solver adding result to Dest, possibly using Dest as initial guess
  virtual void applyAdd ( const VectorType &Arg, VectorType &Dest ) const {
    VectorType result ( Dest, aol::DEEP_COPY ); // so that initial guess can be used
    apply ( Arg, result );
    Dest += result;
  }

  //! Apply solver, possibly using Dest as initial guess
  virtual void apply ( const VectorType &/*Arg*/, VectorType &/*Dest*/ ) const = 0;

  int getMaxIterations ( ) const {
    return ( _infoPtr->getMaxIterations() );
  }

protected:
  ostream& getOstream ( ) const {
    return ( _infoPtr->getOstream() );
  }

};



/**
 * \brief Class for approximation of inverse for cg solver.
 * \ingroup iterativeSolver
 * \author Droske
 */
template < typename VectorType, typename OpType = Op<VectorType> >
class CGInverse : public IterativeInverseOp<VectorType, OpType> {
  typedef typename IterativeInverseOp<VectorType, OpType>::DataType DataType;

public:
  CGInverse ( const OpType       &Op,
              const DataType     Epsilon = 1e-16,
              const int          MaxIter = 1000,
              const StoppingMode Stop = STOPPING_UNSET,
              ostream&           Out = cerr )
      : IterativeInverseOp<VectorType, OpType> ( Op, Epsilon, MaxIter, Stop, false, Out ) {}

  CGInverse ( const OpType       &Op,
              SolverInfo<DataType> & info )
      : IterativeInverseOp<VectorType, OpType> ( Op, info ) {}

  virtual ~CGInverse () {}

  //! Return the squared current residuum
  DataType getResSqr ( ) const {
    return this->_infoPtr->getFinalResidual();
  }

  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    DataType spa = 0, spn, q, quad;

    VectorType r ( Arg, aol::STRUCT_COPY );
    VectorType p ( Arg, aol::STRUCT_COPY );
    VectorType h ( Arg, aol::STRUCT_COPY );

    this->_op.apply ( Dest, h );

    //! \todo think about whether it makes sense to write a method for
    //!       aol::{Multi}Vector to set vector to linear combination
    //!       of h and arg to avoid double loop
    r = h;
    r -= Arg;

    p = Arg;
    p -= h;

    spn = r * r;

    this->_infoPtr->startIterations ( Arg.normSqr(), spn, "cg", "l_2 norm ^2" );

    while ( ! ( this->_infoPtr->stoppingCriterionIsFulfilled() ) && ! ( this->_infoPtr->maxIterIsReached() ) && ! ( this->_infoPtr->currentResidualIsNaN() ) ) {
      this->_infoPtr->startStep();

      // case starting with second iteration
      if ( this->_infoPtr->getIterationCount() > 1 ) {
        const DataType e = spn / spa;
        p *= e;
        p -= r;
      }

      // basic iteration step
      this->_op.apply ( p, h );

      quad = p * h;
      q    = spn / quad;

      Dest.addMultiple ( p, q );
      r.addMultiple ( h, q );

      spa = spn;

      // compute new residuum
      spn = r * r;

      this->_infoPtr->finishStep ( spn );
    }

    this->_op.apply ( Dest, h );
    r = h;
    r -= Arg;
    spn = r * r;

    this->_infoPtr->finishIterations ( spn );
  }

  // end class CGInverse
};


//! Preconditioned Conjugate Gradient solver
/*!
 * Class for preconditioned conjugate gradient solver.
 * Specify approximate inverse of the Operator \b Op by \b ApproxInverseOp, e. g. DiagonalPreconditioner
 * \ingroup iterativeSolver
 * \author Droske
 */
template < typename VectorType, typename OpType = Op<VectorType>, typename iOpType = Op<VectorType> >
class PCGInverse : public IterativeInverseOp<VectorType, OpType> {
  typedef typename IterativeInverseOp<VectorType, OpType>::DataType DataType;

protected:
  const iOpType &_approxInverseOp;

public:
  PCGInverse ( const OpType &Op,
               const iOpType &ApproxInverseOp,
               const DataType Epsilon = 1e-16,
               const int MaxIter = 50,
               const StoppingMode Stop = STOPPING_UNSET,
               ostream &Out = cerr )
      : IterativeInverseOp< VectorType, OpType > ( Op, Epsilon, MaxIter, Stop, false, Out )
      , _approxInverseOp ( ApproxInverseOp ) {}

  PCGInverse ( const OpType &Op,
               const iOpType &ApproxInverseOp,
               SolverInfo<DataType> & info )
      : IterativeInverseOp< VectorType, OpType > ( Op, info )
      , _approxInverseOp ( ApproxInverseOp ) {}

  virtual ~PCGInverse () {}

  DataType getResSqr ( ) const {
    return this->_infoPtr->getFinalResidual();
  }

  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    DataType alpha_numer, alpha_denom, beta_numer, beta_denom, spn;

    VectorType g ( Arg, aol::STRUCT_COPY );
    VectorType d ( Arg, aol::STRUCT_COPY );
    VectorType h ( Arg, aol::STRUCT_COPY ); // automatically cleared

    this->_op.apply ( Dest, h );

#ifdef VERBOSE
    this->_infoPtr->getOstream() << "Solver h norm " << h.norm() << endl;
#endif

    g = h;
    g -= Arg;

    h.setZero();
    _approxInverseOp.apply ( g, h );

#ifdef VERBOSE
    this->_infoPtr->getOstream() << "h*h " << h*h << endl;
#endif

    spn = g * g;

#ifdef VERBOSE
    this->_infoPtr->getOstream() << " spn = " << spn << endl;
#endif
    d -= h;

#ifdef VERBOSE
    this->_infoPtr->getOstream() << "h*h " << h*h << endl;
#endif

    this->_infoPtr->startIterations ( Arg.normSqr(), spn, "p-cg", "l_2 norm ^2" );

    while ( ! ( this->_infoPtr->stoppingCriterionIsFulfilled() ) && ! ( this->_infoPtr->maxIterIsReached() ) && ! ( this->_infoPtr->currentResidualIsNaN() ) ) {
      this->_infoPtr->startStep();

      beta_denom = alpha_numer = g * h;

#ifdef VERBOSE
      this->_infoPtr->getOstream() << "beta_denom = " << beta_denom << endl
      << "alpha_numer = " << alpha_numer << endl;
#endif

      h.setZero();
      this->_op.apply ( d, h );
      alpha_denom = d * h;

#ifdef VERBOSE
      this->_infoPtr->getOstream() << "h * h = " << h*h << endl
      << "alpha_denom = " << alpha_denom << endl;
#endif

      Dest.addMultiple ( d, alpha_numer / alpha_denom );

      g.addMultiple ( h, alpha_numer / alpha_denom );

#ifdef VERBOSE
      this->_infoPtr->getOstream() << "g * g = " << g*g << endl;
#endif

      h.setZero();
      _approxInverseOp.apply ( g, h );
      beta_numer = g * h;

#ifdef VERBOSE
      this->_infoPtr->getOstream() << "h * h = " << h*h << endl
      << "beta_numer = " << beta_numer << endl;
#endif

      d *= ( beta_numer / beta_denom );
      d -= h;

#ifdef VERBOSE
      this->_infoPtr->getOstream() << "d * d = " << d*d << endl;
#endif

      spn = g * g;
#ifdef VERBOSE
      VectorType dummy ( Arg, aol::STRUCT_COPY );
      this->_op.apply ( Dest, dummy );
      dummy -= Arg;
      this->_infoPtr->getOstream() << "Computed residuum = " << dummy * dummy << endl;
#endif
      this->_infoPtr->finishStep ( spn );
    }

    this->_op.apply ( Dest, h );
    g = h;
    g -= Arg;
    spn = g * g;

    this->_infoPtr->finishIterations ( spn );
  }

  // end class PCGInverse
};


//! BiCG for nonsymmetric matrices.
/**
 * Class for approximation of inverse for bicg solver working also for
 * non-symmetric operators, but the transpose of the operator has to be known.
 * \ingroup iterativeSolver
 * \author Droske
 */
template < typename VectorType, typename OpType = Op<VectorType> >
class BiCGInverse : public IterativeInverseOp<VectorType, OpType> {
  typedef typename IterativeInverseOp<VectorType, OpType>::DataType DataType;

protected:
  const OpType &_transposeOp;

public:
  BiCGInverse ( const OpType &Op,
                const OpType &TransposeOp,
                const DataType Epsilon = 1e-16,
                const int MaxIter = 50,
                ostream &Out = cerr )
      : IterativeInverseOp< VectorType, OpType > ( Op, Epsilon, MaxIter, STOPPING_UNSET, false, Out )
      , _transposeOp ( TransposeOp ) {}

  BiCGInverse ( const OpType &Op,
                const OpType &TransposeOp,
                SolverInfo<DataType> & info )
      : IterativeInverseOp< VectorType, OpType > ( Op, info )
      , _transposeOp ( TransposeOp ) {}

  virtual ~BiCGInverse () {}

  virtual void apply ( const VectorType &MArg, VectorType &MDest ) const {

    VectorType r  ( MArg, aol::STRUCT_COPY );
    VectorType rt ( MArg, aol::STRUCT_COPY );
    VectorType p  ( MArg, aol::STRUCT_COPY );
    VectorType pt ( MArg, aol::STRUCT_COPY );
    VectorType q  ( MArg, aol::STRUCT_COPY );
    VectorType qt ( MArg, aol::STRUCT_COPY );

    this->_op.apply ( MDest, q );

    r = MArg;
    r -= q;

    rt = r;

    DataType rho_new, rho_old = 0.0, alpha, beta;

    this->_infoPtr->startIterations ( MArg.normSqr(), r * r, "bicg", "l_2 norm ^2" );

    while ( ! ( this->_infoPtr->stoppingCriterionIsFulfilled() ) && ! ( this->_infoPtr->maxIterIsReached() ) && ! ( this->_infoPtr->currentResidualIsNaN() ) ) {
      this->_infoPtr->startStep();

      rho_new = r * rt;

      if ( rho_new == 0. ) {
        throw Exception ( "bicg fails because rho = 0.", __FILE__, __LINE__ );
      }

      if ( this->_infoPtr->getIterationCount() == 1 ) {
        p = r;
        pt = rt;
      } else {
        beta = rho_new / rho_old;

        p  *= beta;
        p  += r;
        pt *= beta;
        pt += rt;
      }
      this->_op.apply ( p, q );

      _transposeOp.apply ( pt, qt );

      alpha = rho_new / ( pt * q );

      MDest.addMultiple ( p, alpha );

      r.addMultiple ( q, -alpha );
      rt.addMultiple ( qt, -alpha );

      DataType resSqr =  r * r;

      rho_old = rho_new;
      this->_infoPtr->finishStep ( resSqr );
    }
    this->_infoPtr->finishIterations();
  }

  // end of class BiCGInverse
};


//! Preconditioned BiCG.
/**
 * Class for approximation of inverse for preconditioned bicg solver working also for
 * non-symmetric operators, but the transpose of the operator has to be known.
 * \author Droske
 * \ingroup iterativeSolver
 */
template < typename VectorType, typename OpType = Op<VectorType>, typename iOpType = Op<VectorType> >
class PBiCGInverse : public IterativeInverseOp<VectorType, OpType> {
  typedef typename IterativeInverseOp<VectorType, OpType>::DataType DataType;

protected:
  const OpType &_transposeOp;
  const iOpType &_approxInverseOp;
  const iOpType &_approxInverseTransposeOp;

public:
  PBiCGInverse ( const OpType &Op,
                 const OpType &TransposeOp,
                 const iOpType &ApproxInverseOp,
                 const iOpType &ApproxInverseTransposeOp,
                 const DataType Epsilon = 1e-16,
                 const int MaxIter = 50,
                 ostream& Out = cerr )
      : IterativeInverseOp< VectorType, OpType > ( Op, Epsilon, MaxIter, STOPPING_UNSET, false, Out )
      , _transposeOp ( TransposeOp ), _approxInverseOp ( ApproxInverseOp )
      , _approxInverseTransposeOp ( ApproxInverseTransposeOp ) {}

  PBiCGInverse ( const OpType &Op,
                 const OpType &TransposeOp,
                 const iOpType &ApproxInverseOp,
                 const iOpType &ApproxInverseTransposeOp,
                 SolverInfo<DataType> & info )
      : IterativeInverseOp< VectorType, OpType > ( Op, info )
      , _transposeOp ( TransposeOp ), _approxInverseOp ( ApproxInverseOp )
      , _approxInverseTransposeOp ( ApproxInverseTransposeOp ) {}

  virtual ~PBiCGInverse () {}

  virtual void apply ( const VectorType &MArg, VectorType &MDest ) const {
    VectorType r  ( MArg, aol::STRUCT_COPY );
    VectorType rt ( MArg, aol::STRUCT_COPY );
    VectorType z  ( MArg, aol::STRUCT_COPY );
    VectorType zt ( MArg, aol::STRUCT_COPY );
    VectorType p  ( MArg, aol::STRUCT_COPY );
    VectorType pt ( MArg, aol::STRUCT_COPY );
    VectorType q  ( MArg, aol::STRUCT_COPY );
    VectorType qt ( MArg, aol::STRUCT_COPY );

    this->_op.apply ( MDest, q );

    r = MArg;
    r -= q;

    rt = r;

    DataType rho_new, rho_old = 1, alpha, beta;

    this->_infoPtr->startIterations ( MArg.normSqr(), r * r, "p-bicg", "l_2 norm ^2" );

    while ( ! ( this->_infoPtr->stoppingCriterionIsFulfilled() ) && ! ( this->_infoPtr->maxIterIsReached() ) && ! ( this->_infoPtr->currentResidualIsNaN() ) ) {
      this->_infoPtr->startStep();

      _approxInverseOp.apply ( r, z );

      _approxInverseTransposeOp.apply ( rt, zt );

      rho_new = z * rt;

      if ( rho_new == 0. ) {
        throw Exception ( "p-bicg fails because rho = 0.", __FILE__, __LINE__ );
      }

      if ( this->_infoPtr->getIterationCount() == 1 ) {
        p = z;
        pt = zt;
      } else {
        beta = rho_new / rho_old;

        p  *= beta;
        p  += z;
        pt *= beta;
        pt += zt;
      }

      this->_op.apply ( p, q );

      _transposeOp.apply ( pt, qt );

      alpha = rho_new / ( pt * q );

      MDest.addMultiple ( p, alpha );

      r.addMultiple ( q, -alpha );
      rt.addMultiple ( qt, -alpha );

      rho_old = rho_new;

      DataType resSqr =  r * r;

      this->_infoPtr->finishStep ( resSqr );
    }

    this->_op.apply ( MDest, z );
    z -= MArg;
    this->_infoPtr->finishIterations ( z.normSqr() );
  }
  // end of class PBiCGInverse
};


//! BiCG-stab solver with preconditioning.
/**
 * Class for approximation of inverse for preconditioned bicg-stab solver working also for
 * non-symmetric operators.
 * \ingroup iterativeSolver
 * \author Droske
 */
template < typename VectorType, typename OpType = Op<VectorType>, typename iOpType = Op<VectorType> >
class PBiCGStabInverse : public IterativeInverseOp<VectorType, OpType> {
  typedef typename IterativeInverseOp<VectorType, OpType>::DataType DataType;

protected:
  const iOpType &_approxInverseOp;

public:
  PBiCGStabInverse ( const OpType &Op,
                     const iOpType &ApproxInverseOp,
                     DataType Epsilon = 1e-16,
                     const int MaxIter = 50,
                     const StoppingMode Stop = STOPPING_UNSET,
                     ostream &Out = cerr )
      : IterativeInverseOp< VectorType, OpType > ( Op, Epsilon, MaxIter, Stop, false, Out )
      , _approxInverseOp ( ApproxInverseOp ) {}

  PBiCGStabInverse ( const OpType &Op,
                     const iOpType &ApproxInverseOp,
                     SolverInfo<DataType> & info )
      : IterativeInverseOp< VectorType, OpType > ( Op, info )
      , _approxInverseOp ( ApproxInverseOp ) {}

  virtual ~PBiCGStabInverse () {}

  virtual void apply ( const VectorType &MArg, VectorType &MDest ) const {

    DataType residSqr;
    DataType rho1, rho2 = 0, alpha = 0, beta, omega = 0;

    VectorType p      ( MArg,  aol::STRUCT_COPY );
    VectorType phat   ( MDest, aol::STRUCT_COPY ); // This is important if arg and dest space have same dimension but not same structure
    VectorType s      ( MArg,  aol::STRUCT_COPY ); // e.g. block-ops with non-symmetric block sizes
    VectorType shat   ( MDest, aol::STRUCT_COPY ); // ATTN: In this case the preconditioner has the transposed structure of the op!
    VectorType t      ( MArg,  aol::STRUCT_COPY );
    VectorType v      ( MArg,  aol::STRUCT_COPY );
    VectorType r      ( MArg,  aol::STRUCT_COPY );
    VectorType rtilde ( MArg,  aol::STRUCT_COPY );

    this->_op.apply ( MDest, rtilde );
    r = MArg;
    r -= rtilde;
    rtilde = r;
    residSqr = r.normSqr();

    this->_infoPtr->startIterations ( MArg.normSqr(), residSqr, "p-bicg-stab", "l_2 norm ^2" );

    while ( ! ( this->_infoPtr->maxIterIsReached() ) && ! ( this->_infoPtr->currentResidualIsNaN() ) ) {
      this->_infoPtr->startStep();

      rho1 = rtilde * r;
      if ( rho1 == 0. )
        break;

      if ( this->_infoPtr->getIterationCount() == 1 ) {
        p = r;
      } else {
        beta = ( rho1 / rho2 ) * ( alpha / omega );
        // p = r + beta * (p - omega * v);
        p.addMultiple ( v, -omega );
        p *= beta;
        p += r;
      }

      _approxInverseOp.apply ( p, phat );
      this->_op.apply ( phat, v );

      alpha = rho1 / ( rtilde * v );
      // s = r - alpha(0) * v;
      s = r;
      s.addMultiple ( v, -alpha );
      residSqr = s.normSqr ();
      this->_infoPtr->setCurrentResidual ( residSqr );
      if ( this->_infoPtr->stoppingCriterionIsFulfilled() ) {
        MDest.addMultiple ( phat, alpha );
        this->_infoPtr->finishStep ( residSqr );
        break;
      }

      _approxInverseOp.apply ( s, shat );
      this->_op.apply ( shat, t );
      omega = ( t * s ) / ( t * t );
      MDest.addMultiple ( phat, alpha );
      MDest.addMultiple ( shat, omega );
      r = s;
      r.addMultiple ( t, -omega );

      rho2 = rho1;

      residSqr = r.normSqr ();

      this->_infoPtr->finishStep ( residSqr );
      if ( omega == 0 )
        break;
    }

    this->_op.apply ( MDest, r );
    r -= MArg;
    this->_infoPtr->finishIterations ( r.normSqr() );
  }
  // end of class PBiCGStab
};


//! computes the inverse by GMRES(m) algorithm. works on non-symmetric matrices, but is slower.
/*!
 * Generalized Mean Residuals.
 * Class for approximation of inverse of an operator, which
 * also works for non-SPD-operators.
 * The class operates by applying a GMRES(m) algorithm.
 * TODO: implement relativeStoppingCriterion + print convergence warning.
 * \ingroup iterativeSolver
 * \author Droske
 */
template < typename VectorType, typename OpType = aol::Op<VectorType> >
class GMRESInverse : public IterativeInverseOp< VectorType, OpType > {
  typedef typename IterativeInverseOp< VectorType, OpType >::DataType DataType;

protected:
  const int _maxInnerIter;

public:
  // constructor
  GMRESInverse ( const OpType &Op,
                 const DataType Epsilon = 1e-16,
                 const int MaxInnerIter = 10,
                 const int MaxIter = 50,
                 ostream &Out = cerr )
      : IterativeInverseOp< VectorType, OpType > ( Op, Epsilon, MaxIter, aol::STOPPING_ABSOLUTE, false, Out )
      , _maxInnerIter ( MaxInnerIter ) {}

  GMRESInverse ( const OpType &Op,
                 SolverInfo<DataType> & info,
                 const int MaxInnerIter = 10 )
      : IterativeInverseOp< VectorType, OpType > ( Op, info )
      , _maxInnerIter ( MaxInnerIter ) {}

  virtual ~GMRESInverse () {}

  virtual void apply ( const VectorType &MArg, VectorType &MDest ) const {

    if ( this->getStopping() != aol::STOPPING_ABSOLUTE )
      throw aol::UnimplementedCodeException ( "aol::GMRESInverse::apply not implemented for stopping mode other than STOPPING_ABSOLUTE.", __FILE__, __LINE__ );

    DataType beta;

    VectorType z ( MArg, aol::STRUCT_COPY );
    std::vector<VectorType *> v;

    for ( int i = 0; i < _maxInnerIter + 1; ++i ) {
      v.push_back ( new VectorType ( MArg, aol::STRUCT_COPY ) );
    }

    aol::Vector<DataType> c ( _maxInnerIter + 1 );
    FullMatrix<DataType> h ( _maxInnerIter + 1, _maxInnerIter + 1 );

    // z = Ax - b
    this->_op.apply ( MDest, z );
    z -= MArg;

    this->_infoPtr->startIterations ( MArg.normSqr(), z.normSqr(), "GMRES", "l_2 norm ^2" );

    while ( ! ( this->_infoPtr->maxIterIsReached() ) && ! ( this->_infoPtr->currentResidualIsNaN() ) ) {
      this->_infoPtr->startStep();

      h.setZero();
      c.setZero();

      beta = sqrt ( z * z );
      *v[ 0 ] = z;
      *v[ 0 ] *= 1. / beta;

      for ( int k = 0; k < this->_maxInnerIter; ++k ) {
        z.setZero();
        this->_op.apply ( *v[ k ], z );

        for ( int i = 0; i <= k; ++i ) {
          h.set ( i, k, ( *v[i] ) *z );
        }

        *v[ k+1 ] = z;
        for ( int i = 0; i <= k; ++i ) {
          v[ k+1 ]->addMultiple ( *v[ i ], -h.get ( i, k ) );
        }

        h.set ( k + 1, k, sqrt ( ( *v[k+1] ) * ( *v[k+1] ) ) );
        *v[ k+1 ] *= 1. / h.get ( k + 1, k );

        // do the QR-decomposition now
        FullMatrix<DataType> H ( k + 2, k + 1 );
        for ( int i = 0; i < k + 2; ++i ) {
          for ( int j = 0; j < k + 1; ++j ) {
            H.set ( i, j, h.get ( i, j ) );
          }
        }

        QRDecomposeHouseholderMatrix<DataType> decomp;
        FullMatrix<DataType> Qmat ( k + 2, k + 2 );
        FullMatrix<DataType> R ( k + 2, k + 1 );

        decomp.transform ( H, R, Qmat );
        Qmat.transpose();

        aol::Vector<DataType> d ( k + 2 );

        this->_infoPtr->setCurrentResidual ( aol::Sqr ( beta * Qmat.get ( k + 1, 0 ) ) );
        this->_infoPtr->printStats();

        for ( int i = 0; i < k + 2; ++i )
          d[ i ] = Qmat.get ( i, 0 );

        d *= -beta;
        backSolve ( R, c, d );

        if ( ( k == this->_maxInnerIter - 1 ) || ( this->_infoPtr->stoppingCriterionIsFulfilled() ) )  {
          // compute solution vector.
          for ( int i = 0; i < /*maxInnerIter = */ k + 1; ++i )
            MDest.addMultiple ( *v[ i ], c[ i ] );
          break;
        }
      } // end inner loop

      // z = Ax - b, exact residual used at the beginning of the outer loop.
      this->_op.apply ( MDest, z );
      z -= MArg;

      // finishStep before finishIterations, because of estimated <~> exact residual check in finishIterations.
      // Estimated residual: If inner loop finished, the algorithm should be finished.
      this->_infoPtr->finishStep ( this->_infoPtr->getCurrentResidual () );

      if ( this->_infoPtr->stoppingCriterionIsFulfilled() ) {
        break;
      }

    } // end outer loop

    // delete the allocated basis vectors, before finishing...
    for ( typename vector<VectorType *>::iterator it = v.begin(); it != v.end(); ++it )
      delete *it;

    // calculate residual once again
    this->_op.apply ( MDest, z );
    z -= MArg;
    this->_infoPtr->finishIterations ( z.normSqr() );
  }

protected:

  void backSolve ( const FullMatrix<DataType> &R, Vector<DataType> &Solution, const Vector<DataType> &RHS ) const {
    if ( Solution.size() < R.getNumCols() || RHS.size() < R.getNumCols() ) {
      throw Exception ( "GMRESInverse::backSolve: incompatible dimensions.\n", __FILE__, __LINE__ );
    }
    const int n = R.getNumCols();
    for ( int i = n - 1; i >= 0; --i ) {
      DataType s = RHS[ i ];
      for ( int j = i + 1; j < n; ++j ) {
        s -= Solution[ j ] * R.get ( i, j );
      }
      Solution[ i ] = s / R.get ( i, i );
    }
  }
  // end of class GMRESInverse
};

/*!
 * \brief Preconditioned GMRES(m) algorithm
 * \details This class implements a preconditioned version of the GMRES(m) algorithm (GMRES with restarts).
 *          The algorithm is very similiar to the one given in "Numerik linearer Gleichungssysteme: Direkte und iterative Verfahren, C. Kanzow, 2005, Springer",
 *          especially in the naming of the variables. However, the termination criterion of the inner loop is different to imitate the behavior of the older aol::GMRESInverse.
 *          Also the H in this algorithm is called \f$\bar{H}\f$ in the book.
 * \remark One iteration as displayed by aol::SolverInfo corresponds to one complete preconditioned GMRES run with at most maxInnerIter iterations!
 * \author Toelkes
 */
template < typename VectorType, typename OpType = aol::Op < VectorType >, typename PrecondType = aol::Op < VectorType > >
class PGMRESInverse : public IterativeInverseOp < VectorType, OpType > {
  typedef typename IterativeInverseOp < VectorType, OpType >::DataType DataType;

  //! The preconditioner
  const PrecondType &_precond;

  //! Maximum number of inner iterations (=GMRES iterations)
  int _maxIter;

  //! \brief Calculates some variables needed at the start of a GMRES run.
  inline void calculateInitialVariables ( const VectorType &arg, const VectorType &dest, VectorType &temp, VectorType &q, VectorType &v0, DataType &z0 ) const {
    // t = b - Ax^0 (arg = b, dest = x^0)
    this->_op.apply ( dest, temp );
    temp *= -1.0;
    temp += arg;

    // Pq = b - Ax^0 (P preconditioner)
    _precond.apply ( temp, q );

    // z^0 = beta = ||q^0||
    z0 = q.norm ();

    // v^0 = q^0/beta
    v0 = q;
    v0 /= z0;
  }

public:
  /*!
   * \brief Constructor
   * \param[in] op The operator of the linear system
   * \param[in] precond The preconditioner
   * \param[in] epsilon The parameter for the stoppin criterion
   * \param[in] stoppingMode Stopping mode to use
   * \param[in] maxInnerIter Maximum number of iterations in the GMRES algorithm
   * \param[in] maxRestarts Maximum number of restarts, i.e. maximum number of GMRES runs
   * \param[in] out Stream to print output to
   */
  PGMRESInverse ( const OpType &op,
      const PrecondType &precond,
      const DataType epsilon = 1e-8,
      const aol::StoppingMode stoppingMode = aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM,
      const int maxInnerIter = 25,
      const int maxRestarts = 50,
      ostream &out = std::cerr )
   // One iteration as counted by solverInfo is one complete run of the preconditioned GMRES method with at most maxIter iterations,
   // ~> use maxRestarts as bound for the maximum number of iterations.
   : IterativeInverseOp < VectorType, OpType > ( op, epsilon, maxRestarts, stoppingMode, false, out ),
     _precond ( precond ), _maxIter ( maxInnerIter ) {}

  /*!
   * \brief Constructor using the aol::SolverInfo class
   * \param[in] op The operator of the linear system
   * \param[in] precond The preconditioner
   * \param[in] info SolverInfo class
   * \param[in] maxInnerIter Maximum number of iterations in the GMRES algorithm
   */
  PGMRESInverse ( const OpType &op,
      const PrecondType &precond,
      SolverInfo < DataType > &info,
      const int maxInnerIter = 25 )
  : IterativeInverseOp < VectorType, OpType > ( op, info ),
    _precond ( precond ), _maxIter ( maxInnerIter ) {}

  virtual ~PGMRESInverse () {}

  virtual void apply ( const VectorType &arg, VectorType &dest ) const {
    FullMatrix < DataType > h ( _maxIter + 1, _maxIter + 1 );

    VectorType t ( arg, aol::STRUCT_COPY );
    VectorType q ( arg, aol::STRUCT_COPY );
    VectorType w ( arg, aol::STRUCT_COPY );
    aol::Vector < DataType > c ( this->_maxIter );
    aol::Vector < DataType > s ( this->_maxIter );
    aol::Vector < DataType > z ( this->_maxIter + 1 );
    aol::Vector < DataType > y ( this->_maxIter );

    aol::RandomAccessContainer < VectorType* > v ( this->_maxIter + 1 );
    for ( int i = 0; i < v.size (); ++i ) {
      v[i] = new VectorType ( arg, aol::STRUCT_COPY );
    }

    DataType tau;
    DataType nu;

    this->calculateInitialVariables ( arg, dest, t, q, *v[0], z[0] );

    // z[0] is set to ||q^0|| in calculateInitialVariables
    DataType beta = z[0];

    this->_infoPtr->startIterations ( arg.norm (), beta, "preconditioned GMRES", "l^2 norm" );

    while ( !( this->_infoPtr->maxIterIsReached () ) && !( this->_infoPtr->currentResidualIsNaN () ) ) {
      this->_infoPtr->startStep ();

      // inner iteration counter
      int k;
      for ( k = 0; k < this->_maxIter; ++k ) {
        // modified Arnoldi algorithm
        // Pw^k = Av^k
        this->_op.apply ( *v[k], t );
        _precond.apply ( t, w );

        for ( int i = 0; i <= k; ++i ) {
          h.set ( i, k, *v[i] * w );
          w.addMultiple ( *v[i], -h.get ( i, k ) );
        }

        h.set ( k + 1, k, w.norm () );
        *v[k + 1] = w;
        *v[k + 1] /= h.get ( k + 1, k );

        // apply old Givens rotation on the k-th column of H_k
        for ( int i = 0; i < k; ++i ) {
          DataType newik = c[i] * h.get ( i, k ) + s[i] * h.get ( i + 1, k );
          DataType newip1k = -s[i] * h.get ( i, k ) + c[i] * h.get ( i + 1, k );

          h.set ( i, k, newik );
          h.set ( i + 1, k, newip1k );
        }

        // new Givens rotation for the elimination of H_k(k+1, k)
        tau = aol::Abs ( h.get ( k, k ) ) + aol::Abs ( h.get ( k+1, k ) );
        nu = tau * sqrt ( aol::Sqr ( h.get ( k, k ) / tau ) + aol::Sqr ( h.get ( k + 1, k ) / tau ) );
        c[k] = h.get ( k, k ) / nu;
        s[k] = h.get ( k + 1, k ) / nu;

        // apply new Givens rotation to H_k
        h.set ( k, k, nu );
        h.set ( k + 1, k, 0.0 );

        // apply new Givens rotation to the RHS
        z[k + 1] = -s[k] * z[k];
        z[k] = c[k] * z[k];

        this->_infoPtr->setCurrentResidual ( aol::Abs ( z[k + 1] ) );
        this->_infoPtr->printStats ();

        // Check for termination criterion
        if ( this->_infoPtr->stoppingCriterionIsFulfilled () )
          break;
      }

      // If k = maxIter, decrement k by 1 to work with the last set entry
      if ( k == this->_maxIter )
        --k;

      // back solve
      y[k] = z[k] / h.get ( k, k );

      for ( int i = k - 1; i >= 0; --i ) {
        y[i] = z[i];
        for ( int j = i + 1; j <= k; ++j ) {
          y[i] -= h.get ( i, j ) * y[j];
        }
        y[i] /= h.get ( i, i );
      }

      // Calculate approximate solution
      // dest = x[0] already
      for ( int i = 0; i <= k; ++i ) {
        dest.addMultiple ( *v[i], y[i] );
      }

      this->_infoPtr->finishStep ( aol::Abs ( z[k + 1] ) );

      // Check for termination criterion
      if ( this->_infoPtr->stoppingCriterionIsFulfilled () )
        break;

      // Reinitialize variables for the next GMRES run
      this->calculateInitialVariables ( arg, dest, t, q, *v[0], z[0] );
    }

    // Recalculate the (exact) residuum
    this->_op.apply ( dest, t );
    t -= arg;

    // Give the exact residuum to finishIterations, where it will be compared to the estimated residuum.
    this->_infoPtr->finishIterations ( t.norm () );

    for ( int i = 0; i < v.size (); ++i )
      delete v[i];
  }
};


/** \brief Jacobi solver for operators that support makeRowEntries
 *  \ingroup iterativeSolver
 *  \author Schwen (Droske)
 */
template < typename VectorType, typename OpType >
class JacobiInverse : public IterativeInverseOp<VectorType, OpType> {
  typedef typename IterativeInverseOp<VectorType, OpType>::DataType DataType;

protected:
  DataType _relax;

public:
  JacobiInverse ( const OpType &Op,
                  const DataType Epsilon = 1e-16,
                  const int MaxIter = 1000,
                  const DataType Relax = 1.0,
                  const StoppingMode Stop = STOPPING_UNSET,
                  ostream& Out = cerr )
      : IterativeInverseOp<VectorType, OpType> ( Op, Epsilon, MaxIter, Stop, false, Out )
      , _relax ( Relax ) {}

  JacobiInverse ( const OpType &Op,
                  SolverInfo<DataType> & info,
                  const DataType Relax = 1.0 )
      : IterativeInverseOp<VectorType, OpType> ( Op, info )
      , _relax ( Relax ) {}

  virtual ~JacobiInverse () {}

  DataType getResSqr ( ) const {
    return this->_infoPtr->getFinalResidual();
  }

  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    std::vector<typename Row<DataType>::RowEntry > vec;

    Vector<DataType> dummy ( Dest, aol::STRUCT_COPY );

    this->_op.apply ( Dest, dummy );
    dummy -= Arg;
    this->_infoPtr->startIterations ( Arg.normSqr(), dummy.normSqr(), "Jacobi", "l_2 norm ^2" );

    while ( !this->_infoPtr->stoppingCriterionIsFulfilled() && ! ( this->_infoPtr->maxIterIsReached() ) && ! ( this->_infoPtr->currentResidualIsNaN() ) ) {
      this->_infoPtr->startStep();

      this->_op.apply ( Dest, dummy );

      DataType res = 0.0;

      for ( int i = 0; i < Dest.size(); ++i ) {
        this->_op.makeRowEntries ( vec, i );

        for ( typename vector<typename Row<DataType>::RowEntry >::iterator it = vec.begin(); it != vec.end(); ++it ) {
          if ( it->col == i ) {                                             // is this efficient??
            Dest[i] += _relax * ( Arg[i] - dummy[i] ) / it->value;
            break;
          }
        }
        res += aol::Sqr ( Arg[i] - dummy[i] );
      }
      this->_infoPtr->finishStep ( res );
    }
    this->_infoPtr->finishIterations();

  }
  // end of class JacobiInverse
};

/** \brief Gauss-Seidel solver for operators that support makeRowEntries. Default behavior is symmetric Gauss-Seidel due to historical reasons.
 *  \ingroup iterativeSolver
 *  \author Schwen (Droske)
 */
template < typename VectorType, typename OpType >
class GaussSeidelInverse : public IterativeInverseOp<VectorType, OpType> {
  typedef typename IterativeInverseOp<VectorType, OpType>::DataType DataType;

protected:
  DataType _relax;
  aol::GaussSeidelSweepingMode _gss;
  aol::GaussSeidelSweeper<DataType, VectorType, OpType> _sweeper;

public:
  GaussSeidelInverse ( const OpType &Op,
                       const DataType Epsilon = 1e-16,
                       const int MaxIter = 1000,
                       const DataType Relax = 1.0,
                       const aol::GaussSeidelSweepingMode Gss = aol::GAUSS_SEIDEL_SYMMETRIC,
                       const StoppingMode Stop = STOPPING_UNSET,
                       ostream &Out = cerr )
      : IterativeInverseOp<VectorType, OpType> ( Op, Epsilon, MaxIter, Stop, false, Out )
      , _relax ( Relax )
      , _gss ( Gss )
      , _sweeper ( Op, Relax ) {}

  GaussSeidelInverse ( const OpType &Op,
                       SolverInfo<DataType> & info,
                       const DataType Relax = 1.0,
                       const aol::GaussSeidelSweepingMode Gss = aol::GAUSS_SEIDEL_SYMMETRIC )
      : IterativeInverseOp<VectorType, OpType> ( Op, info )
      , _relax ( Relax )
      , _gss ( Gss )
      , _sweeper ( Op, Relax ) {}

  virtual ~GaussSeidelInverse () {}

public:
  void setGaussSeidelSweepingMode ( const aol::GaussSeidelSweepingMode Gss ) {
    _gss = Gss;
  }

  DataType getResSqr () const  {
    return this->_infoPtr->getFinalResidual();
  }

  void setRelax ( const DataType Relax ) {
    _relax = Relax;
  }

  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    VectorType dummy ( Dest, aol::STRUCT_COPY );

    this->_op.apply ( Dest, dummy );
    dummy -= Arg;
    this->_infoPtr->startIterations ( Arg.normSqr(), dummy.normSqr(), "Gauss-Seidel", "l_2 norm ^2" );

    while ( !this->_infoPtr->stoppingCriterionIsFulfilled() && ! ( this->_infoPtr->maxIterIsReached() ) && ! ( this->_infoPtr->currentResidualIsNaN() ) ) {
      this->_infoPtr->startStep();

      switch ( _gss ) {
        case aol::GAUSS_SEIDEL_FORWARD:
          _sweeper.apply ( Arg, Dest, aol::GAUSS_SEIDEL_FORWARD );
          break;
        case aol::GAUSS_SEIDEL_SYMMETRIC:
          _sweeper.apply ( Arg, Dest, ( this->_infoPtr->getIterationCount() & 2 ? aol::GAUSS_SEIDEL_BACKWARD : aol::GAUSS_SEIDEL_FORWARD ) );
          break;
        case aol::GAUSS_SEIDEL_RED_BLACK:
          _sweeper.apply ( Arg, Dest, aol::GAUSS_SEIDEL_EVEN_ONLY );
          _sweeper.apply ( Arg, Dest, aol::GAUSS_SEIDEL_ODD_ONLY );
          break;
        default:
          throw aol::Exception ( "aol::GaussSeidelInverse::apply: No Gauss-Seidel sweeping mode selected", __FILE__, __LINE__ );
          break;
      };

      // compute new residuum
      this->_op.apply ( Dest, dummy );
      dummy -= Arg;
      this->_infoPtr->finishStep ( dummy.normSqr() );
    }
    this->_infoPtr->finishIterations();
  }

private:
  // this class has pointer members, so the following should be implemented:

  GaussSeidelInverse ( ) {
    throw aol::UnimplementedCodeException ( "aol::GaussSeidelInverse standard constructor not implemented", __FILE__, __LINE__ );
  }

  GaussSeidelInverse ( const GaussSeidelInverse<VectorType, OpType>& ) {
    throw aol::UnimplementedCodeException ( "aol::GaussSeidelInverse copy constructor not implemented", __FILE__, __LINE__ );
  }

  GaussSeidelInverse<VectorType, OpType>& operator= ( const GaussSeidelInverse<VectorType, OpType>& ) {
    throw aol::UnimplementedCodeException ( "aol::GaussSeidelInverse::operator= not implemented", __FILE__, __LINE__ );
  }
  // end of class GaussSeidelInverse
};


/** \brief Simple interface for using direct inversion methods of matrix classes as solvers
 *  \ingroup directSolver
 *  \author Geihe
 */
template < typename VectorType, typename MatType >
class DirectSolver
: public MatType, public InverseOp<VectorType>
{
public:
  typedef typename VectorType::DataType DataType;

  DirectSolver ()
  : MatType ( )
  , InverseOp<VectorType> ( )
  {}

  explicit DirectSolver ( const MatType & mat )
  : MatType ( mat.inverse() )
  , InverseOp<VectorType> ( )
  {}

  void setInverseMatrix( const MatType & mat )  {
    MatType::operator= ( mat );
  }

  virtual ~DirectSolver () {}

  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    MatType::apply( Arg, Dest );
  }

  virtual void applyAdd ( const VectorType &Arg, VectorType &Dest ) const {
    MatType::applyAdd( Arg, Dest );
  }
};


/**
 * Solves (R^tR)X=RHS, if R is an upper triangular matrix.
 * \todo Find a more suitable place for this function.
 * \author Berkels
 */

template <typename RealType>
void solveRTransposedR ( FullMatrix<RealType> &R, const aol::Vector<RealType> &RHS, aol::Vector<RealType>& X ) {
#ifdef BOUNDS_CHECK
  if ( X.size() != R.getNumRows() )
    throw Exception ( "Length of X not compatible to R!.", __FILE__, __LINE__ );
  if ( RHS.size() != R.getNumRows() )
    throw Exception ( "Length of RightHandSide not compatible to R!.", __FILE__, __LINE__ );
  if ( R.getNumRows() != R.getNumCols() )
    throw Exception ( "Matrix R is not quadratic!.", __FILE__, __LINE__ );
#endif
  aol::Vector<RealType> Y ( R.getNumRows() );
  int i, j; // May not be unsigned due to >=0 comparison below
  RealType val;

  // first step: solve R^t*y = RHS
  for ( i = 0; i < R.getNumRows(); ++i ) {
    val = RHS[ i ];
    for ( j = 0; j < i; ++j ) {
      val -= Y[j] * R.get ( j, i );
    }
    Y[i] = val  / R.get ( i, i );
  }

  // second step: solve R*x = y
  for ( i = R.getNumRows() - 1; i >= 0; --i ) {
    val = Y[i];
    for ( j = i + 1; j < R.getNumRows() ; ++j ) {
      val -= X[j] * R.get ( i, j );
    }
    X[i] = val / R.get ( i, i );
  }
}



//! Fixpoint iteration, untested.
template <class VectorType, class InverseType>
class SemiimplicitFixpointInverse : public Op<VectorType> {
  typedef typename VectorType::Field DataType;

private:
  const Op<VectorType>& _impInvOp;
  const Op<VectorType>& _expOp;
  int _maxIt;
  DataType _eps;

public:

  SemiimplicitFixpointInverse ( const Op<VectorType>& impOp, const Op<VectorType>& expOp, DataType eps = 1E-16, int maxIt = 1000 )
      : Op<VectorType> (), _impInvOp ( InverseType ( impOp ) ), _expOp ( expOp ), _maxIt ( maxIt ), _eps ( eps ) {

    cerr << aol::color::red << "Attention: No relative stopping for SemiimplicitFixpointInverse" << aol::color::reset << endl;
  }

  void apply ( const VectorType& arg, VectorType& dest ) const {
    VectorType temp ( dest, aol::STRUCT_COPY );
    VectorType old ( dest, aol::DEEP_COPY ); // Dest is new step

    DataType step;

    for ( int i = 0; i < _maxIt; ++i ) {

      old *= -1;
      _expOp.apply ( old, temp );
      temp += arg; // b - E x_k

      _impInvOp.apply ( temp, dest ); // I x_k+1

      old += dest;
      step = old.norm (); // Break early?

      cerr << "fixpoint iteration " << aol::mixedFormat ( i ) << ", step size is " << aol::mixedFormat ( step ) << "\r";

      if ( step < _eps ) break;

      old = dest; // For next step
    }

    cerr << endl;
  }

  void applyAdd ( const VectorType& arg, VectorType& dest ) const {
    VectorType temp ( dest );
    apply ( arg, temp );
    dest += temp;
  }

};

/** Successively applies solver, reducing the residuum by a given
 *  factor a given number of times and saving intermediate results.
 *  Note that different solvers may behave differently if restarted.
 *  \author Schwen
 */
template < class SolverType, class VectorType >
class SolutionStepSaver {
protected:
  typedef typename VectorType::DataType DataType;

  SolverType&  _solver;
  const char*  _filenamemask;
  DataType     _epsilonStep;
  int          _numberOfSteps;

public:
  SolutionStepSaver ( SolverType &solver, const char* filenamemask, const DataType epsilonStep = 1.0e-2, const int numberOfSteps = 8 ) :
      _solver ( solver ),
      _filenamemask ( filenamemask ),
      _epsilonStep ( epsilonStep ),
      _numberOfSteps ( numberOfSteps ) {}

public:
  void apply ( const VectorType &arg, VectorType &dest ) {

    VectorType tmp ( arg, aol::DEEP_COPY );
    tmp *= -1.0;
    _solver.getOpReference().applyAdd ( dest, tmp );
    DataType residuum = tmp.normSqr(), accuracy = residuum;

    _solver.setStopping ( aol::STOPPING_ABSOLUTE );

    for ( int i = 0; i < _numberOfSteps ; ++i ) {
      accuracy *= _epsilonStep;
      _solver.setAccuracy ( accuracy );
      _solver.apply ( arg, dest );

      char filename[1024];
      sprintf ( filename, _filenamemask, i );
      dest.save ( filename, qc::PGM_DOUBLE_BINARY );
    }

  }

  // end class SolutionStepSaver
};


/** Utility class for treating equal zero contraints (storing constraints and projection directions, performing projection).
 *  Base class for Projecting Solvers
 *  \author Schwen
 */
template< typename VectorType >
class ProjectEqConstrSolver {
  typedef typename VectorType::DataType DataType;

protected:
  const aol::RandomAccessContainer<VectorType> &_constrVec;
  const aol::RandomAccessContainer<VectorType> &_projDirsVec;
  DataType                                      _projectThreshold;

public:
  ProjectEqConstrSolver ( const aol::RandomAccessContainer<VectorType> &constrVec,
                          const aol::RandomAccessContainer<VectorType> &projDirsVec ) :
      _constrVec ( constrVec ),
      _projDirsVec ( projDirsVec ),
      _projectThreshold ( 1.0e-7 ) {
  }

  void setProjectThreshold ( const DataType newThres ) {
    _projectThreshold = newThres;
  }

  void projectToEqConstr ( VectorType &vec ) const {
    ProjectEqConstrSolver<VectorType>::projectToEqConstr ( vec, _constrVec, _projDirsVec, _projectThreshold );
  }

  DataType checkCorrectionResiduumNeutrality ( const aol::Op< VectorType > &op ) const {
    DataType max = 0.0;
    for ( int i = 0; i < _projDirsVec.size(); ++i ) {
      if ( fabs ( ProjectEqConstrSolver<VectorType>::checkCorrectionResiduumNeutrality ( op, _projDirsVec[i] ) ) > fabs ( max ) )
        max = ProjectEqConstrSolver<VectorType>::checkCorrectionResiduumNeutrality ( op, _projDirsVec[i] );
    }
    return max;
  }

  /** Successive projection for equal-zero constraints \f[ v_i \cdot x = 0 \f]. Note that this method does not use the members of the class.
   *  \param vec  Value       vector to be projected
   *  \param constrVec        constraints \f[ v_i \f]
   *  \param projDirsVec      direction in which vec is projected: vec -= constrVec[i] * vec;
   *  \param projectThreshold only perform projection if fabs ( constrVec[i] * vec ) > projectThreshold (to avoid numerical "drift")
   */
  static void projectToEqConstr ( VectorType &vec, const aol::RandomAccessContainer< VectorType > &constrVec, const aol::RandomAccessContainer< VectorType > &projDirsVec, const DataType projectThreshold ) {
    VectorType vecCopy ( vec, aol::DEEP_COPY );

    for ( int i = 0; i < constrVec.size(); ++i ) {
      const DataType factor = constrVec[i] * vec;
      if ( fabs ( factor ) > projectThreshold ) {
        vecCopy.addMultiple ( projDirsVec[i], - factor );
      }
#ifdef VERBOSE
      cerr << endl << "Projection factor = " << factor << endl;
#endif
    }

    vec = vecCopy;

#ifdef DEBUG
    // check whether projection was successful
    for ( int i = 0; i < constrVec.size(); ++i ) {
      const DataType check = ( constrVec[i] * vec );
      if ( fabs ( check ) > projectThreshold ) {
        cerr << endl << aol::color::red << "Projection failed for " << i << ": " << check << aol::color::reset << endl;
#ifdef VERBOSE
        cerr << constrVec[i] << endl << vec << endl;
#endif
      }
    }
#endif
  }

  static DataType checkCorrectionResiduumNeutrality ( const aol::Op< VectorType > &op, const VectorType &projDirVec ) {
    VectorType dummy ( projDirVec, aol::STRUCT_COPY );
    op.apply ( projDirVec, dummy );
    return ( dummy.norm() / dummy.getTotalSize() );
  }

};


/** CG solver with projection for zero equality constraints:
 *  if fabs ( _constrVec[i] * Dest ) > some threshold, correct iterate by -( _constrVec[i] * Dest ) * _projDirsVec[i]
 *
 *  \author Schwen
 */
template < typename VectorType, typename OpType = Op<VectorType> >
class CGInverseProjectEqConstr : public IterativeInverseOp<VectorType, OpType>, public ProjectEqConstrSolver<VectorType> {
  typedef typename IterativeInverseOp<VectorType, OpType>::DataType DataType;

public:
  CGInverseProjectEqConstr ( const OpType                                 &Op,
                             const aol::RandomAccessContainer<VectorType> &constrVec,
                             const aol::RandomAccessContainer<VectorType> &projDirsVec,
                             const DataType                                Epsilon = 1e-16,
                             const int                                     MaxIter = 1000,
                             const StoppingMode                            Stop = STOPPING_UNSET,
                             ostream&                                      Out = cerr )
      : IterativeInverseOp<VectorType, OpType> ( Op, Epsilon, MaxIter, Stop, false, Out ),
      ProjectEqConstrSolver<VectorType> ( constrVec, projDirsVec ) {
  }

  virtual ~CGInverseProjectEqConstr () {}

  DataType getResSqr ( ) const {
    return this->_infoPtr->getFinalResidual();
  }

  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    DataType spa = 0, spn, q, quad;

    VectorType r ( Arg, aol::STRUCT_COPY );
    VectorType p ( Arg, aol::STRUCT_COPY );
    VectorType h ( Arg, aol::STRUCT_COPY );

    // projection before iterations
    this->projectToEqConstr ( Dest );

    this->_op.apply ( Dest, h );

    r = h;
    r -= Arg;

    p = Arg;
    p -= h;

    spn = r * r;

    this->_infoPtr->startIterations ( Arg.normSqr(), spn, "cg", "l_2 norm ^2" );

    while ( !this->_infoPtr->stoppingCriterionIsFulfilled() && ! ( this->_infoPtr->maxIterIsReached() ) && ! ( this->_infoPtr->currentResidualIsNaN() ) ) {
      this->_infoPtr->startStep();

      // case starting with second iteration
      if ( this->_infoPtr->getIterationCount() > 1 ) {
        const DataType e = spn / spa;
        p *= e;
        p -= r;
      }

      // basic iteration step
      this->_op.apply ( p, h );

      quad = p * h;
      q    = spn / quad;

      Dest.addMultiple ( p, q );

      // projection during iteration
      this->projectToEqConstr ( Dest );

      r.addMultiple ( h, q );

      spa = spn;

      // compute new residuum
      spn = r * r;

      this->_infoPtr->finishStep ( spn );
    }

    this->_op.apply ( Dest, h );
    r = h;
    r -= Arg;
    spn = r * r;

    this->_infoPtr->finishIterations ( spn );
  }

  // end class CGInverseProjectEqConstr
};


/** preconditioned CG solver with projection for zero equality constraints:
 *  if fabs ( _constrVec[i] * Dest ) > some threshold, correct iterate by -( _constrVec[i] * Dest ) * _projDirsVec[i]
 *
 *  \author Schwen
 */
template < typename VectorType, typename OpType = Op<VectorType>, typename iOpType = Op<VectorType> >
class PCGInverseProjectEqConstr : public IterativeInverseOp<VectorType, OpType>, public ProjectEqConstrSolver<VectorType> {
  typedef typename IterativeInverseOp<VectorType, OpType>::DataType DataType;

protected:
  const iOpType &_approxInverseOp;

public:
  PCGInverseProjectEqConstr ( const OpType &Op,
                              const iOpType &ApproxInverseOp,
                              const aol::RandomAccessContainer<VectorType>  &constrVec,
                              const aol::RandomAccessContainer<VectorType>  &projDirsVec,
                              const DataType Epsilon = 1e-16,
                              const int MaxIter = 50,
                              const StoppingMode Stop = STOPPING_UNSET,
                              ostream &Out = cerr )
      : IterativeInverseOp< VectorType, OpType > ( Op, Epsilon, MaxIter, Stop, false, Out ),
      ProjectEqConstrSolver<VectorType> ( constrVec, projDirsVec ),
      _approxInverseOp ( ApproxInverseOp ) {
  }

  virtual ~PCGInverseProjectEqConstr () {}

  DataType getResSqr ( ) const {
    return this->_info.getFinalResidual();
  }

  virtual void apply ( const VectorType &Arg, VectorType &Dest ) const {
    DataType alpha_numer, alpha_denom, beta_numer, beta_denom, spn;

    VectorType g ( Arg, aol::STRUCT_COPY );
    VectorType d ( Arg, aol::STRUCT_COPY );
    VectorType h ( Arg, aol::STRUCT_COPY );

    // projection before iterations
    this->projectToEqConstr ( Dest );

    this->_op.apply ( Dest, h );

#ifdef VERBOSE
    this->getOstream() << "Solver h norm " << h.norm() << endl;
#endif

    g = h;
    g -= Arg;

    h.setZero();
    _approxInverseOp.apply ( g, h );

#ifdef VERBOSE
    this->getOstream() << "h*h " << h*h << endl;
#endif

    spn = g * g;

#ifdef VERBOSE
    this->getOstream() << " spn = " << spn << endl;
#endif
    d -= h;

#ifdef VERBOSE
    this->getOstream() << "h*h " << h*h << endl;
#endif

    this->_infoPtr->startIterations ( Arg.normSqr(), spn, "p-cg", "l_2 norm ^2" );

    while ( !this->_infoPtr->stoppingCriterionIsFulfilled() && ! ( this->_infoPtr->maxIterIsReached() ) && ! ( this->_infoPtr->currentResidualIsNaN() ) ) {
      this->_infoPtr->startStep();

      beta_denom = alpha_numer = g * h;

      h.setZero();
      this->_op.apply ( d, h );
      alpha_denom = d * h;

      Dest.addMultiple ( d, alpha_numer / alpha_denom );

      // projection during iterations
      this->projectToEqConstr ( Dest );

      g.addMultiple ( h, alpha_numer / alpha_denom );

      h.setZero();
      _approxInverseOp.apply ( g, h );
      beta_numer = g * h;

      d *= ( beta_numer / beta_denom );
      d -= h;

      spn = g * g;

      this->_infoPtr->finishStep ( spn );
    }

    this->_op.apply ( Dest, h );
    g = h;
    g -= Arg;
    spn = g * g;

    this->_infoPtr->finishIterations ( spn );
  }

  // end class PCGInverseProjectEqConstr
};

} // namespace aol


#endif
