#ifndef __NEWTON_H
#define __NEWTON_H

#include <ArmijoSearch.h>
#include <imageTools.h>
#include <vectorExtensions.h>
#include <op.h>
#include <NewtonInfo.h>
#include <timestepSaver.h>

#ifdef USE_CPP11
#define AUTOPTR std::unique_ptr
#else
#define AUTOPTR auto_ptr
#endif

namespace aol {

/**
 * Newton iteration to find a root of the function \f$ F \f$.
 *
 * Supports different timestep control types, at least one of them
 * should be globally convergent (NEWTON_OPTIMAL).
 *
 * In derived classes, the abstract method "writeResult" can be implemented, and
 * the solver can be set in the constructor (or via setSolverPointer).
 * If preconditioning is to be used, then "prepareSolver" also has to be overloaded, in which
 * the solver has to be reinitialized (is necessary, since the preconditioner has to know the used matrix).
 *
 * Currently, the class does not work if VectorType != DerivativeType (see e.g. the beginning of apply()),
 * thus this is checked with an AssertAtCompileTimeThatTypesAreEqual in the constructor.
 *
 * Originally, VectorType != DerivativeType used to work. If you need it and repair it,
 * you can of course remove the assertion.
 *
 * \author Berkels
 */
template < typename RealType,
           typename VectorType,
           typename DerivativeType = VectorType,
           typename SecondDerivativeType = qc::FastUniformGridMatrix<RealType,qc::QC_2D> >
class NewtonIterationBase
  : public ArmijoLineSearchUsingOp<RealType, VectorType, VectorType, DerivativeType >
{
protected:
  enum NEWTON_FLAGS {
    //! if you want, the iterative solver's tolerance can be
    //! automatically adapted to the current residual in
    //! each Newton step. See NewtonInfo::computeSolverTOL() for the
    //! explicit algorithm.
    COMPUTE_SOLVER_TOL = 1
  };
  int _configurationFlags;

  //! Function that we want to find a root of.
  const aol::Op<VectorType, DerivativeType> &_F;

  //! Derivative of _F.
  const aol::Op<VectorType, SecondDerivativeType> &_DF;

  //! Stores the derivative of F and is deleted in the destructor.
  mutable SecondDerivativeType *_pMatDF;

  //! _pSolver has to be intialized and deleted by the derived class!
  mutable aol::InverseOp<DerivativeType> *_pSolver;

  mutable DeleteFlagPointer<NewtonInfo<RealType> > _pInfo;

  //! variable that is only needed by LNEQ_RES timestep control
  mutable AUTOPTR<DerivativeType> _NLEQ_ERR_delta_x_bar;

public:
  NewtonIterationBase ( SecondDerivativeType *PMatDF,
                        const aol::Op<VectorType, DerivativeType> &F,
                        const aol::Op<VectorType, SecondDerivativeType> &DF,
                        const int MaxIterations = 50,
                        const RealType StopEpsilon = 1.e-6,
                        const bool WriteTimeSteps = false,
                        const char *BaseSaveName = NULL )
  : ArmijoLineSearchUsingOp<RealType, VectorType, VectorType, DerivativeType >(0.25, MaxIterations, aol::ZOTrait<RealType>::one),
    _configurationFlags ( 0 ),
    _F(F),
    _DF(DF),
    _pMatDF(PMatDF),
    _pSolver ( NULL ),
    _pInfo ( new aol::NewtonInfo<RealType> (StopEpsilon, MaxIterations, aol::Sqr(StopEpsilon), 1000), true ) {
    AssertAtCompileTimeThatTypesAreEqual<VectorType,DerivativeType> (); // Cf. comment in front of class
    this->setWriteTimeSteps( WriteTimeSteps );
    if( BaseSaveName != NULL )
      this->setSaveName( BaseSaveName );
  }

  NewtonIterationBase ( SecondDerivativeType *PMatDF,
                        const aol::Op<VectorType, DerivativeType> &F,
                        const aol::Op<VectorType, SecondDerivativeType> &DF,
                        aol::NewtonInfo<RealType> & info,
                        const bool WriteTimeSteps = false,
                        const char *BaseSaveName = NULL )
  : ArmijoLineSearchUsingOp<RealType, VectorType, VectorType, DerivativeType >(0.25, info.getMaxIterations(), aol::ZOTrait<RealType>::one),
    _configurationFlags ( 0 ),
    _F(F),
    _DF(DF),
    _pMatDF(PMatDF),
    _pSolver(NULL),
    _pInfo ( &info, false ) {
    this->setWriteTimeSteps( WriteTimeSteps );
    if( BaseSaveName != NULL )
      this->setSaveName( BaseSaveName );
  }

  virtual ~NewtonIterationBase(){
    delete _pMatDF;
    // delete NULL has no effect
    delete _pSolver;
  };
  /* NewtonIterationBase takes ownership of the memory pointed at by PSolver,
   * i.e. it will delete the pointer in the destructor.
   */
  virtual void setSolverPointer( aol::InverseOp<DerivativeType> *PSolver ) const {
    if ( _pSolver )
      delete _pSolver;
    _pSolver = PSolver;
  }
  NewtonInfo<RealType> & getNewtonInfo () {
    return *_pInfo;
  }

  const NewtonInfo<RealType> & getNewtonInfo () const {
    return *_pInfo;
  }

  void setNewtonInfoPointer ( aol::NewtonInfo<RealType> * pInfo ) {
    _pInfo.reset ( pInfo, false );
  }

  void setNewtonInfoReference ( aol::NewtonInfo<RealType> & pInfo ) {
    _pInfo.reset ( &pInfo, false );
  }

  //! The class has an inherited member function aol::TimestepSaver::setQuietMode().
  //! Of course, this does not what an NewtonIterationBase::setQuietMode()
  //! is assumed to do, so we override it. Note that the TimestepSaver::setQuietMode()
  //! function is now only accessible via explicit cast.
  virtual void setQuietMode ( const bool quiet ) {
    _pInfo->setQuietMode ( quiet );
  }

  void setConfigurationFlags( int flags ) {
    _configurationFlags = flags;
  }

  //! \short shadow variable from ArmijoLineSearchUsingOp and use the one
  //! from NewtonInfo.
  virtual RealType getTauMin () const {  return getNewtonInfo().getTauMin();  }
  virtual RealType getTauMax () const {  return getNewtonInfo().getTauMax();  }

  virtual void setTauMin ( const RealType tauMin ) {  getNewtonInfo().setTauMin ( tauMin );  }
  virtual void setTauMax ( const RealType tauMax ) {  getNewtonInfo().setTauMax ( tauMax );  }

protected:
  virtual void writeResult( const VectorType &/*Dest*/, const char* /*FileName*/ ) const {
  }
  virtual void writeUpdate ( const DerivativeType &/*descentDir*/, const char * /*FileName*/ ) const {
  }
  virtual RealType fNorm( const DerivativeType &Position ) const{
    return (Position.norm());
  }
  virtual void postProcess( VectorType &/*Position*/, const int /*Iteration*/ ) const {
  }
  virtual void prepareSolver() const {
  }
  virtual void applySolver( const DerivativeType &Rhs, DerivativeType &DescentDir ) const {
    _pSolver->apply( Rhs, DescentDir );
  }
  virtual void smoothDirection( const VectorType &/*CurrentPosition*/, DerivativeType &/*Direction*/ ) const {
  }
  //! \brief compute the current error norm for the stopping criterion
  //! \param[in] x_k current position 
  //! \param[in] delta_x_k difference to current position
  //! \param[in] F_x_k function value at x_k
  //! \param[in] true, if this is the intial norm calculation (before iterations), false otherwise
  //! \returns error norm
  virtual RealType computeErrorNorm ( const VectorType& /*x_k*/, const VectorType& /*delta_x_k*/, const DerivativeType& F_x_k, bool /*initial*/ ) const {
    return fNorm ( F_x_k );
  }
  virtual RealType ArmijoLineSearchHelpFunction_evaluate( const DerivativeType &DescentDir, const VectorType &CurrentPosition, const RealType timestepWidth ) const{
    VectorType newPosition(CurrentPosition);
    DerivativeType temp(DescentDir, aol::STRUCT_COPY);
    newPosition.addMultiple ( DescentDir, timestepWidth );
    _F.apply( newPosition, temp);
    return (temp.normSqr()/2.);
  }
  // Attention: ArmijoLineSearchHelpFunction_evaluateDerivative only works, if _matDF is filled with DF(position)!
  virtual RealType ArmijoLineSearchHelpFunction_evaluateDerivative( const DerivativeType &DescentDir, const VectorType &Position) const{
    DerivativeType tmp (DescentDir, aol::STRUCT_COPY);
    DerivativeType tmp2 (DescentDir, aol::STRUCT_COPY);
    //SecondDerivativeType* pMatDF = newSecondDerivativeTempObject();
    //_DF.apply(Position, *pMatDF);
    _pMatDF->apply(DescentDir, tmp);
    _F.apply(Position, tmp2);
    RealType temp = (tmp)*(tmp2);
    //delete pMatDF;
    return (temp);
  }

  //! extreme simple timestep width control, just ensures, that FNormNew < FNorm
  RealType getSimpleNewtonTimestepWidth( const DerivativeType &DescentDir,
                                         const VectorType &CurrentPosition ) const {
    DerivativeType tmp(DescentDir, aol::STRUCT_COPY);
    DerivativeType tmp2(DescentDir, aol::STRUCT_COPY);
    DerivativeType descentDir(DescentDir);

    tmp = CurrentPosition;

    tmp += DescentDir;

    _F.apply(CurrentPosition, tmp2);
    RealType FNorm = tmp2.norm();
    _F.apply(tmp, tmp2);
    RealType FNormNew = tmp2.norm();
    RealType tau = 1.;
    while( (FNormNew > FNorm) && (tau >= this->_tauMin ) ){
      tau = tau * 0.5;
      descentDir *= 0.5;
      tmp -= descentDir;  // Here v_i = v + 2^{-i}d => v^{i+1} = v^i - 2^{i+1}d
      _F.apply(tmp, tmp2);
      FNormNew = tmp2.norm();
    }

    // FNormNew < FNorm not fulfilled for tau >= _tauMin, so we don't want the step to be done.
    // The stopping criterion also handles this case.
    if ( tau < this->_tauMin )
      tau = 0;

#ifdef VERBOSE
    cerr << "tau = " << tau << endl;
#endif

    return tau;
  }

  /**
   * Calculates the step size based on "Schaback, Werner - Numerische Mathematik, 4. Auflage, Seite 129".
   *
   * Is global convergent, as long as F fulfils certain regularity conditions.
   */
  RealType getNewtonOptimalTimestepWidth( const DerivativeType &DescentDir,
                                          const VectorType &CurrentPosition ) const{
    const RealType DescentDirNormSqr = DescentDir.normSqr();
    if( DescentDirNormSqr > 0. ){
      // initial guess for L: L = ||F(y)-F(x)-F'(x)(y-x)|| / ||y-x||^2, where x = CurrentPosition, y = x + DescentDir
      VectorType newPosition(CurrentPosition);
      newPosition += DescentDir;

      DerivativeType* pTmp = new DerivativeType(DescentDir, aol::STRUCT_COPY);
      _F.apply(CurrentPosition, *pTmp);
      RealType fNorm = pTmp->norm();
      _pMatDF->applyAdd(DescentDir, *pTmp);
      *pTmp *= -1.;
      _F.applyAdd(newPosition, *pTmp);
      RealType L = (*pTmp).norm()/DescentDirNormSqr;
      delete pTmp;

      // If F is locally linear and matches the gradient with numerical perfection take the full step
      if ( L < 1e-15 ) return aol::ZOTrait<RealType>::one;

#ifdef VERBOSE
      cerr << endl << aol::color::brown << "Calculating tau Newton optimal (starting with L = " << L << ")\n";
#endif
      RealType tau = aol::Min(fNorm/(2*L*DescentDirNormSqr), aol::ZOTrait<RealType>::one);
      RealType temp1 = fNorm -
          sqrt(2*NewtonIterationBase<RealType, VectorType, DerivativeType, SecondDerivativeType>::ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, tau ));
      RealType temp2 = tau*(fNorm - L*tau*DescentDirNormSqr);
      RealType tauCandidate;
      if( temp1 < temp2 ){
        do{
          L *= 2.;
          tau = aol::Min(fNorm/(2*L*DescentDirNormSqr), aol::ZOTrait<RealType>::one);
          temp1 = fNorm -
              sqrt(2*NewtonIterationBase<RealType, VectorType, DerivativeType, SecondDerivativeType>::ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, tau ));
          temp2 = tau*(fNorm - L*tau*DescentDirNormSqr);
          // Prevent the while loop from getting stuck if DescentDir is not a descent direction.
          if ( tau < this->_tauMin )
            temp1 = temp2, tau = 0.;
        } while( temp1 < temp2 );
      }
      else{
        do{
          // Compute tauCandidate, temp1 and temp2 with L/2 instead of L to find out if a bigger step size is still feasible.
          tauCandidate = aol::Min(fNorm/(L*DescentDirNormSqr), aol::ZOTrait<RealType>::one);
          if(tauCandidate == aol::ZOTrait<RealType>::one)
            break;
          temp1 = fNorm -
              sqrt(2*NewtonIterationBase<RealType, VectorType, DerivativeType, SecondDerivativeType>::ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, tauCandidate ));
          temp2 = tauCandidate*(fNorm - 0.5*L*tauCandidate*DescentDirNormSqr);

          if ( temp1 >= temp2 ){
            L /= 2.;
            tau = tauCandidate;
          }
        } while( temp1 >= temp2 );
      }

#ifdef VERBOSE
    cerr << "accepted tau = " << tau << ", L = " << L << endl << endl << aol::color::reset;
#endif
      return tau;
    }
    else
      return 0.;
  }

  /**
   *  Calculates stepsize following algorithm NLEQ-RES from Deuflhard's "Newton Methods
   *  for Nonlinear Problems" (Heidelberg 2004), p. 131
   */
  RealType getNewtonTimestep_NLEQ_RES ( const DerivativeType & /* descent direction */ delta_x,
                                        const VectorType & /* current position */ x) const {
    RealType lambda = 1.;
    RealType lambda_prime = 1.;
    RealType mu = 0.;
    static RealType mu_prime = 0.;
    static RealType norm_F_x = 0.;

    DerivativeType F_x (delta_x, aol::STRUCT_COPY);
    DerivativeType F_x_new (delta_x, aol::STRUCT_COPY);
    _F.apply(x, F_x);

    VectorType x_new ( x, aol::STRUCT_COPY );

    RealType new_norm_F_x = F_x.norm();

    if (mu_prime == 0.) // first step:
      norm_F_x = new_norm_F_x;

    mu = norm_F_x / new_norm_F_x;
    norm_F_x = new_norm_F_x;
    lambda = Min ( aol::ZOTrait<RealType>::one, mu );

    while ( lambda >= this->_tauMin ) {
      x_new = x;
      x_new.addMultiple ( delta_x, lambda );
      _F.apply ( x_new, F_x_new );

      RealType Theta = F_x_new.norm() / norm_F_x;
      // tmp = F(x_new) - (1-lambda) F (x)
      DerivativeType tmp (F_x_new);
      tmp.addMultiple ( F_x, -(1. - lambda) );

      mu_prime = 0.5 * norm_F_x * lambda * lambda / tmp.norm();

      if ( Theta >= 1. ) {
        lambda_prime = Min ( mu_prime, static_cast<RealType>( 0.5 ) * lambda );
        lambda = lambda_prime;
        continue;
      }
      else
        lambda_prime = Min ( aol::ZOTrait<RealType>::one, mu_prime );
      if ( lambda_prime >= 4. * lambda ) {
        lambda_prime = lambda;
        continue;
      }
      else
        return lambda;
    }
    // convergence failure:
    return 0.;
  }

  /**
   *  Calculates stepsize following algorithm NLEQ-ERR from Deuflhard's "Newton Methods
   *  for Nonlinear Problems" (Heidelberg 2004), p. 148 sq.
   */
  RealType getNewtonTimestep_NLEQ_ERR ( const DerivativeType & /* descent direction */ delta_x,
                                        const VectorType & /* current position */ x) const {

      static RealType lambda = 0.;
      static RealType norm_delta_x = 0.;
      static RealType norm_delta_x_bar = 0.;

      if ( !_NLEQ_ERR_delta_x_bar.get() ) {
        _NLEQ_ERR_delta_x_bar.reset ( new DerivativeType ( delta_x, STRUCT_COPY ) );
        lambda = 0.;
      }
      DerivativeType & delta_x_bar = *_NLEQ_ERR_delta_x_bar;
      RealType mu;
      RealType new_norm_delta_x = delta_x.norm();
      if ( lambda == 0. )
        mu = 1.;
      else {
        delta_x_bar -= delta_x;
        mu = norm_delta_x * norm_delta_x_bar / ( delta_x_bar.norm() * new_norm_delta_x ) * lambda;
      }
      lambda = Min ( ZOTrait<RealType>::one, mu );
      norm_delta_x = new_norm_delta_x;

      while ( lambda > this->_tauMin ) {
        VectorType x_new ( x );
        x_new.addMultiple ( delta_x, lambda );

        DerivativeType F_new (delta_x, aol::STRUCT_COPY);
        _F.apply( x_new, F_new );
        F_new *= -1.;
        delta_x_bar.setZero();
        this->applySolver( F_new, delta_x_bar );
        norm_delta_x_bar = delta_x_bar.norm();

        RealType Theta = norm_delta_x_bar / norm_delta_x;
        // tmp = delta_x_bar - (1 - lambda) delta_x
        DerivativeType tmp ( delta_x_bar );
        tmp.addMultiple ( delta_x, -(1. - lambda) );
        RealType mu_prime = 0.5 * norm_delta_x * lambda * lambda / tmp.norm();

        RealType lambda_prime;
        if ( Theta >= 1. ) {
          lambda_prime = Min ( mu_prime, static_cast<RealType>( 0.5 ) * lambda );
          lambda = lambda_prime;
          continue;
        }
        else
          lambda_prime = Min ( ZOTrait<RealType>::one, mu_prime );

        if ( lambda_prime > 4. * lambda ) {
          lambda = lambda_prime;
          continue;
        }
        return lambda;
      }
      // convergence failure:
      return 0.;
  }

public:
  virtual void apply( const VectorType &Arg, VectorType &Dest ) const {

    const VectorType & x_0 = Arg;
    VectorType & x_k = Dest;

    DerivativeType F_x_k     (x_k, aol::STRUCT_COPY);
    DerivativeType tmp       (x_k, aol::STRUCT_COPY);
    DerivativeType delta_x_k (x_k, aol::STRUCT_COPY);

    RealType FNorm = 1.0;
    RealType FNormNew = 1.0;
    RealType tau = 1.0;
    x_k = x_0;

    _pInfo->startIterations();
    _F.apply(x_k, F_x_k);
    FNorm = computeErrorNorm ( x_k, delta_x_k, F_x_k, true ); 

    // Since we post process the current solution after each iteration,
    // it shouldn't hurt to do this on the initial data. This way each
    // iteration starts with post processed input.
    postProcess( x_k, 0 );

    while ( FNorm > _pInfo->getAccuracy() && !(_pInfo->maxIterIsReached()) && tau > 0. ) {
      _pInfo->startStep();

      // Newton iteration given by x^{k+1} = x^k - DF(x^k)^{-1}(F(x^k))
      // we store DF in matDF;
      if ( (this->_pInfo->getIterationCount() - 1) % _pInfo->getNumDerivativeHoldingSteps() == 0)
        _DF.apply(x_k, *_pMatDF);
      F_x_k *= -1.;

      // x^k is not a good initial guess for the descent direction, better use zero.
      delta_x_k.setZero();

      prepareSolver();
      if ( _pSolver == NULL )
        throw Exception( "_pSolver is not initialized!", __FILE__, __LINE__);

      // store solver TOL
      RealType solverTOL = getNewtonInfo().getSolverInfo().getAccuracy();
      // if adaption of solver TOL is requested, do so
      if ( (_configurationFlags & COMPUTE_SOLVER_TOL) && _pInfo->getIterationCount() > 0 )
        _pInfo->computeSolverTOL ();

      applySolver( F_x_k, delta_x_k );

      const aol::IterativeInverseOp<DerivativeType>* pIterativeSolver = dynamic_cast<const aol::IterativeInverseOp<DerivativeType>*>( _pSolver );
      if( pIterativeSolver != NULL ){
        // ensure that correct iteration count and residual
        // is saved in NewtonInfo's SolverInfo
        _pInfo->getSolverInfo().getIterationInfoFrom( pIterativeSolver->getSolverInfoReference() );

        // restore previously stored solver TOL
        _pInfo->getSolverInfo().setAccuracy ( solverTOL );
      }
      
      smoothDirection( x_k, delta_x_k );
      
      tau = getStepsize ( delta_x_k, x_k, tau );
      
      RealType normUpdate = delta_x_k.norm();
      x_k.addMultiple ( delta_x_k, tau );

      // Post processing possibly increases the residual, therefore it should be done
      // before FNormNew is calculated.
      postProcess( x_k, this->_pInfo->getIterationCount() );

      _F.apply(x_k, F_x_k);
      FNormNew = computeErrorNorm ( x_k, delta_x_k, F_x_k, false );

#ifdef VERBOSE
    cerr << aol::color::light_grey;
    cerr << "FNorm before iteration:   " << FNorm << endl;
    cerr << "FNorm after iteration :   " << FNormNew << endl;
    cerr << aol::color::reset;
#endif

      FNorm = FNormNew;

      _pInfo->finishStep(FNormNew, normUpdate, tau);
      // output
      if( this->checkSaveConditions( _pInfo->getIterationCount()) ){
        string filename = strprintf( "%s%s_newtonstep_%03d_%.5f",
                                    this->getSaveDirectory(),
                                    this->getSaveName(),
                                    this->_pInfo->getIterationCount(),
                                    static_cast<double>(FNormNew) );
        writeUpdate( delta_x_k, filename.c_str() );
        writeResult( x_k, filename.c_str() );
      }

    } // end while
    _pInfo->finishIterations();
  }

  virtual RealType getStepsize ( const DerivativeType & descentDir, const VectorType & Dest, const RealType tau_before ) const {
    switch( _pInfo->getTimestepController() ) {
      case NewtonInfo<RealType>::ARMIJO:
        return this->getTimestepWidthWithArmijoLineSearch ( descentDir, Dest, tau_before );

      case NewtonInfo<RealType>::NEWTON_OPTIMAL:
        return getNewtonOptimalTimestepWidth( descentDir, Dest );

      case NewtonInfo<RealType>::NLEQ_RES:
        return getNewtonTimestep_NLEQ_RES ( descentDir, Dest );

      case NewtonInfo<RealType>::NLEQ_ERR:
        return getNewtonTimestep_NLEQ_ERR ( descentDir, Dest );

      case NewtonInfo<RealType>::NEWTON_SIMPLE:
        return getSimpleNewtonTimestepWidth ( descentDir, Dest );

      case NewtonInfo<RealType>::NO_TIMESTEP_CONTROL:
        return aol::NumberTrait<RealType>::one;

      case NewtonInfo<RealType>::WOLFE:
        return this->getTimestepWidthWithPowellWolfeLineSearch( descentDir, Dest, tau_before );

      default:
        throw Exception( "Unhandled timestep control mode", __FILE__, __LINE__);
    }
  }
  virtual void applyAdd( const VectorType &/*Arg*/, VectorType &/*Dest*/ ) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
  void setSolverQuietMode( const bool QuietMode ) {
    if ( _pSolver )
      _pSolver->setQuietMode( QuietMode );
  }
  void setTimestepController( typename NewtonInfo<RealType>::TIMESTEP_CONTROLLER TimestepController ) {
    _pInfo->setTimestepController ( TimestepController );
  }
};

/**
 * This class implements Newton's method for minimizing a scalar function.
 * It has to be provided the objective function, its derivative, and the second derivative.
 * In contrast to NewtonIterationBase, the linesearch/backtracking algorithm checks whether
 * the objective function has sufficiently decreased instead of the norm of the derivative.
 *
 * \author Wirth
 */
template < typename RealType, typename VectorType, typename DerivativeType, typename SecondDerivativeType >
class NewtonMinimizationBase
  : public aol::NewtonIterationBase< RealType, VectorType, DerivativeType, SecondDerivativeType > {

protected:
  // the scalar objective function or "energy"
  const aol::Op<VectorType, aol::Scalar<RealType> > &_E;

  // Cached value
  mutable RealType _energyAtPosition;
public:
  NewtonMinimizationBase( const aol::Op<VectorType, aol::Scalar<RealType> > &F,
                          const aol::Op<VectorType, DerivativeType> &DF,
                          const aol::Op< VectorType, SecondDerivativeType > &D2F,
                          SecondDerivativeType *PMatD2F,
                          const int MaxIterations = 50,
                          const RealType StopEpsilon = 1.e-6,
                          const bool WriteTimeSteps = false,
                          const char *BaseSaveName = NULL ) :
    aol::NewtonIterationBase< RealType, VectorType, DerivativeType, SecondDerivativeType >
                                                                      ( PMatD2F, DF, D2F,
                                                                        MaxIterations,
                                                                        StopEpsilon,
                                                                        WriteTimeSteps,
                                                                        BaseSaveName),
    _E( F ) {
  }

  NewtonMinimizationBase( const aol::Op<VectorType, aol::Scalar<RealType> > &F,
                          const aol::Op<VectorType, DerivativeType> &DF,
                          const aol::Op< VectorType, SecondDerivativeType > &D2F,
                          SecondDerivativeType *PMatD2F,
                          aol::NewtonInfo<RealType> & Info,
                          const bool WriteTimeSteps = false,
                          const char *BaseSaveName = NULL ) :
    aol::NewtonIterationBase< RealType, VectorType, DerivativeType, SecondDerivativeType >
                                                                      ( PMatD2F, DF, D2F,
                                                                        Info,
                                                                        WriteTimeSteps,
                                                                        BaseSaveName),
    _E( F ) {
  }

  virtual ~NewtonMinimizationBase () {
  }

protected:
  /**
   * Returns the scalar objective function evaluated at CurrentPosition+DescentDir*timestepWidth.
   */
  virtual RealType ArmijoLineSearchHelpFunction_evaluate( const DerivativeType &DescentDir,
                                                          const VectorType &CurrentPosition,
                                                          const RealType timestepWidth ) const {
    VectorType newPosition( CurrentPosition );
    newPosition.addMultiple ( DescentDir, timestepWidth );
    aol::Scalar<RealType> Energy;
    _E.apply( newPosition, Energy);
    if ( timestepWidth == 0 )
      _energyAtPosition = Energy[0];
    return ( Energy );
  }

  /**
   * Returns the dot product of the energy derivative and the descent direction.
   */
  virtual RealType ArmijoLineSearchHelpFunction_evaluateDerivative( const DerivativeType &DescentDir,
                                                                    const VectorType &Position) const{
    DerivativeType tmp( DescentDir );
    this->_F.apply( Position, tmp );
    return ( tmp * DescentDir );
  }

  /**
   * Returns the dot product of the energy derivative evaluated at CurrentPosition+DescentDir*timestepWidth and the descent direction.
   */
  virtual RealType ArmijoLineSearchHelpFunction_evaluateDerivativeWithTau( const DerivativeType &DescentDir,
                                                                           const VectorType &Position,
                                                                           const RealType timestepWidth ) const{
    VectorType newPosition( Position );
    DerivativeType temp( Position, aol::STRUCT_COPY );
    newPosition.addMultiple( DescentDir, timestepWidth);
    this->_F.apply( newPosition, temp );
    return ( temp * DescentDir );
  }
};


template <typename ConfiguratorType, typename VectorType = aol::Vector<typename ConfiguratorType::RealType>, typename SecondDerivativeType = qc::FastUniformGridMatrix<typename ConfiguratorType::RealType, ConfiguratorType::Dim> >
class NewtonIteration{
};

template <typename ConfiguratorType, typename SecondDerivativeType>
class NewtonIteration<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType>, SecondDerivativeType >
: public NewtonIterationBase< typename ConfiguratorType::RealType,
                              aol::Vector<typename ConfiguratorType::RealType>,
                              aol::Vector<typename ConfiguratorType::RealType>,
                              SecondDerivativeType > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  aol::SSORPreconditioner<aol::Vector<RealType>, SecondDerivativeType> _ssor;

  void initSolver () {
    this->_pInfo->getSolverInfo().setMaxIterations ( 10000 );
    aol::PCGInverse<aol::Vector<RealType> >* pSolver = new aol::PCGInverse<aol::Vector<RealType> >( *(this->_pMatDF), _ssor, this->_pInfo->getSolverInfo() );
    pSolver->setStopping( aol::STOPPING_ABSOLUTE );
    this->_pSolver = pSolver;
  }
public:
  NewtonIteration( const typename ConfiguratorType::InitType &Initializer,
                   const aol::Op<aol::Vector<RealType> > &F,
                   const aol::Op<aol::Vector<RealType>, SecondDerivativeType> &DF,
                   const int MaxIterations = 50,
                   const RealType StopEpsilon = 1.e-6,
                   const bool WriteTimeSteps = false,
                   const char *BaseSaveName = NULL)
  : NewtonIterationBase< RealType,
                         aol::Vector<RealType>,
                         aol::Vector<RealType>,
                         SecondDerivativeType
                       >
                       ( new SecondDerivativeType(Initializer), F, DF, MaxIterations, StopEpsilon, WriteTimeSteps, BaseSaveName),
    _grid(Initializer),
    _ssor(*(this->_pMatDF))
  {
    initSolver ();
  }
  //! The first argument is only a dummy for compatibility reasons with NewtonIterationNB.
  NewtonIteration( const ConfiguratorType &/*Config*/,
                   const typename ConfiguratorType::InitType &Initializer,
                   const aol::Op<aol::Vector<RealType> > &F,
                   const aol::Op<aol::Vector<RealType>, SecondDerivativeType> &DF,
                   const int MaxIterations = 50,
                   const RealType StopEpsilon = 1.e-6,
                   const bool WriteTimeSteps = false,
                   const char *BaseSaveName = NULL)
  : NewtonIterationBase< RealType,
                         aol::Vector<RealType>,
                         aol::Vector<RealType>,
                         SecondDerivativeType
                       >
                       ( new SecondDerivativeType(Initializer), F, DF, MaxIterations, StopEpsilon, WriteTimeSteps, BaseSaveName),
    _grid(Initializer),
    _ssor(*(this->_pMatDF))
  {
    initSolver ();
  }
  virtual ~NewtonIteration() {
  }
  virtual void writeResult( const aol::Vector<RealType> &Dest, const char* FileName ) const {
    qc::writeImage<RealType> ( _grid, Dest, FileName );
  }
};

/**
 * The template SolverType so far only admits solvers that have a constructor
 * with the same argument list as PCGInverse.
 *
 * Recommended SolverType:
 * PCGInverse for symmetric matrices
 * PBiCGStabInverse for non-symmetric matrices
 *
 * Uses DiagonalPreconditioner as preconditioner.
 *
 * \author Berkels
 */
template <typename MatrixType, typename SolverType, int NumComponents>
class NewtonIterationSparseBlockMatrix
: public NewtonIterationBase< typename MatrixType::DataType,
                              aol::MultiVector<typename MatrixType::DataType >,
                              aol::MultiVector<typename MatrixType::DataType >,
                              aol::SparseBlockMatrix<MatrixType> > {
  typedef typename MatrixType::DataType RealType;
private:
  const aol::IdentityOp<aol::MultiVector<RealType> > _identity;
  mutable aol::DiagonalPreconditioner< aol::MultiVector<RealType> >* _pPrecond;
public:
  NewtonIterationSparseBlockMatrix( const aol::Op<aol::MultiVector<RealType> > &F,
                                    const aol::Op<aol::MultiVector<RealType>, aol::SparseBlockMatrix<MatrixType> > &DF,
                                    const int MaxIterations = 50,
                                    const RealType StopEpsilon = 1.e-6,
                                    const bool WriteTimeSteps = false,
                                    const char *BaseSaveName = NULL)
  : NewtonIterationBase< RealType,
                         aol::MultiVector<RealType>,
                         aol::MultiVector<RealType>,
                         aol::SparseBlockMatrix<MatrixType>
                       >
                       (new aol::SparseBlockMatrix<MatrixType>(NumComponents, NumComponents), F, DF, MaxIterations, StopEpsilon, WriteTimeSteps, BaseSaveName),
    _identity ( ),
    _pPrecond(NULL)
  {
    this->_pInfo->getSolverInfo().setMaxIterations ( 10000 );
    this->_pSolver = new SolverType (*this->_pMatDF, _identity, this->_pInfo->getSolverInfo() );
  }
  virtual ~NewtonIterationSparseBlockMatrix() {
    if ( _pPrecond )
      delete _pPrecond;
  }

  virtual void prepareSolver() const {
    delete this->_pSolver;
    if ( _pPrecond ){
      delete _pPrecond;
      _pPrecond = NULL;
    }
    _pPrecond = new aol::DiagonalPreconditioner< aol::MultiVector<RealType> >(*this->_pMatDF);
    SolverType *pSolver = new SolverType (*this->_pMatDF, *_pPrecond, this->_pInfo->getSolverInfo() );
    pSolver->setStopping( STOPPING_ABSOLUTE );
    this->_pSolver = pSolver;
  }
};

/**
 * \brief This class implements a quasi-Newton method (the DFP-method) for minimizing a scalar function.
 *
 * It has to be provided the scalar function and its derivative.
 *
 * \author Wirth
 * \ingroup Optimization
 */
template < typename RealType, typename VectorType, typename DerivativeType >
class QuasiNewtonIteration
  : public aol::NewtonMinimizationBase< RealType, VectorType, DerivativeType, aol::IdentityOp<VectorType> > {
protected:
  //! Stores the vectors Z such that the approximation of the inverse Hesse matrix is sum of cZZ^T and the identity.
  mutable aol::RandomAccessContainer< DerivativeType > _Z;
  //! Stores the scalars c such that the approximation of the inverse Hesse matrix is sum of cZZ^T and the identity.
  mutable std::vector< RealType > _c;
  //! number of steps before resetting the approximate inverse Hesse matrix to the identity
  const int _resetNum;

public:
  QuasiNewtonIteration( const aol::Op<VectorType, aol::Scalar<RealType> > &E,
                        const aol::Op<VectorType, DerivativeType> &F,
                        const int MaxIterations = 50,
                        const RealType StopEpsilon = 1.e-6,
                        const int ResetNum = 10,
                        const bool WriteTimeSteps = false,
                        const char *BaseSaveName = NULL ) :
    aol::NewtonMinimizationBase< RealType, VectorType, DerivativeType, aol::IdentityOp<VectorType> >( E, F, aol::NullOp<VectorType, aol::IdentityOp<VectorType> >(), NULL, MaxIterations, StopEpsilon, WriteTimeSteps, BaseSaveName),
    _resetNum( ResetNum ) {
    this->_pInfo->getSolverInfo().setQuietMode( true );
  }

  QuasiNewtonIteration( const aol::Op<VectorType, aol::Scalar<RealType> > &E,
                        const aol::Op<VectorType, DerivativeType> &F,
                        aol::NewtonInfo<RealType> & Info,
                        const int ResetNum = 10,
                        const bool WriteTimeSteps = false,
                        const char *BaseSaveName = NULL ) :
    aol::NewtonMinimizationBase< RealType, VectorType, DerivativeType, aol::IdentityOp<VectorType> >( E, F, aol::NullOp<VectorType, aol::IdentityOp<VectorType> >(), NULL, Info, WriteTimeSteps, BaseSaveName),
    _resetNum( ResetNum ) {
    this->_pInfo->getSolverInfo().setQuietMode( true );
  }

  virtual ~QuasiNewtonIteration () {}

protected:
  // CurrentPosition contains the currently optimal point,
  // DescentDir the linesearch direction,
  // Derivative the energy derivative at CurrentPosition.
  // All parameters at the beginning contain the values at the initial position,
  // and their new values at the end.
  virtual void getTauAndUpdateDescentDir( DerivativeType &DescentDir,
                                          VectorType &CurrentPosition,
                                          DerivativeType &Derivative,
                                          aol::Scalar<RealType> &Energy,
                                          aol::Vector<RealType> &Tau ) const {
    switch( this->_pInfo->getTimestepController() ) {
      case aol::NewtonInfo<RealType>::ARMIJO: getTauAndUpdateDescentDirArmijo( DescentDir, CurrentPosition, Derivative, Energy, Tau );
        break;
      case aol::NewtonInfo<RealType>::NEWTON_OPTIMAL:
      case aol::NewtonInfo<RealType>::NEWTON_SIMPLE:
      case aol::NewtonInfo<RealType>::NO_TIMESTEP_CONTROL:
        throw Exception( "aol::NewtonInfo<RealType>::NEWTON_OPTIMAL, NEWTON_SIMPLE, or NO_TIMESTEP_CONTROL is not implemented!", __FILE__, __LINE__);
        break;
      case aol::NewtonInfo<RealType>::WOLFE:
      default: getTauAndUpdateDescentDirWolfe( DescentDir, CurrentPosition, Derivative, Energy, Tau ); // only this ensures the approximation to the Hessian to be positive definite and hence the construction of descent directions
    }
  }

  // CurrentPosition contains the currently optimal point,
  // DescentDir the linesearch direction,
  // Derivative the energy derivative at CurrentPosition.
  // All parameters at the beginning contain the values at the initial position,
  // and their new values at the end.
  virtual void getTauAndUpdateDescentDirArmijo( DerivativeType &DescentDir,
                                                VectorType &CurrentPosition,
                                                DerivativeType &Derivative,
                                                aol::Scalar<RealType> &Energy,
                                                aol::Vector<RealType> &Tau ) const {
    // compute stepsize from Armijo condition

    Tau[0] = aol::Min( aol::Max( Tau[0], this->getTauMin() ), this->getTauMax() );

    RealType energyNew = this->ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, Tau[0] ), energyBackup;
    RealType dE = Derivative * DescentDir;
    RealType G = (energyNew - Energy.v)/(dE*Tau[0]);

    if( G >= this->getSigma() && Energy.v >= energyNew ){
      //time step too small
      do{
        Tau[0] *= 2.;
        energyBackup = energyNew;
        energyNew = this->ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, Tau[0] );
        G = (energyNew - Energy.v)/(dE*Tau[0]);
      } while (G >= this->getSigma() && Energy.v >= energyNew && Tau[0] <= this->getTauMax());
      Tau[0] *= 0.5;
    }// if
    else{
      // time step too large
      do{
        if (Tau[0] > this->getTauMin())
          Tau[0] *= 0.5;
        energyBackup = energyNew;
        energyNew = this->ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, Tau[0] );
        G = (energyNew - Energy.v)/(dE*Tau[0]);
      } while ((G < this->getSigma() || Energy.v < energyNew || aol::isNaN( Energy.v )) && (Tau[0] > this->getTauMin()));
      energyBackup = energyNew;
    }//else

    if( Tau[0] <= this->getTauMin() )
      Tau[0] = 0.;
    else {
      DescentDir *= Tau[0];
      CurrentPosition += DescentDir;
      this->_F.apply( CurrentPosition, Derivative );
      Energy.v = energyBackup;
    }
  }

  // auxiliary method for getTauAndUpdateDescentDirWolfe
  virtual void zoomWolfe( DerivativeType &DescentDir,
                          VectorType &CurrentPosition,
                          DerivativeType &Derivative,
                          aol::Scalar<RealType> &Energy,
                          const RealType dE,
                          RealType ELo,
                          RealType TauLo,
                          RealType TauHi,
                          aol::Vector<RealType> &Tau ) const {
    // either: TauLo satisfies the Armijo condition, but is too short for curvature condition, and TauHi > TauLo does not satisfy Armijo condition
    // or: both satisfy Armijo condition; TauLo is too long and TauHi too short for curvature condition

    RealType dENew;
    aol::Scalar<RealType> energyNew;
    VectorType positionNew( CurrentPosition, aol::STRUCT_COPY );

    bool stepSizeFound = false;
    int counter = 0; // for non-termination due to rounding errors, we should fix the maximum number of iterations
    do {
      Tau[0] = ( TauLo + TauHi ) / 2.;
      positionNew = CurrentPosition;
      positionNew.addMultiple( DescentDir, Tau[0] );
      this->_E.apply( positionNew, energyNew );

      if ( energyNew.v > aol::Min( Energy.v + this->getSigma() * Tau[0] * dE, ELo ) || aol::isNaN( energyNew.v ) )
        // Armijo condition violated...
        TauHi = Tau[0];

      else {
        // Armijo condition fulfilled...
        this->_F.apply( positionNew, Derivative );
        dENew = Derivative * DescentDir;

        if ( aol::Abs( dENew ) <= - this->getBeta() * dE )
          // strong Wolfe (Armijo + curvature) condition fulfilled...
          stepSizeFound = true;

        else {
          // curvature condition violated...
          if ( dENew * ( TauHi - TauLo ) >= 0. )
            TauHi = TauLo;
          TauLo = Tau[0];
          ELo = energyNew;
        }
      }
    } while ( !stepSizeFound && ++counter < 30 );

    if ( counter >= 30 ) {
      Tau[0] = 0.;
      this->_F.apply( CurrentPosition, Derivative );
    }
    else {
      CurrentPosition = positionNew;
      Energy = energyNew;
    }
    DescentDir *= Tau[0];
  }

  // CurrentPosition contains the currently optimal point,
  // DescentDir the linesearch direction,
  // Derivative the energy derivative at CurrentPosition.
  // All parameters at the beginning contain the values at the initial position,
  // and their new values at the end.
  // Algorithm taken from: Nocedal, Wright: Numerical Optimization, Alg.3.5-3.6
  virtual void getTauAndUpdateDescentDirWolfe( DerivativeType &DescentDir,
                                               VectorType &CurrentPosition,
                                               DerivativeType &Derivative,
                                               aol::Scalar<RealType> &Energy,
                                               aol::Vector<RealType> &Tau ) const {

    Tau[0] = 1.; // important for Newton-based methods for superlinear convergence, since the first steplength is chosen which satisfies the Wolfe conditions
    RealType dE = Derivative * DescentDir, dENew, energyBackup = Energy.v, TauOld = 0.;
    aol::Scalar<RealType> energyNew;
    VectorType positionNew( CurrentPosition, aol::STRUCT_COPY );

    bool stepSizeFound = false;
    do {
      positionNew = CurrentPosition;
      positionNew.addMultiple( DescentDir, Tau[0] );
      this->_E.apply( positionNew, energyNew );

      if ( energyNew.v > aol::Min( Energy.v + this->getSigma() * Tau[0] * dE, energyBackup ) || aol::isNaN( energyNew.v ) ) {
        // Armijo condition violated...
        zoomWolfe( DescentDir, CurrentPosition, Derivative, Energy, dE, energyBackup, TauOld, Tau[0], Tau );
        stepSizeFound = true;

      } else {
        // Armijo condition fulfilled...
        this->_F.apply( positionNew, Derivative );
        dENew = Derivative * DescentDir;

        if ( aol::Abs( dENew ) <= - this->getBeta() * dE ) {
          // strong Wolfe condition (Armijo condition + curvature condition) fulfilled...
          DescentDir *= Tau[0];
          CurrentPosition = positionNew;
          Energy = energyNew;
          stepSizeFound = true;

        } else {
          // curvature condition violated...
          if ( dENew >= 0. ) {
            // too long step
            zoomWolfe( DescentDir, CurrentPosition, Derivative, Energy, dE, Energy.v, Tau[0], TauOld, Tau );
            stepSizeFound = true;

          } else {
            // too short step
            TauOld = Tau[0];
            energyBackup = energyNew.v;
            Tau[0] *= 2.;
          }
        }
      }
    } while ( !stepSizeFound );
  }

  virtual void writeResult( const VectorType &Dest, const char* FileName ) const {
    Dest.saveASCII( FileName );
  }

public:
  /**
   * Returns the result of the quasi-Newton iteration in Dest, while Arg is the starting value.
   */
  virtual void apply( const VectorType &Arg, VectorType &Dest ) const {

    aol::Scalar<RealType> e;
    VectorType Dx ( Arg, aol::DEEP_COPY );
    DerivativeType f ( Dest, aol::STRUCT_COPY );
    DerivativeType Df ( Dest, aol::STRUCT_COPY );
    DerivativeType tmp ( Dest, aol::STRUCT_COPY );
    DerivativeType descentDir ( Dest, aol::STRUCT_COPY );

    // initialization
    Dest = Arg;
    this->_E.apply( Dest, e );
    this->_F.apply( Dest, f );
    RealType FNorm = this->fNorm( f );
    RealType FNormNew = aol::NumberTrait<RealType>::one;
    aol::Vector<RealType> tau(1);
    tau.setAll( aol::NumberTrait<RealType>::one );

    this->_pInfo->startIterations();

    while ( FNorm > this->_pInfo->getAccuracy() && !(this->_pInfo->maxIterIsReached()) && ( tau.getMaxValue() > 0. || _c.size() > 0 ) ) {
      this->_pInfo->startStep();
      // Quasi-Newton-iteration given by x^{k+1} = x^k - tau^k H^k f(x^k)
      // iteration number k (i.e. find x^{k+1} from current approximation x^k)

      // compute Dx = x^k - x^{k-1} and Df = f(x^k) - f(x^{k-1})
      Dx -= Dest;
      Dx *= -1.;
      Df -= f;
      Df *= -1.;
      // update Z and c such that H^k = \sum_i c_iZ_iZ_i^T,
      // where H^0 = I and H^{k+1} = H^k + (Dx^k\otimes Dx^k)/(Df^k\cdot Dx^k) - ((H^k Df^k)\otimes(H^k Df^k))/(Df^k\cdot H^k Df^k) (DFP-method)
      if ( ( ( this->_pInfo->getIterationCount() - 1 ) % _resetNum > 0 ) && ( tau.getMaxValue() > 0. ) ) {
        DerivativeType HDf( Df, aol::DEEP_COPY );
        for ( unsigned int i = 0; i < _c.size(); i++ )
          HDf.addMultiple( _Z[i], _c[i] * ( _Z[i] * Df ) );
        _Z.pushBack( HDf );
        _c.push_back( - 1. / ( HDf * Df ) );
        _Z.pushBack( Dx );
        _c.push_back( 1. / ( Df * Dx ) );
      } else {
        _Z.clear();
        _c.clear();
        tau.setAll( aol::NumberTrait<RealType>::one );
      }
      // compute -H^k f(x^k) = -f(x^k) - \sum_i c_iZ_iZ_i^Tf(x^k)
      descentDir = f;
      descentDir *= -1;
      for ( unsigned int i = 0; i < _c.size(); i++ )
        descentDir.addMultiple( _Z[i], -_c[i] * ( _Z[i] * f ) );
      // save x and f for next step
      Dx = Dest;
      Df = f;

      this->smoothDirection( Dest, descentDir );
      RealType normUpdate = descentDir.norm();

      // update x^k, f(x^k), and its norm
      getTauAndUpdateDescentDir( descentDir, Dest, f, e, tau );
      FNormNew = this->fNorm(f);

#ifdef VERBOSE
      cerr << "FNorm before iteration " << this->_pInfo->getIterationCount() << " : " << FNorm << endl;
      cerr << "FNorm after iteration " << this->_pInfo->getIterationCount() << "  : " << FNormNew << endl;
#endif

      FNorm = FNormNew;

      this->postProcess( Dest, this->_pInfo->getIterationCount() );

      // output
      if( this->checkSaveConditions( this->_pInfo->getIterationCount()) ){
        string fn;
        fn = aol::strprintf( "%s%s_newtonstep_%03d_%.5f", this->getSaveDirectory(), this->getSaveName(), this->_pInfo->getIterationCount(), FNormNew );
        writeResult( Dest, fn.c_str() );
      }

      this->_pInfo->finishStep( FNormNew, normUpdate, tau.getMinValue() );
    } // end while

    this->_pInfo->finishIterations();
    cerr << endl;
  }

  virtual void applyAdd( const VectorType &/*Arg*/, VectorType &/*Dest*/ ) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }

  void reset ( ) const {
    cerr<<"\n reset z and c in Quasi Newton Method\n";
    _Z.clear();
    _c.clear();
  }
};

/**
 * \brief This class implements a quasi-Newton method (the DFP-method) for minimizing a scalar function.
 *        The scalar function is a function of a MultiVector, whose components are updated using separate step lengths.
 *
 * It has to be provided the scalar function and its derivative.
 *
 * \author Wirth
 * \ingroup Optimization
 */
template < typename RealType >
class QuasiNewtonIterationComponentWiseTimestepControlled
  : public aol::QuasiNewtonIteration< RealType, MultiVector<RealType>, MultiVector<RealType> > {

protected:
  //! components of the MultiVector apply argument can be grouped to share the same tau with _customTauBlocks
  aol::Vector<int> _customTauBlocks;

public:
  QuasiNewtonIterationComponentWiseTimestepControlled( const aol::Op<MultiVector<RealType>, aol::Scalar<RealType> > &E,
                                                       const aol::Op<MultiVector<RealType>, MultiVector<RealType> > &F,
                                                       const aol::Vector<int> &TauBlocks,
                                                       const int MaxIterations = 50,
                                                       const RealType StopEpsilon = 1.e-6,
                                                       const int ResetNum = 10,
                                                       const bool WriteTimeSteps = false,
                                                       const char *BaseSaveName = NULL ) :
    aol::QuasiNewtonIteration< RealType, MultiVector<RealType>, MultiVector<RealType> >( E, F, MaxIterations, StopEpsilon, ResetNum, WriteTimeSteps, BaseSaveName),
    _customTauBlocks( TauBlocks ) {}

  QuasiNewtonIterationComponentWiseTimestepControlled( const aol::Op<MultiVector<RealType>, aol::Scalar<RealType> > &E,
                                                       const aol::Op<MultiVector<RealType>, MultiVector<RealType> > &F,
                                                       const aol::Vector<int> &TauBlocks,
                                                       aol::NewtonInfo<RealType> &Info,
                                                       const int ResetNum = 10,
                                                       const bool WriteTimeSteps = false,
                                                       const char *BaseSaveName = NULL ) :
    aol::QuasiNewtonIteration< RealType, MultiVector<RealType> , MultiVector<RealType> >( E, F, Info, ResetNum, WriteTimeSteps, BaseSaveName),
    _customTauBlocks( TauBlocks ) {}

  virtual ~QuasiNewtonIterationComponentWiseTimestepControlled () {}

protected:
  virtual void getTauAndUpdateDescentDir( aol::MultiVector<RealType> &DescentDir,
                                          aol::MultiVector<RealType> &CurrentPosition,
                                          aol::MultiVector<RealType> &Derivative,
                                          aol::Scalar<RealType> &/*Energy*/,
                                          aol::Vector<RealType> &Tau ) const {
    aol::Vector<int> tauBlocks(0);
    if ( _customTauBlocks.size() ) {
      tauBlocks.resize( _customTauBlocks.size() );
      tauBlocks = _customTauBlocks;
    } else {
      tauBlocks.resize( DescentDir.numComponents() );
      tauBlocks.setAll( 1 );
    }
    const int numComponents = tauBlocks.size();
    if( numComponents != Tau.size() ) {
      Tau.resize( numComponents );
      Tau.setAll( 1. );
    }
    aol::MultiVector<RealType> descentDirComponent( DescentDir, aol::STRUCT_COPY );

    int handledDDComponents = 0;
    for( int component = 0; component < numComponents; component++){
      for( int i = 0; i < tauBlocks[component]; i++)
        descentDirComponent[i+handledDDComponents] = DescentDir[i+handledDDComponents];
      Tau[component] = this->getStepsize(descentDirComponent, CurrentPosition, Tau[component]);
      for( int i = 0; i < tauBlocks[component]; i++){
        DescentDir[i+handledDDComponents] *= Tau[component];
        descentDirComponent[i+handledDDComponents].setZero();
      }
      handledDDComponents += tauBlocks[component];
    }

    CurrentPosition += DescentDir;
    this->_F.apply( CurrentPosition, Derivative );
  }
};

template < typename RealType, typename VectorType>
class BFGSOp {
protected:
  int _counter;
  const int _reset;
  aol::RandomAccessContainer<VectorType> _y;
  aol::RandomAccessContainer<VectorType> _dx;
  aol::Vector<RealType> _dx_y;
  mutable aol::RandomAccessContainer<VectorType> _auxMemory;

public:
  BFGSOp (const int Reset) :
    _counter ( 0 ),
    _reset ( Reset ),
    _dx_y( Reset ) {}

  void updateOpBFGS ( const VectorType &DX,
                      const VectorType &Y ){
    if ( _counter == _reset )
      this->reset();
    _counter++;
    
    if ( _y.size() < _counter ) {
      _y.pushBack ( Y );
      _dx.pushBack( DX );
      _auxMemory.pushBack( Y );
    } else {
      _y[_counter-1] = Y;
      _dx[_counter-1] = DX;
    }
    _dx_y[_counter-1] = DX*Y;
  }

  void reset () {
    _y.clear();
    _dx.clear();
    _auxMemory.clear();
    _counter = 0;
  }

  void apply ( const VectorType &Arg, VectorType &Dest ) const {
    Dest = Arg;
    for ( int i = 0; i < _counter; ++i ) {
      _auxMemory[i] = _dx[i]; // _auxMemory[k] = B_k \Delta dx_k
      for ( int j = 0; j < i; ++j ) {
        _auxMemory[i].addMultiple( _y[j], ( _y[j] * _dx[i] ) / _dx_y[j] );
        _auxMemory[i].addMultiple( _auxMemory[j], - ( _auxMemory[j] * _dx[i] ) / ( _auxMemory[j] * _dx[j] ) );
      }
      Dest.addMultiple( _y[i], ( _y[i] * Arg ) / _dx_y[i] );
      Dest.addMultiple( _auxMemory[i], - ( _auxMemory[i] * Arg ) / ( _auxMemory[i] * _dx[i] ) );
    }
  }

  void applyInverse ( const VectorType &Arg, VectorType &Dest ) const {
    Dest = Arg;
    for ( int i = _counter - 1; i >= 0; --i ) {
      RealType dotProd = ( _dx[i] * Dest ) / _dx_y[i];
      _auxMemory[i] = _dx[i];
      _auxMemory[i] *= dotProd; // _auxMemory[k] = ( dx_k dx_k^T / (y_k^T dx_k)) (I - y_{k+1} dx_{k+1}^T / (y_{k+1}^T dx_{k+1})) ... (I - y_K dx_K^T / (y_K^T dx_K)) Arg
      Dest.addMultiple( _y[i], -dotProd ); // Dest = (I - y_k dx_k^T / (y_k^T dx_k)) ... (I - y_K dx_K^T / (y_K^T dx_K)) Arg
    }
    for ( int i = 0; i < _counter; ++i ) {
      Dest.addMultiple( _dx[i], -( _y[i] * Dest ) / _dx_y[i] ); // aux = (I - y_k dx_k^T / (y_k^T dx_k))^T H_{k-1} (I - y_k dx_k^T / (y_k^T dx_k)) ... (I - y_K dx_K^T / (y_K^T dx_K)) Arg
      Dest += _auxMemory[i];
    }
  }
};

/**
 * \brief This class implements a quasi-Newton method with BFGS-Update of the matrix for minimizing a scalar function.
 *
 * It has to be provided the scalar function and its derivative.
 *
 * \author Teusner
 * \ingroup Optimization
 */
template < typename RealType, typename VectorType, typename DerivativeType>
class QuasiNewtonBFGS
  : public aol::NewtonMinimizationBase< RealType, VectorType, DerivativeType, aol::IdentityOp<VectorType> > {
protected:

  const aol::IdentityOp<VectorType> _identity;
  const int _reset;
  mutable BFGSOp<RealType, VectorType> _BFGSOp;
  bool _energyDecayBasedStopping;
  const StepSaverBase<RealType, VectorType> *_pStepSaver;
  

public:
  QuasiNewtonBFGS ( const aol::Op<VectorType, aol::Scalar<RealType> > &E,
                    const aol::Op<VectorType, DerivativeType> &F,
                    const int MaxIterations = 50,
                    const RealType StopEpsilon = 1.e-6,
                    const int Reset = 10,
                    const bool WriteTimeSteps = false,
                    const char *BaseSaveName = NULL ) :
    aol::NewtonMinimizationBase< RealType, VectorType, DerivativeType, aol::IdentityOp<VectorType> >( E, F, aol::NullOp<VectorType, aol::IdentityOp<VectorType> >(), NULL, MaxIterations, StopEpsilon, WriteTimeSteps, BaseSaveName),
    _identity ( ),
    _reset ( Reset ),
    _BFGSOp ( Reset ),
    _energyDecayBasedStopping ( false ),
    _pStepSaver ( NULL ) {
    this->_pInfo->getSolverInfo().setQuietMode();
  }

  QuasiNewtonBFGS ( const aol::Op<VectorType, aol::Scalar<RealType> > &E,
                    const aol::Op<VectorType, DerivativeType> &F,
                    aol::NewtonInfo<RealType> & Info,
                    const int Reset = 100,
                    const bool WriteTimeSteps = false,
                    const char *BaseSaveName = NULL ) :
    aol::NewtonMinimizationBase< RealType, VectorType, DerivativeType, aol::IdentityOp<VectorType> >( E, F, aol::NullOp<VectorType, aol::IdentityOp<VectorType> >(), NULL, Info, WriteTimeSteps, BaseSaveName),
    _identity ( ),
    _reset ( Reset ),
    _BFGSOp ( Reset ),
    _energyDecayBasedStopping ( false ),
    _pStepSaver ( NULL ) {
    this->_pInfo->getSolverInfo().setQuietMode();
  }

  virtual ~QuasiNewtonBFGS () {
  }

public:
  /**
   * Returns the result of the quasi-Newton iteration in Dest, while Arg is the starting value.
   */
  virtual void apply( const VectorType &Arg, VectorType &Dest ) const {
    Dest = Arg;
    applySingle ( Dest );
  }
  virtual void applySingle( VectorType &Dest ) const {
    VectorType Dx ( Dest, aol::STRUCT_COPY );
    DerivativeType f ( Dest, aol::STRUCT_COPY );
    DerivativeType Df ( Dest, aol::STRUCT_COPY );
    DerivativeType tmp ( Dest, aol::STRUCT_COPY );
    DerivativeType descentDir ( Dest, aol::STRUCT_COPY );

    RealType FNorm = aol::NumberTrait<RealType>::Inf;
    RealType FNormNew = 1.0;
    RealType tau = 1.0;

    RealType E = aol::NumberTrait<RealType>::Inf;
    RealType ENew = aol::NumberTrait<RealType>::Inf;

    this->_pInfo->startIterations();

    bool forcedReset = false;

    if ( _energyDecayBasedStopping && ( this->_pInfo->getTimestepController() != aol::NewtonInfo<RealType>::ARMIJO ) )
      throw aol::UnimplementedCodeException ( "Energy based stopping so far only works with the Armijo rule.", __FILE__, __LINE__ );

    // Save the current time step if _pStepSaver was initialized with setStepSaverReference
    if ( _pStepSaver )
      _pStepSaver->saveStep ( Dest, this->_pInfo->getIterationCount ( ) );

    while ( FNorm > this->_pInfo->getAccuracy() && !(this->_pInfo->maxIterIsReached()) && tau > 0. ) {
      // Quasi-Newton-iteration given by x_{k+1} = x_k - tau_k B_k^{-1} f(x_k)
      // iteration number k (i.e. find x_{k+1} from current approximation x_k)

      // compute f and its norm
      this->_F.apply( Dest, f );

      // save x_k and f(x_k) for the computation of Dx = x_{k+1}-x_k and Df = f(x_{k+1})-f(x_k)
      Dx = Dest;
      Dx *= -1.;
      Df = f;
      Df *= -1.;

      _BFGSOp.applyInverse( f, descentDir );
      descentDir *= -1.;

      this->smoothDirection( Dest, descentDir );

      switch( this->_pInfo->getTimestepController() ) {
      case aol::NewtonInfo<RealType>::ARMIJO: tau = this->getTimestepWidthWithArmijoLineSearch( descentDir, Dest, tau );
        break;
      case aol::NewtonInfo<RealType>::NEWTON_OPTIMAL: throw Exception( "aol::NewtonInfo<RealType>::NEWTON_OPTIMAL is not implemented!", __FILE__, __LINE__);
        break;
      case aol::NewtonInfo<RealType>::NEWTON_SIMPLE:  tau = this->getSimpleNewtonTimestepWidth( descentDir, Dest );
        break;
      case aol::NewtonInfo<RealType>::NO_TIMESTEP_CONTROL: tau = aol::NumberTrait<RealType>::one;
        break;
      case aol::NewtonInfo<RealType>::WOLFE: tau = this->getTimestepWidthWithPowellWolfeLineSearch( descentDir, Dest, tau );
        break;
      default:
        throw ( aol::Exception ( "aol::NewtonInfo<RealType>: Illegal TimestepController", __FILE__, __LINE__ ) );
      }

      if ( ( tau == 0 ) && ( forcedReset == false ) ) {
        _BFGSOp.reset();
        tau = aol::NumberTrait<RealType>::one;
        ENew = aol::NumberTrait<RealType>::Inf;
        forcedReset = true;
        continue;
      }
      forcedReset = false;

      // In the function startStep the counter of steps is raised by one. If we would call this function at the beginning of the
      // while loop we would raise up the counter twice in one step with forced reset.
      this->_pInfo->startStep();

      RealType normUpdate = descentDir.norm();
      descentDir *= tau;
      Dest += descentDir;

      this->_F.apply(Dest, tmp);
      FNormNew = this->fNorm(tmp);

      Dx += Dest;
      Df += tmp;

      // update of B, B_{k+1} = B_k - \frac{B_k Dx Dx^T B_k}{Dx\cdot (B_k Dx)}+\frac{Df Df^T}{Df \cdot Dx}
      _BFGSOp.updateOpBFGS(Dx, Df);

#ifdef VERBOSE
      cerr << "FNorm before iteration " << this->_pInfo->getIterationCount() << " : " << FNorm << endl;
      cerr << "FNorm after iteration " << this->_pInfo->getIterationCount() << "  : " << FNormNew << endl;
#endif

      FNorm = FNormNew;

      this->postProcess( Dest, this->_pInfo->getIterationCount() );

      this->_pInfo->finishStep(FNormNew, normUpdate, tau);

      // Save the current time step if _pStepSaver was initialized with setStepSaverReference
      if ( _pStepSaver )
        _pStepSaver->saveStep ( Dest, this->_pInfo->getIterationCount ( ) );

      // output (this is the old method for saving time steps, which is kept here for compatibility reasons)
      if( this->checkSaveConditions( this->_pInfo->getIterationCount()) ){
        string fn;
        fn = aol::strprintf( "%s%s_newtonstep_%03d_%.5f", this->getSaveDirectory(), this->getSaveName(), this->_pInfo->getIterationCount(), FNormNew);
        this->writeResult( Dest, fn.c_str() );
      }

      E = ENew;
      ENew = this->_energyAtPosition;

      if ( _energyDecayBasedStopping && ( ( E - ENew ) < this->_pInfo->getAccuracy() ) )
        break;
    } // end while

    this->_pInfo->finishIterations();
    if ( !this->_pInfo->getQuietMode ( ) )
      cerr << endl;
  }

  virtual void applyAdd( const VectorType &/*Arg*/, VectorType &/*Dest*/ ) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }

  virtual void prepareSolver() const {}

  void setEnergyDecayBasedStopping ( const bool EnergyDecayBasedStopping ) {
    _energyDecayBasedStopping = EnergyDecayBasedStopping;
  }

  void setStepSaverReference ( const StepSaverBase<RealType, VectorType> &StepSaver ) {
    _pStepSaver = &StepSaver;
  }
};


} // namespace aol

#endif // __NEWTON_H
