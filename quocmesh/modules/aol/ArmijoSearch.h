#ifndef __ARMIJOSEARCH_H
#define __ARMIJOSEARCH_H

#include <scalarArray.h>
#include <fastUniformGridMatrix.h>
#include <solver.h>
#include <preconditioner.h>
#include <linearSmoothOp.h>
#include <mcm.h>
#include <CSRMatrix.h>
#include <timestepSaver.h>
#include <NewtonInfo.h>
#include <vec.h>

namespace aol {

/** @brief The interface of searching based on Armijo rule.
 *
 */
template <typename RealType, typename DomainType, typename RangeType=DomainType, typename DerivativeType=DomainType>
    class ArmijoLineSearchUsingOp : public aol::Op<DomainType, RangeType >, public virtual aol::TimestepSaver<RealType> {
private:
  /** @brief The tolerence parameters.*/
  RealType _sigma;
  //! Additional parameter for getTimestepWidthWithPowellWolfeLineSearch.
  RealType _beta;
protected:
  mutable RealType _startTau;
  RealType _tauMin, _tauMax;
private:
  /** @brief Get the function evauluation at new position. (Needs to be implemented)
  */
  virtual RealType ArmijoLineSearchHelpFunction_evaluate( const DerivativeType &DescentDir,
                                                          const DomainType &CurrentPosition,
                                                          const RealType timestepWidth ) const = 0;

  /** @brief Get the derivative of function at new position. (Needs to be implemented)
  */
  virtual RealType ArmijoLineSearchHelpFunction_evaluateDerivative( const DerivativeType &DescentDir,
                                                                    const DomainType &Position ) const = 0;

  /** @brief Get the derivative of function at new position. (Needs to be implemented)
  */
  virtual RealType ArmijoLineSearchHelpFunction_evaluateDerivativeWithTau( const DerivativeType &/*DescentDir*/,
                                                                           const DomainType &/*Position*/,
                                                                           const RealType /*timestepWidth*/ ) const {
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }

protected:
  /** @brief save the timestep
      @param Data              The data to be saved
      @param Iteration         The iteration (oh really?)
  */
  virtual void writeTimeStep ( const DomainType &/*Data*/, const int /*Iteration*/ ) const
  {
    cerr << "Please overload writeTimeStep." << endl;
  }


  /** @brief The initialization constructor
      @param  Sigma  The tolerance parameter.
  */
  ArmijoLineSearchUsingOp( const RealType Sigma,
                           const int MaxIterations,
                           const RealType StartTau = aol::ZOTrait<RealType>::one,
                           const RealType TauMax = pow(2., 25.),
                           const RealType Beta = 0.9 )
  : _sigma(Sigma),
    _beta(Beta),
    _startTau(StartTau),
    _tauMin( pow(2.,-30.) ),
    _tauMax( TauMax )
    {
      // Just for pretty console output.
      const int tempInt = MaxIterations / 100;
      if (tempInt > 0)
        this->setSaveTimestepOffset( tempInt );
      this->setQuietMode( true );
      this->setStepDigits( aol::countDigitsOfNumber ( MaxIterations ) );

      this->setWriteTimeSteps( false );
    }

  /** @brief Computation of time step by Armijo rule
      @param DescentDir          The direction of descent.
      @param CurrentPosition     The current point of evulation.
      @param OldTau              The previous time step.
  */
  RealType getTimestepWidthWithSimpleLineSearch( const DerivativeType &DescentDir,
                                                 const DomainType &CurrentPosition,
                                                 const RealType OldTau = aol::ZOTrait<RealType>::one ) const{

    // extreme simple timestep width control, just ensures, that fnew < fnew

    RealType tau = OldTau;

    RealType f = ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, 0. );
    RealType fNew = ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, tau );

    while( (fNew >= f) && (tau >= getTauMin() ) ){
        tau = tau * 0.5;
        fNew = ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, tau );
    }

    // No energy descent for tau >= tauMin, so we don't want the step to be done.
    // The stopping criterion also handles this case.
    if ( tau < getTauMin() )
      tau = 0;

#ifdef VERBOSE
    cerr << "tau = " << tau << endl;
#endif

    return tau;
  }

  RealType getTimestepWidthWithArmijoLineSearch( const DerivativeType &DescentDir,
                                                 const DomainType &CurrentPosition,
                                                 const RealType OldTau) const{

#ifdef VERBOSE
    cerr << endl << aol::color::brown << "Calculating tau with Armijo Line Search (starting with tau^0 = " << OldTau << ")" << endl << endl;
#endif
    // Makes sure that the starting value for tau is between tauMin and tauMax.
    // The case OldTau == 0. happens, if one of the vector components of the
    // descent direction in GradientDescent in the last step didn't lead to an
    // energy reduction and the step is not done by setting tau to zero.
    RealType tau = aol::Clamp( OldTau, getTauMin(), getTauMax() );

    const RealType f = ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, 0. );
    RealType fNew = ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, tau );
    const RealType Df = ArmijoLineSearchHelpFunction_evaluateDerivative( DescentDir, CurrentPosition );
    RealType G = (fNew - f)/(Df*tau);

#ifdef VERBOSE
    cerr << "f = f(x)              = " << aol::detailedFormat ( f ) << endl
         << "f_new = f(x + tau p)  = " << aol::detailedFormat ( fNew ) << endl
         << "gradient slope Df     = " << aol::detailedFormat ( Df ) << endl
         << "G = secant slope / Df = " << aol::detailedFormat ( G ) << endl << endl;
#endif

    if ( Df >= 0 )
      return 0;

    if ( isFinite ( f ) == false ) {
#ifdef VERBOSE
      cerr << "Error: Initial f is not finite. Rejecting step." << endl;
#endif
      return 0;
    }

    //check
    if( G >= _sigma && f >= fNew){
      //time step too small
#ifdef VERBOSE
      cerr << "initial stepsize too small. Will increase tau ..." << endl;
#endif
      do{
        tau *= 2.;
        fNew = ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, tau );
        G = (fNew - f)/(Df*tau);
#ifdef VERBOSE
        cerr << "tau = " << aol::detailedFormat ( tau ) << ",\tf_new = " << aol::detailedFormat ( fNew ) << ",\tG = " << aol::detailedFormat ( G ) << endl;
#endif
      } while (G >= _sigma && f >= fNew && tau <= getTauMax());
      tau *= 0.5;
    }// if
    else{
      // time step too large
#ifdef VERBOSE
      cerr << "initial stepsize too large. Will decrease tau ..." << endl;
#endif
      do{
        if (tau > getTauMin())
          tau *= 0.5;
        fNew = ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, tau );
        G = (fNew - f)/(Df*tau);
#ifdef VERBOSE
        cerr << "tau = " << aol::detailedFormat ( tau ) << ",\tf_new = " << aol::detailedFormat ( fNew ) << ",\tG = " << aol::detailedFormat ( G ) << endl;
#endif
      } while ( ( ( isFinite(fNew) == false ) || (G < _sigma || f < fNew) ) && (tau > getTauMin()));
    }//else

    // If tau == tauMin it's very likely that the condition for the step size control
    // is not satisfied. In such a case we don't want the step to be done. The stopping
    // criterion also handles this case.
    if( tau > getTauMin() ) {
#ifdef VERBOSE
      cerr << "accepted tau = " << aol::detailedFormat ( tau ) << endl << endl << aol::color::reset;
#endif
      return tau;
    }
    else {
#ifdef VERBOSE
      cerr << "tau <= tauMin, returning 0" << endl << endl << aol::color::reset;
#endif
      return 0.;
    }
  }

  /** The algorithm was found in the book
  *  Kosmol, Methoden zur numerischen Bahandlung nichtlinearere Gleichungen und Optimierungsaufgaben,
  *  B.G. Teubner Stuttgart, 1989
  */
  RealType getTimestepWidthWithPowellWolfeLineSearch( const DerivativeType &DescentDir,
                                                      const DomainType &CurrentPosition,
                                                      const RealType OldTau) const{
  #ifdef VERBOSE
      cerr << "Calculating tau with Powell-Wolfe Line Search (starting with tau = " << OldTau << ")\n";
  #endif
      // Makes sure that the starting value for tau is between tauMin and tauMax.
      // The case OldTau == 0. happens, if one of the vector components of the
      // descent direction in GradientDescent in the last step didn't lead to an
      // energy reduction and the step is not done by setting tau to zero.
      RealType tau = aol::Min( aol::Max( OldTau, getTauMin() ), getTauMax() );
      RealType s = tau;

      RealType alpha = 0.5;
      RealType r = 0;
      RealType p = 0;

      RealType f = ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, 0. );
      RealType fNew = ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, tau );
      RealType Df = ArmijoLineSearchHelpFunction_evaluateDerivative( DescentDir, CurrentPosition );
      RealType DfNew = ArmijoLineSearchHelpFunction_evaluateDerivativeWithTau( DescentDir, CurrentPosition, s );
      RealType G = (fNew - f)/(Df*tau);
      RealType P = DfNew/Df;

//stepOne:
      if ( G < _sigma ) {
        goto stepFour;
      }

stepTwo:
      DfNew = ArmijoLineSearchHelpFunction_evaluateDerivativeWithTau( DescentDir, CurrentPosition, s );
      P = DfNew/Df;

      if (P <= _beta ) {
        tau = s;
        goto stepEight;
      }
      r = s;

//stepThree:
      // enlagement of step size
      s = r/alpha;
      fNew = ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, s );
      G = (fNew - f)/(Df*s);

      if ( ( G >= _sigma ) && ( isFinite ( fNew ) ) ) {
        goto stepTwo;
      } else {
        goto stepSix;
      }

stepFour:
      // reduction of step size
      r = s * alpha;
      if ( r < getTauMin() ) {
        goto stepEight;
      }
      fNew = ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, r );
      G = (fNew - f)/(Df*r);

      if ( G < _sigma ) {
        s = r;
        goto stepFour;
      }

stepFive:
      DfNew = ArmijoLineSearchHelpFunction_evaluateDerivativeWithTau( DescentDir, CurrentPosition, r );
      P = DfNew/Df;

      if ( P <= _beta ) {
        tau = r;
        goto stepEight;
      }

stepSix:
      if ( std::fabs(r-s) < 1e-16 ) {
        tau = 0;
        goto stepEight;
      }
      p = (r + s)/2.;
      if ( p < getTauMin() ) {
        goto stepEight;
      }
//stepSeven:
      fNew = ArmijoLineSearchHelpFunction_evaluate( DescentDir, CurrentPosition, p );
      G = (fNew - f)/(Df*p);

      if ( G >= _sigma ) {
        r = p;
        goto stepFive;
      } else {
        s = p;
        goto stepSix;
      }

stepEight:
      if( tau > getTauMin() )
        return tau;
      else
        return 0.;
    }
public:
  RealType getStartTau() const{
    return _startTau;
  }

  //! \short CAVE: This function is overloaded in NewtonIterationBase.
  //! All Newton classes do not access this variable, but the one
  //! from their NewtonInfo object.
  virtual RealType getTauMin () const {  return _tauMin;  }
  virtual RealType getTauMax () const {  return _tauMax;  }

  virtual void setTauMin( const RealType TauMin ) {
    if( TauMin < numeric_limits<RealType>::epsilon( ) )
      throw aol::Exception( "Trying to set TauMin out of range", __FILE__, __LINE__);

    _tauMin = TauMin;
  }
  virtual void setTauMax( const RealType TauMax ) { _tauMax = TauMax; }

  void setSigma(const RealType Sigma){
    if( (Sigma <= 0.) || ( Sigma >= 1. ) )
      throw aol::Exception( "Trying to set sigma out of range", __FILE__, __LINE__);
    else
      _sigma = Sigma;
  }
  RealType getSigma()const{
    return _sigma;
  }
  void setBeta(const RealType Beta){
    if( (Beta <= _sigma) || ( Beta >= 1. ) )
      throw aol::Exception( "Trying to set beta out of range", __FILE__, __LINE__);
    else
      _beta = Beta;
  }
  RealType getBeta()const{
    return _beta;
  }
};

} // namespace aol

#endif
