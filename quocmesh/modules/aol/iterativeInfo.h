#ifndef __ITERATIVEINFO_H
#define __ITERATIVEINFO_H

#include <convergenceEstimator.h>

namespace aol {

/****************************************************************************
 *
 *        CLASS IterativeInfo
 */
/**
 *  \brief Base class for SolverInfo and NewtonInfo.
 *
 *  \author von Deylen
 */
template <typename _RealType>
class IterativeInfo {
public:
  typedef _RealType RealType;

  IterativeInfo ( RealType TOL, int maxIter, bool quietMode, ostream & out,
                  int numAveragingSteps, int estimatedConvergenceOrder );
  virtual ~IterativeInfo ();

  // parameter setting routines
  void setAccuracy ( RealType TOL );
  void setMaxIterations ( int maxIter );
  virtual void setQuietMode ( bool quiet = true );
  virtual void setMegaQuietMode ( bool quiet = true );
  virtual bool getMegaQuietMode ( );
  void setNumAveragingSteps ( int n );

  // parameter getting routines
  RealType getAccuracy() const;
  int      getMaxIterations() const;
  bool     getQuietMode() const;
  int      getNumAveragingSteps () const;

  // begin and end of iteration (steps)
  virtual void startIterations();
  virtual void startStep();
  // void finishStep() is protected
  virtual void finishStep ( RealType residual );
  virtual void finishIterations();

  virtual void setCurrentResidual ( RealType residual );
  RealType getCurrentResidual() const;
  void setTablePrintout ( bool tablePrintout = true );

  // current iteration infos
  void printStats() const;
  void printStats ( ostream & out ) const;

  RealType getFinalResidual() const;
  int getIterationCount() const;
  virtual bool maxIterIsExceeded() const;
  virtual bool maxIterIsReached() const;
  bool currentResidualIsNaN() const;

  //!copy iteration count and residuals
  void getIterationInfoFrom ( const IterativeInfo<RealType> & info );

  RealType getLastConvergenceCoeff ( int order ) const;
  RealType getLastConvergenceOrder() const;
  RealType getConvergenceCoeff ( int order ) const;
  RealType getConvergenceOrder() const;

  ostream & getOstream() const;

protected:
  RealType  _TOL;
  int       _maxIter;
  bool      _quiet;
  bool      _megaQuiet;
  ostream & _out;

  RealType _finalResidual;

  //! number of started steps
  int       _stepCount;

  //! is true as long as startStep() was not successed by finishStep().
  bool      _currentlyProceedingStep;
  bool      _residualWasSet;

  ConvergenceEstimator<RealType> _convergenceEstimator;

  virtual void finishStep();

  void ( IterativeInfo<RealType>::* _logBegin ) ( ostream & out ) const;
  void ( IterativeInfo<RealType>::* _logIter )  ( ostream & out ) const;
  void ( IterativeInfo<RealType>::* _logEnd )   ( ostream & out ) const;

  virtual void logBeginStandard ( ostream & out ) const = 0;
  virtual void logIterStandard ( ostream & out ) const = 0;
  virtual void logEndStandard  ( ostream & out ) const = 0;

  virtual void logBeginTable   ( ostream & out ) const = 0;
  virtual void logIterTable    ( ostream & out ) const = 0;
  virtual void logEndTable     ( ostream & out ) const = 0;
};

/**** IMPLEMENTATION ********************************************************
 *
 *
 *
 */



template <typename _RealType>
IterativeInfo<_RealType>::
IterativeInfo ( RealType TOL, int maxIter, bool quietMode, ostream & out,
                int numAveragingSteps, int estimatedConvergenceOrder )
    : _TOL ( TOL )
    , _maxIter ( maxIter )
    , _quiet ( quietMode )
    , _megaQuiet ( false )
    , _out ( out )
    , _finalResidual ( NumberTrait<_RealType>::NaN )
    , _currentlyProceedingStep ( false )
    , _convergenceEstimator ( numAveragingSteps, estimatedConvergenceOrder )
    , _logBegin ( &IterativeInfo<_RealType>::logBeginStandard )
    , _logIter  ( &IterativeInfo<_RealType>::logIterStandard )
    , _logEnd   ( &IterativeInfo<_RealType>::logEndStandard ) {
}


template <typename _RealType>
IterativeInfo<_RealType>::
~IterativeInfo () {}

template <typename _RealType>
void IterativeInfo<_RealType>::setAccuracy ( RealType TOL ) {  _TOL = TOL;  }

template <typename _RealType>
void IterativeInfo<_RealType>::setMaxIterations ( int maxIter )
                                                  {   _maxIter = maxIter;   }

template <typename _RealType>
void IterativeInfo<_RealType>::setQuietMode ( bool quiet )
                                                      {   _quiet = quiet;
                                                          if ( quiet == false) 
                                                            this->setMegaQuietMode( quiet );
                                                      }
                                                      
template <typename _RealType>
void IterativeInfo<_RealType>::setMegaQuietMode ( bool quiet )
                                                      {   _quiet = quiet; _megaQuiet = quiet;   }
                                                      
template <typename _RealType>
bool IterativeInfo<_RealType>::getMegaQuietMode (  )
                                                      {  return _megaQuiet;   }

template <typename _RealType>
_RealType IterativeInfo<_RealType>::getAccuracy () const  {  return _TOL;   }

template <typename _RealType>
int IterativeInfo<_RealType>::getMaxIterations () const { return _maxIter;  }

template <typename _RealType>
bool IterativeInfo<_RealType>::getQuietMode () const   {   return _quiet;   }

template <typename _RealType>
int IterativeInfo<_RealType>::getNumAveragingSteps () const
                {   return _convergenceEstimator.getNumAveragingSteps ();   }

template <typename _RealType>
void IterativeInfo<_RealType>::printStats() const {  printStats ( _out );   }

template <typename _RealType>
void IterativeInfo<_RealType>::printStats ( ostream & out ) const
                                         {   ( this->*_logIter ) ( out );   }

template <typename _RealType>
_RealType IterativeInfo<_RealType>::
getLastConvergenceCoeff ( int order ) const
      {   return _convergenceEstimator.getLastConvergenceCoeff ( order );   }

template <typename _RealType>
_RealType IterativeInfo<_RealType>::getFinalResidual() const
                                               {   return _finalResidual;   }

template <typename _RealType>
_RealType IterativeInfo<_RealType>::getLastConvergenceOrder() const
              {   return _convergenceEstimator.getLastConvergenceOrder();   }

template <typename _RealType>
_RealType IterativeInfo<_RealType>::getConvergenceCoeff ( int order ) const
          {   return _convergenceEstimator.getConvergenceCoeff ( order );   }

template <typename _RealType>
_RealType IterativeInfo<_RealType>::getConvergenceOrder() const
                  {   return _convergenceEstimator.getConvergenceOrder();   }

template <typename _RealType>
ostream & IterativeInfo<_RealType>::getOstream() const   {   return _out;   }

template <typename _RealType>
int IterativeInfo<_RealType>::getIterationCount() const
                                                   {   return _stepCount;   }

template <typename _RealType>
bool IterativeInfo<_RealType>::maxIterIsExceeded() const
                     {   return getIterationCount() > getMaxIterations();   }

template <typename _RealType>
bool IterativeInfo<_RealType>::maxIterIsReached() const
                    {   return getIterationCount() >= getMaxIterations();   }

template <typename _RealType>
bool IterativeInfo<_RealType>::currentResidualIsNaN() const {
  return ( ( ! ( getIterationCount() == 0 ) ) && aol::isNaN ( _convergenceEstimator.getLastConvergenceCoeff ( 0 ) ) ); // this cannot be evaluated in the first iteration
}

template <typename _RealType>
void IterativeInfo<_RealType>::getIterationInfoFrom ( const IterativeInfo<_RealType> & info ) {
  _stepCount = info._stepCount;
  _convergenceEstimator = info._convergenceEstimator;
}

template <typename _RealType>
void IterativeInfo<_RealType>::setTablePrintout ( bool tablePrintout ) {
  this->_logBegin = ( tablePrintout
                      ? &IterativeInfo<_RealType>::logBeginTable
                      : &IterativeInfo<_RealType>::logBeginStandard );
  this->_logIter = ( tablePrintout
                     ? &IterativeInfo<_RealType>::logIterTable
                     : &IterativeInfo<_RealType>::logIterStandard );
  this->_logEnd = ( tablePrintout
                    ? &IterativeInfo<_RealType>::logEndTable
                    : &IterativeInfo<_RealType>::logEndStandard );
}



template <typename _RealType>
void IterativeInfo<_RealType>::setNumAveragingSteps ( int n ) {
  int maxOrder = _convergenceEstimator.getNumConvergenceCoeff();
  _convergenceEstimator = ConvergenceEstimator<_RealType> ( n, maxOrder );
}



template <typename _RealType>
void IterativeInfo<_RealType>::startIterations() {
  _stepCount = 0;
  _finalResidual = NumberTrait<_RealType>::NaN;

  _convergenceEstimator.clear();

  if ( !_quiet )
    ( this->*_logBegin ) ( _out );
}



template <typename _RealType>
void IterativeInfo<_RealType>::startStep() {
  _residualWasSet = false;
  _currentlyProceedingStep = true;
  _stepCount++;
}



template <typename _RealType>
void IterativeInfo<_RealType>::finishStep() {
  // only output ...
  if (    _currentlyProceedingStep // if step was not finished yet,
       && _residualWasSet          // any meaningful residual has been set
       && !_quiet                  // we want to output at all
       && !(_stepCount % _convergenceEstimator.getNumAveragingSteps() )
     )                             // we want output in this step
    ( this->*_logIter ) ( _out );

  _finalResidual = getCurrentResidual();
  _currentlyProceedingStep = false;
}



template <typename _RealType>
void IterativeInfo<_RealType>::finishStep ( RealType residual ) {
  setCurrentResidual ( residual );
  finishStep();
}



template <typename _RealType>
void IterativeInfo<_RealType>::finishIterations() {
  if ( !_quiet )
    ( this->*_logEnd ) ( _out );
}



template <typename _RealType>
void IterativeInfo<_RealType>::setCurrentResidual ( RealType residual ) {
  if ( _residualWasSet )
    _convergenceEstimator.setResidualAgain ( residual );
  else
    _convergenceEstimator.push ( residual );
  _residualWasSet = true;
}



template <typename _RealType>
_RealType IterativeInfo<_RealType>::getCurrentResidual() const {
  if ( this->getIterationCount() == 0 )
    return _finalResidual;

  if ( _currentlyProceedingStep && _residualWasSet )
    return _convergenceEstimator.getLastConvergenceCoeff ( 0 );
  else
    throw Exception ( "IterativeInfo::getCurrentResidual(): Trying to read "
                      "current residual that was not yet set.",
                      __FILE__, __LINE__ );
}



} // end of namespace aol.

#endif
