#ifndef __SOLVERINFO_H
#define __SOLVERINFO_H

#include <iterativeInfo.h>
#include <parameterParser.h>

namespace aol {

//! Type of stopping criterion used in IterativeInverseOps,
//! residual is usually measured in (l_2-norm)^2
enum StoppingMode {
  STOPPING_UNSET,                         //!< Behaves as STOPPING_ABSOLUTE but gives a warning
  STOPPING_ABSOLUTE,                      //!< Stop if \f$ |r_n|^2 < TOL \f$
  STOPPING_RELATIVE_TO_INITIAL_RESIDUUM,  //!< Stop if \f$ |r_n|^2 / |r_0|^2 < TOL \f$
  STOPPING_RELATIVE_TO_RIGHT_HAND_SIDE    //!< Stop if \f$ |r_n|^2 / |b|^2 < TOL \f$
};

/****************************************************************************
 *
 *        CLASS SolverInfo
 */
/**
 *  \brief Collects all administrative infos about an iterative
 *         linear solver.
 *
 *  \author von Deylen
 */
template <typename RealType>
class SolverInfo : public IterativeInfo<RealType> {
public:
  SolverInfo ( RealType TOL,
               int maxIter,
               StoppingMode mode = STOPPING_UNSET,
               bool quietMode = false,
               ostream & _out = clog,
               int numAveragingSteps = 1 );

  virtual ~SolverInfo () {}

  static SolverInfo<RealType> * createFromParameterFile ( string parameterFileName );

  StoppingMode getStoppingMode() const;
  void setStoppingMode ( StoppingMode mode );

  bool stoppingCriterionIsFulfilled() const;

  using IterativeInfo<RealType>::startIterations;
  using IterativeInfo<RealType>::setMegaQuietMode;
  using IterativeInfo<RealType>::getMegaQuietMode;
  void startIterations ( RealType normRhs,
                         RealType normInitialResiduum,
                         string methodName,
                         string residualDescr );

  void finishIterations();
  void finishIterations ( const RealType residual );

  void setNormRhs ( RealType normRhs );
  void setNormInitialResiduum ( RealType normInRes );

protected:
  StoppingMode _stoppingMode;
  RealType     _normRhs;
  RealType     _normInitialResiduum;

  RealType getTolerableError() const;

  virtual void logBeginStandard ( ostream & out ) const;
  virtual void logIterStandard ( ostream & out ) const;
  virtual void logEndStandard  ( ostream & out ) const;

  virtual void logBeginTable   ( ostream & out ) const;
  virtual void logIterTable    ( ostream & out ) const;
  virtual void logEndTable     ( ostream & out ) const;

private:
  string _methodName;
  string _residualDescr;
};

/**** IMPLEMENTATION ********************************************************
 *
 *
 *
 */

template <typename RealType>
SolverInfo<RealType>::
SolverInfo ( RealType TOL,
             int maxIter,
             StoppingMode mode,
             bool quietMode,
             ostream & out,
             int numAveragingSteps )
    : IterativeInfo<RealType> ( TOL, maxIter, quietMode, out, numAveragingSteps, 1 )
    , _stoppingMode ( mode )
    , _normRhs ( -1. )
    , _normInitialResiduum ( -1. ) {}


template <typename RealType>
SolverInfo<RealType> *
SolverInfo<RealType>::createFromParameterFile ( string parameterFileName ) {

  aol::ParameterParser parser ( parameterFileName.c_str() );

  RealType TOL           = parser.getDouble ( "TOL" );
  int maxIter            = parser.getInt ( "maxIter" );
  int numAveragingSteps  = parser.hasVariable ( "numAveragingSteps" ) ?
                           parser.getInt      ( "numAveragingSteps" ) : 1;
  bool quietMode         = parser.hasVariable ( "quietMode" ) ?
                           (parser.getString  ( "quietMode" ) == "true") : false;

  ostream * out = &clog;
  if ( parser.hasVariable ( "outStream" ) ) {
    string streamName = parser.getString ( "outStream" );
    if ( streamName == "cout" ) out = &cout;
    if ( streamName == "cerr" ) out = &cerr;
    if ( streamName == "clog" ) out = &clog;
  }

  StoppingMode stoppingMode = STOPPING_UNSET;
  if ( parser.hasVariable ( "stoppingMode" ) ) {
    string stopName = parser.getString ( "stoppingMode" );
    if ( stopName == "STOPPING_ABSOLUTE" )
      stoppingMode = STOPPING_ABSOLUTE;
    if ( stopName == "STOPPING_RELATIVE_TO_INITIAL_RESIDUUM" )
      stoppingMode = STOPPING_RELATIVE_TO_INITIAL_RESIDUUM;
    if ( stopName == "STOPPING_RELATIVE_TO_RIGHT_HAND_SIDE" )
      stoppingMode = STOPPING_RELATIVE_TO_RIGHT_HAND_SIDE;
  }

  SolverInfo<RealType> * newInfo
    = new SolverInfo<RealType> ( TOL, maxIter, stoppingMode,
                                 quietMode, *out, numAveragingSteps );

  if ( parser.hasVariable ( "logMode" ) )
    if ( parser.getString ( "logMode" ) == "table" )
      newInfo->setTablePrintout ();

  return newInfo;
}


template <typename RealType>
StoppingMode SolverInfo<RealType>::getStoppingMode() const
                                                {   return _stoppingMode;   }

template <typename RealType>
void SolverInfo<RealType>::setStoppingMode ( StoppingMode mode )
                                                {   _stoppingMode = mode;   }

template <typename RealType>
void SolverInfo<RealType>::setNormRhs ( RealType normRhs )
                                                  {   _normRhs = normRhs;   }

template <typename RealType>
void SolverInfo<RealType>::setNormInitialResiduum ( RealType normInRes )
                                    {   _normInitialResiduum = normInRes;   }

template <typename RealType>
RealType SolverInfo<RealType>::getTolerableError() const {
  RealType error = 0.;
  switch ( this->getStoppingMode() ) {
  case STOPPING_UNSET:
  case STOPPING_ABSOLUTE:
    error = this->getAccuracy();
    break;

  case STOPPING_RELATIVE_TO_INITIAL_RESIDUUM:
    error = this->getAccuracy() * ( _normInitialResiduum != aol::NumberTrait<RealType>::zero ? _normInitialResiduum : aol::NumberTrait<RealType>::one );
    break;

  case STOPPING_RELATIVE_TO_RIGHT_HAND_SIDE:
    error = this->getAccuracy() * ( _normInitialResiduum != aol::NumberTrait<RealType>::zero ? _normInitialResiduum : aol::NumberTrait<RealType>::one );
    break;

  default:
    throw aol::UnimplementedCodeException ( "SolverInfo::getTolerableError: unknown StoppingMode", __FILE__, __LINE__ );
  }
  return error;
}


template <typename RealType>
bool SolverInfo<RealType>::stoppingCriterionIsFulfilled() const {

  if ( this->_currentlyProceedingStep)
    return ( this->getCurrentResidual() < getTolerableError() );
  else
    return ( this->getFinalResidual() < getTolerableError() );
}

template <typename RealType>
void SolverInfo<RealType>::startIterations ( RealType normRhs,
                                             RealType normInitialResidual,
                                             string methodName,
                                             string residualDescr ) {
  setNormRhs ( normRhs );
  setNormInitialResiduum ( normInitialResidual );

  this->_residualWasSet = false;
  this->setCurrentResidual ( normInitialResidual );

  _methodName = methodName;
  _residualDescr = residualDescr;

  // stopping criterion set?
  if ( this->getStoppingMode() == STOPPING_UNSET )
    cerr << aol::color::red << "Stopping criterion not set. "
    "Will use absolute stopping." << aol::color::reset << endl;

  IterativeInfo<RealType>::startIterations();
  this->_finalResidual = normInitialResidual;
}


template <typename RealType>
void SolverInfo<RealType>::logBeginStandard ( ostream & out ) const {

  if ( !this->getQuietMode() )
    out << "Starting " << _methodName << " method, "
           "measuring " << _residualDescr << "," << endl
        << "initial residual: "
        << aol::scientificFormat ( _normInitialResiduum ) << endl;

}


template <typename RealType>
void SolverInfo<RealType>::logIterStandard ( ostream & out ) const {
  if ( !this->getQuietMode() ) {
    out << intFormat ( this->getIterationCount() )
        << " iterations, residual: "
        << scientificFormat ( this->getCurrentResidual() ) << "\r";
#ifdef VERBOSE
    out << endl;
#endif
  }
}


template <typename RealType>
void SolverInfo<RealType>::logEndStandard ( ostream & out ) const {
  if ( !this->getQuietMode() )
    // uncomment this line (and delete the following)
    // to erase inbetween output:
    // out << "                                                          \r";
    out << endl;
}


template <typename RealType>
void SolverInfo<RealType>::logBeginTable ( ostream & out ) const {
  if ( this->getQuietMode() )
    return;

  out << "Starting " << _methodName << " method, "
         "measuring " << _residualDescr << "," << endl
      << "initial residual: "
      << aol::scientificFormat ( _normInitialResiduum ) << endl
      << endl
      << "iter \t"
         " residual \t\t"
         "av. rate" << endl;
}


template <typename RealType>
void SolverInfo<RealType>::logIterTable ( ostream & out ) const {
  if ( this->getQuietMode() )
    return;

  out << intFormat ( this->getIterationCount() ) << "\t"
      << scientificFormat ( this->getLastConvergenceCoeff ( 0 ) ) << "\t "
      << shortFormat ( this->getConvergenceCoeff ( 1 ) )
      << "\r";
  if ( ! ( this->getIterationCount() % this->_convergenceEstimator.getNumAveragingSteps() ) )
    out << endl;
}


template <typename RealType>
void SolverInfo<RealType>::logEndTable ( ostream & out ) const {
  if ( this->getQuietMode() )
    return;

  // if necessary, print last row:
  if ( this->getIterationCount() % this->_convergenceEstimator.getNumAveragingSteps() )
  {
    out << intFormat ( this->getIterationCount() ) << "\t"
        << scientificFormat ( this->getLastConvergenceCoeff ( 0 ) ) << "\t "
        << shortFormat ( this->getConvergenceCoeff ( 1 ) )
        << endl;

    // blank line after table:
    out << endl;
  }
}




template <typename RealType>
void SolverInfo<RealType>::finishIterations ( const RealType exactResidual ) {
  const RealType estimatedResidual = this->_finalResidual;
  this->_finalResidual = exactResidual;
  finishIterations();

  if ( !this->getQuietMode() )
    this->_out << "recomputed residual:       "
               << scientificFormat ( exactResidual ) << endl;

  if ( !this->getMegaQuietMode() && ( ( fabs ( estimatedResidual / exactResidual ) > 10.0 ) || ( fabs ( estimatedResidual / exactResidual ) < 0.1 ) ) ) {
    this->_out << aol::color::red << "Estimated residual " << aol::detailedFormat ( estimatedResidual )  << " and recomputed residual " << aol::detailedFormat ( exactResidual ) << " differ by more than one order of magnitude." << aol::color::reset << endl;
  }
  
}


template <typename RealType>
void SolverInfo<RealType>::finishIterations () {
  IterativeInfo<RealType>::finishIterations();

  // It would seem natural to call maxIterIsReached() instead
  // of testing the max iteration number by hand. But in
  // InterruptableIterativeInfo, maxIterIsReached() is
  // overloaded with a slightly different behaviour.
  if ( !this->getMegaQuietMode() ) {
    if ( ( this->getIterationCount() >= this->getMaxIterations() ) && !stoppingCriterionIsFulfilled() ) {
      this->_out << aol::color::red << "Too many iterations: solver has not converged to " << getTolerableError() << " within " << this->getIterationCount() << " iterations." << aol::color::reset << endl;
    }
    if ( this->currentResidualIsNaN() ) {
      this->_out << aol::color::red << "NaN residual found after " << this->getIterationCount() << " iterations." << aol::color::reset << endl;
    }
  }
}

} // end of namespace aol.

#endif
