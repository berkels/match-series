#ifndef __NEWTONINFO_H
#define __NEWTONINFO_H

#include <iterativeInfo.h>
#include <solverInfo.h>
#include <pointerClasses.h>
#include <parameterParser.h>

namespace aol {

/****************************************************************************
 *
 *        CLASS NewtonInfo
 */
/**
 *  \brief Collects all administrative infos about a Newton method.
 *
 *  \author von Deylen
 */

template <typename RealType>
class NewtonInfo : public IterativeInfo<RealType> {
public:
  // constructors
  NewtonInfo ( RealType TOL_Newton, int maxIterNewton,
               SolverInfo<RealType> & solverInfo,
               double tauMin = 1E-10,
               double tauMax = 2.,
               int numAveragingSteps = 1 );

  NewtonInfo ( RealType TOL_Newton,
               int maxIterNewton,
               RealType TOL_solver,
               int maxIterSolver,
               StoppingMode mode = STOPPING_UNSET,
               ostream & out = clog,
               double tauMin = 1E-10,
               double tauMax = 2.,
               int numAveragingSteps = 1 );

  virtual ~NewtonInfo () {
  }

  //! deep copy (i. e. SolverInfo will be copied, too)
  explicit NewtonInfo ( const NewtonInfo<RealType> & info );

  static NewtonInfo<RealType> * createFromParameterFile ( string parameterFileName );

  // *** parameter set routines ***

  void setQuietMode ( bool quiet = true );
  //! you can change from normal to simplified Newton scheme
  //! by setting this variable to anything larger that 1.
  void setNumDerivativeHoldingSteps ( int numDerivativeHoldingSteps );

  // *** member get routines ***
  SolverInfo<RealType> & getSolverInfo();
  const SolverInfo<RealType> & getSolverInfo() const;

  RealType getCurrentStepsize() const;
  RealType getCurrentNormUpdate() const;

  RealType getStepsize() const;
  RealType getNormUpdate() const;

  // *** parameter get routines ***
  int getNumDerivativeHoldingSteps () const;

  RealType getTauMin () const {  return _tauMin;  }
  RealType getTauMax () const {  return _tauMax;  }

  void setTauMin ( RealType tauMin ) {  _tauMin = tauMin;  }
  void setTauMax ( RealType tauMax ) {  _tauMax = tauMax;  }

  enum TIMESTEP_CONTROLLER {
    ARMIJO,
    NEWTON_OPTIMAL,
    NLEQ_RES,
    NLEQ_ERR,
    NEWTON_SIMPLE,
    NO_TIMESTEP_CONTROL,
    WOLFE
  };

  TIMESTEP_CONTROLLER getTimestepController() const;

  // set routines for status variables
  void setSolverInfo ( SolverInfo<RealType> & solverInfo );
  void setSolverInfoDeleteFlag ( bool deleteSolverInfo );
  void setTimestepController ( TIMESTEP_CONTROLLER timestepController );
  void setLastStepsize ( RealType stepsize );
  void setLastNormUpdate ( RealType normUpdate );

  void startStep();
  void finishStep ( RealType residual, RealType normUpdate, RealType stepsize );
  using IterativeInfo<RealType>::finishStep;

  //! you may want the iterative solver tolerance to be computed
  //! for each Newton step separately, such that the iterative
  //! solver will not solver "too exact" when Newton's method
  //! is still far from the optimal point.
  //! Suppose the next Newton step could produce k correct decimal
  //! digits in the next step, then this method sets the solver TOL
  //! to 10^k * _solverTOLSafetyFactor (usually 0.01).
  void computeSolverTOL ();

  void setWriteUpdate ( bool writeUpdate );

protected:
  TIMESTEP_CONTROLLER _timestepController;
  RingBufferFrontInsertCapped<RealType> _stepsizeBuffer;
  RingBufferFrontInsertCapped<RealType> _normUpdateBuffer;
  MixedFormat _stepsizeFormat;
  MixedFormat _normUpdateFormat;

  RealType _tauMin;
  RealType _tauMax;
  RealType _solverTOLSafetyFactor;

  int _numDerivativeHoldingSteps;

  bool _stepsizeWasSet;
  bool _writeUpdate;

  virtual void logBeginStandard ( ostream & out ) const;
  virtual void logIterStandard ( ostream & out ) const;
  virtual void logEndStandard  ( ostream & out ) const;

  virtual void logBeginTable   ( ostream & out ) const;
  virtual void logIterTable    ( ostream & out ) const;
  virtual void logEndTable     ( ostream & out ) const;

  DeleteFlagPointer<SolverInfo<RealType> > _solverInfo;
};

/**** IMPLEMENTATION ********************************************************
 *
 *
 *
 */

//---------------------------------------------------------------------------

template <typename RealType>
NewtonInfo<RealType>::
NewtonInfo ( RealType TOL_Newton, int maxIterNewton,
             SolverInfo<RealType> & solverInfo,
             double tauMin,
             double tauMax,
             int numAveragingSteps )
    : IterativeInfo<RealType> ( TOL_Newton, maxIterNewton,
                                false, solverInfo.getOstream(),
                                numAveragingSteps, 2 )
    , _timestepController ( ARMIJO )
    , _stepsizeBuffer ( numAveragingSteps )
    , _normUpdateBuffer ( numAveragingSteps )
    , _stepsizeFormat ( 4, 4 )
    , _normUpdateFormat ( 4, 4 )
    , _tauMin ( tauMin )
    , _tauMax ( tauMax )
    , _solverTOLSafetyFactor ( 0.01 )
    , _numDerivativeHoldingSteps ( 1 )
    , _writeUpdate ( false )
    , _solverInfo ( &solverInfo, false )
{}
//---------------------------------------------------------------------------

template <typename RealType>
NewtonInfo<RealType>::
NewtonInfo ( RealType TOL_Newton, int maxIterNewton,
             RealType TOL_solver, int maxIterSolver,
             StoppingMode mode,
             ostream & out,
             double tauMin,
             double tauMax,
             int numAveragingSteps )
    : IterativeInfo<RealType> ( TOL_Newton, maxIterNewton,
                                false, out, numAveragingSteps, 2 )
    , _timestepController ( ARMIJO )
    , _stepsizeBuffer ( numAveragingSteps )
    , _normUpdateBuffer ( numAveragingSteps )
    , _stepsizeFormat ( 4, 4 )
    , _normUpdateFormat ( 4, 4 )
    , _tauMin ( tauMin )
    , _tauMax ( tauMax )
    , _solverTOLSafetyFactor ( 0.01 )
    , _numDerivativeHoldingSteps ( 1 )
    , _writeUpdate ( false )
    , _solverInfo ( new SolverInfo<RealType>(TOL_solver, maxIterSolver,
                                             mode, false, out), true )
{}
//---------------------------------------------------------------------------

template <typename RealType>
NewtonInfo<RealType>::
NewtonInfo ( const NewtonInfo<RealType> & info )
    : IterativeInfo<RealType> ( info )
    , _timestepController     ( info._timestepController )
    , _stepsizeBuffer         ( info.getNumAveragingSteps() )
    , _normUpdateBuffer       ( info.getNumAveragingSteps() )
    , _stepsizeFormat         ( info._stepsizeFormat )
    , _normUpdateFormat       ( info._stepsizeFormat )
    , _tauMin                 ( info._tauMin )
    , _tauMax                 ( info._tauMax )
    , _solverTOLSafetyFactor  ( info._solverTOLSafetyFactor )
    , _numDerivativeHoldingSteps ( info._numDerivativeHoldingSteps )
    , _solverInfo ( new SolverInfo<RealType> ( info.getSolverInfo() ), true )
{}

//---------------------------------------------------------------------------

template <typename RealType>
NewtonInfo<RealType> * NewtonInfo<RealType>::
createFromParameterFile ( string parameterFileName ) {

  aol::ParameterParser parser ( parameterFileName.c_str() );

  RealType TOL_Newton = parser.getDouble ( "TOL_Newton" );
  int maxIterNewton   = parser.getInt    ( "maxIterNewton" );
  int maxIterSolver   = parser.getInt    ( "maxIterSolver" );

  RealType TOL_Solver    = parser.hasVariable ( "TOL_Solver" ) ?
                           parser.getDouble   ( "TOL_Solver" ) : TOL_Newton * TOL_Newton * 0.01;
  RealType tauMin        = parser.hasVariable ( "tauMin" ) ?
                           parser.getDouble   ( "tauMin" ) : 1E-10;
  RealType tauMax        = parser.hasVariable ( "tauMax" ) ?
                           parser.getDouble   ( "tauMax" ) : 2.;
  int numAveragingSteps  = parser.hasVariable ( "numAveragingSteps" ) ?
                           parser.getInt      ( "numAveragingSteps" ) : 1;
  RealType solverTOLSafetyFactor = parser.hasVariable ( "solverTOLSafetyFactor" ) ?
                           parser.getDouble   ( "solverTOLSafetyFactor" ) : 0.01;
  bool quietMode         = parser.hasVariable ( "quietMode" ) ?
                           (parser.getString   ( "quietMode" ) == "true") : false;
  bool solverQuietMode   = parser.hasVariable ( "solverQuietMode" ) ?
                           (parser.getString   ( "solverQuietMode" ) == "true") : true;
  bool writeUpdate       = parser.hasVariable ( "writeUpdate" ) ?
                           (parser.getString   ( "writeUpdate" ) == "true") : false;
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

  TIMESTEP_CONTROLLER tsc = ARMIJO;
  if ( parser.hasVariable ( "timestepController" ) ) {
    string tscName = parser.getString ( "timestepController" );
    if ( tscName == "NEWTON_OPTIMAL" )      tsc = NEWTON_OPTIMAL;
    if ( tscName == "NLEQ_ERR" )            tsc = NLEQ_ERR;
    if ( tscName == "NLEQ_RES" )            tsc = NLEQ_RES;
    if ( tscName == "WOLFE" )               tsc = WOLFE;
    if ( tscName == "NEWTON_SIMPLE" )       tsc = NEWTON_SIMPLE;
    if ( tscName == "NO_TIMESTEP_CONTROL" ) tsc = NO_TIMESTEP_CONTROL;
  }

  NewtonInfo<RealType> * newInfo = new NewtonInfo<RealType>
              ( TOL_Newton, maxIterNewton, TOL_Solver, maxIterSolver,
                stoppingMode, *out, tauMin, tauMax, numAveragingSteps );

  newInfo->setQuietMode ( quietMode );
  newInfo->setWriteUpdate ( writeUpdate );
  newInfo->getSolverInfo().setQuietMode ( solverQuietMode );
  newInfo->setTimestepController ( tsc );
  newInfo->_solverTOLSafetyFactor = solverTOLSafetyFactor;

  if ( parser.hasVariable ( "logMode" )
       && parser.getString ( "logMode" ) == "table" )
    newInfo->setTablePrintout ();

  if ( parser.hasVariable ( "numDerivativeHoldingSteps" ) )
    newInfo->_numDerivativeHoldingSteps = parser.getInt ( "numDerivativeHoldingSteps" );

  return newInfo;
}
//---------------------------------------------------------------------------

template <typename RealType>
void NewtonInfo<RealType>::setQuietMode ( bool quiet ) {
  IterativeInfo<RealType>::setQuietMode ( quiet );
  getSolverInfo().setQuietMode ( quiet );
}
//---------------------------------------------------------------------------
template <typename RealType>
void NewtonInfo<RealType>::
setNumDerivativeHoldingSteps ( int numDerivativeHoldingSteps )
              {   _numDerivativeHoldingSteps = numDerivativeHoldingSteps;   }
//---------------------------------------------------------------------------
template <typename RealType>
SolverInfo<RealType> & NewtonInfo<RealType>::
getSolverInfo()                                  {   return *_solverInfo;   }
//---------------------------------------------------------------------------
template <typename RealType>
const SolverInfo<RealType> & NewtonInfo<RealType>::
getSolverInfo() const                            {   return *_solverInfo;   }
//---------------------------------------------------------------------------

template <typename RealType>
void NewtonInfo<RealType>::startStep() {
  IterativeInfo<RealType>::startStep();
  _stepsizeWasSet = false;
}
//---------------------------------------------------------------------------

template <typename RealType>
RealType NewtonInfo<RealType>::
getCurrentStepsize() const {
  if ( _stepsizeWasSet )
    return _stepsizeBuffer[0];
  else
    throw Exception ( "NewtonInfo::getCurrentStepsize(): Trying to read "
                      "current stepsize that was not yet set.",
                      __FILE__, __LINE__ );
}

//---------------------------------------------------------------------------

template <typename RealType>
RealType NewtonInfo<RealType>::
getCurrentNormUpdate() const {
  if ( _stepsizeWasSet )
    return _normUpdateBuffer[0];
  else
    throw Exception ( "NewtonInfo::getCurrentNormUpdate(): Trying to read "
                      "current updaten norm that was not yet set.",
                      __FILE__, __LINE__ );
}

//---------------------------------------------------------------------------
template <typename RealType>
RealType NewtonInfo<RealType>::
getStepsize() const          {   return _stepsizeBuffer.arithmeticMean();   }
//---------------------------------------------------------------------------
template <typename RealType>
RealType NewtonInfo<RealType>::
getNormUpdate() const      {   return _normUpdateBuffer.arithmeticMean();   }
//---------------------------------------------------------------------------
template <typename RealType>
int NewtonInfo<RealType>::getNumDerivativeHoldingSteps () const
                                   {   return _numDerivativeHoldingSteps;   }
//---------------------------------------------------------------------------
template <typename RealType>
typename NewtonInfo<RealType>::TIMESTEP_CONTROLLER NewtonInfo<RealType>::
getTimestepController() const             {   return _timestepController;   }
//---------------------------------------------------------------------------
template <typename RealType>
void NewtonInfo<RealType>::
setSolverInfo ( SolverInfo<RealType> & solverInfo )
                            {   _solverInfo.reset ( &solverInfo, false );   }
//---------------------------------------------------------------------------
template <typename RealType>
void NewtonInfo<RealType>::setSolverInfoDeleteFlag ( bool deleteSolverInfo )
                      {   _solverInfo.setDeleteFlag ( deleteSolverInfo );   }
//---------------------------------------------------------------------------
template <typename RealType>
void NewtonInfo<RealType>::
setTimestepController ( TIMESTEP_CONTROLLER timestepController )
                            {   _timestepController = timestepController;   }
//---------------------------------------------------------------------------

template <typename RealType>
void NewtonInfo<RealType>::
setLastStepsize ( RealType stepsize ) {
  _stepsizeBuffer.push ( stepsize );
  _stepsizeWasSet = true;
}
//---------------------------------------------------------------------------

template <typename RealType>
void NewtonInfo<RealType>::
setLastNormUpdate ( RealType normUpdate ) {
  _normUpdateBuffer.push ( normUpdate );
}
//---------------------------------------------------------------------------

template <typename RealType>
void NewtonInfo<RealType>::
finishStep ( RealType residual, RealType normUpdate, RealType stepsize ) {
  setLastStepsize ( stepsize );
  setLastNormUpdate ( normUpdate );
  finishStep ( residual );
}
//---------------------------------------------------------------------------

template <typename RealType>
void NewtonInfo<RealType>::computeSolverTOL () {

  RealType lastRes = this->_convergenceEstimator.getLastConvergenceCoeff ( 0 );
  RealType expRes = lastRes > 1. ? sqrt ( lastRes ) : Sqr ( lastRes );

  if ( expRes < Sqr ( this->getAccuracy() ) )
    expRes = Sqr ( this->getAccuracy () );

  if ( getSolverInfo().getAccuracy() > Sqr ( expRes * _solverTOLSafetyFactor ) )
    getSolverInfo().setAccuracy ( Sqr ( expRes * _solverTOLSafetyFactor ) );
}
//---------------------------------------------------------------------------
template <typename RealType>
void NewtonInfo<RealType>::setWriteUpdate ( bool writeUpdate )
                                            {  _writeUpdate = writeUpdate;  }
//---------------------------------------------------------------------------
template <typename RealType>
void NewtonInfo<RealType>::logBeginStandard ( ostream & ) const            {}
//---------------------------------------------------------------------------
template <typename RealType>
void NewtonInfo<RealType>::logEndStandard ( ostream & out ) const {
  if ( ! this->getQuietMode() ) {
    out << endl;
  }
}

//---------------------------------------------------------------------------
template <typename RealType>
void NewtonInfo<RealType>::
logIterStandard ( ostream & out ) const {
  if ( ! this->getQuietMode() ) {
    out << intFormat ( this->getIterationCount() )
        << " iterations, stepsize "
        << _stepsizeFormat ( getCurrentStepsize() );
    if ( _writeUpdate ) out
        <<", norm update "
        << _normUpdateFormat ( getCurrentNormUpdate() );
    out << ", residuum: "
        << scientificFormat ( this->getCurrentResidual() );
#ifndef VERBOSE
    if( !getSolverInfo().getQuietMode() )
      out << "\n";
    else
      out << "\r";
#else
    out << endl;
#endif
  }
}

//---------------------------------------------------------------------------

template <typename RealType>
void NewtonInfo<RealType>::logBeginTable ( ostream & out ) const {
  if ( ! this->getQuietMode() ) {
    out << "iter \t"
           "av. stepsize \t";
    if ( _writeUpdate ) out <<
           "av. update\t";
    out << "residual \t"
           "solver steps \t"
           "solver residual"
        << endl;
  }
}
//---------------------------------------------------------------------------

template <typename RealType>
void NewtonInfo<RealType>::logEndTable ( ostream & out ) const {
  if ( ! this->getQuietMode() )
    out << endl;
}

//---------------------------------------------------------------------------

template <typename RealType>
void NewtonInfo<RealType>::
logIterTable ( ostream & out ) const {
  if ( ! this->getQuietMode() ) {
    out << this->getIterationCount() << "\t"
        << _stepsizeFormat ( getStepsize() ) << "\t";
    if ( _writeUpdate ) out
        << _normUpdateFormat ( getNormUpdate() ) << "\t";
    out << scientificFormat ( this->getLastConvergenceCoeff ( 0 ) ) << "\t"
        << intFormat ( getSolverInfo().getIterationCount() ) << "\t"
        << scientificFormat ( sqrt ( getSolverInfo().getFinalResidual() ) ) << "\r";
    if ( ! ( this->getIterationCount()
                     % this->_convergenceEstimator.getNumAveragingSteps() ) )
      out << endl;
  }
}

//---------------------------------------------------------------------------

} // end of namespace aol.

#endif
