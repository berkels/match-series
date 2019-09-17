#ifndef __GRADIENTDESCENT_H
#define __GRADIENTDESCENT_H

#include <ArmijoSearch.h>
#include <pointerClasses.h>
#include <timestepSaver.h>
#ifdef USE_CPP11
#include <chrono>
#endif

namespace aol {

/**
 * \brief General-purpose energy minimization using an explicit, step size controlled gradient descent.
 *
 * \author Berkels
 * \ingroup Optimization
 */
template <typename RealType, typename VectorType, typename DerivativeType = VectorType>
class GradientDescentBase
: public ArmijoLineSearchUsingOp<RealType, VectorType, VectorType, DerivativeType > {
public:
  enum TIMESTEP_CONTROLLER {
    ARMIJO,
    GRADIENT_SIMPLE,
    NO_TIMESTEP_CONTROL,
    POWELL_WOLFE
  };
  //! With these flags the behavior of this class can be configured.
  enum GRADIENT_DESCENT_FLAGS {
    //! Use a nonlinear conjugate gradient method instead of a simple gradient descent.
    USE_NONLINEAR_CG = 1,
    //! Do not call postProcess in apply before starting the gradient descent.
    DO_NOT_POSTPROCESS_STARTING_VALUE = 2,
    //! Do not apply smoothDirection on the gradient when calculation the descent direction.
    //! This is like a gradient flow with the Euclidean metric.
    DO_NOT_SMOOTH_DESCENT_DIRECTION = 4,
    //! Do not output any of the standard info messages with cerr.
    DO_NOT_WRITE_CONSOLE_OUTPUT = 8,
    //! Write the info for each iteration in a new line instead of overwriting the previous info.
    USE_MULTILINE_CONSOLE_OUTPUT = 16,
    //! Use the norm of the gradient instead of the energy change as stopping criterion.
    USE_GRADIENT_BASED_STOPPING = 32,
    //! In case the energy and/or the derivative Op depend on the argument of applySingle,
    //! i.e. the functionals themself change when Dest is altered, this flag has to be set.
    //! Otherwise the stopping criterion will not work properly.
    CALC_STOPPING_ENERGY_BEFORE_UPDATING_CURRENT_POSITION = 64,
    //! Include the norm of the gradient of the position at the beginning of the step in the
    //! standard info that is written to the console.
    LOG_GRADIENT_NORM_AT_OLD_POSITION = 128,
    LOG_SECONDS_AND_EVALS_PER_STEP = 256
  };
protected:
  const aol::Op<VectorType, aol::Scalar<RealType> > &_E;
  const aol::Op<VectorType, DerivativeType> &_DE;
  mutable DerivativeType *_pVecDE;
  // Cached value needed by ArmijoLineSearchHelpFunction_evaluate.
  mutable RealType _energyAtPosition;
  // Values saved by the last two calls of ArmijoLineSearchHelpFunction_evaluate.
  mutable aol::RingBuffer<RealType> _lastEnergyTimestepWidth;
  mutable aol::RingBuffer<RealType> _energyForLastTimestepWidth;
  RealType _stopEpsilon;
  mutable int _iterations;
  int _maxIterations;
  bool _updateFilterWidth;
  mutable RealType _filterWidth;
  RealType _stopFilterWidth;
  RealType _updateFilterWidthTolerance;
  aol::SimpleFormat _iterationsFormat;
  aol::MixedFormat _tauFormat;
  std::ostream *_out;
  aol::MixedFormat _sigmaFormat;
  aol::MixedFormat _energyFormat;
  std::ostream *_outFilterWidth;
  int _maximumNumIterationsPerFilterWidth;
  TIMESTEP_CONTROLLER _timestepController;
  unsigned int _configurationFlags;
  //! postProcess is called every _postProcessingOffset steps (default 1, i.e. post processing after every step).
  unsigned int _postProcessingOffset;
  mutable double _timeOneStep;
  mutable int _countEvals;
  const StepSaverBase<RealType, VectorType> *_pStepSaver;
public:
  GradientDescentBase( const aol::Op<VectorType, aol::Scalar<RealType> > &E,
                       const aol::Op<VectorType, DerivativeType> &DE,
                       const int MaxIterations = 50,
                       const RealType StartTau = aol::ZOTrait<RealType>::one,
                       const RealType StopEpsilon = aol::ZOTrait<RealType>::zero )
  : ArmijoLineSearchUsingOp<RealType, VectorType, VectorType, DerivativeType >
                      (0.5, MaxIterations, StartTau),
    _E(E),
    _DE(DE),
    _pVecDE(NULL),
    _lastEnergyTimestepWidth(2),
    _energyForLastTimestepWidth(2),
    _stopEpsilon(StopEpsilon),
    _iterations(0),
    _maxIterations(MaxIterations),
    _updateFilterWidth(false),
    _filterWidth(1.),
    _stopFilterWidth(0.0001),
    _updateFilterWidthTolerance(0.00001),
    _iterationsFormat(aol::countDigitsOfNumber(_maxIterations), 0, ios::fixed | ios::right),
    _tauFormat ( 4, 4 ),
    _out ( NULL ),
    _sigmaFormat ( 1, 4 ),
    _energyFormat( 1, 12 ),
    _outFilterWidth(NULL),
    _maximumNumIterationsPerFilterWidth(0),
    _timestepController(ARMIJO),
    _configurationFlags ( 0 ),
    _postProcessingOffset ( 1 ),
    _countEvals ( 0 ),
    _pStepSaver ( NULL )
  {
  }
  virtual ~GradientDescentBase()
  {
  }
protected:
  virtual RealType getStepSizeWithSelectedController( DerivativeType &DescentDir,
                                              const VectorType &CurrentPosition,
                                              const RealType OldTau ) const {
    switch ( _timestepController ) {
    case ARMIJO:
      return this->getTimestepWidthWithArmijoLineSearch(DescentDir, CurrentPosition, OldTau);
      break;
    case GRADIENT_SIMPLE:
      return this->getTimestepWidthWithSimpleLineSearch(DescentDir, CurrentPosition);
      break;
    case NO_TIMESTEP_CONTROL:
      return OldTau;
      break;
    case POWELL_WOLFE:
      return this->getTimestepWidthWithPowellWolfeLineSearch( DescentDir, CurrentPosition, OldTau );
      break;
    default:
      throw Exception( "Unhandled timestep control mode", __FILE__, __LINE__);
    }
    return aol::NumberTrait<RealType>::zero;
  }

  virtual void getTauAndUpdateDescentDir( DerivativeType &DescentDir,
                                          const VectorType &CurrentPosition,
                                          aol::Vector<RealType> &Tau ) const {
    Tau[0] = getStepSizeWithSelectedController(DescentDir, CurrentPosition, Tau[0]);
    DescentDir *= Tau[0];
  }
  virtual bool postProcess( VectorType &/*Position*/, const int /*Iteration*/ ) const {
    return false;
  }
  virtual void smoothDirection( const VectorType &CurrentPosition, const RealType Sigma, DerivativeType &Direction ) const = 0;
  // Attention: ArmijoLineSearchHelpFunction_evaluate only works, if _energyAtPosition contains E(CurrentPosition)!
  RealType ArmijoLineSearchHelpFunction_evaluate( const DerivativeType &DescentDir,
                                                  const VectorType &CurrentPosition,
                                                  const RealType timestepWidth ) const
  {
    // In this case we can use the cached value.
    if ( timestepWidth == 0 )
      return _energyAtPosition;
    else {
      VectorType newPosition ( CurrentPosition );
      aol::Scalar<RealType> Energy;
      newPosition.addMultiple ( DescentDir, timestepWidth );
      _E.apply( newPosition, Energy );

      // Save these values, they may save us from evaluating the energy at the new position again.
      _lastEnergyTimestepWidth.push_back ( timestepWidth );
      _energyForLastTimestepWidth.push_back ( Energy[0] );

      _countEvals++;
      return ( Energy[0] );
    }
  }
  // Attention: ArmijoLineSearchHelpFunction_evaluateDerivative only works, if _vecDE is filled with DE(position)!
  RealType ArmijoLineSearchHelpFunction_evaluateDerivative( const DerivativeType &DescentDir,
                                                            const VectorType &/*Position*/) const{
    //aol::MultiVector<RealType> tmp(ConfiguratorType::Dim, _grid.getSize());
    //_DE.apply(Position, tmp);
    //return (tmp*DescentDir);
    return ((*_pVecDE).dotProduct(DescentDir));
  }

  virtual RealType ArmijoLineSearchHelpFunction_evaluateDerivativeWithTau( const DerivativeType &DescentDir,
                                                                           const VectorType &CurrentPosition,
                                                                           const RealType TimestepWidth ) const {
    // In this case, we can use the faster evaluation function.
    if ( TimestepWidth == 0 )
      return ArmijoLineSearchHelpFunction_evaluateDerivative ( DescentDir, CurrentPosition );
    else {
      VectorType newPosition ( CurrentPosition );
      VectorType tmp ( CurrentPosition, aol::STRUCT_COPY );
      newPosition.addMultiple ( DescentDir, TimestepWidth );

      _DE.apply( newPosition, tmp );
      return (tmp*DescentDir);
    }
  }

  virtual void calcTaudAndDescentDir( VectorType &CurrentPosition, aol::Vector<RealType> &Tau, DerivativeType &DescentDir, DerivativeType &OldGradient, DerivativeType &OldDirection, const int Iterations ) const{
#ifdef VERBOSE
    cerr << "Applying DE ";
#endif
    _DE.apply(CurrentPosition, DescentDir);                           // evaluate E'[x^k]
#ifdef VERBOSE
    cerr << "done\n";
#endif
    *_pVecDE = DescentDir;

    DescentDir *= -1.;

#ifdef VERBOSE
    cerr << "Smoothing descent direction ";
#endif
    if ( !( _configurationFlags & DO_NOT_SMOOTH_DESCENT_DIRECTION ) )
      smoothDirection( CurrentPosition, _filterWidth, DescentDir );                   // apply A^{-1} to -E'[x^k]
#ifdef VERBOSE
    cerr << "done\n";
#endif
    // Fletcher-Reeves nonlinear conjugate gradient method.
    if ( _configurationFlags & USE_NONLINEAR_CG ) {
      const RealType nonlinCGBeta = ( ( Iterations - 1 ) % 10 ) ? ( DescentDir.normSqr() / OldGradient.normSqr() ) : aol::NumberTrait<RealType>::zero;
      OldGradient = DescentDir;
      DescentDir.addMultiple ( OldDirection, nonlinCGBeta );
      OldDirection = DescentDir;
    }

    getTauAndUpdateDescentDir( DescentDir, CurrentPosition, Tau );

    // If the nonlinear CG direction is not a descent direction, use the normal descent direction instead.
    if ( ( _configurationFlags & USE_NONLINEAR_CG ) && ( Tau.sum() == 0 ) ) {
      DescentDir = OldGradient;
      OldDirection.setZero();
      getTauAndUpdateDescentDir( DescentDir, CurrentPosition, Tau );
    }
  }
protected:
  bool writeEnergyOfGenEnergyOp ( const aol::Op<VectorType, aol::Scalar<RealType> > &Op, const RealType ScaleFactor = aol::ZOTrait<RealType>::one ) const{
    if ( _out == NULL )
      return false;

    aol::Vector<RealType> energyComponentVector;
    const aol::GenEnergyOp<VectorType, RealType>* pGenEnergyOp = dynamic_cast<const aol::GenEnergyOp<VectorType, RealType>*>( &Op );
    if( pGenEnergyOp != NULL ) {
        pGenEnergyOp->getLastEnergy( energyComponentVector );
        for( int i = 0; i < energyComponentVector.size(); i++ )
          *(_out) << ( ScaleFactor * energyComponentVector[i] ) << " ";
        return true;
      }
      else
        return false;
  }
public:
  virtual void writeEnergy( const RealType CurrentEnergy ) const{
    if( _out != NULL )
      *(_out) << (this->_iterations)-1 << " ";
    if( _out != NULL ){
      if ( writeEnergyOfGenEnergyOp ( _E ) == false ) {
        const aol::LinCombOp<VectorType, aol::Scalar<RealType> > *pLinCombOp = dynamic_cast<const aol::LinCombOp<VectorType, aol::Scalar<RealType> >*>(&_E);
        if ( pLinCombOp ) {
          const list<const aol::Op<VectorType, aol::Scalar<RealType> > *> &opList = pLinCombOp->getOpList();
          typename list<const aol::Op<VectorType, aol::Scalar<RealType> >* >::const_iterator opIt;
          const list<RealType> &coeffList = pLinCombOp->getCoeffList();
          typename list<RealType>::const_iterator coeffIt;
          for ( opIt = opList.begin(), coeffIt = coeffList.begin(); opIt != opList.end(); ++opIt, ++coeffIt )
            writeEnergyOfGenEnergyOp ( *(*opIt), *coeffIt );
        }
        *(_out) << CurrentEnergy;
      }
      *(_out) << endl;
    }
  }
  void writeConsoleOutput( const RealType CurrentEnergy, const aol::Vector<RealType> &Tau ) const {
    cerr << this->_iterationsFormat(this->_iterations) << " steps, tau: ";
    for( int i = 0; i < Tau.size(); i++)
      cerr << this->_tauFormat(Tau[i]);
    if ( _updateFilterWidth == true )
      cerr << ", sigma: " << _sigmaFormat(_filterWidth);
    cerr << " energy: " << _energyFormat ( CurrentEnergy );
    if ( _configurationFlags & LOG_GRADIENT_NORM_AT_OLD_POSITION )
      cerr << " DEnorm: " << _pVecDE->norm();
#ifdef USE_CPP11
    if ( _configurationFlags & LOG_SECONDS_AND_EVALS_PER_STEP )
      cerr << " seconds: " << _timeOneStep << " evals: " << _countEvals << "     ";
#endif
    if ( _configurationFlags & USE_MULTILINE_CONSOLE_OUTPUT )
      cerr << "\n";
    else
      cerr << "\r";
    _countEvals = 0;
  }

public:
  // x^{k+1} = x^k - tau * A^{-1} E'[x^k]
  virtual void applySingle( VectorType &Dest ) const{
    _pVecDE = new DerivativeType(Dest, aol::STRUCT_COPY);                      // = E'[x^k]
    aol::Vector<RealType> tau(1);
    tau[0] = this->_startTau;
    aol::Scalar<RealType> energyScalar;
    // Since we post process the current solution after each iteration,
    // it shouldn't hurt to do this on the initial data. This way each
    // iteration starts with post processed input.
    // Nevertheless, give the user the ability to deactivate it.
    if ( !( _configurationFlags & DO_NOT_POSTPROCESS_STARTING_VALUE ) )
      postProcess( Dest, 0 );
    _E.apply(Dest, energyScalar);                             // just for the stopping criterion
    RealType energy = energyScalar[0];
    RealType energyNew = energyScalar[0];

    if ( !( _configurationFlags & DO_NOT_WRITE_CONSOLE_OUTPUT ) )
      cerr << "Initial energy: " << _energyFormat ( energy ) << endl;

    _iterations = 0;
    int iterationsPerFilterWidth = 0;

    // Also log the initial energy.
    writeEnergy( energyNew );
    
    // Save the current time step if _pStepSaver was initialized with setStepSaverReference
    if ( _pStepSaver )
      _pStepSaver->saveStep ( Dest, _iterations );

    // This is the old method for saving time steps, which is kept here for compatibility reasons
    if (this->checkSaveConditions(_iterations))
      this->writeTimeStep( Dest, _iterations );     // virtual function to call the TimestepSaver

    bool stoppingCriterion = false;

    DerivativeType descentDir ( Dest, aol::STRUCT_COPY );
    aol::DeleteFlagPointer<DerivativeType> oldDirection;
    aol::DeleteFlagPointer<DerivativeType> oldGradient;
    // Only allocate memory for the nonlin CG storage variables if necessary.
    if ( _configurationFlags & USE_NONLINEAR_CG ) {
      oldDirection.reset ( new DerivativeType ( Dest, aol::STRUCT_COPY ), true );
      oldGradient.reset ( new DerivativeType ( Dest, aol::STRUCT_COPY ), true );
    }

    if ( _maxIterations <= 0 ) {
      cerr << "No descent steps were allowed since ( _maxIterations <= 0 ) is true.\n";
      return;
    }

    do{                                                       // iteration loop
        
#ifdef USE_CPP11
      std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
#endif
      _iterations++;
      // Save the energy at the current position.
      energy = energyNew;
      if ( _configurationFlags & CALC_STOPPING_ENERGY_BEFORE_UPDATING_CURRENT_POSITION ) {
        _E.apply ( Dest, energyScalar );
        energy = energyScalar[0];
      }

      // Cached value needed by ArmijoLineSearchHelpFunction_evaluate.
      _energyAtPosition = energy;
      calcTaudAndDescentDir( Dest, tau, descentDir, *oldGradient, *oldDirection, _iterations );

      // Calculate the energy at the new position and update Dest now.
      {
        // writeEnergy assumes that the last evaluated energy is the energy of our new position
        // in case our energy is an aol::GenEnergyOp or an aol::LinCombOp consisting of aol::GenEnergyOps
        // So we don't reuse the energy calculations in case _out is not NULL.
        const bool energyPossiblyReusable = ( _timestepController == ARMIJO ) && ( tau.size() == 1 ) && ( _out == NULL );
        // Setting this to NaN allows us to check whether we could reuse the energy or have to calculate it.
        energyScalar[0] = aol::NumberTrait<RealType>::NaN;

        // The armijo step size control already had to calculate the energy at the new position,
        // either as last or second to last energy evaluation.
        if ( energyPossiblyReusable && ( _lastEnergyTimestepWidth[0] == tau[0] ) )
          energyScalar[0] = _energyForLastTimestepWidth[0];
        else if ( energyPossiblyReusable && ( _lastEnergyTimestepWidth[1] == tau[0] ) )
          energyScalar[0] = _energyForLastTimestepWidth[1];
        // If no step was made, the new energy is the same as the old one.
        else if ( ( tau.sum() == 0 ) && ( _out == NULL ) )
          energyScalar[0] = energy;

        // Update Dest and if necessary calculate the energy at the new position.
        // In case the derivate contained Inf or NaN, descentDir=tau*derivative
        // is possibly not zero even if tau is zero. Thus, we may only add
        // descentDir to Dest in case tau is not zero.
        if ( ( _configurationFlags & CALC_STOPPING_ENERGY_BEFORE_UPDATING_CURRENT_POSITION ) && ( tau.sum() != 0 ) ){
          descentDir += Dest;
          if ( aol::isNaN ( energyScalar[0] ) )
            _E.apply ( descentDir, energyScalar);
          Dest = descentDir;
        }
        else {
          if ( tau.sum() != 0 )
            Dest += descentDir;

          if ( aol::isNaN ( energyScalar[0] ) )
            _E.apply(Dest, energyScalar);
        }
      }
      energyNew = energyScalar[0];
#ifdef USE_CPP11
      std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
      _timeOneStep = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.;
#endif
      
      writeEnergy( energyNew );
      if ( !( _configurationFlags & DO_NOT_WRITE_CONSOLE_OUTPUT ) )
        writeConsoleOutput( energyNew, tau );
      
      // Save the current time step if _pStepSaver was initialized with setStepSaverReference
      if ( _pStepSaver )
        _pStepSaver->saveStep ( Dest, _iterations );

      // This is the old method for saving time steps, which is kept here for compatibility reasons
      if (this->checkSaveConditions(_iterations))
        this->writeTimeStep( Dest, _iterations );     // virtual function to call the TimestepSaver

      // stoppingCriterion == true means the iteration will be continued
      stoppingCriterion = ( ( ( _configurationFlags & USE_GRADIENT_BASED_STOPPING ) ? _pVecDE->norm() : energy - energyNew ) > this->_stopEpsilon ) && ( _iterations < _maxIterations );

      // If we are not updating the filter width and the current tau is zero, we have to stop now.
      // The current position of course will always lead to tau to be zero.
      // It's possible that this case occurs if the user has set _stopEpsilon lower than zero.
      if ( ( _updateFilterWidth == false ) && ( tau.sum() == 0 ) )
        stoppingCriterion = false;

      if ( (_updateFilterWidth == true)
           // If the new _filterWidth would be <= _stopFilterWidth, just stop the refining.
           && ( 0.5 * _filterWidth > _stopFilterWidth )
           && ( (energy - energyNew) < _updateFilterWidthTolerance
                 || ( (++iterationsPerFilterWidth) >= _maximumNumIterationsPerFilterWidth
                       && _maximumNumIterationsPerFilterWidth != 0 ) ) ){
        iterationsPerFilterWidth = 0;
        _filterWidth *= 0.5;
        tau.setAll( 1. );
        if( ( _filterWidth > _stopFilterWidth ) && ( _iterations < _maxIterations) )
          stoppingCriterion = true;
        else
          stoppingCriterion = false;
        if( _outFilterWidth != NULL )
          *(_outFilterWidth) << (_iterations)-1 << " " << energyNew << " " << _filterWidth << endl;
      }

      // After applying postProcess the energy has to be recalculated, otherwise the stopping criterion may fail in the next step.
      // We also need the recalculated energy for the optimization in ArmijoLineSearchHelpFunction_evaluate.
      // Relies on short-circuit evaluation, i.e. postProcess is only called if the first condition is true.
      if( ( ( _iterations % _postProcessingOffset ) == 0 ) && postProcess( Dest, _iterations ) ){
        _E.apply(Dest, energyScalar);                       // for the stopping criterion
        energyNew = energyScalar[0];
      }


    } while ( stoppingCriterion ); // end iteration loop

    if ( !( _configurationFlags & DO_NOT_WRITE_CONSOLE_OUTPUT ) ) {
      cerr << endl;
      cerr << "Descent needed " << _iterations << " step(s).\n";
    }

    // Save the energy at the final position.
    _energyAtPosition = energyNew;

    // Using zero as _startTau doesn't make sense.
    if ( tau[0] > 0 )
      this->_startTau = tau[0];
    delete _pVecDE;
  }
  virtual void apply( const VectorType &Arg, VectorType &Dest ) const {
    Dest = Arg;
    applySingle ( Dest );
  }
  virtual void applyAdd( const VectorType &/*Arg*/, VectorType &/*Dest*/ ) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
  void setStopEpsilon ( const RealType StopEpsilon ) {
    _stopEpsilon = StopEpsilon;
  }
  void setFilterWidth( const RealType FilterWidth ){
    _filterWidth = FilterWidth;
  }
  RealType getFilterWidth( ) const{
    return _filterWidth;
  }
  void setStopFilterWidth( const RealType StopFilterWidth ){
    _stopFilterWidth = StopFilterWidth;
  }
  void setUpdateFilterWidthTolerance( const RealType UpdateFilterWidthTolerance ){
    _updateFilterWidthTolerance = UpdateFilterWidthTolerance;
  }
  void setFilterWidthOutStream(std::ostream &Out){
    _outFilterWidth = &Out;
  }
  void setMaximumNumIterationsPerFilterWidth(const int MaximumNumIterationsPerFilterWidth){
    _maximumNumIterationsPerFilterWidth = MaximumNumIterationsPerFilterWidth;
  }
  void setTimestepController( const TIMESTEP_CONTROLLER TimestepController ) {
    _timestepController = TimestepController;
  }
  int getIterations() const{
    return _iterations;
  }
  int getMaxIterations() const{
    return _maxIterations;
  }
  RealType getLastTimestepWidth() const{
    // strange access due to the fact that push_back of RingBuffer is used and access via [-1] might not work because of modulo operation
    return _lastEnergyTimestepWidth[ _lastEnergyTimestepWidth.size()-1 ];
  }
  void setOutStream(std::ostream &Out){
    _out = &Out;
  }
  void setConfigurationFlags ( const unsigned int ConfigurationFlags ) {
    _configurationFlags = ConfigurationFlags;
  }
  void setPostProcessingOffset ( const unsigned int PostProcessingOffset ) {
    _postProcessingOffset = PostProcessingOffset;
  }
  RealType getEnergyAtLastPosition ( ) const {
    return _energyAtPosition;
  }
  void setStepSaverReference ( const StepSaverBase<RealType, VectorType>& StepSaver ) {
    _pStepSaver = &StepSaver;
  }
protected:
  void setUpdateFilterWidth( const bool UpdateFilterWidth ){
    _updateFilterWidth = UpdateFilterWidth;
  }
};


template <typename ConfiguratorType, typename VectorType>
class QuocGradientDescent
: public GradientDescentBase<typename ConfiguratorType::RealType, VectorType> {
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_grid;
public:
  QuocGradientDescent( const typename ConfiguratorType::InitType &Initializer,
                       const aol::Op<VectorType, aol::Scalar<RealType> > &E,
                       const aol::Op<VectorType> &DE,
                       const int MaxIterations = 50,
                       const RealType StartTau = aol::ZOTrait<RealType>::one,
                       const RealType StopEpsilon = aol::ZOTrait<RealType>::zero )
  : GradientDescentBase<RealType, VectorType>(E, DE, MaxIterations, StartTau, StopEpsilon),
  _grid(Initializer)
  {}
  virtual ~QuocGradientDescent() {}
  virtual void smoothDirection( const VectorType&, const RealType, VectorType&/*Direction*/ ) const {
  }
  virtual void writeTimeStep( const VectorType &Dest, const int Iterations ) const {
    this->saveTimestepVectorOrMultiVector( Iterations, Dest, _grid );
    this->saveTimestepBZ2( Iterations, Dest, _grid );
  }
};

/**
 * \ingroup Optimization
 */
template <typename ConfiguratorType, typename VectorType>
class GradientDescent{
};

template <typename ConfiguratorType>
class GradientDescent<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> >
: public QuocGradientDescent<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  int _smoothFirstNComponents;
public:
  GradientDescent( const typename ConfiguratorType::InitType &Initializer,
                   const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &E,
                   const aol::Op<aol::MultiVector<RealType> > &DE,
                   const int MaxIterations = 50,
                   const RealType StartTau = aol::ZOTrait<RealType>::one,
                   const RealType StopEpsilon = aol::ZOTrait<RealType>::zero )
    : QuocGradientDescent<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> >(Initializer, E, DE, MaxIterations, StartTau, StopEpsilon),
      _smoothFirstNComponents ( std::numeric_limits<int>::max() ) {
    this->setFilterWidth ( sqrt ( 10 * this->_grid.H() ) );
    this->setUpdateFilterWidth( false );
  }
  virtual ~GradientDescent() {}
  virtual void smoothDirection( const aol::MultiVector<RealType>&, const RealType Sigma, aol::MultiVector<RealType> &Direction ) const {
    qc::LinearSmoothOp<RealType> linSmooth;
    linSmooth.setCurrentGrid( this->_grid );
    linSmooth.setSigma( Sigma );

    const int smoothComponents = aol::Min ( Direction.numComponents(), _smoothFirstNComponents );

    for ( int i = 0; i < smoothComponents; i++ ){
      linSmooth.apply( Direction[i], Direction[i] );
    }
  }

  void setSmoothFirstNComponents ( const int SmoothFirstNComponents ) {
    _smoothFirstNComponents = SmoothFirstNComponents;
  }
};

template <typename ConfiguratorType>
class GradientDescent<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType> >
: public QuocGradientDescent<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
public:
  GradientDescent( const typename ConfiguratorType::InitType &Initializer,
                   const aol::Op<aol::Vector<RealType>, aol::Scalar<RealType> > &E,
                   const aol::Op<aol::Vector<RealType> > &DE,
                   const int MaxIterations = 50,
                   const RealType StartTau = aol::ZOTrait<RealType>::one,
                   const RealType StopEpsilon = aol::ZOTrait<RealType>::zero )
    : QuocGradientDescent<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType> >(Initializer, E, DE, MaxIterations, StartTau, StopEpsilon) {
    this->setFilterWidth ( sqrt ( 10 * this->_grid.H() ) );
    this->setUpdateFilterWidth( false );
  }
  virtual ~GradientDescent () {}
  virtual void smoothDirection( const aol::Vector<RealType>&, const RealType Sigma, aol::Vector<RealType> &Direction ) const {
    qc::LinearSmoothOp<RealType> linSmooth;
    linSmooth.setCurrentGrid( this->_grid );
    linSmooth.setSigma( Sigma );
    linSmooth.apply( Direction, Direction );
  }
};


/**
 * \brief Nearly identical to the class GradientDescent except that it uses
 *        qc::CGBasedInverseH1Metric instead of qc::LinearSmoothOp to apply
 *        the H1 metric to get the descent direction. Because of this, it
 *        also works on rectangular grids.
 *
 * Changing the template InverseH1MetricType replaces qc::CGBasedInverseH1Metric
 * with any desired smooth op, like qc::LinearSmoothOp.
 *
 * \author Berkels
 * \ingroup Optimization
 */
template <typename ConfiguratorType, typename VectorType, typename InverseH1MetricType = qc::CGBasedInverseH1Metric<ConfiguratorType> >
class H1GradientDescent
: public QuocGradientDescent<ConfiguratorType, VectorType > {
  typedef typename ConfiguratorType::RealType RealType;
  mutable InverseH1MetricType _smoothOp;
  int _smoothFirstNComponents;
public:
  H1GradientDescent( const typename ConfiguratorType::InitType &Initializer,
                     const aol::Op<VectorType, aol::Scalar<RealType> > &E,
                     const aol::Op<VectorType> &DE,
                     const int MaxIterations = 50,
                     const RealType StartTau = aol::ZOTrait<RealType>::one,
                     const RealType StopEpsilon = aol::ZOTrait<RealType>::zero )
    : QuocGradientDescent<ConfiguratorType, VectorType>(Initializer, E, DE, MaxIterations, StartTau, StopEpsilon),
      _smoothOp ( Initializer ),
      _smoothFirstNComponents ( std::numeric_limits<int>::max() ) {
    this->setFilterWidth ( sqrt ( 10 * this->_grid.H() ) );
    this->setUpdateFilterWidth( false );
  }
  virtual ~H1GradientDescent () {}
  virtual void smoothDirection( const aol::Vector<RealType>&, const RealType Sigma, aol::Vector<RealType> &Direction ) const {
    _smoothOp.setSigma ( Sigma );
    _smoothOp.apply( Direction, Direction );
  }
  virtual void smoothDirection( const aol::MultiVector<RealType>&, const RealType Sigma, aol::MultiVector<RealType> &Direction ) const {
    _smoothOp.setSigma ( Sigma );

    const int smoothComponents = aol::Min ( Direction.numComponents(), _smoothFirstNComponents );
    for ( int i = 0; i < smoothComponents; i++ )
      _smoothOp.apply( Direction[i], Direction[i] );
  }
  void setSmoothFirstNComponents ( const int SmoothFirstNComponents ) {
    _smoothFirstNComponents = SmoothFirstNComponents;
  }
};

/**
 * \ingroup Optimization
 */
template <typename ConfiguratorType, typename VectorType, typename GradientDescentType = GradientDescent<ConfiguratorType, VectorType> >
class GradientDescentWithAutomaticFilterWidth{
};

template <typename ConfiguratorType, typename GradientDescentType>
class GradientDescentWithAutomaticFilterWidth<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType>, GradientDescentType>
: public GradientDescentType {
  typedef typename ConfiguratorType::RealType RealType;
public:
  GradientDescentWithAutomaticFilterWidth( const typename ConfiguratorType::InitType &Initializer,
                                           const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &E,
                                           const aol::Op<aol::MultiVector<RealType> > &DE,
                                           const int MaxIterations = 50,
                                           const RealType StartTau = aol::ZOTrait<RealType>::one,
                                           const RealType StopEpsilon = aol::ZOTrait<RealType>::zero )
  : GradientDescentType(Initializer, E, DE, MaxIterations, StartTau, StopEpsilon)
  {
    this->setFilterWidth ( aol::ZOTrait<RealType>::one );
    this->setUpdateFilterWidth( true );
  }
};

template <typename ConfiguratorType, typename GradientDescentType>
class GradientDescentWithAutomaticFilterWidth<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType>, GradientDescentType>
: public GradientDescentType {
  typedef typename ConfiguratorType::RealType RealType;
public:
  GradientDescentWithAutomaticFilterWidth( const typename ConfiguratorType::InitType &Initializer,
                                           const aol::Op<aol::Vector<RealType>, aol::Scalar<RealType> > &E,
                                           const aol::Op<aol::Vector<RealType> > &DE,
                                           const int MaxIterations = 50,
                                           const RealType StartTau = aol::ZOTrait<RealType>::one,
                                           const RealType StopEpsilon = aol::ZOTrait<RealType>::zero )
  : GradientDescentType(Initializer, E, DE, MaxIterations, StartTau, StopEpsilon)
  {
    this->setFilterWidth ( aol::ZOTrait<RealType>::one );
    this->setUpdateFilterWidth( true );
  }
};

template <typename ConfiguratorType, typename GradientDescentType = GradientDescent<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> > >
class GradientDescentComponentWiseTimestepControlled
: public GradientDescentType {
  typedef typename ConfiguratorType::RealType RealType;
protected:
  //! components of the MultiVector apply argument can be grouped to share the same tau with _customTauBlocks
  aol::Vector<int> _customTauBlocks;
public:
  GradientDescentComponentWiseTimestepControlled( const typename ConfiguratorType::InitType &Initializer,
                                                  const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &E,
                                                  const aol::Op<aol::MultiVector<RealType> > &DE,
                                                  const int MaxIterations = 50,
                                                  const RealType StartTau = aol::ZOTrait<RealType>::one,
                                                  const RealType StopEpsilon = aol::ZOTrait<RealType>::zero )
  : GradientDescentType(Initializer, E, DE, MaxIterations, StartTau, StopEpsilon), _customTauBlocks( 0 )
  {
  }

  GradientDescentComponentWiseTimestepControlled( const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &E,
                                                  const aol::Op<aol::MultiVector<RealType> > &DE,
                                                  const int MaxIterations = 50,
                                                  const RealType StartTau = aol::ZOTrait<RealType>::one,
                                                  const RealType StopEpsilon = aol::ZOTrait<RealType>::zero )
  : GradientDescentType(E, DE, MaxIterations, StartTau, StopEpsilon), _customTauBlocks( 0 )
  {
  }

  GradientDescentComponentWiseTimestepControlled( const typename ConfiguratorType::InitType &Initializer,
                                                  const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &E,
                                                  const aol::Op<aol::MultiVector<RealType> > &DE,
                                                  const aol::Vector<int> &TauBlocks,
                                                  const int MaxIterations = 50,
                                                  const RealType StartTau = aol::ZOTrait<RealType>::one,
                                                  const RealType StopEpsilon = aol::ZOTrait<RealType>::zero )
  : GradientDescentType(Initializer, E, DE, MaxIterations, StartTau, StopEpsilon), _customTauBlocks( TauBlocks ) {}

  virtual ~GradientDescentComponentWiseTimestepControlled () {}

  virtual void getTauAndUpdateDescentDir( aol::MultiVector<RealType> &DescentDir,
                                          const aol::MultiVector<RealType> &CurrentPosition,
                                          aol::Vector<RealType> &Tau ) const {
    aol::Vector<int> tauBlocks(0);
    if ( _customTauBlocks.size() )
    {
      tauBlocks.resize( _customTauBlocks.size() );
      tauBlocks = _customTauBlocks;
    }
    else
    {
      tauBlocks.resize( DescentDir.numComponents() );
      tauBlocks.setAll( 1 );
    }
    const int numComponents = tauBlocks.size();
    if( numComponents != Tau.size() ){
      Tau.resize( numComponents );
      Tau.setAll( 1. );
    }
    aol::MultiVector<RealType> descentDirComponent( DescentDir, aol::STRUCT_COPY );

    int handledDDComponents = 0;
    for( int component = 0; component < numComponents; component++){
      for( int i = 0; i < tauBlocks[component]; i++)
        descentDirComponent[i+handledDDComponents] = DescentDir[i+handledDDComponents];
      Tau[component] = this->getStepSizeWithSelectedController(descentDirComponent, CurrentPosition, Tau[component]);
      for( int i = 0; i < tauBlocks[component]; i++){
        DescentDir[i+handledDDComponents] *= Tau[component];
        descentDirComponent[i+handledDDComponents].setZero();
      }
      handledDDComponents += tauBlocks[component];
    }
  }
};

/**
 * \author Berkels
 * \ingroup Optimization
 */
template <typename ConfiguratorType, typename VectorType>
class L2GradientDescentDirichletBCs
: public QuocGradientDescent<ConfiguratorType, VectorType> {
  typedef typename ConfiguratorType::RealType RealType;
  // In case no Dirichlet mask is supplied in the contructor, this pointer
  // will be used to store an automatically generated Dirichlet mask.
  aol::DeleteFlagPointer<typename qc::BitArray<ConfiguratorType::Dim> > _pDirichletMaskStorage;
  const typename qc::BitArray<ConfiguratorType::Dim> &_dirichletMask;
  typename ConfiguratorType::MatrixType _massMatMasked;
  bool _applyMetric;
public:
  L2GradientDescentDirichletBCs ( const typename ConfiguratorType::InitType &Initializer,
                                  const aol::Op<VectorType, aol::Scalar<RealType> > &E,
                                  const aol::Op<VectorType> &DE,
                                  const typename qc::BitArray<ConfiguratorType::Dim> &DirichletMask,
                                  const int MaxIterations = 50,
                                  const RealType StartTau = aol::ZOTrait<RealType>::one,
                                  const RealType StopEpsilon = aol::ZOTrait<RealType>::zero )
  : QuocGradientDescent<ConfiguratorType, VectorType>(Initializer, E, DE, MaxIterations, StartTau, StopEpsilon),
    _pDirichletMaskStorage ( ),
    _dirichletMask ( DirichletMask ),
    _massMatMasked ( Initializer ),
    _applyMetric ( true ) {
    aol::MassOp<ConfiguratorType> massOp ( Initializer );
    massOp.assembleAddMatrix ( _massMatMasked, &_dirichletMask, true );
  }

  //! All boundary nodes of the grid will be treated as Dirichlet nodes.
  L2GradientDescentDirichletBCs ( const typename ConfiguratorType::InitType &Initializer,
                                  const aol::Op<VectorType, aol::Scalar<RealType> > &E,
                                  const aol::Op<VectorType> &DE,
                                  const int MaxIterations = 50,
                                  const RealType StartTau = aol::ZOTrait<RealType>::one,
                                  const RealType StopEpsilon = aol::ZOTrait<RealType>::zero )
  : QuocGradientDescent<ConfiguratorType, VectorType>(Initializer, E, DE, MaxIterations, StartTau, StopEpsilon),
    _pDirichletMaskStorage ( new typename qc::BitArray<ConfiguratorType::Dim> ( Initializer ), true ),
    _dirichletMask ( *_pDirichletMaskStorage ),
    _massMatMasked ( Initializer ),
    _applyMetric ( true ) {

    for ( qc::GridStructure::AllBoundaryNodeIterator it = this->_grid; it.notAtEnd(); ++it )
      _pDirichletMaskStorage->set ( *it, true );

    aol::MassOp<ConfiguratorType> massOp ( Initializer );
    massOp.assembleAddMatrix ( _massMatMasked, &_dirichletMask, true );
  }

  virtual ~L2GradientDescentDirichletBCs () {
  }

  void setApplyMetric ( const bool ApplyMetric ) {
    _applyMetric = ApplyMetric;
  }

  virtual void smoothDirection( const aol::Vector<RealType>&, const RealType, aol::Vector<RealType> &Direction ) const {
    for ( int i = 0; i < _dirichletMask.size(); ++i ) {
      if ( _dirichletMask[i] == true ) {
        Direction[i] = 0;
      }
    }
    if ( _applyMetric ) {
      aol::CGInverse< aol::Vector<RealType> > solver ( _massMatMasked );
      solver.setStopping ( aol::STOPPING_ABSOLUTE );
      solver.setQuietMode ( true );
      aol::Vector<RealType> temp ( Direction, aol::STRUCT_COPY );
      solver.apply ( Direction, temp );
      Direction = temp;
    }
  }

  virtual void smoothDirection( const aol::MultiVector<RealType>&, const RealType, aol::MultiVector<RealType> &Direction ) const {
    for ( int i = 0; i < _dirichletMask.size(); ++i ) {
      if ( _dirichletMask[i] == true ) {
        for ( int j = 0; j < Direction.numComponents(); ++j )
          Direction[j][i] = 0;
      }
    }
    if ( _applyMetric ) {
      aol::CGInverse< aol::Vector<RealType> > solver ( _massMatMasked );
      solver.setStopping ( aol::STOPPING_ABSOLUTE );
      solver.setQuietMode ( true );
      aol::Vector<RealType> temp ( Direction[0], aol::STRUCT_COPY );

      for ( int j = 0; j < Direction.numComponents(); ++j ) {
        solver.apply ( Direction[j], temp );
        Direction[j] = temp;
      }
    }
  }
};

/**
 * \ingroup Optimization
 */
template <typename RealType, typename VectorType, typename DerivativeType = VectorType>
class GridlessGradientDescent
: public GradientDescentBase<RealType, VectorType, DerivativeType> {
public:
  GridlessGradientDescent( const aol::Op<VectorType, aol::Scalar<RealType> > &E,
                           const aol::Op<VectorType, DerivativeType> &DE,
                           const int MaxIterations = 50,
                           const RealType StartTau = aol::ZOTrait<RealType>::one,
                           const RealType StopEpsilon = aol::ZOTrait<RealType>::zero )
  : GradientDescentBase<RealType, VectorType, DerivativeType>(E, DE, MaxIterations, StartTau, StopEpsilon) {}

  //! Contructor with an added dummy argument that allows this class to be used as template in qc::StandardRegistration.
  template <typename GridType>
  GridlessGradientDescent( const GridType &/*Dummy*/,
                           const aol::Op<VectorType, aol::Scalar<RealType> > &E,
                           const aol::Op<VectorType, DerivativeType> &DE,
                           const int MaxIterations = 50,
                           const RealType StartTau = aol::ZOTrait<RealType>::one,
                           const RealType StopEpsilon = aol::ZOTrait<RealType>::zero )
  : GradientDescentBase<RealType, VectorType, DerivativeType>(E, DE, MaxIterations, StartTau, StopEpsilon) {}
  virtual ~GridlessGradientDescent() {}
  virtual void smoothDirection( const VectorType&, const RealType, DerivativeType &/*Direction*/ ) const {}
};

template <typename ConfiguratorType, typename VectorType>
class SimpleGradientDescent{
};

template <typename ConfiguratorType>
class SimpleGradientDescent<ConfiguratorType, aol::MultiVector<typename ConfiguratorType::RealType> >
: public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
  const typename ConfiguratorType::InitType &_grid;
  typedef typename ConfiguratorType::RealType RealType;
  const aol::Op<aol::MultiVector<RealType> > &_DE;
  const int _maxIterations;
public:
  SimpleGradientDescent(const typename ConfiguratorType::InitType &Initializer,
                        const aol::Op<aol::MultiVector<RealType> > &DE,
                        const int MaxIterations = 50,
                        const RealType /*StartTau*/ = aol::ZOTrait<RealType>::one)
  : _grid(Initializer),
    _DE(DE),
    _maxIterations(MaxIterations){
  }
  virtual ~SimpleGradientDescent () {
  }
  virtual void apply( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const{
    aol::MultiVector<RealType> descentDir(ConfiguratorType::Dim, _grid.getSize());
    Dest = Arg;
    int iterations = 0;

    do{
      iterations++;

      // old code (use iterations < 1 to replicate old behaviour)
      aol::MultiVector<RealType> mtmp(ConfiguratorType::Dim, _grid.getSize());
      _DE.apply(Dest, mtmp);

      // inverse mass matrix here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      const aol::LumpedMassOp<ConfiguratorType> lumpedMassInv( _grid, true );
      for ( int i = 0; i < ConfiguratorType::Dim; i++ ){
        lumpedMassInv.apply( mtmp[i], descentDir[i] );
      }
      descentDir *= -1.;

      qc::LinearSmoothOp<RealType> linSmooth;
      linSmooth.setCurrentGrid( _grid );
      linSmooth.setTau( 5 * _grid.H() );

      for ( int i = 0; i < ConfiguratorType::Dim; i++ ){
        linSmooth.apply( descentDir[i], descentDir[i] );
      }

      aol::Vector<RealType> maxVec(2*ConfiguratorType::Dim);
      for ( int i = 0; i < ConfiguratorType::Dim; i++ ){
        maxVec[i] = aol::Abs( descentDir[i].getMaxValue() );
        maxVec[i+ConfiguratorType::Dim] = aol::Abs( descentDir[i].getMinValue() );
      }
      RealType max = maxVec.getMaxValue();

      descentDir *= _grid.H() / max;
      Dest += descentDir;
      // end old code
    } while (( iterations < _maxIterations));
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
};

template <typename ConfiguratorType>
class SimpleGradientDescent<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType> >
: public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {
  const typename ConfiguratorType::InitType &_grid;
  typedef typename ConfiguratorType::RealType RealType;
  const aol::Op<aol::Vector<RealType> > &_DE;
  const int _maxIterations;
public:
  SimpleGradientDescent(const typename ConfiguratorType::InitType &Initializer,
                        const aol::Op<aol::Vector<RealType> > &DE,
                        const int MaxIterations = 50,
                        const RealType /*StartTau*/ = aol::ZOTrait<RealType>::one)
  : _grid(Initializer),
    _DE(DE),
    _maxIterations(MaxIterations){
  }
  virtual ~SimpleGradientDescent () {}
  virtual void apply( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const{
    aol::Vector<RealType> descentDir(_grid.getSize());
    Dest = Arg;
    int iterations = 0;

    do{
      iterations++;

      // old code (use iterations < 1 to replicate old behaviour)
      aol::Vector<RealType> mtmp(_grid.getSize());
      _DE.apply(Dest, mtmp);

      // inverse mass matrix here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      const aol::LumpedMassOp<ConfiguratorType> lumpedMassInv( _grid, true );
      lumpedMassInv.apply( mtmp, descentDir );
      descentDir *= -1.;

      qc::LinearSmoothOp<RealType> linSmooth;
      linSmooth.setCurrentGrid( _grid );
      linSmooth.setTau( 5 * _grid.H() );

      linSmooth.apply( descentDir, descentDir );

      aol::Vector<RealType> maxVec(2);
      maxVec[0] = aol::Abs( descentDir.getMaxValue() );
      maxVec[1] = aol::Abs( descentDir.getMinValue() );
      RealType max = maxVec.getMaxValue();

      descentDir *= _grid.H() / max;
      Dest += descentDir;
      // end old code
    } while (( iterations < _maxIterations));
  }
  virtual void applyAdd( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
};

template <typename RealType, typename DomainType, typename DerivativeType=DomainType>
class DifferentialQuotientTestBase{
private:
  virtual void shiftIthDOFByEpsilon( const RealType Epsilon, const int I, DomainType &Position) const = 0;
  virtual RealType getDerivativeWRTIthDOF( const int I, const DomainType &Position) const = 0;
public:
  const aol::Op<DomainType, aol::Scalar<RealType> > &_E;
  const aol::Op<DomainType, DerivativeType > &_DE;
  const int _numDofs;
  mutable DerivativeType _vecDE;
public:
  DifferentialQuotientTestBase(const aol::Op<DomainType, aol::Scalar<RealType> > &E,
                               const aol::Op<DomainType, DerivativeType > &DE,
                               const int NumDofs,
                               const DerivativeType &VecDE)
  : _E(E),
    _DE(DE),
    _numDofs(NumDofs),
    _vecDE(VecDE)
  {}
  virtual ~DifferentialQuotientTestBase(){}

  void test( const RealType Epsilon, const int I, const DomainType &Position ) const{
    DomainType shiftedPosition(Position);
    aol::Scalar<RealType> energyScalar;
    _E.apply(Position, energyScalar);
    const RealType energyAtPosition = energyScalar[0];
    shiftIthDOFByEpsilon( Epsilon, I, shiftedPosition );
    _E.apply(shiftedPosition, energyScalar);
    const RealType energyAtShiftedPosition = energyScalar[0];
    const RealType diffQuotient = (energyAtShiftedPosition - energyAtPosition)/Epsilon;
    const RealType derivative = getDerivativeWRTIthDOF( I, Position );
    cerr << "DQ    = " << diffQuotient << endl;
    cerr << "Deriv = " << derivative << endl;
    cerr << "Diff  = " << diffQuotient - derivative << endl;
    aol::callSystemPauseIfNecessaryOnPlatform();
  }
  void testAllDofs( const RealType Epsilon, const DomainType &Position ) const{
    for( int i = 0; i < _numDofs; i++){
      cerr << "Current Dof " << i << endl;
      test( Epsilon, i, Position);
    }
  }
};

template <typename RealType, typename DomainType, typename DerivativeType=DomainType>
class DifferentialQuotientTest{
};


template <typename RealType>
class DifferentialQuotientTest<RealType, aol::MultiVector<RealType> >
: public DifferentialQuotientTestBase<RealType, aol::MultiVector<RealType> > {
public:
  DifferentialQuotientTest(const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &E,
                           const aol::Op<aol::MultiVector<RealType> > &DE,
                           const int NumComponents,
                           const int NumPoints)
  : DifferentialQuotientTestBase<RealType, aol::MultiVector<RealType> >(E, DE, NumComponents*NumPoints, aol::MultiVector<RealType>(NumComponents, NumPoints) )
  {
  }
  virtual ~DifferentialQuotientTest ()
  {
  }
  virtual void shiftIthDOFByEpsilon( const RealType Epsilon, const int I, aol::MultiVector<RealType> &Position) const{
    const int i = I%Position.numComponents();
    const int j = I/Position.numComponents();
    Position[i][j] += Epsilon;
  }
  virtual RealType getDerivativeWRTIthDOF( const int I, const aol::MultiVector<RealType> &Position) const{
    const int i = I%Position.numComponents();
    const int j = I/Position.numComponents();
    this->_DE.apply( Position, this->_vecDE);
    return this->_vecDE[i][j];
  };
};

template <typename RealType>
 class DifferentialQuotientTest<RealType, aol::Vector<RealType> >
 : public aol::DifferentialQuotientTestBase<RealType, aol::Vector<RealType> > {
 public:
   DifferentialQuotientTest(const aol::Op<aol::Vector<RealType>, aol::Scalar<RealType> > &E,
                            const aol::Op<aol::Vector<RealType> > &DE,
                            const int NumPoints)
   : aol::DifferentialQuotientTestBase<RealType, aol::Vector<RealType> >(E, DE, NumPoints, aol::Vector<RealType>(NumPoints) )
   {
   }
   virtual ~DifferentialQuotientTest() {
   }
   virtual void shiftIthDOFByEpsilon( const RealType Epsilon, const int I, aol::Vector<RealType> &Position) const{
     Position[I] += Epsilon;
   }
   virtual RealType getDerivativeWRTIthDOF( const int I, const aol::Vector<RealType> &Position) const{
     this->_DE.apply( Position, this->_vecDE);
     return this->_vecDE[I];
   };
 };

} // namespace aol

#endif // __GRADIENTDESCENT_H
