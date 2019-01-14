#ifndef __LINEARSMOOTHOP_H
#define __LINEARSMOOTHOP_H

#include <quoc.h>
#include <multilevelArray.h>
#include <solver.h>
#include <preconditioner.h>
#include <suiteSparseSolver.h>

// forward declaration of nb::NarrowGridBase
namespace nb {
  template <typename FullGridType> class NarrowBandGridBase;
}

namespace qc {

/**
 * \author Berkels
 */
template <typename _RealType>
class DyadicGridTrait {
public:
  typedef _RealType RealType;
  typedef qc::GridDefinition GridType;
  typedef qc::ProlongOp<RealType> ProlongOpType;
  typedef qc::RestrictOp<RealType, qc::STD_QUOC_RESTRICT> RestrictOpType;
  typedef qc::MultilevelArray<RealType, qc::Array<RealType>, ProlongOpType, RestrictOpType, GridType> MultilevelArrayType;
};

/** class for multigrid linear smoother via the heat equation
 * \author Droske
 */
template <typename RealType, typename GridTrait = qc::DyadicGridTrait<RealType> >
class LinearSmoothOp : public aol::Op< aol::Vector<RealType> > {
  typedef typename GridTrait::GridType GridType;
  typedef typename GridTrait::MultilevelArrayType MultilevelArrayType;
public:
  LinearSmoothOp() : _grid ( NULL, false ), _x ( NULL ), _rhs ( NULL ), _tau ( 1.0 ) {}

  template <typename InputGridType>
  LinearSmoothOp( const InputGridType &Grid ) : _grid ( NULL, false ), _x ( NULL ), _rhs ( NULL ), _tau ( 1.0 ) {
    setCurrentGrid ( Grid );
  }

  virtual ~LinearSmoothOp() {
    delete _x;
    delete _rhs;
  }

  void setSigma ( RealType Sigma ) {
    _tau = 0.5f * aol::Sqr ( Sigma );
  }

  void setTau ( RealType Tau ) {
    _tau = Tau;
  }

  void setCurrentGrid ( const GridType &Grid ) {
    _grid.reset ( &Grid, false );
    resetMultilevelArrays ();
  }

  void resetMultilevelArrays () {
    delete _x;
    _x = new MultilevelArrayType ( _grid->getGridDepth(), _grid->getDimOfWorld() );

    delete _rhs;
    _rhs = new MultilevelArrayType ( _grid->getGridDepth(), _grid->getDimOfWorld() );
  }

  /**
   * Tries to convert the given qc::RectangularGrid to GridType. If this is
   * possible, this function can be used exactly like the setCurrentGrid version that gets a
   * GridType.
   *
   * Otherwise it still serves as dummy function to allow code involving this class to be
   * compiled for qc::RectangularGrid. For example useful for classes derived from
   * aol::GradientDescent that don't need smoothing in their implementation of smoothDirection.
   */
  template <qc::Dimension Dim>
  void setCurrentGrid ( const qc::RectangularGrid<Dim> &Grid ) {
    _grid.reset ( new GridType ( Grid.getSize() ), true );
    resetMultilevelArrays ();
  }

  /**
   * This function can be used exactly like the setCurrentGrid version that gets a
   * qc::GridDefinition.
   * Furthermore, it serves as dummy function to allow code involving this class to be
   * compiled for qc::simplex::GridStructure. For example useful for classes derived from
   * aol::GradientDescent that don't need smoothing in their implementation of smoothDirection.
   */
  template <qc::Dimension Dim>
  void setCurrentGrid ( const qc::simplex::GridStructure<qc::GridDefinition,Dim> &Grid ) {
    _grid.reset ( &Grid.getCubicGrid(), true );
    resetMultilevelArrays ();
  }

  /**
   * Dummy function to allow code involving this class to be compiled for nb::NarrowGrid.
   * For example useful for classes derived from aol::GradientDescent that don't need
   * smoothing in their implementation of smoothDirection.
   */
  template <typename FullGridType>
  void setCurrentGrid ( const nb::NarrowBandGridBase<FullGridType> &/*Grid*/ ) {
    throw aol::Exception ( "LinearSmoothOp doesn't support nb::NarrowBandGrid", __FILE__, __LINE__ );
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp ( Dest, aol::DEEP_COPY );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const;

  void applySingle ( aol::Vector<RealType> &Arg ) const {
    aol::Vector<RealType> tmp ( Arg, aol::STRUCT_COPY );
    apply ( Arg, tmp );
    Arg = tmp;
  }

  void apply ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    // The apply for aol::Vector is not thread safe, so we can't just inherit the apply for aol::MultiVector from aol::BiOp
    for ( int i = 0; i < Arg.numComponents(); ++i )
      apply ( Arg[i], Dest[i] );
  }

  void applySingle ( aol::MultiVector<RealType> &Arg ) const {
    for ( int i = 0; i < Arg.numComponents(); ++i )
      applySingle ( Arg[i] );
  }

protected:
  aol::DeleteFlagPointer<const GridType> _grid;
  MultilevelArrayType *_x, *_rhs;
  RealType _tau;
};


/** Class to compute (approximately) one step of heat conduction \f$ (M + \tau L)^{-1} \f$
 *  using a PCG solver. In contrast to the LinearSmoothOp, it also works for non-cubic grids
 *  \author Schwen
 */
template <typename ConfiguratorType>
class GeneralLinearSmoothOp : public aol::BiOp< aol::Vector<typename ConfiguratorType::RealType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  const typename ConfiguratorType::InitType &_grid;
  const aol::StiffOp<ConfiguratorType> _stiffOp;
  typename ConfiguratorType::MatrixType _massMat, _systemMat;
  aol::DiagonalPreconditioner< aol::Vector<RealType> >* _prec;
  RealType _tau, _solverAccuracy;
  int _solverSteps;

public:
  // this version must be used with adaptive grids
  GeneralLinearSmoothOp( const ConfiguratorType& Config, const RealType Sigma )
    : _grid ( Config.getInitializer () ),
      _stiffOp ( Config, _grid ),
      _massMat ( _grid ),
      _systemMat ( _grid ),
      _prec ( NULL ),
      _tau ( 0 ),
      _solverAccuracy ( 1.0e-16 ),
      _solverSteps ( 500 ) {
    const aol::LumpedMassOp<ConfiguratorType> massOp ( Config, aol::DO_NOT_INVERT );
    massOp.assembleAddMatrix ( _massMat );
    setSigma ( Sigma ); // also assembles matrices
  }

  GeneralLinearSmoothOp( const typename ConfiguratorType::InitType &Initializer, const RealType Sigma )
    : _grid ( Initializer ),
      _stiffOp ( _grid ),
      _massMat ( _grid ),
      _systemMat ( _grid ),
      _prec ( NULL ),
      _tau ( 0 ),
      _solverAccuracy ( 1.0e-16 ),
      _solverSteps ( 500 ) {
    const aol::LumpedMassOp<ConfiguratorType> massOp ( _grid, aol::DO_NOT_INVERT );
    massOp.assembleAddMatrix ( _massMat );
    setSigma ( Sigma ); // also assembles matrices
  }

  GeneralLinearSmoothOp( const typename ConfiguratorType::InitType &Initializer )
    : _grid ( Initializer ),
      _stiffOp ( _grid ),
      _massMat ( _grid ),
      _systemMat ( _grid ),
      _prec ( NULL ),
      _tau ( aol::NumberTrait<RealType>::NaN ),
      _solverAccuracy ( 1.0e-16 ),
      _solverSteps ( 500 ) {
    const aol::LumpedMassOp<ConfiguratorType> massOp ( _grid, aol::DO_NOT_INVERT );
    massOp.assembleAddMatrix ( _massMat );
  }

  virtual ~GeneralLinearSmoothOp() {
    if ( _prec != NULL )
      delete ( _prec );
  }

  void setSigma ( const RealType Sigma ) {
    setTau ( aol::Sqr ( Sigma ) / 2 );
  }

  void setTau ( const RealType Tau ) {
    _tau = Tau;
    assembleSystemMatrix();
  }

  void setSolverParameters ( const RealType Eps, const int Steps ) {
    _solverAccuracy = Eps;
    _solverSteps = Steps;
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp ( Dest, aol::DEEP_COPY );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::PCGInverse< aol::Vector<RealType> > inv ( _systemMat, *_prec, _solverAccuracy, _solverSteps );
    inv.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    //    inv.setQuietMode ( true );
    // Makes sure that one can use apply ( a, a ).
    aol::Vector<RealType> rhs ( Arg, aol::STRUCT_COPY );
    _massMat.apply ( Arg, rhs );
    Dest = Arg; // better initial guess than zero ...
    inv.apply ( rhs, Dest );
  }

  void applySingle ( aol::Vector<RealType> &Dest ) const {
    apply ( Dest, Dest );
  }

  using aol::BiOp< aol::Vector<RealType> >::applyAdd;
  using aol::BiOp< aol::Vector<RealType> >::apply;
  using aol::BiOp< aol::Vector<RealType> >::applySingle;

protected:
  void assembleSystemMatrix () {
    _systemMat = _massMat;
    _stiffOp.assembleAddMatrix( _systemMat, _tau );
    if ( _prec != NULL ) delete ( _prec );
    _prec = new aol::DiagonalPreconditioner< aol::Vector<RealType> > ( _systemMat );
  }

private:
  GeneralLinearSmoothOp ( const GeneralLinearSmoothOp<ConfiguratorType> &other );
  GeneralLinearSmoothOp<ConfiguratorType>& operator= ( const GeneralLinearSmoothOp<ConfiguratorType> &other );
};


/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class InverseH1MetricBase : public aol::BiOp< aol::Vector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
public:
  InverseH1MetricBase( const typename ConfiguratorType::InitType &Initializer )
    : _grid ( Initializer ),
      _massOp ( _grid ),
      _stiffOp ( _grid ),
      _systemMat ( _grid ),
      _tau ( 1.0 ) {
    assembleSystemMatrix();
  }

  virtual ~InverseH1MetricBase() {}

  void setSigma ( const RealType Sigma ) {
    setTau ( 0.5f * aol::Sqr ( Sigma ) );
  }

  void setTau ( const RealType Tau ) {
    if ( _tau != Tau ) {
      _tau = Tau;
      assembleSystemMatrix();
    }
  }

protected:
  virtual void assembleSystemMatrix () {
    _systemMat.setZero();
    _massOp.assembleAddMatrix ( _systemMat );
    _stiffOp.assembleAddMatrix( _systemMat, _tau );
  }

  const typename ConfiguratorType::InitType &_grid;
  const aol::MassOp<ConfiguratorType> _massOp;
  const aol::StiffOp<ConfiguratorType> _stiffOp;
  typename ConfiguratorType::MatrixType _systemMat;
  RealType _tau;
};

/**
 * Class to apply \f$ (M + \tau L)^{-1} \f$ approximately with a CG solver and a relatively
 * big stopping epsilon.
 *
 * \note apply is not equivalent to one heat equation step with stepsize tau, because
 * the argument is not multiplied with M.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class CGBasedInverseH1Metric : public qc::InverseH1MetricBase<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;
public:
  CGBasedInverseH1Metric( const typename ConfiguratorType::InitType &Initializer )
    : qc::InverseH1MetricBase<ConfiguratorType> ( Initializer ) {}

  virtual ~CGBasedInverseH1Metric() {}

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp ( Dest, aol::DEEP_COPY );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Arg.size() != this->_grid.getNumberOfNodes() )
      throw aol::Exception ( "CGBasedInverseH1Metric::apply: Size of the argument doesn't match the size of the grid.", __FILE__, __LINE__ );

    aol::CGInverse<aol::Vector<RealType> > inv ( this->_systemMat, 1e-8 );
    inv.setStopping ( aol::STOPPING_ABSOLUTE );
    inv.setQuietMode ( true );
    // Makes sure that one can use apply ( a, a ).
    aol::Vector<RealType> rhs ( Arg, aol::DEEP_COPY );
    inv.apply ( rhs, Dest );
  }

  using aol::BiOp<aol::Vector<RealType> >::applyAdd;
  using aol::BiOp<aol::Vector<RealType> >::apply;
};


/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class CholeskyBasedInverseH1Metric : public qc::InverseH1MetricBase<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;
  aol::BlockOp<RealType, typename ConfiguratorType::MatrixType> _blockSystemMat;
  aol::CholeskyBlockInverseOp<RealType, typename ConfiguratorType::MatrixType> _choleskySolver;
  aol::MVecToVecOp<RealType> _inv;
public:
  CholeskyBasedInverseH1Metric( const typename ConfiguratorType::InitType &Initializer )
    : qc::InverseH1MetricBase<ConfiguratorType> ( Initializer ),
      _blockSystemMat( 1, 1 ),
      _inv ( _choleskySolver ) {
    _blockSystemMat.setReference( 0, 0, this->_systemMat );
  }

  virtual void assembleSystemMatrix () {
    qc::InverseH1MetricBase<ConfiguratorType>::assembleSystemMatrix ( );
    _choleskySolver.setMatrix ( _blockSystemMat, this->_systemMat.maxNumNonZeroesPerRow() );
  }

  virtual ~CholeskyBasedInverseH1Metric() {}

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    _inv.applyAdd ( Arg, Dest );
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    // Makes sure that one can use apply ( a, a ).
    aol::Vector<RealType> rhs ( Arg, aol::DEEP_COPY );
    _inv.apply ( rhs, Dest );
  }

  using aol::BiOp<aol::Vector<RealType> >::applyAdd;
  using aol::BiOp<aol::Vector<RealType> >::apply;
};

}
#endif
