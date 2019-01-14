#ifndef __GRADIENTFLOW_H
#define __GRADIENTFLOW_H

#include <ctrlCCatcher.h>
#include <op.h>
#include <gridBase.h>
#include <gnuplotter.h>
#include <progressBar.h>
#include <FEOpInterface.h>

#ifdef _OPENMP
#include "omp.h"
#endif

namespace aol {

/**
 * Class to help validating derivatives by plotting a energy and the
 * linear approximation of this energy using the derivative.
 *
 * \author Berkels, Droske
 */
template <typename RealType, typename VecType, typename DerivativeType = VecType>
class DescentDirValidator {
protected:
  aol::Vector<RealType> _positions;
  const aol::Op<DerivativeType , DerivativeType> &_massOp;
  const aol::Op<VecType, aol::Scalar<RealType> > &_energyOp;
  bool _quiet;
  const bool _useMassOp;

  RealType getScaledError( const aol::Vector<RealType> &Gateaux,
                           const aol::Vector<RealType> &Energy ) const {
    aol::Vector<RealType> errorVector (_positions.size());
    errorVector = Energy;
    errorVector -= Gateaux;
    errorVector.addToAll ( -errorVector.getMinValue() );
    const RealType energyRange = Energy.getMaxValue() - Energy.getMinValue();
    if ( energyRange != 0. )
      errorVector /= energyRange;
    return (errorVector.norm() / errorVector.size());
  }

  void writePlotFiles ( const aol::Vector<RealType> &Positions,
                        const aol::Vector<RealType> &Gateaux,
                        const aol::Vector<RealType> &Energy,
                        ofstream &Gateaux_out,
                        ofstream &En_out,
                        bool PlotOnlyGateaux ) const {
    En_out << setprecision ( 30 );
    Gateaux_out << setprecision ( 30 );

    for ( int i = 0; i < Positions.size(); i++ ) {
      if (!PlotOnlyGateaux)
        En_out << Positions[i] << " " << Energy[i] << endl;
      Gateaux_out << Positions[i]  << " " << Gateaux[i] << endl;
    }
  }

public:
  DescentDirValidator ( const RealType StepSize,
                        int NumSteps,
                        const aol::Op<DerivativeType, DerivativeType> &MassOp,
                        const aol::Op<VecType, aol::Scalar<RealType> > &EnergyOp,
                        const bool UseMassOp = true )
  : _positions ( NumSteps ),
    _massOp ( MassOp ),
    _energyOp ( EnergyOp ),
    _quiet ( true ),
    _useMassOp( UseMassOp ) {
    for (int i = 0; i < NumSteps; ++i)
      _positions[i] = StepSize * i;
  }

  DescentDirValidator ( const aol::Vector<RealType> & Positions,
                        const aol::Op<DerivativeType, DerivativeType> &MassOp,
                        const aol::Op<VecType, aol::Scalar<RealType> > &EnergyOp,
                        const bool UseMassOp = true )
  : _positions ( Positions ),
    _massOp ( MassOp ),
    _energyOp ( EnergyOp ),
    _quiet ( true ),
    _useMassOp( UseMassOp ) {}

  void setQuietMode ( const bool QuietMode ) {
    _quiet = QuietMode;
  }

  void validate ( const VecType &Var, const VecType &DescentDir, const DerivativeType &Gradient ) const {
    ofstream gateaux_out ( "gateaux.dat" );
    ofstream en_out ( "energy.dat" );
    validate ( Var, DescentDir, Gradient, gateaux_out, en_out );
    gateaux_out.close();
    en_out.close();
  }

  /** \brief compute evaluations and linearization of \a EnergyOp
   *         (given in the constructor) for a sequence of points \f$ x + \tau_i dx \f$
   *         (\f$ \tau_i \f$ also given in the constructor).
   *
   *  \param Var point at which the derivative is to be tested (\f$ x \f$)
   *  \param DescentDir direction in which the derivative is to be tested (\f$ dx \f$)
   *  \param Gradient suggested full gradient \f$ DE \f$ of \a EnergyOp at \f$ x \f$
   *  \param Positions will contain step sizes \f$ \tau_i \f$ after function call
   *  \param Gateaux will contain linearization \f$ E(x) + \tau_i DE \cdot dx \f$ after function call
   *  \param Energy will contain evaluations \f$ E(x + \tau_i dx) \f$ of \a EnergyOp \f$ F \f$
   *         after function call.
   *  \param negRate ratio of negative step sizes. Works only for linearly spaced step size
   *         vector! Example: If you pass 0.333 as \a negRate, one third of the positions
   *         (step sizes) will be negative.
   *
   *  \pre \a Positions, \a Gateaux and \a Energy have at least size \a NumSteps of
   *       first \c DescentDirValidator constructor or \c Positions.size() of
   *       second \c DescentDirValidator constructor, depending on which one was called.
   */
  void validate ( const VecType &Var,
                  const VecType &DescentDir,
                  const DerivativeType &Gradient,
                  aol::Vector<RealType> &Positions,
                  aol::Vector<RealType> &Gateaux,
                  aol::Vector<RealType> &Energy,
                  RealType negRate = aol::ZOTrait<RealType>::zero ) const {
    if ( Gradient.checkForNANsAndINFs() )
      cerr << "DescentDirValidator::validate Gradient contains NAN or INF!\n";
    DerivativeType gr ( Gradient );
    if( _useMassOp )
      _massOp.apply ( Gradient, gr );

    Positions = _positions;
    const int negRateIndexOffset = static_cast<int> ( -negRate*Positions.size() );
    if (negRateIndexOffset)
      Positions.addToAll ( -negRate * Positions[Positions.size() - 1] );

    aol::Scalar<RealType> energy;
    RealType gateauxDerivative =  gr * DescentDir;
    _energyOp.apply ( Var, energy );
    RealType initialEnergy = energy[0];
    VecType pos ( Var );

    for ( int i = 0; i < Positions.size(); i++ ) {
      pos = Var;
      pos.addMultiple ( DescentDir, Positions[i] );
      _energyOp.apply ( pos, energy );
      Energy[i] = energy[0];
      Gateaux[i] = (initialEnergy + Positions[i] * gateauxDerivative);
      if ( !_quiet ) {
        cerr << "step size = " << Positions[i] << "\tenergy = " << Energy[i] << "\tderivative = " << Gateaux[i] << "\tdiff. = " << Energy[i] - Gateaux[i] << endl;
      }
    }
  }

  // negRate is the rate of the test in negative DescentDir-direction (between 0 and 1)
  void validate ( const VecType &Var,
                  const VecType &DescentDir,
                  const DerivativeType &Gradient,
                  ofstream &Gateaux_out,
                  ofstream &En_out,
                  RealType negRate = aol::ZOTrait<RealType>::zero ) const {
    int numSteps = _positions.size();
    aol::Vector<RealType> positions( numSteps );
    aol::Vector<RealType> gateaux( numSteps );
    aol::Vector<RealType> energy( numSteps );
    validate ( Var, DescentDir, Gradient, positions, gateaux, energy, negRate );

    writePlotFiles ( positions, gateaux, energy, Gateaux_out, En_out );
  }

  void validateAndPlotToPNG ( const VecType &Var,
                              const VecType &DescentDir,
                              const DerivativeType &Gradient,
                              const char *BaseFileName,
                              RealType negRate = aol::ZOTrait<RealType>::zero,
                              const bool SkipPlottingWithSmallError = false,
                              const RealType SkippingThreshold = 1e-8,
                              bool PlotDifference = false ) const {
    int numSteps = _positions.size();
    aol::Vector<RealType> positions( numSteps );
    aol::Vector<RealType> gateaux( numSteps );
    aol::Vector<RealType> energy( numSteps );
    validate ( Var, DescentDir, Gradient, positions, gateaux, energy, negRate );

    if ( !SkipPlottingWithSmallError || getScaledError( gateaux, energy ) > SkippingThreshold ){
      if (PlotDifference)
        gateaux -= energy;

      char gateauxOutTempFileName[1024];
      sprintf ( gateauxOutTempFileName, "gateaux.datXXXXXX" );
      char enOutTempFileName[1024];
      sprintf ( enOutTempFileName, "energy.datXXXXXX" );
      ofstream gateaux_out;
      generateTemporaryFile ( gateauxOutTempFileName, gateaux_out );
      ofstream en_out;
      generateTemporaryFile ( enOutTempFileName, en_out );
      writePlotFiles ( positions, gateaux, energy, gateaux_out, en_out, PlotDifference );
      gateaux_out.close();
      en_out.close();
      if (PlotDifference) {
        aol::Plotter<RealType> plotter ( gateauxOutTempFileName );
        plotter.set_outfile_base_name ( BaseFileName );
        plotter.genPNG();
      }
      else {
        aol::Plotter<RealType> plotter ( gateauxOutTempFileName, enOutTempFileName, "Gateaux derivative", "Energy" );
        plotter.set_outfile_base_name ( BaseFileName );
        plotter.genPNG();
      }
      remove ( gateauxOutTempFileName );
      remove ( enOutTempFileName );
    }
  }

  /** \brief compute differences between directional derivative and difference quotient
   *         of \a EnergyOp (given in the constructor) for a sequence of points \f$ x + \tau_i dx \f$
   *         (\f$ \tau_i \f$ also given in the constructor).
   *
   *  \param Var point at which the derivative is to be tested (\f$ x \f$)
   *  \param DescentDir direction in which the derivative is to be tested (\f$ dx \f$)
   *  \param Gradient suggested full gradient \f$ DE \f$ of \a EnergyOp at \f$ x \f$
   *  \param Positions will contain step sizes \f$ \tau_i \f$ after function call
   *  \param DerivativeMinusDiffQuot will contain differences between directional
   *         derivative and forward difference quotient \f[ \frac{E(x + \tau_i dx)
   *         - E(x)}{\tau_i} - DE \cdot dx \f] after function call
   *
   *  uses validate().
   *
   *  \author von Deylen
   */
  void directionalDerivMinusDiffQuot ( const VecType &Var,
                                       const VecType &DescentDir,
                                       const DerivativeType &Gradient,
                                       aol::Vector<RealType> &Positions,
                                       aol::Vector<RealType> &DerivativeMinusDiffQuot ) {
    aol::Vector<RealType> & linearization = DerivativeMinusDiffQuot;
    aol::Vector<RealType> energyValues ( Positions.size() );

    validate ( Var, DescentDir, Gradient, Positions, linearization, energyValues );

    aol::Scalar<RealType> energy;
    _energyOp.apply ( Var, energy );
    RealType energyAtVar = energy[0];

    for (int i = 0; i < Positions.size(); ++i) {
      // linearization[i] contains E(x) + tau_i DE * dx
      // energyValues[i]  contains E(x + tau_i dx)
      // we compute dir. deriv.   DE * dx
      // and diff. quot.          ( E(x + tau_i dx) - E(x) ) / tau_i
      // from this:
      RealType tau = Positions[i];
      RealType directionalDeriv = (linearization[i] - energyAtVar) / tau;
      RealType diffQuot = (energyValues[i] - energyAtVar) / tau;
      DerivativeMinusDiffQuot[i] = directionalDeriv - diffQuot;
    }
  }

};

/**
 * gets a vector-valued Op \f$ F \f$ in the constructor, returns
 * \f$ 1/2 ||F(x)||^2 \f$ when being applied to \f$ x \f$.
 *
 * \author Berkels
 */
template <typename VecType>
class ComposeWithNormSqrOp : public Op<VecType, Scalar<typename VecType::DataType> > {
protected:
  const aol::Op<VecType, VecType> &_vectorValuedOp;
public:
  ComposeWithNormSqrOp ( const aol::Op<VecType, VecType> &VectorValuedOp )
    : _vectorValuedOp ( VectorValuedOp ) {}

  virtual ~ComposeWithNormSqrOp () {}

  virtual void applyAdd ( const VecType &Arg, Scalar<typename VecType::DataType> &Dest ) const {
    VecType temp ( Arg, aol::STRUCT_COPY );
    _vectorValuedOp.apply ( Arg, temp );
    Dest[0] += 0.5*( temp.normSqr() );
  }
};

/**
 * gets a vector-valued Op \f$ F \f$ and its derivative \f$ DF \f$
 * in the constructor, returns \f$ \partial_x 1/2 ||F(x)||^2
 * = DF(x) F(x) \f$ when being applied to \f$ x \f$.
 *
 * \author Berkels
 */
template <typename VecType, typename SecondDerivativeType>
class ComposeWithNormSqrDerivOp : public Op<VecType> {
protected:
  const aol::Op<VecType, VecType> &_vectorValuedOp;
  const aol::Op<VecType, SecondDerivativeType> &_matrixValuedOp;
  mutable SecondDerivativeType *_pMat;
public:
  ComposeWithNormSqrDerivOp ( SecondDerivativeType *PMat,
                         const aol::Op<VecType, VecType> &VectorValuedOp,
                         const aol::Op<VecType, SecondDerivativeType> &MatrixValuedOp )
    : _vectorValuedOp ( VectorValuedOp ),
      _matrixValuedOp ( MatrixValuedOp ),
      _pMat ( PMat ) {}
  virtual ~ComposeWithNormSqrDerivOp () {}
  virtual void applyAdd ( const VecType &Arg, VecType &Dest ) const {
    _matrixValuedOp.apply ( Arg, *_pMat );
    VecType temp ( Arg, aol::STRUCT_COPY );
    _vectorValuedOp( Arg, temp );
    _pMat->applyAdd ( temp, Dest );
  }
};

template <typename VecType, typename DerivativeType = VecType>
class FirstDerivativeTester : public Op<VecType, Scalar<typename VecType::DataType> > {
protected:
  const aol::Op<VecType, DerivativeType> &_vectorValuedOp;
  const VecType _direction;
public:
  FirstDerivativeTester ( const aol::Op<VecType, DerivativeType> &VectorValuedOp,
                          const VecType &Direction )
      : _vectorValuedOp ( VectorValuedOp ),
      _direction ( Direction ) {}
  virtual ~FirstDerivativeTester () {}
  virtual void apply ( const VecType &Arg, Scalar<typename VecType::DataType> &Dest ) const {
    DerivativeType temp ( Arg );
    _vectorValuedOp.apply ( Arg, temp );
    Dest[0] = ( temp * _direction );
  }
  virtual void applyAdd ( const VecType &Arg, Scalar<typename VecType::DataType> &Dest ) const {
    DerivativeType temp ( Arg );
    _vectorValuedOp.apply ( Arg, temp );
    Dest[0] += ( temp * _direction );
  }
};

template <typename VecType, typename SecondDerivativeType>
class SecondDerivativeTester : public Op<VecType> {
protected:
  const aol::Op<VecType, SecondDerivativeType> &_matrixValuedOp;
  const VecType _direction;
  mutable SecondDerivativeType *_pMat;
public:
  SecondDerivativeTester ( SecondDerivativeType *PMat,
                           const aol::Op<VecType, SecondDerivativeType> &MatrixValuedOp,
                           const VecType &Direction )
      : _matrixValuedOp ( MatrixValuedOp ),
      _direction ( Direction ),
      _pMat ( PMat ) {}
  virtual ~SecondDerivativeTester() {
  }
  virtual void apply ( const VecType &Arg, VecType &Dest ) const {
    _matrixValuedOp.apply ( Arg, *_pMat );
    _pMat->apply ( _direction, Dest );
  }
  virtual void applyAdd ( const VecType &Arg, VecType &Dest ) const {
    _matrixValuedOp.apply ( Arg, *_pMat );
    _pMat->applyAdd ( _direction, Dest );
  }
};

/**
 * \brief Method to test the first derivative of an energy.
 * This method plots the graph of the energy function E, evaluating it
 * at the points Position+j*StepSize*Direction for j all integers
 * running from -NegRate*NumSteps to (1-NegRate)*NumSteps.
 * Additionally, it plots E(Position)+j*StepSize*(DE*Direction),
 * where the bracketed term denotes a scalar product (i.e. DE is to be
 * interpreted as the gradient of E wrt the scalar product of VecType).
 * The plot is saved under the name provided.
 */
template<typename VecType, typename RealType>
void testFirstDerivative ( const VecType &Position,
                           const aol::Op<VecType, aol::Scalar<RealType> > &E,
                           const aol::Op<VecType> &DE,
                           const VecType &Direction,
                           const char *BasePlotFileName,
                           const RealType StepSize = 1e-2,
                           const int NumSteps = 120,
                           const RealType NegRate = 0 ) {
  VecType derivative ( Direction );
  DE.apply ( Position, derivative );

  aol::IdentityOp<VecType> identity;

  aol::DescentDirValidator<RealType, VecType> validator ( StepSize, NumSteps, identity, E, false );
  validator.validateAndPlotToPNG ( Position, Direction, derivative, BasePlotFileName, NegRate );
}

template<typename ConfiguratorType>
void testFirstDerivativeVector ( const typename ConfiguratorType::InitType &Grid,
                                 const aol::Vector<typename ConfiguratorType::RealType> &Position,
                                 const aol::Op<aol::Vector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > &E,
                                 const aol::Op<aol::Vector<typename ConfiguratorType::RealType> > &DE,
                                 const aol::Vector<typename ConfiguratorType::RealType> &Direction,
                                 const char *BasePlotFileName,
                                 typename ConfiguratorType::RealType StepSize = 0.0001 ) {
  typedef typename ConfiguratorType::RealType RealType;
  aol::Vector<RealType> gr ( Grid.getNumberOfNodes() );        // gradient
  aol::Vector<RealType> derivative ( Grid.getNumberOfNodes() );
  DE.apply ( Position, derivative );

  aol::LumpedMassOp<ConfiguratorType> lumpedMass ( Grid, aol::DO_NOT_INVERT );
  aol::LumpedMassOp<ConfiguratorType> lumpedMassInv ( Grid, aol::INVERT );

  lumpedMassInv.apply ( derivative, gr );

  aol::DescentDirValidator<RealType, aol::Vector<RealType> > validator ( StepSize * Grid.H() * Grid.H(), 120, lumpedMass, E );
  validator.validateAndPlotToPNG ( Position, Direction, gr, BasePlotFileName );
}

template<typename ConfiguratorType>
void testFirstDerivativeVector ( const typename ConfiguratorType::InitType &Grid,
                                 const aol::Vector<typename ConfiguratorType::RealType> &Position,
                                 const aol::Op<aol::Vector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > &E,
                                 const aol::Vector<typename ConfiguratorType::RealType> &DE,
                                 const char *BasePlotFileName,
                                 typename ConfiguratorType::RealType StepSize = 0.0001 ) {
  typedef typename ConfiguratorType::RealType RealType;
  aol::Vector<RealType> gr ( Grid.getNumberOfNodes() );      // gradient
  aol::Vector<RealType> dd ( DE );                  // descent direction

  aol::LumpedMassOp<ConfiguratorType> lumpedMass ( Grid, aol::DO_NOT_INVERT );
  aol::LumpedMassOp<ConfiguratorType> lumpedMassInv ( Grid, aol::INVERT );

  lumpedMassInv.apply ( dd, gr );                   // descent direction is the negative gradient
  dd = gr;
  dd *= -1.;

  aol::DescentDirValidator<RealType, aol::Vector<RealType> > validator ( StepSize * Grid.H() * Grid.H(), 120, lumpedMass, E );
  validator.validateAndPlotToPNG ( Position, dd, gr, BasePlotFileName );
}

template<typename ConfiguratorType>
void testFirstDerivativeVector ( const typename ConfiguratorType::InitType &Grid,
                                 const aol::Vector<typename ConfiguratorType::RealType> &Position,
                                 const aol::Op<aol::Vector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > &E,
                                 const aol::Op<aol::Vector<typename ConfiguratorType::RealType> > &DE,
                                 const char *BasePlotFileName,
                                 typename ConfiguratorType::RealType StepSize = 0.0001 ) {
  typedef typename ConfiguratorType::RealType RealType;
  aol::Vector<RealType> dd ( Position, aol::STRUCT_COPY );

  DE.apply ( Position, dd );
  testFirstDerivativeVector<ConfiguratorType>( Grid, Position, E, dd, BasePlotFileName, StepSize );
}

template<typename ConfiguratorType>
void testFirstDerivativeVectorAllDirections ( const typename ConfiguratorType::InitType &Grid,
                                              const aol::Vector<typename ConfiguratorType::RealType> &Position,
                                              const aol::Op < aol::Vector<typename ConfiguratorType::RealType>,
                                                              aol::Scalar<typename ConfiguratorType::RealType> > &E,
                                              const aol::Op<aol::Vector<typename ConfiguratorType::RealType> > &DE,
                                              const char *BasePlotFileName,
                                              typename ConfiguratorType::RealType StepSize = 0.0001 ) {
  typedef typename ConfiguratorType::RealType RealType;

  aol::Vector<RealType> dd ( Position, aol::STRUCT_COPY );   // test direction
  aol::Vector<RealType> gr ( Position, aol::STRUCT_COPY );

  aol::LumpedMassOp<ConfiguratorType> lumpedMass ( Grid, aol::DO_NOT_INVERT );

  DE.apply ( Position, gr );
  aol::DescentDirValidator<RealType, aol::Vector<RealType> > validator ( StepSize * Grid.H() * Grid.H(), 120, lumpedMass, E, false );
  char filename[1024];
  for ( int i = 0; i < Position.size(); i++ ) {
    sprintf ( filename, "%s_%03d", BasePlotFileName, i );
    dd.setAll ( 0. );
    dd[i] = 1.;
    validator.validateAndPlotToPNG ( Position, dd, gr, filename );
  }
}

template<typename ConfiguratorType>
void testFirstDerivativeMultiVector ( const typename ConfiguratorType::InitType &Grid,
                                      const aol::MultiVector<typename ConfiguratorType::RealType> &Position,
                                      const aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > &E,
                                      const aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > &DE,
                                      const aol::MultiVector<typename ConfiguratorType::RealType> &Direction,
                                      const char *BasePlotFileName,
                                      typename ConfiguratorType::RealType StepSize = 0.0001 ) {
  typedef typename ConfiguratorType::RealType RealType;

  aol::MultiVector<RealType> derivative ( Direction, aol::STRUCT_COPY );
  DE.apply ( Position, derivative );

  aol::LumpedMassOp<ConfiguratorType> lumpedMass ( Grid, aol::DO_NOT_INVERT );
  aol::DiagonalBlockOp<RealType> blockLumpedMass ( lumpedMass );
  aol::LumpedMassOp<ConfiguratorType> lumpedMassInv ( Grid, aol::INVERT );
  aol::DiagonalBlockOp<RealType> blockLumpedMassInv ( lumpedMassInv );

  aol::MultiVector<RealType> gradient ( Position, aol::STRUCT_COPY );
  blockLumpedMassInv.apply ( derivative, gradient );

  aol::DescentDirValidator<RealType, aol::MultiVector<RealType> > validator ( StepSize * Grid.H() * Grid.H(), 120, blockLumpedMass, E );
  validator.validateAndPlotToPNG ( Position, Direction, gradient, BasePlotFileName );
}

template<typename ConfiguratorType>
void testFirstDerivativeMultiVector ( const typename ConfiguratorType::InitType &Grid,
                                      const aol::MultiVector<typename ConfiguratorType::RealType> &Position,
                                      const aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > &E,
                                      const aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > &DE,
                                      const char *BasePlotFileName,
                                      typename ConfiguratorType::RealType StepSize = 0.0001 ) {
  typedef typename ConfiguratorType::RealType RealType;

  aol::MultiVector<RealType> dd ( Position, aol::STRUCT_COPY );   // test direction
  aol::MultiVector<RealType> gr ( Position, aol::STRUCT_COPY );

  aol::LumpedMassOp<ConfiguratorType> lumpedMass ( Grid, aol::DO_NOT_INVERT );
  aol::DiagonalBlockOp<RealType> blockLumpedMass ( lumpedMass );
  aol::LumpedMassOp<ConfiguratorType> lumpedMassInv ( Grid, aol::INVERT );
  aol::DiagonalBlockOp<RealType> blockLumpedMassInv ( lumpedMassInv );

  DE.apply ( Position, dd );
  blockLumpedMassInv.apply ( dd, gr );
  dd = gr;
  dd *= -1.;
  aol::DescentDirValidator<RealType, aol::MultiVector<RealType> > validator ( StepSize * Grid.H() * Grid.H(), 120, blockLumpedMass, E );
  validator.validateAndPlotToPNG ( Position, dd, gr, BasePlotFileName );
}

template<typename ConfiguratorType>
void testFirstDerivativeMultiVectorAllDirections ( const typename ConfiguratorType::InitType &Grid,
                                                   const aol::MultiVector<typename ConfiguratorType::RealType> &Position,
                                                   const aol::Op < aol::MultiVector<typename ConfiguratorType::RealType>,
                                                   aol::Scalar<typename ConfiguratorType::RealType> > &E,
                                                   const aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > &DE,
                                                   const char *BasePlotFileName,
                                                   typename ConfiguratorType::RealType StepSize = 0.0001 ) {
  typedef typename ConfiguratorType::RealType RealType;

  aol::MultiVector<RealType> dd ( Position, aol::STRUCT_COPY );   // test direction
  aol::MultiVector<RealType> gr ( Position, aol::STRUCT_COPY );

  aol::LumpedMassOp<ConfiguratorType> lumpedMass ( Grid, aol::DO_NOT_INVERT );
  aol::DiagonalBlockOp<RealType> blockLumpedMass ( lumpedMass );

  DE.apply ( Position, gr );
  aol::DescentDirValidator<RealType, aol::MultiVector<RealType> > validator ( StepSize * Grid.H() * Grid.H(), 120, blockLumpedMass, E, false );
  char filename[1024];
  for ( int i = 0; i < Position.numComponents(); i++ ) {
    for ( int j = 0; j < Position[i].size(); j++ ) {
      sprintf ( filename, "%s_%d_%03d", BasePlotFileName, i, j );
      dd.setAll ( 0. );
      dd[i][j] = 1.;
      validator.validateAndPlotToPNG ( Position, dd, gr, filename );
    }
  }
}

template<typename ConfiguratorType>
void testFirstDerivativeMultiVectorSingleComponentAllDirections
 ( const typename ConfiguratorType::InitType &Grid,
   const aol::MultiVector<typename ConfiguratorType::RealType> &Position,
   const aol::Op < aol::MultiVector<typename ConfiguratorType::RealType>,
   aol::Scalar<typename ConfiguratorType::RealType> > &E,
   const aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > &DE,
   const int Component,
   const char *BasePlotFileName,
   typename ConfiguratorType::RealType StepSize = 0.0001 )
{
  typedef typename ConfiguratorType::RealType RealType;

  aol::MultiVector<RealType> dd ( Position, aol::STRUCT_COPY );   // test direction
  aol::MultiVector<RealType> gr ( Position, aol::STRUCT_COPY );

  aol::LumpedMassOp<ConfiguratorType> lumpedMass ( Grid, aol::DO_NOT_INVERT );
  aol::DiagonalBlockOp<RealType> blockLumpedMass ( lumpedMass );

  DE.apply ( Position, gr );
  aol::DescentDirValidator<RealType, aol::MultiVector<RealType> > validator ( StepSize * Grid.H() * Grid.H(), 120, blockLumpedMass, E, false );
  char filename[1024];
  for ( int j = 0; j < Position[Component].size(); j++ ) {
    sprintf ( filename, "%s_%03d", BasePlotFileName, j );
    dd.setAll ( 0. );
    dd[Component][j] = 1.;
    validator.validateAndPlotToPNG ( Position, dd, gr, filename );
  }
}

template<typename ConfiguratorType, typename SecondDerivativeType>
void testSecondDerivativeVector ( const typename ConfiguratorType::InitType &Grid,
                                  const aol::Vector<typename ConfiguratorType::RealType> &Position,
                                  const aol::Vector<typename ConfiguratorType::RealType> &Direction,
                                  const aol::Op<aol::Vector<typename ConfiguratorType::RealType> > &F,
                                  const aol::Op<aol::Vector<typename ConfiguratorType::RealType>, SecondDerivativeType > &DF,
                                  const char *BasePlotFileName ) {
  typedef typename ConfiguratorType::RealType RealType;
  aol::FirstDerivativeTester<aol::Vector<RealType> > b ( F, Direction );
  aol::SecondDerivativeTester<aol::Vector<RealType>, SecondDerivativeType > c ( new SecondDerivativeType ( Grid ), DF, Direction );
  aol::Vector<RealType> dd ( Grid.getNumberOfNodes() );
  aol::Vector<RealType> gr ( Grid.getNumberOfNodes() );

  aol::LumpedMassOp<ConfiguratorType> lumpedMass ( Grid, aol::DO_NOT_INVERT );
  aol::LumpedMassOp<ConfiguratorType> lumpedMassInv ( Grid, aol::INVERT );

  c.apply ( Position, dd );

  lumpedMassInv.apply ( dd, gr );
  dd = gr;
  dd *= -1.;

  aol::DescentDirValidator<RealType, aol::Vector<RealType> > validator ( 0.0001 * Grid.H() * Grid.H(), 120, lumpedMass, b );
  validator.validateAndPlotToPNG ( Position, dd, gr, BasePlotFileName );
}

template<typename ConfiguratorType, typename SecondDerivativeType>
void testSecondDerivativeMultiVector ( const typename ConfiguratorType::InitType &Grid,
                                       const aol::MultiVector<typename ConfiguratorType::RealType> &Position,
                                       const aol::MultiVector<typename ConfiguratorType::RealType> &Direction,
                                       const aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > &F,
                                       const aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, SecondDerivativeType > &DF,
                                       const char *BasePlotFileName,
                                       SecondDerivativeType* PMatDF = new SecondDerivativeType ( ConfiguratorType::Dim, ConfiguratorType::Dim ),
                                       typename ConfiguratorType::RealType StepSize = 0.0001,
                                       const bool DeleteMatrixPointer = true ) {
  typedef typename ConfiguratorType::RealType RealType;
  aol::FirstDerivativeTester<aol::MultiVector<RealType> > b ( F, Direction );
  aol::SecondDerivativeTester<aol::MultiVector<RealType>, SecondDerivativeType > c ( PMatDF, DF, Direction );
  aol::MultiVector<RealType> dd ( Position, aol::STRUCT_COPY );
  aol::MultiVector<RealType> gr ( Position, aol::STRUCT_COPY );

  aol::LumpedMassOp<ConfiguratorType> lumpedMass ( Grid, aol::DO_NOT_INVERT );
  aol::DiagonalBlockOp<RealType> blockLumpedMass ( lumpedMass );
  aol::LumpedMassOp<ConfiguratorType> lumpedMassInv ( Grid, aol::INVERT );
  aol::DiagonalBlockOp<RealType> blockLumpedMassInv ( lumpedMassInv );

  c.apply ( Position, dd );

  blockLumpedMassInv.apply ( dd, gr );
  dd = gr;
  dd *= -1.;

  aol::DescentDirValidator<RealType, aol::MultiVector<RealType> > validator ( StepSize * Grid.H() * Grid.H(), 120, blockLumpedMass, b );
  validator.validateAndPlotToPNG ( Position, dd, gr, BasePlotFileName );

  if ( DeleteMatrixPointer )
    delete PMatDF;
}

template<typename ConfiguratorType, typename SecondDerivativeType>
void testSecondDerivativeMultiVectorAllDirections ( const typename ConfiguratorType::InitType &Grid,
                                                    const aol::MultiVector<typename ConfiguratorType::RealType> &Position,
                                                    const aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > &F,
                                                    const aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, SecondDerivativeType > &DF,
                                                    const char *BasePlotFileName,
                                                    SecondDerivativeType* PMatDF = new SecondDerivativeType ( ConfiguratorType::Dim, ConfiguratorType::Dim ),
                                                    typename ConfiguratorType::RealType StepSize = 0.0001,
                                                    const bool DeleteMatrixPointer = true ) {
  typedef typename ConfiguratorType::RealType RealType;
  aol::MultiVector<RealType> mdirection( Position, aol::STRUCT_COPY );

  for( int i = 0; i < Position.numComponents(); i++){
    mdirection.setZero();
    for( int j = 0; j < Grid.getNumberOfNodes(); j++){
      mdirection[i].setZero();
      mdirection[i][j] = 1.;
      char fn[1024];
      sprintf( fn, "%s_%d_%03d", BasePlotFileName, i, j );
      aol::testSecondDerivativeMultiVector<ConfiguratorType, SecondDerivativeType>
        ( Grid, Position, mdirection, F, DF, fn, PMatDF, StepSize, false );
    }
  }
  if ( DeleteMatrixPointer )
    delete PMatDF;
}

/**
 * \brief Class that collects the code shared by FirstDerivativeValidator and SecondDerivativeValidator.
 *
 * \author Berkels
 */
template <typename RealType, typename VecType >
class DerivativeValidatorBase {
public:
  enum StepSizeFunctionType { LINEAR, EXPONENTIAL };

  DerivativeValidatorBase ( const RealType H,
                            StepSizeFunctionType stepsizeFctType,
                            const RealType StepSize,
                            const RealType StepSizeMin,
                            const RealType StepSizeMax )
    : _h(H),
      _stepSizeFctType ( stepsizeFctType ),
      _stepSize ( StepSize ),
      _stepSizeMin ( StepSizeMin ),
      _stepSizeMax ( StepSizeMax ),
      _skipPlottingWithSmallError ( false ),
      _skippingThreshold ( 1e-8 ),
      _maskPtr ( NULL ),
      _catchCtrlC ( false ),
      _negRate ( 0. ),
      _numberOfSamples ( 120 ) {
    this->buildPositionsVector ( _positions );
  }

  void setSkipPlottingWithSmallError ( const bool SkipPlottingWithSmallError ) {
    _skipPlottingWithSmallError = SkipPlottingWithSmallError;
  }

  void setSkippingThreshold ( const RealType SkippingThreshold ) {
    _skippingThreshold = SkippingThreshold;
  }

  //! All derivative tests that loop over multiple directions
  //! can be interrupted by Ctrl-C if you have set \a catchCtrlC to
  //! true. The normal Ctrl-C handler is replaced by the QuocMesh handler
  //! during the iteration. When you press Ctrl-C, the test for the
  //! current derivative direction is finished, and afterwards a message
  //! is printed you really want to interrupt. There you can choose
  //! between terminating the program by pressing Ctrl-C again,
  //! cancelling the current derivative test and continuing the
  //! program or even resuming the derivative test.
  void setCatchCtrlC ( bool catchCtrlC ) {
    _catchCtrlC = catchCtrlC;
  }

  //! if you do not want all entries to be tested, you may set
  //! a mask before calling testAllDirections() or similar fct's.
  void setMask ( const BitVector & mask ) {
    _maskPtr = &mask;
  }

  void unsetMask () {
    _maskPtr = NULL;
  }

  const BitVector & getMask () const {
    return *_maskPtr;
  }

  //! set parameter negRate which is passed to validateAndPlotToPNG
  void setNegRate ( const RealType negRate ) {
    _negRate = negRate;
  }

  //! set number of samples in LINEAR case
  void setNumberOfSamples ( const int numberOfSamples ) {
    _numberOfSamples = numberOfSamples;
    this->buildPositionsVector ( _positions );
  }


protected:
  const RealType _h;
  const StepSizeFunctionType _stepSizeFctType;
  const RealType _stepSize;
  const RealType _stepSizeMin, _stepSizeMax;
  bool _resizePositionVector;
  bool _skipPlottingWithSmallError;
  RealType _skippingThreshold;
  aol::Vector<RealType> _positions;
  const BitVector * _maskPtr; //< only vector entries marked "true" will be tested.
  bool _catchCtrlC;
  RealType _negRate;
  int _numberOfSamples;
  mutable sigfunc _previousCtrlCHandler;


  /** \brief fill passed vector with values at which the difference quotients
   *  will be evaluated.
   *
   *  If \a StepSizeFctType in the constructor was passed as LINEAR (default),
   *  a vector of size 120 with entries \a i * \a StepSize * \a h
   *  will be constructed.
   *
   *  In case of \a StepSizeFctType = EXPONENTIAL, entries of \a Positions
   *  will be \f$ \mbox{StepSizeMin} \cdot \mbox{StepSize}^i \f$ for
   *  \f$ i = 1,\dots,\log(\mbox{StepSizeMin} / \mbox{StepSizeMax}) /
   *  \log \mbox{StepSize} \f$. If the same value of \a StepSizeMin and
   *  \a StepSizeMax was passed in the constructor, a vector with one
   *  entry is built.
   *
   *  \return length of \a Positions after filling.
   *
   *  \author von Deylen
   */
  int buildPositionsVector ( aol::Vector<RealType> & Positions ) const {
    switch ( _stepSizeFctType ) {
      case LINEAR:
        Positions.resize ( _numberOfSamples );
        for ( int i = 0; i < _numberOfSamples; ++i )
          Positions[i] = i * _stepSize * Sqr ( _h );
        return _numberOfSamples;
      case EXPONENTIAL: {
        int n;
        if ( log ( 1. / _stepSize ) < 1E-15 )
          n = 1;
        else
          n = static_cast<int> ( log ( _stepSizeMax / _stepSizeMin ) / log ( 1. / _stepSize ) );
        Positions.resize ( n );
        if ( n ) Positions[0] = _stepSizeMax;
        for ( int i = 1; i < n; ++i )
          Positions[i] = Positions[i - 1] * _stepSize;
        return n;
      }
      default:
        throw ( aol::UnimplementedCodeException ( "Unimplemented StepSizeFunctionType value.", __FILE__, __LINE__ ) );
    }
  }

  aol::Vec< VecType::Depth, int> setIthComponent ( const int I, const RealType Value, VecType &Vec ) const {
    return aol::Vec< VecType::Depth, int> ( Vec.setIthComponent ( I, Value ) );
  }

  string getSuffixString ( const aol::Vec< VecType::Depth, int>& idxVec ) const {
    char fn[1024];
    switch ( VecType::Depth ){
      case 1: sprintf( fn, "%03d", idxVec[0] ); break;
      case 2: sprintf( fn, "%d_%03d", idxVec[0], idxVec[1] ); break;
      case 3: sprintf( fn, "%d_%03d_%d", idxVec[0], idxVec[1], idxVec[2] );  break;
      default: throw aol::Exception ( "aol::DerivativeValidatorBase<>::getSuffixString: unvalid VecTyp.", __FILE__, __LINE__ );
    }
    string temp = fn;
    return temp;
  }

  //! checks if mask is stored. If so, checks I'th entry.
  bool wantsToCheckEntry ( const int I ) const {
    return !_maskPtr || (*_maskPtr)[I];
  }

  void setCtrlCHandler () const {
    if (_catchCtrlC)
      _previousCtrlCHandler = signal ( InterruptSignal, ctrlCHandler );
  }

  void unsetCtrlCHandler () const {
    if (_catchCtrlC)
      signal ( InterruptSignal, _previousCtrlCHandler );
  }

  bool wantsInterrupt() const {
    if (!_catchCtrlC || !getCtrlCState())
      return false;

    signal ( InterruptSignal, DefaultHandler );
    cerr << endl << endl << aol::color::error
         << "Do you really want to interrupt the derivative test "
            "(y/n, Ctrl-C to kill program)? ";
    cerr << aol::color::reset;
    string yes_no;
    cin >> yes_no;
    signal ( InterruptSignal, ctrlCHandler );

    if (yes_no == "y" || yes_no == "yes")
      return true;
    else {
      resetCtrlCState();
      return false;
    }
  }

};

/**
 * \brief Class to validate first derivatives or derivatives of vector valued functions.
 *
 * \note Meant to supersede all testFirstDerivative* functions.
 *
 * \author Berkels
 */
template <typename VecType, typename DerivativeType = VecType, typename RealType = typename VecType::DataType>
class FirstDerivativeValidator : public DerivativeValidatorBase<RealType, VecType> {
  const aol::Op<VecType, aol::Scalar<RealType> > &_E;
  const aol::Op<VecType, DerivativeType> &_DE;
  typedef typename DerivativeValidatorBase<RealType, VecType>::StepSizeFunctionType StepSizeFunctionType;
public:
  FirstDerivativeValidator ( const aol::Op<VecType, aol::Scalar<RealType> > &E,
                             const aol::Op<VecType, DerivativeType> &DE,
                             const RealType H = 1.,
                             // Workaround for older GCC versions. 0 is supposed to be LINEAR.
                             StepSizeFunctionType stepSizeFctType = static_cast<StepSizeFunctionType> ( 0 ),
                             const RealType StepSize = 0.0001,
                             const RealType StepSizeMin = 1E-12,
                             const RealType StepSizeMax = 1E-2 )
    : DerivativeValidatorBase<RealType, VecType> ( H, stepSizeFctType, StepSize, StepSizeMin, StepSizeMax ),
      _E(E),
      _DE(DE) {
  }

  ~FirstDerivativeValidator() {}

  void testDirection ( const VecType &Position,
                       const VecType &Direction,
                       const DerivativeType &Gradient,
                       const char* BasePlotFileName,
                       bool PlotDifference = false ) const {
    aol::IdentityOp<DerivativeType> identity;
    aol::DescentDirValidator<RealType, VecType, DerivativeType> validator ( this->_positions, identity, _E );
    validator.validateAndPlotToPNG ( Position, Direction, Gradient, BasePlotFileName, this->_negRate, this->_skipPlottingWithSmallError, this->_skippingThreshold, PlotDifference );
  }

  void testDirection ( const VecType &Position,
                       const char* BasePlotFileName,
                       bool PlotDifference = false ) const {
    VecType descentDirection ( Position, aol::STRUCT_COPY );
    DerivativeType gradient ( Position, aol::STRUCT_COPY );

    _DE.apply ( Position, gradient );

    descentDirection = gradient;
    descentDirection *= -1.;

    testDirection ( Position, descentDirection, gradient, BasePlotFileName, PlotDifference );
  }

  void testDirection ( const VecType &Position,
                       const VecType &Direction,
                       const char* BasePlotFileName,
                       bool PlotDifference = false ) const {
    DerivativeType gradient ( Position, aol::STRUCT_COPY );
    _DE.apply ( Position, gradient );
    testDirection ( Position, Direction, gradient, BasePlotFileName, PlotDifference );
  }

  void testAllDirections ( const VecType &Position,
                           const char* BasePlotFileName,
                           bool PlotDifference = false ) const {
    VecType descentDirection ( Position, aol::STRUCT_COPY );
    DerivativeType gradient ( Position, aol::STRUCT_COPY );

    _DE.apply ( Position, gradient );

    aol::ProgressBar<> pb ( "Validating all first derivatives" );
    pb.start ( Position.getTotalSize() );
    this->setCtrlCHandler();
    for( int i = 0; i < Position.getTotalSize() && !this->wantsInterrupt(); i++, pb++ ){
      if (!this->wantsToCheckEntry(i)) continue;
      descentDirection.setZero();
      aol::Vec< VecType::Depth, int>  idxVec = this->setIthComponent ( i, 1., descentDirection );

      char fn[1024];
      sprintf( fn, "%s_%s", BasePlotFileName, this->getSuffixString( idxVec ).c_str() );
      testDirection( Position, descentDirection, gradient, fn, PlotDifference );
    }
    this->unsetCtrlCHandler();
    cerr << endl;
  }
  
  void testAllDirectionsParallel ( const VecType &Position,
                                   const char* BasePlotFileName,
                                   bool PlotDifference = false ) const {
    DerivativeType gradient ( Position, aol::STRUCT_COPY );
    _DE.apply ( Position, gradient );

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for( int i = 0; i < Position.getTotalSize(); i++ ){
      if (!this->wantsToCheckEntry(i)) continue;
      VecType descentDirection ( Position, aol::STRUCT_COPY );
      descentDirection.setZero();
      aol::Vec< VecType::Depth, int>  idxVec = this->setIthComponent ( i, 1., descentDirection );

      char fn[1024];
      sprintf( fn, "%s_%s", BasePlotFileName, this->getSuffixString( idxVec ).c_str() );
      testDirection( Position, descentDirection, gradient, fn, PlotDifference );
    }
  }
  
    void testRandomDirections ( const VecType &Position,
			      const int NumOfDirections,
                              const char* BasePlotFileName,
			      bool PlotDifference = false ) const {
    DerivativeType gradient ( Position, aol::STRUCT_COPY );
    _DE.apply ( Position, gradient );
    
    aol::RandomGenerator rGen;
    rGen.randomize();
    aol::ProgressBar<> pb ( "Validating first derivatives in randomly chosen directions " );
    pb.start ( NumOfDirections );
    this->setCtrlCHandler();
    for( int i = 0; i < NumOfDirections; i++ ){
     
      VecType descentDirection ( Position, aol::STRUCT_COPY );
      descentDirection.setZero();
      int directionNum = rGen.rInt( Position.getTotalSize() ); 
      if (!this->wantsToCheckEntry(directionNum)) continue;
      
      aol::Vec< VecType::Depth, int>  idxVec = this->setIthComponent ( directionNum, 1., descentDirection );

      char fn[1024];
      sprintf( fn, "%s_%s", BasePlotFileName, this->getSuffixString( idxVec ).c_str() );
      testDirection( Position, descentDirection, gradient, fn, PlotDifference );
      pb++;
    }
    this->unsetCtrlCHandler();
    cerr << endl;
  }
  

  RealType directionalDerivMinusDiffQuotient ( const VecType &Position,
                                               const VecType &Direction,
                                               const DerivativeType &Gradient,
                                               RealType StepSize = 1E-6) {
    aol::IdentityOp<DerivativeType> identity;
    aol::Vector<RealType> positions ( 1 ),
                          differences ( 1 );
    positions[0] = StepSize;
    aol::DescentDirValidator<RealType, VecType, DerivativeType> validator ( positions, identity, _E );
    validator.directionalDerivMinusDiffQuot ( Position, Direction, Gradient, positions, differences );
    return differences[0];
  }

  void derivMinusDiffQuotientsAllDirections ( const VecType & Position,
                                              VecType & derivMinusDiffQuots,
                                              RealType StepSize = 1E-6 ) {
    VecType descentDirection ( Position, aol::STRUCT_COPY );
    DerivativeType gradient ( Position, aol::STRUCT_COPY );

    _DE.apply ( Position, gradient );

    aol::ProgressBar<> pb ( "computing difference quotients for all components" );
    pb.start ( Position.getTotalSize() );
    this->setCtrlCHandler();
    for( int i = 0; i < Position.getTotalSize() && !this->wantsInterrupt(); i++, pb++ ){
      if (!this->wantsToCheckEntry(i)) continue;
      this->setIthComponent ( i, 1., descentDirection );
      this->setIthComponent ( i, directionalDerivMinusDiffQuotient ( Position, descentDirection, gradient, StepSize ), derivMinusDiffQuots );
      this->setIthComponent ( i, 0., descentDirection );
    }
    this->unsetCtrlCHandler();
    cerr << endl;
  }
};

/**
 * \brief Class to validate second derivatives or derivatives of vector valued functions.
 *        Input is a vector-valued fct \f$ F \f$, maybe the gradient of a scalar
 *        function \f$ \nabla E \f$, and \f$ DF = D^2 E \f$.
 *
 * \note  Meant to supersede all testSecondDerivative* functions.
 *
 * We can only test scalar-valued functions. Therefore, consider
 * an "outer direction" \f$ o \f$ and the objective function
 * \f[ b_o(x) := F(x) \cdot o. \f] In the test of a single
 * direction pair \f$ (o,i) \f$, we compare \f$ b_o(x + \tau i) \f$
 * to \f$ b_o(x) + \tau \partial_i b_o(x) = b_o(x) + \tau DF^T(x) i \f$.
 *
 * \warning This class assumes that domain and range of \f$ F \f$ have the same size!
 *
 * \author Berkels
 */
template <typename VecType, typename SecondDerivativeType>
class SecondDerivativeValidator : public DerivativeValidatorBase<typename VecType::DataType, VecType> {
  typedef typename VecType::DataType RealType;
  const aol::Op<VecType> &_F;
  const aol::Op<VecType, SecondDerivativeType > &_DF;
  mutable SecondDerivativeType* _pMatDF;
  const bool _deleteMatrixPointer;
  typedef typename DerivativeValidatorBase<RealType, VecType>::StepSizeFunctionType StepSizeFunctionType;
public:
  SecondDerivativeValidator ( const aol::Op<VecType> &F,
                              const aol::Op<VecType, SecondDerivativeType > &DF,
                              const RealType H,
                              SecondDerivativeType* PMatDF,
                              const bool DeleteMatrixPointer = true,
                              // Workaround for older GCC versions. 0 is supposed to be LINEAR.
                              StepSizeFunctionType stepSizeFctType = static_cast<StepSizeFunctionType> ( 0 ),
                              const RealType StepSize = 0.0001,
                              const RealType StepSizeMin = 0.,
                              const RealType StepSizeMax = 1. )
    : DerivativeValidatorBase<RealType, VecType> ( H, stepSizeFctType, StepSize, StepSizeMin, StepSizeMax ),
      _F(F),
      _DF(DF),
      _pMatDF(PMatDF),
      _deleteMatrixPointer(DeleteMatrixPointer) {
  }

  ~SecondDerivativeValidator() {
    if ( _deleteMatrixPointer )
      delete _pMatDF;
  }

  void testDirection ( const VecType &Position,
                       const VecType &OuterDirection,
                       const VecType &InnerDirection,
                       const VecType &GradientWRTOuterDirection,
                       const char* BasePlotFileName ) const {
    aol::FirstDerivativeTester<VecType> b ( _F, OuterDirection );
    aol::IdentityOp<VecType> identity;
    aol::DescentDirValidator<RealType, VecType> validator ( this->_positions, identity, b );
    validator.validateAndPlotToPNG ( Position, InnerDirection, GradientWRTOuterDirection, BasePlotFileName, this->_negRate, this->_skipPlottingWithSmallError, this->_skippingThreshold );
  }
  
  void testDirection ( const VecType &Position,
                       const VecType &OuterDirection,
                       const VecType &InnerDirection,
                       const char* BasePlotFileName,
                       const bool applyDF = true,
                       const bool secondDerivTypeIsSymmetric = true  ) const {
                           
    if( applyDF )
      _DF.apply ( Position, *_pMatDF );
    if( !secondDerivTypeIsSymmetric )
      _pMatDF->transpose();
    
    VecType gradientWRTOuterDirection( Position, aol::STRUCT_COPY );
    _pMatDF->apply ( OuterDirection, gradientWRTOuterDirection );    
    
    aol::FirstDerivativeTester<VecType> b ( _F, OuterDirection );
    aol::IdentityOp<VecType> identity;
    
    aol::DescentDirValidator<RealType, VecType> validator ( this->_positions, identity, b );
    validator.validateAndPlotToPNG ( Position, InnerDirection, gradientWRTOuterDirection, BasePlotFileName, this->_negRate, this->_skipPlottingWithSmallError, this->_skippingThreshold );
  }

  /** gets one direction \f$ v \f$, tests the pair
   *  \f$ (v,DF^T(x) v) \f$. The inner direction is chose
   *  as \f$ DF^T(x) v = \nabla b_v(x) \f$ as analogue
   *  to the first-derivative testing case where often
   *  only the direction \f$ \nabla E \f$ is tested.
   */
  void testDirection ( const VecType &Position,
                       const VecType &Direction,
                       const char* BasePlotFileName,
                       const bool applyDF = true,
                       const bool secondDerivTypeIsSymmetric = true ) const {
    if( applyDF )
      _DF.apply ( Position, *_pMatDF );
    if( !secondDerivTypeIsSymmetric )
      _pMatDF->transpose();
    
    VecType descentDirection ( Position, aol::STRUCT_COPY );
    VecType gradient ( Position, aol::STRUCT_COPY );

    _pMatDF->apply( Direction, descentDirection );

    gradient = descentDirection;
    descentDirection *= -1.;

    testDirection ( Position, Direction, descentDirection, gradient, BasePlotFileName );
  }

  /** performs a loop over all outer directions and proceeds
   *  with them as described in void testDirection ( const VecType &,
   *  const VecType &, const char * ) const.
   */
  void testAllOuterDirections ( const VecType &Position,
                                const char* BasePlotFileName,
                                const bool applyDF = true,
                                const bool secondDerivTypeIsSymmetric = true ) const {
    if( applyDF )
     _DF.apply ( Position, *_pMatDF );
    if( !secondDerivTypeIsSymmetric )
      _pMatDF->transpose();

    VecType outerDirection( Position, aol::STRUCT_COPY );
    VecType innerDirection( Position, aol::STRUCT_COPY );
    VecType gradientWRTOuterDirection( Position, aol::STRUCT_COPY );

    aol::ProgressBar<> pb ( "Validating outer second derivatives" );
    pb.start ( Position.getTotalSize() );
    this->setCtrlCHandler();
    for( int i = 0; i < Position.getTotalSize() && !this->wantsInterrupt(); i++, pb++ ){
      if (!this->wantsToCheckEntry(i)) continue;
      outerDirection.setZero();
      aol::Vec< VecType::Depth, int>  idxVec = this->setIthComponent ( i, 1., outerDirection );
      _pMatDF->apply ( outerDirection, gradientWRTOuterDirection );
      innerDirection = gradientWRTOuterDirection;
      innerDirection *= -1;

      char fn[1024];
      sprintf( fn, "%s_%s", BasePlotFileName, this->getSuffixString( idxVec ).c_str() );
      testDirection( Position, outerDirection, innerDirection, gradientWRTOuterDirection, fn );
    }
    this->unsetCtrlCHandler();
    cerr << endl;
  }

  void testAllDirections ( const VecType &Position,
                           const char* BasePlotFileName,
                           const bool applyDF = true,
                           const bool secondDerivTypeIsSymmetric = true ) const {
    if( applyDF )
      _DF.apply ( Position, *_pMatDF );
    if( !secondDerivTypeIsSymmetric )
      _pMatDF->transpose();

    VecType outerDirection( Position, aol::STRUCT_COPY );
    VecType innerDirection( Position, aol::STRUCT_COPY );
    VecType gradientWRTOuterDirection( Position, aol::STRUCT_COPY );

    aol::ProgressBar<> pb ( "Validating all second derivatives" );
    pb.start ( aol::Sqr(Position.getTotalSize()) );
    this->setCtrlCHandler();
    for( int outerDirectionNum = 0; outerDirectionNum < Position.getTotalSize(); outerDirectionNum++ ){
      outerDirection.setZero();
      aol::Vec< VecType::Depth, int>  idxVecOuter = this->setIthComponent ( outerDirectionNum, 1., outerDirection );

      _pMatDF->apply ( outerDirection, gradientWRTOuterDirection );

      for( int innerDirectionNum = 0; innerDirectionNum < Position.getTotalSize() && !this->wantsInterrupt(); innerDirectionNum++, pb++ ){
        if (!this->wantsToCheckEntry(outerDirectionNum) || !this->wantsToCheckEntry(innerDirectionNum)) continue;
        innerDirection.setZero();
        aol::Vec< VecType::Depth, int>  idxVecInner = this->setIthComponent ( innerDirectionNum, 1., innerDirection );
        char fn[1024];
        sprintf( fn, "%s_%s_%s", BasePlotFileName, this->getSuffixString( idxVecOuter ).c_str(), this->getSuffixString( idxVecInner ).c_str() );
        testDirection( Position, outerDirection, innerDirection, gradientWRTOuterDirection, fn );
      }
      if (this->_catchCtrlC && getCtrlCState())
        break;
    }
    this->unsetCtrlCHandler();
    cerr << endl;
  }
  
  void testRandomDirections ( const VecType &Position,
			      const int NumOfDirections,
                              const char* BasePlotFileName,
                              const bool applyDF = true,
                              const bool secondDerivTypeIsSymmetric = true,
			      const bool testNonNegEntriesOfHessianOnly = true ) const {
    if( applyDF )
      _DF.apply ( Position, *_pMatDF );
    if( !secondDerivTypeIsSymmetric )
      _pMatDF->transpose();

    VecType outerDirection( Position, aol::STRUCT_COPY );
    VecType innerDirection( Position, aol::STRUCT_COPY );
    VecType gradientWRTOuterDirection( Position, aol::STRUCT_COPY );

    aol::RandomGenerator rGen;
    rGen.randomize();
    aol::ProgressBar<> pb ( "Validating second derivatives in randomly chosen directions " );
    pb.start ( NumOfDirections );
    this->setCtrlCHandler();
    for( int i = 0; i < NumOfDirections; i++ ){      
      
      int outerDirectionNum = rGen.rInt( Position.getTotalSize() );
      int innerDirectionNum = rGen.rInt( Position.getTotalSize() );
      if (!this->wantsToCheckEntry(outerDirectionNum) || !this->wantsToCheckEntry(innerDirectionNum)) continue;
      
      outerDirection.setZero();
      innerDirection.setZero();
      
      aol::Vec< VecType::Depth, int>  idxVecOuter = this->setIthComponent ( outerDirectionNum, 1., outerDirection );
      aol::Vec< VecType::Depth, int>  idxVecInner = this->setIthComponent ( innerDirectionNum, 1., innerDirection );

      if( testNonNegEntriesOfHessianOnly && (abs( _pMatDF->getReference(idxVecOuter[0], idxVecInner[0]).get(idxVecOuter[1], idxVecInner[1]) ) < 1.e-10) ){
	i--;
	continue;
      }
      
      _pMatDF->apply ( outerDirection, gradientWRTOuterDirection );        
        
      char fn[1024];
      sprintf( fn, "%s_%s_%s", BasePlotFileName, this->getSuffixString( idxVecOuter ).c_str(), this->getSuffixString( idxVecInner ).c_str() );
      testDirection( Position, outerDirection, innerDirection, gradientWRTOuterDirection, fn );
      pb++;

      if (this->_catchCtrlC && getCtrlCState())
        break;
    }
    this->unsetCtrlCHandler();
    cerr << endl;
  }
};

} // end namespace aol

#endif
