#ifndef __CHANVESE_H
#define __CHANVESE_H

#include <FEOpInterface.h>
#include <RudinOsherFatemi.h>
#include <mcm.h>

namespace aol {

/**
 * \brief \f$ \left(\int_\Omega \frac{1}{H^\prime(\phi(x))}\varphi_i \varphi_j\right)_{ij}\f$
 * \author Berkels
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, typename HeavisideFunctionType>
class HeavisideFunctionMassOp : public aol::FELinScalarWeightedMassInterface<ConfiguratorType, HeavisideFunctionMassOp<ConfiguratorType, HeavisideFunctionType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrLevelsetFunction;
  const HeavisideFunctionType &_heavisideFunction;
public:
  HeavisideFunctionMassOp( const typename ConfiguratorType::InitType &Initializer,
                           const aol::Vector<RealType> &LevelsetFunction,
                           const HeavisideFunctionType &HeavisideFunction,
                           aol::OperatorType OpType = aol::ONTHEFLY )
    : aol::FELinScalarWeightedMassInterface<ConfiguratorType, HeavisideFunctionMassOp<ConfiguratorType, HeavisideFunctionType> >( Initializer, OpType ),
      _discrLevelsetFunction( Initializer, LevelsetFunction ),
      _heavisideFunction( HeavisideFunction )
  {
  }
  inline RealType getCoeff( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::VecType& ) const {
    return 1./(_heavisideFunction.evaluateDerivative(_discrLevelsetFunction.evaluateAtQuadPoint( El, QuadPoint)));
  }
};

/**
 * \brief Lumped version of HeavisideFunctionMassOp.
 * \author Berkels
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, typename HeavisideFunctionType>
class HeavisideFunctionLumpedMassOp : public aol::LumpedMassOpInterface<ConfiguratorType, HeavisideFunctionLumpedMassOp<ConfiguratorType, HeavisideFunctionType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrLevelsetFunction;
  const HeavisideFunctionType &_heavisideFunction;
public:
  HeavisideFunctionLumpedMassOp( const typename ConfiguratorType::InitType &Initializer,
                                 const aol::Vector<RealType> &LevelsetFunction,
                                 const HeavisideFunctionType &HeavisideFunction,
                                 const LUMPED_MASS_OP_MODE Invert )
    : aol::LumpedMassOpInterface<ConfiguratorType, HeavisideFunctionLumpedMassOp<ConfiguratorType, HeavisideFunctionType> >( Initializer, Invert ),
      _discrLevelsetFunction( Initializer, LevelsetFunction ),
      _heavisideFunction( HeavisideFunction )
  {
  }
  inline RealType getCoeff( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::VecType& ) const {
    return 1./(_heavisideFunction.evaluateDerivative(_discrLevelsetFunction.evaluateAtQuadPoint( El, QuadPoint)));
  }
};

/**
 * \brief \f$ \left(\int_\Omega H(\phi(x))\varphi_i \varphi_j\right)_{ij}\f$
 * \author Berkels
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, typename HeavisideFunctionType>
class HeavisideFunctionWeightedMassOp : public aol::FELinScalarWeightedMassInterface<ConfiguratorType, HeavisideFunctionWeightedMassOp<ConfiguratorType, HeavisideFunctionType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrLevelsetFunction;
  const HeavisideFunctionType &_heavisideFunction;
public:
  HeavisideFunctionWeightedMassOp( const typename ConfiguratorType::InitType &Initializer,
                                   const aol::Vector<RealType> &LevelsetFunction,
                                   const HeavisideFunctionType &HeavisideFunction,
                                   aol::OperatorType OpType = aol::ONTHEFLY )
    : aol::FELinScalarWeightedMassInterface<ConfiguratorType, HeavisideFunctionWeightedMassOp<ConfiguratorType, HeavisideFunctionType> >( Initializer, OpType ),
      _discrLevelsetFunction( Initializer, LevelsetFunction ),
      _heavisideFunction( HeavisideFunction )
  {
  }
  inline RealType getCoeff( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::VecType& ) const {
    return (_heavisideFunction.evaluate(_discrLevelsetFunction.evaluateAtQuadPoint( El, QuadPoint)));
  }
};

/**
 * \author Berkels
 * \warning The default value for \a Epsilon differs from the default value used in
 *          VariationOfHeavisideLevelsetLengthEnergy. This will be corrected
 *          in the future, but for legacy reasons is not yet done.
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int numberOfLevelsetFunctions = 1>
class HeavisideLevelsetLengthEnergy
: public aol::FENonlinIntegrationVectorInterface<ConfiguratorType,
                                                   HeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>,
                                                   numberOfLevelsetFunctions > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const HeavisideFunctionType &_heavisideFunction;
  const RealType _epsilonSqr;
public:
  HeavisideLevelsetLengthEnergy( const typename ConfiguratorType::InitType &Initializer,
                                 const HeavisideFunctionType &HeavisideFunction,
                                 const RealType Epsilon = 0. )
    : aol::FENonlinIntegrationVectorInterface< ConfiguratorType,
                                               HeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType, numberOfLevelsetFunctions>,
                                               numberOfLevelsetFunctions > ( Initializer ),
      _heavisideFunction ( HeavisideFunction ),
      _epsilonSqr ( aol::Sqr( Epsilon ) ) {
  }

  RealType evaluateIntegrand( const aol::auto_container<numberOfLevelsetFunctions,aol::DiscreteFunctionDefault<ConfiguratorType> > &discrFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    typename ConfiguratorType::VecType gradPhi;
    RealType integrand = 0.;
    for( int i = 0; i < numberOfLevelsetFunctions; i++ ){
      discrFuncs[i].evaluateGradientAtQuadPoint( El, QuadPoint, gradPhi );
      integrand += ( _heavisideFunction.evaluateDerivative(discrFuncs[i].evaluateAtQuadPoint( El, QuadPoint))
                     * sqrt ( gradPhi.normSqr() + _epsilonSqr ) );
    }
    return integrand;
  }
};

/**
 * Calulates the variation of HeavisideLevelsetLengthEnergy without the derivative
 * of the Heaviside function part. This has to be handled by the metric of the
 * gradient flow.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, int numberOfLevelsetFunctions = 1>
class VariationOfHeavisideLevelsetLengthEnergy
: public Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_grid;
  const RealType _epsilon;
  const RealType _scale;
public:
  /**
   * \param Epsilon parameter to regularize the absolute value
   * \param Scale scales the whole result MultiVector
   */
  VariationOfHeavisideLevelsetLengthEnergy( const typename ConfiguratorType::InitType &Initializer,
                                            const RealType Epsilon = 0.1,
                                            const RealType Scale = 1.)
    : _grid( Initializer ),
      _epsilon( Epsilon ),
      _scale( Scale )
  {
  }

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    qc::MCMStiffOp<ConfiguratorType> mcmStiffOp( _grid, aol::ONTHEFLY, _epsilon, _scale );
    for( int i = 0; i < numberOfLevelsetFunctions; i++ ){
      mcmStiffOp.setImageReference( MArg[i] );
      mcmStiffOp.applyAdd( MArg[i], MDest[i] );
    }
  }
};

/**
 * Computes \f$ \int H''(\phi(x))) \|\nabla\phi(x)\|\varphi_j dx \f$, where arg[0]=phi and HeavisideFunction=H.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType>
class FirstPartOfFullVariationOfHeavisideLevelsetLengthEnergy
  : public aol::FENonlinOpInterface< ConfiguratorType, FirstPartOfFullVariationOfHeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const HeavisideFunctionType &_heavisideFunction;
  const RealType _epsilonSqr;
public:
  FirstPartOfFullVariationOfHeavisideLevelsetLengthEnergy( const typename ConfiguratorType::InitType &Initializer,
                                                           const HeavisideFunctionType &HeavisideFunction,
                                                           const RealType Epsilon )
    : aol::FENonlinOpInterface< ConfiguratorType,
                                 FirstPartOfFullVariationOfHeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType> > ( Initializer ),
      _heavisideFunction( HeavisideFunction ),
      _epsilonSqr( aol::Sqr(Epsilon) ) {
  }

  void getNonlinearity( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/,
                        typename ConfiguratorType::RealType &NL ) const {
    const RealType hPrimePrime = _heavisideFunction.evaluateSecondDerivative(DiscFunc.evaluateAtQuadPoint( El, QuadPoint));

    typename ConfiguratorType::VecType gradPhi;
    DiscFunc.evaluateGradientAtQuadPoint( El, QuadPoint, gradPhi );

    NL = hPrimePrime * sqrt( gradPhi.normSqr() + _epsilonSqr );
  }
};

/**
 * Computes \f$ \int H'(\phi(x))) \frac{\nabla\phi(x)}{\|\nabla\phi(x)\|_\epsilon}\cdot \nabla\varphi_j dx \f$, where arg[0]=phi and HeavisideFunction=H.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType>
class SecondPartOfFullVariationOfHeavisideLevelsetLengthEnergy
: public aol::FENonlinDiffOpInterface<ConfiguratorType, SecondPartOfFullVariationOfHeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const HeavisideFunctionType &_heavisideFunction;
  const RealType _epsilonSqr;
public:
  SecondPartOfFullVariationOfHeavisideLevelsetLengthEnergy( const typename ConfiguratorType::InitType &Initializer,
                                                            const HeavisideFunctionType &HeavisideFunction,
                                                            const RealType Epsilon )
  : aol::FENonlinDiffOpInterface< ConfiguratorType,
                                   SecondPartOfFullVariationOfHeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType> > ( Initializer ),
      _heavisideFunction( HeavisideFunction ),
      _epsilonSqr( aol::Sqr(Epsilon) ) {
  }

  inline void getNonlinearity( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint,
                               const typename ConfiguratorType::DomVecType &/*RefCoord*/,
                               typename ConfiguratorType::VecType &NL ) const {
    const RealType hPrime = _heavisideFunction.evaluateDerivative(DiscFunc.evaluateAtQuadPoint( El, QuadPoint));
    DiscFunc.evaluateGradientAtQuadPoint( El, QuadPoint, NL );
    NL *= hPrime / (sqrt(NL.normSqr()+_epsilonSqr));
  }
};

/**
 * Calulates the variation of HeavisideLevelsetLengthEnergy including the derivative
 * of the Heaviside function part.
 *
 * This involves the second derivate of the Heaviside function. Therefore the modelling
 * should try to avoid to use this. Currently it's only meant to be used for testing
 * purposes.
 *
 * Instead of using this class, let the metric of the gradient flow handle the Heaviside
 * function part of the variation and use VariationOfHeavisideLevelsetLengthEnergy.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int numberOfLevelsetFunctions = 1>
class FullVariationOfHeavisideLevelsetLengthEnergy
: public Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_grid;
  const HeavisideFunctionType &_heavisideFunction;
  const RealType _epsilon;
public:
  /**
   * \param Epsilon parameter to regularize the absolute value
   */
  FullVariationOfHeavisideLevelsetLengthEnergy( const typename ConfiguratorType::InitType &Initializer,
                                                const HeavisideFunctionType &HeavisideFunction,
                                                const RealType Epsilon = 0.1 )
    : _grid( Initializer ),
      _heavisideFunction( HeavisideFunction ),
      _epsilon( Epsilon )
  {
  }

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {

    FirstPartOfFullVariationOfHeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType>
      part1( _grid, _heavisideFunction, _epsilon );
    SecondPartOfFullVariationOfHeavisideLevelsetLengthEnergy<ConfiguratorType, HeavisideFunctionType>
      part2( _grid, _heavisideFunction, _epsilon );

    for( int i = 0; i < numberOfLevelsetFunctions; i++ ){
      part1.applyAdd( MArg[i], MDest[i] );
      part2.applyAdd( MArg[i], MDest[i] );
    }
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, typename TargetOpType>
class HeavisideOpBase {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrLevelsetFunction;
  const HeavisideFunctionType &_heavisideFunction;
  const TargetOpType &_targetOp;
public:
  HeavisideOpBase( const typename ConfiguratorType::InitType &Initializer,
                   const aol::Vector<RealType> &LevelsetFunction,
                   const HeavisideFunctionType &HeavisideFunction,
                   const TargetOpType &TargetOp )
    : _discrLevelsetFunction( Initializer, LevelsetFunction ),
      _heavisideFunction( HeavisideFunction ),
      _targetOp( TargetOp ){
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType, typename TargetOpType>
class HeavisideFENonlinIntegrationVectorInterface
: public aol::FENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType>, TargetOpType::NumOfComponents >,
         aol::HeavisideOpBase<ConfiguratorType, HeavisideFunctionType, TargetOpType> {
public:
  typedef typename ConfiguratorType::RealType RealType;

  HeavisideFENonlinIntegrationVectorInterface( const typename ConfiguratorType::InitType &Initializer,
                                               const aol::Vector<RealType> &LevelsetFunction,
                                               const HeavisideFunctionType &HeavisideFunction,
                                               const TargetOpType &TargetOp )
  : aol::FENonlinIntegrationVectorInterface< ConfiguratorType,
                                               HeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType>,
                                               TargetOpType::NumOfComponents > ( Initializer ),
    aol::HeavisideOpBase<ConfiguratorType, HeavisideFunctionType, TargetOpType> ( Initializer, LevelsetFunction, HeavisideFunction, TargetOp ) {
  }

  RealType evaluateIntegrand( const aol::auto_container<TargetOpType::NumOfComponents,aol::DiscreteFunctionDefault<ConfiguratorType> > &discrFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint, const typename ConfiguratorType::VecType &RefCoord ) const {
    return ( this->_heavisideFunction.evaluate ( this->_discrLevelsetFunction.evaluateAtQuadPoint( El, QuadPoint) )
             * this->_targetOp.evaluateIntegrand( discrFuncs, El, QuadPoint, RefCoord ) );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, typename TargetOpType>
class HeavisideFENonlinIntegrationScalarInterface
: public aol::FENonlinIntegrationScalarInterface<ConfiguratorType, HeavisideFENonlinIntegrationScalarInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType> >,
         aol::HeavisideOpBase<ConfiguratorType, HeavisideFunctionType, TargetOpType> {
public:
  typedef typename ConfiguratorType::RealType RealType;

  HeavisideFENonlinIntegrationScalarInterface ( const typename ConfiguratorType::InitType &Initializer,
                                                const aol::Vector<RealType> &LevelsetFunction,
                                                const HeavisideFunctionType &HeavisideFunction,
                                                const TargetOpType &TargetOp )
    : aol::FENonlinIntegrationScalarInterface< ConfiguratorType,
                                               HeavisideFENonlinIntegrationScalarInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType> > ( Initializer ),
      aol::HeavisideOpBase<ConfiguratorType, HeavisideFunctionType, TargetOpType> ( Initializer, LevelsetFunction, HeavisideFunction, TargetOp ) { }

  RealType evaluateIntegrand( const aol::DiscreteFunctionDefault< ConfiguratorType > &DiscFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    return ( this->_heavisideFunction.evaluate ( this->_discrLevelsetFunction.evaluateAtQuadPoint( El, QuadPoint) )
             * this->_targetOp.evaluateIntegrand( DiscFuncs, El, QuadPoint, RefCoord ) );
  }
};

/**
 * \brief Multiplies the integrand of ops derived from aol::FELinScalarWeightedStiffInterface with \f$H(\phi(x))\f$.
 * \author Berkels
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, typename TargetOpType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class HeavisideFELinScalarWeightedStiffInterface
: public aol::FELinScalarWeightedStiffInterface<ConfiguratorType, HeavisideFELinScalarWeightedStiffInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType, IndexMode>, IndexMode >,
         aol::HeavisideOpBase<ConfiguratorType, HeavisideFunctionType, TargetOpType> {

public:
  typedef typename ConfiguratorType::RealType RealType;

  HeavisideFELinScalarWeightedStiffInterface ( const typename ConfiguratorType::InitType &Initializer,
                                               const aol::Vector<RealType> &LevelsetFunction,
                                               const HeavisideFunctionType &HeavisideFunction,
                                               const TargetOpType &TargetOp,
                                               OperatorType OpType = ONTHEFLY )
    : aol::FELinScalarWeightedStiffInterface<ConfiguratorType, HeavisideFELinScalarWeightedStiffInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType, IndexMode>, IndexMode > ( Initializer, OpType ),
      aol::HeavisideOpBase<ConfiguratorType, HeavisideFunctionType, TargetOpType> ( Initializer, LevelsetFunction, HeavisideFunction, TargetOp ) { }

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::VecType &RefCoord ) const {
    return ( this->_heavisideFunction.evaluate ( this->_discrLevelsetFunction.evaluateAtQuadPoint( El, QuadPoint) )
             * this->_targetOp.getCoeff( El, QuadPoint, RefCoord ) );
  }
};

/**
 * \note The corresponding energy can be constructed using aol::IsoEnergyOp and aol::HeavisideFENonlinIntegrationScalarInterface.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType>
class VariationOfHeavisideIsoEnergyOp : public aol::VariationOfIsoEnergyOp<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::ArrayType &_levelsetFunction;
  const HeavisideFunctionType &_heavisideFunction;
public:
  VariationOfHeavisideIsoEnergyOp ( const typename ConfiguratorType::InitType &Initializer,
                                    const typename ConfiguratorType::ArrayType &LevelsetFunction,
                                    const HeavisideFunctionType &HeavisideFunction,
                                    const RealType Delta )
  : aol::VariationOfIsoEnergyOp<ConfiguratorType> ( Initializer, Delta ),
    _levelsetFunction ( LevelsetFunction ),
    _heavisideFunction ( HeavisideFunction ) { }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    qc::MCMStiffOp<ConfiguratorType> mcmStiffOp( this->_grid, aol::ONTHEFLY, this->_delta );
    mcmStiffOp.setImageReference( Arg );
    aol::HeavisideFELinScalarWeightedStiffInterface<ConfiguratorType, HeavisideFunctionType, qc::MCMStiffOp<ConfiguratorType> > hMcmStiffOp ( this->_grid, _levelsetFunction, _heavisideFunction, mcmStiffOp, aol::ONTHEFLY );
    hMcmStiffOp.applyAdd(Arg, Dest);
  }

};

template <typename ConfiguratorType, typename HeavisideFunctionType, typename TargetOpType>
class OneMinusHeavisideFENonlinIntegrationVectorInterface
: public aol::FENonlinIntegrationVectorInterface<ConfiguratorType, OneMinusHeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType>, TargetOpType::NumOfComponents > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrLevelsetFunction;
  const HeavisideFunctionType &_heavisideFunction;
  const TargetOpType &_targetOp;
public:
  OneMinusHeavisideFENonlinIntegrationVectorInterface( const typename ConfiguratorType::InitType &Initializer,
                                                         const aol::Vector<RealType> &LevelsetFunction,
                                                         const HeavisideFunctionType &HeavisideFunction,
                                                         const TargetOpType &TargetOp )
  : aol::FENonlinIntegrationVectorInterface< ConfiguratorType,
                                               OneMinusHeavisideFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType>,
                                               TargetOpType::NumOfComponents > ( Initializer ),
      _discrLevelsetFunction( Initializer, LevelsetFunction ),
      _heavisideFunction( HeavisideFunction ),
      _targetOp( TargetOp ){
  }

  RealType evaluateIntegrand( const aol::auto_container<TargetOpType::NumOfComponents,aol::DiscreteFunctionDefault<ConfiguratorType> > &discrFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint, const typename ConfiguratorType::VecType &RefCoord ) const {
    return ( (1. -_heavisideFunction.evaluate(_discrLevelsetFunction.evaluateAtQuadPoint( El, QuadPoint)))
             * _targetOp.evaluateIntegrand( discrFuncs, El, QuadPoint, RefCoord ) );
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType, typename TargetOpType>
class HeavisideAndOneMinusHFENonlinIntegrationVectorInterface
: public aol::FENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideAndOneMinusHFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType>, TargetOpType::NumOfComponents > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrLevelsetFunction;
  const HeavisideFunctionType &_heavisideFunction;
  const TargetOpType &_hsTargetOp;
  const TargetOpType &_oneMinusHSTargetOp;
public:
  HeavisideAndOneMinusHFENonlinIntegrationVectorInterface( const typename ConfiguratorType::InitType &Initializer,
                                                             const aol::Vector<RealType> &LevelsetFunction,
                                                             const HeavisideFunctionType &HeavisideFunction,
                                                             const TargetOpType &HSTargetOp,
                                                             const TargetOpType &OneMinusHSTargetOp )
  : aol::FENonlinIntegrationVectorInterface< ConfiguratorType,
                                               HeavisideAndOneMinusHFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType>,
                                               TargetOpType::NumOfComponents > ( Initializer ),
      _discrLevelsetFunction( Initializer, LevelsetFunction ),
      _heavisideFunction( HeavisideFunction ),
      _hsTargetOp( HSTargetOp ),
      _oneMinusHSTargetOp( OneMinusHSTargetOp ){
  }

  RealType evaluateIntegrand( const aol::auto_container<TargetOpType::NumOfComponents,aol::DiscreteFunctionDefault<ConfiguratorType> > &discrFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint, const typename ConfiguratorType::VecType &RefCoord ) const {
    const RealType h = _heavisideFunction.evaluate(_discrLevelsetFunction.evaluateAtQuadPoint( El, QuadPoint));
    const RealType hsOp = _hsTargetOp.evaluateIntegrand( discrFuncs, El, QuadPoint, RefCoord );
    const RealType oneMinusHSOp = _oneMinusHSTargetOp.evaluateIntegrand( discrFuncs, El, QuadPoint, RefCoord );
    return ( h * hsOp + ( 1. - h ) * oneMinusHSOp );
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType, typename TargetOpType>
class OneMinusHeavisideFENonlinOpInterface : public aol::FENonlinOpInterface< ConfiguratorType, OneMinusHeavisideFENonlinOpInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType> > {
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrLevelsetFunction;
  const HeavisideFunctionType &_heavisideFunction;
  const TargetOpType &_targetOp;
public:
  OneMinusHeavisideFENonlinOpInterface( const typename ConfiguratorType::InitType &Initializer,
                                        const aol::Vector<RealType> &LevelsetFunction,
                                        const HeavisideFunctionType &HeavisideFunction,
                                        const TargetOpType &TargetOp )
    : aol::FENonlinOpInterface< ConfiguratorType,
                                 OneMinusHeavisideFENonlinOpInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType> > ( Initializer ),
      _discrLevelsetFunction( Initializer, LevelsetFunction ),
      _heavisideFunction( HeavisideFunction ),
      _targetOp( TargetOp ){
  }

  void getNonlinearity( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,
                        typename ConfiguratorType::RealType &NL ) const {
    _targetOp.getNonlinearity( DiscFunc, El, QuadPoint, RefCoord, NL );
    NL *= (1. -_heavisideFunction.evaluate(_discrLevelsetFunction.evaluateAtQuadPoint( El, QuadPoint)));
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType, typename TargetOpType>
class HeavySideOneMinusHFENonlinOpInterface : public aol::FENonlinOpInterface< ConfiguratorType, HeavySideOneMinusHFENonlinOpInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType> > {
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrLevelsetFunction;
  const HeavisideFunctionType &_heavisideFunction;
  const TargetOpType &_hsTargetOp;
  const TargetOpType &_oneMinusHSTargetOp;
public:
  HeavySideOneMinusHFENonlinOpInterface( const typename ConfiguratorType::InitType &Initializer,
                                         const aol::Vector<RealType> &LevelsetFunction,
                                         const HeavisideFunctionType &HeavisideFunction,
                                         const TargetOpType &HSTargetOp,
                                         const TargetOpType &OneMinusHSTargetOp )
    : aol::FENonlinOpInterface< ConfiguratorType,
                                 HeavySideOneMinusHFENonlinOpInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType> > ( Initializer ),
      _discrLevelsetFunction( Initializer, LevelsetFunction ),
      _heavisideFunction( HeavisideFunction ),
      _hsTargetOp( HSTargetOp ),
      _oneMinusHSTargetOp( OneMinusHSTargetOp ){
  }

  void getNonlinearity( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,
                        typename ConfiguratorType::RealType &NL ) const {
    const RealType h = _heavisideFunction.evaluate(_discrLevelsetFunction.evaluateAtQuadPoint( El, QuadPoint));
    RealType hsOp = 0.;
    _hsTargetOp.getNonlinearity( DiscFunc, El, QuadPoint, RefCoord, hsOp );
    RealType oneMinusHSOp = 0.;
    _oneMinusHSTargetOp.getNonlinearity( DiscFunc, El, QuadPoint, RefCoord, oneMinusHSOp );
    NL = ( h * hsOp + ( 1. - h ) * oneMinusHSOp );
  }
};

template <typename ConfiguratorType, typename HeavisideFunctionType, typename TargetOpType>
class HeavisideAndOneMinusHFENonlinDiffOpInterface
: public aol::FENonlinDiffOpInterface<ConfiguratorType, HeavisideAndOneMinusHFENonlinDiffOpInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrLevelsetFunction;
  const HeavisideFunctionType &_heavisideFunction;
  const TargetOpType &_hsTargetOp;
  const TargetOpType &_oneMinusHSTargetOp;
public:
  HeavisideAndOneMinusHFENonlinDiffOpInterface( const typename ConfiguratorType::InitType &Initializer,
                                                 const aol::Vector<RealType> &LevelsetFunction,
                                                 const HeavisideFunctionType &HeavisideFunction,
                                                 const TargetOpType &HSTargetOp,
                                                 const TargetOpType &OneMinusHSTargetOp )
  : aol::FENonlinDiffOpInterface< ConfiguratorType,
                                   HeavisideAndOneMinusHFENonlinDiffOpInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType> > ( Initializer ),
      _discrLevelsetFunction( Initializer, LevelsetFunction ),
      _heavisideFunction( HeavisideFunction ),
      _hsTargetOp( HSTargetOp ),
      _oneMinusHSTargetOp( OneMinusHSTargetOp ){
  }

  inline void getNonlinearity( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint,
                               const typename ConfiguratorType::VecType &RefCoord,
                               typename ConfiguratorType::VecType &NL ) const {
    const RealType h = _heavisideFunction.evaluate(_discrLevelsetFunction.evaluateAtQuadPoint( El, QuadPoint));
    typename ConfiguratorType::VecType temp;

    _hsTargetOp.getNonlinearity( DiscFunc, El, QuadPoint, RefCoord, temp );
    temp *= h;
    NL = temp;

    _oneMinusHSTargetOp.getNonlinearity( DiscFunc, El, QuadPoint, RefCoord, temp );
    temp *= (1.-h);
    NL += temp;
  }
};

template <typename ConfiguratorType, int NumCompArg, int NumCompDest, typename HeavisideFunctionType, typename TargetOpType>
class OneMinusHeavisideFENonlinVectorDiffOpInterface : public aol::FENonlinVectorDiffOpInterface< ConfiguratorType, NumCompArg, NumCompDest, OneMinusHeavisideFENonlinVectorDiffOpInterface<ConfiguratorType, NumCompArg, NumCompDest, HeavisideFunctionType, TargetOpType> > {
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrLevelsetFunction;
  const HeavisideFunctionType &_heavisideFunction;
  const TargetOpType &_targetOp;
public:
  OneMinusHeavisideFENonlinVectorDiffOpInterface ( const typename ConfiguratorType::InitType &Initializer,
                                                   const aol::Vector<RealType> &LevelsetFunction,
                                                   const HeavisideFunctionType &HeavisideFunction,
                                                   const TargetOpType &TargetOp )
    : aol::FENonlinVectorDiffOpInterface< ConfiguratorType,
                                          NumCompArg,
                                          NumCompDest,
                                          OneMinusHeavisideFENonlinVectorDiffOpInterface<ConfiguratorType, NumCompArg, NumCompDest, HeavisideFunctionType, TargetOpType> > ( Initializer ),
      _discrLevelsetFunction( Initializer, LevelsetFunction ),
      _heavisideFunction( HeavisideFunction ),
      _targetOp( TargetOp ){
  }

  void getNonlinearity ( const auto_container<NumCompArg, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,
                         aol::Mat<NumCompDest, NumCompDest, typename ConfiguratorType::RealType> &NL ) const {
    _targetOp.getNonlinearity( DiscFuncs, El, QuadPoint, RefCoord, NL );
    NL *= (1. -_heavisideFunction.evaluate(_discrLevelsetFunction.evaluateAtQuadPoint( El, QuadPoint)));
  }
};

/**
 * Computes \f$ \int (1-H(\phi(x))) f(x) dx \f$, where arg[0]=f, LevelsetFunction=phi and HeavisideFunction=H.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType>
class OneMinusHeavisideIntegrateFEFunction
  : public FENonlinIntegrationScalarInterface<ConfiguratorType, OneMinusHeavisideIntegrateFEFunction<ConfiguratorType, HeavisideFunctionType> > {

  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrLevelsetFunction;
  const HeavisideFunctionType &_heavisideFunction;

public:
  OneMinusHeavisideIntegrateFEFunction ( const typename ConfiguratorType::InitType &Initializer,
                                         const aol::Vector<RealType> &LevelsetFunction,
                                         const HeavisideFunctionType &HeavisideFunction )
    : FENonlinIntegrationScalarInterface<ConfiguratorType, OneMinusHeavisideIntegrateFEFunction<ConfiguratorType, HeavisideFunctionType> > ( Initializer ),
      _discrLevelsetFunction( Initializer, LevelsetFunction ),
      _heavisideFunction( HeavisideFunction ) { }

  RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    const RealType h = _heavisideFunction.evaluate(_discrLevelsetFunction.evaluateAtQuadPoint( El, QuadPoint));
    return (ZOTrait<RealType>::one-h) * DiscFunc.evaluateAtQuadPoint ( El, QuadPoint );
  }
};

/**
 * \brief Computes \f$ - \int H'(\phi(x)) f(x) \varphi_j dx \f$, where arg[0]=f, LevelsetFunction=phi and HeavisideFunction=H.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType>
class FullVariationOfOneMinusHeavisideIntegrateFEFunction
  : public aol::FENonlinOpInterface< ConfiguratorType, FullVariationOfOneMinusHeavisideIntegrateFEFunction<ConfiguratorType, HeavisideFunctionType> > {

  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrLevelsetFunction;
  const HeavisideFunctionType &_heavisideFunction;

public:
  FullVariationOfOneMinusHeavisideIntegrateFEFunction ( const typename ConfiguratorType::InitType &Initializer,
                                                        const aol::Vector<RealType> &LevelsetFunction,
                                                        const HeavisideFunctionType &HeavisideFunction )
    : aol::FENonlinOpInterface< ConfiguratorType, FullVariationOfOneMinusHeavisideIntegrateFEFunction<ConfiguratorType, HeavisideFunctionType> > ( Initializer ),
      _discrLevelsetFunction( Initializer, LevelsetFunction ),
      _heavisideFunction( HeavisideFunction ) { }

  void getNonlinearity( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
                        typename ConfiguratorType::RealType &NL ) const {
    const RealType hPrime = _heavisideFunction.evaluateDerivative(_discrLevelsetFunction.evaluateAtQuadPoint( El, QuadPoint));
    NL = -1 * hPrime * DiscFunc.evaluateAtQuadPoint ( El, QuadPoint );
  }
};

/**
 * \brief \f$f(x)=x^2\f$
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template <typename RealType>
class SquareFunction{
public:
  static RealType evaluate( const RealType x ) {
    return aol::Sqr ( x );
  }
  static RealType evaluateDerivative( const RealType x ) {
    return 2 * x;
  }
};

/**
 * \brief \f$f(x)=x\f$. Heaviside function replacement to convert the Chan Vese model to the one of Chan Esedoglu Nikolova with exact penalty.
 *
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template <typename RealType, typename _ScalingFunctionType = aol::SquareFunction<RealType> >
class IdentityFunction{
public:
  typedef _ScalingFunctionType ScalingFunctionType;

  IdentityFunction( const RealType /*Epsilon*/ = 0. ) {
  }
  static RealType evaluate( const RealType x ) {
    return x;
  }
  static RealType evaluateDerivative( const RealType /*x*/ ) {
    return aol::NumberTrait<RealType>::one;
  }
  //! In a proper Chan Vese type model, you should not need the second derivative of the Heaviside function.
  static RealType evaluateSecondDerivative( const RealType /*x*/ ) {
    return aol::NumberTrait<RealType>::zero;
  }
};

/**
 * \brief \f$ H_\epsilon(x)=\frac12+\frac{1}{\pi}\arctan\left(\frac{x}{\epsilon}\right)\f$.
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template <typename RealType>
class ArcTanHeavisideFunction{
  const RealType _epsilon;
  const RealType _oneOverPi;
public:
  typedef IdentityFunction<RealType> ScalingFunctionType;

  ArcTanHeavisideFunction( const RealType Epsilon )
   : _epsilon(Epsilon),
    _oneOverPi( 1./aol::NumberTrait<RealType>::getPi() )
  {
  }
  RealType evaluate( const RealType x ) const{
    return 0.5 + _oneOverPi*atan( x/_epsilon );
  }
  RealType evaluateDerivative( const RealType x ) const{
    return _oneOverPi/(_epsilon +  aol::Sqr(x)/_epsilon );
  }
  //! In a proper Chan Vese type model, you should not need the second derivative of the Heaviside function.
  RealType evaluateSecondDerivative( const RealType x ) const{
    return -2*x*_oneOverPi/(_epsilon*aol::Sqr(_epsilon +  aol::Sqr(x)/_epsilon ));
  }
};

/**
 * \brief \f$ f(x)=1-H_\epsilon(x)=\frac12-\frac{1}{\pi}\arctan\left(\frac{x}{\epsilon}\right)\f$.
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template <typename RealType>
class OneMinusArcTanHeavisideFunction{
  const RealType _epsilon;
  const RealType _oneOverPi;
public:
  OneMinusArcTanHeavisideFunction( const RealType Epsilon )
   : _epsilon(Epsilon),
     _oneOverPi( 1./aol::NumberTrait<RealType>::getPi() )
  {
  }
  RealType evaluate( const RealType x ) const{
    return 0.5 - _oneOverPi*atan( x/_epsilon );
  }
  RealType evaluateDerivative( const RealType x ) const{
    return -1.*_oneOverPi/(_epsilon +  aol::Sqr(x)/_epsilon );
  }
};

/**
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template <typename RealType>
class PolynomialHeavisideFunction{
  const RealType _epsilon;
public:
  PolynomialHeavisideFunction( const RealType Epsilon )
   : _epsilon( Epsilon )
  {
  }
  RealType evaluate( const RealType x ) const{
    if( x > _epsilon*0.5){
      return 1.;
    }
    else{
      if( x < (-1.)*_epsilon*0.5 ){
        return 0.;
      }
      else{
        const RealType position = -x/_epsilon+0.5;
        return aol::Sqr(position)*(2.*position-3.)+1.;
      }
    }
  }
  RealType evaluateDerivative( const RealType x ) const{
    if( fabs(x) > _epsilon*0.5){
      return 0.;
    }
    else{
      const RealType position = -x/_epsilon+0.5;
      return -6.*position*(position-1.)/_epsilon;
    }
  }
};

/**
 * \brief Smooth approximation (two times continuously differentiable) of the Heaviside function
 *        based on a polynomial of degree 5.
 *
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template <typename RealType>
class C2PolynomialHeavisideFunction{
  const RealType _epsilon;
public:
  C2PolynomialHeavisideFunction( const RealType Epsilon )
    : _epsilon( Epsilon ) { }

  RealType evaluate( const RealType x ) const {
    if ( x > _epsilon*0.5 ) {
      return 1.;
    }
    else {
      if ( x < (-1.)*_epsilon*0.5 ) {
        return 0.;
      }
      else{
        const RealType position = x / _epsilon + 0.5;
        return aol::Cub ( position ) * ( ( 6 * position - 15 ) * position + 10 );
      }
    }
  }

  RealType evaluateDerivative( const RealType x ) const {
    if ( fabs(x) > _epsilon*0.5 ) {
      return 0.;
    }
    else {
      const RealType position = x / _epsilon + 0.5;
      return 30 * aol::Sqr ( position ) * ( ( position - 2 ) * position + 1 ) / _epsilon;
    }
  }

  RealType evaluateSecondDerivative( const RealType x ) const{
    if ( fabs(x) > _epsilon*0.5 ) {
      return 0.;
    }
    else {
      const RealType position = x / _epsilon + 0.5;
      return 60 * position * ( ( 2 * position - 3 ) * position + 1 ) / aol::Sqr ( _epsilon );
    }
  }
};

//! Prints central differential quotioent and implemented derivative of the Heaviside Funtion.
//! If they don't match, it's likely that there is a bug in the Heaviside Funtion.
//! Also calculates the integral over the derivative of the Heaviside Function.
//! This should be one, since the derivatives approximates the delta distribution.
template <typename RealType, typename HeavisideFunctionType>
void consistencyCheckOfHeavisideFuntion( HeavisideFunctionType &H ){
  const RealType epsilon = 0.0001;
  aol::MixedFormat format ( 4, 12 );
  for( int i = -100; i < 100; i++ ){
    const RealType x = static_cast<RealType>(i)/1000.;
    const RealType differentialQuotient = (H.evaluate( x + epsilon ) - H.evaluate( x - epsilon ))/ (2.*epsilon);
    const RealType derivative = H.evaluateDerivative( x );
    cerr << format(differentialQuotient) << " " << format(derivative) << " " << format(differentialQuotient-derivative) << endl;
  }
  const RealType width = 1000.;
  const int N = 1000000;
  RealType integral = 0.;
  const RealType h = 2.*width/static_cast<RealType>(N-1);
  for( int i = 0; i < N; i++ ){
    const RealType x = -width + static_cast<RealType>(i)*h;
    integral += h * H.evaluateDerivative( x );
  }
  cerr << "Integral of derivative = " << integral << endl;
}


/**
 * \brief Heaviside function replacement to convert the Chan Vese model to the one of Chan Esedoglu Nikolova.
 *
 * \attention You can't use the Armijo rule in a gradient descent for an energy using this function:
 * evaluate is clamped, evaluateDerivative is not, in this sense the derivative doesn't match the energy
 * and therefore the Armijo rule will reject the negative gradient as descent direction.
 *
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template <typename RealType>
class ClampedIdentityFunction{
public:
  ClampedIdentityFunction( const RealType /*Epsilon*/ = 0. ) {
  }
  static RealType evaluate( const RealType x ) {
    return aol::Clamp<RealType>( x, 0, 1 );
  }
  static RealType evaluateDerivative( const RealType /*x*/ ) {
    return aol::NumberTrait<RealType>::one;
  }
  //! In a proper Chan Vese type model, you should not need the second derivative of the Heaviside function.
  static RealType evaluateSecondDerivative( const RealType /*x*/ ) {
    return aol::NumberTrait<RealType>::zero;
  }
};

/**
 * \brief Regularized version of \f$ \nu(x)=\max\{0,2|x-\frac{1}{2}|-1\}\f$.
 *
 * Used as penalty term in the Chan Esedoglu Nikolova model.
 *
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template <typename RealType>
class ZeroOneIntervalPenaltyFunction {
  const RealType _oneOverFourEpsilon, _oneOverTwoEpsilon, _epsilonMinusOneHalf, _oneHalfMinusEpsilon, _oneHalfPlusEpsilon;
public:
  ZeroOneIntervalPenaltyFunction ( const RealType Epsilon )
   : _oneOverFourEpsilon( aol::NumberTrait<RealType>::one/(4*Epsilon) ),
     _oneOverTwoEpsilon( aol::NumberTrait<RealType>::one/(2*Epsilon) ),
     _epsilonMinusOneHalf ( Epsilon - 0.5 ),
     _oneHalfMinusEpsilon ( 0.5 - Epsilon ),
     _oneHalfPlusEpsilon ( 0.5 + Epsilon )
  {
  }
  RealType evaluate( const RealType x ) const {
    const RealType y = aol::Abs( x - 0.5 );
    if ( y < _oneHalfMinusEpsilon )
      return aol::NumberTrait<RealType>::zero;
    else if ( y > _oneHalfPlusEpsilon )
      return ( y - 0.5 );
    else
      return _oneOverFourEpsilon * aol::Sqr( y + _epsilonMinusOneHalf );
  }
  RealType evaluateDerivative( const RealType x ) const{
    const RealType y = aol::Abs( x - 0.5 );
    if ( y < _oneHalfMinusEpsilon )
      return aol::NumberTrait<RealType>::zero;
    else if ( y > _oneHalfPlusEpsilon )
      return ( aol::signum( x - 0.5 ) * aol::NumberTrait<RealType>::one );
    else
      return aol::signum( x - 0.5 ) * _oneOverTwoEpsilon * ( y + _epsilonMinusOneHalf );
  }
};

/**
 * \brief \f$ \nu(x)=\max\{0,-x|x|, (x-1)|x-1|\} \f$.
 *
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template <typename RealType>
class QuadraticZeroOneIntervalPenaltyFunction {
public:
  QuadraticZeroOneIntervalPenaltyFunction ( const RealType /*Epsilon*/ ) {}

  RealType evaluate( const RealType x ) const {
    const RealType y = aol::Abs( x - 0.5 );
    if ( y < 0.5 )
      return aol::NumberTrait<RealType>::zero;
    else
      return ( aol::Sqr ( y - 0.5 ) );
  }
  RealType evaluateDerivative( const RealType x ) const{
    const RealType y = aol::Abs( x - 0.5 );
    if ( y < 0.5 )
      return aol::NumberTrait<RealType>::zero;
    else
      return aol::signum( x - 0.5 ) * 2 * ( y -0.5 );
  }
};

/**
 * \brief \f$ \nu(x)=\max\{0,-x^3, (x-1)^3\} \f$.
 *
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template <typename RealType>
class CubicZeroOneIntervalPenaltyFunction {
public:
  CubicZeroOneIntervalPenaltyFunction ( const RealType /*Epsilon*/ ) {}
  RealType evaluate( const RealType x ) const {
    const RealType y = aol::Abs( x - 0.5 );
    if ( y < 0.5 )
      return aol::NumberTrait<RealType>::zero;
    else
      return ( aol::Cub ( y - 0.5 ) );
  }
  RealType evaluateDerivative( const RealType x ) const{
    const RealType y = aol::Abs( x - 0.5 );
    if ( y < 0.5 )
      return aol::NumberTrait<RealType>::zero;
    else
      return aol::signum( x - 0.5 ) * 3 * aol::Sqr( y -0.5 );
  }
  RealType evaluateSecondDerivative( const RealType x ) const {
    const RealType y = aol::Abs( x - 0.5 );
    if ( y < 0.5 )
      return aol::NumberTrait<RealType>::zero;
    else
      return aol::signum( x - 0.5 ) * 6 * ( y -0.5 );
  }
};

/**
 * \brief \f$ \nu(x)=max\{0,-x^3\} \f$.
 *
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template <typename RealType>
class CubicPositivityPenaltyFunction {
public:
  CubicPositivityPenaltyFunction () {}

  RealType evaluate( const RealType x ) const {
    if ( x > 0 )
      return aol::NumberTrait<RealType>::zero;
    else
      return ( - aol::Cub ( x ) );
  }

  RealType evaluateDerivative( const RealType x ) const{
    if ( x > 0 )
      return aol::NumberTrait<RealType>::zero;
    else
      return - 3 * aol::Sqr( x );
  }

  RealType evaluateSecondDerivative( const RealType x ) const {
    if ( x > 0 )
      return aol::NumberTrait<RealType>::zero;
    else
      return - 6 * ( x );
  }
};

/**
 * \brief \f$ \nu(x)=max\{0,x^3\} \f$.
 *
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template <typename RealType>
class CubicNegativityPenaltyFunction {
public:
  CubicNegativityPenaltyFunction () {}

  RealType evaluate( const RealType x ) const {
    if ( x < 0 )
      return aol::NumberTrait<RealType>::zero;
    else
      return ( aol::Cub ( x ) );
  }

  RealType evaluateDerivative( const RealType x ) const{
    if ( x < 0 )
      return aol::NumberTrait<RealType>::zero;
    else
      return 3 * aol::Sqr( x );
  }

  RealType evaluateSecondDerivative( const RealType x ) const {
    if ( x < 0 )
      return aol::NumberTrait<RealType>::zero;
    else
      return 6 * ( x );
  }
};

/**
 * \brief \f$ \nu(u,v,c)=p(u/v-c) \f$.
 *
 * \author Berkels
 * \ingroup AnalyticFunctions
 */
template <typename RealType>
class CubicRatioPenaltyFunction {
  const RealType _c;
  const aol::CubicPositivityPenaltyFunction<RealType> _penaltyFunc;
public:
  CubicRatioPenaltyFunction ( const RealType C ) : _c ( C ) {}

  RealType evaluate( const RealType u, const RealType v ) const {
    return _penaltyFunc.evaluate ( u / v - _c );
  }

  RealType evaluateUDerivative( const RealType u, const RealType v ) const{
    return _penaltyFunc.evaluateDerivative ( u / v - _c ) / v;
  }

  RealType evaluateVDerivative( const RealType u, const RealType v ) const{
    return - _penaltyFunc.evaluateDerivative ( u / v - _c ) * ( u / aol::Sqr ( v ) );
  }

  RealType evaluateUUDerivative( const RealType u, const RealType v ) const{
    return _penaltyFunc.evaluateSecondDerivative ( u / v - _c ) / aol::Sqr ( v );
  }

  RealType evaluateUVDerivative( const RealType u, const RealType v ) const{
    return - _penaltyFunc.evaluateSecondDerivative ( u / v - _c ) * ( u / aol::Cub ( v ) )
           - _penaltyFunc.evaluateDerivative ( u / v - _c ) / aol::Sqr ( v );
  }

  RealType evaluateVVDerivative( const RealType u, const RealType v ) const{
    return _penaltyFunc.evaluateSecondDerivative ( u / v - _c ) * ( aol::Sqr ( u ) / aol::Qrt ( v ) )
           + 2 * _penaltyFunc.evaluateDerivative ( u / v - _c ) * ( u / aol::Cub ( v ) );
  }
};

/**
 * \brief calculates \f$ \int_\Omega \nu(u(x)) dx \f$, where \f$ u \f$ is given by Arg and \f$ \nu \f$ by PenaltyFunction.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename PenaltyFunctionType>
class PointwiseConstraintPenaltyEnergy
  : public FENonlinIntegrationScalarInterface<ConfiguratorType, PointwiseConstraintPenaltyEnergy<ConfiguratorType, PenaltyFunctionType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const PenaltyFunctionType &_penaltyFunction;
public:
  PointwiseConstraintPenaltyEnergy ( const typename ConfiguratorType::InitType &Initializer,
                                     const PenaltyFunctionType &PenaltyFunction )
    : FENonlinIntegrationScalarInterface<ConfiguratorType, PointwiseConstraintPenaltyEnergy<ConfiguratorType, PenaltyFunctionType> > ( Initializer ),
      _penaltyFunction( PenaltyFunction ) { }

  RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return _penaltyFunction.evaluate ( DiscFunc.evaluateAtQuadPoint( El, QuadPoint ) );
  }
};

/**
 * \brief calculates \f$ \int_\Omega \nu^\prime(u(x)) \varphi_i dx \f$, where \f$ u \f$ is given by Arg and \f$ \nu \f$ by PenaltyFunction.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename PenaltyFunctionType>
class VariationOfPointwiseConstraintPenaltyEnergy
  : public aol::FENonlinOpInterface< ConfiguratorType, VariationOfPointwiseConstraintPenaltyEnergy<ConfiguratorType, PenaltyFunctionType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const PenaltyFunctionType &_penaltyFunction;
public:
  VariationOfPointwiseConstraintPenaltyEnergy( const typename ConfiguratorType::InitType &Initializer,
                                               const PenaltyFunctionType &PenaltyFunction )
    : aol::FENonlinOpInterface< ConfiguratorType,
                                VariationOfPointwiseConstraintPenaltyEnergy<ConfiguratorType, PenaltyFunctionType> > ( Initializer ),
      _penaltyFunction( PenaltyFunction ) {
  }

  void getNonlinearity( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
                        typename ConfiguratorType::RealType &NL ) const {
    NL = _penaltyFunction.evaluateDerivative ( DiscFunc.evaluateAtQuadPoint( El, QuadPoint ) );
  }
};

/**
 * \brief calculates \f$ \int_\Omega \nu^{\prime\prime}(u(x)) \varphi_i(x) \varphi_j(x) dx \f$, where \f$ u \f$ is given by UDofs and \f$ \nu \f$ by PenaltyFunction.
 *
 * \author Berkels
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, typename PenaltyFunctionType>
class SecondVariationOfPointwiseConstraintPenaltyEnergy : public aol::FELinScalarWeightedMassInterface<ConfiguratorType, SecondVariationOfPointwiseConstraintPenaltyEnergy<ConfiguratorType, PenaltyFunctionType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrU;
  const PenaltyFunctionType &_penaltyFunction;
public:
  SecondVariationOfPointwiseConstraintPenaltyEnergy( const typename ConfiguratorType::InitType &Initializer,
                                                     const aol::Vector<RealType> &UDofs,
                                                     const PenaltyFunctionType &PenaltyFunction,
                                                     aol::OperatorType OpType = aol::ONTHEFLY )
    : aol::FELinScalarWeightedMassInterface<ConfiguratorType, SecondVariationOfPointwiseConstraintPenaltyEnergy<ConfiguratorType, PenaltyFunctionType> >( Initializer, OpType ),
      _discrU( Initializer, UDofs ),
      _penaltyFunction( PenaltyFunction )
  {
  }
  inline RealType getCoeff( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::VecType& ) const {
    return _penaltyFunction.evaluateSecondDerivative ( _discrU.evaluateAtQuadPoint ( El, QuadPoint ) );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename PenaltyFunctionType, typename MatrixType = qc::FastUniformGridMatrix<typename ConfiguratorType::RealType,ConfiguratorType::Dim> >
class SecondVariationOfEsedogluSegmentationModelWithPenalty : public aol::Op<aol::Vector<typename ConfiguratorType::RealType>, MatrixType > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const PenaltyFunctionType &_penaltyFunction;
  const RealType _lengthFactor, _penaltyFactor, _epsilon;
public:
  SecondVariationOfEsedogluSegmentationModelWithPenalty ( const typename ConfiguratorType::InitType &Initializer,
                                                          const PenaltyFunctionType &PenaltyFunction,
                                                          const RealType LengthFactor,
                                                          const RealType PenaltyFactor,
                                                          const RealType Epsilon )
  : _grid ( Initializer ),
    _penaltyFunction ( PenaltyFunction ),
    _lengthFactor ( LengthFactor ),
    _penaltyFactor ( PenaltyFactor ),
    _epsilon ( Epsilon )
  {
  }
  virtual void apply( const aol::Vector<RealType> &Arg, MatrixType &Dest ) const {
    Dest.setZero();
    aol::SecondVariationOfPointwiseConstraintPenaltyEnergy<ConfiguratorType, PenaltyFunctionType> DDEPenalty ( _grid, Arg, _penaltyFunction );
    DDEPenalty.assembleAddMatrix ( Dest, _penaltyFactor );

    aol::IsotropicMatrixStiffOp<ConfiguratorType> isotropicMatrixStiffOp ( _grid, Arg, _epsilon );
    isotropicMatrixStiffOp.assembleAddMatrix( Dest,_lengthFactor );
  }
  virtual void applyAdd( const aol::Vector<RealType> &/*Arg*/, MatrixType &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__);
  }
};

/**
 * calculates \f$ \sum_j\int_\Omega \nu(u_j(x)) dx \f$, where \f$ u_j \f$ is given by MArg and \f$ \nu \f$ by PenaltyFunction.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename PenaltyFunctionType, int numberOfFunctions>
class PointwiseConstraintPenaltyEnergyMulti
  : public aol::FENonlinIntegrationVectorInterface< ConfiguratorType,
                                                    PointwiseConstraintPenaltyEnergyMulti<ConfiguratorType, PenaltyFunctionType, numberOfFunctions>,
                                                    numberOfFunctions > {
  typedef typename ConfiguratorType::RealType RealType;
  const PenaltyFunctionType &_penaltyFunction;
public:
  PointwiseConstraintPenaltyEnergyMulti ( const typename ConfiguratorType::InitType &Initializer,
                                          const PenaltyFunctionType &PenaltyFunction )
    : aol::FENonlinIntegrationVectorInterface< ConfiguratorType,
                                               PointwiseConstraintPenaltyEnergyMulti<ConfiguratorType, PenaltyFunctionType, numberOfFunctions>,
                                               numberOfFunctions > ( Initializer ),
      _penaltyFunction( PenaltyFunction ) {
  }

  RealType evaluateIntegrand( const aol::auto_container<numberOfFunctions,aol::DiscreteFunctionDefault<ConfiguratorType> > &discrFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/ ) const {
    RealType integrand = 0.;
    for( int i = 0; i < numberOfFunctions; i++ ){
      integrand += _penaltyFunction.evaluate ( discrFuncs[i].evaluateAtQuadPoint( El, QuadPoint ) );
    }
    return integrand;
  }
};

/**
 * calculates \f$ \int_\Omega \nu^\prime(u_j(x)) \varphi_i dx \f$, where \f$ u_j \f$ is given by MArg and \f$ \nu \f$ by PenaltyFunction.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename PenaltyFunctionType, int numberOfFunctions>
class VariationOfPointwiseConstraintPenaltyEnergyMulti
: public Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const VariationOfPointwiseConstraintPenaltyEnergy<ConfiguratorType, PenaltyFunctionType> _DEPenaltySingle;
public:
  VariationOfPointwiseConstraintPenaltyEnergyMulti ( const typename ConfiguratorType::InitType &Initializer,
                                                     const PenaltyFunctionType &PenaltyFunction )
    : _DEPenaltySingle ( Initializer, PenaltyFunction )
  {
  }

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    for( int i = 0; i < numberOfFunctions; i++ ){
      _DEPenaltySingle.applyAdd( MArg[i], MDest[i] );
    }
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename PenaltyFunctionType, int numberOfFunctions, typename MatrixType = qc::FastUniformGridMatrix<typename ConfiguratorType::RealType,ConfiguratorType::Dim> >
class SecondVariationOfEsedogluSegmentationModelWithPenaltyMulti
: public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::SparseBlockMatrix<MatrixType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::SecondVariationOfEsedogluSegmentationModelWithPenalty<ConfiguratorType, PenaltyFunctionType, MatrixType> _DFSingle;
public:
  SecondVariationOfEsedogluSegmentationModelWithPenaltyMulti( const typename ConfiguratorType::InitType &Initializer,
                                                              const PenaltyFunctionType &PenaltyFunction,
                                                              const RealType LengthFactor,
                                                              const RealType PenaltyFactor,
                                                              const RealType Epsilon )
    : _grid ( Initializer ),
      _DFSingle ( Initializer, PenaltyFunction, LengthFactor, PenaltyFactor, Epsilon )
  {
  }

  virtual void apply( const aol::MultiVector<RealType> &MArg,  aol::SparseBlockMatrix<MatrixType> &MDest ) const{
    MDest.deleteBlockEntries();

    for( int i = 0; i < numberOfFunctions; i++ ){
      MatrixType &diagBlockEntryI = MDest.allocateMatrix( i, i, _grid );
      _DFSingle.apply( MArg[i], diagBlockEntryI );
    }
  }

  virtual void applyAdd( const aol::MultiVector<RealType> &/*MArg*/, aol::SparseBlockMatrix<MatrixType> &/*MDest*/ ) const{
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__);
  }
};

/**
 * \brief Iterates over all possible vectors of length NumberOfComponents with component values of zero or one.
 *
 * A vector of length l with component values of zero or one can be considered as an element of the power set
 * of a set with l elements, therefore the decision for the name of this class.
 *
 * \author Berkels
 */
class PowerSetIterator{
  const int _numberOfComponents;
  const int _numberOfModifiers;
  int _currentPostion;
public:
  PowerSetIterator( int NumberOfComponents )
    : _numberOfComponents( NumberOfComponents ),
      _numberOfModifiers( 1<<_numberOfComponents ),
      _currentPostion( 0 )
  {
  }
  //! Step to the next vector
  void increment(){
    if( _currentPostion < _numberOfModifiers )
      _currentPostion++;
  }
  //! Returns true, if the iterations over all possible vectors is finished.
  bool end() const{
    return ( _currentPostion == _numberOfModifiers );
  }
  //! returns the I-th component of the current vector
  int getComponent( int I ) const{
    // 2^I is equal to the bit shift 1<<I
    // Division of an integer by 2^I is equal to the bit shift >>I
    return ((_currentPostion & 1<<I )>>I);
  }
  int getCurrentPosition() const{
    return _currentPostion;
  }
  /**
   * Returns the position number given an element of the power set, i.e.
   * an int vector with component values of zero or one.
   *
   * This is needed for example, if you want to know which parameter is
   * active at a certain position in the computational domain from the
   * signs of the levelset functions.
   */
  static int getPositionNumberFromVector( const aol::Vector<int> &Vec ){
    int position = 0;
    for( int i = 0; i < Vec.size(); i++ )
      position += Vec[i]<<i;
    return position;
  }
  /**
   * Returns the position number given a global index and a set of
   * levelset functions.
   *
   * This is needed for example, if you want to know which parameter is
   * active at a certain position in the computational domain specified
   * by the global index based on the signs of the levelset functions.
   */
  template <typename RealType>
  static int getPositionNumberFromLevelsetFunctions( const aol::MultiVector<RealType> &LevelsetFunctions,
                                                     const int GlobalIndex,
                                                     const RealType Isovalue = aol::NumberTrait<RealType>::zero ){
    static aol::Vector<int> intVec;
    intVec.resize( LevelsetFunctions.numComponents() );
    for ( int i = 0; i < LevelsetFunctions.numComponents(); i++ )
      intVec[i] = (LevelsetFunctions[i].get(GlobalIndex)>= Isovalue) ? 1 : 0;
    return getPositionNumberFromVector(intVec);
  }
};

/**
 * \brief A helper class for Chan Vese with multiple levelset functions.
 *
 * Has to be constructed (or initialized with initValues) for each point where you want to evaluate
 * the product of the Heaviside funtions of the levelset functions.
 *
 * Note: In a proper Chan Vese model evaluateDerivative is not necessary,
 * but caching the necessary values for this function causes a noticeably
 * performance hit. Terefore all code for this is commented out at the
 * moment.
 *
 * If HeavisideFunctionType::ScalingFunctionType is equal to IdentityFunction,
 * this class fits to the classical Chan Vese model.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType>
class HeavisideFunctionProduct {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  typedef typename HeavisideFunctionType::ScalingFunctionType ScalingFunctionType;
  const HeavisideFunctionType &_heavisideFunction;
  const int _numLevelsetFunctions;
  aol::Vector<RealType> _hOfPhiValues;
  //aol::Vector<RealType> _hPrimeOfPhiValues;
  vector< aol::DiscreteFunctionDefault<ConfiguratorType>* > _discrLevesetFunctionsVec;

  //! May only be called once! Calling it more than once will lead to memory leaks.
  void initDiscrLevesetFunctionsVec ( const typename ConfiguratorType::InitType &Initializer,
                                      const aol::MultiVector<RealType> &LevelsetFunctions ) {
    _discrLevesetFunctionsVec.resize( _numLevelsetFunctions );
    for( int i = 0; i < _numLevelsetFunctions; i++ ){
      _discrLevesetFunctionsVec[i] = new aol::DiscreteFunctionDefault<ConfiguratorType> ( Initializer, LevelsetFunctions[i] );
    }
  }
public:
  void initValues ( const typename ConfiguratorType::ElementType &El,
                    const int QuadPoint ) {
    for( int i = 0; i < _numLevelsetFunctions; i++ ){
      const RealType phiValue = _discrLevesetFunctionsVec[i]->evaluateAtQuadPoint( El, QuadPoint);
      _hOfPhiValues[i] = _heavisideFunction.evaluate( phiValue );
      //_hPrimeOfPhiValues[i] = _heavisideFunction.evaluateDerivative( phiValue );
    }
  }
  HeavisideFunctionProduct( const typename ConfiguratorType::InitType &Initializer,
                            const HeavisideFunctionType &HeavisideFunction,
                            const aol::MultiVector<RealType> &LevelsetFunctions )
  : _heavisideFunction( HeavisideFunction ),
    _numLevelsetFunctions( LevelsetFunctions.numComponents() ),
    _hOfPhiValues( _numLevelsetFunctions )//,
    //_hPrimeOfPhiValues( _numLevelsetFunctions )
  {
    initDiscrLevesetFunctionsVec ( Initializer, LevelsetFunctions );
  }
  ~HeavisideFunctionProduct () {
    for( unsigned int i = 0; i < _discrLevesetFunctionsVec.size(); i++ ){
      delete _discrLevesetFunctionsVec[i];
    }
  }
  //! Note: Calling this constructor on every quad point will be slow!
  HeavisideFunctionProduct( const typename ConfiguratorType::InitType &Initializer,
                            const HeavisideFunctionType &HeavisideFunction,
                            const aol::MultiVector<RealType> &LevelsetFunctions,
                            const typename ConfiguratorType::ElementType &El,
                            const int QuadPoint )
  : _heavisideFunction( HeavisideFunction ),
    _numLevelsetFunctions( LevelsetFunctions.numComponents() ),
    _hOfPhiValues( _numLevelsetFunctions )//,
    //_hPrimeOfPhiValues( _numLevelsetFunctions )
  {
    initDiscrLevesetFunctionsVec ( Initializer, LevelsetFunctions );
    initValues ( El, QuadPoint );
  }
  //! Note: Calling this constructor on every quad point will be slow!
  HeavisideFunctionProduct( const typename ConfiguratorType::InitType &Initializer,
                            const HeavisideFunctionType &HeavisideFunction,
                            const aol::MultiVector<RealType> &LevelsetFunctions,
                            const typename ConfiguratorType::ElementType &El,
                            const typename ConfiguratorType::VecType& RefCoord )
  : _heavisideFunction( HeavisideFunction ),
    _numLevelsetFunctions( LevelsetFunctions.numComponents() ),
    _hOfPhiValues( _numLevelsetFunctions )//,
    //_hPrimeOfPhiValues( _numLevelsetFunctions )
  {
    for( int i = 0; i < _numLevelsetFunctions; i++ ){
      aol::DiscreteFunctionDefault<ConfiguratorType> discrFunc( Initializer, LevelsetFunctions[i] );
      const RealType phiValue = discrFunc.evaluate( El, RefCoord );
      _hOfPhiValues[i] = _heavisideFunction.evaluate( phiValue );
      //_hPrimeOfPhiValues[i] = _heavisideFunction.evaluateDerivative( phiValue );
    }
  }
  /**
   * Helper function to decide whether to use "x" or "1-x" without using an if statement.
   * f(x,1) = x
   * f(x,0) = 1-x
   */
  inline RealType f( const RealType x, const int b ) const {
    return ( 2*(b-0.5f)*(x - 0.5f) + 0.5f );
  }
  /**
   * Helper function to decide whether to use "x" or "-x" without using an if statement.
   * g(x,1) = x
   * g(x,0) = -x
   */
  inline RealType g( const RealType x, const int b ) const {
    return ( 2*(b-0.5f)*x );
  }
  RealType evaluate( const aol::PowerSetIterator &Iterator ) const {
    RealType value = aol::ZOTrait<RealType>::one;
    for( int i = 0; i < _numLevelsetFunctions; i++ ){
      value *= ScalingFunctionType::evaluate ( f( _hOfPhiValues[i], Iterator.getComponent(i) ) );
    }
    return value;
  }
  RealType evaluateSkippingOneComponent( const int ComponentToSkip, const aol::PowerSetIterator &Iterator ) const {
    RealType value = aol::ZOTrait<RealType>::one;

    for( int i = 0; i < ComponentToSkip; i++ ){
      value *= ScalingFunctionType::evaluate ( f( _hOfPhiValues[i], Iterator.getComponent(i) ) );
    }
    for( int i = ComponentToSkip+1; i < _numLevelsetFunctions; i++ ){
      value *= ScalingFunctionType::evaluate ( f( _hOfPhiValues[i], Iterator.getComponent(i) ) );
    }
    return value;
  }
  /**
   * ATTENTION: If ScalingFunctionType::evaluateDerivative doesn't always return one,
   * evaluateSkippingOneComponentWithSign doesn't exactly do what the name suggests.
   */
  RealType evaluateSkippingOneComponentWithSign( const int ComponentToSkip, const aol::PowerSetIterator &Iterator ) const {
    const RealType sign = 2*(Iterator.getComponent(ComponentToSkip)-0.5f);
    return (sign * ScalingFunctionType::evaluateDerivative ( f( _hOfPhiValues[ComponentToSkip], Iterator.getComponent(ComponentToSkip) ) ) * evaluateSkippingOneComponent( ComponentToSkip, Iterator ));
  }
  RealType evaluateSkippingTwoComponents( const int ComponentToSkipOne, const int ComponentToSkipTwo, const aol::PowerSetIterator &Iterator ) const {
    RealType value = aol::ZOTrait<RealType>::one;

    for( int i = 0; i < _numLevelsetFunctions; i++ ){
      if ( ( i != ComponentToSkipOne ) && ( i != ComponentToSkipTwo ) )
        value *= ScalingFunctionType::evaluate ( f( _hOfPhiValues[i], Iterator.getComponent(i) ) );
    }
    return value;
  }
  RealType evaluateMixedSecondDerivativeForCEMModel( const int ComponentOne, const int ComponentTwo, const aol::PowerSetIterator &Iterator ) const {
    const RealType sign1 = 2*(Iterator.getComponent(ComponentOne)-0.5f);
    const RealType sign2 = 2*(Iterator.getComponent(ComponentTwo)-0.5f);
    return ( sign1 * ScalingFunctionType::evaluateDerivative ( f( _hOfPhiValues[ComponentOne], Iterator.getComponent(ComponentOne) ) )
             * sign2 * ScalingFunctionType::evaluateDerivative ( f( _hOfPhiValues[ComponentTwo], Iterator.getComponent(ComponentTwo) ) )
             * evaluateSkippingTwoComponents( ComponentOne, ComponentTwo, Iterator ));
  }
  /*
  RealType evaluateDerivative( const int ComponentToDerive, const aol::PowerSetIterator &Iterator ) const {
    const RealType value = g( _hPrimeOfPhiValues[ComponentToDerive], Iterator.getComponent(ComponentToDerive) );
    return (value * evaluateSkippingOneComponent( ComponentToDerive, Iterator ));
  }
  */
};

/**
 * \brief Given n levelset functions and 2^n ops, this class calculates the sum
 *        over these ops, where the integrands of the ops are multiplied with
 *        the proper combination of heaviside functions of the levelset
 *        functions as needed in the Chan Vese model for more than two segments.
 *        The ops have to be derived from FENonlinIntegrationVectorInterface.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, typename TargetOpType>
class HeavisideProductFENonlinIntegrationVectorInterface
: public aol::FENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideProductFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType>, TargetOpType::NumOfComponents > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::MultiVector<RealType> &_levelsetFunctions;
  const HeavisideFunctionType &_heavisideFunction;
  const std::vector<const TargetOpType*> _targetOpVec;
public:
  HeavisideProductFENonlinIntegrationVectorInterface( const typename ConfiguratorType::InitType &Initializer,
                                                        const aol::MultiVector<RealType> &LevelsetFunctions,
                                                        const HeavisideFunctionType &HeavisideFunction,
                                                        const std::vector<const TargetOpType*> &TargetOpVec )
    : aol::FENonlinIntegrationVectorInterface< ConfiguratorType,
                                                 HeavisideProductFENonlinIntegrationVectorInterface<ConfiguratorType, HeavisideFunctionType, TargetOpType>,
                                                 TargetOpType::NumOfComponents > ( Initializer ),
      _levelsetFunctions( LevelsetFunctions ),
      _heavisideFunction( HeavisideFunction ),
      _targetOpVec( TargetOpVec ){
    // Given n levelset functions, we have 2^n possible segments and also need that many target ops.
    if ( 1<<_levelsetFunctions.numComponents() != static_cast<int>(_targetOpVec.size()) )
      throw aol::Exception( "pow( 2, _levelsetFunctions.numComponents()) != _targetOpVec.length()", __FILE__, __LINE__ );
  }

  RealType evaluateIntegrand( const aol::auto_container<TargetOpType::NumOfComponents,aol::DiscreteFunctionDefault<ConfiguratorType> > &discrFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint, const typename ConfiguratorType::VecType &RefCoord ) const {
    // Construct HeavisideFunctionProduct at the point, where we want to evaluate it.
    HeavisideFunctionProduct<ConfiguratorType, HeavisideFunctionType>
      heavisideFunctionProduct( this->_initializer, _heavisideFunction, _levelsetFunctions, El, QuadPoint );

    RealType integrand = aol::ZOTrait<RealType>::zero;
    // Do the summation over all combinations of Levelset function signs.
    aol::PowerSetIterator iterator(_levelsetFunctions.numComponents());
    do
    {
      integrand += ( _targetOpVec[iterator.getCurrentPosition()]->evaluateIntegrand( discrFuncs, El, QuadPoint, RefCoord )
                     * heavisideFunctionProduct.evaluate( iterator ) );
      iterator.increment();
    } while ( !iterator.end() );

    return integrand;
  }
};

/**
 * Checks if the MultiVector is of the correct structure to be an
 * argument of one of the multi levelset Chan Vese operators.
 * Throws an expection if the check fails.
 *
 * \author Berkels
 */
template <typename RealType>
void checkMultiChanVeseStructure( const aol::MultiVector<RealType> &LevelsetfunctionsAndParameters, const int NumberOfLevelsetFunctions, const int NumberOfAdditionalVectors = 0 ){
  if ( LevelsetfunctionsAndParameters.numComponents() != NumberOfLevelsetFunctions + (1<<NumberOfLevelsetFunctions) + NumberOfAdditionalVectors )
    throw aol::Exception( "LevelsetfunctionsAndParameters.numComponents() != NumberOfLevelsetFunctions + 1<<numberOfLevelsetFunctions", __FILE__, __LINE__ );
}

/**
 * General Interface for multi levelset Chan Vese type energies.
 * Only evaluateIndicator has to be provided. For optimizations
 * cacheIndicatorParameters can be overloaded.
 *
 * Apply expects MArg to contain the n levelset functions in the
 * first n components and the indicator paramters in the next
 * 2^n components.
 *
 * If MArg contains more than the above mentioned n+2^n components,
 * _numOfAdditionalVectorComponentsInArgument has to be set
 * accordingly.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int ParameterDimension, typename Imp>
class ChanVeseEnergyInterface
  : public Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  mutable ConfiguratorType _config;
  const typename ConfiguratorType::InitType &_initializer;
  const HeavisideFunctionType &_heavisideFunction;
  int _numOfAdditionalVectorComponentsInArgument;
public:
  ChanVeseEnergyInterface ( const typename ConfiguratorType::InitType &Initializer,
                            const HeavisideFunctionType &HeavisideFunction ) :
    _config ( Initializer ),
    _initializer ( Initializer ),
    _heavisideFunction( HeavisideFunction ),
    _numOfAdditionalVectorComponentsInArgument( 0 )
  {}

  virtual ~ChanVeseEnergyInterface() {}

  virtual void cacheIndicatorParameters( const aol::MultiVector<RealType> &/*Parameters*/ ) const {}

  void setNumOfAdditionalVectorComponentsInArgument ( const int NumOfAdditionalVectorComponentsInArgument ){
    _numOfAdditionalVectorComponentsInArgument = NumOfAdditionalVectorComponentsInArgument;
  }

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {

    checkMultiChanVeseStructure<RealType>(MArg, NumberOfLevelsetFunctions, _numOfAdditionalVectorComponentsInArgument);

    // Group the levelset functions contained in MArg into a new MultiVector.
    aol::MultiVector<RealType> levelsetFunctions( 0, 0 );
    for( int i = 0; i < NumberOfLevelsetFunctions; i++ )
      levelsetFunctions.appendReference( MArg[i] );

    // Group the parameter vectors contained in MArg into a new MultiVector.
    aol::MultiVector<RealType> parameters( 0, 0 );
    for( int i = 0; i < (1<<NumberOfLevelsetFunctions); i++ )
      parameters.appendReference( MArg[i+NumberOfLevelsetFunctions] );

    // Since the indicator parameters are fixed during this calculation
    // it may be possible to speed up evaluateIndicatorDerivative, by
    // pre-calculating all dependencies on indicator parameters here.
    cacheIndicatorParameters( parameters );

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    RealType res = aol::ZOTrait<RealType>::zero;

    HeavisideFunctionProduct<ConfiguratorType, HeavisideFunctionType>
      heavisideFunctionProduct( _initializer, _heavisideFunction, levelsetFunctions );

    for ( IteratorType it = _config.begin(); it != _config.end(); ++it ) {
      typedef typename ConfiguratorType::QuadType QType;
      RealType a = aol::ZOTrait<RealType>::zero;
      for ( int q = 0; q < QType::numQuadPoints; q++ ) {

        // Init HeavisideFunctionProduct at the point, where we want to evaluate it.
        heavisideFunctionProduct.initValues ( *it, q );

        // Loop over all combinations of Levelset function signs, each represents one of the indicator parameters.
        aol::PowerSetIterator iterator( levelsetFunctions.numComponents() );
        const RealType weight = _config.getBaseFunctionSet ( *it ).getWeight ( q );
        const typename ConfiguratorType::VecType &refCoord = _config.getBaseFunctionSet ( *it ).getRefCoord ( q );
        do
        {
          const int i = iterator.getCurrentPosition();
          a += ( this->asImp().evaluateIndicator( MArg[NumberOfLevelsetFunctions+i], i, *it, q, refCoord )
                 * heavisideFunctionProduct.evaluate( iterator ) )
                 * weight;
          iterator.increment();
        } while ( !iterator.end() );

      }

      a *= _config.vol ( *it );

      res += a;
    }
    Dest += res;
  }

  //! interface function, has to be provided in derived classes.
  RealType evaluateIndicator ( const aol::Vector<RealType> &IndicatorParameter,
                               const int IndicatorParameterIndex,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint,
                               const typename ConfiguratorType::VecType &RefCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().evaluateIndicator ( IndicatorParameter, IndicatorParameterIndex, El, QuadPoint, RefCoord );
  }

protected:
  // barton-nackman
  inline Imp& asImp() {
    return static_cast<Imp&> ( *this );
  }
  inline const Imp& asImp() const {
    return static_cast<const Imp&> ( *this );
  }

};


/**
 * General Interface for the derivate of multi levelset Chan Vese type
 * energies with respect to the indicator parameters. Only
 * evaluateIndicatorDerivative has to be provided. For optimizations
 * cacheIndicatorParameters can be overloaded.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int ParameterDimension, typename Imp>
class ChanVeseEnergyParameterDerivativeInterface
      : public Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  mutable ConfiguratorType _config;
  const typename ConfiguratorType::InitType &_initializer;
  const aol::MultiVector<RealType> &_levelsetFunctions;
  const HeavisideFunctionType &_heavisideFunction;
public:
  ChanVeseEnergyParameterDerivativeInterface ( const typename ConfiguratorType::InitType &Initializer,
                                               const aol::MultiVector<RealType> &LevelsetFunctions,
                                               const HeavisideFunctionType &HeavisideFunction ) :
    _config ( Initializer ),
    _initializer ( Initializer ),
    _levelsetFunctions( LevelsetFunctions ),
    _heavisideFunction( HeavisideFunction )
  {}

  virtual ~ChanVeseEnergyParameterDerivativeInterface() {}

  virtual void cacheIndicatorParameters( const aol::MultiVector<RealType> &/*Parameters*/ ) const {}

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {

    // Since the indicator parameters are fixed during this calculation
    // it may be possible to speed up evaluateIndicatorDerivative, by
    // pre-calculating all dependencies on indicator parameters here.
    cacheIndicatorParameters( Arg );

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;

    const int NumberOfSegments = 1<<NumberOfLevelsetFunctions;

    aol::Mat<NumberOfSegments, ParameterDimension, RealType> res;
    aol::Mat<NumberOfSegments, ParameterDimension, RealType> a;
    aol::Vec<ParameterDimension, RealType> tempVec;

    HeavisideFunctionProduct<ConfiguratorType, HeavisideFunctionType>
      heavisideFunctionProduct( _initializer, _heavisideFunction, _levelsetFunctions );

    for ( IteratorType it = _config.begin(); it != _config.end(); ++it ) {
      typedef typename ConfiguratorType::QuadType QType;
      a.setZero();
      for ( int q = 0; q < QType::numQuadPoints; q++ ) {

        // Init HeavisideFunctionProduct at the point, where we want to evaluate it.
        heavisideFunctionProduct.initValues( *it, q );

        // Loop over all combinations of Levelset function signs, each represents one of the indicator parameters.
        aol::PowerSetIterator iterator(_levelsetFunctions.numComponents());
        const RealType weight = _config.getBaseFunctionSet ( *it ).getWeight ( q );
        do
        {
          const int i = iterator.getCurrentPosition();
          this->asImp().evaluateIndicatorDerivative( Arg[i], i, *it, q, _config.getBaseFunctionSet ( *it ).getRefCoord ( q ), tempVec );
          tempVec *= heavisideFunctionProduct.evaluate( iterator ) * weight;
          a[i] += tempVec;
          iterator.increment();
        } while ( !iterator.end() );

      }

      a *= _config.vol ( *it );

      res += a;
    }
    for ( int i = 0; i < NumberOfSegments; i++ ) {
      for ( int j = 0; j < ParameterDimension; j++ ) {
        Dest[i][j] += res[i][j];
      }
    }
  }

  //! interface function, has to be provided in derived classes.
  void evaluateIndicatorDerivative ( const aol::Vector<RealType> &IndicatorParameter,
                                     const int IndicatorParameterIndex,
                                     const typename ConfiguratorType::ElementType &El,
                                     int QuadPoint,
                                     const typename ConfiguratorType::VecType &RefCoord,
                                     aol::Vec<ParameterDimension, RealType> &Derivative) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    this->asImp().evaluateIndicatorDerivative ( IndicatorParameter, IndicatorParameterIndex, El, QuadPoint, RefCoord, Derivative );
  }

protected:
  // barton-nackman
  inline Imp& asImp() {
    return static_cast<Imp&> ( *this );
  }
  inline const Imp& asImp() const {
    return static_cast<const Imp&> ( *this );
  }

};


/**
 * General Interface for the derivate of multi levelset Chan Vese type
 * energies with respect to the levelset functions. Only
 * evaluateIndicator has to be provided. For optimizations
 * cacheIndicatorParameters can be overloaded.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, typename Imp>
class ChanVeseEnergyLevelsetDerivativeInterface
      : public Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  mutable ConfiguratorType _config;
  const typename ConfiguratorType::InitType &_initializer;
  const aol::MultiVector<RealType> &_parameters;
  const HeavisideFunctionType &_heavisideFunction;
public:
  ChanVeseEnergyLevelsetDerivativeInterface ( const typename ConfiguratorType::InitType &Initializer,
                                              const aol::MultiVector<RealType> &Parameters,
                                              const HeavisideFunctionType &HeavisideFunction ) :
    _config ( Initializer ),
    _initializer ( Initializer ),
    _parameters( Parameters ),
    _heavisideFunction( HeavisideFunction )
  {}

  virtual ~ChanVeseEnergyLevelsetDerivativeInterface() {}

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    aol::Vec<NumberOfLevelsetFunctions, RealType> a;
    aol::FullMatrix<RealType> indicatorCache(_config.maxNumQuadPoints(), _parameters.numComponents() );

    std::vector<HeavisideFunctionProduct<ConfiguratorType, HeavisideFunctionType>*> heavisideFunctionProductCache(_config.maxNumQuadPoints());
    for ( int q = 0; q < _config.maxNumQuadPoints(); q++ ) {
      heavisideFunctionProductCache[q] = new HeavisideFunctionProduct<ConfiguratorType, HeavisideFunctionType>( _initializer, _heavisideFunction, MArg );
    }

    for ( IteratorType it = _config.begin(); it != _config.end(); ++it ) {
      a.setZero();
      const int numLocalDofs = _config.getNumLocalDofs ( *it );

      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      // Cache the indicator values at all quad points. We do this outside the "dof" loop
      // since the values don't depend on dof.
      // The extra scope deletes iterator after the iteration.
      {
        aol::PowerSetIterator iterator( NumberOfLevelsetFunctions );
        do
        {
          const int i = iterator.getCurrentPosition();
          for ( int q = 0; q < numQuadPoints; q++ ) {
            indicatorCache.set(q,i,this->asImp().evaluateIndicator( _parameters[i], i, *it, q, _config.getBaseFunctionSet ( *it ).getRefCoord ( q ) ));
          }
          iterator.increment();
        } while ( !iterator.end() );
      }

      // Init HeavisideFunctionProduct at each quad point to cache it's values there.
      // We do this outside the "dof" loop since the values don't depend on dof.
      for ( int q = 0; q < numQuadPoints; q++ ) {
        heavisideFunctionProductCache[q]->initValues ( *it, q );
      }

      for ( int dof = 0; dof < numLocalDofs; dof++ ) {

        for ( int q = 0; q < numQuadPoints; q++ ) {

          // Loop over all combinations of Levelset function signs, each represents one of the indicator parameters.
          aol::PowerSetIterator iterator( NumberOfLevelsetFunctions );
          const RealType weight = bfs.getWeight ( q );
          const RealType basefunctionValue = bfs.evaluate ( dof, q );
          do
          {
            const int i = iterator.getCurrentPosition();
            const RealType temp = indicatorCache.get(q,i)*weight*basefunctionValue;
            for( int j = 0; j < NumberOfLevelsetFunctions; j++ )
              a[j] += temp * heavisideFunctionProductCache[q]->evaluateSkippingOneComponentWithSign( j, iterator );
            // Perhaps evaluateSkippingOneComponentWithSign should be cached somehow?
            iterator.increment();
          } while ( !iterator.end() );

        }

        a *= _config.vol ( *it );

        const int globalIndex = _config.localToGlobal ( *it, dof );
        for( int j = 0; j < NumberOfLevelsetFunctions; j++ ) {
          MDest[j][globalIndex] += a[j];
        }
      }
    }
    for ( int q = 0; q < _config.maxNumQuadPoints(); q++ ) {
      delete heavisideFunctionProductCache[q];
    }
  }

  //! interface function, has to be provided in derived classes.
  RealType evaluateIndicator ( const aol::Vector<RealType> &IndicatorParameter,
                               const int IndicatorParameterIndex,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint,
                               const typename ConfiguratorType::VecType &RefCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().evaluateIndicator ( IndicatorParameter, IndicatorParameterIndex, El, QuadPoint, RefCoord );
  }

protected:
  // barton-nackman
  inline Imp& asImp() {
    return static_cast<Imp&> ( *this );
  }
  inline const Imp& asImp() const {
    return static_cast<const Imp&> ( *this );
  }

};


/**
 * General Interface for the second derivate of multi levelset
 * quadratic Chan Esedoglu Nikolova type energies with respect
 * to one of the levelset functions. Only evaluateIndicator has
 * to be provided. For optimizations cacheIndicatorParameters
 * can be overloaded.
 *
 * Assumes that LevelsetFunctions and Parameters don't change
 * after this class is created, because cacheIndicatorParameters
 * is only called once in the constructor and nowhere else.
 *
 * Only meant to be used with
 * HeavisideFuncType == aol::IdentityFunction<typename ConfiguratorType::RealType>
 * Probably the template argument should be removed completely.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, typename Imp>
class QuadraticCEMEnergySecondLevelsetDerivativeInterface
  : public aol::FELinScalarWeightedMassInterface<ConfiguratorType, QuadraticCEMEnergySecondLevelsetDerivativeInterface<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, Imp> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::MultiVector<RealType> &_levelsetFunctions;
  const aol::MultiVector<RealType> &_parameters;
  const HeavisideFunctionType &_heavisideFunction;
  mutable HeavisideFunctionProduct<ConfiguratorType, HeavisideFunctionType> _heavisideFunctionProduct;
  const int _levelsetFunctionNumber;
public:
  QuadraticCEMEnergySecondLevelsetDerivativeInterface ( const typename ConfiguratorType::InitType &Initializer,
                                                        const aol::MultiVector<RealType> &LevelsetFunctions,
                                                        const aol::MultiVector<RealType> &Parameters,
                                                        const HeavisideFunctionType &HeavisideFunction,
                                                        const int LevelsetFunctionNumber,
                                                        aol::OperatorType OpType = ONTHEFLY )
  : aol::FELinScalarWeightedMassInterface<ConfiguratorType, QuadraticCEMEnergySecondLevelsetDerivativeInterface<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, Imp> > ( Initializer, OpType ),
    _levelsetFunctions ( LevelsetFunctions ),
    _parameters( Parameters ),
    _heavisideFunction( HeavisideFunction ),
    _heavisideFunctionProduct( Initializer, HeavisideFunction, LevelsetFunctions ),
    _levelsetFunctionNumber ( LevelsetFunctionNumber )
  {
    if ( HeavisideFunctionType::ScalingFunctionType::evaluate != aol::SquareFunction<RealType>::evaluate )
      throw aol::Exception( "QuadraticCEMEnergySecondLevelsetDerivativeInterface is only meant to be used with HeavisideFuncType == aol::IdentityFunction<typename ConfiguratorType::RealType> ", __FILE__, __LINE__ );

    // Since the indicator parameters are assumed to be fixed here
    // it may be possible to speed up evaluateIndicatorDerivative, by
    // pre-calculating all dependencies on indicator parameters here.
    cacheIndicatorParameters( Parameters );
  }

  virtual ~QuadraticCEMEnergySecondLevelsetDerivativeInterface() {}

  virtual void cacheIndicatorParameters( const aol::MultiVector<RealType> &/*Parameters*/ ) const {}

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType& RefCoord ) const {

    _heavisideFunctionProduct.initValues ( El, QuadPoint );

    RealType a = aol::ZOTrait<RealType>::zero;
    // Loop over all combinations of Levelset function signs, each represents one of the indicator parameters.
    aol::PowerSetIterator iterator( NumberOfLevelsetFunctions );
    do
    {
      const int i = iterator.getCurrentPosition();
      a += ( this->asImp().evaluateIndicator( _parameters[i], i, El, QuadPoint, RefCoord )
        * _heavisideFunctionProduct.evaluateSkippingOneComponent( _levelsetFunctionNumber, iterator ) )
        * 2;
      iterator.increment();
    } while ( !iterator.end() );

    return a;
  }

  //! interface function, has to be provided in derived classes.
  RealType evaluateIndicator ( const aol::Vector<RealType> &IndicatorParameter,
                               const int IndicatorParameterIndex,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint,
                               const typename ConfiguratorType::VecType &RefCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().evaluateIndicator ( IndicatorParameter, IndicatorParameterIndex, El, QuadPoint, RefCoord );
  }

protected:
  // barton-nackman
  inline Imp& asImp() {
    return static_cast<Imp&> ( *this );
  }
  inline const Imp& asImp() const {
    return static_cast<const Imp&> ( *this );
  }

};

/**
 * General Interface for the mixed second derivates of multi levelset
 * quadratic Chan Esedoglu Nikolova type energies with respect
 * to two different of the levelset functions. Only evaluateIndicator has
 * to be provided. For optimizations cacheIndicatorParameters
 * can be overloaded.
 *
 * Assumes that LevelsetFunctions and Parameters don't change
 * after this class is created, because cacheIndicatorParameters
 * is only called once in the constructor and nowhere else.
 *
 * Only meant to be used with
 * HeavisideFuncType == aol::IdentityFunction<typename ConfiguratorType::RealType>
 * Probably the template argument should be removed completely.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, typename Imp>
class QuadraticCEMEnergyMixedSecondLevelsetDerivativeInterface
  : public aol::FELinScalarWeightedMassInterface<ConfiguratorType, QuadraticCEMEnergyMixedSecondLevelsetDerivativeInterface<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, Imp> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::MultiVector<RealType> &_levelsetFunctions;
  const aol::MultiVector<RealType> &_parameters;
  const HeavisideFunctionType &_heavisideFunction;
  mutable HeavisideFunctionProduct<ConfiguratorType, HeavisideFunctionType> _heavisideFunctionProduct;
  const int _levelsetFunctionNumberOne, _levelsetFunctionNumberTwo;
public:
  QuadraticCEMEnergyMixedSecondLevelsetDerivativeInterface ( const typename ConfiguratorType::InitType &Initializer,
                                                             const aol::MultiVector<RealType> &LevelsetFunctions,
                                                             const aol::MultiVector<RealType> &Parameters,
                                                             const HeavisideFunctionType &HeavisideFunction,
                                                             const int LevelsetFunctionNumberOne,
                                                             const int LevelsetFunctionNumberTwo,
                                                             aol::OperatorType OpType = ONTHEFLY )
  : aol::FELinScalarWeightedMassInterface<ConfiguratorType, QuadraticCEMEnergyMixedSecondLevelsetDerivativeInterface<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, Imp> > ( Initializer, OpType ),
    _levelsetFunctions ( LevelsetFunctions ),
    _parameters( Parameters ),
    _heavisideFunction( HeavisideFunction ),
    _heavisideFunctionProduct( Initializer, HeavisideFunction, LevelsetFunctions ),
    _levelsetFunctionNumberOne ( LevelsetFunctionNumberOne ),
    _levelsetFunctionNumberTwo ( LevelsetFunctionNumberTwo )
  {
    if ( HeavisideFunctionType::ScalingFunctionType::evaluate != aol::SquareFunction<RealType>::evaluate )
      throw aol::Exception( "QuadraticCEMEnergyMixedSecondLevelsetDerivativeInterface is only meant to be used with HeavisideFuncType == aol::IdentityFunction<typename ConfiguratorType::RealType> ", __FILE__, __LINE__ );

    // Since the indicator parameters are assumed to be fixed here
    // it may be possible to speed up evaluateIndicatorDerivative, by
    // pre-calculating all dependencies on indicator parameters here.
    cacheIndicatorParameters( Parameters );
  }

  virtual ~QuadraticCEMEnergyMixedSecondLevelsetDerivativeInterface() {}

  virtual void cacheIndicatorParameters( const aol::MultiVector<RealType> &/*Parameters*/ ) const {}

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType& RefCoord ) const {

    _heavisideFunctionProduct.initValues ( El, QuadPoint );

    RealType a = aol::ZOTrait<RealType>::zero;
    // Loop over all combinations of Levelset function signs, each represents one of the indicator parameters.
    aol::PowerSetIterator iterator( NumberOfLevelsetFunctions );
    do
    {
      const int i = iterator.getCurrentPosition();
      a += ( this->asImp().evaluateIndicator( _parameters[i], i, El, QuadPoint, RefCoord )
        * _heavisideFunctionProduct.evaluateMixedSecondDerivativeForCEMModel( _levelsetFunctionNumberOne, _levelsetFunctionNumberTwo, iterator ) );
      iterator.increment();
    } while ( !iterator.end() );

    return a;
  }

  //! interface function, has to be provided in derived classes.
  RealType evaluateIndicator ( const aol::Vector<RealType> &IndicatorParameter,
                               const int IndicatorParameterIndex,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint,
                               const typename ConfiguratorType::VecType &RefCoord ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().evaluateIndicator ( IndicatorParameter, IndicatorParameterIndex, El, QuadPoint, RefCoord );
  }

protected:
  // barton-nackman
  inline Imp& asImp() {
    return static_cast<Imp&> ( *this );
  }
  inline const Imp& asImp() const {
    return static_cast<const Imp&> ( *this );
  }

};

/**
 * Given a Chan Vese energy op and a parameter vector, this class
 * builds an op, which can be applied on MultiVectors of levelset
 * functions and calculates the energy for the levelset functions
 * with the given parameter vector.
 *
 * \author Berkels
 */
template <typename RealType>
class ChanVeseEnergyLevelsetPart : public aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > {
protected:
  const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &_chanVeseEnergyOp;
  const aol::MultiVector<RealType> &_parameters;
public:
  ChanVeseEnergyLevelsetPart ( const aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > &ChanVeseEnergyOp,
                               const aol::MultiVector<RealType> &Parameters )
    : _chanVeseEnergyOp ( ChanVeseEnergyOp ),
      _parameters ( Parameters ) {}
  virtual void apply ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    aol::MultiVector<RealType> mtemp( MArg, aol::FLAT_COPY );
    mtemp.appendReference( _parameters );
    _chanVeseEnergyOp.apply ( mtemp, Dest );
  }
  virtual void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    aol::MultiVector<RealType> mtemp( MArg, aol::FLAT_COPY );
    mtemp.appendReference( _parameters );
    _chanVeseEnergyOp.applyAdd ( mtemp, Dest );
  }
};

/**
 * Calulates the classical Chan Vese energy for multi levelset gray value segmentation
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions>
class ClassicalChanVeseEnergyMulti
: public aol::ChanVeseEnergyInterface< ConfiguratorType,
                                       HeavisideFunctionType,
                                       NumberOfLevelsetFunctions,
                                       1,
                                       ClassicalChanVeseEnergyMulti<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrU0;
public:
  ClassicalChanVeseEnergyMulti( const typename ConfiguratorType::InitType &Initializer,
                                const HeavisideFunctionType &HeavisideFunction,
                                const aol::Vector<RealType> &U0 )
  : aol::ChanVeseEnergyInterface< ConfiguratorType,
                                  HeavisideFunctionType,
                                  NumberOfLevelsetFunctions,
                                  1,
                                  ClassicalChanVeseEnergyMulti<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> >
                                ( Initializer, HeavisideFunction ),
      _discrU0( Initializer, U0 )
  {
  }
  ~ClassicalChanVeseEnergyMulti () {}

  RealType evaluateIndicator ( const aol::Vector<RealType> &IndicatorParameter,
                               const int /*IndicatorParameterIndex*/,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint,
                               const typename ConfiguratorType::VecType &/*RefCoord*/ ) const {
    return aol::Sqr(_discrU0.evaluateAtQuadPoint( El, QuadPoint) - IndicatorParameter[0] );
  }
};

/**
 * Calulates the variation of ClassicalChanVeseEnergyMulti with respect to the levelset functions.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions>
class ClassicalChanVeseEnergyMultiPhiVariation
: public aol::ChanVeseEnergyLevelsetDerivativeInterface< ConfiguratorType,
                                                         HeavisideFunctionType,
                                                         NumberOfLevelsetFunctions,
                                                         ClassicalChanVeseEnergyMultiPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrU0;
public:
  ClassicalChanVeseEnergyMultiPhiVariation( const typename ConfiguratorType::InitType &Initializer,
                                            const HeavisideFunctionType &HeavisideFunction,
                                            const aol::Vector<RealType> &U0,
                                            const aol::MultiVector<RealType> &Parameters )
  : aol::ChanVeseEnergyLevelsetDerivativeInterface< ConfiguratorType,
                                                    HeavisideFunctionType,
                                                    NumberOfLevelsetFunctions,
                                                    ClassicalChanVeseEnergyMultiPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> >
                                                   ( Initializer, Parameters, HeavisideFunction ),
      _discrU0( Initializer, U0 )
  {
  }
  ~ClassicalChanVeseEnergyMultiPhiVariation () {}

  RealType evaluateIndicator ( const aol::Vector<RealType> &IndicatorParameter,
                               const int /*IndicatorParameterIndex*/,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint,
                               const typename ConfiguratorType::VecType &/*RefCoord*/ ) const {
    return aol::Sqr(_discrU0.evaluateAtQuadPoint( El, QuadPoint) - IndicatorParameter[0] );
  }
};

/**
 * Calulates the classical Chan Vese energy of vector valued images for multi levelset gray value segmentation
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int ImageDimension>
class ClassicalChanVeseVectorEnergyMulti
: public aol::ChanVeseEnergyInterface< ConfiguratorType,
                                       HeavisideFunctionType,
                                       NumberOfLevelsetFunctions,
                                       ImageDimension,
                                       ClassicalChanVeseVectorEnergyMulti<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, ImageDimension> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  aol::auto_container< ImageDimension, const aol::DiscreteFunctionDefault<ConfiguratorType> > _discrU0;
public:
  ClassicalChanVeseVectorEnergyMulti( const typename ConfiguratorType::InitType &Initializer,
                                      const HeavisideFunctionType &HeavisideFunction,
                                      const aol::MultiVector<RealType> &U0 )
  : aol::ChanVeseEnergyInterface< ConfiguratorType,
                                  HeavisideFunctionType,
                                  NumberOfLevelsetFunctions,
                                  ImageDimension,
                                  ClassicalChanVeseVectorEnergyMulti<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, ImageDimension> >
                                ( Initializer, HeavisideFunction )
  {
    for ( int c=0; c<ImageDimension; c++ ) {
      aol::DiscreteFunctionDefault<ConfiguratorType> temp( Initializer, U0[c] );
      _discrU0.set_copy( c, temp );
    }
  }
  ~ClassicalChanVeseVectorEnergyMulti () {}

  RealType evaluateIndicator ( const aol::Vector<RealType> &IndicatorParameter,
                               const int /*IndicatorParameterIndex*/,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint,
                               const typename ConfiguratorType::VecType &/*RefCoord*/ ) const {
    RealType indicator = 0.;
    for ( int i = 0; i < ImageDimension; i++ )
      indicator += aol::Sqr(_discrU0[i].evaluateAtQuadPoint( El, QuadPoint) - IndicatorParameter[i] );
    return indicator;
  }
};

/**
 * Calulates the variation of ClassicalChanVeseVectorEnergyMulti with respect to the levelset functions.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int ImageDimension>
class ClassicalChanVeseVectorEnergyMultiPhiVariation
: public aol::ChanVeseEnergyLevelsetDerivativeInterface< ConfiguratorType,
                                                         HeavisideFunctionType,
                                                         NumberOfLevelsetFunctions,
                                                         ClassicalChanVeseVectorEnergyMultiPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, ImageDimension> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  aol::auto_container< ImageDimension, const aol::DiscreteFunctionDefault<ConfiguratorType> > _discrU0;
public:
  ClassicalChanVeseVectorEnergyMultiPhiVariation( const typename ConfiguratorType::InitType &Initializer,
                                            const HeavisideFunctionType &HeavisideFunction,
                                            const aol::MultiVector<RealType> &U0,
                                            const aol::MultiVector<RealType> &Parameters )
  : aol::ChanVeseEnergyLevelsetDerivativeInterface< ConfiguratorType,
                                                    HeavisideFunctionType,
                                                    NumberOfLevelsetFunctions,
                                                    ClassicalChanVeseVectorEnergyMultiPhiVariation<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, ImageDimension> >
                                                   ( Initializer, Parameters, HeavisideFunction )
  {
    for ( int c=0; c<ImageDimension; c++ ) {
      aol::DiscreteFunctionDefault<ConfiguratorType> temp( Initializer, U0[c] );
      _discrU0.set_copy( c, temp );
    }
  }
  ~ClassicalChanVeseVectorEnergyMultiPhiVariation () {}

  RealType evaluateIndicator ( const aol::Vector<RealType> &IndicatorParameter,
                               const int /*IndicatorParameterIndex*/,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint,
                               const typename ConfiguratorType::VecType &/*RefCoord*/ ) const {
    RealType indicator = 0.;
    for ( int i = 0; i < ImageDimension; i++ )
      indicator += aol::Sqr(_discrU0[i].evaluateAtQuadPoint( El, QuadPoint) - IndicatorParameter[i] );
    return indicator;
  }
};

/**
 * Calculates the volume of each segment given by the levelset functions.
 *
 * Argument in apply/applyAdd is ignored.
 *
 * Is derived from ChanVeseEnergyParameterDerivativeInterface even though
 * this is not a derivative. This can be done nevertheless.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions>
class MultiLevelsetVolumes
: public aol::ChanVeseEnergyParameterDerivativeInterface< ConfiguratorType,
                                                          HeavisideFunctionType,
                                                          NumberOfLevelsetFunctions,
                                                          1,
                                                          MultiLevelsetVolumes<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  MultiLevelsetVolumes( const typename ConfiguratorType::InitType &Initializer,
                        const aol::MultiVector<RealType> &LevelsetFunctions,
                        const HeavisideFunctionType &HeavisideFunction )
  : aol::ChanVeseEnergyParameterDerivativeInterface< ConfiguratorType,
                                                     HeavisideFunctionType,
                                                     NumberOfLevelsetFunctions,
                                                     1,
                                                     MultiLevelsetVolumes<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions> >
                                                   ( Initializer, LevelsetFunctions, HeavisideFunction )
  {
  }
  void evaluateIndicatorDerivative ( const aol::Vector<RealType> &/*IndicatorParameter*/,
                                     const int /*IndicatorParameterIndex*/,
                                     const typename ConfiguratorType::ElementType &/*El*/,
                                     int /*QuadPoint*/,
                                     const typename ConfiguratorType::VecType &/*RefCoord*/,
                                     aol::Vec<1, RealType> &Derivative ) const {
    Derivative[0] = 1.;
  }
};

/**
 * Calculates the volume of each segment given by the levelset functions,
 * weighted by the gray value of the image U0. If you divide the results
 * of this class component-wise by the ones of MultiLevelsetVolumes, you
 * get the vector containing the average gray values of the image in the
 * different segments.
 *
 * This also works for vector valued images U0: In this case you have to
 * divide the vector valued components of the result by the appropriate
 * scalar values (vectors of length one) obtained from MultiLevelsetVolumes.
 *
 * Argument in apply/applyAdd is ignored.
 *
 * Is derived from ChanVeseEnergyParameterDerivativeInterface even though
 * this is not a derivative. This can be done nevertheless.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int ImageDimension>
class MultiLevelsetWeightedVolumes
: public aol::ChanVeseEnergyParameterDerivativeInterface< ConfiguratorType,
                                                          HeavisideFunctionType,
                                                          NumberOfLevelsetFunctions,
                                                          ImageDimension,
                                                          MultiLevelsetWeightedVolumes<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, ImageDimension> > {
  aol::auto_container< ImageDimension, const aol::DiscreteFunctionDefault<ConfiguratorType> > _discrU0;
public:
  typedef typename ConfiguratorType::RealType RealType;
  MultiLevelsetWeightedVolumes( const typename ConfiguratorType::InitType &Initializer,
                                const aol::MultiVector<RealType> &LevelsetFunctions,
                                const HeavisideFunctionType &HeavisideFunction,
                                const aol::MultiVector<RealType> &U0 )
  : aol::ChanVeseEnergyParameterDerivativeInterface< ConfiguratorType,
                                                     HeavisideFunctionType,
                                                     NumberOfLevelsetFunctions,
                                                     ImageDimension,
                                                     MultiLevelsetWeightedVolumes<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, ImageDimension> >
                                                   ( Initializer, LevelsetFunctions, HeavisideFunction )
  {
    for ( int c=0; c<ImageDimension; c++ ) {
      aol::DiscreteFunctionDefault<ConfiguratorType> temp( Initializer, U0[c] );
      _discrU0.set_copy( c, temp );
    }
  }
  void evaluateIndicatorDerivative ( const aol::Vector<RealType> &/*IndicatorParameter*/,
                                     const int /*IndicatorParameterIndex*/,
                                     const typename ConfiguratorType::ElementType &El,
                                     int QuadPoint,
                                     const typename ConfiguratorType::VecType &/*RefCoord*/,
                                     aol::Vec<ImageDimension, RealType> &Derivative) const {
    for ( int i = 0; i < ImageDimension; i++ )
      Derivative[i] = _discrU0[i].evaluateAtQuadPoint( El, QuadPoint);
  }
};


/**
 * Calculates the average mean values (weighted according to LevelsetFunctions
 * and HeavisideFunction) in the segments for the classical piecewise constant
 * Chan Vese segmentation.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename HeavisideFunctionType, int NumberOfLevelsetFunctions, int ImageDimension>
class ClassicalChanVeseMeanValueUpdater {
private:
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const HeavisideFunctionType &_heavisideFunction;
  const aol::MultiVector<RealType> &_imageMVec;
public:
  ClassicalChanVeseMeanValueUpdater( const typename ConfiguratorType::InitType &Initializer,
                                     const HeavisideFunctionType &HeavisideFunction,
                                     const aol::MultiVector<RealType> &ImageMVec )
    : _grid ( Initializer ),
      _heavisideFunction ( HeavisideFunction ),
      _imageMVec( ImageMVec )
  {
  }
  void update( const aol::MultiVector<RealType> &LevelsetFunctions, aol::MultiVector<RealType> &MeanGrayValues ) const {
    // Calculate the new average values based on the current levelset functions.
    aol::MultiLevelsetVolumes<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions>
      volE( _grid, LevelsetFunctions, _heavisideFunction );

    aol::MultiLevelsetWeightedVolumes<ConfiguratorType, HeavisideFunctionType, NumberOfLevelsetFunctions, ImageDimension>
      volEweight( _grid, LevelsetFunctions, _heavisideFunction, _imageMVec );

    aol::MultiVector<RealType> tempMVec( MeanGrayValues.numComponents(), 1 );
    aol::MultiVector<RealType> tempGrayValues( MeanGrayValues, aol::STRUCT_COPY );
    volE.apply( tempMVec, tempMVec );
    // MultiLevelsetWeightedVolumes only uses the structure of the first argument, not the actual content.
    volEweight.apply( MeanGrayValues, tempGrayValues );
    for( int i = 0; i < MeanGrayValues.numComponents(); i++ ) {
      // We can only calculate a meaningful mean gray value for non-empty segments.
      if ( appeqAbsolute ( tempMVec[i][0], aol::ZOTrait<RealType>::zero ) == false ) {
        MeanGrayValues[i] = tempGrayValues[i];
        MeanGrayValues[i] /= tempMVec[i][0];
      }
    }
  }
};

} // namespace aol

#endif
