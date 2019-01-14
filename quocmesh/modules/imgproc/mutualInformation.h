#ifndef __MUTUALINFORMATION_H
#define __MUTUALINFORMATION_H

#include <gradientDescent.h>
#include <registration.h>
#include <jointHistogram.h>
#include <configurators.h>
#include <gnuplotter.h>
#include <deformations.h>
#include <paramReg.h>

namespace im {

template <typename ConfiguratorType, typename ConfiguratorTypeFeatureDomain> class MIRegistrationEnergyWithRegardToPhi;
template <typename ConfiguratorType> class VariationOfMIRegistrationEnergyWithRegardToPhi;

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class MIRegistrationConfigurator {
  typedef typename ConfiguratorType::RealType RealType;
public:
  typedef typename ConfiguratorType::ArrayType ImageDOFType;
private:
  const int _numberOfIntensityLevels;
  const RealType _beta;

public:
  MIRegistrationConfigurator ( const aol::ParameterParser &Parser )
  : _numberOfIntensityLevels ( Parser.getInt ( "numberOfIntensityLevels" ) ),
    _beta ( static_cast<RealType> ( Parser.getDouble ( "beta_hist" ) ) ) {
  }

  template <typename RegistrationMultilevelDescentType>
  MIRegistrationConfigurator ( const RegistrationMultilevelDescentType &RegisMLD )
  : _numberOfIntensityLevels ( RegisMLD.getParserReference().getInt ( "numberOfIntensityLevels" ) ),
    _beta ( static_cast<RealType> ( RegisMLD.getParserReference().getDouble ( "beta_hist" ) ) ) {
  }

  MIRegistrationConfigurator ( const int NumberOfIntensityLevels, const RealType Beta )
  : _numberOfIntensityLevels ( NumberOfIntensityLevels ),
    _beta ( Beta ) {
  }

  void checkInput ( const typename ConfiguratorType::ArrayType &ReferenceImage, const typename ConfiguratorType::ArrayType &TemplateImage ) const {
    if ( ( ReferenceImage.checkRange() == false ) || ( TemplateImage.checkRange() == false ) )
      throw aol::OutOfBoundsException ( "MIRegistrationConfigurator::checkInput:\nValues of the input images are not in [0,1]!\n", __FILE__, __LINE__ );
  }

  int getNumberOfIntensityLevels ( ) const {
    return _numberOfIntensityLevels;
  }

  RealType getBeta ( ) const {
    return _beta;
  }

  template <typename RegistrationMultilevelDescentType>
  void writeCurrentHistogram ( const RegistrationMultilevelDescentType &RegisMLD,
                               const aol::MultiVector<RealType> &Phi,
                               const char* FileName ) const {
    typename RegistrationMultilevelDescentType::ArrayType deformedTemplateArray ( RegisMLD.getCurrentGrid () );
    qc::DeformImage<ConfiguratorType> ( RegisMLD.getTemplImageMLAReference().current(), RegisMLD.getCurrentGrid (), deformedTemplateArray, Phi );
    im::JointHistogram<RealType> _hist ( RegisMLD.getRefImageMLAReference().current(), deformedTemplateArray, getNumberOfIntensityLevels(), getBeta() );
    qc::ScalarArray<RealType, qc::QC_2D> curTable ( _hist.getJointHistogram(), aol::FLAT_COPY );

    // output the table
    curTable.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0., curTable.getMaxValue() / 255 );
    curTable.savePNG ( FileName );
  }

  template <typename RegistrationMultilevelDescentType>
  void writeInitialRegistration ( const RegistrationMultilevelDescentType &RegisMLD,
                                  const qc::MultiArray<RealType, ConfiguratorType::Dim> &Phi ) const {
    RegisMLD.writeInitialRegistration( Phi );
    string filename = aol::strprintf ( "%sjointHistogram_before_%02d.png", RegisMLD.getSaveDirectory(), RegisMLD.getLevel() );
    writeCurrentHistogram ( RegisMLD, Phi, filename.c_str() );
  }

  template <typename RegistrationMultilevelDescentType>
  void writeCurrentRegistration ( const RegistrationMultilevelDescentType &RegisMLD,
                                  const qc::MultiArray<RealType, ConfiguratorType::Dim> &Phi ) const {
    string filename = aol::strprintf ( "%sjointHistogram%02d.png", RegisMLD.getSaveDirectory(), RegisMLD.getLevel() );
    writeCurrentHistogram ( RegisMLD, Phi, filename.c_str() );
    RegisMLD.writeCurrentRegistration ( Phi );
  }

  typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType, qc::QC_2D, 3> > ConfTypeFeatureDomain;
  typedef MIRegistrationEnergyWithRegardToPhi<ConfiguratorType, ConfTypeFeatureDomain> Energy;
  typedef VariationOfMIRegistrationEnergyWithRegardToPhi<ConfiguratorType> EnergyVariation;
};

template <typename ConfiguratorType>
class MIForce: public aol::FENonlinVectorOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, MIForce<ConfiguratorType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const int _numberOfIntensityValues;
  const qc::ScalarArray<RealType, qc::QC_2D> &_mutualInformation;
  const typename ConfiguratorType::InitType &_grid;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _t;
public:
  MIForce ( const typename ConfiguratorType::InitType &Grid,
            const aol::Vector<RealType> &Template,
            im::JointHistogram<RealType> &Hist )
    : aol::FENonlinVectorOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, MIForce<ConfiguratorType> > ( Grid ),
      _numberOfIntensityValues ( Hist.getNumberOfIntensityValues() ),
      _mutualInformation ( Hist.getMutualInformation() ),
      _grid ( Grid ),
      _r ( Grid, Hist.getReference() ),
      _t ( Grid, Template ) {}

  void getNonlinearity ( aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint,
                         const typename ConfiguratorType::VecType &RefCoord,
                         aol::Vec<ConfiguratorType::Dim, typename ConfiguratorType::RealType> &NL ) const {
    typename ConfiguratorType::VecType transformed_local_coord;
    qc::Element transformed_el;
    if ( !qc::transformCoord<ConfiguratorType> ( _grid, DiscFuncs, El, QuadPoint, RefCoord, transformed_el, transformed_local_coord ) ) {
      NL.setZero();
      return;
    }

    const RealType Td = _t.evaluate ( transformed_el, transformed_local_coord );
    const RealType R = _r.evaluateAtQuadPoint ( El, QuadPoint );
    const int ir = aol::Min ( static_cast<int> ( ( _numberOfIntensityValues - 1 ) * R ), _numberOfIntensityValues - 1 );
    const int it = aol::Min ( static_cast<int> ( ( _numberOfIntensityValues - 1 ) * Td ), _numberOfIntensityValues - 1 );

    const RealType factor = _mutualInformation.get ( ir, it );
    _t.evaluateGradient ( transformed_el, transformed_local_coord, NL );

    NL *= factor;
  }
};

template <typename ConfiguratorType,
typename ConfiguratorTypeFeatureDomain = qc::QuocConfiguratorTraitMultiLin<typename ConfiguratorType::RealType, qc::QC_2D, aol::GaussQuadrature<typename ConfiguratorType::RealType, qc::QC_2D, 3> > >
class MIRegistrationEnergyWithRegardToPhi: public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  const typename ConfiguratorType::InitType &_grid;

  const qc::Array<RealType> &_reference;
  const qc::Array<RealType> &_template;
  mutable qc::Array<RealType> _deformedTemplate;

  const MIRegistrationConfigurator<ConfiguratorType> _regisConfig;

public:

  MIRegistrationEnergyWithRegardToPhi ( const typename ConfiguratorType::InitType &Initializer,
                                        const qc::Array<RealType> &Reference,
                                        const qc::Array<RealType> &Template,
                                        const MIRegistrationConfigurator<ConfiguratorType> &RegisConfig ) :
      _grid ( Initializer ),
      _reference ( Reference ),
      _template ( Template ),
      _deformedTemplate ( _grid ),
      _regisConfig ( RegisConfig ) {}
  virtual void apply ( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    RealType energy = 0.;

    qc::DeformImage<ConfiguratorType> ( _template, _grid, _deformedTemplate, Arg );

    im::JointHistogram<RealType> hist ( _reference, _deformedTemplate, _regisConfig.getNumberOfIntensityLevels(), _regisConfig.getBeta() );

    const qc::ScalarArray<RealType, qc::QC_2D> &jointHistogram = hist.getJointHistogram();
/*
    //-----------------------------------------------------------------------

    const typename ConfiguratorTypeFeatureDomain::InitType &histoGrid ( hist.getHistoGrid() );

    const aol::DiscreteFunctionDefault<ConfiguratorTypeFeatureDomain> discrJointHistogram( histoGrid, jointHistogram );

    qc::ScalarArray<RealType, qc::QC_2D> p1Array(histoGrid);
    qc::ScalarArray<RealType, qc::QC_2D> p2Array(histoGrid);

    hist.marginalize(p1Array, 0);
    hist.marginalize(p2Array, 1);

    const aol::DiscreteFunctionDefault<ConfiguratorTypeFeatureDomain> discrP1( histoGrid, p1Array );
    const aol::DiscreteFunctionDefault<ConfiguratorTypeFeatureDomain> discrP2( histoGrid, p2Array );

    ConfiguratorTypeFeatureDomain config(histoGrid);
    typedef typename ConfiguratorTypeFeatureDomain::QuadType QType;

    typedef typename ConfiguratorTypeFeatureDomain::ElementIteratorType IteratorType;
    for ( IteratorType it = config.begin(); it != config.end(); ++it ) {
      RealType a = 0.;
      for ( int q = 0; q < QType::numQuadPoints; q++ ) {
        const RealType p1p2 = discrP1.evaluateAtQuadPoint( *it, q ) * discrP2.evaluateAtQuadPoint( *it, q );
        if ( p1p2 != 0 ) {
          RealType value = discrJointHistogram.evaluateAtQuadPoint( *it, q );
          if ( aol::appeqAbsolute ( value, aol::ZOTrait<RealType>::zero ) == false ) {
            value = -1 * value * log ( value / p1p2 );
            a += value * config.getBaseFunctionSet(*it).getWeight( q );
          }
        }
      }
      a *= aol::Sqr ( config.H( *it ) );
      energy += a;
    }

    energy = 0.;
    //-----------------------------------------------------------------------
*/
    const int numberOfIntensityValues = hist.getNumberOfIntensityValues();
    aol::Vector<RealType> p1 ( numberOfIntensityValues );
    aol::Vector<RealType> p2 ( numberOfIntensityValues );

    hist.marginalizeVector ( p1, 0 );
    hist.marginalizeVector ( p2, 1 );

    for ( int j = 0; j < numberOfIntensityValues; ++j ) {
      if ( p2[j] != 0 ) {
        for ( int i = 0; i < numberOfIntensityValues; ++i ) {
          if ( p1[i] != 0 ) {
            const RealType value = jointHistogram.get ( i, j );
            // The limit x -> 0 of x log ( a x ) for a > 0 is 0. To work around numerical errors
            // ( log ( epsilon ) numerically is -inf for epsilon > 0 small enough )
            // we explicitly catch this case here and don't need to add anything then.
            if ( aol::appeqAbsolute ( value, aol::ZOTrait<RealType>::zero ) == false ) {
              energy += -1 * value * log ( value / ( p1[i] * p2[j] ) );
            }
          }
        }
      }
    }

    // BB: The scaling is necessary, otherwise timestep control doesn't work, but I can't find the reason for the scaling
    energy *= aol::Sqr ( numberOfIntensityValues - 1 );

    Dest[0] = energy;
  }
  virtual void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    aol::Scalar<RealType> tmp;
    apply( MArg, tmp );
    Dest += tmp;
  }

  const typename ConfiguratorType::InitType &getGridReference ( ) const {
    return _grid;
  }
};


template <typename ConfiguratorType>
class VariationOfMIRegistrationEnergyWithRegardToPhi: public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  const typename ConfiguratorType::InitType &_grid;

  const qc::Array<RealType> &_reference;
  const qc::Array<RealType> &_template;
  mutable qc::Array<RealType> _deformedTemplate;

  const MIRegistrationConfigurator<ConfiguratorType> _regisConfig;

public:

  VariationOfMIRegistrationEnergyWithRegardToPhi ( const typename ConfiguratorType::InitType &Initializer,
                                                   const qc::Array<RealType> &Reference,
                                                   const qc::Array<RealType> &Template,
                                                   const MIRegistrationConfigurator<ConfiguratorType> &RegisConfig ) :
      _grid ( Initializer ),
      _reference ( Reference ),
      _template ( Template ),
      _deformedTemplate ( _grid ),
      _regisConfig ( RegisConfig ) {}
  virtual void apply ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {

    qc::DeformImage<ConfiguratorType> ( _template, _grid, _deformedTemplate, Arg );
    im::JointHistogram<RealType> _hist ( _reference, _deformedTemplate, _regisConfig.getNumberOfIntensityLevels(), _regisConfig.getBeta() );
    MIForce<ConfiguratorType> miForceOp ( _grid, _template, _hist );

    miForceOp.apply ( Arg, Dest );
  }
  virtual void applyAdd ( const aol::MultiVector<RealType> &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const {
    throw aol::Exception ( "Not implemented", __FILE__, __LINE__ );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricMIForce
: public aol::VectorFENonlinIntegrationVectorInterface < ConfiguratorType,
                                                         ParametricMIForce<ConfiguratorType, ParametricDeformationType>,
                                                         1, ParametricDeformationType::NumberOfDeformParameters > {
  typedef typename ConfiguratorType::RealType RealType;
  const int _numberOfIntensityValues;
  const qc::ScalarArray<RealType, qc::QC_2D> &_mutualInformation;
  const typename ConfiguratorType::InitType &_grid;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const ParametricDeformationType &_parDef;
public:
  ParametricMIForce ( const typename ConfiguratorType::InitType &Grid,
                      im::JointHistogram<RealType> &Hist,
                      const ParametricDeformationType &ParDef )
    : aol::VectorFENonlinIntegrationVectorInterface < ConfiguratorType,
                                                      ParametricMIForce<ConfiguratorType, ParametricDeformationType>,
                                                      1, ParametricDeformationType::NumberOfDeformParameters >  ( Grid ),
      _numberOfIntensityValues ( Hist.getNumberOfIntensityValues() ),
      _mutualInformation ( Hist.getMutualInformation() ),
      _grid ( Grid ),
      _r ( Grid, Hist.getReference() ),
      _parDef ( ParDef ) {}

  void evaluateIntegrand ( const aol::auto_container<1, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrFuncs, // DiscrFuncs[0] is supposed to contain T.
                          const typename ConfiguratorType::ElementType &El,
                          int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,
                          aol::Vec<ParametricDeformationType::NumberOfDeformParameters, RealType> &Integrand ) const {

    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;
    if ( !_parDef.template evaluateDeformation<false> ( El, RefCoord, transformedEl, transformedLocalCoord) ) {
      Integrand.setZero();
      return;
    }

    aol::Mat< ParametricDeformationType::NumberOfDeformParameters, ConfiguratorType::Dim, RealType > jacobian;
    _parDef.evaluateDerivativeDeformation ( El, RefCoord, jacobian );

    typename ConfiguratorType::VecType gradT;
    DiscrFuncs[0].evaluateGradient ( transformedEl, transformedLocalCoord, gradT );
    jacobian.mult ( gradT, Integrand );

    const RealType Td = DiscrFuncs[0].evaluate ( transformedEl, transformedLocalCoord );
    const RealType R = _r.evaluateAtQuadPoint ( El, QuadPoint );
    const int ir = aol::Min ( static_cast<int> ( ( _numberOfIntensityValues - 1 ) * R ), _numberOfIntensityValues - 1 );
    const int it = aol::Min ( static_cast<int> ( ( _numberOfIntensityValues - 1 ) * Td ), _numberOfIntensityValues - 1 );
    Integrand *= _mutualInformation.get ( ir, it );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricMIEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
  const typename ConfiguratorType::ArrayType _r;
  const typename ConfiguratorType::ArrayType _t;
  const im::MIRegistrationConfigurator<ConfiguratorType> _regisConfig;
  im::MIRegistrationEnergyWithRegardToPhi<ConfiguratorType> _MIEnergy;
public:
  typedef typename ConfiguratorType::RealType RealType;

  ParametricMIEnergy ( const typename ConfiguratorType::InitType &Grid,
                       const aol::Vector<RealType> &ImR,
                       const aol::Vector<RealType> &ImT )
    : _r ( ImR, Grid, aol::FLAT_COPY ),
      _t ( ImT, Grid, aol::FLAT_COPY ),
      _regisConfig ( 6, 1 ),
      _MIEnergy ( Grid, _r, _t, _regisConfig ) {}

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    // This implementation is not very efficient but is hopefully sufficient to test how well
    // MI works in the parametric case.
    qc::DataGenerator<ConfiguratorType> generator ( _MIEnergy.getGridReference() );
    qc::MultiArray<RealType, ConfiguratorType::Dim> deformation ( _MIEnergy.getGridReference() );
    ParametricDeformationType parDef ( _MIEnergy.getGridReference(), MArg );
    generator.template generateDeformationFromParametricDeformation<ParametricDeformationType, true> ( parDef, deformation );
    _MIEnergy.applyAdd ( deformation, Dest );
  }

  void applyDerivative ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    // This implementation is not very efficient but is hopefully sufficient to test how well
    // MI works in the parametric case.
    typename ConfiguratorType::ArrayType deformedTemplate ( _MIEnergy.getGridReference() );
    qc::ParametricDeformImage<ConfiguratorType, ParametricDeformationType> ( _t, _MIEnergy.getGridReference(), deformedTemplate, MArg );
    im::JointHistogram<RealType> hist ( _r, deformedTemplate, _regisConfig.getNumberOfIntensityLevels(), _regisConfig.getBeta() );

    aol::Vector<RealType> destVec ( ParametricDeformationType::NumberOfDeformParameters );
    ParametricDeformationType parDef (  _MIEnergy.getGridReference(), MArg );
    ParametricMIForce<ConfiguratorType, ParametricDeformationType> force ( _MIEnergy.getGridReference(), hist, parDef );
    force.apply ( _t, destVec );
    MDest.copySplitFrom ( destVec );
  }
};

} // end namespace qc

#endif // __MUTUALINFORMATION_H
