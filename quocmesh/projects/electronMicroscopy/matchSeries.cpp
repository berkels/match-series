#include <matchSeries.h>
#include <cellCenteredGrid.h>
#include <paramReg.h>
#include <mutualInformation.h>

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
//typedef qc::QuocConfiguratorTraitMultiLin<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;
typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3>, qc::CellCenteredCubicGrid<DimensionChoice> > ConfType;

int main( int argc, char **argv ) {
  try {
    aol::ParameterParser parser( argc > 1 ? argv[1] : "matchSeries.par" );

    const int deformationModel = parser.getIntOrDefault ( "deformationModel", 0 );
    switch ( deformationModel ) {
      case 0:
        {
          typedef qc::StandardRegistrationMultilevelDescent<ConfType, qc::NCCRegistrationConfigurator<ConfType> > RegisType;
          matchSeries<ConfType, RegisType>( parser, argc, argv );
        }
        break;
      case 1:
        {
          typedef qc::ParametricRigidBodyMotion2D<ConfType> ParametricDefType;
          typedef qc::ParametricRegistrationMultilevelDescent<ConfType, ParametricDefType, qc::ParametricNCCEnergy<ConfType, ParametricDefType> > RegisType;
          matchSeries<ConfType, RegisType>( parser, argc, argv );
        }
        break;
      case 2:
        {
          typedef qc::ParametricTranslation<ConfType> ParametricDefType;
          typedef qc::ParametricRegistrationMultilevelDescent<ConfType, ParametricDefType, qc::ParametricNCCEnergy<ConfType, ParametricDefType> > RegisType;
          matchSeries<ConfType, RegisType>( parser, argc, argv );
        }
        break;
      case 3:
        {
          typedef qc::StandardRegistrationMultilevelDescent<ConfType, im::MIRegistrationConfigurator<ConfType> > RegisType;
          matchSeries<ConfType, RegisType>( parser, argc, argv );
        }
        break;
      case 4:
        {
          typedef qc::StandardRegistrationMultilevelDescent<ConfType, im::MIRegistrationConfigurator<ConfType>, qc::DirichletLaplaceConfigurator<ConfType> > RegisType;
          matchSeries<ConfType, RegisType>( parser, argc, argv );
        }
        break;
      case 5:
      {
        typedef qc::ParametricRigidBodyMotion2D<ConfType> ParametricDefType;
        typedef qc::ParametricRegistrationMultilevelDescent<ConfType, ParametricDefType, im::ParametricMIEnergy<ConfType, ParametricDefType> > RegisType;
        matchSeries<ConfType, RegisType>( parser, argc, argv );
      }
        break;
      case 6:
        {
          typedef qc::ParametricTranslation<ConfType> ParametricDefType;
          typedef qc::ParametricRegistrationMultilevelDescent<ConfType, ParametricDefType, im::ParametricMIEnergy<ConfType, ParametricDefType> > RegisType;
          matchSeries<ConfType, RegisType>( parser, argc, argv );
        }
        break;
      case 7:
        {
          // Prolongation of deformations so far only works with 2^d+1 grids without possibly generating
          // infinite energy after prolongation.
          typedef qc::QuocConfiguratorTraitMultiLin<RType, DimensionChoice, aol::Quadrature2D<RType,aol::SimpsonQuadrature<RType> > > HyperConfType;
          //typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::Quadrature2D<RType,aol::SimpsonQuadrature<RType> > , qc::CellCenteredCubicGrid<DimensionChoice> > HyperConfType;
          typedef qc::StandardRegistrationMultilevelDescent<HyperConfType, qc::NCCRegistrationConfigurator<HyperConfType>, qc::HyperelasticEnergyConfigurator<HyperConfType> > RegisType;
          matchSeries<HyperConfType, RegisType>( parser, argc, argv );
        }
        break;
      case 8:
        {
          typedef qc::FullyAffineTransformation2D<ConfType> ParametricDefType;
          typedef qc::ParametricRegistrationMultilevelDescent<ConfType, ParametricDefType, qc::ParametricNCCEnergy<ConfType, ParametricDefType> > RegisType;
          matchSeries<ConfType, RegisType>( parser, argc, argv );
        }
        break;
      case 9:
        {
          typedef qc::FullyAffineTransformation2D<ConfType> ParametricDefType;
          typedef qc::ParametricRegistrationMultilevelDescent<ConfType, ParametricDefType, qc::ParametricSSDEnergy<ConfType, ParametricDefType> > RegisType;
          matchSeries<ConfType, RegisType>( parser, argc, argv );
        }
        break;
      case 10:
        {
          typedef qc::StandardRegistrationMultilevelDescent<ConfType, qc::NCCRegistrationConfigurator<ConfType>, qc::LaplaceConfigurator<ConfType> > RegisType;
          matchSeries<ConfType, RegisType>( parser, argc, argv );
        }
        break;

      case 11:
        {
          typedef qc::RectangularGridConfigurator<RType, qc::QC_1D, aol::GaussQuadrature<RType, qc::QC_1D, 3> > ConfType1D;
          typedef XShiftRegistrationMultilevelDescent<ConfType, ConfType1D, qc::NCCRegistrationConfigurator<ConfType1D> > RegisType;
          matchSeries<ConfType, RegisType>( parser, argc, argv );
        }
        break;

      case 12:
        {
          typedef qc::ParametricTranslation<ConfType> ParametricDefType;
          typedef qc::ParametricRegistrationMultilevelDescent<ConfType, ParametricDefType, qc::ParametricSSDEnergy<ConfType, ParametricDefType> > RegisType;
          matchSeries<ConfType, RegisType>( parser, argc, argv );
        }
        break;
      case 13:
        {
          typedef qc::StandardRegistrationMultilevelDescent<ConfType, qc::CCRegistrationConfigurator<ConfType> > RegisType;
          matchSeries<ConfType, RegisType>( parser, argc, argv );
        }
        break;
      case 14:
        {
          typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > NonDyadicConfType;
          typedef NonDyadicRegistrationMultilevelDescent<NonDyadicConfType, qc::NCCRegistrationConfigurator, qc::DirichletRegularizationConfigurator, ConfType,
                                                         qc::StandardRegistrationMultilevelDescent<ConfType, qc::NCCRegistrationConfigurator<ConfType>, qc::DirichletRegularizationConfigurator<ConfType> >,
                                                         aol::H1GradientDescent<NonDyadicConfType, aol::MultiVector<RType>, qc::CholeskyBasedInverseH1Metric<NonDyadicConfType> >
                                                         > RegisType;
          matchSeries<ConfType, RegisType>( parser, argc, argv );
        }
        break;
      case 15:
        {
          typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::Quadrature2D<RType,aol::SimpsonQuadrature<RType> > > NonDyadicHyperConfType;
          typedef qc::QuocConfiguratorTraitMultiLin<RType, DimensionChoice, aol::Quadrature2D<RType,aol::SimpsonQuadrature<RType> > > HyperConfType;
          typedef NonDyadicRegistrationMultilevelDescent<NonDyadicHyperConfType, qc::NCCRegistrationConfigurator, qc::HyperelasticEnergyConfigurator, HyperConfType,
                                                         qc::StandardRegistrationMultilevelDescent<HyperConfType, qc::NCCRegistrationConfigurator<HyperConfType>, qc::HyperelasticEnergyConfigurator<HyperConfType> >,
                                                         aol::H1GradientDescent<NonDyadicHyperConfType, aol::MultiVector<RType>, qc::CholeskyBasedInverseH1Metric<NonDyadicHyperConfType> >
                                                        > RegisType;
          matchSeries<HyperConfType, RegisType>( parser, argc, argv );
        }
        break;
      case 16:
        {
          typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::Quadrature2D<RType,aol::SimpsonQuadrature<RType> > > NonDyadicHyperConfType;
          typedef qc::QuocConfiguratorTraitMultiLin<RType, DimensionChoice, aol::Quadrature2D<RType,aol::SimpsonQuadrature<RType> > > HyperConfType;
          typedef NonDyadicRegistrationMultilevelDescent<NonDyadicHyperConfType, qc::NCCRegistrationConfigurator, qc::HyperelasticEnergyConfigurator, HyperConfType,
                                                         qc::StandardRegistrationMultilevelDescent<HyperConfType, qc::NCCRegistrationConfigurator<HyperConfType>, qc::HyperelasticEnergyConfigurator<HyperConfType> >,
                                                         aol::H1GradientDescent<NonDyadicHyperConfType, aol::MultiVector<RType> >
                                                         > RegisType;
          matchSeries<HyperConfType, RegisType>( parser, argc, argv );
        }
        break;
      case 17:
        {
          typedef qc::StandardRegistrationMultilevelDescent<ConfType, qc::NCCRegistrationConfigurator<ConfType>, qc::DirichletLaplaceConfigurator<ConfType> > RegisType;
          matchSeries<ConfType, RegisType>( parser, argc, argv );
        }
        break;
      case 18:
        {
          typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::Quadrature2D<RType,aol::SimpsonQuadrature<RType> > > NonDyadicHyperConfType;
          typedef qc::QuocConfiguratorTraitMultiLin<RType, DimensionChoice, aol::Quadrature2D<RType,aol::SimpsonQuadrature<RType> > > HyperConfType;
          typedef NonDyadicRegistrationMultilevelDescent<NonDyadicHyperConfType, qc::SSDRegistrationConfigurator, qc::HyperelasticEnergyConfigurator, HyperConfType,
                                                         qc::StandardRegistrationMultilevelDescent<HyperConfType, qc::SSDRegistrationConfigurator<HyperConfType>, qc::HyperelasticEnergyConfigurator<HyperConfType> >,
                                                         aol::H1GradientDescent<NonDyadicHyperConfType, aol::MultiVector<RType>, qc::CholeskyBasedInverseH1Metric<NonDyadicHyperConfType> >
                                                        > RegisType;
          matchSeries<HyperConfType, RegisType>( parser, argc, argv );
        }
      case 19:
        {
          typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::Quadrature2D<RType,aol::SimpsonQuadrature<RType> > > NonDyadicHyperConfType;
          typedef qc::QuocConfiguratorTraitMultiLin<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > HyperConfType;
          typedef NonDyadicRegistrationMultilevelDescent<NonDyadicHyperConfType, qc::SSDRegistrationConfigurator, qc::HyperelasticEnergyConfigurator, HyperConfType,
                                                         qc::StandardRegistrationMultilevelDescent<HyperConfType, qc::SSDRegistrationConfigurator<HyperConfType>, qc::HyperelasticEnergyConfigurator<HyperConfType> >,
                                                         aol::H1GradientDescent<NonDyadicHyperConfType, aol::MultiVector<RType>, qc::CholeskyBasedInverseH1Metric<NonDyadicHyperConfType> >
                                                        > RegisType;
          matchSeries<HyperConfType, RegisType>( parser, argc, argv );
        }
        break;
      default:
        throw aol::Exception ( "Unknown mode", __FILE__, __LINE__ );
    }
  }
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
