#include <matchSeismicSeries.h>
#include <cellCenteredGrid.h>
#include <paramReg.h>
#include <mutualInformation.h>

typedef double RType;
const qc::Dimension DimensionChoice = qc::QC_2D;
typedef qc::RectangularGridConfigurator<RType, DimensionChoice, aol::GaussQuadrature<RType,DimensionChoice,3> > ConfType;

int main( int argc, char **argv ) {
  try {
    aol::ParameterParser parser( argc > 1 ? argv[1] : "matchSeries.par" );
    typedef qc::RectangularGridConfigurator<RType, qc::QC_1D, aol::GaussQuadrature<RType, qc::QC_1D, 3> > DyadicConfiguratorType1D;
    typedef qc::RectangularGridConfigurator<RType, qc::QC_2D, aol::GaussQuadrature<RType, qc::QC_2D, 3>, qc::CellCenteredCubicGrid<qc::QC_2D> > DyadicConfiguratorType;
    typedef aol::H1GradientDescent<DyadicConfiguratorType, aol::Vector<RType>, typename qc::MultilevelArrayTrait<RType, typename DyadicConfiguratorType::InitType>::LinSmoothType> DyadicGradientDescentType;
    //typedef aol::GridlessGradientDescent<RType, aol::Vector<RType> > DyadicGradientDescentType;
    typedef XShiftRegistrationMultilevelDescent<DyadicConfiguratorType, DyadicConfiguratorType1D, qc::NCCRegistrationConfigurator<DyadicConfiguratorType1D>, DyadicGradientDescentType> DyadicRegistrationMultilevelDescentType;
    typedef aol::H1GradientDescent<ConfType, aol::Vector<RType>, qc::CholeskyBasedInverseH1Metric<ConfType> > NonDyadicGradientDescentType;
    //typedef aol::GridlessGradientDescent<RType, aol::Vector<RType> > NonDyadicGradientDescentType;
    typedef XShiftNonDyadicRegistrationMultilevelDescent<ConfType, qc::NCCRegistrationConfigurator, qc::DirichletRegularizationConfigurator, DyadicConfiguratorType1D, DyadicConfiguratorType, DyadicRegistrationMultilevelDescentType, NonDyadicGradientDescentType> RegisType;

    if ( parser.checkAndGetBool ( "addParamsToSafeDirName" ) ) {
      string saveDirectory = parser.getStringExpandTilde( "saveDirectory" );
      if ( saveDirectory.back() == '/' )
        saveDirectory = saveDirectory.substr ( 0, saveDirectory.size() - 1 );

      std::ostringstream strs;
      strs << saveDirectory;
      if ( parser.checkAndGetBool ( "saveNamedDeformedTemplatesUsingNearestNeighborInterpolation" ) )
        strs << "-nointerp";
      strs << "-" <<  parser.getDouble ( "lambda" )
           << "-" << parser.getDouble ( "yRegFactor" )
           << "-" << parser.getDouble ( "laplaceWeight" )
           << "-" << parser.getDouble ( "averageXDerivNormThreshold" )
           << "-" << parser.getDouble ( "extraStagesLambdaFactor" );
      if ( parser.checkAndGetBool ( "forceNonDyadicXShiftYConst" ) )
        strs << "-yconst";
      strs << "/";
      saveDirectory = strs.str();
      parser.changeVariableValue ( "saveDirectory", saveDirectory.c_str() );
    }

    if ( parser.hasVariable( "initialAverage" ) == false ) {
      std::vector<std::string> templateFileNames;
      aol::createTemplateFileNameList ( parser, templateFileNames );

      const bool clamp = parser.hasVariable( "clampRange" );
      aol::Vec<2, RType> clampRange;
      if ( clamp )
        clampRange = parser.getRealVec<2, RType> ( "clampRange" );

      ConfType::ArrayType referenceImage ( templateFileNames[0] );
      if ( clamp )
        referenceImage.clamp ( clampRange[0], clampRange[1] );
      const int numTemplates = parser.getInt ( "refLineDistance" ) - 1;

      if ( numTemplates >= static_cast<int> ( templateFileNames.size() ) )
        throw aol::Exception ( "refLineDistance too big for the amount of supplied template images", __FILE__, __LINE__ );

      aol::RandomAccessContainer<ConfType::ArrayType> templateImages ( numTemplates );
      for ( int i = 0; i < numTemplates; ++i ) {
        templateImages[i].load ( templateFileNames[i+1].c_str() );
        if ( clamp )
          templateImages[i].clamp ( clampRange[0], clampRange[1] );
      }
      const ConfType::InitType grid ( qc::GridSize<ConfType::Dim>::createFrom ( referenceImage ) );
      ConfType::ArrayType interleavedImage ( grid );

      interleaveImages<ConfType, RegisType::DyadicConfiguratorType1D, RegisType::DyadicRegistrationMultilevelDescentType::RegistrationConfiguratorType1D> ( referenceImage, templateImages, 1, interleavedImage, parser.getBool ( "justStackIniAverage" ) == false, parser.getBoolOrDefault ( "justCopyLastRefLineForIniAverage", false ) );

      qc::DefaultArraySaver<RType, ConfType::Dim> saver;
      saver.initFromParser ( parser, true );

      const string initialAverageFileName = saver.createSaveName ( "", qc::getDefaultArraySuffix ( ConfType::Dim ), -1, "initialAverage" );
      interleavedImage.save ( initialAverageFileName.c_str(), qc::PGM_DOUBLE_BINARY );
      parser.addVariable( "initialAverage", initialAverageFileName.c_str() );
    }

    matchSeries<ConfType, RegisType>( parser, argc, argv );
  }
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
