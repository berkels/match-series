#ifndef __MATCHSERIES_H
#define __MATCHSERIES_H

#include <registration.h>
#include <dm3Import.h>
#include <paramReg.h>
#include <multiStreambuf.h>
#include <IQFactor.h>
#include <micrographQuantifier.h>
#ifdef USE_INTERNAL_MODULES_IMGPROC
#include <noiseAnalysis.h>
#endif

/**
 * \author Berkels
 */
template <typename RegistrationType>
class SeriesMatching {
  typedef typename RegistrationType::ConfiguratorType ConfiguratorType;
  typedef typename RegistrationType::TransformationDOFType TransformationDOFType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef typename ConfiguratorType::RealType RealType;
  static const qc::Dimension Dim = ConfiguratorType::Dim;

  const aol::ParameterParser &_parser;
  RegistrationType _registrationAlgo;

  const bool _loadStageOneResults;
  const bool _reduceDeformations;
  const bool _useAltStartLevel;
  const int _stage;
public:
  enum ACTION {
    MATCH_AND_AVERAGE_SERIES,
    ONLY_AVERAGE_SERIES,
    ANALYZE_DEFORMATTIONS,
    APPLY_DEFORMATTION
  };

  SeriesMatching ( const aol::ParameterParser &Parser, const bool LoadStageOneResults = false )
    : _parser ( Parser ),
      _registrationAlgo ( Parser ),
      _loadStageOneResults ( LoadStageOneResults ),
      _reduceDeformations ( _parser.checkAndGetBool ( "reduceDeformations" ) ),
      _useAltStartLevel ( _parser.hasVariable ( "altStartLevel" ) && ( _parser.getInt ( "altStartLevel" ) != _parser.getInt ( "startLevel" ) ) ),
      _stage ( _parser.hasVariable ( "stage" ) ? _parser.getInt ( "stage" ) : 1 ) { }

  int getNumTemplates ( ) const {
    // We don't need these filenames here, but this is a convenient way to calculate
    // the number of templates we have.
    std::vector<std::string> templateFileNames;
    aol::createTemplateFileNameList ( _parser, templateFileNames );
    return static_cast<int> ( templateFileNames.size() );
  }

  string createDeformationBaseFileName ( const char *InputDirectory, const int DeformationNumber ) const {
    if ( _useAltStartLevel && ( DeformationNumber == 0 ) )
      return aol::strprintf ( "%sdeformation-%03d", InputDirectory, DeformationNumber );
    else {
      const bool refinedDeformation = ( ( _parser.checkAndGetBool ( "dontAccumulateDeformation" ) == false ) && ( DeformationNumber > 0 ) );
      return aol::strprintf ( "%s%d%s/deformation_%02d", InputDirectory, DeformationNumber, refinedDeformation ? "-r" : "", refinedDeformation ? _registrationAlgo.getMaxGridDepth() : _parser.getInt ( "stopLevel" ) );
    }
  }

  string createDeformationFileName ( const char *InputDirectory, const int DeformationNumber ) const {
    return ( createDeformationBaseFileName ( InputDirectory, DeformationNumber ) + _registrationAlgo.getDeformationFileNameSuffix() );
  }

  void writeCurrentDefParams ( const int I, std::ofstream &DefParamsFile ) const {
    if ( RegistrationType::IsParametric ) {
      TransformationDOFType params ( _registrationAlgo.getTransformationDOFInitializer() );
      _registrationAlgo.getTransformation ( params );
      DefParamsFile << I << " ";
      if ( _parser.checkAndGetBool ( "saveTranslationInPixels" ) && ( params.numComponents() == 1 ) && ( params[0].size() == 2 ) ) {
        for ( int j = 0; j < params[0].size(); ++j )
          DefParamsFile << params[0][j] / _registrationAlgo.H() << " ";
      }
      else {
        for ( int i = 0; i < params.numComponents(); ++i )
          for ( int j = 0; j < params[i].size(); ++j )
            DefParamsFile << params[i][j] << " ";
      }
      DefParamsFile << endl;
    }
  }

  void matchSeries ( ) {
    // We will temporarily change the save directory, so store the original value.
    const string origSaveDir = _registrationAlgo.getSaveDirectory();

    std::ofstream defNormsFile;
    if ( _loadStageOneResults == false )
      defNormsFile.open ( ( origSaveDir + "/defNorms.txt" ).c_str() );

    std::ofstream energiesFile;
    energiesFile.open ( ( origSaveDir + "/energies.txt" ).c_str() );

    std::ofstream defParamsFile;
    if ( RegistrationType::IsParametric )
      defParamsFile.open ( ( origSaveDir + "/defParams.txt" ).c_str() );

    std::vector<std::string> templateFileNames;
    aol::createTemplateFileNameList ( _parser, templateFileNames );
    const int numTemplateImages = static_cast<int> ( templateFileNames.size() );
    const ArrayType reference ( _registrationAlgo.getRefImageReference(), aol::DEEP_COPY );
    const bool reverseRoles = _parser.checkAndGetBool ( "reverseRolesInSeriesMatching" );
    const bool calcInverseDeformation = _parser.checkAndGetBool ( "calcInverseDeformation" );
    if ( reverseRoles )
      _registrationAlgo.loadRefOrTemplate ( _parser.getStringExpandTilde ( "reference" ).c_str(), qc::TEMPLATE, _parser.checkAndGetBool ( "dontNormalizeInputImages" ) );

    // Start with a clean deformation.
    _registrationAlgo.setTransformationToZero ( );
    TransformationDOFType accumulatedTransformation ( _registrationAlgo.getTransformationDOFInitializer() );
    // The "zero transformation" (aka the identity mapping) doesn't necessarily have zero DOFs.
    // So make sure it's initialized properly by using the clean deformation in _registrationAlgo.
    _registrationAlgo.getTransformation ( accumulatedTransformation );

    for ( int i = 0; i < numTemplateImages; ++i ) {
      // First match the current template with the last template (assuming that the first template is the reference image).

      // In case we try to reduce the deformations, we have to adjust the initial deformation used
      // for the first pairing.
      if ( _reduceDeformations && ( reverseRoles == false ) && ( i == 0 ) && ( _stage > 1 ) ) {
        _registrationAlgo.loadTransformation ( aol::strprintf ( "%s../stage%d/reduceDef%s", _parser.getStringExpandTilde ( "saveDirectory" ).c_str(), _stage - 1, _registrationAlgo.getDeformationFileNameSuffix().c_str() ).c_str() );
        // In Stage 3+, this needs to be combined with the deformation calculated in the previous stage.
        if ( _stage > 2 ) {
          TransformationDOFType reducedTransformation ( _registrationAlgo.getTransformationDOFInitializer() );
          _registrationAlgo.getTransformation ( reducedTransformation );

          _registrationAlgo.loadTransformation ( aol::strprintf ( "%s../stage%d/0/deformation_%02d%s", _parser.getStringExpandTilde ( "saveDirectory" ).c_str(), _stage-1, _registrationAlgo.getMaxGridDepth(), _registrationAlgo.getDeformationFileNameSuffix().c_str() ).c_str() );
          TransformationDOFType transformationToLastAverage ( _registrationAlgo.getTransformationDOFInitializer() );
          _registrationAlgo.getTransformation ( transformationToLastAverage );

          _registrationAlgo.setTransformationToComposition ( transformationToLastAverage, reducedTransformation );
        }

      }
      else // Start with a clean displacement (this pairing is new).
        _registrationAlgo.setTransformationToZero();

      _registrationAlgo.loadRefOrTemplate ( templateFileNames[i].c_str(), reverseRoles ? qc::REFERENCE : qc::TEMPLATE, _parser.checkAndGetBool ( "dontNormalizeInputImages" ) );
      if ( ( _loadStageOneResults == false ) || ( i == 0 ) ) {
        _registrationAlgo.makeAndSetSaveDirectory ( aol::strprintf ( "%s%d/", _parser.getStringExpandTilde ( "saveDirectory" ).c_str(), i ).c_str() );

        if ( _parser.checkAndGetBool ( "useCorrelationToInitTranslation" ) ) {
          qc::CoordType signedShift = qc::PhaseCorrelationRegistration<ConfiguratorType>::registerImages ( _registrationAlgo.getTemplImageReference(), _registrationAlgo.getRefImageReference(), _parser.getInt ( "maxCorrShift" ) );
        
          aol::Vec<Dim, RealType> trans;

          for ( int k = 0; k < ConfiguratorType::Dim; ++k )
            trans[k] = signedShift[k] * _registrationAlgo.getInitializerRef().H();
          _registrationAlgo.setTransformationToTranslation ( trans );
        }

        _registrationAlgo.solveAndProlongToMaxDepth();
        if ( _loadStageOneResults == false )
          defNormsFile << i << " " << _registrationAlgo.getTransformationNorm () << " # " << templateFileNames[i] << endl;

        if ( _useAltStartLevel ) {
          const RealType firstTryEnergy = _registrationAlgo.getEnergyOfLastSolution();

          TransformationDOFType savedTransformation ( _registrationAlgo.getTransformationDOFInitializer() );
          _registrationAlgo.getTransformation ( savedTransformation );

          ArrayType temp ( _registrationAlgo.getInitializerRef() );
          _registrationAlgo.applyCurrentTransformation ( _registrationAlgo.getTemplImageReference(), temp );
          const int numOutOfDomainFirstTry = temp.numOccurence ( aol::NumberTrait<RealType>::Inf );

          _registrationAlgo.setTransformationToZero();
          _registrationAlgo.makeAndSetSaveDirectory ( aol::strprintf ( "%s%d-alt/", _parser.getStringExpandTilde ( "saveDirectory" ).c_str(), i ).c_str() );
          _registrationAlgo.solveAndProlongToMaxDepth(  _parser.getInt ( "altStartLevel" ) );

          _registrationAlgo.applyCurrentTransformation ( _registrationAlgo.getTemplImageReference(), temp );
          const int numOutOfDomainSecondTry = temp.numOccurence ( aol::NumberTrait<RealType>::Inf );

          const bool firstTryBetterWRTDomain = ( numOutOfDomainFirstTry < numOutOfDomainSecondTry );
          const bool firstTryBetterWRTEnergy = ( firstTryEnergy < _registrationAlgo.getEnergyOfLastSolution() );
          const bool energyAsAltCriterion = _parser.checkAndGetBool( "energyAsAltCriterion" );

          bool usingAltStartLevel = true;
          cerr << "Using ";
          if ( ( energyAsAltCriterion && firstTryBetterWRTEnergy ) || ( !energyAsAltCriterion && firstTryBetterWRTDomain ) ) {
            cerr << "\"startLevel\"";
            _registrationAlgo.setTransformation ( savedTransformation );
            usingAltStartLevel = false;
          }
          else
            cerr << "\"altStartLevel\"";
          cerr << " result\n";
          cerr << ( firstTryBetterWRTDomain ? "\"startLevel\"" : "\"altStartLevel\"" ) << " has higher overlap\n";
          cerr << ( firstTryBetterWRTEnergy ? "\"startLevel\"" : "\"altStartLevel\"" ) << " has lower energy (" << firstTryEnergy << " vs. " << _registrationAlgo.getEnergyOfLastSolution() << ")\n";

          if ( i == 0 )
            energiesFile << i << " " << ( usingAltStartLevel ? _registrationAlgo.getEnergyOfLastSolution () : firstTryEnergy ) << endl;

          _registrationAlgo.setLevel ( _registrationAlgo.getMaxGridDepth() );
          _registrationAlgo.saveTransformation ( aol::strprintf ( "%s/deformation-%03d%s", _parser.getStringExpandTilde ( "saveDirectory" ).c_str(), i, _registrationAlgo.getDeformationFileNameSuffix().c_str() ).c_str() );
        }
        // The first deformation will not be refined, so save its energy directly.
        else if ( i == 0 )
          energiesFile << i << " " << _registrationAlgo.getEnergyOfLastSolution () << endl;

        if ( i == 0 )
          writeCurrentDefParams ( i, defParamsFile );
      }
      else {
        if ( _useAltStartLevel )
          _registrationAlgo.loadTransformation ( aol::strprintf ( "%s../stage1/deformation-%03d%s", _parser.getStringExpandTilde ( "saveDirectory" ).c_str(), i-1, _registrationAlgo.getDeformationFileNameSuffix().c_str() ).c_str() );
        else
          _registrationAlgo.loadTransformation ( aol::strprintf ( "%s../stage1/%d/deformation_%02d%s", _parser.getStringExpandTilde ( "saveDirectory" ).c_str(), i-1, _registrationAlgo.getMaxGridDepth(), _registrationAlgo.getDeformationFileNameSuffix().c_str() ).c_str() );
      }

      if ( _parser.checkAndGetBool ( "dontAccumulateDeformation" ) == false ) {
        TransformationDOFType currentTransformation ( _registrationAlgo.getTransformationDOFInitializer() );
        _registrationAlgo.getTransformation ( currentTransformation );
        TransformationDOFType tempTransformation ( accumulatedTransformation );
        _registrationAlgo.composeDeformations ( currentTransformation, tempTransformation, accumulatedTransformation );

        const int refineStartLevel = _parser.getInt ( "refineStartLevel" );
        const int refineStopLevel = _parser.getIntOrDefault( "refineStopLevel", _registrationAlgo.getMaxGridDepth() );

        // Now match the current template to the first template aka the reference.
        if ( i > 0 ) {
          // Start with the accumulated displacement.
          _registrationAlgo.setTransformation ( accumulatedTransformation );

          if ( reverseRoles )
            _registrationAlgo.setTemplate ( reference );
          else
            _registrationAlgo.setReference ( reference );

          _registrationAlgo.makeAndSetSaveDirectory ( aol::strprintf ( "%s%d-r/", _parser.getStringExpandTilde ( "saveDirectory" ).c_str(), i ).c_str() );
          _registrationAlgo.solveAndProlongToMaxDepth ( refineStartLevel, refineStopLevel );
          // In case we didn't calculate the deformation on the finest level, save the prolongated solution.
          if ( refineStopLevel < _registrationAlgo.getMaxGridDepth() ) {
            _registrationAlgo.saveTransformation ( aol::strprintf ( "%s/deformation_%02d%s", _registrationAlgo.getSaveDirectory(), _registrationAlgo.getMaxGridDepth(), _registrationAlgo.getDeformationFileNameSuffix().c_str() ).c_str() );
          }

          // Save the energy of the refined deformation.
          energiesFile << i << " " << _registrationAlgo.getEnergyOfLastSolution () << endl;
          writeCurrentDefParams ( i, defParamsFile );

          // Store the refined transformation from the current template to the reference.
          _registrationAlgo.getTransformation ( accumulatedTransformation );
        }

        if ( calcInverseDeformation ) {
          const ArrayType savedTemplate ( _registrationAlgo.getTemplImageReference(), aol::DEEP_COPY );
          _registrationAlgo.setTemplate ( _registrationAlgo.getRefImageReference() );
          _registrationAlgo.setReference ( savedTemplate );
          TransformationDOFType inverseTransformation ( _registrationAlgo.getTransformationDOFInitializer() );
          qc::approxInverseDeformation<ConfiguratorType> ( _registrationAlgo.getInitializerRef(), accumulatedTransformation, inverseTransformation );
          _registrationAlgo.setTransformation ( inverseTransformation );
          _registrationAlgo.makeAndSetSaveDirectory ( aol::strprintf ( "%s%d-r-inv/", _parser.getStringExpandTilde ( "saveDirectory" ).c_str(), i ).c_str() );
          _registrationAlgo.solveAndProlongToMaxDepth ( refineStartLevel, refineStopLevel );
          _registrationAlgo.setReference ( _registrationAlgo.getTemplImageReference() );
          _registrationAlgo.setTemplate ( savedTemplate );
        }

        if ( reverseRoles )
          _registrationAlgo.setTemplate ( _registrationAlgo.getRefImageReference() );
        else
          _registrationAlgo.setReference ( _registrationAlgo.getTemplImageReference() );
      }
    }

    // We temporarily messed with the save directory. Restore the original value.
    _registrationAlgo.setSaveDirectory ( origSaveDir.c_str() );
  }

  void saveAverageMedianAndNumSamples ( const ArrayType &Average, const ArrayType &Median, const ArrayType &NumSamples, const int Iter = -1 ) const {
    qc::DefaultArraySaver<RealType, Dim> saver ( true, false, false, false );
    saver.setSaveDirectory ( _parser.getStringExpandTilde ( "saveDirectory" ).c_str() );
    saver.saveStep ( Average, Iter, "average" );
    saver.saveStep ( Median, Iter, "median" );
    saver.saveStep ( NumSamples, Iter, "numSamples" );

    saveAverageMedianAndNumSamplesHelper ( Median, saver, Iter );
  }

private:
  template <int Dim>
  void saveAverageMedianAndNumSamplesHelper ( const typename qc::ScalarArrayTrait<RealType, Dim>::ArrayType &/*Median*/, qc::DefaultArraySaver<RealType, Dim> &/*Saver*/, const int /*Iter*/ ) const {
  }

  void saveAverageMedianAndNumSamplesHelper ( const qc::ScalarArray<RealType, qc::QC_2D> &Median, qc::DefaultArraySaver<RealType, qc::QC_2D> &Saver, const int Iter ) const {
    if ( _parser.checkAndGetBool ( "calcIQ" ) ) {
      im::IQFactor<RealType> IQ ( Median,
                                  _parser.getInt ( "IQPatchSize" ),
                                  _parser.getReal<RealType> ( "IQThreshold" ),
                                  _parser.getReal<RealType> ( "IQOnePixelInAngstrom" ),
                                  true, false );
      IQ.plotIQFactors ( Saver.createSaveName ( "", "", Iter, "median-IQplot" ).c_str() );
    }
    if ( _parser.checkAndGetBool ( "calcPrecision" ) ) {
      Saver.setSaveName ( "median" );
      im::AtomFinder<RealType> atomFinder;
      aol::MultiVector<RealType> atomPositions;
      aol::MultiVector<RealType> gaussianParameters;
      // A proper choice of gamma is critical for a proper segmentation,
      // which in turn is the basis for the number of atoms.
      atomFinder.setGamma ( _parser.getReal<RealType> ( "precSegGamma" ) );
      atomFinder.getSingleAtomPositions ( atomPositions, gaussianParameters, Median );
      atomFinder.getSegmentedConstRef ( ).save ( Saver.createSaveName ( "", "precCenterMask.pgm", Iter ).c_str() );

      const int numCenters = atomPositions.numComponents();
      aol::MultiVector<RealType> sigma ( 2, numCenters );
      for ( int i = 0; i < numCenters; ++i ) {
        sigma[0][i] = gaussianParameters[i][3]; // Width
        sigma[1][i] = gaussianParameters[i][4]; // Height
      }
      const aol::Vec2<RealType> sigmaMean ( sigma[0].getMeanValue(), sigma[1].getMeanValue() );
      const aol::Vec2<RealType> sigmaStdDev ( sigma[0].getStdDev(), sigma[1].getStdDev() );
      cerr << "Sigma_x = " << sigmaMean[0] << " (+-" << sigmaStdDev[0] << ")" << endl;
      cerr << "Sigma_y = " << sigmaMean[1] << " (+-" << sigmaStdDev[1] << ")" << endl;

      aol::Plotter<RealType> plotter;
      plotter.set_outfile_base_name( Saver.createSaveName ( "", "centers", Iter ) );
      aol::PlotDataFileHandler<RealType> plotHandler;
      plotHandler.generateCurvePlot( atomPositions, true );
      plotter.addPlotCommand( plotHandler.getDataFileNames()[0].c_str(), plotHandler.getDataFilePlotTypes()[0], -2 );
      plotter.genPlot( aol::GNUPLOT_PDF );

      cerr << "median - ";
      im::analyzePrecisionFromCenters<RealType> ( atomPositions, _parser, Saver.createSaveName ( "", ".png", Iter ).c_str(), Saver.createSaveName ( "", "prec" ).c_str() );
      // Potentially output more information about the fit.
      // << sigmaMean[0] << " " << sigmaStdDev[0] << " " << sigmaMean[1] << " " << sigmaStdDev[1] << endl;
    }
  }

  template <int Dim>
  void saveNamedDeformedTemplateHelper ( const typename qc::ScalarArrayTrait<RealType, Dim>::ArrayType &DefTemplateArray, const qc::DefaultArraySaver<RealType, Dim> &Saver, const char *TemplateFileName ) const {
    Saver.saveStep ( DefTemplateArray, -1, aol::getBaseFileName( TemplateFileName ).c_str() );
  }

  void saveNamedDeformedTemplateHelper ( const qc::ScalarArray<RealType, qc::QC_2D>  &DefTemplateArray, const qc::DefaultArraySaver<RealType, qc::QC_2D> &Saver, const char *TemplateFileName ) const {
    if ( _parser.checkAndGetBool ( "saveNamedDeformedDMXTemplatesAsDMX" )
        && ( aol::fileNameEndsWith( TemplateFileName, "dm3" ) || aol::fileNameEndsWith( TemplateFileName, "dm4" ) ) ) {
      qc::DM3Reader dmreader( TemplateFileName );
      dmreader.saveQuocDataInDM3Container ( DefTemplateArray, Saver.createSaveName ( "", "", -1, aol::getBaseFileName( TemplateFileName ).c_str() ) );
    }
    else
      Saver.saveStep ( DefTemplateArray, -1, aol::getBaseFileName( TemplateFileName ).c_str() );
  }
public:

  void saveNamedDeformedTemplate ( const char *TemplateFileName, const ArrayType &DefTemplateArray ) const {
    qc::DefaultArraySaver<RealType, Dim> saver ( false, false, true );
    saver.setSaveDirectory ( _parser.getStringExpandTilde ( "saveDirectory" ).c_str() );
    saveNamedDeformedTemplateHelper ( DefTemplateArray, saver, TemplateFileName );
  }

  void averageSeries ( const char *InputDirectory ) const {
    std::vector<std::string> templateFileNames;
    aol::createTemplateFileNameList ( _parser, templateFileNames );
    const int numTemplateImages = static_cast<int> ( templateFileNames.size() );

    ArrayType curTemplate ( _registrationAlgo.getInitializerRef() );

    ArrayType average ( _registrationAlgo.getInitializerRef() );
    ArrayType median ( _registrationAlgo.getInitializerRef() );
    // Keep track of how many samples were available when averaging a certain pixel.
    ArrayType numSamples ( _registrationAlgo.getInitializerRef() );

    if ( _parser.checkAndGetBool ( "reverseRolesInSeriesMatching" ) == false ) {
      // In stage one, the first frame is used as reference and not deformed. Nevetheless,
      // save it if the user wants to save the deformed templates.
      if ( ( _stage == 1 ) && _parser.checkAndGetBool ( "saveNamedDeformedTemplates" ) ) {
        _registrationAlgo.loadAndPrepareImage ( _parser.getStringExpandTilde ( "reference" ).c_str(), curTemplate, true, false, true );
        saveNamedDeformedTemplate ( _parser.getStringExpandTilde ( "reference" ).c_str(), curTemplate );
      }

      aol::MultiVector<RealType> deformedTemplates ( numTemplateImages, _registrationAlgo.getInitializerRef().getNumberOfNodes() );
      aol::MultiVector<RealType> calculatedDeformedTemplates;

      const int averageSaveIncrement = _parser.hasVariable ( "averageSaveIncrement" ) ? _parser.getInt ( "averageSaveIncrement" ) : numTemplateImages;
      if ( averageSaveIncrement <= 0 )
        throw aol::Exception ( "'averageSaveIncrement' must be bigger than zero.", __FILE__, __LINE__ );

      TransformationDOFType reducedPhi;
      if ( _reduceDeformations )
        reducedPhi.load ( aol::strprintf ( "%sreduceDef%s", getSaveDirectory(), _registrationAlgo.getDeformationFileNameSuffix().c_str() ).c_str() );

      // First generate all deformed template images.
      for ( int i = 0; i < numTemplateImages; ++i ) {
        _registrationAlgo.loadAndPrepareImage ( templateFileNames[i].c_str(), curTemplate, true, false, true );
        const string defFileNameBase = createDeformationBaseFileName ( InputDirectory, i );
        if ( _reduceDeformations ) {
          _registrationAlgo.applyTransformation ( reducedPhi, curTemplate, average );
          _registrationAlgo.applySavedTransformation ( defFileNameBase.c_str(), average, deformedTemplates[i], _parser.checkAndGetBool ( "deformForAverageExtendedWithMean" ) ? average.getMeanValue() : aol::NumberTrait<RealType>::Inf, _parser.checkAndGetBool ( "deformForAverageUsingNearestNeighborInterpolation" ) );
        }
        else
          _registrationAlgo.applySavedTransformation ( defFileNameBase.c_str(), curTemplate, deformedTemplates[i], _parser.checkAndGetBool ( "deformForAverageExtendedWithMean" ) ? curTemplate.getMeanValue() : aol::NumberTrait<RealType>::Inf, _parser.checkAndGetBool ( "deformForAverageUsingNearestNeighborInterpolation" ), _parser.getDoubleOrDefault ( "averageXDerivNormThreshold", 0 ) );
        calculatedDeformedTemplates.appendReference ( deformedTemplates[i] );
        if ( _parser.checkAndGetBool ( "saveNamedDeformedTemplates" ) ) {
          ArrayType defTemplateArray ( _registrationAlgo.getInitializerRef() );
          _registrationAlgo.applySavedTransformation ( defFileNameBase.c_str(), curTemplate, defTemplateArray, _parser.checkAndGetBool ( "saveNamedDeformedTemplatesExtendedWithMean" ) ? curTemplate.getMeanValue() : 0, _parser.checkAndGetBool ( "saveNamedDeformedTemplatesUsingNearestNeighborInterpolation" ), _parser.getDoubleOrDefault ( "averageXDerivNormThreshold", 0 ) );
          saveNamedDeformedTemplate ( templateFileNames[i].c_str(), defTemplateArray );
        }
        if ( _parser.checkAndGetBool ( "saveNamedDeformations" ) ) {
          TransformationDOFType phi;
          _registrationAlgo.loadTransformationTo ( defFileNameBase.c_str(), phi );
          _registrationAlgo.saveTransformationTo ( phi, ( _parser.getStringExpandTilde ( "saveDirectory" ) + aol::getBaseFileName( templateFileNames[i] ) + "_def" ).c_str() );
        }

        if ( ( calculatedDeformedTemplates.numComponents() == numTemplateImages ) || ( ( calculatedDeformedTemplates.numComponents() % averageSaveIncrement ) == 0 ) ) {
          const RealType lastTemplateMean = curTemplate.getMeanValue ();

          // Now calculate the average and median.
          for ( int j = 0; j < average.size(); ++j ) {
            average[j] = 0;
            int numSamplesInDomain = 0;
            for ( int i = 0; i < calculatedDeformedTemplates.numComponents(); ++i ) {
              if ( calculatedDeformedTemplates[i][j] != aol::NumberTrait<RealType>::Inf ) {
                average[j] += calculatedDeformedTemplates[i][j];
                ++numSamplesInDomain;
              }
            }
            if ( numSamplesInDomain > 0 )
              average[j] /= numSamplesInDomain;
            else
              average[j] = lastTemplateMean;

            numSamples[j] = numSamplesInDomain;
          }
          
          if ( _parser.checkAndGetBool ( "performNoiseAnalysis" ) ) {
#ifdef USE_INTERNAL_MODULES_IMGPROC
            if ( !_parser.checkAndGetBool ( "deformForAverageUsingNearestNeighborInterpolation" ) ) throw aol::Exception ( "Noise analysis requires nearest neighbor interpolation!", __FILE__, __LINE__ );
            aol::VectorContainer<aol::Vector<RealType> > samples ( calculatedDeformedTemplates[0].size ( ) );
            for ( int j = 0; j < calculatedDeformedTemplates[0].size(); ++j )
              for ( int i = 0; i < calculatedDeformedTemplates.numComponents(); ++i )
                if ( aol::isFinite<RealType> ( calculatedDeformedTemplates[i][j] ) ) samples[j].pushBack ( calculatedDeformedTemplates[i][j] );
            RealType a = 0, b = 0;
            im::NoiseAnalyzer<RealType>::getCMOSVarianceParams ( samples, a, b, false, _parser.getStringExpandTilde ( "saveDirectory" ), calculatedDeformedTemplates.numComponents() );
            std::cerr << "Noise analysis: a = " << a << "; b = " << b << std::endl;
#else
            throw aol::Exception ( "Option \"performNoiseAnalysis\" needs INTERNAL_MODULES_IMGPROC", __FILE__, __LINE__ );
#endif
          }

          calculatedDeformedTemplates.getMedianVecOverComponents ( median, lastTemplateMean );
          saveAverageMedianAndNumSamples ( average, median, numSamples, ( calculatedDeformedTemplates.numComponents() != numTemplateImages ) ? calculatedDeformedTemplates.numComponents() : -1 );
        }
      }
      if ( numTemplateImages == 0 ) {
        _registrationAlgo.loadAndPrepareImage ( _parser.getStringExpandTilde ( "reference" ).c_str(), curTemplate, true, false, true );
        numSamples.setAll ( 1 );
        saveAverageMedianAndNumSamples ( curTemplate, curTemplate, numSamples, -1 );
      }
    }
    else {
      qc::MultiArray<RealType, Dim> phi ( _registrationAlgo.getInitializerRef() );
      qc::GridSize<Dim> targetSize ( _registrationAlgo.getInitializerRef() );
      const int factor = _parser.getIntOrDefault ( "reverseRolesSRFactor", 1 );
      targetSize *= factor;
      const qc::GridStructure targetGrid ( targetSize );
      if ( factor != 1 ) {
        numSamples.reallocate ( targetGrid );
        average.reallocate ( targetGrid );
        median.reallocate ( targetGrid );
      }
      qc::AArray<aol::Vector<RealType>, Dim> samples ( targetGrid );
      qc::AArray<aol::Vector<RealType>, Dim> weights ( targetGrid );
      const RealType postionDistanceThreshold = _parser.hasVariable ( "reverseRolesPostionDistanceThreshold" ) ? _parser.getReal<RealType> ( "reverseRolesPostionDistanceThreshold" ) : 1;
      const bool reverseRolesWeightMedianAndMean = _parser.checkAndGetBool ( "reverseRolesWeightMedianAndMean" );
      // The code used to calculate middlePointPos makes the assumption that the distance
      // from exactDeformedPos to deformedPos is big enough.
      if ( postionDistanceThreshold < 0.251 )
        throw aol::Exception ( "SeriesMatching::averageSeries: Invalid reverseRolesPostionDistanceThreshold value", __FILE__, __LINE__ );

      aol::ProgressBar<> progressBar ( "Sampling data" );
      progressBar.start ( numTemplateImages );
      progressBar.display ( cerr );
      // First generate all deformed template images.
      for ( int i = 0; i < numTemplateImages; ++i ) {
        _registrationAlgo.loadAndPrepareImage ( templateFileNames[i].c_str(), curTemplate, true );
        const string defFileNameBase = createDeformationFileName ( InputDirectory, i );
        phi.load ( defFileNameBase.c_str() );

        for ( qc::RectangularIterator<Dim> it ( _registrationAlgo.getInitializerRef() ); it.notAtEnd(); ++it ) {
          typename ConfiguratorType::VecType exactDeformedPos;
          qc::CoordType deformedPos;
          typename ConfiguratorType::VecType distToExactDeformedPos;
          for ( int j = 0; j < Dim; ++j ) {
            exactDeformedPos[j] = ( (*it)[j] + phi[j].get ( *it ) / _registrationAlgo.getInitializerRef().H() ) * factor;
            deformedPos[j] = reverseRolesWeightMedianAndMean ? static_cast<short> ( std::floor ( exactDeformedPos[j] ) ) : aol::Rint ( exactDeformedPos[j] );
            distToExactDeformedPos[j] = exactDeformedPos[j] - deformedPos[j];
          }
          if ( reverseRolesWeightMedianAndMean ) {
            // The following only works in 2D.
            for ( int k = 0; k <= 1; ++k ) {
              const RealType xWeight = ( k == 0 ) ? ( 1 - distToExactDeformedPos[0] ) : ( distToExactDeformedPos[0] );
              if ( xWeight > 0 ) {
                qc::CoordType node;
                node[0] = deformedPos[0] + k;
                for ( int l = 0; l <= 1; ++l ) {
                  const RealType yWeight = ( l == 0 ) ? ( 1 - distToExactDeformedPos[1] ) : ( distToExactDeformedPos[1] );
                  if ( yWeight > 0 ) {
                    node[1] = deformedPos[1] + l;
                    if ( targetGrid.isAdmissibleNode ( node ) )
                    {
                      samples.getRef ( node ).pushBack ( curTemplate.get ( *it ) );
                      weights.getRef ( node ).pushBack ( xWeight * yWeight );
                    }
                  }
                }
              }
            }
          }
          else if ( distToExactDeformedPos.getMinAbsValue() > postionDistanceThreshold ) {
            typename ConfiguratorType::VecType middlePointPos;
            for ( int j = 0; j < Dim; ++j )
              middlePointPos[j] = aol::Rint ( 2 * distToExactDeformedPos[j] ) * 0.5 + deformedPos[j];
            // The following only works in 2D.
            for ( int k = -1; k <= 1; k = k + 2 ) {
              deformedPos[0] = aol::Rint ( middlePointPos[0] + k * 0.5 );
              for ( int l = -1; l <= 1; l = l + 2 ) {
                deformedPos[1] = aol::Rint ( middlePointPos[1] + l * 0.5 );
                if ( targetGrid.isAdmissibleNode ( deformedPos ) )
                  samples.getRef ( deformedPos ).pushBack ( curTemplate.get ( *it ) );
              }
            }
          }
          else if ( targetGrid.isAdmissibleNode ( deformedPos ) ) {
            samples.getRef ( deformedPos ).pushBack ( curTemplate.get ( *it ) );
          }
        }
        progressBar++;
      }
      progressBar.finish();

      // Take the value of the reference image as sample for those pixels where we don't have any samples yet.
      {
        ArrayType referenceArray ( _registrationAlgo.getInitializerRef() );
        _registrationAlgo.loadAndPrepareImage ( _parser.getStringExpandTilde ( "reference" ).c_str(), referenceArray, true );

        for ( qc::RectangularIterator<Dim> it ( samples ); it.notAtEnd(); ++it ) {
          if ( samples.getRef ( *it ).size() == 0 ) {
            if ( factor == 1 )
              samples.getRef ( *it ).pushBack ( referenceArray.get ( *it ) );
            else {
              aol::Vec<Dim, RealType> pos;
              for ( int j = 0; j < Dim; ++j )
                pos[j] = aol::Min ( (*it)[j] / static_cast<RealType> ( factor ), static_cast<RealType> ( referenceArray.getSize()[j] - 1 ) );

              samples.getRef ( *it ).pushBack ( referenceArray.interpolate ( pos ) );
            }
            if ( reverseRolesWeightMedianAndMean )
              weights.getRef ( *it ).pushBack ( 1 );
          }
        }
      }

      for ( int i = 0; i < average.size(); ++i ) {
        numSamples[i] = samples[i].size();
        average[i] = reverseRolesWeightMedianAndMean ? samples[i].getWeightedMeanValue ( weights[i] ) : samples[i].getMeanValue();
        median[i] = reverseRolesWeightMedianAndMean ? samples[i].getWeightedMedianValue( weights[i] ) : samples[i].getMedianValue();
      }

      saveAverageMedianAndNumSamples ( average, median, numSamples );
    }
  }

  void matchAndAverageSeries ( ) {
    matchSeries();
    if ( _reduceDeformations )
      reduceDeformations ( _registrationAlgo.getSaveDirectory() );
    averageSeries ( _registrationAlgo.getSaveDirectory() );
  }

private:
  template <qc::Dimension Dim>
  void analyzeDeformationsHelper ( const qc::ScalarArray<RealType, Dim> &/*Temp*/, const char * /*InputDirectory*/ ) const {
    throw aol::UnimplementedCodeException ( "Dimension independent analyzeDeformationsHelper not implemented.", __FILE__, __LINE__ );
  }

  void analyzeDeformationsHelper ( qc::ScalarArray<RealType, qc::QC_2D> &Temp, const char *InputDirectory ) const {
    std::vector<std::string> templateFileNames;
    aol::createTemplateFileNameList ( _parser, templateFileNames );
    const int numDeformations = static_cast<int> ( templateFileNames.size() );

    qc::MultiArray<RealType, Dim> phi ( _registrationAlgo.getInitializerRef() );

    aol::RandomAccessContainer<aol::Vec<Dim, RealType> > translations ( numDeformations );

    RealType averageNonRigidComponentMax = 0;
    RealType averageNonRigidComponent = 0;

    for ( int i = 0; i < numDeformations; ++i ) {
      phi.load ( createDeformationFileName ( InputDirectory, i ).c_str() );

      if ( _parser.checkAndGetBool ( "saveDeformedTemplates" ) ) {
        ArrayType curTemplate ( _registrationAlgo.getInitializerRef() );
        _registrationAlgo.loadAndPrepareImage ( templateFileNames[i].c_str(), curTemplate, true );
        _registrationAlgo.applyTransformation ( phi, curTemplate, Temp, aol::NumberTrait<RealType>::Inf, _parser.checkAndGetBool ( "saveDeformedTemplatesUsingNearestNeighborInterpolation" ) );

        qc::DefaultArraySaver<RealType, Dim> saver ( true, true );
        saver.setSaveDirectory ( getSaveDirectory() );
        saver.setSaveTimestepOffset ( 1 );

        if ( _parser.hasVariable ( "enhanceContrastSaturationPercentage" ) )
          saver.setEnhanceContrastSaturationPercentage ( _parser.getDouble ( "enhanceContrastSaturationPercentage" ) );

        Temp.save( saver.createSaveName ( "", ".dat.bz2", i, "defTempl" ).c_str(), qc::PGM_DOUBLE_BINARY );

        if ( Dim == qc::QC_2D ) {
          // Replace Inf by the mean when saving as PNG.
          const RealType curTemplateMean = curTemplate.getMeanValue();
          for ( int j = 0; j < curTemplate.size(); ++j ) {
            if ( aol::isInf ( Temp[j] ) )
              Temp[j] = curTemplateMean;
          }
          saver.saveStep ( Temp, i, "defTempl" );
        }
      }

      translations[i] = phi.getMeanValue();

      for ( int c = 0; c < Dim; ++c )
        phi[c].addToAll ( -translations[i][c] );

      qc::ScalarArray<RealType, Dim> phiNorms ( _registrationAlgo.getInitializerRef() );
      phi.getPointWiseNorm( phiNorms );
      const RealType phiMaxOfPointWiseNorm = phiNorms.getMaxValue();
      averageNonRigidComponentMax += phiMaxOfPointWiseNorm;
      averageNonRigidComponent += phiNorms.getMeanValue();

      cerr << "Translation " << translations[i] << ", non-rigid component max " << phiMaxOfPointWiseNorm / translations[i].norm() * 100 << "%, " << phiMaxOfPointWiseNorm / _registrationAlgo.getInitializerRef().H() << " pixels" <<endl;

      if ( _parser.checkAndGetBool ( "saveXDerivativeOfDef" ) ) {
        qc::ScalarArrayHelper<ArrayType>::getPointWiseXDerivNorm ( phi, Temp );
        qc::DefaultArraySaver<RealType, Dim> saver ( true, false, true );
        saver.setSaveDirectory ( getSaveDirectory() );
        saver.setSaveTimestepOffset ( 1 );
        saver.saveStep ( Temp, i, "def" );
      }
    }

    averageNonRigidComponentMax /= numDeformations;
    averageNonRigidComponent /= numDeformations;

    RealType averageTranslation = 0;
    for ( int i = 1; i < numDeformations; ++i ) {
      averageTranslation += ( translations[i] - translations[i-1] ).norm();
    }
    averageTranslation /= ( numDeformations - 1 );
    cerr << "average translation = " << averageTranslation << ", " << averageTranslation / _registrationAlgo.getInitializerRef().H() << " pixels" << endl;
    cerr << "average non-rigid component = " << averageNonRigidComponent << ", " << averageNonRigidComponent / _registrationAlgo.getInitializerRef().H() << " pixels" << endl;
    cerr << "ratio  = " << averageNonRigidComponent / averageTranslation << endl;
    cerr << "average non-rigid component max = " << averageNonRigidComponentMax << ", " << averageNonRigidComponentMax / _registrationAlgo.getInitializerRef().H() << " pixels" << endl;
    cerr << "ratio  = " << averageNonRigidComponentMax / averageTranslation << endl;

    aol::PlotDataFileHandler<RealType> plotDataFileHandler;
    plotDataFileHandler.generateCurvePlot ( translations );
    aol::Plotter<RealType> plotter;
    plotter.addPlotCommandsFromHandler( plotDataFileHandler );
    plotter.set_outfile_base_name ( aol::strprintf ( "%strans", getSaveDirectory() ).c_str() );
    plotter.genPlot( aol::GNUPLOT_PNG );

    plotDataFileHandler.clear();
    plotter.clearAdditionalPlotCommands();
    aol::Vector<RealType> translationNorms ( numDeformations );
    for ( int i = 0; i < translations.size(); ++i )
      translationNorms[i] = translations[i].norm();
    plotDataFileHandler.generateFunctionPlot ( translationNorms, 1, translations.size() );
    plotter.addPlotCommandsFromHandler( plotDataFileHandler );
    plotter.set_outfile_base_name ( aol::strprintf ( "%snorms", getSaveDirectory() ).c_str() );
    plotter.genPlot( aol::GNUPLOT_PNG );
  }

public:
  void analyzeDeformations ( const char *InputDirectory ) {
    ArrayType temp ( _registrationAlgo.getInitializerRef() );
    analyzeDeformationsHelper ( temp, InputDirectory );
  }

  void reduceDeformations ( const char *InputDirectory ) {
    const int numDeformations = getNumTemplates();

    aol::RandomAccessContainer< qc::MultiArray<RealType, Dim> > deformations ( numDeformations, _registrationAlgo.getInitializerRef() );
    aol::RandomAccessContainer< const qc::ConsistencyEnergyOp<ConfiguratorType> > energies;
    aol::RandomAccessContainer< const qc::VariationOfConsistencyEnergyOp<ConfiguratorType> > derivatives;

    aol::LinCombOp<aol::MultiVector<RealType>, aol::Scalar<RealType> > E;
    aol::LinCombOp<aol::MultiVector<RealType> > DE;

    for ( int i = 0; i < numDeformations; ++i ) {
      deformations[i].load ( createDeformationFileName ( InputDirectory, i ).c_str() );
      energies.constructDatumAndPushBack ( _registrationAlgo.getInitializerRef(), deformations[i] );
      E.appendReference ( energies[i] );
      derivatives.constructDatumAndPushBack ( _registrationAlgo.getInitializerRef(), deformations[i] );
      DE.appendReference ( derivatives[i] );
    }

    typedef typename RegistrationType::GradientDescentType GDType;
    // typedef aol::H1GradientDescent<ConfiguratorType, aol::MultiVector<RealType>, typename qc::MultilevelArrayTrait<RealType, typename ConfiguratorType::InitType>::LinSmoothType > GDType;
    GDType solver ( _registrationAlgo.getInitializerRef(), E, DE, _parser.getIntOrDefault ( "maxGDIterations", 1000 ), 1, _parser.getDoubleOrDefault ( "stopEpsilon", 5e-7 ) );
    solver.setConfigurationFlags ( GDType::USE_NONLINEAR_CG|GDType::LOG_GRADIENT_NORM_AT_OLD_POSITION|GDType::USE_GRADIENT_BASED_STOPPING );

    qc::MultiArray<RealType, Dim> phi ( _registrationAlgo.getInitializerRef() );
    solver.applySingle ( phi );
    phi.save ( aol::strprintf ( "%sreduceDef%s", getSaveDirectory(), _registrationAlgo.getDeformationFileNameSuffix().c_str() ).c_str(), qc::PGM_DOUBLE_BINARY );
  }

  void applyDeformation ( ) {
    std::vector<std::string> templateFileNames;
    aol::createTemplateFileNameList ( _parser, templateFileNames );
    const int numTemplateImages = static_cast<int> ( templateFileNames.size() );

    if ( _stage == 1 ) {
      std::vector<std::string>::iterator it;
      it = templateFileNames.begin();
      it = templateFileNames.insert ( it , _parser.getStringExpandTilde ( "reference" ) );
    }

    ArrayType curTemplate ( _registrationAlgo.getInitializerRef() );
    ArrayType defTemplateArray ( _registrationAlgo.getInitializerRef() );
    const int numDeformations = _parser.getInt ( "numDeformations" );
    TransformationDOFType deformationDofs ( _registrationAlgo.getTransformationDOFInitializer() );
    _registrationAlgo.loadTransformationTo ( aol::strprintf ( _parser.getString ( "deformationBaseNamePattern" ).c_str(), numDeformations - 1, numDeformations - 1 ).c_str(), deformationDofs );

    for ( int i = numDeformations - 2; i >= 0; --i ) {
      _registrationAlgo.loadTransformation ( ( aol::strprintf ( _parser.getString ( "deformationBaseNamePattern" ).c_str(), i, i ) + _registrationAlgo.getDeformationFileNameSuffix() ).c_str() );
      _registrationAlgo.addTransformationTo ( deformationDofs );
    }

    for ( int i = 0; i < numTemplateImages; ++i ) {
      _registrationAlgo.loadAndPrepareImage ( templateFileNames[i].c_str(), curTemplate, true, false, true );

      _registrationAlgo.applyTransformation ( deformationDofs, curTemplate, defTemplateArray, _parser.checkAndGetBool ( "saveNamedDeformedTemplatesExtendedWithMean" ) ? curTemplate.getMeanValue() : 0, _parser.checkAndGetBool ( "saveNamedDeformedTemplatesUsingNearestNeighborInterpolation" ) );
      saveNamedDeformedTemplate ( templateFileNames[i].c_str(), defTemplateArray );
    }
  }

  void discardBrokenFramesFromSeries ( ) const {
    std::vector<std::string> templateFileNames;
    aol::createTemplateFileNameList ( _parser, templateFileNames );
    const int numTemplateImages = templateFileNames.size();

    ArrayType currentFrame ( templateFileNames[0] );
    const aol::Vector<RealType> lastLineOfCurrentFrame ( currentFrame.getData(), currentFrame.getNumX(), aol::FLAT_COPY );

    int numGoodTemplates = 1;
    currentFrame.save ( aol::strprintf ( "%sgood-%03d.dat.bz2", getSaveDirectory(), numGoodTemplates-1 ).c_str(), qc::PGM_DOUBLE_BINARY );

    for ( int frameNum = 1; frameNum < numTemplateImages; ++frameNum ) {
      const ArrayType nextFrame ( templateFileNames[frameNum] );
      const aol::Vector<RealType> lastLineOfNextFrame ( nextFrame.getData(), nextFrame.getNumX(), aol::FLAT_COPY );

      if ( lastLineOfCurrentFrame != lastLineOfNextFrame ) {
        currentFrame = nextFrame;
        ++numGoodTemplates;
        currentFrame.save ( aol::strprintf ( "%sgood-%03d.dat.bz2", getSaveDirectory(), numGoodTemplates-1 ).c_str(), qc::PGM_DOUBLE_BINARY );
      }
    }
  }

  void calcMedianOfUnregistredTemplates ( ) const {
    std::vector<std::string> templateFileNames;
    aol::createTemplateFileNameList ( _parser, templateFileNames );

    const int numTemplateImages = templateFileNames.size();

    ArrayType median ( _registrationAlgo.getInitializerRef() );

    aol::MultiVector<RealType> templates ( numTemplateImages, _registrationAlgo.getInitializerRef().getNumberOfNodes() );

    for ( int i = 0; i < numTemplateImages; ++i ) {
      ArrayType templateArray ( templates[i], _registrationAlgo.getInitializerRef(), aol::FLAT_COPY );
      _registrationAlgo.loadAndPrepareImage ( templateFileNames[i].c_str(), templateArray, true );
    }

    templates.getMedianVecOverComponents ( median );
    median.save ( aol::strprintf ( "%smedianOfTemplates.dat.bz2", getSaveDirectory() ).c_str(), qc::PGM_DOUBLE_BINARY );
  }

  const char* getSaveDirectory () const {
    return _registrationAlgo.getSaveDirectory();
  }

  void doAction ( const ACTION ActionType ) {
    switch ( ActionType ) {
      case MATCH_AND_AVERAGE_SERIES :
        matchAndAverageSeries();
        break;
      case ONLY_AVERAGE_SERIES:
        averageSeries( getSaveDirectory() );
        break;
      case ANALYZE_DEFORMATTIONS:
        analyzeDeformations( getSaveDirectory() );
        break;
      case APPLY_DEFORMATTION:
        applyDeformation( );
        break;
      default:
        throw aol::Exception ( "SeriesMatching::doAction: Invalid ActionType", __FILE__, __LINE__ );
    } 
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ConfiguratorType1D, typename RegistrationConfiguratorType1D, typename GradientDescentType>
typename ConfiguratorType::RealType findXShift ( const typename ConfiguratorType::InitType &Grid,
                                                 const typename ConfiguratorType::ArrayType &Reference,
                                                 const typename ConfiguratorType::ArrayType &Template,
                                                 const typename ConfiguratorType::RealType Lambda,
                                                 const typename ConfiguratorType::RealType YRegFactor,
                                                 const typename ConfiguratorType::RealType LaplaceWeight,
                                                 qc::ScalarArray<typename ConfiguratorType::RealType, qc::QC_2D>&PhiX,
                                                 const char *EnergyPlotFile = NULL,
                                                 const int MaxGDIterations = 500 ) {

  typedef typename ConfiguratorType::RealType RealType;
  qc::FELineShiftRegisEnergy<ConfiguratorType, ConfiguratorType1D, RegistrationConfiguratorType1D> EData ( Grid, Reference, Template );
  const aol::DerivativeWrapper<RealType, qc::FELineShiftRegisEnergy<ConfiguratorType, ConfiguratorType1D, RegistrationConfiguratorType1D>, aol::Vector < RealType > > DEData ( EData );

  aol::DeleteFlagPointer<aol::Op<aol::Vector<RealType> > > pRegOp;
  const aol::FEOpMixedDerivative<ConfiguratorType> stiffX ( Grid, 0, 0, aol::ONTHEFLY );
  const aol::FEOpMixedDerivative<ConfiguratorType> stiffY ( Grid, 1, 1, aol::ONTHEFLY );

  if ( LaplaceWeight <= 0 ) {
    typename ConfiguratorType::MatrixType *pStiffMat = new typename ConfiguratorType::MatrixType ( Grid );
    stiffX.assembleAddMatrix ( *pStiffMat );
    stiffY.assembleAddMatrix ( *pStiffMat, YRegFactor );
    pRegOp.reset ( pStiffMat, true );
  }
  else {
    aol::SparseMatrix<RealType> tempMat ( Grid );
    qc::assembleLaplaceSquareRegMatrix<ConfiguratorType> ( Grid, tempMat );
    tempMat *= LaplaceWeight;
    stiffX.assembleAddMatrix ( tempMat );
    stiffY.assembleAddMatrix ( tempMat, YRegFactor );

    aol::CSR_Matrix<> regMat ( tempMat );

    pRegOp.reset ( new aol::CSR_Matrix<> ( tempMat ), true );
  }

  const aol::QuadraticFormOp<aol::Vector<RealType> > ERegStiff ( *pRegOp );

  aol::LinCombOp<aol::Vector<RealType>, aol::Scalar<RealType> > E;
  E.appendReference ( EData );
  E.appendReference ( ERegStiff, Lambda );

  aol::LinCombOp<aol::Vector<RealType> > DE;
  DE.appendReference ( DEData );
  DE.appendReference ( *pRegOp, Lambda );

  typedef aol::GradientDescentWithAutomaticFilterWidth<ConfiguratorType, aol::Vector<RealType>, GradientDescentType> GDType;
  GDType gradientDescentSolver ( Grid, E, DE, MaxGDIterations );
  gradientDescentSolver.setConfigurationFlags ( GDType::USE_NONLINEAR_CG );

  std::ofstream out;
  if ( EnergyPlotFile ) {
    out.open ( EnergyPlotFile );
    gradientDescentSolver.setOutStream ( out );
  }

  const RealType scaleFac = static_cast<RealType> ( PhiX.getSize()[0] ) / PhiX.getSize().getMaxValue();
  // The 1D line integrals use h_x to interpret the deformation.
  PhiX /= scaleFac;
  gradientDescentSolver.applySingle ( PhiX );
  PhiX *= scaleFac;
  return gradientDescentSolver.getEnergyAtLastPosition();
}

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ConfiguratorType1D, typename RegistrationConfiguratorType1D>
typename ConfiguratorType::RealType findXShiftYConst ( const typename ConfiguratorType::InitType &Grid,
                                                       const typename ConfiguratorType::ArrayType &Reference,
                                                       const typename ConfiguratorType::ArrayType &Template,
                                                       qc::ScalarArray<typename ConfiguratorType::RealType, qc::QC_2D> &PhiX,
                                                       const char *EnergyPlotFile = NULL,
                                                       const int MaxGDIterations = 500 ) {
  typedef typename ConfiguratorType::RealType RealType;
  typename ConfiguratorType::ArrayType temSmoothed ( Grid );
  
  aol::Vector<RealType> tmp ( Grid.getNumY() );
  PhiX.sumInXDirectionTo ( tmp );
  tmp /= Grid.getNumX();
  RealType energyVal = aol::NumberTrait<RealType>::Inf;
  for ( int i = 3; i >= 0; --i ) {
    
    if ( i > 0 ) {
      qc::ModifiableKernel2d<RealType> xAveraginKernel ( 4*i + 1 );
      const int offset = xAveraginKernel.getOffset();
      for ( int X = -offset; X <= offset; X++ )
        xAveraginKernel.setValue ( X, 0, 1. / ( 2 * offset + 1 ) );
      Template.applyLinearFilterTo ( xAveraginKernel, temSmoothed );
    }
    else
      temSmoothed = Template;
    
    const qc::FELineShiftRegisEnergy<ConfiguratorType, ConfiguratorType1D, RegistrationConfiguratorType1D> EData ( Grid, Reference, temSmoothed );
    const aol::DerivativeWrapper<RealType, qc::FELineShiftRegisEnergy<ConfiguratorType, ConfiguratorType1D, RegistrationConfiguratorType1D>, aol::Vector < RealType > > DEData ( EData );
    
    const qc::PixelWiseToRowWiseOp<ConfiguratorType> E ( Grid, EData, DEData );
    const aol::DerivativeWrapper<RealType, qc::PixelWiseToRowWiseOp<ConfiguratorType>, aol::Vector < RealType > > DE ( E );
    
    typedef aol::GridlessGradientDescent<RealType, aol::Vector<RealType> > GradientDescentType;
    GradientDescentType gdSolver ( Grid, E, DE, MaxGDIterations );
    gdSolver.setConfigurationFlags ( GradientDescentType::USE_NONLINEAR_CG );
    
    std::ofstream out;
    if ( EnergyPlotFile ) {
      out.open ( EnergyPlotFile );
      gdSolver.setOutStream ( out );
    }
    
    gdSolver.applySingle ( tmp );
    energyVal = gdSolver.getEnergyAtLastPosition();
  }
  PhiX.setInXDirectionFrom ( tmp );
  return energyVal;
}

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ConfiguratorType1D, typename _RegistrationConfiguratorType1D, typename _GradientDescentType = aol::H1GradientDescent<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType>, typename qc::MultilevelArrayTrait<typename ConfiguratorType::RealType, typename ConfiguratorType::InitType>::LinSmoothType> >
class XShiftRegistrationMultilevelDescent : public qc::RegistrationMultilevelDescentInterface<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef _RegistrationConfiguratorType1D RegistrationConfiguratorType1D;
  typedef _GradientDescentType GradientDescentType;
  static const bool IsParametric = false;

  RealType _energyOfLastSolution;

public:
  XShiftRegistrationMultilevelDescent ( const aol::ParameterParser &Parser ) :
  qc::RegistrationMultilevelDescentInterface<ConfiguratorType> ( Parser ),
  _energyOfLastSolution ( aol::NumberTrait<RealType>::NaN ) {}

  virtual ~XShiftRegistrationMultilevelDescent( ) {}

  void setLevel ( const int Level ) {
    qc::RegistrationMultilevelDescentInterface<ConfiguratorType>::setLevel ( Level );
    this->_transformation.setCurLevel ( Level );
  }

  void prolongate( ) {
    if ( this->_curLevel < this->getMaxGridDepth() ) {
      if ( this->getParserReference().checkAndGetBool ( "resampleInsteadOfProlongateDeformation" ) ) {
        for ( int i = 0; i < this->_transformation.numComponents(); ++i )
          this->_transformation.getArray ( i, this->_curLevel + 1 ).resampleFrom ( this->_transformation.getArray ( i, this->_curLevel ) );
        this->_transformation.setCurLevel ( this->_curLevel + 1 );
      }
      else
        this->_transformation.levProlongate( );
      setLevel ( this->_curLevel + 1 );
    }
  }

  void descentOnCurrentGrid() {
    cerr << "\n--------------------------------------------------------\n";
    cerr << "Registration on level " << this->_curLevel << " started";
    cerr << "\n";
    cerr << "--------------------------------------------------------\n\n";

    qc::ScalarArray<typename ConfiguratorType::RealType, qc::QC_2D>& phiX = this->_transformation.getArray( 0 );

    // To handle saving, create the full deformation (y-component should always be zero) and a registration configurator.
    const qc::MultiArray<RealType, ConfiguratorType::Dim> phi ( this->_transformation, aol::FLAT_COPY );
    qc::BaseRegistrationConfigurator<ConfiguratorType> regisConfig;
    regisConfig.writeInitialRegistration( *this, phi );

    _energyOfLastSolution = findXShift<ConfiguratorType, ConfiguratorType1D, RegistrationConfiguratorType1D, GradientDescentType> ( this->getCurrentGrid(), this->_org_reference[ this->_curLevel ], this->_org_template[ this->_curLevel ], this->getParserReference().getDouble ( "lambda" ), this->getParserReference().getDouble ( "yRegFactor" ), this->getParserReference().getDouble ( "laplaceWeight" ), phiX, aol::strprintf ( "%senergy_%02d.txt", this->getSaveDirectory(), this->_curLevel ).c_str(), this->getParserReference().getIntOrDefault ( "maxGDIterations", 500 ) );

    regisConfig.writeCurrentRegistration( *this, phi );
  }

  RealType getEnergyOfLastSolution ( ) const {
    return _energyOfLastSolution;
  }
};

/**
 * \author Berkels
 */
template <typename _ConfiguratorType, template<class> class RegistrationConfiguratorType,
          template<class> class RegularizationConfiguratorType = qc::DirichletRegularizationConfigurator,
          typename DyadicConfiguratorType = qc::QuocConfiguratorTraitMultiLin<typename _ConfiguratorType::RealType, _ConfiguratorType::Dim, aol::GaussQuadrature<typename _ConfiguratorType::RealType, _ConfiguratorType::Dim,3> >,
          typename DyadicRegistrationMultilevelDescentType = qc::StandardRegistrationMultilevelDescent<DyadicConfiguratorType, RegistrationConfiguratorType<DyadicConfiguratorType>, RegularizationConfiguratorType<DyadicConfiguratorType> >,
          typename _GradientDescentType = aol::H1GradientDescent<_ConfiguratorType, aol::MultiVector<typename _ConfiguratorType::RealType> > >
class NonDyadicRegistrationMultilevelDescent : public qc::RegistrationInterface<_ConfiguratorType> {
public:
  typedef _ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  typedef typename ConfiguratorType::InitType InitType;
  typedef _GradientDescentType GradientDescentType;
  typedef qc::MultiArray<RealType, ConfiguratorType::Dim> TransformationDOFType;
  static const bool IsParametric = false;

  const aol::ParameterParser &_parser;
  aol::ParameterParser _mldParser;
  aol::DeleteFlagPointer<DyadicRegistrationMultilevelDescentType> _pMld;
  qc::MultiArray<RealType, ConfiguratorType::Dim> _phi;
  RealType _energyOfLastSolution;
public:
  NonDyadicRegistrationMultilevelDescent ( const aol::ParameterParser &Parser )
    : qc::RegistrationInterface<ConfiguratorType> ( Parser.getStringExpandTilde ( "reference" ).c_str(), Parser.getStringExpandTilde ( "template" ).c_str(), Parser.getDouble ( "lambda" ), Parser.checkAndGetBool ( "dontNormalizeInputImages" ) ),
      _parser ( Parser ),
      _mldParser ( Parser ),
      _phi ( this->_grid ),
      _energyOfLastSolution ( aol::NumberTrait<RealType>::NaN ) {
    const int depth = Parser.getIntOrDefault ( "dyadicDepth", qc::logBaseTwo ( this->_grid.getSize().getMaxValue() - 1 ) );
    _mldParser.changeVariableValue ( "stopLevel", depth );
    _mldParser.changeVariableValue ( "precisionLevel", depth );
    if ( _mldParser.hasVariable ( "dontResizeOrCropReference" ) )
      _mldParser.changeVariableValue ( "dontResizeOrCropReference", 0 );
    if ( _mldParser.hasVariable ( "resizeInput" ) )
      _mldParser.changeVariableValue ( "resizeInput", 1 );
    else
      _mldParser.addVariable ( "resizeInput", 1 );

    if ( Parser.hasVariable ( "dyadicMinimizationAlgo" ) )
      _mldParser.changeVariableValue ( "minimizationAlgo", Parser.getString ( "dyadicMinimizationAlgo" ).c_str() );

    _pMld.reset ( new DyadicRegistrationMultilevelDescentType ( _mldParser ), true );

    // We still need to load the input at full resolution.
    loadAndPrepareImage ( Parser.getStringExpandTilde ( "reference" ).c_str(), this->getNonConstRefImageReference(), Parser.checkAndGetBool ( "dontNormalizeInputImages" ) );
    loadAndPrepareImage ( Parser.getStringExpandTilde ( "template" ).c_str(), this->getNonConstTemplImageReference(), Parser.checkAndGetBool ( "dontNormalizeInputImages" ) );

    if ( _parser.getInt( "stopLevel" ) != getMaxGridDepth() )
      throw aol::Exception( aol::strprintf ( "NonDyadicRegistrationMultilevelDescent requires the parameter \"stopLevel\" to be set in accordance to the input data (expected %d, got %d).", getMaxGridDepth(), _parser.getInt( "stopLevel" ) ).c_str(), __FILE__, __LINE__);
  }

  virtual ~NonDyadicRegistrationMultilevelDescent () {}

  RealType H() const {
    return this->_grid.H();
  }

  void setLevel ( const int /*Level*/ ) {}

  void loadAndPrepareImage ( const char* Filename, ArrayType &Dest, const bool NoScaling, const bool /*NoResizeOrCrop*/ = false, const bool /*NoSmoothing*/ = false ) const {
    Dest.load ( Filename );
    if ( _parser.hasVariable( "clampRange" ) ) {
      aol::Vec<2, RealType> clampRange = _parser.getRealVec<2, RealType> ( "clampRange" );
      Dest.clamp ( clampRange[0], clampRange[1] );
    }
    if ( _parser.hasVariable ( "arcTanArgScalingFactor" ) ) {
      const RealType scaleFac = _parser.getReal<RealType> ( "arcTanArgScalingFactor" );
      for ( int i = 0; i < Dest.size(); ++i )
        Dest[i] = atan ( scaleFac * Dest[i] );
    }
    if ( NoScaling == false )
      Dest.scaleValuesTo01();
  }

  virtual void loadRefOrTemplate ( const char* Filename, const qc::REGISTRATION_INPUT_TYPE Input, const bool NoScaling ) {
    _pMld->loadRefOrTemplate( Filename, Input, NoScaling );
    loadAndPrepareImage ( Filename, ( Input == qc::REFERENCE ) ? this->getNonConstRefImageReference() : this->getNonConstTemplImageReference(), NoScaling );
  }

  virtual void setTemplate ( const ArrayType &Template ) {
    this->getNonConstTemplImageReference() = Template;
    _pMld->getTemplImageMLAReference()[_pMld->getMaxGridDepth()].resampleFrom ( Template );
    _pMld->getTemplImageMLAReference().levRestrict ( 0, _pMld->getMaxGridDepth() );
  }

  virtual void setReference ( const ArrayType &Reference ) {
    this->getNonConstRefImageReference() = Reference;
    _pMld->getRefImageMLAReference()[_pMld->getMaxGridDepth()].resampleFrom ( Reference );
    _pMld->getRefImageMLAReference().levRestrict ( 0, _pMld->getMaxGridDepth() );
  }

  const InitType& getTransformationDOFInitializer ( ) const {
    return this->getInitializerRef();
  }

  string getDeformationFileNameSuffix ( ) const {
    return "_%d.dat.bz2";
  }

  void setTransformation ( const qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation ) {
    _phi = Transformation;
  }

  void getTransformation ( qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation ) const {
    Transformation = _phi;
  }

  const qc::MultiArray<RealType, ConfiguratorType::Dim> &getTransformationReference ( ) const {
    return _phi;
  }

  qc::MultiArray<RealType, ConfiguratorType::Dim> &getTransformationReference ( ) {
    return _phi;
  }

  RealType getTransformationNorm ( ) const {
    return _phi.norm();
  }

  void addTransformationTo ( qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation ) const {
    Transformation += _phi;
  }

  void setTransformationToZero ( ) {
    _phi.setZero( );
  }

  void setTransformationToTranslation ( const aol::Vec<ConfiguratorType::Dim, RealType> &Translation ) {
    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      _phi[i].setAll ( Translation[i] );
  }

  void composeDeformations ( const qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation1,  const qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation2, qc::MultiArray<RealType, ConfiguratorType::Dim> &Composition ) {
    qc::ConcatAndSmoothlyExtendDeformations<ConfiguratorType> ( Transformation1, Transformation2, this->_grid, Composition );
  }

  void setTransformationToComposition ( const qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation1,  const qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation2 ) {
    qc::ConcatAndSmoothlyExtendDeformations<ConfiguratorType> ( Transformation1, Transformation2, this->_grid, _phi );
  }

  void applyTransformation ( const qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation, const aol::Vector<RealType> &InputImage, aol::Vector<RealType> &DeformedImage, const RealType ExtensionConstant = aol::NumberTrait<RealType>::Inf, const bool NearestNeighborInterpolation = false ) const {
    qc::DeformImage<ConfiguratorType> ( InputImage, this->_grid, DeformedImage, Transformation, true, ExtensionConstant, NearestNeighborInterpolation );
  }

  void applyCurrentTransformation ( const aol::Vector<RealType> &InputImage, aol::Vector<RealType> &DeformedImage ) const {
    applyTransformation ( _phi, InputImage, DeformedImage );
  }

  void applySavedTransformation ( const char *DefBaseName, const aol::Vector<RealType> &InputImage, aol::Vector<RealType> &DeformedImage, const RealType ExtensionConstant = aol::NumberTrait<RealType>::Inf, const bool NearestNeighborInterpolation = false, const RealType xDerivNormThreshold = 0 ) const {
    qc::MultiArray<RealType, ConfiguratorType::Dim> phi ( this->_grid );
    loadTransformationTo ( DefBaseName, phi );
    applyTransformation ( phi, InputImage, DeformedImage, ExtensionConstant, NearestNeighborInterpolation );
    if ( xDerivNormThreshold > 0 ) {
      // \todo Use qc::applyStretchMute instead.
      ArrayType xDerivNorm ( this->_grid );
      qc::ScalarArrayHelper<ArrayType>::getPointWiseXDerivNorm ( phi, xDerivNorm );
      typename qc::BitArray<ConfiguratorType::Dim> maskedPositions ( this->_grid );
      maskedPositions.thresholdFrom ( xDerivNorm, xDerivNormThreshold );
      maskedPositions.erodeByOne();
      maskedPositions.dilateByOne();
      for ( int k = 0; k < xDerivNorm.size(); ++k )
        if ( maskedPositions[k] )
          DeformedImage[k] = ExtensionConstant;
    }
  }

  void loadTransformationTo ( const char *DefBaseName, qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation ) const {
    Transformation.load ( aol::strprintf ( "%s%s", DefBaseName, getDeformationFileNameSuffix().c_str() ).c_str() );
  }

  void saveTransformationTo ( const qc::MultiArray<RealType, ConfiguratorType::Dim> &Transformation, const char *DefBaseName ) const {
    Transformation.save ( aol::strprintf ( "%s%s", DefBaseName, getDeformationFileNameSuffix().c_str() ).c_str(), qc::PGM_DOUBLE_BINARY );
  }

  void saveTransformation ( const char *FileName ) const {
    _phi.save ( FileName, qc::PGM_DOUBLE_BINARY );
  }

  void loadTransformation ( const char *FileName ) {
    qc::MultiArray<RealType, ConfiguratorType::Dim> phi ( this->_grid );
    phi.load ( FileName );
    setTransformation ( phi );
  }

  const char* getSaveDirectory ( ) const {
    return _pMld->getSaveDirectory();
  }

  void setSaveDirectory ( const char *SaveDirectory ) {
    _pMld->setSaveDirectory( SaveDirectory );
  }

  void makeAndSetSaveDirectory ( const char *SaveDirectory ) {
    _pMld->makeAndSetSaveDirectory ( SaveDirectory );
  }

  int getMaxGridDepth( ) const {
    return _pMld->getMaxGridDepth() + ( ( _pMld->getInitializerRef().getSize() != this->_grid.getSize() ) ? 1 : 0 );
  }

  void solve ( const int StartLevel = -1 ) {
    if ( this->getRefImageReference().getSize() != this->getTemplImageReference().getSize() )
    {
      cerr << this->getRefImageReference().getSize() << " != " << this->getTemplImageReference().getSize() << endl;
      throw aol::Exception( "NonDyadicRegistrationMultilevelDescent requires the input images to be of the same size\n", __FILE__, __LINE__);
    }

    typename DyadicRegistrationMultilevelDescentType::TransformationDOFType phiMld ( _pMld->getTransformationDOFInitializer() );
    phiMld.resampleFrom ( _phi );
    // Account for switching between quadratic and non-quadratic grids.
    const int maxWidth = _phi.getSize().getMaxValue();
    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      phiMld[i] /= _phi.getSize()[i] / static_cast<RealType> ( maxWidth );
    _pMld->setTransformation ( phiMld );
    _pMld->solve ( StartLevel );

    if ( _pMld->getMaxGridDepth() != getMaxGridDepth() ) {
      cerr << "\n--------------------------------------------------------\n";
      cerr << "Registration on full resolution started\n";
      cerr << "--------------------------------------------------------\n\n";

      if ( _pMld->getParserReference().checkAndGetBool ( "saveRefAndTempl" ) ) {
        this->getTemplImageReference().save ( aol::strprintf ( "%stemplateFine%s", _pMld->getSaveDirectory(), qc::getDefaultArraySuffix ( ConfiguratorType::Dim ) ).c_str(), qc::PGM_DOUBLE_BINARY );
        this->getRefImageReference().save ( aol::strprintf ( "%sreferenceFine%s", _pMld->getSaveDirectory(), qc::getDefaultArraySuffix ( ConfiguratorType::Dim ) ).c_str(), qc::PGM_DOUBLE_BINARY );
      }
      qc::MultiArray<RealType, ConfiguratorType::Dim> approxPhi ( _pMld->getTransformationReference(), aol::FLAT_COPY );
      _phi.resampleFrom ( approxPhi );
      // Account for switching between quadratic and non-quadratic grids.
      for ( int i = 0; i < ConfiguratorType::Dim; ++i )
        _phi[i] *= _phi.getSize()[i] / static_cast<RealType> ( maxWidth );

      // Only try to refine the deformation if at least one pixel is not NaN in both images.
      aol::Vector<RealType> tmp ( this->getRefImageReference() );
      tmp += this->getTemplImageReference();
      if ( aol::VectorHelper<aol::Vector<RealType> >::numNonNaNs ( tmp ) > 0 )
        _energyOfLastSolution = updateNonDyadicDeformation ( _phi );
      qc::RegistrationStepSaver<ConfiguratorType>
      stepSaver( this->_grid, this->getRefImageReference(),  this->getTemplImageReference(), _parser.getInt ( "checkboxWidth" ) );
      stepSaver.setSaveName ( aol::strprintf ( "_%02d", _pMld->getMaxGridDepth() + 1 ).c_str() );
      stepSaver.setSaveDirectory ( _pMld->getSaveDirectory() );
      stepSaver.setSaveImagesAsDouble ( _pMld->getParserReference().checkAndGetBool ( "saveAsDouble" ) );
      stepSaver.saveStep ( _phi, -1 );
    }
  }

  virtual RealType updateNonDyadicDeformation ( qc::MultiArray<RealType, ConfiguratorType::Dim> &Phi ) const {
    const RegistrationConfiguratorType<ConfiguratorType> regisConfig ( _parser );
    const RegularizationConfiguratorType<ConfiguratorType> regulConfig ( this->_grid, _parser );
    qc::StandardRegistration<ConfiguratorType, ArrayType, RegistrationConfiguratorType<ConfiguratorType>, RegularizationConfiguratorType<ConfiguratorType>, GradientDescentType> stdRegistration ( this->getRefImageReference(), this->getTemplImageReference(), regisConfig, regulConfig, this->_lambda );

    stdRegistration.setSettingsFromParser ( _parser );
    return stdRegistration.findTransformation ( Phi );
  }

  RealType getEnergyOfLastSolution ( ) const {
    return _energyOfLastSolution;
  }

  void solveAndProlongToMaxDepth ( const int StartLevel = -1, const int /*StopLevel*/ = -1 ) {
    solve ( StartLevel );
  }
};

template <typename ConfiguratorType, typename RegistrationType>
void matchSeries ( aol::ParameterParser &Parser, int argc, char **argv ) {
  aol::StopWatch watch;
  watch.start();

  typedef typename ConfiguratorType::RealType RealType;
  const int numExtraStages = Parser.getInt ( "numExtraStages" );
  const RealType extraStageslambdaFactor = Parser.getDouble ( "extraStagesLambdaFactor" );

  // This creates the save directory and dumps the paramters to a file in that dir.
  qc::DefaultArraySaver<RealType, ConfiguratorType::Dim> saver;
  saver.initFromParser ( Parser, true );
  aol::AdditionalOutputToFile addOut ( saver.createSaveName ( "", ".txt", -1, "log" ).c_str() );

  if ( Parser.hasVariable ( "templateNamePattern" ) ) {
    if ( Parser.hasVariable ( "reference" ) )
      throw aol::Exception( "Paramater \"reference\" may not be specified if \"templateNamePattern\" is used.", __FILE__, __LINE__);
    else
    {
      if ( Parser.hasVariable ( "initialAverage" ) == false ) {
        // Use the first of the templates as reference image
        Parser.addVariable ( "reference", aol::strprintf ( Parser.getString ( "templateNamePattern" ).c_str(), Parser.getInt ( "templateNumOffset" ) ).c_str() );
        Parser.changeVariableValue ( "templateNumOffset", aol::strprintf ( "%d", Parser.getInt ( "templateNumOffset" ) + Parser.getInt ( "templateNumStep" ) ).c_str() );
        Parser.changeVariableValue ( "numTemplates", aol::strprintf ( "%d", Parser.getInt ( "numTemplates" ) - 1 ).c_str() );
      }
      else
        Parser.addVariable ( "reference", Parser.getStringExpandTilde ( "initialAverage" ).c_str() );
    }
  }

  // qc::RegistrationMultilevelDescentInterfaceBase expects the parser to have a valid "template" parameter.
  // SeriesMatching itself doesn't need it, so just add the current "reference" value as valid dummy value for "template".
  if ( Parser.hasVariable ( "template" ) == false )
    Parser.addVariable ( "template", Parser.getString ( "reference" ).c_str() );

  aol::ParameterParser parserStage1 ( Parser );
  parserStage1.changeVariableValue ( "saveDirectory", ( string ( saver.getSaveDirectory() ) + "stage1/" ).c_str() );

  SeriesMatching<RegistrationType> mldS1 ( parserStage1 );
  const typename SeriesMatching<RegistrationType>::ACTION actionType = ( argc > 2 ) ? static_cast<typename SeriesMatching<RegistrationType>::ACTION> ( atoi ( argv[2] ) ) : SeriesMatching<RegistrationType>::MATCH_AND_AVERAGE_SERIES;
  if ( Parser.checkAndGetBool ( "skipStage1" ) == false )
    mldS1.doAction( actionType );
  const string averageName = string ( Parser.checkAndGetBool ( "useMedianAsNewTarget" ) ? "median" : "average" ) + qc::getDefaultArraySuffix ( ConfiguratorType::Dim );
  string lastAverageFileName = string ( mldS1.getSaveDirectory() ) + averageName;

  if ( numExtraStages > 0 ) {
    aol::ParameterParser parserStage2 ( Parser );
    if ( Parser.hasVariable ( "initialAverage" ) == false ) {
      parserStage2.changeVariableValue ( "templateNumOffset", aol::strprintf ( "%d", Parser.getInt ( "templateNumOffset" ) - Parser.getInt ( "templateNumStep" ) ).c_str() );
      parserStage2.changeVariableValue ( "numTemplates", aol::strprintf ( "%d", Parser.getInt ( "numTemplates" ) + 1 ).c_str() );
    }
    parserStage2.changeVariableValue ( "lambda", aol::strprintf ( "%f", Parser.getDouble ( "lambda" ) * extraStageslambdaFactor ).c_str() );
    if ( Parser.checkAndGetBool ( "cropInput" ) || Parser.checkAndGetBool ( "resizeInput" ) )
      parserStage2.changeVariableValue ( "dontResizeOrCropReference", 1 );
    if ( Parser.checkAndGetBool ( "useCorrelationToInitTranslation" ) )
      parserStage2.changeVariableValue ( "useCorrelationToInitTranslation", 0 );

    for ( int i = 0; i < numExtraStages; ++i ) {
      aol::ParameterParser parserStageN ( parserStage2 );
      parserStageN.changeVariableValue ( "saveDirectory", aol::strprintf( "%sstage%d/", saver.getSaveDirectory(), i+2 ).c_str() );
      parserStageN.changeVariableValue ( "reference", lastAverageFileName.c_str() );
      parserStageN.addVariable ( "stage", i+2 );

      SeriesMatching<RegistrationType> mldSN ( parserStageN, Parser.checkAndGetBool ( "reuseStage1Results" ) );
      mldSN.doAction( actionType );
      lastAverageFileName = string ( mldSN.getSaveDirectory() ) + averageName;
    }
  }

  watch.stop();
  watch.printReport ( cerr );
}

#endif // __MATCHSERIES_H
