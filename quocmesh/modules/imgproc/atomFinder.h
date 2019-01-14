#ifndef __ATOMFINDER_H
#define __ATOMFINDER_H

#include <aol.h>
#include <connectedComponents.h>
#include <configurators.h>
#include <segmentation.h>
#include <primalDualSegmentation.h>
#include <imageTools.h>
#include <multiArray.h>
#include <parameterParser.h>
#include <scalarArray.h>
#include <bumpFit.h>
#include <regression.h>


namespace im {


template <typename _RealType,
          typename _MatrixType = aol::FullMatrix<_RealType>,
          typename _LinearRegressionType = aol::LinearRegressionQR<_RealType>,
          typename _ScalarPictureType = qc::ScalarArray<_RealType, qc::QC_2D>,
          typename _ColoredPictureType = qc::MultiArray<_RealType, qc::QC_2D, 3> >
class AtomFinder {
  typedef _RealType RealType;
  typedef _MatrixType MatrixType;
  typedef _LinearRegressionType LinearRegressionType;
  typedef _ScalarPictureType PictureType;
  typedef _ColoredPictureType ColoredPictureType;
  typedef qc::MultiArray<RealType, qc::QC_2D, 1> ArrayType;
  typedef qc::RectangularGridConfigurator<RealType, qc::QC_2D, aol::GaussQuadrature<RealType, qc::QC_2D, 3> > ConfiguratorType;
  typedef typename ConfiguratorType::InitType InitType;
  typedef typename qc::ComponentsCollection<RealType>::NonEmptyComponentsIterator NonEmptyComponentsIterator;
  typedef aol::BoxProjector<RealType, aol::Vector<RealType> > ProjectorType;
public:
  static RealType gammaDefault, epsilonDefault;
  static int maxItDefault, numOuterIterationsDefault;
protected:
  std::string _outputDir;
  bool _quietMode, _diskOutput;
  aol::ProgressBar<> *_progressBar;
  bool _catchCtrlC;
  mutable sigfunc _previousCtrlCHandler;
  
  RealType _gamma, _epsilon;
  int _maxIt, _numOuterIterations;
  bool _lmVerbose;
  qc::BitArray<qc::QC_2D> _segmented;
public:
  AtomFinder ( const std::string &OutputDir = "", const bool Quiet = true )
    : _outputDir ( OutputDir ), _quietMode ( Quiet ), _diskOutput ( OutputDir != "" ), _progressBar ( NULL ), _catchCtrlC ( false ),
      _gamma ( gammaDefault ), _epsilon ( epsilonDefault ), _maxIt ( maxItDefault ), _numOuterIterations ( numOuterIterationsDefault ),
      _lmVerbose ( false ) { }
  
  void getAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                          aol::MultiVector<RealType> &ApproximateCenters,
                          const PictureType &Data ) {
    setCtrlCHandler ( );
    aol::MultiVector<RealType> *approximateDumbbellCenters = new aol::MultiVector<RealType> ( 0, 0 );
    aol::MultiVector<int> approximateDimensions, approximateDumbbellDimensions;
    aol::Vector<RealType> dumbbellSeparations, dumbbellOrientations;
    if ( _progressBar != NULL ) _progressBar->setText ( "AtomFinder: segmenting and finding approximate centers (step 1/3)" );
    getApproximateAtomPositions ( ApproximateCenters, approximateDimensions, *approximateDumbbellCenters, approximateDumbbellDimensions, dumbbellSeparations, dumbbellOrientations, Data );
    aol::MultiVector<RealType> *dumbbellCenters = new aol::MultiVector<RealType> ( 0, 0 ), *dumbbellGaussianParams = new aol::MultiVector<RealType> ( 0, 0 );
    if ( ApproximateCenters.numComponents ( ) > 0 ) {
      if ( _progressBar != NULL ) _progressBar->setText ( "AtomFinder: Gaussian fitting single atoms (step 2/3)" );
      getRefinedAtomPositions ( Centers, GaussianParams, ApproximateCenters, approximateDimensions, Data );
    }
    if ( approximateDumbbellCenters->numComponents ( ) > 0 ) {
      if ( _progressBar != NULL ) _progressBar->setText ( "AtomFinder: Gaussian fitting dumbbells (step 3/3)" );
      getRefinedDumbbellAtomPositions ( *dumbbellCenters, *dumbbellGaussianParams, *approximateDumbbellCenters, approximateDumbbellDimensions, dumbbellSeparations, dumbbellOrientations, Data );
    }
    ApproximateCenters.appendReference ( *approximateDumbbellCenters, true );
    Centers.appendReference ( *dumbbellCenters, true );
    GaussianParams.appendReference ( *dumbbellGaussianParams, true );
  
    unsetCtrlCHandler ( );
  }
  
  void getAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                          const PictureType &Data ) {
    aol::MultiVector<RealType> approximateCenters;
    getAtomPositions ( Centers, GaussianParams, approximateCenters, Data );
  }
  
  void getSingleAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                                aol::MultiVector<RealType> &ApproximateCenters,
                                const PictureType &Data ) {
    setCtrlCHandler ( );
    aol::MultiVector<int> approximateDimensions;
    if ( _progressBar != NULL ) _progressBar->setText ( "AtomFinder: segmenting and finding approximate centers (step 1/2)" );
    getApproximateSingleAtomPositions ( ApproximateCenters, approximateDimensions, Data );
    if ( _progressBar != NULL ) _progressBar->setText ( "AtomFinder: Gaussian fitting single atoms (step 2/2)" );
    getRefinedAtomPositions ( Centers, GaussianParams, ApproximateCenters, approximateDimensions, Data );
    
    unsetCtrlCHandler ( );
  }
  
  void getSingleAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                                const PictureType &Data ) {
    aol::MultiVector<RealType> approximateCenters;
    getSingleAtomPositions ( Centers, GaussianParams, approximateCenters, Data );
  }
  
  void getDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                                  aol::MultiVector<RealType> &ApproximateCenters,
                                  const PictureType &Data ) {
    setCtrlCHandler ( );
    aol::MultiVector<int> approximateDimensions;
    aol::Vector<RealType> dumbbellSeparations, dumbbellOrientations;
    if ( _progressBar != NULL ) _progressBar->setText ( "AtomFinder: segmenting and finding approximate centers (step 1/2)" );
    getApproximateDumbbellAtomPositions ( ApproximateCenters, approximateDimensions, dumbbellSeparations, dumbbellOrientations, Data );
    if ( _progressBar != NULL ) _progressBar->setText ( "AtomFinder: Gaussian fitting dumbbells (step 2/2)" );
    getRefinedDumbbellAtomPositions ( Centers, GaussianParams, ApproximateCenters, approximateDimensions, dumbbellSeparations, dumbbellOrientations, Data );
    
    unsetCtrlCHandler ( );
  }
  
  void getDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                                 const PictureType &Data ) {
    aol::MultiVector<RealType> approximateCenters;
    getDumbbellAtomPositions ( Centers, GaussianParams, approximateCenters, Data );
  }
  
  void getAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                          aol::MultiVector<RealType> &ApproximateCenters,
                          const PictureType &Data, const aol::ParameterParser &Parser ) {
    if ( Parser.hasVariable ( "minCenterValue" ) ) {
      RealType smoothSigma, minCenterValue;
      int eraseInfinityRadius, maxNumCenters;
      readParameters ( Parser, smoothSigma, eraseInfinityRadius, minCenterValue, maxNumCenters );
      aol::Vec2<short> atomSize;
      RealType atomSeparation, atomAngle;
      bool singleFit, dumbbellFit;
      readParameters ( Parser, atomSize, atomSeparation, atomAngle, singleFit, dumbbellFit );
      if ( Parser.checkAndGetBool ( "initializeAtomCentersOnGroundTruth" ) ) {
        qc::ScalarArray<RealType, qc::QC_2D> u ( Parser.getString ( "groundTruthPath" ).c_str ( ) );
        getApproximateAtomPositions ( ApproximateCenters, u, smoothSigma, eraseInfinityRadius, minCenterValue, maxNumCenters );
      } else
        getApproximateAtomPositions ( ApproximateCenters, Data, smoothSigma, eraseInfinityRadius, minCenterValue, maxNumCenters );
      
      if ( dumbbellFit ) {
        const RealType periodDelta = Parser.getDoubleOrDefault ( "periodDelta", 5 );
        const RealType angleDelta = Parser.getDoubleOrDefault ( "dumbbellAngleDelta", Parser.getDoubleOrDefault ( "angleDelta", 5 ) );
        getApproximateDumbbellAtomPositions ( ApproximateCenters, atomSeparation, atomAngle, periodDelta, angleDelta, Data );
        getRefinedDumbbellAtomPositions ( Centers, GaussianParams, ApproximateCenters, atomSize, atomSeparation, atomAngle, Data );
      } else
        getRefinedAtomPositions ( Centers, GaussianParams, ApproximateCenters, atomSize, Data );
    } else {
      aol::Vec2<short> atomSize;
      RealType atomSeparation, atomAngle;
      bool singleFit, dumbbellFit;
      readParameters ( Parser, atomSize, atomSeparation, atomAngle, singleFit, dumbbellFit );
      if ( singleFit ) getAtomPositions ( Centers, GaussianParams, ApproximateCenters, Data, atomSize );
      else if ( dumbbellFit ) getDumbbellAtomPositions ( Centers, GaussianParams, ApproximateCenters, Data, atomSize, atomSeparation, atomAngle );
      else getAtomPositions ( Centers, GaussianParams, ApproximateCenters, Data );
    }
  }
  
  void getAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                          const PictureType &Data, const aol::ParameterParser &Parser ) {
    aol::MultiVector<RealType> approximateCenters;
    getAtomPositions ( Centers, GaussianParams, approximateCenters, Data, Parser );
  }
  
  void analyzeAtoms ( aol::MultiVector<RealType> &Centers, aol::Vector<RealType> &Intensities,
                      aol::MultiVector<RealType> &CentersInitial, aol::Vector<RealType> &IntensitiesInitial,
                      const PictureType &Data, const aol::ParameterParser &Parser ) {
    aol::MultiVector<RealType> gaussianParams;
    getAtomPositions ( Centers, gaussianParams, CentersInitial, Data, Parser );
    getAtomIntensities ( CentersInitial, Data, IntensitiesInitial );
    getAtomIntensities ( gaussianParams, Intensities );
  }
  
  void getAtomPositions ( aol::MultiVector<RealType> &Centers, const PictureType &Data, const aol::ParameterParser &Parser ) {
    aol::MultiVector<RealType> gaussianParams;
    getAtomPositions ( Centers, gaussianParams, Data, Parser );
  }
  
  void getAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                          aol::MultiVector<RealType> &ApproximateCenters,
                          const PictureType &Data, const aol::Vec2<short> &AtomSize ) {
    getApproximateAtomPositions ( ApproximateCenters, Data );
    getRefinedAtomPositions ( Centers, GaussianParams, ApproximateCenters, AtomSize, Data );
  }
  
  void getAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                          const PictureType &Data, const aol::Vec2<short> &AtomSize ) {
    aol::MultiVector<RealType> approximateCenters;
    getAtomPositions ( Centers, GaussianParams, approximateCenters, Data, AtomSize );
  }
  
  void getAtomPositions ( aol::MultiVector<RealType> &Centers,
                         const PictureType &Data, const aol::Vec2<short> &AtomSize ) {
    aol::MultiVector<RealType> gaussianParams;
    getAtomPositions ( Centers, gaussianParams, Data, AtomSize );
  }
  
  void getDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                                  aol::MultiVector<RealType> &ApproximateCenters,
                                  const PictureType &Data, const aol::Vec2<short> &AtomSize,
                                  const RealType AtomSeparation, const RealType AtomAngle ) {
    getApproximateAtomPositions ( ApproximateCenters, Data );
    getRefinedDumbbellAtomPositions ( Centers, GaussianParams, ApproximateCenters, AtomSize, AtomSeparation, AtomAngle, Data );
  }
  
  void getDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                                  const PictureType &Data, const aol::Vec2<short> &AtomSize,
                                  const RealType AtomSeparation, const RealType AtomAngle ) {
    aol::MultiVector<RealType> approximateCenters;
    getDumbbellAtomPositions ( Centers, GaussianParams, approximateCenters, Data, AtomSize, AtomSeparation, AtomAngle );
  }
  
  void getDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers,
                                 const PictureType &Data, const aol::Vec2<short> &AtomSize,
                                 const RealType AtomSeparation, const RealType AtomAngle ) {
    aol::MultiVector<RealType> gaussianParams;
    getDumbbellAtomPositions ( Centers, gaussianParams, Data, AtomSize, AtomSeparation, AtomAngle );
  }
  
  void getApproximateAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<int> &Dimensions,
                                     aol::MultiVector<RealType> &DumbbellCenters, aol::MultiVector<int> &DumbbellDimensions,
                                     aol::Vector<RealType> &DumbbellSeparations, aol::Vector<RealType> &DumbbellOrientations,
                                     const PictureType &Data );
  
  void getApproximateAtomPositions ( aol::MultiVector<RealType> &Centers, const PictureType &Data,
                                     const RealType SmoothSigma, const int EraseInfinityRadius, const RealType MinCenterValue = 0.5, const int MaxNumCenters = 0 );
  
  void getApproximateSingleAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<int> &Dimensions,
                                           const PictureType &Data );
  
  void getApproximateDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<int> &Dimensions,
                                             aol::Vector<RealType> &DumbbellSeparations, aol::Vector<RealType> &DumbbellOrientations,
                                             const PictureType &Data );
  
  void getApproximateDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers,
                                             const RealType AtomSeparation, const RealType AtomAngle, const RealType PeriodDelta, const RealType AngleDelta,
                                             const PictureType &Data );
  
  void getApproximateAtomPositions ( aol::MultiVector<RealType> &Centers, const PictureType &Data );
  
  void getRefinedAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                                 const aol::MultiVector<RealType> &ApproximateCenters, const aol::MultiVector<int> &ApproximateDimensions,
                                 const PictureType &Data, bool ClearBoundaryAtoms = true );
  
  void getRefinedAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                                 const aol::MultiVector<RealType> &ApproximateCenters, const aol::Vec2<short> &AtomSize,
                                 const PictureType &Data, bool ClearBoundaryAtoms = true ) {
    aol::MultiVector<int> approximateDimensions ( ApproximateCenters.numComponents ( ), 2 );
    for ( int k=0; k<approximateDimensions.numComponents ( ) ; ++k ) {
      approximateDimensions[k][0] = AtomSize[0];
      approximateDimensions[k][1] = AtomSize[1];
    }
    getRefinedAtomPositions ( Centers, GaussianParams, ApproximateCenters, approximateDimensions, Data, ClearBoundaryAtoms );
  }
  
  void getRefinedDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                                         const aol::MultiVector<RealType> &ApproximateCenters, const aol::MultiVector<int> &ApproximateDimensions,
                                         const aol::Vector<RealType> &Separations, const aol::Vector<RealType> &Orientations,
                                         const PictureType &Data, bool ClearBoundaryAtoms = true );
  
  void getRefinedDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                                         const aol::MultiVector<RealType> &ApproximateCenters, const aol::Vec2<short> &AtomSize,
                                         const RealType AtomSeparation, const RealType AtomAngle,
                                         const PictureType &Data, bool ClearBoundaryAtoms = true ) {
    aol::MultiVector<int> approximateDimensions ( ApproximateCenters.numComponents ( ), 2 );
    for ( int k=0; k<approximateDimensions.numComponents ( ) ; ++k ) {
      approximateDimensions[k][0] = AtomSize[0];
      approximateDimensions[k][1] = AtomSize[1];
    }
    aol::Vector<RealType> separations ( ApproximateCenters.numComponents ( ) ), orientations ( ApproximateCenters.numComponents ( ) );
    orientations.setAll ( AtomAngle * aol::NumberTrait<RealType>::pi / 180.0 );
    separations.setAll ( AtomSeparation );
    getRefinedDumbbellAtomPositions ( Centers, GaussianParams, ApproximateCenters, approximateDimensions, separations, orientations, Data, ClearBoundaryAtoms );
  }

  void readCentersFromCSV ( aol::MultiVector<RealType> &Centers, const std::string &Path ) const;
  
  void setGamma ( const RealType Gamma ) {
    _gamma = Gamma;
  }
  
  void setEpsilon ( const RealType Epsilon ) {
    _epsilon = Epsilon;
  }
  
  void setMaxIt ( const int MaxIt ) {
    _maxIt = MaxIt;
  }
  
  void setNumOuterIterations ( const int NumOuterIterations ) {
    _numOuterIterations = NumOuterIterations;
  }
  
  const qc::BitArray<qc::QC_2D>& getSegmentedConstRef ( ) {
    return _segmented;
  }
  
  void setOutputDir ( const std::string &OutputDir ) {
    _outputDir = OutputDir;
    _diskOutput = OutputDir != "";
  }
  
  void setQuietMode ( const bool Quiet = true ) {
    _quietMode = Quiet;
  }
  
  static void getAtomIntensities ( const aol::MultiVector<RealType> &Centers, const PictureType &Data, aol::Vector<RealType> &Intensities ) {
    Intensities.reallocate ( 0 );
    for ( int k=0; k<Centers.numComponents ( ) ; ++k )
      Intensities.pushBack ( Data.interpolate ( Centers[k][0], Centers[k][1] ) );
  }
  
  static void getAtomIntensities ( const aol::MultiVector<RealType> &GaussianParams, aol::Vector<RealType> &Intensities ) {
    Intensities.reallocate ( 0 );
    for ( int k=0; k<GaussianParams.numComponents ( ); ++k ) {
      if ( GaussianParams[k].size ( ) == AsymmetricGaussianBumpFunction<RealType>::NumberOfParameters )
        Intensities.pushBack ( GaussianParams[k][2] + GaussianParams[k][6] );
      else if ( GaussianParams[k].size ( ) == AsymmetricGaussianDoubleBumpFunction<RealType>::NumberOfParameters ) {
        Intensities.pushBack ( GaussianParams[k][4] + GaussianParams[k][12] );
        Intensities.pushBack ( GaussianParams[k][5] + GaussianParams[k][12] );
      }
    }
  }
  
  static void getBumpFunctionImage ( const aol::MultiVector<RealType> &GaussianParameters, PictureType &BumpFunctionImg ) {
    BumpFunctionImg.setZero ( );
    for ( int k=0; k<GaussianParameters.numComponents ( ) ; ++k ) {
      if ( GaussianParameters[k].size ( ) == AsymmetricGaussianBumpFunction<RealType>::NumberOfParameters ) {
        const AsymmetricGaussianBumpFunction<RealType> bumpFunc ( GaussianParameters[k] );
        for ( int x=0; x<BumpFunctionImg.getNumX ( ) ; ++x )
          for ( int y=0; y<BumpFunctionImg.getNumY ( ) ; ++y )
            BumpFunctionImg.add ( x, y, bumpFunc.evaluate ( aol::Vec2<RealType> ( x, y ) ) - GaussianParameters[k][6] );
        BumpFunctionImg.addToAll ( 1.0 / static_cast<RealType> ( GaussianParameters.numComponents ( ) ) * GaussianParameters[k][6] );
      } else if ( GaussianParameters[k].size ( ) == AsymmetricGaussianDoubleBumpFunction<RealType>::NumberOfParameters ) {
        const AsymmetricGaussianDoubleBumpFunction<RealType> bumpFunc ( GaussianParameters[k] );
        for ( int x=0; x<BumpFunctionImg.getNumX ( ) ; ++x )
          for ( int y=0; y<BumpFunctionImg.getNumY ( ) ; ++y )
            BumpFunctionImg.add ( x, y, bumpFunc.evaluate ( aol::Vec2<RealType> ( x, y ) ) - GaussianParameters[k][12] );
        BumpFunctionImg.addToAll ( 1.0 / static_cast<RealType> ( GaussianParameters.numComponents ( ) ) * GaussianParameters[k][12] );
      }
    }
  }
  
  static void getAtomCentersImage ( const PictureType &Data, const aol::MultiVector<RealType> &Centers,
                                    ColoredPictureType &AtomCentersImg ) {
    AtomCentersImg.reallocate ( Data );
    for ( short c=0; c<3 ; ++c ) AtomCentersImg[c] = Data;
    AtomCentersImg.scaleValuesTo01 ( );
    aol::Vec2<short> pos;
    for ( int i=0; i<Centers.numComponents ( ) ; ++i ) {
      pos.set ( aol::Rint ( Centers[i][0] ), aol::Rint ( Centers[i][1] ) );
      AtomCentersImg[0].set ( pos, 1 );
      AtomCentersImg[1].set ( pos, 0 );
      AtomCentersImg[2].set ( pos, 0 );
    }
    AtomCentersImg.setOverflowHandlingToCurrentValueRange ( );
  }
  
  void addDumbbellOrientationLinesToAtomCentersImage ( const int kDumbbell, aol::MultiVector<RealType> &DumbbellCenters,
                                                       aol::Vector<RealType> &DumbbellSeparations, aol::Vector<RealType> &DumbbellOrientations ) const {
    const std::string geometricCentersPath = aol::strprintf ( "%s/centers_initial.png", _outputDir.c_str ( ) );
    ColoredPictureType geometricCentersImg ( geometricCentersPath.c_str ( ) );
    aol::Vec2<short> pos;
    RealType maxVal = geometricCentersImg.getMaxValue ( );
    for ( int k=0; k<kDumbbell ; ++k ) {
      for ( short i=-DumbbellSeparations[k]; i<=DumbbellSeparations[k] ; ++i ) {
        pos.set ( aol::Rint ( DumbbellCenters[k][0] + cos ( DumbbellOrientations[k] ) * i ), aol::Rint ( DumbbellCenters[k][1] + sin ( DumbbellOrientations[k] ) * i ) );
        if ( pos[0] >=0 && pos[0] < geometricCentersImg[0].getNumX ( ) && pos[1] >=0 && pos[1] < geometricCentersImg[0].getNumY ( ) ) {
          geometricCentersImg[0].set ( pos, 0 );
          geometricCentersImg[1].set ( pos, maxVal );
          geometricCentersImg[2].set ( pos, 0 );
        }
      }
      pos.set ( DumbbellCenters[k][0], DumbbellCenters[k][1] );
      geometricCentersImg[0].set ( pos, maxVal );
      geometricCentersImg[1].set ( pos, 0 );
    }
    geometricCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    geometricCentersImg.savePNG ( geometricCentersPath.c_str ( ) );
  }
  
  void setProgressBar ( aol::ProgressBar<> *ProgressBar ) {
    _progressBar = ProgressBar;
  }
  
  void setCatchCtrlC ( bool catchCtrlC ) {
    _catchCtrlC = catchCtrlC;
  }
protected:
  void getComponentsCollection ( qc::ComponentsCollection<RealType> &ComponentsCollection, const PictureType &Data );
  
  void readParameters ( const aol::ParameterParser &Parser,
                        aol::Vec2<short> &AtomSize, RealType &AtomSeparation, RealType &AtomAngle,
                        bool &SingleFit, bool &DumbbellFit ) {
    readParameters ( Parser );
    _gamma = Parser.getDoubleOrDefault ( "gamma", gammaDefault );
    _maxIt = Parser.getIntOrDefault ( "maxIt", maxItDefault );
    _epsilon = Parser.getDoubleOrDefault ( "epsilon", epsilonDefault );
    _numOuterIterations = Parser.getDoubleOrDefault ( "numOuterIterations", numOuterIterationsDefault );
    AtomSize.set ( Parser.getIntOrDefault ( "atomWidth", 0 ), Parser.getIntOrDefault ( "atomHeight", 0 ) );
    AtomSeparation = Parser.getDoubleOrDefault ( "atomSeparation", 0 );
    AtomAngle = Parser.getDoubleOrDefault ( "atomAngle", 0 );
    SingleFit = ( AtomSize[0] > 0 && AtomSize[1] > 0 && AtomSize[0] % 2 != 0 && AtomSize[1] % 2 != 0 );
    DumbbellFit = ( Parser.checkAndGetBool ( "dumbbellFit" ) && SingleFit && AtomSeparation > 0 && Parser.hasVariable ( "atomAngle" ) );
    if ( DumbbellFit ) SingleFit = false;
  }
  
  void readParameters ( const aol::ParameterParser &Parser, RealType &SmoothSigma, int &EraseInfinityRadius, RealType &MinCenterValue, int &MaxNumCenters ) {
    readParameters ( Parser );
    SmoothSigma = Parser.getDoubleOrDefault ( "smoothSigma", 0 );
    EraseInfinityRadius = Parser.getDoubleOrDefault ( "eraseInfinityRadius", 11 );
    MinCenterValue = Parser.getDouble ( "minCenterValue" );
    MaxNumCenters = Parser.getIntOrDefault ( "maxNumCenters", 0 );
  }
  
  void readParameters ( const aol::ParameterParser &Parser ) {
    _lmVerbose = Parser.checkAndGetBool ( "lmVerbose" );
  }
  
  void setCtrlCHandler () const {
    if (_catchCtrlC)
      _previousCtrlCHandler = signal ( InterruptSignal, aol::ctrlCHandler );
  }
  
  void unsetCtrlCHandler () const {
    if (_catchCtrlC)
      signal ( InterruptSignal, _previousCtrlCHandler );
  }
  
  bool wantsInterrupt() const {
    if (!_catchCtrlC || !aol::getCtrlCState())
      return false;
    else
      return true;
  }
};

template <typename RealType, typename MatrixType, typename LinearRegressionType, typename ScalarPictureType, typename ColoredPictureType>
RealType AtomFinder<RealType, MatrixType, LinearRegressionType, ScalarPictureType, ColoredPictureType>::gammaDefault = 1e-4;

template <typename RealType, typename MatrixType, typename LinearRegressionType, typename ScalarPictureType, typename ColoredPictureType>
RealType AtomFinder<RealType, MatrixType, LinearRegressionType, ScalarPictureType, ColoredPictureType>::epsilonDefault = 1e-3;

template <typename RealType, typename MatrixType, typename LinearRegressionType, typename ScalarPictureType, typename ColoredPictureType>
int AtomFinder<RealType, MatrixType, LinearRegressionType, ScalarPictureType, ColoredPictureType>::maxItDefault = 10000;

template <typename RealType, typename MatrixType, typename LinearRegressionType, typename ScalarPictureType, typename ColoredPictureType>
int AtomFinder<RealType, MatrixType, LinearRegressionType, ScalarPictureType, ColoredPictureType>::numOuterIterationsDefault = 1;
  
} // end namespace


#endif
