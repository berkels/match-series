#include "atomFinder.h"

namespace im {

template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
                _ColoredPictureType>::getComponentsCollection ( qc::ComponentsCollection<RealType> &ComponentsCollection, const PictureType &Data ) {
  if ( !_quietMode )
    std::cerr << "Approximating atom positions as geometric centers of connected components.." << std::endl;
  
  PictureType u ( Data );
  u /= u.getMaxValue ( );
  u.setOverflowHandling( aol::CLIP_THEN_SCALE, 0, 1 );
  
  InitType grid ( qc::GridSize<ConfiguratorType::Dim> ( u.getSize ( ) ) );
  
  ArrayType uArr ( u, aol::FLAT_COPY );
  im::PiecewiseConstantTwoPhaseMSSegmentor<ConfiguratorType, 1, FirstOrderPrimalDualTwoPhaseMSSegmentor<ConfiguratorType> > segmentor ( grid, _gamma, uArr, true, true );
  segmentor.setQuietMode ( this->_quietMode );
  segmentor.setCatchCtrlC ( true );
  PictureType segmentation ( grid );
  segmentor.setMaxIterations ( _maxIt );
  segmentor.setStopEpsilon ( _epsilon );
  segmentor.setOuterIterations ( _numOuterIterations );
  segmentor.segmentAndAdjustGrayValues ( segmentation );
  
  if ( _diskOutput ) {
    segmentation.setOverflowHandlingToCurrentValueRange ( );
    segmentation.savePNG ( aol::strprintf ( "%s/seg.png", _outputDir.c_str ( ) ).c_str ( ) );
  }

  segmentation.threshold ( 0.5, 0, 1 );
  
  if ( _diskOutput ) {
    segmentation.setOverflowHandlingToCurrentValueRange ( );
    segmentation.savePNG ( aol::strprintf ( "%s/hardSeg.png", _outputDir.c_str ( ) ).c_str ( ) );
  }
  
  qc::BitArray<qc::QC_2D> mask ( grid );
  for ( int i=0; i<segmentation.getNumX ( ); ++i )
    for ( int j=0; j<segmentation.getNumY ( ); ++j )
      mask.set ( i, j, ( segmentation.get ( i, j ) == 1 ) );

  // We assume that atoms intensities are higher than the void intensities.
  // So we have to make sure that the mask corresponds to regions of higher intensities.
  if ( segmentor.getMeanValuesReference()[0][0] > segmentor.getMeanValuesReference()[1][0] )
    mask.invert();

  _segmented.reallocate ( mask.getNumX ( ), mask.getNumY ( ) );
  _segmented = mask;
  
  ComponentsCollection.initializeFrom ( mask );
  
  if ( !_quietMode )
    std::cerr << "# of atoms after segmentation: " << ComponentsCollection.getNumNonEmptyComponents ( ) << std::endl;
  
  if ( _diskOutput ) {
    PictureType connectedRegions ( grid );
    ComponentsCollection.createPictureBWComponents ( connectedRegions );
    connectedRegions.save ( aol::strprintf ( "%s/connectedRegions%s", _outputDir.c_str ( ), getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    
    ColoredPictureType geometricCentersImg ( grid );
    ComponentsCollection.createPictureBWComponentsRedGeometricCenters ( geometricCentersImg );
    geometricCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    geometricCentersImg.savePNG ( aol::strprintf ( "%s/centers_initial.png", _outputDir.c_str ( ) ).c_str ( ) );
  }
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
_ColoredPictureType>::getApproximateAtomPositions ( aol::MultiVector<RealType> &Centers, const PictureType &Data,
                                                    const RealType SmoothSigma, const int EraseInfinityRadius, const RealType MinCenterValue, const int MaxNumCenters ) {
  typedef qc::RectangularGridConfigurator<RealType, qc::QC_2D, aol::GaussQuadrature<RealType,qc::QC_2D,3> > ConfiguratorType;
  BumpFitCenterEstimator<ConfiguratorType> centerEstimator ( SmoothSigma, EraseInfinityRadius, MinCenterValue, MaxNumCenters );
  centerEstimator.guessCenters ( Data, Centers );
  
  if ( _diskOutput ) {
    ColoredPictureType atomCentersImg ( Data.getNumX ( ), Data.getNumY ( ) );
    getAtomCentersImage ( Data, Centers, atomCentersImg );
    atomCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    atomCentersImg.savePNG ( aol::strprintf ( "%s/centers_initial.png", _outputDir.c_str ( ) ).c_str ( ) );
  }
}
  

template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
_ColoredPictureType>::getApproximateAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<int> &Dimensions,
                                                    aol::MultiVector<RealType> &DumbbellCenters, aol::MultiVector<int> &DumbbellDimensions,
                                                    aol::Vector<RealType> &DumbbellSeparations, aol::Vector<RealType> &DumbbellOrientations,
                                                    const PictureType &Data ) {
  qc::ComponentsCollection<RealType> componentsCollection;
  getComponentsCollection ( componentsCollection, Data );
  
  // Split components into dumbbells and single atoms and collect orientation and separation for dumbbells
  int numNonEmptyComps = componentsCollection.getNumNonEmptyComponents ( );
  Centers.reallocate ( numNonEmptyComps, 2 );
  Dimensions.reallocate ( numNonEmptyComps, 2 );
  DumbbellCenters.reallocate ( numNonEmptyComps, 2 );
  DumbbellDimensions.reallocate ( numNonEmptyComps, 2 );
  DumbbellOrientations.reallocate ( numNonEmptyComps );
  DumbbellSeparations.reallocate ( numNonEmptyComps );
  int kSingle = 0, kDumbbell = 0;
  for ( NonEmptyComponentsIterator it ( componentsCollection ); it.notAtEnd ( ) ; ++it ) {
    const qc::RandomAccessStencil<short> boundaryIndices = (*it).getBoundaryIndices ( );
    aol::Vec2<short> geometricCenter = (*it).getGeometricCenter ( );
    const short width = 1 + 2 * aol::Rint ( 1.25 * aol::Max<short> ( aol::Abs<short> ( boundaryIndices.W ( ) - geometricCenter[0] ), aol::Abs<short> ( boundaryIndices.E ( ) - geometricCenter[0] ) ) );
    const short height = 1 + 2 * aol::Rint ( 1.25 * aol::Max<short> ( aol::Abs<short> ( boundaryIndices.N ( ) - geometricCenter[1] ), aol::Abs<short> (  boundaryIndices.S ( ) - geometricCenter[1] ) ) );
    aol::Vec<4,RealType> minMaxDiamAngle = (*it).getMinMaxDiameterAndAngle ( );
    if ( minMaxDiamAngle[0] > 1.5 * ( minMaxDiamAngle[2] ) ) {
      DumbbellCenters[kDumbbell][0] = geometricCenter[0];
      DumbbellCenters[kDumbbell][1] = geometricCenter[1];
      DumbbellDimensions[kDumbbell][0] = width;
      DumbbellDimensions[kDumbbell][1] = height;
      DumbbellOrientations[kDumbbell] = minMaxDiamAngle[1];
      DumbbellSeparations[kDumbbell] = 0.5 * minMaxDiamAngle[0];
      ++kDumbbell;
    } else {
      Centers[kSingle][0] = geometricCenter[0];
      Centers[kSingle][1] = geometricCenter[1];
      Dimensions[kSingle][0] = width;
      Dimensions[kSingle][1] = height;
      ++kSingle;
    }
  }
  DumbbellCenters.resize ( kDumbbell, 2 );
  DumbbellDimensions.resize ( kDumbbell, 2 );
  DumbbellOrientations.resize ( kDumbbell );
  DumbbellSeparations.resize ( kDumbbell );
  Centers.resize ( kSingle, 2 );
  Dimensions.resize ( kSingle, 2 );
  
  if ( _diskOutput )
    addDumbbellOrientationLinesToAtomCentersImage ( kDumbbell, DumbbellCenters, DumbbellSeparations, DumbbellOrientations );
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
                _ColoredPictureType>::getApproximateSingleAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<int> &Dimensions,
                                                                          const PictureType &Data ) {
  qc::ComponentsCollection<RealType> componentsCollection;
  getComponentsCollection ( componentsCollection, Data );
  
  // Split components into dumbbells and single atoms and collect orientation and separation for dumbbells
  int numNonEmptyComps = componentsCollection.getNumNonEmptyComponents ( );
  Centers.reallocate ( numNonEmptyComps, 2 );
  Dimensions.reallocate ( numNonEmptyComps, 2 );
  int k = 0;
  for ( NonEmptyComponentsIterator it ( componentsCollection ); it.notAtEnd ( ) ; ++it ) {
    const qc::RandomAccessStencil<short> boundaryIndices = (*it).getBoundaryIndices ( );
    aol::Vec2<short> geometricCenter = (*it).getGeometricCenter ( );
    const short width = 1 + 2 * aol::Rint ( 1.25 * aol::Max<short> ( aol::Abs<short> ( boundaryIndices.W ( ) - geometricCenter[0] ), aol::Abs<short> ( boundaryIndices.E ( ) - geometricCenter[0] ) ) );
    const short height = 1 + 2 * aol::Rint ( 1.25 * aol::Max<short> ( aol::Abs<short> ( boundaryIndices.N ( ) - geometricCenter[1] ), aol::Abs<short> (  boundaryIndices.S ( ) - geometricCenter[1] ) ) );
    Centers[k][0] = geometricCenter[0];
    Centers[k][1] = geometricCenter[1];
    Dimensions[k][0] = width;
    Dimensions[k][1] = height;
    ++k;
  }
  Centers.resize ( k, 2 );
  Dimensions.resize ( k, 2 );
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
_ColoredPictureType>::getApproximateDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<int> &Dimensions,
                                                            aol::Vector<RealType> &DumbbellSeparations, aol::Vector<RealType> &DumbbellOrientations,
                                                            const PictureType &Data ) {
  qc::ComponentsCollection<RealType> componentsCollection;
  getComponentsCollection ( componentsCollection, Data );
  
  // Split components into dumbbells and single atoms and collect orientation and separation for dumbbells
  int numNonEmptyComps = componentsCollection.getNumNonEmptyComponents ( );
  Centers.reallocate ( numNonEmptyComps, 2 );
  Dimensions.reallocate ( numNonEmptyComps, 2 );
  DumbbellOrientations.reallocate ( numNonEmptyComps );
  DumbbellSeparations.reallocate ( numNonEmptyComps );
  int k = 0;
  for ( NonEmptyComponentsIterator it ( componentsCollection ); it.notAtEnd ( ) ; ++it ) {
    const qc::RandomAccessStencil<short> boundaryIndices = (*it).getBoundaryIndices ( );
    aol::Vec2<short> geometricCenter = (*it).getGeometricCenter ( );
    const short width = 1 + 2 * aol::Rint ( 1.25 * aol::Max<short> ( aol::Abs<short> ( boundaryIndices.W ( ) - geometricCenter[0] ), aol::Abs<short> ( boundaryIndices.E ( ) - geometricCenter[0] ) ) );
    const short height = 1 + 2 * aol::Rint ( 1.25 * aol::Max<short> ( aol::Abs<short> ( boundaryIndices.N ( ) - geometricCenter[1] ), aol::Abs<short> (  boundaryIndices.S ( ) - geometricCenter[1] ) ) );
    aol::Vec<4,RealType> minMaxDiamAngle = (*it).getMinMaxDiameterAndAngle ( );
    Centers[k][0] = geometricCenter[0];
    Centers[k][1] = geometricCenter[1];
    Dimensions[k][0] = width;
    Dimensions[k][1] = height;
    DumbbellOrientations[k] = minMaxDiamAngle[1];
    DumbbellSeparations[k] = 0.5 * minMaxDiamAngle[0];
    ++k;
  }
  Centers.resize ( k, 2 );
  Dimensions.resize ( k, 2 );
  DumbbellOrientations.resize ( k );
  DumbbellSeparations.resize ( k );
  
  if ( _diskOutput )
    addDumbbellOrientationLinesToAtomCentersImage ( k, Centers, DumbbellSeparations, DumbbellOrientations );
}
  
  
template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
_ColoredPictureType>::getApproximateDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, const RealType AtomSeparation, const RealType AtomAngle,
                                                            const RealType PeriodDelta, const RealType AngleDelta,
                                                            const PictureType &Data ) {
  aol::MultiVector<RealType> centers ( Centers );
  Centers.clear ( );
  
  while ( centers.numComponents ( ) > 1 ) {
    aol::Vector<RealType> dists ( centers.numComponents ( ), aol::NumberTrait<RealType>::Inf );
    for ( int i=1 ; i<centers.numComponents ( ) ; ++i ) {
      const RealType dx = centers[0][0] - centers[i][0];
      const RealType dy = centers[0][1] - centers[i][1];
      const RealType angle = atan2 ( dy, dx ) * 180 / aol::NumberTrait<RealType>::pi;
      const RealType dist = aol::Vec2<RealType> ( dx, dy ).norm ( );
      if ( aol::Abs<RealType> ( dist - AtomSeparation ) < PeriodDelta && aol::Min<RealType> ( aol::Abs<RealType> ( angle - AtomAngle ), aol::Abs<RealType> ( angle + AtomAngle ) ) < AngleDelta )
        dists[i] = dist;
    }
    const int iOpt = dists.getMinIndexAndValue ( ).first;
    if ( iOpt > 0 ) {
      Centers.resize ( Centers.numComponents ( ) + 1, 2 );
      Centers[Centers.numComponents ( ) - 1][0] = 0.5 * ( centers[0][0] + centers[iOpt][0] );
      Centers[Centers.numComponents ( ) - 1][1] = 0.5 * ( centers[0][1] + centers[iOpt][1] );
      centers.eraseComponent ( iOpt );
    }
    centers.eraseComponent ( 0 );
  }
  
  if ( _diskOutput ) {
    ColoredPictureType atomCentersImg ( Data.getNumX ( ), Data.getNumY ( ) );
    getAtomCentersImage ( Data, Centers, atomCentersImg );
    atomCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    atomCentersImg.savePNG ( aol::strprintf ( "%s/centers_dumbbell_initial.png", _outputDir.c_str ( ) ).c_str ( ) );
  }
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
                _ColoredPictureType>::getApproximateAtomPositions ( aol::MultiVector<RealType> &Centers, const PictureType &Data ) {
  qc::ComponentsCollection<RealType> componentsCollection;
  getComponentsCollection ( componentsCollection, Data );
  
  if ( !_quietMode )
    std::cerr << "# of atoms after segmentation: " << componentsCollection.getNumNonEmptyComponents ( ) << std::endl;
  
  if ( _diskOutput ) {
    PictureType connectedRegions ( Data, aol::STRUCT_COPY );
    componentsCollection.createPictureBWComponents ( connectedRegions );
    connectedRegions.save ( aol::strprintf ( "%s/connectedRegions%s", _outputDir.c_str ( ), getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    
    ColoredPictureType geometricCentersImg ( Data.getNumX ( ), Data.getNumY ( ) );
    componentsCollection.createPictureBWComponentsRedGeometricCenters ( geometricCentersImg );
    geometricCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    geometricCentersImg.savePNG ( aol::strprintf ( "%s/centers_initial.png", _outputDir.c_str ( ) ).c_str ( ) );
  }
  
  short k = 0;
  Centers.reallocate ( componentsCollection.getNumNonEmptyComponents ( ), 2 );
  for ( NonEmptyComponentsIterator it ( componentsCollection ); it.notAtEnd ( ) ; ++it ) {
    aol::Vec2<short> geometricCenter = (*it).getGeometricCenter ( );
    Centers[k][0] = geometricCenter[0];
    Centers[k][1] = geometricCenter[1];
    ++k;
  }
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
_ColoredPictureType>::getRefinedAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                                                const aol::MultiVector<RealType> &ApproximateCenters, const aol::MultiVector<int> &ApproximateDimensions,
                                                const PictureType &Data, bool ClearBoundaryAtoms ) {
  if ( !_quietMode )
    std::cerr << "Refining single atom positions via 2D Gaussian fit.." << std::endl;
  
  Centers.reallocate ( ApproximateCenters.numComponents ( ), 2 );
  const short numVar = AsymmetricGaussianBumpFunction<RealType>::NumberOfParameters;
  GaussianParams.reallocate ( ApproximateCenters.numComponents ( ), numVar );
  RealType approxAtomRadius = 0;
  if ( _progressBar != NULL ) _progressBar->start ( ApproximateCenters.numComponents ( ) );
  for ( short k=0; k<ApproximateCenters.numComponents ( ) && !wantsInterrupt ( ) ; ++k ) {
    aol::Vec2<short> maxAtomOffset ( ( ApproximateDimensions[k][0] - 1 ) / 2, ( ApproximateDimensions[k][1] - 1 ) / 2 ), atomOffset, atomSize, approxCenterInt;
    approxCenterInt.set ( aol::Rint ( ApproximateCenters[k][0] ), aol::Rint ( ApproximateCenters[k][1] ) );
    atomOffset.set ( aol::Min<short> ( aol::Min<short> ( maxAtomOffset[0], approxCenterInt[0] ), Data.getNumX ( ) - 1 - approxCenterInt[0] ),
                    aol::Min<short> ( aol::Min<short> ( maxAtomOffset[1], approxCenterInt[1] ), Data.getNumY ( ) - 1 - approxCenterInt[1] ) );
    atomSize.set ( 2 * atomOffset[0] + 1, 2 * atomOffset[1] + 1 );
    aol::Vec2<short> corner ( approxCenterInt[0] - atomOffset[0], approxCenterInt[1] - atomOffset[1] );
    PictureType atomData ( atomSize[0], atomSize[1] );
    for ( short dx=0; dx<atomSize[0] ; ++dx )
      for ( short dy=0; dy<atomSize[1] ; ++dy )
        atomData.set ( dx, dy, Data.get ( corner[0] + dx, corner[1] + dy ) );
    SingleAsymmetricBumpFitTargetFunctional<RealType> gaussNewtonF ( atomData );
    SingleAsymmetricBumpFitTargetJacobian<RealType> gaussNewtonDF ( atomData );
    aol::Vector<RealType> lowerBounds ( numVar );
    aol::Vector<RealType> upperBounds ( numVar );
    lowerBounds[0] = 1;
    upperBounds[0] = atomSize[0]-2;
    lowerBounds[1] = 1;
    upperBounds[1] = atomSize[1]-2;
    lowerBounds[2] = atomData.getMaxValue ( ) - atomData.getMeanValue ( );
    upperBounds[2] = atomData.getMaxValue ( ) - atomData.getMinValue ( );
    lowerBounds[3] = 0;
    upperBounds[3] = 0.5 * atomSize[0];
    lowerBounds[4] = 0;
    upperBounds[4] = 0.5 * atomSize[1];
    lowerBounds[6] = atomData.getMinValue ( );
    upperBounds[6] = atomData.getMeanValue ( );
    aol::BitVector constrainedDirections ( numVar );
    constrainedDirections.set ( 0, true );
    constrainedDirections.set ( 1, true );
    constrainedDirections.set ( 2, true );
    constrainedDirections.set ( 3, true );
    constrainedDirections.set ( 4, true );
    constrainedDirections.set ( 6, true );
    ProjectorType boxProjector ( lowerBounds, upperBounds, constrainedDirections );
    aol::ConstrainedLevenbergMarquardtAlgorithm<RealType, MatrixType, LinearRegressionType, ProjectorType> levenbergMarquardtAlg ( atomData.size ( ), gaussNewtonF, gaussNewtonDF, boxProjector, 50, 1, 0.2, 0.8, 1e-6, 1e-6, 1e-6, _lmVerbose );
    aol::Vector<RealType> Arg ( numVar ), Dest ( numVar );
    Arg[0] = atomOffset[0];                                                 // Center x
    Arg[1] = atomOffset[1];                                                 // Center y
    Arg[2] = 0.5 * ( atomData.getMaxValue ( ) - atomData.getMinValue ( ) ); // Intensity
    Arg[3] = 0.25 * atomSize[0];                                            // Width
    Arg[4] = 0.25 * atomSize[1];                                            // Height
    Arg[5] = 0;                                                             // Rotation
    Arg[6] = atomData.getMinValue ( );                                      // Offset
    if ( atomOffset.getMinValue ( ) > 0 )
      levenbergMarquardtAlg.apply ( Arg, Dest );
    else
      Dest = Arg;
    Dest[0] += corner[0];
    Dest[1] += corner[1];
    Centers[k][0] = Dest[0];
    Centers[k][1] = Dest[1];
    
    for ( int j=0; j<numVar ; ++j )
      GaussianParams[k][j] = Dest[j];
    
    approxAtomRadius += Dest[3] + Dest[4];
    
    if ( _progressBar != NULL ) (*_progressBar)++;
  }
  if ( _progressBar != NULL ) _progressBar->finish ( );
  
  approxAtomRadius /= 2 * ApproximateCenters.numComponents ( );
  
  // Clear boundary atoms
  if ( ClearBoundaryAtoms ) {
    short k = 0;
    while ( k < Centers.numComponents ( ) ) {
      if ( Centers[k][0] - approxAtomRadius < 0 || Centers[k][0] + approxAtomRadius >= Data.getNumX ( )
          || Centers[k][1] - approxAtomRadius < 0 || Centers[k][1] + approxAtomRadius >= Data.getNumY ( ) ) {
        Centers.eraseComponent ( k );
        GaussianParams.eraseComponent ( k );
      } else ++k;
    }
  }
  
  if ( !_quietMode )
    std::cerr << "# of (inner) atoms after bump fitting: " << Centers.numComponents ( ) << std::endl;
  
  if ( _diskOutput && GaussianParams.numComponents ( ) > 0 ) {
    PictureType bumpFunctionImg ( Data, aol::STRUCT_COPY );
    getBumpFunctionImage ( GaussianParams, bumpFunctionImg );
    bumpFunctionImg.save ( aol::strprintf ( "%s/bumpFuncs_refined%s", _outputDir.c_str ( ), getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    
    ColoredPictureType atomCentersImg ( Data.getNumX ( ), Data.getNumY ( ) );
    getAtomCentersImage ( Data, Centers, atomCentersImg );
    atomCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    atomCentersImg.savePNG ( aol::strprintf ( "%s/centers_refined.png", _outputDir.c_str ( ) ).c_str ( ) );
  }
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
_ColoredPictureType>::getRefinedDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                                                        const aol::MultiVector<RealType> &ApproximateCenters, const aol::MultiVector<int> &ApproximateDimensions,
                                                        const aol::Vector<RealType> &Separations, const aol::Vector<RealType> &Orientations,
                                                        const PictureType &Data, bool ClearBoundaryAtoms ) {
  if ( !_quietMode )
    std::cerr << "Refining positions within dumbbells via 2D Gaussian fit.." << std::endl;
  
  Centers.reallocate ( 2 * ApproximateCenters.numComponents ( ), 2 );
  const short numVar = AsymmetricGaussianDoubleBumpFunction<RealType>::NumberOfParameters;
  GaussianParams.reallocate ( ApproximateCenters.numComponents ( ), numVar );
  RealType approxAtomRadius = 0;
  if ( _progressBar != NULL ) _progressBar->start ( ApproximateCenters.numComponents ( ) );
  for ( short k=0; k<ApproximateCenters.numComponents ( ) && !wantsInterrupt ( ) ; ++k ) {
    const RealType atomAngleRads = Orientations[k], atomSeparation = Separations[k];
    aol::Vec2<short> maxAtomOffset ( ( ApproximateDimensions[k][0] - 1 ) / 2, ( ApproximateDimensions[k][1] - 1 ) / 2 ), atomOffset, atomSize, approxCenterInt;
    approxCenterInt.set ( aol::Rint ( ApproximateCenters[k][0] ), aol::Rint ( ApproximateCenters[k][1] ) );
    atomOffset.set ( aol::Min<short> ( aol::Min<short> ( maxAtomOffset[0], approxCenterInt[0] ), Data.getNumX ( ) - 1 - approxCenterInt[0] ),
                    aol::Min<short> ( aol::Min<short> ( maxAtomOffset[1], approxCenterInt[1] ), Data.getNumY ( ) - 1 - approxCenterInt[1] ) );
    atomSize.set ( 2 * atomOffset[0] + 1, 2 * atomOffset[1] + 1 );
    aol::Vec2<short> corner ( approxCenterInt[0] - atomOffset[0], approxCenterInt[1] - atomOffset[1] );
    PictureType atomData ( atomSize[0], atomSize[1] );
    for ( short dx=0; dx<atomSize[0] ; ++dx )
      for ( short dy=0; dy<atomSize[1] ; ++dy )
        atomData.set ( dx, dy, Data.get ( corner[0] + dx, corner[1] + dy ) );
    SingleAsymmetricDoubleBumpFitTargetFunctional<RealType> gaussNewtonF ( atomData );
    SingleAsymmetricDoubleBumpFitTargetJacobian<RealType> gaussNewtonDF ( atomData );
    aol::Vector<RealType> lowerBounds ( numVar );
    aol::Vector<RealType> upperBounds ( numVar );
    lowerBounds[0] = 1;
    upperBounds[0] = atomSize[0]-2;
    lowerBounds[1] = 1;
    upperBounds[1] = atomSize[1]-2;
    lowerBounds[2] = 1;
    upperBounds[2] = atomSize[0]-2;
    lowerBounds[3] = 1;
    upperBounds[3] = atomSize[1]-2;
    lowerBounds[4] = atomData.getMaxValue ( ) - atomData.getMeanValue ( );
    upperBounds[4] = atomData.getMaxValue ( ) - atomData.getMinValue ( );
    lowerBounds[5] = atomData.getMaxValue ( ) - atomData.getMeanValue ( );
    upperBounds[5] = atomData.getMaxValue ( ) - atomData.getMinValue ( );
    lowerBounds[6] = 0;
    upperBounds[6] = 0.5 * atomSize.getMinValue ( );
    lowerBounds[7] = 0;
    upperBounds[7] = 0.5 * atomSize.getMinValue ( );
    lowerBounds[8] = 0;
    upperBounds[8] = 0.5 * atomSize.getMinValue ( );
    lowerBounds[9] = 0;
    upperBounds[9] = 0.5 * atomSize.getMinValue ( );
    lowerBounds[12] = atomData.getMinValue ( );
    upperBounds[12] = atomData.getMeanValue ( );
    aol::BitVector constrainedDirections ( numVar );
    constrainedDirections.set ( 0, true );
    constrainedDirections.set ( 1, true );
    constrainedDirections.set ( 2, true );
    constrainedDirections.set ( 3, true );
    constrainedDirections.set ( 4, true );
    constrainedDirections.set ( 5, true );
    constrainedDirections.set ( 6, true );
    constrainedDirections.set ( 7, true );
    constrainedDirections.set ( 8, true );
    constrainedDirections.set ( 9, true );
    constrainedDirections.set ( 12, true );
    ProjectorType boxProjector ( lowerBounds, upperBounds, constrainedDirections );
    aol::ConstrainedLevenbergMarquardtAlgorithm<RealType, MatrixType, LinearRegressionType, ProjectorType> levenbergMarquardtAlg ( atomData.size ( ), gaussNewtonF, gaussNewtonDF, boxProjector, 50, 1, 0.2, 0.8, 1e-6, 1e-6, 1e-6, _lmVerbose );
    aol::Vector<RealType> Arg ( numVar ), Dest ( numVar );
    Arg[0] = 0.5 * atomSize[0] - std::cos ( atomAngleRads ) * 0.5 * atomSeparation;    // Center 1 x
    Arg[1] = 0.5 * atomSize[1] - std::sin ( atomAngleRads ) * 0.5 * atomSeparation;    // Center 1 y
    Arg[2] = 0.5 * atomSize[0] + std::cos ( atomAngleRads ) * 0.5 * atomSeparation;    // Center 2 x
    Arg[3] = 0.5 * atomSize[1] + std::sin ( atomAngleRads ) * 0.5 * atomSeparation;    // Center 2 y
    Arg[4] = 0.5 * ( atomData.getMaxValue ( ) - atomData.getMinValue ( ) );            // Intensity 1
    Arg[5] = 0.5 * ( atomData.getMaxValue ( ) - atomData.getMinValue ( ) );            // Intensity 2
    Arg[6] = 0.25 * atomSize.getMinValue ( );                                          // Width 1
    Arg[7] = 0.25 * atomSize.getMinValue ( );                                          // Height 1
    Arg[8] = 0.25 * atomSize.getMinValue ( );                                          // Width 2
    Arg[9] = 0.25 * atomSize.getMinValue ( );                                          // Height 2
    Arg[10] = 0;                                                                       // Rotation 1
    Arg[11] = 0;                                                                       // Rotation 2
    Arg[12] = atomData.getMinValue ( );                                                // Offset
    if ( atomOffset.getMinValue ( ) > 0 )
      levenbergMarquardtAlg.apply ( Arg, Dest );
    else
      Dest = Arg;
    Dest[0] += corner[0];
    Dest[1] += corner[1];
    Dest[2] += corner[0];
    Dest[3] += corner[1];
    Centers[2*k][0] = Dest[0];
    Centers[2*k][1] = Dest[1];
    Centers[2*k+1][0] = Dest[2];
    Centers[2*k+1][1] = Dest[3];
    
    for ( int j=0; j<numVar ; ++j )
      GaussianParams[k][j] = Dest[j];
    
    approxAtomRadius += Dest[6] + Dest[7] + Dest[8] + Dest[9];
    
    if ( _progressBar != NULL ) (*_progressBar)++;
  }
  if ( _progressBar != NULL ) _progressBar->finish ( );
  
  approxAtomRadius /= 4 * ApproximateCenters.numComponents ( );
  // Clear boundary atoms
  if ( ClearBoundaryAtoms ) {
    short k = 0;
    while ( k < Centers.numComponents ( ) - 1 ) {
      if ( Centers[k][0] - approxAtomRadius < 0 || Centers[k][0] + approxAtomRadius >= Data.getNumX ( )
          || Centers[k][1] - approxAtomRadius < 0 || Centers[k][1] + approxAtomRadius >= Data.getNumY ( )
          || Centers[k+1][0] - approxAtomRadius < 0 || Centers[k+1][0] + approxAtomRadius >= Data.getNumX ( )
          || Centers[k+1][1] - approxAtomRadius < 0 || Centers[k+1][1] + approxAtomRadius >= Data.getNumY ( ) ) {
        Centers.eraseComponent ( k );
        Centers.eraseComponent ( k ); // component k+1 is now component k
        GaussianParams.eraseComponent ( k / 2 );
      } else k += 2;
    }
  }
  
  if ( !_quietMode )
    std::cerr << "# of (inner) dumbbells after bump fitting: " << Centers.numComponents ( ) / 2 << std::endl;
  
  if ( _diskOutput && GaussianParams.numComponents ( ) > 0 ) {
    PictureType bumpFunctionImg ( Data, aol::STRUCT_COPY );
    getBumpFunctionImage ( GaussianParams, bumpFunctionImg );
    bumpFunctionImg.save ( aol::strprintf ( "%s/bumpFuncs_dumbbell_refined%s", _outputDir.c_str ( ), getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    
    ColoredPictureType atomCentersImg ( Data.getNumX ( ), Data.getNumY ( ) );
    getAtomCentersImage ( Data, Centers, atomCentersImg );
    atomCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    atomCentersImg.savePNG ( aol::strprintf ( "%s/centers_dumbbell_refined.png", _outputDir.c_str ( ) ).c_str ( ) );
  }
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
                _ColoredPictureType>::readCentersFromCSV ( aol::MultiVector<RealType> &Centers, const std::string &Path ) const {
  std::ifstream file ( Path.c_str() );
  std::string value;
  std::list<std::string> values;
  while ( file.good() ) {
    getline ( file, value, ',' );
    if ( value.find('\n') != std::string::npos ) {
      size_t pos = 0;
      while ( ( pos = value.find ( "\n", pos + 1 ) ) != std::string::npos ) {
        std::string p = value.substr(0, pos);
        values.push_back(p);
        value = value.substr(pos + 1);
      }
      if (!value.empty()) {
        values.push_back(value);
      }
    } else {
      values.push_back(value);
    }
  }
  for ( int i=0; i<2; ++i ) // skip header
    values.pop_front ( );
  RealType xPos = 0, yPos = 0, xMax = 0, yMax = 0;
  int k = 0;
  for ( std::list<std::string>::const_iterator it = values.begin ( ); it != values.end ( ) ; it++ ) {
    const RealType val = strtod ( (*it).c_str(), NULL );
    if ( k % 2 == 0 ) {
      xPos = val;
      if ( xPos > xMax )
        xMax = xPos;
    }
    if ( k % 2 == 1 ) {
      yPos = val;
      if ( yPos > yMax )
        yMax = yPos;
      
      if ( !_quietMode )
        std::cerr << xPos << "; " << yPos << std::endl;
      
      Centers.resize ( Centers.numComponents ( ) + 1, 2 );
      Centers[Centers.numComponents ( )-1][0] = xPos;
      Centers[Centers.numComponents ( )-1][1] = yPos;
    }
    ++k;
  }
  file.close ( );
  
  if ( _diskOutput ) {
    short nx = aol::Rint ( xMax * 1.1 ), ny = aol::Rint ( yMax * 1.1 );
    ColoredPictureType atomCentersImg ( nx, ny );
    PictureType data ( nx, ny );
    getAtomCentersImage ( data, Centers, atomCentersImg );
    atomCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    atomCentersImg.savePNG ( aol::strprintf ( "%s/atomCenters.png", _outputDir.c_str ( ) ).c_str ( ) );
  }
}

template class AtomFinder<double, aol::FullMatrix<double>, aol::LinearRegressionQR<double>,
                          qc::ScalarArray<double, qc::QC_2D>, qc::MultiArray<double, qc::QC_2D, 3> >;

} // end namespace