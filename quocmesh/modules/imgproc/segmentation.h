#ifndef __SEGMENTATION_H
#define __SEGMENTATION_H

#include <quoc.h>
#include <finiteDifferences.h>
#ifdef USE_MODULES_QT
#include <customPlotHandler.h>
#endif

namespace im {

/**
 * \brief Thresholds an image and determines the threshold automatically with the isodata algorithm.
 *
 * \author Berkels
 */
template <typename RealType>
void isodataThreshold ( qc::ScalarArray<double, qc::QC_2D> &quocArray, const int numBins ) {
  aol::Vector<int> histo;
  quocArray.createHistogramOfValues ( histo, numBins );

  RealType threshold = 0.5 * numBins;
  int iter = 0;
  do {
    ++iter;
    RealType oldThreshold = threshold;

    int thresholdInt = aol::Rint ( threshold );
    threshold = static_cast<RealType> ( qc::computeFirstMoment ( histo, 0, thresholdInt - 1 ) ) / qc::computeMean ( histo, 0, thresholdInt - 1 );
    threshold += static_cast<RealType> ( qc::computeFirstMoment ( histo, thresholdInt, numBins - 1 ) ) / qc::computeMean ( histo, thresholdInt - 1, numBins - 1 );
    threshold *= 0.5;

    if ( ( aol::Abs ( threshold - oldThreshold ) < 0.5 ) || ( iter > 1000 ) )
      break;
  } while ( true );

  // The threshold was determined on the histogram and thus needs to be rescaled to the original range of the input image.
  const aol::Vec2<RealType> minMax = quocArray.getMinMaxValue();
  const RealType scaledThreshold = ( threshold / numBins ) * ( minMax[1] - minMax[0] ) + minMax[0];

  cerr << "Determined " << scaledThreshold << " as threshold in " << iter << " iterations\n";

  quocArray.threshold ( scaledThreshold, 0, 255 );
}

/**
 * Saves a ScalarArray<int, qc::QC_2D> as colored PNG, where non-negative integer entries are converted to gray scale,
 * negative entries are printed in blue and non-finite entries are printed in yellow.
 *
 * \author Mevenkamp
 */
template <typename RealType>
void saveColoredSegmentation ( const qc::ScalarArray<int, qc::QC_2D> &Segmentation, const std::string &baseFileName ) {
  qc::MultiArray<RealType, qc::QC_2D, 3> coloredSeg ( Segmentation.getNumX ( ), Segmentation.getNumY ( ) );
  const int maxVal = Segmentation.getMaxValue ( );
  for ( int i=0; i<Segmentation.size ( ) ; ++i ) {
    if ( aol::isFinite<RealType> ( Segmentation[i] ) ) {
      if ( Segmentation[i] >= 0 ) for ( int d=0; d<3; ++d ) coloredSeg[d][i] = Segmentation[i];
      else coloredSeg[2][i] = maxVal;
    } else {
      coloredSeg[0][i] = maxVal;
      coloredSeg[1][i] = maxVal;
    }
  }
  coloredSeg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, maxVal );
  coloredSeg.savePNG ( aol::strprintf ( "%s.png", baseFileName.c_str ( ) ).c_str ( ) );
}

/**
 * \brief Abstract base class for Mumford-Shah segmentation
 *
 * The interface function generateIndicatorFunction, in which the region indicators \f$ f_i \f$ are
 * defined, needs to be implemented in the derived class.
 *
 * \author Mevenkamp
 * \ingroup Segmentation
 */
template <typename ConfiguratorType>
class MSSegmentorBase : public qc::TVAlgorithmBase<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
protected:
  int _numSegments;
  std::vector<ArrayType*> _pIndicator;
  RealType _tau;
  bool _unknownRegion;
  std::string _outputDir;
public:
  MSSegmentorBase ( const typename ConfiguratorType::InitType &Initializer,
                    const RealType Gamma,
                    const int NumSegments,
                    const int MaxIterations,
                    const RealType StopEpsilon,
                    const std::string OutputDir = "" )
  : qc::TVAlgorithmBase<ConfiguratorType> ( Initializer, Gamma, MaxIterations, StopEpsilon ),
    _numSegments ( NumSegments ),
    _pIndicator ( NumSegments, NULL ),
    _tau ( 0.25 ),
    _unknownRegion ( false ),
    _outputDir ( OutputDir ) { }
  
  virtual ~MSSegmentorBase () {}
protected:
  //! Interface for a generic region indicator function
  virtual void generateIndicatorFunction ( const int IndicatorNumber, ArrayType &IndicatorFunction ) const {
    if ( _pIndicator[IndicatorNumber] == NULL )
      throw ( aol::Exception ( aol::strprintf ( "Indicator %d not set and generateIndicatorFunction() not overloaded.", IndicatorNumber ), __FILE__, __LINE__ ) );
    
    IndicatorFunction = *( _pIndicator[IndicatorNumber] );
  }
  
  virtual void generateIndicatorFunctions ( aol::VectorContainer<ArrayType> &IndicatorFunctions ) const {
    const short unknownRegionOffset = ( _unknownRegion ? 1 : 0 );
    IndicatorFunctions.reallocate ( this->_numSegments + unknownRegionOffset, ArrayType ( this->_grid ) );
    
    for ( int l=0; l<_numSegments ; ++l ) generateIndicatorFunction ( l, IndicatorFunctions[l + unknownRegionOffset] );
    if ( _unknownRegion ) addUnknownRegionConstantIndicatorFunction ( IndicatorFunctions );
    
    if ( this->_outputDir != "" ) {
      for ( int l=0; l<IndicatorFunctions.size ( ) ; ++l )
        IndicatorFunctions[l].save ( aol::strprintf ( "%s/indicator%d%s", _outputDir.c_str ( ), l, qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    }
  }
public:
  //! Indicator function is set point-wise to minimum of all indicator functions
  static void generateMinimumIndicatorFunction ( ArrayType &MinIndicatorFunction,
                                                 const aol::VectorContainer<ArrayType> &IndicatorFunctions ) {
    const int N = IndicatorFunctions[0].size ( );
    if ( MinIndicatorFunction.size ( ) != N ) throw aol::Exception ( "Indicator sizes do not match!", __FILE__, __LINE__ );
    aol::Vector<RealType> indicators ( IndicatorFunctions.size ( ) );
    for ( int i=0; i<N ; ++i ) {
      for ( int l=0; l<IndicatorFunctions.size ( ); ++l ) indicators[l] = IndicatorFunctions[l][i];
      MinIndicatorFunction[i] = indicators.getMinValue ( );
    }
  }
  
  //! Indicator function is set to constant equal to 0.5 * global max of point-wise minimum of all indicator functions
  static void addUnknownRegionConstantIndicatorFunction ( aol::VectorContainer<ArrayType> &IndicatorFunctions,
                                                          const int MinSize = 25 ) {
//    IndicatorFunctions[0].setAll ( aol::NumberTrait<RealType>::Inf );
//    ArrayType minIndicator ( IndicatorFunctions[0], aol::STRUCT_COPY );
//    generateMinimumIndicatorFunction ( minIndicator, IndicatorFunctions );
//    IndicatorFunctions[0].setAll ( 0.5 * minIndicator.getMaxValue ( ) );
    
    aol::Vector<RealType> indicatorThresholds ( IndicatorFunctions.size ( ) );
    for ( int l=0; l<IndicatorFunctions.size ( ) ; ++l ) {
      aol::Vector<RealType> valsInSmallRegionAroundMinVal;
      std::pair<int, RealType> minIndVal = IndicatorFunctions[l].getMinIndexAndValue ( );
      qc::FastILexMapper<ConfiguratorType::Dim> mapper ( IndicatorFunctions[l] );
      qc::CoordType center = mapper.splitGlobalIndex ( minIndVal.first );
      const short size = floor ( sqrt ( static_cast<RealType> ( MinSize ) ) ), radius = floor ( 0.5 * ( size - 1 ) );
      for ( qc::LocalLInfBoxIterator<ConfiguratorType::Dim> it ( center, radius, aol::Vec3<int> ( 0 ), IndicatorFunctions[l].getSize ( ) ); it.notAtEnd ( ) ; ++it )
        valsInSmallRegionAroundMinVal.pushBack ( IndicatorFunctions[l].get ( *it ) );
      indicatorThresholds[l] = valsInSmallRegionAroundMinVal.getMaxValue ( );
    }
    IndicatorFunctions[0].setAll ( indicatorThresholds.getMaxValue ( ) );
  }
  
  //! Indicator function is set to constant equal to ( 1.0 + Delta ) * max. of all known region indicators
  //! over indices specified by ThresholdSubset
  static void addUnknownRegionConstantIndicatorFunction ( aol::VectorContainer<ArrayType> &IndicatorFunctions,
                                                          const aol::Vector<int> &ThresholdSubset,
                                                          const RealType Delta = 0.0 ) {
    RealType max = -aol::NumberTrait<RealType>::Inf;
    for ( int i=0; i<ThresholdSubset.size ( ) ; ++i )
      for ( int l=1; l<IndicatorFunctions.size ( ); ++l )
        max = aol::Max<RealType> ( max, IndicatorFunctions[l][ThresholdSubset[i]] );
    IndicatorFunctions[0].setAll ( max * ( 1.0 + Delta ) );
  }
  
  //! Returns a connected region exceeding a certain size (e.g. 100 px) where the minimum of all indicators is peaked
  //! This is intended to be used for the creation of a new indicator for that region
  void getMaxMinIndicatorConnectedRegion ( aol::Vector<int> &IndicesOfNewRegion,
                                           const aol::VectorContainer<ArrayType> &IndicatorFunctions,
                                           const int MinSize = 100 ) const {
    IndicesOfNewRegion.reallocate ( 0 );
    
    ArrayType minIndicator ( IndicatorFunctions[0], aol::STRUCT_COPY );
    generateMinimumIndicatorFunction ( minIndicator, IndicatorFunctions );
    
    std::pair<int, RealType> maxIndVal = minIndicator.getMaxIndexAndValue ( );
    qc::FastILexMapper<qc::QC_2D> mapper ( this->_grid );
    qc::CoordType center = mapper.splitGlobalIndex ( maxIndVal.first );
    const short size = floor ( sqrt ( MinSize ) ), offset = floor ( 0.5 * ( size - 1 ) );
    for ( int dy=-offset; dy<=offset ; ++dy )
      for ( int dx=-offset; dx<=offset ; ++dx )
        IndicesOfNewRegion.pushBack ( mapper.getGlobalIndex ( center[0] + dx, center[1] + dy ) );
    
//    qc::BitArray<qc::QC_2D> mask ( minIndicator.getSize ( ) );
//    qc::FastILexMapper<qc::QC_2D> mapper ( mask );
//    int size = 0, n = 0;
//    while ( size < MinSize && n < mask.size ( ) ) {
//      std::pair<int, RealType> maxIndVal = minIndicator.getMaxIndexAndValue ( );
//      mask.set ( mapper.splitGlobalIndex ( maxIndVal.first ), true );
//      minIndicator[maxIndVal.first] = 0;
//      if ( n >= MinSize ) {
//        qc::ComponentsCollection<RealType> componentsCollection ( mask );
//        size = componentsCollection.getLargestComponent ( ).size ( );
//      }
//      ++n;
//    }
//    
//    qc::ComponentsCollection<RealType> componentsCollection ( mask );
//    const qc::Component<RealType> &largestComp = componentsCollection.getLargestComponent ( );
//    for ( std::set<aol::Vec2<short> >::const_iterator it=largestComp.begin ( ); it != largestComp.end ( ) ; ++it )
//      IndicesOfNewRegion.pushBack ( mapper.getGlobalIndex ( (*it)[0], (*it)[1] ) );
//    
//    if ( this->_outputDir != "" ) {
//      qc::ScalarArray<int, qc::QC_2D> finalMask ( this->_grid );
//      for ( int i=0; i<mask.size ( ) ; ++i ) if ( mask[i] ) finalMask[i] = 1;
//      for ( int i=0; i<IndicesOfNewRegion.size ( ) ; ++i ) finalMask[IndicesOfNewRegion[i]] = 2;
//      finalMask.setOverflowHandlingToCurrentValueRange ( );
//      finalMask.savePNG ( aol::strprintf ( "%s/finalMask.png", this->_outputDir.c_str ( ) ).c_str ( ) );
//    }
  }
  
  //! Returns maximum indicator function value inside inner part of regions specified by HardSegmentation and ErosionSize
  static RealType getMaxIndicatorValueInsideKnownRegions ( const ArrayType &IndicatorFunction,
                                                           const qc::ScalarArray<int, qc::QC_2D> &HardSegmentation,
                                                           const int ErosionSize = 0 ) {
    ArrayType erodedIndicator ( IndicatorFunction, aol::STRUCT_COPY );
    erodedIndicator.setAll ( -aol::NumberTrait<RealType>::Inf );
    qc::RectangularGrid<qc::QC_2D> grid ( ( qc::GridSize<qc::QC_2D> ( HardSegmentation ) ) );
    qc::BitArray<qc::QC_2D> mask ( grid );
    for ( qc::RectangularIterator<qc::QC_2D> it ( grid ); it.notAtEnd ( ) ; ++it )
      if ( HardSegmentation.get ( *it ) != 0 ) mask.set ( *it, true );
    mask.erodeBy ( ErosionSize );
    for ( qc::RectangularIterator<qc::QC_2D> it ( grid ); it.notAtEnd ( ) ; ++it )
      if ( mask.get ( *it ) ) erodedIndicator.set ( *it, IndicatorFunction.get ( *it ) );
    return erodedIndicator.getMaxValue ( );
  }
  
  //! Returns size of the unknown region (specified by entries of hard segmentation with value zero)
  static int getUnknownRegionSize ( const qc::ScalarArray<int, qc::QC_2D> &HardSegmentation ) {
    qc::RectangularGrid<qc::QC_2D> grid ( ( qc::GridSize<qc::QC_2D> ( HardSegmentation ) ) );
    qc::BitArray<qc::QC_2D> mask ( grid );
    for ( qc::RectangularIterator<qc::QC_2D> it ( grid ); it.notAtEnd ( ) ; ++it )
      if ( HardSegmentation.get ( *it ) == 0 ) mask.set ( *it, true );
    return mask.numTrue ( );
  }
  
  //! Returns a hard segmentation where each entry is 0 if the corresponding soft segmentation is below a specified threshold, and 1 else
  static void getHardSegmentation ( qc::ScalarArray<int, ConfiguratorType::Dim> &HardSegmentation, const ArrayType &SoftSegmentation,
                                   const RealType Threshold = 0.5 ) {
    if ( SoftSegmentation.size ( ) != HardSegmentation.size ( ) )
      throw aol::Exception ( "Dimensions of soft segmentation and hard segmentation do not match!", __FILE__, __LINE__ );
    
    for ( int k=0; k<HardSegmentation.size ( ) ; ++k ) HardSegmentation[k] = ( SoftSegmentation[k] > Threshold ? 1 : 0 );
  }
  
  //! Returns a hard segmentation where each entry is set to index of the soft segmentation with the largest value at that position
  static void getHardSegmentation ( qc::ScalarArray<int, ConfiguratorType::Dim> &HardSegmentation, const aol::VectorContainer<ArrayType> &SoftSegmentation ) {
    if ( SoftSegmentation.size ( ) == 0 || SoftSegmentation[0].size ( ) != HardSegmentation.size ( ) )
      throw aol::Exception ( "Dimensions of components of soft segmentation and hard segmentation do not match!", __FILE__, __LINE__ );
    
    for ( int k=0; k<HardSegmentation.size ( ) ; ++k ) {
      int maxInd = 0;
      for ( int l=1; l<SoftSegmentation.size ( ) ; ++l ) {
        if ( SoftSegmentation[l][k] > SoftSegmentation[maxInd][k] )
          maxInd = l;
      }
      HardSegmentation[k] = maxInd;
    }
  }
public:
  void setTau( const RealType Tau ) {
    _tau = Tau;
  }
  
  void setIndicatorReference ( const int IndicatorNumber, ArrayType &IndicatorFunction ) {
    _pIndicator[IndicatorNumber] = &IndicatorFunction;
  }
  
  int getNumSegments ( ) const {
    return _numSegments;
  }
  
  virtual void setNumSegments ( const int NumSegments ) {
    _numSegments = NumSegments;
    _pIndicator.resize ( NumSegments, NULL );
  }
  
  void setUnknownRegion ( const bool UnknownRegion = true ) {
    _unknownRegion = UnknownRegion;
  }
  
  void setOutputDirectory ( const std::string OutputDir ) {
    _outputDir = OutputDir;
  }
};


/**
 * \brief Abstract base class for two phase Mumford Shah segmentation.
 *
 * Minimizes the quadratic Esedoglu like model
 * \f[ \min_{u} \gamma\int_\Omega g|\nabla u|dx + \int_\Omega f_1u^2+f_2(1-u)^2 dx, \f]
 * where \f$ f_1,f_2 \f$ are indicator functions for the two different regions to be segmented.
 *
 * Uses a finite difference scheme closely based on {An Algorithm
 * for Total Variation Minimization and Applications} by Antonin Chambolle:
 * \f[ p^{k+1} = \left(p^k+\tau h^2\nabla\frac{\mathrm{div}p^k-2f_2/\gamma}{2(f_1+f_2)}\right)
 *    /\left(1+\frac{\tau h^2}{g}\left|\nabla\frac{\mathrm{div}p^k-2f_2/\gamma}{2(f_1+f_2)}\right|\right) \f]
 * where \f$ u^k=\frac{2f_2-\gamma\mathrm{div}p^k}{2(f_1+f_2)} \f$,
 * \f$ h \f$ is the grid size, and \f$ \tau \f$ can be taken as \f$ \frac18 \f$
 * as shown by Chambolle.
 *
 * The edge weighting function \f$ g \f$ can be defined
 * in the derived class by implementing the interface function generateEdgeWeight.
 *
 * \author Berkels
 * \ingroup Segmentation
 */
template <typename ConfiguratorType>
class TwoPhaseMSSegmentor : public MSSegmentorBase<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
protected:
  virtual void prepareIndicatorFunctionGeneration ( ) const {}
  
  virtual void generateEdgeWeight ( ArrayType & edgeWeight ) const {
    edgeWeight.setAll ( aol::ZOTrait<RealType>::one );
  }
public:
  TwoPhaseMSSegmentor ( const typename ConfiguratorType::InitType &Initializer,
                        const RealType Gamma,
                        const std::string OutputDir = "" )
    : MSSegmentorBase<ConfiguratorType> ( Initializer, Gamma, 2, 1000, 0.01, OutputDir ) { }
  
  virtual ~TwoPhaseMSSegmentor () {}
  
  void calcPrimalFromDual ( const qc::MultiArray<RealType, ConfiguratorType::Dim> &Dual,
                           ArrayType &Primal,
                           const ArrayType &Indicator1Plus2,
                           const ArrayType &Indicator2 ) const {
    qc::calculateBackwardFDDivergence<RealType, ConfiguratorType::Dim> ( Dual, Primal );
    Primal *= -0.5 * this->_gamma / this->_grid.H();
    Primal += Indicator2;
    Primal /= Indicator1Plus2;
  }
  
  void segment ( ArrayType &Segmentation, qc::MultiArray<RealType, ConfiguratorType::Dim> * PDual = NULL ) const {
    prepareIndicatorFunctionGeneration ( );
    doSegment ( Segmentation, PDual );
  }
  
  virtual void setNumSegments ( const int NumSegments ) {
    if ( NumSegments > 2 ) throw aol::Exception ( "Number of segments has to be either 2 or 1 (adding an unknown segment in that case)!", __FILE__, __LINE__ );
    
    if ( NumSegments == 1 ) {
      if ( !this->_unknownRegion && !this->_quietMode ) std::cerr << "Number of segments set to 1 (expected two segments). Adding an unknown region instead." << std::endl;
      this->setUnknownRegion ( true );
    } else if ( NumSegments == 2 ) {
      if ( this->_unknownRegion && !this->_quietMode ) std::cerr << "Number of (known) segments set to 2 (max. handled by this class). Disabling unknown region detection." << std::endl;
      this->setUnknownRegion ( false );
    }
    
    MSSegmentorBase<ConfiguratorType>::setNumSegments ( NumSegments );
  }
  
  void setUnknownRegion ( const bool UnknownRegion = true ) {
    if ( UnknownRegion ) {
      if ( this->_numSegments == 2 && !this->_quietMode ) std::cerr << "Total number of regions must be two. Thus, setting number of (known) regions to 1." << std::endl;
      MSSegmentorBase<ConfiguratorType>::setNumSegments ( 1 );
    } else {
      if ( this->_numSegments == 1 && !this->_quietMode ) std::cerr << "Total number of regions must be two. Thus, setting number of (known) regions to 2." << std::endl;
      MSSegmentorBase<ConfiguratorType>::setNumSegments ( 2 );
    }
    MSSegmentorBase<ConfiguratorType>::setUnknownRegion ( UnknownRegion );
  }
private:
  //! Can assume that prepareIndicatorFunctionGeneration ( ) has just been called.
  virtual void doSegment ( ArrayType &Segmentation, qc::MultiArray<RealType, ConfiguratorType::Dim> *PDual ) const {
    aol::VectorContainer<ArrayType> indicators;
    this->generateIndicatorFunctions ( indicators );
    indicators[0] += indicators [1];
    ArrayType indicator2OverGamma ( indicators[1] );
    indicator2OverGamma *= static_cast<RealType> ( this->_grid.H() ) / this->_gamma;
    
    ArrayType edgeWeight ( this->_grid );
    generateEdgeWeight ( edgeWeight );
    
    qc::MultiArray<RealType, ConfiguratorType::Dim> pOld ( this->_grid );
    if ( PDual != NULL )
      pOld = *PDual;
    qc::MultiArray<RealType, ConfiguratorType::Dim> pNew ( this->_grid );
    ArrayType temp ( this->_grid );
    const int numPrimalDofs = temp.size();
    qc::MultiArray<RealType, ConfiguratorType::Dim> gradient ( this->_grid );
    
    aol::ProgressBar<> progressBar ( "Segmenting" );
    if ( !this->_quietMode ) {
      progressBar.start ( this->_maxIterations );
      progressBar.display ( cerr );
    }
    
    for ( int fixpointIterations = 0; fixpointIterations < this->_maxIterations; fixpointIterations++ ) {
      qc::calculateBackwardFDDivergence<RealType, ConfiguratorType::Dim> ( pOld, temp );
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int j = 0; j < numPrimalDofs; ++j )
        temp[j] = ( temp[j] * 0.5 - indicator2OverGamma[j] ) / indicators[0][j];
      qc::calculateForwardFDGradient<RealType, ConfiguratorType::Dim> ( temp, gradient );
      
      typename ConfiguratorType::VecType gradientVec;
#ifdef _OPENMP
#pragma omp parallel for firstprivate ( gradientVec )
#endif
      for ( int j = 0; j < numPrimalDofs; ++j ) {
        for ( int i = 0; i < ConfiguratorType::Dim; ++i )
          gradientVec[i] = gradient[i][j];
        const RealType gradientVecNorm = gradientVec.norm();
        for ( int i = 0; i < ConfiguratorType::Dim; ++i ) {
          pNew[i][j] = ( pOld[i][j] + this->_tau * gradientVec[i] ) / ( 1 + this->_tau * gradientVecNorm / edgeWeight[j] );
        }
      }
      pOld -= pNew;
      const RealType change = pOld.norm();
      pOld = pNew;
      
      if ( this->_pStepSaver ) {
        calcPrimalFromDual ( pNew, temp, indicators[0], indicators[1] );
        this->_pStepSaver->saveStep ( temp, fixpointIterations );
      }
      
      // If the change is small enough, we consider the gradient descent to have converged.
      if ( change < this->_stopEpsilon )
        break;
      
      if ( !this->_quietMode ) progressBar++;
    }
    if ( !this->_quietMode ) progressBar.finish();
    calcPrimalFromDual ( pOld,  Segmentation, indicators[0], indicators[1] );
    if ( PDual != NULL )
      *PDual = pOld;
  }
};

template <typename ConfiguratorType>
class WeightedTwoPhaseMSSegmentor : public TwoPhaseMSSegmentor<ConfiguratorType> {
  
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  
  WeightedTwoPhaseMSSegmentor ( const typename ConfiguratorType::InitType &grid,
                               const RealType gamma,
                               const ArrayType & edgeWeight )
  : TwoPhaseMSSegmentor<ConfiguratorType> ( grid, gamma ),
  _edgeWeight ( edgeWeight )
  {}
  
  virtual void generateEdgeWeight ( ArrayType & edgeWeight ) const {
    edgeWeight = _edgeWeight;
  }
  
protected:
  const ArrayType & _edgeWeight;
};


/**
 * This struct + static function construction is a workaround to C++ limitations.
 * It is used to allow template specialization, in this particular case for the ImageDimension template parameter
 *
 * \author Mevenkamp
 */
template <typename ConfiguratorType, int ImageDimension>
struct doUpdateGrayValues {
  typedef typename ConfiguratorType::RealType RealType;
  
  static void apply ( const aol::Vector<RealType> &CurrentSegmentation,
                     const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &ImageMVec,
                     aol::MultiVector<RealType> &MeanValues,
                     bool QuietMode ) {
    aol::IdentityFunction<RealType>  heavisideFunction;
    aol::ClassicalChanVeseMeanValueUpdater<ConfiguratorType, aol::IdentityFunction<RealType>, 1, ImageDimension> grayValueUpdater ( Grid, heavisideFunction, ImageMVec );
    
    aol::MultiVector<RealType> levelsetFunctions ( 1, CurrentSegmentation.size() );
    levelsetFunctions[0] = CurrentSegmentation;
    
    // The "shift" in generateIndicatorFunction causes a severe loss of contrast
    // in CurrentSegmentation (according to the theory the 0.5 levelset that we are
    // interested in should be unaffected though). To properly calculate the gray
    // values we need to account for this here: The values start in [0,1].
    // First, shift to [-0.5,0.5]
    levelsetFunctions.addToAll ( -0.5 );
    // Rescale such that -0.5 and/or 0.5 is attained (increasing contrast but leaving 0 unaffected).
    levelsetFunctions /= 2*levelsetFunctions.getMaxAbsValue();
    // Finally, shift back to [0,1].
    levelsetFunctions.addToAll ( 0.5 );
    
    // This is an alternative to the code above (neither always better nor always worse).
    //levelsetFunctions.threshold( 0.5, 0., 1. );
    
    aol::MultiVector<RealType> newGrayValuesMVec ( 2, ImageDimension );
    
    grayValueUpdater.update( levelsetFunctions, newGrayValuesMVec );
    
    for ( int j = 0; j < ImageDimension; j++ ) {
      MeanValues[0][j] = newGrayValuesMVec[0][j];
      MeanValues[1][j] = newGrayValuesMVec[1][j];
    }
    if ( !QuietMode ) cerr << MeanValues << endl;
  }
};

/**
 * This struct + static function construction is a workaround to C++ limitations.
 * It is used here to specialize the ImageDimension template parameter
 *
 * \author Mevenkamp
 */
template <typename ConfiguratorType>
struct doUpdateGrayValues<ConfiguratorType, 0> {
  typedef typename ConfiguratorType::RealType RealType;
  
  static void apply ( const aol::Vector<RealType> &CurrentSegmentation,
                     const typename ConfiguratorType::InitType &/*Grid*/, const aol::MultiVector<RealType> &ImageMVec,
                     aol::MultiVector<RealType> &MeanValues,
                     bool QuietMode ) {
    aol::Vector<RealType> u ( CurrentSegmentation );
    u.addToAll ( -0.5 );
    u /= 2 * u.getMaxAbsValue ( );
    u.addToAll ( 0.5 );
    
    const int imageDim = ImageMVec.numComponents ( );
    const int numPixels = ImageMVec[0].size ( );
    
    MeanValues.setZero ( );
    RealType c1Norm = 0, c2Norm = 0, uSqr, uMinusOneSqr;
    for ( int i = 0; i < numPixels ; ++i ) {
      uSqr = aol::Sqr<RealType> ( u[i] );
      uMinusOneSqr = aol::Sqr<RealType> ( 1 - u[i] );
      c1Norm += uSqr;
      c2Norm += uMinusOneSqr;
      for ( int j = 0; j < imageDim ; ++j ) {
        MeanValues[0][j] += uSqr * ImageMVec[j][i];
        MeanValues[1][j] += uMinusOneSqr * ImageMVec[j][i];
      }
    }
    if ( c1Norm > 0 ) MeanValues[0] /= c1Norm;
    if ( c2Norm > 0 ) MeanValues[1] /= c2Norm;
    
    if ( !QuietMode ) cerr << MeanValues << endl;
  }
};


template <typename ConfiguratorType> class FirstOrderPrimalDualTwoPhaseMSSegmentor;

/**
 * Piecewise constant two-phase Mumford Shah image segmentation with support for scalar
 * and vector values images, e.g. color images or vector fields.
 *
 * \author Berkels
 * \ingroup Segmentation
 */
template <typename ConfiguratorType, int ImageDimension, typename BaseClass = FirstOrderPrimalDualTwoPhaseMSSegmentor<ConfiguratorType> >
class PiecewiseConstantTwoPhaseMSSegmentor : public BaseClass {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
protected:
  aol::MultiVector<RealType> &_imageMVec;
  aol::MultiVector<RealType> _meanValues;
protected:
  int _outerIterations;
  aol::ProgressBar<> *_progressBar;
  bool _catchCtrlC;
  mutable sigfunc _previousCtrlCHandler;
  std::string _baseOutputDir;
public:
  PiecewiseConstantTwoPhaseMSSegmentor ( const typename ConfiguratorType::InitType &Initializer,
                                         const RealType Gamma,
                                         aol::MultiVector<RealType> &ImageMVec,
                                         const bool InitializeGrayValues = true,
                                         const bool Quiet = false,
                                         const std::string OutputDir = "" )
  : BaseClass ( Initializer, Gamma ),
    _imageMVec ( ImageMVec ),
    _meanValues ( 2, ImageMVec.numComponents ( ) ),
    _outerIterations ( 5 ),
    _catchCtrlC ( false ) {
    this->setQuietMode ( Quiet );
    this->setOutputDirectory ( OutputDir );
    if ( InitializeGrayValues ) initializeGrayValues ( );
  }
  
  virtual ~PiecewiseConstantTwoPhaseMSSegmentor() {}
  
  virtual void generateIndicatorFunction ( const int IndicatorNumber, ArrayType &IndicatorFunction ) const {
    if ( ( IndicatorNumber >= 2 ) || ( IndicatorNumber < 0 ) )
      throw ( aol::OutOfBoundsException ( "PiecewiseConstantTwoPhaseMSSegmentor only defines two indicator functions.", __FILE__, __LINE__ ) );
    
    const int imageDim = _imageMVec.numComponents ( );
    const RealType shift = 0.5;
    for ( int i = 0; i < IndicatorFunction.size(); ++i) {
      RealType indicator = 0.;
      for ( int j = 0; j < imageDim; j++ )
        // The mean values are assigned to the segments consistently to
        // aol::ClassicalChanVeseEnergyMulti (when used for binary segmentation)
        indicator += aol::Sqr( _imageMVec[j][i] - _meanValues[1-IndicatorNumber][j] );
      IndicatorFunction[i] = indicator + shift;
    }
  }
protected:
  virtual void initializeGrayValues ( ) {
    aol::KMeansClusterer<RealType> kMeansClusterer;
    aol::MultiVector<RealType> clusters;
    kMeansClusterer.apply ( _imageMVec, clusters, 2 );
    const int imageDim = _imageMVec.numComponents ( );
    for ( int j=0; j<imageDim ; ++j ) {
      _meanValues[0][j] = clusters[j][0];
      _meanValues[1][j] = clusters[j][1];
    }
    if ( !this->_quietMode ) cerr << this->_meanValues << endl;
  }
  
  virtual void updateGrayValues ( const ArrayType &Segmentation ) {
    doUpdateGrayValues<ConfiguratorType, ImageDimension>::apply ( Segmentation, this->_grid, _imageMVec, _meanValues, this->_quietMode );
  }
public:
  void segmentAndAdjustGrayValues ( ArrayType &Segmentation, qc::MultiArray<RealType, ConfiguratorType::Dim> * PDual = NULL ) {
    setCtrlCHandler ( );
    
    _baseOutputDir = this->_outputDir;
    for ( int it = 0; it < _outerIterations && !wantsInterrupt ( ) ; ++it ) {
      if ( _baseOutputDir != "" ) {
        this->_outputDir = aol::strprintf ( "%s/it%d", _baseOutputDir.c_str ( ), it );
        aol::makeDirectory ( this->_outputDir.c_str ( ), false );
      }
      
      this->segment ( Segmentation, PDual );
      updateGrayValues ( Segmentation );
    }
    this->_outputDir = _baseOutputDir;
    
    unsetCtrlCHandler ( );
  }
  
  void addUnknownRegion ( qc::ScalarArray<int, qc::QC_2D> &HardSegmentation ) {
    const std::string outputDir = this->_outputDir;
    const std::string unknownRegionSegOutputDir = ( outputDir != "" ) ? aol::strprintf ( "%s/unknownRegionSeg", outputDir.c_str ( ) ) : "";
    if ( unknownRegionSegOutputDir != "" ) aol::makeDirectory ( unknownRegionSegOutputDir.c_str ( ) );
    this->_outputDir = unknownRegionSegOutputDir;
    
    if ( this->_unknownRegion && !this->_quietMode ) std::cerr << "Disabling detection of unknown region DURING segmentation" << std::endl;
    this->setUnknownRegion ( false );
    
    aol::VectorContainer<ArrayType> indicators;
    this->generateIndicatorFunctions ( indicators );
    ArrayType minIndicator ( this->_grid );
    this->generateMinimumIndicatorFunction ( minIndicator, indicators );
    if ( this->_outputDir != "" ) minIndicator.save ( aol::strprintf ( "%s/minIndicator%s", this->_outputDir.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    
    minIndicator.scaleValuesTo01 ( );
    aol::MultiVector<RealType> imageMVec ( _imageMVec );
    aol::MultiVector<RealType> meanValues ( _meanValues );
    _imageMVec.reallocate ( 1, minIndicator.size ( ) );
    _imageMVec[0] = minIndicator;
    _meanValues.reallocate ( 2, 1 );
    _meanValues[0][0] = 0;
    _meanValues[1][0] = 1;
    
    ArrayType minIndicatorSegmentation ( this->_grid );
    this->segment ( minIndicatorSegmentation );
    
    _imageMVec.reallocate ( imageMVec );
    _imageMVec = imageMVec;
    _meanValues.reallocate ( meanValues );
    _meanValues = meanValues;
    
    qc::ScalarArray<int, qc::QC_2D> minIndicatorHardSegmentation ( this->_grid );
    this->getHardSegmentation ( minIndicatorHardSegmentation, minIndicatorSegmentation );
    if ( this->_outputDir != "" ) {
      minIndicatorSegmentation.save ( aol::strprintf ( "%s/softSeg%s", this->_outputDir.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
      minIndicatorHardSegmentation.savePNG ( aol::strprintf ( "%s/hardSeg.png", this->_outputDir.c_str ( ) ).c_str ( ) );
    }
    
    HardSegmentation.addToAll ( 1 );
    if ( minIndicatorHardSegmentation.getMinValue ( ) != minIndicatorHardSegmentation.getMaxValue ( ) ) {
      for ( int i=0; i<minIndicatorHardSegmentation.size ( ) ; ++i )
        if ( minIndicatorHardSegmentation[i] == 1 ) HardSegmentation[i] = 0;
    }
    
    this->_outputDir = outputDir;
  }
  
  const aol::MultiVector<RealType>& getMeanValuesReference () const {
    return _meanValues;
  }
  
  aol::MultiVector<RealType>& getMeanValuesReference () {
    return _meanValues;
  }
  
  void setMeanValues ( const aol::Vector<int> &Indices ) {
    for ( int l=0; l<Indices.size ( ) ; ++l )
      for ( int j=0; j<_imageMVec.numComponents ( ) ; ++j )
        _meanValues[l][j] = _imageMVec[j][Indices[l]];
  }
  
  void setMeanValuesFromInitialSegmentation ( const ArrayType &InitialSegmentation ) {
    updateGrayValues ( InitialSegmentation );
  }
  
  void setNumGhostCells ( const int /*NumGhostCells*/ ) { }
  
  void setOuterIterations ( const int OuterIterations ) {
    _outerIterations = OuterIterations;
  }
  
  void setCatchCtrlC ( bool catchCtrlC ) {
    _catchCtrlC = catchCtrlC;
  }
protected:
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



/**
 * \brief Abstract base class for multi phase Mumford Shah segmentation.
 *
 * Minimizes the multi-phase Mumford-Shah model
 * \f[ \min_{u} \gamma\int_\Omega g|\nabla u|dx + \sum_{l=1}^k \int_\Omega f_l u_l dx, \f]
 * where \f$ f_1,\dots,f_n \f$ are indicator functions for the different regions to be segmented.
 *
 * The interface function generateIndicatorFunctions, in which \f$ f_1,\dots,f_n \f$ are
 * defined, needs to be implemented in the derived class.
 *
 * \author Mevenkamp
 * \ingroup Segmentation
 */
template <typename ConfiguratorType>
class MultiPhaseMSSegmentor : public MSSegmentorBase<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
protected:
  virtual void prepareIndicatorFunctionGeneration ( ) const {}
public:
  MultiPhaseMSSegmentor ( const typename ConfiguratorType::InitType &Initializer,
                          const RealType Gamma,
                          const int NumSegments = 0,
                          const int /*NumGhostCells*/ = 0,
                          const std::string OutputDir = "" )
    : MSSegmentorBase<ConfiguratorType> ( Initializer, Gamma, NumSegments, 1000, 0.01, OutputDir ) { }
  
  virtual ~MultiPhaseMSSegmentor () {}
  
  void calcPrimalFromDual ( const aol::VectorContainer<qc::MultiArray<RealType, ConfiguratorType::Dim> > &/*Dual*/,
                           aol::VectorContainer<ArrayType> &/*Primal*/,
                           const aol::VectorContainer<ArrayType> &/*Indicators*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void segment ( aol::VectorContainer<ArrayType> &Segmentation, aol::VectorContainer<qc::MultiArray<RealType, ConfiguratorType::Dim> > * PDual = NULL ) const {
    prepareIndicatorFunctionGeneration ( );
    doSegment ( Segmentation, PDual );
  }
  
  void setNumGhostCells ( const int /*NumGhostCells*/ ) { }
private:
  //! Can assume that prepareIndicatorFunctionGeneration ( ) has just been called.
  virtual void doSegment ( aol::VectorContainer<ArrayType> &/*Segmentation*/, aol::VectorContainer<qc::MultiArray<RealType, ConfiguratorType::Dim> > * /*PDual*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};


template <typename ConfiguratorType> class FirstOrderPrimalDualMultiPhaseMSSegmentor;

/**
 * Piecewise constant multi-phase Mumford Shah image segmentation with support for scalar
 * and vector values images, e.g. color images or vector fields.
 *
 * \author Mevenkamp
 * \ingroup Segmentation
 */
template <typename ConfiguratorType, typename BaseClass = FirstOrderPrimalDualMultiPhaseMSSegmentor<ConfiguratorType> >
class PiecewiseConstantMultiPhaseMSSegmentor : public BaseClass {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType IndicatorFunctionType;
  typedef aol::VectorContainer<IndicatorFunctionType> ArrayType;
protected:
  aol::MultiVector<RealType> &_imageMVec;
  aol::MultiVector<RealType> _meanValues;
protected:
  int _outerIterations;
  aol::ProgressBar<> *_progressBar;
  bool _catchCtrlC;
  mutable sigfunc _previousCtrlCHandler;
  std::string _baseOutputDir;
public:
  PiecewiseConstantMultiPhaseMSSegmentor ( const typename ConfiguratorType::InitType &Initializer,
                                          const RealType Gamma,
                                          aol::MultiVector<RealType> &ImageMVec,
                                          const bool InitializeGrayValues = true,
                                          const bool Quiet = false,
                                          const std::string OutputDir = "",
                                          const int NumSegments = 0,
                                          const int NumGhostCells = 0 )
  : BaseClass ( Initializer, Gamma, NumSegments, NumGhostCells ),
    _imageMVec ( ImageMVec ),
    _meanValues ( this->_numSegments, ImageMVec.numComponents ( ) ),
    _outerIterations ( 1 ),
    _catchCtrlC ( false ) {
    this->setQuietMode ( Quiet );
    this->setOutputDirectory ( OutputDir );
    if ( InitializeGrayValues ) initializeGrayValues ( );
  }
  
  virtual ~PiecewiseConstantMultiPhaseMSSegmentor() {}
  
  virtual void generateIndicatorFunction ( const int IndicatorNumber, IndicatorFunctionType &IndicatorFunction ) const {
    if ( IndicatorNumber < 0 || IndicatorNumber >= this->_numSegments )
      throw aol::Exception ( "Passed indicator number does not match number of segments!", __FILE__, __LINE__ );
    
    const int imageDim = _imageMVec.numComponents ( );
    for ( int i = 0; i < IndicatorFunction.size(); ++i ) {
      RealType indicator = 0.;
      for ( int j = 0; j < imageDim; j++ )
        indicator += aol::Sqr( _imageMVec[j][i] - _meanValues[IndicatorNumber][j] );
      IndicatorFunction[i] = indicator;
    }
  }
protected:
  virtual void initializeGrayValues ( ) {
    aol::KMeansClusterer<RealType> kMeansClusterer;
    aol::MultiVector<RealType> clusters;
    kMeansClusterer.apply ( _imageMVec, clusters, this->_numSegments, aol::KMEANS_INIT_METHOD::MEC );
    const int imageDim = _imageMVec.numComponents ( );
    for ( int j=0; j<imageDim ; ++j )
      for ( int l=0; l<clusters[j].size ( ) ; ++l )
        this->_meanValues[l][j] = clusters[j][l];
    if ( !this->_quietMode ) cerr << this->_meanValues << endl;
  }
  virtual void updateGrayValues ( const ArrayType &CurrentSegmentation ) {
    const int imageDim = _imageMVec.numComponents ( );
    const int numPixels = _imageMVec[0].size ( );
    
    _meanValues.setZero ( );
    aol::Vector<RealType> cNorms ( this->_numSegments );
    for ( int l = 0; l < this->_numSegments ; ++l ) {
      for ( int i = 0; i < numPixels ; ++i ) {
        cNorms[l] += CurrentSegmentation[l][i];
        for ( int j = 0; j < imageDim ; ++j )
          this->_meanValues[l][j] += CurrentSegmentation[l][i] * _imageMVec[j][i];
      }
      if ( cNorms[l] > 0 ) _meanValues[l] /= cNorms[l];
    }
    
    if ( !this->_quietMode ) cerr << _meanValues << endl;
  }
  virtual void initializeSegmentation ( ArrayType &Segmentation ) const {
    // Best regularity, very low fidelity
    // Initialize all soft segmentations constant with equal values everywhere
    //    Segmentation.setAll ( 1.0 / static_cast<RealType> ( Segmentation.size ( ) ) );
    
    // Best fidelity, very irregular
    // In each pixel, set the soft segmentation with index of the mean value nearest to the input image to one and all others to zero
    Segmentation.setAll ( 0.0 );
    const int imageDim = _imageMVec.numComponents ( );
    const int numPixels = _imageMVec[0].size ( );
    for ( int i=0; i<numPixels ; ++i ) {
      RealType minDist = aol::NumberTrait<RealType>::Inf;
      int nearestMeanIdx = 0;
      for ( int l=0; l<Segmentation.size ( ) ; ++l ) {
        RealType dist = 0;
        for ( int j=0; j<imageDim ; ++j )
          dist += aol::Sqr<RealType> ( _imageMVec[j][i] - _meanValues[l][j] );
        if ( dist < minDist ) {
          minDist = dist;
          nearestMeanIdx = l;
        }
      }
      Segmentation[nearestMeanIdx][i] = 1.0;
    }
  }
public:
  void segmentAndAdjustGrayValues ( ArrayType &Segmentation, aol::VectorContainer<qc::MultiArray<RealType, ConfiguratorType::Dim> > * PDual = NULL ) {
    initializeSegmentation ( Segmentation );
    setCtrlCHandler ( );
    
    _baseOutputDir = this->_outputDir;
    for ( int it = 0; it < _outerIterations && !wantsInterrupt ( ) ; ++it ) {
      if ( _baseOutputDir != "" ) {
        this->_outputDir = aol::strprintf ( "%s/it%d", _baseOutputDir.c_str ( ), it );
        aol::makeDirectory ( this->_outputDir.c_str ( ), false );
      }
      
      this->segment ( Segmentation, PDual );
      this->updateGrayValues ( Segmentation );
      
      if ( !this->_quietMode && this->_outputDir != "" ) {
        for ( int l=0; l<Segmentation.size ( ) ; ++l )
          Segmentation[l].save ( aol::strprintf ( "%s/segmentation_seg%d%s", this->_outputDir.c_str ( ), l, qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
      }
    }
    this->_outputDir = _baseOutputDir;
    
    unsetCtrlCHandler ( );
  }
  
  void addUnknownRegion ( qc::ScalarArray<int, qc::QC_2D> &HardSegmentation ) {
    const std::string outputDir = this->_outputDir;
    const std::string unknownRegionSegOutputDir = ( outputDir != "" ) ? aol::strprintf ( "%s/unknownRegionSeg", outputDir.c_str ( ) ) : "";
    if ( unknownRegionSegOutputDir != "" ) aol::makeDirectory ( unknownRegionSegOutputDir.c_str ( ) );
    
    this->_outputDir = "";
    aol::VectorContainer<IndicatorFunctionType> indicators;
    this->generateIndicatorFunctions ( indicators );
    this->_outputDir = unknownRegionSegOutputDir;
    IndicatorFunctionType minIndicator ( this->_grid );
    this->generateMinimumIndicatorFunction ( minIndicator, indicators );
    if ( this->_outputDir != "" ) minIndicator.save ( aol::strprintf ( "%s/minIndicator%s", this->_outputDir.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    
    minIndicator.scaleValuesTo01 ( );
    qc::MultiArray<RealType, qc::QC_2D, 1> minIndicatorsArr ( minIndicator, aol::DEEP_COPY );
    IndicatorFunctionType minIndicatorSegmentation ( this->_grid );
    PiecewiseConstantTwoPhaseMSSegmentor<ConfiguratorType, 1> segmentor ( this->_grid, this->_gamma, minIndicatorsArr, true, true, unknownRegionSegOutputDir );
    segmentor.setQuietMode ( this->_quietMode );
    segmentor.setCatchCtrlC ( _catchCtrlC );
    segmentor.setMaxIterations ( this->_maxIterations );
    segmentor.setStopEpsilon ( this->_stopEpsilon );
    segmentor.getMeanValuesReference ( )[0][0] = 0;
    segmentor.getMeanValuesReference ( )[1][0] = 1;
    segmentor.segment ( minIndicatorSegmentation );
    qc::ScalarArray<int, qc::QC_2D> minIndicatorHardSegmentation ( this->_grid );
    segmentor.getHardSegmentation ( minIndicatorHardSegmentation, minIndicatorSegmentation );
    if ( this->_outputDir != "" ) {
      minIndicatorSegmentation.save ( aol::strprintf ( "%s/softSeg%s", this->_outputDir.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
      minIndicatorHardSegmentation.savePNG ( aol::strprintf ( "%s/hardSeg.png", this->_outputDir.c_str ( ) ).c_str ( ) );
    }
    
    HardSegmentation.addToAll ( 1 );
    if ( minIndicatorHardSegmentation.getMinValue ( ) != minIndicatorHardSegmentation.getMaxValue ( ) ) {
      for ( int i=0; i<minIndicatorHardSegmentation.size ( ) ; ++i )
        if ( minIndicatorHardSegmentation[i] == 1 ) HardSegmentation[i] = 0;
    }
    
    this->_outputDir = outputDir;
  }
  
  void getMeanValuesImage ( aol::MultiVector<RealType> &MeanValuesImage, const qc::ScalarArray<int, ConfiguratorType::Dim> &HardSegmentation ) const {
    const int imageDim = this->_imageMVec.numComponents ( );
    const int numPixels = this->_imageMVec[0].size ( );
    const int numSegments = HardSegmentation.getMaxValue ( ) + 1;
    
    if ( HardSegmentation.size ( ) != numPixels ) throw aol::Exception ( "Dimension of hard segmentation does not match initial input image!", __FILE__, __LINE__ );
    
    aol::MultiVector<RealType> means ( numSegments, _imageMVec.numComponents ( ) );
    aol::Vector<int> numPixelsPerSegment ( numSegments );
    for ( int i = 0; i < numPixels ; ++i ) {
      for ( int j = 0; j < imageDim; j++ )
        means[HardSegmentation[i]][j] += _imageMVec[j][i];
      ++numPixelsPerSegment[HardSegmentation[i]];
    }
    for ( int l = 0; l < numSegments ; ++l )
      means[l] /= static_cast<RealType> ( numPixelsPerSegment[l] );
    for ( int i = 0; i < numPixels ; ++i )
      for ( int j = 0; j < imageDim ; ++j )
        MeanValuesImage[j][i] = means[HardSegmentation[i]][j];
  }
  
  const aol::MultiVector<RealType>& getMeanValuesReference () const {
    return _meanValues;
  }
  
  aol::MultiVector<RealType>& getMeanValuesReference () {
    return _meanValues;
  }
  
  void setMeanValues ( const aol::Vector<int> &Indices ) {
    for ( int l=0; l<Indices.size ( ) ; ++l )
      for ( int j=0; j<_imageMVec.numComponents ( ) ; ++j )
        _meanValues[l][j] = _imageMVec[j][Indices[l]];
  }
  
  void setMeanValuesFromInitialSegmentation ( const ArrayType &InitialSegmentation ) {
    updateGrayValues ( InitialSegmentation );
  }
  
  void setNumGhostCells ( const int NumGhostCells ) {
    this->_numGhostCells = NumGhostCells;
  }
  
  void setOuterIterations ( const int OuterIterations ) {
    _outerIterations = OuterIterations;
  }
  
  void setCatchCtrlC ( bool catchCtrlC ) {
    _catchCtrlC = catchCtrlC;
  }
protected:
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


/**
 * Piecewise sub-space multi-phase Mumford Shah image segmentation with support for scalar
 * and vector values images, e.g. color images or vector fields.
 *
 * \author Mevenkamp
 * \ingroup Segmentation
 */
template <typename ConfiguratorType, typename BaseClass = FirstOrderPrimalDualMultiPhaseMSSegmentor<ConfiguratorType> >
class PiecewiseSubSpaceMultiPhaseMSSegmentor : public BaseClass {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType IndicatorFunctionType;
  typedef aol::VectorContainer<IndicatorFunctionType> ArrayType;
protected:
  aol::MultiVector<RealType> &_imageMVec;
  aol::MultiVector<RealType> _meanValues;
  aol::RandomAccessContainer<aol::MultiVector<RealType> > _subSpaces;
  const int _subSpaceDim;
protected:
  int _outerIterations;
  aol::ProgressBar<> *_progressBar;
  bool _catchCtrlC;
  mutable sigfunc _previousCtrlCHandler;
  std::string _baseOutputDir;
public:
  PiecewiseSubSpaceMultiPhaseMSSegmentor ( const typename ConfiguratorType::InitType &Initializer,
                                           const RealType Gamma,
                                           aol::MultiVector<RealType> &ImageMVec,
                                           const int NumSegments,
                                           const int SubSpaceDim = 3,
                                           const int NumGhostCells = 0 )
  : BaseClass ( Initializer, Gamma, NumSegments, NumGhostCells ),
    _imageMVec ( ImageMVec ),
    _meanValues ( NumSegments, ImageMVec.numComponents ( ) ),
    _subSpaces ( NumSegments, aol::MultiVector<RealType> ( SubSpaceDim, ImageMVec.numComponents ( ) ) ),
    _subSpaceDim ( SubSpaceDim ),
    _outerIterations ( 1 ),
    _catchCtrlC ( false ) { }
  
  virtual ~PiecewiseSubSpaceMultiPhaseMSSegmentor() {}
  
  virtual void generateIndicatorFunction ( const int IndicatorNumber, IndicatorFunctionType &IndicatorFunction ) const {
    const int imageDim = _imageMVec.numComponents ( );
    aol::Vector<RealType> v ( imageDim ), proj ( v );
    for ( int i = 0; i < IndicatorFunction.size ( ); ++i ) {
      // Calculate projection of v[i] onto subspace of segment l
      _imageMVec.getTo ( i, v );
      v -= _meanValues[IndicatorNumber];
      proj.setZero ( );
      for ( int k = 0; k < _subSpaceDim ; ++k ) proj.addMultiple ( _subSpaces[IndicatorNumber][k], v.dotProduct ( _subSpaces[IndicatorNumber][k] ) );
      v -= proj;
      IndicatorFunction[i] = v.normSqr ( );
      IndicatorFunction[i] = log ( 1 + IndicatorFunction[i] );
    }
  }
protected:
  virtual void updateSubSpaces ( const ArrayType &CurrentSegmentation ) {
    const int numPixels = _imageMVec[0].size ( );
    
    qc::ScalarArray<int, ConfiguratorType::Dim> hardSegmentation ( this->_grid );
    getHardSegmentation ( hardSegmentation, CurrentSegmentation, 0.25 );
    if ( this->_outputDir != "" ) im::saveColoredSegmentation<RealType> ( hardSegmentation, aol::strprintf ( "%s/hardSegmentation", this->_outputDir.c_str ( ) ) );
    
    for ( int l=0; l<this->_numSegments ; ++l ) {
      // Assemble points
      aol::Vector<int> points;
      for ( int i=0; i<numPixels ; ++i ) {
        if ( hardSegmentation[i] == l + ( this->_unknownRegion ? 1 : 0 ) )
          points.pushBack ( i );
      }
      if ( !this->_quietMode ) std::cerr << "#samples for sub-space #" << l << ": " << points.size ( ) << std::endl;
      
      getSubSpace ( points, _meanValues[l], _subSpaces[l] );
      
      if ( this->_outputDir != "" ) visualizeSubSpace ( l, points );
    }
  }
  
  virtual void visualizeSubSpace ( const int SubSpaceNumber, const aol::Vector<int> &Points = aol::Vector<int> ( ) ) const {
    if ( this->_outputDir == "" ) throw aol::Exception ( "Cannot store visualization because output directory was not specified!", __FILE__, __LINE__ );
    
#ifdef USE_MODULES_QT
    { // Plot mean spectrum and some exemplary spectra from subspace
      CustomPlotHandler<RealType> qcpHandler ( aol::strprintf ( "Spectra (region/subspace %d)", SubSpaceNumber ), true );
      if ( Points.size ( ) >= 3 ) {
        aol::Vector<RealType> v ( _imageMVec.numComponents ( ) );
        for ( int i=0; i<3 ; ++i ) {
          _imageMVec.getTo ( Points[i], v );
          v -= _meanValues[SubSpaceNumber];
          qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( v, aol::strprintf ( "example spectrum #%d", i ) );
        }
      }
      qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( _meanValues[SubSpaceNumber], "mean spectrum" );
      qcpHandler.saveToFile ( aol::strprintf ( "%s/spectra_seg%d.q1cp", this->_outputDir.c_str ( ), SubSpaceNumber ).c_str ( ) );
    }
    { // Plot eigen vectors of subspace
      CustomPlotHandler<RealType> qcpHandler ( aol::strprintf ( "Principal components of subspace %d", SubSpaceNumber ), true );
      for ( int k=0; k<_subSpaceDim ; ++k )
        qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( _subSpaces[SubSpaceNumber][k], aol::strprintf ( "comp %d", k ) );
      qcpHandler.saveToFile ( aol::strprintf ( "%s/pca_seg%d.q1cp", this->_outputDir.c_str ( ), SubSpaceNumber ).c_str ( ) );
    }
#endif
    
    qc::ScalarArray<int, qc::QC_2D> mask ( this->_grid );
    for ( int i=0; i<Points.size ( ) ; ++i ) mask[Points[i]] = 1;
    mask.setOverflowHandlingToCurrentValueRange ( );
    mask.savePNG ( aol::strprintf ( "%s/subSpaceMask_seg%d.png", this->_outputDir.c_str ( ), SubSpaceNumber ).c_str ( ) );
  }
  
  virtual void getSubSpace ( const aol::Vector<int> &Points,
                             aol::Vector<RealType> &MeanVal,
                             aol::MultiVector<RealType> &SubSpace ) const {
    // Assemble multi-vector of spectra extracted at given points
    const int imageDim = _imageMVec.numComponents ( );
    aol::MultiVector<RealType> segVecs ( imageDim, Points.size ( ) );
    aol::Vector<RealType> v ( imageDim );
    for ( int i=0; i<Points.size ( ) ; ++i ) {
        _imageMVec.getTo ( Points[i], v );
        segVecs.set ( i, v );
    }
    
    // Compute mean values
    segVecs.getMeanComponents ( MeanVal );
    
    // Compute PCA
    aol::Vector<RealType> eigenVals;
    aol::getPCAEigenVecs<RealType> ( SubSpace, eigenVals, segVecs, _subSpaceDim, 0.0, true );
  }
  
  virtual void initializeSegmentation ( ArrayType &Segmentation ) const {
    // Best regularity, very low fidelity
    // Initialize all soft segmentations constant with equal values everywhere
    Segmentation.setAll ( 1.0 / static_cast<RealType> ( Segmentation.size ( ) ) );
    
    // Best fidelity, very irregular
    // In each pixel, set the soft segmentation with index of the mean value nearest to the input image to one and all others to zero
    // TODO
  }
public:
  void segmentAndAdjustSubSpaces ( ArrayType &Segmentation, aol::VectorContainer<qc::MultiArray<RealType, ConfiguratorType::Dim> > * PDual = NULL,
                                  const bool InitializeSegmentation = false ) {
    if ( InitializeSegmentation ) initializeSegmentation ( Segmentation );
    setCtrlCHandler ( );
    
    _baseOutputDir = this->_outputDir;
    for ( int it = 0; it < _outerIterations && !wantsInterrupt ( ) ; ++it ) {
      if ( _baseOutputDir != "" ) {
        this->_outputDir = aol::strprintf ( "%s/it%d", _baseOutputDir.c_str ( ), it );
        aol::makeDirectory ( this->_outputDir.c_str ( ) );
      }
      
      this->segment ( Segmentation, PDual );
      this->updateSubSpaces ( Segmentation );
      
      if ( this->_outputDir != "" ) {
        for ( int l=0; l<Segmentation.size ( ) ; ++l )
          Segmentation[l].save ( aol::strprintf ( "%s/segmentation_seg%d%s", this->_outputDir.c_str ( ), l, qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
      }
    }
    this->_outputDir = _baseOutputDir;
    
    unsetCtrlCHandler ( );
  }
  
  void getHardSegmentation ( qc::ScalarArray<int, ConfiguratorType::Dim> &HardSegmentation, const ArrayType &SoftSegmentation, const RealType Threshold = 0.0 ) const {
    if ( SoftSegmentation.size ( ) == 0 || SoftSegmentation[0].size ( ) != HardSegmentation.size ( ) )
      throw aol::Exception ( "Dimensions of components of soft segmentation and hard segmentation do not match!", __FILE__, __LINE__ );
    
    for ( int k=0; k<SoftSegmentation[0].size ( ) ; ++k ) {
      int maxInd = -1;
      RealType maxVal = 0;
      for ( int l=0; l<SoftSegmentation.size ( ) ; ++l ) {
        if ( SoftSegmentation[l][k] >= maxVal + Threshold ) {
          maxInd = l;
          maxVal = SoftSegmentation[l][k];
        }
      }
      HardSegmentation[k] = maxInd;
    }
  }
  
  const aol::RandomAccessContainer<aol::MultiVector<RealType> >& getSubSpacesReference () const {
    return _subSpaces;
  }
  
  aol::RandomAccessContainer<aol::MultiVector<RealType> >& getSubSpacesReference () {
    return _subSpaces;
  }
  
  void setSubSpaces ( const aol::RandomAccessContainer<aol::MultiVector<RealType> > &SubSpaces ) {
    _subSpaces.clear ( );
    _subSpaces.pushBack ( SubSpaces );
  }
  
  void setSubSpacesViaSubSpaceClustering ( const bool Update = false ) {
    const std::string baseOutputDir = this->_outputDir;
    this->_outputDir = ( baseOutputDir != "" ) ? aol::strprintf ( "%s/initializeSubSpacesViaClustering", this->_outputDir.c_str ( ) ) : "";
    if ( this->_outputDir != "" ) aol::makeDirectory ( this->_outputDir.c_str ( ), false );
    
    bool unknownRegion = this->_unknownRegion;
    this->_unknownRegion = false;
    
    aol::Vector<int> clusterLabels;
    aol::KSubSpaceClusterer<RealType> subSpaceClusterer;
    subSpaceClusterer.setQuietMode ( this->_quietMode );
    subSpaceClusterer.setOutputDirectory ( this->_outputDir );
    subSpaceClusterer.setLabelDimensions2D ( this->_grid.getNumX ( ), this->_grid.getNumY ( ) );
    subSpaceClusterer.apply ( _imageMVec, this->_numSegments, _subSpaceDim, clusterLabels, aol::KSUBSPACE_INIT_METHOD::RNG, 1 );
    ArrayType segmentation ( this->_numSegments, IndicatorFunctionType ( this->_grid ) );
    for ( int i=0; i<clusterLabels.size ( ) ; ++i ) segmentation[clusterLabels[i]][i] = 1;
    updateSubSpaces ( segmentation );
    
    this->_unknownRegion = unknownRegion;
    
    if ( Update && this->_unknownRegion ) {
      const int outerIterations = _outerIterations;
      _outerIterations = 1;
      ArrayType segmentation ( this->getNumSegments ( ), IndicatorFunctionType ( this->_grid ) );
      segmentAndAdjustSubSpaces ( segmentation );
      _outerIterations = outerIterations;
    }
    
    this->_outputDir = baseOutputDir;
  }
  
  void setSubSpacesFromInitialSegmentation ( const ArrayType &InitialSegmentation, const bool Update = false ) {
    const std::string baseOutputDir = this->_outputDir;
    this->_outputDir = ( baseOutputDir != "" ) ? aol::strprintf ( "%s/initializeSubSpacesFromSegmentation", this->_outputDir.c_str ( ) ) : "";
    if ( this->_outputDir != "" ) aol::makeDirectory ( this->_outputDir.c_str ( ), false );
    
    bool unknownRegion = this->_unknownRegion;
    this->_unknownRegion = false;
    
    updateSubSpaces ( InitialSegmentation );
    
    this->_unknownRegion = unknownRegion;
    
    if ( Update && this->_unknownRegion ) {
      const int outerIterations = _outerIterations;
      _outerIterations = 1;
      ArrayType segmentation ( InitialSegmentation );
      segmentAndAdjustSubSpaces ( segmentation );
      _outerIterations = outerIterations;
    }
    
    this->_outputDir = baseOutputDir;
  }
  
  void addSubSpaceFromHighestMinIndicatorRegion ( const bool Update = false ) {
    const int numSegments = this->getNumSegments ( );
    const int imageDim = _imageMVec.numComponents ( );
    
    const std::string baseOutputDir = this->_outputDir;
    this->_outputDir = ( baseOutputDir != "" ) ? aol::strprintf ( "%s/addSubSpaceFromHighestMinIndiactorRegion_%d", this->_outputDir.c_str ( ), numSegments ) : "";
    if ( this->_outputDir != "" ) aol::makeDirectory ( this->_outputDir.c_str ( ), false );
    
    bool unknownRegion = this->_unknownRegion;
    this->_unknownRegion = false;
    
    ArrayType indicatorFunctions;
    this->generateIndicatorFunctions ( indicatorFunctions );
    aol::Vector<int> indicesOfNewRegion;
    this->getMaxMinIndicatorConnectedRegion ( indicesOfNewRegion, indicatorFunctions );
    aol::Vector<RealType> *newMeanVal = new aol::Vector<RealType> ( imageDim );
    aol::MultiVector<RealType> newSubSpace ( _subSpaceDim, imageDim );
    getSubSpace ( indicesOfNewRegion, *newMeanVal, newSubSpace );
    _meanValues.appendReference ( *newMeanVal, true );
    _subSpaces.pushBack ( newSubSpace );
    BaseClass::setNumSegments ( numSegments + 1 );
    
    if ( !this->_quietMode ) std::cerr << "#samples for new sub-space: " << indicesOfNewRegion.size ( ) << std::endl;
    
    if ( this->_outputDir != "" ) {
      IndicatorFunctionType minIndicator ( this->_grid );
      this->generateMinimumIndicatorFunction ( minIndicator, indicatorFunctions );
      minIndicator.save ( aol::strprintf ( "%s/minIndicator_before%s", this->_outputDir.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
      this->generateIndicatorFunctions ( indicatorFunctions );
      this->generateMinimumIndicatorFunction ( minIndicator, indicatorFunctions );
      minIndicator.save ( aol::strprintf ( "%s/minIndicator_after%s", this->_outputDir.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
      visualizeSubSpace ( numSegments, indicesOfNewRegion );
    }
    
    this->_unknownRegion = unknownRegion;
    
    if ( Update && this->_unknownRegion ) {
      const int outerIterations = _outerIterations;
      _outerIterations = 1;
      ArrayType segmentation ( this->getNumSegments ( ), IndicatorFunctionType ( this->_grid ) );
      segmentAndAdjustSubSpaces ( segmentation );
      _outerIterations = outerIterations;
    }
    
    this->_outputDir = baseOutputDir;
  }
  
  void setNumSegments ( const int NumSegments ) {
    BaseClass::setNumSegments ( NumSegments );
    
    _meanValues.resize ( NumSegments, this->_imageMVec.numComponents ( ) );
    _subSpaces.reallocate ( NumSegments );
    for ( int l=0; l<NumSegments ; ++l ) _subSpaces[l].reallocate ( _subSpaceDim, this->_imageMVec.numComponents ( ) );
  }
  
  void setNumGhostCells ( const int NumGhostCells ) {
    this->_numGhostCells = NumGhostCells;
  }
  
  void setOuterIterations ( const int OuterIterations ) {
    _outerIterations = OuterIterations;
  }
  
  void setCatchCtrlC ( bool catchCtrlC ) {
    _catchCtrlC = catchCtrlC;
  }
protected:
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
  
} // namespace im

#endif
