#ifndef __FEATURESEGMENTATION_H
#define __FEATURESEGMENTATION_H

#include <stack>

#include <quoc.h>
#include <configurators.h>
#include <convolution.h>
#include <segmentation.h>
#include <primalDualSegmentation.h>
#include <ctrlCCatcher.h>
#include <progressBar.h>
#include <gnuplotter.h>
#include <connectedComponents.h>
#include <clustering.h>
#ifdef USE_LIB_HUNGARIAN
#include <Hungarian.h>
#endif

namespace im {
  
template <typename RealType>
static inline void cropSegment ( qc::ScalarArray<RealType, qc::QC_2D> &ArgDest, const qc::ScalarArray<int, qc::QC_2D> &Segmentation, const int SegmentIdx ) {
  const int nx = ArgDest.getNumX ( ), ny = ArgDest.getNumY ( );
  
  if ( Segmentation.getNumX ( ) != nx || Segmentation.getNumY ( ) != ny ) throw aol::Exception ( "Input and segmentation dimensions do not match!", __FILE__, __LINE__ );
  if ( SegmentIdx < 0 || SegmentIdx > Segmentation.getMaxValue ( ) ) throw aol::Exception ( "Segment index invalid!", __FILE__, __LINE__ );
  
  // Determine crop start & end
  aol::Vec<2,int> cropStart; cropStart[0] = -1; cropStart[1] = -1;
  aol::Vec<2,int> cropEnd; cropEnd[0] = nx; cropEnd[1] = ny;
  qc::ScalarArray<int, qc::QC_1D> lineX ( nx ), lineY ( ny );
  do { ++cropStart[0]; Segmentation.getLine ( qc::QC_X, cropStart[0], lineY ); } while ( cropStart[0] < nx-1 && lineY.numOccurence ( SegmentIdx ) == 0 );
  do { ++cropStart[1]; Segmentation.getLine ( qc::QC_Y, cropStart[1], lineX ); } while ( cropStart[1] < ny-1 && lineX.numOccurence ( SegmentIdx ) == 0 );
  do { --cropEnd[0]; Segmentation.getLine ( qc::QC_X, cropEnd[0], lineY ); } while ( cropEnd[0] > 0 && lineY.numOccurence ( SegmentIdx ) == 0 );
  do { --cropEnd[1]; Segmentation.getLine ( qc::QC_Y, cropEnd[1], lineX ); } while ( cropEnd[1] > 0 && lineX.numOccurence ( SegmentIdx ) == 0 );

  aol::Vec<2,int> cropSize; cropSize[0] = cropEnd[0] - cropStart[0] + 1; cropSize[1] = cropEnd[1] - cropStart[1] + 1;
  ArgDest.crop ( cropStart, cropSize );
}
  
template <typename RealType>
static inline bool borderContainsNonFinite ( const qc::ScalarArray<RealType, qc::QC_2D> &Data, const qc::CoordType &Center, const int Radius ) {
  RealType sum = 0;
  for ( int x=Center[0]-Radius; x<=Center[0]+Radius ; ++x ) sum += Data.get ( x, Center[1]-Radius );
  for ( int x=Center[0]-Radius; x<=Center[0]+Radius ; ++x ) sum += Data.get ( x, Center[1]+Radius );
  for ( int y=Center[1]-Radius+1; y<=Center[1]+Radius-1 ; ++y ) sum += Data.get ( Center[0]-Radius, y );
  for ( int y=Center[1]-Radius+1; y<=Center[1]+Radius-1 ; ++y ) sum += Data.get ( Center[0]+Radius, y );
  return !aol::isFinite<RealType> ( sum );
}
  
template <typename RealType>
static inline RealType getMaxSquareRadius ( const qc::ScalarArray<RealType, qc::QC_2D> &Data, const qc::CoordType &Center ) {
  int maxRadius = aol::Min<int> ( aol::Min<int> ( aol::Min<int> ( Center[0], Center[1] ), Data.getNumX ( ) - Center[0] - 1 ), Data.getNumY ( ) - Center[1] - 1 );
  while ( borderContainsNonFinite ( Data, Center, maxRadius ) ) --maxRadius;
  return maxRadius;
}
  
template <typename RealType>
static inline void rotateCropToSquareAndResampleSegment ( qc::ScalarArray<RealType, qc::QC_2D> &ArgDest, RealType &Angle ) {
  Angle = 0; // Rotation and resampling currently not supported
  
  // Determine crop start & end
  const int nx = ArgDest.getNumX ( ), ny = ArgDest.getNumY ( );
  const int minRadius = 10;
  aol::Vec<2,int> cropStart;
  int maxRadius = 0;
  for ( int y=minRadius; y<ny-minRadius ; ++y ) {
    for ( int x=minRadius; x<nx-minRadius ; ++x ) {
      const int radius = getMaxSquareRadius ( ArgDest, qc::CoordType ( x, y ) );
      if ( radius > maxRadius ) {
        maxRadius = radius;
        cropStart[0] = x;
        cropStart[1] = y;
      }
    }
  }
  cropStart[0] -= maxRadius; cropStart[1] -= maxRadius;
  aol::Vec<2,int> cropSize; cropSize[0] = 2 * maxRadius + 1; cropSize[1] = 2 * maxRadius + 1;
  cropSize[0] = aol::Min<int> ( cropSize[0], cropSize[1] );
  cropSize[1] = aol::Min<int> ( cropSize[0], cropSize[1] );
  ArgDest.crop ( cropStart, cropSize );
}

static inline float getMaxMinSideArea ( const int A, const int B ) {
  const int min = aol::Min<int> ( A, B );
  const int max = aol::Max<int> ( A, B );
  return aol::Sqr<float> ( min + ( 1 - 1.0 / static_cast<float> ( max ) ) );
}
  
static inline void getMaxMinSideHistogram ( aol::Vector<int> &Hist, int &Start, aol::Vec<2,int> &Size ) {
  Start = 0;
  Size.setZero ( );
  float maxArea = 0;
  std::stack<std::pair<int, int> > stack;
  
  int x = 0;
  for ( x = 0; x < Hist.size ( ); x++ ) {
    int start = x;
    int height = Hist[x];
    while ( true ) {
      if ( stack.empty ( ) || height > stack.top ( ).second ) {
        stack.push ( std::pair<int, int> ( start, height ) );
      } else if ( height < stack.top ( ).second ) {
        const float tempArea = getMaxMinSideArea ( stack.top ( ).second, x - stack.top ( ).first );
        if ( tempArea > maxArea ) {
          Size[0] = x - stack.top ( ).first;
          Size[1] = stack.top ( ).second;
          Start = stack.top ( ).first;
          maxArea = tempArea;
        }
        
        std::pair<int, int> popped = stack.top ( );
        stack.pop ( );
        start = popped.first;
        continue;
      }
      
      break;
    }
  }
  
  while ( !stack.empty ( ) ) {
    std::pair<int, int> data = stack.top ( );
    stack.pop ( );
    const float tempArea = getMaxMinSideArea ( data.second, x - data.first );
    if ( tempArea > maxArea ) {
      Size[0] = x - data.first;
      Size[1] = data.second;
      Start = data.first;
      maxArea = tempArea;
    }
  }
}

static inline void getMaxMinSideRectangle ( qc::BitArray<qc::QC_2D> &Mask, aol::Vec<2,int> &Start, aol::Vec<2,int> &Size ) {
  // STEP 1: build a seed histogram using the first row of grid points
  aol::Vector<int> hist ( Mask.getNumX ( ) );
  for ( int x=0; x<hist.size ( ) ; ++x )
    hist[x] = ( Mask.get ( x, 0 ) ) ? 1 : 0;
  
  // STEP 2: get a starting max area from the seed histogram we created above.
  int start = 0;
  aol::Vec<2,int> size;
  getMaxMinSideHistogram ( hist, start, size );
  int maxArea = getMaxMinSideArea ( size[0], size[1] );
  Start[0] = start;
  Start[1] = 0;
  Size = size;
  
  // STEP 3: build histograms for each additional row, re-testing for new possible max rectangluar areas
  for ( int y = 1; y < Mask.getNumY ( ) ; ++y ) {
    // build a new histogram for this row. the values of this row are
    // 0 if the current grid point is occupied; otherwise, it is 1 + the value
    // of the previously found historgram value for the previous position.
    // What this does is effectly keep track of the height of continous avilable spaces.
    for  ( int x = 0; x < Mask.getNumX ( ) ; ++x )
      hist[x] = ( Mask.get ( x, y ) ) ? hist[x] + 1 : 0;
      
    // find the maximum size of the current histogram. If it happens to be larger
    // than the currently recorded max size, then it is the new max size.
    getMaxMinSideHistogram ( hist, start, size );
    
    const float tempArea = getMaxMinSideArea ( size[0], size[1] );
    if ( tempArea > maxArea ) {
      maxArea = tempArea;
      Size = size;
      Start[0] = start;
      Start[1] = y - size[1] + 1;
    }
  }
}
  
template <typename RealType>
static inline void rotateAndResampleFromLargestRectangle ( qc::ScalarArray<RealType, qc::QC_2D> &ArgDest, RealType &Angle ) {
  Angle = 0;
  aol::Vec<2,int> cropStart, cropSize;
  
  if ( !aol::isFinite<RealType> ( ArgDest.sum ( ) ) ) {
    qc::BitArray<qc::QC_2D> mask ( ArgDest.getSize ( ) );
    for ( int y=0; y<ArgDest.getNumY ( ) ; ++y )
      for ( int x=0; x<ArgDest.getNumX ( ) ; ++x )
        mask.set ( x, y, aol::isFinite<RealType> ( ArgDest.get ( x, y ) ) );
    
    getMaxMinSideRectangle ( mask, cropStart, cropSize );
  } else {
    cropSize[0] = ArgDest.getNumX ( );
    cropSize[1] = ArgDest.getNumY ( );
  }
  
  // Make crop size odd
  for ( int d=0; d<2 ; ++d ) {
    if ( cropSize[d] % 2 == 0 ) --cropSize[d];
  }
  ArgDest.crop ( cropStart, cropSize );
}

static inline void transformGrayValuesToLabels ( qc::ScalarArray<int, qc::QC_2D> &Segmentation ) {
  int label = -1;
  while ( Segmentation.getMaxValue ( ) >= 0 ) {
    int val = Segmentation.getMaxValue ( );
    for ( int i=0; i<Segmentation.size ( ) ; ++i ) {
      if ( Segmentation[i] == val )
        Segmentation[i] = label;
    }
    --label;
  }
  Segmentation *= -1;
  Segmentation.addToAll ( -1 );
}

static inline void swapLabels ( qc::ScalarArray<int, qc::QC_2D> &Segmentation, const int l1, const int l2 ) {
  qc::ScalarArray<int, qc::QC_2D> segmentation ( Segmentation );
  for ( int i=0; i<Segmentation.size ( ) ; ++i ) {
    if ( Segmentation[i] == l1 ) segmentation[i] = l2;
    if ( Segmentation[i] == l2 ) segmentation[i] = l1;
  }
  Segmentation = segmentation;
}

template <typename RealType, typename PictureType>
static inline void plotSegmentationBoundariesOntoImage ( const char* path, const PictureType &Input, const qc::ScalarArray<int, qc::QC_2D> &Segmentation ) {
  typedef qc::RectangularGridConfigurator<RealType, qc::QC_2D, aol::GaussQuadrature<RealType,qc::QC_2D,3> > ConfiguratorType;
  
  qc::ScalarArray<int, qc::QC_2D> seg ( Segmentation );
  transformGrayValuesToLabels ( seg );
  aol::Plotter<RealType> plotter;
  aol::PlotDataFileHandler<RealType> plotDataFileHandler;
  plotDataFileHandler.generateBackgroundPNGData ( Input );
  plotter.addPlotCommandsFromHandler ( plotDataFileHandler );
  
  aol::MultiVector<RealType> levelsetFunctions ( seg.getMaxValue ( ) + 1, seg.size ( ) );
  for ( int i=0; i<seg.size ( ) ; ++i ) levelsetFunctions[seg[i]][i] = 1.0;
  typename ConfiguratorType::InitType grid ( qc::GridSize<qc::QC_2D> ( Input.getNumX ( ), Input.getNumY ( ) ) );
  ConfiguratorType confDummy ( grid );
  plotDataFileHandler.generateIsolineData ( levelsetFunctions, grid, confDummy, 0.5, false, true );
  const RealType rgb[3] = {1, 0, 0};
  for ( unsigned int i=1; i<plotDataFileHandler.getDataFileNames().size ( ) ; ++i )
    plotter.addLinePlot ( plotDataFileHandler.getDataFileNames ( )[i].c_str ( ), NULL, 3, rgb );
  
  plotter.set_outfile_base_name ( path );
  plotter.setBackgroundPNGSettings ( Input.getNumX ( ), Input.getNumY ( ) );
  plotter.genPlot ( aol::GNUPLOT_PDF );
}

/**
 *  \brief Hungarian Method (Kuhn-Munkres algorithm) for finding a perfect (edge weight minimizing) matching in a bipartite graph
 *
 *  Input: adjacency matrix of a bipartite graph (edge weights)
 *         use the value NEG_INF to designate non-existing edges
 *  Output: a vector that assigns each node in the left part of the graph (index of the vector entry) at most one node in the right part of the graph (value of the vector entry)
 *          vector entries equal to -1 indicate that the corresponding node in the left part of the graph is not connected to any node in the right part of the graph
 *
 *  This is just a wrapper for the Hungarian library found at: http://www.frc.ri.cmu.edu/~lantao/codes/hungarian.php
 *
 *  Lantao Liu, Dylan Shell. "Assessing Optimal Assignment under Uncertainty: An Interval-based Algorithm".
 *  International Journal of Robotics Research (IJRR). vol. 30, no. 7, pp 936-953. Jun 2011.
 *
 *  \author mevenkamp
 *  \ingroup Segmentation
 */
template <typename RealType>
static inline void HungarianMethod ( const aol::FullMatrix<RealType> &AdjacencyMatrix, aol::Vector<int> &Memberships ) {
#ifdef USE_LIB_HUNGARIAN
  const int rows = AdjacencyMatrix.getNumRows ( ), cols = AdjacencyMatrix.getNumCols ( );
  const int k = aol::Max<int> ( rows, cols );
  HungarianMatrix m ( k, std::vector<Edge> ( k ) );
  for ( int i=0; i<k ; ++i )
    for ( int j=0; j<k ; ++j )
      m[i][i].SetWeight ( NEG_INF );
  for ( int i=0; i<rows; ++i )
    for ( int j=0; j<cols ; ++j )
      m[i][j].SetWeight ( -AdjacencyMatrix.get ( i, j ) );

  BipartiteGraph bg(m);
  Hungarian h(bg);
  h.HungarianAlgo();

  Memberships.setAll ( -1 );
  BipartiteGraph* bgFinal = h.GetBG();
  for(unsigned int i=0; i<bgFinal->GetNumAgents(); i++){
    for(unsigned int j=0; j<bgFinal->GetNumTasks(); j++){
      if(bgFinal->GetMatrix(i,j)->GetMatchedFlag() && m[i][j].GetWeight ( ) != NEG_INF ) Memberships[i] = j;
    }
  }
#else
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( AdjacencyMatrix );
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( Memberships );
  throw aol::Exception ( "Hungarian method requires library hungarian. Compile with -BUILD_AND_USE_HUNGARIAN", __FILE__, __LINE__ );
#endif
}

/**
 * \brief Provides different measures for comparing a 2D segmentation (integer label in each pixel) to ground truth.
 * \author mevenkamp
 * \ingroup Segmentation
 */
template <typename RealType>
class SegmentationQualityQuantifier {
  typedef qc::ScalarArray<int, qc::QC_2D> PictureType;
protected:
  PictureType _gt, _seg;
  const int _numPixels;
  int _N, _M, _K;
  aol::FullMatrix<int> _errorMatrix;
  aol::Vector<int> _iHat, _i;
  const RealType _overlapAcceptanceThreshold;
public:
  SegmentationQualityQuantifier ( const PictureType &GroundTruth, const PictureType &Segmentation )
    : _gt ( GroundTruth, aol::DEEP_COPY ), _seg ( Segmentation, aol::DEEP_COPY ), _numPixels ( GroundTruth.size ( ) ),
      _overlapAcceptanceThreshold ( 0.75 ) {

    transformGrayValuesToLabels ( _gt );
    transformGrayValuesToLabels ( _seg );
      
    _N = _gt.getMaxValue ( ) + 1;
    _M = _seg.getMaxValue ( ) + 1;
    _K = aol::Max<int> ( _N, _M );
        
    // Compute error matrix (n_{i,j}) (i.e. number of pixels interpreted as i-th class, but belonging to j-th class
    _errorMatrix.reallocate ( _K, _K );
    for ( int i=0; i<_M ; ++i ) {
      for ( int j=0; j<_N ; ++j ) {
        for ( int k=0; k<_numPixels ; ++k ) {
          if ( _seg[k] == i && _gt[k] == j )
            _errorMatrix.add ( i, j, 1 );
        }
      }
    }
        
    // Compute memberships (i.e. mapping of ground truth segments to classified labels)
    aol::FullMatrix<int> adjacencyMatrix ( _N, _M );
    for ( int gtLabel=0; gtLabel<_N ; ++gtLabel ) {
      for ( int segLabel=0; segLabel<_M ; ++segLabel ) {
        for ( int k=0; k<_numPixels ; ++k ) {
          if ( ( _gt[k] == gtLabel && _seg[k] != segLabel ) || ( _seg[k] == segLabel && _gt[k] != gtLabel ) )
            adjacencyMatrix.add ( gtLabel, segLabel, 1 );
        }
      }
    }
    aol::Vector<int> memberships ( _K );
    HungarianMethod ( adjacencyMatrix, memberships );

    _iHat.reallocate ( _K ), _i.reallocate ( _K );
    for ( int i=0; i<_K ; ++i ) {
      if ( memberships[i] >= 0 ) {
        _iHat[i] = memberships[i];
        _i[_iHat[i]] = i;
      }
    }
  }
  
  void saveSegmentationsWithMatchingLabels ( const char* outputdir ) {
    PictureType seg ( _seg );
    for ( int k=0; k<_numPixels ; ++k ) {
      if ( _i[_seg[k]] < _N )
        seg[k] = _i[_seg[k]];
    }
    
    _gt.setOverflowHandlingToCurrentValueRange ( );
    _gt.savePNG ( aol::strprintf ( "%s/gt.png", outputdir ).c_str ( ) );
    
    seg.setOverflowHandlingToCurrentValueRange ( );
    seg.savePNG ( aol::strprintf ( "%s/seg.png", outputdir ).c_str ( ) );
  }
  
  void printReport ( const char* path ) {
    std::ofstream txtFile ( path );
    txtFile << "Report on segmentation quality" << std::endl << std::endl;
    
    txtFile << "Pixel-wise criteria" << std::endl;
    txtFile << "CS  = " << getCorrectDetection ( ) * 100.0 << "%" << std::endl;
    txtFile << "O   = " << getOmissionError ( ) * 100.0 << "%" << std::endl;
    txtFile << "C   = " << getCommissionError ( ) * 100.0 << "%" << std::endl;
    txtFile << "CA  = " << getClassAccuracy ( ) * 100 << "%" << std::endl;
    txtFile << "CO  = " << getCorrectAssignment ( ) * 100 << "%" << std::endl;
    txtFile << "CC  = " << getObjectAccuracy ( ) * 100 << "%" << std::endl;
    txtFile << std::endl;
  }
  
  /* Begin: region-wise criteria */
  
  RealType getCorrectDetection ( ) const {
    RealType cs = 0;
    for ( int i=0; i<_N ; ++i ) {
      const int numOverlap = getNumOverlap ( i );
      if ( numOverlap >= _overlapAcceptanceThreshold * getNIDot ( _iHat[i] ) && numOverlap >= _overlapAcceptanceThreshold * getNDotI ( i ) )
        ++cs;
    }
    cs /= static_cast<RealType> ( _N );
    return cs;
  }
  
  /* End: region-wise criteria */
  
  /* Begin: pixel-wise criteria */
  
  RealType getOmissionError ( ) const {
    aol::Vector<RealType> omissionErrors ( _N );
    for ( int i=0; i<_N ; ++i ) {
      RealType ndoti = getNDotI ( i );
      omissionErrors[i] = ( ndoti - _errorMatrix.get ( _iHat[i], i ) ) / static_cast<RealType> ( ndoti );
    }
    return omissionErrors.getMedianValue ( );
  }
  
  RealType getCommissionError ( ) const {
    aol::Vector<RealType> commissionErrors ( _M );
    for ( int iHat=0; iHat<_M ; ++iHat ) {
      const RealType nidot = getNIDot ( iHat );
      commissionErrors[iHat] = ( nidot - _errorMatrix.get ( iHat, _i[iHat] ) ) / static_cast<RealType> ( nidot );
    }
    return commissionErrors.getMedianValue ( );
  }
  
  RealType getClassAccuracy ( ) const {
    RealType ca = 0;
    for ( int i=0; i<_K ; ++i ) {
      const RealType ndoti = getNDotI ( i );
      const RealType nidot = getNIDot ( _iHat[i] );
      ca += _errorMatrix.get ( _iHat[i], i ) * ndoti / ( ndoti + nidot - _errorMatrix.get ( _iHat[i], i ) );
    }
    return ca / static_cast<RealType> ( _numPixels );
  }
  
  RealType getCorrectAssignment ( ) const {
    RealType co = 0;
    for ( int i=0; i<_K ; ++i ) co += _errorMatrix.get ( _iHat[i], i );
    return co / static_cast<RealType> ( _numPixels );
  }
  
  RealType getObjectAccuracy ( ) const {
    RealType cc = 0;
    for ( int i=0; i<_K ; ++i ) {
      const RealType nidot = getNIDot ( _iHat[i] );
      cc += ( nidot > 0 ) ? _errorMatrix.get ( _iHat[i], i ) * getNDotI ( i ) / nidot : 0.0;
    }
    return cc / static_cast<RealType> ( _numPixels );
  }
  
  /* End: pixel-wise criteria */
protected:
  RealType getNIDot ( const int I ) const {
    RealType nidot = 0;
    for ( int j=0; j<_N ; ++j ) nidot += _errorMatrix.get ( I, j );
    return nidot;
  }
  
  RealType getNDotI ( const int I ) const {
    RealType ndoti = 0;
    for ( int j=0; j<_M ; ++j ) ndoti += _errorMatrix.get ( j, I );
    return ndoti;
  }
  
  int getNumOverlap ( const int I ) const {
    int numOverlap = 0;
    for ( int k=0; k<_numPixels ; ++k ) {
      if ( _gt[k] == I && _seg[k] == _iHat[I] )
        ++numOverlap;
    }
    return numOverlap;
  }
};

/*
 * BEGIN: Operators suitable for segmentation
 */

/**
 * \brief Converts a 2D image into a MultiVector of Discrete Orthogonal S-Transformed windows (real and imaginary part are stored)
 * \author mevenkamp
 * \ingroup Segmentation
 */
template <typename _RealType, typename _PictureType>
class DOSTOp : public aol::Op<_PictureType, aol::MultiVector<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const int _size, _sizeSqr;
  const int _blockOffset;
public:
  DOSTOp ( const int Size )
    : _size ( Size ), _sizeSqr ( Size * Size ),
      _blockOffset ( ( Size - 1 ) / 2 ) { }
  
  void apply ( const PictureType &Arg, aol::MultiVector<RealType> &Dest ) const {
    Dest.resize ( 2 * _sizeSqr, Arg.size ( ) );
    qc::FastILexMapper<qc::QC_2D> mapper ( Arg.getNumX ( ), Arg.getNumY ( ) );
    qc::MultiArray<RealType, 2, 2> function ( _size, _size ), transform ( _size, _size );
    for ( int y=0; y<Arg.getNumY ( ) ; ++y ) {
      for ( int x=0; x<Arg.getNumX ( ) ; ++x ) {
        for ( int dy=-_blockOffset ; dy<=_blockOffset ; ++dy )
          for ( int dx=-_blockOffset ; dx<=_blockOffset ; ++dx )
            function[0].set ( dx + _blockOffset, dy + _blockOffset, Arg.getReflection ( x + dx, y + dy ) );
        qc::StockwellTransform<RealType> ( function, transform );
        transform /= static_cast<RealType> ( _sizeSqr );
        for ( int i=0; i<_sizeSqr ; ++i )
          for ( int k=0; k<2 ; ++k )
            Dest[k * _sizeSqr + i][mapper.getGlobalIndex ( x, y )] = transform[k][i];
      }
    }
  }
  
  virtual void applyAdd ( const PictureType &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};

/**
 * \brief Converts a 2D image into a MultiVector of the moduli of Discrete Orthogonal S-Transformed windows
 * \author mevenkamp
 * \ingroup Segmentation
 */
template <typename _RealType, typename _PictureType>
class DOSTModulusOp : public aol::Op<_PictureType, aol::MultiVector<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const int _size, _sizeSqr;
  const int _blockOffset;
public:
  DOSTModulusOp ( const int Size )
    : _size ( Size ), _sizeSqr ( Size * Size ),
      _blockOffset ( ( Size - 1 ) / 2 ) { }
  
  void apply ( const PictureType &Arg, aol::MultiVector<RealType> &Dest ) const {
    Dest.resize ( _sizeSqr, Arg.size ( ) );
    qc::FastILexMapper<qc::QC_2D> mapper ( Arg.getNumX ( ), Arg.getNumY ( ) );
    qc::MultiArray<RealType, 2, 2> function ( _size, _size ), transform ( _size, _size );
    for ( int y=0; y<Arg.getNumY ( ) ; ++y ) {
      for ( int x=0; x<Arg.getNumX ( ) ; ++x ) {
        for ( int dy=-_blockOffset ; dy<=_blockOffset ; ++dy )
          for ( int dx=-_blockOffset ; dx<=_blockOffset ; ++dx )
            function[0].set ( dx + _blockOffset, dy + _blockOffset, Arg.getReflection ( x + dx, y + dy ) );
        qc::StockwellTransform<RealType> ( function, transform );
        transform /= static_cast<RealType> ( _sizeSqr );
        for ( int i=0; i<_sizeSqr ; ++i )
          Dest[i][mapper.getGlobalIndex ( x, y )] = sqrt ( aol::Sqr<RealType> ( transform[0][i] ) + aol::Sqr<RealType> ( transform[1][i] ) );
      }
    }
  }
  
  virtual void applyAdd ( const PictureType &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};

/**
 * \brief Converts a 2D image into a MultiVector of the moduli of Fourier transformed windows
 * \author mevenkamp
 * \ingroup Segmentation
 */
template <typename _RealType, typename _PictureType>
class FFTModulusOp : public aol::Op<_PictureType, aol::MultiVector<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const int _size, _sizeSqr;
  const int _blockOffset;
public:
  FFTModulusOp ( const int Size )
    : _size ( Size ), _sizeSqr ( Size * Size ),
      _blockOffset ( ( Size - 1 ) / 2 ) { }
  
  void apply ( const PictureType &Arg, aol::MultiVector<RealType> &Dest ) const {
    Dest.resize ( _sizeSqr, Arg.size ( ) );
    qc::FastILexMapper<qc::QC_2D> mapper ( Arg.getNumX ( ), Arg.getNumY ( ) );
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int y=0; y<Arg.getNumY ( ) ; ++y ) {
      for ( int x=0; x<Arg.getNumX ( ) ; ++x ) {
        PictureType block ( _size ), modulus ( _size );
        for ( int dy=-_blockOffset ; dy<=_blockOffset ; ++dy )
          for ( int dx=-_blockOffset ; dx<=_blockOffset ; ++dx )
            block.set ( dx + _blockOffset, dy + _blockOffset, Arg.getReflection ( x + dx, y + dy ) );
        qc::computeLogFFTModulus<RealType> ( block, modulus, 0, false );
        modulus /= static_cast<RealType> ( _sizeSqr );
        for ( int k=0; k<modulus.size ( ) ; ++k )
          Dest[k][mapper.getGlobalIndex ( x, y )] = modulus[k];
      }
    }
  }
  
  virtual void applyAdd ( const PictureType &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};

/**
 * \brief Abstract class for the convolution of a 2D image with a filter kernel. Requires implementation of the filter kernel matrix.
 * \author mevenkamp
 * \ingroup Segmentation
 */
template <typename _PictureType>
class FilterConvolutionOp : public aol::Op<_PictureType> {
  typedef _PictureType PictureType;
  typedef typename PictureType::RealType RealType;
protected:
  const int _size, _offset;
  aol::FullMatrix<RealType> _M;
public:
  FilterConvolutionOp ( const int Size ) : _size ( Size ), _offset ( ( Size - 1 ) / 2 ), _M ( Size, Size ) { }
  
  virtual void apply ( const PictureType &Arg, PictureType &Dest ) const {
    Dest.setZero ( );
    for ( int y=0; y<Arg.getNumY ( ) ; ++y )
      for ( int x=0; x<Arg.getNumX ( ) ; ++x )
        for ( int dy=-_offset; dy<=_offset ; ++dy )
          for ( int dx=-_offset; dx<=_offset ; ++dx )
            Dest.add ( x, y, _M.get ( dy + _offset , dx + _offset ) * Arg.getReflection ( x + dx, y + dy ) );
  }
  
  virtual void applyAdd ( const PictureType &/*Arg*/, PictureType &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  const aol::FullMatrix<RealType>& getFilterKernel ( ) const {
    return _M;
  }
};

/**
 * \author mevenkamp
 * \ingroup Segmentation
 */
template <typename _RealType, typename _PictureType> class FilterResponseLocalHistogramOp;

/**
 * \brief Converts a 2D scalar image into a MultiVector of histograms of filtered windows.
 * \author mevenkamp
 * \ingroup Segmentation
 */
template <typename _RealType>
class FilterResponseLocalHistogramOp<_RealType, qc::ScalarArray<_RealType, qc::QC_2D> > : public aol::Op<qc::ScalarArray<_RealType, qc::QC_2D>, aol::MultiVector<_RealType> > {
  typedef _RealType RealType;
  typedef qc::ScalarArray<_RealType, qc::QC_2D> PictureType;
protected:
  const FilterConvolutionOp<PictureType> *_filter;
  const int _integrationScale, _integrationScaleSqr, _offset;
  const int _numBins;
  RealType _minVal, _maxVal, _numBinsInvMinMaxDiff;
public:
  FilterResponseLocalHistogramOp ( const FilterConvolutionOp<PictureType> *Filter, const int IntegrationScale, const int NumBins )
    : _filter ( Filter ),
      _integrationScale ( IntegrationScale ),
      _integrationScaleSqr ( aol::Sqr<RealType> ( IntegrationScale ) ),
      _offset ( ( IntegrationScale - 1 ) / 2 ), _numBins ( NumBins ) { }
  
  void apply ( const PictureType &Arg, aol::MultiVector<RealType> &Dest ) const {
    Dest.resize ( _numBins, Arg.size ( ) );
    qc::FastILexMapper<qc::QC_2D> mapper ( Arg.getNumX ( ), Arg.getNumY ( ) );
    PictureType response ( Arg, aol::STRUCT_COPY );
    _filter->apply ( Arg, response );
    RealType minVal = response.getMinValue ( ), maxVal = response.getMaxValue ( );
    RealType numBinsInvMinMaxDiff = static_cast<RealType> ( _numBins ) / ( maxVal - minVal );
    aol::Vector<RealType> histogram ( _numBins );
    for ( int y=0; y<response.getNumY ( ) ; ++y ) {
      for ( int x=0; x<response.getNumX ( ) ; ++x ) {
        histogram.setZero ( );
        int binIdx;
        for ( int dy=-_offset ; dy<=_offset ; ++dy ) {
          for ( int dx=-_offset ; dx<=_offset ; ++dx ) {
            binIdx = floor ( ( response.getClip ( x + dx, y + dy ) - minVal ) * numBinsInvMinMaxDiff );
            if ( binIdx < 0 ) binIdx = 0;
            if ( binIdx >= _numBins ) binIdx = _numBins - 1;
            histogram[binIdx] += 1.0;
          }
        }
        histogram /= static_cast<RealType> ( _integrationScaleSqr );
        for ( int k=0; k<_numBins ; ++k )
          Dest[k][mapper.getGlobalIndex ( x, y )] = histogram[k];
      }
    }
  }
  
  int getNumBins ( ) const {
    return _numBins;
  }
  
  virtual void applyAdd ( const PictureType &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};

/**
 * \brief Converts a 2D multi-dimensional image (e.g. color) into a MultiVector of histograms of filtered windows.
 * \author mevenkamp
 * \ingroup Segmentation
 */
template <typename _RealType>
class FilterResponseLocalHistogramOp<_RealType, qc::MultiArray<_RealType, qc::QC_2D, 3> > : public aol::Op<qc::MultiArray<_RealType, qc::QC_2D, 3>, aol::MultiVector<_RealType> > {
  typedef _RealType RealType;
  typedef qc::MultiArray<RealType, qc::QC_2D, 3> PictureType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> ScalarPictureType;
protected:
  const FilterResponseLocalHistogramOp<RealType, ScalarPictureType> _scalarFilterResponceLocalHistogramOp;
  const bool _convertToGray;
public:
  FilterResponseLocalHistogramOp ( const FilterConvolutionOp<ScalarPictureType> *Filter, const int IntegrationScale, const int NumBins,
                                   const bool ConvertToGray = true )
    : _scalarFilterResponceLocalHistogramOp ( Filter, IntegrationScale, NumBins ),
      _convertToGray ( ConvertToGray ) { }
  
  void apply ( const PictureType &Arg, aol::MultiVector<RealType> &Dest ) const {
    const int numBins = _scalarFilterResponceLocalHistogramOp.getNumBins ( );
    if ( _convertToGray ) {
      ScalarPictureType gray ( Arg[0] ); // Assume CIELab color space, where first coordinate is luminance/gray-scale
      Dest.resize ( numBins, gray.size ( ) );
      _scalarFilterResponceLocalHistogramOp.apply ( gray, Dest );
    } else {
      Dest.resize ( 3 * numBins, Arg[0].size ( ) );
      aol::MultiVector<RealType> tmp ( numBins, Arg[0].size ( ) );
      for ( int c=0; c<3 ; ++ c ) {
        _scalarFilterResponceLocalHistogramOp.apply ( Arg[c], tmp );
        for ( int i=0; i<numBins ; ++i ) Dest[c * numBins + i] = tmp[i];
      }
    }
  }
  
  virtual void applyAdd ( const PictureType &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};

/**
 * \brief Base class for conversion of a 2D (scalar or multi-dimensional) image into a MultiVector of multiple histograms of windows filtered with different kernels (spectral histograms).
 * \author mevenkamp
 * \ingroup Segmentation
 */
template <typename _RealType, typename _PictureType>
class BaseLocalSpectralHistogramOp : public aol::Op<_PictureType, aol::MultiVector<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  std::vector<FilterResponseLocalHistogramOp<_RealType, _PictureType>*> _histOps;
public:
  BaseLocalSpectralHistogramOp ( ) { }
  
  void apply ( const PictureType &Arg, aol::MultiVector<RealType> &Dest ) const {
    for ( int i=0; i<static_cast<int>( _histOps.size ( ) ) ; ++i ) {
      aol::MultiVector<RealType> *dest = new aol::MultiVector<RealType> ( );
      _histOps[i]->apply ( Arg, *dest );
      Dest.appendReference ( *dest, true );
    }
  }
  
  virtual void applyAdd ( const PictureType &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void pushBack ( const FilterConvolutionOp<PictureType> *Filter, const int IntegrationScale, const int NumBins ) {
    _histOps.push_back ( new FilterResponseLocalHistogramOp<RealType, PictureType> ( Filter, IntegrationScale, NumBins ) );
  }
};

/**
 * \brief Local spectral histogram operator for 2D scalar images.
 * \author mevenkamp
 * \ingroup Segmentation
 */
template <typename _RealType, typename _PictureType>
class LocalSpectralHistogramOp : public BaseLocalSpectralHistogramOp<_RealType, _PictureType> { };

/**
 * \brief Local spectral histogram operator for 2D multi-dimensional images (e.g. color)
 * \author mevenkamp
 * \ingroup Segmentation
 */
template <typename _RealType>
class LocalSpectralHistogramOp<_RealType, qc::MultiArray<_RealType, qc::QC_2D, 3> > : public BaseLocalSpectralHistogramOp<_RealType, qc::MultiArray<_RealType, qc::QC_2D, 3> > {
public:
  void pushBack ( const FilterConvolutionOp<qc::ScalarArray<_RealType, qc::QC_2D> > *Filter, const int IntegrationScale, const int NumBins,
                  const bool ConvertToGray = true ) {
    BaseLocalSpectralHistogramOp<_RealType, qc::MultiArray<_RealType, qc::QC_2D, 3> >::_histOps.push_back ( new FilterResponseLocalHistogramOp<_RealType, qc::MultiArray<_RealType, qc::QC_2D, 3> > ( Filter, IntegrationScale, NumBins, ConvertToGray ) );
  }
};

/**
 * \author mevenkamp
 * \ingroup Segmentation
 */
template <typename RealType>
static void histogramEqualization ( qc::ScalarArray<RealType, qc::QC_2D> &ArgDest ) {
  const int argMax = ArgDest.getMaxValue ( );
  aol::Vector<int> frequencies ( argMax + 1 ), lookUp ( argMax + 1 );
  for ( int k=0; k<ArgDest.size ( ) ; ++k ) ++frequencies[static_cast<int> ( ArgDest[k] )];
  for ( int i=0; i<=argMax ; ++i ) {
    for ( int j=0; j<=i ; ++j ) lookUp[i] += frequencies[j];
    lookUp[i] *= argMax / static_cast<RealType> ( ArgDest.size ( ) );
  }
  for ( int k=0; k<ArgDest.size ( ) ; ++k ) ArgDest[k] = lookUp[static_cast<int> ( ArgDest[k] )];
}

/**
 * \brief Implements the intensity filter kernel (identity matrix)
 * \author mevenkamp
 * \ingroup Segmentation
 */
template <typename _PictureType>
class IntensityFilterOp : public FilterConvolutionOp<_PictureType> {
  typedef _PictureType PictureType;
  typedef typename PictureType::RealType RealType;
public:
  IntensityFilterOp ( ) : FilterConvolutionOp<PictureType> ( 1 ) {
    this->_M.set ( 0, 0, 1.0 );
  }
};

/**
 * \brief Implements the Gabor filter kernel: \f$g(x,y;\sigma,\theta,\psi,\gamma,\lambda) = \frac{1}{2\pi \sigma^2} \exp{-\frac{{x^\prime}^2 + \gamma^2 {y^\prime}^2}{2 \sigma^2}} \cos( 2 \pi \lambda x^\prime + \psi ), x^\prime = x \cos \theta + y \sin \theta, y^\prime = -x \sin\theta + y \cos\theta, \sigma = 0.25\text{Size}, \theta = \text{Theta}, \psi = 0, \gamma = 1, \lambda = (2 \sigma)^{-1}\f$
 * \author mevenkamp
 * \ingroup Segmentation
 */
template <typename _PictureType>
class GaborFilterOp : public FilterConvolutionOp<_PictureType> {
  typedef _PictureType PictureType;
  typedef typename PictureType::RealType RealType;
protected:
  const RealType _s;
  const RealType _theta;
  const RealType _phase;
  const RealType _gamma;
  const RealType _f;
public:
  GaborFilterOp ( const int Size, const RealType Theta )
    : FilterConvolutionOp<PictureType> ( Size ),
      _s ( 0.5 * this->_offset ), _theta ( Theta ), _phase ( 0.0 ), _gamma ( 1.0 ),
      _f ( 1.0 / ( 2.0 * _s ) ) {
    RealType sSqr = aol::Sqr<RealType> ( _s );
    RealType xPrime, yPrime;
    for ( int y=-this->_offset; y<=this->_offset ; ++y ) {
      for ( int x=-this->_offset; x<=this->_offset ; ++x ) {
        xPrime = x * cos(_theta) + y * sin(_theta);
        yPrime = y * cos(_theta) - x * sin(_theta);
        this->_M.set ( x+this->_offset, y+this->_offset, 1.0 / ( 2.0 * aol::NumberTrait<RealType>::pi * sSqr ) * exp ( -0.5 * ( aol::Sqr<RealType> ( xPrime )
                                                         + aol::Sqr<RealType> ( yPrime * _gamma ) ) / sSqr ) * cos ( 2.0 * aol::NumberTrait<RealType>::pi * _f * xPrime + _phase ) );
      }
    }
  }
};

/**
 * \brief Implements the Gaussian filter kernel: \f$g(x,y;\sigma) = \frac{1}{2\pi \sigma^2} \exp{-\frac{x^2 + y^2}{2 \sigma^2}}\f$
 * \author mevenkamp
 * \ingroup Segmentation
 */
template <typename _PictureType>
class GaussianFilterOp : public FilterConvolutionOp<_PictureType> {
  typedef _PictureType PictureType;
  typedef typename PictureType::RealType RealType;
protected:
  const RealType _sigma, _minusInvTwoSigmaSqr;
public:
  GaussianFilterOp ( const int Size, const RealType Sigma )
    : FilterConvolutionOp<PictureType> ( Size ),
      _sigma ( Sigma ), _minusInvTwoSigmaSqr ( -1.0 / ( 2.0 * aol::Sqr<RealType> ( Sigma ) ) ) {
    RealType sum = 0;
    for ( int y=-this->_offset; y<=this->_offset ; ++y ) {
      for ( int x=-this->_offset; x<=this->_offset ; ++x ) {
        this->_M.set ( x+this->_offset, y+this->_offset, exp ( ( aol::Sqr<RealType> ( x ) + aol::Sqr<RealType> ( y ) ) * _minusInvTwoSigmaSqr ) );
        sum += this->_M.get ( x+this->_offset, y+this->_offset );
      }
    }
    this->_M *= 1.0 / sum;
  }
  
  RealType getMinusInvTwoSigmaSqr ( ) const {
    return _minusInvTwoSigmaSqr;
  }
};

/**
 * \brief Implements the Laplacian-of-Gaussian (LoG) filter kernel
 * \author mevenkamp
 * \ingroup Segmentation
 */
template <typename _PictureType>
class LaplacianOfGaussianFilterOp : public FilterConvolutionOp<_PictureType> {
  typedef _PictureType PictureType;
  typedef typename PictureType::RealType RealType;
protected:
  const RealType _sigma, _sigmaSqr;
public:
  LaplacianOfGaussianFilterOp ( const int Size, const RealType Sigma )
    : FilterConvolutionOp<PictureType> ( Size ),
      _sigma ( Sigma ), _sigmaSqr ( aol::Sqr<RealType> ( Sigma ) ) {
    GaussianFilterOp<PictureType> gaussianFilterOp ( this->_size, _sigma );
    for ( int y=-this->_offset; y<=this->_offset ; ++y )
      for ( int x=-this->_offset; x<=this->_offset ; ++x )
        this->_M.set ( x+this->_offset, y+this->_offset, gaussianFilterOp.getFilterKernel ( ).get ( x+this->_offset, y+this->_offset )
                                                         * ( 1 + ( aol::Sqr<RealType> ( x ) + aol::Sqr<RealType> ( y ) ) * gaussianFilterOp.getMinusInvTwoSigmaSqr ( ) ) );
    this->_M *= -1.0 / ( aol::NumberTrait<RealType>::pi * aol::Sqr<RealType> ( _sigmaSqr ) );
  }
};

/*
 * END: Local spectral histogram ops
 */


template <typename RealType>
static void saveOptimalRegionIndicators ( const char* baseName, const qc::ScalarArray<int, qc::QC_2D> &GroundTruth, const aol::MultiVector<RealType> &Indicator ) {
  qc::ScalarArray<int, qc::QC_2D> groundTruth ( GroundTruth, aol::DEEP_COPY );
  transformGrayValuesToLabels ( groundTruth );
  for ( int l=0; l<=groundTruth.getMaxValue ( ) ; ++l ) {
    int n = 0;
    aol::Vector<RealType> mean ( Indicator.numComponents ( ) );
    for ( int i=0; i<groundTruth.size ( ) ; ++i ) {
      if ( groundTruth[i] == l ) {
        for ( int c=0; c<Indicator.numComponents ( ) ; ++c )
          mean[c] += Indicator[c][i];
        ++n;
      }
    }
    mean /= static_cast<RealType> ( n );

    qc::ScalarArray<RealType, qc::QC_2D> errorImg ( groundTruth.getNumX ( ), groundTruth.getNumY ( ) );
    aol::Vector<RealType> error ( Indicator.numComponents ( ) );
    for ( int i=0; i<groundTruth.size ( ) ; ++i ) {
      error = mean;
      for ( int c=0; c<Indicator.numComponents ( ) ; ++c )
        error[c] -= Indicator[c][i];
      errorImg[i] = error.norm ( );
    }
    errorImg.setOverflowHandlingToCurrentValueRange ( );
    errorImg.savePNG ( aol::strprintf ( "%s_%d.png", baseName, l ).c_str ( ) );
  }
}


/**
 * A class for PCA-based segmentation of high-dimensional piece-wise constant vector fields
 * e.g. arising in texture segmentation based on local spectral histograms or unitary transforms (e.g. Stockwell)
 *
 * As BaseClass, use either PiecewiseConstantMultiPhaseMSSegmentor<ConfiguratorType, FirstOrderPrimalDualMultiPhaseMSSegmentor<ConfiguratorType> > (multi phase)
 *                       or PiecewiseConstantTwoPhaseMSSegmentor<ConfiguratorType, FirstOrderPrimalDualTwoPhaseMSSegmentor<ConfiguratorType> > (two phase)
 *
 * \author Mevenkamp
 * \ingroup Segmentation
 */
template <typename ConfiguratorType, typename BaseClass>
class PCAMSSegmentor : public BaseClass {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
protected:
  const aol::MultiVector<RealType> _highDimImageMVec;
  aol::MultiVector<RealType> _eigenVecs;
  aol::Vector<RealType> _eigenVals;
  const int _numEvals;
  const RealType _omega;
  const bool _orthonormalizeMeanValues;
  std::string _outputDir;
  qc::ScalarArray<int, qc::QC_2D> _groundTruth;
public:
  PCAMSSegmentor ( const typename ConfiguratorType::InitType &Initializer,
                   const RealType Gamma,
                   aol::MultiVector<RealType> &ImageMVec,
                   const bool InitializeGrayValues = true,
                   const int NumSegments = 0,
                   const int NumEvals = 0,
                   const RealType Omega = 0.05,
                   const bool OrthonormalizeMeanValues = false,
                   const int NumGhostCells = 0,
                   const std::string &OutputDir = "",
                   const bool Quiet = false )
    : BaseClass ( Initializer, Gamma, ImageMVec, false, Quiet, OutputDir ),
      _highDimImageMVec ( ImageMVec ), _numEvals ( ( NumEvals == 0 ) ? NumSegments : NumEvals ), _omega ( Omega ),
      _orthonormalizeMeanValues ( OrthonormalizeMeanValues ), _outputDir ( OutputDir ) {
    this->setNumGhostCells ( NumGhostCells );
    setEigenVecs ( );
    this->setNumSegments ( ( NumSegments > 0 ) ? NumSegments : _eigenVecs.numComponents ( ) );
    if ( NumSegments == 0 && !this->_quietMode ) std::cerr << "PCA detected " << this->_numSegments << " segments." << std::endl;
    setProjectedCoefficients ( );
    this->_meanValues.resize ( this->_numSegments, this->_imageMVec.numComponents ( ) );
    if ( InitializeGrayValues ) initializeGrayValues ( );
        
    if ( _outputDir != "" && _groundTruth.size ( ) > 0 )
      saveOptimalRegionIndicators ( aol::strprintf ( "%s/pcaRegionIndicator", _outputDir.c_str ( ) ).c_str ( ), _groundTruth, this->_imageMVec );
  }
  
  static void getEigenVecs ( aol::MultiVector<RealType> &EigenVecs, aol::Vector<RealType> &EigenVals, const aol::MultiVector<RealType> &HighDimInput, const int NumEvals, const RealType Omega = 0,
                             const bool Quiet = false ) {
    getPCAEigenVecs ( EigenVecs, EigenVals, HighDimInput, NumEvals, Omega, Quiet, false, true );
  }
  
  static void getProjectedCoefficients ( aol::MultiVector<RealType> &LowDimOutput, const aol::MultiVector<RealType> &HighDimInput, const aol::MultiVector<RealType> &EigenVecs,
                                         const bool Quiet = false ) {
    getPCAProjectedCoefficients ( LowDimOutput, HighDimInput, EigenVecs, Quiet );
  }
protected:
  void setEigenVecs ( ) {
    getEigenVecs ( _eigenVecs, _eigenVals, _highDimImageMVec, _numEvals, _omega, this->_quietMode );
  }
  
  void setProjectedCoefficients ( ) const {
    getProjectedCoefficients ( this->_imageMVec, _highDimImageMVec, _eigenVecs, this->_quietMode );
  }
  
  void orthoNormalizeGrayValues ( ) {
    if ( !this->_quietMode ) std::cerr << "Orthonormalizing basis representation of mean values.." << std::endl;
    // Assemble matrix A containing the coefficient representation of mean values in eigen vector basis
    aol::FullMatrix<RealType> A ( _eigenVecs.numComponents ( ), this->_numSegments );
    for ( int l=0; l<this->_numSegments ; ++l )
      for ( int ev=0; ev<_eigenVecs.numComponents ( ) ; ++ev )
        A.set ( ev, l, this->_meanValues[l][ev] );
    
    // Perform orthonormalization of mean values
    aol::FullMatrix<RealType> Q ( _eigenVecs.numComponents ( ), this->_numSegments ), R ( this->_numSegments, this->_numSegments );
    aol::QRDecomposeModifiedGramSchmidt<RealType> qrGramSchmidt;
    qrGramSchmidt.transform ( A, R, Q );
    for ( int l=0; l<this->_numSegments ; ++l )
      for ( int ev=0; ev<_eigenVecs.numComponents ( ) ; ++ev )
        this->_meanValues[l][ev] = Q.get ( ev, l );
    
    // Apply change of basis to eigen vectors
    aol::FullMatrix<RealType> QT ( Q.getNumCols ( ), Q.getNumRows ( ) );
    Q.transposeTo ( QT );
    aol::FullMatrix<RealType> AQT ( A.getNumRows ( ), A.getNumRows ( ) );
    AQT.makeProduct ( A, QT );
    aol::QRInverse<RealType> qrInv ( AQT );
    aol::FullMatrix<RealType> B = qrInv.getFM ( );
    aol::Vector<RealType> vec ( _eigenVecs.numComponents ( ) ), transformedVec ( _eigenVecs.numComponents ( ) );
    for ( int k=0; k<_eigenVecs[0].size ( ) ; ++k ) {
      _eigenVecs.getTo ( k, vec );
      B.apply ( vec, transformedVec );
      _eigenVecs.set ( k, transformedVec );
    }
    
    getProjectedCoefficients ( this->_imageMVec, _highDimImageMVec, _eigenVecs );
  }

  virtual void initializeGrayValues ( ) {
    aol::KMeansClusterer<RealType> kMeansClusterer;
    aol::MultiVector<RealType> clusters;
    kMeansClusterer.applyMultipleRNG ( this->_imageMVec, clusters, this->_numSegments );
    const int imageDim = this->_imageMVec.numComponents ( );
    for ( int l=0; l<this->_numSegments ; ++l )
      for ( int j=0; j<imageDim ; ++j )
        this->_meanValues[l][j] = clusters[j][l];
    if ( !this->_quietMode ) cerr << this->_meanValues << endl;
    
    if ( _orthonormalizeMeanValues ) orthoNormalizeGrayValues ( );
  }
  virtual void updateGrayValues ( const typename BaseClass::ArrayType &CurrentSegmentation ) {
    BaseClass::updateGrayValues ( CurrentSegmentation );
    if ( _orthonormalizeMeanValues ) orthoNormalizeGrayValues ( );
  }
public:
  void setMeanValues ( const aol::Vector<int> &Indices ) {
    BaseClass::setMeanValues ( Indices );
    if ( _orthonormalizeMeanValues ) orthoNormalizeGrayValues ( );
  }
  
  void setGroundTruth ( const qc::ScalarArray<int, qc::QC_2D> &GroundTruth ) {
    _groundTruth.reallocate ( GroundTruth.getSize ( ) );
    _groundTruth = GroundTruth;
  }
};

enum LOCALFEATUREMSS_METHOD {
  PCA_PIECEWISECONSTANT,
  SUBSPACE
};
  
enum UNKNOWNREGION {
  UNKNOWNREGION_NONE,
  UNKNOWNREGION_DURINGSEG,
  UNKNOWNREGION_POSTSEG
};

template<typename _RealType, typename _PictureType>
class LocalFeatureMSSegmentationOptions {
public:
  // I/O
  const _PictureType &input;
  std::string outputDir;
  bool quietMode;
  aol::ProgressBar<> *progressBar;
  qc::ScalarArray<int, qc::QC_2D> groundTruth;
  
  LOCALFEATUREMSS_METHOD method;
  
  // General parameters
  short resampleFactor;
  int numSegments;
  int numGhostCells;
  bool includeBoundary;
  bool forceMultiPhaseSegmentation;
  UNKNOWNREGION unknownRegion;
  _RealType lowEdgenessThreshold;
  aol::MultiVector<_RealType> seedPositions;
  
  // Parameters for PCA
  int numEvals;
  _RealType omega;
  bool orthonormalizeMeanValues;

  // Parameters for Clustering
  int numDominantOperators;
  int transformScale;
  
  bool forceGroundTruthClusters; // for testing purposes
  
  // Parameters for Mumford-Shah
  _RealType gamma;
  int maxIt;
  _RealType epsilon;
  int numOuterIterations;
  
  // Operators
  std::vector<aol::Op<_PictureType, aol::MultiVector<_RealType> >*> operators;
  aol::Vector<_RealType> operatorWeights;
  
  LocalFeatureMSSegmentationOptions ( const _PictureType &Input, const std::string OutputDir = "", const bool Quiet = true )
    : input ( Input ), outputDir ( OutputDir ), quietMode ( Quiet ), progressBar ( NULL ),
      method ( PCA_PIECEWISECONSTANT ),
      resampleFactor ( 1 ), numSegments ( 1 ), numGhostCells ( 0 ),
      includeBoundary ( true ), forceMultiPhaseSegmentation ( false ), unknownRegion ( UNKNOWNREGION_NONE ),
      lowEdgenessThreshold ( 0.5 ), seedPositions ( ),
      numEvals ( 0 ), omega ( 0.05 ), orthonormalizeMeanValues ( false ),
      numDominantOperators ( 1 ), transformScale ( 21 ),
      forceGroundTruthClusters ( false ),
      gamma ( gammaDefault ), maxIt ( maxItDefault ), epsilon ( epsilonDefault ), numOuterIterations ( 3 ) { }
  
  LocalFeatureMSSegmentationOptions ( const LocalFeatureMSSegmentationOptions<_RealType, _PictureType> &Options )
    : input ( Options.input ), outputDir ( Options.outputDir ), quietMode ( Options.quietMode ), progressBar ( Options.progressBar ), groundTruth ( Options.groundTruth ),
      method ( Options.method ),
      resampleFactor ( Options.resampleFactor ), numSegments ( Options.numSegments ), numGhostCells ( Options.numGhostCells ),
      includeBoundary ( Options.includeBoundary ), forceMultiPhaseSegmentation ( Options.forceMultiPhaseSegmentation ), unknownRegion ( Options.unknownRegion ),
      lowEdgenessThreshold ( Options.lowEdgenessThreshold ), seedPositions ( Options.seedPositions ),
      numEvals ( Options.numEvals ), omega ( Options.omega ), orthonormalizeMeanValues ( Options.orthonormalizeMeanValues ),
      numDominantOperators ( Options.numDominantOperators ), transformScale ( Options.transformScale ),
      forceGroundTruthClusters ( Options.forceGroundTruthClusters ),
      gamma ( Options.gamma ), maxIt ( Options.maxIt ), epsilon ( Options.epsilon ), numOuterIterations ( Options.numOuterIterations ) { }
  
  static _RealType gammaDefault, epsilonDefault;
  static int maxItDefault;
};
template<typename _RealType, typename _PictureType> _RealType LocalFeatureMSSegmentationOptions<_RealType, _PictureType>::gammaDefault = 1e-4;
template<typename _RealType, typename _PictureType> _RealType LocalFeatureMSSegmentationOptions<_RealType, _PictureType>::epsilonDefault = 1e-3;
template<typename _RealType, typename _PictureType> int LocalFeatureMSSegmentationOptions<_RealType, _PictureType>::maxItDefault = 10000;

  
/**
 * A class for Mumford-Shah and feature based segmentation 2D images. Utilizes PCA and clustering to reduce the dimension and initialize the segmentation.
 *
 * \author Mevenkamp
 * \ingroup Segmentation
 */
template <typename _RealType, typename _PictureType>
class LocalFeatureMSSegmentor {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> SegmentationType;
  typedef qc::MultiArray<RealType, qc::QC_2D, 3> ColoredPictureType;
  typedef qc::RectangularGridConfigurator<RealType, qc::QC_2D, aol::GaussQuadrature<RealType, qc::QC_2D, 3> > ConfiguratorType;
  typedef typename ConfiguratorType::InitType InitType;
  typedef typename qc::ComponentsCollection<RealType>::NonEmptyComponentsIterator NonEmptyComponentsIterator;
  typedef LocalFeatureMSSegmentationOptions<RealType, PictureType> OptionsType;
  typedef PiecewiseConstantTwoPhaseMSSegmentor<ConfiguratorType, 0> TwoPhaseSegmentorType;
  typedef PiecewiseConstantMultiPhaseMSSegmentor<ConfiguratorType> MultiPhaseSegmentorType;
protected:
  std::string _outputDir;
  aol::ProgressBar<> *_progressBar;
  bool _quietMode;
  bool _catchCtrlC;
  mutable sigfunc _previousCtrlCHandler;
  int _numDominantComponents;
public:
  LocalFeatureMSSegmentor ( ) : _outputDir ( "" ), _quietMode ( false ), _catchCtrlC ( false ) { }
  
  virtual ~LocalFeatureMSSegmentor ( ) { }
  
  void apply ( OptionsType &Options, qc::ScalarArray<int, qc::QC_2D> &Segmentation ) {
    Segmentation.reallocate ( Options.input.getNumX ( ), Options.input.getNumY ( ) );
    initialize ( Options );
    getSegmentationFromLocalFeatures ( Options, Segmentation );
  }
  
  void getLocalFeatures ( const OptionsType &Options, aol::MultiVector<RealType> &LocalFeatures ) {
    InitType *grid = NULL;
    getLocalFeatures ( Options, LocalFeatures, grid );
  }
  
  void setQuietMode ( bool Quiet ) {
    _quietMode = Quiet;
  }
  
  void setCatchCtrlC ( bool catchCtrlC ) {
    _catchCtrlC = catchCtrlC;
  }
protected:
  void initialize ( const OptionsType &Options ) {
    _outputDir = Options.outputDir;
    _quietMode = Options.quietMode;
    _progressBar = Options.progressBar;
    
    // Copy input to destination
    if ( _outputDir != "" ) {
      PictureType u ( Options.input );
      u.scaleValuesTo01 ( );
      u.setOverflowHandling ( aol::SCALE, 0, 1 );
      u.savePNG ( aol::strprintf ( "%s/input.png", _outputDir.c_str ( ) ).c_str ( ) );
    }
  }
  
  void getSegmentationFromLocalFeatures ( OptionsType &Options, qc::ScalarArray<int, qc::QC_2D> &Dest ) {
    aol::MultiVector<RealType> localFeatures;
    InitType *grid = NULL;
    getLocalFeatures ( Options, localFeatures, grid );
    if ( Options.method == PCA_PIECEWISECONSTANT ) getHardSegmentationPCAPiecewiseConstant ( Options, Dest, *grid, localFeatures );
    else if ( Options.method == SUBSPACE ) getHardSegmentationSubSpace ( Options, Dest, *grid, localFeatures );
    if ( !Options.includeBoundary ) extendSegmentationToBoundary ( Options, Dest );
    upsampleHardSegmentation ( Dest, aol::Vec2<short> ( Options.input.getNumX ( ), Options.input.getNumY ( ) ) );
    
    if ( _outputDir != "" )
      plotSegmentationBoundariesOntoImage<RealType, PictureType> ( aol::strprintf ( "%s/boundaries", _outputDir.c_str ( ) ).c_str ( ), Options.input, Dest );
  }
  
  void getLocalFeatures ( const OptionsType &Options, aol::MultiVector<RealType> &LocalFeatures, InitType *&Grid ) {
    // Resample image
    PictureType u ( Options.input.getNumX ( ) / Options.resampleFactor, Options.input.getNumY ( ) / Options.resampleFactor );
    u.resampleFrom ( Options.input );
    
    // Post-process operator weights
    aol::Vector<RealType> operatorWeights ( static_cast<int>( Options.operators.size ( ) ) );
    RealType sum = Options.operatorWeights.sum ( ) + Options.operators.size ( ) - Options.operatorWeights.size ( );
    operatorWeights.setAll ( 1.0 / sum );
    for ( int k=0; k<Options.operatorWeights.size ( ) ; ++k ) operatorWeights[k] = Options.operatorWeights[k] / sum;
    
    // Apply specified operators and assemble MultiVector of local features
    qc::FastILexMapper<qc::QC_2D> mapper ( u.getNumX ( ), u.getNumY ( ) );
    int offset = ( Options.includeBoundary ) ? 0 : ( Options.transformScale - 1 ) / 2;
    Grid = new InitType ( qc::GridSize<qc::QC_2D> ( u.getNumX ( ) - 2 * offset, u.getNumY ( ) - 2 * offset ) );
    aol::MultiVector<RealType> localFeatures;
    if ( this->_progressBar != NULL ) {
      this->_progressBar->setText ( "Computing local descriptors" );
      this->_progressBar->start ( Options.operators.size ( ) );
    }
    _numDominantComponents = 0;
    for ( int k=0; k<static_cast<int>(Options.operators.size ( )) ; ++k ) {
      Options.operators[k]->apply ( u, localFeatures );
      localFeatures *= operatorWeights[k];
      if ( k < Options.numDominantOperators ) _numDominantComponents += localFeatures.numComponents ( );
      int c0 = LocalFeatures.numComponents ( );
      LocalFeatures.resize ( c0 + localFeatures.numComponents ( ), Grid->getNumberOfNodes ( ) );
      int j = 0;
      for ( int y=offset; y<u.getNumY ( )-offset ; ++y ) {
        for ( int x=offset; x<u.getNumX ( )-offset ; ++x ) {
          for ( int c=0; c<localFeatures.numComponents ( ) ; ++c )
            LocalFeatures[c0+c][j] = localFeatures[c][mapper.getGlobalIndex ( x, y )];
          ++j;
        }
      }
      if ( this->_progressBar != NULL ) (*this->_progressBar)++;
    }
    if ( this->_progressBar != NULL ) this->_progressBar->finish ( );
    
    if ( _outputDir != "" ) {
      u.scaleValuesTo01 ( );
      u.setOverflowHandling ( aol::SCALE, 0, 1 );
      u.savePNG ( aol::strprintf ( "%s/input_resampled.png", _outputDir.c_str ( ) ).c_str ( ) );
    }
    
    if ( _outputDir != "" && Options.groundTruth.size ( ) > 0 )
      saveOptimalRegionIndicators ( aol::strprintf ( "%s/highDimFeatureIndicator", _outputDir.c_str ( ) ).c_str ( ), Options.groundTruth, LocalFeatures );
  }
  
  void getHardSegmentationPCAPiecewiseConstant ( OptionsType &Options, qc::ScalarArray<int, qc::QC_2D> &HardSegmentation, const InitType &Grid, aol::MultiVector<RealType> &Array ) const {
    const bool initialMeanPositions = Options.numSegments > 0 && Options.seedPositions.numComponents ( ) == Options.numSegments && Options.seedPositions[0].size ( ) == 2;
    aol::Vector<int> indices ( Options.numSegments );
    int offset = ( Options.includeBoundary ) ? 0 : ( Options.transformScale - 1 ) / 2;
    if ( initialMeanPositions ) {
      qc::FastILexMapper<qc::QC_2D> mapper ( Grid );
      int x, y;
      for ( int l=0; l<Options.numSegments ; ++l ) {
        x = Options.seedPositions[l][0] / Options.resampleFactor - offset;
        y = Options.seedPositions[l][1] / Options.resampleFactor - offset;
        if ( x < 0 || y < 0 ) throw aol::Exception ( "Seed position outside of support! If boundary is not included, choose seed positions away from boundary.", __FILE__, __LINE__ );
        indices[l] = mapper.getGlobalIndex ( x, y );
      }
    } else if ( Options.seedPositions.numComponents ( ) > 0 ) std::cerr << "WARNING: seed positions were specified but have wrong dimensions! Will use automatic seeding instead.." << std::endl;
    
    HardSegmentation.reallocate ( Grid );
    if ( !Options.forceMultiPhaseSegmentation && ( Options.numSegments == 2 || ( Options.numSegments == 1 && Options.unknownRegion == UNKNOWNREGION_DURINGSEG ) ) ) {
      // Compute soft segmentation
      SegmentationType segmentation ( Grid );
      PCAMSSegmentor<ConfiguratorType, TwoPhaseSegmentorType> segmentor ( Grid, Options.gamma, Array, !initialMeanPositions, Options.numSegments,
                                                                          Options.numEvals, Options.omega, Options.orthonormalizeMeanValues, Options.numGhostCells, _outputDir, _quietMode );
      segmentor.setCatchCtrlC ( _catchCtrlC );
      segmentor.setMaxIterations ( Options.maxIt );
      segmentor.setStopEpsilon ( Options.epsilon );
      segmentor.setOuterIterations ( Options.numOuterIterations );
      if ( initialMeanPositions ) segmentor.setMeanValues ( indices );
      if ( Options.unknownRegion == UNKNOWNREGION_DURINGSEG ) segmentor.setUnknownRegion ( true );
      segmentor.segmentAndAdjustGrayValues ( segmentation );
      
      if (_outputDir != "" ) segmentation.save ( aol::strprintf ( "%s/softSeg%s", _outputDir.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
      
      segmentor.getHardSegmentation ( HardSegmentation, segmentation );
      
      if ( Options.unknownRegion == UNKNOWNREGION_POSTSEG ) segmentor.addUnknownRegion ( HardSegmentation );
    } else {
      // Compute initial segmentation
      aol::VectorContainer<SegmentationType> initialSegmentation ( aol::Max<int> ( 1, Options.numSegments ), SegmentationType ( Grid ) );
      if ( !initialMeanPositions ) getInitialSegmentationFromLowEdgenessDominantFeatureClustering ( Options, initialSegmentation, Grid, Array );
      
      // Compute soft segmentation
      aol::VectorContainer<SegmentationType> segmentation ( Options.numSegments, SegmentationType ( Grid ) );
      PCAMSSegmentor<ConfiguratorType, MultiPhaseSegmentorType> segmentor ( Grid, Options.gamma, Array, false, Options.numSegments,
                                                                            Options.numEvals, Options.omega, Options.orthonormalizeMeanValues, Options.numGhostCells, _outputDir, _quietMode );
      segmentor.setOutputDirectory ( _outputDir );
      segmentor.setCatchCtrlC ( _catchCtrlC );
      segmentor.setMaxIterations ( Options.maxIt );
      segmentor.setStopEpsilon ( Options.epsilon );
      segmentor.setOuterIterations ( Options.numOuterIterations );
      if ( initialMeanPositions ) segmentor.setMeanValues ( indices );
      else segmentor.setMeanValuesFromInitialSegmentation ( initialSegmentation );
      if ( Options.unknownRegion == UNKNOWNREGION_DURINGSEG ) segmentor.setUnknownRegion ( true );
      segmentor.segmentAndAdjustGrayValues ( segmentation );
      
      if ( _outputDir != "" ) {
        for ( int l=0; l<segmentation.size ( ) ; ++l )
          segmentation[l].save ( aol::strprintf ( "%s/softSeg%d%s", _outputDir.c_str ( ), l, qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
      }
      
      segmentor.getHardSegmentation ( HardSegmentation, segmentation );
      
      if ( Options.unknownRegion == UNKNOWNREGION_POSTSEG ) segmentor.addUnknownRegion ( HardSegmentation );
    }
  }
  
  void getHardSegmentationSubSpace ( OptionsType &Options, qc::ScalarArray<int, qc::QC_2D> &HardSegmentation, const InitType &Grid, aol::MultiVector<RealType> &Array ) const {
    // Compute initial segmentation
    aol::VectorContainer<SegmentationType> initialSegmentation ( aol::Max<int> ( 1, Options.numSegments ), SegmentationType ( Grid ) );
    // TODO
    
    // Compute soft segmentation
    aol::VectorContainer<SegmentationType> segmentation ( Options.numSegments, SegmentationType ( Grid ) );
    PiecewiseSubSpaceMultiPhaseMSSegmentor<ConfiguratorType, FirstOrderPrimalDualMultiPhaseMSSegmentor<ConfiguratorType> > segmentor ( Grid, Options.gamma, Array,
                                                                                                                                           Options.numSegments, Options.numEvals, Options.numGhostCells );
    segmentor.setOutputDirectory ( _outputDir );
    segmentor.setCatchCtrlC ( _catchCtrlC );
    segmentor.setMaxIterations ( Options.maxIt );
    segmentor.setStopEpsilon ( Options.epsilon );
    segmentor.setOuterIterations ( Options.numOuterIterations );
    segmentor.setSubSpacesFromInitialSegmentation ( initialSegmentation );
    segmentor.segmentAndAdjustSubSpaces ( segmentation );
    
    if ( _outputDir != "" ) {
      for ( int l=0; l<Options.numSegments ; ++l )
        segmentation[l].save ( aol::strprintf ( "%s/softSeg%d%s", _outputDir.c_str ( ), l, qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    }
    
    segmentor.getHardSegmentation ( HardSegmentation, segmentation );
  }

  void extendSegmentationToBoundary ( const OptionsType &Options, qc::ScalarArray<int, qc::QC_2D> &HardSegmentation ) const {
    int offset = ( Options.includeBoundary ) ? 0 : ( Options.transformScale - 1 ) / 2;
    qc::ScalarArray<int, qc::QC_2D> hardSegmentation ( HardSegmentation );
    HardSegmentation.reallocate ( Options.input.getNumX ( ) / Options.resampleFactor, Options.input.getNumY ( ) / Options.resampleFactor );
    for ( int y=0; y<HardSegmentation.getNumY ( ) ; ++y )
      for ( int x=0; x<HardSegmentation.getNumX ( ) ; ++x )
        HardSegmentation.set ( x, y, hardSegmentation.getClip ( x - offset, y - offset ) );
  }
  
  void upsampleHardSegmentation ( qc::ScalarArray<int, qc::QC_2D> &HardSegmentation, const aol::Vec2<short> &Size ) const {
    qc::ScalarArray<int, qc::QC_2D> hardSegmentation ( HardSegmentation );
    HardSegmentation.reallocate ( Size[0], Size[1] );
    for ( int y=0; y<HardSegmentation.getNumY ( ) ; ++y )
      for ( int x=0; x<HardSegmentation.getNumX ( ) ; ++x )
        HardSegmentation.set ( x, y, hardSegmentation.get ( round ( x * hardSegmentation.getNumX ( ) / Size[0] ), round ( y * hardSegmentation.getNumY ( ) / Size[1] ) ) );
  }

  void getInitialSegmentationFromLowEdgenessDominantFeatureClustering ( OptionsType &Options, aol::VectorContainer<SegmentationType> &Segmentation, const InitType &Grid, aol::MultiVector<RealType> &Array ) const {
    // Reduce array to the components created by the dominant operators
    // Then perform PCA on it
    aol::MultiVector<RealType> reducedArray ( _numDominantComponents, Array[0].size ( ) );
    for ( int i=0; i<Array[0].size ( ) ; ++i )
      for ( int c=0; c<reducedArray.numComponents ( ); ++c )
        reducedArray[c][i] = Array[c][i];
    aol::MultiVector<RealType> eigenVecs, reducedPCAArray;
    aol::Vector<RealType> eigenVals;
    PCAMSSegmentor<ConfiguratorType, MultiPhaseSegmentorType>::getEigenVecs ( eigenVecs, eigenVals, reducedArray, Options.numEvals, Options.omega );
    PCAMSSegmentor<ConfiguratorType, MultiPhaseSegmentorType>::getProjectedCoefficients ( reducedPCAArray, reducedArray, eigenVecs );
    if ( Options.numSegments == 0 ) {
      Options.numSegments = eigenVecs.numComponents ( );
      if ( !this->_quietMode ) std::cerr << "PCA detected " << Options.numSegments << " segments." << std::endl;
    }
    
    if ( _outputDir != "" ) {
      eigenVals.saveASCII ( aol::strprintf ( "%s/eigenValues_full.txt", this->_outputDir.c_str ( ) ).c_str ( ) );
      std::ofstream txtFile ( aol::strprintf ( "%s/numSamples_full.txt", this->_outputDir.c_str ( ) ).c_str ( ) );
      txtFile << reducedArray[0].size ( ) << std::endl;
      txtFile.close ( );
    }
    
    // Determine edge-ness indicators
    if ( Grid.getNumberOfNodes ( ) != reducedPCAArray[0].size ( ) ) throw aol::Exception ( "Grid dimensions and PCA array size do not match!", __FILE__, __LINE__ );
    qc::FastILexMapper<qc::QC_2D> mapper ( Grid );
    int x, y, h = ( Options.transformScale - 1 ) / 2;
    int i1, i2;
    aol::Vector<RealType> f1 ( reducedPCAArray.numComponents ( ) ), f2 ( f1, aol::STRUCT_COPY );
    aol::Vector<RealType> indicators ( reducedPCAArray[0].size ( ) );
    RealType meanIndicator = 0;
    int numValidIndicators = 0;
    for ( int i=0; i<reducedPCAArray[0].size ( ) ; ++i ) {
      mapper.splitGlobalIndex ( i, x, y );
      if ( x - h >= 0 && x + h < Grid.getNumX ( ) && y - h >= 0 && y + h < Grid.getNumY ( ) ) {
        i1 = mapper.getGlobalIndex ( x - h, y );
        i2 = mapper.getGlobalIndex ( x + h, y );
        reducedPCAArray.getTo ( i1, f1 );
        reducedPCAArray.getTo ( i2, f2 );
        f1 -= f2;
        indicators[i] += f1.normSqr ( );
        i1 = mapper.getGlobalIndex ( x, y - h );
        i2 = mapper.getGlobalIndex ( x, y + h );
        reducedPCAArray.getTo ( i1, f1 );
        reducedPCAArray.getTo ( i2, f2 );
        f1 -= f2;
        indicators[i] += f1.normSqr ( );
        
        meanIndicator += indicators[i];
        ++numValidIndicators;
      } else indicators[i] = -1;
    }
    meanIndicator /= static_cast<RealType> ( numValidIndicators );
    
    // Construct array with low edge-ness
    aol::MultiVector<RealType> reducedPCAArrayLowEdgeness ( reducedPCAArray.numComponents ( ), reducedPCAArray[0].size ( ) );
    aol::Vector<int> indexCorrespondences ( Array[0].size ( ) );
    int j = 0;
    for ( int i = 0; i<reducedPCAArray[0].size ( ); ++i ) {
      if ( indicators[i] > 0 && indicators[i] < Options.lowEdgenessThreshold * meanIndicator ) {
        for ( int c = 0; c<reducedPCAArray.numComponents ( ); ++c )
          reducedPCAArrayLowEdgeness[c][j] = reducedPCAArray[c][i];
        indexCorrespondences[j] = i;
        ++j;
      }
    }
    reducedPCAArrayLowEdgeness.resize ( reducedPCAArray.numComponents ( ), j );
    
    // Perform clustering on reduced array with low edge-ness
    aol::KMeansClusterer<RealType> kMeansClusterer;
    aol::MultiVector<RealType> clusters;
    aol::Vector<int> clusterLabels;
    
    RealType residual = 0;
    if ( Options.forceGroundTruthClusters && Options.groundTruth.size ( ) > 0 ) {
      // If clustering result is bad, use this to check if it could be improved by a better initialization of k-means
      qc::ScalarArray<int, qc::QC_2D> gt ( Options.groundTruth, aol::DEEP_COPY );
      transformGrayValuesToLabels ( gt );
      aol::MultiVector<RealType> trueCenters ( reducedPCAArrayLowEdgeness.numComponents ( ), Options.numSegments );
      for ( int l=0; l<=gt.getMaxValue ( ) ; ++l ) {
        qc::BitArray<qc::QC_2D> bitMask ( gt.getNumX ( ), gt.getNumY ( ) );
        for ( int y=1; y<gt.getNumY ( )-1 ; ++y )
          for ( int x=1; x<gt.getNumX ( )-1 ; ++x )
            if ( gt.get ( x, y ) == l ) bitMask.set ( x, y, true );
        bitMask.invert ( );
        for ( int j=0; j<12 ; ++j ) bitMask.dilateByOne ( );
        bitMask.invert ( );
        
        int n = 0;
        for ( int i=0; i<bitMask.size ( ) ; ++i ) {
          if ( indicators[i] > 0 && indicators[i] < Options.lowEdgenessThreshold * meanIndicator && bitMask[i] ) {
            for ( int c=0; c<reducedPCAArray.numComponents ( ) ; ++c )
              trueCenters[c][l] += reducedPCAArray[c][i];
            ++n;
          }
        }
        for ( int c=0; c<reducedPCAArray.numComponents ( ) ; ++c )
          trueCenters[c][l] /= static_cast<RealType> ( n );
      }
      kMeansClusterer.setInitialCenters ( trueCenters );
      residual = kMeansClusterer.apply ( reducedPCAArrayLowEdgeness, clusters, Options.numSegments, clusterLabels, aol::KMEANS_INIT_METHOD::MAN );
    } else
      residual = kMeansClusterer.applyMultipleRNG ( reducedPCAArrayLowEdgeness, clusters, Options.numSegments, clusterLabels, 1000 );
    if ( !_quietMode ) {
      std::cerr << "Residual after clustering: " << residual << std::endl;
      std::cerr << clusters << std::endl;
    }
    
    // Create segmentation labels based on clustering
    for ( int l=Segmentation.size ( ) ; l<Options.numSegments ; ++l ) Segmentation.pushBack ( Segmentation[0] ); // add required number of segments if automatically detected
    for ( int j=0; j<clusterLabels.size ( ) ; ++j ) Segmentation[clusterLabels[j]][indexCorrespondences[j]] = 1.0;
    
    if ( _outputDir != "" ) {
      for ( int l=0; l<Options.numSegments ; ++l )
        Segmentation[l].save ( aol::strprintf ( "%s/initialSegmentation_%d%s", _outputDir.c_str ( ), l, qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    }
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

template<typename _RealType, typename _PictureType>
class PiecewisePeriodicSegmentationOptions : public LocalFeatureMSSegmentationOptions<_RealType, _PictureType> {
public:
  PiecewisePeriodicSegmentationOptions ( const _PictureType &Input, const int BlockSize, const std::string OutputDir = "", const bool Quiet = true )
    : LocalFeatureMSSegmentationOptions<_RealType, _PictureType> ( Input, OutputDir, Quiet ) {
    this->operators.push_back ( new FFTModulusOp<_RealType, _PictureType> ( BlockSize ) );
    this->includeBoundary = false;
    this->operatorWeights.pushBack ( 1.0 );
    this->numDominantOperators = 1;
    this->transformScale = BlockSize;
  }
};
  
} // namespace im

#endif
