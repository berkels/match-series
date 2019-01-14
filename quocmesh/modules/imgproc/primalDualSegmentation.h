#ifndef __PRIMALDUALSEGMENTATION_H
#define __PRIMALDUALSEGMENTATION_H

#include <quoc.h>
#include <firstOrderTVAlgos.h>
#include <segmentation.h>

namespace im {

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class FirstOrderPrimalDualTwoPhaseMSSegmentorEngine : public qc::FirstOrderChambollePockTVAlgorithmType2<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  
private:
  const ArrayType &_indicator1Plus2;
  const ArrayType &_indicator2;
  
protected:
  void applyResolventOfDataTermSingle ( const RealType TauOverGamma, ArrayType &ArgDest ) const {
    const int dofs = ArgDest.size();
    const RealType twoTauOverGamma = 2*TauOverGamma;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < dofs; ++i )
      ArgDest[i] = ( ArgDest[i] + twoTauOverGamma*_indicator2[i] ) / ( 1 + twoTauOverGamma*_indicator1Plus2[i] );
  }
  
public:
  FirstOrderPrimalDualTwoPhaseMSSegmentorEngine ( const typename ConfiguratorType::InitType &Initializer,
                                                 const RealType Gamma,
                                                 const ArrayType &Indicator1Plus2,
                                                 const ArrayType &Indicator2 )
  : qc::FirstOrderChambollePockTVAlgorithmType2<ConfiguratorType> ( Initializer, Gamma ),
  _indicator1Plus2 ( Indicator1Plus2 ),
  _indicator2 ( Indicator2 ) {}
};

/**
 * \author Berkels
 * \ingroup Segmentation
 */
template <typename ConfiguratorType>
class FirstOrderPrimalDualTwoPhaseMSSegmentor : public TwoPhaseMSSegmentor<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  
  FirstOrderPrimalDualTwoPhaseMSSegmentor ( const typename ConfiguratorType::InitType &Initializer,
                                           const RealType Gamma )
  : TwoPhaseMSSegmentor<ConfiguratorType> ( Initializer, Gamma ) {}
  
private:
  void doSegment ( ArrayType &Segmentation, qc::MultiArray<RealType, ConfiguratorType::Dim> *PDual ) const {
    aol::VectorContainer<ArrayType> indicators;
    this->generateIndicatorFunctions ( indicators );
    indicators[0] += indicators[1];
    
    FirstOrderPrimalDualTwoPhaseMSSegmentorEngine<ConfiguratorType> engine ( this->_grid, this->_gamma, indicators[0], indicators[1] );
    engine.setMaxIterations ( this->getMaxIterations ( ) );
    engine.setStopEpsilon ( this->getStopEpsilon( ) );
    engine.setQuietMode ( this->_quietMode );
    if ( this->getStepSaverPointer ( ) )
      engine.setStepSaverReference ( *(this->getStepSaverPointer ( )) );
    engine.minimize ( Segmentation, PDual );
  }
};
  
/**
 * \author Mevenkamp
 */
template <typename ConfiguratorType>
class FirstOrderPrimalDualMultiPhaseMSSegmentorEngine : public qc::FirstOrderChambollePockTVAlgorithmType1<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  
private:
  const aol::VectorContainer<ArrayType> &_indicators;
  
protected:
  void applyResolventOfDataTermSingle ( const RealType TauOverGamma, aol::VectorContainer<ArrayType> &ArgDest ) const {
    if ( ArgDest.size ( ) != _indicators.size ( ) || ArgDest[0].size ( ) != _indicators[0].size ( ) )
      throw aol::Exception ( "Resolvent argument dimension does not match indicator dimension!", __FILE__, __LINE__ );
    
    const int dofs = ArgDest[0].size();
    aol::CanonicalSimplexProjector<RealType, aol::Vector<RealType> > canonicalSimplexProjector;
    aol::Vector<RealType> projArg ( ArgDest.size ( ) ), projDest ( ArgDest.size ( ) );
#ifdef _OPENMP
#pragma omp parallel for firstprivate ( projArg, projDest )
#endif
    for ( int i = 0; i < dofs; ++i ) {
      for ( int k = 0; k < ArgDest.size ( ) ; ++k )
        projArg[k] = ArgDest[k][i] - TauOverGamma*_indicators[k][i];
      canonicalSimplexProjector.apply ( projArg, projDest );
      for ( int k = 0; k < ArgDest.size ( ) ; ++k )
        ArgDest[k][i] = projDest[k];
    }
  }
  
public:
  FirstOrderPrimalDualMultiPhaseMSSegmentorEngine ( const typename ConfiguratorType::InitType &Initializer,
                                                   const RealType Gamma,
                                                   const aol::VectorContainer<ArrayType> &Indicators )
  : qc::FirstOrderChambollePockTVAlgorithmType1<ConfiguratorType> ( Initializer, Gamma ),
  _indicators ( Indicators ) {}
};

/**
 * \author Mevenkamp
 * \ingroup Segmentation
 */
template <typename ConfiguratorType>
class FirstOrderPrimalDualMultiPhaseMSSegmentor : public MultiPhaseMSSegmentor<ConfiguratorType> {
protected:
  int _numGhostCells;
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  
  FirstOrderPrimalDualMultiPhaseMSSegmentor ( const typename ConfiguratorType::InitType &Initializer,
                                             const RealType Gamma,
                                             const int NumSegments = 0,
                                             const int NumGhostCells = 0 )
  : MultiPhaseMSSegmentor<ConfiguratorType> ( Initializer, Gamma, NumSegments ), _numGhostCells ( NumGhostCells ) {}
  
private:
  void doSegment ( aol::VectorContainer<ArrayType> &Segmentation, aol::VectorContainer<qc::MultiArray<RealType, ConfiguratorType::Dim> > *PDual ) const {
    if ( this->_unknownRegion && Segmentation.size ( ) == this->_numSegments ) {
      Segmentation.pushBack ( ArrayType ( this->_grid ) );
      for ( int l=0; l<this->_numSegments ; ++l ) Segmentation[l+1] = Segmentation[l];
      Segmentation[0].setZero ( );
    }
    
    if ( _numGhostCells == 0 ) {
      aol::VectorContainer<ArrayType> indicators;
      this->generateIndicatorFunctions ( indicators );
      
      // Compute segmentation
      FirstOrderPrimalDualMultiPhaseMSSegmentorEngine<ConfiguratorType> engine ( this->_grid, this->_gamma, indicators );
      engine.setMaxIterations ( this->getMaxIterations ( ) );
      engine.setStopEpsilon ( this->getStopEpsilon( ) );
      engine.setQuietMode ( this->_quietMode );
      if ( this->getStepSaverPointer ( ) )
        engine.setStepSaverReference ( *(this->getStepSaverPointer ( )) );
      engine.minimize ( Segmentation, PDual );
    } else doSegmentWithGhostCells ( Segmentation, PDual );
  }
  
  void doSegmentWithGhostCells ( aol::VectorContainer<ArrayType> &Segmentation, aol::VectorContainer<qc::MultiArray<RealType, ConfiguratorType::Dim> > *PDual ) const {
    Segmentation.reallocate ( this->_numSegments + ( this->_unknownRegion ? 1 : 0 ), ArrayType ( this->_grid ) );
    
    aol::VectorContainer<ArrayType> indicators;
    this->generateIndicatorFunctions ( indicators );
    
    // Zero-extend indicators, dual and segmentation to a larger grid (according to the number of specified ghost cells)
    std::cerr << "Zero-extending indicators by " << _numGhostCells << " ghost cells." << std::endl;
    typename ConfiguratorType::InitType extendedGrid ( qc::GridSize<ConfiguratorType::Dim> ( this->_grid.getWidth ( ) + 2 * _numGhostCells, this->_grid.getHeight ( ) + 2 * _numGhostCells ) );
    aol::VectorContainer<ArrayType> extendedIndicators ( this->_numSegments, ArrayType ( extendedGrid ) );
    aol::VectorContainer<qc::MultiArray<RealType, ConfiguratorType::Dim> > extendedPDual ( this->_numSegments, qc::MultiArray<RealType, ConfiguratorType::Dim> ( extendedGrid ) );
    aol::VectorContainer<ArrayType> extendedSegmentation ( this->_numSegments, ArrayType ( extendedGrid ) );
    for ( int l=0; l<this->_numSegments ; ++l ) {
      extendedIndicators[l].padFrom ( indicators[l] );
      if ( PDual != NULL ) {
        for ( int i=0; i<ConfiguratorType::Dim ; ++i )
          extendedPDual[l][i].padFrom ( (*PDual)[l][i] );
      }
    }
    
    // Compute segmentation
    FirstOrderPrimalDualMultiPhaseMSSegmentorEngine<ConfiguratorType> engine ( extendedGrid, this->_gamma, extendedIndicators );
    engine.setMaxIterations ( this->getMaxIterations ( ) );
    engine.setStopEpsilon ( this->getStopEpsilon( ) );
    engine.setQuietMode ( this->_quietMode );
    if ( this->getStepSaverPointer ( ) )
      engine.setStepSaverReference ( *(this->getStepSaverPointer ( )) );
    engine.minimize ( extendedSegmentation, &extendedPDual );
    
    // Extract segmentation and dual from their extended versions
    for ( int l=0; l<this->_numSegments ; ++l ) {
      Segmentation[l].padFrom ( extendedSegmentation[l] );
      if ( PDual != NULL ) {
        for ( int i=0; i<ConfiguratorType::Dim ; ++i )
          (*PDual)[l][i].padFrom ( extendedPDual[l][i] );
      }
    }
  }
};
  
} // namespace im

#endif
