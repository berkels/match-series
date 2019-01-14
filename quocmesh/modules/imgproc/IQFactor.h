#ifndef __IQFACTOR_H
#define __IQFACTOR_H

#include <atomFinder.h>

namespace im {

/**
 * \author Berkels
 */
template <typename RealType>
class IQFactor {
  typedef qc::RectangularGridConfigurator<RealType, qc::QC_2D, aol::GaussQuadrature<RealType,qc::QC_2D,3> > BumpFitConfType;
  std::vector<std::pair<RealType, RealType> > _IQFactors;
  qc::ScalarArray<RealType, qc::QC_2D> _modulus;
  aol::RandomAccessContainer<aol::Vec2<short> > _maxima;
  aol::BitVector _maximaSelected;
  const int _patchSize;
  RealType _IQValueThreshold;
  const RealType _IQFrequencyThreshold;
  const RealType _onePixelInAngstrom;
  const bool _subtractMeanBeforeModulus;
  const int _eraseInfinityRadius;
  const bool _bumpFitPeaks;
  const RealType _minPeakIntensityRatio;
  const RealType _minPeakBackgroundRatio;
public:
  IQFactor ( const qc::ScalarArray<RealType, qc::QC_2D> &Image,
             const int PatchSize,
             const RealType IQValueThreshold,
             const RealType OnePixelInAngstrom,
             const bool SubtractMeanBeforeModulus,
             const bool BumpFitPeaks,
             const RealType IQFrequencyThreshold = aol::NumberTrait<RealType>::Inf,
             const int EraseInfinityRadius = 30,
             const RealType MinPeakIntensityRatio = 0.01,
             const RealType MinPeakBackgroundRatio = 1 )
    : _modulus ( Image, aol::STRUCT_COPY ),
      _patchSize ( PatchSize ),
      _IQValueThreshold ( IQValueThreshold ),
      _IQFrequencyThreshold ( IQFrequencyThreshold ),
      _onePixelInAngstrom ( OnePixelInAngstrom ),
      _subtractMeanBeforeModulus ( SubtractMeanBeforeModulus ),
      _eraseInfinityRadius ( EraseInfinityRadius ),
      _bumpFitPeaks ( BumpFitPeaks ),
      _minPeakIntensityRatio ( MinPeakIntensityRatio ),
      _minPeakBackgroundRatio ( MinPeakBackgroundRatio ) {
    qc::computeLogFFTModulus<RealType> ( Image, _modulus, 0, false, SubtractMeanBeforeModulus );
    findModulusPeaks ( );

    // The modulus of a real image is symmetric, so we can drop maxima in the left quadrants.
    calcIQFactorsPatchBased ( true );
  }

  void rerunUsingOldPeakEstimate ( const qc::ScalarArray<RealType, qc::QC_2D> &Image ) {
    qc::computeLogFFTModulus<RealType> ( Image, _modulus, 0, false, _subtractMeanBeforeModulus );
    if ( _bumpFitPeaks == false )
      refinePeaksByLocalMaxima ( );
    calcIQFactorsPatchBased ( true );
  }

  const std::vector<std::pair<RealType, RealType> > &getIQFactorsRef ( ) const {
    return _IQFactors;
  }

  const qc::ScalarArray<RealType, qc::QC_2D> &getModulusRef ( ) const {
    return _modulus;
  }

  const aol::RandomAccessContainer<aol::Vec2<short> > &getMaximaRef ( ) const {
    return _maxima;
  }

  const aol::BitVector &getMaximaSelectedRef ( ) const {
    return _maximaSelected;
  }

  void setMaxima ( const aol::RandomAccessContainer<aol::Vec2<short> > &Maxima, const aol::BitVector &MaximaSelected ) {
    _maxima = Maxima;
    _maximaSelected = MaximaSelected;
  }

  void setIQValueThreshold ( const RealType IQValueThreshold ) {
    _IQValueThreshold = IQValueThreshold;
  }

  void plotIQFactors  ( const char *BaseOutName ) const {
    aol::Plotter<RealType> plotter;
    plotter.set_outfile_base_name( BaseOutName );
    aol::PlotDataFileHandler<RealType> plotHandler;
    plotHandler.generateFunctionPlot( _IQFactors );
    plotter.addPlotCommandsFromHandler( plotHandler );
    plotter.genPlot( aol::GNUPLOT_PNG );
  }

  void exportIQFactorsAsASCII ( const char *Filename ) const {
    std::ofstream outdat ( Filename );
    for ( unsigned int i = 0; i < _IQFactors.size(); ++i )
      outdat << _IQFactors[i].first << " " << _IQFactors[i].second << endl;
    outdat.close();
  }

  RealType getIQValueSum ( ) const {
    RealType sum = 0;
    for ( unsigned int i = 0; i < _IQFactors.size(); ++i )
      sum += _IQFactors[i].second;
    return sum;
  }

  void visualizePeaksOnModulus ( qc::MultiArray<double, qc::QC_2D, 3> &Image ) const {
    Image[0] = _modulus;
    aol::Vec2<RealType> minMax = Image[0].getSaturatedMinMaxValue ( 0.15 );
    Image[0].clamp ( minMax[0], minMax[1] );
    for ( int c = 1; c < 3; ++c )
      Image[c] = _modulus;

    const RealType maxVal = Image[0].getMaxValue();

    for ( int i = 0; i < _maxima.size(); ++i ) {
      for ( int y = aol::Max ( _maxima[i][1] - 2, 0 ); y <= aol::Min ( _maxima[i][1] + 2, Image.getNumY() - 1 ); ++y ) {
        for ( int x = aol::Max ( _maxima[i][0] - 2, 0 ); x <= aol::Min ( _maxima[i][0] + 2, Image.getNumX() - 1 ); ++x ) {
          Image[_maximaSelected[i] ? 0 : 2].set ( x, y, maxVal );
          Image[_maximaSelected[i] ? 2 : 0].set ( x, y, 0 );
          Image[1].set ( x, y, 0 );
        }
      }
    }
  }

private:
  void findModulusPeaks ( ) {
    _maxima.clear();
    if ( ( _bumpFitPeaks == false ) && ( _eraseInfinityRadius == 0 ) ) {
      qc::findMaxima<RealType> ( _modulus, _patchSize, _maxima );
    }
    else {
      estimatePeaks ( );
    }
  }

  void estimatePeaks ( ) {
    aol::MultiVector<RealType> estimatedCenters;
    const RealType smoothSigma = 0;
    const RealType minCenterValue = _minPeakIntensityRatio * ( _modulus.getMaxValue() - _modulus.getMinValue() ) + _modulus.getMinValue();
    im::BumpFitCenterEstimator<BumpFitConfType> centerEstimator ( smoothSigma, _eraseInfinityRadius, minCenterValue );
    centerEstimator.guessCenters ( _modulus, estimatedCenters );

    _maxima.clear();
    cerr << "Initial guess contains " << estimatedCenters.numComponents() << " centers." << endl;
    qc::ScalarArray<RealType, qc::QC_2D> temp ( 2*_patchSize, 2*_patchSize );
    for ( int i = 0; i < estimatedCenters.numComponents(); ++i ) {
      aol::Vec2<short> peakPos ( static_cast<short> ( estimatedCenters[i][0] ),
                                 static_cast<short> ( estimatedCenters[i][1] ) );

      if ( ( peakPos.getMinValue() - _patchSize ) < 0 )
        continue;
      if ( ( peakPos[0] + _patchSize ) >= _modulus.getNumX() )
        continue;
      if ( ( peakPos[1] + _patchSize ) >= _modulus.getNumY() )
        continue;

      _modulus.copyBlockTo ( aol::Vec2<short> ( peakPos[0] - _patchSize, peakPos[1] - _patchSize ), temp );
      if ( _modulus.get ( peakPos ) > _minPeakBackgroundRatio * temp.getMeanValue() )
        _maxima.pushBack ( peakPos );
    }
    _maximaSelected.reallocate( _maxima.size() );
    _maximaSelected.setAll( true );
  }

  void refinePeaksByLocalMaxima ( ) {
    if ( _eraseInfinityRadius == 0 )
      throw aol::ParameterRangeException( "im::IQFactor::refinePeaksByLocalMaxima only makes sense for positive erase radius.", __FILE__, __LINE__ );

    for ( int k = 0; k < _maxima.size(); ++k ) {
      if ( _maximaSelected.get ( k ) == false )
        continue;

      const aol::Vec2<short> pos = _maxima[k];
      RealType maxVal = _modulus.get ( pos );

      for ( int i = aol::Max ( 0, pos[0] - _eraseInfinityRadius/2 ); i < aol::Min ( _modulus.getNumX(), pos[0]+_eraseInfinityRadius/2 ); ++i ) {
        for ( int j = aol::Max ( 0, pos[1] - _eraseInfinityRadius/2 ); j < aol::Min ( _modulus.getNumY(), pos[1]+_eraseInfinityRadius/2 ); ++j ) {
          const RealType val = _modulus.get ( i, j );
          if ( val > maxVal ) {
            maxVal = val;
            _maxima[k][0] = i;
            _maxima[k][1] = j;
          }
        }
      }
    }
  }

  void calcIQFactorsPatchBased ( const bool DropLeftQuadrants ) {
    aol::MultiVector<RealType> gaussianParams;
    aol::MultiVector<RealType> centers;
    if ( _bumpFitPeaks ) {
      aol::MultiVector<RealType> estimatedCenters ( _maxima.size(), 2 );
      for ( int i = 0; i < _maxima.size(); ++i ) {
        estimatedCenters[i][0] = _maxima[i][0];
        estimatedCenters[i][1] = _maxima[i][1];
      }

      im::AtomFinder<RealType> atomFinder;
      aol::MultiVector<int> approximateDimensions ( estimatedCenters.numComponents(), 2 );
      approximateDimensions.setAll( _eraseInfinityRadius );
      atomFinder.getRefinedAtomPositions ( centers, gaussianParams, estimatedCenters, approximateDimensions, _modulus, false );
    }

    qc::ScalarArray<RealType, qc::QC_2D> temp ( _patchSize, _patchSize );
    const RealType hx = aol::ZOTrait<RealType>::one / ( _modulus.getNumX() - 1 );
    const RealType hy = aol::ZOTrait<RealType>::one / ( _modulus.getNumY() - 1 );
    _IQFactors.clear();
    for ( int i = 0; i < _maxima.size(); ++i ) {
      // Skip maxima we already ignored earlier.
      if ( _maximaSelected.get ( i ) == false )
        continue;

      // There are several reasons to discard a maxima, but only one place
      // where it is finally accepted. We mark the accepted ones there.
      _maximaSelected.set ( i, false );

      aol::Vector<RealType> backgroundValues;
      backgroundValues.reserve ( 4 );

      const aol::Vec2<short> peakPos = _bumpFitPeaks ? aol::Vec2<short>  ( static_cast<short> ( centers[i][0] ),
                                                                           static_cast<short> ( centers[i][1] ) ) : _maxima[i];

      // Ignore maxima that are so close to the boundary that one of the background
      // sample areas is not in the image domain anymore.
      if ( ( peakPos.getMinValue() - _patchSize * 3 / 2 ) < 0 )
        continue;
      if ( ( peakPos[0] + _patchSize * 3 / 2 ) >= _modulus.getNumX() )
        continue;
      if ( ( peakPos[1] + _patchSize * 3 / 2 ) >= _modulus.getNumY() )
        continue;

      for ( int k = -1; k <= 1; k = k + 2 ) {
        for ( int l = -1; l <= 1; l = l + 2 ) {
          _modulus.copyBlockTo ( aol::Vec2<short> ( peakPos[0] - _patchSize/2 + k * _patchSize, peakPos[1] - _patchSize/2 + l * _patchSize ), temp );
          backgroundValues.pushBack ( temp.getMeanValue() );
        }
      }
      // Here one could try to remove bumps that are too elongated.
      //if ( gaussianParams[i][3] / gaussianParams[i][4] > 0.1 )
      const RealType IQ = _bumpFitPeaks ? ( gaussianParams[i][2] / gaussianParams[i][6] ) : ( _modulus.get ( peakPos ) / backgroundValues.getMeanValue() - 1 );
      // In case the left quadrants are dropped, try not to drop the "center" peak, which may be one pixel off.
      if ( ( IQ > _IQValueThreshold ) && ( ( DropLeftQuadrants == false ) || ( _modulus.getNumX() / 2 - 1 <= peakPos[0] ) ) ) {
        const RealType frequency = aol::Vec2<RealType> ( ( _modulus.getNumX() / 2 - peakPos[0] ) * hx / _onePixelInAngstrom, ( _modulus.getNumY() / 2 - peakPos[1] ) * hy / _onePixelInAngstrom ).norm();
        if ( frequency < _IQFrequencyThreshold ) {
          _IQFactors.push_back ( std::pair<RealType, RealType> ( frequency, IQ ) );
          _maximaSelected.set ( i, true );
        }
      }
    }
    std::sort( _IQFactors.begin(), _IQFactors.end() );
  }

};

} // end of namespace im.

#endif // __IQFACTOR_H
