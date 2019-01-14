#ifndef __DENOISING_H
#define __DENOISING_H

#include <quoc.h>
#include <firstOrderTVAlgos.h>

namespace im {

/**
 * Based on http://stackoverflow.com/questions/5695865/bilateral-filter
 *
 * \author Berkels
 */
template <typename RealType>
void bilateralFilter ( const qc::ScalarArray<RealType, qc::QC_2D> &Input,
                       qc::ScalarArray<RealType, qc::QC_2D> &Output,
                       const int KernelRadius,
                       const RealType Sigma,
                       const RealType SigmaIntensity ) {
  const int numX = Input.getNumX();
  const int numY = Input.getNumY();
  for ( int y = 0; y < numY; ++y ) {
    for ( int x = 0; x < numX; ++x ) {

      RealType sumWeight = 0;
      RealType sum = 0;

      const RealType ctrPix = Input.get( x, y );

      const int kernelStartX = aol::Max ( x - KernelRadius, 0 );
      const int kernelEndX   = aol::Min ( x + KernelRadius, numX-1 );
      const int kernelStartY = aol::Max ( y - KernelRadius, 0 );
      const int kernelEndY   = aol::Min ( y + KernelRadius, numY-1 );

      for ( int j = kernelStartY; j <= kernelEndY; ++j ) {
        for ( int i = kernelStartX; i <= kernelEndX; ++i ) {

          const RealType curPix = Input.get( i, j );
          const RealType imageDist = sqrt ( static_cast<RealType> ( aol::Sqr ( i - x ) + aol::Sqr ( j - y ) ) );
          const RealType colorDist = sqrt ( static_cast<RealType> ( aol::Sqr ( curPix - ctrPix ) ) );

          const RealType currWeight = exp ( - aol::Sqr ( imageDist / Sigma ) * 0.5 ) * exp ( - aol::Sqr ( colorDist / SigmaIntensity ) * 0.5 );
          sumWeight += currWeight;

          sum += currWeight * curPix;
        }
      }
      Output.set ( x, y, sum / sumWeight );
    }
  }
}

/**
 * 1D version of the 2D bilateral filter above
 *
 * \author Berkels
 */
template <typename RealType>
void bilateralFilter ( const aol::Vector<RealType> &Input,
                       aol::Vector<RealType> &Output,
                       const int KernelRadius,
                       const RealType Sigma,
                       const RealType SigmaIntensity ) {
  const int length = Input.size();
  for ( int x = 0; x < length; ++x ) {
    
    RealType sumWeight = 0;
    RealType sum = 0;
    
    const RealType ctrPix = Input[x];
    
    const int kernelStartX = aol::Max ( x - KernelRadius, 0 );
    const int kernelEndX   = aol::Min ( x + KernelRadius, length-1 );
  
    for ( int i = kernelStartX; i <= kernelEndX; ++i ) {
      
      const RealType curPix = Input[i];
      const RealType imageDist = sqrt ( static_cast<RealType> ( aol::Sqr ( i - x ) ) );
      const RealType colorDist = sqrt ( static_cast<RealType> ( aol::Sqr ( curPix - ctrPix ) ) );
      
      const RealType currWeight = exp ( - aol::Sqr ( imageDist / Sigma ) * 0.5 ) * exp ( - aol::Sqr ( colorDist / SigmaIntensity ) * 0.5 );
      sumWeight += currWeight;
      
      sum += currWeight * curPix;
    }
    Output[x] = sum / sumWeight;
  }
}
  
  
/**
 * 1D median filter
 *
 * \author Mevenkamp
 */
template <typename RealType>
void medianFilter ( const aol::Vector<RealType> &Input,
                    aol::Vector<RealType> &Output,
                    const int FilterSize = 3 ) {
  const int Off = ( FilterSize - 1 ) >> 1;
  const int length = Input.size();
  for ( int x = 0; x < length; ++x ) {
    const int XMin = aol::Max ( 0, x - Off );
    const int XMax = aol::Min ( length - 1, x + Off );
    aol::Vector<RealType> tmp ( XMax - XMin + 1 );
    for ( int xx=XMin; xx<=XMax ; ++xx )
      tmp[xx-XMin] = Input[xx];
    Output[x] = tmp.getMedianValue ( );
  }
}
  

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class FirstOrderPrimalDualROFMinimizer : public qc::FirstOrderChambollePockTVAlgorithmType2<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  
private:
  const ArrayType &_image;
  
  void applyResolventOfDataTermSingle ( const RealType TauOverGamma, ArrayType &ArgDest ) const {
    ArgDest.addMultiple ( _image, TauOverGamma );
    ArgDest /= ( 1 + TauOverGamma );
  }
  
public:
  FirstOrderPrimalDualROFMinimizer ( const typename ConfiguratorType::InitType &Initializer,
                                    const RealType Gamma,
                                    const ArrayType &Image,
                                    const int MaxIterations = 1000,
                                    const RealType StopEpsilon = 0 )
  : qc::FirstOrderChambollePockTVAlgorithmType2<ConfiguratorType> ( Initializer, Gamma, MaxIterations, StopEpsilon ),
  _image ( Image ) {}
};

/**
 * \brief Inpaint an image by trying to detect and recover broken pixels
 * 
 * Up to four kinds of broken pixels can be considered:
 *  1) Pixels with unusually small values (inpaintsmall)
 *  2) Pixels with unusually large values (inpaintlarge)
 *  3) Pixels with a value of zero (inpaintzero)
 *  4) Pixels with a negative value (inpaintnegative)
 * 
 * \author Berkels, Doberstein
 */
template <typename RealType>
void inpaintImage ( const qc::ScalarArray<RealType, qc::QC_2D>& input,
                    qc::ScalarArray<RealType, qc::QC_2D>& inpaintedImage,
                    qc::BitArray<qc::QC_2D>& brokenPixels,
                    const bool inpaintsmall = false,
                    const RealType inpaint_small_factor = 0.8,
                    const bool inpaintlarge = true,
                    const RealType inpaint_large_factor = 1.5,
                    const bool inpaintzero = false,
                    const bool inpaintnegative = false,
                    const int radius = 1,
                    const bool median = true,
                    const bool verbose = false,
                    const qc::ScalarArray<RealType, qc::QC_2D>& inputMask = qc::ScalarArray<RealType, qc::QC_2D> ( ) ) {
  // Initialize inpaintedImage
  if ( inpaintedImage.getNumX ( ) != input.getNumX ( ) || inpaintedImage.getNumY ( ) != input.getNumY ( ) )
    inpaintedImage.reallocate ( input.getNumX ( ), input.getNumY ( ) );
  inpaintedImage = input;
  
  // Reallocate bit array and set all entries to zero
  brokenPixels.reallocate ( input.getNumX ( ), input.getNumY ( ) );
  
  // Assemble broken pixel mask
  if ( inputMask.size ( ) != 0 ) {
    brokenPixels.thresholdFrom ( inputMask, inputMask.getMinMaxValue().getMeanValue() );
  } else {
    const aol::Vec2<RealType> saturatedMinMax = input.getSaturatedMinMaxValue ( 1 );
    
    for ( int i = 0; i < input.size(); ++i ) {
      // If desired, pixels with very low values are considered to be "dead".
      if ( inpaintsmall && ( input[i] < inpaint_small_factor * saturatedMinMax[0] ) )
        brokenPixels.set ( i, true );
      
      // If desired, pixels with very large values are considered to be "hot".
      if ( inpaintlarge && ( input[i] > inpaint_large_factor * saturatedMinMax[1] ) )
        brokenPixels.set ( i, true );
      
      // If desired, pixels with value 0 are considered to be "dead".
      if ( inpaintzero && ( input[i] == 0 ) )
        brokenPixels.set ( i, true );

      // If desired, also mark pixels with a negative value as broken.
      if ( inpaintnegative && ( input[i] < 0 ) )
        brokenPixels.set ( i, true );
    }
  }
  
  // Inpainting
  const int max_radius = aol::Max ( input.getNumX ( ), input.getNumY ( ) );
  for ( int y = 0; y < input.getNumY(); ++y ) {
    for ( int x = 0; x < input.getNumX(); ++x ) {
      if ( brokenPixels.get ( x, y ) == false )
        continue;

      int current_radius = radius;
      while ( current_radius < max_radius ) {
        aol::Vector<RealType> neighborhoodValues;
        for ( int j = aol::Max ( 0, y - current_radius ); j <= aol::Min ( input.getNumY() - 1, y + current_radius ); ++j ) {
          for ( int i = aol::Max ( 0, x - current_radius ); i <= aol::Min ( input.getNumX() - 1, x + current_radius ); ++i ) {
            if ( brokenPixels.get ( i, j ) == false )
              neighborhoodValues.pushBack ( input.get ( i, j ) );
          }
        }
        
        if ( neighborhoodValues.size ( ) == 0 ) {
          ++current_radius;
          if ( verbose )
            cerr << "Warning: Increased radius for pixel (" << x << ", " << y << ") from " << current_radius - 1 << " to " << current_radius << "." << endl;
          continue;
        }
        
        inpaintedImage.set ( x, y, median ?  neighborhoodValues.getMedianValue() : neighborhoodValues.getMeanValue() );
        break;
      }

    }
  }
  
}


} // namespace im

#endif
