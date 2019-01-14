#include <convolution.h>

namespace qc {
  
#ifdef USE_LIB_FFTW

//! 1D Fourier transform (complex-to-complex)
//! Unnecessarily complicated because Vector<complex> is not possible
//! (otherwise copying would not be necessary)
template <typename RealType>
void FourierTransform ( const aol::MultiVector<RealType>& function, aol::MultiVector<RealType>& transform, enum FourierTransformDirection direction, bool normalize ) {
  typedef ConvolutionTrait<qc::QC_1D,RealType> ConvType;
  typedef typename ConvType::FFTWComplex FFTWComplex;
  
  
  // Prepare transformation
  if ( function.numComponents ( ) != 2 || transform.numComponents ( ) != 2 )
    throw aol::Exception ( "MultiVector components not equal to two in FourierTransform", __FILE__, __LINE__ );
  int size = function [0].size ();
  if ( function [1].size () != size || transform [0].size () != size ||  transform [1].size () != size )
    throw aol::Exception ( "Array sizes not equal in FourierTransform", __FILE__, __LINE__ );
  FFTWComplex* f = static_cast<FFTWComplex*> ( ConvType::fftwMalloc ( sizeof ( FFTWComplex ) * size ) );
  FFTWComplex* t = static_cast<FFTWComplex*> ( ConvType::fftwMalloc ( sizeof ( FFTWComplex ) * size ) );
  typename ConvType::FFTWPlan plan = ConvType::fftwPlan_dft ( size, f, t, direction, FFTW_ESTIMATE );
  
  // Copy data, transform and copy back
  for ( int i = 0; i < size; ++i )
    for ( int k = 0; k < 2; ++k ) f [i] [k] = function [k] [i];
  ConvType::fftwExecute ( plan );
  for ( int i = 0; i < size; ++i )
    for ( int k = 0; k < 2; ++k ) transform [k] [i] = t [i] [k];
  
  // Cleanup
  ConvType::fftwDestroy_plan ( plan );
  ConvType::fftwFree ( f );
  ConvType::fftwFree ( t );
  
  if ( normalize ) transform /= sqrt ( transform[0].size ( ) * 2.0 );
}
  
//! 2D Fourier transform (complex-to-complex)./
//! Unnecessarily complicated because ScalarArray<complex> is not possible
//! (otherwise copying would not be necessary)
template <typename RealType>
void FourierTransform ( const qc::MultiArray<RealType, 2, 2>& function, qc::MultiArray<RealType, 2, 2>& transform, enum FourierTransformDirection direction, bool normalize ) {
  typedef ConvolutionTrait<qc::QC_2D,RealType> ConvType;
  typedef typename ConvType::FFTWComplex FFTWComplex;


  // Prepare transformation
  int numX = function [0].getNumX (), numY = function [0].getNumY ();
  if ( function [1].getNumX () != numX || transform [0].getNumX () != numX ||  transform [1].getNumX () != numX ||
       function [1].getNumY () != numY || transform [0].getNumY () != numY ||  transform [1].getNumY () != numY )
    throw aol::Exception ( "Array sizes not equal in FourierTransform", __FILE__, __LINE__ );
  FFTWComplex* f = static_cast<FFTWComplex*> ( ConvType::fftwMalloc ( sizeof ( FFTWComplex ) * numX * numY ) );
  FFTWComplex* t = static_cast<FFTWComplex*> ( ConvType::fftwMalloc ( sizeof ( FFTWComplex ) * numX * numY ) );
  typename ConvType::FFTWPlan plan = ConvType::fftwPlan_dft ( aol::Vec2<int> ( numX, numY ), f, t, direction, FFTW_ESTIMATE );

  // Copy data, transform and copy back
  for ( int j = 0; j < numY; ++j ) {
    for ( int i = 0; i < numX; ++i ) {
      const int ind = qc::ILexCombine2 ( i, j, numX );
      const aol::Vec2<short int> pos ( i, j );
      aol::Vec2<RealType> val = function.get ( pos );
      for ( int k = 0; k < 2; ++k ) f [ind] [k] = val [k];
    }
  }
  ConvType::fftwExecute ( plan );
  for ( int j = 0; j < numY; ++j ) {
    for ( int i = 0; i < numX; ++i ) {
      const int ind = qc::ILexCombine2 ( i, j, numX );
      const aol::Vec2<short int> pos ( i, j );
      aol::Vec2<RealType> val ( t [ind] [0], t [ind] [1] );
      transform.set ( pos, val );
    }
  }

  // Cleanup
  ConvType::fftwDestroy_plan ( plan );
  ConvType::fftwFree ( f );
  ConvType::fftwFree ( t );
  
  if ( normalize ) transform /= sqrt ( transform[0].size ( ) * 4.0 );
}
  
//! 1D Cosine transform (real-to-real)
template <typename RealType>
void CosineTransform ( const aol::Vector<RealType>& function, aol::Vector<RealType>& transform, enum FourierTransformDirection direction, bool normalize ) {
  FastCosineTransform<QC_1D,RealType>::apply ( function, transform, direction, normalize );
}
  
//! 2D Cosine transform (real-to-real)
template <typename RealType>
void CosineTransform ( const qc::ScalarArray<RealType, qc::QC_2D>& function, qc::ScalarArray<RealType, qc::QC_2D>& transform, enum FourierTransformDirection direction, bool normalize ) {
  FastCosineTransform<QC_2D,RealType>::apply ( function, transform, direction, normalize );
}
  
#elif defined ( USE_KISSFFT )

//! 1D Fourier transform (complex-to-complex)
//! Unnecessarily complicated because ScalarArray<complex> is not possible
//! (otherwise copying would not be necessary)
template<>
void FourierTransform<double> ( const aol::MultiVector<double>& function, aol::MultiVector<double>& transform, enum FourierTransformDirection direction, bool normalize ) {
  typedef ConvolutionTrait<qc::QC_1D,double> ConvType;
  typedef ConvType::KISSFFTComplex KISSFFTComplex;
  
  
  // Prepare transformation
  if ( function.numComponents ( ) != 2 || transform.numComponents ( ) != 2 )
    throw aol::Exception ( "MultiVector components not equal to two in FourierTransform", __FILE__, __LINE__ );
  int size = function [0].size ();
  if ( function [1].size () != size || transform [0].size () != size ||  transform [1].size () != size )
    throw aol::Exception ( "Array sizes not equal in FourierTransform", __FILE__, __LINE__ );
  KISSFFTComplex* f = static_cast<KISSFFTComplex*> ( ConvType::kissfftMalloc ( sizeof ( KISSFFTComplex ) * size ) );
  KISSFFTComplex* t = static_cast<KISSFFTComplex*> ( ConvType::kissfftMalloc ( sizeof ( KISSFFTComplex ) * size ) );
  ConvType::KISSFFTPlan plan = ConvType::kissfftPlan_dft ( size, direction );
  
  // Copy data, transform and copy back
  for ( int i = 0; i < size; ++i ) {
    f [i].r = function [0] [i];
    f [i].i = function [1] [i];
  }
  ConvType::kissfftExecute ( plan, f, t );
  for ( int i = 0; i < size; ++i ) {
    transform [0] [i] = t [i].r;
    transform [1] [i] = t [i].i;
  }
  
  // Cleanup
  ConvType::kissfftDestroy_plan ( plan );
  ConvType::kissfftFree ( f );
  ConvType::kissfftFree ( t );
}
  
//! 2D Fourier transform (complex-to-complex)
//! Unnecessarily complicated because ScalarArray<complex> is not possible
//! (otherwise copying would not be necessary)
template<>
void FourierTransform<double> ( const qc::MultiArray<double, 2, 2>& function, qc::MultiArray<double, 2, 2>& transform, enum FourierTransformDirection direction, bool normalize ) {
  typedef ConvolutionTrait<qc::QC_2D,double> ConvType;
  typedef ConvType::KISSFFTComplex KISSFFTComplex;
  
  
  // Prepare transformation
  int numX = function [0].getNumX (), numY = function [0].getNumY ();
  if ( function [1].getNumX () != numX || transform [0].getNumX () != numX ||  transform [1].getNumX () != numX ||
      function [1].getNumY () != numY || transform [0].getNumY () != numY ||  transform [1].getNumY () != numY )
    throw aol::Exception ( "Array sizes not equal in FourierTransform", __FILE__, __LINE__ );
  KISSFFTComplex* f = static_cast<KISSFFTComplex*> ( ConvType::kissfftMalloc ( sizeof ( KISSFFTComplex ) * numX * numY ) );
  KISSFFTComplex* t = static_cast<KISSFFTComplex*> ( ConvType::kissfftMalloc ( sizeof ( KISSFFTComplex ) * numX * numY ) );
  ConvType::KISSFFTPlan plan = ConvType::kissfftPlan_dft ( aol::Vec2<int> ( numX, numY ), direction );
  
  // Copy data, transform and copy back
  for ( int j = 0; j < numY; ++j ) {
    for ( int i = 0; i < numX; ++i ) {
      const int ind = qc::ILexCombine2 ( i, j, numX );
      const aol::Vec2<short int> pos ( i, j );
      aol::Vec2<double> val = function.get ( pos );
      f [ind].r = val [0];
      f [ind].i = val [1];
    }
  }
  ConvType::kissfftExecute ( plan, f, t );
  for ( int j = 0; j < numY; ++j ) {
    for ( int i = 0; i < numX; ++i ) {
      const int ind = qc::ILexCombine2 ( i, j, numX );
      const aol::Vec2<short int> pos ( i, j );
      aol::Vec2<double> val ( t [ind].r, t [ind].i );
      transform.set ( pos, val );
    }
  }
  
  // Cleanup
  ConvType::kissfftDestroy_plan ( plan );
  ConvType::kissfftFree ( f );
  ConvType::kissfftFree ( t );
}

template<>
void FourierTransform<float> ( const aol::MultiVector<float>& /*function*/, aol::MultiVector<float>& /*transform*/, enum FourierTransformDirection /*direction*/, bool /*normalize*/ )
{
  throw aol::Exception ( "KissFFT does not support float type! Compile with -DUSE_LIB_FFTW", __FILE__, __LINE__ );
}
  
template<>
void FourierTransform<float> ( const qc::MultiArray<float, 2, 2>& /*function*/, qc::MultiArray<float, 2, 2>& /*transform*/, enum FourierTransformDirection /*direction*/, bool /*normalize*/ )
{
  throw aol::Exception ( "KissFFT does not support float type! Compile with -DUSE_LIB_FFTW", __FILE__, __LINE__ );
}
  

//! 1D Cosine transform (real-to-real)
//! Since KissFFT does not have an inherent fast DCT, we use a fast DCT based on FFT as described here: http://fourier.eng.hmc.edu/e161/lectures/dct/node2.html
template<>
void CosineTransform<double> ( const aol::Vector<double>& function, aol::Vector<double>& transform, enum FourierTransformDirection direction, bool normalize ) {
  typedef ConvolutionTrait<qc::QC_1D,double> ConvType;
  typedef ConvType::KISSFFTScalar KISSFFTScalar;
  typedef ConvType::KISSFFTComplex KISSFFTComplex;
  
  int size = function.size ();
  if ( transform.size () != size )
    throw aol::Exception ( "Array sizes not equal in CosineTransform", __FILE__, __LINE__ );
  
  if ( direction == FTForward ) {
    KISSFFTScalar* f = static_cast<KISSFFTScalar*> ( ConvType::kissfftMalloc ( sizeof ( KISSFFTScalar ) * size ) );
    KISSFFTComplex* t = static_cast<KISSFFTComplex*> ( ConvType::kissfftMalloc ( sizeof ( KISSFFTComplex ) * size ) );
    ConvType::KISSFFTRPlan plan = ConvType::kissfftrPlan_dft ( size, direction );
    
    // Generate symmetric real sequence
    for ( int i = 0; i < size / 2; ++i ) {
      f [i] = function [2*i];
      f [size-1-i] = function [2*i+1];
    }
    
    // Transform sequence
    ConvType::kissfftrExecute ( plan, f, t );
    
    // Copy real part of complex phase shifted transformed sequence
    throw aol::UnimplementedCodeException ( "Implementation incomplete!", __FILE__, __LINE__ );
    for ( int i = 0; i < size; ++i ) {
      // TODO apply phase shift by exp{-Iipi/(2size)} and extract real part
    }
    
    // Cleanup
    ConvType::kissfftrDestroy_plan ( plan );
    ConvType::kissfftFree ( f );
    ConvType::kissfftFree ( t );
  } else {
    KISSFFTComplex* f = static_cast<KISSFFTComplex*> ( ConvType::kissfftMalloc ( sizeof ( KISSFFTComplex ) * size ) );
    KISSFFTScalar* t = static_cast<KISSFFTScalar*> ( ConvType::kissfftMalloc ( sizeof ( KISSFFTScalar ) * size ) );
    ConvType::KISSFFTRPlan plan = ConvType::kissfftrPlan_dft ( size, direction );
    
    // Generate symmetric real sequence
    throw aol::UnimplementedCodeException ( "Implementation incomplete!", __FILE__, __LINE__ );
    for ( int i = 0; i < size / 2; ++i ) {
      // TODO apply phase shift by exp{Iipi/(2size)} to real input
    }
    
    // Transform sequence
    ConvType::kissfftriExecute ( plan, f, t );
    
    // Copy real part of complex phase shifted transformed sequence
    for ( int i = 0; i < size / 2; ++i ) {
      transform [2*i] = t [i];
      transform [2*i+1] = t [size-1-i];
    }
    
    // Cleanup
    ConvType::kissfftrDestroy_plan ( plan );
    ConvType::kissfftFree ( f );
    ConvType::kissfftFree ( t );
  }
}
  
//! 2D Cosine transform (real-to-real)
//! Since KissFFT does not have an inherent fast DCT, we use a fast DCT based on FFT as described here: http://fourier.eng.hmc.edu/e161/lectures/dct/node2.html
template<>
void CosineTransform<double> ( const qc::ScalarArray<double, qc::QC_2D>& function, qc::ScalarArray<double, qc::QC_2D>& transform, enum FourierTransformDirection direction, bool normalize ) {
  throw aol::UnimplementedCodeException ( "Not implemented", __FILE__, __LINE__ );
}
 
  
template<>
void CosineTransform<float> ( const aol::Vector<float>& /*function*/, aol::Vector<float>& /*transform*/, enum FourierTransformDirection /*direction*/, bool /*normalize*/ ) {
  throw aol::Exception ( "KissFFT does not support float type! Compile with -DUSE_LIB_FFTW", __FILE__, __LINE__ );
}

template<>
void CosineTransform<float> ( const qc::ScalarArray<float, qc::QC_2D>& /*function*/, qc::ScalarArray<float, qc::QC_2D>& /*transform*/, enum FourierTransformDirection /*direction*/, bool /*normalize*/ ) {
  throw aol::Exception ( "KissFFT does not support float type! Compile with -DUSE_LIB_FFTW", __FILE__, __LINE__ );
}
  
#else

//! 1D Fourier transform (complex-to-complex)
template <typename RealType>
void FourierTransform ( const aol::MultiVector<RealType>& /*function*/, aol::MultiVector<RealType>& /*transform*/, enum FourierTransformDirection /*direction*/, bool /*normalize*/ )
// this (unnamed parameter but default value) seems to work with gcc-4.2.0 and VC++ 2005. if it does not work with other compilers, write separate function with two parameters.
{
  throw aol::Exception ( "FourierTransform needs either libfftw or KissFFT! Compile with either -DUSE_LIB_FFTW or -DBUILD_AND_USE_KISSFFT", __FILE__, __LINE__ );
}
  
//! 2D Fourier transform (complex-to-complex)
template <typename RealType>
void FourierTransform ( const qc::MultiArray<RealType, 2, 2>& /*function*/, qc::MultiArray<RealType, 2, 2>& /*transform*/, enum FourierTransformDirection /*direction*/, bool /*normalize*/ )
// this (unnamed parameter but default value) seems to work with gcc-4.2.0 and VC++ 2005. if it does not work with other compilers, write separate function with two parameters.
{
  throw aol::Exception ( "FourierTransform needs either libfftw or KissFFT! Compile with either -DUSE_LIB_FFTW or -DBUILD_AND_USE_KISSFFT", __FILE__, __LINE__ );
}
                                      
                                      
//! 1D Cosine transform (real-to-real)
template <typename RealType>
void CosineTransform ( const aol::Vector<RealType>& /*function*/, aol::Vector<RealType>& /*transform*/, enum FourierTransformDirection /*direction*/, bool /*normalize*/ )
// this (unnamed parameter but default value) seems to work with gcc-4.2.0 and VC++ 2005. if it does not work with other compilers, write separate function with two parameters.
{
  throw aol::Exception ( "CosineTransform needs either libfftw or KissFFT! Compile with either -DUSE_LIB_FFTW or -DBUILD_AND_USE_KISSFFT", __FILE__, __LINE__ );
}
  
//! 2D Cosine transform (real-to-real)
template <typename RealType>
void CosineTransform ( const qc::ScalarArray<RealType, qc::QC_2D>& /*function*/, qc::ScalarArray<RealType, qc::QC_2D>& /*transform*/, enum FourierTransformDirection /*direction*/, bool /*normalize*/ )
// this (unnamed parameter but default value) seems to work with gcc-4.2.0 and VC++ 2005. if it does not work with other compilers, write separate function with two parameters.
{
  throw aol::Exception ( "CosineTransform needs either libfftw or KissFFT! Compile with either -DUSE_LIB_FFTW or -DBUILD_AND_USE_KISSFFT", __FILE__, __LINE__ );
}
  
#endif
  
  
const char* const DWT::dwts[] = {
  "haar", "db1", "db2", "db3", "db4", "db5", "db6", "db7", "db8", "db9", "db10",
  "db11", "db12", "db13", "db14", "db15",                                                 // Daubechies
  
  "bior1.1", "bior1.3", "bior1.5", "bior2.2", "bior2.4", "bior2.6", "bior2.8",
  "bior3.1", "bior3.3", "bior3.5", "bior3.7", "bior3.9", "bior4.4", "bior5.5", "bior6.8", // Bi-orthogonal wavelets
  
  "coif1", "coif2", "coif3", "coif4", "coif5"                                             // Coiflets
};
  
#ifdef USE_LIB_WAVELET
  
//! 1D Wavelet transform (real-to-real)
template <>
void WaveletTransform<double> ( const aol::Vector<double>& function, aol::Vector<double>& transform,
                                const std::string& name, enum FourierTransformDirection direction ) {
  int signalLength = function.size ( );
  if ( !aol::isPowerOfTwo ( signalLength ) ) throw aol::Exception ( "Signal length must be a power of two!", __FILE__, __LINE__ );
  if ( ! DWT::isDWT ( name ) ) throw aol::Exception ( "Invalid transform name in WavletTransform!", __FILE__, __LINE__ );
  
  if ( signalLength == 1 ) transform[0] = function[0];
  else {
    std::vector<double> signal ( signalLength ), dwtOutput, flags;
    for ( int i=0; i<signalLength ; ++i ) signal[i] = function[i];
      
    if ( direction == qc::FTForward ) dwt ( signal, log2 ( signalLength ), name, dwtOutput, flags );
    else {
      flags.resize ( 2 );
      flags[1] = log2 ( signalLength );
      idwt ( signal, flags, name, dwtOutput );
    }
      
    for ( int i=0; i<signalLength ; ++i ) transform[i] = dwtOutput[i];
  }
}
  
//! 2D Wavelet transform (real-to-real)
template <>
void WaveletTransform<double> ( const qc::ScalarArray<double, qc::QC_2D>& function, qc::ScalarArray<double, qc::QC_2D>& transform,
                                const std::string& name, enum FourierTransformDirection direction ) {
  if ( ! DWT::isDWT ( name ) ) throw aol::Exception ( "Invalid transform name in WavletTransform!", __FILE__, __LINE__ );
  
  int nx = function.getNumX ( ), ny = function.getNumY ( );
  if ( transform.getNumX ( ) != nx || transform.getNumY ( ) != ny )
    throw aol::Exception ( "Array sizes not equal in WaveletTransform", __FILE__, __LINE__ );
  if ( !aol::isPowerOfTwo ( nx ) || !aol::isPowerOfTwo ( ny ) )
    throw aol::Exception ( "Array sizes must be powers of two!", __FILE__, __LINE__ );

  // 2D wavelet transforms in wavelet1d do not seem to work
  // Instead, apply 1D transform to rows and columns consecutively
  aol::Vector<double> rowSig ( nx ), rowOut ( nx ), colSig ( ny ), colOut ( ny );
  for ( int y=0; y<ny ; ++y ) {
    for ( int x=0; x<nx ; ++x ) rowSig[x] = function.get ( x, y );
    WaveletTransform<double> ( rowSig, rowOut, name, direction );
    for ( int x=0; x<nx ; ++x ) transform.set ( x, y, rowOut[x] );
  }
  for ( int x=0; x<nx ; ++x ) {
    for ( int y=0; y<ny ; ++y ) colSig[y] = transform.get ( x, y );
    WaveletTransform<double> ( colSig, colOut, name, direction );
    for ( int y=0; y<ny ; ++y ) transform.set ( x, y, colOut[y] );
  }
}
  
template<>
void WaveletTransform<float> ( const aol::Vector<float>& function, aol::Vector<float>& transform,
                               const std::string& name, enum FourierTransformDirection direction )
{
  aol::Vector<double> functionDouble ( function.size ( ) ), transformDouble ( transform.size ( ) );
  for ( int i=0; i<function.size ( ) ; ++i )
    functionDouble[i] = static_cast<double> ( function[i] );
  WaveletTransform<double> ( functionDouble, transformDouble, name, direction );
  for ( int i=0; i<transform.size ( ) ; ++i )
    transform[i] = static_cast<float> ( transformDouble[i] );
}

template<>
void WaveletTransform<float> ( const qc::ScalarArray<float, qc::QC_2D>& function, qc::ScalarArray<float, qc::QC_2D>& transform,
                               const std::string &name, enum FourierTransformDirection direction )
{
  qc::ScalarArray<double, qc::QC_2D> functionDouble ( function.getNumX ( ), function.getNumY ( ) ), transformDouble ( functionDouble, aol::STRUCT_COPY );
  for ( int k=0; k<function.size ( ) ; ++k )
    functionDouble[k] = static_cast<float> ( function[k] );
  WaveletTransform<double> ( functionDouble, transformDouble, name, direction );
  for ( int k=0; k<transform.size ( ) ; ++k )
    transform[k] = static_cast<float> ( transformDouble[k] );
}
  
#else
          
//! 1D Wavelet transform (real-to-real)
template <typename RealType>
void WaveletTransform ( const aol::Vector<RealType>& /*function*/, aol::Vector<RealType>& /*transform*/,
                        const std::string &/*name*/, enum FourierTransformDirection /*direction*/ )
// this (unnamed parameter but default value) seems to work with gcc-4.2.0 and VC++ 2005. if it does not work with other compilers, write separate function with two parameters.
{
  throw aol::Exception ( "WaveletTransform needs wavelet1d! Compile with -DBUILD_AND_USE_WAVELET", __FILE__, __LINE__ );
}
  
//! 2D Wavelet transform (real-to-real)
template <typename RealType>
void WaveletTransform ( const qc::ScalarArray<RealType, qc::QC_2D>& /*function*/, qc::ScalarArray<RealType, qc::QC_2D>& /*transform*/,
                        const std::string &/*name*/, enum FourierTransformDirection /*direction*/ )
// this (unnamed parameter but default value) seems to work with gcc-4.2.0 and VC++ 2005. if it does not work with other compilers, write separate function with two parameters.
{
  throw aol::Exception ( "WaveletTransform needs wavelet1d! Compile with -DBUILD_AND_USE_WAVELET", __FILE__, __LINE__ );
}
  
template void WaveletTransform<float> ( const aol::Vector<float>&, aol::Vector<float>&, const std::string&, enum FourierTransformDirection );
template void WaveletTransform<double> ( const aol::Vector<double>&, aol::Vector<double>&, const std::string&, enum FourierTransformDirection );
template void WaveletTransform<float> ( const qc::ScalarArray<float, qc::QC_2D>&, qc::ScalarArray<float, qc::QC_2D>&, const std::string&, enum FourierTransformDirection );
template void WaveletTransform<double> ( const qc::ScalarArray<double, qc::QC_2D>&, qc::ScalarArray<double, qc::QC_2D>&, const std::string&, enum FourierTransformDirection );
  
#endif

#if !defined ( USE_KISSFFT )
template void FourierTransform<float> ( const aol::MultiVector<float>&, aol::MultiVector<float>&, enum FourierTransformDirection, bool );
template void FourierTransform<double> ( const aol::MultiVector<double>&, aol::MultiVector<double>&, enum FourierTransformDirection, bool );
template void FourierTransform<float> ( const qc::MultiArray<float, 2, 2>&, qc::MultiArray<float, 2, 2>&, enum FourierTransformDirection, bool );
template void FourierTransform<double> ( const qc::MultiArray<double, 2, 2>&, qc::MultiArray<double, 2, 2>&, enum FourierTransformDirection, bool );
// fftw3l is not available everywhere.
#if defined ( HAVE_LIB_FFTW_LONGDOUBLE )
template void FourierTransform<long double> ( const aol::MultiVector<long double>&, aol::MultiVector<long double>&, enum FourierTransformDirection, bool );
template void FourierTransform<long double> ( const qc::MultiArray<long double, 2, 2>&, qc::MultiArray<long double, 2, 2>&, enum FourierTransformDirection, bool );
#endif

template void CosineTransform<float> ( const aol::Vector<float>&, aol::Vector<float>&, enum FourierTransformDirection, bool );
template void CosineTransform<double> ( const aol::Vector<double>&, aol::Vector<double>&, enum FourierTransformDirection, bool );
template void CosineTransform<float> ( const qc::ScalarArray<float, qc::QC_2D>&, qc::ScalarArray<float, qc::QC_2D>&, enum FourierTransformDirection, bool );
template void CosineTransform<double> ( const qc::ScalarArray<double, qc::QC_2D>&, qc::ScalarArray<double, qc::QC_2D>&, enum FourierTransformDirection, bool );
#endif

void addMotionBlurToArray ( const aol::Vec2<double> &Velocity, const qc::ScalarArray<double, qc::QC_2D> &Arg, qc::ScalarArray<double, qc::QC_2D> &Dest ) {
  qc::Convolution<qc::QC_2D> conv ( aol::Vec2<int>( Arg.getNumX(), Arg.getNumY() ) );
  qc::ScalarArray<double, qc::QC_2D> kernel ( Arg, aol::STRUCT_COPY );
  qc::generateMotionBlurKernel<double> ( Velocity, kernel );
  conv.convolve ( Arg, kernel, Dest );
}


} // end namespace
