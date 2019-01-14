#include <kernel1d.h>

template<typename RealType>
void qc::Kernel1d<RealType>::dump() const {
  for ( int X = -this->offset; X <= this->offset; X++ ) {
    cerr << setw ( 6 ) << getValue ( X ) << " ";
  }
  cerr << endl;
}

template<typename RealType>
void qc::GaussKernel1d<RealType>::makeKernel() {
  for ( int X = -this->offset; X <= this->offset; X++ ) {
    RealType val = exp ( -0.5f * ( X * X ) / ( sigma * sigma ) );
    this->setValue ( X, val );
  }
  this->normalize();
}


template<typename RealType>
void qc::GaussDiffKernel1d<RealType>::makeKernel() {

  const RealType scale = 1 / ( 2*aol::NumberTrait<RealType>::pi * aol::Sqr(sigma) );
  const RealType s = 2 * aol::Sqr(sigma);

  for ( int X = -this->offset; X <= this->offset; X++ ) {
      RealType val = scale * exp ( - ( X * X ) / s );

      switch ( comp ) {
      case DIFF_X:
        val *= 2 * static_cast< RealType > ( X ) / s;
        break;
      case DIFF_XX:
        val *= 4 * aol::Sqr ( ( static_cast< RealType > ( X ) ) / s ) - 2 / s;
        break;
      default:
        cerr << comp << " invalid kernel type.\n";
        return;
      }
      this->setValue ( X, val );
  }

  // Diff kernels should have a mean value of zero.
  this->addToAll ( - this->getMeanValue() );
}


template class qc::Kernel1d<long double>;
template class qc::Kernel1d<double>;
template class qc::Kernel1d<float>;
template class qc::GaussDiffKernel1d<long double>;
template class qc::GaussDiffKernel1d<double>;
template class qc::GaussDiffKernel1d<float>;
template class qc::GaussKernel1d<long double>;
template class qc::GaussKernel1d<double>;
template class qc::GaussKernel1d<float>;

