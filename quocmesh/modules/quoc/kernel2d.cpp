#include <kernel2d.h>

template<typename RealType>
void qc::Kernel2d<RealType>::dump() const {
  for ( int Y = this->offset; Y >= -this->offset; Y-- ) {
    for ( int X = -this->offset; X <= this->offset; X++ ) {
      cerr << setw ( 6 ) << getValue ( X, Y ) << " ";
    }
    cerr << endl;
  }
}

template<typename RealType>
void qc::GaussKernel2d<RealType>::makeKernel() {
  for ( int X = -this->offset; X <= this->offset; X++ ) {
    for ( int Y = -this->offset; Y <= this->offset; Y++ ) {
      RealType val = exp ( -0.5f * ( X * X + Y * Y ) / ( sigma * sigma ) );
      this->setValue ( X, Y, val );
    }
  }
  this->normalize();
}


template<typename RealType>
void qc::GaussDiffKernel2d<RealType>::makeKernel() {

  const RealType scale = 1 / ( 2*aol::NumberTrait<RealType>::pi * aol::Sqr(sigma) );
  const RealType s = 2 * aol::Sqr(sigma);

  for ( int X = -this->offset; X <= this->offset; X++ ) {
    for ( int Y = -this->offset; Y <= this->offset; Y++ ) {
      RealType val = scale * exp ( - ( X * X + Y * Y ) / s );

      switch ( comp ) {
      case DIFF_X:
        val *= 2 * static_cast< RealType > ( X ) / s;
        break;
      case DIFF_Y:
        val *= 2 * static_cast< RealType > ( Y ) / s;
        break;
      case DIFF_XX:
        val *= 4 * aol::Sqr ( ( static_cast< RealType > ( X ) ) / s ) - 2 / s;
        break;
      case DIFF_YY:
        val *= 4 * aol::Sqr ( ( static_cast< RealType > ( Y ) ) / s ) - 2 / s;
        break;
      case DIFF_XY:
        val *= 4 * static_cast< RealType > ( Y )  * static_cast< RealType > ( X ) / aol::Sqr ( s );
        break;
      default:
        cerr << comp << " invalid kernel type.\n";
        return;
      }
      this->setValue ( X, Y, val );
    }
  }

  // Diff kernels should have a mean value of zero.
  this->addToAll ( - this->getMeanValue() );
}


template class qc::Kernel2d<long double>;
template class qc::Kernel2d<double>;
template class qc::Kernel2d<float>;
template class qc::GaussDiffKernel2d<long double>;
template class qc::GaussDiffKernel2d<double>;
template class qc::GaussDiffKernel2d<float>;
template class qc::GaussKernel2d<long double>;
template class qc::GaussKernel2d<double>;
template class qc::GaussKernel2d<float>;

template<typename RealType>
void qc::DiscreteGaussKernel2d<RealType>::makeKernel() {
  // First we calculate the line in Pascal's triangle that is as long as our desired filter width.
  aol::Vector<int> triangleLine ( 1 );
  aol::Vector<int> nextTriangleLine;
  triangleLine[0]=1;

  for ( int i = 1; i < this->size; ++i ) {
    nextTriangleLine.resize ( triangleLine.size() + 1 );
    nextTriangleLine[0] = 1;
    nextTriangleLine[ nextTriangleLine.size() - 1 ] = 1;
    for ( int j = 1; j < triangleLine.size(); ++j )
      nextTriangleLine[j] = triangleLine[j-1]+triangleLine[j];
    triangleLine.growBy(1);
    triangleLine = nextTriangleLine;
  }

  // The kernel is the tensor product of the triangle line with itself ...
  for ( int X = -this->offset; X <= this->offset; X++ ) {
    for ( int Y = -this->offset; Y <= this->offset; Y++ ) {
      this->setValue ( X, Y, triangleLine[X+this->offset]*triangleLine[Y+this->offset] );
    }
  }
  // ... normalized
  this->normalize();
}

template class qc::DiscreteGaussKernel2d<long double>;
template class qc::DiscreteGaussKernel2d<double>;
template class qc::DiscreteGaussKernel2d<float>;

template<typename RealType>
void qc::CircleAverageKernel2d<RealType>::makeKernel() {
  const int offsetSqr = aol::Sqr ( this->offset );
  for ( int X = -this->offset; X <= this->offset; X++ ) {
    for ( int Y = -this->offset; Y <= this->offset; Y++ ) {
      if ( aol::Sqr ( X ) + aol::Sqr ( Y ) <= offsetSqr )
        this->setValue ( X, Y, 1 );
    }
  }
  this->normalize();
}

template class qc::CircleAverageKernel2d<long double>;
template class qc::CircleAverageKernel2d<double>;
template class qc::CircleAverageKernel2d<float>;
