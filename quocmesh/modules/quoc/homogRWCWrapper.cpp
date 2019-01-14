#include <homogRWCWrapper.h>

namespace qc {

template< typename RealType >
HomogRWCMapper<RealType, qc::QC_3D>::HomogRWCMapper ( const aol::Vec3<RealType> &Offset, const aol::Vec3<RealType> &Resolutions, const qc::GridSize<qc::QC_3D> &GridSize ) : _gridSize ( GridSize ) {
  // set matrices
  for ( short i = 0; i < 3; ++i ) {
    _toInternalConvertMat[i][i] = 1.0 / Resolutions[i];
    _toInternalConvertMat[i][3] = - Offset[i] / Resolutions[i];
  }
  _toInternalConvertMat[3][3] = 1.0;

  for ( short i = 0; i < 3; ++i ) {
    _toWorldConvertMat[i][i] = Resolutions[i];
    _toWorldConvertMat[i][3] = Offset[i];
  }
  _toWorldConvertMat[3][3] = 1.0;
}


template< typename RealType >
void HomogRWCMapper<RealType, qc::QC_3D>::toInternal ( const aol::Vec3<RealType> &Arg, aol::Vec3<RealType> &Dest ) const {
  aol::TransformCartesianCoordinatesByHomogeneousMapping ( Arg, _toInternalConvertMat, Dest );
}

template< typename RealType >
void HomogRWCMapper<RealType, qc::QC_3D>::toRWC ( const aol::Vec3<RealType> &Arg, aol::Vec3<RealType> &Dest ) const {
  aol::TransformCartesianCoordinatesByHomogeneousMapping ( Arg, _toWorldConvertMat, Dest );
}


template< typename RealType >
RealType HomogRWCMapper<RealType, qc::QC_3D>::getVoxelVolume ( ) const {
  return ( _toWorldConvertMat.det() );
}


template class HomogRWCMapper<float, qc::QC_3D>;
template class HomogRWCMapper<double, qc::QC_3D>;
template class HomogRWCMapper<long double, qc::QC_3D>;

}
