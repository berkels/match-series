#include "micrographQuantifier.h"

namespace im {
  
  
template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void MicrographQuantifier<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
                          _ColoredPictureType>::getInterAtomicDistances ( aol::MultiVector<RealType> &Distances, const aol::MultiVector<RealType> &Centers,
                                                                          const RealType PeriodX, const RealType PeriodY, const RealType PeriodDelta,
                                                                          const RealType AngleX, const RealType AngleY, const RealType AngleDelta,
                                                                          const char *GnuplotXMatches, const char *GnuplotYMatches ) const {
  Distances.reallocate ( 6, 0 );
  
  // Iterate over all atoms, find corresponding neighbors (atoms within the specified x-, y-distance and angle)
  // and collect the corresponding distances separately for x and y
  std::ofstream outX, outY;
  if ( GnuplotXMatches )
    outX.open ( GnuplotXMatches );
  if ( GnuplotYMatches )
    outY.open ( GnuplotYMatches );

  for ( int i=0; i<Centers.numComponents ( ) ; ++i ) {
    short numXNeighbors = 0, numYNeighbors = 0;
    aol::MultiVector<RealType> centersXNeighbors ( Centers.numComponents ( ), 2 ), centersYNeighbors ( Centers.numComponents ( ), 2 );
    for ( int j=0; j<Centers.numComponents ( ) ; ++j ) {
      if ( i != j ) {
        const RealType dx = Centers[j][0] - Centers[i][0];
        const RealType dy = Centers[j][1] - Centers[i][1];
        const RealType angle = atan2 ( dy, dx ) * 180 / aol::NumberTrait<RealType>::pi;
        const RealType dist = aol::Vec2<RealType> ( dx, dy ).norm ( );
        if ( dx > 0 && aol::Abs<RealType> ( dist - PeriodX ) < PeriodDelta && aol::Abs<RealType> ( angle - AngleX ) < AngleDelta ) {
          Distances[0].pushBack ( dist );
          Distances[2].pushBack ( dx );
          Distances[4].pushBack ( dy );
          centersXNeighbors[numXNeighbors][0] = Centers[j][0];
          centersXNeighbors[numXNeighbors][1] = Centers[j][1];
          ++numXNeighbors;
        }
        if ( dy > 0 && aol::Abs<RealType> ( dist - PeriodY ) < PeriodDelta && aol::Abs<RealType> ( angle - AngleY ) < AngleDelta ) {
          Distances[1].pushBack ( dist );
          Distances[3].pushBack ( dx );
          Distances[5].pushBack ( dy );
          centersYNeighbors[numYNeighbors][0] = Centers[j][0];
          centersYNeighbors[numYNeighbors][1] = Centers[j][1];
          ++numYNeighbors;
        }
      }
    }
    if ( numXNeighbors > 1 && !_quietMode ) {
      std::cerr << "Warning: found multiple neighbors in x-direction for atom at (" << Centers[i][0] << "," << Centers[i][1] << ")!" << std::endl;
      std::cerr << "Centers of the neighbors:" << std::endl;
      for ( int k=0; k<numXNeighbors ; ++k ) std::cerr << centersXNeighbors[k] << std::endl;
      std::cerr << std::endl;
    }
    if ( numYNeighbors > 1 && !_quietMode ) {
      std::cerr << "Warning: found multiple neighbors in y-direction for atom at (" << Centers[i][0] << "," << Centers[i][1] << ")!" << std::endl;
      std::cerr << "Centers of the neighbors:" << std::endl;
      for ( int k=0; k<numYNeighbors ; ++k ) std::cerr << centersYNeighbors[k] << std::endl;
      std::cerr << std::endl;
    }
    if ( GnuplotXMatches ) {
      for ( int k = 0; k < numXNeighbors; ++k ) {
        outX << Centers[i][0] << " " << Centers[i][1] << " ";
        outX << centersXNeighbors[k][0]-Centers[i][0] << " " << centersXNeighbors[k][1]-Centers[i][1] << endl;
      }
    }
    if ( GnuplotYMatches ) {
      for ( int k = 0; k < numYNeighbors; ++k ) {
        outY << Centers[i][0] << " " << Centers[i][1] << " ";
        outY << centersYNeighbors[k][0]-Centers[i][0] << " " << centersYNeighbors[k][1]-Centers[i][1] << endl;
      }
    }
  }
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
_RealType MicrographQuantifier<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
                               _ColoredPictureType>::getFidelity ( const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate ) const {
  aol::Vector<short> correspondences;
  getCorrespondences ( correspondences, CentersRef, CentersEstimate );
  RealType res = 0;
  short N = 0;
  for ( short i=0; i<correspondences.size ( ) ; ++i ) {
    if ( correspondences[i] >= 0 ) {
      res += aol::Vec2<RealType> ( CentersRef[i][0] - CentersEstimate[correspondences[i]][0], CentersRef[i][1] - CentersEstimate[correspondences[i]][1] ).normSqr ( );
      ++N;
    }
  }
  res = ( N > 0 ) ? res / N : aol::NumberTrait<RealType>::Inf;
  return sqrt ( res );
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void MicrographQuantifier<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
                          _ColoredPictureType>::getCorrespondences ( aol::Vector<short> &Correspondences,
                                                                     const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate ) const {
  Correspondences.reallocate ( CentersRef.numComponents ( ) );
  RealType dist, minInterAtomicDist = aol::NumberTrait<RealType>::Inf;
                                   
  for ( short i=0; i<CentersRef.numComponents ( ) ; ++i )
    for ( short j=0; j<CentersRef.numComponents ( ) ; ++j ) {
      if ( i != j ) {
        dist = aol::Vec2<RealType> ( CentersRef[i][0] - CentersRef[j][0], CentersRef[i][1] - CentersRef[j][1] ).norm ( );
        if ( dist < minInterAtomicDist )
          minInterAtomicDist = dist;
    }
  }
                                   
  RealType minDist;
  short correspondence;
  std::set<short> matchedCentersEstimate;
  for ( short i=0; i<CentersRef.numComponents ( ) ; ++i ) {
    correspondence = -1;
    minDist = aol::NumberTrait<RealType>::Inf;
    for ( short j=0; j<CentersEstimate.numComponents ( ) ; ++j ) {
      if ( matchedCentersEstimate.find ( j ) == matchedCentersEstimate.end ( ) ) {
        dist = aol::Vec2<RealType> ( CentersRef[i][0] - CentersEstimate[j][0], CentersRef[i][1] - CentersEstimate[j][1] ).norm ( );
        if ( dist < minDist ) {
          minDist = dist;
          correspondence = j;
        }
      }
    }
    if ( minDist < 0.5 * minInterAtomicDist ) {
      Correspondences[i] = correspondence;
      matchedCentersEstimate.insert ( correspondence );
    } else Correspondences[i] = -1;
  }
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
int MicrographQuantifier<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
                         _ColoredPictureType>::getNumCorrespondences ( const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate ) const {
  aol::Vector<short> correspondences;
  getCorrespondences ( correspondences, CentersRef, CentersEstimate );
  short N = 0;
  for ( short i=0; i<correspondences.size ( ) ; ++i )
    if ( correspondences[i] >= 0 ) ++N;
  return N;
}


template class MicrographQuantifier<double, aol::FullMatrix<double>, aol::LinearRegressionQR<double>,
                                    qc::ScalarArray<double, qc::QC_2D>, qc::MultiArray<double, qc::QC_2D, 3> >;

} // namespace im
