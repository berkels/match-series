#include <geom.h>

namespace aol {

template < typename RealType >
RealType aol::LineSegment<RealType, qc::QC_2D>::dist ( const aol::Vec2<RealType> &Pt ) const {
  aol::Vec2<RealType> a, b;
  a = _pt2; a -= Pt;
  b = _pt1; b -= _pt2;
  RealType lambda = - ( a * b ) / ( b * b );
  if ( lambda >= 0. && lambda <= 1. ) {
    b *= lambda; a += b;
    return a.norm();
  } else {
    RealType d1, d2;
    a = _pt2; a -= Pt;
    d1 = a.norm();
    a = _pt1; a -= Pt;
    d2 = a.norm();
    return aol::Min ( d1, d2 );
  }
}


template < typename RealType >
RealType aol::LineSegment<RealType, qc::QC_3D>::dist ( const aol::Vec3<RealType> &Pt ) const {
  aol::Vec3<RealType> a, b;
  a = _pt2; a -= Pt;
  b = _pt1; b -= _pt2;
  RealType lambda = - ( a * b ) / ( b * b );
  if ( lambda >= 0. && lambda <= 1. ) {
    b *= lambda; a += b;
    return a.norm();
  } else {
    RealType d1, d2;
    a = _pt2; a -= Pt;
    d1 = a.norm();
    a = _pt1; a -= Pt;
    d2 = a.norm();
    return aol::Min ( d1, d2 );
  }
}


template < typename RealType >
RealType aol::LineSegment<RealType, qc::QC_3D>::dist ( const aol::LineSegment<RealType, qc::QC_3D> &OtherLS ) const {
  const aol::Vec3<RealType> dir0 = ( _pt2 - _pt1 ), dir1 = ( OtherLS._pt2 - OtherLS._pt1 );
  if ( aol::appeqAbsolute ( ( dir0.crossProduct ( dir1 ) ).normSqr(), aol::ZTrait<RealType>::zero ) ) { // parallel case
    return ( aol::Min ( this->dist ( OtherLS._pt1 ), this->dist ( OtherLS._pt2 ) ) );
  } else {
    // points on lines are parametrized as below with lm and mu
    // find minimum of Const + B lm^2 + C mu^2 + D lm + E mu + f lm mu as root of gradient
    // this can probably be done more efficiently ...

    const RealType
      B = dir0.normSqr(),
      C = dir1.normSqr(),
      D = 2 * ( ( _pt1 - OtherLS._pt1 ) * dir0 ),
      E = -2 * ( ( _pt1 - OtherLS._pt1 ) * dir1 ),
      F = -2 * ( dir0 * dir1 );

    const aol::Matrix22<RealType> sysMat ( 2 * B, F, F, 2 * C );
    const aol::Vec2<RealType> sysRhs ( - D, - E );
    const aol::Matrix22<RealType> sysMatInv = sysMat.inverse();
    const aol::Vec2<RealType> sysSoln = sysMatInv * sysRhs;

    const RealType
      lm = sysSoln[0],
      mu = sysSoln[1];

    if ( ( ( lm < 0  ) && ( mu < 0  ) ) || ( ( lm < 0  ) && ( mu > 1  ) ) ||
         ( ( lm > 1  ) && ( mu < 0  ) ) || ( ( lm > 1  ) && ( mu > 1  ) ) ) {
      // note that lm < 0 and mu < 0 does not imply that the first one is the minimum!
      return ( aol::Min ( aol::Min ( ( _pt1 - OtherLS._pt1 ).norm(), ( _pt1 - OtherLS._pt2 ).norm()), aol::Min ( ( _pt2 - OtherLS._pt1 ).norm(), ( _pt2 - OtherLS._pt2 ).norm()) ) );
    } else if ( ( lm < 0  ) && ( mu <= 1 ) ) {
      return ( OtherLS.dist ( this->_pt1 ) );
    } else if ( ( lm <= 1 ) && ( mu < 0  ) ) {
      return ( this->dist ( OtherLS._pt1 ) );
    } else if ( ( lm <= 1 ) && ( mu <= 1 ) ) {
      return ( ( ( _pt1 + lm * dir0 ) - ( OtherLS._pt1 + mu * dir1 ) ).norm() );
    } else if ( ( lm <= 1 ) && ( mu > 1  ) ) {
      return ( this->dist ( OtherLS._pt2 ) );
    } else if ( ( lm > 1  ) && ( mu <= 1 ) ) {
      return ( OtherLS.dist ( this->_pt2 ) );
    } else {
      throw aol::Exception ( "This cannot happen" );
    }
  }
}


template < typename RealType >
void aol::LineSegment<RealType, qc::QC_3D>::projectTo ( const aol::Vec3<RealType> &Pt, aol::Vec3<RealType> &ProjPt ) const {
  aol::Vec3<RealType> a, b;
  a = _pt2; a -= Pt;
  b = _pt1; b -= _pt2;
  RealType lambda = - ( a * b ) / ( b * b );
  if ( lambda >= 0. && lambda <= 1. ) {
    ProjPt = b;
    ProjPt *= lambda;
    ProjPt += _pt2;
  } else if ( lambda < 0. ) {
    ProjPt = _pt2;
  } else {
    ProjPt = _pt1;
  }
}


template class aol::LineSegment<float, qc::QC_2D>;
template class aol::LineSegment<double, qc::QC_2D>;
template class aol::LineSegment<long double, qc::QC_2D>;

template class aol::LineSegment<float, qc::QC_3D>;
template class aol::LineSegment<double, qc::QC_3D>;
template class aol::LineSegment<long double, qc::QC_3D>;


template< typename RealType >
RealType CylinderSegment<RealType>::distWithSphericalCaps ( const aol::CylinderSegment<RealType> & Other ) const {
  RealType ctrDist;
  // treat degenerate cases that are currently not treated correctly in the aol::LineSegment

  if ( ( this->beginPoint() - this->endPoint() ).normSqr() < 1e-20 ) {
    ctrDist = Other._lineSeg.dist ( this->beginPoint() );
  } else if ( ( Other.beginPoint() - Other.endPoint() ).normSqr() < 1e-20 ) {
    ctrDist = this->_lineSeg.dist ( Other.beginPoint() );
  } else {
    ctrDist = this->_lineSeg.dist ( Other._lineSeg );
  }
  return ( aol::Max ( ctrDist - this->_radius - Other._radius, aol::NumberTrait<RealType>::zero ) );
}


template class aol::CylinderSegment<float>;
template class aol::CylinderSegment<double>;
template class aol::CylinderSegment<long double>;

}
