#ifndef __GEOM_H
#define __GEOM_H

#include <aol.h>
#include <quoc.h>
#include <triangle.h>
#include <smallMat.h>
#include <qmElement.h>

namespace aol {

enum PROJ_TYPE {
  PT_INSIDE, P1, P2
};

template <class PRIMITIVE_A, class PRIMITIVE_B>
class Intersector { };

template <class RealType> class Plane;
//template <class RealType> class Triang;
template <class RealType> class AlignedQuad;
template <class RealType> class Parallelogram;


template <typename RealType, qc::Dimension EmbedDim>
class LineSegment { };

template <typename RealType>
class LineSegment<RealType, qc::QC_2D> {
  friend class Intersector<Plane<RealType>, Triangle<RealType> >;
  friend class Intersector<AlignedQuad<RealType>, LineSegment<RealType, qc::QC_2D> >;
public:
  LineSegment(  ) : _pt1 ( aol::ZOTrait<RealType>::zero ), _pt2 ( aol::ZOTrait<RealType>::zero ) { }

  LineSegment ( const aol::Vec2<RealType> &Pt1, const aol::Vec2<RealType> &Pt2 ) :
      _pt1 ( Pt1 ), _pt2 ( Pt2 ) {}

  void getBoundingBox ( aol::Vec2<RealType> &min, aol::Vec2<RealType> &max ) const {
    for ( int i = 0; i < 2; i++ ) {
      min[i] = aol::Min ( _pt1[i], _pt2[i] );
      max[i] = aol::Max ( _pt1[i], _pt2[i] );
    }
  }


  void getBoundingBox ( aol::Vec2<int> &min, aol::Vec2<int> &max ) const {
    aol::Vec2<RealType> rmin, rmax;
    getBoundingBox ( rmin, rmax );
    for ( int i = 0; i < 2; i++ ) {
      min[i] = static_cast<int> ( floor ( rmin[i] ) );
      max[i] = static_cast<int> ( ceil ( rmax[i] ) );
    }
  }

  /** compute distance of a point (ignoring the third coordinate of the point) */
  RealType dist ( const aol::Vec3<RealType> &Pt ) const {
    return dist ( aol::Vec2<RealType> ( Pt[0], Pt[1] ) );
  }

  /** compute distance to a point */
  RealType dist ( const aol::Vec2<RealType> &Pt ) const;

  /** compute distance to linesegment */
  PROJ_TYPE projectTo ( const aol::Vec2<RealType> &Pt, aol::Vec2<RealType> &ProjPt ) const {
    aol::Vec2<RealType> a, b;
    a = _pt2; a -= Pt;
    b = _pt1; b -= _pt2;
    RealType lambda = - ( a * b ) / ( b * b );
    if ( lambda >= 0. && lambda <= 1. ) {
      ProjPt = b;
      ProjPt *= lambda;
      ProjPt += _pt2;
      return PT_INSIDE;
    } else if ( lambda < 0. ) {
      ProjPt = _pt2;
      return P2;
    } else {
      ProjPt = _pt1;
      return P2;
    }
  }

  void normal ( aol::Vec2<RealType> &n ) const {
    for ( int c = 0; c < 2; c++ ) {
      n[c] = _pt2[c] - _pt1[c];
    }
    n /= sqrt ( n * n );
    RealType tmp = -n[0];
    n[0] = n[1];
    n[1] = tmp;
  }

  RealType length( ) const {
    return sqrt ( aol::Sqr ( _pt1[0] - _pt2[0] ) +
                  aol::Sqr ( _pt1[1] - _pt2[1] ) );
  }

  aol::Vec2<RealType>& beginPoint() { return _pt1; }
  aol::Vec2<RealType>& endPoint() { return _pt2; }

  const aol::Vec2<RealType>& beginPoint() const { return _pt1; }
  const aol::Vec2<RealType>& endPoint() const { return _pt2; }
protected:
  aol::Vec2<RealType> _pt1, _pt2;
};


/**
 * utility class, representing a line segment, i.e. the connection between two points.
 * @ingroup Line
 */
template <typename RealType>
class LineSegment<RealType, qc::QC_3D> {
  friend class Intersector<Plane<RealType>, Triangle<RealType> >;
public:
  LineSegment(  ) {
    _pt1.setZero();
    _pt2.setZero();
  }

  LineSegment ( const aol::Vec3<RealType> &Pt1, const aol::Vec3<RealType> &Pt2 ) :
      _pt1 ( Pt1 ), _pt2 ( Pt2 ) {}

  /** compute distance to a point */
  RealType dist ( const aol::Vec3<RealType> &Pt ) const;

  /** compute distance to another line segment */
  RealType dist ( const aol::LineSegment<RealType, qc::QC_3D> &OtherLS ) const;

  /** compute distance to linesegment */
  void projectTo ( const aol::Vec3<RealType> &Pt, aol::Vec3<RealType> &ProjPt ) const;

  RealType length( ) const {
    return sqrt ( aol::Sqr ( _pt1[0] - _pt2[0] ) +
                  aol::Sqr ( _pt1[1] - _pt2[1] ) +
                  aol::Sqr ( _pt1[2] - _pt2[2] ) );
  }

  void getBoundingBox ( aol::Vec3<RealType> &BBoxMin, aol::Vec3<RealType> &BBoxMax ) const {
    BBoxMin = aol::CompWiseMin ( beginPoint(), endPoint() );
    BBoxMax = aol::CompWiseMax ( beginPoint(), endPoint() );
  }

  aol::Vec3<RealType>& beginPoint() { return _pt1; }
  aol::Vec3<RealType>& endPoint() { return _pt2; }

  const aol::Vec3<RealType>& beginPoint() const { return _pt1; }
  const aol::Vec3<RealType>& endPoint() const { return _pt2; }
protected:
  aol::Vec3<RealType> _pt1, _pt2;
};


template< typename RealType, qc::Dimension Dim >
ostream& operator<< ( ostream &Os, const aol::LineSegment<RealType, Dim> &LineSeg ) {
  Os << LineSeg.beginPoint() << " - " << LineSeg.endPoint();
  return ( Os );
}


template <class RealType>
class Plane {
public:
  aol::Vec3<RealType> _normal;
  RealType _alpha;

  Plane ( const aol::Vec3<RealType> &normal, RealType alpha )
      : _normal ( normal ), _alpha ( alpha ) {}
};


template <class RealType>  class AlignedCube {
public:
  aol::Vec3<RealType> _min, _max;

  AlignedCube ( const aol::Vec3<RealType> &min, const aol::Vec3<RealType> &max )
      : _min ( min ), _max ( max ) {}

  AlignedCube ( const RealType minX,
                const RealType minY,
                const RealType minZ,
                const RealType maxX,
                const RealType maxY,
                const RealType maxZ )
    : _min ( minX, minY, minZ ),
      _max ( maxX, maxY, maxZ ) {}
};

template <class RealType>  class AlignedQuad {
public:
  aol::Vec2<RealType> _min, _max;

  AlignedQuad ( const aol::Vec2<RealType> &min, const aol::Vec2<RealType> &max )
      : _min ( min ), _max ( max ) {}

};
  
template <class RealType>  class Parallelogram {
public:
  aol::Vec<2, RealType> _origin, _v1, _v2;
  
  Parallelogram ( const aol::Vec<2, RealType> &origin, const aol::Vec<2, RealType> &v1, const aol::Vec<2, RealType> &v2 )
    : _origin ( origin ), _v1 ( v1 ), _v2 ( v2 ) { }
  
  bool contains ( const aol::Vec<2, RealType> &p ) {
    const RealType d = _v1[0] * _v2[1] - _v1[1] * _v2[0];
    aol::Vec<2, RealType> pp ( p );
    pp -= _origin;
    const RealType det1 = -( pp[0] * _v1[1] - pp[1] * _v1[0] ) / d;
    const RealType det2 = ( pp[0] * _v2[1] - pp[1] * _v2[0] ) / d;
    return 0 <= det1 && det1 <= 1 && 0 <= det2 && det2 <= 1;
  }
};

/** 3D Cylinder segments between two points with radius
 *
 *  \author Schwen (MEVIS)
 */
template< typename RealType >
class CylinderSegment {
public:
  aol::LineSegment< RealType, qc::QC_3D > _lineSeg;
  RealType _radius;

  CylinderSegment ( const aol::Vec3<RealType> &PtA, const aol::Vec3<RealType> &PtB, const RealType Radius ) : _lineSeg ( PtA, PtB ), _radius ( Radius ) {
  }

  virtual ~CylinderSegment ( ) { /* do nothing */ }

  // upper bound on distance between two cylinder segments, currently this uses a simplified view interpreting cylinder segments as ending in a hemisphere with radius _radius
  RealType distWithSphericalCaps ( const aol::CylinderSegment<RealType> & Other ) const;

  aol::Vec3<RealType>& beginPoint() { return ( _lineSeg.beginPoint() ); }
  aol::Vec3<RealType>& endPoint() { return ( _lineSeg.endPoint() ); }

  const aol::Vec3<RealType>& beginPoint() const { return ( _lineSeg.beginPoint() ); }
  const aol::Vec3<RealType>& endPoint() const { return ( _lineSeg.endPoint() ); }

};

template< typename RealType >
ostream& operator<< ( ostream &Os, const aol::CylinderSegment<RealType> &CylSeg ) {
  Os << CylSeg.beginPoint() << " - " << CylSeg.endPoint() << ", " << aol::shortFormat ( CylSeg._radius );
  return ( Os );
}


template <class RealType>
class Intersector<Plane<RealType>, Triangle<RealType> > {
public:

  // returns number of intersections of t's edges with p.
  // see ref/geom1.tex
  static int cut ( const Plane<RealType> &p, const Triangle<RealType> &t, LineSegment<RealType, qc::QC_3D> &s ) {
    aol::Vec3<RealType> d1, d2;
    d1 = t[0]; d1 -= t[2];
    d2 = t[1]; d2 -= t[2];
    const RealType beta = p._alpha - p._normal * t[2];

    const RealType mu1 = d1 * p._normal, mu2 = d2 * p._normal;

    int num = 0;
    aol::Vec3<RealType> *pt = & ( s._pt1 );

    // lambda1 = 0
    RealType lambda2 = beta / mu2;
    if ( mu2 != 0. && lambda2 >= 0. && lambda2 <= 1. ) {
      *pt = t[1];
      *pt -= t[2];
      *pt *= lambda2;
      *pt += t[2];
      num ++;
      if ( num == 1 ) {
        pt = & ( s._pt2 );
      }
    }

    // lambda2 = 0
    RealType lambda1 = beta / mu1;
    if ( mu1 != 0. && lambda1 >= 0. && lambda1 <= 1. ) {
      *pt = t[0];
      *pt -= t[2];
      *pt *= lambda1;
      *pt += t[2];

      num++;
      if ( num == 1 ) {
        pt = & ( s._pt2 );
      }
    }

    // lambda3 = 0
    lambda2 = ( beta - mu1 ) / ( mu2 - mu1 );
    if ( lambda2 >= 0. && lambda2 <= 1. ) {
      *pt = t[1];
      *pt -= t[0];
      *pt *= lambda2;
      *pt += t[0];
      num ++;
    }

    if ( num > 2 ) {
      cerr << "error: num > 2!\n";
    }
    return num;
  }
};




template <class RealType>
class Intersector<LineSegment<RealType, qc::QC_2D>, LineSegment<RealType, qc::QC_2D> > {

public:
  static int cut ( const LineSegment<RealType, qc::QC_2D> &s1,
                   const LineSegment<RealType, qc::QC_2D> &s2, aol::Vec2<RealType> &pt  ) {

    aol::Matrix22<RealType> mat, inv;
    aol::Vec2<RealType> d1, d2, d3;
    aol::Vec2<RealType> rhs, alpha;

    d1 = s1.beginPoint(); d1 -= s1.endPoint();
    d2 = s2.endPoint(); d2 -= s2.beginPoint();
    d3 = s1.endPoint(); d3 -= s2.endPoint();

    mat[0][0] = d1 * d1;
    mat[1][0] = mat[0][1] = d1 * d2;
    mat[1][1] = d2 * d2;

    rhs[0] = -d1 * d3;
    rhs[1] = -d2 * d3;

    if ( mat.det() == 0. ) {
      return 0;
      // TODO: this could also be infinity, i.e. overlapping of the two segments
    }
    
    try {
      mat.make_inverse ( inv );
    } catch ( aol::Exception &ex ) {
      cerr << mat
      << " " << s1.beginPoint() << " " << s1.endPoint()
      << " " << s2.beginPoint() << " " << s2.endPoint() << endl;
      ex.dump();
      return 0;
    }
    inv.mult ( rhs, alpha );

    pt = s1.beginPoint();
    pt -= s1.endPoint();
    pt *= alpha[0];
    pt += s1.endPoint();

    if ( alpha[0] >= 0. && alpha[0] <= 1. &&
         alpha[1] >= 0. && alpha[1] <= 1. ) {
      return 1;
    } else {
      return 0;
    }

  }

};



template <class RealType>
class Intersector<LineSegment<RealType, qc::QC_3D>, LineSegment<RealType, qc::QC_3D> > {

public:
  static int cut ( const LineSegment<RealType, qc::QC_3D> &s1,
                   const LineSegment<RealType, qc::QC_3D> &s2, aol::Vec3<RealType> &pt  ) {

    aol::Matrix22<RealType> mat, inv;
    aol::Vec3<RealType> d1, d2, d3;
    aol::Vec2<RealType> rhs, alpha;

    d1 = s1.beginPoint(); d1 -= s1.endPoint();
    d2 = s2.endPoint(); d2 -= s2.beginPoint();
    d3 = s1.endPoint(); d3 -= s2.endPoint();

    mat[0][0] = d1 * d1;
    mat[1][0] = mat[0][1] = d1 * d2;
    mat[1][1] = d2 * d2;

    rhs[0] = -d1 * d3;
    rhs[1] = -d2 * d3;

    try {
      mat.make_inverse ( inv );
    } catch ( aol::Exception &ex ) {
      cerr << mat
      << " " << s1.beginPoint() << " " << s1.endPoint()
      << " " << s2.beginPoint() << " " << s2.endPoint() << endl;
      ex.dump();
      return 0;
    }
    inv.mult ( rhs, alpha );

    if ( mat.det() == 0. ) {
      return 0;
      // TODO: this could also be infinity, i.e. overlapping of the two segments
    }

    pt = s1.beginPoint();
    pt -= s1.endPoint();
    pt *= alpha[0];
    pt += s1.endPoint();

    if ( alpha[0] >= 0. && alpha[0] <= 1. &&
         alpha[1] >= 0. && alpha[1] <= 1. ) {
      return 1;
    } else {
      return 0;
    }

  }

};




template <class RealType>
class Intersector<AlignedCube<RealType>, Triangle<RealType> > {
public:

  static bool has_intersection ( const AlignedCube<RealType> &c,
                                 const Triangle<RealType> &t )  {
    aol::Vec3<RealType> bbmin, bbmax;
    t.getBoundingBox ( bbmin, bbmax );

    bool completely_inside = true;
    for ( int i = 0; i < 3; i++ ) {
      if ( ! ( c._min[i] <= bbmin[i] &&
               c._max[i] >= bbmax[i] ) ) {
        completely_inside = false;
        break;
      }
    }

    if ( completely_inside ) {
      //cerr << "completely inside: "
      //   << c._min << " " << c._max << " "
      //   << bbmin << " " << bbmax << " " << endl;

      return true;
    }


    if ( bbmin[0] <= c._min[0] && bbmax[0] >= c._min[0] ) {
      Plane<RealType> p ( aol::Vec3<RealType> ( 1., 0., 0. ), c._min[0] );
      aol::LineSegment<RealType, qc::QC_3D> seg;
      int num = Intersector<Plane<RealType>, Triangle<RealType> >::cut ( p, t, seg );

      aol::Vec2<RealType> min2d ( c._min[1], c._min[2] );
      aol::Vec2<RealType> max2d ( c._max[1], c._max[2] );
      aol::Vec2<RealType> seg2d1 ( seg.beginPoint() [1], seg.beginPoint() [2] );
      aol::Vec2<RealType> seg2d2 ( seg.endPoint() [1], seg.endPoint() [2] );

      if ( num == 2 &&
           ( is_inside ( min2d, max2d, seg2d1 ) ||
             is_inside ( min2d, max2d, seg2d2 ) ||
             rect_segment_intersect ( min2d, max2d, seg2d1, seg2d2 ) ) ) {
        return true;
      }
    }

    if ( bbmin[0] <= c._max[0] && bbmax[0] >= c._max[0] ) {
      Plane<RealType> p ( aol::Vec3<RealType> ( 1., 0., 0. ), c._max[0] );
      aol::LineSegment<RealType, qc::QC_3D> seg;
      int num = Intersector<Plane<RealType>, Triangle<RealType> >::cut ( p, t, seg );

      aol::Vec2<RealType> min2d ( c._min[1], c._min[2] );
      aol::Vec2<RealType> max2d ( c._max[1], c._max[2] );
      aol::Vec2<RealType> seg2d1 ( seg.beginPoint() [1], seg.beginPoint() [2] );
      aol::Vec2<RealType> seg2d2 ( seg.endPoint  () [1], seg.endPoint  () [2] );

      if ( num == 2 &&
           ( is_inside ( min2d, max2d, seg2d1 ) ||
             is_inside ( min2d, max2d, seg2d2 ) ||
             rect_segment_intersect ( min2d, max2d, seg2d1, seg2d2 ) ) ) {
        return true;
      }
    }


    if ( bbmin[1] <= c._min[1] && bbmax[1] >= c._min[1] ) {
      Plane<RealType> p ( aol::Vec3<RealType> ( 0., 1., 0. ), c._min[1] );
      aol::LineSegment<RealType, qc::QC_3D> seg;
      int num = Intersector<Plane<RealType>, Triangle<RealType> >::cut ( p, t, seg );

      aol::Vec2<RealType> min2d ( c._min[0], c._min[2] );
      aol::Vec2<RealType> max2d ( c._max[0], c._max[2] );
      aol::Vec2<RealType> seg2d1 ( seg.beginPoint() [0], seg.beginPoint() [2] );
      aol::Vec2<RealType> seg2d2 ( seg.endPoint() [0], seg.endPoint() [2] );

      if ( num == 2 &&
           ( is_inside ( min2d, max2d, seg2d1 ) ||
             is_inside ( min2d, max2d, seg2d2 ) ||
             rect_segment_intersect ( min2d, max2d, seg2d1, seg2d2 ) ) ) {
        return true;
      }
    }


    if ( bbmin[1] <= c._max[1] && bbmax[1] >= c._max[1] ) {
      Plane<RealType> p ( aol::Vec3<RealType> ( 0., 1., 0. ), c._max[1] );
      aol::LineSegment<RealType, qc::QC_3D> seg;
      int num = Intersector<Plane<RealType>, Triangle<RealType> >::cut ( p, t, seg );

      aol::Vec2<RealType> min2d ( c._min[0], c._min[2] );
      aol::Vec2<RealType> max2d ( c._max[0], c._max[2] );
      aol::Vec2<RealType> seg2d1 ( seg.beginPoint() [0], seg.beginPoint() [2] );
      aol::Vec2<RealType> seg2d2 ( seg.endPoint() [0], seg.endPoint() [2] );

      if ( num == 2 &&
           ( is_inside ( min2d, max2d, seg2d1 ) ||
             is_inside ( min2d, max2d, seg2d2 ) ||
             rect_segment_intersect ( min2d, max2d, seg2d1, seg2d2 ) ) ) {
        return true;
      }
    }


    if ( bbmin[2] <= c._min[2] && bbmax[2] >= c._min[2] ) {
      Plane<RealType> p ( aol::Vec3<RealType> ( 0., 0., 1. ), c._min[2] );
      aol::LineSegment<RealType, qc::QC_3D> seg;
      int num = Intersector<Plane<RealType>, Triangle<RealType> >::cut ( p, t, seg );

      aol::Vec2<RealType> min2d ( c._min[0], c._min[1] );
      aol::Vec2<RealType> max2d ( c._max[0], c._max[1] );
      aol::Vec2<RealType> seg2d1 ( seg.beginPoint() [0], seg.beginPoint() [1] );
      aol::Vec2<RealType> seg2d2 ( seg.endPoint() [0], seg.endPoint() [1] );

      if ( num == 2 &&
           ( is_inside ( min2d, max2d, seg2d1 ) ||
             is_inside ( min2d, max2d, seg2d2 ) ||
             rect_segment_intersect ( min2d, max2d, seg2d1, seg2d2 ) ) ) {
        return true;
      }
    }

    if ( bbmin[2] <= c._max[2] && bbmax[2] >= c._max[2] ) {
      Plane<RealType> p ( aol::Vec3<RealType> ( 0., 0., 1. ), c._max[2] );
      aol::LineSegment<RealType, qc::QC_3D> seg;
      int num = Intersector<Plane<RealType>, Triangle<RealType> >::cut ( p, t, seg );

      aol::Vec2<RealType> min2d ( c._min[0], c._min[1] );
      aol::Vec2<RealType> max2d ( c._max[0], c._max[1] );
      aol::Vec2<RealType> seg2d1 ( seg.beginPoint() [0], seg.beginPoint() [1] );
      aol::Vec2<RealType> seg2d2 ( seg.endPoint() [0], seg.endPoint() [1] );

      if ( num == 2 &&
           ( is_inside ( min2d, max2d, seg2d1 ) ||
             is_inside ( min2d, max2d, seg2d2 ) ||
             rect_segment_intersect ( min2d, max2d, seg2d1, seg2d2 ) ) ) {
        return true;
      }
    }
    return false;
  }

  inline static bool is_inside ( const aol::Vec2<RealType> &min, const aol::Vec2<RealType> &max,
                                 const aol::Vec2<RealType> &p ) {
    return ( min[0] <= p[0] && min[1] <= p[1] && p[0] <= max[0] && p[1] <= max[1] );
  }

  inline static bool rect_segment_intersect ( const aol::Vec2<RealType> &min, const aol::Vec2<RealType> &max,
                                              const aol::Vec2<RealType> &p1, const aol::Vec2<RealType> &p2 ) {
    return ( segments_intersect ( p1, p2, aol::Vec2<RealType> ( min[0], min[1] ),  aol::Vec2<RealType> ( min[0], max[1] ) ) ||
             segments_intersect ( p1, p2, aol::Vec2<RealType> ( min[0], min[1] ),  aol::Vec2<RealType> ( max[0], min[1] ) ) ||
             segments_intersect ( p1, p2, aol::Vec2<RealType> ( min[0], max[1] ),  aol::Vec2<RealType> ( max[0], max[1] ) ) ||
             segments_intersect ( p1, p2, aol::Vec2<RealType> ( max[0], min[1] ),  aol::Vec2<RealType> ( max[0], max[1] ) ) );
  }

  inline static bool segments_intersect ( const aol::Vec2<RealType> &p00, const aol::Vec2<RealType> &p01,
                                          const aol::Vec2<RealType> &p10, const aol::Vec2<RealType> &p11 ) {
    aol::Vec2<RealType> alpha;
    return segments_intersect ( p00, p01, p10, p11, alpha );
  }


  inline static bool segments_intersect ( const aol::Vec2<RealType> &p00, const aol::Vec2<RealType> &p01,
                                          const aol::Vec2<RealType> &p10, const aol::Vec2<RealType> &p11,
                                          aol::Vec2<RealType> &Alpha ) {
    aol::Matrix22<RealType> mat, inv;
    aol::Vec2<RealType> rhs;

    mat[0] = p00;
    mat[0] -= p01;

    mat[1] = p11;
    mat[1] -= p10;

    mat.transpose();

    rhs = p11;
    rhs -= p01;

    try {
      mat.make_inverse ( inv );
    } catch ( aol::Exception &ex ) {
      ex.consume();
      aol::Vec2<RealType> v, w;
      v  = p01;
      v -= p00;

      w  = p10;
      w -= p00;


      RealType alpha = ( v * w ) / ( w * w );

      v *= alpha;
      v += p00;

      if ( alpha >= 0. && alpha <= 1. && euclidianDist ( v, p10 ) < 1e-10 ) {
        return true;
      }

      w = p11;
      w -= p00;


      alpha = ( v * w ) / ( w * w );

      v *= alpha;
      v += p00;


      if ( alpha >= 0. && alpha <= 1. && euclidianDist ( v, p11 ) < 1e-10 ) {
        return true;
      }

      return false;
    }

    inv.mult ( rhs, Alpha );

    if ( Alpha[0] >= 0. && Alpha[0] <= 1. &&
         Alpha[1] >= 0. && Alpha[1] <= 1. ) {
      return true;
    } else {
      return false;
    }
  }
};



template <class RealType>
class Intersector<AlignedQuad<RealType>, LineSegment<RealType, qc::QC_2D> > {
public:

  static bool has_intersection ( const AlignedQuad<RealType> &c,
                                 const LineSegment<RealType, qc::QC_2D> &s )  {

    aol::Vec2<RealType> bbmin, bbmax, cut;
    s.getBoundingBox ( bbmin, bbmax );

    bool completely_inside = true;
    for ( int i = 0; i < 2; i++ ) {
      if ( ! ( c._min[i] <= bbmin[i] &&
               c._max[i] >= bbmax[i] ) ) {
        completely_inside = false;
        break;
      }
    }

    if ( completely_inside ) {
      return true;
    }

    if ( bbmin[0] <= c._max[0] && bbmax[0] >= c._max[0] ) {
      aol::Vec2<RealType> a ( c._max[0], c._min[1] );
      aol::Vec2<RealType> b ( c._max[0], c._max[1] );
      aol::LineSegment<RealType, qc::QC_2D> qs ( a, b );
      if ( Intersector <
           aol::LineSegment<RealType, qc::QC_2D>,
           aol::LineSegment<RealType, qc::QC_2D> >:: cut ( qs, s, cut ) ) {
        return true;
      }
    }

    if ( bbmin[0] <= c._min[0] && bbmax[0] >= c._min[0] ) {
      aol::Vec2<RealType> a ( c._min[0], c._min[1] );
      aol::Vec2<RealType> b ( c._min[0], c._max[1] );
      aol::LineSegment<RealType, qc::QC_2D> qs ( a, b );
      if ( Intersector <
           aol::LineSegment<RealType, qc::QC_2D>,
           aol::LineSegment<RealType, qc::QC_2D> >:: cut ( qs, s, cut ) ) {
        return true;
      }
    }


    if ( bbmin[1] <= c._max[1] && bbmax[1] >= c._max[1] ) {
      aol::Vec2<RealType> a ( c._min[0], c._max[1] );
      aol::Vec2<RealType> b ( c._max[0], c._max[1] );
      aol::LineSegment<RealType, qc::QC_2D> qs ( a, b );
      if ( Intersector <
           aol::LineSegment<RealType, qc::QC_2D>,
           aol::LineSegment<RealType, qc::QC_2D> >:: cut ( qs, s, cut ) ) {
        return true;
      }
    }

    if ( bbmin[1] <= c._min[1] && bbmax[1] >= c._min[1] ) {
      aol::Vec2<RealType> a ( c._min[0], c._min[1] );
      aol::Vec2<RealType> b ( c._max[0], c._min[1] );
      aol::LineSegment<RealType, qc::QC_2D> qs ( a, b );
      if ( Intersector <
           aol::LineSegment<RealType, qc::QC_2D>,
           aol::LineSegment<RealType, qc::QC_2D> >:: cut ( qs, s, cut ) ) {
        return true;
      }
    }
    return false;
  }
  
  // Sets the LineSegment intersec as the intersection between c and s, if there is an intersection (otherwise, sets begin and end point to infinitiy)
  // Returns 0 in case of no intersections, 1 if one end point of s is in c and the other outside, otherwise 2
  // If intersection is only a point, this is represented as a LineSegment with equal begin and end point
  static int cut ( const AlignedQuad<RealType> &c,
                   const LineSegment<RealType, qc::QC_2D> &s,
                   LineSegment<RealType, qc::QC_2D> &intersec ) {
    MultiVector<RealType> pts;
    
    aol::Vec2<RealType> bbmin, bbmax, cut;
    s.getBoundingBox ( bbmin, bbmax );
    
    // Check if the begin or end point of the LineSegment (or both) are inside the AlignedQuad (if so, add them as end points of intersecting line segment)
    if ( s.beginPoint ( )[0] > c._min[0] && s.beginPoint ( )[0] < c._max[0]
      && s.beginPoint ( )[1] > c._min[1] && s.beginPoint ( )[1] < c._max[1] ) {
      pts.resize ( pts.numComponents ( ) + 1, 2 );
      pts[pts.numComponents ( )-1][0] = s.beginPoint( )[0]; pts[pts.numComponents ( )-1][1] = s.beginPoint( )[1];
    }
    if ( s.endPoint ( )[0] > c._min[0] && s.endPoint ( )[0] < c._max[0]
        && s.endPoint ( )[1] > c._min[1] && s.endPoint ( )[1] < c._max[1] ) {
      pts.resize ( pts.numComponents ( ) + 1, 2 );
      pts[pts.numComponents ( )-1][0] = s.endPoint( )[0]; pts[pts.numComponents ( )-1][1] = s.endPoint( )[1];
    }
    
    // Cut with left side of AlignedQuad
    if ( bbmin[0] <= c._min[0] && bbmax[0] >= c._min[0] ) {
      aol::Vec2<RealType> a ( c._min[0], c._min[1] );
      aol::Vec2<RealType> b ( c._min[0], c._max[1] );
      aol::LineSegment<RealType, qc::QC_2D> qs ( a, b );
      if ( Intersector <
          aol::LineSegment<RealType, qc::QC_2D>,
          aol::LineSegment<RealType, qc::QC_2D> >:: cut ( qs, s, cut ) ) {
        pts.resize ( pts.numComponents ( ) + 1, 2 );
        pts[pts.numComponents ( )-1][0] = cut[0]; pts[pts.numComponents ( )-1][1] = cut[1];
      }
    }
    
    // Cut with right side of AlignedQuad
    if ( bbmin[0] <= c._max[0] && bbmax[0] >= c._max[0] ) {
      aol::Vec2<RealType> a ( c._max[0], c._min[1] );
      aol::Vec2<RealType> b ( c._max[0], c._max[1] );
      aol::LineSegment<RealType, qc::QC_2D> qs ( a, b );
      if ( Intersector <
          aol::LineSegment<RealType, qc::QC_2D>,
          aol::LineSegment<RealType, qc::QC_2D> >:: cut ( qs, s, cut ) ) {
        pts.resize ( pts.numComponents ( ) + 1, 2 );
        pts[pts.numComponents ( )-1][0] = cut[0]; pts[pts.numComponents ( )-1][1] = cut[1];
      }
    }
    
    // Cut with bottom side of AlignedQuad
    if ( bbmin[1] <= c._max[1] && bbmax[1] >= c._max[1] ) {
      aol::Vec2<RealType> a ( c._min[0], c._max[1] );
      aol::Vec2<RealType> b ( c._max[0], c._max[1] );
      aol::LineSegment<RealType, qc::QC_2D> qs ( a, b );
      if ( Intersector <
          aol::LineSegment<RealType, qc::QC_2D>,
          aol::LineSegment<RealType, qc::QC_2D> >:: cut ( qs, s, cut ) ) {
        pts.resize ( pts.numComponents ( ) + 1, 2 );
        pts[pts.numComponents ( )-1][0] = cut[0]; pts[pts.numComponents ( )-1][1] = cut[1];
      }
    }
    
    // Cut with top side of AlignedQuad
    if ( bbmin[1] <= c._min[1] && bbmax[1] >= c._min[1] ) {
      aol::Vec2<RealType> a ( c._min[0], c._min[1] );
      aol::Vec2<RealType> b ( c._max[0], c._min[1] );
      aol::LineSegment<RealType, qc::QC_2D> qs ( a, b );
      if ( Intersector <
          aol::LineSegment<RealType, qc::QC_2D>,
          aol::LineSegment<RealType, qc::QC_2D> >:: cut ( qs, s, cut ) ) {
        pts.resize ( pts.numComponents ( ) + 1, 2 );
        pts[pts.numComponents ( )-1][0] = cut[0]; pts[pts.numComponents ( )-1][1] = cut[1];
      }
    }
    
    // Remove possible duplicate cuts in corners
    for ( int i=0; i<pts.numComponents ( )-1 ; ++i )
      for ( int j=i+1; j<pts.numComponents ( ) ; ++j )
        if ( pts[i][0] == pts[j][0] && pts[i][1] == pts[j][1] )
          pts.eraseComponent ( j );
    
    if ( pts.numComponents ( ) >= 1 ) {
      intersec._pt1[0] = pts[0][0];
      intersec._pt1[1] = pts[0][1];
      
      if ( pts.numComponents ( ) == 2 ) {
        intersec._pt2[0] = pts[1][0];
        intersec._pt2[1] = pts[1][1];
      } else {
        intersec._pt2[0] = pts[0][0];
        intersec._pt2[1] = pts[0][1];
      }
    }
    
    return pts.numComponents ( );
  }
  
  static RealType cutLength ( const AlignedQuad<RealType> &c,
                              const LineSegment<RealType, qc::QC_2D> &s ) {
    LineSegment<RealType, qc::QC_2D> intersec;
    cut ( c, s, intersec );
    return intersec.length ( );
  }
};


template <class RealType>
class Intersector<AlignedCube<RealType>, LineSegment<RealType, qc::QC_3D> > {
public:

  /**
   * Heavily optimized code, taken from http://www.3dkingdoms.com/weekly/bbox.cpp,
   * CBBox::IsLineInBox( const CVec3& L1, const CVec3& L2 )
   *
   * \author Berkels
   */
  static bool has_intersection ( const AlignedCube<RealType> &Cube,
                                 const LineSegment<RealType, qc::QC_3D> &Segment )  {
    // Translate everything so that the center of the cube is the origin.
    const aol::Vec3<RealType> cubeExtent = (Cube._max - Cube._min) / 2.0f;
    const aol::Vec3<RealType> cubeMid = (Cube._min + Cube._max) * 0.5f;
    const aol::Vec3<RealType> translatesSegmentBegin = Segment.beginPoint() - cubeMid;
    const aol::Vec3<RealType> translatesSegmentEnd = Segment.endPoint() - cubeMid;

    // Get line midpoint and extent
    const aol::Vec3<RealType> sMid = (translatesSegmentBegin + translatesSegmentEnd) * 0.5f;
    const aol::Vec3<RealType> s = (translatesSegmentBegin - sMid);
    const aol::Vec3<RealType> sExt = aol::Vec3<RealType>( aol::Abs(s[0]), aol::Abs(s[1]), aol::Abs(s[2]) );

    // Use Separating Axis Test
    // Separation vector from box center to line center is sMid, since the line is in box space
    if ( aol::Abs( sMid[0] ) > cubeExtent[0] + sExt[0] ) return false;
    if ( aol::Abs( sMid[1] ) > cubeExtent[1] + sExt[1] ) return false;
    if ( aol::Abs( sMid[2] ) > cubeExtent[2] + sExt[2] ) return false;
    // Crossproducts of line and each axis
    if ( aol::Abs( sMid[1] * s[2] - sMid[2] * s[1])  >  (cubeExtent[1] * sExt[2] + cubeExtent[2] * sExt[1]) ) return false;
    if ( aol::Abs( sMid[0] * s[2] - sMid[2] * s[0])  >  (cubeExtent[0] * sExt[2] + cubeExtent[2] * sExt[0]) ) return false;
    if ( aol::Abs( sMid[0] * s[1] - sMid[1] * s[0])  >  (cubeExtent[0] * sExt[1] + cubeExtent[1] * sExt[0]) ) return false;
    // No separating axis, the line intersects
    return true;
  }

  /**
   * Interprets the element as cube given by (Cube[0], Cube[1], Cube[2]) and (Cube[0] + 1, Cube[1] + 1, Cube[2] + 1)
   *
   * \author Berkels
   */
  static bool has_intersection ( const qc::Element &Cube,
                                 const LineSegment<RealType, qc::QC_3D> &Segment )  {
    aol::AlignedCube<RealType> cube ( Cube[0], Cube[1], Cube[2],
                                      Cube[0] + 1, Cube[1] + 1, Cube[2] + 1 );
    return has_intersection ( cube, Segment );
  }

  static bool inline GetIntersection( RealType fDst1, RealType fDst2, aol::Vec3<RealType> P1, aol::Vec3<RealType> P2, aol::Vec3<RealType> &Hit) {
    if ( (fDst1 * fDst2) >= 0.0f) return 0;
    if ( fDst1 == fDst2) return 0;
    Hit = P1 + (P2-P1) * ( -fDst1/(fDst2-fDst1) );
    return 1;
  }

  static bool inline InBox( aol::Vec3<RealType> Hit, aol::Vec3<RealType> B1, aol::Vec3<RealType> B2, const int Axis) {
    if ( Axis==1 && Hit[2] > B1[2] && Hit[2] < B2[2] && Hit[1] > B1[1] && Hit[1] < B2[1]) return 1;
    if ( Axis==2 && Hit[2] > B1[2] && Hit[2] < B2[2] && Hit[0] > B1[0] && Hit[0] < B2[0]) return 1;
    if ( Axis==3 && Hit[0] > B1[0] && Hit[0] < B2[0] && Hit[1] > B1[1] && Hit[1] < B2[1]) return 1;
    return 0;
  }

  /**
   * Returns true if line (L1, L2) intersects with the box (B1, B2).
   * If the line intersects with the box, an intersection point is returned
   * in Hit. In case the line is completely inside the box, Hit is set to L1,
   * otherwise it is a point on the boundary of the box.
   *
   * Code from http://www.3dkingdoms.com/weekly/weekly.php?a=3.
   *
   * \author Berkels
   *
   * \note Should only be used if one actually needs the intersection point as
   *       the version of has_intersection that doesn't calculate the point is
   *       considerably faster.
   */
  static bool has_intersection ( const AlignedCube<RealType> &Cube,
                                 const LineSegment<RealType, qc::QC_3D> &Segment,
                                 aol::Vec3<RealType> &Hit )  {
    const aol::Vec3<RealType> &B1 = Cube._min;
    const aol::Vec3<RealType> &B2 = Cube._max;
    const aol::Vec3<RealType> &L1 = Segment.beginPoint();
    const aol::Vec3<RealType> &L2 = Segment.endPoint();

    for ( int i = 0; i < 3; ++ i ) {
      if (L2[i] < B1[i] && L1[i] < B1[i]) return false;
      if (L2[i] > B2[i] && L1[i] > B2[i]) return false;
    }

    // Is the line completely inside the box?
    if ( L1[0] > B1[0] && L1[0] < B2[0] &&
         L1[1] > B1[1] && L1[1] < B2[1] &&
         L1[2] > B1[2] && L1[2] < B2[2]) {
      Hit = L1;
      return true;
    }

    if ( (GetIntersection( L1[0]-B1[0], L2[0]-B1[0], L1, L2, Hit) && InBox( Hit, B1, B2, 1 ))
      || (GetIntersection( L1[1]-B1[1], L2[1]-B1[1], L1, L2, Hit) && InBox( Hit, B1, B2, 2 ))
      || (GetIntersection( L1[2]-B1[2], L2[2]-B1[2], L1, L2, Hit) && InBox( Hit, B1, B2, 3 ))
      || (GetIntersection( L1[0]-B2[0], L2[0]-B2[0], L1, L2, Hit) && InBox( Hit, B1, B2, 1 ))
      || (GetIntersection( L1[1]-B2[1], L2[1]-B2[1], L1, L2, Hit) && InBox( Hit, B1, B2, 2 ))
      || (GetIntersection( L1[2]-B2[2], L2[2]-B2[2], L1, L2, Hit) && InBox( Hit, B1, B2, 3 )))
      return true;

    return false;
  }

};


//! Conversion from Cartesian to homogeneous coordinates (copy entries and write 1 to last component)
template< typename RealType, int Dim >
inline void CartesianToHomogeneousCoordinates ( const aol::Vec<Dim,RealType> &arg, aol::Vec<Dim+1,RealType> &dest ) {
  for ( int d = 0; d < Dim; ++d ) {
    dest[d] = arg[d];
  }
  dest[Dim] = aol::NumberTrait<RealType>::one;
}

//! Conversion from homogeneous to Cartesian coordinates (division by last component)
template< typename RealType, int Dim >
inline void HomogeneousToCartesianCoordinates ( const aol::Vec<Dim+1,RealType> &arg, aol::Vec<Dim,RealType> &dest ) {
  for ( int d = 0; d < Dim; ++d ) {
    dest[d] = arg[d] / arg[Dim];
  }
}

//! perform a transformation given for homoegeneous coordinates to vectors in Cartesian coordinates
template< typename RealType, int Dim >
void TransformCartesianCoordinatesByHomogeneousMapping ( const aol::Vec<Dim, RealType> &arg, const aol::Mat<Dim+1, Dim+1, RealType> &transMat, aol::Vec<Dim, RealType> &dest ) {
  aol::Vec<Dim+1,RealType> argH, destH;
  aol::CartesianToHomogeneousCoordinates ( arg, argH );
  transMat.mult ( argH, destH );
  aol::HomogeneousToCartesianCoordinates ( destH, dest );
}

//! compute angle (in radians) between two Vecs
template< int Dim, typename RealType >
RealType AngleBetween ( const aol::Vec<Dim, RealType> &A, const aol::Vec<Dim, RealType> &B ) {
  return ( acos ( ( A * B ) / ( A.norm() * B.norm() ) ) );
}
  
//! check if a point lies within a line segment / rectangle / cuboid
template< qc::Dimension Dim, typename IndexType>
bool InsideQuad ( const aol::Vec<Dim, IndexType> &Lower, const aol::Vec<Dim, IndexType> &Upper, const aol::Vec<Dim, IndexType> &X, const bool ExclusiveUpper = true ) {
  for ( int d=0; d<Dim ; ++d ) {
    if ( X[d] < Lower[d] || X[d] > Upper[d] ) return false;
    if ( ExclusiveUpper && X[d] == Upper[d] ) return false;
  }
  return true;
}
  
template< qc::Dimension Dim>
bool InsideQuad ( const qc::CoordType &Lower, const qc::CoordType &Upper, const qc::CoordType &X, const bool ExclusiveUpper = true, const bool IgnoreExcessiveComponents = false ) {
  if ( !IgnoreExcessiveComponents ) {
    for ( int d=Dim ; d<3 ; ++d )
      if ( X[d] != 0 ) return false;
  }
  for ( int d=0; d<Dim ; ++d ) {
    if ( X[d] < Lower[d] || X[d] > Upper[d] ) return false;
    if ( ExclusiveUpper && X[d] == Upper[d] ) return false;
  }
  return true;
}


} // end namespace aol

#endif
