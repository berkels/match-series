#ifndef __SIMPLEXBASEFUNCTIONSET_H
#define __SIMPLEXBASEFUNCTIONSET_H

#include <quoc.h>
#include <baseFunctionSet.h>
#include <smallVec.h>
#include <simplexLookup.h>

namespace qc {

namespace simplex {

//=============================================================================================================================

template <typename RealType, Dimension Dim>
class MidpointQuadrature {};

//-----------------------------------------------------------------------------------------------------------------------------

//! midpoint quadrature on triangle.
//! \author von Deylen
template <typename RealType>
class MidpointQuadrature<RealType, qc::QC_2D> {
public:
  typedef aol::Vec2<RealType> VecType;
  typedef aol::BarCoord<qc::QC_2D, RealType> DomVecType;

  enum { numQuadPoints = 1 };
  inline static const DomVecType & getRefCoord ( int /*QuadPoint*/ );
  inline static RealType getWeight ( int /*QuadPoint*/ );

protected:
  static const DomVecType _barycentricPoint;
};

//-----------------------------------------------------------------------------------------------------------------------------

// definition of static members
template <typename RealType>
const aol::BarCoord<qc::QC_2D, RealType>
MidpointQuadrature<RealType, qc::QC_2D>::_barycentricPoint ( 1./3., 1./3., 1./3. );

// class functions implementation
template <typename RealType>
inline const aol::BarCoord<qc::QC_2D, RealType> &
MidpointQuadrature<RealType, qc::QC_2D>::getRefCoord ( int /*QuadPoint*/ ) {
  return _barycentricPoint;
}

template <typename RealType>
inline RealType
MidpointQuadrature<RealType, qc::QC_2D>::getWeight ( int /*QuadPoint*/ ) {
  return 1.;
}

//=============================================================================================================================

//! midpoint quadrature on tetrahedra
//! \author von Deylen
template <typename _RealType>
class MidpointQuadrature<_RealType, qc::QC_3D> {
public:
  typedef _RealType RealType;
  static const qc::Dimension Dim=qc::QC_3D;
  static const int Order = 1;
  typedef aol::Vec3<RealType> VecType;
  typedef aol::BarCoord<qc::QC_3D, RealType> DomVecType;

  enum { numQuadPoints = 1 };
  inline static const DomVecType & getRefCoord ( int /*QuadPoint*/ );
  inline static RealType getWeight ( int /*QuadPoint*/ );

protected:
  static const DomVecType _barycentricPoint;
};

//-----------------------------------------------------------------------------------------------------------------------------

// definition of static members
template <typename RealType>
const aol::BarCoord<qc::QC_3D, RealType>
MidpointQuadrature<RealType, qc::QC_3D>::_barycentricPoint ( 0.25, 0.25, 0.25, 0.25 );

// class functions implementation
template <typename RealType>
inline const aol::BarCoord<qc::QC_3D, RealType> &
MidpointQuadrature<RealType, qc::QC_3D>::getRefCoord ( int /*QuadPoint*/ ) {
  return _barycentricPoint;
}

template <typename RealType>
inline RealType
MidpointQuadrature<RealType, qc::QC_3D>::getWeight ( int /*QuadPoint*/ ) {
  return 1.;
}

//=============================================================================================================================

/**
 * Small helper function to create barycentric coordinates in 2D from 2 values.
 *
 * \author Effland
 */
template <typename RealType>
aol::BarCoord<qc::QC_2D, RealType> createBarycentricCoord ( const RealType A, const RealType B ) {
  return aol::BarCoord<qc::QC_2D, RealType> ( A, B, aol::ZOTrait<RealType>::one - A - B );
}

/**
 * Small helper function to create barycentric coordinates in 3D from 3 values.
 *
 * \author Berkels
 */
template <typename RealType>
aol::BarCoord<qc::QC_3D, RealType> createBarycentricCoord ( const RealType A, const RealType B, const RealType C ) {
  return aol::BarCoord<qc::QC_3D, RealType> ( A, B, C, aol::ZOTrait<RealType>::one - A - B - C );
}

/**
 * Second order quadrature on tetrahedra.
 *
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
class SecondOrderQuadrature {};

template <typename RealType>
class SecondOrderQuadrature<RealType, qc::QC_3D> {
public:
  enum { numQuadPoints = 8 };

  inline static const aol::BarCoord<qc::QC_3D, RealType> & getRefCoord ( int QuadPoint ) {
    return _points[ QuadPoint ];
  }
  inline static RealType getWeight ( int QuadPoint ) {
    return _weights[QuadPoint];
  }
protected:
  static aol::BarCoord<qc::QC_3D, RealType> _points [ numQuadPoints ];
  static RealType _weights[ numQuadPoints ];
};

template <typename RealType>
RealType SecondOrderQuadrature<RealType, qc::QC_3D>::_weights[] = { 6 * 0.0091694,
                                                                    6 * 0.0160270,
                                                                    6 * 0.0211570,
                                                                    6 * 0.0369799,
                                                                    6 * 0.0091694,
                                                                    6 * 0.0160270,
                                                                    6 * 0.0211570,
                                                                    6 * 0.0369799 };

template <typename RealType>
aol::BarCoord<qc::QC_3D, RealType> SecondOrderQuadrature<RealType, qc::QC_3D>::_points[] = {
  createBarycentricCoord ( 0.5441518, 0.2939988, 0.0342028 ),
  createBarycentricCoord ( 0.5441518, 0.0706797, 0.0813957 ),
  createBarycentricCoord ( 0.1225148, 0.5659332, 0.0658387 ),
  createBarycentricCoord ( 0.1225148, 0.1360550, 0.1566826 ),
  createBarycentricCoord ( 0.5441518, 0.2939988, 0.1276466 ),
  createBarycentricCoord ( 0.5441518, 0.0706797, 0.3037728 ),
  createBarycentricCoord ( 0.1225148, 0.5659332, 0.2457133 ),
  createBarycentricCoord ( 0.1225148, 0.1360550, 0.5847476 ) };

/**
 * Second order quadrature on triangles.
 *
 * \author Franken
 */

template <typename RealType>
class SecondOrderQuadrature<RealType, qc::QC_2D> {
public:
  enum { numQuadPoints = 3 };

  inline static const aol::BarCoord<qc::QC_2D, RealType> & getRefCoord ( int QuadPoint ) {
    return _points[ QuadPoint ];
  }
  inline static RealType getWeight ( int QuadPoint ) {
    return _weights[QuadPoint];
  }
protected:
  static aol::BarCoord<qc::QC_2D, RealType> _points [ numQuadPoints ];
  static RealType _weights[ numQuadPoints ];
};

template <typename RealType>
RealType SecondOrderQuadrature<RealType, qc::QC_2D>::_weights[] = { 1./3.,
                                                                    1./3.,
                                                                    1./3. };

template <typename RealType>
aol::BarCoord<qc::QC_2D, RealType> SecondOrderQuadrature<RealType, qc::QC_2D>::_points[] = {
  aol::BarCoord<qc::QC_2D, RealType> ( 0.5, 0.5, 0.0 ),
  aol::BarCoord<qc::QC_2D, RealType> ( 0.0, 0.5, 0.5 ),
  aol::BarCoord<qc::QC_2D, RealType> ( 0.5, 0.0, 0.5 ) };

//=============================================================================================================================

/**
 * Second order quadrature on tetrahedra
 * (see "MODERATE DEGREE CUBATURE FORMULAS FOR 3-D TETRAHEDRAL FINITE-ELEMENT APPROXIMATIONS", Gellert and Harbord)
 *
 * \author Effland
 */

template <typename RealType, qc::Dimension Dim>
class SecondOrderQuadratureImproved {};

template <typename RealType>
class SecondOrderQuadratureImproved<RealType, qc::QC_3D> {
public:
  enum { numQuadPoints = 4 };

  inline static const aol::BarCoord<qc::QC_3D, RealType> & getRefCoord ( int QuadPoint ) {
    return _points[ QuadPoint ];
  }
  inline static RealType getWeight ( int QuadPoint ) {
    return _weights[QuadPoint];
  }
protected:
  static aol::BarCoord<qc::QC_3D, RealType> _points [ numQuadPoints ];
  static RealType _weights[ numQuadPoints ];
};

template <typename RealType>
RealType SecondOrderQuadratureImproved<RealType, qc::QC_3D>::_weights[] = { 0.25,
                                                                            0.25,
                                                                            0.25,
                                                                            0.25 };

template <typename RealType>
aol::BarCoord<qc::QC_3D, RealType> SecondOrderQuadratureImproved<RealType, qc::QC_3D>::_points[] = {
  createBarycentricCoord ( 0.5854101966249685, 0.1381966011250105, 0.1381966011250105 ),
  createBarycentricCoord ( 0.1381966011250105, 0.5854101966249685, 0.1381966011250105 ),
  createBarycentricCoord ( 0.1381966011250105, 0.1381966011250105, 0.5854101966249685 ),
  createBarycentricCoord ( 0.1381966011250105, 0.1381966011250105, 0.1381966011250105 ) };

/**
 * Third order quadrature on tetrahedra
 * (see "MODERATE DEGREE CUBATURE FORMULAS FOR 3-D TETRAHEDRAL FINITE-ELEMENT APPROXIMATIONS", Gellert and Harbord)
 *
 * \author Effland
 */

template <typename RealType, qc::Dimension Dim>
class ThirdOrderQuadrature {};

template <typename RealType>
class ThirdOrderQuadrature<RealType, qc::QC_3D> {
public:
  enum { numQuadPoints = 5 };

  inline static const aol::BarCoord<qc::QC_3D, RealType> & getRefCoord ( int QuadPoint ) {
    return _points[ QuadPoint ];
  }
  inline static RealType getWeight ( int QuadPoint ) {
    return _weights[QuadPoint];
  }
protected:
  static aol::BarCoord<qc::QC_3D, RealType> _points [ numQuadPoints ];
  static RealType _weights[ numQuadPoints ];
};

template <typename RealType>
RealType ThirdOrderQuadrature<RealType, qc::QC_3D>::_weights[] = { -0.8,
                                                                    0.45,
                                                                    0.45,
                                                                    0.45,
                                                                    0.45 };

template <typename RealType>
aol::BarCoord<qc::QC_3D, RealType> ThirdOrderQuadrature<RealType, qc::QC_3D>::_points[] = {
  createBarycentricCoord ( 0.25, 0.25, 0.25 ),
  createBarycentricCoord ( 0.5, 0.1666666666666667, 0.1666666666666667 ),
  createBarycentricCoord ( 0.1666666666666667, 0.5, 0.1666666666666667 ),
  createBarycentricCoord ( 0.1666666666666667, 0.1666666666666667, 0.5 ),
  createBarycentricCoord ( 0.1666666666666667, 0.1666666666666667, 0.1666666666666667 ) };

/**
 * Fourth order quadrature on triangle
 * (see "Several new quadrature formulas for polynomial integration in the triangle", Taylor, Wingate and Bos)
 *
 * \author Effland
 */

template <typename RealType, qc::Dimension Dim>
class FourthOrderQuadrature {};

template <typename RealType>
class FourthOrderQuadrature<RealType, qc::QC_2D> {
public:
  enum { numQuadPoints = 6 };

  inline static const aol::BarCoord<qc::QC_2D, RealType> & getRefCoord ( int QuadPoint ) {
    return _points[ QuadPoint ];
  }
  inline static RealType getWeight ( int QuadPoint ) {
    return _weights[QuadPoint];
  }
protected:
  static aol::BarCoord<qc::QC_2D, RealType> _points [ numQuadPoints ];
  static RealType _weights[ numQuadPoints ];
};

template <typename RealType>
RealType FourthOrderQuadrature<RealType, qc::QC_2D>::_weights[] = { 0.10995174365,
                                                                    0.10995174365,
                                                                    0.10995174365,
                                                                    0.22338158967,
                                                                    0.22338158967,
                                                                    0.22338158967 };

template <typename RealType>
aol::BarCoord<qc::QC_2D, RealType> FourthOrderQuadrature<RealType, qc::QC_2D>::_points[] = {
  createBarycentricCoord ( 0.0915762135098, 0.0915762135098 ),
  createBarycentricCoord ( 0.8168475729805, 0.0915762135098 ),
  createBarycentricCoord ( 0.0915762135098, 0.8168475729805 ),
  createBarycentricCoord ( 0.1081030181681, 0.4459484909160 ),
  createBarycentricCoord ( 0.4459484909160, 0.1081030181681 ),
  createBarycentricCoord ( 0.4459484909160, 0.4459484909160 ) };


/**
 * Fourth order quadrature on tetrahedra
 * (see "MODERATE DEGREE CUBATURE FORMULAS FOR 3-D TETRAHEDRAL FINITE-ELEMENT APPROXIMATIONS", Gellert and Harbord)
 *
 * \author Effland
 */

template <typename RealType>
class FourthOrderQuadrature<RealType, qc::QC_3D> {
public:
  enum { numQuadPoints = 11 };

  inline static const aol::BarCoord<qc::QC_3D, RealType> & getRefCoord ( int QuadPoint ) {
    return _points[ QuadPoint ];
  }
  inline static RealType getWeight ( int QuadPoint ) {
    return _weights[QuadPoint];
  }
protected:
  static aol::BarCoord<qc::QC_3D, RealType> _points [ numQuadPoints ];
  static RealType _weights[ numQuadPoints ];
};

template <typename RealType>
RealType FourthOrderQuadrature<RealType, qc::QC_3D>::_weights[] = { -0.07893333333333333,
                                                                     0.04573333333333333,
                                                                     0.04573333333333333,
                                                                     0.04573333333333333,
                                                                     0.04573333333333333,
                                                                     0.14933333333333333,
                                                                     0.14933333333333333,
                                                                     0.14933333333333333,
                                                                     0.14933333333333333,
                                                                     0.14933333333333333,
                                                                     0.14933333333333333 };

template <typename RealType>
aol::BarCoord<qc::QC_3D, RealType> FourthOrderQuadrature<RealType, qc::QC_3D>::_points[] = {
  createBarycentricCoord ( 0.25, 0.25, 0.25 ),
  createBarycentricCoord ( 0.7857142857142857, 0.0714285714285714, 0.0714285714285714 ),
  createBarycentricCoord ( 0.0714285714285714, 0.7857142857142857, 0.0714285714285714 ),
  createBarycentricCoord ( 0.0714285714285714, 0.0714285714285714, 0.7857142857142857 ),
  createBarycentricCoord ( 0.0714285714285714, 0.0714285714285714, 0.0714285714285714 ),
  createBarycentricCoord ( 0.3994035761667992, 0.3994035761667992, 0.1005964238332008 ),
  createBarycentricCoord ( 0.3994035761667992, 0.1005964238332008, 0.3994035761667992 ),
  createBarycentricCoord ( 0.3994035761667992, 0.1005964238332008, 0.1005964238332008 ),
  createBarycentricCoord ( 0.1005964238332008, 0.3994035761667992, 0.3994035761667992 ),
  createBarycentricCoord ( 0.1005964238332008, 0.3994035761667992, 0.1005964238332008 ),
  createBarycentricCoord ( 0.1005964238332008, 0.1005964238332008, 0.3994035761667992 ) };

/**
 * Fifth order quadrature on tetrahedra
 * (see "MODERATE DEGREE CUBATURE FORMULAS FOR 3-D TETRAHEDRAL FINITE-ELEMENT APPROXIMATIONS", Gellert and Harbord)
 *
 * \author Effland
 */

template <typename RealType, qc::Dimension Dim>
class FifthOrderQuadrature {};

template <typename RealType>
class FifthOrderQuadrature<RealType, qc::QC_3D> {
public:
  enum { numQuadPoints = 14 };

  inline static const aol::BarCoord<qc::QC_3D, RealType> & getRefCoord ( int QuadPoint ) {
    return _points[ QuadPoint ];
  }
  inline static RealType getWeight ( int QuadPoint ) {
    return _weights[QuadPoint];
  }
protected:
  static aol::BarCoord<qc::QC_3D, RealType> _points [ numQuadPoints ];
  static RealType _weights[ numQuadPoints ];
};

template <typename RealType>
RealType FifthOrderQuadrature<RealType, qc::QC_3D>::_weights[] = { 0.1126879257180159,
                                                                   0.1126879257180159,
                                                                   0.1126879257180159,
                                                                   0.1126879257180159,
                                                                   0.07349304311636195,
                                                                   0.07349304311636195,
                                                                   0.07349304311636195,
                                                                   0.07349304311636195,
                                                                   0.04254602077708147,
                                                                   0.04254602077708147,
                                                                   0.04254602077708147,
                                                                   0.04254602077708147,
                                                                   0.04254602077708147,
                                                                   0.04254602077708147 };

template <typename RealType>
aol::BarCoord<qc::QC_3D, RealType> FifthOrderQuadrature<RealType, qc::QC_3D>::_points[] = {
  createBarycentricCoord ( 0.06734224221009817, 0.3108859192633006, 0.3108859192633006 ),
  createBarycentricCoord ( 0.3108859192633006, 0.06734224221009817, 0.3108859192633006 ),
  createBarycentricCoord ( 0.3108859192633006, 0.3108859192633006, 0.06734224221009817 ),
  createBarycentricCoord ( 0.3108859192633006, 0.3108859192633006, 0.3108859192633006 ),
  createBarycentricCoord ( 0.7217942490673263, 0.09273525031089123, 0.09273525031089123 ),
  createBarycentricCoord ( 0.09273525031089123, 0.7217942490673263, 0.09273525031089123 ),
  createBarycentricCoord ( 0.09273525031089123, 0.09273525031089123, 0.7217942490673263 ),
  createBarycentricCoord ( 0.09273525031089123, 0.09273525031089123, 0.09273525031089123 ),
  createBarycentricCoord ( 0.4544962958743504, 0.4544962958743504, 0.04550370412564965 ),
  createBarycentricCoord ( 0.4544962958743504, 0.04550370412564965, 0.4544962958743504 ),
  createBarycentricCoord ( 0.4544962958743504, 0.04550370412564965, 0.04550370412564965 ),
  createBarycentricCoord ( 0.04550370412564965, 0.4544962958743504, 0.4544962958743504 ),
  createBarycentricCoord ( 0.04550370412564965, 0.4544962958743504, 0.04550370412564965 ),
  createBarycentricCoord ( 0.04550370412564965, 0.04550370412564965, 0.4544962958743504 ) };

/**
 * Seventh order quadrature on triangle
 * (see "Several new quadrature formulas for polynomial integration in the triangle", Taylor, Wingate and Bos)
 *
 * \author Effland
 */

template <typename RealType, qc::Dimension Dim>
class SeventhOrderQuadrature {};

template <typename RealType>
class SeventhOrderQuadrature<RealType, qc::QC_2D> {
public:
  enum { numQuadPoints = 15 };

  inline static const aol::BarCoord<qc::QC_2D, RealType> & getRefCoord ( int QuadPoint ) {
    return _points[ QuadPoint ];
  }
  inline static RealType getWeight ( int QuadPoint ) {
    return _weights[QuadPoint];
  }
protected:
  static aol::BarCoord<qc::QC_2D, RealType> _points [ numQuadPoints ];
  static RealType _weights[ numQuadPoints ];
};

template <typename RealType>
RealType SeventhOrderQuadrature<RealType, qc::QC_2D>::_weights[] = { 0.0051279087046,
                                                                    0.0051279087046,
                                                                    0.0051279087046,
                                                                    0.05580235233235,
                                                                    0.05580235233235,
                                                                    0.05580235233235,
                                                                    0.05580235233235,
                                                                    0.05580235233235,
                                                                    0.05580235233235,
                                                                    0.08398877976675,
                                                                    0.08398877976675,
                                                                    0.08398877976675,
                                                                    0.1326119401973,
                                                                    0.1326119401973,
                                                                    0.1326119401973 };

template <typename RealType>
aol::BarCoord<qc::QC_2D, RealType> SeventhOrderQuadrature<RealType, qc::QC_2D>::_points[] = {
  createBarycentricCoord ( aol::ZOTrait<RealType>::one, aol::ZTrait<RealType>::zero ),
  createBarycentricCoord ( aol::ZTrait<RealType>::zero, aol::ZTrait<RealType>::zero ),
  createBarycentricCoord ( aol::ZTrait<RealType>::zero, aol::ZOTrait<RealType>::one ),
  createBarycentricCoord ( 0.7839656651012, 0.0421382841642 ),
  createBarycentricCoord ( 0.1738960507345, 0.7839656651012 ),
  createBarycentricCoord ( 0.1738960507345, 0.0421382841642 ),
  createBarycentricCoord ( 0.0421382841642, 0.1738960507345 ),
  createBarycentricCoord ( 0.7839656651012, 0.1738960507345 ),
  createBarycentricCoord ( 0.0421382841642, 0.7839656651012 ),
  createBarycentricCoord ( 0.4743880861752, 0.4743880861752 ),
  createBarycentricCoord ( 0.4743880861752, 0.0512238276497 ),
  createBarycentricCoord ( 0.0512238276497, 0.4743880861752 ),
  createBarycentricCoord ( 0.2385615300181, 0.5228769399639 ),
  createBarycentricCoord ( 0.5228769399639, 0.2385615300181 ),
  createBarycentricCoord ( 0.2385615300181, 0.2385615300181 ) };


template < typename RealType, qc::Dimension Dim, unsigned short int QuadChoice >
struct SimplexCubatureTrait {
};

template < typename RealType >
struct SimplexCubatureTrait < RealType, qc::QC_2D, 1 > {
  typedef MidpointQuadrature < RealType, qc::QC_2D > QuadType;
};

template < typename RealType >
struct SimplexCubatureTrait < RealType, qc::QC_2D, 2 > {
  typedef SecondOrderQuadrature < RealType, qc::QC_2D > QuadType;
};

template < typename RealType >
struct SimplexCubatureTrait < RealType, qc::QC_2D, 4 > {
  typedef FourthOrderQuadrature < RealType, qc::QC_2D > QuadType;
};

template < typename RealType >
struct SimplexCubatureTrait < RealType, qc::QC_2D, 7 > {
  typedef SeventhOrderQuadrature < RealType, qc::QC_2D > QuadType;
};

template < typename RealType >
struct SimplexCubatureTrait < RealType, qc::QC_3D, 1 > {
  typedef MidpointQuadrature < RealType, qc::QC_3D > QuadType;
};

template < typename RealType >
struct SimplexCubatureTrait < RealType, qc::QC_3D, 2 > {
  typedef SecondOrderQuadratureImproved < RealType, qc::QC_3D > QuadType;
};

template < typename RealType >
struct SimplexCubatureTrait < RealType, qc::QC_3D, 3 > {
  typedef ThirdOrderQuadrature < RealType, qc::QC_3D > QuadType;
};

template < typename RealType >
struct SimplexCubatureTrait < RealType, qc::QC_3D, 4 > {
  typedef FourthOrderQuadrature < RealType, qc::QC_3D > QuadType;
};

template < typename RealType >
struct SimplexCubatureTrait < RealType, qc::QC_3D, 5 > {
  typedef FifthOrderQuadrature < RealType, qc::QC_3D > QuadType;
};

//=============================================================================================================================

template <typename RealType, qc::Dimension Dim, typename QuadRuleType>
class BaseFunctionSetLin {};

//-----------------------------------------------------------------------------------------------------------------------------

//! linear basis functions on triangles
//! \author von Deylen (july 2008)
template <typename RealType, class QuadRuleType>
class BaseFunctionSetLin<RealType, QC_2D, QuadRuleType> :
public aol::BaseFunctionSetInterface<RealType,
                     typename aol::VecDimTrait<RealType, QC_2D>::VecType,
                     typename aol::VecDimTrait<RealType, QC_2D>::BarCoordType,
                     3, QuadRuleType,
                     BaseFunctionSetLin<RealType, QC_2D, QuadRuleType> > {
public:
  typedef aol::BaseFunctionSetInterface<RealType,
                     typename aol::VecDimTrait<RealType, QC_2D>::VecType,
                     typename aol::VecDimTrait<RealType, QC_2D>::BarCoordType,
                     3, QuadRuleType,
                     BaseFunctionSetLin<RealType, QC_2D, QuadRuleType> > Base;
  typedef typename aol::VecDimTrait<RealType, QC_2D>::BarCoordType       DomVecType;
  typedef typename aol::VecDimTrait<RealType, QC_2D>::VecType            VecType;

  BaseFunctionSetLin ( RealType H, int SimplexNumber )
  : _h ( H )
  , _simplexNumber ( SimplexNumber ) {
    this->initializeQuadCache();
  }

  using Base::evaluate;
  using Base::evaluateGradient;

  inline RealType evaluate ( int BaseFuncNum, const DomVecType & RefCoord ) const {
    return RefCoord[BaseFuncNum];
  }
  inline void evaluateGradient ( int BaseFuncNum, const DomVecType &, VecType & grad ) const {
    grad = Lookup<RealType, qc::QC_2D>::gradients[_simplexNumber][BaseFuncNum];
    grad /= _h;
  }

protected:
  RealType _h;
  short _simplexNumber;
};

//-----------------------------------------------------------------------------------------------------------------------------

//! linear basis functions on tetrahedra
//! \author von Deylen (july 2008)
template <typename RealType, class QuadRuleType>
class BaseFunctionSetLin<RealType, QC_3D, QuadRuleType> :
  public aol::BaseFunctionSetInterface<RealType,
                     typename aol::VecDimTrait<RealType, QC_3D>::VecType,
                     typename aol::VecDimTrait<RealType, QC_3D>::BarCoordType,
                     4, QuadRuleType,
                     BaseFunctionSetLin<RealType, QC_3D, QuadRuleType> > {
public:
  typedef aol::BaseFunctionSetInterface<RealType,
                     typename aol::VecDimTrait<RealType, QC_3D>::VecType,
                     typename aol::VecDimTrait<RealType, QC_3D>::BarCoordType,
                     4, QuadRuleType,
                     BaseFunctionSetLin<RealType, QC_3D, QuadRuleType> > Base;
  typedef typename aol::VecDimTrait<RealType, QC_3D>::BarCoordType       DomVecType;
  typedef typename aol::VecDimTrait<RealType, QC_3D>::VecType            VecType;

  BaseFunctionSetLin ( RealType H, int SimplexNumber )
  : _h ( H )
  , _simplexNumber ( SimplexNumber ) {
    this->initializeQuadCache();
  }

  using Base::evaluate;
  using Base::evaluateGradient;

  inline RealType evaluate ( int BaseFuncNum, const DomVecType & RefCoord ) const {
    return RefCoord[BaseFuncNum];
  }
  inline void evaluateGradient ( int BaseFuncNum, const DomVecType &, VecType & grad ) const {
    grad = Lookup<RealType, qc::QC_3D>::gradients[_simplexNumber][BaseFuncNum];
    grad /= _h;
  }

protected:
  RealType _h;
  short _simplexNumber;
};

//-----------------------------------------------------------------------------------------------------------------------------

} // end of namespace simplex.

} // end of namespace qc.

#endif
