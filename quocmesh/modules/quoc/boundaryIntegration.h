#ifndef __BOUNDARYINTEGRATION_H
#define __BOUNDARYINTEGRATION_H

#include <aol.h>
#include <quoc.h>
#include <op.h>
#include <discreteFunction.h>
#include <iterators.h>
#include <baseFunctionSet.h>
#include <simplexBaseFunctionSet.h>
#include <simplexConfigurators.h>

namespace qc {


// --------------------------------- PART 1: INTERFACES -----------------------------------------------------------------


//! Interface class for computing \f$ \int_{\partial \Omega} f \psi \, da \f$.
/** The function \f$ f \f$
 *  has to be implemented in the derived class by the method integrandAtQuadPoint( arg, normal, el, edgePoint ).
 *  The domain \f$ \Omega \f$ is a 2d or 3d-quocmesh-grid.
 *  \ingroup FENonlinOpInt
 */
template <typename ConfiguratorType, typename VecType, qc::Dimension, typename QuadType, typename Imp>
class BoundaryIntegrationInterface : public aol::FEOpInterface<ConfiguratorType, VecType, VecType> { };

//! Interface class for computing \f$ \int_{\partial \Omega} f \psi \, da \f$.
/** The function \f$ f \f$
 *  has to be implemented in the derived class by the method integrandAtQuadPoint( arg, normal, el, edgePoint ).
 *  The domain \f$ \Omega \f$ is a 2d or 3d-quocmesh-grid.
 *  2d-specification
 *  QuadType has to be a 1D quadrature type
 *  \ingroup FENonlinOpInt
 */
template <typename ConfiguratorType, typename VecType, typename QuadType, typename Imp>
class BoundaryIntegrationInterface<ConfiguratorType, VecType, qc::QC_2D, QuadType, Imp> :
      public aol::FEOpInterface<ConfiguratorType, VecType, VecType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::InitType GridType;

protected:
  const GridType &_grid;

public:
  explicit BoundaryIntegrationInterface ( const GridType &grid )
      : aol::FEOpInterface<ConfiguratorType, VecType, VecType> ( grid ), _grid ( grid ) {}

  explicit BoundaryIntegrationInterface ( const GridType &grid, const ConfiguratorType &configurator )
      : aol::FEOpInterface<ConfiguratorType, VecType, VecType> ( configurator ), _grid ( grid ) {}

  //! this method computes \f$ f \f$ and has to be implemented in the derived class
  RealType integrandAtQuadPoint ( const VecType &arg, const aol::Vec2<RealType> &normal, const ElementType &/*el*/, const aol::Vec2<RealType> edgePoint ) const {
    return asImp().integrandAtQuadPoint ( arg, normal, edgePoint );
  }

  void applyAdd ( const VecType &arg, VecType &dest ) const {
    // This lookUp-table represents the four corners of the the unit-square and is needed
    // for the computation of the boundary-vertex.
    static RealType lookUp[4][2] = { {0, 0}, {1, 0}, {0, 1}, {1, 1} };

    // OldBoundaryHyperfaceIterator2D is an element-iterator that traverses the boundary of the unit-square
    for ( typename GridType::OldBoundaryHyperfaceIterator2D eit = _grid.ebegin(); eit; ++eit ) {
      // _locNode stores, which nodes of the element belong to the boundary
      aol::Vec2<RealType> edgePt1 ( lookUp[eit->_locNode[0]][0], lookUp[eit->_locNode[0]][1] );
      aol::Vec2<RealType> edgePt2 ( lookUp[eit->_locNode[1]][0], lookUp[eit->_locNode[1]][1] );

      // computation of the normal
      aol::Vec2<RealType> edgePt ( edgePt1 );
      edgePt += edgePt2;
      edgePt /= 2.;

      aol::Vec2<RealType> normal ( edgePt[0] - 0.5, edgePt[1] - 0.5 );
      normal *= 2.;

      QuadType quad;    // this is 1d-quadrature-rule as the boundary are 1d-structures
      for ( int q = 0; q < QuadType::numQuadPoints; q++ ) {
        /*** interpolate ***/
        RealType rc = quad.getRefCoord ( q ) [0];

        // compute the interpolation point
        for ( int i = 0; i < 2; i++ ) {
          edgePt[i] = edgePt1[i] * ( 1. - rc ) + rc * edgePt2[i];
        }
        // get the integrand at this point
        RealType integrand = asImp().integrandAtQuadPoint ( arg, normal, *eit, edgePt );

        // traverse the local dof's compute the complete integrand incl. basis-function
//         for ( int b = 0; b < _config.getNumLocalDofs ( *eit ); b++ ) {
        for ( int nodeIndex = 0; nodeIndex < 2; nodeIndex++ ) {
          int b = eit->_locNode[nodeIndex];
          dest[ this->getConfigurator().localToGlobal ( *eit, b ) ] += quad.getWeight ( q )
                                                       * integrand * this->getConfigurator().H ( *eit )
                                                       * this->getConfigurator().getBaseFunctionSet ( *eit ).evaluate ( b, edgePt );
        }
      }
    }   // end element-loop
  }

protected:
  /** barton-nackman **/
  Imp& asImp( ) { return static_cast<Imp&> ( *this ); }
  const Imp& asImp( ) const { return static_cast<const Imp&> ( *this ); }

};


// class for computing boundary integrals
// requires the grid to supply an RectangularBoundaryFaceElementIterator (iterator over boundary face elements)
// ConfiguratorType specifies your elements (e.g. trilinear, iterator through elements), quadrature, index mapping.
//! 3d-specification
//! \ingroup FENonlinOpInt
template <typename ConfiguratorType, typename VecType, typename QuadType, typename Imp>
class BoundaryIntegrationInterface<ConfiguratorType, VecType, qc::QC_3D, QuadType, Imp> :
      public aol::FEOpInterface<ConfiguratorType, VecType, VecType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::InitType GridType;

protected:
  const GridType &_grid;

public:

  BoundaryIntegrationInterface ( const ConfiguratorType &config, const GridType &grid )
      : aol::FEOpInterface<ConfiguratorType, VecType, VecType> ( config ), _grid ( grid ) {}

  RealType integrandAtQuadPoint ( const VecType &arg, const aol::Vec3<RealType> &normal, const ElementType &/*el*/, const aol::Vec3<RealType> edgePoint ) const {
    return asImp().integrandAtQuadPoint ( arg, normal, edgePoint );
  }

  void applyAdd ( const VecType &arg, VecType &dest ) const {

    for ( RectangularBoundaryFaceElementIterator<ConfiguratorType> bfeIterator ( _grid ); bfeIterator.notAtEnd (); ++bfeIterator ) {
      aol::Vec3<RealType> normal ( bfeIterator->getNormal() );

      QuadType quad;
      for ( int q = 0; q < QuadType::numQuadPoints; q++ ) {
        /*** interpolate ***/
        const aol::Vec2<RealType> &rc = quad.getRefCoord ( q );
        aol::Vec3<RealType> refCoord ( bfeIterator->getRefCoordOnElementFromRefCoordOnFace ( rc ) );

        // get the integrand
        RealType integrand = asImp().integrandAtQuadPoint ( arg, normal, *bfeIterator, refCoord );

        // traverse the local dof's and compute the complete integrand incl. basis function
//         for ( int b = 0; b < _config.getNumLocalDofs ( eit->_el ); b++ ) {
        for ( int nodeIndex = 0; nodeIndex < 4; nodeIndex++ ) {
          int b = bfeIterator->getLocalNodeOnBoundary ( nodeIndex );
          dest[ this->getConfigurator().localToGlobal ( *bfeIterator, b ) ] += quad.getWeight ( q ) * integrand * aol::Sqr ( this->getConfigurator().H ( *bfeIterator ) )
                                                               * this->getConfigurator().getBaseFunctionSet ( *bfeIterator ).evaluate ( b, refCoord );
        }
      }
    }
  }

protected:
  /** barton-nackman **/
  Imp& asImp( ) { return static_cast<Imp&> ( *this ); }
  const Imp& asImp( ) const { return static_cast<const Imp&> ( *this ); }

};




//! \brief Interface class for computing \f$ \int_{\partial \Omega} f \nabla \psi \cdot \nu \, da \f$. The function \f$ f \f$
//!        has to be implemented in the derived class by the method integrandAtQuadPoint( arg, normal, el, edgePoint ).
//!
//! The domain \f$ \Omega \f$ is a 2d or 3d-quocmesh-grid.
//! \ingroup FENonlinOpInt
template <typename ConfiguratorType, typename VecType, qc::Dimension, typename QuadType, typename Imp>
class BoundaryDiffIntegrationInterface : public aol::Op<VecType, VecType> { };

//! \brief Interface class for computing \f$ \int_{\partial \Omega} f \nabla \psi \cdot \nu \, da \f$. The function \f$ f \f$
//!        has to be implemented in the derived class by the method integrandAtQuadPoint( arg, normal, el, edgePoint ).
//!
//! The domain \f$ \Omega \f$ is a 2d or 3d-quocmesh-grid.
//! 2d-specification
//! \ingroup FENonlinOpInt
template <typename ConfiguratorType, typename VecType, typename QuadType, typename Imp>
class BoundaryDiffIntegrationInterface<ConfiguratorType, VecType, qc::QC_2D, QuadType, Imp> :
      public aol::Op<VecType, VecType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::InitType GridType;

protected:
  const ConfiguratorType _config;
  const GridType &_grid;

public:
  BoundaryDiffIntegrationInterface ( const GridType &grid )
      : _config ( grid ), _grid ( grid ) {}

  //! this method computes \f$ f \f$ and has to be implemented in the derived class
  RealType integrandAtQuadPoint ( const VecType &arg, const aol::Vec2<RealType> &normal, const ElementType &/*el*/, const aol::Vec2<RealType> edgePoint ) const {
    return asImp().integrandAtQuadPoint ( arg, normal, edgePoint );
  }

  void applyAdd ( const VecType &arg, VecType &dest ) const {
    cerr << aol::color::red << "This operator uses an additional multiplication with h which seems to be correct but has to be justified!\n";
    cerr << aol::color::reset;
    // This lookUp-table represents the four corners of the the unit-square and is needed
    // for the computation of the boundary-vertex.
    static RealType lookUp[4][2] = { {0, 0}, {1, 0}, {0, 1}, {1, 1} };

    // OldBoundaryHyperfaceIterator2D is an element-iterator that traverses the boundary of the unit-square
    for ( typename GridType::OldBoundaryHyperfaceIterator2D eit = _grid.ebegin(); eit; ++eit ) {
      // _locNode stores, which nodes of the element belong to the boundary
      aol::Vec2<RealType> edgePt1 ( lookUp[eit->_locNode[0]][0], lookUp[eit->_locNode[0]][1] );
      aol::Vec2<RealType> edgePt2 ( lookUp[eit->_locNode[1]][0], lookUp[eit->_locNode[1]][1] );

      // computation of the normal
      aol::Vec2<RealType> edgePt ( edgePt1 );
      edgePt += edgePt2;
      edgePt /= 2.;

      aol::Vec2<RealType> normal ( edgePt[0] - 0.5, edgePt[1] - 0.5 );
      normal *= 2.;

      QuadType quad;                // this is 1d-quadrature-rule as the boundary are 1d-structures
      aol::Vec2<RealType> grad;     // gradient of the basis-function

      for ( int q = 0; q < QuadType::numQuadPoints; q++ ) {
        /*** interpolate ***/
        RealType rc = quad.getRefCoord ( q );

        // compute the interpolation point
        for ( int i = 0; i < 2; i++ ) {
          edgePt[i] = edgePt1[i] * ( 1. - rc ) + rc * edgePt2[i];
        }
        // evaluate the integrand, that means the function f, at this point
        RealType integrand = asImp().integrandAtQuadPoint ( arg, normal, *eit, edgePt );

        // traverse the local dof's and compute the complete integrand incl. basis-function
//         for ( int b = 0; b < _config.getNumLocalDofs ( *eit ); b++ ) {
        for ( int nodeIndex = 0; nodeIndex < 2; nodeIndex++ ) {
          int b = eit->_locNode[nodeIndex];
          // get the integrand at this point
          _config.getBaseFunctionSet( *eit ).evaluateGradient ( b, edgePt, grad );
          grad *= _config.H( *eit );      // HACK: But this seems to be correct, exact reason therefore?
//           cerr << grad << ", ";

          dest[ _config.localToGlobal ( *eit, b ) ] += ( normal * grad ) * quad.getWeight ( q )* _config.H ( *eit ) * integrand;
        }   // end localDof-loop
      }   // end quad-point-loop
    }   // end element iterator
  }   // end applyAdd

protected:
  /** barton-nackman **/
  Imp& asImp( ) { return static_cast<Imp&> ( *this ); }
  const Imp& asImp( ) const { return static_cast<const Imp&> ( *this ); }

};

//! \brief Interface class for computing \f$ \int_{\partial \Omega} f \nabla \psi \cdot \nu \, da \f$. The function \f$ f \f$
//!        has to be implemented in the derived class by the method integrandAtQuadPoint( arg, normal, el, edgePoint ).
//!
//! The domain \f$ \Omega \f$ is a 2d or 3d-quocmesh-grid.
//! 3d-specification
//! \ingroup FENonlinOpInt
template <typename ConfiguratorType, typename VecType, typename QuadType, typename Imp>
class BoundaryDiffIntegrationInterface<ConfiguratorType, VecType, qc::QC_3D, QuadType, Imp> :
      public aol::Op<VecType, VecType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::InitType GridType;

protected:
  const ConfiguratorType _config;
  const GridType &_grid;

public:
  BoundaryDiffIntegrationInterface ( const GridType &grid )
      : _config ( grid ), _grid ( grid ) {}

  //! this method computes \f$ f \f$ and has to be implemented in the derived class
  RealType integrandAtQuadPoint ( const VecType &arg, const aol::Vec3<RealType> &normal, const ElementType &/*el*/, const aol::Vec3<RealType> edgePoint ) const {
    return asImp().integrandAtQuadPoint ( arg, normal, edgePoint );
  }

  void applyAdd ( const VecType &arg, VecType &dest ) const {

    for ( RectangularBoundaryFaceElementIterator<ConfiguratorType> bfeIterator ( _grid ); bfeIterator.notAtEnd (); ++bfeIterator ) {
      aol::Vec3<RealType> normal ( bfeIterator->getNormal() );

      QuadType quad;                // this is 2d-quadrature-rule as the boundary are 2d-structures
      aol::Vec3<RealType> grad;     // gradient of the basis-function

      for ( int q = 0; q < QuadType::numQuadPoints; q++ ) {
        /*** interpolate ***/
        const aol::Vec2<RealType> &rc = quad.getRefCoord ( q );
        aol::Vec3<RealType> refCoord(bfeIterator->getRefCoordOnElementFromRefCoordOnFace(rc)) /*( bfeIterator->getRefCoordIn3DFromRefCoordIn2D ( rc ) )*/;

        // evaluate the integrand, that means the function f, at this point
        RealType integrand = asImp().integrandAtQuadPoint ( arg, normal, *bfeIterator, refCoord );

        // traverse the local dof's and compute the complete integrand incl. basis-function
        //         for ( int b = 0; b < _config.getNumLocalDofs ( *eit ); b++ ) {
        for ( int nodeIndex = 0; nodeIndex < 8; nodeIndex++ ) {
          // get the integrand at this point
          _config.getBaseFunctionSet( *bfeIterator ).evaluateGradient ( nodeIndex, refCoord, grad );

          dest[ _config.localToGlobal ( *bfeIterator, nodeIndex ) ] += ( normal * grad ) * quad.getWeight ( q )* aol::Sqr ( _config.H ( *bfeIterator ) ) * integrand;
        }   // end localDof-loop
      }   // end quad-point-loop
    }  // end element iterator
  }   // end applyAdd

protected:
  /** barton-nackman **/
  Imp& asImp( ) { return static_cast<Imp&> ( *this ); }
  const Imp& asImp( ) const { return static_cast<const Imp&> ( *this ); }

};


//! Interface class for computing \f$ \int_{\partial \Omega} f(...) \, A(...) \nabla \psi \cdot \nu \, da \f$.
/** The function \f$ f \f$
 *  has to be implemented in the derived class by the method getCoeffAtQuadPoint( arg, normal, el, edgePoint ), as well as the
 *  matrix A(...) by the function getCoeffMatrix( arg, normal, e1, edgePoint ).
 *  The domain \f$ \Omega \f$ is a 2d or 3d-quocmesh-grid.
 *  \ingroup FENonlinOpInt
 */
template <typename ConfiguratorType, typename VecType, qc::Dimension, typename QuadType, typename Imp>
class BoundaryWeightedMatrixDiffIntegrationInterface : public aol::Op<VecType, VecType> { };

//! Interface class for computing \f$ \int_{\partial \Omega} f(...) \, A(...) \nabla \psi \cdot \nu \, da \f$.
/** The function \f$ f \f$
 *  has to be implemented in the derived class by the method getCoeff( arg, normal, el, edgePoint ), as well as the
 *  matrix A(...) by the function getCoeffMatrix( arg, normal, e1, edgePoint ).
 *  The domain \f$ \Omega \f$ is a 2d or 3d-quocmesh-grid.
 *  2d-specification
 *  \ingroup FENonlinOpInt
 */
template <typename ConfiguratorType, typename VecType, typename QuadType, typename Imp>
class BoundaryWeightedMatrixDiffIntegrationInterface<ConfiguratorType, VecType, qc::QC_2D, QuadType, Imp> :
      public aol::Op<VecType, VecType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::InitType GridType;
//   typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

protected:
  const ConfiguratorType _config;
  const GridType &_grid;

public:
  BoundaryWeightedMatrixDiffIntegrationInterface ( const GridType &grid )
      : _config ( grid ), _grid ( grid ) {}

  const ConfiguratorType & getConfigurator () const {
    return _config;
  }

  //! this method computes \f$ f \f$ and has to be implemented in the derived class
  RealType getCoeff ( const VecType &arg, const aol::Vec2<RealType> &normal, const ElementType &/*el*/, const aol::Vec2<RealType> edgePoint ) const {
    return asImp().integrandAtQuadPoint ( arg, normal, edgePoint );
  }

  //! this function computes the matrix \f$ A \f$  in the derived class
  inline void getCoeffMatrix ( const VecType &arg, const aol::Vec2<RealType> &normal, const ElementType &el,
                               const aol::Vec2<RealType> edgePoint, MatType &Matrix ) const {
    this->asImp().getCoeffMatrix ( arg, normal, el, edgePoint, Matrix );
  }

  void applyAdd ( const VecType &arg, VecType &dest ) const {
    cerr << aol::color::red << "This operator uses an additional multiplication with h which seems to be correct but has to be justified!\n";
    cerr << aol::color::reset;
    // This lookUp-table represents the four corners of the the unit-square and is needed
    // for the computation of the boundary-vertex.
    static RealType lookUp[4][2] = { {0, 0}, {1, 0}, {0, 1}, {1, 1} };

    // OldBoundaryHyperfaceIterator2D is an element-iterator that traverses the boundary of the unit-square
    for ( typename GridType::OldBoundaryHyperfaceIterator2D eit = _grid.ebegin(); eit; ++eit ) {
      // _locNode stores, which nodes of the element belong to the boundary
      aol::Vec2<RealType> edgePt1 ( lookUp[eit->_locNode[0]][0], lookUp[eit->_locNode[0]][1] );
      aol::Vec2<RealType> edgePt2 ( lookUp[eit->_locNode[1]][0], lookUp[eit->_locNode[1]][1] );

      // computation of the normal
      aol::Vec2<RealType> edgePt ( edgePt1 );
      edgePt += edgePt2;
      edgePt /= 2.;

      aol::Vec2<RealType> normal ( edgePt[0] - 0.5, edgePt[1] - 0.5 );
      normal *= 2.;

      QuadType quad;                // this is 1d-quadrature-rule as the boundary are 1d-structures
      aol::Vec2<RealType> grad;     // gradient of the basis-function
      aol::Vec2<RealType> matGrad;  // A*gradient
      MatType mat;

      // traverse the quadpoints
      for ( int q = 0; q < QuadType::numQuadPoints; q++ ) {
        /*** interpolate ***/
        RealType rc = quad.getRefCoord ( q );

        // compute the interpolation point
        for ( int i = 0; i < 2; i++ ) {
          edgePt[i] = edgePt1[i] * ( 1. - rc ) + rc * edgePt2[i];
        }

        // traverse the local dof's and compute the complete integrand incl. basis-function
//         for ( int b = 0; b < _config.getNumLocalDofs ( *eit ); b++ ) {
        for ( int nodeIndex = 0; nodeIndex < 2; nodeIndex++ ) {
          int b = eit->_locNode[nodeIndex];
          // get the integrand at this point
          RealType coeff = asImp().getCoeff( arg, normal, *eit, edgePt );
          asImp().getCoeffMatrix( arg, normal, *eit, edgePt, mat );
          _config.getBaseFunctionSet( *eit ).evaluateGradient ( b, edgePt, grad );

          grad *= _config.H( *eit );              // HACK

          mat.mult( grad, matGrad );        // A*grad

          dest[ _config.localToGlobal ( *eit, b ) ] += ( normal * matGrad ) * quad.getWeight ( q )* _config.H ( *eit ) * coeff;
        }   // end localDof-loop
      }   // end quad-point-loop
    }   // end element iterator
  }   // end applyAdd

protected:
  /** barton-nackman **/
  Imp& asImp( ) { return static_cast<Imp&> ( *this ); }
  const Imp& asImp( ) const { return static_cast<const Imp&> ( *this ); }

};




// --------------------------------- PART 2: APPLICABLE OPERATORS -----------------------------------------------------------------

/**
 * This class computes the following boundary integral:
 * \f[ \int_{\partial\Omega} \nabla u \cdot n \, \varphi da \f]
 *
 */
template <typename ConfiguratorType, typename BoundaryQuadratureType>
class IntegrateGradUOverBoundary
  : public BoundaryIntegrationInterface < ConfiguratorType,
                                          aol::Vector<typename ConfiguratorType::RealType>,
                                          ConfiguratorType::Dim,
                                          BoundaryQuadratureType, IntegrateGradUOverBoundary< ConfiguratorType, BoundaryQuadratureType > > {

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;

public:
  IntegrateGradUOverBoundary( const typename ConfiguratorType::InitType &Initializer )
    : BoundaryIntegrationInterface< ConfiguratorType,
                                    aol::Vector<RealType>,
                                    ConfiguratorType::Dim,
                                    BoundaryQuadratureType,
                                    IntegrateGradUOverBoundary< ConfiguratorType, BoundaryQuadratureType > >( Initializer ) {
  }

  // the only interesting function here: evalute gradU * normal
  RealType integrandAtQuadPoint( const aol::Vector<RealType> &arg, const typename ConfiguratorType::VecType &normal, const ElementType &el, const typename ConfiguratorType::VecType edgePoint ) const {
    aol::DiscreteFunctionDefault<ConfiguratorType> discrU( this->_grid, arg );

    // compute gradient of boundaryfunction @ edgePt
    typename ConfiguratorType::VecType gradU;
    discrU.evaluateGradient( el, edgePoint, gradU );

    return gradU * normal;
  }
};


/**
 * This class computes the following boundary integral:
 * \f[ \int_{\partial\Omega} u \, \nabla \vartheta \cdot n da \f]
 *
 */
template <typename ConfiguratorType, typename BoundaryQuadratureType>
class IntegrateUDiffBasisFunctionOverBoundary
  : public BoundaryDiffIntegrationInterface < ConfiguratorType,
                                              aol::Vector<typename ConfiguratorType::RealType>,
                                              ConfiguratorType::Dim,
                                              BoundaryQuadratureType, IntegrateUDiffBasisFunctionOverBoundary< ConfiguratorType, BoundaryQuadratureType > > {

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;

public:
  IntegrateUDiffBasisFunctionOverBoundary( const typename ConfiguratorType::InitType &Initializer )
    : BoundaryDiffIntegrationInterface< ConfiguratorType,
                                        aol::Vector<RealType>,
                                        ConfiguratorType::Dim,
                                        BoundaryQuadratureType,
                                        IntegrateUDiffBasisFunctionOverBoundary< ConfiguratorType, BoundaryQuadratureType > >( Initializer ) {
  }

  // the only interesting function here: evaluate u at the current coordinates
  RealType integrandAtQuadPoint( const aol::Vector<RealType> &arg, const typename ConfiguratorType::VecType &/*normal*/, const ElementType &el, const typename ConfiguratorType::VecType edgePoint ) const {
    aol::DiscreteFunctionDefault<ConfiguratorType> discrU( this->_grid, arg );

    return discrU.evaluate( el, edgePoint );
  }
};

template < class QuadType, typename RealType /*= QuadType::RealType*/, qc::Dimension Dim /*= QuadType::Dim*/, int Order /*= QuadType::Order*/ >
class Codim1QuadratureTrait {};

template < typename RealType, qc::Dimension Dim, int Order >
class Codim1QuadratureTrait<aol::GaussQuadrature<RealType,Dim,Order>,RealType,Dim,Order> {
public:
  typedef aol::GaussQuadrature<RealType,static_cast<qc::Dimension>( Dim - qc::QC_1D ),Order> QuadType;
};

template < typename RealType, qc::Dimension Dim, int Order >
class Codim1QuadratureTrait<qc::simplex::MidpointQuadrature<RealType, Dim>,RealType,Dim,Order> {
public:
  typedef qc::simplex::MidpointQuadrature<RealType,static_cast<qc::Dimension>( Dim - qc::QC_1D )> QuadType;
};

//! Interface for linear FE operators which only operate on the domain boundary (implemented for a rectangular grid).
/**
 *  \author Wirth
 */
template <typename ConfiguratorType, typename Imp, typename IteratorType = RectangularBoundaryFaceElementIterator<ConfiguratorType> >
class FEBoundaryOpInterface :
  public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  aol::DeleteFlagPointer<const ConfiguratorType> _config;
  mutable typename ConfiguratorType::MatrixType *_mat;
  aol::OperatorType _opType;

public:
  explicit FEBoundaryOpInterface ( const typename ConfiguratorType::InitType &Grid, const aol::OperatorType OpType = aol::ONTHEFLY ) :
    _config ( new ConfiguratorType( Grid ), true ), _mat ( NULL ), _opType ( OpType ) {}

  explicit FEBoundaryOpInterface ( const ConfiguratorType * Conf, const aol::OperatorType OpType = aol::ONTHEFLY ) :
    _config ( Conf, false ), _mat ( NULL ), _opType ( OpType ) {}

  virtual ~FEBoundaryOpInterface( ) {
    delete _mat;
  }

  void reset( ) {
    if ( _mat )
      delete _mat;
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    switch ( _opType ) {
    case aol::ONTHEFLY:
      multiplyOnTheFly ( Arg, Dest );
      break;
    case aol::ASSEMBLED:
      if ( !_mat )
        assembleMatrix();
      _mat->applyAdd ( Arg, Dest );
      break;
    default:
      throw aol::UnimplementedCodeException ( "FEBoundaryOpInterface::applyAdd: unknown opType", __FILE__, __LINE__ );
    };
  }

  typename ConfiguratorType::MatrixType& getMatrix( ) {
    if ( !_mat ) {
      _mat = _config->createNewMatrix( );
      assembleMatrix( );
    }
    //return dynamic_cast<typename ConfiguratorType::MatrixType&>(*_mat);
    return *_mat;
  }

protected:
  void multiplyOnTheFly ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {

    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;
    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];

    // traverse the boundary faces of the boundary elements of the grid and on each face the quadrature points
    for ( IteratorType bfeIterator ( _config->getInitializer() ); bfeIterator.notAtEnd (); ++bfeIterator ) {
      // assemble the local matrix belonging to the current element
      this->asImp().prepareLocalMatrix ( *bfeIterator, localMatrix );

      // get the global indices to the current Dofs
      const int numLocalDofs = _config->getNumLocalDofs ( *bfeIterator );
      for ( int i = 0; i < numLocalDofs; ++i )
        globalDofs[ i ] = _config->localToGlobal ( *bfeIterator, i );

      // add the locally computed value to the global result
      for ( int i = 0; i < numLocalDofs; ++i )
        for ( int j = 0; j < numLocalDofs; ++j )
          Dest[ globalDofs[ i ] ] += localMatrix [ i ][ j ] * Arg[ globalDofs[ j ] ] ;
    }
  }

  void assembleMatrix( ) const {
    if ( _mat )
      delete _mat;
    _mat = _config->createNewMatrix( );
    assembleAddMatrix ( *_mat );
  }

public:
  template <typename MatrixType>
void assembleAddMatrix ( MatrixType &Mat, const RealType Factor = aol::NumberTrait<RealType>::one ) const {          //Diese wird aufgerufen

    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;
    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];

    // traverse the boundary faces of the boundary elements of the grid and on each face the quadrature points
    for ( IteratorType bfeIterator ( _config->getInitializer() ); bfeIterator.notAtEnd (); ++bfeIterator ) {
      // assemble the local matrix belonging to the current element
      this->asImp().prepareLocalMatrix ( *bfeIterator, localMatrix );

      // get the global indices to the current Dofs
      const int numLocalDofs = _config->getNumLocalDofs ( *bfeIterator );

      for ( int i = 0; i < numLocalDofs; ++i )
        globalDofs[ i ] = _config->localToGlobal ( *bfeIterator, i );

      // add the locally assembled matrix to the globally assembled matrix
      for ( int i = 0; i < numLocalDofs; ++i )
        for ( int j = 0; j < numLocalDofs; ++j )
          Mat.add ( globalDofs[ i ], globalDofs[ j ], Factor * localMatrix [ i ][ j ] );
    }
  }

  template <typename MatrixType >
  void assembleAddMatrix ( const aol::BitVector &Mask, MatrixType &Mat, const RealType Factor = aol::NumberTrait<RealType>::one, const bool SetMaskedRowsToIdentity = true ) const {

    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;
    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];

    // traverse the boundary faces of the boundary elements of the grid and on each face the quadrature points
    for ( IteratorType bfeIterator ( _config->getInitializer() ); bfeIterator.notAtEnd (); ++bfeIterator ) {
      // assemble the local matrix belonging to the current element
      this->asImp().prepareLocalMatrix ( *bfeIterator, localMatrix );

      // get the global indices to the current Dofs
      const int numLocalDofs = _config->getNumLocalDofs ( *bfeIterator );
      for ( int i = 0; i < numLocalDofs; ++i )
        globalDofs[ i ] = _config->localToGlobal ( *bfeIterator, i );

      // add the locally assembled matrix to the globally assembled matrix,
      // but write the ith row only if the ith node is not a Dirichlet node and the jth colume only if the jth node is not a Dirichlet node
      for ( int i = 0; i < numLocalDofs; ++i )
        if ( !Mask[ globalDofs[ i ] ] )
          for ( int j = 0; j < numLocalDofs; ++j )
            if ( !Mask[ globalDofs[ j ] ] )
              Mat.add ( globalDofs[ i ], globalDofs[ j ], Factor * localMatrix [ i ][ j ] );
    }

    // set ones on the diagonal for Dirichlet nodes
    if ( SetMaskedRowsToIdentity )
      for ( int i = 0; i < Mask.size(); i++ )
        if ( Mask[i] )
          Mat.add( i, i, Factor );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

//! Interface to compute \f$(\int_{\partial\Omega} w(x)\varphi_i(x)\varphi_j(x) dx)_{ij}\f$,
/** where \f$\partial\Omega\f$ is the domain boundary, \f$\varphi_i\f$ is the \f$i\f$th FE basis function,
 *  and the scalar weight \f$w\f$ has to be provided in derived classes.
 *
 *  \author Wirth
 *  \ingroup FELinOpInt
 */
template <typename ConfiguratorType, typename Imp, typename IteratorType = RectangularBoundaryFaceElementIterator<ConfiguratorType> >
class FELinBoundaryScalarWeightedMassInterface :
  public FEBoundaryOpInterface<ConfiguratorType,FELinBoundaryScalarWeightedMassInterface<ConfiguratorType,Imp,IteratorType>,IteratorType > {

protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

public:
  explicit FELinBoundaryScalarWeightedMassInterface ( const typename ConfiguratorType::InitType &Grid, const aol::OperatorType OpType = aol::ONTHEFLY )
  : FEBoundaryOpInterface<ConfiguratorType,FELinBoundaryScalarWeightedMassInterface<ConfiguratorType,Imp,IteratorType>,IteratorType > ( Grid, OpType )
  {}

  explicit FELinBoundaryScalarWeightedMassInterface ( const ConfiguratorType * Conf, const aol::OperatorType OpType = aol::ONTHEFLY )
  : FEBoundaryOpInterface<ConfiguratorType,FELinBoundaryScalarWeightedMassInterface<ConfiguratorType,Imp,IteratorType>,IteratorType > ( Conf, OpType )
  {}

  //! This function has to be provided in the implementation (derived class) of the interface.
  inline RealType getCoeff ( const qc::BoundaryFaceElement<RealType,ConfiguratorType::Dim> &El, const DomVecType &Local3DCoord ) const {
    return this->asImp().getCoeff ( El, Local3DCoord );
  }

  //! Performs the numerical quadrature of the bilinear form and saves the values locally.
void prepareLocalMatrix ( const typename IteratorType::IteratedType &El, aol::Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> &LocalMatrix ) const {
    LocalMatrix.setZero();
    aol::Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> localMatrix;

    const int numDofs = this->_config->getNumLocalDofs ( El );
    const typename ConfiguratorType::BaseFuncSetType &bfs = this->_config->getBaseFunctionSet ( El );
    aol::Vec<ConfiguratorType::maxNumLocalDofs,RealType> basisFunctionValues;

    // perform the quadrature
    typedef typename qc::Codim1QuadratureTrait<typename ConfiguratorType::QuadType,typename ConfiguratorType::QuadType::RealType,ConfiguratorType::QuadType::Dim,ConfiguratorType::QuadType::Order>       ::QuadType        QuadType;
    QuadType quadrature;
    for ( int q = 0; q < QuadType::numQuadPoints; q++ ) {
      // compute the weight and the basis function values
      const DomVecType local3DCoord ( El.getRefCoordOnElementFromRefCoordOnFace ( quadrature.getRefCoord ( q ) ) );
      const RealType coeff = this->asImp().getCoeff ( El, local3DCoord );
      for ( int i = 0; i < numDofs; i++ )
        basisFunctionValues[i] = bfs.evaluate( this->_config->localOnFaceToLocal(El,i), local3DCoord );

      // assemble the local mass matrix
      localMatrix.makeTensorProduct( basisFunctionValues, basisFunctionValues );
      localMatrix *= coeff * quadrature.getWeight ( q ) * ( this->_config->vol ( El ) / this->_config->H( El ) );
      LocalMatrix += localMatrix;
    }
  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }
};

//! Computes \f$(\int_{\partial\Omega} \varphi_i(x)\varphi_j(x) dx)_{ij}\f$,
/** where \f$\partial\Omega\f$ is the domain boundary, \f$\varphi_i\f$ is the \f$i\f$th FE basis function.
 *
 *  \author Wirth
 *  \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, typename IteratorType = RectangularBoundaryFaceElementIterator<ConfiguratorType> >
class FELinBoundaryMassOp :
  public FELinBoundaryScalarWeightedMassInterface<ConfiguratorType,FELinBoundaryMassOp<ConfiguratorType,IteratorType>,IteratorType > {

  typedef typename ConfiguratorType::RealType RealType;

public:
  explicit FELinBoundaryMassOp ( const typename ConfiguratorType::InitType &Grid, const aol::OperatorType OpType = aol::ONTHEFLY )
  : FELinBoundaryScalarWeightedMassInterface<ConfiguratorType,FELinBoundaryMassOp<ConfiguratorType,IteratorType>,IteratorType > ( Grid, OpType )
  {}

  explicit FELinBoundaryMassOp ( const ConfiguratorType * Conf, const aol::OperatorType OpType = aol::ONTHEFLY )
  : FELinBoundaryScalarWeightedMassInterface<ConfiguratorType,FELinBoundaryMassOp<ConfiguratorType,IteratorType>,IteratorType > ( Conf, OpType )
  {}

  inline RealType getCoeff ( const typename IteratorType::IteratedType &/*El*/, const typename ConfiguratorType::DomVecType &/*Local3DCoord*/ ) const {
    return 1.;
  }
};

//! Computes \f$(\int_{\partial\Omega} \chi(x) \varphi_i(x)\varphi_j(x) dx)_{ij}\f$,
/** where \f$\partial\Omega\f$ is the domain boundary, \f$\chi(x)\f$ is a characteristic function selecting parts of the boundary, \f$\varphi_i\f$ is the \f$i\f$th FE basis function.
 *
 *  \author Geihe
 *  \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, typename IteratorType = RectangularBoundaryFaceElementIterator<ConfiguratorType> >
class FELinPartialBoundaryMassOp :
public FELinBoundaryScalarWeightedMassInterface<ConfiguratorType,FELinPartialBoundaryMassOp<ConfiguratorType,IteratorType>,IteratorType > {

  typedef typename ConfiguratorType::RealType RealType;

  public:
    FELinPartialBoundaryMassOp ( const typename ConfiguratorType::InitType &Grid, const aol::BitVector & FaceMask, const aol::OperatorType OpType = aol::ONTHEFLY ) :
    FELinBoundaryScalarWeightedMassInterface<ConfiguratorType,FELinPartialBoundaryMassOp<ConfiguratorType,IteratorType>,IteratorType > ( Grid, OpType ),
    _faceMask( FaceMask )
    {}

    FELinPartialBoundaryMassOp ( const ConfiguratorType * Conf, const aol::BitVector & FaceMask, const aol::OperatorType OpType = aol::ONTHEFLY ) :
    FELinBoundaryScalarWeightedMassInterface<ConfiguratorType,FELinPartialBoundaryMassOp<ConfiguratorType,IteratorType>,IteratorType > ( Conf, OpType ),
    _faceMask( FaceMask )
    {}

    inline RealType getCoeff ( const typename IteratorType::IteratedType & El, const typename ConfiguratorType::DomVecType &/*Local3DCoord*/ ) const {
      if ( _faceMask[ El.getBoundaryFaceType() ] )
        return 0.;
      else
        return 1.;
      return 0.; // to avoid warning
    }

  protected:
    const aol::BitVector & _faceMask;
};

} // end namespace qc

#endif
