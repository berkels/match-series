#ifndef __IMPLICITSURFACEFEOPS_H
#define __IMPLICITSURFACEFEOPS_H

#include <FEOpInterface.h>
#include <baseFunctionSet.h>
#include <qmException.h>
#include <discreteFunction.h>
#include <boundaryIntegration.h>

namespace qc {

// ------------------------------ OPERATORS ON THE WHOLE GRID ---------------------------------------------------

/** \brief This class computes the following weighted mass op:
 *         \f[ \int_{\Omega} \| \nabla \Phi \|_\epsilon \vartheta_i \vartheta_j \, dx \f]
 *
 * \todo This is essentially aol::SquaredDiffWeightMassOp without the square. Rename to aol::DiffWeightMassOp?
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class WeightedMassOp : public aol::FELinScalarWeightedMassInterface<ConfiguratorType, WeightedMassOp<ConfiguratorType, IndexMode>, IndexMode > {
public:
  typedef typename ConfiguratorType::RealType RealType;

protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrImg;
  RealType _epsSqr;

public:
  WeightedMassOp ( const typename ConfiguratorType::InitType &Initializer, aol::OperatorType OpType = aol::ONTHEFLY, RealType Epsilon = 1. )
      : aol::FELinScalarWeightedMassInterface<ConfiguratorType, WeightedMassOp<ConfiguratorType, IndexMode>, IndexMode > ( Initializer, OpType ),
      _discrImg ( NULL ), _epsSqr ( Epsilon*Epsilon )  {}

  WeightedMassOp ( const ConfiguratorType & Config, const typename ConfiguratorType::InitType &Initializer, aol::OperatorType OpType = aol::ONTHEFLY, RealType Epsilon = 1. )
      : aol::FELinScalarWeightedMassInterface<ConfiguratorType, WeightedMassOp<ConfiguratorType, IndexMode>, IndexMode > ( Config, Initializer, OpType ),
      _discrImg ( NULL ), _epsSqr ( Epsilon*Epsilon )  {}

  ~WeightedMassOp( ) {
    if ( _discrImg ) delete _discrImg;
  }

  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( _discrImg )
      delete _discrImg;
    _discrImg = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator(), Image );
    this->reset();
  }

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType& ) const {
    if ( !_discrImg ) {
      throw aol::Exception ( "WeightedMassOp::getCoeff: set image first!", __FILE__, __LINE__ );
    }
    typename ConfiguratorType::VecType grad, v;
    _discrImg->evaluateGradientAtQuadPoint ( El, QuadPoint, grad );
    return sqrt ( grad.normSqr() + _epsSqr );
  }
};


/** \brief This class computes the following weighted mass op:
 *         \f[ \int_{\Omega} G \| \nabla \Phi \| \vartheta_i \vartheta_j \, dx \f]
 *         where G is an additional FE-function
 *
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class WeightedAdditionalFEFctMassOp : public aol::FELinScalarWeightedMassInterface < ConfiguratorType,
      WeightedAdditionalFEFctMassOp<ConfiguratorType, IndexMode>, IndexMode > {
public:
  typedef typename ConfiguratorType::RealType RealType;

protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrPhi;    // for computing \| \nabla \Phi \|
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrFEFct;  // the additional FE-function
  RealType _epsSqr;

public:
  WeightedAdditionalFEFctMassOp ( const typename ConfiguratorType::InitType &Initializer,
                                  aol::OperatorType OpType = aol::ONTHEFLY, RealType Epsilon = 1. )
      : aol::FELinScalarWeightedMassInterface<ConfiguratorType, WeightedAdditionalFEFctMassOp<ConfiguratorType, IndexMode>, IndexMode > ( Initializer, OpType ),
      _discrPhi ( NULL ), _discrFEFct ( NULL ), _epsSqr ( Epsilon*Epsilon )  {}

  ~WeightedAdditionalFEFctMassOp( ) {
    if ( _discrPhi ) delete _discrPhi;
    if ( _discrFEFct ) delete _discrFEFct;
  }

  void setPhiReference ( const aol::Vector<RealType> &Image ) {
    _discrPhi = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator(), Image );
  }

  void setFEFctReference ( const aol::Vector<RealType> &Image ) {
    _discrFEFct = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator(), Image );
  }

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::VecType& ) const {
    if ( !_discrPhi ) {
      throw aol::Exception ( "WeightedAdditionalFEFctMassOp::getCoeff: set Phi first!", __FILE__, __LINE__ );
    }
    if ( !_discrFEFct ) {
      throw aol::Exception ( "WeightedAdditionalFEFctMassOp::getCoeff: set additional FEFct first!", __FILE__, __LINE__ );
    }
    typename ConfiguratorType::VecType grad, v;
    _discrPhi->evaluateGradientAtQuadPoint ( El, QuadPoint, grad );

    RealType G = _discrFEFct->evaluateAtQuadPoint ( El, QuadPoint );

    return G * sqrt ( grad.normSqr() + _epsSqr );
  }
};


/** \brief This class computes the following op which is the weak formulation of the LaplaceBeltrami-operator:
 *         \f[ \int_{\Omega} \| \nabla \Phi \| P[\Phi] \nabla \vartheta_i \nabla \vartheta_j \, dx \f]
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class LaplaceBeltramiOp : public aol::FELinMatrixWeightedStiffInterface < ConfiguratorType,
      LaplaceBeltramiOp<ConfiguratorType, IndexMode>, IndexMode > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrImg;
  const RealType _epsSqr;

public:
  LaplaceBeltramiOp ( const typename ConfiguratorType::InitType &Initializer,
                      aol::OperatorType OpType, RealType Epsilon )
      : aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, LaplaceBeltramiOp<ConfiguratorType, IndexMode>, IndexMode >
      ( Initializer, OpType ), _discrImg ( NULL ), _epsSqr ( Epsilon*Epsilon )  {  }

  LaplaceBeltramiOp ( const ConfiguratorType & Config, const typename ConfiguratorType::InitType &Initializer,
                      aol::OperatorType OpType, RealType Epsilon )
      : aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, LaplaceBeltramiOp<ConfiguratorType, IndexMode>, IndexMode >
      ( Config, Initializer, OpType ), _discrImg ( NULL ), _epsSqr ( Epsilon*Epsilon )  {  }

  ~LaplaceBeltramiOp( ) {
    if ( _discrImg ) delete _discrImg;
  }

  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( static_cast< int > ( Image.size() ) != this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "Array.size() = " << Image.size() << " getConfigurator().getNumGlobalDofs() = " << this->getConfigurator().getNumGlobalDofs()  << endl;
      throw aol::Exception ( "array phi has wrong size\n", __FILE__, __LINE__ );
    }
    if ( _discrImg )
      delete _discrImg;
    _discrImg = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator(), Image );
    this->reset();
  }

  // implement the projection-matrix:
  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                               const typename ConfiguratorType::DomVecType& /*RefCoord*/, typename ConfiguratorType::MatType &mat ) const {
    if ( !_discrImg ) {
      throw aol::Exception ( "LaplaceBeltramiOp::getCoeffMatrix: set image first!", __FILE__, __LINE__ );
    }

    typename ConfiguratorType::VecType grad;
    _discrImg->evaluateGradientAtQuadPoint ( El, QuadPoint, grad );

    RealType norm = sqrt ( grad.normSqr() + _epsSqr );
    grad /= norm;

    // compute the projection matrix
    for ( int i = 0; i < ConfiguratorType::VecType::dim; i++ )
      for ( int j = 0; j < ConfiguratorType::VecType::dim; j++ )
        mat[i][j] = -grad[i] * grad[j];
    for ( int i = 0; i < ConfiguratorType::VecType::dim; i++ )
      mat[i][i] += 1.;

    mat *= norm;    // \| \nabla \Phi \|
  }

};




/** This class computes an anisotropic stiffness matrix. The gradients
 * are projected onto the tangent space and weighted different in directions \f$ v,v^{\perp}\f$
 * These directions are given by a vector-field which is projected onto the tangent space.
 * The complete expression is:
 * \f[ \int_{\Omega} (\alpha^2 P[v] + \beta^2 P[v^{\perp}]) \nabla \Phi \cdot \nabla \Psi \, dx \f]
 * \attention This operator is only available in 3D!
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class AnisoTangentialProjectionStiffOp : public aol::FELinMatrixWeightedStiffInterface < ConfiguratorType,
      AnisoTangentialProjectionStiffOp<ConfiguratorType, IndexMode>, IndexMode > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrPhi;    // the levelset-function
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrAlpha;  // weight in the main smoothing direction
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrBeta;   // weight in direction orthogonal to main dir. but in TxM
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType, ConfiguratorType::Dim> *_discrVectorField;
//   RealType _alphaSqr, _betaSqr;                                       // The anisotropic weights of the tangent directions
  const RealType _epsSqr;
  const RealType _zeroBound;                                          // bound to test whether a vector is zero or not

public:
  AnisoTangentialProjectionStiffOp (  const typename ConfiguratorType::InitType &Initializer,
                                      aol::OperatorType OpType, RealType Epsilon )
      : aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, AnisoTangentialProjectionStiffOp<ConfiguratorType, IndexMode>, IndexMode >
      ( Initializer, OpType ), _discrPhi ( NULL ), _discrAlpha ( NULL ), _discrBeta( NULL ), _discrVectorField ( NULL ),
      _epsSqr ( Epsilon*Epsilon ), _zeroBound ( _epsSqr )  {
    if ( ConfiguratorType::Dim != 3 )
      throw aol::Exception ( "AnisoTangentialProjectionStiffOp only available in 3D!", __FILE__, __LINE__ );
//     _alphaSqr  = 1.;
//     _betaSqr   = 1.;
  }

  ~AnisoTangentialProjectionStiffOp( ) {
    if ( _discrPhi )          delete _discrPhi;                       // necessary to compute the normal
    if ( _discrAlpha )        delete _discrAlpha;                     // weight in the main smoothing direction
    if ( _discrBeta )         delete _discrBeta;                      // weight in direction orthogonal to main dir. but in TxM
    if ( _discrVectorField )  delete _discrVectorField;               // the given main directions
  }

//   void setAlphaBeta ( const RealType Alpha, const RealType Beta ) {
//     _alphaSqr = Alpha * Alpha;
//     _betaSqr  = Beta * Beta;
//   }

  void setVectorFieldReference ( const aol::MultiVector<RealType> &VectorField ) {
    _discrVectorField = new aol::DiscreteVectorFunctionDefault<ConfiguratorType, ConfiguratorType::Dim> ( this->getConfigurator().getInitializer(), VectorField );
  }

  // level set function, necessary to compute the normal on the surface, this is necessary to compute v x n
  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( static_cast< int > ( Image.size() ) != this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "Array.size() = " << Image.size() << " getConfigurator().getNumGlobalDofs() = " << this->getConfigurator().getNumGlobalDofs()  << endl;
      throw aol::Exception ( "array phi has wrong size\n", __FILE__, __LINE__ );
    }
    if ( _discrPhi )
      delete _discrPhi;
    _discrPhi = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator(), Image );
    this->reset();
  }

  void setWeightReferences ( const aol::Vector<RealType> &Alpha, const aol::Vector<RealType> &Beta ) {
    if (    static_cast< int > ( Alpha.size() ) != this->getConfigurator().getNumGlobalDofs()
         || static_cast< int > ( Beta.size() )  != this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "Alpha.size(): "<< Alpha.size() << ", Beta.size(): " << Beta.size() << endl;
      cerr << "Weight.size() (Alpha or Beta) != config.getNumGlobalDofs() = " << this->getConfigurator().getNumGlobalDofs()  << endl;
      throw aol::Exception ( "array Alpha or Beta has wrong size\n", __FILE__, __LINE__ );
    }
    _discrAlpha = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator(), Alpha );
    _discrBeta  = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator(), Beta );
  }

  // create the projection matrix onto the tangent space
  inline void addProjectionMatrix ( const typename ConfiguratorType::VecType& v, typename ConfiguratorType::MatType &mat ) const {
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        mat[i][j] += v[i] * v[j];
  }

  // implement the projection-matrix: ONLY available in 3D!!
  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                               const typename ConfiguratorType::VecType& /*RefCoord*/, typename ConfiguratorType::MatType &mat ) const {
    if ( !_discrPhi )
      throw aol::Exception ( "AnisoTangentialProjectionStiffOp::getCoeffMatrix: set levelset-function Phi first!", __FILE__, __LINE__ );
    if ( !_discrVectorField )
      throw aol::Exception ( "AnisoTangentialProjectionStiffOp::getCoeffMatrix: set vectorField first!", __FILE__, __LINE__ );

    // first project the vector field onto the tangent space
    typename ConfiguratorType::DomVecType gradPhi;              // normal on the surface
    typename ConfiguratorType::DomVecType v;                    // the vector from the vector-field (main smoothing direction)
    typename ConfiguratorType::DomVecType vTangent;             // vectorfield projected on the tangent space
    typename ConfiguratorType::DomVecType vOrthogonal;          // vectorfield projected on the tangent space but orthogonal to v

    // evaluate the normal and the vector field and normalize it (epsilon-norm)
    _discrPhi->evaluateGradientAtQuadPoint ( El, QuadPoint, gradPhi );
    double normGradPhi = sqrt ( gradPhi.normSqr() + _epsSqr );
    gradPhi /= normGradPhi;
    _discrVectorField->evaluateAtQuadPoint ( El, QuadPoint, v );
    v /= sqrt ( v.normSqr() + _epsSqr );

    // evaluate the anisotropic weights alpha and beta
    RealType alphaSqr = aol::Sqr( _discrAlpha->evaluateAtQuadPoint ( El, QuadPoint ) );
    RealType betaSqr  = aol::Sqr( _discrBeta->evaluateAtQuadPoint  ( El, QuadPoint ) );

    // declare the projection matrices
    typename ConfiguratorType::MatType Pn;                      // projection on the tangent space (Id - n x n)
    typename ConfiguratorType::MatType PvTangent;               // projection on v (v x v)
    typename ConfiguratorType::MatType PvOrthogonal;            // projection on v^{\perp}

    // compute the projection matrix that maps on the tangent space
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        Pn[i][j] = - gradPhi[i] * gradPhi[j];
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      Pn[i][i] += 1.;

    Pn.mult ( v, vTangent );                                    // project v to the tangent space
    // if vTangent is zero define it to be arbitrary in the tangent space
    if ( vTangent.norm() < _zeroBound ) {
//       cerr << "AnisoTangentialProjectionStiffOp::getCoeffMatrix: vTangent is nearly zero!\n";
//       cerr << gradPhi << ", " << v << endl;
      vTangent.setZero();
      if ( abs ( gradPhi[0] ) > _zeroBound ) {                  // use components that are big enough
        if ( abs ( gradPhi[1] ) > _zeroBound ) {
          vTangent[1] = -gradPhi[0];
          vTangent[0] = gradPhi[1];
          vTangent.normalize();
        } else {
          vTangent[1] = 1.;
        }
      } else {
        vTangent[0] = 1.;
      }
    }

    vTangent /= vTangent.norm();                                // normalize vTangent

    // now define v^{\perp} by the cross product of n and v
    // if vTangent and gradPhi are normalized, then the cross-product is normalized too.
    vOrthogonal[0] = vTangent[1] * gradPhi[2] - vTangent[2] * gradPhi[1];
    vOrthogonal[1] = vTangent[2] * gradPhi[0] - vTangent[0] * gradPhi[2];
    vOrthogonal[2] = vTangent[0] * gradPhi[1] - vTangent[1] * gradPhi[0];

    // define the projection on v and v^{\perp}
    mat.setZero();
    addProjectionMatrix ( vTangent, mat );
    mat *= alphaSqr / betaSqr;
    addProjectionMatrix ( vOrthogonal, mat );
    mat *= betaSqr;

  }

};

/** This class computes the following operator:
 * \f[ \int_{\Omega} \| \nabla \Phi \| \vartheta \, dx \f]
 */
template <typename ConfiguratorType>
class GradPhiTestFctOp : public aol::FENonlinOpInterface< ConfiguratorType, GradPhiTestFctOp<ConfiguratorType> > {
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrImg;
  const RealType _epsSqr;

public:
  GradPhiTestFctOp ( const typename ConfiguratorType::InitType &Initializer, RealType Epsilon  )
      : aol::FENonlinOpInterface< ConfiguratorType, GradPhiTestFctOp<ConfiguratorType> > ( Initializer ),
      _discrImg ( NULL ), _epsSqr ( Epsilon*Epsilon ) { }

  ~GradPhiTestFctOp( ) {
    if ( _discrImg ) delete _discrImg;
  }

  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( static_cast< int > ( Image.size() ) != this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "Array.size() = " << Image.size() << " getConfigurator().getNumGlobalDofs() = " << this->getConfigurator().getNumGlobalDofs()  << endl;
      throw aol::Exception ( "array phi has wrong size\n", __FILE__, __LINE__ );
    }
    if ( _discrImg )
      delete _discrImg;
    _discrImg = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator(), Image );
    this->reset();
  }

  // just evaluate \f$ \| \nabla \Phi \| \f$
  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &/*DiscFunc*/,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
                         RealType &NL ) const {
    if ( !_discrImg ) {
      throw aol::Exception ( "GradPhiTestFctOp::getNonlinearity: set image first!", __FILE__, __LINE__ );
    }
    typename ConfiguratorType::DomVecType grad;
    _discrImg->evaluateGradientAtQuadPoint ( El, QuadPoint, grad );

    NL = sqrt ( grad.normSqr() + _epsSqr );
  }
};




// ------------------------------ BOUNDARY OPERATORS ---------------------------------------------------

/**
 * This class computes the following boundary integral:
 * \f[ \int_{\partial\Omega} \nabla u \cdot n \, \varphi da \f]
 * ATTENTION: ONLY WORKING IN 2D!
 *
 */
template <typename ConfiguratorType, typename BoundaryQuadratureType>
class IntegrateProjectionUOverBoundary
      : public BoundaryIntegrationInterface < ConfiguratorType,
      aol::Vector<typename ConfiguratorType::RealType>,
      ConfiguratorType::Dim,
      BoundaryQuadratureType, IntegrateProjectionUOverBoundary< ConfiguratorType, BoundaryQuadratureType > > {
  typename ConfiguratorType::RealType _epsSqr;
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrPhi;

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

public:
  IntegrateProjectionUOverBoundary ( const typename ConfiguratorType::InitType &Initializer, RealType Epsilon )
      : BoundaryIntegrationInterface < ConfiguratorType,
      aol::Vector<RealType>,
      ConfiguratorType::Dim,
      BoundaryQuadratureType,
      IntegrateProjectionUOverBoundary< ConfiguratorType, BoundaryQuadratureType > > ( Initializer ),
      _epsSqr ( Epsilon*Epsilon ), _discrPhi ( NULL ) {}

  // set the reference to phi, the level set function where we project to
  void setImageReference ( const aol::Vector<RealType> &Phi ) {
    if ( static_cast< int > ( Phi.size() ) != this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "Array.size() = " << Phi.size() << " getConfigurator().getNumGlobalDofs() = " << this->getConfigurator().getNumGlobalDofs()  << endl;
      throw aol::Exception ( "IntegrateProjectionUOverBoundary: Array Phi has wrong size\n", __FILE__, __LINE__ );
    }
    if ( _discrPhi )
      delete _discrPhi;
    _discrPhi = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator(), Phi );
    this->reset();
  }

  //! computes \f[ P[\Phi] \nabla U \f]
  RealType integrandAtQuadPoint ( const aol::Vector<RealType> &arg, const typename ConfiguratorType::VecType &normal,
                                  const ElementType &el, const typename ConfiguratorType::VecType edgePoint ) const {
    if ( !_discrPhi ) {
      throw aol::Exception ( "IntegrateProjectionUOverBoundary::integrandAtQuadPoint: set image first!", __FILE__, __LINE__ );
    }

    // compute gradient of boundaryfunction @ edgePt
    aol::DiscreteFunctionDefault<ConfiguratorType> discrU ( this->getConfigurator(), arg );
    VecType gradPhi, gradU, PGradU;
    discrU.evaluateGradient ( el, edgePoint, gradU );
    _discrPhi->evaluateGradient ( el, edgePoint, gradPhi );

    RealType norm = sqrt ( gradPhi.normSqr() + _epsSqr );       // sqrt( nabla Phi^2 + eps^2 )
    gradPhi /= norm;

    MatType mat;

    // compute the projection: Id - (nabla Phi)/||...|| \times (nabla Phi)/||...||
    for ( int i = 0; i < ConfiguratorType::VecType::dim; i++ )
      for ( int j = 0; j < ConfiguratorType::VecType::dim; j++ )
        mat[i][j] = -gradPhi[i] * gradPhi[j];
    for ( int i = 0; i < ConfiguratorType::VecType::dim; i++ )
      mat[i][i] += 1.;

    mat.mult ( gradU, PGradU );

    return PGradU * normal;
  }
};



}   // end namespace

#endif
