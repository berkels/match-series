#ifndef __WILLMORE_H
#define __WILLMORE_H

#include <scalarArray.h>
#include <gridOp.h>
#include <FEOpInterface.h>
#include <anisotropies.h>
#include <solver.h>
#include <discreteFunction.h>
#include <boundaryIntegration.h>
#include <mcm.h>

namespace qc {

// ---------------------------------------  INTEGRATION OPERATORS ON THE WHOLE DOMAIN --------------------------------------------

//! \brief operator from willmore flow with p-regularization.
//! \ingroup MatrixFEOp
template <typename ConfiguratorType, typename AnisotropyType>
class WillmoreStiffExpOp : public aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, WillmoreStiffExpOp<ConfiguratorType, AnisotropyType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;

protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrImg;
  RealType _eps;
  const AnisotropyType &_anisotropy;

public:
  WillmoreStiffExpOp ( const typename ConfiguratorType::InitType &Initializer, aol::OperatorType OpType, RealType Epsilon, const AnisotropyType &Anisotropy )
      : aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, WillmoreStiffExpOp<ConfiguratorType, AnisotropyType> > ( Initializer, OpType ), _discrImg ( NULL ), _eps ( Epsilon ), _anisotropy ( Anisotropy )  {}

  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( static_cast< int > ( Image.size() ) != this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "Array.size() = " << Image.size() << " getConfigurator().getNumGlobalDofs() = " << this->getConfigurator().getNumGlobalDofs()  << endl;
      throw aol::Exception ( "array phi has wrong size\n", __FILE__, __LINE__ );
    }
    if ( _discrImg )
      delete _discrImg;
    _discrImg = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator().getInitializer(), Image );
    this->reset();
  }

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::VecType& /*RefCoord*/, typename ConfiguratorType::MatType &mat ) const {
    if ( !_discrImg ) throw aol::Exception ( "qc::WillmoreStiffExpOp: Set Image first!\n", __FILE__, __LINE__ );
    typename ConfiguratorType::DomVecType grad;
    _discrImg->evaluateGradientAtQuadPoint ( El, QuadPoint, grad );

    _anisotropy.explicitPart ( grad, mat );
  }

};


/**
 * \brief implementation of \f[ \int_{\Omega} \frac{P[\nabla \Phi]}{|\nabla \Phi|} \nabla \psi_i \nabla \psi_j \, dx \f],
 *        with \f$ P[\nabla \Phi] = Id - \frac{\nabla \Phi}{|\nabla \Phi|} \otimes \frac{\nabla \Phi}{|\nabla \Phi|} \f$.
 * \attention Just for the ISOTROPIC willmore flow!
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType>
class WillmoreProjectionStiffOp : public aol::FELinMatrixWeightedStiffInterface < ConfiguratorType,
      WillmoreProjectionStiffOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrImg;
  RealType _eps;

public:
  WillmoreProjectionStiffOp ( const typename ConfiguratorType::InitType &Initializer,
                                aol::OperatorType OpType, RealType Epsilon )
      : aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, WillmoreProjectionStiffOp<ConfiguratorType> >
      ( Initializer, OpType ), _discrImg ( NULL ), _eps ( Epsilon ) {  }

  ~WillmoreProjectionStiffOp( ) {
    if ( _discrImg )
      delete _discrImg;
  }

  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( static_cast< int > ( Image.size() ) != this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "Array.size() = " << Image.size() << " getConfigurator().getNumGlobalDofs() = " << this->getConfigurator().getNumGlobalDofs()  << endl;
      throw aol::Exception ( "array phi has wrong size\n", __FILE__, __LINE__ );
    }
    if ( _discrImg )
      delete _discrImg;
    _discrImg = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator().getInitializer(), Image );
    this->reset();
  }

  // implement the projection-matrix:
  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                               const typename ConfiguratorType::VecType& /*RefCoord*/, typename ConfiguratorType::MatType &mat ) const {

    typename ConfiguratorType::DomVecType grad;
    _discrImg->evaluateGradientAtQuadPoint ( El, QuadPoint, grad );

    RealType norm = sqrt ( grad.normSqr() + aol::Sqr ( _eps ) );
    grad /= norm;

    for ( int i = 0; i < ConfiguratorType::VecType::dim; i++ )
      for ( int j = 0; j < ConfiguratorType::VecType::dim; j++ )
        mat[i][j] = -grad[i] * grad[j];
    for ( int i = 0; i < ConfiguratorType::VecType::dim; i++ )
      mat[i][i] += 1.;

    mat /= norm;
  }

};


/**
 * \brief implementation of \f[ \int_{\Omega} ( \gamma_{zz}(\nabla \Phi) \nabla \psi_i \nabla \psi_j \, dx \f],
 * \attention For the ANISOTROPIC willmore flow!
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, typename AnisotropyType>
class WillmoreStiffImpOp : public aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, WillmoreStiffImpOp<ConfiguratorType, AnisotropyType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrImg;
  const AnisotropyType &_anisotropy;
public:

  WillmoreStiffImpOp ( const typename ConfiguratorType::InitType &Initializer, aol::OperatorType OpType, const AnisotropyType &Anisotropy )
      : aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, WillmoreStiffImpOp<ConfiguratorType, AnisotropyType> > ( Initializer, OpType ),
      _discrImg ( NULL ), _anisotropy ( Anisotropy )  {}


  //! set image (by reference!)
  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( _discrImg )
      delete _discrImg;
    _discrImg = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator().getInitializer(), Image );
    this->reset();
  }

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::VecType& /*RefCoord*/, typename ConfiguratorType::MatType &mat ) const {
    if ( !_discrImg ) throw aol::Exception ( "qc::WillmoreStiffImpOp: Set Image first!\n", __FILE__, __LINE__ );
    typename ConfiguratorType::DomVecType grad;
    _discrImg->evaluateGradientAtQuadPoint ( El, QuadPoint, grad );

    _anisotropy.implicitPart ( grad, mat );

  }
};

/**
 * \brief implementation of \f[ \int_{\Omega} ( \gamma_{zz}(\nabla \Phi) \nabla \psi_i \nabla \psi_j \, dx \f],
 *        Difference to non-general version is that the anisotropy is called with an element and a quadpoint.
 * \attention For the ANISOTROPIC willmore flow!
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, typename AnisotropyType>
class GeneralWillmoreStiffImpOp : public aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, GeneralWillmoreStiffImpOp<ConfiguratorType, AnisotropyType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrImg;
  const AnisotropyType &_anisotropy;
public:

  GeneralWillmoreStiffImpOp ( const typename ConfiguratorType::InitType &Initializer,
                                aol::OperatorType OpType,
                                const AnisotropyType &Anisotropy )
    : aol::FELinMatrixWeightedStiffInterface<ConfiguratorType,
                              GeneralWillmoreStiffImpOp<ConfiguratorType, AnisotropyType> > ( Initializer, OpType ),
      _discrImg ( NULL ), _anisotropy ( Anisotropy )  {}


  //! set image (by reference!)
  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( _discrImg )
      delete _discrImg;
    _discrImg = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator().getInitializer(), Image );
    this->reset();
  }

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::VecType& /*RefCoord*/, typename ConfiguratorType::MatType &mat ) const {
    if ( !_discrImg ) throw aol::Exception ( "qc::WillmoreStiffImpOp: Set Image first!\n", __FILE__, __LINE__ );
    typename ConfiguratorType::DomVecType grad;
    _discrImg->evaluateGradientAtQuadPoint ( El, QuadPoint, grad );

    _anisotropy.implicitPart ( El, QuadPoint, grad, mat );

  }
};



//! \brief operator from willmore flow with p-regularization.
//! \ingroup MatrixFEOp
template <typename ConfiguratorType, typename AnisotropyType>
class WillmoreStiffAnisoOp : public aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, WillmoreStiffAnisoOp<ConfiguratorType, AnisotropyType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrImg;
  RealType _eps;
  const AnisotropyType &_anisotropy;
public:

  WillmoreStiffAnisoOp ( const typename ConfiguratorType::InitType &Initializer, aol::OperatorType OpType, RealType Epsilon, const AnisotropyType &Anisotropy )
      : aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, WillmoreStiffAnisoOp<ConfiguratorType, AnisotropyType> > ( Initializer, OpType ), _discrImg ( NULL ), _eps ( Epsilon ), _anisotropy ( Anisotropy )  {}


  //! set image (by reference!)
  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( _discrImg )
      delete _discrImg;
    _discrImg = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator().getInitializer(), Image );
    this->reset();
  }

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::VecType& /*RefCoord*/, typename ConfiguratorType::MatType &mat ) const {
    if ( !_discrImg ) throw aol::Exception ( "qc::WillmoreStiffAnisoOp: Set Image first!\n", __FILE__, __LINE__ );
    typename ConfiguratorType::DomVecType grad;
    _discrImg->evaluateGradientAtQuadPoint ( El, QuadPoint, grad );

    typename ConfiguratorType::MatType matExp;
    _anisotropy.implicitPart ( grad, mat );
    _anisotropy.explicitPart ( grad, matExp );
    mat -= matExp;
  }

};

/**
 * \brief implementation of \f[ \int_{\Omega} \frac12 \frac{W^2}{\| \nabla \phi \|^3_{\epsilon} } \nabla \psi_i \nabla \psi_j \, dx \f]
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType>
class WillmoreStiffWOp : public aol::FELinScalarWeightedStiffInterface<ConfiguratorType, WillmoreStiffWOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;

  WillmoreStiffWOp ( const typename ConfiguratorType::InitType &Initializer, aol::OperatorType OpType, RealType Epsilon )
      : aol::FELinScalarWeightedStiffInterface<ConfiguratorType, WillmoreStiffWOp<ConfiguratorType> > ( Initializer, OpType ),
      _discrImg ( NULL ), _discrW ( NULL ), _epsSqr ( Epsilon*Epsilon ) {}

  //! set reference to the image
  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( Image.size() != this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "array.size() = " << Image.size() << " getConfigurator().getNumGlobalDofs() = " <<  this->getConfigurator().getNumGlobalDofs() << endl;
      throw aol::Exception ( "array phi has wrong size\n", __FILE__, __LINE__ );
    }
    if ( _discrImg )
      delete _discrImg;
    _discrImg = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator().getInitializer(), Image );
    this->reset();
  }
  //! set reference to the curvature concentration W
  void setWReference ( const aol::Vector<RealType> &W ) {
    if ( W.size() != this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "array.size() = " << W.size() << " getConfigurator().getNumGlobalDofs() = " <<  this->getConfigurator().getNumGlobalDofs() << endl;
      throw aol::Exception ( "array W has wrong size\n", __FILE__, __LINE__ );
    }
    if ( _discrW )
      delete _discrW;
    _discrW = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator().getInitializer(), W );
    this->reset();
  }

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::VecType& /*RefCoord*/ ) const {
    if ( !_discrW || !_discrImg ) throw aol::Exception ( "qc::WillmoreStiffWOp: Set Image and W first!\n", __FILE__, __LINE__ );
    typename ConfiguratorType::DomVecType grad;
    _discrImg->evaluateGradientAtQuadPoint ( El, QuadPoint, grad );

    return 0.5 * pow ( grad.normSqr() + _epsSqr, -1.5 ) * aol::Sqr ( _discrW->evaluateAtQuadPoint ( El, QuadPoint ) );
  }

protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrImg;
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrW;
  RealType _epsSqr;
};




/** ********************************************************************************************
// qc::AnisotropyIntegrationOp
// implementation of \f[ \int_{\Omega} \gamma_z(\nabla \phi) \nabla \psi_j \, dx \f],
// ******************************************************************************************** */
template <typename ConfiguratorType, typename AnisoType>
class AnisotropyIntegrationOp : public aol::FENonlinDiffOpInterface< ConfiguratorType, AnisotropyIntegrationOp<ConfiguratorType, AnisoType> > {
  const AnisoType &_anisotropy;
public:
  typedef typename ConfiguratorType::RealType RealType;

  AnisotropyIntegrationOp ( const typename ConfiguratorType::InitType &Initializer,
                              const AnisoType &Anisotropy ) :
      aol::FENonlinDiffOpInterface < ConfiguratorType,
      AnisotropyIntegrationOp<ConfiguratorType, AnisoType> > ( Initializer ),
      _anisotropy ( Anisotropy ) {}


  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
                         typename ConfiguratorType::VecType &NL ) const {
    typename ConfiguratorType::VecType Gradient;
    DiscFunc.evaluateGradientAtQuadPoint ( El, QuadPoint, Gradient );

    _anisotropy.gammaFirstDerivative ( Gradient, NL );
  }
};

/** ********************************************************************************************
// qc::AnisotropyIntegrationOp
// implementation of \f[ \int_{\Omega} \gamma_z(\nabla \phi) \nabla \psi_j \, dx \f],
// Difference to non-general version is that the anisotropy is called with an element and a quadpoint.
// ******************************************************************************************** */
template <typename ConfiguratorType, typename AnisoType>
class GeneralAnisotropyIntegrationOp : public aol::FENonlinDiffOpInterface< ConfiguratorType, GeneralAnisotropyIntegrationOp<ConfiguratorType, AnisoType> > {
  const AnisoType &_anisotropy;
public:
  typedef typename ConfiguratorType::RealType RealType;

  GeneralAnisotropyIntegrationOp ( const typename ConfiguratorType::InitType &Initializer,
                                     const AnisoType &Anisotropy ) :
      aol::FENonlinDiffOpInterface < ConfiguratorType,
      GeneralAnisotropyIntegrationOp<ConfiguratorType, AnisoType> > ( Initializer ),
      _anisotropy ( Anisotropy ) {}


  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
                         typename ConfiguratorType::VecType &NL ) const {
    typename ConfiguratorType::VecType Gradient;
    DiscFunc.evaluateGradientAtQuadPoint ( El, QuadPoint, Gradient );

    _anisotropy.gammaFirstDerivative ( El, QuadPoint, Gradient, NL );
  }
};


/** ********************************************************************************************
// qc::AnisotropyMassIntegrationOp
// implementation of \f[ \int_{\Omega} \gamma(\nabla \phi) \psi_j \, dx \f],
// ******************************************************************************************** */
template <typename ConfiguratorType, typename AnisoType>
class AnisotropyMassIntegrationOp : public aol::FENonlinOpInterface< ConfiguratorType, AnisotropyMassIntegrationOp<ConfiguratorType, AnisoType> > {
  const AnisoType &_anisotropy;
public:
  typedef typename ConfiguratorType::RealType RealType;

  AnisotropyMassIntegrationOp ( const typename ConfiguratorType::InitType &Initializer,
                                  const AnisoType &Anisotropy ) :
      aol::FENonlinOpInterface < ConfiguratorType,
      AnisotropyMassIntegrationOp<ConfiguratorType, AnisoType> > ( Initializer ),
      _anisotropy ( Anisotropy ) {}


  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
                         typename ConfiguratorType::RealType &NL ) const {
    typename ConfiguratorType::VecType Gradient;
    DiscFunc.evaluateGradientAtQuadPoint ( El, QuadPoint, Gradient );

    // Gradient normieren
    if ( Gradient.norm() != 0. )  Gradient /= Gradient.norm();

    NL = _anisotropy.gammaNorm ( Gradient );
  }
};



/** ------------------ qcLoaclRotatedAnisotropyIntegrationOp -------------------------------------
// anisotropy depends local on the element, a vector-field with the rotations-directions
// of the anisotropy must be given. Before evaluating the NonLinearity, the anisotropy
// is rotated by R into the direction given from the vector-field, so it's the implementation of
// \f[ \int_{\Omega} \gamma_z( R \nabla \phi) \nabla \psi_j \, dx \f],
// ---------------------------------------------------------------------------------------------- */

template <typename ConfiguratorType, typename AnisoType>
class LocalRotatedAnisotropyIntegrationOp : public aol::FENonlinDiffOpInterface < ConfiguratorType,
      LocalRotatedAnisotropyIntegrationOp<ConfiguratorType, AnisoType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VecType VecType;

private:
  AnisoType &_anisotropy;
  const ScalarArray<RealType, qc::QC_3D>* _rotDirectx;     // the components of the vector-field
  const ScalarArray<RealType, qc::QC_3D>* _rotDirecty;
  const ScalarArray<RealType, qc::QC_3D>* _rotDirectz;

public:
  // constructor
  LocalRotatedAnisotropyIntegrationOp ( const typename ConfiguratorType::InitType &Initializer,
                                          AnisoType &Anisotropy,
                                          const ScalarArray<RealType, qc::QC_3D> *rotDx,
                                          const ScalarArray<RealType, qc::QC_3D> *rotDy,
                                          const ScalarArray<RealType, qc::QC_3D> *rotDz ) :
      aol::FENonlinDiffOpInterface < ConfiguratorType,
      LocalRotatedAnisotropyIntegrationOp<ConfiguratorType, AnisoType> > ( Initializer ),
      _anisotropy ( Anisotropy ) {
    _rotDirectx = rotDx;
    _rotDirecty = rotDy;
    _rotDirectz = rotDz;
  }


  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const VecType &RefCoord,
                         VecType &NL ) const {
    VecType Gradient, coord;
    DiscFunc.evaluateGradientAtQuadPoint ( El, QuadPoint, Gradient );

    coord[0] = El.x() + RefCoord[0];
    coord[1] = El.y() + RefCoord[1];
    if ( ConfiguratorType::VecType::dim == 3 ) {
      coord[2] = El.z() + RefCoord[2];
    }

    // now, before evaluating the derivative of the anisotropy, rotate it
    // into the direction given from the vector-field:
    aol::Vec3<RealType> rotDirect;

    rotDirect[0] = _rotDirectx->interpolate ( coord[0], coord[1], coord[2] );
    rotDirect[1] = _rotDirecty->interpolate ( coord[0], coord[1], coord[2] );
    rotDirect[2] = _rotDirectz->interpolate ( coord[0], coord[1], coord[2] );

    _anisotropy.setRotateDirection ( rotDirect );

    _anisotropy.gammaFirstDerivative ( coord, Gradient, NL );
  }
};


//! Operator for computing \f$ \int_{\Omega} \gamma(\nabla u) \, dx. \f$
/*!
 * \author Nemitz
*/
template <typename ConfiguratorType, typename AnisoType>
class GammaIntegrationOp : public aol::FENonlinIntegrationVectorGeneralInterface < ConfiguratorType,
      GammaIntegrationOp<ConfiguratorType, AnisoType> > {
  const AnisoType &_anisotropy;
public:
  typedef typename ConfiguratorType::RealType RealType;

  GammaIntegrationOp ( const typename ConfiguratorType::InitType &Initializer,
                         const AnisoType &Anisotropy ) :
      aol::FENonlinIntegrationVectorGeneralInterface < ConfiguratorType,
      GammaIntegrationOp<ConfiguratorType, AnisoType> > ( Initializer ),
      _anisotropy ( Anisotropy ) {}

  RealType evaluateIntegrand ( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFunc,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/ ) const {
    typename ConfiguratorType::VecType Gradient;
    DiscFunc[0].evaluateGradientAtQuadPoint ( El, QuadPoint, Gradient );

    return _anisotropy.gammaNorm ( Gradient );
  }
};


// ---------------------------------------  BOUNDARY INTEGRATION OPERATORS --------------------------------------------

/**
 * This class computes the following boundary integral:
 * \f[ \int_{\partial\Omega} \gamma_z( \nabla u ) \cdot n \, \varphi da \f]
 */
template <typename ConfiguratorType, typename BoundaryQuadratureType, typename AnisoType>
class IntegrateGammaZOverBoundary
      : public BoundaryIntegrationInterface < ConfiguratorType,
      aol::Vector<typename ConfiguratorType::RealType>,
      ConfiguratorType::Dim,
      BoundaryQuadratureType,
      IntegrateGammaZOverBoundary< ConfiguratorType, BoundaryQuadratureType, AnisoType > > {

  const AnisoType &_anisotropy;

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::VecType VecType;

public:
  IntegrateGammaZOverBoundary ( const typename ConfiguratorType::InitType &Initializer, const AnisoType &Anisotropy  )
      : BoundaryIntegrationInterface < ConfiguratorType,
      aol::Vector<RealType>,
      ConfiguratorType::Dim,
      BoundaryQuadratureType,
      IntegrateGammaZOverBoundary< ConfiguratorType, BoundaryQuadratureType, AnisoType > > ( Initializer ),
      _anisotropy ( Anisotropy ) {  }

  //! evaluate \f[ \gamma_z( \nabla U ) \cdot n \f]
  RealType integrandAtQuadPoint ( const aol::Vector<RealType> &arg, const VecType &normal, const ElementType &el, const VecType edgePoint ) const {
    aol::DiscreteFunctionDefault<ConfiguratorType> discrU ( this->_grid, arg );

    // compute gradient of boundaryfunction @ edgePt and gammaZ ( grad U )
    VecType gradU;
    VecType gammaZ;

    discrU.evaluateGradient ( el, edgePoint, gradU );
    _anisotropy.gammaFirstDerivative ( gradU, gammaZ );   // gamma_z( grad U )

    return gammaZ * normal;
  }
};

/**
 * This class computes the following boundary integral:
 * \f[ \int_{\partial\Omega} (\gamma_{zz}( \nabla u ) \nabla \varphi) \cdot n \, w da \f]
 * \f[ u \f] has to be provided by setImageReference, \f[ w \f] is the argument of apply.
 */
template <typename ConfiguratorType, typename BoundaryQuadratureType, typename AnisoType>
class IntegrateGammaZZWOverBoundary
      : public BoundaryWeightedMatrixDiffIntegrationInterface < ConfiguratorType,
      aol::Vector<typename ConfiguratorType::RealType>,
      ConfiguratorType::Dim,
      BoundaryQuadratureType,
      IntegrateGammaZZWOverBoundary< ConfiguratorType, BoundaryQuadratureType, AnisoType > > {

  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrImg;
  const AnisoType &_anisotropy;

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

public:
  IntegrateGammaZZWOverBoundary ( const typename ConfiguratorType::InitType &Initializer, const AnisoType &Anisotropy  )
      : BoundaryWeightedMatrixDiffIntegrationInterface < ConfiguratorType,
      aol::Vector<RealType>,
      ConfiguratorType::Dim,
      BoundaryQuadratureType,
      IntegrateGammaZZWOverBoundary< ConfiguratorType, BoundaryQuadratureType, AnisoType > > ( Initializer ),
      _discrImg ( NULL ), _anisotropy ( Anisotropy ) {  }

  //! set image (by reference!)
  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( _discrImg )
      delete _discrImg;
    _discrImg = new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->getConfigurator().getInitializer(), Image );
  }

  //! evaluate W
  RealType getCoeff ( const aol::Vector<RealType> &W, const VecType &/*normal*/, const ElementType &el, const VecType edgePoint ) const {
    aol::DiscreteFunctionDefault<ConfiguratorType> discrW ( this->_grid, W );

    return discrW.evaluate ( el, edgePoint );
  }

  //! evaluate \f[ gamma_{zz}( \nabla U ) \f]
  void getCoeffMatrix ( const aol::Vector<RealType> &/*W*/, const VecType &/*normal*/, const ElementType &el,
                        const VecType edgePoint, MatType &Matrix ) const {
    if ( !_discrImg ) throw aol::Exception ( "qc::IntegrateGammaZZWOverBoundary: Set Image first!\n", __FILE__, __LINE__ );
    VecType grad;
    _discrImg->evaluateGradient ( el, edgePoint, grad );

    _anisotropy.implicitPart ( grad, Matrix );     // gamma_{zz} ( grad u )
  }

};



// ---------------------------------------  ENERGIES --------------------------------------------

//! Computes \f$ E[\Phi] = \int_{\Omega} \left( div \gamma_z (\nabla \Phi) \right)^2 \|\nabla \Phi\|^2 \, dx \f$,
//! which is just the energy of the anisotropic willmore flow.
template <typename ConfiguratorType, typename BoundaryQuadratureType, typename AnisoType>
class AnisoWillmoreEnergyOp : public aol::Op < aol::Vector<typename ConfiguratorType::RealType>,
      aol::Scalar<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_initializer;
  const AnisoType &_anisotropy;
  AnisotropyIntegrationOp< ConfiguratorType, AnisoType > _anisoIntegrationOp;
  IntegrateGammaZOverBoundary< ConfiguratorType, BoundaryQuadratureType, AnisoType > _boundaryGammaZOp;
  mutable MCMLumpedMassOp<ConfiguratorType> _invMassOp;

public:
  AnisoWillmoreEnergyOp (  const typename ConfiguratorType::InitType &Initializer,
                           const AnisoType &Anisotropy, const RealType Epsilon = 1. )
      : _initializer ( Initializer ), _anisotropy ( Anisotropy ),
      _anisoIntegrationOp ( _initializer, _anisotropy ),
      _boundaryGammaZOp ( _initializer, _anisotropy ),
      _invMassOp ( _initializer, Epsilon, aol::INVERT ) {}

  ~AnisoWillmoreEnergyOp( ) {  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {

    _invMassOp.setImageReference ( Arg );

    aol::Vector<RealType> W ( Arg.size() );
    aol::Vector<RealType> temp ( Arg.size() );

    _boundaryGammaZOp.apply ( Arg, temp );
    temp *= -1.;
    _anisoIntegrationOp.applyAdd ( Arg, temp );       // temp = A_gamma[Phi] - B Phi
    _invMassOp.apply ( temp, W );                     // W = M^-1 (A_gamma[Phi] - B Phi)

    Dest[0] = temp * W;                               // (A_gamma[Phi] - B Phi) * M^-1 (A_gamma[Phi] - B Phi)
    Dest[0] *= 0.5;
  }

  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Scalar<RealType> &/*Dest*/ ) const {
    throw aol::Exception ( "Not implemented", __FILE__, __LINE__ );
  }

};

}    // end namespace qc


#endif
