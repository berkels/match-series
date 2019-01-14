#ifndef __ANISOSTIFFOPS_H
#define __ANISOSTIFFOPS_H

#include <quoc.h>
#include <aol.h>
#include <op.h>
#include <gridBase.h>
#include <linearSmoothOp.h>
#include <FEOpInterface.h>
#include <cellCenteredGrid.h>

namespace qc {

/**
 * Class to evaluate \f$g(s)=\frac{1}{\alpha+\beta{s^2}}\f$.
 *
 * \author Berkels
 */
template <typename RealType>
class PeronaMalikWeightingFunction {
  const RealType _alpha;
  const RealType _beta;
public:
  PeronaMalikWeightingFunction ( const RealType Alpha, const RealType Beta )
    : _alpha ( Alpha ),
      _beta ( Beta )
  {
  }
  RealType evaluate( const RealType x ) const {
    return aol::NumberTrait<RealType>::one/(_alpha + _beta*aol::Sqr( x ) );
  }
};

/** \brief class for anisotropic diffusion filtering in 3D according to
 *         Perona Malik (NOT Level Set case)
 *
 * \f$  \partial_t \Phi - \mathrm{div} \left( g(\|\nabla \Phi\|) \nabla \Phi \right) =0,
      \quad g(s)=\frac{1}{1+\left(\frac{s}{\lambda}\right)^2} \f$
 *
 * \author Droske
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType>
class PeronaMalikStiffOp : public aol::FELinScalarWeightedStiffInterface<ConfiguratorType, PeronaMalikStiffOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VecType  VecType;

  PeronaMalikStiffOp ( const typename ConfiguratorType::InitType &Grid, aol::OperatorType OpType = aol::ONTHEFLY, RealType Lambda = 1. )
      : aol::FELinScalarWeightedStiffInterface<ConfiguratorType, PeronaMalikStiffOp<ConfiguratorType> > ( Grid, OpType ), _pDiscrImg ( ), _lambda ( Lambda )  {}

  void setImageReference ( const aol::Vector<RealType> &Image ) {
    _pDiscrImg.reset ( new aol::DiscreteFunctionDefault<ConfiguratorType> ( this->_grid, Image ), true );
    this->reset();
  }

  void setLambda ( RealType Lambda ) {
    _lambda = Lambda;
  }

  inline RealType getCoeff ( const qc::Element &El, int QuadPoint, const VecType& /*RefCoord*/ ) const {
    if ( _pDiscrImg.get() == NULL ) {
      throw aol::Exception ( "PeronaMalikStiffOp::getCoeff: set image first!", __FILE__, __LINE__ );
    }

    VecType grad;
    _pDiscrImg->evaluateGradientAtQuadPoint ( El, QuadPoint, grad );
    return ( 1. / ( 1. + grad.normSqr() / aol::Sqr ( _lambda ) ) );
  }

protected:
  aol::DeleteFlagPointer<const aol::DiscreteFunctionDefault<ConfiguratorType> > _pDiscrImg;
  RealType _lambda;
};

/** \brief class for anisotropic diffusion filtering in 3D according to Perona Malik (Level Set case)
 *
 * \f$  \partial_t \Phi - \|\nabla \Phi\|
      \mathrm{div} \left( g(\|\nabla \Phi\|) \frac{ \nabla \Phi}{\| \nabla \Phi\|} \right) =0,
      \quad g(s)=\frac{1}{1+\left(\frac{s}{\lambda}\right)^2} \f$
 *
 * \author Nemitz
 * (last date: 17.03.2006)
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType>
class PeronaMalikLevelSetStiffOp : public aol::FELinScalarWeightedStiffInterface<ConfiguratorType, PeronaMalikLevelSetStiffOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VecType  VecType;

  PeronaMalikLevelSetStiffOp ( const typename ConfiguratorType::InitType &Grid, aol::OperatorType OpType = aol::ONTHEFLY,
                               RealType Lambda = 1., RealType Epsilon = 1. )
      : aol::FELinScalarWeightedStiffInterface<ConfiguratorType, PeronaMalikLevelSetStiffOp<ConfiguratorType> > ( Grid, OpType ),
      _img ( NULL ), _lambda ( Lambda ), _eps ( Epsilon ) {}

  void setImageReference ( const aol::Vector<RealType> &Image ) {
    _img = &Image;
    this->reset();
  }

  void setLambda ( RealType Lambda ) {
    _lambda = Lambda;
  }

  inline RealType getCoeff ( const qc::Element &El, int QuadPoint, const VecType& /*RefCoord*/ ) const {
    if ( !_img ) {
      throw aol::Exception ( "PeronaMalikLevelSetStiffOp::getCoeff: set image first!", __FILE__, __LINE__ );
    }

    VecType grad, v;
    grad.setZero();

    for ( int b = 0; b < this->getConfigurator().getNumLocalDofs ( El ); b++ ) {
      v = this->getBaseFunctionSet ( El ).evaluateGradient ( b, QuadPoint );
      v *= ( *_img ) [ this->getConfigurator().localToGlobal ( El, b ) ];
      grad += v;
    }
    grad /= this->_grid.H();

    RealType norm = sqrt ( grad.normSqr() + aol::Sqr ( _eps ) );

    cerr << "PeronaMalikLevelSetStiffOp seems to be broken" << endl; // should the term 1 + s^2 / lambda^2 be in the numerator???
    return ( 1. / norm * ( 1. + grad.normSqr() / aol::Sqr ( _lambda ) ) );
  }

protected:
  const aol::Vector<RealType> *_img;
  RealType _lambda;
  RealType _eps;
};



/**
 * \brief class for anisotropic diffusion filtering for 2d and 3d and arbitrary base functions
 * \author Droske
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType>
class TVStiffOp : public aol::FELinScalarWeightedStiffInterface<ConfiguratorType, TVStiffOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;

  TVStiffOp ( const qc::GridDefinition &Grid, aol::OperatorType OpType = aol::ONTHEFLY, RealType Epsilon = 1. )
      : aol::FELinScalarWeightedStiffInterface<ConfiguratorType, TVStiffOp<ConfiguratorType> > ( Grid, OpType ), _img ( NULL ), _epsilon ( Epsilon )  {}

  void setImageReference ( const aol::MultiVector<RealType> &Image ) {
    _img = &Image;
    this->reset();
  }

  void setEpsilon ( RealType Epsilon ) {
    _epsilon = Epsilon;
  }

  inline RealType getCoeff ( const qc::Element &El, int QuadPoint, const aol::Vec2<RealType>& /*RefCoord*/ ) const {
    if ( !_img ) {
      throw aol::Exception ( "TVStiffOp::getCoeff: set image first!", __FILE__, __LINE__ );
    }

    aol::Vec2<RealType> grad, v;
    grad.setZero();

    RealType gSqr = 0.;

    for ( int c = 0; c < _img->numComponents(); c++ ) {
      grad.setZero();
      for ( int b = 0; b < this->getConfigurator().getNumLocalDofs ( El ); b++ ) {
        v = this->getBaseFunctionSet ( El ).evaluateGradient ( b, QuadPoint );
        v *= ( *_img ) [ c ][ this->getConfigurator().localToGlobal ( El, b ) ];
        grad += v;
      }
      gSqr += grad.normSqr();
    }
    gSqr /= aol::Sqr ( this->_grid.H() );

    return 1. / sqrt ( gSqr + aol::Sqr ( _epsilon ) );
  }
protected:
  const aol::MultiVector<RealType> *_img;
  RealType _epsilon;
};



/**
 * \brief class for anisotropic diffusion filtering in 2D
 * \author Droske
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType>
class AnisotropicDiffusion2DStiffOp : public aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, AnisotropicDiffusion2DStiffOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  RealType _lambda, _struct_sigma, _struct_rho;
  qc::ScalarArray<RealType, qc::QC_2D> _struct1, _struct2, _struct3, _filter_im;
  qc::LinearSmoothOp<RealType, typename qc::MultilevelArrayTrait<RealType, typename ConfiguratorType::InitType>::GridTraitType> linSmooth;
public:
  AnisotropicDiffusion2DStiffOp ( const typename ConfiguratorType::InitType &Grid, aol::OperatorType OpType = aol::ONTHEFLY, const RealType Lambda = 1. )
      : aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, AnisotropicDiffusion2DStiffOp<ConfiguratorType> > ( Grid, OpType ),
      _lambda ( Lambda ),
      _struct_sigma ( 0.5 * Grid.H() ), _struct_rho ( 4.0 * Grid.H() ),
      _struct1 ( Grid ), _struct2 ( Grid ), _struct3 ( Grid ), _filter_im ( Grid ) {
    linSmooth.setCurrentGrid ( Grid );
  }


  // controls sensitivity of diffusion in estimated normal direction with respect to gradient magnitude
  void setLambda ( RealType Lambda ) {
    _lambda = Lambda;
  }

  /*!
   * Set presmoothing parameters
   * @param[in] Sigma Presmooting width of the image before structure tensor is computed
   * @param[in] Rho   Presmooting width of the structure tensor components
   * Parameters are not multiplied by the discretization parameter
   */
  void setStructSmooth ( RealType Sigma, RealType Rho ) {
    _struct_sigma = Sigma;
    _struct_rho = Rho;
  }

  inline void getCoeffMatrix ( const qc::Element &El, int /*QuadPoint*/,
                               const typename ConfiguratorType::VecType& RefCoord,
                               typename ConfiguratorType::MatType &Matrix ) const {
    RealType s1 = _struct1.interpolate ( El.x() + RefCoord[0], El.y() + RefCoord[1] );
    RealType s2 = _struct2.interpolate ( El.x() + RefCoord[0], El.y() + RefCoord[1] );
    RealType s3 = _struct3.interpolate ( El.x() + RefCoord[0], El.y() + RefCoord[1] );

    // See section 2.2 of Joachim Weickert's book Anisotropic Diffusion in Image Processing for details on the structure tensor.
    RealType v = sqrt ( aol::Sqr ( s1 - s3 ) + 4. * s2 * s2 );

    // Eigenvalues of the structure tensor.
    RealType mu1 = 0.5 * ( s1 + s3 + v );
    RealType mu2 = 0.5 * ( s1 + s3 - v );

    // Eigenvectors of the structure tensor.
    aol::Vec2<RealType> ev1, ev2;
    ev1 [ 0 ] = 1.;
    if ( aol::Abs ( s2 ) > aol::Abs ( s3 - mu1 ) ) {
      ev1 [ 1 ] =  ( mu1 - s1 ) / s2;
    } else {
      ev1 [ 1 ] =   s2 / ( mu1 - s3 );
    }

    ev2 [ 0 ] = 1.;
    if ( aol::Abs ( s2 ) > aol::Abs ( s3 - mu2 ) ) {
      ev2 [ 1 ] =  ( mu2 - s1 ) / s2;
    } else {
      ev2 [ 1 ] =   s2 / ( mu2 - s3 );
    }

    ev1.normalize();
    ev2.normalize();

    // tangential auf der kante
    RealType lambda2 = ( aol::Abs ( mu1 - mu2 ) < 1e-16 ) ? _lambda : ( _lambda + ( 1. - _lambda ) * exp ( -1. / aol::Sqr ( mu1 - mu2 ) ) );
    RealType lambda1 = 1. / ( 1. + mu1 ); // orthogonal zur kante

    // Matrix = lamda1 ev1 \otimes ev1 + lambda2 ev2 \otimes ev2
    Matrix.set ( 0, 0, lambda1 * aol::Sqr ( ev1 [ 0 ] ) + lambda2 * aol::Sqr ( ev2 [ 0 ] ) );
    Matrix.set ( 1, 0, lambda1 * ev1 [ 1 ] * ev1 [ 0 ] + lambda2 * ev2 [ 0 ] * ev2 [ 1 ] );
    Matrix.set ( 0, 1, lambda1 * ev1 [ 1 ] * ev1 [ 0 ] + lambda2 * ev2 [ 0 ] * ev2 [ 1 ] );
    Matrix.set ( 1, 1, lambda1 * aol::Sqr ( ev1 [ 1 ] ) + lambda2 * aol::Sqr ( ev2 [ 1 ] ) );
  }


  void computeStructureTensor ( const qc::ScalarArray<RealType, qc::QC_2D> &Img ) {
    cerr << "Computing structure tensor...";

    this->reset();

    linSmooth.setSigma ( _struct_sigma );
    linSmooth.apply ( Img, _filter_im );

    typename ConfiguratorType::InitType::OldAllNodeIterator it;
    for ( it = this->_grid._nBeginIt; it != this->_grid._nEndIt; ++it ) {
      double dx = _filter_im.dxFD ( it->x(), it->y() ) / this->_grid.H();
      double dy = _filter_im.dyFD ( it->x(), it->y() ) / this->_grid.H();
      _struct1.set ( *it, dx*dx );
      _struct2.set ( *it, dx*dy );
      _struct3.set ( *it, dy*dy );
    }
    linSmooth.setSigma ( _struct_rho );
    linSmooth.apply ( _struct1, _struct1 );
    linSmooth.apply ( _struct2, _struct2 );
    linSmooth.apply ( _struct3, _struct3 );
    cerr << "done.\n";
  }
};


/** \brief class for anisotropic diffusion filtering in 3D (NOT Level Set case)
 *
 * \f$  \partial_t \Phi -
      \mathrm{div} \left( A(\nabla \Phi) \nabla \Phi \right) =0,
      \quad A = Q \left( \begin{array}{ccc} g(\lambda_1) &&  \\ & \ddots & \\ && g(\lambda_n)
      \end{array} \right) Q^T \f$
 *
 * \author Nemitz (with regard to Marc's template)
 * \note Before apply, you always have to compute the structure tensor!
 * \ingroup MatrixFEOp
 */
template < typename ConfiguratorType, typename LinearSmootherType = qc::LinearSmoothOp< typename ConfiguratorType::RealType > >
class AnisotropicDiffusion3DStiffOp : public aol::FELinMatrixWeightedStiffInterface< ConfiguratorType, AnisotropicDiffusion3DStiffOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;

protected:
  RealType _lambda, _sigma, _rho;
  qc::ScalarArray<RealType, qc::QC_3D> _struct00, _struct10, _struct11, _struct20, _struct21, _struct22, _imgSmooth;
  LinearSmootherType _linSmooth;

public:
  AnisotropicDiffusion3DStiffOp ( const typename ConfiguratorType::InitType &Grid, aol::OperatorType OpType = aol::ONTHEFLY )
      : aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, AnisotropicDiffusion3DStiffOp<ConfiguratorType> > ( Grid, OpType ),
      _lambda ( 1.0 ),
      _sigma ( 0.5 * Grid.H() ), _rho ( 4.0 * Grid.H() ),
      _struct00 ( Grid ), _struct10 ( Grid ), _struct11 ( Grid ), _struct20 ( Grid ), _struct21 ( Grid ), _struct22 ( Grid ), _imgSmooth ( Grid ),
      _linSmooth ( Grid ) {
  }

  void setLambda ( RealType Lambda ) {
    _lambda = Lambda;
  }

  void setStructSmooth ( RealType Sigma, RealType Rho ) {
    _sigma = Sigma;
    _rho = Rho;
  }

  inline RealType g ( const RealType s ) const {
    return ( 1. / ( 1. + ( ( s*s ) / ( _lambda*_lambda ) ) ) );
  }

  //! this method computes the anisotropic weighted coefficient matrix
  inline void getCoeffMatrix ( const qc::Element &El, int /*QuadPoint*/,
                               const typename ConfiguratorType::VecType& RefCoord,
                               typename ConfiguratorType::MatType &Matrix ) const {
    const RealType s00 = _struct00.interpolate ( El.x() + RefCoord[0], El.y() + RefCoord[1], El.z() + RefCoord[2] );
    const RealType s10 = _struct10.interpolate ( El.x() + RefCoord[0], El.y() + RefCoord[1], El.z() + RefCoord[2] );
    const RealType s11 = _struct11.interpolate ( El.x() + RefCoord[0], El.y() + RefCoord[1], El.z() + RefCoord[2] );
    const RealType s20 = _struct20.interpolate ( El.x() + RefCoord[0], El.y() + RefCoord[1], El.z() + RefCoord[2] );
    const RealType s21 = _struct21.interpolate ( El.x() + RefCoord[0], El.y() + RefCoord[1], El.z() + RefCoord[2] );
    const RealType s22 = _struct22.interpolate ( El.x() + RefCoord[0], El.y() + RefCoord[1], El.z() + RefCoord[2] );

    // first build the symmetric 3x3-matrix of the tensor product
    aol::Matrix33Symm<RealType> Tensor ( s00, s10, s20,
                                         s10, s11, s21,
                                         s20, s21, s22 );

    // now compute the eigenvalues and -vectors of this matrix
    aol::Matrix33<RealType> eigenVectors;
    aol::Vec3<RealType> eigenValues;

    Tensor.eigenVectors ( eigenValues, eigenVectors );

    // apply g to the eigenValues
    for ( int i = 0; i < 3; i++ )
      eigenValues[i] = g ( eigenValues[i] );

    Matrix.setZero();

    // Finally build the coefficient matrix
    for ( int i = 0; i < 3; i++ )   // tensor rows
      for ( int j = 0; j < 3; j++ )   // tensor columns
        for ( int k = 0; k < 3; k++ )   // the sum in each entry
          Matrix.add ( i, j, eigenValues[k] * eigenVectors[i][k] * eigenVectors[j][k] );
  }

  //! this method computes the structure tensor and smoothes it component-wise
  void computeStructureTensor ( const qc::ScalarArray<RealType, qc::QC_3D> &Img ) {
    cerr << "Computing structure tensor...";

    this->reset();

    // first smooth the image (generate u_sigma)
    _linSmooth.setSigma ( _sigma );
    _linSmooth.apply ( Img, _imgSmooth );

    // compute the tensor product of the smoothed gradient
    for ( qc::RectangularIterator<qc::QC_3D> it ( this->_grid ); it.notAtEnd(); ++it ) { // FullNodeIterator cannot be used here
      const RealType dx = _imgSmooth.dxFD ( it->x(), it->y(), it->z() ) / this->_grid.H();
      const RealType dy = _imgSmooth.dyFD ( it->x(), it->y(), it->z() ) / this->_grid.H();
      const RealType dz = _imgSmooth.dzFD ( it->x(), it->y(), it->z() ) / this->_grid.H();
      _struct00.set ( *it, dx*dx );
      _struct10.set ( *it, dx*dy );
      _struct11.set ( *it, dy*dy );
      _struct20.set ( *it, dx*dz );
      _struct21.set ( *it, dy*dz );
      _struct22.set ( *it, dz*dz );
    }

    // finally smooth the components of the tensor product
    _linSmooth.setSigma ( _rho );
    _linSmooth.apply ( _struct00, _struct00 );
    _linSmooth.apply ( _struct10, _struct10 );
    _linSmooth.apply ( _struct11, _struct11 );
    _linSmooth.apply ( _struct20, _struct20 );
    _linSmooth.apply ( _struct21, _struct21 );
    _linSmooth.apply ( _struct22, _struct22 );
    cerr << "done.\n";
  }

};


/**  class for anisotropic diffusion filtering in 3D (Level Set case)
 *
 * \f$  \partial_t \Phi - \|\nabla \Phi\|
      \mathrm{div} \left( A(\nabla \Phi) \frac{ \nabla \Phi}{\| \nabla \Phi\|} \right) =0,
      \quad A = Q \left( \begin{array}{ccc} g(\lambda_1) &&  \\ & \ddots & \\ && g(\lambda_n)
      \end{array} \right) Q^T \f$
 *
 * \author Nemitz (with regard to Marc's template)
 * Notice: Before apply, you always have to compute the structure tensor!
 */
template < typename ConfiguratorType, typename LinearSmootherType = qc::LinearSmoothOp< typename ConfiguratorType::RealType > >
class AnisotropicDiffusion3DLevelSetStiffOp : public AnisotropicDiffusion3DStiffOp< ConfiguratorType, LinearSmootherType > {
public:
  typedef typename ConfiguratorType::RealType RealType;

protected:
  const aol::Vector<RealType> *_img;
  RealType _eps;

public:
  AnisotropicDiffusion3DLevelSetStiffOp ( const typename ConfiguratorType::InitType &Grid, aol::OperatorType OpType = aol::ONTHEFLY ) :
    AnisotropicDiffusion3DStiffOp<ConfiguratorType, LinearSmootherType> ( Grid, OpType ) {
  }

  // -------- Zeiger auf den letzten Zeitschritt setzen (semi-implizites Verfahren) --------
  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( static_cast<int> ( Image.size() ) != this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "Array.size() = " << Image.size() << " getConfigurator().getNumGlobalDofs() = " << this->getConfigurator().getNumGlobalDofs()  << endl;
      throw aol::Exception ( "Array phi has wrong size.", __FILE__, __LINE__ );
    }
    _img = &Image;
    this->reset();
  }

  //! this method computes the anisotropic weighted coefficient matrix
  inline void getCoeffMatrix ( const qc::Element &El, int QuadPoint,
                               const typename ConfiguratorType::VecType& RefCoord,
                               typename ConfiguratorType::MatType &Matrix ) const {
    AnisotropicDiffusion3DStiffOp<ConfiguratorType,LinearSmootherType>::getCoeffMatrix ( El, QuadPoint, RefCoord, Matrix );

    typename ConfiguratorType::VecType grad;
    gradient ( El, QuadPoint, grad );

    RealType norm = sqrt ( grad.normSqr() + aol::Sqr ( _eps ) );

    Matrix /= norm;
  }

  // -------- Gradienten des letzten Zeitschritts berechnen -------------------------
  inline void gradient ( const qc::Element &El, int QuadPoint, typename ConfiguratorType::VecType &Grad ) const {
    if ( !_img ) {
      throw aol::Exception ( "AnisotropicDiffusion3DLevelSetStiffOp::getCoeff: set image first!", __FILE__, __LINE__ );
    }

    typename ConfiguratorType::VecType v;
    Grad.setZero();
    for ( int b = 0; b < this->getConfigurator().getNumLocalDofs ( El ); b++ ) {
      v = this->getBaseFunctionSet ( El ).evaluateGradient ( b, QuadPoint );
      v *= ( *_img ) [ this->getConfigurator().localToGlobal ( El, b ) ];
      Grad += v;
    }
    Grad /= this->getConfigurator().H ( El );
  }

};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename PeronaMalikFunctionType>
class ImageDrivenAnisoTVVecNorm2D : public aol::FENonlinIntegrationVectorInterface < ConfiguratorType,  ImageDrivenAnisoTVVecNorm2D<ConfiguratorType, PeronaMalikFunctionType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrImage;
  const PeronaMalikFunctionType &_peronaMalikFunction;
  const RealType _lambda;
  const RealType _epsilonSqr;
protected:
public:
  ImageDrivenAnisoTVVecNorm2D ( const typename ConfiguratorType::InitType &Initializer,
                                const aol::Vector<RealType> &ImageDofs,
                                const PeronaMalikFunctionType &PeronaMalikFunction,
                                const RealType Lambda,
                                const RealType Epsilon )
    : aol::FENonlinIntegrationVectorInterface < ConfiguratorType,
                                                  ImageDrivenAnisoTVVecNorm2D<ConfiguratorType, PeronaMalikFunctionType> > ( Initializer ),
      _discrImage ( Initializer, ImageDofs ),
      _peronaMalikFunction ( PeronaMalikFunction ),
      _lambda ( Lambda ),
      _epsilonSqr ( aol::Sqr( Epsilon ) )
  {}

  RealType evaluateIntegrand ( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &discrFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/ ) const {

    typename ConfiguratorType::VecType gradU;
    _discrImage.evaluateGradientAtQuadPoint( El, QuadPoint, gradU );

    const RealType normGradU = sqrt ( gradU.normSqr() + _epsilonSqr );
    typename ConfiguratorType::VecType normal = gradU;
    normal /= normGradU;

    typename ConfiguratorType::VecType normalOrthogonal;
    normalOrthogonal[0] = normal[1];
    normalOrthogonal[1] = -normal[0];

    const RealType g = _peronaMalikFunction.evaluate( normGradU );

    typename ConfiguratorType::VecType gradVX;
    typename ConfiguratorType::VecType gradVY;

    discrFuncs[0].evaluateGradientAtQuadPoint( El, QuadPoint, gradVX );
    discrFuncs[1].evaluateGradientAtQuadPoint( El, QuadPoint, gradVY );

    return ( g * sqrt ( gradVX.normSqr() + gradVY.normSqr() + _epsilonSqr )
           + ( 1. - g ) * sqrt ( aol::Sqr(_lambda*normal*gradVX) + aol::Sqr(normalOrthogonal*gradVX)
                                 + aol::Sqr(_lambda*normal*gradVY) + aol::Sqr(normalOrthogonal*gradVY) + _epsilonSqr ) );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename PeronaMalikFunctionType>
class VariationOfImageDrivenAnisoTVVecNorm2D : public aol::FENonlinVectorDiffOpInterface< ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, VariationOfImageDrivenAnisoTVVecNorm2D<ConfiguratorType, PeronaMalikFunctionType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrImage;
  const PeronaMalikFunctionType &_peronaMalikFunction;
  const RealType _lambdaSqr;
  const RealType _epsilonSqr;
public:
  VariationOfImageDrivenAnisoTVVecNorm2D ( const typename ConfiguratorType::InitType &Initializer,
                                           const aol::Vector<RealType> &ImageDofs,
                                           const PeronaMalikFunctionType &PeronaMalikFunction,
                                           const RealType Lambda,
                                           const RealType Epsilon )
    :aol::FENonlinVectorDiffOpInterface< ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, VariationOfImageDrivenAnisoTVVecNorm2D<ConfiguratorType, PeronaMalikFunctionType> >( Initializer ),
    _discrImage( Initializer, ImageDofs ),
    _peronaMalikFunction( PeronaMalikFunction ),
    _lambdaSqr( aol::Sqr(Lambda) ),
    _epsilonSqr( aol::Sqr(Epsilon) ) {
  }


  void getNonlinearity ( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
                         aol::Mat<ConfiguratorType::Dim, ConfiguratorType::Dim, typename ConfiguratorType::RealType> &NL ) const {

    typename ConfiguratorType::VecType gradU;
    _discrImage.evaluateGradientAtQuadPoint( El, QuadPoint, gradU );

    const RealType normGradU = sqrt ( gradU.normSqr() + _epsilonSqr );
    typename ConfiguratorType::VecType normal = gradU;
    normal /= normGradU;

    typename ConfiguratorType::VecType normalOrthogonal;
    normalOrthogonal[0] = normal[1];
    normalOrthogonal[1] = -normal[0];

    const RealType g = _peronaMalikFunction.evaluate( normGradU );

    typename ConfiguratorType::VecType gradV0;
    DiscFuncs[0].evaluateGradientAtQuadPoint( El, QuadPoint, gradV0 );
    typename ConfiguratorType::VecType gradV1;
    DiscFuncs[1].evaluateGradientAtQuadPoint( El, QuadPoint, gradV1 );

    NL.setRow( 0, gradV0 );
    NL.setRow( 1, gradV1 );
    NL *= g / sqrt( gradV0.normSqr() + gradV1.normSqr() + _epsilonSqr );

    const RealType normalTimesV0 = normal * gradV0;
    const RealType normalOrthogonalTimesV0 = normalOrthogonal * gradV0;
    const RealType normalTimesV1 = normal * gradV1;
    const RealType normalOrthogonalTimesV1 = normalOrthogonal * gradV1;

    RealType gamma = sqrt( _lambdaSqr*aol::Sqr(normalTimesV0) + aol::Sqr(normalOrthogonalTimesV0)
                           + _lambdaSqr*aol::Sqr(normalTimesV1) + aol::Sqr(normalOrthogonalTimesV1) + _epsilonSqr );

    typename ConfiguratorType::MatType tempMat;

    typename ConfiguratorType::VecType temp;
    temp.addMultiple( normal, _lambdaSqr*normalTimesV0 );
    temp.addMultiple( normalOrthogonal, normalOrthogonalTimesV0 );

    tempMat.setRow( 0, temp );

    temp.setZero();
    temp.addMultiple( normal, _lambdaSqr*normalTimesV1 );
    temp.addMultiple( normalOrthogonal, normalOrthogonalTimesV1 );
    tempMat.setRow( 1, temp );
    tempMat *= (1.-g)/gamma;

    NL += tempMat;
  }
};


}

#endif
