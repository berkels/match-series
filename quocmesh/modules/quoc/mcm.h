#ifndef __MCM_H
#define __MCM_H

#include <FEOpInterface.h>
#include <baseFunctionSet.h>
#include <scalarArray.h>
#include <fastUniformGridMatrix.h>
#include <solver.h>
#include <preconditioner.h>
#include <qmException.h>
#include <maskedVector.h>
#include <discreteFunction.h>
#include <quocMatrices.h>

namespace qc {

/**
 * \brief calculates \f$ \int_\Omega \frac{1}{|\nabla u|_\epsilon}\nabla\varphi_i\cdot\nabla\varphi_j dx \f$,
 *        where \f$ u \f$ is the image specified by setImageReference.
 *
 * \author Droske
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class MCMStiffOp : public aol::FELinScalarWeightedStiffInterface<ConfiguratorType, MCMStiffOp<ConfiguratorType, IndexMode>, IndexMode > {

public:
  typedef typename ConfiguratorType::RealType RealType;

  MCMStiffOp ( const typename ConfiguratorType::InitType &Initializer, aol::OperatorType OpType = aol::ONTHEFLY,
               RealType Epsilon = 1., RealType Scale = 1.  )
      : aol::FELinScalarWeightedStiffInterface<ConfiguratorType, MCMStiffOp<ConfiguratorType, IndexMode>, IndexMode > ( Initializer, OpType ),
        _discrImg ( NULL ), _epsSqr ( Epsilon*Epsilon ), _scale( Scale )  {}

  ~MCMStiffOp( ) {
    if ( _discrImg ) delete _discrImg;
  }

  //! set image (by reference!)
  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( _discrImg )
      delete _discrImg;
    _discrImg = new aol::DiscreteFunctionDefault<ConfiguratorType>( this->getConfigurator().getInitializer(), Image );
    this->reset();
  }

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType& ) const {
    if ( !_discrImg ) {
      throw aol::Exception ( "qcMCMStiffOp::getCoeff: set image first!", __FILE__, __LINE__ );
    }

    typename ConfiguratorType::VecType grad;
    _discrImg->evaluateGradientAtQuadPoint( El, QuadPoint, grad );

    return _scale / sqrt ( grad.normSqr() + _epsSqr );
  }

protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrImg;
  const RealType _epsSqr;
  const RealType _scale;

};



/**
 * \brief calculates \f$ \int_\Omega \frac{1}{|\nabla u|_\epsilon}\varphi_i\varphi_j dx \f$,
 *        where \f$ u \f$ is the image specified by setImageReference.
 *
 * \author Droske
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class MCMMassOp : public aol::FELinScalarWeightedMassInterface<ConfiguratorType, MCMMassOp<ConfiguratorType, IndexMode>, IndexMode > {
protected:
public:
  typedef typename ConfiguratorType::RealType RealType;

  MCMMassOp ( const typename ConfiguratorType::InitType &Initializer, aol::OperatorType OpType = aol::ONTHEFLY,
              RealType Epsilon = 1., RealType Scale = 1. )
      : aol::FELinScalarWeightedMassInterface<ConfiguratorType, MCMMassOp<ConfiguratorType, IndexMode>, IndexMode > ( Initializer, OpType ),
        _discrImg ( NULL ), _eps ( Epsilon ), _scale( Scale )  {}


  ~MCMMassOp( ) {
    if ( _discrImg ) delete _discrImg;
  }


  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( _discrImg )
      delete _discrImg;
    _discrImg = new aol::DiscreteFunctionDefault<ConfiguratorType>( this->getConfigurator().getInitializer(), Image );
    this->reset();
  }


  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::VecType& ) const {
    if ( !_discrImg ) {
      throw aol::Exception ( "qcMCMMassOp::getCoeff: set image first!", __FILE__, __LINE__ );
    }
    typename ConfiguratorType::VecType grad;
    _discrImg->evaluateGradientAtQuadPoint( El, QuadPoint, grad );
    return _scale / sqrt ( grad.normSqr() + _eps*_eps );
  }
protected:

  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrImg;
  RealType _eps;
  RealType _scale;

};


template <typename ConfiguratorType>
class MCMLumpedMassOp : public  aol::LumpedMassOpInterface<ConfiguratorType, MCMLumpedMassOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;

  MCMLumpedMassOp ( const typename ConfiguratorType::InitType &Initializer, RealType Epsilon, aol::LUMPED_MASS_OP_MODE Invert )
  : aol::LumpedMassOpInterface<ConfiguratorType, MCMLumpedMassOp<ConfiguratorType> > ( Initializer, Invert ),
        _discrImg ( NULL ), _epsSqr ( Epsilon*Epsilon ) {}

  ~MCMLumpedMassOp( ) {
    if ( _discrImg ) delete _discrImg;
  }

  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( _discrImg )
      delete _discrImg;
    _discrImg = new aol::DiscreteFunctionDefault<ConfiguratorType>( this->getConfigurator().getInitializer(), Image );
    this->reset();
  }

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::VecType& ) const {
    if ( !_discrImg ) {
      throw aol::Exception ( "qc::MCMLumpedMassOp::getCoeff: set image first!", __FILE__, __LINE__ );
    }

    typename ConfiguratorType::DomVecType grad;
    _discrImg->evaluateGradientAtQuadPoint( El, QuadPoint, grad );

    return 1. / sqrt ( grad.normSqr() + _epsSqr );
  }

protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrImg;
  RealType _epsSqr;
};




// ------------------------ AnisotropyStiffOp --------------------------------------

/**
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, typename AnisoType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class AnisotropyStiffOp : public aol::FELinScalarWeightedStiffInterface < ConfiguratorType,
      AnisotropyStiffOp<ConfiguratorType, AnisoType, IndexMode>, IndexMode >   {
public:
  typedef typename ConfiguratorType::RealType RealType;

  AnisotropyStiffOp ( const typename ConfiguratorType::InitType &Initializer,
                      aol::OperatorType OpType,
                      const AnisoType &Anisotropy, RealType epsilon )
      : aol::FELinScalarWeightedStiffInterface<ConfiguratorType, AnisotropyStiffOp<ConfiguratorType, AnisoType, IndexMode>, IndexMode > (
        Initializer, OpType ),  _anisotropy ( Anisotropy ),
      _discrImg ( NULL ), _epsSqr ( epsilon*epsilon ) {}


  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( static_cast<int> ( Image.size() ) != this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "array.size() = " << Image.size() << "  getConfigurator().getNumGlobalDofs() = ";
      cerr << this->getConfigurator().getNumGlobalDofs() << endl;
      throw aol::Exception ( "array phi has wrong size\n", __FILE__, __LINE__ );
    }
    if ( _discrImg )
      delete _discrImg;
    _discrImg = new aol::DiscreteFunctionDefault<ConfiguratorType>( this->getConfigurator().getInitializer(), Image );
    this->reset();
  }

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                             const typename ConfiguratorType::VecType& /*RefCoord*/ ) const {

    typename ConfiguratorType::VecType grad;
    _discrImg->evaluateGradientAtQuadPoint( El, QuadPoint, grad );
    return _anisotropy.gammaNorm ( grad ) / sqrt ( grad.normSqr() + _epsSqr );
  }

protected:
  const AnisoType &_anisotropy;
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrImg;
  RealType _epsSqr;

};


//! the difference to the AnisotropyStiffOp is just, that this op calls the anisotropy
//! with the Element and QuadPoint as additional arguments.
//! \ingroup MatrixFEOp
template <typename ConfiguratorType, typename AnisoType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class AnisotropyStiffOpPositionDepending : public aol::FELinScalarWeightedStiffInterface < ConfiguratorType,
      AnisotropyStiffOpPositionDepending<ConfiguratorType, AnisoType, IndexMode>, IndexMode >   {
public:
  typedef typename ConfiguratorType::RealType RealType;

  AnisotropyStiffOpPositionDepending ( const typename ConfiguratorType::InitType &Initializer,
                                       aol::OperatorType OpType,
                                       aol::Vector<RealType> &/*Dofs*/,
                                       const AnisoType &Anisotropy, RealType epsilon )
          : aol::FELinScalarWeightedStiffInterface<ConfiguratorType,
            AnisotropyStiffOpPositionDepending<ConfiguratorType, AnisoType, IndexMode>, IndexMode > (
        Initializer, OpType ), _img ( NULL ), _anisotropy ( Anisotropy ),
      _discrImg ( NULL ), _eps ( epsilon ) {}


  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( static_cast< int > ( Image.size() ) != this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "array.size() = " << Image.size() << "  getConfigurator().getNumGlobalDofs() = ";
      cerr << this->getConfigurator().getNumGlobalDofs() << endl;
      throw aol::Exception ( "array phi has wrong size\n", __FILE__, __LINE__ );
    }
    _img = &Image;
    if ( _discrImg )
      delete _discrImg;
    _discrImg = new aol::DiscreteFunctionDefault<ConfiguratorType>( this->getConfigurator().getInitializer(), Image );
    this->reset();
  }

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                             const typename ConfiguratorType::VecType& /*RefCoord*/ ) const {

    typename ConfiguratorType::VecType grad;
    _discrImg->evaluateGradientAtQuadPoint( El, QuadPoint, grad );
    return _anisotropy.gammaNorm ( El, QuadPoint, grad ) / sqrt ( grad.normSqr() + _eps*_eps );
  }

protected:
  const aol::Vector<RealType> *_img;
  const AnisoType &_anisotropy;
  const aol::DiscreteFunctionDefault<ConfiguratorType> *_discrImg;
  RealType _eps;

};




// ------------------------ LocalRotatedAnisotropyStiffOp ---------------------------------------
// anisotropy depends local on the element, a vector-field with the rotations-directions
// of the anisotropy must be given. Before evaluating the weight, the anisotropy
// is rotated into the direction given from the vector-field.
// ----------------------------------------------------------------------------------------------

//! \ingroup MatrixFEOp
template <typename ConfiguratorType, typename AnisoType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class LocalRotatedAnisotropyStiffOp : public aol::FELinScalarWeightedStiffInterface < ConfiguratorType,
      LocalRotatedAnisotropyStiffOp<ConfiguratorType, AnisoType, IndexMode>, IndexMode >   {

public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VecType VecType;

private:

  const qc::ScalarArray<RealType, qc::QC_3D>* _rotDirectx;     // the components of the vector-field
  const qc::ScalarArray<RealType, qc::QC_3D>* _rotDirecty;
  const qc::ScalarArray<RealType, qc::QC_3D>* _rotDirectz;

public:

  // constructor
  LocalRotatedAnisotropyStiffOp ( const typename ConfiguratorType::InitType &Initializer,
                                  aol::OperatorType OpType,
                                  aol::Vector<RealType> &Dofs,
                                  AnisoType &Anisotropy, RealType epsilon,
                                  const qc::ScalarArray<RealType, qc::QC_3D> *rotDx,
                                  const qc::ScalarArray<RealType, qc::QC_3D> *rotDy,
                                  const qc::ScalarArray<RealType, qc::QC_3D> *rotDz )
      : aol::FELinScalarWeightedStiffInterface < ConfiguratorType,
      LocalRotatedAnisotropyStiffOp<ConfiguratorType, AnisoType, IndexMode>, IndexMode > (
        Initializer, OpType ), _img ( NULL ), _anisotropy ( Anisotropy ),
      _discFunc ( Initializer, Dofs ), _eps ( epsilon ) {
    _rotDirectx = rotDx;
    _rotDirecty = rotDy;
    _rotDirectz = rotDz;
  }

  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( static_cast< int > ( Image.size() ) != this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "array.size() = " << Image.size() << "  getConfigurator().getNumGlobalDofs() = ";
      cerr << this->getConfigurator().getNumGlobalDofs() << endl;
      throw aol::Exception ( "array phi has wrong size\n", __FILE__, __LINE__ );
    }
    _img = &Image;
    this->reset();
  }

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, VecType RefCoord ) const {
    VecType grad;
    gradient ( El, QuadPoint, grad );
    return evaluateWeight ( El, QuadPoint, RefCoord ) /  sqrt ( grad.normSqr() + aol::Sqr ( _eps ) );
  }

  inline void gradient ( const qc::Element &El, int QuadPoint, VecType &Grad )
  const {
    if ( !_img ) {
      throw aol::Exception ( "LocalRotatedAnisotropyStiffOp::getCoeff: set image first!", __FILE__, __LINE__ );
    }

    VecType v;
    Grad.setZero();
    for ( int b = 0; b < this->getConfigurator().getNumLocalDofs ( El ); b++ ) {
      v = this->getBaseFunctionSet ( El ).evaluateGradient ( b, QuadPoint );
      v *= ( *_img ) [ this->getConfigurator().localToGlobal ( El, b ) ];
      Grad += v;
    }
  }

  // this method returns the anisotropie on the rhs, that means gamma_z(n)
  inline RealType evaluateWeight ( const qc::Element &El, int QuadPoint,
                                   const VecType RefCoord ) const {
    if ( !_img ) {
      throw aol::Exception ( "LocalRotatedAnisotropyStiffOp::getCoeff: set image first!", __FILE__, __LINE__ );
    }

    VecType coord;
    coord[0] = El.x() + RefCoord[0];
    coord[1] = El.y() + RefCoord[1];
    if ( ConfiguratorType::VecType::dim == 3 ) {
      coord[2] = El.z() + RefCoord[2];
    }

    typename ConfiguratorType::VecType Gradient;
    _discFunc.evaluateGradientAtQuadPoint ( El, QuadPoint, Gradient );

    aol::Vec3<RealType> rotDirect;
    rotDirect[0] = _rotDirectx->interpolate ( coord[0], coord[1], coord[2] );
    rotDirect[1] = _rotDirecty->interpolate ( coord[0], coord[1], coord[2] );
    rotDirect[2] = _rotDirectz->interpolate ( coord[0], coord[1], coord[2] );

    _anisotropy.setRotateDirection ( rotDirect );

    return _anisotropy.gamma ( Gradient );

  }

protected:
  const aol::Vector<RealType> *_img;
  AnisoType &_anisotropy;
  aol::DiscreteFunctionDefault<ConfiguratorType> _discFunc;
  RealType _eps;

};



// ---------------------------------------------------------------------------------


//! \ingroup MatrixFEOp
template <typename ConfiguratorType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class MCMWeightedStiffOp : public aol::FELinScalarWeightedStiffInterface<ConfiguratorType, MCMWeightedStiffOp<ConfiguratorType, IndexMode>, IndexMode > {
public:
  typedef typename ConfiguratorType::RealType RealType;

  MCMWeightedStiffOp ( const typename ConfiguratorType::InitType &Initializer, aol::OperatorType OpType, RealType Epsilon )
      : aol::FELinScalarWeightedStiffInterface<ConfiguratorType, MCMWeightedStiffOp<ConfiguratorType, IndexMode>, IndexMode >
        ( Initializer, OpType ), _img ( NULL ), _weight ( NULL ), _eps ( Epsilon ) {}

  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( static_cast< int > ( Image.size() ) != this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "array.size() = " << Image.size() << "  getConfigurator().getNumGlobalDofs() = " <<   this->getConfigurator().getNumGlobalDofs() << endl;
      throw aol::Exception ( "array phi has wrong size\n", __FILE__, __LINE__ );
    }
    _img = &Image;
    this->reset();
  }

  void setWeightFunction ( const aol::Vector<RealType> &Weight ) {
    if ( static_cast< int > ( Weight.size() ) !=  this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "array.size() = " << Weight.size() << "  getConfigurator().getNumGlobalDofs() = " <<   this->getConfigurator().getNumGlobalDofs() << endl;
      throw aol::Exception ( "array phi has wrong size\n", __FILE__, __LINE__ );
    }
    _weight = &Weight;
  }

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::VecType& /*RefCoord*/ ) const {
    typename ConfiguratorType::VecType grad;
    gradient ( El, QuadPoint, grad );
    return evaluateWeight ( El, QuadPoint ) /  sqrt ( grad.normSqr() + aol::Sqr ( _eps ) );
  }

  inline void gradient ( const typename ConfiguratorType::ElementType &El, int QuadPoint, typename ConfiguratorType::VecType &Grad ) const {
    if ( !_img || !_weight ) {
      throw aol::Exception ( "MCMWeightedStiffOp::getCoeff: set image first!", __FILE__, __LINE__ );
    }

    typename ConfiguratorType::VecType v;
    Grad.setZero();
    for ( int b = 0; b < this->getConfigurator().getNumLocalDofs ( El ); b++ ) {
      v = this->getBaseFunctionSet ( El ).evaluateGradient ( b, QuadPoint );
      v *= ( *_img ) [ this->getConfigurator().localToGlobal ( El, b ) ];
      Grad += v;
    }
  }

  inline RealType evaluateWeight ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    if ( !_img || !_weight ) {
      throw aol::Exception ( "MCMWeightedStiffOp::getCoeff: set image first!", __FILE__, __LINE__ );
    }

    RealType w = 0.;
    for ( int b = 0; b < this->getConfigurator().getNumLocalDofs ( El ); b++ ) {
      w += ( *_weight ) [ this->getConfigurator().localToGlobal ( El, b ) ] * this->getBaseFunctionSet ( El ).evaluate ( b, QuadPoint );
    }
    return w;
  }

protected:
  const aol::Vector<RealType> *_img;
  const aol::Vector<RealType> *_weight;
  RealType _eps;

};

/**
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class MCMWeightedMassOp : public aol::FELinScalarWeightedMassInterface<ConfiguratorType, MCMWeightedMassOp<ConfiguratorType, IndexMode>, IndexMode > {
public:
  typedef typename ConfiguratorType::RealType RealType;

  MCMWeightedMassOp ( const typename ConfiguratorType::InitType &Initializer, aol::OperatorType OpType, RealType Epsilon )
      : aol::FELinScalarWeightedMassInterface<ConfiguratorType, MCMWeightedMassOp<ConfiguratorType, IndexMode>, IndexMode >
        ( Initializer, OpType ), _img ( NULL ), _weight ( NULL ), _eps ( Epsilon ) {}

  void setImageReference ( const aol::Vector<RealType> &Image ) {
    if ( static_cast< int > ( Image.size() ) !=  this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "array.size() = " << Image.size() << "  getConfigurator().getNumGlobalDofs() = " <<   this->getConfigurator().getNumGlobalDofs() << endl;
      throw aol::Exception ( "array phi has wrong size\n", __FILE__, __LINE__ );
    }
    _img = &Image;
    this->reset();
  }

  void setWeightFunction ( const aol::Vector<RealType> &Weight ) {
    if ( static_cast< int > ( Weight.size() ) !=  this->getConfigurator().getNumGlobalDofs() ) {
      cerr << "array.size() = " << Weight.size() << "  getConfigurator().getNumGlobalDofs() = " <<   this->getConfigurator().getNumGlobalDofs() << endl;
      throw aol::Exception ( "array phi has wrong size\n", __FILE__, __LINE__ );
    }
    _weight = &Weight;
  }

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::VecType& /*RefCoord*/ ) const {
    typename ConfiguratorType::VecType grad;
    gradient ( El, QuadPoint, grad );
    return evaluateWeight ( El, QuadPoint ) / sqrt ( grad.normSqr() + aol::Sqr ( _eps ) );
  }

  inline void gradient ( const typename ConfiguratorType::ElementType &El, int QuadPoint, typename ConfiguratorType::VecType &Grad ) const {
    if ( !_img || !_weight ) {
      throw aol::Exception ( "MCMMassOp::getCoeff: set image first!", __FILE__, __LINE__ );
    }

    typename ConfiguratorType::VecType v;
    Grad.setZero();
    for ( int b = 0; b < this->getConfigurator().getNumLocalDofs ( El ); b++ ) {
      v = this->getBaseFunctionSet ( El ).evaluateGradient ( b, QuadPoint );
      v *= ( *_img ) [ this->getConfigurator().localToGlobal ( El, b ) ];
      Grad += v;
    }
  }

  inline RealType evaluateWeight ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    if ( !_img || !_weight ) {
      throw aol::Exception ( "MCMMassOp::getCoeff: set image first!", __FILE__, __LINE__ );
    }

    RealType w = 0.;
    for ( int b = 0; b < this->getConfigurator().getNumLocalDofs ( El ); b++ ) {
      w += ( *_weight ) [ this->getConfigurator().localToGlobal ( El, b ) ] * this->getBaseFunctionSet ( El ).evaluate ( b, QuadPoint );
    }
    return w;
  }

protected:
  const aol::Vector<RealType> *_img;
  const aol::Vector<RealType> *_weight;
  RealType _eps;
};

/***
 * class takes over the whole management and calculation of one timestep of the
 * isotropic mean curvature flow. No Boundary conditions are specified (i.e.
 * Neumann-conditions are received).
 * \author Droske
 */
template <typename ConfType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class MCMSmoothOp : public aol::Op<aol::Vector<typename ConfType::RealType> > {
public:
  typedef typename ConfType::RealType RealType;
protected:
  const qc::GridDefinition &_grid;
  mutable qc::FastUniformGridMatrix<RealType, ConfType::Dim > _mat;
  RealType _tau, _epsilon;
public:
  MCMSmoothOp ( const qc::GridDefinition &grid, RealType tau, RealType epsilon )
      : _grid ( grid ), _mat ( grid ), _tau ( tau ), _epsilon ( epsilon ) {}

  virtual ~MCMSmoothOp() { }

  void setTau ( RealType tau ) {
    _tau = tau;
  }

  void applyAdd ( const aol::Vector<typename ConfType::RealType> &arg, aol::Vector<typename ConfType::RealType> &dest ) const {
    if ( &arg == &dest ) {
      throw aol::Exception ( "dest and arg are not allowed to be the same..\n", __FILE__, __LINE__ );
    }

    qc::MCMStiffOp<ConfType, IndexMode> stiffOp ( _grid, aol::ONTHEFLY, _epsilon );
    qc::MCMMassOp<ConfType, IndexMode> massOp ( _grid, aol::ONTHEFLY, _epsilon );

    stiffOp.setImageReference ( arg );
    massOp.setImageReference ( arg );

    _mat.setZero();
    stiffOp.assembleAddMatrix ( _mat );
    _mat *= _tau;
    massOp.assembleAddMatrix ( _mat );

    aol::Vector<RealType> rhs ( dest );

    aol::SSORPreconditioner<aol::Vector<RealType>, qc::FastUniformGridMatrix<RealType, ConfType::Dim> > precond ( _mat );
    aol::PCGInverse<aol::Vector<RealType> > solver ( _mat, precond, 1e-18, 1000 );

    massOp.apply ( arg, rhs );
    solver.applyAdd ( rhs, dest );
  }
};



/**
 * MCMDirichletSmoothOp
 * class takes over the whole management and calculation of one timestep of the
 * isotropic mean curvature flow. One can specify Dirichlet conditions. Therefore
 * one has to deliver a DofMask (0 = inner node, 1 = boundary node) and a vector
 * with Dirichlet-values.
 * NOTICE: the data should be in [0,1], otherwise there might not exist solutions!
 * Author: Marc Droske, boundary stuff by Oliver Nemitz.
 */
template <typename ConfType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class MCMDirichletSmoothOp : public aol::Op<aol::Vector<typename ConfType::RealType> > {
public:
  typedef typename ConfType::RealType RealType;

protected:
  const qc::GridDefinition &_grid;
  mutable qc::UniformGridSparseMatrix<RealType> _mat;
  const aol::DofMask &_boundaryMask;
  aol::Vector<RealType> _dirichletValues;
  RealType _tau, _epsilon;

public:
  MCMDirichletSmoothOp ( const qc::GridDefinition &grid, const aol::DofMask &boundaryMask,
                         const aol::Vector<RealType> &dirichletValues, RealType tau, RealType epsilon )
      : _grid ( grid ), _mat ( grid ), _boundaryMask ( boundaryMask ), _dirichletValues ( dirichletValues ),
      _tau ( tau ), _epsilon ( epsilon ) {
    _dirichletValues.setZero();
    for ( aol::DofMask::iterator it = _boundaryMask.begin(); it != _boundaryMask.end(); ++it )
      _dirichletValues[ ( *it ) ] = dirichletValues[ ( *it ) ];
  }

virtual ~MCMDirichletSmoothOp() { }

  void setTau ( RealType tau ) {
    _tau = tau;
  }

  void applyAdd ( const aol::Vector<typename ConfType::RealType> &arg, aol::Vector<typename ConfType::RealType> &dest ) const {
    if ( &arg == &dest ) {
      throw aol::Exception ( "dest and arg are not allowed to be the same..\n", __FILE__, __LINE__ );
    }

    // the normal MCM stuff
    qc::MCMStiffOp<ConfType, IndexMode> stiffOp ( _grid, aol::ONTHEFLY, _epsilon );
    qc::MCMMassOp<ConfType, IndexMode> massOp ( _grid, aol::ONTHEFLY, _epsilon );

    stiffOp.setImageReference ( arg );
    massOp.setImageReference ( arg );

    _mat.setZero();
    stiffOp.assembleAddMatrix ( _mat );
    _mat *= _tau;
    massOp.assembleAddMatrix ( _mat );

    aol::Vector<RealType> rhs ( dest );

    massOp.apply ( arg, rhs );

    // now take the Dirichlet values into account:
    // apply the matrix to the dirichlet vector and subtract it from the rhs (Lu^{tilde} = f - Lg)
    aol::Vector<RealType> tmp ( rhs.size() );
    _mat.apply ( _dirichletValues, tmp );
    tmp *= -1.;
    rhs += tmp;

    // then traverse the boundary elements, delete the rows (columns aren't necessary) in the matrix
    // belonging to the nodes (except diagonal = 1) and put the Dirichlet value (0)
    // into the rhs-vector. (u^{tilde} = 0 on the boundary)
    for ( aol::DofMask::iterator it = _boundaryMask.begin(); it != _boundaryMask.end(); ++it ) {
      int index = ( *it );
      _mat.setRowToZero ( index );                  // clear the row and set diagonal to 1
      _mat.set ( index, index, 1. );

      rhs[index] = 0.;                          // put 0-dirichlet value into the rhs
    }

    aol::SSORPreconditioner<aol::Vector<RealType>, qc::UniformGridSparseMatrix<RealType> > precond ( _mat );
    aol::PCGInverse<aol::Vector<RealType> > solver ( _mat, precond, 1e-18, 1000 );

    solver.applyAdd ( rhs, dest );

    // finally add the dirichlet-values to the computed solution (u = u^{tilde} + g)
    for ( aol::DofMask::iterator it = _boundaryMask.begin(); it != _boundaryMask.end(); ++it )
      dest[ ( *it ) ] += _dirichletValues[ ( *it ) ];

  }
};






} // end namespace qc

#endif
