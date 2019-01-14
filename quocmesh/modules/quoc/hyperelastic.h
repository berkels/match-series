#ifndef __HYPERELASTIC_H
#define __HYPERELASTIC_H

#include <op.h>
#include <vec.h>
#include <multiVector.h>
#include <linearSmoothOp.h>
#include <FEOpInterface.h>
#include <deformations.h>
#include <regression.h>

namespace qc {

template <typename ConfiguratorType, qc::Dimension Dim>
class VolumeGradient
      : public aol::FENonlinVectorDiffOpInterface<ConfiguratorType, Dim, Dim, VolumeGradient<ConfiguratorType, Dim> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  const qc::GridDefinition &_grid;
protected:
public:
  VolumeGradient ( const qc::GridDefinition &Grid )
      : aol::FENonlinVectorDiffOpInterface<ConfiguratorType, Dim, Dim, VolumeGradient<ConfiguratorType, Dim> > ( Grid ),
      _grid ( Grid ) {};

  RealType gammaPrime ( RealType det ) const {
    return ( 2. * det - 2. );
  }

  void getNonlinearity ( const aol::auto_container<Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
                         aol::Mat<Dim,ConfiguratorType::Dim,RealType> &NL ) const {
    typename ConfiguratorType::VecType grad;
    typename ConfiguratorType::MatType dphi;

    NL.setZero();
    dphi.setZero();

    for ( int i = 0; i < Dim; i++ ) {
      DiscFuncs[i].evaluateGradientAtQuadPoint ( El, QuadPoint, grad );
      for ( int j = 0; j < Dim; j++ ) {
        dphi[i][j] = grad[j];
      }
    }

    for ( int i = 0; i < Dim; i++ ) {
      dphi[i][i] += 1.;
    }

    typename ConfiguratorType::MatType aux;
    aux.makeCofactorMatrix ( dphi );
    NL = aux;

    RealType det = dphi.det();
    NL *= gammaPrime ( det );
  }
};

template <typename ConfiguratorType, qc::Dimension Dim>
class VolumeEnergy
      : public aol::FENonlinIntegrationVectorInterface < ConfiguratorType,
      VolumeEnergy<ConfiguratorType, Dim> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  const qc::GridDefinition &_grid;
protected:
public:
  VolumeEnergy ( const qc::GridDefinition &Grid )
      : aol::FENonlinIntegrationVectorInterface < ConfiguratorType,
      VolumeEnergy<ConfiguratorType, Dim> > ( Grid ),
      _grid ( Grid ) {};

  RealType evaluateIntegrand ( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &discrFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/ ) const {
    typename ConfiguratorType::VecType grad;
    typename ConfiguratorType::MatType dphi;

    dphi.setZero();

    for ( int i = 0; i < Dim; i++ ) {
      discrFuncs[i].evaluateGradientAtQuadPoint ( El, QuadPoint, grad );
      for ( int j = 0; j < Dim; j++ ) {
        dphi[i][j] = grad[j];
      }
    }

    for ( int j = 0; j < Dim; j++ ) {
      dphi[j][j] += 1.;
    }

    RealType det = dphi.det();
    return gamma ( det );
  }
protected:
  RealType gamma ( RealType det ) const {
    return aol::Sqr ( det - 1. );
  }
};


template <typename ConfiguratorType, qc::Dimension Dim>
class LengthGradient
      : public aol::FENonlinVectorDiffOpInterface<ConfiguratorType, Dim, Dim, LengthGradient<ConfiguratorType, Dim> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  const qc::GridDefinition &_grid;
protected:
public:
  LengthGradient ( const qc::GridDefinition &Grid )
      : aol::FENonlinVectorDiffOpInterface<ConfiguratorType, Dim, Dim, LengthGradient<ConfiguratorType, Dim> > ( Grid ),
      _grid ( Grid ) {};

  void getNonlinearity ( const aol::auto_container<Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/,
                         aol::Mat<Dim,ConfiguratorType::Dim,RealType> &NL ) const {
    typename ConfiguratorType::VecType grad;

    for ( int i = 0; i < Dim; i++ ) {
      DiscFuncs[i].evaluateGradientAtQuadPoint ( El, QuadPoint, grad );
      for ( int j = 0; j < Dim; j++ ) {
        NL[i][j] = grad[j];
      }
    }
  }
};


template <typename ConfiguratorType>
class LengthEnergy
      : public aol::FENonlinIntegrationVectorInterface < ConfiguratorType,
      LengthEnergy<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
protected:
public:
  LengthEnergy ( const typename ConfiguratorType::InitType &Grid )
      : aol::FENonlinIntegrationVectorInterface < ConfiguratorType,
      LengthEnergy<ConfiguratorType> > ( Grid ),
      _grid ( Grid ) {};

  RealType evaluateIntegrand ( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &discrFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::VecType &/*RefCoord*/ ) const {
    typename ConfiguratorType::VecType grad;

    RealType energy = 0.;

    for ( int i = 0; i < ConfiguratorType::Dim; ++i ) {
      discrFuncs[i].evaluateGradientAtQuadPoint ( El, QuadPoint, grad );
      energy += grad.normSqr();
    }
    return 0.5 * energy;
  }
};

/**
 * Does the same as LengthEnergy, but uses a StiffOp (has to be supplied as argument) to do so.
 * instead of using aol::FENonlinIntegrationVectorInterface. If the supplied StiffOp is assembled,
 * applying DisplacementLengthEnergy is much faster than applying LengthEnergy.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class DisplacementLengthEnergy: public aol::StandardGenEnergyOp<aol::MultiVector<typename ConfiguratorType::RealType> >
{
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiagonalBlockOp<RealType> _blockStiff;
public:

  DisplacementLengthEnergy( const aol::StiffOp<ConfiguratorType> &Stiff )
    : _blockStiff ( Stiff ) {
  }
  virtual ~DisplacementLengthEnergy() {}
  virtual void apply( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::MultiVector <RealType> mtmp( Arg, aol::STRUCT_COPY );
    _blockStiff.apply( Arg, mtmp); // E = 1/2 L_dPhi*Phi
    Dest[0] = (0.5)*(Arg * mtmp);
    this->_lastEnergy = Dest[0];
  }
  virtual void applyAdd( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest) const {
    aol::Scalar<RealType> tmp ( Dest );
    apply ( Arg, tmp );
    Dest += tmp;
    this->_lastEnergy = tmp;
  }
};

/**
 * \brief To use the Dirichlet energy in Gauss-Newton approaches. Only provides one component of the gradient.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename MatrixType = aol::SparseMatrix<typename ConfiguratorType::RealType> >
class DirichtletEnergyComponentGNFunc : public aol::LinearGNFuncBase<typename ConfiguratorType::RealType, MatrixType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const ConfiguratorType _config;

  void assembleDirichtletEnergyComponentMatrix ( MatrixType &Mat, const int GradComp, const RealType Lambda ) {
    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const typename IteratorType::EndType end_it = _config.end();
    int i=0;
    for ( IteratorType it = _config.begin(); it != end_it; ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );
      const int numLocalDofs = _config.getNumLocalDofs ( *it );
      const RealType vol = _config.vol ( *it );
      for ( int dof = 0; dof < numLocalDofs; dof++ )
        for ( int q = 0; q < numQuadPoints; ++q )
          Mat.set(i+q, _config.localToGlobal ( *it, dof ), Lambda * sqrt( bfs.getWeight ( q ) * vol ) * (bfs.evaluateGradient (dof, q))[GradComp] );
      i+=numQuadPoints;
    }
  }
public:
  DirichtletEnergyComponentGNFunc ( const typename ConfiguratorType::InitType &Grid,
                                    const int GradComp,
                                    const RealType Lambda = aol::ZOTrait<RealType>::one )
    : aol::LinearGNFuncBase<RealType, MatrixType> ( 0, 0 ),
      _config ( Grid ) {
    this->_mat.reallocate ( aol::totalNumberOfQuadPoints ( _config ), Grid.getNumberOfNodes() );
    assembleDirichtletEnergyComponentMatrix ( this->_mat, GradComp, Lambda );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
void assembleLaplaceSquareRegMatrix ( const typename ConfiguratorType::InitType &Grid, aol::SparseMatrix<typename ConfiguratorType::RealType> &RegMat ) {

  aol::SparseMatrix<typename ConfiguratorType::RealType> tempMatrix1( Grid ); // going to be L
  aol::SparseMatrix<typename ConfiguratorType::RealType> tempMatrix2( Grid ); // going to be M^{-1}L

  aol::LumpedMassOp<ConfiguratorType> lumpedMassInv( Grid, aol::INVERT );
  aol::StiffOp<ConfiguratorType> stiff( Grid );
  stiff.assembleAddMatrix( tempMatrix1 );
  tempMatrix2 = tempMatrix1;

  const int numNodes = Grid.getNumberOfNodes();

  for ( int i = 0; i < numNodes; ++i )
    tempMatrix2.scaleRow ( i, lumpedMassInv.getMatrix( )[ i ] );

  RegMat.reallocate( numNodes, numNodes );
  RegMat.addMatrixProduct( tempMatrix1, tempMatrix2 ); // regMat = LM^{-1}L
}

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class LaplaceEnergy : public aol::StandardGenEnergyOp<aol::MultiVector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  aol::SparseMatrix<RealType> _mat;
  aol::DiagonalBlockOp<RealType> _diagBlockOp;
  aol::QuadraticFormOp<aol::MultiVector<RealType> > _q;
public:
  LaplaceEnergy ( const typename ConfiguratorType::InitType &Initializer )
    : _mat ( Initializer ),
      _diagBlockOp ( _mat ),
      _q ( _diagBlockOp ) {
    qc::assembleLaplaceSquareRegMatrix<ConfiguratorType> ( Initializer, _mat );
  }

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    this->_lastEnergy = Dest[0];
    _q.applyAdd ( MArg, Dest );
    this->_lastEnergy = Dest[0] - this->_lastEnergy;
  }

  void applyDerivative ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    _diagBlockOp.apply ( MArg, MDest );
  }
};

/**
 * Like qc::LaplaceEnergy but uses an aol::CSR_Matrix<> matrix instead of an aol::SparseMatrix<RealType>
 * to store the assembled operator. Thus, it is considerably faster than qc::LaplaceEnergy but only works
 * with double precision.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class LaplaceEnergyCSR : public aol::StandardGenEnergyOp<aol::MultiVector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  aol::DeleteFlagPointer<aol::CSR_Matrix<> >  _pMat;
  aol::DeleteFlagPointer<aol::DiagonalBlockOp<RealType> > _pDiagBlockOp;
  aol::DeleteFlagPointer<aol::QuadraticFormOp<aol::MultiVector<RealType> > > _pQ;
public:
  LaplaceEnergyCSR ( const typename ConfiguratorType::InitType &Initializer ) {
    aol::SparseMatrix<RealType> tempMat ( Initializer );
    qc::assembleLaplaceSquareRegMatrix<ConfiguratorType> ( Initializer, tempMat );
    _pMat.reset ( new aol::CSR_Matrix<> ( tempMat ), true );
    _pDiagBlockOp.reset ( new aol::DiagonalBlockOp<RealType> ( *_pMat ), true );
    _pQ.reset ( new aol::QuadraticFormOp<aol::MultiVector<RealType> > ( *_pDiagBlockOp ), true );
  }
    
  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    this->_lastEnergy = Dest[0];
    _pQ->applyAdd ( MArg, Dest );
    this->_lastEnergy = Dest[0] - this->_lastEnergy;
  }
    
  void applyDerivative ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    _pDiagBlockOp->apply ( MArg, MDest );
  }

  const aol::CSR_Matrix<> &getMatrixRef ( ) const {
    return *_pMat;
  }
};

/**
 * \brief To use the Laplace energy in Gauss-Newton approaches.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename MatrixType = aol::SparseMatrix<typename ConfiguratorType::RealType> >
class LaplaceEnergyGNFunc : public aol::LinearGNFuncBase<typename ConfiguratorType::RealType, MatrixType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const typename ConfiguratorType::InitType &_grid;
public:
  LaplaceEnergyGNFunc ( const typename ConfiguratorType::InitType &Grid,
                        const RealType Lambda = aol::ZOTrait<RealType>::one )
    : aol::LinearGNFuncBase<RealType, MatrixType> ( Grid.getNumberOfNodes(), Grid.getNumberOfNodes() ),
      _grid ( Grid ) {
    const aol::LumpedMassOp<ConfiguratorType> lumpedMassOp ( _grid, aol::INVERT );
    aol::DiagonalMatrix<RealType> MInv ( _grid );
    lumpedMassOp.assembleAddMatrix ( MInv );
    const aol::StiffOp<ConfiguratorType> stiffOp ( _grid, aol::ONTHEFLY );
    stiffOp.assembleAddMatrix ( this->_mat );
    for ( int i = 0; i < this->_mat.getNumRows(); i++)
      this->_mat.scaleRow ( i, Lambda * MInv.getDiag(i) );
  }
};

/**
  * \brief This class computes via "apply(...)" or "applyAdd(...)" a hyperelastic energy \f$ \int_\Omega W(I_1,I_2,I_3,x) dx \f$,
  * where the hyperelastic invariants are given by \f$ I_1=||\nabla\phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi)||_F^2 \f$,
  * \f$ I_1=det(\nabla\phi) \f$ and the function \f$ W(I_1,I_2,I_3,x) \f$ is passed to the constructor (e.g. HyerelasticEnergyDensityDefault).
  * In its first argument, "apply" and "applyAdd" expect a Multivector which represents the displacements \f$ d_i \f$
  * in the different Cartesian directions; the deformation \f$ \phi \f$ is then automatically computed via \f$ \phi=d+\f$ identity.
  *
  * \author Wirth
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class HyperelasticEnergy
      : public aol::FENonlinIntegrationVectorInterface < ConfiguratorType, HyperelasticEnergy<ConfiguratorType,HyperelasticEnergyDensityType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;
public:
  HyperelasticEnergy ( const typename ConfiguratorType::InitType &Grid, const HyperelasticEnergyDensityType &HyperelasticEnergyDensity )
      : aol::FENonlinIntegrationVectorInterface < ConfiguratorType, HyperelasticEnergy<ConfiguratorType,HyperelasticEnergyDensityType> > ( Grid ),
        _hyperelasticEnergyDensity( HyperelasticEnergyDensity ) {}

  /**
    * \brief Returns the hyperelastic energy \f$ W(I_1,I_2,I_3) \f$ for the invariants \f$ I_1=||\nabla\phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi)||_F^2\f$, \f$ I_3=det(\nabla\phi) \f$.
    */
  RealType evaluateIntegrand ( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El, int QuadPoint,
                               const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    // compute the deformation gradient $\nabla\phi$: for each displacment component do...
    typename ConfiguratorType::MatType dphi;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      // evaluate gradient of the displacement $d_i$
      DiscFuncs[i].evaluateGradientAtQuadPoint ( El, QuadPoint, dphi[i] );
      // compute corresponding deformation gradient of the $i$th component
      dphi[i][i] += 1.;
    }
    return energyDensity( dphi, _hyperelasticEnergyDensity, El, QuadPoint );
  }

  /**
    * \brief Returns the hyperelastic energy density \f$ W(I_1,I_2,I_3) \f$ for the invariants
    * \f$ I_1=||\nabla\phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi)||_F^2\f$, \f$ I_3=det(\nabla\phi) \f$,
    * where \f$ \nabla\phi \f$, \f$ W \f$, and the spatial position are passed as arguments.
    */
  static inline RealType energyDensity ( const typename ConfiguratorType::MatType &DPhi,
                                         const HyperelasticEnergyDensityType &EnergyDensity,
                                         const typename ConfiguratorType::ElementType &El, int QuadPoint ) {
    // compute the hyperelastic invariants $||\nabla\phi||_F^2,||cof(\nabla\phi)||_F^2,det(\nabla\phi)$
    RealType
      I1 = DPhi.normSqr(),
      I2,
      I3 = DPhi.det();
    if ( ConfiguratorType::Dim != qc::QC_2D ) {
      typename ConfiguratorType::MatType cof;
      cof.makeCofactorMatrix( DPhi );
      I2 = cof.normSqr();
    } else
      // in 2d, $||cof(\nabla\phi)||_F$ equals $||\nabla\phi||_F$
      I2 = I1;
    return EnergyDensity.evaluateAtQuadPoint ( I1, I2, I3, El, QuadPoint );
  }
};

/**
  * \brief This class computes via "apply(...)" or "applyAdd(...)" a hyperelastic energy gradient with respect to the displacement,
  * i.e. \f$ (<\frac{\partial \int_\Omega W(I_1,I_2,I_3,x) dx}{\partial \phi},\psi_i>)_i \f$, where the \f$ \psi_i \f$ are the FE basis
  * functions, the hyperelastic invariants are given by \f$ I_1=||\nabla\phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi)||_F^2 \f$,
  * \f$ I_3=det(\nabla\phi)\f$ and the function \f$ W(I_1,I_2,I_3,x) \f$ is passed to the constructor (e.g. HyerelasticEnergyDensityDefault).
  * In its first argument, "apply" and "applyAdd" expect a Multivector which represents the displacements \f$ d_i \f$
  * in the different Cartesian directions; the deformation \f$ \phi \f$ is then automatically computed via \f$ \phi=d+\f$ identity.
  *
  * \author Wirth
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class HyperelasticGradient
      : public aol::FENonlinVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, HyperelasticGradient<ConfiguratorType,HyperelasticEnergyDensityType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;
public:
  HyperelasticGradient ( const typename ConfiguratorType::InitType &Grid, const HyperelasticEnergyDensityType &HyperelasticEnergyDensity )
      : aol::FENonlinVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, HyperelasticGradient<ConfiguratorType,HyperelasticEnergyDensityType> > ( Grid ),
        _hyperelasticEnergyDensity( HyperelasticEnergyDensity ) {}

  /**
    * \brief Returns \f$ 2\frac{\partial W}{\partial I_1}\nabla\phi+2\frac{\partial W}{\partial I_2}\frac{\partial cof(\nabla\phi)}{\partial \phi}+\frac{\partial W}{\partial I_3}cof(\nabla\phi) \f$,
    * where \f$ \frac{\partial cof(\nabla\phi)}{\partial \phi} \f$ shall be short for the matrix \f$ M \f$ such that
    * \f$ tr(M^T\nabla\psi) = tr[cof(\nabla\phi)cof(\nabla\phi)^T(tr((\nabla\phi)^{-1}\nabla\psi)-\nabla\psi(\nabla\phi)^{-1})] \f$
    * for any \f$ \psi \f$.
    */
  void getNonlinearity ( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/,
                         aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim,RealType> &NL ) const {
    // compute the deformation gradient $\nabla\phi$: for each displacment component do...
    typename ConfiguratorType::MatType dPhi;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      // evaluate gradient of the displacement $d_i$
      DiscFuncs[i].evaluateGradientAtQuadPoint ( El, QuadPoint, dPhi[i] );
      // compute corresponding deformation gradient of the $i$th component
      dPhi[i][i] += 1.;
    }

    firstPiolaKirchhoffStress( dPhi, _hyperelasticEnergyDensity, El, QuadPoint, NL );
  }

  /**
    * \brief Computes Stress=\f$ \frac{\partial W}{\partial\nabla\phi}(I_1,I_2,I_3) \f$ for the invariants
    * \f$ I_1=||\nabla\phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi)||_F^2\f$, \f$ I_3=det(\nabla\phi) \f$,
    * where \f$ \nabla\phi \f$, \f$ W \f$, and the spatial position are passed as arguments.
    */
  static inline void firstPiolaKirchhoffStress ( const typename ConfiguratorType::MatType &DPhi,
                                                 const HyperelasticEnergyDensityType &EnergyDensity,
                                                 const typename ConfiguratorType::ElementType &El, int QuadPoint,
                                                 aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim,RealType> &Stress ) {
    Stress.setZero();

    // compute the cofactor matrix of the deformation gradient, $cof(\nabla\phi)$
    typename ConfiguratorType::MatType dPhi( DPhi ), cof;
    cof.makeCofactorMatrix( DPhi );
    // compute the hyperelastic invariants $||\nabla\phi||_F^2,||cof(\nabla\phi)||_F^2,det(\nabla\phi)$
    RealType
      I1 = DPhi.normSqr(),
      I2 = cof.normSqr(),
      I3 = DPhi.det();

    // compute the outer derivative of the hyperelastic energy
    aol::Vec3<RealType> outerDeriv ( EnergyDensity.evaluateDerivativeAtQuadPoint ( I1, I2, I3, El, QuadPoint ) );

    // compute the first Piola-Kirchhoff stress
    // the term belonging to $I_2$...
    if ( ConfiguratorType::Dim == qc::QC_3D ) {
      Stress[0][0] = DPhi[1][1] * cof[2][2] + DPhi[2][2] * cof[1][1] - DPhi[1][2] * cof[2][1] - DPhi[2][1] * cof[1][2];
      Stress[0][1] = DPhi[1][2] * cof[2][0] + DPhi[2][0] * cof[1][2] - DPhi[1][0] * cof[2][2] - DPhi[2][2] * cof[1][0];
      Stress[0][2] = DPhi[2][1] * cof[1][0] + DPhi[1][0] * cof[2][1] - DPhi[2][0] * cof[1][1] - DPhi[1][1] * cof[2][0];
      Stress[1][0] = DPhi[0][2] * cof[2][1] + DPhi[2][1] * cof[0][2] - DPhi[0][1] * cof[2][2] - DPhi[2][2] * cof[0][1];
      Stress[1][1] = DPhi[0][0] * cof[2][2] + DPhi[2][2] * cof[0][0] - DPhi[0][2] * cof[2][0] - DPhi[2][0] * cof[0][2];
      Stress[1][2] = DPhi[2][0] * cof[0][1] + DPhi[0][1] * cof[2][0] - DPhi[2][1] * cof[0][0] - DPhi[0][0] * cof[2][1];
      Stress[2][0] = DPhi[0][1] * cof[1][2] + DPhi[1][2] * cof[0][1] - DPhi[0][2] * cof[1][1] - DPhi[1][1] * cof[0][2];
      Stress[2][1] = DPhi[1][0] * cof[0][2] + DPhi[0][2] * cof[1][0] - DPhi[1][2] * cof[0][0] - DPhi[0][0] * cof[1][2];
      Stress[2][2] = DPhi[0][0] * cof[1][1] + DPhi[1][1] * cof[0][0] - DPhi[0][1] * cof[1][0] - DPhi[1][0] * cof[0][1];
      Stress *= outerDeriv[1] * 2;
      // alternative implementation
      /*if ( I3 != 0. ) {
        Stress.makeProductABtransposed( cof, cof );
        RealType trace = Stress.tr();
        for ( int i = 0; i < ConfiguratorType::Dim; i++ )
          Stress[i][i] -= trace;
        Stress *= cof;
        Stress *= - outerDeriv[1] * 2 / I3;
      }*/
    } else
      // in 2d, the term belonging to $I_2$ equals the one belonging to $I_1$
      outerDeriv[0] += outerDeriv[1];
    // the term belonging to $I_1$...
    dPhi *= 2 * outerDeriv[0];
    Stress += dPhi;
    // the term belonging to $I_3$...
    cof *= outerDeriv[2];
    Stress += cof;
  }
};

/**
  * \brief This class computes via "apply(...)" or "applyAdd(...)" parts of a hyperelastic energy Hessian with respect to the displacement,
  * i.e. \f$ (<\frac{\partial^2 \int_\Omega W(I_1,I_2,I_3,x) dx}{\partial \phi^2},\psi_i,\psi_j>)_{i,j} \f$,
  * where the \f$ \psi_i \f$ and \f$ \psi_i \f$ are those vector-valued FE basis functions that are only nonzero in the kth and lth component,
  * the hyperelastic invariants are given by \f$ I_1=||\nabla\phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi)||_F^2 \f$,
  * \f$ I_3=det(\nabla\phi)\f$, and the function \f$ W(I_1,I_2,I_3,x) \f$ is passed to the constructor (e.g. HyerelasticEnergyDensityDefault).
  * In its first argument, "apply" and "applyAdd" expect a Multivector which represents the displacements \f$ d_i \f$
  * in the different Cartesian directions; the deformation \f$ \phi \f$ is then automatically computed via \f$ \phi=d+\f$ identity.
  *
  * \author Wirth
  * \ingroup MatrixFEOp
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class HyperelasticSubHessian
      : public aol::FELinAsymMatrixWeightedStiffInterface<ConfiguratorType, HyperelasticSubHessian<ConfiguratorType,HyperelasticEnergyDensityType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;
  // the displacement, at which the sub-Hessian shall be evaluated
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement;
  // the indices of the sub-Hessian
  const int _k, _l;
public:
  HyperelasticSubHessian ( const typename ConfiguratorType::InitType &Grid,
                           const HyperelasticEnergyDensityType &HyperelasticEnergyDensity,
                           const aol::MultiVector<RealType> &Displacement,
                           const int K,
                           const int L )
      : aol::FELinAsymMatrixWeightedStiffInterface<ConfiguratorType, HyperelasticSubHessian<ConfiguratorType,HyperelasticEnergyDensityType> > ( Grid ),
        _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
        _displacement( Grid, Displacement ),
        _k( K ),
        _l( L ) {}

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                               const typename ConfiguratorType::DomVecType& /*DomRefCoord*/, typename ConfiguratorType::MatType &Matrix ) const {
    // compute the deformation gradient $\nabla\phi$
    typename ConfiguratorType::MatType dPhi;
    _displacement.evaluateGradientAtQuadPoint ( El, QuadPoint, dPhi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dPhi[i][i] += 1.;

    // compute the Hessian matrix for the kth component multiplied from the left and the lth from the right
    elasticityTensor( dPhi, _hyperelasticEnergyDensity, El, QuadPoint, _k, _l, Matrix );
  }

  /**
    * \brief Computes ElasticityTensor=\f$ \frac{\partial^2W}{\partial\nabla\phi_K\partial\nabla\phi_L}(I_1,I_2,I_3) \f$ for the invariants
    * \f$ I_1=||\nabla\phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi)||_F^2\f$, \f$ I_3=det(\nabla\phi) \f$,
    * where \f$ \nabla\phi \f$, \f$ W \f$, and the spatial position are passed as arguments.
    * The result matrix \f$ C= \f$ ElasticityTensor is such that
    * \f$ <\frac{\partial^2W}{\partial\nabla\phi_K\partial\nabla\phi_L},\theta_K,\psi_L>=\theta_K\cdot(C\psi_L) \f$
    * for a variation \f$ \theta_K \f$ of \f$ \nabla\phi_K \f$ and \f$ \psi_L \f$ of \f$ \nabla\phi_L \f$.
    */
  static inline void elasticityTensor( const typename ConfiguratorType::MatType &DPhi,
                                       const HyperelasticEnergyDensityType &EnergyDensity,
                                       const typename ConfiguratorType::ElementType &El, int QuadPoint,
                                       const int K, const int L,
                                       aol::Mat< ConfiguratorType::Dim, ConfiguratorType::Dim, RealType > &ElasticityTensor ) {
    // the cofactor matrix of the deformation gradient
    typename ConfiguratorType::MatType cof;
    // auxiliary vectors and matrices
    typename ConfiguratorType::VecType auxVec1;
    typename ConfiguratorType::MatType auxMat1, auxMat2;

    ElasticityTensor.setZero();

    // compute the cofactor matrix of the deformation gradient, $cof(\nabla\phi)$
    cof.makeCofactorMatrix( DPhi );
    // compute the hyperelastic invariants $||\nabla\phi||_F^2,||cof(\nabla\phi)||_F^2,det(\nabla\phi)$
    RealType
      I1 = DPhi.normSqr(),
      I2 = cof.normSqr(),
      I3 = DPhi.det();

    // compute the outer first and second derivative of the hyperelastic energy
    aol::Vec3<RealType> outerDeriv ( EnergyDensity.evaluateDerivativeAtQuadPoint ( I1, I2, I3, El, QuadPoint ) );
    aol::Matrix33<RealType> outerSecondDeriv ( EnergyDensity.evaluateSecondDerivativeAtQuadPoint ( I1, I2, I3, El, QuadPoint ) );

    // compute auxiliary matrices h, H
    typename ConfiguratorType::MatType H( cof.transposed() ), h, hCof;
    h.makeProduct( cof, H );
    RealType auxScal1 = h.tr();
    auxMat1 = h;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      auxMat1[i][i] -= auxScal1;
    H.makeProduct( auxMat1, cof );
    H /= - I3;
    hCof.makeProduct( h, cof );

    // compute the Hessian matrix for the kth component multiplied from the left and the lth from the right
    if ( K == L )
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        ElasticityTensor[i][i] = 2 * outerDeriv[0]; // note: \epsilon:\varespilon = \sum_{k,l=1..d}\delta_{kl}\epsilon_{k:}I\varepsilon_{l:}^T
    auxVec1.setMultiple( DPhi[K], 4 * outerSecondDeriv[0][0] );
    auxVec1.addMultiple( H[K], 4 * outerSecondDeriv[0][1] );
    auxVec1.addMultiple( cof[K], 2 * outerSecondDeriv[0][2] );
    auxMat1.makeTensorProduct( auxVec1, DPhi[L] );
    ElasticityTensor += auxMat1; // note: (A:\epsilon)(B:\varepsilon) = \sum_{k,l=1..d}\epsilon_{k:}A_{k:}^T B_{l:}\varepsilon_{l:}^T

    auxVec1.setMultiple( cof[L], - outerDeriv[2] / I3 );
    auxMat1.makeTensorProduct( auxVec1, cof[K] );
    ElasticityTensor += auxMat1; // note: tr(A\epsilon^T B\varepsilon^T) = \sum_{k,l=1..d}\epsilon_{k:}A_{l:}^T B_{k:}\varepsilon_{l:}^T
    auxVec1.setMultiple( cof[K], outerDeriv[2] / I3 );
    auxVec1.addMultiple( DPhi[K], 2 * outerSecondDeriv[2][0] );
    auxVec1.addMultiple( H[K], 2 * outerSecondDeriv[2][1] );
    auxVec1.addMultiple( cof[K], outerSecondDeriv[2][2] );
    auxMat1.makeTensorProduct( auxVec1, cof[L] );
    ElasticityTensor += auxMat1;

    auxVec1.setMultiple( cof[K], - 2 * outerDeriv[1] / aol::Sqr( I3 ) );
    auxMat1.makeTensorProduct( auxVec1, hCof[L] );
    ElasticityTensor += auxMat1;
    auxVec1.setMultiple( cof[L], 2 * outerDeriv[1] / aol::Sqr( I3 ) );
    auxMat1.makeTensorProduct( auxVec1, hCof[K] );
    ElasticityTensor += auxMat1;
    auxMat2.transposeFrom( cof );
    auxMat1.makeProduct( auxMat2, cof );
    auxMat1 *= 2 * outerDeriv[1] / aol::Sqr( I3 ) * h[L][K];
    ElasticityTensor += auxMat1; // note: tr(A\epsilon B\varepsilon^T) = \sum_{k,l=1..d}\epsilon_{k:}A_{lk}B\varepsilon_{l:}^T
    auxVec1.setMultiple( cof[K],  2 * outerDeriv[1] / aol::Sqr( I3 ) * h.tr() );
    auxVec1.addMultiple( hCof[K], - 4 * outerDeriv[1] / aol::Sqr( I3 ) );
    auxMat1.makeTensorProduct( auxVec1, cof[L] );
    ElasticityTensor += auxMat1;
    auxVec1.setMultiple( H[L], - 2 * outerDeriv[1] / I3 );
    auxMat1.makeTensorProduct( auxVec1, cof[K] );
    ElasticityTensor += auxMat1;
    auxVec1.setMultiple( cof[K], 2 * outerDeriv[1] / I3 );
    auxVec1.addMultiple( DPhi[K], 4 * outerSecondDeriv[1][0] );
    auxVec1.addMultiple( H[K], 4 * outerSecondDeriv[1][1] );
    auxVec1.addMultiple( cof[K], 2 * outerSecondDeriv[1][2] );
    auxMat1.makeTensorProduct( auxVec1, H[L] );
    ElasticityTensor += auxMat1;
  }
};

/**
  * \brief This class computes via "apply(...)" or "applyAdd(...)" a hyperelastic energy Hessian with respect to the displacement,
  * i.e. \f$ (<\frac{\partial^2 \int_\Omega W(I_1,I_2,I_3,x) dx}{\partial \phi^2},\psi_i,\psi_j>)_{i,j} \f$, where the \f$ \psi_i \f$ are the FE basis
  * functions, the hyperelastic invariants are given by \f$ I_1=||\nabla\phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi)||_F^2 \f$,
  * \f$ I_3=det(\nabla\phi)\f$ and the function \f$ W(I_1,I_2,I_3,x) \f$ is passed to the constructor (e.g. HyerelasticEnergyDensityDefault).
  * In its first argument, "apply" and "applyAdd" expect a Multivector which represents the displacements \f$ d_i \f$
  * in the different Cartesian directions; the deformation \f$ \phi \f$ is then automatically computed via \f$ \phi=d+\f$ identity.
  *
  * \author Wirth
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType, typename SubMatrixType = typename ConfiguratorType::MatrixType, typename BlockOpType = aol::BlockOpBase<typename ConfiguratorType::RealType, SubMatrixType> >
class HyperelasticHessian
      : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, BlockOpType > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the underlying FE grid
  const typename ConfiguratorType::InitType &_grid;
  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;

public:
  HyperelasticHessian ( const typename ConfiguratorType::InitType &Grid,
                        const HyperelasticEnergyDensityType &HyperelasticEnergyDensity )
      : _grid ( Grid ),
        _hyperelasticEnergyDensity( HyperelasticEnergyDensity ) {}

  void applyAdd( const aol::MultiVector<RealType> &Arg, BlockOpType &Dest ) const {
    for ( int k = 0; k < ConfiguratorType::Dim; k++ )
      for ( int l = 0; l < ConfiguratorType::Dim; l++ )
        HyperelasticSubHessian<ConfiguratorType,HyperelasticEnergyDensityType>( _grid, _hyperelasticEnergyDensity, Arg, k, l ).assembleAddMatrix( Dest.getReference( k, l ) );
  }
};

/**
  * \brief This class computes via "apply(...)" or "applyAdd(...)" a hyperelastic energy \f$ \int_\Omega W(I_1,I_2,I_3,x) dx \f$,
  * where the hyperelastic invariants are given by \f$ I_1=||\nabla\phi\nabla\Phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi\nabla\Phi)||_F^2 \f$,
  * \f$ I_1=det(\nabla\phi\nabla\Phi) \f$ and the function \f$ W(I_1,I_2,I_3,x) \f$ is passed to the constructor (e.g. HyerelasticEnergyDensityDefault).
  * \f$ \Phi \f$ is a deformation passed to the constructor.
  * In its first argument, "apply" and "applyAdd" expect a Multivector which represents the displacements \f$ d_i \f$
  * in the different Cartesian directions; the deformation \f$ \phi \f$ is then automatically computed via \f$ \phi=d+\f$ identity.
  *
  * \author Wirth
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class HyperelasticDeformEnergy
      : public FENonlinDeformIntegrationVectorInterface < ConfiguratorType, HyperelasticDeformEnergy<ConfiguratorType,HyperelasticEnergyDensityType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;
  // the deformation \Phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement;
public:
  HyperelasticDeformEnergy ( const typename ConfiguratorType::InitType &Grid,
                             const HyperelasticEnergyDensityType &HyperelasticEnergyDensity,
                             const aol::MultiVector<RealType> &Displacement ) :
    FENonlinDeformIntegrationVectorInterface < ConfiguratorType, HyperelasticDeformEnergy<ConfiguratorType,HyperelasticEnergyDensityType> > ( Grid, Displacement ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacement( Grid, Displacement ) {}

  RealType evaluateIntegrand ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &/*LocalCoord*/,
                               const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord ) const {
    // compute the deformation gradient $\nabla\Phi$
    typename ConfiguratorType::MatType dPhi;
    _displacement.evaluateGradientAtQuadPoint( El, QuadPoint, dPhi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dPhi[i][i] += 1.;

    // compute the deformation gradient $\nabla\phi$
    typename ConfiguratorType::MatType dphi;
    DiscFuncs.evaluateGradient( TransformedEl, TransformedLocalCoord, dphi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dphi[i][i] += 1.;

    dphi *= dPhi;
    return HyperelasticEnergy<ConfiguratorType,HyperelasticEnergyDensityType>::energyDensity( dphi, _hyperelasticEnergyDensity, El, QuadPoint );
  }
};

/**
  * \brief This class computes via "apply(...)" or "applyAdd(...)" a hyperelastic energy gradient with respect to the displacement,
  * i.e. \f$ (<\frac{\partial \int_\Omega W(I_1,I_2,I_3,x) dx}{\partial \phi},\psi_i>)_i \f$, where the \f$ \psi_i \f$ are the FE basis
  * functions, the hyperelastic invariants are given by \f$ I_1=||\nabla\phi\nabla\Phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi\nabla\Phi)||_F^2 \f$,
  * \f$ I_3=det(\nabla\phi\nabla\Phi)\f$ and the function \f$ W(I_1,I_2,I_3,x) \f$ is passed to the constructor (e.g. HyerelasticEnergyDensityDefault).
  * \f$ \Phi \f$ is some prestressing deformation.
  * In its first argument, "apply" and "applyAdd" expect a Multivector which represents the displacements \f$ d_i \f$
  * in the different Cartesian directions; the deformation \f$ \phi \f$ is then automatically computed via \f$ \phi=d+\f$ identity.
  *
  * \author Wirth
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class HyperelasticDeformGradient
      : public FENonlinDeformVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, HyperelasticDeformGradient<ConfiguratorType,HyperelasticEnergyDensityType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;
  // the deformation \Phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement;
public:
  HyperelasticDeformGradient ( const typename ConfiguratorType::InitType &Grid,
                               const HyperelasticEnergyDensityType &HyperelasticEnergyDensity,
                               const aol::MultiVector<RealType> &Displacement ) :
    FENonlinDeformVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, HyperelasticDeformGradient<ConfiguratorType,HyperelasticEnergyDensityType> > ( Grid, Displacement ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacement( Grid, Displacement ) {}

  void getNonlinearity ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &/*LocalCoord*/,
                         const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                         aol::Mat<ConfiguratorType::Dim, ConfiguratorType::Dim, RealType> &NL ) const {
    // compute the deformation gradient $\nabla\Phi$
    typename ConfiguratorType::MatType dPhi;
    _displacement.evaluateGradientAtQuadPoint( El, QuadPoint, dPhi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dPhi[i][i] += 1.;

    // compute the deformation gradient $\nabla\phi$
    typename ConfiguratorType::MatType dphi;
    DiscFuncs.evaluateGradient( TransformedEl, TransformedLocalCoord, dphi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dphi[i][i] += 1.;

    dphi *= dPhi;
    HyperelasticGradient<ConfiguratorType,HyperelasticEnergyDensityType>::firstPiolaKirchhoffStress( dphi, _hyperelasticEnergyDensity, El, QuadPoint, NL );
  }
};

/**
  * \brief This class computes via "apply(...)" or "applyAdd(...)" parts of a hyperelastic energy Hessian with respect to the displacement,
  * where the material is prestressed by a deformation \f$ \Phi \f$,
  * i.e. \f$ (<\frac{\partial^2 \int_\Omega W(I_1,I_2,I_3,x) dx}{\partial\phi^2},\psi_i\circ\Phi,\psi_j\circ\Phi>)_{i,j} \f$,
  * where the \f$ \psi_i \f$ and \f$ \psi_j \f$ are those vector-valued FE basis functions that are only nonzero in the kth and lth component,
  * the hyperelastic invariants are given by \f$ I_1=||\nabla\phi\nabla\Phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi\nabla\Phi)||_F^2 \f$,
  * \f$ I_3=det(\nabla\phi\nabla\Phi)\f$, and the function \f$ W(I_1,I_2,I_3,x) \f$ is passed to the constructor (e.g. HyerelasticEnergyDensityDefault).
  * In its first argument, "apply" and "applyAdd" expect a Multivector which represents the displacements \f$ d_i \f$
  * in the different Cartesian directions; the deformation \f$ \phi \f$ is then automatically computed via \f$ \phi=d+\f$ identity.
  *
  * \author Wirth
  * \ingroup MatrixFEOp
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class HyperelasticDeformSubHessian
      : public qc::FELinDeformAsymMatrixWeightedStiffInterface<ConfiguratorType, HyperelasticDeformSubHessian<ConfiguratorType,HyperelasticEnergyDensityType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;
  // the deformations \Phi and \phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement, _argDisplacement;
  // the indices of the sub-Hessian
  const int _k, _l;

public:
  HyperelasticDeformSubHessian ( const typename ConfiguratorType::InitType &Grid,
                                 const HyperelasticEnergyDensityType &HyperelasticEnergyDensity,
                                 const aol::MultiVector<RealType> &Displacement,
                                 const aol::MultiVector<RealType> &ArgDisplacement,
                                 const int K,
                                 const int L ) :
    qc::FELinDeformAsymMatrixWeightedStiffInterface<ConfiguratorType, HyperelasticDeformSubHessian<ConfiguratorType,HyperelasticEnergyDensityType> > ( Grid, Displacement ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacement( Grid, Displacement ),
    _argDisplacement( Grid, ArgDisplacement ),
    _k( K ),
    _l( L ){}

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &/*LocalCoord*/,
                               const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                               typename ConfiguratorType::MatType &Matrix ) const {
    // compute the deformation gradient $\nabla\Phi$
    typename ConfiguratorType::MatType dPhi;
    _displacement.evaluateGradientAtQuadPoint( El, QuadPoint, dPhi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dPhi[i][i] += 1.;

    // compute the deformation gradient $\nabla\phi$
    typename ConfiguratorType::MatType dphi;
    _argDisplacement.evaluateGradient( TransformedEl, TransformedLocalCoord, dphi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dphi[i][i] += 1.;

    // compute the Hessian matrix for the kth component multiplied from the left and the lth from the right
    dphi *= dPhi;
    HyperelasticSubHessian<ConfiguratorType,HyperelasticEnergyDensityType>::elasticityTensor( dphi, _hyperelasticEnergyDensity, El, QuadPoint, _k, _l, Matrix );
  }
};

/**
  * \brief This class computes via "apply(...)" or "applyAdd(...)" a hyperelastic energy Hessian with respect to the displacement,
  * where the material is prestressed by a deformation \f$ \Phi \f$,
  * i.e. \f$ (<\frac{\partial^2 \int_\Omega W(I_1,I_2,I_3,x) dx}{\partial \phi^2},\psi_i\circ\Phi,\psi_j\circ\Phi>)_{i,j} \f$,
  * where the \f$ \psi_i \f$ and \f$ \psi_j \f$ are the vector-valued FE basis functions,
  * the hyperelastic invariants are given by \f$ I_1=||\nabla\phi\nabla\Phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi\nabla\Phi)||_F^2 \f$,
  * \f$ I_3=det(\nabla\phi\nabla\Phi)\f$, and the function \f$ W(I_1,I_2,I_3,x) \f$ is passed to the constructor (e.g. HyerelasticEnergyDensityDefault).
  * In its first argument, "apply" and "applyAdd" expect a Multivector which represents the displacements \f$ d_i \f$
  * in the different Cartesian directions; the deformation \f$ \phi \f$ is then automatically computed via \f$ \phi=d+\f$ identity.
  *
  * \author Wirth
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType, typename SubMatrixType = typename ConfiguratorType::MatrixType>
class HyperelasticDeformHessian
      : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::BlockMatrix<SubMatrixType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the underlying FE grid
  const typename ConfiguratorType::InitType &_grid;
  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;
  // the prestressing displacement
  const aol::MultiVector<RealType> &_displacement;

public:
  HyperelasticDeformHessian ( const typename ConfiguratorType::InitType &Grid,
                              const HyperelasticEnergyDensityType &HyperelasticEnergyDensity,
                              const aol::MultiVector<RealType> &Displacement ) :
    _grid ( Grid ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacement( Displacement ){}

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::BlockMatrix<SubMatrixType> &Dest ) const {
    for ( int k = 0; k < ConfiguratorType::Dim; k++ )
      for ( int l = 0; l < ConfiguratorType::Dim; l++ )
        HyperelasticDeformSubHessian<ConfiguratorType,HyperelasticEnergyDensityType>( _grid, _hyperelasticEnergyDensity, _displacement, Arg, k, l ).assembleAddMatrix( Dest.getReference( k, l ) );
  }

  void applyAdd( const aol::MultiVector<RealType> &Arg, aol::BlockMatrix<SubMatrixType> &Dest, const aol::BitVector &Mask ) const {
    for ( int k = 0; k < ConfiguratorType::Dim; k++ )
      for ( int l = 0; l < ConfiguratorType::Dim; l++ )
        HyperelasticDeformSubHessian<ConfiguratorType,HyperelasticEnergyDensityType>( _grid, _hyperelasticEnergyDensity, _displacement, Arg, k, l ).assembleAddMatrix( Mask, Dest.getReference( k, l ), aol::NumberTrait<RealType>::one, k == l );
  }
};

/**
  * \brief This class computes via "apply(...)" or "applyAdd(...)" a hyperelastic energy \f$ \int_\Omega W(I_1,I_2,I_3,x)det(\nabla\phi_r) dx \f$,
  * where the hyperelastic invariants are given by \f$ I_1=||\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r||_F^2 \f$, \f$ I_2=||cof(\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r||_F^2 \f$,
  * \f$ I_1=det(\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r) \f$ and the function \f$ W(I_1,I_2,I_3,x) \f$ is passed to the constructor (e.g. HyerelasticEnergyDensityDefault).
  * \f$ \phi \f$ is a deformation passed to the constructor.
  * In its first argument, "apply" and "applyAdd" expect a Multivector which represents the displacements \f$ d_i \f$ of \f$ \phi_l \f$ and \f$ \phi_r \f$
  * in the different Cartesian directions.
  *
  * \author Wirth
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class HyperelasticPrePostDeformEnergy
      : public FENonlinDeformIntegrationVectorInterface < ConfiguratorType, HyperelasticPrePostDeformEnergy<ConfiguratorType,HyperelasticEnergyDensityType>, 2*ConfiguratorType::Dim > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;
  // the deformation \phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement;
public:
  HyperelasticPrePostDeformEnergy ( const typename ConfiguratorType::InitType &Grid,
                                    const HyperelasticEnergyDensityType &HyperelasticEnergyDensity,
                                    const aol::MultiVector<RealType> &Displacement ) :
    FENonlinDeformIntegrationVectorInterface < ConfiguratorType, HyperelasticPrePostDeformEnergy<ConfiguratorType,HyperelasticEnergyDensityType>, 2*ConfiguratorType::Dim > ( Grid, Displacement ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacement( Grid, Displacement ) {}

  RealType evaluateIntegrand ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,2*ConfiguratorType::Dim> &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &LocalCoord,
                               const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord ) const {
    // compute the deformation gradient $\nabla\phi$
    typename ConfiguratorType::MatType dPhi;
    _displacement.evaluateGradientAtQuadPoint( El, QuadPoint, dPhi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dPhi[i][i] += 1.;

    // compute the deformation gradient $\nabla\phi_l$ and $\nabla\phi_r$
    typename ConfiguratorType::MatType dPhiL, dPhiR;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ){
      DiscFuncs[i].evaluateGradient( TransformedEl, TransformedLocalCoord, dPhiL[i] );
      DiscFuncs[ConfiguratorType::Dim+i].evaluateGradient( El, LocalCoord, dPhiR[i] );
      dPhiL[i][i] += 1.;
      dPhiR[i][i] += 1.;
    }

    dPhiL *= dPhi;
    dPhiL *= dPhiR.inverse();
    return aol::Abs( dPhiR.det() ) * HyperelasticEnergy<ConfiguratorType,HyperelasticEnergyDensityType>::energyDensity( dPhiL, _hyperelasticEnergyDensity, El, QuadPoint );
  }
};

/**
  * \brief This class computes via "apply(...)" or "applyAdd(...)" a hyperelastic energy gradient with respect to one of two displacements,
  * i.e. \f$ (<\frac{\partial \int_\Omega W(I_1,I_2,I_3,x)det(\nabla\phi_r) dx}{\partial \phi_l},\psi_i>)_i \f$,
  * where the \f$ \psi_i \f$ are the vector-valued FE basis functions,
  * the hyperelastic invariants are given by \f$ I_1=||\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r||_F^2 \f$, \f$ I_2=||cof(\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r)||_F^2 \f$,
  * \f$ I_1=det(\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r) \f$ and the function \f$ W(I_1,I_2,I_3,x) \f$ is passed to the constructor (e.g. HyerelasticEnergyDensityDefault).
  * \f$ \phi \f$ is some prestressing deformation.
  * In its first argument, "apply" and "applyAdd" expect a Multivector which represents the displacements \f$ d_i \f$ of \f$ \phi_l \f$ and \f$ \phi_r \f$
  * in the different Cartesian directions.
  *
  * \author Wirth
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class HyperelasticPrePostDeformGradientPostComp
      : public FENonlinDeformVectorDiffOpInterface<ConfiguratorType, 2*ConfiguratorType::Dim, ConfiguratorType::Dim, HyperelasticPrePostDeformGradientPostComp<ConfiguratorType,HyperelasticEnergyDensityType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;
  // the deformation \phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement;
public:
  HyperelasticPrePostDeformGradientPostComp ( const typename ConfiguratorType::InitType &Grid,
                                              const HyperelasticEnergyDensityType &HyperelasticEnergyDensity,
                                              const aol::MultiVector<RealType> &Displacement ) :
    FENonlinDeformVectorDiffOpInterface<ConfiguratorType, 2*ConfiguratorType::Dim, ConfiguratorType::Dim, HyperelasticPrePostDeformGradientPostComp<ConfiguratorType,HyperelasticEnergyDensityType> > ( Grid, Displacement ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacement( Grid, Displacement ) {}

  void getNonlinearity ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,2*ConfiguratorType::Dim> &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &LocalCoord,
                         const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                         aol::Mat<ConfiguratorType::Dim, ConfiguratorType::Dim, RealType> &NL ) const {
    // compute the deformation gradient $\nabla\phi$
    typename ConfiguratorType::MatType dPhi;
    _displacement.evaluateGradientAtQuadPoint( El, QuadPoint, dPhi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dPhi[i][i] += 1.;

    // compute the deformation gradient $\nabla\phi_l$ and $\nabla\phi_r$
    typename ConfiguratorType::MatType dPhiL, dPhiR;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ){
      DiscFuncs[i].evaluateGradient( TransformedEl, TransformedLocalCoord, dPhiL[i] );
      DiscFuncs[ConfiguratorType::Dim+i].evaluateGradient( El, LocalCoord, dPhiR[i] );
      dPhiL[i][i] += 1.;
      dPhiR[i][i] += 1.;
    }

    // compute \partial_{,A}W(I_1,I_2,I_3)
    dPhiL *= dPhi;
    dPhiL *= dPhiR.inverse();
    HyperelasticGradient<ConfiguratorType,HyperelasticEnergyDensityType>::firstPiolaKirchhoffStress( dPhiL, _hyperelasticEnergyDensity, El, QuadPoint, NL );

    // compute NL = det(\nabla\phi_r) \partial_{,A}W(I_1,I_2,I_3) \nabla\phi_r^{-T}
    typename ConfiguratorType::MatType cofDPhiR;
    cofDPhiR.makeCofactorMatrix( dPhiR );
    NL *= cofDPhiR;
  }
};

/**
  * \brief This class computes via "apply(...)" or "applyAdd(...)" a hyperelastic energy gradient with respect to one of two displacements,
  * i.e. \f$ (<\frac{\partial \int_\Omega W(I_1,I_2,I_3,x)det(\nabla\phi_r) dx}{\partial \phi_r},\psi_i>)_i \f$,
  * where the \f$ \psi_i \f$ are the vector-valued FE basis functions,
  * the hyperelastic invariants are given by \f$ I_1=||\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r||_F^2 \f$, \f$ I_2=||cof(\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r)||_F^2 \f$,
  * \f$ I_1=det(\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r) \f$ and the function \f$ W(I_1,I_2,I_3,x) \f$ is passed to the constructor (e.g. HyerelasticEnergyDensityDefault).
  * \f$ \phi \f$ is some prestressing deformation.
  * In its first argument, "apply" and "applyAdd" expect a Multivector which represents the displacements \f$ d_i \f$ of \f$ \phi_l \f$ and \f$ \phi_r \f$
  * in the different Cartesian directions.
  *
  * \author Wirth
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class HyperelasticPrePostDeformGradientPreComp
      : public aol::FENonlinVectorDiffOpInterface<ConfiguratorType, 2*ConfiguratorType::Dim, ConfiguratorType::Dim, HyperelasticPrePostDeformGradientPreComp<ConfiguratorType,HyperelasticEnergyDensityType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;
  // the deformation \phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement;
public:
  HyperelasticPrePostDeformGradientPreComp ( const typename ConfiguratorType::InitType &Grid,
                                             const HyperelasticEnergyDensityType &HyperelasticEnergyDensity,
                                             const aol::MultiVector<RealType> &Displacement ) :
    aol::FENonlinVectorDiffOpInterface<ConfiguratorType, 2*ConfiguratorType::Dim, ConfiguratorType::Dim, HyperelasticPrePostDeformGradientPreComp<ConfiguratorType,HyperelasticEnergyDensityType> > ( Grid ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacement( Grid, Displacement ) {}

  void getNonlinearity ( const aol::auto_container<2*ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &LocalCoord,
                         aol::Mat<ConfiguratorType::Dim, ConfiguratorType::Dim, RealType> &NL ) const {
    // compute the deformation gradient $\nabla\phi$
    typename ConfiguratorType::MatType dPhi;
    _displacement.evaluateGradientAtQuadPoint( El, QuadPoint, dPhi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dPhi[i][i] += 1.;

    // compute $\phi(x)$
    typename ConfiguratorType::ElementType transformedEl;
    typename ConfiguratorType::DomVecType  transformedCoord;
    qc::transformAndClipCoord<ConfiguratorType>( *this->_config, _displacement, El, QuadPoint, LocalCoord, transformedEl, transformedCoord );

    // compute the deformation gradient $\nabla\phi_l$ and $\nabla\phi_r$
    typename ConfiguratorType::MatType dPhiL, dPhiR;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ){
      DiscFuncs[i].evaluateGradient( transformedEl, transformedCoord, dPhiL[i] );
      DiscFuncs[ConfiguratorType::Dim+i].evaluateGradient( El, LocalCoord, dPhiR[i] );
      dPhiL[i][i] += 1.;
      dPhiR[i][i] += 1.;
    }

    // compute W(I_1,I_2,I_3) and \partial_{,A}W(I_1,I_2,I_3)
    dPhiL *= dPhi;
    dPhiL *= dPhiR.inverse();
    const RealType density = HyperelasticEnergy<ConfiguratorType,HyperelasticEnergyDensityType>::energyDensity( dPhiL, _hyperelasticEnergyDensity, El, QuadPoint );
    HyperelasticGradient<ConfiguratorType,HyperelasticEnergyDensityType>::firstPiolaKirchhoffStress( dPhiL, _hyperelasticEnergyDensity, El, QuadPoint, NL );

    // compute NL = W(I_1,I_2,I_3) cof(\nabla\phi_r) + det(\nabla\phi_r) \partial_{,A}W(I_1,I_2,I_3)
    typename ConfiguratorType::MatType cofDPhiR, aux;
    cofDPhiR.makeCofactorMatrix( dPhiR );
    aux.makeProduct( NL, cofDPhiR );
    NL.makeProductAtransposedB( dPhiL, aux );
    NL *= -1.;
    cofDPhiR *= density;
    NL += cofDPhiR;
  }
};

/**
  * \brief This class computes via "apply(...)" or "applyAdd(...)" part of the hyperelastic energy Hessian with respect to one of two displacements,
  * i.e. \f$ (<\frac{\partial^2\int_\Omega W(I_1,I_2,I_3,x)det(\nabla\phi_r) dx}{\partial \phi_l^2},\psi_i,\psi_j>)_{ij} \f$,
  * where the \f$ \psi_i \f$ are the vector-valued FE basis functions,
  * the hyperelastic invariants are given by \f$ I_1=||\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r||_F^2 \f$, \f$ I_2=||cof(\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r)||_F^2 \f$,
  * \f$ I_1=det(\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r) \f$ and the function \f$ W(I_1,I_2,I_3,x) \f$ is passed to the constructor (e.g. HyerelasticEnergyDensityDefault).
  * \f$ \phi \f$ is some prestressing deformation.
  * In its first argument, "apply" and "applyAdd" expect a Multivector which represents the displacements \f$ d_i \f$ of \f$ \phi_l \f$ and \f$ \phi_r \f$
  * in the different Cartesian directions.
  *
  * \author Wirth
  * \ingroup MatrixFEOp
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class HyperelasticPrePostDeformSubHessianPostComp
      : public qc::FELinDeformAsymMatrixWeightedStiffInterface<ConfiguratorType, HyperelasticPrePostDeformSubHessianPostComp<ConfiguratorType,HyperelasticEnergyDensityType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;
  // the deformations \phi, \phi_r and \phi_l
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement, _preDisplacement, _postDisplacement;
  // the indices of the sub-Hessian
  const int _k, _l;
public:
  HyperelasticPrePostDeformSubHessianPostComp ( const typename ConfiguratorType::InitType &Grid,
                                                const HyperelasticEnergyDensityType &HyperelasticEnergyDensity,
                                                const aol::MultiVector<RealType> &Displacement,
                                                const aol::MultiVector<RealType> &PreDisplacement,
                                                const aol::MultiVector<RealType> &PostDisplacement,
                                                const int K,
                                                const int L ) :
    qc::FELinDeformAsymMatrixWeightedStiffInterface<ConfiguratorType, HyperelasticPrePostDeformSubHessianPostComp<ConfiguratorType,HyperelasticEnergyDensityType> > ( Grid, Displacement ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacement( Grid, Displacement ),
    _preDisplacement( Grid, PreDisplacement ),
    _postDisplacement( Grid, PostDisplacement ),
    _k( K ),
    _l( L ){}

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &LocalCoord,
                               const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                               typename ConfiguratorType::MatType &Matrix ) const {
    // compute the deformation gradient $\nabla\Phi$
    typename ConfiguratorType::MatType dPhi;
    _displacement.evaluateGradientAtQuadPoint( El, QuadPoint, dPhi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dPhi[i][i] += 1.;

    // compute the deformation gradient $\nabla\phi_l$ and $\nabla\phi_r$
    typename ConfiguratorType::MatType dPhiL, dPhiR;
    _preDisplacement.evaluateGradient( El, LocalCoord, dPhiR );
    _postDisplacement.evaluateGradient( TransformedEl, TransformedLocalCoord, dPhiL );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ){
      dPhiL[i][i] += 1.;
      dPhiR[i][i] += 1.;
    }

    // compute the Hessian matrix for the kth component multiplied from the left and the lth from the right
    typename ConfiguratorType::MatType dPhiRInv, cofDPhiR, aux;
    cofDPhiR.makeCofactorMatrix( dPhiR );
    dPhiRInv = dPhiR.inverse();
    dPhiL *= dPhi;
    dPhiL *= dPhiRInv;
    HyperelasticSubHessian<ConfiguratorType,HyperelasticEnergyDensityType>::elasticityTensor( dPhiL, _hyperelasticEnergyDensity, El, QuadPoint, _k, _l, Matrix );
    aux.makeProduct( Matrix, cofDPhiR );
    Matrix.makeProduct( dPhiRInv, aux );
  }
};

/**
  * \brief This class computes via "apply(...)" or "applyAdd(...)" part of the hyperelastic energy Hessian with respect to one of two displacements,
  * i.e. \f$ (<\frac{\partial^2\int_\Omega W(I_1,I_2,I_3,x)det(\nabla\phi_r) dx}{\partial \phi_r^2},\psi_i,\psi_j>)_{ij} \f$,
  * where the \f$ \psi_i \f$ are the vector-valued FE basis functions,
  * the hyperelastic invariants are given by \f$ I_1=||\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r||_F^2 \f$, \f$ I_2=||cof(\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r)||_F^2 \f$,
  * \f$ I_1=det(\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r) \f$ and the function \f$ W(I_1,I_2,I_3,x) \f$ is passed to the constructor (e.g. HyerelasticEnergyDensityDefault).
  * \f$ \phi \f$ is some prestressing deformation.
  * In its first argument, "apply" and "applyAdd" expect a Multivector which represents the displacements \f$ d_i \f$ of \f$ \phi_l \f$ and \f$ \phi_r \f$
  * in the different Cartesian directions.
  *
  * \author Wirth
  * \ingroup MatrixFEOp
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class HyperelasticPrePostDeformSubHessianPreComp
      : public aol::FELinAsymMatrixWeightedStiffInterface<ConfiguratorType, HyperelasticPrePostDeformSubHessianPreComp<ConfiguratorType,HyperelasticEnergyDensityType> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;
  // the deformations \phi, \phi_r and \phi_l
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement, _preDisplacement, _postDisplacement;
  // the indices of the sub-Hessian
  const int _k, _l;
public:
  HyperelasticPrePostDeformSubHessianPreComp( const typename ConfiguratorType::InitType &Grid,
                                              const HyperelasticEnergyDensityType &HyperelasticEnergyDensity,
                                              const aol::MultiVector<RealType> &Displacement,
                                              const aol::MultiVector<RealType> &PreDisplacement,
                                              const aol::MultiVector<RealType> &PostDisplacement,
                                              const int K,
                                              const int L ) :
    aol::FELinAsymMatrixWeightedStiffInterface<ConfiguratorType, HyperelasticPrePostDeformSubHessianPreComp<ConfiguratorType,HyperelasticEnergyDensityType> > ( Grid ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacement( Grid, Displacement ),
    _preDisplacement( Grid, PreDisplacement ),
    _postDisplacement( Grid, PostDisplacement ),
    _k( K ),
    _l( L ) {}

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint,
                               const typename ConfiguratorType::DomVecType& DomRefCoord, typename ConfiguratorType::MatType &Matrix ) const {
    // compute the deformation gradient $\nabla\phi$
    typename ConfiguratorType::MatType dPhi;
    _displacement.evaluateGradientAtQuadPoint( El, QuadPoint, dPhi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dPhi[i][i] += 1.;

    // compute $\phi(x)$
    typename ConfiguratorType::ElementType transformedEl;
    typename ConfiguratorType::DomVecType  transformedCoord;
    qc::transformAndClipCoord<ConfiguratorType>( *this->_config, _displacement, El, QuadPoint, DomRefCoord, transformedEl, transformedCoord );

    // compute the deformation gradient $\nabla\phi_l$ and $\nabla\phi_r$
    typename ConfiguratorType::MatType dPhiL, dPhiR;
    _postDisplacement.evaluateGradient( transformedEl, transformedCoord, dPhiL );
    _preDisplacement.evaluateGradient( El, DomRefCoord, dPhiR );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ){
      dPhiL[i][i] += 1.;
      dPhiR[i][i] += 1.;
    }

    // compute D\phi_l D\phi D\phi_r^{-1} and \partial_{,A}W
    typename ConfiguratorType::MatType defGrad( dPhiL ), cofDPhiR, invDPhiR( dPhiR.inverse() ), dPhiInvDPhiR( dPhi ), dW;
    cofDPhiR.makeCofactorMatrix( dPhiR );
    dPhiInvDPhiR *= invDPhiR;
    defGrad *= dPhiInvDPhiR;
    HyperelasticGradient<ConfiguratorType,HyperelasticEnergyDensityType>::firstPiolaKirchhoffStress( defGrad, _hyperelasticEnergyDensity, El, QuadPoint, dW );

    // compute the Hessian matrix; denote variations of \phi_r by \psi and \chi; they only act on the component "_k" and "_l"
    // \psi shall belong to vectors multiplied from the left, \chi from the right
    // use tr(A D\psi B D\chi) = ( D\psi_{k:} B_{:l} )( D\chi_{l:} A_{:k} ) = D\psi_{k:} ( B_{:l} A_{:k}^T ) D\chi_{l:}
    typename ConfiguratorType::MatType aux1, aux2;
    typename ConfiguratorType::VecType vec1;
    // compute M s.t. (D\psi)_{k}M(D\chi)_{l}^T = [-\partial_{,A}W:(D\phi_l D\phi D\phi_r^{-1} D\psi D\phi_r^{-1})][cof(D\phi_r):D\chi]
    aux1.makeProductAtransposedB( defGrad, dW );
    aux2.makeProductABtransposed( aux1, invDPhiR );
    Matrix.makeTensorProduct( aux2[_k], cofDPhiR[_l] );
    Matrix *= -1.;
    // compute M s.t. (D\psi)_{k}M(D\chi)_{l}^T = [-\partial_{,A}W:(D\phi_l D\phi D\phi_r^{-1} D\chi D\phi_r^{-1})][cof(D\phi_r):D\psi]
    aux1.makeTensorProduct( cofDPhiR[_k], aux2[_l] );
    Matrix -= aux1;
    // compute M s.t. (D\psi)_{k}M(D\chi)_{l}^T = W<\partial_{,AA}det D\phi_r,D\psi,D\chi>
    //   = [(cof D\phi_r;D\psi)(cof D\phi_r;D\chi)-tr(cof D\phi_r D\psi^T cof D\phi_r D\chi^T)] / det D\phi_r
    aux1.makeTensorProduct( cofDPhiR[_k], cofDPhiR[_l] );
    aux2.makeTensorProduct( cofDPhiR[_l], cofDPhiR[_k] );
    aux1 -= aux2;
    Matrix.addMultiple( aux1, HyperelasticEnergy<ConfiguratorType,HyperelasticEnergyDensityType>::energyDensity( defGrad, _hyperelasticEnergyDensity, El, QuadPoint ) / dPhiR.det() );
    // compute M s.t. (D\psi)_{k}M(D\chi)_{l}^T = det D\phi_r<\partial_{,AA}W,D\phi_l D\phi D\phi_r^{-1} D\psi D\phi_r^{-1},D\phi_l D\phi D\phi_r^{-1} D\chi D\phi_r^{-1}>
    aux1.setZero();
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        HyperelasticSubHessian<ConfiguratorType,HyperelasticEnergyDensityType>::elasticityTensor( defGrad, _hyperelasticEnergyDensity, El, QuadPoint, i, j, aux2 );
        aux1.addMultiple( aux2, defGrad[j][_l] * defGrad[i][_k] );
      }
    aux1 *= cofDPhiR;
    aux2.makeProduct( invDPhiR, aux1 );
    Matrix += aux2;
    // compute M s.t. (D\psi)_{k}M(D\chi)_{l}^T = det D\phi_r[\partial_{,A}W:(D\phi_l D\phi D\phi_r^{-1} [D\psi D\phi_r^{-1} D\chi + D\chi D\phi_r^{-1} D\psi] D\phi_r^{-1})]
    aux1.makeProductAtransposedB( defGrad, dW );
    aux1 *= cofDPhiR;
    invDPhiR.getColumn( _l, vec1 );
    aux2.makeTensorProduct( vec1, aux1[_k] );
    Matrix += aux2;
    invDPhiR.getColumn( _k, vec1 );
    aux2.makeTensorProduct( aux1[_l], vec1 );
    Matrix += aux2;
  }
};

/**
  * \brief This class computes via "apply(...)" or "applyAdd(...)" part of the hyperelastic energy Hessian with respect to two displacements,
  * i.e. \f$ (<\frac{\partial^2\int_\Omega W(I_1,I_2,I_3,x)det(\nabla\phi_r) dx}{\partial\phi_l\partial\phi_r},\psi_i,\psi_j>)_{ij} \f$,
  * where the \f$ \psi_i \f$ are the vector-valued FE basis functions,
  * the hyperelastic invariants are given by \f$ I_1=||\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r||_F^2 \f$, \f$ I_2=||cof(\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r)||_F^2 \f$,
  * \f$ I_1=det(\nabla(\phi_l\circ\phi\circ\phi_r^{-1})\circ\phi_r) \f$ and the function \f$ W(I_1,I_2,I_3,x) \f$ is passed to the constructor (e.g. HyerelasticEnergyDensityDefault).
  * \f$ \phi \f$ is some prestressing deformation.
  * In its first argument, "apply" and "applyAdd" expect a Multivector which represents the displacements \f$ d_i \f$ of \f$ \phi_l \f$ and \f$ \phi_r \f$
  * in the different Cartesian directions.
  * If ``PhiRVariationFromRight''=true, then the resulting matrix is the Hessian with respect to a variation of the ``L''th component
  * of \f$ \phi_r \f$ from the right and a vatiation of the ``K''th component of \f$ \phi_l \f$ from the left;
  * else \f$ \phi_l \f$ and \f$ \phi_r \f$ are exchanged.
  *
  * \author Wirth
  * \ingroup MatrixFEOp
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class HyperelasticPrePostDeformSubHessianMixedPart
      : public qc::FELinDeformAsymMatrixWeightedStiffInterface<ConfiguratorType, HyperelasticPrePostDeformSubHessianMixedPart<ConfiguratorType,HyperelasticEnergyDensityType>, false> {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::MatType    MatType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;
  // the deformations \phi, \phi_r and \phi_l
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement, _preDisplacement, _postDisplacement;
  // the indices of the sub-Hessian
  const int _k, _l;
  // flag indicating whether the test direction for phi_r is multiplied from the right or the left
  const bool _phiRVariationFromRight;
public:
  HyperelasticPrePostDeformSubHessianMixedPart( const typename ConfiguratorType::InitType &Grid,
                                                const HyperelasticEnergyDensityType &HyperelasticEnergyDensity,
                                                const aol::MultiVector<RealType> &Displacement,
                                                const aol::MultiVector<RealType> &PreDisplacement,
                                                const aol::MultiVector<RealType> &PostDisplacement,
                                                const int K,
                                                const int L,
                                                const bool PhiRVariationFromRight ) :
    qc::FELinDeformAsymMatrixWeightedStiffInterface<ConfiguratorType, HyperelasticPrePostDeformSubHessianMixedPart<ConfiguratorType,HyperelasticEnergyDensityType>, false> ( Grid, Displacement, aol::ASSEMBLED, PhiRVariationFromRight ? qc::LEFT : qc::RIGHT ),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacement( Grid, Displacement ),
    _preDisplacement( Grid, PreDisplacement ),
    _postDisplacement( Grid, PostDisplacement ),
    _k( K ),
    _l( L ),
    _phiRVariationFromRight( PhiRVariationFromRight ){}

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const DomVecType &LocalCoord,
                               const typename ConfiguratorType::ElementType &TransformedEl, const DomVecType &TransformedLocalCoord,
                               MatType &Matrix ) const {
    // compute the deformation gradient $\nabla\phi$
    typename ConfiguratorType::MatType dPhi;
    _displacement.evaluateGradientAtQuadPoint( El, QuadPoint, dPhi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dPhi[i][i] += 1.;

    // compute the deformation gradient $\nabla\phi_l$ and $\nabla\phi_r$
    typename ConfiguratorType::MatType dPhiL, dPhiR;
    _postDisplacement.evaluateGradient( TransformedEl, TransformedLocalCoord, dPhiL );
    _preDisplacement.evaluateGradient( El, LocalCoord, dPhiR );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ){
      dPhiL[i][i] += 1.;
      dPhiR[i][i] += 1.;
    }

    // compute D\phi_l D\phi D\phi_r^{-1} and \partial_{,A}W
    typename ConfiguratorType::MatType defGrad( dPhiL ), cofDPhiR, invDPhiR( dPhiR.inverse() ), dPhiInvDPhiR( dPhi ), dW;
    cofDPhiR.makeCofactorMatrix( dPhiR );
    dPhiInvDPhiR *= invDPhiR;
    defGrad *= dPhiInvDPhiR;
    HyperelasticGradient<ConfiguratorType,HyperelasticEnergyDensityType>::firstPiolaKirchhoffStress( defGrad, _hyperelasticEnergyDensity, El, QuadPoint, dW );

    // compute the Hessian matrix, pretending the gradient of the variation of \phi_r is multiplied from the right
    // denote variations of \phi_r and \phi_l by \psi_r and \psi_l; these variations only act on the component "compPhiR" and "compPhiL"
    int compPhiR = _phiRVariationFromRight ? _l : _k,
        compPhiL = _phiRVariationFromRight ? _k : _l;
    typename ConfiguratorType::MatType aux1, aux2;
    typename ConfiguratorType::VecType vec1;
    // compute M s.t. (D\psi_l)_{"compPhiL"}M(D\psi_r)_{"compPhiR"}^T  = [\partial_{,A}W:(D\psi_l D\phi D\phi_r^{-1})][cof(D\phi_r):D\psi_r]
    aux1.makeProductABtransposed( dW, dPhiInvDPhiR );
    Matrix.makeTensorProduct( aux1[compPhiL], cofDPhiR[compPhiR] );
    // compute M s.t. (D\psi_l)_{"compPhiL"}M(D\psi_r)_{"compPhiR"}^T  = -det(D\phi_r)[\partial_{,A}W:(D\psi_l D\phi D\phi_r^{-1} D\psi_r D\phi_r^{-1})]
    aux1.makeProduct( dW, cofDPhiR );
    dPhiInvDPhiR.getColumn( compPhiR, vec1 );
    aux2.makeTensorProduct( vec1, aux1[compPhiL] );
    Matrix -= aux2;
    // compute M s.t. (D\psi_l)_{"compPhiL"}M(D\psi_r)_{"compPhiR"}^T  = -det(D\phi_r)[\partial_{,AA}W(D\psi_l D\phi D\phi_r^{-1})]:(D\phi_l D\phi D\phi_r^{-1} D\psi_r D\phi_r^{-1})
    aux1.setZero();
    for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
      HyperelasticSubHessian<ConfiguratorType,HyperelasticEnergyDensityType>::elasticityTensor( defGrad, _hyperelasticEnergyDensity, El, QuadPoint, compPhiL, j, aux2 );
      aux1.addMultiple( aux2, defGrad[j][compPhiR] );
    }
    aux2.makeProduct( dPhiInvDPhiR, aux1 );
    aux2 *= cofDPhiR;
    Matrix -= aux2;
    // adapt local matrix if gradient of the variation of \phi_r is multiplied from the left
    if ( !_phiRVariationFromRight )
      Matrix.transpose();
  }
};

/**
  * \brief This class computes via "apply(...)" or "applyAdd(...)" the full hyperelastic energy Hessian with respect to two displacements,
  * Note that due to the fact that there are deformed and undeformed test functions in qc::HyperelasticPrePostDeformSubHessianMixedPart<>,
  * we really need to ensure that each block entry of the Hessian is  aol::SparseMatrix<> (and not ConfiguratorType::MatrixType, which is a FastUniformGridMatrix).
  * 
  * \author Heeren
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class HyperelasticPrePostDeformHessianMixedPart 
     : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::BlockMatrix<aol::SparseMatrix<typename ConfiguratorType::RealType> >  > {
protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::MatType    MatType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::InitType   GridType;

  const GridType& _grid;
  // the hyperelastic energy density $W$ as function of the three deformation gradient invariants
  const HyperelasticEnergyDensityType &_hyperelasticEnergyDensity;
  // the deformations \phi, \phi_r and \phi_l
  const aol::MultiVector<RealType> _displacement, _preDisplacement, _postDisplacement;
  // flag indicating whether the test direction for phi_r is multiplied from the right or the left
  const bool _phiRVariationFromRight;
public:
  HyperelasticPrePostDeformHessianMixedPart( const GridType &Grid,
                                             const HyperelasticEnergyDensityType &HyperelasticEnergyDensity,
                                             const aol::MultiVector<RealType> &Displacement,
                                             const aol::MultiVector<RealType> &PreDisplacement,
                                             const aol::MultiVector<RealType> &PostDisplacement,
                                             const bool PhiRVariationFromRight ) :
    _grid(Grid),
    _hyperelasticEnergyDensity( HyperelasticEnergyDensity ),
    _displacement( Displacement ),
    _preDisplacement( PreDisplacement ),
    _postDisplacement( PostDisplacement ),
    _phiRVariationFromRight( PhiRVariationFromRight ){}

    void assembleAddMatrix( aol::BlockMatrix< aol::SparseMatrix<RealType> >& Hessian ) const {
      for ( int k = 0; k < ConfiguratorType::Dim; k++ )
        for ( int l = 0; l < ConfiguratorType::Dim; l++ )
          qc::HyperelasticPrePostDeformSubHessianMixedPart<ConfiguratorType,HyperelasticEnergyDensityType>( _grid, _hyperelasticEnergyDensity, _displacement, _preDisplacement, _postDisplacement, k, l, _phiRVariationFromRight ).assembleAddMatrix( Hessian.getReference( k, l ) );
    }

    void assembleMatrix( aol::BlockMatrix< aol::SparseMatrix<RealType> >& Hessian ) const {
      Hessian.setZero();
      assembleAddMatrix( Hessian );
    }
    
    void applyAdd( const aol::MultiVector<RealType>& /*Arg*/, aol::BlockMatrix< aol::SparseMatrix<RealType> >& Hessian ) const {
      assembleAddMatrix( Hessian );
    }

};
      
/**
  * \brief This class represents a standard implementation of the hyperelastic energy density \f$ W(I_1,I_2,I_3) \f$
  * with respect to the invariants \f$ I_1=||\nabla\phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi)||_F^2 \f$, \f$ I_3=det(\nabla\phi) \f$
  * of the deformation gradient \f$\nabla\phi\f$. In addition to \f$I_1,I_2,I_3\f$, the methods
  * "evaluate(...)" and "evaluateDerivative(...)" have to be passed a spatial position, since
  * in principle the energy density may be space dependent.
  *
  * \author Wirth
  */
template <typename ConfiguratorType>
class HyperelasticEnergyDensityDefault{
protected:
  typedef typename ConfiguratorType::RealType RealType;
  // weights of the different energy components
  const RealType _weightLengthEnergy;
  const RealType _weightSurfaceEnergy;
  const RealType _weightVolumeEnergy;
  // Maximum compression allowed. If this compression is reached, energy is set to a high value.
  // If the threshold I3t is reached, a gradient descent step with minimum step length should still be Armijo-rule-
  // compatible, i.e. in this case ln(I3)-ln(I3t)>[1/I3t*(I3-I3t)]/2 with I3=I3t+\tau_{min}*1/I3t, hence I3t\circ1e-5.
  const RealType _compressionThreshold;
  const RealType _valueAtThreshold;
public:
  HyperelasticEnergyDensityDefault( const RealType WeightLengthEnergy, const RealType WeightSurfaceEnergy, const RealType WeightVolumeEnergy, const RealType ValueAtThreshold = aol::NumberTrait<RealType>::Inf ) :
    _weightLengthEnergy ( WeightLengthEnergy ),
    _weightSurfaceEnergy ( WeightSurfaceEnergy ),
    _weightVolumeEnergy ( WeightVolumeEnergy ),
    _compressionThreshold ( 1.e-10 ),
    _valueAtThreshold ( ValueAtThreshold ) {}

  //! returns \f$\frac\mu2I_1+\frac\lambda4I_3^2-(\mu+\frac\lambda2)\log I_3-d\frac\mu2-\frac\lambda4\f$ for parameters \f$\mu,\lambda\f$, dimension \f$d\f$
  inline RealType evaluate ( const RealType I1, const RealType /*I2*/, const RealType I3,
                             const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return _weightLengthEnergy * I1 / 2
           + _weightVolumeEnergy * aol::Sqr( I3 ) / 4
           - ( _weightLengthEnergy + _weightVolumeEnergy / 2 ) * ( I3 <= _compressionThreshold ? log( _compressionThreshold ) + ( I3 - _compressionThreshold ) / _compressionThreshold  -_valueAtThreshold : log( I3 ) )
           - ( ConfiguratorType::Dim / 2. ) * _weightLengthEnergy - _weightVolumeEnergy / 4;
  }

  //! returns \f$\frac\mu2I_1+\frac\lambda4I_3^2-(\mu+\frac\lambda2)\log I_3-d\frac\mu2-\frac\lambda4\f$ for parameters \f$\mu,\lambda\f$, dimension \f$d\f$
  inline RealType evaluateAtQuadPoint ( const RealType I1, const RealType /*I2*/, const RealType I3,
                                        const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    return _weightLengthEnergy * I1 / 2
           + _weightVolumeEnergy * aol::Sqr( I3 ) / 4
           - ( _weightLengthEnergy + _weightVolumeEnergy / 2 ) * ( I3 <= _compressionThreshold ? log( _compressionThreshold ) + ( I3 - _compressionThreshold ) / _compressionThreshold  -_valueAtThreshold : log( I3 ) )
           - ( ConfiguratorType::Dim / 2. ) * _weightLengthEnergy - _weightVolumeEnergy / 4;
  }

  //! returns \f$(\frac\mu2,0,\frac\lambda2I_3-(\mu+\frac\lambda2)\frac1{I_3})\f$
  inline aol::Vec3<RealType> evaluateDerivative ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                                                  const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return aol::Vec3<RealType>( _weightLengthEnergy / 2,
                                0,
                                _weightVolumeEnergy * I3 / 2 - ( _weightLengthEnergy + _weightVolumeEnergy / 2 ) / aol::Max( I3, _compressionThreshold ) );
  }

  //! returns \f$(\frac\mu2,0,\frac\lambda2I_3-(\mu+\frac\lambda2)\frac1{I_3})\f$
  inline aol::Vec3<RealType> evaluateDerivativeAtQuadPoint ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                                                             const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    return aol::Vec3<RealType>( _weightLengthEnergy / 2,
                                0,
                                _weightVolumeEnergy * I3 / 2 - ( _weightLengthEnergy + _weightVolumeEnergy / 2 ) / aol::Max( I3, _compressionThreshold ) );
  }

  //! returns diag\f$(0,0,\frac\lambda2+(\mu+\frac\lambda2)\frac1{I_3^2})\f$
  inline aol::Matrix33<RealType> evaluateSecondDerivative ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                                                            const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return aol::Matrix33<RealType>( 0, 0, 0, 0, 0, 0, 0, 0,
                                    _weightVolumeEnergy / 2 + ( I3 <= _compressionThreshold ? 0. : ( _weightLengthEnergy + _weightVolumeEnergy / 2 ) / aol::Sqr( I3 ) ) );
  }

  //! returns diag\f$(0,0,\frac\lambda2+(\mu+\frac\lambda2)\frac1{I_3^2})\f$
  inline aol::Matrix33<RealType> evaluateSecondDerivativeAtQuadPoint ( const RealType /*I1*/, const RealType /*I2*/, const RealType I3,
                                                                       const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    return aol::Matrix33<RealType>( 0, 0, 0, 0, 0, 0, 0, 0,
                                    _weightVolumeEnergy / 2 + ( I3 <= _compressionThreshold ? 0. : ( _weightLengthEnergy + _weightVolumeEnergy / 2 ) / aol::Sqr( I3 ) ) );
  }
};

/**
  * \brief Implements the St Venant-Kirchhoff energy density \f$ W(I_1,I_2,I_3)=\frac\lambda2(I_1-d)^2+\frac\mu4((I_1-1)^2+d-1-2\zeta)+\nu\Gamma(I_3) \f$
  * with \f$ I_1=||\nabla\phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi)||_F^2 \f$, \f$ I_3=det(\nabla\phi) \f$, \f$ \zeta=I_3^2 \f$ in 2D and \f$ \zeta=I_2 \f$ in 3D,
  * \f$ \Gamma(I_3)=\frac1{I_3^2}+I_3^2-2 \f$, and the deformation gradient \f$\nabla\phi\f$, where \f$ \Gamma \f$ in fact does not belong to the law and
  * only prohibits complete compression.
  *
  * \author Wirth
  */
template <typename ConfiguratorType>
class HyperelasticEnergyDensityStVenantKirchhoff{
protected:
  typedef typename ConfiguratorType::RealType RealType;
  // Lame constants and weight of term preventing complete volume compression (which in fact does not belong to St Venant-Kirchhoff)
  const RealType _lambdaQuarter, _lambdaEighth;
  const RealType _mu, _muQuarter;
  const RealType _weightAntiCompression;
  const RealType _dimMinusOne;
public:
  HyperelasticEnergyDensityStVenantKirchhoff( const RealType Lambda, const RealType Mu, const RealType WeightAntiCompression ) :
    _lambdaQuarter( Lambda / 4 ),
    _lambdaEighth( Lambda / 8 ),
    _mu( Mu ),
    _muQuarter( Mu / 4 ),
    _weightAntiCompression( WeightAntiCompression ),
    _dimMinusOne( ConfiguratorType::Dim - 1 ) {}

  inline RealType evaluate ( const RealType I1, const RealType I2, const RealType I3,
                             const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return _lambdaEighth * aol::Sqr( I1 - ConfiguratorType::Dim )
           + _muQuarter * ( aol::Sqr( I1 - 1 ) + _dimMinusOne - 2 * ( ConfiguratorType::Dim == qc::QC_2D ? aol::Sqr( I3 ) : I2 ) )
           + _weightAntiCompression * ( 1 / aol::Sqr( I3 ) + aol::Sqr( I3 ) - 2 );
  }

  inline RealType evaluateAtQuadPoint ( const RealType I1, const RealType I2, const RealType I3,
                                        const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    return _lambdaEighth * aol::Sqr( I1 - ConfiguratorType::Dim )
           + _muQuarter * ( aol::Sqr( I1 - 1 ) + _dimMinusOne - 2 * ( ConfiguratorType::Dim == qc::QC_2D ? aol::Sqr( I3 ) : I2 ) )
           + _weightAntiCompression * ( 1 / aol::Sqr( I3 ) + aol::Sqr( I3 ) - 2 );
  }

  inline aol::Vec3<RealType> evaluateDerivative ( const RealType I1, const RealType /*I2*/, const RealType I3,
                                                  const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return aol::Vec3<RealType>( _lambdaQuarter * ( I1 - ConfiguratorType::Dim ) + _mu / 2 * ( I1 - 1 ),
                                ( ConfiguratorType::Dim == qc::QC_2D ? 0 : - _mu / 2 ),
                                ( ConfiguratorType::Dim == qc::QC_2D ? - _mu * I3 : 0 ) + _weightAntiCompression * 2 * ( I3 - 1 / aol::Cub( I3 ) ) );
  }

  inline aol::Vec3<RealType> evaluateDerivativeAtQuadPoint ( const RealType I1, const RealType /*I2*/, const RealType I3,
                                                             const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    return aol::Vec3<RealType>( _lambdaQuarter * ( I1 - ConfiguratorType::Dim ) + _mu / 2 * ( I1 - 1 ),
                                ( ConfiguratorType::Dim == qc::QC_2D ? 0 : - _mu / 2 ),
                                ( ConfiguratorType::Dim == qc::QC_2D ? - _mu * I3 : 0 ) + _weightAntiCompression * 2 * ( I3 - 1 / aol::Cub( I3 ) ) );
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivative ( const RealType /*I1*/, const RealType /*I2*/, const RealType /*I3*/,
                                                            const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return aol::Matrix33<RealType>( _lambdaQuarter + _mu / 2, 0, 0, 0, 0, 0, 0, 0, ( ConfiguratorType::Dim == qc::QC_2D ? - _mu : 0 ) );
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivativeAtQuadPoint ( const RealType /*I1*/, const RealType /*I2*/, const RealType /*I3*/,
                                                                       const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    return aol::Matrix33<RealType>( _lambdaQuarter + _mu / 2, 0, 0, 0, 0, 0, 0, 0, ( ConfiguratorType::Dim == qc::QC_2D ? - _mu : 0 ) );
  }
};

/**
  * \brief Implements the 2D linear energy density \f$ W(I_1,I_2,I_3)=\alpha(I_1+2-2\sqrt{I_1+2\det I_3})=\alpha dist(\nabla\phi,SO(2))^2 \f$
  * with \f$ I_1=||\nabla\phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi)||_F^2 \f$, \f$ I_3=det(\nabla\phi) \f$
  *
  * \author Wirth
  */
template <typename ConfiguratorType>
class HyperelasticEnergyDensity2DLinear{
protected:
  typedef typename ConfiguratorType::RealType RealType;
  const RealType _alpha;

public:
  HyperelasticEnergyDensity2DLinear( const RealType Alpha ) :
    _alpha( Alpha ) {
    if ( ConfiguratorType::Dim != qc::QC_2D )
      throw aol::Exception ( "qc::HyperelasticEnergyDensity2DLinear currently only implemented for 2D", __FILE__, __LINE__ );
  }

  inline RealType evaluate ( const RealType I1, const RealType /*I2*/, const RealType I3,
                             const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return _alpha * ( I1 + 2 - 2 * sqrt( I1 + 2 * I3 ) );
  }

  inline RealType evaluateAtQuadPoint ( const RealType I1, const RealType /*I2*/, const RealType I3,
                                        const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    return _alpha * ( I1 + 2 - 2 * sqrt( I1 + 2 * I3 ) );
  }

  inline aol::Vec3<RealType> evaluateDerivative ( const RealType I1, const RealType /*I2*/, const RealType I3,
                                                  const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    RealType auxTerm = - 1. / sqrt( I1 + 2 * I3 );
    return aol::Vec3<RealType>( _alpha * ( 1 + auxTerm ), 0, _alpha * 2 * auxTerm );
  }

  inline aol::Vec3<RealType> evaluateDerivativeAtQuadPoint ( const RealType I1, const RealType /*I2*/, const RealType I3,
                                                             const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    RealType auxTerm = - 1. / sqrt( I1 + 2 * I3 );
    return aol::Vec3<RealType>( _alpha * ( 1 + auxTerm ), 0, _alpha * 2 * auxTerm );
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivative ( const RealType I1, const RealType /*I2*/, const RealType I3,
                                                            const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    RealType auxTerm = aol::Cub( 1. / sqrt( I1 + 2 * I3 ) );
    return aol::Matrix33<RealType>( _alpha / 2. * auxTerm, 0, _alpha * auxTerm,
                                    0, 0, 0,
                                    _alpha * auxTerm, 0, _alpha * 2 * auxTerm );
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivativeAtQuadPoint ( const RealType I1, const RealType /*I2*/, const RealType I3,
                                                                       const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    RealType auxTerm = aol::Cub( 1. / sqrt( I1 + 2 * I3 ) );
    return aol::Matrix33<RealType>( _alpha / 2. * auxTerm, 0, _alpha * auxTerm,
                                    0, 0, 0,
                                    _alpha * auxTerm, 0, _alpha * 2 * auxTerm );
  }
};

/**
  * \brief Implements a generalized Mooney-Rivlin energy density \f$ W(I_1,I_2,I_3)=aI_1^{\frac{p}2}+bI_2^{\frac{q}2}+c\Gamma(I_3)-ad^{\frac{p}2}-bd^{\frac{q}2} \f$
  * with \f$ I_1=||\nabla\phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi)||_F^2 \f$, \f$ I_3=det(\nabla\phi) \f$,
  * \f$ \Gamma(I_3)=\frac1{I_3^s}+\frac{s}rI_3^r-1-\frac{s}r \f$, and the deformation gradient \f$\nabla\phi\f$.
  *
  * \author Wirth
  */
template <typename ConfiguratorType>
class HyperelasticEnergyDensityGeneralizedMooneyRivlin{
protected:
  typedef typename ConfiguratorType::RealType RealType;
  // material parameters
  const RealType _pHalf, _qHalf;
  const RealType _lengthFac, _areaFac, _volumeFac;
  const RealType _r, _s;
public:
  HyperelasticEnergyDensityGeneralizedMooneyRivlin( const RealType P, const RealType Q, const RealType LengthFac, const RealType AreaFac, const RealType VolumeFac, const RealType R, const RealType S ) :
    _pHalf( P / 2 ),
    _qHalf( Q / 2 ),
    _lengthFac( LengthFac ),
    _areaFac( AreaFac ),
    _volumeFac( VolumeFac ),
    _r( R ),
    _s( S ) {}

  inline RealType evaluate ( const RealType I1, const RealType I2, const RealType I3,
                             const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return _lengthFac * ( pow( I1, _pHalf ) - pow( ConfiguratorType::Dim, _pHalf ) )
           + _areaFac * ( pow( I2, _qHalf ) - pow( ConfiguratorType::Dim, _qHalf ) )
           + _volumeFac * ( pow( I3, -_s ) - 1 + ( _s / _r ) * ( pow( I3, _r ) - 1 ) );
  }

  inline RealType evaluateAtQuadPoint ( const RealType I1, const RealType I2, const RealType I3,
                                        const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    return _lengthFac * ( pow( I1, _pHalf ) - pow( ConfiguratorType::Dim, _pHalf ) )
           + _areaFac * ( pow( I2, _qHalf ) - pow( ConfiguratorType::Dim, _qHalf ) )
           + _volumeFac * ( pow( I3, -_s ) - 1 + ( _s / _r ) * ( pow( I3, _r ) - 1 ) );
  }

  inline aol::Vec3<RealType> evaluateDerivative ( const RealType I1, const RealType I2, const RealType I3,
                                                  const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return aol::Vec3<RealType>( _lengthFac * _pHalf * pow( I1, _pHalf - 1 ),
                                _areaFac * _qHalf * pow( I2, _qHalf - 1 ),
                                _volumeFac * _s * ( pow( I3, _r - 1 ) - pow( I3, -_s - 1 ) ) );
  }

  inline aol::Vec3<RealType> evaluateDerivativeAtQuadPoint ( const RealType I1, const RealType I2, const RealType I3,
                                                             const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    return aol::Vec3<RealType>( _lengthFac * _pHalf * pow( I1, _pHalf - 1 ),
                                _areaFac * _qHalf * pow( I2, _qHalf - 1 ),
                                _volumeFac * _s * ( pow( I3, _r - 1 ) - pow( I3, -_s - 1 ) ) );
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivative ( const RealType I1, const RealType I2, const RealType I3,
                                                            const typename ConfiguratorType::ElementType &/*El*/, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    return aol::Matrix33<RealType>( _lengthFac * _pHalf * ( _pHalf - 1 ) * pow( I1, _pHalf - 2 ), 0, 0,
                                    0, _areaFac * _qHalf * ( _qHalf - 1 ) * pow( I2, _qHalf - 2 ), 0,
                                    0, 0, _volumeFac * _s * ( ( _r - 1 ) * pow( I3, _r - 2 ) + ( _s + 1 ) * pow( I3, -_s - 2 ) ) );
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivativeAtQuadPoint ( const RealType I1, const RealType I2, const RealType I3,
                                                                       const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/ ) const {
    return aol::Matrix33<RealType>( _lengthFac * _pHalf * ( _pHalf - 1 ) * pow( I1, _pHalf - 2 ), 0, 0,
                                    0, _areaFac * _qHalf * ( _qHalf - 1 ) * pow( I2, _qHalf - 2 ), 0,
                                    0, 0, _volumeFac * _s * ( ( _r - 1 ) * pow( I3, _r - 2 ) + ( _s + 1 ) * pow( I3, -_s - 2 ) ) );
  }
};

/**
  * \brief This class represents a spatially weighted hyperelastic energy density \f$ w(x)W(I_1,I_2,I_3) \f$
  * with respect to the invariants \f$ I_1=||\nabla\phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi)||_F^2 \f$, \f$ I_3=det(\nabla\phi) \f$
  * of the deformation gradient \f$\nabla\phi\f$, where the hyperelastic energy \f$W\f$ and the weight \f$w\f$ have to be
  * passed to the constructor. \f$W\f$ can e.g. be "HyperelasticEnergyDensityDefault".
  *
  * \author Wirth
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class WeightedHyperelasticEnergyDensity {
private:
  typedef typename ConfiguratorType::RealType RealType;
  const HyperelasticEnergyDensityType &_energy;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _weight;
public:
  WeightedHyperelasticEnergyDensity( const HyperelasticEnergyDensityType &Energy, const typename ConfiguratorType::InitType &Grid, const aol::Vector<RealType> &Weight ) :
    _energy( Energy ),
    _weight( Grid, Weight ) {}

  inline RealType evaluate ( const RealType I1, const RealType I2, const RealType I3,
                             const typename ConfiguratorType::ElementType &El, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    return _weight.evaluate( El, RefCoord ) * _energy.evaluate( I1, I2, I3, El, RefCoord );
  }

  inline RealType evaluateAtQuadPoint ( const RealType I1, const RealType I2, const RealType I3,
                                        const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    return _weight.evaluateAtQuadPoint( El, QuadPoint ) * _energy.evaluateAtQuadPoint( I1, I2, I3, El, QuadPoint );
  }

  inline aol::Vec3<RealType> evaluateDerivative ( const RealType I1, const RealType I2, const RealType I3,
                                                  const typename ConfiguratorType::ElementType &El, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    return _weight.evaluate( El, RefCoord ) * _energy.evaluateDerivative( I1, I2, I3, El, RefCoord );
  }

  inline aol::Vec3<RealType> evaluateDerivativeAtQuadPoint ( const RealType I1, const RealType I2, const RealType I3,
                                                             const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    return _weight.evaluateAtQuadPoint( El, QuadPoint ) * _energy.evaluateDerivativeAtQuadPoint( I1, I2, I3, El, QuadPoint );
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivative ( const RealType I1, const RealType I2, const RealType I3,
                                                            const typename ConfiguratorType::ElementType &El, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    aol::Matrix33<RealType> result( _energy.evaluateSecondDerivative( I1, I2, I3, El, RefCoord ) );
    result *= _weight.evaluate( El, RefCoord );
    return result;
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivativeAtQuadPoint ( const RealType I1, const RealType I2, const RealType I3,
                                                                       const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    aol::Matrix33<RealType> result( _energy.evaluateSecondDerivativeAtQuadPoint( I1, I2, I3, El, QuadPoint ) );
    result *= _weight.evaluateAtQuadPoint( El, QuadPoint );
    return result;
  }
};

/**
  * \brief This class represents the hyperelastic energy density \f$ W(I_1,I_2,I_3) + f(I_3) \f$
  * with respect to the invariants \f$ I_1=||\nabla\phi||_F^2 \f$, \f$ I_2=||cof(\nabla\phi)||_F^2 \f$, \f$ I_3=det(\nabla\phi) \f$
  * of the deformation gradient \f$\nabla\phi\f$, where the hyperelastic energy \f$W\f$ has to be
  * passed to the constructor and \f$ f(x)=-\log(3x)+\frac5{12} \f$ for \f$ x<\frac13 \f$,
  * \f$ f(x)=-\frac92(\frac92x+1)(x-\frac23)^3 \f$ for \f$ \frac13< x <\frac23 \f$, and
  * \f$ f(x)=0 \f$ for \f$ x>\frac23 \f$.
  *
  * \author Wirth
  */
template <typename ConfiguratorType, typename HyperelasticEnergyDensityType>
class AntiCompressionHyperelasticEnergyDensity {
private:
  typedef typename ConfiguratorType::RealType RealType;
  const HyperelasticEnergyDensityType &_energy;
public:
  AntiCompressionHyperelasticEnergyDensity( const HyperelasticEnergyDensityType &Energy ) :
    _energy( Energy ) {}

  inline RealType evaluate ( const RealType I1, const RealType I2, const RealType I3,
                             const typename ConfiguratorType::ElementType &El, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    return _energy.evaluate( I1, I2, I3, El, RefCoord )
           + ( I3 < 0 ? aol::NumberTrait<RealType>::Inf
           : ( 3 * I3 < 1 ? -log( 3 * I3 ) + 5. / 12.
           : ( 3 * I3 < 2 ? -4.5 * ( 4.5 * I3 + 1 ) * aol::Cub( I3 - 2. / 3. ) : 0 ) ) );
  }

  inline RealType evaluateAtQuadPoint ( const RealType I1, const RealType I2, const RealType I3,
                                        const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    return _energy.evaluateAtQuadPoint( I1, I2, I3, El, QuadPoint )
           + ( I3 < 0 ? aol::NumberTrait<RealType>::Inf
           : ( 3 * I3 < 1 ? -log( 3 * I3 ) + 5. / 12.
           : ( 3 * I3 < 2 ? -4.5 * ( 4.5 * I3 + 1 ) * aol::Cub( I3 - 2. / 3. ) : 0 ) ) );
  }

  inline aol::Vec3<RealType> evaluateDerivative ( const RealType I1, const RealType I2, const RealType I3,
                                                  const typename ConfiguratorType::ElementType &El, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    aol::Vec3<RealType> result = _energy.evaluateDerivative( I1, I2, I3, El, RefCoord );
    result[2] += ( I3 < 0 ? -1
               : ( 3 * I3 < 1 ? -1. / I3
               : ( 3 * I3 < 2 ? -4.5 * ( 4.5 * aol::Cub( I3 - 2. / 3. ) + 3 * ( 4.5 * I3 + 1 ) * aol::Sqr( I3 - 2. / 3. ) ) : 0 ) ) );
    return result;
  }

  inline aol::Vec3<RealType> evaluateDerivativeAtQuadPoint ( const RealType I1, const RealType I2, const RealType I3,
                                                             const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    aol::Vec3<RealType> result = _energy.evaluateDerivativeAtQuadPoint( I1, I2, I3, El, QuadPoint );
    result[2] += ( I3 < 0 ? -1
               : ( 3 * I3 < 1 ? -1. / I3
               : ( 3 * I3 < 2 ? -4.5 * ( 4.5 * aol::Cub( I3 - 2. / 3. ) + 3 * ( 4.5 * I3 + 1 ) * aol::Sqr( I3 - 2. / 3. ) ) : 0 ) ) );
    return result;
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivative ( const RealType I1, const RealType I2, const RealType I3,
                                                            const typename ConfiguratorType::ElementType &El, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    aol::Matrix33<RealType> result( _energy.evaluateSecondDerivative( I1, I2, I3, El, RefCoord ) );
    result[2][2] += ( I3 < 0 ? 0
                  : ( 3 * I3 < 1 ? 1. / aol::Sqr( I3 )
                  : ( 3 * I3 < 2 ? -4.5 * ( 27 * aol::Sqr( I3 - 2. / 3. ) + 6 * ( 4.5 * I3 + 1 ) * ( I3 - 2. / 3. ) ) : 0 ) ) );
    return result;
  }

  inline aol::Matrix33<RealType> evaluateSecondDerivativeAtQuadPoint ( const RealType I1, const RealType I2, const RealType I3,
                                                                       const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    aol::Matrix33<RealType> result( _energy.evaluateSecondDerivativeAtQuadPoint( I1, I2, I3, El, QuadPoint ) );
    result[2][2] += ( I3 < 0 ? 0
                  : ( 3 * I3 < 1 ? 1. / aol::Sqr( I3 )
                  : ( 3 * I3 < 2 ? -4.5 * ( 27 * aol::Sqr( I3 - 2. / 3. ) + 6 * ( 4.5 * I3 + 1 ) * ( I3 - 2. / 3. ) ) : 0 ) ) );
    return result;
  }
};


template <typename RealType, typename VecType>
class RegularizedLinCombDescentDirOp : public aol::Op<VecType, VecType> { };

// specialization
template <class RealType>
class RegularizedLinCombDescentDirOp<RealType, aol::MultiVector<RealType> > : public aol::Op<aol::MultiVector<RealType> > {
protected:
  typedef aol::MultiVector<RealType> VecType;
  const qc::GridDefinition &_grid;
  RealType _smoothSigma;
  const aol::Op<VecType, VecType> &_massOpInv;
  aol::LinCombOp<VecType> _linCombOpReg;
  aol::LinCombOp<VecType> _linCombOpUnreg;

  mutable VecType *_gradient;

public:
  RegularizedLinCombDescentDirOp ( const qc::GridDefinition &grid,
                                   RealType smoothSigma,
                                   const aol::Op<VecType, VecType> &massOpInv )
      : _grid ( grid ), _smoothSigma ( smoothSigma ), _massOpInv ( massOpInv ), _gradient ( NULL ) {}

  virtual ~RegularizedLinCombDescentDirOp() {
    delete _gradient;
  }

  const VecType& getGradient() const {
    if ( !_gradient ) {
      throw aol::Exception ( "gradient not computed yet.", __FILE__, __LINE__ );
    }
    return *_gradient;
  }

  void appendReg ( const aol::Op<aol::MultiVector<RealType> > &op, RealType factor ) {
    if ( factor != 0. ) {
      _linCombOpReg.appendReference ( op, factor );
    }
  }

  void appendUnreg ( const aol::Op<aol::MultiVector<RealType> > &op, RealType factor ) {
    if ( factor != 0. ) {
      _linCombOpUnreg.appendReference ( op, factor );
    }
  }

  virtual void applyAdd ( const aol::MultiVector<RealType> &arg, aol::MultiVector<RealType> &dest ) const {
    VecType  tmp ( arg );
    if ( !_gradient ) {
      _gradient = new VecType ( arg );
    }
    VecType &gradient = *_gradient;

    // collect all gradients and reverse sign
    _linCombOpReg.apply ( arg, tmp );

    cerr << sqrt ( tmp.normSqr() ) << endl;
    // multiply with inverse lumped mass matrix
    _massOpInv.apply ( tmp, gradient );
    cerr << sqrt ( gradient.normSqr() ) << "   max x " << gradient[0].getMaxValue() << "   max y " << gradient[1].getMaxValue() << endl;


    qc::LinearSmoothOp<RealType> _linSmooth;
    _linSmooth.setCurrentGrid ( _grid );
    _linSmooth.setSigma ( _smoothSigma );
    for ( int c = 0; c < dest.numComponents(); c++ ) {
      _linSmooth.applyAdd ( gradient[c], dest[c] );
    }

    _linCombOpUnreg.apply ( arg, tmp );
    _massOpInv.applyAdd ( tmp, gradient );
    _massOpInv.applyAdd ( tmp, dest );

    dest *= -1.;

  }

};

//! Computes \f$ \frac{1}{2}\int |Du+Du^T|^2 \f$, where marg=u.
template <typename ConfiguratorType>
class SymmetricLengthEnergy : public aol::FENonlinIntegrationVectorInterface < ConfiguratorType,
                                                                               SymmetricLengthEnergy<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;

  SymmetricLengthEnergy ( const typename ConfiguratorType::InitType &Grid )
    : aol::FENonlinIntegrationVectorInterface < ConfiguratorType,
                                                SymmetricLengthEnergy<ConfiguratorType> > ( Grid ) {};

  RealType evaluateIntegrand ( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    typename ConfiguratorType::VecType grad;
    typename ConfiguratorType::MatType matrix;

    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      DiscrFuncs[i].evaluateGradientAtQuadPoint ( El, QuadPoint, grad );
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        matrix[i][j] += grad[j];
        matrix[j][i] += grad[j];
      }
    }

    return 0.5 * matrix.normSqr();
  }
};

//! Computes the first variation of SymmetricLengthEnergy, e.q. \f$ 2 int (Du+Du^T):D\zeta \f$, where marg=u.
template <typename ConfiguratorType>
class VariationOfSymmetricLengthEnergy
  : public aol::FENonlinVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, VariationOfSymmetricLengthEnergy<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;

  VariationOfSymmetricLengthEnergy ( const typename ConfiguratorType::InitType &Grid )
    : aol::FENonlinVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, VariationOfSymmetricLengthEnergy<ConfiguratorType> > ( Grid ) {};

  void getNonlinearity ( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/,
                         typename aol::Mat<ConfiguratorType::Dim, ConfiguratorType::Dim, RealType> &NL ) const {
    typename ConfiguratorType::VecType grad;
    NL.setZero();
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      DiscFuncs[i].evaluateGradientAtQuadPoint ( El, QuadPoint, grad );
      for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
        NL[i][j] += grad[j];
        NL[j][i] += grad[j];
      }
    }
    NL *= 2;
  }
};

/**
 * Computes \f$ \frac{1}{2}\int (\partial_iu_j-\partial_ju_i \f$, where marg=u.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class SkewSymmetricMean : public aol::FENonlinIntegrationVectorInterface < ConfiguratorType,
                                                                           SkewSymmetricMean<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const int _i, _j;
public:
  SkewSymmetricMean ( const typename ConfiguratorType::InitType &Grid,
                      const int I,
                      const int J )
    : aol::FENonlinIntegrationVectorInterface < ConfiguratorType,
                                                  SkewSymmetricMean<ConfiguratorType> > ( Grid ),
      _i(I),
      _j(J) {};

  RealType evaluateIntegrand ( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    typename ConfiguratorType::VecType grad;
    DiscrFuncs[_j].evaluateGradientAtQuadPoint ( El, QuadPoint, grad );
    RealType temp = grad[_i];
    DiscrFuncs[_i].evaluateGradientAtQuadPoint ( El, QuadPoint, grad );
    temp -= grad[_j];
    return 0.5 * temp;
  }
};

/**
 * Projects a deformation such that the mean value of the skew symmetric part of the
 * Jacobian of the deformation is zero. The projection is done in such a way that the
 * deformation at the domain's center of mass is not changed.
 *
 * \note Only works in 2D and assumes that the center of mass is (0.5, 0.5) and the
 *       volume of the domain is 1.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class SymmetricProjector{
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType InitType;
  const InitType &_grid;
  qc::GridDefinition::OldFullNodeIterator _fnit;
  qc::FastILexMapper<ConfiguratorType::Dim> _mapper;
  const RealType _h;
  qc::SkewSymmetricMean<ConfiguratorType> _skewSymmetricMean01;
  const bool _quietMode;
public:
  SymmetricProjector( const InitType &Initializer, const bool QuietMode = true )
   : _grid( Initializer ),
     _mapper( Initializer ),
     _h( Initializer.H() ),
    _skewSymmetricMean01( Initializer, 0, 1 ),
    _quietMode( QuietMode )
  {
    if ( _grid.getDimOfWorld() != 2 )
      throw aol::UnimplementedCodeException( "SymmetricProjector not implemented for Dim != 2", __FILE__, __LINE__ );
  }
  void project( aol::MultiVector<RealType> &Dest ){
    aol::Scalar<RealType> scalarTemp;
    _skewSymmetricMean01.apply( Dest, scalarTemp );
    const RealType s = scalarTemp[0];
    if( !_quietMode ){
      cerr << endl;
      cerr << "before s = " << s << endl;
    }
    for ( _fnit = _grid.begin(); _fnit != _grid.end(); ++_fnit ) {
      const RealType x = (*_fnit)[0] * _h - 0.5;
      const RealType y = (*_fnit)[1] * _h - 0.5;
      Dest[0][ _mapper.getGlobalIndex(*_fnit) ] = Dest[0][ _mapper.getGlobalIndex(*_fnit) ] + s*y;
      Dest[1][ _mapper.getGlobalIndex(*_fnit) ] = Dest[1][ _mapper.getGlobalIndex(*_fnit) ] - s*x;
    }

    if( !_quietMode ) {
      _skewSymmetricMean01.apply( Dest, scalarTemp );
      cerr << "after  s = " << scalarTemp[0] << endl;
    }
  }
};

/**
 * \brief Computes via the method "apply(...)" \f$\frac1{h^d}(\int_{\Omega}(\|D\phi\|,\|cof(D\phi)\|,det(D\phi))^T\psi_idx)_i\f$,
 * where the domain \f$\Omega\f$ is given by a grid (passed to the constructor), \f$\psi_i\f$ denotes the finite element functions,
 * and \f$D\phi\f$ denotes the deformation gradient of a deformation \f$\phi=identity+d\f$, where the displacement \f$d\f$
 * is passed to "apply(...)". Expressed differently, the three invariants of the deformation gradient are evaluated on the grid.
 *
 * \author Wirth
 */
template <typename ConfiguratorType>
class DeformationGradientInvariants :
  public aol::FENonlinVectorOpInterface<ConfiguratorType, ConfiguratorType::Dim, 3, DeformationGradientInvariants<ConfiguratorType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const RealType _gridCellVolume;

public:
  DeformationGradientInvariants( const qc::GridDefinition &Grid ) :
    // load the grid
    aol::FENonlinVectorOpInterface<ConfiguratorType, ConfiguratorType::Dim, 3, DeformationGradientInvariants<ConfiguratorType> >( Grid ),
    _gridCellVolume( std::pow( Grid.H(), ConfiguratorType::Dim ) ) {
  }

  /**
   * \brief Computes \f$(\|D\phi\|,\|cof(D\phi)\|,det(D\phi))^T\f$.
   */
  void getNonlinearity( const aol::auto_container< ConfiguratorType::Dim, aol::DiscreteFunctionDefault< ConfiguratorType > > &DiscFuncs,
                        const typename ConfiguratorType::ElementType &El,
                        int QuadPoint,
                        const typename ConfiguratorType::DomVecType /*&LocalCoord*/,
                        aol::Vec< 3, RealType > &NL ) const {
    // compute the transposed displacement gradient and put it into "dphi"
    typename ConfiguratorType::MatType dphi;
    for ( int j=0; j<ConfiguratorType::Dim; j++ )
      DiscFuncs[j].evaluateGradientAtQuadPoint ( El, QuadPoint, dphi[j] );
    // compute the transposed deformation gradient "dphi"
    for ( int j=0; j<ConfiguratorType::Dim; j++ )
      dphi[j][j] += 1.;
    // return norm, cofactor-norm, and determinant, weighted by the inverse volume/area of a grid cell
    typename ConfiguratorType::MatType cof;
    cof.makeCofactorMatrix( dphi );
    NL[0] = dphi.norm();
    NL[1] = cof.norm();
    NL[2] = dphi.det();
    NL /= _gridCellVolume;
  }
};

} // end namespace qc

#endif
