#ifndef __ELASTOPS_H
#define __ELASTOPS_H

#include <FEOpInterface.h>
#include <tensor.h>

namespace qc {

template <typename ConfiguratorType, qc::Dimension Dim, typename SparseType = qc::FastUniformGridMatrix<typename ConfiguratorType::RealType, qc::QC_2D> >
class InfinitesimalElastOp {};

/**
 * \f[ \int_\Omega \epsilon(u):\epsilon(v) dx \f]
 */
template <typename ConfiguratorType, typename SparseType>
class InfinitesimalElastOp<ConfiguratorType, qc::QC_2D, SparseType >
      : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  SparseType _mat11, _mat12, _mat21, _mat22;
protected:
  aol::BlockOp<RealType> _blockOp;
public:
  InfinitesimalElastOp ( const typename ConfiguratorType::InitType &Grid )
      : _mat11 ( Grid ), _mat12 ( Grid ), _mat21 ( Grid ), _mat22 ( Grid ), _blockOp ( 2, 2 ) {

    aol::FEOpMixedDerivative<ConfiguratorType> xx ( Grid, 0, 0 );
    aol::FEOpMixedDerivative<ConfiguratorType> xy ( Grid, 0, 1 );
    aol::FEOpMixedDerivative<ConfiguratorType> yx ( Grid, 1, 0 );
    aol::FEOpMixedDerivative<ConfiguratorType> yy ( Grid, 1, 1 );

    _mat11.clear();
    _mat12.clear();
    _mat21.clear();
    _mat22.clear();

    yy.assembleAddMatrix ( _mat11 );
    _mat11.scale ( 0.25 );
    xx.assembleAddMatrix ( _mat11 );

    xx.assembleAddMatrix ( _mat22 );
    _mat22.scale ( 0.25 );
    yy.assembleAddMatrix ( _mat22 );

    xy.assembleAddMatrix ( _mat12 );
    _mat12.scale ( 0.25 );

    yx.assembleAddMatrix ( _mat21 );
    _mat21.scale ( 0.25 );

    _blockOp.set ( 0, 0, _mat11 );
    _blockOp.set ( 0, 1, _mat12 );
    _blockOp.set ( 1, 0, _mat21 );
    _blockOp.set ( 1, 1, _mat22 );
  }

  virtual ~InfinitesimalElastOp() { }

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    _blockOp.applyAdd ( Arg, Dest );
  }

};

/** GradDivergence operator
 * \f[ \int_\Omega \mathrm{div}(u)\cdot\mathrm(v) dx \f]
 */
template <typename ConfiguratorType, qc::Dimension Dim>
class GradDivergenceOp : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  aol::FEOpMixedDerivative<ConfiguratorType>* _mixedOps[ Dim*Dim ];
public:
  aol::BlockOp<RealType> _blockOp;
  GradDivergenceOp ( const typename ConfiguratorType::InitType &Grid,
                     aol::OperatorType OpType = aol::ONTHEFLY ) {
    for ( int i = 0; i < Dim; i++ ) {
      for ( int j = 0; j < Dim; j++ ) {
        _mixedOps[ i*Dim+j ] = new aol::FEOpMixedDerivative<ConfiguratorType> ( Grid, j, i, OpType );
        _blockOp.set ( i, j, _mixedOps[ i*Dim+j ] );
      }
    }
  }

  virtual ~GradDivergenceOp() {
    for ( int i = 0; i < Dim; i++ ) {
      for ( int j = 0; j < Dim; j++ ) {
        delete _mixedOps[ i*Dim+j ];
      }
    }
  }

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    _blockOp.applyAdd ( Arg, Dest );
  }

};

/**
 * Computes \f$ \frac{1}{2} (Du+Du^T) \f$, where DiscrFuncs=u.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
void calcSymmetrizedJacobian ( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               const int QuadPoint,
                               typename ConfiguratorType::MatType &SymJacobian ) {
  typename ConfiguratorType::VecType grad;
  SymJacobian.setZero();

  for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
    DiscrFuncs[i].evaluateGradientAtQuadPoint ( El, QuadPoint, grad );
    for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
      SymJacobian[i][j] += grad[j];
      SymJacobian[j][i] += grad[j];
    }
  }
  SymJacobian *= 0.5;
}

/**
 * Computes \f$ \frac{1}{2}\int C\epsilon(u):\epsilon(u) \f$, where marg=u.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class LinearizedElasticEnergy : public aol::FENonlinIntegrationVectorInterface < ConfiguratorType,
                                                                                 LinearizedElasticEnergy<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::ElasticTensor<RealType> &_c;
public:
  LinearizedElasticEnergy ( const typename ConfiguratorType::InitType &Grid, const aol::ElasticTensor<RealType> &C )
    : aol::FENonlinIntegrationVectorInterface < ConfiguratorType,
                                                LinearizedElasticEnergy<ConfiguratorType> > ( Grid ),
      _c ( C ) {};

  RealType evaluateIntegrand ( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    typename ConfiguratorType::VecType grad;
    typename ConfiguratorType::MatType matrix, matrix2;

    calcSymmetrizedJacobian<ConfiguratorType> ( DiscrFuncs, El, QuadPoint, matrix );
    _c.apply ( matrix, matrix2 );
    return 0.5 * matrix2.ddprod ( matrix );
  }
};

/**
 * Computes the first variation of LinearizedElasticEnergy, e.q. \f$ int C\epsilon(u):D\zeta \f$, where marg=u.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class VariationOfLinearizedElasticEnergy
  : public aol::FENonlinVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, VariationOfLinearizedElasticEnergy<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::ElasticTensor<RealType> &_c;
public:
  VariationOfLinearizedElasticEnergy ( const typename ConfiguratorType::InitType &Grid, const aol::ElasticTensor<RealType> &C )
    : aol::FENonlinVectorDiffOpInterface<ConfiguratorType, ConfiguratorType::Dim, ConfiguratorType::Dim, VariationOfLinearizedElasticEnergy<ConfiguratorType> > ( Grid ),
      _c ( C ) {};

  void getNonlinearity ( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                         const typename ConfiguratorType::ElementType &El,
                         int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/,
                         typename aol::Mat<ConfiguratorType::Dim, ConfiguratorType::Dim, RealType> &NL ) const {
    typename ConfiguratorType::MatType matrix;
    calcSymmetrizedJacobian<ConfiguratorType> ( DiscFuncs, El, QuadPoint, matrix );
    _c.apply ( matrix, NL );
  }
};

} // end namespace qc

#endif
