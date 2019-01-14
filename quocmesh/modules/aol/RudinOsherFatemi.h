#ifndef __RUDINOSHERFATEMI_H
#define __RUDINOSHERFATEMI_H

#include <mcm.h>

namespace aol {

//! Computes \f$ \int | \nabla u |_{2,\delta} \f$, where arg=u.
template <typename ConfiguratorType>
class IsoEnergyOp
: public aol::FENonlinIntegrationScalarInterface<ConfiguratorType, IsoEnergyOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const RealType _deltaSqr;
public:
  IsoEnergyOp( const typename ConfiguratorType::InitType &Initializer,
               const RealType Delta )
  : aol::FENonlinIntegrationScalarInterface< ConfiguratorType,
                                             IsoEnergyOp<ConfiguratorType> > ( Initializer ),
    _deltaSqr(aol::Sqr(Delta)){
  }

  RealType evaluateIntegrand( const aol::DiscreteFunctionDefault< ConfiguratorType > &DiscFuncs,
                              const typename ConfiguratorType::ElementType &El,
                              int QuadPoint, const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    typename ConfiguratorType::VecType gradU;
    DiscFuncs.evaluateGradientAtQuadPoint( El, QuadPoint, gradU );

    return sqrt(gradU.normSqr() + _deltaSqr);
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename StiffOpType = qc::MCMStiffOp<ConfiguratorType> >
class VariationOfIsoEnergyOp : public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_grid;
  const RealType _delta;
public:
  VariationOfIsoEnergyOp ( const typename ConfiguratorType::InitType &Initializer,
                           const RealType Delta )
  : _grid ( Initializer ),
    _delta ( Delta ) {}

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    StiffOpType mcmStiffOp( _grid, aol::ONTHEFLY, _delta );
    mcmStiffOp.setImageReference( Arg );
    mcmStiffOp.applyAdd(Arg, Dest);
  }
};

/**
 * \todo Essentially the same as qc::WillmoreProjectionStiffOp. Drop one of these classes and
 *       find a better name for the remaining class.
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType>
class IsotropicMatrixStiffOp : public aol::FELinMatrixWeightedStiffInterface<ConfiguratorType,
                                           IsotropicMatrixStiffOp<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _discrImage;
  const RealType _eps;
public:
  IsotropicMatrixStiffOp( const typename ConfiguratorType::InitType &Initializer,
                          const aol::Vector<RealType> &Image,
                          const RealType Epsilon,
                          aol::OperatorType OpType = aol::ONTHEFLY )
    : aol::FELinMatrixWeightedStiffInterface<ConfiguratorType, IsotropicMatrixStiffOp<ConfiguratorType> >
      ( Initializer, OpType ),
      _discrImage( Initializer, Image ),
      _eps( Epsilon )  {  }

  // implement the projection-matrix:
  inline void getCoeffMatrix( const qc::Element &El, int QuadPoint,
      const typename ConfiguratorType::VecType& /*RefCoord*/, typename ConfiguratorType::MatType &mat ) const {

    typename ConfiguratorType::VecType grad;
    _discrImage.evaluateGradientAtQuadPoint( El, QuadPoint, grad );

    const RealType norm = sqrt( grad.normSqr() + aol::Sqr( _eps ) );
    const RealType normCub = aol::Cub( norm );

    for( int i = 0; i < ConfiguratorType::Dim; i++){
      mat[i][i] = 1./norm - aol::Sqr(grad[i]) / normCub;
      for( int j = i+1; j < ConfiguratorType::Dim; j++){
        mat[i][j] = - grad[i] * grad[j] / normCub;
      }
      for( int j = 0; j < i; j++){
        mat[i][j] = mat[j][i];
      }
    }
  }
};


//! Computes \f$ \int_\Omega \frac{\beta_1}{2}( u - u_0 )^2 + \beta_2 | \nabla u |_{2,\delta}dx \f$
template <typename ConfiguratorType, typename MatrixType = typename ConfiguratorType::MatrixType, typename RegOpType = IsoEnergyOp<ConfiguratorType> >
class IsotropicROFEnergy : public aol::Op<aol::Vector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const MatrixType &_massMatrix;
  aol::Vector<RealType> _u0constantVec; // stores -beta_1*M*u_0
  RealType _u0constant;
  const RealType _beta1;
  const RealType _beta2;
  const RealType _delta;
public:
  IsotropicROFEnergy ( const typename ConfiguratorType::InitType &Initializer,
                       const aol::Vector<RealType> &U0,
                       const RealType Beta1,
                       const RealType Beta2,
                       const RealType Delta,
                       const MatrixType &MassMatrix):
  _grid(Initializer),
  _massMatrix(MassMatrix),
  _u0constantVec(U0, aol::STRUCT_COPY),
  _beta1(Beta1),
  _beta2(Beta2),
  _delta(Delta)
  {
    MassMatrix.apply(U0, _u0constantVec);
    _u0constantVec *= -2.;
    _u0constant = (-1.)*(_u0constantVec * U0);
  }
  virtual void apply( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const{
    //! Isotropic BV-term
    RegOpType _isoEnergyOp( _grid, _delta );
    _isoEnergyOp.apply( Arg, Dest );
    Dest *= _beta2;

    //! Fidelity-term
    aol::Vector<RealType> temp( Arg, aol::STRUCT_COPY );
    RealType fidelityEnergy = _u0constant;
    fidelityEnergy += _u0constantVec*Arg;
    _massMatrix.apply( Arg, temp );
    fidelityEnergy += Arg*temp;

    fidelityEnergy = _beta1*0.5*fidelityEnergy;
    Dest[0] += fidelityEnergy;

  }
  virtual void applyAdd( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const{
    aol::Scalar<RealType> tmp;
    apply( Arg, tmp );
    Dest += tmp;
  }
};

template <typename ConfiguratorType, typename MatrixType = typename ConfiguratorType::MatrixType, typename VarOfRegOpType = VariationOfIsoEnergyOp<ConfiguratorType> >
class VariationOfIsotropicROF : public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const MatrixType &_massMatrix;
  aol::Vector<RealType> _u0constantVec; // stores -beta_1*M*u_0
  const RealType _beta1;
  const RealType _beta2;
  const RealType _delta;
  typename qc::BitArray<ConfiguratorType::Dim> *_pDirichletMask;
public:
  /* MassMatrix needs to be the full mass matrix, regardless if DirichletMask
   * is supplied or not.
   */
  VariationOfIsotropicROF ( const typename ConfiguratorType::InitType &Initializer,
                            const aol::Vector<RealType> &U0,
                            const RealType Beta1,
                            const RealType Beta2,
                            const RealType Delta,
                            const MatrixType &MassMatrix,
                            typename qc::BitArray<ConfiguratorType::Dim> *DirichletMask = NULL ):
  _grid(Initializer),
  _massMatrix(MassMatrix),
  _u0constantVec(U0, aol::STRUCT_COPY),
  _beta1(Beta1),
  _beta2(Beta2),
  _delta(Delta),
  _pDirichletMask(DirichletMask)
  {
    MassMatrix.apply(U0, _u0constantVec);
    _u0constantVec *= -1.;
  }
  virtual void apply( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const{
    Dest = _u0constantVec; // Dest=-M*u_0
    _massMatrix.applyAdd(Arg, Dest); // Dest = beta1/beta2*(M*u-M*u_0)
    Dest *= _beta1/_beta2;

    VarOfRegOpType tvDE ( _grid, _delta); // Dest = beta1(M*u-M*u_0)+beta2*L_{mcm}u
    tvDE.applyAdd ( Arg, Dest );
    Dest *= _beta2;

    if ( _pDirichletMask ){
      for ( int i = 0; i < _pDirichletMask->size(); ++i ) {
        if ( (*_pDirichletMask)[i] == true ) {
          Dest[i] = 0;
        }
      }
    }
  }
  virtual void applyAdd( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const{
    aol::Vector<RealType> tmp( Arg, aol::STRUCT_COPY );
    apply( Arg, tmp );
    Dest += tmp;
  }
};




template <typename ConfiguratorType, typename MatrixType = typename ConfiguratorType::MatrixType>
class SecondVariationOfIsotropicROF : public aol::Op<aol::Vector<typename ConfiguratorType::RealType>, MatrixType > {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const MatrixType &_massMatrix;
  const RealType _beta1;
  const RealType _beta2;
  const RealType _delta;
  typename qc::BitArray<ConfiguratorType::Dim> *_pDirichletMask;
public:
  /* If you supply DirichletMask, the MassMatrixPossiblyMasked needs to already respect the DirichletMask.
   * Without DirichletMask, MassMatrixPossiblyMasked should be the standard mass matrix.
   */
  SecondVariationOfIsotropicROF ( const typename ConfiguratorType::InitType &Initializer,
                                  const RealType Beta1,
                                  const RealType Beta2,
                                  const RealType Delta,
                                  const MatrixType &MassMatrixPossiblyMasked,
                                  typename qc::BitArray<ConfiguratorType::Dim> *DirichletMask = NULL ):
  _grid(Initializer),
  _massMatrix(MassMatrixPossiblyMasked),
  _beta1(Beta1),
  _beta2(Beta2),
  _delta(Delta),
  _pDirichletMask(DirichletMask)
  {
  }
  virtual void apply( const aol::Vector<RealType> &Arg, MatrixType &Dest ) const{
    Dest.setZero(); //Dest = beta2/beta1*L[..]
    IsotropicMatrixStiffOp<ConfiguratorType> isotropicMatrixStiffOp(_grid, Arg, _delta);
    if ( _pDirichletMask == NULL )
      isotropicMatrixStiffOp.assembleAddMatrix(Dest);
    else
      isotropicMatrixStiffOp.assembleAddMatrix(Dest, _pDirichletMask);
    Dest *= (_beta2/_beta1);
    Dest += _massMatrix; //Dest = beta1*M + beta2L[..]
    Dest *= (_beta1);
  }
  virtual void applyAdd( const aol::Vector<RealType> &/*Arg*/, MatrixType &/*Dest*/ ) const{
    throw aol::Exception( "Not implemented", __FILE__, __LINE__);
  }
};

} // end namespace aol

#endif // __RUDINOSHERFATEMI_H
