#ifndef __FIRSTORDERTVALGOS_H
#define __FIRSTORDERTVALGOS_H

#include <finiteDifferences.h>
#include <projectors.h>

namespace qc {

/**
 * \brief Minimizing algorithm for problems of the form \f$ F(x)+\gamma G(\nabla_h x) \f$, where \f$ F,G \f$ are convex,
 *        one of them strictly convex and \f$ \nabla_h \f$ is the gradient approximated with forward finite-differences.
 *        
 *        For an example on how to use this algorithm, see FirstOrderPrimalDualROFMinimizer in modules/imgproc/denoising.h
 *        or FirstOrderPrimalDualTwoPhaseMSSegmentor and FirstOrderPrimalDualTwoPhaseMSSegmentorEngine in modules/imgproc/primalDualSegmentation.h
 *
 * \note \arg Gamma is the weight of \f$ G \f$, not the \f$ \gamma \f$ from the Chambolle Pock algorithm.
 *
 * \author Berkels
 * \ingroup Optimization
 */
template <typename ConfiguratorType>
class FirstOrderChambollePockTVAlgorithmType2 : public TVAlgorithmBase<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  FirstOrderChambollePockTVAlgorithmType2 ( const typename ConfiguratorType::InitType &Initializer,
                                            const RealType Gamma,
                                            const int MaxIterations = 1000,
                                            const RealType StopEpsilon = 0 )
    : TVAlgorithmBase<ConfiguratorType> ( Initializer, Gamma, MaxIterations, StopEpsilon ) {}

  virtual ~FirstOrderChambollePockTVAlgorithmType2 () {}

protected:
  virtual void applyResolventOfDataTermSingle ( const RealType TauOverGamma, ArrayType &ArgDest ) const = 0;

  virtual void applyResolventOfAdjointRegTermSingle ( const RealType /*Sigma*/, qc::MultiArray<RealType, ConfiguratorType::Dim> &ArgDest ) const {
    const int primalDOFs = this->_grid.getNumberOfNodes();
    typename ConfiguratorType::VecType pVec;
#ifdef _OPENMP
#pragma omp parallel for firstprivate ( pVec )
#endif
    for ( int j = 0; j < primalDOFs; ++j ) {
      for ( int i = 0; i < ConfiguratorType::Dim; ++i )
        pVec[i] = ArgDest[i][j];
      const RealType tempVal = aol::Max ( pVec.norm(), aol::ZOTrait<RealType>::one );
      for ( int i = 0; i < ConfiguratorType::Dim; ++i ) {
        ArgDest[i][j] = ArgDest[i][j] / tempVal;
      }
    }
  }

  virtual void initializeDualVariable ( qc::MultiArray<RealType, ConfiguratorType::Dim> &P ) const {
    P.setZero();
  }

  virtual void saveStep ( const ArrayType &/*U*/, const int /*Iterations*/ ) const { }

public:
  void minimize ( ArrayType &U, qc::MultiArray<RealType, ConfiguratorType::Dim> *PDual = NULL ) const {
    qc::MultiArray<RealType, ConfiguratorType::Dim> p ( this->_grid );
    if ( PDual == NULL )
      initializeDualVariable ( p );
    else
      p = *PDual;
    qc::MultiArray<RealType, ConfiguratorType::Dim> divTemp ( this->_grid );
    ArrayType uOld ( U, aol::FLAT_COPY );
    ArrayType uNew ( this->_grid );
    RealType tau = aol::ZOTrait<RealType>::one / 8;
    RealType sigma = tau;
    RealType theta = 0.0;
    const RealType gammaHDependent = this->_gamma / this->_grid.H();
    const RealType pockGamma =  0.7 / gammaHDependent;

    aol::ProgressBar<> progressBar ( "Minimizing" );
    if ( !this->_quietMode ) {
      progressBar.start ( this->_maxIterations );
      progressBar.display ( cerr );
    }

    for ( int iterations = 0; iterations < this->_maxIterations; ++iterations ) {
      // Update p
      qc::calculateForwardFDGradient<RealType, ConfiguratorType::Dim> ( uOld, divTemp );
      p.addMultiple ( divTemp, sigma );

      applyResolventOfAdjointRegTermSingle ( sigma, p );

      // Calc uNew
      qc::calculateBackwardFDDivergence<RealType, ConfiguratorType::Dim> ( p, uNew );
      uNew.scaleAndAdd ( tau, uOld );

      applyResolventOfDataTermSingle ( tau / gammaHDependent, uNew );

      // Update parameters
      theta = 1 / sqrt ( 1 + 2 * pockGamma *tau );
      tau *= theta;
      sigma /= theta;

      uNew.scaleAndAddMultiple ( ( 1 + theta ), uOld, -theta );

      uOld -= uNew;
      const RealType change = uOld.norm();

      uOld = uNew;

// TODO: implement step saving for VectorContainer
//      if ( this->_pStepSaver ) {
//        this->_pStepSaver->saveStep ( uOld, iterations );
//      }

      saveStep ( uOld, iterations );

      // If the change is small enough, we consider the algorithm to be converged
      if ( change < this->_stopEpsilon ) {
        if ( !this->_quietMode ) cerr << "\nStopping after " << iterations + 1 << " iterations.\n";
        break;
      }

      if ( !this->_quietMode ) progressBar++;
    }
    if ( !this->_quietMode ) progressBar.finish();

    if ( PDual != NULL )
      *PDual = p;
  }
};
  
/**
 * \brief For examples on how to use this algorithm, see FirstOrderPrimalDualMultiPhaseMSSegmentor and FirstOrderPrimalDualMultiPhaseMSSegmentorEngine in modules/imgproc/primalDualSemgentation.h
 *
 * \author Mevenkamp
 * \ingroup Optimization
 */
template <typename ConfiguratorType>
class FirstOrderChambollePockTVAlgorithmType1 : public TVAlgorithmBase<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  
  FirstOrderChambollePockTVAlgorithmType1 ( const typename ConfiguratorType::InitType &Initializer,
                                            const RealType Gamma,
                                            const int MaxIterations = 1000,
                                            const RealType StopEpsilon = 0 )
    : TVAlgorithmBase<ConfiguratorType> ( Initializer, Gamma, MaxIterations, StopEpsilon ) { }
  
  virtual ~FirstOrderChambollePockTVAlgorithmType1 () {}
  
protected:
  virtual void applyResolventOfDataTermSingle ( const RealType TauOverGamma, aol::VectorContainer<ArrayType> &ArgDest ) const = 0;
  
  virtual void applyResolventOfAdjointRegTermSingle ( const RealType /*Sigma*/, aol::VectorContainer<qc::MultiArray<RealType, ConfiguratorType::Dim> > &ArgDest ) const {
    const int primalDOFs = this->_grid.getNumberOfNodes();
    typename ConfiguratorType::VecType pVec;
#ifdef _OPENMP
#pragma omp parallel for firstprivate ( pVec )
#endif
    for ( int k = 0; k<ArgDest.size ( ) ; ++k ) {
      for ( int j = 0; j < primalDOFs; ++j ) {
        for ( int i = 0; i < ConfiguratorType::Dim; ++i )
          pVec[i] = ArgDest[k][i][j];
        const RealType tempVal = aol::Max ( pVec.norm(), aol::ZOTrait<RealType>::one );
        for ( int i = 0; i < ConfiguratorType::Dim; ++i )
          ArgDest[k][i][j] = ArgDest[k][i][j] / tempVal;
      }
    }
  }
  
  virtual void initializeDualVariable ( aol::VectorContainer<qc::MultiArray<RealType, ConfiguratorType::Dim> > &P ) const {
    P.setZero();
  }
  
  virtual void saveStep ( const aol::VectorContainer<ArrayType> &/*U*/, const int /*Iterations*/ ) const { }
  
public:
  void minimize ( aol::VectorContainer<ArrayType> &U, aol::VectorContainer<qc::MultiArray<RealType, ConfiguratorType::Dim> > *PDual = NULL ) const {
    aol::VectorContainer<qc::MultiArray<RealType, ConfiguratorType::Dim> > p ( U.size ( ), qc::MultiArray<RealType, ConfiguratorType::Dim> ( this->_grid ) );
    if ( PDual == NULL )
      initializeDualVariable ( p );
    else
      p = *PDual;
    aol::VectorContainer<qc::MultiArray<RealType, ConfiguratorType::Dim> > divTemp ( U.size ( ), qc::MultiArray<RealType, ConfiguratorType::Dim> ( this->_grid ) );
    aol::VectorContainer<ArrayType> uOld ( U );
    
    aol::VectorContainer<ArrayType> uNew ( U.size ( ), ArrayType ( this->_grid ) );
    RealType tau = aol::ZOTrait<RealType>::one / 8;
    RealType sigma = tau;
    RealType theta = 0.7;
    const RealType gammaHDependent = this->_gamma / this->_grid.H();
    
    aol::ProgressBar<> progressBar ( "Minimizing" );
    if ( !this->_quietMode ) {
      progressBar.start ( this->_maxIterations );
      progressBar.display ( cerr );
    }
    
    for ( int iterations = 0; iterations < this->_maxIterations; ++iterations ) {
      // Update p
      for ( int k=0; k<uOld.size ( ) ; ++k )
        qc::calculateForwardFDGradient<RealType, ConfiguratorType::Dim> ( uOld[k], divTemp[k] );
      p.addMultiple ( divTemp, sigma );
      
      applyResolventOfAdjointRegTermSingle ( sigma, p );
      
      // Calc uNew
      for ( int k=0; k<uNew.size ( ) ; ++k )
        qc::calculateBackwardFDDivergence<RealType, ConfiguratorType::Dim> ( p[k], uNew[k] );
      uNew.scaleAndAdd ( tau, uOld );
      
      applyResolventOfDataTermSingle ( tau / gammaHDependent, uNew );
      
      uNew.scaleAndAddMultiple ( ( 1 + theta ), uOld, -theta );
      
      uOld -= uNew;
      const RealType change = uOld.norm();
      
      uOld = uNew;
      
// TODO: extend StepSaver to Vectors
//      if ( this->_pStepSaver ) {
//        this->_pStepSaver->saveStep ( uOld, iterations );
//      }
      
      saveStep ( uOld, iterations );
      
      // If the change is small enough, we consider the algorithm to be converged
      if ( change < this->_stopEpsilon ) {
        if ( !this->_quietMode ) cerr << "\nStopping after " << iterations + 1 << " iterations.\n";
        break;
      }
      
      if ( !this->_quietMode ) progressBar++;
    }
    if ( !this->_quietMode ) progressBar.finish();
    
    if ( PDual != NULL )
      *PDual = p;
    
    U = uOld; // This is required, since aol::VectorContainer does not support FLAT_COPY
  }
};
  
} // namespace qc

#endif // __FIRSTORDERTVALGOS_H
