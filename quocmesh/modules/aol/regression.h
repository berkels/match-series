#ifndef __REGRESSION_H
#define __REGRESSION_H

#include <matrixInverse.h>
#include <preconditioner.h>
#include <ArmijoSearch.h>
#include <eigenWrapper.h>


namespace aol {

  
template <typename _RealType>
class LinearRegression {
  typedef _RealType RealType;
  typedef aol::Vector<_RealType> VectorType;
  typedef aol::FullMatrix<_RealType> MatrixType;
protected:
  MatrixType _A;
public:
  LinearRegression ( const MatrixType &SystemMatrix )
  : _A ( SystemMatrix ) { }
  
  virtual ~LinearRegression ( ) { }
  
  virtual void apply ( const VectorType &/*RHS*/, VectorType &/*X*/ ) = 0;
};


template <typename _RealType>
class LinearRegressionQR : public LinearRegression<_RealType> {
  typedef _RealType RealType;
  typedef aol::Vector<_RealType> VectorType;
  typedef aol::FullMatrix<_RealType> MatrixType;
protected:
  const aol::QRInverse<RealType> _qrInv;
public:
  LinearRegressionQR ( const MatrixType &SystemMatrix )
  : LinearRegression<RealType> ( SystemMatrix ), _qrInv ( SystemMatrix ) { }
  
  void apply ( const VectorType &RHS, VectorType &X ) {
    if ( RHS.size ( ) != this->_A.getNumRows ( ) )
      throw aol::Exception ( "Number of equations does not match between system matrix and right-hand side!", __FILE__, __LINE__ );
    
    if ( X.size ( ) != this->_A.getNumCols ( ) )
      throw aol::Exception ( "Vector size does not match number of columns of the system matrix!", __FILE__, __LINE__ );
    
    _qrInv ( RHS, X );
    
    if ( X.checkForNANsAndINFs() )
      throw aol::Exception ( "Linear regression failed! NaN / Inf entries detected in the solution!", __FILE__, __LINE__ );
  }
};


template <typename _RealType>
class LinearRegressionNormalEquations : public LinearRegression<_RealType> {
  typedef _RealType RealType;
  typedef aol::Vector<_RealType> VectorType;
  typedef aol::FullMatrix<_RealType> MatrixType;
protected:
  MatrixType _ATA;
  VectorType _ATb;
  aol::DiagonalMatrix<RealType> _P;
public:
  LinearRegressionNormalEquations ( const MatrixType &SystemMatrix )
  : LinearRegression<RealType> ( SystemMatrix ), _ATA ( SystemMatrix.getNumCols ( ), SystemMatrix.getNumCols ( ) ), _ATb ( SystemMatrix.getNumCols ( ) ),
  _P ( SystemMatrix.getNumCols ( ) ) {
    // Calculate _A^T * _A
    for ( int i=0; i<this->_A.getNumCols ( ) ; ++i ) {
      for ( int j=0; j<this->_A.getNumCols ( ) ; ++j ) {
        _ATA.set ( i, j, 0 );
        for ( int k=0; k<this->_A.getNumRows ( ) ; ++k )
          _ATA.add ( i, j, this->_A.get ( k, i ) * this->_A.get ( k, j ) );
      }
    }
    
    // Apply preconditioner to system matrix
    aol::DiagonalPreconditioner<aol::Vector<RealType> > preconditioner ( _ATA );
    _P.resize ( _ATA.getNumRows ( ) );
    preconditioner.getPreconditionerMatrix ( _P );
    for ( int i=0; i<_ATA.getNumRows ( ) ; ++i ) {
      for ( int j=0; j<_ATA.getNumCols ( ) ; ++j )
        _ATA.set ( i, j, _ATA.get ( i, j ) * _P.get ( i, i ) );
    }
  }
  
  void apply ( const VectorType &RHS, VectorType &X ) {
    if ( RHS.size ( ) != this->_A.getNumRows ( ) )
      throw aol::Exception ( "Number of equations does not match between system matrix and right-hand side!", __FILE__, __LINE__ );
    
    if ( X.size ( ) != this->_A.getNumCols ( ) )
      throw aol::Exception ( "Vector size does not match number of columns of the system matrix!", __FILE__, __LINE__ );
    
    // Calculate _A^T * _b and apply preconditioner
    _ATb.setZero ( );
    this->_A.applyAddTranspose ( RHS, _ATb );
    this->_P.apply ( _ATb, _ATb );
    
    // Solve P ATA X = P ATb using QR decomposition
    aol::QRInverse<RealType> qrInverse ( _ATA );
    qrInverse.apply ( _ATb, X );
    
    if ( X.checkForNANsAndINFs() )
      throw aol::Exception ( "Linear regression failed! NaN / Inf entries detected in the solution!", __FILE__, __LINE__ );
  }
};
  
  

  
  
/**
 *  \brief Implementation of the Levenberg-Marquardt algorithm:
 *         K. Levenberg: A Method for the Solution of Certain Problems in Least Squares, Quart. Appl. Math. 2, 164-168, 1944.
 *         D. Marquardt: An Algorithm for Least-Squares Estimation of Nonlinear Parameters, SIAM J. Appl. Math. 11, 431-441, 1963.
 *  \author mevenkamp
 *  \ingroup Optimization
 */
template <typename _RealType, typename _MatrixType = aol::FullMatrix<_RealType>, typename _LinearRegressionType = LinearRegressionQR<_RealType> >
class LevenbergMarquardtAlgorithm {
  typedef _RealType RealType;
  typedef _MatrixType MatrixType;
  typedef _LinearRegressionType LinearRegressionType;
  
protected:
  //! stores the evaluation of DF and is deleted in the destructor.
  mutable aol::DeleteFlagPointer<MatrixType> _jacobian;
  
  const int _dimRangeF;
  const aol::Op<aol::Vector<RealType> > &_F;
  const aol::Op<aol::Vector<RealType>, MatrixType> &_DF;
  const int _maxIterations;
  const RealType _mu0, _muMin, _muMax, _rho0, _rho1;
  const RealType _epsDeltaX, _epsF, _epsGradRealTargetFunc;
  const bool _verbose;
  
public:
  LevenbergMarquardtAlgorithm ( const int DimRangeF,
                               const aol::Op<aol::Vector<RealType> > &F,
                               const aol::Op<aol::Vector<RealType>, MatrixType> &DF,
                               const int MaxIterations = 50,
                               const RealType Mu0 = 1,
                               const RealType Rho0 = 0.2,
                               const RealType Rho1 = 0.8,
                               const RealType EpsDeltaX = 1e-6,
                               const RealType EpsF = 1e-6,
                               const RealType EpsGradRealTargetFunc = 1e-6,
                               const bool Verbose = true )
  : _jacobian ( NULL ),
  _dimRangeF ( DimRangeF ),
  _F ( F ),
  _DF ( DF ),
  _maxIterations ( MaxIterations ),
  _mu0 ( Mu0 ), _muMin ( 1e-6 ), _muMax ( 1e6 ),
  _rho0 ( Rho0 ), _rho1 ( Rho1 ),
  _epsDeltaX ( EpsDeltaX ), _epsF ( EpsF ), _epsGradRealTargetFunc ( EpsGradRealTargetFunc ),
  _verbose ( Verbose ) { }
  
  virtual ~LevenbergMarquardtAlgorithm ( ) { }
  
  void apply( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( _jacobian.get ( ) == NULL )
      _jacobian.reset ( new MatrixType ( _dimRangeF, Arg.size ( ) ), true );
    
    if ( _verbose )
      std::cerr << "LevenbergMarquardt: Optimizing..." << std::endl;
    
    RealType mu = _mu0;
    int iteration = 0;
    aol::Vector<RealType> gradientRealValuedTargetFunction ( _jacobian->getNumCols ( ) );
    
    aol::Vector<RealType> x0 ( Arg ), x1 ( Arg );
    aol::Vector<RealType> f0 ( _dimRangeF ), f1 ( _dimRangeF );
    aol::Vector<RealType> direction ( Dest, aol::STRUCT_COPY );
    
    _F.apply ( x0, f0 );
    f1 = f0;
    
    if ( _verbose )
      std::cerr << "LevenbergMarquardt: Iteration " << iteration << ": x=" << x1 << "; res=" << f1.norm ( ) << std::endl;
    
    do {
      x0 = x1;
      f0 = f1;
      
      _DF.apply ( x0, *_jacobian );
      
      // Setup right-hand side
      aol::Vector<RealType> rhs ( _jacobian->getNumRows ( ) + _jacobian->getNumCols ( ) );
      for ( int i=0; i<f0.size ( ) ; ++i )
        rhs[i] = -f0[i];
      
      do {
        if ( mu < _muMin || mu > _muMax ) {
          if ( _verbose )
            std::cerr << "LevenbergMarquardt: mu outside bounds [ " << _muMin << ", " << _muMax << " ]. Breaking inner iteration." << std::endl;
          break;
        }
          
        // Setup system matrix
        aol::FullMatrix<RealType> systemMatrix ( _jacobian->getNumRows ( ) + _jacobian->getNumCols ( ), _jacobian->getNumCols ( ) );
        for ( int i=0; i<_jacobian->getNumRows ( ) ; ++i ) {
          for ( int j=0; j<_jacobian->getNumCols ( ) ; ++j ) {
            systemMatrix.set ( i, j, _jacobian->get ( i, j ) );
            systemMatrix.set ( _jacobian->getNumRows ( ) + j, j, mu );
          }
        }
        
        LinearRegressionType linearRegression ( systemMatrix );
        linearRegression.apply ( rhs, direction );
        
        step ( x0, direction, x1 );
        _F.apply ( x1, f1 );
        
        // Test the step (and possibly correct mu)
        aol::Vector<RealType> linEstimate ( f0 );
        _jacobian->applyAdd ( direction, linEstimate );
        RealType linDiffNormSqr = f0.normSqr ( ) - linEstimate.normSqr ( );
        if ( linDiffNormSqr < 0 ) {
          if ( _verbose )
            std::cerr << "LevenbergMarquardt: Linear residual is negative. Breaking inner iteration." << std::endl;
          break;
        }
        RealType rhoMu = aol::NumberTrait<RealType>::Inf;
        rhoMu = ( f0.normSqr ( ) - f1.normSqr ( ) ) / linDiffNormSqr;
        
        if ( !aol::isFinite ( rhoMu ) ) {
          mu *= 10;
          if ( _verbose )
            std::cerr << "LevenbergMarquardt: Quotient between linear and non-linear Residual is NaN or Inf. Increasing to mu=" << mu << " and recalculating step." << std::endl;
          x1 = x0;
          _F.apply ( x1, f1 );
          continue;
        }
        
        if ( rhoMu <= _rho0 ) {
          mu *= 10;
          if ( _verbose )
            std::cerr << "LevenbergMarquardt: Adapting damping parameter: Increasing to mu=" << mu << " and recalculating step." << std::endl;
        } else if ( _rho0 < rhoMu && rhoMu < _rho1 ) {
          if ( _verbose )
            std::cerr << "LevenbergMarquardt: Accepting step without changing the damping parameter." << std::endl;
          break;
        } else {
          mu /= 10;
          if ( _verbose )
            std::cerr << "LevenbergMarquardt: Adapting damping parameter: Decreasing to mu=" << mu << " and accepting step." << std::endl;
          break;
        }
      } while ( true );
      
      ++iteration;
      if ( _verbose )
        std::cerr << "LevenbergMarquardt: Iteration " << iteration << ": x=" << x1 << "; res=" << f1.norm ( ) << std::endl;
      
      // Calculate gradient of scalar target function phi(x) = 1/2 * ( F(x) ).normSqr ( )
      // Use: grad phi(x) = ( F'(x) )^T F(x)
      for ( int i=0; i<_jacobian->getNumCols ( ) ; ++i ) {
        gradientRealValuedTargetFunction[i] = 0;
        for ( int k=0; k<_jacobian->getNumRows ( ) ; ++k )
          gradientRealValuedTargetFunction[i] += _jacobian->get ( k, i ) * f1[k];
      }
    } while ( iteration < _maxIterations && f1.norm ( ) >= _epsF
             && gradientRealValuedTargetFunction.norm ( ) >= _epsGradRealTargetFunc
             && direction.norm ( ) / x0.norm ( ) >= _epsDeltaX );
    
    if ( _verbose ) {
      std::cerr << "LevenebrgMarquardt: Terminating iteration. Reason(s):";
      if ( iteration == _maxIterations )
        std::cerr << " maximum number of iterations reached";
      if ( f1.norm ( ) < _epsF )
        std::cerr << " residual below threshold";
      if ( gradientRealValuedTargetFunction.norm ( ) < _epsGradRealTargetFunc )
        std::cerr << " gradient of target function below threshold";
      if ( direction.norm ( ) / x0.norm ( ) < _epsDeltaX )
        std::cerr << " step size below threshold";
      std::cerr << "." << std::endl;
    }
    
    if ( _verbose ) {
      std::cerr << "LevenbergMarquardt: Optimization finished." << std::endl;
      std::cerr << "LevenbergMarquardt: Residual: " << f1.norm ( ) << std::endl;
      std::cerr << "LevenbergMarquardt: Optimal parameters: " << x1 << std::endl;
    }
    
    Dest = x1;
  }
  
  void applyAdd( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__);
  }
  
protected:
  virtual void step ( const aol::Vector<RealType> &X0, const aol::Vector<RealType> &Direction, aol::Vector<RealType> &X1 ) const {
    X1 = X0;
    X1 += Direction;
  }
};

/**
 *  \brief Combination of Levenberg-Marquardt and projections for solving non-linear least-squares problems with convex constraints
 *
 *  In each iteration solves unconstrained linearized regression problem \f$\|DF(x_k) s_k + F(x_k)\|^2 + \mu^2 \|x_k\|^2 \rightarrow \min\f$.
 *  The projection of \f$x_k+1 = x_k + s_k\f$ onto the convex space \f$X = \{ x \in \mathbb{R}^n | A x <= B \}\f$ is then used as the update
 *
 *  \author mevenkamp
 *  \ingroup Optimization
 */
template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ProjectorType>
class ConstrainedLevenbergMarquardtAlgorithm : public LevenbergMarquardtAlgorithm<_RealType, _MatrixType, _LinearRegressionType> {
  typedef _RealType RealType;
  typedef aol::Vector<RealType> VectorType;
  typedef _MatrixType MatrixType;
  typedef _LinearRegressionType LinearRegressionType;
  typedef _ProjectorType ProjectorType;
protected:
  const ProjectorType _projector;
public:
  ConstrainedLevenbergMarquardtAlgorithm (  const int DimRangeF,
                                          const aol::Op<aol::Vector<RealType> > &F,
                                          const aol::Op<aol::Vector<RealType>, MatrixType> &DF,
                                          const ProjectorType &Projector,
                                          const int MaxIterations = 50,
                                          const RealType Mu0 = 1,
                                          const RealType Rho0 = 0.2,
                                          const RealType Rho1 = 0.8,
                                          const RealType EpsDeltaX = 1e-6,
                                          const RealType EpsF = 1e-6,
                                          const RealType EpsGradRealTargetFunc = 1e-6,
                                          const bool Verbose = true )
  : LevenbergMarquardtAlgorithm<RealType, MatrixType, LinearRegressionType> ( DimRangeF, F, DF, MaxIterations, Mu0, Rho0, Rho1, EpsDeltaX, EpsF, EpsGradRealTargetFunc, Verbose ),
  _projector ( Projector ) { }
protected:
  void step ( const aol::Vector<RealType> &X0, const aol::Vector<RealType> &Direction, aol::Vector<RealType> &X1 ) const {
    LevenbergMarquardtAlgorithm<RealType, MatrixType, LinearRegressionType>::step ( X0, Direction, X1 );
    _projector.apply ( X1, X1 );
  }
};


/*
 *  \brief Solves the linear complementary problem:
 *         $w = M * z + q$, s.t. $w,z \geq 0$, $w \cdot z = 0$
 *
 *  Implements the principal pivoting method due to
 *     The principal pivoting method of quadratic programming
 *     RW Cottle - Mathematics of Decision Sciences, Part, 1968
 *
 *  Code was ported from a MATLAB implementation by
 *     Andreas Almqvist, Andrew Spencer and Peter Wall
 *     http://www.mathworks.com/matlabcentral/fileexchange/41485-a-pivoting-algorithm-solving-linear-complementarity-problems
 *
 *  \author mevenkamp
 */
template <typename _RealType, typename _MatrixType = aol::FullMatrix<_RealType> >
class LinearComplementaryProblemSolver {
  typedef _RealType RealType;
  typedef aol::Vector<RealType> VectorType;
  typedef _MatrixType MatrixType;
protected:
  const MatrixType &_M;
  const RealType _pivTol;
  const int _maxIt;
  MatrixType _tableau;
  VectorType _pivot;
public:
  LinearComplementaryProblemSolver ( const MatrixType &M,
                                    const RealType PivTol = 1e-8,
                                    const int MaxIt = 1e4 )
  : _M ( M ), _pivTol ( PivTol ), _maxIt ( MaxIt ) {
    if ( M.getNumCols ( ) != M.getNumRows ( ) ) throw aol::Exception ( "Matrix must be square!", __FILE__, __LINE__ );
    
    const int n = M.getNumRows ( ), m = 2 * ( n + 1 );
    _tableau.reallocate ( n, m );
    _pivot.reallocate ( m );
  }
  
  bool apply ( const VectorType &Q, VectorType &Z ) {
    const int n = _M.getNumCols ( ), m = 2 * ( n + 1 );
    if ( Q.size ( ) != n ) throw aol::Exception ( "Input vector has wrong size!", __FILE__, __LINE__ );
    if ( Z.size ( ) != n ) throw aol::Exception ( "Output vector has wrong size!", __FILE__, __LINE__ );
    
    bool rayTerm = false;
    int loopCount = 0;
    if ( Q.getMinValue ( ) >= 0 ) {
      Z.setZero ( ); // trivial solution
      return true;
    } else {
      // Create initial tableau
      _tableau.setZero ( );
      for ( int i=0; i<n ; ++i ) {
        _tableau.set ( i, i, 1 );
        for ( int j=0; j<n ; ++j )
          _tableau.set ( i, j+n, -_M.get ( i, j ) );
      }
      for ( int i=0; i<n ; ++i )
        _tableau.set ( i, 2*n, -1 );
      for ( int i=0; i<n ; ++i )
        _tableau.set ( i, 2*n+1, Q[i] );
      
      // Let artificial variable enter the basis
      aol::Vector<int> basis ( n );
      for ( int i=0; i<n ; ++i ) basis[i] = i;
      std::pair<int,RealType> indVal = Q.getMinIndexAndValue ( );
      basis[indVal.first] = 2*n;
      int cand = indVal.first + n;
      for ( int j=0; j<m ; ++j ) {
        _pivot[j] = -_tableau.get ( indVal.first, j );
        for ( int i=0; i<n ; ++i )
          _tableau.ref ( i, j ) += _pivot[j];
        _tableau.set ( indVal.first, j, _pivot[j] );
      }
      
      // Perform complementary pivoting
      aol::Vector<RealType> tableauCandCol ( n );
      while ( basis.getMaxValue ( ) == 2 * n && loopCount < _maxIt ) {
        ++loopCount;
        
        // Check if at least one element is not missing
        int numMissing = 0;
        int minQuotIdx = 0;
        RealType minQuot = aol::NumberTrait<RealType>::Inf, quot;
        for ( int i=0; i<n ; ++i ) {
          if ( _tableau.get ( i, cand ) <= 0 ) ++numMissing;
          else {
            quot = _tableau.get ( i, 2*n+1 ) / _tableau.get ( i, cand );
            if ( quot < minQuot ) {
              minQuot = quot;
              minQuotIdx = i;
            }
          }
        }
        
        if ( numMissing != n && aol::Abs<RealType> ( _tableau.get ( minQuotIdx, cand ) ) > _pivTol ) {
          // Reduce tableau
          for ( int i=0; i<n ; ++i ) tableauCandCol[i] = _tableau.get ( i, cand );
          for ( int j=0; j<m ; ++j ) {
            _pivot[j] = _tableau.get ( minQuotIdx, j ) / tableauCandCol[minQuotIdx];
            for ( int i=0; i<n ; ++i ) _tableau.ref ( i, j ) -= tableauCandCol[i] * _pivot[j];
            _tableau.set ( minQuotIdx, j, _pivot[j] );
          }
          int oldVar = basis[minQuotIdx];
          
          // New variable enters the basis
          basis[minQuotIdx] = cand;
          
          // Select next candidate for entering the basis
          if ( oldVar > n ) cand = oldVar - n;
          else cand = oldVar + n;
        } else {
          rayTerm = true;
          break;
          
          std::cerr << "LCP solver: encountered ray termination!" << std::endl;
        }
      }
      
      std::cerr << "LCP solver: terminated after " << loopCount << " loops" << std::endl;
      
      Z.setZero ( );
      aol::Vector<RealType> vars ( 2*n+1 );
      for ( int i=0; i<n ; ++i ) vars[basis[i]] = _tableau.get ( i, 2*n+1 );
      for ( int i=0; i<n ; ++i ) Z[i] = vars[i+n];
      return !rayTerm;
    }
  }
};
  
  
/*
 *  \brief given a noisy signal and its approximate standard deviation, returns the smoothest signal that lies within a standard deviation of the noisy one
 *
 *  for given $y \in \mathbb{R}^n, \sigma > 0$, solves $\min_{x \in \mathbb{R}^n} {\frac{1}{2} \|b - A x\|^2 s.t. y_i - \sigma \leq x_i \leq y_i + \sigma \forall i=1,\dots,n}
 *  where $A=(e_1^T,D)$ and $D$ is a finite difference matrix with an additional $e_1$ in the first row and $b = (y_1,0,\dots,0)$
 *
 *  \author mevenkamp
 */
template <typename RealType>
static void getSmoothMeanVals ( const aol::Vector<RealType> &Input, const RealType StdDev, aol::Vector<RealType> &SmoothMeanVals ) {
  // Assemble matrix for residual minimization (FD + constraint for constant shift) and right-hand side
  // This implementation allows this part to be altered for testing
  const int n = Input.size ( );
  aol::FullMatrix<RealType> A ( n, n );
  A.set ( 0, 0, 1 );
  for ( int i=1; i<n-1 ; ++i ) {
    A.set ( i, i-1, -1 );
    A.set ( i, i, 2 );
    A.set ( i, i+1, -1 );
  }
  A.set ( n-1, n-1, 1 );
  
  const int m = A.getNumRows ( );
  aol::Vector<RealType> b ( m );
  b[0] = Input[0];
  b[n-1] = Input[n-1];
  
  // Assemble inequality constraints and matrices used for the linear complementary problem
  aol::FullMatrix<RealType> G ( 2 * n, n );
  for ( int i=0; i<n ; ++i ) {
    G.set ( i, i, 1 );
    G.set ( i+n, i, -1 );
  }
  
  aol::Vector<RealType> h ( 2 * n );
  for ( int i=0; i<n ; ++i ) {
    h[i] = Input[i] - StdDev;
    h[i+n] = - ( Input[i] + StdDev );
  }
  
  aol::FullMatrix<RealType> AT ( A.getNumCols ( ), A.getNumRows ( ) );
  A.transposeTo ( AT );
  aol::FullMatrix<RealType> ATA ( n, n );
  ATA.makeProduct ( AT, A );
  aol::FullMatrix<RealType> ATAinv ( n, n );
  ATAinv.makeInverse ( ATA );
  
  aol::FullMatrix<RealType> GT ( G.getNumCols ( ), G.getNumRows ( ) );
  G.transposeTo ( GT );
  
  aol::FullMatrix<RealType> T ( n, 2 * n );
  T.makeProduct ( ATAinv, GT );
  
  aol::FullMatrix<RealType> M ( 2 * n, 2 * n );
  M.makeProduct ( G, T );
  
  
  aol::Vector<RealType> ATb ( n );
  AT.apply ( b, ATb );
  aol::Vector<RealType> q ( 2 * n );
  aol::Vector<RealType> ATAinvATb ( n );
  ATAinv.apply ( ATb, ATAinvATb );
  G.apply ( ATAinvATb, q );
  q -= h;
  
  // Solve linear complementary problem: M z + q = w, w,z >= 0, w * z = 0
  aol::Vector<RealType> z ( 2 * n );
  LinearComplementaryProblemSolver<RealType> lcpSolver ( M );
  lcpSolver.apply ( q, z );
  
  // Retrieve smooth mean vals from z
  T.apply ( z, SmoothMeanVals );
  SmoothMeanVals += ATAinvATb;
}


/*
 *  \brief fast version of getSmoothMeanVals with hard-coded matrix and vector assembling for a specific FD matrix
 *
 *  \author mevenkamp
 */
template <typename RealType>
static void fastGetSmoothMeanVals ( const aol::Vector<RealType> &Input, const RealType StdDev, aol::Vector<RealType> &SmoothMeanVals ) {
  const int n = Input.size ( );
  
  // Assemble q
  aol::Vector<RealType> q ( 2*n );
  const RealType a = Input[0], b = Input[n-1];
  RealType ti;
  for ( int i=0; i<n ; ++i ) {
    ti = i / static_cast<RealType> ( n - 1 );
    q[i] = ( 1 - ti ) * a + ti * b;
    q[i+n] = -q[i];
    q[i] -= Input[i] - StdDev;
    q[i+n] -= - ( Input[i] + StdDev );
  }
  
  // Assemble M
  aol::FullMatrix<RealType> M ( 2*n, 2*n );
  RealType mij;
  int64_t s11, s12, s13, s1, s21, s22, s23, s2, s31, s32, s33, s34, s35, s3, sum;
  const RealType nm1Sqr = static_cast<RealType> ( ( n-1 ) * ( n-1 ) );
  for ( int i=1; i<=n ; ++i ) {
    for ( int j=i; j<=n ; ++j ) {
      s11 = (n-i)*(n-j);
      s12 = (i-2)*(i-1);
      s13 = 2*i-3;
      s1 = s11 * ( 6 + s12 * s13 );
      
      s21 = (i-1)*(j-1);
      s22 = (n-j)*(n-j+1);
      s23 = 2*(n-j)+1;
      s2 = s21 * ( 6 + s22 * s23 );
      
      s31 = (i-1)*(n-j);
      s32 = (i-j);
      s33 = 2*(aol::Sqr<int64_t>(i+j)-i*j+2);
      s34 = 9*n;
      s35 = -3*(2+n)*(i+j);
      s3 = s31 * s32 * ( s33 + s34 + s35 );
  
      sum = s1 + s2 + s3;
      mij = sum / nm1Sqr / 6.0;

      M.set ( i-1, j-1, mij );
      M.set ( j-1, i-1, mij );
    }
  }
  for ( int i=0; i<n ; ++i ) {
    for ( int j=0; j<n ; ++j ) {
      mij = M.get ( i, j );
      M.set ( n+i, n+j, mij );
      M.set ( n+i, j, -mij );
      M.set ( i, n+j, -mij );
    }
  }

  // Solve linear complementary problem: M z + q = w, w,z >= 0, w * z = 0
  aol::Vector<RealType> z ( 2*n );
  LinearComplementaryProblemSolver<RealType> lcpSolver ( M );
  lcpSolver.apply ( q, z );
  
  // Fast multiplication of (ATA)^-1 * G^T with z
  SmoothMeanVals.setZero ( );
  for ( int j=0; j<2*n ; ++j ) {
    if ( z[j] != 0 ) {
      for ( int i=0; i<n ; ++i )
        SmoothMeanVals[i] += z[j] * M.get ( i, j );
    }
  }
  
  // Fast addition of (A^TA)^-1 A^T b to x
  for ( int i=0; i<n ; ++i ) {
    ti = i / static_cast<RealType> ( n - 1 );
    SmoothMeanVals[i] += ( 1 - ti ) * a + ti * b;
  }
}
  

/*
 * \brief fits data $(x_i,y_i)$ with a linear model $y(x) = P_0 * x + P_1$
 *
 * \author mevenkamp
 */
template <typename RealType>
static void getLinearFit ( const std::vector<std::pair<RealType, RealType> > &Data, aol::Vec2<RealType> &Params ) {
  const int N = static_cast<int> ( Data.size ( ) );
  RealType sumXSqr = 0, sumY = 0, sumX = 0, sumXY = 0, x, y;
  for ( int i=0; i<N ; ++i ) {
    x = Data[i].first;
    y = Data[i].second;
    sumXSqr += aol::Sqr<RealType> ( x );
    sumY += y;
    sumX += x;
    sumXY += x * y;
  }
  Params[0] = ( N * sumXY - sumX * sumY ) / ( N * sumXSqr - aol::Sqr<RealType> ( sumX ) );
  Params[1] = ( sumXSqr * sumY - sumX * sumXY ) / ( N * sumXSqr - aol::Sqr<RealType> ( sumX ) );
}
  
  
/*
 *  \brief fits data $(x_i,y_i)$ with a model of the form $y(x) = P_0 * \exp{P_1 * x}$
 * 
 *  A linear model through log transformation that minimizes a weighted least squares functional is used, according to:
 *  http://mathworld.wolfram.com/LeastSquaresFittingExponential.html
 *  Positive input $x_i,y_i > 0$ is assumed. This is enforced by clipping values smaller than 1e-8
 *
 *  \author mevenkamp
 */
template <typename RealType>
static void getExponentialFit ( const std::vector<std::pair<RealType, RealType> > &Data, aol::Vec2<RealType> &Params ) {
  RealType sumXSqrY = 0, sumYlnY = 0, sumXY = 0, sumXYlnY = 0, sumY = 0, x, y;
  for ( int i=0; i<Data.size ( ) ; ++i ) {
    x = aol::Max ( 1e-8, Data[i].first );
    y = aol::Max ( 1e-8, Data[i].second );
    sumXSqrY += aol::Sqr<RealType> ( x ) * y;
    sumYlnY += y * log ( y );
    sumXY += x * y;
    sumXYlnY += x * y * log ( y );
    sumY += y;
  }
  Params[0] = exp ( ( sumXSqrY * sumYlnY - sumXY * sumXYlnY ) / ( sumY * sumXSqrY - aol::Sqr<RealType> ( sumXY ) ) );
  Params[1] = ( sumY * sumXYlnY - sumXY * sumYlnY ) / ( sumY * sumXSqrY - aol::Sqr<RealType> ( sumXY ) );
}

/*
 *  \brief fits data $(x_i,y_i)$ with a model of the form $y(x) = P_0 * x^{P_1}$
 *
 *  Algorithm is based on: http://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html
 *  Positive input $x_i,y_i > 0$ is assumed. This is enforced by clipping values smaller than 1e-8
 *
 *  \author mevenkamp
 */
template <typename RealType>
static void getPowerLawFit ( const std::vector<std::pair<RealType, RealType> > &Data, aol::Vec2<RealType> &Params ) {
  RealType sumlnXlnY = 0, sumlnX = 0, sumlnY = 0, sumlnXSqr = 0, x, y;
  for ( unsigned int i=0; i<Data.size ( ) ; ++i ) {
    x = aol::Max<RealType> ( 1e-8, Data[i].first );
    y = aol::Max<RealType> ( 1e-8, Data[i].second );
    sumlnXlnY += log ( x ) * log ( y );
    sumlnX += log ( x );
    sumlnY += log ( y );
    sumlnXSqr += aol::Sqr<RealType> ( log ( x ) );
  }
  
  const RealType n = static_cast<RealType> ( Data.size ( ) );
  Params[1] = ( n * sumlnXlnY - sumlnX * sumlnY ) / ( n * sumlnXSqr - aol::Sqr<RealType> ( sumlnX ) );
  Params[0] = exp ( ( sumlnY - Params[1] * sumlnX ) / n );
}

/**
 * \brief Helper class for aol::GaussNewtonAlgorithm to solve the linear least squares problem that occurs in each Gauss-Newton step.
 *
 * \author Berkels
 */
template <typename RealType>
class GaussNewtonDenseQRSolver {
#ifdef USE_EXTERNAL_EIGEN
  eig::FullMatrix<RealType> A;
  eig::Vector<RealType> rhs;
  eig::Vector<RealType> sol;
#endif
  aol::Vector<RealType> f;
  aol::Vector<RealType> direction;

public:
  GaussNewtonDenseQRSolver ( const int NumRows, const int NumCols ) :
#ifdef USE_EXTERNAL_EIGEN
    A ( NumRows, NumCols ),
    rhs ( NumRows ),
    sol ( NumCols ),
    f ( rhs.data(), NumRows, aol::FLAT_COPY ),
    direction ( sol.data(), NumCols, aol::FLAT_COPY )
#else
    f ( NumRows ),
    direction ( NumCols )
#endif
  {
#ifndef USE_EXTERNAL_EIGEN
    cerr << "aol::QRInversePivot is an order of magnitude slower than the algorithms in the Eigen librarby. It is highly recommended to use the Eigen external for this class.\n";
#endif
  }

  aol::Vector<RealType> &getFReference ( ) {
    return f;
  }

  aol::Vector<RealType> &getDirectionReference ( ) {
    return direction;
  }

  template <typename MatrixType>
  void solve ( const MatrixType &DF ) {
#ifdef USE_EXTERNAL_EIGEN
    // The performance of colPivHouseholderQr is best if the matrix uses column-major storage.
    // The easiest way to achieve this is to simply copy the data from our matrix structure to
    // an Eigen matrix. Otherwise _DF needs to be implemented in a way that works with
    // column-major storage, but currently BumpFitGaussNewtonTargetJacobian is written with
    // row-major storage in mind.
    A.copyFrom ( DF );
    sol = A.colPivHouseholderQr().solve(rhs);
#else
    aol::QRInversePivot<RealType> qrInverse ( DF );
    qrInverse.apply ( f, direction );
#endif
  }
};

/**
 * \author Berkels
 */
template <typename RealType>
class GaussNewtonSparseEigenBlockSolverBase {
protected:
#if defined(USE_EXTERNAL_EIGEN)
  eig::SparseMatrix<RealType> A;
  eig::Vector<RealType> rhs;
  eig::Vector<RealType> sol;
#endif
  aol::MultiVector<RealType> f;
  aol::MultiVector<RealType> direction;
public:
  GaussNewtonSparseEigenBlockSolverBase ( const aol::Vector<int> &NumRows, const aol::Vector<int> &NumCols )
#if defined(USE_EXTERNAL_EIGEN)
    : A ( NumRows.sum(), NumCols.sum() ),
      rhs ( NumRows.sum() ),
      sol ( NumCols.sum() ),
      f ( rhs.data(), NumRows, aol::FLAT_COPY ),
      direction ( sol.data(), NumCols, aol::FLAT_COPY )
#endif
  { }

  aol::MultiVector<RealType> &getFReference ( ) {
    return f;
  }

  aol::MultiVector<RealType> &getDirectionReference ( ) {
    return direction;
  }
};

#if defined(USE_EXTERNAL_EIGEN) && defined(USE_EXTERNAL_SUITESPARSE)

/**
 * \author Berkels
 */
template <typename RealType, typename EigenQRSolverType = Eigen::SPQR<Eigen::SparseMatrix<RealType> > >
class GaussNewtonSparseQRBlockSolver : public GaussNewtonSparseEigenBlockSolverBase<RealType> {
public:
  GaussNewtonSparseQRBlockSolver ( const aol::Vector<int> &NumRows, const aol::Vector<int> &NumCols )
    : GaussNewtonSparseEigenBlockSolverBase<RealType> ( NumRows, NumCols ) { }

  template <typename MatrixType>
  void solve ( const MatrixType &DF ) {
    this->A.copyFrom ( DF );
    EigenQRSolverType solver;
    solver.compute ( this->A );
    if ( solver.info() != Eigen::Success ) {
      cerr << "decomposition failed\n";
      return;
    }
    this->sol = solver.solve ( this->rhs );
    if ( solver.info() != Eigen::Success ) {
      cerr << "solving failed\n";
      return;
    }
  }
};

/**
 * \author Berkels
 */
template <typename RealType, typename EigenSolverType = Eigen::CholmodDecomposition<Eigen::SparseMatrix<RealType> > >
class GaussNewtonSparseNormalEquationsBlockSolver : public GaussNewtonSparseEigenBlockSolverBase<RealType> {
public:
  GaussNewtonSparseNormalEquationsBlockSolver ( const aol::Vector<int> &NumRows, const aol::Vector<int> &NumCols )
    : GaussNewtonSparseEigenBlockSolverBase<RealType> ( NumRows, NumCols ) { }

  template <typename MatrixType>
  void solve ( const MatrixType &DF ) {
    this->A.copyFrom ( DF );
    EigenSolverType solver;
    solver.compute ( this->A.transpose() * this->A );
    if ( solver.info() != Eigen::Success ) {
      cerr << "decomposition failed\n";
      return;
    }
    eig::Vector<RealType> ATrhs;
    ATrhs = this->A.transpose() * this->rhs;
    this->sol = solver.solve ( ATrhs );
    if ( solver.info() != Eigen::Success ) {
      cerr << "solving failed\n";
      return;
    }
  }
};
#else

/**
 * \author Berkels
 */
template <typename RealType>
class GaussNewtonSparseNormalEquationsBlockSolver : public GaussNewtonSparseEigenBlockSolverBase<RealType> {
public:
  GaussNewtonSparseNormalEquationsBlockSolver ( const aol::Vector<int> &NumRows, const aol::Vector<int> &NumCols )
    : GaussNewtonSparseEigenBlockSolverBase<RealType> ( NumRows, NumCols ) { }

  template <typename MatrixType>
  void solve ( const MatrixType &/*DF*/ ) {
    throw aol::Exception ( "aol::GaussNewtonSparseNormalEquationsBlockSolver needs eigen and suitesparse! Compile with -USE_EXTERNAL_EIGEN and -DUSE_SUITESPARSE", __FILE__, __LINE__ );
  }
};

#endif

/**
 * \brief <a href="https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm">Gauss-Newton algorithm</a> to solve non-linear
 *        least squares problems of type \f$\Vert F(x)\Vert_2^2\rightarrow \min\f$, where \f$ F : \mathbb{R}^n \rightarrow \mathbb{R}^m\f$.
 *
 * \author Berkels
 * \ingroup Optimization
 */
template <typename VectorType, typename MatrixType, typename LinLSSolverType = GaussNewtonDenseQRSolver<typename VectorType::RealType> >
class GaussNewtonAlgorithm : public ArmijoLineSearchUsingOp<typename VectorType::RealType, VectorType > {
  typedef typename VectorType::RealType RealType;
  typedef typename VectorInitTrait<VectorType>::InitType VectorInitTypeType;

  //! stores the evaluation of DF and is deleted in the destructor.
  mutable DeleteFlagPointer<MatrixType> _pMatDF;

  mutable RealType _fNormSqrAtLastPosition;

  const VectorInitTypeType _dimRangeF;
  const VectorInitTypeType _dimDomF;
  const aol::Op<VectorType> &_F;
  const aol::Op<VectorType, MatrixType> &_DF;
  const int _maxIterations;
  const RealType _stopEpsilon;

  const StepSaverBase<RealType, VectorType> *_pStepSaver;

  RealType ArmijoLineSearchHelpFunction_evaluate( const VectorType &DescentDir, const VectorType &CurrentPosition, const RealType timestepWidth ) const {
    VectorType temp ( CurrentPosition );
    temp.addMultiple ( DescentDir, timestepWidth );
    VectorType f ( _dimRangeF );
    _F.apply ( temp, f );
    return f.normSqr();
  }
  
  // Attention: ArmijoLineSearchHelpFunction_evaluateDerivative only works, if _pMatDF is filled with _DF(Position)!
  RealType ArmijoLineSearchHelpFunction_evaluateDerivative( const VectorType &DescentDir, const VectorType &Position ) const {
    VectorType tmp ( _dimRangeF );
    VectorType tmp2 ( _dimRangeF );
    _pMatDF->apply ( DescentDir, tmp );
    _F.apply ( Position, tmp2 );
    return ( 2*(tmp*tmp2) );
  }

public:
  GaussNewtonAlgorithm ( const VectorInitTypeType &DimRangeF,
                         const VectorInitTypeType &DimDomF,
                         const aol::Op<VectorType> &F,
                         const aol::Op<VectorType, MatrixType> &DF,
                         const int MaxIterations = 50,
                         const RealType StopEpsilon = 0 )
    : ArmijoLineSearchUsingOp<RealType, VectorType> ( 0.5, MaxIterations ),
      _pMatDF ( NULL ),
      _dimRangeF ( DimRangeF ),
      _dimDomF ( DimDomF ),
      _F ( F ),
      _DF ( DF ),
      _maxIterations ( MaxIterations ),
      _stopEpsilon ( StopEpsilon ),
      _pStepSaver ( NULL ) { }

  void apply( const VectorType &Arg, VectorType &Dest ) const {
    Dest = Arg;
    applySingle ( Dest );
  }

  void applySingle( VectorType &Dest ) const {
    if ( _pMatDF.get() == NULL )
      _pMatDF.reset ( new MatrixType ( _dimRangeF, _dimDomF ), true );

    int iteration = 0;
    RealType fNormSqrOld = aol::NumberTrait<RealType>::Inf;
    RealType fNormSqr = aol::NumberTrait<RealType>::Inf;

    LinLSSolverType solver ( _dimRangeF, _dimDomF );
    VectorType &f = solver.getFReference();
    VectorType &direction = solver.getDirectionReference();

    _F.apply ( Dest, f );
    fNormSqrOld = f.normSqr();

    cerr << "Initial fNormSqr " << fNormSqrOld << endl;

    RealType tau = 1;

    do {
      _DF.apply ( Dest, *_pMatDF );

      solver.solve ( *_pMatDF );

      if ( direction.checkForNANsAndINFs() ) {
        cerr << "Errr: LinLSSolverType failed.\n";
        break;
      }

      Dest -= direction;

      _F.apply ( Dest, f );
      fNormSqr = f.normSqr();

      // If the target functional did not decrease with the update, try to find a smaller step
      // so that it does. This step size control is extremely simple and not very efficient, but
      // it's certainly better than letting the algorithm diverge.
      if ( fNormSqr >= fNormSqrOld ) {
        Dest += direction;
        direction *= -1;
        // getTimestepWidthWithSimpleLineSearch doesn't support "widening", so let it start with 2*tau.
        tau = this->getTimestepWidthWithSimpleLineSearch ( direction, Dest, aol::Min ( 2 * tau, aol::ZOTrait<RealType>::one) );
        Dest.addMultiple ( direction, tau );
        _F.apply ( Dest, f );
        fNormSqr = f.normSqr();
      }
      else
        tau = 1;

      cerr << "Step " << iteration + 1 << ", tau " << tau << ", fNormSqr " << fNormSqr << ", diff " << fNormSqrOld - fNormSqr << endl;

      const bool stop = ( ( fNormSqrOld - fNormSqr ) ) <= _stopEpsilon * fNormSqr;
      fNormSqrOld = fNormSqr;
      ++iteration;

      if ( _pStepSaver )
        _pStepSaver->saveStep ( Dest, iteration );

      if ( stop )
        break;
    } while ( ( iteration < _maxIterations ) && ( appeqAbsolute ( fNormSqr, aol::ZOTrait<RealType>::zero ) == false ) );
    
    _fNormSqrAtLastPosition = fNormSqrOld;
  }

  void applyAdd( const VectorType &/*Arg*/, VectorType &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__);
  }

  void setStepSaverReference ( const StepSaverBase<RealType, VectorType> &StepSaver ) {
    _pStepSaver = &StepSaver;
  }

  RealType getFNormSqrAtLastPosition ( ) const {
    return _fNormSqrAtLastPosition;
  }
};

/**
 * \author Berkels
 * \ingroup Optimization
 */
template <typename RealType, typename MatrixType = aol::SparseMatrix<RealType> >
class LinearGNFuncBase : public aol::Op<aol::Vector<RealType> > {
protected:
  MatrixType _mat;
public:
  LinearGNFuncBase ( const int DimRangeF, const int DimDomF )
    : _mat ( DimRangeF, DimDomF ) { }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    _mat.apply ( Arg, Dest );
  }

  void applyAdd ( const aol::Vector<RealType> &, aol::Vector<RealType> & ) const {
    throw aol::UnimplementedCodeException ( "applyAdd not implemented...", __FILE__, __LINE__ );
  }

  void applyDerivative ( const aol::Vector<RealType> &/*Arg*/, MatrixType &MatDest ) const {
    MatDest = _mat;
  }

  int getDimRangeF ( ) const {
    return _mat.getNumRows();
  }

  const MatrixType &getMatrixRef ( ) const {
    return _mat;
  }
};

} // end namespace
  

#endif
