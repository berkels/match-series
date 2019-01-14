#ifndef __CONVERGENCEESTIMATOR_H
#define __CONVERGENCEESTIMATOR_H

#include <ringBufferFrontInsertCapped.h>

namespace aol {

/****************************************************************************
 *
 *        CLASS ConvergenceEstimator
 */
/**
 *  \brief Administrates the storage of previous iterations' residuals
 *         and computation of convergence coefficents and order
 *         estimates.
 *
 *  Use this class when computing entries of a sequence that
 *  should tend to zero with a certain convergence rate and
 *  order. It takes care of storing the last residuals and
 *  computes estimates for a given number of convergence
 *  coefficients.
 *
 *  That means: Suppose you feed an instance of this class
 *  with values \f$ x_i \f$. Convergence of order \f$ t \f$
 *  then means that \f$ x_i = C^{(t)} x_{i-1}^t \f$ should hold.
 *  Thus, for given \f$ t_{\textrm{max}} \f$, the class computes
 *  \f$ C^{(t)}_i = x_i / x_{i-1}^t \f$ for all \f$ t=1, \dots, t_{\textrm{max}} \f$.
 *  Additionally, a (real valued) convergence order estimation
 *  \f$ \kappa_i \f$ is computed such that \f$ (x_i / x_{i-1})^{\kappa_i}
 *  = (x_{i+1} / x_i) \f$.
 *
 *  Note that \f$ C^{(0)} \f$ contains the passed residuals.
 *
 *  You can request either (a) these quantities computed only from
 *  the last and before-last step (and the third-last step for
 *  \f$ \kappa_i \f$) by calling getLast...() or (b) values that are
 *  averaged over a certain number of steps (given in constructor).
 *
 *  \author von Deylen
 *          Created 2008-04-03
 *          Last changed 2008-04-04
 */

template <typename RealType>
class ConvergenceEstimator {
public:
  ConvergenceEstimator ( int numAveragingSteps = 5, int numConvergenceCoeff = 1 );

  void clear();

  //! use push_back() to exprime that residual belongs to
  //! a new step.
  void push ( RealType residual );

  //! use setResidualAgain() to reset the residual for the current
  //! step when push_back() has already been called.
  void setResidualAgain ( RealType residual );

  RealType getLastConvergenceCoeff ( int order ) const;
  RealType getLastConvergenceOrder() const;
  RealType getConvergenceCoeff ( int order ) const;
  RealType getConvergenceOrder() const;

  RealType getPotentialConvergenceCoeff() const;

  int getNumComputedSteps() const;
  int getNumConvergenceCoeff() const;

  int getNumAveragingSteps() const;

protected:
  //! computes convergence order and coefficients
  //! depending on the residuals stored in
  //! _convergenceCoeff[0].
  void computeStep();

  RingBufferFrontInsertCapped<RealType> _convergenceOrder;
  vector<RingBufferFrontInsertCapped<RealType> > _convergenceCoeff;

};

/**** IMPLEMENTATION ********************************************************
 *
 *
 *
 */

//---------------------------------------------------------------------------
template <typename RealType>
ConvergenceEstimator<RealType>::
ConvergenceEstimator ( int numAveragingSteps, int numConvergenceCoeff )
    : _convergenceOrder ( numAveragingSteps )
    , _convergenceCoeff ( numConvergenceCoeff + 1, _convergenceOrder ) {}
//---------------------------------------------------------------------------
template <typename RealType>
void ConvergenceEstimator<RealType>::clear() {
  _convergenceOrder.clear();
  for ( size_t i = 0; i < _convergenceCoeff.size(); ++i )
    _convergenceCoeff[i].clear();
}
//---------------------------------------------------------------------------
template <typename RealType>
void ConvergenceEstimator<RealType>::push ( RealType residual ) {
  _convergenceCoeff[0].push ( residual );
  computeStep();
}
//---------------------------------------------------------------------------
template <typename RealType>
void ConvergenceEstimator<RealType>::setResidualAgain ( RealType residual ) {
  _convergenceCoeff[0][0] = residual;
  computeStep();
}
//---------------------------------------------------------------------------

template <typename RealType>
void ConvergenceEstimator<RealType>::computeStep() {
// if this is at least the second step, compute
  // convergence coefficients.
  if ( getNumComputedSteps() > 1 ) {
    RealType coeff = _convergenceCoeff[0][0];
    RealType oldResidual = _convergenceCoeff[0][1];
    for ( size_t i = 1; i < _convergenceCoeff.size(); ++i ) {
      coeff /= oldResidual;
      _convergenceCoeff[i].push ( coeff );
    }
  }

  // if this is at least the third step, compute
  // convergence order, too.
  if ( getNumComputedSteps() > 2 ) {
    RealType quot = _convergenceCoeff[0][0] / _convergenceCoeff[0][1];
    RealType oldQuot = _convergenceCoeff[0][1] / _convergenceCoeff[0][2];
    RealType order = log ( quot ) / log ( oldQuot );
    _convergenceOrder.push ( order );
  }
}

//---------------------------------------------------------------------------

template <typename RealType>
RealType ConvergenceEstimator<RealType>::
getPotentialConvergenceCoeff() const {
  RealType x_i_pot = pow ( _convergenceCoeff[0][1], _convergenceOrder[0] );
  return _convergenceCoeff[0][0] / x_i_pot;
}

//---------------------------------------------------------------------------
template <typename RealType>
RealType ConvergenceEstimator<RealType>::
getLastConvergenceCoeff ( int order ) const
                                  {   return _convergenceCoeff[order][0];   }
//---------------------------------------------------------------------------
template <typename RealType>
RealType ConvergenceEstimator<RealType>::
getLastConvergenceOrder() const          {   return _convergenceOrder[0];   }
//---------------------------------------------------------------------------
template <typename RealType>
RealType ConvergenceEstimator<RealType>::
getConvergenceCoeff ( int order ) const
                    {   return _convergenceCoeff[order].arithmeticMean();   }
//---------------------------------------------------------------------------
template <typename RealType>
RealType ConvergenceEstimator<RealType>::
getConvergenceOrder() const {  return _convergenceOrder.arithmeticMean();   }
//---------------------------------------------------------------------------
template <typename RealType>
int ConvergenceEstimator<RealType>::
getNumComputedSteps() const       {   return _convergenceCoeff[0].size();   }
//---------------------------------------------------------------------------
template <typename RealType>
int ConvergenceEstimator<RealType>::
getNumConvergenceCoeff() const
                {   return static_cast<int> ( _convergenceCoeff.size() );   }
//---------------------------------------------------------------------------
template <typename RealType>
int ConvergenceEstimator<RealType>::
getNumAveragingSteps() const      {   return _convergenceOrder.maxSize();   }
//---------------------------------------------------------------------------

} // end of namespace aol.

#endif
