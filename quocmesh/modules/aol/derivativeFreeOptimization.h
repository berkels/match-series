#ifndef __DERIVATIVEFREEOPTIMIZATION_H
#define __DERIVATIVEFREEOPTIMIZATION_H

#include <vec.h>
#include <op.h>
#include <matrix.h>

namespace aol {

//-----------------------------------------------------------------------------------------------

template <typename RealType>
class FunctorFromOp {

public:
  FunctorFromOp ( const Op<RealType, Scalar<RealType> > & Operator )
    : _operator ( Operator )
  {}

  RealType operator () ( RealType Arg ) const {
    Scalar<RealType> Dest;
    _operator.apply ( Arg, Dest );
    return Dest;
  }

protected:
  const Op<RealType, Scalar<RealType> > & _operator;
};

//-----------------------------------------------------------------------------------------------

template <typename RealType, typename FunctorType>
RealType findOptimumByBisection ( const FunctorType & scalarOp,
                                  RealType posMin,
                                  RealType posMax,
                                  RealType TOL = 1E-4 ) {

  RealType en_upper, en_lower;

  RealType delta = posMax - posMin;
  RealType g = (3. - sqrt(5.)) / 2.;
  RealType pos_upper = posMax;
  RealType pos_lower = posMin;

  en_upper = scalarOp ( pos_upper );
  en_lower = scalarOp ( pos_lower );

  string tab = "\t";

  clog << "lower bd" << tab << "evaluation" << tab << tab << "upper bd" << tab << "evaluation" << endl
       << pos_lower << tab << en_lower << tab << tab << pos_upper << tab << en_upper << endl;

  clog.precision ( 5 );

  while ( delta > TOL ) {
    RealType pos_middle = pos_lower + g * delta;

    if ( en_lower < en_upper ) {
      pos_upper = pos_middle;
      en_upper = scalarOp ( pos_upper );
    }
    else {
      pos_lower = pos_middle;
      en_lower = scalarOp ( pos_lower );
    }
    clog << pos_lower << tab << en_lower << tab << tab << pos_upper << tab << en_upper << endl;
    delta = pos_upper - pos_lower;
  }
  RealType pos_middle = pos_lower + g * delta;
  clog << endl << "Bracketed pos in an interval of length " << delta << "." << endl
       << "Will return " << pos_middle << " as optimal position." << endl;

  return pos_middle;
}

//-----------------------------------------------------------------------------------------------

template <typename RealType, typename FunctorType>
RealType findRootByBisection ( const FunctorType & scalarOp,
                               RealType posMin,
                               RealType posMax,
                               RealType TOL = 1E-4 ) {

  aol::Scalar<RealType> en_upper, en_middle, en_lower;

  RealType delta = posMax - posMin;
  RealType g = (3. - sqrt(5.)) / 2.;
  RealType pos_upper = posMax;
  RealType pos_lower = posMin;

  en_upper = scalarOp ( pos_upper );
  en_lower = scalarOp ( pos_lower );

  if ( en_upper * en_lower > 0. ) {
    cerr << "Error: evaluation of given operator has equal sign" << endl
         << "at position " << posMin << " and " << posMax << ". Will return 1." << endl;
    return 1.;
  }

  if ( en_lower > 0. ) {
    swap ( pos_lower, pos_upper );
    swap ( en_lower, en_upper );
  }

  string tab = "\t";

  clog << "neg bd" << tab << "evaluation" << tab << tab << "pos bd" << tab << "evaluation" << endl
       << pos_lower << tab << en_lower << tab << tab << pos_upper << tab << en_upper << endl;

  clog.precision ( 5 );

  RealType pos_middle = pos_lower + g * delta;
  en_middle = scalarOp ( pos_middle );

  while ( delta > TOL ) {
    if ( en_middle > 0. ) {
      pos_upper = pos_middle;
      en_upper = en_middle;
    }
    else {
      pos_lower = pos_middle;
      en_lower = en_middle;
    }
    clog << pos_lower << tab << en_lower << tab << tab << pos_upper << tab << en_upper << endl;

    delta = pos_upper - pos_lower;
    pos_middle = pos_lower + g * delta;
    en_middle = scalarOp ( pos_middle );
  }
  clog << endl << "Bracketed pos in an interval of length " << delta << "." << endl
       << "Will return " << pos_middle << " as optimal position." << endl;

  return pos_middle;
}

//-----------------------------------------------------------------------------------------------

//! Compute derivative of a scalar function by central difference quotients
//! multivector arg contains position and direction,
//! used by SimpleDescent and subclasses to approximate directional derivative
template <class DataType = double, class IndexType = int> class CDQEnergyDerivativeOp
      : public aol::Op<aol::MultiVector<DataType>, aol::Scalar<DataType> > {
public:
  CDQEnergyDerivativeOp ( const aol::Op<aol::Vector<DataType>, aol::Scalar<DataType> > & Energy, DataType h = 1E-8 )
      : _energy ( Energy ), _h ( h ) {}
  void applyAdd ( const aol::MultiVector<DataType>& arg, aol::Scalar<DataType>& dest ) const {
    if ( arg.numComponents () != 2 )
      throw ( aol::DimensionMismatchException ( "CDQEnergyDerivativeOp::applyAdd, arg must contain point and direction", __FILE__, __LINE__ ) );
    aol::Vector<DataType> xph ( arg [0] ), xmh ( arg [0] ), hd ( arg [1] );
    hd *= _h; xph += hd; xmh -= hd;
    dest += ( energy ( xph ) - energy ( xmh ) ) / ( 2 * _h );
  }
  void seth ( DataType h ) const {
    _h = h;
  }
private:
  DataType energy ( const aol::Vector<DataType>& arg ) const {
    aol::Scalar<DataType> temp;
    _energy.apply ( arg, temp );
    return temp;
  }
  const aol::Op<aol::Vector<DataType>, aol::Scalar<DataType> > & _energy;
  mutable DataType _h;
};

//! Compute derivative of a scalar function by forward difference quotients
//! multivector arg contains position and direction,
//! used by SimpleDescent and subclasses to compute slope in Armijo rule
template <class DataType = double, class IndexType = int> class DQEnergyDerivativeOp
      : public aol::Op<aol::MultiVector<DataType>, aol::Scalar<DataType> > {
public:
  DQEnergyDerivativeOp ( const aol::Op<aol::Vector<DataType>, aol::Scalar<DataType> > & Energy, DataType h = 1E-8 )
      : _energy ( Energy ), _h ( h ) {}
  void applyAdd ( const aol::MultiVector<DataType>& arg, aol::Scalar<DataType>& dest ) const {
    if ( arg.numComponents () != 2 )
      throw ( aol::DimensionMismatchException ( "DQEnergyDerivativeOp::applyAdd, arg must contain point and direction", __FILE__, __LINE__ ) );
    aol::Vector<DataType> x ( arg [0] ), xph ( arg [1] );
    xph *= _h; xph += x;
    dest += ( energy ( xph ) - energy ( x ) ) / _h;
  }
  void seth ( DataType h ) const {
    _h = h;
  }
private:
  DataType energy ( const aol::Vector<DataType>& arg ) const {
    aol::Scalar<DataType> temp;
    _energy.apply ( arg, temp );
    return temp;
  }
  const aol::Op<aol::Vector<DataType>, aol::Scalar<DataType> > & _energy;
  mutable DataType _h;
};

//! Compute gradient of a scalar function by central difference quotients
//! apply maps position to gradient direction,
//! used by SimpleDescent and subclasses to approximate gradient
template <class DataType = double, class IndexType = int> class CDQEnergyGradientOp
      : public aol::Op<aol::Vector<DataType>, aol::Vector<DataType> > {
public:
  CDQEnergyGradientOp ( const aol::Op<aol::Vector<DataType>, aol::Scalar<DataType> > & Energy,
                       const aol::Vector<IndexType>& dirs, DataType h = 1E-8 )
      : _energy ( Energy ), _dirs ( dirs ), _h ( h ) {}
  void applyAdd ( const aol::Vector<DataType>& arg, aol::Vector<DataType>& dest ) const {

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for ( IndexType i = 0; i < _dirs.size (); ++i ) {
      aol::Vector<DataType> point ( arg );
      point [_dirs [i]] -= _h;
      DataType em = energy ( point );
      point [_dirs [i]] += 2 * _h;
      DataType ep = energy ( point );
      dest [_dirs [i]] += ( ep - em ) / ( 2 * _h );
    }
  }
  void seth ( DataType h ) const {
    _h = h;
  }
private:
  DataType energy ( const aol::Vector<DataType>& arg ) const {
    aol::Scalar<DataType> temp;
    _energy.apply ( arg, temp );
    return temp;
  }
  const aol::Op<aol::Vector<DataType>, aol::Scalar<DataType> > & _energy;
  const aol::Vector<int>& _dirs;
  mutable DataType _h;
};

  /** \brief Generic l2 energy descent with conjugate directions,
   *         uses difference quotients by default
   *  \author Lenz
   *  \ingroup Optimization
   */
template <class DataType, class IndexType>
class SimpleDescent {
public:
  /** Constructor.
   *  \param e     objective operator to be minimized
   *  \param optimizeExcept directions to be left out during minimization
   *  \param eps   lower bound on gradient (=direction) size
   *  \param h     distance for finite difference
   *  \param hMin  lower bound for h
   *  \param hDec  factor by which h is multiplied
   *  \param beta  for Armijo rule:
   *  \param sigma for Armijo rule:
   *  \param start initial step size = direction norm * start
   *  \param hfac  minimal Armijo step size must be hfac times h
   *  \param reset reset conjugate directions every reset steps
   *  \param maxit maximum number of iterations
   */
  SimpleDescent ( const aol::Op<aol::Vector<DataType>, aol::Scalar<DataType> > & e,
            const aol::Vector<IndexType>& optimizeExcept,
            DataType eps = 1E-3, DataType h = 1E-4, DataType hMin = 1E-6, DataType hDec = 0.5, DataType beta = 0.5,
            DataType sigma = 0.2, DataType start = 0.01, DataType hfac = 2, int reset = 50, int maxit = 1000 )
      : _e ( &e ), _optvars ( 0 ), _cge ( NULL ), _de ( e, h ), _cde ( e, h ), _eps ( eps ), _h ( h ), _hMin ( hMin ), _hDec ( hDec ),
      _beta ( beta ), _sigma ( sigma ), _start ( start ), _hfac ( hfac ), _reset ( reset ), _maxit ( maxit ) {
    int p = 0;
    for ( int i = 0; i < optimizeExcept.size (); ++i )
      if ( optimizeExcept [i] == 0 ) {
        _optvars.resize ( p + 1 );
        _optvars [p] = i;
        ++p;
      }

    _cge = new CDQEnergyGradientOp<DataType, IndexType> ( *_e, _optvars, _h );
  }
  void adjust ( const aol::Vector<DataType>& point, const aol::Vector<DataType>& dir, DataType& step ) const {
    aol::MultiVector<DataType> arg ( 0, 1 );
    arg.appendReference ( point ); arg.appendReference ( dir );

    aol::Scalar<DataType> m_tan;
    aol::Scalar<DataType> m_sec;

    DataType norm = dir.norm ();
    if ( !step ) {
      if ( norm ) step = _start / norm;
      else step = _start;
    }

    _cde.apply ( arg, m_tan );

    while ( step * norm > _start ) step *= _beta;
    do {
      step /= _beta;
      _de.seth ( step );
      _de.apply ( arg, m_sec );
    } while ( m_sec / m_tan > _sigma && step * norm < _start );
    do {
      step *= _beta;
      _de.seth ( step );
      _de.apply ( arg, m_sec );
    } while ( m_sec / m_tan < _sigma && step * norm > _hfac * _h );
  }
  void solve ( aol::Vector<DataType>& point ) const {
    aol::Vector<DataType> dir ( point.size () ), grad ( point.size () ), old_point ( point.size () );
    aol::Scalar<DataType> val_old, val_new;
    DataType normS_old, beta_old = 1, step = 0, realstep;
    int i = 1;
    _cge->apply ( point, dir ); dir *= -1;
    _e->apply ( point, val_new );
    normS_old = dir.normSqr ();

    do {

      adjust ( point, dir, step );

      old_point = point;

      if ( step != 0 ) {

        dir *= step;
        point += dir;
        dir /= step;
      } else {
        dir.setZero ();
      }

      val_old = val_new;
      _e->apply ( point, val_new );
      if ( val_new > val_old ) {
        point = old_point;
        val_new = val_old;
        dir.setZero ();
        _h *= _hDec;
        if ( _h < _hMin ) _h = _hMin;
        _cge->seth ( _h );
        _cde.seth ( _h );
      }

      _cge->apply ( point, grad );

      DataType normS_new = grad.normSqr ();
      DataType beta = normS_new / normS_old;
      
      old_point -= point;
      realstep =  old_point.norm ();

      if (!(i % _reset)) { _h *= _hDec; beta = 0; if (_h < _hMin) _h = _hMin; }

      dir *= beta;
      dir -= grad;

      normS_old = normS_new;
      beta_old = beta;

      cout << aol::color::green << aol::mixedFormat ( sqrt ( normS_old ) )
           << "  " << aol::mixedFormat ( step ) << "  " << aol::mixedFormat ( realstep ) << "  " << aol::mixedFormat ( _h )
           << "  " << aol::detailedFormat ( static_cast<DataType> ( val_old ) )
           << "  " << aol::detailedFormat ( static_cast<DataType> ( val_new ) )
           << "  " << aol::intFormat ( i ) << aol::color::reset << endl;

      cout << point << endl;

      // Save boundary
      aol::Scalar<DataType> nanval = std::numeric_limits<DataType>::infinity ();
      _e->apply ( point, nanval );

    } while ( sqrt ( normS_old ) > _eps && ++i < _maxit );
  }
protected:
  const aol::Op<aol::Vector<DataType>, aol::Scalar<DataType> > * _e;
  aol::Vector<IndexType> _optvars;
  const CDQEnergyGradientOp<DataType, IndexType>* _cge;
  const DQEnergyDerivativeOp<DataType, IndexType> _de;
  const CDQEnergyDerivativeOp<DataType, IndexType> _cde;
  const DataType _eps;
  mutable DataType _h;
  const DataType _hMin, _hDec;
  const DataType _beta, _sigma, _start, _hfac;
  const int _reset, _maxit;
};

//! \brief Generic energy descent using an arbitrary metric tensor with conjugate directions, projecting back into some allowed set in each step
//!        uses difference quotients by default
//! \ingroup Optimization
template <class DataType, class IndexType>
class ProjectingMetricDescent : public SimpleDescent<DataType,IndexType> {
public:

  ProjectingMetricDescent ( const aol::Op<aol::Vector<DataType>, aol::Scalar<DataType> > & e,
          const aol::Op<aol::Vector<DataType>, aol::Vector<DataType> > & proj, // Projecting point back into allowed area
          const aol::Op<aol::MultiVector<DataType>, aol::Vector<DataType> > & gradproj, // Projecting gradient at boundary onto tangent space, if pointing inwards
          const aol::Op<aol::MultiVector<DataType>, aol::Vector<DataType> > & invmetric, // Inverse of metric, mapping differential to gradient
          const aol::Vector<IndexType>& optimizeExcept,
          DataType eps = 1E-3, DataType h = 1E-4, DataType hMin = 1E-6, DataType hDec = 0.5, DataType beta = 0.5,
          DataType sigma = 0.2, DataType start = 0.01, DataType hfac = 2, int reset = 50, int maxit = 1000 )
    : SimpleDescent<DataType,IndexType> ( e, optimizeExcept, eps, h, hMin, hDec, beta, sigma, start, hfac, reset, maxit ), _proj ( proj ), _gradproj ( gradproj ), _invmetric ( invmetric ) {}

  void solve ( aol::Vector<DataType>& point ) const {
    aol::Vector<DataType> dir ( point.size () ), tmppoint ( point.size () ), diff ( point.size () ), grad ( point.size () ), old_point ( point.size () );
    aol::MultiVector<DataType> mvarg; mvarg.appendReference ( point ); mvarg.appendReference ( diff );
    aol::Scalar<DataType> val_old, val_new;
    DataType normS_old, beta_old = 1, step = 0, realstep;
    int i = 1;
    this->_cge->apply ( point, dir ); dir *= -1;
    this->_e->apply ( point, val_new );
    normS_old = dir.normSqr ();

    do {
      mvarg [1] = dir;
      this->_gradproj.apply ( mvarg, dir );
      this->adjust ( point, dir, step );

      //cerr << point << endl << dir << endl;

      old_point = point;

      if ( step != 0 ) {

        dir *= step;
        point += dir;
        dir /= step;
      } else {
        dir.setZero ();
      }

      tmppoint = point;
      this->_proj.apply ( tmppoint, point );

      val_old = val_new;
      this->_e->apply ( point, val_new );

      if ( val_new > val_old ) {
        point = old_point;
        val_new = val_old;
        dir.setZero ();
        this->_h *= this->_hDec;
        if ( this->_h < this->_hMin ) this->_h = this->_hMin;
        this->_cge->seth ( this->_h );
        this->_cde.seth ( this->_h );
      }

      this->_cge->apply ( point, diff );

      //cerr << diff << endl;

      this->_invmetric.apply ( mvarg, grad );

      //cerr << grad << endl;

      DataType normS_new = grad.normSqr ();
      DataType beta = normS_new / normS_old;

      old_point -= point;
      realstep =  old_point.norm ();
      if ( !realstep ) {
  if ( this->_h == this->_hMin && beta_old == 0 ) break;
        this->_h *= this->_hDec; beta = 0; if ( this->_h < this->_hMin ) this->_h = this->_hMin;
      }

      if ( ! ( i % this->_reset ) ) {
        this->_h *= this->_hDec; beta = 0; if ( this->_h < this->_hMin ) this->_h = this->_hMin;
      }

      dir *= beta;
      dir -= grad;

      normS_old = normS_new;
      beta_old = beta;

      cout << aol::color::green << aol::mixedFormat ( sqrt ( normS_old ) )
      << "  " << aol::mixedFormat ( step ) << "  " << aol::mixedFormat ( realstep ) << "  " << aol::mixedFormat ( this->_h )
      << "  " << aol::detailedFormat ( static_cast<DataType> ( val_old ) )
      << "  " << aol::detailedFormat ( static_cast<DataType> ( val_new ) ) << aol::color::reset << endl;

      cout << point << endl;

      // Save boundary
      aol::Scalar<DataType> nanval = std::numeric_limits<DataType>::infinity ();
      this->_e->apply ( point, nanval );

    } while ( sqrt ( normS_old ) > this->_eps && ++i < this->_maxit );
  }

 private:
  const aol::Op<aol::Vector<DataType>, aol::Vector<DataType> > & _proj;
  const aol::Op<aol::MultiVector<DataType>, aol::Vector<DataType> > & _gradproj;
  const aol::Op<aol::MultiVector<DataType>, aol::Vector<DataType> > & _invmetric;
};

//! \brief Gradient descent (no conjugate directions!) for an energy e + d, where e is smooth and d has some kinks.
//!
//! Assumptions: d has a kink at x[i]=x0[i] for the i's listed in dissipationDirections
//! (x0 being the starting point of the descent),
//! with a subdifferential of the form [-m,m] for some m in each of those directions.
//! dd gives derivative of d away from the kink and +m or -m at the kink.
//!
//! For all dissipation directions, store whether the current position is "<=", "=" or ">=" the kink position.
//!
//! If one is at the kink in some directions:
//!
//! Take negative gradient of e+d as a potential descent direction and compute directional derivative in this direction (including dd).
//! If this is positive, set components of the direction (starting from largest) to zero until it is negative.
//! (Positive components can only be "=" kink position, because otherwise a real gradient exists which ensures descent.)
//!
//! Start a line search in the modified direction.
//! Set position to ">=" / "<="  if one starts to move away from the kink in some direction.
//! If line search hits the kink in a ">=" / "<=" direction, stop and set this direction to "=".
//!
//! Repeat until descent direction is close to zero.
//! \ingroup Optimization
template <class DataType, class IndexType>
class RateIndependentDescent : public SimpleDescent<DataType,IndexType> {
public:
  RateIndependentDescent ( const aol::Op<aol::Vector<DataType>, aol::Scalar<DataType> > & e, // Energy
         const aol::Op<aol::Vector<DataType>, aol::Scalar<DataType> > & d, // Dissipated Energy
         const aol::Op<aol::Vector<DataType>, aol::Vector<DataType> > & dd, // Derivative of Dissipation
         const aol::Op<aol::Vector<DataType>, aol::Vector<DataType> > & proj, // Projection for Constraints
         const aol::Vector<IndexType>& dissipationDirections,
         const aol::Vector<IndexType>& optimizeExcept,
         const aol::Vector<DataType>& startpoint,
         DataType eps = 1E-3, DataType h = 1E-4, DataType hMin = 1E-6, DataType hDec = 0.5, DataType beta = 0.5,
                           DataType sigma = 0.2, DataType start = 0.01, DataType hfac = 2, int maxit = 1000, bool debug = false )
    : SimpleDescent<DataType,IndexType> ( e, optimizeExcept, eps, h, hMin, hDec, beta, sigma, start, hfac, 1, maxit ),
    _d ( d ), _dd ( dd ), _projection ( proj ), _dissipationDirections ( dissipationDirections ), _startpoint (startpoint), _debug ( debug ) {}

  //! Full energy
  DataType fullEnergy ( const aol::Vector<DataType>& point ) const {
    aol::Scalar<DataType> val;
    this->_e->apply ( point, val );
    this->_d.applyAdd ( point, val );
    return val;
  }

  //! Full gradient, assume zero dissipation gradient at kink (medgrad) or dissipation opposite to energy gradient (badgrad)
  void fullGradient ( const aol::Vector<DataType>& point, const aol::Vector<signed int>& sign, aol::Vector<DataType>& badgrad, aol::Vector<DataType>& medgrad ) const {

    aol::Vector<DataType> dissmedgrad ( badgrad.size () );
    aol::Vector<DataType> dissbadgrad ( badgrad.size () );
    aol::Vector<DataType> egrad ( badgrad.size () );

    this->_cge->apply ( point, egrad );
    this->_dd.apply ( point, dissbadgrad );
    dissmedgrad = dissbadgrad;

    if ( this->_debug ) {
      cerr << " point " << point << endl;
      cerr << " egrad " << egrad << endl;
      cerr << " dgrad " << dissbadgrad << endl;
    }

    for (int j = 0; j < _dissipationDirections.size (); ++j) {
      int i = _dissipationDirections[j];
      if (sign[j] == 0) {
        dissmedgrad[i]  = 0;
        dissbadgrad[i] *= aol::signum1at0 (-egrad[i]);
      }
    }

    medgrad = egrad; medgrad += dissmedgrad;
    badgrad = egrad; badgrad += dissbadgrad;

    if ( this->_debug ) {
      cerr << " mgrad " << medgrad  << endl;
      cerr << " bgrad " << badgrad  << endl;
    }
  }

  // Forward difference quotient for full energy
  DataType fullDifferenceQuotient ( const aol::Vector<DataType>& point, const aol::Vector<DataType>& dir, double h ) const {
    aol::Vector<DataType> forward_point ( dir ); forward_point *= h; forward_point += point;
    return ( fullEnergy ( forward_point ) - fullEnergy ( point ) ) / h;
  }

  //! Find a descent direction by projecting away dissipation directions
  void findDirection ( const aol::Vector<DataType>& point, aol::Vector<DataType>& dir, aol::Vector<DataType>& badgrad, aol::Vector<signed int>& sign ) const {

    fullGradient ( point, sign, badgrad, dir ); dir *= -1;

    if ( this->_debug ) cerr << " --> testing  " << dir  << endl;
    if ( this->_debug ) cerr << " --> signs    " << sign << endl;

    // dir_deriv contains the full directional derivative in direction dir, before summing
    // (which makes it easy to project out directions by setting components zero)
    aol::Vector<DataType> dir_deriv ( badgrad );
    aol::BitVector settozero ( dir.size () );
    settozero.setAll ( false );
    for (int i = 0; i < dir.size (); ++i)
      dir_deriv[i] *= dir[i];

    // Project away movement directions with too much dissipation, starting with highest dir_deriv component,
    // until total slope gets negative
    while ( true ) {
      std::pair<int,DataType> max = dir_deriv.getMaxIndexAndValue ();
      if ( max.second <= 0 ) break;
      int i = max.first;
      dir_deriv[i] = 0;
      dir[i] = 0;
      settozero.set(i, true);
    }

    // Start to move in directions that have not been projected away
    for ( int j = 0; j < _dissipationDirections.size (); ++j ) {
      int i = _dissipationDirections[j];
      if (sign[j] == 0 && !settozero[i])
        sign[j] = aol::signum (dir[i]);
    }

    if ( this->_debug ) cerr << " <-- modified " << dir  << endl;
    if ( this->_debug ) cerr << " <-- signs    " << sign << endl;
  }

  // In which directions did I hit the kink again?
  void stopAtKink ( const aol::Vector<DataType>& kink_point, const aol::Vector<DataType>& old_point, const aol::Vector<DataType>& dir, aol::Vector<DataType>& new_point, aol::Vector<signed int>& sign ) const {
    if ( this->_debug ) cerr << " --> goto     " << old_point << " + " << dir << endl;
    aol::Vector<DataType> old_diff ( old_point ); old_diff -= kink_point; // Distance old point to kink
    aol::Vector<DataType> new_diff ( old_diff ); new_diff += dir; // Distance new point to kink  - should remain same sign as sign variable
    aol::Vector<DataType> overshoot ( _dissipationDirections.size () );
    for (int j = 0; j < _dissipationDirections.size (); ++j ) {
      int i = _dissipationDirections [j];
      overshoot [j] = - new_diff[i]*sign[j];
    }
    std::pair<int,DataType> max = overshoot.getMaxIndexAndValue ();
    if (max.second > 0) {
      int j = max.first;
      int i = _dissipationDirections [j];
      double factor = - old_diff[i] / dir[i];
      aol::Vector<DataType> fdir ( dir ); fdir *= factor;
      new_point = old_point; new_point += fdir;
      sign[j] = 0;
      new_point[i] = kink_point[i]; // Should be (roughly) nop.
      if ( this->_debug ) cerr << " <-- scale by " << factor << " because of index " << i << endl;
    } else {
      new_point = old_point; new_point += dir;
    }
    if ( this->_debug ) cerr << " <-- goto     " << new_point << endl;
  }

  // Armijo line search
  void adjust ( const aol::Vector<DataType>& point, const aol::Vector<DataType>& dir, const aol::Vector<DataType>& grad, DataType& step ) const {

    DataType m_sec, m_tan = grad * dir;

    DataType norm = dir.norm ();
    if ( !step ) {
      if ( norm ) step = this->_start / norm;
      else step = this->_start;
    }

    if ( this->_debug ) cerr << " --> step     " << step << endl;
    while ( step * norm > this->_start ) step *= this->_beta;
    do {
      step /= this->_beta;
      m_sec = fullDifferenceQuotient ( point, dir, step );
      if (this->_debug) cerr << "+";
    } while ( m_sec / m_tan > this->_sigma && step * norm < this->_start );
    do {
      step *= this->_beta;
      m_sec = fullDifferenceQuotient ( point, dir, step );
      if (this->_debug) cerr << "-";
    } while ( m_sec / m_tan < this->_sigma && step * norm > this->_hfac * this->_h );
    if ( this->_debug ) cerr << " <-- step     " << step << endl;
  }

  void solve ( aol::Vector<DataType>& point ) const {

    int i = 0;
    DataType val_new = fullEnergy ( point ), val_old, step = 0;
    aol::Vector<DataType> kink_point ( _startpoint );
    aol::Vector<DataType> dir ( point.size () );
    aol::Vector<DataType> badgrad ( point.size () );
    aol::Vector<DataType> newpoint ( point.size () );

    // For all directions in which dissipation happens, store whether the current position is greater (+1), equal (0) or less than (-1) the starting position
    aol::Vector<signed int> sign ( _dissipationDirections.size () );
    
    for (int j = 0; j < _dissipationDirections.size (); ++j) {
      int di = _dissipationDirections[i];
      if ( point[di] > _startpoint[di] ) sign[i] =  1;
      if ( point[di] < _startpoint[di] ) sign[i] = -1;
    }

    do {

      findDirection ( point, dir, badgrad, sign );
      adjust ( point, dir, badgrad, step );
      aol::Vector<DataType> diff ( dir ); diff *= step;
      stopAtKink ( kink_point, point, diff, newpoint, sign );

      val_old = val_new;
      val_new = fullEnergy ( newpoint );

      if ( val_new < val_old ) {
        // If energy decreased: do step
        _projection.apply (newpoint, point);
      } else {
        // If energy increased: do not step but refine difference quotient
        val_new = val_old;
        this->_h *= this->_hDec;
        if ( this->_h < this->_hMin ) this->_h = this->_hMin;
        this->_cge->seth ( this->_h );
        this->_cde.seth ( this->_h );
      }

      // If step too small: do it but refine difference quotient
      if ( false && diff.norm () < this->_eps ) {
        if ( this->_h == this->_hMin ) { throw aol::Exception ( "early break", __FILE__, __LINE__, "RateIndependentDescent::solve" ); }
        this->_h *= this->_hDec; if ( this->_h < this->_hMin ) this->_h = this->_hMin;
      }

      cout << aol::color::green << aol::mixedFormat ( badgrad.norm () )
           << "  " << aol::mixedFormat ( diff.norm () ) << "  " << aol::mixedFormat ( this->_h )
           << "  " << aol::detailedFormat ( static_cast<DataType> ( val_old ) )
           << "  " << aol::detailedFormat ( static_cast<DataType> ( val_new ) )
           << "  " << aol::intFormat ( i ) << aol::color::reset << endl;

      cout << "**" << point << endl;

      // Save boundary
      aol::Scalar<DataType> nanval = std::numeric_limits<DataType>::infinity ();
      this->_e->apply ( point, nanval );

    } while ( dir.norm () > this->_eps && ++i < this->_maxit );
  }

private:
  const aol::Op<aol::Vector<DataType>, aol::Scalar<DataType> > & _d;
  const aol::Op<aol::Vector<DataType>, aol::Vector<DataType> > & _dd;
  const aol::Op<aol::Vector<DataType>, aol::Vector<DataType> > & _projection;
  const aol::Vector<IndexType>& _dissipationDirections;
  const aol::Vector<DataType> _startpoint;
  const bool _debug;
};
  

//! \brief Nelder-Mead downhill simplex minimization method for continuous functions \f$ f : \mathbb{R}^n \rightarrow \mathbb{R}\f$
//!
//! Reference: Nelder, John A.; R. Mead (1965). "A simplex method for function minimization". Computer Journal 7: 308–313. doi:10.1093/comjnl/7.4.308.
//!
//! See also: https://en.wikipedia.org/wiki/Nelder–Mead_method
//!           http://www.math.uni-hamburg.de/home/oberle/skripte/optimierung/optim04.pdf (german)
//! \ingroup Optimization
template <typename RealType>
class NelderMeadDownhillSimplexAlgorithm {
  const aol::Op<aol::Vector<RealType>, aol::Scalar<RealType> > *_e;
  const RealType tol;
  const int maxIt;
  bool _quietMode;
public:
  NelderMeadDownhillSimplexAlgorithm ( aol::Op<aol::Vector<RealType>, aol::Scalar<RealType> > &e, const RealType tol = 1e-6, const int maxIt = 100, const bool QuietMode = false )
    : _e ( &e ), tol ( tol ), maxIt ( maxIt ), _quietMode ( QuietMode ) { }
  
  void solve ( aol::Vector<RealType>& point ) const {
    const int dim = point.size ( );
    
    // Initialize simplex
    aol::MultiVector<RealType> x ( dim + 1, dim );
    x.setAll ( point );
    for ( int j=0; j<dim ; ++j ) {
      if ( x[j+1][j] != 0 ) x[j+1][j] *= 1.05;
      else x[j+1][j] += point.getMeanValue ( ) * 0.05;
    }
    
    aol::FullMatrix<RealType> dists ( dim+1, dim+1 );
    aol::Vector<RealType> centroid ( dim ), reflected ( dim ), expanded ( dim ), contracted ( dim );
    aol::Scalar<RealType> energy, reflectedEnergy, expandedEnergy, contractedEnergy;
    aol::Vector<RealType> energies ( dim+1 );
    const RealType alpha = 1.0, gamma = 2.0, rho = -0.5, sigma = 0.5;
    int it = 0;
    aol::ProgressBar<> progressBar ( "Nelder-Mead optimization", std::cerr );
    if ( !_quietMode ) progressBar.start ( maxIt );
    do {
      // 1. Compute all energies (in the Wikipedia version all points are ordererd w.r.t. energy here)
      for ( int j=0; j<dim+1 ; ++j ) {
        _e->apply ( x[j], energy );
        energies[j] = energy;
      }
      
      // 2. Compute centroid
      centroid.setZero ( );
      std::pair<int, RealType> worstEnergyIndVal = energies.getMaxIndexAndValue ( );
      for ( int j=0; j<dim+1 ; ++j ) if ( j != worstEnergyIndVal.first ) centroid += x[j];
      centroid /= dim;
      
      // 3. Reflection
      reflected = centroid;
      reflected.addMultiple ( centroid, alpha );
      reflected.addMultiple ( x[worstEnergyIndVal.first], -alpha );
      _e->apply ( reflected, reflectedEnergy );
      const std::pair<int, RealType> bestEnergyIndVal = energies.getMinIndexAndValue ( );
      energies[bestEnergyIndVal.first] = aol::NumberTrait<RealType>::Inf;
      const std::pair<int, RealType> secondBestEnergyIndVal = energies.getMinIndexAndValue ( );
      if ( reflectedEnergy >= bestEnergyIndVal.second && reflectedEnergy < secondBestEnergyIndVal.second ) {
        x[worstEnergyIndVal.first] = reflected;
      } else {
        // 4. Expansion
        if ( reflectedEnergy < bestEnergyIndVal.second ) {
          expanded = centroid;
          expanded.addMultiple ( centroid, gamma );
          expanded.addMultiple ( x[worstEnergyIndVal.first], -gamma );
          _e->apply ( expanded, expandedEnergy );
          if ( expandedEnergy < reflectedEnergy ) {
            x[worstEnergyIndVal.first] = expanded;
          } else {
            x[worstEnergyIndVal.first] = reflected;
          }
        } else {
          // 5. Contraction
          contracted = centroid;
          contracted.addMultiple ( centroid, rho );
          contracted.addMultiple ( x[worstEnergyIndVal.first], -rho );
          _e->apply ( contracted, contractedEnergy );
          if ( contractedEnergy < worstEnergyIndVal.second ) {
            x[worstEnergyIndVal.first] = contracted;
          } else {
            // 6. Reduction
            for ( int j=0; j<dim+1 ; ++j ) {
              if ( j != bestEnergyIndVal.first ) {
                x[j] *= ( 1 + sigma );
                x[j].addMultiple ( x[bestEnergyIndVal.first], 1 - sigma );
              }
            }
          }
        }
      }
      ++it;
      if ( !_quietMode ) progressBar++;
    } while ( getMaxSimplexNodeDistance ( x ) > tol && getCentroidVariance ( x, centroid ) > tol && getEnergyVariance ( x ) > tol && it < maxIt );
    if ( !_quietMode ) progressBar.finish ( );
    
    // Compute all energies and set result to point with least energy
    for ( int j=0; j<dim+1 ; ++j ) {
      _e->apply ( x[j], energy );
      energies[j] = energy;
    }
    point = x[energies.getMinIndexAndValue ( ).first];
  }
protected:
  RealType getMaxSimplexNodeDistance ( const aol::MultiVector<RealType> &x ) const {
    RealType res = 0.0, dist;
    aol::Vector<RealType> diff ( x[0].size ( ) );
    for ( int j=0; j<x.numComponents ( ) ; ++j ) {
      for ( int i=j+1; i<x.numComponents ( ) ; ++i ) {
        diff = x[j];
        diff -= x[i];
        dist = diff.normSqr ( );
        if ( dist > res ) res = dist;
      }
    }
    return res;
  }
  
  RealType getCentroidVariance ( aol::MultiVector<RealType> &x, aol::Vector<RealType> &centroid ) const {
    RealType res = 0.0;
    aol::Scalar<RealType> centroidEnergy, energy;
    _e->apply ( centroid, centroidEnergy );
    for ( int j=0; j<x.numComponents ( ) ; ++j ) {
      _e->apply ( x[j], energy );
      res += aol::Sqr<RealType> ( energy - centroidEnergy );
    }
    return res / static_cast<RealType> ( x.numComponents ( ) );
  }
    
  RealType getEnergyVariance ( aol::MultiVector<RealType> &x ) const {
    aol::Vector<RealType> energies ( x.numComponents ( ) );
    aol::Scalar<RealType> energy;
    for ( int j=0; j<x.numComponents ( ) ; ++j ) {
      _e->apply ( x[j], energy );
      energies[j] = energy;
    }
    energies.addToAll ( -energies.getMeanValue ( ) );
    return energies.getVariance ( );
  }
};

} // end of namespace aol.

#endif
