#ifndef __PROJECTORS_H
#define __PROJECTORS_H

#include <aol.h>
#include <matrix.h>


namespace aol {


template <typename _RealType, typename _VectorType, typename _MatrixType>
class Projector {
  typedef _RealType RealType;
  typedef _VectorType VectorType;
  typedef _MatrixType MatrixType;
public:
  Projector ( ) { }

  virtual ~Projector ( ) { }

  virtual void apply ( const VectorType &/*Arg*/, VectorType &/*Dest*/ ) const = 0;

  virtual bool isFeasible ( const VectorType &/*X*/ ) const = 0;
};


/**
 *  \brief Provides a class that projects vectors onto a box (possibly only in part of the dimensions)
 *
 *  \f[x \in \mathbb{R}^n is projected onto B = \{ x \in \mathbb{R}^n : l_i \leq x_i \leq u_i \forall i \text{ with } c_i=1 \}\f],
 *  where \f[l \in \mathbb{R}^n\f] are the lower bounds, \f[u \in \mathbb{R}^n\f] the upper bounds
 *  and \f[c \in \{0,1\}^n\f] is a BitVector indicating for each dimension if it should be constrained or not.
 *
 *  \author mevenkamp
 */
template <typename _RealType, typename _VectorType>
class BoxProjector : public Projector<_RealType, _VectorType, aol::FullMatrix<_RealType> > {
  typedef _RealType RealType;
  typedef _VectorType VectorType;
protected:
  VectorType _lowerBounds, _upperBounds;
  aol::BitVector _constrainedDirections;
public:
  BoxProjector ( const VectorType &LowerBounds, const VectorType &UpperBounds )
    : Projector<RealType, VectorType, aol::FullMatrix<_RealType> > ( ), _lowerBounds ( LowerBounds ), _upperBounds ( UpperBounds ),
      _constrainedDirections ( LowerBounds.Size ( ), true ) {
    if ( LowerBounds.size ( ) != UpperBounds.size ( ) )
      throw aol::Exception ( "Lower and upper bounds dimensions do not match!", __FILE__, __LINE__ );
  }

  BoxProjector ( const VectorType &LowerBounds, const VectorType &UpperBounds, const aol::BitVector &ConstrainedDirections )
    : Projector<RealType, VectorType, aol::FullMatrix<_RealType> > ( ), _lowerBounds ( LowerBounds ), _upperBounds ( UpperBounds ),
      _constrainedDirections ( ConstrainedDirections ) {
    if ( LowerBounds.size ( ) != UpperBounds.size ( ) )
      throw aol::Exception ( "Lower and upper bounds dimensions do not match!", __FILE__, __LINE__ );

    if ( ConstrainedDirections.size ( ) != LowerBounds.size ( ) )
      throw aol::Exception ( "Dimension of vector specifying the constrained directions does not match lower/upper variable bounds!", __FILE__, __LINE__ );
  }

  void apply ( const VectorType &Arg, VectorType &Dest ) const {
    if ( Arg.size ( ) != Dest.size ( ) || Arg.size ( ) != _lowerBounds.size ( ) )
      throw aol::Exception ( "Given variable dimensions do not match the dimensions of previously specified bounds!", __FILE__, __LINE__ );

    Dest = Arg;
    for ( int i=0; i<Arg.size ( ) ; ++i ) {
      if ( _constrainedDirections[i] ) {
        if ( Arg[i] < _lowerBounds[i] )
          Dest[i] = _lowerBounds[i];
        else if ( Arg[i] > _upperBounds[i] )
          Dest[i] = _upperBounds[i];
      }
    }
  }

  bool isFeasible ( const VectorType &Arg ) const {
    if ( Arg.size ( ) != _lowerBounds.size( ) )
      throw aol::Exception ( "Given variable dimension does not match the dimensions of previously specified bounds!", __FILE__, __LINE__ );

    for ( int i=0; i<Arg.size ( ) ; ++i ) {
      if ( Arg[i] < _lowerBounds[i] || Arg[i] > _upperBounds[i] )
        return false;
    }
    return true;
  }
};

/**
 *  \brief Provides a class that projects vectors onto the canonical simplex
 * 
 *  \f[x \in \mathbb{R}^n is projected onto U = \{ x \in \mathbb{R}^n : x_i \geq 0 \forall i, \sum_i x_i = 1 \}\f]
 *
 *  Algorithm implemented based on the following paper:
 *    C. Michelot: A Finite Algorithm for Finding the Projection of a Point onto the Canonical Simplex of \f[\mathbb{R}^n\f].
 *    J. Optim. Theory Appl., 50(1):195--200, 1986.
 *
 *  \author mevenkamp
 */
template <typename _RealType, typename _VectorType>
class CanonicalSimplexProjector : public Projector<_RealType, _VectorType, aol::FullMatrix<_RealType> > {
  typedef _RealType RealType;
  typedef _VectorType VectorType;
public:
  CanonicalSimplexProjector ( ) : Projector<RealType, VectorType, aol::FullMatrix<_RealType> > ( ) { }
  
  void apply ( const VectorType &Arg, VectorType &Dest ) const {
    VectorType x ( Arg );
    Dest.setAll ( -1.0 );
    std::set<unsigned int> indices;
    while ( Dest.getMinValue ( ) < 0 ) {
      RealType xSum = x.sum ( );
      for ( int i=0; i<Dest.size ( ) ; ++i ) {
        if ( indices.find ( i ) == indices.end ( ) ) Dest[i] = x[i] - ( xSum - 1 ) / ( Arg.size ( ) - indices.size ( ) );
        else Dest[i] = 0.0;
      }
      if ( Dest.getMinValue ( ) < 0 ) {
        for ( int i=0; i<Dest.size ( ) ; ++i ) {
          if ( Dest[i] < 0 ) {
            indices.insert ( i );
            x[i] = 0.0;
          } else x[i] = Dest[i];
        }
      }
    }
  }
  
  bool isFeasible ( const VectorType &X ) const {
    return ( X.getMinValue ( ) >= 0 && X.sum ( ) == 1 );
  }
};
  
  
} // end namespace


#endif
