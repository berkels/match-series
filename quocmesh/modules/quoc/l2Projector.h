#ifndef __L2PROJECTOR_H
#define __L2PROJECTOR_H

#include <quoc.h>
#include <matrix.h>
#include <op.h>

namespace qc {

const int QC_SELECT_BASIS_2D_CONSTANT = 1;
const int QC_SELECT_BASIS_2D_X        = ( 1 << 1 );
const int QC_SELECT_BASIS_2D_Y        = ( 1 << 2 );
const int QC_SELECT_BASIS_2D_XY       = ( 1 << 3 );
const int QC_SELECT_BASIS_2D_XX       = ( 1 << 4 );
const int QC_SELECT_BASIS_2D_YY       = ( 1 << 5 );
const int QC_SELECT_BASIS_2D_ALL_LINEAR = QC_SELECT_BASIS_2D_X | QC_SELECT_BASIS_2D_Y;
const int QC_SELECT_BASIS_2D_ALL_SECOND = QC_SELECT_BASIS_2D_XY | QC_SELECT_BASIS_2D_XX | QC_SELECT_BASIS_2D_YY;
const int QC_SELECT_BASIS_2D_ALL = QC_SELECT_BASIS_2D_CONSTANT | QC_SELECT_BASIS_2D_ALL_LINEAR | QC_SELECT_BASIS_2D_ALL_SECOND;

template <typename RealType>
class BasisFunction2D {
protected:
  static RealType _1 ( RealType  , RealType ) {
    return 1.;
  }
  static RealType _x ( RealType X, RealType ) {
    return X;
  }
  static RealType _y ( RealType  , RealType Y ) {
    return Y;
  }
  static RealType _xy ( RealType X, RealType Y ) {
    return X*Y;
  }
  static RealType _xx ( RealType X, RealType ) {
    return X*X;
  }
  static RealType _yy ( RealType  , RealType Y ) {
    return Y*Y;
  }
  typedef RealType ( *BASIS_FUNC_TYPE ) ( RealType, RealType );
public:

  class IncompBasisEx {
  public:
    IncompBasisEx() {}
    int j;
  };

  virtual ~BasisFunction2D() {} // subclasses have virtual functions

  BasisFunction2D() : _func ( NULL ) {}

  virtual aol::ExitStatus scalarProduct ( BasisFunction2D<RealType> &B, RealType &Prod ) = 0;

  RealType operator() ( RealType X, RealType Y ) {
    return _func ( X, Y );
  }

protected:
  BASIS_FUNC_TYPE _func;
};

template <typename RealType>
class Basis2DConst : public qc::BasisFunction2D<RealType> {
public:
  Basis2DConst() : qc::BasisFunction2D<RealType>() {
    this->_func = this->_1;
  }

  virtual aol::ExitStatus scalarProduct ( BasisFunction2D<RealType> &B, RealType &Prod );
};

template <typename RealType>
class Basis2DX : public qc::BasisFunction2D<RealType> {
public:
  Basis2DX() : qc::BasisFunction2D<RealType>() {
    this->_func = this->_x;
  }

  virtual aol::ExitStatus scalarProduct ( BasisFunction2D<RealType> &B, RealType &Prod );
};

template <typename RealType>
class Basis2DY : public qc::BasisFunction2D<RealType> {
public:
  Basis2DY() : qc::BasisFunction2D<RealType>() {
    this->_func = this->_y;
  }

  virtual aol::ExitStatus scalarProduct ( BasisFunction2D<RealType> &B, RealType &Prod );
};


template <typename RealType>
class Basis2DXY : public qc::BasisFunction2D<RealType> {
public:
  Basis2DXY() : qc::BasisFunction2D<RealType>() {
    this->_func = this->_xy;
  }

  virtual aol::ExitStatus scalarProduct ( BasisFunction2D<RealType> &B, RealType &Prod );
};

template <typename RealType>
class Basis2DXX : public qc::BasisFunction2D<RealType> {
public:
  Basis2DXX() : qc::BasisFunction2D<RealType>() {
    this->_func = this->_xx;
  }

  virtual aol::ExitStatus scalarProduct ( BasisFunction2D<RealType> &B, RealType &Prod );
};

template <typename RealType>
class Basis2DYY : public qc::BasisFunction2D<RealType> {
public:
  Basis2DYY() : qc::BasisFunction2D<RealType>() {
    this->_func = this->_yy;
  }

  virtual aol::ExitStatus scalarProduct ( BasisFunction2D<RealType> &B, RealType &Prod );
};

/** Class L2Projector is used for performing L2-projections of functions onto
 *  discrete function spaces
 *  \author Droske
 */
template <typename RealType>
class L2Projector  {
public:
  enum QuadraturOrder { FIRST_ORDER, SECOND_ORDER };

  /** Constructor for Projection
   * @param PatchSize sample size of the patch on which projection will be performed
   * @param Dim   defines whether the class should work in 2D or 3D */
  L2Projector ( const int PatchSize, Dimension Dim )
      : _patchSize ( PatchSize ), _dim ( Dim ),
      _m ( NULL ), _minv ( NULL ), _eval ( NULL ), _rhs ( NULL ),
      _w ( PatchSize - aol::ZOTrait<RealType>::one ), _order ( SECOND_ORDER ) {
        _eval = new aol::FullMatrix<RealType> ( ( PatchSize - 1 ) *2, ( PatchSize - 1 ) *2 );
  }

  ~L2Projector() {
    if ( _m ) delete _m;
    if ( _minv ) delete _minv;
    if ( _rhs ) delete _rhs;
    if ( _eval ) delete _eval;
  }

  /** Member function to add a basis function to the discrete Function space
   * @param B an object derived from BasisFunction2D
   */
  void addBasisFunction ( BasisFunction2D<RealType> &B ) {
    _basisFunctions.push_back ( &B );
  }

  /** Member function for the initialization of the operator.
   *  MUST BE CALLED BEFORE SUBSEQUENT CALLS OF APPLY.
   */
  void initBasis();

  template <class EVALUATOR_TYPE>
  void computeProjection ( aol::Vector<RealType> &Dest, EVALUATOR_TYPE &LocalEvaluator ) {
    if ( static_cast< int > ( Dest.size() ) != static_cast< int > ( _basisFunctions.size() ) ) {
      throw aol::Exception ( "Destination vector size unequal to number of basis functions", __FILE__, __LINE__ );
    }
    buildRHS ( LocalEvaluator );
    _minv->mult ( *_rhs, Dest );
  }

protected:

  template <class EVALUATOR_TYPE>
  void buildRHS ( EVALUATOR_TYPE &LocalEvaluator ) {

    const RealType alpha[2] = {
                                0.21132487, 0.78867513
                              };

    for ( int i = 0; i < _patchSize - 1; i++ ) {
      for ( int j = 0; j < _patchSize - 1; j++ ) {
        for ( int k = 0; k < 2; k++ ) {
          for ( int l = 0; l < 2; l++ ) {
            _eval->set ( i * 2 + k, j * 2 + l, LocalEvaluator ( i, j, alpha[ k ], alpha[ l ] ) );
          }
        }
      }
    }
    switch ( _order ) {
    case SECOND_ORDER: {
      const int numBasis = _basisFunctions.size();
      for ( int b = 0; b < numBasis; b++ ) {
        RealType v = 0.;
        BasisFunction2D<RealType> &func = *_basisFunctions[ b ];
        for ( int i = 0; i < _patchSize - 1; i++ ) {
          for ( int j = 0; j < _patchSize - 1; j++ ) {
            for ( int k = 0; k < 2; k++ ) {
              for ( int l = 0; l < 2; l++ ) {
                //float x = (i + alpha[ k ])/_w;
                //float y = (j + alpha[ l ])/_w;
                // v += LocalEvaluator( i, j, alpha[ k ], alpha[ l ] ) * func( (i + alpha[ k ])/_w, (j + alpha[ l ])/_w );
                v += _eval->get ( i * 2 + k, j * 2 + l ) * func ( ( i + alpha[ k ] ) / _w, ( j + alpha[ l ] ) / _w );
              }
            }
          }
        }
        v /= _w * _w * 4.;
        ( *_rhs ) [ b ] = v;
      }
    }
    break;
    default:
      throw aol::Exception ( "only second order integration implemented for L2Projector", __FILE__, __LINE__ );
      break;
    }
  }

  vector<BasisFunction2D<RealType>* > _basisFunctions;
  const int _patchSize;
  const Dimension _dim;
  aol::FullMatrix<RealType> *_m, *_minv, *_eval;
  aol::Vector<RealType> *_rhs;

  RealType _w;
  QuadraturOrder _order;
};


template <typename RealType>
void qc::L2Projector<RealType>::initBasis() {
  int numBasis = static_cast<int> ( _basisFunctions.size() );

  if ( _m ) delete _m;
  _m = new aol::FullMatrix<RealType> ( numBasis, numBasis );

  if ( _minv ) delete _minv;
  _minv = new aol::FullMatrix<RealType> ( numBasis, numBasis );

  if ( _rhs ) delete _rhs;
  _rhs = new aol::Vector<RealType> ( numBasis );

  for ( int i = 0; i < numBasis; i++ ) {
    for ( int j = i; j < numBasis; j++ ) {
      RealType v = 0;
      if ( _basisFunctions[ i ]->scalarProduct ( *_basisFunctions[ j ], v ) == aol::FAILURE &&
           _basisFunctions[ j ]->scalarProduct ( *_basisFunctions[ i ], v ) == aol::FAILURE ) {
        throw aol::Exception ( "basisFunctions don't know how to integrate their product", __FILE__, __LINE__ );
      }
      _m->set ( i, j, v );
      _m->set ( j, i, v );
    }
  }

  // cerr << "M:\n";
  //   _m->dump();

  _minv->makeInverse ( *_m );

  //  cerr << "M inverse:\n";
  //_minv->dump();

  //aol::FullMatrix<RealType> id( numBasis, numBasis );
  //id = *_m;
  //id *= *_minv;
}

} // end namespace qc

#endif
