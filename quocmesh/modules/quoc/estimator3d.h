#ifndef __ESTIMATOR3D_H
#define __ESTIMATOR3D_H

#include <scalarArray.h>
#include <estimator2d.h>

namespace qc {

/*!
  Implemenation of a class containing saturated error indicators.
 */
template <typename DataType>
class Estimator3d : public ScalarArray<DataType, qc::QC_3D> {
public:
  Estimator3d ( int Size );

  ~Estimator3d();

  /** Return the depth of this hierarchical data
   */
  int          getMaxLevel() const {
    return maxLevel;
  }

  /** Return the dimension of this data
   */
  Dimension  getDimOfWorld() const {
    return QC_3D;
  }

  //! Saturate the indicator to be robust downwards.
  /*! Downwards saturation means
    \f[ \eta(E) \geq \eta(E_C) \quad \forall E_C \in C(E)  \f] */
  void         makeSaturationDown();

  //! Saturate the indicator to be robust upwards.
  void         makeSaturationUp();

  //! Set the current threshold for the estimator
  void         setThreshold ( DataType Threshold ) {
    threshold = Threshold;
  }

  //! method walks across the data-array and overwrites it
  //! with the saturated error-array
  void makeSaturatedErrorArray ( ScalarArray<DataType, qc::QC_3D>* );

  // auxiliary functions needed for saturation
  inline DataType     maxChildNodeValue ( int X, int Y, int Z, int Level ) const;

  inline DataType     minChildNodeValue ( int X, int Y, int Z, int Level ) const;

  inline int   checkElement ( const Element &El ) const {
    return checkElement ( El, threshold );
  }

  inline int   checkElement ( const Element &El, DataType Thres ) const {
    return checkElement ( El.x(), El.y(), El.z(), El.level(), Thres );
  }

  //! Returns boolean. Checks if the estimator is ok for the
  //! Element at X, Y, Z at Level according to Thres.
  inline int   checkElement ( int X, int Y, int Z, int Level, DataType Thres ) const;

  inline int   checkElement ( int X, int Y, int Z, int Level ) const {
    return checkElement ( X, Y, Z, Level, threshold );
  }

  virtual inline bool   checkElementDown ( int X, int Y, int Z, int Level, DataType Thres ) const;

  virtual inline bool   checkElementUp ( int X, int Y, int Z, int Level, DataType Thres ) const;

  void reset() {
    this->setZero();
    sat_type &= ~ ( EST_SAT_UP | EST_SAT_DOWN );
  }

protected:

  void         makeSaturationDownNode();

  void         makeSaturationUpNode();

  void         makeSaturationDownElement();

  void         makeSaturationUpElement();

  int          maxLevel;

  // stores flags to indicate saturation state
  // (upwards, downwards, elementwise, nodal)
  int          sat_type;

protected:

  ScalarArray<DataType, qc::QC_3D> **elSatArrays;

  DataType     threshold;
};

template <typename DataType>
DataType Estimator3d<DataType>::maxChildNodeValue ( int X, int Y, int Z, int halfStep ) const {
  int fullStep = halfStep << 1;

  DataType max, v;
  max = this->get ( X + halfStep, Y + halfStep, Z + halfStep );

  for ( int lX = X; lX <= X + fullStep; lX += halfStep ) {
    for ( int lY = Y; lY <= Y + fullStep; lY += halfStep ) {
      for ( int lZ = Z; lZ <= Z + fullStep; lZ += halfStep ) {
        if ( lX % fullStep != 0
             || lY % fullStep != 0
             || lZ % fullStep != 0 ) {
          v = this->get ( lX, lY, lZ );
          if ( v > max ) max = v;
        }
      }
    }
  }
  return max;
}

template <typename DataType>
DataType Estimator3d<DataType>::minChildNodeValue ( int X, int Y, int Z, int halfStep ) const {
  int fullStep = halfStep << 1;

  DataType min, v;
  min = this->get ( X + halfStep, Y + halfStep, Z + halfStep );

  for ( int lX = X; lX <= X + fullStep; lX += halfStep ) {
    for ( int lY = Y; lY <= Y + fullStep; lY += halfStep ) {
      for ( int lZ = Z; lZ <= Z + fullStep; lZ += halfStep ) {
        if ( lX % fullStep != 0
             || lY % fullStep != 0
             || lZ % fullStep != 0 ) {
          v = this->get ( lX, lY, lZ );
          if ( v < min ) min = v;
        }
      }
    }
  }
  return min;
}

template <typename DataType>
int Estimator3d<DataType>::checkElement ( int X, int Y, int Z, int Level, DataType Thres ) const {

  //if ( Level < maxLevel ) return 0;
  //else return 1;

  switch ( sat_type & ( EST_SAT_UP | EST_SAT_DOWN ) ) {
  case EST_SAT_DOWN:
    return checkElementDown ( X, Y, Z, Level, Thres );
  case EST_SAT_UP:
    return checkElementUp ( X, Y, Z, Level, Thres );
  case EST_SAT_NONE:
    cerr << "ERROR in qcEstimator::checkElement: not saturated yet!\n";
    /* FALLTHRU */
  default:
    return 0;
  }
}

template <class T>
bool Estimator3d<T>::checkElementDown ( int X, int Y, int Z, int Level, T Thres ) const {
  if ( sat_type & EST_SAT_NODE ) {
    if ( maxChildNodeValue ( X, Y, Z, Level ) > Thres ) return 0;
    return 1;
  } else {
    int Shift = ( maxLevel - Level );
    return ( elSatArrays[ Level ]->get ( X >> Shift, Y >> Shift, Z >> Shift ) <= Thres );
  }
}

template <class T>
bool Estimator3d<T>::checkElementUp ( int X, int Y, int Z, int Level, T Thres ) const {

  if ( sat_type & EST_SAT_NODE ) {
    if ( minChildNodeValue ( X, Y, Z, Level ) < Thres ) return 0;
    return 1;
  } else {
    int Shift = ( maxLevel - Level );
    return ( elSatArrays[ Level ]->get ( X >> Shift, Y >> Shift, Z >> Shift ) >= Thres );
  }
}

}

#endif
