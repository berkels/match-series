#ifndef __ESTIMATOR2D_H
#define __ESTIMATOR2D_H

#include <scalarArray.h>

namespace qc {

  static const signed char EST_SAT_NONE = 0;
  static const signed char EST_SAT_UP   = 1;
  static const signed char EST_SAT_DOWN = 2;
  static const signed char EST_SAT_NODE = 4;

//! Estimator class 2d.
template <class T>
class Estimator2d : public ScalarArray<T, qc::QC_2D> {
public:
  /** The constuctor.
   *  @param Size is the size of the data in one dimension, e.g. 129 or 128
   *  depending on whether you hand over a value of the form \f$ 2^n \f$ or
   *  \f$ 2^n +1 \f$, the class will choose a nodal or an element-based
   *  saturation.
   */
  Estimator2d ( int Size );

  ~Estimator2d();

  /** Return the depth of this hierarchical data
   */
  int          getMaxLevel() const {
    return maxLevel;
  }

  /** Return the dimension of this data
   */
  Dimension  getDimOfWorld() const {
    return QC_2D;
  }

  void         makeSaturationDown();

  void         makeSaturationUp();

  void         setThreshold ( T Threshold ) {
    threshold = Threshold;
  }

  //! method walks across the data-array and overwrites it
  //! with the saturated error-array
  void makeSaturatedErrorArray ( ScalarArray<T, qc::QC_2D>* );

  // auxiliary functions needed for saturation
  inline T     maxChildNodeValue ( int X, int Y, int Level );

  inline T     minChildNodeValue ( int X, int Y, int Level );


  //! Returns boolean. Checks if the estimator is ok for the
  //! Element at X, Y at Level according to Thres.
  inline bool  checkElement ( int X, int Y, int Z, int Level, T Thres ) const;

  inline bool  checkElement ( int X, int Y, int Z, int Level ) const {
    return checkElement ( X, Y, Z, Level, threshold );
  }

  //! Checks estimator for element El with threshold Thres.
  inline bool  checkElement ( const Element &El, T Thres ) const {
    return checkElement ( El.x(), El.y(), El.z(), El.level(), Thres );
  }

  inline bool  checkElement ( const Element &El ) const {
    return checkElement ( El, threshold );
  }

  virtual inline bool checkElementDown ( int X, int Y,
                                         int Level, T Thres ) const ;

  virtual inline bool  checkElementUp ( int X, int Y,
                                        int Level, T Thres ) const ;

  void reset() {
    this->setZero();
    sat_type &= ~ ( EST_SAT_UP | EST_SAT_DOWN );
  }

protected:

  void         makeSaturationDownNode();

  void         makeSaturationUpNode();

  void         makeSaturationDownElement();

  void         makeSaturationUpElement();

  int      maxLevel;

  // stores flags to indicate saturation state
  // (upwards, downwards, elementwise, nodal)
  int      sat_type;

  ScalarArray<T, qc::QC_2D> **elSatArrays;
  T threshold;                 //!< The actual threshold

};

template <class T>
T Estimator2d<T>::maxChildNodeValue ( int X, int Y, int halfStep ) {
  int fullStep = halfStep << 1;

  T max, v;
  max = this->get ( X + halfStep, Y + halfStep );

  for ( int lX = X; lX <= X + fullStep; lX += halfStep ) {
    for ( int lY = Y; lY <= Y + fullStep; lY += halfStep ) {
      if ( ( lX % fullStep ) != 0
           || ( lY % fullStep ) != 0 ) {
        v = this->get ( lX, lY );
        if ( v > max ) max = v;
      }
    }
  }
  return max;
}

template <class T>
T Estimator2d<T>::minChildNodeValue ( int X, int Y, int halfStep ) {
  int fullStep = halfStep << 1;

  T min, v;
  min = this->get ( X + halfStep, Y + halfStep );

  for ( int lX = X; lX <= X + fullStep; lX += halfStep ) {
    for ( int lY = Y; lY <= Y + fullStep; lY += halfStep ) {
      if ( lX % fullStep != 0
           || lY % fullStep != 0 ) {
        v = this->get ( lX, lY );
        if ( v < min ) min = v;

      }
    }
  }
  return min;
}

template <class T>
bool Estimator2d<T>::checkElement ( int X, int Y, int, int Level, T Thres ) const {

  switch ( sat_type & ( EST_SAT_UP | EST_SAT_DOWN ) ) {
  case EST_SAT_DOWN:
    return checkElementDown ( X, Y, Level, Thres );
  case EST_SAT_UP:
    return checkElementUp ( X, Y, Level, Thres );
  case EST_SAT_NONE:
    throw aol::Exception ( "Estimator2d::checkElement: not saturated yet",
                           __FILE__, __LINE__ );
  default:
    return 0;
  }
}

template <class T>
bool Estimator2d<T>::checkElementDown ( int X, int Y, int Level, T Thres ) const {
  if ( sat_type & EST_SAT_NODE ) {
    // for nodal saturation check all values on nodes
    int fullStep = 1 << ( maxLevel - Level );
    int halfStep = fullStep >> 1;

    //cerr << "level is " << Level << "/" << maxLevel << endl;
    //cerr << X << " | " << Y << " (" << fullStep << ", " << halfStep << ")\n";

    if ( Thres < this->get ( X + fullStep, Y + halfStep ) ) return 0;
    if ( Thres < this->get ( X + halfStep, Y + fullStep ) ) return 0;
    if ( Thres < this->get ( X           , Y + halfStep ) ) return 0;
    if ( Thres < this->get ( X + halfStep, Y ) ) return 0;
    if ( Thres < this->get ( X + halfStep, Y + halfStep ) ) return 0;

    return 1;
  } else {
    int Shift = ( maxLevel - Level );
    return ( elSatArrays[ Level ]->get ( X >> Shift, Y >> Shift ) <= Thres );
  }
}

template <class T>
bool Estimator2d<T>::checkElementUp ( int X, int Y, int Level, T Thres ) const {

  if ( sat_type & EST_SAT_NODE ) {
    // for nodal saturation check all values on nodes
    int fullStep = 1 << ( maxLevel - Level );
    int halfStep = fullStep >> 1;

    if ( Thres > this->get ( X + fullStep, Y + halfStep ) ) return 0;
    if ( Thres > this->get ( X + halfStep, Y + fullStep ) ) return 0;
    if ( Thres > this->get ( X           , Y + halfStep ) ) return 0;
    if ( Thres > this->get ( X + halfStep, Y ) ) return 0;
    if ( Thres > this->get ( X + halfStep, Y + halfStep ) ) return 0;

    return 1;
  } else {
    int Shift = ( maxLevel - Level );
    return ( elSatArrays[ Level ]->get ( X >> Shift, Y >> Shift ) >= Thres );
  }
}

}

#endif
