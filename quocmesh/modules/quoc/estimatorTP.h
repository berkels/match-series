#ifndef __ESTIMATORTP_H
#define __ESTIMATORTP_H

#include <scalarArray.h>


/**
* Estimator class 2d.
*
*  \author Paetz
*/
template <class T>
class EstimatorTP2d : public qc::ScalarArray<T, qc::QC_2D> {
public:
  /** The constuctor.
   *  @param Size is the size of the data in one dimension, e.g. 129 or 128
   *  depending on whether you hand over a value of the form \f$ 2^n \f$ or
   *  \f$ 2^n +1 \f$, the class will choose a nodal or an element-based
   *  saturation.
   */
  EstimatorTP2d ( int Size );

  ~EstimatorTP2d();

  /** Return the depth of this hierarchical data
   */
  int          getMaxLevel() const {
    return maxLevel;
  }

  /** Return the dimension of this data
   */
  qc::Dimension  getDimOfWorld() const {
    return qc::QC_2D;
  }


  void         setThreshold ( T Threshold ) {
    threshold = Threshold;
  }

  T getThreshold() {
    return threshold;
  }

  //! method walks across the data-array and overwrites it
  //! with the saturated error-array
  void makeSaturatedErrorArray ( const qc::ScalarArray<T, qc::QC_2D> &meshData );


  //! method walks across the data-array and overwrites it
  //! with the saturated error-array
  void setMinOfBothArrays ( const qc::ScalarArray<T, qc::QC_2D> &oldError ) {
    for ( int i = 0; i < oldError.getNumX(); ++i ) { // x-direction
      for ( int j = 0; j < oldError.getNumY(); ++j ) { // y-direction
        if ( oldError.get ( i, j ) < this->get ( i, j ) ) {
          this->set ( i, j, oldError.get ( i, j ) );
        }
      }
    }
  }

  //! method walks across the data-array and overwrites it
  //! with the saturated error-array
  void setErrorArray ( const qc::ScalarArray<T, qc::QC_2D> &error ) {
    for ( int i = 0; i < error.getNumX(); ++i ) { // x-direction
      for ( int j = 0; j < error.getNumY(); ++j ) { // y-direction
        this->set ( i, j, error.get ( i, j ) );
      }
    }
  }

  T dx ( const qc::ScalarArray<T, qc::QC_2D> &data, const int &x, const int &y ) const {
    return data.getClip ( x + 1, y ) - data.get ( x, y );
  }

  T dy ( const qc::ScalarArray<T, qc::QC_2D> &data, const int &x, const int &y ) const {
    return data.getClip ( x , y + 1 ) - data.get ( x, y );
  }


  //! Returns boolean. Checks if the estimator is ok for the
  //! Element at X, Y at Level according to Thres.
  //inline bool  checkElement ( int X, int Y, int Z, int Level, T Thres ) const;

  bool checkElement ( int X, int Y, int, int Level, T Thres ) const {
    int fullStep = 1 << ( maxLevel - Level );
    //    int halfStep = fullStep >> 1;

    if ( Thres > this->get ( X + fullStep, Y + fullStep ) ) return 0;
    if ( Thres > this->get ( X , Y + fullStep ) ) return 0;
    if ( Thres > this->get ( X           , Y  ) ) return 0;
    if ( Thres > this->get ( X + fullStep, Y ) ) return 0;

    return 1;
  }


  inline bool  checkElement ( int X, int Y, int Z, int Level ) const {
    return checkElement ( X, Y, Z, Level, threshold );
  }

  //! Checks estimator for element El with threshold Thres.
  inline bool  checkElement ( const qc::Element &El, T Thres ) const {
    return checkElement ( El.x(), El.y(), El.z(), El.level(), Thres );
  }

  inline bool  checkElement ( const qc::Element &El ) const {
    return checkElement ( El, threshold );
  }

  void reset() {
    this->setZero();
  }

protected:

  int      maxLevel;

  T threshold;                 //!< The actual threshold

};

/*!
  Implemenation of a class containing saturated error indicators.
 */
template <typename DataType>
class EstimatorTP3d : public qc::ScalarArray<DataType, qc::QC_3D> {
public:
  EstimatorTP3d ( int Size );

  ~EstimatorTP3d();

  /** Return the depth of this hierarchical data
   */
  int          getMaxLevel() const {
    return maxLevel;
  }

  /** Return the dimension of this data
   */
  qc::Dimension  getDimOfWorld() const {
    return qc::QC_3D;
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
  void setMinOfBothArrays ( const qc::ScalarArray<DataType, qc::QC_3D> &oldError ) {
    for ( int i = 0; i < oldError.getNumX(); ++i ) { // x-direction
      for ( int j = 0; j < oldError.getNumY(); ++j ) { // y-direction
        for ( int k = 0; k < oldError.getNumZ(); ++k ) { // z-direction
          if ( oldError.get ( i, j, k ) < this->get ( i, j, k ) ) {
            this->set ( i, j, k, oldError.get ( i, j, k ) );
          }
        }
      }
    }
  }

  //! method walks across the data-array and overwrites it
  //! with the saturated error-array
  void setErrorArray ( const qc::ScalarArray<DataType, qc::QC_3D> &error ) {
    for ( int i = 0; i < error.getNumX(); ++i ) { // x-direction
      for ( int j = 0; j < error.getNumY(); ++j ) { // y-direction
        for ( int k = 0; k < error.getNumZ(); ++k ) { // z-direction
          this->set ( i, j, k, error.get ( i, j, k ) );
        }
      }
    }
  }


  //! method walks across the data-array and overwrites it
  //! with the saturated error-array
  void makeSaturatedErrorArray ( qc::ScalarArray<DataType, qc::QC_3D>& );

  // auxiliary functions needed for saturation
  DataType maxChildNodeValue ( int X, int Y, int Z, int halfStep ) const {
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

  DataType minChildNodeValue ( int X, int Y, int Z, int halfStep ) const {
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

  inline int   checkElement ( const qc::Element &El ) const {
    return checkElement ( El, threshold );
  }

  inline int   checkElement ( const qc::Element &El, DataType Thres ) const {
    return checkElement ( El.x(), El.y(), El.z(), El.level(), Thres );
  }

  //! Returns boolean. Checks if the estimator is ok for the
  //! Element at X, Y, Z at Level according to Thres.
  int checkElement ( int /*X*/, int /*Y*/, int /*Z*/, int /*Level*/, DataType /*Thres*/ ) const {
    cerr << "has to be implemented" << endl;
    return 1;
  }

  inline int   checkElement ( int X, int Y, int Z, int Level ) const {
    return checkElement ( X, Y, Z, Level, threshold );
  }

  void reset() {
    this->setZero();
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

  qc::ScalarArray<DataType, qc::QC_3D> **elSatArrays;

  DataType     threshold;
};

#endif
