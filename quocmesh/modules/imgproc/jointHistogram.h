#ifndef __JOINTHISTOGRAM_H
#define __JOINTHISTOGRAM_H

#include <FEOpInterface.h>
#include <linearSmoothOp.h>

namespace im {

/**
 * Computes the joint probability of two images. The input images have to be scaled to [0,1]
 *
 * \author Berkels, Han
 */
template <typename RealType>
class JointHistogram {
protected:
  //! Reference image.
  const qc::Array<RealType> &_r;

  //! Template image.
  const qc::Array<RealType> &_t;

  //! The level of bin in each dimension. Number of bin = 2^_intensityLevel;
  const int _intensityLevel;
  const int _numberOfIntensityValues;
  //! The standard variance of parzen windows (Gaussian windows)
  const RealType _beta;

  //! How many percents of grid points are sampled
  const RealType _percentSample;

  //! The pointer that refers to the 2D grid of histogram.
  const qc::GridDefinition _histoGrid;

  bool _jointHistogramIsComputed;
  bool _mutualInformationIsComputed;

  //! table of joint histgram
  qc::ScalarArray<RealType, qc::QC_2D> _histoTable;
  qc::ScalarArray<RealType, qc::QC_2D> _mutualInformation;

  //! Sampling data from images. Currently, Sampling covers all the voxels in images.
  //! \todo Random sampling should be implemented
  void sampleData();

  //! Smoothing the histogram
  void smooth ( qc::ScalarArray<RealType, qc::QC_2D> &Array );

public:
  //! \warning The values of Reference and Template are supposed to be in [0,1]. If they are not, the class crashes when sampleData() is called!
  JointHistogram ( const qc::Array<RealType> &Reference,
                   const qc::Array<RealType> &Template,
                   const int IntensityLevel,
                   const RealType Beta = 1.0,
                   const RealType PercentSample = 1.0 ) :
      _r ( Reference ),
      _t ( Template ),
      _intensityLevel ( IntensityLevel ),
      _numberOfIntensityValues ( static_cast<int> ( pow ( ( 2.0 ), _intensityLevel ) + 1 ) ),
      _beta ( Beta ),
      _percentSample ( PercentSample ),
      _histoGrid ( _intensityLevel, qc::QC_2D ),
      _jointHistogramIsComputed ( false ),
      _mutualInformationIsComputed ( false ),
      _histoTable ( _histoGrid ),
      _mutualInformation ( _histoGrid ) {
    // check
    if ( IntensityLevel < 1 ) {
      throw aol::Exception ( "IntensityLevel needs to be one or bigger", __FILE__, __LINE__ );
    }
    //cerr << "_intensityLevel = " << _intensityLevel << " " << "_numberOfIntensityValues = " << _numberOfIntensityValues << endl;
  }
  virtual ~JointHistogram() {}

  /** @brief Compute the marginal probablities
      @param  fr    [in] The flag of marginalization. 0: reference. 1: template.
      @return           The array of marginal probabilities
   */
  void marginalizeVector ( aol::Vector<RealType> &MarProb, int fr );

  virtual void computeJointHistogram();

  //! Compute the joint histogram with Parzen Windowing.
  //! \note Reference implementation for testing purposes only. Should not be used, because it's REALLY slow!
  void computeJointHistogramWithParzenWindowing();

  /** @brief Transform the joint probablity to L. This function is important for \
      MI based registration. The example is given in MI_reg.cpp.
      1) A joint histogram object is initiated like:
             JointHistogram<ConfType,qc::QC_2D> Joint(grid,r0,t,level,beta);
      2) The joint histogram is computed:
             Joint.computeJointHistogram();
      3) The computation of MI information is done as:
             Joint.computeMutualInformation();
   */
  virtual void computeMutualInformation();

  RealType computeEntropyOfJointHistogram();

  void marginalize ( qc::ScalarArray<RealType, qc::QC_2D> &MargTable, int fr ) {
    // check
    if ( fr < 0 || fr > 1 ) {
      throw aol::Exception ( "fr needs to be equal to zero or one", __FILE__, __LINE__ );
    }

    // define the marginal histogram table
    MargTable.setZero();

    if ( fr == 0 ) {
      // Pr
      // store the marginal probabities on first column
      for ( int i = 0; i < _numberOfIntensityValues; i++ ) {
        for ( int j = 0; j < _numberOfIntensityValues; j++ ) {
          MargTable.add ( i, 0, _histoTable.get ( i, j ) );
        }
      }
      // copy to other column
      for ( int i = 0; i < _numberOfIntensityValues; i++ ) {
        for ( int j = 1; j < _numberOfIntensityValues; j++ ) {
          MargTable.set ( i, j, MargTable.get ( i, 0 ) );
        }
      }
    } else {
      // t marginal
      // store the marginal probabities on first row
      for ( int i = 0; i < _numberOfIntensityValues; i++ ) {
        for ( int j = 0; j < _numberOfIntensityValues; j++ ) {
          MargTable.add ( 0, j, _histoTable.get ( i, j ) );
        }
      }
      // copy to other rows
      for ( int i = 1; i < _numberOfIntensityValues; i++ ) {
        for ( int j = 0; j < _numberOfIntensityValues; j++ ) {
          MargTable.set ( i, j, MargTable.get ( 0, j ) );
        }
      }
    }
  }

  int getNumberOfIntensityValues() const {
    return _numberOfIntensityValues;
  }

  const qc::ScalarArray<RealType, qc::QC_2D>& getJointHistogram() {
    if ( _jointHistogramIsComputed == false )
      this->computeJointHistogram();
    const qc::ScalarArray<RealType, qc::QC_2D> &temp = _histoTable;
    return temp;
  }
  const qc::ScalarArray<RealType, qc::QC_2D>& getMutualInformation() {
    if ( _mutualInformationIsComputed == false )
      this->computeMutualInformation();
    const qc::ScalarArray<RealType, qc::QC_2D> &temp = _mutualInformation;
    return temp;
  }
  const qc::GridDefinition& getHistoGrid() const {
    const qc::GridDefinition &temp = _histoGrid;
    return temp;
  }
  const aol::Vector<RealType>& getReference() const {
    const aol::Vector<RealType> &temp = _r;
    return temp;
  }
  const aol::Vector<RealType>& getTemplate() const {
    const aol::Vector<RealType> &temp = _t;
    return temp;
  }
private:
  RealType temp ( const int i, const int j, const RealType MPj, const RealType dMP ) const {
    RealType value = 0.;
    if ( MPj != 0. )
      value = dMP / MPj;
    if ( _histoTable.get ( i, j ) != 0. )
      value -= _histoTable.dyFD ( i, j ) / _histoTable.get ( i, j );
    value *= ( _numberOfIntensityValues - 1 );
    return value;
  }

};
template <typename RealType>
RealType JointHistogram<RealType>::computeEntropyOfJointHistogram() {
  if ( _jointHistogramIsComputed == false )
    this->computeJointHistogram();
  RealType entropy = 0.;
  RealType temp = 0.;
  for ( int i = 0; i < _numberOfIntensityValues; i++ ) {
    for ( int j = 0; j < _numberOfIntensityValues; j++ ) {
      temp = _histoTable.get ( i, j );
      if ( temp != 0 )
        entropy -= temp * log ( temp );
    }
  }
  return entropy;
}
template <typename RealType>
void JointHistogram<RealType>::computeMutualInformation() {
  if ( _jointHistogramIsComputed == false )
    this->computeJointHistogram();
  aol::Vector<RealType> MP ( _numberOfIntensityValues );
  this->marginalizeVector ( MP, 1 );
  // compute the L

  _mutualInformation.setZero();

  RealType dMP = 0.;

  for ( int j = 0; j < 1; j++ ) {
    dMP = MP[j+1] - MP[j];
    for ( int i = 0; i < _numberOfIntensityValues; i++ ) {
      _mutualInformation.set ( i, j, temp ( i, j, MP[j], dMP ) );
    }
  }
  for ( int j = 1; j < _numberOfIntensityValues - 1; j++ ) {
    dMP = 0.5 * ( MP[j+1] - MP[j-1] );
    for ( int i = 0; i < _numberOfIntensityValues; i++ ) {
      _mutualInformation.set ( i, j, temp ( i, j, MP[j], dMP ) );
    }
  }
  for ( int j = _numberOfIntensityValues - 1;j < _numberOfIntensityValues; j++ ) {
    dMP = MP[j] - MP[j-1];
    for ( int i = 0; i < _numberOfIntensityValues; i++ ) {
      _mutualInformation.set ( i, j, temp ( i, j, MP[j], dMP ) );
    }
  }
  // smoothing
  this->smooth ( _mutualInformation );
  _mutualInformationIsComputed = true;
}
template <typename RealType>
void JointHistogram<RealType>::marginalizeVector ( aol::Vector<RealType> &MarProb, int fr ) {
  if ( _jointHistogramIsComputed == false )
    this->computeJointHistogram();

  // marginal probability
  MarProb.setZero();

  // check
  if ( fr < 0 || fr > 1 ) {
    throw aol::Exception ( "fr needs to be equal to zero or one", __FILE__, __LINE__ );
  }
  if ( fr == 0 ) {
    // r marginal
    for ( int i = 0; i < _numberOfIntensityValues; i++ ) {
      for ( int j = 0; j < _numberOfIntensityValues; j++ ) {
        MarProb[i] +=  _histoTable.get ( i, j );
      }
    }
  } else {
    // t marginal
    for ( int i = 0 ; i < _numberOfIntensityValues; i++ ) {
      for ( int j = 0; j < _numberOfIntensityValues; j++ ) {
        MarProb[j] +=  _histoTable.get ( i, j );
      }
    }
  }
}

template <typename RealType>
void JointHistogram<RealType>::computeJointHistogram() {
  // build up discrete histogram
  this->sampleData();

  // smooth the histogram
  this->smooth ( _histoTable );

  // normalization
  _histoTable /= _histoTable.sum();

  _jointHistogramIsComputed = true;
}

template <typename RealType>
void JointHistogram<RealType>::computeJointHistogramWithParzenWindowing() {
  cerr << "Reference implementation for testing purposes only!\n";
  const RealType betaScaled = _beta * _histoGrid.H() * _histoGrid.H();
  const RealType scaling = 1. / ( 2 * aol::NumberTrait<RealType>::pi * betaScaled * ( _r.getNumX() - 1. ) * ( _r.getNumY() - 1. ) * aol::Max ( 1., _r.getNumZ() - 1. ) );
  for ( int i = 0; i < _numberOfIntensityValues; i++ ) {
    for ( int j = 0; j < _numberOfIntensityValues; j++ ) {
      const RealType iScaled = static_cast<RealType> ( i ) / static_cast<RealType> ( _numberOfIntensityValues - 1 );
      const RealType jScaled = static_cast<RealType> ( j ) / static_cast<RealType> ( _numberOfIntensityValues - 1 );
      RealType integral = 0.;
      for ( int x = 0; x < _r.getNumX(); x++ ) {
        for ( int y = 0; y < _r.getNumY(); y++ ) {
          for ( int z = 0; z < _r.getNumZ(); z++ ) {
            if ( ( x >= _r.getNumX() ) || ( y >= _r.getNumY() ) || ( z >= _r.getNumZ() ) )
              cerr << " Obacht! ";
            integral += exp ( -1. * ( aol::Sqr ( _r.get ( x, y, z ) - iScaled ) + aol::Sqr ( _t.get ( x, y, z ) - jScaled ) ) / ( 2 * betaScaled ) );
          }
        }
      }
      _histoTable.set ( i, j, scaling*integral );
    }
  }
  // normalization
  _histoTable /= _histoTable.sum();

  _jointHistogramIsComputed = true;
}
template <typename RealType>
void JointHistogram<RealType>::smooth ( qc::ScalarArray<RealType, qc::QC_2D> &Array ) {
  qc::LinearSmoothOp<RealType> linSmooth;
  linSmooth.setCurrentGrid ( _histoGrid );
  linSmooth.setSigma ( _beta * ( _histoGrid.H() ) );
  linSmooth.apply ( Array ,  Array );
  return;
}

template <typename RealType>
void JointHistogram<RealType>::sampleData() {
  //clear histogram table
  _histoTable.setZero();

  // length of vector
  const int length = _r.size();

  // fill the histogram table
  for ( int i = 0; i < length; ++i ) {
    // The gray values are supposed to be in [0,1]. If they are not, this will crash!
    _histoTable.add ( static_cast<int> ( ( _numberOfIntensityValues - 1 ) * _r[i] ), static_cast<int> ( ( _numberOfIntensityValues - 1 ) * _t[i] ), 1.0 );
  }
}

} // end namespace qc

#endif // __JOINTHISTOGRAM_H
