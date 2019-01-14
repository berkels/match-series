#ifndef __TIMESTEPPING_H
#define __TIMESTEPPING_H

#include <multiVector.h>
#include <vectorExtensions.h>

namespace aol {

/** Abstract basis class for time stepping schemes.
 *  \author Schwen (MEVIS)
 */
template< typename VectorType, typename TimeDepType >
class TimeStepBase {
public:
  typedef typename VectorType::RealType RealType;

protected:
  TimeStepBase ( const RealType TimeStep, const TimeDepType &TimeDep )
    : _timeStep ( TimeStep ), _timeDep ( TimeDep ) { }

  // default constructor and assignment operator do not work due to const ref member
  // synthesized copy constructor OK

  virtual ~TimeStepBase ( ) {
    // do nothing
  }

public:
  virtual void doTimeStep ( const VectorType &ProfileOld,
                            const RealType Time,
                            VectorType &ProfileNew ) const = 0;

  void doTimeStep ( const RealType Time, VectorType &Profile ) {
    VectorType profileTmp ( Profile, aol::STRUCT_COPY );
    this->doTimeStep ( Profile, Time, profileTmp );
    Profile = profileTmp;
  }

protected:
  RealType _timeStep;
  const TimeDepType  &_timeDep;
};


/** Abstract basis class for time stepping schemes with automatically adapted internal time step size.
 *  \author Schwen (MEVIS)
 */
template< typename VectorType, typename TimeDepType >
class TimeStepAdaptiveBase : public TimeStepBase < VectorType, TimeDepType > {
public:
  typedef typename VectorType::RealType RealType;

  TimeStepAdaptiveBase ( const RealType TimeStep, const TimeDepType &TimeDep, const RealType RelativeTolerance = 1.0e-6, const RealType IncreaseThreshold = 1.0e-2 )
    : TimeStepBase < VectorType, TimeDepType > ( TimeStep, TimeDep ),
      _lastInternalTimeStep ( TimeStep ), _relTol ( RelativeTolerance ), _incFac ( IncreaseThreshold ), _verbose ( false ), _internalTempVectors ( ), _internalTempVectorsInitialized ( false ) { }

  // copy constructor OK
  // assignment not possible due to const reference members

  void setTimeStep ( const RealType TimeStepNew ) {
    this->_timeStep = TimeStepNew;
  }

  void setVerbose ( const bool Verbose ) {
    _verbose = Verbose;
  }

  void setTimeStepAndLastTimeStep ( const RealType TimeStepNew ) {
    this->_timeStep = TimeStepNew;
    _lastInternalTimeStep = TimeStepNew;
  }

  RealType getLastInternalTimeStep ( ) const {
    return ( _lastInternalTimeStep );
  }

protected:
  //! to be implemented in derived classes
  virtual RealType errorForTimestep ( const VectorType &ProfileOld,
                                      const RealType Time,
                                      VectorType &ProfileNew,
                                      const RealType TimeStep,
                                      aol::RandomAccessContainer< VectorType > &InternalTempVectors ) const = 0;

public:
  using aol::TimeStepBase<VectorType, TimeDepType>::doTimeStep;

  void doTimeStep ( const VectorType &ProfileOld,
                    const RealType Time,
                    VectorType &ProfileNew ) const {

    VectorType updateOld ( ProfileOld, aol::DEEP_COPY ), updateTmp ( ProfileOld, aol::STRUCT_COPY_UNINIT );
    if ( ( ! _internalTempVectorsInitialized ) || ( _internalTempVectors[0].compareDim ( ProfileOld ) == true ) ) {
      _internalTempVectors.clear();
      for ( int i = 0; i < 9; ++i ) {
        _internalTempVectors.constructDatumAndPushBack ( ProfileOld, aol::STRUCT_COPY_UNINIT );
      }
      _internalTempVectorsInitialized = true;
    }

    uint64_t curtUN = 0, tUN = 1, tUD = std::max ( static_cast<uint64_t> ( floor ( ( this->_timeStep / _lastInternalTimeStep ) + 0.5 ) ), static_cast<uint64_t> ( 1 ) ); // numerator and denominator for internal time step
    RealType curInternalTime = 0.0;
    RealType curInternalTimeStep = aol::NumberTrait<RealType>::NaN;

    while ( curtUN < tUD ) {
      tUN = std::min ( tUN, tUD - curtUN ); // don't go beyond _timeStep
      curInternalTimeStep = ( tUN * this->_timeStep ) / tUD;
      RealType curError = this->errorForTimestep ( updateOld, Time + curInternalTime, updateTmp, curInternalTimeStep, _internalTempVectors );

      const RealType tolerance = std::min ( aol::NumberTrait<RealType>::one, updateTmp.norm() ) * _relTol; // this mix of relative and absolute tolerance can probably be improved

      if ( _verbose ) {
        cerr << ". initially  " << aol::detailedFormat ( curInternalTimeStep ) << " at " << aol::detailedFormat ( curInternalTime ) << " " << aol::detailedFormat ( curError ) << " " << aol::detailedFormat ( tolerance ) << endl;
      }

      // if error sufficiently small, increase time step
      while ( ( curError < ( _incFac * tolerance ) ) && ( ( curtUN + tUN ) < tUD ) ) {
        if ( tUN >= ( static_cast<uint64_t> ( 1 ) << 62 ) ) {
          throw aol::Exception ( "TimeStepAdaptiveBase::doTimeStep: integer overflow expected", __FILE__, __LINE__ );
        }
        tUN *= 2;
        tUN = std::min ( tUN, tUD - curtUN ); // don't go beyond _timeStep
        curInternalTimeStep = ( tUN * this->_timeStep ) / tUD;
        curError = this->errorForTimestep ( updateOld, Time + curInternalTime, updateTmp, curInternalTimeStep, _internalTempVectors );
        if ( _verbose ) {
          cerr << "+ increasing " << aol::detailedFormat ( curInternalTimeStep ) << " " << tUD << " " << tUN << " " << curtUN << " at " << aol::detailedFormat ( curInternalTime ) << " " << aol::detailedFormat ( curError ) << " " << aol::detailedFormat ( tolerance ) << endl;
        }
      }

      // if error too large, decrease time step
      while ( curError > tolerance ) {
        if ( ( tUD >= ( static_cast<uint64_t> ( 1 ) << 62 ) ) || ( curtUN >= ( static_cast<uint64_t> ( 1 ) << 62 ) ) ){
          throw aol::Exception ( "TimeStepAdaptiveBase::doTimeStep: integer overflow expected", __FILE__, __LINE__ );
        }
        tUD *= 2;
        curtUN *= 2; // ratio curtUN / tUD not changed.
        // adaption of tUN not necessary
        curInternalTimeStep = ( tUN * this->_timeStep ) / tUD;
        curError = this->errorForTimestep ( updateOld, Time + curInternalTime, updateTmp, curInternalTimeStep, _internalTempVectors );
        if ( _verbose ) {
          cerr << "- decreasing " << aol::detailedFormat ( curInternalTimeStep ) << " " << tUD << " " << tUN << " " << curtUN << " at " << aol::detailedFormat ( curInternalTime ) << " " << aol::detailedFormat ( curError ) << " " << aol::detailedFormat ( tolerance ) << endl;
        }
      }

      // If all numbers are even, we can cancel common factors 2. Otherwise, integer overflows may occur.
      while ( ( tUN % 2 == 0 ) && ( curtUN % 2 == 0 ) && ( tUD % 2 == 0 ) ) {
        tUN /= 2;
        curtUN /= 2;
        tUD /= 2;
        if ( _verbose ) {
          cerr << "/ dividing" << endl;
        }
      }

      // this should be an appropriate time step now. update things accordingly.
      curtUN += tUN;
      curInternalTime = ( curtUN * this->_timeStep ) / tUD;
      updateOld = updateTmp;
    }
    if ( _verbose ) {
      cerr << "went until " << aol::detailedFormat ( curInternalTime ) << " of " << aol::detailedFormat ( this->_timeStep ) << endl;
    }
    _lastInternalTimeStep = curInternalTimeStep;
    ProfileNew = updateOld;
  }


protected:
  mutable RealType _lastInternalTimeStep;
  const RealType _relTol, _incFac;
  bool _verbose;
  mutable aol::RandomAccessContainer< VectorType > _internalTempVectors;
  mutable bool _internalTempVectorsInitialized;
};


/** Runge-Kutta-Fehlbarg 4th/5th order time stepping (automatically adapting time step),
 *  implemented according to E. Fehlberg: "Klassische Runge-Kutta-Formeln vierter und niedrigerer Ordnung mit Schrittweiten-Kontrolle und ihre Anwendung auf Wärmeleitungsprobleme", Computing 6, pp. 61-70, 1970. DOI 10.1007/BF02241732.
 *  \author Schwen (MEVIS)
 */
template< typename VectorType, typename TimeDepType >
class TimeStepRKF45 : public TimeStepAdaptiveBase < VectorType, TimeDepType > {
public:
  typedef typename VectorType::RealType RealType;

  TimeStepRKF45 ( const RealType TimeStep, const TimeDepType &Reaction, const RealType RelativeTolerance = 1.0e-6, const RealType IncreaseThreshold = 1.0e-2 )
    : TimeStepAdaptiveBase < VectorType, TimeDepType > ( TimeStep, Reaction, RelativeTolerance, IncreaseThreshold ) {
  }

protected:
  /* virtual */
  RealType errorForTimestep ( const VectorType &ProfileOld,
                              const RealType Time,
                              VectorType &ProfileNew,
                              const RealType TimeStep,
                              aol::RandomAccessContainer< VectorType > &InternalTempVectors ) const {

    VectorType
      &PI = InternalTempVectors[0],
      &K0 = InternalTempVectors[1],
      &K1 = InternalTempVectors[2],
      &K2 = InternalTempVectors[3],
      &K3 = InternalTempVectors[4],
      &K4 = InternalTempVectors[5],
      &K5 = InternalTempVectors[6],
      &DT4 = ProfileNew,
      &DT5 = InternalTempVectors[7],
      &difference = InternalTempVectors[8];
    RealType t = 0.0;

    // constants in this method are those from the paper cited above

    PI = ProfileOld;
    t = Time;
    this->_timeDep.evaluateDerivative ( PI, t, K0 );
    K0 *= TimeStep;

    setLinComb1 ( PI, ProfileOld, 1./4., K0 );
    t = Time + 1. / 4. * TimeStep;
    this->_timeDep.evaluateDerivative ( PI, t, K1 );
    K1 *= TimeStep;

    setLinComb2 ( PI, ProfileOld, 3./32., K0, 9./32., K1 );
    t = Time + 3. / 8. * TimeStep;
    this->_timeDep.evaluateDerivative ( PI, t, K2 );
    K2 *= TimeStep;

    setLinComb3 ( PI, ProfileOld, 1932./2197., K0, -7200./2197., K1, 7296./2197., K2 );
    t = Time + 12. / 13. * TimeStep;
    this->_timeDep.evaluateDerivative ( PI, t, K3 );
    K3 *= TimeStep;

    setLinComb4 ( PI, ProfileOld, 439./216., K0, -8., K1, 3680./513., K2, -845./4104., K3 );
    t = Time + TimeStep;
    this->_timeDep.evaluateDerivative ( PI, t, K4 );
    K4 *= TimeStep;

    setLinComb5 ( PI, ProfileOld, -8./27., K0, 2., K1, -3544./2565., K2, 1859./4104., K3, -11./40., K4 );
    t = Time + 1. / 2. * TimeStep;
    this->_timeDep.evaluateDerivative ( PI, t, K5 );
    K5 *= TimeStep;

    DT4 = ProfileOld; DT5 = ProfileOld;
    setLinComb4 ( DT4, ProfileOld, 25./216., K0, 1408./2565., K2, 2197./4104., K3, -1./5., K4 );
    setLinComb5 ( DT5, ProfileOld, 16./135., K0, 6656./12825., K2, 28561./56430., K3, -9./50., K4, 2./55., K5 );

    difference = DT4;
    difference -= DT5;
    return ( ( difference.norm() / difference.getTotalSize() ) / TimeStep );
  }

protected:
  static inline void doAddMultiple ( aol::Vector<RealType> &Dest, const aol::Vector<RealType> &Arg, const RealType Factor ) {
    Dest.addMultiple ( Arg, Factor );
  }

  static inline void doAddMultiple ( aol::MultiVector<RealType> &Dest, const aol::MultiVector<RealType> &Arg, const RealType Factor ) {
    Dest.addMultipleParallel ( Arg, Factor );
  }


  static inline void setLinComb1 ( aol::Vector<RealType> &Dest, const aol::Vector<RealType> &Arg0, const RealType Fac1, const aol::Vector<RealType> &Arg1 ) {
    for ( int i = 0; i < Dest.size(); ++i ) {
      Dest[i] = Arg0[i] + Fac1 * Arg1[i];
    }
  }

  static inline void setLinComb2 ( aol::Vector<RealType> &Dest, const aol::Vector<RealType> &Arg0, const RealType Fac1, const aol::Vector<RealType> &Arg1, const RealType Fac2, const aol::Vector<RealType> &Arg2 ) {
    for ( int i = 0; i < Dest.size(); ++i ) {
      Dest[i] = Arg0[i] + Fac1 * Arg1[i] + Fac2 * Arg2[i];
    }
  }

  static inline void setLinComb3 ( aol::Vector<RealType> &Dest, const aol::Vector<RealType> &Arg0, const RealType Fac1, const aol::Vector<RealType> &Arg1, const RealType Fac2, const aol::Vector<RealType> &Arg2,
                                   const RealType Fac3, const aol::Vector<RealType> &Arg3 ) {
    for ( int i = 0; i < Dest.size(); ++i ) {
      Dest[i] = Arg0[i] + Fac1 * Arg1[i] + Fac2 * Arg2[i] + Fac3 * Arg3[i];
    }
  }

  static inline void setLinComb4 ( aol::Vector<RealType> &Dest, const aol::Vector<RealType> &Arg0, const RealType Fac1, const aol::Vector<RealType> &Arg1, const RealType Fac2, const aol::Vector<RealType> &Arg2,
                                   const RealType Fac3, const aol::Vector<RealType> &Arg3, const RealType Fac4, const aol::Vector<RealType> &Arg4 ) {
    for ( int i = 0; i < Dest.size(); ++i ) {
      Dest[i] = Arg0[i] + Fac1 * Arg1[i] + Fac2 * Arg2[i] + Fac3 * Arg3[i] + Fac4 * Arg4[i];
    }
  }

  static inline void setLinComb5 ( aol::Vector<RealType> &Dest, const aol::Vector<RealType> &Arg0, const RealType Fac1, const aol::Vector<RealType> &Arg1, const RealType Fac2, const aol::Vector<RealType> &Arg2,
                            const RealType Fac3, const aol::Vector<RealType> &Arg3, const RealType Fac4, const aol::Vector<RealType> &Arg4, const RealType Fac5, const aol::Vector<RealType> &Arg5 ) {
    for ( int i = 0; i < Dest.size(); ++i ) {
      Dest[i] = Arg0[i] + Fac1 * Arg1[i] + Fac2 * Arg2[i] + Fac3 * Arg3[i] + Fac4 * Arg4[i] + Fac5 * Arg5[i];
    }
  }


  static inline void setLinComb1 ( aol::MultiVector<RealType> &Dest, const aol::MultiVector<RealType> &Arg0, const RealType Fac1, const aol::MultiVector<RealType> &Arg1 ) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < Dest.numComponents(); ++i ) {
      for ( int j = 0; j < Dest[i].size(); ++j ) {
        Dest[i][j] = Arg0[i][j] + Fac1 * Arg1[i][j];
      }
    }
  }

  static inline void setLinComb2 ( aol::MultiVector<RealType> &Dest, const aol::MultiVector<RealType> &Arg0, const RealType Fac1, const aol::MultiVector<RealType> &Arg1, const RealType Fac2, const aol::MultiVector<RealType> &Arg2 ) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < Dest.numComponents(); ++i ) {
      for ( int j = 0; j < Dest[i].size(); ++j ) {
        Dest[i][j] = Arg0[i][j] + Fac1 * Arg1[i][j] + Fac2 * Arg2[i][j];
      }
    }
  }

  static inline void setLinComb3 ( aol::MultiVector<RealType> &Dest, const aol::MultiVector<RealType> &Arg0, const RealType Fac1, const aol::MultiVector<RealType> &Arg1, const RealType Fac2, const aol::MultiVector<RealType> &Arg2,
                                   const RealType Fac3, const aol::MultiVector<RealType> &Arg3 ) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < Dest.numComponents(); ++i ) {
      for ( int j = 0; j < Dest[i].size(); ++j ) {
        Dest[i][j] = Arg0[i][j] + Fac1 * Arg1[i][j] + Fac2 * Arg2[i][j] + Fac3 * Arg3[i][j];
      }
    }
  }

  static inline void setLinComb4 ( aol::MultiVector<RealType> &Dest, const aol::MultiVector<RealType> &Arg0, const RealType Fac1, const aol::MultiVector<RealType> &Arg1, const RealType Fac2, const aol::MultiVector<RealType> &Arg2,
                                   const RealType Fac3, const aol::MultiVector<RealType> &Arg3, const RealType Fac4, const aol::MultiVector<RealType> &Arg4 ) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < Dest.numComponents(); ++i ) {
      for ( int j = 0; j < Dest[i].size(); ++j ) {
        Dest[i][j] = Arg0[i][j] + Fac1 * Arg1[i][j] + Fac2 * Arg2[i][j] + Fac3 * Arg3[i][j] + Fac4 * Arg4[i][j];
      }
    }
  }

  static inline void setLinComb5 ( aol::MultiVector<RealType> &Dest, const aol::MultiVector<RealType> &Arg0, const RealType Fac1, const aol::MultiVector<RealType> &Arg1, const RealType Fac2, const aol::MultiVector<RealType> &Arg2,
                                   const RealType Fac3, const aol::MultiVector<RealType> &Arg3, const RealType Fac4, const aol::MultiVector<RealType> &Arg4, const RealType Fac5, const aol::MultiVector<RealType> &Arg5 ) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < Dest.numComponents(); ++i ) {
      for ( int j = 0; j < Dest[i].size(); ++j ) {
        Dest[i][j] = Arg0[i][j] + Fac1 * Arg1[i][j] + Fac2 * Arg2[i][j] + Fac3 * Arg3[i][j] + Fac4 * Arg4[i][j] + Fac5 * Arg5[i][j];
      }
    }
  }

public:
  /** Class used in the selfTest
   *  Test example from Granville Sewell: The Numerical Solution of Ordinary and Partial Differential Equations, Academic Press 1998, ISBN 0-12-637475-9, p. 62/63
   */
  struct RKFTestReaction {
    void evaluateDerivative ( const aol::Vector<RealType> &/*ProfileOld*/,
                              const RealType Time,
                              aol::Vector<RealType> &ProfileTDeriv ) const {
      ProfileTDeriv[0] = 100.0 / ( 1.0 + 10000 * aol::Sqr ( Time ) );
    }
  };

  static bool isSelfTestOK ( ) {
    aol::Vector<double> uCur ( 1 );
    uCur[0] = atan ( -100.0 );

    const double tau = 0.25; // time step

    RKFTestReaction reac;
    TimeStepRKF45< aol::Vector<double>, RKFTestReaction > ts ( tau, reac );

    for ( int t = -4; t < 4; ++t ) {
      ts.doTimeStep ( t * tau, uCur );
      // cerr << aol::shortFormat ( ( t + 1 ) * tau ) << " " << aol::detailedFormat ( uCur[0] ) << endl;
    }

    return ( fabs ( uCur[0] - atan ( 100.0 ) ) < 1.0e-6 );
  }

};

}

#endif
