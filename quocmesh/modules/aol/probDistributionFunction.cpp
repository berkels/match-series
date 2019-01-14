#include "probDistributionFunction.h"

#include <arrayExtensions.h>
#include <discreteFunction.h>

namespace aol {

template< typename RealType >
RealType ProbDistFuncHelper<RealType>::KolmogorovProb ( const RealType z ) {
  // code adapted from the CERN root library TMash::KolmogorovProb(Double_t z)
  const RealType u = fabs ( z );

  if ( u < 0.2 ) {
    return ( 1.0 );
  } else if ( u < 0.755 ) {
    const RealType
      v = 1.0 / aol::Sqr ( u ),
      w = 2.50662827,
      c1 = - aol::Sqr ( aol::NumberTrait<RealType>::pi ) / 8.0,
      c2 = 9.0 * c1,
      c3 = 25.0 * c1;
    return ( 1.0 - w * ( exp ( c1*v ) + exp ( c2*v ) + exp ( c3*v ) ) / u );
  } else if ( u < 6.8116 ) {
    const RealType
      fj[4] = { -2.0, -8.0, -18.0, -32.0 },
      v = aol::Sqr ( u );
      RealType r[4] = {0.0, 0.0, 0.0, 0.0 };
      int maxj = aol::Max ( 1, static_cast<int> ( 3.0 / u + 0.5 ) );
      for ( int j = 0; j < maxj; ++j ) {
        r[j] = exp ( fj[j] * v );
      }
      return ( 2.0 * ( r[0] - r[1] + r[2] - r[3] ) );
  } else {
    return ( 0.0 );
  }
}


template< typename RealType >
RealType ProbDistFuncHelper<RealType>::KolmogorovProbTwoSmallSamples ( const RealType X, const unsigned int N0, const unsigned int N1 ) {
  RealType ret = 0.0;
  if ( ! aol::isNaN ( X ) ) {

    const unsigned int
      m = std::min ( N0, N1 ),
      n = std::max ( N0, N1 );

    const long double
      md = m,
      oneOverM = 1.0 / md,
      nd = n,
      oneOverN = 1.0 / nd,
      q = ( 0.5 + floor ( X * m * n - 1.0e-10 ) ) / ( m * n );

    aol::Vector<long double> u ( n+1 );
    for( unsigned int j = 0; j <= n; ++j ) {
      u[j] = ( ( ( j * oneOverN ) > q ) ? 0.0 : 1.0 );
    }

    for ( unsigned int i = 1; i <= m; ++i ) {
      u[0] = ( ( ( i * oneOverM ) > q ) ? 0.0 : u[0] );
      for( unsigned int j = 1; j <= n; ++j ) {
        u[j] = ( ( fabs( i * oneOverM - j * oneOverN ) > q ) ? 0.0 : ( u[j] + u[j-1] ) );
      }
      u *= i / ( i + nd );
    }
    ret = static_cast<RealType> ( 1.0 - u[n] );

  } else {
    ret = aol::NumberTrait<RealType>::NaN;
  }
  return ( ret );
}


template< typename RealType >
RealType ProbDistFuncHelper<RealType>::CramerVonMisesProb ( const RealType z, const unsigned int n0, const unsigned int n1 ) {
  // code adapted from a matlab script by Juan Cardelino (cmtest2.m); values from Anderson, Darling: Asymptotic theory of certain "Goodness of fit" criteria based on stochastic processes, The Annals of Mathematical Statistics 23, 1952.
  RealType zValArray[101] = { 0.00000, 0.02480, 0.02878, 0.03177, 0.03430, 0.03656, 0.03865, 0.04061, 0.04247, 0.04427,
                              0.04601, 0.04772, 0.04939, 0.05103, 0.05265, 0.05426, 0.05586, 0.05746, 0.05904, 0.06063,
                              0.06222, 0.06381, 0.06541, 0.06702, 0.06863, 0.07025, 0.07189, 0.07354, 0.07521, 0.07690,
                              0.07860, 0.08032, 0.08206, 0.08383, 0.08562, 0.08744, 0.08928, 0.09115, 0.09306, 0.09499,
                              0.09696, 0.09896, 0.10100, 0.10308, 0.10520, 0.10736, 0.10956, 0.11182, 0.11412, 0.11647,
                              0.11888, 0.12134, 0.12387, 0.12646, 0.12911, 0.13183, 0.13463, 0.13751, 0.14046, 0.14350,
                              0.14663, 0.14986, 0.15319, 0.15663, 0.16018, 0.16385, 0.16765, 0.17159, 0.17568, 0.17992,
                              0.18433, 0.18892, 0.19371, 0.19870, 0.20392, 0.20939, 0.21512, 0.22114, 0.22748, 0.23417,
                              0.24124, 0.24874, 0.25670, 0.26520, 0.27429, 0.28406, 0.29460, 0.30603, 0.31849, 0.33217,
                              0.34730, 0.36421, 0.38331, 0.40520, 0.43077, 0.46136, 0.49929, 0.54885, 0.61981, 0.74346, 1.16786 };
  // At least under MinGW zValArray is not guaranteed to be 16 byte aligned and thus
  // can't be used as FLAT_COPY in a Vector.
  aol::Vector<RealType> zVal ( zValArray, 101, aol::DEEP_COPY ), pVal ( 101 );
  for ( short i = 0; i < 100; ++i ) {
    pVal[i] = 0.01 * i;
  }
  pVal[100] = 0.999;

  aol::DiscreteValueInterpolator< RealType, RealType > interp ( zVal, pVal );

  const RealType
    n0d = static_cast<RealType> ( n0 ),
    n1d = static_cast<RealType> ( n1 ), // we need to prevent integer overflows!
    meanT = 1.0 / 6.0 + 1.0 / ( 6.0 * ( n0d + n1d ) ),
    varTFactor0 = ( 1.0 / 45.0 ) * ( ( n0d + n1d + 1.0 ) / ( aol::Sqr ( n0d + n1d ) ) ),
    varTFactor1 = ( 4.0 * n0d * n1d * ( n0d + n1d ) - 3.0 * ( aol::Sqr ( n0d ) + aol::Sqr ( n1d ) ) - 2 * n0d * n1d ) / ( 4.0 * n0d * n1d ),
    varT = varTFactor0 * varTFactor1,
    CvMAdapted = ( z - meanT ) / sqrt ( 45.0 * varT ) + 1.0 / 6.0;

  RealType pValue = aol::NumberTrait<RealType>::NaN;

  if ( CvMAdapted < zVal[0] ) {
    pValue = 0.0;
  } else if ( CvMAdapted < zVal[100] ) {
    pValue = interp.evaluate ( CvMAdapted );
  } else {
    pValue = 1.0;
  }

  return ( 1.0 - pValue );
}


template< typename RealType >
ProbDistributionFunctionAnyD<RealType>::ProbDistributionFunctionAnyD ( ) : _nSamples ( 0 ), _L2Dist ( aol::NumberTrait<RealType>::NaN ), _LInfDist ( aol::NumberTrait<RealType>::NaN ), _CvMDist ( aol::NumberTrait<RealType>::NaN ) {
  // do nothing else
}


template< typename RealType >
RealType ProbDistributionFunctionAnyD<RealType>::getScaledKSDistanceTo ( const ProbDistributionFunctionAnyD &other ) const {
  const RealType factor = sqrt ( ( static_cast<RealType> ( this->_nSamples ) * static_cast<RealType> ( other._nSamples ) ) / (  this->_nSamples + other._nSamples ) );
  return ( factor * _LInfDist );
}



template< typename RealType >
ProbDistributionFunction1D<RealType>::ProbDistributionFunction1D ( ) {
  // do nothing else
}


template< typename RealType >
void ProbDistributionFunction1D<RealType>::initialize ( const std::map< RealType, unsigned int > &Histogram ) {
  unsigned int hSum = 0;
  for ( typename std::map< RealType, unsigned int >::const_iterator it = Histogram.begin(); it != Histogram.end(); ++it ) {
    hSum += it->second;
  }
  this->_nSamples = hSum;

  RealType sumSoFar = 0;
  for ( typename std::map< RealType, unsigned int >::const_iterator it = Histogram.begin(); it != Histogram.end(); ++it ) {
    sumSoFar += it->second / static_cast<RealType> ( hSum );
    _pdf[ it->first ] += sumSoFar;
  }
}


template< typename RealType >
void ProbDistributionFunction1D<RealType>::dump ( ostream& out ) const {
  for ( typename std::map< RealType, RealType >::const_iterator it = _pdf.begin(); it != _pdf.end(); ++it ) {
    out << aol::longScientificFormat ( it->first ) << " " << aol::longScientificFormat ( it->second ) << endl;
  }
}


template< typename RealType >
RealType ProbDistributionFunction1D<RealType>::getDomainScaledL2DistanceTo ( const ProbDistributionFunction1D &other ) const {
  if ( this->_pdf.empty() || other._pdf.empty() ) {
    return ( aol::NumberTrait<RealType>::NaN );
  } else {
    const RealType
      minVal = aol::Min ( _pdf.begin()->first,  other._pdf.begin()->first  ),
      maxVal = aol::Max ( _pdf.rbegin()->first, other._pdf.rbegin()->first );
    return ( this->_L2Dist / sqrt ( maxVal - minVal ) );
  }
}


template< typename RealType >
RealType ProbDistributionFunction1D<RealType>::getScaledL2DistanceTo ( const ProbDistributionFunction1D &other ) const {
  const RealType factor = sqrt ( ( static_cast<RealType> ( this->_nSamples ) * static_cast<RealType> ( other._nSamples ) ) / (  this->_nSamples + other._nSamples ) );
  return ( factor * this->getDomainScaledL2DistanceTo ( other ) );
}


template< typename RealType >
RealType ProbDistributionFunction1D<RealType>::getScaledCvMDistanceTo ( const ProbDistributionFunction1D &other ) const {
  const RealType factor = ( static_cast<RealType> ( this->_nSamples ) * static_cast<RealType> ( other._nSamples ) ) / (  this->_nSamples + other._nSamples );
  return ( factor * this->_CvMDist );
}


template< typename RealType >
void ProbDistributionFunction1D<RealType>::computePDFdistTo ( const ProbDistributionFunction1D &other ) {
  if ( ! ( this->_pdf.empty() ) && ! ( other._pdf.empty() ) ) {

    std::map< RealType, PDFDiffStep1D > steps;
    {
      RealType valOld = 0;
      for ( typename std::map< RealType, RealType >::const_iterator it = this->_pdf.begin(); it != this->_pdf.end(); ++it ) {
        steps[ it->first ]._val[0]  = it->second;
        steps[ it->first ]._step[0] = it->second - valOld;
        valOld = it->second;
      }
      valOld = 0;

      for ( typename std::map< RealType, RealType >::const_iterator it = other._pdf.begin(); it != other._pdf.end(); ++it ) {
        steps[ it->first ]._val[1]  = it->second;
        steps[ it->first ]._step[1] = it->second - valOld;
        valOld = it->second;
      }
    }

#ifdef VERBOSE
    // dump (with holes)
    for ( typename std::map< RealType, PDFDiffStep1D >::const_iterator it = steps.begin(); it != steps.end(); ++it ) {
      cerr << it->first << ": " << it->second._val << ", " << it->second._step << endl;
    }
#endif

    // fill values
    aol::Vec2<RealType> oldVal;
    for ( typename std::map< RealType, PDFDiffStep1D >::iterator it = steps.begin(); it != steps.end(); ++it ) {
      for ( short d = 0; d < 2; ++d ) {
        if ( it->second._val[d] == 0 ) {
          it->second._val[d] = oldVal[d];
        }
        oldVal[d] = it->second._val[d];
      }
    }

#ifdef VERBOSE
    // dump
    for ( typename std::map< RealType, PDFDiffStep1D >::const_iterator it = steps.begin(); it != steps.end(); ++it ) {
      cerr << it->first << ": " << it->second._val << ", " << it->second._step << endl;
    }
#endif

    RealType l2sum = 0, cvmsum = 0, linfdist = -aol::NumberTrait<RealType>::Inf;
    typename std::map< RealType, PDFDiffStep1D >::const_iterator oldIt = steps.begin();
    for ( typename std::map< RealType, PDFDiffStep1D >::const_iterator it = steps.begin(); it != steps.end(); ++it ) {
      typename std::map< RealType, PDFDiffStep1D >::const_iterator nextIt = it; ++nextIt;

      if ( nextIt == steps.end() ) {
        if ( fabs ( it->second._val[1] - it->second._val[0] ) > 1.0e-4 ) { // THINK ABOUT THIS SOME MORE!
          cerr << __FILE__ << ": " << __LINE__ << ": difference detected: " << aol::longScientificFormat ( it->second._val[1] ) << " " << aol::longScientificFormat ( it->second._val[0] ) << endl;
        }
        continue; // and then for loop should be finished.
      }

      const RealType vdiff = it->second._val[1] - it->second._val[0];
      linfdist = aol::Max ( fabs ( vdiff ), linfdist );

#ifdef VERBOSE
      cerr << it->first << " " << vdiff << endl;
#endif

      {
        const RealType xdiff = nextIt->first - it->first;
        l2sum += aol::Sqr ( vdiff ) * xdiff;
      }

      cvmsum += it->second._step[0] * aol::Sqr ( vdiff ) * this->_nSamples;
      cvmsum += it->second._step[1] * aol::Sqr ( vdiff ) * other._nSamples;

      oldIt = it;
    }

    this->_L2Dist = sqrt ( l2sum );
    this->_LInfDist = linfdist;
    this->_CvMDist = cvmsum / ( this->_nSamples + other._nSamples );
  }
}


template< typename RealType >
void ProbDistributionFunction1D<RealType>::printPDFforGnuplot ( ostream & gpout ) const {
  RealType valBefore = 0;
  for ( typename std::map < RealType, RealType >::const_iterator it = this->_pdf.begin(); it != this->_pdf.end(); ++it ) {
    gpout << aol::longScientificFormat ( it->first ) << " " << aol::longScientificFormat ( valBefore ) << endl
          << aol::longScientificFormat ( it->first ) << " " << aol::longScientificFormat ( it->second ) << endl;
    valBefore = it->second;
  }
  gpout << endl;
}



template< typename RealType >
ProbDistributionFunction2D<RealType>::ProbDistributionFunction2D ( ) {
  // do nothing
}


template< typename RealType >
void ProbDistributionFunction2D<RealType>::initialize ( const typename std::map< typename ProbDistributionFunction2D::Pt2d, unsigned int > &Histogram ) {
  const int sz = static_cast<int> ( Histogram.size() );

  _xyCo[0].reserve ( Histogram.size() );
  _xyCo[1].reserve ( Histogram.size() );

  qc::ScalarArray< unsigned int, qc::QC_2D > dHisto ( sz, sz );

  {
    std::set< std::pair<RealType, int> > xyInd[2];
    unsigned short hCntr = 0;
    for ( typename std::map< Pt2d, unsigned int >::const_iterator it = Histogram.begin(); it != Histogram.end(); ++it, ++hCntr ) {
      xyInd[0].insert ( std::pair<RealType, int> ( (it->first)[0], hCntr ) );
      xyInd[1].insert ( std::pair<RealType, int> ( (it->first)[1], hCntr ) );
    }

    std::map< int, int > xyRank[2];
    for ( unsigned short d = 0; d < 2; ++d ) {
      for ( typename std::set< std::pair<RealType, int> >::const_iterator it = xyInd[d].begin(); it != xyInd[d].end(); ++it ) {
        xyRank[d][ it->second ] = _xyCo[d].size();
        _xyCo[d].pushBack ( it->first );
      }
    }

    hCntr = 0;
    for ( typename std::map< Pt2d, unsigned int >::const_iterator it = Histogram.begin(); it != Histogram.end(); ++it, ++hCntr ) {
      dHisto.set ( xyRank[0][ hCntr ], xyRank[1][ hCntr ], it->second );
    }

#ifdef VERBOSE
    for ( unsigned short d = 0; d < 2; ++d ) {
      cerr << "xyInd" << endl;
      for (  typename std::set< std::pair<RealType, int> >::const_iterator it = xyInd[d].begin(); it != xyInd[d].end(); ++it ) {
        cerr << it->first << " " << it->second << endl;
      }
      cerr << endl << "_xyCo" << endl;
      for ( int i = 0; i < sz; ++i ) {
        cerr << i << " " << _xyCo[d][i] << endl;
      }
      cerr << endl << "xyRank" << endl;
      for ( typename std::map<int, int>::const_iterator it = xyRank[d].begin(); it != xyRank[d].end(); ++it ) {
        cerr << it->first << " " << it->second << endl;
      }
      cerr << endl;
    }

    cerr << "dHisto" << endl
         << dHisto << endl;
#endif
  }

  this->_nSamples = dHisto.sum();

  // sum of all {x,y} values for indices {smaller (Lower sum), larger (Upper sum)} than current index (excluded)
  qc::ScalarArray< unsigned int, qc::QC_2D > lSumX ( sz, sz ), lSumY ( sz, sz ), uSumX ( sz, sz ), uSumY ( sz, sz );

  for ( unsigned short xI = 0; xI < sz - 1; ++xI ) {
    for ( unsigned short xJ = 0; xJ < sz; ++xJ ) {
      lSumX.set ( xI + 1, xJ, lSumX.get ( xI, xJ ) + dHisto.get ( xI, xJ ) );
    }
  }
  for ( unsigned short xI = sz - 1; xI > 0; --xI ) {
    for ( unsigned short xJ = 0; xJ < sz; ++xJ ) {
      uSumX.set ( xI - 1, xJ, uSumX.get ( xI, xJ ) + dHisto.get ( xI, xJ ) );
    }
  }
  for ( unsigned short xI = 0; xI < sz; ++xI ) {
    for ( unsigned short xJ = 0; xJ < sz - 1; ++xJ ) {
      lSumY.set ( xI, xJ + 1, lSumY.get ( xI, xJ ) + dHisto.get ( xI, xJ ) );
    }
  }
  for ( unsigned short xI = 0; xI < sz; ++xI ) {
    for ( unsigned short xJ = sz - 1; xJ > 0; --xJ ) {
      uSumY.set ( xI, xJ - 1, uSumY.get ( xI, xJ ) + dHisto.get ( xI, xJ ) );
    }
  }

  // set _dPdfs to correct size
  for ( short di = 0; di < 2; ++di ) {
    for ( short dj = 0; dj < 2; ++dj ) {
      _dPdf[di][dj].reallocate ( sz, sz );
    }
  }

  // and compute their content
  // 0, 0 is <=, <=
  for ( unsigned short xI = 0; xI < sz; ++xI ) {
    for ( unsigned short xJ = 0; xJ < sz; ++xJ ) {
      // 0, 0 is <=, <=
      _dPdf[0][0].set ( xI, xJ, lSumX.get ( xI, xJ ) + lSumY.get ( xI, xJ ) + dHisto.get ( xI, xJ ) );
      if ( ( xI > 0 ) && ( xJ > 0 ) ) {
        _dPdf[0][0].add ( xI, xJ, _dPdf[0][0].get ( xI - 1, xJ - 1 ) );
      }
    }
  }
  _dPdf[0][0] /= static_cast<RealType> ( this->_nSamples );

  // 0, 1 is <=, >
  for ( unsigned short xI = 0; xI < sz; ++xI ) {
    for ( unsigned short xJ = sz - 1; xJ < sz; --xJ ) {
      _dPdf[0][1].set ( xI, xJ, lSumX.get ( xI, xJ ) + uSumY.get ( xI, xJ ) + dHisto.get ( xI, xJ ) );
      if ( ( xI > 0 ) && ( xJ < sz - 1 ) ) {
        _dPdf[0][1].add ( xI, xJ, _dPdf[0][1].get ( xI - 1, xJ + 1 ) );
      }
    }
  }
  _dPdf[0][1] /= static_cast<RealType> ( this->_nSamples );

  // 1, 0 is >, <=
  for ( unsigned short xI = sz - 1; xI < sz; --xI ) {
    for ( unsigned short xJ = 0; xJ < sz; ++xJ ) {
      _dPdf[1][0].set ( xI, xJ, uSumX.get ( xI, xJ ) + lSumY.get ( xI, xJ ) + dHisto.get ( xI, xJ ) );
      if ( ( xI < sz - 1 ) && ( xJ > 0 ) ) {
        _dPdf[1][0].add ( xI, xJ, _dPdf[1][0].get ( xI + 1, xJ - 1 ) );
      }
    }
  }
  _dPdf[1][0] /= static_cast<RealType> ( this->_nSamples );

  // 1, 1 is >, >
  for ( unsigned short xI = sz - 1; xI < sz; --xI ) {
    for ( unsigned short xJ = sz - 1; xJ < sz; --xJ ) {
      _dPdf[1][1].set ( xI, xJ, uSumX.get ( xI, xJ ) + uSumY.get ( xI, xJ ) + dHisto.get ( xI, xJ ) );
      if ( ( xI < sz - 1 ) && ( xJ < sz - 1 ) ) {
        _dPdf[1][1].add ( xI, xJ, _dPdf[1][1].get ( xI + 1, xJ + 1 ) );
      }
    }
  }
  _dPdf[1][1] /= static_cast<RealType> ( this->_nSamples );

#ifdef VERBOSE
  // dump all the stuff computed above
  cerr << "lSumX" << endl;
  for ( unsigned short xI = 0; xI < sz; ++xI ) {
    for ( unsigned short xJ = 0; xJ < sz; ++xJ ) {
      cerr << aol::shortFormat ( lSumX.get ( xI, xJ ) ) << " ";
    }
    cerr << endl;
  }
  cerr << endl;

  cerr << "lSumY" << endl;
  for ( unsigned short xI = 0; xI < sz; ++xI ) {
    for ( unsigned short xJ = 0; xJ < sz; ++xJ ) {
      cerr << aol::shortFormat ( lSumY.get ( xI, xJ ) ) << " ";
    }
    cerr << endl;
  }
  cerr << endl;

  cerr << "uSumX" << endl;
  for ( unsigned short xI = 0; xI < sz; ++xI ) {
    for ( unsigned short xJ = 0; xJ < sz; ++xJ ) {
      cerr << aol::shortFormat ( uSumX.get ( xI, xJ ) ) << " ";
    }
    cerr << endl;
  }
  cerr << endl;

  cerr << "uSumY" << endl;
  for ( unsigned short xI = 0; xI < sz; ++xI ) {
    for ( unsigned short xJ = 0; xJ < sz; ++xJ ) {
      cerr << aol::shortFormat ( uSumY.get ( xI, xJ ) ) << " ";
    }
    cerr << endl;
  }
  cerr << endl;

  for ( unsigned short di = 0; di < 2; ++di ) {
    for ( unsigned short dj = 0; dj < 2; ++dj ) {
      cerr << "dPdf " << di << " " << dj << endl
           << _dPdf[di][dj] << endl;
    }
  }
#endif
}


template< typename RealType >
void ProbDistributionFunction2D<RealType>::dump ( ostream& out ) const {
  for ( qc::RectangularIterator< qc::QC_2D, aol::Vec2<short> > rit ( _dPdf[0][0] ); rit.notAtEnd(); ++rit ) {
    const Pt2d coord = getCoord ( *rit );
    out << aol::longScientificFormat ( coord[0] ) << " "
        << aol::longScientificFormat ( coord[1] );
    for ( unsigned short di = 0; di < 2; ++di ) {
      for ( unsigned short dj = 0; dj < 2; ++dj ) {
        out << " " << aol::longScientificFormat ( _dPdf[di][dj].get ( *rit ) );
      }
    }
    out << endl;
  }
}


template< typename RealType >
RealType ProbDistributionFunction2D<RealType>::getScaledL2DistanceTo ( const ProbDistributionFunction2D &other ) const {
  const Pt2d
    minVal ( aol::Min ( this->_xyCo[0][0], other._xyCo[0][0] ), aol::Min ( this->_xyCo[1][0], other._xyCo[1][0] ) ),
    maxVal ( aol::Max ( this->_xyCo[0][ this->_xyCo[0].size() - 1 ], other._xyCo[0][ other._xyCo[0].size() - 1 ] ), aol::Max ( this->_xyCo[1][ this->_xyCo[1].size() - 1 ], other._xyCo[1][ other._xyCo[1].size() - 1 ] ) );
  const RealType factor = sqrt ( ( static_cast<RealType> ( this->_nSamples ) * static_cast<RealType> ( other._nSamples ) ) / (  this->_nSamples + other._nSamples ) );
  return ( factor / sqrt ( ( maxVal - minVal ).prod() ) * this->_L2Dist );
}


template< typename RealType >
RealType ProbDistributionFunction2D<RealType>::getScaledCvMDistanceTo ( const ProbDistributionFunction2D &/*other*/ ) const {
  return ( this->_CvMDist );
}


template< typename RealType >
void ProbDistFuncHelper<RealType>::doCompute2DPDFdistTo ( const aol::MultiVector< RealType > &thisXYCo, const aol::MultiVector< RealType > &otherXYCo,
                                                          const qc::ScalarArray< RealType, qc::QC_2D > &thisDPdf, const qc::ScalarArray< RealType, qc::QC_2D > &otherDPdf,
                                                          const unsigned int thisNSamples, const unsigned int otherNSamples,
                                                          RealType &L2Dist, RealType& LInfDist, RealType &CvMDist ) {

  std::vector<RealType> stepCoord[2];
  qc::AArray < typename ProbDistributionFunction2D<RealType>::PDFDiffStep2D, qc::QC_2D > dSteps;

  std::map < RealType, aol::Vec3<unsigned short> > stepInd[2];
  for ( unsigned short d = 0; d < 2; ++d ) {
    unsigned short cntr = 0;
    for ( int i = 0; i < thisXYCo[d].size(); ++i, ++cntr ) {
      stepInd[d][ ( thisXYCo[d] )[i] ] += aol::Vec3<unsigned short> ( cntr, 0, 1 ); // in this
    }
    cntr = 0;
    for ( int i = 0; i < otherXYCo[d].size(); ++i, ++cntr ) {
      stepInd[d][ ( otherXYCo[d] )[i] ] += aol::Vec3<unsigned short> ( 0, cntr, 2 ); // in other; last index is 3 if it occurs in both
    }
  }

#ifdef VERBOSE
  for ( unsigned short d = 0; d < 2; ++d ) {
    cerr << "dumping steps " << d << endl;
    for ( typename std::map < RealType, aol::Vec3<unsigned short> >::const_iterator it = stepInd[d].begin(); it != stepInd[d].end(); ++it ) {
      cerr << aol::detailedFormat ( it->first  ) << " " << it->second << endl;
    }
    cerr << endl;
  }
#endif

  typename std::map< unsigned short, unsigned short > stepRankThis[2], stepRankOther[2];
  for ( unsigned short d = 0; d < 2; ++d ) {
    stepCoord[d].reserve( thisXYCo[d].size() + otherXYCo[d].size() );

    for ( typename std::map< RealType, aol::Vec3<unsigned short> >::const_iterator it = stepInd[d].begin(); it != stepInd[d].end(); ++it ) {
      if ( (it->second)[2] & 1 ) { // occurs in this (but not necessarily exclusively there)
        stepRankThis [d][ (it->second)[0] ] = static_cast< unsigned short > ( stepCoord[d].size() );
      }
      if ( (it->second)[2] & 2 ) { // occurs in other
        stepRankOther[d][ (it->second)[1] ] = static_cast< unsigned short > ( stepCoord[d].size() );
      }
      stepCoord[d].push_back( it->first );
    }
  }

#ifdef VERBOSE
  for ( unsigned short d = 0; d < 2; ++d ) {
    cerr << "dumping step ranks this " << d << endl;
    for ( typename std::map< unsigned short, unsigned short >::const_iterator it = stepRankThis[d].begin(); it != stepRankThis[d].end(); ++it ) {
      cerr << it->first << " " << it->second << endl;
    }
    cerr << endl;
  }
  for ( unsigned short d = 0; d < 2; ++d ) {
    cerr << "dumping step ranks other " << d << endl;
    for ( typename std::map< unsigned short, unsigned short >::const_iterator it = stepRankOther[d].begin(); it != stepRankOther[d].end(); ++it ) {
      cerr << it->first << " " << it->second << endl;
    }
    cerr << endl;
  }
#endif

  // compute dSteps
  dSteps.reallocate ( stepCoord[0].size(), stepCoord[1].size() ); // coordinates may be duplicate for the two pdf, this is taken care of by mapping in the loops below

  for ( unsigned short xI = 0; xI < static_cast< unsigned short > ( (thisXYCo)[0].size() ); ++xI ){
    RealType valOld = 0;
    for ( unsigned short xJ = 0; xJ < static_cast< unsigned short > ( (thisXYCo)[1].size() ); ++xJ ){
      const RealType val = thisDPdf.get ( xI, xJ );
      dSteps.getRef ( stepRankThis[0][ xI ], stepRankThis[1][ xJ ] )._val[0] = val;
      dSteps.getRef ( stepRankThis[0][ xI ], stepRankThis[1][ xJ ] )._step[0][0] = val - valOld;
      valOld = val;
    }
  }

  for ( unsigned short xJ = 0; xJ < static_cast< unsigned short > ( (thisXYCo)[1].size() ); ++xJ ){
    RealType valOld = 0;
    for ( unsigned short xI = 0; xI < static_cast< unsigned short > ( (thisXYCo)[0].size() ); ++xI ){
      const RealType val = thisDPdf.get ( xI, xJ );
      // val already set above
      dSteps.getRef ( stepRankThis[0][ xI ], stepRankThis[1][ xJ ] )._step[0][1] = val - valOld;
      valOld = val;
    }
  }

  for ( unsigned short xI = 0; xI < static_cast< unsigned short > ( (otherXYCo)[0].size() ); ++xI ){
    RealType valOld = 0;
    for ( unsigned short xJ = 0; xJ < static_cast< unsigned short > ( (otherXYCo)[1].size() ); ++xJ ){
      const RealType val = otherDPdf.get ( xI, xJ );
      dSteps.getRef ( stepRankOther[0][ xI ], stepRankOther[1][ xJ ] )._val[1] = val;
      dSteps.getRef ( stepRankOther[0][ xI ], stepRankOther[1][ xJ ] )._step[1][0] = val - valOld;
      valOld = val;
    }
  }

  for ( unsigned short xJ = 0; xJ < static_cast< unsigned short > ( (otherXYCo)[1].size() ); ++xJ ){
    RealType valOld = 0;
    for ( unsigned short xI = 0; xI < static_cast< unsigned short > ( (otherXYCo)[0].size() ); ++xI ){
      const RealType val = otherDPdf.get ( xI, xJ );
      // val already set above
      dSteps.getRef ( stepRankOther[0][ xI ], stepRankOther[1][ xJ ] )._step[1][1] = val - valOld;
      valOld = val;
    }
  }

#ifdef VERBOSE
  cerr << "dumping steps and values before filling" << endl;
  for ( qc::RectangularIterator< qc::QC_2D > rit ( dSteps ); rit.notAtEnd(); ++rit ) {
    cerr << (*rit)[0] << " " << (*rit)[1] << " "
         << dSteps.getRef ( *rit )._val << " "
         << dSteps.getRef ( *rit )._step[0] << " " << dSteps.getRef ( *rit )._step[1] << endl;
  }
  cerr << endl;
#endif

  // fill dSteps
  for ( unsigned short xI = 0; xI < static_cast< unsigned short > ( stepCoord[0].size() ); ++xI ) {
    aol::Vec2<RealType> oldVal;
    for ( unsigned short xJ = 0; xJ < static_cast< unsigned short > ( stepCoord[1].size() ); ++xJ ) {
      for ( short d = 0; d < 2; ++d ) {
        if ( dSteps.getRef ( xI, xJ )._val[d] == 0 ) {
          dSteps.getRef ( xI, xJ )._val[d] = oldVal[d];
        }
        oldVal[d] = dSteps.getRef ( xI, xJ )._val[d];
      }
    }
  }

  for ( unsigned short xJ = 0; xJ < static_cast< unsigned short > ( stepCoord[1].size() ); ++xJ ) {
    aol::Vec2<RealType> oldVal;
    for ( unsigned short xI = 0; xI < static_cast< unsigned short > ( stepCoord[0].size() ); ++xI ) {
      for ( short d = 0; d < 2; ++d ) {
        if ( dSteps.getRef ( xI, xJ )._val[d] == 0 ) {
          dSteps.getRef ( xI, xJ )._val[d] = oldVal[d];
        }
        oldVal[d] = dSteps.getRef ( xI, xJ )._val[d];
      }
    }
  }

#ifdef VERBOSE
  cerr << "dumping steps and values after filling" << endl;
  for ( qc::RectangularIterator< qc::QC_2D > rit ( dSteps ); rit.notAtEnd(); ++rit ) {
    cerr << (*rit)[0] << " " << (*rit)[1] << " "
         << dSteps.getRef ( *rit )._val << " "
         << dSteps.getRef ( *rit )._step[0] << " " << dSteps.getRef ( *rit )._step[1] << endl;
  }
  cerr << endl;

  cerr << "dumping steps and values after filling, once again" << endl;
  for ( short i = 0; i < dSteps.getNumX(); ++i ) {
    for ( short j = 0; j < dSteps.getNumX(); ++j ) {
      cerr << aol::detailedFormat ( dSteps.getRef ( i, j )._val[0] - dSteps.getRef ( i, j )._val[1] ) << "  ";
    }
    cerr << endl;
  }
  cerr << endl;
#endif

  // compute differences
  RealType l2sum = 0, cvmsum = 0, linfdist = -aol::NumberTrait<RealType>::Inf;
  for ( unsigned short xI = 0; xI < static_cast< unsigned short > ( stepCoord[0].size() ); ++xI ) {
    for ( unsigned short xJ = 0; xJ < static_cast< unsigned short > ( stepCoord[1].size() ); ++xJ ) {
      const aol::Vec2<RealType> &values = dSteps.getRef( xI, xJ )._val;
      const RealType vdiff = values[1] - values[0];
      linfdist = aol::Max ( fabs ( vdiff ), linfdist );

      if ( ( xI != static_cast< unsigned short > ( stepCoord[0].size() ) - 1 ) && ( xJ != static_cast< unsigned short > ( stepCoord[1].size() ) - 1 ) ) {
        {
          aol::Vec2<RealType> xydiff ( stepCoord[0][ xI + 1 ] - stepCoord[0][ xI ], stepCoord[1][ xJ + 1 ] - stepCoord[1][ xJ ] );
          l2sum += aol::Sqr ( vdiff ) * xydiff.prod();
        }

        const aol::Matrix22<RealType> &step = dSteps.getRef( xI, xJ )._step; // the last step does not matter because the difference is zero

        cvmsum += ( step[0][0] + step[0][1] ) * aol::Sqr ( vdiff ) * thisNSamples;
        cvmsum += ( step[1][0] + step[1][1] ) * aol::Sqr ( vdiff ) * otherNSamples;
      }
    }
  }

  L2Dist = sqrt ( l2sum );
  LInfDist = linfdist;
  CvMDist = cvmsum / ( thisNSamples + otherNSamples );
}



template< typename RealType >
void ProbDistributionFunction2D<RealType>::computePDFdistTo ( const ProbDistributionFunction2D &other ) {
  if ( ( this->_xyCo[0].size() != 0 ) && ( other._xyCo[0].size() != 0 ) ) {
  qc::ScalarArray<RealType, qc::QC_2D> L2DistVals(2,2), LInfDistVals(2,2), CvMDistVals(2,2);

  aol::Vector<RealType> thisFliXYCo[2], otherFliXYCo[2];
  for ( unsigned short dd = 0; dd < 2; ++dd ) {
    thisFliXYCo[dd].resize ( this->_xyCo[dd].size() );
    thisFliXYCo[dd].revertOrderFrom ( this->_xyCo[dd] );
    thisFliXYCo[dd].revert();

    otherFliXYCo[dd].resize ( other._xyCo[dd].size() );
    otherFliXYCo[dd].revertOrderFrom ( other._xyCo[dd] );
    otherFliXYCo[dd].revert();
  }

  aol::MultiVector<RealType> thisFlippedXYCo[2][2], otherFlippedXYCo[2][2];
  thisFlippedXYCo[0][0].appendReference ( this->_xyCo[0] );    thisFlippedXYCo[0][0].appendReference ( this->_xyCo[1] );
  thisFlippedXYCo[0][1].appendReference ( this->_xyCo[0] );    thisFlippedXYCo[0][1].appendReference ( thisFliXYCo[1] );
  thisFlippedXYCo[1][0].appendReference ( thisFliXYCo[0] );    thisFlippedXYCo[1][0].appendReference ( this->_xyCo[1] );
  thisFlippedXYCo[1][1].appendReference ( thisFliXYCo[0] );    thisFlippedXYCo[1][1].appendReference ( thisFliXYCo[1] );

  otherFlippedXYCo[0][0].appendReference ( other._xyCo[0] );   otherFlippedXYCo[0][0].appendReference ( other._xyCo[1]  );
  otherFlippedXYCo[0][1].appendReference ( other._xyCo[0] );   otherFlippedXYCo[0][1].appendReference ( otherFliXYCo[1] );
  otherFlippedXYCo[1][0].appendReference ( otherFliXYCo[0] );  otherFlippedXYCo[1][0].appendReference ( other._xyCo[1]  );
  otherFlippedXYCo[1][1].appendReference ( otherFliXYCo[0] );  otherFlippedXYCo[1][1].appendReference ( otherFliXYCo[1] );


  qc::ScalarArray< RealType, qc::QC_2D > thisFlippedDPdf[2][2], otherFlippedDPdf[2][2];
  {
    const int Ni = this->_dPdf[0][0].getNumX(), Nj = this->_dPdf[0][0].getNumY(); // all have same size
    for ( short di = 0; di < 2; ++di ) {
      for ( short dj = 0; dj < 2; ++dj ) {
        thisFlippedDPdf[di][dj].reallocate ( Ni, Nj );
      }
    }

    for ( int j = 0; j < Nj; ++j ) {
      for ( int i = 0; i < Ni; ++i ) {
        thisFlippedDPdf[0][0].set ( i, j, this->_dPdf[0][0].get (          i,          j ) );
        thisFlippedDPdf[0][1].set ( i, j, this->_dPdf[0][1].get (          i, Nj - 1 - j ) );
        thisFlippedDPdf[1][0].set ( i, j, this->_dPdf[1][0].get ( Ni - 1 - i, j          ) );
        thisFlippedDPdf[1][1].set ( i, j, this->_dPdf[1][1].get ( Ni - 1 - i, Nj - 1 - j ) );
      }
    }
  }

  {
    const int Ni = other._dPdf[0][0].getNumX(), Nj = other._dPdf[0][0].getNumY(); // all have same size
    for ( short di = 0; di < 2; ++di ) {
      for ( short dj = 0; dj < 2; ++dj ) {
        otherFlippedDPdf[di][dj].reallocate ( Ni, Nj );
      }
    }

    for ( int j = 0; j < Nj; ++j ) {
      for ( int i = 0; i < Ni; ++i ) {
        otherFlippedDPdf[0][0].set ( i, j, other._dPdf[0][0].get (          i,          j ) );
        otherFlippedDPdf[0][1].set ( i, j, other._dPdf[0][1].get (          i, Nj - 1 - j ) );
        otherFlippedDPdf[1][0].set ( i, j, other._dPdf[1][0].get ( Ni - 1 - i,          j ) );
        otherFlippedDPdf[1][1].set ( i, j, other._dPdf[1][1].get ( Ni - 1 - i, Nj - 1 - j ) );
      }
    }

#ifdef VERBOSE
    for ( unsigned short di = 0; di < 2; ++di ) {
      for ( unsigned short dj = 0; dj < 2; ++dj ) {
        cerr << "thisFlippedDPdf " << di << " " << dj << endl
             << thisFlippedXYCo[di][dj] << endl
             << thisFlippedDPdf[di][dj] << endl
             << "otherFlippedDPdf " << di << " " << dj << endl
             << otherFlippedXYCo[di][dj] << endl
             << otherFlippedDPdf[di][dj] << endl;
      }
    }
#endif
  }


#ifdef _OPENMP
#pragma omp parallel for
#endif
  for ( short i = 0; i < 4; ++i ) {
    const unsigned short di = i / 2, dj = i % 2; // for simpler parallelization, both just run from 0 to 1
    RealType l2=0.0, li=0.0, cvm=0.0;
    ProbDistFuncHelper<RealType>::doCompute2DPDFdistTo ( thisFlippedXYCo[di][dj], otherFlippedXYCo[di][dj], thisFlippedDPdf[di][dj], otherFlippedDPdf[di][dj], this->_nSamples, other._nSamples, l2, li, cvm );
    L2DistVals.set(di,dj, l2);
    LInfDistVals.set(di,dj, li);
    CvMDistVals.set(di,dj, cvm);
  }

#ifdef VERBOSE
  cerr << "L2   distance values ";   L2DistVals.dump();
  cerr << "Linf distance values ";   LInfDistVals.dump();
  cerr << "CvM  distance values ";   CvMDistVals.dump();
#endif

  // is this the correct thing to do?
  this->_L2Dist   = L2DistVals.getMeanValue();
  this->_LInfDist = LInfDistVals.getMeanValue();
  this->_CvMDist  = CvMDistVals.getMeanValue();
  }
}


template class ProbDistFuncHelper<float>;
template class ProbDistributionFunctionAnyD<float>;
template class ProbDistributionFunction1D<float>;
template class ProbDistributionFunction1DForSample<float>;
template class ProbDistributionFunction1DForDiscreteHisto<float>;
template class ProbDistributionFunction1DForVecHisto<float>;
template class ProbDistributionFunction2D<float>;
template class ProbDistributionFunction2DForSample<float>;

template class ProbDistFuncHelper<double>;
template class ProbDistributionFunctionAnyD<double>;
template class ProbDistributionFunction1D<double>;
template class ProbDistributionFunction1DForSample<double>;
template class ProbDistributionFunction1DForDiscreteHisto<double>;
template class ProbDistributionFunction1DForVecHisto<double>;
template class ProbDistributionFunction2D<double>;
template class ProbDistributionFunction2DForSample<double>;

template class ProbDistFuncHelper<long double>;
template class ProbDistributionFunctionAnyD<long double>;
template class ProbDistributionFunction1D<long double>;
template class ProbDistributionFunction1DForSample<long double>;
template class ProbDistributionFunction1DForDiscreteHisto<long double>;
template class ProbDistributionFunction1DForVecHisto<long double>;
template class ProbDistributionFunction2D<long double>;
template class ProbDistributionFunction2DForSample<long double>;

}
