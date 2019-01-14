#ifndef __PATTERNANALYSIS_H
#define __PATTERNANALYSIS_H

#include <aol.h>
#include <ctrlCCatcher.h>
#include <multiArray.h>
#include <convolution.h>
#ifdef USE_MODULES_QT
#include <customPlotHandler.h>
#endif
#include <FEOpInterface.h>
#include <paramReg.h>
#include <scalarArray.h>
#include <solver.h>
#include <fastUniformGridMatrix.h>
#include <linearSmoothOp.h>
#include <configurators.h>
#include <generator.h>
#include <regression.h>
#include <cellCenteredGrid.h>
#include <linearSmoothOp.h>
#include <clustering.h>
#include <bumpFit.h>
#include <atomFinder.h>
#include <gradientDescent.h>
#include <probDistributionFunction.h>
#include <featureSegmentation.h>
#include <derivativeFreeOptimization.h>

namespace im {

template <typename RealType>
void getMeanFilteredImage ( const qc::ScalarArray<RealType, qc::QC_2D> &Data, qc::ScalarArray<RealType, qc::QC_2D> &FilteredData,
                            bool Crop = true, short FilterSize = 7 ) {
  const short filterSize = FilterSize, filterOffset = ( filterSize - 1 ) / 2;
  if ( Crop ) FilteredData.reallocate ( Data.getNumX ( ) - filterSize + 1, Data.getNumY ( ) - filterSize + 1 );
  else {
    FilteredData.reallocate ( Data.getNumX ( ), Data.getNumY ( ) );
    FilteredData = Data;
  }
  qc::ScalarArray<RealType, qc::QC_2D> block ( filterSize, filterSize );
  const int pasteOffset = Crop ? -filterOffset : 0;
  for ( int x=filterOffset; x<Data.getNumX ( )-filterOffset ; ++x ) {
    for ( int y=filterOffset; y<Data.getNumY ( )-filterOffset ; ++ y) {
      Data.copyBlockTo ( x-filterOffset, y-filterOffset, block );
      if ( aol::isFinite<RealType> ( block.sum ( ) ) ) FilteredData.set ( x + pasteOffset, y + pasteOffset, block.getMeanValue ( ) );
      else FilteredData.set ( x + pasteOffset, y + pasteOffset, Data.get ( x, y ) );
    }
  }
}
  
template <typename RealType>
void getMeanFilteredImage ( const qc::ScalarArray<RealType, qc::QC_3D> &Data, qc::ScalarArray<RealType, qc::QC_3D> &FilteredData, bool Crop = true ) {
  qc::ScalarArray<RealType, qc::QC_2D> slice ( Data.getNumX ( ), Data.getNumY ( ) ), filteredSlice ( slice );
  for ( int z=0; z<Data.getNumZ ( ) ; ++z ) {
    Data.getSlice ( qc::QC_Z, z, slice );
    getMeanFilteredImage ( slice, filteredSlice, Crop );
    FilteredData.putSlice ( qc::QC_Z, z, filteredSlice );
  }
}

template <typename RealType>
void getMeanFilteredVectorPeriodic ( const aol::Vector<RealType> &Data, aol::Vector<RealType> &FilteredData ) {
  qc::ScalarArray<RealType, qc::QC_1D> data ( Data );
  const short filterSize = 5, filterOffset = ( filterSize - 1 ) / 2;
  FilteredData.reallocate ( Data.size ( ) );
  aol::Vector<RealType> block ( filterSize );
  for ( int x=0; x<Data.size ( ) ; ++x ) {
    for ( int dx=-filterOffset; dx<=filterOffset ; ++dx ) {
      block[dx+filterOffset] = data.getPeriodic ( x + dx );
      if ( aol::isFinite<RealType> ( block.sum ( ) ) ) FilteredData.set ( x, block.getMeanValue ( ) );
      else FilteredData.set ( x, Data[x] );
    }
  }
}

template <typename RealType>
void getMeanFilteredVector ( const aol::Vector<RealType> &Data, aol::Vector<RealType> &FilteredData, bool Crop = true ) {
  const short filterSize = 3, filterOffset = ( filterSize - 1 ) / 2;
  if ( Crop ) FilteredData.reallocate ( Data.size ( ) - filterSize + 1 );
  else {
    FilteredData.reallocate ( Data.size ( ) );
    FilteredData = Data;
  }
  qc::ScalarArray<RealType, qc::QC_1D> data ( Data );
  aol::Vector<RealType> block ( filterSize );
  const int pasteOffset = Crop ? -filterOffset : 0;
  for ( int x=filterOffset; x<Data.size ( )-filterOffset ; ++x ) {
    for ( int dx=-filterOffset; dx<=filterOffset ; ++dx ) block[filterOffset+dx] = Data[x+dx];
    if ( aol::isFinite<RealType> ( block.sum ( ) ) ) FilteredData.set ( x-filterOffset, block.getMeanValue ( ) );
    else FilteredData.set ( x + pasteOffset, Data[x] );
  }
}
  
template <typename RealType, const qc::Dimension Dim, typename PictureType>
qc::CoordType getAtomCenterReferenceCoordinate ( const PictureType &Data,
                                                 const qc::CoordType &X0 = qc::CoordType ( 0, 0 ), const qc::CoordType &XEnd = qc::CoordType ( -1, -1 ),
                                                 const bool ForceToImgCenter = true, const bool MeanFilterImg = true ) {
  qc::CoordType xEnd ( XEnd );
  PictureType data ( Data );
  qc::GridSize<Dim> size ( data.getSize ( ) );
  if ( XEnd[0] < 0 ) for ( int d=0; d<Dim ; ++d ) xEnd[d] = size[d];
  if ( ForceToImgCenter ) {
    data.setAll ( -aol::NumberTrait<RealType>::Inf );
    const int padding = 21;
    aol::Vec<Dim, RealType> center;
    int n = 0;
    for ( qc::RectangularIterator<Dim> it ( data ); it.notAtEnd ( ) ; ++it ) {
      if ( aol::isFinite<RealType> ( Data.get ( *it ) ) ) {
        for ( int d=0; d<Dim ; ++d ) center[d] += (*it)[d];
        ++n;
      }
    }
    center /= static_cast<RealType> ( n );
    qc::CoordType lower, upper;
    for ( int d=0; d<Dim ; ++d ) {
      lower[d] = aol::Clamp<short> ( round ( center[d] ) - padding, 0, size[d] );
      upper[d] = aol::Clamp<short> ( round ( center[d] ) + padding, 0, size[d] );
    }
    for ( qc::RectangularIterator<Dim> it ( lower, upper ); it.notAtEnd ( ) ; ++it )
      if ( aol::isFinite<RealType> ( Data.get ( *it ) ) ) data.set ( *it, Data.get ( *it ) );
  }
  PictureType meanFilteredData ( data );
  if ( MeanFilterImg ) getMeanFilteredImage<RealType> ( data, meanFilteredData, false );
  qc::CoordType xRef;
  for ( int d=0; d<Dim ; ++d ) xRef[d] = -1;
  const qc::FastILexMapper<Dim> mapper ( size );
  std::pair<int, RealType> maxIndVal ( -1, 0 );
  while ( aol::isFinite<RealType> ( maxIndVal.second ) && !aol::InsideQuad<Dim> ( X0, xEnd, xRef ) ) {
    maxIndVal = meanFilteredData.getFiniteMaxIndexAndValue ( );
    for ( int d=0; d<Dim ; ++d ) xRef[d] = mapper.splitGlobalIndex ( maxIndVal.first )[d];
    meanFilteredData.set ( xRef, -aol::NumberTrait<RealType>::Inf );
  }
  return xRef;
}
  
template <typename RealType, typename PictureType>
void getPNGYellowNaNs ( const PictureType &Arg, qc::MultiArray<RealType, qc::QC_2D, 3> &Dest ) {
  PictureType arg ( Arg );
  arg.addToAll ( -arg.getFiniteMinValue ( ) );
  arg /= arg.getFiniteMaxValue ( );
  Dest.reallocate ( arg.getNumX ( ), arg.getNumY ( ) );
  for ( int i=0; i<3 ; ++i ) Dest[i] = arg;
  for ( int k=0; k<Arg.size ( ) ; ++k ) {
    if ( !aol::isFinite<RealType> ( Arg[k] ) ) {
      Dest[0][k] = 1.0;
      Dest[1][k] = 1.0;
      Dest[2][k] = 0.0;
    }
  }
  Dest.setOverflowHandlingToCurrentValueRange ( );
}
  

template <typename RealType, typename PictureType>
void setInteriorPixels ( aol::RandomAccessContainer<aol::Vec2<int> > &InteriorPixels, const PictureType &Data, const int Margin ) {
  InteriorPixels.clear ( );
  
  PictureType nanGrid ( Data, aol::STRUCT_COPY );
  for ( int y=0; y<Data.getNumY ( ) ; ++y ) {
    for ( int x=0; x<Data.getNumX ( ) ; ++x ) {
      if ( aol::isNaN<RealType> ( Data.get ( x, y ) ) ) {
        for ( int dx=-Margin; dx<=Margin ; ++dx ) {
          for ( int dy=-Margin; dy<=Margin ; ++dy ) {
            if ( x + dx >= 0 && x + dx < Data.getNumX ( ) && y + dy >= 0 && y + dy < Data.getNumY ( ) )
              nanGrid.set ( x+dx, y+dy, aol::NumberTrait<RealType>::NaN );
          }
        }
      }
    }
  }
  for ( int y=Margin; y<Data.getNumY ( )-Margin ; ++y ) {
    for ( int x=Margin; x<Data.getNumX ( )-Margin ; ++x ) {
      if ( !aol::isNaN<RealType> ( nanGrid.get ( x, y ) ) )
        InteriorPixels.pushBack ( aol::Vec2<int> ( x, y ) );
    }
  }
}

template <typename _RealType>
class SumOfSinesLeastSquaresEnergy : public aol::Op<aol::Vector<_RealType>, aol::Scalar<_RealType> > {
  typedef _RealType RealType;
protected:
  const std::vector<std::pair<RealType, RealType> > &_data;
  const int _numTerms;
public:
  SumOfSinesLeastSquaresEnergy ( const std::vector<std::pair<RealType, RealType> > &Data, const int NumTerms )
  : _data ( Data ), _numTerms ( NumTerms ) { }
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    if ( Arg.size ( ) != 3 * _numTerms )
      throw aol::Exception ( "Arguments size does not match 3 * number of sin terms!", __FILE__, __LINE__ );
    
    for ( unsigned int i=0; i<_data.size ( ) ; ++i )
      Dest += aol::Sqr<RealType> ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
    Dest *= 0.5;
  }
  
  RealType sumOfSines ( const RealType X, const aol::Vector<RealType> &Arg ) const {
    RealType res = 0;
    for ( short k=0; k<_numTerms ; ++k )
      res += Arg[3*k] * std::sin ( Arg[3*k+1] * X + Arg[3*k+2] );
    
    return res;
  }
  
  static RealType sumOfSines ( const RealType X, const aol::Vector<RealType> &Arg, const short NumTerms ) {
    RealType res = 0;
    for ( short k=0; k<NumTerms ; ++k )
      res += Arg[3*k] * std::sin ( Arg[3*k+1] * X + Arg[3*k+2] );
    
    return res;
  }
};

template <typename _RealType>
class SumOfSinesLeastSquaresEnergyDerivative : public aol::Op<aol::Vector<_RealType> > {
  typedef _RealType RealType;
protected:
  const std::vector<std::pair<RealType, RealType> > &_data;
  const int _numTerms;
public:
  SumOfSinesLeastSquaresEnergyDerivative ( const std::vector<std::pair<RealType, RealType> > &Data, const int NumTerms )
  : _data ( Data ), _numTerms ( NumTerms ) { }
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Arg.size ( ) != 3 * _numTerms )
      throw aol::Exception ( "Arguments size does not match 3 * number of sin terms!", __FILE__, __LINE__ );
    
    if ( Dest.size ( ) != Arg.size ( ) )
      throw aol::Exception ( "Destination vector (gradient) size does not match argument size!", __FILE__, __LINE__ );
    
    for ( unsigned int i=0; i<_data.size ( ) ; ++i ) {
      const RealType x = _data[i].first;
      for ( short k=0; k<_numTerms ; ++k ) {
        Dest[3*k] += ( sumOfSines ( _data[i].first, Arg ) - _data[i].second ) * std::sin ( Arg[3*k+1] * x + Arg[3*k+2] );
        Dest[3*k+1] += ( sumOfSines ( _data[i].first, Arg ) - _data[i].second ) * Arg[3*k] * x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] );
        Dest[3*k+2] += ( sumOfSines ( _data[i].first, Arg ) - _data[i].second ) * Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] );
      }
    }
  }
  
  RealType sumOfSines ( const RealType X, const aol::Vector<RealType> &Arg ) const {
    RealType res = 0;
    for ( short k=0; k<_numTerms ; ++k )
      res += Arg[3*k] * std::sin ( Arg[3*k+1] * X + Arg[3*k+2] );
    
    return res;
  }
};

template <typename _RealType>
class SumOfSinesLeastSquaresEnergySecondDerivative : public aol::Op<aol::Vector<_RealType>, aol::FullMatrix<_RealType> > {
  typedef _RealType RealType;
protected:
  const std::vector<std::pair<RealType, RealType> > &_data;
  const int _numTerms;
public:
  SumOfSinesLeastSquaresEnergySecondDerivative ( const std::vector<std::pair<RealType, RealType> > &Data, const int NumTerms )
  : _data ( Data ), _numTerms ( NumTerms ) { }
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::FullMatrix<RealType> &Dest ) const {
    if ( Arg.size ( ) != 3 * _numTerms )
      throw aol::Exception ( "Arguments size does not match 3 * number of sin terms!", __FILE__, __LINE__ );
    
    if ( Dest.getNumRows ( ) != Arg.size ( ) || Dest.getNumCols ( ) != Arg.size ( ) )
      throw aol::Exception ( "Destination matrix (Hessian) row and/or column size does not match arguments size!", __FILE__, __LINE__ );
    
    for ( unsigned int i=0; i<_data.size ( ) ; ++i ) {
      const RealType x = _data[i].first;
      for ( short k=0; k<_numTerms ; ++k ) {
        for ( short l=0; l<_numTerms ; ++l ) {
          if ( k == l ) {
            Dest.ref ( 3*k, 3*k ) += aol::Sqr<RealType> ( std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) );
            Dest.ref ( 3*k, 3*k+1 ) += std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) * Arg[3*k] * x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] )
            + x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
            Dest.ref ( 3*k, 3*k+2 ) += std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) * Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] )
            + std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
            Dest.ref ( 3*k+1, 3*k ) += Arg[3*k] * x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * std::sin ( Arg[3*k+1] * x + Arg[3*k+2] )
            + x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
            Dest.ref ( 3*k+1, 3*k+1 ) += aol::Sqr<RealType> ( Arg[3*k] * x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) )
            - Arg[3*k] * aol::Sqr<RealType> ( x ) * std::sin ( Arg[3*k+1] * x + Arg[3*k+2] )
            * ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
            Dest.ref ( 3*k+1, 3*k+2 ) += x * aol::Sqr<RealType> ( Arg[3*k] *  std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) )
            - Arg[3*k] * x * std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) * ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
            Dest.ref ( 3*k+2, 3*k ) += Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * std::sin ( Arg[3*k+1] * x + Arg[3*k+2] )
            + std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
            Dest.ref ( 3*k+2, 3*k+1 ) += x * aol::Sqr<RealType> ( Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) )
            - Arg[3*k] * x * std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) * ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
            Dest.ref ( 3*k+2, 3*k+2 ) += aol::Sqr<RealType> ( Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) )
            - Arg[3*k] * std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) * ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
          } else {
            Dest.ref ( 3*k, 3*l ) += std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) * std::sin ( Arg[3*l+1] * x + Arg[3*l+2] );
            Dest.ref ( 3*k, 3*l+1 ) += std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) * Arg[3*l] * x * std::cos ( Arg[3*l+1] * x + Arg[3*l+2] );
            Dest.ref ( 3*k, 3*l+2 ) += std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) * Arg[3*l] * std::cos ( Arg[3*l+1] * x + Arg[3*l+2] );
            Dest.ref ( 3*k+1, 3*l ) += Arg[3*k] * x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * std::sin ( Arg[3*l+1] * x + Arg[3*l+2] );
            Dest.ref ( 3*k+1, 3*l+1 ) += Arg[3*k] * x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * Arg[3*l] * x * std::cos ( Arg[3*l+1] * x + Arg[3*l+2] );
            Dest.ref ( 3*k+1, 3*l+2 ) += Arg[3*k] * x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * Arg[3*l] * std::cos ( Arg[3*l+1] * x + Arg[3*l+2] );
            Dest.ref ( 3*k+2, 3*l ) += Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * std::sin ( Arg[3*l+1] * x + Arg[3*l+2] );
            Dest.ref ( 3*k+2, 3*l+1 ) += Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * Arg[3*l] * x * std::cos ( Arg[3*l+1] * x + Arg[3*l+2] );
            Dest.ref ( 3*k+2, 3*l+2 ) += Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * Arg[3*l] * std::cos ( Arg[3*l+1] * x + Arg[3*l+2] );
          }
        }
      }
    }
  }
  
  RealType sumOfSines ( const RealType X, const aol::Vector<RealType> &Arg ) const {
    RealType res = 0;
    for ( short k=0; k<_numTerms ; ++k )
      res += Arg[3*k] * std::sin ( Arg[3*k+1] * X + Arg[3*k+2] );
    
    return res;
  }
};


template <typename _RealType>
class SumOfSinesTargetFunctional : public aol::Op<aol::Vector<_RealType> > {
  typedef _RealType RealType;
protected:
  const std::vector<std::pair<RealType, RealType> > &_data;
  const int _numTerms;
public:
  SumOfSinesTargetFunctional ( const std::vector<std::pair<RealType, RealType> > &Data, const int NumTerms )
  : _data ( Data ), _numTerms ( NumTerms ) { }
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Arg.size ( ) != 3 * _numTerms )
      throw aol::Exception ( "Arguments size does not match 3 * number of sin terms!", __FILE__, __LINE__ );
    
    if ( Dest.size ( ) != static_cast<int> ( _data.size ( ) ) )
      throw aol::Exception ( "Destination vector does not match number of data points!", __FILE__, __LINE__ );
    
    for ( unsigned int i=0; i<_data.size ( ) ; ++i )
      Dest[i] += sumOfSines ( _data[i].first, Arg ) - _data[i].second;
  }
  
  RealType sumOfSines ( const RealType X, const aol::Vector<RealType> &Arg ) const {
    RealType res = 0;
    for ( short k=0; k<_numTerms ; ++k )
      res += Arg[3*k] * std::sin ( Arg[3*k+1] * X + Arg[3*k+2] );
    
    return res;
  }
};

template <typename _RealType, typename _MatrixType>
class SumOfSinesTargetJacobian : public aol::Op<aol::Vector<_RealType>, _MatrixType> {
  typedef _RealType RealType;
protected:
  const std::vector<std::pair<RealType, RealType> > &_data;
  const int _numTerms;
public:
  SumOfSinesTargetJacobian ( const std::vector<std::pair<RealType, RealType> > &Data, const int NumTerms )
  : _data ( Data ), _numTerms ( NumTerms ) { }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::FullMatrix<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::FullMatrix<RealType> &Dest ) const {
    if ( Arg.size ( ) != 3 * _numTerms )
      throw aol::Exception ( "Arguments size does not match 3 * number of sin terms!", __FILE__, __LINE__ );
    
    if ( static_cast<unsigned int> ( Dest.getNumRows ( ) ) != _data.size ( ) || Dest.getNumCols ( ) != 3 * _numTerms )
      throw aol::Exception ( "Destination dimensions do not fit number of data points and parameters!", __FILE__, __LINE__ );
    
    for ( unsigned int i=0; i<_data.size ( ) ; ++i ) {
      const RealType x = _data[i].first;
      for ( short k=0; k<_numTerms ; ++k ) {
        Dest.set ( i, 3*k, std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) );
        Dest.set ( i, 3*k+1, Arg[3*k] * x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) );
        Dest.set ( i, 3*k+2, Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) );
      }
    }
  }
};


template <typename _RealType>
class Gaussian1DTargetFunctional : public aol::Op<aol::Vector<_RealType> > {
  typedef _RealType RealType;
protected:
  const std::vector<std::pair<RealType, RealType> > &_data;
public:
  Gaussian1DTargetFunctional ( const std::vector<std::pair<RealType, RealType> > &Data ) : _data ( Data ) { }
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Arg.size ( ) != 3 )
      throw aol::Exception ( "Arguments size is inequal to 3 (mean, variance, intensity)!", __FILE__, __LINE__ );
    
    if ( static_cast<unsigned int> ( Dest.size ( ) ) != _data.size ( ) )
      throw aol::Exception ( "Destination vector does not match number of data points!", __FILE__, __LINE__ );
    
    for ( unsigned int i=0; i<_data.size ( ) ; ++i )
      Dest[i] += Arg[2] * aol::NormalDistribution<RealType>::PDF ( _data[i].first, Arg[0], Arg[1] ) - _data[i].second;
  }
};

template <typename _RealType, typename _MatrixType>
class Gaussian1DTargetJacobian : public aol::Op<aol::Vector<_RealType>, _MatrixType> {
  typedef _RealType RealType;
protected:
  const std::vector<std::pair<RealType, RealType> > &_data;
public:
  Gaussian1DTargetJacobian ( const std::vector<std::pair<RealType, RealType> > &Data ) : _data ( Data ) { }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::FullMatrix<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::FullMatrix<RealType> &Dest ) const {
    if ( Arg.size ( ) != 3 )
      throw aol::Exception ( "Arguments size is inequal to 3 (mean, variance, intensity)!", __FILE__, __LINE__ );
    
    if ( static_cast<unsigned int> ( Dest.getNumRows ( ) ) != _data.size ( ) || Dest.getNumCols ( ) != 3 )
      throw aol::Exception ( "Destination dimensions do not fit number of data points and parameters!", __FILE__, __LINE__ );
    
    for ( unsigned int i=0; i<_data.size ( ) ; ++i ) {
      const RealType x = _data[i].first;
      Dest.set ( i, 0, Arg[2] * ( ( x - Arg[0] ) / Arg[1] * aol::NormalDistribution<RealType>::PDF ( _data[i].first, Arg[0], Arg[1] ) ) );
      Dest.set ( i, 1, Arg[2] * ( 0.5 * ( aol::Sqr<RealType> ( x - Arg[0] ) * aol::NormalDistribution<RealType>::PDF ( _data[i].first, Arg[0], Arg[1] ) / aol::Sqr<RealType> ( Arg[1] )
                                         - aol::NormalDistribution<RealType>::PDF ( _data[i].first, Arg[0], Arg[1] ) / Arg[1] ) ) );
      Dest.set ( i, 2, aol::NormalDistribution<RealType>::PDF ( _data[i].first, Arg[0], Arg[1] ) );
    }
  }
};

  
  
template <typename RealType, typename PictureType,
          typename ConfiguratorType = qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType, qc::QC_2D, 3> > >
class PeriodicityTargetFunctionalFE : public aol::FELeastSquaresFunctionalInterface<ConfiguratorType, PeriodicityTargetFunctionalFE<RealType, PictureType, ConfiguratorType>, aol::FullMatrix<RealType> > {
  typedef qc::ParametricTranslation<ConfiguratorType> ParamTransfType;
protected:
  const qc::BitArray<qc::QC_2D> &_mask;
  aol::Vector<RealType> _latticeVector;
public:
  PeriodicityTargetFunctionalFE ( const qc::BitArray<qc::QC_2D> &Mask, const aol::Vector<RealType> &LatticeVector )
    : aol::FELeastSquaresFunctionalInterface<ConfiguratorType, PeriodicityTargetFunctionalFE<RealType, PictureType, ConfiguratorType>,
                                             aol::FullMatrix<RealType> > ( typename ConfiguratorType::InitType ( qc::GridSize<qc::QC_2D> ( Mask ) ) ),
      _mask ( Mask ), _latticeVector ( LatticeVector ) {
    _latticeVector /= Mask.getNumX ( );
  }
  
  RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrFunct,
                               const typename ConfiguratorType::ElementType &El, int /*QuadPoint*/, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    if ( _mask.elementTrue ( El ) ) {
      qc::Element TransformedEl;
      typename ConfiguratorType::DomVecType TransformedLocalCoord;
      const ParamTransfType parDef ( qc::GridDefinition ( qc::GridSize<qc::QC_2D> ( _mask ) ), _latticeVector );
      parDef.template evaluateDeformation<false> ( El, RefCoord, TransformedEl, TransformedLocalCoord );
      
      return DiscrFunct.evaluate ( TransformedEl, TransformedLocalCoord ) - DiscrFunct.evaluate ( El, RefCoord );
    } else return 0;
  }
};

template <typename RealType, typename PictureType,
          typename ConfiguratorType = qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType, qc::QC_2D, 3> > >
class PeriodicityTargetJacobianFE : public aol::FELeastSquaresFunctionalInterface<ConfiguratorType, PeriodicityTargetJacobianFE<RealType, PictureType, ConfiguratorType>, aol::FullMatrix<RealType> > {
  typedef qc::ParametricTranslation<ConfiguratorType> ParamTransfType;
protected:
  const qc::BitArray<qc::QC_2D> &_mask;
  aol::Vector<RealType> _latticeVector;
  const int _gradComp;
public:
  PeriodicityTargetJacobianFE ( const qc::BitArray<qc::QC_2D> &Mask, const aol::Vector<RealType> &LatticeVector, const int GradientComponent )
    : aol::FELeastSquaresFunctionalInterface<ConfiguratorType, PeriodicityTargetJacobianFE<RealType, PictureType, ConfiguratorType>,
                                             aol::FullMatrix<RealType> > ( typename ConfiguratorType::InitType ( qc::GridSize<qc::QC_2D> ( Mask ) ) ),
      _mask ( Mask ), _latticeVector ( LatticeVector ), _gradComp ( GradientComponent ) {
    _latticeVector /= Mask.getNumX ( );
  }
  
  RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscrFunct,
                               const typename ConfiguratorType::ElementType &El, int /*QuadPoint*/, const typename ConfiguratorType::DomVecType &RefCoord ) const {
    if ( _mask.elementTrue ( El ) ) {
      qc::Element TransformedEl;
      typename ConfiguratorType::DomVecType TransformedLocalCoord;
      const ParamTransfType parDef ( qc::GridDefinition ( qc::GridSize<qc::QC_2D> ( _mask ) ), _latticeVector );
      parDef.template evaluateDeformation<true> ( El, RefCoord, TransformedEl, TransformedLocalCoord );
      
      aol::Vec<2, RealType> grad;
      DiscrFunct.evaluateGradient ( TransformedEl, TransformedLocalCoord, grad );
      aol::Mat<ParamTransfType::NumberOfDeformParameters, qc::QC_2D, RealType> jacobian;
      parDef.evaluateDerivativeDeformation ( TransformedEl, TransformedLocalCoord, jacobian );
      
      return grad[0] * jacobian[_gradComp][0] + grad[1] * jacobian[_gradComp][1];
    } else return 0;
  }
};
  
template <typename _RealType, typename _PictureType>
class PeriodicityTargetFunctional : public aol::Op<aol::Vector<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const PictureType &_data;
  qc::BitArray<qc::QC_2D> _mask;
public:
  PeriodicityTargetFunctional ( const PictureType &Data, const int Margin ) : _data ( Data ), _mask ( qc::GridSize<qc::QC_2D> ( Data ) ) {
    aol::RandomAccessContainer<aol::Vec2<int> > interiorPixels;
    setInteriorPixels<RealType, PictureType> ( interiorPixels, Data, Margin );
    for ( int i=0; i<interiorPixels.size ( ) ; ++i ) _mask.set ( interiorPixels[i], true );
  }
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Arg.size ( ) != 2 )
      throw aol::Exception ( "Arguments size must be equal to 2!", __FILE__, __LINE__ );
    
    if ( Dest.size ( ) != size ( ) )
      throw aol::Exception ( "Destination vector size does not match #pixels * #Gauss quadrature points!", __FILE__, __LINE__ );
    
    aol::Vector<RealType> dest ( size ( ) );
    PeriodicityTargetFunctionalFE<RealType, PictureType> f ( _mask, Arg );
    f.apply ( _data, dest );
    Dest += dest;
  }
  
  int size ( ) const {
    return _data.size ( ) * 4; // #pixels * #Gauss quadrature points
  }
};

template <typename _RealType, typename _MatrixType, typename _PictureType>
class PeriodicityTargetJacobian : public aol::Op<aol::Vector<_RealType>, _MatrixType> {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const PictureType &_data;
  qc::BitArray<qc::QC_2D> _mask;
public:
  PeriodicityTargetJacobian ( const PictureType &Data, const int Margin ) : _data ( Data ), _mask ( qc::GridSize<qc::QC_2D> ( Data ) ) {
    aol::RandomAccessContainer<aol::Vec2<int> > interiorPixels;
    setInteriorPixels<RealType, PictureType> ( interiorPixels, Data, Margin );
    for ( int i=0; i<interiorPixels.size ( ) ; ++i ) _mask.set ( interiorPixels[i], true );
  }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::FullMatrix<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::FullMatrix<RealType> &Dest ) const {
    if ( Arg.size ( ) != 2 )
      throw aol::Exception ( "Arguments size must be equal to 2!", __FILE__, __LINE__ );
    
    if ( Dest.getNumRows ( ) != size ( ) || Dest.getNumCols ( ) != 2 )
      throw aol::Exception ( "Destination dimensions do not match #pixels * #Gauss quadrature points and parameters!", __FILE__, __LINE__ );
    
    for ( int j=0; j<2 ; ++j ) {
      aol::Vector<RealType> dest ( size ( ) );
      PeriodicityTargetJacobianFE<RealType, PictureType> f ( _mask, Arg, j );
      f.apply ( _data, dest );
      Dest.setSubColumn ( 0, j, dest );
    }
  }
  
  int size ( ) const {
    return _data.size ( ) * 4; // #pixels * #Gauss quadrature points
  }
};
  
  
template <typename RealType, typename PictureType>
class PeriodicityEnergy : public aol::Op<aol::Vector<RealType>, aol::Scalar<RealType> > {
protected:
  const PictureType &_data;
  const int _radius;
  const aol::Vector<RealType> &_initialGuess;
  aol::RandomAccessContainer<aol::Vec2<int> > _interiorPixels;
public:
  PeriodicityEnergy ( const PictureType &Data, const int Margin, const aol::Vector<RealType> &InitialGuess )
    : _data ( Data ), _radius ( aol::Min<RealType> ( 1, aol::Abs<RealType> ( Margin - InitialGuess.norm ( ) ) ) ), _initialGuess ( InitialGuess ) {
    setInteriorPixels<RealType, PictureType> ( _interiorPixels, Data, Margin );
  }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    Dest = 0;
    aol::Vector<RealType> diff ( _initialGuess );
    diff -= Arg;
    if ( diff.norm ( ) < _radius ) { // ensure that solution stays within radius of initial guess (avoids shifted position running out of frame and convergence towards trivial minimum)
      for ( int i=0; i<_interiorPixels.size ( ) ; ++i ) {
        aol::Vec2<RealType> xpv ( _interiorPixels[i][0] + Arg[0], _interiorPixels[i][1] + Arg[1] );
        Dest += aol::Sqr<RealType> ( _data.get ( _interiorPixels[i] ) - _data.interpolate ( xpv ) );
      }
    } else
      Dest = aol::NumberTrait<RealType>::Inf;
  }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Scalar<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException ( "Not implemented", __FILE__, __LINE__ );
  }
};
  

template <typename RealType, typename PictureType>
RealType periodicityEnergy1D ( const PictureType &Data, const RealType &Angle, const RealType Period ) {
  RealType energy = 0;
  int numPoints = 0;
  aol::Vec2<RealType> shiftedPos;
  const RealType periodCosAngle = Period * cos ( Angle );
  const RealType periodSinAngle = Period * sin ( Angle );
  const RealType numX = Data.getNumX ( );
  const RealType numY = Data.getNumY ( );
  for ( int y=0; y<numY; ++y ) {
    for ( int x=0; x<numX; ++x ) {
      shiftedPos.set ( x + periodCosAngle, y + periodSinAngle );
      if ( shiftedPos[0] >= 0 && shiftedPos[0] < numX && shiftedPos[1] >= 0 && shiftedPos[1] < numY && aol::isFinite<RealType> ( Data.get ( x, y ) ) ) {
        const RealType value = Data.interpolate ( shiftedPos );
        if ( aol::isFinite<RealType> ( value ) ) {
          energy += aol::Sqr<RealType> ( Data.get ( x, y ) - value );
          ++numPoints;
        }
      }
    }
  }
  return energy / static_cast<RealType> ( numPoints );
}


enum PatternAnalysisLatticeAngleOptimizationType {
  FourierPeaks,
  ProjectiveStdDevPeaks
};
  
enum PatternAnalysisFundamentalPeriodOptimizationType {
  SineFit,
  PeriodicityEnergyMinimization
};

  
/**
 * Unsupervised recognition of primitive unit cell dimensions from atomic-scale electron micrographs.
 *
 * \author Mevenkamp
 * \ingroup ElectronMicroscopy
 */
template <typename _RealType, typename _PictureType>
class PatternAnalyzer {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
  typedef qc::MultiArray<RealType, qc::QC_2D, 3> ColoredPictureType;
  typedef qc::RectangularGridConfigurator<RealType, qc::QC_1D, aol::GaussQuadrature<RealType,qc::QC_1D,3> > ConfType;
  typedef qc::QuocConfiguratorTraitMultiLin<RealType, qc::QC_2D, aol::GaussQuadrature<RealType,qc::QC_2D,3> > GDConfType;
protected:
  const std::string &_outputDir;
  const bool _quietMode;
  aol::ProgressBar<> *_progressBar;
  int _maxNumFFTPeaks, _maxNumPSDPeaks;
  mutable sigfunc _previousCtrlCHandler;
public:
  PatternAnalyzer ( const std::string &OutputDir = "", const bool QuietMode = true,
                    aol::ProgressBar<> *ProgressBar = NULL )
    : _outputDir ( OutputDir ), _quietMode ( QuietMode ), _progressBar ( ProgressBar ), _maxNumFFTPeaks ( 4 ), _maxNumPSDPeaks ( 6 ) { }
  
  void getLatticeVectors ( aol::MultiVector<RealType> &LatticeVectors, const PictureType &Data,
                           const PatternAnalysisLatticeAngleOptimizationType &AngleOptimizationType = FourierPeaks,
                           const PatternAnalysisFundamentalPeriodOptimizationType &FundamentalPeriodOptimizationType = PeriodicityEnergyMinimization,
                           const bool RefineByPeriodicityFunctionalMinimization = false,
                           const int MeanFilterWindowSize = 7 ) {
    // Preprocess data (periodicity should not be afflicted, but analysis becomes more robust)
    PictureType preprocessedData;
    getMeanFilteredImage<RealType> ( Data, preprocessedData, true, MeanFilterWindowSize );
    if ( _outputDir.size ( ) > 0 ) {
      preprocessedData.setOverflowHandlingToCurrentValueRange ( );
      preprocessedData.save ( aol::strprintf ( "%s/preprocessedData%s", _outputDir.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    }
    
    LatticeVectors.resize ( 2, 2 );
    aol::Vector<RealType> latticeAngles, fundamentalPeriods;
    if ( AngleOptimizationType == FourierPeaks ) getLatticeAnglesFromFourierPeaks ( latticeAngles, preprocessedData, Data );
    else getLatticeAnglesFromProjectiveStdDevPeaks ( latticeAngles, preprocessedData );
    if ( FundamentalPeriodOptimizationType == SineFit ) getFundamentalPeriodsFromSineFit ( fundamentalPeriods, latticeAngles, preprocessedData );
    else getFundamentalPeriodsFromPeriodicityEnergyMinimization ( fundamentalPeriods, latticeAngles, preprocessedData );
    getPrimitiveUnitCellVectors ( LatticeVectors, fundamentalPeriods, latticeAngles, Data );
    if ( RefineByPeriodicityFunctionalMinimization ) refinePrimitiveUnitCellVectorsByPeriodicityFunctionalMinimization ( LatticeVectors, preprocessedData, Data );
  }
  
  void getLatticeAnglesAndPeriods ( aol::Vec2<RealType> &LatticeAngles, aol::Vec2<RealType> &FundamentalPeriods, const PictureType &Data,
                                    const PatternAnalysisLatticeAngleOptimizationType &AngleOptimizationType = FourierPeaks,
                                    const PatternAnalysisFundamentalPeriodOptimizationType &FundamentalPeriodOptimizationType = PeriodicityEnergyMinimization,
                                    const bool RefineByPeriodicityFunctionalMinimization = false ) {
    aol::MultiVector<RealType> latticeVectors;
    getLatticeVectors ( latticeVectors, Data, AngleOptimizationType, FundamentalPeriodOptimizationType, RefineByPeriodicityFunctionalMinimization );
    for ( int i=0; i<2 ; ++i ) {
      LatticeAngles[i] = atan ( latticeVectors[i][1] / latticeVectors[i][0] );
      FundamentalPeriods[i] = latticeVectors[i].norm ( );
    }
  }
  
  void refinePrimitiveUnitCellVectorsByPeriodicityFunctionalMinimization ( aol::MultiVector<RealType> &LatticeVectors, const PictureType &PreprocessedData, const PictureType &Data ) {
    // Functional requires square domain and cannot cope with NaN entries, thus crop to largest inner square containing only finite values
    PictureType data ( PreprocessedData );
    RealType angle = 0;
    rotateAndResampleFromLargestRectangle<RealType> ( data, angle );

    for ( int i=0; i<2 ; ++i ) {
      int margin = ceil ( LatticeVectors.norm ( ) ) + 3;

      PeriodicityEnergy<RealType, PictureType> E ( data, margin, LatticeVectors[i] );
      aol::NelderMeadDownhillSimplexAlgorithm<RealType> nelderMeadAlg ( E, 1e-16, 100, false );
      nelderMeadAlg.solve ( LatticeVectors[i] );
    }
    
    if ( !_quietMode ) std::cerr << "Refined primitive unit cell vectors:" << std::endl << LatticeVectors << std::endl;
  
    if ( _outputDir.size ( ) > 0 ) {
      aol::Vector<RealType> latticeAngles ( 2 );
      for ( int i=0; i<2 ; ++i ) latticeAngles[i] = atan2 ( LatticeVectors[i][1], LatticeVectors[i][0] );
      saveDataPeriodicityAxesImg ( aol::strprintf ( "%s/refined_latticeSymmetryAxes.png", _outputDir.c_str ( ) ).c_str ( ), Data, latticeAngles );
      saveDataPeriodicPatternImg ( aol::strprintf ( "%s/refined_crystalLattice.png", _outputDir.c_str ( ) ).c_str ( ), Data, LatticeVectors );
      saveDataLatticeVectorsFigure ( aol::strprintf ( "%s/refined_latticeVectors", _outputDir.c_str ( ) ).c_str ( ), Data, LatticeVectors );
      LatticeVectors.saveToFile ( aol::strprintf ( "%s/refined_latticeVectors.dat", _outputDir.c_str ( ) ).c_str ( ) );
    }
  }
  
  void setMaxNumFFTPeaks ( const int MaxNumFFTPeaks ) {
    _maxNumFFTPeaks = MaxNumFFTPeaks;
  }
  
  static void saveDataPeriodicPatternImg ( const char* Path, const PictureType &Data, const aol::MultiVector<RealType> &LatticeVectors,
                                           const aol::Vec2<short> &Origin = aol::Vec2<short> ( -1, -1 ), const short SearchWindowSize = 1 ) {
    qc::CoordType origin ( Origin[0], Origin[1] );
    if ( !aol::InsideQuad<qc::QC_2D> ( qc::CoordType ( 0, 0 ), qc::CoordType ( Data.getNumX ( ), Data.getNumY ( ) ), origin ) )
      origin = getAtomCenterReferenceCoordinate<RealType, qc::QC_2D> ( Data, qc::CoordType ( 0, 0 ), qc::CoordType ( Data.getNumX ( ), Data.getNumY ( ) ) );
    
    const short searchWindowOffset = ( SearchWindowSize - 1 ) / 2;
    
    ColoredPictureType periodicityAxesImg;
    getPNGYellowNaNs ( Data, periodicityAxesImg );
    
    aol::Vec2<int> pos1, pos2, pos3;
    pos1.set ( origin[0], origin[1] );
    for ( int i=-Data.getNumX ( ); i<Data.getNumX ( ) ; ++i ) {
      pos1.set ( round ( origin[0] + i * LatticeVectors[0][0] ), round ( origin[1] + i * LatticeVectors[0][1] ) );
      for ( int j=-Data.getNumX ( ); j<Data.getNumX ( ) ; ++j ) {
        pos2.set ( round ( pos1[0] + j * LatticeVectors[1][0] ), round ( pos1[1] + j * LatticeVectors[1][1] ) );
        for ( int dx=-searchWindowOffset; dx<=searchWindowOffset ; ++dx ) {
          for ( int dy=-searchWindowOffset; dy<=searchWindowOffset ; ++dy ) {
            pos3.set ( pos2[0] + dx, pos2[1] + dy );
            if ( pos3[0] >= 0 && pos3[0] < Data.getNumX ( ) && pos3[1] >= 0 && pos3[1] < Data.getNumY ( ) ) {
              periodicityAxesImg[0].set ( pos3, 1.0 );
              periodicityAxesImg[1].set ( pos3, 0 );
              periodicityAxesImg[2].set ( pos3, 0 );
            }
          }
        }
      }
    }
    periodicityAxesImg[0].set ( origin, 0 );
    periodicityAxesImg[1].set ( origin, 1.0 );
    periodicityAxesImg.savePNG ( Path );
  }
  
  static void saveDataPeriodicityAxesImg ( const char* Path, const PictureType &Data, const aol::Vector<RealType> &LatticeAnglesRadians,
                                           const aol::Vec2<short> &Origin = aol::Vec2<short> ( -1, -1 ) ) {
    qc::CoordType origin ( Origin[0], Origin[1] );
    if ( !aol::InsideQuad<qc::QC_2D> ( qc::CoordType ( 0, 0 ), qc::CoordType ( Data.getNumX ( ), Data.getNumY ( ) ), origin ) )
      origin = getAtomCenterReferenceCoordinate<RealType, qc::QC_2D> ( Data, qc::CoordType ( 0, 0 ), qc::CoordType ( Data.getNumX ( ), Data.getNumY ( ) ) );
    
    ColoredPictureType periodicityAxesImg;
    getPNGYellowNaNs ( Data, periodicityAxesImg );
    aol::Vec2<short> pos;
    for ( short k=0; k<LatticeAnglesRadians.size ( ) ; ++k ) {
      for ( short i=-Data.getNumX ( ); i<Data.getNumX ( ) ; ++i ) {
        pos.set ( origin[0] + cos ( LatticeAnglesRadians[k] ) * i, origin[1] + sin ( LatticeAnglesRadians[k] ) * i );
        if ( pos[0] >= 0 && pos[0] < Data.getNumX ( ) && pos[1] >= 0 && pos[1] < Data.getNumY ( ) ) {
          periodicityAxesImg[0].set ( pos, 0 );
          periodicityAxesImg[1].set ( pos, 0 );
          periodicityAxesImg[2].set ( pos, 1.0 );
        }
      }
    }
    periodicityAxesImg.savePNG ( Path );
  }
  
  static void saveDataLatticeVectorsFigure ( const char* BasePath, const PictureType &Data, const aol::MultiVector<RealType> &LatticeVectors,
                                             const aol::Vec2<short> &Origin = aol::Vec2<short> ( -1, -1 ) ) {
    qc::CoordType origin ( Origin[0], Origin[1] );
    if ( !aol::InsideQuad<qc::QC_2D> ( qc::CoordType ( 0, 0 ), qc::CoordType ( Data.getNumX ( ), Data.getNumY ( ) ), origin ) )
      origin = getAtomCenterReferenceCoordinate<RealType, qc::QC_2D> ( Data, qc::CoordType ( 0, 0 ), qc::CoordType ( Data.getNumX ( ), Data.getNumY ( ) ) );
    
    std::string backGroundFileName = aol::strprintf ( "%s_gnuBackground.png", BasePath );
    ColoredPictureType background;
    getPNGYellowNaNs ( Data, background );
    background.savePNG ( backGroundFileName.c_str ( ) );
    
    std::string outLatVecFileName = aol::strprintf ( "%s_gnuLatVec.dat", BasePath );
    std::ofstream outLatVec;
    outLatVec.open ( outLatVecFileName.c_str ( ) );
    for ( int k = 0; k < 2; ++k ) {
      outLatVec << origin[0] << " " << Data.getNumY ( )-1-origin[1] << " ";
      outLatVec << LatticeVectors[k][0] << " " << -LatticeVectors[k][1] << endl;
    }
    outLatVec.close ( );
    std::vector<std::string> latVecFiles; latVecFiles.push_back ( outLatVecFileName );
    aol::plotUnitCells<RealType, PictureType> ( aol::strprintf ( "%s_gnuplot.dat", BasePath ), aol::strprintf ( "%s.eps", BasePath ), backGroundFileName, latVecFiles );
  }
protected:
  void getPrimitiveUnitCellVectors ( aol::MultiVector<RealType> &LatticeVectors,
                                     const aol::Vector<RealType> &FundamentalPeriods, const aol::Vector<RealType> &LatticeAngles,
                                     const PictureType &Data ) {
    aol::Vector<RealType> fundamentalPeriods ( FundamentalPeriods );
    aol::Vector<RealType> latticeAnglesMinPeriod ( 2 );
    for ( int i=0; i<2 ; ++i ) {
      std::pair<int, RealType> minIndVal = fundamentalPeriods.getMinIndexAndValue ( );
      
      latticeAnglesMinPeriod[i] = LatticeAngles[minIndVal.first];
      LatticeVectors[i][0] = minIndVal.second * cos ( latticeAnglesMinPeriod[i] );
      LatticeVectors[i][1] = minIndVal.second * sin ( latticeAnglesMinPeriod[i] );
      
      fundamentalPeriods[minIndVal.first] = fundamentalPeriods.getMaxValue ( ) + 1;
    }
    
    if ( !_quietMode ) std::cerr << "Primitive unit cell vectors:" << std::endl << LatticeVectors << std::endl;
    if ( _outputDir.size ( ) > 0 ) {
      saveDataPeriodicityAxesImg ( aol::strprintf ( "%s/latticeSymmetryAxes.png", _outputDir.c_str ( ) ).c_str ( ), Data, latticeAnglesMinPeriod );
      saveDataPeriodicPatternImg ( aol::strprintf ( "%s/crystalLattice.png", _outputDir.c_str ( ) ).c_str ( ), Data, LatticeVectors );
      saveDataLatticeVectorsFigure ( aol::strprintf ( "%s/latticeVectors", _outputDir.c_str ( ) ).c_str ( ), Data, LatticeVectors );
      LatticeVectors.saveToFile ( aol::strprintf ( "%s/latticeVectors.dat", _outputDir.c_str ( ) ).c_str ( ) );
    }
  }
    
  void getLatticeAnglesFromProjectiveStdDevPeaks ( aol::Vector<RealType> &LatticeAngles, const PictureType &Data ) {
    LatticeAngles.resize ( 0 );
    
    // Calculate projective standard deviations for each angle with 1 degree increments
    const int l0 = 0.1 * sqrt ( Data.size ( ) );
    aol::Vector<RealType> projectiveStandardDeviations ( 180 ), projectedAverageIntensities;
    for ( int angleDegrees=0; angleDegrees<180 ; ++angleDegrees ) {
      getProjectedAverageIntensities ( projectedAverageIntensities, Data, angleDegrees, l0 ) ;
      projectiveStandardDeviations[angleDegrees] = projectedAverageIntensities.getStdDev ( );
    }
    projectiveStandardDeviations.addToAll ( -projectiveStandardDeviations.getMinValue ( ) );
    aol::Vector<RealType> projectiveStandardDeviationsNoisy ( projectiveStandardDeviations );
    getMeanFilteredVectorPeriodic<RealType> ( projectiveStandardDeviationsNoisy, projectiveStandardDeviations );
    
#ifdef USE_MODULES_QT
    if ( _outputDir.size ( ) > 0 ) {
      CustomPlotHandler<RealType> qcpHandler ( "Projective standard deviation" );
      qcpHandler.setAxesLabels ( "angle", "psd" );
      qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( projectiveStandardDeviations );
      qcpHandler.saveToFile ( aol::strprintf ( "%s/psd.q1cp", _outputDir.c_str ( ) ).c_str ( ) );
    }
#endif
    
    // Find and refine all significant peaks of projective standard deviation
    const RealType threshold = 0.1 * projectiveStandardDeviations.getMaxValue ( );
    RealType angle;
    // Keep extracting significant peaks until no proper candidates are left
    while ( LatticeAngles.size ( ) < _maxNumPSDPeaks && projectiveStandardDeviations.getMaxValue ( ) > 0 ) {
      if ( getAndRemovePeakAndAngle ( Data, angle, projectiveStandardDeviations, LatticeAngles.size ( ), threshold ) )
        LatticeAngles.pushBack ( to180Degrees ( angle ) * aol::NumberTrait<RealType>::pi / 180.0 );
    }
  }
  
  RealType to180Degrees ( const RealType Degrees ) const {
    RealType res = Degrees;
    while ( res < 0 ) res += 180;
    while ( res >= 180 ) res -= 180;
    return res;
  }
  
  RealType to360Degrees ( const RealType Degrees ) const {
    RealType res = Degrees;
    while ( res < 0 ) res += 360;
    while ( res >= 360 ) res -= 360;
    return res;
  }
  
  bool getAndRemovePeakAndAngle ( const PictureType &Data, RealType &Angle, aol::Vector<RealType> &ProjectiveStandardDeviations, const int I,
                                  const RealType Threshold = 0 ) const {
    if ( ProjectiveStandardDeviations.getMaxValue ( ) == 0 ) return false;
    
    // Find peak
    std::pair<int, RealType> indVal = ProjectiveStandardDeviations.getMaxIndexAndValue ( );
    ProjectiveStandardDeviations[indVal.first] = 0;
    const RealType localMinVal = removePeak ( indVal.first, ProjectiveStandardDeviations );
    if ( aol::Abs<RealType> ( indVal.second - localMinVal ) < Threshold )
      return false;
    
    if ( !_quietMode ) std::cerr << "#" << I+1 << ": peak = " << indVal.first;
    
    // Refine peak
    Angle = getLocalMaxProjectiveStandardDeviationAngle ( Data, indVal.first, 1.0 );
    if ( !_quietMode ) std::cerr << "; fitted angle = " << Angle << std::endl;
    removePeak ( Angle, ProjectiveStandardDeviations );
    
#ifdef USE_MODULES_QT
    if ( _outputDir.size ( ) > 0 ) {
      CustomPlotHandler<RealType> qcpHandler ( aol::strprintf ( "Projective standard deviation (-peak %d)", I+1 ) );
      qcpHandler.setAxesLabels ( "angle", "psd" );
      qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( ProjectiveStandardDeviations );
      qcpHandler.saveToFile ( aol::strprintf ( "%s/psd_-peak%d.q1cp", _outputDir.c_str ( ), I+1 ).c_str ( ) );
    }
#endif
    
    Angle = to180Degrees ( Angle + 90 );
    
    return true;
  }
  
  RealType removePeak ( const RealType Angle, aol::Vector<RealType> &ProjectiveStandardDeviations ) const {
    int k = 0;
    int angleDegrees = round ( Angle )-1;
    aol::Vec2<RealType> localMinVals ( ProjectiveStandardDeviations.getMaxValue ( ), ProjectiveStandardDeviations.getMaxValue ( ) );
    while ( k < 3 ) {
      ++k;
      if ( ProjectiveStandardDeviations[to180Degrees ( angleDegrees - 1 )] < ProjectiveStandardDeviations[to180Degrees ( angleDegrees )] ) k = 0;
      if ( ProjectiveStandardDeviations[to180Degrees ( angleDegrees )] > 0 )
        localMinVals[0] = aol::Min<RealType> ( localMinVals[0], ProjectiveStandardDeviations[to180Degrees ( angleDegrees )] );
      ProjectiveStandardDeviations[to180Degrees ( angleDegrees )] = 0;
      --angleDegrees;
    }
    k = 0;
    angleDegrees = round ( Angle )+1;
    while ( k < 3 ) {
      ++k;
      if ( ProjectiveStandardDeviations[to180Degrees ( angleDegrees + 1 )] < ProjectiveStandardDeviations[to180Degrees ( angleDegrees )] ) k = 0;
      if ( ProjectiveStandardDeviations[to180Degrees ( angleDegrees )] > 0 )
        localMinVals[1] = aol::Min<RealType> ( localMinVals[1], ProjectiveStandardDeviations[to180Degrees ( angleDegrees )] );
      ProjectiveStandardDeviations[to180Degrees ( angleDegrees )] = 0;
      ++angleDegrees;
    }
    
    return localMinVals.getMaxValue ( );
  }
  
  RealType getLocalMaxProjectiveStandardDeviationAngle ( const PictureType &Data, const RealType Angle, const RealType Delta, const int I = -1 ) const {
    const int l0 = 0.1 * sqrt ( Data.size ( ) );
    const int n = 20;
    const RealType minAngle = Angle - Delta, maxAngle = Angle + Delta;
    aol::Vector<RealType> psd ( n ), angleDegrees ( n ), projectedAverageIntensities;
    for ( int i=0; i<n ; ++i ) {
      angleDegrees[i] = minAngle + i / static_cast<RealType> ( n - 1 ) * ( maxAngle - minAngle );
      getProjectedAverageIntensities ( projectedAverageIntensities, Data, angleDegrees[i], l0 ) ;
      psd[i] = projectedAverageIntensities.getStdDev ( );
    }
    aol::Vector<RealType> psdNoisy ( psd );
    getMeanFilteredVector<RealType> ( psdNoisy, psd, false );
    
#ifdef USE_MODULES_QT
    if ( I >= 0 && _outputDir != "" ) {
      CustomPlotHandler<RealType> qcpHandler ( "Refined projective standard deviation" );
      qcpHandler.setAxesLabels ( "angle", "psd" );
      qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( angleDegrees, psd );
      qcpHandler.saveToFile ( aol::strprintf ( "%s/psd_refined%d.q1cp", _outputDir.c_str ( ), I+1 ).c_str ( ) );
    }
#else
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( I );
#endif
    
    return angleDegrees[psd.getMaxIndexAndValue ( ).first];
  }
  
  void fitGaussianToProjectiveStandardDeviationPeak ( const PictureType &Data, RealType &Angle, const int I ) const {
    const int l0 = 0.1 * sqrt ( Data.size ( ) );
    const int n = 50;
    const RealType minAngle = Angle - 2.5, maxAngle = Angle + 2.5;
    aol::Vector<RealType> psd ( n ), angleDegrees ( n ), projectedAverageIntensities;
    for ( int i=0; i<n ; ++i ) {
      angleDegrees[i] = minAngle + i / static_cast<RealType> ( n - 1 ) * ( maxAngle - minAngle );
      getProjectedAverageIntensities ( projectedAverageIntensities, Data, angleDegrees[i], l0 ) ;
      psd[i] = projectedAverageIntensities.getStdDev ( );
    }
    aol::Vector<RealType> psdNoisy ( psd );
    getMeanFilteredVectorPeriodic<RealType> ( psdNoisy, psd );
    
    std::vector<std::pair<RealType, RealType> > dataPairs;
    for ( int i=0; i<n ; ++i ) dataPairs.push_back ( std::pair<RealType, RealType> ( angleDegrees[i], psd[i] ) );
    Gaussian1DTargetFunctional<RealType> F ( dataPairs );
    Gaussian1DTargetJacobian<RealType, aol::FullMatrix<RealType> > DF ( dataPairs );
    aol::LevenbergMarquardtAlgorithm<RealType, aol::FullMatrix<RealType>, aol::LinearRegressionQR<RealType> > levenbergMarquardtAlg ( dataPairs.size ( ),
                                                                                                                                      F, DF, 50, 1, 0.2, 0.8, 1e-6, 1e-6, 1e-6, false );
    aol::Vector<RealType> arg ( 3 ), dest ( 3 );
    arg[0] = Angle;
    arg[1] = 1.0;
    arg[2] = psd.getMaxValue ( );
    levenbergMarquardtAlg.apply ( arg, dest );
    
    Angle = dest[0];
  }
  
  void getProjectedAverageIntensities ( aol::Vector<RealType> &ProjectedAverageIntensities, const PictureType &Data, const RealType AngleDegrees, const int L0 ) const {
    RealType angleRadians = AngleDegrees * aol::NumberTrait<RealType>::pi / 180.0;
    std::map<int, RealType> aDelta, nDelta;
    
    // Project image onto line with origin (0,0) and angle angleRadians
    for ( int y=0; y<Data.getNumY ( ) ; ++y ) {
      for ( int x=0; x<Data.getNumX ( ) ; ++x ) {
        if ( !aol::isNaN<RealType> ( Data.get ( x, y ) ) ) {
          const int p = floor ( getProjectedPosition ( x, y, angleRadians ) );
          nDelta[p] = nDelta[p] + 1;
          aDelta[p] = aDelta[p] + Data.get ( x, y );
        }
      }
    }
    
    // Normalize and threshold bins, then convert projected average intensities to vector
    ProjectedAverageIntensities.resize ( 0 );
    for ( typename std::map<int, RealType>::iterator it=aDelta.begin ( ); it != aDelta.end ( ); ++it ) {
      if ( nDelta[it->first] >= L0 ) aDelta[it->first] = aDelta[it->first] / nDelta[it->first];
      else aDelta[it->first] = 0;
      ProjectedAverageIntensities.pushBack ( aDelta[it->first] );
    }
    RealType meanVal = ProjectedAverageIntensities.getMeanValue ( );
    for ( int i=0; i<ProjectedAverageIntensities.size ( ) ; ++i ) {
      if ( ProjectedAverageIntensities[i] == 0 ) ProjectedAverageIntensities[i] = meanVal;
    }
  }
  
  RealType getProjectedPosition ( const int X, const int Y, const RealType AngleRadians ) const {
    return X * cos ( AngleRadians ) + Y * sin ( AngleRadians );
  }
  
    
  void getFundamentalPeriodsFromPeriodicityEnergyMinimization ( aol::Vector<RealType> &FundamentalPeriods, aol::Vector<RealType> &LatticeAngles, const PictureType &Data ) {
    if ( _outputDir.size ( ) > 0 ) {
      qc::MultiArray<RealType, qc::QC_2D, 3> dataMulti;
      getPNGYellowNaNs ( Data, dataMulti );
      dataMulti.savePNG ( aol::strprintf ( "%s/dataYellowNaNs.png", _outputDir.c_str ( ) ).c_str ( ) );
    }
    
    FundamentalPeriods.reallocate ( LatticeAngles.size ( ) );
#ifdef USE_MODULES_QT
    CustomPlotHandler<RealType> qcpHandler ( "Periodicity energies", true );
    qcpHandler.setAxesLabels ( "period", "energy" );
#endif
    for ( int i=0; i<LatticeAngles.size ( ) ; ++i ) {
      // Compute energies of difference between image and shifted image for a range of possible periods
      aol::Vector<RealType> energies;
      for ( int period=0; period<0.5*aol::Min<int> ( Data.getNumX ( ), Data.getNumY ( ) ) ; ++period ) energies.pushBack ( periodicityEnergy1D ( Data, LatticeAngles[i], static_cast<RealType> ( period ) ) );
      energies /= energies.getMaxValue ( );
      aol::Vector<RealType> preprocessedEnergies;
      getMeanFilteredVector ( energies, preprocessedEnergies, false );
      
#ifdef USE_MODULES_QT
      if ( _outputDir.size ( ) > 0 )
        qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( preprocessedEnergies,
                                                                         aol::strprintf ( "angle = %f", LatticeAngles[i] * 180 / aol::NumberTrait<RealType>::pi ).c_str ( ),
                                                                         QuocQCPGraphStyle<RealType> ( QuocQCPGraphStyle<RealType> ( quocQCPColorName ( i ).toStdString ( ) ) ) );
#endif
      
      // Perform periodicity analysis on energies
      
      // Step 1: Find all local minima
      aol::Vector<RealType> localMinima;
      aol::Vector<int> localMinimaSpacings;
      for ( int x=1; x<preprocessedEnergies.size ( )-1 ; ++x ) {
        if ( preprocessedEnergies[x] < preprocessedEnergies[x-1] && preprocessedEnergies[x] < preprocessedEnergies[x+1] ) {
          localMinima.pushBack ( preprocessedEnergies[x] );
          localMinimaSpacings.pushBack ( x );
        }
      }
      
      if ( !_quietMode ) {
        std::cerr << "Local Minima: " << localMinima << std::endl;
        std::cerr << "Spacings: " << localMinimaSpacings << std::endl;
      }
      
      if ( localMinima.size ( ) > 0 ) {
        // Step 2: Throw away all local minima whose energy is not below the threshold (absolute minimum + delta)
        const RealType delta = 0.2;
        const RealType threshold = aol::Min<RealType> ( localMinima.getMinValue ( ) + delta, 0.5 * preprocessedEnergies.getMaxValue ( ) );
        int j=0;
        while ( j < localMinima.size ( ) ) {
          if ( localMinima[j] > threshold ) {
            localMinima.erase ( j );
            localMinimaSpacings.erase ( j );
          } else
            ++j;
        }
        
        // Step 3: Determine local minimum that corresponds to smallest spacing and which belongs to the same class of local minima as the global minimum
        if ( localMinima.size ( ) == 1 ) {
          FundamentalPeriods[i] = localMinimaSpacings[0];
        } else {
          if ( !_quietMode ) std::cerr << "Local minima #" << i+1 << ": " << localMinima << std::endl;
          
          // Step 3.1 Linear correction of energy penalty due to errors in angle estimation
          aol::Vec2<RealType> linearFitParams;
          std::vector<std::pair<RealType, RealType> > linearFitData;
          aol::Vector<RealType> localMinimaTmp ( localMinima );
          while ( true ) {
            std::pair<int, RealType> minIndVal = localMinimaTmp.getFiniteMinIndexAndValue ( );
            linearFitData.push_back ( std::pair<RealType, RealType> ( localMinimaSpacings[minIndVal.first], localMinimaTmp[minIndVal.first] ) );
            if ( minIndVal.first == localMinimaTmp.size ( ) - 1 ) break;
            for ( int j=0; j<=minIndVal.first ; ++j ) localMinimaTmp[j] = aol::NumberTrait<RealType>::NaN;
          }
          if ( linearFitData.size ( ) >= 4 ) {
            getLinearFit ( linearFitData, linearFitParams );
            std::pair<int, RealType> minIndVal = localMinima.getMinIndexAndValue ( );
            aol::Vector<RealType> constraints;
            for ( int j = minIndVal.first + 1; j<localMinima.size ( ); ++j ) constraints.pushBack ( ( localMinima[j] - minIndVal.second ) / ( localMinimaSpacings[j] - localMinimaSpacings[minIndVal.first] ) );
            const RealType slope = aol::Clamp<RealType> ( linearFitParams[0], 0, 0.99 * constraints.getMinValue ( ) );
            for ( int j = minIndVal.first + 1; j<localMinima.size ( ); ++j ) localMinima[j] -= slope * ( localMinimaSpacings[j] - localMinimaSpacings[minIndVal.first] );
            
#ifdef USE_MODULES_QT
            if ( _outputDir.size ( ) > 0 ) {
              aol::Vector<RealType> preprocessedEnergiesLinCorr ( preprocessedEnergies.size ( ) );
              for ( int j=0; j<preprocessedEnergies.size ( ) ; ++j ) preprocessedEnergiesLinCorr[j] = preprocessedEnergies[j] - slope * ( j - localMinimaSpacings[minIndVal.first] );
              qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( preprocessedEnergiesLinCorr,
                                                                               aol::strprintf ( "angle = %f (lin. corr.)", LatticeAngles[i] * 180 / aol::NumberTrait<RealType>::pi ).c_str ( ),
                                                                               QuocQCPGraphStyle<RealType> ( quocQCPColorName ( i ).toStdString ( ), -1, 2 ) );
            }
#endif
          }
          
          // Step 3.2: Cluster local Minima using k-Means and determining k based on Akaike's Information Criterion
          aol::PMeansClusterer<RealType> pMeansClusterer ( localMinima.size ( ) / 2, 0.05 );
          pMeansClusterer.setQuietMode ( _quietMode );
          aol::Vector<RealType> clusters;
          aol::Vector<int> clusterLabels;
          pMeansClusterer.apply ( localMinima, clusters, clusterLabels, aol::AIC, aol::MIN_NUMCLUSTERS_SIGNIFICANT_CRITERIA );
          if ( !_quietMode ) {
            std::cerr << "Significant relative AIC criterion identified " << clusters.size ( ) << " cluster(s) for energy " << i+1 << std::endl;
            std::cerr << "Clusters #" << i+1 << ": " << clusters << std::endl;
          }
          
          // Step 3.3: Select cluster with smallest mean (i.e. the one corresponding to the global minimum)
          std::pair<int, RealType> minClustersIndVal = clusters.getMinIndexAndValue ( );
          aol::Vector<RealType> clusterSpacings, clusterEnergies;
          for ( int j=0; j<localMinima.size ( ) ; ++j ) {
            if ( clusterLabels[j] == minClustersIndVal.first ) {
              clusterSpacings.pushBack ( localMinimaSpacings[j] );
              clusterEnergies.pushBack ( localMinima[j] );
            }
          }
          
          // Step 3.4: From cluster with smallest mean, select local minimum that corresponds to smallest period among all that are within the threshold of the absolute minimum
          std::pair<int, RealType> minClusterIndVal = clusterEnergies.getMinIndexAndValue ( );
          int minSpacingIndex = minClusterIndVal.first;
          for ( int j=0; j<clusterSpacings.size ( ) ; ++j ) {
            if ( clusterEnergies[j] < minClusterIndVal.second + delta && clusterSpacings[j] < clusterSpacings[minSpacingIndex] )
              minSpacingIndex = j;
          }
          
          FundamentalPeriods[i] = clusterSpacings[minSpacingIndex];
        }
      
        // Step 4: Refine local minimum
        const int N = 50;
        aol::Vector<RealType> energiesRefined;
        aol::Vector<RealType> periodsRefined;
        for ( int t=-N; t<=N ; ++t ) {
          const RealType period = FundamentalPeriods[i] + 2 * static_cast<RealType> ( t ) / static_cast<RealType> ( N );
          const RealType energy = periodicityEnergy1D ( Data, LatticeAngles[i], period );
          energiesRefined.pushBack ( energy );
          periodsRefined.pushBack ( period );
        }
        aol::Vector<RealType> preprocessedEnergiesRefined;
        getMeanFilteredVector ( energiesRefined, preprocessedEnergiesRefined, false );
        
#ifdef USE_MODULES_QT
        if ( _outputDir.size ( ) > 0 ) {
          CustomPlotHandler<RealType> plotHandler ( aol::strprintf ( "Refined periodicity energy (angle = %f)", LatticeAngles[i] * 180 / aol::NumberTrait<RealType>::pi ) );
          plotHandler.setAxesLabels ( "period", "energy" );
          plotHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( periodsRefined, preprocessedEnergiesRefined );
          plotHandler.saveToFile ( aol::strprintf ( "%s/energies_refined%d.q1cp", _outputDir.c_str ( ), i ).c_str ( ) );
        }
#endif
        
        FundamentalPeriods[i] = periodsRefined[preprocessedEnergiesRefined.getMinIndexAndValue ( ).first];
      } else
        FundamentalPeriods[i] = aol::NumberTrait<RealType>::Inf;
      
      if ( !_quietMode ) std::cerr << "Fundamental period #" << i+1 << ": " << FundamentalPeriods[i] << std::endl;
    }
    
#ifdef USE_MODULES_QT
    if ( _outputDir.size ( ) > 0 )
      qcpHandler.saveToFile ( aol::strprintf ( "%s/energies.q1cp", _outputDir.c_str ( ) ).c_str ( ) );
#endif
  }
  
  
  void getLatticeAnglesFromFourierPeaks ( aol::Vector<RealType> &LatticeAngles, const PictureType &PreprocessedData, const PictureType &Data, bool RefineByProjectiveStdDevPeaks = true ) {
    // Find largest (possibly rotated) square that fits into image without containing NaNs and resample
    if ( _outputDir.size ( ) > 0 ) PreprocessedData.save ( aol::strprintf ( "%s/fftInput%s", _outputDir.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    PictureType data ( PreprocessedData );
    RealType angle = 0;
    rotateAndResampleFromLargestRectangle<RealType> ( data, angle );
    if ( _outputDir.size ( ) > 0 ) data.save ( aol::strprintf ( "%s/fftInput_rect%s", _outputDir.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    
    // Calculate fourier power coefficients
    const RealType mean = data.getMeanValue ( );
    qc::MultiArray<RealType, 2, 2> dataEliminatedMean ( data.getNumX ( ), data.getNumY ( ) );
    for ( short i=0; i<data.getNumX ( ) ; ++i )
      for ( short j=0; j<data.getNumY ( ) ; ++j )
        dataEliminatedMean[0].set ( i, j, data.get ( i, j ) - mean );
    qc::ScalarArray<RealType, qc::QC_2D> modulus ( data.getNumX ( ), data.getNumY ( ) );
    qc::computeLogFFTModulus<RealType> ( dataEliminatedMean[0], modulus, 0, false );
    modulus.scaleValuesTo01 ( );
    if ( _outputDir.size ( ) > 0 ) modulus.save ( aol::strprintf ( "%s/fftOutput%s", _outputDir.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    
    // Find fourier power peaks
    aol::MultiVector<RealType> peaks ( _maxNumFFTPeaks, 2 );
    PictureType fourierPowerPeaks ( modulus );
    qc::FastILexMapper<qc::QC_2D> mapper ( fourierPowerPeaks.getNumX ( ), fourierPowerPeaks.getNumY ( ) );
    std::pair<int, RealType> maxIndVal;
    aol::Vec2<RealType> peakPos, center ( ( data.getNumX ( ) - 1 ) / 2, ( data.getNumY ( ) - 1 ) / 2 );
    int peakIdx = 0;
    while ( peakIdx < peaks.numComponents ( ) ) {
      maxIndVal = fourierPowerPeaks.getMaxIndexAndValue ( );
      const int x = mapper.splitGlobalIndex ( maxIndVal.first )[0];
      const int y = mapper.splitGlobalIndex ( maxIndVal.first )[1];
      for ( short dx=-1; dx<=1 ; ++dx )
        for ( short dy=-1; dy<=1 ; ++dy )
          if ( x + dx >= 0 && x + dx < fourierPowerPeaks.getNumX ( ) && y >= 0 && y < fourierPowerPeaks.getNumY ( ) ) fourierPowerPeaks.set ( x + dx, y + dy, 0 );
      peakPos.set ( mapper.splitGlobalIndex ( maxIndVal.first )[0], mapper.splitGlobalIndex ( maxIndVal.first )[1] );
      peakPos -= center;
      bool angleTooSmall = false;
      for ( short k=0; k<peakIdx ; ++k ) {
        const RealType p1dotp2 = peakPos.dotProduct ( aol::Vec2<RealType> ( peaks[k][0], peaks[k][1] ) ) / ( peakPos.norm ( ) * peaks[k].norm ( ) );
        if ( p1dotp2 < -0.9 || p1dotp2 > 0.9 )
          angleTooSmall = true;
      }
      if ( !angleTooSmall ) {
        peaks[peakIdx][0] = peakPos[0]; peaks[peakIdx][1] = peakPos[1];
        fourierPowerPeaks[maxIndVal.first] = -1;
        ++peakIdx;
      }
    }
    
    LatticeAngles.reallocate ( peaks.numComponents ( ) );
    for ( short i=0; i<peaks.numComponents ( ) ; ++i ) {
      LatticeAngles[i] = atan2 ( static_cast<RealType> ( -peaks[i][0] ) / static_cast<RealType> ( modulus.getNumX ( ) ),
                                 static_cast<RealType> ( peaks[i][1] ) / static_cast<RealType> ( modulus.getNumY ( ) ) );
    }
    LatticeAngles.addToAll ( -angle );
    
    if ( _outputDir.size ( ) > 0 ) {
      saveFourierCoefficients ( aol::strprintf ( "%s/periodicityFourierPeaks.png", _outputDir.c_str ( ) ).c_str ( ), modulus, true,
                                true, center, peaks,
                                true, aol::Vec2<RealType> ( LatticeAngles[0], LatticeAngles[1] ) );
      peaks.saveASCII ( aol::strprintf ( "%s/fftPeakPositions.txt", _outputDir.c_str ( ) ).c_str ( ) );
      saveDataPeriodicityAxesImg ( aol::strprintf ( "%s/latticeSymmetryAxes_Fourier.png", _outputDir.c_str ( ) ).c_str ( ), Data, LatticeAngles );
    }
    
    if ( RefineByProjectiveStdDevPeaks ) {
      for ( int j=0; j<LatticeAngles.size ( ) ; ++j ) {
        LatticeAngles[j] = to180Degrees ( LatticeAngles[j] * 180.0 / aol::NumberTrait<RealType>::pi - 90.0 );
        LatticeAngles[j] = getLocalMaxProjectiveStandardDeviationAngle ( Data, LatticeAngles[j], 1.5, j );
        LatticeAngles[j] = to180Degrees ( ( LatticeAngles[j] + 90.0 ) * aol::NumberTrait<RealType>::pi / 180.0 );
      }
    }
  }
  
  void getFundamentalPeriodsFromSineFit ( aol::Vector<RealType> &FundamentalPeriods, aol::Vector<RealType> &LatticeAngles, const PictureType &Data ) {
    LatticeAngles.resize ( 2 );
    FundamentalPeriods.reallocate ( 2 );
    
    // Find brightest peak and see which of the two axes offers the larger intersection with the image
    const qc::FastILexMapper<qc::QC_2D> mapperReduced ( Data.getNumX ( ), Data.getNumY ( ) );
    aol::Vec2<short> peak ( mapperReduced.splitGlobalIndex ( Data.getMaxIndexAndValue ( ).first )[0],
                           mapperReduced.splitGlobalIndex ( Data.getMaxIndexAndValue ( ).first )[1] );
    aol::RandomAccessContainer<std::vector<std::pair<RealType, RealType> > > intensitiesContainer ( 2 );
    for ( int k=0; k<2 ; ++k ) setIntensitiesAlongAxis ( intensitiesContainer[k], Data, LatticeAngles[k], peak );
    const short kFirst = ( intensitiesContainer[0].size ( ) > intensitiesContainer[1].size ( ) ) ? 0 : 1;
    
    // Get periodicity spacing along primary axis kFirst
    RealType energy;
    FundamentalPeriods[kFirst] = getPeriodicitySpacing ( intensitiesContainer[kFirst], energy, 1 );
    
    // Move origin of secondary axis along primary axis from brightest peak in periodicity steps,
    // and extract intensities from where the intersection of secondary axis with the image is largest
    std::vector<std::pair<RealType, RealType> > intensities;
    setIntensitiesAlongAxis ( intensities, Data, LatticeAngles[1-kFirst], peak );
    short direction = 1;
    aol::Vec2<RealType> pos ( peak[0], peak[1] );
    aol::Vec2<RealType> stepVector ( direction * FundamentalPeriods[kFirst] * cos ( LatticeAngles[kFirst] ),
                                    direction * FundamentalPeriods[kFirst] * sin ( LatticeAngles[kFirst] ) );
    unsigned short maxIntensitiesSize = 0;
    aol::Vec2<short> maxIntersectionOrigin;
    while ( true ) {
      pos += stepVector;
      if ( pos[0] < 1 || pos[0] >= Data.getNumX ( ) -1 || pos[1] < 1 || pos[1] >= Data.getNumY ( ) - 1 ) {
        if ( direction == -1 ) break;
        else {
          direction = -1;
          stepVector *= direction;
          pos.set ( peak[0], peak[1] );
          pos += stepVector;
        }
      }
      RealType localMax = 0;
      aol::Vec2<short> localPos, localMaxPos;
      for ( int dx=-1; dx<=1 ; ++dx ) {
        for ( int dy=-1; dy<=1 ; ++dy ) {
          localPos.set ( round ( pos[0] ) + dx, round ( pos[1] ) + dy );
          if ( localPos[0] >= 0 && localPos[0] < Data.getNumX ( ) && localPos[1] >= 0 && localPos[1] < Data.getNumY ( )
              && Data.get ( localPos ) > localMax ) {
            localMax = Data.get ( round ( pos[0] ) + dx, round ( pos[1] ) + dy );
            localMaxPos.set ( localPos );
          }
        }
      }
      pos.set ( localMaxPos[0], localMaxPos[1] );
      setIntensitiesAlongAxis ( intensities, Data, LatticeAngles[1-kFirst], localMaxPos );
      if ( intensities.size ( ) > maxIntensitiesSize ) {
        maxIntensitiesSize = intensities.size ( );
        maxIntersectionOrigin.set ( localMaxPos );
      }
    }
    setIntensitiesAlongAxis ( intensities, Data, LatticeAngles[1-kFirst], maxIntersectionOrigin );
    
    // Get periodicity spacing along secondary axis 1-kFirst
    FundamentalPeriods[1-kFirst] = getPeriodicitySpacing ( intensities, energy, 2 );
    
    if ( !_quietMode ) std::cerr << "Fundamental periods: " << std::endl << FundamentalPeriods << std::endl;
    if ( _outputDir.size ( ) > 0 ) {
      ofstream txtFile ( aol::strprintf ( "%s/periodicityAnalysis.txt", _outputDir.c_str ( ) ).c_str ( ) );
      txtFile << "Estimated grid parameters" << std::endl;
      txtFile << std::endl;
      txtFile << "Delta x_1 = " << FundamentalPeriods[0] << " pixels" << std::endl;
      txtFile << "Delta x_2 = " << FundamentalPeriods[1] << " pixels" << std::endl;
      txtFile << std::endl;
      txtFile << "alpha_1 = " << LatticeAngles[0] * 180.0 / aol::NumberTrait<RealType>::pi << " degrees" << std::endl;
      txtFile << "alpha_2 = " << LatticeAngles[1] * 180.0 / aol::NumberTrait<RealType>::pi << " degrees" << std::endl;
      txtFile << "alpha_1 = " << LatticeAngles[0] << " radians" << std::endl;
      txtFile << "alpha_2 = " << LatticeAngles[1] << " radians" << std::endl;
      txtFile.close ( );
    }
  }
  
  void setIntensitiesAlongAxis ( std::vector<std::pair<RealType, RealType> > &Intensities,
                                const PictureType &Data, const RealType AngleRadians,
                                const aol::Vec2<short> &Origin = aol::Vec2<short> ( -1, -1 ) ) {
    Intensities.clear ( );
    const qc::FastILexMapper<qc::QC_2D> mapper ( Data.getNumX ( ), Data.getNumY ( ) );
    
    aol::Vec2<short> origin ( Origin );
    if ( origin[0] < 0 || origin[0] >= Data.getNumX ( ) || origin[1] < 0 || origin[1] >= Data.getNumY ( ) )
      origin.set ( mapper.splitGlobalIndex ( Data.getMaxIndexAndValue ( ).first )[0], mapper.splitGlobalIndex ( Data.getMaxIndexAndValue ( ).first )[1] );
    
    aol::Vec2<short> pos;
    for ( short sign=-1; sign<=1 ; sign+=2 ) {
      int i = ( sign == -1 ) ? 1 : 0;
      pos.set ( origin[0], origin[1] );
      while ( pos[0] >= 0 && pos[0] < Data.getNumX ( ) && pos[1] >= 0 && pos[1] < Data.getNumY ( ) ) {
        Intensities.insert ( ( sign == -1 ) ? Intensities.begin ( ) : Intensities.end ( ),
                            std::pair<RealType, RealType> ( sign * i * aol::Vec2<RealType> ( cos ( AngleRadians ), sin ( AngleRadians ) ).norm ( ),
                                                           Data.get ( pos ) ) );
        ++i;
        pos.set ( origin[0] + sign * cos ( AngleRadians ) * i, origin[1] + sign * sin ( AngleRadians ) * i );
      }
    }
  }
  
  RealType getPeriodicitySpacing ( const std::vector<std::pair<RealType, RealType> > &Intensities, RealType &Energy, const int OutputNr = 0 ) {
    // Calculate mean value and create intensities vector with zero mean
    aol::Vector<RealType> intensities;
    for ( unsigned int i=0; i<Intensities.size ( ) ; ++i )
      intensities.pushBack ( Intensities[i].second );
    const RealType mean = intensities.getMeanValue ( ), maxZeroMean = intensities.getMaxValue ( ) - mean;
    std::vector<std::pair<RealType, RealType> > intensitiesZeroMean;
    for ( unsigned int i=0; i<Intensities.size ( ) ; ++i )
      intensitiesZeroMean.push_back ( std::pair<RealType, RealType> ( Intensities[i].first, Intensities[i].second - mean ) );
    
    // Define energy, derivative and second derivative of a non-linear least squares sum of sines fit target functional
    const int numTerms = 1;
    SumOfSinesTargetFunctional<RealType> F ( intensitiesZeroMean, numTerms );
    SumOfSinesTargetJacobian<RealType, aol::FullMatrix<RealType> > DF ( intensitiesZeroMean, numTerms );
    
    // Define a Levenberg-Marquardt method to find optimal parameters (especially the frequency of the sin functions)
    aol::LevenbergMarquardtAlgorithm<RealType, aol::FullMatrix<RealType>, aol::LinearRegressionQR<RealType> > levenbergMarquardtAlg ( intensitiesZeroMean.size ( ), F, DF, 50, 1, 0.2, 0.8, 1e-6, 1e-6, 1e-6, false );
    aol::Vector<RealType> arg ( 3 * numTerms ), dest ( 3 * numTerms );
    arg[0] = maxZeroMean;
    
    // Try different initial values (trust region method does not seem to converge globally)
    aol::Vector<RealType> energies, frequencies, optDest ( dest.size ( ) );
    aol::Vector<RealType> diffs ( intensitiesZeroMean.size ( ) );
    RealType energy = 0;
    
    setCtrlCHandler ( );
    if ( _progressBar != NULL ) _progressBar->setText ( aol::strprintf ( "PatternAnalysis: computing spacing along axis #%d (step %d/2)", OutputNr, OutputNr ).c_str() );
    if ( _progressBar != NULL ) _progressBar->start ( intensitiesZeroMean.size ( ) - 1 );
    for ( unsigned int dx=1; dx<intensitiesZeroMean.size ( ) && !wantsInterrupt ( ) ; ++dx ) {
      arg[1] = 2 * aol::NumberTrait<RealType>::pi / dx;
      levenbergMarquardtAlg.apply ( arg, dest );
      F.apply ( dest, diffs );
      energy = diffs.normSqr ( );
      energies.pushBack ( energy );
      frequencies.pushBack ( dest[1] );
      if ( energy == energies.getMinValue ( ) )
        optDest = dest;
      if ( _progressBar != NULL ) (*_progressBar)++;
    }
    if ( _progressBar != NULL ) _progressBar->finish ( );
    unsetCtrlCHandler ( );
    
    if ( !_quietMode ) {
      std::cerr << "Optimal initial frequency=" << 2 * aol::NumberTrait<RealType>::pi / ( energies.getMinIndexAndValue ( ).first + 1)
      << "; optimal frequency=" << frequencies[energies.getMinIndexAndValue ( ).first]
      << "; minimal energy=" << energies[energies.getMinIndexAndValue ( ).first] << std::endl;
    }
    
    RealType periodicitySpacing = aol::Abs<RealType> ( 2 * aol::NumberTrait<RealType>::pi / frequencies[energies.getMinIndexAndValue ( ).first] );
    
#ifdef USE_MODULES_QT
    if ( _outputDir.size ( ) > 0 ) {
      std::vector<std::pair<RealType, RealType> > sumOfSines;
      for ( unsigned int i=0; i<intensitiesZeroMean.size ( ) ; ++i )
        sumOfSines.push_back ( std::pair<RealType, RealType> ( intensitiesZeroMean[i].first,
                                                              SumOfSinesLeastSquaresEnergy<RealType>::sumOfSines ( intensitiesZeroMean[i].first, optDest, numTerms ) ) );
      CustomPlotHandler<RealType> qcpHandler ( "Sum of sines fit", true );
      qcpHandler.setAxesLabels ( "distance", "intensity" );
      qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( intensitiesZeroMean, "normalized image intensities" );
      qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( sumOfSines, "fitted sum of sines" );
      qcpHandler.saveToFile ( aol::strprintf ( "%s/sumOfSinesFit_%d.q1cp", _outputDir.c_str ( ), OutputNr ).c_str ( ) );
    }
#endif
    
    if ( !_quietMode )
      std::cerr << "Periodicity spacing: " << periodicitySpacing << " pixels." << std::endl;
    
    Energy = energies.getMinValue ( );
    return periodicitySpacing;
  }
  
  
  
  void saveFourierCoefficients ( const char* Path, const PictureType &Coefficients, const bool LogScale = false,
                                 const bool RedPeaks = false, const aol::Vec2<RealType> &Center = aol::Vec2<RealType> ( -1, -1 ), const aol::MultiVector<RealType> &Peaks = aol::MultiVector<RealType> ( ),
                                 const bool BlueAxes = false, const aol::Vec2<RealType> &LatticeAngles = aol::Vec2<RealType> ( -1, -1 ) ) const {
    PictureType u ( Coefficients );
    RealType uMax = u.getMaxValue ( );
    if ( LogScale ) {
      for ( int k = 0; k<u.size ( ) ; ++k )
        u[k] = log ( 1 + u[k] / uMax * 255 );
      uMax = u.getMaxValue ( );
    }
    
    if ( RedPeaks ) {
      ColoredPictureType v ( u.getNumX ( ), u.getNumY ( ) );
      v[0] = u;
      v[1] = u;
      v[2] = u;
      for ( int k=0; k<Peaks.numComponents ( ) ; ++k ) {
        v[0].set ( round ( Peaks[k][0] + Center[0] ), round ( Peaks[k][1] + Center[1] ), uMax );
        v[1].set ( round ( Peaks[k][0] + Center[0] ), round ( Peaks[k][1] + Center[1] ), 0 );
        v[2].set ( round ( Peaks[k][0] + Center[0] ), round ( Peaks[k][1] + Center[1] ), 0 );
      }
      v.setOverflowHandling ( aol::CLIP_THEN_SCALE, v.getMinValue ( ), v.getMaxValue ( ) );
      v.savePNG ( Path );
    } else if ( BlueAxes ) {
      ColoredPictureType v ( u.getNumX ( ), u.getNumY ( ) );
      v[0] = u;
      v[1] = u;
      v[2] = u;
      aol::Vec2<short> pos;
      for ( short k=0; k<2 ; ++k ) {
        for ( short i=-Coefficients.getNumX ( ); i<Coefficients.getNumX ( ) ; ++i ) {
          pos.set ( Center[0] + cos ( LatticeAngles[k] ) * i, Center[1] + sin ( LatticeAngles[k] ) * i );
          if ( pos[0] >= 0 && pos[0] < Coefficients.getNumX ( ) && pos[1] >= 0 && pos[1] < Coefficients.getNumY ( ) ) {
            v[0].set ( pos, 0 );
            v[1].set ( pos, 0 );
            v[2].set ( pos, uMax );
          }
        }
      }
      v.setOverflowHandling ( aol::CLIP_THEN_SCALE, v.getMinValue ( ), v.getMaxValue ( ) );
      v.savePNG ( Path );
    } else u.save ( Path, qc::PGM_DOUBLE_BINARY );
  }
  
  void saveFourierCoefficientsRedPeaks ( const char* Path, const PictureType Coefficients,
                                         const aol::Vec2<RealType> &Center, const aol::MultiVector<RealType> &Peaks,
                                         const bool LogScale = false ) const {
    aol::Vec2<RealType> latticeAngles;
    saveFourierCoefficients ( Path, Coefficients, LogScale, true, Center, Peaks, false, latticeAngles );
  }
  
  void saveFourierCoefficientsBlueAxes ( const char* Path, const PictureType Coefficients,
                                         const aol::Vec2<RealType> &LatticeAngles,
                                         const bool LogScale = false ) const {
    aol::Vec2<RealType> center;
    aol::MultiVector<RealType> peaks;
    saveFourierCoefficients ( Path, Coefficients, LogScale, false, center, peaks, true, LatticeAngles );
  }
  
  void saveIntensityPlotAlongPeriodicAxis ( const char* Path, const PictureType &Data, const aol::Vec2<RealType> &LatticeAnglesRadians,
                                            const short Axis, const aol::Vec2<short> &Origin = aol::Vec2<short> ( -1, -1 ) ) const {
#ifdef USE_MODULES_QT
    if ( Axis < 0 || Axis > 1 )
      throw aol::Exception ( "Argument \"Axis\" has to be between 0 and 1!", __FILE__, __LINE__ );
    
    aol::Vec2<short> origin ( Origin );
    if ( origin[0] < 0 || origin[0] >= Data.getNumX ( ) || origin[1] < 0 || origin[1] >= Data.getNumY ( ) ) {
      qc::FastILexMapper<qc::QC_2D> mapper ( Data.getNumX ( ), Data.getNumY ( ) );
      origin.set ( mapper.splitGlobalIndex ( Data.getFiniteMaxIndexAndValue ( ).first )[0], mapper.splitGlobalIndex ( Data.getFiniteMaxIndexAndValue ( ).first )[1] );
    }
    
    std::vector<std::pair<RealType, RealType> > intensities;
    aol::Vec2<short> pos;
    for ( short sign=-1; sign<=1 ; sign+=2 ) {
      int i = ( sign == -1 ) ? 1 : 0;
      pos.set ( origin[0], origin[1] );
      while ( pos[0] >= 0 && pos[0] < Data.getNumX ( ) && pos[1] >= 0 && pos[1] < Data.getNumY ( ) ) {
        intensities.insert ( ( sign == -1 ) ? intensities.begin ( ) : intensities.end ( ),
                             std::pair<RealType, RealType> ( sign * i * aol::Vec2<RealType> ( cos ( LatticeAnglesRadians[Axis] ),
                                                                                              sin ( LatticeAnglesRadians[Axis] ) ).norm ( ),
                                                                                              Data.get ( pos ) ) );
        ++i;
        pos.set ( origin[0] + sign * cos ( LatticeAnglesRadians[Axis] ) * i, origin[1] + sign * sin ( LatticeAnglesRadians[Axis] ) * i );
      }
    }
    CustomPlotHandler<RealType> qcpHandler ( aol::strprintf ( "Intensities along line [origin = (%d,%d), angle = %f]", origin[0], origin[1], LatticeAnglesRadians[Axis] * 180 / aol::NumberTrait<RealType>::pi ) );
    qcpHandler.setAxesLabels ( "raw image intensity", "distance" );
    qcpHandler.template addPlottable<QuocQCPGraphStyle<RealType> > ( intensities );
    qcpHandler.saveToFile ( aol::strprintf ( "%s.q1cp", Path ).c_str ( ) );
#else
    std::cerr << "PatternAnalyzer: plot requires external library QCustomPlot!" << std::endl;
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( Path );
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( Data );
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( LatticeAnglesRadians );
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( Axis );
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( Origin );
#endif
  }
  
  void setCtrlCHandler () const {
    _previousCtrlCHandler = signal ( InterruptSignal, aol::ctrlCHandler );
  }
  
  void unsetCtrlCHandler () const {
    signal ( InterruptSignal, _previousCtrlCHandler );
  }
  
  bool wantsInterrupt() const {
    if (!aol::getCtrlCState())
      return false;
    else
      return true;
  }
};

  
  
template <typename ConfiguratorType>
class PiecewisePeriodicTwoPhaseMSSegmentor : public FirstOrderPrimalDualTwoPhaseMSSegmentor<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> PictureType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
protected:
  const aol::MultiVector<RealType> &_imageMVec;
  const aol::Vec2<short> _spatialDims;
  aol::MultiVector<RealType> _latticeVectors;
private:
  int _outerIterations;
  aol::ProgressBar<> *_progressBar;
  bool _catchCtrlC;
  mutable sigfunc _previousCtrlCHandler;
public:
  
  PiecewisePeriodicTwoPhaseMSSegmentor ( const typename ConfiguratorType::InitType &Initializer,
                                        const RealType Gamma,
                                        aol::MultiVector<RealType> &ImageMVec )
    : FirstOrderPrimalDualTwoPhaseMSSegmentor<ConfiguratorType> ( Initializer, Gamma ),
  _imageMVec ( ImageMVec ),
  _spatialDims ( Initializer.getNumX ( ), Initializer.getNumY ( ) ),
  _latticeVectors ( 4, 2 ),
  _outerIterations ( 5 ),
  _catchCtrlC ( false ) {
    // TODO: approximate lattice vectors
  }
  
  virtual ~PiecewisePeriodicTwoPhaseMSSegmentor ( ) { }
  
protected:
  virtual void generateIndicatorFunction ( const int IndicatorNumber, ArrayType &IndicatorFunction ) const {
    if ( ( IndicatorNumber >= 2 ) || ( IndicatorNumber < 0 ) )
      throw ( aol::OutOfBoundsException ( "PiecewisePeriodicTwoPhaseMSSegmentor only defines two indicator functions.", __FILE__, __LINE__ ) );
    
    qc::FastILexMapper<qc::QC_2D> mapper ( _spatialDims );
    const int imageDim = _imageMVec.numComponents ( );
    const RealType shift = 0.5;
    int x,y;
    for ( int i = 0; i < IndicatorFunction.size(); ++i) {
      RealType indicator = 0.;
      for ( int vecIdx=0; vecIdx<2 ; ++vecIdx ) {
        mapper.splitGlobalIndex ( i, x, y );
        x += _latticeVectors[2*(1-IndicatorNumber)+vecIdx][0];
        y += _latticeVectors[2*(1-IndicatorNumber)+vecIdx][1];
        if ( x >=0 && x < _spatialDims[0] && y >= 0 && y < _spatialDims[1] ) {
          int iShifted = mapper.getGlobalIndex ( x, y );
          for ( int j = 0; j < imageDim; j++ )
            indicator += aol::Sqr( _imageMVec[j][i] - _imageMVec[j][iShifted] );
        }
      }
      IndicatorFunction[i] = indicator + shift;
    }
  }
  virtual void updateLatticeVectors ( const aol::Vector<RealType> &CurrentSegmentation, const aol::MultiVector<RealType> &ImageMVec ) {
    if ( ImageMVec.numComponents ( ) != 1 )
      throw aol::UnimplementedCodeException( "Update of lattice vectors undefined for non-scalar images!", __FILE__, __LINE__ );
    
    aol::Vector<RealType> u ( CurrentSegmentation );
    u.addToAll ( -0.5 );
    u /= 2 * u.getMaxAbsValue ( );
    u.addToAll ( 0.5 );
    u.threshold ( 0.5, 0, 1 );
    
    PictureType segment ( _spatialDims );
    for ( int i=0; i<2 ; ++i ) {
      segment.setAll ( aol::NumberTrait<RealType>::NaN );
      for ( int k=0; k<ImageMVec[0].size ( ) ; ++k )
        if ( u[k] == 1-i ) segment[k] = ImageMVec[0][k];
      aol::MultiVector<RealType> latticeVectors ( 2 , 2 );
      for ( int vecIdx=0; vecIdx<1 ; ++vecIdx )
        for ( int vecEntryIdx=0; vecEntryIdx<2 ; ++vecEntryIdx )
          latticeVectors[vecIdx][vecEntryIdx] = _latticeVectors[2*i+vecIdx][vecEntryIdx];
      im::PatternAnalyzer<RealType, PictureType> patternAnalyzer ( "", false );
      patternAnalyzer.refineLatticeVectorsByPeriodicityFunctionalMinimization ( latticeVectors, segment );
      for ( int vecIdx=0; vecIdx<1 ; ++vecIdx )
        for ( int vecEntryIdx=0; vecEntryIdx<2 ; ++vecEntryIdx )
          _latticeVectors[2*i+vecIdx][vecEntryIdx] = latticeVectors[vecIdx][vecEntryIdx];
    }
  }
public:
  void segmentAndAdjustGrayValues ( ArrayType &Segmentation, qc::MultiArray<RealType, ConfiguratorType::Dim> * PDual = NULL ) {
    setCtrlCHandler ( );
    for ( int i = 0; i < _outerIterations && !wantsInterrupt ( ) ; ++i ) {
      this->segment ( Segmentation, PDual );
      this->updateLatticeVectors ( Segmentation );
    }
    unsetCtrlCHandler ( );
  }
  
  const aol::MultiVector<RealType>& getLatticeVectorsReference () const {
    return _latticeVectors;
  }
  
  aol::MultiVector<RealType>& getLatticeVectorsReference () {
    return _latticeVectors;
  }
  
  void setOuterIterations ( const int OuterIterations ) {
    _outerIterations = OuterIterations;
  }
  
  void setCatchCtrlC ( bool catchCtrlC ) {
    _catchCtrlC = catchCtrlC;
  }
protected:
  void setCtrlCHandler () const {
    if (_catchCtrlC)
      _previousCtrlCHandler = signal ( InterruptSignal, aol::ctrlCHandler );
  }
  
  void unsetCtrlCHandler () const {
    if (_catchCtrlC)
      signal ( InterruptSignal, _previousCtrlCHandler );
  }
  
  bool wantsInterrupt() const {
    if (!_catchCtrlC || !aol::getCtrlCState())
      return false;
    else
      return true;
  }
};
  
  
  

template <typename _RealType, typename _PictureType>
class PatternAnalyzerBergmann {
  typedef _RealType RealType;
  typedef aol::Matrix22<int> MatrixType;
  typedef _PictureType PictureType;
  typedef qc::MultiArray<RealType, qc::QC_2D, 3> ColoredPictureType;
protected:
  const PictureType &_data;
  const std::string &_outputDir;
  const bool _quietMode;
  const qc::FastILexMapper<qc::QC_2D> _mapper;
  PictureType _fourierPowerCoefficients, _fourierPowerPeaks;
  MatrixType _patternNormalForm;
  aol::RandomAccessContainer<aol::Vec2<short> > _pattern;
public:
  PatternAnalyzerBergmann ( const PictureType &Data, const std::string &OutputDir = "", const bool QuietMode = true )
  : _data ( Data ), _outputDir ( OutputDir ), _quietMode ( QuietMode ), _mapper ( Data.getNumX ( ), Data.getNumY ( ) ),
    _fourierPowerCoefficients ( Data.getNumX ( ), Data.getNumY ( ) ), _fourierPowerPeaks ( Data.getNumX ( ), Data.getNumY ( ) ),
    _patternNormalForm ( ) {
    setFourierPowerCoefficients ( );
    setFourierPowerPeaks ( );
    setPatternNormalForm ( );
    setPattern ( );
  }
  
  static MatrixType getPatternNormalForm ( const MatrixType &Matrix ) {
    if ( Matrix.det ( ) == 0 )
      throw aol::Exception ( "The specified matrix is not regular!", __FILE__, __LINE__ );
    
    MatrixType pnf ( Matrix );
    // Form upper triangular matrix
    gcdOnRows ( pnf, 1, 0 );
    // Make diagonal positive
    for ( short row=0; row<2 ; ++row ) {
      if ( pnf.get ( row, row ) < 0 ) {
        for ( short col=0; col<2 ; ++col )
          pnf.set ( row, col, -pnf.get ( row, col ) );
      }
    }
    // Make upper non-zero values of a column lie between 0 and
    int f = 0;
    if ( pnf.get ( 0, 1 ) < 0 || pnf.get ( 0, 1 ) >= pnf.get ( 1, 1 ) )
      f = -floor ( pnf.get ( 0, 1 ) / pnf.get ( 1, 1 ) );
    for ( short col=0; col<2 ; ++col )
      pnf.set ( 0, col, pnf.get ( 0, col ) + f * pnf.get ( 1, col ) );
    
    return pnf;
  }
  
  static void setPattern ( aol::RandomAccessContainer<aol::Vec2<short> > &Positions, const MatrixType &M, const aol::Vec2<short> &TorusSize = aol::Vec2<short> ( 0, 0 ) ) {
    Positions.clear ( );
    const RealType stepSize1 = 1.0 / abs ( M.get ( 1, 1 ) );
    aol::Vector<RealType> steps1;
    RealType curStep = -0.5;
    do {
      steps1.pushBack ( curStep );
      curStep += stepSize1;
    } while ( curStep < 0.5-stepSize1 );
    const RealType stepSize2 = 1.0 / abs ( M.get ( 0, 0 ) );
    aol::Vector<RealType> steps2;
    curStep = -0.5;
    do {
      steps2.pushBack ( curStep );
      curStep += stepSize2;
    } while ( curStep < 0.5-stepSize2 );
    aol::Vector<RealType> tSums ( steps1 );
    tSums *= M.get ( 0, 1 );
    aol::RandomAccessContainer<aol::Vec2<RealType> > positions;
    for ( int i=0; i<steps1.size ( ) ; ++i ) {
      aol::Vec2<RealType> newPos;
      for ( int j=abs ( M.get ( 0, 0 ) ) * i; j<abs ( M.get ( 0, 0 ) ) * ( i + 1 ) ; ++j ) {
        newPos[0] = steps2[j - abs ( M.get ( 0, 0 ) ) * i] + stepSize2 * ( ceil ( tSums[i] ) - tSums[i] );
        newPos[1] = steps1[i];
        positions.pushBack ( newPos );
      }
    }
    aol::Vec2<short> torusSize ( TorusSize );
    if ( torusSize.norm ( ) == 0 )
      torusSize.set ( M.det ( ), M.det ( ) );
    for ( int i=0; i<positions.size ( ) ; ++i ) {
      int x = round ( positions[i][0] * torusSize[0] ), y = round ( positions[i][1] * torusSize[1] );
      x %= static_cast<int> ( torusSize[0] );
      y %= static_cast<int> ( torusSize[1] );
      x = ( x < 0 ) ? x + torusSize[0] : x;
      y = ( y < 0 ) ? y + torusSize[1] : y;
      Positions.pushBack ( aol::Vec2<short> ( x , y ) );
    }
  }
  
  void savePatternDataImage ( const aol::Vec2<short> &Point = aol::Vec2<short> ( -1, -1 ), const std::string &OutputDir = "" ) {
    if ( _outputDir.size ( ) == 0 && OutputDir.size ( ) == 0 )
      throw aol::Exception ( "No output directory specified!", __FILE__, __LINE__ );
    const std::string outputDir = ( _outputDir.size ( ) > 0 ) ? _outputDir : OutputDir;
    
    aol::Vec2<short> point ( Point );
    if ( point[0] < 0 || point[0] >= _data.getNumX ( ) || point[1] < 0 || point[1] >= _data.getNumY ( ) )
      point.set ( _mapper.splitGlobalIndex ( _data.getMaxIndexAndValue ( ).first )[0], _mapper.splitGlobalIndex ( _data.getMaxIndexAndValue ( ).first )[1] );
    translatePatternTo ( point );
    
    ColoredPictureType patternDataImg ( _data.getNumX ( ), _data.getNumY ( ) );
    for ( short i=0; i<3 ; ++i )
      patternDataImg[i] = _data;
    for ( short i=0; i<_pattern.size ( ) ; ++i ) {
      patternDataImg[0].set ( _pattern[i], _data.getMaxValue ( ) );
      patternDataImg[1].set ( _pattern[i], 0 );
      patternDataImg[2].set ( _pattern[i], 0 );
    }
    patternDataImg[0].set ( _pattern[0], 0 );
    patternDataImg[1].set ( _pattern[0], _data.getMaxValue ( ) );
    
    std::stringstream ss;
    ss << outputDir << "/patternDataImg.png";
    patternDataImg.savePNG ( ss.str ( ).c_str ( ) );
  }
  
private:
  void setFourierPowerCoefficients ( ) {
    qc::MultiArray<RealType, 2, 2> function ( _data.getNumX ( ), _data.getNumY ( ) ), transform ( _data.getNumX ( ), _data.getNumY ( ) );
    function[0] = _data;
    qc::FourierTransform ( function, transform );
    for ( short x=0; x<_data.getNumX ( ) ; ++x )
      for ( short y=0; y<_data.getNumY ( ) ; ++y )
        _fourierPowerCoefficients.set ( x, y, aol::Vec2<RealType> ( transform[0].get ( x, y ), transform[1].get ( x, y ) ).norm ( ) );
  }
  
  void setFourierPowerPeaks ( ) {
    eraseFourierPowerPeak ( 0, 0 );
    if ( _outputDir.size ( ) > 0 ) {
      _fourierPowerCoefficients.setOverflowHandlingToCurrentValueRange ( );
      _fourierPowerCoefficients.save ( aol::strprintf ( "%s/fourierPowerCoefficients_meanErased.pgm", _outputDir.c_str ( ) ).c_str ( ), qc::PGM_UNSIGNED_CHAR_BINARY );
    }
    
    for ( int i=0; i<6 ; ++i ) {
      const std::pair<int, RealType> maxIndexAndValue = _fourierPowerCoefficients.getMaxIndexAndValue ( );
      const aol::Vec2<short> maxIndex ( _mapper.splitGlobalIndex ( maxIndexAndValue.first )[0], _mapper.splitGlobalIndex ( maxIndexAndValue.first )[1] );
      _fourierPowerPeaks.set ( maxIndex, 1 );
      eraseFourierPowerPeak ( maxIndex );
    }
    if ( _outputDir.size ( ) > 0 ) {
      _fourierPowerPeaks.setOverflowHandlingToCurrentValueRange ( );
      _fourierPowerPeaks.save ( aol::strprintf ( "%s/fourierPowerPeaks.pgm", _outputDir.c_str ( ) ).c_str ( ), qc::PGM_UNSIGNED_CHAR_BINARY );
    }
  }
  
  void setPatternNormalForm ( ) {
    // 1. Search for two peaks that are closest to the origin
    MatrixType patternMatrix;
    PictureType fourierPowerPeakDistances ( _fourierPowerPeaks );
    RealType maxDist = aol::Vec2<short> ( _fourierPowerPeaks.getNumX ( ) + 1, _fourierPowerPeaks.getNumY ( ) + 1 ).norm ( );
    for ( int x=0; x<_fourierPowerPeaks.getNumX ( ) ; ++x ) {
      for ( int y=0; y<_fourierPowerPeaks.getNumY ( ) ; ++y ) {
        RealType val = _fourierPowerPeaks.get ( x, y );
        fourierPowerPeakDistances.set ( x, y, ( val > 0 ) ? aol::Vec2<short> ( x, y ).norm ( ) : maxDist );
      }
    }
    std::pair<int, RealType> minIndexAndValue = fourierPowerPeakDistances.getMinIndexAndValue ( );
    patternMatrix.set ( 0, 0, _mapper.splitGlobalIndex ( minIndexAndValue.first )[1] );
    patternMatrix.set ( 1, 0, _mapper.splitGlobalIndex ( minIndexAndValue.first )[0] );
    fourierPowerPeakDistances.set ( minIndexAndValue.first, maxDist );
    minIndexAndValue = fourierPowerPeakDistances.getMinIndexAndValue ( );
    patternMatrix.set ( 0, 1, _mapper.splitGlobalIndex ( minIndexAndValue.first )[1] );
    patternMatrix.set ( 1, 1, _mapper.splitGlobalIndex ( minIndexAndValue.first )[0] );
    _patternNormalForm = getPatternNormalForm ( patternMatrix );
  }
  
  void translatePatternTo ( const aol::Vec2<short> &Point ) {
    aol::Vec2<short> translation ( Point );
    translation -= _pattern[0];
    for ( int i=0; i<_pattern.size ( ) ; ++i ) {
      _pattern[i] += translation;
      _pattern[i][0] %= _data.getNumX ( ) - 1;
      _pattern[i][1] %= _data.getNumY ( ) - 1;
      _pattern[i][0] = ( _pattern[i][0] < 0 ) ? _pattern[i][0] + _data.getNumX ( ) - 1 : _pattern[i][0];
      _pattern[i][1] = ( _pattern[i][1] < 0 ) ? _pattern[i][1] + _data.getNumY ( ) - 1 : _pattern[i][1];
    }
  }
  
  void eraseFourierPowerPeak ( const short X, const short Y, const short RubberSize = 3 ) {
    const short offset = ( RubberSize - 1 ) / 2;
    for ( short dx=-offset; dx<=offset ; ++dx ) {
      for ( short dy=-offset; dy<=offset ; ++dy ) {
        if ( X+dx >= 0 && X+dx < _data.getNumX ( ) && Y+dy >= 0 && Y+dy < _data.getNumY ( ) )
          _fourierPowerCoefficients.set ( X+dx, Y+dy, 0 );
      }
    }
  }
  
  void eraseFourierPowerPeak ( const aol::Vec2<short> &Pos, const short RubberSize = 3 ) {
    eraseFourierPowerPeak ( Pos[0], Pos[1], RubberSize );
  }
  
  static void gcdOnRows ( MatrixType &M, const short Ri, const short Ci ) {
    if ( M.get ( Ri, Ci ) != 0 ) {
      // Modify by row addition, such that M[Ci,Ci] is non-negative
      if ( M.get ( Ci, Ci ) < 0 ) {
        for ( short col=0; col<2 ; ++col )
          M.set ( Ci, col, M.get ( Ci, col ) - floor ( M.get ( Ci, Ci ) / M.get ( Ri, Ci ) ) * M.get ( Ri, col ) );
      }
      // Make M[Ci,Ci] positive
      if ( M.get ( Ci, Ci ) == 0 ) {
        for ( short col=0; col<2 ; ++col )
          M.set ( Ci, col, M.get ( Ci, col ) + sign ( M.get ( Ri, Ci ) ) * M.get ( Ri, col ) );
      }
      // Make second entry in that column positive as well
      if ( M.get ( Ri, Ci ) < 0 ) {
        int f = ceil ( M.get ( Ri, Ci ) / M.get ( Ci, Ci ) );
        if ( f == 0 )
          ++f;
        for ( short col=0; col<2 ; ++col )
          M.set ( Ri, col, M.get ( Ri, col ) - f * M.get ( Ci, col ) );
      }
      // Euclidian algorithm on rows in order to get M[Ri,Ci] to zero
      while ( M.get ( Ri, Ci ) != 0 ) {
        if ( abs ( M.get ( Ci, Ci ) ) > abs ( M.get ( Ri, Ci ) ) ) {
          int f = floor ( M.get ( Ci, Ci ) / M.get ( Ri, Ci ) );
          if ( M.get ( Ci, Ci ) % M.get ( Ri, Ci ) == 0 )
            f = f - sign ( M.get ( Ri, Ci ) ) * sign ( M.get ( Ci, Ci ) );
          for ( short col=0; col<2 ; ++col )
            M.set ( Ci, col, M.get ( Ci, col ) - f * M.get ( Ri, col ) );
        } else {
          int f = floor ( M.get ( Ri, Ci ) / M.get ( Ci, Ci ) );
          for ( short col=0; col<2 ; ++col )
            M.set ( Ri, col, M.get ( Ri, col ) - f * M.get ( Ci, col ) );
        }
      }
    }
  }
  
  void setPattern ( ) {
    setPattern ( _pattern, _patternNormalForm, aol::Vec2<short> ( _data.getNumX ( ) - 1, _data.getNumY (  ) - 1 ) );
  }
  
  static int sign (int val) {
    return ( val < 0 ) ? -1 : ( ( val > 0 ) ? 1 : 0 );
  }
};


template <typename _RealType, typename _PictureType>
class BravaisLatticeImageCreator {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const std::string &_outputDir;
public:
  BravaisLatticeImageCreator ( const std::string &OutputDir = "" ) : _outputDir ( OutputDir ) { }
  
  void getBravaisLatticeImage ( const aol::MultiVector<RealType> &LatticeVectors, const aol::MultiVector<RealType> &GaussianParameters, PictureType &Result ) {
    aol::MultiVector<RealType> gaussianParams ( GaussianParameters );
    gaussianParams.resize ( gaussianParams.numComponents ( ), im::AsymmetricGaussianBumpFunction<RealType>::NumberOfParameters );
    PictureType unitCellImg ( Result.getNumX ( ), Result.getNumY ( ) );
    im::AtomFinder<RealType>::getBumpFunctionImage ( gaussianParams, unitCellImg );
    
    // Determine reasonable boundaries for unit cell image
    aol::Vector<RealType> heights ( gaussianParams.numComponents ( ) );
    gaussianParams.getTo ( 3, heights );
    RealType threshold = 1e-4 * heights.getMaxValue ( );
    PictureType block ( 1, Result.getNumY ( ) );
    int nx = unitCellImg.getNumX ( ) - 1;
    RealType maxVal = 0;
    while ( maxVal < threshold && nx > 1 ) {
      unitCellImg.copyBlockTo ( nx, 0, block );
      maxVal = block.getMaxValue ( );
      --nx;
    }
    nx += 2;
    block.reallocate ( Result.getNumX ( ), 1 );
    int ny = unitCellImg.getNumY ( ) - 1;
    maxVal = 0;
    while ( maxVal < threshold && ny > 1 ) {
      unitCellImg.copyBlockTo ( 0, ny, block );
      maxVal = block.getMaxValue ( );
      --ny;
    }
    ny += 2;
    unitCellImg.crop ( aol::Vec2<int> ( 0, 0 ), aol::Vec2<int> ( nx, ny ) );
    
    if ( _outputDir != "" )
      unitCellImg.save ( aol::strprintf ( "%s/unitCellImg%s", _outputDir.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    
    aol::MultiVector<RealType> latticePoints;
    getBravaisLattice ( LatticeVectors, aol::Vec2<RealType> ( 0, 0 ), aol::Vec2<RealType> ( Result.getNumX ( ), Result.getNumY ( ) ), latticePoints );
    getBravaisLatticeImage ( latticePoints, unitCellImg, Result );
  }
  
  static void getBravaisLattice ( const aol::MultiVector<RealType> &LatticeVectors,
                                  const aol::Vec2<RealType> &X0, const aol::Vec2<RealType> &XEnd,
                                  aol::MultiVector<RealType> &LatticePoints,
                                  const aol::Vec2<RealType> &Origin = aol::Vec2<RealType> ( 0, 0 ) ) {
    LatticePoints.resize ( 0, 2 );
    RealType diam = ( XEnd[0] - X0[0] + 1 ) + ( XEnd[1] - X0[1] + 1 );
    int k1Bound = diam / LatticeVectors[0].norm ( ), k2Bound = diam / LatticeVectors[1].norm ( );
    aol::Vector<RealType> pos ( 2 );
    for ( int k1=-k1Bound; k1<=k1Bound ; ++k1 ) {
      for ( int k2=-k2Bound; k2<=k2Bound ; ++k2 ) {
        pos[0] = Origin[0];
        pos[1] = Origin[1];
        pos.addMultiple ( LatticeVectors[0], k1 );
        pos.addMultiple ( LatticeVectors[1], k2 );
        if ( aol::InsideQuad<qc::QC_2D, RealType> ( X0, XEnd, aol::Vec2<RealType> ( pos[0], pos[1] ) ) ) {
          LatticePoints.resize ( LatticePoints.numComponents ( ) + 1, 2 );
          LatticePoints[LatticePoints.numComponents ( ) - 1][0] = pos[0];
          LatticePoints[LatticePoints.numComponents ( ) - 1][1] = pos[1];
        }
      }
    }
  }
  
  static void getBravaisLatticeImage ( const aol::MultiVector<RealType> &LatticePoints, const PictureType &UnitCell,
                                       PictureType &Result ) {
    for ( int i=0; i<LatticePoints.numComponents ( ) ; ++i ) {
      for ( int dy=0; dy<UnitCell.getNumY ( ) ; ++dy ) {
        for ( int dx=0; dx<UnitCell.getNumX ( ) ; ++dx ) {
          int x = round ( LatticePoints[i][0] ) + dx, y = round ( LatticePoints[i][1] ) + dy;
          if ( x >= 0 && x < Result.getNumX ( ) && y >=0 && y < Result.getNumY ( ) )
            Result.add ( x, y, UnitCell.interpolate ( x - LatticePoints[i][0], y - LatticePoints[i][1] ) );
        }
      }
    }
  }
  
  static void getBravaisLatticeImage ( const aol::MultiVector<RealType> &LatticePoints,
                                       const RealType AtomDiameter, const RealType AtomBrightness,
                                       PictureType &Result ) {
    for ( int i=0; i<LatticePoints.numComponents ( ) ; ++i )
      drawGaussianBell ( Result, aol::Vec2<RealType> ( LatticePoints[i][0], LatticePoints[i][1] ), AtomDiameter, AtomBrightness );
  }
  
  static void removeBoundaryAtoms ( const aol::Vec2<RealType> &X0, const aol::Vec2<RealType> &XEnd,
                                    const RealType Padding,
                                    aol::MultiVector<RealType> &LatticePoints ) {
    aol::MultiVector<RealType> latticePoints ( LatticePoints );
    LatticePoints.resize ( 0, 2 );
    for ( int i=0; i<latticePoints.numComponents ( ) ; ++i ) {
      if ( latticePoints[i][0] - Padding >= X0[0] && latticePoints[i][0] + Padding < XEnd[0] && latticePoints[i][1] - Padding >= X0[1] && latticePoints[i][1] + Padding < XEnd[1] ) {
        LatticePoints.resize ( LatticePoints.numComponents ( ) + 1, 2 );
        LatticePoints[LatticePoints.numComponents ( ) - 1] = latticePoints[i];
      }
    }
  }
  
  static void getPiecewiseBravaisLattice ( const aol::MultiVector<RealType> &BaseLatticeVectors,
                                           const qc::ScalarArray<int, qc::QC_2D> &Mask,
                                           aol::RandomAccessContainer<aol::MultiVector<RealType> > &LatticeVectors,
                                           aol::RandomAccessContainer<aol::MultiVector<RealType> > &LatticePoints,
                                           const aol::RandomAccessContainer<aol::Vec2<RealType> > &Origins = aol::RandomAccessContainer<aol::Vec2<RealType> > ( ) ) {
    LatticeVectors.reallocate ( Mask.getMaxValue ( ) + 1 );
    aol::RandomGenerator rng;
    for ( int l=0; l<= Mask.getMaxValue ( ) ; ++l ) {
      LatticeVectors[l].resize ( 2, 2 );
      aol::MultiVector<RealType> centers;
      const RealType angle = rng.rReal<RealType> ( ) * 2.0 * aol::NumberTrait<RealType>::pi;
      for ( int i=0; i<2 ; ++i ) {
        LatticeVectors[l][i][0] = cos ( angle ) * BaseLatticeVectors[i][0] + sin ( angle ) * BaseLatticeVectors[i][1];
        LatticeVectors[l][i][1] = -sin ( angle ) * BaseLatticeVectors[i][0] + cos ( angle ) * BaseLatticeVectors[i][1];
      }
      getBravaisLattice ( LatticeVectors[l], aol::Vec2<RealType> ( 0, 0 ), aol::Vec2<RealType> ( Mask.getNumX ( ), Mask.getNumY ( ) ),
                                 centers, ( l < Origins.size ( ) ) ? Origins[l] : aol::Vec2<RealType> ( 0, 0 ) );
      LatticePoints.pushBack ( centers );
    }
    removeBoundaryAtoms ( Mask, 0.5, LatticePoints );
  }
  
  static void getPiecewiseBravaisLatticeImage ( const aol::RandomAccessContainer<aol::MultiVector<RealType> > &LatticePoints,
                                                const RealType AtomDiameter, const RealType AtomBrightness,
                                                PictureType &Result ) {
    Result.setZero ( );
    for ( int l=0; l<LatticePoints.size ( ) ; ++l ) {
      PictureType res ( Result, aol::STRUCT_COPY );
      getBravaisLatticeImage ( LatticePoints[l], AtomDiameter, AtomBrightness, res );
      Result += res;
    }
  }
  
  static void getMosaikMask ( const int NumSegments, qc::ScalarArray<int, qc::QC_2D> &Mask ) {
    aol::RandomGenerator rng;
//    rng.randomize ( );    turn this on to enable true randomness
    aol::Vector<int> min ( 2 ), max ( 2 );
    max[0] = Mask.getNumX ( );
    max[1] = Mask.getNumY ( );
    aol::MultiVector<int> regionCentroids ( 2, NumSegments );
    rng.rIntMultiVecPairwiseDifferent ( regionCentroids, min, max );
    for ( int y=0; y<Mask.getNumY ( ) ; ++y ) {
      for ( int x=0; x<Mask.getNumX ( ) ; ++x ) {
        aol::Vector<RealType> distances ( NumSegments );
        for ( int l=0; l<NumSegments ; ++l ) distances[l] = aol::Vec2<RealType> ( x - regionCentroids[0][l], y - regionCentroids[1][l] ).norm ( );
        Mask.set ( x, y, distances.getMinIndexAndValue ( ).first );
      }
    }
  }
  
  static void removeBoundaryAtoms ( const qc::ScalarArray<int, qc::QC_2D> &Mask,
                                    const RealType Padding,
                                    aol::RandomAccessContainer<aol::MultiVector<RealType> > &LatticePoints ) {
    for ( int l=0; l<LatticePoints.size ( ) ; ++l ) {
      removeBoundaryAtoms ( aol::Vec2<RealType> ( 0, 0 ), aol::Vec2<RealType> ( Mask.getNumX ( ), Mask.getNumY (  ) ), Padding, LatticePoints[l] );
      aol::MultiVector<RealType> latticePoints ( LatticePoints[l] );
      LatticePoints[l].resize ( 0, 2 );
      for ( int i=0; i<latticePoints.numComponents ( ) ; ++i ) {
        bool latticePointIsInterior = true;
        for ( int y=0; y<Mask.getNumY ( ) && latticePointIsInterior ; ++y ) {
          for ( int x=0; x<Mask.getNumX ( ) && latticePointIsInterior ; ++x ) {
            if ( Mask.get ( x, y ) != l && aol::Abs<RealType> ( latticePoints[i][0] - x ) < Padding && aol::Abs<RealType> ( latticePoints[i][1] - y ) < Padding )
              latticePointIsInterior = false;
          }
        }
        if ( latticePointIsInterior ) {
          LatticePoints[l].resize ( LatticePoints[l].numComponents ( ) + 1, 2 );
          LatticePoints[l][LatticePoints[l].numComponents ( ) - 1] = latticePoints[i];
        }
      }
    }
  }
  
  static void removeBoundaryAtoms ( const qc::ScalarArray<int, qc::QC_2D> &Mask,
                                    const RealType Padding,
                                    aol::MultiVector<RealType> &LatticePoints ) {
    removeBoundaryAtoms ( aol::Vec2<RealType> ( 0, 0 ), aol::Vec2<RealType> ( Mask.getNumX ( ), Mask.getNumY (  ) ), Padding, LatticePoints );
    
    aol::MultiVector<RealType> boundaryPoints ( 0, 2 );
    for ( int y=1; y<Mask.getNumY ( )-1 ; ++y ) {
      for ( int x=1; x<Mask.getNumX ( )-1 ; ++x ) {
        bool posIsBoundaryPoint = false;
        for ( int dy=-1; dy<=1 && !posIsBoundaryPoint ; ++dy )
          for ( int dx=-1; dx<=1 && !posIsBoundaryPoint ; ++dx )
            if ( Mask.get ( x, y ) != Mask.get ( x + dx, y + dy ) ) posIsBoundaryPoint = true;
        if ( posIsBoundaryPoint ) {
          boundaryPoints.resize ( boundaryPoints.numComponents ( ) + 1, 2 );
          boundaryPoints[boundaryPoints.numComponents ( ) - 1][0] = x;
          boundaryPoints[boundaryPoints.numComponents ( ) - 1][1] = y;
        }
      }
    }
    
    aol::MultiVector<RealType> latticePoints ( LatticePoints );
    LatticePoints.resize ( 0, 2 );
    for ( int i=0; i<latticePoints.numComponents ( ) ; ++i ) {
      bool latticePointIsInterior = true;
      for ( int j=0; j<boundaryPoints.numComponents ( ) && latticePointIsInterior ; ++j ) {
        if ( aol::Abs<RealType> ( latticePoints[i][0] - boundaryPoints[j][0] ) < Padding && aol::Abs<RealType> ( latticePoints[i][1] - boundaryPoints[j][1] ) < Padding )
          latticePointIsInterior = false;
      }
      if ( latticePointIsInterior ) {
        LatticePoints.resize ( LatticePoints.numComponents ( ) + 1, 2 );
        LatticePoints[LatticePoints.numComponents ( ) - 1] = latticePoints[i];
      }
    }
  }
  
  static void getLatticePointsUnion ( const aol::RandomAccessContainer<aol::MultiVector<RealType> > &LatticePoints,
                                      aol::MultiVector<RealType> &LatticePointsUnion ) {
    LatticePointsUnion.resize ( 0, 2 );
    int k=0;
    for ( int l=0; l<LatticePoints.size ( ) ; ++l ) {
      LatticePointsUnion.resize ( LatticePointsUnion.numComponents ( ) + LatticePoints[l].numComponents ( ), 2 );
      for ( int i=0; i<LatticePoints[l].numComponents ( ) ; ++i, ++k ) LatticePointsUnion[k] = LatticePoints[l][i];
    }
  }
  
  static void getHexBravaisLattice ( const RealType AtomSeparation, const RealType Angle, const aol::Vec2<RealType> &X0, const aol::Vec2<RealType> &XEnd, aol::MultiVector<RealType> &LatticePoints ) {
    aol::MultiVector<RealType> hexLatticeVectors ( 2, 2 );
    hexLatticeVectors[0][0] = 1;
    hexLatticeVectors[0][1] = 0;
    hexLatticeVectors[1][0] = cos ( 2.0 * aol::NumberTrait<RealType>::pi / 3.0 );
    hexLatticeVectors[1][1] = sin ( 2.0 * aol::NumberTrait<RealType>::pi / 3.0 );
    hexLatticeVectors *= AtomSeparation;
    
    aol::MultiVector<RealType> rotatedHexLatticeVectors ( 2, 2 );
    for ( int i=0; i<2 ; ++i ) {
      rotatedHexLatticeVectors[i][0] = cos ( Angle ) * hexLatticeVectors[i][0] + sin ( Angle ) * hexLatticeVectors[i][1];
      rotatedHexLatticeVectors[i][1] = -sin ( Angle ) * hexLatticeVectors[i][0] + cos ( Angle ) * hexLatticeVectors[i][1];
    }
    getBravaisLattice ( rotatedHexLatticeVectors, X0, XEnd, LatticePoints );
  }
  
  static void getPiecewiseHexBravaisLattice ( const int NumSegments, const RealType AtomSeparation,
                                              const aol::Vec2<RealType> &X0, const aol::Vec2<RealType> &XEnd,
                                              aol::RandomAccessContainer<aol::MultiVector<RealType > > &LatticePoints ) {
    LatticePoints.reallocate ( NumSegments );
    aol::RandomGenerator rng;
    for ( int l=0; l<NumSegments ; ++l ) getHexBravaisLattice ( AtomSeparation, rng.rReal<RealType> ( ) * 2.0 * aol::NumberTrait<RealType>::pi, X0, XEnd, LatticePoints[l] );
  }
protected:
  static void drawGaussianBell ( PictureType &Picture, const aol::Vec2<RealType> &Pos, const RealType AtomDiameter, const RealType AtomBrightness ) {
    im::AsymmetricGaussianBumpFunction<RealType> bumpFunction ( Pos, AtomBrightness, AtomDiameter, AtomDiameter, 0.0, 0.0 );
    for ( int y=aol::Max<int> ( 0, floor ( Pos[1] - 3 * AtomDiameter ) ); y<aol::Min<int> ( Picture.getNumY ( ), ceil ( Pos[1] + 3 * AtomDiameter ) ) ; ++y )
      for ( int x=aol::Max<int> ( 0, floor ( Pos[0] - 3 * AtomDiameter ) ); x<aol::Min<int> ( Picture.getNumX ( ), ceil ( Pos[0] + 3 * AtomDiameter ) ) ; ++x )
        Picture.add ( x, y, bumpFunction.evaluate ( aol::Vec2<RealType> ( x, y ) ) );
  }
};
  
} // namespace im

#endif
