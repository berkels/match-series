#ifndef __GENERATOR_H
#define __GENERATOR_H

#include <indexMapper.h>
#include <multiArray.h>
#include <ChanVese.h>
#include <hashSet.h>
#include <configurators.h>

namespace qc {

/**
 * \author Berkels
 */
template <typename RealType>
class ArtificialPFCFrontFunction{
  const RealType _domainScaleFactor;
  const RealType _frontPosition;
  const RealType _transitionWidth;
  const bool _circularFront;
public:
  ArtificialPFCFrontFunction ( const RealType DomainScaleFactor,
                               const RealType FrontPosition,
                               const RealType TransitionWidth,
                               const bool CircularFront )
    : _domainScaleFactor ( DomainScaleFactor ),
      _frontPosition ( FrontPosition ),
      _transitionWidth ( TransitionWidth ),
      _circularFront ( CircularFront ) { }

  RealType evaluate( const aol::Vec2<RealType> &X ) const{
    const RealType q = aol::NumberTrait<long double>::pi * sqrt ( 3. );
    const RealType scaledX = _domainScaleFactor * X[0];
    const RealType scaledY = _domainScaleFactor * X[1];
    const RealType periodic = cos(q*scaledX) * cos(q*scaledY/sqrt(3.)) - 0.5 * cos(2*q*scaledY/sqrt(3.));
    const RealType temp = _circularFront ? X.norm() : X[0] ;
    const RealType cut = 0.5*(1.-tanh((temp-_frontPosition)/_transitionWidth));
    return cut*(periodic+1) + 0.2;
  }
};

/**
 * \brief A class for generating some standard deformations and levelset functions.
 *
 * Note: Most of the functions only work properly in 2D, but still compile in 3D.
 * If you use DataGenerator in 3D, verify that the functions you use are working as
 * intended.
 * \see ShapeLevelsetGenerator
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class DataGenerator {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::InitType InitType;
  typedef typename ConfiguratorType::VectorType VectorType;
  const ConfiguratorType _config;
  const InitType &_grid;
  typename ConfiguratorType::InitType::OldAllNodeIterator _fnit;
  const RealType _h;
public:
  DataGenerator ( const InitType &Grid )
    : _config ( Grid ),
      _grid ( Grid ),
      _h ( static_cast<RealType> ( Grid.H() ) ) {};
  void generateDeformation ( const RealType Alpha, aol::MultiVector<RealType> &Deformation ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      Deformation[0][ _config.localToGlobal ( *_fnit, 0 )  ] = Alpha * ( ( *_fnit ) [0] * _h - 0.5 ) * fabs ( ( *_fnit ) [1] * _h - 0.5 );
    }
  }
  void generateSkewSymmetricDeformation ( const RealType Alpha, aol::MultiVector<RealType> &Deformation ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      Deformation[0][ _config.localToGlobal ( *_fnit, 0 )  ] = -1. * Alpha * ( ( *_fnit ) [1] * _h - 0.5 );
      Deformation[1][ _config.localToGlobal ( *_fnit, 0 )  ] = Alpha * ( ( *_fnit ) [0] * _h - 0.5 );
    }
  }
  void generateSymmetricDeformation ( const RealType Alpha,
                                      aol::MultiVector<RealType> &Deformation,
                                      const typename ConfiguratorType::VecType &Center ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        Deformation[i][ _config.localToGlobal ( *_fnit, 0 )  ] = Alpha * ( ( *_fnit ) [i] * _h - Center[i] );
    }
  }
  void generateSymmetricDeformation ( const RealType Alpha,
                                      aol::MultiVector<RealType> &Deformation,
                                      const RealType CenterX = 0.5,
                                      const RealType CenterY = 0.5 ) {
    aol::Vec2<RealType> center ( CenterX, CenterY );
    generateSymmetricDeformation ( Alpha, Deformation, center );
  }
  void generateShearDeformation ( const RealType Alpha,
                                  aol::MultiVector<RealType> &Deformation,
                                  const RealType CenterX = 0.5 ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      Deformation[0][ _config.localToGlobal ( *_fnit, 0 )  ] = Alpha * ( ( *_fnit ) [1] * _h - CenterX );
    }
  }
  void generateNonLinearDeformation ( const RealType Alpha, aol::MultiVector<RealType> &Deformation ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const RealType x = ( *_fnit ) [0] * _h;
      const RealType y = ( *_fnit ) [1] * _h;
      Deformation[0][ _config.localToGlobal ( *_fnit, 0 )  ] = x * ( -1. * Alpha * ( y - 0.5 ) ) + ( 1 - x ) * ( Alpha * ( x - 0.5 ) );
      Deformation[1][ _config.localToGlobal ( *_fnit, 0 )  ] = x * ( Alpha * ( x - 0.5 ) ) + ( 1 - x ) * ( Alpha * ( y - 0.5 ) );
    }
  }
  void generateStretchDeformation ( const RealType StretchStart,
                                    const RealType StretchEnd,
                                    const RealType StrechFactor,
                                    aol::MultiVector<RealType> &Deformation ) {
    Deformation[1].setZero();
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const RealType x = ( *_fnit ) [0] * _h;
      const RealType displacement = ( 1 / StrechFactor - 1 ) * aol::Clamp<RealType> ( ( x - StretchStart ), 0, ( StretchEnd - StretchStart ) );

      Deformation[0][ _config.localToGlobal ( *_fnit, 0 )  ] = displacement;
    }
  }
  void generateIdentity ( aol::MultiVector<RealType> &Identity ) {
    typename ConfiguratorType::VecType zeroVec;
    generateSymmetricDeformation ( 1., Identity, zeroVec );
  }
  template <typename ParametricDeformationType, bool ClipCoord>
  void generateDeformationFromParametricDeformation ( const ParametricDeformationType &ParDef,
                                                      aol::MultiVector<RealType> &Deformation ) {
    typename ConfiguratorType::VecType transformedLocalCoord, zeroVec;
    qc::Element transformedEl;

    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      ParDef.template evaluateDeformation<ClipCoord> ( *_fnit, zeroVec, transformedEl, transformedLocalCoord );
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        Deformation[i][ _config.localToGlobal ( *_fnit, 0 ) ] = ( transformedEl[i] + transformedLocalCoord[i] - ( *_fnit ) [i] ) * _h;
    }
  }
  void generateLineLevelset ( VectorType &LevelsetFunction, const int Axis = 0, const RealType Offset = 0.5 ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 )  ] = 1 * ( ( *_fnit ) [Axis] * _h - Offset );

    }
  }
  void generatePlaneLevelset ( VectorType &LevelsetFunction, const typename ConfiguratorType::VecType &Normal, const typename ConfiguratorType::VecType &Point ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      RealType temp = 0;
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        temp += ( ( *_fnit ) [i] *_h - Point[i] ) * Normal[i];
      LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 )  ] = temp;
    }
  }
  void generateDiagonalLevelset ( VectorType &LevelsetFunction, const int DiagonalNumber = 0 ) {
    aol::Vec2<RealType> lineNormal (1.,1);
    if ( DiagonalNumber == 1 )
      lineNormal[1] = -1.;
    lineNormal /= lineNormal.norm();
    aol::Vec2<RealType> point ( 0.5, 0.5 );
    generatePlaneLevelset ( LevelsetFunction, lineNormal, point );
  }

  /**
   * \author Olischlaeger
   */
  void generateDiagonalLevelsetPlusCirclesInTheCorners ( VectorType &LevelsetFunction ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      RealType x = ( *_fnit ) [0] *_h;
      RealType y = ( *_fnit ) [1] *_h;

      RealType absToDiagonal = sqrt( aol::Sqr(0.5*(1-y-x)) +aol::Sqr(0.5*(x-1+y)) );
      RealType skp = x-1+y;
      RealType absToOrigin = sqrt( aol::Sqr(x) +aol::Sqr(y) );
      RealType absToOne = sqrt( aol::Sqr(x-1) +aol::Sqr(y-1) );

      if (skp>0){
        if(absToOne<0.3){
          LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 ) ] = -(0.3-absToOne); //-1
        }
        else{
          if( (absToOne-0.3)<absToDiagonal  ){
            LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 ) ] =absToOne-0.3;  // 1
          }
          else{
            LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 ) ] = absToDiagonal;
          }
        }
      }
      else{
        if(absToOrigin<0.3){
          LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 ) ] = 0.3-absToOrigin; //1
        }
        else{
          if( (absToOrigin-0.3)<absToDiagonal  ){
            LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 ) ] = -(absToOrigin-0.3);  // -1
          }
          else{
            LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 ) ] = -absToDiagonal;
          }
        }
      }
    }
  }

  /**
   * \author Olischlaeger
   */
  void generateThreeDiagonalsPiecewiseConstant ( VectorType &LevelsetFunction ) {
    VectorType Initialization(LevelsetFunction.size());
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      RealType x = ( *_fnit ) [0] *_h;
      RealType y = ( *_fnit ) [1] *_h;

      //RealType absToDiagonal = sqrt( aol::Sqr(0.5*(1-y-x)) +aol::Sqr(0.5*(x-1+y)) );
      RealType skp = x-1+y;
      RealType absToOrigin = sqrt( aol::Sqr(x) +aol::Sqr(y) );
      RealType absToOne = sqrt( aol::Sqr(x-1) +aol::Sqr(y-1) );

      if (skp>0){
        if(absToOne<0.3){
          LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 )  ] = -1;
        }
        else{
          LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 )  ] = 1;
        }
      }
      else{
        if(absToOrigin<0.3){
          LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 )  ] = 1;
        }
        else{
          LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 )  ] = -1;
        }
      }
    }
  }
  void generateDoubleLineLevelset ( VectorType &LevelsetFunction, const int Axis = 0, const RealType Center = 0.5, const RealType Offset = 0.25 ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 )  ] = 1 * fabs( ( *_fnit ) [Axis] * _h - Center ) - Offset;
    }
  }
  void generateEllipsoidLevelset ( VectorType &LevelsetFunction, const typename ConfiguratorType::VecType &Center, const typename ConfiguratorType::VecType &Scaling, const RealType Offset = 0.5 ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      RealType tmp = 0;
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        tmp += aol::Sqr ( Scaling[i]*(( *_fnit ) [i] * _h - Center[i]) );
      LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 )  ] = sqrt ( tmp ) - Offset;
    }
  }
  void generateSphereLevelset ( VectorType &LevelsetFunction, const typename ConfiguratorType::VecType &Center, const RealType Offset = 0.5 ) {
    typename ConfiguratorType::VecType scaling;
    scaling.setAll( aol::NumberTrait<RealType>::one );
    generateEllipsoidLevelset ( LevelsetFunction, Center, scaling, Offset );
  }
  void generateEllipseLevelset ( VectorType &LevelsetFunction, const RealType Offset = 0.5, const RealType CenterX = 0.5, const RealType CenterY = 0.5, const RealType ScaleX = 1.,  const RealType ScaleY = 1. ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 )  ] = sqrt ( ( aol::Sqr ( ScaleX*(( *_fnit ) [0] * _h - CenterX) ) + aol::Sqr ( ScaleY*(( *_fnit ) [1] * _h - CenterY) ) ) ) - Offset;
    }
  }
  void generateL1SphereLevelset ( VectorType &LevelsetFunction, const typename ConfiguratorType::VecType &Center, const RealType Offset = 0.5 ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      RealType tmp = 0;
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        tmp += aol::Abs (( *_fnit ) [i] * _h - Center[i]);
      LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 )  ] = tmp - Offset;
    }
  }
  void generateRegularizedL1SphereLevelset ( VectorType &LevelsetFunction, const typename ConfiguratorType::VecType &Center, const RealType Offset = 0.5, const RealType Epsilon = 0. ) {
    const RealType epsilonSqr = aol::Sqr ( Epsilon );
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      RealType tmp = 0;
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        tmp += sqrt ( aol::Sqr ( aol::Abs (( *_fnit ) [i] * _h - Center[i]) ) + epsilonSqr );
      LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 )  ] = tmp - Offset;
    }
  }
  void generateScaledLInfinitySphereLevelset ( VectorType &LevelsetFunction, const typename ConfiguratorType::VecType &Center, const typename ConfiguratorType::VecType &Scaling, const RealType Offset = 0.5 ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      RealType tmp = 0;
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        tmp = aol::Max ( Scaling[i] * aol::Abs (( *_fnit ) [i] * _h - Center[i]), tmp );
      LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 )  ] = tmp - Offset;
    }
  }
  void generateLInfinitySphereLevelset ( VectorType &LevelsetFunction, const typename ConfiguratorType::VecType &Center, const RealType Offset = 0.5 ) {
    typename ConfiguratorType::VecType scaling;
    scaling.setAll ( 1. );
    generateScaledLInfinitySphereLevelset ( LevelsetFunction, Center, scaling, Offset );
  }
  //! Generates levelset function whose zero-levelset is a cylinder aligned to the Dim-axis, centered around Center,
  //! whose cross-section is cylindrical and which extends from the center in direction i by Offset/Scaling[i]
  void generateCylindricalLevelset ( VectorType &LevelsetFunction, const typename ConfiguratorType::VecType &Center, const typename ConfiguratorType::VecType &Scaling, const RealType Offset = 0.5 ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      RealType tmp = 0;
      for ( int i = 0; i < ConfiguratorType::Dim - 1; i++ )
        tmp += aol::Sqr ( Scaling[i]*(( *_fnit ) [i] * _h - Center[i]) );
      tmp = aol::Max ( Scaling[ConfiguratorType::Dim - 1] * aol::Abs (( *_fnit ) [ConfiguratorType::Dim - 1] * _h - Center[ConfiguratorType::Dim - 1]), sqrt ( tmp ) );
      LevelsetFunction[ _config.localToGlobal ( *_fnit, 0 )  ] = tmp - Offset;
    }
  }
  void generateCircleLevelset ( VectorType &LevelsetFunction, const RealType Offset = 0.5, const RealType CenterX = 0.5, const RealType CenterY = 0.5 ) {
    generateEllipseLevelset( LevelsetFunction, Offset, CenterX, CenterY, 1., 1. );
  }
  void addLevelset ( const VectorType &LevelsetFunctionArg, VectorType &LevelsetFunctionDest ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const int index = _config.localToGlobal ( *_fnit, 0 ) ;
      if ( LevelsetFunctionArg[index] < LevelsetFunctionDest[index] )
        LevelsetFunctionDest[index] = LevelsetFunctionArg[index];
    }
  }
  void subtractLevelset ( const VectorType &LevelsetFunctionArg, VectorType &LevelsetFunctionDest ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const int index = _config.localToGlobal ( *_fnit, 0 ) ;
      if ( -LevelsetFunctionArg[index] > LevelsetFunctionDest[index] )
        LevelsetFunctionDest[index] = -LevelsetFunctionArg[index];
    }
  }
  void generateCircleLevelsets ( VectorType &LevelsetFunction, const int Number ) {
    const RealType diameter = 0.4/Number;
    const RealType offset = 0.5/Number;
    generateCircleLevelset ( LevelsetFunction, diameter, offset, offset );
    VectorType tempCircle( LevelsetFunction, aol::STRUCT_COPY );
    for ( int i = 0; i < Number; i++ ){
      for ( int j = 0; j < Number; j++ ){
        if ( i != 0 || j != 0 ) {
          generateCircleLevelset ( tempCircle, diameter, offset*(2*i+1), offset*(2*j+1) );
          addLevelset( tempCircle, LevelsetFunction );
        }
      }
    }
  }
  void generatePeriodicSphereLevelset( VectorType &Levelset, const RealType Radius = 0.5 ) {
    typename ConfiguratorType::ArrayType  levelset( Levelset, _grid, aol::FLAT_COPY ), aux( levelset, aol::STRUCT_COPY );
    typename ConfiguratorType::VecType center;
    center.setAll( .5 );
    int Offset[3] = { levelset.getNumX() >> 1, levelset.getNumY() >> 1, levelset.getNumZ() >> 1 };
    generateSphereLevelset( aux, center, Radius );
    for ( int z = 0; z < levelset.getNumZ(); z++ )
      for ( int y = 0; y < levelset.getNumY(); y++ )
        for ( int x = 0; x < levelset.getNumX(); x++ )
          ( qc::Array<RealType>( levelset, aol::FLAT_COPY ) ).set( x, y, z, ( qc::Array<RealType>( aux, aol::FLAT_COPY ) ).getPeriodic( x + Offset[0], y + Offset[1], z + Offset[2] ) );
  }
  void generatePeriodicSphereCharacFunc( VectorType &CharacFunc, const RealType Radius = 0.5 ) {
    generatePeriodicSphereLevelset( CharacFunc, Radius );
    CharacFunc *= -1. / _h;
    CharacFunc.clamp( -.5, .5 );
    CharacFunc.addToAll( .5 );
  }
  void generateLineLevelsets ( VectorType &LevelsetFunction, const int Axis = 0, const int Depth = 1 ) {
    const RealType offset = 0.5/Depth;
    generateDoubleLineLevelset ( LevelsetFunction, Axis, offset, offset*0.5 );
    VectorType tempCircle( LevelsetFunction, aol::STRUCT_COPY );
    for ( int i = 1; i < Depth; i++ ){
      generateDoubleLineLevelset ( tempCircle, Axis, offset*(2*i+1), offset*0.5 );
      addLevelset( tempCircle, LevelsetFunction );
    }
  }
  void generateChanVeseInitialization ( aol::MultiVector<RealType> &LevelsetFunction ) {
    for( int i = 0; i < LevelsetFunction.numComponents(); i++ )
    {
      if ( i < 2 ){
        generateLineLevelset( LevelsetFunction[i], i%2, 0.5 * ( ( i == 1 ) ? ( _grid.getNumY() - 1 ) * _grid.H() : 1. ));
      }
      else{
        generateLineLevelsets( LevelsetFunction[i], i%2, i/2 );
      }
    }
  }
  void generatePiecewiseAffineFlowField ( aol::MultiVector<RealType> &FlowField,
                                          const aol::MultiVector<RealType> &LevelsetFunctions,
                                          const aol::MultiVector<RealType> &DeformationParameters,
                                          const RealType Isovalue = aol::NumberTrait<RealType>::zero ) {
    aol::Vector<int> intVec( LevelsetFunctions.numComponents() );
    const int dim = FlowField.numComponents();
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const int globalIndex = _config.localToGlobal ( *_fnit, 0 ) ;
      const int segmentIndex = aol::PowerSetIterator::getPositionNumberFromLevelsetFunctions<RealType>(LevelsetFunctions, globalIndex, Isovalue);
      for ( int i = 0; i < dim; i++ ){
        RealType temp = 0.;
        for ( int j = 0; j < dim; j++ ){
          temp += DeformationParameters[segmentIndex][dim*i+j]*( *_fnit )[j] * _h;
        }
        FlowField[i][globalIndex] = temp + DeformationParameters[segmentIndex][aol::Sqr(dim)+i];
      }
    }
  }
  /**
   * Discretizes the affine transformation \f$ Ax+b \f$ of a position x with matrix \f$ A \f$ and vector \f$ b \f$ passed as arguments.
   * The corresponding discretized displacement is put into "Displacement".
   *
   * \author Wirth
   */
  void generateAffineDisplacement ( const typename ConfiguratorType::MatType DeformMatrix, const typename ConfiguratorType::VecType Shift, aol::MultiVector<RealType> &Displacement ) {
    // for all grid nodes
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      // compute the displacement
      typename ConfiguratorType::VecType position, zeroCoords;
      _config.getGlobalCoords( *_fnit, zeroCoords, position );
      typename ConfiguratorType::VecType displacement = DeformMatrix * position + Shift - position;
      // save the result into the displacement MultiVector
      const int globalIndex = _config.localToGlobal ( *_fnit, 0 ) ;
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        Displacement[i][globalIndex] = displacement[i];
    }
  }

  /**
   * \todo Check if this works.
   *
   * \author Olischlaeger
   */
  void generatePiecewiseConstantFlowField ( aol::MultiVector<RealType> &FlowField, const aol::MultiVector<RealType> &LevelsetFunctions, const aol::MultiVector<RealType> &DeformationParameters ) {
    aol::Vector<int> intVec( LevelsetFunctions.numComponents() );
    const int dim = FlowField.numComponents();
    cerr<<"dim: "<<dim<<endl;
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const int globalIndex = _config.localToGlobal ( *_fnit, 0 );
      const int segmentIndex = aol::PowerSetIterator::getPositionNumberFromLevelsetFunctions<RealType>(LevelsetFunctions, globalIndex);
      for ( int i = 0; i < dim; i++ ){
        FlowField[i][globalIndex] =  DeformationParameters[segmentIndex][i];
      }
    }
  }

  void generateColoredSegmentation ( const int NumberOfLevelsetFunctions,
                                     const aol::MultiVector<RealType> &LevelsetFunctions,
                                     const VectorType &BackgroundImage,
                                     qc::MultiArray<RealType, 2, 3> &OutputImageArray,
                                     const RealType Isovalue = aol::NumberTrait<RealType>::zero ) {
    aol::Vector<int> intVec( NumberOfLevelsetFunctions );

    if ( (1<<NumberOfLevelsetFunctions) <= 16 ){
      RealType colorMap[16][3] = { {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}, {1., 1., 0.},
                                   {1., 0., 1.}, {0., 1., 1.}, {1., 1., 1.}, {1., 0.5, 0.},
                                   {1., 0., 0.5}, {0.5, 1., 0.}, {0.5, 0., 1.}, {0., 0.5, 1.},
                                   {0., 1., 0.5}, {1., 1., 0.5}, {1., 0.5, 1.}, {0.5, 1., 1.}
                                 };

      for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
        for ( int i = 0; i < NumberOfLevelsetFunctions; i++ )
          intVec[i] = (LevelsetFunctions[i].get ( _config.localToGlobal ( *_fnit, 0 ) )>= Isovalue ) ? 1 : 0;
        for ( int i = 0; i < 3; i++ )
          OutputImageArray[i].set ( *_fnit, BackgroundImage[_config.localToGlobal ( *_fnit, 0 )]*colorMap[aol::PowerSetIterator::getPositionNumberFromVector(intVec)][i] );
      }
    }
    else
      throw aol::UnimplementedCodeException( "Not implemented for more than 16 segments.\n", __FILE__, __LINE__ );
  }

  //! Sets all values in the stripes to 1., the rest is untouched.
  void generateStripes ( VectorType &Image, const int Width ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      if ( ( ( *_fnit ) [0] / Width ) % 2 == 0 )
        Image[ _config.localToGlobal ( *_fnit, 0 )  ] = 1.;
    }
  }
  //! Sets all values in white checkerboard fields to 1., the rest is untouched.
  void generateCheckerBoard ( VectorType &Image, const int Width ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      int checkValue = 0;
      for ( int i = 0; i < ConfiguratorType::Dim; ++i )
        checkValue += ( *_fnit ) [i] / Width;

      if ( ( checkValue ) % 2 == 0 )
        Image[ _config.localToGlobal ( *_fnit, 0 )  ] = 1.;
    }
  }

  //! Generates a slice or checkerboard view of two images.
  void generateCheckView ( VectorType &CheckViewImage, const VectorType &ImageA, const VectorType &ImageB, const int Width, const bool SliceView = true ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      int checkValue = 0;

      if( SliceView ){
        for ( int i = 0; i < ConfiguratorType::Dim; ++i )
          checkValue += ( *_fnit ) [i];
        checkValue /= Width;
      }
      else{
        for ( int i = 0; i < ConfiguratorType::Dim; ++i )
          checkValue += ( *_fnit ) [i] / Width;
      }
      const int globalIndex = _config.localToGlobal ( *_fnit, 0 );
      CheckViewImage[globalIndex] = ( ( checkValue ) % 2 == 0 ) ? ImageA[globalIndex] : ImageB[globalIndex];
    }
  }

  void generateCheckView ( aol::MultiVector<RealType> &CheckViewImage, const aol::MultiVector<RealType> &ImageA, const aol::MultiVector<RealType> &ImageB, const int Width, const bool SliceView = true ) {
    for ( int i = 0; i < CheckViewImage.numComponents(); ++i )
      generateCheckView ( CheckViewImage[i], ImageA[i], ImageB[i], Width, SliceView );
  }

  //! Sets all values inside the ellipse to Intensity, the rest is untouched.
  void generateEllipse ( VectorType &Image, const RealType Intensity = 1., const RealType a = 0.3, const RealType b = 0.5, const RealType CenterX = 0.5, const RealType CenterY = 0.5 ) {
    typename ConfiguratorType::ElementIteratorType it;
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      if ( ( aol::Sqr ( ( *_fnit ) [0] * _h - CenterX )/ aol::Sqr(a) + aol::Sqr ( ( *_fnit ) [1] * _h - CenterY )/ aol::Sqr(b) ) < 1.0 )
        Image[ _config.localToGlobal ( *_fnit, 0 )  ] = Intensity;
    }
  }
  
  //! Sets all values inside the circle to Intensity, the rest is untouched.
  void generateCircle ( VectorType &Image, const RealType Intensity = 1., const RealType Radius = 0.5, const RealType CenterX = 0.5, const RealType CenterY = 0.5 ) {
    generateEllipse(Image, Intensity, Radius, Radius, CenterX, CenterY);
  }

  //! Sets all values inside the rectangular to Intensity, the rest is untouched.
  void generateRectangular ( VectorType &Image, const RealType Intensity = 1., const RealType StartX = 0.25, const RealType StartY = 0.25, const RealType EndX = 0.75, const RealType EndY = 0.75 ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const RealType x = ( *_fnit ) [0] * _h;
      const RealType y = ( *_fnit ) [1] * _h;
      if ( x >= StartX && x <= EndX && y >=  StartY && y <= EndY )
        Image[ _config.localToGlobal ( *_fnit, 0 )  ] = Intensity;
    }
  }

  //! Sets Image to 0 on the boundary of [0,1]^d, to 1 in [TransitionWidth,1-TransitionWidth]^d and creates a smooth transition in between.
  void generateBlendingKernel ( VectorType &Image, const RealType TransitionWidth = 0.1 ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      typename ConfiguratorType::VecType pos;
      for ( int i = 0; i < ConfiguratorType::Dim; ++i )
        pos[i] = ( *_fnit ) [i] * _h;
      // Calculate the distance of pos to the boundary of the unit cube [0,1]^d
      const RealType dist = aol::Min ( pos.getMinValue(), 1 - pos.getMaxValue() );
      Image[ _config.localToGlobal ( *_fnit, 0 ) ] = aol::Clamp ( dist, aol::ZOTrait<RealType>::zero, TransitionWidth ) / TransitionWidth;
    }
  }

  //! Generates a phase field function to the boundary of a rectangular, centered in \f[\lbrack 0,1\rbrack^2\f]
  void generateRectangularPhaseField ( VectorType &Image, const RealType OffsetX, const RealType OffsetY, const RealType TransitionLayerWidth ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const RealType x = ( *_fnit ) [0] * _h;
      const RealType y = ( *_fnit ) [1] * _h;

      Image[ _config.localToGlobal ( *_fnit, 0 )  ] = std::tanh ( aol::Max(aol::Abs(x-0.5)-OffsetX, aol::Abs(y-0.5)-OffsetY) / TransitionLayerWidth );
    }
  }

  //! Generates a phase field function to the boundary of a cuboid, centered in \f[\lbrack 0,1\rbrack^3\f]
  void generateCuboidPhaseField ( VectorType &Image, const RealType OffsetX, const RealType OffsetY, const RealType OffsetZ, const RealType TransitionLayerWidth ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const RealType x = ( *_fnit ) [0] * _h;
      const RealType y = ( *_fnit ) [1] * _h;
      const RealType z = ( *_fnit ) [2] * _h;

      Image[ _config.localToGlobal ( *_fnit, 0 )  ] = std::tanh ( aol::Max(aol::Max(aol::Abs(x-0.5)-OffsetX, aol::Abs(y-0.5)-OffsetY), aol::Abs(z-0.5)-OffsetZ) / TransitionLayerWidth );
    }
  }

  //! Generates a phase field function to the boundary of a cuboid, centered in \f[\lbrack 0,1\rbrack^3\f]
  void generateDiscPhaseField ( VectorType &Image, const RealType OffsetR, const RealType OffsetZ, const RealType TransitionLayerWidth ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const RealType x = ( *_fnit ) [0] * _h;
      const RealType y = ( *_fnit ) [1] * _h;
      const RealType z = ( *_fnit ) [2] * _h;

      Image[ _config.localToGlobal ( *_fnit, 0 )  ] = std::tanh ( aol::Max(std::sqrt(aol::Sqr(x-0.5)+aol::Sqr(y-0.5))-OffsetR, aol::Abs(z-0.5)-OffsetZ) / TransitionLayerWidth );
    }
  }

  //! Generates a phase field function to the boundary of two circles, the center of the first circle is \f[(CenterX1,0.5)\f] and of the second one \f[(CenterX2,0.5)\f].
  void generateTwoCirclesPhaseField ( VectorType &Image, const RealType CenterX1, const RealType CenterX2, const RealType Radius, const RealType TransitionLayerWidth ) {
  RealType tmp = 0.;
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const RealType x = ( *_fnit ) [0] * _h;
      const RealType y = ( *_fnit ) [1] * _h;

      tmp = aol::Min(std::sqrt(aol::Sqr(x-CenterX1)+aol::Sqr(y-0.5))-Radius, std::sqrt(aol::Sqr(x-CenterX2)+aol::Sqr(y-0.5))-Radius);
      Image[ _config.localToGlobal ( *_fnit, 0 )  ] = std::tanh ( tmp / TransitionLayerWidth );
    }
  }

  //! Generates a phase field function to the boundary of two balls, the center of the first ball is \f[(CenterX1,0.5, 0.5)\f] and of the second one \f[(CenterX2,0.5,0.5)\f].
  void generateTwoBallsPhaseField ( VectorType &Image, const RealType CenterX1, const RealType CenterX2, const RealType Radius, const RealType TransitionLayerWidth ) {
  RealType tmp = 0.;
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const RealType x = ( *_fnit ) [0] * _h;
      const RealType y = ( *_fnit ) [1] * _h;
      const RealType z = ( *_fnit ) [2] * _h;

      tmp = aol::Min(std::sqrt(aol::Sqr(x-CenterX1)+aol::Sqr(y-0.5)+aol::Sqr(z-0.5))-Radius, std::sqrt(aol::Sqr(x-CenterX2)+aol::Sqr(y-0.5)+aol::Sqr(z-0.5))-Radius);
      Image[ _config.localToGlobal ( *_fnit, 0 )  ] = std::tanh ( tmp / TransitionLayerWidth );
    }
  }

  void generateSignedDistanceFunctionOfTorus ( VectorType &Image,
                                               const typename ConfiguratorType::VecType &Center,
                                               const RealType RadiusCenterline,
                                               const RealType Radius ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const RealType x = ( *_fnit ) [0] * _h - Center[0];
      const RealType y = ( *_fnit ) [1] * _h - Center[1];
      const RealType z = ( *_fnit ) [2] * _h - Center[2];

      Image[ _config.localToGlobal ( *_fnit, 0 )  ] = std::sqrt( aol::Sqr( RadiusCenterline - std::sqrt( x*x + y*y ) ) + aol::Sqr( z ) )-Radius;
    }
  }

  void generateSignedDistanceFunctionOfRotatedTorus ( VectorType &Image,
                                                      const typename ConfiguratorType::VecType &Center,
                                                      const RealType RadiusCenterline,
                                                      const RealType Radius,
                                                      const RealType Angle ) {
    RealType cosAngle = cos( Angle );
    RealType sinAngle = sin( Angle );
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const RealType x = ( *_fnit ) [0] * _h - Center[0];
      const RealType y = ( *_fnit ) [1] * _h - Center[1];
      const RealType z = ( *_fnit ) [2] * _h - Center[2];
      const RealType yr = cosAngle * y - sinAngle * z;
      const RealType zr = sinAngle * y + cosAngle * z;

      Image[ _config.localToGlobal ( *_fnit, 0 )  ] = std::sqrt( aol::Sqr( RadiusCenterline - std::sqrt( x*x + yr*yr ) ) + aol::Sqr( zr ) )-Radius;
    }
  }

  //! Generates a mask for each segment identified by the signs of the levelset functions.
  void generateMaskFromLevelsetfunctions ( std::vector<qc::BitArray<qc::QC_2D>*> &MaskVector,
                                           const aol::MultiVector<RealType> &LevelsetFunctions,
                                           const RealType Isovalue = aol::NumberTrait<RealType>::zero ) {
    aol::Vector<int> intVec( LevelsetFunctions.numComponents() );
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      for ( int i = 0; i < LevelsetFunctions.numComponents(); i++ )
        intVec[i] = (LevelsetFunctions[i].get ( _config.localToGlobal ( *_fnit, 0 )  )>= Isovalue ) ? 1 : 0;
      MaskVector[aol::PowerSetIterator::getPositionNumberFromVector(intVec)]->set ( *_fnit, true );
    }
  }

  //! Calculates the norm of the gradient of InputImage at every node using finite differences.
  void generateGradientNormFD ( const qc::ScalarArray<RealType, qc::QC_2D> &InputImage, VectorType &GradientNorm ) {
    aol::Vec2<RealType> grad;
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const aol::Vec2<int> pos ( ( *_fnit ) [0], ( *_fnit ) [1] );
      InputImage.gradientFD ( pos, grad );
      grad /= _h;
      GradientNorm[ _config.localToGlobal ( *_fnit, 0 ) ] = grad.norm();
    }
  }

  void generateArtificialPFCFront ( VectorType &Image,
                                    const RealType DomainScaleFactor,
                                    const RealType FrontPosition,
                                    const RealType TransitionWidth,
                                    const bool CircularFront ) {

    const ArtificialPFCFrontFunction<RealType> f ( DomainScaleFactor, FrontPosition, TransitionWidth, CircularFront );

    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const aol::Vec2<RealType> x ( ( *_fnit ) [0] * _h, ( *_fnit ) [1] * _h );
      Image[ _config.localToGlobal ( *_fnit, 0 ) ] = f.evaluate ( x );
    }
  }

  template <typename FunctionType, typename ParametricDeformationType>
  void generateDeformedFunctionImage ( VectorType &Image, const FunctionType &F, const ParametricDeformationType &ParDef ) {
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      const aol::Vec2<RealType> x ( ( *_fnit ) [0] * _h, ( *_fnit ) [1] * _h );
      aol::Vec2<RealType> defX;
      ParDef.evaluateDeformationOn01 ( ParDef, x, defX );
      Image[ _config.localToGlobal ( *_fnit, 0 ) ] = F.evaluate ( defX );
    }
  }

  void generateHashSetFromSubLevelset ( const VectorType &LevelsetFunction,
                                        const RealType IsoValue,
                                        aol::HashSet<typename ConfiguratorType::ElementType> &Hash ) {
    Hash.clear( );
    for ( _fnit = _grid._nBeginIt; _fnit != _grid._nEndIt; ++_fnit ) {
      if ( (LevelsetFunction[_config.localToGlobal ( *_fnit, 0 )] ) < IsoValue )
        Hash.insert ( *_fnit );
    }
  }

};

/**
 * \brief Creates the characteristic function corresponding to the subgraph of a given function.
 *
 * Assumes the height of the domain to be 1 and that the function values are relative to 1.
 *
 * \note Doesn't fit in qc::DataGenerator.
 *
 * \author Berkels
 */
template <typename RealType>
void generateSubgraphVolume ( const qc::ScalarArray<RealType, qc::QC_2D> &Function, qc::ScalarArray<RealType, qc::QC_3D> &SubgraphVolume, const bool AntiAlias = true ) {
  if ( ( Function.getNumX() != SubgraphVolume.getNumX() ) || ( Function.getNumY() != SubgraphVolume.getNumY() ) )
      throw aol::DimensionMismatchException( "Sizes of the arguments don't match.\n", __FILE__, __LINE__ );

  const int numX = SubgraphVolume.getNumX();
  const int numY = SubgraphVolume.getNumY();
  const int numZ = SubgraphVolume.getNumZ();

  SubgraphVolume.setZero();
  for ( int y = 0; y < numY; ++y ) {
    for ( int x = 0; x < numX; ++x ) {
      const int zMax = aol::Min( numZ, static_cast<int> ( floor ( Function.get ( x, y ) * numZ ) ) );
      for ( int z = 0; z < zMax; ++z )
        SubgraphVolume.set ( x, y, z, 1 );

      if ( AntiAlias && ( zMax < numZ - 1 ) && ( zMax >= 0 ) )
        SubgraphVolume.set ( x, y, zMax, Function.get ( x, y ) * numZ -zMax );
    }
  }
}

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class SmoothedImage {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const typename ConfiguratorType::InitType &_grid;
  typename ConfiguratorType::ArrayType _image;
  typename ConfiguratorType::ArrayType _imageDx;
  typename ConfiguratorType::ArrayType _imageDy;
  const RealType _sigma;
  const int _kernelSize;
public:
  SmoothedImage ( const typename ConfiguratorType::InitType &Grid,
                  const aol::Vector<RealType> &Image,
                  const bool ComputeSmoothDerivative,
                  const RealType Sigma = 1.5,
                  const int KernelSize = 5 )
    : _grid ( Grid ),
      _image ( Grid ),
      _imageDx ( ),
      _imageDy ( ),
      _sigma ( Sigma ),
      _kernelSize ( KernelSize ) {
    const typename ConfiguratorType::ArrayType image ( Image, Grid, aol::FLAT_COPY );
    qc::GaussKernel2d<RealType> kernel ( _kernelSize, _sigma );
    image.applyLinearFilterTo ( kernel, _image );
    if ( ComputeSmoothDerivative ) {
      _imageDx.reallocate ( Grid );
      _imageDy.reallocate ( Grid );
      qc::GaussDiffKernel2d<RealType> kernelX ( _kernelSize, _sigma, qc::DIFF_X );
      qc::GaussDiffKernel2d<RealType> kernelY ( _kernelSize, _sigma, qc::DIFF_Y );
      image.applyLinearFilterTo ( kernelX, _imageDx );
      image.applyLinearFilterTo ( kernelY, _imageDy );
    }
  }

  const typename ConfiguratorType::ArrayType &getSmoothImRef ( ) const {
    return _image;
  }

  const typename ConfiguratorType::ArrayType &getSmoothDXImRef ( ) const {
    return _imageDx;
  }

  const typename ConfiguratorType::ArrayType &getSmoothDYImRef ( ) const {
    return _imageDy;
  }

  RealType interpolate ( const typename ConfiguratorType::VecType &Pos ) const {
    return _image.interpolate ( Pos );
  }

  RealType interpolateDX ( const typename ConfiguratorType::VecType &Pos ) const {
    return _imageDx.interpolate ( Pos );
  }

  RealType interpolateDY ( const typename ConfiguratorType::VecType &Pos ) const {
    return _imageDy.interpolate ( Pos );
  }
};

/**
 * \author Berkels
 */
template <typename ScalarArrayType>
class ScalarArrayHelper {
  typedef typename ScalarArrayType::DataType DataType;
  typedef typename ScalarArrayType::RealType RealType;
  static const qc::Dimension Dim = ScalarArrayType::Dim;
public:
  static void getPointWiseXDerivNorm ( const qc::MultiArray<DataType, Dim> &MArray, ScalarArrayType &XDNorm ) {
    qc::RectangularGrid<Dim> grid ( qc::GridSize<Dim>::createFrom ( MArray ) );
    for ( qc::RectangularIterator<Dim> it ( grid ); it.notAtEnd(); ++it ) {
      aol::Vec<Dim, RealType> tempVec;
      for ( int comp = 0; comp < Dim; ++comp )
        tempVec[comp] = MArray[comp].dxFD( *it ) / grid.H();
      XDNorm.set ( *it, tempVec.norm() );
    }
  }

  static void getXDerivative ( const ScalarArrayType &Array, ScalarArrayType &XDeriv ) {
    qc::RectangularGrid<Dim> grid ( qc::GridSize<Dim>::createFrom ( Array ) );
    for ( qc::RectangularIterator<Dim> it ( grid ); it.notAtEnd(); ++it ) {
      XDeriv.set ( *it, Array.dxFD( *it ) / grid.H() );
    }
  }

  static void getYDerivative ( const ScalarArrayType &Array, ScalarArrayType &YDeriv ) {
    qc::RectangularGrid<Dim> grid ( qc::GridSize<Dim>::createFrom ( Array ) );
    for ( qc::RectangularIterator<Dim> it ( grid ); it.notAtEnd(); ++it ) {
      YDeriv.set ( *it, Array.dyFD( *it ) / grid.H() );
    }
  }

  static void getSmoothedGradient ( const qc::ScalarArray<RealType, qc::QC_2D> &Array,
                                    qc::MultiArray<RealType, qc::QC_2D> &Gradient,
                                    const RealType Sigma = 5,
                                    const int KernelSize = 9 ) {
    qc::ScalarArray<RealType, qc::QC_2D> extendedData ( 2 * Array.getNumX(), 2 * Array.getNumY() );
    qc::ScalarArrayHelper<qc::ScalarArray<RealType, qc::QC_2D> >::extendMirrored ( Array, extendedData );
    qc::ScalarArray<RealType, qc::QC_2D> temp ( extendedData, aol::STRUCT_COPY );
    temp.shiftByOffset ( extendedData.getNumX() / 2, extendedData.getNumY() / 2, extendedData );

    typedef qc::RectangularGridConfigurator<RealType, qc::QC_2D, aol::GaussQuadrature<RealType,qc::QC_2D,3> > ConfType;
    const typename ConfType::InitType grid ( qc::GridSize2d::createFrom ( extendedData ) );
    qc::SmoothedImage<ConfType> smoothedData ( grid, temp, true, Sigma, KernelSize );
    smoothedData.getSmoothDXImRef ( ).copyBlockTo ( extendedData.getNumX() / 2, extendedData.getNumY() / 2, Gradient[0] );
    smoothedData.getSmoothDYImRef ( ).copyBlockTo ( extendedData.getNumX() / 2, extendedData.getNumY() / 2, Gradient[1] );
    const typename ConfType::InitType dataGrid ( qc::GridSize2d::createFrom ( Array ) );
    Gradient /=  dataGrid.H();
  }

  //! \note straightforward, but not very efficient implementation
  static void erode ( const qc::ScalarArray<DataType, qc::QC_2D> &InputArray, qc::ScalarArray<DataType, qc::QC_2D> &ErodedArray, const int InfRadius ) {
    for ( int y = 0; y < InputArray.getNumY(); ++y ) {
      for ( int x = 0; x < InputArray.getNumX(); ++x ) {
        const int XMin = aol::Max ( 0, x - InfRadius );
        const int YMin = aol::Max ( 0, y - InfRadius );
        const int XMax = aol::Min ( InputArray.getNumX() - 1, x + InfRadius );
        const int YMax = aol::Min ( InputArray.getNumY() - 1, y + InfRadius );

        qc::ScalarArray<DataType, qc::QC_2D> block ( ( XMax - XMin + 1 ), ( YMax - YMin + 1 ) );
        InputArray.copyBlockTo( XMin, YMin, block );
        ErodedArray.set ( x, y, block.getMinValue() );
      }
    }
  }

  //! \note straightforward, but not very efficient implementation
  static void dilate ( const qc::ScalarArray<DataType, qc::QC_2D> &InputArray, qc::ScalarArray<DataType, qc::QC_2D> &ErodedArray, const int InfRadius ) {
    for ( int y = 0; y < InputArray.getNumY(); ++y ) {
      for ( int x = 0; x < InputArray.getNumX(); ++x ) {
        const int XMin = aol::Max ( 0, x - InfRadius );
        const int YMin = aol::Max ( 0, y - InfRadius );
        const int XMax = aol::Min ( InputArray.getNumX() - 1, x + InfRadius );
        const int YMax = aol::Min ( InputArray.getNumY() - 1, y + InfRadius );

        qc::ScalarArray<DataType, qc::QC_2D> block ( ( XMax - XMin + 1 ), ( YMax - YMin + 1 ) );
        InputArray.copyBlockTo( XMin, YMin, block );
        ErodedArray.set ( x, y, block.getMaxValue() );
      }
    }
  }

  static void extendMirrored ( const qc::ScalarArray<DataType, qc::QC_2D> &InputArray, qc::ScalarArray<DataType, qc::QC_2D> &ExtendedArray ) {
    const int numX = InputArray.getNumX();
    const int numY = InputArray.getNumY();

    QUOC_ASSERT ( ExtendedArray.getNumX() == 2 * numX );
    QUOC_ASSERT ( ExtendedArray.getNumY() == 2 * numY );

    for ( int j = 0; j < numY; ++j ){
      for ( int i = 0; i < numX; ++i ){
        const DataType value = InputArray.get( i, j );
        ExtendedArray.set( i, j, value );
        ExtendedArray.set( 2 * numX - 1 - i, j, value );
        ExtendedArray.set( i, 2 * numY - 1 - j, value );
        ExtendedArray.set( 2 * numX - 1 - i, 2 * numY - 1 - j, value );
      }
    }
  }
};


} // end namespace qc

#endif // __GENERATOR_H
