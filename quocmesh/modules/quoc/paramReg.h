#ifndef __PARAMREG_H
#define __PARAMREG_H

#include <registration.h>
#include <convolution.h>
#include <tensor.h>
#include <regression.h>

namespace qc {

/** 
 * \author Berkels
 */
template <bool AllowScaling, bool AllowTranslation>
class ParametricRigidBodyMotion2DHelper {};

template <>
class ParametricRigidBodyMotion2DHelper<false, true> {
public:
  static const int NumberOfScalingParameters = 0;
  static aol::Vec<2, int> getDeformParametersSize ( ) {
    return aol::Vec2<int> ( 2, 1 );
  }
};

template <>
class ParametricRigidBodyMotion2DHelper<false, false> {
public:
  static const int NumberOfScalingParameters = 0;
  static aol::Vec<1, int> getDeformParametersSize ( ) {
    return aol::Vec<1, int> ( 1 );
  }
};

template <>
class ParametricRigidBodyMotion2DHelper<true, true> {
public:
  static const int NumberOfScalingParameters = 1;
  static aol::Vec<3, int> getDeformParametersSize ( ) {
    return aol::Vec3<int> ( 2, 1, 1 );
  }
};
  
template <>
class ParametricRigidBodyMotion2DHelper<true, false> {
public:
  static const int NumberOfScalingParameters = 1;
  static aol::Vec<2, int> getDeformParametersSize ( ) {
    return aol::Vec2<int> ( 1, 1 );
  }
};

/**
  * \author Tatano
 */
template <bool AllowScaling, bool NumberScalingFactorIsOne>
class ParametricRigidBodyMotion3DHelper {};
  
template <>
class ParametricRigidBodyMotion3DHelper<false, false> {
public:
  static const int NumberOfScalingParameters = 0;
  static aol::Vec<2, int> getDeformParametersSize ( ) {
    return aol::Vec2<int> ( 3, 3 );
  }
};

template <>//just repetition to have it compatible with tamplate
class ParametricRigidBodyMotion3DHelper<false, true> {
public:
  static const int NumberOfScalingParameters = 0;
  static aol::Vec<2, int> getDeformParametersSize ( ) {
    return aol::Vec2<int> ( 3, 3 );
  }
};
  
template <>
class ParametricRigidBodyMotion3DHelper<true, true> {
public:
  static const int NumberOfScalingParameters = 1;
  static aol::Vec<3, int> getDeformParametersSize ( ) {
    return aol::Vec3<int> ( 3, 3, 1 );
  }
};
  
template <>
class ParametricRigidBodyMotion3DHelper<true, false> {
public:
  static const int NumberOfScalingParameters = 3;
  static aol::Vec<3, int> getDeformParametersSize ( ) {
    return aol::Vec3<int> ( 3, 3, 3 );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class ParametricDeformationBase {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_grid;
public:
  ParametricDeformationBase ( const typename ConfiguratorType::InitType &Initializer )
    : _grid ( Initializer ) { }

  //! \note Inefficient implementation, don't use for anything where performance is critical.
  template <typename ParametricDeformationType>
  void evaluateDeformationOn01 ( const ParametricDeformationType &ParDef, const typename ConfiguratorType::VecType &Position, typename ConfiguratorType::VecType &DeformedPosition ) const {
    typename ConfiguratorType::ElementType el;
    typename ConfiguratorType::DomVecType localCoord;
    qc::getLocalCoordsRegularRectangularGrid<ConfiguratorType> ( Position, _grid, el, localCoord );
    typename ConfiguratorType::ElementType transformedEl;
    typename ConfiguratorType::VecType transformedLocalCoord;
    ParDef.template evaluateDeformation<false> ( el, localCoord, transformedEl, transformedLocalCoord );

    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      DeformedPosition[i] = ( transformedEl[i] + transformedLocalCoord[i] ) * _grid.H();
  }

  const typename ConfiguratorType::InitType &getGridRef () const {
    return _grid;
  }
};

/**
 * \note With AllowScaling == true the name is misleading: The class is a rigid body motion plus scaling in this case.
 *
 * \author Berkels
 */
template <typename _ConfiguratorType, bool AllowScaling = false, bool AllowTranslation = true>
class ParametricRigidBodyMotion2D {
public:
  typedef _ConfiguratorType ConfiguratorType;
  static const int Dim = ConfiguratorType::Dim;
  static const int NumberOfDeformParameters = 1 + ( AllowTranslation ? ConfiguratorType::Dim : 0 ) + ParametricRigidBodyMotion2DHelper<AllowScaling, AllowTranslation>::NumberOfScalingParameters;
  static const int NumOfDeformParametersComponents = 1 + ( AllowTranslation ? 1 : 0 ) + ParametricRigidBodyMotion2DHelper<AllowScaling, AllowTranslation>::NumberOfScalingParameters;
  typedef typename ConfiguratorType::RealType RealType;
private:
  static const int RotationVecIndex = AllowTranslation ? ConfiguratorType::Dim : 0;
  static const int RotationMVecIndex = AllowTranslation ? 1 : 0;
  const typename ConfiguratorType::InitType &_grid;
  RealType _rotationAngle;
  aol::Matrix22<RealType> _rotationMatrix;
  typename ConfiguratorType::VecType _translation;
  RealType _scaling;
public:
  ParametricRigidBodyMotion2D ( const typename ConfiguratorType::InitType &Initializer, const aol::MultiVector<RealType> &DeformParameters )
    : _grid ( Initializer ),
      _rotationAngle ( DeformParameters[RotationMVecIndex][0] ),
      _scaling ( AllowScaling ? DeformParameters[RotationMVecIndex+1][0] : 1 ) {
    _rotationMatrix.makeRotation ( _rotationAngle );
    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      _translation[i] = ( AllowTranslation ? DeformParameters[0][i] : 0 ) + 0.5;
  }

  ParametricRigidBodyMotion2D ( const typename ConfiguratorType::InitType &Initializer, const aol::Vector<RealType> &DeformParametersVec )
    : _grid ( Initializer ) {
    setDeformationParametersAsVec ( DeformParametersVec );
  }

  ParametricRigidBodyMotion2D ( const typename ConfiguratorType::InitType &Initializer )
    : _grid ( Initializer ),
      _rotationAngle ( 0 ),
      _scaling ( 1 ) {
    _rotationMatrix.makeRotation ( _rotationAngle );
    _translation.setAll ( 0.5 );
  }

  void setDeformationParametersAsVec ( const aol::Vector<RealType> &DeformParametersVec ) {
    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      _translation[i] = ( AllowTranslation ? DeformParametersVec[i] : 0 ) + 0.5;
    _rotationAngle = DeformParametersVec[RotationVecIndex];
    _scaling = ( AllowScaling ? DeformParametersVec[RotationVecIndex+1] : 1 );
    _rotationMatrix.makeRotation ( _rotationAngle );
  }

  static void setIdentityDeformationParameters ( aol::MultiVector<RealType> &DeformParameters ) {
    DeformParameters.setZero();
    if ( AllowScaling )
      DeformParameters[RotationMVecIndex + 1][0] = 1;
  }

  static void setIdentityDeformationParametersAsVec ( aol::Vector<RealType> &DeformParametersAsVec ) {
    DeformParametersAsVec.setZero();
    if ( AllowScaling )
      DeformParametersAsVec[RotationVecIndex + 1] = 1;
  }

  // DestDeform := ArgDeform \circ DestDeform
  static void concatenateDeformationParameters ( const aol::MultiVector<RealType> &ArgDeformPars, aol::MultiVector<RealType> &DestDeformPars ) {
    if ( AllowScaling )
      DestDeformPars[RotationMVecIndex+1][0] *= ArgDeformPars[RotationMVecIndex+1][0];

    if ( AllowTranslation ) {
      aol::Matrix22<RealType > _rotationMatrix;
      aol::Vec2<RealType> tmp ( ArgDeformPars[0][0], ArgDeformPars[0][1] );
      _rotationMatrix.makeRotation ( ArgDeformPars[RotationMVecIndex][0] );
      _rotationMatrix.multAdd ( DestDeformPars[0], tmp );
      if ( AllowScaling )
        tmp *= ArgDeformPars[RotationMVecIndex+1][0];

      for ( int i = 0; i < 2; ++i )
        DestDeformPars[0][i] = tmp[i];
    }

    DestDeformPars[RotationMVecIndex][0] += ArgDeformPars[RotationMVecIndex][0];
  }

  static void getDeformationParametersOfInverse ( const aol::MultiVector<RealType> &DeformPars, aol::MultiVector<RealType> &InverseDeformPars )  {
    InverseDeformPars[RotationMVecIndex][0] = -DeformPars[RotationMVecIndex][0];
    if ( AllowScaling )
      InverseDeformPars[RotationMVecIndex + 1][0] = aol::ZOTrait<RealType>::one / DeformPars[RotationMVecIndex + 1][0];

    if ( AllowTranslation ) {
      aol::Matrix22<RealType> inverseRotationMatrix;
      inverseRotationMatrix.makeRotation ( InverseDeformPars[RotationMVecIndex][0] );
      aol::Vec2<RealType> tempVec, inverseTranslation;

      for ( int i = 0; i < ConfiguratorType::Dim; ++i )
        tempVec[i] = -DeformPars[0][i];

      inverseRotationMatrix.mult ( tempVec, inverseTranslation );

      if ( AllowScaling )
        inverseTranslation *= InverseDeformPars[RotationMVecIndex + 1][0];

      for ( int i = 0; i < ConfiguratorType::Dim; ++i )
        InverseDeformPars[0][i] = inverseTranslation[i];
    }
  }

  static aol::Vec<NumOfDeformParametersComponents, int> getDeformParametersSize ( ) {
    return ParametricRigidBodyMotion2DHelper<AllowScaling, AllowTranslation>::getDeformParametersSize();
  }

  template <bool ClipCoord>
  bool evaluateDeformation ( const typename ConfiguratorType::ElementType &El,
                             const typename ConfiguratorType::VecType &RefCoord,
                             qc::Element &TransformedEl, typename ConfiguratorType::VecType &TransformedLocalCoord ) const {

    typename ConfiguratorType::VecType coord, transformedCoord;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      coord[i] = El[i] + RefCoord[i] - 0.5 / _grid.H();
    }
    _rotationMatrix.mult ( coord, transformedCoord );
    if ( AllowScaling )
      transformedCoord *= _scaling;
    return qc::transformCoord<ConfiguratorType, ClipCoord> ( _grid, transformedCoord, _translation, TransformedEl, TransformedLocalCoord );
  }

  void evaluateDeformationOn01 ( const ParametricRigidBodyMotion2D<ConfiguratorType, AllowScaling, AllowTranslation> &/*ParDef*/, const typename ConfiguratorType::VecType &Position, typename ConfiguratorType::VecType &DeformedPosition ) const {
    typename ConfiguratorType::VecType tmp ( Position );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      tmp[i] -= 0.5;
    _rotationMatrix.mult ( tmp, DeformedPosition );
    if ( AllowScaling )
      DeformedPosition *= _scaling;
    DeformedPosition += _translation;
  }

  void evaluateDeformationOn0N ( typename aol::Vec<ConfiguratorType::Dim, RealType> &DeformedPosition ) const {
    typename ConfiguratorType::VecType tmp ( DeformedPosition );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      tmp[i] -= 0.5 / this->_grid.H();
    _rotationMatrix.mult ( tmp, DeformedPosition );
    if ( AllowScaling )
      DeformedPosition *= _scaling;
    DeformedPosition.addMultiple ( _translation, 1 / this->_grid.H() );
  }

  void evaluateDerivativeDeformationOn01 ( const typename ConfiguratorType::VecType &Position,
                                           aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    const RealType cosAlpha = _rotationMatrix[0][0];
    const RealType sinAlpha = _rotationMatrix[1][0];
    Jacobian[RotationVecIndex][0] = ( AllowScaling ? _scaling : 1 ) * ( -sinAlpha * ( Position[0] - 0.5 ) - cosAlpha * ( Position[1] - 0.5 ) );
    Jacobian[RotationVecIndex][1] = ( AllowScaling ? _scaling : 1 ) * (  cosAlpha * ( Position[0] - 0.5 ) - sinAlpha * ( Position[1] - 0.5 ) );

    if ( AllowTranslation ) {
      Jacobian[0][0] = Jacobian[1][1] = 1;
      Jacobian[1][0] = Jacobian[0][1] = 0;
    }

    if ( AllowScaling ) {
      Jacobian[RotationVecIndex+1][0] = ( cosAlpha * ( Position[0] - 0.5 ) - sinAlpha * ( Position[1] - 0.5 ) );
      Jacobian[RotationVecIndex+1][1] = ( sinAlpha * ( Position[0] - 0.5 ) + cosAlpha * ( Position[1] - 0.5 ) );
    }
  }

  void evaluateDerivativeDeformation ( const typename ConfiguratorType::ElementType &El,
                                       const typename ConfiguratorType::VecType &RefCoord,
                                       aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    typename ConfiguratorType::VecType x;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      x[i] = ( El[i] + RefCoord[i] ) * _grid.H();
    }
    evaluateDerivativeDeformationOn01 ( x, Jacobian );
  }

  void evaluateDerivativeDeformationOn0N ( const typename ConfiguratorType::VecType &Position,
                                          aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    const RealType cosAlpha = _rotationMatrix[0][0];
    const RealType sinAlpha = _rotationMatrix[1][0];
    const RealType center0N = 0.5 / this->_grid.H();
    Jacobian[RotationVecIndex][0] = ( AllowScaling ? _scaling : 1 ) * ( -sinAlpha * ( Position[0] - center0N ) - cosAlpha * ( Position[1] - center0N ) );
    Jacobian[RotationVecIndex][1] = ( AllowScaling ? _scaling : 1 ) * (  cosAlpha * ( Position[0] - center0N ) - sinAlpha * ( Position[1] - center0N ) );

    if ( AllowTranslation ) {
      Jacobian[0][0] = Jacobian[1][1] = aol::ZOTrait<RealType>::one / this->_grid.H();
      Jacobian[1][0] = Jacobian[0][1] = 0;
    }

    if ( AllowScaling ) {
      Jacobian[RotationVecIndex+1][0] = ( cosAlpha * ( Position[0] - center0N ) - sinAlpha * ( Position[1] - center0N ) );
      Jacobian[RotationVecIndex+1][1] = ( sinAlpha * ( Position[0] - center0N ) + cosAlpha * ( Position[1] - center0N ) );
    }
  }


  void evaluateSecondDerivativeDeformationOn0N ( const typename ConfiguratorType::VecType &Position,
                                                 aol::Tensor3rdOrder< NumberOfDeformParameters, NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Hessian ) const {
    const RealType cosAlpha = _rotationMatrix[0][0];
    const RealType sinAlpha = _rotationMatrix[1][0];
    const RealType center0N = 0.5 / this->_grid.H();

    for ( int i = 0; i < Dim; i++ ) {
      if ( AllowTranslation )
        Hessian[i].setZero();
      Hessian[RotationVecIndex][i].setZero();
    }

    Hessian[RotationVecIndex][RotationVecIndex][0] = ( -cosAlpha * ( Position[0] - center0N ) + sinAlpha * ( Position[1] - center0N ) );
    Hessian[RotationVecIndex][RotationVecIndex][1] = ( -sinAlpha * ( Position[0] - center0N ) - cosAlpha * ( Position[1] - center0N ) );

    if ( AllowScaling ) {
      for ( int i = 0; i < Dim; i++ )
        Hessian[RotationVecIndex+1][i].setZero();
      Hessian[RotationVecIndex][RotationVecIndex+1][0] = Hessian[RotationVecIndex+1][RotationVecIndex][0] = ( -sinAlpha * ( Position[0] - center0N ) - cosAlpha * ( Position[1] - center0N ) );
      Hessian[RotationVecIndex][RotationVecIndex+1][1] = Hessian[RotationVecIndex+1][RotationVecIndex][1] = (  cosAlpha * ( Position[0] - center0N ) - sinAlpha * ( Position[1] - center0N ) );
      Hessian[RotationVecIndex+1][RotationVecIndex+1].setZero();
    }
  }

  void evaluateSpatialJacobianOn0N ( const typename ConfiguratorType::VecType &/*Position*/,
                                     aol::Mat< ConfiguratorType::Dim, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    Jacobian = _rotationMatrix;
    if ( AllowScaling )
      Jacobian *= _scaling;
  }

  void evaluateMixedHessianOn0N ( const typename ConfiguratorType::VecType &/*Position*/,
                                  aol::Tensor3rdOrder< NumberOfDeformParameters, ConfiguratorType::Dim, ConfiguratorType::Dim, RealType > &Hessian ) const {
    const RealType cosAlpha = _rotationMatrix[0][0];
    const RealType sinAlpha = _rotationMatrix[1][0];

    if ( AllowTranslation ) {
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        Hessian[i].setZero();
    }

    Hessian[RotationVecIndex][0][0] = -sinAlpha;
    Hessian[RotationVecIndex][0][1] = -cosAlpha;
    Hessian[RotationVecIndex][1][0] =  cosAlpha;
    Hessian[RotationVecIndex][1][1] = -sinAlpha;

    if ( AllowScaling ) {
      Hessian[RotationVecIndex] *= _scaling;
      Hessian[RotationVecIndex+1] = _rotationMatrix;
    }
  }

  void evaluateSpatialHessianOn0N ( const typename ConfiguratorType::VecType &/*Position*/,
                                    aol::Tensor3rdOrder< ConfiguratorType::Dim, ConfiguratorType::Dim, ConfiguratorType::Dim, RealType > &Hessian ) const {
    Hessian.setZero();
  }
};

/**
 * \note With AllowScaling == true the name is misleading: The class is a rigid body motion plus scaling in this case. If NumberScalingFactorIsOne == true, one scaling factor is used, othewise, three scaling factors are considered (s_x, s_y, s_z).
 *
 * \author Berkels
 */
template <typename ConfiguratorType, bool AllowScaling = false, bool NumberScalingFactorIsOne = true>
class ParametricRigidBodyMotion3D : public ParametricDeformationBase<ConfiguratorType> {
public:
  static const int NumberOfDeformParameters = 3 + ConfiguratorType::Dim + ParametricRigidBodyMotion3DHelper<AllowScaling, NumberScalingFactorIsOne>::NumberOfScalingParameters;
  static const int NumOfDeformParametersComponents = 2 + ((AllowScaling) ? 1 : 0);
  typedef typename ConfiguratorType::RealType RealType;
private:
  const RealType _yaw;
  const RealType _pitch;
  const RealType _roll;
  aol::Matrix33<RealType> _rotationYaw;
  aol::Matrix33<RealType> _rotationPitch;
  aol::Matrix33<RealType> _rotationRoll;
  typename ConfiguratorType::VecType _translation;
  typename ConfiguratorType::VecType _scalingVec;
  RealType _scaling;
public:
  ParametricRigidBodyMotion3D ( const typename ConfiguratorType::InitType &Initializer, const aol::MultiVector<RealType> &DeformParameters )
    : ParametricDeformationBase<ConfiguratorType> ( Initializer ),
      _yaw ( DeformParameters[1][0] ),
      _pitch ( DeformParameters[1][1] ),
      _roll ( DeformParameters[1][2] ),
      _scaling ( ( AllowScaling && NumberScalingFactorIsOne ) ? DeformParameters[2][0] : 1 ) {
      _rotationYaw.setRotationAboutZ ( _yaw );
      _rotationPitch.setRotationAboutY ( _pitch );
      _rotationRoll.setRotationAboutX ( _roll );
      for ( int i = 0; i < ConfiguratorType::Dim; ++i ){
        _scalingVec[i] = ( AllowScaling && (NumberScalingFactorIsOne == false) ) ? DeformParameters[2][i] : 1;
        _translation[i] = DeformParameters[0][i] + 0.5;
      }
      
  }

  static void setIdentityDeformationParameters ( aol::MultiVector<RealType> &DeformParameters ) {
    DeformParameters.setZero();
    if ( AllowScaling && NumberScalingFactorIsOne )
      DeformParameters[2][0] = 1;
    else if ( AllowScaling && (NumberScalingFactorIsOne == false) )
      for ( int i = 0; i < ConfiguratorType::Dim; ++i )
        DeformParameters[2][i] = 1;
  }

  static void concatenateDeformationParameters ( const aol::MultiVector<RealType> &/*ArgDeformPars*/, aol::MultiVector<RealType> &/*DestDeformPars*/ ) {
    throw aol::UnimplementedCodeException ( "not implemented", __FILE__, __LINE__ );
  }

  static aol::Vec<NumOfDeformParametersComponents, int> getDeformParametersSize ( ) {
    return ParametricRigidBodyMotion3DHelper<AllowScaling, NumberScalingFactorIsOne>::getDeformParametersSize();
  }

  template <bool ClipCoord>
  bool evaluateDeformation ( const typename ConfiguratorType::ElementType &El,
                             const typename ConfiguratorType::VecType &RefCoord,
                             qc::Element &TransformedEl, typename ConfiguratorType::VecType &TransformedLocalCoord ) const {

    typename ConfiguratorType::VecType coord, transformedCoord;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      coord[i] = El[i] + RefCoord[i] - 0.5 / this->_grid.H();
    }
    _rotationYaw.mult ( coord, transformedCoord );
    _rotationPitch.mult ( transformedCoord, coord );
    _rotationRoll.mult ( coord, transformedCoord );
    if ( AllowScaling && NumberScalingFactorIsOne )
      transformedCoord *= _scaling;
    else if ( AllowScaling && (NumberScalingFactorIsOne == false) )
      for ( int i = 0; i < ConfiguratorType::Dim; ++i )
        transformedCoord[i] *= _scalingVec[i];
    return qc::transformCoord<ConfiguratorType, ClipCoord> ( this->_grid, transformedCoord, _translation, TransformedEl, TransformedLocalCoord );
  }
  void evaluateDeformationOn01 ( const ParametricRigidBodyMotion3D<ConfiguratorType, AllowScaling, NumberScalingFactorIsOne> &/*ParDef*/, const typename ConfiguratorType::VecType &Position, typename ConfiguratorType::VecType &DeformedPosition ) const {
    typename ConfiguratorType::VecType tmp ( Position );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      tmp[i] -= 0.5;
    _rotationYaw.mult ( tmp, DeformedPosition );
    _rotationPitch.mult ( DeformedPosition, tmp );
    _rotationRoll.mult ( tmp, DeformedPosition );
    if ( AllowScaling && NumberScalingFactorIsOne ){
      DeformedPosition *= _scaling;
    }
    else if ( AllowScaling && (NumberScalingFactorIsOne == false) ){
      for ( int i = 0; i < ConfiguratorType::Dim; ++i )
        DeformedPosition[i] *= _scalingVec[i];
    }
    DeformedPosition += _translation;
  }
  
  void evaluateDerivativeDeformationOn01 ( const typename ConfiguratorType::VecType &Position,
                                          aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    typename ConfiguratorType::VecType temp;
    typename ConfiguratorType::VecType x ( Position );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      x[i] -= 0.5;
    // Yaw derivative
    {
      const RealType cosYaw = _rotationYaw[0][0];
      const RealType sinYaw = _rotationYaw[1][0];
      Jacobian[ConfiguratorType::Dim][0] = -sinYaw * x[0] - cosYaw * x[1];
      Jacobian[ConfiguratorType::Dim][1] =  cosYaw * x[0] - sinYaw * x[1];
      Jacobian[ConfiguratorType::Dim][2] =  0;
      _rotationPitch.mult ( Jacobian[ConfiguratorType::Dim], temp );
      _rotationRoll.mult ( temp, Jacobian[ConfiguratorType::Dim] );
      if ( AllowScaling && NumberScalingFactorIsOne )
        Jacobian[ConfiguratorType::Dim] *= _scaling;
      else if ( AllowScaling && (NumberScalingFactorIsOne == false) )
        for ( int i = 0; i < ConfiguratorType::Dim; ++i )
          Jacobian[ConfiguratorType::Dim][i] *= _scalingVec[i];
    }

    // Pitch derivative
    {
      const RealType cosPitch = _rotationPitch[0][0];
      const RealType sinPitch = _rotationPitch[2][0];
      _rotationYaw.mult ( x, Jacobian[ConfiguratorType::Dim+1] );
      temp[0] = -sinPitch * Jacobian[ConfiguratorType::Dim+1][0] - cosPitch * Jacobian[ConfiguratorType::Dim+1][2];//-sinPitch * x[0] - cosPitch * x[2];
      temp[1] = 0;
      temp[2] = cosPitch * Jacobian[ConfiguratorType::Dim+1][0] - sinPitch * Jacobian[ConfiguratorType::Dim+1][2];//cosPitch * x[0] - sinPitch * x[2];
      _rotationRoll.mult ( temp, Jacobian[ConfiguratorType::Dim+1] );
      if ( AllowScaling && NumberScalingFactorIsOne )
        Jacobian[ConfiguratorType::Dim+1] *= _scaling;
      else if ( AllowScaling && (NumberScalingFactorIsOne == false) )
        for ( int i = 0; i < ConfiguratorType::Dim; ++i )
          Jacobian[ConfiguratorType::Dim+1][i] *= _scalingVec[i];
    }

    // Roll derivative
    {
      const RealType cosRoll = _rotationRoll[1][1];
      const RealType sinRoll = _rotationRoll[2][1];
      _rotationYaw.mult ( x, Jacobian[ConfiguratorType::Dim+2] );
      _rotationPitch.mult ( Jacobian[ConfiguratorType::Dim+2], temp );
      Jacobian[ConfiguratorType::Dim+2][0] = 0;
      Jacobian[ConfiguratorType::Dim+2][1] = -sinRoll * temp[1] - cosRoll * temp[2];//-sinRoll * x[1] - cosRoll * x[2];
      Jacobian[ConfiguratorType::Dim+2][2] =  cosRoll * temp[1] - sinRoll * temp[2];//cosRoll * x[1] - sinRoll * x[2];
      if ( AllowScaling && NumberScalingFactorIsOne )
        Jacobian[ConfiguratorType::Dim+2] *= _scaling;
      else if ( AllowScaling && (NumberScalingFactorIsOne == false) )
        for ( int i = 0; i < ConfiguratorType::Dim; ++i )
          Jacobian[ConfiguratorType::Dim+2][i] *= _scalingVec[i];
    }

    Jacobian[0][0] = Jacobian[1][1] = Jacobian[2][2] = 1;
    Jacobian[1][0] = Jacobian[0][1] = Jacobian[2][0] = Jacobian[0][2] = Jacobian[2][1] = Jacobian[1][2] = 0;
    
    if ( AllowScaling && NumberScalingFactorIsOne ) {
      _rotationYaw.mult ( x, temp );
      _rotationPitch.mult ( temp, x );
      _rotationRoll.mult ( x, Jacobian[ConfiguratorType::Dim+3] );
    } else if ( AllowScaling && (NumberScalingFactorIsOne == false) ){
      _rotationYaw.mult ( x, temp );
      _rotationPitch.mult ( temp, x );
      _rotationRoll.mult ( x, temp);
      Jacobian[ConfiguratorType::Dim+3][0] = temp[0]; Jacobian[ConfiguratorType::Dim+3][1] = 0.; Jacobian[ConfiguratorType::Dim+3][2] = 0.;
      Jacobian[ConfiguratorType::Dim+4][0] = 0.; Jacobian[ConfiguratorType::Dim+4][1] = temp[1]; Jacobian[ConfiguratorType::Dim+4][2] = 0.;
      Jacobian[ConfiguratorType::Dim+5][0] = 0.; Jacobian[ConfiguratorType::Dim+5][1] = 0.; Jacobian[ConfiguratorType::Dim+5][2] = temp[2];
    }
  }
  
  void evaluateDerivativeDeformation ( const typename ConfiguratorType::ElementType &El,
                                      const typename ConfiguratorType::VecType &RefCoord,
                                      aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    typename ConfiguratorType::VecType x;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      x[i] = ( El[i] + RefCoord[i] ) * this->_grid.H();
    }
    evaluateDerivativeDeformationOn01 ( x, Jacobian );
  }
  
  void evaluateSpatialJacobian ( const typename ConfiguratorType::VecType &/*Position*/,
                                    aol::Mat< ConfiguratorType::Dim, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    aol::Matrix33<RealType> rotationMatrix;
    aol::Matrix33<RealType> temp;
    temp.makeProduct(_rotationRoll, _rotationPitch);
    rotationMatrix.makeProduct(temp,_rotationYaw);
    Jacobian = rotationMatrix;
    if ( AllowScaling && NumberScalingFactorIsOne )
      Jacobian *= _scaling;
    else if ( AllowScaling && (NumberScalingFactorIsOne == false) )
      for ( int i = 0; i < ConfiguratorType::Dim; ++i )
        Jacobian[i] *= _scalingVec[i];
  }

};

/** 
 * \author Berkels
 */
template <typename _ConfiguratorType>
class ParametricTranslation : public ParametricDeformationBase<_ConfiguratorType> {
public:
  typedef _ConfiguratorType ConfiguratorType;
  static const int NumberOfDeformParameters = ConfiguratorType::Dim;
  static const int NumOfDeformParametersComponents = 1;
  typedef typename ConfiguratorType::RealType RealType;
private:
  typename ConfiguratorType::VecType _translation;
public:
  ParametricTranslation ( const typename ConfiguratorType::InitType &Initializer, const aol::MultiVector<RealType> &DeformParameters )
    : ParametricDeformationBase<ConfiguratorType> ( Initializer ) {
    DeformParameters.copyTo ( _translation );
  }

  ParametricTranslation ( const typename ConfiguratorType::InitType &Initializer, const aol::Vector<RealType> &DeformParametersVec )
    : ParametricDeformationBase<ConfiguratorType> ( Initializer ) {
    setDeformationParametersAsVec ( DeformParametersVec );
  }

  ParametricTranslation ( const typename ConfiguratorType::InitType &Initializer, const aol::Vec<ConfiguratorType::Dim, RealType> &Translation )
    : ParametricDeformationBase<ConfiguratorType> ( Initializer ) {
    _translation = Translation;
  }

  ParametricTranslation ( const typename ConfiguratorType::InitType &Initializer )
    : ParametricDeformationBase<ConfiguratorType> ( Initializer ) {
    _translation.setZero();
  }

  void setTranslation ( const aol::Vec<ConfiguratorType::Dim, RealType> &Translation ) {
    _translation = Translation;
  }

  void setDeformationParametersAsVec ( const aol::Vector<RealType> &DeformParametersVec ) {
    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      _translation[i] = DeformParametersVec[i];
  }

  static void setIdentityDeformationParameters ( aol::MultiVector<RealType> &DeformParameters ) {
    DeformParameters.setZero();
  }

  static void setIdentityDeformationParametersAsVec ( aol::Vector<RealType> &DeformParametersAsVec ) {
    DeformParametersAsVec.setZero();
  }

  static void concatenateDeformationParameters ( const aol::MultiVector<RealType> &ArgDeformPars, aol::MultiVector<RealType> &DestDeformPars ) {
    DestDeformPars += ArgDeformPars;
  }

  static void getDeformationParametersOfInverse ( const aol::MultiVector<RealType> &DeformPars, aol::MultiVector<RealType> &InverseDeformPars )  {
    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      InverseDeformPars[0][i] = -DeformPars[0][i];
  }

  static aol::Vec<1, int> getDeformParametersSize ( ) {
    return aol::Vec<1, int> ( ConfiguratorType::Dim );
  }

  template <bool ClipCoord>
  bool evaluateDeformation ( const typename ConfiguratorType::ElementType &El,
                             const typename ConfiguratorType::VecType &RefCoord,
                             qc::Element &TransformedEl, typename ConfiguratorType::VecType &TransformedLocalCoord ) const {
    return qc::transformCoord<ConfiguratorType, ClipCoord> ( this->_grid, El, RefCoord, _translation, TransformedEl, TransformedLocalCoord );
  }
  void evaluateDeformationOn01 ( const aol::Vec<ConfiguratorType::Dim, RealType> &Position, typename aol::Vec<ConfiguratorType::Dim, RealType> &DeformedPosition ) const {
    DeformedPosition = Position;
    DeformedPosition += _translation;
  }
  void evaluateDeformationOn01 ( const ParametricTranslation<ConfiguratorType> &/*ParDef*/, const typename ConfiguratorType::VecType &Position, typename ConfiguratorType::VecType &DeformedPosition ) const {
    evaluateDeformationOn01 ( Position, DeformedPosition );
  }
  void evaluateDeformationOn0N ( typename aol::Vec<ConfiguratorType::Dim, RealType> &DeformedPosition ) const {
    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      DeformedPosition[i] += _translation[i] / this->_grid.H();
  }
  void evaluateDerivativeDeformationOn01 ( const typename ConfiguratorType::VecType &/*Position*/,
                                           aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    Jacobian.setIdentity();
  }
  void evaluateDerivativeDeformation ( const typename ConfiguratorType::ElementType &/*El*/,
                                       const typename ConfiguratorType::VecType &/*RefCoord*/,
                                       aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    Jacobian.setIdentity();
  }
  void evaluateDerivativeDeformationOn0N ( const typename ConfiguratorType::VecType &/*Position*/,
                                           aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    Jacobian.setZero();
    Jacobian.addToDiagonal ( aol::ZOTrait<RealType>::one / this->_grid.H() );
  }

  void evaluateSecondDerivativeDeformationOn0N ( const typename ConfiguratorType::VecType &/*Position*/,
                                                 aol::Tensor3rdOrder< NumberOfDeformParameters, NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Hessian ) const {
    Hessian.setZero();
  }
  void evaluateSpatialJacobianOn01 ( const typename ConfiguratorType::VecType &/*Position*/,
                                     aol::Mat< ConfiguratorType::Dim, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    Jacobian.setIdentity();
  }
  void evaluateSpatialJacobianOn0N ( const typename ConfiguratorType::VecType &/*Position*/,
                                     aol::Mat< ConfiguratorType::Dim, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    Jacobian.setIdentity();
  }
  void evaluateMixedHessianOn01 ( const typename ConfiguratorType::VecType &/*Position*/,
                                  aol::Tensor3rdOrder< NumberOfDeformParameters, ConfiguratorType::Dim, ConfiguratorType::Dim, RealType > &Hessian ) const {
    Hessian.setZero();
  }
  void evaluateMixedHessianOn0N ( const typename ConfiguratorType::VecType &/*Position*/,
                                  aol::Tensor3rdOrder< NumberOfDeformParameters, ConfiguratorType::Dim, ConfiguratorType::Dim, RealType > &Hessian ) const {
    Hessian.setZero();
  }
  void evaluateSpatialHessianOn0N ( const typename ConfiguratorType::VecType &/*Position*/,
                                    aol::Tensor3rdOrder< ConfiguratorType::Dim, ConfiguratorType::Dim, ConfiguratorType::Dim, RealType > &Hessian ) const {
    Hessian.setZero();
  }
};

/**
 * \author Berkels
 */
template <typename _ConfiguratorType>
class ParametricIdentity : public ParametricDeformationBase<_ConfiguratorType> {
public:
  typedef _ConfiguratorType ConfiguratorType;
  static const int NumberOfDeformParameters = 0;
  static const int NumOfDeformParametersComponents = 0;
  typedef typename ConfiguratorType::RealType RealType;
private:
  typename ConfiguratorType::VecType _translation;
public:
  ParametricIdentity ( const typename ConfiguratorType::InitType &Initializer, const aol::MultiVector<RealType> &/*DeformParameters*/ )
    : ParametricDeformationBase<ConfiguratorType> ( Initializer ) { }

  ParametricIdentity ( const typename ConfiguratorType::InitType &Initializer, const aol::Vector<RealType> &/*DeformParametersVec*/ )
    : ParametricDeformationBase<ConfiguratorType> ( Initializer ) { }

  ParametricIdentity ( const typename ConfiguratorType::InitType &Initializer )
    : ParametricDeformationBase<ConfiguratorType> ( Initializer ) { }

  void setDeformationParametersAsVec ( const aol::Vector<RealType> &/*DeformParametersVec*/ ) { }

  static void setIdentityDeformationParameters ( aol::MultiVector<RealType> &/*DeformParameters*/ ) { }

  static void setIdentityDeformationParametersAsVec ( aol::Vector<RealType> &/*DeformParametersAsVec*/ ) { }

  static void concatenateDeformationParameters ( const aol::MultiVector<RealType> &/*ArgDeformPars*/, aol::MultiVector<RealType> &/*DestDeformPars*/ ) { }
  static void getDeformationParametersOfInverse ( const aol::MultiVector<RealType> &/*DeformPars*/, aol::MultiVector<RealType> &/*InverseDeformPars*/ )  { }
  static aol::Vec<0, int> getDeformParametersSize ( ) {
    return aol::Vec<0, int> ();
  }

  template <bool ClipCoord>
  bool evaluateDeformation ( const typename ConfiguratorType::ElementType &El,
                             const typename ConfiguratorType::VecType &RefCoord,
                             qc::Element &TransformedEl, typename ConfiguratorType::VecType &TransformedLocalCoord ) const {
    TransformedEl = El;
    TransformedLocalCoord = RefCoord;
    return true;
  }

  void evaluateDeformationOn01 ( const aol::Vec<ConfiguratorType::Dim, RealType> &Position, typename aol::Vec<ConfiguratorType::Dim, RealType> &DeformedPosition ) const {
      DeformedPosition = Position;
  }

  void evaluateDeformationOn01 ( const ParametricIdentity<ConfiguratorType> &/*ParDef*/, const typename ConfiguratorType::VecType &Position, typename ConfiguratorType::VecType &DeformedPosition ) const {
      evaluateDeformationOn01 ( Position, DeformedPosition );
  }

  void evaluateDeformationOn0N ( typename aol::Vec<ConfiguratorType::Dim, RealType> &/*DeformedPosition*/ ) const { }

  void evaluateDerivativeDeformationOn01 ( const typename ConfiguratorType::VecType &/*Position*/,
                                           aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &/*Jacobian*/ ) const {  }

  void evaluateDerivativeDeformation ( const typename ConfiguratorType::ElementType &/*El*/,
                                       const typename ConfiguratorType::VecType &/*RefCoord*/,
                                       aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &/*Jacobian*/ ) const { }

  void evaluateDerivativeDeformationOn0N ( const typename ConfiguratorType::VecType &/*Position*/,
                                           aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &/*Jacobian*/ ) const { }

  void evaluateSecondDerivativeDeformationOn0N ( const typename ConfiguratorType::VecType &/*Position*/,
                                                 aol::Tensor3rdOrder< NumberOfDeformParameters, NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Hessian ) const {
    Hessian.setZero();
  }

  void evaluateSpatialJacobianOn01 ( const typename ConfiguratorType::VecType &/*Position*/,
                                     aol::Mat< ConfiguratorType::Dim, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    Jacobian.setIdentity();
  }

  void evaluateSpatialJacobianOn0N ( const typename ConfiguratorType::VecType &/*Position*/,
                                     aol::Mat< ConfiguratorType::Dim, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    Jacobian.setIdentity();
  }

  void evaluateMixedHessianOn01 ( const typename ConfiguratorType::VecType &/*Position*/,
                                  aol::Tensor3rdOrder< NumberOfDeformParameters, ConfiguratorType::Dim, ConfiguratorType::Dim, RealType > &Hessian ) const {
    Hessian.setZero();
  }

  void evaluateMixedHessianOn0N ( const typename ConfiguratorType::VecType &/*Position*/,
                                   aol::Tensor3rdOrder< NumberOfDeformParameters, ConfiguratorType::Dim, ConfiguratorType::Dim, RealType > &Hessian ) const {
    Hessian.setZero();
  }

  void evaluateSpatialHessianOn0N ( const typename ConfiguratorType::VecType &/*Position*/,
                                    aol::Tensor3rdOrder< ConfiguratorType::Dim, ConfiguratorType::Dim, ConfiguratorType::Dim, RealType > &Hessian ) const {
    Hessian.setZero();
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
inline bool parametricTransformAndClipCoord ( const ConfiguratorType &Config,
                                              const typename ConfiguratorType::ElementType &El,
                                              const typename ConfiguratorType::DomVecType &RefCoord,
                                              const ParametricDeformationType &ParDef,
                                              typename ConfiguratorType::ElementType &TransformedEl,
                                              typename ConfiguratorType::DomVecType &TransformedLocalCoord ) {
  bool insideFlag = true;
  typename ConfiguratorType::VecType position, coord;
  Config.getGlobalCoords( El, RefCoord, position );
  ParDef.evaluateDeformationOn01 ( ParDef, position, coord );
  for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
    if ( coord[i] < 0. ) {
      coord[i] = 0.;
      insideFlag = false;
    } else if ( coord[i] >= 1. ) {
      coord[i] = 1. - std::numeric_limits<typename ConfiguratorType::RealType>::epsilon(); // largest number < 1.0 (for clipping to [0,1[)
      insideFlag = false;
    }
  }
  Config.getLocalCoords( coord, TransformedEl, TransformedLocalCoord );
  return insideFlag;
}

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
inline void parametricTransformAndClipCoord ( const ConfiguratorType &Config,
                                              const typename ConfiguratorType::ElementType &El,
                                              const typename ConfiguratorType::DomVecType &RefCoord,
                                              const ParametricDeformationType &ParDef,
                                              typename ConfiguratorType::ElementType &TransformedEl,
                                              typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                                              aol::Vec<ConfiguratorType::Dim,bool> &CoordinateWithinLimits ) {
  CoordinateWithinLimits.setAll( true );
  typename ConfiguratorType::VecType position, coord;
  Config.getGlobalCoords( El, RefCoord, position );
  ParDef.evaluateDeformationOn01 ( ParDef, position, coord );
  for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
    if ( coord[i] < 0. ) {
      coord[i] = 0.;
      CoordinateWithinLimits[i] = false;
    } else if ( coord[i] >= 1. ) {
      coord[i] = 1. - std::numeric_limits<typename ConfiguratorType::RealType>::epsilon(); // largest number < 1.0 (for clipping to [0,1[)
      CoordinateWithinLimits[i] = false;
    }
  }
  Config.getLocalCoords( coord, TransformedEl, TransformedLocalCoord );
}

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricDeformedPosition {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::VecType _pos01;
  typename ConfiguratorType::VecType _deformedPos;
  bool _transformPositionInDomain;
public:
  ParametricDeformedPosition ( const int X, const int Y,
                              const ParametricDeformationType &ParDef,
                              const aol::Vec3<int> &GridSize,
                              const RealType H )
  : _pos01 ( X*H, Y*H ),
  _transformPositionInDomain ( true ) {

    ParDef.evaluateDeformationOn01 ( ParDef, _pos01, _deformedPos );
    _deformedPos /= H;
    for ( int c = 0; c < ConfiguratorType::Dim; c++ ) {
      const RealType temp = aol::Clamp ( _deformedPos[c], aol::ZOTrait<RealType>::zero, static_cast<RealType> ( GridSize[c] - 1 ) );
      if ( temp != _deformedPos[c] )
        _transformPositionInDomain = false;

      _deformedPos[c] = temp;
    }
  }

  const typename ConfiguratorType::VecType &getUndefPos01 ( ) const {
    return _pos01;
  }

  const typename ConfiguratorType::VecType &getDefPos ( ) const {
    return _deformedPos;
  }

  bool isDefPosInDomain ( ) const {
    return _transformPositionInDomain;
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType, qc::Dimension Dim>
struct doFastParametricDeformImage {
static void apply ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                    const typename ConfiguratorType::InitType &Grid,
                    aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                    const ParametricDeformationType &ParDef );
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
struct doFastParametricDeformImage<ConfiguratorType, ParametricDeformationType, qc::QC_2D> {
  static void apply ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                      const typename ConfiguratorType::InitType &Grid,
                      aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                      const ParametricDeformationType &ParDef ) {
    typedef typename ConfiguratorType::RealType RealType;
    const typename ConfiguratorType::ArrayType imageArray ( Image, Grid, aol::FLAT_COPY );

    const aol::Vec3<int> gridSize = Grid.getSize();
    const RealType h = static_cast<RealType> ( Grid.H() );

    for ( int j = 0; j < gridSize[1]; ++j ) {
      for ( int i = 0; i < gridSize[0]; ++i ) {
        const int index = qc::ILexCombine2( i, j, Grid.getNumX() );
        qc::ParametricDeformedPosition<ConfiguratorType, ParametricDeformationType> defPos ( i, j, ParDef, gridSize, h );
        DeformedImage[index] = imageArray.interpolate ( defPos.getDefPos() );
      }
    }
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
struct doFastParametricDeformImage<ConfiguratorType, ParametricDeformationType, qc::QC_3D> {
  static void apply ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                      const typename ConfiguratorType::InitType &Grid,
                      aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                      const ParametricDeformationType &ParDef ) {
    typedef typename ConfiguratorType::RealType RealType;
    const typename ConfiguratorType::ArrayType imageArray ( Image, Grid, aol::FLAT_COPY );

    const aol::Vec3<int> gridSize = Grid.getSize();
    const RealType h = static_cast<RealType> ( Grid.H() );

    for ( int k = 0; k < gridSize[2]; ++k ) {
      for ( int j = 0; j < gridSize[1]; ++j ) {
        for ( int i = 0; i < gridSize[0]; ++i ) {
          const int index = qc::ILexCombine3( i, j, k, Grid.getNumX(), Grid.getNumY() );
          const typename ConfiguratorType::VecType pos ( i*h, j*h, k*h );
          typename ConfiguratorType::VecType ds;
          ParDef.evaluateDeformationOn01 ( ParDef, pos, ds );
          for ( int c = 0; c < ConfiguratorType::Dim; c++ )
            ds[c] = aol::Clamp ( ds[c] / h, aol::ZOTrait<RealType>::zero, static_cast<RealType> ( gridSize[c] - 1 ) );
          DeformedImage[index] = imageArray.interpolate ( ds );
        }
      }
    }
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType, qc::Dimension Dim>
void FastParametricDeformImage ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                                 const typename ConfiguratorType::InitType &Grid,
                                 aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                                 const ParametricDeformationType &ParDef ) {
  doFastParametricDeformImage<ConfiguratorType, ParametricDeformationType, Dim>::apply ( Image, Grid, DeformedImage, ParDef );
}

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType, qc::Dimension Dim>
void FastParametricInvDeformImage ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                                    const typename ConfiguratorType::InitType &Grid,
                                    aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                                    const aol::MultiVector<typename ConfiguratorType::RealType> &ParDefParam ) {
  aol::MultiVector<typename ConfiguratorType::RealType> inverseParDefParam ( ParametricDeformationType::getDeformParametersSize() );
  ParametricDeformationType::getDeformationParametersOfInverse ( ParDefParam, inverseParDefParam );
  ParametricDeformationType inverseParDef ( Grid, inverseParDefParam );
  qc::FastParametricDeformImage<ConfiguratorType, ParametricDeformationType, Dim> ( Image, Grid, DeformedImage, inverseParDef );
}

/**
 * \f$ \frac{1}{2}\int ((T\circ\Phi)(x)-R(x))^2 dx \f$ where \f$\Phi\f$ is a parametric deformation
 * whose parameters are given in the constructor and T is given as argument in apply(Add).
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricSSDEnergyFromImage : public aol::FENonlinIntegrationScalarInterface<ConfiguratorType, ParametricSSDEnergyFromImage<ConfiguratorType, ParametricDeformationType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const ParametricDeformationType &_parDef;
public:

  ParametricSSDEnergyFromImage ( const typename ConfiguratorType::InitType &Grid, 
                                 const aol::Vector<RealType> &ImR,
                                 const ParametricDeformationType &ParDef )
  : aol::FENonlinIntegrationScalarInterface<ConfiguratorType, ParametricSSDEnergyFromImage<ConfiguratorType, ParametricDeformationType> > ( Grid ),
    _r( Grid, ImR ),
    _parDef ( ParDef ) {}

  RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord ) const {

    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;
    if ( !_parDef.template evaluateDeformation<true> ( El, RefCoord, transformedEl, transformedLocalCoord) ) {
      return 0;
    }

    return 0.5 * aol::Sqr ( DiscFuncs.evaluate(transformedEl, transformedLocalCoord) - _r.evaluateAtQuadPoint(El, QuadPoint) );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class VariationOfParametricSSDEnergyWRTParamsFromImage
  : public aol::VectorFENonlinIntegrationVectorInterface < ConfiguratorType,
                                                           VariationOfParametricSSDEnergyWRTParamsFromImage<ConfiguratorType, ParametricDeformationType>,
                                                           1, ParametricDeformationType::NumberOfDeformParameters > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
private:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const ParametricDeformationType &_parDef;

public:
  VariationOfParametricSSDEnergyWRTParamsFromImage ( const typename ConfiguratorType::InitType &Initializer,
                                                     const aol::Vector<RealType> &ImR,
                                                     const ParametricDeformationType &ParDef )
    : aol::VectorFENonlinIntegrationVectorInterface < ConfiguratorType,
                                                      VariationOfParametricSSDEnergyWRTParamsFromImage<ConfiguratorType, ParametricDeformationType>,
                                                      1, ParametricDeformationType::NumberOfDeformParameters > ( Initializer ),
      _r( Initializer, ImR ),
      _parDef ( ParDef ) {}

  void evaluateIntegrand ( const aol::auto_container<1, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrFuncs, // DiscrFuncs[0] is supposed to contain T.
                           const typename ConfiguratorType::ElementType &El,
                           int QuadPoint, const typename ConfiguratorType::VecType &RefCoord, aol::Vec<ParametricDeformationType::NumberOfDeformParameters, RealType> &Integrand ) const {

    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;
    if ( !_parDef.template evaluateDeformation<true> ( El, RefCoord, transformedEl, transformedLocalCoord) ) {
      Integrand.setZero();
      return;
    }

    aol::Mat< ParametricDeformationType::NumberOfDeformParameters, ConfiguratorType::Dim, RealType > jacobian;
    _parDef.evaluateDerivativeDeformation ( El, RefCoord, jacobian );

    typename ConfiguratorType::VecType gradT;
    DiscrFuncs[0].evaluateGradient ( transformedEl, transformedLocalCoord, gradT );

    jacobian.mult ( gradT, Integrand );
    Integrand *= ( DiscrFuncs[0].evaluate ( transformedEl, transformedLocalCoord ) - _r.evaluateAtQuadPoint ( El, QuadPoint ) );
  }
};

/**
 * \brief  \f$ \frac{1}{2}\int ((T\circ\Phi)(x)-R(x))^2 dx \f$ where \f$\Phi\f$ is a parametric deformation
 * whose parameters are given as argument in apply(Add).
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricSSDEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> _r;
  const aol::Vector<RealType> _t;
public:
  ParametricSSDEnergy ( const typename ConfiguratorType::InitType &Grid, 
                        const aol::Vector<RealType> &ImR,
                        const aol::Vector<RealType> &ImT )
  : _grid ( Grid ),
    _r( ImR ), 
    _t( ImT ) {}

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    const ParametricDeformationType parDef ( _grid, MArg );
    ParametricSSDEnergyFromImage<ConfiguratorType, ParametricDeformationType> E ( _grid, _r, parDef );
    E.applyAdd ( _t, Dest );
  }

  void applyDerivative ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    const ParametricDeformationType parDef ( _grid, MArg );

    aol::Vector<RealType> dest ( ParametricDeformationType::NumberOfDeformParameters );
    VariationOfParametricSSDEnergyWRTParamsFromImage<ConfiguratorType, ParametricDeformationType> E ( _grid, _r, parDef );
    E.applyAdd ( _t, dest );
    MDest.copySplitFrom ( dest );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricCrossCorrelationForce
  : public aol::VectorFENonlinIntegrationVectorInterface < ConfiguratorType,
                                                           ParametricCrossCorrelationForce<ConfiguratorType, ParametricDeformationType>,
                                                           1, ParametricDeformationType::NumberOfDeformParameters > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const ParametricDeformationType &_parDef;

  ParametricCrossCorrelationForce ( const typename ConfiguratorType::InitType &Initializer,
                                    const aol::Vector<RealType> &ImR,
                                    const ParametricDeformationType &ParDef )
    : aol::VectorFENonlinIntegrationVectorInterface < ConfiguratorType,
                                                      ParametricCrossCorrelationForce<ConfiguratorType, ParametricDeformationType>,
                                                      1, ParametricDeformationType::NumberOfDeformParameters > ( Initializer ),
      _r( Initializer, ImR ),
      _parDef ( ParDef ) {}

  void evaluateIntegrand ( const aol::auto_container<1, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrFuncs, // DiscrFuncs[0] is supposed to contain T.
                           const typename ConfiguratorType::ElementType &El,
                           int QuadPoint, const typename ConfiguratorType::VecType &RefCoord, aol::Vec<ParametricDeformationType::NumberOfDeformParameters, RealType> &Integrand ) const {

    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;
    if ( !_parDef.template evaluateDeformation<true> ( El, RefCoord, transformedEl, transformedLocalCoord) ) {
      Integrand.setZero();
      return;
    }

    aol::Mat< ParametricDeformationType::NumberOfDeformParameters, ConfiguratorType::Dim, RealType > jacobian;
    _parDef.evaluateDerivativeDeformation ( El, RefCoord, jacobian );

    typename ConfiguratorType::VecType gradT;
    DiscrFuncs[0].evaluateGradient ( transformedEl, transformedLocalCoord, gradT );

    jacobian.mult ( gradT, Integrand );
    Integrand *= - _r.evaluateAtQuadPoint ( El, QuadPoint );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricNCCEnergy : public qc::NormalizedCrossCorrelationEnergy<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;

  ParametricNCCEnergy ( const typename ConfiguratorType::InitType &Grid, 
                        const aol::Vector<RealType> &ImR,
                        const aol::Vector<RealType> &ImT )
  : qc::NormalizedCrossCorrelationEnergy<ConfiguratorType> ( Grid, ImR, ImT ) {}

  // By overloading deformImage, we can use apply(Add) from the base class to calculate the energy.
  void deformImage ( const typename ConfiguratorType::ArrayType &Image, typename ConfiguratorType::ArrayType &DeformedImage , const aol::MultiVector<RealType> &Phi) const {
    ParametricDeformationType parDef ( this->_grid, Phi );
    qc::FastParametricDeformImage<ConfiguratorType, ParametricDeformationType, ConfiguratorType::Dim> ( Image, this->_grid, DeformedImage, parDef );
  }

  void applyDerivative ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    aol::Scalar<RealType> energy;
    this->apply ( MArg, energy );
    aol::Vector<RealType> temp ( this->_normalizedR );
    temp.addMultiple ( this->_lastNormalizedDeformedT, energy[0] );
    temp /= this->_varOfLastDeformedT;

    aol::Vector<RealType> destVec ( ParametricDeformationType::NumberOfDeformParameters );
    ParametricDeformationType parDef ( this->_grid, MArg );
    ParametricCrossCorrelationForce<ConfiguratorType, ParametricDeformationType> force ( this->_grid, temp, parDef );
    force.apply ( this->_t, destVec );
    MDest.copySplitFrom ( destVec );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class DeformationMatrix : protected aol::SparseMatrix<typename ConfiguratorType::RealType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const typename ConfiguratorType::InitType &_grid;
  const ParametricDeformationType _parDef;
public:
  DeformationMatrix ( const typename ConfiguratorType::InitType &Grid,
                     const aol::MultiVector<RealType> &DeformParams )
  : aol::SparseMatrix<RealType> ( Grid ),
  _grid ( Grid ),
  _parDef ( Grid, DeformParams ) {

    const aol::Vec3<int> gridSize = _grid.getSize();
    const RealType h = static_cast<RealType> ( _grid.H() );

    for ( int j = 0; j < gridSize[1]; ++j ) {
      for ( int i = 0; i < gridSize[0]; ++i ) {
        const int index = qc::ILexCombine2( i, j, _grid.getNumX() );
        qc::ParametricDeformedPosition<ConfiguratorType, ParametricDeformationType> defPos ( i, j, _parDef, gridSize, h );

        RealType X = defPos.getDefPos()[0];
        RealType Y = defPos.getDefPos()[1];
        int xL = static_cast<int> ( X ), yL = static_cast<int> ( Y );

        if ( X >= static_cast<RealType> ( gridSize[0] - 1 ) ) xL = gridSize[0] - 2;
        if ( Y >= static_cast<RealType> ( gridSize[1] - 1 ) ) yL = gridSize[1] - 2;
        if ( X < 0. ) xL = 0;
        if ( Y < 0. ) yL = 0;

        X -= xL;
        Y -= yL;

        const int xU = xL + 1, yU = yL + 1;

        this->set ( index, qc::ILexCombine2( xL, yL, _grid.getNumX() ), ( 1 - X ) * ( 1 - Y ) );
        this->set ( index, qc::ILexCombine2( xU, yL, _grid.getNumX() ), X * ( 1 - Y ) );
        this->set ( index, qc::ILexCombine2( xL, yU, _grid.getNumX() ), ( 1 - X ) * Y );
        this->set ( index, qc::ILexCombine2( xU, yU, _grid.getNumX() ), X * Y );
      }
    }
  }

  const aol::SparseMatrix<typename ConfiguratorType::RealType> &getSparseMatrixRef ( ) const {
    return *this;
  }

  using aol::SparseMatrix<RealType>::apply;
  using aol::SparseMatrix<RealType>::applyAddTransposed;
  using aol::SparseMatrix<RealType>::transposeTo;
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType, typename Imp, typename MatrixType>
class FELeastSquaresFunctionalDeformationDerivativeInterface : public aol::FEOpInterface<ConfiguratorType, aol::Vector<typename ConfiguratorType::RealType>, MatrixType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_grid;
public:
  explicit FELeastSquaresFunctionalDeformationDerivativeInterface ( const typename ConfiguratorType::InitType &Initializer )
    : aol::FEOpInterface<ConfiguratorType, aol::Vector<RealType>, MatrixType> (Initializer),
      _grid ( Initializer ) { }

  FELeastSquaresFunctionalDeformationDerivativeInterface ( const ConfiguratorType & Config,
                                                           const typename ConfiguratorType::InitType &Initializer )
    : aol::FEOpInterface<ConfiguratorType, aol::Vector<RealType>, MatrixType> (Config, Initializer),
      _grid ( Initializer ) {}

  virtual ~FELeastSquaresFunctionalDeformationDerivativeInterface( ) {}

  void apply ( const aol::Vector<RealType> &Arg, MatrixType &Dest ) const {
    const ParametricDeformationType parDef ( _grid, Arg );

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    aol::Mat<ParametricDeformationType::NumberOfDeformParameters, ConfiguratorType::Dim, RealType> jacobian;
    aol::Vec<ConfiguratorType::Dim, RealType> spatGrad;
    aol::Vec<ParametricDeformationType::NumberOfDeformParameters, RealType> parGrad;

    const typename IteratorType::EndType end = this->getConfigurator().end();
    int i = 0;
    for ( IteratorType it = this->getConfigurator().begin(); it != end; ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );
      const RealType vol = this->getConfigurator().vol ( *it );
      for ( int q = 0; q < numQuadPoints; ++q ) {
        parDef.evaluateDerivativeDeformation ( *it, bfs.getRefCoord ( q ), jacobian );
        this->asImp().evaluateIntegrandDerivative( parDef, *it, q, bfs.getRefCoord ( q ), spatGrad );
        jacobian.mult ( spatGrad, parGrad );
        for ( int c = 0; c < ParametricDeformationType::NumberOfDeformParameters; ++c ) {
          Dest.set ( i+q, c, parGrad[c] * sqrt( bfs.getWeight ( q ) * vol ) );
        }
      }
      i += numQuadPoints;
    }
  }

  void applyAdd ( const aol::Vector<RealType> &, MatrixType & ) const {
    throw aol::UnimplementedCodeException ( "applyAdd not implemented...", __FILE__, __LINE__ );
  }

  //! interface function, has to be provided in derived classes.
  void evaluateIntegrandDerivative ( const ParametricDeformationType &ParDef,
                                     const typename ConfiguratorType::ElementType &El,
                                     int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord,
                                     aol::Vec<ConfiguratorType::Dim, RealType> &Derivative ) const {
    throw aol::Exception ( "called the interface function", __FILE__, __LINE__ );
    return this->asImp().evaluateIntegrandDerivative ( ParDef, El, QuadPoint, RefCoord, Derivative );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricSSDGNFuncFromImage : public aol::FELeastSquaresFunctionalInterface<ConfiguratorType, ParametricSSDGNFuncFromImage<ConfiguratorType, ParametricDeformationType>, aol::FullMatrix<typename ConfiguratorType::RealType>, true> {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const ParametricDeformationType &_parDef;
public:
  ParametricSSDGNFuncFromImage ( const typename ConfiguratorType::InitType &Grid,
                                 const aol::Vector<RealType> &ImR,
                                 const ParametricDeformationType &ParDef )
  : aol::FELeastSquaresFunctionalInterface<ConfiguratorType, ParametricSSDGNFuncFromImage<ConfiguratorType, ParametricDeformationType>, aol::FullMatrix<typename ConfiguratorType::RealType>, true> ( Grid ),
    _r( Grid, ImR ),
    _parDef ( ParDef ) {}

  RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord ) const {

    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;
    qc::parametricTransformAndClipCoord<ConfiguratorType, ParametricDeformationType> ( this->getConfigurator(), El, RefCoord, _parDef, transformedEl, transformedLocalCoord );

    return ( DiscFuncs.evaluate(transformedEl, transformedLocalCoord) - _r.evaluateAtQuadPoint(El, QuadPoint) );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class VarOfParametricSSDGNFunc : public FELeastSquaresFunctionalDeformationDerivativeInterface<ConfiguratorType, ParametricDeformationType, VarOfParametricSSDGNFunc<ConfiguratorType, ParametricDeformationType>, aol::FullMatrix<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _t;
public:

  VarOfParametricSSDGNFunc ( const typename ConfiguratorType::InitType &Grid,
                             const aol::Vector<RealType> &ImR,
                             const aol::Vector<RealType> &ImT )
    : FELeastSquaresFunctionalDeformationDerivativeInterface<ConfiguratorType, ParametricDeformationType, VarOfParametricSSDGNFunc<ConfiguratorType, ParametricDeformationType>, aol::FullMatrix<typename ConfiguratorType::RealType> > ( Grid ),
      _r( Grid, ImR ),
      _t ( Grid, ImT ) {}

  void evaluateIntegrandDerivative ( const ParametricDeformationType &ParDef,
                                     const typename ConfiguratorType::ElementType &El,
                                     int /*QuadPoint*/, const typename ConfiguratorType::DomVecType &RefCoord,
                                     aol::Vec<ConfiguratorType::Dim, RealType> &Derivative ) const {

    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;
    aol::Vec<ConfiguratorType::Dim, bool> coordinateWithinLimits;
    qc::parametricTransformAndClipCoord<ConfiguratorType, ParametricDeformationType> ( this->getConfigurator(), El, RefCoord, ParDef, transformedEl, transformedLocalCoord, coordinateWithinLimits );

    _t.evaluateGradient ( transformedEl, transformedLocalCoord, Derivative );

    for ( int i = 0; i < ConfiguratorType::Dim; ++i ) {
      if ( coordinateWithinLimits[i] == false )
        Derivative[i] = 0;
    }
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
  class ParametricSSDGNFunc : public aol::Op<aol::Vector<typename ConfiguratorType::RealType>, aol::Vector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename aol::FullMatrix<RealType> DerivativeType;
private:
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_r;
  const aol::Vector<RealType> &_t;
  const qc::VarOfParametricSSDGNFunc<ConfiguratorType, ParametricDeformationType> _DF;
public:
  ParametricSSDGNFunc ( const typename ConfiguratorType::InitType &Grid,
                        const aol::Vector<RealType> &ImR,
                        const aol::Vector<RealType> &ImT )
    : _grid ( Grid ),
      _r( ImR ),
      _t( ImT ),
      _DF ( Grid, ImR, ImT ) {}

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    const ParametricDeformationType parDef ( _grid, Arg );
    ParametricSSDGNFuncFromImage<ConfiguratorType, ParametricDeformationType> F ( _grid, _r, parDef );
    F.apply ( _t, Dest );
  }
    
  void applyAdd ( const aol::Vector<RealType> &, aol::Vector<RealType> & ) const {
    throw aol::UnimplementedCodeException ( "applyAdd not implemented...", __FILE__, __LINE__ );
  }

  void applyDerivative ( const aol::Vector<RealType> &Arg, DerivativeType &MatDest ) const {
    _DF.apply ( Arg, MatDest );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricNCCGNFuncFromImage : public aol::FELeastSquaresFunctionalInterface<ConfiguratorType, ParametricNCCGNFuncFromImage<ConfiguratorType, ParametricDeformationType>, aol::FullMatrix<typename ConfiguratorType::RealType>, true> {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _normalizedR;
  const RealType _meanOfDeformedT;
  const RealType _oneOverVarOfDeformedT;
  const ParametricDeformationType &_parDef;
public:
  ParametricNCCGNFuncFromImage ( const typename ConfiguratorType::InitType &Grid,
                                 const aol::Vector<RealType> &NormalizedR,
                                 const aol::Vec2<RealType> &MeanAndVarOfDeformedT,
                                 const ParametricDeformationType &ParDef )
    : aol::FELeastSquaresFunctionalInterface<ConfiguratorType, ParametricNCCGNFuncFromImage<ConfiguratorType, ParametricDeformationType>, aol::FullMatrix<typename ConfiguratorType::RealType>, true> ( Grid ),
      _normalizedR( Grid, NormalizedR ),
      _meanOfDeformedT ( MeanAndVarOfDeformedT[0] ),
      _oneOverVarOfDeformedT ( aol::Reciprocal ( MeanAndVarOfDeformedT[1] ) ),
      _parDef ( ParDef ) {}

  RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord ) const {

    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;
    qc::parametricTransformAndClipCoord<ConfiguratorType, ParametricDeformationType> ( this->getConfigurator(), El, RefCoord, _parDef, transformedEl, transformedLocalCoord );

    return ( _oneOverVarOfDeformedT *( DiscFuncs.evaluate(transformedEl, transformedLocalCoord) - _meanOfDeformedT ) - _normalizedR.evaluateAtQuadPoint(El, QuadPoint) );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricNCCGNForcePart1
  : public aol::VectorFENonlinIntegrationVectorInterface < ConfiguratorType,
                                                           ParametricNCCGNForcePart1<ConfiguratorType, ParametricDeformationType>,
                                                           1, ParametricDeformationType::NumberOfDeformParameters > {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const ParametricDeformationType &_parDef;
public:
  ParametricNCCGNForcePart1 ( const typename ConfiguratorType::InitType &Initializer,
                              const ParametricDeformationType &ParDef )
    : aol::VectorFENonlinIntegrationVectorInterface < ConfiguratorType,
                                                      ParametricNCCGNForcePart1<ConfiguratorType, ParametricDeformationType>,
                                                      1, ParametricDeformationType::NumberOfDeformParameters > ( Initializer ),
      _parDef ( ParDef ) {}

  void evaluateIntegrand ( const aol::auto_container<1, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrFuncs, // DiscrFuncs[0] is supposed to contain T.
                           const typename ConfiguratorType::ElementType &El,
                           int /*QuadPoint*/, const typename ConfiguratorType::VecType &RefCoord, aol::Vec<ParametricDeformationType::NumberOfDeformParameters, RealType> &Integrand ) const {

    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;
    aol::Vec<ConfiguratorType::Dim, bool> coordinateWithinLimits;
    qc::parametricTransformAndClipCoord<ConfiguratorType, ParametricDeformationType> ( this->_config, El, RefCoord, _parDef, transformedEl, transformedLocalCoord, coordinateWithinLimits );

    aol::Mat< ParametricDeformationType::NumberOfDeformParameters, ConfiguratorType::Dim, RealType > jacobian;
    _parDef.evaluateDerivativeDeformation ( El, RefCoord, jacobian );

    typename ConfiguratorType::VecType gradT;
    DiscrFuncs[0].evaluateGradient ( transformedEl, transformedLocalCoord, gradT );

    for ( int i = 0; i < ConfiguratorType::Dim; ++i ) {
      if ( coordinateWithinLimits[i] == false )
        gradT[i] = 0;
    }

    jacobian.mult ( gradT, Integrand );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricNCCGNForcePart2
  : public aol::VectorFENonlinIntegrationVectorInterface < ConfiguratorType,
                                                           ParametricNCCGNForcePart2<ConfiguratorType, ParametricDeformationType>,
                                                           1, ParametricDeformationType::NumberOfDeformParameters > {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const RealType _meanOfDeformedT;
  const RealType _oneOverVarOfDeformedT;
  const ParametricDeformationType &_parDef;
public:
  ParametricNCCGNForcePart2 ( const typename ConfiguratorType::InitType &Initializer,
                              const aol::Vec2<RealType> &MeanAndVarOfDeformedT,
                              const ParametricDeformationType &ParDef )
    : aol::VectorFENonlinIntegrationVectorInterface < ConfiguratorType,
                                                      ParametricNCCGNForcePart2<ConfiguratorType, ParametricDeformationType>,
                                                      1, ParametricDeformationType::NumberOfDeformParameters > ( Initializer ),
      _meanOfDeformedT ( MeanAndVarOfDeformedT[0] ),
      _oneOverVarOfDeformedT ( aol::Reciprocal ( MeanAndVarOfDeformedT[1] ) ),
      _parDef ( ParDef ) {}

  void evaluateIntegrand ( const aol::auto_container<1, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrFuncs, // DiscrFuncs[0] is supposed to contain T.
                           const typename ConfiguratorType::ElementType &El,
                           int /*QuadPoint*/, const typename ConfiguratorType::VecType &RefCoord, aol::Vec<ParametricDeformationType::NumberOfDeformParameters, RealType> &Integrand ) const {

    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;
    aol::Vec<ConfiguratorType::Dim, bool> coordinateWithinLimits;
    qc::parametricTransformAndClipCoord<ConfiguratorType, ParametricDeformationType> ( this->_config, El, RefCoord, _parDef, transformedEl, transformedLocalCoord, coordinateWithinLimits );

    aol::Mat< ParametricDeformationType::NumberOfDeformParameters, ConfiguratorType::Dim, RealType > jacobian;
    _parDef.evaluateDerivativeDeformation ( El, RefCoord, jacobian );

    typename ConfiguratorType::VecType gradT;
    DiscrFuncs[0].evaluateGradient ( transformedEl, transformedLocalCoord, gradT );

    for ( int i = 0; i < ConfiguratorType::Dim; ++i ) {
      if ( coordinateWithinLimits[i] == false )
        gradT[i] = 0;
    }

    jacobian.mult ( gradT, Integrand );
    Integrand *= _oneOverVarOfDeformedT * ( DiscrFuncs[0].evaluate ( transformedEl, transformedLocalCoord ) - _meanOfDeformedT );
  }
};

/**
 * \author Berkels
 *
 * \note Note fully tested yet.
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricNCCGNFunc : public aol::Op<aol::Vector<typename ConfiguratorType::RealType>, aol::Vector<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename aol::FullMatrix<RealType> DerivativeType;
private:
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_r;
  const aol::Vector<RealType> &_t;
  // This takes care of normalizing R and provides a MassOp.
  const qc::ParametricNCCEnergy<ConfiguratorType, ParametricDeformationType> _E;
public:
  ParametricNCCGNFunc ( const typename ConfiguratorType::InitType &Grid,
                        const aol::Vector<RealType> &ImR,
                        const aol::Vector<RealType> &ImT )
    : _grid ( Grid ),
      _r( ImR ),
      _t( ImT ),
      _E ( Grid, ImR, ImT ) {}

  aol::Vec2<RealType> computeMeanAndVarOfDeformedImage ( const aol::Vector<RealType> &Image, const ParametricDeformationType &ParDef ) const {

    aol::Vec2<RealType> meanAndVarOfDeformedImage;
    aol::Vector<RealType> imageCopy ( Image );

    qc::DataGenerator<ConfiguratorType> generator ( this->_grid );
    qc::MultiArray<RealType, ConfiguratorType::Dim> deformation ( this->_grid );
    generator.template generateDeformationFromParametricDeformation<ParametricDeformationType, true> ( ParDef, deformation );

    const qc::DeformMassOp<ConfiguratorType> MDefRight ( this->_grid, deformation, aol::ONTHEFLY, RIGHT );
    const qc::DeformMassOp<ConfiguratorType> MDefBoth ( this->_grid, deformation, aol::ONTHEFLY, BOTH );

    aol::Vector<RealType> temp ( imageCopy, aol::STRUCT_COPY );
    MDefRight.apply ( Image, temp );
    meanAndVarOfDeformedImage[0] = temp.sum();

    imageCopy.addToAll ( -meanAndVarOfDeformedImage[0] );
    MDefBoth.apply ( imageCopy, temp );
    meanAndVarOfDeformedImage[1] = sqrt ( imageCopy * temp );
    return meanAndVarOfDeformedImage;
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    const ParametricDeformationType parDef ( _grid, Arg );
    qc::ParametricNCCGNFuncFromImage<ConfiguratorType, ParametricDeformationType> F ( _grid, _E.getNormalizedRRef(), computeMeanAndVarOfDeformedImage ( _t, parDef ), parDef );
    F.apply ( _t, Dest );
  }

  void applyAdd ( const aol::Vector<RealType> &, aol::Vector<RealType> & ) const {
    throw aol::UnimplementedCodeException ( "applyAdd not implemented...", __FILE__, __LINE__ );
  }
    
  void applyDerivative ( const aol::Vector<RealType> &Arg, DerivativeType &MatDest ) const {
    const ParametricDeformationType parDef ( _grid, Arg );
    const aol::Vec2<RealType> meanAndVarOfDeformedT ( computeMeanAndVarOfDeformedImage ( _t, parDef ) );
    const qc::VarOfParametricSSDGNFunc<ConfiguratorType, ParametricDeformationType> DF ( _grid, _r, _t );
    const RealType oneOverVarOfDeformedT = aol::Reciprocal ( meanAndVarOfDeformedT[1] );
    DF.apply ( Arg, MatDest );
    MatDest *= oneOverVarOfDeformedT;

    aol::Vector<RealType> defNormedTUndef ( _t );
    defNormedTUndef.addToAll ( -meanAndVarOfDeformedT[0] );
    defNormedTUndef *= oneOverVarOfDeformedT;

    const aol::DiscreteFunctionDefault<ConfiguratorType> discrDefNormedTUndef ( DF.getConfigurator(), defNormedTUndef );

    const qc::ParametricNCCGNForcePart1<ConfiguratorType, ParametricDeformationType> dfPart1 ( _grid, parDef );
    aol::Vector<RealType> part1 ( ParametricDeformationType::NumberOfDeformParameters );
    dfPart1.apply ( _t, part1 );

    const qc::ParametricNCCGNForcePart2<ConfiguratorType, ParametricDeformationType> dfPart2 ( _grid, meanAndVarOfDeformedT, parDef );
    aol::Vector<RealType> part2 ( ParametricDeformationType::NumberOfDeformParameters );
    dfPart2.apply ( _t, part2 );

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;

    const typename IteratorType::EndType end = DF.getConfigurator().end();
    int i = 0;
    for ( IteratorType it = DF.getConfigurator().begin(); it != end; ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = DF.getConfigurator().getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );
      const RealType vol = DF.getConfigurator().vol ( *it );
      for ( int q = 0; q < numQuadPoints; ++q ) {
        const RealType weight = sqrt( bfs.getWeight ( q ) * vol );
        qc::parametricTransformAndClipCoord<ConfiguratorType, ParametricDeformationType> ( DF.getConfigurator(), *it, bfs.getRefCoord( q ), parDef, transformedEl, transformedLocalCoord );

        const RealType normedDeforedTEval = discrDefNormedTUndef.evaluate ( transformedEl, transformedLocalCoord );

        for ( int c = 0; c < ParametricDeformationType::NumberOfDeformParameters; ++c ) {
          MatDest.add ( i+q, c, - oneOverVarOfDeformedT * ( part1[c] + normedDeforedTEval * part2[c] ) * weight );
        }
      }
      i += numQuadPoints;
    }
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class SimpleParametricSSDEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const typename ConfiguratorType::InitType &_grid;
  const qc::SmoothedImage<ConfiguratorType> _r;
  const qc::SmoothedImage<ConfiguratorType> _t;
public:
  SimpleParametricSSDEnergy ( const typename ConfiguratorType::InitType &Grid,
                              const aol::Vector<RealType> &ImR,
                              const aol::Vector<RealType> &ImT )
    : _grid ( Grid ),
      _r ( Grid, ImR, false ),
      _t ( Grid, ImT, true ) { }

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    const ParametricDeformationType parDef ( _grid, MArg );
    typename ConfiguratorType::ArrayType deformedT ( _grid );
    qc::FastParametricDeformImage<ConfiguratorType, ParametricDeformationType, ConfiguratorType::Dim>
      ( _t.getSmoothImRef(), _grid, deformedT, parDef );

    deformedT -= _r.getSmoothImRef();
    Dest[0] += 0.5*deformedT.normSqr();
  }

  void applyDerivative ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    const ParametricDeformationType parDef ( _grid, MArg );

    aol::Mat< ParametricDeformationType::NumberOfDeformParameters, ConfiguratorType::Dim, RealType > jacobian;

    const aol::Vec3<int> gridSize = _grid.getSize();
    const RealType h = static_cast<RealType> ( _grid.H() );

    aol::Vector<RealType> dest ( ParametricDeformationType::NumberOfDeformParameters );

    for ( int j = 0; j < gridSize[1]; ++j ) {
      for ( int i = 0; i < gridSize[0]; ++i ) {
        const int index = qc::ILexCombine2( i, j, _grid.getNumX() );
        ParametricDeformedPosition<ConfiguratorType, ParametricDeformationType> defPos ( i, j, parDef, gridSize, h );

        if ( defPos.isDefPosInDomain() == false )
          continue;

        parDef.evaluateDerivativeDeformationOn01 ( defPos.getUndefPos01(), jacobian );
        const typename ConfiguratorType::VecType gradT ( _t.interpolateDX ( defPos.getDefPos() ), _t.interpolateDY ( defPos.getDefPos() ) );
        aol::Vec<ParametricDeformationType::NumberOfDeformParameters, RealType> tmp;

        jacobian.mult ( gradT, tmp );
        tmp *= ( _t.interpolate ( defPos.getDefPos() ) - _r.getSmoothImRef()[index] );
        for ( int k = 0; k < ParametricDeformationType::NumberOfDeformParameters; ++k )
          dest[k] += tmp[k];
      }
    }
    MDest.copySplitFrom ( dest );
  }
};

/**
 * \author Berkels
 */
template <typename RealType>
class SquaredL2DistanceEnergy : public aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > {
private:
  const aol::MultiVector<RealType> &_origin;
public:
  SquaredL2DistanceEnergy ( const aol::MultiVector<RealType> &Origin )
    : _origin ( Origin ) {}

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    Dest[0] *= 2;
    for ( int i = 0; i < _origin.numComponents(); ++i )
      for ( int j = 0; j < _origin[i].size(); ++j )
        Dest[0] += aol::Sqr ( MArg[i][j] - _origin[i][j] );
    Dest[0] /= 2;
  }

  void applyDerivative ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    MDest.setSum ( MArg, _origin, -1 );
  }
};
  

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
typename ConfiguratorType::RealType updateDeformParameters ( const typename ConfiguratorType::InitType &Grid,
                              const aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > &E,
                              const aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::MultiVector<typename ConfiguratorType::RealType> > &DE,
                              aol::MultiVector<typename ConfiguratorType::RealType> &DeformParameters,
                              const int MaxGradientDescentSteps = 1000,
                              const bool UseComponentWiseTimestep = true ) {
  typedef typename ConfiguratorType::RealType RealType;

  /*
  aol::FirstDerivativeValidator<aol::MultiVector<RealType> > tester ( E, DE, Grid.H(), aol::FirstDerivativeValidator<aol::MultiVector<RealType> >::LINEAR, 0.01 );
  tester.testDirection ( DeformParameters, "test/test" );
  tester.testAllDirections( DeformParameters, "test/test" );
  */

  aol::MultiVector<RealType> temp ( DeformParameters );
  aol::DeleteFlagPointer< aol::GradientDescentBase<RealType, aol::MultiVector<RealType> > > pGradientDescentSolver;
  typedef aol::GridlessGradientDescent<RealType, aol::MultiVector<RealType> > GradientDescentType;
  if ( UseComponentWiseTimestep ) {
    // Since GridlessGradientDescent doesn't need a grid, the Grid argument here is completely ignored
    // It's only necessary for the syntax.
    // Furthermore, we use a negative stopping epsilon here since GradientDescentComponentWiseTimestepControlled
    // only calculates tau for the different components approximately, possibly leading to a premature stopping.
    pGradientDescentSolver.reset ( new aol::GradientDescentComponentWiseTimestepControlled<ConfiguratorType, GradientDescentType> ( Grid, E, DE, MaxGradientDescentSteps, 1, -1000 ), true );
  }
  else {
    pGradientDescentSolver.reset ( new GradientDescentType ( Grid, E, DE, MaxGradientDescentSteps ), true );
  }
  pGradientDescentSolver->apply ( DeformParameters, temp );
  DeformParameters = temp;
  return pGradientDescentSolver->getEnergyAtLastPosition();
}

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
void ParametricDeformImage ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                             const typename ConfiguratorType::InitType &Grid,
                             aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                             const aol::MultiVector<typename ConfiguratorType::RealType> &DeformParameters,
                             const bool ExtendWithConstant = true,
                             const typename ConfiguratorType::RealType ExtensionConstant = 0,
                             const bool NearestNeighborInterpolation = false ) {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::ArrayType imageArray ( Image, Grid, aol::FLAT_COPY );

  qc::DataGenerator<ConfiguratorType> generator ( Grid );
  qc::MultiArray<RealType, ConfiguratorType::Dim> deformation ( Grid );
  ParametricDeformationType parDef ( Grid, DeformParameters );
  generator.template generateDeformationFromParametricDeformation<ParametricDeformationType, false> ( parDef, deformation );
  qc::DeformImage<ConfiguratorType> ( Image, Grid, DeformedImage, deformation, ExtendWithConstant, ExtensionConstant, NearestNeighborInterpolation );
}
  
/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricEnergyType>
typename ConfiguratorType::RealType updateDeformParameters ( const typename ConfiguratorType::InitType &Grid,
                              const aol::Vector<typename ConfiguratorType::RealType> &ImR,
                              const aol::Vector<typename ConfiguratorType::RealType> &ImT,
                              aol::MultiVector<typename ConfiguratorType::RealType> &DeformParameters,
                              const int MaxGradientDescentSteps = 1000,
                              const bool UseComponentWiseTimestep = true,
                              const typename ConfiguratorType::RealType ParameterPenaltyWeight = 0 ) {
  typedef typename ConfiguratorType::RealType RealType;
  ParametricEnergyType dataE ( Grid, ImR, ImT );
  aol::DerivativeWrapper<RealType, ParametricEnergyType, aol::MultiVector<RealType> > dataDE( dataE );
  if ( ParameterPenaltyWeight <= 0 )
    return updateDeformParameters<ConfiguratorType> ( Grid, dataE, dataDE, DeformParameters, MaxGradientDescentSteps, UseComponentWiseTimestep );
  else {
    SquaredL2DistanceEnergy<RealType> penaltyE ( DeformParameters );
    aol::DerivativeWrapper<RealType, SquaredL2DistanceEnergy<RealType>, aol::MultiVector<RealType> > penaltyDE( penaltyE );
    aol::LinCombOp<aol::MultiVector<RealType>, aol::Scalar<RealType> > E;
    E.appendReference ( dataE );
    E.appendReference ( penaltyE, ParameterPenaltyWeight );
    aol::LinCombOp<aol::MultiVector<RealType> > DE;
    DE.appendReference ( dataDE );
    DE.appendReference ( penaltyDE, ParameterPenaltyWeight );
    return updateDeformParameters<ConfiguratorType> ( Grid, E, DE, DeformParameters, MaxGradientDescentSteps, UseComponentWiseTimestep );
  }
}

/**
 * \author Berkels
 */
template <typename _ConfiguratorType, typename ParametricDeformationType>
class ParametricRegistrationMultilevelDescentBase : public qc::RegistrationMultilevelDescentInterfaceBase<_ConfiguratorType> {
public:
  typedef _ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef aol::MultiVector<RealType> TransformationDOFType;
  static const bool IsParametric = true;

protected:
  aol::MultiVector<RealType> _deformParameters;
  int _checkboxWidth, _maxGradientDescentSteps;
  bool _saveSteps;
  bool _useComponentWiseTimestep;
  RealType _parameterPenaltyWeight;
  RealType _energyOfLastSolution;

public:
  ParametricRegistrationMultilevelDescentBase ( const aol::ParameterParser &Parser, const bool TryToInitializeTranslation = false )
    : qc::RegistrationMultilevelDescentInterfaceBase<ConfiguratorType> ( Parser ),
      _deformParameters ( ParametricDeformationType::getDeformParametersSize() ),
      _checkboxWidth( Parser.getInt("checkboxWidth") ),
      _maxGradientDescentSteps( Parser.getInt("MaxGradientDescentSteps") ),
      _saveSteps( true ),
      _useComponentWiseTimestep ( Parser.getInt ( "UseComponentWiseTimestep" ) != 0 ),
      _parameterPenaltyWeight ( Parser.getRealOrDefault<RealType> ( "parameterPenaltyWeight", 0 ) ),
      _energyOfLastSolution ( aol::NumberTrait<RealType>::NaN ) {
    ParametricDeformationType::setIdentityDeformationParameters ( _deformParameters );

    if ( TryToInitializeTranslation ) {
      const aol::Vec<ConfiguratorType::Dim, RealType> centerRef = qc::getCenterOfMassOfArray<RealType, ConfiguratorType::Dim> ( this->getRefImageReference() ) / this->_grid.getNumX();
      const aol::Vec<ConfiguratorType::Dim, RealType> centerTem = qc::getCenterOfMassOfArray<RealType, ConfiguratorType::Dim> ( this->getTemplImageReference() ) / this->_grid.getNumX();
      const aol::Vec<ConfiguratorType::Dim, RealType> offset = ( centerTem - centerRef );

      for ( int i = 0; i < ConfiguratorType::Dim; ++i )
        _deformParameters[0][i] = offset[i]; 
    }
  }

  ParametricRegistrationMultilevelDescentBase ( const int MaxDepth, bool saveSteps = true )
    : qc::RegistrationMultilevelDescentInterfaceBase<ConfiguratorType> ( MaxDepth ),
      _deformParameters ( ParametricDeformationType::getDeformParametersSize() ),
      _checkboxWidth( 32 ),
      _maxGradientDescentSteps( 1000 ),
      _saveSteps( saveSteps ),
      _useComponentWiseTimestep ( true ),
      _parameterPenaltyWeight ( 0 ) {
    ParametricDeformationType::setIdentityDeformationParameters ( _deformParameters );
  }

  virtual ~ParametricRegistrationMultilevelDescentBase( ) {}

  void prolongate( ) {
    if ( this->_curLevel < this->getMaxGridDepth() )
      this->setLevel ( this->_curLevel + 1 );
  }

  aol::Vec<ParametricDeformationType::NumOfDeformParametersComponents, int> getTransformationDOFInitializer ( ) const {
    return ParametricDeformationType::getDeformParametersSize();
  }

  void setTransformation ( const aol::MultiVector<RealType> &DeformParameters ) {
    _deformParameters = DeformParameters;
  }

  void getTransformation ( aol::MultiVector<RealType> &DeformParameters ) const {
    DeformParameters = _deformParameters;
  }

  RealType getTransformationNorm ( ) const {
    return _deformParameters.norm();
  }

  void addTransformationTo ( aol::MultiVector<RealType> &DeformParameters ) const {
    ParametricDeformationType::concatenateDeformationParameters ( _deformParameters, DeformParameters );
  }

  void setTransformationToZero ( ) {
    ParametricDeformationType::setIdentityDeformationParameters ( _deformParameters );
  }

  void setTransformationToTranslation ( const aol::Vec<ConfiguratorType::Dim, RealType> &/*Translation*/ ) {
    throw aol::UnimplementedCodeException ( "not implemented", __FILE__, __LINE__ );
  }

  void composeDeformations ( const aol::MultiVector<RealType> &DeformParameters1, const aol::MultiVector<RealType> &DeformParameters2, aol::MultiVector<RealType> &CompositionDeformParameters ) {
    CompositionDeformParameters = DeformParameters2;
    ParametricDeformationType::concatenateDeformationParameters ( DeformParameters1, CompositionDeformParameters );
  }

  void setTransformationToComposition ( const aol::MultiVector<RealType> &/*DeformParameters1*/, const aol::MultiVector<RealType> &/*DeformParameters2*/ ) {
    throw aol::UnimplementedCodeException ( "not implemented", __FILE__, __LINE__ );
  }

  void applyTransformation ( const aol::MultiVector<RealType> &DeformParameters, const aol::Vector<RealType> &InputImage, aol::Vector<RealType> &DeformedImage, const RealType ExtensionConstant = aol::NumberTrait<RealType>::Inf, const bool NearestNeighborInterpolation = false ) const {
    qc::ParametricDeformImage<ConfiguratorType, ParametricDeformationType> ( InputImage, this->_grid, DeformedImage, DeformParameters, true, ExtensionConstant, NearestNeighborInterpolation );
  }

  void applySavedTransformation ( const char *DefBaseName, const aol::Vector<RealType> &InputImage, aol::Vector<RealType> &DeformedImage, const RealType ExtensionConstant = aol::NumberTrait<RealType>::Inf, const bool NearestNeighborInterpolation = false, const RealType xDerivNormThreshold = 0 ) const {
    aol::MultiVector<RealType> deformParameters;
    loadTransformationTo ( DefBaseName, deformParameters );
    applyTransformation ( deformParameters, InputImage, DeformedImage, ExtensionConstant, NearestNeighborInterpolation );
    if ( xDerivNormThreshold > 0 )
      throw aol::UnimplementedCodeException ( "not implemented", __FILE__, __LINE__ );
  }

  void loadTransformationTo ( const char *DefBaseName, aol::MultiVector<RealType> &DeformParameters ) const {
    DeformParameters.load ( aol::strprintf ( "%s%s", DefBaseName, getDeformationFileNameSuffix().c_str() ).c_str() );
  }

  void saveTransformationTo ( const aol::MultiVector<RealType> &DeformParameters, const char *DefBaseName ) const {
    DeformParameters.save ( aol::strprintf ( "%s%s", DefBaseName, getDeformationFileNameSuffix().c_str() ).c_str() );
  }

  int getMaxGradientSteps ( ) const {
    return _maxGradientDescentSteps;
  }

  void setMaxGradientSteps ( int maxGradientDescentSteps ) {
    _maxGradientDescentSteps = maxGradientDescentSteps;
  }

  virtual RealType updateDeformationOnOnCurrentGrid ( ) = 0;

  void descentOnCurrentGrid( ) {
    cerr << "Registration on level " << this->_curLevel << " started\n";

    if( this->_saveSteps )
      this->saveCurrentDeformation ( "_before", false );

    this->_energyOfLastSolution = updateDeformationOnOnCurrentGrid();

    cerr << "Detected parameters are \n" << this->_deformParameters << endl;

    if( this->_saveSteps )
      this->saveCurrentDeformation ( );
  }

  void setCheckboxWidth ( int checkboxWidth ) {
    _checkboxWidth = checkboxWidth;
  }

  string getDeformationFileNameSuffix ( ) const {
    return ".dat";
  }

  void saveCurrentDeformation ( const char *FileNameAddition = NULL, const bool SaveDisplacement = true ) const {
    qc::DataGenerator<ConfiguratorType> generator ( this->_grid );
    qc::MultiArray<RealType, ConfiguratorType::Dim> deformation ( this->_grid );
    ParametricDeformationType parDef ( this->_grid, _deformParameters );
    generator.template generateDeformationFromParametricDeformation<ParametricDeformationType, false> ( parDef, deformation );

    qc::RegistrationStepSaver<ConfiguratorType, typename ConfiguratorType::ArrayType, false, false>
      stepSaver( this->_grid, this->getRefImageReference(), this->getTemplImageReference(), _checkboxWidth  );
    stepSaver.setSaveName ( aol::strprintf ( "%s_%02d", FileNameAddition ? FileNameAddition : "", this->_curLevel ).c_str() );
    stepSaver.setSaveDirectory ( this->getSaveDirectory() );
    if ( this->getParserReference().checkAndGetBool ( "onlySaveDisplacement" ) == false )
      stepSaver.saveStep ( deformation, -1 );
    if ( SaveDisplacement )
      _deformParameters.save ( stepSaver.createSaveName ( "deformation", aol::strprintf ( "%s.dat", FileNameAddition ? FileNameAddition : "" ).c_str(), -1 ).c_str() );
  }

  RealType getEnergyOfLastSolution ( ) const {
    return _energyOfLastSolution;
  }

  void applyCurrentTransformation ( const aol::Vector<RealType> &InputImage, aol::Vector<RealType> &DeformedImage ) const {
    applyTransformation ( _deformParameters, InputImage, DeformedImage );
  }

  void saveTransformation ( const char *FileName ) const {
    _deformParameters.save ( FileName );
  }

  void loadTransformation ( const char *FileName ) {
    _deformParameters.load ( FileName );
  }

};

/**
 * \author Berkels
 */
template <typename _ConfiguratorType, typename ParametricDeformationType, typename ParametricEnergyType>
class ParametricRegistrationMultilevelDescent : public ParametricRegistrationMultilevelDescentBase<_ConfiguratorType, ParametricDeformationType> {
public:
  typedef _ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef aol::MultiVector<RealType> TransformationDOFType;
  static const bool IsParametric = true;

public:
  ParametricRegistrationMultilevelDescent ( const aol::ParameterParser &Parser, const bool TryToInitializeTranslation = false )
    : ParametricRegistrationMultilevelDescentBase<_ConfiguratorType, ParametricDeformationType> ( Parser, TryToInitializeTranslation ) { }

  ParametricRegistrationMultilevelDescent ( const int MaxDepth, bool saveSteps = true )
    : ParametricRegistrationMultilevelDescentBase<_ConfiguratorType, ParametricDeformationType> ( MaxDepth, saveSteps ) {}

  virtual ~ParametricRegistrationMultilevelDescent( ) {}

  RealType updateDeformationOnOnCurrentGrid ( ) {
    return updateDeformParameters<ConfiguratorType, ParametricEnergyType> ( this->getCurrentGrid(), this->_org_reference[ this->_curLevel ], this->_org_template[ this->_curLevel ], this->_deformParameters, this->_maxGradientDescentSteps, this->_useComponentWiseTimestep, this->_parameterPenaltyWeight );
  }
};

/**
 * \author Berkels
 */
template <typename _ConfiguratorType, typename ParametricDeformationType, typename ParametricObjectiveType>
class ParametricRegistrationMultilevelDescentGN : public ParametricRegistrationMultilevelDescentBase<_ConfiguratorType, ParametricDeformationType> {
public:
  typedef _ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ParametricObjectiveType::DerivativeType DerivativeType;
  typedef aol::MultiVector<RealType> TransformationDOFType;
  static const bool IsParametric = true;

public:
  ParametricRegistrationMultilevelDescentGN ( const aol::ParameterParser &Parser, const bool TryToInitializeTranslation = false )
    : ParametricRegistrationMultilevelDescentBase<_ConfiguratorType, ParametricDeformationType> ( Parser, TryToInitializeTranslation ) { }

  ParametricRegistrationMultilevelDescentGN ( const int MaxDepth, bool saveSteps = true )
    : ParametricRegistrationMultilevelDescentBase<_ConfiguratorType, ParametricDeformationType> ( MaxDepth, saveSteps ) {}

  virtual ~ParametricRegistrationMultilevelDescentGN( ) {}

  RealType updateDeformationOnOnCurrentGrid ( ) {
    const int dimRangeF = this->getCurrentGrid().getNumberOfElements() * ConfiguratorType::QuadType::numQuadPoints;

    ParametricObjectiveType F ( this->getCurrentGrid(), this->_org_reference[ this->_curLevel ], this->_org_template[ this->_curLevel ] );
    aol::DerivativeWrapper<RealType, ParametricObjectiveType, aol::Vector<RealType>, DerivativeType> DF ( F );

    aol::Vector<RealType> deformParametersVec ( this->_deformParameters.getTotalSize() );
    deformParametersVec.copyUnblockedFrom ( this->_deformParameters );
    aol::GaussNewtonAlgorithm<aol::Vector<RealType>, DerivativeType> gnSolver ( dimRangeF, deformParametersVec.size(), F, DF, this->_maxGradientDescentSteps, this->getParserReference().getDoubleOrDefault ( "stopEpsilon", 0 ) );
    gnSolver.applySingle ( deformParametersVec );
    this->_deformParameters.copySplitFrom ( deformParametersVec );
    return gnSolver.getFNormSqrAtLastPosition();
  }
};

/**
 * \brief Registers two images using phase correlation to find the optimal integer shift between the two images.
 *
 * \author Berkels
 */
template<typename ConfiguratorType>
class PhaseCorrelationRegistration {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  static const qc::Dimension Dim = ConfiguratorType::Dim;
public:
  //! \brief Retuns the translation vector \f$b\f$ (in pixels) such that \f$R(\cdot+b)\approx T\f$.
  //! \note If the allowed shift is limited, the normal correlation is used instead of the phase correlation.
  static qc::CoordType registerImages ( const ArrayType &Reference, const ArrayType &Template, const int MaxInfAbsShift = 0 ) {
    const typename ConfiguratorType::InitType grid ( Reference.getSize() );
    const aol::Vec3<int> imageSize = Reference.getSize();
    ArrayType phaseCorrelation ( grid );

    const bool limitShift = ( MaxInfAbsShift > 0 );

    qc::Convolution<Dim, RealType> conv ( qc::GridSize<Dim> ( grid ).getSizeAsVecDim() );
    conv.phaseCorrelation  ( Reference, Template, phaseCorrelation, !limitShift );

    qc::CoordType shift;
    if ( !limitShift ) {
      qc::OTFILexMapper<Dim> mapper ( grid );
      shift = mapper.splitGlobalIndex ( phaseCorrelation.getMaxIndexAndValue().first );
    }
    else {
      shift.setZero();
      RealType maxCorr = phaseCorrelation.get ( shift );
      typedef typename aol::VecDimTrait<int, ConfiguratorType::Dim>::VecType VecType;
      for ( qc::RectangularIterator<ConfiguratorType::Dim> it ( VecType ( -MaxInfAbsShift ), VecType ( MaxInfAbsShift + 1 )  ); it.notAtEnd(); ++it ) {

        VecType pos;
        for ( int i = 0; i < Dim; ++i )
          pos[i] = ( (*it)[i] + imageSize[i] ) % imageSize[i];

        if ( phaseCorrelation.get ( pos ) > maxCorr ) {
          maxCorr = phaseCorrelation.get ( pos );
          for ( int i = 0; i < Dim; ++i )
            shift[i] = pos[i];
        }
      }
    }

    qc::CoordType signedShift;
    for ( int i = 0; i < Dim; ++i )
      signedShift[i] = ( ( shift[i] <= imageSize[i] / 2 ) ? shift[i] : shift[i] - imageSize[i] );
    return signedShift;
  }
};
  
/**
 *\brief Class for affine transformation. It allows translation, rotation, scaling and shearing
 *\author Tatano
 *
 */
template <typename ConfiguratorType>
class ParametricAffineTransformation2D {
public:
  static const int NumberOfDeformParameters = 2 + 2 * ConfiguratorType::Dim ;
  static const int NumOfDeformParametersComponents = 4;
  typedef typename ConfiguratorType::RealType RealType;
private:
  const typename ConfiguratorType::InitType &_grid;
  RealType _rotationAngle;
  aol::Matrix22<RealType> _rotationMatrix;
  aol::Matrix22<RealType> _shearingMatrix;
  typename ConfiguratorType::VecType _translation;
  typename ConfiguratorType::VecType _shearing;
  RealType _scaling;
public:
  ParametricAffineTransformation2D ( const typename ConfiguratorType::InitType &Initializer, const aol::MultiVector<RealType> &DeformParameters )
  : _grid ( Initializer ),
  _rotationAngle ( DeformParameters[1][0] ),
  _scaling ( DeformParameters[2][0] ){
    _rotationMatrix.makeRotation ( _rotationAngle );
    _shearingMatrix.setIdentity();
    for ( int i = 0; i < ConfiguratorType::Dim; ++i ){
      _translation[i] = DeformParameters[0][i] + 0.5;
      _shearing[i] = DeformParameters[3][i];
    }
    _shearingMatrix[1][0] = _shearing[1];
    _shearingMatrix[0][1] = _shearing[0];
  }
    
  ParametricAffineTransformation2D ( const typename ConfiguratorType::InitType &Initializer, const aol::Vector<RealType> &DeformParametersVec )
    : _grid ( Initializer ) {
    setDeformationParametersAsVec ( DeformParametersVec );
  }

  void setDeformationParametersAsVec ( const aol::Vector<RealType> &DeformParametersVec ) {
    _rotationAngle = DeformParametersVec[2];
    _scaling = DeformParametersVec[3];
    _rotationMatrix.makeRotation ( _rotationAngle );
    _shearingMatrix.setIdentity();
    for ( int i = 0; i < ConfiguratorType::Dim; ++i ){
      _translation[i] = DeformParametersVec[i] + 0.5;
      _shearing[i] = DeformParametersVec[4+i];
    }
    _shearingMatrix[1][0] = _shearing[1];
    _shearingMatrix[0][1] = _shearing[0];
  }

  static void setIdentityDeformationParameters ( aol::MultiVector<RealType> &DeformParameters ) {
    DeformParameters.setZero();
    DeformParameters[ConfiguratorType::Dim][0] = 1;
  }
    
  // DestDeform := ArgDeform \circ DestDeform
  static void concatenateDeformationParameters ( const aol::MultiVector<RealType> &ArgDeformPars, aol::MultiVector<RealType> &DestDeformPars ) {
    throw aol::Exception("Not implemented yet!", __FILE__, __LINE__);
  }
  
  static aol::Vec<NumOfDeformParametersComponents, int> getDeformParametersSize ( ) {
    aol::Vec<NumOfDeformParametersComponents, int> size;
    size[0] = 2;
    size[1] = 1;
    size[2] = 1;
    size[3] = 2;
    return size;
  }
    
  template <bool ClipCoord>
  bool evaluateDeformation ( const typename ConfiguratorType::ElementType &El,
                              const typename ConfiguratorType::VecType &RefCoord,
                              qc::Element &TransformedEl, typename ConfiguratorType::VecType &TransformedLocalCoord ) const {
      
    typename ConfiguratorType::VecType coord, transformedCoord, tmp;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      coord[i] = El[i] + RefCoord[i] - 0.5 / _grid.H();
    }
    _shearingMatrix.mult ( coord, tmp );
    _rotationMatrix.mult ( tmp, transformedCoord );
    transformedCoord *= _scaling;
      
    return qc::transformCoord<ConfiguratorType, ClipCoord> ( _grid, transformedCoord, _translation, TransformedEl, TransformedLocalCoord );
  }
    
  void evaluateDeformationOn01 ( const ParametricAffineTransformation2D<ConfiguratorType> &/*ParDef*/, const typename ConfiguratorType::VecType &Position, typename ConfiguratorType::VecType &DeformedPosition ) const {
    typename ConfiguratorType::VecType tmp ( Position ), tmp2(tmp, aol::STRUCT_COPY);
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      tmp[i] -= 0.5;
    _shearingMatrix.apply ( tmp, tmp2 );
    _rotationMatrix.apply ( tmp2, DeformedPosition );
    DeformedPosition *= _scaling;
    DeformedPosition += _translation;
  }
    
  void evaluateDerivativeDeformationOn01 ( const typename ConfiguratorType::VecType &Position, aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    const RealType cosAlpha = _rotationMatrix[0][0];
    const RealType sinAlpha = _rotationMatrix[1][0];
      
    const RealType tempX = (Position[0] - 0.5) + _shearing[0] * (Position[1] - 0.5);
    const RealType tempY = (Position[1] - 0.5) + _shearing[1] * (Position[0] - 0.5);
      
    Jacobian[ConfiguratorType::Dim][0] = _scaling * ( -sinAlpha * tempX - cosAlpha * tempY );
    Jacobian[ConfiguratorType::Dim][1] = _scaling * (  cosAlpha * tempX - sinAlpha * tempY );
      
    Jacobian[0][0] = Jacobian[1][1] = 1;
    Jacobian[1][0] = Jacobian[0][1] = 0;
      
    Jacobian[ConfiguratorType::Dim+1][0] = ( cosAlpha * tempX - sinAlpha * tempY );
    Jacobian[ConfiguratorType::Dim+1][1] = ( sinAlpha * tempX + cosAlpha * tempY );
      
    Jacobian[ConfiguratorType::Dim+2][0] =   _scaling * cosAlpha * (Position[1] - 0.5);
    Jacobian[ConfiguratorType::Dim+2][1] = - _scaling * sinAlpha * (Position[0] - 0.5);
    Jacobian[ConfiguratorType::Dim+3][1] =   _scaling * sinAlpha * (Position[1] - 0.5);
    Jacobian[ConfiguratorType::Dim+3][0] =   _scaling * cosAlpha * (Position[0] - 0.5);
  }
    
  void evaluateDerivativeDeformation ( const typename ConfiguratorType::ElementType &El,
                                        const typename ConfiguratorType::VecType &RefCoord,
                                        aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    typename ConfiguratorType::VecType x;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      x[i] = ( El[i] + RefCoord[i] ) * _grid.H();
    }
    evaluateDerivativeDeformationOn01 ( x, Jacobian );
  }
};
  
/**
 *\brief Class for fully affine transformation
 *\author Tatano
 *
 */
template <typename ConfiguratorType>
class FullyAffineTransformation2D {
public:
  static const int NumberOfDeformParameters = 6 ;
  static const int NumOfDeformParametersComponents = 2;
  typedef typename ConfiguratorType::RealType RealType;
private:
  const typename ConfiguratorType::InitType &_grid;
  aol::Matrix22<RealType> _affineMatrix;
  typename ConfiguratorType::VecType _translation;
public:
  FullyAffineTransformation2D ( const typename ConfiguratorType::InitType &Initializer, const aol::MultiVector<RealType> &DeformParameters )
  : _grid ( Initializer ){
    _affineMatrix[0][0] = DeformParameters[1][0];
    _affineMatrix[0][1] = DeformParameters[1][1];
    _affineMatrix[1][0] = DeformParameters[1][2];
    _affineMatrix[1][1] = DeformParameters[1][3];
    _translation[0] = DeformParameters[0][0] + 0.5;
    _translation[1] = DeformParameters[0][1] + 0.5;
  }
    
  static void setIdentityDeformationParameters ( aol::MultiVector<RealType> &DeformParameters ) {
    DeformParameters.setZero();
    DeformParameters[1][0] = 1;
    DeformParameters[1][3] = 1;
  }
  
  static void setDeformationParametersToShrinked ( aol::MultiVector<RealType> &DeformParameters ) {
    DeformParameters.setZero();
    DeformParameters[1][0] = 1./2;
    DeformParameters[1][3] = 1./2;
  }
    
  // DestDeform := ArgDeform \circ DestDeform
  static void concatenateDeformationParameters ( const aol::MultiVector<RealType> &ArgDeformPars, aol::MultiVector<RealType> &DestDeformPars ) {

    aol::Matrix22<RealType> argMat, destMat;
    argMat.fill ( ArgDeformPars[1] );
    destMat.fill ( DestDeformPars[1] );

    aol::Vec2<RealType> b, destB;
    for ( int i = 0; i < ConfiguratorType::Dim; ++i ) {
      b[i] = ArgDeformPars[0][i];
      destB[i] = DestDeformPars[0][i];
    }

    argMat.multAdd ( destB, b );
    argMat *= destMat;

    DestDeformPars[1][0] = argMat[0][0];
    DestDeformPars[1][1] = argMat[0][1];
    DestDeformPars[1][2] = argMat[1][0];
    DestDeformPars[1][3] = argMat[1][1];

    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      DestDeformPars[0][i] = b[i];
  }
  
  static aol::Vec<NumOfDeformParametersComponents, int> getDeformParametersSize ( ) {
    return aol::Vec2<int> (2, 4);
  }
    
  template <bool ClipCoord>
  bool evaluateDeformation ( const typename ConfiguratorType::ElementType &El,
                              const typename ConfiguratorType::VecType &RefCoord,
                              qc::Element &TransformedEl, typename ConfiguratorType::VecType &TransformedLocalCoord ) const {
      
    typename ConfiguratorType::VecType coord, transformedCoord;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      coord[i] = El[i] + RefCoord[i] - 0.5 / _grid.H();
    }
    _affineMatrix.apply(coord, transformedCoord);
    return qc::transformCoord<ConfiguratorType, ClipCoord> ( _grid, transformedCoord, _translation, TransformedEl, TransformedLocalCoord );
  }
    
  void evaluateDeformationOn01 ( const FullyAffineTransformation2D<ConfiguratorType> &/*ParDef*/, const typename ConfiguratorType::VecType &Position, typename ConfiguratorType::VecType &DeformedPosition ) const {
    typename ConfiguratorType::VecType tmp ( Position );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      tmp[i] -= 0.5;
    _affineMatrix.apply ( tmp, DeformedPosition );
    DeformedPosition += _translation;
  }
    
  void evaluateDerivativeDeformationOn01 ( const typename ConfiguratorType::VecType &Position, aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    Jacobian[0][0] = 1;                       Jacobian[0][1] = 0;
    Jacobian[1][0] = 0;                       Jacobian[1][1] = 1;
    Jacobian[2][0] = Position[0] - 0.5;       Jacobian[2][1] = 0;
    Jacobian[3][0] = Position[1] - 0.5;       Jacobian[3][1] = 0;
    Jacobian[4][0] = 0;                       Jacobian[4][1] = Position[0] - 0.5;
    Jacobian[5][0] = 0;                       Jacobian[5][1] = Position[1] - 0.5;
  }
    
  void evaluateDerivativeDeformation ( const typename ConfiguratorType::ElementType &El,
                                        const typename ConfiguratorType::VecType &RefCoord,
                                        aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    typename ConfiguratorType::VecType x;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      x[i] = ( El[i] + RefCoord[i] ) * _grid.H();
    }
    evaluateDerivativeDeformationOn01 ( x, Jacobian );
  }
};

  /**
   *\brief Class for fully affine transformation
   *\author Tatano
   *
   */
  template <typename ConfiguratorType>
  class FullyAffineTransformation3D {
  public:
    static const int NumberOfDeformParameters = 12 ;
    static const int NumOfDeformParametersComponents = 2;
    typedef typename ConfiguratorType::RealType RealType;
  private:
    const typename ConfiguratorType::InitType &_grid;
    aol::Matrix33<RealType> _affineMatrix;
    typename ConfiguratorType::VecType _translation;
  public:
    FullyAffineTransformation3D ( const typename ConfiguratorType::InitType &Initializer, const aol::MultiVector<RealType> &DeformParameters )
    : _grid ( Initializer ){
      _affineMatrix[0][0] = DeformParameters[1][0];
      _affineMatrix[0][1] = DeformParameters[1][1];
      _affineMatrix[0][2] = DeformParameters[1][2];
      _affineMatrix[1][0] = DeformParameters[1][3];
      _affineMatrix[1][1] = DeformParameters[1][4];
      _affineMatrix[1][2] = DeformParameters[1][5];
      _affineMatrix[2][0] = DeformParameters[1][6];
      _affineMatrix[2][1] = DeformParameters[1][7];
      _affineMatrix[2][2] = DeformParameters[1][8];
      _translation[0] = DeformParameters[0][0] + 0.5;
      _translation[1] = DeformParameters[0][1] + 0.5;
      _translation[2] = DeformParameters[0][2] + 0.5;
    }
    
    static void setIdentityDeformationParameters ( aol::MultiVector<RealType> &DeformParameters ) {
      DeformParameters.setZero();
      DeformParameters[1][0] = 1;
      DeformParameters[1][4] = 1;
      DeformParameters[1][8] = 1;
    }
    
    static void setDeformationParametersToShrinked ( aol::MultiVector<RealType> &DeformParameters ) {
      DeformParameters.setZero();
      DeformParameters[1][0] = 1./2;
      DeformParameters[1][4] = 1./2;
      DeformParameters[1][8] = 1./2;
    }
    
    // DestDeform := ArgDeform \circ DestDeform
    static void concatenateDeformationParameters ( const aol::MultiVector<RealType> &ArgDeformPars, aol::MultiVector<RealType> &DestDeformPars ) {
      //TODO
      throw aol::Exception("Not implemented yet!", __FILE__, __LINE__);
    }
    
    static aol::Vec<NumOfDeformParametersComponents, int> getDeformParametersSize ( ) {
      return aol::Vec2<int> (3, 9);
    }
    
    template <bool ClipCoord>
    bool evaluateDeformation ( const typename ConfiguratorType::ElementType &El,
                              const typename ConfiguratorType::VecType &RefCoord,
                              qc::Element &TransformedEl, typename ConfiguratorType::VecType &TransformedLocalCoord ) const {
      
      typename ConfiguratorType::VecType coord, transformedCoord;
      for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
        coord[i] = El[i] + RefCoord[i] - 0.5 / _grid.H();
      }
      _affineMatrix.apply(coord, transformedCoord);
      return qc::transformCoord<ConfiguratorType, ClipCoord> ( _grid, transformedCoord, _translation, TransformedEl, TransformedLocalCoord );
    }
    
    void evaluateDeformationOn01 ( const FullyAffineTransformation3D<ConfiguratorType> &/*ParDef*/, const typename ConfiguratorType::VecType &Position, typename ConfiguratorType::VecType &DeformedPosition ) const {
      typename ConfiguratorType::VecType tmp ( Position );
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        tmp[i] -= 0.5;
      _affineMatrix.apply ( tmp, DeformedPosition );
      DeformedPosition += _translation;
    }
    
    void evaluateDerivativeDeformationOn01 ( const typename ConfiguratorType::VecType &Position, aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
      Jacobian[0][0] = 1;                       Jacobian[0][1] = 0;                   Jacobian[0][2] = 0;
      Jacobian[1][0] = 0;                       Jacobian[1][1] = 1;                   Jacobian[1][2] = 0;
      Jacobian[2][0] = 0;                       Jacobian[2][1] = 0;                   Jacobian[2][2] = 1;
      Jacobian[3][0] = Position[0] - 0.5;       Jacobian[3][1] = 0;                   Jacobian[3][2] = 0;
      Jacobian[4][0] = Position[1] - 0.5;       Jacobian[4][1] = 0;                   Jacobian[4][2] = 0;
      Jacobian[5][0] = Position[2] - 0.5;       Jacobian[5][1] = 0;                   Jacobian[5][2] = 0;
      Jacobian[6][0] = 0;                       Jacobian[6][1] = Position[0] - 0.5;   Jacobian[6][2] = 0;
      Jacobian[7][0] = 0;                       Jacobian[7][1] = Position[1] - 0.5;   Jacobian[7][2] = 0;
      Jacobian[8][0] = 0;                       Jacobian[8][1] = Position[2] - 0.5;   Jacobian[8][2] = 0;
      Jacobian[9][0] = 0;                       Jacobian[9][1] = 0;                   Jacobian[9][2] = Position[0] - 0.5;
      Jacobian[10][0] = 0;                      Jacobian[10][1] = 0;                  Jacobian[10][2] = Position[1] - 0.5;
      Jacobian[11][0] = 0;                      Jacobian[11][1] = 0;                  Jacobian[11][2] = Position[2] - 0.5;
    }
    
    void evaluateDerivativeDeformation ( const typename ConfiguratorType::ElementType &El,
                                        const typename ConfiguratorType::VecType &RefCoord,
                                        aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
      typename ConfiguratorType::VecType x;
      for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
        x[i] = ( El[i] + RefCoord[i] ) * _grid.H();
      }
      evaluateDerivativeDeformationOn01 ( x, Jacobian );
    }
  };

/**
 * \author Berkels
 */
template<typename ConfiguratorType, typename ParametricDeformationType>
void parametricSeriesRegis ( const aol::RandomAccessContainer<typename ConfiguratorType::ArrayType> &Images,
                             const aol::Vector<typename ConfiguratorType::RealType> &PairParDefParams,
                             aol::MultiVector<typename ConfiguratorType::RealType> &ParDefParams,
                             const bool UseCorrelationToInitTranslation = false,
                             const int MaxCorrShift = 0 ) {
  const typename ConfiguratorType::InitType grid ( qc::GridSize<ConfiguratorType::Dim>::createFrom ( Images[0] ) );

  typedef typename ConfiguratorType::RealType RealType;
  aol::MultiVector<RealType> pairParDefParamMVec ( ParametricDeformationType::getDeformParametersSize() );
  pairParDefParamMVec.copySplitFrom ( PairParDefParams );

  if ( PairParDefParams.size() != ParametricDeformationType::NumberOfDeformParameters )
    throw aol::Exception ( "pairParDefParams size mismatch", __FILE__, __LINE__ );

  for ( int i = 0; i < Images.size() - 1; ++i ) {
    aol::MultiVector<RealType> parDefParamMVec ( pairParDefParamMVec );

    if ( i > 0 ) {
      aol::MultiVector<RealType> translationMVec ( ParametricDeformationType::getDeformParametersSize() );
      translationMVec.copySplitFrom ( ParDefParams[i-1] );
      ParametricDeformationType::concatenateDeformationParameters ( translationMVec, parDefParamMVec );
    }

    if ( UseCorrelationToInitTranslation ) {
      // Remove the old translation.
      parDefParamMVec[0].setZero();

      // Estimate translation using the optimal integer shift obtained with phase correlation.
      typename ConfiguratorType::ArrayType deformedImage ( grid );
      const ParametricDeformationType parDef ( grid, parDefParamMVec );
      qc::FastParametricInvDeformImage<ConfiguratorType, ParametricDeformationType, ConfiguratorType::Dim>
      ( Images[i+1], grid, deformedImage, parDefParamMVec );

      qc::CoordType signedShift = qc::PhaseCorrelationRegistration<ConfiguratorType>::registerImages ( Images[0], deformedImage, MaxCorrShift );

      for ( int k = 0; k < ConfiguratorType::Dim; ++k )
        parDefParamMVec[0][k] = signedShift[k]*grid.H();
    }

    qc::updateDeformParameters<ConfiguratorType,  qc::ParametricNCCEnergy<ConfiguratorType, ParametricDeformationType> > ( grid, Images[i+1], Images[0], parDefParamMVec, 1000, false );

    ParDefParams[i].copyUnblockedFrom ( parDefParamMVec );
  }
}

/**
 * \author Berkels
 */
template<typename ConfiguratorType, typename ParametricDeformationType>
void parametricSeriesRegis ( const aol::RandomAccessContainer<typename ConfiguratorType::ArrayType> &Images,
                             aol::MultiVector<typename ConfiguratorType::RealType> &ParDefParams,
                             const bool UseCorrelationToInitTranslation = false,
                             const int MaxCorrShift = 0 ) {
  typedef typename ConfiguratorType::RealType RealType;
  aol::MultiVector<RealType> indentityParDefParamMVec ( ParametricDeformationType::getDeformParametersSize() );
  ParametricDeformationType::setIdentityDeformationParameters ( indentityParDefParamMVec );
  aol::Vector<RealType> indentityParDefParams ( ParametricDeformationType::NumberOfDeformParameters );
  indentityParDefParams.copyUnblockedFrom( indentityParDefParamMVec );

  qc::parametricSeriesRegis<ConfiguratorType, ParametricDeformationType>
    ( Images, indentityParDefParams, ParDefParams, UseCorrelationToInitTranslation, MaxCorrShift );
}

/**
 * \brief Works with the output from qc::parametricSeriesRegis.
 *
 * \author Berkels
 */
template<typename ConfiguratorType, typename ParametricDeformationType>
void averageSeries ( const aol::RandomAccessContainer<typename ConfiguratorType::ArrayType> &Images,
                     aol::MultiVector<typename ConfiguratorType::RealType> &ParDefParams,
                     typename ConfiguratorType::ArrayType &Average ) {
  const typename ConfiguratorType::InitType grid ( qc::GridSize<ConfiguratorType::Dim>::createFrom ( Images[0] ) );
  typedef typename ConfiguratorType::RealType RealType;
  Average.reallocate ( grid );
  aol::RandomAccessContainer<typename ConfiguratorType::ArrayType> deformedImages ( Images.size() - 1, grid );

  for ( int i = 0; i < Images.size() - 1; ++i ) {
    aol::MultiVector<RealType> parDefParamMVec ( ParametricDeformationType::getDeformParametersSize() );
    parDefParamMVec.copySplitFrom ( ParDefParams[i] );
    aol::MultiVector<RealType> inverseParDefParam ( ParametricDeformationType::getDeformParametersSize() );
    ParametricDeformationType::getDeformationParametersOfInverse ( parDefParamMVec, inverseParDefParam );
    qc::ParametricDeformImage<ConfiguratorType, ParametricDeformationType> ( Images[i+1], grid, deformedImages[i], inverseParDefParam, true, aol::NumberTrait<RealType>::Inf );
  }

  for ( int i = 0; i < Average.size(); ++i ) {
    Average[i] = Images[0][i];
    int numSamplesInDomain = 1;
    for ( int j = 0; j < deformedImages.size(); ++j ) {
      if ( deformedImages[j][i] != aol::NumberTrait<RealType>::Inf ) {
        Average[i] += deformedImages[j][i];
        ++numSamplesInDomain;
      }
    }
    Average[i] /= numSamplesInDomain;
  }
}

/**
 * \author Berkels, Wirth
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
void getCentralPixelOfCommonSupportInAllImages ( const typename ConfiguratorType::InitType &Grid,
                                                 const aol::MultiVector<typename ConfiguratorType::RealType> &ParDefParams,
                                                 aol::RandomAccessContainer<aol::Vec2<typename ConfiguratorType::RealType> > &Centers ) {
  typedef typename ConfiguratorType::RealType RealType;
  const int numImages = ParDefParams.numComponents() + 1;

  // maximum and minimum common pixels of all images, indexed by position in 1st image
  RealType minX = 0, minY = 0, maxX = Grid.getNumX()-1, maxY = Grid.getNumY()-1;
  for ( int i = 0; i < numImages-1; ++i ) {
    // corners of (i+1)th image
    aol::Vec2<RealType> origin( 0, 0 ), topRight( Grid.getNumX()-1, Grid.getNumY()-1 );
    // deformation
    ParametricDeformationType parDef ( Grid, ParDefParams[i] );
    parDef.evaluateDeformationOn0N ( origin );
    parDef.evaluateDeformationOn0N ( topRight );
    // identify pixel position
    minX = aol::Max( minX, aol::Min( origin[0], topRight[0] ) );
    minY = aol::Max( minY, aol::Min( origin[1], topRight[1] ) );
    maxX = aol::Min( maxX, aol::Max( origin[0], topRight[0] ) );
    maxY = aol::Min( maxY, aol::Max( origin[1], topRight[1] ) );
  }
  // identify central pixel of the common pixels of all input images
  Centers.reallocate ( numImages, aol::Vec2<RealType> ( ( minX + maxX ) / 2, ( minY + maxY ) / 2 ) );
  for ( int i = 0; i < numImages-1; ++i ) {
    // compute central common pixel in (i+1)th image by applying the inverse deformation
    aol::MultiVector<RealType> inverseParDefParam ( ParametricDeformationType::getDeformParametersSize() );
    aol::MultiVector<RealType> parDefParam ( ParametricDeformationType::getDeformParametersSize() );
    parDefParam.copySplitFrom ( ParDefParams[i] );
    ParametricDeformationType::getDeformationParametersOfInverse ( parDefParam, inverseParDefParam );
    ParametricDeformationType inverseParDef ( Grid, inverseParDefParam );

    inverseParDef.evaluateDeformationOn0N ( Centers[i+1] );
  }
}
/**
 * \author Berkels, Wirth
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
void computeSeriesCropParameters ( const typename ConfiguratorType::InitType &Grid,
                                   const aol::MultiVector<typename ConfiguratorType::RealType> &ParDefParams,
                                   aol::MultiVector<int> &CropStart,
                                   int &CropSizeX ) {
  typedef typename ConfiguratorType::RealType RealType;
  const int numImages = ParDefParams.numComponents() + 1;

  aol::RandomAccessContainer<aol::Vec2<RealType> > centers;
  getCentralPixelOfCommonSupportInAllImages<ConfiguratorType, ParametricDeformationType> ( Grid, ParDefParams, centers );

  // find crop start values
  RealType cropWidth = 0;
  if ( CropSizeX > 0 )
    cropWidth = static_cast<RealType>( CropSizeX ) / 2.;
  else {
    aol::Vec2<int> minCenter ( Grid.getNumX(), Grid.getNumY() );
    aol::Vec2<int> maxCenter ( 0, 0 );
    for ( int i = 0; i < numImages; ++i ) {
      for ( int j = 0; j < 2; ++j ) {
        minCenter[j] = aol::Min ( minCenter[j], aol::Rint( centers[i][j] ) );
        maxCenter[j] = aol::Max ( maxCenter[j], aol::Rint( centers[i][j] ) );
      }
    }
    cropWidth = aol::Min ( minCenter.getMinValue(), Grid.getNumX() - 1 - maxCenter.getMaxValue() );
    CropSizeX = 2*cropWidth;
  }
  for ( int i = 0; i < numImages; ++i ) {
    CropStart[i][0] = aol::Rint( centers[i][0] - cropWidth );
    CropStart[i][1] = aol::Rint( centers[i][1] - cropWidth );
  }
  if ( ( CropStart.getMinValue() < 0 ) || ( CropStart.getMaxValue() >= Grid.getNumX() - 2 * cropWidth ) )
    throw aol::Exception ( "cropwidth too large", __FILE__, __LINE__ );
}

/**
 * \author Berkels, Wirth
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
void computeSeriesCropParameters ( const typename ConfiguratorType::InitType &Grid,
                                   const aol::MultiVector<typename ConfiguratorType::RealType> &ParDefParams,
                                   aol::MultiVector<int> &CropStart,
                                   aol::Vec2<int> &CropSize ) {
  typedef typename ConfiguratorType::RealType RealType;
  const int numImages = ParDefParams.numComponents() + 1;

  aol::RandomAccessContainer<aol::Vec2<RealType> > centers;
  qc::getCentralPixelOfCommonSupportInAllImages<ConfiguratorType, ParametricDeformationType> ( Grid, ParDefParams, centers );

  // find crop start values
  aol::Vec2<RealType> cropWidth;

  aol::Vec2<int> minCenter ( Grid.getNumX(), Grid.getNumY() );
  aol::Vec2<int> maxCenter ( 0, 0 );
  for ( int i = 0; i < numImages; ++i ) {
    for ( int j = 0; j < 2; ++j ) {
      minCenter[j] = aol::Min ( minCenter[j], aol::Rint( centers[i][j] ) );
      maxCenter[j] = aol::Max ( maxCenter[j], aol::Rint( centers[i][j] ) );
    }
  }

  for ( int j = 0; j < 2; ++j ) {
    cropWidth[j] = aol::Min ( minCenter[j], Grid.getSize()[j] - 1 - maxCenter[j] );
    CropSize[j] = 2*cropWidth[j];
  }

  for ( int i = 0; i < numImages; ++i ) {
    for ( int j = 0; j < 2; ++j ) {
      CropStart[i][j] = aol::Rint( centers[i][j] - cropWidth[j] );
    }
  }
}

} // end namespace qc

#endif // __PARAMREG_H
