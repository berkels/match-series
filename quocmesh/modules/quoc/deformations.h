#ifndef __DEFORMATIONS_H
#define __DEFORMATIONS_H

#include <generator.h>
#include <Newton.h>
#include <multiArray.h>

namespace qc {

template <typename ConfiguratorType>
void deformImageWithCoarseDeformation ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                                        const typename ConfiguratorType::InitType &Finegrid,
                                        const typename ConfiguratorType::InitType &Coarsegrid,
                                        aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                                        const aol::MultiVector<typename ConfiguratorType::RealType> &Phidofs ) {
  typedef typename ConfiguratorType::RealType RealType;

  qc::MultiArray<RealType, ConfiguratorType::Dim> phi ( Coarsegrid, Phidofs );
  typename ConfiguratorType::ArrayType aImage ( Image, Finegrid );
  typename ConfiguratorType::ArrayType aDeformedImage ( DeformedImage, Finegrid );

  const RealType h = Finegrid.H();
  const qc::GridSize<ConfiguratorType::Dim> size ( Finegrid );

  typename ConfiguratorType::InitType::OldAllNodeIterator fnit;
  for ( fnit = Finegrid._nBeginIt; fnit != Finegrid._nEndIt; ++fnit ) {
    RealType scaling =  ( static_cast<RealType> ( Coarsegrid.getWidth() - 1 ) ) / ( static_cast<RealType> ( Finegrid.getWidth() - 1 ) );
    typename ConfiguratorType::VecType ds;
    typename ConfiguratorType::VecType position;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      position[i] = ( *fnit ) [i];
    }
    position *= scaling;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      ds[i] = ( *fnit ) [i] + phi[i].interpolate ( position ) / h;
      ds[i] = aol::Clamp ( ds[i], aol::ZOTrait<RealType>::zero, static_cast<RealType> ( size[i] - 1 ) );
    }
    const RealType value = aImage.interpolate ( ds );
    aDeformedImage.set ( *fnit, value );
  }
  return;


}

template <typename ConfiguratorType>
void deformImageWithCoarseDeformation ( const aol::MultiVector<typename ConfiguratorType::RealType> &Image,
                                        const typename ConfiguratorType::InitType &Finegrid,
                                        const typename ConfiguratorType::InitType &Coarsegrid,
                                        aol::MultiVector<typename ConfiguratorType::RealType> &DeformedImage,
                                        const aol::MultiVector<typename ConfiguratorType::RealType> &Phidofs ) {
  for ( int i = 0; i < Image.numComponents(); ++i )
    deformImageWithCoarseDeformation<ConfiguratorType> ( Image[i], Finegrid, Coarsegrid, DeformedImage[i], Phidofs );
}

/**
 * \brief DeformedImage = Image circ Phi
 *        In case ExtendWithConstant is false, points outside of the image domain are projected back orthogonally (normal constant extension)
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
void DeformImage ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                   const typename ConfiguratorType::InitType &Grid,
                   aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                   const aol::MultiVector<typename ConfiguratorType::RealType> &Phi,
                   const bool ExtendWithConstant = true,
                   const typename ConfiguratorType::RealType ExtensionConstant = 0,
                   const bool NearestNeighborInterpolation = false ) {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::ArrayType imageArray ( Image, Grid, aol::FLAT_COPY );
  typename ConfiguratorType::ArrayType deformedImageArray ( DeformedImage, Grid, aol::FLAT_COPY );
  const qc::MultiArray<RealType, ConfiguratorType::Dim> phiMArray ( Grid, Phi, aol::FLAT_COPY );

  const aol::Vec3<int> gridSize = Grid.getSize();
  const RealType h = static_cast<RealType> ( Grid.H() );

  typename ConfiguratorType::InitType::OldAllNodeIterator fnit;
  for ( fnit = Grid._nBeginIt; fnit != Grid._nEndIt; ++fnit ) {
    typename ConfiguratorType::VecType ds;

    bool transformPositionInDomain = true;

    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      ds[i] = ( *fnit ) [i] + phiMArray[i].get ( *fnit ) / h;
      if ( ExtendWithConstant == true ) {
        if ( ds[i] < aol::ZOTrait<RealType>::zero || ds[i] > static_cast<RealType> ( gridSize[i] - 1 ) || aol::isNaN ( ds[i] ) )
          transformPositionInDomain = false;
      } else
        ds[i] = aol::Clamp ( ds[i], aol::ZOTrait<RealType>::zero, static_cast<RealType> ( gridSize[i] - 1 ) );

      if ( NearestNeighborInterpolation )
        ds[i] = aol::Rint ( ds[i] );
    }

    RealType value = ExtensionConstant;
    if ( transformPositionInDomain == true )
      value = imageArray.interpolate ( ds );
    deformedImageArray.set ( *fnit, value );
  }
}

/**
 * This struct + static function construction is a workaround to C++ limitations.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, qc::Dimension Dim>
struct doFastDeformImage {};

template <typename ConfiguratorType>
struct doFastDeformImage<ConfiguratorType, qc::QC_1D> {
  static void apply ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                      const typename ConfiguratorType::InitType &Grid,
                      aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                      const aol::MultiVector<typename ConfiguratorType::RealType> &Phi ) {


    static bool doFastDeformImage1DWarningPrinted = false;
    if ( !doFastDeformImage1DWarningPrinted ) {
      cerr << aol::color::error << "qc::doFastDeformImage<ConfiguratorType, qc::QC_1D> has not been tested yet! Please don't remove this warning before testing it properly!" << aol::color::reset << endl;;
      doFastDeformImage1DWarningPrinted = true;
    }

    typedef typename ConfiguratorType::RealType RealType;
    const typename ConfiguratorType::ArrayType imageArray ( Image, Grid, aol::FLAT_COPY );

    const aol::Vec3<int> gridSize = Grid.getSize();
    const RealType h = static_cast<RealType> ( Grid.H() );

    for ( int i = 0; i < gridSize[0]; ++i ) {
      typename ConfiguratorType::VecType ds ( i + Phi[0][i] / h );
      for ( int c = 0; c < ConfiguratorType::Dim; c++ )
        ds[c] = aol::Clamp ( ds[c], aol::ZOTrait<RealType>::zero, static_cast<RealType> ( gridSize[c] - 1 ) );
      DeformedImage[i] = imageArray.interpolate ( ds );
    }
  }
};

template <typename ConfiguratorType>
struct doFastDeformImage<ConfiguratorType, qc::QC_2D> {
  static void apply ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                      const typename ConfiguratorType::InitType &Grid,
                      aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                      const aol::MultiVector<typename ConfiguratorType::RealType> &Phi ) {
    typedef typename ConfiguratorType::RealType RealType;
    const typename ConfiguratorType::ArrayType imageArray ( Image, Grid, aol::FLAT_COPY );

    const aol::Vec3<int> gridSize = Grid.getSize();
    const RealType h = static_cast<RealType> ( Grid.H() );

    for ( int j = 0; j < gridSize[1]; ++j ) {
      for ( int i = 0; i < gridSize[0]; ++i ) {
        const int index = qc::ILexCombine2( i, j, Grid.getNumX() );
        typename ConfiguratorType::VecType ds ( i + Phi[0][index] / h, j + Phi[1][index] / h );
        for ( int c = 0; c < ConfiguratorType::Dim; c++ )
          ds[c] = aol::Clamp ( ds[c], aol::ZOTrait<RealType>::zero, static_cast<RealType> ( gridSize[c] - 1 ) );
        DeformedImage[index] = imageArray.interpolate ( ds );
      }
    }
  }
};

template <typename ConfiguratorType>
struct doFastDeformImage<ConfiguratorType, qc::QC_3D> {
  static void apply ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                      const typename ConfiguratorType::InitType &Grid,
                      aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                      const aol::MultiVector<typename ConfiguratorType::RealType> &Phi ) {
    typedef typename ConfiguratorType::RealType RealType;
    const typename ConfiguratorType::ArrayType imageArray ( Image, Grid, aol::FLAT_COPY );

    const aol::Vec3<int> gridSize = Grid.getSize();
    const RealType h = static_cast<RealType> ( Grid.H() );

    for ( int k = 0; k < gridSize[2]; ++k ) {
      for ( int j = 0; j < gridSize[1]; ++j ) {
        for ( int i = 0; i < gridSize[0]; ++i ) {
          const int index = qc::ILexCombine3( i, j, k, Grid.getNumX(), Grid.getNumY() );
          typename ConfiguratorType::VecType ds ( i + Phi[0][index] / h, j + Phi[1][index] / h, k + Phi[2][index] / h  );
          for ( int c = 0; c < ConfiguratorType::Dim; c++ )
            ds[c] = aol::Clamp ( ds[c], aol::ZOTrait<RealType>::zero, static_cast<RealType> ( gridSize[c] - 1 ) );
          DeformedImage[index] = imageArray.interpolate ( ds );
        }
      }
    }
  }
};


/**
 * \brief Specialized version of qc:DeformImage that is less flexible, but noticeably faster.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, qc::Dimension Dim>
void FastDeformImage ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                       const typename ConfiguratorType::InitType &Grid,
                       aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                       const aol::MultiVector<typename ConfiguratorType::RealType> &Phi ) {
  doFastDeformImage<ConfiguratorType, Dim>::apply ( Image, Grid, DeformedImage, Phi );
}

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
void DeformImage ( const aol::MultiVector<typename ConfiguratorType::RealType> &Image,
                   const typename ConfiguratorType::InitType &Grid,
                   aol::MultiVector<typename ConfiguratorType::RealType> &DeformedImage,
                   const aol::MultiVector<typename ConfiguratorType::RealType> &Phi,
                   const bool ExtendWithConstant = true,
                   const typename ConfiguratorType::RealType ExtensionConstant = 0,
                   const bool NearestNeighborInterpolation = false ) {
  for ( int i = 0; i < Image.numComponents(); ++i )
    qc::DeformImage<ConfiguratorType>( Image[i], Grid, DeformedImage[i], Phi, ExtendWithConstant, ExtensionConstant, NearestNeighborInterpolation );
}

/**
 * @brief DeformedImage = Image circ Phi,
 * where deformed image is extended by a second image in regions, where the deformed image is not defined.
 *
 * @author Wirth
 */
template <typename ConfiguratorType>
void DeformImage ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                   const typename ConfiguratorType::InitType &Grid,
                   aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                   const aol::MultiVector<typename ConfiguratorType::RealType> &Phi,
                   const aol::Vector<typename ConfiguratorType::RealType> &ExtendImage ) {
  typedef typename ConfiguratorType::RealType RealType;
  typename ConfiguratorType::ArrayType imageArray ( Image, Grid, aol::FLAT_COPY );
  typename ConfiguratorType::ArrayType deformedImageArray ( DeformedImage, Grid, aol::FLAT_COPY );
  typename ConfiguratorType::ArrayType extendImageArray ( ExtendImage, Grid, aol::FLAT_COPY );

  const aol::Vec3<int> gridSize = Grid.getSize();
  const RealType h = Grid.H();

  qc::FastILexMapper<ConfiguratorType::Dim> mapper ( Grid );
  for ( qc::RectangularIterator<ConfiguratorType::Dim> fnit ( Grid ); fnit.notAtEnd(); ++fnit ) {
    typename ConfiguratorType::VecType ds;

    bool transformPositionInDomain = true;

    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      ds[i] = ( *fnit ) [i] + Phi[i].get ( mapper.getGlobalIndex ( *fnit ) ) / h;
      if ( ds[i] < aol::ZOTrait<RealType>::zero || ds[i] > static_cast<RealType> ( gridSize[i] - 1 ) )
        transformPositionInDomain = false;
    }
    if ( transformPositionInDomain == true )
      deformedImageArray.set ( *fnit, imageArray.interpolate ( ds ) );
    else
      deformedImageArray.set ( *fnit, extendImageArray.get ( *fnit ) );
  }
  return;
}

template <typename ConfiguratorType>
void DeformImageFromMVecPart ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                               const typename ConfiguratorType::InitType &Grid,
                               aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                               const aol::MultiVector<typename ConfiguratorType::RealType> &MultiVec,
                               const bool ExtendWithZero,
                               const int StartIndex ) {
  aol::MultiVector<typename ConfiguratorType::RealType> phi ( 0, 0 );
  for ( int i = StartIndex; i < StartIndex + ConfiguratorType::Dim; i++ )
    phi.appendReference ( MultiVec[i] );
  DeformImage<ConfiguratorType> ( Image, Grid, DeformedImage, phi, ExtendWithZero );
}

/**
 * @brief DeformedImage = Image\f$\circ\phi^{-1}\f$, where \f$\phi\f$=D+identity.
 * "VectorType" may be either aol::Vector or aol::MultiVector.
 *
 * @author Wirth
 */
template <typename ConfiguratorType, typename VectorType>
void InvDeformImage ( const VectorType &Image,
                      const typename ConfiguratorType::InitType &Grid,
                      VectorType &DeformedImage,
                      const aol::MultiVector<typename ConfiguratorType::RealType> &D,
                      // VC++ 2008 somehow doesn't see here that ConfiguratorType::Dim is of type qc::Dimension,
                      // thus we have to use qc::BitArrayTrait which can use a simple int.
                      const typename qc::BitArrayTrait<ConfiguratorType::Dim>::ArrayType *Mask = NULL ) {

  // deform the image
  qc::TransformFunction<typename ConfiguratorType::RealType,ConfiguratorType::Dim> transformFunction( Grid, Mask );
  transformFunction.setDeformation( D );
  transformFunction.apply( Image, DeformedImage );
}

/**
 * @brief DeformedImage = Image\f$\circ\phi^{-1}\f$, where \f$\phi\f$=D+identity.
 * The deformed image is extended by "ExtendImage", where it is not defined.
 * "VectorType" may be either aol::Vector or aol::MultiVector.
 *
 * @author Wirth
 */
template <typename ConfiguratorType, typename VectorType>
void InvDeformImage ( const VectorType &Image,
                      const typename ConfiguratorType::InitType &Grid,
                      VectorType &DeformedImage,
                      const aol::MultiVector<typename ConfiguratorType::RealType> &D,
                      const VectorType &ExtendImage ) {

  // deform the image
  qc::TransformFunction<typename ConfiguratorType::RealType,ConfiguratorType::Dim> transformFunction( Grid );
  transformFunction.setDeformation( D );
  typename qc::BitArray<ConfiguratorType::Dim> mask;
  transformFunction.transform( Image, DeformedImage, mask, ExtendImage );
}

/**
 * @brief DeformedImage = Image\f$\circ\phi^{-1}\f$, where \f$\phi\f$=D+identity.
 * The deformed image is extended by "ExtendImage", where it is not defined.
 * "VectorType" may be either aol::Vector or aol::MultiVector.
 *
 * @author Wirth
 */
template <typename ConfiguratorType, typename VectorType>
void InvDeformImage ( const VectorType &Image,
                      const typename ConfiguratorType::InitType &Grid,
                      VectorType &DeformedImage,
                      const aol::MultiVector<typename ConfiguratorType::RealType> &D,
                      const VectorType &ExtendImage,
                      typename qc::BitArray<ConfiguratorType::Dim> &InvDispDefined ) {

  // deform the image
  qc::TransformFunction<typename ConfiguratorType::RealType,ConfiguratorType::Dim> transformFunction( Grid );
  transformFunction.setDeformation( D );
  transformFunction.transform( Image, DeformedImage, InvDispDefined, ExtendImage );
}


/**
 * \brief Inverse = Phi^{-1}
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
void approxInverseDeformation ( const typename ConfiguratorType::InitType &Grid,
                                const aol::MultiVector<typename ConfiguratorType::RealType> &Phi,
                                aol::MultiVector<typename ConfiguratorType::RealType> &Inverse ) {
  aol::MultiVector<typename ConfiguratorType::RealType> identity( Grid );
  qc::DataGenerator<ConfiguratorType>( Grid ).generateIdentity( identity );
  // At least VC++ doesn't seem to be able to discern the different qc::InvDeformImage versions, so we can't call it here.
  qc::TransformFunction<typename ConfiguratorType::RealType,ConfiguratorType::Dim> transformFunction( Grid );
  transformFunction.setDeformation( Phi );
  qc::BitArray<ConfiguratorType::Dim> valuesSet( GridSize<ConfiguratorType::Dim>::createFrom ( Grid ) );
  transformFunction.transform( identity, Inverse, valuesSet );
  Inverse -= identity;
  // Try to set the values not covered by qc::TransformFunction to something reasonable.
  for ( int i = 0; i < ConfiguratorType::Dim; ++i )
    Inverse[i].setAllMasked ( 0.5 * ( Inverse[i].getMinValueMasked ( valuesSet ) + Inverse[i].getMaxValueMasked ( valuesSet ) ), valuesSet, true );
}

/**
 * @brief Copies Image to ExtendedImage, but replaces values at all nodes where ValueDefined is false such that L2-norm of gradient is minimized.
 * "VectorType" may be either aol::Vector or aol::MultiVector.
 *
 * @author Wirth
 */
template <typename ConfiguratorType, typename VectorType>
void SmoothlyExtendImage ( const VectorType &Image,
                           const typename ConfiguratorType::InitType &Grid,
                           VectorType &ExtendedImage,
                           // VC++ 2008 somehow doesn't see here that ConfiguratorType::Dim is of type qc::Dimension,
                           // thus we have to use qc::BitArrayTrait which can use a simple int.
                           const typename qc::BitArrayTrait<ConfiguratorType::Dim>::ArrayType &ValueDefined ) {
  typedef typename ConfiguratorType::RealType RealType;

  aol::MultiVector<RealType> image, extendedImage;
  image.appendReference( Image );
  extendedImage.appendReference( ExtendedImage );

  // smoothly extend Image $u$ to ExtendedImage $v$ by minimizing $\int_\Omega|\nabla v|^2dx$
  // (i.e.solving $\Delta v=0$) under boundary conditions $v=u$ where ValueDefined = true

  // compute boundary conditions
  aol::MultiVector<RealType> bc( image, aol::DEEP_COPY );
  for ( int i = 0; i < bc.numComponents(); i++ )
    for ( int j = 0; j < ValueDefined.size(); j++ )
      if ( !ValueDefined[j] )
        bc[i][j] = 0.;
  // compute right hand side
  aol::StiffOp<ConfiguratorType> stiffOp( Grid );
  aol::MultiVector<RealType> rhs( image, aol::STRUCT_COPY );
  for ( int i = 0; i < rhs.numComponents(); i++ ){
    stiffOp.apply( bc[i], rhs[i] );
    rhs[i] *= -1.;
    // build in boundary conditions
    for ( int j = 0; j < ValueDefined.size(); j++ )
      if ( ValueDefined[j] )
        rhs[i][j] = bc[i][j];
  }
  // find the solution
  typename ConfiguratorType::MatrixType systemMatrix( Grid );
  stiffOp.assembleAddMatrix( systemMatrix, &ValueDefined, true );
  aol::DiagonalPreconditioner<aol::Vector<RealType> > preCond( systemMatrix );
  aol::PCGInverse<aol::Vector<RealType> > inv( systemMatrix, preCond, 1.e-16, 1000, aol::STOPPING_ABSOLUTE );
  inv.setQuietMode( true );
  for ( int j = 0; j < rhs.numComponents(); j++ )
    inv.apply( rhs[j], extendedImage[j] );
}

/**
 * @brief DeformedImage = \f$L(\f$Image\f$\circ\phi^{-1})\f$, where \f$\phi\f$=D+identity.
 * The deformed image is smoothly extended via the operator \f$L\f$, where the inverse deformation is not defined.
 * \f$L\f$ minimizes the L2-norm of the gradient.
 * "VectorType" may be either aol::Vector or aol::MultiVector.
 *
 * @author Wirth
 */
template <typename ConfiguratorType, typename VectorType>
void InvDeformAndSmoothlyExtendImage ( const VectorType &Image,
                                       const typename ConfiguratorType::InitType &Grid,
                                       VectorType &DeformedImage,
                                       const aol::MultiVector<typename ConfiguratorType::RealType> &D ) {

  // deform the image
  typename qc::BitArray<ConfiguratorType::Dim> invDispDefined;
  VectorType deformationResult( Image, aol::STRUCT_COPY );
  qc::TransformFunction<typename ConfiguratorType::RealType,ConfiguratorType::Dim> transformFunction( Grid );
  transformFunction.setDeformation( D );
  transformFunction.transform( Image, deformationResult, invDispDefined );
  // smoothly extend the image, where the inverse deformation was not defined
  qc::SmoothlyExtendImage<ConfiguratorType,VectorType>( deformationResult, Grid, DeformedImage, invDispDefined );
}

/**
 * @brief identity + ResultDisplacement = \f$L(\phi_1\circ\phi_2^{-1})\f$, where \f$\phi_i\f$=Displacementi+identity.
 * The resulting displacement is smoothly extended via the operator \f$L\f$ where the displacement is not defined.
 * \f$L\f$ minimizes the L2-norm of the gradient.
 *
 * @author Wirth
 */
template <typename ConfiguratorType>
void InvConcatAndSmoothlyExtendDeformations ( const aol::MultiVector<typename ConfiguratorType::RealType> &Displacement1,
                                              const aol::MultiVector<typename ConfiguratorType::RealType> &Displacement2,
                                              const typename ConfiguratorType::InitType &Grid,
                                              aol::MultiVector<typename ConfiguratorType::RealType> &ResultDisplacement ) {
  // compute the deformation \phi_1
  aol::MultiVector<typename ConfiguratorType::RealType> identity( Grid ), deformation1( Displacement1 ), resultDisplacement( Grid );
  qc::DataGenerator<ConfiguratorType>( Grid ).generateIdentity( identity );
  deformation1 += identity;

  // concatenate \phi_1 with \phi_2^{-1}
  typename qc::BitArray<ConfiguratorType::Dim> invDispDefined;
  qc::TransformFunction<typename ConfiguratorType::RealType,ConfiguratorType::Dim> transformFunction( Grid );
  transformFunction.setDeformation( Displacement2 );
  transformFunction.transform( deformation1, resultDisplacement, invDispDefined );

  // compute the corresponding displacement
  resultDisplacement -= identity;

  // smoothly extend the displacement where it is not yet defined
  qc::SmoothlyExtendImage<ConfiguratorType,aol::MultiVector<typename ConfiguratorType::RealType> >( resultDisplacement, Grid, ResultDisplacement, invDispDefined );
}

/**
 * @brief DeformedImage = \f$L(\f$Image\f$\circ\phi)\f$, where \f$\phi\f$=D+identity.
 * The deformed image is smoothly extended via the operator \f$L\f$, where the image is not defined.
 * \f$L\f$ minimizes the L2-norm of the gradient.
 * "VectorType" may be either aol::Vector or aol::MultiVector.
 *
 * @author Wirth
 */
template <typename ConfiguratorType, typename VectorType>
void DeformAndSmoothlyExtendImage ( const VectorType &Image,
                                    const typename ConfiguratorType::InitType &Grid,
                                    VectorType &DeformedImage,
                                    const aol::MultiVector<typename ConfiguratorType::RealType> &D ) {
  aol::MultiVector<typename ConfiguratorType::RealType> image, deformedImage;
  image.appendReference( Image );
  deformedImage.appendReference( DeformedImage );

  // deform the image
  aol::MultiVector<typename ConfiguratorType::RealType> deformationResult( image, aol::STRUCT_COPY );
  aol::Vector<typename ConfiguratorType::RealType> extendImage( Grid.getNumberOfNodes() );
  extendImage.setAll( aol::NumberTrait<typename ConfiguratorType::RealType>::NaN );
  for ( int i = 0; i < image.numComponents(); i++ )
    DeformImage<ConfiguratorType>( image[i], Grid, deformationResult[i], D, extendImage );
  typename qc::BitArray<ConfiguratorType::Dim> imageDefined( Grid );
  for ( int i = 0; i < imageDefined.size(); i++ )
    imageDefined.set( i, !aol::isNaN( deformationResult[0][i] ) );

  // smoothly extend the image where it is not yet defined
  qc::SmoothlyExtendImage<ConfiguratorType,aol::MultiVector<typename ConfiguratorType::RealType> >( deformationResult, Grid, deformedImage, imageDefined );
}

/**
 * @brief identity + ResultDisplacement = \f$L(\phi_1\circ\phi_2)\f$, where \f$\phi_i\f$=Displacementi+identity.
 * The resulting displacement is smoothly extended via the operator \f$L\f$ where the displacement is not defined.
 * \f$L\f$ minimizes the L2-norm of the gradient.
 *
 * @author Wirth
 */
template <typename ConfiguratorType>
void ConcatAndSmoothlyExtendDeformations ( const aol::MultiVector<typename ConfiguratorType::RealType> &Displacement1,
                                           const aol::MultiVector<typename ConfiguratorType::RealType> &Displacement2,
                                           const typename ConfiguratorType::InitType &Grid,
                                           aol::MultiVector<typename ConfiguratorType::RealType> &ResultDisplacement ) {
  aol::MultiVector<typename ConfiguratorType::RealType> resultDisplacement( Displacement1, aol::STRUCT_COPY );

  // concatenate the displacement with the deformation
  aol::Vector<typename ConfiguratorType::RealType> extendImage( Grid.getNumberOfNodes() );
  extendImage.setAll( aol::NumberTrait<typename ConfiguratorType::RealType>::NaN );
  for ( int i = 0; i < Displacement1.numComponents(); i++ )
    DeformImage<ConfiguratorType>( Displacement1[i], Grid, resultDisplacement[i], Displacement2, extendImage );

  // compute the correct deformation and then displacement
  resultDisplacement += Displacement2;

  // smoothly extend the displacement where it is not yet defined
  typename qc::BitArray<ConfiguratorType::Dim> imageDefined( qc::GridSize<ConfiguratorType::Dim>::createFrom( Grid ) );
  for ( int i = 0; i < imageDefined.size(); i++ )
    imageDefined.set( i, !aol::isNaN( resultDisplacement[0][i] ) );

  qc::SmoothlyExtendImage<ConfiguratorType,aol::MultiVector<typename ConfiguratorType::RealType> >( resultDisplacement, Grid, ResultDisplacement, imageDefined );
}

/**
 * \note NearestNeighborInterpolation only has an effect if UseInverseDeformation is false.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
void deformAndSaveColoredImage ( const char *InputImageFileName, const qc::MultiArray<typename ConfiguratorType::RealType, qc::QC_2D> &Deformation, const char *OutputFileName, const bool UseInverseDeformation = false, const bool NearestNeighborInterpolation = false ){
  typedef typename ConfiguratorType::RealType RealType;
  qc::MultiArray< RealType, 2, 3 > deformedImage( Deformation[0].getNumX(), Deformation[0].getNumY() );
  if ( aol::fileNameEndsWith ( InputImageFileName, ".pgm" ) ) {
    deformedImage[0].load ( InputImageFileName );
    for ( int i = 1; i < 3; ++i )
      deformedImage[i] = deformedImage[0];
  }
  else
    deformedImage.loadPNG( InputImageFileName );

  typename ConfiguratorType::InitType grid ( deformedImage[0].getSize() );

  qc::ScalarArray<RealType,qc::QC_2D> tmp( deformedImage[0], aol::STRUCT_COPY );
  for ( int i = 0; i < 3; i++ ){
    if ( UseInverseDeformation == false )
      qc::DeformImage<ConfiguratorType>( deformedImage[i], grid, tmp, Deformation, true, 0, NearestNeighborInterpolation );
    else
      qc::InvDeformImage<ConfiguratorType>( deformedImage[i], grid, tmp, Deformation);
    deformedImage[i] = tmp;
  }

  deformedImage.setQuietMode( true );
  deformedImage.savePNG( OutputFileName );
}

/**
 * Loads image (possibly colored) named "InputImageFileName"
 * deforms it with "Identity + (InputFileNameDefX, InputFileNameDefY)"
 * and writes the result to OutputFileName.
 *
 * \note NearestNeighborInterpolation only has an effect if UseInverseDeformation is false.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
void deformAndSaveColoredImage ( const char *InputImageFileName, const char *InputFileNameDefX, const char *InputFileNameDefY, const char *OutputFileName, const bool UseInverseDeformation = false, const bool NearestNeighborInterpolation = false ){
  typedef typename ConfiguratorType::RealType RealType;
  qc::ScalarArray<RealType,qc::QC_2D> defX( InputFileNameDefX );
  qc::ScalarArray<RealType,qc::QC_2D> defY( InputFileNameDefY );
  aol::MultiVector<RealType> deformation ( 0, 0 );
  deformation.appendReference ( defX );
  deformation.appendReference ( defY );

  qc::RectangularGrid<qc::QC_2D> grid ( defX.getSize() );
  qc::MultiArray<RealType, qc::QC_2D> deformationMArray ( grid, deformation, aol::FLAT_COPY );
  deformAndSaveColoredImage<ConfiguratorType> ( InputImageFileName, deformationMArray, OutputFileName, UseInverseDeformation, NearestNeighborInterpolation );
}

template <typename ConfiguratorType>
void concatenateTwoDeformations ( const qc::GridDefinition &Grid,
                                  const aol::MultiVector<typename ConfiguratorType::RealType> &Deformation1,
                                  const aol::MultiVector<typename ConfiguratorType::RealType> &Deformation2,
                                  aol::MultiVector<typename ConfiguratorType::RealType> &Dest ) {
  typedef typename ConfiguratorType::RealType RealType;

  qc::MultiArray<RealType, ConfiguratorType::Dim> deformation2Array ( Grid, Deformation2 );
  qc::MultiArray<RealType, ConfiguratorType::Dim> destArray ( Grid, Dest );

  const int num = Grid.getWidth();
  RealType h = Grid.H();

  qc::GridDefinition::OldFullNodeIterator fnit;
  qc::FastILexMapper<ConfiguratorType::Dim> mapper ( Grid );
  for ( fnit = Grid.begin(); fnit != Grid.end(); ++fnit ) {
    typename ConfiguratorType::VecType ds;

    bool transformPositionInDomain = true;

    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      ds[i] = ( *fnit ) [i] + Deformation1[i].get ( mapper.getGlobalIndex ( *fnit ) ) / h;
      if ( ds[i] < aol::ZOTrait<RealType>::zero || ds[i] > static_cast<RealType> ( num - 1 ) )
        transformPositionInDomain = false;
    }

    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      RealType value = 0.;
      if ( transformPositionInDomain == true ) {
        value = ( deformation2Array[i].interpolate ( ds ) ) + Deformation1[i].get ( mapper.getGlobalIndex ( *fnit ) );
      }
      destArray[i].set ( *fnit, value );
    }
  }
  return;
}

/**
 * Shifts coordinates by Offset and checks whether the transformed coordinates are still in [0,1]^d.
 * ClipCoord == false -> Returns true if in [0,1]^d and false if not.
 * ClipCoord == true -> Clips the transformed coordinates into [0,1[^d and always returns true.
 * Offset contains the offset relative to [0,1]^d, i.e. the pixel offset is Offset / H.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, bool ClipCoord>
inline bool transformCoord ( const typename ConfiguratorType::InitType &Grid,
                             const typename ConfiguratorType::VecType &Coord,
                             const typename ConfiguratorType::VecType &Offset,
                             qc::Element &TransformedEl,
                             typename ConfiguratorType::VecType &TransformedLocalCoord ) {
  typename ConfiguratorType::VecType transformedCoord;
  for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
    transformedCoord[i] = Coord[i] + Offset[i] / Grid.H();

    if ( ClipCoord && ( transformedCoord[i] >= ( (Grid.getSize())[i] - 1 ) ) )
      transformedCoord[i] = ( (Grid.getSize())[i] - 1 ) - std::numeric_limits<typename ConfiguratorType::RealType>::epsilon();

    TransformedEl[i] = static_cast<short> ( transformedCoord[i] );
    TransformedLocalCoord[i] = transformedCoord[i] - TransformedEl[i];

    if ( ( TransformedEl[i] == ( (Grid.getSize())[i] - 1 ) ) && ( TransformedLocalCoord[i] == 0 ) ) {
      TransformedEl[i] = (Grid.getSize())[i] - 2;
      TransformedLocalCoord[i] = 1;
    }

    if ( ( transformedCoord[i] < 0. ) || ( TransformedEl[i] < 0 ) || ( TransformedEl[i] >= ( (Grid.getSize())[i] - 1 ) ) ) {
      if ( ClipCoord == false )
        return false;
      else {
        TransformedEl[i] = aol::Clamp ( static_cast<int> ( TransformedEl[i] ), 0, Grid.getSize()[i] - 2 );
        if ( transformedCoord[i] < 0. )
          transformedCoord[i] = 0;
      }
    }
  }
  return true;
}

template <typename ConfiguratorType, bool ClipCoord>
inline bool transformCoord ( const typename ConfiguratorType::InitType &Grid,
                             const typename ConfiguratorType::ElementType &El,
                             const typename ConfiguratorType::VecType &RefCoord,
                             const typename ConfiguratorType::VecType &Offset,
                             qc::Element &TransformedEl,
                             typename ConfiguratorType::VecType &TransformedLocalCoord ) {
  typename ConfiguratorType::VecType coord;
  for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
    coord[i] = El[i] + RefCoord[i];
  }
  return qc::transformCoord<ConfiguratorType, ClipCoord> ( Grid, coord, Offset, TransformedEl, TransformedLocalCoord );
}

// Dummy function to emulate a default for the second template argument for the function above.
template <typename ConfiguratorType>
inline bool transformCoord ( const typename ConfiguratorType::InitType &Grid,
                             const typename ConfiguratorType::ElementType &El,
                             const typename ConfiguratorType::VecType &RefCoord,
                             const typename ConfiguratorType::VecType &Offset,
                             qc::Element &TransformedEl,
                             typename ConfiguratorType::VecType &TransformedLocalCoord ) {
  return qc::transformCoord<ConfiguratorType, false> ( Grid, El, RefCoord, Offset, TransformedEl, TransformedLocalCoord );
}

template <typename ConfiguratorType, bool ClipCoord>
inline bool transformCoord ( const typename ConfiguratorType::InitType &Grid,
                             const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrTransformation,
                             const typename ConfiguratorType::ElementType &El,
                             const int QuadPoint,
                             const typename ConfiguratorType::VecType &RefCoord,
                             qc::Element &TransformedEl,
                             typename ConfiguratorType::VecType &TransformedLocalCoord ) {
  typename ConfiguratorType::VecType offset;
  for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
    offset[i] = DiscrTransformation[i].evaluateAtQuadPoint ( El, QuadPoint );
  }
  return qc::transformCoord<ConfiguratorType, ClipCoord> ( Grid, El, RefCoord, offset, TransformedEl, TransformedLocalCoord );
}

template <typename ConfiguratorType, bool ClipCoord>
inline bool transformCoord ( const ConfiguratorType &Config,
                             const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrTransformation,
                             const typename ConfiguratorType::ElementType &El,
                             const int QuadPoint,
                             const typename ConfiguratorType::DomVecType &RefCoord,
                             typename ConfiguratorType::ElementType &TransformedEl,
                             typename ConfiguratorType::DomVecType &TransformedLocalCoord ) {
  typename ConfiguratorType::VecType offset;
  for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
    offset[i] = DiscrTransformation[i].evaluateAtQuadPoint ( El, QuadPoint );
  }
  return Config.template transformCoord<ClipCoord> ( El, RefCoord, offset, TransformedEl, TransformedLocalCoord );
}

// Dummy function to emulate a default for the second template argument for the function above.
template <typename ConfiguratorType>
inline bool transformCoord ( const typename ConfiguratorType::InitType &Grid,
                             const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrTransformation,
                             const typename ConfiguratorType::ElementType &El,
                             const int QuadPoint,
                             const typename ConfiguratorType::VecType &RefCoord,
                             qc::Element &TransformedEl,
                             typename ConfiguratorType::VecType &TransformedLocalCoord ) {
  return qc::transformCoord<ConfiguratorType, false> ( Grid, DiscrTransformation, El, QuadPoint, RefCoord, TransformedEl, TransformedLocalCoord );
}

template <typename ConfiguratorType, bool ClipCoord>
inline bool transformCoord ( const typename ConfiguratorType::InitType &Grid,
                             const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> &DiscrTransformation,
                             const typename ConfiguratorType::ElementType &El,
                             const int QuadPoint,
                             const typename ConfiguratorType::VecType &RefCoord,
                             qc::Element &TransformedEl,
                             typename ConfiguratorType::VecType &TransformedLocalCoord ) {
  typename ConfiguratorType::VecType offset;
  DiscrTransformation.evaluateAtQuadPoint ( El, QuadPoint, offset );
  return qc::transformCoord<ConfiguratorType, ClipCoord> ( Grid, El, RefCoord, offset, TransformedEl, TransformedLocalCoord );
}

template <typename ConfiguratorType, bool ClipCoord>
inline bool transformCoord ( const ConfiguratorType &Config,
                             const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> &DiscrTransformation,
                             const typename ConfiguratorType::ElementType &El,
                             const int QuadPoint,
                             const typename ConfiguratorType::DomVecType &RefCoord,
                             typename ConfiguratorType::ElementType &TransformedEl,
                             typename ConfiguratorType::DomVecType &TransformedLocalCoord ) {
  typename ConfiguratorType::VecType offset;
  DiscrTransformation.evaluateAtQuadPoint ( El, QuadPoint, offset );
  return Config.template transformCoord<ClipCoord> ( El, RefCoord, offset, TransformedEl, TransformedLocalCoord );
}

// Dummy function to emulate a default for the second template argument for the function above.
template <typename ConfiguratorType>
inline bool transformCoord ( const typename ConfiguratorType::InitType &Grid,
                             const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> &DiscrTransformation,
                             const typename ConfiguratorType::ElementType &El,
                             const int QuadPoint,
                             const typename ConfiguratorType::VecType &RefCoord,
                             qc::Element &TransformedEl,
                             typename ConfiguratorType::VecType &TransformedLocalCoord ) {
  return qc::transformCoord<ConfiguratorType, false> ( Grid, DiscrTransformation, El, QuadPoint, RefCoord, TransformedEl, TransformedLocalCoord );
}

//! Shifts a point (given by element El and local coordinates RefCoord) by Offset and computes corresponding point (TransformedEl, TransformedLocalCoord).
/** If the point is shifted outside [0,1[^d, its coordinates are clipped to the interval [0,1[, and the method returns false.
 */
template <typename ConfiguratorType>
inline bool transformAndClipCoord ( const ConfiguratorType &Config,
                                    const typename ConfiguratorType::ElementType &El,
                                    const typename ConfiguratorType::DomVecType &RefCoord,
                                    const typename ConfiguratorType::VecType &Offset,
                                    typename ConfiguratorType::ElementType &TransformedEl,
                                    typename ConfiguratorType::DomVecType &TransformedLocalCoord ) {
  bool insideFlag = true;
  typename ConfiguratorType::VecType coord;
  Config.getGlobalCoords( El, RefCoord, coord );
  coord += Offset;
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

//! Shifts a point (given by element El and local coordinates RefCoord) by Offset and computes corresponding point (TransformedEl, TransformedLocalCoord).
/** If the point is shifted outside [0,1[^d, its coordinates are clipped to the interval [0,1[, and CoordinateWithinLimits is false for the clipped coordinates.
 */
template <typename ConfiguratorType>
inline void transformAndClipCoord ( const ConfiguratorType &Config,
                                    const typename ConfiguratorType::ElementType &El,
                                    const typename ConfiguratorType::DomVecType &RefCoord,
                                    const typename ConfiguratorType::VecType &Offset,
                                    typename ConfiguratorType::ElementType &TransformedEl,
                                    typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                                    aol::Vec<ConfiguratorType::Dim,bool> &CoordinateWithinLimits ) {
  CoordinateWithinLimits.setAll( true );
  typename ConfiguratorType::VecType coord;
  Config.getGlobalCoords( El, RefCoord, coord );
  coord += Offset;
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

//! Shifts a point (given by element El and local coordinates RefCoord) by the transformation DiscrTransformation and computes corresponding point (TransformedEl, TransformedLocalCoord).
/** If the point is shifted outside [0,1]^d, its coordinates are clipped to the interval [0,1], and the method returns false.
 */
template <typename ConfiguratorType>
inline bool transformAndClipCoord ( const ConfiguratorType &Config,
                                    const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrTransformation,
                                    const typename ConfiguratorType::ElementType &El,
                                    const int QuadPoint,
                                    const typename ConfiguratorType::DomVecType &RefCoord,
                                    typename ConfiguratorType::ElementType &TransformedEl,
                                    typename ConfiguratorType::DomVecType &TransformedLocalCoord ) {
  typename ConfiguratorType::VecType offset;
  for ( int i = 0; i < ConfiguratorType::Dim; i++ )
    offset[i] = DiscrTransformation[i].evaluateAtQuadPoint ( El, QuadPoint );
  return qc::transformAndClipCoord<ConfiguratorType> ( Config, El, RefCoord, offset, TransformedEl, TransformedLocalCoord );
}

//! Shifts a point (given by element El and local coordinates RefCoord) by the transformation DiscrTransformation and computes corresponding point (TransformedEl, TransformedLocalCoord).
/** If the point is shifted outside [0,1]^d, its coordinates are clipped to the interval [0,1], and CoordinateWithinLimits is false for the clipped coordinates.
 */
template <typename ConfiguratorType>
inline void transformAndClipCoord ( const ConfiguratorType &Config,
                                    const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrTransformation,
                                    const typename ConfiguratorType::ElementType &El,
                                    const int QuadPoint,
                                    const typename ConfiguratorType::DomVecType &RefCoord,
                                    typename ConfiguratorType::ElementType &TransformedEl,
                                    typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                                    aol::Vec<ConfiguratorType::Dim,bool> &CoordinateWithinLimits ) {
  typename ConfiguratorType::VecType offset;
  for ( int i = 0; i < ConfiguratorType::Dim; i++ )
    offset[i] = DiscrTransformation[i].evaluateAtQuadPoint ( El, QuadPoint );
  qc::transformAndClipCoord<ConfiguratorType> ( Config, El, RefCoord, offset, TransformedEl, TransformedLocalCoord, CoordinateWithinLimits );
}

//! Shifts a point (given by element El and local coordinates RefCoord) by the transformation DiscrTransformation and computes corresponding point (TransformedEl, TransformedLocalCoord).
/** If the point is shifted outside [0,1]^d, its coordinates are clipped to the interval [0,1], and the method returns false.
 */
template <typename ConfiguratorType>
inline bool transformAndClipCoord ( const ConfiguratorType &Config,
                                    const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> &DiscrTransformation,
                                    const typename ConfiguratorType::ElementType &El,
                                    const int QuadPoint,
                                    const typename ConfiguratorType::DomVecType &RefCoord,
                                    typename ConfiguratorType::ElementType &TransformedEl,
                                    typename ConfiguratorType::DomVecType &TransformedLocalCoord ) {
  typename ConfiguratorType::VecType offset;
  DiscrTransformation.evaluateAtQuadPoint ( El, QuadPoint, offset );
  return qc::transformAndClipCoord<ConfiguratorType> ( Config, El, RefCoord, offset, TransformedEl, TransformedLocalCoord );
}

//! Shifts a point (given by element El and local coordinates RefCoord) by the transformation DiscrTransformation and computes corresponding point (TransformedEl, TransformedLocalCoord).
/** If the point is shifted outside [0,1]^d, its coordinates are clipped to the interval [0,1], and CoordinateWithinLimits is false for the clipped coordinates.
 */
template <typename ConfiguratorType>
inline void transformAndClipCoord ( const ConfiguratorType &Config,
                                    const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> &DiscrTransformation,
                                    const typename ConfiguratorType::ElementType &El,
                                    const int QuadPoint,
                                    const typename ConfiguratorType::DomVecType &RefCoord,
                                    typename ConfiguratorType::ElementType &TransformedEl,
                                    typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                                    aol::Vec<ConfiguratorType::Dim,bool> &CoordinateWithinLimits ) {
  typename ConfiguratorType::VecType offset;
  DiscrTransformation.evaluateAtQuadPoint ( El, QuadPoint, offset );
  qc::transformAndClipCoord<ConfiguratorType> ( Config, El, RefCoord, offset, TransformedEl, TransformedLocalCoord, CoordinateWithinLimits );
}

//! Shifts coordinates by TranslationOffset and then deforms the shifted coordinates with DiscrTransformation.
//! Returns true if the transformed coordinates are still in [0,1]^d after each step and false if not.
//! TranslationOffset contains the offset relative to [0,1]^d, i.e. the pixel offset is Offset / H
template <typename ConfiguratorType>
inline bool translateAndTransformCoord ( const typename ConfiguratorType::InitType &Grid,
                                         const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrTransformation,
                                         const typename ConfiguratorType::ElementType &El,
                                         const typename ConfiguratorType::DomVecType &RefCoord,
                                         const typename ConfiguratorType::VecType &TranslationOffset,
                                         typename ConfiguratorType::ElementType &TransformedEl,
                                         typename ConfiguratorType::DomVecType &TransformedLocalCoord ) {
  typename ConfiguratorType::VecType translated_local_coord;
  qc::Element translated_el;
  // translation by TranslationOffset -> store result in translated_el, translated_local_coord
  if ( qc::transformCoord<ConfiguratorType> ( Grid, El, RefCoord, TranslationOffset, translated_el, translated_local_coord ) ) {
    // deformation by DiscreteFunctionDefault
    typename ConfiguratorType::VecType offset;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      offset[i] = DiscrTransformation[i].evaluate ( translated_el, translated_local_coord );
    }
    return qc::transformCoord<ConfiguratorType> ( Grid, translated_el, translated_local_coord, offset, TransformedEl, TransformedLocalCoord );
  } else
    return false;
}

/**
 *  @brief Save the MultiVector to file
 *  @author Jingfeng Han
**/
template<typename ConfiguratorType>
bool SaveMultiVector ( const qc::GridDefinition &Grid,
                       const aol::MultiVector<typename ConfiguratorType::RealType> &Phi,
                       const char *Filename ) {
  typedef typename ConfiguratorType::RealType RealType;
  ofstream of ( Filename, ios::binary );
  if ( !of )
    return false;

  int width = Grid.getWidth ();
  int depth = Grid.getGridDepth();
  int NumOfBytes = sizeof ( RealType );
  int NumOfNodes = Grid.getNumberOfNodes();

  if ( NumOfBytes == sizeof ( float ) ) {
    of << "F" << "\n" << width << "\n" << depth << "\n" << NumOfNodes << "\n";
    cerr << "F" << "\n" << width << "\n" << depth << "\n" << NumOfNodes << "\n";
  } else {
    of << "D" << "\n" << width << "\n" << depth << "\n" << NumOfNodes << "\n";
    cerr << "D" << "\n" << width << "\n" << depth << "\n" << NumOfNodes << "\n";
  }
  RealType* dummy = NULL;
  dummy = new RealType[ConfiguratorType::Dim*NumOfNodes];
  if ( !dummy ) {
    return false;
  }


  //ofstream tmpOf("P:/local/tmp/save.txt");

  qc::GridDefinition::OldFullNodeIterator fnit;
  qc::FastILexMapper<ConfiguratorType::Dim> mapper ( Grid );
  int i = 0;
  for ( fnit = Grid.begin(); fnit != Grid.end(); ++fnit ) {

    for ( int n = 0;n < ConfiguratorType::Dim;n++ ) {
      dummy[i++] = Phi[n][mapper.getGlobalIndex ( *fnit ) ];
      //tmpOf<< dummy[i-1]<<" ";
    }

  }

  of.write ( reinterpret_cast<char*> ( dummy ), ConfiguratorType::Dim*NumOfNodes*sizeof ( RealType ) );
  of.close();
  delete[] dummy;

  return true;
}

/**
 *  @brief load the MultiVector from file
 *  @author Jingfeng Han
**/
template<typename ConfiguratorType>
bool LoadMultiVector ( int &LevelOfDeformation,
                       aol::MultiVector<typename ConfiguratorType::RealType>* &Phi,
                       const char *Filename ) {
  typedef typename ConfiguratorType::RealType RealType;
  ifstream in ( Filename, ios::binary );
  if ( !in ) {
    return false;
  }
  char          tmp[256];
  in.get ( tmp, 255, '\n' );
  int NumOfBytes = sizeof ( RealType );
  if ( tmp[0] != 'F' && tmp[0] != 'D' ) {
    return false;
  }
  if ( tmp[0] == 'F' && ( NumOfBytes != sizeof ( float ) ) ) {
    return false;
  }
  if ( tmp[0] == 'D' && ( NumOfBytes != sizeof ( double ) ) ) {
    return false;
  }
  cerr << "type: " << tmp << endl;
  int inWidth = 0;
  int NumOfNodes = 0;

  aol::READ_COMMENTS ( in );
  in >> inWidth;
  cerr << "width: " << inWidth << endl;

  aol::READ_COMMENTS ( in );
  in >> LevelOfDeformation;
  cerr << "inDepth: " << LevelOfDeformation << endl;
  if ( inWidth != static_cast<int> ( pow ( 2., LevelOfDeformation ) + 1. ) ) {
    cerr << "end " << endl;
    return false;
  }

  aol::READ_COMMENTS ( in );
  in >> NumOfNodes;
  cerr << "NumOfNodes: " << NumOfNodes << endl;
  if ( NumOfNodes != static_cast<int> ( pow ( pow ( 2., LevelOfDeformation ) + 1., ConfiguratorType::Dim ) ) ) {
    cerr << "end " << endl;
    return false;
  }
  in.ignore();

  Phi = new aol::MultiVector<RealType> ( ConfiguratorType::Dim, NumOfNodes );

  RealType* dummy = NULL;
  dummy = new RealType[ConfiguratorType::Dim*NumOfNodes];
  in.read ( reinterpret_cast<char*> ( dummy ),  ConfiguratorType::Dim*NumOfNodes * sizeof ( RealType ) );
  in.close();

  qc::GridDefinition::OldFullNodeIterator fnit;
  qc::GridDefinition grid ( LevelOfDeformation, ConfiguratorType::Dim );
  qc::FastILexMapper<ConfiguratorType::Dim> mapper ( grid );
  int i = 0;

  for ( fnit = grid.begin(); fnit != grid.end(); ++fnit ) {

    for ( int n = 0; n < ConfiguratorType::Dim; n++ ) {
      ( *Phi ) [n][mapper.getGlobalIndex ( *fnit ) ] = dummy[i++];
    }

  }

  delete[] dummy;
  return true;
}

/**
 * \brief This class computes via "apply(...)": \f$ \int_\Omega w(x)(d-Ax-b)^2 dx \f$ in "Dest",
 * where the domain \f$ \Omega \f$, the weight \f$ w \f$, and the displacement \f$ d \f$ are passed to the constructor, and
 * \f$ b \f$ and \f$ A \f$ are passed to "apply(...)".
 *
 * \author Wirth
 */
template <typename ConfiguratorType>
class AffineDisplacementProjectionEnergy :
  public aol::Op< aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim+1,typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType VectorType;

  // image grid
  const typename ConfiguratorType::InitType &_grid;
  // auxiliary operator and corresponding matrix
  const aol::WeightedMassOp<ConfiguratorType> _massOp;
  qc::FastUniformGridMatrix<RealType, ConfiguratorType::Dim> _massMatrix;
  // displacement
  const aol::MultiVector<RealType> &_d;

public:
  AffineDisplacementProjectionEnergy( const typename ConfiguratorType::InitType &Grid,
                                      const aol::MultiVector<RealType> &D,
                                      const aol::Vector<RealType> &W ) :
    _grid( Grid ),
    _massOp( Grid, W, aol::ASSEMBLED ),
    _massMatrix( Grid ),
    _d( D ) {
    _massOp.assembleAddMatrix( _massMatrix, 1.0 );
  }

  /**
   * \brief Returns \f$ \int_\Omega w(x)(d-Ax-b)^2 dx \f$ in "Dest",
   * where the domain \f$ \Omega \f$, the weight \f$ w \f$, and the displacement \f$ d \f$ were passed to the constructor, and
   * \f$ b \f$ is the last column of "Arg" and \f$ A \f$ the rest.
   */
  void applyAdd( const aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim+1,RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    // obtain matrix $A+I$ and vector $b$
    typename ConfiguratorType::VecType shift;
    typename ConfiguratorType::MatType transformMatrix;
    Arg.getCol( ConfiguratorType::Dim, shift );
    Arg.getSubMatrix( 0, 0, transformMatrix );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      transformMatrix[i][i] += 1;
    // compute $-(d-Ax-b)$
    aol::MultiVector<RealType> displacement( _d, aol::STRUCT_COPY );
    (qc::DataGenerator<ConfiguratorType>( _grid )).generateAffineDisplacement( transformMatrix, shift, displacement );
    displacement -= _d;
    // return $\int_\Omega w(x)(d-Ax-b)^2 dx$
    aol::MultiVector<RealType> tmp( _d, aol::STRUCT_COPY );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      _massMatrix.apply( displacement[i], tmp[i] );
    Dest += displacement * tmp;
  }
};

/**
 * \brief This class computes via "apply(...)" the gradient of \f$ \int_\Omega w(x)(d-Ax-b)^2 dx \f$ wrt \f$ A \f$ and \f$ b \f$ in "Dest",
 * where the domain \f$ \Omega \f$ and the displacement \f$ d \f$ are passed to the constructor, and
 * \f$ b \f$ and \f$ A \f$ are passed to "apply(...)".
 *
 * \author Wirth
 */
template <typename ConfiguratorType>
class AffineDisplacementProjectionGradient :
  public aol::Op< aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim+1,typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType VectorType;

  // image grid
  const typename ConfiguratorType::InitType &_grid;
  // auxiliary operator and corresponding matrix
  const aol::WeightedMassOp<ConfiguratorType> _massOp;
  qc::FastUniformGridMatrix<RealType, ConfiguratorType::Dim> _massMatrix;
  // displacement
  const aol::MultiVector<RealType> &_d;

public:
  AffineDisplacementProjectionGradient( const typename ConfiguratorType::InitType &Grid,
                                        const aol::MultiVector<RealType> &D,
                                        const aol::Vector<RealType> &W ) :
    _grid( Grid ),
    _massOp( Grid, W, aol::ASSEMBLED ),
    _massMatrix( Grid ),
    _d( D ) {
    _massOp.assembleAddMatrix( _massMatrix, 1.0 );
  }

  /**
   * \brief Returns variation of \f$ \int_\Omega w(x)(d-Ax-b)^2 dx \f$ in "Dest",
   * where the domain \f$ \Omega \f$ and the displacement \f$ d \f$ were passed to the constructor, and
   * \f$ b \f$ is the last column of "Arg" and \f$ A \f$ the rest.
   */
  void applyAdd( const aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim+1,RealType> &Arg, aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim+1,RealType> &Dest ) const {
    // obtain matrix $A+I$ and vector $b$
    typename ConfiguratorType::VecType shift;
    aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim,RealType> transformMatrix;
    Arg.getCol( ConfiguratorType::Dim, shift );
    Arg.getSubMatrix( 0, 0, transformMatrix );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      transformMatrix[i][i] += 1;
    // compute $-(d-Ax-b)$
    aol::MultiVector<RealType> displacement( _d, aol::STRUCT_COPY );
    (qc::DataGenerator<ConfiguratorType>( _grid )).generateAffineDisplacement( transformMatrix, shift, displacement );
    displacement -= _d;
    // return $\int_\Omega -2w(x)(d-Ax-b) dx$ and $\int_\Omega -2w(x)(d-Ax-b) x^T dx$
    aol::Vector<RealType> ones( _d[0].size() );
    aol::MultiVector<RealType> id( _d, aol::STRUCT_COPY );
    ones.setAll( 1 );
    (qc::DataGenerator<ConfiguratorType>( _grid )).generateIdentity( id );
    aol::Vector<RealType> tmp( _d[0].size() );
    _massMatrix.apply( ones, tmp );
    for ( int j = 0; j < ConfiguratorType::Dim; j++ )
      Dest[j][ConfiguratorType::Dim] += tmp * displacement[j];
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      _massMatrix.apply( id[i], tmp );
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        Dest[j][i] += tmp * displacement[j];
    }
    Dest *= 2;
  }
};

/**
 * \brief This class computes via "apply(...)": \f$ \int_\Omega w(x)(\phi-R(\alpha_i)(x-x_0)-b-x_0)^2 dx \f$ in "Dest",
 * where the domain \f$ \Omega \f$, the weight \f$ w \f$, and the displacement \f$ \phi-id \f$ are passed to the constructor, and
 * \f$ b \f$ and the angles \f$ \alpha_i \f$ for the rotation matrix \f$ R(\alpha_i) \f$ are passed to "apply(...)".
 * If the offset \f$ x_0 \f$ was not passed to the constructor it is taken to be zero.
 *
 * \author Wirth
 */
template <typename ConfiguratorType>
class RigidDisplacementProjectionEnergy :
  public aol::Op< aol::Vector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const AffineDisplacementProjectionEnergy<ConfiguratorType> _affineProjector;
  // the offset $x_0$
  const typename ConfiguratorType::VecType _offset;

public:
  RigidDisplacementProjectionEnergy( const typename ConfiguratorType::InitType &Grid,
                                     const aol::MultiVector<RealType> &D,
                                     const aol::Vector<RealType> &W ) :
    _affineProjector( Grid, D, W ) {
  }

public:
  RigidDisplacementProjectionEnergy( const typename ConfiguratorType::InitType &Grid,
                                     const aol::MultiVector<RealType> &D,
                                     const aol::Vector<RealType> &W,
                                     const typename ConfiguratorType::VecType &Offset ) :
    _affineProjector( Grid, D, W ),
    _offset( Offset ) {
  }

  /**
   * \brief Returns \f$ \int_\Omega w(x)(\phi-R(\alpha_i)(x-x_0)-b-x_0)^2 dx \f$ in "Dest",
   * where the domain \f$ \Omega \f$, the weight \f$ w \f$, and the displacement \f$ \phi-id \f$ were passed to the constructor, and
   * \f$ b \f$ and the angles \f$ \alpha_i \f$ for the rotation matrix \f$ R(\alpha_i) \f$ are passed as "Arg".
   * The first dim entries of "Arg" constitute the vector \f$ b \f$, and the rest are rotation angles for rotation around
   * the z-, y- and x- axis (in 3D).
   * If the offset \f$ x_0 \f$ was not passed to the constructor it is taken to be zero.
   */
  void applyAdd( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    // obtain translation b and put it into the last column of a rectangular matrix
    aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim+1,RealType> totalMatrix;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      totalMatrix[i][ConfiguratorType::Dim] = Arg[i];
    // obtain rotation R
    int index1[3] = { 0, 0, 1 };
    int index2[3] = { 1, 2, 2 };
    int numberAngles = Arg.size() - ConfiguratorType::Dim;
    typename ConfiguratorType::MatType transformMatrix, rotationMatrix;
    transformMatrix.setIdentity();
    for ( int i = 0; i < numberAngles; i++ ) {
      RealType alpha = Arg[ConfiguratorType::Dim + i];
      rotationMatrix.setIdentity();
      rotationMatrix[index1[i]][index1[i]] = cos( alpha );
      rotationMatrix[index2[i]][index2[i]] = cos( alpha );
      rotationMatrix[index1[i]][index2[i]] = -sin( alpha );
      rotationMatrix[index2[i]][index1[i]] = sin( alpha );
      transformMatrix *= rotationMatrix;
    }
    // put A=R-I into the first columns of the rectangular matrix
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      transformMatrix[i][i] -= 1;
    totalMatrix.setSubMatrix( 0, 0, transformMatrix );
    // add $x_0-Rx_0$ to the translation
    typename ConfiguratorType::VecType offset;
    transformMatrix.mult( _offset, offset );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      totalMatrix[i][ConfiguratorType::Dim] -= offset[i];
    // compute the energy
    _affineProjector.applyAdd( totalMatrix, Dest );
  }
};

/**
 * \brief This class computes via "apply(...)" the gradient of \f$ \int_\Omega w(x)(\phi-R(\alpha_i)(x-x_0)-b-x_0)^2 dx \f$ wrt \f$ b \f$ and the \f$ \alpha_i \f$ in "Dest",
 * where the domain \f$ \Omega \f$, the weight \f$ w \f$, and the displacement \f$ \phi-id \f$ are passed to the constructor, and
 * \f$ b \f$ and the angles \f$ \alpha_i \f$ for the rotation matrix \f$ R(\alpha_i) \f$ are passed to "apply(...)".
 * If the offset \f$ x_0 \f$ was not passed to the constructor it is taken to be zero.
 *
 * \author Wirth
 */
template <typename ConfiguratorType>
class RigidDisplacementProjectionGradient :
  public aol::Op< aol::Vector<typename ConfiguratorType::RealType>, aol::Vector<typename ConfiguratorType::RealType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  const AffineDisplacementProjectionGradient<ConfiguratorType> _affineProjectorGradient;
  // the offset $x_0$
  const typename ConfiguratorType::VecType _offset;

public:
  RigidDisplacementProjectionGradient( const typename ConfiguratorType::InitType &Grid,
                                       const aol::MultiVector<RealType> &D,
                                       const aol::Vector<RealType> &W ) :
    _affineProjectorGradient( Grid, D, W ) {
  }

  RigidDisplacementProjectionGradient( const typename ConfiguratorType::InitType &Grid,
                                       const aol::MultiVector<RealType> &D,
                                       const aol::Vector<RealType> &W,
                                       const typename ConfiguratorType::VecType &Offset ) :
    _affineProjectorGradient( Grid, D, W ),
    _offset( Offset ) {
  }

  /**
   * \brief Returns gradient of \f$ \int_\Omega w(x)(\phi-R(\alpha_i)(x-x_0)-b-x_0)^2 dx \f$ wrt \f$ b \f$ and the \f$ \alpha_i \f$ in "Dest",
   * where the domain \f$ \Omega \f$, the weight \f$ w \f$, and the displacement \f$ \phi-id \f$ were passed to the constructor, and
   * \f$ b \f$ and the angles \f$ \alpha_i \f$ for the rotation matrix \f$ R(\alpha_i) \f$ are passed as "Arg".
   * The first dim entries of "Arg" constitute the vector \f$ b \f$, and the rest are rotation angles for rotation around
   * the z-, y- and x- axis (in 3D).
   * If the offset \f$ x_0 \f$ was not passed to the constructor it is taken to be zero.
   */
  void applyAdd( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    // obtain translation b and put it into the last column of a rectangular matrix
    aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim+1,RealType> totalMatrix;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      totalMatrix[i][ConfiguratorType::Dim] = Arg[i];
    // obtain rotation R
    int index1[3] = { 0, 0, 1 };
    int index2[3] = { 1, 2, 2 };
    int numberAngles = Arg.size() - ConfiguratorType::Dim;
    typename ConfiguratorType::MatType transformMatrix, derivMatrix;
    typename ConfiguratorType::VecType derivVector;
    std::vector<typename ConfiguratorType::MatType> rotationMatrices;
    rotationMatrices.resize( numberAngles );
    transformMatrix.setIdentity();
    for ( int i = 0; i < numberAngles; i++ ) {
      RealType alpha = Arg[ConfiguratorType::Dim + i];
      rotationMatrices[i].setIdentity();
      rotationMatrices[i][index1[i]][index1[i]] = cos( alpha );
      rotationMatrices[i][index2[i]][index2[i]] = cos( alpha );
      rotationMatrices[i][index1[i]][index2[i]] = -sin( alpha );
      rotationMatrices[i][index2[i]][index1[i]] = sin( alpha );
      transformMatrix *= rotationMatrices[i];
    }
    // put A=R-I into the first columns of the rectangular matrix
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      transformMatrix[i][i] -= 1;
    totalMatrix.setSubMatrix( 0, 0, transformMatrix );
    // add $x_0-Rx_0$ to the translation
    typename ConfiguratorType::VecType offset;
    transformMatrix.mult( _offset, offset );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      totalMatrix[i][ConfiguratorType::Dim] -= offset[i];
    // obtain the derivative of the projection energy wrt A and b
    aol::Mat<ConfiguratorType::Dim,ConfiguratorType::Dim+1,RealType> deriv;
    _affineProjectorGradient.apply( totalMatrix, deriv );
    // put the derivative wrt b into the result vector
    deriv.getCol( ConfiguratorType::Dim, derivVector );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      Dest[i] += derivVector[i];
    // compute the derivative wrt the rotation matrix R
    deriv.getSubMatrix( 0, 0, derivMatrix );
    transformMatrix.makeTensorProduct( derivVector, _offset );
    derivMatrix -= transformMatrix;
    // compute the derivative wrt the rotation angles and put it into the result vector
    for ( int i = 0; i < numberAngles; i++ ) {
      std::vector<typename ConfiguratorType::MatType> derivRotMatrices( rotationMatrices );
      RealType alpha = Arg[ConfiguratorType::Dim + i];
      derivRotMatrices[i].setZero();
      derivRotMatrices[i][index1[i]][index1[i]] = -sin( alpha );
      derivRotMatrices[i][index2[i]][index2[i]] = -sin( alpha );
      derivRotMatrices[i][index1[i]][index2[i]] = -cos( alpha );
      derivRotMatrices[i][index2[i]][index1[i]] = cos( alpha );
      transformMatrix.setIdentity();
      for ( int j = 0; j < numberAngles; j++ )
        transformMatrix *= derivRotMatrices[j];
      Dest[ConfiguratorType::Dim + i] += transformMatrix.dotProduct( derivMatrix );
    }
  }
};

/**
 * \brief Minimizes \f$ \int_\Omega w(x)(\phi-R(x-x_0)-b-x_0)^2 dx \f$ wrt \f$ b \f$ for rotation \f$ R \f$ and translation \f$ b \f$,
 * where the domain \f$ \Omega \f$, the weight \f$ w \f$, the center \f$ x_0 \f$, and the displacement \f$ \phi-id \f$
 * are passed to the constructor.
 *
 * \author Wirth
 */
template <typename ConfiguratorType>
void projectOntoRigidDisplacement ( const aol::MultiVector<typename ConfiguratorType::RealType> &Displacement,
                                    const aol::Vector<typename ConfiguratorType::RealType> &Weight,
                                    const typename ConfiguratorType::InitType &Grid,
                                    const typename ConfiguratorType::VecType &Center,
                                    typename ConfiguratorType::VecType &Translation,
                                    typename ConfiguratorType::MatType &Rotation,
                                    const bool QuietMode = true ) {
  typedef typename ConfiguratorType::RealType RealType;

  // find the best fit rigid body motion
  const int numberAngles = ( ConfiguratorType::Dim == QC_3D ? 3 : 1 );
  RigidDisplacementProjectionEnergy<ConfiguratorType> e( Grid, Displacement, Weight, Center );
  RigidDisplacementProjectionGradient<ConfiguratorType> de( Grid, Displacement, Weight, Center );
  aol::NewtonInfo<RealType> newtonInfo ( 1E-10, 300, 1E-20, 1000, aol::STOPPING_ABSOLUTE );
  aol::QuasiNewtonIteration< RealType, aol::Vector<RealType>, aol::Vector<RealType> > descent( e, de, newtonInfo );
  descent.setQuietMode( QuietMode );
  aol::Vector<RealType> initialX( ConfiguratorType::Dim + numberAngles ), x( ConfiguratorType::Dim + numberAngles );
  descent.apply( initialX, x );

  // obtain matrix $R$ and vector $b$
  int index1[3] = { 0, 0, 1 };
  int index2[3] = { 1, 2, 2 };
  typename ConfiguratorType::MatType rotationMatrix;
  Rotation.setIdentity();
  for ( int i = 0; i < numberAngles; i++ ) {
    RealType alpha = x[ConfiguratorType::Dim + i];
    rotationMatrix.setIdentity();
    rotationMatrix[index1[i]][index1[i]] = cos( alpha );
    rotationMatrix[index2[i]][index2[i]] = cos( alpha );
    rotationMatrix[index1[i]][index2[i]] = -sin( alpha );
    rotationMatrix[index2[i]][index1[i]] = sin( alpha );
    Rotation *= rotationMatrix;
  }
  for ( int i = 0; i < ConfiguratorType::Dim; i++ )
    Translation[i] = x[i];
}

/**
 * \brief This class represents the weighted mass matrix \f$ \left(\int_\Omega w(\phi(x))\varphi_i(x)\varphi_j(x)dx\right)_{ij} \f$.
 *
 * Here, the \f$ \varphi \f$ represent the finite element base functions, and the domain \f$ \Omega \f$ is represented by the grid
 * which is passed to the constructor. Furthermore, \f$ w \f$ and \f$ d \f$ are passed to the constructor as aol::Vector and
 * aol::MultiVector respectively, where \f$ d \f$ is the displacement such that the deformation \f$ \phi \f$ is given by \f$ d \f$+identity.
 *
 * Same usage as aol::MassOp.
 *
 * \author Wirth
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class DeformedWeightMassOp :
  public aol::FELinScalarWeightedMassInterface<ConfiguratorType, DeformedWeightMassOp<ConfiguratorType, IndexMode>, IndexMode> {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  // the displacement $d$ in all Cartesian directions
  const aol::DiscreteVectorFunctionDefault< ConfiguratorType, ConfiguratorType::Dim > _d;
  // the weight $w$
  const aol::DiscreteFunctionDefault<ConfiguratorType> _w;

public:
  DeformedWeightMassOp( const typename ConfiguratorType::InitType &Grid,
                        const aol::Vector<RealType> &W,
                        const aol::MultiVector<RealType> &D,
                        aol::OperatorType OpType = aol::ONTHEFLY ) :
    // initialise the grid
    aol::FELinScalarWeightedMassInterface<ConfiguratorType, DeformedWeightMassOp<ConfiguratorType, IndexMode>, IndexMode>( Grid, OpType ),
    _d( Grid, D ),
    _w( Grid, W ) {}

  /**
   * Returns \f$ w(\phi(x)) \f$ at the point $x$ specified by element, quadrature point, and local coordinates.
   */
  inline RealType getCoeff( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType& RefCoord ) const {
    // compute $\phi(x)$
    typename ConfiguratorType::DomVecType transformedLocalCoord;
    typename ConfiguratorType::ElementType transformedEl;
    if( !qc::transformAndClipCoord<ConfiguratorType> ( *this->_config, _d, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord ) )
      // if point $x$ was displaced outside the grid, return 0
      return 0.;
    // return $w(\phi(x))$
    return _w.evaluate( transformedEl, transformedLocalCoord );
  }
};

/**
 * \brief This class represents the weighted mass matrix \f$ \left(\int_\Omega w(\phi(x))^2\varphi_i(x)\varphi_j(x)dx\right)_{ij} \f$.
 *
 * Here, the \f$ \varphi \f$ represent the finite element base functions, and the domain \f$ \Omega \f$ is represented by the grid
 * which is passed to the constructor. Furthermore, \f$ w \f$ and \f$ d \f$ are passed to the constructor as aol::Vector and
 * aol::MultiVector respectively, where \f$ d \f$ is the displacement such that the deformation \f$ \phi \f$ is given by \f$ d \f$+identity.
 *
 * Same usage as aol::MassOp.
 *
 * \author Wirth
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class SquaredDeformedWeightMassOp :
  public aol::FELinScalarWeightedMassInterface<ConfiguratorType, SquaredDeformedWeightMassOp<ConfiguratorType, IndexMode>, IndexMode> {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  // the displacement $d$ in all Cartesian directions
  const aol::DiscreteVectorFunctionDefault< ConfiguratorType, ConfiguratorType::Dim > _d;
  // the weight $w$ to be squared
  const aol::DiscreteFunctionDefault<ConfiguratorType> _w;

public:
  SquaredDeformedWeightMassOp( const typename ConfiguratorType::InitType &Grid,
                               const aol::Vector<RealType> &W,
                               const aol::MultiVector<RealType> &D,
                               aol::OperatorType OpType = aol::ONTHEFLY ) :
    // initialise the grid
    aol::FELinScalarWeightedMassInterface<ConfiguratorType, SquaredDeformedWeightMassOp<ConfiguratorType, IndexMode>, IndexMode>( Grid, OpType ),
    _d( Grid, D ),
    _w( Grid, W ) {}

  /**
   * Returns \f$ w(\phi(x))^2 \f$ at the point $x$ specified by element, quadrature point, and local coordinates.
   */
  inline RealType getCoeff( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::VecType& RefCoord ) const {
    // compute $\phi(x)$
    typename ConfiguratorType::DomVecType transformedLocalCoord;
    typename ConfiguratorType::ElementType transformedEl;
    if( !qc::transformCoord<ConfiguratorType> ( this->_grid, _d, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord ) )
      // if point $x$ was displaced outside the grid, return 0
      return 0.;
    // return $w(\phi(x))^2$
    return aol::Sqr( _w.evaluate( transformedEl, transformedLocalCoord ) );
  }
};

/**
 * \brief This class represents the weighted mass matrix \f$ \left(\int_\Omega w(\phi^{-1}(x))\varphi_i(x)\varphi_j(x)|det(\nabla\phi^{-1}(x))|dx\right)_{ij} \f$.
 *
 * Here, the \f$ \varphi \f$ represent the finite element base functions, and the domain \f$ \Omega \f$ is represented by the grid
 * which is passed to the constructor. Furthermore, \f$ w \f$ and \f$ d \f$ are passed to the constructor as aol::Vector and
 * aol::MultiVector respectively, where \f$ d \f$ is the displacement such that the deformation \f$ \phi \f$ is given by \f$ d \f$+identity.
 *
 * Same usage as aol::MassOp.
 *
 * \author Wirth
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class InvDeformedWeightMassOp :
  public aol::FELinScalarWeightedMassInterface<ConfiguratorType, InvDeformedWeightMassOp<ConfiguratorType, IndexMode>, IndexMode> {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  // the displacement $d$ in all Cartesian directions
  const aol::DiscreteVectorFunctionDefault< ConfiguratorType, ConfiguratorType::Dim > _d;
  // the weight $w$
  const aol::DiscreteFunctionDefault<ConfiguratorType> _w;
  // the inverse displacement, $\phi^{-1}$-identity, in all Cartesian directions as vector of the dofs and as multilinear interpolation
  aol::MultiVector<RealType> _invDispDOFs;
  const aol::DiscreteVectorFunctionDefault< ConfiguratorType, ConfiguratorType::Dim > _invDisp;
  // array encoding whether $\phi^{-1}$ is defined at the corresponding node
  typename qc::BitArray<ConfiguratorType::Dim> _invPhiDefined;

public:
  InvDeformedWeightMassOp( const typename ConfiguratorType::InitType &Grid,
                           const aol::Vector<RealType> &W,
                           const aol::MultiVector<RealType> &D,
                           aol::OperatorType OpType = aol::ONTHEFLY ) :
    // initialise the grid
    aol::FELinScalarWeightedMassInterface<ConfiguratorType, InvDeformedWeightMassOp<ConfiguratorType, IndexMode>, IndexMode>( Grid, OpType ),
    _d( Grid, D ),
    _w( Grid, W ),
    _invDispDOFs( D, aol::STRUCT_COPY ),
    _invDisp( Grid, _invDispDOFs ),
    _invPhiDefined( Grid ) {
    // generate the identity function
    aol::MultiVector<RealType> identity( D, aol::STRUCT_COPY );
    qc::DataGenerator<ConfiguratorType>( Grid ).generateIdentity( identity );
    // compute $\phi^{-1}$-identity
    qc::TransformFunction<RealType, ConfiguratorType::Dim> transformOp( Grid );
    transformOp.setDeformation( D );
    transformOp.transform( identity, _invDispDOFs, _invPhiDefined );
    _invDispDOFs -= identity;
  }

  /**
   * Returns \f$ w\circ\phi^{-1}|det(\nabla\phi^{-1})| \f$ evaluated at the point specified by the element "El" and quadrature point "QuadPoint".
   */
  inline RealType getCoeff( const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType& RefCoord ) const {
    // if $\phi^{-1}$ is not defined at this position, return 0
    for ( int i = (*this->_config).getNumLocalDofs( El ) - 1; i >= 0; i-- )
      if ( !_invPhiDefined[(*this->_config).localToGlobal( El, i )] )
        return 0.;

    // compute $\phi^{-1}(x)$
    typename ConfiguratorType::DomVecType transformedLocalCoord;
    typename ConfiguratorType::ElementType transformedEl;
    if ( !qc::transformAndClipCoord<ConfiguratorType> ( *this->_config, _invDisp, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord ) )
      return 0.;

    // compute $D\phi(\phi^{-1}(x))$
    typename ConfiguratorType::MatType dPhi;
    _d.evaluateGradient( transformedEl, transformedLocalCoord, dPhi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dPhi[i][i] += 1.;

    // compute $w(\phi^{-1}(x))|det(D\phi^{-1}(x))|$
    return _w.evaluate( transformedEl, transformedLocalCoord ) / aol::Abs( dPhi.det() );
  }
};

/**
 * \brief This class represents the weighted mass matrix \f$ \left(\int_\Omega w(\phi^{-1}(x))^2\varphi_i(x)\varphi_j(x)|det(\nabla\phi^{-1}(x))|dx\right)_{ij} \f$.
 *
 * Here, the \f$ \varphi \f$ represent the finite element base functions, and the domain \f$ \Omega \f$ is represented by the grid
 * which is passed to the constructor. Furthermore, \f$ w \f$ and \f$ d \f$ are passed to the constructor as aol::Vector and
 * aol::MultiVector respectively, where \f$ d \f$ is the displacement such that the deformation \f$ \phi \f$ is given by \f$ d \f$+identity.
 *
 * Same usage as aol::MassOp.
 *
 * \author Wirth
 * \ingroup MatrixFEOp
 */
template <typename ConfiguratorType, aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE>
class SquaredInvDeformedWeightMassOp :
  public aol::FELinScalarWeightedMassInterface<ConfiguratorType, SquaredInvDeformedWeightMassOp<ConfiguratorType, IndexMode>, IndexMode> {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  // the displacement $d$ in all Cartesian directions
  const aol::DiscreteVectorFunctionDefault< ConfiguratorType, ConfiguratorType::Dim > _d;
  // the weight $w$
  const aol::DiscreteFunctionDefault<ConfiguratorType> _w;
  // the inverse displacement, $\phi^{-1}$-identity, in all Cartesian directions as vector of the dofs and as multilinear interpolation
  aol::MultiVector<RealType> _invDispDOFs;
  const aol::DiscreteVectorFunctionDefault< ConfiguratorType, ConfiguratorType::Dim > _invDisp;
  // array encoding whether $\phi^{-1}$ is defined at the corresponding node
  typename qc::BitArray<ConfiguratorType::Dim> _invPhiDefined;

public:
  SquaredInvDeformedWeightMassOp( const typename ConfiguratorType::InitType &Grid,
                                  const aol::Vector<RealType> &W,
                                  const aol::MultiVector<RealType> &D,
                                  aol::OperatorType OpType = aol::ONTHEFLY ) :
    // initialise the grid
    aol::FELinScalarWeightedMassInterface<ConfiguratorType, SquaredInvDeformedWeightMassOp<ConfiguratorType, IndexMode>, IndexMode>( Grid, OpType ),
    _d( Grid, D ),
    _w( Grid, W ),
    _invDispDOFs( D, aol::STRUCT_COPY ),
    _invDisp( Grid, _invDispDOFs ),
    _invPhiDefined( Grid ) {
    // generate the identity function
    aol::MultiVector<RealType> identity( D, aol::STRUCT_COPY );
    qc::DataGenerator<ConfiguratorType>( Grid ).generateIdentity( identity );
    // compute $\phi^{-1}$-identity
    qc::TransformFunction<RealType, ConfiguratorType::Dim> transformOp( Grid );
    transformOp.setDeformation( D );
    transformOp.transform( identity, _invDispDOFs, _invPhiDefined );
    _invDispDOFs -= identity;
  }

  /**
   * Returns \f$ (w\circ\phi^{-1})^2|det(\nabla\phi^{-1})| \f$ evaluated at the point specified by the element "El" and quadrature point "QuadPoint".
   */
  inline RealType getCoeff( const qc::Element &El, int QuadPoint, const typename ConfiguratorType::VecType& RefCoord ) const {
    // if $\phi^{-1}$ is not defined at this position, return 0
    if ( !_invPhiDefined.elementTrue( El ) )
      return 0.;

    // compute $\phi^{-1}(x)$
    typename ConfiguratorType::DomVecType transformedLocalCoord;
    typename ConfiguratorType::ElementType transformedEl;
    if ( !qc::transformCoord<ConfiguratorType> ( this->_grid, _invDisp, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord ) )
      return 0.;

    // compute $D\phi(\phi^{-1}(x))$
    typename ConfiguratorType::MatType dPhi;
    _d.evaluateGradient( transformedEl, transformedLocalCoord, dPhi );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      dPhi[i][i] += 1.;

    // compute $w(\phi^{-1}(x))|det(D\phi^{-1}(x))|$
    return aol::Sqr( _w.evaluate( transformedEl, transformedLocalCoord ) / aol::Abs( dPhi.det() ) );
  }
};

//! Interface to compute \f$(\int_\Omega f(\phi,\Phi(x),x) \varphi_i\circ\Phi(x) dx)_i\f$,
/** where \f$\Phi\f$ is a deformation on \f$\Omega\f$, \f$\varphi_i\f$ is the \f$i\f$th FE basis function,
 *  and \f$\phi\f$ is the argument of the operator. If \f$\Phi(x)\f$ comes to lie outside \f$\Omega\f$,
 *  it is projected onto the nearest point in \f$\partial\Omega\f$.
 *
 *  \author Wirth
 */
template <typename ConfiguratorType, typename Imp, bool ClipCoord = true>
class FENonlinDeformOpInterface :
  public aol::FENonlinOpInterfaceWithMovedTestFunctions<ConfiguratorType, FENonlinDeformOpInterface<ConfiguratorType, Imp, ClipCoord> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  // the deformation \Phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement;

public:
  explicit FENonlinDeformOpInterface( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &Displacement ) :
    aol::FENonlinOpInterfaceWithMovedTestFunctions<ConfiguratorType, FENonlinDeformOpInterface<ConfiguratorType, Imp, ClipCoord> >( Grid ),
    _displacement( Grid, Displacement ) {}

  //! Represents \f$f(\phi,\Phi(x),x)\f$ and has to be implemented in the derived class.
  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &LocalCoord,
                         const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                         typename ConfiguratorType::RealType &NL ) const {
    this->asImp().getNonlinearity ( DiscFunc, El, QuadPoint, LocalCoord, TransformedEl, TransformedLocalCoord, NL );
  }

  //! Returns \f$\Phi(x)\f$.
  bool getTransformedPosition ( const aol::DiscreteFunctionDefault<ConfiguratorType> &/*DiscFunc*/,
                                const typename ConfiguratorType::ElementType &El,
                                int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord,
                                typename ConfiguratorType::ElementType &TransformedEl,
                                typename ConfiguratorType::DomVecType &TransformedLocalCoord ) const {
    // It should be possible to use qc::transformCoord for both cases. This construction
    // is just kept for the moment to guarantee perfect backwards compatibility.
    if ( ClipCoord ) {
      qc::transformAndClipCoord<ConfiguratorType>( *this->_config, _displacement, El, QuadPoint, RefCoord, TransformedEl, TransformedLocalCoord );
      return true;

    }
    else
      return qc::transformCoord<ConfiguratorType, ClipCoord>( this->getConfigurator(), _displacement, El, QuadPoint, RefCoord, TransformedEl, TransformedLocalCoord );
  }

protected:
  // Barton-Nackman trick
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

//! Interface to compute \f$\int_\Omega f(\phi,\Phi(x),x)dx\f$,
/** where \f$\Phi\f$ is a deformation on \f$\Omega\f$ and \f$\phi\f$ is the argument of the operator.
 *  If \f$\Phi(x)\f$ comes to lie outside \f$\Omega\f$,
 *  it is projected onto the nearest point in \f$\partial\Omega\f$.
 *
 *  \author Wirth
 */
template <typename ConfiguratorType, typename Imp>
class FENonlinDeformIntegrationInterface
      : public aol::Op<aol::Vector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  aol::DeleteFlagPointer<const ConfiguratorType> _config;
  // the deformation \Phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement;

public:
  explicit FENonlinDeformIntegrationInterface ( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &Displacement ) :
    _config ( new ConfiguratorType ( Grid ), true ), _displacement( Grid, Displacement ) {}

  explicit FENonlinDeformIntegrationInterface ( const ConfiguratorType &Config, const aol::MultiVector<RealType> &Displacement ) :
    _config ( &Config, false ), _displacement( _config->getInitializer(), Displacement ) {}

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {

    const aol::DiscreteFunctionDefault<ConfiguratorType> arg( _config->getInitializer(), Arg );

    typename ConfiguratorType::ElementType transformedEl;
    typename ConfiguratorType::DomVecType  transformedCoord;

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const typename IteratorType::EndType end_it = _config->end();

    // traverse the elements of the grid and on each element the quadrature points
    for ( IteratorType it = _config->begin(); it != end_it; ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config->getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      RealType elIntegral = aol::ZTrait<RealType>::zero;
      for ( int q = 0; q < numQuadPoints; q++ ) {
        // compute \Phi(x) ("transformedEl", "transformedCoord")
        typename ConfiguratorType::DomVecType localCoord( bfs.getRefCoord ( q ) );
        qc::transformAndClipCoord<ConfiguratorType>( *_config, _displacement, *it, q, localCoord, transformedEl, transformedCoord );

        // compute the integrand at the quadrature points
        elIntegral += this->asImp().evaluateIntegrand ( arg, *it, q, localCoord, transformedEl, transformedCoord ) * bfs.getWeight( q );
      }

      Dest += elIntegral * _config->vol( *it );
    }
  }

  //! Represents \f$f(\phi,\Phi(x),x)\f$ and has to be implemented in the derived class.
  RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                               const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &LocalCoord,
                               const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord ) const {
    return this->asImp().evaluateIntegrand ( DiscFunc, El, QuadPoint, LocalCoord, TransformedEl, TransformedLocalCoord );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

//! Interface to compute \f$\int_\Omega f(\phi,\Phi(x),x)dx\f$,
/** where \f$\Phi\f$ is a deformation on \f$\Omega\f$ and \f$\phi\f$ is the argument of the operator.
 *  If \f$\Phi(x)\f$ comes to lie outside \f$\Omega\f$,
 *  it is projected onto the nearest point in \f$\partial\Omega\f$.
 *
 *  \author Wirth
 */
template <typename ConfiguratorType, typename Imp, int NumComponents = ConfiguratorType::Dim >
class FENonlinDeformIntegrationVectorInterface
      : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  aol::DeleteFlagPointer<const ConfiguratorType> _config;
  // the deformation \Phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement;

public:
  explicit FENonlinDeformIntegrationVectorInterface ( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &Displacement ) :
    _config ( new ConfiguratorType ( Grid ), true ), _displacement( Grid, Displacement ) {}

  explicit FENonlinDeformIntegrationVectorInterface ( const ConfiguratorType &Config, const aol::MultiVector<RealType> &Displacement ) :
    _config ( &Config, false ), _displacement( _config->getInitializer(), Displacement ) {}

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {

    const aol::DiscreteVectorFunctionDefault<ConfiguratorType,NumComponents> arg( _config->getInitializer(), Arg );

    typename ConfiguratorType::ElementType transformedEl;
    typename ConfiguratorType::DomVecType  transformedCoord;

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const typename IteratorType::EndType end_it = _config->end();

    // traverse the elements of the grid and on each element the quadrature points
    for ( IteratorType it = _config->begin(); it != end_it; ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config->getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      RealType elIntegral = aol::ZTrait<RealType>::zero;
      for ( int q = 0; q < numQuadPoints; q++ ) {
        // compute \Phi(x) ("transformedEl", "transformedCoord")
        typename ConfiguratorType::DomVecType localCoord( bfs.getRefCoord ( q ) );
        qc::transformAndClipCoord<ConfiguratorType>( *_config, _displacement, *it, q, localCoord, transformedEl, transformedCoord );

        // compute the integrand at the quadrature points
        elIntegral += this->asImp().evaluateIntegrand ( arg, *it, q, localCoord, transformedEl, transformedCoord ) * bfs.getWeight( q );
      }

      Dest += elIntegral * _config->vol( *it );
    }
  }

  //! Represents \f$f(\phi,\Phi(x),x)\f$ and has to be implemented in the derived class.
  RealType evaluateIntegrand ( const aol::DiscreteVectorFunctionDefault<ConfiguratorType,NumComponents> &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &LocalCoord,
                               const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord ) const {
    return this->asImp().evaluateIntegrand ( DiscFuncs, El, QuadPoint, LocalCoord, TransformedEl, TransformedLocalCoord );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

//! Interface to compute \f$(\int_\Omega f(\phi,\Phi(x),x):\nabla(\varphi_i\circ\Phi(x)) dx)_i\f$
//! or \f$(\int_\Omega f(\phi,\Phi(x),x):(\nabla\varphi_i)(\Phi(x)) dx)_i\f$ (for second, set "DeformDiffMode"=false),
/** where \f$\Phi\f$ is a deformation on \f$\Omega\f$, \f$\varphi_i\f$ is the \f$i\f$th vector-valued FE basis function,
 *  and \f$\phi\f$ is the argument of the operator. If \f$\Phi(x)\f$ comes to lie outside \f$\Omega\f$,
 *  it is projected onto the nearest point in \f$\partial\Omega\f$.
 *
 *  \author Wirth
 */
template <typename ConfiguratorType, int NumCompArg, int NumCompDest, typename Imp, bool DeformDiffMode = true>
class FENonlinDeformVectorDiffOpInterface : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  aol::DeleteFlagPointer<const ConfiguratorType> _config;
  // the deformation \Phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement;

public:
  explicit FENonlinDeformVectorDiffOpInterface ( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &Displacement ) :
    _config ( new ConfiguratorType ( Grid ), true ), _displacement( Grid, Displacement ) {}

  explicit FENonlinDeformVectorDiffOpInterface ( const ConfiguratorType &Config, const aol::MultiVector<RealType> &Displacement ) :
    _config ( &Config, false ), _displacement( _config->getInitializer(), Displacement ) {}

  void applyAdd ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {

    const aol::DiscreteVectorFunctionDefault<ConfiguratorType,NumCompArg> arg( _config->getInitializer(), Arg );

    typename ConfiguratorType::ElementType transformedEl;
    typename ConfiguratorType::DomVecType  transformedCoord;
    typename ConfiguratorType::MatType     dPhi;
    typename ConfiguratorType::VecType     grad;
    aol::Vec<NumCompDest, typename ConfiguratorType::RealType> integrand;

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const typename IteratorType::EndType end_it = _config->end();

    // traverse the elements of the grid and on each element the quadrature points
    for ( IteratorType it = _config->begin(); it != end_it; ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config->getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; q++ ) {
        // compute \Phi(x) ("transformedEl", "transformedCoord") and \nabla\Phi (if needed)
        typename ConfiguratorType::DomVecType localCoord( bfs.getRefCoord ( q ) );
        qc::transformAndClipCoord<ConfiguratorType>( *_config, _displacement, *it, q, localCoord, transformedEl, transformedCoord );
        if ( DeformDiffMode ) {
          _displacement.evaluateGradientAtQuadPoint ( *it, q, dPhi );
          for ( int i = 0; i < ConfiguratorType::Dim; i++ )
            dPhi[i][i] += 1.;
        }

        // compute f(\phi,\Phi(x),x) and its product with \nabla\Phi^T (if needed)
        aol::Mat<NumCompDest, ConfiguratorType::Dim, RealType> f;
        this->asImp().getNonlinearity ( arg, *it, q, localCoord, transformedEl, transformedCoord, f );
        if ( DeformDiffMode ) {
          dPhi.transpose();
          f *= dPhi;
        }

        // add f\nabla\Phi^T\nabla\phi_i (or f\nabla\phi_i) to the ith position in the result
        const int numLocalDofs = _config->getNumLocalDofs ( transformedEl );
        const typename ConfiguratorType::BaseFuncSetType &tbfs = _config->getBaseFunctionSet ( transformedEl );
        for ( int dof = 0; dof < numLocalDofs; ++dof ) {
          tbfs.evaluateGradient ( dof, transformedCoord, grad );
          f.mult( grad, integrand );
          integrand *= bfs.getWeight( q ) * _config->vol( *it );
          for ( int j = 0; j < NumCompDest; ++j )
            Dest[j][_config->localToGlobal( transformedEl, dof )] += integrand[j];
        }
      }
    }
  }

  //! Represents \f$f(\phi,\Phi(x),x)\f$ and has to be implemented in the derived class.
  void getNonlinearity ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFunc,
                         const typename ConfiguratorType::ElementType &El, int QuadPoint, const typename ConfiguratorType::DomVecType &LocalCoord,
                         const typename ConfiguratorType::ElementType &TransformedEl, const typename ConfiguratorType::DomVecType &TransformedLocalCoord,
                         aol::Mat<NumCompDest, ConfiguratorType::Dim, RealType> &NL ) const {
    this->asImp().getNonlinearity ( DiscFunc, El, QuadPoint, LocalCoord, TransformedEl, TransformedLocalCoord, NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

//! type signifying whether in the definition of a linear FE operator
//! all basis functions are deformed or none or only the ones belonging to vectors that are multiplied from the right or left
enum FELinOpBaseFunDeformModeType { BOTH = 0, RIGHT = 1, LEFT = 2, NONE = 3 };

inline bool leftBasisFunDeformed( const FELinOpBaseFunDeformModeType defMode ) { return !( defMode % 2 ); }

inline bool rightBasisFunDeformed( const FELinOpBaseFunDeformModeType defMode ) { return ( defMode < LEFT ); }

//! Interface for linear FE operators with basis functions deformed by a deformation \f$\Phi\f$.
/**
 *  \note If DeformationMode == qc::LEFT (qc::RIGHT) the matrix type ConfiguratorType::MatrixType is possibly not appropiate for assembleAddMatrix<> !
 * 
 *  \author Wirth
 */
template <typename ConfiguratorType, typename Imp, typename _MatrixType = typename ConfiguratorType::MatrixType>
class FELinDeformOpInterface : public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  aol::DeleteFlagPointer<const ConfiguratorType> _config;
  mutable _MatrixType *_mat;
  aol::OperatorType _opType;
  FELinOpBaseFunDeformModeType _deformationMode;
  // the displacement corresponding to the deformation \Phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement;

public:
  explicit FELinDeformOpInterface ( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &Displacement, aol::OperatorType OpType = aol::ONTHEFLY, FELinOpBaseFunDeformModeType DeformationMode = BOTH ) :
    _config ( new ConfiguratorType(Grid), true ),
    _mat ( NULL ),
    _opType ( OpType ),
    _deformationMode( DeformationMode ),
    _displacement( Grid, Displacement ) {}

  virtual ~FELinDeformOpInterface( ) {
    delete _mat;
  }

  //! Clears the assembled matrix.
  void reset( ) {
    if ( _mat )
      delete _mat;
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    switch ( _opType ) {
    case aol::ONTHEFLY:
      multiplyOnTheFly ( Arg, Dest );
      break;
    case aol::ASSEMBLED:
      if ( !_mat )
        assembleMatrix();
      _mat->applyAdd ( Arg, Dest );
      break;
    default:
      throw aol::UnimplementedCodeException ( "FELinDeformOpInterface::applyAdd: unknown opType", __FILE__, __LINE__ );
    };
  }

  _MatrixType& getMatrix( ) {
    if ( !_mat ) {
      _mat = _config->createNewMatrix( );
      assembleMatrix( );
    }
    //return dynamic_cast<typename ConfiguratorType::MatrixType&>(*_mat);
    return *_mat;
  }

protected:
  void multiplyOnTheFly ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {

    typename ConfiguratorType::ElementType transformedEl;
    typename ConfiguratorType::DomVecType  transformedCoord;
    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;
    int globalDofsR[ ConfiguratorType::maxNumLocalDofs ],
        globalDofsL[ ConfiguratorType::maxNumLocalDofs ];

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const typename IteratorType::EndType end_it = _config->end();

    // traverse the elements of the grid and on each element the quadrature points
    for ( IteratorType it = _config->begin(); it != end_it; ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config->getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; q++ ) {
        // compute \Phi(x) ("transformedEl", "transformedCoord")
        typename ConfiguratorType::DomVecType localCoord( bfs.getRefCoord ( q ) );
        qc::transformAndClipCoord<ConfiguratorType>( *_config, _displacement, *it, q, localCoord, transformedEl, transformedCoord );

        // assemble the local matrix belonging to the current quadrature point
        this->asImp().prepareLocalMatrix ( *it, q, localCoord, transformedEl, transformedCoord, localMatrix );

        // get the global indices to the current DoFs belonging to vectors multiplied from the right
        const typename ConfiguratorType::ElementType& elR( rightBasisFunDeformed( _deformationMode ) ? transformedEl : *it );
        const int numLocalDofsR = _config->getNumLocalDofs ( elR );
        for ( int i = 0; i < numLocalDofsR; ++i )
          globalDofsR[ i ] = _config->localToGlobal ( elR, i );

        // get the global indices to the current DoFs belonging to vectors multiplied from the left
        const typename ConfiguratorType::ElementType& elL( leftBasisFunDeformed( _deformationMode ) ? transformedEl : *it );
        const int numLocalDofsL = _config->getNumLocalDofs ( elL );
        for ( int i = 0; i < numLocalDofsL; ++i )
          globalDofsL[ i ] = _config->localToGlobal ( elL, i );

        // add the locally computed value to the global result
        for ( int i = 0; i < numLocalDofsL; ++i )
          for ( int j = 0; j < numLocalDofsR; ++j )
            Dest[ globalDofsL[ i ] ] += localMatrix[ i ][ j ] * Arg[ globalDofsR[ j ] ] ;
      }
    }
  }

  void assembleMatrix( ) const {
    if ( _mat )
      delete _mat;
    _mat = _config->createNewMatrix( );
    assembleAddMatrix ( *_mat );
  }

public:
  /** (this assembled matrix * Factor) is added to Mat  */
  template <typename MatrixType>
  void assembleAddMatrix ( MatrixType &Mat, const RealType Factor = aol::NumberTrait<RealType>::one ) const {

    typename ConfiguratorType::ElementType transformedEl;
    typename ConfiguratorType::DomVecType  transformedCoord;
    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;
    int globalDofsR[ ConfiguratorType::maxNumLocalDofs ],
        globalDofsL[ ConfiguratorType::maxNumLocalDofs ];

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const typename IteratorType::EndType end_it = _config->end();

    // traverse the elements of the grid and on each element the quadrature points
    for ( IteratorType it = _config->begin(); it != end_it; ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config->getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; q++ ) {
        // compute \Phi(x) ("transformedEl", "transformedCoord")
        typename ConfiguratorType::DomVecType localCoord( bfs.getRefCoord ( q ) );
        qc::transformAndClipCoord<ConfiguratorType>( *_config, _displacement, *it, q, localCoord, transformedEl, transformedCoord );

        // assemble the local matrix belonging to the current quadrature point
        this->asImp().prepareLocalMatrix ( *it, q, localCoord, transformedEl, transformedCoord, localMatrix );

        // get the global indices to the current DoFs belonging to vectors multiplied from the right
        const typename ConfiguratorType::ElementType& elR( rightBasisFunDeformed( _deformationMode ) ? transformedEl : *it );
        const int numLocalDofsR = _config->getNumLocalDofs ( elR );
        for ( int i = 0; i < numLocalDofsR; ++i )
          globalDofsR[ i ] = _config->localToGlobal ( elR, i );

        // get the global indices to the current DoFs belonging to vectors multiplied from the left
        const typename ConfiguratorType::ElementType& elL( leftBasisFunDeformed( _deformationMode ) ? transformedEl : *it );
        const int numLocalDofsL = _config->getNumLocalDofs ( elL );
        for ( int i = 0; i < numLocalDofsL; ++i )
          globalDofsL[ i ] = _config->localToGlobal ( elL, i );

        // add the locally assembled matrix to the globally assembled matrix
        for ( int i = 0; i < numLocalDofsL; ++i )
          for ( int j = 0; j < numLocalDofsR; ++j )
            Mat.add ( globalDofsL[ i ], globalDofsR[ j ], Factor * localMatrix [ i ][ j ] );
      }
    }
  }

  template <typename MatrixType>
  void assembleAddMatrix ( const aol::BitVector &Mask, MatrixType &Mat, const RealType Factor = aol::NumberTrait<RealType>::one, const bool SetMaskedRowsToIdentity = true ) const {

    typename ConfiguratorType::ElementType transformedEl;
    typename ConfiguratorType::DomVecType  transformedCoord;
    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;
    int globalDofsR[ ConfiguratorType::maxNumLocalDofs ],
        globalDofsL[ ConfiguratorType::maxNumLocalDofs ];

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const typename IteratorType::EndType end_it = _config->end();

    // traverse the elements of the grid and on each element the quadrature points
    for ( IteratorType it = _config->begin(); it != end_it; ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config->getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; q++ ) {
        // compute \Phi(x) ("transformedEl", "transformedCoord")
        typename ConfiguratorType::DomVecType localCoord( bfs.getRefCoord ( q ) );
        qc::transformAndClipCoord<ConfiguratorType>( *_config, _displacement, *it, q, localCoord, transformedEl, transformedCoord );

        // assemble the local matrix belonging to the current quadrature point
        this->asImp().prepareLocalMatrix ( *it, q, localCoord, transformedEl, transformedCoord, localMatrix );

        // get the global indices to the current DoFs belonging to vectors multiplied from the right
        const typename ConfiguratorType::ElementType& elR( rightBasisFunDeformed( _deformationMode ) ? transformedEl : *it );
        const int numLocalDofsR = _config->getNumLocalDofs ( elR );
        for ( int i = 0; i < numLocalDofsR; ++i )
          globalDofsR[ i ] = _config->localToGlobal ( elR, i );

        // get the global indices to the current DoFs belonging to vectors multiplied from the left
        const typename ConfiguratorType::ElementType& elL( leftBasisFunDeformed( _deformationMode ) ? transformedEl : *it );
        const int numLocalDofsL = _config->getNumLocalDofs ( elL );
        for ( int i = 0; i < numLocalDofsL; ++i )
          globalDofsL[ i ] = _config->localToGlobal ( elL, i );

        // add the locally assembled matrix to the globally assembled matrix,
        // but write the ith row only if the ith node is not a Dirichlet node and the jth column only if the jth node is not a Dirichlet node
        for ( int i = 0; i < numLocalDofsL; ++i )
          if ( !Mask[ globalDofsL[ i ] ] )
            for ( int j = 0; j < numLocalDofsR; ++j )
              if ( !Mask[ globalDofsR[ j ] ] )
                Mat.add ( globalDofsL[ i ], globalDofsR[ j ], Factor * localMatrix [ i ][ j ] );
      }
    }

    // set ones on the diagonal for Dirichlet nodes
    if ( SetMaskedRowsToIdentity )
      for ( int i = 0; i < Mask.size(); i++ )
        if ( Mask[i] )
          Mat.add( i, i, Factor );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

//! Interface to compute \f$(\int_\Omega w(x,\Phi(x))(\varphi_i\circ\Phi)(\varphi_j\circ\Phi) dx)_{ij}\f$,
/** where \f$\Phi\f$ is a deformation on \f$\Omega\f$, \f$\varphi_i\f$ is the \f$i\f$th FE basis function,
 *  and the scalar weight \f$w\f$ has to be provided in derived classes. If \f$\Phi(x)\f$ comes to lie outside \f$\Omega\f$,
 *  it is projected onto the nearest point in \f$\partial\Omega\f$.
 *
 *  \author Wirth
 */
template <typename ConfiguratorType, typename Imp, typename MatrixType = typename ConfiguratorType::MatrixType>
class FELinDeformScalarWeightedMassInterface :
  public FELinDeformOpInterface< ConfiguratorType, FELinDeformScalarWeightedMassInterface<ConfiguratorType, Imp, MatrixType>, MatrixType > {

protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

public:
  explicit FELinDeformScalarWeightedMassInterface ( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &Displacement, aol::OperatorType OpType = aol::ONTHEFLY, FELinOpBaseFunDeformModeType DeformationMode = BOTH ) :
    FELinDeformOpInterface<ConfiguratorType, FELinDeformScalarWeightedMassInterface<ConfiguratorType, Imp, MatrixType>, MatrixType > ( Grid, Displacement, OpType, DeformationMode ) {}

  //! This function has to be provided in the implementation (derived class) of the interface.
  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const DomVecType &LocalCoord,
                             const typename ConfiguratorType::ElementType &TransformedEl, const DomVecType &TransformedLocalCoord ) const {
    return this->asImp().getCoeff ( El, QuadPoint, LocalCoord, TransformedEl, TransformedLocalCoord );
  }

  //! Performs the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const DomVecType &LocalCoord,
                                   const typename ConfiguratorType::ElementType &TransformedEl, const DomVecType &TransformedLocalCoord,
                                   aol::Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> &LocalMatrix ) const {

    const RealType coeff = getCoeff ( El, QuadPoint, LocalCoord, TransformedEl, TransformedLocalCoord );
    aol::Vec<ConfiguratorType::maxNumLocalDofs,RealType> basisFunctionValuesR, basisFunctionValuesL;

    // compute the basis function values belonging to vectors multiplied from the right
    const typename ConfiguratorType::ElementType& elR( rightBasisFunDeformed( this->_deformationMode ) ? TransformedEl : El );
    const int numDofsR = this->_config->getNumLocalDofs ( elR );
    const typename ConfiguratorType::BaseFuncSetType &tbfsR = this->_config->getBaseFunctionSet ( elR );
    for ( int i = 0; i < numDofsR; i++ )
      basisFunctionValuesR[i] = tbfsR.evaluate( i, rightBasisFunDeformed( this->_deformationMode ) ? TransformedLocalCoord : LocalCoord );

    // compute the basis function values belonging to vectors multiplied from the left
    const typename ConfiguratorType::ElementType& elL( leftBasisFunDeformed( this->_deformationMode ) ? TransformedEl : El );
    const int numDofsL = this->_config->getNumLocalDofs ( elL );
    const typename ConfiguratorType::BaseFuncSetType &tbfsL = this->_config->getBaseFunctionSet ( elL );
    for ( int i = 0; i < numDofsL; i++ )
      basisFunctionValuesL[i] = tbfsL.evaluate( i, leftBasisFunDeformed( this->_deformationMode ) ? TransformedLocalCoord : LocalCoord );

    // assemble the local mass matrix
    const typename ConfiguratorType::BaseFuncSetType &bfs = this->_config->getBaseFunctionSet ( El );
    LocalMatrix.makeTensorProduct( basisFunctionValuesL, basisFunctionValuesR );
    LocalMatrix *= coeff * bfs.getWeight ( QuadPoint ) * this->_config->vol ( El );
  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }
};

//! Interface to compute \f$(\int_\Omega A(x,\Phi(x))\nabla(\varphi_j\circ\Phi)\cdot\nabla(\varphi_i\circ\Phi) dx)_{ij}\f$,
//! or \f$(\int_\Omega A(x,\Phi(x))(\nabla\varphi_j)\circ\Phi\cdot(\nabla\varphi_i)\circ\Phi dx)_{ij}\f$ (for second, set "DeformDiffMode"=false),
/** where \f$\Phi\f$ is a deformation on \f$\Omega\f$, \f$\varphi_i\f$ is the \f$i\f$th FE basis function,
 *  and the (asymmetric) matrix \f$A\f$ has to be provided in derived classes. If \f$\Phi(x)\f$ comes to lie outside \f$\Omega\f$,
 *  it is projected onto the nearest point in \f$\partial\Omega\f$.
 *
 *  \author Wirth
 */
template <typename ConfiguratorType, typename Imp, bool DeformDiffMode = true, typename MatrixType = typename ConfiguratorType::MatrixType>
class FELinDeformAsymMatrixWeightedStiffInterface :
  public FELinDeformOpInterface< ConfiguratorType,FELinDeformAsymMatrixWeightedStiffInterface<ConfiguratorType, Imp, DeformDiffMode, MatrixType>, MatrixType > {

protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::MatType    MatType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

public:
  explicit FELinDeformAsymMatrixWeightedStiffInterface ( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &Displacement, aol::OperatorType OpType = aol::ONTHEFLY, FELinOpBaseFunDeformModeType DeformationMode = BOTH ) :
    FELinDeformOpInterface<ConfiguratorType, FELinDeformAsymMatrixWeightedStiffInterface<ConfiguratorType, Imp, DeformDiffMode> > ( Grid, Displacement, OpType, DeformationMode ) {}

  //! This function has to be provided in the implementation (derived class) of the interface. It returns the matrix \f$A\f$.
  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const DomVecType &LocalCoord,
                               const typename ConfiguratorType::ElementType &TransformedEl, const DomVecType &TransformedLocalCoord,
                               MatType &Matrix ) const {
    this->asImp().getCoeffMatrix ( El, QuadPoint, LocalCoord, TransformedEl, TransformedLocalCoord, Matrix );
  }

  //! Performs the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const DomVecType &LocalCoord,
                                   const typename ConfiguratorType::ElementType &TransformedEl, const DomVecType &TransformedLocalCoord,
                                   aol::Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> &LocalMatrix ) const {

    MatType mat, dPhi;
    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::Dim,RealType> gradL, gradR, auxMat2;

    // compute the matrices A and \nabla\Phi and the product \nabla\Phi A \nabla\Phi^T (or corresponding product for _deformationMode!=BOTH), if needed
    getCoeffMatrix ( El, QuadPoint, LocalCoord, TransformedEl, TransformedLocalCoord, mat ); // A
    if ( DeformDiffMode ) {
      this->_displacement.evaluateGradientAtQuadPoint ( El, QuadPoint, dPhi );               // \nabla\Phi
      for ( int i = 0; i < ConfiguratorType::Dim; i++ )
        dPhi[i][i] += 1.;
      MatType auxMat( mat );
      if ( leftBasisFunDeformed( this->_deformationMode ) )
        auxMat.makeProduct( dPhi, mat );
      mat = auxMat;
      if ( rightBasisFunDeformed( this->_deformationMode ) )
        mat *= dPhi.transposed();                                                            // the product
    }

    // compute the gradients of the basis functions belonging to vectors multiplied from the right
    const typename ConfiguratorType::ElementType& elR( rightBasisFunDeformed( this->_deformationMode ) ? TransformedEl : El );
    const int numDofsR = this->_config->getNumLocalDofs ( elR );
    const typename ConfiguratorType::BaseFuncSetType &tbfsR = this->_config->getBaseFunctionSet ( elR );
    for ( int i = 0; i < numDofsR; i++ )
      tbfsR.evaluateGradient( i, rightBasisFunDeformed( this->_deformationMode ) ? TransformedLocalCoord : LocalCoord, static_cast<VecType&>( gradR[i] ) );

    // compute the gradients of the basis functions belonging to vectors multiplied from the left
    const typename ConfiguratorType::ElementType& elL( leftBasisFunDeformed( this->_deformationMode ) ? TransformedEl : El );
    const int numDofsL = this->_config->getNumLocalDofs ( elL );
    const typename ConfiguratorType::BaseFuncSetType &tbfsL = this->_config->getBaseFunctionSet ( elL );
    for ( int i = 0; i < numDofsL; i++ )
      tbfsL.evaluateGradient( i, leftBasisFunDeformed( this->_deformationMode ) ? TransformedLocalCoord : LocalCoord, static_cast<VecType&>( gradL[i] ) );

    // assemble the local stiffness matrix
    const typename ConfiguratorType::BaseFuncSetType &bfs = this->_config->getBaseFunctionSet ( El );
    auxMat2.makeProduct( gradL, mat );
    LocalMatrix.makeProductABtransposed( auxMat2, gradR );
    LocalMatrix *= bfs.getWeight ( QuadPoint ) * this->_config->vol ( El );
  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }
};

//! Class to compute \f$(\int_\Omega (\varphi_i\circ\Phi)(\varphi_j\circ\Phi) dx)_{ij}\f$,
/** where \f$\Phi\f$ is a deformation on \f$\Omega\f$, \f$\varphi_i\f$ is the \f$i\f$th FE basis function.
 *  If \f$\Phi(x)\f$ comes to lie outside \f$\Omega\f$, it is projected onto the nearest point in \f$\partial\Omega\f$.
 *
 *  \author Wirth
 */
template <typename ConfiguratorType, typename MatrixType = typename ConfiguratorType::MatrixType>
class DeformMassOp :
  public FELinDeformScalarWeightedMassInterface<ConfiguratorType, DeformMassOp<ConfiguratorType, MatrixType>, MatrixType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

public:
  explicit DeformMassOp ( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &Displacement, aol::OperatorType OpType = aol::ONTHEFLY, FELinOpBaseFunDeformModeType DeformationMode = BOTH ) :
    FELinDeformScalarWeightedMassInterface<ConfiguratorType, DeformMassOp<ConfiguratorType, MatrixType>, MatrixType> ( Grid, Displacement, OpType, DeformationMode ) {}

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/, const DomVecType &/*LocalCoord*/,
                             const typename ConfiguratorType::ElementType &/*TransformedEl*/, const DomVecType &/*TransformedLocalCoord*/ ) const {
    return 1;
  }
};

//! Class to compute \f$(\int_\Omega w(x)^2(\varphi_i\circ\Phi)(\varphi_j\circ\Phi) dx)_{ij}\f$,
/** where \f$\Phi\f$ is a deformation on \f$\Omega\f$, \f$\varphi_i\f$ is the \f$i\f$th FE basis function,
 *  and the scalar weight \f$w\f$ is given as a discretized function. If \f$\Phi(x)\f$ comes to lie outside \f$\Omega\f$,
 *  it is projected onto the nearest point in \f$\partial\Omega\f$.
 *
 *  \author Wirth
 */
template <typename ConfiguratorType, typename MatrixType = typename ConfiguratorType::MatrixType>
class SquaredWeightDeformMassOp :
  public FELinDeformScalarWeightedMassInterface<ConfiguratorType, SquaredWeightDeformMassOp<ConfiguratorType,MatrixType>, MatrixType> {

protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

  // the weight to be squared
  const aol::DiscreteFunctionDefault<ConfiguratorType> _weight;

public:
  explicit SquaredWeightDeformMassOp ( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &Displacement, const aol::Vector<RealType> &Weight, aol::OperatorType OpType = aol::ONTHEFLY, FELinOpBaseFunDeformModeType DeformationMode = BOTH ) :
    FELinDeformScalarWeightedMassInterface<ConfiguratorType, SquaredWeightDeformMassOp<ConfiguratorType> > ( Grid, Displacement, OpType, DeformationMode ),
    _weight( Grid, Weight ) {}

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const DomVecType &/*LocalCoord*/,
                             const typename ConfiguratorType::ElementType &/*TransformedEl*/, const DomVecType &/*TransformedLocalCoord*/ ) const {
    return aol::Sqr( _weight.evaluateAtQuadPoint( El, QuadPoint ) );
  }
};


//! Class to compute \f$ E[u] = \int_\Omega |\det( 1 + Du)| dx = |\phi(\Omega)| \f$,
/** where \f$ u : \Omega \subset R^n \rightarrow R^n \f$,  \f$ u(x) = \phi(x) - x \f$ is a displacement.
 *
 *  \author Wirth
 */
template <typename ConfiguratorType>
class DeformedVolumeIntegrator :
  public aol::FENonlinIntegrationVectorInterface<ConfiguratorType,DeformedVolumeIntegrator<ConfiguratorType>,ConfiguratorType::Dim> {

private:
  typedef typename ConfiguratorType::RealType RealType;

public:
  DeformedVolumeIntegrator( const typename ConfiguratorType::InitType &Grid ) :
    aol::FENonlinIntegrationVectorInterface<ConfiguratorType,DeformedVolumeIntegrator<ConfiguratorType>,ConfiguratorType::Dim>( Grid ) {}

  RealType evaluateIntegrand( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                              const typename ConfiguratorType::ElementType &El, int QuadPoint,
                              const typename ConfiguratorType::DomVecType &/*RefCoord*/ ) const {
    // compute the deformation gradient $\nabla\phi$: for each displacment component do...
    typename ConfiguratorType::MatType dphi;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      // evaluate gradient of the displacement $d_i$
      DiscFuncs[i].evaluateGradientAtQuadPoint ( El, QuadPoint, dphi[i] );
      // compute corresponding deformation gradient of the $i$th component
      dphi[i][i] += 1.;
    }

    return aol::Abs( dphi.det() );
  }
};

/** This class adds a (linearized or nonlinear) rigid body motion onto a given displacement such that mean displacement and moment are zero.
 *
 *  \author Wirth
 */
template <typename ConfiguratorType>
class meanZeroShiftAndRotationProjector :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const bool _linearized;
  aol::Vector<RealType> _integrator;
  aol::MultiVector<RealType> _momentIntegrator;
  aol::MultiVector<RealType> _identity;
  RealType _volume;
  typename ConfiguratorType::VecType _center;
  typename ConfiguratorType::MatType _moments;

  class Displacement3DMomentNormSqr :
    public aol::Op< aol::Vec3<RealType>, aol::Scalar<RealType> > {
  private:
    aol::Matrix33<RealType> _outerProd;
  public:
    explicit Displacement3DMomentNormSqr( const aol::Matrix33<RealType> OuterProd ) :
      _outerProd( OuterProd ) {}
    void applyAdd( const aol::Vec3<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
      // compute the rotation matrix from the given three angles
      aol::Matrix33<RealType> rotMat, auxMat;
      rotMat.setRotationAboutAxis( Arg[0], 0 );
      for ( int i = 1; i < 3; i++ ) {
        auxMat.setRotationAboutAxis( Arg[i], i );
        rotMat *= auxMat;
      }
      // compute the corresponding moment
      int index1[3] = { 1, 2, 0 };
      int index2[3] = { 2, 0, 1 };
      aol::Vec3<RealType> moment;
      for ( int i = 0; i < qc::QC_3D; i++ )
        for ( int j = 0; j < qc::QC_3D; j++ )
          moment[i] += _outerProd[index1[i]][j] * rotMat[index2[i]][j] - _outerProd[index2[i]][j] * rotMat[index1[i]][j];
      // return the square of the moment
      Dest += moment.normSqr();
    }
  };

  class Displacement3DMomentNormSqrGradient :
    public aol::Op< aol::Vec3<RealType>, aol::Vec3<RealType> > {
  private:
    aol::Matrix33<RealType> _outerProd;
  public:
    explicit Displacement3DMomentNormSqrGradient( const aol::Matrix33<RealType> OuterProd ) :
      _outerProd( OuterProd ) {}
    void applyAdd( const aol::Vec3<RealType> &Arg, aol::Vec3<RealType> &Dest ) const {
      // compute the rotation matrix from the given three angles
      aol::Matrix33<RealType> rotMat;
      rotMat.setIdentity();
      std::vector< aol::Matrix33<RealType> > axRotMat( 3 ), axRotMatDeriv( 3 );
      for ( int i = 0; i < 3; i++ ) {
        axRotMat[i].setRotationAboutAxis( Arg[i], i );
        rotMat *= axRotMat[i];
        axRotMatDeriv[i].setRotationAboutAxis( Arg[i] + ( aol::NumberTrait<long double>::pi / 2. ), i );
        for ( int j = 0; j < qc::QC_3D; j++ )
          axRotMatDeriv[i][j][j] = 0.;
      }
      // compute the corresponding moment
      int index1[3] = { 1, 2, 0 };
      int index2[3] = { 2, 0, 1 };
      aol::Vec3<RealType> moment;
      for ( int i = 0; i < qc::QC_3D; i++ )
        for ( int j = 0; j < qc::QC_3D; j++ )
          moment[i] += _outerProd[index1[i]][j] * rotMat[index2[i]][j] - _outerProd[index2[i]][j] * rotMat[index1[i]][j];
      // compute and return the moment derivatives
      for( int k = 0; k < qc::QC_3D; k++ ) {
        // compute derivative of rotation matrix wrt kth angle
        aol::Matrix33<RealType> rotMat;
        rotMat.setIdentity();
        for ( int i = 0; i < qc::QC_3D; i++ )
          if ( i != k )
            rotMat *= axRotMat[i];
          else
            rotMat *= axRotMatDeriv[i];
        aol::Vec3<RealType> momentDeriv;
        for ( int i = 0; i < qc::QC_3D; i++ )
          for ( int j = 0; j < qc::QC_3D; j++ )
            momentDeriv[i] += _outerProd[index1[i]][j] * rotMat[index2[i]][j] - _outerProd[index2[i]][j] * rotMat[index1[i]][j];
        Dest[k] += 2 * ( moment * momentDeriv );
      }
    }
  };

public:
  explicit meanZeroShiftAndRotationProjector( const typename ConfiguratorType::InitType &Grid, const aol::Vector<RealType> &Weight, const bool Linearized = true ) :
    _grid( Grid ),
    _linearized( Linearized ),
    _integrator( Grid.getNumberOfNodes() ),
    _momentIntegrator( ConfiguratorType::Dim, Grid.getNumberOfNodes() ),
    _identity( ConfiguratorType::Dim, Grid.getNumberOfNodes() ) {

    // compute the integration operator, \int_\Omega w ... dx
    aol::MassOp<ConfiguratorType>( Grid ).apply( Weight, _integrator );

    // compute the moment integration operator, \int_O w x ... dx
    qc::DataGenerator<ConfiguratorType>( Grid ).generateIdentity( _identity );
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      aol::WeightedMassOp<ConfiguratorType>( Grid, Weight ).apply( _identity[i], _momentIntegrator[i] );

    // compute volume, \int_\Omega w dx
    aol::Vector<RealType> ones( Grid.getNumberOfNodes() );
    ones.setAll( 1. );
    _volume = _integrator * ones;

    // compute the center of mass, x_c = \int_\Omega w x dx / \int_\Omega w dx
    for ( int j = 0; j < ConfiguratorType::Dim; j++ )
      _center[j] = _momentIntegrator[j] * ones;
    _center /= _volume;

    // compute second moments, \int_\Omega w (x-x_c)_i (x-x_c)_j dx
    for ( int k = 0; k < ConfiguratorType::Dim; k++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        _moments[k][j] = _momentIntegrator[k] * _identity[j] - _volume * _center[k] * _center[j];
  }

  //! computes total shift, \f$ \int_\Omega w u dx \f$
  typename ConfiguratorType::VecType totalShift( const aol::MultiVector<RealType> &Arg ) const {
    typename ConfiguratorType::VecType shift;
    for ( int j = 0; j < ConfiguratorType::Dim; j++ )
      shift[j] = _integrator * Arg[j];
    return shift;
  }

  //! computes total moment, \f$ \int_\Omega w x \times u dx \f$
  aol::Vec3<RealType> totalMoment( const aol::MultiVector<RealType> &Arg ) const {
    aol::Vec3<RealType> moment;
    moment[2] = _momentIntegrator[0] * Arg[1] - _momentIntegrator[1] * Arg[0];
    if ( ConfiguratorType::Dim == qc::QC_3D ) {
      moment[0] = _momentIntegrator[1] * Arg[2] - _momentIntegrator[2] * Arg[1];
      moment[1] = _momentIntegrator[2] * Arg[0] - _momentIntegrator[0] * Arg[2];
    }
    return moment;
  }

  //! computes a skew-symmetric matrix \f$ R \f$ such that \f$ \int_\Omega w x \times [R(x-x_c)] dx = m \f$,
  //! where \f$ x_c \f$ is \f$ \int_\Omega w x dx \f$ and \f$ m \f$ the given moment
  aol::Matrix33<RealType> momentMatrixLinearized( const aol::Vec3<RealType> &Moment ) const {
    aol::Matrix33<RealType> matrix, auxMat;
    if ( ConfiguratorType::Dim == qc::QC_2D ) {
      matrix[0][1] = Moment[2] / ( _moments[0][0] + _moments[1][1] );
      matrix[1][0] = -matrix[0][1];
    } else {
      auxMat[0][0] =  _moments[0][2];
      auxMat[0][1] = -_moments[0][1];
      auxMat[0][2] = -_moments[1][1] - _moments[2][2];
      auxMat[1][0] =  _moments[1][2];
      auxMat[1][1] =  _moments[0][0] + _moments[2][2];
      auxMat[1][2] =  _moments[0][1];
      auxMat[2][0] = -_moments[0][0] - _moments[1][1];
      auxMat[2][1] = -_moments[1][2];
      auxMat[2][2] =  _moments[0][2];
      aol::Vec3<RealType> abc;
      abc = auxMat.inverse() * Moment;
      matrix[0][1] =  abc[0];
      matrix[0][2] =  abc[1];
      matrix[1][0] = -abc[0];
      matrix[1][2] =  abc[2];
      matrix[2][0] = -abc[1];
      matrix[2][1] = -abc[2];
    }
    return matrix;
  }

  //! computes a rotation matrix \f$ R \f$ such that \f$ \int_\Omega w x \times (R d) dx = 0 \f$,
  //! where \f$ d = x + u - x_c \f$ is given in ``ShiftedDeform''
  typename ConfiguratorType::MatType momentMatrixNonlinear( const aol::MultiVector<RealType> &ShiftedDeform ) const {
    typename ConfiguratorType::MatType matrix;
    aol::Matrix33<RealType> outerProd;
    // create (\int_\Omega w x_k(u_j(x)+x_j-x_{cj}) dx)_{kj}
    for ( int k = 0; k < ConfiguratorType::Dim; k++ )
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        outerProd[k][j] = _momentIntegrator[k] * ShiftedDeform[j];
    // compute rotation matrix such that \int_\Omega w x \times [R(x+u-x_c)+x_c-x] dx = \int_\Omega w x \times [R(x+u-x_c)] dx = 0
    if ( ConfiguratorType::Dim == qc::QC_2D ) {
      RealType aux1 = outerProd[0][0] + outerProd[1][1],
               aux2 = outerProd[0][1] - outerProd[1][0];
      RealType sin = sqrt( 1. / ( 1. + aol::Sqr( aux1 / aux2 ) ) ) * ( aux1 * aux2 < 0 ? 1 : -1 ),
               cos = sqrt( 1. - aol::Sqr( sin ) ); // note: \int_\Omega (w x \times u)_3 dx = \int_\Omega (w x \times ( u + x - x_c) ) dx
      matrix[0][1] = -sin;
      matrix[1][0] = sin;
      matrix[0][0] = cos;
      matrix[1][1] = cos;
    } else {
      aol::Vec3<RealType> initialAngles, angles;
      Displacement3DMomentNormSqr e( outerProd );
      Displacement3DMomentNormSqrGradient de( outerProd );
      aol::QuasiNewtonIteration< RealType, aol::Vec3<RealType>, aol::Vec3<RealType> > descent( e, de, 1000, 1.e-20, 10, false, NULL );
      descent.getNewtonInfo().setQuietMode( true );
      descent.setTimestepController( aol::NewtonInfo<RealType>::WOLFE );
      descent.apply( initialAngles, angles );
      aol::Matrix33<RealType> auxMat, resultMat;
      resultMat.setRotationAboutAxis( angles[0], 0 );
      for ( int i = 1; i < 3; i++ ) {
        auxMat.setRotationAboutAxis( angles[i], i );
        resultMat *= auxMat;
      }
      matrix.setSubMatrix( 0, 0, resultMat );
    }
    return matrix;
  }

  void apply( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest ) const {
    if ( &Dest != &Arg )
      Dest = Arg;

    // compute mean shift and make it zero
    typename ConfiguratorType::VecType shift( totalShift( Arg ) );
    for ( int j = 0; j < ConfiguratorType::Dim; j++ )
      Dest[j].addToAll( - shift[j] / _volume );

    // compute the moment of the displacement (with zero total shift), \int_\Omega w x \times u dx, and make it zero
    aol::Vec3<RealType> moment( totalMoment( Dest ) );
    aol::MultiVector<RealType> momentDisp( Arg, aol::DEEP_COPY );
    if ( _linearized ) {
      aol::Matrix33<RealType> matrix( momentMatrixLinearized( moment ) );
      typename ConfiguratorType::MatType momentMat;
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        for ( int k = 0; k < ConfiguratorType::Dim; k++ )
          momentMat[j][k] = - matrix[j][k];
      shift = momentMat * _center;
      shift *= -1.;
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        momentMat[j][j] = 1.;
      qc::DataGenerator<ConfiguratorType>( _grid ).generateAffineDisplacement( momentMat, shift, momentDisp );
    } else {
      // use the following: given the center x_c of \Omega and a deformation \phi and displacement u=\phi-id with zero mean,
      // find rotation \psi or rotation matrix R such that \phi_0=R(\phi-x_c)+x_c with displacement u_0=u+(R-I)(u+id-x_c) has zero moment
      // create u+id-x_c
      momentDisp += _identity;
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        momentDisp[j].addToAll( - _center[j] );
      // find R and compute R-I
      typename ConfiguratorType::MatType rotMat( momentMatrixNonlinear( momentDisp ) );
      for ( int j = 0; j < ConfiguratorType::Dim; j++ )
        rotMat[j][j] -= 1.;
      // compute (R-I)(u+id-x_c)
      for ( int k = 0; k < momentDisp[0].size(); k++ ) {
        typename ConfiguratorType::VecType vec, newVec;
        for ( int j = 0; j < ConfiguratorType::Dim; j++ )
          vec[j] = momentDisp[j][k];
        newVec = rotMat * vec;
        for ( int j = 0; j < ConfiguratorType::Dim; j++ )
          momentDisp[j][k] = newVec[j];
      }
    }
    Dest += momentDisp;
  }

  void applyAdd( const aol::MultiVector<RealType> &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const {
    throw aol::Exception( "applyAdd makes no sense in meanZeroShiftAndRotationProjector" );
  }
};

//! Class to compute the energy \f$ E[\phi] := \int_\Omega (\phi(u_t) - u_f)^2 dx \f$,
/** where \f$\phi\f$ is a deformation on \f$\Omega\f$, \f$u_t\f$ and \f$u_f\f$ are images.
 *
 *  \author Wirth
 */
template <typename ConfiguratorType>
class MismatchEnergyWRTDisp :
  public aol::FENonlinIntegrationVectorInterface<ConfiguratorType,MismatchEnergyWRTDisp<ConfiguratorType>,ConfiguratorType::Dim> {

private:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _UFixed, _UTemplate;

public:
  MismatchEnergyWRTDisp( const typename ConfiguratorType::InitType &Grid,
                         const aol::Vector<RealType> &UFixed,
                         const aol::Vector<RealType> &UTemplate ) :
    aol::FENonlinIntegrationVectorInterface<ConfiguratorType,MismatchEnergyWRTDisp<ConfiguratorType>,ConfiguratorType::Dim>( Grid ),
    _UFixed( Grid, UFixed ),
    _UTemplate( Grid, UTemplate ) {}

  RealType evaluateIntegrand( const aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                              const typename ConfiguratorType::ElementType &El, int QuadPoint,
                              const typename ConfiguratorType::DomVecType &RefCoord ) const {
    typename ConfiguratorType::DomVecType transformedLocalCoord;
    typename ConfiguratorType::ElementType transformedEl;
    // attention! Clipping is important to extend the image "continuously" outside the computational domain and thus prevent energy discontinuities wrt displacements.
    qc::transformAndClipCoord<ConfiguratorType> ( *this->_config, DiscFuncs, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord );

    return aol::Sqr( _UTemplate.evaluate( transformedEl, transformedLocalCoord ) - _UFixed.evaluateAtQuadPoint( El, QuadPoint ) );
  }
};

//! Class to compute the variation \f$ \partial_\phi E[\phi]\f$ of \f$ E[\phi] := \int_\Omega (\phi(u_t) - u_f)^2 dx \f$,
/** where \f$\phi\f$ is a deformation on \f$\Omega\f$, \f$u_f\f$ and \f$u_t\f$ are images.
 *
 *  \author Wirth
 */
template <typename ConfiguratorType>
class MismatchEnergyVariationWRTDisp :
  public aol::FENonlinVectorOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim,MismatchEnergyVariationWRTDisp<ConfiguratorType> > {

private:
  typedef typename ConfiguratorType::RealType RealType;

  // the fixed and deformable images to be matched
  const aol::DiscreteFunctionDefault<ConfiguratorType> _UFixed, _UTemplate;

public:
  MismatchEnergyVariationWRTDisp( const typename ConfiguratorType::InitType &Grid,
                                  const aol::Vector<RealType> &UFixed,
                                  const aol::Vector<RealType> &UTemplate ) :
    aol::FENonlinVectorOpInterface<ConfiguratorType,ConfiguratorType::Dim,ConfiguratorType::Dim,MismatchEnergyVariationWRTDisp<ConfiguratorType> >( Grid ),
    _UFixed( Grid, UFixed ),
    _UTemplate( Grid, UTemplate ) {}

  void getNonlinearity( aol::auto_container<ConfiguratorType::Dim, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscFuncs,
                        const typename ConfiguratorType::ElementType &El, int QuadPoint,
                        const typename ConfiguratorType::DomVecType &RefCoord,
                        aol::Vec<ConfiguratorType::Dim,RealType > &NL) const {
    typename ConfiguratorType::DomVecType transformedLocalCoord;
    typename ConfiguratorType::ElementType transformedEl;
    aol::Vec<ConfiguratorType::Dim,bool> coordinateWithinLimits;
    qc::transformAndClipCoord<ConfiguratorType> ( *this->_config, DiscFuncs, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord, coordinateWithinLimits );

    RealType uTemplate = _UTemplate.evaluate( transformedEl, transformedLocalCoord );
    _UTemplate.evaluateGradient( transformedEl, transformedLocalCoord, NL );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      if ( coordinateWithinLimits[i] == false )
        NL[i] = 0.;

    NL *= 2. * ( uTemplate - _UFixed.evaluateAtQuadPoint( El, QuadPoint ) );
  }
};

//! Helping class to compute the second variation \f$ \partial_\phi^2 E[\phi]\f$ of \f$ E[\phi] := \int_\Omega (\phi(u_t) - u_f)^2 dx \f$,
/** 
 *  \author Wirth
 *  \ingroup MatrixFEOp
 */
template <typename ConfiguratorType>
class MismatchEnergySecondVariationWRTDispAssembly :
  public aol::FELinMatrixWeightedVectorMassInterface<ConfiguratorType,MismatchEnergySecondVariationWRTDispAssembly<ConfiguratorType> > {

private:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::MatType    MatType;
  typedef typename ConfiguratorType::VecType    VecType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

  // the fixed and deformable images to be matched
  const aol::DiscreteFunctionDefault<ConfiguratorType> _UFixed, _UTemplate;
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _gradUTemplate, _displacement;

public:
  MismatchEnergySecondVariationWRTDispAssembly( const typename ConfiguratorType::InitType &Grid,
                                                const aol::Vector<RealType> &UFixed,
                                                const aol::Vector<RealType> &UTemplate,
                                                const aol::MultiVector<RealType> &GradUTemplate,
                                                const aol::MultiVector<RealType> &Displacement ) :
    aol::FELinMatrixWeightedVectorMassInterface<ConfiguratorType,MismatchEnergySecondVariationWRTDispAssembly<ConfiguratorType> >( Grid ),
    _UFixed( Grid, UFixed ),
    _UTemplate( Grid, UTemplate ),
    _gradUTemplate( Grid, GradUTemplate ),
    _displacement( Grid, Displacement ) {}

  inline void getCoeffMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const DomVecType &RefCoord, MatType &Matrix ) const {
    typename ConfiguratorType::DomVecType transformedLocalCoord;
    typename ConfiguratorType::ElementType transformedEl;
    aol::Vec<ConfiguratorType::Dim,bool> coordinateWithinLimits;
    qc::transformAndClipCoord<ConfiguratorType> ( *this->_config, _displacement, El, QuadPoint, RefCoord, transformedEl, transformedLocalCoord, coordinateWithinLimits );

    RealType uTemplate = _UTemplate.evaluate( transformedEl, transformedLocalCoord );
    VecType gradUTemplate;
    _UTemplate.evaluateGradient( transformedEl, transformedLocalCoord, gradUTemplate );

    _gradUTemplate.evaluateGradient( transformedEl, transformedLocalCoord, Matrix );
    Matrix *= ( uTemplate - _UFixed.evaluateAtQuadPoint( El, QuadPoint ) );

    MatType auxMat;
    auxMat.makeTensorProduct( gradUTemplate, gradUTemplate );
    Matrix += auxMat;

    // Matrix is not quite symmetric, since \nabla[(\nabla G_\varepsilon)*u] is numerically only approximately symmetric
    // hence, instead of multiplying with 2, add the transpose
    Matrix += Matrix.transposed();

    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      if ( coordinateWithinLimits[i] == false )
        for ( int j = 0; j < ConfiguratorType::Dim; j++ ) {
          Matrix[i][j] = 0.;
          Matrix[j][i] = 0.;
        }
  }
};

//! Class to compute the second variation \f$ \partial_\phi^2 E[\phi]\f$ of \f$ E[\phi] := \int_\Omega (\phi(u_t) - u_f)^2 dx \f$,
/** 
 *  \author Wirth
 */
template <typename ConfiguratorType, typename BlockOpType = aol::BlockMatrix<typename ConfiguratorType::MatrixType> >
class MismatchEnergySecondVariationWRTDisp :
  public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>,BlockOpType> {

private:
  typedef typename ConfiguratorType::RealType   RealType;
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> &_UFixed, &_UTemplate;
  const aol::MultiVector<RealType> &_gradUTemplate;

public:
  MismatchEnergySecondVariationWRTDisp( const typename ConfiguratorType::InitType &Grid,
                                        const aol::Vector<RealType> &UFixed,
                                        const aol::Vector<RealType> &UTemplate,
                                        const aol::MultiVector<RealType> &GradUTemplate ) :
    _grid( Grid ),
    _UFixed( UFixed ),
    _UTemplate( UTemplate ),
    _gradUTemplate( GradUTemplate ) {}

  void applyAdd( const aol::MultiVector<RealType> &Arg, BlockOpType &Dest ) const {
    qc::MismatchEnergySecondVariationWRTDispAssembly<ConfiguratorType>( _grid,  _UFixed, _UTemplate, _gradUTemplate, Arg ).assembleAddMatrix( Dest );
  }

  void applyAddMultiple( const aol::MultiVector<RealType> &Arg, BlockOpType &Dest, const RealType Factor = aol::NumberTrait<RealType>::One ) const {
    qc::MismatchEnergySecondVariationWRTDispAssembly<ConfiguratorType>( _grid,  _UFixed, _UTemplate, _gradUTemplate, Arg ).assembleAddMatrix( Dest, Factor );
  }
};

//! Interface for linear FE operators with basis functions deformed by two different deformations \f$\Phi_1, \Phi_2\f$.
/**
 *
 *  \author Effland
 */
template <typename ConfiguratorType, typename Imp, typename _MatrixType = typename ConfiguratorType::MatrixType>
class FELinDeformOpInterfaceDifferentDeformations : public aol::Op<aol::Vector<typename ConfiguratorType::RealType> > {

protected:
  typedef typename ConfiguratorType::RealType RealType;

  aol::DeleteFlagPointer<const ConfiguratorType> _config;
  mutable _MatrixType *_mat;
  aol::OperatorType _opType;
  // the displacement corresponding to the deformation \Phi
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement1;
  const aol::DiscreteVectorFunctionDefault<ConfiguratorType,ConfiguratorType::Dim> _displacement2;

public:
  explicit FELinDeformOpInterfaceDifferentDeformations ( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &Displacement1, const aol::MultiVector<RealType> &Displacement2, aol::OperatorType OpType = aol::ONTHEFLY ) :
  _config ( new ConfiguratorType(Grid), true ),
  _mat ( NULL ),
  _opType ( OpType ),
  _displacement1( Grid, Displacement1 ),
  _displacement2( Grid, Displacement2 ) {}

  virtual ~FELinDeformOpInterfaceDifferentDeformations( ) {
    delete _mat;
  }

  //! Clears the assembled matrix.
  void reset( ) {
    if ( _mat )
      delete _mat;
  }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    switch ( _opType ) {
    case aol::ONTHEFLY:
      multiplyOnTheFly ( Arg, Dest );
      break;
    case aol::ASSEMBLED:
      if ( !_mat )
        assembleMatrix();
      _mat->applyAdd ( Arg, Dest );
      break;
    default:
      throw aol::UnimplementedCodeException ( "FELinDeformOpInterface::applyAdd: unknown opType", __FILE__, __LINE__ );
    };
  }

  _MatrixType& getMatrix( ) {
    if ( !_mat ) {
      _mat = _config->createNewMatrix( );
      assembleMatrix( );
    }
    //return dynamic_cast<typename ConfiguratorType::MatrixType&>(*_mat);
    return *_mat;
  }

protected:
  void multiplyOnTheFly ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {

    typename ConfiguratorType::ElementType transformedEl1, transformedEl2;
    typename ConfiguratorType::DomVecType  transformedCoord1, transformedCoord2;
    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;
    int globalDofs1[ ConfiguratorType::maxNumLocalDofs ], globalDofs2[ ConfiguratorType::maxNumLocalDofs ];

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const typename IteratorType::EndType end_it = _config->end();

    // traverse the elements of the grid and on each element the quadrature points
    for ( IteratorType it = _config->begin(); it != end_it; ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config->getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; q++ ) {
        // compute \Phi(x) ("transformedEl", "transformedCoord")
        typename ConfiguratorType::DomVecType localCoord( bfs.getRefCoord ( q ) );
        qc::transformAndClipCoord<ConfiguratorType>( *_config, _displacement1, *it, q, localCoord, transformedEl1, transformedCoord1 );
        qc::transformAndClipCoord<ConfiguratorType>( *_config, _displacement2, *it, q, localCoord, transformedEl2, transformedCoord2 );

        // assemble the local matrix belonging to the current quadrature point
        this->asImp().prepareLocalMatrix ( *it, q, localCoord, transformedEl1, transformedCoord1, transformedEl2, transformedCoord2, localMatrix );

        const int numLocalDofs1 = _config->getNumLocalDofs ( transformedEl1 );
        for ( int i = 0; i < numLocalDofs1; ++i )
          globalDofs1[ i ] = _config->localToGlobal ( transformedEl1, i );

        const int numLocalDofs2 = _config->getNumLocalDofs ( transformedEl2 );
          for ( int i = 0; i < numLocalDofs2; ++i )
            globalDofs2[ i ] = _config->localToGlobal ( transformedEl2, i );

        // add the locally computed value to the global result
        for ( int i = 0; i < numLocalDofs1; ++i )
          for ( int j = 0; j < numLocalDofs2; ++j )
            Dest[ globalDofs1[ i ] ] += localMatrix[ i ][ j ] * Arg[ globalDofs2[ j ] ] ;
      }
    }
  }

  void assembleMatrix( ) const {
    if ( _mat )
      delete _mat;
    _mat = _config->createNewMatrix( );
    assembleAddMatrix ( *_mat );
  }

public:
  /** (this assembled matrix * Factor) is added to Mat  */
  template <typename MatrixType>
  void assembleAddMatrix ( MatrixType &Mat, const RealType Factor = aol::NumberTrait<RealType>::one ) const {

    typename ConfiguratorType::ElementType transformedEl1, transformedEl2;
    typename ConfiguratorType::DomVecType  transformedCoord1, transformedCoord2;
    aol::Mat<ConfiguratorType::maxNumLocalDofs,ConfiguratorType::maxNumLocalDofs,RealType> localMatrix;
    int globalDofs1[ ConfiguratorType::maxNumLocalDofs ], globalDofs2[ ConfiguratorType::maxNumLocalDofs ];

    typedef typename ConfiguratorType::ElementIteratorType IteratorType;
    const typename IteratorType::EndType end_it = _config->end();

    // traverse the elements of the grid and on each element the quadrature points
    for ( IteratorType it = _config->begin(); it != end_it; ++it ) {
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config->getBaseFunctionSet ( *it );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; q++ ) {
        typename ConfiguratorType::DomVecType localCoord( bfs.getRefCoord ( q ) );
        qc::transformAndClipCoord<ConfiguratorType>( *_config, _displacement1, *it, q, localCoord, transformedEl1, transformedCoord1 );
        qc::transformAndClipCoord<ConfiguratorType>( *_config, _displacement2, *it, q, localCoord, transformedEl2, transformedCoord2 );

        // assemble the local matrix belonging to the current quadrature point
        this->asImp().prepareLocalMatrix ( *it, q, localCoord, transformedEl1, transformedCoord1, transformedEl2, transformedCoord2, localMatrix );

        const int numLocalDofs1 = _config->getNumLocalDofs ( transformedEl1 );
        for ( int i = 0; i < numLocalDofs1; ++i )
          globalDofs1[ i ] = _config->localToGlobal ( transformedEl1, i );

        const int numLocalDofs2 = _config->getNumLocalDofs ( transformedEl2 );
          for ( int i = 0; i < numLocalDofs2; ++i )
            globalDofs2[ i ] = _config->localToGlobal ( transformedEl2, i );

        // add the locally assembled matrix to the globally assembled matrix
        for ( int i = 0; i < numLocalDofs1; ++i )
          for ( int j = 0; j < numLocalDofs2; ++j )
            Mat.add ( globalDofs1[ i ], globalDofs2[ j ], Factor * localMatrix [ i ][ j ] );
      }
    }
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

/**
 * Same as FELinDeformScalarWeightedMassInterface, but with two different deformations
 *
 *  \author Effland
 */
template <typename ConfiguratorType, typename Imp, typename MatrixType = typename ConfiguratorType::MatrixType>
class FELinDeformScalarWeightedMassInterfaceDifferentDeformations :
    public FELinDeformOpInterfaceDifferentDeformations< ConfiguratorType, FELinDeformScalarWeightedMassInterfaceDifferentDeformations<ConfiguratorType, Imp, MatrixType>, MatrixType > {

    protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

    public:
  explicit FELinDeformScalarWeightedMassInterfaceDifferentDeformations ( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &Displacement1, const aol::MultiVector<RealType> &Displacement2, aol::OperatorType OpType = aol::ONTHEFLY ) :
  FELinDeformOpInterfaceDifferentDeformations<ConfiguratorType, FELinDeformScalarWeightedMassInterfaceDifferentDeformations<ConfiguratorType, Imp, MatrixType>, MatrixType > ( Grid, Displacement1, Displacement2, OpType ) {}

  //! This function has to be provided in the implementation (derived class) of the interface.
  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const DomVecType &LocalCoord,
      const typename ConfiguratorType::ElementType &TransformedEl1, const DomVecType &TransformedLocalCoord1,
      const typename ConfiguratorType::ElementType &TransformedEl2, const DomVecType &TransformedLocalCoord2 ) const {
    return this->asImp().getCoeff ( El, QuadPoint, LocalCoord, TransformedEl1, TransformedLocalCoord1, TransformedEl2, TransformedLocalCoord2 );
  }

  //! Performs the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, int QuadPoint, const DomVecType &LocalCoord,
      const typename ConfiguratorType::ElementType &TransformedEl1, const DomVecType &TransformedLocalCoord1,
      const typename ConfiguratorType::ElementType &TransformedEl2, const DomVecType &TransformedLocalCoord2,
      aol::Mat<ConfiguratorType::maxNumLocalDofs, ConfiguratorType::maxNumLocalDofs, RealType> &LocalMatrix ) const {

    const RealType coeff = getCoeff ( El, QuadPoint, LocalCoord, TransformedEl1, TransformedLocalCoord1, TransformedEl2, TransformedLocalCoord2 );
    aol::Vec<ConfiguratorType::maxNumLocalDofs,RealType> basisFunctionValues1, basisFunctionValues2;

    const int numDofs1 = this->_config->getNumLocalDofs ( TransformedEl1 );
    const typename ConfiguratorType::BaseFuncSetType &tbfs1 = this->_config->getBaseFunctionSet ( TransformedEl1 );
    for ( int i = 0; i < numDofs1; i++ )
      basisFunctionValues1[i] = tbfs1.evaluate( i, TransformedLocalCoord1 );

    // compute the basis function values belonging to vectors multiplied from the right
    const int numDofs2 = this->_config->getNumLocalDofs ( TransformedEl2 );
    const typename ConfiguratorType::BaseFuncSetType &tbfs2 = this->_config->getBaseFunctionSet ( TransformedEl2 );
    for ( int i = 0; i < numDofs2; i++ )
      basisFunctionValues2[i] = tbfs2.evaluate( i, TransformedLocalCoord2 );

    // assemble the local mass matrix
    const typename ConfiguratorType::BaseFuncSetType &bfs = this->_config->getBaseFunctionSet ( El );
    LocalMatrix.makeTensorProduct( basisFunctionValues1, basisFunctionValues2 );
    LocalMatrix *= coeff * bfs.getWeight ( QuadPoint ) * this->_config->vol ( El );
  }

    protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }
};

/**
 * Same as DeformMassOp, but with two different deformations
 *
 *  \author Effland
 */
template <typename ConfiguratorType, typename MatrixType = typename ConfiguratorType::MatrixType>
class DeformMassDifferentDeformationsOp :
    public FELinDeformScalarWeightedMassInterfaceDifferentDeformations<ConfiguratorType, DeformMassDifferentDeformationsOp<ConfiguratorType, MatrixType>, MatrixType> {

    protected:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

    public:
  explicit DeformMassDifferentDeformationsOp ( const typename ConfiguratorType::InitType &Grid, const aol::MultiVector<RealType> &Displacement1, const aol::MultiVector<RealType> &Displacement2, aol::OperatorType OpType = aol::ONTHEFLY ) :
  FELinDeformScalarWeightedMassInterfaceDifferentDeformations<ConfiguratorType, DeformMassDifferentDeformationsOp<ConfiguratorType, MatrixType>, MatrixType> ( Grid, Displacement1, Displacement2, OpType ) {}

  inline RealType getCoeff ( const typename ConfiguratorType::ElementType &/*El*/, int /*QuadPoint*/, const DomVecType &/*LocalCoord*/,
      const typename ConfiguratorType::ElementType &/*TransformedEl1*/, const DomVecType &/*TransformedLocalCoord1*/,
      const typename ConfiguratorType::ElementType &/*TransformedEl2*/, const DomVecType &/*TransformedLocalCoord2*/) const {
    return aol::ZOTrait<RealType>::one;
  }
};

} // end namespace qc

#endif // __DEFORMATIONS
