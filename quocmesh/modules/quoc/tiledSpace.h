#ifndef __TILEDSPACE_H
#define __TILEDSPACE_H

#include <op.h>
#include <aol.h>
#include <scalarArray.h>
#include <smallMat.h>

//---------------------------------------------------------------------------

namespace aol {

//   CLASS :  AffineTransformation
/**
 *   kapselt Matrix und Offset einer affinen Transformation.
 *
 *   \author von Deylen
 *
 *   bisher benutzt und getestet: nur im TiledSpace
 */

template <typename MatrixType, typename DomainType, typename RangeType>
class AffineTransformation : public Op<DomainType, RangeType> {
public:
  MatrixType matrix;
  RangeType  offset;

  // constructors:
  // (a) matrix and offset
  AffineTransformation ( const MatrixType & M, const RangeType & b ) :
      matrix ( M ),
      offset ( b ) {}

  // (b) only matrix (offset = 0 implicitly)
  AffineTransformation ( const MatrixType & M ) :
      matrix ( M ) {}

  // (c) nothing (transformation = identity)
  AffineTransformation() :
      matrix ( aol::ZOTrait<MatrixType>::one ) {}

  // applyAdd has to be re-implemented, as for every Op:
  void applyAdd ( const DomainType & x, RangeType & y ) const {
    matrix.multAdd ( x, y );
    y += offset;
  }

  bool isRigidBodyMotion () const {
    typename MatrixType::DataType det = matrix.det();
    return (fabs (det * det - ZOTrait<typename MatrixType::DataType>::one) < 1E-14 );
  }
};

} // end of namespace aol

//---------------------------------------------------------------------------

namespace qc {

//   CLASS :  TiledSpace
/**
 *   "gekachelter Raum", das heisst 2D- oder 3D-Raum, auf dem die Textur,
 *   die durch ein ScalarArray(2d/3d) gegeben ist,
 *   (a) periodisch fortgesetzt und
 *   (b) affin transformiert ist.
 *
 *   Die Textur sowie die anzuwendende affine Tranformation sind
 *   oeffentliche Variablen und koennen also vom Benutzer selbstaendig
 *   gesetzt werden. Er wird benutzt durch
 *
 *   get(x) = image(affineTransformation(x)).
 *
 *   \author von Deylen
 *   \date 2007-12-19
 *
 *   bisher benutzt und getestet: von Deylen: Projekt "lengthTest".
 */

template <typename ArrayType, typename MatrixType, typename VectorType>
class TiledSpace {
public:
  typedef aol::AffineTransformation<MatrixType, VectorType, VectorType> TransformType;
  typedef typename VectorType::DataType RealType;

  // public member variabled:
  ArrayType image;
  TransformType affineTransformation;

  // constructor:
  TiledSpace ( const ArrayType & _image, const TransformType & transf ) :
      image ( _image ),
      affineTransformation ( transf ) {}

  TiledSpace ( string fileName, const MatrixType & M ) :
      image ( fileName ),
      affineTransformation ( M ) {}

  TiledSpace ( const ArrayType & _image, RealType scaleFactor = aol::ZOTrait<RealType>::one ) :
      image ( _image ) {
    affineTransformation.matrix *= scaleFactor;
  }

  TiledSpace ( string fileName, RealType scaleFactor = aol::ZOTrait<RealType>::one ) :
      image ( fileName ) {
    affineTransformation.matrix *= scaleFactor;
  }

  // get function:
  RealType get ( const VectorType & x ) const {
      VectorType y;
      affineTransformation.applyAdd ( x, y );
      return image.interpolate_on01_periodic ( y );
    }
};

//---------------------------------------------------------------------------

} // end of namespace qc

#endif
