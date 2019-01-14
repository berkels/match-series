#include <auxiliary.h>

#include <scalarArray.h>
#include <bzipiostream.h>

#include <multiArray.h>

template <typename RealType>
void qc::TransformFunction<RealType, qc::QC_2D>::transform ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest, qc::BitArray<qc::QC_2D> &ValuesSet ) const {
  aol::MultiVector<RealType> extendImage( Arg, aol::STRUCT_COPY );
  //extendImage.setZero();
  transform( Arg, Dest, ValuesSet, extendImage );
}

template <typename RealType>
void qc::TransformFunction<RealType, qc::QC_2D>::transform ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest, qc::BitArray<qc::QC_2D> &ValuesSet ) const {
  aol::MultiVector<RealType> arg, dest;
  arg.appendReference( Arg );
  dest.appendReference( Dest );
  transform( arg, dest, ValuesSet );
}

template <typename RealType>
void qc::TransformFunction<RealType, qc::QC_2D>::transform ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest, qc::BitArray<qc::QC_2D> &ValuesSet, const aol::Vector<RealType> &ExtendImage ) const {
  aol::MultiVector<RealType> arg, dest, extendImage;
  arg.appendReference( Arg );
  dest.appendReference( Dest );
  extendImage.appendReference( ExtendImage );
  transform( arg, dest, ValuesSet, extendImage );
}

template <typename RealType>
void qc::TransformFunction<RealType, qc::QC_2D>::transform ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest, qc::BitArray<qc::QC_2D> &ValuesSet, const aol::MultiVector<RealType> &ExtendImage ) const {
  ValuesSet.reallocate( _grid.getNumX(), _grid.getNumY() );

  if ( !_def ) {
    throw aol::Exception ( "set deformation first\n", __FILE__, __LINE__ );
  }

  Dest = ExtendImage;

  const int numComponents = Arg.numComponents();

  std::vector<qc::Array<RealType>* > ADestVec ( numComponents );
  for ( int i = 0; i < numComponents; i++ ) {
    ADestVec[i] = new qc::Array<RealType> ( Dest[i], _grid );
  }

  std::vector<const qc::Array<RealType>* > AArgVec ( numComponents );
  for ( int i = 0; i < numComponents; i++ ) {
    AArgVec[i] = new const qc::Array<RealType> ( Arg[i], _grid );
  }

  qc::Array<RealType> defx ( ( *_def ) [0], _grid );
  qc::Array<RealType> defy ( ( *_def ) [1], _grid );

  qc::RectangularGrid<qc::QC_2D>::OldAllElementIterator it;

  int offsets[2][4][2] = { { {0, 0}, {1, 0}, {0, 1}, {1, 1} }, { {0, 1}, {0, 0}, {1, 1}, {1, 0} } };
  // The value of bit selects either criss or cross triangularization of the rectangular grid elements when computing the inverse.
  int bit = 0;
  int pts[4][2];
  RealType dx[4], dy[4];
  aol::Vec2<RealType> coords[3];
  RealType *values = new RealType[numComponents*3];
  for ( it = _grid.begin_it; it != _grid.end_it; ++it ) {
    //bit = aol::Abs( 1 - bit );
    for ( int i = 0; i < 4; i++ ) {
      pts[i][0] = it->x() + offsets[bit][i][0];
      pts[i][1] = it->y() + offsets[bit][i][1];

      dx[i] = static_cast< RealType > ( pts[i][0] ) + defx.get ( pts[i][0], pts[i][1] ) / static_cast< RealType > (_grid.H());
      dy[i] = static_cast< RealType > ( pts[i][1] ) + defy.get ( pts[i][0], pts[i][1] ) / static_cast< RealType > (_grid.H());

      dx[i] = aol::Min ( dx[i], static_cast< RealType > ( _grid.getNumX() - 1 ) );
      dy[i] = aol::Min ( dy[i], static_cast< RealType > ( _grid.getNumY() - 1 ) );
      dx[i] = aol::Max ( dx[i], static_cast< RealType > ( 0.0 ) );
      dy[i] = aol::Max ( dy[i], static_cast< RealType > ( 0.0 ) );
    }

    for( int i = 0; i < 3; i++ ){
      coords[i][0] = dx[i];
      coords[i][1] = dy[i];
    }

    if ( ( _mask == NULL ) || ( _mask->get( pts[1][0], pts[1][1] ) && _mask->get( pts[2][0], pts[2][1] ) ) ) {
      if ( ( _mask == NULL ) || ( _mask->get( pts[0][0], pts[0][1] ) ) ) {
        for ( int k = 0; k < numComponents; k++ ) {
          values[k*3] = AArgVec[k]->get ( pts[0][0], pts[0][1] );
          values[k*3+1] = AArgVec[k]->get ( pts[1][0], pts[1][1] );
          values[k*3+2] = AArgVec[k]->get ( pts[2][0], pts[2][1] );
        }
        interpolateTriangle ( coords, values, ADestVec, numComponents, ValuesSet );
      }

      coords[0][0] = dx[3];
      coords[0][1] = dy[3];
      if ( ( _mask == NULL ) || ( _mask->get( pts[3][0], pts[3][1] ) ) ) {
        for ( int k = 0; k < numComponents; k++ ) {
          values[k*3] = AArgVec[k]->get ( pts[3][0], pts[3][1] );
        }
        interpolateTriangle ( coords, values, ADestVec, numComponents, ValuesSet );
      }
    }
  }
  for ( int i = 0; i < numComponents; i++ ) {
    delete ADestVec[i];
    delete AArgVec[i];
  }
  delete[] values;
}

template <typename RealType>
void qc::TransformFunction<RealType, qc::QC_2D>::interpolateTriangle ( const aol::Vec2<RealType> Coords[3], const RealType *Values, std::vector<qc::Array<RealType>* > &Array, const int NumComponents, qc::BitArray<qc::QC_2D> &ValuesSet ) const {
  int miny = _grid.getNumY(), maxy = -1;
  int minx = _grid.getNumX(), maxx = -1;

  for ( int loc = 0; loc < 3; loc++ ) {
    int qx = static_cast< int > ( Coords[loc][0] );
    int qy = static_cast< int > ( Coords[loc][1] );

    maxy = aol::Max ( qy, maxy );
    if ( qy < Coords[loc][1] ) qy++;
    miny = aol::Min ( qy, miny );

    maxx = aol::Max ( qx, maxx );
    if ( qx < Coords[loc][0] ) qx++;
    minx = aol::Min ( qx, minx );
  }

  // check bounds
  minx = aol::Max ( 0, minx );
  miny = aol::Max ( 0, miny );
  maxx = aol::Min ( _grid.getNumX() - 1, maxx );
  maxy = aol::Min ( _grid.getNumY() - 1, maxy );

  aol::Matrix22<RealType> mat, inv;
  mat[0][0] = Coords[0][0] - Coords[2][0];
  mat[1][0] = Coords[0][1] - Coords[2][1];

  mat[0][1] = Coords[1][0] - Coords[2][0];
  mat[1][1] = Coords[1][1] - Coords[2][1];

  const RealType matdet = mat.det();

  if ( matdet == 0 || ! aol::isFinite ( matdet ) )
    return;

  mat.make_inverse( inv );

  aol::Vec2<RealType> lambda;
  aol::Vec2<RealType> r;
  for ( int x = minx; x <= maxx; x++ ) {
    for ( int y = miny; y <= maxy; y++ ) {
      if ( ValuesSet.get ( x, y ) ) {} else {
        const RealType sx = static_cast<RealType> ( x );
        const RealType sy = static_cast<RealType> ( y );
        r.set ( sx - Coords[2][0], sy - Coords[2][1] );
        inv.mult ( r, lambda );

        if ( lambda[0] < -1.e-6 || lambda[1] < -1.e-6 || ( lambda[0] + lambda[1] ) > 1. + 1.e-6 || ! aol::isFinite ( mat.det() ) ) {} else {
          if ( ! ( lambda[0] < -1.e-6 || lambda[1] < -1.e-6 || ( lambda[0] + lambda[1] ) > 1. + 1.e-6 ) ) {
            ValuesSet.set ( x, y, true );
          }
          for ( int k = 0; k < NumComponents; k++ ) {
            RealType v = lambda[0] * Values[k*3] + lambda[1] * Values[k*3+1] + ( aol::ZOTrait<RealType>::one - lambda[0] - lambda[1] ) * Values[k*3+2];
            if ( aol::isFinite ( v ) ) {
              Array[k]->set ( x, y, v );
            }
          }
        }
      }
    }
  }
}

/**
 * Does the same as "transform", only uses multilinear instead of only affine interpolation, which is more appropriate for the quocmeshes.
 */
template <typename RealType>
void qc::TransformFunction<RealType, qc::QC_2D>::transformMultiLin ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest, qc::BitArray<qc::QC_2D> &ValuesSet, const aol::MultiVector<RealType> &ExtendImage ) const {

  if ( !_def )
    throw aol::Exception ( "set deformation first\n", __FILE__, __LINE__ );

  // initialize the results
  Dest = ExtendImage;
  ValuesSet.reallocate( _grid.getNumX(), _grid.getNumY() );

  // produce arrays of ScalarArray<QC_2D> from Arg and Dest and the deformation (so that e.g. interpolation methods can be used on them)
  int numComponents = Arg.numComponents();
  std::vector<qc::ScalarArray<RealType, qc::QC_2D>* > destPointerVec ( numComponents );
  std::vector<const qc::ScalarArray<RealType, qc::QC_2D>* > argPointerVec ( numComponents );
  for ( int i = 0; i < numComponents; i++ ) {
    destPointerVec[i] = new qc::ScalarArray<RealType, qc::QC_2D> ( Dest[i], _grid, aol::FLAT_COPY );
    argPointerVec[i] = new const qc::ScalarArray<RealType, qc::QC_2D> ( Arg[i], _grid, aol::FLAT_COPY );
  }
  qc::MultiArray<RealType,qc::QC_2D,qc::QC_2D> def ( _grid, *_def );

  // compute the inverse deformation
  int offsets[4][2] = { {0, 0}, {1, 0}, {0, 1}, {1, 1} };
  int elementNodes[4][2];
  aol::MultiVector<RealType> preImages( 4, 2 ), images( 4, 2 );
  qc::RectangularGrid<qc::QC_2D>::OldAllElementIterator it;
  // iterate through all images of the finite elements
  for ( it = _grid.begin_it; it != _grid.end_it; ++it ) {
    for ( int i = 0; i < 4; i++ ) {
      elementNodes[i][0] = it->x() + offsets[i][0];
      elementNodes[i][1] = it->y() + offsets[i][1];

      preImages[i][0] = static_cast< RealType > ( elementNodes[i][0] );
      preImages[i][1] = static_cast< RealType > ( elementNodes[i][1] );

      images[i][0] = static_cast< RealType > ( elementNodes[i][0] ) + def[0].get ( elementNodes[i][0], elementNodes[i][1] ) / static_cast< RealType > ( _grid.H() );
      images[i][1] = static_cast< RealType > ( elementNodes[i][1] ) + def[1].get ( elementNodes[i][0], elementNodes[i][1] ) / static_cast< RealType > ( _grid.H() );
    }

    interpolateQuad ( preImages, images, argPointerVec, destPointerVec, ValuesSet );
  }

  // free storage
  for ( int i = 0; i < numComponents; i++ ) {
    delete destPointerVec[i];
    delete argPointerVec[i];
  }
}

/**
 * For the four vertices of a square quocmesh element, passed as "NodalPreImages", and their positions after deformation, passed
 * as "NodalImages", this method interpolates the original position (i.e. before deformation) of "Image" and puts it into "PreImage".
 */
template <typename RealType>
inline void qc::TransformFunction<RealType, qc::QC_2D>::getLocalPreImage ( const aol::Vec2<RealType> &Image,
                                                                           const aol::MultiVector<RealType> &NodalPreImages,
                                                                           const aol::MultiVector<RealType> &NodalImages,
                                                                           aol::Vec2<RealType> &PreImage ) const {

  // let \f$ f=(\xi-\xi_0)(\eta-\eta_0)(x,y)_{11}^T+(\xi-\xi_0)(\eta-\eta_1)(x,y)_{10}^T+(\xi-\xi_1)(\eta-\eta_0)(x,y)_{01}^T+(\xi-\xi_1)(\eta-\eta_1)(x,y)_{00}^T-(x,y)^T\f$,
  // where we try to find \f$(\xi,\eta)\f$ such that \f$f=0\f$, and where Image=\f$(x,y)\f$, NodalPreImages=\f$[(\xi_0,\eta_0),(\xi_1,\eta_0),(\xi_0,\eta_1),(\xi_1,\eta_1)]\f$,
  // NodalImages=\f$[(x,y)_{00},(x,y)_{10},(x,y)_{01},(x,y)_{11}]\f$

  // perform few steps of Newton's method
  aol::Vec2<RealType> f;
  do {
    // abbreviations
    RealType xiL = PreImage[0] - NodalPreImages[0][0],
             xiU = PreImage[0] - NodalPreImages[3][0],
             etaL = PreImage[1] - NodalPreImages[0][1],
             etaU = PreImage[1] - NodalPreImages[3][1];

    // compute the derivative of \f$f\f$ wrt \f$(\xi,\eta)\f$ and its inverse
    aol::Matrix22<RealType> df( etaL * ( NodalImages[3][0] - NodalImages[2][0] ) + etaU * ( NodalImages[0][0] - NodalImages[1][0] ),
                                xiL * ( NodalImages[3][0] - NodalImages[1][0] ) + xiU * ( NodalImages[0][0] - NodalImages[2][0] ),
                                etaL * ( NodalImages[3][1] - NodalImages[2][1] ) + etaU * ( NodalImages[0][1] - NodalImages[1][1] ),
                                xiL * ( NodalImages[3][1] - NodalImages[1][1] ) + xiU * ( NodalImages[0][1] - NodalImages[2][1] ) ),
                            invDf;
    invDf.makeInverse( df );

    // compute \f$f\f$
    f[0] = xiL * ( etaL * NodalImages[3][0] - etaU * NodalImages[1][0] ) + xiU * ( etaU * NodalImages[0][0] - etaL * NodalImages[2][0] ) - Image[0];
    f[1] = xiL * ( etaL * NodalImages[3][1] - etaU * NodalImages[1][1] ) + xiU * ( etaU * NodalImages[0][1] - etaL * NodalImages[2][1] ) - Image[1];

    // compute \f$(\xi,\eta)_{i+1}^T=(\xi,\eta)_i^T-Df^{-1}f\f$
    PreImage -= invDf * f;
  } while ( f.normSqr() > 1.e-12 );
}

/**
 * For the four vertices of a square quocmesh element, passed as "PreImages", and their positions after deformation, passed
 * as "Images", this method multilinearly interpolates the original position (i.e. before deformation) of each grid node of
 * "Dest" and gives it the corresponding value from "Arg", if the original position lies within the square defined by "PreImages".
 * Also, the corresponding entry in "ValuesSet" is then set to true to signify that the value has been set.
 */
template <typename RealType>
void qc::TransformFunction<RealType, qc::QC_2D>::interpolateQuad ( const aol::MultiVector<RealType> &PreImages,
                                                                   const aol::MultiVector<RealType> &Images,
                                                                   const std::vector<const qc::ScalarArray<RealType, qc::QC_2D>* > &Arg,
                                                                   const std::vector<qc::ScalarArray<RealType, qc::QC_2D>* > &Dest,
                                                                   qc::BitArray<qc::QC_2D> &ValuesSet ) const {

  // compute the range [minX,maxX]x[minY,maxY], in which the image of the convex hull of PreImages lies
  int minX = static_cast< int > ( ceil( aol::Min( aol::Min( Images[0][0], Images[1][0] ), aol::Min( Images[2][0], Images[3][0] ) ) ) ),
      maxX = static_cast< int > ( floor( aol::Max( aol::Max( Images[0][0], Images[1][0] ), aol::Max( Images[2][0], Images[3][0] ) ) ) ),
      minY = static_cast< int > ( ceil( aol::Min( aol::Min( Images[0][1], Images[1][1] ), aol::Min( Images[2][1], Images[3][1] ) ) ) ),
      maxY = static_cast< int > ( floor( aol::Max( aol::Max( Images[0][1], Images[1][1] ), aol::Max( Images[2][1], Images[3][1] ) ) ) );

  // check bounds
  minX = aol::Max ( 0, minX );
  minY = aol::Max ( 0, minY );
  maxX = aol::Min ( _grid.getNumX() - 1, maxX );
  maxY = aol::Min ( _grid.getNumY() - 1, maxY );

  // iterate over all nodes in the range, in which the image of the convex hull of PreImages lies
  for ( int x = minX; x <= maxX; x++ )
    for ( int y = minY; y <= maxY; y++ )
      if ( ValuesSet.get ( x, y ) ) {} else {

        // obtain the preimage of (x,y)
        aol::Vec2<RealType> image( static_cast<RealType> ( x ), static_cast<RealType> ( y ) ),
                            preImage( PreImages[0][0], PreImages[0][1] );
        getLocalPreImage( image, PreImages, Images, preImage );

        // if the preimage of (x,y) lies within the convex hull of PreImages, compute the corresponding interpoland values
        RealType delta = .01;
        if ( ( preImage[0] > PreImages[0][0] - delta ) && ( preImage[0] < PreImages[3][0] + delta )
            && ( preImage[1] > PreImages[0][1] - delta ) && ( preImage[1] < PreImages[3][1] + delta ) ) {
          ValuesSet.set( x, y, true );
          for ( unsigned int k = 0; k < Arg.size(); k++ )
            Dest[k]->set( x, y, Arg[k]->interpolate( preImage ) );
        }
      }
}

template <typename RealType>
void qc::TransformFunction<RealType, qc::QC_3D>::transform ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest, qc::BitArray<qc::QC_3D> &ValuesSet ) const {
  aol::MultiVector<RealType> extendImage( Arg, aol::STRUCT_COPY );
  //extendImage.setZero();
  transform( Arg, Dest, ValuesSet, extendImage );
}

template <typename RealType>
void qc::TransformFunction<RealType, qc::QC_3D>::transform ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest, qc::BitArray<qc::QC_3D> &ValuesSet ) const {
  aol::MultiVector<RealType> arg, dest;
  arg.appendReference( Arg );
  dest.appendReference( Dest );
  transform( arg, dest, ValuesSet );
}

template <typename RealType>
void qc::TransformFunction<RealType, qc::QC_3D>::transform ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest, qc::BitArray<qc::QC_3D> &ValuesSet, const aol::Vector<RealType> &ExtendImage ) const {
  aol::MultiVector<RealType> arg, dest, extendImage;
  arg.appendReference( Arg );
  dest.appendReference( Dest );
  extendImage.appendReference( ExtendImage );
  transform( arg, dest, ValuesSet, extendImage );
}

template <typename RealType>
void qc::TransformFunction<RealType, qc::QC_3D>::transform ( const aol::MultiVector<RealType> &Arg, aol::MultiVector<RealType> &Dest, qc::BitArray<qc::QC_3D> &ValuesSet, const aol::MultiVector<RealType> &ExtendImage ) const {
  ValuesSet.reallocate( _grid.getNumX(), _grid.getNumY(), _grid.getNumZ() );

  if ( !_def ) {
    throw aol::Exception ( "set deformation first\n", __FILE__, __LINE__ );
  }

  Dest = ExtendImage;

  const int numComponents = Arg.numComponents();

  std::vector<qc::Array<RealType>* > ADestVec ( numComponents );
  for ( int i = 0; i < numComponents; i++ ) {
    ADestVec[i] = new qc::Array<RealType> ( Dest[i], _grid );
  }

  std::vector<const qc::Array<RealType>* > AArgVec ( numComponents );
  for ( int i = 0; i < numComponents; i++ ) {
    AArgVec[i] = new const qc::Array<RealType> ( Arg[i], _grid );
  }

  qc::Array<RealType> defx ( ( *_def ) [0], _grid );
  qc::Array<RealType> defy ( ( *_def ) [1], _grid );
  qc::Array<RealType> defz ( ( *_def ) [2], _grid );

  qc::RectangularGrid<qc::QC_3D>::OldAllElementIterator it;

  int offsets[8][3] = { {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1} };
  int pts[8][3];
  RealType dx[8], dy[8], dz[8];
  aol::Vec3<RealType> coords[4];
  RealType *values = new RealType[numComponents*4];
  for ( it = _grid.begin_it; it != _grid.end_it; ++it ) {
    for ( int i = 0; i < 8; i++ ) {
      pts[i][0] = it->x() + offsets[i][0];
      pts[i][1] = it->y() + offsets[i][1];
      pts[i][2] = it->z() + offsets[i][2];

      dx[i] = static_cast< RealType > ( pts[i][0] ) + defx.get ( pts[i][0], pts[i][1], pts[i][2] ) / static_cast< RealType > (_grid.H());
      dy[i] = static_cast< RealType > ( pts[i][1] ) + defy.get ( pts[i][0], pts[i][1], pts[i][2] ) / static_cast< RealType > (_grid.H());
      dz[i] = static_cast< RealType > ( pts[i][2] ) + defz.get ( pts[i][0], pts[i][1], pts[i][2] ) / static_cast< RealType > (_grid.H());

      dx[i] = aol::Min ( dx[i], static_cast< RealType > ( _grid.getNumX() - 1 ) );
      dy[i] = aol::Min ( dy[i], static_cast< RealType > ( _grid.getNumY() - 1 ) );
      dz[i] = aol::Min ( dz[i], static_cast< RealType > ( _grid.getNumZ() - 1 ) );
      dx[i] = aol::Max ( dx[i], static_cast< RealType > ( 0.0 ) );
      dy[i] = aol::Max ( dy[i], static_cast< RealType > ( 0.0 ) );
      dz[i] = aol::Max ( dz[i], static_cast< RealType > ( 0.0 ) );
    }

    // The nodes of the cube are numbered from 0 to 7, the cube is cut into 6 tetrahedrons
    int nodenumbers[6][4] = { {0, 1, 2, 4}, {1, 2, 3, 4}, {2, 3, 4, 6}, {1, 3, 4, 5}, {3, 4, 5, 6}, {3, 5, 6, 7} };

    for ( int i = 0; i < 6; i++ ) {
      bool mask = true;
      for ( int j = 0; j < 4; j++ ) {
        const int index = nodenumbers[i][j];
        coords[j].set ( dx[index], dy[index], dz[index] );
        for ( int k = 0; k < numComponents; k++ ) {
          values[k*4+j] = AArgVec[k]->get ( it->x() + offsets[index][0], it->y() + offsets[index][1], it->z() + offsets[index][2] );
        }
        mask = ( _mask == NULL ) || ( mask && _mask->get( it->x() + offsets[index][0], it->y() + offsets[index][1], it->z() + offsets[index][2] ) );
      }
      if ( mask )
        interpolateTetrahedron ( coords, values, ADestVec, numComponents, ValuesSet );
    }
  }
  for ( int i = 0; i < numComponents; i++ ) {
    delete ADestVec[i];
    delete AArgVec[i];
  }
  delete[] values;
}

template <typename RealType>
void qc::TransformFunction<RealType, qc::QC_3D>::interpolateTetrahedron ( const aol::Vec3<RealType> Coords[4], const RealType *Values, std::vector<qc::Array<RealType>* > &Array, const int NumComponents, qc::BitArray<qc::QC_3D> &ValuesSet ) const {
  int minz = _grid.getNumZ(), maxz = -1;
  int miny = _grid.getNumY(), maxy = -1;
  int minx = _grid.getNumX(), maxx = -1;

  for ( int loc = 0; loc < 4; loc++ ) {
    int qx = static_cast< int > ( Coords[loc][0] );
    int qy = static_cast< int > ( Coords[loc][1] );
    int qz = static_cast< int > ( Coords[loc][2] );

    maxz = aol::Max ( qz, maxz );
    if ( qz < Coords[loc][2] ) qz++;
    minz = aol::Min ( qz, minz );

    maxy = aol::Max ( qy, maxy );
    if ( qy < Coords[loc][1] ) qy++;
    miny = aol::Min ( qy, miny );

    maxx = aol::Max ( qx, maxx );
    if ( qx < Coords[loc][0] ) qx++;
    minx = aol::Min ( qx, minx );
  }

  // check bounds
  minx = aol::Max ( 0, minx );
  miny = aol::Max ( 0, miny );
  minz = aol::Max ( 0, minz );
  maxx = aol::Min ( _grid.getNumX() - 1, maxx );
  maxy = aol::Min ( _grid.getNumY() - 1, maxy );
  maxz = aol::Min ( _grid.getNumZ() - 1, maxz );

  aol::Matrix33<RealType> mat, inv;
  for ( int i = 0; i < 3; i ++ ) {
    for ( int j = 0; j < 3; j ++ ) {
      mat[i][j] = Coords[j][i] - Coords[3][i];
    }
  }
  inv.makeCofactorMatrix ( mat );
  inv.transpose();
  const RealType matdet = mat.det();

  if ( matdet == 0. || ! aol::isFinite ( mat.det() ) )
    return;

  inv /= matdet;

  for ( int x = minx; x <= maxx; x++ ) {
    for ( int y = miny; y <= maxy; y++ ) {
      for ( int z = minz; z <= maxz; z++ ) {
        if ( ValuesSet.get ( x, y, z ) ) {} else {
          RealType sx = static_cast<RealType> ( x );
          RealType sy = static_cast<RealType> ( y );
          RealType sz = static_cast<RealType> ( z );
          aol::Vec3<RealType> lambda;
          aol::Vec3<RealType> r ( sx - Coords[3][0], sy - Coords[3][1], sz - Coords[3][2] );

          inv.mult ( r, lambda );

          // If the transformation is not invertible, it seems to be possible, that there are points, which aren't contained in one of the tetrahedrons. So it's necessary to set values also in the neighborhood of a given tetrahedron and mark those points, which are found in a given tetrahedron, to prevent overwriting of correct values.
          if ( lambda[0] < -5.e-1 || lambda[1] < -5.e-1 || lambda[2] < -5.e-1 || ( lambda[0] + lambda[1] + lambda[2] ) > 1. + 5.e-1 ) {} else {
            if ( ! ( lambda[0] < -1.e-6 || lambda[1] < -1.e-6 || lambda[2] < -1.e-6 || ( lambda[0] + lambda[1] + lambda[2] ) > 1. + 1.e-6 ) ) {
              ValuesSet.set ( x, y, z, true );
            }
            for ( int k = 0; k < NumComponents; k++ ) {
              RealType v = lambda[0] * Values[k*4] + lambda[1] * Values[k*4+1] + lambda[2] * Values[k*4+2] + ( aol::ZOTrait<RealType>::one - lambda[0] - lambda[1] - lambda[2] ) * Values[k*4+3];
              if ( aol::isFinite ( v ) ) {
                Array[k]->set ( x, y, z, v );
              }
            }
          }
        }
      }
    }
  }
}

int qc::getGridLevelFromArrayFile ( const string &ArrayFileName ) {
  aol::Bzipifstream file ( ArrayFileName.c_str() );
  ArrayHeader header;
  ReadArrayHeader ( file, header );
  switch ( header.numX ) {
  case 3: {
    return 1;
    break;
  }
  case 5: {
    return 2;
    break;
  }
  case 9: {
    return 3;
    break;
  }
  case 17: {
    return 4;
    break;
  }
  case 33: {
    return 5;
    break;
  }
  case 65: {
    return 6;
    break;
  }
  case 129: {
    return 7;
    break;
  }
  case 257: {
    return 8;
    break;
  }
  case 513: {
    return 9;
    break;
  }
  case 1025: {
    return 10;
    break;
  }
  case 2049: {
    return 11;
    break;
  }
  default: {
    char errorMessage[1024];
    sprintf( errorMessage, "Read header from %s:\nArray width of %d not equal to 2^d+1\nNote: qc::getGridLevelFromArrayFile can't be used on PNGs.\n", ArrayFileName.c_str(), header.numX );
    throw aol::Exception ( errorMessage, __FILE__, __LINE__ );
    return 0;
  }
  }
  // this point cannot be reached, return to prevent warning
  return 0;
}

aol::Vec3<int> qc::getSizeFromArrayFile ( const string &ArrayFileName ) {
  aol::Bzipifstream file ( ArrayFileName.c_str() );
  ArrayHeader header;
  ReadArrayHeader ( file, header );
  return aol::Vec3<int> ( header.numX, header.numY, header.numZ );
}

qc::Dimension qc::getDimensionFromArrayFile ( const string &ArrayFileName ) {
  // If the suffix is ".png" just assume that ArrayFileName actually is a PNG file.
  if ( aol::fileNameEndsWith ( ArrayFileName.c_str(), ".png" ) )
    return qc::QC_2D;

  aol::Bzipifstream file ( ArrayFileName.c_str() );
  ArrayHeader header;
  ReadArrayHeader ( file, header );

  switch ( header.magic[0] ) {
    case 'O':
      return qc::QC_1D;
      break;

    case 'P':
      return qc::QC_2D;
      break;

    case 'Q':
      return qc::QC_3D;
      break;

    default:
      throw aol::Exception( "Unexpected magic char in header encountered.\n", __FILE__, __LINE__);
      return qc::QC_1D;
  }
  // this point cannot be reached, return to prevent warning
  return qc::QC_1D;
}

template class qc::TransformFunction<float, qc::QC_2D>;
template class qc::TransformFunction<double, qc::QC_2D>;

template class qc::TransformFunction<float, qc::QC_3D>;
template class qc::TransformFunction<double, qc::QC_3D>;
