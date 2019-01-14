#include <quoc.h>
#include "estimatorTP.h"




template <class T>
EstimatorTP2d<T>::EstimatorTP2d ( int Size )
  : qc::ScalarArray<T, qc::QC_2D> ( Size, Size ) {

  maxLevel = qc::logBaseTwo ( Size );
}

template <class T>
EstimatorTP2d<T>::~EstimatorTP2d() {
}



/** *******************************************************
 * Method generates a saturated Error Array, by first
 * calculating the error and second saturate it.
 * ********************************************************/
//! method walks across the data-array and fills the own array
//! with the saturated error-entries
template <class T>
void EstimatorTP2d<T>::makeSaturatedErrorArray ( const qc::ScalarArray<T, qc::QC_2D>& meshData ) {
  int sizeOfMeshData = meshData.getNumX();
  int maxLevelOfMeshData = qc::logBaseTwo ( sizeOfMeshData );

  int level, levelDist;
  int i, j;
  int mx, my;
  T error;
  T max;

  // the error on the finest level is zero
  for ( i = 0; i < sizeOfMeshData; ++i ) {
    for ( j = 0; j < sizeOfMeshData; ++j ) {
      this->set ( i, j, T ( 0.0 ) );
    }
  }

  // calc the number of qc::Elements in one line in the max. level
  int numberOfElements = 1 << ( maxLevelOfMeshData - 1 );
  levelDist = 1;


  for ( i = 0; i < sizeOfMeshData; i++ ) { // x-direction
    for ( j = 0; j < sizeOfMeshData; j++ ) { // y-direction
      error = sqrt ( double ( aol::Sqr ( this->dx ( meshData, i, j ) ) + aol::Sqr ( this->dy ( meshData, i, j ) ) ) );
      this->set ( i, j, error );
    } // j
  } // i

  // -------------------- 2. Saturation ---------------------
  // now saturate the grid by adding the maximum of the five
  // circumjacent vertices of the finer level.

  numberOfElements = 1 << ( maxLevelOfMeshData - 1 );
  levelDist = 2;

  // again walk through the levels, maxLevelOfMeshData is the finest level
  for ( level = maxLevelOfMeshData - 1; level >= 0; --level ) {
    // walk through the vertices
    for ( i = 0; i < numberOfElements; ++i ) { // x-direction
      for ( j = 0; j < numberOfElements; ++j ) { // y-direction
        // get the coords of the point of interest
        mx = i * levelDist;
        my = j * levelDist;

        // now get the maximum of the five inner points of each
        // adjacent element and add the maximum of them to the
        // value of the vertex of interest
        // there are at most 4 adjacent elements, 1) left upper
        // 2) right upper, 3) left lower, 4) right lower
        max = 0;

        if ( this->get ( mx + levelDist / 2, my ) > max ) { max = this->get ( mx + levelDist / 2, my );}
        if ( this->get ( mx, my + levelDist / 2 ) > max ) { max = this->get ( mx, my + levelDist / 2 );}
        if ( this->get ( mx + levelDist / 2, my + levelDist / 2 ) > max ) { max = this->get ( mx + levelDist / 2, my + levelDist / 2 );}
        if ( this->get ( mx + levelDist, my + levelDist / 2 ) > max ) { max = this->get ( mx + levelDist, my + levelDist / 2 );}
        if ( this->get ( mx + levelDist / 2, my + levelDist ) > max ) { max = this->get ( mx + levelDist / 2, my + levelDist );}

        if ( this->get ( mx, my ) < max ) {this->set ( mx, my, max ); }
        if ( this->get ( mx + levelDist, my ) < max ) { this->set ( mx + levelDist, my, max ); }
        if ( this->get ( mx, my + levelDist ) < max ) { this->set ( mx, my + levelDist, max ); }
        if ( this->get ( mx + levelDist, my + levelDist ) < max ) { this->set ( mx + levelDist, my + levelDist, max ); }
      } // j
    } // i

    // in the next coarser level there are only the half number of elements,
    // but the steps are 2 times bigger
    numberOfElements /= 2;
    levelDist *= 2;
  } // level
}

template class EstimatorTP2d<float>;
template class EstimatorTP2d<double>;
template class EstimatorTP2d<long double>;
template class EstimatorTP2d<int>;
template class EstimatorTP2d<unsigned char>;










template <class T>
EstimatorTP3d<T>::EstimatorTP3d ( int Size )
  : qc::ScalarArray<T, qc::QC_3D> ( Size, Size, Size ) {

  maxLevel = qc::logBaseTwo ( Size );
}

template <class T>
EstimatorTP3d<T>::~EstimatorTP3d() {
}


template <class T>
void EstimatorTP3d<T>::makeSaturationDownNode() {
  int count = 0, total = 0, Level;
  for ( Level = maxLevel - 1; Level >= 0; Level-- ) {
    total |= ( 1 << Level );
  }

  for ( Level = maxLevel - 1; Level >= 0; Level-- ) {

    int fullStep = 1 << ( maxLevel - Level );
    int halfStep = fullStep >> 1;
    for ( int X = 0; X < this->numX - 1; X += fullStep ) {
      cout << ( count * 100 / total ) << "% done.     \r";
      count++;
      for ( int Y = 0; Y < this->numY - 1; Y += fullStep ) {
        for ( int Z = 0; Z < this->numZ - 1; Z += fullStep ) {
          T max = maxChildNodeValue ( X, Y, Z, halfStep );

          this->setMax ( X           , Y           , Z, max );
          this->setMax ( X           , Y + fullStep, Z, max );
          this->setMax ( X + fullStep, Y           , Z, max );
          this->setMax ( X + fullStep, Y + fullStep, Z, max );

          this->setMax ( X           , Y           , Z + fullStep, max );
          this->setMax ( X           , Y + fullStep, Z + fullStep, max );
          this->setMax ( X + fullStep, Y           , Z + fullStep, max );
          this->setMax ( X + fullStep, Y + fullStep, Z + fullStep, max );
        }
      }
    }
  }
  cout << endl;
}


template <class T>
void EstimatorTP3d<T>::makeSaturationUpNode() {
  int total = 0, count = 0, Level;
  for ( Level = maxLevel - 1; Level >= 0; Level-- ) {
    total |= ( 1 << Level );
  }
  for ( Level = maxLevel - 1; Level >= 0; Level-- ) {
    int fullStep = 1 << ( maxLevel - Level );
    int halfStep = fullStep >> 1;
    for ( int X = 0; X < this->numX - 1; X += fullStep ) {
      cout << ( count * 100 / total ) << "% done.   \r";
      count++;
      for ( int Y = 0; Y < this->numY - 1; Y += fullStep ) {
        for ( int Z = 0; Z < this->numZ - 1; Z += fullStep ) {
          T min = minChildNodeValue ( X, Y, Z, halfStep );

          this->setMin ( X           , Y           , Z, min );
          this->setMin ( X           , Y + fullStep, Z, min );
          this->setMin ( X + fullStep, Y           , Z, min );
          this->setMin ( X + fullStep, Y + fullStep, Z, min );

          this->setMin ( X           , Y           , Z + fullStep, min );
          this->setMin ( X           , Y + fullStep, Z + fullStep, min );
          this->setMin ( X + fullStep, Y           , Z + fullStep, min );
          this->setMin ( X + fullStep, Y + fullStep, Z + fullStep, min );

        }
      }
    }
  }
  cout << endl;
}


template <class T>
void EstimatorTP3d<T>::makeSaturationUpElement() {
  int total = 0, count = 0, cs, level;
  int FullStep, X, Y, Z, XShift, YShift, ZShift;

  for ( level = maxLevel - 1; level > 0; level-- ) {
    total += ( 1 << level ) * ( 1 << level ) * ( 1 << level );
  }

  for ( level = maxLevel - 1; level >= 0; level -- ) {
    int HalfShift = maxLevel - ( level + 1 );
    int FullShift = maxLevel - level;
    FullStep = 1 << ( maxLevel - level );

    cs = ( 1 << level ) * ( 1 << level );

    for ( X = 0; X < this->numX ; X += FullStep ) {
      cerr << ( count * 100 ) / ( total ) << "% done.   \r";
      count += cs;

      XShift = X >> FullShift;
      for ( Y = 0; Y < this->numY ; Y += FullStep ) {
        YShift = Y >> FullShift;
        for ( Z = 0; Z < this->numZ ; Z += FullStep ) {
          ZShift = Z >> FullShift;
          T min = elSatArrays[ level + 1 ]->getElementSaturationMin ( X >> HalfShift,
                                                                      Y >> HalfShift,
                                                                      Z >> HalfShift );
          elSatArrays[ level ]->set ( XShift, YShift, ZShift, min );
        }
      }
    }
  }
  cout << endl;
}

template <class T>
void EstimatorTP3d<T>::makeSaturationDownElement() {
  int total = 0, count = 0, cs, level;
  int FullStep, X, Y, Z, XShift, YShift, ZShift;

  for ( level = maxLevel - 1; level > 0; level-- ) {
    total += ( 1 << level ) * ( 1 << level ) * ( 1 << level );
  }

  for ( level = maxLevel - 1; level >= 0; level -- ) {
    int HalfShift = maxLevel - ( level + 1 );
    int FullShift = maxLevel - level;
    FullStep = 1 << ( maxLevel - level );

    cs = ( 1 << level ) * ( 1 << level );

    for ( X = 0; X < this->numX ; X += FullStep ) {
      cerr << ( count * 100 ) / ( total ) << "% done.   \r";
      count += cs;

      XShift = X >> FullShift;
      for ( Y = 0; Y < this->numY ; Y += FullStep ) {
        YShift = Y >> FullShift;
        for ( Z = 0; Z < this->numZ ; Z += FullStep ) {
          ZShift = Z >> FullShift;
          T max = elSatArrays[ level + 1 ]->getElementSaturationMax ( X >> HalfShift,
                                                                      Y >> HalfShift,
                                                                      Z >> HalfShift );
          elSatArrays[ level ]->set ( XShift, YShift, ZShift, max );
        }
      }
    }
  }
  cout << endl;
}



/** *******************************************************
 * Method generates a saturated Error Array, by first
 * calculating the error and second saturate it.
 * ********************************************************/
//! method walks across the data-array and fills the own array
//! with the saturated error-entries
template <class T>
void EstimatorTP3d<T>::makeSaturatedErrorArray ( qc::ScalarArray<T, qc::QC_3D>& meshData ) {
  int sizeOfMeshData = meshData.getNumX();
  int maxLevelOfMeshData = qc::logBaseTwo ( sizeOfMeshData );

  int level, levelDist;
  int i, j, k;
  int mx, my, mz, lx, ly, lz;
  T fx, fy, fz, kompl;
  T max, a, b, c, d, error;

  // the error on the finest level is zero
  for ( i = 0; i < sizeOfMeshData; ++i )
    for ( j = 0; j < sizeOfMeshData; ++j )
      for ( k = 0; k < sizeOfMeshData; ++k )
        this->set ( i, j, k, 0 );

  // calc the number of qc::Elements in one line in the max. level
  int numberOfElements = 1 << ( maxLevelOfMeshData - 1 );
  levelDist = 1;

  // -------------------- 1. Error-calculation --------------
  // the first thing to do is calculating the error between a
  // vertex on level n and on the finer level n+1.
  // Now walk through the levels, maxLevelOfMeshData is the finest level
  for ( level = maxLevelOfMeshData; level > 0; --level ) {
    // walk through the vertices
    for ( i = 0; i < numberOfElements; ++i ) { // x-direction
      for ( j = 0; j < numberOfElements; ++j ) { // y-direction
        for ( k = 0; k < numberOfElements; ++k ) { // z-direction
          // calc the position of the middle vertex of the el. on the finer level
          mx = ( i * 2 + 1 ) * levelDist;
          my = ( j * 2 + 1 ) * levelDist;
          mz = ( k * 2 + 1 ) * levelDist;

          // now calc the max distance from the coarser approx. by
          // bilinear interpolation
          for ( lx = -1; lx <= 1; ++lx ) {
            for ( ly = -1; ly <= 1; ++ly ) {
              for ( lz = -1; lz <= 1; ++lz ) {
                // no vertices of the coarser cube
                if ( ( aol::Abs ( lx ) != 1 ) || ( aol::Abs ( ly ) != 1 ) || ( aol::Abs ( lz ) ) != 1 ) {
                  // calc the local coordinates of the actual point
                  fx = static_cast<T> ( lx + 1 ) / 2;
                  fy = static_cast<T> ( ly + 1 ) / 2;
                  fz = static_cast<T> ( lz + 1 ) / 2;
                  kompl = -fx + 1; // 1. - fx

                  // now trilinear interpolation of the data at this point
                  a = fx * meshData.get ( mx + levelDist, my + levelDist, mz + levelDist )
                      + kompl * meshData.get ( mx - levelDist, my + levelDist, mz + levelDist );
                  b = fx * meshData.get ( mx + levelDist, my - levelDist, mz + levelDist )
                      + kompl * meshData.get ( mx - levelDist, my - levelDist, mz + levelDist );
                  c = fx * meshData.get ( mx + levelDist, my + levelDist, mz - levelDist )
                      + kompl * meshData.get ( mx - levelDist, my + levelDist, mz - levelDist );
                  d = fx * meshData.get ( mx + levelDist, my - levelDist, mz - levelDist )
                      + kompl * meshData.get ( mx - levelDist, my - levelDist, mz - levelDist );

                  kompl = -fy + 1; //1.-fy;
                  a = fy * a + kompl * b;
                  b = fy * c + kompl * d;

                  a = fz * a + ( -fz + 1 ) * b; // (1.-fz)*b;

                  error = aol::Abs ( meshData.get ( mx + lx * levelDist, my + ly * levelDist, mz + lz * levelDist ) - a );
                  this->set ( mx + lx * levelDist, my + ly * levelDist, mz + lz * levelDist, error );
                } // if
              } // lz
            } // ly
          } // lx

        } // k (z)
      } // j (y)
    } // i (x)
    // in the next coarser level there are only the half number of elements,
    // but the steps are 2 times bigger
    numberOfElements /= 2;
    levelDist *= 2;
  } // level


  // -------------------- 2. Saturation ---------------------
  // now saturate the grid by adding the maximum of the five
  // circumjacent vertices of the finer level.

  numberOfElements = 1 << ( maxLevelOfMeshData - 1 );
  levelDist = 2;

  // again walk through the levels, maxLevelOfMeshData is the finest level
  for ( level = maxLevelOfMeshData - 1; level >= 0; --level ) {
    // walk through the vertices
    for ( i = 0; i <= numberOfElements; ++i ) { // x-direction
      for ( j = 0; j <= numberOfElements; ++j ) { // y-direction
        for ( k = 0; k <= numberOfElements; ++k ) { // z-direction
          // get the coords of the point of interest
          mx = i * levelDist;
          my = j * levelDist;
          mz = k * levelDist;

          // now get the maximum of the five inner points of each
          // adjacent element and add the maximum of them to the
          // value of the vertex of interest
          // there are at most 4 adjacent elements, 1) left upper
          // 2) right upper, 3) left lower, 4) right lower
          max = 0;
          a = 0;
          if ( ( mx > 0 ) && ( my > 0 ) && ( mz > 0 ) )
            a = maxChildNodeValue ( mx - levelDist, my - levelDist, mz - levelDist, levelDist / 2 );
          if ( max < a ) max = a;
          if ( ( mx < sizeOfMeshData - 1 ) && ( my > 0 ) && ( mz > 0 ) )
            a = maxChildNodeValue ( mx, my - levelDist, mz - levelDist, levelDist / 2 );
          if ( max < a ) max = a;
          if ( ( mx > 0 ) && ( my < sizeOfMeshData - 1 ) && ( mz > 0 ) )
            a = maxChildNodeValue ( mx - levelDist, my, mz - levelDist, levelDist / 2 );
          if ( max < a ) max = a;
          if ( ( mx < sizeOfMeshData - 1 ) && ( my < sizeOfMeshData - 1 ) && ( mz > 0 ) )
            a = maxChildNodeValue ( mx, my, mz - levelDist, levelDist / 2 );
          if ( max < a ) max = a;

          if ( ( mx > 0 ) && ( my > 0 ) && ( mz < sizeOfMeshData - 1 ) )
            a = maxChildNodeValue ( mx - levelDist, my - levelDist, mz, levelDist / 2 );
          if ( max < a ) max = a;
          if ( ( mx < sizeOfMeshData - 1 ) && ( my > 0 ) && ( mz < sizeOfMeshData - 1 ) )
            a = maxChildNodeValue ( mx, my - levelDist, mz, levelDist / 2 );
          if ( max < a ) max = a;
          if ( ( mx > 0 ) && ( my < sizeOfMeshData - 1 ) && ( mz < sizeOfMeshData - 1 ) )
            a = maxChildNodeValue ( mx - levelDist, my, mz, levelDist / 2 );
          if ( max < a ) max = a;
          if ( ( mx < sizeOfMeshData - 1 ) && ( my < sizeOfMeshData - 1 ) && ( mz < sizeOfMeshData - 1 ) )
            a = maxChildNodeValue ( mx, my, mz, levelDist / 2 );
          if ( max < a ) max = a;

          // now add this maximum to the actual point
          a = this->get ( mx, my, mz );
          this->set ( mx, my, mz, max + a );

        } // k
      } // j
    } // i

    // in the next coarser level there are only the half number of elements,
    // but the steps are 2 times bigger
    numberOfElements /= 2;
    levelDist *= 2;
  } // level

}

template class EstimatorTP3d<unsigned char>;
template class EstimatorTP3d<float>;
template class EstimatorTP3d<double>;
