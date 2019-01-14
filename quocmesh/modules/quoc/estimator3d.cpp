#include <estimator3d.h>

template <class T>
qc::Estimator3d<T>::Estimator3d ( int Size )
    : qc::ScalarArray<T, qc::QC_3D> ( Size, Size, Size ),
    sat_type ( EST_SAT_NONE ) {

  maxLevel = qc::logBaseTwo ( Size );

  if ( Size & 1 ) {
    sat_type |= EST_SAT_NODE;
    elSatArrays = NULL;
  } else {
    // ELEMENT saturation needs more arrays for storing values
    typedef qc::ScalarArray<T, qc::QC_3D>* ScalarArrayPtr;
    elSatArrays = new ScalarArrayPtr[ maxLevel + 1 ];
    for ( int l = 0; l < maxLevel; l++ ) {
      elSatArrays[ l ] = new qc::ScalarArray<T, qc::QC_3D> ( Size >> ( maxLevel - l ),
                                                    Size >> ( maxLevel - l ),
                                                    Size >> ( maxLevel - l ) );
    }
    elSatArrays[ maxLevel ] = this;
  }
}

template <class T>
qc::Estimator3d<T>::~Estimator3d() {
  if ( elSatArrays ) {
    for ( int l = 0; l < maxLevel; l++ ) {
      delete elSatArrays[ l ];
    }
    delete[] elSatArrays;
  }
}

template <class T>
void qc::Estimator3d<T>::makeSaturationDown() {
  if ( sat_type & EST_SAT_NODE ) {
    makeSaturationDownNode();
  } else {
    makeSaturationDownElement();
  }
}

template <class T>
void qc::Estimator3d<T>::makeSaturationUp() {
  if ( sat_type & EST_SAT_NODE ) {
    makeSaturationUpNode();
  } else {
    makeSaturationUpElement();
  }
}


template <class T>
void qc::Estimator3d<T>::makeSaturationDownNode() {
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
  sat_type = EST_SAT_DOWN | EST_SAT_NODE;
  cout << endl;
}


template <class T>
void qc::Estimator3d<T>::makeSaturationUpNode() {
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
  sat_type = EST_SAT_UP | EST_SAT_NODE;
}


template <class T>
void qc::Estimator3d<T>::makeSaturationUpElement() {
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
  sat_type = EST_SAT_UP;
}

template <class T>
void qc::Estimator3d<T>::makeSaturationDownElement() {
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
  sat_type = EST_SAT_DOWN;
}



/** *******************************************************
 * Method generates a saturated Error Array, by first
 * calculating the error and second saturate it.
 * ********************************************************/
//! method walks across the data-array and fills the own array
//! with the saturated error-entries
template <class T>
void qc::Estimator3d<T>::makeSaturatedErrorArray ( qc::ScalarArray<T, qc::QC_3D>* meshData ) {
  int sizeOfMeshData = meshData->getNumX();
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
    for ( i = 0; i < numberOfElements; ++i ) // x-direction
    {
      for ( j = 0; j < numberOfElements; ++j ) // y-direction
      {
        for ( k = 0; k < numberOfElements; ++k ) // z-direction
        {
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
                  a = fx * meshData->get ( mx + levelDist, my + levelDist, mz + levelDist )
                      + kompl * meshData->get ( mx - levelDist, my + levelDist, mz + levelDist );
                  b = fx * meshData->get ( mx + levelDist, my - levelDist, mz + levelDist )
                      + kompl * meshData->get ( mx - levelDist, my - levelDist, mz + levelDist );
                  c = fx * meshData->get ( mx + levelDist, my + levelDist, mz - levelDist )
                      + kompl * meshData->get ( mx - levelDist, my + levelDist, mz - levelDist );
                  d = fx * meshData->get ( mx + levelDist, my - levelDist, mz - levelDist )
                      + kompl * meshData->get ( mx - levelDist, my - levelDist, mz - levelDist );

                  kompl = -fy + 1; //1.-fy;
                  a = fy * a + kompl * b;
                  b = fy * c + kompl * d;

                  a = fz * a + ( -fz + 1 ) * b; // (1.-fz)*b;

                  error = aol::Abs ( meshData->get ( mx + lx * levelDist, my + ly * levelDist, mz + lz * levelDist ) - a );
                  this->set ( mx + lx*levelDist, my + ly*levelDist, mz + lz*levelDist, error );
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
    for ( i = 0; i <= numberOfElements; ++i ) // x-direction
    {
      for ( j = 0; j <= numberOfElements; ++j ) // y-direction
      {
        for ( k = 0; k <= numberOfElements; ++k ) // z-direction
        {
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

template class qc::Estimator3d<unsigned char>;
template class qc::Estimator3d<float>;
template class qc::Estimator3d<double>;
