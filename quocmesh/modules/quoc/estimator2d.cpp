#include <quoc.h>
#include <estimator2d.h>


template <class T>
qc::Estimator2d<T>::Estimator2d ( int Size )
    : qc::ScalarArray<T, qc::QC_2D> ( Size, Size ),
    sat_type ( EST_SAT_NONE ) {

  maxLevel = qc::logBaseTwo ( Size );

  if ( Size & 1 ) {
    sat_type |= EST_SAT_NODE;
    elSatArrays = NULL;
  } else {
    // ELEMENT saturation needs more arrays for storing values
    typedef qc::ScalarArray<T, qc::QC_2D>* ScalarArrayPtr;
    elSatArrays = new ScalarArrayPtr[ maxLevel + 1 ];
    for ( int l = 0; l < maxLevel; l++ ) {
      elSatArrays[ l ] = new qc::ScalarArray<T, qc::QC_2D> ( Size >> ( maxLevel - l ),
                                                    Size >> ( maxLevel - l ) );
    }
    elSatArrays[ maxLevel ] = this;
  }
}

template <class T>
qc::Estimator2d<T>::~Estimator2d() {
  if ( elSatArrays ) {
    for ( int l = 0; l < maxLevel; l++ ) {
      delete elSatArrays[ l ];
    }
    delete[] elSatArrays;
  }
}

template <class T>
void qc::Estimator2d<T>::makeSaturationDown() {
  if ( sat_type & EST_SAT_NODE ) {
    makeSaturationDownNode();
  } else {
    makeSaturationDownElement();
  }
}

template <class T>
void qc::Estimator2d<T>::makeSaturationUp() {
  if ( sat_type & EST_SAT_NODE ) {
    makeSaturationUpNode();
  } else {
    makeSaturationUpElement();
  }
}

template <class T>
void qc::Estimator2d<T>::makeSaturationDownNode() {
  for ( int Level = maxLevel - 1; Level >= 0; Level-- ) {

    int fullStep = 1 << ( maxLevel - Level );
    int halfStep = fullStep >> 1;

    for ( int X = 0; X < this->numX - 1; X += fullStep ) {
      for ( int Y = 0; Y < this->numY - 1; Y += fullStep ) {
        T max = maxChildNodeValue ( X, Y, halfStep );

        this->setMax ( X           , Y           , max );
        this->setMax ( X           , Y + fullStep, max );
        this->setMax ( X + fullStep, Y           , max );
        this->setMax ( X + fullStep, Y + fullStep, max );

      }
    }
  }
  sat_type |= EST_SAT_DOWN;
  cout << endl;
}


template <class T>
void qc::Estimator2d<T>::makeSaturationUpNode() {
  for ( int Level = maxLevel - 1; Level >= 0; Level-- ) {
    int fullStep = 1 << ( maxLevel - Level );
    int halfStep = fullStep >> 1;
    for ( int X = 0; X < this->numX - 1; X += fullStep ) {
      for ( int Y = 0; Y < this->numY - 1; Y += fullStep ) {
        T min = minChildNodeValue ( X, Y, halfStep );

        this->setMin ( X           , Y           , min );
        this->setMin ( X           , Y + fullStep, min );
        this->setMin ( X + fullStep, Y           , min );
        this->setMin ( X + fullStep, Y + fullStep, min );
      }
    }
  }
  sat_type |= EST_SAT_UP;
}


template <class T>
void qc::Estimator2d<T>::makeSaturationUpElement() {

  cerr << "makeSaturationUpElement( ) \n";

  int FullStep, X, Y, XShift, YShift;
  for ( int level = maxLevel - 1; level >= 0; level -- ) {
    int HalfShift = maxLevel - ( level + 1 );
    int FullShift = maxLevel - level;
    FullStep = 1 << ( maxLevel - level );

    for ( X = 0; X < this->numX ; X += FullStep ) {
      XShift = X >> FullShift;
      for ( Y = 0; Y < this->numY ; Y += FullStep ) {
        YShift = Y >> FullShift;
        T min = elSatArrays[ level + 1 ]->getElementSaturationMin ( X >> HalfShift,
                                                                    Y >> HalfShift );
        elSatArrays[ level ]->set ( XShift, YShift, min );
      }
    }
  }
  sat_type |= EST_SAT_UP;
}

template <class T>
void qc::Estimator2d<T>::makeSaturationDownElement() {


  cerr << "makeSaturationDownElement( ) \n";

  int FullStep, X, Y, XShift, YShift;
  for ( int level = maxLevel - 1; level >= 0; level -- ) {
    int HalfShift = maxLevel - ( level + 1 );
    int FullShift = maxLevel - level;
    FullStep = 1 << ( maxLevel - level );

    for ( X = 0; X < this->numX ; X += FullStep ) {
      XShift = X >> FullShift;
      for ( Y = 0; Y < this->numY ; Y += FullStep ) {
        YShift = Y >> FullShift;
        T max = elSatArrays[ level + 1 ]->getElementSaturationMax ( X >> HalfShift,
                                                                    Y >> HalfShift );
        elSatArrays[ level ]->set ( XShift, YShift, max );
      }
    }
  }
  sat_type |= EST_SAT_DOWN;
}


/** *******************************************************
 * Method generates a saturated Error Array, by first
 * calculating the error and second saturate it.
 * ********************************************************/
//! method walks across the data-array and fills the own array
//! with the saturated error-entries
template <class T>
void qc::Estimator2d<T>::makeSaturatedErrorArray ( qc::ScalarArray<T,qc::QC_2D>* meshData ) {
  int sizeOfMeshData = meshData->getNumX();
  int maxLevelOfMeshData = qc::logBaseTwo ( sizeOfMeshData );

  int level, levelDist;
  int i, j, k;
  int mx, my;
  T error;
  T max, a, b;

  // the error on the finest level is zero
  for ( i = 0; i < sizeOfMeshData; ++i )
    for ( j = 0; j < sizeOfMeshData; ++j )
      this->set ( i, j, 0 );

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
        // calc the position of the middle vertex of the el. on the finer level
        mx = ( i * 2 + 1 ) * levelDist;
        my = ( j * 2 + 1 ) * levelDist;
        // now calc the max distance from the coarser approx.
        // to the finest, first the horizontal lines
        for ( k = -1; k <= 1; k += 2 ) {
          a = meshData->get ( mx + levelDist, my + k * levelDist );
          b = meshData->get ( mx - levelDist, my + k * levelDist );

          error = aol::Abs ( meshData->get ( mx, my + k * levelDist ) - ( a + b ) / static_cast<T> ( 2 ) );
          this->set ( mx, my + k*levelDist, error );
        }
        // now the vertical lines
        for ( k = -1; k <= 1; k += 2 ) {
          a = meshData->get ( mx + k * levelDist, my + levelDist );
          b = meshData->get ( mx + k * levelDist, my - levelDist );

          error = aol::Abs ( meshData->get ( mx + k * levelDist, my ) - ( a + b ) / 2 );
          this->set ( mx + k*levelDist, my, error );
        }
        // still missing: the middle element
        // here the error is calculated by bilinear interpolation

        a = ( meshData->get ( mx - levelDist, my - levelDist ) + meshData->get ( mx + levelDist, my - levelDist ) ) / 2;
        b = ( meshData->get ( mx - levelDist, my + levelDist ) + meshData->get ( mx + levelDist, my + levelDist ) ) / 2;

        error = aol::Abs ( meshData->get ( mx, my ) - ( a + b ) / 2 );
        this->set ( mx, my, error );
      } // j
    } // i
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
        // get the coords of the point of interest
        mx = i * levelDist;
        my = j * levelDist;

        // now get the maximum of the five inner points of each
        // adjacent element and add the maximum of them to the
        // value of the vertex of interest
        // there are at most 4 adjacent elements, 1) left upper
        // 2) right upper, 3) left lower, 4) right lower
        max = 0;
        a = 0;
        if ( ( mx > 0 ) && ( my > 0 ) ) a = maxChildNodeValue ( mx - levelDist, my - levelDist, levelDist / 2 );
        if ( max < a ) max = a;
        if ( ( mx < sizeOfMeshData - 1 ) && ( my > 0 ) ) a = maxChildNodeValue ( mx, my - levelDist, levelDist / 2 );
        if ( max < a ) max = a;
        if ( ( mx > 0 ) && ( my < sizeOfMeshData - 1 ) ) a = maxChildNodeValue ( mx - levelDist, my, levelDist / 2 );
        if ( max < a ) max = a;
        if ( ( mx < sizeOfMeshData - 1 ) && ( my < sizeOfMeshData - 1 ) ) a = maxChildNodeValue ( mx, my, levelDist / 2 );
        if ( max < a ) max = a;

        // now add this maximum to the actual point
        a = this->get ( mx, my );
        this->set ( mx, my, max + a );

      } // j
    } // i

    // in the next coarser level there are only the half number of elements,
    // but the steps are 2 times bigger
    numberOfElements /= 2;
    levelDist *= 2;
  } // level
}

template class qc::Estimator2d<float>;
template class qc::Estimator2d<double>;
template class qc::Estimator2d<long double>;
template class qc::Estimator2d<int>;
template class qc::Estimator2d<unsigned char>;

