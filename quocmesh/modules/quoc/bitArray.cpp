#include <bitArray.h>
#include <scalarArray.h>

void qc::BitArray<qc::QC_2D>::save ( const char *FileName ) const {
  qc::ScalarArray<unsigned char, qc::QC_2D> saveArray ( getNumX(), getNumY() );
  for ( int i = 0; i < this->size(); i++ ) {
    if ( this->get ( i ) )
      saveArray[i] = 255;
  }
  saveArray.save ( FileName, qc::PGM_UNSIGNED_CHAR_ASCII );
}

void qc::BitArray<qc::QC_2D>::dilateByOne ( ) {
  qc::ScalarArray<unsigned char, qc::QC_2D> temp ( getNumX(), getNumY() );
  temp.assignFrom ( *this );

  // Inspired by http://ostermiller.org/dilate_and_erode.html.
  for ( int i = 0; i < getNumX(); ++i ) {
    for ( int j = 0; j < getNumY(); ++j ) {
      if ( temp.get ( i, j ) == 1 ) {
        if ( i > 0 && temp.get ( i-1, j ) == 0 )
          temp.set ( i-1, j, 2 );
        if ( j > 0 && temp.get ( i, j-1 ) == 0 )
          temp.set ( i, j-1, 2 );
        if ( i+1 < getNumX() && temp.get ( i+1, j ) == 0 )
          temp.set ( i+1, j, 2 );
        if ( j+1 < getNumY() && temp.get ( i, j+1 ) == 0 )
          temp.set ( i, j+1, 2 );
      }
    }
  }

  for ( int i = 0; i < getNumX(); ++i ) {
    for ( int j = 0; j < getNumY(); ++j ) {
      if ( temp.get ( i, j ) == 2 ) {
        set ( i, j, true );
      }
    }
  }
}


void qc::BitArray<qc::QC_2D>::erodeByOne ( ) {
  qc::ScalarArray<unsigned char, qc::QC_2D> temp ( getNumX(), getNumY() );
  temp.assignFrom ( *this );
  
  // Inspired by http://ostermiller.org/dilate_and_erode.html.
  for ( int i = 0; i < getNumX(); ++i ) {
    for ( int j = 0; j < getNumY(); ++j ) {
      if ( temp.get ( i, j ) == 0 ) {
        if ( i > 0 && temp.get ( i-1, j ) == 1 )
          temp.set ( i-1, j, 2 );
        if ( j > 0 && temp.get ( i, j-1 ) == 1 )
          temp.set ( i, j-1, 2 );
        if ( i+1 < getNumX() && temp.get ( i+1, j ) == 1 )
          temp.set ( i+1, j, 2 );
        if ( j+1 < getNumY() && temp.get ( i, j+1 ) == 1 )
          temp.set ( i, j+1, 2 );
      }
    }
  }
  
  for ( int i = 0; i < getNumX(); ++i ) {
    for ( int j = 0; j < getNumY(); ++j ) {
      if ( temp.get ( i, j ) == 2 ) {
        set ( i, j, false );
      }
    }
  }
}



void qc::BitArray<qc::QC_2D>::saveToFile ( const char *filename ) const {
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( filename );
  throw aol::UnimplementedCodeException ( "qc::BitArray<QC_2D>::saveToFile not implemented yet", __FILE__, __LINE__ );
}


void qc::BitArray<qc::QC_2D>::loadFromFile ( const char *filename ) {
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( filename );
  throw aol::UnimplementedCodeException ( "qc::BitArray<QC_2D>::loadFromFile not implemented yet", __FILE__, __LINE__ );
}


void qc::BitArray<qc::QC_3D>::save ( const char *FileName ) const {
  qc::ScalarArray<unsigned char, qc::QC_3D> saveArray ( getNumX(), getNumY(), getNumZ() );
  for ( int i = 0; i < this->size(); i++ ) {
    if ( this->get ( i ) )
      saveArray[i] = 255;
  }
  saveArray.save ( FileName, qc::PGM_UNSIGNED_CHAR_ASCII );
}

void qc::BitArray<qc::QC_3D>::prolongateTo ( qc::BitArray<qc::QC_3D> & finerArray ) const {

  QUOC_ASSERT ( finerArray.getNumX() == 2 * this->getNumX() - 1 );
  QUOC_ASSERT ( finerArray.getNumY() == 2 * this->getNumY() - 1 );
  QUOC_ASSERT ( finerArray.getNumZ() == 2 * this->getNumZ() - 1 );

  int n_x = this->getNumX(),
      n_y = this->getNumY(),
      n_z = this->getNumZ();

  // compute all other entries
  for (int i = 0; i < n_x; ++i)
    for (int j = 0; j < n_y; ++j)
      for (int k = 0; k < n_z; ++k) {

        // point that exists in both arrays
        finerArray.set ( 2 * i, 2 * j, 2 * k, this->get(i, j, k) );

        // midpoint of lines
        if ( i + 1 < n_x )
          finerArray.set ( 2 * i + 1, 2 * j,     2 * k,     this->get(i, j, k) && this->get(i+1, j,   k)   );
        if ( j + 1 < n_y )
          finerArray.set ( 2 * i,     2 * j + 1, 2 * k,     this->get(i, j, k) && this->get(i,   j+1, k)   );
        if ( k + 1 < n_z )
          finerArray.set ( 2 * i,     2 * j,     2 * k + 1, this->get(i, j, k) && this->get(i,   j,   k+1) );

        // midpoint of surfaces
        if ( i + 1 < n_x && j + 1 < n_y )
          finerArray.set ( 2 * i + 1, 2 * j + 1, 2 * k,        this->get(i,   j,   k) && this->get(i,   j+1, k  )
                                                            && this->get(i+1, j,   k) && this->get(i+1, j+1, k  ) );
        if (i + 1 < n_x && k + 1 < n_z )
          finerArray.set ( 2 * i + 1, 2 * j,     2 * k + 1,    this->get(i,   j,   k) && this->get(i,   j,   k+1)
                                                            && this->get(i+1, j,   k) && this->get(i+1, j,   k+1) );
        if ( j + 1 < n_y && k + 1 < n_z )
          finerArray.set ( 2 * i,     2 * j + 1, 2 * k + 1,    this->get(i,   j,   k) && this->get(i,   j,   k+1)
                                                            && this->get(i,   j+1, k) && this->get(i,   j+1, k+1) );

        // midpoint of element
        if (i + 1 < n_x && j + 1 < n_y && k + 1 < n_z )
          finerArray.set ( 2 * i + 1, 2 * j + 1, 2 * k + 1,    this->get(i,   j,   k  ) && this->get(i+1, j,   k  )
                                                            && this->get(i,   j+1, k  ) && this->get(i,   j,   k+1)
                                                            && this->get(i+1, j+1, k  ) && this->get(i+1, j,   k+1)
                                                            && this->get(i+1, j,   k+1) && this->get(i+1, j+1, k+1) );
      }
}

void qc::BitArray<qc::QC_3D>::restrictTo ( qc::BitArray<qc::QC_3D> & coarserArray ) const {

  QUOC_ASSERT ( 2 * coarserArray.getNumX() - 1 == this->getNumX() );
  QUOC_ASSERT ( 2 * coarserArray.getNumY() - 1 == this->getNumY() );
  QUOC_ASSERT ( 2 * coarserArray.getNumZ() - 1 == this->getNumZ() );

  for (int i = 0; i < coarserArray.getNumX(); ++i)
    for (int j = 0; j < coarserArray.getNumY(); ++j)
      for (int k = 0; k < coarserArray.getNumZ(); ++k)
        coarserArray.set ( i, j, k, this->get(2 * i, 2 * j, 2 * k) );
}



void qc::BitArray<qc::QC_3D>::saveToFile ( const char *filename ) const {
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( filename );
  throw aol::UnimplementedCodeException ( "qc::BitArray<QC_3D>::saveToFile not implemented yet", __FILE__, __LINE__ );
}


void qc::BitArray<qc::QC_3D>::loadFromFile ( const char *filename ) {
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( filename );
  throw aol::UnimplementedCodeException ( "qc::BitArray<QC_3D>::loadFromFile not implemented yet", __FILE__, __LINE__ );
}
