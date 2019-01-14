#include <shapeLevelsetGenerator.h>

#include <iterators.h>
#include <randomGenerator.h>

namespace qc {

template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateColumnLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const DataType radius ) {
  generateEllipticColumnLevelset ( levelset, radius, radius );
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateEllipticColumnLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const DataType x_rad, const DataType y_rad ) {
  const int width = levelset.getNumXYZ();

  const DataType ctr_x = 0.5, ctr_y = 0.5;

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const int i = ( *bit ) [0], j = ( *bit ) [1];
    const DataType  x = ( 1.0 * i ) / ( width - 1.0 ),  y = ( 1.0 * j ) / ( width - 1.0 );
    const DataType value = aol::Sqr ( x - ctr_x ) / aol::Sqr ( x_rad ) + aol::Sqr ( y - ctr_y ) / aol::Sqr ( y_rad ) - 1;

    levelset.set ( *bit, value );
  }
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateConeLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const DataType low_rad, const DataType high_rad ) {
  const int width = levelset.getNumXYZ();

  const DataType ctr_x = 0.5, ctr_y = 0.5;

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const int i = ( *bit ) [0], j = ( *bit ) [1], k = ( *bit ) [2];
    const DataType rad = low_rad + ( high_rad - low_rad ) * ( static_cast<DataType> ( k ) / static_cast<DataType> ( width - 1 ) );
    const DataType x = ( 1.0 * i ) / ( width - 1.0 ),   y = ( 1.0 * j ) / ( width - 1.0 );
    const DataType value = aol::Sqr ( x - ctr_x ) / aol::Sqr ( rad ) + aol::Sqr ( y - ctr_y ) / aol::Sqr ( rad ) - 1;

    levelset.set ( *bit, value );
  }

}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateBallLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const DataType radius ) {
  /*
  // something is wrong if we do the following (OK, the values are different, but the reconstructed level set should topologically be the same - still the selfTest fails)
  ShapeLevelsetGenerator<DataType>::generateEllipsoidLevelset ( levelset, aol::Vec3<DataType> ( radii, radii, radii), aol::Vec3<DataType> ( 0.5, 0.5, 0.5 ) );
  */

  const DataType h = 1.0 / aol::Max ( levelset.getNumX() - 1, levelset.getNumY() - 1 , levelset.getNumZ() - 1 );
  const DataType ctr_x = 0.5, ctr_y = 0.5, ctr_z = 0.5;

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const int i = ( *bit ) [0], j = ( *bit ) [1], k = ( *bit ) [2];
    const DataType x = h * i, y = h * j, z = h * k;
    const DataType value = aol::Sqr ( x - ctr_x ) + aol::Sqr ( y - ctr_y ) + aol::Sqr ( z - ctr_z ) - aol::Sqr ( radius );
    levelset.set ( *bit, value );
  }
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateEllipsoidLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const aol::Vec3<DataType> radii, const aol::Vec3<DataType> center ) {
  const DataType h = 1.0 / aol::Max ( levelset.getNumX() - 1, levelset.getNumY() - 1 , levelset.getNumZ() - 1 );

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const int i = ( *bit ) [0], j = ( *bit ) [1], k = ( *bit ) [2];
    const DataType x = h * i, y = h * j, z = h * k;
    const DataType value = aol::Sqr ( x - center[0] ) / aol::Sqr ( radii[0] ) + aol::Sqr ( y - center[1] ) / aol::Sqr ( radii[1] )  + aol::Sqr ( z - center[2] ) / aol::Sqr ( radii[2] ) - 1;
    levelset.set ( *bit, value );
  }
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generatePlateRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods, const DataType rad ) {
  const int width = levelset.getNumXYZ();
  aol::ProgressBar<> pb ( "Creating Plate/Rods levelset" );
  pb.start ( aol::Cub ( width ) );

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const int i = ( *bit ) [0], j = ( *bit ) [1], k = ( *bit ) [2];
    pb++;

    const DataType x = ( 1.0 * i ) / ( width - 1.0 ), y = ( 1.0 * j ) / ( width - 1.0 ), z = ( 1.0 * k ) / ( width - 1.0 );

    DataType value = 1e20;

    if ( ( aol::Abs ( z ) - 0.1 ) < value ) {                             // bottom plate
      value = ( aol::Abs ( z ) - 0.1 );
    }

    if ( ( aol::Abs ( 1.0 - z ) - 0.1 ) < value ) {                       // top plate
      value = ( aol::Abs ( 1.0 - z ) - 0.1 );
    }

    for ( int ir = 0; ir < n_rods; ++ir ) {
      for ( int jr = 0; jr < n_rods; ++jr ) {
        const DataType A = 0.1 + 0.8 * ( 2.0 * ir + 1.0 ) / ( 2.0 * n_rods );
        const DataType B = 0.1 + 0.8 * ( 2.0 * jr + 1.0 ) / ( 2.0 * n_rods );
        const DataType val = sqrt ( ( x - A ) * ( x - A ) + ( y - B ) * ( y - B ) ) - rad ;

        if ( val < value )
          value = val;
      }
    }

    levelset.set ( *bit, value );

  }

  pb.finish();
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generatePlateRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods ) {
  generatePlateRodsLevelset ( levelset, n_rods, 0.12 / ( 1.0 * n_rods ) );
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generate3DPlaterodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods, const DataType rad ) {
  const int width = levelset.getNumXYZ();
  aol::ProgressBar<> pb ( "Creating 3DPlateRods levelset" );
  pb.start ( aol::Cub ( width ) );

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const int i = ( *bit ) [0], j = ( *bit ) [1], k = ( *bit ) [2];
    pb++;

    const DataType x = ( 1.0 * i ) / ( width - 1.0 ), y = ( 1.0 * j ) / ( width - 1.0 ), z = ( 1.0 * k ) / ( width - 1.0 );

    DataType value = 1e20;

    const DataType val_bottom = ( aol::Abs ( z ) - 0.1 );
    const DataType val_top    = ( aol::Abs ( 1.0 - z ) - 0.1 );
    const DataType val_left   = ( aol::Abs ( x ) - 0.1 );
    const DataType val_right  = ( aol::Abs ( 1.0 - x ) - 0.1 );
    const DataType val_front  = ( aol::Abs ( y ) - 0.1 );
    const DataType val_back   = ( aol::Abs ( 1.0 - y ) - 0.1 );

    if ( val_bottom < value )   value = val_bottom;
    if ( val_top    < value )   value = val_top   ;
    if ( val_left   < value )   value = val_left  ;
    if ( val_right  < value )   value = val_right ;
    if ( val_front  < value )   value = val_front ;
    if ( val_back   < value )   value = val_back  ;

    for ( int ir = 0; ir < n_rods; ++ir ) {
      for ( int jr = 0; jr < n_rods; ++jr ) {
        const DataType A = 0.1 + 0.8 * ( 2.0 * ir + 1.0 ) / ( 2.0 * n_rods );
        const DataType B = 0.1 + 0.8 * ( 2.0 * jr + 1.0 ) / ( 2.0 * n_rods );
        const DataType val_x = sqrt ( ( y - A ) * ( y - A ) + ( z - B ) * ( z - B ) ) - rad ;
        const DataType val_y = sqrt ( ( x - A ) * ( x - A ) + ( z - B ) * ( z - B ) ) - rad ;
        const DataType val_z = sqrt ( ( x - A ) * ( x - A ) + ( y - B ) * ( y - B ) ) - rad ;

        if ( val_x < value )   value = val_x;
        if ( val_y < value )   value = val_y;
        if ( val_z < value )   value = val_z;
      }
    }

    levelset.set ( *bit, value );

  }

  pb.finish();
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generate3DPlaterodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods ) {
  generate3DPlaterodsLevelset ( levelset, n_rods, 0.12 / ( 1.0 * n_rods ) );
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generate3DRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods, const DataType rad ) {
  const int width = levelset.getNumXYZ();
  aol::ProgressBar<> pb ( "Creating 3DRods levelset" );
  pb.start ( aol::Cub ( width ) );

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const int i = ( *bit ) [0], j = ( *bit ) [1], k = ( *bit ) [2];
    pb++;

    const DataType x = ( 1.0 * i ) / ( width - 1.0 );
    const DataType y = ( 1.0 * j ) / ( width - 1.0 );
    const DataType z = ( 1.0 * k ) / ( width - 1.0 );

    DataType value = 1.0e20;

    for ( int ir = 0; ir < n_rods; ++ir ) {
      for ( int jr = 0; jr < n_rods; ++jr ) {
        const DataType A = 0.0 + 1.0 * ( 2.0 * ir + 1.0 ) / ( 2.0 * n_rods );
        const DataType B = 0.0 + 1.0 * ( 2.0 * jr + 1.0 ) / ( 2.0 * n_rods );
        const DataType val_x = sqrt ( ( y - A ) * ( y - A ) + ( z - B ) * ( z - B ) ) - rad ;
        const DataType val_y = sqrt ( ( x - A ) * ( x - A ) + ( z - B ) * ( z - B ) ) - rad ;
        const DataType val_z = sqrt ( ( x - A ) * ( x - A ) + ( y - B ) * ( y - B ) ) - rad ;

        if ( val_x < value )   value = val_x;
        if ( val_y < value )   value = val_y;
        if ( val_z < value )   value = val_z;
      }
    }

    levelset.set ( *bit, value );

  }

  pb.finish();
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generate3DRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods ) {
  generate3DRodsLevelset ( levelset, n_rods, 0.12 / ( 1.0 * n_rods ) );
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateAnisoSelected3DRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods,
                                                                             const DataType x_rad, const DataType y_rad, const DataType z_rad,
                                                                             const qc::BitArray<qc::QC_3D> &selMask_x, const qc::BitArray<qc::QC_3D> &selMask_y, const qc::BitArray<qc::QC_3D> &selMask_z ) {

  const int width = levelset.getNumXYZ(), n_cells = n_rods + 1;
  aol::ProgressBar<> pb ( "Selected_3DRods levelset" );
  pb.start ( aol::Cub ( width ) );

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const int i = ( *bit ) [0], j = ( *bit ) [1], k = ( *bit ) [2];
    pb++;

    const DataType x = ( 1.0 * i ) / ( width - 1.0 ), y = ( 1.0 * j ) / ( width - 1.0 ), z = ( 1.0 * k ) / ( width - 1.0 );

    DataType dist_value = 1.0;//1.0e20;

    // rods in x direction
    {
      const int ix = ( i * n_rods + width / 2 ) / width; // integer division
      for ( int jy = 0; jy < n_cells - 1; ++jy ) {
        for ( int kz = 0; kz < n_cells - 1; ++kz ) {
          if ( selMask_x.get ( ix, jy, kz ) == true ) {
            const DataType B = ( 2.0 * jy + 1.0 ) / ( 2.0 * n_rods );
            const DataType C = ( 2.0 * kz + 1.0 ) / ( 2.0 * n_rods );
            const DataType val_x = aol::Sqr ( y - B ) + aol::Sqr ( z - C ) - aol::Sqr ( x_rad );

            if ( val_x < dist_value )
              dist_value = val_x;

          }
        }
      }
    }

    // rods in y direction
    int ix;
    for ( ix = 0; ix < n_cells - 1; ++ix ) {
      const int jy = ( j * n_rods + width / 2 ) / width; { // integer division
        for ( int kz = 0; kz < n_cells - 1; ++kz ) {
          if ( selMask_y.get ( ix, jy, kz ) == true ) {
            const DataType A = ( 2.0 * ix + 1.0 ) / ( 2.0 * n_rods );
            const DataType C = ( 2.0 * kz + 1.0 ) / ( 2.0 * n_rods );
            const DataType val_y = aol::Sqr ( x - A ) + aol::Sqr ( z - C ) - aol::Sqr ( y_rad );

            if ( val_y < dist_value )
              dist_value = val_y;
          }
        }
      }
    }

    // rods in z direction
    for ( ix = 0; ix < n_cells - 1; ++ix ) {
      for ( int jy = 0; jy < n_cells - 1; ++jy ) {
        const int kz = ( k * n_rods + width / 2 ) / width; { // integer division
          if ( selMask_z.get ( ix, jy, kz ) == true ) {
            const DataType A = ( 2.0 * ix + 1.0 ) / ( 2.0 * n_rods );
            const DataType B = ( 2.0 * jy + 1.0 ) / ( 2.0 * n_rods );
            const DataType val_z = aol::Sqr ( x - A ) + aol::Sqr ( y - B ) - aol::Sqr ( z_rad );

            if ( val_z < dist_value )
              dist_value = val_z;
          }
        }
      }
    }

    levelset.set ( *bit, dist_value );
  }

  pb.finish();
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateSelected3DRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods, const DataType rad,
                                                                        const qc::BitArray<qc::QC_3D> &selMask_x, const qc::BitArray<qc::QC_3D> &selMask_y, const qc::BitArray<qc::QC_3D> &selMask_z ) {
  generateAnisoSelected3DRodsLevelset ( levelset, n_rods, rad, rad, rad, selMask_x, selMask_y, selMask_z );
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateRandomRodsRemoved3DRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods, const DataType rad, const DataType removeFraction ) {
  qc::BitArray<qc::QC_3D> selMask[3];

  aol::RandomGenerator rg;

  aol::ProgressBar<> pb ( "Selecting rods to remove" );
  pb.start ( 3 * static_cast<int> ( removeFraction * n_rods * n_rods * ( n_rods + 1 ) ) );
  int num_removed = 0;

  for ( int i = 0 ; i < 3; ++i ) {

    selMask[i].reallocate ( n_rods + 2, n_rods + 2, n_rods + 2 );
    selMask[i].setAll ( true );

    for ( int j = 0; j < static_cast<int> ( removeFraction * n_rods * n_rods * ( n_rods + 1 ) ); ++j ) {
      pb++;
      qc::CoordType pos ( rg.rInt ( n_rods + 1 ), rg.rInt ( n_rods + 1 ), rg.rInt ( n_rods + 1 ) );
      while ( selMask[i].get ( pos ) == false ) { // want to make sure that entry is not set yet so that we don't remove less entries than removeFraction (birthday paradoxon!)
        pos = qc::CoordType ( rg.rInt ( n_rods + 0 ), rg.rInt ( n_rods + 0 ), rg.rInt ( n_rods + 0 ) );
      }
      selMask[i].set ( pos, false );
    }

    num_removed += selMask[i].numFalse();

  }

  pb.finish();

  cerr << num_removed << " rods to be removed" << endl;

  generateSelected3DRodsLevelset ( levelset, n_rods, rad, selMask[0], selMask[1], selMask[2] );
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateAnisoRandomRodsRemoved3DRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods,
                                                                                      const DataType x_rad, const DataType y_rad, const DataType z_rad,
                                                                                      const DataType x_removeFraction, const DataType y_removeFraction, const DataType z_removeFraction,
                                                                                      const unsigned int seed ) {
  qc::BitArray<qc::QC_3D> selMask[3];
  aol::Vec3<DataType> removeFraction ( x_removeFraction, y_removeFraction, z_removeFraction );

  aol::RandomGenerator rg ( seed );

  int num_removed = 0;

  for ( int i = 0 ; i < 3; ++i ) {
    const int size_x = n_rods + ( i == 0 ? 1 : 0 ), size_y = n_rods + ( i == 1 ? 1 : 0 ), size_z = n_rods + ( i == 2 ? 1 : 0 );

    selMask[i].reallocate ( size_x, size_y, size_z );
    selMask[i].setAll ( true );
    for ( int j = 0; j < static_cast<int> ( removeFraction[i] * n_rods * n_rods * ( n_rods + 1 ) ); ++j ) {
      qc::CoordType pos ( rg.rInt ( size_x ), rg.rInt ( size_y ), rg.rInt ( size_z ) );
      while ( selMask[i].get ( pos ) == false ) { // want to make sure that entry is not set yet so that we don't remove less entries than removeFraction (birthday paradoxon!)
        pos = qc::CoordType ( rg.rInt ( size_x ), rg.rInt ( size_y ), rg.rInt ( size_z ) );
      }
#ifdef VERBOSE
      cerr << "Direction " << i << ", removing " << pos << endl;
#endif
      selMask[i].set ( pos, false );
    }

    num_removed += selMask[i].numFalse();

  }

  cerr << num_removed << " rods to be removed" << endl;

  generateAnisoSelected3DRodsLevelset ( levelset, n_rods, x_rad, y_rad, z_rad, selMask[0], selMask[1], selMask[2] );

}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generatePeriodicAnisoRandom3DRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods,
                                                                                   const DataType x_rad, const DataType y_rad, const DataType z_rad,
                                                                                   const DataType x_removeFraction, const DataType y_removeFraction, const DataType z_removeFraction,
                                                                                   const unsigned int seed ) {
  qc::BitArray<qc::QC_3D> selMask[3];
  aol::Vec3<DataType> removeFraction ( x_removeFraction, y_removeFraction, z_removeFraction );

  aol::RandomGenerator rg ( seed );

  int num_removed = 0;

  for ( int i = 0 ; i < 3; ++i ) {
    selMask[i].reallocate ( n_rods, n_rods, n_rods );
    selMask[i].setAll ( true );
    for ( int j = 0; j < static_cast<int> ( removeFraction[i] * aol::Cub ( n_rods ) ); ++j ) {
      qc::CoordType pos ( rg.rInt ( n_rods ), rg.rInt ( n_rods ), rg.rInt ( n_rods ) );
      while ( selMask[i].get ( pos ) == false ) { // want to make sure that entry is not set yet so that we don't remove less entries than removeFraction (birthday paradoxon!)
        pos = qc::CoordType ( rg.rInt ( n_rods ), rg.rInt ( n_rods ), rg.rInt ( n_rods ) );
      }
#ifdef VERBOSE
      cerr << "Direction " << i << ", removing " << pos << endl;
#endif
      selMask[i].set ( pos, false );
    }

    num_removed += selMask[i].numFalse();

  }

  selMask[0].resize ( n_rods + 1, n_rods    , n_rods     );
  selMask[1].resize ( n_rods    , n_rods + 1, n_rods     );
  selMask[2].resize ( n_rods    , n_rods    , n_rods + 1 );

  for ( int m = 0; m < n_rods; ++m ) {
    for ( int n = 0; n < n_rods; ++n ) {
      if ( selMask[0].get ( 0, m, n ) ) {
        selMask[0].set ( n_rods, m, n, true );
      }

      if ( selMask[1].get ( m, 0, n ) ) {
        selMask[1].set ( m, n_rods, n, true );
      }
      if ( selMask[2].get ( m, n, 0 ) ) {
        selMask[2].set ( m, n, n_rods, true );
      }
    }
  }

  cerr << "tpcfe::generatePeriodicAnisoRandom3DRodsLevelset: " << num_removed << " rods to be removed" << endl;

  generateAnisoSelected3DRodsLevelset ( levelset, n_rods, x_rad, y_rad, z_rad, selMask[0], selMask[1], selMask[2] );

}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateSolidcubeLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset ) {
  levelset.setAll ( -1.0 );
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateHalfcubeLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const DataType threshold, const qc::Comp direction ) {

  const DataType thres = threshold * ( levelset.getSize() [ direction ] - 1.0 );

#ifdef VERBOSE
  cerr << "halfcube: " << levelset.getSize() [ direction ] << " " << thres << endl;
#endif

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    levelset.set ( *bit, 1.0 * ( *bit ) [ direction ] - thres );
  }

}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateRotatedHalfcubeLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const DataType threshold, const double angle ) {

#ifdef VERBOSE
  cerr << "halfcube: " << levelset.getNumY() << " " << threshold << endl;
#endif

  aol::Matrix33<DataType> rotMat;
  rotMat.setRotationAboutX ( angle );
  const aol::Vec3<DataType> ctr ( 0.5, 0.5, 0.5 );
  const int width = levelset.getNumY();

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    aol::Vec3<DataType> pos ( ( *bit ) [0] / ( width - 1.0 ), ( *bit ) [1] / ( width - 1.0 ), ( *bit ) [2] / ( width - 1.0 ) ), rotPos;
    pos -= ctr;
    rotPos = rotMat * pos;
    rotPos += ctr;

    levelset.set ( *bit, rotPos[1] - threshold );
  }

}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateSlotLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const DataType width, const DataType depth ) {

  const DataType wi = width * levelset.getNumZ(), zmid = 0.5 * ( levelset.getNumZ() - 1.0 ), de = depth * levelset.getNumX() - wi ;

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const int i = ( *bit ) [0], k = ( *bit ) [2];
    DataType value = wi / 2;
    if ( i > de ) {
      value -= ( i - de ) ;
    }
    value -= aol::Abs ( k - zmid );
    levelset.set ( *bit,  value );
  }

}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateLaminateLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int number_of_layers, const DataType alpha_deg, const DataType beta_deg, const bool periodic ) {
  const int  width = levelset.getNumXYZ();
  const DataType N = static_cast<DataType> ( width );
  const DataType alpha = aol::DegreesToRadians ( alpha_deg );
  const DataType beta  = aol::DegreesToRadians (  beta_deg );

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const int i = ( *bit ) [0], j = ( *bit ) [1], k = ( *bit ) [2];

    aol::Vec3<DataType> pos ( 1.0 * i, 1.0 * j, 1.0 * k ), pos_rot ( 0.0, 0.0, 0.0 );
    pos /= N ; // pos is in world coordinates

    aol::Matrix33<DataType>
      rot_x ( cos ( alpha ), -sin ( alpha ), 0., sin ( alpha ), cos ( alpha ), 0., 0., 0., 1. ),
      rot_y ( 1., 0., 0., 0., sin ( beta ), cos ( beta ), 0., -cos ( beta ), sin ( beta ) ),
      rot = rot_x;

    rot *= rot_y;

    pos_rot = rot * pos;

    const DataType val = cos ( number_of_layers * aol::NumberTrait<DataType>::pi * ( periodic ?  pos_rot[0] - 0.5 : pow ( aol::Abs ( static_cast<double>(pos_rot[0]) - 0.5 ), 1.5 ) ) );
    levelset.set ( *bit, val );

  }
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateHalfcubeRippleLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const DataType threshold,
                                                                        const aol::Vec3<DataType> a, const aol::Vec3<DataType> k ) {
  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const DataType x = ( *bit ) [0] / ( levelset.getNumX() - 1.0 ), y = ( *bit ) [1] / ( levelset.getNumY() - 1.0 ), z = ( *bit ) [2] / ( levelset.getNumZ() - 1.0 );
    const DataType pi = aol::NumberTrait<DataType>::pi;
    const DataType value = y -  a[0] * ( sin ( k[0] * 2 * pi * x ) + sin ( k[0] * 2 * pi * z ) ) - a[1] * ( sin ( k[1] * 2 * pi * x ) + sin ( k[1] * 2 * pi * z ) ) - a[2] * ( sin ( k[2] * 2 * pi * x ) + sin ( k[2] * 2 * pi * z ) ) - threshold;
    levelset.set ( *bit, value );
  }
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateTwoPlaneRippleLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset,
                                                                        const aol::Vec3<DataType> a0, const aol::Vec3<DataType> a1, const aol::Vec3<DataType> k0, const aol::Vec3<DataType> k1 ) {
  qc::ScalarArray<DataType, qc::QC_3D> ls1 ( levelset, aol::STRUCT_COPY ), ls2 ( levelset, aol::STRUCT_COPY );
  qc::ShapeLevelsetGenerator<DataType>::generateHalfcubeRippleLevelset ( ls1, 1./3., a0, k0 );
  qc::ShapeLevelsetGenerator<DataType>::generateHalfcubeRippleLevelset ( ls2, 2./3., a1, k1 );
  for ( int i = 0; i < levelset.size(); ++i ) {
    levelset[i] = - aol::Min ( ls1[i], - ls2[i] );
  }
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateFourPlaneRippleLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset,
                                                                         const aol::Vec3<DataType> a0, const aol::Vec3<DataType> a1, const aol::Vec3<DataType> a2, const aol::Vec3<DataType> a3,
                                                                         const aol::Vec3<DataType> k0, const aol::Vec3<DataType> k1, const aol::Vec3<DataType> k2, const aol::Vec3<DataType> k3 ) {
  qc::ScalarArray<DataType, qc::QC_3D> ls1 ( levelset, aol::STRUCT_COPY ), ls2 ( levelset, aol::STRUCT_COPY ), ls3 ( levelset, aol::STRUCT_COPY ), ls4 ( levelset, aol::STRUCT_COPY );
  qc::ShapeLevelsetGenerator<DataType>::generateHalfcubeRippleLevelset ( ls1, 1./7., a0, k0 );
  qc::ShapeLevelsetGenerator<DataType>::generateHalfcubeRippleLevelset ( ls2, 2./7., a1, k1 );
  qc::ShapeLevelsetGenerator<DataType>::generateHalfcubeRippleLevelset ( ls3, 5./7., a2, k2 );
  qc::ShapeLevelsetGenerator<DataType>::generateHalfcubeRippleLevelset ( ls4, 6./7., a3, k3 );

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    if ( (*bit)[1] < levelset.getNumY() / 2 ) {
      levelset.set( *bit, - aol::Min ( ls1.get(*bit), - ls2.get(*bit) ) );
    } else {
      levelset.set( *bit, - aol::Min ( ls3.get(*bit), - ls4.get(*bit) ) );
    }
  }
}



template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateSwissCheeseLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_holes, const DataType radius ) {

  const int width = levelset.getNumXYZ();
  aol::ProgressBar<> pb ( "Creating SwissCheese levelset" );
  pb.start ( aol::Cub ( width ) );

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const int i = ( *bit ) [0], j = ( *bit ) [1], k = ( *bit ) [2];
    pb++;

    const DataType x = ( 1.0 * i ) / ( width - 1.0 ), y = ( 1.0 * j ) / ( width - 1.0 ), z = ( 1.0 * k ) / ( width - 1.0 );

    DataType value = -1e20;

    for ( qc::RectangularIterator<qc::QC_3D> holeit ( aol::Vec3<short> ( 0, 0, 0 ), aol::Vec3<short> ( n_holes, n_holes, n_holes ) ); holeit.notAtEnd(); ++holeit ) {
      const int ir = ( *holeit ) [0], jr = ( *holeit ) [1], kr = ( *holeit ) [2];
      const DataType A = 0.1 + 0.8 * ( 2.0 * ir + 1.0 ) / ( 2.0 * n_holes );
      const DataType B = 0.1 + 0.8 * ( 2.0 * jr + 1.0 ) / ( 2.0 * n_holes );
      const DataType C = 0.1 + 0.8 * ( 2.0 * kr + 1.0 ) / ( 2.0 * n_holes );
      const DataType val = radius - sqrt ( ( x - A ) * ( x - A ) + ( y - B ) * ( y - B ) + ( z - C ) * ( z - C ) ) ;

      if ( val > value )
        value = val;
    }

    levelset.set ( *bit, value );

  }

  pb.finish();

}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateAustrianCheeseLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_holes, const DataType radius ) {

  const int width = levelset.getNumXYZ();
  aol::ProgressBar<> pb ( "Creating AustrianCheese levelset" );
  pb.start ( aol::Cub ( width ) );

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const int i = ( *bit ) [0], j = ( *bit ) [1], k = ( *bit ) [2];
    pb++;

    const DataType x = ( 1.0 * i ) / ( width - 1.0 ), y = ( 1.0 * j ) / ( width - 1.0 ), z = ( 1.0 * k ) / ( width - 1.0 );

    DataType value = -1e20;

    for ( qc::RectangularIterator<qc::QC_3D> holeit ( aol::Vec3<short> ( 0, 0, 0 ), aol::Vec3<short> ( n_holes, n_holes, n_holes ) ); holeit.notAtEnd(); ++holeit ) {
      const int ir = ( *holeit ) [0], jr = ( *holeit ) [1], kr = ( *holeit ) [2];
      const DataType A = ( 1.0 * ir ) / ( n_holes - 1.0 ),  B = ( 1.0 * jr ) / ( n_holes - 1.0 ), C = ( 1.0 * kr ) / ( n_holes - 1.0 );
      const DataType val = radius - sqrt ( ( x - A ) * ( x - A ) + ( y - B ) * ( y - B ) + ( z - C ) * ( z - C ) ) ;

      if ( val > value )
        value = val;
    }

    levelset.set ( *bit, value );

  }

  pb.finish();

}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateRandomRadiiSwissCheeseLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_holes ) {

  const int width = levelset.getNumXYZ();
  aol::ProgressBar<> pb ( "Creating SwissCheese levelset" );
  pb.start ( aol::Cub ( width ) );

  qc::ScalarArray<DataType, qc::QC_3D> radii ( n_holes, n_holes, n_holes );

  {

    aol::RandomGenerator rg;
    for ( qc::RectangularIterator<qc::QC_3D> holeit ( aol::Vec3<short> ( 0, 0, 0 ), aol::Vec3<short> ( n_holes, n_holes, n_holes ) ); holeit.notAtEnd(); ++holeit ) {
      const int ir = ( *holeit ) [0], jr = ( *holeit ) [1], kr = ( *holeit ) [2];
      const DataType radius = rg.rReal ( 2.0 / ( width - 1.0 ), 0.3 / ( 1.0 * n_holes ) * 1.0 );
      radii.set ( ir, jr, kr, radius );
#ifdef VERBOSE
      cerr << radius << endl;
#endif
    }

  }

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const int i = ( *bit ) [0], j = ( *bit ) [1], k = ( *bit ) [2];
    pb++;

    const DataType x = ( 1.0 * i ) / ( width - 1.0 ), y = ( 1.0 * j ) / ( width - 1.0 ), z = ( 1.0 * k ) / ( width - 1.0 );

    DataType value = -1e20;

    for ( qc::RectangularIterator<qc::QC_3D> holeit ( aol::Vec3<short> ( 0, 0, 0 ), aol::Vec3<short> ( n_holes, n_holes, n_holes ) ); holeit.notAtEnd(); ++holeit ) {
      const int ir = ( *holeit ) [0], jr = ( *holeit ) [1], kr = ( *holeit ) [2];
      const DataType A = 0.1 + 0.8 * ( 2.0 * ir + 1.0 ) / ( 2.0 * n_holes );
      const DataType B = 0.1 + 0.8 * ( 2.0 * jr + 1.0 ) / ( 2.0 * n_holes );
      const DataType C = 0.1 + 0.8 * ( 2.0 * kr + 1.0 ) / ( 2.0 * n_holes );
      const DataType radius = radii.get ( ir, jr, kr );
      const DataType val = radius - sqrt ( ( x - A ) * ( x - A ) + ( y - B ) * ( y - B ) + ( z - C ) * ( z - C ) ) ;

      if ( val > value )
        value = val;
    }

    levelset.set ( *bit, value );

  }

  pb.finish();

}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateSwissCheeseLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_holes ) {
  generateSwissCheeseLevelset ( levelset, n_holes,  0.3 / ( 1.0 * n_holes ) );
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::superimposeBoundingBox ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int numBL, const DataType LSValueThere ) {
  const int nX = levelset.getNumX(), nY = levelset.getNumY(), nZ = levelset.getNumZ();
  for ( qc::RectangularIterator<qc::QC_3D> it ( levelset ); !it.atEnd(); ++it ) {
    qc::CoordType pos = *it;
    if ( pos[0] < numBL || pos[1] < numBL || pos[2] < numBL || pos[0] > nX - numBL - 1 || pos[1] > nY - numBL - 1 || pos[2] > nZ - numBL - 1 ) {
      levelset.set ( pos, LSValueThere );
    }
  }
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateZLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset ) {
  const int width = levelset.getNumXYZ();
  const DataType h = 1.0 / ( width - 1.0 );

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const aol::Vec3<DataType>
    pos ( h * ( *bit ) [0], h * ( *bit ) [1], h * ( *bit ) [2] ),
        proj0 ( 0.5 - 0.6 * pos[1], pos[1], 0.5 ),
        proj1 ( pos[0], 0.5, 0.5 ),
        proj2 ( 1.1 - 0.6 * pos[1], pos[1], 0.5 );

    const DataType dist0 = ( pos[1] < 0.5 ?  ( pos - proj0 ).norm()  :  aol::NumberTrait<DataType>::Inf ) - 0.1;
    const DataType dist1 = sqrt ( 2 * aol::Sqr ( pos[1] - proj1[1] ) + aol::Sqr ( pos[2] - proj1[2] ) ) - 0.1;
    const DataType dist2 = ( pos[1] > 0.5 ?  ( pos - proj2 ).norm()  :  aol::NumberTrait<DataType>::Inf ) - 0.1;

    levelset.set ( *bit, aol::Min ( dist0, dist1, dist2 ) );
  }
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateZ3Levelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset ) {
  const int width = levelset.getNumXYZ();
  const DataType h = 1.0 / ( width - 1.0 );

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const aol::Vec3<DataType>
    pos ( h * ( *bit ) [0], h * ( *bit ) [1], h * ( *bit ) [2] ),
        proj0 ( 0.5 - 0.6 * pos[1], pos[1], pos[2] ),
        proj1 ( pos[0], 0.5, pos[2] ),
        proj2 ( 1.1 - 0.6 * pos[1], pos[1], pos[2] );

    const DataType dist0 = ( pos[1] < 0.5 ?  ( pos - proj0 ).norm()  :  aol::NumberTrait<DataType>::Inf );
    const DataType dist1 = sqrt ( 2 * aol::Sqr ( pos[1] - proj1[1] ) + aol::Sqr ( pos[2] - proj1[2] ) );
    const DataType dist2 = ( pos[1] > 0.5 ?  ( pos - proj2 ).norm()  :  aol::NumberTrait<DataType>::Inf );

    levelset.set ( *bit, aol::Min ( dist0, dist1, dist2 ) - 0.11 );

  }
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateDrunkenZebraLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset ) {
  const int width = levelset.getNumXYZ();
  const DataType h = 1.0 / ( width - 1.0 );

  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
    const aol::Vec3<DataType>
    pos ( h * ( *bit ) [0], h * ( *bit ) [1], h * ( *bit ) [2] ),
        proj0 ( 0.5 - pos[1], pos[1], pos[2] ),
        proj1 ( 1.5 - pos[1], pos[1], pos[2] );

    const DataType dist0 = ( pos - proj0 ).norm();
    const DataType dist1 = ( pos - proj1 ).norm();

    levelset.set ( *bit, aol::Min ( dist0, dist1 ) - 0.11 );

  }
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateHoneycombLevelset ( qc::ScalarArray<DataType, qc::QC_3D> & levelset, const int n, const DataType theta ) {
  levelset.setAll ( aol::NumberTrait<DataType>::Inf );
  const int m = static_cast<int> ( floor ( 2 * n / sqrt ( 3. ) ) );
  if ( m % 2 == 0 ) {
    const DataType h = aol::NumberTrait<DataType>::one / ( levelset.getNumXYZ() - 1 );
    const DataType yoffset = ( 1 - m * sqrt ( 3. ) / ( 2 * n ) ) / 2;

    aol::ProgressBar<> pb ( "Generating Honeycomb Levelset" );
    pb.start ( ( n + 1 ) * ( m + 1 ) + n * m );

    for ( int i = 0; i <= n; ++i ) {
      for ( int j = 0; j <= m; ++j, pb++ ) {
        for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
          const aol::Vec3<DataType> pointOffset ( h * ( *bit ) [0] - static_cast<DataType> ( i ) / n, h * ( *bit ) [1] - ( sqrt ( 3. ) * j ) / n - yoffset, aol::NumberTrait<DataType>::zero );

          const DataType value = - aol::NumberTrait<DataType>::one / ( 2 * n ) + theta / 2 + pointOffset.norm();

          if ( value < levelset.get ( *bit ) )
            levelset.set ( *bit, value );

        }
      }
    }

    for ( int i = 0; i < n; ++i ) {
      for ( int j = 0; j < m; ++j, pb++ ) {
        for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
          const aol::Vec3<DataType> pointOffset ( h * ( *bit ) [0] - static_cast<DataType> ( 2 * i + 1 ) / ( 2 * n ), h * ( *bit ) [1] - ( sqrt ( 3. ) * ( 2 * j + 1 ) ) / ( 2 * n ) - yoffset, aol::NumberTrait<DataType>::zero );

          const DataType value = - aol::NumberTrait<DataType>::one / ( 2 * n ) + theta / 2 + pointOffset.norm();

          if ( value < levelset.get ( *bit ) )
            levelset.set ( *bit, value );

        }
      }
    }
    pb.finish();

    levelset *= -1.0;
  } else {
    throw aol::Exception ( "tpcfe::generateHoneycombLevelset: this number of holes in one direction will not produce a periodic pattern", __FILE__, __LINE__ );
  }
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateGeomOrthoHoneycombLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n, const aol::Mat<2,2,DataType> &thetaMat ) {
  levelset.setAll ( aol::NumberTrait<DataType>::Inf );

  const int m = static_cast<int> ( floor ( 2 * n / sqrt ( 3. ) ) );
  if ( m % 2 == 0 ) {
    const DataType h = aol::NumberTrait<DataType>::one / ( levelset.getNumXYZ() - 1 );
    const DataType yoffset = ( 1 - m * sqrt ( 3. ) / ( 2 * n ) ) / 2;

    aol::ProgressBar<> pb ( "Generating Honeycomb Levelset" );
    pb.start ( (n+1)*(m+1) + n*m );

    for ( int i = 0; i <= n; ++i ) {
      for ( int j = 0; j <= m; ++j, pb++ ) {
        for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
          const aol::Vec3<DataType> pointOffset ( h * ( *bit ) [0] - static_cast<DataType> ( i ) / n, h * ( *bit ) [1] - ( sqrt ( 3. ) * j ) / n - yoffset, aol::NumberTrait<DataType>::zero );

          const DataType theta = thetaMat[ ( (i%2) + (j%2) ) % 2 ][ 0 ];
          const DataType value = - aol::NumberTrait<DataType>::one / ( 2 * n ) + theta / 2 + pointOffset.norm();

          if ( value < levelset.get ( *bit ) )
            levelset.set ( *bit, value );

        }
      }
    }

    for ( int i = 0; i < n; ++i ) {
      for ( int j = 0; j < m; ++j, pb++ ) {
        for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit ) {
          const aol::Vec3<DataType> pointOffset ( h * ( *bit ) [0] - static_cast<DataType> ( 2*i + 1 ) / ( 2*n ), h * ( *bit ) [1] - ( sqrt ( 3. ) * ( 2*j + 1 ) ) / ( 2*n ) - yoffset, aol::NumberTrait<DataType>::zero );

          const DataType theta = thetaMat[ ( (i%2) + (j%2) ) % 2 ][ 1 ];
          const DataType value = - aol::NumberTrait<DataType>::one / ( 2 * n ) + theta / 2 + pointOffset.norm();

          if ( value < levelset.get ( *bit ) )
            levelset.set ( *bit, value );

        }
      }
    }
    pb.finish();

    levelset *= -1.0;
  } else {
    throw aol::Exception ( "tpcfe::generateHoneycombLevelset: this number of holes in one direction will not produce a periodic pattern", __FILE__, __LINE__ );
  }
}


template< typename DataType >
void ShapeLevelsetGenerator<DataType>::generateXRotated3DRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods, const aol::Vec3<DataType> &radii, const DataType angle ) {
  aol::Matrix33<DataType> rotMat;
  rotMat.setRotationAboutX ( angle );

  rotMat /= cos ( angle );

  const int width = levelset.getNumXYZ();
  aol::ProgressBar<> pb ( "Creating Rotated 3DRods levelset" );
  pb.start ( aol::Cub ( width ) );

  const aol::Vec3<DataType> ctr ( 0.5, 0.5, 0.5 );
  for ( qc::RectangularIterator<qc::QC_3D> bit ( levelset ); bit.notAtEnd(); ++bit, pb++ ) {
    aol::Vec3<DataType> pos ( ( *bit ) [0] / ( width - 1.0 ), ( *bit ) [1] / ( width - 1.0 ), ( *bit ) [2] / ( width - 1.0 ) ), rotPos;
    pos -= ctr;
    rotPos = rotMat * pos;
    rotPos += ctr;

    DataType value = 1e20;

    for ( int ir = - n_rods; ir < 2 * n_rods; ++ir ) {
      for ( int jr = - n_rods; jr < 2 * n_rods; ++jr ) {
        const DataType A = 0.0 + 1.0 * ( 2.0 * ir + 1.0 ) / ( 2.0 * n_rods );
        const DataType B = 0.0 + 1.0 * ( 2.0 * jr + 1.0 ) / ( 2.0 * n_rods );
        const DataType val_x = aol::Sqr ( rotPos[1] - A ) + aol::Sqr ( rotPos[2] - B ) - aol::Sqr ( radii[0] );
        const DataType val_y = aol::Sqr ( rotPos[0] - A ) + aol::Sqr ( rotPos[2] - B ) - aol::Sqr ( radii[1] );
        const DataType val_z = aol::Sqr ( rotPos[0] - A ) + aol::Sqr ( rotPos[1] - B ) - aol::Sqr ( radii[2] );

        if ( val_x < value )   value = val_x;
        if ( val_y < value )   value = val_y;
        if ( val_z < value )   value = val_z;
      }
    }

    levelset.set ( *bit, value );

  }

  pb.finish();
}

template class ShapeLevelsetGenerator<float>;
template class ShapeLevelsetGenerator<double>;
template class ShapeLevelsetGenerator<long double>;

// end namespace
}
