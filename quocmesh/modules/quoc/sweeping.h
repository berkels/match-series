#ifndef __SWEEPING_H
#define __SWEEPING_H

#include <array.h>
#include <gridBase.h>
#include <quoc.h>
#include <scalarArray.h>

namespace qc {

template< typename DataType >
class DistanceSweeper2d {

protected:

  const BitArray<qc::QC_2D>            &_isSeed;      // true if corresponding node is seed node; local data structure

  ScalarArray<DataType, qc::QC_2D>     _phiOld, _phiNew;

  const DataType              _h;

  const DataType              _eps;
  const DataType              _smallValue;   // in updateOnSide
  bool                        _verbose;

  const GridDefinition        &_grid;

public:
  //! Constructor where only the seed points are specified and all points are implicitely assumed to be relevant
  DistanceSweeper2d ( const BitArray<qc::QC_2D> &isSeed, const GridDefinition &grid )
      : _isSeed ( isSeed ), _phiOld( ), _phiNew( ), _h ( grid.H() ), _eps ( 1.0e-5 ), _smallValue ( 1.0e-10 ), _verbose ( false ), _grid ( grid ) {}


  ~DistanceSweeper2d ( ) {}

public:
  /** Compute distances for all nodes in the grid
   *  At seed points, initial distances must be given.
   *  At non-seed points, ??? -> same below.
   *  \param delta       stopping criterion for the sweeping algorithm
   *  \param max_steps   maximum number of sweeps
   */
  void computeDistances ( ScalarArray<DataType, qc::QC_2D> &distanceField, const DataType delta = 1e-12, const int max_steps = 500 ) {
    _phiOld.reallocate ( _grid.getWidth(), _grid.getHeight() ); // automatically set to zero
    _phiNew.reallocate ( _grid.getWidth(), _grid.getHeight() );
    _phiNew = distanceField;



    DataType difference = aol::NumberTrait<DataType>::one; // is this the correct initialization?
    int counter = 0;

    DataType theta = 0.0;

    while ( difference > delta && counter < max_steps ) {

      _phiOld = _phiNew;

      for ( int s1 = -1; s1 <= 1;s1 += 2 ) {
        for ( int s2 = -1; s2 <= 1; s2 += 2 ) {

          for ( int i = ( s1 < 0 ? ( _grid.getWidth() - 1 ) : 0 ); ( s1 < 0 ? i >= 0 : i <= ( _grid.getWidth() - 1 ) ); i += s1 ) {
            for ( int j = ( s2 < 0 ? ( _grid.getHeight() - 1 ) : 0 ); ( s2 < 0 ? j >= 0 : j <= ( _grid.getHeight() - 1 ) ); j += s2 ) {

              if ( ! ( _isSeed.get ( i, j ) ) ) {

                for ( int sx = -1;  sx <= 1 ;sx += 2 ) {
                  for ( int sy = -1;  sy <= 1 ;sy += 2 ) {

                    if ( checkBx ( i - sx ) && checkBy ( j - sy ) ) {

                      int sign;
                      DataType phi_tmp = aol::NumberTrait<DataType>::Inf;

                      if ( _phiNew.get ( i, j - sy ) < 0 || _phiNew.get ( i - sx, j ) < 0 ) {
                        sign = -1;
                      } else {
                        sign = 1;
                      }

                      DataType m = sx * sy * ( std::abs ( _phiNew.get ( i, j - sy ) ) - std::abs ( _phiNew.get ( i - sx, j ) ) ) / _h;

                      if ( sx*sy != 0 && -1 < m && m < 1 ) {

                        const DataType
                          dummy1 = -sx * sy + m * sqrt ( 2 - m * m ),
                          dummy2 = -sx * sy - m * sqrt ( 2 - m * m ),
                          denom = m * m - sx * sx;

                        if ( aol::Abs ( m ) > 1.e-15 && aol::Abs ( denom ) > _smallValue ) {

                          const DataType
                            frac1 = static_cast<DataType> ( dummy1 / denom ),
                            frac2 = static_cast<DataType> ( dummy2 / denom ),

                            theta1 = atan ( frac1 ) + ( 1. - sx ) * aol::NumberTrait<DataType>::pi * 0.5,
                            theta2 = atan ( frac2 ) + ( 1. - sx ) * aol::NumberTrait<DataType>::pi * 0.5,

                            T1 = -sx * sin ( theta1 ) + cos ( theta1 ) * sy,

                            test =  m - T1;

                          if ( test < _smallValue ) {
                            theta = theta1;
                          } else {
                            theta = theta2;
                          }

                        }// end if( abs(m)>1.e-15 && abs(denom)>1.e-15  )

                        if ( aol::Abs ( denom ) < _smallValue ) {
                          theta = 0. + ( 1 - sx ) * aol::NumberTrait<DataType>::pi * 0.5;
                        }

                        if ( aol::Abs ( m ) < _smallValue ) {
                          const DataType
                            frac = static_cast<DataType> ( sy / sx ),
                            sum_1 = atan ( frac ),
                            sum_2 = ( 1 - sx ) * aol::NumberTrait<DataType>::pi * 0.5;

                          theta = sum_1 + sum_2;
                        }

                        // Update phi

                        phi_tmp = ( sx * cos ( theta ) * std::abs ( _phiNew.get ( i - sx, j ) ) + sy * sin ( theta ) * std::abs ( _phiNew.get ( i, j - sy ) )  + _h ) / ( abs ( cos ( theta ) ) +  abs ( sin ( theta ) ) );

                      }// end  if(sx*sy!=0 && -1<m && m<1){

                      if ( checkBx ( i - sx ) ) {
                        const DataType K_0 = std::abs ( _phiNew.get ( i - sx, j ) ) + _h ;
                        if ( phi_tmp > K_0 ) {
                          phi_tmp = K_0;
                        }
                      }
                      if ( checkBy ( j - sy ) ) {
                        const DataType K_pi_2 = std::abs ( _phiNew.get ( i, j - sy ) )  + _h;
                        if ( K_pi_2 < phi_tmp ) {
                          phi_tmp =  K_pi_2;
                        }
                      }

                      if ( phi_tmp < std::abs ( _phiNew.get ( i, j ) ) ) {
                        _phiNew.set ( i, j, phi_tmp*sign );
                      }

                      phi_tmp = aol::NumberTrait<DataType>::Inf;

                    }// end if(checkBx(i-sx) && checkBy(j-sy)   )

                  }// end for(sy=-1;  sy <=1 ;sy +=2)
                }// end  for(sx=-1;  sx <=1 ;sx +=2){

              }//end if( ! ( _isSeed.get(i,j) ) )

            }//end for(j=(s2<0?(1/_h):0); (s2<0?j>=0:j<=(1/_h)); j+=s2 )
          }// end for(i=(s1<0?(1/_h):0); (s1<0?i>=0:i<=(1/_h)); i+=s1 )
        }// end  for( int s2=-1; s2<=1; s2+=2)
      }//end for( int s1 =-1; s1 <=1;s1 +=2){


#ifdef DEBUG
      for ( int i = 0; i < _grid.getWidth() ; ++i )
        for ( int j = 0; j < _grid.getHeight() ; ++j )
          if ( _phiNew.get ( i, j ) == aol::NumberTrait<DataType>::Inf )
            cerr << i << " " << j << "contains inf. " << endl;
      cerr << endl;
#endif

      _phiOld -= _phiNew;

      difference = _phiOld.norm();

      ++counter;

       }//end while (difference > delta)

    distanceField = _phiNew;

  }

  void setVerboseMode ( bool mode = true ) {
      _verbose = mode;
  }

 private:
  DistanceSweeper2d ( );

  DistanceSweeper2d ( const DistanceSweeper2d< DataType > &/*other*/ );

  DistanceSweeper2d< DataType >& operator= ( const DistanceSweeper2d< DataType > &/*other*/ );

 protected:
  inline bool checkBx ( const int x ) const {
    return ( x >= 0 && x < _grid.getWidth() );
  }

  inline bool checkBy ( const int y ) const {
    return ( y >= 0 && y < _grid.getHeight() );
  }


#ifdef VERBOSE
  inline void printIfInfinite ( const char* const nameOfVariable, const DataType &Value ) const {
    if ( !aol::isFinite ( Value ) ) {
      cerr << nameOfVariable << " is infinite" << endl;
    }
  }
#else
  inline void printIfInfinite ( const char* const, const DataType& ) const {}
#endif

  // }

};



/** Computation of distances based on the sweeping method, either on a full ScalarArray (default) or on local boxes
 *  Default behavior is full ScalarArray for which the two "Local" can be ignored, for local boxes, use qc::RectangularContainer
 *  The Sweeper keeps a reference to the BitArray of seeds and a pointer to the BitArray of relevant points (if applicable)
 *
 *  \author Schwen, Teusner
 */

template< typename DataType, typename LocalBitArrayType = BitArray<qc::QC_3D>, typename LocalDistanceArrayType = ScalarArray<DataType, qc::QC_3D>, typename _GridType = qc::GridDefinition >
class DistanceSweeper3d {
protected:
  typedef _GridType GridType;

  const LocalBitArrayType     &_isSeed;      // true if corresponding node is seed node; local data structure
  const BitArray<qc::QC_3D>* const     _pIsRelevant;  // NULL if all points are implicitely relevant; true if corresponding node is an interior node (thus relevant for distance computation); global data structure (for now, may need to be local)

  LocalDistanceArrayType      _phiOld, _phiNew;

  const DataType              _h;

  const DataType              _eps;
  const DataType              _smallValue;   // in updateOnSide
  bool                        _verbose;

  const GridType &_grid;

public:
  //! Constructor where only the seed points are specified and all points are implicitely assumed to be relevant
  DistanceSweeper3d ( const LocalBitArrayType &isSeed, const GridType &grid )
      : _isSeed ( isSeed ), _pIsRelevant ( NULL ), _phiOld( ), _phiNew( ), _h ( grid.H() ), _eps ( 1.0e-5 ), _smallValue ( 1.0e-10 ), _verbose ( false ), _grid ( grid ) {}

  //! Constructor where both seed points and relevant points are specified
  DistanceSweeper3d ( const LocalBitArrayType &isSeed, const BitArray<qc::QC_3D>* const pIsRelevant, const GridType &grid )
      : _isSeed ( isSeed ), _pIsRelevant ( pIsRelevant ), _phiOld( ), _phiNew( ), _h ( grid.H() ), _eps ( 1.0e-5 ), _smallValue ( 1.0e-10 ), _verbose ( false ), _grid ( grid ) {}

  ~DistanceSweeper3d ( ) {}

public:
  /** Compute distances for all nodes in the grid
   *  At seed points, initial distances must be given.
   *  At non-seed points, ??? -> same below.
   *  \param delta       stopping criterion for the sweeping algorithm
   *  \param max_steps   maximum number of sweeps
   */
  void computeDistances ( LocalDistanceArrayType &distanceField, const DataType delta = 1e-12, const int max_steps = 500 ) {
    computeDistances ( distanceField, delta, max_steps, CoordType ( 0, 0, 0 ), CoordType ( _grid.getNumX(), _grid.getNumY(), _grid.getNumZ() ) );
  }

  /** Compute distances for all nodes in the cuboid given by two vertices min_vertex and max_vertex (lower left front vertex included and upper right back vertex not included).
   *  At seed points, initial distances must be given.
   */
  void computeDistances ( LocalDistanceArrayType &distanceField, const CoordType &min_vertex, const CoordType &max_vertex ) {
    computeDistances ( distanceField, 1e-12, 500, min_vertex, max_vertex );
  }

  void computeDistances ( LocalDistanceArrayType &distanceField, const DataType delta, const int max_steps, const CoordType &min_vertex, const CoordType &max_vertex ) {
    _phiOld.reallocate ( distanceField ); // automatically set to zero
    _phiNew.reallocate ( distanceField );
    _phiNew = distanceField;

    DataType difference = aol::NumberTrait<DataType>::one, noninf_difference = aol::NumberTrait<DataType>::one; // is this the correct initialization?
    int counter = 0, num_inf = _phiNew.numOccurence ( aol::NumberTrait<DataType>::Inf ), num_inf_difference = 1;

    while ( difference > delta && ( num_inf_difference != 0 || noninf_difference > delta )  && counter < max_steps ) {

      if ( _verbose ) {
        aol::SimpleFormat format ( 4, 0 );
        cerr << endl << "STEP " << format ( counter ) << " ... " << endl;
      }

      _phiOld = _phiNew;

      aol::ProgressBar<> pb( "Sweeping" );
      pb.start( 2 * 2 * 2 * ( max_vertex[0] - min_vertex[0] ) * ( max_vertex[1] - min_vertex[1] ) * ( max_vertex[2] - min_vertex[2] ), 2, 1 ); // at least approximately
      // these nested loops are used for looping forward and backward.
      for ( signed char s0 = -1; s0 <= 1; s0 += 2 ) {
        for ( signed char s1 = -1; s1 <= 1; s1 += 2 ) {
          for ( signed char s2 = -1; s2 <= 1; s2 += 2 ) {

            for ( int i = ( s0 < 0 ? max_vertex[0] - 1 : min_vertex[0] ); ( s0 < 0 ? i >= min_vertex[0] : i < max_vertex[0] ); i += s0 ) {
              for ( int j = ( s1 < 0 ? max_vertex[1] - 1 : min_vertex[1] ); ( s1 < 0 ? j >= min_vertex[1] : j < max_vertex[1] ); j += s1 ) {
                for ( int k = ( s2 < 0 ? max_vertex[2] - 1 : min_vertex[2] ); ( s2 < 0 ? k >= min_vertex[2] : k < max_vertex[2] ); k += s2, pb++ ) {

                  if ( ( ( _pIsRelevant == NULL ) || ( _pIsRelevant->get ( i, j, k ) ) ) && ! ( _isSeed.get ( i, j, k ) ) ) {
                    treatPoint ( i, j, k, min_vertex, max_vertex );
                  }

                }
              }
            }

          }
        }
      }
      pb.finish();

#ifdef DEBUG
      for ( int i = min_vertex[0]; i < max_vertex[0] ; ++i )
        for ( int j = min_vertex[1]; j < max_vertex[1] ; ++j )
          for ( int k = min_vertex[2]; k < max_vertex[2] ; ++k )
            if ( _phiNew.get ( i, j, k ) == aol::NumberTrait<DataType>::Inf )
              cerr << i << " " << j << " " << k << "contains inf. " << endl;
      cerr << endl;
#endif

      _phiOld -= _phiNew;

      const int temp_num_inf = _phiNew.numOccurence ( aol::NumberTrait<DataType>::Inf );
      if ( temp_num_inf > 0 ) {
        DataType dummy = aol::NumberTrait<DataType>::zero;
        for ( int i = 0; i < _phiNew.size(); ++i ) {
          if ( aol::isFinite ( _phiNew[i] ) ) { // distance still infinity: do not consider
            dummy += aol::Sqr ( _phiOld[i] );  // _phiOld here contains difference
          }
        }
        noninf_difference = sqrt ( dummy );

        difference = aol::NumberTrait<DataType>::Inf;
        // in this case, _phiOld contains Inf - Inf = NaN entries, hence norm is NaN, but should be Inf
      } else {
        difference = _phiOld.norm();
      }

      ++counter;

      num_inf_difference = num_inf - temp_num_inf;
      num_inf = temp_num_inf;

      if ( _verbose ) {
        aol::SimpleFormat format ( 4, 0 );
        cerr << "difference = " << aol::detailedFormat ( difference ) << ", noninf_difference = " << aol::detailedFormat ( noninf_difference ) << ", num_inf = " << format ( num_inf ) << endl;
      }

    }

    distanceField = _phiNew;
  }


  void treatPoint ( const int i, const int j, const int k, const CoordType &min_vertex, const CoordType &max_vertex ) {

    for ( signed char sx = -1; sx <= 1 ; sx += 2 ) { // loop over neighborhood
      for ( signed char sy = -1; sy <= 1 ; sy += 2 ) {
        for ( signed char sz = -1; sz <= 1; sz += 2 ) {

          if ( checkBx ( i - sx, min_vertex, max_vertex ) && checkBy ( j - sy, min_vertex, max_vertex ) && checkBz ( k - sz, min_vertex, max_vertex ) ) {

            const signed char sign_lsf = ( ( ( ( _pIsRelevant == NULL ) || ( _pIsRelevant->get ( i - sx, j, k ) ) ) && _phiNew.get ( i - sx, j, k ) < 0 ) ||
                                           ( ( ( _pIsRelevant == NULL ) || ( _pIsRelevant->get ( i, j - sy, k ) ) ) && _phiNew.get ( i, j - sy, k ) < 0 ) ||
                                           ( ( ( _pIsRelevant == NULL ) || ( _pIsRelevant->get ( i, j, k - sz ) ) ) && _phiNew.get ( i, j, k - sz ) < 0 )     ?  -1  :  1 );

            aol::Vec3<int> y, z, w;
            aol::Vec3<DataType> vec_x, vec_y, vec_z, vec_w;

            signed char s_y = 0, s_z = 0, s_w = 0;

            choose_yzw ( y, z, w, i, j, k, sx, sy, sz, s_y, s_z, s_w );

            if ( HLUpdatePossible ( y, z, w ) ) {

              //update within the cube
              updateFunction ( i, j, k, s_y, s_z, s_w, vec_x, vec_y, vec_z, vec_w, y, z, w, sign_lsf );

            }

            const DataType
              //updates calculated on 3 sides of the cube
              phi_tmp0 = updateOnSide ( sx, sy, getOrInf ( i, j - sy, k ), getOrInf ( i - sx, j, k ) ),
              phi_tmp1 = updateOnSide ( sx, sz, getOrInf ( i, j, k - sz ), getOrInf ( i - sx, j, k ) ),
              phi_tmp2 = updateOnSide ( sy, sz, getOrInf ( i, j, k - sz ), getOrInf ( i, j - sy, k ) ),
              phi_tmpA = aol::Min ( phi_tmp0, phi_tmp1, phi_tmp2 ),

              //update calculated on 3 edges of the cube
              phi_tmp3 = aol::Abs ( getOrInf ( i - sx, j, k ) ) + _h,
              phi_tmp4 = aol::Abs ( getOrInf ( i, j - sy, k ) ) + _h,
              phi_tmp5 = aol::Abs ( getOrInf ( i, j, k - sz ) ) + _h,
              phi_tmpB = aol::Min ( phi_tmp3, phi_tmp4, phi_tmp5 ),

              phi_tmp = aol::Min ( phi_tmpA, phi_tmpB );

              // the minimal phi_tmp will be saved in phiNEW
              // note that Inf is not < Inf, hence we do not update by Inf.
              if ( phi_tmp < aol::Abs ( getOrInf ( i, j, k ) ) ) {
                printIfInfinite ( "phi_tmp", phi_tmp );
                setOrError ( i, j, k, phi_tmp * sign_lsf );
              }

          } // if in interior

        } // for sz
      } // for sy
    } // for sx

  }

  void setVerboseMode ( bool mode = true ) {
    _verbose = mode;
  }

private:
  DistanceSweeper3d ( );

  DistanceSweeper3d ( const DistanceSweeper3d< DataType, LocalBitArrayType, LocalDistanceArrayType > &/*other*/ );

  DistanceSweeper3d< DataType, LocalBitArrayType, LocalDistanceArrayType >& operator= ( const DistanceSweeper3d< DataType, LocalBitArrayType, LocalDistanceArrayType > &/*other*/ );

protected:
  inline bool checkBx ( const int x, const CoordType &min_vertex, const CoordType &max_vertex ) const {
    return ( x >= min_vertex[0] && x < max_vertex[0] );
  }

  inline bool checkBy ( const int y, const CoordType &min_vertex, const CoordType &max_vertex ) const {
    return ( y >= min_vertex[1] && y < max_vertex[1] );
  }

  inline bool checkBz ( const int z, const CoordType &min_vertex, const CoordType &max_vertex ) const {
    return ( z >= min_vertex[2] && z < max_vertex[2] );
  }

#ifdef VERBOSE
  inline void printIfInfinite ( const char* const nameOfVariable, const DataType &Value ) const {
    if ( !aol::isFinite ( Value ) ) {
      cerr << nameOfVariable << " is infinite" << endl;
    }
  }
#else
  inline void printIfInfinite ( const char* const, const DataType& ) const {}
#endif

  //! get the specified value from _phiNew or return infinity if the point is not relevant
  inline DataType getOrInf ( const int i, const int j, const int k ) const {
    return ( ( ( _pIsRelevant == NULL ) || ( _pIsRelevant->get ( i, j, k ) ) ) ? _phiNew.get ( i, j, k ) : aol::NumberTrait<DataType>::Inf );
  }

  //! get the specified value from _phiNew or return infinity if the point is not relevant
  inline DataType getOrInf ( const aol::Vec3<int> &pos ) const {
    return ( ( ( _pIsRelevant == NULL ) || ( _pIsRelevant->get ( pos ) ) ) ? _phiNew.get ( aol::Vec3<short>( pos ) ) : aol::NumberTrait<DataType>::Inf );
  }

  //! set value in _phiNew or print error message if point not relevant
  inline void setOrError ( const int i, const int j, const int k, const DataType Value ) {
#ifdef VERBOSE
    if ( ( _pIsRelevant == NULL ) || ( _pIsRelevant->get ( i, j, k ) ) ) {
      _phiNew.set ( i, j, k, Value );
    }
    else {
      cerr << "Irrelevant position " << i << " " << j << " " << k << endl;
    }
#else
    _phiNew.set ( i, j, k, Value );
#endif
  }



  //! function for sorting the three vertices which are needed for the update
  void choose_yzw ( aol::Vec3<int>& y, aol::Vec3<int>& z, aol::Vec3<int>& w,
                    const int i, const int j, const int k,
                    const signed char sx, const signed char sy, const signed char sz,
                    signed char& s_y, signed char& s_z, signed char& s_w ) const {

    if ( aol::Abs ( getOrInf ( i - sx, j, k ) ) <= aol::Abs ( getOrInf ( i, j - sy, k ) ) ) {

      if ( aol::Abs ( getOrInf ( i - sx, j, k ) ) <= aol::Abs ( getOrInf ( i, j, k - sz ) ) ) {

        y[0] = i - sx;   y[1] = j;        y[2] = k;        s_y = sx;
        z[0] = i;        z[1] = j - sy;   z[2] = k;        s_z = sy;
        w[0] = i;        w[1] = j;        w[2] = k - sz;   s_w = sz;

      } else {

        y[0] = i;        y[1] = j;        y[2] = k - sz;   s_y = sz;
        z[0] = i - sx;   z[1] = j;        z[2] = k;        s_z = sx;
        w[0] = i;        w[1] = j - sy;   w[2] = k;        s_w = sy;

      }

    } else {

      if ( aol::Abs ( getOrInf ( i, j - sy, k ) ) <= aol::Abs ( getOrInf ( i, j, k - sz ) ) ) {

        y[0] = i;        y[1] = j - sy;   y[2] = k;        s_y = sy;
        z[0] = i;        z[1] = j;        z[2] = k - sz;   s_z = sz;
        w[0] = i - sx;   w[1] = j;        w[2] = k;        s_w = sx;

      } else {

        y[0] = i;        y[1] = j;        y[2] = k - sz;   s_y = sz;
        z[0] = i - sx;   z[1] = j;        z[2] = k;        s_z = sx;
        w[0] = i;        w[1] = j - sy;   w[2] = k;        s_w = sy;

      }

    }
  }


  //! function for checking the condition whether it is possible to calulate the update by "rotation"
  bool HLUpdatePossible ( const aol::Vec3<int>& y, const aol::Vec3<int>& z, const aol::Vec3<int>& w ) const {
    aol::Vec3<int> p1, p2;
    const DataType
      phi_z = aol::Abs ( getOrInf ( z ) ),
      phi_w = aol::Abs ( getOrInf ( w ) );

    if ( phi_z >= phi_w ) {
      p1 = z;
      p2 = w;
    } else {
      p1 = w;
      p2 = z;
    }

    const bool cond1 = aol::Abs ( getOrInf ( p1 ) ) - aol::Abs ( getOrInf ( y ) ) < sqrt ( 2. ) * _h;

    const DataType cond2_rhs = ( aol::Abs ( getOrInf ( y ) ) + aol::Abs ( getOrInf ( p1 ) ) ) / 2. + sqrt ( 3 * ( _h * _h / 2. - aol::Sqr ( aol::Abs ( getOrInf ( p1 ) ) - aol::Abs ( getOrInf ( y ) ) ) / 4. ) );

    return ( cond1 && aol::Abs ( getOrInf ( p2 ) ) <= cond2_rhs );

  }

  //! update function
  void updateFunction ( const int i, const int j, const int k,
                        const signed char s_y, const signed char s_z, const signed char s_w,
                        aol::Vec3<DataType>& vec_x, aol::Vec3<DataType>& vec_y, aol::Vec3<DataType>& vec_z, aol::Vec3<DataType>& vec_w,
                        const aol::Vec3<int>& y, const aol::Vec3<int>& z, const aol::Vec3<int>& w, const int sign_lsf ) {

    aol::Vec3<DataType> vec_n, vec_isp, vec_bar;

    calcCoords_yzw ( y, z, w, vec_y, vec_z, vec_w );

    calcCoords_x ( vec_y, vec_z, vec_w, vec_n, vec_x );

    calcIntersectionPoint ( vec_x, vec_y, vec_z, vec_w, vec_isp );

    if ( calcBarycentricCoords ( vec_x, vec_y, vec_z, vec_w, vec_isp, vec_bar ) == false ) {
      cerr << " error: bar. coords are wrong! \n";
    }


    if ( vec_bar[0] > 0 && vec_bar[1] > 0 && vec_bar[2] > 0 ) {             // intersection point within the triangle
      if ( vec_x[2] < aol::Abs ( getOrInf ( i, j, k ) ) ) {
        printIfInfinite ( "vec_x[2]", vec_x[2] );
        setOrError ( i, j, k, vec_x[2] * sign_lsf );
      }

    } else if ( vec_bar[0] <= 0 && vec_bar[1] > 0 && vec_bar[2] > 0 ) {     // update has to be calculated on one side of the tetrahedron

      const DataType phi_tmp = updateOnSide ( s_z, s_w, aol::Abs ( getOrInf ( w ) ), aol::Abs ( getOrInf ( z ) ) );
      if ( phi_tmp < aol::Abs ( getOrInf ( i, j, k ) ) ) {
        printIfInfinite ( "phi_tmp", phi_tmp );
        setOrError ( i, j, k, phi_tmp * sign_lsf );
      }

    } else if ( vec_bar[1] <= 0 && vec_bar[0] > 0 && vec_bar[2] > 0 ) {

      const DataType phi_tmp = updateOnSide ( s_y, s_w, aol::Abs ( getOrInf ( w ) ), aol::Abs ( getOrInf ( y ) ) );
      if ( phi_tmp < aol::Abs ( getOrInf ( i, j, k ) ) ) {
        printIfInfinite ( "phi_tmp", phi_tmp );
        setOrError ( i, j, k, phi_tmp * sign_lsf );
      }

    } else if ( vec_bar[2] <= 0 && vec_bar[0] > 0 && vec_bar[1] > 0 ) {

      const DataType phi_tmp = updateOnSide ( s_y, s_z, aol::Abs ( getOrInf ( z ) ), aol::Abs ( getOrInf ( y ) ) );
      if ( phi_tmp < aol::Abs ( getOrInf ( i, j, k ) ) ) {
        printIfInfinite ( "phi_tmp", phi_tmp );
        setOrError ( i, j, k, phi_tmp * sign_lsf );
      }

    } else if ( vec_bar[0] <= 0 && vec_bar[1] <= 0 && vec_bar[2] > 0 ) {      // update has to be calculated on one edge of the tetrahedron

      if ( vec_w[2] + _h < aol::Abs ( getOrInf ( i, j, k ) ) ) {
        printIfInfinite ( "vec_w[2] + _h", vec_w[2] + _h );
        setOrError ( i, j, k, ( vec_w[2] + _h ) * sign_lsf );
      }

    } else if ( vec_bar[1] <= 0 && vec_bar[2] <= 0 && vec_bar[0] > 0 ) {

      if ( vec_y[2] + _h < aol::Abs ( getOrInf ( i, j, k ) ) ) {
        printIfInfinite ( "vec_y[2] + _h", vec_y[2] + _h );
        setOrError ( i, j, k, ( vec_y[2] + _h ) * sign_lsf );
      }

    } else if ( vec_bar[0] <= 0 && vec_bar[2] <= 0 && vec_bar[1] > 0 ) {

      if ( vec_z[2] + _h < aol::Abs ( getOrInf ( i, j, k ) ) ) {
        printIfInfinite ( "vec_z[2] + _h", vec_z[2] + _h );
        setOrError ( i, j, k, ( vec_z[2] + _h ) * sign_lsf );
      }
    }
  }


  //! function for calculating the coordinates of y,z and w
  void calcCoords_yzw ( const aol::Vec3<int>& y, const aol::Vec3<int>& z, const aol::Vec3<int>& w, aol::Vec3<DataType>& vec_y, aol::Vec3<DataType>& vec_z, aol::Vec3<DataType>& vec_w ) const {

    vec_y[0] = 0;
    vec_y[1] = 0;
    vec_y[2] = aol::Abs ( getOrInf ( y ) );

    vec_w[0] = 0;
    vec_w[1] = sqrt ( aol::Abs ( 2 * _h * _h - aol::Sqr ( aol::Abs ( getOrInf ( w ) ) - aol::Abs ( getOrInf ( y ) ) ) ) );
    vec_w[2] = aol::Abs ( getOrInf ( w ) );

    vec_z[0] = sqrt ( aol::Abs ( 3. / 2. * _h * _h - aol::Sqr ( aol::Abs ( getOrInf ( z ) ) - aol::Abs ( getOrInf ( y ) ) ) ) );
    vec_z[1] = vec_w[1] / 2. ;            // this depends on vec_w, so calculate z later
    vec_z[2] = aol::Abs ( getOrInf ( z ) );

    // TODO: static_casts may be necessary here!

  }

  //! function for calculating the coordinates of x which will be updated later
  void calcCoords_x ( const aol::Vec3<DataType>& vec_y, const aol::Vec3<DataType>& vec_z, const aol::Vec3<DataType>& vec_w, aol::Vec3<DataType>& vec_n, aol::Vec3<DataType>& vec_x ) const {
    const aol::Vec3<DataType> vec_zy = vec_z - vec_y;
    const aol::Vec3<DataType> vec_wy = vec_w - vec_y;

    vec_n.makeCrossProduct ( vec_zy, vec_wy );
    vec_n.normalize();

    const aol::Vec3<DataType> vec_m_yw = vec_y + ( vec_wy / 2 );
    const aol::Vec3<DataType> vec_m    = vec_z + static_cast<DataType> ( 2. / 3. ) * ( vec_m_yw - vec_z );

    vec_x = vec_m + static_cast<DataType> ( _h / sqrt ( 3.0 ) ) * vec_n;

  }

  //! function for calculating the intersection point of normal and plane yzw
  void calcIntersectionPoint ( const aol::Vec3<DataType>& vec_x, const aol::Vec3<DataType>& vec_y, const aol::Vec3<DataType>& vec_z, const aol::Vec3<DataType>& vec_w, aol::Vec3<DataType>& vec_isp ) const {

    aol::Matrix33<DataType> mat_isp;
    mat_isp.setCol ( 0, aol::Vec3<DataType> ( 0, 0, 1 ) );
    mat_isp.setCol ( 1, vec_z - vec_y );
    mat_isp.setCol ( 2, vec_w - vec_y );

    aol::Vec3<DataType> vec_par;

    aol::Matrix33<DataType> QRI = mat_isp.inverse();

    QRI.mult ( vec_y - vec_x , vec_par );

    vec_isp = vec_x;
    vec_isp[2] += vec_par[0];

    // now corresponds to handwritten notes but differs from previous code!
  }

  //! function for calculating the barycentric coordinates
  bool calcBarycentricCoords ( const aol::Vec3<DataType>& vec_x, const aol::Vec3<DataType>& vec_y, const aol::Vec3<DataType>& vec_z, const aol::Vec3<DataType>& vec_w, aol::Vec3<DataType> &vec_isp, aol::Vec3<DataType>& vec_bar ) const {

    aol::Matrix33<DataType> mat_bar;
    mat_bar.setCol ( 0, vec_y - vec_x );
    mat_bar.setCol ( 1, vec_z - vec_x );
    mat_bar.setCol ( 2, vec_w - vec_x );

    vec_isp -= vec_x;

    aol::Matrix33<DataType> QRI = mat_bar.inverse();

    QRI.mult ( vec_isp, vec_bar );

    return ( 1 - _eps < vec_bar[0] + vec_bar[1] + vec_bar[2]  &&  vec_bar[0] + vec_bar[1] + vec_bar[2] < 1 + _eps );

  }


  //! function for calculating the update of a vertex on one side of the cube
  DataType updateOnSide ( const signed char sx, const signed char sy, const DataType phi_i_sy, const DataType phi_sx_j ) {

    DataType theta = 0, phi_tmp = 0;

    const DataType m = sx * sy * ( aol::Abs ( phi_i_sy ) - aol::Abs ( phi_sx_j ) ) / _h;

    if ( phi_i_sy == aol::NumberTrait<DataType>::Inf && phi_sx_j  == aol::NumberTrait<DataType>::Inf )
      return ( aol::NumberTrait<DataType>::Inf );

    if ( ( sx*sy != 0 ) && ( -1 < m ) && ( m < 1 ) ) {

      const DataType
      dummy1 = -sx * sy + m * sqrt ( 2 - m * m ),
               dummy2 = -sx * sy - m * sqrt ( 2 - m * m ),
                        denom = m * m - sx * sx;

      if ( aol::Abs ( m ) > 1.e-15 && aol::Abs ( denom ) > _smallValue  ) {

        const DataType
          frac1 = static_cast<DataType> ( ( 1.0 * dummy1 ) / denom ),
          frac2 = static_cast<DataType> ( ( 1.0 * dummy2 ) / denom ),

          theta1 = atan ( frac1 ) + ( 1.0 - sx ) * aol::NumberTrait<DataType>::pi / 2,
          theta2 = atan ( frac2 ) + ( 1.0 - sx ) * aol::NumberTrait<DataType>::pi / 2,

          T1 = -sx * sin ( theta1 ) + cos ( theta1 ) * sy,

          test =  m - T1;

        if ( test < _smallValue ) {
          theta = theta1;
        } else {
          theta = theta2;
        }

      }

      if ( aol::Abs ( denom ) < _smallValue  ) {
        theta = ( 1 - sx ) * aol::NumberTrait<long double>::pi * 0.5;
      }

      if ( aol::Abs ( m ) < _smallValue ) {
        const DataType
          fraction = static_cast<DataType> ( ( 1.0 * sy ) / sx ),
          sum_1 = atan ( fraction ),
          sum_2 = ( 1 - sx ) * aol::NumberTrait<long double>::pi * 0.5;
        theta = sum_1 + sum_2;
      }

      phi_tmp = ( sx * cos ( theta ) * std::abs ( phi_sx_j ) + sy * sin ( theta ) * std::abs ( phi_i_sy ) + _h ) / ( aol::Abs ( cos ( theta ) ) +  aol::Abs ( sin ( theta ) )  );

      return phi_tmp;

    } else {

      return ( aol::NumberTrait<DataType>::Inf );

    }

  }

};

}
#endif
