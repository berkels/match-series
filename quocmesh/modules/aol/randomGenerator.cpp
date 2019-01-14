#include <randomGenerator.h>

void aol::RandomGenerator::randomize ( ) {
  time_t wct;
  unsigned short wct_millisec;
  aol::getWallClockTime ( wct, wct_millisec );
  _seed = static_cast<unsigned int> ( wct ) ^ wct_millisec ;
  initializeY();
}


void aol::RandomGenerator::initializeY ( ) {
  if ( sizeof ( unsigned int ) != 4 )
    throw aol::Exception ( "Cannot use aol::RandomGenerator if sizeof( unsigned int ) != 4", __FILE__, __LINE__ );

  _y[0] = _seed ^ 0xffffffff;
  for ( unsigned int i = 1; i < 624; ++i ) {
    _y[i] = ( 69069 * _y[i-1] ) & 0xffffffff;
  }

  // 7 times refilling the array of random numbers by the Mersenne twister should be sufficient to remove any artifacts by the initialization
  for ( unsigned int i = 0; i < 7 * 624; ++i ) {
    getNextRnd();
  }
}

unsigned int aol::RandomGenerator::getNextRnd() {
  static const unsigned int HI = 0x80000000, LO = 0x7fffffff;
  static const unsigned int A[2] = { 0, 0x9908b0df };
  static const unsigned int M = 397, N = 624;

  unsigned int e;

#ifdef _OPENMP
#pragma omp critical ( aol_RandomGenerator_getNextRnd )
#endif
  {
    if ( _index == N ) {
      unsigned int h;
      for ( unsigned int k = 0 ; k < N - M ; ++k ) {
        h = ( _y[k] & HI ) | ( _y[k+1] & LO );
        _y[k] = _y[k+M] ^ ( h >> 1 ) ^ A[h & 1];
      }
      for ( unsigned int k = N - M ; k < N - 1 ; ++k ) {
        h = ( _y[k] & HI ) | ( _y[k+1] & LO );
        _y[k] = _y[k+ ( M-N ) ] ^ ( h >> 1 ) ^ A[h & 1];
      }
      h = ( _y[N-1] & HI ) | ( _y[0] & LO );
      _y[N-1] = _y[M-1] ^ ( h >> 1 ) ^ A[h & 1];
      _index = 0;
    }

    e = _y[_index++];
    // tempering:
    e ^= ( e >> 11 );
    e ^= ( e << 7 ) & 0x9d2c5680;
    e ^= ( e << 15 ) & 0xefc60000;
    e ^= ( e >> 18 );
  }
  return ( e );
}



