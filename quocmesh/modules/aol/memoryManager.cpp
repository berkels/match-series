#include <memoryManager.h>

void aol::MemoryManager::deleteUnlocked ( int NumToRetain ) {
#ifdef DO_NOT_USE_MEMORYMANAGER
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( NumToRetain );
#ifdef VERBOSE
  cerr << "aol::MemoryManager is turned off, deleteUnlocked has no effect" << endl;
#endif

#else

#ifdef _OPENMP
#pragma omp critical ( aol_MemoryManager )
#endif
  {
    list< MMStore >::iterator it;
    const list< MMStore >::iterator end = _storage.end ();

#ifdef VERBOSE
    cerr << "aol::MemoryManager: deleting from " << _numStored << " unlocked blocks, retaining " << NumToRetain << endl;
#endif

    it = _storage.begin();
    for ( int retained = 0; ( retained < NumToRetain ) && ( it != end ); ++retained ) {
      ++it;
    }

    for ( ; it != end; it = _storage.erase ( it ) ) {
#ifdef USE_SSE
      aol::aligned_memory_deallocation ( (*it).pBlock );
#else
#ifdef USE_DUMA
      free ( (*it).pBlock );
#else
      operator delete ( (*it).pBlock );
#endif
#endif
      --_numStored;
      _memusage -= ( (*it).size );
    }
#ifdef VERBOSE
    cerr << "aol::MemoryManager: after deleteUnlocked: MemoryManager stores " << _memusage / ( 1024*1024 ) << " MiB in " << _numStored << " blocks." << endl;;
#endif
  }

#endif
}


int64_t aol::MemoryManager::memoryManagerMemoryUsage () {
  return _memusage;
}


void aol::MemoryManager::setMaxRetain ( const int MaxRetain ) {
#ifdef DO_NOT_USE_MEMORYMANAGER
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( MaxRetain );
#ifdef VERBOSE
  cerr << "aol::MemoryManager is turned off, setMaxRetain has no effect" << endl;
#endif
#else
#ifdef _OPENMP
#pragma omp critical ( aol_MemoryManager )
#endif
  {
    _maxRetain = MaxRetain;
  }
#endif
}


void aol::MemoryManager::setMemusageLimit ( const int64_t MaxMemusage ) {
#ifdef DO_NOT_USE_MEMORYMANAGER
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( MaxMemusage );
#ifdef VERBOSE
  cerr << "aol::MemoryManager is turned off, setMemusageLimit has no effect" << endl;
#endif
#else
#ifdef _OPENMP
#pragma omp critical ( aol_MemoryManager )
#endif
  {
    _memusageLimit = MaxMemusage;
  }
#endif
}


void* aol::MemoryManager::allocateAtLeast ( int& Length, const size_t PointeeSize ) {
  // If the user wants a memory block of length 0, just return a NULL pointer. Such a pointer
  // can't be used for anything anyway and this saves us from doing more checks further below
  // since allocating a memory block of lenght 0 will give a NULL pointer under some platforms.
  if ( Length == 0 )
    return NULL;

#ifdef BOUNDS_CHECK
  if ( Length < 0 )
    throw aol::Exception ( "aol::MemoryManager: Cannot allocate negative amount of memory", __FILE__, __LINE__ );
#endif

#ifndef DO_NOT_USE_MEMORYMANAGER
  void* ptr = NULL;
#ifdef _OPENMP
#pragma omp critical ( aol_MemoryManager )
#endif
  {
    int64_t maxlen = ( Length * PointeeSize ) << 1;

    // Find vector of correct length
    const list<MMStore>::iterator end = _storage.end ();

    for ( list<MMStore>::iterator it = _storage.begin (); it != end; ++it ) {
#ifdef VERBOSE
      cerr << "aol::MemoryManager: searching, found length " << (*it).size << endl;
#endif

      if ( ( (*it).size >= static_cast<int64_t> ( Length * PointeeSize ) ) && ( (*it).size < maxlen ) && ( (*it).size % PointeeSize == 0 ) ) {
        Length = (*it).size / PointeeSize; // integer division on purpose, must be multiple
#ifdef VERBOSE
        cerr << "aol::MemoryManager: recycling memory with " << ( Length * PointeeSize ) << " bytes at address " << (*it).pBlock << " with length " << (*it).size << endl;
#endif

        ptr = (*it).pBlock;
        --_numStored;
        _memusage -= Length * PointeeSize;
        _storage.erase ( it );
        break;
      }
    }
  }
  if ( ptr != NULL )
    return ptr;
  // Otherwise allocate new data

#endif //DO_NOT_USE_MEMORYMANAGER

#ifdef USE_SSE
  void *new_vec = aol::aligned_memory_allocation ( Length * PointeeSize, 16 );
#else
#ifdef USE_DUMA
  void *new_vec = malloc ( Length * PointeeSize );
#else
  void *new_vec = operator new ( Length * PointeeSize );
#endif
#endif // USE_SSE

#ifndef DO_NOT_USE_MEMORYMANAGER
#ifdef VERBOSE
  cerr << "aol::MemoryManager: allocating new memory at address " << new_vec << " with length " << Length * PointeeSize << endl;
#endif
#endif

  if ( !new_vec )
    throw aol::OutOfMemoryException ( "aol::MemoryManager: Allocate not successful.", __FILE__, __LINE__ );

  return new_vec;
}


void aol::MemoryManager::deallocate ( void *Ptr, const int Length, const size_t PointeeSize ) {
#ifdef DO_NOT_USE_MEMORYMANAGER
  // no-vectormanager case:
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( Length + PointeeSize );
#ifdef USE_SSE
  aol::aligned_memory_deallocation ( Ptr );
#else
#ifdef USE_DUMA
  free( Ptr );
#else
  operator delete ( Ptr );
#endif
#endif

#else

#ifdef USE_SSE
  if ( ( reinterpret_cast<const uintptr_t> ( Ptr ) % 16 ) != 0 ) {
    throw Exception ( "aol::MemoryManager::deallocate: With USE_SSE, the data pointer must be 16 byte aligned.\n", __FILE__, __LINE__ );
  }
#endif

#ifdef VERBOSE
  cerr << "aol::MemoryManager: unlocking memory at address " << Ptr << " with length " << Length * PointeeSize << endl;
#endif

#ifdef _OPENMP
#pragma omp critical ( aol_MemoryManager )
#endif
  {
    _storage.push_front ( MMStore ( Length * PointeeSize, Ptr ) );
    ++_numStored;
    _memusage += Length * PointeeSize;

    while ( ( _numStored > _maxRetain ) || ( _memusage > _memusageLimit ) ) {
#ifdef VERBOSE
      cerr << "aol::MemoryManager: removing memory of length " << ( _storage.back() ).size << " bytes" << endl;
#endif

#ifdef USE_SSE
      aol::aligned_memory_deallocation ( (_storage.back()).pBlock );
#else
#ifdef USE_DUMA
      free ( (_storage.back()).pBlock );
#else
      operator delete ( (_storage.back()).pBlock );
#endif
#endif

      --_numStored;
      _memusage -= ( (_storage.back()).size );
      _storage.pop_back ();

#ifdef VERBOSE
      cerr << "aol::MemoryManager: was too large, deleted last entry" << endl;
#endif

    }

#ifdef VERBOSE
    cerr << "aol::MemoryManager stores " << _memusage / ( 1024*1024 ) << " MiB in " << _numStored << " blocks." << endl;
#endif
  }

#endif
}


list< aol::MemoryManager::MMStore > aol::MemoryManager::_storage = std::list< MMStore >();
int aol::MemoryManager::_numStored = 0;
int64_t aol::MemoryManager::_memusage = 0;
int aol::MemoryManager::_maxRetain = 256;
int64_t aol::MemoryManager::_memusageLimit = 512 * ( 1 << 20 ); // 512 MiB


#ifndef DO_NOT_USE_MEMORYMANAGER
namespace {
  static const int DUMMY ( atexit ( aol::MemoryManager::clearOnExit ) );
}
#endif
