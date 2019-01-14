#ifndef __MEMORYMANAGER_H
#define __MEMORYMANAGER_H

#include <aol.h>

namespace aol {

/** Class for allocating and deallocating memory (16-byte aligned in
 *  case SSE is used).  A certain amount of memory is kept here rather
 *  than returned to the system to improve performance (unless
 *  explicitly disabled by definig DO_NOT_USE_MEMORYMANAGER).  Newest
 *  blocks of memory deallocated are recycled first, oldest blocks are
 *  dropped first.
 *
 *  \todo consistent use of size_t and int data types
 *
 *  \author Schwen (MEVIS), based on older code
 */
class MemoryManager {
public:
  struct MMStore {
    int64_t size;
    void* pBlock;
    MMStore ( const int64_t Size, void* PBlock ) : size ( Size ), pBlock ( PBlock ) { };
  };

private:
  static list< MMStore > _storage;

  static int _numStored;         //!< number of memory blocks being stored
  static int64_t _memusage;      //!< memory being used
  static int _maxRetain;         //!< Maximum number of memory blocks to retain in manager
  static int64_t _memusageLimit; //!< maximum memory usage

public:
  //! Return how much memory is used by the MemoryManager
  static int64_t memoryManagerMemoryUsage ();

  //! Return a pointer to memory of at least size Length * PointeeSize, write actual length to Length
  //! \warning The user has to store the actual length returned and pass this value to deallocate later, otherwise memory will get lost.
  static void* allocateAtLeast ( int& Length, const size_t PointeeSize );

  //! Vector manager, store for reuse
  static void deallocate ( void *Ptr, const int Length, const size_t PointeeSize );

  //! Clean reuse storage
  //! \todo rename method
  static void deleteUnlocked ( const int NumToRetain = 0 );

  //! set maximum number of vectors to be kept, does not affect current size
  static void setMaxRetain ( const int MaxRetain );

  //! set maximum amount of memory to be kept, does not affect current size
  static void setMemusageLimit ( const int64_t MaxMemusage );

  //! to free used memory on exit to prevent memory leak, should not be called directly except by atexit
  static void clearOnExit ( ) {
    deleteUnlocked();
  }
};

}

#endif
