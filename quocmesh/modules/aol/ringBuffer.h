#ifndef __RINGBUFFER_H
#define __RINGBUFFER_H

// use aol.h to include stl vector definition:
#include <aol.h>

namespace aol {

/****************************************************************************
 *
 *        CLASS RingBuffer
 */
/**
 *  \brief Container that can overwrite its oldest item when
 *         inserting a new one.
 *
 *  It is internally designed as a container with additional knowledge
 *  about which item is the first. Pushing back one item to the front
 *  means overwriting this item, pushing to the back of the list
 *  overwrites the item behind the first one. Additionally, as the ring
 *  buffer is currently internally implemented as an std::vector, random
 *  access via operator[] is provided.
 *
 *  The used may choose for himself if the indexing runs from oldest
 *  to newest or vice versa (or even mixed) by pushing to front or
 *  back. Derived classes may decide and provide additional access
 *  to a well-defined "newest" and "oldest" item.
 *
 *  Resizing is done via resize() or successive growing of the
 *  ring buffer. Note that only building up the entire list with
 *  grow_back() asserts that every enlargement is O(1). In general,
 *  both grow_front() and grow_back() are O(n), and even starting
 *  from an empty list and constructing via grow_front() takes
 *  in summa O(n^2/2) operations.
 *
 *  Using only the RingBuffer standard constructor and grow_back()
 *  or grow_front(), only copy constructor and assignment operator
 *  of T are needed. Constructing a RingBuffer of given size or
 *  resizing also needs a default constructor for T.
 *
 *  \author von Deylen
 */
template <typename T>
class RingBuffer {
public:
  RingBuffer();
  explicit RingBuffer ( int len );

  // access functions
  T & operator[] ( int i );
  const T & operator[] ( int i ) const;

  // read-only functions
  int size() const;

  T sum() const;
  T getMaxValue() const;
  T arithmeticMean() const;

  // setting functions
  void push_front ( const T & entry );
  void push_back ( const T & entry );
  void grow_front ( const T & entry );
  void grow_back ( const T & entry );

  void setAll ( const T & entry );

  void shift ( int offset );

  //! CAVE: all entries are invalidated after a
  //! call of resize().
  void resize ( int len );
  void clear();

protected:
  int map_index ( int index_outside ) const;

private:
  vector<T> _data;
  int       _start;
};

/**** IMPLEMENTATION ********************************************************
 *
 *
 *
 */

//---------------------------------------------------------------------------

template <typename T>
RingBuffer<T>::RingBuffer() :
    _start ( 0 ) {}

//---------------------------------------------------------------------------

template <typename T>
RingBuffer<T>::RingBuffer ( int len ) :
    _data  ( len ),
    _start ( 0 ) {}

//---------------------------------------------------------------------------
template <typename T>
const T & RingBuffer<T>::operator[] ( int i ) const
                                      {   return _data[map_index ( i ) ];   }
//---------------------------------------------------------------------------
template <typename T>
T & RingBuffer<T>::operator[] ( int i ) {  return _data[map_index ( i ) ];  }
//---------------------------------------------------------------------------

template <typename T>
void RingBuffer<T>::push_front ( const T & entry ) {
  int last_index = static_cast<int> ( _data.size() ) - 1;
  ( *this ) [last_index] = entry;
  _start--;

  // unfortunately, the sign of (neg) % (pos) is
  // implementation specific, so we cannot write
  // start %= size() but have to implement
  // a conditioned addition by hand.
  if ( _start < 0 )
    _start += size();
}

//---------------------------------------------------------------------------

template <typename T>
void RingBuffer<T>::push_back ( const T & entry ) {
  ( *this ) [0] = entry;
  _start++;
  _start %= size();
}

//---------------------------------------------------------------------------

template <typename T>
void RingBuffer<T>::grow_front ( const T & entry ) {
  // use push_back to avoid need of a standard
  // constructor:
  _data.push_back ( entry );

  // if _start points to the very beginning
  // of _data, we simply can change it to the
  // back of _data and move no elements at all
  if ( _start == 0 ) {
    _start = size() - 1;
    return;
  }

  // move all elements behind _start:
  for ( int i = size() - 1; i > _start; --i )
    _data[i] = _data[i - 1];

  // write given entry to its correct position
  // (front of buffer)
  _data[_start] = entry;
}

//---------------------------------------------------------------------------

template <typename T>
void RingBuffer<T>::grow_back ( const T & entry ) {
  // use push_back to avoid need of a standard
  // constructor:
  _data.push_back ( entry );

  // if _start points to the very beginning
  // of _data, inserting at the back of _data
  // was just where entry should be.
  if ( _start == 0 )
    return;

  // move all elements behind _start:
  for ( int i = _start; i < size() - 1; ++i )
    _data[i + 1] = _data[i];

  // write given entry to its correct position
  // (front of buffer) and set _start just behind
  _data[_start] = entry;
  _start++;
}

//---------------------------------------------------------------------------
template <typename T>
void RingBuffer<T>::shift ( int offset )   {   _start += offset & size();   }
//---------------------------------------------------------------------------
template <typename T>
int RingBuffer<T>::size() const
                            {   return static_cast<int> ( _data.size() );   }
//---------------------------------------------------------------------------
template <typename T>
void RingBuffer<T>::clear()                             {   resize ( 0 );   }
//---------------------------------------------------------------------------

template <typename T>
void RingBuffer<T>::resize ( int len ) {
  _data.resize ( len );
  _start = 0;
}

//---------------------------------------------------------------------------

template <typename T>
void RingBuffer<T>::setAll ( const T & entry ) {
  for ( int i = 0; i < size(); ++i )
    ( *this ) [i] = entry;
}

//---------------------------------------------------------------------------

template <typename T>
T RingBuffer<T>::sum() const {
  T ret = aol::ZOTrait<T>::zero;
  for ( size_t i = 0; i < _data.size(); ++i )
    ret += _data[i];
  return ret;
}

//---------------------------------------------------------------------------

template <typename T>
T RingBuffer<T>::getMaxValue() const {
  QUOC_ASSERT ( size() > 0 );
  T ret = _data[0];
  for ( size_t i = 1; i < _data.size(); ++i )
    ret = aol::Max( _data[i], ret );
  return ret;
}

//---------------------------------------------------------------------------

template <typename T>
T RingBuffer<T>::arithmeticMean() const {
  return sum() / size();
}

//---------------------------------------------------------------------------

template <typename T>
int RingBuffer<T>::map_index ( int index_outside ) const {
  return ( _start + index_outside ) % size();
}

//---------------------------------------------------------------------------

} // end of namespace aol.

#endif
