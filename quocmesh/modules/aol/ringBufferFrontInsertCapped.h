#ifndef __RINGBUFFERFRONTINSERTCAPPED_H
#define __RINGBUFFERFRONTINSERTCAPPED_H

#include <ringBuffer.h>

namespace aol {

/****************************************************************************
 *
 *        CLASS RingBufferFrontInsertCapped
 */
/**
 *  \brief Ring buffer that inserts in front, start with size zero
 *         and grows up to its initially given maximal size.
 *
 *  See documentation of RingBuffer for details.
 *
 *  \author von Deylen
 */

template <typename T>
class RingBufferFrontInsertCapped {
public:
  RingBufferFrontInsertCapped ( int maxLen );

  // access functions
  T & operator[] ( int i );
  const T & operator[] ( int i ) const;

  // read-only functions
  int size() const;
  int maxSize() const;

  T sum() const;
  T arithmeticMean() const;
  T getMaxValue() const;

  // set function
  void push ( const T & entry );
  void setAll ( const T & entry );
  void clear();

protected:
  int _maxLen;

private:
  RingBuffer<T> _buffer;
};

/**** IMPLEMENTATION ********************************************************
 *
 *
 *
 */

//---------------------------------------------------------------------------

template <typename T>
RingBufferFrontInsertCapped<T>::RingBufferFrontInsertCapped ( int maxLen )
    : _maxLen ( maxLen ) {}

//---------------------------------------------------------------------------
template <typename T>
T & RingBufferFrontInsertCapped<T>::operator[] ( int i )
                                                   {   return _buffer[i];   }
//---------------------------------------------------------------------------
template <typename T>
const T & RingBufferFrontInsertCapped<T>::operator[] ( int i ) const
                                                   {   return _buffer[i];   }
//---------------------------------------------------------------------------
template <typename T>
int RingBufferFrontInsertCapped<T>::size() const {  return _buffer.size();  }
//---------------------------------------------------------------------------
template <typename T>
int RingBufferFrontInsertCapped<T>::maxSize() const   {   return _maxLen;   }
//---------------------------------------------------------------------------
template <typename T>
T RingBufferFrontInsertCapped<T>::sum() const   {   return _buffer.sum();   }
//---------------------------------------------------------------------------
template <typename T>
T RingBufferFrontInsertCapped<T>::arithmeticMean() const
                                     {   return _buffer.arithmeticMean();   }
//---------------------------------------------------------------------------
template <typename T>
T RingBufferFrontInsertCapped<T>::getMaxValue() const {
  return _buffer.getMaxValue();
}

//---------------------------------------------------------------------------

template <typename T>
void RingBufferFrontInsertCapped<T>::push ( const T & entry ) {
  if ( size() < _maxLen )
    _buffer.grow_front ( entry );
  else
    _buffer.push_front ( entry );
}

//---------------------------------------------------------------------------
template <typename T>
void RingBufferFrontInsertCapped<T>::setAll ( const T & entry )
                                            {   _buffer.setAll ( entry );   }
//---------------------------------------------------------------------------
template <typename T>
void RingBufferFrontInsertCapped<T>::clear()         {   _buffer.clear();   }
//---------------------------------------------------------------------------

} // end of namespace aol.

#endif
