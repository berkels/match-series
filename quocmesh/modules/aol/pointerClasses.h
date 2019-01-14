#ifndef __POINTERCLASSES_H
#define __POINTERCLASSES_H

#include <aol.h>

namespace aol {

/**
 *  \brief
 *  This class lets you choose whether it is responsible for
 *  destroying the held object or not.
 *
 *  When creating copies via copy constructor or assignment,
 *  the rhs operand will never be responsible for destroying
 *  the object passed from the lhs operand.
 *
 *  \author von Deylen
 *          Created 2008-04-22
 */
template <typename T>
class DeleteFlagPointer {

public:
  explicit DeleteFlagPointer( T* p = NULL, bool deleteFlag = false);
  ~DeleteFlagPointer();

        T * operator -> ();
  const T * operator -> () const;

        T & operator* ();
  const T & operator* () const;

  T * get() const;
  void reset(T * p = NULL, bool deleteFlag = false);

  bool getDeleteFlag() const;
  void setDeleteFlag ( bool deleteFlag );

private:
  DeleteFlagPointer (const DeleteFlagPointer<T>& alius);
  DeleteFlagPointer<T>& operator= (const DeleteFlagPointer<T> & alius);

  T * _pointer;
  bool _deleteFlag;
};

//---------------------------------------------------------------------------

template <typename T> DeleteFlagPointer<T>::
DeleteFlagPointer( T* p, bool deleteFlag )
  : _pointer ( p )
  , _deleteFlag ( deleteFlag ) {
}

//---------------------------------------------------------------------------

template <typename T> DeleteFlagPointer<T>::
~DeleteFlagPointer() {
  reset();
}

//---------------------------------------------------------------------------

template <typename T> DeleteFlagPointer<T>::
DeleteFlagPointer (const DeleteFlagPointer<T>& alius) {
  *this = alius;
}

//---------------------------------------------------------------------------

template <typename T> DeleteFlagPointer<T>& DeleteFlagPointer<T>::
operator= (const DeleteFlagPointer<T> & alius) {
  reset ( alius.get(), false );
  return *this;
}

//---------------------------------------------------------------------------
template <typename T> T * DeleteFlagPointer<T>::
operator -> ()                                       {   return _pointer;   }
//---------------------------------------------------------------------------
template <typename T> const T * DeleteFlagPointer<T>::
operator -> () const                                 {   return _pointer;   }
//---------------------------------------------------------------------------
template <typename T> T & DeleteFlagPointer<T>::
operator*  ()                                       {   return *_pointer;   }
//---------------------------------------------------------------------------
template <typename T> const T & DeleteFlagPointer<T>::
operator*  () const                                 {   return *_pointer;   }
//---------------------------------------------------------------------------
template <typename T> T * DeleteFlagPointer<T>::
get() const                                          {   return _pointer;   }
//---------------------------------------------------------------------------
template <typename T> bool DeleteFlagPointer<T>::
getDeleteFlag() const                             {   return _deleteFlag;   }
//---------------------------------------------------------------------------
template <typename T> void DeleteFlagPointer<T>::
setDeleteFlag ( bool deleteFlag )           {   _deleteFlag = deleteFlag;   }
//---------------------------------------------------------------------------

template <typename T> void DeleteFlagPointer<T>::
reset(T * p, bool deleteFlag ) {
  if (p != get() )
    if ( getDeleteFlag() )
      delete get();
  _pointer = p;
  _deleteFlag = deleteFlag;
}

//---------------------------------------------------------------------------

} // end of namespace aol.

#endif
