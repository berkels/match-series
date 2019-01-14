#ifndef __QMHEAP_H
#define __QMHEAP_H

#include <array.h> // for IndexedHeap

namespace qc {

// uses Barton and Nackman trick

template <class T, class T_leaftype>
class HeapBase {
public:
  int append ( const T &Element ) {
    this->data.push_back ( Element );
    return static_cast<int> ( this->data.size() ) - 1;
  }

  void push ( const T &Element ) {
    asLeaf()->push ( Element );
  }

  void pop ( T &Element ) {
    asLeaf()->pop ( Element );
  }

  int isEmpty( ) const {
    return ( data.size( ) == 0 );
  }

  int upHeap ( int Index ) {
    while ( Index > 0  && ( data[ ( Index - 1 ) >> 1 ] > data[ Index ] ) ) {
      swap ( ( Index - 1 ) >> 1, Index );
      Index >>= 1;
    }
    return Index;
  }

  void downHeap ( int Index ) {
    if ( data.size( ) < 2 ) return;

    while ( firstChild ( Index ) + 1 <  static_cast<int> ( data.size( ) ) &&
            ( data[ Index ] > data[ firstChild ( Index ) ]
              || data[ Index ] > data[ firstChild ( Index ) + 1 ] ) ) {
      if ( data[ firstChild ( Index ) ] < data[ firstChild ( Index ) + 1] ) {
        swap ( Index, firstChild ( Index ) );
        Index = firstChild ( Index );
      } else  {
        swap ( Index, firstChild ( Index ) + 1 );
        Index = firstChild ( Index ) + 1;
      }
    }
    if ( firstChild ( Index ) + 1 ==  static_cast<int> ( data.size( ) )
         && data[ Index ] > data[ firstChild ( Index ) ] ) {
      swap ( Index, firstChild ( Index )  );
    }
  }

  void touch ( int Index ) {
    downHeap ( upHeap ( Index ) );
  }

  T &operator[] ( const int Index ) { return data[ Index ]; }

  const T &operator[] ( const int Index ) const { return data[ Index ]; }

  typename vector<T>::iterator begin( ) { return data.begin( ); }

  typename vector<T>::iterator end( ) { return data.end( ); }

  void erase ( const typename vector<T>::iterator &first, const typename vector<T>::iterator &last ) {
    data.erase ( first, last );
  }

  void erase( ) {
    data.erase ( data.begin( ), data.end( ) );
  }

  int size( ) const  { return static_cast<int>(data.size( )); }

protected:
  vector<T> data;

  inline int firstChild ( int Index ) {
    return ( Index << 1 ) + 1;
  }

  void swap_def ( int I1, int I2 ) {
    T tmp = data[ I1 ];
    data[ I1 ] = data[ I2 ];
    data[ I2 ] = tmp;
  }

  void swap ( int I1, int I2 ) {
    asLeaf()->swap ( I1, I2 );
  }

  T_leaftype* asLeaf( ) {
    return static_cast<T_leaftype*> ( this );
  }
};


template <class T>
class Heap : public HeapBase<T, Heap<T> > {
public:


  void push ( const T &Element ) {
    this->data.push_back ( Element );
    this->upHeap ( static_cast<int> ( this->data.size( ) ) - 1 );
  }

  void pop ( T &Element ) {
    if ( !this->isEmpty( ) ) {
      Element = this->data[ 0 ];
      this->data[ 0 ] = this->data.back( );
      this->data.pop_back( );
      this->downHeap ( 0 );
    }
  }

  void swap ( int I1, int I2 ) {
    this->swap_def ( I1, I2 );
  }
};


template <class T>
class IndexedHeap : public HeapBase<T, IndexedHeap<T> > {
public:
  IndexedHeap( ) :
      indexField ( NULL ) {}

  void setIndexField ( Array<int> *IndexField ) {
    indexField = IndexField;
  }

  void peekMin ( T &Element ) {
    if ( !this->isEmpty( ) ) {
      Element = this->data[ 0 ];
    }
  }

  void dump( ) {
    for ( int i = 0; i < this->size( ); i++ ) {
      this->data[ i ].dump( );
    }
  }

  void append ( const T &El ) {
    this->data.push_back ( El );
    if ( indexField->get ( El.x(), El.y(), El.z() ) != -1 ) {
      cerr << "WARNING: IndexedHeap::append: indexField already set for appended element.\n";
      // cerr << "Element = " << qc::Element << endl;
      //cerr << "heap[ index ] = " << data[ indexField->get( Element.x(), qc::Element.y(), qc::Element.z() ) ] << endl;
    }
    indexField->set ( El.x(), El.y(), El.z(), static_cast<int>(this->data.size( ) - 1 ));
  }

  void push ( const T &El ) {
    append ( El );
    this->upHeap (static_cast<int> ( this->data.size( ) ) - 1 );
  }

  void pop ( T &El ) {
    if ( !this->isEmpty( ) ) {
      El = this->data[ 0 ];
      indexField->set ( El.x(), El.y(), El.z(), -1 );


      this->data[ 0 ] = this->data.back( );
      indexField->set ( this->data[ 0 ].x(), this->data[ 0 ].y(), this->data[ 0 ].z(), 0 );

      this->data.pop_back( );

      this->downHeap ( 0 );
    }
  }

  void swap ( int I1, int I2 ) {
    this->swap_def ( I1, I2 );
    indexField->set ( this->data[ I1 ].x() , this->data[ I1 ].y() , this->data[ I1 ].z() , I1 );
    indexField->set ( this->data[ I2 ].x() , this->data[ I2 ].y() , this->data[ I2 ].z() , I2 );
  }

protected:

  Array<int> *indexField;

};

}

#endif // __QMHEAP_H
