#ifndef __ROWS_H
#define __ROWS_H

#include <vec.h>
#include <gridBase.h>

namespace aol {

/** abstract base class for rows in sparseMatrices
 */
template <typename DataType>
class Row {
protected:
public:

  struct RowEntry {
    RowEntry ( int Col, DataType Value )
        : col ( Col ), value ( Value ) {}
    RowEntry() { }

    int col;
    DataType value;
  };

  // For makeSortedRowEntries
  class RowEntryComp {
  public:
    bool operator () ( const RowEntry& a, const RowEntry& b ) const {
      return a.col < b.col;
    }
  };

  virtual ~Row() {}

  virtual DataType get ( int I, int J ) const = 0;
  virtual void set ( int I, int J, DataType Value ) = 0;
  virtual void add ( int I, int J, DataType Value ) = 0;
  //! Scales row by factor, telling row it is the $I$th row.
  virtual void scale ( int I, DataType factor ) = 0;
  virtual DataType mult ( const Vector<DataType> &V, const int Row ) const = 0;

  typedef DataType ( Row<DataType>::* MultMaskedFctPtrType ) ( const Vector<DataType> &, int Row, const BitVector & );
  virtual DataType multMaskedFunctorTrue     ( const Vector<DataType> &Src, int Row, const BitVector & Mask ) = 0;
  virtual DataType multMaskedFunctorFalse    ( const Vector<DataType> &Src, int Row, const BitVector & Mask ) = 0;
  virtual DataType multMaskedFunctorIdentity ( const Vector<DataType> &Src, int Row, const BitVector & Mask ) = 0;
  virtual DataType multMaskedFunctorNegate   ( const Vector<DataType> &Src, int Row, const BitVector & Mask ) = 0;

  virtual int numNonZeroes() const = 0;

  virtual int numStoredEntries() const {
    throw aol::UnimplementedCodeException ( "aol::Row<DataType>::numStoredEntries not implemented, implement on derived subclasses!", __FILE__, __LINE__ );
    return ( -1 );
  }

  virtual void setZero() = 0;

  //! Scales row by factor. Overload this function if scaling cannot be done without knowing the index of the row.
  void scale ( DataType factor ) {
    scale ( -1, factor );
  }

  //! Add multiple of other row to this row. This is a generic, slow implementation.
  void addMultiple ( const int I, const Row &AddedRow, const DataType Factor = aol::NumberTrait<DataType>::one ) {
    vector<RowEntry> addedRowEntries;
    AddedRow.makeRowEntries ( addedRowEntries, I );

    for ( typename vector<RowEntry>::iterator it = addedRowEntries.begin(); it != addedRowEntries.end(); ++it ) {
      add ( I, it->col, Factor * it->value );
    }
  }

  //! compares in l^inf norm, prints first entry where difference is larger than epsilon
  bool isApproxEqual ( const int I, const Row &ComparedRow, const DataType epsilon ) const {
    vector<RowEntry> entriesA, entriesB;
    this->makeSortedRowEntries ( entriesA, I );
    ComparedRow.makeSortedRowEntries ( entriesB, I );

    unsigned int ia, ib;
    for ( ia = 0, ib = 0; ia < entriesA.size () && ib < entriesB.size(); ) {
      if ( entriesA[ia].col < entriesB[ib].col ) {
        if ( fabs( entriesA[ia].value ) > epsilon ) {
          #ifdef VERBOSE
          std::cerr << "Data differs by " << fabs ( entriesA[ia].value ) << " at pos (" << I << ", " << entriesA[ia].col << "): "
                    << get ( I, entriesA[ia].col ) << " != " << ComparedRow.get ( I, entriesA[ia].col ) << std::endl;
          #endif
          return false;
        }
        ++ia;
      }
      if ( entriesA[ia].col > entriesB[ib].col ) {
        if ( fabs( entriesB[ib].value ) > epsilon ) {
          #ifdef VERBOSE
          std::cerr << "Data differs by " << fabs ( entriesB[ib].value ) << " at pos (" << I << ", " << entriesB[ib].col << "): "
                    << get ( I, entriesB[ib].col ) << " != " << ComparedRow.get ( I, entriesB[ib].col ) << std::endl;
          #endif
          return false;
        }
        ++ib;
      }
      if ( entriesA[ia].col == entriesB[ib].col ) {
        if ( fabs( entriesA[ia].value - entriesB[ib].value ) > epsilon ) {
          #ifdef VERBOSE
          std::cerr << "Data differs by " << fabs ( entriesA[ia].value - entriesB[ib].value ) << " at pos (" << I << ", " << entriesA[ia].col << "): "
                    << get ( I, entriesA[ia].col ) << " != " << ComparedRow.get ( I, entriesA[ia].col ) << std::endl;
          #endif
          return false;
        }
        ++ia; ++ib;
      }
    }

    for ( ; ia < entriesA.size (); ++ia ) {
      if ( fabs( entriesA[ia].value ) > epsilon ) {
        #ifdef VERBOSE
        std::cerr << "Data differs by " << fabs ( entriesA[ia].value ) << " at pos (" << I << ", " << entriesA[ia].col << "): "
                  << get ( I, entriesA[ia].col ) << " != " << ComparedRow.get ( I, entriesA[ia].col ) << std::endl;
        #endif
        return false;
      }
    }

    for ( ; ib < entriesB.size (); ++ib ) {
      if ( fabs( entriesB[ib].value ) > epsilon ) {
        #ifdef VERBOSE
        std::cerr << "Data differs by " << fabs ( entriesB[ib].value ) << " at pos (" << I << ", " << entriesB[ib].col << "): "
                  << get ( I, entriesB[ib].col ) << " != " << ComparedRow.get ( I, entriesB[ib].col ) << std::endl;
        #endif
        return false;
      }
    }

    return true;
  }

  virtual DataType sum ( int I ) = 0;

  //! This is the standard interface to the row's entries, which has to be implemented in each subclass
  virtual void makeRowEntries ( vector<RowEntry > &vec, const int RowNum ) const = 0;

  //! This is the standard implementation using makeRowEntries and std::sort
  //! this function could be overloaded in derived classes.
  virtual void makeSortedRowEntries ( vector<RowEntry > &vec, const int RowNum ) const {
    makeRowEntries ( vec, RowNum );
    sort ( vec.begin (), vec.end (), RowEntryComp () );
    cerr << "Row::makeSortedRowEntries is slow, please overload!" << endl;
  }

  virtual bool checkForNANsAndINFs() const = 0;

  void print ( ostream &/*out*/ ) {}

  Row<DataType>& operator= ( const Row<DataType>& ) {
    throw aol::Exception ( " aol::Row::operator= should not be called.", __FILE__, __LINE__ );
    return *this;
  }
};

template <class T> class MatrixEntry;

template < typename DataType > class SparseMatrixRowIterator;
template < typename DataType > class SparseMatrix;

/** class for general sparse rows that do not necessarily have a special structure
 */
template <typename DataType>
class SparseRow : public Row<DataType> {
  friend class SparseMatrixRowIterator < DataType >;
  friend class SparseMatrix< DataType >;

private:
  //typedef MatrixEntry<DataType> qcCurMatrixEntry;
  typedef typename Row<DataType>::RowEntry qcCurMatrixEntry;

public:
  class ConstIterator : public std::vector <qcCurMatrixEntry>::const_iterator {
    typedef typename std::vector <qcCurMatrixEntry>::const_iterator ParentClass;
    std::vector <qcCurMatrixEntry> &_row;

  public:
    ConstIterator ( aol::SparseRow < DataType > &SparseRow )
      : _row ( SparseRow.row ) {
      ParentClass::operator= ( _row.begin () );
    }

    bool atEnd () const {
      return *this == _row.end ();
    }

    bool notAtEnd () const {
      return *this != _row.end ();
    }
  };

  SparseRow() {}

  virtual ~SparseRow() {}

  DataType get ( int , int J ) const {
    return ( get ( J ) );
  }

  DataType get ( int J ) const {
    typename vector<qcCurMatrixEntry>::const_iterator it;
    for ( it = row.begin(); it != row.end() && it->col <= J; ++it ) {
      if ( ( *it ).col == J ) {
        return ( *it ).value;
      }
    }
    return 0;
  }

  void set ( int , int J, DataType Value ) {
    set ( J, Value );
  }

  void set ( int J, DataType Value ) {
    typename vector<qcCurMatrixEntry>::iterator it;
    if ( Value != 0.0 ) {
      for ( it = row.begin(); it != row.end() && it->col <= J; ++it ) {
        if ( ( *it ).col == J ) {
          ( *it ).value = Value;
          return;
        }
      }

      // Not found, insert before it
      row.insert ( it, qcCurMatrixEntry ( J, Value ) );
    } else {
      for ( it = row.begin(); it != row.end() && it->col <= J; ++it ) {
        if ( ( *it ).col == J ) {
          row.erase ( it );
          break;
        }
      }
    }
  }

  //! return number of nonzero entries in the row (which is NOT the number of stored entries)
  int numNonZeroes ( ) const {
    typename vector<qcCurMatrixEntry>::const_iterator it;
    int nNonZeroes = 0;
    for ( it = row.begin(); it != row.end(); ++it ) {
      if ( it->value != aol::NumberTrait<DataType>::zero )
        ++nNonZeroes;
    }
    return nNonZeroes;
  }

  int numStoredEntries() const {
    return static_cast<int> ( row.size() );
  }

  /** row-vector scalar multiplication
   */
  DataType mult ( const Vector<DataType> &Src, const int /* Row */ ) const {
    typename vector<qcCurMatrixEntry>::const_iterator it;

    DataType dst = ZOTrait<DataType>::zero;
    for ( it = row.begin(); it != row.end(); ++it )
      dst += ( *it ).value * Src.get ( ( *it ).col );
    return ( dst );
  }

  DataType mult ( const Vector<DataType> &src ) const {
    return ( mult ( src, 0 ) ); // or any other index instead of 0
  }

  virtual DataType multMaskedFunctorTrue ( const Vector<DataType> &Src, int Row, const BitVector & Mask ) {
    return multMasked<BitMaskFunctorTrue> ( Src, Row, Mask );
  }
  virtual DataType multMaskedFunctorFalse ( const Vector<DataType> &Src, int Row, const BitVector & Mask ) {
    return multMasked<BitMaskFunctorFalse> ( Src, Row, Mask );
  }
  virtual DataType multMaskedFunctorIdentity ( const Vector<DataType> &Src, int Row, const BitVector & Mask ) {
    return multMasked<BitMaskFunctorIdentity> ( Src, Row, Mask );
  }
  virtual DataType multMaskedFunctorNegate ( const Vector<DataType> &Src, int Row, const BitVector & Mask ) {
    return multMasked<BitMaskFunctorNegate> ( Src, Row, Mask );
  }

  /** masked row-vector scalar mulitplication (if desired only include nodes that are (not) masked
    */
  template <typename BitMaskFunctorType>
  DataType multMasked ( const Vector<DataType> &Src, int,     // 2nd parameter is ignored
                        const BitVector & Mask ) {

    BitMaskFunctorType maskFunctor;
    DataType dst = ZOTrait<DataType>::zero;
    typename vector<qcCurMatrixEntry>::const_iterator it;

    for ( it = row.begin(); it != row.end(); ++it )
      if ( maskFunctor ( Mask[ ( *it ).col] ) )     // i.e. this node is to be included
        dst += ( *it ).value * Src.get ( ( *it ).col );
    return ( dst );
  }


  /** row is mapped to scalar * row
   */
  void scale ( DataType factor ) {
    typename vector<qcCurMatrixEntry>::iterator it;
    for ( it = row.begin(); it != row.end(); ++it ) {
      ( *it ).value *= factor;
    }
  }

  /** computes row sum
   */
  DataType sum ( int ) {
    // I is ignored
    typename vector<qcCurMatrixEntry>::iterator it;
    DataType result = 0;
    for ( it = row.begin(); it != row.end(); ++it ) {
      result += ( *it ).value;
    }
    return ( result );
  }

  void scale ( int , DataType factor ) {
    scale ( factor );
  }

  void add ( int , int J, DataType Value ) {
    add ( J, Value );
  }

  void add ( int J, DataType Value ) {
    typename vector<qcCurMatrixEntry>::iterator it;
    if ( Value != 0.0 ) {
      for ( it = row.begin(); it != row.end() && it->col <= J; ++it ) {
        if ( ( *it ).col == J ) {
          ( *it ).value += Value;
          return;
        }
      }
      // Not found, insert before it
      row.insert ( it, qcCurMatrixEntry ( J, Value ) );
    }
  }

  void setZero() {
    row.erase ( row.begin(), row.end() );
  }

  virtual void makeRowEntries ( vector<typename Row<DataType>::RowEntry > &vec, const int /*RowNum*/ ) const {
    vec = row;
  }

  //! SparseMatrix has already sorted rows
  virtual void makeSortedRowEntries ( vector<typename Row<DataType>::RowEntry > &vec, const int /*RowNum*/ ) const {
    vec = row;
  }

  //! ATTENTION: the name is misleading!
  bool checkForNANsAndINFs() const {
    typename vector<qcCurMatrixEntry>::const_iterator it;
    for ( it = row.begin(); it != row.end(); ++it ) {
      if ( !aol::isFinite ( it->value ) ) {
        return true;
      }
    }
    return false;
  }


  //using Row<DataType>::operator=; // this was protected
  SparseRow<DataType>& operator= ( const SparseRow<DataType> &from ) {
    this->row = from.row;
    return ( *this );
  }


  void eraseZeroEntries ( ) {
    for ( typename vector<qcCurMatrixEntry>::iterator it = row.begin(); it != row.end(); ) {
      if ( it->value == aol::NumberTrait<DataType>::zero ) {
        it = row.erase ( it );
      } else {
        ++it;
      }
    }
  }

protected:

  vector<qcCurMatrixEntry> row;
};



/** A row that only has the diagonal entry.
 *  This is useful if a GenSparseMatrix has many identity rows because they can all be represented by the same instance of DiagonalRow (because the DiagonalRow does not know which row it is).
 *  \author Schwen
 */

template <typename DataType>
class DiagonalRow : public aol::Row< DataType > {

protected:
  DataType _entry;

  mutable bool _entryCurrentlyReturned;

public:

  DiagonalRow ( ) : _entry ( 0 ), _entryCurrentlyReturned ( false ) {}

  explicit DiagonalRow ( const DataType &value ) : _entry ( value ), _entryCurrentlyReturned ( false ) {}

  DiagonalRow<DataType>& operator= ( const DiagonalRow<DataType> &other ) {
    _entry = other._entry;
    _entryCurrentlyReturned = other._entryCurrentlyReturned;
  }

  ~DiagonalRow() {}

  DataType get ( int I, int J ) const {
      return ( I == J ? _entry : aol::NumberTrait<DataType>::zero );
    }

  void set ( int I, int J, DataType Value ) {
    if ( I != J && Value != aol::NumberTrait<DataType>::zero )
      throw aol::Exception ( "aol::DiagonalRow cannot set off-diagonal entry", __FILE__, __LINE__ );

    if ( I == J )
      _entry = Value;
  }

  /** row-vector scalar mulitplication
   */
  DataType mult ( const aol::Vector<DataType> &src, const int rowIndex ) const {
    return ( _entry * src[rowIndex] );
  }

  virtual DataType multMaskedFunctorTrue ( const Vector<DataType> &Src, int Row, const BitVector & Mask ) {
    return multMasked<BitMaskFunctorTrue> ( Src, Row, Mask );
  }
  virtual DataType multMaskedFunctorFalse ( const Vector<DataType> &Src, int Row, const BitVector & Mask ) {
    return multMasked<BitMaskFunctorFalse> ( Src, Row, Mask );
  }
  virtual DataType multMaskedFunctorIdentity ( const Vector<DataType> &Src, int Row, const BitVector & Mask ) {
    return multMasked<BitMaskFunctorIdentity> ( Src, Row, Mask );
  }
  virtual DataType multMaskedFunctorNegate ( const Vector<DataType> &Src, int Row, const BitVector & Mask ) {
    return multMasked<BitMaskFunctorNegate> ( Src, Row, Mask );
  }

  /** masked row-multiplication for the apply(Add)Masked-method
   */
  template <typename BitMaskFunctorType>
  DataType multMasked ( const Vector<DataType> & /*Src*/, int /*Row*/,
                        const BitVector & /*Mask*/ ) {
    throw UnimplementedCodeException ( "DiagonalRow: multMasked is not implemented yet.", __FILE__, __LINE__ );
  }

  /** row is mapped to scalar * row
   */
  void scale ( DataType factor ) {
    _entry *= factor;
  }

  /** compute row sum
   */
  DataType sum ( int /*I*/ ) {
    return ( _entry );
  }

  void scale ( int /*I*/, DataType factor ) {  // 1st parameter is ignored
    scale ( factor );
  }

  void add ( int I, int J, DataType Value ) {
    if ( I != J && Value != aol::NumberTrait<DataType>::zero )
      throw aol::Exception ( "aol::DiagonalRow cannot set off-diagonal entry", __FILE__, __LINE__ );
    _entry += Value;
  }


  void addMultiple ( const int /*rowNum*/, const aol::DiagonalRow<DataType>& other, const DataType factor ) {
    _entry += factor * other._entry;
  }

  using aol::Row<DataType>::addMultiple;

  void setZero() {
    _entry = aol::NumberTrait<DataType>::zero;
  }

  //! ATTENTION: the name is misleading!
  bool checkForNANsAndINFs() const {
    return ( !isFinite ( _entry ) );
  }

  int numNonZeroes() const {
    throw aol::Exception ( "aol::DiagonalRow<DataType>::numNonZeroes() does not return what the method name suggests", __FILE__, __LINE__ );
    return ( 1 );
  }

protected:
  virtual void makeRowEntries ( vector<typename aol::Row<DataType>::RowEntry > &vec, const int RowNum ) const {
    vec.resize ( 1 );
    vec[0].col   = RowNum;
    vec[0].value = _entry;
  }

  // sorting a single entry is not too difficult ...
  virtual void makeSortedRowEntries ( vector<typename aol::Row<DataType>::RowEntry > &vec, const int RowNum ) const {
    vec.resize ( 1 );
    vec[0].col   = RowNum;
    vec[0].value = _entry;
  }
};

}

#endif

