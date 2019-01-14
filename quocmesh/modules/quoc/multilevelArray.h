#ifndef __MULTILEVELARRAY_H
#define __MULTILEVELARRAY_H

#include <quoc.h>
#include <scalarArray.h>
#include <gridBase.h>
#include <op.h>
#include <prolongation.h>
#include <restriction.h>
#include <vectorExtensions.h>

namespace qc {

/** @ingroup multigrid
 */

template <typename DataType, class _ArrayType = qc::Array<DataType>, class ProlongOpType = qc::ProlongOp<DataType>, class RestrictOpType = qc::RestrictOp<DataType, qc::STD_QUOC_RESTRICT>, typename GridType = qc::GridDefinition>
class MultilevelArray {
public:
  typedef _ArrayType ArrayType;
private:
  void init ( ) {
    for ( int level = 0; level <= _fineDepth; ++level ) {
      GridType temp ( level, _dim );
      arrays.push_back ( new ArrayType ( temp ) );
    }
  }
public:
  MultilevelArray ( int FineDepth, Dimension Dim )
      : _fineDepth ( FineDepth ), _cur_level ( FineDepth ), _dim ( Dim ) {
    init();
  }

  explicit MultilevelArray ( const GridType &Grid )
      : _fineDepth ( Grid.getGridDepth() ), _cur_level ( Grid.getGridDepth() ), _dim ( Grid.getDimOfWorld() ) {
    init();
  }

  virtual ~MultilevelArray() {
    typename vector<ArrayType* >::const_iterator it;
    for ( it = arrays.begin(); it != arrays.end(); ++it ) {
      if ( *it ) delete *it;
    }
  }

  virtual void ascent() {
    if ( _cur_level > 0 ) _cur_level--;
  }
  virtual void descent() {
    if ( _cur_level < _fineDepth ) _cur_level++;
  }

  virtual ArrayType* operator() () {
    return arrays[_fineDepth];
  }

  virtual ArrayType &operator[] ( int Level ) {
    return *arrays[ Level ];
  }

  virtual const ArrayType &operator[] ( int Level ) const {
    return *arrays[ Level ];
  }

  /**
   * Clear all arrays of this hierarchy
   * \author Preusser
   */
  virtual void clearAll() {
    cerr << "Clearing ml array" << endl;
    for ( int i = 0; i < _fineDepth; i++ ) {
      arrays[i]->setZero();
    }
  }

  virtual void levProlongate() {
    if ( _cur_level < _fineDepth ) levProlongate ( _cur_level, _cur_level + 1 );
    descent();
  }

  virtual void levProlongate ( int CoarseLevel, int FineLevel ) {
    if ( CoarseLevel < 0 || FineLevel > _fineDepth ) {
      throw aol::Exception ( "ERROR in prolongation!\n" );
    }
    for ( int l = CoarseLevel; l < FineLevel; l++ ) {
      GridType coarse ( l, _dim ), fine ( l + 1, _dim );
      /*       qc::ProlongOp<aol::Vector<DataType> > p( coarse, fine ); */
      ProlongOpType p ( coarse, fine );
      arrays[l+1]->setZero();
      p.apply ( *arrays[ l ], *arrays[ l + 1 ] );
    }
  }

  void levRestrict() {
    levRestrict ( _cur_level - 1, _cur_level );
    ascent();
  }

  virtual void levRestrict ( int CoarseLevel, int FineLevel ) {
    if ( CoarseLevel < 0 || FineLevel > _fineDepth || CoarseLevel >= FineLevel ) {
      throw aol::Exception ( "ERROR in restriction!", __FILE__, __LINE__ );
    }
    for ( int l = FineLevel; l > CoarseLevel; l-- ) {
      GridType coarse ( l - 1, _dim ), fine ( l, _dim );
      /*       qc::RestrictOp<aol::Vector<DataType> > r( coarse, fine ); */
      RestrictOpType r ( coarse, fine );
      arrays[l-1]->setZero();
      r.apply ( *arrays[ l ], *arrays[ l - 1 ] );
    }
  }

  virtual ArrayType& current() {
    return *arrays[ _cur_level ];
  }

  virtual const ArrayType& current() const {
    return *arrays[ _cur_level ];
  }

  virtual ArrayType* operator->() {
    return arrays[ _cur_level ];
  }

  int getCurLevel() {
    return _cur_level;
  }

  void setCurLevel ( int Level ) {
    _cur_level = Level;
  }

  int getDepth() {
    return _fineDepth;
  }

  Dimension getDimension() {
    return _dim;
  }

protected:
  vector<ArrayType* > arrays;
  int _fineDepth, _cur_level;
  Dimension _dim;
};


template<typename DataType, class ArrayType = qc::Array<DataType>, class ProlongOpType = qc::ProlongOp<DataType>, class RestrictOpType = qc::RestrictOp<DataType, qc::STD_QUOC_RESTRICT>, typename GridType = qc::GridDefinition>
class MultiDimMultilevelArray {
public:
  MultiDimMultilevelArray ( int FineDepth, Dimension Dim, int Components = 1 )
      : _cur_level ( FineDepth ), _num_comp ( Components ), _fineDepth ( FineDepth ) {
    for ( int i = 0; i < _num_comp; i++ ) {
      _ml_arrays.push_back ( new MultilevelArray<DataType, ArrayType, ProlongOpType, RestrictOpType, GridType> ( FineDepth, Dim ) );
    }
  }

  explicit MultiDimMultilevelArray ( const GridType &Grid, int Components = 1 )
      : _cur_level ( Grid.getGridDepth() ), _num_comp ( Components ), _fineDepth ( Grid.getGridDepth() ) {
    for ( int i = 0; i < _num_comp; i++ ) {
      _ml_arrays.push_back ( new MultilevelArray<DataType, ArrayType, ProlongOpType, RestrictOpType, GridType> ( Grid ) );
    }
  }

  virtual ~MultiDimMultilevelArray() {
    for ( int i = 0; i < _num_comp; i++ ) {
      delete _ml_arrays[ i ];
    }
  }

  void clear() {
    for ( int i = 0; i < _num_comp; i++ ) {
      _ml_arrays[ i ]->current().setZero();
    }
  }

  ArrayType &getArray ( int Comp ) {
    return ( *_ml_arrays[ Comp ] ) [ _cur_level ];
  }

  const ArrayType &getArray ( int Comp ) const {
    return ( *_ml_arrays[ Comp ] ) [ _cur_level ];
  }

  ArrayType &getArray ( int Comp, int Level ) {
    return ( *_ml_arrays[ Comp ] ) [ Level ];
  }

  const ArrayType &getArray ( int Comp, int Level ) const {
    return ( *_ml_arrays[ Comp ] ) [ Level ];
  }

  MultilevelArray<DataType, ArrayType, ProlongOpType, RestrictOpType, GridType> &getMultilevelArray ( int Comp ) {
    return *_ml_arrays[ Comp ];
  }

  void getArrayReferences ( aol::RandomAccessContainer<ArrayType> &Arrays ) const {
    Arrays.reallocate ( 0 );
    for ( int i = 0; i < _num_comp; i++ )
      Arrays.constructDatumAndPushBack( getArray(i), aol::FLAT_COPY );
  }

  void getArrays ( int Level, vector<ArrayType* > &Arrays ) {
    if ( Arrays.size() == _num_comp ) {
      for ( int i = 0; i < _num_comp; i++ ) {
        Arrays[ i ] = & ( ( *_ml_arrays[ i ] ) [ Level ] );
      }
    } else {
      Arrays.erase ( Arrays.begin(), Arrays.end() );
      for ( int i = 0; i < _num_comp; i++ ) {
        Arrays.push_back ( & ( ( *_ml_arrays[ i ] ) [ Level ] ) );
      }
    }
  }

  void appendReferencesTo ( aol::MultiVector<DataType> &MDest ) {
    for ( int i = 0; i < _num_comp; i++ ){
      MDest.appendReference( getArray( i ) );
    }
  }

  int getCurLevel() const {
    return _cur_level;
  }

  void setCurLevel ( int Level ) {
    _cur_level = Level;
    for ( int i = 0; i < _num_comp; i++ ) {
      _ml_arrays[ i ]->setCurLevel ( Level );
    }
  }

  void levRestrict() {
    if ( _cur_level > 0 ) {
      _cur_level--;
      for ( int i = 0; i < _num_comp; i++ ) {
        _ml_arrays[ i ]->levRestrict();
      }
    } else {
      // error
    }
  }



  virtual void levRestrict ( int CoarseLevel, int FineLevel ) {
    for ( int i = 0; i < _num_comp; i++ ) {
      _ml_arrays[ i ]->levRestrict ( CoarseLevel, FineLevel );
    }
    _cur_level = CoarseLevel;
  }

  void levProlongate() {
    if ( _cur_level < _fineDepth ) {
      _cur_level++;
      for ( int i = 0; i < _num_comp; i++ ) {
        _ml_arrays[ i ]->levProlongate();
      }
    } else {
      // error
    }
  }

  virtual void levProlongate ( int CoarseLevel, int FineLevel ) {
    for ( int i = 0; i < _num_comp; i++ )
      _ml_arrays[ i ]->levProlongate( CoarseLevel, FineLevel );
    _cur_level = FineLevel;
  }

  int getDepth() const {
    return _fineDepth;
  }

  int numComponents() const {
    return _num_comp;
  }

protected:
  vector<MultilevelArray<DataType, ArrayType, ProlongOpType, RestrictOpType, GridType>* > _ml_arrays;
  int _cur_level;
  int _num_comp;
  int _fineDepth;
};

}

#endif
