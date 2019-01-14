#ifndef __MASKEDVECTOR_H
#define __MASKEDVECTOR_H

#include <vec.h>
#include <op.h>
#include <bitVector.h>

namespace aol {

/**
 * class for managing a subset of indices for masked vectors
 * using run-length-encoding. mainly provides an iterator to all used dofindices.
 * \author Droske
 * \see MaskedVector
 */
class DofMask {
protected:
  // mask stores alternatingly the number of used/unused dofs
  // the first entry is number of used dofs
  vector<int> _mask;

public:
  DofMask ( int size ) : end_it ( *this ) {
    _mask.push_back ( size );
    end_it._maskIt = _mask.end();
    end_it._dofIndex = -1;
  }

public:
  DofMask() : end_it ( *this ) {
    end_it._maskIt = _mask.end();
    end_it._dofIndex = -1;
  }

public:
  DofMask ( const aol::BitVector &bf ) : end_it ( *this ) {
    end_it._maskIt = _mask.end();
    end_it._dofIndex = -1;
    createDofMaskFromBitVector ( bf );
  }

public:
  // copy-constructor
  DofMask ( const DofMask &dm )
      : _mask ( dm._mask ), end_it ( ( *this ) ) {}

public:
  /**
   * adds a pair of used and unused dof numbers. call this function to successively to initialize
   * the dofmask.
   * \param nonzeros number of used dofs
   * \param zeros number of unused dofs
   */
  void add ( int nonzeros, int zeros ) {
    _mask.push_back ( nonzeros );
    _mask.push_back ( zeros );
    end_it._maskIt = _mask.end();
  }

  class iterator {
    friend class DofMask;

    int _dofIndex;
    int _lastIndexInSequence;
    const DofMask &_dofMask;
    vector<int>::const_iterator _maskIt;

  public:
    iterator ( const DofMask &dofMask )
        : _dofMask ( dofMask ) {}

  public:
    iterator& operator++ ( ) {
      if ( ++_dofIndex >= _lastIndexInSequence ) {
        _lastIndexInSequence += * ( ++_maskIt );
        _dofIndex += *_maskIt;
        ++_maskIt;
        if ( _maskIt != _dofMask._mask.end() )
          _lastIndexInSequence += * ( _maskIt ); // else set to undefined value ...
      }
      return *this;
    }

    iterator& incrementUntilIndex ( int i ) {
      if ( i <= _dofIndex ) return *this;
      _dofIndex = i;
      while ( _dofIndex >= _lastIndexInSequence ) {
        _lastIndexInSequence += * ( ++_maskIt );   // maskIt points to non-zero count, so here we add the count for the zeros
        _dofIndex = max ( _dofIndex, _lastIndexInSequence );  // skip to next non-zero
        ++_maskIt;                                  // point to non-zero
        if ( _maskIt != _dofMask._mask.end() )
          _lastIndexInSequence += * ( _maskIt );
        else
          _lastIndexInSequence = _dofIndex + 1;
      }
      return *this;
    }

    bool operator== ( iterator it ) const {
      return ( ( ( _maskIt == _dofMask._mask.end() ) && ( _maskIt == it._maskIt ) ) ||
               ( _dofIndex == it._dofIndex ) );
    }

    bool operator!= ( iterator it ) const {
      return !operator== ( it );
    }

    int operator*() const {
      return _dofIndex;
    }

    const int *operator->() const {
      return &_dofIndex;
    }

  };

  iterator begin() const {
    iterator it ( *this );

    it._maskIt = _mask.begin();

    if ( * ( it._maskIt ) == 0 ) {
      it._maskIt++;
      it._dofIndex = * ( it._maskIt );
      it._maskIt++;
    } else {
      it._dofIndex = 0;
    }
    it._lastIndexInSequence = it._dofIndex + * ( it._maskIt );
    return it;
  }

  const iterator& end() const {
    return end_it;
  }

private:
  iterator end_it;

public:

  template <typename DataType>
  void createDofMaskFromVector ( const aol::Vector<DataType> &vec, bool verbose = true ) {
    _mask.clear();
    int numDofs = 0;
    int nonzeros = 0, zeros = 0;
    bool count_nonzeros = true;
    for ( int i = 0; i < vec.size(); i++ ) {
      if ( count_nonzeros ) {
        if ( vec[i] != 0 ) {
          nonzeros++;
          numDofs++;
        } else {
          count_nonzeros = false;
          zeros++;
        }
      } else { // count zeros
        if ( vec[i] != 0 ) {
          numDofs++;
          add ( nonzeros, zeros );
          count_nonzeros = true;
          nonzeros = 1;
          zeros = 0;
        } else {
          zeros++;
        }
      }
    }
    add ( nonzeros, zeros );
    if ( verbose ) cerr << "number of active dofs = " << numDofs << endl;
  }
  void createDofMaskFromBitVector ( const aol::BitVector &bf, bool verbose = true ) {
    _mask.clear();
    int numDofs = 0;
    int nonzeros = 0, zeros = 0;
    bool count_nonzeros = true;
    for ( int i = 0; i < bf.size(); i++ ) {
      if ( count_nonzeros ) {
        if ( bf[i] != 0 ) {
          nonzeros++;
          numDofs++;
        } else {
          count_nonzeros = false;
          zeros++;
        }
      } else { // count zeros
        if ( bf[i] != 0 ) {
          numDofs++;
          add ( nonzeros, zeros );
          count_nonzeros = true;
          nonzeros = 1;
          zeros = 0;
        } else {
          zeros++;
        }
      }
    }
    add ( nonzeros, zeros );
    if ( verbose ) cerr << "number of active dofs = " << numDofs << endl;

  }

};

/**
 * vector class which operates only on a subset of all global indices.
 * \author Droske
 */
template <typename DataType>
class MaskedVector : public Vector<DataType> {

protected:
  const DofMask &_dofMask;

public:
  MaskedVector ( int size, const DofMask &dofMask )
      : Vector<DataType> ( size ), _dofMask ( dofMask ) {}

  MaskedVector ( DataType *values, int size, const DofMask &dofMask )
      : Vector<DataType> ( values, size ), _dofMask ( dofMask ) {}

public:
  // copy constructor
  MaskedVector ( const MaskedVector<DataType> &v, CopyFlag copyFlag = DEEP_COPY )
      : Vector<DataType> ( v, copyFlag ), _dofMask ( v._dofMask ) {}

public:
  const DofMask& getDofMask() const {
    return _dofMask;
  }

  MaskedVector<DataType>& operator-= ( const Vector<DataType> &vec ) {
    if ( vec.size() == this->_size ) {
      for ( DofMask::iterator it = _dofMask.begin(); it != _dofMask.end(); ++it ) {
        ( *this ) [*it] -= vec[*it];
      }
    } else {
      throw aol::Exception ( "aol::MaskedVector<DataType>::operator-= : incompatible vector lengths!", __FILE__, __LINE__ );
    }

    return *this;
  }

  MaskedVector<DataType>& addToAll ( DataType scalar ) {
    for ( DofMask::iterator it = _dofMask.begin(); it != _dofMask.end(); ++it ) {
      ( *this ) [*it] += scalar;
    }
    return *this;
  }


  MaskedVector<DataType>& operator+= ( const Vector<DataType> &vec ) {
    if ( vec.size() == this->_size ) {
      for ( DofMask::iterator it = _dofMask.begin(); it != _dofMask.end(); ++it ) {
        ( *this ) [*it] += vec[*it];
      }
    } else {
      throw aol::Exception ( "aol::MaskedVector<DataType>::operator+= : incompatible vector lengths!", __FILE__, __LINE__ );
    }

    return *this;
  }

  MaskedVector<DataType>& operator*= ( const DataType scalar ) {
    for ( DofMask::iterator it = _dofMask.begin(); it != _dofMask.end(); ++it ) {
      ( *this ) [*it] *= scalar;
    }
    return *this;
  }

  MaskedVector<DataType>& operator/= ( const DataType scalar ) {
    for ( DofMask::iterator it = _dofMask.begin(); it != _dofMask.end(); ++it ) {
      ( *this ) [*it] /= scalar;
    }
    return *this;
  }
  using aol::Vector<DataType>::operator/= ;

  bool operator== ( const Vector<DataType> &vec ) const {
    if ( vec.size() == this->_size ) {
      for ( DofMask::iterator it = _dofMask.begin(); it != _dofMask.end(); ++it ) {
        if ( ( *this ) [*it] != vec[*it] ) return false ;
      }
    } else {
      throw aol::Exception ( "aol::MaskedVector<DataType>::operator== : incompatible vector lengths!", __FILE__, __LINE__ );
    }

    return true;
  }

  bool operator!= ( const Vector<DataType> &vec ) const {
    return ! ( *this == vec );
  }

  MaskedVector<DataType>& operator= ( const Vector<DataType> &vec ) {
    if ( vec.size() == this->_size ) {
      for ( DofMask::iterator it = _dofMask.begin(); it != _dofMask.end(); ++it ) {
        ( *this ) [*it] = vec[*it];
      }
    } else {
      throw aol::Exception ( "aol::MaskedVector<DataType>::operator= : incompatible vector lengths!", __FILE__, __LINE__ );
    }

    return *this;
  }

  MaskedVector<DataType>& operator= ( const MaskedVector<DataType> &vec ) {
    if ( vec.size() == this->_size ) {
      for ( DofMask::iterator it = _dofMask.begin(); it != _dofMask.end(); ++it ) {
        ( *this ) [*it] = vec[*it];
      }
    } else {
      throw aol::Exception ( "aol::MaskedVector<DataType>::operator= : incompatible vector lengths!", __FILE__, __LINE__ );
    }

    return *this;
  }


  MaskedVector<DataType>& addMultiple ( const Vector<DataType> &vec, DataType scalar ) {
    if ( vec.size() == this->_size ) {
      for ( DofMask::iterator it = _dofMask.begin(); it != _dofMask.end(); ++it ) {
        ( *this ) [*it] += scalar * vec[*it];
      }
    } else {
      throw aol::Exception ( "aol::MaskedVector<DataType>::addMultiple : incompatible vector lengths!", __FILE__, __LINE__ );
    }

    return *this;
  }

  DataType operator* ( const Vector<DataType> &vec ) const {
    DataType dot = 0.;
    if ( vec.size() == this->_size ) {
      for ( DofMask::iterator it = _dofMask.begin(); it != _dofMask.end(); ++it ) {
        dot += ( *this ) [*it] * vec[*it];
      }
    } else {
      throw aol::Exception ( "aol::MaskedVector<DataType>::operator* : incompatible vector lengths!", __FILE__, __LINE__ );
    }

    return dot;
  }
};

/**
 * allows to convert a Op<VecType> to a Op<DerivedVecType>
 * which calls the apply of the Op<VecType>-operator.
 * used for technical reasons.
 * \author Droske
 */
template <typename VecType, typename DerivedVecType>
class CastOp : public Op<DerivedVecType> {
protected:
  const Op<VecType> &_op;

public:
  CastOp ( const Op<VecType> &op ) : _op ( op ) {}

public:
  virtual void apply ( const DerivedVecType &Arg, DerivedVecType &Dest ) const {
    _op.apply ( Arg, Dest );
  }

  virtual void applyAdd ( const DerivedVecType &Arg, DerivedVecType &Dest ) const {
    _op.applyAdd ( Arg, Dest );
  }
};

} // end namespace aol

#endif
