#ifndef __BANDMATRIX_H
#define __BANDMATRIX_H

#include <op.h>
#include <matrix.h>
#include <gridBase.h>

namespace aol {
/** \brief A class for general band matrices.
 *
 *  Needs an aol::Vector<int> of the offsets of the nonzero diagonals from the main diagonal,
 *  for example [-1, 0, 1] is a tridiagonal matrix.
 *  These diagonals need not be right next to the main diagonal, e. g. [-10, 0, 10].
 *  Entries will be stored diagonal-wise.
 *  Focus is on speed, not minimal memory usage.
 *  Write/Use subclasses if you know the structure of the matrix beforehand!
 *  Use of virtual may not be optimal yet. TODO: think about this.
 *  \author Schwen
 *  \ingroup Matrix
 */
template < typename _DataType >
class GenBandMatrix : public Matrix<_DataType> {
public:
  typedef _DataType DataType;

protected:
  aol::Vector<int> _startApply, _endApply, _globalToLocal, _localToGlobal;
  DataType*  _pData;
  int _nDiags, _sizeReserved;
  static const int UNDEFINED_ENTRY;

  //! internal function to compute index in data vector
  inline int map_index ( const int diag_no, const int row ) const {
    return ( row * _nDiags + diag_no );    // row-wise: for apply, this should be better for cache usage than diagonal-wise storage
  }

public:
  //! constructor
  GenBandMatrix ( const int Rows, const int Cols, const aol::Vector<int> &Offsets );

  //! constructor with no structure information given: assumes zero matrix of given size, no bands. Can define band structure later on via reallocate.
  GenBandMatrix ( const int Rows, const int Cols );

  //! constructor from GridDefinition: only uses number of nodes and boundary mask
  explicit GenBandMatrix ( const qc::GridStructure & Grid );

  //! default constructor
  GenBandMatrix ( ) ;

  //! copy constructor
  explicit GenBandMatrix ( const GenBandMatrix<DataType> &Other );

  //! destructor
  virtual ~GenBandMatrix();

  //! reallocate Matrix
  virtual void reallocate ( const int Rows, const int Cols, const aol::Vector<int>& Offsets );

  using aol::Matrix<DataType>::reallocate;

  //! Returns a matrix entry.
  virtual DataType get ( int I, int J ) const {

#ifdef BOUNDS_CHECK
      this->boundsCheck ( I, J, "aol::GenBandMatrix::get: Index out of bounds", __FILE__, __LINE__ );
#endif

      if ( _globalToLocal[J - I + this->getNumCols() - 1] == UNDEFINED_ENTRY ) {
        return ( aol::NumberTrait<DataType>::zero );
      } else {
        return ( _pData[ map_index ( _globalToLocal[J - I + this->getNumCols() - 1] , I ) ] );
      }
    }

  //! Sets a band matrix entry. Non-existent entries are fixed to 0.
  virtual void set ( int I, int J, DataType Value ) {

#ifdef BOUNDS_CHECK
    this->boundsCheck ( I, J, "aol::GenBandMatrix::set: Index out of bounds", __FILE__, __LINE__ );

    if ( ( _globalToLocal[J - I + this->getNumCols() - 1] == UNDEFINED_ENTRY ) && ( Value != aol::NumberTrait<DataType>::zero ) ) {
      cerr << "trying to set " << I << ", " << J << ", value " << Value << "; this entry does not exist!" << endl;
      throw aol::Exception ( "GenBandMatrix::set: trying to set illegal position to value != 0", __FILE__, __LINE__ );
    }
#endif

    _pData[ map_index ( _globalToLocal[J - I + this->getNumCols() - 1] , I ) ] = Value;

  }

  //! Adds to a band matrix entry. Non-existent entries are fixed to 0.
  virtual void add ( int I, int J, DataType Value ) {

#ifdef BOUNDS_CHECK
    this->boundsCheck ( I, J, "aol::GenBandMatrix::add: Index out of bounds", __FILE__, __LINE__ );

    if ( ( _globalToLocal[J - I + this->getNumCols() - 1] == UNDEFINED_ENTRY ) && ( Value != aol::NumberTrait<DataType>::zero ) ) {
      cerr << "trying to add to " << I << ", " << J << ", value " << Value << "; this entry does not exist!" << endl;
      throw aol::Exception ( "GenBandMatrix::add: trying to set illegal position to value != 0", __FILE__, __LINE__ );
    }
#endif

    _pData[ map_index ( _globalToLocal[J - I + this->getNumCols() - 1] , I ) ] += Value;

  }


  //! Multiply matrix by scalar
  virtual GenBandMatrix<DataType>& operator*= ( const DataType Alpha );

  /** Approximate comparison
   */
  bool isApproxEqual ( const GenBandMatrix<DataType> &other, DataType Epsilon );

  GenBandMatrix<DataType>& addMultiple ( const GenBandMatrix<DataType> &Mat, DataType Factor );

  using aol::Matrix<DataType>::operator+=;

  GenBandMatrix<DataType>& operator+= ( const GenBandMatrix<DataType> &Mat ) {
    return addMultiple ( Mat, aol::NumberTrait<DataType>::one );
  }

  GenBandMatrix<DataType>& operator= ( const GenBandMatrix<DataType> &Mat );

  //! Matrix-vector multiplication.
  //! This works rather efficiently for the general case but you might want to write a subclass for specially structured band matrices and implement a special apply method tailored to this structure.
  virtual void apply ( const Vector<DataType> &src, Vector<DataType> &dst ) const;

  //! Matrix-vector multiplication with masking functionality.
  //! Differently from apply, this function is not re-implemented
  //! in subclasses.
  void applyMasked ( const Vector<DataType> &Arg, Vector<DataType> &Dest,
                     const BitVector & Mask, IncludeWriteMode applyMode ) const;

  template <typename PreimBitMaskFunctorType, typename ImageBitMaskFunctorType>
  void applyMasked ( const Vector<DataType> &Arg, Vector<DataType> &Dest,
                     const BitVector & Mask ) const {
    // might want to compare sizes
    PreimBitMaskFunctorType preimMaskFunctor;
    ImageBitMaskFunctorType imageMaskFunctor;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i <  this->getNumRows() ; ++i )
      if ( imageMaskFunctor ( Mask[i] ) ) {
        DataType result = aol::ZOTrait<DataType>::zero;
        for ( int j = _startApply[i]; j < _endApply[i]; ++j )
          if ( preimMaskFunctor ( Mask[i + _localToGlobal[j]] ) )
            // note: first indices have different order!
            result += _pData[map_index ( j, i ) ] * Arg[i + _localToGlobal[j]];
        Dest[i] = result;
      }
  }

  //! This works rather efficiently for the general case but you might want
  //! to write a subclass for specially structured band matrices and implement
  //! a special apply method tailored to this structure.
  virtual void applyAdd ( const Vector<DataType> &src, Vector<DataType> &dst ) const;

  //! Matrix-vector multiplication with masking functionality.
  //! Differently from applyAdd, this function is not re-implemented
  //! in subclasses.
  void applyAddMasked ( const Vector<DataType> &Arg, Vector<DataType> &Dest,
                        const BitVector & Mask, IncludeWriteMode applyMode ) const;

  template <typename PreimBitMaskFunctorType, typename ImageBitMaskFunctorType>
  void applyAddMasked ( const Vector<DataType> &Arg, Vector<DataType> &Dest,
                        const BitVector & Mask ) const {
    // might want to compare sizes
    PreimBitMaskFunctorType preimMaskFunctor;
    ImageBitMaskFunctorType imageMaskFunctor;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i <  this->getNumRows() ; ++i )
      if ( imageMaskFunctor ( Mask[i] ) ) {
        DataType result = aol::ZOTrait<DataType>::zero;
        for ( int j = _startApply[i]; j < _endApply[i]; ++j ) {
          if ( preimMaskFunctor ( Mask[i + _localToGlobal[j]] ) )
            // note: first indices have different order!
            result += _pData[map_index ( j, i ) ] * Arg[i + _localToGlobal[j]];
        }
        Dest[i] += result;
      }
  }

  //! Sets to all entries to zero
  void setZero ( ) {
    memset ( _pData, 0, _nDiags *  this->getNumRows()  * sizeof ( DataType ) );
  }


  void makeSortedRowEntries ( vector<typename Row<DataType>::RowEntry > &vec, const int RowNum ) const {
    makeRowEntries ( vec, RowNum );
  }


  void makeRowEntries ( vector<typename Row<DataType>::RowEntry > &vec, const int RowNum ) const {
    vec.resize ( _endApply[ RowNum ] - _startApply[ RowNum ] ); // and all of these will be set
    for ( int k = 0, j = _startApply[RowNum]; j < _endApply[RowNum]; ++j, ++k ) {
      vec[k].col = RowNum + _localToGlobal[j];
      vec[k].value = _pData[ map_index ( j, RowNum ) ];
    }
  }

  //! Comparison of band structure
  bool hasSameStructureAs ( const aol::GenBandMatrix<DataType> &Mat ) const;

  //! efficient reimplementation of setRowToZero
  void setRowToZero ( const int i );

  //! efficient reimplementation of setColToZero
  void setColToZero ( const int i );

  int get_first_offset() {
    return ( _localToGlobal[0] ) ;
  }

  int get_last_offset() {
    return ( _localToGlobal[_nDiags-1] );
  }

};


/** \brief A class for band matrices that have adjacent bands, e. g. centered at the main diagonal.
 *  \author Schwen
 *  \ingroup Matrix
 */
template<typename DataType>
class BandMatrix : public GenBandMatrix<DataType> {

public:
  //! Constructor adapted to this matrix type
  BandMatrix ( const int M, const int N, const int FirstOffset, const int NBands );

public:
  //! optimized apply for a matrix of this structure
  virtual void apply ( const Vector<DataType> &src, Vector<DataType> &dst ) const;

private:
  //! optimized applyAdd for a matrix of this structure
  virtual void applyAdd ( const Vector<DataType> &/*src*/, Vector<DataType> &/*dst*/ ) const {
    throw aol::Exception ( "aol::BandMatrix::applyAdd not implemented yet", __FILE__, __LINE__ );
    // work this out
    // is this not just dst[i] += result instead of = result??
  }

protected:
  int _firstOffset, _nBands, _firstComplete, _lastComplete;
};



/** \brief A class for tridiagonal matrices.
 *  \author Schwen
 *  \ingroup Matrix
 */
template<typename DataType>
class TriBandMatrix : public GenBandMatrix<DataType> {

public:
  //! Constructor adapted to this matrix type
  explicit TriBandMatrix ( const int N );

  explicit TriBandMatrix ( const qc::GridDefinition & Grid );

  explicit TriBandMatrix ( const qc::GridSize<qc::QC_1D> &Size );

  //! Constructor for compatibility
  TriBandMatrix ( const int Nx, const int Ny );

public:
  //! optimized apply for a matrix of this structure
  virtual void apply ( const Vector<DataType> &Src, Vector<DataType> &Dst ) const;

  //! optimized applyAdd for a matrix of this structure
  virtual void applyAdd ( const Vector<DataType> &Src, Vector<DataType> &Dst ) const;

};


/** \brief A class for band matrices with two lower, main and one upper diagonals.
 *  \author Schwen
 *  \ingroup Matrix
 */
template<typename DataType>
class LQuadBandMatrix : public GenBandMatrix<DataType> {
public:
  //! Constructor adapted to this matrix type
  explicit LQuadBandMatrix ( const int N );

public:
  //! optimized apply for a matrix of this structure
  virtual void apply ( const Vector<DataType> &Src, Vector<DataType> &Dst ) const;

  //! optimized applyAdd for a matrix of this structure
  virtual void applyAdd ( const Vector<DataType> &Src, Vector<DataType> &Dst ) const;

};


/** \brief A class for band matrices with one lower, main and two upper diagonals.
 *  \author Schwen
 *  \ingroup Matrix
 */
template<typename DataType>
class RQuadBandMatrix : public GenBandMatrix<DataType> {
public:
  //! Constructor adapted to this matrix type
  explicit RQuadBandMatrix ( const int N );

public:
  //! optimized apply for a matrix of this structure
  virtual void apply ( const Vector<DataType> &Src, Vector<DataType> &Dst ) const;

  //! optimized applyAdd for a matrix of this structure
  virtual void applyAdd ( const Vector<DataType> &Src, Vector<DataType> &Dst ) const;

};

}

#endif
