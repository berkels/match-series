#ifndef __INDEXTRANSLATIONMATRIX_H
#define __INDEXTRANSLATIONMATRIX_H

#include <matrix.h>

#ifdef USE_CPP11
#define AUTOPTR std::unique_ptr
#else
#define AUTOPTR auto_ptr
#endif

namespace aol {

/*    CONTENTS OF THIS FILE:
 *
 *    class TransparentIndexTranslationMatrix: 1.1 declaration
 *                                             2.1 implementation
 *
 *    class IndexTranslationMatrix: 1.2 declaration
 *                                  2.2 implementation
 */

/****************************************************************************
 *
 *        1.1 CLASS TransparentIndexTranslationMatrix
 */
/**
 *  \brief Wrapper for Matrices that are intended to work only on a
 *         small set of vector components. Small matrix is visible
 *         through the wrapper.
 *
 *  The TransparentIndexTranslationMatrix<MatrixType> is derived from
 *  MatrixType.
 *  This underlying object stores a small matrix and can, due to public
 *  inheritance, be accessed directly via apply() and applyAdd(). This
 *  class only extends functionality by storing a translation
 *  vector (more precisely: one translation for preimage and one for image
 *  domain vectors). It has translation functions translateSmallToLarge()
 *  and translateLargeToSmall() in usual apply() signature.

 *  The procedure
 *  of translation large->small, apply on small vector and back translation
 *  small->large can be done by apply{Small|Large}to{Small|Large}.
 *  It gets arg and dest parameters as usual and performs the obvious
 *  translations. You may choose if for all vector entries that do not
 *  belong to the translated components' set, you want to apply identity
 *  (that means, vector entries from arg to dest are copied) or leave
 *  them unchanged. For the translation procedure, additional
 *  vectors are constructed in the constructor. Therefore, translation
 *  during apply{Small|Large}to{Small|Large} calls does not require
 *  additional memory.
 *
 *  If you want to hide the apply() and applyAdd() from the underlying
 *  matrix, use an IndexTranslationMatrix instead of
 *  TransparentIndexTranslationMatrix.
 *
 *  \author von Deylen
 */

template < typename MatrixType,
           typename LargeDomainType = typename MatrixType::DomainType,
           typename LargeRangeType  = typename MatrixType::RangeType >
class TransparentIndexTranslationMatrix : public MatrixType {
public:
  typedef typename MatrixType::DataType   DataType;
  typedef typename MatrixType::DomainType DomainType;
  typedef typename MatrixType::RangeType  RangeType;
  typedef typename MatrixType::MaskType   MaskType;

  enum ApplyModeOnUnusedNodes {
    LEAVE_UNCHANGED,
    APPLY_IDENTITY
  };

  // *** CONSTRUCTORS ***

  TransparentIndexTranslationMatrix ( int dimSmall, int dimLarge, const Vector<int> & smallToLarge,
                                      ApplyModeOnUnusedNodes mode = APPLY_IDENTITY );

  TransparentIndexTranslationMatrix ( int dimSmallRows, int dimSmallCols,
                                      int dimLargeRows, int dimLargeCols,
                                      const Vector<int> & smallToLargePreimage,
                                      const Vector<int> & smallToLargeImage,
                                      ApplyModeOnUnusedNodes mode = APPLY_IDENTITY );

  template<class LargeMatrixType>
  void FillFromLarge ( const LargeMatrixType & large );

  // *** VECTOR TRANSLATIONS ***

  void translateLargeToSmallPreimage ( const LargeDomainType & large, DomainType & small ) const;
  void translateLargeToSmallImage    ( const LargeRangeType &  large, RangeType &  small ) const;
  void translateSmallToLargePreimage ( const DomainType & small, LargeDomainType & large ) const;
  void translateSmallToLargeImage    ( const RangeType &  small, LargeRangeType &  large ) const;

  void translateLargeToSmall         ( const LargeDomainType & large, DomainType & small ) const;
  void translateSmallToLarge         ( const DomainType & large, LargeDomainType & small ) const;

  // *** APPLY ROUTINES ***

  //! If the ApplyModeOnUnusedNodes is set to APPLY_IDENTITY, the small vector
  //! is virtually enlarged with zeros (ie, dest will be set to zero on
  //! unused nodes).
  void applySmallToLarge ( const DomainType & arg, LargeRangeType dest ) const;
  void applySmallToSmall ( const DomainType & arg, RangeType dest ) const;
  void applyLargeToSmall ( const LargeDomainType & arg, RangeType dest ) const;
  void applyLargeToLarge ( const LargeDomainType & arg, LargeRangeType dest ) const;

  void applyAddSmallToLarge ( const DomainType & arg, LargeRangeType dest ) const;
  void applyAddSmallToSmall ( const DomainType & arg, RangeType dest ) const;
  void applyAddLargeToSmall ( const LargeDomainType & arg, RangeType dest ) const;
  void applyAddLargeToLarge ( const LargeDomainType & arg, LargeRangeType dest ) const;

  // *** GET DIMENSION FUNCTIONS ***

  int getNumColsSmall() const    {  return this->getNumCols();  }
  int getNumRowsSmall()    const    {  return this->getNumRows();     }
  int getNumColsLarge() const    {  return _dimLargeCols;          }
  int getNumRowsLarge()    const    {  return _dimLargeRows;          }

  int getSizeSmall() const;
  int getSizeLarge() const;

  const Vector<int> getSmallToLargePreimage() const;
  const Vector<int> getSmallToLargeImage() const;
  const Vector<int> getSmallToLarge () const;

  ApplyModeOnUnusedNodes getApplyModeOnUnusedNodes() const;
  void setApplyModeOnUnusedNodes ( ApplyModeOnUnusedNodes mode );

  const MaskType & getUsedNodesMaskPreimage() const;
  const MaskType & getUsedNodesMaskImage() const;

  bool isTranslationSymmetric() const;

private:
  Vector<int> _preimSmallToLarge;
  Vector<int> _imSmallToLarge;

  DomainType _preimTempVector;
  RangeType  _imTempVector;

  mutable AUTOPTR<MaskType> _preimUsedNodesMaskPtr;
  mutable AUTOPTR<MaskType> _imUsedNodesMaskPtr;

  int _dimLargeRows;
  int _dimLargeCols;

  ApplyModeOnUnusedNodes _unusedNodesMode;

  bool _translationIsSymmetric;
}; // end of class TransparentIndexTranslationMatrix.

/****************************************************************************
 *
 *        1.2 CLASS IndexTranslationMatrix
 */
/**
 *  \brief Wrapper for Matrices that are intended to work only on a
 *         small set of vector components. In contrast to the behaviour of
 *         TransparentIndexTranslationMatrix, the underlying small
 *         matrix ist not visible.
 *
 *  For the concept of index translation matrices, please refer to the
 *  TransparentIndexTranslationMatrix. This class is derived therefrom.
 *  It overloads the apply() and applyAdd() method such that with an
 *  additional method setIndexTranslationMode(), the user can choose
 *  which of the routines apply{Small|Large}to|Small|Large} will be
 *  used (and, correspondingly, which size of arg and dest vectors will
 *  be expected).
 *
 *  Use this class if you want to discard the index translation process,
 *  that means, if you want to make a small matrix look like a larger one
 *  (either in preimage or in image space) via index translation.
 *
 *  You cannot specify other domain and range types that
 *  MatrixType::DomainType and MatrixType::RangeType (in contrast
 *  to TransparentIndexTranslationMatrix) because the apply() and
 *  applyAdd() arguments have to be of the same type for all
 *  smallOrLargeApplyMode's.
 *
 *  Up to this moment, the class is only implemented for symmetric
 *  translation (i. e. same translation on image and preimage).
 *
 *  \author von Deylen
 */

template <typename MatrixType>
class IndexTranslationMatrix
      : public TransparentIndexTranslationMatrix<MatrixType> {

public:
  typedef typename MatrixType::DataType DataType;
  typedef typename MatrixType::DomainType DomainType;
  typedef typename MatrixType::RangeType RangeType;

  enum SmallOrLargeApplyMode {
    SMALL_TO_SMALL,
    SMALL_TO_LARGE,
    LARGE_TO_SMALL,
    LARGE_TO_LARGE
  };

  IndexTranslationMatrix ( int dimSmall, int dimLarge, Vector<int> smallToLarge,
                           typename TransparentIndexTranslationMatrix<MatrixType>::ApplyModeOnUnusedNodes unusedNodesMode
                             = TransparentIndexTranslationMatrix<MatrixType>::APPLY_IDENTITY,
                           SmallOrLargeApplyMode smallOrLargeApplyMode = LARGE_TO_LARGE );

  void apply ( const DomainType & arg, RangeType & dest ) const;
  void applyAdd ( const DomainType & arg, RangeType & dest ) const;

  SmallOrLargeApplyMode  getSmallOrLargeApplyMode  () const;
  void setSmallOrLargeApplyMode  ( SmallOrLargeApplyMode  smallOrLargeApplyMode );

protected:
private:
  SmallOrLargeApplyMode _smallOrLargeApplyMode;
};

// ============================= IMPLEMENTATION =============================

//              2.1 Implementation of class TransparentIndexTranslationMatrix

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
TransparentIndexTranslationMatrix ( int dimSmall, int dimLarge,
                                    const Vector<int> & smallToLarge,
                                    ApplyModeOnUnusedNodes mode ) :
    MatrixType ( dimSmall, dimSmall ),
    _preimSmallToLarge ( smallToLarge ),
    _imSmallToLarge ( smallToLarge ),
    _preimTempVector ( dimSmall ),
    _imTempVector ( dimSmall ),
    _dimLargeRows ( dimLarge ),
    _dimLargeCols ( dimLarge ),
    _unusedNodesMode ( mode ),
    _translationIsSymmetric ( true ) {}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
TransparentIndexTranslationMatrix ( int dimSmallRows, int dimSmallCols,
                                    int dimLargeRows, int dimLargeCols,
                                    const Vector<int> & smallToLargePreimage,
                                    const Vector<int> & smallToLargeImage,
                                    ApplyModeOnUnusedNodes mode ) :
    MatrixType ( dimSmallRows, dimSmallCols ),
    _preimSmallToLarge ( smallToLargePreimage ),
    _imSmallToLarge ( smallToLargeImage ),
    _preimTempVector ( dimSmallCols ),
    _imTempVector ( dimSmallRows ),
    _dimLargeRows ( dimLargeRows ),
    _dimLargeCols ( dimLargeCols ),
    _unusedNodesMode ( mode ) {

  _translationIsSymmetric = isTranslationSymmetric();
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
template<class LargeMatrixType>
void TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
FillFromLarge ( const LargeMatrixType & large ) {
  int n = getNumRowsSmall();
  int m = getNumColsSmall();
  for ( int i = 0; i < n; ++i )
    for ( int j = 0; j < m; ++j )
      this->set ( i, j, large.get ( _imSmallToLarge[i], _preimSmallToLarge[j] ) );
}

//---------------------------------------------------------------------------

// *** VECTOR TRANSLATIONS ***

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
void TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
translateLargeToSmallPreimage ( const LargeDomainType & large, DomainType & small ) const {
  QUOC_ASSERT ( small.size() == this->getNumCols() );
  QUOC_ASSERT ( large.size() == _dimLargeCols );
  for ( int i = 0; i < _preimSmallToLarge.size(); ++i )
    small[i] = large[_preimSmallToLarge[i]];
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
void TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
translateLargeToSmallImage ( const LargeRangeType & large, RangeType & small ) const {
  QUOC_ASSERT ( small.size() == this->getNumRows() );
  QUOC_ASSERT ( large.size() == _dimLargeRows );
  for ( int i = 0; i < _imSmallToLarge.size(); ++i )
    small[i] = large[_imSmallToLarge[i]];
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
void TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
translateSmallToLargePreimage ( const DomainType & small, LargeDomainType & large ) const {
  QUOC_ASSERT ( small.size() == this->getNumCols() );
  QUOC_ASSERT ( large.size() == _dimLargeCols );
  for ( int i = 0; i < _preimSmallToLarge.size(); ++i )
    large[_preimSmallToLarge[i]] = small[i];
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
void TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
translateSmallToLargeImage ( const RangeType & small, LargeRangeType & large ) const {
  QUOC_ASSERT ( small.size() == this->getNumRows() );
  QUOC_ASSERT ( large.size() == _dimLargeRows );
  for ( int i = 0; i < _imSmallToLarge.size(); ++i )
    large[_imSmallToLarge[i]] = small[i];
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
void TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
translateLargeToSmall ( const LargeDomainType & large, DomainType & small ) const {
  QUOC_ASSERT ( _translationIsSymmetric );
  translateLargeToSmallImage ( large, small );
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
void TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
translateSmallToLarge ( const DomainType & large, LargeDomainType & small ) const {
  QUOC_ASSERT ( _translationIsSymmetric );
  translateSmallToLargeImage ( large, small );
}

//---------------------------------------------------------------------------

// *** APPLY ROUTINES ***

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
void TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
applySmallToSmall ( const DomainType & arg, RangeType dest ) const {
  QUOC_ASSERT ( arg .size() == getNumColsSmall() );
  QUOC_ASSERT ( dest.size() == getNumRowsSmall() );

  MatrixType::apply ( arg, dest );
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
void TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
applyLargeToSmall ( const LargeDomainType & arg, RangeType dest ) const {
  QUOC_ASSERT ( arg .size() == getNumColsLarge() );
  QUOC_ASSERT ( dest.size() == getNumRowsSmall() );

  translateLargeToSmallPreimage ( arg, _preimTempVector );
  MatrixType::apply ( _preimTempVector, dest );
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
void TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
applySmallToLarge ( const DomainType & arg, LargeRangeType dest ) const {
  QUOC_ASSERT ( arg .size() == getNumColsSmall() );
  QUOC_ASSERT ( dest.size() == getNumRowsLarge() );

  MatrixType::apply ( arg, _imTempVector );
  translateSmallToLargeImage ( _imTempVector, dest );

  if ( _unusedNodesMode == APPLY_IDENTITY )
    dest.setAllMasked ( 0., getUsedNodesMaskImage(), /* invertMask = */ true );
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
void TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
applyLargeToLarge ( const LargeDomainType & arg, LargeRangeType dest ) const {
  QUOC_ASSERT ( arg .size() == getNumColsLarge() );
  QUOC_ASSERT ( dest.size() == getNumRowsLarge() );

  translateLargeToSmallPreimage ( arg, _preimTempVector );
  MatrixType::apply ( _preimTempVector, _imTempVector );
  translateSmallToLargeImage ( _imTempVector, dest );

  if ( _unusedNodesMode == APPLY_IDENTITY ) {
    QUOC_ASSERT ( _translationIsSymmetric );
    dest.assignMasked ( arg, getUsedNodesMaskPreimage(), /* invertMask = */ true );
  }
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
void TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
applyAddSmallToSmall ( const DomainType & arg, RangeType dest ) const {
  QUOC_ASSERT ( arg .size() == getNumColsSmall() );
  QUOC_ASSERT ( dest.size() == getNumRowsSmall() );

  MatrixType::applyAdd ( arg, dest );
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
void TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
applyAddLargeToSmall ( const LargeDomainType & arg, RangeType dest ) const {
  QUOC_ASSERT ( arg .size() == getNumColsLarge() );
  QUOC_ASSERT ( dest.size() == getNumRowsSmall() );

  translateLargeToSmallPreimage ( arg, _preimTempVector );
  MatrixType::applyAdd ( _preimTempVector, dest );
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
void TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
applyAddSmallToLarge ( const DomainType & arg, LargeRangeType dest ) const {
  QUOC_ASSERT ( arg .size() == getNumColsSmall() );
  QUOC_ASSERT ( dest.size() == getNumRowsLarge() );

  MatrixType::applyAdd ( arg, _imTempVector );
  translateSmallToLargeImage ( _imTempVector, dest );
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
void TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
applyAddLargeToLarge ( const LargeDomainType & arg, LargeRangeType dest ) const {
  QUOC_ASSERT ( arg .size() == getNumColsLarge() );
  QUOC_ASSERT ( dest.size() == getNumRowsLarge() );

  translateLargeToSmallPreimage ( arg, _preimTempVector );
  MatrixType::applyAdd ( _preimTempVector, _imTempVector );
  translateSmallToLargeImage ( _imTempVector, dest );

  if ( _unusedNodesMode == APPLY_IDENTITY ) {
    QUOC_ASSERT ( _translationIsSymmetric );
    dest.addMasked ( arg, getUsedNodesMaskPreimage(), /* invertMask = */ true );
  }
}

//---------------------------------------------------------------------------

// *** GET DIMENSION FUNCTIONS ***

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
int TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
getSizeSmall() const {
  QUOC_ASSERT ( _translationIsSymmetric );
  return getNumColsSmall();
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
int TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
getSizeLarge() const {
  QUOC_ASSERT ( _translationIsSymmetric );
  return getNumColsLarge();
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
const Vector<int> TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
getSmallToLargePreimage() const {
  return _preimSmallToLarge;
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
const Vector<int> TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
getSmallToLargeImage() const {
  return _imSmallToLarge;
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
const Vector<int> TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
getSmallToLarge () const {
  QUOC_ASSERT ( _translationIsSymmetric );
  return getSmallToLargeImage();
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
const typename TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::MaskType &
TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
getUsedNodesMaskPreimage() const {
  if ( !_preimUsedNodesMaskPtr.get() ) {
    _preimUsedNodesMaskPtr.reset ( new MaskType ( getNumColsSmall() ) );
    for ( int i = 0; i < _preimSmallToLarge.size(); ++i )
      _preimUsedNodesMaskPtr->set ( _preimSmallToLarge[i], true );
  }
  return *_preimUsedNodesMaskPtr;
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
const typename TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::MaskType &
TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
getUsedNodesMaskImage() const {
  // for symmetric translations, store only
  // one bit array.
  if ( _translationIsSymmetric )
    return *_preimUsedNodesMaskPtr;

  if ( !_imUsedNodesMaskPtr.get() ) {
    _imUsedNodesMaskPtr.reset ( new MaskType ( getNumColsSmall() ) );
    for ( int i = 0; i < _preimSmallToLarge.size(); ++i )
      _imUsedNodesMaskPtr->set ( _imSmallToLarge[i], true );
  }
  return *_imUsedNodesMaskPtr;
}

//---------------------------------------------------------------------------
template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
typename TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::ApplyModeOnUnusedNodes
TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
getApplyModeOnUnusedNodes() const            {   return _unusedNodesMode;   }
//---------------------------------------------------------------------------
template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
void TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
setApplyModeOnUnusedNodes ( ApplyModeOnUnusedNodes mode )
                                             {   _unusedNodesMode = mode;   }
//---------------------------------------------------------------------------

template <typename MatrixType, typename LargeDomainType, typename LargeRangeType>
bool TransparentIndexTranslationMatrix<MatrixType, LargeDomainType, LargeRangeType>::
isTranslationSymmetric() const {
  bool ret = (    _dimLargeRows         == _dimLargeCols
                  && this->getNumCols() == this->getNumRows()
                  && _preimSmallToLarge    == _imSmallToLarge );
  if ( ret )
    _imUsedNodesMaskPtr.reset();

  return ret;
}

// ==========================================================================

//                         2.2 Implementation of class IndexTranslationMatrix

template <typename MatrixType>
IndexTranslationMatrix<MatrixType>::
IndexTranslationMatrix ( int dimSmall, int dimLarge, Vector<int> smallToLarge,
      typename TransparentIndexTranslationMatrix<MatrixType>::ApplyModeOnUnusedNodes unusedNodesMode,
      SmallOrLargeApplyMode smallOrLargeApplyMode ) :

    TransparentIndexTranslationMatrix<MatrixType> ( dimSmall, dimLarge,
                                                    smallToLarge, unusedNodesMode ),
  _smallOrLargeApplyMode ( smallOrLargeApplyMode ) {}

//---------------------------------------------------------------------------

template <typename MatrixType>
void IndexTranslationMatrix<MatrixType>::
apply ( const DomainType & arg, RangeType & dest ) const {
  switch ( _smallOrLargeApplyMode ) {
  case SMALL_TO_SMALL:
    applySmallToSmall ( arg, dest );
    break;
  case SMALL_TO_LARGE:
    applySmallToLarge ( arg, dest );
    break;
  case LARGE_TO_SMALL:
    applyLargeToSmall ( arg, dest );
    break;
  case LARGE_TO_LARGE:
    applyLargeToLarge ( arg, dest );
    break;
  }
}

//---------------------------------------------------------------------------

template <typename MatrixType>
void IndexTranslationMatrix<MatrixType>::
applyAdd ( const DomainType & arg, RangeType & dest ) const {
  switch ( _smallOrLargeApplyMode ) {
  case SMALL_TO_SMALL:
    applyAddSmallToSmall ( arg, dest );
    break;
  case SMALL_TO_LARGE:
    applyAddSmallToLarge ( arg, dest );
    break;
  case LARGE_TO_SMALL:
    applyAddLargeToSmall ( arg, dest );
    break;
  case LARGE_TO_LARGE:
    applyAddLargeToLarge ( arg, dest );
    break;
  }
}

//---------------------------------------------------------------------------
template <typename MatrixType>
void IndexTranslationMatrix<MatrixType>::
setSmallOrLargeApplyMode  ( SmallOrLargeApplyMode  smallOrLargeApplyMode )
                      {   _smallOrLargeApplyMode = smallOrLargeApplyMode;   }
//---------------------------------------------------------------------------
template <typename MatrixType>
typename IndexTranslationMatrix<MatrixType>::SmallOrLargeApplyMode
IndexTranslationMatrix<MatrixType>::
getSmallOrLargeApplyMode  () const     {   return _smallOrLargeApplyMode;   }
//---------------------------------------------------------------------------

} // end of namespace aol.

#endif
