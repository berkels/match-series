#ifndef __MASKEDOP_H
#define __MASKEDOP_H

#include <op.h>
#include <gridBase.h>

namespace aol {

//   CLASS MaskedOp
/**
 *   This is a simple wrapper-class that calls for each apply(Add) of the
 *   nested operator the apply(Add)Masked with the defined includeWriteMode.
 *   It is needed, because we want to be able to use an apply(Add)-call with
 *   two parameters (e.g. needed in the solvers).
 *
 *   \author Nemitz, von Deylen
 */
template <typename OpType>
class MaskedOp : public OpType {
public:
  typedef typename OpType::DomainType DomainType;
  typedef typename OpType::RangeType  RangeType;
  typedef typename OpType::MaskType   MaskType;

  template <typename OpInitType>
  MaskedOp ( const OpInitType & initializer, const MaskType & mask, const IncludeWriteMode includeWriteMode )
      : OpType ( initializer ), _maskPtr ( &mask ), _includeWriteMode ( includeWriteMode ) {}

  template <typename OpInitType1, typename OpInitType2>
  MaskedOp ( const OpInitType1 & initializer1, const OpInitType2 initializer2, const MaskType & mask, const IncludeWriteMode includeWriteMode )
      : OpType ( initializer1, initializer2 ), _maskPtr ( &mask ), _includeWriteMode ( includeWriteMode ) {}

  aol::IncludeWriteMode getIncludeWriteMode () const {
    return _includeWriteMode;
  }

  //! sets include-write-mode to passed value, returns
  //! include-write-mode that was stored before.
  IncludeWriteMode setIncludeWriteMode ( const IncludeWriteMode includeWriteMode ) {
    IncludeWriteMode ret = getIncludeWriteMode();
    _includeWriteMode = includeWriteMode;
    return ret;
  }

  const MaskType & getMaskReference () const {
    return *_maskPtr;
  }
  const MaskType * getMaskPointer () const {
    return _maskPtr;
  }

  //! sets mask to passed value, returns
  //! pointer to mask that was stored before.
  const MaskType * setMaskReference ( const MaskType & mask ) {
    const MaskType * ret = _maskPtr;
    _maskPtr = &mask;
    return ret;
  }
  //! sets mask to passed value, returns
  //! pointer to mask that was stored before.
  const MaskType * setMaskPointer ( const MaskType * maskPtr ) {
    const MaskType * ret = _maskPtr;
    _maskPtr = maskPtr;
    return ret;
  }

  // call the masked applyAdd method of the base operator
  virtual void applyAdd ( const DomainType &Arg, RangeType &Dest ) const {
    OpType::applyAddMasked ( Arg, Dest, *_maskPtr, _includeWriteMode );
  }

  // call the masked apply method of the base operator
  virtual void apply ( const DomainType &Arg, RangeType &Dest ) const {
    OpType::applyMasked ( Arg, Dest, *_maskPtr, _includeWriteMode );
  }

private:
  //! which nodes should be included and to which should be written?
  const MaskType * _maskPtr;
  IncludeWriteMode _includeWriteMode;
};

} // end of namespace aol.

#endif
