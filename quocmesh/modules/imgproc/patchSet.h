#ifndef __PATCHSET_H
#define __PATCHSET_H

#include <multiArray.h>
#include <convolution.h>

namespace im {

/**
 * \author Berkels
 */
template <typename RealType, const int NumComponents>
void extractPatchesFrom2DArray ( const qc::MultiArray<RealType, qc::QC_2D, NumComponents> &InMArray, const int PatchWidth, aol::MultiVector<RealType> &Patches,
                                 aol::Vec2<int> *PNumPatches = NULL,
                                 const bool FullOverlap = false ) { 
  const int entrySize = NumComponents * aol::Sqr ( PatchWidth );
  const int patchesNumX = FullOverlap ? ( InMArray[0].getNumX() - PatchWidth + 1 ) : ( InMArray[0].getNumX() / PatchWidth );
  const int patchesNumY = FullOverlap ? ( InMArray[0].getNumY() - PatchWidth + 1 ) : ( InMArray[0].getNumY() / PatchWidth );
  const int patchOffset = FullOverlap ? 1 : PatchWidth;
  Patches.reallocate ( patchesNumX*patchesNumY, entrySize );
  if ( PNumPatches != NULL ) {
    (*PNumPatches)[0] = patchesNumX;
    (*PNumPatches)[1] = patchesNumY;
  }

  for ( int yPatch = 0; yPatch < patchesNumY; ++yPatch ) {
    for ( int xPatch = 0; xPatch < patchesNumX; ++xPatch ) {
      for ( int comp = 0; comp < NumComponents; ++comp ) {
        qc::ScalarArray<RealType, qc::QC_2D> temp ( PatchWidth, PatchWidth, Patches[ qc::ILexCombine2 ( xPatch, yPatch, patchesNumX ) ].getData() + comp * ( entrySize / NumComponents ), aol::FLAT_COPY );
        InMArray[comp].copyBlockTo ( xPatch * patchOffset, yPatch * patchOffset, temp );
      }
    }
  }
}

/**
 * \author Berkels
 */
template <typename RealType>
void extractPatchesFrom2DArray ( const qc::ScalarArray<RealType, qc::QC_2D> &InArray, const int PatchWidth, aol::MultiVector<RealType> &Patches,
                                 aol::Vec2<int> *PNumPatches = NULL,
                                 const bool FullOverlap = false ) {
  const qc::MultiArray<RealType, qc::QC_2D, 1> mArray ( InArray, aol::FLAT_COPY );
  extractPatchesFrom2DArray ( mArray, PatchWidth, Patches, PNumPatches, FullOverlap );
}

/**
 * \author Berkels
 */
template <typename RealType, const int NumComponents>
void build2DArrayFromPatches ( const aol::MultiVector<RealType> &Patches,
                               const int PatchWidth,
                               const int PatchesNumX,
                               const int PatchesNumY,
                               qc::MultiArray<RealType, qc::QC_2D, NumComponents> &OutMArray,
                               const int Padding = 0,
                               const RealType PadValue = 0.5,
                               const bool ReallocateArray = true,
                               const RealType Lambda = 0,
                               const qc::MultiArray<RealType, qc::QC_2D, NumComponents> *InMArray = NULL ) { 
  const int entryCompSize = Patches[0].size() / NumComponents;

  if ( ReallocateArray )
    OutMArray.reallocate ( PatchesNumX*PatchWidth+Padding*(PatchesNumX-1), PatchesNumY*PatchWidth+Padding*(PatchesNumY-1) );

  if ( Padding > 0 )
    OutMArray.setAll ( PadValue );

  qc::ScalarArray<RealType, qc::QC_2D> oneArray;
  qc::ScalarArray<RealType, qc::QC_2D> countArray;
  if ( Padding < 0 ) {
    oneArray.reallocate ( PatchWidth, PatchWidth );
    oneArray.setAll ( 1 );
    countArray.reallocate ( OutMArray[0].getNumX(), OutMArray[0].getNumY() );
    if ( ReallocateArray == false )
      OutMArray.setZero();
  }

  for ( int yPatch = 0; yPatch < PatchesNumY; ++yPatch ) {
    for ( int xPatch = 0; xPatch < PatchesNumX; ++xPatch ) {
      const int xPos = xPatch * ( PatchWidth + Padding );
      const int yPos = yPatch * ( PatchWidth + Padding );
      for ( int comp = 0; comp < NumComponents; ++comp ) {
        qc::ScalarArray<RealType, qc::QC_2D> temp ( PatchWidth, PatchWidth, Patches[ qc::ILexCombine2 ( xPatch, yPatch, PatchesNumX ) ].getData() + comp * entryCompSize, aol::FLAT_COPY );
        if ( Padding < 0 )
          OutMArray[comp].pasteAddFrom ( temp, xPos, yPos );
        else
          OutMArray[comp].pasteFrom ( temp, xPos, yPos );
      }
      if ( Padding < 0 )
        countArray.pasteAddFrom ( oneArray, xPos, yPos );
    }
  }

  if ( Padding < 0 ) {
    if ( ( Lambda > 0 ) && ( InMArray != 0 ) ) {
      OutMArray.addMultiple ( *InMArray, Lambda );
      countArray.addToAll ( Lambda );
    }

    for ( int comp = 0; comp < NumComponents; ++comp ) {
      OutMArray[comp] /= countArray;
    }
  }
}

/**
 * \author Berkels
 */
template <typename RealType>
void build2DArrayFromPatches ( const aol::MultiVector<RealType> &Patches,
                               const int PatchWidth,
                               const int PatchesNumX,
                               const int PatchesNumY,
                               qc::ScalarArray<RealType, qc::QC_2D> &OutArray,
                               const int Padding = 0,
                               const RealType PadValue = 0.5 ) { 
  // We can't call reallocate on flat copies, so we need to reallocate OutArray here.
  OutArray.reallocate ( PatchesNumX*PatchWidth+Padding*(PatchesNumX-1), PatchesNumY*PatchWidth+Padding*(PatchesNumY-1) );
  qc::MultiArray<RealType, qc::QC_2D, 1> mArray ( OutArray, aol::FLAT_COPY );
  build2DArrayFromPatches ( Patches, PatchWidth, PatchesNumX, PatchesNumY, mArray, Padding, PadValue, false );
}

/**
 * \author Berkels
 */
template <typename RealType, int NumComponents>
void buildAndSavePNGFromPatches ( const aol::MultiVector<RealType> &Patches,
                                  const int PatchWidth,
                                  const int PatchesNumX,
                                  const int PatchesNumY,
                                  const char *FileName,
                                  const int Padding = 0,
                                  const RealType PadValue = 0.5 ) { 
  qc::MultiArray<RealType, qc::QC_2D, NumComponents> outMArray;
  build2DArrayFromPatches<RealType, NumComponents> ( Patches, PatchWidth, PatchesNumX, PatchesNumY, outMArray, Padding, PadValue );
  outMArray.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
  outMArray.savePNG ( FileName );
}

/**
 * \author Berkels
 */
template <typename RealType, int NumComponents>
class PatchSet2D : protected aol::MultiVector<RealType> {
  const int _patchWidth;
  int _overlap;
  aol::Vec2<int> _size;

  PatchSet2D () {}

public:
  PatchSet2D ( const int PatchWidth )
    : aol::MultiVector<RealType> ( 0, 0 ),
      _patchWidth ( PatchWidth ),
      _overlap ( 0 ),
      _size () {}

  PatchSet2D ( const PatchSet2D<RealType, NumComponents> &Other, aol::CopyFlag copyFlag = aol::DEEP_COPY )
    : aol::MultiVector<RealType> ( Other, copyFlag ),
      _patchWidth ( Other.getPatchWidth() ),
      _overlap ( Other.getOverlap() ),
      _size ( Other.getSize() ) {}

  //! Creates patches filled with zero.
  PatchSet2D ( const int PatchWidth, const int Overlap, const aol::Vec2<int> &Size )
    : aol::MultiVector<RealType> ( Size[0]*Size[1], NumComponents * aol::Sqr ( PatchWidth ) ),
      _patchWidth ( PatchWidth ),
      _overlap ( Overlap ),
      _size ( Size ) {}

  PatchSet2D ( const int PatchWidth, const int Overlap, const aol::Vec2<int> &Size, const aol::MultiVector<RealType> &Patches, aol::CopyFlag copyFlag = aol::FLAT_COPY )
    : aol::MultiVector<RealType> ( Patches, copyFlag ),
      _patchWidth ( PatchWidth ),
      _overlap ( Overlap ),
      _size ( Size ) {}

  void clear ( ) {
    this->reallocate ( 0, 0 );
    _size.setZero();
    _overlap = 0;
  }

  void extractPatchesFromArray ( const qc::MultiArray<RealType, qc::QC_2D, NumComponents> &InMArray, const bool FullOverlap ) {
    if ( getNumPatches() > 0 )
      this->reallocate ( 0, 0 );

    im::extractPatchesFrom2DArray<RealType, NumComponents> ( InMArray, _patchWidth, *this, &_size, FullOverlap );
    _overlap = FullOverlap ? ( 1 - _patchWidth) : 0;
  }

  void extractPatchesFromArray ( const qc::ScalarArray<RealType, qc::QC_2D> &InArray, const bool FullOverlap ) {
    const qc::MultiArray<RealType, qc::QC_2D, NumComponents> inMArray ( InArray, aol::FLAT_COPY );
    extractPatchesFromArray ( inMArray, FullOverlap );
  }

  void extractPatchesFromArray ( const char* InImageFilename, const bool FullOverlap ) {
    qc::MultiArray<RealType, qc::QC_2D, NumComponents> mArray;

    if ( aol::fileNameEndsWith ( InImageFilename, ".png" ) ) {
      mArray.loadPNG ( InImageFilename );
    } else {
      if ( NumComponents != 1 )
        throw aol::UnimplementedCodeException ( "PatchSet2D::extractPatchesFromArray: Loading of non-png files only implemented for NumComponents == 1.", __FILE__, __LINE__ );

      mArray.reallocate ( qc::GridSize<qc::QC_2D> ( qc::getSizeFromArrayFile ( InImageFilename ) ) );
      mArray.load ( InImageFilename );
    }
    mArray /= mArray.getMaxValue();
    extractPatchesFromArray ( mArray, FullOverlap );
  }

  void multPatchesWithGaussian ( const RealType Sigma ) {
    qc::ScalarArray<RealType, qc::QC_2D> gaussian ( _patchWidth, _patchWidth );
    qc::generateGaussKernel<RealType> ( Sigma, gaussian, false, false );

    const int entryCompSize = aol::Sqr ( _patchWidth );

    for ( int i = 0; i < getNumPatches(); ++i ) {
      for ( int comp = 0; comp < NumComponents; ++comp ) {
        aol::Vector<RealType> temp ( (*this)[i].getData() + comp * entryCompSize, entryCompSize, aol::FLAT_COPY );
        temp.entrywiseMultiplyFrom ( gaussian );
      }
    }
  }

  void saveAsPNG ( const char* Filename ) const {
    im::buildAndSavePNGFromPatches<RealType, NumComponents> ( *this, _patchWidth, _size[0], _size[1], Filename, _overlap );
  }

  const aol::MultiVector<RealType>& getMultiVectorReference ( ) const {
    return *this;
  }

  aol::MultiVector<RealType>& getMultiVectorReference ( ) {
    return *this;
  }

  const aol::Vector<RealType>& getPatch ( const int XIndex, const int YIndex ) const {
    return (*this)[ qc::ILexCombine2 ( XIndex, YIndex, _size[0] ) ];
  }

  const aol::Vector<RealType>& getPatch ( const aol::Vec<2, short> &Coord ) const {
    return getPatch ( Coord[0], Coord[1] );
  }

  int getPatchWidth ( ) const {
    return _patchWidth;
  }

  const aol::Vec2<int>& getSize ( ) const {
    return _size;
  }

  int getEntrySize ( ) const {
    return getPatch ( 0, 0 ).size();
  }

  //! If this patch set was created from an array with extractPatchesFromArray and ( FullOverlap == true ),  
  // a patch at position (X,Y) corresponds to the pixel (X+getXYOffset(),Y+getXYOffset) in the original array.
  int getXYOffset ( ) const {
    return static_cast<int> ( ceil ( getPatchWidth() / 2. ) );
  }

  int getOverlap ( ) const {
    return _overlap;
  }

  int getNumPatches() const {
   return this->numComponents();
  }
};

} // namespace qc

#endif // __PATCHSET_H
