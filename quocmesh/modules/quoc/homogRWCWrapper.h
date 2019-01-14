#ifndef __HOMOGRWCWRAPPER_H
#define __HOMOGRWCWRAPPER_H

#include <aol.h>
#include <quoc.h>
#include <smallMat.h>
#include <geom.h>
#include <gridSize.h>

namespace qc {

template < typename RealType, qc::Dimension Dim >
class HomogRWCMapper;


/** Mapping of 3D "world coordinates" (in the usual quoc understanding that the 3D structure starts at the origin and just fits into the unit cube) to "real world coordinates" lying in a cuboid aligned with the coordinate axes.
 *  This allows for an anisotropic voxel size and an offset, but no rotation.
 *  \author Schwen (MEVIS)
 */
template < typename RealType >
class HomogRWCMapper < RealType, qc::QC_3D > {
protected:
  aol::Matrix44<RealType> _toInternalConvertMat, _toWorldConvertMat;
  const qc::GridSize<qc::QC_3D> _gridSize;

public:
  // default constructor, default destructor, copy constructor, assignment operator OK

  HomogRWCMapper ( const aol::Vec3<RealType> &Offset, const aol::Vec3<RealType> &Resolutions, const qc::GridSize<qc::QC_3D> &GridSize );

  void toInternal ( const aol::Vec3<RealType> &Arg, aol::Vec3<RealType> &Dest ) const;

  void toRWC ( const aol::Vec3<RealType> &Arg, aol::Vec3<RealType> &Dest ) const;

  const aol::Matrix44<RealType>& getToWorldConvertMat ( ) const {
    return ( _toWorldConvertMat );
  }

  const aol::Matrix44<RealType>& getToInternalConvertMat ( ) const {
    return ( _toInternalConvertMat );
  }

  RealType getVoxelVolume ( ) const;

};


template< typename WrappedType, typename WrappedContainedType, typename RealType, qc::Dimension Dim >
class HomogRWCWrapper;


/** Wrapping of 3D data structures so that they are interpreted as being shifted and (anisotropically) scaled (but not rotated).
 *  \author Schwen (MEVIS)
 */
template < typename WrappedType, typename WrappedContainedType, typename RealType >
class HomogRWCWrapper < WrappedType, WrappedContainedType, RealType, qc::QC_3D > {
protected:
  WrappedType& _wrapped;
  const HomogRWCMapper<RealType, qc::QC_3D> _mapper;

private:
  HomogRWCWrapper(); // do not implement
  HomogRWCWrapper& operator= ( const HomogRWCWrapper < WrappedType, WrappedContainedType, RealType, qc::QC_3D > &other ); // do not implement

public:
  // default destructor OK
  // default copy constructor OK

  HomogRWCWrapper ( WrappedType &Wrapped, const aol::Vec3<RealType> &Offset, const aol::Vec3<RealType> &Resolutions ) : _wrapped ( Wrapped ), _mapper ( Offset, Resolutions, Wrapped.getSize() ) {
  }

  HomogRWCWrapper ( WrappedType &Wrapped, const HomogRWCMapper<RealType, qc::QC_3D> &Mapper ) : _wrapped ( Wrapped ), _mapper ( Mapper ) {
  }

  //! get value at position given in real world coordinates (notice the rounding necessary, this is mainly useful if interpolation does not make sense)
  WrappedContainedType getAtRWC ( const aol::Vec3<RealType> &RWCPos ) const {
    return ( _wrapped.get ( this->getIntPos ( RWCPos ) ) );
  }

  const WrappedContainedType& getRefAtRWC ( const aol::Vec3<RealType> &RWCPos ) const {
    return ( _wrapped.getRef ( this->getIntPos ( RWCPos ) ) );
  }

  WrappedContainedType& getRefAtRWC ( const aol::Vec3<RealType> &RWCPos ) {
    return ( _wrapped.getRef ( this->getIntPos ( RWCPos ) ) );
  }

  //! get value at position given in real world coordinates, returning default value for positions outside
  WrappedContainedType getAtRWCDefaultOutside ( const aol::Vec3<RealType> &RWCPos, const WrappedContainedType DefaultValueOutside ) const {
    aol::Vec3<RealType> intPosR;
    _mapper.toInternal ( RWCPos, intPosR );

    const aol::Vec3<int> intPos ( static_cast<int> ( intPosR[0] ), static_cast<int> ( intPosR[1] ), static_cast<int> ( intPosR[2] ) );

    if ( this->checkBounds ( intPos[0], intPos[1], intPos[2] ) ) {
      return ( _wrapped.get ( intPos ) );
    } else {
      return ( DefaultValueOutside );
    }
  }

  //! interpolate value at position given in real world coordinates
  WrappedContainedType interpolateAtRWC ( const aol::Vec3<RealType> &RWCPos ) const {
    aol::Vec3<RealType> intPosR;
    _mapper.toInternal ( RWCPos, intPosR );
    this->boundsCheck ( intPosR[0], intPosR[1], intPosR[2], "HomogRWCWrapper::interpolateAtRWC" );
    return ( _wrapped.interpolate ( intPosR ) );
  }

  //! interpolate value at position given in real world coordinates, returning default value for positions outside
  WrappedContainedType interpolateAtRWCDefaultOutside ( const aol::Vec3<RealType> &RWCPos, const WrappedContainedType DefaultValueOutside ) const {
    aol::Vec3<RealType> intPosR;
    _mapper.toInternal ( RWCPos, intPosR );
    if ( this->checkBounds ( intPosR[0], intPosR[1], intPosR[2] ) ) {
      return ( _wrapped.interpolate ( intPosR ) );
    } else {
      return ( DefaultValueOutside );
    }
  }

  //! compute the bounding box of the represented volume
  RealType computeBoundingBoxRWC ( aol::Vec3<RealType> &boxLower, aol::Vec3<RealType> &boxUpper ) const {
    aol::Vec3<RealType> intLower ( 0, 0, 0 ), intUpper ( _wrapped.getNumX(), _wrapped.getNumY(), _wrapped.getNumZ() );
    _mapper.toRWC ( intLower, boxLower );
    _mapper.toRWC ( intUpper, boxUpper );
    return ( ( boxUpper - boxLower ).getMaxValue() );
  }

  const aol::Matrix44<RealType>& getToWorldConvertMat ( ) const {
    return ( _mapper.getToWorldConvertMat() );
  }

  const HomogRWCMapper<RealType, qc::QC_3D> & getMapperRef ( ) const {
    return ( _mapper );
  }

  inline void convertCoordsToInternal ( const aol::Vec3<RealType> &Arg, aol::Vec3<RealType> &Dest ) const {
    _mapper.toInternal ( Arg, Dest );
  }

  inline void convertCoordsToRW ( const aol::Vec3<RealType> &Arg, aol::Vec3<RealType> &Dest ) const {
    _mapper.toRWC ( Arg, Dest );
  }

  aol::Vec3<RealType> getConvertedCoordsToRW ( const qc::CoordType &Arg ) const {
    aol::Vec3<RealType> ArgR ( Arg[0], Arg[1], Arg[2] ), dest;
    this->convertCoordsToRW ( ArgR, dest );
    return ( dest );
  }

  const WrappedType& getWrappedRef ( ) const {
    return ( _wrapped );
  }

  aol::Vec3<int> getSize() const {
    return ( _wrapped.getSize() );
  }

  RealType getVoxelVolume ( ) const {
    return ( _mapper.getVoxelVolume() );
  }


private:
  inline aol::Vec3<int> getIntPos ( const aol::Vec3<RealType> &RWCPos ) {
    aol::Vec3<RealType> intPosR;
    _mapper.toInternal ( RWCPos, intPosR );

    const aol::Vec3<int> ret ( static_cast<int> ( intPosR[0] ), static_cast<int> ( intPosR[1] ), static_cast<int> ( intPosR[2] ) );
    this->boundsCheck ( ret[0], ret[1], ret[2], "HomogRWCWrapper::getIntPos" );
    return ( ret );
  }

  template< typename Num >
  inline bool checkBounds ( const Num ix, const Num iy, const Num iz ) const {
    return ( ( ix >= 0 ) && ( iy >= 0 ) && ( iz >= 0 ) && ( ix <= _wrapped.getNumX() - 1 ) && ( iy <= _wrapped.getNumY() - 1 ) && ( iz <= _wrapped.getNumZ() - 1 ) );
    // e.g. _wrapped.getNumX() - 0.5 is already outside
  }

  inline void boundsCheck ( const int ix, const int iy, const int iz, const char* where ) const {
    if ( ! ( checkBounds ( ix, iy, iz ) ) ) {
      char errmsg[1024];
      sprintf( errmsg, "qc::HomogRWCWrapper:%s %d %d %d are out of bounds (%d %d %d)!", where, ix, iy, iz, _wrapped.getNumX(), _wrapped.getNumY(), _wrapped.getNumZ() );
      cerr << errmsg << endl;
      throw aol::OutOfBoundsException( errmsg, __FILE__, __LINE__ );
    }
  }

  inline void boundsCheck ( const RealType ix, const RealType iy, const RealType iz, const char* where ) const {
    if ( ! ( checkBounds ( ix, iy, iz ) ) ) {
      char errmsg[1024];
      sprintf( errmsg, "qc::HomogRWCWrapper:%s %f %f %f are out of bounds (%d %d %d)!", where, static_cast<double>( ix ), static_cast<double>( iy ), static_cast<double>( iz ), _wrapped.getNumX(), _wrapped.getNumY(), _wrapped.getNumZ() );
      cerr << errmsg << endl;
      throw aol::OutOfBoundsException( errmsg, __FILE__, __LINE__ );
    }
  }

};

}

#endif
