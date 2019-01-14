#ifndef __KERNEL2D_H
#define __KERNEL2D_H

#include <quoc.h>
#include <vec.h>

namespace qc {

//! Abstract class for the implementation of an arbitrary 2d filtering kernel.
template<typename RealType = float>
class Kernel2d : public aol::Vector<RealType> {
public:
  Kernel2d ( int Size ) : aol::Vector<RealType> ( Size * Size ),
      size ( Size ) {
    offset = ( size >> 1 );
    QUOC_ASSERT ( ( ( offset << 1 ) + 1 ) == size );
  }

  virtual ~Kernel2d() { }

  RealType getValue ( int OffX, int OffY ) const {
    return this->_pData[ ( OffY + offset ) *size + OffX + offset ];
  }

  int getSize() const {
    return size;
  }

  int getOffset() const {
    return offset;
  }

  //! Abstract virtual function for the initialization of the kernel.
  virtual void makeKernel() = 0;

  void dump() const;

  Kernel2d<RealType>& operator= ( const Kernel2d<RealType>& rhs ) {
    // Beware of self-assignment
    if ( this == &rhs ) return *this;
    aol::Vector<RealType>::operator = ( rhs );
    size = rhs.size;
    offset = rhs.offset;
    return *this;
  }

protected:

  void setValue ( int OffX, int OffY, RealType Val ) {
    this->_pData[ ( OffY + offset ) *size + OffX + offset ] = Val;
  }

  void normalize() {
    ( *this ) /= this->sum();
  }

  int size, offset;
};

//! A 2D Gaussian Filter Kernel.

template<typename RealType>
class GaussKernel2d : public qc::Kernel2d<RealType> {
public:
  GaussKernel2d ( int Size, RealType Sigma ) :
      Kernel2d<RealType> ( Size ), sigma ( Sigma ) {
    makeKernel();
  }

  virtual void makeKernel();

  void setSigma ( RealType Sigma ) {
    sigma = Sigma;
    makeKernel();
  }

  RealType getSigma() const {
    return sigma;
  }

protected:
  RealType sigma;
};

//! A 2D Gaussian Differentiation Filter Kernel.
template<typename RealType>
class GaussDiffKernel2d : public qc::Kernel2d<RealType> {
public:

  GaussDiffKernel2d ( int Size, RealType Sigma, DiffVarType Comp ) :
      Kernel2d<RealType> ( Size ), sigma ( Sigma ), comp ( Comp ) {
    makeKernel();
  }

  virtual void makeKernel();

  void setSigma ( RealType Sigma ) {
    sigma = Sigma;
    makeKernel();
  }

  RealType getSigma() const {
    return sigma;
  }

protected:
  RealType sigma;
  DiffVarType comp;
};

//! Special discrete approximation of the Gauss kernel based on Pascal's triangle.
//! \author Berkels
template<typename RealType>
class DiscreteGaussKernel2d : public qc::Kernel2d<RealType> {
public:
  DiscreteGaussKernel2d ( int Size ) : Kernel2d<RealType> ( Size ) {
    makeKernel();
  }

  virtual void makeKernel();
};

//! 2D circle kernel used for the moving average filter.
//! \author Berkels
template<typename RealType>
class CircleAverageKernel2d : public qc::Kernel2d<RealType> {
public:
  CircleAverageKernel2d ( int Size ) : Kernel2d<RealType> ( Size ) {
    makeKernel();
  }

  virtual void makeKernel();
};

//! Blank 2D kernel that can be filled with setValue.
//! \author Berkels
template<typename RealType>
class ModifiableKernel2d : public qc::Kernel2d<RealType> {
public:
  ModifiableKernel2d ( int Size ) : qc::Kernel2d<RealType> ( Size ) {}

  virtual void makeKernel() {}

  using qc::Kernel2d<RealType>::setValue;
  using aol::Vector<RealType>::sum;
};

}

#endif
