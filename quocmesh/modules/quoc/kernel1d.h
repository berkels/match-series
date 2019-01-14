#ifndef __KERNEL1D_H
#define __KERNEL1D_H

#include <quoc.h>
#include <vec.h>

namespace qc {

//! Abstract class for the implementation of an arbitrary 3d filtering kernel.
template<typename RealType = float>
class Kernel1d : public aol::Vector<RealType> {
public:
  Kernel1d ( int Size ) : aol::Vector<RealType> ( Size ),
      size ( Size ) {
    offset = ( size >> 1 );
    QUOC_ASSERT ( ( ( offset << 1 ) + 1 ) == size );
  }

  virtual ~Kernel1d() { }

  RealType getValue ( int OffX ) const {
    return this->_pData[ OffX + offset ];
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

  Kernel1d<RealType>& operator= ( const Kernel1d<RealType>& rhs ) {
    // Beware of self-assignment
    if ( this == &rhs ) return *this;
    aol::Vector<RealType>::operator = ( rhs );
    size = rhs.size;
    offset = rhs.offset;
    return *this;
  }

protected:

  void setValue ( int OffX, RealType Val ) {
    this->_pData[ OffX + offset ] = Val;
  }

  void normalize() {
    ( *this ) /= this->sum();
  }

  int size, offset;
};

//! A 1D Gaussian Filter Kernel.

template<typename RealType>
class GaussKernel1d : public qc::Kernel1d<RealType> {
public:
  GaussKernel1d ( int Size, RealType Sigma ) :
      Kernel1d<RealType> ( Size ), sigma ( Sigma ) {
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

//! A 1D Gaussian Differentiation Filter Kernel.
template<typename RealType>
class GaussDiffKernel1d : public qc::Kernel1d<RealType> {
public:

  GaussDiffKernel1d ( int Size, RealType Sigma, DiffVarType Comp ) :
      Kernel1d<RealType> ( Size ), sigma ( Sigma ), comp ( Comp ) {
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

}

#endif
