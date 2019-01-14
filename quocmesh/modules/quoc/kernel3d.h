#ifndef __KERNEL3D_H
#define __KERNEL3D_H

#include <vec.h>
#include <quoc.h>
#include <kernel2d.h>

namespace qc {

//! Abstract class for an arbitrary convolution kernel.
template<typename T_RealType = float>
class Kernel3d : public aol::Vector<T_RealType> {
public:
  Kernel3d ( int Size ) :
      aol::Vector<T_RealType> ( Size * Size * Size ),
      size ( Size ) {
    if ( size % 2 == 0 ) throw aol::Exception ( "Kernel must be odd", __FILE__, __LINE__ );
    offset = ( size >> 1 );
  }

virtual ~Kernel3d() { }

  T_RealType getValue ( int OffX, int OffY, int OffZ ) const {
    return this->_pData[ ( OffZ + offset ) * size * size + ( OffY + offset ) *size + OffX + offset ];
  }

  int getSize() const {
    return size;
  }

  virtual void makeKernel() = 0;

protected:

  void setValue ( int OffX, int OffY, int OffZ, T_RealType Val ) {
    this->_pData[ ( OffZ + offset ) *size*size + ( OffY + offset ) *size + OffX + offset ] = Val;
  }

  //! Normalize the kernel such that all values add up to 1.
  void normalize() {
    ( *this ) /= this->sum();
  }

  int size, offset;
};

//! A Gaussian filtering kernel.
/*! This class is the implementation of a filtering kernel of the
  form
  \f[ F(x)=\frac{1}{(2 \pi \sigma^2)^\frac{3}{2}}e^{-\frac{x^2}{2\sigma^2}} \f]
*/
template<typename T_RealType>
class GaussKernel3d : public qc::Kernel3d<T_RealType> {
public:
  GaussKernel3d ( int Size, T_RealType Sigma ) :
      Kernel3d<T_RealType> ( Size ), sigma ( Sigma ) {
    makeKernel();
  }

  //! Make this kernel a Gaussian.
  virtual void makeKernel();

  //! Set the smoothing parameter \f$\sigma=\sqrt{2t}\f$
  void setSigma ( T_RealType Sigma ) {
    sigma = Sigma;
    makeKernel();
  }

  //! Returns \f$\sigma\f$
  T_RealType getSigma() const {
    return sigma;
  }

protected:
  T_RealType sigma;
};

//! Kernel for Derivatives on Gaussian smoothed data.
template<typename T_RealType>
class GaussDiffKernel3d : public qc::Kernel3d<T_RealType> {
public:

  GaussDiffKernel3d ( int Size, T_RealType Sigma, DiffVarType Comp ) :
      Kernel3d<T_RealType> ( Size ), sigma ( Sigma ), comp ( Comp ) {
    makeKernel();
  }

  virtual void makeKernel();

  void setSigma ( T_RealType Sigma ) {
    sigma = Sigma;
    makeKernel();
  }

  T_RealType getSigma() const {
    return sigma;
  }

protected:
  T_RealType sigma;
  DiffVarType comp;

};

template <class RealType, int dim>
class KernelTrait {};

template <class RealType>
class KernelTrait<RealType, 2> {
public:
  typedef qc::Kernel2d<RealType> KernelType;
};

template <class RealType>
class KernelTrait<RealType, 3> {
public:
  typedef qc::Kernel3d<RealType> KernelType;
};

template <class RealType, int Dim>
class AveragingKernel :
  public KernelTrait<RealType, Dim>::KernelType {
public:
  AveragingKernel( int Size ) :
    KernelTrait<RealType, Dim>::KernelType( Size ) {
    makeKernel();
  }
  void makeKernel() {
    this->setAll( 1 );
    this->normalize();
  }
};

}

#endif
