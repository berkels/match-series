#include <l2Projector.h>

template <typename RealType>
aol::ExitStatus qc::Basis2DConst<RealType>::scalarProduct ( qc::BasisFunction2D<RealType> &B, RealType &Prod ) {
  if ( dynamic_cast<qc::Basis2DConst<RealType>* > ( &B ) )  {
    Prod = 1.0;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DX<RealType>* > ( &B ) )  {
    Prod = 0.5;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DY<RealType>* > ( &B ) )  {
    Prod = 0.5;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DXY<RealType>* > ( &B ) ) {
    Prod = 0.25;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DXX<RealType>* > ( &B ) ) {
    Prod = 1. / 3.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DYY<RealType>* > ( &B ) ) {
    Prod = 1. / 3.;
    return aol::SUCCESS;
  }
  return aol::FAILURE;
}

template <typename RealType>
aol::ExitStatus qc::Basis2DX<RealType>::scalarProduct ( qc::BasisFunction2D<RealType> &B, RealType &Prod ) {
  if ( dynamic_cast<qc::Basis2DX<RealType>* > ( &B ) )  {
    Prod = 1. / 3.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DY<RealType>* > ( &B ) )  {
    Prod = 0.25;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DXY<RealType>* > ( &B ) ) {
    Prod = 1. / 6.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DXX<RealType>* > ( &B ) ) {
    Prod = 1. / 4.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DYY<RealType>* > ( &B ) ) {
    Prod = 1. / 6.;
    return aol::SUCCESS;
  }
  return aol::FAILURE;
}

template <typename RealType>
aol::ExitStatus qc::Basis2DY<RealType>::scalarProduct ( qc::BasisFunction2D<RealType> &B, RealType &Prod ) {
  if ( dynamic_cast<qc::Basis2DY<RealType>* > ( &B ) )  {
    Prod = 1. / 3.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DXY<RealType>* > ( &B ) ) {
    Prod = 1. / 6.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DXX<RealType>* > ( &B ) ) {
    Prod = 1. / 6.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DYY<RealType>* > ( &B ) ) {
    Prod = 1. / 4.;
    return aol::SUCCESS;
  }
  return aol::FAILURE;
}

template <typename RealType>
aol::ExitStatus qc::Basis2DXY<RealType>::scalarProduct ( qc::BasisFunction2D<RealType> &B, RealType &Prod ) {
  if ( dynamic_cast<qc::Basis2DXY<RealType>* > ( &B ) )  {
    Prod = 1. / 9.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DXX<RealType>* > ( &B ) ) {
    Prod = 1. / 8.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DYY<RealType>* > ( &B ) ) {
    Prod = 1. / 8.;
    return aol::SUCCESS;
  }
  return aol::FAILURE;
}

template <typename RealType>
aol::ExitStatus qc::Basis2DXX<RealType>::scalarProduct ( qc::BasisFunction2D<RealType> &B, RealType &Prod ) {
  if ( dynamic_cast<qc::Basis2DXX<RealType>* > ( &B ) ) {
    Prod = 1. / 5.;
    return aol::SUCCESS;
  }
  if ( dynamic_cast<qc::Basis2DYY<RealType>* > ( &B ) ) {
    Prod = 1. / 9.;
    return aol::SUCCESS;
  }
  return aol::FAILURE;
}

template <typename RealType>
aol::ExitStatus qc::Basis2DYY<RealType>::scalarProduct ( qc::BasisFunction2D<RealType> &B, RealType &Prod ) {
  if ( dynamic_cast<qc::Basis2DYY<RealType>* > ( &B ) ) {
    Prod = 1. / 5.;
    return aol::SUCCESS;
  }
  return aol::FAILURE;
}

template class qc::Basis2DConst<float>;
template class qc::Basis2DX<float>;
template class qc::Basis2DY<float>;
template class qc::Basis2DXY<float>;
template class qc::Basis2DXX<float>;
template class qc::Basis2DYY<float>;

template class qc::L2Projector<float>;


template class qc::Basis2DConst<double>;
template class qc::Basis2DX<double>;
template class qc::Basis2DY<double>;
template class qc::Basis2DXY<double>;
template class qc::Basis2DXX<double>;
template class qc::Basis2DYY<double>;

template class qc::L2Projector<double>;
