#ifndef __GRIDSIZE_H
#define __GRIDSIZE_H

#include <aol.h>
#include <quoc.h>

// --------------------------------------------------------------------------

namespace qc {

template <Dimension _Dim>
class GridSize : public CoordType {
public:
  typedef CoordType::DataType IndexType;

  static const Dimension Dim = _Dim;

  // ***** Constructors *****

  //! constructs GridSize with equal widths in all _Dim directions.
  explicit GridSize ( IndexType NumXYZ );
  //! Constructor for 2d. Not implemented for 3d,
  //! thus raising a linker error there
  GridSize ( IndexType NumX, IndexType NumY );
  //! Constructor for 3d. Not implemented for 2d,
  //! thus raising a linker error there
  GridSize ( IndexType NumX, IndexType NumY, IndexType NumZ );

  //! Copy constructur. Needs to be implemented, the default one won't work since it'll possibly call :operator [] (2)
  //! on the argument even if (_Dim == QC_2D). This happens when calling the default copy constructor of the parent class CoordType.
  GridSize ( const GridSize<_Dim> &Size );

  template <typename IntType> GridSize ( const aol::Vec2<IntType> & Size );
  template <typename IntType> GridSize ( const aol::Vec3<IntType> & Size );
  template <int d, typename IntType>
  GridSize ( const aol::Vec<d, IntType> & Size );

  // creates a GridSize object from getNum{X|Y|Z}
  template <typename InitType>
  explicit GridSize ( const InitType & initObj );

  // creates a GridSize object from getNum{X|Y|Z}
  template <typename InitType>
  static GridSize<_Dim> createFrom ( const InitType & initObj );

  IndexType getNumX() const;
  IndexType getNumY() const;
  IndexType getNumZ() const;

  int getNumberOfNodes () const;

  // number of Dofs is the same as number of nodes
  // this is introduced to get compatibility with adaptive grids
  int getNumberOfDofs () const;

  typename aol::VecDimTrait<int,_Dim>::VecType getSizeAsVecDim ( ) const;

  IndexType operator [] (int i) const;
  IndexType & operator [] (int i);

  void quadraticOrDie () const;
  
  void setAll ( IndexType NumXYZ );
};

typedef GridSize<QC_2D> GridSize2d;
typedef GridSize<QC_3D> GridSize3d;

// --------------------------------------------------------------------------
template <>
template <typename InitType>
GridSize<QC_1D> GridSize<QC_1D>::createFrom ( const InitType & initObj ) {
  return GridSize<QC_1D> ( static_cast<short> ( initObj.getNumX() ) );
}
// --------------------------------------------------------------------------
template <>
template <typename InitType>
GridSize<QC_2D> GridSize<QC_2D>::createFrom ( const InitType & initObj ) {
  return GridSize<QC_2D> ( initObj.getNumX(), initObj.getNumY() );
}
// --------------------------------------------------------------------------
template <>
template <typename InitType>
GridSize<QC_3D> GridSize<QC_3D>::createFrom ( const InitType & initObj ) {
  return GridSize<QC_3D> ( initObj.getNumX(), initObj.getNumY(),
                           initObj.getNumZ() );
}

// --------------------------------------------------------------------------
template <>
template <typename InitType>
GridSize<QC_1D>::GridSize ( const InitType & initObj ) {
  (*this)[0] = initObj.getNumX();
}
// --------------------------------------------------------------------------
template <>
template <typename InitType>
GridSize<QC_2D>::GridSize ( const InitType & initObj ) {
  (*this)[0] = initObj.getNumX();
  (*this)[1] = initObj.getNumY();
}
// --------------------------------------------------------------------------
template <>
template <typename InitType>
GridSize<QC_3D>::GridSize ( const InitType & initObj ) {
  (*this)[0] = initObj.getNumX();
  (*this)[1] = initObj.getNumY();
  (*this)[2] = initObj.getNumZ();
}
// --------------------------------------------------------------------------
template <Dimension _Dim>
GridSize<_Dim>::GridSize ( typename GridSize<_Dim>::IndexType NumXYZ ) {
  for (int i = 0; i < Dim; ++i)
    (*this)[i] = NumXYZ;
  for (int i = Dim; i < dim; ++i)
    CoordType::operator [] (i) = 1;
}
// --------------------------------------------------------------------------
template <Dimension _Dim>
template <typename IntType>
GridSize<_Dim>::GridSize ( const aol::Vec2<IntType> & Size ) {
  if ( Dim > QC_2D )
    throw aol::DimensionMismatchException ( "GridSize<d> constructor called "
              "with d > 2, but only a Vec2 was given.", __FILE__, __LINE__ );

  (*this)[0] = static_cast<DataType> ( Size[0] );
  (*this)[1] = static_cast<DataType> ( Size[1] );
  CoordType::operator [] (2) = 1;
}
// --------------------------------------------------------------------------
template <Dimension _Dim>
template <typename IntType>
GridSize<_Dim>::GridSize ( const aol::Vec3<IntType> & Size ) {
  // initialize ourselves and assert that additional
  // components of Size are all one:
  for (int i = 0; i < Dim; ++i)
    (*this)[i] = static_cast<DataType> ( Size[i] );
  for (int i = Dim; i < 3; ++i)
    if ( Size[i] != 1 )
      throw aol::DimensionMismatchException ( aol::strprintf( "GridSize"
                     "<QC_%iD> constructor for Vec3 Size with Size[%i] = %i",
                     Dim, i, Size[i] ), __FILE__, __LINE__ );
    else
      CoordType::operator [] (i) = 1;
}
// --------------------------------------------------------------------------
template <Dimension _Dim>
template <int d, typename IntType>
GridSize<_Dim>::GridSize ( const aol::Vec<d, IntType> & Size ) {
  // do not allow initialization from a vector smaller than
  // out own size:
  if ( d < Dim )
    throw aol::DimensionMismatchException ( aol::strprintf ( "GridSize<%i> "
              "constructor called with a Vec<%i>.", Dim, d ),
              __FILE__, __LINE__ );

  // initialize ourselves and assert that additional
  // components of Size are all one:
  for (int i = 0; i < Dim; ++i)
    (*this)[i] = static_cast<DataType> ( Size[i] );
  for (int i = Dim; i < d; ++i)
    if ( Size[i] != 1 )
      throw aol::DimensionMismatchException ( aol::strprintf("GridSize"
                  "<QC_%iD> constructor for Vec<%i> Size with Size[%i] = %i",
                         Dim, d, i, Size[i] ), __FILE__, __LINE__ );
    else
      CoordType::operator [] (i) = 1;
}

// --------------------------------------------------------------------------
template <Dimension _Dim>
typename aol::VecDimTrait<int,_Dim>::VecType GridSize<_Dim>::getSizeAsVecDim ( ) const {
  typename aol::VecDimTrait<int,_Dim>::VecType sizeVec;
  for ( int i = 0; i < _Dim; ++i )
    sizeVec[i] = (*this)[i];
  return sizeVec;
}
// --------------------------------------------------------------------------
template <Dimension _Dim>
void GridSize<_Dim>::quadraticOrDie () const {
  bool isQuadratic = true;
  for ( int i = 1; i < _Dim; ++i )
    isQuadratic &= ( (*this)[i] == (*this)[0] );
  if ( !isQuadratic )
    throw aol::Exception ( "GridSize should be quadratic here, but is not.",
                           __FILE__, __LINE__ );
}
// --------------------------------------------------------------------------
template <Dimension _Dim>
void GridSize<_Dim>::setAll ( typename GridSize<_Dim>::IndexType NumXYZ ) {
  for (int i = 0; i < Dim; ++i)
    (*this)[i] = NumXYZ;
  for (int i = Dim; i < dim; ++i)
    CoordType::operator [] (i) = 1;
}
// --------------------------------------------------------------------------

} // end of namespace qc.

#endif
