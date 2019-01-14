#ifndef __CONFIGURATORS_H
#define __CONFIGURATORS_H

#include <gridBase.h>
#include <indexMapper.h>
#include <baseFunctionSet.h>
#include <aol.h>
#include <rectangularGrid.h>
#include <FEOpInterface.h>
#include <matrix.h>
#include <pointerClasses.h>
#include <quocMatrices.h>
#include <diagBandMatrix.h>

namespace qc {

// Forward declaration because we can't include deformations.h here.
template <typename ConfiguratorType, bool ClipCoord>
inline bool transformCoord ( const typename ConfiguratorType::InitType &Grid,
                             const typename ConfiguratorType::ElementType &El,
                             const typename ConfiguratorType::VecType &RefCoord,
                             const typename ConfiguratorType::VecType &Offset,
                             qc::Element &TransformedEl,
                             typename ConfiguratorType::VecType &TransformedLocalCoord );

/**
 * Helper function to reduce code duplication in the configurators.
 *
 * \todo Use this function to get rid of the outrageous code duplication in QuocConfiguratorTraitMultiLin.
 *
 * \author Berkels, Droske
 */
template < typename ConfiguratorType >
bool getLocalCoordsRegularRectangularGrid ( const typename ConfiguratorType::VecType &Coord, const typename ConfiguratorType::InitType &Grid, qc::Element &El, typename ConfiguratorType::DomVecType &LocalCoord ) {
  for ( int c = 0; c < ConfiguratorType::Dim; ++c ) {
    const typename ConfiguratorType::RealType sc = Coord[c] / Grid.H();
    El[c] = static_cast<short> ( sc );
    LocalCoord[c] = sc - El[c];

    // We are exactly on the "right" boundary, adjust for this here.
    if ( ( El[c] == ( (Grid.getSize())[c] - 1 ) ) && ( aol::appeqAbsolute ( LocalCoord[c], aol::ZOTrait<typename ConfiguratorType::RealType>::zero ) ) ) {
      El[c] = (Grid.getSize())[c] - 2;
      LocalCoord[c] = 1;
    }

    // The position is not in the domain.
    if ( ( Coord[c] < 0. ) || ( El[c] < 0 ) || ( El[c] >= ( (Grid.getSize())[c] - 1 ) ) )
      return false;
  }
  return true;
}

/**
 * Helper class to provide dummy functions necessary to compile code involving DT_GRID_INDEX_MODE.
 *
 * \author Berkels
 */
class DTGridStubs {
public:
  inline void localToGlobal ( const qc::Element &/*El*/, const int /*localIndex0*/, const int /*localIndex1*/, aol::Vec2<int> &/*glob*/ ) const {
    throw aol::UnimplementedCodeException( "Special function only necessary for DT_GRID_INDEX_MODE", __FILE__, __LINE__);
  }
};

//! This serves as a basic default configuration for qc::Ops.
/*!
 *
 * It supplies begin and end iterator for the standard quocmesh grids, and
 * cares for default lexicographic mapping of the degrees of freedom.
 * \author Droske
 */
template <typename RealType, qc::Dimension Dim>
class QuocConfiguratorTraitBase {
public:

  typedef qc::Element                           ElementType;            //!< use the quoc element here
  typedef qc::GridDefinition::OldAllElementIterator  ElementIteratorType;    //!< use the standard element iterator of qc::GridDefinition
  // typedef UniformGridSparseMatrix<RealType> MatrixType;          //!< for uniform grids, the UniformSparseMatrices are most efficient
  // typedef qc::FastUniformGridMatrix<RealType,Dim>   MatrixType;
  typedef RealType Real;
  typedef qc::GridDefinition                    InitType;               //!< that's the type, that is needed by the constructor of the configurator
  typedef typename qc::BitArray<Dim>            MaskType;

  static const qc::Dimension DimOfWorld = Dim;

  explicit QuocConfiguratorTraitBase ( const InitType &Grid ) :
      _grid ( Grid ) {}

  //! returns the begin iterator of the grid
  const ElementIteratorType &begin( ) const {
    return _grid.begin_it;
  }

  //! returns the end iterator of the grid
  inline const ElementIteratorType &end( ) const {
    return _grid.end_it;
  }

  RealType H ( const qc::Element& ) const {
    return _grid.H();
  }

  const InitType& getInitializer( ) const { return _grid; }

  //! method needed in boundaryIntegration.h
  //! although the name suggests something different this implementation reflects the current behavior
  inline int localOnFaceToLocal ( const qc::Element &/*El*/, const int localIndex ) const {
    return localIndex;
  }

protected:
  const InitType &_grid;   //!< memorize reference to the grid
  //qc::IndexMapper _mapper;         //!< that's the index mapper
};

template < typename RealType,
           qc::Dimension Dim,
           typename QuadType,
           typename _MatrixType = qc::FastUniformGridMatrix< RealType, Dim> >
class QuocConfiguratorTraitMultiLin {};

/**
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename _QuadType, typename _MatrixType >
class QuocConfiguratorTraitMultiLin<_RealType, qc::QC_1D, _QuadType, _MatrixType> : public QuocConfiguratorTraitBase<_RealType, qc::QC_1D> {
  aol::BaseFunctionSetMultiLin<_RealType, qc::QC_1D, _QuadType> _baseFuncSet;
public:
  explicit QuocConfiguratorTraitMultiLin ( const qc::GridDefinition &Grid )
      : QuocConfiguratorTraitBase<_RealType, qc::QC_1D> ( Grid ),  _baseFuncSet ( Grid.H() ), _volEl ( Grid.H() ), _mapper ( Grid ) {}

  typedef aol::Vector<_RealType>       VectorType;
  typedef aol::Vec<1, _RealType>       VecType;
  typedef aol::Vec<1, _RealType>       DomVecType;
  typedef aol::Mat<1, 1, _RealType>    MatType;
  typedef _MatrixType                  MatrixType;
  typedef aol::BaseFunctionSetMultiLin<_RealType, qc::QC_1D, _QuadType> BaseFuncSetType;
  typedef _QuadType                    QuadType;
  typedef _RealType                    RealType;
  typedef qc::ScalarArray<_RealType, qc::QC_1D> ArrayType;
  typedef aol::FullMatrix<RealType>    FullMatrixType;
  typedef qc::QuocConfiguratorTraitMultiLin<_RealType, qc::QC_1D, _QuadType, _MatrixType> FullGridConfiguratorType;

  static const int maxNumLocalDofs = 2;

  static const qc::Dimension Dim = qc::QC_1D;
  static const qc::Dimension DomDim = qc::QC_1D;

  const RealType _volEl;

  inline bool getLocalCoords ( const VecType &Coord, qc::Element &El, DomVecType &LocalCoord ) const {
    return getLocalCoordsRegularRectangularGrid<QuocConfiguratorTraitMultiLin<_RealType, qc::QC_1D, _QuadType, _MatrixType> > ( Coord, this->_grid, El, LocalCoord );
  }

  int maxNumQuadPoints( ) const {
    return _QuadType::numQuadPoints;
  }

  inline void getGlobalCoords ( const qc::Element &El, const DomVecType &LocalCoord, VecType &Coord ) const {
    for ( int c = 0; c < 1; c++ ) {
      Coord[c] = ( static_cast<RealType> ( El[c] ) + LocalCoord[c] ) * this->_grid.H();
    }
  }

  inline int getNumLocalDofs ( const qc::Element & ) const {
    return 2;
  }

  int getNumGlobalDofs( ) const {
    return this->_grid.getNumberOfNodes();
  }

  qc::Element getEmptyElement() const {
    return qc::Element();
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element& /*El*/ ) const {
    return _baseFuncSet;
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
    return _mapper.localToGlobal ( El, localIndex );
  }

  inline void localToGlobal ( const qc::Element &El, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal( El, localIndex0 );
    glob[1] = localToGlobal( El, localIndex1 );
  }

  RealType vol ( const qc::Element & ) const {
    return _volEl;
  }

  RealType H ( const qc::Element & ) const {
    return this->_grid.H();
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( qc::GridSize<qc::QC_1D>::createFrom ( this->_grid ) );
    // mat->clearRows( );
    return mat;
  }

  const qc::FastILexMapper<Dim>& getIndexMapper() const {
    return _mapper;
  }

  using QuocConfiguratorTraitBase<_RealType, qc::QC_1D>::localOnFaceToLocal;

  template <bool ClipCoord>
  inline bool transformCoord ( const qc::Element &El, const VecType &RefCoord, const VecType &Offset, qc::Element &TransformedEl, VecType &TransformedLocalCoord ) const {
    return qc::transformCoord<QuocConfiguratorTraitMultiLin<_RealType, qc::QC_1D, _QuadType, _MatrixType>, ClipCoord> ( this->_grid, El, RefCoord, Offset, TransformedEl, TransformedLocalCoord );
  }
protected:
  qc::FastILexMapper<Dim> _mapper;         //!< that's the index mapper
};

/**
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename _QuadType, typename _MatrixType >
class QuocConfiguratorTraitMultiLin<_RealType, qc::QC_2D, _QuadType, _MatrixType> : public QuocConfiguratorTraitBase<_RealType, qc::QC_2D> {
  aol::BaseFunctionSetMultiLin<_RealType, qc::QC_2D, _QuadType> _baseFuncSet;
public:
  explicit QuocConfiguratorTraitMultiLin ( const qc::GridDefinition &Grid )
      : QuocConfiguratorTraitBase<_RealType, qc::QC_2D> ( Grid ),  _baseFuncSet ( Grid.H() ), _volEl ( aol::Sqr ( Grid.H() ) ), _mapper ( Grid ) {}

  typedef aol::Vector<_RealType>       VectorType;
  typedef aol::Vec2<_RealType>         VecType;
  typedef aol::Vec2<_RealType>         DomVecType;
  typedef aol::Matrix22<_RealType>     MatType;
  typedef _MatrixType                  MatrixType;
  typedef aol::BaseFunctionSetMultiLin<_RealType, qc::QC_2D, _QuadType> BaseFuncSetType;
  typedef _QuadType                    QuadType;
  typedef _RealType                    RealType;
  typedef qc::ScalarArray<_RealType, qc::QC_2D> ArrayType;
  typedef qc::BitArray<qc::QC_2D>               MaskType;
  typedef aol::FullMatrix<RealType>    FullMatrixType;
  typedef qc::QuocConfiguratorTraitMultiLin<_RealType, qc::QC_2D, _QuadType, _MatrixType> FullGridConfiguratorType;

  static const int maxNumLocalDofs = 4;

  static const qc::Dimension Dim = qc::QC_2D;
  static const qc::Dimension DomDim = qc::QC_2D;
  static const aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE;

  const RealType _volEl;

  inline bool getLocalCoords ( const VecType &Coord, qc::Element &El, DomVecType &LocalCoord ) const {
    return getLocalCoordsRegularRectangularGrid<QuocConfiguratorTraitMultiLin<_RealType, qc::QC_2D, _QuadType, _MatrixType> > ( Coord, this->_grid, El, LocalCoord );
  }


  int maxNumQuadPoints( ) const {
    return _QuadType::numQuadPoints;
  }

  inline void getGlobalCoords ( const qc::Element &El, const DomVecType &LocalCoord, VecType &Coord ) const {
    for ( int c = 0; c < 2; c++ ) {
      Coord[c] = ( static_cast<RealType> ( El[c] ) + LocalCoord[c] ) * this->_grid.H();
    }
  }

  inline int getNumLocalDofs ( const qc::Element & ) const {
    return 4;
  }

  int getNumGlobalDofs( ) const {
    return this->_grid.getNumberOfNodes();
  }

  qc::Element getEmptyElement() const {
    return qc::Element();
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element& /*El*/ ) const {
    return _baseFuncSet;
  }

  //! get Number of Element ( = Number of ll dof ) (BG)
  int getElementNumber ( const qc::Element & El ) const {
    return _mapper.localToGlobal ( El, 0 );
  }

  //! compute consecutive element numbers (BG)
  int getConsecutiveElementNumber ( const qc::Element & El ) const {
    int elNum = getElementNumber( El );
    return elNum - ( elNum / this->_grid.getNumX() );
  }
  
  //! get element from consecutive number
  void getElementfromConsecutiveNumber ( int number, qc::Element & El ) const {
    number += number / (this->_grid.getNumX()-1);
    _mapper.splitGlobalIndex ( number, El ) ;
    El.setLevel( this->_grid.getGridDepth() );
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
    return _mapper.localToGlobal ( El, localIndex );
  }

  //! returns a vec2 of global indices i,j of the dofs with number localIndex_i,j
  inline void localToGlobal ( const qc::Element &El, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal( El, localIndex0 );
    glob[1] = localToGlobal( El, localIndex1 );
  }


  RealType vol ( const qc::Element & ) const {
    return _volEl;
  }

  RealType H ( const qc::Element & ) const {
    return this->_grid.H();
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( qc::GridSize<qc::QC_2D>::createFrom ( this->_grid ) );
    // mat->clearRows( );
    return mat;
  }

  const qc::FastILexMapper<Dim>& getIndexMapper() const {
    return _mapper;
  }

  using QuocConfiguratorTraitBase<_RealType, qc::QC_2D>::localOnFaceToLocal;

  template <bool ClipCoord>
  inline bool transformCoord ( const qc::Element &El, const VecType &RefCoord, const VecType &Offset, qc::Element &TransformedEl, VecType &TransformedLocalCoord ) const {
    return qc::transformCoord<QuocConfiguratorTraitMultiLin<_RealType, qc::QC_2D, _QuadType, _MatrixType>, ClipCoord> ( this->_grid, El, RefCoord, Offset, TransformedEl, TransformedLocalCoord );
  }
protected:
  qc::FastILexMapper<Dim> _mapper;         //!< that's the index mapper
};


/**
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename _QuadType, typename _MatrixType >
class QuocConfiguratorTraitMultiLin<_RealType, qc::QC_3D, _QuadType, _MatrixType> : public QuocConfiguratorTraitBase<_RealType, qc::QC_3D> {
  aol::BaseFunctionSetMultiLin<_RealType, qc::QC_3D, _QuadType> _baseFuncSet;
public:
  explicit QuocConfiguratorTraitMultiLin ( const qc::GridDefinition &Grid )
    : QuocConfiguratorTraitBase<_RealType, qc::QC_3D> ( Grid ), _baseFuncSet( Grid.H() ), _volEl ( aol::Cub ( Grid.H() ) ), _mapper ( Grid )  {}

  typedef qc::QuocConfiguratorTraitMultiLin<_RealType, qc::QC_3D, _QuadType, _MatrixType> Self;

  typedef aol::Vector<_RealType>       VectorType;
  typedef aol::Vec3<_RealType>         VecType;
  typedef aol::Vec3<_RealType>         DomVecType;
  typedef aol::Matrix33<_RealType>     MatType;
  typedef _MatrixType                  MatrixType;
  typedef aol::BaseFunctionSetMultiLin<_RealType, qc::QC_3D, _QuadType> BaseFuncSetType;
  typedef _QuadType                    QuadType;
  typedef _RealType                    RealType;
  typedef qc::ScalarArray<_RealType, qc::QC_3D> ArrayType;
  typedef qc::BitArray<qc::QC_3D>               MaskType;
  typedef aol::FullMatrix<RealType>    FullMatrixType;
  typedef Self                         FullGridConfiguratorType;

  static const int maxNumLocalDofs = 8;

  static const qc::Dimension Dim = qc::QC_3D;
  static const qc::Dimension DomDim = qc::QC_3D;
  static const aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE;

  const RealType _volEl;

  inline int getNumLocalDofs ( const qc::Element & ) const {
    return 8;
  }

  int getNumGlobalDofs( ) const {
    return this->_grid.getNumberOfNodes();
  }

  int maxNumQuadPoints( ) const {
    return _QuadType::numQuadPoints;
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element& ) const {
    return _baseFuncSet;
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
    return _mapper.localToGlobal ( El, localIndex );
  }

  //! returns a vec2 of global indices i,j of the dofs with number localIndex_i,j
  inline void localToGlobal ( const qc::Element &El, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal( El, localIndex0 );
    glob[1] = localToGlobal( El, localIndex1 );
  }

  inline bool getLocalCoords ( const VecType &Coord, qc::Element &El, DomVecType &LocalCoord ) const {
    return getLocalCoordsRegularRectangularGrid<QuocConfiguratorTraitMultiLin<_RealType, qc::QC_3D, _QuadType, _MatrixType> > ( Coord, this->_grid, El, LocalCoord );
  }

  inline void getGlobalCoords ( const qc::Element &El, const DomVecType &LocalCoord, VecType &Coord ) const {
    for ( int c = 0; c < 3; c++ ) {
      Coord[c] = ( static_cast<RealType> ( El[c] ) + LocalCoord[c] ) * this->_grid.H();
    }
  }

  qc::Element getEmptyElement() const {
    return qc::Element();
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( qc::GridSize<qc::QC_3D>::createFrom ( this->_grid ) );
    // mat->clearRows( );
    return mat;
  }

  RealType vol ( const qc::Element & ) const {
    return _volEl;
  }

  RealType H ( const qc::Element & ) const {
    return this->_grid.H();
  }

  //! returns consecutive Element Number (ST)
  int getConsecutiveElementNumber (const qc::Element & El) const {
    int ElNum = localToGlobal(El, 0);
    return ElNum - ElNum / this->_grid.getNumX() - ElNum / ( this->_grid.getNumX() * this->_grid.getNumY() ) * ( this->_grid.getNumX() - 1 );
  }

  using QuocConfiguratorTraitBase<_RealType, qc::QC_3D>::localOnFaceToLocal;

  template <bool ClipCoord>
  inline bool transformCoord ( const qc::Element &El, const VecType &RefCoord, const VecType &Offset, qc::Element &TransformedEl, VecType &TransformedLocalCoord ) const {
    return qc::transformCoord<QuocConfiguratorTraitMultiLin<_RealType, qc::QC_3D, _QuadType, _MatrixType>, ClipCoord> ( this->_grid, El, RefCoord, Offset, TransformedEl, TransformedLocalCoord );
  }
protected:
  qc::FastILexMapper<qc::QC_3D> _mapper;         //!< that's the index mapper

};


template <typename RealType, qc::Dimension Dim, typename QuadType, typename _MatrixType = aol::SparseMatrix<RealType> >
class QuocConfiguratorTraitMultiQuad {
};

/**
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename _QuadType, typename _MatrixType >
  class QuocConfiguratorTraitMultiQuad<_RealType, qc::QC_2D, _QuadType, _MatrixType> : public QuocConfiguratorTraitBase<_RealType, qc::QC_2D> {
  aol::BaseFunctionSetMultiQuad<_RealType, qc::QC_2D, _QuadType> _baseFuncSet;
public:
  explicit QuocConfiguratorTraitMultiQuad ( const qc::GridDefinition &Grid )
    : QuocConfiguratorTraitBase<_RealType, qc::QC_2D> ( Grid ), _baseFuncSet( Grid.H() ), _volEl( aol::Sqr( Grid.H() )), _mapper ( Grid ) {}

  typedef qc::QuocConfiguratorTraitMultiQuad<_RealType, qc::QC_2D, _QuadType, _MatrixType> Self;

  typedef aol::Vector<_RealType>       VectorType;
  typedef aol::Vec2<_RealType>         VecType;
  typedef aol::Vec2<_RealType>         DomVecType;
  typedef aol::Matrix22<_RealType>     MatType;
  typedef _MatrixType                  MatrixType;
  typedef aol::BaseFunctionSetMultiQuad<_RealType, qc::QC_2D, _QuadType> BaseFuncSetType;
  typedef _QuadType                    QuadType;
  typedef _RealType                    RealType;
  typedef qc::BitArray<qc::QC_2D>      MaskType;
  typedef ScalarArray<RealType, QC_2D> ArrayType;
  typedef aol::FullMatrix<RealType>    FullMatrixType;

  static const int maxNumLocalDofs = 9;
  static const aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE;

  static const qc::Dimension Dim = qc::QC_2D;
  static const qc::Dimension DomDim = qc::QC_2D;

  const RealType _volEl;

  inline int getNumLocalDofs ( const qc::Element & ) const {
    return 9;
  }

  int getNumGlobalDofs( ) const {
    return aol::Sqr ( ( this->_grid.getWidth() - 1 ) *2 + 1 );
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element &/*El*/ ) const {
    return _baseFuncSet;
  }

  int maxNumQuadPoints( ) const {
    return _QuadType::numQuadPoints;
  }

  RealType vol ( const qc::Element & ) const {
    return _volEl;
  }

  RealType H ( const qc::Element & ) const {
    return this->_grid.H();
  }

  //! get Number of Element ( = Number of ll dof ) (BG)
  int getElementNumber ( const qc::Element & El ) const {
    return _mapper.localToGlobal ( El, 0 );
  }

  //! compute consecutive element numbers (BG)
  int getConsecutiveElementNumber ( const qc::Element & El ) const {
    int elNum = getElementNumber( El );
    return ( elNum - ( elNum / (2*this->_grid.getWidth()-1) ) )                     // eliminate last node index in each row
           / 2                                                                      // forget about virtual nodes in between
           - ( elNum / (4*this->_grid.getWidth()-2) ) * (this->_grid.getWidth()-1); // subtract intermediate rows of virtual nodes
  }

  //! get element from consecutive number
  void getElementfromConsecutiveNumber ( int number, qc::Element & El ) const {
    // element numbers only refer to real dofs
    number += number / (this->_grid.getNumX()-1);
    // do manual splitting wrt reals dofs
    El.set( number % this->_grid.getNumX(), number / this->_grid.getNumX(), 0, this->_grid.getGridDepth() );
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
    return _mapper.localToGlobal ( El, localIndex );
  }

  //! returns a vec2 of global indices i,j of the dofs with number localIndex_i,j
  inline void localToGlobal ( const qc::Element &El, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal( El, localIndex0 );
    glob[1] = localToGlobal( El, localIndex1 );
  }

  inline void getGlobalCoords ( const qc::Element &El, const DomVecType &LocalCoord, VecType &Coord ) const {
    for ( int c = 0; c < 2; c++ ) {
      Coord[c] = ( static_cast<RealType> ( El[c] ) + LocalCoord[c] ) * this->_grid.H();
    }
  }

  inline bool getLocalCoords ( const VecType &Coord, qc::Element &El, DomVecType &LocalCoord ) const {
    return getLocalCoordsRegularRectangularGrid<QuocConfiguratorTraitMultiQuad<_RealType, qc::QC_2D, _QuadType, _MatrixType> > ( Coord, this->_grid, El, LocalCoord );
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( getNumGlobalDofs( ), getNumGlobalDofs( ) );
    mat->setZero( );
    return mat;
  }

  using QuocConfiguratorTraitBase<_RealType, qc::QC_2D>::localOnFaceToLocal;

protected:
  qc::FastQuadILexMapper<qc::QC_2D> _mapper;         //!< that's the index mapper
};


template <typename RealType, qc::Dimension Dim, typename QuadType, typename _MatrixType = aol::SparseMatrix<RealType> >
class QuocConfiguratorTraitMultiQuart {
};

/**
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename _QuadType, typename _MatrixType >
class QuocConfiguratorTraitMultiQuart<_RealType, qc::QC_2D, _QuadType, _MatrixType> : public QuocConfiguratorTraitBase<_RealType, qc::QC_2D> {
  aol::BaseFunctionSetMultiQuart<_RealType, qc::QC_2D, _QuadType> _baseFuncSet;
public:
  explicit QuocConfiguratorTraitMultiQuart ( const qc::GridDefinition &Grid )
  : QuocConfiguratorTraitBase<_RealType, qc::QC_2D> ( Grid ), _baseFuncSet( Grid.H() ), _volEl( aol::Sqr( Grid.H() )), _mapper ( Grid ) {}

  typedef qc::QuocConfiguratorTraitMultiQuart<_RealType, qc::QC_2D, _QuadType, _MatrixType> Self;

  typedef aol::Vector<_RealType>       VectorType;
  typedef aol::Vec2<_RealType>         VecType;
  typedef aol::Vec2<_RealType>         DomVecType;
  typedef aol::Matrix22<_RealType>     MatType;
  typedef _MatrixType                  MatrixType;
  typedef aol::BaseFunctionSetMultiQuart<_RealType, qc::QC_2D, _QuadType> BaseFuncSetType;
  typedef _QuadType                    QuadType;
  typedef _RealType                    RealType;
  typedef qc::BitArray<qc::QC_2D>      MaskType;
  typedef ScalarArray<RealType, QC_2D> ArrayType;
  typedef aol::FullMatrix<RealType>    FullMatrixType;

  static const int maxNumLocalDofs = 25;
  static const aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE;

  static const qc::Dimension Dim = qc::QC_2D;
  static const qc::Dimension DomDim = qc::QC_2D;

  const RealType _volEl;

  inline int getNumLocalDofs ( const qc::Element & ) const {
    return 25;
  }

  int getNumGlobalDofs( ) const {
    return aol::Sqr ( ( this->_grid.getWidth() - 1 ) * 4 + 1 );
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element &/*El*/ ) const {
    return _baseFuncSet;
  }

  int maxNumQuadPoints( ) const {
    return _QuadType::numQuadPoints;
  }

  RealType vol ( const qc::Element & ) const {
    return _volEl;
  }

  RealType H ( const qc::Element & ) const {
    return this->_grid.H();
  }

  //! get Number of Element ( = Number of ll dof ) (BG)
  int getElementNumber ( const qc::Element & El ) const {
    return _mapper.localToGlobal ( El, 0 );
  }

  //! compute consecutive element numbers (BG)
  int getConsecutiveElementNumber ( const qc::Element & El ) const {
    int elNum = getElementNumber( El );
    return ( elNum - ( elNum / (4*this->_grid.getWidth()-3) ) )                           // eliminate last node index in each row
           / 4                                                                            // forget about virtual nodes in between
           - ( elNum / (16*this->_grid.getWidth()-12) * 3 ) * (this->_grid.getWidth()-1); // subtract intermediate rows of virtual nodes
  }

  //! get element from consecutive number
  void getElementfromConsecutiveNumber ( int number, qc::Element & El ) const {
    // element numbers only refer to real dofs
    number += number / (this->_grid.getNumX()-1);
    // do manual splitting wrt reals dofs
    El.set( number % this->_grid.getNumX(), number / this->_grid.getNumX(), 0, this->_grid.getGridDepth() );
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
    return _mapper.localToGlobal ( El, localIndex );
  }

  //! returns a vec2 of global indices i,j of the dofs with number localIndex_i,j
  inline void localToGlobal ( const qc::Element &El, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal( El, localIndex0 );
    glob[1] = localToGlobal( El, localIndex1 );
  }

  inline void getGlobalCoords ( const qc::Element &El, const DomVecType &LocalCoord, VecType &Coord ) const {
    for ( int c = 0; c < 2; c++ ) {
      Coord[c] = ( static_cast<RealType> ( El[c] ) + LocalCoord[c] ) * this->_grid.H();
    }
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( getNumGlobalDofs( ), getNumGlobalDofs( ) );
    mat->setZero( );
    return mat;
  }

protected:
  qc::FastQuartILexMapper<qc::QC_2D> _mapper;         //!< that's the index mapper
};


template <typename _RealType, qc::Dimension Dim, typename RectangularGridType = RectangularGrid<Dim> >
class RectangularGridConfiguratorBase : public DTGridStubs {
public:
  typedef typename RectangularGrid<Dim>::OldAllElementIterator ElementIteratorType;    //!< use the standard element iterator of qc::RectangularGrid
  typedef qc::Element                             ElementType;            //!< use the quoc element here
  typedef RectangularGridType                     InitType;               //!< that's the type, that is needed by the constructor of the configurator
  typedef _RealType                               RealType;
  typedef qc::MultilinFEBandMatrix<_RealType,Dim> MatrixType;

  static const qc::Dimension DimOfWorld = Dim;

  static const aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE;

protected:
  const InitType &_grid;   //!< memorize reference to the grid

public:
  const aol::Vec3<int> &_size;   //!< memorize grid size

  RectangularGridConfiguratorBase ( const InitType &Grid ) :
      _grid ( Grid ), _size ( Grid.getSize() ) {
#ifdef VERBOSE
    cerr << "constructor = " << this->_size;
#endif
  }
  //! returns the begin iterator of the grid
  const typename RectangularGrid<Dim>::OldAllElementIterator &begin( ) const {
    return this->_grid.begin_it;
  }

  //! returns the end iterator of the grid
  inline const typename RectangularGrid<Dim>::OldAllElementIterator &end( ) const {
    return this->_grid.end_it;
  }

  RealType H ( const qc::Element& ) const {
    return this->_grid.H();
  }


  const InitType& getInitializer( ) const { return this->_grid; }

  int getNumGlobalDofs( ) const {
    return this->_grid.getNumberOfNodes();
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( GridSize<Dim> ( _grid ) );
    return mat;
  }

  int getConsecutiveElementNumber ( const qc::Element &/*El*/ ) const {
    throw aol::UnimplementedCodeException( "RectangularGridConfiguratorBase::getConsecutiveElementNumber not implemented", __FILE__, __LINE__);
    return -1;
  }

  //! method needed in boundaryIntegration.h
  //! copy and paste from qc::QuocConfiguratorTraitBase
  inline int localOnFaceToLocal ( const qc::Element &/*El*/, const int localIndex ) const {
    return localIndex;
  }
};


template <typename _RealType, qc::Dimension Dim, typename _QuadType, typename RectangularGridType = RectangularGrid<Dim> >
class RectangularGridConfigurator { };

// special stuff for 1d
/**
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename _QuadType, typename RectangularGridType>
class RectangularGridConfigurator<_RealType, qc::QC_1D, _QuadType, RectangularGridType>
      : public RectangularGridConfiguratorBase<_RealType, qc::QC_1D, RectangularGridType> {
public:
  RectangularGridConfigurator ( const RectangularGridType &Grid )
      : RectangularGridConfiguratorBase<_RealType, qc::QC_1D, RectangularGridType> ( Grid ), _baseFuncSet ( static_cast<_RealType> ( Grid.H() ) ), _volEl ( static_cast<_RealType> ( Grid.H() ) ) {}

  typedef aol::Vector<_RealType>       VectorType;
  typedef _RealType                    RealType;
  typedef aol::Vec<1, _RealType>       VecType;
  typedef aol::Vec<1, _RealType>       DomVecType;
  typedef aol::Mat<1, 1, _RealType>    MatType;
  typedef aol::DiagBandMatrix<_RealType, 1, 1> MatrixType;
  typedef aol::BaseFunctionSetMultiLin<_RealType, qc::QC_1D, _QuadType> BaseFuncSetType;
  typedef _QuadType                    QuadType;
  typedef qc::ScalarArray<_RealType, qc::QC_1D> ArrayType;
  typedef aol::FullMatrix<_RealType>   FullMatrixType;


  BaseFuncSetType _baseFuncSet;

  const RealType _volEl;

  static const int maxNumLocalDofs = 2;

  static const qc::Dimension Dim = qc::QC_1D;
  static const qc::Dimension DomDim = qc::QC_1D;

  inline bool getLocalCoords ( const VecType &Coord, qc::Element &El, DomVecType &LocalCoord ) const {
    return getLocalCoordsRegularRectangularGrid <RectangularGridConfigurator<_RealType, qc::QC_1D, _QuadType, RectangularGridType> > ( Coord, this->_grid, El, LocalCoord );
  }

  inline int getNumLocalDofs ( const qc::Element & ) const {
    return 2;
  }

  inline void getGlobalCoords ( const qc::Element &El, const DomVecType &LocalCoord, VecType &Coord ) const {
    Coord[0] = ( static_cast<RealType> ( El[0] ) + LocalCoord[0] ) * this->_grid.H();
  }

  inline RealType vol ( const qc::Element & ) const {
    return _volEl;
  }

  int maxNumQuadPoints( ) const {
    return _QuadType::numQuadPoints;
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element &/*El*/ ) const {
    return _baseFuncSet;
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
#ifdef VERBOSE
    cerr << this->_size;
#endif
    int r = ( El.x() + ( ( localIndex & 1 ) != 0 ) );
#ifdef VERBOSE
    cerr << " r = " << r << " localIndex = " << localIndex << endl;
#endif
    return r;
  }

  using DTGridStubs::localToGlobal;

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( this->_grid.getNumberOfNodes() );
    return mat;
  }

  template <bool ClipCoord>
  inline bool transformCoord ( const qc::Element &El, const VecType &RefCoord, const VecType &Offset, qc::Element &TransformedEl, VecType &TransformedLocalCoord ) const {
    return qc::transformCoord<RectangularGridConfigurator<_RealType, qc::QC_1D, _QuadType, RectangularGridType>, ClipCoord> ( this->_grid, El, RefCoord, Offset, TransformedEl, TransformedLocalCoord );
  }
};

// special stuff for 2d
/**
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename _QuadType, typename RectangularGridType>
class RectangularGridConfigurator<_RealType, qc::QC_2D, _QuadType, RectangularGridType>
      : public RectangularGridConfiguratorBase<_RealType, qc::QC_2D, RectangularGridType> {
public:
  RectangularGridConfigurator ( const RectangularGridType &Grid )
      : RectangularGridConfiguratorBase<_RealType, qc::QC_2D, RectangularGridType> ( Grid ), _baseFuncSet ( static_cast<_RealType> ( Grid.H() ) ), _volEl ( static_cast<_RealType> ( aol::Sqr ( Grid.H() ) ) ) {}

  typedef aol::Vector<_RealType>       VectorType;
  typedef _RealType                    RealType;
  typedef aol::Vec2<_RealType>         VecType;
  typedef aol::Vec2<_RealType>         DomVecType;
  typedef aol::Matrix22<_RealType>     MatType;
  typedef qc::MultilinFEBandMatrix<_RealType,qc::QC_2D>  MatrixType;
  typedef aol::BaseFunctionSetMultiLin<_RealType, qc::QC_2D, _QuadType> BaseFuncSetType;
  typedef _QuadType                    QuadType;
  typedef qc::ScalarArray<_RealType, qc::QC_2D> ArrayType;
  typedef qc::BitArray<qc::QC_2D>               MaskType;
  typedef aol::FullMatrix<_RealType>   FullMatrixType;


  BaseFuncSetType _baseFuncSet;

  const RealType _volEl;

  int getConsecutiveElementNumber ( const qc::Element & El ) const {
    int elNum = localToGlobal ( El, 0 );
    return elNum - ( elNum / this->_grid.getNumX() );
  }
  static const int maxNumLocalDofs = 4;

  static const qc::Dimension Dim = qc::QC_2D;
  static const qc::Dimension DomDim = qc::QC_2D;

  inline bool getLocalCoords ( const VecType &Coord, qc::Element &El, DomVecType &LocalCoord ) const {
    return getLocalCoordsRegularRectangularGrid <RectangularGridConfigurator<_RealType, qc::QC_2D, _QuadType, RectangularGridType> > ( Coord, this->_grid, El, LocalCoord );
  }

  inline void getGlobalCoords ( const qc::Element &El, const DomVecType &LocalCoord, VecType &Coord ) const {
    for ( int c = 0; c < 2; c++ ) {
      Coord[c] = ( static_cast<RealType> ( El[c] ) + LocalCoord[c] ) * this->_grid.H();
    }
  }

  inline int getNumLocalDofs ( const qc::Element & ) const {
    return 4;
  }

  inline RealType vol ( const qc::Element & ) const {
    return _volEl;
  }

  int maxNumQuadPoints( ) const {
    return _QuadType::numQuadPoints;
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element &/*El*/ ) const {
    return _baseFuncSet;
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
#ifdef VERBOSE
    cerr << this->_size;
#endif
    int r =
      ( El.x() + ( ( localIndex & 1 ) != 0 ) ) * 1       +
      ( El.y() + ( ( localIndex & 2 ) != 0 ) ) * this->_size[0];
#ifdef VERBOSE
    cerr << " r = " << r << " localIndex = " << localIndex << endl;
#endif
    return r;
  }

  using DTGridStubs::localToGlobal;

  template <bool ClipCoord>
  inline bool transformCoord ( const qc::Element &El, const VecType &RefCoord, const VecType &Offset, qc::Element &TransformedEl, VecType &TransformedLocalCoord ) const {
    return qc::transformCoord<RectangularGridConfigurator<_RealType, qc::QC_2D, _QuadType, RectangularGridType>, ClipCoord> ( this->_grid, El, RefCoord, Offset, TransformedEl, TransformedLocalCoord );
  }
};

// special stuff for 3d
/**
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename _QuadType, typename RectangularGridType>
class RectangularGridConfigurator<_RealType, qc::QC_3D, _QuadType, RectangularGridType>
      : public RectangularGridConfiguratorBase<_RealType, qc::QC_3D, RectangularGridType> {
public:
  RectangularGridConfigurator ( const RectangularGridType &Grid )
      : RectangularGridConfiguratorBase<_RealType, qc::QC_3D, RectangularGridType> ( Grid ), _baseFuncSet( Grid.H()), _volEl ( aol::Cub ( Grid.H() ) ) {}
public:
  typedef aol::Vector<_RealType>       VectorType;
  typedef qc::Element                  ElementType;            //!< use the quoc element here
  typedef _RealType                    RealType;
  typedef aol::Vec3<_RealType>         DomVecType;
  typedef aol::Vec3<_RealType>         VecType;
  typedef aol::Matrix33<_RealType>     MatType;
  typedef qc::MultilinFEBandMatrix<_RealType,qc::QC_3D>  MatrixType;
  typedef aol::BaseFunctionSetMultiLin<_RealType, qc::QC_3D, _QuadType> BaseFuncSetType;
  typedef qc::ScalarArray<_RealType, qc::QC_3D> ArrayType;
  typedef qc::BitArray<qc::QC_3D>               MaskType;
  typedef aol::FullMatrix<_RealType>   FullMatrixType;

  typedef _QuadType QuadType;

  BaseFuncSetType _baseFuncSet;

  const RealType _volEl;


 int getConsecutiveElementNumber (const qc::Element & El) const {
   int ElNum = localToGlobal(El, 0);
   return ElNum - ElNum / this->_grid.getNumX() - ElNum / ( this->_grid.getNumX() * this->_grid.getNumY() ) * ( this->_grid.getNumX() - 1 );
 }
  static const int maxNumLocalDofs = 8;

  static const qc::Dimension Dim = qc::QC_3D;
  static const qc::Dimension DomDim = qc::QC_3D;

  inline RealType vol ( const qc::Element & ) const {
    return _volEl;
  }

  inline bool getLocalCoords ( const VecType &Coord, qc::Element &El, DomVecType &LocalCoord ) const {
    return getLocalCoordsRegularRectangularGrid <RectangularGridConfigurator<_RealType, qc::QC_3D, _QuadType, RectangularGridType> > ( Coord, this->_grid, El, LocalCoord );
  }

  inline int getNumLocalDofs ( const qc::Element & ) const {
    return 8;
  }

  int maxNumQuadPoints( ) const {
    return _QuadType::numQuadPoints;
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element &/*El*/ ) const {
    return _baseFuncSet;
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
#ifdef VERBOSE
    cerr << this->_size;
#endif
    int r =
      ( El.x() + ( ( localIndex & 1 ) != 0 ) ) * 1       +
      ( El.y() + ( ( localIndex & 2 ) != 0 ) ) * this->_size[0] +
      ( El.z() + ( ( localIndex & 4 ) != 0 ) ) * this->_size[0] * this->_size[1];
#ifdef VERBOSE
    cerr << " r = " << r << " localIndex = " << localIndex << endl;
#endif
    return r;
  }

  using DTGridStubs::localToGlobal;

  template <bool ClipCoord>
  inline bool transformCoord ( const qc::Element &El, const VecType &RefCoord, const VecType &Offset, qc::Element &TransformedEl, VecType &TransformedLocalCoord ) const {
    return qc::transformCoord<RectangularGridConfigurator<_RealType, qc::QC_3D, _QuadType, RectangularGridType>, ClipCoord> ( this->_grid, El, RefCoord, Offset, TransformedEl, TransformedLocalCoord );
  }
};


// compatibility classes
template <typename _RealType, typename _QuadType>
  class RectangularGridConfigurator3D : public qc::RectangularGridConfigurator<_RealType, qc::QC_3D, _QuadType> { };

template <typename _RealType, typename _QuadType>
  class RectangularGridConfigurator2D : public qc::RectangularGridConfigurator<_RealType, qc::QC_2D, _QuadType> { };


template <typename _RealType, qc::Dimension Dim, typename _InitType>
class RegularSimplexConfiguratorBase {
public:
  typedef aol::Vector<_RealType>       VectorType;
  typedef typename _InitType::OldFullElementIterator
                                                  ElementIteratorType; //!< use the standard element iterator of _InitType
  typedef typename _InitType::BeginIterType BeginIterType;
  typedef typename _InitType::EndIterType EndIterType;
  typedef qc::Element                  ElementType;         //!< use the quoc element here
  typedef _InitType                    InitType;            //!< that's the type, that is needed by the constructor of the configurator
  typedef _RealType                    RealType;
  typedef aol::SparseMatrix<_RealType> MatrixType;

protected:
  const InitType &_grid;                                     //!< memorize reference to the grid

public:
  const aol::Vec3<int> &_size;   //!< memorize grid size
  RegularSimplexConfiguratorBase ( const InitType &Grid ) :
      _grid ( Grid ), _size ( Grid.getSize() ) {
#ifdef VERBOSE
    cerr << "constructor = " << this->_size;
#endif
  }
  //! returns the begin iterator of the grid
  const BeginIterType &begin( ) const {
    return this->_grid.begin();
  }

  //! returns the end iterator of the grid
  inline EndIterType end( ) const {
    return this->_grid.end();
  }

  RealType H ( const qc::Element& ) const {
    return this->_grid.H();
  }


  const InitType& getInitializer( ) const { return this->_grid; }

  int getNumGlobalDofs( ) const {
    return this->_grid.getNumberOfNodes();
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( getNumGlobalDofs(), getNumGlobalDofs() );
    return mat;
  }

};


template <typename _RealType, qc::Dimension, typename _QuadType, typename _InitType>
class RegularSimplexConfigurator { };


// special stuff for 2d
/**
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename _QuadType, typename _InitType>
class RegularSimplexConfigurator<_RealType, qc::QC_2D, _QuadType, _InitType>
      : public RegularSimplexConfiguratorBase<_RealType, qc::QC_2D, _InitType> {
public:
  RegularSimplexConfigurator ( const _InitType &Grid )
      : RegularSimplexConfiguratorBase<_RealType, qc::QC_2D, _InitType> ( Grid ), _baseFuncSet ( Grid.H() ), _volEl ( aol::Sqr ( Grid.H() ) ) {}

  typedef aol::Vector<_RealType>       VectorType;
  typedef _RealType                    RealType;
  typedef aol::Vec2<_RealType>         VecType;
  typedef aol::Vec2<_RealType>         DomVecType;
  typedef aol::Matrix22<_RealType>     MatType;
  typedef aol::SparseMatrix<_RealType> MatrixType;
  typedef aol::BaseFunctionSetLinear<_RealType, qc::QC_2D, _QuadType> BaseFuncSetType;
  typedef _QuadType                    QuadType;
  typedef qc::ScalarArray<_RealType, qc::QC_2D> ArrayType;
  typedef qc::BitArray<qc::QC_2D>               BitArrayType;

  static const int maxNumLocalDofs = 4;

  static const qc::Dimension Dim = qc::QC_2D;
  static const qc::Dimension DomDim = qc::QC_2D;

  inline int getNumLocalDofs ( const qc::Element & ) const {
    return 4;
  }

  inline RealType vol ( const qc::Element & ) const {
    return _volEl;
  }

  int maxNumQuadPoints( ) const {
    return _QuadType::numQuadPoints;
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element &/*El*/ ) const {
    return _baseFuncSet;
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
#ifdef VERBOSE
    cerr << this->_size;
#endif
    int r =
      ( El.x() + ( ( localIndex & 1 ) != 0 ) ) * 1       +
      ( El.y() + ( ( localIndex & 2 ) != 0 ) ) * this->_size[0];
#ifdef VERBOSE
    cerr << " r = " << r << " localIndex = " << localIndex << endl;
#endif
    return r;
  }

  //! returns a vec2 of global indices i,j of the dofs with number localIndex_i,j
  inline void localToGlobal ( const qc::Element &El, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal( El, localIndex0 );
    glob[1] = localToGlobal( El, localIndex1 );
  }

protected:
  BaseFuncSetType _baseFuncSet;
  const RealType  _volEl;

};


// special stuff for 3d
/**
 * \ingroup FEConfigurator
 */
template <typename _RealType, typename _QuadType, typename _InitType>
class RegularSimplexConfigurator<_RealType, qc::QC_3D, _QuadType, _InitType>
      : public RegularSimplexConfiguratorBase<_RealType, qc::QC_3D, _InitType> {
public:
  RegularSimplexConfigurator ( const _InitType &Grid )
      : RegularSimplexConfiguratorBase<_RealType, qc::QC_3D, _InitType> ( Grid ), _baseFuncSet( Grid.H()), _volEl ( aol::Cub ( Grid.H() ) ) {}
public:
  typedef aol::Vector<_RealType>       VectorType;
  typedef qc::Element                  ElementType;            //!< use the quoc element here
  typedef _RealType                    RealType;
  typedef aol::Vec3<_RealType>         DomVecType;
  typedef aol::Vec3<_RealType>         VecType;
  typedef aol::Matrix33<_RealType>     MatType;
  typedef aol::SparseMatrix<_RealType> MatrixType;
  typedef aol::BaseFunctionSetLinear<_RealType, qc::QC_3D, _QuadType> BaseFuncSetType;
  typedef qc::ScalarArray<_RealType, qc::QC_3D> ArrayType;
  typedef qc::BitArray<qc::QC_3D>               BitArrayType;
  typedef _QuadType                    QuadType;
  typedef BitArrayType                 MaskType;

  static const int maxNumLocalDofs = 8;

  static const qc::Dimension Dim = qc::QC_3D;
  static const qc::Dimension DomDim = qc::QC_3D;

  inline RealType vol ( const qc::Element & ) const {
    return _volEl;
  }

  inline int getNumLocalDofs ( const qc::Element & ) const {
    return 8;
  }

  int maxNumQuadPoints( ) const {
    return _QuadType::numQuadPoints;
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element &/*El*/ ) const {
    return _baseFuncSet;
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
#ifdef VERBOSE
    cerr << this->_size;
#endif
    int r =
      ( El.x() + ( ( localIndex & 1 ) != 0 ) ) * 1       +
      ( El.y() + ( ( localIndex & 2 ) != 0 ) ) * this->_size[0] +
      ( El.z() + ( ( localIndex & 4 ) != 0 ) ) * this->_size[0] * this->_size[1];
#ifdef VERBOSE
    cerr << " r = " << r << " localIndex = " << localIndex << endl;
#endif
    return r;
  }

  //! returns a vec2 of global indices i,j of the dofs with number localIndex_i,j
  inline void localToGlobal ( const qc::Element &El, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal( El, localIndex0 );
    glob[1] = localToGlobal( El, localIndex1 );
  }

protected:
  BaseFuncSetType _baseFuncSet;
  const RealType  _volEl;

};

}

#endif
