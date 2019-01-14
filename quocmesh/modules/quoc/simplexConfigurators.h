#ifndef __SIMPLEXCONFIGURATORS_H
#define __SIMPLEXCONFIGURATORS_H

#include <scalarArray.h>
#include <indexMapper.h>
#include <quocMatrices.h>
#include <simplexLookup.h>
#include <simplexGrid.h>
#include <simplexBaseFunctionSet.h>
#include <simplexBaseFuncSetTFE.h>
#include <gridSize.h>
#include <configurators.h>

#ifdef USE_CPP11
#define AUTOPTR std::unique_ptr
#else
#define AUTOPTR auto_ptr
#endif

namespace qc {

namespace simplex {

// --------------------------------------------------------------------------

//! inverse lexicographic mapper (fast, because using a lookup table)
//! for simplex elements
//! \author von Deylen (july 2008)
template < aol::GridGlobalIndexMode IndexMode,
           typename _CubicGridType,
           Dimension Dim>
class FastILexMapper;

// --------------------------------------------------------------------------

//! template specialization of simplex::FastILexMapper for Quoc grids
//! using a FastILexMapper for cubic elements' mapping.
//! \author von Deylen (july 2008)
template <typename _CubicGridType, Dimension Dim>
class FastILexMapper<aol::QUOC_GRID_INDEX_MODE, _CubicGridType, Dim> {
public:
  typedef typename GridStructure<_CubicGridType, Dim>::ElementType ElementType;

  FastILexMapper ( const GridStructure<_CubicGridType, Dim> & Grid )
  : _cubicIndexMapper ( Grid.getCubicGrid() )
  {}

  FastILexMapper ( const GridSize<Dim> & GridSize )
  : _cubicIndexMapper ( GridSize )
  {}

  int getElementNumber ( const ElementType &El ) const {
    // cubic Element is identified by LL dof
    int cubeNum = _cubicIndexMapper.localToGlobal ( El.getCubicElement(), 0 );
    return TopologyLookup<Dim>::numSimplexesPerCube * cubeNum + El.getSimplexNumber();
  }

  int localToGlobal ( const ElementType &El, const int localIndex ) const {
    int simNum = El.getSimplexNumber();
    int localCubicCoord = TopologyLookup<Dim>::
                               localIndicesSimplexToCube[simNum][localIndex];
    return _cubicIndexMapper.localToGlobal
                                   ( El.getCubicElement(), localCubicCoord );
  }

  void splitGlobalIndex ( int globalIndex, int & x, int & y, int & z ) const {
    _cubicIndexMapper.splitGlobalIndex ( globalIndex, x, y, z );
  }

  qc::FastILexMapper<Dim> _cubicIndexMapper;
};

// --------------------------------------------------------------------------

//! template specialization of simplex::FastILexMapper for DT grids
//! using DT grid's index mapping (which is delegated by
//! DTGrid to the DTMatrix).
//! \author von Deylen (july 2008)
template <typename _CubicGridType, Dimension Dim>
class FastILexMapper<aol::DT_GRID_INDEX_MODE, _CubicGridType, Dim> {
public:
  typedef typename GridStructure<_CubicGridType, Dim>::ElementType ElementType;

  FastILexMapper ( const GridStructure<_CubicGridType, Dim> & /*Grid*/ )
  {}

  template <typename CubicElementType>
  int localToGlobal ( const ElementType &El, const int localIndex ) const {
    // int simNum = El.getSimplexNumber();
    // int localCubicCoord = TopologyLookup<Dim>::localIndicesSimplexToCube[simNum][localIndex];
    return _CubicGridType::ConfiguratorTraitBaseType::localToGlobal
                                        ( El.getCubicElement(), localIndex );
  }
};

// --------------------------------------------------------------------------

//! base class for simplex grids' configurator
//! \author von Deylen (july 2008)
template <typename _RealType, typename _InitType>
class ConfiguratorTraitBase {
public:
  typedef _RealType                                    RealType;
  typedef _InitType                                    InitType;
  typedef aol::Vector<RealType>                        VectorType;
  typedef typename InitType::ElementType               ElementType;
  typedef typename InitType::OldFullElementIterator    ElementIteratorType;
  typedef aol::FullMatrix<RealType>                    FullMatrixType;
  typedef typename InitType::BeginIterType             BeginIterType;
  typedef typename InitType::EndIterType               EndIterType;

  ConfiguratorTraitBase ( const InitType & Grid )
  : _grid ( Grid )
  {}

  virtual ~ConfiguratorTraitBase () {}

  const BeginIterType & begin( ) const {
    return getInitializer().begin();
  }

  EndIterType end( ) const {
    return getInitializer().end();
  }

  const InitType & getInitializer () const {
    return _grid;
  }

  RealType H ( const ElementType& ) const {
    return getInitializer().H();
  }

  //! method needed in boundaryIntegration.h
  //! although the name suggests something different this implementation reflects the current behavior
  inline int localOnFaceToLocal ( const ElementType &/*El*/, const int localIndex ) const {
    return localIndex;
  }

protected:
  const InitType & _grid;
};

// --------------------------------------------------------------------------

template <typename RealType,
          qc::Dimension Dim,
          typename QuadType,
          typename _InitType>
class ConfiguratorTraitLinear {};

// --------------------------------------------------------------------------

//! \brief configurator for linear basis functions on triangles
//! \author von Deylen (july 2008)
//! \ingroup FEConfigurator
template <typename _RealType, typename _QuadType, typename _InitType>
class ConfiguratorTraitLinear<_RealType, QC_2D, _QuadType, _InitType>:
  public ConfiguratorTraitBase<_RealType, _InitType> {
public:
  typedef ConfiguratorTraitBase<_RealType, _InitType> Base;
  typedef _QuadType                          QuadType;
  typedef typename Base::RealType            RealType;
  typedef typename Base::InitType            InitType;
  typedef typename Base::VectorType          VectorType;
  typedef typename Base::ElementType         ElementType;
  typedef typename Base::ElementIteratorType ElementIteratorType;
  typedef typename Base::FullMatrixType      FullMatrixType;
  typedef typename qc::BitArray<QC_2D>       MaskType;
  typedef aol::Vec2<RealType>                VecType;
  typedef aol::BarCoord<2, RealType>         DomVecType;
  typedef aol::Matrix22<RealType>            MatType;
  typedef FastUniformGridMatrix<RealType, QC_2D> MatrixType;
  typedef typename qc::ScalarArray<RealType, qc::QC_2D> ArrayType;
  typedef BaseFunctionSetLin<RealType, QC_2D, QuadType>
                                             BaseFuncSetType;
  typedef ConfiguratorTraitLinear<RealType, QC_2D, QuadType, InitType>
                                             FullGridConfiguratorType;
  typedef ConfiguratorTraitLinear<_RealType, QC_2D, _QuadType, _InitType>
                                             Self;

  typedef qc::QuocConfiguratorTraitMultiLin<_RealType, qc::QC_2D, _QuadType >
                                             CubicConfType;                                 //arbitrary cubic configurator using the right grid

  static const bool isAdaptive = false;
  static const int maxNumLocalDofs = 3;
  static const Dimension Dim = QC_2D;
  static const Dimension DomDim = QC_2D;
  static const aol::GridGlobalIndexMode IndexMode = aol::QUOC_GRID_INDEX_MODE;

  static const Dimension DimOfWorld = QC_2D;

  explicit ConfiguratorTraitLinear ( const InitType & Grid )
  : ConfiguratorTraitBase<_RealType, _InitType> ( Grid )
  , _indexMapper ( GridSize<QC_2D> ( Grid ) )
  {
    short numSimpl = TopologyLookup<Dim>::numSimplexesPerCube;

    _volEl = aol::Sqr ( Grid.H() ) / numSimpl;

    _bfs.reserve ( numSimpl );
    for (int i = 0; i < numSimpl; ++i)
      _bfs.push_back ( BaseFuncSetType ( Grid.H(), i ) );
  }

  int maxNumQuadPoints () const {
    return QuadType::numQuadPoints;
  }

  inline static int getNumLocalDofs ( const ElementType & ) {
    return maxNumLocalDofs;
  }

  int getNumGlobalDofs () const {
    return this->getInitializer().getNumberOfNodes();
  }

  const BaseFuncSetType &
  getBaseFunctionSet ( const ElementType & El ) const {
    return _bfs [ El.getSimplexNumber() ];
  }

  const FastILexMapper<_InitType::IndexMode, typename InitType::CubicGridType,Dim>& getIndexMapper() const {
    return _indexMapper;
  }

  int getElementNumber ( const ElementType & El ) const {
    return _indexMapper.getElementNumber ( El );
  }

  int getConsecutiveElementNumber ( const ElementType & El ) const {
    int elNum = getElementNumber( El );
    return elNum - 2 * ( elNum / ( 2 * this->_grid.getNumX()) );
  }

  void getElementfromConsecutiveNumber ( int /*number*/, ElementType & /*El*/ ) const {
    throw aol::UnimplementedCodeException ( "ConfiguratorTraitLinear::getElementfromConsecutiveNumber", __FILE__, __LINE__ );
  }

  int localToGlobal ( const ElementType & El, int LocalIndex ) const {
    return _indexMapper.localToGlobal ( El, LocalIndex );
  }

  inline void localToGlobal ( const ElementType &El,
                              int localIndex0, int localIndex1,
                              aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal( El, localIndex0 );
    glob[1] = localToGlobal( El, localIndex1 );
  }

  inline void getLocalCoords ( const VecType &GlobalCoords,
                               ElementType &El,
                               DomVecType &LocalCoords ) const {
    VecType localCoords;

    // compute the cubic element
    for ( int i = 0; i < QC_2D; ++i ) {
      localCoords[i] = GlobalCoords[i] / this->_grid.H();
      El.getCubicElement()[i] = static_cast<short>( localCoords[i] );
      localCoords[i] -= El.getCubicElement()[i];
    }

    // compute simplex
    El.set( El.getCubicElement(), localCoords[0] + localCoords[1] > 1. ? 1 : 0 );

    getBarycentricCoords ( localCoords, El, LocalCoords );
  }

  inline void getBarycentricCoords ( const VecType& CartesianCoords,
                                     const ElementType& El,
                                     DomVecType &LocalCoords ) const {
    // compute barycentric coordinates
     if ( El.getSimplexNumber() == 0 ) {
       LocalCoords [0] = 1 - CartesianCoords[0] - CartesianCoords[1];
       LocalCoords [1] = CartesianCoords[0];
       LocalCoords [2] = CartesianCoords[1];
     } else {
       LocalCoords [0] = 1 - CartesianCoords[1];
       LocalCoords [1] = 1 - CartesianCoords[0];
       LocalCoords [2] = CartesianCoords[0] + CartesianCoords[1] - 1;
     }
  }

  inline void getGlobalCoords ( const ElementType &El,
                                const DomVecType &LocalCoords,
                                VecType &GlobalCoords ) const {
    // retrieve coordinates of cubic element
    for ( int i = 0; i < 2; i++ )
      GlobalCoords[i] = static_cast<RealType>( El.getCubicElement()[i] );

    // retrieve coordinates within a simplex
    if ( El.getSimplexNumber() == 0 ) {
      GlobalCoords[0] += LocalCoords [1];
      GlobalCoords[1] += LocalCoords [2];
    } else {
      GlobalCoords[0] += 1 - LocalCoords [1];
      GlobalCoords[1] += 1 - LocalCoords [0];
    }

    // scale by grid size
    GlobalCoords *= this->_grid.H();
  }

  RealType vol ( const ElementType & ) const {
    return _volEl;
  }

  MatrixType* createNewMatrix( ) const {
    MatrixType *mat =
      new MatrixType ( this->getInitializer().getCubicGrid() );
    return mat;
  }

  void writeNodeExistMaskTo ( MaskType & existMask ) const {
    existMask.setAll ( false );
    for ( ElementIteratorType iter = this->begin(); iter != this->end(); ++iter ) {
      for ( int i = 0; i < getNumLocalDofs ( *iter ); ++i )
        existMask.set ( localToGlobal ( *iter, i ), true );
    }
  }

  using ConfiguratorTraitBase<_RealType, _InitType>::localOnFaceToLocal;

  template <bool ClipCoord>
  inline bool transformCoord ( const ElementType &/*El*/, const DomVecType &/*RefCoord*/, const VecType &/*Offset*/, ElementType &/*TransformedEl*/, DomVecType &/*TransformedLocalCoord*/ ) const {
    throw aol::UnimplementedCodeException ( "ConfiguratorTraitLinear::transformCoord", __FILE__, __LINE__ );
    return false;
  }
protected:
  FastILexMapper<_InitType::IndexMode,
                 typename InitType::CubicGridType,Dim>    _indexMapper;
  std::vector<BaseFuncSetType>                            _bfs;
  RealType                                                _volEl;
};

// --------------------------------------------------------------------------

//! \brief configurator for linear basis functions on tetrahedra
//! \author von Deylen (july 2008)
//! \ingroup FEConfigurator
template <typename _RealType, typename _QuadType, typename _InitType>
class ConfiguratorTraitLinear<_RealType, QC_3D, _QuadType, _InitType>
  : public ConfiguratorTraitBase<_RealType, _InitType> {
public:
  typedef ConfiguratorTraitBase<_RealType, _InitType> Base;
  typedef _QuadType                          QuadType;
  typedef typename Base::RealType            RealType;
  typedef typename Base::InitType            InitType;
  typedef typename Base::VectorType          VectorType;
  typedef typename Base::ElementType         ElementType;
  typedef typename Base::ElementIteratorType ElementIteratorType;
  typedef typename Base::FullMatrixType      FullMatrixType;
  typedef typename qc::BitArray<QC_3D>       MaskType;
  typedef aol::Vec3<RealType>                VecType;
  typedef aol::BarCoord<3, RealType>         DomVecType;
  typedef aol::Matrix33<RealType>            MatType;
  typedef FastUniformGridMatrix<RealType, QC_3D> MatrixType;
  typedef typename qc::ScalarArray<RealType, qc::QC_3D> ArrayType;
  typedef BaseFunctionSetLin<RealType, QC_3D, QuadType>
                                             BaseFuncSetType;
  typedef ConfiguratorTraitLinear<RealType, QC_3D, QuadType, InitType>
                                             FullGridConfiguratorType;
  typedef ConfiguratorTraitLinear<_RealType, QC_3D, _QuadType, _InitType>
                                             Self;

  typedef qc::QuocConfiguratorTraitMultiLin<_RealType, qc::QC_3D, _QuadType >
                                             CubicConfType;                                 //arbitrary cubic configurator using the right grid

  static const bool isAdaptive = false;
  static const int maxNumLocalDofs = 4;
  static const Dimension Dim = QC_3D;
  static const Dimension DimOfWorld = QC_3D;
  static const Dimension DomDim = QC_3D;
  static const aol::GridGlobalIndexMode IndexMode = InitType::IndexMode;


  explicit ConfiguratorTraitLinear ( const InitType & Grid )
  : ConfiguratorTraitBase<_RealType, _InitType> ( Grid )
  , _indexMapper ( GridSize<QC_3D> ( Grid ) )
  {
    short numSimpl = TopologyLookup<Dim>::numSimplexesPerCube;

    _volEl = aol::Cub ( Grid.H() ) / numSimpl;

    _bfs.reserve ( numSimpl );
    for (int i = 0; i < numSimpl; ++i)
      _bfs.push_back ( BaseFuncSetType ( Grid.H(), i ) );
  }

  int maxNumQuadPoints () const {
    return QuadType::numQuadPoints;
  }

  inline static int getNumLocalDofs ( const ElementType & ) {
    return maxNumLocalDofs;
  }

  inline static int getNumLocalDofs ( const SimplexBoundaryFaceElement<typename InitType::CubicGridType::ElementType, RealType, qc::QC_3D> & ) {
    return maxNumLocalDofs;
  }

  int getNumGlobalDofs () const {
    return this->getInitializer().getNumberOfNodes();
  }

  const BaseFuncSetType &
  getBaseFunctionSet ( const ElementType & El ) const {
    return _bfs [ El.getSimplexNumber() ];
  }

  int localToGlobal ( const SimplexBoundaryFaceElement<typename InitType::CubicGridType::ElementType, RealType, qc::QC_3D>& El, int LocalIndex ) const {
    return _indexMapper.localToGlobal ( El, LocalIndex);
    }

  int localToGlobal ( const ElementType & El, int LocalIndex ) const {
    return _indexMapper.localToGlobal ( El, LocalIndex );
  }

  inline void localToGlobal ( const ElementType &El,
                              int localIndex0, int localIndex1,
                              aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal( El, localIndex0 );
    glob[1] = localToGlobal( El, localIndex1 );
  }

  inline void getLocalCoords ( const VecType &GlobalCoords,
                               ElementType &El,
                               DomVecType &LocalCoords ) const {
    // compute the cubic element
    VecType localCoords;
    for ( int i = 0; i < QC_3D; ++i ) {
      localCoords[i] = GlobalCoords[i] / this->_grid.H();
      El.getCubicElement()[i] = static_cast<short>( localCoords[i] );
      localCoords[i] -= El.getCubicElement()[i];
    }

    // compute simplex
    short simplexCode = localCoords[0] + localCoords[1] > 1. ? 4 : 0;
    simplexCode |= localCoords[1] > localCoords[2] ? 2 : 0;
    simplexCode |= localCoords[0] + localCoords[1] - localCoords[2] > simplexCode >> 2 ? 1 : 0;
    const static short simplices[8] = { 2, 1, -1, 0, 5, -1, 4, 3 };
    El.set( El.getCubicElement(), simplices[simplexCode] );

    getBarycentricCoords ( localCoords, El, LocalCoords );
  }

  inline void getBarycentricCoords ( const VecType& CartesianCoords,
                                     const ElementType& El,
                                     DomVecType &LocalCoords ) const {
    // compute barycentric coordinates
    switch ( El.getSimplexNumber() ) {
    case 0:
      LocalCoords [0] = 1 - CartesianCoords[0] - CartesianCoords[1];
      LocalCoords [1] = CartesianCoords[0];
      LocalCoords [2] = CartesianCoords[1] - CartesianCoords[2];
      LocalCoords [3] = CartesianCoords[2];
      break;
    case 1:
      LocalCoords [0] = CartesianCoords[0] + CartesianCoords[1] - CartesianCoords[2];
      LocalCoords [1] = 1 - CartesianCoords[0] - CartesianCoords[1];
      LocalCoords [2] = CartesianCoords[2] - CartesianCoords[1];
      LocalCoords [3] = CartesianCoords[1];
      break;
    case 2:
      LocalCoords [0] = 1 - CartesianCoords[2];
      LocalCoords [1] = CartesianCoords[2] - CartesianCoords[0] - CartesianCoords[1];
      LocalCoords [2] = CartesianCoords[0];
      LocalCoords [3] = CartesianCoords[1];
      break;
    case 3:
      LocalCoords [0] = 1 - CartesianCoords[0];
      LocalCoords [1] = 1 - CartesianCoords[1];
      LocalCoords [2] = CartesianCoords[0] + CartesianCoords[1] - CartesianCoords[2] - 1;
      LocalCoords [3] = CartesianCoords[2];
      break;
    case 4:
      LocalCoords [0] = 1 - CartesianCoords[1];
      LocalCoords [1] = CartesianCoords[1] - CartesianCoords[2];
      LocalCoords [2] = 1 + CartesianCoords[2] - CartesianCoords[0] - CartesianCoords[1];
      LocalCoords [3] = CartesianCoords[0] + CartesianCoords[1] - 1;
      break;
    case 5:
      LocalCoords [0] = CartesianCoords[2] - CartesianCoords[1];
      LocalCoords [1] = 1 - CartesianCoords[2];
      LocalCoords [2] = 1 - CartesianCoords[0];
      LocalCoords [3] = CartesianCoords[0] + CartesianCoords[1] - 1;
      break;
    default:
      throw aol::Exception ( "Wrong simplex number!", __FILE__, __LINE__ );
      break;
    }
  }

  inline void getGlobalCoords ( const ElementType &El,
                                const DomVecType &LocalCoords,
                                VecType &GlobalCoords ) const {
    // retrieve coordinates of cubic element
    for ( int i = 0; i < 3; i++ )
      GlobalCoords[i] = static_cast<RealType>( El.getCubicElement()[i] );

    // retrieve coordinates within a simplex
    switch ( El.getSimplexNumber() ) {
    case 0:
      GlobalCoords[0] += LocalCoords[1];
      GlobalCoords[1] += LocalCoords[2] + LocalCoords[3];
      GlobalCoords[2] += LocalCoords[3];
      break;
    case 1:
      GlobalCoords[0] += LocalCoords[0] + LocalCoords[2];
      GlobalCoords[1] += LocalCoords[3];
      GlobalCoords[2] += LocalCoords[2] + LocalCoords[3];
      break;
    case 2:
      GlobalCoords[0] += LocalCoords[2];
      GlobalCoords[1] += LocalCoords[3];
      GlobalCoords[2] += 1 - LocalCoords[0];
      break;
    case 3:
      GlobalCoords[0] += 1 - LocalCoords[0];
      GlobalCoords[1] += 1 - LocalCoords[1];
      GlobalCoords[2] += LocalCoords[3];
      break;
    case 4:
      GlobalCoords[0] += LocalCoords[0] + LocalCoords[3];
      GlobalCoords[1] += 1 - LocalCoords[0];
      GlobalCoords[2] += LocalCoords[2] + LocalCoords[3];
      break;
    case 5:
      GlobalCoords[0] += 1 - LocalCoords[2];
      GlobalCoords[1] += LocalCoords[2] + LocalCoords[3];
      GlobalCoords[2] += 1 - LocalCoords[1];
      break;
    default:
      throw aol::Exception ( "Wrong simplex number!", __FILE__, __LINE__ );
      break;
    }

    // scale by grid size
    GlobalCoords *= this->_grid.H();
  }

  RealType vol ( const ElementType & ) const {
    return _volEl;
  }

  const FastILexMapper<_InitType::IndexMode, typename InitType::CubicGridType,Dim>& getIndexMapper() const {
    return _indexMapper;
  }

  MatrixType* createNewMatrix( ) const {
    MatrixType *mat =
      new MatrixType ( this->getInitializer().getCubicGrid() );
    return mat;
  }

  void writeNodeExistMaskTo ( MaskType & existMask ) const {
    existMask.setAll ( false );
    for ( ElementIteratorType iter = this->begin(); iter != this->end(); ++iter ) {
      for ( int i = 0; i < getNumLocalDofs ( *iter ); ++i )
        existMask.set ( localToGlobal ( *iter, i ), true );
    }
  }

  //! returns consecutive Element Number (ST)
  int getConsecutiveElementNumber ( const ElementType & El ) const {
    int cElNum = _indexMapper._cubicIndexMapper.localToGlobal(El.getCubicElement(), 0);                   //not yet consecutive!
    int pFactor = cElNum / ( this->_grid.getNumX() * this->_grid.getNumY() ) * ( this->_grid.getNumX() - 1 );         //n-1 elements per row

    int consecutiveCubicElNum = cElNum - cElNum / this->_grid.getNumX() - pFactor;                        //now it is

    return consecutiveCubicElNum * TopologyLookup<Dim>::numSimplexesPerCube + El.getSimplexNumber();
  }

  using ConfiguratorTraitBase<_RealType, _InitType>::localOnFaceToLocal;

protected:
  FastILexMapper<_InitType::IndexMode,
                 typename InitType::CubicGridType, Dim>   _indexMapper;
  std::vector<BaseFuncSetType>                            _bfs;
  RealType                                                _volEl;
};

// --------------------------------------------------------------------------

template <typename RealType,
          qc::Dimension Dim,
          typename QuadType,
          typename InitType>
class ConfiguratorTraitTFE {};

// --------------------------------------------------------------------------

//! \brief configurator for cut-off linear basis functions on tetrahedra
//! \author von Deylen (july 2008)
//! \ingroup FEConfigurator
template <typename _RealType, typename _QuadType, typename _InitType>
class ConfiguratorTraitTFE<_RealType, QC_3D, _QuadType, _InitType>
  : public ConfiguratorTraitBase<_RealType, _InitType> {
public:
  typedef ConfiguratorTraitBase<_RealType, _InitType>   Base;
  typedef _QuadType                                     QuadType;
  typedef typename Base::RealType                       RealType;
  typedef typename Base::InitType                       InitType;
  typedef typename Base::VectorType                     VectorType;
  typedef typename Base::ElementType                    ElementType;
  typedef typename Base::ElementIteratorType            ElementIteratorType;
  typedef typename Base::FullMatrixType                 FullMatrixType;
  typedef typename qc::BitArray<QC_3D>                  MaskType;
  typedef aol::Vec3<RealType>                           VecType;
  typedef aol::BarCoord<QC_3D, RealType>                DomVecType;
  typedef aol::Matrix33<RealType>                       MatType;
  typedef MultilinFEBandMatrix<RealType, QC_3D>         MatrixType;
  typedef typename qc::ScalarArray<RealType, qc::QC_3D> ArrayType;
  typedef BaseFunctionSetTFE<RealType, QC_3D, QuadType> BaseFuncSetType;
  typedef ConfiguratorTraitTFE<RealType, QC_3D, QuadType, InitType>
                                                        FullGridConfiguratorType;
  typedef FastILexMapper<_InitType::IndexMode,
      typename InitType::CubicGridType, QC_3D>          IndexMapperType;
  typedef BaseFunctionCollection<RealType, QuadType, QC_3D> BFSCollection;

  static const int maxNumLocalDofs  = 4;
  static const Dimension Dim        = QC_3D;
  static const Dimension DimOfWorld = QC_3D;
  static const Dimension DomDim     = QC_3D;
  static const aol::GridGlobalIndexMode IndexMode = InitType::IndexMode;

  //! usual constructor
  ConfiguratorTraitTFE ( const InitType & Grid, const VectorType & LevelValues,
                         RealType BandRadius )
      : ConfiguratorTraitBase<_RealType, _InitType> ( Grid )
      , _indexMapper ( GridSize<QC_3D> ( Grid ) )
      , _bfsCollection ( Grid, LevelValues, BandRadius, _indexMapper )
  {}

  //! copy constructor
  //!
  //! copies everything like the implicit copy constructor,
  //! only the Mask pointer is not copied and will be newly
  //! created when accessed on the copy for the first time.
  //! \note Be aware that copying this costructor is quite
  //! expensive, due to its base function set collection.
  explicit ConfiguratorTraitTFE ( const ConfiguratorTraitTFE<RealType, Dim, QuadType, InitType> & other )
      : ConfiguratorTraitBase<_RealType, _InitType> ( other )
      , _indexMapper ( other._indexMapper )
      , _bfsCollection ( other._bfsCollection )
  {}

  int maxNumQuadPoints () const {
    return 3;
  }

  inline static int getNumLocalDofs ( const ElementType & ) {
    return maxNumLocalDofs;
  }

  int getNumGlobalDofs () const {
    return this->getInitializer().getNumberOfNodes();
  }

  void splitGlobalIndex ( int globalIndex, int & x, int & y, int & z ) const {
    _indexMapper.splitGlobalIndex ( globalIndex, x, y, z );
  }


  const BaseFuncSetType &
  getBaseFunctionSet ( const ElementType & El ) const {
    return _bfsCollection.getBaseFunctionSet
            ( this->getInitializer().getElementIndex ( El ) );
  }

  int localToGlobal ( const ElementType & El, int LocalIndex ) const {
    return _indexMapper.localToGlobal ( El, LocalIndex );
  }

  inline void localToGlobal ( const ElementType &El,
                              int localIndex0, int localIndex1,
                              aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal( El, localIndex0 );
    glob[1] = localToGlobal( El, localIndex1 );
  }

  RealType vol ( const ElementType & El ) const {
    return getBaseFunctionSet( El ).getVolume();
  }

  MatrixType* createNewMatrix( ) const {
    MatrixType *mat =
      new MatrixType ( this->getInitializer().getCubicGrid() );
    return mat;
  }

  void writeNodeExistMaskTo ( MaskType & existMask ) const {
    existMask.setAll ( false );
    for ( ElementIteratorType iter = this->begin(); iter != this->end(); ++iter ) {
      for ( int i = 0; i < getNumLocalDofs ( *iter ); ++i )
        existMask.set ( localToGlobal ( *iter, i ), true );
    }
  }

  const MaskType & getNodeExistMask () const {
    _nodeExistMask.reset ( new MaskType ( this->getInitializer().getCubicGrid() ) );
    writeNodeExistMaskTo ( *_nodeExistMask );
    return *_nodeExistMask;
  }

protected:
  IndexMapperType                _indexMapper;
  BFSCollection                  _bfsCollection;
  mutable AUTOPTR<MaskType>      _nodeExistMask;
};

// --------------------------------------------------------------------------

} // end of namespace simplex.

} // end of namespace qc.

#endif
