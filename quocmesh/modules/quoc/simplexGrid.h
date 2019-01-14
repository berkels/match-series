#ifndef __SIMPLEXGRID_H
#define __SIMPLEXGRID_H

#include <quoc.h>
#include <gridBase.h>
#include <gridSize.h>
#include <simplexLookup.h>
#include <iterators.h>


namespace qc {

namespace simplex {

// ------------------------------------------------------------------------------------------------

/** simplicial element.
 *
 *  The class is templatized over the underlying cubic element,
 *  storing a cubic element and a simplex number.
 *
 *  \author von Deylen (july 2008)
 */
template <typename CubicElementType>
class Element {
public:
  Element() : _simplexNumber ( 0 ) {}
  Element ( const CubicElementType & CubicElement, short SimplexNumber )
  : _cubicEl ( CubicElement )
  , _simplexNumber ( SimplexNumber )
  {}

  Element<CubicElementType> &operator= ( const Element<CubicElementType> & El ) {
    _cubicEl = El.getCubicElement();
    _simplexNumber = El.getSimplexNumber();
    return *this;
  }
  const CubicElementType & getCubicElement () const {
    return _cubicEl;
  }

  CubicElementType & getCubicElement () {
    return _cubicEl;
  }

  short getSimplexNumber () const {
    return _simplexNumber;
  }

  void set ( const CubicElementType & CubicElement, short SimplexNumber ) {
    _cubicEl.set ( CubicElement );
    _simplexNumber = SimplexNumber;
  }

  short & simplexNumber ( ) {
    return _simplexNumber;
  }

  bool operator== ( const Element<CubicElementType> & El ) const {
    return (   getCubicElement() == El.getCubicElement()
            && getSimplexNumber() == El.getSimplexNumber() );
  }

  bool operator!= ( const Element<CubicElementType> & El ) const {
    return (   getCubicElement() != El.getCubicElement()
            || getSimplexNumber() != El.getSimplexNumber() );
  }

  unsigned char type() const {
    throw aol::Exception ( "qc::simplex::ElementType::type() is not implemented!", __FILE__, __LINE__ );
    return 23;
  }

  short level() const {
    throw aol::Exception ( "qc::simplex::ElementType::level() is not implemented!", __FILE__, __LINE__ );
    return -1;
  }

  short operator[] ( const int I ) const {
    return getCubicElement()[I];
  }

protected:
  CubicElementType _cubicEl;
  short _simplexNumber;
};

template < typename CubicElementType >
std::ostream& operator<< ( std::ostream& out, const Element< CubicElementType >& El ) {
  out << El.getCubicElement() << "\tsimplex number " << El.getSimplexNumber();
  return out;
}

// ------------------------------------------------------------------------------------------------

//! \brief Simplicial face element
template<typename CubicElementType, typename RealType, qc::Dimension Dim>
class SimplexBoundaryFaceElement {};

 //! \brief Simplicial face element implementation for 3D
template<typename CubicElementType, typename RealType>
class SimplexBoundaryFaceElement<CubicElementType, RealType, qc::QC_3D> : public Element<CubicElementType>
{
public:
  //! \brief SimplexBoundaryFaceElement::X_LOWER_BOUNDARY denotes the y-z-plane with x = 0, SimplexBoundaryFaceElement::X_UPPER_BOUNDARY the y-z-plane with x = 1 and so on
  enum BoundaryFaceType { X_LOWER_BOUNDARY = 0, ///< y-z-plane with x = 0,
                          X_UPPER_BOUNDARY = 1, ///< y-z-plane with x = 1
                          Y_LOWER_BOUNDARY = 2,
                          Y_UPPER_BOUNDARY = 3,
                          Z_LOWER_BOUNDARY = 4,
                          Z_UPPER_BOUNDARY = 5,
                          UNDEF_BOUNDARY = 42   };

  // Lookup bfType, Triangle, Node
  // assigns local 2D nodes to local 3D nodes
  static const short simplexNodeLookup[6][2][3];

  // Lookup bfType, simplexFaceNumber
  // assigns 2D elements to 3D elements
  static const short fullElementSimplexNumberLookup[6][2];

  //! standart constructor
  SimplexBoundaryFaceElement()
  : _cubicFace(NULL), _bfType(UNDEF_BOUNDARY), _simplexFaceNumber(-1)
    {}

  //! \brief Construct simplicial face element from cubic face element
  SimplexBoundaryFaceElement(qc::BoundaryFaceElement<RealType, qc::QC_3D>& cubicFace, BoundaryFaceType bfType, short simplexFaceNumber)
   : _cubicFace(cubicFace), _bfType(bfType), _simplexFaceNumber(simplexFaceNumber)
  {
    this->_cubicEl.set(cubicFace);
    this->_simplexNumber = fullElementSimplexNumberLookup[static_cast<int>(this->_bfType)][_simplexFaceNumber];
  }

  void set(const qc::BoundaryFaceElement<RealType, qc::QC_3D>& cubicFace, short simplexFaceNumber)
  {
    this->_cubicEl.set(cubicFace);
    _cubicFace = &cubicFace;
    _bfType = static_cast<BoundaryFaceType>(_cubicFace->getBoundaryFaceType());
    _simplexFaceNumber = simplexFaceNumber;
    this->_simplexNumber = fullElementSimplexNumberLookup[static_cast<int>(this->_bfType)][_simplexFaceNumber];
  }

  aol::Vec3<RealType> getNormal()
  {
    return _cubicFace->getNormal();
  }

  aol::BarCoord<3, RealType> getRefCoordOnElementFromRefCoordOnFace(const aol::BarCoord<2, RealType>& barCoord) const
  {
    aol::BarCoord<3, RealType> barCoord3D;
    barCoord3D.setZero();
    for(int i=0; i<3; ++i)
    {
      barCoord3D[ simplexNodeLookup[static_cast<int>(_bfType)][_simplexFaceNumber][i] ] = barCoord[i];
    }

    return barCoord3D;
  }

  //! Converts local boundary face nodes to nodes on the element.
  short getLocalNodeOnBoundary(short i) const
  {
    if( (i < 0) || (i > 2) )
      throw aol::Exception("wrong choice of i, it must be 0 <= i <= 2", __FILE__, __LINE__);

    return simplexNodeLookup[static_cast<int>(this->_bfType)][_simplexFaceNumber][i];
  }

  //! Sets this to the next simplex Element on the same face. Only works, if there is a next element.
  void setToNextSimplexFace()
  {
    if(_simplexFaceNumber > 0)
      {
  throw aol::Exception("Cannot go to next simplex face, as this is already the last one.", __FILE__, __LINE__);
      }

    else
      {
  ++_simplexFaceNumber;
  this->set(*_cubicFace, _simplexFaceNumber);
      }
  }

  short getSimplexNumber() const
  {
    return this->_simplexNumber;
  }

  short getSimplexFaceNumber() const
  {
    return _simplexFaceNumber;
  }

  BoundaryFaceType getBoundaryFaceType()
  {
    return _bfType;
  }

  int level() const
  {
    return _cubicFace->level();
  }

 protected:
  const qc::BoundaryFaceElement<RealType, qc::QC_3D>* _cubicFace;
  BoundaryFaceType _bfType;
  short _simplexFaceNumber;
};

template <typename CubicElementType, typename RealType>
const short SimplexBoundaryFaceElement<CubicElementType, RealType, qc::QC_3D>::simplexNodeLookup[6][2][3] = {  { {0, 2, 3}, {0, 1, 3} },
                                                                               { {1, 2, 3}, {0, 1, 3} },
                                                                           { {0, 1, 2}, {0, 1, 2} },
                                                                           { {0, 2, 3}, {1, 2, 3} },
                                                                           { {0, 1, 2}, {0, 1, 2} },
                                                                           { {1, 2, 3}, {0, 2, 3} }  };


template <typename CubicElementType, typename RealType>
const short SimplexBoundaryFaceElement<CubicElementType, RealType, qc::QC_3D>::fullElementSimplexNumberLookup[6][2] = { {0, 2},
                              {3, 5},
                              {1, 2},
                              {3, 4},
                              {0, 3},
                              {2, 5} };

// ------------------------------------------------------------------------------------------------

/** simplicial grid structure without iterators.
 *
 *  A simplicial grid is always constructed on top of a cubic
 *  element grid, dividing each cube in the appropriate number
 *  of simplices.
 *  Inheritance from _CubicGridType is protected to model an
 *  "is implemented as" structure, but no "is a" relation.
 *
 *  \author von Deylen (july 2008)
 */
template <typename _CubicGridType>
class GridStructureBase : protected _CubicGridType {
public:
  static const aol::GridGlobalIndexMode IndexMode = _CubicGridType::IndexMode;

  explicit GridStructureBase ( const _CubicGridType & cubicGrid )
  : _CubicGridType ( cubicGrid )
  {}

  template <qc::Dimension DimOfWorld>
  explicit GridStructureBase ( const GridSize<DimOfWorld> & gridSize )
  : _CubicGridType ( gridSize )
  {}

  const _CubicGridType & getCubicGrid () const {
    return static_cast<const _CubicGridType &> (*this);
  }

  _CubicGridType & getCubicGrid () {
    return static_cast<_CubicGridType &> (*this);
  }

  // functions from qc::GridStructure
  using _CubicGridType::getSize;
  using _CubicGridType::getNumX;
  using _CubicGridType::getNumY;
  using _CubicGridType::getNumZ;
  using _CubicGridType::getDimOfWorld;

  // functions from qc::GridDefinition
  using _CubicGridType::H;
  using _CubicGridType::getGridDepth;
  using _CubicGridType::getNumberOfNodes;
  using _CubicGridType::getNumberOfDofs;
  using _CubicGridType::getNumberOfBoundaryNodes;

  using _CubicGridType::isAdaptive;

  using _CubicGridType::checkForHangingNode;

  virtual int checkForHangingNode ( const qc::simplex::Element< qc::Element > &, int ) const {
    return -1;
  }

};

// ------------------------------------------------------------------------------------------------

/** base class for simplicial iterators
 *
 *  all simplicial iterators follow the iterator idea "initialization
 *  via grid, termination via EndElement object" (just as
 *  GridDefinition::OldFullElementIterator, opposite to
 *  GridDefinition::iterator).
 *
 *  \author von Deylen (july 2008)
 */
template <typename IteratedType, typename _CubicGridType, typename Imp>
class Iterator {

public:
  typedef Iterator<IteratedType, _CubicGridType, Imp> Self;

  Self & operator= ( const _CubicGridType & Grid ) {
    asImp() = Grid;
  }
  bool operator== ( const EndElement & End ) const {
    return (asImp() == End );
  }
  bool operator!= ( const EndElement & End ) const {
    return !(*this == End);
  }

  Self & operator++ () {
    ++(asImp());
  }
  Self & operator++ (int) {
    (asImp())++;
  }

  IteratedType & operator* () {
    return _cur;
  }
  IteratedType * operator->() {
    return &_cur;
  }

protected:
  Imp & asImp () {
    return static_cast<Imp&> (*this);
  }

  const Imp & asImp () const {
    return static_cast<const Imp&> (*this);
  }

  IteratedType _cur;
};

// ------------------------------------------------------------------------------------------------

/** iterator over all simplices in a cubic element
 *
 *  As simplicial grids are "refined cubic grids", a natural
 *  idea for traversing all elements is to run (a) through the
 *  cubic elements, (b) in each cube through all simplices.
 *
 *  \author von Deylen (july 2008)
 */
template <typename CubicElIterator, typename _CubicGridType, qc::Dimension Dim>
class SimplexInCubeIterator :
 public Iterator< Element<typename CubicElIterator::IteratedType>,
                   _CubicGridType,
                   SimplexInCubeIterator<CubicElIterator, _CubicGridType, Dim> > {
public:
  typedef SimplexInCubeIterator<CubicElIterator,
                                 _CubicGridType, Dim> Self;
  typedef TopologyLookup<Dim>                         Lookup;
  typedef typename CubicElIterator::IteratedType      CubicElType;
  typedef Element<CubicElType>                        IteratedType;
  typedef _CubicGridType                              BeginType;
  typedef EndElement                                  EndType;

  SimplexInCubeIterator ( ) {}

  SimplexInCubeIterator ( const BeginType & Grid ) {
    *this = Grid;
  }

  Self & operator= ( const BeginType & ) {
    return *this;
  }

  Self & operator= ( const CubicElType & CubicEl ) {
    this->_cur.set ( CubicEl, 0 );
    return *this;
  }

  bool operator== ( const EndType & ) const {
    return atEnd();
  }

  bool atEnd () const {
    return ( this->_cur.getSimplexNumber() == Lookup::numSimplexesPerCube );
  }

  Self & operator++ () {
    ++(this->_cur.simplexNumber());
    return *this;
  }

  Self operator++ (int) {
    Self tmp ( *this );
    ++(*this);
    return tmp;
  }
};

// ------------------------------------------------------------------------------------------------

/** simplicial grids' element iterator
 *
 *  This iterator combines a cube iterator (known from the
 *  template argument "CubicElIterator") and a SimplexInCubeIterator.
 *  Therefore, no assumptions on the elements' order in the cubic
 *  grid is needed, but delegated to the cubic grid's element
 *  iterator.
 *
 *  \author von Deylen (july 2008)
 */
template <typename CubicElIterator, typename _CubicGridType, qc::Dimension Dim>
class FullSimplexIterator :
  public Iterator< Element<typename CubicElIterator::IteratedType>,
                   _CubicGridType,
                   FullSimplexIterator<CubicElIterator, _CubicGridType, Dim> > {

public:
  typedef FullSimplexIterator<CubicElIterator,
                               _CubicGridType, Dim>   Self;
  typedef typename CubicElIterator::IteratedType      CubicElType;
  typedef Element<CubicElType>                        IteratedType;
  typedef _CubicGridType                              BeginType;
  typedef EndElement                                  EndType;

  FullSimplexIterator ( ) {}

  FullSimplexIterator ( const BeginType & Grid ) {
    *this = Grid;
  }

  Self & operator= ( const BeginType & Grid ) {
    _cubicIter = Grid;
    _simplexIter = *_cubicIter;
    return *this;
  }

  bool operator== ( const EndType & ) const {
    bool cubicAtEnd = _cubicIter.atEnd();
    bool simplexAtEnd = _simplexIter.atEnd();
    return ( cubicAtEnd && simplexAtEnd );
  }

  Self & operator++ () {
    ++_simplexIter;
    if ( _simplexIter.atEnd() ) {
      ++_cubicIter;
      if ( !_cubicIter.atEnd() )
        _simplexIter = *_cubicIter;
    }
    return *this;
  }

  Self operator++ (int) {
    Self tmp ( *this );
    ++(*this);
    return tmp;
  }

  IteratedType & operator* () {
    return _simplexIter.operator*();
  }

  IteratedType * operator->() {
    return _simplexIter.operator->();
  }

protected:
  CubicElIterator _cubicIter;
  SimplexInCubeIterator<CubicElIterator, _CubicGridType, Dim> _simplexIter;
};

// ------------------------------------------------------------------------------------------------

//! \brief Iterates over all simplices of a cubic boundary face
//! \warning Only 3D!
template <typename CubicBoundaryFaceElementIterator, typename _GridType, /*typename _CubicConfType,*/ qc::Dimension Dim, typename RealType = double>
class SimplexInFaceIterator
: public Iterator< SimplexBoundaryFaceElement<typename _GridType::CubicGridType::ElementType, RealType, Dim>, _GridType,
  SimplexInFaceIterator<CubicBoundaryFaceElementIterator, _GridType, /*_CubicConfType, */Dim, RealType> >
{
public:
  typedef SimplexInFaceIterator <CubicBoundaryFaceElementIterator, _GridType, /*_CubicConfType,*/ Dim, RealType>
                                            Self;
  typedef TopologyLookup<Dim>                Lookup;
  typedef _GridType                          BeginType;
  typedef EndElement                         EndType;

  SimplexInFaceIterator()
  : _numSimplexesPerFace(2)
  {}


  SimplexInFaceIterator(const BeginType& grid)
  : _numSimplexesPerFace(2)
  {
    *this = grid;
  }

  void setGrid(const BeginType& grid)
  {
    *this = grid;
  }

  SimplexInFaceIterator(const BoundaryFaceElement<RealType, Dim>& cubicBoundaryElement)
  : _numSimplexesPerFace(2)
  {
    this->_cur.set(cubicBoundaryElement, 0);
  }

  Self& operator=(const BeginType&)
  {
    return *this;
  }

  Self& operator=(const BoundaryFaceElement<RealType, Dim>& cubicBoundaryElement)
  {
    this->_cur.set(cubicBoundaryElement, 0);

    return *this;
  }

  bool operator==(const EndType&) const
  {
    return atEnd();
  }

  bool atEnd() const
  {
    return ( this->_cur.getSimplexFaceNumber() == (_numSimplexesPerFace - 1) );
  }

  Self& operator++()
  {
    this->_cur.setToNextSimplexFace();

    return *this;
  }

  Self operator++(int)
  {
    Self tmp(*this);
    ++(*this);

    return tmp;
  }

private:
  short _numSimplexesPerFace;
};


template <typename CubicBoundaryFaceElementIterator, typename _GridType, typename _SimplexConfType,  qc::Dimension Dim>
class SimplexBoundaryFaceElementIterator :
public Iterator< SimplexBoundaryFaceElement<typename _GridType::CubicGridType::ElementType, typename _SimplexConfType::RealType, Dim>, _GridType,
  SimplexBoundaryFaceElementIterator<CubicBoundaryFaceElementIterator, _GridType, _SimplexConfType, Dim> >
  {
public:
  typedef typename _SimplexConfType::RealType                             RealType;
  typedef SimplexBoundaryFaceElementIterator<CubicBoundaryFaceElementIterator, _GridType, _SimplexConfType, Dim>
                                                                              Self;
  typedef TopologyLookup<Dim>                                             Lookup;
  typedef SimplexBoundaryFaceElement<typename _GridType::CubicGridType::ElementType, /*_CubicConfType::RealType*/RealType, qc::QC_3D> IteratedType;   //warning: Usage of RealType not really consistent with declaration of SimplexBoundaryFaceElement
  typedef _GridType                                                       BeginType;
  typedef EndElement                                                      EndType;

  SimplexBoundaryFaceElementIterator()
  : _cubicFaceIter( NULL )
  {}

  SimplexBoundaryFaceElementIterator(const BeginType& grid)
  : _simplexFaceIter(grid)
  {
    _cubicFaceIter = NULL;
    *this = grid;
    _cubicFaceIter = new CubicBoundaryFaceElementIterator(grid);
    _simplexFaceIter = _cubicFaceIter->operator*();
  }

  ~SimplexBoundaryFaceElementIterator()
   {
     if(_cubicFaceIter != NULL)
       delete _cubicFaceIter;
   }

  Self& operator=(const BeginType& grid)
  {
    if(_cubicFaceIter != NULL)
      delete _cubicFaceIter;

    _cubicFaceIter = NULL;
    _cubicFaceIter = new CubicBoundaryFaceElementIterator(grid);
    _simplexFaceIter = _cubicFaceIter->operator*();
    return *this;
  }

  bool operator==(const EndType&) const
  {
    if(_cubicFaceIter == NULL)
    {
      throw aol::Exception("SimplexBoundaryFaceElementIterator::operator==: _cubicFaceIter == NULL but needed.", __FILE__, __LINE__);
      return true;               //at end
    }

    bool cubicAtEnd = _cubicFaceIter->atEnd();
    bool simplexAtEnd = _simplexFaceIter.atEnd();
    return ( cubicAtEnd && simplexAtEnd );
  }

  bool notAtEnd() const
  {
    if(_cubicFaceIter == NULL)
    {
      throw aol::Exception("SimplexBoundaryFaceElementIterator::notAtEnd==: _cubicFaceIter == NULL but needed.", __FILE__, __LINE__);
      return false;             //at end
    }

    bool cubicAtEnd = _cubicFaceIter->atEnd();
    bool simplexAtEnd = _simplexFaceIter.atEnd();
    return !( cubicAtEnd && simplexAtEnd );
  }

   Self& operator++()
   {
     if(_cubicFaceIter == NULL)
     {
       throw aol::Exception("SimplexBoundaryFaceElementIterator::operator++: _cubicFaceIter == NULL but needed.", __FILE__, __LINE__);
       return *this;          //do nothing
     }

     if(_simplexFaceIter.atEnd())
     {
       ++(*_cubicFaceIter);

       if(!_cubicFaceIter->atEnd())
       {
         _simplexFaceIter = _cubicFaceIter->operator*();
       }
     }

     else
     {
       ++_simplexFaceIter;
     }

     return *this;
  }

  Self operator++(int)
  {
    Self tmp(*this);
    ++(*this);

    return tmp;
  }

  IteratedType& operator*()
  {
    return _simplexFaceIter.operator*();
  }

  IteratedType* operator->()
  {
    return _simplexFaceIter.operator->();
  }

protected:
  CubicBoundaryFaceElementIterator* _cubicFaceIter;
  SimplexInFaceIterator<CubicBoundaryFaceElementIterator, _GridType, Dim, RealType> _simplexFaceIter;
};

//! node iterator for simplicial grid
//! (which actually only iterates over all cubic grid elements, including virtual, non-existent elements at boundary nodes)
//! \author wirth
template<typename CubicGridType, Dimension Dim>
class FullNodeIterator {
  typedef FullNodeIterator Self;

protected:
  typedef Element<typename CubicGridType::ElementType> ElementType;
  ElementType _cur;
  qc::GridSize<Dim> _size;

public:
  FullNodeIterator() :
    _size( aol::Vec<Dim,short>() ) {}

  explicit FullNodeIterator( const qc::GridSize<Dim> &Size ) :
    _size ( Size ) {}

  FullNodeIterator( const qc::GridSize<Dim> &Size, const ElementType &El ) :
    _size ( Size ),
    _cur( El ) {}

  FullNodeIterator &operator= ( const FullNodeIterator& Other ) {
    _size = Other._size;
    _cur = Other._cur;
    return *this;
  }

  inline FullNodeIterator& operator++() {
    _cur.getCubicElement().xref() ++;
    if ( _cur.getCubicElement().x() == ( _size.getNumX() ) ) {
      _cur.getCubicElement().xref() = 0;
      _cur.getCubicElement().yref() ++;
      if ( _cur.getCubicElement().y() == ( _size.getNumY() ) ) {
        _cur.getCubicElement().yref() = 0;
        _cur.getCubicElement().zref() ++;
      }
    }
    return *this;
  }

  ElementType& operator*() {
    return _cur;
  }

  const ElementType& operator*() const {
    return _cur;
  }

  const ElementType* operator->() const {
    return &_cur;
  }

  inline bool operator!= ( const Self& Other ) const {
    return ( _cur != Other._cur );
  }

  inline bool operator== ( const Self& Other ) const {
    return ( _cur == Other._cur );
  }
};


// ------------------------------------------------------------------------------------------------

//! simplicial grid with iterators
//! \author von Deylen (july 2008)
template <typename _CubicGridType, Dimension Dim>
class GridStructure : public GridStructureBase<_CubicGridType> {
public:
  typedef FullSimplexIterator<typename _CubicGridType::OldFullElementIterator,
                               _CubicGridType, Dim> OldFullElementIterator;
  typedef typename _CubicGridType::OldFullNodeIterator OldFullNodeIterator;
  typedef FullNodeIterator<_CubicGridType,Dim> OldAllNodeIterator;
  typedef Element<typename _CubicGridType::ElementType> ElementType;
  typedef _CubicGridType                       CubicGridType;
  typedef typename _CubicGridType::ElementType CubicElementType;
  typedef _CubicGridType                       BeginIterType;
  typedef EndElement                           EndIterType;

  OldAllNodeIterator _nBeginIt;
  OldAllNodeIterator _nEndIt;

  static const Dimension DimOfWorld = Dim;

  explicit GridStructure ( const _CubicGridType & CubicGrid )
  : GridStructureBase<_CubicGridType> ( CubicGrid ),
    _nBeginIt( qc::GridSize<DimOfWorld>::createFrom( CubicGrid ) ),
    _nEndIt( qc::GridSize<DimOfWorld>::createFrom( CubicGrid ) ) {
    if ( DimOfWorld == qc::QC_1D )
      (*_nEndIt).set( CubicElementType( 0, 1, 0, 0 ), 0 );
    else if ( DimOfWorld == qc::QC_2D )
      (*_nEndIt).set( CubicElementType( 0, 0, 1, 0 ), 0 );
    else if ( DimOfWorld == qc::QC_3D )
      (*_nEndIt).set( CubicElementType( 0, 0, CubicGrid.getNumZ(), 0 ), 0 );
    else
      throw aol::Exception ( "unsupported dimension" );
  }

  template <qc::Dimension DimOfWorld>
  explicit GridStructure ( const GridSize<DimOfWorld> &GridSize )
  : GridStructureBase<_CubicGridType> ( GridSize ),
    _nBeginIt( GridSize ),
    _nEndIt( GridSize ) {
    if ( DimOfWorld == qc::QC_1D )
      (*_nEndIt).set( CubicElementType( 0, 1, 0, 0 ), 0 );
    else if ( DimOfWorld == qc::QC_2D )
      (*_nEndIt).set( CubicElementType( 0, 0, 1, 0 ), 0 );
    else if ( DimOfWorld == qc::QC_3D )
      (*_nEndIt).set( CubicElementType( 0, 0, GridSize.getNumZ(), 0 ), 0 );
    else
      throw aol::Exception ( "unsupported dimension" );
  }

  const BeginIterType & begin () const {
    return this->getCubicGrid();
  }

  const EndIterType & end () const {
    return this->getCubicGrid().end();
  }

  int getGridDepth() const {
    return this->getCubicGrid().getGridDepth();
  }

  int getNumberOfElements () const {
    short numSimpl = TopologyLookup<Dim>::numSimplexesPerCube;
    return this->getCubicGrid().getNumberOfElements() * numSimpl;
  }

  const aol::Vec3<int>& getSize () const {
    return this->getCubicGrid().getSize();
  }

  inline Dimension getDimOfWorld() const {
    return DimOfWorld;
  }

  int getElementIndex ( const ElementType & El ) const {
    return this->getCubicGrid().getElementIndex ( El.getCubicElement() )
                                  * TopologyLookup<Dim>::numSimplexesPerCube
      + El.getSimplexNumber();
  }

  const GridStructure<_CubicGridType, Dim> & getFullGrid () const {
    return *this;
  }
};

// ------------------------------------------------------------------------------------------------

} // end of namespace simplex.

} // end of namespace qc.

#endif
