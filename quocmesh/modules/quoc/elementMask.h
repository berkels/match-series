#ifndef __ELEMENTMASK_H
#define __ELEMENTMASK_H

#include <gridBase.h>
#include <simplexGrid.h>
#include <quoc.h>

namespace qc {

/*! \brief Masking class that divides a cubicGrid into two regions.
 *  \author toelkes
 *  \warning Only regions with rectangular inner region are implemented.
 *  \warning Only implemented for non-adaptive grids.
 *
 *  This class is designed to divide the unit cube into an inner and outer region of arbitrary (reasonable) sizes and state, if a given element is in the inner or outer region.
 *  Element numbering for inner elements is provided.
 */
template<typename GridType, typename ConfiguratorType, qc::Dimension Dim>
class ElementMask {
public:
  template < typename ElementType >
  bool isInDesiredRegion ( const ElementType &/*el*/ ) const {
    throw aol::Exception ( "Element mask not implemented for this element type", __FILE__, __LINE__, __FUNCTION__ );
  }
};

/*! \brief 2D implementation of qc::ElementMask
 *  \author toelkes
 */
template <typename ConfiguratorType>
class ElementMask<qc::GridDefinition, ConfiguratorType, qc::QC_2D> {
  typedef qc::GridDefinition GridType;
public:
  enum RegionType { INNER_REGION = 0,
                    OUTER_REGION = 1
                  };

  typedef typename ConfiguratorType::ElementType ElementType;

  //! \brief Constructor for centered quadratic inner region.
  //! \param a width of the outer region (in elements)
  //! \author toelkes
  ElementMask ( const GridType& grid, const ConfiguratorType& configurator, const unsigned int a )
    : _grid ( grid ), _configurator ( configurator ), _a ( a ), _b ( a ), _innerIsDesired ( true ) {

    //The number of elements per cell is the number of elements divided by the number of cell. Do a consistency check first..
    if ( ( _grid.getNumberOfElements() % ( ( _grid.getNumX() - 1 ) * ( _grid.getNumY() - 1 ) ) ) != 0 )
      throw aol::Exception ( "ElementMask: Could not calculate the number of elements per cell!", __FILE__, __LINE__ );

    //...and calculate.
    _numElementsPerCell = ( _grid.getNumberOfElements() / ( ( _grid.getNumX() - 1 ) * ( _grid.getNumY() - 1 ) ) );
  }

  //! \brief Return true if el is in the inner region.
  //! \author toelkes
  bool isInInnerRegion ( const ElementType &el ) const {
    unsigned int elNum = _configurator.getConsecutiveElementNumber ( el );
    unsigned int row = elNum / ( ( _grid.getNumX() - 1 ) * _numElementsPerCell );       //Calculate the row of el.

    if ( ( row < _a ) || ( row >= ( _grid.getNumY() - _a - 1 ) ) )                      //If the element is below or above the inner region, return false.
      return false;

    unsigned int col = elNum % ( ( _grid.getNumX() - 1 ) * _numElementsPerCell );       //Calculate the column of el.
    col /= _numElementsPerCell;

    if ( ( col < _b ) || ( col >= ( _grid.getNumX() - _b - 1 ) ) )                      //If the element is right or left of the inner region, return false.
      return false;

    return true;                                                                //If none of the above is true, the element has to be in the inner region. Return true.
  }

  //! \brief Return true if el is in the desired region (determined by _innerIsDesired).
  //! \author toelkes
  bool isInDesiredRegion ( const ElementType &el ) const {
    if ( _innerIsDesired )
      return isInInnerRegion ( el );
    else
      return ! ( isInInnerRegion ( el ) );
  }

  //! \brief Return true if el is in the desired region (determined by _innerIsDesired). Not implemented for general elements.
  template < typename ElementType >
  bool isInDesiredRegion ( const ElementType &el ) const {
    throw aol::Exception ( "Element mask not implemented for this element type", __FILE__, __LINE__, __FUNCTION__ );
  }

  //! \brief Set desired region.
  //! \author toelkes
  void setDesiredRegion ( RegionType rt ) {
    if ( rt == INNER_REGION )
      _innerIsDesired = true;
    else if ( rt == OUTER_REGION )
      _innerIsDesired = false;
    else
      throw aol::Exception ( "ElementMask: Undefined region!", __FILE__, __LINE__ );
  }

  //! \brief Swaps desired region
  //! \author toelkes
  void swapDesiredRegion() {
    _innerIsDesired = !_innerIsDesired;
  }

  //! \brief Consecutively numbers elements in one (the inner) region.
  //! \warning Currently only inner elements have regional element numbers
  //! \author toelkes
  unsigned int getRegionalElementNumber ( const ElementType &el ) const {
    if ( !_innerIsDesired ) {
      throw aol::Exception ( "ElementMask: Outer region element numbering not implemented.", __FILE__, __LINE__ );
      return 0;
    }

    unsigned int elNum = _configurator.getConsecutiveElementNumber ( el );
    unsigned int row = elNum / ( ( _grid.getNumX() - 1 ) * _numElementsPerCell );       //Calculate the row of el.
    unsigned int col = elNum % ( ( _grid.getNumX() - 1 ) * _numElementsPerCell );       //Calculate the column of el.
    unsigned int simplexNumber = col % _numElementsPerCell;                             //Calculate if el is the first, second, ... element of the cell
    col /= _numElementsPerCell;                                                         //finish calculation of col

    //Return consecutive "inner element number"
    return ( ( row - _a ) * ( _grid.getNumX() - 2 * _b - 1 ) + ( col - _b ) ) * _numElementsPerCell + simplexNumber;
  }

  //! \brief Returns the number of inner elements
  //! \author toelkes
  unsigned int getNumInnerElements() {
    return _numElementsPerCell * ( ( _grid.getNumX() - 2 * _b - 1 ) * ( _grid.getNumY() - 2 * _a - 1 ) );
  }

protected:
  // Do not use standart constructor!
  ElementMask();

  const GridType& _grid;
  const ConfiguratorType& _configurator;
  //! _a is the distance of the inner region (in elements) from the upper and lower bound of the rectangular grid, _b is the distance from the right and left bound.
  unsigned int _a, _b;
  unsigned int _numElementsPerCell;
  bool _innerIsDesired;
};

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief 3D implementation of qc::ElementMask
 *  \author toelkes
 */
template <typename ConfiguratorType>
class ElementMask<qc::GridDefinition, ConfiguratorType, qc::QC_3D> {
  typedef qc::GridDefinition GridType;
public:
  enum RegionType { INNER_REGION = 0,
                    OUTER_REGION = 1
                  };

  typedef typename ConfiguratorType::ElementType ElementType;

  //! \brief Constructor for centered cubic inner region.
  //! \param a width of the outer region (in elements)
  //! \author toelkes
  ElementMask ( const GridType& grid, const ConfiguratorType& configurator, const unsigned int a )
    : _grid ( grid ), _configurator ( configurator ), _a ( a ), _b ( a ), _c ( a ), _innerIsDesired ( true ) {

    //The number of elements per cell is the number of elements divided by the number of cell. Do a consistency check first..
    if ( ( _grid.getNumberOfElements() % ( ( _grid.getNumX() - 1 ) * ( _grid.getNumY() - 1 ) * ( _grid.getNumZ() - 1 ) ) ) != 0 )
      throw aol::Exception ( "ElementMask: Could not calculate the number of elements per cell!", __FILE__, __LINE__ );

    //...and calculate.
    _numElementsPerCell = ( _grid.getNumberOfElements() / ( ( _grid.getNumX() - 1 ) * ( _grid.getNumY() - 1 ) * ( _grid.getNumZ() - 1 ) ) );
  }

  //! \brief Return true if el is in the inner region.
  //! \author toelkes
  bool isInInnerRegion ( const ElementType &el ) const {
    unsigned int elNum = _configurator.getConsecutiveElementNumber ( el );
    unsigned int row = ( elNum / ( ( _grid.getNumX() - 1 ) * _numElementsPerCell ) ) % ( _grid.getNumY() - 1 );     //Calculate the row of el.

    if ( ( row < _a ) || ( row >= ( _grid.getNumY() - _a - 1 ) ) )                                                  //If the element is below or above the inner region, return false.
      return false;

    unsigned int col = elNum % ( ( _grid.getNumX() - 1 ) * _numElementsPerCell );                                   //Calculate the column of el.
    col /= _numElementsPerCell;

    if ( ( col < _b ) || ( col >= ( _grid.getNumX() - _b - 1 ) ) )                                                  //If the element is right or left of the inner region, return false.
      return false;

    //Calculate the layer (=row in z direction) of el.
    unsigned int layer = ( elNum / ( ( _grid.getNumX() - 1 ) * ( _grid.getNumY() - 1 ) * _numElementsPerCell ) );

    if ( ( layer < _c ) || ( layer >= ( _grid.getNumZ() - _c - 1 ) ) )                                               //If the element is in front of or behind the inner region, return false.
      return false;

    return true;               //If none of the above is true, the element has to be in the inner region. Return true.
  }

  //! \brief Return true if el is in the desired region (determined by _innerIsDesired).
  //! \author toelkes
  bool isInDesiredRegion ( const ElementType &el ) const {
    if ( _innerIsDesired )
      return isInInnerRegion ( el );
    else
      return ! ( isInInnerRegion ( el ) );
  }

  //! \brief Return true if el is in the desired region (determined by _innerIsDesired). Not implemented for general elements.
  template < typename ElementType >
  bool isInDesiredRegion ( const ElementType &el ) const {
    throw aol::Exception ( "Element mask not implemented for this element type", __FILE__, __LINE__, __FUNCTION__ );
  }

  //! \brief Set desired region.
  //! \author toelkes
  void setDesiredRegion ( RegionType rt ) {
    if ( rt == INNER_REGION )
      _innerIsDesired = true;
    else if ( rt == OUTER_REGION )
      _innerIsDesired = false;
    else
      throw aol::Exception ( "ElementMask: Undefined region!", __FILE__, __LINE__ );
  }

  //! \brief Swaps desired region
  //! \author toelkes
  void swapDesiredRegion() {
    _innerIsDesired = !_innerIsDesired;
  }

  //! \brief Consecutively numbers elements in one (the inner) region.
  //! \warning Currently only inner elements have regional element numbers
  //! \author toelkes
  unsigned int getRegionalElementNumber ( const ElementType &el ) const {
    if ( !_innerIsDesired ) {
      throw aol::Exception ( "ElementMask: Outer region element numbering not implemented.", __FILE__, __LINE__ );
      return 0;
    }

    unsigned int elNum = _configurator.getConsecutiveElementNumber ( el );
    unsigned int row = ( elNum / ( ( _grid.getNumX() - 1 ) * _numElementsPerCell ) ) % ( _grid.getNumY() - 1 );     //Calculate the row of el.
    unsigned int col = elNum % ( ( _grid.getNumX() - 1 ) * _numElementsPerCell );                                   //Calculate the column of el.
    unsigned int simplexNumber = col % _numElementsPerCell;                                                         //Calculate if el is the first, second, ... element of the cell
    col /= _numElementsPerCell;                                                                                     //finish calculation of col

    //Calculate the layer (=row in z direction) of el.
    unsigned int layer = ( elNum / ( ( _grid.getNumX() - 1 ) * ( _grid.getNumY() - 1 ) * _numElementsPerCell ) );

    //Return consecutive "inner element number"
    return ( ( layer - _c ) * ( ( _grid.getNumX() - 2 * _b - 1 ) * ( _grid.getNumY() - 2 * _a - 1 ) ) + ( row - _a ) * ( _grid.getNumX() - 2 * _b - 1 ) + ( col - _b ) ) * _numElementsPerCell + simplexNumber;
  }

  //! \brief Returns the number of inner elements
  //! \author toelkes
  unsigned int getNumInnerElements() {
    return _numElementsPerCell * ( ( _grid.getNumX() - 2 * _b - 1 ) * ( _grid.getNumY() - 2 * _a - 1 ) * ( _grid.getNumZ() - 2 * _c - 1 ) );
  }

protected:
  // Do not use standart constructor!
  ElementMask();

  const GridType& _grid;
  const ConfiguratorType& _configurator;
  //! _a is the distance of the inner region (in elements) from the upper and lower bound of the rectangular grid, _b is the distance from the right and left bound, _c is the distance from the front and back bound.
  unsigned int _a, _b, _c;
  unsigned int _numElementsPerCell;
  bool _innerIsDesired;
};

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief 3D implementation of qc::ElementMask for simplicial grids
 *  \warning code duplication from above
 */
template <typename ConfiguratorType>
class ElementMask<qc::simplex::GridStructure<qc::GridDefinition,qc::QC_3D>, ConfiguratorType, qc::QC_3D> {
  typedef qc::simplex::GridStructure<qc::GridDefinition,qc::QC_3D> GridType;
public:
  enum RegionType { INNER_REGION = 0,
    OUTER_REGION = 1
  };

  typedef typename ConfiguratorType::ElementType ElementType;

  //! \brief Constructor for centered cubic inner region.
  //! \param a width of the outer region (in elements)
  //! \author toelkes
  ElementMask ( const GridType& grid, const ConfiguratorType& configurator, const unsigned int a )
  : _grid ( grid ), _configurator ( configurator ), _a ( a ), _b ( a ), _c ( a ), _innerIsDesired ( true ) {

    //The number of elements per cell is the number of elements divided by the number of cell. Do a consistency check first..
    if ( ( _grid.getNumberOfElements() % ( ( _grid.getNumX() - 1 ) * ( _grid.getNumY() - 1 ) * ( _grid.getNumZ() - 1 ) ) ) != 0 )
      throw aol::Exception ( "ElementMask: Could not calculate the number of elements per cell!", __FILE__, __LINE__ );

    //...and calculate.
    _numElementsPerCell = ( _grid.getNumberOfElements() / ( ( _grid.getNumX() - 1 ) * ( _grid.getNumY() - 1 ) * ( _grid.getNumZ() - 1 ) ) );
  }

  //! \brief Return true if el is in the inner region.
  //! \author toelkes
  bool isInInnerRegion ( const ElementType &el ) const {
    unsigned int elNum = _configurator.getConsecutiveElementNumber ( el );
    unsigned int row = ( elNum / ( ( _grid.getNumX() - 1 ) * _numElementsPerCell ) ) % ( _grid.getNumY() - 1 );     //Calculate the row of el.

    if ( ( row < _a ) || ( row >= ( _grid.getNumY() - _a - 1 ) ) )                                                  //If the element is below or above the inner region, return false.
      return false;

    unsigned int col = elNum % ( ( _grid.getNumX() - 1 ) * _numElementsPerCell );                                   //Calculate the column of el.
    col /= _numElementsPerCell;

    if ( ( col < _b ) || ( col >= ( _grid.getNumX() - _b - 1 ) ) )                                                  //If the element is right or left of the inner region, return false.
      return false;

    //Calculate the layer (=row in z direction) of el.
    unsigned int layer = ( elNum / ( ( _grid.getNumX() - 1 ) * ( _grid.getNumY() - 1 ) * _numElementsPerCell ) );

    if ( ( layer < _c ) || ( layer >= ( _grid.getNumZ() - _c - 1 ) ) )                                               //If the element is in front of or behind the inner region, return false.
      return false;

    return true;               //If none of the above is true, the element has to be in the inner region. Return true.
  }

  //! \brief Return true if el is in the desired region (determined by _innerIsDesired).
  //! \author toelkes
  bool isInDesiredRegion ( const ElementType &el ) const {
    if ( _innerIsDesired )
      return isInInnerRegion ( el );
    else
      return ! ( isInInnerRegion ( el ) );
  }

  //! \brief Return true if el is in the desired region (determined by _innerIsDesired). Not implemented for general elements.
  template < typename ElementType >
  bool isInDesiredRegion ( const ElementType &el ) const {
    throw aol::Exception ( "Element mask not implemented for this element type", __FILE__, __LINE__, __FUNCTION__ );
  }

  //! \brief Set desired region.
  //! \author toelkes
  void setDesiredRegion ( RegionType rt ) {
    if ( rt == INNER_REGION )
      _innerIsDesired = true;
    else if ( rt == OUTER_REGION )
      _innerIsDesired = false;
    else
      throw aol::Exception ( "ElementMask: Undefined region!", __FILE__, __LINE__ );
  }

  //! \brief Swaps desired region
  //! \author toelkes
  void swapDesiredRegion() {
    _innerIsDesired = !_innerIsDesired;
  }

  //! \brief Consecutively numbers elements in one (the inner) region.
  //! \warning Currently only inner elements have regional element numbers
  //! \author toelkes
  unsigned int getRegionalElementNumber ( const ElementType &el ) const {
    if ( !_innerIsDesired ) {
      throw aol::Exception ( "ElementMask: Outer region element numbering not implemented.", __FILE__, __LINE__ );
      return 0;
    }

    unsigned int elNum = _configurator.getConsecutiveElementNumber ( el );
    unsigned int row = ( elNum / ( ( _grid.getNumX() - 1 ) * _numElementsPerCell ) ) % ( _grid.getNumY() - 1 );     //Calculate the row of el.
    unsigned int col = elNum % ( ( _grid.getNumX() - 1 ) * _numElementsPerCell );                                   //Calculate the column of el.
    unsigned int simplexNumber = col % _numElementsPerCell;                                                         //Calculate if el is the first, second, ... element of the cell
    col /= _numElementsPerCell;                                                                                     //finish calculation of col

    //Calculate the layer (=row in z direction) of el.
    unsigned int layer = ( elNum / ( ( _grid.getNumX() - 1 ) * ( _grid.getNumY() - 1 ) * _numElementsPerCell ) );

    //Return consecutive "inner element number"
    return ( ( layer - _c ) * ( ( _grid.getNumX() - 2 * _b - 1 ) * ( _grid.getNumY() - 2 * _a - 1 ) ) + ( row - _a ) * ( _grid.getNumX() - 2 * _b - 1 ) + ( col - _b ) ) * _numElementsPerCell + simplexNumber;
  }

  //! \brief Returns the number of inner elements
  //! \author toelkes
  unsigned int getNumInnerElements() {
    return _numElementsPerCell * ( ( _grid.getNumX() - 2 * _b - 1 ) * ( _grid.getNumY() - 2 * _a - 1 ) * ( _grid.getNumZ() - 2 * _c - 1 ) );
  }

protected:
  // Do not use standart constructor!
  ElementMask();

  const GridType& _grid;
  const ConfiguratorType& _configurator;
  //! _a is the distance of the inner region (in elements) from the upper and lower bound of the rectangular grid, _b is the distance from the right and left bound, _c is the distance from the front and back bound.
  unsigned int _a, _b, _c;
  unsigned int _numElementsPerCell;
  bool _innerIsDesired;
};

}

#endif
