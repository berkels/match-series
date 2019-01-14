#ifndef __ADAPTIVEGRIDTP_H
#define __ADAPTIVEGRIDTP_H

#include <quoc.h>
#include <gridBase.h>
#include <op.h>
#include <indexMapper.h>
#include <scalarArray.h>

/**
* Class AdaptiveGrid is the representation of adaptive grids in quocmesh
*
*  \author Paetz
*/
template <typename EstType, typename RealType>
class AdaptiveGridTP : public qc::GridDefinition {
public:
  class ElementIterator {
  public:
    typedef ElementIterator _Self;
    typedef _Self    BeginType;
    typedef qc::EndElement    EndType;

    ElementIterator() : _grid ( NULL ), _estimator ( NULL ), _grid_depth ( -1 ), _more ( true ) {
      _it_stack.reserve ( 2000 );
      _it_stack.push_back ( qc::Element ( 0, 0, 0, 0 ) );
    }

    ElementIterator ( const AdaptiveGridTP<EstType, RealType> &Grid ) : _estimator ( Grid.getEstimator() ), _grid ( &Grid ), _more ( true ) {
      _grid_depth = Grid.getGridDepth();
      _estimator = Grid.getEstimator();
      _it_stack.reserve ( 2000 );
      _it_stack.push_back ( qc::Element ( 0, 0, 0, 0 ) );
      _more = true;
      operator++ ( 1 );
    }

    qc::Element& operator*() {
      return _cur;
    }

    _Self& operator= ( const AdaptiveGridTP<EstType, RealType> &Grid ) {
      _grid = &Grid;
      _grid_depth = _grid->getGridDepth();
      _estimator = _grid->getEstimator();
      _more = true;
      _it_stack.erase ( _it_stack.begin(), _it_stack.end() );
      _it_stack.push_back ( qc::Element ( 0, 0, 0, 0 ) );
      operator++ ( 1 );
      return *this;
    }

    void restart() {
      _more = true;
      _it_stack.erase ( _it_stack.begin(), _it_stack.end() );
      _it_stack.push_back ( qc::Element ( 0, 0, 0, 0 ) );
      operator++ ( 1 );
    }


    /**
     * Compare current iterator element with another iterator.
     */
    bool operator== ( const _Self &__x ) const {
      return _cur == __x._cur;
    }

    /**
     * Compare current iterator element with another iterator.
     */
    bool operator!= ( const _Self &__x ) const {
      return _cur != __x._cur;
    }
    /**
     * check, whether iterator has reached the end
     */
    bool operator== ( const qc::EndElement & ) const {
      return !_it_stack.size();
    }

    /**
     * check, whether iterator has not reached the end
     */
    bool operator!= ( const qc::EndElement & ) const {
      return _more;
    } //{ return _it_stack.size(); }

    qc::Element* operator->() {
      return &_cur;
    }

    inline ElementIterator& operator++ ( int );

    inline ElementIterator& operator++ ( ) {
      return this->operator++ ( 1 );
    }

    EstType *_estimator;

    void dump() {
      typename vector<qc::Element>::iterator it;
      for ( it = _it_stack.begin(); it != _it_stack.end(); ++it ) {
        cout << *it << endl;
      }
    }
  private:

    vector<qc::Element> _it_stack;
    qc::Element _cur;
    const AdaptiveGridTP<EstType, RealType> *_grid;
    bool _more;
    int _grid_depth;
  };

  class NodeIterator {
  public:

    typedef NodeIterator _Self;

    NodeIterator() : _grid ( NULL ), _array ( NULL ), localNodeCounter ( -1 ) {}

    ~NodeIterator() {
      if ( _array ) delete _array;
    }

    inline _Self& operator= ( AdaptiveGridTP<EstType, RealType> &Grid );

    inline _Self& operator++ ( int ) {

      while ( 1 ) {
        if ( ++localNodeCounter ) {}
      };
      return *this;
    }

  protected:
    ElementIterator it;
    qc::Array<unsigned char> *_array;
    aol::Vec3<int> localNodes[ 8 ];
    int localNodeCounter;
    AdaptiveGridTP<EstType, RealType> *_grid;
  };

protected:
  int *h;           /* spatial gridstep IN IMAGES, i.e. = 1 on finest grid*/
  float *hSqr;      /* square of spatial gridstep */
  float *hSqr_2;    /* square of spatial gridstep divided by 2 */
  float *spaceH;
  //lookuptables for local to global dof mapping
  int **nmap,       /* maps the nodes of one element from the lower left node */
      *map,         /* maps a y coordinate to dof */
      **cmap,       /* maps the nodes of children of an el from the lower left node */
      ***hmap;      /* maps the constraining nodes of one hanging node */


  /** Constructor for adaptive grids
   *  @param Estimator estimator that defines the adaptive 3D grid
   */
public:

  ElementIterator& makeElementIterator() const {
    elemIt->restart();
    return *elemIt;
  }

  ElementIterator *elemIt;
  const qc::EndElement _endEl;

  explicit AdaptiveGridTP ( EstType *Estimator ) :
    GridDefinition ( Estimator->getMaxLevel(),
                     Estimator->getDimOfWorld() ), _endEl ( qc::EndElement() ), _estimator ( Estimator ) {
    nodeIsVertexOfElement = new aol::BitVector ( width * width );
    nodeIsHangingNode = new aol::BitVector ( width * width );

    init();
  }

  explicit AdaptiveGridTP ( const aol::Vec3<int> &size ) : GridDefinition ( qc::logBaseTwo ( size.x() ),
                                                                            qc::QC_2D ), _estimator ( NULL ), _endEl ( qc::EndElement() ) {
    if ( size.z() != 1 ) {
      cerr << "Adaptive grid not implemented for 3D!" << endl;
      exit ( 1 );
    }
    nodeIsVertexOfElement = new aol::BitVector ( width * width );
    nodeIsHangingNode = new aol::BitVector ( width * width );
    init();
  }


  explicit AdaptiveGridTP ( const qc::GridSize<qc::QC_2D> &gridSize ) : GridDefinition ( qc::logBaseTwo ( gridSize.getNumX() ),
                                                                                         qc::QC_2D ), _estimator ( NULL ) {
    if ( gridSize.getNumZ() != 1 ) {
      cerr << "Adaptive grid not implemented for 3D!" << endl;
      exit ( 1 );
    }
    nodeIsVertexOfElement = new aol::BitVector ( width * width );
    nodeIsHangingNode = new aol::BitVector ( width * width );
    init();
  }

  void init() {
    int i, j, w;

    nmap = new int*[gridDepth + 1];
    cmap = new int*[gridDepth + 1];
    h = new int[gridDepth + 1];
    hSqr = new float[gridDepth + 1];
    spaceH = new float[gridDepth + 1];
    hSqr_2 = new float[gridDepth + 1];
    map = new int[width];
    hmap = new int**[gridDepth + 1];

    w = width - 1;
    spaceH[0] = w;

    for ( i = 0; i <= gridDepth; i++ ) {
      nmap[i] = new int[4];
      cmap[i] = new int[5];
      hmap[i] = new int*[4];
      for ( j = 0; j < 4; j++ ) {
        hmap[i][j] = new int[8];
      }

      nmap[i][0] = 0;
      nmap[i][1] = w;
      nmap[i][2] = w + width * w;
      nmap[i][3] = width * w;

      cmap[i][0] = w >> 1;
      cmap[i][1] = w + width * ( w >> 1 );
      cmap[i][2] = ( w >> 1 ) + width * ( w >> 1 );
      cmap[i][3] = width * ( w >> 1 );
      cmap[i][4] = ( w >> 1 ) + width * w;

      h[i] = w;
      if ( i > 0 ) spaceH[i] = spaceH[i - 1] / 2.;
      hSqr[i] = spaceH[i] * spaceH[i];
      hSqr_2[i] = hSqr[i] / 2.;

      if ( i > 0 ) {
        hmap[i][0][0] = hmap[i][0][1] = 0;
        hmap[i][0][2] = 0;
        hmap[i][0][3] = nmap[i - 1][1];
        hmap[i][0][4] = hmap[i][0][5] = 0;
        hmap[i][0][6] = 0;
        hmap[i][0][7] = nmap[i - 1][3];

        hmap[i][1][0] = -nmap[i][1];
        hmap[i][1][1] = nmap[i][1];
        hmap[i][1][2] = hmap[i][1][3] = 0;
        hmap[i][1][4] = nmap[i][1];
        hmap[i][1][5] = nmap[i - 1][2] - nmap[i][1];
        hmap[i][1][6] = hmap[i][1][7] = 0;

        hmap[i][2][0] = hmap[i][2][1] = 0;
        hmap[i][2][2] = nmap[i][2];
        hmap[i][2][3] = nmap[i][1] - nmap[i][3];
        hmap[i][2][4] = hmap[i][2][5] = 0;
        hmap[i][2][6] = nmap[i][2];
        hmap[i][2][7] = -nmap[i][1] + nmap[i][3];

        hmap[i][3][0] = -nmap[i][3];
        hmap[i][3][1] = nmap[i][3];
        hmap[i][3][2] = hmap[i][3][3] = 0;
        hmap[i][3][4] = nmap[i][1] + nmap[i][2];
        hmap[i][3][5] = nmap[i][3];
        hmap[i][3][6] = hmap[i][3][7] = 0;
      }

      w >>= 1;
    }

    for ( i = 0; i < width; i++ ) map[i] = width * i;

    elemIt = new ElementIterator ( *this );

  }




  //! Checks whether the passed Element is valid.
  bool elementInGrid ( const qc::Element &El ) const {
    if ( !isAdaptive() ) {
      if ( !nodeInGrid ( El, El.level() ) ) return false;
    }
    int elSize = 1 << ( gridDepth - El.level() );

    switch ( this->getDimOfWorld() ) {
      case qc::QC_1D:
        return ( El.x() >= 0 &&
                 El.x() + elSize < this->width );
        break;
      case qc::QC_2D:
        if ( !isAdaptive() ) {
          return ( El.x() >= 0 &&
                   El.x() + elSize < this->width &&
                   El.y() >= 0 &&
                   El.y() + elSize < this->width );
        } else {
          int ind1 = ( El.x() + ( ( 0 & 1 ) != 0 ) * elSize ) * 1 + ( El.y() + ( ( 0 & 2 ) != 0 ) * elSize ) * this->width;
          int ind2 = ( El.x() + ( ( 1 & 1 ) != 0 ) * elSize ) * 1 + ( El.y() + ( ( 1 & 2 ) != 0 ) * elSize ) * this->width;
          int ind3 = ( El.x() + ( ( 2 & 1 ) != 0 ) * elSize ) * 1 + ( El.y() + ( ( 2 & 2 ) != 0 ) * elSize ) * this->width;
          int ind4 = ( El.x() + ( ( 3 & 1 ) != 0 ) * elSize ) * 1 + ( El.y() + ( ( 3 & 2 ) != 0 ) * elSize ) * this->width;
          if ( nodeIsVertexOfElement->get ( ind1 ) && nodeIsVertexOfElement->get ( ind2 )
               && nodeIsVertexOfElement->get ( ind3 ) && nodeIsVertexOfElement->get ( ind4 ) ) {
            return true;
          } else {
            return false;
          }
        }
        break;
      case qc::QC_3D:
        return ( El.x() >= 0 &&
                 El.x() + elSize < this->width &&
                 El.y() >= 0 &&
                 El.y() + elSize < this->width &&
                 El.z() >= 0 &&
                 El.z() + elSize < this->width  );
        break;
    }
    return 0;
  }





  /* interpolateResult interpolates the result vector on all inactive
  * nodes
  */
  void interpolateResult ( aol::Vector<double> &result ) const {
    int x, y, d;

    for ( d = 0; d < gridDepth; d++ ) {
      for ( x = 0; x < this->getNumX() - 1; x += h[d] ) {
        for ( y = 0; y < this->getNumX() - 1; y += h[d] ) {
          int k = x + map[y];

          int level = d;

          if ( level == gridDepth ) return;

          if ( !isActive ( k + cmap[level][0] ) )
            ( result ) [k + cmap[level][0]] = ( ( result ) [k] +
                                                ( result ) [k + nmap[level][1]] ) / 2.;

          if ( !isActive ( k + cmap[level][3] ) )
            ( result ) [k + cmap[level][3]] = ( ( result ) [k] +
                                                ( result ) [k + nmap[level][3]] ) / 2.;

          if ( !isActive ( k + cmap[level][1] ) )
            ( result ) [k + cmap[level][1]] =
              ( ( result ) [k + nmap[level][1]] + ( result ) [k + nmap[level][2]] ) / 2.;

          if ( !isActive ( k + cmap[level][4] ) )
            ( result ) [k + cmap[level][4]] =
              ( ( result ) [k + nmap[level][3]] + ( result ) [k + nmap[level][2]] ) / 2.;

          if ( !isActive ( k + cmap[level][2] ) )
            ( result ) [k + cmap[level][2]] =
              ( ( result ) [k + cmap[level][3]] + ( result ) [k + cmap[level][1]] ) / 2.;
        }
      }
    }
  }


  void adoptToAdaptiveGrid ( aol::Vector<double> &result ) const {
    qc::Element cEl;
    int numY = this->getNumY();

    for ( int d = gridDepth; d < 0; d-- ) {
      const int elSize = 1 << ( gridDepth - d );
      const int elSizeBig = 1 << ( gridDepth - ( d - 1 ) );
      const int maxNodes = ( 1 << gridDepth ) /*- 1*/ - elSize;
      for ( int x = 0; x < this->getNumX() - 1; x += h[d] ) {
        for ( int y = 0; y < this->getNumX() - 1; y += h[d] ) {
          int pos = x + y * this->getNumX();
          int type = 3;
          if ( ( x % elSizeBig ) == 0 ) { type -= 2; }
          if ( ( y % elSizeBig ) == 0 ) { type -= 1; }
          cEl.set ( x, y, 0, d, type );
          if ( elementInGrid ( cEl ) ) {
            continue;
          } else {
            if ( type == 0 ) {
              result[pos] += result[pos + elSize];
              result[pos + elSizeBig] += result[pos + elSize];

              result[pos] += result[pos + elSize * numY];
              result[pos + elSizeBig * numY] += result[pos + elSize * numY];

              result[pos] += result[pos + elSize + elSize * numY];
              result[pos + elSizeBig] += result[pos + elSize + elSize * numY];
              result[pos + elSizeBig * numY] += result[pos + elSize + elSize * numY];
              result[pos + elSizeBig + elSizeBig * numY] += result[pos + elSize + elSize * numY];
            }

            else if ( type == 1 ) {
              result[pos - elSize] += result[pos];
              result[pos + elSize] += result[pos];

              result[pos + elSize] += result[pos + elSize + elSize * numY];
              result[pos + elSize + elSizeBig * numY] += result[pos + elSize + elSize * numY];

              result[pos - elSize] += result[pos + elSize * numY];
              result[pos + elSize] += result[pos + elSize * numY];
              result[pos - elSize + elSizeBig * numY] += result[pos + elSize * numY];
              result[pos + elSize + elSizeBig * numY] += result[pos + elSize * numY];
            }

            else if ( type == 2 ) {
              result[pos - elSize] += result[pos];
              result[pos + elSize] += result[pos];

              result[pos + elSize ] += result[pos + elSize + elSize * numY];
              result[pos + elSizeBig + elSize * numY] += result[pos + elSize + elSize * numY];

              result[pos - elSize] += result[pos + elSize];
              result[pos + elSize] += result[pos + elSize];
              result[pos + elSizeBig - elSize * numY] += result[pos + elSize];
              result[pos + elSizeBig + elSize * numY] += result[pos + elSize];
            }

            else if ( type == 3 ) {
              result[pos - elSize + elSize * numY] += result[pos + elSize * numY];
              result[pos + elSize + elSize * numY] += result[pos + elSize * numY];

              result[pos + elSize + elSize * numY] += result[pos + elSize];
              result[pos + elSize - elSize * numY] += result[pos + elSize];

              result[pos - elSize - elSize * numY] += result[pos];
              result[pos + elSize - elSize * numY] += result[pos];
              result[pos - elSize + elSize * numY] += result[pos];
              result[pos + elSize + elSize * numY] += result[pos];
            } else {
              cerr << "BIG PROBLEM in adoptToAdaptiveGrid!!!" << endl;
            }
          }
        }
      }
    }
  }







  /* interpolateResult interpolates the result vector on all inactive
  * nodes
  */
  void addVirtualNodeToConstrainingNodes ( aol::Vector<double> &result ) const {
    int x, y, d;

    for ( d = 0; d < gridDepth; d++ ) {
      for ( x = 0; x < this->getNumX() - 1; x += h[d] ) {
        for ( y = 0; y < this->getNumX() - 1; y += h[d] ) {
          int k = x + map[y];

          int level = d;

          if ( level == gridDepth ) return;

          if ( isHangingNode ( k + cmap[level][0] ) ) {
            ( result ) [k] += ( result ) [k + cmap[level][0]] / 2.0;
            ( result ) [k + nmap[level][1]] += ( result ) [k + cmap[level][0]] / 2.0;
          }

          if ( isHangingNode ( k + cmap[level][3] ) ) {
            ( result ) [k] += ( result ) [k + cmap[level][3]] / 2.0;
            ( result ) [k + nmap[level][3]] += ( result ) [k + cmap[level][3]] / 2.0;
          }

          if ( isHangingNode ( k + cmap[level][1] ) ) {
            ( result ) [k + nmap[level][1]] += ( result ) [k + cmap[level][1]] / 2.0;
            ( result ) [k + nmap[level][2]] += ( result ) [k + cmap[level][1]] / 2.0;
          }

          if ( isHangingNode ( k + cmap[level][4] ) ) {
            ( result ) [k + nmap[level][3]] += ( result ) [k + cmap[level][4]] / 2.0;
            ( result ) [k + nmap[level][2]] += ( result ) [k + cmap[level][4]] / 2.0;
          }

          if ( isHangingNode ( k + cmap[level][2] ) ) {
            ( result ) [k + cmap[level][3]] += ( result ) [k + cmap[level][2]] / 2.0;
            ( result ) [k + cmap[level][1]] += ( result ) [k + cmap[level][2]] / 2.0;
          }
        }
      }
    }
  }



  void interpolateHangingNodes ( aol::Vector<double> &result ) const {
    int x, y, d;

    for ( d = 0; d < gridDepth; d++ ) {
      for ( x = 0; x < this->getNumX() - 1; x += h[d] ) {
        for ( y = 0; y < this->getNumX() - 1; y += h[d] ) {
          int k = x + map[y];

          int level = d;

          if ( level == gridDepth ) return;

          if ( isHangingNode ( k + cmap[level][0] ) )
            ( result ) [k + cmap[level][0]] = ( ( result ) [k] +
                                                ( result ) [k + nmap[level][1]] ) / 2.;

          if ( isHangingNode ( k + cmap[level][3] ) )
            ( result ) [k + cmap[level][3]] = ( ( result ) [k] +
                                                ( result ) [k + nmap[level][3]] ) / 2.;

          if ( isHangingNode ( k + cmap[level][1] ) )
            ( result ) [k + cmap[level][1]] =
              ( ( result ) [k + nmap[level][1]] + ( result ) [k + nmap[level][2]] ) / 2.;

          if ( isHangingNode ( k + cmap[level][4] ) )
            ( result ) [k + cmap[level][4]] =
              ( ( result ) [k + nmap[level][3]] + ( result ) [k + nmap[level][2]] ) / 2.;

          if ( isHangingNode ( k + cmap[level][2] ) )
            ( result ) [k + cmap[level][2]] =
              ( ( result ) [k + cmap[level][3]] + ( result ) [k + cmap[level][1]] ) / 2.;
        }
      }
    }
  }


  int getNumberOfNodes() const {
    return this->getNumX() * this->getNumY() * this->getNumZ();
  }


  /** Destructor
   */
  ~AdaptiveGridTP() {

    if ( nodeIsVertexOfElement != NULL ) delete nodeIsVertexOfElement;
    if ( nodeIsHangingNode != NULL ) delete nodeIsHangingNode;

    for ( int i = 0; i <= gridDepth; i++ ) {
      for ( int j = 0; j < 4; j++ ) {
        delete[] hmap[i][j];
      }
      delete[] nmap[i];
      delete[] cmap[i];
      delete[] hmap[i];
    }

    delete[] nmap;
    delete[] cmap;
    delete[] h;
    delete[] hSqr;
    delete[] spaceH;
    delete[] hSqr_2;
    delete[] map;
    delete[] hmap;
    delete elemIt;
  }

  AdaptiveGridTP<EstType, RealType> &begin() {
    return *this;
  }

  virtual bool isActive ( int k ) const {
    return ( ! ( nodeIsHangingNode->get ( k ) ) ) && ( nodeIsVertexOfElement->get ( k ) );
  }


  EstType* getEstimator() const {
    return _estimator;
  }

  /** Set error estimator that defines adaptivity for this 2D grid
   */
  void setEstimator ( EstType &e ) {
    _estimator = &e;
  }

  /** Return whether error estimator exists or not
   */
  bool isAdaptive() const {
    return ( _estimator != NULL );
  }

  bool isHangingNode ( int pos ) const {
    return nodeIsHangingNode->get ( pos );
  }

  void getActivityVector ( aol::Vector<RealType> &vec ) const {
    for ( int i = 0; i < vec.size(); i++ ) {
      if ( this->isActive ( i ) ) {
        vec.set ( i, 1.0 );
      } else {
        vec.set ( i, 0.0 );
      }
    }
  }


  void getElementVertexVector ( aol::Vector<RealType> &vec ) const {
    for ( int i = 0; i < vec.size(); i++ ) {
      if ( nodeIsVertexOfElement->get ( i ) ) {
        vec.set ( i, 1.0 );
      } else {
        vec.set ( i, 0.0 );
      }
    }
  }

  void getHangingNodeVector ( aol::Vector<RealType> &vec ) const {
    for ( int i = 0; i < vec.size(); i++ ) {
      if ( nodeIsHangingNode->get ( i ) ) {
        vec.set ( i, 1.0 );
      } else {
        vec.set ( i, 0.0 );
      }
    }
  }






  template <typename OpType>
  void getConstrainingNodes ( const OpType &op, aol::Vector<double> &constrainingNodes ) const {
    const typename ElementIterator::EndType end_it = op.getConfigurator().end();
    int globalDofs[ 4 ];
    int globalDofsConstraint[ 4 ];
    int globalDofsConstraintOtherElem[ 4 ];

    for ( ElementIterator it = op.getConfigurator().begin(); it != end_it; it++ ) {

      const int numLocalDofs = 4;
      for ( int i = 0; i < numLocalDofs; ++i ) {
        globalDofs[ i ] = op.getConfigurator().localToGlobal ( *it, i );
        globalDofsConstraintOtherElem[i] = globalDofs[i];
        globalDofsConstraint[i] = globalDofs[i];
        if ( op.getConfigurator().getGrid().isAdaptive() ) {
          const int elSize = 1 << ( op.getConfigurator().getGrid().getGridDepth() - ( *it ).level() );
          int hanging = op.getConfigurator().getGrid().checkForHangingNode ( *it, i );
          if ( hanging == i ) {
            globalDofsConstraintOtherElem[i] = globalDofs[i];
          } else {
            // constraining nodes in x direction
            if ( ( ( ( *it ).type() == 0 ) && ( i == 1 ) ) || ( ( ( *it ).type() == 1 ) && ( i == 0 ) ) || ( ( ( *it ).type() == 2 ) && ( i == 3 ) ) || ( ( ( *it ).type() == 3 ) && ( i == 2 ) ) ) {
              globalDofsConstraint[ i ] = op.getConfigurator().localToGlobal ( *it, hanging );
              if ( ( ( *it ).type() == 0 ) || ( ( *it ).type() == 2 ) ) {
                globalDofsConstraintOtherElem[i] = globalDofs[i] + 1 * elSize;
              } else {
                globalDofsConstraintOtherElem[i] = globalDofs[i] - 1 * elSize;
              }
            } else if ( ( ( ( *it ).type() == 0 ) && ( i == 2 ) ) || ( ( ( *it ).type() == 1 ) && ( i == 3 ) ) || ( ( ( *it ).type() == 2 ) && ( i == 0 ) ) || ( ( ( *it ).type() == 3 ) && ( i == 1 ) ) ) {
              globalDofsConstraint[ i ] = op.getConfigurator().localToGlobal ( *it, hanging );
              if (  ( *it ).type() == 0 || ( *it ).type() == 1 ) {
                globalDofsConstraintOtherElem[i] = globalDofs[i] + 1 * op.getConfigurator().getGrid().getNumX() * elSize;
              } else {
                globalDofsConstraintOtherElem[i] = globalDofs[i] - 1 * op.getConfigurator().getGrid().getNumX() * elSize;
              }
            } else {
              cerr << endl << endl << " BIG PROBLEM IN ASSEMBLEADDMATRIX!!!!!!!!!!!!!" << endl;
            }
          }
        }
      }

      // finally add the locally computed values to the matrix
      for ( int i = 0; i < numLocalDofs; ++i ) {


        int hanging = op.getConfigurator().getGrid().checkForHangingNode ( *it, i );
        if ( hanging != i ) {
          constrainingNodes.set ( globalDofsConstraint[i], 1.0 );
          constrainingNodes.set ( globalDofsConstraintOtherElem[i], 1.0 );
        }

      }
    }
  }






  double H() const {
    cerr << "H called on adaptive grid!!!" << endl;
    return 0.0;
  }


  template <typename MatrixType, typename VectorType>
  void removeInactiveNodesFromMatrixAndVector ( MatrixType &Matr, VectorType &vec ) const {
    int num = 0;
    for ( int i = 0; i < Matr.getNumCols(); i++ ) {
      if ( ! ( this->isActive ( i ) ) ) {
        Matr.set ( i, i, 1.0 );
        vec.set ( i, 0.0 );
        num++;
      }
    }
  }

  template <typename MatrixType>
  void removeInactiveNodesFromMatrix ( MatrixType &Matr ) const {
    for ( int i = 0; i < Matr.getNumCols(); i++ ) {
      if ( ! ( this->isActive ( i ) ) ) {
        Matr.setRowToZero ( i );
      }
    }
  }

  template <typename MatrixType>
  void setInactiveNodesFromMatrixToUnity ( MatrixType &Matr ) const {
    for ( int i = 0; i < Matr.getNumCols(); i++ ) {
      if ( ! ( this->isActive ( i ) ) ) {
        Matr.set ( i, i, 1.0 );
      }
    }
  }

  template <typename VectorType>
  void removeInactiveNodesFromVector ( VectorType &vec ) const {
    for ( int i = 0; i < vec.size(); i++ ) {
      if ( ! ( this->isActive ( i ) ) ) {
        vec.set ( i, 0.0 );
      }
    }
  }


  template <typename OpType>
  void initVertexArrays ( const OpType &op ) {
    int numElem = 0;
    int numHanging = 0;
    nodeIsVertexOfElement->setAll ( false );
    nodeIsHangingNode->setAll ( false );
    const typename ElementIterator::EndType end_it = op.getConfigurator().end();
    for ( ElementIterator it = op.getConfigurator().begin(); it != end_it; it++, numElem++ ) {
      const int numLocalDofs = 4;
      for ( int i = 0; i < numLocalDofs; ++i ) {
        int globalDof = op.getConfigurator().localToGlobal ( *it, i );
        nodeIsVertexOfElement->set ( globalDof, true );
      }
    }

    for ( ElementIterator it = op.getConfigurator().begin(); it != end_it; it++, numElem++ ) {
      const int numLocalDofs = 4;
      for ( int i = 0; i < numLocalDofs; ++i ) {
        int hanging = op.getConfigurator().getGrid().checkForHangingNode ( *it, i );
        int globalDof = op.getConfigurator().localToGlobal ( *it, i );
        if ( hanging != i ) {
          nodeIsHangingNode->set ( globalDof, true );
          numHanging++;
        }
      }
    }
  }

  int getNumNoVertexNodes() {
    int res = 0;
    int vertexNodes = 0;
    for ( int i = 0; i < nodeIsVertexOfElement->size(); i++ ) {
      if ( ! ( nodeIsVertexOfElement->get ( i ) ) ) {
        res++;
      } else {
        vertexNodes++;
      }
    }
    return res;
  }

  template <typename OpType>
  void drawMesh ( const OpType &op, qc::ScalarArray<int, qc::QC_2D> &img, int scale ) {

    img.setAll ( 255 );
    const typename ElementIterator::EndType end_it = op.getConfigurator().end();
    for ( ElementIterator it = op.getConfigurator().begin(); it != end_it; it++ ) {
      drawElement ( *it, img, scale );
    }

  }

  void drawElement ( const qc::Element &el, qc::ScalarArray<int, qc::QC_2D> &img, int scale ) {

    int i /*, k=el.x()+map[el.y()] */;

    /* Draw element first */
    for ( i = el.x() * scale; i < ( el.x() + h[el.level()] ) *scale; i++ ) {
      if ( img[el.y() * ( width - 1 ) *scale / scale + i] != 0 ) {
        img[el.y() * ( width - 1 ) *scale * scale + i] = 0;
      }
    }

    for ( i = el.x() * scale; i < ( el.x() + h[el.level()] ) *scale; i++ ) {
      if ( img[ ( el.y() + h[el.level()] ) * ( width - 1 ) *scale * scale + i] != 0 ) {
        img[ ( el.y() + h[el.level()] ) * ( width - 1 ) *scale * scale + i] = 0;
      }
    }

    for ( i = el.y() * scale * scale * ( width - 1 ); i < ( el.y() + h[el.level()] ) *scale * scale * ( width - 1 ); i += ( width - 1 ) * scale ) {
      if ( img[i + el.x() *scale] != 0 ) {
        img[i + el.x() *scale] =  0;
      }
    }

    for ( i = el.y() * scale * scale * ( width - 1 ); i < ( el.y() + h[el.level()] ) *scale * scale * ( width - 1 ); i += ( width - 1 ) * scale ) {
      if ( img[i + ( el.x() + h[el.level()] ) *scale] != 0 ) { /*/==0) */
        img[i + ( el.x() + h[el.level()] ) *scale] =  0;
      }
    }
  }


  virtual int checkForHangingNode ( const qc::Element &El, int VertInd ) const ;

protected:
  EstType     *_estimator;
  aol::BitVector *nodeIsVertexOfElement;
  aol::BitVector *nodeIsHangingNode;

};


template <class EstType, typename RealType>
int AdaptiveGridTP<EstType, RealType>::checkForHangingNode ( const qc::Element &El,
                                                             int VertInd ) const {
  qc::Element cEl;
  if ( !_estimator )
    throw aol::Exception ( "Set Error estimator in grid before "
                           "invoking checkForHangingNode!",
                           __FILE__,  __LINE__ );


  // Elements on level zero do never hang!
  if ( El.level() == 0 ) { return VertInd; }

  const int elSize = 1 << ( gridDepth - El.level() );

  const int maxNodes = ( 1 << gridDepth ) - elSize;

  switch ( this->getDimOfWorld() ) {
    case qc::QC_2D:

      // Check element type
      switch ( El.type() ) {
        case 0:
          if ( VertInd == 0 || VertInd == 3 ) return VertInd;
          if ( VertInd == 1 ) {
            cEl.set ( El.x(), aol::Max ( 0, El.y() - elSize ), 0, El.level(), 2 );
            if ( ! ( elementInGrid ( cEl ) ) )
              return 0;
          }
          if ( VertInd == 2 ) {
            cEl.set ( aol::Max ( 0, El.x() - elSize ), El.y(), 0,  El.level(), 1 );
            if ( ! ( elementInGrid ( cEl ) ) )
              return 0;
          }
          break;

        case 1:
          if ( VertInd == 1 || VertInd == 2 ) return VertInd;
          if ( VertInd == 0 ) {
            cEl.set ( El.x(), aol::Max ( 0, El.y() - elSize ), 0,  El.level(), 3 );
            if ( ! ( elementInGrid ( cEl ) ) )
              return 1;
          }
          if ( VertInd == 3 ) {
            cEl.set ( aol::Min ( maxNodes, El.x() + elSize ), El.y(), 0, El.level(), 0 );
            if ( ! ( elementInGrid ( cEl ) ) )
              return 1;
          }
          break;

        case 2:
          if ( VertInd == 1 || VertInd == 2 ) return VertInd;
          if ( VertInd == 0 ) {
            cEl.set ( aol::Max ( 0, El.x() - elSize ), El.y(), 0,  El.level(), 3 );
            if ( ! ( elementInGrid ( cEl ) ) )
              return  2;
          }
          if ( VertInd == 3 ) {
            cEl.set ( El.x(), aol::Min ( maxNodes, El.y() + elSize ), 0, El.level(), 0 );
            if ( ! ( elementInGrid ( cEl ) ) )
              return  2;
          }
          break;

        case 3:
          if ( VertInd == 0 || VertInd == 3 ) return VertInd;
          if ( VertInd == 1 ) {
            cEl.set ( aol::Min ( maxNodes, El.x() + elSize ), El.y(), 0, El.level(), 2 );
            if ( ! ( elementInGrid ( cEl ) ) )
              return  3;
          }
          if ( VertInd == 2 ) {
            cEl.set ( El.x(), aol::Min ( maxNodes, El.y() + elSize ), 0, El.level(), 1 );
            if ( ! ( elementInGrid ( cEl ) ) )
              return  3;
          }
          break;

        default:
          cerr << "unknown element type" << endl;
          ;
      }

      break;

    case qc::QC_3D:
      throw aol::Exception ( "checkForHangingNode not implemented for 3D!", __FILE__,
                             __LINE__ );

      break;

    default:
      throw aol::Exception ( "AdaptiveGridTP::checkForHangingNodes not implemented for dimension other than 2 or 3", __FILE__, __LINE__ );
  }
  return VertInd;

}


template <class EstType, typename RealType>
inline typename AdaptiveGridTP<EstType, RealType>::ElementIterator&
AdaptiveGridTP<EstType, RealType>::ElementIterator::operator++ ( int ) {

  if ( _it_stack.empty() ) {
    _more = false;
    return *this;
  }

  qc::Element cur = _it_stack.back();
  _it_stack.pop_back();
  if ( cur.level() < _grid_depth && ( _estimator->checkElement ( cur ) ) /*|| (cur.level() +2  < _grid_depth)*/ ) {
    const int nl = cur.level() + 1;
    const int x = cur.x(), y = cur.y(), z = cur.z();
    const int X = x + _grid->h[nl];
    const int Y = y + _grid->h[nl];
    //const int Z = z + _grid->h[nl];

    cur.set ( x, y, z, nl, 0 );
    _it_stack.push_back ( cur );
    cur[0] = X;
    cur.setType ( 1 );
    _it_stack.push_back ( cur );
    cur[1] = Y;
    cur.setType ( 3 );
    _it_stack.push_back ( cur );
    cur[0] = x;

    /*
    if(_grid->getDimOfWorld() == qc::QC_3D){
      cur.setType ( 2 );
      _it_stack.push_back ( cur );
      cur[2] = Z;
      cur.setType ( 6 );
      _it_stack.push_back ( cur );
      cur[1] = y;
      cur.setType ( 4 );
      _it_stack.push_back ( cur );
      cur[0] = X;
      cur.setType ( 5 );
      _it_stack.push_back ( cur );
      cur[1] = Y;
    }
    else{
    */
    cur.setType ( 2 );
    _it_stack.push_back ( cur );
    this->operator++ ( 1 );
    return *this;
    //}
  }
  _cur = cur;
  return *this;
}


template <class EstType, typename RealType>
inline typename AdaptiveGridTP<EstType, RealType>::NodeIterator&
AdaptiveGridTP<EstType, RealType>::NodeIterator::operator= ( AdaptiveGridTP<EstType, RealType> &Grid ) {
  _grid = &Grid;
  const int width = _grid->getWidth();
  if ( !_array || _array->getNumX() != width ) {
    if ( _array ) delete _array;
    _array = new qc::Array<unsigned char> ( width, width, width );
  }
  int num = 0;
  cerr << endl;
  // initialize the array by traversing local element nodes
  for ( ElementIterator eit = Grid.begin(); eit != Grid.end(); ++eit ) {
    const int &h = _grid->h[ eit->level() ];
    const int x = eit->x(), y = eit->y(), z = eit->z();
    const int X = x + h, Y = y + h, Z = z + h;
    _array->set ( x, y, z, 1 );
    _array->set ( X, y, z, 1 );
    _array->set ( x, Y, z, 1 );
    _array->set ( X, Y, z, 1 );
    _array->set ( x, y, Z, 1 );
    _array->set ( X, y, Z, 1 );
    _array->set ( x, Y, Z, 1 );
    _array->set ( X, Y, Z, 1 );
    num++;
  }
  cerr << "NodeField initialized.num = " << num << "\n";
  return *this;
}

#endif
