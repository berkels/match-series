#ifndef __SIMPLEXBASEFUNCSETTFE_H
#define __SIMPLEXBASEFUNCSETTFE_H

#include <vec.h>
#include <smallMat.h>
#include <matrix.h>
#include <simplexLookup.h>
#include <simplexBaseFunctionSet.h>

namespace qc {

namespace simplex {

//-----------------------------------------------------------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType, qc::Dimension Dim>
class BaseFuncSetTFEConstructor;

//=============================================================================================================================

template <typename RealType, typename VecType, typename DomVecType, int NumBaseFuncs, typename Imp>
class VariableNumberOfQuadPointsBaseFunctionSet {
public:
  static const int numBaseFunctions = NumBaseFuncs;

  int numQuadPoints( ) const {
    return _numQuadPoints;
  }

  inline RealType getWeight ( int QuadPoint ) const {
    return _weights[QuadPoint];
  }

  inline const DomVecType & getRefCoord ( int QuadPoint ) const {
    return _quadPoints[QuadPoint];
  }

  inline RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return _basisQuadValues[BaseFuncNum][QuadPoint];
  }

  inline const VecType & evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return _basisQuadGradients[BaseFuncNum][QuadPoint];
  }

protected:
  void initializeQuadCache () {
    for (int b = 0; b < Imp::numBaseFunctions; ++b) {
      vector<RealType> thisBFSQuadValues ( _numQuadPoints );
      vector<VecType>  thisBFSQuadGradients ( _numQuadPoints );
      for (int q = 0; q < _numQuadPoints; ++q) {
        thisBFSQuadValues[q] = asImp().evaluate(b, _quadPoints[q] );
        asImp().evaluateGradient ( b, _quadPoints[q], thisBFSQuadGradients[q] );
      }
      _basisQuadValues.push_back ( thisBFSQuadValues );
      _basisQuadGradients.push_back ( thisBFSQuadGradients );
    }
  }

  const Imp & asImp () const {
    return static_cast<const Imp &> ( *this );
  }

  /* _quadPoints and _weights are vectors of length _numQuadPoints.
     They have to be filled by the derived class before
     initializeQuadCache() is called.                               */
  int                _numQuadPoints;
  vector<DomVecType> _quadPoints;
  vector<RealType>   _weights;

  vector<vector<RealType> > _basisQuadValues;
  vector<vector<VecType> >  _basisQuadGradients;
};

//=============================================================================================================================

template <typename RealType, qc::Dimension Dim, typename QuadRuleType>
class BaseFunctionSetTFE
  : public VariableNumberOfQuadPointsBaseFunctionSet < RealType,
                                                       typename aol::VecDimTrait<RealType, Dim>::VecType,
                                                       aol::BarCoord<Dim, RealType>, /* NumBaseFunctions = */ Dim + 1,
                                                       BaseFunctionSetTFE<RealType, Dim, QuadRuleType> > {

  friend class BaseFuncSetTFEConstructor<RealType, QuadRuleType, Dim>;

  typedef VariableNumberOfQuadPointsBaseFunctionSet < RealType,
            typename aol::VecDimTrait<RealType, Dim>::VecType,
            aol::BarCoord<Dim, RealType>, Dim + 1,
            BaseFunctionSetTFE<RealType, Dim, QuadRuleType> >    Base;

public:
  typedef typename aol::VecDimTrait<RealType, Dim>::VecType      VecType;
  typedef aol::BarCoord<Dim, RealType>                           DomVecType;

  //! public constructor for non-cut elements
  BaseFunctionSetTFE ( RealType H, int SimplexNumber )
  : _uncutBFS ( H, SimplexNumber ) {
    this->_numQuadPoints = QuadRuleType::numQuadPoints;
    this->_quadPoints.resize ( this->_numQuadPoints );
    this->_weights.resize ( this->_numQuadPoints );

    for (int i = 0; i < this->_numQuadPoints; ++i) {
      this->_quadPoints[i] = QuadRuleType::getRefCoord ( i );
      this->_weights[i] = QuadRuleType::getWeight ( i );
    }

    // as long as no interface is set explicitly, we
    // assume that the whole element lies inside the band.
    for (int i = 0; i < Dim + 1; ++i)
      _interface[i] = - aol::ZOTrait<RealType>::one;

    _volume = ( Dim == QC_2D ? aol::Sqr ( H ) : aol::Cub ( H ) )
              / static_cast<RealType>(TopologyLookup<Dim>::numSimplexesPerCube);
  }

  using Base::evaluate;
  using Base::evaluateGradient;

  inline RealType evaluate ( int BaseFuncNum, const DomVecType & RefCoord ) const {
    if (_interface * RefCoord <= 0.)
      return _uncutBFS.evaluate ( BaseFuncNum, RefCoord );
    else
      return aol::ZOTrait<RealType>::zero;
  }

  inline void evaluateGradient ( int BaseFuncNum, const DomVecType & RefCoord, VecType & Gradient ) const {
    if (_interface * RefCoord <= 0.)
      _uncutBFS.evaluateGradient ( BaseFuncNum, RefCoord, Gradient );
    else
      Gradient.setZero();
  }

  RealType getVolume () const {
    return _volume;
  }

protected:
  //! constructor to be called by BaseFuncSetCFEConstructor
  BaseFunctionSetTFE ( RealType H, int SimplexNumber, const aol::Vec<Dim + 1, RealType> & Interface )
  : _interface ( Interface )
  , _volume ( 0. )
  , _uncutBFS ( H, SimplexNumber ) {}

  //! the basis functions are only nonzero in the region
  //! {lambda | interface * lambda <= 0}.
  aol::Vec<Dim + 1, RealType>                     _interface;
  RealType                                        _volume;
  BaseFunctionSetLin<RealType, Dim, QuadRuleType> _uncutBFS;
};

//=============================================================================================================================

template <typename RealType, typename QuadRuleType, qc::Dimension Dim>
class BaseFuncSetTFEConstructor {};

//-----------------------------------------------------------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
class BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D> {

public:
  static const Dimension Dim = QC_3D;
  typedef BaseFunctionSetTFE<RealType, Dim, QuadRuleType> BFS;
  typedef Lookup<RealType, Dim>                           LookupType;
  typedef aol::BarCoord<Dim, RealType>                    DomVecType;

  static const int numSimplexes = LookupType::numSimplexesPerCube;

  BaseFuncSetTFEConstructor ( RealType H );
  BFS * getNewBFS ( const aol::Vec<4, RealType> & LevelValues, short SimplexNumber ) const;

protected:
  int signature ( const aol::Vec<4, RealType> & LevelValues ) const;
  int edgeBetween ( int node_a, int node_b ) const;
  void computeCutRelations ( const aol::Vec<4, RealType> & LevelValues, int sigma,
                             RealType cutRelations[], bool interfaced[] ) const;
  int findNextUncutEdge ( bool interfaced[], int start ) const;
  const DomVecType & convexComb ( int a, int b, RealType lambda ) const;

  void addTetra ( BFS * bfs,
                  const DomVecType & e0, const DomVecType & e1,
                  const DomVecType & e2, const DomVecType & e3,
                  const aol::Mat<3, 4, RealType> & M ) const;

  inline void addTetra ( BFS * bfs, const DomVecType * e[],
                         const aol::Mat<3, 4, RealType> & M ) const;

  int                              _allSet;
  RealType                         _h;
  aol::Mat<Dim, Dim + 1, RealType> _cartesianNodeCoords[numSimplexes];
};

//-----------------------------------------------------------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D>::BaseFuncSetTFEConstructor( RealType H )
: _h ( H ) {
  aol::Vec<Dim + 1, RealType> allNeg;
  allNeg.setAll ( -aol::ZOTrait<RealType>::one );
  _allSet = signature ( allNeg );

  // fill _cartesianNodeCoords with cartesian node coordinates
  // of the simplex'es nodes. Use therefore the list of cube node
  // coordinates and an index translation simplex->cube.
  for (int s = 0; s < numSimplexes; ++s)
    for (int n = 0; n < Dim + 1; ++n)
      for (int d = 0; d < Dim; ++d)
        _cartesianNodeCoords[s][d][n] = LookupType::cartesianCubeNodeCoords[LookupType::localIndicesSimplexToCube[s][n]][d];
}

//-----------------------------------------------------------------------------------------------------------------------------

//! create a new CFE base function set
//! \param LevelValues values of such a level function that the sub-zero level
//!                    is thought to lie inside the band
//! \param SimplexNumber number of simplex in cube (0..1 for 2D, 0..5 for 3D)
//! \return pointer to base function set object on the heap (has to be
//!         destructed by caller!)
template <typename RealType, typename QuadRuleType>
typename BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D>::BFS *
BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D>::
getNewBFS ( const aol::Vec<4, RealType> & LevelValues, short SimplexNumber ) const {
  int sigma = signature ( LevelValues );

  // this function should not be called if the element lies
  // outside the band
  if (sigma == 0)
    throw aol::Exception ( "BaseFuncSetCFEConstructor::getNewBFS() "
                           "was called for an element outside the band.",
                           __FILE__, __LINE__ );


  // if all nodes lie inside, return standard-constructed BFS
  if (sigma == _allSet) {
    BFS * bfs = new BFS ( _h, SimplexNumber );
    bfs->initializeQuadCache();
    return bfs;
  }

  // call protected constructor which does not initialize
  // the quad point coordinates and weights
  BFS * bfs = new BFS ( _h, SimplexNumber, LevelValues );

  // choose appropriate cartesian node coordinate collection M
  const aol::Mat<3, 4, RealType> & M = _cartesianNodeCoords[SimplexNumber];

  // determine which edges are cut and compute corresponding cut relations
  RealType cutRelations[LookupType::maxEdgeBit + 1];
  bool     interfaced[LookupType::maxEdgeBit + 1];
  computeCutRelations ( LevelValues, sigma, cutRelations, interfaced );

  int a, b, c, d, indexSum = Dim * (Dim + 1) / 2;

  // save nodes of first not-cut edge
  int e = findNextUncutEdge ( interfaced, 0 );
  a = LookupType::edges[e][0];
  b = LookupType::edges[e][1];

  e = findNextUncutEdge ( interfaced, ++e );

  // if one node of the second non-cut edges has already been saved,
  // all non-cut edges have the same sign (tetra-prism situation).
  // Because the LookupType::edges is sorted in the node inidices, it
  // suffices to test the first node of e.
  if (LookupType::edges[e][0] == a || LookupType::edges[e][0] == b) {
    c = LookupType::edges[e][1];
    d = indexSum - (a + b + c);

    // create barycentric coordinate vectors for all nodes
    // (BarCoord constructor with int argument sets the corresponding
    // coordinate to one, all others to zero)
    const DomVecType barCoord_a ( a ),
                     barCoord_b ( b ),
                     barCoord_c ( c ),
                     barCoord_d ( d ),
                     barCoord_da = convexComb ( a, d, cutRelations[edgeBetween(a, d)] ),
                     barCoord_db = convexComb ( b, d, cutRelations[edgeBetween(b, d)] ),
                     barCoord_dc = convexComb ( c, d, cutRelations[edgeBetween(c, d)] );

    // if node d is negative, (a, b, c) are positive, situation + + + -
    if (sigma & 1 << d) {
      addTetra ( bfs, barCoord_d, barCoord_da, barCoord_db, barCoord_dc, M );
    }
    // otherwise: situation - - - +
    else {
      addTetra ( bfs, barCoord_a,  barCoord_b,  barCoord_db, barCoord_dc, M );
      addTetra ( bfs, barCoord_b,  barCoord_c,  barCoord_a,  barCoord_dc, M );
      addTetra ( bfs, barCoord_a,  barCoord_da, barCoord_db, barCoord_dc, M );
    }
  }

  // if the two non-intersected edges have no node in common,
  // there are only these two and have different sign
  // (prism-prism situation).
  else {
    c = LookupType::edges[e][0];
    d = indexSum - (a + b + c); /* why not LookupType::edges[e][1] ?? */

    // create barycentric coordinate vectors for all nodes
    // (BarCoord constructor with int argument sets the corresponding
    // coordinate to one, all others to zero)
    const DomVecType barCoord_a ( a ),
                     barCoord_b ( b ),
                     barCoord_c ( c ),
                     barCoord_d ( d ),
                     barCoord_ac = convexComb ( a, c, cutRelations[edgeBetween(a, c)] ),
                     barCoord_ad = convexComb ( a, d, cutRelations[edgeBetween(a, d)] ),
                     barCoord_bc = convexComb ( b, c, cutRelations[edgeBetween(b, c)] ),
                     barCoord_bd = convexComb ( b, d, cutRelations[edgeBetween(b, d)] );

    // if nodes a and b are negative, c and d are positive, situation - - + +
    if (sigma & edgeBetween(a, b)) {
      addTetra ( bfs, barCoord_a,  barCoord_b,  barCoord_bd, barCoord_bc, M );
      addTetra ( bfs, barCoord_a,  barCoord_ac, barCoord_ad, barCoord_bc, M );
      addTetra ( bfs, barCoord_a,  barCoord_bd, barCoord_ad, barCoord_bc, M );
    }

    // otherwise: inverse situation + + - -
    else {
      addTetra ( bfs, barCoord_c,  barCoord_d,  barCoord_ad, barCoord_bd, M );
      addTetra ( bfs, barCoord_c,  barCoord_bd, barCoord_ad, barCoord_bc, M );
      addTetra ( bfs, barCoord_c,  barCoord_ac, barCoord_ad, barCoord_bc, M );
    }
  }

  // reset bfs->_numQuadPoints appropriately
  bfs->_numQuadPoints = static_cast<int>(bfs->_quadPoints.size());

  // finally, all quad point weights have to be divided by bfs->_volume to obtain
  // weights relative to this volume (as the FEOps expect)
  for (size_t i = 0; i < bfs->_weights.size(); ++i)
    bfs->_weights[i] /= bfs->_volume;

  bfs->_volume *= aol::Cub ( _h );

  bfs->initializeQuadCache();
  return bfs;
}

//-----------------------------------------------------------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
void BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D>::
computeCutRelations ( const aol::Vec<4, RealType> & LevelValues, int sigma,
                      RealType cutRelations[], bool interfaced[] ) const {
  for (int e = 0; e < LookupType::numEdges; ++e) {
    int edgeAndSignature = LookupType::edgesAsBits[e] & sigma;
    // edge is cut if its nodes have different sign, i. e. in this AND-ing
    // exactly one of the two egde's nodes' bits is set.
    if (edgeAndSignature != 0 && edgeAndSignature != LookupType::edgesAsBits[e]) {
      RealType a = LevelValues[LookupType::edges[e][0]],
               b = LevelValues[LookupType::edges[e][1]];
      cutRelations[LookupType::edgesAsBits[e]] = a / (a - b);
      interfaced[LookupType::edgesAsBits[e]] = true;
    }
    else {
      cutRelations[LookupType::edgesAsBits[e]] = - aol::ZOTrait<RealType>::one;
      interfaced[LookupType::edgesAsBits[e]] = false;
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
int BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D>::
signature ( const aol::Vec<4, RealType> & LevelValues ) const {
  short ret = 0;
  for (int i = 0; i < Dim + 1; ++i)
    if (LevelValues[i] <= 0)
      ret += (1 << i);
  return ret;
}

//-----------------------------------------------------------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
int BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D>::
edgeBetween ( int node_a, int node_b ) const {
  return (1 << node_a) + (1 << node_b);
}

//-----------------------------------------------------------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
int BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D>::
findNextUncutEdge ( bool interfaced[], int start ) const {
  while ( ( start < LookupType::numEdges ) && interfaced[LookupType::edgesAsBits[start++]]) ;
  return --start;
}
//-----------------------------------------------------------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
const typename BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D>::DomVecType &
BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D>::
convexComb ( int a, int b, RealType lambda ) const {
  int iSmall = (a < b ? a : b);
  int iLarge = (a < b ? b : a);
  static DomVecType ret;
  ret.setZero();
  ret[iSmall] = aol::ZOTrait<RealType>::one - lambda;
  ret[iLarge] = lambda;
  return ret;
}

//-----------------------------------------------------------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
void BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D>::
addTetra ( BFS * bfs, const typename BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D>::DomVecType & e0,
                      const typename BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D>::DomVecType & e1,
                      const typename BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D>::DomVecType & e2,
                      const typename BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D>::DomVecType & e3,
           const aol::Mat<3, 4, RealType> & M ) const {

  const DomVecType * e[] = { &e0, &e1, &e2, &e3 };
  addTetra ( bfs, e, M );
}

//-----------------------------------------------------------------------------------------------------------------------------

template <typename RealType, typename QuadRuleType>
void BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D>::
addTetra ( BFS * bfs, const typename BaseFuncSetTFEConstructor<RealType, QuadRuleType, QC_3D>::DomVecType * e[],
           const aol::Mat<3, 4, RealType> & M ) const {
  aol::Vec<Dim, RealType> Me0 = M * *(e[0]);
  typename aol::MatDimTrait<RealType, Dim>::MatType cartCoord;
  for (int i = 0; i < Dim; ++i) {
    cartCoord[i] = M * *(e[i + 1]);
    cartCoord[i] -= Me0;
  }

  RealType det = cartCoord.det();
  RealType vol = fabs ( det ) / static_cast<RealType>(TopologyLookup<QC_3D>::numSimplexesPerCube);

  DomVecType quadPt;
  for (int i = 0; i < QuadRuleType::numQuadPoints; ++i) {
    quadPt.setZero();
    const DomVecType & regQuadPt = QuadRuleType::getRefCoord(i);
    for (int d = 0; d < Dim + 1; ++d)
      quadPt.addMultiple ( *(e[d]), regQuadPt[d] );
    bfs->_quadPoints.push_back ( quadPt );
    bfs->_weights.push_back ( vol * QuadRuleType::getWeight ( i ) );
  }
  bfs->_volume += vol;
}

//=============================================================================================================================

template <typename RealType, typename QuadRuleType, Dimension Dim>
class BaseFunctionCollection {
public:
  typedef BaseFunctionSetTFE<RealType, Dim, QuadRuleType> BFS;

  //! constructor
  template <typename GridType, typename IndexMapperType>
  BaseFunctionCollection ( const GridType & Grid, const aol::Vector<RealType> & LevelValues, RealType BandRadius,
                           const IndexMapperType & IndexMapper ) {
    BaseFuncSetTFEConstructor<RealType, QuadRuleType, Dim> bfsConstructor ( Grid.H() );

    // numberOfNodes times numSimplexesPerCube is a bit too large, but we do not care.
    _bfs.resize ( Grid.getNumberOfNodes() * TopologyLookup<Dim>::numSimplexesPerCube, NULL );

    // run over all elements, create appropriate small level values vector
    // and construct BFS with it.
    for (typename GridType::OldFullElementIterator iter = Grid.begin(); iter != Grid.end(); ++iter) {
      aol::Vec<Dim + 1, RealType> elementLevelValues;
      for (int i = 0; i < Dim + 1; ++i)
        elementLevelValues[i] = LevelValues[IndexMapper.localToGlobal ( *iter, i )];

      careForNearlyOutlyingNodes ( elementLevelValues, BandRadius );

      int type = signedType ( elementLevelValues, BandRadius );

      // if no bit in "type" is set, the element has level values below -BandRadius
      // and above BandRadius, thus the given band radius is illegal
      if (type == 0)
        throw aol::Exception ( aol::strprintf("BaseFunctionSetCFE: found element with level values "
                                "below -r and above r, band radius r = %lf", BandRadius),
                                __FILE__, __LINE__ );

      // do not create a base function set for this element if all
      // nodes lie outside the band (if so, all nodes will lie
      // either above or below the band's level values)
      if (type == tooLarge || type == tooSmall)
        continue;

      // if the second bit of "type" is set, the element is crossed by
      // by the -BandRadius level set. To obtain level values that are
      // larger than BandRadius outside and smaller inside, we replace
      // invert the level values.
      if (type == 2)
        for (int i = 0; i < Dim + 1; ++i)
          elementLevelValues[i] = - elementLevelValues[i];

      // now substract BandRadius from the level values to obtain
      // level values <= 0 inside and > 0 outside the band,
      // as expected by the BaseFuncSetCFEConstructor:
      for (int i = 0; i < Dim + 1; ++i)
        elementLevelValues[i] -= BandRadius;

      _bfs[Grid.getElementIndex(*iter)] = bfsConstructor.getNewBFS ( elementLevelValues, iter->getSimplexNumber() );
    }
  }

  //! this object is only shallow-copied. Copied object will behave well only
  //! until the initial object or the copy is destructed (or re-initialized,
  //! which is not implemented yet). When the second destructor is called,
  //! your program will probably crash.
  BaseFunctionCollection ( const BaseFunctionCollection<RealType, QuadRuleType, Dim> & other ) {
    cerr << aol::color::red << "You shall not have any BFS before me." << aol::color::reset << endl;
    this->_bfs = other._bfs;
  }

  ~BaseFunctionCollection () {
    for (size_t i = 0; i < _bfs.size(); ++i)
      delete _bfs[i];
  }

  const BFS & getBaseFunctionSet ( int Index ) const {
    return *_bfs[Index];
  }

protected:
  int signedType ( const aol::Vec<Dim + 1, RealType> & LevelValues, RealType BandRadius ) const {
    int ret = 7; // bits 1 1 1
    for (int i = 0; i < Dim + 1; ++i)
      if (LevelValues[i] >= BandRadius)
        ret &= tooLarge;
      else if (LevelValues[i] <= -BandRadius)
        ret &= tooSmall;
      else
        ret &= inbetween;
    return ret;
  }

  //! pushes values very near to +-BandRadius into the band
  void careForNearlyOutlyingNodes ( aol::Vec<Dim + 1, RealType> & x, RealType r, RealType TOL = 1E-12 ) const {
    for (int i = 0; i < Dim + 1; ++i)
      if ( x[i] < r + TOL && x[i] > r - TOL)
        x[i] = r - TOL;
      else if (x[i] > -(r + TOL) && x[i] < -(r - TOL) )
        x[i] = -(r - TOL);
  }

  static const int tooLarge  = 5, // bits 1 0 1
                   inbetween = 3, // bits 0 1 1
                   tooSmall  = 6; // bits 1 1 0

  // vector of BFS pointers. Each entry has to be
  // deleted when resizing the vector or destructing
  // the BaseFunctionCollection object. The entries
  // are indexed over return values of Grid.getElementIndex().
  vector<BFS *> _bfs;
};

//-----------------------------------------------------------------------------------------------------------------------------

} // end of namespace simplex.

} // end of namespace qc.

#endif
