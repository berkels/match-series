#include <linearSmoothOp.h>

#include <gridBase.h>
#include <multilevelArray.h>
#include <cellCenteredGrid.h>

namespace {
// nameless namespace only visible in this cpp file

// some old finite-difference based helper class for the linearSmoothOp
template <typename RealType>
class HeatEquation2DFD : public aol::Op<aol::Vector<RealType> > {
protected:
  const short    _width;
  RealType      _tau;
  short         _a2d_off1[4];

  int indexOffset ( int dx, int dy ) const {
    return dx + dy * _width;
  }

public:

  template <typename GridType>
  HeatEquation2DFD ( const GridType &Grid, RealType Tau )
      : _width ( Grid.getWidth() ),
        _tau ( Tau ) {
    init();
  }

  virtual ~HeatEquation2DFD() {}

  /** Set the time step for this heat equation
   */
  void setTimeStep ( RealType Tau ) {
    _tau = Tau;
    init();
  }

  /** Get the time step of this heat equation
   */
  RealType getTimeStep() const {
    return _tau;
  }

  void init() {
    _a2d_off1[ 0 ] = indexOffset ( + 1,  0 );
    _a2d_off1[ 1 ] = indexOffset ( -1,  0 );
    _a2d_off1[ 2 ] = indexOffset ( 0, + 1 );
    _a2d_off1[ 3 ] = indexOffset ( 0, -1 );
  }

  virtual void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp ( Dest.size() );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  virtual void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    // If we are on a 1x1 grid, we can't do anything.
    if ( _width == 1 ) {
      Dest = Arg;
      return;
    }

    const RealType _hsqrtau = _tau * aol::Sqr ( _width - 1 );


    const qc::Array<RealType> ArgIm ( Arg, _width, _width );
    qc::Array<RealType> DestIm ( Dest, _width, _width );

    RealType a, b;
    int X, Y, i;

    const int off[ 4 ] = { _a2d_off1[0], _a2d_off1[1], _a2d_off1[2], _a2d_off1[3] };

    // traverse the interior nodes.
    // minimize number of multiplication and index computations
    for ( X = 1; X < _width - 1; X++ ) {
      for ( Y = 1; Y < _width - 1; Y++ ) {
        int index = ArgIm.index ( X, Y );
        a  = ArgIm.get ( index ) * 4;

        b = 0.0;
        for ( i = 0; i < 4; i++ ) {
          b -= ArgIm.get ( index + off[ i ] );
        }
        DestIm.set ( X, Y, ArgIm.get ( index ) + _hsqrtau * ( a + b ) );
      }
    }
    const short e = _width - 1;
    const short corners[4][2] = { { 0, 0 }, {e, 0}, {0, e}, {e, e} };
    const short co_offs[4][2] = { { 1, 1 }, { -1, 1}, {1, -1}, { -1, -1} };

    for ( i = 0; i < 4; i++ ) {
      const short *c = corners[i];
      const short *o = co_offs[i];

      a = ArgIm.get ( c[0],  c[1] ) * 4;

      b = - ArgIm.get ( c[0] + o[0],  c[1] )
          - ArgIm.get ( c[0],  c[1] + o[1] );

      DestIm.set ( c[0], c[1], ArgIm.get ( c[0], c[1] ) + _hsqrtau * ( a + 2 * b ) );
    }

    // traverse the 2 edges parallel to X-axis;
    for ( X = 1; X < _width - 1; X++ ) {
      a  = ArgIm.get ( X,  e ) * 4;

      b = -2 * ( ArgIm.get ( X,  e - 1 ) );
      b -= ArgIm.get ( X + 1,  e );
      b -= ArgIm.get ( X - 1,  e );

      DestIm.set ( X, e, ArgIm.get ( X, e ) + _hsqrtau * ( a + b ) );

      a  = ArgIm.get ( X,  0 ) * 4;

      b = -2 * ( ArgIm.get ( X,  1 ) );
      b -= ArgIm.get ( X + 1,  0 );
      b -= ArgIm.get( X - 1,  0 );

      DestIm.set ( X, 0, ArgIm.get ( X,  0 ) + _hsqrtau * ( a + b ) );
    }

    // traverse the 2 edges parallel to Y-axis;
    for ( X = 1; X < _width - 1; X++ ) {

      a  = ArgIm.get ( e,  X ) * 4;

      b = -2 * ( ArgIm.get ( e - 1,  X ) );
      b -= ArgIm.get ( e,  X + 1 );
      b -= ArgIm.get ( e,  X - 1 );

      DestIm.set ( e, X, ArgIm.get ( e, X ) + _hsqrtau * ( a + b ) );

      a  = ArgIm.get ( 0,  X ) * 4;

      b = -2 * ( ArgIm.get ( 1,  X ) );
      b -= ArgIm.get ( 0,  X + 1 );
      b -= ArgIm.get ( 0,  X - 1 );

      DestIm.set ( 0, X, ArgIm.get ( 0, X ) + _hsqrtau * ( a + b ) );

    }
  }
};

// some old finite-difference based helper class for the linearSmoothOp
template <typename RealType>
class HeatEquation3DFD : public aol::Op<aol::Vector<RealType> > {
protected:
  const int         _width;
  RealType          _tau;
  int _a3d_off1[6];

  int indexOffset ( int dx, int dy, int dz ) const {
    return dx + dy * _width + dz * aol::Sqr ( _width );
  }

public:

  template <typename GridType>
  HeatEquation3DFD ( const GridType &Grid, RealType Tau )
    : _width ( Grid.getWidth() ),
      _tau ( Tau ) {
    init();
  }

  virtual ~HeatEquation3DFD() {}

  /** Set the time step for this heat equation
   */
  void setTimeStep ( RealType Tau ) {
    _tau = Tau;
    init();
  }

  /** Get the time step of this heat equation
   */
  RealType getTimeStep() const {
    return _tau;
  }

  void init() {
    _a3d_off1[ 0 ] = indexOffset ( + 1,  0,  0 );
    _a3d_off1[ 1 ] = indexOffset ( -1,  0,  0 );
    _a3d_off1[ 2 ] = indexOffset ( 0, + 1,  0 );
    _a3d_off1[ 3 ] = indexOffset ( 0, -1,  0 );
    _a3d_off1[ 4 ] = indexOffset ( 0,  0, + 1 );
    _a3d_off1[ 5 ] = indexOffset ( 0,  0, -1 );
  }

  virtual void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp ( Dest.size() );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  virtual void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    // If we are on a 1x1x1 grid, we can't do anything.
    if ( _width == 1 ) {
      Dest = Arg;
      return;
    }

    const RealType _hsqrtau = _tau * aol::Sqr ( _width - 1 );

    const qc::Array<RealType> ArgIm ( Arg, _width, _width, _width );
    qc::Array<RealType> DestIm ( Dest, _width, _width, _width );

    RealType a, b;
    int X, Y, Z, i;

    const int off[ 6 ] = {
                           _a3d_off1[0], _a3d_off1[1], _a3d_off1[2],
                           _a3d_off1[3], _a3d_off1[4], _a3d_off1[5]
                         };

    // traverse the interior nodes.
    // minimize number of multiplication and index computations
    for ( X = 1; X < _width - 1; X++ ) {
      for ( Y = 1; Y < _width - 1; Y++ ) {
        for ( Z = 1; Z < _width - 1; Z++ ) {
          int index = ArgIm.index ( X, Y, Z );
          a  = ArgIm.get ( index ) * 6;

          b = 0.0;
          for ( i = 0; i < 6; i++ ) {
            b -= ArgIm.get ( index + off[ i ] );
          }
          DestIm.set ( X, Y, Z, ArgIm.get ( index ) + _hsqrtau * ( a + b ) );
        }
      }
    }


    const int w = _width - 1;
    const int corners[8][3] = { { 0, 0, 0 }, {w, 0, 0}, {0, w, 0}, {w, w, 0}, {0, 0, w}, {w, 0, w}, {0, w, w}, {w, w, w} };
    const int co_offs[8][3] = { { 1, 1, 1 }, { -1, 1, 1}, {1, -1, 1}, { -1, -1, 1}, {1, 1, -1}, { -1, 1, -1}, {1, -1, -1}, { -1, -1, -1} };
    for ( i = 0; i < 8; i++ ) {
      const int *c = corners[i];
      const int *o = co_offs[i];

      a = ArgIm.get ( c[0], c[1],  c[2] ) * 6;

      b = - ArgIm.get ( c[0] + o[0], c[1],  c[2] )
          - ArgIm.get ( c[0], c[1] + o[1],  c[2] )
          - ArgIm.get ( c[0], c[1],  c[2] + o[2] );

      DestIm.set ( c[0], c[1], c[2], ArgIm.get ( c[0], c[1], c[2] ) + _hsqrtau * ( a + 2 * b ) );
    }

    const int edges [4][2] = { { 0, 0}, {0, w}, {w, 0}, {w, w} };
    const int e_offs[4][2] = { { 1, 1}, {1, -1}, { -1, 1}, { -1, -1} };

    // traverse the 4 edges parallel to x-axis;
    for ( i = 0; i < 4; i++ ) {
      const int *e = edges[i];
      const int *o = e_offs[i];
      for ( X = 1; X < _width - 1; X++ ) {
        a  = ArgIm.get ( X, e[0],  e[1] ) * 6;

        b = -2 * ( ArgIm.get ( X, e[0] + o[0], e[1] ) + ArgIm.get ( X, e[0], e[1] + o[1] ) );
        b -= ArgIm.get ( X + 1, e[0],  e[1] );
        b -= ArgIm.get ( X - 1, e[0],  e[1] );

        DestIm.set ( X, e[0], e[1], ArgIm.get ( X, e[0], e[1] ) + _hsqrtau * ( a + b ) );
      }
    }

    // traverse the 4 edges parallel to y-axis;
    for ( i = 0; i < 4; i++ ) {
      const int *e = edges[i];
      const int *o = e_offs[i];
      for ( Y = 1; Y < _width - 1; Y++ ) {
        a  = ArgIm.get ( e[0], Y,  e[1] ) * 6;

        b = -2 * ( ArgIm.get ( e[0] + o[0], Y, e[1] ) + ArgIm.get ( e[0], Y, e[1] + o[1] ) );
        b -= ArgIm.get ( e[0], Y + 1,  e[1] );
        b -= ArgIm.get ( e[0], Y - 1,  e[1] );

        DestIm.set ( e[0], Y, e[1], ArgIm.get ( e[0], Y, e[1] ) + _hsqrtau * ( a + b ) );
      }
    }

    // traverse the 4 edges parallel to z-axis;
    for ( i = 0; i < 4; i++ ) {
      const int *e = edges[i];
      const int *o = e_offs[i];
      for ( Z = 1; Z < _width - 1; Z++ ) {
        a  = ArgIm.get ( e[0], e[1],  Z ) * 6;

        b = -2 * ( ArgIm.get ( e[0] + o[0], e[1], Z ) + ArgIm.get ( e[0], e[1] + o[1], Z ) );
        b -= ArgIm.get ( e[0], e[1],  Z + 1 );
        b -= ArgIm.get ( e[0], e[1],  Z - 1 );

        DestIm.set ( e[0], e[1], Z, ArgIm.get ( e[0], e[1], Z ) + _hsqrtau * ( a + b ) );
      }
    }

    const int faces[2] = {
                           0, w
                         };
    const int f_offs[2] = {
                            1, -1
                          };

    // walk on yz-aligned boundary face
    for ( i = 0; i < 2; i++ ) {
      const int f = faces[i];
      const int o = f_offs[i];
      for ( Y = 1; Y < _width - 1; Y++ ) {
        for ( Z = 1; Z < _width - 1; Z++ ) {
          a  = ArgIm.get ( f, Y,  Z ) * 6;

          b = -2 * ArgIm.get ( f + o, Y,  Z );
          b -= ArgIm.get ( f, Y + 1, Z ) + ArgIm.get ( f, Y - 1, Z );
          b -= ArgIm.get ( f, Y, Z + 1 ) + ArgIm.get ( f, Y, Z - 1 );

          DestIm.set ( f, Y, Z, ArgIm.get ( f, Y, Z ) + _hsqrtau * ( a + b ) );
        }
      }
    }

    // walk on xz-aligned boundary face
    for ( i = 0; i < 2; i++ ) {
      const int f = faces[i];
      const int o = f_offs[i];
      for ( X = 1; X < _width - 1; X++ ) {
        for ( Z = 1; Z < _width - 1; Z++ ) {
          a  = ArgIm.get ( X, f,  Z ) * 6;

          b = -2 * ArgIm.get ( X, f + o,  Z );
          b -= ArgIm.get ( X + 1, f, Z ) + ArgIm.get ( X - 1, f, Z );
          b -= ArgIm.get ( X, f, Z + 1 ) + ArgIm.get ( X, f, Z - 1 );

          DestIm.set ( X, f, Z, ArgIm.get ( X, f, Z ) + _hsqrtau * ( a + b ) );
        }
      }
    }

    // walk on xy-aligned boundary face
    for ( i = 0; i < 2; i++ ) {
      const int f = faces[i];
      const int o = f_offs[i];
      for ( X = 1; X < _width - 1; X++ ) {
        for ( Y = 1; Y < _width - 1; Y++ ) {
          a  = ArgIm.get ( X, Y,  f ) * 6;

          b = -2 * ArgIm.get ( X, Y,  f + o );
          b -= ArgIm.get ( X + 1, Y, f ) + ArgIm.get ( X - 1, Y, f );
          b -= ArgIm.get ( X, Y + 1, f ) + ArgIm.get ( X, Y - 1, f );

          DestIm.set ( X, Y, f, ArgIm.get ( X, Y, f ) + _hsqrtau * ( a + b ) );
        }
      }
    }
  }
};


/** Abstract class implementing a very crude but fast version of multigrid.
 *  Not for general use, only used in the linearSmoothOp
 *  \author Droske
 */
template <typename GridTrait>
class MDMultigrid {
public:
  typedef typename GridTrait::RealType RealType;
private:
  typedef typename GridTrait::GridType GridType;
  typedef typename GridTrait::ProlongOpType ProlongOpType;
  typedef typename GridTrait::RestrictOpType RestrictOpType;
  typedef typename GridTrait::MultilevelArrayType MultilevelArrayType;
public:
  MDMultigrid ( const GridType &Grid, int pre = 2, int post = 2 ) :
      _depth ( Grid.getGridDepth() ), _x ( NULL ), _rhs ( NULL ), _dummy ( Grid ), _residuum ( Grid ) {
    preSmoothSteps = pre;
    postSmoothSteps = post;
  }

  /** Delete the MultilevelArrays which have been allocated
   */
  virtual ~MDMultigrid() {}

  /** Run the multigrid.
   */
  void solve ( MultilevelArrayType &X, MultilevelArrayType &Rhs ) const {

    _x = &X;
    _rhs = &Rhs;
    solve ( _x->getDepth() );
  }

protected:
  /** Run the multigrid
   */
  inline void solve ( int level ) const;

  virtual void applyOperator ( const GridType &Grid, const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const = 0;

  /** Interpolate from one a coarse grid to a fine grid
   *  @param Coarse The Vector containing the coarse data
   *  @param Fine The vector in which the prolongated fine data is stored
   *  @param CoarseGrid The definition of the coarse grid
   *  @param FineGrid The definition of the fine grid
   */
  virtual void mgProlongate ( const aol::Vector<RealType> &Coarse,
                              aol::Vector<RealType> &Fine,
                              const GridType &CoarseGrid,
                              const GridType &FineGrid ) const {
    ProlongOpType p ( CoarseGrid, FineGrid );
    p.apply ( Coarse, Fine );
  }

  /** Restrict from one a fine grid to a coarse grid
   *  @param Coarse The Vector containing the coarse data
   *  @param Fine The vector in which the prolongated fine data is stored
   *  @param CoarseGrid The definition of the coarse grid
   *  @param FineGrid The definition of the fine grid
   */

  virtual void mgRestrict ( aol::Vector<RealType> &Coarse,
                            const aol::Vector<RealType> &Fine,
                            const GridType &CoarseGrid,
                            const GridType &FineGrid ) const {
    RestrictOpType r ( CoarseGrid, FineGrid );
    r.apply ( Fine, Coarse );
  }


  /** The smoother. Runs a few cycles of the iterative solver on the Grid
   *  for the system of equations determined by X and RHS.
   */
  virtual void jacobi ( const GridType &Grid,
                        aol::Vector<RealType> &X,
                        const aol::Vector<RealType> &RHS ) const = 0;

  int _depth;            //!< Depth of the grid hierarchy
  mutable MultilevelArrayType *_x;         //!< The hierarchy of solution vectors
  mutable MultilevelArrayType *_rhs;       //!< The hierarchy of right hand sides
  mutable MultilevelArrayType _dummy;      //!< A hierarchy of dummy vectors
  mutable MultilevelArrayType _residuum;   //!< The hierarchy of residual vectors

  int preSmoothSteps;    //!< Number of presmoothig steps of the smoother
  int postSmoothSteps;   //!< Number of postsmoothing steps of the smoother
};


template <typename GridTrait>
void MDMultigrid<GridTrait>::solve ( int level ) const {
  int i;
  GridType grid ( level, _x->getDimension() );
  GridType coarse_grid ( level - 1, _x->getDimension() );

  // If on coarses level, call smoother and return afterwards
  if ( level == 0 ) {
    for ( i = 0; i < 100; i++ ) {
      jacobi ( grid, ( *_x ) [level], ( *_rhs ) [level] );
    }
    return;
  }

  // Step 1: Presmoothing
  for ( i = 0; i < preSmoothSteps; i++ ) {
    jacobi ( grid, ( *_x ) [level], ( *_rhs ) [level] );
  }

  // Step 2: Compute residuum
  applyOperator ( grid, ( *_x ) [level], ( _dummy ) [level] );

  const int size = ( *_x ) [level].size();
  for ( int l = 0; l < size; l++ ) {
    ( _residuum ) [level][l] = ( *_rhs ) [level][l] - ( _dummy ) [level][l];
  }

  // Step 3: Restrict the residuum
  mgRestrict ( ( ( *_rhs ) ) [level-1], ( _residuum ) [level], coarse_grid, grid );

  // Step 4: Ascent one level
  ( ( *_x ) ) [level-1].setZero();
  solve ( level - 1 );

  // Step 5: Add correction
  mgProlongate ( ( ( *_x ) ) [level-1], ( _dummy ) [level], coarse_grid, grid );

  ( ( *_x ) ) [level] += ( _dummy ) [level];

  // Step 6: Postsmoothing
  for ( i = 0; i < postSmoothSteps; i++ ) {
    jacobi ( grid, ( *_x ) [level], ( *_rhs ) [level] );
  }
}


/* ********************************************************************** */

/** Class implementing a multigrid algorithm for the heat-equation in 2D and
 *  3D. Very crude but fast. Not for general use, only for use in linearSmoothOp.h
 *  \author Droske
 */

template <typename GridTrait>
class HeatMultigrid : public MDMultigrid<GridTrait> {
  typedef typename GridTrait::RealType RealType;
  typedef typename GridTrait::GridType GridType;
public:
  /** Constructor which initializes the multigrid for the hea-equation to
   *  run on the given grid.
   */
  HeatMultigrid ( const GridType &Grid )
      : MDMultigrid<GridTrait> ( Grid ) {
    _tau = 0.01;
  }

  void setTau ( RealType Tau ) {
    _tau = Tau;
  }

protected:
  /** Thie method defines how the matrix of the system of equations for the
   *  heat equation is working. In particular this method computes
   *  \f$ Ax = b\f$, where \f$x\f$ is the argument vector and \f$b\f$ is
   *  the destination vector
   *  @param Grid The grid on which the data is defined
   *  @param Arg The argument vector
   *  @param Dest The destination vector
   */

  void applyOperator ( const GridType &Grid,
                       const aol::Vector<RealType> &Arg,
                       aol::Vector<RealType> &Dest ) const {
    switch ( Grid.getDimOfWorld() ) {
    case qc::QC_2D: {
        HeatEquation2DFD<RealType> he ( Grid, _tau );
        he.apply ( Arg, Dest );
        break;
      }
    case qc::QC_3D: {
        HeatEquation3DFD<RealType> he ( Grid, _tau );
        he.apply ( Arg, Dest );
        break;
      }
    default:
      throw aol::Exception ( "invalid dimension", __FILE__, __LINE__ );
    }
  }

  /** Jacobi smoother for the system of equations \f$ A\;X = \rm RHS \f$.
   *  @param Grid The definition of the grid on which the data is defined
   *  @param X The solution vector of the linear system
   *  @param RHS The right-hand-side vector of the linear system
   */

  void jacobi ( const GridType &Grid,
                aol::Vector<RealType> &X,
                const aol::Vector<RealType> &RHS ) const {
    RealType diag;

    if ( Grid.getDimOfWorld() == qc::QC_2D ) {
      diag = 1 + 4 * _tau * aol::Sqr ( Grid.getWidth() - 1 );
    } else if ( Grid.getDimOfWorld() == qc::QC_3D ) {
      diag = 1 + 6 * _tau * aol::Sqr ( Grid.getWidth() - 1 );
    } else {
      throw aol::Exception ( "invalid dimension", __FILE__, __LINE__ );
    }

    int level = Grid.getGridDepth();
    const int size = Grid.getNumberOfNodes();

    RealType relax = 1.0;

    applyOperator ( Grid, X, ( this->_dummy ) [level] );

    for ( int l = 0; l < size; l++ ) {
      RealType v = RHS.get ( l ) - ( this->_dummy ) [level].get ( l );
      X[l] += relax * v / diag;
    }
  }

  RealType _tau;
};


} // end nameless namespace


namespace qc {

template< typename RealType, typename GridTrait>
void LinearSmoothOp<RealType, GridTrait>::apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
  if ( !(_grid.get()) ) {
    throw aol::Exception ( "Set grid prior to calling apply for LinearSmoothOp", __FILE__, __LINE__ );
  }

  if ( Arg.size() != _grid->getNumberOfNodes() || Dest.size() != _grid->getNumberOfNodes() ) {
    throw aol::Exception ( "LinearSmoothOp: Size of Arg or Dest incompatible to gridsize.", __FILE__, __LINE__ );
  }

  const qc::Array<RealType> ArgIm ( Arg, *_grid );
  qc::Array<RealType> DestIm ( Dest, *_grid );

  _x->current().setZero();
  _rhs->current() = ArgIm;

  HeatMultigrid<GridTrait> mg ( *_grid );
  mg.setTau ( _tau );
  mg.solve ( *_x, *_rhs );

  DestIm = _x->current();
}

template class LinearSmoothOp<float, qc::DyadicGridTrait<float> >;
template class LinearSmoothOp<double, qc::DyadicGridTrait<double> >;
template class LinearSmoothOp<long double, qc::DyadicGridTrait<long double> >;
template class LinearSmoothOp<float, qc::CellCenteredGridTrait<float, qc::QC_2D> >;
template class LinearSmoothOp<double, qc::CellCenteredGridTrait<double, qc::QC_2D> >;
template class LinearSmoothOp<long double, qc::CellCenteredGridTrait<long double, qc::QC_2D> >;
template class LinearSmoothOp<float, qc::CellCenteredGridTrait<float, qc::QC_3D> >;
template class LinearSmoothOp<double, qc::CellCenteredGridTrait<double, qc::QC_3D> >;
template class LinearSmoothOp<long double, qc::CellCenteredGridTrait<long double, qc::QC_3D> >;

} // end namespace qc
