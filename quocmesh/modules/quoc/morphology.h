#ifndef __MORPHOLOGY_H
#define __MORPHOLOGY_H

// some basic morphological operators that are handy from time to time
// they are based on upwind-schemes
// WARNING: 3D not tested yet


#include <levelSet.h>

//! Basic morphological operators
namespace morph {

template <typename RealType, qc::Dimension Dim>
class Dilation { };

template <typename RealType>
class Dilation<RealType, qc::QC_2D> : public aol::Op<qc::ScalarArray<RealType, qc::QC_2D> > {
protected:
  const qc::GridDefinition &_grid;
  RealType _tau, _radius;
  const RealType _velocity;
  int _order;

public:
  class Dilate2d : public qc::LevelSet2dInt<RealType, Dilate2d> {
  protected:
    const RealType _v;
  public:
    Dilate2d ( const RealType v = -1. ) : qc::LevelSet2dInt<RealType, Dilate2d> (),  _v ( v ) {}
    RealType velocity ( const qc::ScalarArray<RealType, qc::QC_2D> &, int, int ) const {
      return _v;
    }
  };


  Dilation ( const qc::GridDefinition &grid, int order = 3 /* of the ENO */, RealType velocity = -1. ) :
    _grid ( grid ), _tau ( 0.5 * grid.H() ), _radius ( grid.H() ), _velocity ( velocity ), _order ( order ) {
    if ( grid.getDimOfWorld() != qc::QC_2D ) {
      throw aol::Exception ( "According to the template argument, you should give me a 2d grid, please.", __FILE__, __LINE__ );
    }
  };

  void setRadius ( RealType radius ) {
    _radius = radius;
  }

  void setTau ( RealType tau ) {
    _tau = tau;
  }

  void apply ( const qc::ScalarArray<RealType, qc::QC_2D> &arg, qc::ScalarArray<RealType, qc::QC_2D> &dest ) const {
    Dilate2d dilate ( _velocity );
    dest = arg;
    dilate.setData ( &dest );

    // because speed is 1.
    // let's be generous and allow some round errors in case the fraction is not an integer
    const int steps = static_cast<int> ( _radius / _tau );

    for ( int i = 0; i < steps; i++ ) {
      dilate.timeStepENO ( _tau, _order );
    }
  }

  void applyAdd ( const qc::ScalarArray<RealType, qc::QC_2D> &arg, qc::ScalarArray<RealType, qc::QC_2D> &dest ) const {
    qc::ScalarArray<RealType, qc::QC_2D> tmp ( dest );
    apply ( arg, tmp );
    dest += tmp;
  }
};

template <typename RealType>
class Dilation<RealType, qc::QC_3D> : public aol::Op<qc::ScalarArray<RealType, qc::QC_3D> > {
protected:
  const qc::GridDefinition &_grid;
  RealType _tau, _radius;
  const RealType _velocity;
  int _order;

public:
  class Dilate3d : public qc::LevelSet3dInt<RealType, Dilate3d> {
  protected:
    const RealType _v;
  public:
    Dilate3d ( const RealType v = -1. ) : qc::LevelSet3dInt<RealType, Dilate3d> (),  _v ( v ) {}
    RealType velocity ( const qc::ScalarArray<RealType, qc::QC_3D> &, int, int, int ) const {
      return _v;
    }
  };


  Dilation ( const qc::GridDefinition &grid, int order = 1 /* of the ENO */, RealType velocity = -1. ) :
    _grid ( grid ), _tau ( 0.5 * grid.H() ), _radius ( grid.H() ), _velocity ( velocity ), _order ( order ) {

    cerr << "WARNING: 3D version of dilation has not been tested yet..." << endl;

    if ( grid.getDimOfWorld() != qc::QC_3D ) {
      throw aol::Exception ( "According to the template argument, you should give me a 3d grid, please.", __FILE__, __LINE__ );
    }
    if ( order != 1 ) {
      throw aol::Exception ( "The 3d upwind scheme does not support higher order so far.", __FILE__, __LINE__ );
    }
  };

  void setRadius ( RealType radius ) {
    _radius = radius;
  }

  void setTau ( RealType tau ) {
    _tau = tau;
  }

  void apply ( const qc::ScalarArray<RealType, qc::QC_3D> &arg, qc::ScalarArray<RealType, qc::QC_3D> &dest ) const {
    Dilate3d dilate ( _velocity );
    dest = arg;
    dilate.setData ( &dest );

    // because speed is 1.
    // let's be generous and allow some round errors in case the fraction is not an integer
    const int steps = static_cast<int> ( _radius / _tau );

    for ( int i = 0; i < steps; i++ ) {
      dilate.timeStepEO ( _tau );
    }
  }

  void applyAdd ( const qc::ScalarArray<RealType, qc::QC_3D> &arg, qc::ScalarArray<RealType, qc::QC_3D> &dest ) const {
    qc::ScalarArray<RealType, qc::QC_3D> tmp ( dest );
    apply ( arg, tmp );
    dest += tmp;
  }
};



template <typename RealType, qc::Dimension Dim>
class Erosion : public Dilation<RealType, Dim> {
public:
  Erosion ( const qc::GridDefinition &grid, int order = 3 /* of the ENO */ )
    : Dilation<RealType, Dim> ( grid, order, 1. ) { }
};

template <typename RealType, qc::Dimension Dim>
class Closing : public aol::Op<qc::ScalarArray<RealType, qc::QC_2D> > {
protected:
  const qc::GridDefinition &_grid;
  RealType _tau, _radius;
  int _order;
public:
  Closing ( const qc::GridDefinition &grid, int order = 3 ) : _grid ( grid ), _order ( order ) { }

  void setRadius ( RealType radius ) {
    _radius = radius;
  }

  void setTau ( RealType tau ) {
    _tau = tau;
  }

  void apply ( const qc::ScalarArray<RealType, qc::QC_2D> &arg, qc::ScalarArray<RealType, qc::QC_2D> &dest ) const {
    qc::ScalarArray<RealType, qc::QC_2D> tmp ( dest );

    Dilation<RealType, Dim> dilation ( _grid, _order );
    dilation.setTau ( _tau );
    dilation.setRadius ( _radius );

    dilation.apply ( arg, tmp );

    Erosion<RealType, Dim> erosion ( _grid, _order );
    erosion.setTau ( _tau );
    erosion.setRadius ( _radius );

    erosion.apply ( tmp, dest );
  }

  void applyAdd ( const qc::ScalarArray<RealType, qc::QC_2D> &arg, qc::ScalarArray<RealType, qc::QC_2D> &dest ) const {
    qc::ScalarArray<RealType, qc::QC_2D> tmp ( dest );
    apply ( arg, tmp );
    dest += tmp;
  }
};


template <typename RealType, qc::Dimension Dim>
class Opening : public aol::Op<qc::ScalarArray<RealType, qc::QC_2D> > {
protected:
  const qc::GridDefinition &_grid;
  RealType _tau, _radius;
  const RealType _velocity;
  int _order;
public:
  Opening ( const qc::GridDefinition &grid, int order = 3 ) : _grid ( grid ), _order ( order ) { }

  void setRadius ( RealType radius ) {
    _radius = radius;
  }

  void setTau ( RealType tau ) {
    _tau = tau;
  }

  void apply ( const qc::ScalarArray<RealType, qc::QC_2D> &arg, qc::ScalarArray<RealType, qc::QC_2D> &dest ) const {
    qc::ScalarArray<RealType, qc::QC_2D> tmp ( dest );

    Erosion<RealType, Dim> erosion ( _grid, _order );
    erosion.setTau ( _tau );
    erosion.setRadius ( _radius );

    erosion.apply ( arg, tmp );

    Dilation<RealType, Dim> dilation ( _grid, _order );
    dilation.setTau ( _tau );
    dilation.setRadius ( _radius );

    dilation.apply ( tmp, dest );
  }

  void applyAdd ( const qc::ScalarArray<RealType, qc::QC_2D> &arg, qc::ScalarArray<RealType, qc::QC_2D> &dest ) const {
    qc::ScalarArray<RealType, qc::QC_2D> tmp ( dest );
    apply ( arg, tmp );
    dest += tmp;
  }
};




} // end namespace morph;


#endif
