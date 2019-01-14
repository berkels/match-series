#ifndef __ISOLINEITERATOR2D_H
#define __ISOLINEITERATOR2D_H

#include <quoc.h>
#include <scalarArray.h>
#include <geom.h>

namespace qc {

/** Represents the class which is returned by the isolineiterator.
 * \author Droske
 */
template <typename RealType, const int d = 2>
class IsoFragmentInf {
public:
  IsoFragmentInf()
      : segment ( points[0], points[1] ) {}

  void makePointsFromWeights() {
    for ( int i = 0; i < d; i++ ) {
      for ( int j = 0; j < d; j++ ) {
        points[i][j] = weights[i] * ( static_cast< RealType > ( intersect[i][1][j] - intersect[i][0][j] ) ) + static_cast<RealType> ( intersect[i][0][j] );
      }
    }
  }

  IsoFragmentInf<RealType, d> &operator= ( const qc::IsoFragmentInf<RealType, d> &Frag ) {
    el = Frag.el;
    for ( int i = 0; i < d; i++ ) {
      intersect[i][0] = Frag.intersect[i][0];
      intersect[i][1] = Frag.intersect[i][1];
      weights[i] = Frag.weights[i];
      points[i] = Frag.points[i];
    }
    return *this;
  }

  // all following coordinates are global and run in each direction from zero up to the grid width (e.g. up to 257 for a 257³ grid)
  // the element which the line (the IsoFragment) passes through
  Element el;
  // the d sides of the element el (each described by start and end point) which the line intersects
  CoordType intersect[d][2];
  // the fraction of the line from intersect[d][0] to intersect[d][1], where the line intersects
  RealType weights[d];
  // intersection points in the reference element
  aol::Vec3<RealType> points[d];
  aol::LineSegment<RealType,qc::QC_3D> segment;

};


template <typename ConfiguratorType> class IsoLineManager2d;

/** Iterator traversing through all linesegments of a discrete isoline
 * of a discrete function represented by a Array
 * \author Droske
 */
template <typename ConfiguratorType>
class IsoLineIterator2d {
  friend class IsoLineManager2d<ConfiguratorType>;
protected:
  const ConfiguratorType _config;
  const typename ConfiguratorType::InitType &_grid;
  typedef typename ConfiguratorType::RealType RealType;
  const Array<RealType> &_func;
  RealType _isoValue;
  bool _reached_end;
  typename ConfiguratorType::ElementIteratorType _elit;
  vector<IsoFragmentInf<RealType, 2> > _curFragments;
  bool _printAmbiguousCaseWarning;

public:

  IsoLineIterator2d ( const typename ConfiguratorType::InitType &Grid, const qc::Array<RealType> &Func )
      : _config ( Grid ), _grid ( Grid ), _func ( Func ), _isoValue ( 0. ), _printAmbiguousCaseWarning ( true ) {
    _curFragments.resize ( 1 );
  }

  IsoLineIterator2d ( qc::IsoLineIterator2d<ConfiguratorType> &It, const bool PrintAmbiguousCaseWarning = true )
      : _config ( It._grid ), _grid ( It._grid ), _func ( It._func ), _isoValue ( It._isoValue ), _reached_end ( It._reached_end ), _printAmbiguousCaseWarning ( PrintAmbiguousCaseWarning ) {
    _elit = It._elit;
    _curFragments = It._curFragments;
  }

  void makeWeights ( IsoFragmentInf<RealType, 2> &FragInf ) {
    for ( int i = 0; i < 2; i++ ) {

      // doch aus versehen falsches vorzeichen?
      if ( ( _isoValue >= _func.get ( FragInf.intersect[i][0] ) ) == ( _isoValue >= _func.get ( FragInf.intersect[i][1] ) ) ) {
        cout << "Error!\n";
      }

      FragInf.weights[i] = ( _func.get ( FragInf.intersect[i][0] ) - _isoValue )
                           / ( _func.get ( FragInf.intersect[i][0] ) - _func.get ( FragInf.intersect[i][1] ) );

      if ( FragInf.weights[i] < 0. || FragInf.weights[i] > 1. ) {
        cout << "weight " << FragInf.weights[i] << " out of bounds.\n";
      }

    }
  }

  bool operator!= ( const IsoLineIterator2d<ConfiguratorType> &It ) {
    return ( It._reached_end != _reached_end );
  }

  bool operator== ( const IsoLineIterator2d<ConfiguratorType> &It ) {
    return ( It._reached_end == _reached_end );
  }

  IsoLineIterator2d<ConfiguratorType>& operator++ ( int ) {
    findNextFragment();
    return *this;
  }

  IsoFragmentInf<RealType, 2> &operator*() {
    return _curFragments[0];
  }

  IsoFragmentInf<RealType, 2>* operator->() {
    return &_curFragments[0];
  }

  void setPrintAmbiguousCaseWarning ( const bool PrintAmbiguousCaseWarning ) {
    _printAmbiguousCaseWarning = PrintAmbiguousCaseWarning;
  }

protected:
  void findNextFragment() {
    // Berkels: ATTENTION!!! For some reason the original author
    // of this class decided to store references to _curFragments[0].
    // YOU MAY NOT DO ANYTHING WHICH BREAKS THIS REFERENCE! Unfortunately
    // the original author didn't take this into account: By using
    // the push_back in initAmbiguousCase the reference is changed!
    // To prevent this I call reserve here.
    _curFragments.reserve ( 3 );
    if ( _curFragments.size() > 1 ) {
      _curFragments[0] = _curFragments.back();
      _curFragments.pop_back();
    } else {
      IsoFragmentInf<RealType, 2> &frag = _curFragments[0];
      for ( ;_elit != _config.end(); ++_elit ) {
        short x = _elit->x(), y = _elit->y();
        int numpos = 0;
        if ( _func.get ( x, y ) > _isoValue ) numpos++;
        if ( _func.get ( x + 1, y ) > _isoValue ) numpos++;
        if ( _func.get ( x + 1, y + 1 ) > _isoValue ) numpos++;
        if ( _func.get ( x, y + 1 ) > _isoValue ) numpos++;

        if ( numpos == 1 || numpos == 3 ) {
          frag.el = *_elit;
          if ( numpos == 1 ) {
            if ( _func.get ( x, y ) > _isoValue ) {
              frag.intersect[0][0].set ( x, y, 0 );
              frag.intersect[0][1].set ( x + 1, y, 0 );
              frag.intersect[1][0].set ( x, y, 0 );
              frag.intersect[1][1].set ( x, y + 1, 0 );
            } else if ( _func.get ( x + 1, y ) > _isoValue ) {
              frag.intersect[0][0].set ( x  , y, 0 );
              frag.intersect[0][1].set ( x + 1, y, 0 );
              frag.intersect[1][0].set ( x + 1, y, 0 );
              frag.intersect[1][1].set ( x + 1, y + 1, 0 );
            } else if ( _func.get ( x + 1, y + 1 ) > _isoValue ) {
              frag.intersect[0][0].set ( x  , y + 1, 0 );
              frag.intersect[0][1].set ( x + 1, y + 1, 0 );
              frag.intersect[1][0].set ( x + 1, y, 0 );
              frag.intersect[1][1].set ( x + 1, y + 1, 0 );
            } else {
              frag.intersect[0][0].set ( x, y, 0 );
              frag.intersect[0][1].set ( x, y + 1, 0 );
              frag.intersect[1][0].set ( x, y + 1, 0 );
              frag.intersect[1][1].set ( x + 1, y + 1, 0 );
            }
          } else {
            if ( _func.get ( x, y ) <= _isoValue ) {
              frag.intersect[0][0].set ( x, y, 0 );
              frag.intersect[0][1].set ( x + 1, y, 0 );
              frag.intersect[1][0].set ( x, y, 0 );
              frag.intersect[1][1].set ( x, y + 1, 0 );
            } else if ( _func.get ( x + 1, y ) <= _isoValue ) {
              frag.intersect[0][0].set ( x, y, 0 );
              frag.intersect[0][1].set ( x + 1, y, 0 );
              frag.intersect[1][0].set ( x + 1, y, 0 );
              frag.intersect[1][1].set ( x + 1, y + 1, 0 );
            } else if ( _func.get ( x + 1, y + 1 ) <= _isoValue ) {
              frag.intersect[0][0].set ( x  , y + 1, 0 );
              frag.intersect[0][1].set ( x + 1, y + 1, 0 );
              frag.intersect[1][0].set ( x + 1, y  , 0 );
              frag.intersect[1][1].set ( x + 1, y + 1, 0 );
            } else {
              frag.intersect[0][0].set ( x  , y, 0 );
              frag.intersect[0][1].set ( x  , y + 1, 0 );
              frag.intersect[1][0].set ( x  , y + 1, 0 );
              frag.intersect[1][1].set ( x + 1, y + 1, 0 );
            }
          }
          makeWeights ( frag );
          frag.makePointsFromWeights();
          _elit++;
          return;
        } // endif ( numpos == 1 || numpos == 3 )
        else if ( numpos == 2 ) {
          frag.el = *_elit;
          RealType v11 = _func.get ( x, y );
          RealType v21 = _func.get ( x + 1, y );
          RealType v22 = _func.get ( x + 1, y + 1 );
          if ( ( v11 > _isoValue ) == ( v22 > _isoValue ) ) {
            // the ambiguous case: values of opposite corners have same sign
            if ( _printAmbiguousCaseWarning )
              cerr << "\n\n\n ambiguous case!!!! \n\n\n";
            CoordType coord ( x, y, 0 );
            initAmbiguousCase ( *_elit );
          } else {
            if ( ( v11 > _isoValue ) == ( v21 > _isoValue ) ) {
              frag.intersect[0][0].set ( x, y  , 0 );
              frag.intersect[0][1].set ( x, y + 1, 0 );
              frag.intersect[1][0].set ( x + 1, y, 0 );
              frag.intersect[1][1].set ( x + 1, y + 1, 0 );
            } else {
              frag.intersect[0][0].set ( x, y, 0 );
              frag.intersect[0][1].set ( x + 1, y, 0 );
              frag.intersect[1][0].set ( x, y + 1, 0 );
              frag.intersect[1][1].set ( x + 1, y + 1, 0 );
            }
          }
          makeWeights ( frag );
          frag.makePointsFromWeights();
          _elit++;
          return;
        } // endif numpos == 2
      } // for
      _reached_end = true;
    }
  }

  void initAmbiguousCase ( const Element &el ) {
    const short x = el.x();
    const short y = el.y();
    const RealType v11 = _func.get ( x, y );
    const RealType v21 = _func.get ( x + 1, y );
    const RealType v22 = _func.get ( x + 1, y + 1 );
    const RealType v12 = _func.get ( x, y + 1 );

    RealType len1, len2;
    const RealType y1 = ( _isoValue - v11 ) / ( v12 - v11 );
    const RealType y2 = ( _isoValue - v21 ) / ( v22 - v21 );
    const RealType x1 = ( _isoValue - v11 ) / ( v21 - v11 );
    const RealType x2 = ( _isoValue - v12 ) / ( v22 - v12 );

    // case 1: [\\]
    aol::Vec3<RealType> pt1 ( 0., y1, 0. );
    aol::Vec3<RealType> pt2 ( x1, 0., 0. );
    aol::LineSegment<RealType,qc::QC_3D> segment ( pt1, pt2 );
    len1 = segment.length();

    pt1.set ( x2, 1., 0. );
    pt2.set ( 1., y2, 0. );
    len1 += segment.length();

    // case 2: [//]
    pt1.set ( x1, 0., 0. );
    pt2.set ( 1., y2, 0. );
    len2 = segment.length();

    pt1.set ( 0., y1, 0. );
    pt2.set ( x2, 1., 0. );
    len2 += segment.length();

    IsoFragmentInf<RealType, 2> &frag1 = _curFragments[0];
    ;
    IsoFragmentInf<RealType, 2> frag2;

    frag2.el = el;

    if ( len1 < len2 ) {
      // [\\]
      frag1.intersect[0][0].set ( x, y, 0 );
      frag1.intersect[0][1].set ( x + 1, y, 0 );
      frag1.intersect[1][0].set ( x, y, 0 );
      frag1.intersect[1][1].set ( x, y + 1, 0 );
      frag1.weights[0] = x1;
      frag1.weights[1] = y1;

      frag2.intersect[0][0].set ( x + 1, y, 0 );
      frag2.intersect[0][1].set ( x + 1, y + 1, 0 );
      frag2.intersect[1][0].set ( x, y + 1, 0 );
      frag2.intersect[1][1].set ( x + 1, y + 1, 0 );
      frag2.weights[0] = y2;
      frag2.weights[1] = x2;
    } else {
      // [//]
      frag1.intersect[0][0].set ( x  , y, 0 );
      frag1.intersect[0][1].set ( x + 1, y, 0 );
      frag1.intersect[1][0].set ( x + 1, y, 0 );
      frag1.intersect[1][1].set ( x + 1, y + 1, 0 );
      frag1.weights[0] = x1;
      frag1.weights[1] = y2;

      frag2.intersect[0][0].set ( x, y, 0 );
      frag2.intersect[0][1].set ( x, y + 1, 0 );
      frag2.intersect[1][0].set ( x, y + 1, 0 );
      frag2.intersect[1][1].set ( x + 1, y + 1, 0 );
      frag2.weights[0] = y1;
      frag2.weights[1] = x2;
    }
    frag1.makePointsFromWeights();
    frag2.makePointsFromWeights();
    // Berkels: ATTENTION!!! This could break the reference
    // to _curFragments[0] stored in IsoFragmentInf<RealType, 2> &frag
    // (used in void findNextFragment())
    _curFragments.push_back ( frag2 );
  }


};

/** Class needed by the isolineiterator to specify begin and end iterators
 * and to manage the discrete function as well as the isovalue
 * \author Droske
 */
template <typename ConfiguratorType>
class IsoLineManager2d {
protected:
  const ConfiguratorType _config;
  const typename ConfiguratorType::InitType &_grid;
  typedef typename ConfiguratorType::RealType RealType;
  const Array<RealType> &_func;
  RealType _isoValue;
  IsoLineIterator2d<ConfiguratorType> _begin_it, _end_it;
public:
  IsoLineManager2d ( const typename ConfiguratorType::InitType &Grid,
                     const Array<RealType> &Func )
      : _config( Grid), _grid ( Grid ), _func ( Func ), _isoValue ( 0. ),
      _begin_it ( Grid, Func ), _end_it ( Grid, Func ) {

    _begin_it._elit = _config.begin();
    _end_it._elit = _config.begin();
    _begin_it._isoValue = _end_it._isoValue = _isoValue;
    _begin_it.findNextFragment();
    _end_it._reached_end = true;
    _begin_it._reached_end = false;

  }

  IsoLineIterator2d<ConfiguratorType> &begin() {
    return _begin_it;
  }

  IsoLineIterator2d<ConfiguratorType> &end() {
    return _end_it;
  }

  void setIsoValue ( RealType IsoValue ) {
    _isoValue = IsoValue;
    _begin_it._isoValue = _end_it._isoValue = _isoValue;
    _begin_it._elit = _config.begin();
    _begin_it.findNextFragment();
  }

protected:

};

}

#endif
