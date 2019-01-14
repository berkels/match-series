#ifndef __DISCRETEFUNCTION_H
#define __DISCRETEFUNCTION_H

#include <vec.h>
#include <multiVector.h>
#include <smallMat.h>
#include <pointerClasses.h>

namespace aol {

/** Piecewise linear interpolation of discrete data given at not necessarily equidistant sample points
 *  \author Schwen
 */
template < typename DomType, typename DiscValType, typename ContValType = DiscValType >
class DiscreteValueInterpolator {
public:
  typedef std::map < DomType, DiscValType > MapType;

protected:
  MapType _vals;

public:
  // destructor, default copy constructor and assignment operators are OK

  DiscreteValueInterpolator ( ) {
  }

  DiscreteValueInterpolator ( const aol::Vector<DomType> &domPoints, const aol::Vector<DiscValType> &valPoints ) {
    if ( domPoints.size() != valPoints.size() ) {
      throw aol::Exception ( "vector sizes do not match", __FILE__, __LINE__ );
    }

    for ( int i = 0; i < domPoints.size(); ++i ) {
      _vals[ domPoints[i] ] = valPoints[i];
    }
  }

  ContValType evaluate ( const DomType &point ) const {
    if ( point < _vals.begin()->first ) {
      cerr << point << " " << _vals.begin()->first << endl;
      throw aol::Exception ( "Value too small, not in interpolatable range.", __FILE__, __LINE__ );
    } else if ( point > _vals.rbegin()->first ) {
      cerr << point << " " << _vals.rbegin()->first << endl;
      throw aol::Exception ( "Value too large, not in interpolatable range.", __FILE__, __LINE__ );
    }

    const typename MapType::const_iterator found = _vals.find ( point );
    if ( found != _vals.end() ) {
      return ( static_cast< ContValType > ( found->second ) );
    }

    typename MapType::const_iterator ubd1 = _vals.upper_bound ( point ), ubd0 = ubd1; --ubd0;
    const double v = static_cast<double> ( point - ubd1->first ) / ( ubd0->first - ubd1->first );

    if ( ( v < 0 ) || ( v > 1 ) ) {
      cerr << point << " " << ubd0->first << " " << ubd0->second << " " << ubd1->first << " " << ubd1->second << " " << v << endl;
      throw aol::Exception ( "barycentric coordinate outside [0,1]", __FILE__, __LINE__ );
    }

    return ( v * ubd0->second + ( 1.0 - v ) * ubd1->second );
  }

  void setVals ( const MapType &Vals ) {
    _vals = Vals;
  }

  void invertIfMonotonic ( ) {
    aol::Vector<DiscValType> differences;
    for ( typename MapType::const_iterator it = _vals.begin(); it != _vals.end(); ++it ) {
      typename MapType::const_iterator itp = it; ++itp;
      if ( itp != _vals.end() ) {
        differences.pushBack ( itp->second - it->second );
      }
    }
    for ( int i = 0; i < differences.size() - 1; ++i ) {
      if ( aol::signum ( differences[i+1] ) * aol::signum ( differences[i] ) != aol::NumberTrait<DiscValType>::one ) {
        throw aol::Exception ( "cannot invert DiscreteValueInterpolator that is not strictly monotonic", __FILE__, __LINE__ );
      }
    }

    MapType inverted;
    for ( typename MapType::const_iterator it = _vals.begin(); it != _vals.end(); ++it ) {
      inverted[ it->second ] = it->first;
    }

    _vals = inverted;
  }

  void dump ( ostream &out = cout ) const {
    for ( typename MapType::const_iterator it = _vals.begin(); it != _vals.end(); ++it ) {
      out << it->first << " " << it->second << endl;
    }
  }

  void clear ( ) {
    _vals.clear();
  }

  MapType& getInternalMapRef ( ) {
    return ( _vals );
  }

  const MapType& getInternalMapRef ( ) const {
    return ( _vals );
  }
};


/**
 * interface class for discrete functions (FE-functions)
 * \author Droske
 */
template <typename ConfiguratorType, typename Imp>
class DiscreteFunctionInterface {

public:
  typedef typename ConfiguratorType::RealType RealType;

protected:
  DeleteFlagPointer<const ConfiguratorType> _configPtr;
  const typename ConfiguratorType::InitType &_initializer;

  // barton-nackman
  typedef typename ConfiguratorType::ElementType ElementType;
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

public:
  //! constructor from grid ( = initializer)
  DiscreteFunctionInterface ( const typename ConfiguratorType::InitType &Initializer ) :
      _configPtr ( new ConfiguratorType(Initializer), /* deleteFlag = */ true ),
      _initializer ( Initializer ) {}

  //! constructor from configurator
  DiscreteFunctionInterface ( const ConfiguratorType & config ) :
      _configPtr ( &config, /* deleteFlag = */ false ),
      _initializer ( config.getInitializer() ) {}

  //! copy constructor
  //! We have to make sure that the copy survives, after the copied object is destroyed.
  DiscreteFunctionInterface ( const DiscreteFunctionInterface<ConfiguratorType, Imp> & other ) :
      _configPtr ( other._configPtr.getDeleteFlag() ? new ConfiguratorType( *(other._configPtr.get()) ) : other._configPtr.get(),
                   /* deleteFlag = */ other._configPtr.getDeleteFlag() ),
      _initializer ( other._initializer ) {}

  typedef typename ConfiguratorType::VecType  VecType;
  typedef typename ConfiguratorType::DomVecType  DomVecType;

  const ConfiguratorType& getConfigurator( ) const { return *_configPtr; }

  /*************************************************************************
   ***** interface begins here >>> *****************************************
   *************************************************************************/
  RealType evaluateAtQuadPoint ( const ElementType &El, int QuadPoint ) const {
    return asImp().evaluateAtQuadPoint ( El, QuadPoint );
  }

  RealType evaluate ( const ElementType &El, const DomVecType& RefCoord ) const {
    return asImp().evaluate ( El, RefCoord );
  }

  void evaluateGradientAtQuadPoint ( const ElementType &El, int QuadPoint, VecType& Grad ) const {
    asImp().evaluateGradientAtQuadPoint ( El, QuadPoint, Grad );
  }

  void evaluateGradient ( const ElementType &El, const DomVecType& RefCoord, VecType& Grad ) const {
    asImp().evaluateGradient ( El, RefCoord, Grad );
  }

#if defined ( USE_CPP11 )
  template < typename MatType >
  void evaluateSecondDerivativeAtQuadPoint ( const ElementType &El, int QuadPoint, MatType& Hessian ) const {
    asImp().evaluateSecondDerivativeAtQuadPoint ( El, QuadPoint, Hessian );
  }

  template < typename MatType >
  void evaluateSecondDerivative ( const ElementType &El, const DomVecType& RefCoord, MatType& Hessian ) const {
    asImp().evaluateSecondDerivative ( El, RefCoord, Hessian );
  }
#endif

  /*************************************************************************
   ***** <<< interface ends here *******************************************
   *************************************************************************/

  RealType evaluateGlobal ( const DomVecType& Coord ) const {
    DomVecType localCoord;
    ElementType el = getConfigurator().getEmptyElement();
    getConfigurator().getLocalCoords ( Coord, el, localCoord );

    return evaluate ( el, localCoord );
  }

  void evaluateGradientGlobal ( const DomVecType& Coord, VecType& Grad ) const {
    DomVecType localCoord;
    ElementType el = getConfigurator().getEmptyElement();
    getConfigurator().getLocalCoords ( Coord, el, localCoord );
    evaluateGradient ( el, localCoord, Grad );
  }
};

/**
 * special implementation for discrete functions, that make use of the basefunctionset
 * \author Droske, Effland
 */
template <typename ConfiguratorType, typename DofType = aol::Vector < typename ConfiguratorType::RealType > >
class DiscreteFunctionDefault : public DiscreteFunctionInterface<ConfiguratorType, DiscreteFunctionDefault<ConfiguratorType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VecType  VecType;
  typedef typename ConfiguratorType::DomVecType  DomVecType;
  typedef typename ConfiguratorType::ElementType ElementType;

  const DofType & _dofs;

  typedef ConfiguratorType ConfType;

  //! constructor from grid ( = initializer)
  DiscreteFunctionDefault ( const typename ConfiguratorType::InitType &Initializer, const DofType &Dofs )
      : DiscreteFunctionInterface<ConfiguratorType, DiscreteFunctionDefault<ConfiguratorType> > ( Initializer ), _dofs ( Dofs ) {}

  //! constructor from configurator
  DiscreteFunctionDefault ( const ConfiguratorType & config, const DofType & Dofs )
    : DiscreteFunctionInterface<ConfiguratorType, DiscreteFunctionDefault<ConfiguratorType> > ( config ), _dofs ( Dofs ) {}

  // implicit copy constructor does the right thing.
  // DiscreteFunctionDefault ( const DiscreteFunctionDefault<ConfiguratorType> &DiscrFunc );

  RealType evaluate ( const ElementType &El, const DomVecType& RefCoord ) const {
    const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( El );
    RealType w = 0.;
    for ( int b = 0; b < this->getConfigurator().getNumLocalDofs ( El ); b++ ) {
      w += this->_dofs[ this->getConfigurator().localToGlobal ( El, b ) ] * bfs.evaluate ( b, RefCoord );
    }
    return w;
  }

  RealType evaluateAtQuadPoint ( const ElementType &El, int QuadPoint ) const {
    const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( El );
    RealType w = 0.;
    for ( int b = 0; b < static_cast<int> ( this->getConfigurator().getNumLocalDofs ( El ) ); b++ ) {
      w += this->_dofs[ this->getConfigurator().localToGlobal ( El, b ) ] * bfs.evaluate ( b, QuadPoint );
    }
    return w;
  }

  void evaluateGradientAtQuadPoint ( const ElementType &El, int QuadPoint, aol::Vec<ConfiguratorType::Dim, RealType>& Grad ) const {
    const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( El );
    VecType v;
    Grad.setZero();
    for ( int b = 0; b < this->getConfigurator().getNumLocalDofs ( El ); b++ ) {
      v = bfs.evaluateGradient ( b, QuadPoint );
      v *= this->_dofs[ this->getConfigurator().localToGlobal ( El, b ) ];
      Grad += v;
    }
  }

  void evaluateGradient ( const ElementType &El, const DomVecType& RefCoord, aol::Vec<ConfiguratorType::Dim, RealType>& Grad ) const {
    const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( El );
    VecType v;
    Grad.setZero();

    for ( int b = 0; b < this->getConfigurator().getNumLocalDofs ( El ); b++ ) {
      bfs.evaluateGradient ( b, RefCoord, v );
      v *= this->_dofs[ this->getConfigurator().localToGlobal ( El, b ) ];
      Grad += v;
    }
  }

  // note: not every base function set allows the evaluation of second derivatives
  // in these cases, calling this function results in a compilation error
  template < typename MatType >
  void evaluateSecondDerivativeAtQuadPoint ( const ElementType &El, int QuadPoint, MatType& Hessian ) const {
    const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( El );
    MatType v;
    Hessian.setZero();
    for ( int b = 0; b < this->getConfigurator().getNumLocalDofs ( El ); ++b ) {
      v = bfs.evaluateSecondDerivative ( b, QuadPoint );
      v *= _dofs[this->getConfigurator().localToGlobal ( El, b )];
      Hessian += v;
    }
  }

  // note: not every base function set allows the evaluation of second derivatives
  // in these cases, calling this function results in a compilation error
  template < typename MatType >
  void evaluateSecondDerivative ( const ElementType &El, const DomVecType& RefCoord, MatType& Hessian ) const {
    const typename ConfiguratorType::BaseFuncSetType &bfs = this->getConfigurator().getBaseFunctionSet ( El );
    MatType v;
    Hessian.setZero();
    for ( int b = 0; b < this->getConfigurator().getNumLocalDofs ( El ); ++b ) {
      bfs.evaluateSecondDerivative ( b, RefCoord, v );
      v *= _dofs[this->getConfigurator().localToGlobal ( El, b )];
      Hessian += v;
    }
  }

  const DofType & getDofs() const {
    return _dofs;
  }

  ostream& print ( ostream &os ) const {
    bool oldPrettyFormat = _dofs.getPrettyFormat();
    _dofs.setPrettyFormat ( false );
    os << "Dimension #Dofs: " << ConfiguratorType::Dim << ' ' << _dofs.size() << '\n';
    os << _dofs;
    _dofs.setPrettyFormat ( oldPrettyFormat );
    return os;
  }
};

template < typename ConfiguratorType >
ostream& operator<< ( ostream &os, DiscreteFunctionDefault < ConfiguratorType > &dfd ) {
  return dfd.print ( os );
}

template < typename ConfiguratorType, typename DofType >
inline void readDiscreteFunctionDofsFromFile ( const ConfiguratorType &configurator, string &filename, DofType &dofs ) {
  readDiscreteFunctionDofsFromFile ( configurator.Dim, filename, dofs );
}

template < typename DofType >
void readDiscreteFunctionDofsFromFile ( const qc::Dimension dim, string &filename, DofType &dofs ) {
  ifstream is ( filename.c_str () );

  char input[17];
  int dimension;
  int numDofs;
  is.read ( input, 17 );
  if ( strncmp ( input, "Dimension #Dofs: ", 17 ) != 0 ) {
    throw aol::Exception ( "DiscreteFunctionDefault: Wrong input format!", __FILE__, __LINE__ );
  }

  is >> dimension;
  if ( dimension != dim )
    throw aol::Exception ( "DiscreteFunctionDefault: Wrong dimension!", __FILE__, __LINE__ );

  is >> numDofs;

  dofs.resize ( numDofs );

  is >> dofs;

  is.close ();
}

/**
 * special implementation for discrete functions, that make use of the basefunctionset
 * \author Droske
 */
template <typename ConfiguratorType, typename DiscreteFunctionType, typename Imp>
class ComposedDiscreteFunctionInterface : public DiscreteFunctionInterface < ConfiguratorType,
      ComposedDiscreteFunctionInterface<ConfiguratorType, DiscreteFunctionType, Imp> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VecType  VecType;
  typedef typename ConfiguratorType::DomVecType  DomVecType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef ConfiguratorType ConfType;

  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

  const DiscreteFunctionType &_discreteFunction;

  ComposedDiscreteFunctionInterface ( const typename ConfiguratorType::InitType &Initializer, const DiscreteFunctionType &DiscreteFunction )
      : DiscreteFunctionInterface<ConfiguratorType, ComposedDiscreteFunctionInterface<ConfiguratorType, DiscreteFunctionType, Imp> > ( Initializer ),
      _discreteFunction ( DiscreteFunction ) {
	throw aol::Exception ( "Check whether gradient aol::Vec<> in evaluateGradient[AtQuadPoint] has correct dimension and remove this exception!", __FILE__, __LINE__ );
      }

  //! copy constructor
  ComposedDiscreteFunctionInterface ( const ComposedDiscreteFunctionInterface<ConfiguratorType, DiscreteFunctionType, Imp> &DiscrFunc )
      : DiscreteFunctionInterface<ConfiguratorType, DiscreteFunctionDefault<ConfiguratorType> > ( DiscrFunc._initializer ),
      _discreteFunction ( DiscrFunc._discreteFunction ) {
	throw aol::Exception ( "Check whether gradient aol::Vec<> in evaluateGradient[AtQuadPoint] has correct dimension and remove this exception!", __FILE__, __LINE__ );
      }

  RealType evaluateCompositionFunction ( RealType s ) const {
    return asImp().evaluateCompositionFunction ( s );
  }

  RealType evaluateCompositionFunctionDerivative ( RealType s ) const {
    return asImp().evaluateCompositionFunctionDerivative ( s );
  }

  RealType evaluateGlobal ( const DomVecType& Coord ) const {
    return evaluateCompositionFunction ( _discreteFunction.evaluateGlobal ( Coord ) );
  }

  RealType evaluate ( const ElementType &El, const DomVecType& RefCoord ) const {
    return evaluateCompositionFunction ( _discreteFunction.evaluate ( El, RefCoord ) );
  }

  RealType evaluateAtQuadPoint ( const ElementType &El, int QuadPoint ) const {
    return evaluateCompositionFunction ( _discreteFunction.evaluate ( El, QuadPoint ) );
  }

  void evaluateGradientAtQuadPoint ( const ElementType &El, int QuadPoint, aol::Vec<ConfiguratorType::DomDim, RealType>& Grad ) const {
    _discreteFunction.evaluateGradientAtQuadPointWeight ( El, QuadPoint, Grad );
    Grad *= evaluateCompositionFunctionDerivative ( _discreteFunction.evaluateAtQuadPoint ( El, QuadPoint ) );
  }

  void evaluateGradient ( const ElementType &El, const DomVecType& RefCoord, aol::Vec<ConfiguratorType::DomDim, RealType>& Grad ) const {
    _discreteFunction.evaluateGradient ( El, RefCoord, Grad );
    Grad *= evaluateCompositionFunctionDerivative ( _discreteFunction.evaluate ( El, RefCoord ) );
  }
};

template <typename ConfiguratorType, typename DiscFuncType, int DimRange>
class DiscreteVectorFunction {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VecType  VecType;
  typedef typename ConfiguratorType::DomVecType  DomVecType;
  typedef typename ConfiguratorType::ElementType ElementType;

  mutable aol::auto_container<DimRange, DiscFuncType> _discrFuncs;
  const aol::MultiVector<RealType> & _dofsReference;
protected:
public:
  //! constructor from grid
  DiscreteVectorFunction ( const typename ConfiguratorType::InitType &Initializer,
                           const aol::MultiVector<RealType> &Dofs,
                           bool allowMoreComponentsThanNeeded ) :
      _dofsReference ( Dofs ) {
    if ( !allowMoreComponentsThanNeeded ) {
      if ( Dofs.numComponents() != DimRange )
        throw aol::Exception ( "you have to pass a multivector with the correct amount of components.", __FILE__, __LINE__ );
    } else {
      if ( Dofs.numComponents() < DimRange )
        throw aol::Exception ( "you have to pass a multivector with at least as much components as the discrete function has.", __FILE__, __LINE__ );
    }
    for ( int c = 0; c < DimRange; c++ ) {
      _discrFuncs.take_over ( c, new DiscFuncType ( Initializer, Dofs[c] ) );
    }
  }

  //! constructor from configurator
  DiscreteVectorFunction ( const ConfiguratorType &Configurator,
                           const aol::MultiVector<RealType> &Dofs,
                           bool allowMoreComponentsThanNeeded ) :
      _dofsReference ( Dofs ) {
    if ( !allowMoreComponentsThanNeeded ) {
      if ( Dofs.numComponents() != DimRange )
        throw aol::Exception ( "you have to pass a multivector with the correct amount of components.", __FILE__, __LINE__ );
    } else {
      if ( Dofs.numComponents() < DimRange )
        throw aol::Exception ( "you have to pass a multivector with at least as much components as the discrete function has.", __FILE__, __LINE__ );
    }
    for ( int c = 0; c < DimRange; c++ ) {
      _discrFuncs.take_over ( c, new DiscFuncType ( Configurator, Dofs[c] ) );
    }
  }


  void evaluate ( const ElementType &El, const DomVecType& RefCoord, aol::Vec<DimRange, RealType> &Value ) const {
    for ( int c = 0; c < DimRange; c++ ) {
      Value[c] = _discrFuncs[c].evaluate ( El, RefCoord );
    }
  }

  void evaluateAtQuadPoint ( const ElementType &El, int QuadPoint, aol::Vec<DimRange, RealType> &Value ) const {
    for ( int c = 0; c < DimRange; c++ ) {
      Value[c] = _discrFuncs[c].evaluateAtQuadPoint ( El, QuadPoint );
    }
  }

  void evaluateGradient ( const ElementType &El, const DomVecType& RefCoord, Mat<DimRange, ConfiguratorType::Dim, RealType> &Mat ) const {
    VecType v;
    for ( int c = 0; c < DimRange; c++ ) {
      _discrFuncs[c].evaluateGradient ( El, RefCoord, v );
      Mat.setRow ( c, v );
    }
  }

  void evaluateGradientAtQuadPoint ( const ElementType &El, int QuadPoint, Mat<DimRange, ConfiguratorType::Dim, RealType> &Mat ) const {
    VecType v;
    for ( int c = 0; c < DimRange; c++ ) {
      _discrFuncs[c].evaluateGradientAtQuadPoint ( El, QuadPoint, v );
      Mat.setRow ( c, v );
    }
  }


  const DiscFuncType& operator[] ( int i ) const {
    return _discrFuncs[i];
  }

  DiscFuncType& operator[] ( int i ) {
    return _discrFuncs[i];
  }


  const aol::MultiVector<RealType> & getDofs() const {
    return _dofsReference;
  }
};


template <typename ConfiguratorType, int DimRange>
class DiscreteVectorFunctionDefault : public DiscreteVectorFunction<ConfiguratorType, DiscreteFunctionDefault<ConfiguratorType>, DimRange > {

public:
  typedef typename ConfiguratorType::RealType RealType;

  //! constructor from grid
  DiscreteVectorFunctionDefault ( const typename ConfiguratorType::InitType &Initializer,
                                  const aol::MultiVector<RealType> &Dofs,
                                  bool allowMoreDofsThanNeeded = false )
      : DiscreteVectorFunction<ConfiguratorType, DiscreteFunctionDefault<ConfiguratorType>, DimRange > ( Initializer, Dofs, allowMoreDofsThanNeeded ) { }

  //! constructor from configurator
  DiscreteVectorFunctionDefault ( const ConfiguratorType &Configurator,
                                  const aol::MultiVector<RealType> &Dofs,
                                  bool allowMoreDofsThanNeeded = false )
      : DiscreteVectorFunction<ConfiguratorType, DiscreteFunctionDefault<ConfiguratorType>, DimRange > ( Configurator, Dofs, allowMoreDofsThanNeeded ) { }
};

} // end namespace aol

#endif
