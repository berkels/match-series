#ifndef __LOOKUPMAP_H
#define __LOOKUPMAP_H

#include <aol.h>

namespace aol {

/** Extension of STL set. Provides method for checking whether an element is contained.
 *  \author Schwen (MEVIS)
 */
template < typename StoredType >
class LookupSet : public std::set < StoredType > {
public:
  // default standard constructor, copy constructor, assignment operator and destructor are correct

  // provide access to iterator typedefs
  using typename std::set < StoredType >::iterator;
  using typename std::set < StoredType >::const_iterator;

  //! return whether Key is already contained
  bool contains ( const StoredType &Key ) const {
    typename std::set < StoredType >::const_iterator kit = this->find ( Key );
    return ( kit != this->end() );
  }

  //! remove KeyOld if it exists and insert KeyNew (actually add it if it is not contained yet)
  std::pair < typename std::set<StoredType>::iterator, bool > replace ( const StoredType &KeyOld, const StoredType &KeyNew ) {
    typename std::set < StoredType >::const_iterator kit = this->find ( KeyOld );
    if ( kit != this->end() ) {
      this->erase ( kit );
    }
    return ( this->insert ( KeyNew ) );
  }

  //! dump contents
  ostream& dump ( ostream &Out = cout ) const {
    for ( typename std::set < StoredType >::const_iterator it = this->begin(); it != this->end(); ++it ) {
      Out << (*it) << endl;
    }
    return ( Out );
  }
};


template < typename StoredType >
ostream &operator<< ( ostream &os, const LookupSet<StoredType> &LSet ) {
  return ( LSet.dump ( os ) );
}


/** Extension of STL map. LookupMap allows setting entries only once, get checks whether entries exist and is not limited by non-const access method.
 *  \author Schwen (MEVIS)
 */
template < typename FromType, typename ToType >
class LookupMap : public std::map < FromType, ToType > {
public:
  // default standard constructor, copy constructor, assignment operator and destructor are correct

  // provide access to iterator typedefs
  using typename std::map < FromType, ToType >::iterator;
  using typename std::map < FromType, ToType >::const_iterator;

  //! set a new entry, throw exception if the entry exists already
  void set ( const FromType &NewKey, const ToType &NewVal ) {
    typename std::map < FromType, ToType >::const_iterator kit = this->find ( NewKey );
    if ( kit != this->end() ) {
      cerr << NewKey << " exists." << endl;
      throw aol::Exception ( "LookupMap::set: entry exists", __FILE__, __LINE__ );
    }
    std::map < FromType, ToType >::operator[] ( NewKey ) = NewVal;
  }

  //! reset entry if it exists, otherwise throw exception
  void reset ( const FromType &Key, const ToType &NewVal ) {
    typename std::map < FromType, ToType >::const_iterator kit = this->find ( Key );
    if ( kit == this->end() ) {
      cerr << Key << " does not exists." << endl;
      throw aol::Exception ( "LookupMap::reset: entry does not exist", __FILE__, __LINE__ );
    }
    std::map < FromType, ToType >::operator[] ( Key ) = NewVal;
  }

  //! get reference to entry or throw exception if entry is not contained
  const ToType & getRef ( const FromType &Key ) const {
    typename std::map < FromType, ToType >::const_iterator kit = this->find ( Key );
    if ( kit != this->end() ) {
      return ( kit->second );
    } else {
      cerr << Key << " not found." << endl;
      throw aol::Exception ( "LookupMap::getRef: entry not found", __FILE__, __LINE__ );
    }
  }

  ToType & getRef ( const FromType &Key ) {
    typename std::map < FromType, ToType >::iterator kit = this->find ( Key );
    if ( kit != this->end() ) {
      return ( kit->second );
    } else {
      cerr << Key << " not found." << endl;
      throw aol::Exception ( "LookupMap::getRef: entry not found", __FILE__, __LINE__ );
    }
  }

  //! get reference to entry or return a given default entry
  const ToType & getRefOrReturnDefault ( const FromType &Key, const ToType &DefaultValue ) const {
    typename std::map < FromType, ToType >::const_iterator kit = this->find ( Key );
    if ( kit != this->end() ) {
      return ( kit->second );
    } else {
      return ( DefaultValue );
    }
  }

  //! get reference to contained entry or newly created entry as in map::operator[]
  ToType & getRefOrCreate ( const FromType &Key ) {
    return ( std::map< FromType, ToType >::operator[] ( Key ) );
  }

  //! check whether entry is contained
  bool contains ( const FromType &Key ) const {
    typename std::map < FromType, ToType >::const_iterator kit = this->find ( Key );
    return ( kit != this->end() );
  }

  //! fill other LookupMap with the inverse lookup, throw an exception if map is not injective
  void invertTo ( aol::LookupMap < ToType, FromType > & Other ) const {
    for ( typename std::map < FromType, ToType >::const_iterator it = this->begin(); it != this->end(); ++it ) {
      Other.set ( it->second, it->first );
    }
  }

  //! fill this LookupMap with the inverse lookup of Other, throw an exception if map is not injective
  void invertFrom ( const aol::LookupMap < ToType, FromType > & Other ) {
    for ( typename std::map < ToType, FromType >::const_iterator it = Other.begin(); it != Other.end(); ++it ) {
      this->set ( it->second, it->first );
    }
  }

  //! returns whether there is at least one value to which more than one key is mapped
  bool isInjective ( ) const {
    aol::LookupMap< ToType, FromType > wouldBeInverse;
    try {
      this->invertTo ( wouldBeInverse );
    } catch ( aol::Exception &Ex ) {
      Ex.consume();
    }
    return ( this->size() == wouldBeInverse.size() );
  }

  //! dump contents
  ostream& dump ( ostream &Out = cout ) const {
    for ( typename std::map < FromType, ToType >::const_iterator it = this->begin(); it != this->end(); ++it ) {
      Out << it->first << " " << it->second << endl;
    }
    return ( Out );
  }

  //! erase all but the Num first entries (i.e. all but the Num smallest entries), do nothing if the map is smaller
  void eraseButFirstNumEntries ( const unsigned int Num ) {
    unsigned int currentSize = this->size();
    if ( Num < currentSize ) {
      typename std::map < FromType, ToType >::iterator startErase = this->begin();
      for ( unsigned int n = 0; n < Num; ++n ) {
        ++startErase;
      }
      this->erase ( startErase, this->end() );
    }
  }

  //! copy the first (up to) Num values to an std::vector
  void copyFirstValuesTo ( const unsigned int Num, std::vector< ToType > &Dest ) const {
    Dest.clear();
    Dest.reserve ( Num );
    unsigned short cntr = 0;
    for ( typename std::map < FromType, ToType >::const_iterator it = this->begin(); it != this->end() && cntr < Num; ++it, ++cntr ) {
      Dest.push_back ( it->second );
    }
  }

  static bool isSelfTestOK ( ); // only implemented for <int, int>, but does not depend on types anyway.

};


template < typename FromType, typename ToType >
/*static*/ bool aol::LookupMap<FromType,ToType>::isSelfTestOK ( ) {
  bool success = true;
  {
    aol::LookupMap< int, int > lookupMap;
    for ( int i = 0; i < 20; ++i ) {
      lookupMap [ (i * i) % 23 ] = i;
    }
    lookupMap.set ( 5, 42 ); // 5 does not exist yet
    lookupMap.reset ( 13, -17 ); // 13 already exists
    lookupMap.eraseButFirstNumEntries ( 7 );
    success &= ( lookupMap.size() == 7 );
  }
  {
    aol::LookupMap<short, unsigned short> injMap, nonInjMap, injMapInvInv;
    aol::LookupMap<unsigned short, short> injMapInv;
    injMap.set ( 1, 1 );
    injMap.set ( 2, 2 );
    nonInjMap.set ( 1, 1 );
    nonInjMap.set ( 2, 1 );
    injMap.invertTo ( injMapInv );
    injMapInvInv.invertFrom ( injMapInv );
    success &= ( ( injMap.isInjective() == true ) && ( nonInjMap.isInjective() == false ) && ( injMap == injMapInvInv ) );
  }
  return ( success );
}


template < typename FromType, typename ToType >
ostream &operator<< ( ostream &os, const LookupMap < FromType, ToType > &LMap ) {
  return ( LMap.dump ( os ) );
}

}

#endif
