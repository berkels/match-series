#ifndef __HASHSET_H
#define __HASHSET_H

#include <aol.h>
#include <quoc.h>
#include <simplexGrid.h>

#ifdef  _MSC_VER
#include <hash_set>
#elif defined(_LIBCPP_VERSION)
#include <unordered_set>
#else
#if (__GNUC__ >= 4)
#include <tr1/unordered_set>
#else
#include <ext/hash_set>
#endif
#endif

/**
 * The hash_set implementation of MSC is very different from
 * the one of GCC. The following things are platform dependent
 * preparations that are necessary before we can derive a class
 * from hash_set.
 */

#ifndef  _MSC_VER

#if defined(_LIBCPP_VERSION)
#define STD_SUB_NAMESPACE __1
#elif (__GNUC__ >= 4)
#define STD_SUB_NAMESPACE tr1
#endif

//! GNU C++ Extensions (contains template scpecialization for hash)
#if (__GNUC__ >= 4)
namespace std {
// When using C++11 and libc++, the STD_SUB_NAMESPACE is inline.
#if ( defined (_LIBCPP_VERSION) ) && ( defined ( USE_CPP11 ) )
inline
#endif
namespace STD_SUB_NAMESPACE {
#else
namespace __gnu_cxx {
#endif

template <> struct hash<qc::Element> {
  size_t operator() ( const qc::CoordType& t ) const {
    return t.x() + ( t.z() * 769 + t.y() ) * 1543;
  }
};

template <> struct hash<qc::simplex::Element<qc::Element> > {
  size_t operator() ( const qc::simplex::Element<qc::Element> & El ) const {
    const qc::Element & cubicEl = El.getCubicElement();
    hash<qc::Element> cubicHash;
    size_t cubicHashValue = cubicHash ( cubicEl );
    return cubicHashValue * 6 + El.getSimplexNumber();
  }
};

#if (__GNUC__ >= 4)
} // end of namespace tr1.
} // end of namespace std.
#else
} // end of __gnu_cxx
#endif

#else

namespace {

struct ElementHashFn
{
  size_t operator() ( const qc::Element & El ) const {
    return El.x() + ( El.z() * 769 + El.y() ) * 1543;
  }

  size_t operator() ( const qc::simplex::Element<qc::Element> & El ) const {
    size_t cubicHashValue = (*this) ( El.getCubicElement() );
    return cubicHashValue * 6 + El.getSimplexNumber();
  }
};

struct ElementHashTraits
{
  static const size_t bucket_size = 4;
  static const size_t min_buckets = 8;

  template <typename T>
  size_t operator()(const T & key) const {
    return hash_func(key);
  }

  template <typename T>
  bool operator()(const T & lhs, const T & rhs) const {
    return hash_func(lhs) < hash_func (rhs);
  }

private:
  ElementHashFn hash_func;
};
} // end nameless namespace

#endif

namespace aol {

template <typename KeyType>
class HashSet;

/**
 * This class is a hash_set implementation for qc::Element that encapsulates
 * all platform dependent code. If a hash_set for other types than qc::Elemtent
 * is needed, a new template specialization needs to be added to accomplish this.
 *
 * \author Berkels
 */
template <>
class HashSet<qc::Element> : public
#ifdef  _MSC_VER
  stdext::hash_set<qc::Element, ElementHashTraits>
#else
#if (__GNUC__ >= 4)
  std::STD_SUB_NAMESPACE::unordered_set<qc::Element, std::STD_SUB_NAMESPACE::hash<qc::Element> >
#else
  hash_set<qc::Element, hash<qc::Element> >
#endif
#endif
{
#ifndef  _MSC_VER
  enum { numBuckets = 10000 };
#endif
public:
  HashSet () :
#ifdef  _MSC_VER
  stdext::hash_set<qc::Element, ElementHashTraits>()
#else
#if (__GNUC__ >= 4)
  std::STD_SUB_NAMESPACE::unordered_set<qc::Element, std::STD_SUB_NAMESPACE::hash<qc::Element> >( numBuckets ) // create a hash with a default of numBuckets buckets
#else
  hash_set<qc::Element, hash<qc::Element> >( numBuckets ) // create a hash with a default of numBuckets buckets
#endif
#endif
  {}
};

template <>
class HashSet<qc::simplex::Element<qc::Element> > : public
#ifdef  _MSC_VER
  stdext::hash_set<qc::simplex::Element<qc::Element>, ElementHashTraits>
#else
#if (__GNUC__ >= 4)
  std::STD_SUB_NAMESPACE::unordered_set<qc::simplex::Element<qc::Element>, std::STD_SUB_NAMESPACE::hash<qc::simplex::Element<qc::Element> > >
#else
  hash_set<qc::simplex::Element<qc::Element>, hash<qc::simplex::Element<qc::Element> > >
#endif
#endif
{
#ifndef  _MSC_VER
  enum { numBuckets = 10000 };
#endif
public:
  HashSet () :
#ifdef  _MSC_VER
  stdext::hash_set<qc::simplex::Element<qc::Element>, ElementHashTraits>()
#else
#if (__GNUC__ >= 4)
  std::STD_SUB_NAMESPACE::unordered_set<qc::simplex::Element<qc::Element>, std::STD_SUB_NAMESPACE::hash<qc::simplex::Element<qc::Element> > > ( numBuckets ) // create a hash with a default of numBuckets buckets
#else
  hash_set<qc::simplex::Element<qc::Element>, hash<qc::simplex::Element<qc::Element> > >( numBuckets ) // create a hash with a default of numBuckets buckets
#endif
#endif
  {}
};

} // end of namespace aol.

#endif // __HASHSET_H
