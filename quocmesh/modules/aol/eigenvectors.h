#ifndef __EIGENVECTORS_H
#define __EIGENVECTORS_H

#include <op.h>
#include <matrix.h>
#include <matrixInverse.h>
#include <progressBar.h>
#include <vectorExtensions.h>
#include <preconditioner.h>
#include <solver.h>
#ifdef USE_EXTERNAL_EIGEN
#include <eigenIncludes.h>
#endif

namespace aol {

//**************************************************************************
/**
 *  \brief base class for ordering operators
 *
 *  By calling apply(), the "dest" argument is reordered such that an
 *  analogous reordering of "arg" would give an ascendingly ordered vector.
 *
 *  Thus, when you pass the identity as permutation matrix P or an int-vector
 *  i with entries i[j] = j, P * arg or reorder ( i, arg ) gives an
 *  ascendingly ordered vector.
 *
 *  apply() can be called with a permutation matrix or an int-vector.
 *  Subclasses will only implement the int-vector version.
 *
 *  \author von Deylen
 */
template <typename RealType>
class OrderOp : public Op<Vector<RealType>, PermutationMatrix<RealType> > {
public:
  //! subclasses do not need to use the PermutationMatrix interface,
  //! but can directly act on a int-vector.
  virtual void apply ( const Vector<RealType> & arg, Vector<int> & dest ) const = 0;

  virtual void apply ( const Vector<RealType> & arg, PermutationMatrix<RealType> & dest ) const {
    Vector<int> dest_vec ( dest.getPermutationVectorRef() );
    apply ( arg, dest_vec );
    dest.setPermutationVector ( dest_vec );
  }

  virtual void applyAdd ( const Vector<RealType> &, PermutationMatrix<RealType> & ) const {
    throw aol::Exception ( "OrderOp::applyAdd(...) is useless.", __FILE__, __LINE__ );
  }

  void reorder ( const PermutationMatrix<RealType> & orderMatrix, Vector<RealType> & toReorder ) {
    reorder ( orderMatrix.getPermutationVectorRef(), toReorder );
  }

  void reorder ( const Vector<int> & order, Vector<RealType> & toReorder ) {
    Vector<RealType> temp ( toReorder );
    PermutationMatrix<RealType> ( order ).apply ( temp, toReorder );
  }
};

//**************************************************************************
/**
 *  \brief Ordering operator using heapsort
 *
 *  for explaination of this algorithm implementation see Numerical Recipes
 *  chap. 8.3 "Heapsort".
 *
 *  \author von Deylen
 */
template <typename RealType>
class HeapsortOp : public OrderOp<RealType> {
public:
  typedef Vector<RealType> VectorType;
  typedef Vector<int> OrderType;

  using OrderOp<RealType>::apply;

  void apply ( const VectorType & arg, OrderType & dest ) const {
    QUOC_ASSERT ( arg.size() == dest.size() );
    int n = arg.size();
    for ( int i = n / 2 - 1; i >= 0; --i )
      siftDown ( arg, dest, i, n - 1 );
    for ( int i = n - 1; i > 0; --i ) {
      swap ( dest[0], dest[i] );
      siftDown ( arg, dest, 0, i - 1 );
    }
  }

protected:
  //! bewegt das Element dest[anfang] solange im Heap nach unten,
  //! bis die Heap-Bedingung wieder erfuellt ist.
  //! \pre Heap-Bedingung ist fuer alle j > anfang erfuellt.
  //! \post Heap-Bedingung ist fuer fuer alle j >= anfang erfuellt.
  void siftDown ( const VectorType & values, OrderType & keys, int anfang, int ende ) const {

    int a_index = keys[anfang];
    int j_old = anfang;
    int j     = 2 * anfang + 1;
    // solange wir nicht am unteren Ende des Baums sind:
    while ( j <= ende ) {
      // wenn rechtes Element groesser: Gehe zum rechten Teilbaum ueber
      if ( j < ende && values[keys[j]] < values[keys[j+1]] )
        j++;
      // jetzt ist values[keys[j]] der groessere der beiden Werte auf
      // dieser Baumebene.
      // Wenn bisher nicht eingepasstes Element groesser ist als
      // alle Elemente auf dieser Baumebene: fertig.
      if ( values[a_index] >= values[keys[j]] )
        break;
      // schiebe Baumelemente nach oben
      keys[j_old] = keys[j];
      j_old = j;
      // setze Index j auf linkes Kind von bisherigem j
      j = 2 * j + 1;
    }
    keys[j_old] = a_index;
  }
};

//**************************************************************************
/**
 *  \brief Parent class of all eigenvector and -value computing operators.
 *
 *  The eigenvector operators compute their result via
 *  apply(const MatrixType & matrix, MultiVector & dest).
 *  It always stores its result as follows: dest[0] contains
 *  the eigenvalues. dest may, but does not have to, have
 *  additional entrys (at most dest[0].size() ones), which are
 *  the eigenvectors corresponding to the eigenvalues in dest[0].
 *  So, if dest.size() is at least k+2, in dest[k+1] is stored
 *  the eigenvector corresponding to dest[0][k]:
 *
 *     matrix * dest[k+1] = dest[0][k] * dest[k+1].
 *
 *  It depends on the implementation of a specific EigenvectorOp
 *  child class how many eigenvalues are computed and if also
 *  eigenvectors are computed. When calling apply(...), the destination
 *  MultiVector is appropriately resized.
 *
 *  \author von Deylen
 */
template <typename MatrixType>
class EigenvectorOp : public Op<MatrixType, MultiVector<typename MatrixType::DataType> > {
public:
  typedef typename MatrixType::DataType DataType;
  void applyAdd ( const MatrixType &, MultiVector<typename MatrixType::DataType> & ) const {
    throw Exception ( "EigenvectorOp::applyAdd(...) called.", __FILE__, __LINE__ );
  }

protected:
  template <typename FullMatrixType>
  void orderEigenvaluesAndCopyEigenvectors ( const FullMatrixType & vectors,
                                             MultiVector<typename MatrixType::DataType> & dest ) const {

    int n = dest.numComponents() - 1;
    PermutationMatrix<DataType> indexReorderMatrix ( n );

    HeapsortOp<DataType> heapsort;
    heapsort.apply ( dest[0], indexReorderMatrix );
    // bring eigenvalues in dest[0] into ascending order
    // use Heapsort::reorder because Matrix::apply
    // does not support Arg == Dest.
    heapsort.reorder ( indexReorderMatrix, dest[0] );

    // write eigenvectors from vectors into dest[i+1]
    // in eigenvalues' ascending order
    Vector<int> indexOrder ( indexReorderMatrix.getPermutationVectorRef() );
    for ( int i = 0; i < n; ++i )
      vectors.getColumn ( indexOrder[i], dest[i + 1] );
  }
};


//**************************************************************************
/**
 *  \brief Computes the leading eigenvalues and eigenvectors of a matrix.
 *
 *  \author Paetz
 */
template <typename MatrixType>
class DeflationEigenvectorOp : public EigenvectorOp<MatrixType> {
public:
  typedef typename MatrixType::DataType RealType;

  DeflationEigenvectorOp ( int numberEigenvalues = 1, bool Quiet = false )
    : _quietMode ( Quiet ),
      _numberEigenvalues ( numberEigenvalues ),
      _relativeAccuracy ( 0.0 ),
      _threshold ( 0.0 ),
      _sumEvals ( 0.0 ) { }

  void apply ( const MatrixType & arg, MultiVector<RealType> & dest ) const {
    int n = arg.getNumCols();
    int nEvals = aol::Max<int> ( _numberEigenvalues, 1 );
    
    if ( nEvals > n ) throw aol::Exception ( aol::strprintf ( "User requested %d eigen values, but the specified matrix is only of size %dx%d", nEvals, n, n ).c_str ( ) , __FILE__, __LINE__ );
    
    if ( _numberEigenvalues == 0 && ( _sumEvals == 0 || _threshold == 0 ) ) throw aol::Exception ( "Neither number of eigenvalues nor threshold and sum of eigenvalues was specified!", __FILE__, __LINE__ );
    
    dest.reallocate ( nEvals + 1, n );
    Vector<RealType> newVec ( n ), old ( n ), deflationTmp ( n ), tmp ( n );
    newVec.setAll ( 1.0 );
    old.setAll ( 1.0 );
    ProgressBar<> pb ( "Computing eigenvalues " );
    if ( _numberEigenvalues > 0 && !_quietMode ) pb.start ( _numberEigenvalues );
    for ( int i = 0; i < nEvals; i++ ) {
      RealType res = 1.0, resOld;
      int numNonDecreasingIterations = 0;
      while ( res > 1e-12 && numNonDecreasingIterations < 10 ) {
        resOld = res;
        arg.apply ( old, newVec );
        for ( int d = 0; d < i; d++ ) {
          RealType t = old.dotProduct ( dest[d+1] ) * dest[0][d];
          deflationTmp = dest[d+1];
          deflationTmp *= t;
          newVec -= deflationTmp;
        }
        newVec /= newVec.norm();
        tmp = old;
        tmp -= newVec;
        res = tmp.norm();
        old = newVec;
        if ( res >= resOld ) ++numNonDecreasingIterations;
      }
      arg.apply ( newVec, old );
      for ( int d = 0; d < i; d++ ) {
        RealType t = newVec.dotProduct ( dest[d+1] ) * dest[0][d];
        deflationTmp = dest[d+1];
        deflationTmp *= t;
        old -= deflationTmp;
      }
      dest[0][i] = newVec.dotProduct ( old );
      dest[i+1] = newVec;
      
      if ( _numberEigenvalues == 0 ) {
        RealType error = _sumEvals;
        for ( int j=0; j<=i ; ++j ) error -= dest[0][j];
        
        if ( !_quietMode ) std::cerr << error << " > " << _threshold << std::endl;
        
        if ( error < _threshold ) break;
        else {
          ++nEvals;
          dest.resize ( nEvals + 1, n );
        }
      } else if ( _relativeAccuracy > 0 && aol::Abs<RealType> ( dest[0][i] / dest[0][0] ) < _relativeAccuracy ) {
        dest.resize ( i+1, n );
        break;
      }
      
      if ( _numberEigenvalues > 0 && !_quietMode ) pb++;
    }
    if ( _numberEigenvalues > 0 && !_quietMode ) pb.finish();
  }

  void setNumberEigenvalues ( int numberEigenvalues ) {
    _numberEigenvalues = numberEigenvalues;
  }
  
  void setRelativeAccuracy ( RealType relativeAccuracy ) {
    _relativeAccuracy = relativeAccuracy;
  }
  
  void setThreshold ( const RealType Threshold, const RealType SumEvals ) {
    _threshold = Threshold;
    _sumEvals = SumEvals;
  }
  
  void setQuietMode ( bool Quiet = true ) {
    _quietMode = Quiet;
  }

protected:
  bool _quietMode;
  int _numberEigenvalues;
  RealType _relativeAccuracy;
  RealType _threshold, _sumEvals;
};


//**************************************************************************
/**
 *  \brief Computes all eigenvalues, and all eigenvectors,
 *         of a symmetrix matrix. The template matrix
 *         class has to support set to every possible entry,
 *         thus full matrix classes are prefereable.
 *
 *  Typical matrices require 6 to 10 sweeps to achieve convergence,
 *  or 3n^2 to 5n^2 Jacobi rotations. Each rotation requires of order
 *  4n operations, so total afford is of order 12n^3 to 20n^3.
 *  Calculation of the eigenvectors as well as eigenvalues is only
 *  50% overhead.
 *
 *  See Numerical Recipes, 11.1 "Jacobi Transformation" for detailled
 *  description.
 *
 *  \author von Deylen
 */
template <typename MatrixType>
class JacobiEigenvectorOp : public EigenvectorOp<MatrixType> {
public:
  typedef typename MatrixType::DataType RealType;

  JacobiEigenvectorOp() : _maxNumberOfSweeps ( 100 ) {}

  void apply ( const MatrixType & arg, MultiVector<RealType> & dest ) const {
    int n = arg.getNumCols();

    // initialize
    dest.reallocate ( n + 1, n );
    // allocate additional space to have
    // eigenvectors in matrix order.
    // After finishing, they are copied
    // into dest.
    MatrixType eigenvectors ( arg );
    eigenvectors.setIdentity();

    MatrixType a ( arg );
    Vector<RealType> & eigenvalues = dest[0];
    for ( int i = 0; i < n; i++ )
      eigenvalues[i] = a.getDiag ( i );
    Vector<RealType> b ( eigenvalues );
    Vector<RealType> z ( n );

    int iter;
    for ( iter = 0; iter < _maxNumberOfSweeps; ++iter ) { // begin rotation sequence

      // check for convergence:
      RealType sm = 0.;
      for ( int i = 0; i < n - 1; ++i )
        for ( int j = i + 1; j < n; ++j )
          sm += fabs ( a.get ( i, j ) );

      if ( sm == 0. )
        break;

      RealType tresh;
      if ( iter < 3 ) tresh = sm / ( 5. * n * n ); // first 3 sweeps
      else            tresh = 0.;

      for ( int i = 0; i < n - 1; ++i )
        for ( int j = i + 1; j < n; ++j ) {
          RealType g = 100. * fabs ( a.get ( i, j ) );

          // after 4 sweeps
          if ( iter > 3 && ( fabs ( eigenvalues[i] ) + g ) == fabs ( eigenvalues[i] )
               && ( fabs ( eigenvalues[j] ) + g ) == fabs ( eigenvalues[j] ) )
            a.set ( i, j, 0. );
          else if ( fabs ( a.get ( i, j ) ) > tresh ) {
            RealType h = eigenvalues[j] - eigenvalues[i];
            RealType t;
            if ( ( fabs ( h ) + g ) == fabs ( h ) )
              t = ( a.get ( i, j ) ) / h;
            else {
              RealType theta = h / ( 2. * a.get ( i, j ) );
              t = 1. / ( fabs ( theta ) + sqrt ( 1. + theta * theta ) );
              if ( theta < 0 ) t = -t;
            }
            RealType c = 1. / sqrt ( 1 + t * t );
            RealType s = t * c;
            RealType tau = s / ( 1. + c );
            h = t * a.get ( i, j );
            z[i] -= h;
            z[j] += h;
            eigenvalues[i] -= h;
            eigenvalues[j] += h;
            a.set ( i, j, 0. );

            // i already shifted left by 1 unit
            for ( int k = 0; k < i; ++k )
              rotate ( a, k, i, k, j, g, h, s, tau );

            // i and j already shifted left by 1 unit
            for ( int k = i + 1; k < j; ++k )
              rotate ( a, i, k, k, j, g, h, s, tau );

            // j already shifted left by 1 unit
            for ( int k = j + 1; k < n; ++k )
              rotate ( a, i, k, j, k, g, h, s, tau );

            for ( int k = 0; k < n; ++k )
              rotate ( eigenvectors, k, i, k, j, g, h, s, tau );
          }
        }
      b += z;
      eigenvalues = b;
      z.setAll ( 0. );
    }

    if ( iter >= _maxNumberOfSweeps )
      throw Exception ( "JacobiEVOp::appy(...) has exceeded max number of sweeps.", __FILE__, __LINE__ );

    this->orderEigenvaluesAndCopyEigenvectors ( eigenvectors, dest );
  }

  void setMaxNumberOfSweeps ( int maxSweeps ) {
    _maxNumberOfSweeps = maxSweeps;
  }

protected:
  inline void rotate ( MatrixType & M,
                       int &i, int &j, int &k, int &l,
                       RealType &g, RealType &h,
                       RealType &s, RealType &tau ) const {
    g = M.get ( i, j );
    h = M.get ( k, l );
    M.set ( i, j, g - s * ( h + g * tau ) );
    M.set ( k, l, h + s * ( g - h * tau ) );
  }

  int _maxNumberOfSweeps;
};

//**************************************************************************
/**
 *  \brief QL algorith with implicit shift. Determines eigenvalues
 *         (and eigenvectors, if desired) of a real symmetric
 *         tridiagonal matrix.
 *
 *  Other classes (as QREigenvectorOp) use a combination of
 *  QRGivensTridiag and this class to compute eigenvalues and -vectors
 *  of non-tridiagonal matrices.
 *
 *  Upon calling apply(...) with any matrix, it is converted into a TriBandMatrix
 *  to check whether it really was tridiagonal. Afterwards, symmetry
 *  is checked. Both can be turned of by the second constructor
 *  argument.
 *
 *  Workload ca. 3n^3, could be 30n^2 for only eigenvalue computation
 *  (code pieces not to perform are marked, but interface lacks).
 *
 *  For further explainations on the algorithm, confer Numerical
 *  Recipes sec. 11.3 "QL Algorithm with Implicit Shifts".
 *
 *  \author von Deylen
 */
template <typename BandMatrixType, typename OrthogonalMatrixType>
class QLTridiagEigenvectorOp : public EigenvectorOp<BandMatrixType> {
public:
  typedef typename BandMatrixType::DataType DataType;

  QLTridiagEigenvectorOp ( const OrthogonalMatrixType & O, bool checkTridiagSymm = true ) :
      _maxIter ( 300 ),
      _O ( O ),
      _checkTridiagSymm ( checkTridiagSymm ) {}

  void setMaxIter ( int maxIter ) {
    _maxIter = maxIter;
  }

  int getMaxIter() const {
    return _maxIter;
  }

  void apply ( const BandMatrixType & arg, MultiVector<DataType> & dest ) const {
    QUOC_ASSERT ( arg.getNumCols() == arg.getNumRows() );
    int n = arg.getNumRows();
    if ( _checkTridiagSymm )
      isTridiagAndSymm ( arg );

    Vector<DataType> diag ( n ),
    subdiag ( n );

    for ( int i = 0; i < n - 1; ++i ) {
      diag[i] = arg.getDiag ( i );
      subdiag[i] = arg.get ( i + 1, i );
    }
    diag[n - 1] = arg.getDiag ( n - 1 );

    OrthogonalMatrixType eigenvectors ( _O );
    performQLAlgorithm ( diag, subdiag, eigenvectors );
    dest.reallocate ( n + 1, n );
    dest[0] = diag;
    this->orderEigenvaluesAndCopyEigenvectors ( eigenvectors, dest );
  }

  bool isTridiagAndSymm ( const BandMatrixType & m ) const {
    int n = m.getNumRows();

    // m is tridiagonal if operator+= does not try to
    // to write to non-existing entries in triBandMatrix.
    // Thus, tridiagonal check is o.k. if no exception
    // is thrown.
    TriBandMatrix<DataType> triBandMatrix ( n );
    try {
      triBandMatrix += m;
    } catch ( ... ) {
      throw Exception ( "QLTridiagEigenvectorOp: passed matrix is not tridiagonal." );
    }

    // symmetry check
    DataType diff;
    for ( int i = 1; i < n; ++i )
      if ( fabs ( diff = triBandMatrix.get ( i, i - 1 ) - triBandMatrix.get ( i - 1, i ) ) > 1E-10 ) {
        string msg = strprintf ( "QLTridiagEigenvectorOp: passed matrix is not symmetric: diff(%i,%i) = %f.", i, i - 1, diff );
        throw Exception ( msg, __FILE__, __LINE__ );
      }
    return true;
  }

  void performQLAlgorithm ( Vector<DataType> & d, Vector<DataType> & e, OrthogonalMatrixType & z ) const {
    int n = d.size();
    // on entry, last entry subdiag is not neccessarily
    // set to anything.
    e[n - 1] = 0.;
    // ProgressBar<false, true, int, 50> progressBar ( "eigenvalue computation " ); // does not seem to work
    ProgressBar<> progressBar ( "eigenvalue computation " );
    progressBar.start ( n );
    for ( int l = 0; l < n; ++l ) {
      progressBar++;
      int iter = 0;
      int m;
      do {
        // Look for a single small subdiagonal element
        // to split the matrix.
        for ( m = l; m < n - 1; ++m ) {
          DataType dd = fabs ( d[m] ) + fabs ( d[m+1] );
          if ( fabs ( e[m] ) + dd == dd )
            break;
        }
        if ( m != l ) {
          if ( ++iter == _maxIter )
            throw Exception ( "Maximum number of iterations exceeded in "
                              "QLTridiagEigenvectorOp::performQLAlgorithm()." );
          // Form shift
          DataType g = ( d[l + 1] - d[l] ) / ( 2. * e[l] );
          DataType r = sqrt ( g * g + 1. );
          // d_m - k_s
          g = d[m] - d[l] + e[l] / ( g + sign_prod ( r, g ) );
          DataType s = 1.,
            c = 1.,
            p = 0.;
          int i;
          for ( i = m - 1; i >= l; --i ) {
            // a plane rotation as in the original QL,
            // followed by Givens rotations to restore
            // tridiagonal form.
            DataType f = s * e[i],
              b = c * e[i];
            e[i + 1] = r = sqrt ( f * f + g * g );
            // recover from underflow:
            if ( r == 0. ) {
              d[i + 1] -= p;
              e[m] = 0.;
              break;
            }
            s = f / r;
            c = g / r;
            g = d[i + 1] - p;
            r = ( d[i] - g ) * s + 2. * c * b;
            d[i + 1] = g + ( p = s * r );
            g = c * r - b;
            // next loop can be omitted if eigenvectors
            // not wanted:
            for ( int k = 0; k < n; ++k ) {
              f = z.get ( k, i + 1 );
              z.set ( k, i + 1, s * z.get ( k, i ) + c * f );
              z.set ( k, i, c * z.get ( k, i ) - s * f );
            }
          }
          if ( r == 0 && i >= l )
            continue;
          d[l] -= p;
          e[l] = g;
          e[m] = 0.;
        }
      } while ( m != l );
    }
    progressBar.finish();
  }

  template <typename T>
  T sign_prod ( const T & a, const T & b ) const {
    return ( b > 0 ? ( a >= 0 ? a : -a ) : ( a >= 0 ? -a : a ) );
  }

  //! maximum number of iterations performed for finding
  //! each eigenvalue (iterations needed for first few
  //! eigenvalues approx. 4-5, later ca. 1.3-1.6)
  int _maxIter;

  //! orthogonal transformation from initial matrix
  //! to tridiagonal form
  OrthogonalMatrixType _O;

  //! save whether argument of apply(...) has to be checked
  //! for tridiagonal form and symmetry
  bool _checkTridiagSymm;
};

//**************************************************************************
/**
 *  \brief Composition of QRGivensTridiag and QLTridiagEigenvectorOp.
 *         Computes eigenvalues by first producing tridiagonal
 *         form and then finding eigenvalues of this tridiagonal
 *         matrix.
 *
 *  For reference on the algorithms, see documentations of the
 *  used classes QRGivensTridiag and QLTridiagEigenvectorOp.
 *
 *  \author von Deylen
 */
template <class MatrixType>
class QREigenvectorOp : public EigenvectorOp<MatrixType> {
public:
  typedef typename MatrixType::DataType DataType;

  void apply ( const MatrixType & arg, MultiVector<DataType> & dest ) const {
    TriBandMatrix<DataType> tridiag ( arg.getNumRows() );
    MatrixType Q ( arg );

    // perform tridiagonalization in constructor
    QRGivensTridiag<MatrixType> tridiagonlizer ( arg );
    tridiagonlizer.getQ ( Q );
    tridiagonlizer.getTridiagBand ( tridiag );

    // perform QL iteration
    QLTridiagEigenvectorOp<TriBandMatrix<DataType>, MatrixType > tridiagEVOp ( Q );
    tridiagEVOp.apply ( tridiag, dest );
  }
};

/**
 * \brief For a number of samples (e.g. vectors or matrices), this operator computes the arithmetic mean
 * and (if desired) the input samples with the mean subtracted.
 */
template <typename SampleType>
class MeanOp :
      public aol::Op<aol::RandomAccessContainer<SampleType>, SampleType> {

public:
  MeanOp() {}

  /**
   * Returns the arithmetic mean.
   */
  void applyAdd ( const aol::RandomAccessContainer<SampleType> &Arg, SampleType &Dest ) const {

    // number of samples/probes
    int n = Arg.size();

    // compute the arithmetic mean
    SampleType mean ( Arg[0], aol::STRUCT_COPY );
    for ( int i = 0; i < n; i++ )
      mean += Arg[i];
    mean /= n;

    Dest += mean;
  }

  /**
   * Returns the arithmetic mean and the input samples minus the mean.
   */
  void meanAndVariation ( const aol::RandomAccessContainer<SampleType> &Samples, SampleType &Mean, aol::RandomAccessContainer<SampleType> &Variation ) const {

    // number of samples/probes
    int n = Samples.size();

    // compute the arithmetic mean
    for ( int i = 0; i < n; i++ )
      Mean += Samples[i];
    Mean /= n;

    // subtract mean from samples
    Variation = Samples;
    for ( int i = 0; i < n; i++ )
      Variation[i] -= Mean;
  }
};

/**
 * \brief For a number of samples (e.g. vectors or matrices) with zero mean, this operator computes the covariance matrix.
 * If the mean is nonzero, then the method "applyForNonZeroMean" has to be used, which first subtracts the mean from all
 * input samples. The dot product of the space has to be ptovided.
 */
template <typename RealType, typename SampleType>
class CovarianceMatrixOp :
      public aol::Op<aol::RandomAccessContainer<SampleType>, aol::FullMatrix<RealType> > {

private:
  // the symmetric linear operator D so that the dot product <x,y> is defined as x^TDy
  const aol::Op<SampleType> &_dotProd;

public:
  /**
   * \brief The dotproduct to be used for the covariance matrix has to be specified,
   * i.e. the symmetric linear operator D so that the dot product <x,y> is defined as x^TDy.
   */
  CovarianceMatrixOp ( const aol::Op<SampleType> &DotProd ) :
      _dotProd ( DotProd ) {}

  /**
   * \brief "Arg" contains all given data points (e.g. of type "aol::Vector" or "aol::MultiVector"),
   * that must have mean zero (else use "applyForNonzeroMean").
   * "Dest" will contain the corresponding covariance matrix.
   */
  void apply ( const aol::RandomAccessContainer<SampleType> &Arg, aol::FullMatrix<RealType> &Dest ) const {

    // number of samples/probes
    int n = Arg.size();

    // compute the covariance matrix, defined as <x_i,x_j>
    Dest.reallocate ( n, n );
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < n; i++ ) {
      SampleType aux ( Arg[0], aol::STRUCT_COPY );
      _dotProd.apply ( Arg[i], aux );
      Dest.set ( i, i, aux * Arg[i] );
      for ( int j = i + 1; j < n; j++ ) {
        Dest.set ( i, j, aux * Arg[j] );
        Dest.set ( j, i, Dest.get ( i, j ) );
      }
    }
  }

  /**
   * \brief "Arg" contains all given data points (e.g. of type "aol::Vector" or "aol::MultiVector"),
   * "Dest" will contain the corresponding covariance matrix.
   */
  void applyForNonzeroMean ( const aol::RandomAccessContainer<SampleType> &Arg, aol::FullMatrix<RealType> &Dest ) const {

    // compute the arithmetic mean and the Variations
    SampleType mean ( Arg[0], aol::STRUCT_COPY );
    aol::RandomAccessContainer<SampleType> variations ( Arg );
    ( MeanOp<SampleType>() ).meanAndVariation ( Arg, mean, variations );

    // compute the covariance matrix, defined as <x_i,x_j>
    apply ( variations, Dest );
  }

  void applyAdd ( const aol::RandomAccessContainer<SampleType> &, aol::FullMatrix<RealType> & ) const {
    throw aol::Exception ( "CovarianceMatrixOp::applyAdd(...) called.", __FILE__, __LINE__ );
  }
};

/**
 * \brief Class for performing a principal component analysis on data of the type "aol::Vector" or "aol::MultiVector"
 * (or any type providing the same class methods and fitting into a "aol::RandomAccessContainer").
 * The dot product used for the covariance matrix can be freely chosen.
 *
 * \author Wirth
 */
template <typename RealType, typename SampleType>
class PCAOp :
      public aol::Op<aol::VectorContainer<SampleType> > {

private:
  // the symmetric linear operator D so that the dot product <x,y> is defined as x^TDy
  const aol::Op<SampleType> &_dotProd;
  // number of modes to be returned
  const int _numModes;

public:
  /**
   * \brief The number of modes of variation that shall be returned and the dotproduct to be used for the covariance matrix have to be specified.
   */
  PCAOp ( const aol::Op<SampleType> &DotProd, int NumModes ) :
      _dotProd ( DotProd ),
      _numModes ( NumModes ) {}

  /**
   * \brief "Samples" contains all given data points (e.g. "aol::Vector" or "aol::MultiVector"),
   * "Modes" will contain the first modes of variation and "Variances" the corresponding variances.
   */
  void modesAndVariances ( const aol::VectorContainer<SampleType> &Samples, aol::VectorContainer<SampleType> &Modes, aol::Vector<RealType> &Variances ) const {

    // number of samples/probes
    int n = Samples.size();

    // compute the mean and variations of the input data
    aol::VectorContainer<SampleType> variations;
    SampleType mean ( Samples[0], aol::STRUCT_COPY );
    ( MeanOp<SampleType>() ).meanAndVariation ( Samples, mean, variations );

    // compute the covariance matrix of the input data
    aol::FullMatrix<RealType> covMatrix;
    ( CovarianceMatrixOp<RealType, SampleType> ( _dotProd ) ).apply ( variations, covMatrix );

    // decompose covariance matrix (first row of decompResult contains eigenvalues in decreasing order, further rows contain corresponding eigenvectors)
    aol::MultiVector<RealType> decompResult;
    ( aol::QREigenvectorOp<aol::FullMatrix<RealType> >() ).apply ( covMatrix, decompResult );

    // compute modes of variation (note that eigenvalues from eigenvalue operator are in increasing absolute value order)
    int numModes = aol::Min ( n, _numModes );
    Modes.reallocate ( numModes );
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < numModes; i++ ) {
      Modes[i].reallocate ( variations[0] );
      for ( int j = 0; j < n; j++ )
        Modes[i].addMultiple ( variations[j], decompResult[n-i][j] );
      Modes[i] /= sqrt ( decompResult[0][n-i-1] );
    }

    // compute variances
    Variances.reallocate ( numModes );
    for ( int i = 0; i < numModes; i++ )
      Variances[i] = decompResult[0][n-i-1];
  }

  /**
   * \brief "Arg" contains all given data points (e.g. "aol::Vector" or "aol::MultiVector"),
   * "Dest" will contain the first modes of variation.
   */
  void apply ( const aol::VectorContainer<SampleType> &Arg, aol::VectorContainer<SampleType> &Dest ) const {
    aol::Vector<RealType> variances;
    modesAndVariances ( Arg, Dest, variances );
  }

  void applyAdd (  const aol::VectorContainer<SampleType> &, aol::VectorContainer<SampleType> & ) const {
    throw aol::Exception ( "ModesOfVariation::applyAdd(...) called.", __FILE__, __LINE__ );
  }

  //! \author Paetz
  void eigenvaluesAndEigenvectors ( const aol::MultiVector<RealType> &Samples, aol::MultiVector<RealType> &eigenvectors, aol::Vector<RealType> &eigenvalues, int numberOfModes ) const {
    // number of samples
    int n = Samples.size();
    aol::SymmetricMatrix<RealType> covMatrix ( Samples[0].size() );
    covMatrix.setZero();
    ProgressBar<> pb ( "Assemble covariance matrix " );
    pb.start ( n );
    for ( int k = 0; k < n; k++ ) {
      pb++;
      for ( int i = 0; i < Samples[0].size(); i++ ) {
        for ( int j = 0; j < Samples[0].size(); j++ ) {
          covMatrix.add ( i, j, Samples[k].get ( i ) *Samples[k].get ( j ) );
        }
      }
    }
    pb.finish();
    covMatrix *= 1.0 / ( n - 1.0 );

    // decompose covariance matrix (first row of decompResult contains eigenvalues in decreasing order, further rows contain corresponding eigenvectors)
    aol::MultiVector<RealType> decompResult;
    ( aol::DeflationEigenvectorOp<aol::SymmetricMatrix<RealType> > ( numberOfModes ) ).apply ( covMatrix, decompResult );
    eigenvectors.reallocate ( numberOfModes, Samples[0].size()  );

    for ( int i = 0; i < numberOfModes; i++ ) {
      eigenvectors[i] = decompResult[i+1];
    }

    // set eigenvalues
    eigenvalues.reallocate ( numberOfModes );
    for ( int i = 0; i < numberOfModes; i++ ) {
      eigenvalues[i] = decompResult[0][i];
      cerr << "Eigenvalue " << i << " : " << eigenvalues[i] << endl;
    }
  }


  void getDiagOfCovariance ( const aol::MultiVector<RealType> &Samples, aol::Vector<RealType> &diag ) const {
    int n = Samples.numComponents();
    diag.setZero();
    for ( int k = 0; k < n; k++ ) {
      for ( int i = 0; i < Samples[0].size(); i++ ) {
        diag.add ( i, Samples[k].get ( i ) *Samples[k].get ( i ) );
      }
    }
  }


  void getRowOfCovariance ( int rowNumber, const aol::MultiVector<RealType> &Samples, aol::Vector<RealType> &row ) const {
    int n = Samples.numComponents();
    row.setZero();
    for ( int k = 0; k < n; k++ ) {
      for ( int j = 0; j < Samples[0].size(); j++ ) {
        row.add ( j, Samples[k].get ( rowNumber ) *Samples[k].get ( j ) );
      }
    }
  }

//! return std::pair of index and value of maximal entry above index with possible permutation
  std::pair< int64_t, RealType > getMaxIndexAndValueCroppedWithPermutation ( Vector<RealType> &data, Vector<int> &permutation, int index, int &posInVec ) const {
    std::pair< int64_t, RealType > IndVal ( -1, 0.0 );
    for ( int64_t i = index; i < data.size(); ++i ) {
      if ( data.get ( permutation.get ( i ) ) >= IndVal.second ) {
        IndVal.first = permutation.get ( i );
        posInVec = i;
        IndVal.second = data.get ( permutation.get ( i ) );
      }
    }
    return ( IndVal );
  }




// This function is based on www.simtech.uni-stuttgart.de/publikationen/prints.php?ID=166
//! \author Paetz
  void eigenvaluesAndEigenvectorsPivotedCholesky ( const aol::MultiVector<RealType> &Samples, aol::MultiVector<RealType> &eigenvectors, aol::Vector<RealType> &eigenvalues, int numberOfModes ) const {
    int n = Samples[0].size();
    aol::Vector<RealType> diag ( n );
    aol::Vector<RealType> row ( n );
    aol::FullMatrix<RealType> lMatrix ( numberOfModes, Samples[0].size() );
    aol::Vector<int> permutation ( n );

    // Compute pivoted Cholesky decomposition
    getDiagOfCovariance ( Samples, diag );
    for ( int i = 0; i < n; i++ ) {
      permutation.set ( i, i );
    }
    for ( int m = 0; m < numberOfModes; m++ ) {
      std::pair<int, RealType> posAndValue;
      int posInVec = -1;
      posAndValue = getMaxIndexAndValueCroppedWithPermutation ( diag, permutation, m, posInVec );
      int temp = permutation.get ( m );
      permutation.set ( m, posAndValue.first );
      permutation.set ( posInVec, temp );
      lMatrix.set ( m, permutation.get ( m ), sqrt ( posAndValue.second ) );
      getRowOfCovariance ( permutation.get ( m ), Samples, row );
      for ( int i = m + 1; i < n; i++ ) {
        RealType value = row.get ( permutation.get ( i ) );
        for ( int j = 0; j < m; j++ ) {
          value -= lMatrix.get ( j, permutation.get ( m ) ) * lMatrix.get ( j, permutation.get ( i ) );
        }
        value /= lMatrix.get ( m, permutation.get ( m ) );
        lMatrix.set ( m, permutation.get ( i ), value );
        diag.add ( permutation.get ( i ), -lMatrix.get ( m, permutation.get ( i ) ) * lMatrix.get ( m, permutation.get ( i ) ) );
      }
    }

    // Compute low dimensional matrix with same leading eigenvalues
    aol::FullMatrix<RealType> smallMatrix ( numberOfModes, numberOfModes );
    for ( int i = 0; i < numberOfModes; i++ ) {
      for ( int j = 0; j < numberOfModes; j++ ) {
        for ( int k = 0; k < n; k++ ) {
          // low rank approximation
          smallMatrix.add ( i, j, lMatrix.get ( i, k ) * lMatrix.get ( j, k ) );
        }
      }
    }

    // Compute eigenvalues and eigenvectors of low dimensional matrix
    aol::MultiVector<RealType> decompResult;
    aol::JacobiEigenvectorOp<aol::FullMatrix<RealType> >  eigenOp;
    eigenOp.apply ( smallMatrix, decompResult );

    // decompose covariance matrix (first row of decompResult contains eigenvalues in decreasing order, further rows contain corresponding eigenvectors)

    eigenvectors.reallocate ( numberOfModes, Samples[0].size()  );

    for ( int i = 0; i < numberOfModes; i++ ) {
      eigenvectors[i].setZero();
      lMatrix.applyAddTranspose ( decompResult[i+1], eigenvectors[i] );
      eigenvectors[i] /= eigenvectors[i].norm();
    }

    // set eigenvalues
    eigenvalues.reallocate ( numberOfModes );
    for ( int i = 0; i < numberOfModes; i++ ) {
      eigenvalues[i] = decompResult[0][i];
    }
  }


};


/**
 *  Calculates the condition number of a matrix as largest/smallest eigenvalue by vector iteration and inverse vector iteration.
 *  Matrix should be spd, else solver will probably fail.
 *  \author Schwen
 */
template< typename MatrixType, typename RealType >
RealType computeConditionNumberByVectorIteration ( const MatrixType &mat, const RealType symmetryThreshold = aol::NumberTrait<RealType>::zero, const bool randomize = false ) {

  if ( ! ( aol::isMatrixSymmetric< MatrixType, RealType > ( mat, symmetryThreshold ) ) ) {
    throw aol::Exception ( "Cannot compute condition number of non-symmetric matrix by vector iteration and inverse vector iteration", __FILE__, __LINE__ );
  }

  // vector iteration and inverse vector iteration
  RealType lambdaMax = 0, lambdaMin = 0;

  aol::Vector<RealType> vOld ( mat.getNumRows() ), vNew ( mat.getNumRows() );
  aol::NoiseOperator<RealType> nO ( aol::NoiseOperator<RealType>::EQUALLY_DISTRIBUTED );

  if ( randomize ) {
    nO.randomize();
  }

  {
    nO.applySingle ( vNew );

    RealType lambdaNew = vNew.norm(), lambdaOld = 0;
    while ( fabs ( lambdaNew - lambdaOld ) / lambdaNew > 1.e-12 ) {
      lambdaOld = lambdaNew;
      vOld = vNew;
      vOld /= lambdaOld;

      mat.apply ( vOld, vNew );
      lambdaNew = vNew.norm();
    }

    lambdaMax = lambdaNew;
#ifdef VERBOSE
    cerr << "Largest Eigenvalue = " << lambdaMax << endl;
#endif
  }

  {
    aol::DiagonalPreconditioner< aol::Vector<RealType> > prec ( mat );
    aol::PCGInverse< aol::Vector<RealType> > matInv ( mat, prec, 1.0e-16, 10000 );
    matInv.setStopping ( aol::STOPPING_RELATIVE_TO_INITIAL_RESIDUUM );
    matInv.setQuietMode ( true );

    nO.applySingle ( vOld );
    mat.apply ( vOld, vNew );

    RealType lambdaNew = vNew.norm(), lambdaOld = 0;
    while ( fabs ( lambdaNew - lambdaOld ) / lambdaNew > 1.e-12 ) {
      lambdaOld = lambdaNew;
      vOld = vNew;
      vOld /= lambdaOld;

      vNew.setZero();
      matInv.apply ( vOld, vNew );
      lambdaNew = vNew.norm();
      // cerr << lambdaNew << endl;
    }

    lambdaMin = 1. / lambdaNew;
#ifdef VERBOSE
    cerr << "Smallest Eigenvalue = " << lambdaMin << endl;
#endif
  }

  return ( lambdaMax / lambdaMin );
}

template <typename RealType, typename MatrixType>
class EigenLibraryInterfaceOp {
#ifdef USE_EXTERNAL_EIGEN
  typedef Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic> EigenMatrixType;
#endif
  
public:
  static void getEigenVectorsRealSymmetric ( MatrixType &Arg, aol::MultiVector<RealType> &EigenVectors, aol::Vector<RealType> &EigenValues ) {
#ifdef USE_EXTERNAL_EIGEN
    if ( Arg.getNumRows ( ) != Arg.getNumCols ( ) ) throw aol::Exception ( "Matrix not square! Expected real-valued symmetric matrix.", __FILE__, __LINE__ );
    const int numColRows = Arg.getNumCols ( );
    const int numEigenVecs = numColRows;
    const int numEvals = numColRows;
    EigenMatrixType A ( numColRows, numColRows );
    for ( int i=0; i<numColRows ; ++i )
      for ( int j=0; j<numColRows ; ++j )
        A ( i, j ) = Arg.get ( i, j );
    Eigen::SelfAdjointEigenSolver<EigenMatrixType> es ( A );
    EigenVectors.reallocate ( numEigenVecs, numColRows );
    for ( int i=0; i<numEigenVecs ; ++i )
      for ( int j=0; j<numColRows ; ++j )
        EigenVectors[numEigenVecs-i-1][j] = ( es.eigenvectors ( ) ) ( i, j );
    EigenValues.reallocate ( numEvals );
    for ( int i=0; i<numEvals ; ++i )
      EigenValues[numEvals-i-1] = ( es.eigenvalues ( ) ) ( i );
#else
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( Arg );
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( EigenVectors );
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( EigenValues );
    throw aol::Exception ( "Eigen library required! Compile with -DUSE_EIGEN", __FILE__, __LINE__ );
#endif
  }
};
  
  
template <typename RealType>
void getPCAEigenVecs ( aol::MultiVector<RealType> &EigenVecs, aol::Vector<RealType> &EigenVals, const aol::MultiVector<RealType> &HighDimInput, const int NumEvals,
                       const RealType Omega = 0, const bool Quiet = false, const bool UseEigenLibrary = false, const bool UseRandomSubset = false ) {
  // Setup system matrix
  const int K = HighDimInput.numComponents ( );
  int N = HighDimInput[0].size ( );
  if ( UseRandomSubset ) N *= log ( static_cast<RealType> ( K ) ) / static_cast<RealType> ( K );
  
  // Compute sample mean from data
  aol::Vector<RealType> sampleMean ( K );
  HighDimInput.getMeanComponents ( sampleMean );
  
  // Assemble covariance matrix from N random samples
  aol::MultiVector<RealType> samples ( K, N );
  if ( UseRandomSubset ) {
    aol::RandomGenerator randomGenerator;
    aol::Vector<int> indices ( N );
    randomGenerator.rIntVecPairwiseDifferent ( indices, 0, HighDimInput[0].size ( ) );
    for ( int row=0; row<N ; ++row )
      for ( int col=0; col<K ; ++col )
        samples[col][row] = HighDimInput[col][indices[row]] - sampleMean[col];
  } else {
    for ( int row=0; row<N ; ++row )
      for ( int col=0; col<K ; ++col )
        samples[col][row] = HighDimInput[col][row] - sampleMean[col];
  }
  
  // Assemble AA^T
  aol::ProgressBar<> progressBar ( "Assembling ATA" );
  if ( !Quiet ) progressBar.start ( ceil ( K * ( K + 1 ) / 2 ) );
  aol::FullMatrix<RealType> AAT ( K, K );
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for ( int i=0; i<K ; ++i ) {
    for ( int j=0; j<=i ; ++j ) {
      AAT.set ( i, j, samples[i].dotProduct ( samples[j] ) );
      if ( !Quiet ) progressBar++;
    }
  }
  if ( !Quiet ) progressBar.finish ( );
  for ( int i=0; i<K ; ++i )
    for ( int j=i+1 ; j<K ; ++j )
      AAT.set ( i, j, AAT.get ( j, i ) );
  
  // Compute eigen vectors
  if ( UseEigenLibrary ) {
    aol::EigenLibraryInterfaceOp<RealType, aol::FullMatrix<RealType> >::getEigenVectorsRealSymmetric ( AAT, EigenVecs, EigenVals );
    
    // Determine number of eigen values / vectors to be used (if necessary)
    if ( NumEvals > 0 ) EigenVecs.resize ( NumEvals, K );
    else {
      int numEvals = EigenVals.size ( );
      RealType lse = 0;
      while ( numEvals > 1 && lse < Omega * N ) {
        --numEvals;
        lse += EigenVals[numEvals];
      }
      EigenVecs.resize ( numEvals+1, K );
      
      std::cerr << lse << std::endl;
      std::cerr << N << std::endl;
      
      std::cerr << "PCAMSSegmentor: detected " << numEvals + 1 << " segments." << std::endl;
    }
  } else {
    aol::MultiVector<RealType> eig;
    aol::DeflationEigenvectorOp<aol::FullMatrix<RealType> > evalOp ( NumEvals, Quiet );
    RealType sumSingularValsSqr = 0;
    for ( int row = 0; row < N; ++row )
      for ( int col = 0; col < K; ++col )
        sumSingularValsSqr += aol::Sqr<RealType> ( samples[col][row] );
    if ( NumEvals == 0 )
      evalOp.setThreshold ( Omega * N, sumSingularValsSqr );
    evalOp.apply ( AAT, eig );
    EigenVecs.resize ( eig.numComponents ( ) - 1, K );
    for ( int ev = 0; ev<eig.numComponents ( ) - 1; ++ev )
      for ( int i = 0; i<K; ++i )
        EigenVecs[ev][i] = eig[ev + 1][i];
    EigenVals.resize ( eig.numComponents ( ) ); // DeflationEigenvectorOp does not return eigenvalues
    EigenVals[0] = sumSingularValsSqr;
    for ( int ev = 1; ev < eig.numComponents ( ); ++ev )
      EigenVals[ev] = eig[0][ev];
  }
}
  
template <typename RealType>
void getPCAProjectedCoefficients ( aol::MultiVector<RealType> &LowDimOutput, const aol::MultiVector<RealType> &HighDimInput, const aol::MultiVector<RealType> &EigenVecs,
                                   const bool Quiet = false ) {
  const int K = HighDimInput.numComponents ( );
  
  // Compute sample mean from data
  aol::Vector<RealType> sampleMean ( K );
  HighDimInput.getMeanComponents ( sampleMean );
  
  // Project mean-centralized high-dimensional input onto specified eigen-vector basis
  aol::ProgressBar<> progressBar ( "Calculating projected coefficients" );
  if ( !Quiet ) progressBar.start ( HighDimInput[0].size ( ) );
  LowDimOutput.reallocate ( EigenVecs.numComponents ( ), HighDimInput[0].size ( ) );
  aol::Vector<RealType> highDimDatum ( HighDimInput.numComponents ( ) );
  for ( int k=0; k<HighDimInput[0].size ( ) ; ++k ) {
    HighDimInput.getTo ( k, highDimDatum );
    highDimDatum -= sampleMean;
    for ( int ev=0; ev<EigenVecs.numComponents ( ) ; ++ev )
      LowDimOutput[ev][k] = highDimDatum.dotProduct ( EigenVecs[ev] );
    if ( !Quiet ) progressBar++;
  }
  if ( !Quiet ) progressBar.finish ( );
}


} // end of namespace aol.

#endif
