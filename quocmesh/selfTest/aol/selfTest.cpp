// include all aol header files; in alphabetical order
#include <adaptiveTriangMesh.h>
#include <aol.h>
#include <ArmijoSearch.h>
#include <bandMatrix.h>
#include <baseFunctionSet.h>
#include <bitVector.h>
#include <boostIncludes.h>
#include <bzipiostream.h>
#include <ChanVese.h>
#include <clustering.h>
#include <complexUtils.h>
#include <convergenceEstimator.h>
#include <CSRMatrix.h>
#include <ctrlCCatcher.h>
#include <derivativeFreeOptimization.h>
#include <diagBandMatrix.h>
#include <discreteFunction.h>
#include <eigenvectors.h>
#include <eigenWrapper.h>
#include <epswriter.h>
#include <FEOpInterface.h>
#include <functionSpaceProjection.h>
#include <GaussSeidel.h>
#include <geom.h>
#include <gnuplotter.h>
#include <gradientDescent.h>
#include <gradientflow.h>
#include <hashSet.h>
#include <indexTranslationMatrix.h>
#include <interruptableIterativeInfo.h>
#include <iterativeInfo.h>
#include <lookupMap.h>
#include <maskedOp.h>
#include <maskedVector.h>
#include <matrix.h>
#include <matrixInverse.h>
#include <memoryManager.h>
#include <meshWithData.h>
#include <multiStreambuf.h>
#include <multiVector.h>
#include <Newton.h>
#include <NewtonInfo.h>
#include <op.h>
#include <parameterParser.h>
#include <platformDependent.h>
#include <pointerClasses.h>
#include <polarCoords.h>
#include <polynomial.h>
#include <preconditioner.h>
#include <probDistributionFunction.h>
#include <projectors.h>
#include <progressBar.h>
#include <qmException.h>
#include <QRDecomposition.h>
#include <randomGenerator.h>
#include <rgbColorMap.h>
#include <ringBuffer.h>
#include <ringBufferFrontInsertCapped.h>
#include <rows.h>
#include <RudinOsherFatemi.h>
#include <simultaneouslyAssembled.h>
#include <smallMat.h>
#include <smallVec.h>
#include <solver.h>
#include <solverInfo.h>
#include <sparseMatrices.h>
#include <sparseMatrixRowIterator.h>
#include <stats.h>
#include <suiteSparseSolver.h>
#include <tensor.h>
#include <timeStepping.h>
#include <timestepSaver.h>
#include <triangle.h>
#include <triangMesh.h>
#include <triangMeshConfigurators.h>
#include <trustRegionMethod.h>
#include <vec.h>
#include <vectorExtensions.h>
#include <verySparseMatrix.h>

// this quoc-header needs to be included to provide the gridDefinition
#include <gridBase.h>
#include <configurators.h>
#include <fastUniformGridMatrix.h>

using namespace aol;

WARNING_OFF ( type-limits )

// VC++'s C4146 warnings are unwarranted here.
#if defined ( _MSC_VER )
#pragma warning ( push )
#pragma warning ( disable : 4146 )
#endif

template <typename Type> void testSqrCubQrt ( Type (*powerfunc) (Type), int exponent, const std::string type, const Type a, const bool overflowExpected, const bool verbose = false ) {

  bool exception = false;
  Type result = 0;

  try  { result = (*powerfunc) (a); }
  catch ( aol::Exception& e ) { e.consume(); exception = true; };

  if ( verbose ) { cout << "(" << type << ") " << a << " ^ " << exponent << " = " << result; }
  if ( verbose ) { cout << "   -->   "; }
  if ( overflowExpected ) {
    if ( !exception ) {
      if ( verbose ) { cout << aol::color::red << "overflow expected but not found" << aol::color::reset << endl; }
      throw aol::Exception ("testSqrCubQrt", __FILE__, __LINE__);
    }
    else {
      if ( verbose ) { cout << "overflow as expected" << endl; }
    }
  }
  else {
    if ( !exception ) {
      if ( verbose ) { cout << "fits as expected" << endl; }
    }
    else {
      if ( verbose ) { cout << aol::color::red << "unexpected overflow" << aol::color::reset << endl; }
      throw aol::Exception ("testSqrCubQrt", __FILE__, __LINE__);
    }
  }
}

template <typename Type> void testSqrCubQrt ( Type (*powerfunc) (Type), int exponent, const std::string type ) {
  const Type max = std::numeric_limits<Type>::max();

  const bool inte = std::numeric_limits<Type>::is_integer;
  const bool sign = std::numeric_limits<Type>::is_signed;

  const Type root = pow (static_cast<long double> (max), static_cast<long double> (1.0/exponent));

  testSqrCubQrt<Type> (powerfunc, exponent, type, root / 2, false );
  testSqrCubQrt<Type> (powerfunc, exponent, type, root * 2, inte );
  if (sign) testSqrCubQrt<Type> (powerfunc, exponent, type, - root / 2, false );
  if (sign) testSqrCubQrt<Type> (powerfunc, exponent, type, - root * 2, inte );
}

template <typename Type> void testSqrCubQrt ( const std::string type ) {
  testSqrCubQrt<Type> ( aol::Sqr, 2, type );
  testSqrCubQrt<Type> ( aol::Cub, 3, type );
  testSqrCubQrt<Type> ( aol::Qrt, 4, type );
}

void testSqrCubQrt () {
  testSqrCubQrt<unsigned short> ("unsigned short");
  testSqrCubQrt<signed short> ("signed short");
  testSqrCubQrt<unsigned int> ("unsigned int");
  testSqrCubQrt<signed int> ("signed int");
  testSqrCubQrt<uint64_t> ("uint64_t");
  testSqrCubQrt<int64_t> ("int64_t");

  testSqrCubQrt<float> ("float");
  testSqrCubQrt<double> ("double");
  testSqrCubQrt<long double> ("long double");
}

template <typename Type> void testSumWillFit ( const std::string type, const Type a, const Type b, const bool overflowExpected, const bool verbose = false ) {

  if ( verbose ) { cout << "(" << type << ") " << a << " + " << b << " = " << static_cast<Type> (a + b); }
  if ( verbose ) { cout << "   -->   "; }
  if ( overflowExpected ) {
    if ( sumWillFit<Type> (a, b) ) {
      if ( verbose ) { cout << aol::color::red << "overflow expected but not found" << aol::color::reset << endl; }
      throw aol::Exception ("testSumWillFit", __FILE__, __LINE__);
    }
    else {
      if ( verbose ) { cout << "overflow as expected" << endl; }
    }
  }
  else {
    if ( sumWillFit<Type> (a, b) ) {
      if ( verbose ) { cout << "fits as expected" << endl; }
    }
    else {
      if ( verbose ) { cout << aol::color::red << "unexpected overflow" << aol::color::reset << endl; }
      throw aol::Exception ("testSumWillFit", __FILE__, __LINE__);
    }
  }
}

template <typename Type> void testProductWillFit ( const std::string type, const Type a, const Type b, const bool overflowExpected, const bool verbose = false ) {

  if ( verbose ) { cout << "(" << type << ") " << a << " * " << b << " = " << static_cast<Type> (a * b); }
  if ( verbose ) { cout << "   -->   "; }
  if ( overflowExpected ) {
    if ( productWillFit<Type> (a, b) ) {
      if ( verbose ) { cout << aol::color::red << "overflow expected but not found" << aol::color::reset << endl; }
      throw aol::Exception ("testProductWillFit", __FILE__, __LINE__);
    }
    else {
      if ( verbose ) { cout << "overflow as expected" << endl; }
    }
  }
  else {
    if ( productWillFit<Type> (a, b) ) {
      if ( verbose ) { cout << "fits as expected" << endl; }
    }
    else {
      if ( verbose ) { cout << aol::color::red << "unexpected overflow" << aol::color::reset << endl; }
      throw aol::Exception ("testProductWillFit", __FILE__, __LINE__);
    }
  }
}

template <typename Type> void testSumWillFit ( const std::string type ) {
  const Type max = std::numeric_limits<Type>::max();
  const Type min = std::numeric_limits<Type>::is_integer ?
    std::numeric_limits<Type>::min() : -std::numeric_limits<Type>::max();

  const bool sign = std::numeric_limits<Type>::is_signed;

  const Type offseta = static_cast<Type> ( std::numeric_limits<Type>::is_integer ? 42 : max * 0.1 * (1 - 5*std::numeric_limits<Type>::epsilon()) );
  const Type offsetb = static_cast<Type> ( std::numeric_limits<Type>::is_integer ? 43 : max * 0.1 * (1 + 5*std::numeric_limits<Type>::epsilon()) );
  // 3 is smallest factor that works

  testSumWillFit<Type> ( type, max - offseta, offseta, false );
  testSumWillFit<Type> ( type, max - offseta, offsetb, true  );
  testSumWillFit<Type> ( type, offseta, max - offseta, false );
  testSumWillFit<Type> ( type, offsetb, max - offseta, true  );
  if (sign) testSumWillFit<Type> ( type, min + offseta, -offseta, false );
  if (sign) testSumWillFit<Type> ( type, min + offseta, -offsetb, true  );
  if (sign) testSumWillFit<Type> ( type, -offseta, min + offseta, false );
  if (sign) testSumWillFit<Type> ( type, -offsetb, min + offseta, true  );
}

template <typename Type> void testProductWillFit ( const std::string type ) {
  const Type max = std::numeric_limits<Type>::max();

  const bool sign = std::numeric_limits<Type>::is_signed;

  Type smallfac = static_cast<Type> ( floor ( sqrt ( static_cast<long double> ( max ) ) ) );
  Type largefac = static_cast<Type> ( ceil  ( sqrt ( static_cast<long double> ( max ) ) ) );

  if ( !std::numeric_limits<Type>::is_integer) {
    smallfac *= (1 - 2*std::numeric_limits<Type>::epsilon());
    largefac *= (1 + 2*std::numeric_limits<Type>::epsilon());
    // even works with 1 -/+ epsilon
  }

  testProductWillFit<Type> ( type, smallfac, smallfac, false );
  testProductWillFit<Type> ( type, largefac, largefac, true  );
  if (sign) testProductWillFit<Type> ( type, -smallfac,  smallfac, false );
  if (sign) testProductWillFit<Type> ( type, -largefac,  largefac, true  );
  if (sign) testProductWillFit<Type> ( type,  smallfac, -smallfac, false );
  if (sign) testProductWillFit<Type> ( type,  largefac, -largefac, true  );
  if (sign) testProductWillFit<Type> ( type, -smallfac, -smallfac, false );
  if (sign) testProductWillFit<Type> ( type, -largefac, -largefac, true  );
}

// Turn on the warning C4146 again.
#if defined ( _MSC_VER )
#pragma warning ( pop )
#endif

void testSumAndProductWillFit () {

  testSumWillFit<unsigned short> ("unsigned short");
  testSumWillFit<signed short> ("signed short");
  testSumWillFit<unsigned int> ("unsigned int");
  testSumWillFit<signed int> ("signed int");
  testSumWillFit<uint64_t> ("uint64_t");
  testSumWillFit<int64_t> ("int64_t");

  testProductWillFit<unsigned short> ("unsigned short");
  testProductWillFit<signed short> ("signed short");
  testProductWillFit<unsigned int> ("unsigned int");
  testProductWillFit<signed int> ("signed int");
#ifndef _MSC_VER
  testProductWillFit<uint64_t> ("uint64_t");
#endif
  testProductWillFit<int64_t> ("int64_t");

  testSumWillFit<signed short> ("signed short", 20000, 20000, true);
  testProductWillFit<signed short> ("signed short", 256, 128, true);
}

//! Tests random number distribution, taking n ints in [0,k)
double chiSquare ( int n, int k, unsigned int seed ) {
  aol::Vector<int> count ( k );
  aol::RandomGenerator rg( seed );
  for ( int i = 0; i < n; ++i )
    ++ count [ rg.rInt ( 0, k ) ];

  double prob = 1.0 / k; double expected = prob * n;
  double chisquare = 0;
  for (int i = 0; i < k; ++i )
    chisquare += ( count [i] - expected ) * ( count [i] - expected ) / expected;

  return chisquare;
}

aol::Vec2<double> chiSquareNinetyPercentInterval ( int k ) {
  // Five percent fall below lower bound, five percent above upper bound
  // Approximation from Abramowitz, Handbook of Mathematical Functions, 26.4.17.
  // Checked against tabulated values for k = 30 and 100, against octave chi2inv for k = 1000
  double xp = 1.64; // Gaussian distribution is 0.95
  return aol::Vec2<double> (k * pow (1 - 2./(9*k) - xp * sqrt (2./(9*k)), 3),
          k * pow (1 - 2./(9*k) + xp * sqrt (2./(9*k)), 3));
}

double chiSquareInNinetyPercentInterval ( int n, int k, unsigned int firstseed, int repeat ) {
  int count = 0;

  for ( int i = 0; i < repeat; ++i ) {
    aol::Vec2<double> interval = chiSquareNinetyPercentInterval ( k );
    double chisquare = chiSquare ( n, k, firstseed + i );
    if ( chisquare > interval [0] && chisquare < interval [1] )
      ++ count;
  }
  return 100.0 * count / repeat;
}

template<typename RealType>
bool VectorOperatorStarTest(){
  bool testFailed = false;
  const int n = 155;
  aol::Vector<RealType> vec(n);
  for( int i = 0; i < n; i++){
    vec[i] = i+aol::ZOTrait<RealType>::one;
  }
  const RealType sum = static_cast<RealType>(n*(n+1)*(2*n+1))/6.;
  if( vec.normSqr() != sum )
    testFailed = true;
  return testFailed;
}

template <typename Operator, typename MatrixType>
bool testEigenvectorOp(const MatrixType & M, string opName)
{
  typedef typename MatrixType::DataType RealType;

  cerr << "Testing " << opName << "...";
  aol::MultiVector<RealType> ev;

  Operator evOp;

  evOp.apply(M, ev);

  int n = M.getNumRows();
  for (int i = 0; i < n; ++i) {
    aol::Vector<RealType> lambda_v (ev[i + 1]);
    lambda_v *= ev[0][i];
    aol::Vector<RealType> M_v(n);
    M.apply(ev[i+1], M_v);
    lambda_v -= M_v;
    RealType error = lambda_v.norm();
    if (error > 5E-13) {
      cerr << aol::color::red << opName << " failed, error " << error << aol::color::reset << endl;
      return false;
    }
  }
  cerr << aol::color::green << " passed." << aol::color::reset << endl;
  return true;
}

bool isappok ()
{
  bool appok = true;
  double reps = 2 * numeric_limits<double>::epsilon();
  double aeps = 5E-16;

  if ( ! aol::appleqRelative<double> ( -1,      2 ) ) appok = false;
  if ( ! aol::appleqRelative<double> (  1,      2 ) ) appok = false;
  if ( ! aol::appleqRelative<double> (  2-reps, 2 ) ) appok = false;
  if ( ! aol::appleqRelative<double> (  2,      2 ) ) appok = false;
  if ( ! aol::appleqRelative<double> (  2+reps, 2 ) ) appok = false;
  if (   aol::appleqRelative<double> (  3,      2 ) ) appok = false;

  if (   aol::appleqRelative<double> (  1,      -2 ) ) appok = false;
  if (   aol::appleqRelative<double> ( -1,      -2 ) ) appok = false;
  if ( ! aol::appleqRelative<double> ( -2+reps, -2 ) ) appok = false;
  if ( ! aol::appleqRelative<double> ( -2,      -2 ) ) appok = false;
  if ( ! aol::appleqRelative<double> ( -2-reps, -2 ) ) appok = false;
  if ( ! aol::appleqRelative<double> ( -3,      -2 ) ) appok = false;

  if ( ! aol::appleqAbsolute<double> ( -1,      2 ) ) appok = false;
  if ( ! aol::appleqAbsolute<double> (  1,      2 ) ) appok = false;
  if ( ! aol::appleqAbsolute<double> (  2-aeps, 2 ) ) appok = false;
  if ( ! aol::appleqAbsolute<double> (  2,      2 ) ) appok = false;
  if ( ! aol::appleqAbsolute<double> (  2+aeps, 2 ) ) appok = false;
  if (   aol::appleqAbsolute<double> (  3,      2 ) ) appok = false;

  if (   aol::appleqAbsolute<double> (  1,      -2 ) ) appok = false;
  if (   aol::appleqAbsolute<double> ( -1,      -2 ) ) appok = false;
  if ( ! aol::appleqAbsolute<double> ( -2+aeps, -2 ) ) appok = false;
  if ( ! aol::appleqAbsolute<double> ( -2,      -2 ) ) appok = false;
  if ( ! aol::appleqAbsolute<double> ( -2-aeps, -2 ) ) appok = false;
  if ( ! aol::appleqAbsolute<double> ( -3,      -2 ) ) appok = false;


  if (   aol::appgeqRelative<double> ( -1,      2 ) ) appok = false;
  if (   aol::appgeqRelative<double> (  1,      2 ) ) appok = false;
  if ( ! aol::appgeqRelative<double> (  2-reps, 2 ) ) appok = false;
  if ( ! aol::appgeqRelative<double> (  2,      2 ) ) appok = false;
  if ( ! aol::appgeqRelative<double> (  2+reps, 2 ) ) appok = false;
  if ( ! aol::appgeqRelative<double> (  3,      2 ) ) appok = false;

  if ( ! aol::appgeqRelative<double> (  1,      -2 ) ) appok = false;
  if ( ! aol::appgeqRelative<double> ( -1,      -2 ) ) appok = false;
  if ( ! aol::appgeqRelative<double> ( -2+reps, -2 ) ) appok = false;
  if ( ! aol::appgeqRelative<double> ( -2,      -2 ) ) appok = false;
  if ( ! aol::appgeqRelative<double> ( -2-reps, -2 ) ) appok = false;
  if (   aol::appgeqRelative<double> ( -3,      -2 ) ) appok = false;

  if (   aol::appgeqAbsolute<double> ( -1,      2 ) ) appok = false;
  if (   aol::appgeqAbsolute<double> (  1,      2 ) ) appok = false;
  if ( ! aol::appgeqAbsolute<double> (  2-aeps, 2 ) ) appok = false;
  if ( ! aol::appgeqAbsolute<double> (  2,      2 ) ) appok = false;
  if ( ! aol::appgeqAbsolute<double> (  2+aeps, 2 ) ) appok = false;
  if ( ! aol::appgeqAbsolute<double> (  3,      2 ) ) appok = false;

  if ( ! aol::appgeqAbsolute<double> (  1,      -2 ) ) appok = false;
  if ( ! aol::appgeqAbsolute<double> ( -1,      -2 ) ) appok = false;
  if ( ! aol::appgeqAbsolute<double> ( -2+aeps, -2 ) ) appok = false;
  if ( ! aol::appgeqAbsolute<double> ( -2,      -2 ) ) appok = false;
  if ( ! aol::appgeqAbsolute<double> ( -2-aeps, -2 ) ) appok = false;
  if (   aol::appgeqAbsolute<double> ( -3,      -2 ) ) appok = false;


  if (   aol::appeqRelative<double> ( -1,      2 ) ) appok = false;
  if (   aol::appeqRelative<double> (  1,      2 ) ) appok = false;
  if ( ! aol::appeqRelative<double> (  2-reps, 2 ) ) appok = false;
  if ( ! aol::appeqRelative<double> (  2,      2 ) ) appok = false;
  if ( ! aol::appeqRelative<double> (  2+reps, 2 ) ) appok = false;
  if (   aol::appeqRelative<double> (  3,      2 ) ) appok = false;

  if (   aol::appeqRelative<double> (  1,      -2 ) ) appok = false;
  if (   aol::appeqRelative<double> ( -1,      -2 ) ) appok = false;
  if ( ! aol::appeqRelative<double> ( -2+reps, -2 ) ) appok = false;
  if ( ! aol::appeqRelative<double> ( -2,      -2 ) ) appok = false;
  if ( ! aol::appeqRelative<double> ( -2-reps, -2 ) ) appok = false;
  if (   aol::appeqRelative<double> ( -3,      -2 ) ) appok = false;

  if (   aol::appeqAbsolute<double> ( -1,      2 ) ) appok = false;
  if (   aol::appeqAbsolute<double> (  1,      2 ) ) appok = false;
  if ( ! aol::appeqAbsolute<double> (  2-aeps, 2 ) ) appok = false;
  if ( ! aol::appeqAbsolute<double> (  2,      2 ) ) appok = false;
  if ( ! aol::appeqAbsolute<double> (  2+aeps, 2 ) ) appok = false;
  if (   aol::appeqAbsolute<double> (  3,      2 ) ) appok = false;

  if (   aol::appeqAbsolute<double> (  1,      -2 ) ) appok = false;
  if (   aol::appeqAbsolute<double> ( -1,      -2 ) ) appok = false;
  if ( ! aol::appeqAbsolute<double> ( -2+aeps, -2 ) ) appok = false;
  if ( ! aol::appeqAbsolute<double> ( -2,      -2 ) ) appok = false;
  if ( ! aol::appeqAbsolute<double> ( -2-aeps, -2 ) ) appok = false;
  if (   aol::appeqAbsolute<double> ( -3,      -2 ) ) appok = false;
  return appok;
}

/**
 * \brief This class implements the energy
 * \f$\sum_{i=1}^N\left(x_3exp\left(-\frac{(t_i-x_1)^2}{x_2}\right)-y_i\right)^2\f$,
 * where the \f$t_i\f$ and \f$y_i\f$ are passed to the constructor and \f$(x_i)_i\f$ are passed to "apply".
 */
class TestFittingEnergy
  : public aol::Op< aol::Vector<double>, aol::Scalar<double> > {
private:
  const aol::Vector<double> &_samplingNodes;
  const aol::Vector<double> &_samplingValues;
public:
  TestFittingEnergy( const aol::Vector<double> &SamplingNodes,
                     const aol::Vector<double> &SamplingValues ) :
    _samplingNodes( SamplingNodes ),
    _samplingValues( SamplingValues ) {}

  void applyAdd( const aol::Vector<double> &Arg, aol::Scalar<double> &Dest ) const {
    for ( int i = 0; i < _samplingNodes.size(); i++ )
      Dest += aol::Sqr( Arg[2] * exp( - aol::Sqr( _samplingNodes[i] - Arg[0] ) / Arg[1] ) - _samplingValues[i] );
  }
};

/**
 * \brief This class implements the derivative of the energy
 * \f$\sum_{i=1}^N\left(x_3exp\left(-\frac{(t_i-x_1)^2}{x_2}\right)-y_i\right)^2\f$ with respect to the \f$x_i\f$,
 * where the \f$t_i\f$ and \f$y_i\f$ are passed to the constructor and (x_i)_i are passed to "apply".
 */
class TestFittingEnergyVariation
  : public aol::Op< aol::Vector<double>, aol::Vector<double> > {
private:
  const aol::Vector<double> &_samplingNodes;
  const aol::Vector<double> &_samplingValues;
public:
  TestFittingEnergyVariation( const aol::Vector<double> &SamplingNodes,
                              const aol::Vector<double> &SamplingValues ) :
    _samplingNodes( SamplingNodes ),
    _samplingValues( SamplingValues ) {}

  void applyAdd( const aol::Vector<double> &Arg, aol::Vector<double> &Dest ) const {
    aol::Vector<double> aux( 3 );
    for ( int i = 0; i < _samplingNodes.size(); i++ ) {
      double diff = _samplingNodes[i] - Arg[0];
      double pot = exp( - aol::Sqr( diff ) / Arg[1] );
      aux[0] = Arg[2] * 2 * diff / Arg[1];
      aux[1] = Arg[2] * aol::Sqr( diff / Arg[1] );
      aux[2] = 1;
      Dest.addMultiple( aux, 2 * pot * ( Arg[2] * pot - _samplingValues[i] ) );
    }
  }
};

/**
 * \brief This class implements the second derivative of the energy
 * \f$\sum_{i=1}^N\left(x_3exp\left(-\frac{(t_i-x_1)^2}{x_2}\right)-y_i\right)^2\f$ with respect to the \f$x_i\f$,
 * where the \f$t_i\f$ and \f$y_i\f$ are passed to the constructor and (x_i)_i are passed to "apply".
 */
class TestFittingEnergySecondVariation
  : public aol::Op< aol::Vector<double>, aol::FullMatrix<double> > {
private:
  const aol::Vector<double> &_samplingNodes;
  const aol::Vector<double> &_samplingValues;
public:
  TestFittingEnergySecondVariation( const aol::Vector<double> &SamplingNodes,
                                    const aol::Vector<double> &SamplingValues ) :
    _samplingNodes( SamplingNodes ),
    _samplingValues( SamplingValues ) {}

  void applyAdd( const aol::Vector<double> &Arg, aol::FullMatrix<double> &Dest ) const {
    aol::Vector<double> auxV( 3 );
    aol::FullMatrix<double> auxM( 3, 3 );
    for ( int i = 0; i < _samplingNodes.size(); i++ ) {
      aol::Vec<2,int> index;
      double diff = _samplingNodes[i] - Arg[0];
      double pot = exp( - aol::Sqr( diff ) / Arg[1] );
      auxV[0] = Arg[2] * 2 * diff / Arg[1];
      auxV[1] = Arg[2] * aol::Sqr( diff / Arg[1] );
      auxV[2] = 1;
      for ( int j = 0; j < 3; j++ )
        for ( int k = 0; k < 3; k++ )
          auxM.set( j, k, 2 * auxV[j] * auxV[k] * aol::Sqr( pot ) );
      Dest += auxM;
      auxM.set( 0, 0, -2 * Arg[2] / Arg[1]  + auxV[0] * auxV[0] / Arg[2] );
      auxM.set( 0, 1, -auxV[1] / Arg[1]     + auxV[1] * auxV[0] / Arg[2] );
      auxM.set( 0, 2, auxV[0] / Arg[2] );
      auxM.set( 1, 1, -2 * auxV[1] / Arg[1] + auxV[1] * auxV[1] / Arg[2]);
      auxM.set( 1, 2, auxV[1] / Arg[2]);
      auxM.set( 2, 2, 0 );
      auxM.set( 1, 0, auxM.get( 0, 1 ) );
      auxM.set( 2, 0, auxM.get( 0, 2 ) );
      auxM.set( 2, 1, auxM.get( 1, 2 ) );
      auxM *= 2 * pot * ( Arg[2] * pot - _samplingValues[i] );
      Dest += auxM;
    }
  }
};

bool isQuocAssertOk () {
#ifdef DEBUG
  try {
    // in debug mode, an exception should be thrown
    QUOC_ASSERT ( 41.999 == 42 );
  }
  catch ( Exception & exc ) {
    if (exc.getMessage() == "Assertion failed: 41.999 == 42") {
      exc.consume();
      return true;
    }
  }
  return false;
#else
  // in optimized mode, no exception should be thrown
  QUOC_ASSERT ( 41.999 == 42 );
  return true;
#endif
}

template< typename DataType >
bool vectorSaveTest ( ) {
  aol::Vector<DataType> vec( 4 ), vecL, vecCL;
  vec[0] = static_cast<DataType> ( 23 );
  vec[2] = static_cast<DataType> ( -42 );

  vec.saveToFile ( "test.quo.bz2" ); // test compressed save&load
  vecL.loadFromFile ( "test.quo.bz2" );
  vecL -= vec;

  vecCL.template loadConvertFromFile< DataType > ( "test.quo.bz2" ); // template keyword is necessary here
  vecCL -= vec;

  aol::MultiVector<DataType> mVec( 2, 4 ), mVecL, mVecCL;
  mVec[0][1] = static_cast<DataType> ( -23 );
  mVec[0][3] = static_cast<DataType> (  -1 );
  mVec[1][1] = static_cast<DataType> (   1 );
  mVec[1][2] = static_cast<DataType> (  42 );

  mVec.saveToFile ( "test.quo" ); // also test uncompressed save&load
  mVecL.loadFromFile ( "test.quo" );
  mVecL -= mVec;

  mVecCL.template loadConvertFromFile<DataType> ( "test.quo" );
  mVecCL -= mVec;

  remove ( "test.quo" );
  remove ( "test.quo.bz2" );
  // could this fail due to floating point noise? need to change criterion to vec.norm() < 1.0e-16 or similar?
  return ( ( vecL.norm() == 0 ) && ( vecCL.norm() == 0 ) && ( mVecL.norm() == 0 ) && ( mVecCL.norm() == 0 ) );
}

template< typename DataType >
bool typeSizeTest ( int ExpectedSize, const string DataTypeName ) {
  cerr << "--- Testing ( sizeof ( " << DataTypeName << " ) == " << ExpectedSize << " ) ... ";
  const bool success = ( sizeof ( DataType ) == ExpectedSize );
  if ( success )
    cerr << "OK" << endl;
  else
    cerr << "FAILED" << endl;
  return success;
}

// Note: A define is the only way to use T as template and to convert it to a string.
#define TYPESIZETEST(T,N) typeSizeTest<T> ( N, #T )

// For the aol::TimeStepRKF45 test. MinGW doesn't seem to allow this to be defined inside main.
struct DummyReaction {
  void evaluateDerivative ( const aol::Vector<double> &,
                            const double ,
                            aol::Vector<double> & ) const {}
};

int main( int, char** ) {

  try {
    bool failed = false;

    failed = !( FullMatrix<double>::isSelfTestOk() );

#if (!defined(WIN32) && !defined(__APPLE__))
    {
      cerr << "--- Testing aol::memusage() ... "<< endl;

      bool memusagefailed = false;

      cerr << "Virtual memory:        " << aol::memusage (        VIRTUAL_MEMORY ) << endl;
      cerr << "Resident memory:       " << aol::memusage (       RESIDENT_MEMORY ) << endl;
      cerr << "Swap memory:           " << aol::memusage (           SWAP_MEMORY ) << endl;
      cerr << "Data memory:           " << aol::memusage (           DATA_MEMORY ) << endl;
      cerr << "Stack memory:          " << aol::memusage (          STACK_MEMORY ) << endl;
      cerr << "Max. Virtual memory:   " << aol::memusage (    MAX_VIRTUAL_MEMORY ) << endl;
      cerr << "Max. Resident memory:  " << aol::memusage (   MAX_RESIDENT_MEMORY ) << endl;
      cerr << "Memory manager memory: " << aol::memusage ( MEMORY_MANAGER_MEMORY ) << endl;
      cerr << "Memory from malloc():  " << aol::memusage (       MALLOCED_MEMORY ) << endl;
      cerr << "Memory from mmap():    " << aol::memusage (        MMAPPED_MEMORY ) << endl;

      int n = 1000000;

#ifdef DO_NOT_USE_MEMORYMANAGER
      if ( areWeUsingEFence() ) {
        cerr << "using efence malloc debugger, skipping aol::memusage() test" << endl;
      }
      else {
        { Vector<int> vec (n); }
        uint64_t mem1 = aol::memusage ();
        Vector<int> vec (n);
        uint64_t mem2 = aol::memusage ();
        uint64_t memdiff = mem2 - mem1;
        uint64_t minmemdiff = 0.9 * n * sizeof (int) - 9000, maxmemdiff = 1.1 * n * sizeof (int) + 9000;
        cerr << "Vector<int> (" << n << ") needs " << memdiff << " bytes, should be between " << minmemdiff << " and " << maxmemdiff << "." << endl;
        if ( memdiff < minmemdiff || memdiff > maxmemdiff ) memusagefailed = true;

        cerr << "Resident memory:       " << aol::memusage (       RESIDENT_MEMORY ) << endl;
        cerr << "Swap memory:           " << aol::memusage (           SWAP_MEMORY ) << endl;
      }
#else
      Vector<int>* vec = new Vector<int> (n);
      int64_t mem1 = aol::memusage ( MEMORY_MANAGER_MEMORY );
      delete vec;
      int64_t mem2 = aol::memusage ( MEMORY_MANAGER_MEMORY );
      int memdiff = mem2 - mem1;
      cerr << "Vector<int> (" << n << ") returned " << memdiff << " bytes to vector manager, should be " << n * sizeof (int) << "." << endl;
      if ( memdiff != static_cast<int>( n * sizeof (int) ) ) memusagefailed = true;
#endif

      if ( memusagefailed ) {
        failed = true;
        cerr << "--- Testing aol::memusage() ... FAILED." << endl;
      }
      else
        cerr << "--- Testing aol::memusage() ... OK." << endl;
    }
#endif

    { // solver test
      int size = 33;
      cerr << "--- Testing solvers ... (norm of solution - approx solution < 1e-8) " << endl;
      aol::Vector<double> orig(size), prod(size), soln(size);
      aol::SparseMatrix<double> M(size,size), Mt(size,size);

      for(int i = 0; i < size; i++){
        orig.set(i, 1.5 * i); // fill with some data
      }

      // set M, do not use FE Op
      M.set(0, 0, 4.0);
      M.set(0, 1, 1.0);
      for(int i = 1; i < (size - 1); i++){
        M.set(i, i-1, 1.0);
        M.set(i, i  , 4.0);
        M.set(i, i+1, 1.0);
      }
      M.set(size-1, size-2, 1.0);
      M.set(size-1, size-1, 4.0);

      M.transposeTo( Mt );

      M.apply(orig, prod);

      aol::DiagonalPreconditioner< aol::Vector<double> >D ( M );

      {
        aol::BiCGInverse< aol::Vector<double> >                                        bi_cg_solver( M, Mt );
        bi_cg_solver.setStopping ( aol::STOPPING_ABSOLUTE );
        bi_cg_solver.apply(prod, soln);
        soln -= orig;
        cerr << "BiCG:                         difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      {
        aol::CGInverse< aol::Vector<double> >                                          cg_solver( M );
        cg_solver.setStopping ( aol::STOPPING_ABSOLUTE );
        cg_solver.apply(prod, soln);
        soln -= orig;
        cerr << "CG:                           difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      {
        aol::JacobiInverse< aol::Vector<double>, aol::SparseMatrix<double> >           jacobi_solver( M );
        jacobi_solver.setStopping ( aol::STOPPING_ABSOLUTE );
        jacobi_solver.apply(prod, soln);
        soln -= orig;
        cerr << "Jacobi:                       difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      { // this is the default
        aol::GaussSeidelInverse< aol::Vector<double>, aol::SparseMatrix<double> >      gauss_seidel_solver( M );
        gauss_seidel_solver.setStopping ( aol::STOPPING_ABSOLUTE );
        gauss_seidel_solver.apply(prod, soln);
        soln -= orig;
        cerr << "Gauss-Seidel: (symmetric)     difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      {
        aol::GaussSeidelInverse< aol::Vector<double>, aol::SparseMatrix<double> >      gauss_seidel_solver( M );
        gauss_seidel_solver.setStopping ( aol::STOPPING_ABSOLUTE );
        gauss_seidel_solver.setGaussSeidelSweepingMode( aol::GAUSS_SEIDEL_FORWARD );
        gauss_seidel_solver.apply(prod, soln);
        soln -= orig;
        cerr << "Gauss-Seidel: (forward)       difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      {
        aol::GaussSeidelInverse< aol::Vector<double>, aol::SparseMatrix<double> >      gauss_seidel_solver( M );
        gauss_seidel_solver.setStopping ( aol::STOPPING_ABSOLUTE );
        gauss_seidel_solver.setGaussSeidelSweepingMode( aol::GAUSS_SEIDEL_RED_BLACK );
        gauss_seidel_solver.apply(prod, soln);
        soln -= orig;
        cerr << "Gauss-Seidel: (red-black)     difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      {
        aol::GMRESInverse< aol::Vector<double> >                                       gmres_solver( M );
        gmres_solver.apply(prod, soln);
        soln -= orig;
        cerr << "GMRES:                        difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      {
        aol::PGMRESInverse< aol::Vector<double> >                                      pgmres_solver( M, D, 1e-8, aol::STOPPING_ABSOLUTE, 10, 50 );
        pgmres_solver.apply( prod, soln );
        soln -= orig;
        cerr << "PGMRES:                        difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      {
        aol::PBiCGInverse< aol::Vector<double> >                                       p_bi_cg_solver( M, Mt, D, D );
        p_bi_cg_solver.setStopping ( aol::STOPPING_ABSOLUTE );
        p_bi_cg_solver.apply(prod, soln);
        soln -= orig;
        cerr << "PBiCG:                        difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      {
        aol::PBiCGStabInverse< aol::Vector<double> >                                   p_bi_cg_stab_solver( M, D );
        p_bi_cg_stab_solver.setStopping ( aol::STOPPING_ABSOLUTE );
        p_bi_cg_stab_solver.apply(prod, soln);
        soln -= orig;
        cerr << "PBiCGStab:                    difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      {
        aol::PCGInverse< aol::Vector<double> >                                         pcg_solver( M, D );
        pcg_solver.setStopping ( aol::STOPPING_ABSOLUTE );
        pcg_solver.apply(prod, soln);
        soln  -= orig;
        cerr << "PCG: (diagonal)               difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      {
        aol::GeometricScalingPreconditioner< aol::Vector<double> > GSP ( M );
        aol::PCGInverse< aol::Vector<double> >                                         pcg_solver( M, GSP );
        pcg_solver.setStopping ( aol::STOPPING_ABSOLUTE );
        pcg_solver.apply(prod, soln);
        soln  -= orig;
        cerr << "PCG: (geometric scaling)      difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      {
        aol::SSORPreconditioner<aol::Vector<double>, aol::SparseMatrix<double> > SSOR( M );
        aol::PCGInverse< aol::Vector<double> >                                         pcg_solver( M, SSOR );
        pcg_solver.setStopping ( aol::STOPPING_ABSOLUTE );
        pcg_solver.apply(prod, soln);
        soln  -= orig;
        cerr << "PCG: (SSOR)                   difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      {
        aol::GaussSeidelPreconditioner<aol::Vector<double>, aol::SparseMatrix<double> > GSPrec( M );
        aol::PCGInverse< aol::Vector<double> >                                         pcg_solver( M, GSPrec );
        pcg_solver.setStopping ( aol::STOPPING_ABSOLUTE );
        pcg_solver.apply(prod, soln);
        soln  -= orig;
        cerr << "PCG: (GaussSeidel)            difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      {
        aol::ILU0Preconditioner< double, aol::SparseMatrix<double> > ILU0( M );
        aol::PCGInverse< aol::Vector<double> >                                         pcg_solver( M, ILU0 );
        pcg_solver.setStopping ( aol::STOPPING_ABSOLUTE );
        pcg_solver.apply(prod, soln);
        soln  -= orig;
        cerr << "PCG: (ILU0)                   difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      {
        aol::FullMatrix < double > FM ( M );
        aol::LUInverse < double >                                                      lu_inverse ( M );
        lu_inverse.apply ( prod, soln );
        soln -= orig;
        cerr << "LU:                           difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      {
        aol::FullMatrix < double > FM ( M );
        aol::QRInverse < double >                                                      qr_inverse ( M );
        qr_inverse.apply ( prod, soln );
        soln -= orig;
        cerr << "QR:                           difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

      {
        aol::FullMatrix < double > FM ( M );
        aol::QRInversePivot < double >                                                 qr_inverse_pivot ( M );
        qr_inverse_pivot.apply ( prod, soln );
        soln -= orig;
        cerr << "pivoted QR:                   difference = " << soln.norm() << endl;
        failed = failed || ( soln.norm() > 1e-8 );
      }

#ifdef USE_EXTERNAL_SUITESPARSE
      {
        aol::BlockOp< double, aol::SparseMatrix<double> > BM( 1, 1 );
        aol::MultiVector< double > prodMV, solnMV;
        prodMV.appendReference( prod );
        solnMV.appendReference( soln );
        BM.setReference( 0, 0, M );

        {
          aol::CholeskyBlockInverseOp< double, aol::SparseMatrix<double> >               cholesky_solver( 3 );
          cholesky_solver.setMatrix ( BM );
          cholesky_solver.apply(prodMV, solnMV);
          soln -= orig;
          cerr << "Cholesky:                     difference = " << soln.norm() << endl;
          failed = failed || ( soln.norm() > 1e-8 );
        }

        {
          aol::UMFPACKBlockInverseOp< double, aol::SparseMatrix<double> >                umfpack_solver( BM );
          umfpack_solver.apply(prodMV, solnMV);
          soln -= orig;
          cerr << "UMFPACK-LU:                   difference = " << soln.norm() << endl;
          failed = failed || ( soln.norm() > 1e-8 );
        }
      }
#endif // USE_EXTERNAL_SUITESPARSE

      // missing: multigrid solver

      // missing: block solvers and preconditioners

      if( !failed )
        cerr << "Solver tests OK." << endl;
      else
        cerr << "Solver tests FAILED!" << endl;
    }

    { // iterative minimization test
      cerr << "--- Testing iterative minimizers ... (norm of solution - approx solution < 1e-7) " << endl;

      // define the sampling nodes and sampling values of a best fit problem
      double nodes[50] = { -1.9000000e+00, -1.8000000e+00, -1.7000000e+00, -1.6000000e+00, -1.5000000e+00, -1.4000000e+00,
                           -1.3000000e+00, -1.2000000e+00, -1.1000000e+00, -1.0000000e-00, -9.0000000e-01, -8.0000000e-01, -7.0000000e-01, -6.0000000e-01, -5.0000000e-01, -4.0000000e-01, -3.0000000e-01,
                           -2.0000000e-01, -1.0000000e-01, 2.2204460e-16, 1.0000000e-01, 2.0000000e-01, 3.0000000e-01, 4.0000000e-01, 5.0000000e-01, 6.0000000e-01, 7.0000000e-01, 8.0000000e-01,
                           9.0000000e-01, 1.0000000e+00, 1.1000000e+00, 1.2000000e+00, 1.3000000e+00, 1.4000000e+00, 1.5000000e+00, 1.6000000e+00, 1.7000000e+00, 1.8000000e+00, 1.9000000e+00,
                           2.0000000e+00, 2.1000000e+00, 2.2000000e+00, 2.3000000e+00, 2.4000000e+00, 2.5000000e+00, 2.6000000e+00, 2.7000000e+00, 2.8000000e+00, 2.9000000e+00, 3.0000000e+00 };

      double values[50] = { -4.0050688e-01, -7.8172281e-01, 5.2514201e-01, -3.7216753e-02, 3.5063169e-01, 3.7429662e-01,
                            3.6361322e-01, 4.3030301e-02, 1.7239420e-01, 5.2690474e-01, 3.4801780e-01, 5.6713233e-01, 1.7677224e+00, 1.4124840e+00, 1.1363354e+00, 1.8600646e+00, 1.6964492e+00,
                            1.8940503e+00, 2.6259315e+00, 3.4140395e+00, 3.3863418e+00, 3.8933322e+00, 3.4463949e+00, 3.9921090e+00, 4.3075085e+00, 4.1303209e+00, 4.2522126e+00, 5.2734804e+00,
                            4.9795604e+00, 4.7419417e+00, 5.2973539e+00, 4.9936438e+00, 4.3840835e+00, 5.1514160e+00, 4.5282853e+00, 4.7679179e+00, 4.3687339e+00, 3.3570898e+00, 2.8181096e+00,
                            3.0034828e+00, 2.5981326e+00, 2.0963102e+00, 2.3468947e+00, 2.4719517e+00, 1.4046720e+00, 1.0514832e+00, 1.0801958e+00, 1.2547032e+00, 4.8069333e-01, 1.9615049e-01 };

      // using FLAT_COPY here may result in problems when using SSE because the arrays defined above are not 16-byte aligned
      aol::Vector<double> samplingNodes( nodes, 50, aol::DEEP_COPY ), samplingValues( values, 50, aol::DEEP_COPY );

      // define the best fit problem
      TestFittingEnergy e( samplingNodes, samplingValues );
      TestFittingEnergyVariation de( samplingNodes, samplingValues );
      TestFittingEnergySecondVariation d2e( samplingNodes, samplingValues );

      // set the starting parameter values for the optimization as well as the exact result
      aol::Vector<double> initialX( 3 ), x( 3 ), solution( 3 );
      initialX[0] = 1.0;      initialX[1] = 2.0;      initialX[2] = 5.0;
      solution[0] = 1.028970500; solution[1] = 1.839322574; solution[2] = 5.046554164;

      // prepare output for the Newton methods
      aol::NewtonInfo<double> newtonInfo ( 1E-10, 100, 1E-20, 1000, aol::STOPPING_ABSOLUTE );
      newtonInfo.setQuietMode ( true );
      newtonInfo.getSolverInfo().setQuietMode ( true );
      //newtonInfo.setTablePrintout();

      { // gradient descent method
        aol::GridlessGradientDescent< double, aol::Vector<double> > gradientDescent( e, de, 200 );
        gradientDescent.apply( initialX, x );
        x -= solution;
        cerr << "gradient descent:                                  difference = " << x.norm() << endl;
        failed = failed || ( x.norm() > 1e-7 );

        x.setZero();
        gradientDescent.setConfigurationFlags ( aol::GridlessGradientDescent< double, aol::Vector<double> >::USE_NONLINEAR_CG );
        gradientDescent.apply( initialX, x );
        x -= solution;
        cerr << "gradient descent with nonlinear CG:                difference = " << x.norm() << endl;
        failed = failed || ( x.norm() > 1e-7 );

        x.setZero();
        gradientDescent.setConfigurationFlags ( 0 );
        gradientDescent.setTimestepController ( aol::GridlessGradientDescent< double, aol::Vector<double> >::POWELL_WOLFE );
        gradientDescent.apply( initialX, x );
        x -= solution;
        cerr << "gradient descent with Powel-Wolfe                  difference = " << x.norm() << endl;
        failed = failed || ( x.norm() > 1e-7 );
      }

      { // quasi Newton method (DFP)
        aol::QuasiNewtonIteration< double, aol::Vector<double>, aol::Vector<double> > quasiNewtonDescent( e, de, newtonInfo );
        quasiNewtonDescent.apply( initialX, x );
        x -= solution;
        cerr << "quasi Newton method:                               difference = " << x.norm() << endl;
        failed = failed || ( x.norm() > 1e-7 );
      }

      { // quasi Newton method (BFGS-Update of the Matrix)
        aol::QuasiNewtonBFGS< double, aol::Vector<double>, aol::Vector<double> > quasiNewtonBFGS ( e, de, newtonInfo );
        quasiNewtonBFGS.apply( initialX, x );
        x -= solution;
        cerr << "quasi Newton method (BFGS):                        difference = " << x.norm() << endl;
        failed = failed || ( x.norm() > 1e-7 );
      }

      { // Newton method with objective function linesearch
        aol::FullMatrix<double>* pMat = new aol::FullMatrix<double>( 3, 3 );
        aol::GMRESInverse< aol::Vector<double>, aol::Op<aol::Vector<double>,aol::Vector<double> > >* pSol = new aol::GMRESInverse< aol::Vector<double>, aol::Op<aol::Vector<double>,aol::Vector<double> > >( *pMat, newtonInfo.getSolverInfo() );
        aol::NewtonMinimizationBase< double, aol::Vector<double>, aol::Vector<double>, aol::FullMatrix<double> > newtonDescent( e, de, d2e, pMat, newtonInfo );
        newtonDescent.setSolverPointer( pSol );
        newtonDescent.setQuietMode( false );
        pSol->setQuietMode( true );
        newtonDescent.apply( initialX, x );
        x -= solution;
        cerr << "Newton method with objective function linesearch:  difference = " << x.norm() << endl;
        failed = failed || ( x.norm() > 1e-7 );
      }

      { // Newton method with objective function linesearch
        aol::FullMatrix<double>* pMat = new aol::FullMatrix<double>( 3, 3 );
        aol::GMRESInverse< aol::Vector<double>, aol::Op<aol::Vector<double>,aol::Vector<double> > >* pSol = new aol::GMRESInverse< aol::Vector<double>, aol::Op<aol::Vector<double>,aol::Vector<double> > >( *pMat, newtonInfo.getSolverInfo() );
        aol::NewtonIterationBase< double, aol::Vector<double>, aol::Vector<double>, aol::FullMatrix<double> > newtonIteration( pMat, de, d2e, newtonInfo );
        newtonIteration.setSolverPointer( pSol );
        newtonIteration.setQuietMode( false );
        pSol->setQuietMode( true );
        newtonIteration.apply( initialX, x );
        x -= solution;
        cerr << "Newton method with gradient norm linesearch:       difference = " << x.norm() << endl;
        failed = failed || ( x.norm() > 1e-7 );

        x.setZero();
        newtonIteration.setTimestepController( NewtonInfo<double>::NEWTON_OPTIMAL );
        newtonIteration.apply( initialX, x );
        x -= solution;
        cerr << "Newton method with NEWTON_OPTIMAL linesearch:      difference = " << x.norm() << endl;
        failed = failed || ( x.norm() > 1e-7 );
      }

#ifdef USE_EXTERNAL_SUITESPARSE
      { // Trust region Newton method with Cholesky factorization
        aol::FullMatrix<double>* pMat = new aol::FullMatrix<double>( 3, 3 );
        aol::TrustRegionMethod<double,aol::Vector<double>,aol::FullMatrix<double> > trustRegionMethod( e, de, d2e, pMat, 100, 1E-10 );
        trustRegionMethod.apply( initialX, x );
        x -= solution;
        cerr << "trust region method with Cholesky factorization:   difference = " << x.norm() << endl;
        failed = failed || ( x.norm() > 1e-7 );
      }

      { // Trust region Newton method with Lanczos iteration
        aol::FullMatrix<double>* pMat = new aol::FullMatrix<double>( 3, 3 );
        aol::TrustRegionMethod<double,aol::Vector<double>,aol::FullMatrix<double>,aol::DiagonalPreconditioner<aol::Vector<double> > > trustRegionMethod( e, de, d2e, pMat, 100, 1E-10 );
        trustRegionMethod.setIterativeMode();
        trustRegionMethod.apply( initialX, x );
        x -= solution;
        cerr << "trust region method with Lanczos iteration:        difference = " << x.norm() << endl;
        failed = failed || ( x.norm() > 1e-7 );
      }
#endif

      // missing: some specializations

      if( !failed )
        cerr << "Iterative minimization tests OK." << endl;
      else
        cerr << "Iterative minimization tests FAILED!" << endl;
    }

    {
      cerr << "--- Testing aol::Vector copy constructors and resize... ";

      aol::Vector<float> vec_orig(3);
      vec_orig.setZero();
      vec_orig[1] = 42.0;

      {
        aol::Vector<float> vec_default(vec_orig); // default behavior of copy constructor should be deep copy.
        vec_default[1] = 23.0;

        failed = failed || ( ( abs( vec_orig[1] - 42.0 ) + abs( vec_default[1] - 23.0) ) > 1e-6 );

        aol::Vector<float> vec_deep(vec_orig, aol::DEEP_COPY );
        vec_deep[1] = 23.0;

        failed = failed || ( ( abs( vec_orig[1] - 42.0 ) + abs( vec_deep[1] - 23.0) ) > 1e-6 );

        aol::Vector<float> vec_flat(vec_orig, aol::FLAT_COPY );
        vec_flat[1] = 23.0;

        failed = failed || ( ( abs( vec_orig[1] - 23.0 ) + abs( vec_flat[1] - 23.0) ) > 1e-6 );

        aol::Vector<float> vec_struct(vec_orig, aol::STRUCT_COPY ); // new Vector with same size
        vec_struct[1] = 17.0;

        failed = failed || ( ( abs( vec_orig[1] - 23.0 ) + abs( vec_struct[1] - 17.0) ) > 1e-6 );

      } // End of Scope - destroy all vector copies

      aol::MemoryManager::deleteUnlocked(); // delete retained vectors from vector manager

      failed = failed || ( abs( vec_orig[1] - 23.0 )  > 1e-6 );

      // now try different variants of resizing the vector and see whether new entries are set to zero correctly
      vec_orig.resize( 23 );
      failed = failed || ( abs ( vec_orig[15] - 0.0 ) > 1.0e-6 );
      vec_orig[15] = 3.1415;

      vec_orig.resize( 3 );
      vec_orig.resize( 20 ); // shorter than before
      failed = failed || ( abs ( vec_orig[15] - 0.0 ) > 1.0e-6 );
      vec_orig[15] = 3.1415;

      vec_orig.resize( 3 );
      vec_orig.resize( 30 ); // shorter than before
      failed = failed || ( abs ( vec_orig[15] - 0.0 ) > 1.0e-6 );

      vec_orig = *&vec_orig; // check whether self-assignment is treated correctly

      // check whether eraseFirstOccurence works properly
      aol::Vector<int> intVec ( 3 );
      for ( int i = 0; i < intVec.size(); ++i ) {
        intVec[i] = i;
      }
      intVec.eraseFirstOccurence ( 1 );
      failed |= ( intVec.crc32OfData() != 819259421 );

      if( !failed )
        cerr << "OK." << endl;
      else
        cerr << "FAILED!" << endl;
    }

    {
      cerr << "--- Testing aol::MultiVector<float> copy constructors ... ";

      aol::MultiVector<float> mvec_orig(3,3);
      mvec_orig.setZero();
      mvec_orig[1][1] = 42.0;

      {
        aol::MultiVector<float> mvec_default(mvec_orig); // default behavior of copy constructor should be deep copy.
        mvec_default[1][1] = 23.0;

        failed = failed || ( ( abs( mvec_orig[1][1] - 42.0 ) + abs( mvec_default[1][1] - 23.0) ) > 1e-6 );

        aol::MultiVector<float> mvec_deep(mvec_orig, aol::DEEP_COPY );
        mvec_deep[1][1] = 23.0;

        failed = failed || ( ( abs( mvec_orig[1][1] - 42.0 ) + abs( mvec_deep[1][1] - 23.0) ) > 1e-6 );

        aol::MultiVector<float> mvec_flat(mvec_orig, aol::FLAT_COPY );
        mvec_flat[1][1] = 23.0;

        failed = failed || ( ( abs( mvec_orig[1][1] - 23.0 ) + abs( mvec_flat[1][1] - 23.0) ) > 1e-6 );

        aol::MultiVector<float> mvec_struct( mvec_orig, aol::STRUCT_COPY ); // new multiVector with same structure
        mvec_struct[1][1] = 17.0;

        failed = failed || ( ( abs( mvec_orig[1][1] - 23.0 ) + abs( mvec_struct[1][1] - 17.0) ) > 1e-6 );
      }

      aol::MemoryManager::deleteUnlocked();

      failed = failed || ( abs( mvec_orig[1][1] - 23.0 ) > 1e-6 );

      mvec_orig = *&mvec_orig;

      if( !failed )
        cerr << "OK." << endl;
      else
        cerr << "FAILED!" << endl;

    }

    {
      cerr << "--- Testing aol::BitVector ... ";
      aol::BitVector bf(20);
      bf.setAll( true );
      failed = failed || ( bf.numTrue() != 20 );

      bf.resize ( 200 );
      failed = failed || ( bf.numTrue() != 20 );

      bf.resize ( 10 );
      failed = failed || ( bf.numTrue() != 10 );

      aol::BitVector bf2 ( bf );
      failed = failed || ( bf != bf2 );

      bf2.set ( 9, false );
      failed = failed || ( bf == bf2 );

      bf2.set ( 9, true );
      bf2.set ( 2, false );
      failed = failed || ( bf == bf2 );

      bf = *&bf;

      if( !failed )
        cerr << "OK." << endl;
      else
        cerr << "FAILED!" << endl;
    }

    {
      cerr << "--- Testing aol::Vector<float> operator * ... ";
      if( !VectorOperatorStarTest<float>() )
        cerr << "OK." << endl;
      else{
        cerr << "FAILED!" << endl;
        failed = true;
      }

      cerr << "--- Testing aol::Vector<double> operator * ... ";
      if( !VectorOperatorStarTest<double>() )
        cerr << "OK." << endl;
      else{
        cerr << "FAILED!" << endl;
        failed = true;
      }
    }

    {
      cerr << "--- Testing aol::SparseMatrix<float>::resize ... ";
      aol::SparseMatrix<float> smat( 3, 4);
      smat.resize( 3, 5 );
      smat.resize( 4, 5 );
      smat.resize( 3, 5 );
      // this is not implemented yet.      smat.resize( 3, 4 );
      cerr << "OK." << endl;
    }

    {
      cerr << "--- Testing aol::Matrix<double>::operator+=/operator-= ... ";
      aol::FullMatrix<double> M0 ( 5, 5 );
      aol::SparseMatrix<double> M1 ( 5, 5 );
      for ( int i = 0; i < 5; ++i ) {
        for ( int j = 0; j < 5; ++j ) {
          M0.set ( i, j, i*j );
          M1.set ( i, j, i*j );
        }
      }
      M0 += M1;
      M0 *= 0.5;
      M0 -= M1;

      if( M0.getFrobeniusNormSqr() < 1.0e-14 )
        cerr << "OK." << endl;
      else{
        cerr << "FAILED!" << endl;
        failed = true;
      }
    }

    {
      cerr << "--- Testing aol::SimultaneouslyAssembledLinCombOp ... ";
      qc::GridDefinition grid( 4, qc::QC_2D );      // 2d grid of level 4 (16^2 cells)
      int N = grid.getNumberOfNodes();
      typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature<double,qc::QC_2D,3> > ConfigType;
      ConfigType config( grid );
      double tau = grid.H();

      aol::MassOp<ConfigType> massOp( grid );
      aol::StiffOp<ConfigType> stiffOp( grid );

      // the new op
      aol::SimultaneouslyAssembledLinCombOp<ConfigType> SAOp( grid );
      SAOp.appendGetCoeffOpReference( massOp, aol::SCALED_MASS_INT );
      SAOp.appendGetCoeffOpReference( stiffOp, aol::SCALED_STIFF_INT, tau );

      // instantiate and rebuild the matrices that are fixed during the evolution ...
      qc::FastUniformGridMatrix< double, qc::QC_2D > matSimulAssembled( grid );    // matrix for the new op
      qc::FastUniformGridMatrix< double, qc::QC_2D > matUsualAssembled( grid );    // matrix for the usual assembly

      matUsualAssembled.setZero();
      stiffOp.assembleAddMatrix( matUsualAssembled );
      matUsualAssembled *= tau;
      massOp.assembleAddMatrix( matUsualAssembled );

      matSimulAssembled.setZero();
      SAOp.assembleAddMatrix( matSimulAssembled );

      bool comparisonSuccess = aol::compareOps<double>( matSimulAssembled, matUsualAssembled, N, N, 1e-16 );

      if ( comparisonSuccess )  cerr << "OK." << endl;
      else {
        failed = true;
        cerr << aol::color::red << "<---------------- SimultaneouslyAssembledLinCombOp failed!\n" << aol::color::reset;
      }
    }

    {
      // random generator should produce the same sequence of random data for a fixed seed!
      cerr << "--- Testing aol::RandomGenerator reproducibility ... " ;

      aol::RandomGenerator rg( 31415 ); // specify a random seed: this way, different but reproducible sequences of random data can be generated

      if ( rg.rInt( 23, 42 ) != 36 || fabs ( rg.rReal( 1.414, 3.141 ) - 1.75776 ) > 1e-3 || (rg.rBool() != true) ) {
        failed = true;
        cerr << aol::color::red << "<---------------- RandomGenerator failed!\n" << aol::color::reset;
      } else {
        cerr << "OK" << endl;
      }
    }

    {
      cerr << "--- Testing aol::RandomGenerator randomness ... " ;
      double percent = chiSquareInNinetyPercentInterval ( 10000, 1000, 8442431, 100 );
      if ( percent < 80 ) {
        failed = true;
        cerr << aol::color::red << "<---------------- RandomGenerator failed!\n" << aol::color::reset;
      } else {
        cerr << "OK" << endl;
      }
      cerr << "    " << percent << " % of chi^2 tests were in the 90 % interval." << endl;
    }

    {
      cerr << "--- Testing aol::solveRTransposedR<double> ... ";
      const int k = 10;
      aol::FullMatrix<double> R( k, k );
      aol::Vector<double> x( k ), rhs ( k ), temp ( k );

      // fill R and rhs with something.
      for ( int i = 0; i < R.getNumRows(); i++ ) {
        rhs[i] = 1. + i;
        for ( int j = 0; j <= i; j++ ) {
          R.set ( j, i, i+j+1. );
        }
      }

      // Solve R^tRx=rhs.
      aol::solveRTransposedR<double>( R, rhs, x );

      // Calculate the norm of (R^tRx-rhs)
      aol::FullMatrix<double> L (R);
      L.transpose();
      R.apply( x, temp );
      L.apply( temp, x );
      x -= rhs;

      if( x.norm() < 0.0001 )
        cerr << "OK." << endl;
      else{
        cerr << "FAILED!" << endl;
        failed = true;
      }
    }

    {
      cerr << "--- Testing aol::isNaN<double> ... ";
        double zero = 0.;
        bool isNaNok = true;
        isNaNok = isNaNok && ( aol::isNaN( zero/zero ) == true );
        isNaNok = isNaNok && ( aol::isNaN( aol::NumberTrait<float>::NaN ) == true );
        isNaNok = isNaNok && ( aol::isNaN( 1.0/zero ) == false );
        isNaNok = isNaNok && ( aol::isNaN( -11.0/zero ) == false );
        isNaNok = isNaNok && ( aol::isNaN( -11.0 ) == false );
        isNaNok = isNaNok && ( aol::isNaN( aol::NumberTrait<float>::Inf ) == false );

        if( isNaNok )
          cerr << "OK." << endl;
        else{
          cerr << "FAILED!" << endl;
          failed = true;
        }
    }

    {
      cerr << "--- Testing aol::Matrix22<int> transposition ... ";
      aol::Matrix22<int> mat22a, mat22b, mat22c, mat22d;
      mat22a[0][1] = 42;
      mat22b[1][0] = 42;
      mat22a.transposeTo ( mat22c );
      mat22d.transposeFrom ( mat22b );

      failed |= ! ( ( mat22b == mat22c ) && ( mat22a == mat22d ) );

      mat22d.transpose();
      mat22c = mat22a.transposed();

      failed |= ! ( ( mat22d == mat22b ) && ( mat22c == mat22b ) );

      if( !failed )
        cerr << "OK." << endl;
      else
        cerr << "FAILED!" << endl;
    }

    {
      cerr << "--- Testing aol::Matrix33Symm<double>::eigenVectors ... ";
      aol::Matrix33Symm<double> mat ( 1, 2, 3,
                                      2, 3, 4,
                                      3, 4, 5 );

      aol::Vec3<double> eigVals;
      aol::Matrix33<double> eigVecs;

      mat.eigenVectors ( eigVals, eigVecs );

      aol::Matrix33<double> diag ( eigVals [0], 0, 0,
                                   0, eigVals [1], 0,
                                   0, 0, eigVals [2] ), q = eigVecs, qt = q;
      qt.transpose ();

      q *= diag; q *= qt;

      aol::Matrix33<double> diff ( mat );
      diff -= q;
      diff *= diff;

      if( diff.tr () < 1e-8 )
        cerr << "OK." << endl;
      else{
        cerr << "FAILED!" << endl;
        failed = true;
      }
    }

    {
      cerr << "--- Testing aol::SparseBlockMatrix<aol::Matrix<double> >::maxNumNonZeroesPerRow ... ";
      aol::SparseBlockMatrix<aol::Matrix<double> > blockMat (2,2);
      bool firstTest = false;
      if ( blockMat.maxNumNonZeroesPerRow() == 0 )
        firstTest = true;

      qc::GridDefinition grid( 4, qc::QC_2D );      // 2d grid of level 4 (16^2 cells)
      aol::DiagonalMatrix<double> diagonalDummyMat ( grid );

      bool secondTest = false;
      aol::Matrix<double> &blockMat01 = blockMat.allocateMatrix ( 0, 1, grid, diagonalDummyMat );
      blockMat01.setAll(1);
      aol::Matrix<double> &blockMat10 = blockMat.allocateMatrix ( 1, 0, grid, diagonalDummyMat );
      blockMat10.setAll(1);
      aol::Matrix<double> &blockMat11 = blockMat.allocateMatrix ( 1, 1, grid, diagonalDummyMat );
      blockMat11.setAll(1);

      if ( blockMat.maxNumNonZeroesPerRow() == 2 )
        secondTest = true;

      if( firstTest && secondTest )
        cerr << "OK." << endl;
      else{
        cerr << "FAILED!" << endl;
        failed = true;
      }
    }

    // *** EigenvectorOp derived classes test ***
    aol::FullMatrix<double> ev_matrix ( 64, 64 );
    for ( int i = 0; i < 64; ++i ) {
      ev_matrix.set ( i, i, -4 );
      if ( ( i % 8 ) != 0 ) ev_matrix.set ( i, i - 1, 1. );
      if ( ( i % 8 ) != 8 - 1 ) ev_matrix.set ( i, i + 1, 1. );
      if ( i >= 8 ) ev_matrix.set ( i, i - 8, 1. );
      if ( i + 8 < 64 ) ev_matrix.set ( i, i + 8, 1. );
    }

    cerr << endl;
    failed = failed || !testEigenvectorOp<aol::JacobiEigenvectorOp<aol::FullMatrix<double> > > (ev_matrix, "JacobiEigenvectorOp");
    failed = failed || !testEigenvectorOp<aol::QREigenvectorOp<aol::FullMatrix<double> > > (ev_matrix, "QREigenvectorOp");

    failed = failed || !aol::Triangle<double>::isSelfTestOK();


    {
      cerr << "--- Testing matrix symmetry check ... ";

      aol::SparseMatrix<double> symmMat ( 3, 3 ), asymmMat ( 3, 3 );
      symmMat.setIdentity();
      asymmMat.setIdentity(); asymmMat.set ( 0, 1, 1.0 );

      if ( aol::isMatrixSymmetric<aol::SparseMatrix<double>, double>( symmMat ) && ! ( aol::isMatrixSymmetric<aol::SparseMatrix<double>, double>( asymmMat ) ) )
        cerr << "OK." << endl;
      else {
        cerr << "FAILED!" << endl;
        failed =  true;
      }
    }


    {
      cerr << "--- Testing QUOC_ASSERT ... ";

      if ( isQuocAssertOk () )
        cerr << "OK." << endl;
      else {
        cerr << "FAILED!" << endl;
        failed =  true;
      }
    }

    try {
      cerr << "--- Testing overflow detection ... ";

      testSumAndProductWillFit ();

#ifdef DEBUG
      testSqrCubQrt (); // checks are only enabled in debug mode
#endif

      cerr << "OK." << endl;
    }
    catch ( aol::Exception &ex ) {
      ex.consume ();
      cerr << "FAILED!" << endl;
      failed = true;
    }

    {
      cerr << "--- Testing aol::Erf ... ";
      bool erfok = aol::isErfTestOK();
      failed |= ! ( erfok );
      if ( !erfok )
        cerr << "FAILED" << endl;
      else
        cerr << "OK" << endl;
    }

    {
      cerr << "--- Testing aol::appleq, appeq, appgeq ... ";
      bool appok = isappok();

      failed |= ! ( appok );
      if ( !appok )
        cerr << "FAILED" << endl;
      else
        cerr << "OK" << endl;
    }

    {
      cerr << "--- Testing aol::Vector::sortValues ... ";
      bool sortok = aol::isErfTestOK();

      aol::RandomGenerator rg;
      aol::Vector<int> tVec ( 10 );
      for ( int i = 0; i < tVec.size(); ++i ) {
        tVec[i] = rg.rInt( 1024 );
      }
      tVec.sortValues();
      sortok &= ( tVec.crc32OfData() == 447719242 );
      failed |= ! ( sortok );
      if ( !sortok )
        cerr << "FAILED" << endl;
      else
        cerr << "OK" << endl;
    }

    {
      cerr << "--- Testing aol::ProbDistributionFunction in 1D ... ";
      bool pdfOK = true;

      aol::RandomGenerator rg ( 4711 );
      const int szA = 20, szB = 22, szC = 18, szD = 65;
      aol::Vector< double > vecA ( szA ), vecB ( szB ), vecC ( szC ), vecD ( szD );

      for ( int i = 0; i < szA; ++i ) {
        vecA[i] = rg.rReal<double>() + rg.rReal<double>();
      }

      for ( int i = 0; i < szB; ++i ) {
        vecB[i] = rg.rReal<double>() + rg.rReal<double>();
      }

      for ( int i = 0; i < szC; ++i ) {
        vecC[i] = rg.rReal<double>() * rg.rReal<double>();
      }

      aol::ProbDistributionFunction1DForSample<double> pdfA ( vecA ), pdfB ( vecB ), pdfC ( vecC );
      pdfB.computePDFdistTo ( pdfA );
      pdfC.computePDFdistTo ( pdfA );

      aol::PRNGForGivenDistr1D<double> rg2 ( vecB );
      for ( int i = 0; i < szD; ++i ) {
        vecD[i] = rg2.rReal();
      }
      aol::ProbDistributionFunction1DForSample<double> pdfD ( vecD );
      pdfD.computePDFdistTo ( pdfB );

      pdfOK &= ( ( aol::ProbDistFuncHelper<double>::CramerVonMisesProb ( pdfB.getScaledCvMDistanceTo ( pdfA ), szB, szA ) > 0.9 ) &&
                 ( aol::ProbDistFuncHelper<double>::KolmogorovProb (     pdfB.getScaledKSDistanceTo (  pdfA ) ) > 0.9 ) &&
                 ( aol::ProbDistFuncHelper<double>::CramerVonMisesProb ( pdfC.getScaledCvMDistanceTo ( pdfA ), szB, szA ) < 0.01 ) &&
                 ( aol::ProbDistFuncHelper<double>::KolmogorovProb (     pdfC.getScaledKSDistanceTo (  pdfA ) ) < 0.01 ) &&
                 ( aol::ProbDistFuncHelper<double>::KolmogorovProb (     pdfD.getScaledKSDistanceTo (  pdfB ) ) > 0.95 ) );

      if ( !pdfOK )
        cerr << "FAILED" << endl;
      else
        cerr << "OK" << endl;
      failed |= !pdfOK;
    }

    {
      cerr << "--- Testing aol::LookupMap ... ";
      if ( aol::LookupMap< int, int >::isSelfTestOK() ) {
        cerr << "OK" << endl;
      } else {
        cerr << "FAILED" << endl;
        failed = true;
      }
    }

    {
      cerr << "--- Testing aol::Row::isApproxEqual ... ";

      aol::SparseRow<double> r1, r2;
      bool rowok = true;

      if ( ! r1.isApproxEqual ( 1, r2, 0.001 ) ) rowok = false;

      r1.set ( 1, 2, 2 );
      if ( r1.isApproxEqual ( 1, r2, 0.001 ) ) rowok = false;
      if ( r2.isApproxEqual ( 1, r1, 0.001 ) ) rowok = false;

      r2.set ( 1, 3, 3 );
      if ( r1.isApproxEqual ( 1, r2, 0.001 ) ) rowok = false;
      if ( r2.isApproxEqual ( 1, r1, 0.001 ) ) rowok = false;

      r1.set ( 1, 3, 3 );
      r2.set ( 1, 2, 2.0001 );
      if ( ! r1.isApproxEqual ( 1, r2, 0.001 ) ) rowok = false;
      if ( ! r2.isApproxEqual ( 1, r1, 0.001 ) ) rowok = false;

      r2.set ( 1, 7, 1 );
      if ( r1.isApproxEqual ( 1, r2, 0.001 ) ) rowok = false;
      if ( r2.isApproxEqual ( 1, r1, 0.001 ) ) rowok = false;

      if ( rowok ) {
        cerr << "OK" << endl;
      } else {
        cerr << "FAILED" << endl;
      }
      failed |= ( !rowok );
    }

    { // sizes of standard data types
      cout << "Sizes of standard data types used in the library" << endl
           << "bool " << sizeof ( bool ) << endl
           << "char " << sizeof ( char ) << endl
           << "short " << sizeof ( short ) << endl
           << "int " << sizeof ( int ) << endl
           << "float " << sizeof ( float ) << endl
           << "double " << sizeof ( double ) << endl
           << "long double " << sizeof ( long double ) << endl
           << "std::complex<float> " << sizeof ( std::complex<float> ) << endl
           << "std::complex<double> " << sizeof ( std::complex<double> ) << endl;

      failed |= sizeof ( int ) < 4; // e.g. for ScalarArray

      failed |= !TYPESIZETEST(int16_t,2);
      failed |= !TYPESIZETEST(uint16_t,2);
      failed |= !TYPESIZETEST(int32_t,4);
      failed |= !TYPESIZETEST(uint32_t,4);
      failed |= !TYPESIZETEST(int64_t,8);
      failed |= !TYPESIZETEST(uint64_t,8);

      // check whether large enums are treated correctly.
      // Note: 2^31-1=0x7FFFFFFF is the largest number the C99 standard guarantees to work.
      enum /* dummyEnum */ { dummyBigValue = 0x7FFFFFFF, dummyNegValue = -42 };
      failed |= ( dummyBigValue >> 29 != 3 ) || ( dummyNegValue != -42 );
    }

    // test aol::Vector and aol::MultiVector save and load
    failed |= ( ! ( vectorSaveTest<signed char>() ) ||
                ! ( vectorSaveTest<unsigned char>() ) ||
                ! ( vectorSaveTest<signed short>() ) ||
                ! ( vectorSaveTest<unsigned short>() ) ||
                ! ( vectorSaveTest<signed int>() ) ||
                ! ( vectorSaveTest<unsigned int>() ) ||
                ! ( vectorSaveTest<int64_t>() ) ||
                ! ( vectorSaveTest<uint64_t>() ) ||
                ! ( vectorSaveTest<float>() ) ||
                ! ( vectorSaveTest<double>() ) ||
                ! ( vectorSaveTest<long double>() ) );

    {
      aol::BitVector bVec ( 42 ), bVecL;
      for ( int i = 0; i < 42 ; ++i ) {
        bVec.set ( i, ( i % 3 == 0 ) );
      }
      bVec.saveToFile ( "bVec.quo" );
      bVecL.loadFromFile ( "bVec.quo" );
      remove ( "bVec.quo" );
      failed |= ! ( bVec == bVecL );
    }

    {
      failed |= ! (aol::TimeStepRKF45< aol::Vector<double>, DummyReaction >::isSelfTestOK() );
    }

    if( !failed ){
      aol::printSelfTestSuccessMessage ( "--                         AOL Self Test Successful                           --" );
      aol::callSystemPauseIfNecessaryOnPlatform();
      return EXIT_SUCCESS;
    }

    {
      using namespace aol;
      // this should not compile if someone is "using namespace qc;" in the headers included.
      aol::Vector<avoid_namespace_collision> test_vector(1);
    }

  }
  catch(std::exception &ex){
    cerr << color::error << "\n\nstd::exception caught:\n" << color::reset;
    cerr << ex.what () << endl;
  }
  catch(aol::Exception &ex){
    cerr << color::error << "\n\naol::Exception caught:\n" << color::reset;
    ex.dump ();
  }
  catch (...){
    cerr << color::error << "\n\nUnknown exception caught.\n" << color::reset;
  }
  aol::printSelfTestFailureMessage ( "!!                       AOL Self Test FAILED                                 !!" );
  aol::callSystemPauseIfNecessaryOnPlatform();
  return EXIT_FAILURE ;
}
