#ifndef __EIGENWRAPPER_H
#define __EIGENWRAPPER_H

#ifdef USE_EXTERNAL_EIGEN

#include <eigenIncludes.h>
#include <configurators.h>

namespace eig {

/**
 * \brief Wrapper to Eigen that aims to be a substitute for aol::Vector in existing code.
 *
 * So far tested to be a valid VectorType for aol::GridlessGradientDescent.
 *
 * \author Berkels
 */
template <typename _DataType>
class Vector : public Eigen::Matrix<_DataType, Eigen::Dynamic, 1> {
protected:
  typedef Eigen::Matrix<_DataType, Eigen::Dynamic, 1> SuperClass;
public:
  typedef _DataType DataType;

  Vector ( ) { }

  explicit Vector ( int Length ) : SuperClass ( Length ) { }

  explicit Vector ( const qc::GridStructure &Grid ) : SuperClass ( Grid.getNumberOfNodes() ) { }

  Vector ( const Vector<DataType> &Vec, aol::CopyFlag copyFlag = aol::DEEP_COPY ) : SuperClass ( ) {
    switch ( copyFlag ) {
      case aol::DEEP_COPY:
        * ( static_cast<SuperClass*> ( this ) ) = Vec;
        break;

      case aol::STRUCT_COPY:
        this->setZero ( Vec.size(), 1 );
        break;

      case aol::STRUCT_COPY_UNINIT:
        this->resize ( Vec.size(), 1 );
        break;

      default:
        throw aol::Exception ( "eig::Vector<DataType>::Vector( DataType* Data, int n, CopyFlag copyFlag ): illegal copy flag", __FILE__, __LINE__ );
        break;
    };

  }

  Vector ( DataType* Data, int N, aol::CopyFlag CopyFlag ) {
    switch ( CopyFlag ) {
      case aol::DEEP_COPY: {
        * ( static_cast<SuperClass*> ( this ) ) = Eigen::Map<SuperClass> ( Data, N );
      }
        break;

      default:
        throw aol::Exception ( "eig::Vector<DataType>::Vector( DataType* Data, int n, CopyFlag copyFlag ): illegal copy flag", __FILE__, __LINE__ );
    }
  }

  Vector<DataType>& operator= ( const SuperClass &Other ) {
    SuperClass::operator= ( Other );
    return *this;
  }

  void reallocate ( const int Length ) {
    this->resize ( Length, 1 );
  }

  void reallocate ( const Vector<DataType> &other ) {
    this->reallocate ( other.size() );
  }

  DataType dotProduct ( const Vector<DataType> &Vec ) const {
    return this->dot ( Vec );
  }

  DataType operator* ( const Vector<DataType> &Vec ) const {
    return this->dot ( Vec );
  }

  DataType normSqr() const {
    return this->squaredNorm();
  }

  Vector<DataType> &addMultiple ( const Vector<DataType>& Vec, DataType Factor ) {
    (*this) += Factor * Vec;
    return *this;
  }

  DataType getMaxValue() const {
    return this->maxCoeff();
  }

  // Unfortunately, qc::ScalarArray allows const access to a non-const data pointer.
  // For compatibility reasons, we have to allow the same, which needs a const_cast.
  DataType* getData() const {
    return const_cast<_DataType*> ( this->data() );
  }
};

/**
 * \author Berkels
 */
template <typename _DataType, qc::Dimension Dim> class ScalarArray;

template <typename _DataType>
class ScalarArray<_DataType, qc::QC_2D> : public eig::Vector<_DataType > {
  aol::Vec2<int> _size;
public:
  explicit ScalarArray ( const std::string &Filename ) {
    qc::ScalarArray<_DataType, qc::QC_2D> quocArray ( Filename );
    _size[0] = quocArray.getNumX();
    _size[1] = quocArray.getNumY();
    eig::Vector<_DataType >::SuperClass::operator= ( Eigen::Map<typename eig::Vector<_DataType >::SuperClass> ( quocArray.getData(), quocArray.size() ) );
  }

  ScalarArray ( const ScalarArray<_DataType, qc::QC_2D> &Input, aol::CopyFlag CopyFlag = aol::DEEP_COPY )
    : eig::Vector<_DataType > ( Input, CopyFlag ),
      _size ( Input._size ) {
  }

  explicit ScalarArray ( const qc::GridStructure &Grid )
    : eig::Vector<_DataType > ( Grid.getNumX() * Grid.getNumY() ),
      _size ( Grid.getNumX(), Grid.getNumY() ) { }

  eig::Vector<_DataType> &asVector() {
    return static_cast<eig::Vector<_DataType>&> ( *this );
  }

  int getNumX() const {
    return _size[0];
  }

  int getNumY() const {
    return _size[1];
  }

  void save ( const char *FileName, qc::SaveType Type, const char *Comment = NULL ) const {
    qc::ScalarArray<_DataType, qc::QC_2D> quocArray ( _size[0], _size[1], this->getData(), aol::FLAT_COPY );
    quocArray.save ( FileName, Type, Comment );
  }

};

/**
 * \author Berkels
 */
template <typename DataType>
class FullMatrix : public Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> {
  typedef Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> SuperClass;
public:
  FullMatrix ( int Rows, int Columns ) : SuperClass ( Rows, Columns ) { }

  void copyFrom ( const aol::Matrix<DataType> &Mat ) {
    const int numRows = this->rows();
    const int numCols = this->cols();
    // Likely not the most efficient way to convert an aol::Matrix to the Eigen format.
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        (*this) ( i, j ) = Mat.get ( i, j );
  }
};

/**
 * \brief Wrapper to Eigen is a suitable MatrixType for aol::FELinOpInterface::assembleAddMatrix.
 *
 * \author Berkels
 */
template <typename DataType>
class TripletList : public std::vector<Eigen::Triplet<DataType> > {
  typedef std::vector<Eigen::Triplet<DataType> > SuperClass;
public:
  TripletList ( ) : SuperClass ( ) { }

  inline void add ( int i, int j, DataType value ) {
    this->push_back ( Eigen::Triplet<DataType> ( i, j, value ) );
  }

  // This doesn't have a size, so resize does nothing.
  void resize ( const int, const int ) { }
};

/**
 * \brief Wrapper to Eigen that aims to be a suitable MatrixType in our FE Code.
 *        It's not possible to assemble into this directly (would be too inefficient),
 *        but it can filled with a TripletList.
 *
 * \author Berkels
 */
template <typename DataType>
class SparseMatrix : public Eigen::SparseMatrix<DataType> {
  typedef Eigen::SparseMatrix<DataType> SuperClass;
public:
  SparseMatrix ( ) : SuperClass ( ) { }

  SparseMatrix ( int Rows, int Columns ) : SuperClass ( Rows, Columns ) { }

  template <typename GridType>
  explicit SparseMatrix ( const GridType &Grid ) : SuperClass ( Grid.getNumberOfNodes(), Grid.getNumberOfNodes() ) { }

  template <typename MatrixType>
  explicit SparseMatrix ( const aol::SparseBlockMatrix<MatrixType> &BlockMat )
    : SuperClass ( BlockMat.getTotalNumRows(), BlockMat.getTotalNumCols() ) {
    copyFrom ( BlockMat );
  }

  template <typename MatrixType>
  void copyFrom ( const aol::SparseBlockMatrix<MatrixType> &BlockMat ) {
    if ( ( this->rows() != BlockMat.getTotalNumRows() ) || ( this->cols() != BlockMat.getTotalNumCols() ) )
      throw aol::DimensionMismatchException ( "Size mismatch.", __FILE__, __LINE__ );

    // Likely not the most efficient way to convert an aol::SparseBlockMatrix to the Eigen format.
    // Still, compared to runtime of Eigen's direct sparse solvers, this overhead seems to be
    // negligible for now.
    eig::TripletList<DataType> triplets;
    BlockMat.addUnblockedMatrixTo ( triplets );
    // setFromTriplets destroy the current contents of *this, so no need to clear anything beforehand.
    this->setFromTriplets(triplets.begin(), triplets.end());
    // the result of setFromTriplets is already in compressed form.
    //this->makeCompressed();
  }

  inline void add ( int i, int j, DataType value ) {
    throw aol::Exception ( "This method is only implemented for compatibility. Don't assemble into this matrix type directly. Instead, ", __FILE__, __LINE__, __FUNCTION__ );
    // This seems to be an extremely inefficient way to create new values.
    this->coeffRef ( i, j ) += value ;
  }

  void applyAdd ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const {
    Dest += (*this) * Arg;
  }

  void apply ( const Vector<DataType> &Arg, Vector<DataType> &Dest ) const {
    static_cast<Eigen::Matrix<DataType, Eigen::Dynamic, 1>&> (Dest) = (*this) * Arg;
  }
};

/**
 * \brief Variant of qc::RectangularGridConfigurator that uses Eigen data structures for VectorType and MatrixType.
 *
 * \author Berkels
 */
template <typename _RealType, qc::Dimension Dim, typename _QuadType, typename RectangularGridType = qc::RectangularGrid<Dim> >
class RectangularGridConfigurator : public qc::RectangularGridConfigurator<_RealType, Dim, _QuadType, RectangularGridType> {

public:
  RectangularGridConfigurator ( const RectangularGridType &Grid )
    : qc::RectangularGridConfigurator<_RealType, Dim, _QuadType, RectangularGridType> ( Grid ) { }

  typedef Vector<_RealType>        VectorType;
  typedef SparseMatrix<_RealType>  MatrixType;
  typedef ScalarArray<_RealType, Dim> ArrayType;

  MatrixType* createNewMatrix( ) const {
    MatrixType *mat = new MatrixType ( qc::GridSize<Dim> ( this->_grid ) );
    const int maxNumNonZeroesPerRow = aol::Pow ( 3, Dim );
    cerr << maxNumNonZeroesPerRow << endl;
    mat->reserve ( Eigen::VectorXi::Constant ( mat->cols(), maxNumNonZeroesPerRow ) );
    return mat;
  }
};

}

#endif // USE_EXTERNAL_EIGEN

#endif // __EIGENWRAPPER_H
