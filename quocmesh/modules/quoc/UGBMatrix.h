#ifndef __UGBMATRIX_H
#define __UGBMATRIX_H

#include <scalarArray.h>
#include <vec.h>
#include <multiVector.h>
#include <sparseMatrices.h>

namespace qc {

/** @brief template class for efficient multiplication of a band matrix
 *         uses sparse banded block multiplication scheme
 *  @param columns defines how many bands this matrix contains
 *  @param localDimension defines the size of the lookup tables for the advanced add/set routines @see setlc and @see addlc
 *  @param blocksize defines the size of the blocks which are used in the multiplication method. blocksize should be choosen such that sizeof(DataType)*columns*blocksize &lt; Size of the largest cache of the processor.
 *  @param DataType defines the data type to be used for computations
 *  @ingroup Matrix
 */
template < int columns, int localDimension, int blocksize, typename _DataType, typename InitType = qc::GridDefinition > class UGBMatrix : public aol::Op< aol::Vector<_DataType> > {
public:
  typedef  _DataType DataType;
  typedef  DataType* DataTypeP;

protected:
  DataTypeP data[columns];
  int       indices[columns];
  int       offset[localDimension];
  int       localIndices[localDimension][localDimension];
  int       dimension;
  int       _columns;
  int       _localDimension;

public:
  /** Delete me!
   */
  virtual ~UGBMatrix() {
    for ( int i = 0; i < columns; ++i )
      delete [] data[i];
  }

  UGBMatrix() : _columns(columns), _localDimension(localDimension) {}

  /** Create a matrix which is suitable to carry discrete operators on the given grid.
   */
  explicit UGBMatrix ( const InitType &grid ) : _columns(columns), _localDimension(localDimension) {
    init_for_grid ( grid );
    setZero();
  }

  template <Dimension Dim>
  explicit UGBMatrix ( const GridSize<Dim> & gridSize ) : _columns(columns), _localDimension(localDimension) {
    InitType tempGrid ( gridSize );
    init_for_grid ( tempGrid );
    setZero();
  }

  /** Create a matrix of the given dimension. The offset indices and the local offset indices
   *  must be given in Indices resp. LocalIndices
   */
  UGBMatrix ( int Dimension,
              const int Indices[columns],
              const int LocalIndices[localDimension][localDimension] ) : dimension ( Dimension ), _columns(columns), _localDimension(localDimension) {
    int i, j;
    for ( i = 0; i < columns; ++i ) {
      indices[i] = Indices[i];
      data[i] = new DataType[Dimension];
    }

    for ( i = 0; i < localDimension; ++i ) {
      for ( j = 0; j < localDimension; ++j ) {
        localIndices[i][j] = LocalIndices[i][j];
      }
    }
    setZero();
  }

  /** Conversion constructor aol::SparseMatrix -> UGBMatrix.
   *  Make sure sp_mat has the correct sparsity structure.
   */
  UGBMatrix ( const InitType &grid,
              const aol::SparseMatrix< DataType > &sp_mat ) : _columns(columns), _localDimension(localDimension) {
    // OS: untested
    init_for_grid ( grid );

    vector< typename aol::Row< DataType >::RowEntry > sp_row_vec;

    for ( int i = 0; i < sp_mat.getNumRows(); ++i ) {
      sp_mat.makeRowEntries ( sp_row_vec, i );

      for ( unsigned int k = 0; k < sp_row_vec.size(); ++k ) {
        set ( i, sp_row_vec[k].col, sp_row_vec[k].value );
        // let set deal with data not fitting into this UGBMatrix.
      }

    }
    // TODO: write selfTest
  }

  //! makeRowEntries for compatibility with solvers like Gauss-Seidel and Jacobi
  // TODO: test. optimize (vec.resize( appropriate size ))
  void makeRowEntries ( vector< typename aol::Row<DataType>::RowEntry > &vec, const int i ) const {
    vec.clear();
    for ( int j = 0; j < columns; ++j ) {
      const int rowindex = i + indices[j];
      if ( rowindex >= 0 && rowindex < dimension ) {
        typename aol::Row<DataType>::RowEntry entry ( rowindex, data[j][i] );
        vec.push_back ( entry );
      }
    }
  }

  int getNumRows() const {
    return dimension;
  }

  int getNumCols() const {
    return dimension;
  }

  UGBMatrix<columns, localDimension, blocksize, DataType, InitType> &operator *= ( const DataType &scalar ) {
    for ( int i = 0; i < columns; ++i ) {
      for ( int j = 0; j < dimension; ++j ) {
        data[i][j] *= scalar;
      }
    }
    return *this;
  }

  UGBMatrix<columns, localDimension, blocksize, DataType, InitType> &operator += ( const UGBMatrix<columns, localDimension, blocksize, DataType, InitType> &other ) {
    for ( int i = 0; i < columns; ++i ) {
      for ( int j = 0; j < dimension; ++j ) {
        data[i][j] += other.data[i][j];
      }
    }
    return *this;
  }

  UGBMatrix<columns, localDimension, blocksize, DataType, InitType> &operator -= ( const UGBMatrix<columns, localDimension, blocksize, DataType, InitType> &other ) {
    for ( int i = 0; i < columns; ++i ) {
      for ( int j = 0; j < dimension; ++j ) {
        data[i][j] -= other.data[i][j];
      }
    }
    return *this;
  }

  void addMultiple ( const UGBMatrix<columns, localDimension, blocksize, DataType, InitType> &other, const DataType factor = aol::NumberTrait<DataType>::one ) {
    for ( int j = 0; j < columns; ++j ) {
      for ( int i = 0; i < dimension; ++i ) {
        data[j][i] += factor * other.data[j][i];
      }
    }
  }

  void add ( UGBMatrix<columns, localDimension, blocksize, DataType, InitType> &other ) {
    addMultiple ( other /*, 1.0 */ );
  }

  /** Write a PGM file of size of the matrix and which shows the non-zero
   *  entries of the matrix as black dot.
   */
  void writeNonZeroPattern ( const char *fileName, const unsigned char zeroValue = 0, const unsigned char nonzeroValue = 255 ) const {
    qc::ScalarArray<unsigned char, qc::QC_2D> img ( dimension, dimension );

    createNonZeroPattern ( img, zeroValue, nonzeroValue );

    img.save ( fileName, PGM_UNSIGNED_CHAR_BINARY );
  }

  /** Store a non zero pattern in the given image
   */
  void createNonZeroPattern ( qc::ScalarArray<unsigned char, qc::QC_2D> &img, const unsigned char zeroValue = 0, const unsigned char nonzeroValue = 255 ) const {
    img.setAll ( zeroValue );

    for ( int l = 0; l < columns; ++l ) {
      for ( int i = 0; i < dimension; ++i ) {
        const int x = i + indices[l];
        const int y = i;

        if ( data[l][i] && x >= 0 && y >= 0 && x < dimension - 1 && y < dimension - 1 ) {
          int k = x + y * dimension;
          img[k] = nonzeroValue;
        }
      }
    }
  }

  /** Print the contents of this matrix
   */
  void dump ( ostream &out ) {
    ios_base::fmtflags flags = out.flags();
    out.flags ( ios_base::scientific );
    for ( int i = 0; i < columns; ++i ) {
      for ( int j = 0; j < dimension; ++j ) {
        out << static_cast< double > ( data[i][j] ) << " ";
        if ( ( j + 1 ) % 10 == 0 ) out << std::endl;
      }
    }
    out.flags ( flags );
  }

  /** Return one entry of the diagonal
   */
  DataType getDiag ( const int i ) const {
    return data[columns/2][i];
  }


  /** print out one line of the matrix
   */
  void dumpLine ( const int which, ostream &out ) {
    ios_base::fmtflags flags = out.flags();
    out.flags ( ios_base::scientific );
    for ( int i = 0; i < columns; ++i ) {
      out << static_cast< double > ( data[i][which] ) << " ";
    }
    out << std::endl;
    out.flags ( flags );
  }

  /** Strict comparison
   */
  bool operator== ( UGBMatrix<columns, localDimension, blocksize, DataType> &other ) {
    for ( int i = 0; i < columns; ++i ) {
      for ( int j = 0; j < dimension; ++j ) {
        if ( data[i][j] != other.data[i][j] ) {
          std::cerr << "Data differs at pos (" << i << ", " << j << "): "
                    << data[i][j] << " != " << other.data[i][j] << std::endl;
          return false;
        }
      }
    }
    return true;
  }

  /** assignment operator
   */
  UGBMatrix<columns, localDimension, blocksize, DataType, InitType> &operator= ( const UGBMatrix<columns, localDimension, blocksize, DataType, InitType> &other ) {

    // Beware of self-assignment
    if ( this == &other ) return *this;

    for ( int i = 0; i < columns; ++i ) {
      memcpy ( data[i], other.data[i], dimension*sizeof ( DataType ) );
    }
    return *this;
  }


  /** Approximate comparison
   */
  bool isApproxEqual ( const UGBMatrix<columns, localDimension, blocksize, DataType, InitType> &other,
                     DataType epsilon ) {
    for ( int i = 0; i < columns; ++i ) {
      for ( int j = 0; j < dimension; ++j ) {
        if ( fabs ( data[i][j] - other.data[i][j] ) > epsilon ) {
          std::cerr << "Data differs by " << fabs ( data[i][j] - other.data[i][j] ) << " at pos ("
                    << i << ", " << j << "): "
                    << data[i][j] << " != " << other.data[i][j] << std::endl;
          return false;
        }
      }
    }
    return true;
  }

  /** Copy constructor
   */
  UGBMatrix ( const UGBMatrix<columns, localDimension, blocksize, DataType, InitType> &m ) {
    int i, j;
    dimension = m.dimension;

    for ( i = 0; i < columns; ++i ) {
      indices[i] = m.indices[i];
      data[i] = new DataType[dimension];
      memcpy ( data[i], m.data[i], dimension*sizeof ( DataType ) );
    }
    for ( i = 0; i < localDimension; ++i ) {
      for ( j = 0; j < localDimension; ++j ) {
        localIndices[i][j] = m.localIndices[i][j];
      }
    }

    std::cerr << "\n\nCOPYING MATRIX!!!\n\n" << std::endl;
  }

  /** Add the contents of the localMatrix for the element at the given
   *  global indices to the matrix
   */
  void add ( const int elIndex, const DataType localMatrix[localDimension][localDimension] ) {

    for ( int i = 0; i < localDimension; ++i ) {
      const int row = elIndex + offset[i];
      for ( int j = 0; j < localDimension; ++j ) {
        const int colIndex = localIndices[i][j];
        data[colIndex][row] += localMatrix[i][j];
        if ( colIndex == 0 && row == 0 ) std::cerr << localMatrix[i][j] << " ";
      }
    }
  }
  // may want to add bounds_check here.

  /** set all matrix entries to zero
   */
  void setZero() {
    for ( int i = 0; i < columns; ++i ) {
      memset ( data[i], 0, sizeof ( DataType ) *dimension );
    }
  }

  /** add value to entry at (row, col)
   *  @see addlc
   */
  void add ( int row, int col, DataType value ) {
#ifdef BOUNDS_CHECK
    if ( ( row < 0 ) || ( row > dimension ) ) {
      cerr << "trying to add to " << row << ", " << col << ", this entry does not exist!" << endl;
      throw aol::Exception ( "UGBMatrix::set: incorrect row index!", __FILE__, __LINE__ );
    }
    // entry_found is not needed, since the function is exited as soon as a valid position has been found!
    //bool entry_found = false; // needed below
#endif
    const int rowIndex = col - row;

    for ( int l = 0; l < columns; ++l ) {
      if ( indices[l] == rowIndex ) {
        data[l][row] += value;
//#ifdef BOUNDS_CHECK
//  entry_found = true;
// #endif
        return;
      }
    }
    // No good index has been found

#ifdef BOUNDS_CHECK
    char message[1024];
    sprintf ( message, "UGBMatrix: Cannot add at row %d, col %d -> bandIndex %d", row, col, rowIndex );
    throw aol::Exception ( message, __FILE__, __LINE__ );
#endif
  }

  /** add value to local indices (li, lj) in (row, *)
   *  @see add
   */
  void addlc ( int row, int li, int lj, DataType value ) {
    const int l = localIndices[li][lj];
    data[l][row] += value;
  }

  void save ( const char *fileName ) {
    FILE *out;

    if ( ! ( out = fopen ( fileName, "wb" ) ) ) {
      throw aol::Exception ( "Can not open file for writing", __FILE__, __LINE__ );
    }
    for ( int l = 0; l < columns; ++l ) {
      fwrite ( data[l], sizeof ( DataType ), dimension, out );
    }
    fclose ( out );
  }

  void load ( const char *fileName ) {
    FILE *in;

    if ( ! ( in = fopen ( fileName, "rb" ) ) ) {
      throw aol::Exception ( "Can not open file for reading", __FILE__, __LINE__ );
    }
    for ( int l = 0; l < columns; ++l ) {
      fread ( data[l], sizeof ( DataType ), dimension, in );
    }
    fclose ( in );
  }

  /** set value to entry at (row, col)
   *  @see setlc
   */
  void set ( int row, int col, DataType value )  {
    // may want to treat value == 0.0 case differently.
#ifdef BOUNDS_CHECK
    if ( ( row < 0 ) || ( row > dimension ) ) {
      cerr << "trying to set " << row << ", " << col << ", this entry does not exist!" << endl;
      throw aol::Exception ( "UGBMatrix::set: incorrect row index!", __FILE__, __LINE__ );
    }
    // entry_found is not needed, since the function is exited as soon as a valid position has been found!
//     bool entry_found = false; // needed below
#endif
    const int rowIndex = col - row;

    for ( int l = 0; l < columns; ++l ) {
      if ( indices[l] == rowIndex ) {
        data[l][row] = value;
// #ifdef BOUNDS_CHECK
//  entry_found = true;
// #endif
        return;
      }
    }

#ifdef BOUNDS_CHECK
    // No good index has been found
    char message[1024];
    sprintf ( message, "UGBMatrix: Cannot add at row %d, col %d -> bandIndex %d", row, col, rowIndex );
    throw aol::Exception ( message, __FILE__, __LINE__ );
#endif
  }

  /** set value at local indices (li, lj) in (row, *)
   *  @see set
   */
  void setlc ( int row, int li, int lj, DataType value ) {
    const int l = localIndices[li][lj];
    data[l][row] = value;
  }

  /** get value at (row, col)
   */
  DataType get ( int row, int col ) {
    // OS: untested
#ifdef BOUNDS_CHECK
    if ( ( row < 0 ) || ( row > dimension ) || ( col < 0 ) || ( col > dimension ) ) {
      cerr << "trying to get " << row << ", " << col << ", this entry does not exist!" << endl;
      throw aol::Exception ( "UGBMatrix::set: incorrect row or column index!", __FILE__, __LINE__ );
    }
#endif
    const int rowIndex = col - row;

    for ( int l = 0; l < columns; ++l ) {
      if ( indices[l] == rowIndex ) {
        return ( data[l][row] );
      }
    }

    return ( static_cast<DataType> ( 0.0 ) ); // if no entry found above
  }

  void info() {
    DataType mi = 1e10;
    DataType ma = -1e10;
    for ( int l = 0; l < columns; ++l ) {
      for ( int r = 0; r < dimension; ++r ) {
        if ( data[l][r] > ma ) ma = data[l][r];
        if ( data[l][r] < mi ) mi = data[l][r];
      }
    }

    cerr << "min/max of UGBMatrix " << mi << "/" << ma << endl;
  }


  /** delete entries in row and column rowCol
   */
  void setRowColToZero ( int rowCol ) {
    for ( int l = 0; l < columns; ++l ) {
      data[l][rowCol] = 0.0;
      int possibleRow = rowCol - indices[l];
      if ( possibleRow >= 0 && possibleRow < dimension ) {
        data[l][possibleRow] = 0.0;
      }
    }
  }

  void setRowToZero ( const int row ) {
    for ( int l = 0; l < columns; ++l ) {
      data[l][row] = 0.0;
    }
  }

  int getDimension() {
    return ( dimension );
  }

  /** perform matrix vector multiplication on multi vectors
   *
   *  @note The method uses the applyAdd method for real type vectors
   *  @note The method assumes there is no overlap between Dest and Arg
   */
  virtual void applyAdd ( const aol::MultiVector<DataType> &Arg, aol::MultiVector<DataType> &Dest ) const {
    const aol::Vector<DataType> &arg = Arg[0];
    aol::Vector<DataType> &dest = Dest[0];
    applyAdd ( arg, dest );
  }

  /** perform matrix vector multiplication
   *  the method used here is a sparse banded block multiplication scheme
   *  which takes advantage of data localization in the cache of the processor
   *  With use of the intel compiler icc/ecc and the pragma directives
   *  #pragma ivdep and efficient vectorization of the code can be achieved by
   *  setting the define VECTORIZE_INTEL during compilation.
   *
   *  @note The method assumes there is no overlap between Dest and Arg
   */
  virtual void applyAdd ( const aol::Vector<DataType> &arg, aol::Vector<DataType> &dest ) const {
    int i, l;
    // Why the hell was dest being cleared here?
    //     dest.clear();

    const DataTypeP destPtr = dest.getData();
    const DataTypeP argPtr = arg.getData();
    DataTypeP dptr[columns], aptr[columns], mptr[columns];
    int startIndex[columns], mi[columns], endIndex[columns];

    // Initialize pointers
    for ( l = 0; l < columns; ++l ) {
      startIndex[l] = aol::Max ( -indices[l], 0 );
      mi[l] = aol::Max ( indices[l], 0 );
      endIndex[l] = dimension - startIndex[l] - mi[l] ;

      mptr[l] = data[l] + startIndex[l];
      dptr[l] = destPtr + startIndex[l];
      aptr[l] = argPtr + startIndex[l] + indices[l];
    }

    const int nb = dimension / ( blocksize );
    const int nr = dimension % ( blocksize );

    if ( arg.size() != dest.size() || arg.size() != dimension ) {
      throw aol::Exception ( "Dimensions do not match", __FILE__, __LINE__ );
    }

    // Loop over nb blocks of size blocksize
    for ( int b = 0; b < nb; ++b ) {
      for ( l = 0; l < columns; ++l ) {
        DataTypeP dp = dptr[l];
        DataTypeP ap = aptr[l];
        DataTypeP mp = mptr[l];

        const int ei = aol::Min ( endIndex[l], blocksize );

#ifdef VECTORIZE_INTEL
#pragma ivdep
        for ( i = 0; i < ei; ++i ) dp[i] += ap[i] * mp[i];
        dptr[l] += ei;
        aptr[l] += ei;
        mptr[l] += ei;
#else
        for ( i = ei; i; i-- ) *dp++ += *ap++ * *mp++;

        dptr[l] = dp;
        aptr[l] = ap;
        mptr[l] = mp;
#endif

        endIndex[l] -= ei;
      }
    }
    if ( nr ) {
      for ( l = 0; l < columns; ++l ) {
        DataTypeP dp = dptr[l];
        DataTypeP ap = aptr[l];
        DataTypeP mp = mptr[l];
        const int ei = endIndex[l];
#ifdef VECTORIZE_INTEL
#pragma ivdep
        for ( i = 0; i < ei; ++i ) dp[i] += ap[i] * mp[i];
#else
        for ( i = ei; i; i-- ) *dp++ += *ap++ * *mp++;
#endif
      }
    }
  }

  //! M_result = M1 * this
  void multiplyFromLeft ( const aol::GenSparseMatrix< DataType > &M1, aol::GenSparseMatrix< DataType > &M_result ) {
    // OS: untested
    if ( M1.getNumRows() != M_result.getNumRows()  ||  M1.getNumCols() != this->getDimension()  ||  this->getDimension() != M_result.getNumCols() ) {
      throw aol::Exception ( "UGBMatrix::multiplyFromLeft: incompatible matrix sizes.\n", __FILE__, __LINE__ );
    }
    M_result.setZero();

    vector<typename aol::Row<DataType>::RowEntry > vec;
    for ( int i = 0; i < M1.getNumRows(); ++i ) {
      M1.getRow ( i ).makeRowEntries ( vec, i );
      for ( int j = 0; j < this->getDimension(); ++j ) {
        for ( typename vector< typename aol::Row<DataType>::RowEntry >::iterator it1 = vec.begin(); it1 != vec.end(); ++it1 ) {
          M_result.add ( i, j, it1->value * this->get ( it1->col, j ) );
        }
      }
    }
  }

private:
  void init_for_grid ( const InitType &grid ) {
    int i, j;
    const int numX = grid.getWidth();
    const int numY = grid.getHeight();

    dimension = grid.getNumberOfNodes();

    if ( grid.getDimOfWorld() == qc::QC_3D ) {  // Grid is 3D
      if ( localDimension != 8 || columns != 27 ) {
        throw aol::Exception ( "localDimension or columns do not match given grid", __FILE__, __LINE__ );
      }
      const int nn = numX * numY;
      indices[9] =  -numX - 1;
      indices[10] =  -numX;
      indices[11] =  -numX + 1;
      indices[12] =  -1;
      indices[13] =  0;
      indices[14] =  1;
      indices[15] =  numX - 1;
      indices[16] =  numX;
      indices[17] =  numX + 1;

      for ( i = 0; i < 9; ++i ) indices[i] = indices[i+9] - nn;
      for ( i = 18; i < 27; ++i ) indices[i] = indices[i-9] + nn;

      offset[0] = 0;
      offset[1] = 1;
      offset[2] = numX;
      offset[3] = numX + 1;

      for ( i = 4; i < 8; ++i ) offset[i] = offset[i-4] + nn;
    } else if ( grid.getDimOfWorld() == qc::QC_2D ) { // Grid is 2D
      if ( localDimension != 4 || columns != 9 ) {
        throw aol::Exception ( "localDimension or columns do not match given grid", __FILE__, __LINE__ );
      }
      indices[0] =  -numX - 1;
      indices[1] =  -numX;
      indices[2] =  -numX + 1;
      indices[3] =  -1;
      indices[4] =  0;
      indices[5] =  1;
      indices[6] =  numX - 1;
      indices[7] =  numX;
      indices[8] =  numX + 1;

      offset[0] = 0;
      offset[1] = 1;
      offset[2] = numX;
      offset[3] = numX + 1;
    } else {
      throw aol::Exception ( "Illegal dimension", __FILE__, __LINE__ );
    }

    // Create data vectors
    for ( i = 0; i < columns; ++i ) {
      data[i] = new DataType[dimension];
    }

    for ( i = 0; i < localDimension; ++i ) {
      const int i1 = offset[i];
      for ( j = 0; j < localDimension; ++j ) {
        const int j1 = offset[j];
        const int rowIndex = j1 - i1;
        for ( int l = 0; l < columns; ++l ) {
          if ( indices[l] == rowIndex ) {
            localIndices[i][j] = l;
            break;
          }
        }
      }
    }

  }

};


/**
 *  \brief  UGB Matrix for graphs, i.e. grids with reduced neighbor relations (left, right, up, down, front, back).
 *
 *  This class generates a UGBMatrix with the sparsity structure needed for graphs.
 *
 *  \author Torben Paetz (JUB)
 *  \ingroup Matrix
 */
template <int blocksize, typename _DataType, typename InitType>
class UGBMatrixForGraph : public qc::UGBMatrix<7, 4, blocksize, _DataType, InitType > {

public:

  //! The data type used for the matrix elements.
  typedef  _DataType DataType;

  //! Creates a matrix which is suitable for calculations on the given graph.
  //! \param graph The graph for which the matrix has to be build.
  explicit UGBMatrixForGraph ( const InitType &graph ) {
    init_for_graph ( graph.getSize() );
    this->setZero();
  }

  //! Creates a matrix which is suitable for calculations on a graph related to the given size.
  //! \param graph The graph for which the matrix has to be build.
  explicit UGBMatrixForGraph ( const qc::GridSize<qc::QC_3D> &graphSize ) {
    init_for_graph ( graphSize );
    this->setZero();
  }

private:

  //! Sets up the sparsity structure needed for graphs.
  //! \param size The size of the graph.
  void init_for_graph ( const qc::GridSize<qc::QC_3D> &size ) {
    int i, j;
    const int numX = size[0];
    const int numY = size[1];

    this->dimension = size[0] * size[1] * size[2];

    const int nn = numX * numY;
    this->indices[0] = -nn;
    this->indices[1] =  -numX;
    this->indices[2] =  -1;
    this->indices[3] =  0;
    this->indices[4] =  1;
    this->indices[5] =  numX;
    this->indices[6] = numX*numY;

    this->offset[0] = 0;
    this->offset[1] = 1;
    this->offset[2] = numX;
    this->offset[3] = nn;

    // Create data vectors
    for ( i = 0; i < this->_columns; ++i ) {
      this->data[i] = new DataType[this->dimension];
    }

    for ( i = 0; i < this->_localDimension; ++i ) {
      const int i1 = this->offset[i];
      for ( j = 0; j < this->_localDimension; ++j ) {
        const int j1 = this->offset[j];
        const int rowIndex = j1 - i1;
        for ( int l = 0; l < this->_columns; ++l ) {
          if ( this->indices[l] == rowIndex ) {
            this->localIndices[i][j] = l;
            break;
          }
        }
      }
    }
  }


};

} // end namespace

#endif
