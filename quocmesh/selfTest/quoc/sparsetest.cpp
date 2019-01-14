#include <scalarArray.h>
#include <bandMatrix.h>
#include <configurators.h>
#include <diagBandMatrix.h>
#include <FEOpInterface.h>
#include <mcm.h>
#include <quoc.h>
#include <solver.h>
#include <quocMatrices.h>
#include <verySparseMatrix.h>


template< typename MatrixType, typename GridSizeType, typename OpType >
bool doCheckMatrixOpReOp ( const GridSizeType &gridSize, const OpType &op, MatrixType &mat ) {
  mat.setZero();
  op.assembleAddMatrix ( mat );
  const bool compOpsOK = compareOps( mat, op, gridSize.getNumberOfNodes(), gridSize.getNumberOfNodes(), 1.0e-18 );
  if ( compOpsOK )
    cerr << "OK ";
  else
    cerr << "FAILED ";

  aol::RowEntryOp< double, MatrixType > REOp ( mat, gridSize.getNumberOfNodes() );
  const bool compReOpOK = compareOps( REOp, op, gridSize.getNumberOfNodes(), gridSize.getNumberOfNodes(), 1.0e-18 );
  if ( compReOpOK )
    cerr << "OK ";
  else
    cerr << "FAILED ";

  return ( compOpsOK && compReOpOK );
}


template <typename MatType>
bool checkDelRowCol( MatType &mat ) {

  bool okay = true;

  const int indexValues[3] = {0, 9, mat.getNumRows() / 2 };
  for ( int k = 0; k < 3; ++k ) {
    const int index = indexValues[k];
    mat.setRowColToZero( index );

    for( int i = 0; i < mat.getNumRows(); ++i ){
      if( mat.get(i, index) != 0.0 || mat.get(index, i) != 0.0 ){
        okay = false;
        cerr << "Error: " << i << " " << index << " " << mat.get(i, index) << "   " << index << " " << i << " " << mat.get(index, i) << endl;
      }
    }
  }

  return( okay );
}


template< typename MatrixType, typename GridSizeType, typename OpType >
bool checkMatrix ( const char* name, const GridSizeType &gridSize, const OpType &op ) {
  cerr << "checking " << name;
  MatrixType mat ( gridSize );
  const bool doCheckOK = doCheckMatrixOpReOp ( gridSize, op, mat );

  const bool drcOK = checkDelRowCol( mat );
  if ( drcOK )
    cerr << "OK" << endl;
  else
    cerr << "FAILED" << endl;

  return ( doCheckOK && drcOK );
}


template< typename MatrixType, typename GridSizeType, typename OpType >
bool checkMatrixR ( const char* name, const GridSizeType &gridSize, const OpType &op ) {
  cerr << "checking " << name;
  MatrixType mat ( gridSize );
  const bool ret = doCheckMatrixOpReOp ( gridSize, op, mat );
  cerr << endl;
  return ( ret );
}


template< typename MatrixType, typename GridSizeType, typename OpType >
bool checkMatrixS ( const char* name, const GridSizeType &gridSize, const OpType &op ) {
  cerr << "checking " << name;
  MatrixType mat ( gridSize.getNumberOfNodes(), gridSize.getNumberOfNodes() );
  const bool ret = doCheckMatrixOpReOp ( gridSize, op, mat );
  cerr << endl;
  return ( ret );
}


int main( ) {

  bool OK = true;

  try {
    const int MY_CACHE_SIZE = 1024 * 1024; // 1024 KB, to be divided by number of bands (rounded up to nearest power of 2) times sizeof(double), used for UGBMatrix

    // ---------------- 2d ----------------------------
    // here we just test a MassOp with identical entries at positions "global index of position" ( x,y ) and "global index of position" ( y, x )
    qc::CubicGrid<qc::QC_2D> grid( 3 );
    qc::GridSize<qc::QC_2D> gridSize = qc::GridSize<qc::QC_2D> ( grid );
    typedef qc::RectangularGridConfigurator< double, qc::QC_2D, aol::GaussQuadrature< double, qc::QC_2D, 3 > > ConfigType;
    const aol::MassOp<ConfigType>   M_2d ( grid, aol::ASSEMBLED );

    //           matrix type used for assembling and comparing operator to            grid
    OK &= checkMatrix < aol::SparseMatrix<double>,                                                                  qc::GridSize<qc::QC_2D>, aol::MassOp<ConfigType> > ( "aol::SparseMatrix...................................", gridSize, M_2d );
    OK &= checkMatrixS< aol::VerySparseMatrix<double>,                                                              qc::GridSize<qc::QC_2D>, aol::MassOp<ConfigType> > ( "aol::VerySparseMatrix...............................", gridSize, M_2d );
    OK &= checkMatrix < qc::UniformGridSparseMatrix<double>,                                                        qc::GridSize<qc::QC_2D>, aol::MassOp<ConfigType> > ( "qc::UniformGridSparseMatrix.........................", gridSize, M_2d );
    OK &= checkMatrix < qc::MultilinFEBandMatrix<double, qc::QC_2D>,                                                qc::GridSize<qc::QC_2D>, aol::MassOp<ConfigType> > ( "qc::MultilinFEBandMatrix............................", gridSize, M_2d );
    OK &= checkMatrixS< aol::DiagBandMatrix<double, 10, 10>,                                                        qc::GridSize<qc::QC_2D>, aol::MassOp<ConfigType> > ( "qc::DiagBandMatrix..................................", gridSize, M_2d );
    OK &= checkMatrix < qc::FastUniformGridMatrix<double, qc::QC_2D>,                                               qc::GridSize<qc::QC_2D>, aol::MassOp<ConfigType> > ( "qc::FastUniformGridMatrix...........................", gridSize, M_2d );
    OK &= checkMatrixR< qc::FastAssembleUniformGridMatrix<double>,                                                  qc::GridSize<qc::QC_2D>, aol::MassOp<ConfigType> > ( "qc::FastAssembleUniformGridMatrix...................", gridSize, M_2d );
    OK &= checkMatrix < qc::UGBMatrix < 9, 4, MY_CACHE_SIZE / ( 16 * 8 ), double, qc::RectangularGrid<qc::QC_2D> >, qc::GridSize<qc::QC_2D>, aol::MassOp<ConfigType> > ( "qc::UGBMatrix.......................................", gridSize, M_2d );
    OK &= checkMatrixR< qc::UniGridCSR_Matrix<qc::QC_2D>,                                                           qc::GridSize<qc::QC_2D>, aol::MassOp<ConfigType> > ( "qc::UniGridCSR_Matrix...............................", gridSize, M_2d );

    // ---------------- 3d ----------------------------
    // the MCMMassOp does not have identical entries at positions "global index of position" ( x, y, z ) and "global index of position" ( permutation  ( x, y, z ) )
    qc::CubicGrid<qc::QC_3D> grid3d( 3 );
    qc::GridSize<qc::QC_3D> gridSize3d = qc::GridSize<qc::QC_3D> ( grid3d );
    typedef qc::MCMMassOp< qc::RectangularGridConfigurator< double, qc::QC_3D, aol::GaussQuadrature< double, qc::QC_3D, 3 > > > MCMMassOpType;
    MCMMassOpType MCMMassOp3d( grid3d, aol::ONTHEFLY, 1.0 );

    qc::ScalarArray<double, qc::QC_3D> img3d( grid3d );
    img3d.load( "../../examples/testdata/volume_9.dat.bz2" );
    img3d /= img3d.getMaxValue();

    MCMMassOp3d.setImageReference( img3d );

    OK &= checkMatrix < aol::SparseMatrix<double>,                                                       qc::GridSize<qc::QC_3D>, MCMMassOpType > ( "aol::SparseMatrix (3D)..............................", gridSize3d, MCMMassOp3d );
    OK &= checkMatrix < qc::UniformGridSparseMatrix<double>,                                             qc::GridSize<qc::QC_3D>, MCMMassOpType > ( "qc::UniformGridSparseMatrix (3D)....................", gridSize3d, MCMMassOp3d );
    OK &= checkMatrix < qc::MultilinFEBandMatrix<double,qc::QC_3D>,                                      qc::GridSize<qc::QC_3D>, MCMMassOpType > ( "qc::MultilinFEBandMatrix............................", gridSize3d, MCMMassOp3d );
    OK &= checkMatrix < qc::FastUniformGridMatrix<double,qc::QC_3D>,                                     qc::GridSize<qc::QC_3D>, MCMMassOpType > ( "qc::FastUniformGridMatrix (3D)......................", gridSize3d, MCMMassOp3d );
    // there is no FastAssembleUniformGridMatrix in 3D
    OK &= checkMatrix < qc::UGBMatrix<27,8,MY_CACHE_SIZE/(32*8),double,qc::RectangularGrid<qc::QC_3D> >, qc::GridSize<qc::QC_3D>, MCMMassOpType > ( "qc::UGBMatrix (3D)..................................", gridSize3d, MCMMassOp3d );
    OK &= checkMatrixR< qc::UniGridCSR_Matrix<qc::QC_3D>,                                                qc::GridSize<qc::QC_3D>, MCMMassOpType > ( "qc::UniGridCSRMatrix (3D)...........................", gridSize3d, MCMMassOp3d );

    // ---------------- symmetry check ----------------
    cerr << "Symmetry check on SparseMatrix: ";
    aol::SparseMatrix<double> mat ( grid );
    M_2d.assembleAddMatrix ( mat );
    const bool symOK = mat.isSymmetric ( 1E-14 );
    cerr << ( symOK ? "OK" : "FAILED" ) << endl;
    OK &= symOK;

    if ( OK ) {
      aol::printSelfTestSuccessMessage ( "--                      QUOC SparseTest Successful                            --" );
      return ( EXIT_SUCCESS );
    }
  } catch ( aol::Exception &ex ) {
    ex.dump();
  }
  aol::printSelfTestFailureMessage ( "!!                      QUOC SparseTest FAILED                                !!" );
  return ( EXIT_FAILURE );
}
