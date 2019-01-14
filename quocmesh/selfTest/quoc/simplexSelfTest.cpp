#include <string>

using namespace std;

#include <simplexGrid.h>
#include <simplexConfigurators.h>
#include <FEOpInterface.h>
#include <quocMatrices.h>
#include <maskedOp.h>
#include <solver.h>
#include <preconditioner.h>

using namespace aol;
using namespace qc;

const   qc::Dimension Dim =                                  QC_2D;
typedef double                                               RealType;
typedef qc::simplex::MidpointQuadrature<RealType, Dim>       QuadType;
typedef qc::GridDefinition                                   CubicGridType;
typedef qc::simplex::GridStructure<CubicGridType, Dim>       GridType;
typedef qc::simplex::ConfiguratorTraitLinear<RealType,
                                    Dim, QuadType, GridType> ConfType;
typedef ConfType::MaskType                                   MaskType;
typedef ConfType::ArrayType                                  ArrayType;
typedef qc::MultilinFEBandMatrix<RealType,qc::QC_2D>         MatrixType; // should this actually be AffineFEBandMatrix, which does not exist in 2D yet?

//---------------------------------------------------------------------------

void solvePoissonProblem ( const GridType & Grid,
                           const ArrayType & BdryValues,
                           const MaskType &  BdryMask,
                           ArrayType &       Sol ) {

  ArrayType rhs ( BdryValues ),
            bdryValues ( BdryValues, DEEP_COPY );

  // assert that bdryValues is only != 0 on bdry nodes
  bdryValues.setAllMasked ( 0., BdryMask );

  // right hand side (let apply write to all nodes, values on
  //                  bdry nodes will be reset separately)
  /*
  typedef qc::RegularSimplexConfigurator<RealType, qc::QC_2D,
                   aol::TriangleIntegration <RealType, qc::QC_2D, 1>,
                   qc::GridDefinition> Conf2dType;
  */
  aol::StiffOp<ConfType> stiffOp ( Grid );
  aol::MaskedOp<MatrixType> stiffMatrix ( GridSize<ConfType::DomDim> ( Grid ), BdryMask, INCLUDE_ALL_WRITE_ALL );
  stiffOp.assembleAddMatrix ( stiffMatrix );
  stiffMatrix.apply( bdryValues, rhs );
  // multiply inner nodes with -1.
  rhs *= -1.;

  // set boundary nodes to zero: This is necessary if we don't use
  // maskedVectors, because the maskedDTMatrix is NOT the identity
  // on masked regions, but skips them. Thus the residual in the
  // cg-method would also include boundary-values.
  rhs.setAllMasked ( 0., BdryMask, true );

  // 3. ******* Solve *******
  stiffMatrix.setIncludeWriteMode(INCLUDE_INT_WRITE_INT);
  // preconditioner
  aol::DiagonalPreconditioner< aol::Vector<RealType> > precond( stiffMatrix );

  aol::PCGInverse< Vector<RealType> > cg ( stiffMatrix, // matrix to invert
                                           precond,     // preconditioner
                                           1e-25,       // stopping residual squared
                                           10000,       // max step number
                                           aol::STOPPING_ABSOLUTE,
                                           clog         // output stream
                                          );
  cg.setQuietMode();
  cg.apply( rhs, Sol );

  // set boundary nodes to bdryValues:
  // third argument is "true" to invert the mask such that
  // bdry node become "true", inner "false", hence only
  // bdry nodes are assigned.
  Sol.assignMaskedFrom(bdryValues, BdryMask, true);
}

//---------------------------------------------------------------------------

void setBoundaryValuesFundSol ( const GridType & grid,
                                ArrayType & bdryValues ) {
  int n = bdryValues.getNumX();
  RealType h = 1. / (n - 1);

  GridType::OldFullNodeIterator iter;

  for (iter = grid.begin(); iter != grid.end(); ++iter) {
    RealType x = (*iter)[0] * h + 1.,
             y = (*iter)[1] * h + 1.,
             r = sqrt ( x*x + y*y );
    bdryValues.set ( *iter, log(r) );
  }
}

//---------------------------------------------------------------------------

void computeExactSol ( const GridType & grid,
                       ArrayType & exactSol ) {

  // if the fundamental solution u(x)=log |x| is given
  // as boundary values, we can give this function
  // as exact solution:
  int n = exactSol.getNumX();
  RealType h = 1. / (n - 1);

  GridType::OldFullNodeIterator iter;

  for (iter = grid.begin(); iter != grid.end(); ++iter) {
    RealType x = (*iter)[0] * h + 1.,
             y = (*iter)[1] * h + 1.,
             r = sqrt ( x*x + y*y );
    exactSol.set ( *iter, log(r) );
  }
}

//---------------------------------------------------------------------------

void setMaskXYBoundary ( MaskType & mask ) {
  int n = mask.getNumX();
  mask.setAll( true );
  for (int i = 0; i < n; ++i) {
      mask.set (   i,   0, false );
      mask.set (   i, n-1, false );
      mask.set (   0,   i, false );
      mask.set ( n-1,   i, false );
    }
}

//---------------------------------------------------------------------------

int main() {

  int numSimpl = qc::simplex::TopologyLookup<Dim>::numSimplexesPerCube;

  try {
    bool success = true;

    CubicGridType cubicGrid ( 4, Dim );
    GridType simplexGrid ( cubicGrid );
    ConfType configurator ( simplexGrid );

    // ----------------------------------------------------------------------
    cerr << "Testing simplex element iterator... ";

    vector<qc::BitArray<qc::QC_2D>*> simplexReached;
    for (int i = 0; i < numSimpl; ++i) {
      simplexReached.push_back ( new qc::BitArray<qc::QC_2D> ( GridSize<QC_2D>::createFrom ( cubicGrid ) ) );
      simplexReached[i]->setAll ( false );
    }

    for (GridType::OldFullElementIterator iter = simplexGrid.begin(); iter != simplexGrid.end(); ++iter)
      if (simplexReached[iter->getSimplexNumber()]->get(iter->getCubicElement())) {
        cerr << "Element " << iter->getCubicElement() << ", BFS " << iter->getSimplexNumber() << " twice iterated." << endl;
        success = false;
      }
      else
        simplexReached[iter->getSimplexNumber()]->set ( iter->getCubicElement(), true );

    for (int i = 0; i < numSimpl; ++i)
      for (CubicGridType::OldFullElementIterator iter = cubicGrid.begin(); iter != cubicGrid.end(); ++iter)
        if (!simplexReached[i]->get(*iter)) {
          cerr << "Element " << *iter << ", BFS " << i << " not iterated." << endl;
          success = false;
        }

    for (int i = 0; i < numSimpl; ++i)
      delete simplexReached[i];

    if ( success )
      cerr << "OK" << endl;
    else
      cerr << "FAILED" << endl;

    // ----------------------------------------------------------------------

    cerr << "Assembling stiffness matrix... ";

    MatrixType stiffMat ( GridSize<ConfType::DomDim>::createFrom ( cubicGrid ) );
    aol::StiffOp<ConfType, ConfType::IndexMode> stiffOp ( simplexGrid );
    stiffOp.assembleAddMatrix ( stiffMat );

    cerr << "OK" << endl;

    // ----------------------------------------------------------------------

    cerr << "Assembling mass matrix... ";

    MatrixType massMat ( GridSize<ConfType::DomDim>::createFrom ( cubicGrid ) );
    aol::MassOp<ConfType, ConfType::IndexMode> massOp ( simplexGrid );
    massOp.assembleAddMatrix ( massMat );
    ConfType::ArrayType ones ( cubicGrid );
    ones.setAll ( aol::ZOTrait<RealType>::one );
    ConfType::ArrayType M_times_one ( cubicGrid );
    massMat.apply ( ones, M_times_one );
    cerr << "M 1 * 1 = " << M_times_one * ones << "\t";

    cerr << "OK" << endl;

    // ----------------------------------------------------------------------


    cerr << "Solving Laplace u = 0 with fundamental solution "
            "as boundary values ... " << endl << endl
         << "level\tL^2 error" << endl;

    for (int level = 2; level < 7; ++level) {
      // solve 2d Poisson problem:
      CubicGridType cubicGrid2 ( level, qc::QC_2D );
      GridType grid ( cubicGrid2 );
      MaskType boundaryMask ( GridSize<QC_2D>::createFrom ( cubicGrid2 ) );
      setMaskXYBoundary ( boundaryMask );

      ArrayType boundaryValues ( grid );
      setBoundaryValuesFundSol ( grid, boundaryValues );

      ArrayType sol ( grid );
      solvePoissonProblem ( grid, boundaryValues, boundaryMask, sol );
      ArrayType exactSol ( grid );
      computeExactSol ( grid, exactSol );

      exactSol -= sol;
      cerr << level << "\t" << exactSol.norm() * Sqr(grid.H()) << endl;
    }
    cerr << endl;

    if(success) {
      aol::printSelfTestSuccessMessage ( "--               QUOC Simplex Self Test Successful                            --" );
      aol::callSystemPauseIfNecessaryOnPlatform();
      return 0;
    }
  } // end try

  catch(std::exception &ex){
    cerr << "\n\nstd::exception caught:\n";
    cerr << ex.what () << endl;
  }
  catch(aol::Exception &ex){
    cerr << "\n\naol::Exception caught:\n";
    ex.dump ();
  }
  catch (...){
    cerr << "\n\nUnknown exception caught.\n";
  }
  aol::printSelfTestFailureMessage ( "!!               QUOC Simplex Self Test FAILED                                !!" );
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
