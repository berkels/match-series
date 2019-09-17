// include all quoc header files; in alphabetical order.
#include <AmbrosioTortorelli.h>
#include <anisoStiffOps.h>
#include <anisotropies.h>
#include <anisotropyVisualization.h>
#include <arrayExtensions.h>
#include <array.h>
#include <auxiliary.h>
#include <bitArray.h>
#include <boundaryIntegration.h>
#include <cellCenteredGrid.h>
#include <colorWheel.h>
#include <configurators.h>
#include <convolution.h>
#include <deformations.h>
#include <dm3Import.h>
#include <elastOps.h>
#include <elementMask.h>
#include <enoConvect.h>
#include <estimator2d.h>
#include <estimator3d.h>
#include <fastUniformGridMatrix.h>
#include <finiteDifferences.h>
#include <firstOrderTVAlgos.h>
#include <generator.h>
#include <gridBase.h>
#include <gridOp.h>
#include <gridSize.h>
#include <homogRWCWrapper.h>
#include <hyperelastic.h>
#include <imageTools.h>
#include <implicitSurfaceFEOps.h>
#include <indexMapper.h>
#include <isolineIterator2d.h>
#include <iterators.h>
#include <kernel2d.h>
#include <kernel3d.h>
#include <l2Projector.h>
#include <levelSetDrawer.h>
#include <levelSet.h>
#include <linearSmoothOp.h>
#include <mcm.h>
#include <morphology.h>
#include <multiDObject.h>
#include <multilevelArray.h>
#include <netCDFInterface.h>
#include <ocTree.h>
#include <paramReg.h>
#include <periodicBC.h>
#include <pngInterface.h>
#include <prolongation.h>
#include <qmElement.h>
#include <qmHeap.h>
#include <quadTree.h>
#include <quoc.h>
#include <quocDescent.h>
#include <quocMatrices.h>
#include <quocOps.h>
#include <quocTimestepSaver.h>
#include <rectangularGrid.h>
#include <registration.h>
#include <restriction.h>
#include <scalarArray.h>
#include <shapeLevelsetGenerator.h>
#include <simplexBaseFuncSetTFE.h>
#include <simplexBaseFunctionSet.h>
#include <simplexConfigurators.h>
#include <simplexGrid.h>
#include <simplexLookup.h>
#include <sweeping.h>
#include <tiledSpace.h>
#include <UGBMatrix.h>
#include <multiArray.h>
#include <Willmore.h>

int main( int, char** ) {

  try {
    // Macro Tests

    cerr << "--- Testing aol::ipfstream and ScalarArray<float, qc::QC_2D>::load ... \n";
    qc::ScalarArray<float, qc::QC_2D> a2d ("../../examples/testdata/image_129.pgm.bz2");
    cerr << "Test of aol::ipfstream ............................................... OK.\n";

    // Micro Tests

    bool success = true;

    {
      cerr << "--- Testing qc::ScalarArray<..., qc::QC_1D>::save and load ... ";

      qc::ScalarArray<unsigned char, qc::QC_1D> array1d ( 20 );
      for ( int i = 0; i < array1d.size(); ++i ) {
        array1d[i] = i*14;
      }
      array1d.save( "savetest.bz2", qc::PGM_UNSIGNED_CHAR_ASCII );
      qc::ScalarArray<unsigned char, qc::QC_1D> array1("savetest.bz2");
      remove( "savetest.bz2" );

      success &= ( array1d == array1 );

      cerr << success;

      array1d.save( "savetest.bz2", qc::PGM_UNSIGNED_CHAR_BINARY );
      qc::ScalarArray<unsigned char, qc::QC_1D> array2("savetest.bz2");
      remove( "savetest.bz2" );

      success &= ( array1d == array2 );

      cerr << success;

      qc::ScalarArray<float, qc::QC_1D> a1d ( a2d.getNumX() );
      for ( int i = 0; i < a1d.size(); ++i ) {
        a1d[i] = a2d.get( i, 0);
      }
      a1d.save( "savetest2.bz2", qc::PGM_FLOAT_ASCII );
      qc::ScalarArray<float, qc::QC_1D> array3("savetest2.bz2");
      remove( "savetest2.bz2" );

      array3 -= a1d;
      success &= ( array3.norm()/array3.size()<1e-9 );

      cerr << success;

      a1d.save( "savetest.bz2", qc::PGM_FLOAT_BINARY );
      qc::ScalarArray<float, qc::QC_1D> array4("savetest.bz2");
      remove( "savetest.bz2" );

      success &= ( a1d == array4 );

      cerr << success;

      qc::ScalarArray<double, qc::QC_1D> testArray ( 20 );
      for ( int i = 0; i < testArray.size(); ++i ) {
        testArray[i] = i*10.25;
      }
      testArray.save( "savetest.bz2", qc::PGM_DOUBLE_BINARY );
      qc::ScalarArray<double, qc::QC_1D> array5("savetest.bz2");
      remove( "savetest.bz2" );

      success &= ( testArray == array5 );

      cerr << success;

      if(success)
        cerr << "..... OK\n";
    }

    {
      cerr << "--- Testing qc::ScalarArray<..., qc::QC_2D>::save and load ... ";

      qc::ScalarArray<unsigned char, qc::QC_2D> array2d ("../../examples/testdata/image_129.pgm.bz2");

      array2d.save( "savetest.bz2", qc::PGM_UNSIGNED_CHAR_ASCII );
      qc::ScalarArray<unsigned char, qc::QC_2D> array1("savetest.bz2");
      remove( "savetest.bz2" );

      success &= ( array2d == array1 );

      cerr << success;

      array2d.save( "savetest.bz2", qc::PGM_UNSIGNED_CHAR_BINARY );
      qc::ScalarArray<unsigned char, qc::QC_2D> array2("savetest.bz2");
      remove( "savetest.bz2" );

      success &= ( array2d == array2 );

      cerr << success;

      a2d.save( "savetest2.bz2", qc::PGM_FLOAT_ASCII );
      qc::ScalarArray<float, qc::QC_2D> array3("savetest2.bz2");
      remove( "savetest2.bz2" );

      array3 -= a2d;
      success &= ( array3.norm()/array3.size()<1e-9 );

      cerr << success;

      a2d.save( "savetest.bz2", qc::PGM_FLOAT_BINARY );
      qc::ScalarArray<float, qc::QC_2D> array4("savetest.bz2");
      remove( "savetest.bz2" );

      success &= ( a2d == array4 );

      cerr << success;

      qc::GridDefinition grid (7, qc::QC_2D);
      qc::ScalarArray<double, qc::QC_2D> testArray (grid);
      typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature<double, qc::QC_2D, 3> > ConfType;
      qc::DataGenerator<ConfType> dataGenerator ( grid );
      dataGenerator.generateDiagonalLevelset ( testArray );
      testArray.save( "savetest.bz2", qc::PGM_DOUBLE_BINARY );
      qc::ScalarArray<double, qc::QC_2D> array5("savetest.bz2");
      remove( "savetest.bz2" );

      success &= ( testArray == array5 );

      cerr << success;

#ifdef USE_LIB_PNG
      a2d.savePNG( "savetest.png" );
      qc::ScalarArray<float, qc::QC_2D> array6;
      array6.setQuietMode ( true );
      // load will automatically call loadPNG because of the ".png" suffix.
      array6.load( "savetest.png" );
      remove( "savetest.png" );

      success &= ( a2d == array6 );

      cerr << success;
#endif
      if(success)
        cerr << "..... OK\n";

      ofstream myOut ( "savetest.txt" );
      {
        qc::ScalarArray<float, qc::QC_2D> rectArray ( 15, 27 ), otherRectArray ( 37, 13 );
        aol::NoiseOperator<float> nOp;
        nOp.applySingle ( rectArray );
        nOp.applySingle ( otherRectArray );
        myOut << rectArray << endl << otherRectArray << endl;
      }
      remove ( "savetest.txt" );
    }

    {
      cerr << "--- Testing qc::ScalarArray<..., qc::QC_3D>::save and load ... ";

      qc::ScalarArray<unsigned char, qc::QC_3D> array3d ( "../../examples/testdata/volume_9.dat.bz2" );

      array3d.save( "savetest.bz2", qc::PGM_UNSIGNED_CHAR_ASCII );
      qc::ScalarArray<unsigned char, qc::QC_3D> array1("savetest.bz2");
      remove( "savetest.bz2" );

      success &= ( array3d == array1 );

      cerr << success;

      array3d.save( "savetest.bz2", qc::PGM_UNSIGNED_CHAR_BINARY );
      qc::ScalarArray<unsigned char, qc::QC_3D> array2("savetest.bz2");
      remove( "savetest.bz2" );

      success &= ( array3d == array2 );

      cerr << success;

      qc::ScalarArray<float, qc::QC_3D> a3d ( "../../examples/testdata/volume_9.dat.bz2" );
      a3d.save( "savetest2.bz2", qc::PGM_FLOAT_ASCII );
      qc::ScalarArray<float, qc::QC_3D> array3("savetest2.bz2");
      remove( "savetest2.bz2" );

      array3 -= a3d;
      success &= ( array3.norm()/array3.size()<1e-9 );

      cerr << success;

      a3d.save( "savetest.bz2", qc::PGM_FLOAT_BINARY );
      qc::ScalarArray<float, qc::QC_3D> array4("savetest.bz2");
      remove( "savetest.bz2" );

      success &= ( a3d == array4 );

      cerr << success;

      qc::GridDefinition grid (5, qc::QC_3D);
      qc::ScalarArray<double, qc::QC_3D> testArray (grid);
      typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double, qc::QC_3D, 3> > ConfType;
      qc::DataGenerator<ConfType> dataGenerator ( grid );
      aol::Vec3<double> center (0.5, 0.5 , 0.5);
      dataGenerator.generateSphereLevelset ( testArray, center );
      testArray.save( "savetest.bz2", qc::PGM_DOUBLE_BINARY );
      qc::ScalarArray<double, qc::QC_3D> array5("savetest.bz2");
      remove( "savetest.bz2" );

      success &= ( testArray == array5 );

      cerr << success;


// does not work yet
/*      qc::ScalarArray<unsigned short, qc::QC_3D> array3dShort ( "../../examples/testdata/volume_9.dat.bz2" );
      for ( int i = 0; i < array3dShort.size(); ++i ) array3dShort[i] = 42023 * testArray[i];

      array3dShort.save( "savetest.bz2", qc::PGM_UNSIGNED_SHORT_BINARY );
      qc::ScalarArray<unsigned short, qc::QC_3D> array6("savetest.bz2");
      remove( "savetest.bz2" );

      success = success && ( array3dShort == array6 );

      cerr << success;
*/
      if(success)
        cerr << "..... OK\n";
    }

    {
      cerr << "Testing saveToFile, loadFromFile ... ";
      qc::ScalarArray<signed short, qc::QC_1D> array1d ( 20 ), array1dL;
      for ( int i = 0; i < array1d.size(); ++i ) {
        array1d[i] = i*14;
      }
      array1d.saveToFile ( "test.quo" );
      array1dL.loadFromFile ( "test.quo" );
      success &= ( array1d == array1dL );

      qc::ScalarArray<unsigned char, qc::QC_2D> array2d ("../../examples/testdata/image_129.pgm.bz2"), array2dL;
      array2d.saveToFile ( "test.quo" );
      array2dL.loadFromFile ( "test.quo" );
      success &= ( array2d == array2dL );

      qc::ScalarArray<double, qc::QC_3D> array3d ( "../../examples/testdata/volume_9.dat.bz2" ), array3dL;
      array3d.saveToFile ( "test.quo" );
      array3dL.loadFromFile ( "test.quo" );
      success &= ( array3d == array3dL );

      qc::MultiArray<double, qc::QC_3D, 3> multiArray2d3d ( qc::GridSize<qc::QC_3D>::createFrom ( array3d  ) ), multiArray2d3dL;
      multiArray2d3d[0] = array3d;
      multiArray2d3d[1] = array3d;
      multiArray2d3d[1] *= 42.0;
      multiArray2d3d.saveToFile ( "test.quo" );
      multiArray2d3dL.loadFromFile ( "test.quo" );
      success &= ( multiArray2d3d == multiArray2d3dL );

      remove ( "test.quo" );

      if(success) {
        cerr << "OK" << endl;
      }
    }

    {
      cerr << "--- Testing qc::Array<float> copy constructor ... " ;
      qc::Array<float> array_orig(3,3,3);
      array_orig.setZero();
      array_orig.set(1, 0, 2, 42.0);

      {
        qc::Array<float> array_default(array_orig);
        array_default.set(1, 0, 2, 23.0);

        success &= ( ( abs( array_orig.get(1, 0, 2) - 42.0 ) + abs( array_default.get(1,0,2) - 23.0) ) < 1e-6 );

        qc::Array<float> array_deep(array_orig, aol::DEEP_COPY );
        array_deep.set(1, 0, 2, 23.0);

        success &= ( ( abs( array_orig.get(1, 0, 2) - 42.0 ) + abs( array_deep.get(1, 0, 2) - 23.0) ) < 1e-6 );

        qc::Array<float> array_flat(array_orig, aol::FLAT_COPY );
        array_flat.set(1, 0, 2, 23.0);

        success &= ( ( abs( array_orig.get(1, 0, 2) - 23.0 ) + abs( array_flat.get(1, 0, 2) - 23.0) ) < 1e-6 );

        qc::Array<float> array_struct(array_orig, aol::STRUCT_COPY );
        array_struct.set(1, 0, 2, 17.0);

        success &= ( ( abs( array_orig.get(1, 0, 2) - 23.0 ) + abs( array_struct.get(1, 0, 2) - 17.0) ) < 1e-6 );

        array_default = *&array_default; // check whether self-assignment is treated correctly

      } // End of Scope - destroy all array copies

      aol::MemoryManager::deleteUnlocked();

      success &= ( abs(array_orig.get(1, 0, 2) - 23.0 ) < 1e-6 );

      if(success)
        cerr << "OK" << endl;
    }

    {
      cerr << "--- Testing qc::ScalarArray<float, qc::QC_2D> copy constructor ... " ;
      qc::ScalarArray<float, qc::QC_2D> array_orig(3,3);
      array_orig.setZero();
      array_orig.set(1, 2, 42.0);

      {
        qc::ScalarArray<float, qc::QC_2D> array_default(array_orig);
        array_default.set(1, 2, 23.0);

        success &= ( ( abs( array_orig.get(1, 2) - 42.0 ) + abs( array_default.get(1, 2) - 23.0) ) < 1e-6 );

        qc::ScalarArray<float, qc::QC_2D> array_deep(array_orig, aol::DEEP_COPY );
        array_deep.set(1, 2, 23.0);

        success &= ( ( abs( array_orig.get(1, 2) - 42.0 ) + abs( array_deep.get(1, 2) - 23.0) ) < 1e-6 );

        qc::ScalarArray<float, qc::QC_2D> array_flat(array_orig, aol::FLAT_COPY );
        array_flat.set(1, 2, 23.0);

        success &= ( ( abs( array_orig.get(1, 2) - 23.0 ) + abs( array_flat.get(1, 2) - 23.0) ) < 1e-6 );

        qc::ScalarArray<float, qc::QC_2D> array_struct(array_orig, aol::STRUCT_COPY );
        array_struct.set(1, 2, 17.0);

        success &= ( ( abs( array_orig.get(1, 2) - 23.0 ) + abs( array_struct.get(1, 2) - 17.0) ) < 1e-6 );

        array_default = *&array_default; // check self-assignment

      } // End of Scope - destroy all array copies

      aol::MemoryManager::deleteUnlocked();

      success &= ( abs(array_orig.get(1, 2) - 23.0 ) < 1e-6 );

      if(success)
        cerr << "OK" << endl;
    }

    { // copy constructors of qc::ScalarArray<float, qc::QC_3D>
      cerr << "--- Testing qc::ScalarArray<float, qc::QC_3D> copy constructor ... " ;
      qc::ScalarArray<float, qc::QC_3D> array_orig(3,3,3);
      array_orig.setZero();
      array_orig.set(1, 0, 2, 42.0);

      {
        qc::ScalarArray<float, qc::QC_3D> array_default(array_orig);
        array_default.set(1, 0, 2, 23.0);

        success &= ( ( abs( array_orig.get(1, 0, 2) - 42.0 ) + abs( array_default.get(1, 0, 2) - 23.0) ) < 1e-6 );

        qc::ScalarArray<float, qc::QC_3D> array_deep(array_orig, aol::DEEP_COPY );
        array_deep.set(1, 0, 2, 23.0);

        success &= ( ( abs( array_orig.get(1, 0, 2) - 42.0 ) + abs( array_deep.get(1, 0, 2) - 23.0) ) < 1e-6 );

        qc::ScalarArray<float, qc::QC_3D> array_flat(array_orig, aol::FLAT_COPY );
        array_flat.set(1, 0, 2, 23.0);

        success &= ( ( abs( array_orig.get(1, 0, 2) - 23.0 ) + abs( array_flat.get(1, 0, 2) - 23.0) ) < 1e-6 );

        qc::ScalarArray<float, qc::QC_3D> array_struct(array_orig, aol::STRUCT_COPY );
        array_struct.set(1, 0, 2, 17.0);

        success &= ( ( abs( array_orig.get(1, 0, 2) - 23.0 ) + abs( array_struct.get(1, 0, 2) - 17.0) ) < 1e-6 );

        array_default = *&array_default; // check self-assignment

      } // End of Scope - destroy all array copies

      aol::MemoryManager::deleteUnlocked();

      success &= ( abs(array_orig.get(1, 0, 2) - 23.0 ) < 1e-6 );

      if(success)
        cerr << "OK" << endl;
    }

    {
      cerr << "--- Testing Array classes' resize/reallocate ... " ;
      qc::ScalarArray<int, qc::QC_3D> array ( 3, 4, 5 );
      array.set ( 1, 2, 3, 42 );
      array.reallocate ( 2, 3, 4 );     success &= ( array.get ( 1, 2, 3 ) ==  0 ) && ( array.norm() == 0  );

      qc::RectangularContainer< int, qc::ScalarArray<int, qc::QC_3D>, qc::QC_3D > brick ( aol::Vec3<int>( 0, 0, 0), aol::Vec3<int>( 3, 4, 5) );
      brick.set ( 1, 2, 3, 42 );

      brick.resize ( aol::Vec3<int>( 0, 0, 0), aol::Vec3<int>( 2, 3, 4 ) );
      success &= ( brick.get ( 1, 2, 3 ) == 42 ) && ( brick.getContainedRef().norm() == 42 );

      brick.resize ( aol::Vec3<int>( 0, 0, 0), aol::Vec3<int>( 4, 5, 6 ) );
      success &= ( brick.get ( 1, 2, 3 ) == 42 ) && ( brick.getContainedRef().norm() == 42 );

      if(success)
        cerr << " OK" << endl;
    }

    { // int compatibility of arrays.
      cerr << "--- Testing qc::Array classes with size > 2^16 ... " ;
      cerr << "sizeof(short) = " << sizeof(short) << ", sizeof(int) = " << sizeof(int);

      qc::Array<float> array_big(70000,2);
      array_big.setZero();
      array_big.set(69000, 1, 42.0);
      array_big.set( 3464, 1, 23.0);
      success &= ( ( abs( array_big.get(69000, 1) - 42.0 ) + abs( array_big.get(3464, 1) - 23.0) ) < 1e-6 );

      qc::ScalarArray<float, qc::QC_2D> array_big2(70000, 2);
      array_big2.setZero();
      array_big2.set(69000, 1, 42.0);
      array_big2.set( 3464, 1, 23.0);
      success &= ( ( abs( array_big2.get(69000, 1) - 42.0 ) + abs( array_big2.get(3464, 1) - 23.0) ) < 1e-6 );

      qc::ScalarArray<float, qc::QC_3D> array_big3(70000, 2, 2);
      array_big3.setZero();
      array_big3.set(69000, 1, 1, 42.0);
      array_big3.set( 3464, 1, 1, 23.0);
      success &= ( ( abs( array_big3.get(69000, 1,1 ) - 42.0 ) + abs( array_big3.get(3464, 1, 1) - 23.0) ) < 1e-6 );

      if(success)
        cerr << " OK" << endl;
    }

    { // interpolate of qc::ScalarArray<float, qc::QC_3D>
      cerr << "--- Testing qc::ScalarArray<float, qc::QC_3D>::interpolate ... " ;
      qc::ScalarArray<float, qc::QC_3D> array(2,2,2);

      array.set( 0, 0, 0, 1.);
      array.set( 0, 0, 1, 2.);
      array.set( 0, 1, 0, 3.);
      array.set( 0, 1, 1, 4.);
      array.set( 1, 0, 0, 5.);
      array.set( 1, 0, 1, 6.);
      array.set( 1, 1, 0, 7.);
      array.set( 1, 1, 1, 8.);

      success &= ( array.interpolate( 0.5, 0.5, 0.5 ) == (1.+2.+3.+4.+5.+6.+7.+8.)/8.);

      if(success)
        cerr << "OK" << endl;
    }

    { // test some anisotropies and their visualization
      cerr << "--- Testing some anisotropies and their visualization ... " ;

      qc::L1Norm2d<double> L1Norm( 0.0000001 );
      qc::Rotated3dAnisotropy<double, qc::L1Norm2d<double> > aniso( L1Norm, 0.02 );

      qc::AnisotropyVisualizer3d<double, qc::Rotated3dAnisotropy<double, qc::L1Norm2d<double> > > visualizer( aniso, 32, 35, 43 );
      visualizer.setVerbose( false );
      visualizer.setSizeOfOutputImg( 3 );
      visualizer.generateWulffShapeByGammaZ3d( 0.1, 0.1 );
      visualizer.generateFrankDiagramByLevelLines3d();

      qc::AnisotropyVisualizer2d<double, qc::L1Norm2d<double> > visualizer2( L1Norm, 200, 100 );
      visualizer2.setVerbose( false );
      visualizer2.setSizeOfOutputImg( 5 );
      visualizer2.generateWulffShapeByGammaZ2d( 0.1, 0.1 );
      visualizer2.generateFrankDiagramByLevelLines2d();

      if (success) cerr << "OK" << endl;
    }

    { // test boundary integration
      cerr << "--- Testing IntegrateGradUOverBoundary ... " ;

      // testing 0 = integral over laplacae u for u(x,y) = x, which should be the same as
      // - integral \nabla u \nabla \vartheta + boundary-integral over \nabla u * normal \vartheta.
      // => First define the function x on a small grid.
      qc::GridDefinition grid( 4, qc::QC_2D );
      typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature<double,qc::QC_2D,7> > ConfType;
      qc::ScalarArray<double, qc::QC_2D> functionX( grid );    // the function ( u(x,y)=x )
      qc::ScalarArray<double, qc::QC_2D> temp( grid );

      for ( int i=0; i<functionX.getNumX(); i++ )     // filling the function with some senseful values
        for ( int j=0; j<functionX.getNumY(); j++ )
          functionX.set( i,j, static_cast<double>( i ) / static_cast<double>( functionX.getNumX() ) );

      // the operators
      aol::StiffOp<ConfType> stiffOp( grid );
      qc::IntegrateGradUOverBoundary< ConfType, aol::GaussQuadrature< double, qc::QC_1D, 7 > > boundaryOp( grid );

      stiffOp.apply( functionX, temp );
      temp *= -1.;
      boundaryOp.applyAdd( functionX, temp );     // temp should now contain laplace u which is 0
      success &= ( temp.norm() < 1e-10 );

      if (success) cerr << "OK" << endl;
    }

    { // AArray and RectangularContainer
      cerr << "--- Testing AArray and RectangularContainer ... ";
      qc::RectangularContainer< double, qc::ScalarArray<double, qc::QC_3D>, qc::QC_3D > pa ( aol::Vec3<int>(2, 2, 2), aol::Vec3<int>(4, 4, 4) );
      aol::Vec3<int> posn( 2, 2, 3 );

      pa.set ( posn, 2.0 );
      success &= fabs ( pa.get( posn ) - 2.0 ) < 1e-6;
      pa.resize( aol::Vec3<int>(1, 1, 1), aol::Vec3<int>(5, 6, 7) );
      success &= fabs ( pa.get( posn ) - 2.0 ) < 1e-6;
      pa.resize( aol::Vec3<int>(2, 2, 2), aol::Vec3<int>(4, 4, 4) );
      success &= fabs ( pa.get( posn ) - 2.0 ) < 1e-6;
      pa.resize( aol::Vec3<int>(8, 8, 8), aol::Vec3<int>(10, 11, 19) );
      pa.resize( aol::Vec3<int>(2, 2, 2), aol::Vec3<int>(4, 4, 4) );
      success &= fabs ( pa.get( posn ) - 0.0 ) < 1e-6;

      qc::AArray< aol::Matrix33<double>, qc::QC_3D > aa ( 4, 5, 6 );
      aa.getRef(2, 2, 3).setIdentity();
      aa.getRef(2, 2, 3).set( 1, 0, 2.0 );
      aa.getRef(2, 2, 1).makeInverse( aa.getRef(2,2,3) );

      success &= fabs ( aa.getRef(2, 2, 1).get(1,0) + 2.0 ) < 1e-6 ;

      if (success) cerr << "OK" << endl;
    }

    {
      cerr << "--- Testing Quoc restriction and prolongation ... ";
      qc::GridDefinition cGrid( 2, qc::QC_3D ), fGrid  ( 3, qc::QC_3D );
      {
        qc::ProlongOp< double > prOp ( cGrid, fGrid, aol::ONTHEFLY );
        qc::ProlongOp< double > prOpA ( cGrid, fGrid, aol::ASSEMBLED );
        success &= compareOps ( prOp, prOpA, cGrid.getNumberOfNodes(), fGrid.getNumberOfNodes() );
        qc::ScalarArray<double, qc::QC_3D> fArr ( fGrid ), cArr ( cGrid );
        cArr.setAll ( 1.0 );
        prOp.apply ( cArr, fArr );
        success &= ( fabs ( 1 - fArr.getMinValue() ) + fabs ( 1 - fArr.getMaxValue() ) < 1e-6 );
        cerr << "Prolong ... ";
      }
      {
        qc::RestrictOp< double, qc::STD_MG_RESTRICT > reOpM ( cGrid, fGrid, aol::ONTHEFLY );
        qc::RestrictOp< double, qc::STD_MG_RESTRICT > reOpMA ( cGrid, fGrid, aol::ASSEMBLED );
        success &= compareOps ( reOpM, reOpMA, fGrid.getNumberOfNodes(), cGrid.getNumberOfNodes() );
        qc::ScalarArray<double, qc::QC_3D> fArr ( fGrid ), fMOne ( fGrid ), cArr ( cGrid ), cOne ( cGrid );
        typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D,3> > ConfigType;
        aol::MassOp<ConfigType> massOp( fGrid );
        fArr.setAll ( 1.0 ); cOne.setAll ( 1.0 );
        massOp.apply ( fArr, fMOne );
        reOpM.apply  ( fMOne, cArr );
        success &= ( fabs ( 1 - cArr * cOne ) < 1e-6 );
        cerr << "StdMGRestict ... ";
      }
      {
        qc::RestrictOp< double, qc::STD_QUOC_RESTRICT > reOpQ ( cGrid, fGrid, aol::ONTHEFLY );
        qc::RestrictOp< double, qc::STD_QUOC_RESTRICT > reOpQA ( cGrid, fGrid, aol::ASSEMBLED );
        success &= compareOps ( reOpQ, reOpQA, fGrid.getNumberOfNodes(), cGrid.getNumberOfNodes() );
        qc::ScalarArray<double, qc::QC_3D> fArr ( fGrid ), cArr ( cGrid );
        fArr.setAll ( 1.0 );
        reOpQ.apply ( fArr, cArr );
        success &= ( fabs ( 1 - cArr.getMinValue() ) + fabs ( 1 - cArr.getMaxValue() ) < 1e-6 );
        cerr << "StdQuocRestrict ...";
      }
      {
        qc::RestrictOp< double, qc::THROW_AWAY_RESTRICT > reOpT ( cGrid, fGrid, aol::ONTHEFLY );
        qc::RestrictOp< double, qc::THROW_AWAY_RESTRICT > reOpTA ( cGrid, fGrid, aol::ASSEMBLED );
        success &= compareOps ( reOpT, reOpTA, fGrid.getNumberOfNodes(), cGrid.getNumberOfNodes() );
        qc::ScalarArray<double, qc::QC_3D> fArr ( fGrid ), cArr ( cGrid );
        fArr.setAll ( 1.0 );
        reOpT.apply ( fArr, cArr );
        success &= ( fabs ( 1 - ( cArr.sum() / cArr.size() ) ) < 1e-6 );
        cerr << "ThrowAwayRestrict ...";
      }
      if ( success )
        cerr << "OK" << endl;
    }

    {
      cerr << "--- Testing OTFILexMapper and FastILexMapper ... ";
      {
        qc::ScalarArray<int, qc::QC_2D> rArr ( 4, 7 );
        qc::FastILexMapper<qc::QC_2D> flm2d ( rArr );
        qc::OTFILexMapper<qc::QC_2D> olm2d ( rArr );
        success &= ( flm2d.getGlobalIndex ( 3, 5 ) == 23 ) && ( olm2d.getGlobalIndex ( 3, 5 ) == 23 );
        int x, y;
        flm2d.splitGlobalIndex ( 23, x, y );
        success &= ( x == 3 && y == 5 );
        olm2d.splitGlobalIndex ( 23, x, y );
        success &= ( x == 3 && y == 5 );
      }
      {
        qc::RectangularGrid<qc::QC_3D> rGrid ( aol::Vec3<int> ( 4, 7, 3 ) );
        qc::FastILexMapper<qc::QC_3D> flm3d ( rGrid );
        qc::OTFILexMapper<qc::QC_3D> olm3d ( rGrid );
        success &= ( flm3d.getGlobalIndex ( 2, 5, 1 ) == 50 ) && ( olm3d.getGlobalIndex ( 2, 5, 1 ) == 50 );
        int x, y, z;
        flm3d.splitGlobalIndex ( 50, x, y, z );
        success &= ( x == 2 && y == 5 && z == 1 );
        olm3d.splitGlobalIndex ( 50, x, y, z );
        success &= ( x == 2 && y == 5 && z == 1 );
      }
      if ( success )
        cerr << "OK" << endl;
    }

    {
      cerr << "--- qc::FastUniformGridMatrix< double, qc::QC_2D >::transposeTo ... ";
      {
        typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature<double,qc::QC_2D, 3> > ConfType;
        qc::GridDefinition grid( 4, qc::QC_2D );
        qc::FastUniformGridMatrix< double, qc::QC_2D > A ( grid );
        qc::FastUniformGridMatrix< double, qc::QC_2D > B ( grid );
        qc::DataGenerator< ConfType > generator( grid  );
        qc::ScalarArray< double, qc::QC_2D > image ( grid );
        generator.generateCircleLevelset( image, 0.25 );

        // Use a non symmetric operator here to generate a non symmetric matrix.
        aol::ImageGradientSemiDiffOp< ConfType > imageGradientSemiDiffOp( grid, image, false );
        imageGradientSemiDiffOp.assembleAddMatrix( A );
        A.transposeTo( B );

        aol::ImageGradientSemiDiffOp< ConfType > imageGradientSemiDiffOpTransposed( grid, image, true );
        A.setZero();
        imageGradientSemiDiffOpTransposed.assembleAddMatrix( A );
        success &= A.isApproxEqual( B, 1E-16 );;
      }

      if ( success )
        cerr << "OK" << endl;
    }

    {
      cerr << "--- RectangularGrid<qc::QC_2D>::MatrixType ... ";

      qc::RectangularGrid<qc::QC_2D> rGrid ( aol::Vec3<int>( 23, 42, 1 ) );
      typedef qc::RectangularGridConfigurator< double, qc::QC_2D, aol::GaussQuadrature<double,qc::QC_2D, 3> > ConfType;
      aol::MassOp< ConfType > massOp ( rGrid );
      ConfType::MatrixType bandMat ( qc::GridSize<qc::QC_2D>::createFrom ( rGrid ) );
      aol::SparseMatrix<double> sparseMat ( rGrid );
      massOp.assembleAddMatrix ( bandMat );
      massOp.assembleAddMatrix ( sparseMat );
      sparseMat -= bandMat;

      success &= ( sparseMat.getFrobeniusNormSqr() < 1e-14 );

      if ( success )
        cerr << "OK" << endl;
    }

    {
      cerr << "--- RectangularGrid<qc::QC_3D>::MatrixType ... ";

      qc::RectangularGrid<qc::QC_3D> rGrid ( aol::Vec3<int>( 13, 20, 5 ) );
      typedef qc::RectangularGridConfigurator< double, qc::QC_3D, aol::GaussQuadrature<double,qc::QC_3D, 3> > ConfType;
      aol::MassOp< ConfType > massOp ( rGrid );
      qc::MultilinFEBandMatrix<double,qc::QC_3D> bandMat ( qc::GridSize<qc::QC_3D>::createFrom ( rGrid ) );
      aol::SparseMatrix<double> sparseMat ( rGrid );
      massOp.assembleAddMatrix ( bandMat );
      massOp.assembleAddMatrix ( sparseMat );
      sparseMat -= bandMat;

      success &= ( sparseMat.getFrobeniusNormSqr() < 1e-14 );

      if ( success )
        cerr << "OK" << endl;
    }

    {
      cerr << "--- CubicGrid<>::Iterators ... ";
      {
        qc::CubicGrid<qc::QC_2D> cgrid ( 2 );
        qc::BitArray<qc::QC_2D> check ( qc::GridSize<qc::QC_2D>::createFrom ( cgrid ) );

        check.setAll ( false );
        for ( qc::CubicGrid<qc::QC_2D>::FullElementIterator elit ( cgrid ); elit.notAtEnd(); ++elit ) {
          check.set ( *elit, true );
        }
        success &= ( check.crc32() == 3529542153UL ); // make sure  this is interpreted as unsigned long integer
        check.setAll ( false );

        for ( qc::CubicGrid<qc::QC_2D>::FullNodeIterator nit ( cgrid ); nit.notAtEnd(); ++nit ) {
          check.set ( *nit, true );
        }
        success &= ( check.crc32() == 1510334235UL );
        check.setAll ( false );

        for ( qc::CubicGrid<qc::QC_2D>::FullBoundaryNodeIterator bnit ( cgrid ); bnit.notAtEnd(); ++bnit ) {
          check.set ( *bnit, true );
        }
        success &= ( check.crc32() == 1219721669UL );
        check.setAll ( false );
      }

      {
        qc::CubicGrid<qc::QC_3D> cgrid ( 2 );
        qc::BitArray<qc::QC_3D> check ( qc::GridSize<qc::QC_3D>::createFrom ( cgrid ) );

        for ( qc::CubicGrid<qc::QC_3D>::FullElementIterator elit ( cgrid ); elit.notAtEnd(); ++elit ) {
          check.set ( *elit, true );
        }
        success &= ( check.crc32() == 413096062UL);
        check.setAll ( false );

        for ( qc::CubicGrid<qc::QC_3D>::FullNodeIterator nit ( cgrid ); nit.notAtEnd(); ++nit ) {
          check.set ( *nit, true );
        }
        success &= ( check.crc32() == 1615256477UL );
        check.setAll ( false );

        for ( qc::CubicGrid<qc::QC_3D>::FullBoundaryNodeIterator bnit ( cgrid ); bnit.notAtEnd(); ++bnit ) {
          check.set ( *bnit, true );
        }
        success &= ( check.crc32() == 857664509UL );
      }

      if ( success )
        cerr << "OK" << endl;
    }

    {
      cerr << "--- Testing SymmetricProjector ... " ;

      typedef qc::QuocConfiguratorTraitMultiLin<double, qc::QC_2D, aol::GaussQuadrature<double,qc::QC_2D,3> > ConfType;
      ConfType::InitType grid( 4, qc::QC_2D );
      qc::DataGenerator<ConfType> dataGenerator ( grid );
      qc::MultiArray<double, 2> phi ( grid );
      dataGenerator.generateSkewSymmetricDeformation ( .6, phi );
      qc::SymmetricProjector<ConfType> projector ( grid );
      projector.project ( phi );
      success &= ( phi.getMaxOfPointWiseNorm() < 1.e-6 );

      if (success) cerr << "OK" << endl;
    }

    {
      cerr << "--- Testing UniformGridMassOp in 2D ... " ;

      typedef qc::RectangularGridConfigurator< double, qc::QC_2D, aol::GaussQuadrature<double,qc::QC_2D, 3> > ConfType;
      const qc::RectangularGrid<qc::QC_2D> grid ( aol::Vec3<int>( 23, 12, 1 ) );
      const qc::UniformGridMassOp<double, qc::QC_2D> uMass ( qc::GridSize2d::createFrom( grid ) );
      const aol::MassOp<ConfType> mass ( grid, aol::ASSEMBLED );
      aol::Vector<double> a ( grid ), b ( grid ), c ( grid );
      // Apply both mass ops to all canonical basis vectors and compare the results.
      for ( int i = 0; i < a.size(); ++i ) {
        a.setZero();
        a[i] = 1;
        uMass.apply ( a, b );
        mass.apply ( a, c );
        c-=b;
        success &= ( c.normSqr() < 1.e-32 );
      }

      if (success) cerr << "OK" << endl;
    }

    {
      cerr << "--- aol::StiffOp on RectangularGrid<qc::QC_2D> ... ";
      typedef qc::RectangularGridConfigurator< double, qc::QC_2D, aol::GaussQuadrature<double,qc::QC_2D, 3> > ConfType;
      ConfType::InitType grid ( aol::Vec3<int> ( 8, 4, 1 ) );
      aol::StiffOp<ConfType> stiffOp( grid );
      aol::StiffOp<ConfType> stiffOpAssem ( grid, aol::ASSEMBLED );
      success &= compareOps( stiffOp, stiffOpAssem, grid.getNumberOfNodes(), grid.getNumberOfNodes(), 1e-8 );

      if (success) cerr << "OK" << endl;
    }

    if(success) {
      aol::printSelfTestSuccessMessage ( "--                       QUOC Self Test Successful                            --" );

      cerr << aol::color::reset;
      aol::callSystemPauseIfNecessaryOnPlatform();
      return 0;
    }

    {
      // this should not compile if someone is "using namespace aol;" in the headers included.
      using namespace qc;
      aol::Vector<avoid_namespace_collision> test_vector(1);
    }

  }
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

  aol::printSelfTestFailureMessage ( "!!                       QUOC Self Test FAILED                                !!" );
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 1;
}
