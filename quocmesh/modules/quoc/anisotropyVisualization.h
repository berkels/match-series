#ifndef __ANISOTROPYVISUALIZATION_H
#define __ANISOTROPYVISUALIZATION_H

#include <anisotropies.h>
#include <scalarArray.h>
#include <aol.h>
#include <qmException.h>
#include <levelSet.h>

namespace qc {

//! Operator for overloading the velocity function of the Enquist-Osher scheme
/*!
 * \author Nemitz
*/
template <typename RealType, typename AnisoType>
    class GrowWulffshapeEO3D : public qc::LevelSet3dInt<RealType, GrowWulffshapeEO3D<RealType, AnisoType> > {
  const qc::GridDefinition &_grid;
  const AnisoType &_anisotropy;

public:
  GrowWulffshapeEO3D( const qc::GridDefinition &Grid, const AnisoType &Anisotropy ) :
    qc::LevelSet3dInt< RealType, GrowWulffshapeEO3D<RealType, AnisoType> >(  ),
    _grid( Grid ), _anisotropy( Anisotropy ) {  }


  RealType velocity( const qc::ScalarArray<RealType, qc::QC_3D> &Data, int x, int y, int z ) const {
        aol::Vec3<RealType> grad;
        Data.gradientFD( x,y,z, grad );
        grad /= sqrt( grad.normSqr() + aol::Sqr(_grid.H()) );
        return _anisotropy.gammaNorm( grad );
  }

};



/**
 * A class with different methods for visualizing the Wulff shape or the Frank diagram
 * to a given anisotropy.
 * ATTENTION: The anisotropy must at least have a method "gammaNorm(...)", furthermore
 * a method "gammaFirstDerivative(...)" if you want to use gammaZSphereToWulffshape(...)
 * (see the anisotropies from anisotropies.h for examples).
 * It is assumed, that the level-lines of \f$ \gamma \f$ are shapes of the Frank diagram,
 * that means \f$ \gamma \f$ is its distance function and at the same time the support-
 * function of the Wulff shape.
 *
 * Author: Nemitz
 */
template <typename RealType, typename AnisoType>
class AnisotropyVisualizer2d {
private:
  const AnisoType &_anisotropy;
  int _nx;                                //! size of the OutputImg-array
  int _ny;
  bool _verbose;                          //! do some output
  ScalarArray<RealType, qc::QC_2D> _outputImg;

public:
  // ------------------------------------------------ methods -----------------------------------------------

  //! The constructor expects an anisotropy
  AnisotropyVisualizer2d ( const AnisoType &Anisotropy ) :
      _anisotropy ( Anisotropy ), _verbose ( true ), _outputImg ( 129, 129 ) {  }

  //! The constructor expects an anisotropy and the size of the output
  AnisotropyVisualizer2d ( const AnisoType &Anisotropy, int Nx, int Ny ) :
      _anisotropy ( Anisotropy ), _verbose ( true ), _outputImg ( Nx, Ny ) { _nx = Nx; _ny = Ny; }

  //! The constructor expects an anisotropy and the level of the output
  AnisotropyVisualizer2d ( const AnisoType &Anisotropy, int Level ) :
      _anisotropy ( Anisotropy ), _verbose ( true ), _outputImg ( 1, 1 ) { setSizeOfOutputImg ( Level ); }

  void setSizeOfOutputImg ( int Level ) {
    _nx = _ny = ( 2 << ( Level - 1 ) ) + 1;
    _outputImg.reallocate ( _nx, _ny );
  }

  void setVerbose ( bool Verbose ) { _verbose = Verbose; }

  //! This function simply saves the outputImg
  void saveImage ( const char* fileName, SaveType type = qc::PGM_UNSIGNED_CHAR_BINARY,
                   char* comment = NULL ) {
    _outputImg.setOverflowHandling ( aol::SCALE, 0., 1. );
    _outputImg.save ( fileName, type, comment );
  }

  //! This function provides the outputImg
  ScalarArray<RealType, qc::QC_2D>& getImageReference() { return _outputImg; }

  void copyImage( ScalarArray<RealType, qc::QC_2D> &copy ) {
    if ( copy.size() != _outputImg.size() )
      throw aol::Exception ( "qc::AnisotropyVisualizer2d: Argument of copyImage must have same size as the visualizer-image!!", __FILE__, __LINE__ );
    for (int i=0; i<_outputImg.size(); i++)
      copy[i] = _outputImg[i];
  }


  //! Visualization-methods: I) ******* The Frank diagram *********

  //! The Frank diagram is given by the level lines of \f$ \gamma \f$.
  //! This function only generates the Frank diagram in a given ScalarArray<QC_2D>, it will NOT be saved!
  void generateFrankDiagramByLevelLines2d ( void ) {
    if ( _verbose ) cerr << "Generating the Frank diagram just by evaluating gamma ...\n";
    // now traverse all nodes of the image
    for ( int i = 0; i < _nx; i++ ) {
      for ( int j = 0; j < _ny; j++ ) {
        // first: scale to [-1,1]^2
        RealType x = 2. * static_cast<RealType> ( i ) / static_cast<RealType> ( _nx ) - 1.;
        RealType y = 2. * static_cast<RealType> ( j ) / static_cast<RealType> ( _ny ) - 1.;

        // evaluate the anisotropy
        aol::Vec2<RealType> arg ( x, y );
        RealType val = _anisotropy.gammaNorm ( arg );
        _outputImg.set ( i, j, val );
      }
    }
  }   // end method


  //! Visualization-methods: II) ******* The Wulff shape *********

  //! Generate the WulffShape by \f$ W=\gamma_z(n), n \in S^1 \f$
  //! NOTICE: The anisotropy must provide a function gammaFirstDerivative( z,v ).
  void generateWulffShapeByGammaZ2d ( double stepGetMax = 0.01, double traverseStep = 0.000001 ) {
    if ( _verbose ) cerr << "Generating the Wulff shape by applying gamma_z to the sphere ...\n";
    _outputImg.setAll ( 0. );


    aol::Vec2<double> z ( 1., 0. );    // argument of gamma_z
    aol::Vec2<double> v ( 0., 0. );    // result = derivative

    // first traversing: Get the maximum of the derivative.
    // This is necessary to fit the Wulff shape into the image.
    if ( _verbose ) cerr << "Compute the max-norm of the derivative...";
    double maxNorm = 0.;
    for ( double phi = 0.; phi < 6.3; phi += stepGetMax ) {
      z[0] = cos ( phi );              // first: compute the normal
      z[1] = sin ( phi );

      // evaluate the derivative of the anisotropy
      _anisotropy.gammaFirstDerivative ( z, v );

      if ( v.norm() > maxNorm ) maxNorm = v.norm();
    }

    // now traverse the sphere to draw the Wulff shape
    if ( _verbose ) cerr << "now traverse the sphere...";
    for ( double phi = 0.; phi < 6.3; phi += traverseStep ) {
      // first: compute the normal
      z[0] = cos ( phi );
      z[1] = sin ( phi );

      // evaluate the derivative of the anisotropy
      _anisotropy.gammaFirstDerivative ( z, v );

      v[0] = static_cast<double> ( _nx ) / 3. * ( v[0] / maxNorm );
      v[1] = static_cast<double> ( _ny ) / 3. * ( v[1] / maxNorm );

      int i = _nx / 2 + static_cast<int> ( v[0] );
      int j = _ny / 2 + static_cast<int> ( v[1] );
      _outputImg.set ( i, j, 1. );
    }
    if ( _verbose ) cerr << "done!\n";
  }


  //! Generate the WulffShape by evaluating \f$ \sup \frac{\langle z,n \rangle }{\gamma_z(n)} \f$
  //! NOTICE: This method takes quite a long time, but needs only the gammaNorm-function and not
  //! the derivative.
  //! Originally written by Berkels, adjusted to this class by Nemitz

  //! auxiliary function for the generating
  RealType evaluateDualMappingOfAnisotropy( const aol::Vec2<RealType> Position,
                                            const qc::ScalarArray<RealType, qc::QC_2D> &WulffShape,
                                            double sizeOfSquare = 1.5 ) const {
    aol::Vec2<RealType> currentPosition;
    RealType value = 0.;
    RealType tmp;
    const int num = WulffShape.getNumX();
    if ( num != WulffShape.getNumY() )
      throw aol::Exception( "qc::AnisotropyVisualization: evaluateDualMappingOfAnisotropy needs square ScalarArray<QC_2D> WulffShape", __FILE__, __LINE__);
    const RealType h = 1./( static_cast<RealType>(num - 1) );
    for ( int X = 0; X < num; X++ ) {
      for ( int Y = 0; Y < num; Y++ ) {
        currentPosition[0] = ( X * h - 0.5 ) * sizeOfSquare;      // Werte von 0..num-1 nach -1,1 skalieren
        currentPosition[1] = ( Y * h - 0.5 ) * sizeOfSquare;
        if ( WulffShape[ X + Y*num ] == 0.){
          tmp = (currentPosition*Position);
          if( value < tmp )
            value = tmp;
        }
      }
    }
    return value;
  }

  //! The generating function
  void generateWulffShapeBySupFormula2d( double sizeOfSquare = 2. ) {
    if ( _verbose ) cerr << "Generating Wulff shape by evaluating the sup-formula. This will take quite some time.\n";
    if ( _nx != _ny )
      throw aol::Exception( "qc::AnisotropyVisualization: generateWulffShapeBySupFormula needs a quadratic image!", __FILE__, __LINE__);
    qc::ScalarArray<RealType, qc::QC_2D> wulffShape( _nx, _ny );
    qc::ScalarArray<RealType, qc::QC_2D> frankDiagram( _nx, _ny );

    // generate the Frank diagram (here called Wulff shape)
    if ( _verbose ) cerr << "Generating Frank diagram...\n";
    const RealType h = 1. / ( static_cast<RealType>( _nx - 1 ) );
    aol::Vec2<RealType> position;
    for ( int X = 0; X < _nx; X++ ) {
      for ( int Y = 0; Y < _ny; Y++ ) {
        position[0] = (X * h - 0.5) * sizeOfSquare;        // Werte von 0..num-1 nach -1,1 skalieren
        position[1] = (Y * h - 0.5) * sizeOfSquare;
        if( _anisotropy.gammaNorm( position) > 1.)
          wulffShape[ X + Y*_nx ] = 1.;
      }
    }
//     wulffShape.save("wulffshape.bz2",qc::PGM_DOUBLE_BINARY);

    // now evaluate the sup-formula for each pixel
    if ( _verbose ) cerr << "Evaluating sup-formula...\n";
    for ( int X = 0; X < _nx; X++ ) {
      for ( int Y = 0; Y < _ny; Y++ ) {
        position[0] = (X * h - 0.5) * sizeOfSquare;        // Werte von 0..num-1 nach -1,1 skalieren
        position[1] = (Y * h - 0.5) * sizeOfSquare;
        if( evaluateDualMappingOfAnisotropy( position, wulffShape, sizeOfSquare ) > 1.)
          _outputImg[ X + Y*_nx ] = 1.;
      }
      if ( _verbose ) cerr << static_cast<RealType>(X)/static_cast<RealType>(_nx-1)*100 << " % complete    \r";
    }
  }

};




// ---------------------------------------------------------------------------------------------------------------
// ---------------------------------- 3D - 3D - 3D - 3D --------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------


template <typename RealType, typename AnisoType>
class AnisotropyVisualizer3d {
private:
  const AnisoType &_anisotropy;
  int _nx;                                //! size of the OutputImg-array
  int _ny;
  int _nz;
  bool _verbose;                          //! do some output
  ScalarArray<RealType, qc::QC_3D> _outputImg;     //! contains the data

public:

  // ------------------------------------------------ methods -----------------------------------------------

  //! The constructor expects an anisotropy
  AnisotropyVisualizer3d ( const AnisoType &Anisotropy ) :
      _anisotropy ( Anisotropy ), _verbose ( true ), _outputImg ( 65, 65, 65 ) {  }

  //! The constructor expects an anisotropy and the size of the output
  AnisotropyVisualizer3d ( const AnisoType &Anisotropy, int Nx, int Ny, int Nz ) :
      _anisotropy ( Anisotropy ), _verbose ( true ), _outputImg ( Nx, Ny, Nz ) { _nx = Nx; _ny = Ny; _nz = Nz; }

  //! The constructor expects an anisotropy and the level of the output
  AnisotropyVisualizer3d ( const AnisoType &Anisotropy, int Level ) :
      _anisotropy ( Anisotropy ), _verbose ( true ), _outputImg ( 1, 1, 1 ) { setSizeOfOutputImg ( Level ); }

  void setSizeOfOutputImg ( int Level ) {
    _nx = _ny = _nz = ( 2 << ( Level - 1 ) ) + 1;
    _outputImg.reallocate ( _nx, _ny, _nz );
  }

  void setVerbose ( bool Verbose ) { _verbose = Verbose; }

  //! This function simply saves the outputImg
  void saveImage ( const char* fileName, SaveType type = qc::PGM_FLOAT_BINARY,
                   char* comment = NULL ) const {
    _outputImg.save ( fileName, type, comment );
  }

  ScalarArray<RealType, qc::QC_3D>& getImageReference() { return _outputImg; }

  void copyImage( ScalarArray<RealType, qc::QC_3D> &copy ) {
    if ( copy.size() != _outputImg.size() )
      throw aol::Exception ( "qc::AnisotropyVisualizer3d: Argument of copyImage must have same size as the visualizer-image!!", __FILE__, __LINE__ );
    for (int i=0; i<_outputImg.size(); i++)
      copy[i] = _outputImg[i];
  }


  //! Visualization-methods: I) ******* The Frank diagram *********

  //! This function only generates the Frank diagram in a given ScalarArray<QC_3D>, it will NOT be saved!
  void generateFrankDiagramByLevelLines3d ( void )  {
    if ( _verbose ) cerr << "Generating the Frank diagram just by evaluating gamma ...\n";
    // traverse all nodes of the image
    for ( int i = 0; i < _nx; i++ ) {
      for ( int j = 0; j < _ny; j++ ) {
        for ( int k = 0; k < _nz; k++ ) {
          // first: scale to [-1,1]^3
          RealType x = 2. * static_cast<RealType> ( i ) / static_cast<RealType> ( _nx ) - 1.;
          RealType y = 2. * static_cast<RealType> ( j ) / static_cast<RealType> ( _ny ) - 1.;
          RealType z = 2. * static_cast<RealType> ( k ) / static_cast<RealType> ( _nz ) - 1.;

          // evaluate the anisotropy
          aol::Vec3<RealType> arg ( x, y, z );
          RealType val = _anisotropy.gammaNorm ( arg );
          _outputImg.set ( i, j, k,  val );
        }
      }
    }
  }   // end method


  //! Visualization-methods: II) ******* The Wulff shape *********

  //! Generate the WulffShape by \f$ W=\gamma_z(n), n \in S^1 \f$
  //! NOTICE: The anisotropy must provide a function gammaFirstDerivative( z,v ).
  void generateWulffShapeByGammaZ3d ( double getMaxStep = 0.01, double traverseStep = 0.0005 ) {
    if ( _verbose ) cerr << "Generating the Wulff shape by applying gamma_z to the sphere ...\n";
    _outputImg.setAll ( 0. );

    // first traversing: Get the maximum of the derivative.
    // This is necessary to fit the Wulff shape into the image.
    if ( _verbose ) cerr << "Compute the max-norm of the derivative...";
    double maxNorm = 0.;
    aol::Vec3<double> v;
    for ( double phi = -aol::NumberTrait<double>::pi; phi <= aol::NumberTrait<double>::pi + 0.1; phi += getMaxStep ) {
      if ( _verbose ) cerr << "\rComputing " << ( phi + aol::NumberTrait<double>::pi ) *100. / ( 2.*aol::NumberTrait<double>::pi + 0.1 ) << " % ...";
      for ( double theta = 0.; theta <= 4.*aol::NumberTrait<double>::pi; theta += getMaxStep ) {
        // first: compute the normal
        aol::Vec3<double> z ( cos ( phi ) *cos ( theta ), cos ( phi ) *sin ( theta ), sin ( phi ) );

        // evaluate the derivative of the anisotropy
        _anisotropy.gammaFirstDerivative ( z, v );

        if ( v.norm() > maxNorm ) maxNorm = v.norm();
      }
    }

    // now traverse the sphere to draw the Wulff shape
    if ( _verbose ) cerr << "now traverse the sphere...\n";
    for ( double phi = -aol::NumberTrait<double>::pi; phi <= aol::NumberTrait<double>::pi + 0.1; phi += traverseStep ) {
      if ( _verbose ) cerr << "\rComputing " << ( phi + aol::NumberTrait<double>::pi ) *100. / ( 2.*aol::NumberTrait<double>::pi + 0.1 ) << " % ...";
      for ( double theta = 0.; theta <= 4.*aol::NumberTrait<double>::pi; theta += traverseStep ) {
        aol::Vec3<double> z ( cos ( phi ) *cos ( theta ), cos ( phi ) *sin ( theta ), sin ( phi ) );


        // evaluate the derivative of the anisotropy
        _anisotropy.gammaFirstDerivative ( z, v );

        v[0] = static_cast<double> ( _nx ) / 3. * ( v[0] / maxNorm );
        v[1] = static_cast<double> ( _ny ) / 3. * ( v[1] / maxNorm );
        v[2] = static_cast<double> ( _nz ) / 3. * ( v[2] / maxNorm );

        int i = _nx / 2 + static_cast<int> ( v[0] );
        int j = _ny / 2 + static_cast<int> ( v[1] );
        int k = _nz / 2 + static_cast<int> ( v[2] );
        _outputImg.set ( i, j, k, 1. );
      }
    }
    if ( _verbose ) cerr << "done!\n";
  }


  //! Generate the WulffShape by \f$ W=\gamma_z(n), n \in S^1 \f$
  //! NOTICE: The anisotropy must provide a function gammaFirstDerivative( z,v ).
  void generateWulffShapeByEnquistOsherScheme3d ( int numSteps = 30, double tau = 0.5 ) {
    if ( _verbose ) cerr << aol::color::red << "\nGenerating the Wulff shape by applying an Enquist Osher scheme ...\n";
    _outputImg.setAll ( 0. );

    if ( _nx != _ny || _nx != _nz || _ny != _nz ) {
      throw aol::Exception ( "qc::AnisotropyVisualizer3d: Enquist Osher scheme input image has to be a cube!", __FILE__, __LINE__ );
    }

    if ( _verbose ) cerr << aol::color::reset << "Generating a sphere as initial image ...";
    for ( int i = 0; i < _nx; i++ ) {
      for ( int j = 0; j < _ny; j++ ) {
        for ( int k = 0; k < _nz; k++ ) {
          // first: scale to [-1,1]^3
          RealType x = static_cast<RealType> ( i ) / static_cast<RealType> ( _nx );
          RealType y = static_cast<RealType> ( j ) / static_cast<RealType> ( _ny );
          RealType z = static_cast<RealType> ( k ) / static_cast<RealType> ( _nz );

          _outputImg.set ( i,j,k, sqrt( sqr(x-0.5) + sqr(y-0.5) + sqr(z-0.5) ) );
        }
      }
    }
    if ( _verbose ) cerr << "ready. Starting the evolution...\n";

    // create a grid and declare the operator
    int d = qc::logBaseTwo (_nx);
    qc::GridDefinition grid( d, qc::QC_3D );

    GrowWulffshapeEO3D<RealType, AnisoType> EnquistOsherSchemeOp( grid, _anisotropy );
    EnquistOsherSchemeOp.setData( &_outputImg );

    // ------------ get the timestep - size -----------------------------
    tau *= grid.H();

    for ( int iter=0; iter<numSteps; ++iter ) {
      if ( _verbose ) cerr << "\rstep " << aol::color::blue << iter << aol::color::reset;
      EnquistOsherSchemeOp.timeStepEO( tau );
    }

    if ( _verbose ) cerr << "\ndone!\n";
  }


};


}   // end namespace

#endif

