#ifndef __FINITEDIFFERENCES_H
#define __FINITEDIFFERENCES_H

#include <multiArray.h>
#include <ctrlCCatcher.h>
#include <ChanVese.h>
#include <quocTimestepSaver.h>
#include <clustering.h>
#include <eigenvectors.h>
#include <QRDecomposition.h>

namespace qc {
  
/**
 * This struct + static function construction is a workaround to C++ limitations.
 *
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
struct doCalculateBackwardFDDivergence {
  static void apply ( const qc::MultiArray<RealType, Dim> &MArg, qc::ScalarArray<RealType, Dim> &Divergence );
};

template <typename RealType>
struct doCalculateBackwardFDDivergence<RealType, qc::QC_2D> {
  static void apply ( const qc::MultiArray<RealType, qc::QC_2D> &MArg, qc::ScalarArray<RealType, qc::QC_2D> &Divergence ) {
    const int numX = Divergence.getNumX();
    const int numY = Divergence.getNumY();
    // inner nodes
    for ( int y = 1; y < numY - 1; ++y ) {
      for ( int x = 1; x < numX - 1; ++x ) {
        Divergence.set ( x, y, MArg[0].get(x,y) -MArg[0].get(x-1,y) + MArg[1].get(x,y) -MArg[1].get(x,y-1) );
      }
    }

    // top and bottom nodes (without corners)
    for ( int y = 1; y < numY - 1; ++y ) {
      Divergence.set ( 0, y, MArg[0].get(0,y) + MArg[1].get(0,y)  -MArg[1].get(0,y-1) );
      Divergence.set ( numX-1, y, -MArg[0].get(numX-2,y) + MArg[1].get(numX-1,y) -MArg[1].get(numX-1,y-1) );
    }

    // left and right nodes (without corners)
    for ( int x = 1; x < numX-1; ++x ) {
      Divergence.set ( x, 0, MArg[0].get(x,0) -MArg[0].get(x-1,0) + MArg[1].get(x,0) );
      Divergence.set ( x, numY-1, MArg[0].get(x,numY-1) -MArg[0].get(x-1,numY-1) -MArg[1].get(x,numY-2) );
    }

    // corner (x,y) = (0,0)
    Divergence.set ( 0, 0, MArg[0].get(0,0) + MArg[1].get(0,0) );

    // corner (x,y) = (0,numY-1)
    Divergence.set ( 0, numY-1, MArg[0].get(0,numY-1) - MArg[1].get(0,numY-2) );

    // corner (x,y) = (numX-1,0)
    Divergence.set ( numX-1, 0, -MArg[0].get(numX-2,0) + MArg[1].get(numX-1,0) );

    // corner (x,y) = (numX-1,numY-1)
    Divergence.set ( numX-1, numY-1, -MArg[0].get(numX-2,numY-1) - MArg[1].get(numX-1,numY-2) );

    // Straightforward and readable, but less efficient implementation:
    /*
    Divergence.setZero();
    for ( int x = 0; x < numX; ++x ) {
      for ( int y = 0; y < numY; ++y ) {
        if ( x < numX - 1 )
          Divergence.add ( x, y, MArg[0].get(x,y) );
        if ( x > 0 )
          Divergence.add ( x, y, -MArg[0].get(x-1,y) );
        if ( y < numY - 1 )
          Divergence.add ( x, y, MArg[1].get(x,y) );
        if ( y > 0 )
          Divergence.add ( x, y, -MArg[1].get(x,y-1) );
      }
    }
    */
  }
};

//! \todo Rewrite to get rid of all the if statements inside the for loops.
template <typename RealType>
struct doCalculateBackwardFDDivergence<RealType, qc::QC_3D> {
  static void apply ( const qc::MultiArray<RealType, qc::QC_3D> &MArg, qc::ScalarArray<RealType, qc::QC_3D> &Divergence ) {
    const int numX = Divergence.getNumX();
    const int numY = Divergence.getNumY();
    const int numZ = Divergence.getNumZ();
    // Straightforward and readable, but not most efficient implementation:
    Divergence.setZero();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int z = 0; z < numZ; ++z ) {
      for ( int y = 0; y < numY; ++y ) {
        for ( int x = 0; x < numX; ++x ) {
          if ( x < numX - 1 )
            Divergence.add ( x, y, z, MArg[0].get(x,y,z) );
          if ( x > 0 )
            Divergence.add ( x, y, z, -MArg[0].get(x-1,y,z) );
          if ( y < numY - 1 )
            Divergence.add ( x, y, z, MArg[1].get(x,y,z) );
          if ( y > 0 )
            Divergence.add ( x, y, z, -MArg[1].get(x,y-1,z) );
          if ( z < numZ - 1 )
            Divergence.add ( x, y, z, MArg[2].get(x,y,z) );
          if ( z > 0 )
            Divergence.add ( x, y, z, -MArg[2].get(x,y,z-1) );
        }
      }
    }
  }
};

/**
 * Calulates the divergence with backward finite differences (not scaled by the grid width),
 * as used in {An Algorithm for Total Variation Minimization and Applications} by Antonin Chambolle.
 *
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
void calculateBackwardFDDivergence ( const qc::MultiArray<RealType, Dim> &MArg, qc::ScalarArray<RealType, Dim> &Divergence ) {
  doCalculateBackwardFDDivergence<RealType, Dim>::apply ( MArg, Divergence );
}

template <typename RealType, qc::Dimension Dim>
struct doCalculateCentralFDDivergencePBC {
  static void apply ( const qc::MultiArray<RealType, Dim> &MArg, qc::ScalarArray<RealType, Dim> &Divergence );
};

//! \todo Rewrite to get rid of all the if statements inside the for loops.
template <typename RealType>
struct doCalculateCentralFDDivergencePBC<RealType, qc::QC_3D> {
  static void apply ( const qc::MultiArray<RealType, qc::QC_3D> &MArg, qc::ScalarArray<RealType, qc::QC_3D> &Divergence ) {
    const int numX = Divergence.getNumX();
    const int numY = Divergence.getNumY();
    const int numZ = Divergence.getNumZ();
    // Straightforward and readable, but not most efficient implementation:
    Divergence.setZero();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int z = 0; z < numZ; ++z ) {
      for ( int y = 0; y < numY; ++y ) {
        for ( int x = 0; x < numX; ++x ) {
          Divergence.add ( x, y, z, MArg[0].get((x+1) %numX,y,z) );
          Divergence.add ( x, y, z, -MArg[0].get((x-1 + numX)%numX,y,z) );
          Divergence.add ( x, y, z, MArg[1].get(x,(y+1) %numY,z) );
          Divergence.add ( x, y, z, -MArg[1].get(x,(y-1 + numY) % numY,z) );
          Divergence.add ( x, y, z, MArg[2].get(x,y,(z+1) %numZ) );
          Divergence.add ( x, y, z, -MArg[2].get(x,y,(z-1 + numZ) % numZ) );
          Divergence.set ( x, y, z, Divergence.get(x,y,z) * 0.5 );
        }
      }
    }
  }
};

/**
 * Calulates the divergence with central finite differences (not scaled by the grid width),
 * and periodic boundary conditions
 *
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
void calculateCentralFDDivergencePBC ( const qc::MultiArray<RealType, Dim> &MArg, qc::ScalarArray<RealType, Dim> &Divergence ) {
  doCalculateCentralFDDivergencePBC<RealType, Dim>::apply ( MArg, Divergence );
}

/**
 * This struct + static function construction is a workaround to C++ limitations.
 *
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
struct doCalculateForwardFDGradient {
  static void apply ( const qc::ScalarArray<RealType, qc::QC_2D> &Arg, qc::MultiArray<RealType, qc::QC_2D> &Gradient );
};

template <typename RealType>
struct doCalculateForwardFDGradient<RealType, qc::QC_2D> {
  static void apply ( const qc::ScalarArray<RealType, qc::QC_2D> &Arg, qc::MultiArray<RealType, qc::QC_2D> &Gradient ) {
    const int numX = Arg.getNumX();
    const int numY = Arg.getNumY();
    for ( int y = 0; y < numY-1; ++y ) {
      for ( int x = 0; x < numX-1; ++x ) {
        Gradient[0].set ( x, y, Arg.get(x+1,y) - Arg.get(x,y) );
        Gradient[1].set ( x, y, Arg.get(x,y+1) - Arg.get(x,y) );
      }
    }

    for ( int x = 0; x < numX-1; ++x ) {
      Gradient[0].set ( x, numY-1, Arg.get(x+1,numY-1) - Arg.get(x, numY-1) );
      Gradient[1].set ( x, numY-1, 0 );
    }

    for ( int y = 0; y < numY-1; ++y ) {
      Gradient[0].set ( numX-1, y, 0 );
      Gradient[1].set ( numX-1, y, Arg.get(numX-1,y+1) - Arg.get(numX-1,y) );
    }

    Gradient[0].set ( numX-1, numY-1, 0 );
    Gradient[1].set ( numX-1, numY-1, 0 );

    // Straightforward and readable, but less efficient implementation:
    /*
    for ( int x = 0; x < numX; ++x ) {
      for ( int y = 0; y < numY; ++y ) {
        if ( x < numX - 1 )
          Gradient[0].set ( x, y, Arg.get(x+1,y) - Arg.get(x,y) );
        else
          Gradient[0].set ( x, y, 0 );
        if ( y < numY - 1 )
          Gradient[1].set ( x, y, Arg.get(x,y+1) - Arg.get(x,y) );
        else
          Gradient[1].set ( x, y, 0 );
      }
    }
    */
  }
};

//! \todo Rewrite to get rid of all the if statements inside the for loops.
template <typename RealType>
struct doCalculateForwardFDGradient<RealType, qc::QC_3D> {
  static void apply ( const qc::ScalarArray<RealType, qc::QC_3D> &Arg, qc::MultiArray<RealType, qc::QC_3D> &Gradient ) {
    const int numX = Arg.getNumX();
    const int numY = Arg.getNumY();
    const int numZ = Arg.getNumZ();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    // Straightforward and readable, but not most efficient implementation:
    for ( int z = 0; z < numZ; ++z ) {
      for ( int y = 0; y < numY; ++y ) {
        for ( int x = 0; x < numX; ++x ) {
          if ( x < numX - 1 )
            Gradient[0].set ( x, y, z, Arg.get(x+1,y,z) - Arg.get(x,y,z) );
          else
            Gradient[0].set ( x, y, z, 0 );
          if ( y < numY - 1 )
            Gradient[1].set ( x, y, z, Arg.get(x,y+1,z) - Arg.get(x,y,z) );
          else
            Gradient[1].set ( x, y, z, 0 );
          if ( z < numZ - 1 )
            Gradient[2].set ( x, y, z, Arg.get(x,y,z+1) - Arg.get(x,y,z) );
          else
            Gradient[2].set ( x, y, z, 0 );
        }
      }
    }
  }
};

/**
 * Calulates the gradient with forward finite differences (not scaled by the grid width),
 * as used in {An Algorithm for Total Variation Minimization and Applications} by Antonin Chambolle.
 *
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
void calculateForwardFDGradient ( const qc::ScalarArray<RealType, Dim> &Arg, qc::MultiArray<RealType, Dim> &Gradient ) {
  doCalculateForwardFDGradient<RealType, Dim>::apply ( Arg, Gradient );
}

template <typename RealType, qc::Dimension Dim>
struct doCalculateCentralFDGradientPBC{
  static void apply ( const qc::ScalarArray<RealType, qc::QC_2D> &Arg, qc::MultiArray<RealType, qc::QC_2D> &Gradient );
};

//! \todo Rewrite to get rid of all the if statements inside the for loops.
template <typename RealType>
struct doCalculateCentralFDGradientPBC<RealType, qc::QC_3D> {
  static void apply ( const qc::ScalarArray<RealType, qc::QC_3D> &Arg, qc::MultiArray<RealType, qc::QC_3D> &Gradient ) {
    const int numX = Arg.getNumX();
    const int numY = Arg.getNumY();
    const int numZ = Arg.getNumZ();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    // Straightforward and readable, but not most efficient implementation:
    for ( int z = 0; z < numZ; ++z ) {
      for ( int y = 0; y < numY; ++y ) {
        for ( int x = 0; x < numX; ++x ) {
            Gradient[0].set ( x, y, z, 0.5 * (Arg.get((x+1) % numX,y,z) - Arg.get((x-1 + numX) % numX,y,z) ));
            Gradient[1].set ( x, y, z, 0.5 * (Arg.get(x,(y+1) %numY,z) - Arg.get(x,(y-1 + numY) % numY,z) ));
            Gradient[2].set ( x, y, z, 0.5 * (Arg.get(x,y,(z+1) % numZ) - Arg.get(x,y,(z-1 + numZ) % numZ) ));
        }
      }
    }
  }
};

/**
 * Calulates the gradient with central finite differences (not scaled by the grid width),
 * and periodic boundary conditions.
 *
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
void calculateCentralFDGradientPBC ( const qc::ScalarArray<RealType, Dim> &Arg, qc::MultiArray<RealType, Dim> &Gradient ) {
  doCalculateCentralFDGradientPBC<RealType, Dim>::apply ( Arg, Gradient );
}

/**
 * Collects very basic stuff for iterative algorithms.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class TVAlgorithmBase {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const RealType _gamma;
  int _maxIterations;
  RealType _stopEpsilon;
  const DefaultArraySaver<RealType, ConfiguratorType::Dim>* _pStepSaver;
  bool _quietMode;
public:

  TVAlgorithmBase ( const typename ConfiguratorType::InitType &Initializer,
                    const RealType Gamma,
                    const int MaxIterations,
                    const RealType StopEpsilon )
    : _grid ( Initializer ),
      _gamma ( Gamma ),
      _maxIterations ( MaxIterations ),
      _stopEpsilon ( StopEpsilon ),
      _pStepSaver ( NULL ),
      _quietMode ( false ) {}

  void setMaxIterations ( const int MaxIterations ) {
    _maxIterations = MaxIterations;
  }

  int getMaxIterations ( ) const {
    return _maxIterations;
  }

  void setStopEpsilon ( const RealType StopEpsilon ) {
    _stopEpsilon = StopEpsilon;
  }

  RealType getStopEpsilon ( ) const {
    return _stopEpsilon;
  }

  void setStepSaverReference ( const DefaultArraySaver<RealType, ConfiguratorType::Dim> &StepSaver ) {
    _pStepSaver = &StepSaver;
  }

  const DefaultArraySaver<RealType, ConfiguratorType::Dim> *getStepSaverPointer ( ) const {
    return _pStepSaver;
  }
  
  void setQuietMode ( bool qmode ) {
    _quietMode = qmode;
  }
};

} // namespace qc

#endif // __FINITEDIFFERENCES_H
