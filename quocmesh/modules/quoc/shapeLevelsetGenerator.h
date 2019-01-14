#ifndef __SHAPELEVELSETGENERATOR_H
#define __SHAPELEVELSETGENERATOR_H

#include <scalarArray.h>

namespace qc {

/** Class to generate voxel datasets whose zero level set represents specific geometric shapes.
 *  The values away from the zero level set are not guaranteed to satisfy specific properties (such as being a signed distance function).
 *  The class only contains static methods and should not be instantiated (there is no need to do so).
 *  \see DataGenerator
 *  \todo improve constness, move to qc module
 *  \author Schwen
 */
template< typename DataType >
class ShapeLevelsetGenerator {
private:
  ShapeLevelsetGenerator ( );
  ShapeLevelsetGenerator ( const ShapeLevelsetGenerator &other );
  ShapeLevelsetGenerator& operator= ( const ShapeLevelsetGenerator &other );

public:

  //! a single column as computational domain
  static void generateColumnLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const DataType radius = 0.3 );

  //! a single elliptic column as computational domain
  static void generateEllipticColumnLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const DataType x_rad = 0.35, const DataType y_rad = 0.2 );

  //! a single cone segment as computational domain
  static void generateConeLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const DataType low_rad = 0.45, const DataType high_rad = 0.15 );

  //! a single ball as computational domain, centered in the domain
  static void generateBallLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const DataType radius = 0.6 );

  //! a single ball as computational domain, centered in the domain
  static void generateEllipsoidLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const aol::Vec3<DataType> radii, const aol::Vec3<DataType> center = aol::Vec3<DataType> ( 0.5, 0.5, 0.5 ) );

  //! plates of thickness 0.1 at bottom and top, n_rods squared rods of radius rad in between as computational domain
  static void generatePlateRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods, const DataType rad );

  //! plate-rods levelset with radius determined by number of rods
  static void generatePlateRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods );

  //! plates of thickness 0.1 on all faces and a 3D array of rods in the interior
  static void generate3DPlaterodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods, const DataType rad );

  //! 3d plate-rods levelset with radii determined by number of rods
  static void generate3DPlaterodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods );

  //! a 3D array of rods, no plates on the outside
  static void generate3DRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods, const DataType rad );

  //! 3d rods levelset with radii determined by number of rods
  static void generate3DRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods );

  //! a 3D array of rods, no plates on the outside, specified rods removed
  static void generateAnisoSelected3DRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods,
                                                    const DataType x_rad, const DataType y_rad, const DataType z_rad,
                                                    const qc::BitArray<qc::QC_3D> &selMask_x, const qc::BitArray<qc::QC_3D> &selMask_y, const qc::BitArray<qc::QC_3D> &selMask_z );

  //! a 3D array of rods, no plates on the outside, specified rods removed
  static void generateSelected3DRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods, const DataType rad,
                                               const qc::BitArray<qc::QC_3D> &selMask_x, const qc::BitArray<qc::QC_3D> &selMask_y, const qc::BitArray<qc::QC_3D> &selMask_z );

  //! a 3D array of rods, no plates on the outside, with approximately the specified fraction of rods removed
  static void generateRandomRodsRemoved3DRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods, const DataType rad, const DataType removeFraction );

  //! a 3D array of rods, no plates on the outside, with approximately the specified fraction of rods removed
  static void generateAnisoRandomRodsRemoved3DRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods,
                                                             const DataType x_rad, const DataType y_rad, const DataType z_rad,
                                                             const DataType x_removeFraction, const DataType y_removeFraction, const DataType z_removeFraction,
                                                             const unsigned int seed = 0 );

  static void generatePeriodicAnisoRandom3DRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods,
                                                          const DataType x_rad, const DataType y_rad, const DataType z_rad,
                                                          const DataType x_removeFraction, const DataType y_removeFraction, const DataType z_removeFraction,
                                                          const unsigned int seed = 0 );

  //! Whole unit cube as computational domain
  static void generateSolidcubeLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset );

  //! brick-shaped part of the unit cube as computational domain; [0, 0, 0] x [numX, threshold * numY, numZ] or in other direction
  static void generateHalfcubeLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const DataType threshold = 0.50001, const qc::Comp direction = qc::QC_Y );

  static void generateRotatedHalfcubeLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const DataType threshold = 0.50001, const double angle = 23 * aol::NumberTrait<DataType>::pi / 180 );

  //! Slot levelset with bratwurst effect.
  /* The slot is centered in z-direction, parallel to the xy plane and has depth from the yz plane.
     Hence the following box is missing:
     [0; depth]  x  [0; 1]  x  [ 0.5 - width/2; 0.5 + width/2]
  */
  static void generateSlotLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const DataType width = 0.1, const DataType depth = 0.8 );

  /** @param number_of_layers     number of layers div 2 is the number of laminates (the others are empty space)
      @param alpha_deg rotation   around z axis in degrees
      @param beta_deg rotation    around y axis in degrees
  */
  static void generateLaminateLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int number_of_layers, const DataType alpha_deg = 0.0, const DataType beta_deg = 0.0, const bool periodic = true );

  //! "Egg-box": halfcube with oscillations
  static void generateHalfcubeRippleLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const DataType threshold,
                                               const aol::Vec3<DataType> a, const aol::Vec3<DataType> k );

  //! "double egg-box"
  static void generateTwoPlaneRippleLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset,
                 const aol::Vec3<DataType> a0, const aol::Vec3<DataType> a1, const aol::Vec3<DataType> k0, const aol::Vec3<DataType> k1 );

  //! "quadruple egg box
  static void generateFourPlaneRippleLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset,
           const aol::Vec3<DataType> a0, const aol::Vec3<DataType> a1, const aol::Vec3<DataType> a2, const aol::Vec3<DataType> a3,
           const aol::Vec3<DataType> k0, const aol::Vec3<DataType> k1, const aol::Vec3<DataType> k2, const aol::Vec3<DataType> k3 );

  //! Swiss cheese has spherical holes in the interior
  static void generateSwissCheeseLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_holes, const DataType radius );

  //! Austrian cheese has spherical holes extending out of the bounding box
  static void generateAustrianCheeseLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_holes, const DataType radius );

  static void generateRandomRadiiSwissCheeseLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_holes );

  static void generateSwissCheeseLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_holes );

  static void superimposeBoundingBox ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int numBL = 3, const DataType LSValueThere = -1 );

  //! essentially 2D Z-shaped levelset
  static void generateZLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset );

  //! 2D Z-shaped x const levelset
  static void generateZ3Levelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset );

  //! 2D rotated parallel x const levelset
  static void generateDrunkenZebraLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset );

  //! Generate 2D honeycomb x cylinder levelset, only works for certain values of n (else is non-periodic), e.g. 2, 4, 6, 7, 9, 11, 14, 16, 18, 21, 23, 25, 26, 28, 30; theta is the wall thickness
  static void generateHoneycombLevelset ( qc::ScalarArray<DataType, qc::QC_3D> & levelset, const int n, const DataType theta );

  //! Generate geometrically orthotropic honeycomb levelset, only works for certain values of n.
  static void generateGeomOrthoHoneycombLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n, const aol::Mat<2,2,DataType> &thetaMat );

  //! a 3D array of rods, rotated
  static void generateXRotated3DRodsLevelset ( qc::ScalarArray<DataType, qc::QC_3D> &levelset, const int n_rods, const aol::Vec3<DataType> &radii, const DataType angle );

// end class
};

// end namespace qc
}

#endif
