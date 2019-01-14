#ifndef __SIMPLEXLOOKUP_H
#define __SIMPLEXLOOKUP_H

#include <quoc.h>
#include <smallVec.h>

namespace qc {

namespace simplex {

// --------------------------------------------------------------------------

template <Dimension Dim>
struct TopologyLookup {};

template <typename RealType, Dimension Dim>
struct Lookup {};

// --------------------------------------------------------------------------

template <>
struct TopologyLookup<QC_2D> {
  //! localIndicesSimplexToCube[i][j] is the local
  //! index wrt square element of the j-th vertex
  //! in the i-th triangle.
  static const short localIndicesSimplexToCube[2][3];
  static const short offsetLocalIndices[4][2];
  static const short edges[3][2];
  static const short edgesAsBits[3];

  static const short numSimplexesPerCube = 2;
  static const short numEdges = 3;
  static const short maxEdgeBit = 6;
};

// --------------------------------------------------------------------------

template <typename RealType>
struct Lookup<RealType, QC_2D> :
  public TopologyLookup<QC_2D> {

  typedef aol::Vec2<RealType> VecType;
  //! gradients[i][j] is the cartesian gradient
  //! of the j-th basis function on the i-th triangle.
  static const VecType gradients[2][3];

  static const VecType cartesianCubeNodeCoords[4];
};

// --------------------------------------------------------------------------

template <>
struct TopologyLookup<QC_3D> {
  //! localIndicesSimplexToCube[i][j] is the local
  //! index wrt cube element of the j-th vertex
  //! in the i-th tetrahedron.
  static const short localIndicesSimplexToCube[6][4];
  static const short offsetLocalIndices[8][3];
  static const short edges[6][2];
  static const short edgesAsBits[6];

  static const short numSimplexesPerCube = 6;
  static const short numEdges = 6;
  static const short maxEdgeBit = 12;
};

// --------------------------------------------------------------------------

template <typename RealType>
struct Lookup<RealType, QC_3D> :
  public TopologyLookup<QC_3D> {

  typedef aol::Vec3<RealType> VecType;
  //! gradients[i][j] is the cartesian gradient
  //! of the j-th basis function on the i-th tetrahedron.
  static const VecType gradients[6][4];

  static const VecType cartesianCubeNodeCoords[8];
};

// --------------------------------------------------------------------------

} // end of namespace simplex.

} // end of namespace qc.

#endif
