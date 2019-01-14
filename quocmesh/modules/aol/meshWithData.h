#ifndef __MESHWITHDATA_H
#define __MESHWITHDATA_H

#include <qmException.h>
#include <bzipiostream.h>
#include <vec.h>
#include <multiVector.h>
//#include <triangMesh.h>

namespace aol {

//! data can either belong to vertices or to faces
enum DataSupp { VERTEX_DATA, FACE_DATA };
enum PLY_FORMAT { UD_PLY, STANFORD_PLY };

//! vector-valued data can be saved (in VTK legacy format) as 3-vectors,
//! normals or texture coordinates (the file format also supports
//! color scalars and lookup tables, which we will not use).
//! PLY saving does not respect this property.
enum VectorSpec { VECTORS, NORMALS, TCOORDS };


/*  Use this class like
 *
 *  aol::MeshWithData<TriMesh<RealType> > ( mesh )
 *           .addData ( result, "color", VERTEX_DATA )
 *           .saveAsPLY ( "result.ply" );
 *
 */

/** This class is a container for a mesh plus data vectors on vertices and
 *  faces. The mesh informations (geometry and topology, that means
 *  vertex coordinates and connectivity) are stored in a mesh that is given
 *  in the constructor. Its type is a template argument of MeshWithData,
 *  which has to provide an ElementIteratorType and a NodeIteratorType.
 *  Inside the QuocMeshes, aol::TriangMesh and om::TriMesh are meaningful
 *  choices. However, there is normally no need to use this class with
 *  aol::TriangMesh explicitely: The TriangMesh can store data vectors
 *  itself and just uses this class as helper for writing to files.
 *
 *  We will not store any of the data vectors here, but only keep pointers
 *  to them.
 *
 *  Please note that PLY files cannot handle vector-valued data. For this
 *  reason, vectorial data is printed to PLY files as multiple scalar
 *  arrays which are named just like your vector data, only with
 *  "0", "1" etc. appended.
 */
template <class MeshType>
class MeshWithData {
public:
  typedef typename MeshType::RealType RealType;

  MeshWithData ( const MeshType & mesh )
      : _mesh ( mesh ), _precision( 8 )
  {}

  MeshWithData ( const MeshType & mesh, int precision )
      : _mesh ( mesh ), _precision( precision )
  {}

  MeshWithData & addData ( const aol::Vector<RealType> & data, string dataDescr, DataSupp supp ) {
  ScalarData entry = { dataDescr, &data };

    switch ( supp ) {
    case VERTEX_DATA:
      if ( data.size() != _mesh.getNumVertices () ) {
        stringstream msg;
        msg << "MeshWithData::addData(): attemp to add vertex data array with " << data.size ()
            << " entries, but mesh has " << _mesh.getNumVertices () << " vertices.";
        throw aol::DimensionMismatchException ( msg.str(), __FILE__, __LINE__ );
      }
      _scalarVertexData.push_back ( entry );
      break;

    case FACE_DATA:
      if ( data.size() != _mesh.getNumFaces () ) {
        stringstream msg;
        msg << "MeshWithData::addData(): attemp to add face data array with " << data.size ()
            << " entries, but mesh has " << _mesh.getNumFaces () << " vertices.";
        throw aol::DimensionMismatchException ( msg.str(), __FILE__, __LINE__ );
      }
      _scalarFaceData.push_back ( entry );
      break;
      
    default:
      throw aol::UnimplementedCodeException ( "MeshWithData::addData: unknown DataSupp", __FILE__, __LINE__ );
    }
    return *this;
  }

  MeshWithData & addData ( const aol::MultiVector<RealType> & data, string dataDescr, DataSupp supp, VectorSpec vSpec = VECTORS ) {

    if ( !data.allDimsEqual () )
      throw aol::DimensionMismatchException ( "MeshWithData::addData(): attemp to add vector-valued vertex data, "
                                              "but entries of the passed aol::MultiVector have differing sizes.",
                                              __FILE__, __LINE__ );
    VectorData entry = { dataDescr, vSpec, &data };

    switch ( supp ) {
    case VERTEX_DATA:
      if ( data[0].size() != _mesh.getNumVertices () ) {
        stringstream msg;
        msg << "MeshWithData::addData(): attemp to add vertex data array with " << data[0].size ()
            << " entries, but mesh has " << _mesh.getNumVertices () << " vertices.";
        throw aol::DimensionMismatchException ( msg.str(), __FILE__, __LINE__ );
      }
      _vectorVertexData.push_back ( entry );
      break;

    case FACE_DATA:
      if ( data[0].size() != _mesh.getNumFaces () ) {
        stringstream msg;
        msg << "MeshWithData::addData(): attemp to add face data array with " << data[0].size ()
            << " entries, but mesh has " << _mesh.getNumFaces () << " vertices.";
        throw aol::DimensionMismatchException ( msg.str(), __FILE__, __LINE__ );
      }
      _vectorFaceData.push_back ( entry );
      break;
      
    default:
      throw aol::UnimplementedCodeException ( "MeshWithData::addData: unknown DataSupp", __FILE__, __LINE__ );
    }
    return *this;
  }

  //! set precision in saving methods
  void setPrecisionTo( int prec ){ _precision = prec; }

  //! save in Stanford ply format
  void saveAsPLY ( string filename, const bool Binary = false ) const {
    saveAsStanfordOrUDPLY ( filename, aol::STANFORD_PLY, Binary );
  }

  //! save in UD ply format
  void saveAsUDPLY ( string filename ) const {
    saveAsStanfordOrUDPLY ( filename, aol::UD_PLY );
  }  

  void saveAsLegacyVTK ( string filename ) const {

    aol::Bzipofstream out ( filename.c_str() );                         // ctor checks if file could be opened

    out << "# vtk DataFile Version 2.0" << endl
        << "written by method MeshWithData::saveAsLegacyVTK" << endl    // paraview supports time = x.x here. useful?
        << "ASCII" << endl
        << "DATASET POLYDATA" << endl
        << "POINTS " << _mesh.getNumVertices() << " float" << endl;

    // vertex coordinates
    for ( typename MeshType::NodeIteratorType nIter = _mesh; nIter.notAtEnd(); ++nIter ) {
      Vec3<RealType> nCoords = nIter.getCoords();
      for ( short i = 0; i < 3; ++i )
        out << ( i == 0 ? "" : " " ) << nCoords[i];
      out << endl;
    }

    out << "POLYGONS " << _mesh.getNumFaces() << " " << 4 * _mesh.getNumFaces() << endl;
    // triangles' vertex indices
    for ( typename MeshType::ElementIteratorType tIter = _mesh; tIter.notAtEnd(); ++tIter ) {
      Vec3<int> vIndices = tIter.getNodeIndices();
      out << "3 ";
      for ( short i = 0; i < 3; ++i )
        out << ( i == 0 ? "" : " " ) << vIndices[i];
      out << endl;
    }

    if ( _scalarVertexData.size() > 0 || _vectorVertexData.size() > 0 )
      out << "POINT_DATA " << _mesh.getNumVertices() << endl;
    // scalar data on vertices
    for (size_t i = 0; i < _scalarVertexData.size(); ++i) {
      out << "SCALARS " << _scalarVertexData[i]._descr << " float" << endl;
      out << "LOOKUP_TABLE default" << endl;
      for ( int vx = 0; vx < (*_scalarVertexData[i]._data).size(); ++vx )
        out << (*_scalarVertexData[i]._data)[vx] << endl;
    }
    // vector data on vertices
    for (size_t i = 0; i < _vectorVertexData.size(); ++i) {
      string spec;
      if ( _vectorVertexData[i]._spec == VECTORS ) spec = "VECTORS";
      if ( _vectorVertexData[i]._spec == NORMALS ) spec = "NORMALS";
      if ( _vectorVertexData[i]._spec == TCOORDS ) spec = "TEXTURE_COORDINATES";
      out << spec << " " << _vectorVertexData[i]._descr << " float" << endl;
      for ( int vx = 0; vx < (*_vectorVertexData[i]._data)[0].size(); ++vx ) {
        for ( int comp = 0; comp < (*_vectorVertexData[i]._data).numComponents(); ++comp )
          out << ( comp == 0 ? "" : " " ) << (*_vectorVertexData[i]._data)[comp][vx];
        out << endl;
      }
    }

    if ( _scalarFaceData.size() > 0 || _vectorFaceData.size() > 0 )
      out << "CELL_DATA " << _mesh.getNumFaces() << endl;
    // scalar data on faces
    for (size_t i = 0; i < _scalarFaceData.size(); ++i) {
      out << "SCALARS " << _scalarFaceData[i]._descr << " float" << endl;
      out << "LOOKUP_TABLE default" << endl;
      for ( int vx = 0; vx < (*_scalarFaceData[i]._data).size(); ++vx )
        out << (*_scalarFaceData[i]._data)[vx] << endl;
    }
    // vector data on faces
    for (size_t i = 0; i < _vectorFaceData.size(); ++i) {
      string spec;
      if ( _vectorFaceData[i]._spec == VECTORS ) spec = "VECTORS";
      if ( _vectorFaceData[i]._spec == NORMALS ) spec = "NORMALS";
      if ( _vectorFaceData[i]._spec == TCOORDS ) spec = "TEXTURE_COORDINATES";
      out << spec << " " << _vectorFaceData[i]._descr << " float" << endl;
      for ( int tx = 0; tx < (*_vectorFaceData[i]._data)[0].size(); ++tx ) {
        for ( int comp = 0; comp < (*_vectorFaceData[i]._data).numComponents(); ++comp )
          out << ( comp == 0 ? "" : " " ) << (*_vectorFaceData[i]._data)[comp][tx];
        out << endl;
      }
    }
  }

  //! save in obj format (caution: face data is not supported!)
  void saveAsOBJ ( string filename, int precision = 8, bool mode2D = false ) const {    
      
    aol::Bzipofstream out ( filename.c_str() );
    out.precision( precision );

    if ( !out )
      throw ( aol::Exception ( "MeshWithData::saveAsOBJ(): could not open file \"" + filename + "\" for writing", __FILE__, __LINE__ ) );     
    
    out << "# written by method MeshWithData::saveAsOBJ" << endl;
    out << "# " << _mesh.getNumVertices () << " vertices, "<< _mesh.getNumFaces () << " faces" <<endl;
    
    // vertex coordinates
    for ( typename MeshType::NodeIteratorType nIter = _mesh; nIter.notAtEnd(); ++nIter ) {
      Vec3<RealType> nCoords = nIter.getCoords();
      out<<"v "<< nCoords[0] <<" "<<nCoords[1];
      if( !mode2D ) out <<" "<<nCoords[2];
      out<<endl;
    }

    // triangles' vertex indices (obj counting starts from 1!)
    for ( typename MeshType::ElementIteratorType tIter = _mesh; tIter.notAtEnd(); ++tIter ) {
      Vec3<int> vIndices = tIter.getNodeIndices();
      out<<"f "<< vIndices[0]+1 <<" "<<vIndices[1]+1<<" "<<vIndices[2]+1<<endl;
    }
    
    // write vertex data: scalar-valued  
    for ( size_t i = 0; i < _scalarVertexData.size(); ++i ){
      for ( typename MeshType::NodeIteratorType nIter = _mesh; nIter.notAtEnd(); ++nIter ) 
        out << "vp " << (*(_scalarVertexData)[i]._data)[nIter.getIndex()];
      out << endl;
    }
    
    // write texture data: vector-valued  
    for ( size_t i = 0; i < _vectorVertexData.size(); ++i ){
      for ( typename MeshType::NodeIteratorType nIter = _mesh; nIter.notAtEnd(); ++nIter ) {
    	out << "vt";
        for ( int j = 0; j < _vectorVertexData[i]._data->numComponents(); ++j )
          out << " " << (*(_vectorVertexData)[i]._data)[j][nIter.getIndex()];
	out << endl;
      }
    }
    
    if ( _scalarFaceData.size() > 0 || _vectorFaceData.size() > 0 )
      cerr << "WARNING in aol::MeshWithData::saveAsOBJ(): obj format does not support any face data!" << endl;

  }

  // VRML 2.0 data format, used by Blender, JavaView, ...
  void saveAsWRL ( string filename ) const {

    aol::Bzipofstream out ( filename.c_str() );                         // ctor checks if file could be opened

    out << "#VRML V2.0 utf8" << endl
        << "# Produced with QuocMesh: written by method MeshWithData::saveAsWRL" << endl
        << "# File Format = WRL" << endl
        << "DATASET POLYDATA" << endl
        << "#     Number of Vertices = " << _mesh.getNumVertices() << " " << endl
        << "#     Number of Elements = " << _mesh.getNumFaces() << " " << endl
        << "# " << endl
        << "# End of Header" << endl
        << "NavigationInfo {" << endl
        << "   headlight TRUE" << endl
        << "   type      [ \" EXAMINE\", \"WALK\", \"ANY\" ]" << endl
        << "}" << endl
        << "Group {" << endl
        << "children [" << endl
        << "Transform {" << endl
        << "rotation    1.0 0.0 0.0 0.0" << endl
        << "scale       1.0 1.0 1.0" << endl
        << "translation 0.0 0.0 0.0" << endl
        << "children [" << endl
        << "Shape {" << endl
        << "geometry IndexedFaceSet {" << endl
        << "solid FALSE" << endl
        << "coord Coordinate {" << endl
        << "point [" << endl;


    // vertex coordinates
    for ( typename MeshType::NodeIteratorType nIter = _mesh; nIter.notAtEnd(); ++nIter ) {
      Vec3<RealType> nCoords = nIter.getCoords();
      for ( short i = 0; i < 3; ++i )
        out << ( i == 0 ? "" : " " ) << nCoords[i];
      if( nIter.getIndex() == _mesh.getNumVertices() - 1 ){
        out << endl;
      }
      else{
        out << "," << endl;        
      }
    }

    out << "]" << endl
        << "}" << endl
        << "coordIndex [" << endl;

    // triangles' vertex indices
    for ( typename MeshType::ElementIteratorType tIter = _mesh; tIter.notAtEnd(); ++tIter ) {
      Vec3<int> vIndices = tIter.getNodeIndices();
      for ( short i = 0; i < 3; ++i )
        out << ( i == 0 ? "" : " " ) << vIndices[i];
      if( tIter.getIndex() == _mesh.getNumFaces() - 1 ){
        out << " -1" << endl;
      }
      else{
        out << " -1," << endl;
      }
    }

    // vector data on vertices
    for (size_t i = 0; i < _vectorVertexData.size(); ++i) {
      out << "]" << endl
          << "colorPerVertex TRUE" << endl
          << "color Color {" << endl
          << "color [" << endl;

      for ( int vx = 0; vx < (*_vectorVertexData[i]._data)[0].size(); ++vx ) {
        for ( int comp = 0; comp < (*_vectorVertexData[i]._data).numComponents(); ++comp )
          out << ( comp == 0 ? "" : " " ) << (*_vectorVertexData[i]._data)[comp][vx];
        if( vx == _mesh.getNumVertices() - 1 ){
          out << endl;
        }
        else{
          out << "," << endl;        
        }
      }
    }

    // vector data on faces
    for (size_t i = 0; i < _vectorFaceData.size(); ++i) {
      out << "]" << endl
          << "colorPerVertex FALSE" << endl
          << "color Color {" << endl
          << "color [" << endl;

      for ( int tx = 0; tx < (*_vectorFaceData[i]._data)[0].size(); ++tx ) {
        for ( int comp = 0; comp < (*_vectorFaceData[i]._data).numComponents(); ++comp )
          out << ( comp == 0 ? "" : " " ) << (*_vectorFaceData[i]._data)[comp][tx];
        if( tx == _mesh.getNumFaces() - 1 ){
          out << endl;
        }
        else{
          out << "," << endl;        
        }
      }
    }


    out << "]" << endl
        << "}" << endl
        << "}" << endl
        << "appearance Appearance {" << endl
        << "material Material {" << endl
        << "diffuseColor 1.0 1.0 1.0" << endl
        << "}" << endl
        << "}" << endl
        << "}" << endl
        << "]" << endl
        << "}" << endl
        << "]" << endl
        << "}" << endl;
    
  }
  
protected:

  void saveAsStanfordOrUDPLY ( string filename, aol::PLY_FORMAT format, const bool Binary = false ) const {

    Bzipofstream out ( filename.c_str() );
    out.precision( _precision );

    if ( !out )
      throw ( aol::Exception ( "MeshWithData::saveAsStanfordOrUDPLY(): could not open file \""
                               + filename + "\" for writing", __FILE__, __LINE__ ) );

    out << "ply" << endl
        << "format " << ( Binary ? "binary_little_endian" : "ascii" ) << " 1.0" << endl
        << "comment written by method MeshWithData::saveAs" << ( format == UD_PLY ? "UDPLY" : "PLY" ) << endl
        << "element vertex " << _mesh.getNumVertices () << endl
        << "property float x" << endl
        << "property float y" << endl
        << "property float z" << endl;

    // PLY format can only store scalar data.
    // We will write our vector-valued data as scalar-values arrays.
    for ( size_t i = 0; i < _scalarVertexData.size(); ++i )
      out << "property float " << _scalarVertexData[i]._descr << endl;

    for ( size_t i = 0; i < _vectorVertexData.size(); ++i )
      for ( int j = 0; j < _vectorVertexData[i]._data->numComponents(); ++j )
        out << "property float " << _vectorVertexData[i]._descr << j << endl;

    out << "element face " << _mesh.getNumFaces () << endl
        << "property list uchar int vertex_index" << endl;

    for ( size_t i = 0; i < _scalarFaceData.size(); ++i )
      out << "property float " << _scalarFaceData[i]._descr << endl;

    for ( size_t i = 0; i < _vectorFaceData.size(); ++i )
      for ( int j = 0; j < _vectorFaceData[i]._data->numComponents(); ++j )
        out << "property float " << _vectorFaceData[i]._descr << j << endl;

    out << "end_header" << endl;

    // vertex coordinates and data
    for ( typename MeshType::NodeIteratorType nIter = _mesh; nIter.notAtEnd(); ++nIter ) {
      // write vertex coordinates
      aol::Vec3<RealType> nCoords = nIter.getCoords ();
      for ( short i = 0; i < 3; ++i )
        if ( Binary )
          aol::writeBinaryData<RealType, float> ( nCoords[i], out );
        else
          out << ( i == 0 ? "" : " " ) << nCoords[i];
      // write vertex data: (a) scalar-valued
      for ( size_t i = 0; i < _scalarVertexData.size(); ++i )
        if ( Binary )
          aol::writeBinaryData<RealType, float> ( (*_scalarVertexData[i]._data)[nIter.getIndex()], out );
        else
          out << " " << (*_scalarVertexData[i]._data)[nIter.getIndex()];
      // (b) vector-valued
      for ( size_t i = 0; i < _vectorVertexData.size(); ++i )
        for ( int j = 0; j < _vectorVertexData[i]._data->numComponents(); ++j )
          if ( Binary )
            aol::writeBinaryData<RealType, float> ( (*_vectorVertexData[i]._data)[j][nIter.getIndex()], out );
          else
            out << " " << (*_vectorVertexData[i]._data)[j][nIter.getIndex()];
      if ( Binary == false )
        out << endl;
    }

    // triangle vertex indices and data
    for ( typename MeshType::ElementIteratorType tIter = _mesh; tIter.notAtEnd(); ++tIter ) {
      // write node inidices
      aol::Vec3<int> nIndices = tIter.getNodeIndices ();
      if (format == aol::STANFORD_PLY) {
        if ( Binary )
          aol::writeBinaryData<int, unsigned char> ( 3, out );
        else
          out << "3 ";
      }
      for ( short i = 0; i < 3; ++i )
        if ( Binary )
          aol::writeBinaryData<int, int> ( nIndices[i], out );
        else
          out << ( i == 0 ? "" : " " ) << nIndices[i];
      // write face  data: (a) scalar-valued
      for ( size_t i = 0; i < _scalarFaceData.size(); ++i )
        if ( Binary )
          aol::writeBinaryData<RealType, float> ( (*_scalarFaceData[i]._data)[tIter.getIndex()], out );
        else
          out << " " << (*_scalarFaceData[i]._data)[tIter.getIndex()];
      // (b) vector-valued
      for ( size_t i = 0; i < _vectorFaceData.size(); ++i )
        for ( int j = 0; j < _vectorFaceData[i]._data->numComponents(); ++j )
          if ( Binary )
            aol::writeBinaryData<RealType, float> ( (*_vectorFaceData[i]._data)[j][tIter.getIndex()], out );
          else
            out << " " << (*_vectorFaceData[i]._data)[j][tIter.getIndex()];
      if ( Binary == false )
        out << endl;
    }
  }

  typedef aol::Vector<RealType>        ScalarDataType;
  typedef aol::MultiVector<RealType>   VectorDataType;

  struct ScalarData {
    string                 _descr;
    const ScalarDataType * _data;
  };

  struct VectorData {
    string                 _descr;
    VectorSpec             _spec;
    const VectorDataType * _data;
  };

  vector<ScalarData> _scalarVertexData;
  vector<VectorData> _vectorVertexData;
  vector<ScalarData> _scalarFaceData;
  vector<VectorData> _vectorFaceData;

  const MeshType _mesh;
  int _precision;      
};

} // end of namespace aol.

#endif
