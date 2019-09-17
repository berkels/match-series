#ifndef __NETCDFINTERFACE_H
#define __NETCDFINTERFACE_H

#include <smallVec.h>
#include <vec.h>
#include <scalarArray.h>

#ifdef USE_LIB_NETCDF
#include <netcdf.h>
#include <stack>
#endif

namespace qc {

class NetCDFReader {
  int _ncid, _ndimsp, _dataVarID, _dataNCID;
  aol::Vec3<int> _dim;
public:
  NetCDFReader ( const char *fileName, const std::vector<std::string> &groupNames, const std::string &dataName ) {
#ifdef USE_LIB_NETCDF
    // Declare ids and open data file
    int grpParentId, grpChildId, retval;
    if ((retval = nc_open(fileName, NC_NOWRITE, &_ncid)))
      throw aol::Exception ( nc_strerror(retval), __FILE__, __LINE__ );

    // Recursively enter groups until the leaf group is reached
    grpParentId = grpChildId = _ncid;
    for ( unsigned int i=0; i<groupNames.size ( ) ; ++i ) {
      if ((retval = nc_inq_ncid(grpParentId, groupNames[i].c_str ( ), &grpChildId)))
        throw aol::Exception ( nc_strerror(retval), __FILE__, __LINE__ );
      grpParentId = grpChildId;
    }

    // Find data within leaf group
    if ((retval = nc_inq_varid(grpChildId, dataName.c_str ( ), &_dataVarID)))
      throw aol::Exception ( nc_strerror(retval), __FILE__, __LINE__ );

    // Check number of data dimensions
    if ((retval = nc_inq_varndims(grpChildId, _dataVarID, &_ndimsp)))
      throw aol::Exception ( nc_strerror(retval), __FILE__, __LINE__ );

    // Determine length of data dimensions and reallocate ScalarArray
    aol::Vector<int> dims ( _ndimsp );
    if ((retval =  nc_inq_vardimid(grpChildId, _dataVarID, dims.getData())))
      throw aol::Exception ( nc_strerror(retval), __FILE__, __LINE__ );

    size_t lenp;
    for ( int i = 0; i < dims.size ( ); ++i ) {
      if ((retval = nc_inq_dimlen(grpChildId, dims[i], &lenp)))
        throw aol::Exception ( nc_strerror(retval), __FILE__, __LINE__ );

      _dim[i] = lenp;
    }

    _dataNCID = grpChildId;
#else
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( fileName );
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( groupNames );
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( dataName );
    throw aol::Exception ( "Loading NetCDF data requires the NetCDF library! Compile with -DUSE_NETCDF", __FILE__, __LINE__ );
#endif
  }

  ~NetCDFReader ( ) {
#ifdef USE_LIB_NETCDF
    int retval;
    if ((retval = nc_close(_ncid)))
      cerr << aol::color::error << "Error in ~NetCDFReader: " << nc_strerror(retval) << endl;
#endif
  }

  int getNumDim () const {
    return _ndimsp;
  }

  const aol::Vec3<int> &getDimRef ( ) const {
    return _dim;
  }

  template <typename DataType>
  void readDataTo ( DataType* DataDest ) {
#ifdef USE_LIB_NETCDF
    int retval;
    if ((retval = NetCDFTrait<DataType>::getVar (_dataNCID, _dataVarID, DataDest ) ))
      throw aol::Exception ( nc_strerror(retval), __FILE__, __LINE__ );
#else
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( DataDest );
#endif
  }
};

#ifdef USE_LIB_NETCDF
template <>
class NetCDFTrait<unsigned short> {
public:
  static int getVar ( int ncid, int varid, unsigned short *p ) {
    return nc_get_var_ushort(ncid, varid, p);
  }

  static int defVar ( int ncid, const char* varname, int ndims, int* dimids, int *varid ) {
    return nc_def_var(ncid, varname, NC_USHORT, ndims, dimids, varid);
  }

  static int putVar ( int ncid, int varid, unsigned short *p ) {
    return nc_put_var_ushort(ncid, varid, p);
  }
};

template <>
class NetCDFTrait<float> {
public:
  static int getVar ( int ncid, int varid, float *p ) {
    return nc_get_var_float(ncid, varid, p);
  }

  static int defVar ( int ncid, const char* varname, int ndims, int* dimids, int *varid ) {
    return nc_def_var(ncid, varname, NC_FLOAT, ndims, dimids, varid);
  }

  static int putVar ( int ncid, int varid, float *p ) {
    return nc_put_var_float(ncid, varid, p);
  }
};

template <>
class NetCDFTrait<double> {
public:
  static int getVar ( int ncid, int varid, double *p ) {
    return nc_get_var_double(ncid, varid, p);
  }

  static int defVar ( int ncid, const char* varname, int ndims, int* dimids, int *varid ) {
    return nc_def_var(ncid, varname, NC_DOUBLE, ndims, dimids, varid);
  }

  static int putVar ( int ncid, int varid, double *p ) {
    return nc_put_var_double(ncid, varid, p);
  }
};
#endif

}

#endif  // __NETCDFINTERFACE_H
