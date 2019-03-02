find_path(NETCDF_INCLUDE_DIR netcdf.mod NETCDF.mod
  DOC "NetCDF include directory (must contain netcdf.mod)"
  PATHS
  $ENV{NETCDF_HOME}/include
  /usr/lib/gfortran/modules
  /usr/lib64/gfortran/modules
  /opt/local/include)
find_library(NETCDF_C_LIB netcdf
  DOC "NetCDF C library"
  PATHS $ENV{NETCDF_HOME}/lib /opt/local/lib)
find_library(NETCDF_FORTRAN_LIB netcdff
  DOC "NetCDF Fortran library"
  PATHS $ENV{NETCDF_HOME}/lib /opt/local/lib)
if(NETCDF_FORTRAN_LIB)
  set(NETCDF_LIBS ${NETCDF_FORTRAN_LIB})
endif()
set(NETCDF_LIBS ${NETCDF_LIBS} ${NETCDF_C_LIB})
include_directories(${NETCDF_INCLUDE_DIR})
