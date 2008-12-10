cmake_minimum_required(VERSION 2.6)
project(PARTMC Fortran)

set(PACKAGE_BUGREPORT "mwest@illinois.edu")
set(PACKAGE_NAME "PartMC")
set(PACKAGE_STRING "PartMC 1.2.0")
set(PACKAGE_TARNAME "partmc")
set(PACKAGE_VERSION "1.2.0")

option(ENABLE_MOSAIC "Enable MOSAIC chemistry support" OFF)
option(ENABLE_MPI "Enable MPI support" OFF)

find_path(NETCDF_INCLUDE_DIR netcdf.mod
  DOC "NetCDF include directory (must contain netcdf.mod)"
  PATHS
  /usr/include
  /usr/local/include
  /opt/include
  $ENV{HOME}/opt/include)
find_library(NETCDF_LIB netcdf
  DOC "NetCDF library"
  PATHS
  /usr/lib
  /usr/local/lib
  /opt/lib
  $ENV{HOME}/opt/lib)

find_path(MOSAIC_INCLUDE_DIR module_data_mosaic_main.mod
  DOC "MOSAIC include directory"
  PATHS
  /usr/lib
  /usr/local/lib
  /opt/lib
  $ENV{HOME}/opt/lib
  $ENV{HOME}/subversion/mosaic/trunk/datamodules
  $ENV{HOME}/svn/mosaic/trunk/datamodules)
find_library(MOSAIC_LIB mosaic
  DOC "MOSAIC library"
  PATHS
  /usr/lib
  /usr/local/lib
  /opt/lib
  $ENV{HOME}/opt/lib
  $ENV{HOME}/subversion/mosaic/trunk
  $ENV{HOME}/svn/mosaic/trunk)

set(EXTRA_COMPILE_FLAGS "")
set(EXTRA_LINK_FLAGS "")

set(OPTIONAL_LIBS)

if(ENABLE_MOSAIC)
  include_directories(${MOSAIC_INCLUDE_DIR})
  set(OPTIONAL_LIBS ${OPTIONAL_LIBS} ${MOSAIC_LIB})
endif(ENABLE_MOSAIC)

include_directories(${NETCDF_INCLUDE_DIR})

if(CMAKE_Fortran_COMPILER MATCHES "gfortran$")
  set(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -x f95-cpp-input")
elseif(CMAKE_Fortran_COMPILER MATCHES "pgf95$")
  set(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -Mpreprocess")
endif(CMAKE_Fortran_COMPILER MATCHES "gfortran$")

add_subdirectory(src)
