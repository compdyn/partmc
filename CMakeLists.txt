CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(PARTMC Fortran)
ADD_SUBDIRECTORY(src)

SET(PACKAGE_BUGREPORT "mwest@illinois.edu")
SET(PACKAGE_NAME "PartMC")
SET(PACKAGE_STRING "PartMC 1.2.0")
SET(PACKAGE_TARNAME "partmc")
SET(PACKAGE_VERSION "1.2.0")

OPTION(ENABLE_MOSAIC "Enable MOSAIC chemistry support" OFF)
OPTION(ENABLE_MPI "Enable MPI support" OFF)

SET(NETCDF_PREFIX "" CACHE PATH "Path to search for NetCDF header and library files")

FIND_PATH(NETCDF_INCLUDE_DIR netcdf.mod
          #{NETCDF_PREFIX}/include
          /usr/include
	  /usr/local/include
          /opt/include
          $ENV{HOME}/opt/include)
FIND_LIBRARY(NETCDF_LIB netcdf
          #{NETCDF_PREFIX}/lib
	  /usr/lib
	  /usr/local/lib
          /opt/lib
          $ENV{HOME}/opt/lib)
