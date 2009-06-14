cmake_minimum_required(VERSION 2.6)
project(PARTMC Fortran)

set(PACKAGE_BUGREPORT "mwest@illinois.edu")
set(PACKAGE_NAME "PartMC")
set(PACKAGE_STRING "PartMC 1.2.0")
set(PACKAGE_TARNAME "partmc")
set(PACKAGE_VERSION "1.2.0")

option(ENABLE_MOSAIC "Enable MOSAIC chemistry support" OFF)

######################################################################

find_path(NETCDF_INCLUDE_DIR netcdf.mod
  DOC "NetCDF include directory (must contain netcdf.mod)"
  PATHS
  /include
  /usr/include
  /usr/local/include
  /opt/include
  /usr/lib64/gfortran/modules
  /usr/lib/gfortran/modules
  $ENV{HOME}/include
  $ENV{HOME}/opt/include
  PATH_SUFFIXES netcdf-3)
find_library(NETCDF_C_LIB netcdf
  DOC "NetCDF C library"
  PATHS
  /usr/lib
  /usr/lib64
  /usr/local/lib
  /usr/local/lib64
  $ENV{HOME}/opt/lib
  $ENV{HOME}/opt/lib64)
find_library(NETCDF_FORTRAN_LIB netcdff
  DOC "NetCDF Fortran library"
  PATHS
  /usr/lib
  /usr/lib64
  /usr/local/lib
  /usr/local/lib64
  $ENV{HOME}/opt/lib
  $ENV{HOME}/opt/lib64)

set(NETCDF_LIBS ${NETCDF_C_LIB})
if(NETCDF_FORTRAN_LIB)
  set(NETCDF_LIBS ${NETCDF_LIBS} ${NETCDF_FORTRAN_LIB})
endif(NETCDF_FORTRAN_LIB)

######################################################################

find_path(MOSAIC_INCLUDE_DIR module_data_mosaic_main.mod
  DOC "MOSAIC include directory"
  PATHS
  /include
  /usr/include
  /usr/local/include
  /opt/include
  $ENV{HOME}/include
  $ENV{HOME}/opt/include
  $ENV{HOME}/subversion/mosaic/trunk/datamodules
  $ENV{HOME}/svn/mosaic/trunk/datamodules
  ${PROJECT_SOURCE_DIR}/../mosaic/datamodules
  ${PROJECT_SOURCE_DIR}/../../mosaic/trunk/datamodules)
find_library(MOSAIC_LIB mosaic
  DOC "MOSAIC library"
  PATHS
  $ENV{HOME}/opt/lib
  $ENV{HOME}/subversion/mosaic/trunk
  $ENV{HOME}/svn/mosaic/trunk
  ${PROJECT_SOURCE_DIR}/../mosaic
  ${PROJECT_SOURCE_DIR}/../../mosaic/trunk)

######################################################################

set(EXTRA_LIBS)

include_directories(${NETCDF_INCLUDE_DIR})
set(EXTRA_LIBS ${EXTRA_LIBS} ${NETCDF_LIBS})

if(ENABLE_MOSAIC)
  include_directories(${MOSAIC_INCLUDE_DIR})
  set(EXTRA_LIBS ${EXTRA_LIBS} ${MOSAIC_LIB})
  add_definitions(-DPMC_USE_MOSAIC)
endif(ENABLE_MOSAIC)

get_filename_component(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# MW 2009-01-20: Below we add -I/usr/include manually because gfortran
# does not search this by default [1], but cmake strips it from paths
# given by include_directories(), assuming that the compiler will
# search it anyway (as gcc does). I reported this as a cmake bug [2].
# [1] http://gcc.gnu.org/bugzilla/show_bug.cgi?id=35707
# [2] http://www.cmake.org/Bug/view.php?id=8408

#if(COMPILER_NAME STREQUAL "gfortran")
#if(CMAKE_COMPILER_IS_GNUG77)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_definitions(-g -x f95-cpp-input -fimplicit-none -W -Wall -Wconversion -Wunderflow -Wimplicit-interface -Wno-unused -fbounds-check -I/include -I/usr/include -I/usr/local/include)
elseif(COMPILER_NAME STREQUAL "pgf95")
  add_definitions(-Mpreprocess)
endif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

######################################################################

enable_testing()
add_custom_target(copy_test ALL ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/test ${CMAKE_BINARY_DIR}/run_test)
set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES run_test)

add_test(bidisperse_1 ${CMAKE_BINARY_DIR}/run_test/bidisperse/run_test_1.sh)
add_test(bidisperse_2 ${CMAKE_BINARY_DIR}/run_test/bidisperse/run_test_2.sh)
add_test(bidisperse_3 ${CMAKE_BINARY_DIR}/run_test/bidisperse/run_test_3.sh)
add_test(brownian_1 ${CMAKE_BINARY_DIR}/run_test/brownian/run_test_1.sh)
add_test(brownian_2 ${CMAKE_BINARY_DIR}/run_test/brownian/run_test_2.sh)
add_test(brownian_3 ${CMAKE_BINARY_DIR}/run_test/brownian/run_test_3.sh)
add_test(brownian_4 ${CMAKE_BINARY_DIR}/run_test/brownian/run_test_4.sh)
add_test(brownian_5 ${CMAKE_BINARY_DIR}/run_test/brownian/run_test_5.sh)
add_test(brownian_6 ${CMAKE_BINARY_DIR}/run_test/brownian/run_test_6.sh)
add_test(emission_1 ${CMAKE_BINARY_DIR}/run_test/emission/run_test_1.sh)
add_test(emission_2 ${CMAKE_BINARY_DIR}/run_test/emission/run_test_2.sh)
add_test(emission_3 ${CMAKE_BINARY_DIR}/run_test/emission/run_test_3.sh)
add_test(emission_4 ${CMAKE_BINARY_DIR}/run_test/emission/run_test_4.sh)
add_test(emission_5 ${CMAKE_BINARY_DIR}/run_test/emission/run_test_5.sh)
add_test(golovin_1 ${CMAKE_BINARY_DIR}/run_test/golovin/run_test_1.sh)
add_test(golovin_2 ${CMAKE_BINARY_DIR}/run_test/golovin/run_test_2.sh)
if(ENABLE_MOSAIC)
  add_test(mosaic_1 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_1.sh)
  add_test(mosaic_2 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_2.sh)
  add_test(mosaic_3 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_3.sh)
  add_test(mosaic_4 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_4.sh)
  add_test(mosaic_5 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_5.sh)
  add_test(mosaic_6 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_6.sh)
  add_test(mosaic_7 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_7.sh)
  add_test(mosaic_8 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_8.sh)
endif(ENABLE_MOSAIC)
add_test(poisson_1 ${CMAKE_BINARY_DIR}/run_test/poisson/run_test_1.sh)
add_test(poisson_2 ${CMAKE_BINARY_DIR}/run_test/poisson/run_test_2.sh)
add_test(poisson_3 ${CMAKE_BINARY_DIR}/run_test/poisson/run_test_3.sh)
add_test(poisson_4 ${CMAKE_BINARY_DIR}/run_test/poisson/run_test_4.sh)
add_test(poisson_5 ${CMAKE_BINARY_DIR}/run_test/poisson/run_test_5.sh)
add_test(sedi_1 ${CMAKE_BINARY_DIR}/run_test/sedi/run_test_1.sh)
add_test(sedi_2 ${CMAKE_BINARY_DIR}/run_test/sedi/run_test_2.sh)

######################################################################

add_executable(partmc src/partmc.f90 src/aero_state.f90
               src/aero_binned.f90 src/bin_grid.f90
               src/condensation.f90 src/constants.f90 src/env_data.f90
               src/env_state.f90 src/aero_mode.f90 src/aero_dist.f90
               src/kernel_golovin.f90 src/kernel_sedi.f90
               src/kernel_constant.f90 src/kernel_brown.f90
               src/kernel_zero.f90 src/aero_data.f90 src/run_exact.f90
               src/run_part.f90 src/util.f90 src/run_sect.f90
               src/output.f90 src/mosaic.f90 src/gas_data.f90
               src/gas_state.f90 src/coagulation.f90 src/kernel.f90
               src/spec_line.f90 src/spec_file.f90 src/rand.f90
               src/aero_particle.f90 src/aero_particle_array.f90
               src/mpi.f90 src/netcdf.f90 src/aero_info.f90
               src/aero_info_array.f90)
target_link_libraries(partmc ${EXTRA_LIBS})

######################################################################

add_executable(test_bidisperse_ode
	       test/bidisperse/test_bidisperse_ode.f90
	       src/kernel_sedi.f90 src/env_data.f90 src/env_state.f90
	       src/constants.f90 src/aero_data.f90 src/util.f90
	       src/gas_data.f90 src/gas_state.f90 src/aero_state.f90
	       src/output.f90 src/bin_grid.f90 src/spec_line.f90
	       src/spec_file.f90 src/aero_mode.f90 src/aero_dist.f90
	       src/aero_binned.f90 src/rand.f90 src/aero_particle.f90
	       src/aero_particle_array.f90 src/mpi.f90 src/netcdf.f90
	       src/aero_info.f90 src/aero_info_array.f90)
target_link_libraries(test_bidisperse_ode ${EXTRA_LIBS})

add_executable(test_bidisperse_extract
               test/bidisperse/test_bidisperse_extract.f90)
target_link_libraries(test_bidisperse_extract ${NETCDF_LIBS})

add_executable(test_poisson_sample test/poisson/test_poisson_sample.f90
               src/util.f90 src/rand.f90 src/constants.f90)
target_link_libraries(test_poisson_sample)

######################################################################

add_executable(extract_aero_particle_mass
	src/extract_aero_particle_mass.f90)
target_link_libraries(extract_aero_particle_mass ${NETCDF_LIBS})

add_executable(extract_aero_size_num
	src/extract_aero_size_num.f90)
target_link_libraries(extract_aero_size_num ${NETCDF_LIBS})

add_executable(extract_aero_size_mass
	src/extract_aero_size_mass.f90)
target_link_libraries(extract_aero_size_mass ${NETCDF_LIBS})

add_executable(extract_aero_species
	src/extract_aero_species.f90)
target_link_libraries(extract_aero_species ${NETCDF_LIBS})

add_executable(extract_aero_total
	src/extract_aero_total.f90)
target_link_libraries(extract_aero_total ${NETCDF_LIBS})

add_executable(extract_gas src/extract_gas.f90)
target_link_libraries(extract_gas ${NETCDF_LIBS})

add_executable(extract_env src/extract_env.f90)
target_link_libraries(extract_env ${NETCDF_LIBS})

######################################################################

add_executable(extract_sectional_aero_size_num
	src/extract_sectional_aero_size_num.f90)
target_link_libraries(extract_sectional_aero_size_num ${NETCDF_LIBS})

add_executable(extract_sectional_aero_size_mass
	src/extract_sectional_aero_size_mass.f90)
target_link_libraries(extract_sectional_aero_size_mass ${NETCDF_LIBS})

add_executable(extract_sectional_aero_species
	src/extract_sectional_aero_species.f90)
target_link_libraries(extract_sectional_aero_species ${NETCDF_LIBS})

add_executable(extract_sectional_aero_total
	src/extract_sectional_aero_total.f90)
target_link_libraries(extract_sectional_aero_total ${NETCDF_LIBS})

######################################################################

add_executable(numeric_diff src/numeric_diff.f90)

add_executable(numeric_average src/numeric_average.f90)

######################################################################
