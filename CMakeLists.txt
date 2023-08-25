cmake_minimum_required(VERSION 2.6.4)
# Require 2.6.4 due to -I/usr/include behavior:
# http://www.cmake.org/Bug/view.php?id=8598
# http://gcc.gnu.org/bugzilla/show_bug.cgi?id=35707
# http://www.cmake.org/Bug/view.php?id=8408
project(PARTMC Fortran C)

set(PACKAGE_BUGREPORT "mwest@illinois.edu")
set(PACKAGE_NAME "PartMC")
set(PACKAGE_STRING "PartMC 2.7.0")
set(PACKAGE_TARNAME "partmc")
set(PACKAGE_VERSION "2.7.0")

######################################################################
# options

option(ENABLE_GSL "Enable GSL library for random number generation" OFF)
option(ENABLE_MOSAIC "Enable MOSAIC chemistry support" OFF)
option(ENABLE_CAMP "Enable CAMP chemistry support" OFF)
option(ENABLE_MPI "Enable MPI parallel support" OFF)
option(ENABLE_SUNDIALS "Enable SUNDIALS solver for condensation support" OFF)
option(ENABLE_C_SORT "Enable C sorting routines" OFF)

######################################################################
# CPack

set(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${PACKAGE_TARNAME}-${PACKAGE_VERSION}")
set(CPACK_SOURCE_IGNORE_FILES "${CPACK_SOURCE_IGNORE_FILES};/.*~$;/[.].*/;/build/;/figures/;/scenarios/[^123456].*/;/doc/condensation/;/doc/deposition/;/doc/ship_track/;/old/;/tool/;/TODO;/CMakeCache.txt")

include(CPack)

######################################################################
# NetCDF

include(netcdf.cmake)

######################################################################
# GSL

if(ENABLE_GSL)
  find_path(GSL_INCLUDE_DIR gsl/gsl_math.h
    DOC "GSL include directory (must have gsl/ subdir)"
    PATHS $ENV{GSL_HOME}/include /opt/local/include)
  find_library(GSL_LIB gsl
    DOC "GSL library"
    PATHS $ENV{GSL_HOME}/lib /opt/local/lib)
  find_library(GSL_CBLAS_LIB gslcblas
    DOC "GSL CBLAS library"
    PATHS $ENV{GSL_HOME}/lib /opt/local/lib)
  find_library(M_LIB m
    DOC "standard C math library")
  set(GSL_SRC src/rand_gsl.c)
  set(GSL_LIBS ${GSL_LIB} ${GSL_CBLAS_LIB} ${M_LIB})
  include_directories(${GSL_INCLUDE_DIR})
  add_definitions(-DPMC_USE_GSL)
endif()

######################################################################
# C sort

if(ENABLE_C_SORT)
  set(C_SORT_SRC src/sort.c)
  add_definitions(-DPMC_USE_C_SORT)
endif()

######################################################################
# MOSAIC

if(ENABLE_MOSAIC)
  find_path(MOSAIC_INCLUDE_DIR module_data_mosaic_main.mod
    DOC "MOSAIC include directory"
    PATHS $ENV{MOSAIC_HOME}/datamodules $ENV{MOSAIC_HOME}/include)
  find_library(MOSAIC_LIB mosaic
    DOC "MOSAIC library"
    PATHS $ENV{MOSAIC_HOME} $ENV{MOSAIC_HOME}/lib)
  include_directories(${MOSAIC_INCLUDE_DIR})
  add_definitions(-DPMC_USE_MOSAIC)
endif()

######################################################################
# CAMP

if(ENABLE_CAMP)
  find_path(CAMP_INCLUDE_DIR camp_camp_core.mod
    DOC "CAMP include directory"
    PATHS $ENV{CAMP_HOME}/include)
  find_library(CAMP_LIB camp
    DOC "CAMP library"
    PATHS $ENV{CAMP_HOME}/lib)
  include_directories(${CAMP_INCLUDE_DIR})
  add_definitions(-DPMC_USE_CAMP)
endif()

######################################################################
# MPI

if(ENABLE_MPI)
  add_definitions(-DPMC_USE_MPI)
endif()

######################################################################
# SUNDIALS

if(ENABLE_SUNDIALS)
  find_path(SUNDIALS_INCLUDE_DIR cvode/cvode.h
    DOC "SUNDIALS include directory (must have cvode/, sundials/, nvector/ subdirs)"
    PATHS $ENV{SUNDIALS_HOME}/include /opt/local/include /usr/local/include)
  find_library(SUNDIALS_NVECSERIAL_LIB sundials_nvecserial
    DOC "SUNDIALS serial vector library"
    PATHS $ENV{SUNDIALS_HOME}/lib /opt/local/lib /usr/local/lib)
  find_library(SUNDIALS_CVODE_LIB sundials_cvode
    DOC "SUNDIALS CVODE library"
    PATHS $ENV{SUNDIALS_HOME}/lib /opt/local/lib /usr/local/lib)
  set(SUNDIALS_LIBS ${SUNDIALS_NVECSERIAL_LIB} ${SUNDIALS_CVODE_LIB})
  set(SUNDIALS_SRC src/condense_solver.c)
  include_directories(${SUNDIALS_INCLUDE_DIR})
  add_definitions(-DPMC_USE_SUNDIALS)
endif()

######################################################################
# tests

enable_testing()
add_custom_target(copy_test ALL ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/test ${CMAKE_BINARY_DIR}/test_run)
set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES test_run)

add_test(test_additive ${CMAKE_BINARY_DIR}/test_run/run_test_directory.sh additive)
add_test(test_average ${CMAKE_BINARY_DIR}/test_run/run_test_directory.sh average)
add_test(test_bidisperse ${CMAKE_BINARY_DIR}/test_run/run_test_directory.sh bidisperse)
add_test(test_brownian ${CMAKE_BINARY_DIR}/test_run/run_test_directory.sh brownian)
if(ENABLE_CAMP)
  add_test(test_camp ${CMAKE_BINARY_DIR}/test_run/run_test_directory.sh camp)
endif()
if(ENABLE_SUNDIALS)
  add_test(test_condense ${CMAKE_BINARY_DIR}/test_run/run_test_directory.sh condense)
endif()
add_test(test_emission ${CMAKE_BINARY_DIR}/test_run/run_test_directory.sh emission)
add_test(test_fractal ${CMAKE_BINARY_DIR}/test_run/run_test_directory.sh fractal)
add_test(test_loss ${CMAKE_BINARY_DIR}/test_run/run_test_directory.sh loss)
add_test(test_nucleate ${CMAKE_BINARY_DIR}/test_run/run_test_directory.sh nucleate)
add_test(test_mixing_state ${CMAKE_BINARY_DIR}/test_run/run_test_directory.sh mixing_state)
if(ENABLE_MOSAIC)
  add_test(test_mosaic ${CMAKE_BINARY_DIR}/test_run/run_test_directory.sh mosaic)
endif()
if(ENABLE_MPI)
  add_test(test_parallel ${CMAKE_BINARY_DIR}/test_run/run_test_directory.sh parallel)
endif()
add_test(test_rand ${CMAKE_BINARY_DIR}/test_run/run_test_directory.sh rand)
add_test(test_sedi ${CMAKE_BINARY_DIR}/test_run/run_test_directory.sh sedi)

######################################################################
# partmc library

add_library(partmclib src/aero_state.F90 src/integer_varray.F90
  src/integer_rmap.F90 src/integer_rmap2.F90 src/aero_sorted.F90
  src/aero_binned.F90 src/bin_grid.F90 src/constants.F90
  src/scenario.F90 src/env_state.F90 src/aero_mode.F90
  src/aero_dist.F90 src/aero_weight.F90 src/aero_weight_array.F90
  src/coag_kernel_additive.F90 src/coag_kernel_sedi.F90
  src/coag_kernel_constant.F90 src/coag_kernel_brown.F90
  src/coag_kernel_zero.F90 src/coag_kernel_brown_free.F90
  src/coag_kernel_brown_cont.F90 src/aero_data.F90 src/run_exact.F90
  src/run_part.F90 src/util.F90 src/stats.F90 src/run_sect.F90 src/output.F90
  src/mosaic.F90 src/gas_data.F90 src/gas_state.F90
  src/coagulation.F90 src/exact_soln.F90 src/coagulation_dist.F90
  src/coag_kernel.F90 src/spec_line.F90 src/spec_file.F90 src/rand.F90
  src/aero_particle.F90 src/aero_particle_array.F90 src/mpi.F90
  src/netcdf.F90 src/aero_info.F90 src/aero_info_array.F90
  src/nucleate.F90 src/condense.F90 src/fractal.F90 src/chamber.F90
  src/camp_interface.F90 src/photolysis.F90 src/sys.F90
  ${SUNDIALS_SRC} ${GSL_SRC}
  ${C_SORT_SRC})

target_link_libraries(partmclib ${NETCDF_LIBS} ${SUNDIALS_LIBS}
  ${MOSAIC_LIB} ${GSL_LIBS} ${CAMP_LIB})

set_target_properties(partmclib PROPERTIES OUTPUT_NAME partmc)

######################################################################
# partmc executable

add_executable(partmc src/partmc.F90)

target_link_libraries(partmc partmclib)

######################################################################
# test_bidisperse_ode

add_executable(test_bidisperse_ode
  test/bidisperse/test_bidisperse_ode.F90)

target_link_libraries(test_bidisperse_ode partmclib)

######################################################################
# test_bidisperse_extract

add_executable(test_bidisperse_extract
  test/bidisperse/test_bidisperse_extract.F90)

target_link_libraries(test_bidisperse_extract ${NETCDF_LIBS})

######################################################################
# test_nucleate_ode

add_executable(test_nucleate_ode test/nucleate/test_nucleate_ode.F90
  src/util.F90 src/constants.F90 src/sys.F90 ${C_SORT_SRC})

######################################################################
# test_poisson_sample

add_executable(test_poisson_sample test/rand/test_poisson_sample.F90
  src/util.F90 src/rand.F90 src/constants.F90 src/mpi.F90 ${GSL_SRC}
  src/sys.F90 ${C_SORT_SRC})

target_link_libraries(test_poisson_sample ${GSL_LIBS})

######################################################################
# test_binomial_sample

add_executable(test_binomial_sample test/rand/test_binomial_sample.F90
  src/util.F90 src/rand.F90 src/constants.F90 src/mpi.F90 ${GSL_SRC}
  src/sys.F90 ${C_SORT_SRC})

target_link_libraries(test_binomial_sample ${GSL_LIBS})

######################################################################
# test_fractal_self_preserve

add_executable(test_fractal_self_preserve
  test/fractal/test_fractal_self_preserve.F90
  src/getopt.F90)

target_link_libraries(test_fractal_self_preserve partmclib)

######################################################################
# test_fractal_dimless_time

add_executable(test_fractal_dimless_time
  test/fractal/test_fractal_dimless_time.F90
  src/getopt.F90)

target_link_libraries(test_fractal_dimless_time partmclib)

######################################################################
# test_fractal_radii_conversion

add_executable(test_fractal_radii_conversion
  test/fractal/test_fractal_radii_conversion.F90)

target_link_libraries(test_fractal_radii_conversion partmclib)

######################################################################
# test_mixing_state_process

add_executable(test_mixing_state_process
  test/mixing_state/test_mixing_state_process.F90)

target_link_libraries(test_mixing_state_process partmclib)

######################################################################
# bin_average_comp

add_executable(bin_average_comp src/bin_average_comp.F90)

target_link_libraries(bin_average_comp partmclib)

######################################################################
# bin_average_size

add_executable(bin_average_size src/bin_average_size.F90)

target_link_libraries(bin_average_size partmclib)

######################################################################
# extract_aero_*

add_executable(extract_aero_particles
  src/extract_aero_particles.F90 src/getopt.F90)
target_link_libraries(extract_aero_particles partmclib)

add_executable(extract_aero_size src/extract_aero_size.F90
  src/getopt.F90)
target_link_libraries(extract_aero_size partmclib)

add_executable(extract_aero_time src/extract_aero_time.F90
  src/getopt.F90)
target_link_libraries(extract_aero_time partmclib)

add_executable(extract_gas src/extract_gas.F90 src/getopt.F90)
target_link_libraries(extract_gas partmclib)

add_executable(extract_env src/extract_env.F90 src/getopt.F90)
target_link_libraries(extract_env partmclib)

######################################################################
# extract_sectional_*

add_executable(extract_sectional_aero_size
  src/extract_sectional_aero_size.F90 src/getopt.F90)
target_link_libraries(extract_sectional_aero_size partmclib)

add_executable(extract_sectional_aero_time
  src/extract_sectional_aero_time.F90 src/getopt.F90)
target_link_libraries(extract_sectional_aero_time partmclib)

######################################################################
# numeric_*

add_executable(numeric_diff src/numeric_diff.F90 src/getopt.F90)
target_link_libraries(numeric_diff partmclib)

add_executable(numeric_average src/numeric_average.F90)

######################################################################
# scenarios/1_urban_plume/urban_plum_process

add_executable(urban_plume_process
  scenarios/1_urban_plume/urban_plume_process.F90)
target_link_libraries(urban_plume_process partmclib)

######################################################################
# scenarios/4_chamber/chamber_process

add_executable(chamber_process
  scenarios/4_chamber/chamber_process.F90)
target_link_libraries(chamber_process partmclib)

######################################################################
