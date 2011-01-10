cmake_minimum_required(VERSION 2.6.4)
# Require 2.6.4 due to -I/usr/include behavior:
# http://www.cmake.org/Bug/view.php?id=8598
# http://gcc.gnu.org/bugzilla/show_bug.cgi?id=35707
# http://www.cmake.org/Bug/view.php?id=8408
project(PARTMC Fortran C)

set(PACKAGE_BUGREPORT "mwest@illinois.edu")
set(PACKAGE_NAME "PartMC")
set(PACKAGE_STRING "PartMC 2.0.0")
set(PACKAGE_TARNAME "partmc")
set(PACKAGE_VERSION "2.0.0")

######################################################################
# options

option(ENABLE_GSL "Enable GSL library for random number generation" OFF)
option(ENABLE_MOSAIC "Enable MOSAIC chemistry support" OFF)
option(ENABLE_MPI "Enable MPI parallel support" OFF)
option(ENABLE_SUNDIALS "Enable SUNDIALS solver for condensation support" OFF)

######################################################################
# CPack

set(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${PACKAGE_TARNAME}-${PACKAGE_VERSION}")
set(CPACK_SOURCE_IGNORE_FILES "${CPACK_SOURCE_IGNORE_FILES};/.*~$;/[.].*/;/build/;/figures/;/scenarios/[^12].*/;/doc/condensation/;/doc/ship_track/;/old/;/tool/")

include(CPack)

######################################################################
# NetCDF

find_path(NETCDF_INCLUDE_DIR netcdf.mod NETCDF.mod
  DOC "NetCDF include directory (must contain netcdf.mod)"
  PATHS $ENV{NETCDF_HOME}/include /opt/local/include)
find_library(NETCDF_C_LIB netcdf
  DOC "NetCDF C library"
  PATHS $ENV{NETCDF_HOME}/lib /opt/local/lib)
find_library(NETCDF_FORTRAN_LIB netcdff
  DOC "NetCDF Fortran library"
  PATHS $ENV{NETCDF_HOME}/lib /opt/local/lib)
set(NETCDF_LIBS ${NETCDF_C_LIB})
if(NETCDF_FORTRAN_LIB)
  set(NETCDF_LIBS ${NETCDF_LIBS} ${NETCDF_FORTRAN_LIB})
endif()
include_directories(${NETCDF_INCLUDE_DIR})

######################################################################
# GSL

if(ENABLE_GSL)
  find_path(GSL_INCLUDE_DIR gsl/gsl_math.h
    DOC "GSL include directory (must have gsl/ subdir)"
    PATHS $ENV{GSL_HOME}/include /opt/local/include)
  find_library(GSL_LIB gsl
    DOC "GSL library"
    PATHS $ENV{GSL_HOME}/lib /opt/local/lib)
  set(GSL_SRC src/rand_gsl.c)
  include_directories(${GSL_INCLUDE_DIR})
  add_definitions(-DPMC_USE_GSL)
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
# MPI

if(ENABLE_MPI)
  add_definitions(-DPMC_USE_MPI)
endif()

######################################################################
# SUNDIALS

if(ENABLE_SUNDIALS)
  find_path(SUNDIALS_INCLUDE_DIR cvode/cvode.h
    DOC "SUNDIALS include directory (must have cvode/, sundials/, nvector/ subdirs)"
    PATHS $ENV{SUNDIALS_HOME}/include /opt/local/include)
  find_library(SUNDIALS_NVECSERIAL_LIB sundials_nvecserial
    DOC "SUNDIALS serial vector library"
    PATHS $ENV{SUNDIALS_HOME}/lib /opt/local/lib)
  find_library(SUNDIALS_CVODE_LIB sundials_cvode
    DOC "SUNDIALS CVODE library"
    PATHS $ENV{SUNDIALS_HOME}/lib /opt/local/lib)
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

add_test(test_additive_1 ${CMAKE_BINARY_DIR}/test_run/additive/test_additive_1.sh)
add_test(test_additive_2 ${CMAKE_BINARY_DIR}/test_run/additive/test_additive_2.sh)

add_test(test_average_01 ${CMAKE_BINARY_DIR}/test_run/average/test_average_01.sh)
add_test(test_average_02 ${CMAKE_BINARY_DIR}/test_run/average/test_average_02.sh)
add_test(test_average_03 ${CMAKE_BINARY_DIR}/test_run/average/test_average_03.sh)
add_test(test_average_04 ${CMAKE_BINARY_DIR}/test_run/average/test_average_04.sh)
add_test(test_average_05 ${CMAKE_BINARY_DIR}/test_run/average/test_average_05.sh)
add_test(test_average_06 ${CMAKE_BINARY_DIR}/test_run/average/test_average_06.sh)
add_test(test_average_07 ${CMAKE_BINARY_DIR}/test_run/average/test_average_07.sh)
add_test(test_average_08 ${CMAKE_BINARY_DIR}/test_run/average/test_average_08.sh)
add_test(test_average_09 ${CMAKE_BINARY_DIR}/test_run/average/test_average_09.sh)
add_test(test_average_10 ${CMAKE_BINARY_DIR}/test_run/average/test_average_10.sh)
add_test(test_average_11 ${CMAKE_BINARY_DIR}/test_run/average/test_average_11.sh)
add_test(test_average_12 ${CMAKE_BINARY_DIR}/test_run/average/test_average_12.sh)
add_test(test_average_13 ${CMAKE_BINARY_DIR}/test_run/average/test_average_13.sh)
add_test(test_average_14 ${CMAKE_BINARY_DIR}/test_run/average/test_average_14.sh)
add_test(test_average_15 ${CMAKE_BINARY_DIR}/test_run/average/test_average_15.sh)
add_test(test_average_16 ${CMAKE_BINARY_DIR}/test_run/average/test_average_16.sh)
add_test(test_average_17 ${CMAKE_BINARY_DIR}/test_run/average/test_average_17.sh)
add_test(test_average_18 ${CMAKE_BINARY_DIR}/test_run/average/test_average_18.sh)
add_test(test_average_19 ${CMAKE_BINARY_DIR}/test_run/average/test_average_19.sh)
add_test(test_average_20 ${CMAKE_BINARY_DIR}/test_run/average/test_average_20.sh)

add_test(test_bidisperse_1 ${CMAKE_BINARY_DIR}/test_run/bidisperse/test_bidisperse_1.sh)
add_test(test_bidisperse_2 ${CMAKE_BINARY_DIR}/test_run/bidisperse/test_bidisperse_2.sh)
add_test(test_bidisperse_3 ${CMAKE_BINARY_DIR}/test_run/bidisperse/test_bidisperse_3.sh)

add_test(test_brownian_1 ${CMAKE_BINARY_DIR}/test_run/brownian/test_brownian_1.sh)
add_test(test_brownian_2 ${CMAKE_BINARY_DIR}/test_run/brownian/test_brownian_2.sh)
add_test(test_brownian_3 ${CMAKE_BINARY_DIR}/test_run/brownian/test_brownian_3.sh)
add_test(test_brownian_4 ${CMAKE_BINARY_DIR}/test_run/brownian/test_brownian_4.sh)
add_test(test_brownian_5 ${CMAKE_BINARY_DIR}/test_run/brownian/test_brownian_5.sh)
add_test(test_brownian_6 ${CMAKE_BINARY_DIR}/test_run/brownian/test_brownian_6.sh)

if(ENABLE_SUNDIALS)
  add_test(test_condense_1 ${CMAKE_BINARY_DIR}/test_run/condense/test_condense_1.sh)
  add_test(test_condense_2 ${CMAKE_BINARY_DIR}/test_run/condense/test_condense_2.sh)
endif()

add_test(test_emission_1 ${CMAKE_BINARY_DIR}/test_run/emission/test_emission_1.sh)
add_test(test_emission_2 ${CMAKE_BINARY_DIR}/test_run/emission/test_emission_2.sh)
add_test(test_emission_3 ${CMAKE_BINARY_DIR}/test_run/emission/test_emission_3.sh)
add_test(test_emission_4 ${CMAKE_BINARY_DIR}/test_run/emission/test_emission_4.sh)
add_test(test_emission_5 ${CMAKE_BINARY_DIR}/test_run/emission/test_emission_5.sh)

add_test(test_nucleate_1 ${CMAKE_BINARY_DIR}/test_run/nucleate/test_nucleate_1.sh)
add_test(test_nucleate_2 ${CMAKE_BINARY_DIR}/test_run/nucleate/test_nucleate_2.sh)
add_test(test_nucleate_3 ${CMAKE_BINARY_DIR}/test_run/nucleate/test_nucleate_3.sh)

if(ENABLE_MOSAIC)
  add_test(test_mosaic_01 ${CMAKE_BINARY_DIR}/test_run/mosaic/test_mosaic_01.sh)
  add_test(test_mosaic_02 ${CMAKE_BINARY_DIR}/test_run/mosaic/test_mosaic_02.sh)
  add_test(test_mosaic_03 ${CMAKE_BINARY_DIR}/test_run/mosaic/test_mosaic_03.sh)
  add_test(test_mosaic_04 ${CMAKE_BINARY_DIR}/test_run/mosaic/test_mosaic_04.sh)
  add_test(test_mosaic_05 ${CMAKE_BINARY_DIR}/test_run/mosaic/test_mosaic_05.sh)
  add_test(test_mosaic_06 ${CMAKE_BINARY_DIR}/test_run/mosaic/test_mosaic_06.sh)
  add_test(test_mosaic_07 ${CMAKE_BINARY_DIR}/test_run/mosaic/test_mosaic_07.sh)
  add_test(test_mosaic_08 ${CMAKE_BINARY_DIR}/test_run/mosaic/test_mosaic_08.sh)
  add_test(test_mosaic_09 ${CMAKE_BINARY_DIR}/test_run/mosaic/test_mosaic_09.sh)
  add_test(test_mosaic_10 ${CMAKE_BINARY_DIR}/test_run/mosaic/test_mosaic_10.sh)
endif()

if(ENABLE_MPI)
  add_test(test_parallel_1 ${CMAKE_BINARY_DIR}/test_run/parallel/test_parallel_1.sh)
  add_test(test_parallel_2 ${CMAKE_BINARY_DIR}/test_run/parallel/test_parallel_2.sh)
  add_test(test_parallel_3 ${CMAKE_BINARY_DIR}/test_run/parallel/test_parallel_3.sh)
  add_test(test_parallel_4 ${CMAKE_BINARY_DIR}/test_run/parallel/test_parallel_4.sh)
endif()

add_test(test_rand_1 ${CMAKE_BINARY_DIR}/test_run/rand/test_rand_1.sh)
add_test(test_rand_2 ${CMAKE_BINARY_DIR}/test_run/rand/test_rand_2.sh)
add_test(test_rand_3 ${CMAKE_BINARY_DIR}/test_run/rand/test_rand_3.sh)
add_test(test_rand_4 ${CMAKE_BINARY_DIR}/test_run/rand/test_rand_4.sh)
add_test(test_rand_5 ${CMAKE_BINARY_DIR}/test_run/rand/test_rand_5.sh)
add_test(test_rand_6 ${CMAKE_BINARY_DIR}/test_run/rand/test_rand_6.sh)
add_test(test_rand_7 ${CMAKE_BINARY_DIR}/test_run/rand/test_rand_7.sh)
add_test(test_rand_8 ${CMAKE_BINARY_DIR}/test_run/rand/test_rand_8.sh)

add_test(test_sedi_1 ${CMAKE_BINARY_DIR}/test_run/sedi/test_sedi_1.sh)
add_test(test_sedi_2 ${CMAKE_BINARY_DIR}/test_run/sedi/test_sedi_2.sh)

######################################################################
# partmc

add_executable(partmc src/partmc.F95 src/aero_state.F95
  src/aero_binned.F95 src/bin_grid.F95 src/constants.F95
  src/env_data.F95 src/env_state.F95 src/aero_mode.F95
  src/aero_dist.F95 src/aero_weight.F95 src/coag_kernel_additive.F95
  src/coag_kernel_sedi.F95 src/coag_kernel_constant.F95
  src/coag_kernel_brown.F95 src/coag_kernel_zero.F95 src/aero_data.F95
  src/run_exact.F95 src/run_part.F95 src/util.F95 src/run_sect.F95
  src/output.F95 src/mosaic.F95 src/gas_data.F95 src/gas_state.F95
  src/coagulation.F95 src/exact_soln.F95 src/coagulation_dist.F95
  src/coag_kernel.F95 src/spec_line.F95 src/spec_file.F95 src/rand.F95
  src/aero_particle.F95 src/aero_particle_array.F95 src/mpi.F95
  src/netcdf.F95 src/aero_info.F95 src/aero_info_array.F95
  src/nucleate.F95 src/condense.F95 ${SUNDIALS_SRC} ${GSL_SRC})

target_link_libraries(partmc ${NETCDF_LIBS} ${SUNDIALS_LIBS}
  ${MOSAIC_LIB} ${GSL_LIB})

######################################################################
# test_bidisperse_ode

add_executable(test_bidisperse_ode
  test/bidisperse/test_bidisperse_ode.F95 src/coag_kernel_sedi.F95
  src/env_data.F95 src/env_state.F95 src/constants.F95
  src/aero_data.F95 src/util.F95 src/gas_data.F95 src/gas_state.F95
  src/aero_state.F95 src/output.F95 src/bin_grid.F95 src/spec_line.F95
  src/spec_file.F95 src/aero_mode.F95 src/aero_dist.F95
  src/aero_weight.F95 src/aero_binned.F95 src/rand.F95
  src/aero_particle.F95 src/aero_particle_array.F95 src/mpi.F95
  src/netcdf.F95 src/aero_info.F95 src/aero_info_array.F95 ${GSL_SRC})

target_link_libraries(test_bidisperse_ode ${NETCDF_LIBS} ${GSL_LIB})

######################################################################
# test_bidisperse_extract

add_executable(test_bidisperse_extract
test/bidisperse/test_bidisperse_extract.F95)

target_link_libraries(test_bidisperse_extract ${NETCDF_LIBS})

######################################################################
# test_nucleate_ode

set(TEST_NUCLEATE_ODE_SRC test/nucleate/test_nucleate_ode.F95
  src/util.F95 src/constants.F95)

add_executable(test_nucleate_ode ${TEST_NUCLEATE_ODE_SRC})

######################################################################
# test_poisson_sample

add_executable(test_poisson_sample test/rand/test_poisson_sample.F95
  src/util.F95 src/rand.F95 src/constants.F95 src/mpi.F95 ${GSL_SRC})

target_link_libraries(test_poisson_sample ${GSL_LIB})

######################################################################
# test_binomial_sample

add_executable(test_binomial_sample test/rand/test_binomial_sample.F95
  src/util.F95 src/rand.F95 src/constants.F95 src/mpi.F95 ${GSL_SRC})

target_link_libraries(test_binomial_sample ${GSL_LIB})

######################################################################
# bin_average_comp

add_executable(bin_average_comp src/bin_average_comp.F95
  src/bin_grid.F95 src/gas_state.F95 src/gas_data.F95 src/util.F95
  src/env_state.F95 src/output.F95 src/aero_state.F95
  src/aero_data.F95 src/netcdf.F95 src/mpi.F95 src/spec_file.F95
  src/aero_binned.F95 src/aero_particle.F95
  src/aero_particle_array.F95 src/aero_info.F95
  src/aero_info_array.F95 src/aero_dist.F95 src/constants.F95
  src/spec_line.F95 src/aero_mode.F95 src/aero_weight.F95 src/rand.F95
  ${GSL_SRC})

target_link_libraries(bin_average_comp ${NETCDF_LIBS} ${GSL_LIB})

######################################################################
# bin_average_size

add_executable(bin_average_size src/bin_average_size.F95
  src/bin_grid.F95 src/gas_state.F95 src/gas_data.F95 src/util.F95
  src/env_state.F95 src/output.F95 src/aero_state.F95
  src/aero_data.F95 src/netcdf.F95 src/mpi.F95 src/spec_file.F95
  src/aero_binned.F95 src/aero_particle.F95
  src/aero_particle_array.F95 src/aero_info.F95
  src/aero_info_array.F95 src/aero_dist.F95 src/constants.F95
  src/spec_line.F95 src/aero_mode.F95 src/aero_weight.F95 src/rand.F95
  ${GSL_SRC})

target_link_libraries(bin_average_size ${NETCDF_LIBS} ${GSL_LIB})

######################################################################
# extract_aero_*

add_executable(extract_aero_particle_mass
  src/extract_aero_particle_mass.F95)
target_link_libraries(extract_aero_particle_mass ${NETCDF_LIBS})

add_executable(extract_aero_size_num
  src/extract_aero_size_num.F95)
target_link_libraries(extract_aero_size_num ${NETCDF_LIBS})

add_executable(extract_aero_size_mass
  src/extract_aero_size_mass.F95)
target_link_libraries(extract_aero_size_mass ${NETCDF_LIBS})

add_executable(extract_aero_species
  src/extract_aero_species.F95)
target_link_libraries(extract_aero_species ${NETCDF_LIBS})

add_executable(extract_aero_total
  src/extract_aero_total.F95)
target_link_libraries(extract_aero_total ${NETCDF_LIBS})

add_executable(extract_gas src/extract_gas.F95)
target_link_libraries(extract_gas ${NETCDF_LIBS})

add_executable(extract_env src/extract_env.F95)
target_link_libraries(extract_env ${NETCDF_LIBS})

######################################################################
# extract_sectional_*

add_executable(extract_sectional_aero_size_num
  src/extract_sectional_aero_size_num.F95)
target_link_libraries(extract_sectional_aero_size_num ${NETCDF_LIBS})

add_executable(extract_sectional_aero_size_mass
  src/extract_sectional_aero_size_mass.F95)
target_link_libraries(extract_sectional_aero_size_mass ${NETCDF_LIBS})

add_executable(extract_sectional_aero_species
  src/extract_sectional_aero_species.F95)
target_link_libraries(extract_sectional_aero_species ${NETCDF_LIBS})

add_executable(extract_sectional_aero_total
  src/extract_sectional_aero_total.F95)
target_link_libraries(extract_sectional_aero_total ${NETCDF_LIBS})

######################################################################
# numeric_*

add_executable(numeric_diff src/numeric_diff.F95)

add_executable(numeric_average src/numeric_average.F95)

######################################################################
