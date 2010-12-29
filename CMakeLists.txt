cmake_minimum_required(VERSION 2.6)
project(PARTMC Fortran C)

set(PACKAGE_BUGREPORT "mwest@illinois.edu")
set(PACKAGE_NAME "PartMC")
set(PACKAGE_STRING "PartMC 2.0.0")
set(PACKAGE_TARNAME "partmc")
set(PACKAGE_VERSION "2.0.0")

option(ENABLE_MOSAIC "Enable MOSAIC chemistry support" OFF)
option(ENABLE_MPI "Enable MPI parallel support" OFF)
option(ENABLE_SUNDIALS "Enable SUNDIALS solver for condensation support" OFF)

######################################################################

find_path(NETCDF_INCLUDE_DIR netcdf.mod NETCDF.mod
  DOC "NetCDF include directory (must contain netcdf.mod)"
  PATHS
  $ENV{NETCDF_HOME}/include
  /include
  /usr/include
  /usr/local/include
  /opt/include
  /opt/local/include
  /usr/lib64/gfortran/modules
  /usr/lib/gfortran/modules
  /sw/include
  /sw/lib/netcdf-gfortran/include
  /sw/lib/netcdf-g95/include
  /usr/include/netcdf
  $ENV{HOME}/include
  $ENV{HOME}/opt/include
  $ENV{HOME}/netcdf/include
  PATH_SUFFIXES netcdf-3)
find_library(NETCDF_C_LIB netcdf
  DOC "NetCDF C library"
  PATHS
  $ENV{NETCDF_HOME}/lib
  $ENV{NETCDF_HOME}/lib64
  /usr/lib
  /usr/lib64
  /usr/local/lib
  /usr/local/lib64
  /opt/lib
  /opt/local/lib
  /sw/lib
  /sw/lib/netcdf-gfortran/lib
  /sw/lib/netcdf-g95/lib
  $ENV{HOME}/lib
  $ENV{HOME}/lib64
  $ENV{HOME}/opt/lib
  $ENV{HOME}/opt/lib64
  $ENV{HOME}/netcdf/lib
  $ENV{HOME}/netcdf/lib64)
find_library(NETCDF_FORTRAN_LIB netcdff
  DOC "NetCDF Fortran library"
  PATHS
  $ENV{NETCDF_HOME}/lib
  $ENV{NETCDF_HOME}/lib64
  /usr/lib
  /usr/lib64
  /usr/local/lib
  /usr/local/lib64
  /opt/lib
  /opt/local/lib
  /sw/lib
  /sw/lib/netcdf-gfortran/lib
  /sw/lib/netcdf-g95/lib
  $ENV{HOME}/lib
  $ENV{HOME}/lib64
  $ENV{HOME}/opt/lib
  $ENV{HOME}/opt/lib64
  $ENV{HOME}/netcdf/lib
  $ENV{HOME}/netcdf/lib64)

set(NETCDF_LIBS ${NETCDF_C_LIB})
if(NETCDF_FORTRAN_LIB)
  set(NETCDF_LIBS ${NETCDF_LIBS} ${NETCDF_FORTRAN_LIB})
endif()

######################################################################

if(ENABLE_SUNDIALS)
  find_path(SUNDIALS_INCLUDE_DIR cvode/cvode.h
    DOC "SUNDIALS include directory (must have cvode/, sundials/, nvector/ subdirs)"
    PATHS
    /include
    /usr/include
    /usr/local/include
    /opt/include
    /usr/lib64/gfortran/modules
    /usr/lib/gfortran/modules
    /opt/include
    /opt/local/include
    /sw/include
    /sw/lib/netcdf-gfortran/include
    /sw/lib/netcdf-g95/include
    /usr/include/netcdf
    $ENV{HOME}/include
    $ENV{HOME}/opt/include
    $ENV{HOME}/sundials/include)
  find_library(SUNDIALS_NVECSERIAL_LIB sundials_nvecserial
    DOC "SUNDIALS serial vector library"
    PATHS
    /usr/lib
    /usr/lib64
    /usr/local/lib
    /usr/local/lib64
    /opt/lib
    /opt/local/lib
    /sw/lib
    /sw/lib/netcdf-gfortran/lib
    /sw/lib/netcdf-g95/lib
    $ENV{HOME}/lib
    $ENV{HOME}/lib64
    $ENV{HOME}/opt/lib
    $ENV{HOME}/opt/lib64
    $ENV{HOME}/sundials/lib
    $ENV{HOME}/sundials/lib64)
  find_library(SUNDIALS_CVODE_LIB sundials_cvode
    DOC "SUNDIALS CVODE library"
    PATHS
    /usr/lib
    /usr/lib64
    /usr/local/lib
    /usr/local/lib64
    /opt/lib
    /opt/local/lib
    /sw/lib
    /sw/lib/netcdf-gfortran/lib
    /sw/lib/netcdf-g95/lib
    $ENV{HOME}/lib
    $ENV{HOME}/lib64
    $ENV{HOME}/opt/lib
    $ENV{HOME}/opt/lib64
    $ENV{HOME}/sundials/lib
    $ENV{HOME}/sundials/lib64)
  
  set(SUNDIALS_LIBS ${SUNDIALS_NVECSERIAL_LIB} ${SUNDIALS_CVODE_LIB})
endif()

######################################################################

if(ENABLE_MOSAIC)
  find_path(MOSAIC_INCLUDE_DIR module_data_mosaic_main.mod
    DOC "MOSAIC include directory"
    PATHS
    ${PROJECT_SOURCE_DIR}/../mosaic/datamodules
    ${PROJECT_SOURCE_DIR}/../../mosaic/trunk/datamodules
    $ENV{HOME}/include
    $ENV{HOME}/opt/include
    $ENV{HOME}/mosaic/include
    $ENV{HOME}/subversion/mosaic/trunk/datamodules
    $ENV{HOME}/svn/mosaic/trunk/datamodules)
  find_library(MOSAIC_LIB mosaic
    DOC "MOSAIC library"
    PATHS
    ${PROJECT_SOURCE_DIR}/../mosaic
    ${PROJECT_SOURCE_DIR}/../../mosaic/trunk
    $ENV{HOME}/lib
    $ENV{HOME}/lib64
    $ENV{HOME}/opt/lib
    $ENV{HOME}/opt/lib64
    $ENV{HOME}/mosaic/lib
    $ENV{HOME}/mosaic/lib64
    $ENV{HOME}/subversion/mosaic/trunk
    $ENV{HOME}/svn/mosaic/trunk)
endif()

######################################################################

set(EXTRA_LIBS)

include_directories(${NETCDF_INCLUDE_DIR})
set(EXTRA_LIBS ${EXTRA_LIBS} ${NETCDF_LIBS})

if(ENABLE_MOSAIC)
  include_directories(${MOSAIC_INCLUDE_DIR})
  set(EXTRA_LIBS ${EXTRA_LIBS} ${MOSAIC_LIB})
  add_definitions(-DPMC_USE_MOSAIC)
endif()

if(ENABLE_MPI)
  add_definitions(-DPMC_USE_MPI)
endif()

if(ENABLE_SUNDIALS)
  include_directories(${SUNDIALS_INCLUDE_DIR})
  set(EXTRA_LIBS ${EXTRA_LIBS} ${SUNDIALS_LIBS})
  add_definitions(-DPMC_USE_SUNDIALS)
endif()

get_filename_component(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# MW 2009-01-20: Below we add -I/usr/include manually because gfortran
# does not search this by default [1], but cmake strips it from paths
# given by include_directories(), assuming that the compiler will
# search it anyway (as gcc does). I reported this as a cmake bug [2].
# [1] http://gcc.gnu.org/bugzilla/show_bug.cgi?id=35707
# [2] http://www.cmake.org/Bug/view.php?id=8408
#
# MW 2009-06-15: The issue is rather tricky from the gfortran side
# (module files should not go in include directories at all). On the
# cmake side the issue has been resolved by removing include directory
# suppression altogether for Fortran [3]. This change apparently went
# into cmake 2.6.4, but we'll keep the explicit include hack below for
# a while until everyone updates.
# [3] http://www.cmake.org/Bug/view.php?id=8598

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(DEFAULT_PARTMC_FFLAGS "-g -std=f2003 -pedantic -x f95-cpp-input -fimplicit-none -W -Wall -Wconversion -Wunderflow -Wno-unused -fbounds-check")
elseif(COMPILER_NAME MATCHES ".*pgf90" OR COMPILER_NAME MATCHES ".*pgf95")
  set(DEFAULT_PARTMC_FFLAGS "-Mpreprocess -DDEFINE_LOCAL_COMMAND_ARGUMENT")
else(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(DEFAULT_PARTMC_FFLAGS "")
endif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

set(PARTMC_CFLAGS "" CACHE STRING "CFLAGS used while compiling PartMC")
set(PARTMC_FFLAGS ${DEFAULT_PARTMC_FFLAGS} CACHE STRING "FFLAGS used while compiling PartMC")

######################################################################

enable_testing()
add_custom_target(copy_test ALL ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/test ${CMAKE_BINARY_DIR}/run_test)
set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES run_test)

add_test(additive_1 ${CMAKE_BINARY_DIR}/run_test/additive/run_test_1.sh)
add_test(additive_2 ${CMAKE_BINARY_DIR}/run_test/additive/run_test_2.sh)

add_test(average_01 ${CMAKE_BINARY_DIR}/run_test/average/run_test_01.sh)
add_test(average_02 ${CMAKE_BINARY_DIR}/run_test/average/run_test_02.sh)
add_test(average_03 ${CMAKE_BINARY_DIR}/run_test/average/run_test_03.sh)
add_test(average_04 ${CMAKE_BINARY_DIR}/run_test/average/run_test_04.sh)
add_test(average_05 ${CMAKE_BINARY_DIR}/run_test/average/run_test_05.sh)
add_test(average_06 ${CMAKE_BINARY_DIR}/run_test/average/run_test_06.sh)
add_test(average_07 ${CMAKE_BINARY_DIR}/run_test/average/run_test_07.sh)
add_test(average_08 ${CMAKE_BINARY_DIR}/run_test/average/run_test_08.sh)
add_test(average_09 ${CMAKE_BINARY_DIR}/run_test/average/run_test_09.sh)
add_test(average_10 ${CMAKE_BINARY_DIR}/run_test/average/run_test_10.sh)
add_test(average_11 ${CMAKE_BINARY_DIR}/run_test/average/run_test_11.sh)
add_test(average_12 ${CMAKE_BINARY_DIR}/run_test/average/run_test_12.sh)
add_test(average_13 ${CMAKE_BINARY_DIR}/run_test/average/run_test_13.sh)
add_test(average_14 ${CMAKE_BINARY_DIR}/run_test/average/run_test_14.sh)
add_test(average_15 ${CMAKE_BINARY_DIR}/run_test/average/run_test_15.sh)
add_test(average_16 ${CMAKE_BINARY_DIR}/run_test/average/run_test_16.sh)
add_test(average_17 ${CMAKE_BINARY_DIR}/run_test/average/run_test_17.sh)
add_test(average_18 ${CMAKE_BINARY_DIR}/run_test/average/run_test_18.sh)
add_test(average_19 ${CMAKE_BINARY_DIR}/run_test/average/run_test_19.sh)
add_test(average_20 ${CMAKE_BINARY_DIR}/run_test/average/run_test_20.sh)

add_test(bidisperse_1 ${CMAKE_BINARY_DIR}/run_test/bidisperse/run_test_1.sh)
add_test(bidisperse_2 ${CMAKE_BINARY_DIR}/run_test/bidisperse/run_test_2.sh)
add_test(bidisperse_3 ${CMAKE_BINARY_DIR}/run_test/bidisperse/run_test_3.sh)

add_test(brownian_1 ${CMAKE_BINARY_DIR}/run_test/brownian/run_test_1.sh)
add_test(brownian_2 ${CMAKE_BINARY_DIR}/run_test/brownian/run_test_2.sh)
add_test(brownian_3 ${CMAKE_BINARY_DIR}/run_test/brownian/run_test_3.sh)
add_test(brownian_4 ${CMAKE_BINARY_DIR}/run_test/brownian/run_test_4.sh)
add_test(brownian_5 ${CMAKE_BINARY_DIR}/run_test/brownian/run_test_5.sh)
add_test(brownian_6 ${CMAKE_BINARY_DIR}/run_test/brownian/run_test_6.sh)

if(ENABLE_SUNDIALS)
  add_test(condense_1 ${CMAKE_BINARY_DIR}/run_test/condense/run_test_1.sh)
  add_test(condense_2 ${CMAKE_BINARY_DIR}/run_test/condense/run_test_2.sh)
endif()

add_test(emission_1 ${CMAKE_BINARY_DIR}/run_test/emission/run_test_1.sh)
add_test(emission_2 ${CMAKE_BINARY_DIR}/run_test/emission/run_test_2.sh)
add_test(emission_3 ${CMAKE_BINARY_DIR}/run_test/emission/run_test_3.sh)
add_test(emission_4 ${CMAKE_BINARY_DIR}/run_test/emission/run_test_4.sh)
add_test(emission_5 ${CMAKE_BINARY_DIR}/run_test/emission/run_test_5.sh)

add_test(nucleate_1 ${CMAKE_BINARY_DIR}/run_test/nucleate/run_test_1.sh)
add_test(nucleate_2 ${CMAKE_BINARY_DIR}/run_test/nucleate/run_test_2.sh)
add_test(nucleate_3 ${CMAKE_BINARY_DIR}/run_test/nucleate/run_test_3.sh)

if(ENABLE_MOSAIC)
  add_test(mosaic_01 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_01.sh)
  add_test(mosaic_02 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_02.sh)
  add_test(mosaic_03 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_03.sh)
  add_test(mosaic_04 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_04.sh)
  add_test(mosaic_05 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_05.sh)
  add_test(mosaic_06 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_06.sh)
  add_test(mosaic_07 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_07.sh)
  add_test(mosaic_08 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_08.sh)
  add_test(mosaic_09 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_09.sh)
  add_test(mosaic_10 ${CMAKE_BINARY_DIR}/run_test/mosaic/run_test_10.sh)
endif()

if(ENABLE_MPI)
  add_test(parallel_1 ${CMAKE_BINARY_DIR}/run_test/parallel/run_test_1.sh)
  add_test(parallel_2 ${CMAKE_BINARY_DIR}/run_test/parallel/run_test_2.sh)
  add_test(parallel_3 ${CMAKE_BINARY_DIR}/run_test/parallel/run_test_3.sh)
  add_test(parallel_4 ${CMAKE_BINARY_DIR}/run_test/parallel/run_test_4.sh)
  add_test(parallel_5 ${CMAKE_BINARY_DIR}/run_test/parallel/run_test_5.sh)
  add_test(parallel_6 ${CMAKE_BINARY_DIR}/run_test/parallel/run_test_6.sh)
endif()

add_test(poisson_1 ${CMAKE_BINARY_DIR}/run_test/poisson/run_test_1.sh)
add_test(poisson_2 ${CMAKE_BINARY_DIR}/run_test/poisson/run_test_2.sh)
add_test(poisson_3 ${CMAKE_BINARY_DIR}/run_test/poisson/run_test_3.sh)
add_test(poisson_4 ${CMAKE_BINARY_DIR}/run_test/poisson/run_test_4.sh)
add_test(poisson_5 ${CMAKE_BINARY_DIR}/run_test/poisson/run_test_5.sh)

add_test(sedi_1 ${CMAKE_BINARY_DIR}/run_test/sedi/run_test_1.sh)
add_test(sedi_2 ${CMAKE_BINARY_DIR}/run_test/sedi/run_test_2.sh)

######################################################################

set(PARTMC_SRC src/partmc.f90 src/aero_state.f90 src/aero_binned.f90
  src/bin_grid.f90 src/constants.f90 src/env_data.f90
  src/env_state.f90 src/aero_mode.f90 src/aero_dist.f90
  src/aero_weight.f90 src/coag_kernel_additive.f90
  src/coag_kernel_sedi.f90 src/coag_kernel_constant.f90
  src/coag_kernel_brown.f90 src/coag_kernel_zero.f90 src/aero_data.f90
  src/run_exact.f90 src/run_part.f90 src/util.f90 src/run_sect.f90
  src/output.f90 src/mosaic.f90 src/gas_data.f90 src/gas_state.f90
  src/coagulation.f90 src/exact_soln.f90 src/coagulation_mpi_equal.f90
  src/coag_kernel.f90 src/spec_line.f90 src/spec_file.f90 src/rand.f90
  src/aero_particle.f90 src/aero_particle_array.f90 src/mpi.f90
  src/netcdf.f90 src/aero_info.f90 src/aero_info_array.f90
  src/nucleate.f90)

if(PARTMC_FFLAGS)
  set_source_files_properties(${PARTMC_SRC}
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

if(ENABLE_SUNDIALS)
  if(PARTMC_CFLAGS)
    set_source_files_properties(src/condense_solver.c
      PROPERTIES COMPILE_FLAGS ${PARTMC_CFLAGS})
  endif()
  if(PARTMC_FFLAGS)
    set_source_files_properties(src/condense.f90
      PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
    set(PARTMC_SRC ${PARTMC_SRC} src/condense.f90 src/condense_solver.c)
  endif(PARTMC_FFLAGS)
endif()

add_executable(partmc ${PARTMC_SRC})
target_link_libraries(partmc ${EXTRA_LIBS})

######################################################################

add_executable(test_bidisperse_ode
  test/bidisperse/test_bidisperse_ode.f90 src/coag_kernel_sedi.f90
  src/env_data.f90 src/env_state.f90 src/constants.f90
  src/aero_data.f90 src/util.f90 src/gas_data.f90 src/gas_state.f90
  src/aero_state.f90 src/output.f90 src/bin_grid.f90 src/spec_line.f90
  src/spec_file.f90 src/aero_mode.f90 src/aero_dist.f90
  src/aero_weight.f90 src/aero_binned.f90 src/rand.f90
  src/aero_particle.f90 src/aero_particle_array.f90 src/mpi.f90
  src/netcdf.f90 src/aero_info.f90 src/aero_info_array.f90)
target_link_libraries(test_bidisperse_ode ${EXTRA_LIBS})
if(PARTMC_FFLAGS)
  set_source_files_properties(test/bidisperse/test_bidisperse_ode.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

add_executable(test_bidisperse_extract
  test/bidisperse/test_bidisperse_extract.f90)
target_link_libraries(test_bidisperse_extract ${NETCDF_LIBS})
if(PARTMC_FFLAGS)
  set_source_files_properties(test/bidisperse/test_bidisperse_extract.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

add_executable(test_nucleate_ode test/nucleate/test_nucleate_ode.f90
  src/util.f90 src/constants.f90)
if(PARTMC_FFLAGS)
  set_source_files_properties(test/nucleate/test_nucleate_ode.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

add_executable(test_poisson_sample test/poisson/test_poisson_sample.f90
  src/util.f90 src/rand.f90 src/constants.f90)
target_link_libraries(test_poisson_sample)
if(PARTMC_FFLAGS)
  set_source_files_properties(test/poisson/test_poisson_sample.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

######################################################################

add_executable(extract_aero_particle_mass
  src/extract_aero_particle_mass.f90)
target_link_libraries(extract_aero_particle_mass ${NETCDF_LIBS})
if(PARTMC_FFLAGS)
  set_source_files_properties(src/extract_aero_particle_mass.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

add_executable(extract_aero_size_num
  src/extract_aero_size_num.f90)
target_link_libraries(extract_aero_size_num ${NETCDF_LIBS})
if(PARTMC_FFLAGS)
  set_source_files_properties(src/extract_aero_size_num.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

add_executable(extract_aero_size_mass
  src/extract_aero_size_mass.f90)
target_link_libraries(extract_aero_size_mass ${NETCDF_LIBS})
if(PARTMC_FFLAGS)
  set_source_files_properties(src/extract_aero_size_mass.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

add_executable(extract_aero_species
  src/extract_aero_species.f90)
target_link_libraries(extract_aero_species ${NETCDF_LIBS})
if(PARTMC_FFLAGS)
  set_source_files_properties(src/extract_aero_species.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

add_executable(extract_aero_total
  src/extract_aero_total.f90)
target_link_libraries(extract_aero_total ${NETCDF_LIBS})
if(PARTMC_FFLAGS)
  set_source_files_properties(src/extract_aero_total.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

add_executable(extract_gas src/extract_gas.f90)
target_link_libraries(extract_gas ${NETCDF_LIBS})
if(PARTMC_FFLAGS)
  set_source_files_properties(src/extract_gas.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

add_executable(extract_env src/extract_env.f90)
target_link_libraries(extract_env ${NETCDF_LIBS})
if(PARTMC_FFLAGS)
  set_source_files_properties(src/extract_env.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

######################################################################

add_executable(bin_average_comp src/bin_average_comp.f90
  src/bin_grid.f90 src/gas_state.f90 src/gas_data.f90 src/util.f90
  src/env_state.f90 src/output.f90 src/aero_state.f90
  src/aero_data.f90 src/netcdf.f90 src/mpi.f90 src/spec_file.f90
  src/aero_binned.f90 src/aero_particle.f90
  src/aero_particle_array.f90 src/aero_info.f90
  src/aero_info_array.f90 src/aero_dist.f90 src/constants.f90
  src/spec_line.f90 src/aero_mode.f90 src/aero_weight.f90
  src/rand.f90)
target_link_libraries(bin_average_comp ${NETCDF_LIBS})
if(PARTMC_FFLAGS)
  set_source_files_properties(src/bin_average_comp.f90 PROPERTIES
    COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

add_executable(bin_average_size src/bin_average_size.f90
  src/bin_grid.f90 src/gas_state.f90 src/gas_data.f90 src/util.f90
  src/env_state.f90 src/output.f90 src/aero_state.f90
  src/aero_data.f90 src/netcdf.f90 src/mpi.f90 src/spec_file.f90
  src/aero_binned.f90 src/aero_particle.f90
  src/aero_particle_array.f90 src/aero_info.f90
  src/aero_info_array.f90 src/aero_dist.f90 src/constants.f90
  src/spec_line.f90 src/aero_mode.f90 src/rand.f90
  src/aero_weight.f90)
target_link_libraries(bin_average_size ${NETCDF_LIBS})
if(PARTMC_FFLAGS)
  set_source_files_properties(src/bin_average_size.f90 PROPERTIES
    COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

######################################################################

add_executable(extract_sectional_aero_size_num
  src/extract_sectional_aero_size_num.f90)
target_link_libraries(extract_sectional_aero_size_num ${NETCDF_LIBS})
if(PARTMC_FFLAGS)
  set_source_files_properties(src/extract_sectional_aero_size_num.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

add_executable(extract_sectional_aero_size_mass
  src/extract_sectional_aero_size_mass.f90)
target_link_libraries(extract_sectional_aero_size_mass ${NETCDF_LIBS})
if(PARTMC_FFLAGS)
  set_source_files_properties(src/extract_sectional_aero_size_mass.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

add_executable(extract_sectional_aero_species
  src/extract_sectional_aero_species.f90)
target_link_libraries(extract_sectional_aero_species ${NETCDF_LIBS})
if(PARTMC_FFLAGS)
  set_source_files_properties(src/extract_sectional_aero_species.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

add_executable(extract_sectional_aero_total
  src/extract_sectional_aero_total.f90)
target_link_libraries(extract_sectional_aero_total ${NETCDF_LIBS})
if(PARTMC_FFLAGS)
  set_source_files_properties(src/extract_sectional_aero_total.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

######################################################################

add_executable(numeric_diff src/numeric_diff.f90)
if(PARTMC_FFLAGS)
  set_source_files_properties(src/numeric_diff.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

add_executable(numeric_average src/numeric_average.f90)
if(PARTMC_FFLAGS)
  set_source_files_properties(src/numeric_average.f90
    PROPERTIES COMPILE_FLAGS ${PARTMC_FFLAGS})
endif(PARTMC_FFLAGS)

######################################################################
