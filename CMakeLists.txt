cmake_minimum_required(VERSION 2.6.4)
# Require 2.6.4 due to -I/usr/include behavior:
# http://www.cmake.org/Bug/view.php?id=8598
# http://gcc.gnu.org/bugzilla/show_bug.cgi?id=35707
# http://www.cmake.org/Bug/view.php?id=8408
project(PARTMC Fortran C)

set(PACKAGE_BUGREPORT "mwest@illinois.edu")
set(PACKAGE_NAME "PartMC")
set(PACKAGE_STRING "PartMC 2.4.0")
set(PACKAGE_TARNAME "partmc")
set(PACKAGE_VERSION "2.4.0")

######################################################################
# options

option(ENABLE_GSL "Enable GSL library for random number generation" OFF)
option(ENABLE_MOSAIC "Enable MOSAIC chemistry support" OFF)
option(ENABLE_MPI "Enable MPI parallel support" OFF)
option(ENABLE_SUNDIALS "Enable SUNDIALS solver for condensation support" OFF)
option(ENABLE_C_SORT "Enable C sorting routines" OFF)
option(ENABLE_JSON "Enable json-fortran support" ON)
option(ENABLE_TESTS "Enable tests" ON)
option(ENABLE_DEBUG "Compile debugging functions" OFF)
option(FAILURE_DETAIL "Output conditions before and after solver failures" OFF)
option(ENABLE_CXX "Enable C++" OFF)
option(ENABLE_GPU "Enable use of GPUs in chemistry solving" OFF)

mark_as_advanced(FORCE ENABLE_DEBUG FAILURE_DETAIL)

######################################################################
# CPack

set(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${PACKAGE_TARNAME}-${PACKAGE_VERSION}")
set(CPACK_SOURCE_IGNORE_FILES "${CPACK_SOURCE_IGNORE_FILES};/.*~$;/[.].*/;/build/;/figures/;/scenarios/[^1234].*/;/doc/condensation/;/doc/deposition/;/doc/ship_track/;/old/;/tool/;/TODO")

include(CPack)

######################################################################
# NetCDF

find_path(NETCDF_INCLUDE_DIR netcdf.mod NETCDF.mod
  DOC "NetCDF include directory (must contain netcdf.mod)"
  PATHS
  $ENV{NETCDF_HOME}/include
  /usr/lib/gfortran/modules
  /usr/lib64/gfortran/modules
  /opt/local/include)
find_library(NETCDF_C_LIB netcdf
  DOC "NetCDF C library"
  PATHS $ENV{NETCDF_HOME}/lib $ENV{NETCDF_HOME}/lib64 /opt/local/lib)
find_library(NETCDF_FORTRAN_LIB netcdff
  DOC "NetCDF Fortran library"
  PATHS $ENV{NETCDF_HOME}/lib $ENV{NETCDF_HOME}/lib64 /opt/local/lib)
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
# MPI

if(ENABLE_MPI)
  add_definitions(-DPMC_USE_MPI)
endif()

######################################################################
# SUNDIALS

if(ENABLE_SUNDIALS)
  find_path(SUITE_SPARSE_INCLUDE_DIR klu.h
    DOC "SuiteSparse include directory (must have klu.h)"
    PATHS $ENV{SUITE_SPARSE_HOME}/include $ENV{SUNDIALS_HOME}/include
          /opt/local/include /usr/local/include)
  find_library(SUITE_SPARSE_KLU_LIB klu
    DOC "SuiteSparse klu library"
    PATHS $ENV{SUITE_SPARSE_HOME}/lib $ENV{SUNDIALS_HOME}/lib
          /opt/local/lib /usr/local/lib)
  find_path(SUNDIALS_INCLUDE_DIR cvode/cvode.h
    DOC "SUNDIALS include directory (must have cvode/, sundials/, nvector/ subdirs)"
    PATHS $ENV{SUNDIALS_HOME}/include /opt/local/include /usr/local/include)
  find_library(SUNDIALS_NVECSERIAL_LIB sundials_nvecserial
    DOC "SUNDIALS serial vector library"
    PATHS $ENV{SUNDIALS_HOME}/lib /opt/local/lib /usr/local/lib)
  find_library(SUNDIALS_CVODE_LIB sundials_cvode
    DOC "SUNDIALS CVODE library"
    PATHS $ENV{SUNDIALS_HOME}/lib /opt/local/lib /usr/local/lib)
  find_library(SUNDIALS_KLU_LIB sundials_sunlinsolklu
    DOC "SUNDIALS KLU library"
    PATHS $ENV{SUITE_SPARSE_HOME}/lib $ENV{SUNDIALS_HOME}/lib
          /opt/local/lib /usr/local/lib)
  find_library(SUNDIALS_SUNMATRIX_SPARSE_LIB sundials_sunmatrixsparse
    DOC "SUNDIALS SUNMatrixSparse library"
    PATHS $ENV{SUITE_SPARSE_HOME}/lib $ENV{SUNDIALS_HOME}/lib
          /opt/local/lib /usr/local/lib)
  set(SUNDIALS_LIBS ${SUNDIALS_NVECSERIAL_LIB} ${SUNDIALS_CVODE_LIB}
    ${SUITE_SPARSE_KLU_LIB} ${SUNDIALS_KLU_LIB} ${SUNDIALS_SUNMATRIX_SPARSE_LIB})
  set(SUNDIALS_SRC src/condense_solver.c)
  include_directories(${SUNDIALS_INCLUDE_DIR} ${SUITE_SPARSE_INCLUDE_DIR})
  add_definitions(-DPMC_USE_SUNDIALS)
endif()

######################################################################
# json-fortran

if (ENABLE_JSON)
  find_path(JSON_INCLUDE_DIR json_module.mod
    DOC "json-fortran include directory (must include json_*.mod files)"
    PATHS $ENV{JSON_FORTRAN_HOME}/lib /opt/local/lib /usr/local/lib)
  find_library(JSON_LIB jsonfortran
    DOC "json-fortran library"
    PATHS $ENV{JSON_FORTRAN_HOME}/lib /opt/local/lib /usr/local/lib)
  include_directories(${JSON_INCLUDE_DIR})
  add_definitions(-DPMC_USE_JSON)
endif()

######################################################################
##CXX

if(ENABLE_CXX)
  add_definitions(-DPMC_USE_CXX)
  enable_language(CXX)
endif()

######################################################################
# GPU support


#Compile GPU tests if CUDA detected
find_package(CUDA)
if (CUDA_FOUND)
    enable_language(CUDA)
    #If enable_gpu, all will be executed on GPU
    if (ENABLE_GPU)
        #Define global variable of use gpu no matter what test to include user main test
        add_definitions(-DPMC_USE_GPU)
    endif()
endif()

if(NOT CUDA_FOUND AND ENABLE_GPU)
    message( WARNING "Not CUDA module found, can't support GPU computation")
endif()

######################################################################
# debugging options

if (ENABLE_DEBUG)
  add_definitions(-DPMC_DEBUG)
endif()
if (FAILURE_DETAIL)
  add_definitions(-DFAILURE_DETAIL)
endif()

######################################################################
# debugging options

if (ENABLE_DEBUG)
  add_definitions(-DPMC_DEBUG)
endif()
if (FAILURE_DETAIL)
  add_definitions(-DFAILURE_DETAIL)
endif()

######################################################################
# Unit test macro

macro(do_unit_test test_name result)
  if(ENABLE_MPI)
    add_test(unit_test_${test_name} mpirun -v -np 2 ${CMAKE_BINARY_DIR}/unit_test_${test_name})
  else()
    add_test(unit_test_${test_name} ${CMAKE_BINARY_DIR}/unit_test_${test_name})
  endif()
  set_tests_properties(unit_test_${test_name}
	  PROPERTIES PASS_REGULAR_EXPRESSION ${result})
endmacro(do_unit_test)

######################################################################
# copy dirs

add_custom_target(copy_test ALL ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/test ${CMAKE_BINARY_DIR}/test_run)
add_custom_target(copy_mechanism ALL ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/mechanisms ${CMAKE_BINARY_DIR}/mechanisms_run)
add_custom_target(copy_scenario ALL ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/scenarios ${CMAKE_BINARY_DIR}/scenarios_run)
set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES test_run mechanisms_run scenarios_run)

######################################################################
# tests

enable_testing()
do_unit_test(property "PASS")
do_unit_test(chem_spec_data "PASS")
do_unit_test(aero_phase_data "PASS")
if(ENABLE_JSON)
  do_unit_test(aero_rep_single_particle "PASS")
  do_unit_test(aero_rep_modal_binned_mass "PASS")
  do_unit_test(camp_core "PASS")
endif()

if (ENABLE_MPI)
  set(MPI_TEST_FLAG MPI)
else()
  set(MPI_TEST_FLAG serial)
endif()

if(CUDA_FOUND)
    add_test(test_rxn_arrhenius_mech_gpu ${CMAKE_BINARY_DIR}/test_run/unit_tests/input_files/run_rxn_arrhenius_gpu.sh ${MPI_TEST_FLAG})
endif()


add_test(test_rxn_arrhenius_mech ${CMAKE_BINARY_DIR}/test_run/unit_tests/input_files/run_rxn_arrhenius.sh ${MPI_TEST_FLAG})

if (ENABLE_TESTS)
    if (ENABLE_JSON AND ENABLE_SUNDIALS)
      # New unit tests (UNDER DEVELOPMENT)
      add_test(test_rxn_arrhenius_mech ${CMAKE_BINARY_DIR}/test_run/unit_tests/input_files/run_rxn_arrhenius.sh ${MPI_TEST_FLAG})

      # Old-style unit tests
      add_test(test_rxn_aqueous_equilibrium ${CMAKE_BINARY_DIR}/test_run/unit_rxn_data/test_aqueous_equilibrium.sh ${MPI_TEST_FLAG})
      add_test(test_rxn_CMAQ_H2O2 ${CMAKE_BINARY_DIR}/test_run/unit_rxn_data/test_CMAQ_H2O2.sh ${MPI_TEST_FLAG})
      add_test(test_rxn_CMAQ_OH_HNO3 ${CMAKE_BINARY_DIR}/test_run/unit_rxn_data/test_CMAQ_OH_HNO3.sh ${MPI_TEST_FLAG})
      add_test(test_rxn_condensed_phase_arrhenius ${CMAKE_BINARY_DIR}/test_run/unit_rxn_data/test_condensed_phase_arrhenius.sh ${MPI_TEST_FLAG})
      add_test(test_rxn_emission ${CMAKE_BINARY_DIR}/test_run/unit_rxn_data/test_emission.sh ${MPI_TEST_FLAG})
      add_test(test_rxn_first_order_loss ${CMAKE_BINARY_DIR}/test_run/unit_rxn_data/test_first_order_loss.sh ${MPI_TEST_FLAG})
      add_test(test_rxn_HL_phase_transfer ${CMAKE_BINARY_DIR}/test_run/unit_rxn_data/test_HL_phase_transfer.sh ${MPI_TEST_FLAG})
      add_test(test_rxn_photolysis ${CMAKE_BINARY_DIR}/test_run/unit_rxn_data/test_photolysis.sh ${MPI_TEST_FLAG})
      add_test(test_rxn_SIMPOL_phase_transfer ${CMAKE_BINARY_DIR}/test_run/unit_rxn_data/test_SIMPOL_phase_transfer.sh ${MPI_TEST_FLAG})
      add_test(test_rxn_troe ${CMAKE_BINARY_DIR}/test_run/unit_rxn_data/test_troe.sh ${MPI_TEST_FLAG})
      add_test(test_rxn_wet_deposition ${CMAKE_BINARY_DIR}/test_run/unit_rxn_data/test_wet_deposition.sh ${MPI_TEST_FLAG})
      add_test(test_sub_model_PDFiTE ${CMAKE_BINARY_DIR}/test_run/unit_sub_model_data/test_PDFiTE.sh ${MPI_TEST_FLAG})
      add_test(test_sub_model_UNIFAC ${CMAKE_BINARY_DIR}/test_run/unit_sub_model_data/test_UNIFAC.sh ${MPI_TEST_FLAG})
      add_test(test_sub_model_ZSR_aerosol_water ${CMAKE_BINARY_DIR}/test_run/unit_sub_model_data/test_ZSR_aerosol_water.sh ${MPI_TEST_FLAG})
      add_test(test_chem_mech_solver ${CMAKE_BINARY_DIR}/test_run/chemistry/test_chemistry_1.sh ${MPI_TEST_FLAG})
      add_test(test_chemistry_cb05cl_ae5 ${CMAKE_BINARY_DIR}/test_run/chemistry/cb05cl_ae5/test_chemistry_cb05cl_ae5.sh ${MPI_TEST_FLAG})
      add_test(test_chemistry_eqsam_v03d ${CMAKE_BINARY_DIR}/test_run/chemistry/eqsam_v03d/test_chemistry_eqsam_v03d.sh ${MPI_TEST_FLAG})
      add_test(test_MONARCH_1 ${CMAKE_BINARY_DIR}/test_run/monarch/test_monarch_1.sh ${MPI_TEST_FLAG})
      add_test(test_MONARCH_2 ${CMAKE_BINARY_DIR}/test_run/monarch/test_monarch_2.sh ${MPI_TEST_FLAG})

    endif()

    add_test(test_additive_1 ${CMAKE_BINARY_DIR}/test_run/additive/test_additive_1.sh)
    add_test(test_additive_2 ${CMAKE_BINARY_DIR}/test_run/additive/test_additive_2.sh)
    set_tests_properties(test_additive_2 PROPERTIES DEPENDS test_additive_1)

    add_test(test_average_01 ${CMAKE_BINARY_DIR}/test_run/average/test_average_01.sh)
    add_test(test_average_02 ${CMAKE_BINARY_DIR}/test_run/average/test_average_02.sh)
    set_tests_properties(test_average_02 PROPERTIES DEPENDS test_average_01)
    add_test(test_average_03 ${CMAKE_BINARY_DIR}/test_run/average/test_average_03.sh)
    set_tests_properties(test_average_03 PROPERTIES DEPENDS test_average_02)
    add_test(test_average_04 ${CMAKE_BINARY_DIR}/test_run/average/test_average_04.sh)
    set_tests_properties(test_average_04 PROPERTIES DEPENDS test_average_03)
    add_test(test_average_05 ${CMAKE_BINARY_DIR}/test_run/average/test_average_05.sh)
    set_tests_properties(test_average_05 PROPERTIES DEPENDS test_average_04)
    add_test(test_average_06 ${CMAKE_BINARY_DIR}/test_run/average/test_average_06.sh)
    set_tests_properties(test_average_06 PROPERTIES DEPENDS test_average_05)
    add_test(test_average_07 ${CMAKE_BINARY_DIR}/test_run/average/test_average_07.sh)
    set_tests_properties(test_average_07 PROPERTIES DEPENDS test_average_06)
    add_test(test_average_08 ${CMAKE_BINARY_DIR}/test_run/average/test_average_08.sh)
    set_tests_properties(test_average_08 PROPERTIES DEPENDS test_average_07)
    add_test(test_average_09 ${CMAKE_BINARY_DIR}/test_run/average/test_average_09.sh)
    set_tests_properties(test_average_09 PROPERTIES DEPENDS test_average_08)
    add_test(test_average_10 ${CMAKE_BINARY_DIR}/test_run/average/test_average_10.sh)
    set_tests_properties(test_average_10 PROPERTIES DEPENDS test_average_09)
    add_test(test_average_11 ${CMAKE_BINARY_DIR}/test_run/average/test_average_11.sh)
    set_tests_properties(test_average_11 PROPERTIES DEPENDS test_average_10)
    add_test(test_average_12 ${CMAKE_BINARY_DIR}/test_run/average/test_average_12.sh)
    set_tests_properties(test_average_12 PROPERTIES DEPENDS test_average_11)
    add_test(test_average_13 ${CMAKE_BINARY_DIR}/test_run/average/test_average_13.sh)
    set_tests_properties(test_average_13 PROPERTIES DEPENDS test_average_12)
    add_test(test_average_14 ${CMAKE_BINARY_DIR}/test_run/average/test_average_14.sh)
    set_tests_properties(test_average_14 PROPERTIES DEPENDS test_average_13)
    add_test(test_average_15 ${CMAKE_BINARY_DIR}/test_run/average/test_average_15.sh)
    set_tests_properties(test_average_15 PROPERTIES DEPENDS test_average_14)
    add_test(test_average_16 ${CMAKE_BINARY_DIR}/test_run/average/test_average_16.sh)
    set_tests_properties(test_average_16 PROPERTIES DEPENDS test_average_15)
    add_test(test_average_17 ${CMAKE_BINARY_DIR}/test_run/average/test_average_17.sh)
    set_tests_properties(test_average_17 PROPERTIES DEPENDS test_average_16)
    add_test(test_average_18 ${CMAKE_BINARY_DIR}/test_run/average/test_average_18.sh)
    set_tests_properties(test_average_18 PROPERTIES DEPENDS test_average_17)
    add_test(test_average_19 ${CMAKE_BINARY_DIR}/test_run/average/test_average_19.sh)
    set_tests_properties(test_average_19 PROPERTIES DEPENDS test_average_18)
    add_test(test_average_20 ${CMAKE_BINARY_DIR}/test_run/average/test_average_20.sh)
    set_tests_properties(test_average_20 PROPERTIES DEPENDS test_average_19)

    add_test(test_bidisperse_1 ${CMAKE_BINARY_DIR}/test_run/bidisperse/test_bidisperse_1.sh)

    add_test(test_brownian_1 ${CMAKE_BINARY_DIR}/test_run/brownian/test_brownian_1.sh)
    add_test(test_brownian_2 ${CMAKE_BINARY_DIR}/test_run/brownian/test_brownian_2.sh)
    set_tests_properties(test_brownian_2 PROPERTIES DEPENDS test_brownian_1)
    add_test(test_brownian_3 ${CMAKE_BINARY_DIR}/test_run/brownian/test_brownian_3.sh)
    set_tests_properties(test_brownian_3 PROPERTIES DEPENDS test_brownian_2)
    add_test(test_brownian_4 ${CMAKE_BINARY_DIR}/test_run/brownian/test_brownian_4.sh)
    set_tests_properties(test_brownian_4 PROPERTIES DEPENDS test_brownian_3)

    if(ENABLE_SUNDIALS)
      add_test(test_condense_1 ${CMAKE_BINARY_DIR}/test_run/condense/test_condense_1.sh)
      add_test(test_condense_2 ${CMAKE_BINARY_DIR}/test_run/condense/test_condense_2.sh)
      set_tests_properties(test_condense_2 PROPERTIES DEPENDS test_condense_1)
    endif()

    add_test(test_emission_1 ${CMAKE_BINARY_DIR}/test_run/emission/test_emission_1.sh)
    add_test(test_emission_2 ${CMAKE_BINARY_DIR}/test_run/emission/test_emission_2.sh)
    set_tests_properties(test_emission_2 PROPERTIES DEPENDS test_emission_1)
    add_test(test_emission_3 ${CMAKE_BINARY_DIR}/test_run/emission/test_emission_3.sh)
    set_tests_properties(test_emission_3 PROPERTIES DEPENDS test_emission_2)

    add_test(test_fractal_1 ${CMAKE_BINARY_DIR}/test_run/fractal/test_fractal_1.sh)
    add_test(test_fractal_2 ${CMAKE_BINARY_DIR}/test_run/fractal/test_fractal_2.sh)
    set_tests_properties(test_fractal_2 PROPERTIES DEPENDS test_fractal_1)
    add_test(test_fractal_3 ${CMAKE_BINARY_DIR}/test_run/fractal/test_fractal_3.sh)
    set_tests_properties(test_fractal_3 PROPERTIES DEPENDS test_fractal_2)

    add_test(test_loss_01 ${CMAKE_BINARY_DIR}/test_run/loss/test_loss_01.sh)
    add_test(test_loss_02 ${CMAKE_BINARY_DIR}/test_run/loss/test_loss_02.sh)
    add_test(test_loss_03 ${CMAKE_BINARY_DIR}/test_run/loss/test_loss_03.sh)
    set_tests_properties(test_loss_02 PROPERTIES DEPENDS test_loss_01)
    set_tests_properties(test_loss_03 PROPERTIES DEPENDS test_loss_01)
    add_test(test_loss_04 ${CMAKE_BINARY_DIR}/test_run/loss/test_loss_04.sh)
    add_test(test_loss_05 ${CMAKE_BINARY_DIR}/test_run/loss/test_loss_05.sh)
    add_test(test_loss_06 ${CMAKE_BINARY_DIR}/test_run/loss/test_loss_06.sh)
    set_tests_properties(test_loss_05 PROPERTIES DEPENDS test_loss_04)
    set_tests_properties(test_loss_06 PROPERTIES DEPENDS test_loss_04)
    add_test(test_loss_07 ${CMAKE_BINARY_DIR}/test_run/loss/test_loss_07.sh)
    add_test(test_loss_08 ${CMAKE_BINARY_DIR}/test_run/loss/test_loss_08.sh)
    add_test(test_loss_09 ${CMAKE_BINARY_DIR}/test_run/loss/test_loss_09.sh)
    set_tests_properties(test_loss_08 PROPERTIES DEPENDS test_loss_07)
    set_tests_properties(test_loss_09 PROPERTIES DEPENDS test_loss_07)
    add_test(test_loss_10 ${CMAKE_BINARY_DIR}/test_run/loss/test_loss_10.sh)
    add_test(test_loss_11 ${CMAKE_BINARY_DIR}/test_run/loss/test_loss_11.sh)
    add_test(test_loss_12 ${CMAKE_BINARY_DIR}/test_run/loss/test_loss_12.sh)
    set_tests_properties(test_loss_11 PROPERTIES DEPENDS test_loss_10)
    set_tests_properties(test_loss_12 PROPERTIES DEPENDS test_loss_10)

    add_test(test_nucleate_1 ${CMAKE_BINARY_DIR}/test_run/nucleate/test_nucleate_1.sh)
    add_test(test_nucleate_2 ${CMAKE_BINARY_DIR}/test_run/nucleate/test_nucleate_2.sh)
    set_tests_properties(test_nucleate_2 PROPERTIES DEPENDS test_nucleate_1)

    if(ENABLE_MOSAIC)
      add_test(test_mosaic_1 ${CMAKE_BINARY_DIR}/test_run/mosaic/test_mosaic_1.sh)
      add_test(test_mosaic_2 ${CMAKE_BINARY_DIR}/test_run/mosaic/test_mosaic_2.sh)
      set_tests_properties(test_mosaic_2 PROPERTIES DEPENDS test_mosaic_1)
      add_test(test_mosaic_3 ${CMAKE_BINARY_DIR}/test_run/mosaic/test_mosaic_3.sh)
      set_tests_properties(test_mosaic_3 PROPERTIES DEPENDS test_mosaic_2)
      add_test(test_mosaic_4 ${CMAKE_BINARY_DIR}/test_run/mosaic/test_mosaic_4.sh)
      set_tests_properties(test_mosaic_4 PROPERTIES DEPENDS test_mosaic_3)
      add_test(test_mosaic_5 ${CMAKE_BINARY_DIR}/test_run/mosaic/test_mosaic_5.sh)
      set_tests_properties(test_mosaic_5 PROPERTIES DEPENDS test_mosaic_4)
      add_test(test_mosaic_6 ${CMAKE_BINARY_DIR}/test_run/mosaic/test_mosaic_6.sh)
      set_tests_properties(test_mosaic_6 PROPERTIES DEPENDS test_mosaic_5)
    endif()

    if(ENABLE_MPI)
      add_test(test_parallel_1 ${CMAKE_BINARY_DIR}/test_run/parallel/test_parallel_1.sh)
      #add_test(test_parallel_2 ${CMAKE_BINARY_DIR}/test_run/parallel/test_parallel_2.sh)
      #set_tests_properties(test_parallel_2 PROPERTIES DEPENDS test_parallel_1)
      add_test(test_parallel_3 ${CMAKE_BINARY_DIR}/test_run/parallel/test_parallel_3.sh)
      set_tests_properties(test_parallel_3 PROPERTIES DEPENDS test_parallel_1)
      #set_tests_properties(test_parallel_3 PROPERTIES DEPENDS test_parallel_2)
      #add_test(test_parallel_4 ${CMAKE_BINARY_DIR}/test_run/parallel/test_parallel_4.sh)
      #set_tests_properties(test_parallel_4 PROPERTIES DEPENDS test_parallel_3)
    endif()

    add_test(test_rand_1 ${CMAKE_BINARY_DIR}/test_run/rand/test_rand_1.sh)
    add_test(test_rand_2 ${CMAKE_BINARY_DIR}/test_run/rand/test_rand_2.sh)
    set_tests_properties(test_rand_2 PROPERTIES DEPENDS test_rand_1)
    add_test(test_rand_3 ${CMAKE_BINARY_DIR}/test_run/rand/test_rand_3.sh)
    set_tests_properties(test_rand_3 PROPERTIES DEPENDS test_rand_2)
    add_test(test_rand_4 ${CMAKE_BINARY_DIR}/test_run/rand/test_rand_4.sh)
    set_tests_properties(test_rand_4 PROPERTIES DEPENDS test_rand_3)
    add_test(test_rand_5 ${CMAKE_BINARY_DIR}/test_run/rand/test_rand_5.sh)
    set_tests_properties(test_rand_5 PROPERTIES DEPENDS test_rand_4)
    add_test(test_rand_6 ${CMAKE_BINARY_DIR}/test_run/rand/test_rand_6.sh)
    set_tests_properties(test_rand_6 PROPERTIES DEPENDS test_rand_5)
    add_test(test_rand_7 ${CMAKE_BINARY_DIR}/test_run/rand/test_rand_7.sh)
    set_tests_properties(test_rand_7 PROPERTIES DEPENDS test_rand_6)
    add_test(test_rand_8 ${CMAKE_BINARY_DIR}/test_run/rand/test_rand_8.sh)
    set_tests_properties(test_rand_8 PROPERTIES DEPENDS test_rand_7)

    add_test(test_sedi_1 ${CMAKE_BINARY_DIR}/test_run/sedi/test_sedi_1.sh)
    add_test(test_sedi_2 ${CMAKE_BINARY_DIR}/test_run/sedi/test_sedi_2.sh)
    set_tests_properties(test_sedi_2 PROPERTIES DEPENDS test_sedi_1)
endif()

######################################################################
# partmc library

set(STD_C_FLAGS "-std=c99")

set(STD_CUDA_FLAGS "-dc -arch=compute_70 -code=sm_70")




if(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set(STD_F_FLAGS "-extend_source -warn truncated_source")
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set(STD_F_FLAGS "-ffree-line-length-none")
endif()

set(REACTIONS_F_SRC
	src/rxns/rxn_aqueous_equilibrium.F90
	src/rxns/rxn_arrhenius.F90
	src/rxns/rxn_CMAQ_H2O2.F90
	src/rxns/rxn_CMAQ_OH_HNO3.F90
	src/rxns/rxn_condensed_phase_arrhenius.F90
	src/rxns/rxn_emission.F90
	src/rxns/rxn_first_order_loss.F90
	src/rxns/rxn_HL_phase_transfer.F90
	src/rxns/rxn_photolysis.F90
        src/rxns/rxn_SIMPOL_phase_transfer.F90
	src/rxns/rxn_troe.F90
        src/rxns/rxn_wet_deposition.F90)

set(REACTIONS_C_SRC
	src/rxns/rxn_aqueous_equilibrium.c
	src/rxns/rxn_arrhenius.c
	src/rxns/rxn_CMAQ_H2O2.c
	src/rxns/rxn_CMAQ_OH_HNO3.c
	src/rxns/rxn_condensed_phase_arrhenius.c
	src/rxns/rxn_emission.c
	src/rxns/rxn_first_order_loss.c
	src/rxns/rxn_HL_phase_transfer.c
	src/rxns/rxn_photolysis.c
        src/rxns/rxn_SIMPOL_phase_transfer.c
	src/rxns/rxn_troe.c
        src/rxns/rxn_wet_deposition.c)

set_source_files_properties(${REACTIONS_F_SRC} PROPERTIES COMPILE_FLAGS
        ${STD_F_FLAGS})

set_source_files_properties(${REACTIONS_C_SRC} PROPERTIES COMPILE_FLAGS
        ${STD_C_FLAGS})

set(REACTIONS_SRC ${REACTIONS_F_SRC} ${REACTIONS_C_SRC})

set(AEROSOL_REPS_F_SRC
        src/aero_reps/aero_rep_modal_binned_mass.F90
        src/aero_reps/aero_rep_single_particle.F90)

set(AEROSOL_REPS_C_SRC
        src/aero_reps/aero_rep_modal_binned_mass.c
        src/aero_reps/aero_rep_single_particle.c)

set_source_files_properties(${AEROSOL_REPS_F_SRC} PROPERTIES COMPILE_FLAGS
        ${STD_F_FLAGS})

set_source_files_properties(${AEROSOL_REPS_C_SRC} PROPERTIES COMPILE_FLAGS
        ${STD_C_FLAGS})

set(AEROSOL_REPS_SRC ${AEROSOL_REPS_F_SRC} ${AEROSOL_REPS_C_SRC})

set(SUB_MODELS_F_SRC
  src/sub_models/sub_model_PDFiTE.F90
        src/sub_models/sub_model_UNIFAC.F90
        src/sub_models/sub_model_ZSR_aerosol_water.F90)

set(SUB_MODELS_C_SRC
  src/sub_models/sub_model_PDFiTE.c
        src/sub_models/sub_model_UNIFAC.c
        src/sub_models/sub_model_ZSR_aerosol_water.c)

set_source_files_properties(${SUB_MODELS_F_SRC} PROPERTIES COMPILE_FLAGS
        ${STD_F_FLAGS})

set_source_files_properties(${SUB_MODELS_C_SRC} PROPERTIES COMPILE_FLAGS
        ${STD_C_FLAGS})

set(SUB_MODELS_SRC ${SUB_MODELS_F_SRC} ${SUB_MODELS_C_SRC})

set(CAMP_C_SRC
        src/camp_solver.c src/rxn_solver.c src/aero_phase_solver.c
        src/aero_rep_solver.c src/sub_model_solver.c)

set_source_files_properties(${CAMP_C_SRC} PROPERTIES COMPILE_FLAGS
        ${STD_C_FLAGS})

set(CAMP_CXX_SRC "")

if (CUDA_FOUND)
#if(ENABLE_GPU)
  set(CAMP_CUDA_SRC
          src/cuda/camp_gpu_solver.cu
          src/cuda/rxns_gpu/rxn_aqueous_equilibrium.cu
          src/cuda/rxns_gpu/rxn_arrhenius.cu
          src/cuda/rxns_gpu/rxn_CMAQ_H2O2.cu
          src/cuda/rxns_gpu/rxn_CMAQ_OH_HNO3.cu
          src/cuda/rxns_gpu/rxn_condensed_phase_arrhenius.cu
          src/cuda/rxns_gpu/rxn_emission.cu
          src/cuda/rxns_gpu/rxn_first_order_loss.cu
          src/cuda/rxns_gpu/rxn_HL_phase_transfer.cu
          src/cuda/rxns_gpu/rxn_photolysis.cu
          src/cuda/rxns_gpu/rxn_SIMPOL_phase_transfer.cu
          src/cuda/rxns_gpu/rxn_troe.cu
          src/cuda/rxns_gpu/rxn_wet_deposition.cu
          src/cuda/aeros_gpu/aero_phase_solver.cu
          src/cuda/aeros_gpu/aero_rep_solver.cu
          src/cuda/aeros_gpu/sub_model_solver.cu
          src/cuda/aero_reps_gpu/aero_rep_modal_binned_mass.cu
          src/cuda/aero_reps_gpu/aero_rep_single_particle.cu)

  set_source_files_properties(${CAMP_CUDA_SRC} PROPERTIES COMPILE_FLAGS
        ${STD_CUDA_FLAGS})

  set_source_files_properties(${CAMP_CUDA_SRC} PROPERTIES LANGUAGE CUDA)
else()
#  set(CAMP_CUDA_SRC "")
endif()

set(CAMP_F_SRC
        src/aero_state.F90 src/integer_varray.F90
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
        src/property.F90 src/chem_spec_data.F90
        src/rxn_data.F90 src/camp_state.F90 src/mechanism_data.F90
        src/camp_core.F90 src/camp_solver_data.F90 src/aero_rep_data.F90
        src/aero_phase_data.F90 src/aero_rep_factory.F90 src/camp_interface.F90
        src/rxn_factory.F90 src/sub_model_data.F90 src/sub_model_factory.F90
        src/solver_stats.F90 src/camp_box_model_data.F90
        )

add_library(partmclib ${CAMP_F_SRC}
  ${CAMP_C_SRC} ${AEROSOL_REPS_SRC} ${SUB_MODELS_SRC} ${REACTIONS_SRC}
  ${SUNDIALS_SRC} ${GSL_SRC} ${C_SORT_SRC} ${CAMP_CUDA_SRC} ${CAMP_CXX_SRC} )

target_link_libraries(partmclib ${NETCDF_LIBS} ${SUNDIALS_LIBS}
  ${MOSAIC_LIB} ${GSL_LIBS} ${JSON_LIB})

set_target_properties(partmclib PROPERTIES OUTPUT_NAME partmc)

if (CUDA_FOUND)
    add_library(partmclib_gpu ${CAMP_F_SRC}
            ${CAMP_C_SRC} ${AEROSOL_REPS_SRC} ${SUB_MODELS_SRC} ${REACTIONS_SRC}
            ${SUNDIALS_SRC} ${GSL_SRC} ${C_SORT_SRC} ${CAMP_CUDA_SRC} ${CAMP_CXX_SRC} )

    target_compile_definitions(partmclib_gpu PRIVATE
            "PMC_USE_GPU=ON")

    target_link_libraries(partmclib_gpu ${NETCDF_LIBS} ${SUNDIALS_LIBS}
            ${MOSAIC_LIB} ${GSL_LIBS} ${JSON_LIB})

    set_target_properties(partmclib_gpu PROPERTIES OUTPUT_NAME partmc_gpu)

    add_executable(partmc_gpu src/partmc.F90)

    target_link_libraries(partmc_gpu partmclib_gpu)
endif()

######################################################################
# partmc executable

add_executable(partmc src/partmc.F90)

target_link_libraries(partmc partmclib)

######################################################################
# camp-chem box model

add_executable(camp_box_model src/camp_box_model.F90)

target_link_libraries(camp_box_model partmclib)

######################################################################
# test_chemistry_cb05cl_ae5

set(CB5_EBI_SOLVER
	test/chemistry/cb05cl_ae5/EXT-CONST_mod.F90
	test/chemistry/cb05cl_ae5/ext-hrg2.F90
	test/chemistry/cb05cl_ae5/ext-hrcalcks.F90
	test/chemistry/cb05cl_ae5/ext-hrinit.F90
	test/chemistry/cb05cl_ae5/ext-hrdata_mod.F90
	test/chemistry/cb05cl_ae5/ext-hrg3.F90
	test/chemistry/cb05cl_ae5/ext-hrprodloss.F90
	test/chemistry/cb05cl_ae5/ext-hrg1.F90
	test/chemistry/cb05cl_ae5/ext-hrg4.F90
	test/chemistry/cb05cl_ae5/ext-hrrates.F90
	test/chemistry/cb05cl_ae5/ext-hrsolver.F90
	test/chemistry/cb05cl_ae5/EXT-RXCM_mod.F90
	)

set(CB5_KPP_SOLVER
	test/chemistry/cb05cl_ae5/cb05cl_ae5_Function.f90
	test/chemistry/cb05cl_ae5/cb05cl_ae5_Global.f90
	test/chemistry/cb05cl_ae5/cb05cl_ae5_Initialize.f90
        test/chemistry/cb05cl_ae5/cb05cl_ae5_Integrator.F90
	test/chemistry/cb05cl_ae5/cb05cl_ae5_Jacobian.f90
	test/chemistry/cb05cl_ae5/cb05cl_ae5_JacobianSP.f90
	test/chemistry/cb05cl_ae5/cb05cl_ae5_LinearAlgebra.f90
	test/chemistry/cb05cl_ae5/cb05cl_ae5_Model.f90
	test/chemistry/cb05cl_ae5/cb05cl_ae5_Monitor.f90
	test/chemistry/cb05cl_ae5/cb05cl_ae5_Parameters.f90
	test/chemistry/cb05cl_ae5/cb05cl_ae5_Precision.f90
	test/chemistry/cb05cl_ae5/cb05cl_ae5_Rates.f90
	test/chemistry/cb05cl_ae5/cb05cl_ae5_Util.f90
	)

if(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
	set(CB5_EBI_FLAGS "-fixed -extend_source -warn truncated_source")
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
	set(CB5_EBI_FLAGS "-ffixed-form -ffixed-line-length-none")
endif()

set_source_files_properties(${CB5_EBI_SOLVER} PROPERTIES COMPILE_FLAGS
	${CB5_EBI_FLAGS})

######################################################################
# test_chemistry_cb05cl_ae5_big

add_executable(test_chemistry_cb05cl_ae5_big
        test/chemistry/cb05cl_ae5/test_cb05cl_ae5_big.F90
        ${CB5_EBI_SOLVER} ${CB5_KPP_SOLVER}
        test/chemistry/cb05cl_ae5/module_BSC_CHEM_DATA.F90)
target_link_libraries(test_chemistry_cb05cl_ae5_big partmclib)

######################################################################
# MONARCH interface

add_executable(mock_monarch test/monarch/mock_monarch.F90
        test/monarch/pmc_monarch_interface.F90)

target_link_libraries(mock_monarch partmclib)


set(UNIT_TEST_SRC
        test/unit_tests/unit_test_data.F90
        test/unit_tests/unit_test_driver.F90)

set(UNIT_TEST_RXN_ARRHENIUS_SRC
        test/unit_tests/rxns/unit_test_rxn_arrhenius.F90 ${UNIT_TEST_SRC})
set_source_files_properties(${UNIT_TEST_RXN_ARRHENIUS_SRC} PROPERTIES COMPILE_FLAGS ${STD_F_FLAGS})
add_executable(unit_test_rxn_arrhenius ${UNIT_TEST_RXN_ARRHENIUS_SRC})
target_compile_definitions(unit_test_rxn_arrhenius PRIVATE
        "UNIT_TEST_MODULE_=pmc_unit_test_rxn_arrhenius \
  ;UNIT_TEST_TYPE_=unit_test_rxn_arrhenius_t()")
target_link_libraries(unit_test_rxn_arrhenius partmclib)

if (CUDA_FOUND)
    add_executable(unit_test_rxn_arrhenius_gpu ${UNIT_TEST_RXN_ARRHENIUS_SRC})
    target_compile_definitions(unit_test_rxn_arrhenius_gpu PRIVATE
            "UNIT_TEST_MODULE_=pmc_unit_test_rxn_arrhenius \
    ;UNIT_TEST_TYPE_=unit_test_rxn_arrhenius_t()")
    target_link_libraries(unit_test_rxn_arrhenius_gpu partmclib_gpu)
endif()

if (ENABLE_TESTS)

######################################################################
# test_chemistry_cb05cl_ae5

add_executable(test_chemistry_cb05cl_ae5
    test/chemistry/cb05cl_ae5/test_cb05cl_ae5.F90
    ${CB5_EBI_SOLVER} ${CB5_KPP_SOLVER}
    test/chemistry/cb05cl_ae5/module_BSC_CHEM_DATA.F90)
target_link_libraries(test_chemistry_cb05cl_ae5 partmclib)

######################################################################
# test_chemistry_eqsam_v03d

set(EQSAM_V03D
	test/chemistry/eqsam_v03d/eqsam_v03d.F90
	)

add_executable(test_chemistry_eqsam_v03d test/chemistry/eqsam_v03d/test_eqsam_v03d.F90
	${EQSAM_V03D})

target_link_libraries(test_chemistry_eqsam_v03d partmclib)

######################################################################
# test_property

add_executable(unit_test_property test/unit_property/test_property.F90)

target_link_libraries(unit_test_property partmclib)

######################################################################
# test_chem_spec_data

add_executable(unit_test_chem_spec_data test/unit_chem_spec_data/test_chem_spec_data.F90)

target_link_libraries(unit_test_chem_spec_data partmclib)

######################################################################
# test_rxn_data

add_executable(test_rxn_aqueous_equilibrium test/unit_rxn_data/test_rxn_aqueous_equilibrium.F90)
target_link_libraries(test_rxn_aqueous_equilibrium partmclib)
add_executable(test_rxn_CMAQ_H2O2 test/unit_rxn_data/test_rxn_CMAQ_H2O2.F90)
target_link_libraries(test_rxn_CMAQ_H2O2 partmclib)
add_executable(test_rxn_CMAQ_OH_HNO3 test/unit_rxn_data/test_rxn_CMAQ_OH_HNO3.F90)
target_link_libraries(test_rxn_CMAQ_OH_HNO3 partmclib)
add_executable(test_rxn_condensed_phase_arrhenius test/unit_rxn_data/test_rxn_condensed_phase_arrhenius.F90)
target_link_libraries(test_rxn_condensed_phase_arrhenius partmclib)
add_executable(test_rxn_emission test/unit_rxn_data/test_rxn_emission.F90)
target_link_libraries(test_rxn_emission partmclib)
add_executable(test_rxn_first_order_loss test/unit_rxn_data/test_rxn_first_order_loss.F90)
target_link_libraries(test_rxn_first_order_loss partmclib)
add_executable(test_rxn_HL_phase_transfer test/unit_rxn_data/test_rxn_HL_phase_transfer.F90)
target_link_libraries(test_rxn_HL_phase_transfer partmclib)
add_executable(test_rxn_photolysis test/unit_rxn_data/test_rxn_photolysis.F90)
target_link_libraries(test_rxn_photolysis partmclib)
add_executable(test_rxn_SIMPOL_phase_transfer test/unit_rxn_data/test_rxn_SIMPOL_phase_transfer.F90)
target_link_libraries(test_rxn_SIMPOL_phase_transfer partmclib)
add_executable(test_rxn_troe test/unit_rxn_data/test_rxn_troe.F90)
target_link_libraries(test_rxn_troe partmclib)
add_executable(test_rxn_wet_deposition test/unit_rxn_data/test_rxn_wet_deposition.F90)
target_link_libraries(test_rxn_wet_deposition partmclib)

######################################################################
# test_sub_model_data

add_executable(test_sub_model_PDFiTE test/unit_sub_model_data/test_sub_model_PDFiTE.F90)
target_link_libraries(test_sub_model_PDFiTE partmclib)
add_executable(test_sub_model_UNIFAC test/unit_sub_model_data/test_sub_model_UNIFAC.F90)
target_link_libraries(test_sub_model_UNIFAC partmclib)
add_executable(test_sub_model_ZSR_aerosol_water
               test/unit_sub_model_data/test_sub_model_ZSR_aerosol_water.c
               test/unit_sub_model_data/test_sub_model_ZSR_aerosol_water.F90)
target_link_libraries(test_sub_model_ZSR_aerosol_water partmclib)

######################################################################
# New unit tests (UNDER DEVELOPMENT)





######################################################################
# test_camp_core

add_executable(unit_test_camp_core test/unit_camp_core/test_camp_core.F90)

target_link_libraries(unit_test_camp_core partmclib)

######################################################################
# test_aero_phase_data

add_executable(unit_test_aero_phase_data test/unit_aero_phase_data/test_aero_phase_data.F90)

target_link_libraries(unit_test_aero_phase_data partmclib)

######################################################################
# test_aero_rep_single_particle

add_executable(unit_test_aero_rep_single_particle
	test/unit_aero_rep_data/test_aero_rep_single_particle.c
	test/unit_aero_rep_data/test_aero_rep_single_particle.F90)

target_link_libraries(unit_test_aero_rep_single_particle partmclib)

######################################################################
# test_aero_rep_modal_binned_mass

add_executable(unit_test_aero_rep_modal_binned_mass
        test/unit_aero_rep_data/test_aero_rep_modal_binned_mass.c
	test/unit_aero_rep_data/test_aero_rep_modal_binned_mass.F90)

target_link_libraries(unit_test_aero_rep_modal_binned_mass partmclib)

######################################################################
# test_chem_mech_solver

add_executable(test_chem_mech_solver test/chemistry/test_chem_mech_solver.F90)

target_link_libraries(test_chem_mech_solver partmclib)

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
  src/util.F90 src/constants.F90 ${C_SORT_SRC})

######################################################################
# test_poisson_sample

add_executable(test_poisson_sample test/rand/test_poisson_sample.F90
  src/util.F90 src/rand.F90 src/constants.F90 src/mpi.F90 ${GSL_SRC}
  ${C_SORT_SRC})

target_link_libraries(test_poisson_sample ${GSL_LIBS})

######################################################################
# test_binomial_sample

add_executable(test_binomial_sample test/rand/test_binomial_sample.F90
  src/util.F90 src/rand.F90 src/constants.F90 src/mpi.F90 ${GSL_SRC}
  ${C_SORT_SRC})

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
endif()
