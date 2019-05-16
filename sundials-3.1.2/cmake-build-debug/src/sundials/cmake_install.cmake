# Install script for directory: /tmp/tmp.aSqUIr3xL5/src/sundials

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  MESSAGE("
Install shared components
")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sundials" TYPE FILE FILES
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_band.h"
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_dense.h"
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_direct.h"
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_fnvector.h"
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_iterative.h"
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_linearsolver.h"
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_math.h"
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_matrix.h"
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_nvector.h"
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_pcg.h"
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_sparse.h"
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_spbcgs.h"
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_spfgmr.h"
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_spgmr.h"
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_sptfqmr.h"
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_types.h"
    "/tmp/tmp.aSqUIr3xL5/include/sundials/sundials_version.h"
    )
endif()

