# Install script for directory: /tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense

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
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/examples/sunlinsol/dense/test_sunlinsol_dense.c;/usr/local/examples/sunlinsol/dense/test_sunlinsol.h;/usr/local/examples/sunlinsol/dense/test_sunlinsol.c;/usr/local/examples/sunlinsol/dense/sundials_matrix.c;/usr/local/examples/sunlinsol/dense/sundials_linearsolver.c;/usr/local/examples/sunlinsol/dense/sundials_nvector.c")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/examples/sunlinsol/dense" TYPE FILE FILES
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/test_sunlinsol_dense.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../test_sunlinsol.h"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../test_sunlinsol.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../../../src/sundials/sundials_matrix.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../../../src/sundials/sundials_linearsolver.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../../../src/sundials/sundials_nvector.c"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/examples/sunlinsol/dense/test_sunlinsol_dense.c;/usr/local/examples/sunlinsol/dense/test_sunlinsol.h;/usr/local/examples/sunlinsol/dense/test_sunlinsol.c;/usr/local/examples/sunlinsol/dense/sundials_matrix.c;/usr/local/examples/sunlinsol/dense/sundials_linearsolver.c;/usr/local/examples/sunlinsol/dense/sundials_nvector.c")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/examples/sunlinsol/dense" TYPE FILE FILES
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/test_sunlinsol_dense.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../test_sunlinsol.h"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../test_sunlinsol.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../../../src/sundials/sundials_matrix.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../../../src/sundials/sundials_linearsolver.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../../../src/sundials/sundials_nvector.c"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/examples/sunlinsol/dense/test_sunlinsol_dense.c;/usr/local/examples/sunlinsol/dense/test_sunlinsol.h;/usr/local/examples/sunlinsol/dense/test_sunlinsol.c;/usr/local/examples/sunlinsol/dense/sundials_matrix.c;/usr/local/examples/sunlinsol/dense/sundials_linearsolver.c;/usr/local/examples/sunlinsol/dense/sundials_nvector.c")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/examples/sunlinsol/dense" TYPE FILE FILES
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/test_sunlinsol_dense.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../test_sunlinsol.h"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../test_sunlinsol.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../../../src/sundials/sundials_matrix.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../../../src/sundials/sundials_linearsolver.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../../../src/sundials/sundials_nvector.c"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/examples/sunlinsol/dense/test_sunlinsol_dense.c;/usr/local/examples/sunlinsol/dense/test_sunlinsol.h;/usr/local/examples/sunlinsol/dense/test_sunlinsol.c;/usr/local/examples/sunlinsol/dense/sundials_matrix.c;/usr/local/examples/sunlinsol/dense/sundials_linearsolver.c;/usr/local/examples/sunlinsol/dense/sundials_nvector.c")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/examples/sunlinsol/dense" TYPE FILE FILES
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/test_sunlinsol_dense.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../test_sunlinsol.h"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../test_sunlinsol.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../../../src/sundials/sundials_matrix.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../../../src/sundials/sundials_linearsolver.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/dense/../../../src/sundials/sundials_nvector.c"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/examples/sunlinsol/dense/CMakeLists.txt")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/examples/sunlinsol/dense" TYPE FILE FILES "/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/sunlinsol/dense/CMakeLists.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/examples/sunlinsol/dense/Makefile")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/examples/sunlinsol/dense" TYPE FILE RENAME "Makefile" FILES "/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/sunlinsol/dense/Makefile_ex")
endif()

