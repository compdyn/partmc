# Install script for directory: /tmp/tmp.aSqUIr3xL5/examples/sunlinsol/sptfqmr/serial

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
   "/usr/local/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial.c;/usr/local/examples/sunlinsol/sptfqmr/serial/test_sunlinsol.h;/usr/local/examples/sunlinsol/sptfqmr/serial/test_sunlinsol.c;/usr/local/examples/sunlinsol/sptfqmr/serial/sundials_linearsolver.c;/usr/local/examples/sunlinsol/sptfqmr/serial/sundials_nvector.c")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/examples/sunlinsol/sptfqmr/serial" TYPE FILE FILES
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/sptfqmr/serial/../../test_sunlinsol.h"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/sptfqmr/serial/../../test_sunlinsol.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/sptfqmr/serial/../../../../src/sundials/sundials_linearsolver.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/sptfqmr/serial/../../../../src/sundials/sundials_nvector.c"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial.c;/usr/local/examples/sunlinsol/sptfqmr/serial/test_sunlinsol.h;/usr/local/examples/sunlinsol/sptfqmr/serial/test_sunlinsol.c;/usr/local/examples/sunlinsol/sptfqmr/serial/sundials_linearsolver.c;/usr/local/examples/sunlinsol/sptfqmr/serial/sundials_nvector.c")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/examples/sunlinsol/sptfqmr/serial" TYPE FILE FILES
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/sptfqmr/serial/../../test_sunlinsol.h"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/sptfqmr/serial/../../test_sunlinsol.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/sptfqmr/serial/../../../../src/sundials/sundials_linearsolver.c"
    "/tmp/tmp.aSqUIr3xL5/examples/sunlinsol/sptfqmr/serial/../../../../src/sundials/sundials_nvector.c"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/examples/sunlinsol/sptfqmr/serial/CMakeLists.txt")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/examples/sunlinsol/sptfqmr/serial" TYPE FILE FILES "/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/sunlinsol/sptfqmr/serial/CMakeLists.txt")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/examples/sunlinsol/sptfqmr/serial/Makefile")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/examples/sunlinsol/sptfqmr/serial" TYPE FILE RENAME "Makefile" FILES "/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/sunlinsol/sptfqmr/serial/Makefile_ex")
endif()

