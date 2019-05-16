# Install script for directory: /tmp/tmp.aSqUIr3xL5

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sundials" TYPE FILE FILES "/tmp/tmp.aSqUIr3xL5/cmake-build-debug/include/sundials/sundials_config.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sundials" TYPE FILE FILES "/tmp/tmp.aSqUIr3xL5/cmake-build-debug/include/sundials/sundials_fconfig.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sundials" TYPE FILE FILES "/tmp/tmp.aSqUIr3xL5/LICENSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/sundials/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/nvec_ser/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/sunmat_dense/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/sunmat_band/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/sunmat_sparse/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/sunlinsol_band/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/sunlinsol_dense/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/sunlinsol_spgmr/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/sunlinsol_spfgmr/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/sunlinsol_spbcgs/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/sunlinsol_sptfqmr/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/sunlinsol_pcg/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/arkode/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/cvode/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/cvodes/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/ida/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/idas/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/src/kinsol/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/arkode/C_serial/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/cvode/serial/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/cvodes/serial/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/ida/serial/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/idas/serial/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/kinsol/serial/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/nvector/serial/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/sunmatrix/dense/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/sunmatrix/band/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/sunmatrix/sparse/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/sunlinsol/band/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/sunlinsol/dense/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/sunlinsol/spgmr/serial/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/sunlinsol/spfgmr/serial/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/sunlinsol/spbcgs/serial/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/sunlinsol/sptfqmr/serial/cmake_install.cmake")
  include("/tmp/tmp.aSqUIr3xL5/cmake-build-debug/examples/sunlinsol/pcg/serial/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/tmp/tmp.aSqUIr3xL5/cmake-build-debug/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
