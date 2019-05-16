mkdir build
mkdir install
mkdir install/examples
cd build
cmake -D CMAKE_BUILD_TYPE=release MPI_ENABLE:BOOL=TRUE -D KLU_ENABLE:BOOL=TRUE -D CMAKE_C_FLAGS_RELEASE="-g -pg" -D CMAKE_CXX_FLAGS_RELEASE="-g -pg"  -D CMAKE_SHARED_LINKER_FLAGS_RELEASE="-g -pg"  -D CMAKE_EXE_LINKER_FLAGS_RELEASE="-g -pg" -D KLU_LIBRARY_DIR=$SUITE_SPARSE_PARTMC_ROOT/lib -D KLU_INCLUDE_DIR=$SUITE_SPARSE_PARTMC_ROOT/include -D CMAKE_INSTALL_PREFIX=$(pwd)/../install -D EXAMPLES_INSTALL_PATH=$(pwd)/../install/examples ..
make install
cd ..

