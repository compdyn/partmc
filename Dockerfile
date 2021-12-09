FROM fedora:33

RUN dnf -y update \
    && dnf -y install \
        gcc-gfortran \
        gcc-c++ \
        make \
        netcdf-fortran-devel \
        sundials-devel \
        openblas-devel \
        metis-devel \
        lapack-devel \
        gsl-devel \
        cmake \
    && dnf clean all

# Build the SuiteSparse libraries for sparse matrix support
# (-k included because of problem with SuiteSparse security certificate - 1 Aug 2021)
RUN curl -kLO http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-5.1.0.tar.gz \
    && tar -zxvf SuiteSparse-5.1.0.tar.gz \
    && export CXX=/usr/bin/cc \
    && cd SuiteSparse \
    && make install INSTALL=/usr/local BLAS="-L/lib64 -lopenblas"

# NOTE: Modify .dockerignore to whitelist files/directories to copy.
COPY . /partmc/

RUN mkdir /build \
    && cd /build \
    && cmake -D CMAKE_BUILD_TYPE=release \
             -D CMAKE_C_FLAGS_RELEASE="-O2 -g -Werror -Wall -Wextra" \
             -D CMAKE_Fortran_FLAGS_RELEASE="-O2 -g -Werror -fimplicit-none -Wall -Wextra -Wconversion -Wunderflow -Wimplicit-interface -Wno-compare-reals -Wno-unused -Wno-unused-parameter -Wno-unused-dummy-argument -fbounds-check" \
             -D ENABLE_GSL:BOOL=TRUE \
             -D ENABLE_SUNDIALS:BOOL=TRUE \
             /partmc \
    && make
