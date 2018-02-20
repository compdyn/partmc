FROM fedora:27

RUN dnf -y update \
    && dnf -y install \
        gcc-gfortran \
        gcc-c++ \
        netcdf-fortran-devel \
        gsl-devel \
        metis-devel \
        lapack-devel \
        openblas-devel \
        cmake \
    && dnf clean all

# Build the SuiteSparse libraries for sparse matrix support
RUN curl -LO http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-5.1.0.tar.gz \
    && tar -zxvf SuiteSparse-5.1.0.tar.gz \
    && export CXX=/usr/bin/cc \
    && cd SuiteSparse \
    && make install INSTALL=/usr/local BLAS='-L/lib64 -lopenblas'

# Install SUNDIALS with sparse matrix functionality
RUN curl -LO https://computation.llnl.gov/projects/sundials/download/cvodes-3.1.0.tar.gz \
    && tar -zxvf cvodes-3.1.0.tar.gz \
    && cd cvodes-3.1.0 \
    && mkdir build \
    && cd build \
    && cmake -D KLU_ENABLE:BOOL=TRUE -D KLU_LIBRARY_DIR=/usr/local/lib -D KLU_INCLUDE_DIR=/usr/local/include .. \
    && make install

# Install json-fortran
RUN curl -LO https://github.com/jacobwilliams/json-fortran/archive/6.1.0.tar.gz \
    && tar -zxvf 6.1.0.tar.gz \
    && cd json-fortran-6.1.0 \
    && export FC=gfortran \
    && mkdir build \
    && cd build \
    && cmake -D SKIP_DOC_GEN:BOOL=TRUE .. \
    && make install

# NOTE: Modify .dockerignore to whitelist files/directories to copy.
COPY . /partmc/

RUN mkdir partmc/build \
    && cd partmc/build \
    && cmake -D ENABLE_SUNDIALS:BOOL=TRUE -D SUNDIALS_CVODE_LIB=/usr/local/lib/libsundials_cvodes.so -D SUNDIALS_INCLUDE_DIR=/usr/local/include /partmc \
    && make
