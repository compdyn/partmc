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
    && make install INSTALL=/usr/local BLAS="-L/lib64 -lopenblas"

# Install SUNDIALS with sparse matrix functionality
# RUN curl -LO https://computation.llnl.gov/projects/sundials/download/sundials-3.1.2.tar.gz \
#    && tar -zxvf sundials-3.1.2.tar.gz \
#    && cd sundials-3.1.2 \
    # Increase number of iterations for initial timestep calculation in cvode
    # to allow for RHS function failures due to concentrations < zero at coarse time steps
#    && sed -i s/"MAX_ITERS[ ]*4"/"MAX_ITERS  4000"/g src/cvode/cvode.c \
    # Add a condition to the cvodes error test to disallow negative concentrations
#    && sed -i s/"realtype dsm;"/"realtype dsm, min_conc; N_VLinearSum(cv_mem->cv_l[0], cv_mem->cv_acor, ONE, cv_mem->cv_zn[0], cv_mem->cv_ftemp); min_conc = N_VMin(cv_mem->cv_ftemp); if (min_conc < ZERO \&\& min_conc > -1.0e-30) {N_VAbs(cv_mem->cv_ftemp, cv_mem->cv_ftemp); N_VLinearSum(-cv_mem->cv_l[0], cv_mem->cv_acor, ONE, cv_mem->cv_ftemp, cv_mem->cv_zn[0]); min_conc = ZERO;}"/g src/cvode/cvode.c \
#    && sed -i s/"dsm <= ONE"/"dsm <= ONE \&\& min_conc >= ZERO"/g src/cvode/cvode.c \
#    && mkdir build \
#    && cd build \
#    && cmake -D CMAKE_BUILD_TYPE=release \
#    && cmake -D CMAKE_BUILD_TYPE=debug -D CMAKE_C_FLAGS_DEBUG="-g -pg" -D CMAKE_EXE_LINKER_FLAGS_DEBUG="-pg" -D CMAKE_MODULE_LINKER_FLAGS_DEBUG="-pg" -D CMAKE_SHARED_LINKER_FLAGS_DEBUG="-pg" \
#     -D KLU_ENABLE:BOOL=TRUE -D KLU_LIBRARY_DIR=/usr/local/lib -D KLU_INCLUDE_DIR=/usr/local/include .. \
#    && make install

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

# Install a modified version of CVODE
RUN tar -zxvf /partmc/cvode-1.0-alpha.tar.gz \
    && cd cvode-1.0-alpha \
    && mkdir build \
    && cd build \
    && cmake -D CMAKE_BUILD_TYPE=debug \
             -D CMAKE_C_FLAGS_DEBUG="-g -pg" \
             -D CMAKE_EXE_LINKER_FLAGS_DEBUG="-pg" \
             -D CMAKE_MODULE_LINKER_FLAGS_DEBUG="-pg" \
             -D CMAKE_SHARED_LINKER_FLAGS_DEBUG="-pg" \
             -D KLU_ENABLE:BOOL=TRUE \
             -D KLU_LIBRARY_DIR=/usr/local/lib \
             -D KLU_INCLUDE_DIR=/usr/local/include \
             .. \
    && make install

RUN mkdir build \
    && cd build \
    && export JSON_FORTRAN_HOME="/usr/local/jsonfortran-gnu-6.1.0" \
    && cmake -D CMAKE_BUILD_TYPE=debug -D CMAKE_C_FLAGS_DEBUG="-g -pg" -D CMAKE_Fortran_FLAGS_DEBUG="-g -pg" -D CMAKE_MODULE_LINKER_FLAGS="-pg" \
    -D ENABLE_SUNDIALS:BOOL=TRUE -D SUNDIALS_CVODE_LIB=/usr/local/lib/libsundials_cvode.so -D SUNDIALS_INCLUDE_DIR=/usr/local/include /partmc \
    && make
