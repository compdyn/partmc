FROM fedora:33

RUN dnf -y update \
    && dnf -y install \
        gcc-gfortran \
        make \
        netcdf-fortran-devel \
        sundials-devel \
        gsl-devel \
        cmake \
    && dnf clean all

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
