FROM fedora:30

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
    && cmake /partmc \
    && make
