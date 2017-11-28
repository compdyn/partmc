FROM fedora:27

RUN dnf -y update \
    && dnf -y install \
        gcc-gfortran \
        netcdf-fortran-devel \
        sundials-devel \
        gsl-devel \
        cmake \
    && dnf clean all

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

RUN mkdir /build \
    && cd /build \
    && cmake /partmc \
    && make
