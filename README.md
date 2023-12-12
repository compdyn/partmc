![PartMC](https://raw.githubusercontent.com/wiki/compdyn/partmc/logo.svg)
======

PartMC: Particle-resolved Monte Carlo code for atmospheric aerosol simulation

[![Latest version](https://img.shields.io/github/tag/compdyn/partmc.svg?label=version)](https://github.com/compdyn/partmc/blob/master/ChangeLog.md) [![Docker build status](https://img.shields.io/docker/automated/compdyn/partmc.svg)](https://hub.docker.com/r/compdyn/partmc/builds/) [![Github Actions Status](https://github.com/compdyn/partmc/workflows/CI/badge.svg?branch=master)](https://github.com/compdyn/partmc/actions/workflows/main.yml) [![License](https://img.shields.io/github/license/compdyn/partmc.svg)](https://github.com/compdyn/partmc/blob/master/COPYING) [![DOI](https://zenodo.org/badge/24058992.svg)](https://zenodo.org/badge/latestdoi/24058992)

Version 2.7.0  
Released 2023-08-11

**Source:** <https://github.com/compdyn/partmc>

**Homepage:** <http://lagrange.mechse.illinois.edu/partmc/>

**Cite as:** M. West, N. Riemer, J. Curtis, M. Michelotti, and J. Tian (2023) PartMC, [![version](https://img.shields.io/github/release/compdyn/partmc.svg?label=version)](https://github.com/compdyn/partmc), [![DOI](https://zenodo.org/badge/24058992.svg)](https://zenodo.org/badge/latestdoi/24058992)

Copyright (C) 2005-2022 Nicole Riemer and Matthew West  
Portions copyright (C) Andreas Bott, Richard Easter, Jeffrey Curtis,
Matthew Michelotti, and Jian Tian  
Licensed under the GNU General Public License version 2 or (at your
option) any later version.  
For details see the file COPYING or
<http://www.gnu.org/licenses/old-licenses/gpl-2.0.html>.

**References:**

   * N. Riemer, M. West, R. A. Zaveri, and R. C. Easter (2009)
     Simulating the evolution of soot mixing state with a
     particle-resolved aerosol model, _J. Geophys. Res._ 114(D09202),
     <http://dx.doi.org/10.1029/2008JD011073>.
   * N. Riemer, M. West, R. A. Zaveri, and R. C. Easter (2010)
     Estimating black carbon aging time-scales with a
     particle-resolved aerosol model, _J. Aerosol Sci._ 41(1),
     143-158, <http://dx.doi.org/10.1016/j.jaerosci.2009.08.009>.
   * R. A. Zaveri, J. C. Barnard, R. C. Easter, N. Riemer, and M. West
     (2010) Particle-resolved simulation of aerosol size, composition,
     mixing state, and the associated optical and cloud condensation
     nuclei activation properties in an evolving urban plume,
     _J. Geophys. Res._ 115(D17210),
     <http://dx.doi.org/10.1029/2009JD013616>.
   * R. E. L. DeVille, N. Riemer, and M. West (2011) Weighted Flow
     Algorithms (WFA) for stochastic particle coagulation,
     _J. Comp. Phys._ 230(23), 8427-8451,
     <http://dx.doi.org/10.1016/j.jcp.2011.07.027>
   * J. Ching, N. Riemer, and M. West (2012) Impacts of black carbon
     mixing state on black carbon nucleation scavenging: Insights from
     a particle-resolved model, _J. Geophys. Res._ 117(D23209),
     <http://dx.doi.org/10.1029/2012JD018269>
   * M. D. Michelotti, M. T. Heath, and M. West (2013) Binning for
     efficient stochastic multiscale particle simulations, _Multiscale
     Model. Simul._ 11(4), 1071-1096,
     <http://dx.doi.org/10.1137/130908038>
   * N. Riemer and M. West (2013) Quantifying aerosol mixing state
     with entropy and diversity measures, _Atmos. Chem. Phys._ 13,
     11423-11439, <http://dx.doi.org/10.5194/acp-13-11423-2013>
   * J. Tian, N. Riemer, M. West, L. Pfaffenberger, H. Schlager, and
     A. Petzold (2014) Modeling the evolution of aerosol particles in
     a ship plume using PartMC-MOSAIC, _Atmos. Chem. Phys._ 14,
     5327-5347, <http://dx.doi.org/10.5194/acp-14-5327-2014>
   * R. M. Healy, N. Riemer, J. C. Wenger, M. Murphy, M. West,
     L. Poulain, A. Wiedensohler, I. P. O'Connor, E. McGillicuddy,
     J. R. Sodeau, and G. J. Evans (2014) Single particle diversity
     and mixing state measurements, _Atmos. Chem. and Phys._ 14,
     6289-6299, <http://dx.doi.org/10.5194/acp-14-6289-2014>
   * J. H. Curtis, M. D. Michelotti, N. Riemer, M. Heath, and M. West
     (2016) Accelerated simulation of stochastic particle removal
     processes in particle-resolved aerosol models, _J. Comp. Phys._
     322, 21-32, <http://dx.doi.org/10.1016/j.jcp.2016.06.029>
   * J. Ching, N. Riemer, and M. West (2016) Black carbon mixing state
     impacts on cloud microphysical properties: Effects of aerosol
     plume and environmental conditions, _J. Geophys. Res._ 121(10),
     5990-6013, <http://dx.doi.org/10.1002/2016JD024851>
   * J. Ching, J. Fast, M. West, and N. Riemer (2017) Metrics to
     quantify the importance of mixing state for CCN activity, _Atmos.
     Chem. and Phys._ 17, 7445-7458,
     <http://dx.doi.org/10.5194/acp-17-7445-2017>
   * J. Tian, B. T. Brem, M. West, T. C. Bond, M. J. Rood, and
     N. Riemer (2017) Simulating aerosol chamber experiments with the
     particle-resolved aerosol model PartMC, _Aerosol Sci. Technol._
     51(7), 856-867, <http://dx.doi.org/10.1080/02786826.2017.1311988>
   * J. H. Curtis, N. Riemer, and M. West (2017) A single-column
     particle-resolved model for simulating the vertical distribution
     of aerosol mixing state: WRF-PartMC-MOSAIC-SCM v1.0,
     _Geosci. Model Dev._ 10, 4057-4079,
     <http://dx.doi.org/10.5194/gmd-10-4057-2017>
   * J. Ching, M. West, and N. Riemer (2018) Quantifying impacts of
     aerosol mixing state on nucleation-scavenging of black carbon
     aerosol particles, _Atmosphere_ 9(1), 17,
     <http://dx.doi.org/10.3390/atmos9010017>
   * M. Hughes, J. K. Kodros, J. R. Pierce, M. West, and N. Riemer
     (2018) Machine learning to predict the global distribution of
     aerosol mixing state metrics, _Atmosphere_ 9(1), 15,
     <http://dx.doi.org/10.3390/atmos9010015>
   * R. E. L. DeVille, N. Riemer, and M. West (2019) Convergence of a
     generalized Weighted Flow Algorithm for stochastic particle
     coagulation, _Journal of Computational Dynamics_ 6(1), 69-94,
     <http://dx.doi.org/10.3934/jcd.2019003>
   * N. Riemer, A. P. Ault, M. West, R. L. Craig, and J. H. Curtis
     (2019) Aerosol mixing state: Measurements, modeling, and impacts,
     _Reviews of Geophysics_ 57(2), 187-249,
     <http://dx.doi.org/10.1029/2018RG000615>
   * C. Shou, N. Riemer, T. B. Onasch, A. J. Sedlacek, A. T. Lambe,
     E. R. Lewis, P. Davidovits, and M. West (2019) Mixing state
     evolution of agglomerating particles in an aerosol chamber:
     Comparison of measurements and particle-resolved simulations,
     _Aerosol Science and Technology_ 53(11), 1229-1243,
     <http://dx.doi.org/10.1080/02786826.2019.1661959>
   * J. T. Gasparik, Q. Ye, J. H. Curtis, A. A. Presto, N. M. Donahue,
     R. C. Sullivan, M. West, and N. Riemer (2020) Quantifying errors
     in the aerosol mixing-state index based on limited particle
     sample size, _Aerosol Science and Technology_ 54(12), 1527-1541,
     <http://dx.doi.org/10.1080/02786826.2020.1804523>
   * Z. Zheng, J. H. Curtis, Y. Yao, J. T. Gasparik, V. G. Anantharaj,
     L. Zhao, M. West, and N. Riemer (2021) Estimating submicron
     aerosol mixing state at the global scale with machine learning
     and earth system modeling, _Earth and Space Science_ 8(2),
     e2020EA001500, <http://dx.doi.org/10.1029/2020EA001500>


Running PartMC with Docker
==========================

This is the fastest way to get running.

* **_Step 1:_** Install [Docker Community Edition](https://www.docker.com/community-edition).
    * On Linux and MacOS this is straightforward. [Download from here](https://store.docker.com/search?type=edition&offering=community).
    * On Windows the best version is [Docker Community Edition for Windows](https://store.docker.com/editions/community/docker-ce-desktop-windows), which requires Windows 10 Pro/Edu.

* **_Step 2:_** (Optional) Run the PartMC test suite with:

        docker run -it --rm compdyn/partmc bash -c 'cd /build; make test'

* **_Step 3:_** Run a scenario like the following. This example uses `partmc/scenarios/4_chamber`. This mounts the current directory (`$PWD`, replace with `%cd%` on Windows) into `/run` inside the container, changes into that directory, and then runs PartMC.

        cd partmc/scenarios/4_chamber
        docker run -it --rm -v $PWD:/run compdyn/partmc bash -c 'cd /run; /build/partmc chamber.spec'

In the above `docker run` command the arguments are:

- `-it`: activates "interactive" mode so Ctrl-C works to kill the command
- `--rm`: remove temporary docker container files after running
- `-v LOCAL:REMOTE`: mount the `LOCAL` directory to the `REMOTE` directory inside the container
- `compdyn/partmc`: the docker image to run
- `bash -c 'COMMAND'`: run `COMMAND` inside the docker container

The directory structure inside the docker container is:

    /partmc           # a copy of the partmc git source code repository
    /build            # the diretory in which partmc was compiled
    /build/partmc     # the compiled partmc executable
    /run              # the default diretory to run in


Dependencies
============

Required dependencies:

   * Fortran 2003 compiler - <https://gcc.gnu.org/fortran/> or similar
   * CMake version 2.6.4 or higher - <http://www.cmake.org/>
   * NetCDF version 4.2 or higher -
     <http://www.unidata.ucar.edu/software/netcdf/>

Optional dependencies:

   * CAMP chemistry code - <https://github.com/open-atmos/camp>
   * MOSAIC chemistry code version 2012-01-25 - Available from Rahul
     Zaveri - <Rahul.Zaveri@pnl.gov>
   * MPI parallel support - <http://www.open-mpi.org/>
   * GSL for random number generators -
     <http://www.gnu.org/software/gsl/>
   * SUNDIALS ODE solver for condensation support -
     <http://www.llnl.gov/casc/sundials/>
   * gnuplot for testcase plotting - <http://www.gnuplot.info/>


Installation
============

1. Install cmake and NetCDF (see above). The NetCDF libraries are
   required to compile PartMC. The `netcdf.mod` Fortran 90 module file
   is required, and it must be produced by the same compiler being
   used to compile PartMC.

2. Unpack PartMC:

        tar xzvf partmc-2.7.0.tar.gz

3. Change into the main PartMC directory (where this README file is
   located):

        cd partmc-2.7.0

4. Make a directory called `build` and change into it:

        mkdir build
        cd build

5. If desired, set environment variables to indicate the install
   locations of supporting libraries. If running `echo $SHELL`
   indicates that you are running `bash`, then you can do something
   like:

        export NETCDF_HOME=/
        export MOSAIC_HOME=${HOME}/mosaic-2012-01-25
        export SUNDIALS_HOME=${HOME}/opt
        export GSL_HOME=${HOME}/opt

   Of course the exact directories will depend on where the libraries
   are installed. You only need to set variables for libraries
   installed in non-default locations, and only for those libraries
   you want to use. Everything except NetCDF is optional.

   If `echo $SHELL` instead is `tcsh` or similar, then the environment
   variables can be set like `setenv NETCDF_HOME /` and similarly.

6. Run cmake with the main PartMC directory as an argument (note the
   double-c):

        ccmake ..

7. Inside ccmake press `c` to configure, edit the values as needed,
   press `c` again, then `g` to generate. Optional libraries can be
   activated by setting the `ENABLE` variable to `ON`. For a parallel
   build, toggle advanced mode with `t` and set the
   `CMAKE_Fortran_COMPILER` to `mpif90`, then reconfigure.

8. Optionally, enable compiler warnings by pressing `t` inside ccmake
   to enable advanced options and then setting `CMAKE_Fortran_FLAGS`
   to:

        -O2 -g -fimplicit-none -W -Wall -Wconversion -Wunderflow -Wimplicit-interface -Wno-compare-reals -Wno-unused -Wno-unused-parameter -Wno-unused-dummy-argument -fbounds-check

8. Compile PartMC and test it as follows.

        make
        make test

9. To run just a single test do something like:

        ctest -R bidisperse   # argument is a regexp for test names

10. To see what make is doing run it like:

        VERBOSE=1 make

11. To run tests with visible output or to make some plots from the
    tests run them as follows. Note that tests often rely on earlier
    tests in the same directory, so always run `test_1`, then
    `test_2`, etc. Tests occasionally fail due to random sampling, so
    re-run the entire sequence after failures. For example:

        cd test_run/emission
        ./test_emission_1.sh
        ./test_emission_2.sh
        ./test_emission_3.sh            # similarly for other tests
        gnuplot -persist plot_species.gnuplot # etc...

12. To run full scenarios, do, for example:

        cd ../scenarios/1_urban_plume
        ./1_run.sh


Usage
=====

The main `partmc` command reads `.spec` files and does the run
specified therein. Either particle-resolved runs, sectional-code runs,
or exact solutions can be generated. A run produces one NetCDF file
per output timestep, containing per-particle data (from
particle-resolved runs) or binned data (from sectional or exact
runs). The `extract_*` programs can read these per-timestep NetCDF
files and output ASCII data (the `extract_sectional_*` programs are
used for sectional and exact model output).

Python bindings
===============

The [PyPartMC](https://github.com/open-atmos/PyPartMC) project offers
pip-installable Python bindings to PartMC. Both source and binary
packages are available and ship with all PartMC dependencies included.
PyPartMC exposes internal components of PartMC (utility routines and 
derived types) which then can serve as building blocks to develop PartMC 
simulations in Python. Time stepping can be performed either using the
internal PartMC time-stepper or externally within a Python loop. The
latter allows to couple the simulation with external Python components
in each timestep. PyPartMC features examples developed as Jupyter notebooks.
Snippets of code provided in the README file depict how to use PyPartMC 
from Julia (using PyCall.jl) and Matlab (using Matlab's built-in Python bridge).
