
PartMC
======

PartMC: Particle-resolved Monte Carlo code for atmospheric aerosol simulation

Version 2.2.0  
Released 2012-02-25

<http://lagrange.mechse.illinois.edu/partmc/>

References:

   * N. Riemer, M. West, R. A. Zaveri, and R. C. Easter (2009),
     Simulating the evolution of soot mixing state with a
     particle-resolved aerosol model, _J. Geophys. Res._ 114(D09202),
     <http://dx.doi.org/10.1029/2008JD011073>.
   * N. Riemer, M. West, R. A. Zaveri, and R. C. Easter (2010),
     Estimating black carbon aging time-scales with a
     particle-resolved aerosol model, _J. Aerosol Sci._ 41(1), 143-158,
     <http://dx.doi.org/10.1016/j.jaerosci.2009.08.009>.
   * R. A. Zaveri, J. C. Barnard, R. C. Easter, N. Riemer, and M. West
     (2010), Particle-resolved simulation of aerosol size,
     composition, mixing state, and the associated optical and cloud
     condensation nuclei activation properties in an evolving urban
     plume, _J. Geophys. Res._ 115(D17210),
     <http://dx.doi.org/10.1029/2009JD013616>.
   * R. E. L. DeVille, N. Riemer, and M. West, Weighted Flow
     Algorithms (WFA) for stochastic particle coagulation,
     _J. Comp. Phys._ 230(23), 8427-8451, 2011,
     <http://dx.doi.org/10.1016/j.jcp.2011.07.027>
   * J. Ching, N. Riemer, and M. West, Impacts of black carbon mixing
     state on black carbon nucleation scavenging: Insights from a
     particle-resolved model, _J. Geophys. Res._ 117(D23209), 2012,
     <http://dx.doi.org/10.1029/2012JD018269>
   * M. D. Michelotti, M. T. Heath, and M. West, Binning for efficient
     stochastic multiscale particle simulations, _Atmos. Chem. Phys._
     11(4), 1071-1096, 2013, <http://dx.doi.org/10.1137/130908038>
   * N. Riemer and M. West, Quantifying aerosol mixing state with
     entropy and diversity measures, _Atmos. Chem. Phys._ 13,
     11423-11439, 2013, <http://dx.doi.org/10.5194/acp-13-11423-2013>
   * J. Tian, N. Riemer, M. West, L. Pfaffenberger, H. Schlager, and
     A. Petzold, Modeling the evolution of aerosol particles in a ship
     plume using PartMC-MOSAIC, _Atmos. Chem. Phys._ 14, 5327-5347,
     2014, <http://dx.doi.org/10.5194/acp-14-5327-2014>

Copyright (C) 2005-2015 Nicole Riemer and Matthew West  
Portions copyright (C) Andreas Bott, Richard Easter, and Jeffrey Curtis  
Licensed under the GNU General Public License version 2 or (at your
option) any later version.  
For details see the file COPYING or
<http://www.gnu.org/licenses/old-licenses/gpl-2.0.html>.


Dependencies
============

Required dependencies:

   * Fortran 2003 compiler - <https://gcc.gnu.org/fortran/> or similar
   * CMake version 2.6.4 or higher - <http://www.cmake.org/>
   * NetCDF version 4.2 or higher -
     <http://www.unidata.ucar.edu/software/netcdf/>

Optional dependencies:

   * MOSAIC chemistry code version 2012-01-25 - Available from Rahul
     Zaveri - <Rahul.Zaveri@pnl.gov>
   * MPI parallel support - <http://www.open-mpi.org/>
   * GSL for random number generators -
     <http://www.gnu.org/software/gsl/>
   * SUNDIALS ODE solver version 2.4 for condensation support -
     <http://www.llnl.gov/casc/sundials/>
   * gnuplot for testcase plotting - <http://www.gnuplot.info/>


Installation
============

1. Install cmake and NetCDF (see above). The NetCDF libraries are
   required to compile PartMC. The `netcdf.mod` Fortran 90 module file
   is required, and it must be produced by the same compiler being
   used to compile PartMC.

2. Unpack PartMC:

        tar xzvf partmc-2.2.0.tar.gz

3. Change into the main PartMC directory (where this README file is
   located):

        cd partmc-2.2.0

4. Make a directory called `build` and change into it:

        mkdir build
        cd build

5. If necessary, set environment variables to indicate the install
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
   you want to use.

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

8. Compile PartMC and test it as follows. Some tests may fail due to
   bad random initial conditions, so re-run the tests a few times to
   see if failures persist.

        make
        make test

9. To run just a single set of tests do something like:

        ctest -R bidisperse   # argument is a regexp for test names

10. To see what make is doing run it like:

        VERBOSE=1 make

11. To run tests with visible output or to make some plots from the
    tests run them as:

        cd test_run/emission
        ./test_emission_1.sh
        ./test_emission_2.sh
        ./test_emission_3.sh            # similarly for other tests
        gnuplot -persist plot_species.gnuplot # etc...

12. To run full scenarios, do, for example:

        cd ../scenarios/1_urban_plume
        ./run.sh


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
