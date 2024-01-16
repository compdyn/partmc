2.7.0 - 2023-08-11

  * Add support for SUNDIALS 6+ (Alex Hirzel).

  * Add wrapper around stop (Sylwester Arabas).

  * Change timestepping structure in run_part (Jeff Curtis).

  * Update kappa value of NaCl (Jeff Curtis).

  * Remove optional communicator in photolysis.F90 (Jeff Curtis).

  * Use relative tolerance for comparing water properties (Jeff Curtis).

2.6.1 - 2022-02-18

  * Update to support SUNDIALS version 5.8.0 (Sylwester Arabas).

2.6.0 - 2021-11-03

  * Clean up README formatting (Matt West).

  * Fix particle ID initialization after loading of `aero_state`
    (Jeff Curtis).

  * Automatically retry failing tests up to 10 times (Matt West).

  * Update Docker build to use Fedora 33 (Matt West).

  * Add multiple groups for mixing state index calculations (Matt West).

  * Add interface to CAMP chemistry library (Matt Dawson and Jeff Curtis).

  * Fix typo in `do_init_equilibrate` (Sylwester Arabas).

  * Add `scenarios/5_coag_brownian` (Zhonghua Zheng and Jeff Curtis).

2.5.0 - 2018-11-17

  * Shift NetCDF CMake rules to `netcdf.cmake` (Matt West).

  * Add `include,exclude` parameters to `aero_state_diameters()`
    (Matt West).

  * Add Docker build and TravisCI support (Matt West).

  * Fix particle sorting in `test/average` (Matt West).

2.4.0 - 2017-02-14

  * Delete unnecessary scenarios and add documentation for remaining
    ones (Matt West).

  * Allow selection of weighting function in spec file (Jeff Curtis).

  * Entropy calculation helper functions added (Nicole Riemer).

  * Add chamber models for wall and settling losses (Jian Tian).

  * Add fractal particle model treatment, following Naumann [2003]
    (Jian Tian).

  * Composition sampling uses a new mean-projection algorithm to give
    better spread (Matt West).

  * Clean up browian coagulation function (Jeff Curtis).

  * Convert code internals to use allocatable arrays rather than
    pointers (Jeff Curtis).

2.3.0 - 2015-03-29

  * Pressure is now specified as a profile rather than a constant value.

  * Add variable-composition aerosol distributions with a mean and
    standard deviation per species.

  * Each particle source is now resolved with equal number of
    computational particles.

  * Output NetCDF files now store per-particle number concentration in
    `aero_num_conc`, which is the inverse of the previous
    `aero_comp_vol` variable.

  * Use Fortran for post-processing rather than Python using `stats.F90`.
    An example is in `scenarios/1_urban_plume/1_urban_plume_process.F90`.

  * Add example Python plotting scripts to `scenarios/1_urban_plume`.

  * Add processing code to calculate particle composition entropy.

  * Add exact calculations of `aero_particle_crit_diameter()` and
    `aero_particle_crit_rel_humid()`.

  * Fix averaging with dry particles to correctly handle weightings.

  * Add dry deposition using an accelerated stochastic particle loss
    algorithm.

  * Fix accelerated coagulation for negative weightings to only activate
    for small/large events.

  * Fix `allow_doubling`/`allow_halving` with weighted particles.

  * Replace Bessel function code with an LGPL-licensed version.

  * Distributed parallel coagulation is disabled due to lack of
    weighting support.

2.2.1 - 2012-03-30

  * De-duplicate error codes.

  * Correct long name of NetCDF variable `aero_mass_concentration`.

  * Weighted Flow Algorithm (WFA) reference added to paper
    <http://dx.doi.org/10.1016/j.jcp.2011.07.027>.

2.2.0 - 2012-02-25

  * No longer need to specify size bins for particle simulations
    (adaptive binning occurs automatically).

  * Coagulation with negative weightings is accelerated.

  * Added sampled input mode type.

  * Parallel `n_part` is now total computational particles over all
    processors, not per-processor.

  * No longer need to specify weight for particle simulations (a mixed
    number/mass weighting is always used and adapts automatically).

2.1.5 - 2012-01-30

  * Fix bug causing erroneous `n_orig_part` entries.

2.1.4 - 2011-08-15

  * Error if running in parallel with `do_parallel` false.

2.1.3 - 2011-07-02

  * Fix `n_source` bug in `aero_data` (broke `condense` test).

2.1.2 - 2011-06-28

  * Update documentation to include source code.

2.1.1 - 2011-06-22

  * Fix documentation about `restart` and `gas_data`/`aerosol_data`.

2.1.0 - 2011-05-17

  * Include `solar_zenith_angle` in output files when using MOSAIC.

  * Added source-oriented capability for aerosol particles.

2.0.0 - 2011-01-13

  * Parallel implementations based on remote particle access.

  * Added restart capability from NetCDF state files.

  * Added partmc.py library for python-based analysis of PartMC output
    files.

  * Reimplemented water condensation to be much faster and also correct.

  * Nucleation added with parameterization due to Kuang, McMurray, et al.

  * Weighted particles, with full support for coagulation and MOSAIC.

  * Changed to Poisson sampling for coagulation test number.

  * Input and output now use diameter everywhere rather than radius.

  * `urban_plume2` test-case added, as described in the paper
     <http://dx.doi.org/10.1029/2009JD013616>.

1.2.0 - 2009-06-15

  * Output of full per-particle data in NetCDF format.

  * `urban_plume` test-case added, as described in the paper
    <http://dx.doi.org/10.1029/2008JD011073>.

  * Build system switched to `cmake`.

  * Automated test suite added (`make test`).

1.1.0 - 2008-02-17

  * Internal reorganization to use Fortran 90 derived types for the
    data structures.

  * Integration with the MOSAIC gas- and aerosol-chemistry code.

  * Output is in binary NetCDF format.

  * Parallel implementation using 1D mixing.

1.0.0 - 2007-02-26

  * First release, including hierarchical coagulation and full water
    condensation.
