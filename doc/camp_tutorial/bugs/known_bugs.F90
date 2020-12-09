! Copyright (C) 2020 Matt Dawson and Christian Guzman
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \page known_bugs Known bugs
!!
!! This page reports some knowned bugs that can affect the program behaviour
!! under some conditions. The user should check this page to understand
!! the impact of these bugs and acts consequently until the bugs are fixed under
!! all the possible conditions.
!! - \ref known_bugs_0

! ***********************************************************************
! ***********************************************************************
! ***********************************************************************

!> \page known_bugs_0 Accuracy divergence under different compilers
!! <a> This bug was found during the merge of the GPU branch develop-128-monarch-multicells
!! into develop_129_merge_cpu_gpu (Github commit: 03d7efe1d345753d22db582cb953f6308d390b35) </a>
!!
!! The results obtained on CAMP and EBI differs from using ICC compiler and GCC.
!! Common configuration:
!!
!! - Test: Test_cb05, 10 Timesteps, tolerance: E-04
!! - CMAKE compiler options: CMAKE_BUILD_TYPE=release ENABLE_MPI=OFF ENABLE_GPU=OFF ENABLE_JSON=ON ENABLE_SUNDIALS=ON
!! - Same "src" and "test" code on both sides
!!
!! Configuration for ICC:
!!
!! - Arquitecture: Marenostrum4
!! - Fortran compiler: Intel 17.0.4.20170411
!! - C compiler: GNU 4.8.5
!! - ICC compiler flags: -extend_source -warn truncated_source
!!
!! Configuration for GCC:
!!
!! - Arquitecture: POWER9
!! - Fortran compiler: GNU 6.4.0
!! - C compiler: GNU 6.4.0
!! - GCC compiler flags: -ffree-line-length-none
!!
!! Sample results from file CAMP_EBI_RESULTS.txt:
!!
!! \image html icc_gcc_results_test_cb05.jpg
!!
!! We can appreciate how the results differ in some species, especially in the case of NO. However,
!! the results above the level of tolerance configured are exactly equal. <b> We recommend increasing
!! the level of tolerance to reduce the effects of this event. </b>
!!
!! The results also varies on the results obtained by the EBI solver. Below we can appreciate
!! multiple warnings from EBI solver that only appears on ICC case.
!!
!! \image html icc_test_cb05_warnings_ebi.jpg
!!
!! On simpler tests like mock_monarch_2 with 4 gas species we saw a very similar behaviour on overall.
!! Only some small differences appears after during the 180 time-steps of the simulation.
!! We can deduce the impact of the compiler depends also on the stiffness of the solver.


