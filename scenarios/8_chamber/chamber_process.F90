! Copyright (C) 2009-2013, 2016 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The process program.

!> Read NetCDF output files and process them.
program process

  use pmc_output
  use pmc_stats

  character(len=PMC_MAX_FILENAME_LEN), parameter :: prefix &
       = "out/chamber"

  character(len=PMC_MAX_FILENAME_LEN) :: in_filename, out_filename
  type(bin_grid_t) :: diam_grid
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  type(env_state_t) :: env_state
  integer :: ncid, index, repeat, i_index, i_repeat, n_index, n_repeat
  real(kind=dp) :: time, del_t, tot_num_conc, tot_mass_conc
  character(len=PMC_UUID_LEN) :: uuid
  real(kind=dp), allocatable :: times(:), mobility_diameters(:), &
       num_concs(:), dry_masses(:), masses(:), num_dist(:), mass_dist(:)
  type(stats_1d_t) :: stats_num_dist, stats_mass_dist, stats_tot_num_conc, &
       stats_tot_mass_conc

  call pmc_mpi_init()

  call input_n_files(prefix, n_repeat, n_index)

  call bin_grid_make(diam_grid, BIN_GRID_TYPE_LOG, 180, 1d-9, 1d-3)

  allocate(times(n_index))

  do i_index = 1,n_index
     do i_repeat = 1,n_repeat
        call make_filename(in_filename, prefix, ".nc", i_index, i_repeat)
        write(*,*) "Processing " // trim(in_filename)
        call input_state(in_filename, index, time, del_t, repeat, &
             uuid, aero_data=aero_data, aero_state=aero_state, &
             env_state=env_state)

        times(i_index) = time

        mobility_diameters = aero_state_mobility_diameters(aero_state, &
             aero_data, env_state)
        num_concs = aero_state_num_concs(aero_state, aero_data)
        num_dist = bin_grid_histogram_1d(diam_grid, mobility_diameters, &
             num_concs) * log(10d0)
        call stats_1d_add(stats_num_dist, num_dist)

        tot_num_conc = sum(num_concs)
        call stats_1d_add_entry(stats_tot_num_conc, tot_num_conc, i_index)

        masses = aero_state_masses(aero_state, aero_data)
        mass_dist = bin_grid_histogram_1d(diam_grid, mobility_diameters, &
             masses * num_concs) * log(10d0)
        call stats_1d_add(stats_mass_dist, mass_dist)

        tot_mass_conc = sum(masses * num_concs)
        call stats_1d_add_entry(stats_tot_mass_conc, tot_mass_conc, i_index)
     end do

     call make_filename(out_filename, prefix, "_num_dist.txt", index)
     write(*,*) "Writing " // trim(out_filename)
     call stats_1d_output_text(stats_num_dist, out_filename, diam_grid%centers)
     call stats_1d_clear(stats_num_dist)

     call make_filename(out_filename, prefix, "_mass_dist.txt", index)
     write(*,*) "Writing " // trim(out_filename)
     call stats_1d_output_text(stats_mass_dist, out_filename, &
          diam_grid%centers)
     call stats_1d_clear(stats_mass_dist)
  end do

  call make_filename(out_filename, prefix, "_tot_num_conc.txt")
     write(*,*) "Writing " // trim(out_filename)
  call stats_1d_output_text(stats_tot_num_conc, out_filename, times)
  call stats_1d_clear(stats_tot_num_conc)

  call make_filename(out_filename, prefix, "_tot_mass_conc.txt")
     write(*,*) "Writing " // trim(out_filename)
  call stats_1d_output_text(stats_tot_mass_conc, out_filename, times)
  call stats_1d_clear(stats_tot_mass_conc)

  call pmc_mpi_finalize()

end program process
