! Copyright (C) 2009 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Process data for test/bidisperse.

program extract_summary

  use netcdf

  integer, parameter :: out_unit = 64
  character(len=*), parameter :: in_filename = "out/bidisperse_mc_0001.nc"
  character(len=*), parameter :: out_filename = "out/bidisperse_mc_data.txt"
  integer, parameter :: i_unit_num = 1
  integer, parameter :: i_unit_mass = 3
  real*8, parameter :: radius_cutoff = 5d-5
  real*8, parameter :: desired_small_init_num_den = 1d9
  real*8, parameter :: desired_large_init_num_den = 1d5
  real*8, parameter :: init_rel_tol = 5d-2

  integer :: ncid
  integer :: dimid_time, dimid_radius, dimid_aero_species, dimid_unit
  integer :: varid_time, varid_radius, varid_aero_species, varid_unit
  integer :: varid_radius_widths, varid_aero
  integer :: n_time, n_radius, n_aero_species, n_unit
  character(len=1000) :: tmp_str, aero_species_names, data_units
  real*8, allocatable :: time(:)
  real*8, allocatable :: radius(:)
  real*8, allocatable :: radius_widths(:)
  real*8, allocatable :: aero(:,:,:,:)
  integer :: xtype, ndims, nAtts
  integer, dimension(nf90_max_var_dims) :: dimids
  integer :: ios, i_time, i_radius, i_spec
  real*8 :: val, large_init_num_den, small_init_num_den
  real*8 :: large_mass_den, small_num_den
  real*8 :: large_rel_error, small_rel_error
  logical :: bad_init

  ! read NetCDF file
  call nc_check(nf90_open(in_filename, NF90_NOWRITE, ncid))

  call nc_check(nf90_inq_dimid(ncid, "time", dimid_time))
  call nc_check(nf90_Inquire_Dimension(ncid, dimid_time, &
       tmp_str, n_time))
  allocate(time(n_time))
  call nc_check(nf90_inq_varid(ncid, "time", varid_time))
  call nc_check(nf90_get_var(ncid, varid_time, time))
  write(*,*) "n_time:", n_time
  write(*,*) "min time:", minval(time)
  write(*,*) "max time:", maxval(time)

  call nc_check(nf90_inq_dimid(ncid, "radius", dimid_radius))
  call nc_check(nf90_Inquire_Dimension(ncid, dimid_radius, &
       tmp_str, n_radius))
  allocate(radius(n_radius))
  call nc_check(nf90_inq_varid(ncid, "radius", varid_radius))
  call nc_check(nf90_get_var(ncid, varid_radius, radius))
  write(*,*) "n_radius:", n_radius
  write(*,*) "min radius:", minval(radius)
  write(*,*) "max radius:", maxval(radius)

  allocate(radius_widths(n_radius))
  call nc_check(nf90_inq_varid(ncid, "radius_widths", varid_radius_widths))
  call nc_check(nf90_get_var(ncid, varid_radius_widths, radius_widths))
  write(*,*) "min radius_width:", minval(radius_widths)
  write(*,*) "max radius_width:", maxval(radius_widths)

  call nc_check(nf90_inq_dimid(ncid, "aero_species", dimid_aero_species))
  call nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_species, &
       tmp_str, n_aero_species))
  call nc_check(nf90_inq_varid(ncid, "aero_species", varid_aero_species))
  call nc_check(nf90_get_att(ncid, varid_aero_species, &
       "names", aero_species_names))
  write(*,*) "n_aero_species:", n_aero_species
  write(*,*) "aero_species_names: ", trim(aero_species_names)

  call nc_check(nf90_inq_dimid(ncid, "unit", dimid_unit))
  call nc_check(nf90_Inquire_Dimension(ncid, dimid_unit, &
       tmp_str, n_unit))
  if (n_unit /= 4) then
     write(*,*) "ERROR: unexpected number of units"
     call exit(1)
  end if
  call nc_check(nf90_inq_varid(ncid, "unit", varid_unit))
  call nc_check(nf90_get_att(ncid, varid_unit, &
       "data_units", data_units))
  write(*,*) "n_unit:", n_unit
  write(*,*) "data_units: ", trim(data_units)

  call nc_check(nf90_inq_varid(ncid, "aero", varid_aero))
  call nc_check(nf90_Inquire_Variable(ncid, varid_aero, tmp_str, &
       xtype, ndims, dimids, nAtts))
  if ((ndims /= 4) &
       .or. (dimids(1) /= dimid_radius) &
       .or. (dimids(2) /= dimid_aero_species) &
       .or. (dimids(3) /= dimid_unit) &
       .or. (dimids(4) /= dimid_time)) then
     write(*,*) "ERROR: unexpected aero dimids"
     call exit(1)
  end if
  allocate(aero(n_radius, n_aero_species, n_unit, n_time))
  call nc_check(nf90_get_var(ncid, varid_aero, aero))

  call nc_check(nf90_close(ncid))

  ! check initial conditions
  small_init_num_den = 0d0
  large_init_num_den = 0d0
  i_time = 1
  do i_radius = 1,n_radius
     do i_spec = 1,n_aero_species
        val = radius_widths(i_radius) &
             * aero(i_radius, i_spec, i_unit_num, i_time)
        if (radius(i_radius) < radius_cutoff) then
           small_init_num_den = small_init_num_den + val
        else
           large_init_num_den = large_init_num_den + val
        end if
     end do
  end do
  small_rel_error = abs((small_init_num_den &
       - desired_small_init_num_den) / desired_small_init_num_den)
  large_rel_error = abs((large_init_num_den &
       - desired_large_init_num_den) / desired_large_init_num_den)
  bad_init = .false.
  if ((small_rel_error > init_rel_tol) &
       .or. (large_rel_error > init_rel_tol)) then
     write(*,*) "WARNING: randomly chosen initial conditions differ", &
          " by too much from the desired values, so results", &
          " may be bad. To fix, simply re-run, or change the rand_init", &
          " value in run_mc.spec until an acceptable run is obtained."
     bad_init = .true.
  end if

  ! output data
  open(unit=out_unit, file=out_filename, iostat=ios)
  if (ios /= 0) then
     write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(out_filename), ' for writing: ', ios
     call exit(1)
  end if
  do i_time = 1,n_time
     small_num_den = 0d0
     large_mass_den = 0d0
     do i_radius = 1,n_radius
        do i_spec = 1,n_aero_species
           if (radius(i_radius) < radius_cutoff) then
              small_num_den = small_num_den &
                   + radius_widths(i_radius) &
                   * aero(i_radius, i_spec, i_unit_num, i_time)
           else
              large_mass_den = large_mass_den &
                   + radius_widths(i_radius) &
                   * aero(i_radius, i_spec, i_unit_mass, i_time)
           end if
        end do
     end do
     write(out_unit,'(e20.10,e20.10,e20.10)') &
          time(i_time), small_num_den, large_mass_den
  end do
  close(out_unit)

  deallocate(time)
  deallocate(radius)
  deallocate(radius_widths)
  deallocate(aero)

  if (bad_init) then
     call exit(1)
  end if

contains

  subroutine nc_check(status)

    !> Status return value.
    integer, intent(in) :: status

    if (status /= NF90_NOERR) then
       write(0,*) nf90_strerror(status)
       call exit(1)
    end if

  end subroutine nc_check

end program extract_summary
