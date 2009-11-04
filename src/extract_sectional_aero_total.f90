! Copyright (C) 2009 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The extract_sectional_aero_total program.

!> Read NetCDF sectional output files and write out the total aerosol
!> number and mass concentrations in text format.
program extract_sectional_aero_total

  use netcdf

  integer, parameter :: dp = kind(0.d0)
  integer, parameter :: out_unit = 64

  character(len=1000) :: in_prefix, in_filename, out_filename
  integer :: ncid
  integer :: dimid_aero_species, dimid_aero_radius
  integer :: varid_time, varid_aero_species
  integer :: varid_aero_radius, varid_aero_radius_widths
  integer :: varid_aero_mass_concentration, varid_aero_number_concentration
  integer :: n_aero_species, n_bin
  character(len=1000) :: tmp_str, aero_species_names
  real(kind=dp) :: time
  real(kind=dp), allocatable :: aero_radius(:)
  real(kind=dp), allocatable :: aero_radius_widths(:)
  real(kind=dp), allocatable :: aero_mass_concentration(:,:)
  real(kind=dp), allocatable :: aero_number_concentration(:)
  integer :: xtype, ndims, nAtts
  integer, dimension(nf90_max_var_dims) :: dimids
  integer :: ios, i_time, i_spec, i_part, status
  integer :: i_bin, n_time
  real(kind=dp) :: num_conc, mass_conc

  ! process commandline arguments
  if (command_argument_count() .ne. 2) then
     write(6,*) 'Usage: extract_sectional_aero_total' &
          // ' <netcdf_state_prefix> <output_filename>'
     stop 2
  endif
  call get_command_argument(1, in_prefix)
  call get_command_argument(2, out_filename)

  ! open output file
  open(unit=out_unit, file=out_filename, status='replace', iostat=ios)
  if (ios /= 0) then
     write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(out_filename), ' for writing: ', ios
     stop 1
  end if

  ! write information
  write(*,*) "Output file array A has:"
  write(*,*) "  A(i, 1) = time(i) (s)"
  write(*,*) "  A(i, 2) = total number concentration at time(i) (#/m^3)"
  write(*,*) "  A(i, 3) = total mass concentration at time(i) (kg/m^3)"

  ! process NetCDF files
  i_time = 0
  n_time = 0
  do while (.true.)
     i_time = i_time + 1
     write(in_filename,'(a,i8.8,a)') trim(in_prefix), i_time, ".nc"
     status = nf90_open(in_filename, NF90_NOWRITE, ncid)
     if (status /= NF90_NOERR) then
        exit
     end if
     n_time = i_time

     ! read time
     call nc_check(nf90_inq_varid(ncid, "time", varid_time))
     call nc_check(nf90_get_var(ncid, varid_time, time))

     ! read aero_species
     call nc_check(nf90_inq_dimid(ncid, "aero_species", dimid_aero_species))
     call nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_species, &
          tmp_str, n_aero_species))
     
     ! read aero_radius dimension
     call nc_check(nf90_inq_dimid(ncid, "aero_radius", dimid_aero_radius))
     call nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_radius, &
          tmp_str, n_bin))

     ! read aero_radius variable
     call nc_check(nf90_inq_varid(ncid, "aero_radius", &
          varid_aero_radius))
     call nc_check(nf90_Inquire_Variable(ncid, varid_aero_radius, &
          tmp_str, xtype, ndims, dimids, nAtts))
     if ((ndims /= 1) &
          .or. (dimids(1) /= dimid_aero_radius)) then
        write(*,*) "ERROR: unexpected aero_radius dimids"
        stop 1
     end if
     allocate(aero_radius(n_bin))
     call nc_check(nf90_get_var(ncid, varid_aero_radius, &
          aero_radius))

     ! read aero_radius_widths variable
     call nc_check(nf90_inq_varid(ncid, "aero_radius_widths", &
          varid_aero_radius_widths))
     call nc_check(nf90_Inquire_Variable(ncid, varid_aero_radius_widths, &
          tmp_str, xtype, ndims, dimids, nAtts))
     if ((ndims /= 1) &
          .or. (dimids(1) /= dimid_aero_radius)) then
        write(*,*) "ERROR: unexpected aero_radius_widths dimids"
        stop 1
     end if
     allocate(aero_radius_widths(n_bin))
     call nc_check(nf90_get_var(ncid, varid_aero_radius_widths, &
          aero_radius_widths))

     ! read aero_number_concentration
     call nc_check(nf90_inq_varid(ncid, "aero_number_concentration", &
          varid_aero_number_concentration))
     call nc_check(nf90_Inquire_Variable(ncid, &
          varid_aero_number_concentration, tmp_str, xtype, ndims, &
          dimids, nAtts))
     if ((ndims /= 1) &
          .or. (dimids(1) /= dimid_aero_radius)) then
        write(*,*) "ERROR: unexpected aero_number_concentration dimids"
        stop 1
     end if
     allocate(aero_number_concentration(n_bin))
     call nc_check(nf90_get_var(ncid, varid_aero_number_concentration, &
          aero_number_concentration))
     
     ! read aero_mass_concentration
     call nc_check(nf90_inq_varid(ncid, "aero_mass_concentration", &
          varid_aero_mass_concentration))
     call nc_check(nf90_Inquire_Variable(ncid, &
          varid_aero_mass_concentration, tmp_str, xtype, ndims, &
          dimids, nAtts))
     if ((ndims /= 2) &
          .or. (dimids(1) /= dimid_aero_radius) &
          .or. (dimids(2) /= dimid_aero_species)) then
        write(*,*) "ERROR: unexpected aero_mass_concentration dimids"
        stop 1
     end if
     allocate(aero_mass_concentration(n_bin, n_aero_species))
     call nc_check(nf90_get_var(ncid, varid_aero_mass_concentration, &
          aero_mass_concentration))
     
     call nc_check(nf90_close(ncid))

     ! compute number and mass concentrations
     num_conc = sum(aero_number_concentration * aero_radius_widths)
     mass_conc = 0d0
     do i_spec = 1,n_aero_species
        mass_conc = mass_conc + sum(aero_mass_concentration(:,i_spec) &
             * aero_radius_widths)
     end do

     deallocate(aero_radius)
     deallocate(aero_radius_widths)
     deallocate(aero_mass_concentration)
     deallocate(aero_number_concentration)

     ! output data
     write(out_unit, '(3e30.15e3)') time, num_conc, mass_conc

  end do

  if (n_time == 0) then
     write(*,'(a,a)') 'ERROR: no input file found matching: ', &
          trim(in_filename)
     stop 1
  end if

  close(out_unit)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check return status of NetCDF function calls.
  subroutine nc_check(status)

    !> Status return value.
    integer, intent(in) :: status

    if (status /= NF90_NOERR) then
       write(0,*) nf90_strerror(status)
       stop 1
    end if

  end subroutine nc_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEFINE_LOCAL_COMMAND_ARGUMENT
  integer function command_argument_count()
    command_argument_count = iargc()
  end function command_argument_count
  subroutine get_command_argument(i, arg)
    integer, intent(in) :: i
    character(len=*), intent(out) :: arg
    call getarg(i, arg)
  end subroutine get_command_argument
#endif
  
end program extract_sectional_aero_total
