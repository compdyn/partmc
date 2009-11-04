! Copyright (C) 2009 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The extract_aero_particle_mass program.

!> Read NetCDF output files and write out the individual particle
!> masses.
program extract_aero_particle_mass

  use netcdf

  integer, parameter :: out_unit = 64
  integer, parameter :: dp = kind(0.d0)

  character(len=1000) :: in_filename, out_filename
  integer :: ncid
  integer :: dimid_aero_species, dimid_aero_particle
  integer :: varid_time, varid_aero_species
  integer :: varid_aero_particle_mass, varid_aero_density
  integer :: varid_aero_comp_vol, varid_aero_id
  integer :: n_aero_species, n_aero_particle
  character(len=1000) :: tmp_str, aero_species_names
  character(len=1000) :: remaining_species
  real(kind=dp) :: time
  real(kind=dp), allocatable :: aero_particle_mass(:,:)
  real(kind=dp), allocatable :: aero_density(:)
  real(kind=dp), allocatable :: aero_comp_vol(:)
  integer, allocatable :: aero_id(:)
  integer :: xtype, ndims, nAtts
  integer, dimension(nf90_max_var_dims) :: dimids
  integer :: ios, i_time, i_spec, i_part, status
  integer :: n_bin, i_bin, n_time, i
  real(kind=dp) :: r_min, r_max, radius, volume, dlnr

  ! process commandline arguments
  if (command_argument_count() .ne. 2) then
     write(6,*) 'Usage: extract_aero_particle_mass ' &
          // '<netcdf_state_file> <output_filename>'
     stop 2
  endif
  call get_command_argument(1, in_filename)
  call get_command_argument(2, out_filename)

  ! open NetCDF file
  call nc_check(nf90_open(in_filename, NF90_NOWRITE, ncid))
  
  ! read time
  call nc_check(nf90_inq_varid(ncid, "time", varid_time))
  call nc_check(nf90_get_var(ncid, varid_time, time))
  
  ! read aero_species
  call nc_check(nf90_inq_dimid(ncid, "aero_species", dimid_aero_species))
  call nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_species, &
       tmp_str, n_aero_species))
  call nc_check(nf90_inq_varid(ncid, "aero_species", varid_aero_species))
  call nc_check(nf90_get_att(ncid, varid_aero_species, &
       "names", aero_species_names))
  
  ! read aero_particle dimension
  status = nf90_inq_dimid(ncid, "aero_particle", dimid_aero_particle)
  if (status == NF90_EBADDIM) then
     ! dimension missing ==> no particles, so skip this time
     write(0,*) 'ERROR: no particles found'
     stop 1
  end if
  call nc_check(status)
  call nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_particle, &
       tmp_str, n_aero_particle))
  
  ! read aero_particle_mass
  call nc_check(nf90_inq_varid(ncid, "aero_particle_mass", &
       varid_aero_particle_mass))
  call nc_check(nf90_Inquire_Variable(ncid, varid_aero_particle_mass, &
       tmp_str, xtype, ndims, dimids, nAtts))
  if ((ndims /= 2) &
       .or. (dimids(1) /= dimid_aero_particle) &
       .or. (dimids(2) /= dimid_aero_species)) then
     write(*,*) "ERROR: unexpected aero_particle_mass dimids"
     stop 1
  end if
  allocate(aero_particle_mass(n_aero_particle, n_aero_species))
  call nc_check(nf90_get_var(ncid, varid_aero_particle_mass, &
       aero_particle_mass))
  
  ! read aero_density
  call nc_check(nf90_inq_varid(ncid, "aero_density", &
       varid_aero_density))
  call nc_check(nf90_Inquire_Variable(ncid, varid_aero_density, &
       tmp_str, xtype, ndims, dimids, nAtts))
  if ((ndims /= 1) &
       .or. (dimids(1) /= dimid_aero_species)) then
     write(*,*) "ERROR: unexpected aero_density dimids"
     stop 1
  end if
  allocate(aero_density(n_aero_species))
  call nc_check(nf90_get_var(ncid, varid_aero_density, &
       aero_density))
  
  ! read aero_comp_vol
  call nc_check(nf90_inq_varid(ncid, "aero_comp_vol", &
       varid_aero_comp_vol))
  call nc_check(nf90_Inquire_Variable(ncid, varid_aero_comp_vol, &
       tmp_str, xtype, ndims, dimids, nAtts))
  if ((ndims /= 1) &
       .or. (dimids(1) /= dimid_aero_particle)) then
     write(*,*) "ERROR: unexpected aero_comp_vol dimids"
     stop 1
  end if
  allocate(aero_comp_vol(n_aero_particle))
  call nc_check(nf90_get_var(ncid, varid_aero_comp_vol, &
       aero_comp_vol))
  
  ! read aero_id
  call nc_check(nf90_inq_varid(ncid, "aero_id", &
       varid_aero_id))
  call nc_check(nf90_Inquire_Variable(ncid, varid_aero_id, &
       tmp_str, xtype, ndims, dimids, nAtts))
  if ((ndims /= 1) &
       .or. (dimids(1) /= dimid_aero_particle)) then
     write(*,*) "ERROR: unexpected aero_id dimids"
     stop 1
  end if
  allocate(aero_id(n_aero_particle))
  call nc_check(nf90_get_var(ncid, varid_aero_id, &
       aero_id))
  
  call nc_check(nf90_close(ncid))
  
  ! write information
  write(*,'(a,a)') "Output file: ", trim(out_filename)
  write(*,'(a,e10.3,a)') "  Output data is for time = ", time, " (s)"
  write(*,'(a)') "  Each row of output is one particle."
  write(*,'(a)') "  The columns of output are:"
  write(*,'(a)') "    column  1: particle ID number"
  write(*,'(a)') "    column  2: computational volume (m^3)"
  write(*,'(a)') "    column  3: particle radius (m)"
  write(*,'(a)') "    column  4: particle total mass (kg)"
  remaining_species = aero_species_names
  do i_spec = 1,n_aero_species
     if (i_spec < n_aero_species) then
        i = index(remaining_species, ',')
     else
        i = index(remaining_species, ' ')
     end if
     write(*,'(a,i2,a,a,a,e10.4,a)') '    column ', i_spec + 4, &
          ': ', remaining_species(1:(i-1)), &
          ' mass (kg) - density = ', aero_density(i_spec), ' (kg/m^3)'
     remaining_species = remaining_species((i+1):)
  end do

  ! open output file
  open(unit=out_unit, file=out_filename, status='replace', iostat=ios)
  if (ios /= 0) then
     write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(out_filename), ' for writing: ', ios
     stop 1
  end if

  ! output data
  do i_part = 1,n_aero_particle
     write(out_unit, '(i15)', advance='no') aero_id(i_part)
     write(out_unit, '(e30.15e3)', advance='no') &
          aero_comp_vol(i_part)
     volume = sum(aero_particle_mass(i_part,:) / aero_density)
     radius = (volume / (4d0 / 3d0 &
          * 3.14159265358979323846d0))**(1d0/3d0)
     write(out_unit, '(e30.15e3)', advance='no') radius
     write(out_unit, '(e30.15e3)', advance='no') &
          sum(aero_particle_mass(i_part,:))
     do i_spec = 1,n_aero_species
        write(out_unit, '(e30.15e3)', advance='no') &
             aero_particle_mass(i_part, i_spec)
     end do
     write(out_unit, *) ''
  end do

  close(out_unit)

  deallocate(aero_particle_mass)
  deallocate(aero_density)
  deallocate(aero_comp_vol)

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
  
end program extract_aero_particle_mass
