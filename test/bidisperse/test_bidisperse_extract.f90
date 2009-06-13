! Copyright (C) 2009 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Process data for test/bidisperse.

program extract_summary

  use netcdf

  integer, parameter :: out_unit = 64
  character(len=*), parameter :: in_prefix = "out/bidisperse_mc_0001_"
  character(len=*), parameter :: out_filename = "out/bidisperse_mc_data.txt"
  real*8, parameter :: radius_cutoff = 5d-5
  real*8, parameter :: desired_small_init_num_conc = 1d9
  real*8, parameter :: desired_large_init_num_conc = 1d5
  real*8, parameter :: init_rel_tol = 5d-2

  character(len=1000) :: in_filename
  integer :: ncid
  integer :: dimid_aero_species, dimid_aero_particle
  integer :: varid_time, varid_aero_species
  integer :: varid_aero_particle_mass, varid_aero_density
  integer :: varid_aero_comp_vol
  integer :: n_aero_species, n_aero_particle
  character(len=1000) :: tmp_str, aero_species_names
  real*8 :: time
  real*8, allocatable :: aero_particle_mass(:,:)
  real*8, allocatable :: aero_density(:)
  real*8, allocatable :: aero_comp_vol(:)
  real*8, allocatable :: aero_dist(:,:)
  integer :: xtype, ndims, nAtts
  integer, dimension(nf90_max_var_dims) :: dimids
  integer :: ios, i_time, i_spec, i_part, status
  real*8 :: val, large_init_num_conc, small_init_num_conc
  real*8 :: large_mass_conc, small_num_conc
  real*8 :: radius, volume
  integer :: large_number
  real*8 :: large_rel_error, small_rel_error
  logical :: bad_large_num

  ! open output
  open(unit=out_unit, file=out_filename, iostat=ios)
  if (ios /= 0) then
     write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(out_filename), ' for writing: ', ios
     call exit(1)
  end if

  ! process NetCDF files
  i_time = 0
  bad_large_num = .false.
  do while (.true.)
     i_time = i_time + 1
     write(in_filename,'(a,i8.8,a)') trim(in_prefix), i_time, ".nc"
     status = nf90_open(in_filename, NF90_NOWRITE, ncid)
     if (status /= NF90_NOERR) then
        exit
     end if

     ! read time
     call nc_check(nf90_inq_varid(ncid, "time", varid_time))
     call nc_check(nf90_get_var(ncid, varid_time, time))

     ! read aero_species dimension
     call nc_check(nf90_inq_dimid(ncid, "aero_species", dimid_aero_species))
     call nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_species, &
          tmp_str, n_aero_species))
     
     ! read aero_particle dimension
     call nc_check(nf90_inq_dimid(ncid, "aero_particle", dimid_aero_particle))
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
        call exit(1)
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
        call exit(1)
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
        call exit(1)
     end if
     allocate(aero_comp_vol(n_aero_particle))
     call nc_check(nf90_get_var(ncid, varid_aero_comp_vol, &
          aero_comp_vol))
     
     call nc_check(nf90_close(ncid))

     ! compute information
     small_num_conc = 0d0
     large_mass_conc = 0d0
     large_number = 0
     do i_part = 1,n_aero_particle
        volume = sum(aero_particle_mass(i_part,:) / aero_density)
        radius = (volume / (4d0 / 3d0 &
             * 3.14159265358979323846d0))**(1d0/3d0)
        if (radius < radius_cutoff) then
           small_num_conc = small_num_conc + 1d0 / aero_comp_vol(i_part)
        else
           large_number = large_number + 1
           large_mass_conc = large_mass_conc &
                + sum(aero_particle_mass(i_part,:)) / aero_comp_vol(i_part)
        end if
     end do

     deallocate(aero_particle_mass)
     deallocate(aero_density)
     deallocate(aero_comp_vol)

     ! write output
     write(out_unit,'(e20.10,e20.10,e20.10)') &
          time, small_num_conc, large_mass_conc

     ! check we only ever have one large particle
     if (large_number /= 1) then
        bad_large_num = .true.
     end if

  end do

  close(out_unit)

  if (bad_large_num) then
     write(*,*) "WARNING: other than one large particle found,", &
          " so results will be bad."
     call exit(1)
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check return status of NetCDF function calls.
  subroutine nc_check(status)

    !> Status return value.
    integer, intent(in) :: status

    if (status /= NF90_NOERR) then
       write(0,*) nf90_strerror(status)
       call exit(1)
    end if

  end subroutine nc_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program extract_summary
