! ---

program extract_freezing
    use netcdf
    implicit none
    integer, parameter :: dp = kind(0.d0)
    integer, parameter :: out_unit = 65
    character(len=*), parameter :: in_dir = "out"
    character(len=*), parameter :: in_prefix = "freezing_part_"
    character(len=*), parameter :: out_prefix = "out/freezing_part_data_"
    character(len=1000) :: in_filename
    character(len=1000) :: out_filename
    integer :: ncid
    integer :: varid_time, varid_temperature, dimid_aero_particle
    integer :: varid_aero_num_conc, varid_aero_frozen
    integer :: n_aero_particle
    integer :: xtype, ndims, nAtts
    character(len=1000) :: tmp_str
    real(kind=dp), allocatable :: aero_num_conc(:)
    integer(kind=dp), allocatable :: aero_frozen(:)
    real(kind=dp) :: time, temperature
    real(kind=dp) :: frozen_fraction
    real(kind=dp) :: frozen_fraction_total, frozen_fraction_mean
    real(kind=dp) :: frozen_fraction_max, frozen_fraction_min
    integer, dimension(nf90_max_var_dims) :: dimids
    integer :: i_time, i_ens, status, ios
    character(len=1000) :: caseName

    call get_command_argument(1, caseName)
    write(out_filename, '(a,a,a)') trim(out_prefix), trim(caseName), ".txt"
    ! open output
    open(unit=out_unit, file=out_filename, iostat=ios)
    if (ios /= 0) then
        write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(out_filename), ' for writing: ', ios
        stop 1
    end if

    i_time = 0
    do while(.true.)
        i_time = i_time + 1
        !do i_ens = 1, 10
        i_ens = 0
        frozen_fraction_total = 0d0
        frozen_fraction_max = -9999d0
        frozen_fraction_min = 9999d0
        do while(.true.)
            i_ens = i_ens + 1
            write(in_filename,'(a,a,a,a,a,i4.4, a, i8.8,a)') &
                 trim(in_dir), "/", trim(caseName), "/", trim(in_prefix), i_ens, "_", i_time, ".nc"
            !print*, in_filename
            status = nf90_open(trim(in_filename), NF90_NOWRITE, ncid)
            if (status /= NF90_NOERR) then
                exit
            end if
            !print*, trim(in_filename)

            call nc_check(nf90_inq_varid(ncid, "time", varid_time))
            !call nc_check(nf90_inq_varid(ncid, "time", varid_time))
            call nc_check(nf90_get_var(ncid, varid_time, time))
            call nc_check(nf90_inq_varid(ncid, "temperature", varid_temperature))
            !call nc_check(nf90_inq_varid(ncid, "temperature", varid_temperature))
            call nc_check(nf90_get_var(ncid, varid_temperature, temperature))

            call nc_check(nf90_inq_dimid(ncid, "aero_particle", &
                        dimid_aero_particle))
            call nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_particle, &
                        tmp_str, n_aero_particle))
            call nc_check(nf90_inq_varid(ncid, "aero_num_conc", &
                        varid_aero_num_conc))
            call nc_check(nf90_Inquire_Variable(ncid, varid_aero_num_conc, &
                        tmp_str, xtype, ndims, dimids, nAtts))
            if ((ndims /= 1) &
                        .or. (dimids(1) /= dimid_aero_particle)) then
                write(*,*) "ERROR: unexpected aero_num_conc dimids"
                stop 1
            end if
            allocate(aero_num_conc(n_aero_particle))
            call nc_check(nf90_get_var(ncid, varid_aero_num_conc, &
                        aero_num_conc))

            call nc_check(nf90_inq_varid(ncid, "aero_frozen", &
                        varid_aero_frozen))
            call nc_check(nf90_Inquire_Variable(ncid, varid_aero_frozen, &
                        tmp_str, xtype, ndims, dimids, nAtts))
            if ((ndims /= 1) &
                        .or. (dimids(1) /= dimid_aero_particle)) then
                write(*,*) "ERROR: unexpected aero_frozen dimids"
                stop 1
            end if
            allocate(aero_frozen(n_aero_particle))
            call nc_check(nf90_get_var(ncid, varid_aero_frozen, &
                        aero_frozen))

            call nc_check(nf90_close(ncid))
            frozen_fraction = sum(aero_frozen * aero_num_conc) / & 
                    sum(aero_num_conc)
            frozen_fraction_total = frozen_fraction_total + &
                    frozen_fraction
            if (frozen_fraction .gt. frozen_fraction_max) then
                frozen_fraction_max = frozen_fraction
            end if
            if (frozen_fraction .lt. frozen_fraction_min) then
                frozen_fraction_min = frozen_fraction
            end if
            
            !print*, time, frozen_fraction, frozen_fraction_max, &
            !    frozen_fraction_min
            deallocate(aero_num_conc)
            deallocate(aero_frozen)
        end do
        if (i_ens .eq. 1) then
            exit
        end if
        frozen_fraction_mean = frozen_fraction_total / (i_ens - 1)
        write(out_unit,'(f12.6,f12.6,e20.10,e20.10,e20.10)') time, temperature,&
                frozen_fraction_mean, frozen_fraction_max, frozen_fraction_min
            
    end do

    close(out_unit)
contains
    !> Check return status of NetCDF function calls.
    subroutine nc_check(status)

        !> Status return value.
        integer, intent(in) :: status

        if (status /= NF90_NOERR) then
           write(0,*) nf90_strerror(status)
           stop 1
        end if

    end subroutine nc_check
end program extract_freezing

   

