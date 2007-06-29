! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Save and restore the state of particle-resolved Monte Carlo
! runs. The state file should contain enough data to restart the
! simulation at the point it was written.
!
! Because it contains the full state of every particle, this is also
! the best way to gain complete access to all statistics of the
! simulation.

module mod_state
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_state_header(state_unit, filename, n_bin, n_spec)

    ! Open the state file, read only the header and then close it
    ! again.
    
    use mod_environ
    use mod_aero_state
    use mod_util
       
    integer, intent(in) :: state_unit   ! unit number to use for state file
    character(len=*), intent(in) :: filename ! input filename
    integer, intent(out) :: n_bin       ! number of bins
    integer, intent(out) :: n_spec      ! number of species
       
    character :: dum*1000

    call open_existing(state_unit, filename)

    read(state_unit, '(a20,i20)') dum, n_bin
    read(state_unit, '(a20,i20)') dum, n_spec

    close(unit=state_unit)

  end subroutine read_state_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_state_bins(state_unit, filename, n_bin, bin_v, dlnr)

    ! Open the state file, read only the bin_v array out and then
    ! close the file again.
    
    use mod_environ
    use mod_aero_state
    use mod_util
       
    integer, intent(in) :: state_unit   ! unit number to use for state file
    character(len=*), intent(in) :: filename ! input filename
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(out) :: bin_v(n_bin) ! volume of particles in bins (m^3)
    real*8, intent(out) :: dlnr         ! bin scale factor
       
    character :: dum*1000
    integer :: i, n_bin_test, dum_int
    real*8 :: dum_real

    call open_existing(state_unit, filename)

    read(state_unit, '(a20,i20)') dum, n_bin_test
    read(state_unit, '(a20,i20)') dum, dum_int
    read(state_unit, '(a20,e20.10)') dum, dum_real
    read(state_unit, '(a20,e20.10)') dum, dum_real
    read(state_unit, '(a20,e20.10)') dum, dum_real
    read(state_unit, '(a20,e20.10)') dum, dum_real
    read(state_unit, '(a20,e20.10)') dum, dum_real
    read(state_unit, '(a20,e20.10)') dum, dlnr

    if (n_bin_test /= n_bin) then
       write(0,*) 'ERROR: n_bin mismatch when reading state'
       call exit(1)
    end if
    
    read(state_unit,'(a)') dum
    do i = 1,n_bin
       read(state_unit,'(i20,e30.20)') dum_int, bin_v(i)
    end do

    close(unit=state_unit)

  end subroutine read_state_bins
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_state(state_unit, filename, n_bin, n_spec, &
       MH, VH, env, time)

    ! Open and read the entire state file (including the header) and
    ! then close it.

    use mod_environ
    use mod_aero_state
    use mod_util
    
    integer, intent(in) :: state_unit   ! unit number to use for state file
    character(len=*), intent(in) :: filename ! input filename
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    integer, intent(out) :: MH(n_bin)   ! number of particles per bin
    type(bin_p), intent(out) :: VH(n_bin) ! particle volumes (m^3)
    type(environ), intent(out) :: env   ! environment state
    real*8, intent(out) :: time         ! current time (s)
    
    character :: dum*1000
    integer :: i, j, k, dum_int_1, dum_int_2, dum_int_3
    integer :: n_bin_test, n_spec_test
    real*8 :: dum_real

    call open_existing(state_unit, filename)
    
    read(state_unit, '(a20,i20)') dum, n_bin_test
    read(state_unit, '(a20,i20)') dum, n_spec_test
    read(state_unit, '(a20,e20.10)') dum, time
    read(state_unit, '(a20,e20.10)') dum, env%T
    read(state_unit, '(a20,e20.10)') dum, env%RH
    read(state_unit, '(a20,e20.10)') dum, env%V_comp
    read(state_unit, '(a20,e20.10)') dum, env%p
    read(state_unit, '(a20,e20.10)') dum, dum_real

    if (n_bin_test /= n_bin) then
       write(0,*) 'ERROR: n_bin mismatch when reading state'
       call exit(1)
    end if
    if (n_spec_test /= n_spec) then
       write(0,*) 'ERROR: n_spec mismatch when reading state'
       call exit(1)
    end if

    read(state_unit,'(a)') dum
    do i = 1,n_bin
       read(state_unit,'(i20,e30.20)') dum_int_1, dum_real
    end do
    
    call zero_array(n_spec, MH, VH)    

    read(state_unit,'(a)') dum
    do i = 1,n_bin
       read(state_unit,'(i20,i20)') dum_int_1, MH(i)
    end do
    
    read(state_unit,'(a)') dum
    do i = 1,n_bin
       do j = 1,MH(i)
          do k = 1,n_spec
             if (j > size(VH(i)%p,1)) then
                call enlarge_bin(VH(i))
             end if
             read(state_unit,'(i12,i12,i12,e30.20)') &
                  dum_int_1, dum_int_2, dum_int_3, VH(i)%p(j,k)
          end do
       end do
    end do
    
    close(unit=state_unit)
    
  end subroutine read_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_state(state_unit, state_prefix, n_bin, n_spec, MH, &
       VH, bin_v, dlnr, env, index, time, i_loop)

    ! Write the current state.
    
    use mod_environ
    use mod_aero_state
    
    integer, intent(in) :: state_unit   ! unit number to use for state file
    character(len=*), intent(in) :: state_prefix ! prefix of state file
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    real*8, intent(in) :: dlnr          ! bin scale factor
    integer, intent(in) :: MH(n_bin)    ! number of particles per bin
    type(bin_p), intent(in) :: VH(n_bin) ! particle volumes (m^3)
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    type(environ), intent(in) :: env    ! environment state
    integer, intent(in) :: index        ! filename index
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: i_loop       ! current loop number
    
    character*300 filename
    integer i, j, k
    
    write(filename, '(a,a,i4.4,a,i8.8,a)') trim(state_prefix), &
         '_', i_loop, '_', index, '.d'
    open(unit=state_unit, file=filename)
    write(state_unit,'(a20,i20)') 'n_bin', n_bin
    write(state_unit,'(a20,i20)') 'n_spec', n_spec
    write(state_unit,'(a20,e20.10)') 'time(s)', time
    write(state_unit,'(a20,e20.10)') 'temp(K)', env%T
    write(state_unit,'(a20,e20.10)') 'RH(1)', env%RH
    write(state_unit,'(a20,e20.10)') 'V_comp(m^3)', env%V_comp
    write(state_unit,'(a20,e20.10)') 'p(Pa)', env%p
    write(state_unit,'(a20,e20.10)') 'dlnr', dlnr

    write(state_unit,'(a1,a19,a30)') '#', 'bin_num', 'vol of bin (m^3)'
    do i = 1,n_bin 
       write(state_unit,'(i20,e30.20)') i, bin_v(i)
    end do  
    write(state_unit,'(a1,a19,a20)') '#', 'bin_num', 'num in bin (#)'
    do i = 1,n_bin
       write(state_unit,'(i20,i20)') i, MH(i)
    end do
    write(state_unit,'(a1,a11,a12,a12,a30)') '#', 'bin_num', 'index', &
         'species', 'volume (m^3)'
    do i = 1,n_bin
       do j = 1,MH(i)
          do k = 1,n_spec
             write(state_unit,'(i12,i12,i12,e30.20)') i, j, k, VH(i)%p(j,k)
          end do
       end do
    end do
    close(unit=state_unit)
    
  end subroutine write_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_state
