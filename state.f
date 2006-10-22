! -*- mode: f90; -*-
! Copyright (C) 2005,2006 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

module mod_state
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_state(filename, n_bin, TDV, n_spec, MH, VH, env, time)

    use mod_environ
    
    character, intent(in) :: filename*100   ! input filename
    integer, intent(in) :: n_bin            ! number of bins
    integer, intent(in) :: TDV              ! trailing dimension of VH      
    integer, intent(in) :: n_spec           ! number of species
    integer, intent(out) :: MH(n_bin)       ! number of particles per bin
    real*8, intent(out) :: VH(n_bin,TDV,n_spec) ! particle volumes (m^3)
    type(environ), intent(out) :: env       ! environment state
    real*8, intent(out) :: time             ! current time (s)
    
    integer, parameter :: f_in = 20
    
    character :: dum*100
    integer :: i, j, k, dum_int_1, dum_int_2, dum_int_3
    integer :: n_bin_test, TDV_test, n_spec_test

    open(f_in, file=filename)
    
    read(f_in, '(a20,e20.10)') dum, time
    read(f_in, '(a20,e20.10)') dum, env%T
    read(f_in, '(a20,e20.10)') dum, env%RH
    read(f_in, '(a20,e20.10)') dum, env%V_comp
    read(f_in, '(a20,e20.10)') dum, env%p
    read(f_in, '(a20,i20)') dum, n_bin_test
    read(f_in, '(a20,i20)') dum, TDV_test
    read(f_in, '(a20,i20)') dum, n_spec_test
    
    if (n_bin .ne. n_bin_test) then
       write(0,*) 'ERROR: n_bin mismatch'
       call exit(1)
    end if
    if (TDV .lt. TDV_test) then
       write(0,*) 'ERROR: TDV mismatch: too small'
       call exit(1)
    end if
    if (n_spec .ne. n_spec_test) then
       write(0,*) 'ERROR: n_spec mismatch'
       call exit(1)
    end if
    
    do i = 1,n_bin
       read(f_in,'(i20,i20)') dum_int_1, MH(i)
    end do
    
    do i = 1,n_bin
       do j = 1,MH(i)
          do k = 1,n_spec
             read(f_in,'(i12,i12,i12,e30.20)') &
                  dum_int_1, dum_int_2, dum_int_2, VH(i,j,k)
          end do
       end do
    end do
    
    close(unit=f_in)
    
  end subroutine read_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_state_hybrid(n_bin, TDV, n_spec, MH, VH, env, &
       index, time)
    
    use mod_environ
    
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: TDV          ! trailing dimension of VH      
    integer, intent(in) :: n_spec       ! number of species
    integer, intent(in) :: MH(n_bin)    ! number of particles per bin
    real*8, intent(in) :: VH(n_bin,TDV,n_spec)  ! particle volumes (m^3)
    type(environ), intent(in) :: env    ! environment state
    integer, intent(in) :: index        ! filename index
    real*8, intent(in) :: time          ! current time (s)
    
    integer, parameter :: funit = 31  ! unit for output
    
    character*50 outname
    integer i, j, k
    
    write(outname, '(a6,i4.4,a2)') 'state_', index, '.d'
    open(unit=funit,file=outname)
    write(funit,'(a20,e20.10)') 'time(s)', time
    write(funit,'(a20,e20.10)') 'temp(K)', env%T
    write(funit,'(a20,e20.10)') 'rh(1)', env%RH
    write(funit,'(a20,e20.10)') 'V_comp(m^3)', env%V_comp
    write(funit,'(a20,e20.10)') 'p(Pa)', env%p
    write(funit,'(a20,i20)') 'n_bin', n_bin
    write(funit,'(a20,i20)') 'TDV', TDV
    write(funit,'(a20,i20)') 'n_spec', n_spec
    do i = 1,n_bin
       write(funit,'(i20,i20)') i, MH(i)
    end do
    do i = 1,n_bin
       do j = 1,MH(i)
          do k = 1,n_spec
             write(funit,'(i12,i12,i12,e30.20)') i, j, k, VH(i,j,k)
          end do
       end do
    end do
    close(unit=funit)
    
  end subroutine write_state_hybrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_state
