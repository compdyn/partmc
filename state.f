! -*- mode: f90;-*-

module mod_state
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_state(n_bin, TDV, n_spec, MH, VH, env, time, basename)

    use mod_environ
    
    integer, intent(in) :: n_bin            ! number of bins
    integer, intent(in) :: TDV              ! trailing dimension of VH      
    integer, intent(in) :: n_spec           ! number of species
    integer, intent(out) :: MH(n_bin)       ! number of particles per bin
    real*8, intent(out) :: VH(n_bin,TDV,n_spec) ! particle volumes (m^3)
    type(environ), intent(out) :: env       ! environment state
    real*8, intent(out) :: time             ! current time (s)
    character, intent(out) :: basename*100   ! basename of the input filename
    
    integer, parameter :: f_in = 20
    
    character :: dum*100
    integer :: i, j, k, dum_int_1, dum_int_2, dum_int_3
    integer :: n_bin_test, TDV_test, n_spec_test
    
    ! check there is exactly one commandline argument
    if (iargc() .ne. 1) then
       write(0,*) 'Usage: process_out <filename.d>'
       call exit(1)
    end if
    
    ! get and check first commandline argument (must be "filename.d")
    call getarg(1, basename)
    i = len_trim(basename)
    if (i .gt. 40) then
       write(0,*) 'ERROR: filename too long'
       call exit(1)
    end if
    if ((basename(i:i) .ne. 'd') .or. &
         (basename((i-1):(i-1)) .ne. '.')) then
       write(0,*) 'ERROR: Filename must end in .d'
       call exit(1)
    end if
    
    open(f_in, file=basename)
    
    ! chop .d off the end of the filename to get the basename
    basename((i-1):i) = '  '
    
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
    if (TDV .ne. TDV_test) then
       write(0,*) 'ERROR: TDV mismatch'
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

end module mod_state
