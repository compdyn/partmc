! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Process the saved state files to obtain summary data. This is pretty
! much a hack at the moment, as lots of the processing routines are
! specific to particles iwth 3 species, of which one is water and the
! other two are of interest.

program process_state

  use mod_bin
  use mod_environ
  use mod_aero_data
  use mod_aero_state
  use mod_output_state

  integer :: n_bin                      ! number of bins
  integer :: n_spec                     ! number of species

  integer, parameter :: state_unit = 33 ! unit to use for state files
  integer, parameter :: scal = 3        ! scale factor for bins
  real*8, parameter :: v_min = 1d-24    ! minimum volume for making grid (m^3)
  integer, parameter :: spec_1 = 1      ! first solute species
  integer, parameter :: spec_2 = 2      ! second solute species
  real*8, parameter :: cutoff_frac = 0.01d0 ! fraction to count as mixed
  integer, parameter :: n_comp = 20     ! number of composition bins

  integer, allocatable :: MH(:)         ! number of particles per bin
  type(bin_p), allocatable :: VH(:)     ! particle volumes (m^3)
  real*8, allocatable :: bin_v(:)       ! volume of particles in bins
  type(bin_grid_t) :: bin_grid          ! bin grid

  type(environ) :: env                  ! environment state
  type(aero_data_t) :: aero_data                 ! aerosol data
  real*8 :: time                        ! current time (s)
  real*8 :: dlnr                        ! bin scale factor
  real*8, allocatable :: bin_g(:)       ! volume in bins
  real*8, allocatable :: bin_gs(:,:)    ! species volume in bins
  integer, allocatable :: bin_n(:)      ! number in bins
  integer, allocatable :: bin_n_2d(:,:) ! 2D species number distribution
  real*8, allocatable :: bin_g_2d(:,:)  ! 2D species volume distribution
  integer, allocatable :: bin_n_mixed(:,:) ! species number by composition
  real*8, allocatable :: bin_g_mixed(:,:) ! species volume by composition
  character :: filename*100             ! input filename
  character :: basename*100             ! basename of the input filename
  integer :: comp_n(n_comp)             ! number in composition bins

  call get_filename(filename, basename)

  call read_state_header(state_unit, filename, n_bin, n_spec)

  allocate(MH(n_bin))
  allocate(VH(n_bin))
  allocate(bin_v(n_bin))
  allocate(bin_g(n_bin))
  allocate(bin_gs(n_bin,n_spec))
  allocate(bin_n(n_bin))
  allocate(bin_g_2d(n_bin,n_bin))
  allocate(bin_n_2d(n_bin,n_bin))
  allocate(bin_n_mixed(n_bin,3))
  allocate(bin_g_mixed(n_bin,3))

  call init_array(n_spec, MH, VH)
  call read_state_bins(state_unit, filename, n_bin, bin_v, dlnr)
  call read_state(state_unit, filename, n_bin, n_spec, MH, VH, env, time)

  call make_bin_grid(n_bin, scal, v_min, bin_grid)
  ! FIXME: eventually delete following two lines
  bin_v = bin_grid%v
  dlnr = bin_grid%dlnr

  call allocate_aero_data(aero_data, n_spec)
  
  call moments(n_bin, n_spec, MH, VH, bin_v, &
       bin_g, bin_gs, bin_n, dlnr)
  call write_moments(basename, n_bin, n_spec, dlnr, env, bin_v, &
       bin_g, bin_gs, bin_n)

  if (n_spec > 1) then
       call moments_2d(n_bin, n_spec, MH, VH, bin_v, aero_data, &
          spec_1, spec_2, bin_n_2d, bin_g_2d)
       call write_moments_2d(basename, n_bin, dlnr, env, bin_v, &
          bin_n_2d, bin_g_2d)

       call moments_composition_2d(n_bin, n_spec, MH, VH, &
          spec_1, spec_2, n_comp, comp_n)
       call write_composition_2d(basename, n_comp, dlnr, env, comp_n)
  
       call moments_mixed_2d(n_bin, n_spec, MH, VH, bin_v, aero_data, &
           spec_1, spec_2, cutoff_frac, bin_n_mixed, bin_g_mixed)
       call write_moments_mixed_2d(basename, n_bin, dlnr, env, bin_v, &
          bin_n_mixed, bin_g_mixed)
       write(6,*) 'volume density in pure species 1: ', sum(bin_g_mixed(:,1))/env%V_comp
       write(6,*) 'volume density in pure species 2: ', sum(bin_g_mixed(:,2))/env%V_comp
       write(6,*) 'volume density in mixed: ', sum(bin_g_mixed(:,3))/env%V_comp
       write(6,*) 'volume density total: ', sum(bin_g)/env%V_comp
       write(6,*) 'volume density in species: ', sum(bin_gs,1)/env%V_comp
   else
       write(6,*) 'volume density total: ', sum(bin_g)/env%V_comp
   endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_filename(filename, basename)
    
    character, intent(out) :: filename*100 ! input filename
    character, intent(out) :: basename*100 ! basename of the input filename

    integer i
    
    ! check there is exactly one commandline argument
    if (iargc() .ne. 1) then
       write(0,*) 'Usage: process_state <filename.d>'
       call exit(1)
    end if
    
    ! get and check first commandline argument (must be "filename.d")
    call getarg(1, filename)
    i = len_trim(filename)
    if (i .gt. 40) then
       write(0,*) 'ERROR: filename too long'
       call exit(1)
    end if
    if ((filename(i:i) .ne. 'd') .or. &
         (filename((i-1):(i-1)) .ne. '.')) then
       write(0,*) 'ERROR: Filename must end in .d'
       call exit(1)
    end if
    
    ! chop .d off the end of the filename to get the basename
    basename = filename
    basename((i-1):i) = '  '
    
  end subroutine get_filename
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine moments_2d(n_bin, n_spec, MH, VH, bin_v, aero_data, &
       spec_1, spec_2, bin_n_2d, bin_g_2d)
    
    use mod_aero_data
    use mod_bin
    
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    integer, intent(in) :: MH(n_bin)    ! number of particles per bin
    type(bin_p), intent(in) :: VH(n_bin) ! particle volumes (m^3)
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    integer, intent(in) :: spec_1       ! first species
    integer, intent(in) :: spec_2       ! second species
    integer, intent(out) :: bin_n_2d(n_bin,n_bin) ! 2D species number dist
    real*8, intent(out) :: bin_g_2d(n_bin,n_bin) ! 2D species volume dist

    integer :: i, j, b1, b2
    
    bin_n_2d = 0
    bin_g_2d = 0d0
    
    do i = 1,n_bin
       do j = 1,MH(i)
          call particle_in_bin(VH(i)%p(j,spec_1), n_bin, bin_v, b1)
          call particle_in_bin(VH(i)%p(j,spec_2), n_bin, bin_v, b2)
          bin_n_2d(b1,b2) = bin_n_2d(b1,b2) + 1
          bin_g_2d(b1,b2) = bin_g_2d(b1,b2) + particle_volume(VH(i)%p(j,:))
       end do
    end do
    
  end subroutine moments_2d
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine moments_composition_2d(n_bin, n_spec, MH, VH, &
       spec_1, spec_2, n_comp, comp_n)
    
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    integer, intent(in) :: MH(n_bin)    ! number of particles per bin
    type(bin_p), intent(in) :: VH(n_bin) ! particle volumes (m^3)
    integer, intent(in) :: spec_1       ! first species
    integer, intent(in) :: spec_2       ! second species
    integer, intent(in) :: n_comp       ! number of composition bins
    integer, intent(out) :: comp_n(n_comp) ! number in composition bins
    
    integer :: i, j, i_comp
    real*8 :: comp
    
    comp_n = 0
    
    do i = 1,n_bin
       do j = 1,MH(i)
          comp = VH(i)%p(j,spec_1) / (VH(i)%p(j,spec_1) + VH(i)%p(j,spec_2))
          i_comp = floor(comp * dble(n_comp)) + 1
          if (i_comp .gt. n_comp) i_comp = n_comp
          comp_n(i_comp) = comp_n(i_comp) + 1
       end do
    end do
    
  end subroutine moments_composition_2d
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine moments_mixed_2d(n_bin, n_spec, MH, VH, bin_v, aero_data, &
       spec_1, spec_2, cutoff_frac, bin_n_mixed, bin_g_mixed)
    
    use mod_aero_data
    use mod_bin
    
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    integer, intent(in) :: MH(n_bin)    ! number of particles per bin
    type(bin_p), intent(in) :: VH(n_bin) ! particle volumes (m^3)
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    integer, intent(in) :: spec_1       ! first species
    integer, intent(in) :: spec_2       ! second species
    real*8, intent(in) :: cutoff_frac   ! fraction to count as mixed
    integer, intent(out) :: bin_n_mixed(n_bin,3) ! species number by composition
    real*8, intent(out) :: bin_g_mixed(n_bin,3) ! species volume by composition
    
    integer :: i, j, b, k
    real*8 :: comp, pv
    
    bin_n_mixed = 0
    bin_g_mixed = 0d0
    
    do i = 1,n_bin
       do j = 1,MH(i)
          comp = VH(i)%p(j,spec_1) / (VH(i)%p(j,spec_1) + VH(i)%p(j,spec_2))
          if (comp .lt. cutoff_frac) then
             k = 2
          else if (comp .gt. (1d0 - cutoff_frac)) then
             k = 1
          else
             k = 3
          end if
          pv = particle_volume(VH(i)%p(j,:))
          call particle_in_bin(pv, n_bin, bin_v, b)
          bin_n_mixed(b,k) = bin_n_mixed(b,k) + 1
          bin_g_mixed(b,k) = bin_g_mixed(b,k) + pv
       end do
    end do
    
  end subroutine moments_mixed_2d
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine write_moments(basename, n_bin, n_spec, dlnr, env, &
             bin_v, bin_g, bin_gs, bin_n)
    
    use mod_util
    
    character, intent(in) :: basename*100 ! basename of the input filename
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    real*8, intent(in) :: dlnr          ! bin scale factor
    type(environ) :: env                ! environment state
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins
    real*8, intent(in) :: bin_g(n_bin)  ! volume in bins
    real*8, intent(in) :: bin_gs(n_bin,n_spec) ! species volume in bins
    integer, intent(in) :: bin_n(n_bin) ! number in bins
    
    integer, parameter :: f_out = 20
    
    character :: filename*100
    integer :: i
    
    filename = basename
    i = len_trim(filename)
    filename((i+1):) = '_moments.d'
    open(f_out, file=filename)
    
    do i = 1,n_bin
       write(f_out,'(i20,e20.10,e20.10,e20.10)') i, vol2rad(bin_v(i)), &
            dble(bin_n(i)) / dlnr / env%V_comp, &
            bin_g(i) / dlnr / env%V_comp
    end do
    
    close(unit=f_out)
    
  end subroutine write_moments
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine write_moments_2d(basename, n_bin, dlnr, env, bin_v, &
             bin_n_2d, bin_g_2d)
    
    use mod_util
    
    character, intent(in) :: basename*100 ! basename of the input filename
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: dlnr          ! bin scale factor
    type(environ) :: env                ! environment state
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins
    integer, intent(in) :: bin_n_2d(n_bin,n_bin) ! 2D species number dist
    real*8, intent(in) :: bin_g_2d(n_bin,n_bin) ! 2D species volume dist
    
    integer, parameter :: f_out = 20
    
    character :: filename*100
    integer :: i, j
    
    filename = basename
    i = len_trim(filename)
    filename((i+1):) = '_moments_2d.d'
    open(f_out, file=filename)
    
    do i = 1,n_bin
       do j = 1,n_bin
          write(f_out,'(i20,i20,e20.10,e20.10,e20.10,e20.10)') &
               i, j, vol2rad(bin_v(i)), vol2rad(bin_v(j)), &
               dble(bin_n_2d(i,j)) / dlnr / env%V_comp, &
               bin_g_2d(i,j) / dlnr / env%V_comp
       end do
    end do
    
    close(unit=f_out)
    
  end subroutine write_moments_2d
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine write_composition_2d(basename, n_comp, dlnr, env, comp_n)
    
    character, intent(in) :: basename*100 ! basename of the input filename
    integer, intent(in) :: n_comp       ! number of composition bins
    real*8, intent(in) :: dlnr          ! bin scale factor
    type(environ) :: env                ! environment state
    integer, intent(in) :: comp_n(n_comp) ! number in composition bins

    integer, parameter :: f_out = 20
    
    character :: filename*100
    integer :: i
    
    filename = basename
    i = len_trim(filename)
    filename((i+1):) = '_composition.d'
    open(f_out, file=filename)
    
    do i = 1,n_comp
       write(f_out,'(i20,e20.10)') i, dble(comp_n(i)) / dlnr / env%V_comp
    end do
    
    close(unit=f_out)
    
  end subroutine write_composition_2d
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine write_moments_mixed_2d(basename, n_bin, dlnr, env, bin_v, &
       bin_n_mixed, bin_g_mixed)
    
    use mod_util
    
    character, intent(in) :: basename*100 ! basename of the input filename
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: dlnr          ! bin scale factor
    type(environ) :: env                ! environment state
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins
    integer, intent(in) :: bin_n_mixed(n_bin,3) ! species number by composition
    real*8, intent(in) :: bin_g_mixed(n_bin,3) ! species volume by composition
    
    integer, parameter :: f_out = 20
    
    character :: filename*100
    integer :: i
    
    filename = basename
    i = len_trim(filename)
    filename((i+1):) = '_moments_comp.d'
    open(f_out, file=filename)
    
    do i = 1,n_bin
       write(f_out,'(i20,e20.10,3e20.10,3e20.10)') i, vol2rad(bin_v(i)), &
            dble(bin_n_mixed(i,:)) / dlnr / env%V_comp, &
            bin_g_mixed(i,:) / dlnr / env%V_comp
    end do
    
    close(unit=f_out)
    
  end subroutine write_moments_mixed_2d
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end program process_state
