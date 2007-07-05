! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Coagulation subroutines.

module mod_coagulation

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine maybe_coag_pair(bin_grid, aero_binned, env, aero_data, &
       aero_state, b1, b2, del_t, k_max, kernel, did_coag, bin_change)
    
    ! Choose a random pair for potential coagulation and test its
    ! probability of coagulation. If it happens, do the coagulation and
    ! update all structures. The probability of a coagulation will be
    ! taken as (kernel / k_max).

    use mod_bin
    use mod_aero_data
    use mod_util
    use mod_environ
    use mod_aero_state
    use mod_aero_binned

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_binned_t), intent(out) :: aero_binned ! binned distributions
    type(environ), intent(inout) :: env ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    integer, intent(in) :: b1           ! bin of first particle
    integer, intent(in) :: b2           ! bin of second particle
    real*8, intent(in) :: del_t         ! timestep
    real*8, intent(in) :: k_max         ! k_max scale factor
    logical, intent(out) :: did_coag    ! whether a coagulation occured
    logical, intent(out) :: bin_change  ! whether bin structure changed
    
    interface
       subroutine kernel(v1, v2, env, k)
         use mod_environ
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(environ), intent(in) :: env 
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
    
    integer s1, s2
    real*8 p, k, pv1, pv2
    
    bin_change = .false.
    did_coag = .false.
    
    if ((aero_state%n(b1) .le. 0) .or. (aero_state%n(b2) .le. 0)) then
       return
    end if
    
    call find_rand_pair(aero_state, b1, b2, s1, s2)
    pv1 = particle_volume(aero_state%v(b1)%p(s1,:))
    pv2 = particle_volume(aero_state%v(b2)%p(s2,:))
    call kernel(pv1, pv2, env, k)
    p = k / k_max
    
    if (util_rand() .lt. p) then
       call coagulate(bin_grid, aero_binned, env, aero_data, aero_state, &
            b1, s1, b2, s2, bin_change)
       did_coag = .true.
    end if
    
  end subroutine maybe_coag_pair
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_rand_pair(aero_state, b1, b2, s1, s2)
    
    ! Given bins b1 and b2, find a random pair of particles (b1, s1)
    ! and (b2, s2).
    
    use mod_util
    use mod_aero_state

    type(aero_state_t), intent(in) :: aero_state ! aerosol state
    integer, intent(in) :: b1           ! bin number of first particle
    integer, intent(in) :: b2           ! bin number of second particle
    integer, intent(out) :: s1          ! first rand particle 1 <= s1 <= M(b1)
    integer, intent(out) :: s2          ! second rand particle 1 <= s2 <= M(b2)
                                        !         (b1,s1) != (b2,s2)

    integer :: n_bin
    logical :: found

    n_bin = size(aero_state%n)
    if ((aero_state%n(b1) < 1) .or. (aero_state%n(b2) < 1) .or. &
         ((b1 == b2) .and. (aero_state%n(b1) < 2))) then
       write(*,*) 'ERROR: find_rand_pair(): insufficient particles in bins', &
            b1, b2
       call exit(1)
    end if
    
    ! FIXME: rand() only returns a REAL*4, so we might not be able to
    ! generate all integers between 1 and M if M is too big.

100 s1 = int(util_rand() * dble(aero_state%n(b1))) + 1
    if ((s1 .lt. 1) .or. (s1 .gt. aero_state%n(b1))) goto 100
101 s2 = int(util_rand() * dble(aero_state%n(b2))) + 1
    if ((s2 .lt. 1) .or. (s2 .gt. aero_state%n(b2))) goto 101
    if ((b1 .eq. b2) .and. (s1 .eq. s2)) goto 101

!    found = .false.
!    do while (.not. found)
!       found = .true.
!       s1 = int(util_rand() * dble(aero_state%n(b1))) + 1
!       if ((s1 .lt. 1) .or. (s1 .gt. aero_state%n(b1))) found = .false.
!       s2 = int(util_rand() * dble(aero_state%n(b2))) + 1
!       if ((s2 .lt. 1) .or. (s2 .gt. aero_state%n(b2))) found = .false.
!       if ((b1 .eq. b2) .and. (s1 .eq. s2)) found = .false.
!    end do
    
  end subroutine find_rand_pair
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine coagulate(bin_grid, aero_binned, env, aero_data, aero_state, &
       b1, s1, b2, s2, bin_change)

    ! Join together particles (b1, s1) and (b2, s2), updating all
    ! particle and bin structures to reflect the change. bin_change is
    ! true if the used bin set changed due to the coagulation (i.e. an
    ! empty bin filled or a filled bin became empty).

    use mod_aero_data
    use mod_bin
    use mod_environ
    use mod_aero_state
    use mod_aero_binned
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_binned_t), intent(out) :: aero_binned ! binned distributions
    type(environ), intent(inout) :: env ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(inout) :: aero_state ! aerosol state
    integer, intent(in) :: b1           ! first particle (bin number)
    integer, intent(in) :: s1           ! first particle (number in bin)
    integer, intent(in) :: b2           ! second particle (bin number)
    integer, intent(in) :: s2           ! second particle (number in bin)
    logical, intent(out) :: bin_change  ! whether an empty bin filled,
                                        ! or a filled bin became empty
    
    integer bn, i, j
    real*8 new_v(aero_data%n_spec), new_v_tot

    bin_change = .false.
    
    ! remove s1 and s2 from bins
    aero_binned%num_den(b1) = aero_binned%num_den(b1) &
         - 1d0 / aero_state%comp_vol
    aero_binned%num_den(b2) = aero_binned%num_den(b2) &
         - 1d0 / aero_state%comp_vol
    aero_binned%vol_den(b1,:) = aero_binned%vol_den(b1,:) &
         - aero_state%v(b1)%p(s1,:) / aero_state%comp_vol
    aero_binned%vol_den(b2,:) = aero_binned%vol_den(b2,:) &
         - aero_state%v(b2)%p(s2,:) / aero_state%comp_vol

    ! do coagulation in aero_state%n, aero_state%v arrays
    ! add particle volumes
    new_v(:) = aero_state%v(b1)%p(s1,:) + aero_state%v(b2)%p(s2,:)
    new_v_tot = sum(new_v)

    call particle_in_bin(new_v_tot, bin_grid, bn)  ! find new bin

    ! handle a tricky corner case
    if ((b1 .eq. b2) .and. (s2 .eq. aero_state%n(b2))) then
       aero_state%v(b1)%p(s1,:) = aero_state%v(b1)%p(aero_state%n(b1)-1,:)
       aero_state%n(b1) = aero_state%n(b1) - 2
    else
       ! shift last particle into empty slot
       aero_state%v(b1)%p(s1,:) = aero_state%v(b1)%p(aero_state%n(b1),:)
       aero_state%n(b1) = aero_state%n(b1) - 1 ! decrease length of array
       aero_state%v(b2)%p(s2,:) = aero_state%v(b2)%p(aero_state%n(b2),:)
       aero_state%n(b2) = aero_state%n(b2) - 1 ! same for second particle
    end if
    if ((aero_state%n(b1) .lt. 0) .or. (aero_state%n(b2) .lt. 0)) then
       write(0,*)'ERROR: invalid aero_state%n'
       call exit(2)
    end if
    aero_state%n(bn) = aero_state%n(bn) + 1 ! increase the length of array
    call enlarge_bin_to(aero_state%v(bn), aero_state%n(bn))
    aero_state%v(bn)%p(aero_state%n(bn),:) = new_v ! add new particle at end
    
    ! add new particle to bins
    aero_binned%num_den(bn) = aero_binned%num_den(bn) &
         + 1d0 / aero_state%comp_vol
    aero_binned%vol_den(bn,:) = aero_binned%vol_den(bn,:) &
         + aero_state%v(bn)%p(aero_state%n(bn),:) / aero_state%comp_vol

    ! did we empty a bin?
    if ((aero_state%n(b1) .eq. 0) .or. (aero_state%n(b2) .eq. 0)) &
         bin_change = .true.
    ! did we start a new bin?
    if ((aero_state%n(bn) .eq. 1) .and. (bn .ne. b1) .and. (bn .ne. b2)) &
         bin_change = .true.

    ! possibly repack memory
    call shrink_bin(aero_state%n(b1), aero_state%v(b1))
    call shrink_bin(aero_state%n(b2), aero_state%v(b2))

  end subroutine coagulate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_coagulation
