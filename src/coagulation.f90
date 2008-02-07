! Copyright (C) 2005-2008 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_coagulation module.

!> Aerosol particle coagulation.
module pmc_coagulation

  use pmc_bin_grid
  use pmc_aero_data
  use pmc_util
  use pmc_env_state
  use pmc_aero_state
  use pmc_aero_binned
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Choose a random pair for potential coagulation and test its
  !> probability of coagulation. If it happens, do the coagulation and
  !> update all structures.
  !!
  !! The probability of a coagulation will be taken as <tt>(kernel /
  !! k_max)</tt>.
  subroutine maybe_coag_pair(bin_grid, aero_binned, env_state, aero_data, &
       aero_state, b1, b2, del_t, k_max, kernel, did_coag)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Binned distributions.
    type(aero_binned_t), intent(out) :: aero_binned
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin of first particle.
    integer, intent(in) :: b1
    !> Bin of second particle.
    integer, intent(in) :: b2
    !> Timestep.
    real*8, intent(in) :: del_t
    !> K_max scale factor.
    real*8, intent(in) :: k_max
    !> Whether a coagulation occured.
    logical, intent(out) :: did_coag
    
    interface
       subroutine kernel(v1, v2, env_state, k)
         use pmc_env_state
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(env_state_t), intent(in) :: env_state   
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
    
    integer s1, s2
    real*8 p, k, pv1, pv2
    
    did_coag = .false.
    
    if ((aero_state%bins(b1)%n_part <= 0) &
         .or. (aero_state%bins(b2)%n_part <= 0)) then
       return
    end if
    
    call find_rand_pair(aero_state, b1, b2, s1, s2)
    pv1 = aero_particle_volume(aero_state%bins(b1)%particle(s1))
    pv2 = aero_particle_volume(aero_state%bins(b2)%particle(s2))
    call kernel(pv1, pv2, env_state, k)
    p = k / k_max
    
    if (pmc_rand() .lt. p) then
       call coagulate(bin_grid, aero_binned, aero_data, aero_state, &
            b1, s1, b2, s2)
       did_coag = .true.
    end if
    
  end subroutine maybe_coag_pair
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Given bins b1 and b2, find a random pair of particles (b1, s1)
  !> and (b2, s2) that are not the same particle particle as each
  !> other.
  subroutine find_rand_pair(aero_state, b1, b2, s1, s2)
    
    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Bin number of first particle.
    integer, intent(in) :: b1
    !> Bin number of second particle.
    integer, intent(in) :: b2
    !> First rand particle.
    integer, intent(out) :: s1
    !> Second rand particle.
    integer, intent(out) :: s2

    if ((aero_state%bins(b1)%n_part < 1) &
         .or. (aero_state%bins(b2)%n_part < 1) &
         .or. ((b1 == b2) .and. (aero_state%bins(b1)%n_part < 2))) then
       write(*,*) 'ERROR: find_rand_pair(): insufficient particles in bins', &
            b1, b2
       call exit(1)
    end if
    
    ! FIXME: rand() only returns a REAL*4, so we might not be able to
    ! generate all integers between 1 and M if M is too big.

100 s1 = int(pmc_rand() * dble(aero_state%bins(b1)%n_part)) + 1
    if ((s1 .lt. 1) .or. (s1 .gt. aero_state%bins(b1)%n_part)) goto 100
101 s2 = int(pmc_rand() * dble(aero_state%bins(b2)%n_part)) + 1
    if ((s2 .lt. 1) .or. (s2 .gt. aero_state%bins(b2)%n_part)) goto 101
    if ((b1 .eq. b2) .and. (s1 .eq. s2)) goto 101

! FIXME: enable this and delete the above junk
!    do
!       s1 = pmc_rand_int(aero_state%bins(b1)%n_part)
!       s2 = pmc_rand_int(aero_state%bins(b2)%n_part)
!       if (.not. ((b1 .eq. b2) .and. (s1 .eq. s2))) then
!          exit
!       end if
!    end do
    
  end subroutine find_rand_pair
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Join together particles (b1, s1) and (b2, s2), updating all
  !> particle and bin structures to reflect the change.
  subroutine coagulate(bin_grid, aero_binned, aero_data, aero_state, &
       b1, s1, b2, s2)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Binned distributions.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> First particle (bin number).
    integer, intent(in) :: b1
    !> First particle (number in bin).
    integer, intent(in) :: s1
    !> Second particle (bin number).
    integer, intent(in) :: b2
    !> Second particle (number in bin).
    integer, intent(in) :: s2
    
    type(aero_particle_t) :: new_particle
    integer :: bn

    call aero_particle_alloc(new_particle, aero_data%n_spec)

    ! coagulate particles
    call aero_particle_coagulate(aero_state%bins(b1)%particle(s1), &
         aero_state%bins(b2)%particle(s2), new_particle)
    bn = aero_particle_in_bin(new_particle, bin_grid)

    ! update binned data
    call aero_binned_remove_particle_in_bin(aero_binned, bin_grid, &
         b1, aero_state%comp_vol, aero_state%bins(b1)%particle(s1))
    call aero_binned_remove_particle_in_bin(aero_binned, bin_grid, &
         b2, aero_state%comp_vol, aero_state%bins(b2)%particle(s2))
    call aero_binned_add_particle_in_bin(aero_binned, bin_grid, &
         bn, aero_state%comp_vol, new_particle)

    ! remove old particles and add new one
    if ((b1 == b2) .and. (s2 > s1)) then
       ! handle a tricky corner case where we have to watch for s2 or
       ! s1 being the last entry in the array and being repacked when
       ! the other one is removed
       call aero_state_remove_particle(aero_state, b2, s2)
       call aero_state_remove_particle(aero_state, b1, s1)
    else
       call aero_state_remove_particle(aero_state, b1, s1)
       call aero_state_remove_particle(aero_state, b2, s2)
    end if
    call aero_state_add_particle(aero_state, bn, new_particle)

    call aero_particle_free(new_particle)
    
  end subroutine coagulate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_coagulation
