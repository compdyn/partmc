! Copyright (C) 2005-2009 Nicole Riemer and Matthew West
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
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Choose a random pair for potential coagulation and test its
  !> probability of coagulation. If it happens, do the coagulation and
  !> update all structures.
  !!
  !! The probability of a coagulation will be taken as <tt>(kernel /
  !! k_max)</tt>.
  subroutine maybe_coag_pair(bin_grid, env_state, aero_data, &
       aero_state, b1, b2, del_t, k_max, kernel, did_coag)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
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
    
#ifndef DOXYGEN_SKIP_DOC
    interface
       subroutine kernel(aero_particle_1, aero_particle_2, aero_data, &
            env_state, k)
         use pmc_aero_particle
         use pmc_aero_data
         use pmc_env_state
         type(aero_particle_t), intent(in) :: aero_particle_1
         type(aero_particle_t), intent(in) :: aero_particle_2
         type(aero_data_t), intent(in) :: aero_data
         type(env_state_t), intent(in) :: env_state  
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
#endif
    
    integer :: s1, s2
    real*8 :: p, k
    
    did_coag = .false.
    
    if ((aero_state%bin(b1)%n_part <= 0) &
         .or. (aero_state%bin(b2)%n_part <= 0)) then
       return
    end if
    
    call find_rand_pair(aero_state, b1, b2, s1, s2)
    call kernel(aero_state%bin(b1)%particle(s1), &
         aero_state%bin(b2)%particle(s2), aero_data, env_state, k)
    p = k / k_max
    
    if (pmc_random() .lt. p) then
       call coagulate(bin_grid, aero_data, aero_state, &
            b1, s1, b2, s2)
       did_coag = .true.
    end if
    
  end subroutine maybe_coag_pair
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

    ! check we have enough particles to avoid being stuck in an
    ! infinite loop below
    call assert(362349482, aero_state%bin(b1)%n_part >= 1)
    call assert(479121681, aero_state%bin(b2)%n_part >= 1)
    if (b1 == b2) then
       call assert(161928491, aero_state%bin(b1)%n_part >= 2)
    end if
    
    ! FIXME: rand() only returns a REAL*4, so we might not be able to
    ! generate all integers between 1 and M if M is too big.

100 s1 = int(pmc_random() * dble(aero_state%bin(b1)%n_part)) + 1
    if ((s1 .lt. 1) .or. (s1 .gt. aero_state%bin(b1)%n_part)) goto 100
101 s2 = int(pmc_random() * dble(aero_state%bin(b2)%n_part)) + 1
    if ((s2 .lt. 1) .or. (s2 .gt. aero_state%bin(b2)%n_part)) goto 101
    if ((b1 .eq. b2) .and. (s1 .eq. s2)) goto 101

! FIXME: enable this and delete the above junk
!    do
!       s1 = pmc_rand_int(aero_state%bin(b1)%n_part)
!       s2 = pmc_rand_int(aero_state%bin(b2)%n_part)
!       if (.not. ((b1 .eq. b2) .and. (s1 .eq. s2))) then
!          exit
!       end if
!    end do
    
  end subroutine find_rand_pair
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Join together particles (b1, s1) and (b2, s2), updating all
  !> particle and bin structures to reflect the change.
  subroutine coagulate(bin_grid, aero_data, aero_state, &
       b1, s1, b2, s2)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
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
    
    type(aero_particle_t), pointer :: particle_1, particle_2
    type(aero_particle_t) :: new_particle
    integer :: bn
    type(aero_info_t) :: aero_info
    logical :: p1_removed, p2_removed

    call aero_particle_allocate_size(new_particle, aero_data%n_spec)
    particle_1 => aero_state%bin(b1)%particle(s1)
    particle_2 => aero_state%bin(b2)%particle(s2)
    call assert(371947172, particle_1%id /= particle_2%id)

    ! coagulate particles
    call aero_particle_coagulate(particle_1, particle_2, new_particle)
    bn = aero_particle_in_bin(new_particle, bin_grid)

    ! remove old particles
    call aero_info_allocate(aero_info)
    if (new_particle%id /= particle_1%id) then
       ! particle_1 is the removed particle
       call assert(361912382, new_particle%id == particle_2%id)
       aero_info%id = particle_1%id
       aero_info%action = AERO_INFO_COAG
       aero_info%other_id = particle_2%id
       p1_removed = .true.
       p2_removed = .false.
    else
       ! particle_2 is the removed particle
       call assert(742917292, new_particle%id /= particle_2%id)
       aero_info%id = particle_2%id
       aero_info%action = AERO_INFO_COAG
       aero_info%other_id = particle_1%id
       p1_removed = .false.
       p2_removed = .true.
    end if
    if ((b1 == b2) .and. (s2 > s1)) then
       ! handle a tricky corner case where we have to watch for s2 or
       ! s1 being the last entry in the array and being repacked when
       ! the other one is removed
       call aero_state_remove_particle(aero_state, b2, s2, &
            p2_removed, aero_info)
       call aero_state_remove_particle(aero_state, b1, s1, &
            p1_removed, aero_info)
    else
       call aero_state_remove_particle(aero_state, b1, s1, &
            p1_removed, aero_info)
       call aero_state_remove_particle(aero_state, b2, s2, &
            p2_removed, aero_info)
    end if
    call aero_info_deallocate(aero_info)

    ! add new particle
    call aero_state_add_particle(aero_state, bn, new_particle)
    call aero_particle_deallocate(new_particle)
    
  end subroutine coagulate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_coagulation
