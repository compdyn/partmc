! Copyright (C) 2005-2010 Nicole Riemer and Matthew West
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
  use pmc_mpi
  use pmc_kernel
#ifdef PMC_USE_MPI
  use mpi
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the number of samples required for the pair of bins.
  subroutine compute_n_samp(ni, nj, same_bin, k_max, comp_vol, &
       del_t, n_samp_real)

    !> Number particles in first bin .
    integer, intent(in) :: ni
    !> Number particles in second bin.
    integer, intent(in) :: nj
    !> Whether first bin is second bin.
    logical, intent(in) :: same_bin
    !> Maximum kernel value.
    real(kind=dp), intent(in) :: k_max
    !> Computational volume (m^3).
    real(kind=dp), intent(in) :: comp_vol
    !> Timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Number of samples per timestep.
    real(kind=dp), intent(out) :: n_samp_real
    
    real(kind=dp) :: r_samp
    real(kind=dp) :: n_possible ! use real(kind=dp) to avoid integer overflow
    ! FIXME: should use integer*8 or integer(kind = 8)
    ! or even better, di = selected_int_kind(18), integer(kind=di)
    ! to represent 10^{-18} to 10^{18}
    
    if (same_bin) then
       n_possible = real(ni, kind=dp) * (real(nj, kind=dp) - 1d0) / 2d0
    else
       n_possible = real(ni, kind=dp) * real(nj, kind=dp) / 2d0
    endif
    
    r_samp = k_max / comp_vol * del_t
    n_samp_real = r_samp * n_possible
    
  end subroutine compute_n_samp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Choose a random pair for potential coagulation and test its
  !> probability of coagulation. If it happens, do the coagulation and
  !> update all structures.
  !!
  !! The probability of a coagulation will be taken as <tt>(kernel /
  !! k_max)</tt>.
  subroutine maybe_coag_pair(bin_grid, env_state, aero_data, &
       aero_weight, aero_state, b1, b2, del_t, k_max, kernel, did_coag)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin of first particle.
    integer, intent(in) :: b1
    !> Bin of second particle.
    integer, intent(in) :: b2
    !> Timestep.
    real(kind=dp), intent(in) :: del_t
    !> K_max scale factor.
    real(kind=dp), intent(in) :: k_max
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
         real(kind=dp), intent(out) :: k
       end subroutine kernel
    end interface
#endif
    
    integer :: s1, s2
    real(kind=dp) :: p, k
    
    did_coag = .false.
    
    call assert(210827476, aero_state%bin(b1)%n_part >= 1)
    call assert(368973460, aero_state%bin(b2)%n_part >= 1)
    if (b1 == b2) then
       call assert(528541565, aero_state%bin(b1)%n_part >= 2)
    end if
    
    call find_rand_pair(aero_state, b1, b2, s1, s2)
    call weighted_kernel(kernel, aero_state%bin(b1)%particle(s1), &
         aero_state%bin(b2)%particle(s2), aero_data, aero_weight, &
         env_state, k)
    p = k / k_max
    
    if (pmc_random() .lt. p) then
       call coagulate(bin_grid, aero_data, aero_weight, aero_state, &
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

100 s1 = int(pmc_random() * real(aero_state%bin(b1)%n_part, kind=dp)) + 1
    if ((s1 .lt. 1) .or. (s1 .gt. aero_state%bin(b1)%n_part)) goto 100
101 s2 = int(pmc_random() * real(aero_state%bin(b2)%n_part, kind=dp)) + 1
    if ((s2 .lt. 1) .or. (s2 .gt. aero_state%bin(b2)%n_part)) goto 101
    if ((b1 .eq. b2) .and. (s1 .eq. s2)) goto 101

! FIXME: enable this and delete the above junk
!    do
!       s1 = pmc_rand_int(aero_state%bin(b1)%n_part)
!       s2 = pmc_rand_int(aero_state%bin(b2)%n_part)
!       if (.not. ((b1 .eq. b2) .and. (s1 .eq. s2))) then
!          ! stop generating if we have two distinct particles
!          exit
!       end if
!    end do
    
  end subroutine find_rand_pair
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Join together particles (b1, s1) and (b2, s2), updating all
  !> particle and bin structures to reflect the change.
  subroutine coagulate(bin_grid, aero_data, aero_weight, aero_state, &
       b1, s1, b2, s2)
 
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
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
    type(aero_particle_t) :: particle_new
    integer :: bn
    type(aero_info_t) :: aero_info_1, aero_info_2
    real(kind=dp) :: weight_1, weight_2, weight_new
    real(kind=dp) :: prob_remove_1, prob_remove_2
    logical :: remove_1, remove_2, id_1_lost, id_2_lost

    call aero_particle_allocate_size(particle_new, aero_data%n_spec)
    particle_1 => aero_state%bin(b1)%particle(s1)
    particle_2 => aero_state%bin(b2)%particle(s2)
    call assert(371947172, particle_1%id /= particle_2%id)

    ! coagulate particles
    call aero_particle_coagulate(particle_1, particle_2, particle_new)
    bn = aero_particle_in_bin(particle_new, bin_grid)

    ! decide which old particles are to be removed
    if (aero_weight%type == AERO_WEIGHT_TYPE_NONE) then
       remove_1 = .true.
       remove_2 = .true.
    else
       weight_1 = aero_weight_value(aero_weight, &
            aero_particle_radius(particle_1))
       weight_2 = aero_weight_value(aero_weight, &
            aero_particle_radius(particle_2))
       weight_new = aero_weight_value(aero_weight, &
            aero_particle_radius(particle_new))
       prob_remove_1 = weight_new / weight_1
       prob_remove_2 = weight_new / weight_2
       remove_1 = (pmc_random() < prob_remove_1)
       remove_2 = (pmc_random() < prob_remove_2)
    end if

    if (remove_1 .and. remove_2) then
       if (aero_particle_volume(particle_1) &
            > aero_particle_volume(particle_2)) then
          particle_new%id = particle_1%id
          id_1_lost = .false.
          id_2_lost = .true.
       else
          particle_new%id = particle_2%id
          id_1_lost = .true.
          id_2_lost = .false.
       end if
    elseif (remove_1 .and. (.not. remove_2)) then
       particle_new%id = particle_1%id
       id_1_lost = .false.
       id_2_lost = .false.
    elseif ((.not. remove_1) .and. remove_2) then
       particle_new%id = particle_2%id
       id_1_lost = .false.
       id_2_lost = .false.
    else ! ((.not. remove_1) .and. (.not. remove_2))
       call aero_particle_new_id(particle_new)
       id_1_lost = .false.
       id_2_lost = .false.
    end if

    ! remove old particles
    call aero_info_allocate(aero_info_1)
    aero_info_1%id = particle_1%id
    aero_info_1%action = AERO_INFO_COAG
    aero_info_1%other_id = particle_new%id
    call aero_info_allocate(aero_info_2)
    aero_info_2%id = particle_2%id
    aero_info_2%action = AERO_INFO_COAG
    aero_info_2%other_id = particle_new%id
    if ((b1 == b2) .and. (s2 > s1)) then
       ! handle a tricky corner case where we have to watch for s2 or
       ! s1 being the last entry in the array and being repacked when
       ! the other one is removed
       if (remove_2) then
          call aero_state_remove_particle(aero_state, b2, s2, &
               id_2_lost, aero_info_2)
       end if
       if (remove_1) then
          call aero_state_remove_particle(aero_state, b1, s1, &
               id_1_lost, aero_info_1)
       end if
    else
       if (remove_1) then
          call aero_state_remove_particle(aero_state, b1, s1, &
               id_1_lost, aero_info_1)
       end if
       if (remove_2) then
          call aero_state_remove_particle(aero_state, b2, s2, &
               id_2_lost, aero_info_2)
       end if
    end if
    call aero_info_deallocate(aero_info_1)
    call aero_info_deallocate(aero_info_2)

    ! add new particle
    call aero_state_add_particle(aero_state, bn, particle_new)
    call aero_particle_deallocate(particle_new)
    
  end subroutine coagulate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_coagulation
