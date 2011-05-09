! Copyright (C) 2005-2011 Nicole Riemer and Matthew West
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
  use pmc_aero_weight
  use pmc_mpi
  use pmc_coag_kernel
  use pmc_aero_sorted
#ifdef PMC_USE_MPI
  use mpi
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do coagulation for time del_t.
  subroutine mc_coag(coag_kernel_type, env_state, aero_data, aero_weight, &
       aero_state, del_t, tot_n_samp, tot_n_coag)

    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Timestep for coagulation.
    real(kind=dp) :: del_t
    !> Total number of samples tested.
    integer, intent(out) :: tot_n_samp
    !> Number of coagulation events.
    integer, intent(out) :: tot_n_coag

    logical :: did_coag
    integer :: i, j, n_samp, i_samp
    real(kind=dp) :: accept_factor

    call aero_state_sort(aero_state)
    if (.not. aero_state%aero_sorted%coag_kernel_bounds_valid) then
       call est_k_minmax_binned(aero_state%aero_sorted%bin_grid, &
            coag_kernel_type, aero_data, aero_weight, env_state, &
            aero_state%aero_sorted%coag_kernel_min, &
            aero_state%aero_sorted%coag_kernel_max)
       aero_state%aero_sorted%coag_kernel_bounds_valid = .true.
    end if
       
    tot_n_samp = 0
    tot_n_coag = 0
    do i = 1,aero_state%aero_sorted%bin_grid%n_bin
       do j = 1,aero_state%aero_sorted%bin_grid%n_bin
          call compute_n_samp(aero_state%aero_sorted%bin(i)%n_entry, &
               aero_state%aero_sorted%bin(j)%n_entry, i == j, &
               aero_state%aero_sorted%coag_kernel_max(i,j), &
               aero_state%comp_vol, del_t, n_samp, accept_factor)
          tot_n_samp = tot_n_samp + n_samp
          do i_samp = 1,n_samp
             ! check we still have enough particles to coagulate
             if ((aero_state%aero_sorted%bin(i)%n_entry < 1) &
                  .or. (aero_state%aero_sorted%bin(j)%n_entry < 1) &
                  .or. ((i == j) &
                  .and. (aero_state%aero_sorted%bin(i)%n_entry < 2))) then
                exit
             end if
             call maybe_coag_pair(env_state, aero_data, aero_weight, &
                  aero_state, i, j, coag_kernel_type, accept_factor, did_coag)
             if (did_coag) tot_n_coag = tot_n_coag + 1
          end do
       end do
    end do

  end subroutine mc_coag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the number of samples required for the pair of bins.
  subroutine compute_n_samp(ni, nj, same_bin, k_max, comp_vol, &
       del_t, n_samp, accept_factor)

    !> Number particles in first bin.
    integer, intent(in) :: ni
    !> Number particles in second bin.
    integer, intent(in) :: nj
    !> Whether first bin is second bin.
    logical, intent(in) :: same_bin
    !> Maximum kernel value (s^{-1}).
    real(kind=dp), intent(in) :: k_max
    !> Computational volume (m^3).
    real(kind=dp), intent(in) :: comp_vol
    !> Timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Number of samples per timestep.
    integer, intent(out) :: n_samp
    !> Scale factor for accept probability (1).
    real(kind=dp), intent(out) :: accept_factor
    
    real(kind=dp) :: r_samp, n_samp_mean
    real(kind=dp) :: n_possible ! use real(kind=dp) to avoid integer overflow
    ! could use integer*8 or integer(kind = 8)
    ! or di = selected_int_kind(18), integer(kind=di)
    ! to represent 10^{-18} to 10^{18}
    
    if (same_bin) then
       ! don't change this to ni * (ni - 1) as the ni/nj distinction
       ! is important for coagulation_dist, which also calls this
       n_possible = real(ni, kind=dp) * (real(nj, kind=dp) - 1d0) / 2d0
    else
       n_possible = real(ni, kind=dp) * real(nj, kind=dp) / 2d0
    end if
    
    r_samp = k_max / comp_vol * del_t
    n_samp_mean = r_samp * n_possible
    n_samp = rand_poisson(n_samp_mean)
    accept_factor = 1d0 / k_max

    ! possible variants:
    ! A: accept_factor = 1d0 / k_max
    ! B: accept_factor = del_t * n_possible &
    !                    / (real(n_samp, kind=dp) * comp_vol)
    ! timings of test suite as of 2010-12-22T17:12:14-0600:
    !   A with n_samp = prob_round(n_samp_mean):
    !       159.82 162.18 156.28
    !   A with n_samp = rand_poisson(n_samp_mean):
    !       151.93 158.38 174.74 157.65
    !   B with n_samp = ceiling(n_samp_mean):
    !       196.06 200.23 192.41
    !   B with n_samp = ceiling(n_samp_mean + 1 * sqrt(n_samp_mean)):
    !       189.78 211.12 195.19
    !   B with n_samp = ceiling(n_samp_mean + 3 * sqrt(n_samp_mean)):
    !       214.60 201.25 203.55
    
  end subroutine compute_n_samp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Choose a random pair for potential coagulation and test its
  !> probability of coagulation. If it happens, do the coagulation and
  !> update all structures.
  !!
  !! The probability of a coagulation will be taken as <tt>(kernel /
  !! k_max)</tt>.
  subroutine maybe_coag_pair(env_state, aero_data, aero_weight, aero_state, &
       b1, b2, coag_kernel_type, accept_factor, did_coag)

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
    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Scale factor for accept probability (1).
    real(kind=dp), intent(in) :: accept_factor
    !> Whether a coagulation occured.
    logical, intent(out) :: did_coag
    
    integer :: s1, s2
    integer :: p1, p2
    real(kind=dp) :: p, k
    
    did_coag = .false.
    
    call assert(210827476, aero_state%aero_sorted%bin(b1)%n_entry >= 1)
    call assert(368973460, aero_state%aero_sorted%bin(b2)%n_entry >= 1)
    if (b1 == b2) then
       call assert(528541565, aero_state%aero_sorted%bin(b1)%n_entry >= 2)
    end if
    
    call find_rand_pair(aero_state%aero_sorted, b1, b2, s1, s2)
    p1 = aero_state%aero_sorted%bin(b1)%entry(s1)
    p2 = aero_state%aero_sorted%bin(b2)%entry(s2)
    call weighted_kernel(coag_kernel_type, aero_state%p%particle(p1), &
         aero_state%p%particle(p2), aero_data, aero_weight, &
         env_state, k)
    p = k * accept_factor

    if (pmc_random() .lt. p) then
       call coagulate(aero_data, aero_weight, aero_state, p1, p2)
       did_coag = .true.
    end if
    
  end subroutine maybe_coag_pair
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Given bins b1 and b2, find a random pair of particles (b1, s1)
  !> and (b2, s2) that are not the same particle particle as each
  !> other.
  subroutine find_rand_pair(aero_sorted, b1, b2, s1, s2)
    
    !> Aerosol sorted data.
    type(aero_sorted_t), intent(in) :: aero_sorted
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
    call assert(619608562, aero_sorted%bin(b1)%n_entry >= 1)
    call assert(271635751, aero_sorted%bin(b2)%n_entry >= 1)
    if (b1 == b2) then
       call assert(956184336, aero_sorted%bin(b1)%n_entry >= 2)
    end if
    
    ! FIXME: don't loop, just do:
    ! if (b1 == b2) then
    !    s2 = pmc_rand_int(aero_sorted%bin(b2)%n_entry - 1)
    !    if (s2 == s1) then
    !       s2 = aero_sorted%bin(b2)%n_entry
    !    end if
    ! else
    !    s2 = pmc_rand_int(aero_sorted%bin(b2)%n_entry)
    ! end if
    do
       s1 = pmc_rand_int(aero_sorted%bin(b1)%n_entry)
       s2 = pmc_rand_int(aero_sorted%bin(b2)%n_entry)
       if ((b1 /= b2) .or. (s1 /= s2)) then
          ! stop generating if we have two distinct particles
          exit
       end if
    end do
    
  end subroutine find_rand_pair
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Actually coagulate particle_1 and particle_2 to form particle_new
  !> and compute weighting effects, including which particles should
  !> be lost and which gained.
  subroutine coagulate_weighting(particle_1, particle_2, particle_new, &
       aero_data, aero_weight, remove_1, remove_2, create_new, &
       id_1_lost, id_2_lost, aero_info_1, aero_info_2)

    !> First coagulating aerosol particle.
    type(aero_particle_t), intent(in) :: particle_1
    !> Second coagulating aerosol particle.
    type(aero_particle_t), intent(in) :: particle_2
    !> Combined aerosol particle resulting from coagulation of particle_1
    !> and particle_2.
    type(aero_particle_t), intent(inout) :: particle_new
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Whether to remove particle_1.
    logical, intent(out) :: remove_1
    !> Whether to remove particle_2.
    logical, intent(out) :: remove_2
    !> Whether to create particle_new.
    logical, intent(out) :: create_new
    !> Whether the ID of particle_1 will be lost due to coagulation.
    logical, intent(out) :: id_1_lost
    !> Whether the ID of particle_2 will be lost due to coagulation.
    logical, intent(out) :: id_2_lost
    !> The removal information associated with particle_1.
    type(aero_info_t), intent(inout) :: aero_info_1
    !> The removal information associated with particle_2.
    type(aero_info_t), intent(inout) :: aero_info_2

    real(kind=dp) :: radius_1, radius_2, radius_new
    real(kind=dp) :: weight_1, weight_2, weight_new, weight_min
    real(kind=dp) :: prob_remove_1, prob_remove_2, prob_create_new
    integer :: info_other_id

    call assert(371947172, particle_1%id /= particle_2%id)

    ! decide which old particles are to be removed and whether to
    ! create the resulting coagulated particle
    if (aero_weight%type == AERO_WEIGHT_TYPE_NONE) then
       remove_1 = .true.
       remove_2 = .true.
       create_new = .true.
    elseif ((aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
         .or. (aero_weight%type == AERO_WEIGHT_TYPE_MFA)) then
       radius_1 = aero_particle_radius(particle_1)
       radius_2 = aero_particle_radius(particle_2)
       radius_new = vol2rad(rad2vol(radius_1) + rad2vol(radius_2))
       weight_1 = aero_weight_value(aero_weight, radius_1)
       weight_2 = aero_weight_value(aero_weight, radius_2)
       weight_new = aero_weight_value(aero_weight, radius_new)
       weight_min = min(weight_1, weight_2, weight_new)
       prob_remove_1 = weight_min / weight_1
       prob_remove_2 = weight_min / weight_2
       prob_create_new = weight_min / weight_new
       remove_1 = (pmc_random() < prob_remove_1)
       if (aero_weight%type == AERO_WEIGHT_TYPE_MFA) then
          remove_2 = .not. remove_1
       else
          remove_2 = (pmc_random() < prob_remove_2)
       end if
       create_new = (pmc_random() < prob_create_new)
    else
       call die_msg(886524113, "unknown aero_weight type: " &
            // trim(integer_to_string(aero_weight%type)))
    end if

    ! figure out what to do about the ID numbers of the various
    ! particles --- we try to preserve particle IDs as much as
    ! possible
    if (create_new) then
       id_1_lost = .false.
       id_2_lost = .false.
       if (remove_1 .and. remove_2) then
          if (aero_particle_volume(particle_1) &
               > aero_particle_volume(particle_2)) then
             id_2_lost = .true.
          else
             id_1_lost = .true.
          end if
       end if
    else
       id_1_lost = remove_1
       id_2_lost = remove_2
    end if

    ! create a new particle and set its ID
    if (create_new) then
       call aero_particle_deallocate(particle_new)
       call aero_particle_allocate_size(particle_new, aero_data%n_spec, &
            aero_data%n_source)
       call aero_particle_coagulate(particle_1, particle_2, particle_new)
       if (remove_1 .and. (.not. id_1_lost)) then
          particle_new%id = particle_1%id
          call assert(975059559, id_2_lost .eqv. remove_2)
       elseif (remove_2 .and. (.not. id_2_lost)) then
          particle_new%id = particle_2%id
          call assert(246529753, id_1_lost .eqv. remove_1)
       else
          call aero_particle_new_id(particle_new)
          call assert(852038606, id_1_lost .eqv. remove_1)
          call assert(254018921, id_2_lost .eqv. remove_2)
       end if
       info_other_id = particle_new%id
    else
       info_other_id = 0
    end if

    aero_info_1%id = particle_1%id
    aero_info_1%action = AERO_INFO_COAG
    aero_info_1%other_id = info_other_id

    aero_info_2%id = particle_2%id
    aero_info_2%action = AERO_INFO_COAG
    aero_info_2%other_id = info_other_id

  end subroutine coagulate_weighting

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Join together particles (b1, s1) and (b2, s2), updating all
  !> particle and bin structures to reflect the change.
  subroutine coagulate(aero_data, aero_weight, aero_state, p1, p2)
 
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> First particle index.
    integer, intent(in) :: p1
    !> Second particle index.
    integer, intent(in) :: p2

    type(aero_particle_t), pointer :: particle_1, particle_2
    type(aero_particle_t) :: particle_new
    integer :: bn
    type(aero_info_t) :: aero_info_1, aero_info_2
    logical :: remove_1, remove_2, create_new, id_1_lost, id_2_lost

    call aero_particle_allocate(particle_new)
    call aero_info_allocate(aero_info_1)
    call aero_info_allocate(aero_info_2)

    particle_1 => aero_state%p%particle(p1)
    particle_2 => aero_state%p%particle(p2)

    call coagulate_weighting(particle_1, particle_2, particle_new, &
         aero_data, aero_weight, remove_1, remove_2, create_new, &
         id_1_lost, id_2_lost, aero_info_1, aero_info_2)

    ! remove old particles
    if (p2 > p1) then
       ! handle a tricky corner case where we have to watch for p2 or
       ! p1 being the last entry in the array and being repacked when
       ! the other one is removed
       if (remove_2) then
          call aero_state_remove_particle(aero_state, p2, &
               id_2_lost, aero_info_2)
       end if
       if (remove_1) then
          call aero_state_remove_particle(aero_state, p1, &
               id_1_lost, aero_info_1)
       end if
    else
       if (remove_1) then
          call aero_state_remove_particle(aero_state, p1, &
               id_1_lost, aero_info_1)
       end if
       if (remove_2) then
          call aero_state_remove_particle(aero_state, p2, &
               id_2_lost, aero_info_2)
       end if
    end if

    ! add new particle
    if (create_new) then
       call aero_state_add_particle(aero_state, particle_new, &
            allow_resort=.false.)
    end if

    call aero_info_deallocate(aero_info_1)
    call aero_info_deallocate(aero_info_2)
    call aero_particle_deallocate(particle_new)
    
  end subroutine coagulate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_coagulation
