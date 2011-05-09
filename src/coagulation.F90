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

    logical :: did_coag, do_accel_coag
    integer :: i_bin, j_bin, n_samp, n_coag, i_samp, j_entry, j_part, new_bin
    real(kind=dp) :: accept_factor, n_samp_mean, n_samp_per_j_part

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
    do i_bin = 1,aero_state%aero_sorted%bin_grid%n_bin
       do j_bin = i_bin,aero_state%aero_sorted%bin_grid%n_bin
          do_accel_coag = .false.
          if ((j_bin > i_bin) &
               .and. (aero_weight%type == AERO_WEIGHT_TYPE_POWER) &
               .and. (aero_weight%exponent < 0d0)) then
             n_samp_per_j_part &
                  = aero_state%aero_sorted%coag_kernel_max(i_bin, j_bin) &
                  / aero_state%comp_vol * del_t &
                  * real(aero_state%aero_sorted%bin(i_bin)%n_entry, kind=dp)
             if (n_samp_per_j_part > 3d0) then
                do_accel_coag = .true.
             end if
          end if
          if (do_accel_coag) then
             do j_entry = 1,aero_state%aero_sorted%bin(j_bin)%n_entry
                call accel_coag(env_state, aero_data, aero_weight, &
                     aero_state, i_bin, j_bin, j_entry, coag_kernel_type, &
                     del_t, n_samp, n_coag)
                tot_n_samp = tot_n_samp + n_samp
                tot_n_coag = tot_n_coag + n_coag
             end do
             ! particles might be in the wrong bins: work backwards to
             ! avoid particle movement issues
             do j_entry = aero_state%aero_sorted%bin(j_bin)%n_entry,1,-1
                j_part = aero_state%aero_sorted%bin(j_bin)%entry(j_entry)
                new_bin &
                     = aero_sorted_particle_in_bin(aero_state%aero_sorted, &
                     aero_state%p%particle(j_part))
                if (new_bin /= j_bin) then
                   call aero_sorted_move_particle(aero_state%aero_sorted, &
                        j_part, new_bin)
                end if
             end do
          else
             call compute_n_samp(aero_state%aero_sorted%bin(i_bin)%n_entry, &
                  aero_state%aero_sorted%bin(j_bin)%n_entry, i_bin == j_bin, &
                  aero_state%aero_sorted%coag_kernel_max(i_bin, j_bin), &
                  aero_state%comp_vol, del_t, n_samp_mean, n_samp, &
                  accept_factor)
             tot_n_samp = tot_n_samp + n_samp
             do i_samp = 1,n_samp
                ! check we still have enough particles to coagulate
                if ((aero_state%aero_sorted%bin(i_bin)%n_entry < 1) &
                     .or. (aero_state%aero_sorted%bin(j_bin)%n_entry < 1) &
                     .or. ((i_bin == j_bin) &
                     .and. (aero_state%aero_sorted%bin(i_bin)%n_entry < 2))) &
                     then
                   exit
                end if
                call maybe_coag_pair(env_state, aero_data, aero_weight, &
                     aero_state, i_bin, j_bin, coag_kernel_type, &
                     accept_factor, did_coag)
             end do
          end if
       end do
    end do

  end subroutine mc_coag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Accelerated coagulation for a single particle.
  subroutine accel_coag(env_state, aero_data, aero_weight, aero_state, &
       i_bin, j_bin, j_entry, coag_kernel_type, del_t, n_samp, n_coag)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Bin of first particle.
    integer, intent(in) :: i_bin
    !> Bin of second particle.
    integer, intent(in) :: j_bin
    !> Entry in bin of second particle.
    integer, intent(in) :: j_entry
    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Number of samples used.
    integer, intent(out) :: n_samp
    !> Number of coagulations performed.
    integer, intent(out) :: n_coag

    real(kind=dp) :: n_samp_mean, k_max, prob_remove_i, prob_remove_i_max
    real(kind=dp) :: prob_coag, prob_coag_tot, prob_coag_mean
    real(kind=dp) :: weight_i, weight_i_min, weight_j, k
    real(kind=dp) :: vol_sq(aero_data%n_spec), vol_mean(aero_data%n_spec)
    real(kind=dp) :: vol_cv(aero_data%n_spec), vol_cv_max, mean_95_conf_cv
    integer :: n_samp_remove, n_samp_extra, n_samp_total, n_avg, i_samp
    integer :: j_part, i_entry, i_part, j_id
    type(aero_particle_t) :: aero_particle_delta
    type(aero_particle_t), pointer :: i_particle
    type(aero_info_t) :: aero_info

    call assert(617651414, aero_state%aero_sorted%bin(i_bin)%n_entry > 0)
    call assert(777800964, aero_state%aero_sorted%bin(j_bin)%n_entry > 0)

    j_part = aero_state%aero_sorted%bin(j_bin)%entry(j_entry)
    weight_j = aero_weight_value(aero_weight, &
         aero_particle_radius(aero_state%p%particle(j_part)))
    j_id = aero_state%p%particle(j_part)%id
    ! don't reuse j_part later, as it may change (j_bin/j_entry are stable)
    j_part = 0

    weight_i_min = aero_weight_value(aero_weight, &
         aero_state%aero_sorted%bin_grid%edge_radius(i_bin))
    prob_remove_i_max = weight_j / weight_i_min
    call assert(653606684, prob_remove_i_max <= 1d0)

    k_max = aero_state%aero_sorted%coag_kernel_max(i_bin, j_bin)
    n_samp_mean =  k_max / aero_state%comp_vol * del_t &
         * real(aero_state%aero_sorted%bin(i_bin)%n_entry, kind=dp)
    n_samp_remove = rand_poisson(prob_remove_i_max * n_samp_mean)
    n_samp_extra = rand_poisson((1d0 - prob_remove_i_max) * n_samp_mean)
    n_samp_total = n_samp_remove + n_samp_extra

    n_avg = 0
    n_samp = 0
    prob_coag_tot = 0d0
    call aero_particle_allocate_size(aero_particle_delta, &
         aero_data%n_spec, aero_data%n_source)
    
    do i_samp = 1,n_samp_total
       if (aero_state%aero_sorted%bin(i_bin)%n_entry == 0) exit
       if ((n_samp > n_samp_remove) .and. (n_avg >= 2)) then
          vol_mean = aero_particle_delta%vol / real(n_avg, kind=dp)
          where(vol_mean > 0d0) &
               vol_cv = sqrt(vol_sq / real(n_avg, kind=dp) - (vol_mean)**2) &
               / vol_mean
          vol_cv_max = maxval(vol_cv)
          mean_95_conf_cv = vol_cv_max / sqrt(real(n_avg, kind=dp)) &
               * student_t_95_coeff(n_avg)
          ! FIXME: We are using just the max of the diagonal of the
          ! covariance matrix. Is this well-justified?
          if (mean_95_conf_cv < 0.1d0) exit
       end if
       n_samp = n_samp + 1
       ! FIXME: We are sampling with replacement. Is this a problem?
       i_entry = pmc_rand_int(aero_state%aero_sorted%bin(i_bin)%n_entry)
       i_part = aero_state%aero_sorted%bin(i_bin)%entry(i_entry)
       i_particle => aero_state%p%particle(i_part)
       ! re-get j_part as particle ordering may be changing
       j_part = aero_state%aero_sorted%bin(j_bin)%entry(j_entry)
       call weighted_kernel(coag_kernel_type, i_particle, &
            aero_state%p%particle(j_part), aero_data, aero_weight, &
            env_state, k)
       ! don't reuse j_part
       j_part = 0
       prob_coag = k / k_max
       prob_coag_tot = prob_coag_tot + prob_coag
       if (pmc_random() < prob_coag) then
          n_avg = n_avg + 1
          call aero_particle_coagulate(aero_particle_delta, i_particle, &
               aero_particle_delta)
          vol_sq = vol_sq + i_particle%vol**2
          if (i_samp <= n_samp_remove) then
             weight_i = aero_weight_value(aero_weight, &
                  aero_particle_radius(i_particle))
             prob_remove_i = weight_j / weight_i
             if (pmc_random() < prob_remove_i / prob_remove_i_max) then
                call aero_info_allocate(aero_info)
                aero_info%id = i_particle%id
                aero_info%action = AERO_INFO_COAG
                aero_info%other_id = j_id
                call aero_state_remove_particle_with_info(aero_state, &
                     i_part, aero_info)
                call aero_info_deallocate(aero_info)
             end if
          end if
       end if
    end do

    if (n_avg == 0) then
       n_coag = 0
       return
    end if

    prob_coag_mean = prob_coag_tot / real(n_avg, kind=dp)
    n_coag = rand_binomial(n_samp_total, prob_coag_mean)
    aero_particle_delta%vol = aero_particle_delta%vol &
         * (real(n_coag, kind=dp) / real(n_avg, kind=dp))
    ! re-get j_part as particle ordering may be changing
    j_part = aero_state%aero_sorted%bin(j_bin)%entry(j_entry)
    call aero_particle_coagulate(aero_state%p%particle(j_part), &
         aero_particle_delta, aero_state%p%particle(j_part))
    aero_state%p%particle(j_part)%id = j_id
    ! adjust particle number to account for weight changes
    call aero_state_reweight_particle(aero_state, aero_weight, j_part, &
         weight_j)
    ! we should only be doing this for decreasing weights
    call assert(654300924, aero_state%p%particle(j_part)%id == j_id)

  end subroutine accel_coag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the number of samples required for the pair of bins.
  subroutine compute_n_samp(ni, nj, same_bin, k_max, comp_vol, &
       del_t, n_samp_mean, n_samp, accept_factor)

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
    !> Mean number of samples per timestep.
    real(kind=dp), intent(out) :: n_samp_mean
    !> Number of samples per timestep.
    integer, intent(out) :: n_samp
    !> Scale factor for accept probability (1).
    real(kind=dp), intent(out) :: accept_factor
    
    real(kind=dp) :: r_samp
    real(kind=dp) :: n_possible ! use real(kind=dp) to avoid integer overflow
    ! could use integer*8 or integer(kind = 8)
    ! or di = selected_int_kind(18), integer(kind=di)
    ! to represent 10^{-18} to 10^{18}
    
    if (same_bin) then
       ! don't change this to ni * (ni - 1) as the ni/nj distinction
       ! is important for coagulation_dist, which also calls this
       n_possible = real(ni, kind=dp) * (real(nj, kind=dp) - 1d0) / 2d0
    else
       n_possible = real(ni, kind=dp) * real(nj, kind=dp)
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
