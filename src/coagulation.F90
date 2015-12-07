! Copyright (C) 2005-2012 Nicole Riemer
! Copyright (C) 2005-2015 Matthew West
! Copyright (C) 2015 Jeffrey Curtis
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_coagulation module.

!> Aerosol particle coagulation.
module pmc_coagulation

  use pmc_bin_grid
  use pmc_aero_data
  use pmc_util
  use pmc_stats
  use pmc_env_state
  use pmc_aero_state
  use pmc_aero_weight
  use pmc_aero_weight_array
  use pmc_mpi
  use pmc_coag_kernel
  use pmc_aero_sorted
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Minimum number of coagulation events per large particle for which
  !> accelerated coagulation is used.
  real(kind=dp), parameter :: COAG_ACCEL_N_EVENT = 1d0
  !> Maximum allowed coefficient-of-variation due to undersampling in
  !> accelerated coagulation.
  real(kind=dp), parameter :: COAG_ACCEL_MAX_CV = 0.1d0
  !> Maximum allowable growth factor of a target particle volume within one
  !> timestep when using accelerated coagulation.
  real(kind=dp), parameter :: MAX_ALLOWABLE_GROWTH_FACTOR = 1.5d0

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do coagulation for time del_t.
  subroutine mc_coag(coag_kernel_type, env_state, aero_data, aero_state, &
       del_t, tot_n_samp, tot_n_coag)

    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Timestep for coagulation.
    real(kind=dp) :: del_t
    !> Total number of samples tested.
    integer, intent(out) :: tot_n_samp
    !> Number of coagulation events.
    integer, intent(out) :: tot_n_coag

    integer :: b1, b2, c1, c2, b2_start

    call aero_state_sort(aero_state, aero_data)
    if (.not. aero_state%aero_sorted%coag_kernel_bounds_valid) then
       call est_k_minmax_binned_unweighted(aero_state%aero_sorted%bin_grid, &
            coag_kernel_type, aero_data, env_state, &
            aero_state%aero_sorted%coag_kernel_min, &
            aero_state%aero_sorted%coag_kernel_max)
       aero_state%aero_sorted%coag_kernel_bounds_valid = .true.
    end if

    tot_n_samp = 0
    tot_n_coag = 0
    do c1 = 1,aero_sorted_n_class(aero_state%aero_sorted)
       do c2 = 1,c1
          do b1 = 1,aero_sorted_n_bin(aero_state%aero_sorted)
             if (integer_varray_n_entry( &
                  aero_state%aero_sorted%size_class%inverse(b1, c1)) == 0) &
                  cycle
             if (c1 == c2) then
                b2_start = b1
             else
                b2_start = 1
             end if
             do b2 = b2_start,aero_sorted_n_bin(aero_state%aero_sorted)
                if (integer_varray_n_entry( &
                     aero_state%aero_sorted%size_class%inverse(b2, c2)) == 0) &
                     cycle
                call mc_coag_for_bin(coag_kernel_type, env_state, aero_data, &
                     aero_state, del_t, tot_n_samp, tot_n_coag, b1, b2, c1, c2)
             end do
          end do
       end do
    end do

  end subroutine mc_coag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do coagulation for time del_t for the given bins.
  subroutine mc_coag_for_bin(coag_kernel_type, env_state, aero_data, &
       aero_state, del_t, tot_n_samp, tot_n_coag, b1, b2, c1, c2)

    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Timestep for coagulation.
    real(kind=dp) :: del_t
    !> Total number of samples tested.
    integer, intent(inout) :: tot_n_samp
    !> Number of coagulation events.
    integer, intent(inout) :: tot_n_coag
    !> First bin number.
    integer, intent(in) :: b1
    !> Second bin number.
    integer, intent(in) :: b2
    !> First weight class.
    integer, intent(in) :: c1
    !> Second weight class.
    integer, intent(in) :: c2

    logical :: per_particle_coag_succeeded
    integer :: cc
    real(kind=dp) :: f_max, k_max

    cc = coag_dest_class(aero_state%awa, aero_data, &
         aero_state%aero_sorted%bin_grid, b1, b2, c1, c2)
    call max_coag_num_conc_factor(aero_state%awa, aero_data, &
         aero_state%aero_sorted%bin_grid, b1, b2, c1, c2, cc, f_max)
    k_max = aero_state%aero_sorted%coag_kernel_max(b1, b2) * f_max

    call try_per_particle_coag(coag_kernel_type, k_max, env_state, aero_data, &
         aero_state, del_t, tot_n_samp, tot_n_coag, b1, b2, c1, c2, cc, &
         per_particle_coag_succeeded)
    if (per_particle_coag_succeeded) return

    call per_set_coag(coag_kernel_type, k_max, env_state, aero_data, &
         aero_state, del_t, tot_n_samp, tot_n_coag, b1, b2, c1, c2, cc)

  end subroutine mc_coag_for_bin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Attempt per-particle coagulation.
  subroutine try_per_particle_coag(coag_kernel_type, k_max, env_state, &
       aero_data, aero_state, del_t, tot_n_samp, tot_n_coag, b1, b2, c1, c2, &
       cc, per_particle_coag_succeeded)

    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Maximum coagulation kernel (s^{-1} m^3).
    real(kind=dp), intent(in) :: k_max
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Timestep for coagulation.
    real(kind=dp) :: del_t
    !> Total number of samples tested.
    integer, intent(inout) :: tot_n_samp
    !> Number of coagulation events.
    integer, intent(inout) :: tot_n_coag
    !> First bin number.
    integer, intent(in) :: b1
    !> Second bin number.
    integer, intent(in) :: b2
    !> First weight class.
    integer, intent(in) :: c1
    !> Second weight class.
    integer, intent(in) :: c2
    !> Coagulated weight class.
    integer, intent(in) :: cc
    !> Whether we succeeded in doing per-particle coag.
    logical, intent(inout) :: per_particle_coag_succeeded

    logical :: correct_weight_ordering
    integer :: target_unif_entry, target_part, n_samp, n_coag, n_remove, bt, bs
    integer :: ct, cs
    real(kind=dp) :: n_source_per_target, accept_factor
    real(kind=dp) :: min_target_vol, max_source_vol, max_new_target_vol
    real(kind=dp) :: max_target_growth_factor
    type(aero_particle_t) :: target_particle, source_particle

    call determine_target_and_source(aero_state%awa, &
         aero_state%aero_sorted%bin_grid, b1, b2, c1, c2, cc, bt, bs, ct, cs, &
         correct_weight_ordering)
    if (.not. correct_weight_ordering) then
       per_particle_coag_succeeded = .false.
       return
    end if

    call compute_n_source(integer_varray_n_entry( &
         aero_state%aero_sorted%size_class%inverse(bs, cs)), k_max, del_t, &
         n_source_per_target, accept_factor)
    if (n_source_per_target < COAG_ACCEL_N_EVENT) then
       per_particle_coag_succeeded = .false.
       return
    end if

    min_target_vol = aero_data_rad2vol(aero_data, &
         aero_state%aero_sorted%bin_grid%edges(bt))
    max_source_vol = aero_data_rad2vol(aero_data, &
         aero_state%aero_sorted%bin_grid%edges(bs+1))
    max_new_target_vol = min_target_vol + max_source_vol * n_source_per_target
    ! check for unacceptably large volume change
    max_target_growth_factor = max_new_target_vol / min_target_vol
    if (max_target_growth_factor > MAX_ALLOWABLE_GROWTH_FACTOR) then
       per_particle_coag_succeeded = .false.
       return
    end if

    ! work backwards to avoid particle movement issues
    do target_unif_entry = integer_varray_n_entry( &
         aero_state%aero_sorted%size_class%inverse(bt, ct)),1,-1
       target_part = aero_state%aero_sorted%size_class%inverse(bt, &
            ct)%entry(target_unif_entry)
       ! need to copy coag_particle as the underlying storage may be
       ! rearranged due to removals
       target_particle = aero_state%apa%particle(target_part)
       call sample_source_particle(aero_state, aero_data, env_state, &
            coag_kernel_type, bs, cs, target_particle, n_source_per_target, &
            accept_factor, n_samp, n_coag, n_remove, source_particle)
       if (n_coag > 0) then
          call coag_target_with_source(aero_state, aero_data, bt, ct, &
               target_unif_entry, source_particle, cc)
       end if
       tot_n_samp = tot_n_samp + n_samp
       tot_n_coag = tot_n_coag + n_coag
       ! we discard n_remove information at present
    end do

    per_particle_coag_succeeded = .true.

  end subroutine try_per_particle_coag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the source and target particle bin/group for
  !> per-particle coagulation, if possible.
  subroutine determine_target_and_source(aero_weight_array, bin_grid, b1, &
       b2, c1, c2, cc, bt, bs, ct, cs, correct_weight_ordering)

    !> Aero weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> First bin number.
    integer, intent(in) :: b1
    !> Second bin number.
    integer, intent(in) :: b2
    !> First weight class.
    integer, intent(in) :: c1
    !> Second weight class.
    integer, intent(in) :: c2
    !> Coagulated weight class.
    integer, intent(in) :: cc
    !> Target bin number.
    integer, intent(out) :: bt
    !> Source bin number.
    integer, intent(out) :: bs
    !> Target weight class.
    integer, intent(out) :: ct
    !> Source weight class.
    integer, intent(out) :: cs
    !> Whether the weight ordering is correct for per-particle coagulation.
    logical, intent(out) :: correct_weight_ordering

    real(kind=dp) :: nc1_min, nc1_max, nc2_min, nc2_max
    logical :: monotone_increasing, monotone_decreasing

    call aero_weight_array_check_monotonicity(aero_weight_array, &
         monotone_increasing, monotone_decreasing)
    if (.not. monotone_decreasing) then
       correct_weight_ordering = .false.
       return
    end if

    call aero_weight_array_minmax_num_conc(aero_weight_array, c1, &
         bin_grid%edges(b1), bin_grid%edges(b1 + 1), nc1_min, nc1_max)
    call aero_weight_array_minmax_num_conc(aero_weight_array, c2, &
         bin_grid%edges(b2), bin_grid%edges(b2 + 1), nc2_min, nc2_max)

    ! we have already confirmed monotone_decreasing weights above
    correct_weight_ordering = .false.
    if ((nc1_max < nc2_min) .and. (cc == c1)) then
       bt = b1
       bs = b2
       ct = c1
       cs = c2
       correct_weight_ordering = .true.
    end if
    if ((nc2_max < nc1_min) .and. (cc == c2)) then
       bt = b2
       bs = b1
       ct = c2
       cs = c1
       correct_weight_ordering = .true.
    end if

  end subroutine determine_target_and_source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compute_n_source(n_part, k_max, del_t, n_source_per_target, &
       accept_factor)

    !> Number of particles available as partners.
    integer, intent(in) :: n_part
    !> Maximum coagulation kernel (s^{-1} m^3).
    real(kind=dp), intent(in) :: k_max
    !> Timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Mean number of source particles per target particle.
    real(kind=dp), intent(out) :: n_source_per_target
    !> Accept factor for samples.
    real(kind=dp), intent(out) :: accept_factor

    n_source_per_target = k_max * del_t * real(n_part, kind=dp)
    accept_factor = 1d0 / k_max

  end subroutine compute_n_source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sample coagulation partners for a single coagulation event.
  subroutine sample_source_particle(aero_state, aero_data, env_state, &
       coag_kernel_type, bs, cs, coag_particle, n_samp_mean, &
       accept_factor, n_samp, n_coag, n_remove, source_particle)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Bin to sample particles from.
    integer, intent(in) :: bs
    !> Weight class to sample particles from.
    integer, intent(in) :: cs
    !> Aerosol particle that coagulation will be with.
    type(aero_particle_t), intent(in) :: coag_particle
    !> Mean number of samples to use.
    real(kind=dp), intent(in) :: n_samp_mean
    !> Probability factor of accepting samples.
    real(kind=dp), intent(in) :: accept_factor
    !> Number of samples used.
    integer, intent(out) :: n_samp
    !> Number of coagulations performed.
    integer, intent(out) :: n_coag
    !> Number of particles removed.
    integer, intent(out) :: n_remove
    !> Sampled average coagulation partner particle.
    type(aero_particle_t), intent(inout) :: source_particle

    real(kind=dp) :: prob_remove_i, prob_remove_source_max
    real(kind=dp) :: prob_coag, prob_coag_tot, prob_coag_mean
    real(kind=dp) :: num_conc_i, num_conc_source_min, num_conc_target, k
    real(kind=dp) :: vol_sq(aero_data_n_spec(aero_data))
    real(kind=dp) :: vol_mean(aero_data_n_spec(aero_data))
    real(kind=dp) :: vol_cv(aero_data_n_spec(aero_data)), vol_cv_max
    real(kind=dp) :: mean_95_conf_cv
    integer :: n_samp_remove, n_samp_extra, n_samp_total, n_avg, i_samp
    integer :: i_unif_entry, i_part, target_id, new_bin, ct
    type(aero_info_t) :: aero_info

    if (integer_varray_n_entry( &
         aero_state%aero_sorted%size_class%inverse(bs, cs)) == 0) then
       n_samp = 0
       n_remove = 0
       n_coag = 0
       return
    end if

    num_conc_target = aero_weight_array_num_conc(aero_state%awa, &
         coag_particle, aero_data)
    target_id = coag_particle%id
    ct = coag_particle%weight_class

    num_conc_source_min = aero_weight_array_num_conc_at_radius( &
         aero_state%awa, cs, aero_state%aero_sorted%bin_grid%edges(bs))
    prob_remove_source_max = num_conc_target / num_conc_source_min
    call assert(653606684, prob_remove_source_max <= 1d0)

    n_samp_remove = rand_poisson(prob_remove_source_max * n_samp_mean)
    n_samp_extra = rand_poisson((1d0 - prob_remove_source_max) * n_samp_mean)
    n_samp_total = n_samp_remove + n_samp_extra

    n_avg = 0
    n_samp = 0
    n_remove = 0
    prob_coag_tot = 0d0
    vol_sq = 0d0

    call aero_particle_zero(source_particle, aero_data)

    ! FIXME: Can't we just do n_samp = 1,n_samp_total and shift tests
    ! to the end?
    do i_samp = 1,n_samp_total
       if (integer_varray_n_entry( &
            aero_state%aero_sorted%size_class%inverse(bs, cs)) == 0) exit
       if ((n_samp > n_samp_remove) .and. (n_avg >= 2)) then
          vol_mean = source_particle%vol / real(n_avg, kind=dp)
          where(vol_mean > 0d0) &
               vol_cv = sqrt(vol_sq / real(n_avg, kind=dp) - (vol_mean)**2) &
               / vol_mean
          vol_cv_max = maxval(vol_cv)
          mean_95_conf_cv = vol_cv_max / sqrt(real(n_avg, kind=dp)) &
               * student_t_95_coeff(n_avg)
          ! FIXME: We are using just the max of the diagonal of the
          ! covariance matrix. Is this well-justified?
          if (mean_95_conf_cv < COAG_ACCEL_MAX_CV) exit
       end if
       n_samp = n_samp + 1
       ! FIXME: We are sampling with replacement. Is this a problem?
       i_unif_entry = pmc_rand_int(integer_varray_n_entry( &
            aero_state%aero_sorted%size_class%inverse(bs, cs)))
       i_part = aero_state%aero_sorted%size_class%inverse(bs, &
            cs)%entry(i_unif_entry)
       ! re-get j_part as particle ordering may be changing
       call num_conc_weighted_kernel(coag_kernel_type, &
            aero_state%apa%particle(i_part), coag_particle, cs, ct, ct, &
            aero_data, aero_state%awa, env_state, k)
       prob_coag = k * accept_factor
       prob_coag_tot = prob_coag_tot + prob_coag
       if (pmc_random() < prob_coag) then
          n_avg = n_avg + 1
          call aero_particle_coagulate(source_particle, &
               aero_state%apa%particle(i_part), source_particle)
          vol_sq = vol_sq + aero_state%apa%particle(i_part)%vol**2
          if (i_samp <= n_samp_remove) then
             num_conc_i = aero_weight_array_num_conc(aero_state%awa, &
                  aero_state%apa%particle(i_part), aero_data)
             prob_remove_i = num_conc_target / num_conc_i
             if (pmc_random() < prob_remove_i / prob_remove_source_max) then
                n_remove = n_remove + 1
                aero_info%id = aero_state%apa%particle(i_part)%id
                aero_info%action = AERO_INFO_COAG
                aero_info%other_id = target_id
                call aero_state_remove_particle_with_info(aero_state, &
                     i_part, aero_info)
             end if
          end if
       end if
    end do

    if (n_avg == 0) then
       n_coag = 0
       return
    end if

    prob_coag_mean = prob_coag_tot / real(n_samp, kind=dp) ! note not n_avg
    call warn_assert_msg(383983415, prob_coag_mean <= 1d0, &
         "kernel upper bound estimation is too tight")
    n_coag = rand_binomial(n_samp_total, prob_coag_mean)
    source_particle%vol = source_particle%vol &
         * (real(n_coag, kind=dp) / real(n_avg, kind=dp))

  end subroutine sample_source_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Coagulate a sampled source particle with a target particle.
  subroutine coag_target_with_source(aero_state, aero_data, bt, ct, &
       target_unif_entry, source_particle, cc)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Bin of coagulating particle.
    integer, intent(in) :: bt
    !> Weight class of coagulating particle.
    integer, intent(in) :: ct
    !> Entry-in-bin of coagulating particle.
    integer, intent(in) :: target_unif_entry
    !> Sampled particle to coagulate with.
    type(aero_particle_t), intent(in) :: source_particle
    !> Weight class for coagulated particle.
    integer, intent(in) :: cc

    integer :: target_part, target_id, new_bin, new_group
    real(kind=dp) :: old_num_conc_target, new_num_conc_target

    target_part = aero_state%aero_sorted%size_class%inverse(bt, &
         ct)%entry(target_unif_entry)
    target_id = aero_state%apa%particle(target_part)%id
    old_num_conc_target = aero_weight_array_num_conc(aero_state%awa, &
         aero_state%apa%particle(target_part), aero_data)
    call aero_particle_coagulate(aero_state%apa%particle(target_part), &
         source_particle, aero_state%apa%particle(target_part))
    aero_state%apa%particle(target_part)%id = target_id
    ! assign to a randomly chosen group
    new_group = aero_weight_array_rand_group(aero_state%awa, cc, &
         aero_particle_radius(aero_state%apa%particle(target_part), &
         aero_data))
    call aero_particle_set_weight(aero_state%apa%particle(target_part), &
         new_group, cc)
    ! fix bin due to composition changes
    new_bin = aero_sorted_particle_in_bin(aero_state%aero_sorted, &
         aero_state%apa%particle(target_part), aero_data)
    if ((new_bin < 1) &
         .or. (new_bin > bin_grid_size(aero_state%aero_sorted%bin_grid))) then
       call die_msg(765620746, "particle outside of bin_grid: " &
            // "try reducing the timestep del_t")
    end if
    call aero_sorted_move_particle(aero_state%aero_sorted, target_part, &
         new_bin, new_group, cc)
    ! Adjust particle number to account for weight changes
    ! bt/target_group/target_entry are invalid, but
    ! target_part is still good. We are treating all particles in all
    ! groups together, and randomly reassigning between groups above,
    ! so here we can't use aero_state_reweight_particle(), as that
    ! assumes we are staying in the same weight group.
    new_num_conc_target = aero_weight_array_num_conc(aero_state%awa, &
         aero_state%apa%particle(target_part), aero_data)
    call aero_state_dup_particle(aero_state, aero_data, target_part, &
         old_num_conc_target / new_num_conc_target, random_weight_group=.true.)
    ! we should only be doing this for decreasing weights
    call assert(654300924, aero_state%apa%particle(target_part)%id &
         == target_id)

  end subroutine coag_target_with_source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do set-wise coagulation.
  subroutine per_set_coag(coag_kernel_type, k_max, env_state, aero_data, &
       aero_state, del_t, tot_n_samp, tot_n_coag, b1, b2, c1, c2, cc)

    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Maximum coagulation kernel (s^{-1} m^3).
    real(kind=dp), intent(in) :: k_max
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Timestep for coagulation.
    real(kind=dp) :: del_t
    !> Total number of samples tested.
    integer, intent(inout) :: tot_n_samp
    !> Number of coagulation events.
    integer, intent(inout) :: tot_n_coag
    !> First bin number.
    integer, intent(in) :: b1
    !> Second bin number.
    integer, intent(in) :: b2
    !> First weight class.
    integer, intent(in) :: c1
    !> Second weight class.
    integer, intent(in) :: c2
    !> Coagulated weight class.
    integer, intent(in) :: cc

    real(kind=dp) :: n_samp_mean, accept_factor
    integer :: i_samp, n_samp, n1, n2
    logical :: did_coag

    n1 = integer_varray_n_entry( &
         aero_state%aero_sorted%size_class%inverse(b1, c1))
    n2 = integer_varray_n_entry( &
         aero_state%aero_sorted%size_class%inverse(b2, c2))
    call compute_n_samp(n1, n2, ((b1 == b2) .and. (c1 == c2)), k_max, del_t, &
         n_samp_mean, n_samp, accept_factor)
    tot_n_samp = tot_n_samp + n_samp

    do i_samp = 1,n_samp
       ! check we still have enough particles to coagulate
       n1 = integer_varray_n_entry( &
            aero_state%aero_sorted%size_class%inverse(b1, c1))
       n2 = integer_varray_n_entry( &
            aero_state%aero_sorted%size_class%inverse(b2, c2))
       if (((n1 < 2) .and. (b1 == b2) .and. (c1 == c2)) &
            .or. (n1 < 1) .or. (n2 < 1)) &
            exit
       call maybe_coag_pair(env_state, aero_data, aero_state, b1, b2, c1, c2, &
            cc, coag_kernel_type, accept_factor, did_coag)
       if (did_coag) tot_n_coag = tot_n_coag + 1
    end do

  end subroutine per_set_coag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the number of samples required for the pair of bins.
  subroutine compute_n_samp(ni, nj, same_bin, k_max, del_t, n_samp_mean, &
       n_samp, accept_factor)

    !> Number particles in first bin.
    integer, intent(in) :: ni
    !> Number particles in second bin.
    integer, intent(in) :: nj
    !> Whether first bin is second bin.
    logical, intent(in) :: same_bin
    !> Maximum kernel value (s^{-1} m^3).
    real(kind=dp), intent(in) :: k_max
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

    r_samp = k_max * del_t
    n_samp_mean = r_samp * n_possible
    n_samp = rand_poisson(n_samp_mean)
    accept_factor = 1d0 / k_max

    ! possible variants:
    ! A: accept_factor = 1d0 / k_max
    ! B: accept_factor = del_t * n_possible / (n_samp * V)
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
  subroutine maybe_coag_pair(env_state, aero_data, aero_state, b1, b2, &
       c1, c2, cc, coag_kernel_type, accept_factor, did_coag)

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
    !> First weight class.
    integer, intent(in) :: c1
    !> Second weight class.
    integer, intent(in) :: c2
    !> Coagulated weight class.
    integer, intent(in) :: cc
    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Scale factor for accept probability (1).
    real(kind=dp), intent(in) :: accept_factor
    !> Whether a coagulation occured.
    logical, intent(out) :: did_coag

    integer :: i1, i2
    integer :: p1, p2
    real(kind=dp) :: p, k

    call find_rand_pair(aero_state%aero_sorted, b1, b2, c1, c2, i1, i2)
    p1 = aero_state%aero_sorted%size_class%inverse(b1, c1)%entry(i1)
    p2 = aero_state%aero_sorted%size_class%inverse(b2, c2)%entry(i2)
    call num_conc_weighted_kernel(coag_kernel_type, &
         aero_state%apa%particle(p1), aero_state%apa%particle(p2), &
         c1, c2, cc, aero_data, aero_state%awa, env_state, k)
    p = k * accept_factor
    call warn_assert_msg(532446093, p <= 1d0, &
         "kernel upper bound estimation is too tight, " &
         // "could be caused by changing env_state")

    did_coag = .false.
    if (pmc_random() .lt. p) then
       call coagulate(aero_data, aero_state, p1, p2, c1, c2, cc)
       did_coag = .true.
    end if

  end subroutine maybe_coag_pair

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Given <tt>(b1, c1)</tt> and <tt>(b2, c2)</tt>, find a random pair
  !> of particles <tt>(b1, c1, i1)</tt> and <tt>(b2, c2, i2)</tt> that
  !> are not the same particle particle as each other.
  subroutine find_rand_pair(aero_sorted, b1, b2, c1, c2, i1, i2)

    !> Aerosol sorted data.
    type(aero_sorted_t), intent(in) :: aero_sorted
    !> Bin number of first particle.
    integer, intent(in) :: b1
    !> Bin number of second particle.
    integer, intent(in) :: b2
    !> Weight class of first particle.
    integer, intent(in) :: c1
    !> Weight class of second particle.
    integer, intent(in) :: c2
    !> First rand particle.
    integer, intent(out) :: i1
    !> Second rand particle.
    integer, intent(out) :: i2

    call assert(619608562, &
         integer_varray_n_entry(aero_sorted%size_class%inverse(b1, c1)) >= 1)
    i1 = pmc_rand_int( &
         integer_varray_n_entry(aero_sorted%size_class%inverse(b1, c1)))

    if ((b1 == b2) .and. (c1 == c2)) then
       call assert(956184336, integer_varray_n_entry( &
            aero_sorted%size_class%inverse(b2, c2)) >= 2)
       i2 = pmc_rand_int( &
            integer_varray_n_entry(aero_sorted%size_class%inverse(b2, c2)) - 1)
       if (i2 == i1) then
          i2 = integer_varray_n_entry(aero_sorted%size_class%inverse(b2, c2))
       end if
    else
       call assert(271635751, integer_varray_n_entry( &
            aero_sorted%size_class%inverse(b2, c2)) >= 1)
       i2 = pmc_rand_int( &
            integer_varray_n_entry(aero_sorted%size_class%inverse(b2, c2)))
    end if

  end subroutine find_rand_pair

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Actually coagulate pt1 and pt2 to form ptc and compute weighting
  !> effects, including which particles should be lost and which
  !> gained.
  subroutine coagulate_weighting(pt1, pt2, ptc, c1, c2, cc, aero_data, &
       aero_weight_array, remove_1, remove_2, create_new, id_1_lost, &
       id_2_lost, aero_info_1, aero_info_2)

    !> First coagulating aerosol particle.
    type(aero_particle_t), intent(in) :: pt1
    !> Second coagulating aerosol particle.
    type(aero_particle_t), intent(in) :: pt2
    !> Combined aerosol particle resulting from coagulation of pt1
    !> and pt2.
    type(aero_particle_t), intent(inout) :: ptc
    !> First weight class.
    integer, intent(in) :: c1
    !> Second weight class.
    integer, intent(in) :: c2
    !> Coagulated weight class.
    integer, intent(in) :: cc
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Whether to remove pt1.
    logical, intent(out) :: remove_1
    !> Whether to remove pt2.
    logical, intent(out) :: remove_2
    !> Whether to create ptc.
    logical, intent(out) :: create_new
    !> Whether the ID of pt1 will be lost due to coagulation.
    logical, intent(out) :: id_1_lost
    !> Whether the ID of pt2 will be lost due to coagulation.
    logical, intent(out) :: id_2_lost
    !> The removal information associated with pt1.
    type(aero_info_t), intent(inout) :: aero_info_1
    !> The removal information associated with pt2.
    type(aero_info_t), intent(inout) :: aero_info_2

    real(kind=dp) :: r1, r2, rc, nc_min, nc1, nc2, ncc
    real(kind=dp) :: prob_remove_1, prob_remove_2, prob_create_new
    integer :: info_other_id, new_group

    call assert(371947172, pt1%id /= pt2%id)

    ! decide which old particles are to be removed and whether to
    ! create the resulting coagulated particle
    r1 = aero_particle_radius(pt1, aero_data)
    r2 = aero_particle_radius(pt2, aero_data)
    rc = aero_data_vol2rad(aero_data, aero_data_rad2vol(aero_data, r1) &
         + aero_data_rad2vol(aero_data, r2))
    nc1 = aero_weight_array_num_conc_at_radius(aero_weight_array, c1, r1)
    nc2 = aero_weight_array_num_conc_at_radius(aero_weight_array, c2, r2)
    ncc = aero_weight_array_num_conc_at_radius(aero_weight_array, cc, rc)
    new_group = aero_weight_array_rand_group(aero_weight_array, cc, rc)
    nc_min = min(nc1, nc2, ncc)
    prob_remove_1 = nc_min / nc1
    prob_remove_2 = nc_min / nc2
    prob_create_new = nc_min / ncc
    remove_1 = (pmc_random() < prob_remove_1)
    ! FIXME
    !if (aero_weight%type == AERO_WEIGHT_TYPE_MFA) then
    !   remove_2 = .not. remove_1
    !else
    remove_2 = (pmc_random() < prob_remove_2)
    !end if
    create_new = (pmc_random() < prob_create_new)

    ! figure out what to do about the ID numbers of the various
    ! particles --- we try to preserve particle IDs as much as
    ! possible
    if (create_new) then
       id_1_lost = .false.
       id_2_lost = .false.
       if (remove_1 .and. remove_2) then
          if (aero_particle_volume(pt1) > aero_particle_volume(pt2)) then
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
       call aero_particle_coagulate(pt1, pt2, ptc)
       call aero_particle_set_weight(ptc, new_group, cc)
       if (remove_1 .and. (.not. id_1_lost)) then
          ptc%id = pt1%id
          call assert(975059559, id_2_lost .eqv. remove_2)
       elseif (remove_2 .and. (.not. id_2_lost)) then
          ptc%id = pt2%id
          call assert(246529753, id_1_lost .eqv. remove_1)
       else
          call aero_particle_new_id(ptc)
          call assert(852038606, id_1_lost .eqv. remove_1)
          call assert(254018921, id_2_lost .eqv. remove_2)
       end if
       info_other_id = ptc%id
    else
       info_other_id = 0
    end if

    aero_info_1%id = pt1%id
    aero_info_1%action = AERO_INFO_COAG
    aero_info_1%other_id = info_other_id

    aero_info_2%id = pt2%id
    aero_info_2%action = AERO_INFO_COAG
    aero_info_2%other_id = info_other_id

  end subroutine coagulate_weighting

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Join together particles \c p1 and \c p2, updating all particle
  !> and bin structures to reflect the change.
  subroutine coagulate(aero_data, aero_state, p1, p2, c1, c2, cc)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> First particle index.
    integer, intent(in) :: p1
    !> Second particle index.
    integer, intent(in) :: p2
    !> First weight class.
    integer, intent(in) :: c1
    !> Second weight class.
    integer, intent(in) :: c2
    !> Coagulated weight class.
    integer, intent(in) :: cc

    type(aero_particle_t) :: ptc
    integer :: bn
    type(aero_info_t) :: aero_info_1, aero_info_2
    logical :: remove_1, remove_2, create_new, id_1_lost, id_2_lost

    call coagulate_weighting(aero_state%apa%particle(p1), &
         aero_state%apa%particle(p2), ptc, c1, c2, cc, aero_data, &
         aero_state%awa, remove_1, remove_2, create_new, id_1_lost, &
         id_2_lost, aero_info_1, aero_info_2)

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
       call aero_state_add_particle(aero_state, ptc, aero_data, &
            allow_resort=.false.)
    end if

  end subroutine coagulate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_coagulation
