! Copyright (C) 2005-2015 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for the.

!> \file
!> The pmc_aero_binned module.

!> The aero_binned_t structure and associated subroutines.
module pmc_aero_binned

  use pmc_bin_grid
  use pmc_aero_particle
  use pmc_spec_file
  use pmc_util
  use pmc_bin_grid
  use pmc_aero_dist
  use pmc_mpi
  use pmc_aero_data
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Aerosol number and volume distributions stored per bin.
  !!
  !! These quantities are densities both in volume (per m^3) and in
  !! radius (per log_width). The total concentration per volume is computed as
  !! sum(aero_binned%%num_conc * bin_grid%%log_width).
  !!
  !! An aero_binned_t is similar to an aero_dist_t in that they both
  !! store binned aerosol distributions. The difference is that an
  !! aero_dist_t has the same composition in every bin, whereas an
  !! aero_binned_t can have aerosol composition that varies per bin.
  !!
  !! By convention, if aero_binned_is_allocated() return \c .false.,
  !! then the aero_binned is treated as zero for all operations on
  !! it. This will be the case for new \c aero_binned_t structures.
  type aero_binned_t
     !> Number concentration per bin (#/m^3/log_width).
     !! Array length is \c bin_grid_size(bin_grid).
     real(kind=dp), allocatable :: num_conc(:)
     !> Volume concentration per bin and per species (m^3/m^3/log_width).
     !! Array size is <tt>bin_grid_size(bin_grid) x
     !! aero_data_n_spec(aero_data)</tt>.
     real(kind=dp), allocatable :: vol_conc(:,:)
  end type aero_binned_t

  !> Moving-center bin remap: whole-bin moves to the destination bin.
  integer, parameter :: AERO_BINNED_REDIST_MOVING_CENTER = 1
  !> Two-moment (linear-discrete) bin remap: splits across two bins.
  integer, parameter :: AERO_BINNED_REDIST_TWO_MOMENT = 2

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine whether the \c aero_binned is correctly allocated.
  logical function aero_binned_is_allocated(aero_binned)

    !> Aerosol binned to check.
    type(aero_binned_t), intent(in) :: aero_binned

    logical :: valid

    valid = .true.
    valid = valid .and. allocated(aero_binned%num_conc)
    valid = valid .and. allocated(aero_binned%vol_conc)
    valid = valid &
         .and. (size(aero_binned%num_conc) == size(aero_binned%num_conc, 1))
    aero_binned_is_allocated = valid

  end function aero_binned_is_allocated

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the number of bins and species in an aero_binned_t.
  subroutine aero_binned_set_sizes(aero_binned, n_bin, n_spec)

    !> Structure to be allocated.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Number of aerosol bins to allocate (typically \c bin_grid%%n_bin).
    integer, intent(in) :: n_bin
    !> Number of aerosol species to allocate (typically
    !> \c aero_data%%n_spec).
    integer, intent(in) :: n_spec

    if (allocated(aero_binned%num_conc)) deallocate(aero_binned%num_conc)
    if (allocated(aero_binned%vol_conc)) deallocate(aero_binned%vol_conc)
    allocate(aero_binned%num_conc(n_bin))
    allocate(aero_binned%vol_conc(n_bin, n_spec))
    call aero_binned_zero(aero_binned)

  end subroutine aero_binned_set_sizes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set all internal data in an aero_binned_t structure to zero.
  subroutine aero_binned_zero(aero_binned)

    !> Structure to zero.
    type(aero_binned_t), intent(inout) :: aero_binned

    if (aero_binned_is_allocated(aero_binned)) then
       aero_binned%num_conc = 0d0
       aero_binned%vol_conc = 0d0
    end if

  end subroutine aero_binned_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add two aero_binned_t structures together.
  !!
  !! Symbolically does aero_binned = aero_binned + aero_binned_delta.
  subroutine aero_binned_add(aero_binned, aero_binned_delta)

    !> Base aero_binned_t structure that will be added to.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Structure to add to aero_binned.
    type(aero_binned_t), intent(in) :: aero_binned_delta

    if (aero_binned_is_allocated(aero_binned_delta)) then
       if (aero_binned_is_allocated(aero_binned)) then
          aero_binned%num_conc = aero_binned%num_conc &
               + aero_binned_delta%num_conc
          aero_binned%vol_conc = aero_binned%vol_conc &
               + aero_binned_delta%vol_conc
       else
          aero_binned%num_conc = aero_binned_delta%num_conc
          aero_binned%vol_conc = aero_binned_delta%vol_conc
       end if
    end if

  end subroutine aero_binned_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a scaled \c aero_binned_t structure to an existing one.
  !!
  !! Symbolically does aero_binned = aero_binned + alpha * aero_binned_delta.
  subroutine aero_binned_add_scaled(aero_binned, aero_binned_delta, alpha)

    !> Base aero_binned_t structure that will be added to.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Structure to add to aero_binned.
    type(aero_binned_t), intent(in) :: aero_binned_delta
    !> Scale factor.
    real(kind=dp), intent(in) :: alpha

    if (aero_binned_is_allocated(aero_binned_delta)) then
       if (aero_binned_is_allocated(aero_binned)) then
          aero_binned%num_conc = aero_binned%num_conc &
               + alpha * aero_binned_delta%num_conc
          aero_binned%vol_conc = aero_binned%vol_conc &
               + alpha * aero_binned_delta%vol_conc
       else
          aero_binned%num_conc = aero_binned_delta%num_conc
          aero_binned%vol_conc = aero_binned_delta%vol_conc
       end if
    end if

  end subroutine aero_binned_add_scaled

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Subtract one aero_binned_t structure from another.
  !!
  !! Symbolically does aero_binned = aero_binned - aero_binned_delta.
  subroutine aero_binned_sub(aero_binned, aero_binned_delta)

    !> Base aero_binned_t structure that will be subtracted from.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Structure to subtract from aero_binned.
    type(aero_binned_t), intent(in) :: aero_binned_delta

    if (aero_binned_is_allocated(aero_binned_delta)) then
       if (aero_binned_is_allocated(aero_binned)) then
          aero_binned%num_conc = aero_binned%num_conc &
               - aero_binned_delta%num_conc
          aero_binned%vol_conc = aero_binned%vol_conc &
               - aero_binned_delta%vol_conc
       else
          aero_binned%num_conc = - aero_binned_delta%num_conc
          aero_binned%vol_conc = - aero_binned_delta%vol_conc
       end if
    end if

  end subroutine aero_binned_sub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scale an aero_binned_t by a real number.
  !!
  !! Symbolically does aero_binned = aero_binned * alpha.
  subroutine aero_binned_scale(aero_binned, alpha)

    !> Base aero_binned to scale.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Scale factor.
    real(kind=dp), intent(in) :: alpha

    if (aero_binned_is_allocated(aero_binned)) then
       aero_binned%num_conc = aero_binned%num_conc * alpha
       aero_binned%vol_conc = aero_binned%vol_conc * alpha
    end if

  end subroutine aero_binned_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scales an aero_binned_t element-wise by an array of reals.
  subroutine aero_binned_scale_by_array(aero_binned, alpha_array)

    !> Base aero_binned_t structure that will be scaled.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Structure to scale aero_binned.
    real(kind=dp), allocatable, intent(in) :: alpha_array(:)

    integer :: i

    do i = 1, size(aero_binned%num_conc)
       aero_binned%num_conc(i) = alpha_array(i)*aero_binned%num_conc(i)
       aero_binned%vol_conc(i,:) = alpha_array(i)*aero_binned%vol_conc(i,:)
    end do

  end subroutine aero_binned_scale_by_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add an aero_dist_t to an aero_binned_t.
  !!
  !! Symbolically does aero_binned = aero_binned + aero_dist.
  subroutine aero_binned_add_aero_dist(aero_binned, bin_grid, &
       aero_data, aero_dist)

    !> Base aero_binned_t structure to add to.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> The aero_dist_t structure to add.
    type(aero_dist_t), intent(in) :: aero_dist

    real(kind=dp) :: dist_num_conc(bin_grid_size(bin_grid))
    real(kind=dp) :: dist_vol_conc(bin_grid_size(bin_grid), &
         aero_data_n_spec(aero_data))

    call aero_dist_num_conc(aero_dist, bin_grid, aero_data, &
         dist_num_conc)
    call aero_dist_vol_conc(aero_dist, bin_grid, aero_data, &
         dist_vol_conc)
    if (aero_binned_is_allocated(aero_binned)) then
       aero_binned%num_conc = aero_binned%num_conc + dist_num_conc
       aero_binned%vol_conc = aero_binned%vol_conc + dist_vol_conc
    else
       aero_binned%num_conc = dist_num_conc
       aero_binned%vol_conc = dist_vol_conc
    end if

  end subroutine aero_binned_add_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remap aerosol back onto the fixed bin grid after growth or
  !> shrinkage, using the moving-center scheme.
  !!
  !! Each bin's entire number and per-species volume are moved to the bin
  !! that contains the bin's current number-mean particle volume (computed
  !! from the bin totals). This is the operator-split companion to a
  !! condensation/evaporation step (e.g. TChem) that changes the per-bin
  !! masses while holding the bin number concentrations fixed, so that the
  !! mean particle size no longer matches the bin it sits in.
  !!
  !! Total number concentration and total per-species volume concentration
  !! are both conserved exactly. A bin whose mean falls below the first bin
  !! or above the last bin is clamped into the first/last bin.
  !!
  !! This is the moving-center variant (whole-bin moves); see
  !! aero_binned_redistribute_two_moment() for the lower-diffusion
  !! linear-discrete alternative.
  subroutine aero_binned_redistribute_moving_center(aero_binned, bin_grid, &
       aero_data)

    !> Binned aerosol distribution to redistribute in place.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Bin grid (radius-based).
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol material data.
    type(aero_data_t), intent(in) :: aero_data

    integer :: n_bin, n_spec, i_bin, i_new
    real(kind=dp) :: num_actual, total_vol, vol_mean, rad_mean
    real(kind=dp) :: vol_actual(aero_data_n_spec(aero_data))
    real(kind=dp) :: new_num(bin_grid_size(bin_grid))
    real(kind=dp) :: new_vol(bin_grid_size(bin_grid), &
         aero_data_n_spec(aero_data))

    if (.not. aero_binned_is_allocated(aero_binned)) return

    n_bin = bin_grid_size(bin_grid)
    n_spec = aero_data_n_spec(aero_data)
    if (n_bin < 1) return

    ! Accumulate into actual concentrations (#/m^3 and m^3/m^3), i.e. the
    ! per-log-width densities multiplied by the bin widths, so that moving
    ! material between bins of (possibly) different width conserves the
    ! totals rather than the densities.
    new_num = 0d0
    new_vol = 0d0

    do i_bin = 1,n_bin
       num_actual = aero_binned%num_conc(i_bin) * bin_grid%widths(i_bin)
       vol_actual = aero_binned%vol_conc(i_bin,:) * bin_grid%widths(i_bin)
       total_vol = sum(vol_actual)

       if ((num_actual > 0d0) .and. (total_vol > 0d0)) then
          ! number-mean single-particle volume, then the bin containing it
          vol_mean = total_vol / num_actual
          rad_mean = aero_data_vol2rad(aero_data, vol_mean)
          i_new = bin_grid_find(bin_grid, rad_mean)
          i_new = max(1, min(n_bin, i_new))
       else
          ! empty or mass-free bin: leave its contents in place
          i_new = i_bin
       end if

       new_num(i_new) = new_num(i_new) + num_actual
       new_vol(i_new,:) = new_vol(i_new,:) + vol_actual
    end do

    ! convert the actual concentrations back to per-log-width densities
    do i_bin = 1,n_bin
       aero_binned%num_conc(i_bin) = new_num(i_bin) / bin_grid%widths(i_bin)
       aero_binned%vol_conc(i_bin,:) = new_vol(i_bin,:) / bin_grid%widths(i_bin)
    end do

  end subroutine aero_binned_redistribute_moving_center

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remap aerosol back onto the fixed bin grid after growth, using the
  !> two-moment (linear-discrete) scheme of Jacobson / MOSAIC move_sections.
  !!
  !! Unlike the moving-center variant (whole-bin moves, which leaves a
  !! comb of spikes), this reconstructs a linear sub-bin number distribution
  !! in each source bin from its pre-growth number and mean volume, then
  !! splits the grown number and volume *fractionally* between the
  !! destination bin and one neighbour -- conserving both number and volume
  !! while greatly reducing the numerical-diffusion / spiking artifact.
  !!
  !! \c aero_binned holds the post-growth (aftgrow) state and is updated in
  !! place; \c aero_binned_pregrow holds the pre-growth state (captured before
  !! the condensation step). Number is unchanged by condensation, so the two
  !! share the same per-bin number; only the volumes differ.
  subroutine aero_binned_redistribute_two_moment(aero_binned, &
       aero_binned_pregrow, bin_grid, aero_data)

    !> Post-growth binned distribution, redistributed in place.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Pre-growth binned distribution (same number, smaller volumes).
    type(aero_binned_t), intent(in) :: aero_binned_pregrow
    !> Bin grid (radius-based).
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol material data.
    type(aero_data_t), intent(in) :: aero_data

    integer :: n_bin, i_bin, dest1, dest2
    real(kind=dp) :: num_actual, vtot_pre, vtot_aft, frac_num, frac_vol
    real(kind=dp) :: vol_actual(aero_data_n_spec(aero_data))
    real(kind=dp) :: vol_edge(bin_grid_size(bin_grid) + 1)
    real(kind=dp) :: new_num(bin_grid_size(bin_grid))
    real(kind=dp) :: new_vol(bin_grid_size(bin_grid), &
         aero_data_n_spec(aero_data))

    if (.not. aero_binned_is_allocated(aero_binned)) return
    n_bin = bin_grid_size(bin_grid)
    if (n_bin < 1) return

    ! single-particle dry volume at each bin edge
    vol_edge = aero_data_rad2vol(aero_data, bin_grid%edges)

    new_num = 0d0
    new_vol = 0d0

    do i_bin = 1,n_bin
       num_actual = aero_binned%num_conc(i_bin) * bin_grid%widths(i_bin)
       vol_actual = aero_binned%vol_conc(i_bin,:) * bin_grid%widths(i_bin)
       vtot_aft = sum(aero_binned%vol_conc(i_bin,:))
       vtot_pre = sum(aero_binned_pregrow%vol_conc(i_bin,:))

       ! per-log-width densities cancel inside the split (it uses ratios), so
       ! pass them directly; deposit the actual (x width) amounts below.
       call aero_binned_two_moment_split(aero_binned%num_conc(i_bin), &
            vtot_pre, vtot_aft, i_bin, n_bin, vol_edge, dest1, dest2, &
            frac_num, frac_vol)

       new_num(dest1) = new_num(dest1) + num_actual * frac_num
       new_vol(dest1,:) = new_vol(dest1,:) + vol_actual * frac_vol
       if (dest2 > 0) then
          new_num(dest2) = new_num(dest2) + num_actual * (1d0 - frac_num)
          new_vol(dest2,:) = new_vol(dest2,:) + vol_actual * (1d0 - frac_vol)
       end if
    end do

    do i_bin = 1,n_bin
       aero_binned%num_conc(i_bin) = new_num(i_bin) / bin_grid%widths(i_bin)
       aero_binned%vol_conc(i_bin,:) = new_vol(i_bin,:) / bin_grid%widths(i_bin)
    end do

  end subroutine aero_binned_redistribute_two_moment

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> For one source bin, determine the destination bin(s) and the number /
  !> volume fractions for the two-moment remap.
  !!
  !! Ports the linear-discrete logic of MOSAIC \c move_sections: reconstruct a
  !! linear sub-bin distribution \c n(z) = aa + bb*z from the pre-growth number
  !! and mean volume, find the destination bin \c dest1 of the grown mean,
  !! map that bin's edges back into pre-growth volume space, and integrate the
  !! 0th/1st moments over the overlap to get the fraction staying in \c dest1
  !! (the rest goes to neighbour \c dest2). Falls back to a whole-bin move
  !! (\c dest2 = 0, fractions = 1) whenever the reconstruction is degenerate.
  subroutine aero_binned_two_moment_split(num, vtot_pre, vtot_aft, n, n_bin, &
       vol_edge, dest1, dest2, frac_num, frac_vol)

    !> Bin number (any consistent units; ratios cancel).
    real(kind=dp), intent(in) :: num
    !> Pre-growth and post-growth total volume in the bin (same units as num).
    real(kind=dp), intent(in) :: vtot_pre, vtot_aft
    !> Source bin index and number of bins.
    integer, intent(in) :: n, n_bin
    !> Single-particle dry volume at the bin edges (m^3), length n_bin+1.
    real(kind=dp), intent(in) :: vol_edge(n_bin + 1)
    !> Primary destination bin and secondary (neighbour) bin (0 if none).
    integer, intent(out) :: dest1, dest2
    !> Fraction of number and of volume going to dest1 (remainder to dest2).
    real(kind=dp), intent(out) :: frac_num, frac_vol

    integer :: nnew, nnew2
    real(kind=dp) :: vbar_aft, vbar_pre, vlo, vhi, vdel, gamma, ratio, aa, bb
    real(kind=dp) :: vtmp, vcutlo, vcuthi, zlo, zhi, d1, d2, d3

    ! default: whole-bin move (moving-center fallback)
    dest2 = 0
    frac_num = 1d0
    frac_vol = 1d0

    ! negligible / empty bin -> leave in place
    if ((num <= 0d0) .or. (vtot_aft <= 0d0)) then
       dest1 = n
       return
    end if

    ! destination bin of the grown number-mean particle volume
    vbar_aft = vtot_aft / num
    if (vbar_aft >= vol_edge(n_bin + 1)) then
       dest1 = n_bin
       return
    else if (vbar_aft <= vol_edge(1)) then
       dest1 = 1
       return
    end if
    nnew = n
    if (vbar_aft > vol_edge(n + 1)) then
       do while ((nnew < n_bin) .and. (vbar_aft > vol_edge(nnew + 1)))
          nnew = nnew + 1
       end do
    else if (vbar_aft < vol_edge(n)) then
       do while ((nnew > 1) .and. (vbar_aft < vol_edge(nnew)))
          nnew = nnew - 1
       end do
    end if
    dest1 = nnew

    ! need a valid pre-growth distribution to do better than moving-center
    if (vtot_pre <= 0d0) return

    vlo = vol_edge(n)
    vhi = vol_edge(n + 1)
    vdel = vhi - vlo
    vbar_pre = vtot_pre / num

    ! pre-growth mean too close to (or outside) the bin edges -> moving-center
    if ((vbar_pre >= vhi - 0.01d0 * vdel) .or. &
        (vbar_pre <= vlo + 0.01d0 * vdel)) return

    ! linear sub-bin reconstruction n(z) = aa + bb*z, z in [0,1] over [vlo,vhi],
    ! matching the bin number and mean volume (with edge-clamping that keeps the
    ! linear density non-negative, as in MOSAIC move_sections)
    gamma = vhi / vlo - 1d0
    ratio = vbar_pre / vlo
    if (ratio <= 1.0001d0 + gamma / 3d0) then
       vtmp = vlo + 3d0 * (vbar_pre - vlo)
       vhi = min(vhi, vtmp)
       vdel = vhi - vlo
       gamma = vhi / vlo - 1d0
       ratio = vbar_pre / vlo
    else if (ratio >= 0.9999d0 + gamma * 2d0 / 3d0) then
       vtmp = vhi + 3d0 * (vbar_pre - vhi)
       vlo = max(vlo, vtmp)
       vdel = vhi - vlo
       gamma = vhi / vlo - 1d0
       ratio = vbar_pre / vlo
    end if
    bb = (ratio - 1d0 - 0.5d0 * gamma) * 12d0 / gamma
    aa = 1d0 - 0.5d0 * bb

    ! destination bin edges mapped back into pre-growth volume space
    vcutlo = vol_edge(nnew)     * (vbar_pre / vbar_aft)
    vcuthi = vol_edge(nnew + 1) * (vbar_pre / vbar_aft)

    ! choose the neighbour bin, or fall back to moving-center if the grown bin
    ! sits entirely within the destination
    if (nnew == 1) then
       if (vhi <= vcuthi) return
       nnew2 = nnew + 1
    else if (nnew == n_bin) then
       if (vlo >= vcutlo) return
       nnew2 = nnew - 1
    else
       if ((vlo >= vcutlo) .and. (vhi <= vcuthi)) return
       if (vlo < vcutlo) then
          nnew2 = nnew - 1
       else
          nnew2 = nnew + 1
       end if
    end if

    ! integrate the linear distribution over the part landing in dest1
    zlo = max(0d0, (vcutlo - vlo) / vdel)
    zhi = min(1d0, (vcuthi - vlo) / vdel)
    d1 = zhi - zlo
    d2 = (zhi**2 - zlo**2) * 0.5d0
    d3 = (zhi**3 - zlo**3) / 3d0
    frac_num = aa * d1 + bb * d2
    frac_vol = (vlo / vbar_pre) &
         * (aa * d1 + (aa * gamma + bb) * d2 + (bb * gamma) * d3)

    if ((frac_num <= 0d0) .or. (frac_vol <= 0d0)) then
       ! all goes to the neighbour
       dest1 = nnew2
       dest2 = 0
       frac_num = 1d0
       frac_vol = 1d0
    else if ((frac_num >= 1d0) .or. (frac_vol >= 1d0)) then
       ! all stays in dest1
       dest2 = 0
       frac_num = 1d0
       frac_vol = 1d0
    else
       dest2 = nnew2
    end if

  end subroutine aero_binned_two_moment_split

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the number of bytes required to pack the structure.
  !!
  !! See pmc_mpi for usage details.
  integer function pmc_mpi_pack_size_aero_binned(val)

    !> Structure to pack.
    type(aero_binned_t), intent(in) :: val

    pmc_mpi_pack_size_aero_binned = &
         pmc_mpi_pack_size_real_array(val%num_conc) &
         + pmc_mpi_pack_size_real_array_2d(val%vol_conc)

  end function pmc_mpi_pack_size_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the structure into the buffer and advance position.
  !!
  !! See pmc_mpi for usage details.
  subroutine pmc_mpi_pack_aero_binned(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Structure to pack.
    type(aero_binned_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_real_array(buffer, position, val%num_conc)
    call pmc_mpi_pack_real_array_2d(buffer, position, val%vol_conc)
    call assert(348207873, &
         position - prev_position <= pmc_mpi_pack_size_aero_binned(val))
#endif

  end subroutine pmc_mpi_pack_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the structure from the buffer and advance position.
  !!
  !! See pmc_mpi for usage details.
  subroutine pmc_mpi_unpack_aero_binned(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Structure to unpack into (must not be allocated).
    type(aero_binned_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_real_array(buffer, position, val%num_conc)
    call pmc_mpi_unpack_real_array_2d(buffer, position, val%vol_conc)
    call assert(878267066, &
         position - prev_position <= pmc_mpi_pack_size_aero_binned(val))
#endif

  end subroutine pmc_mpi_unpack_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of the structure across all processes,
  !> storing the result on the root process.
  subroutine pmc_mpi_reduce_avg_aero_binned(val, val_avg)

    !> Per-process value to average.
    type(aero_binned_t), intent(in) :: val
    !> Averaged result (only valid on root process).
    type(aero_binned_t), intent(inout) :: val_avg

    call pmc_mpi_reduce_avg_real_array(val%num_conc, val_avg%num_conc)
    call pmc_mpi_reduce_avg_real_array_2d(val%vol_conc, val_avg%vol_conc)

  end subroutine pmc_mpi_reduce_avg_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine aero_binned_output_netcdf(aero_binned, ncid, bin_grid, &
       aero_data)

    !> Aero_binned to write.
    type(aero_binned_t), intent(in) :: aero_binned
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> aero_data structure.
    type(aero_data_t), intent(in) :: aero_data

    integer :: dimid_aero_diam, dimid_aero_species
    real(kind=dp) :: mass_conc(bin_grid_size(bin_grid), &
         aero_data_n_spec(aero_data))
    integer :: i_bin

    !> \page output_format_aero_binned Output File Format: Aerosol Binned Sectional State
    !!
    !! The aerosol size distributions (number and mass) are stored on
    !! a logarmithmic grid (see the \ref output_format_diam_bin_grid
    !! section). To compute the total number or mass concentration,
    !! compute the sum over \c i of <tt>aero_number_concentration(i) *
    !! aero_diam_widths(i)</tt>, for example.
    !!
    !! The aerosol binned sectional state uses the \c aero_species
    !! NetCDF dimension as specified in the \ref
    !! output_format_aero_data section, as well as the \c aero_diam
    !! NetCDF dimension specified in the \ref
    !! output_format_diam_bin_grid section.
    !!
    !! The aerosol binned sectional state NetCDF variables are:
    !!   - \b aero_number_concentration (unit 1/m^3, dim \c aero_diam): the
    !!     number size distribution for the aerosol population,
    !!     \f$ dN(r)/d\ln r \f$, per bin
    !!   - \b aero_mass_concentration (unit kg/m^3, dim
    !!     <tt>dimid_aero_diam x dimid_aero_species</tt>): the mass size
    !!     distribution for the aerosol population,
    !!     \f$ dM(r,s)/d\ln r \f$, per bin and per species

    ! output_format_diam_bin_grid is here, as this is the only place it's used

    !> \page output_format_diam_bin_grid Output File Format: Diameter Bin Grid Data
    !!
    !! The aerosol diameter bin grid data NetCDF dimensions are:
    !!   - \b aero_diam: number of bins (grid cells) on the diameter axis
    !!   - \b aero_diam_edges: number of bin edges (grid cell edges) on
    !!     the diameter axis --- always equal to <tt>aero_diam + 1</tt>
    !!
    !! The aerosol diameter bin grid data NetCDF variables are:
    !!   - \b aero_diam (unit m, dim \c aero_diam): aerosol diameter axis
    !!     bin centers --- centered on a logarithmic scale from the edges, so
    !!     that <tt>aero_diam(i) / aero_diam_edges(i) =
    !!     sqrt(aero_diam_edges(i+1) / aero_diam_edges(i))</tt>
    !!   - \b aero_diam_edges (unit m, dim \c aero_diam_edges): aersol
    !!     diameter axis bin edges (there is one more edge than center)
    !!   - \b aero_diam_widths (dimensionless, dim \c aero_diam):
    !!     the base-e logarithmic bin widths --- <tt>aero_diam_widths(i)
    !!     = ln(aero_diam_edges(i+1) / aero_diam_edges(i))</tt>, so
    !!     all bins have the same width
    !!
    !! See also:
    !!   - \ref input_format_diam_bin_grid --- the corresponding input format

    do i_bin = 1,bin_grid_size(bin_grid)
       mass_conc(i_bin,:) = aero_binned%vol_conc(i_bin,:) * aero_data%density
    end do

    call bin_grid_netcdf_dim(bin_grid, ncid, "aero_diam", "m", &
         dimid_aero_diam, "aerosol diameter", scale=2d0)
    call aero_data_netcdf_dim_aero_species(aero_data, ncid, dimid_aero_species)

    call pmc_nc_write_real_1d(ncid, aero_binned%num_conc, &
         "aero_number_concentration", (/ dimid_aero_diam /), &
         unit="1/m^3", &
         long_name="aerosol number size concentration distribution", &
         description="logarithmic number size concentration, " &
         // "d N(r)/d ln D --- multiply by aero_diam_widths(i) " &
         // "and sum over i to compute the total number concentration")
    call pmc_nc_write_real_2d(ncid, mass_conc, &
         "aero_mass_concentration", &
         (/ dimid_aero_diam, dimid_aero_species /), unit="kg/m^3", &
         long_name="aerosol mass size concentration distribution", &
         description="logarithmic mass size concentration per species, " &
         // "d M(r,s)/d ln D --- multiply by aero_diam_widths(i) " &
         // "and sum over i to compute the total mass concentration of " &
         // "species s")

  end subroutine aero_binned_output_netcdf

  ! output_format_diam_bin_grid is here, as this is the only place it's used

  ! this belongs in the subroutine above, but is outside because
  ! Doxygen 1.8.7 doesn't resolve references when multiple \page
  ! blocks are in one subroutine

  !> \page output_format_diam_bin_grid Output File Format: Diameter Bin Grid Data
  !!
  !! The aerosol diameter bin grid data NetCDF dimensions are:
  !!   - \b aero_diam: number of bins (grid cells) on the diameter axis
  !!   - \b aero_diam_edges: number of bin edges (grid cell edges) on
  !!     the diameter axis --- always equal to <tt>aero_diam + 1</tt>
  !!
  !! The aerosol diameter bin grid data NetCDF variables are:
  !!   - \b aero_diam (unit m, dim \c aero_diam): aerosol diameter axis
  !!     bin centers --- centered on a logarithmic scale from the edges, so
  !!     that <tt>aero_diam(i) / aero_diam_edges(i) =
  !!     sqrt(aero_diam_edges(i+1) / aero_diam_edges(i))</tt>
  !!   - \b aero_diam_edges (unit m, dim \c aero_diam_edges): aersol
  !!     diameter axis bin edges (there is one more edge than center)
  !!   - \b aero_diam_widths (dimensionless, dim \c aero_diam):
  !!     the base-e logarithmic bin widths --- <tt>aero_diam_widths(i)
  !!     = ln(aero_diam_edges(i+1) / aero_diam_edges(i))</tt>, so
  !!     all bins have the same width
  !!
  !! See also:
  !!   - \ref input_format_diam_bin_grid --- the corresponding input format

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine aero_binned_input_netcdf(aero_binned, ncid, bin_grid, &
       aero_data)

    !> Aero_binned to write.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> aero_data structure.
    type(aero_data_t), intent(in) :: aero_data

    integer :: i_bin

    call pmc_nc_read_real_1d(ncid, aero_binned%num_conc, &
         "aero_number_concentration")
    call pmc_nc_read_real_2d(ncid, aero_binned%vol_conc, &
         "aero_mass_concentration")
    ! convert mass concentation to volume concentration
    do i_bin = 1,bin_grid_size(bin_grid)
       aero_binned%vol_conc(i_bin,:) = aero_binned%vol_conc(i_bin,:) &
            / aero_data%density
    end do

  end subroutine aero_binned_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_binned
