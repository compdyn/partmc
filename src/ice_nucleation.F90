! Copyright (C) 2024 University of Illinois at Urbana-Champaign
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
! The pmc_ice_nucleation module.

module pmc_ice_nucleation
  use pmc_aero_state
  use pmc_env_state
  use pmc_aero_data
  use pmc_util
  use pmc_aero_particle
  use pmc_constants
  use pmc_rand

  implicit none

  !> Type code for an undefined or invalid immersion freezing scheme.
  integer, parameter :: IMMERSION_FREEZING_SCHEME_INVALID = 0
  !> Type code for constant ice nucleation rate (J_het) immersion freezing scheme.
  integer, parameter :: IMMERSION_FREEZING_SCHEME_CONST = 1
  !> Type code for the singular (INAS) immersion freezing scheme.
  integer, parameter :: IMMERSION_FREEZING_SCHEME_SINGULAR = 2
  !> Type code for the ABIFM immersion freezing scheme.
  integer, parameter :: IMMERSION_FREEZING_SCHEME_ABIFM = 3

  !> Whether to use the binned-tau leaping algorithm for time-dependent scheme.
  logical :: do_speedup = .true.

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Main subroutine for immersion freezing simulation.
  subroutine ice_nucleation_immersion_freezing(aero_state, aero_data, &
       env_state, del_t, immersion_freezing_scheme_type, freezing_rate)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Total time to integrate.
    real(kind=dp), intent(in) :: del_t
    !> Immersion freezing scheme type.
    integer, intent(in) :: immersion_freezing_scheme_type
    !> Freezing rate (only used for the constant rate scheme).
    real(kind=dp), intent(in) :: freezing_rate

    !> Call the immersion freezing subroutine according to the immersion
    !> freezing scheme.
    if (env_state%temp <= const%water_freeze_temp) then
       if ((immersion_freezing_scheme_type == IMMERSION_FREEZING_SCHEME_ABIFM) &
            .OR. (immersion_freezing_scheme_type == IMMERSION_FREEZING_SCHEME_CONST)) then
          if (do_speedup) then
             call ice_nucleation_immersion_freezing_time_dependent( &
                  aero_state, aero_data, env_state, del_t, &
                  immersion_freezing_scheme_type, freezing_rate)
          else
             call ice_nucleation_immersion_freezing_time_dependent_naive( &
                  aero_state, aero_data, env_state, del_t, &
                  immersion_freezing_scheme_type, freezing_rate)
          end if
       else if (immersion_freezing_scheme_type == &
            IMMERSION_FREEZING_SCHEME_SINGULAR) then
          call ice_nucleation_immersion_freezing_singular(aero_state, &
               aero_data, env_state)
       else
          call assert_msg(121370299, .false., &
               'Error type of immersion freezing scheme')
       endif
    endif

  end subroutine ice_nucleation_immersion_freezing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialization for the sigular scheme, sampling the freezing temperature
  !> for each particles.
  subroutine ice_nucleation_singular_initialize(aero_state, aero_data)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    integer :: i_part
    real(kind=dp) :: a_INAS, b_INAS, p, S, T0, temp
    real(kind=dp) :: aerosol_diameter

    T0 = const%water_freeze_temp
    a_INAS = -0.517
    b_INAS = 8.934
    do i_part = 1, aero_state_n_part(aero_state)
       aerosol_diameter = aero_particle_dry_diameter( &
            aero_state%apa%particle(i_part), aero_data)
       S = const%pi * aerosol_diameter **2
       p = pmc_random()
       temp = (log(1 - p) + exp(-S * exp(-a_INAS * T0 + b_INAS))) / (-S)
       aero_state%apa%particle(i_part)%imf_temperature = T0 + (log(temp) &
          - b_INAS) / a_INAS
    end do
  end subroutine ice_nucleation_singular_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Simulation for singular scheme, deciding whether to freeze for each
  !> particle. Run in each time step.
  subroutine ice_nucleation_immersion_freezing_singular(aero_state, aero_data,&
       env_state)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state

    real(kind=dp), allocatable :: H2O_masses(:), total_masses(:), &
        H2O_frac(:)
    integer :: i_part

    ! FIXME: Do this to avoid compiler warning/error, fix it in the future.
    allocate(total_masses(aero_state_n_part(aero_state)))
    allocate(H2O_masses(aero_state_n_part(aero_state)))
    allocate(H2O_frac(aero_state_n_part(aero_state)))

    total_masses = aero_state_masses(aero_state, aero_data)
    H2O_masses = aero_state_masses(aero_state, aero_data, include=["H2O"])
    H2O_frac = H2O_masses / total_masses
    do i_part = 1, aero_state_n_part(aero_state)
       if (aero_state%apa%particle(i_part)%frozen) then
          cycle
       end if
       if (H2O_frac(i_part) < const%imf_water_threshold) then
          cycle
       end if
       if (env_state%temp <= &
            aero_state%apa%particle(i_part)%imf_temperature) then
          aero_state%apa%particle(i_part)%frozen = .TRUE.
          aero_state%apa%particle(i_part)%den_ice = &
               const%reference_ice_density
          aero_state%apa%particle(i_part)%ice_shape_phi = 1d0
       end if
    end do

    deallocate(total_masses)
    deallocate(H2O_masses)
    deallocate(H2O_frac)

  end subroutine ice_nucleation_immersion_freezing_singular

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Simulation for time-dependent scheme (e.g., ABIFM, constant rate),
  !> deciding whether to freeze for each particle. Run in each time step.
  !> This subroutine applys the binned-tau leaping algorithm for speeding up.
  subroutine ice_nucleation_immersion_freezing_time_dependent(aero_state, &
       aero_data, env_state, del_t, immersion_freezing_scheme_type, &
       freezing_rate)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Total time to integrate.
    real(kind=dp), intent(in) :: del_t
    !> Freezing rate (only used for the constant rate scheme).
    real(kind=dp), intent(in) :: freezing_rate

    integer :: i_part, i_bin, i_class, n_bins, n_class
    integer, intent(in) :: immersion_freezing_scheme_type
    real(kind=dp) :: a_w_ice, pis, pvs
    real(kind=dp) :: p_freeze

    real(kind=dp) :: p_freeze_max, radius_max, diameter_max

    integer :: k_th, n_parts_in_bin
    real(kind=dp) :: rand
    real(kind=dp), allocatable :: H2O_masses(:), total_masses(:), &
         H2O_frac(:)
    integer :: i_spec_max
    real(kind=dp) :: j_het_max
    integer :: rand_geo


    allocate(total_masses(aero_state_n_part(aero_state)))
    allocate(H2O_masses(aero_state_n_part(aero_state)))
    allocate(H2O_frac(aero_state_n_part(aero_state)))

    call aero_state_sort(aero_state, aero_data)

    total_masses = aero_state_masses(aero_state, aero_data)
    H2O_masses = aero_state_masses(aero_state, aero_data, include=["H2O"])
    H2O_frac = H2O_masses / total_masses
    pvs = env_state_saturated_vapor_pressure_water(env_state%temp)
    pis = env_state_saturated_vapor_pressure_ice(env_state%temp)
    a_w_ice = pis / pvs
    if (immersion_freezing_scheme_type == IMMERSION_FREEZING_SCHEME_ABIFM) then
       call ABIFM_max_spec(aero_data, a_w_ice, i_spec_max, j_het_max)
    endif

    n_bins = aero_sorted_n_bin(aero_state%aero_sorted)
    n_class = aero_sorted_n_class(aero_state%aero_sorted)

    loop_bins: do i_bin = 1, n_bins
       loop_classes: do i_class = 1, n_class
          n_parts_in_bin = integer_varray_n_entry(&
               aero_state%aero_sorted%size_class%inverse(i_bin, i_class))
          radius_max = aero_state%aero_sorted%bin_grid%edges(i_bin + 1)
          diameter_max = radius_max * 2
          if (immersion_freezing_scheme_type == &
               IMMERSION_FREEZING_SCHEME_ABIFM) then
             p_freeze_max = ABIFM_Pfrz_max(diameter_max, aero_data, j_het_max, &
                   del_t)
          else if (immersion_freezing_scheme_type == &
               IMMERSION_FREEZING_SCHEME_CONST) then
             p_freeze_max = 1 - exp(freezing_rate * del_t)
          endif

          k_th = n_parts_in_bin + 1
          loop_choosed_particles: do while(.TRUE.)
             rand_geo = rand_geometric(p_freeze_max)
             k_th = k_th - rand_geo
             if (k_th <= 0) then
                EXIT loop_choosed_particles
             endif
             i_part = aero_state%aero_sorted%size_class &
                  %inverse(i_bin, i_class)%entry(k_th)
             if (aero_state%apa%particle(i_part)%frozen) then
                cycle
             end if
             if (H2O_frac(i_part) < const%imf_water_threshold) then
                cycle
             end if
             if (immersion_freezing_scheme_type == &
                  IMMERSION_FREEZING_SCHEME_ABIFM) then
                p_freeze = ABIFM_Pfrz_particle( & 
                     aero_state%apa%particle(i_part), aero_data, a_w_ice, del_t)
                call warn_assert_msg(301184565, p_freeze <= p_freeze_max,&
                     "p_freeze > p_freeze_max.")
                rand = pmc_random()
                if (rand < p_freeze / p_freeze_max) then
                   aero_state%apa%particle(i_part)%frozen = .TRUE.
                   aero_state%apa%particle(i_part)%den_ice = &
                        const%reference_ice_density
                   aero_state%apa%particle(i_part)%ice_shape_phi = 1d0

                endif
             else
                aero_state%apa%particle(i_part)%frozen = .TRUE.
                aero_state%apa%particle(i_part)%den_ice = &
                     const%reference_ice_density
                aero_state%apa%particle(i_part)%ice_shape_phi = 1d0
             endif

          enddo loop_choosed_particles
       enddo loop_classes
    enddo loop_bins

    deallocate(total_masses)
    deallocate(H2O_masses)
    deallocate(H2O_frac)

  end subroutine ice_nucleation_immersion_freezing_time_dependent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Simulation for time-dependent scheme (e.g., ABIFM, constant rate),
  !> deciding whether to freeze for each particle. Run in each time step.
  !> This subroutine applies the naive algorithm that checks each particle.
  subroutine ice_nucleation_immersion_freezing_time_dependent_naive( &
       aero_state, aero_data, env_state, del_t, &
       immersion_freezing_scheme_type, freezing_rate)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Total time to integrate.
    real(kind=dp), intent(in) :: del_t
    !> Type of the immersion freezing scheme
    integer, intent(in) :: immersion_freezing_scheme_type
    !> Freezing rate (only used for the constant rate scheme).
    real(kind=dp), intent(in) :: freezing_rate

    integer :: i_part
    real(kind=dp) :: a_w_ice, pis, pvs
    real(kind=dp) :: p_freeze = 0
    real(kind=dp), allocatable :: H2O_masses(:), total_masses(:), &
         H2O_frac(:)
    real(kind=dp) :: rand


    ! FIXME: Do this to avoid compiler warning/error, fix it in the future.
    allocate(total_masses(aero_state_n_part(aero_state)))
    allocate(H2O_masses(aero_state_n_part(aero_state)))
    allocate(H2O_frac(aero_state_n_part(aero_state)))

    total_masses = aero_state_masses(aero_state, aero_data)
    H2O_masses = aero_state_masses(aero_state, aero_data, include=["H2O"])
    H2O_frac = H2O_masses / total_masses

    pvs = env_state_saturated_vapor_pressure_water(env_state%temp)
    pis = env_state_saturated_vapor_pressure_ice(env_state%temp)
    a_w_ice = pis / pvs

    do i_part = 1, aero_state_n_part(aero_state)
       if (aero_state%apa%particle(i_part)%frozen) cycle
       if (H2O_frac(i_part) < const%imf_water_threshold) cycle
       rand = pmc_random()

       if (immersion_freezing_scheme_type == &
            IMMERSION_FREEZING_SCHEME_ABIFM) then
          p_freeze = ABIFM_Pfrz_particle(aero_state%apa%particle(i_part), &
               aero_data, a_w_ice, del_t)
       else if (immersion_freezing_scheme_type == &
            IMMERSION_FREEZING_SCHEME_CONST) then
          p_freeze = 1 - exp(freezing_rate * del_t)
       end if

       if (rand < p_freeze) then
          aero_state%apa%particle(i_part)%frozen = .TRUE.
          aero_state%apa%particle(i_part)%den_ice = &
               const%reference_ice_density
          aero_state%apa%particle(i_part)%ice_shape_phi = 1d0
       end if
    end do

    deallocate(total_masses)
    deallocate(H2O_masses)
    deallocate(H2O_frac)

  end subroutine ice_nucleation_immersion_freezing_time_dependent_naive

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> This subroutine simulates the melting process.(Set frozen = .False.
  !> for each particle once the temperature is higher than water freezing
  !> tempearture)
  subroutine ice_nucleation_melting(aero_state, aero_data, env_state)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state

    integer :: i_part

    if (env_state%temp > const%water_freeze_temp) then
       do i_part = 1, aero_state_n_part(aero_state)
         aero_state%apa%particle(i_part)%frozen = .false.
         aero_state%apa%particle(i_part)%den_ice = const%nan
         aero_state%apa%particle(i_part)%ice_shape_phi = const%nan
       end do
    end if

  end subroutine ice_nucleation_melting

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculating the freezing probability for the particle (i_part) using ABIFM
  !> method (Knopf et al.,2013)
  real(kind=dp) function ABIFM_Pfrz_particle(aero_particle, aero_data, &
       a_w_ice, del_t)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> The water activity w.r.t. ice.
    real(kind=dp), intent(in) :: a_w_ice
    !> Time interval.
    real(kind=dp), intent(in) :: del_t

    real(kind=dp) :: aerosol_diameter
    real(kind=dp) :: immersed_surface_area
    real(kind=dp) :: total_vol
    real(kind=dp) :: surface_ratio
    real(kind=dp) :: abifm_m, abifm_c
    real(kind=dp) :: j_het, j_het_x_area
    integer :: i_spec

    aerosol_diameter =  aero_particle_dry_diameter(aero_particle, aero_data)
    immersed_surface_area = const%pi * aerosol_diameter **2

    total_vol = 0d0
    total_vol = aero_particle_dry_volume(aero_particle, aero_data)

    j_het_x_area = 0d0
    do i_spec = 1,aero_data_n_spec(aero_data)
       if (i_spec == aero_data%i_water) cycle
       abifm_m = aero_data%abifm_m(i_spec)
       abifm_c = aero_data%abifm_c(i_spec)
       j_het = 10 ** (abifm_m  * (1 - a_w_ice) + abifm_c) * 10000
       surface_ratio = aero_particle%vol(i_spec) / total_vol
       j_het_x_area = j_het_x_area + j_het * immersed_surface_area * &
            surface_ratio
    end do

    ABIFM_Pfrz_particle = 1 - exp(-j_het_x_area * del_t)

  end function ABIFM_Pfrz_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculating the heterogeneous ice nucleation rate coefficient
  !> for each species, and fining the species having the largest rate.
  subroutine ABIFM_max_spec(aero_data, a_w_ice, i_spec_max, j_het_max)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> The water activity w.r.t. ice.
    real(kind=dp), intent(in) :: a_w_ice
    !> The index of the maximum J_het species.
    integer, intent(out) :: i_spec_max
    !> The maximum value of J_het among all species.
    real(kind=dp), intent(out) :: j_het_max

    real(kind=dp) :: abifm_m, abifm_c
    real(kind=dp) :: j_het
    integer :: i_spec

    j_het_max = const%nan
    do i_spec = 1,aero_data_n_spec(aero_data)
       if (i_spec == aero_data%i_water) cycle
       abifm_m = aero_data%abifm_m(i_spec)
       abifm_c = aero_data%abifm_c(i_spec)
       j_het = 10 ** (abifm_m  * (1 - a_w_ice) + abifm_c) * 10000
       if (j_het > j_het_max) then
          j_het_max = j_het
          i_spec_max = i_spec
       end if
    end do

  end subroutine ABIFM_max_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculating the maximum freezing probability for particles in
  !> one bin using ABIFM method (Knopf et al.,2013). Only used by
  !> the binned-tau leaping algorithm.
  real(kind=dp) function ABIFM_Pfrz_max(diameter_max, aero_data, j_het_max, &
       del_t)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Maximum diameter.
    real(kind=dp), intent(in) :: diameter_max
    !> Time interval.
    real(kind=dp), intent(in) :: del_t
    !> Maximum J_het among all species.
    real(kind=dp), intent(in) :: j_het_max

    real(kind=dp) :: immersed_surface_area

    immersed_surface_area = const%pi * diameter_max **2
    ABIFM_Pfrz_max = 1 - exp(-j_het_max * immersed_surface_area  * del_t)

  end function ABIFM_Pfrz_max

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the specification for a kernel type from a spec file and
  !> generate it.
  subroutine spec_file_read_immersion_freezing_scheme_type(file, &
       immersion_freezing_scheme_type)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Kernel type.
    integer, intent(out) :: immersion_freezing_scheme_type

    character(len=SPEC_LINE_MAX_VAR_LEN) :: imf_scheme

    !> \page input_format_immersion_freezing_scheme Input File Format: Immersion freezing scheme
    !!
    !! The immersion freezing scheme is specified by the parameter:
    !!   - \b immersion_freezing_scheme (string): the type of immersion freezing scheme
    !!   must be one of: \c sedi for the gravitational sedimentation
    !!   kernel; \c additive for the additive kernel; \c constant
    !!   for the constant kernel; \c brown for the Brownian kernel,
    !!   or \c zero for no immersion freezing
    !!
    !! If \c immersion_freezing_scheme is \c additive, the kernel coefficient needs to be
    !! provided using the \c additive_kernel_coeff parameter
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    call spec_file_read_string(file, 'immersion_freezing_scheme', imf_scheme)
    if (trim(imf_scheme) == 'const') then
       immersion_freezing_scheme_type = IMMERSION_FREEZING_SCHEME_CONST
    elseif (trim(imf_scheme) == 'singular') then
       immersion_freezing_scheme_type = IMMERSION_FREEZING_SCHEME_SINGULAR
    elseif (trim(imf_scheme) == 'ABIFM') then
       immersion_freezing_scheme_type = IMMERSION_FREEZING_SCHEME_ABIFM
    else
       call spec_file_die_msg(920761229, file, &
        "Unknown immersion freezing scheme: " // trim(imf_scheme))
    end if

  end subroutine spec_file_read_immersion_freezing_scheme_type

end module pmc_ice_nucleation
