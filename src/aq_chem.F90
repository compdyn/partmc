! Copyright (C) 2015 Matthew Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aq_chem module.

!> Interface to the aqueous chemistry code
module pmc_aq_chem

  use pmc_env_state
  use pmc_aq_mech_data
  use pmc_aq_spec_data
  use pmc_aq_state
  use pmc_gas_data
  use pmc_gas_state
  use pmc_aero_state
  use pmc_aq_integrate_f

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run aqueous chemistry mechanism over timestep del_t
  subroutine aq_chem_timestep(env_state_initial, env_state_final, aq_mech_data, &
        aq_spec_data, aq_state, aq_state_init, gas_data, gas_state, aero_data, &
        aero_state, del_t)

    !> Environment state at t0
    type(env_state_t), intent(in) :: env_state_initial
    !> Environment state at t0 + del_t
    type(env_state_t), intent(in) :: env_state_final
    !> Aqueous chemistry mechanism
    type(aq_mech_data_t), intent(in), target :: aq_mech_data
    !> Aqueous chemistry related species data
    type(aq_spec_data_t), intent(in), target :: aq_spec_data
    !> Aqueous chemistry state
    type(aq_state_t), intent(inout) :: aq_state
    !> Aqueous chemistry state containing initial and constant species concentrations
    type(aq_state_t), intent(inout) :: aq_state_init
    !> Gas data
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state
    type(gas_state_t), intent(inout) :: gas_state
    !> Aerosol data
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state
    type(aero_state_t), intent(inout) :: aero_state
    !> Timestep (s)
    real(kind=dp), intent(in) :: del_t

    ! current particle index
    integer :: part_index
    ! current sub time step
    integer :: i_time

    integer :: gly_idx, form_idx

    ! Water must be a species to do aqueous chemistry
    if (aero_data%i_water.eq.0) then
        call die_msg(585776938, 'Water must be an aerosol-phase ' &
            // 'species to do aqueous chemistry.')
    endif

    ! Reset the aqueous chemistry species concentrations to start
    ! with initial and constant concentrations
    call aq_state_copy(aq_state_init, aq_state)

    ! TEMP - FIX LATER
    ! Estimate gas-phase glyoxal as 3% of formaldehyde
    form_idx = gas_data_spec_by_name(gas_data, "HCHO")
    gly_idx = gas_data_spec_by_name(gas_data, "GLY")
    if (form_idx.eq.0) then
        call die_msg(586230727, 'Cannot find gas-phase index for formaldehyde.')
    elseif (gly_idx.eq.0) then
        call die_msg(586415604, 'Cannot find gas-phase index for glyoxal.')
    else
        gas_state%mix_rat(gly_idx) = gas_state%mix_rat(form_idx) * 0.03
    endif

    ! Map PMC gas-phase species to aqueous chemistry species
    call aq_chem_gas_spec_from_PMC(aq_spec_data, aq_state, gas_state)

    ! Cycle through each particle and do the aqueous chemistry
    do part_index=1,aero_state%apa%n_part

        ! Skip particles that have too little condensed water
        ! ( 5d-19 m^3 corresponds to ~ 1um diameter particle of water)
        if (aero_state%apa%particle(part_index)%vol(aero_data%i_water).le.5d-19) then
            cycle
        endif

        ! Map PMC aerosol-phase species to aqueous chemistry species
        call aq_chem_aero_spec_from_PMC(aq_spec_data, aq_state, aero_data, &
            aero_state%apa%particle(part_index))

        ! Set the number and size of the particles in this aero_particle_t variable
        aq_state%n_particle = aero_weight_array_num_conc(aero_state%awa, &
                                aero_state%apa%particle(part_index)) / 1d6 ! (#/cc)
        aq_state%radius = aero_particle_radius(aero_state%apa%particle(part_index)) ! (m)

        ! Solve aqueous chemistry mechanism
        call aq_integrate(aq_state, aq_mech_data, aq_spec_data, env_state_initial, &
            env_state_final, del_t)

        ! Map aqueous chemistry species back to PMC aerosol-phase species
        call aq_chem_aero_spec_to_PMC(aq_spec_data, aq_state, aero_data, &
            aero_state%apa%particle(part_index))

    enddo

    ! Map aqueous chemistry species back to PMC gas-phase species
    call aq_chem_gas_spec_to_PMC(aq_spec_data, aq_state, gas_state)

  end subroutine aq_chem_timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Map PMC gas-phase species to aqueous chemistry species
  subroutine aq_chem_gas_spec_from_PMC(aq_spec_data, aq_state, gas_state)

    !> Aqueous chemistry related species data
    type(aq_spec_data_t), intent(in), target :: aq_spec_data
    !> Aqueous chemistry state
    type(aq_state_t), intent(inout) :: aq_state
    !> Gas state
    type(gas_state_t), intent(inout) :: gas_state

    ! Convert gas-phase species from ppb (PMC) to atm (aq. chemistry)
    where (aq_spec_data%pmc_gas_index.ne.0)
        aq_state%mix_rat = gas_state%mix_rat(aq_spec_data%pmc_gas_index) * 1d-9
    endwhere

  end subroutine aq_chem_gas_spec_from_PMC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Map aqueous chemistry species back to PMC gas-phase species
  subroutine aq_chem_gas_spec_to_PMC(aq_spec_data, aq_state, gas_state)

    !> Aqueous chemistry related species data
    type(aq_spec_data_t), intent(in), target :: aq_spec_data
    !> Aqueous chemistry state
    type(aq_state_t), intent(inout) :: aq_state
    !> Gas state
    type(gas_state_t), intent(inout) :: gas_state

    ! Convert gas-phase species from atm (aq. chemistry) to ppb (PMC)
    where (aq_spec_data%pmc_gas_index.ne.0)
         gas_state%mix_rat(aq_spec_data%pmc_gas_index) = aq_state%mix_rat / 1d-9
    endwhere

  end subroutine aq_chem_gas_spec_to_PMC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Map PMC aerosol-phase species to aqueous chemistry species
  subroutine aq_chem_aero_spec_from_PMC(aq_spec_data, aq_state, aero_data, &
                    particle)

    !> Aqueous chemistry related species data
    type(aq_spec_data_t), intent(in), target :: aq_spec_data
    !> Aqueous chemistry state
    type(aq_state_t), intent(inout) :: aq_state
    !> Aerosol data
    type(aero_data_t), intent(in) :: aero_data
    !> Particle state
    type(aero_particle_t), intent(inout) :: particle

    ! aerosol-phase water (L)
    real(kind=dp) :: water_conc

    ! Convert water concentration from m^3 (PMC) to L (aq. chemistry)
    water_conc = particle%vol(aero_data%i_water) * real(1d3, kind=dp)

    ! Convert aerosol-phase species from m^3 (PMC) to mol/L (aq. chemistry)
    where (aq_spec_data%pmc_aero_index.ne.0)
        aq_state%mix_rat = particle%vol(aq_spec_data%pmc_aero_index) &
            * aero_data%density(aq_spec_data%pmc_aero_index) &
            / aero_data%molec_weight(aq_spec_data%pmc_aero_index) &
            / water_conc
    endwhere

  end subroutine aq_chem_aero_spec_from_PMC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Map aqueous chemistry species back to PMC aerosol-phase species
  subroutine aq_chem_aero_spec_to_PMC(aq_spec_data, aq_state, aero_data, &
                    particle)

    !> Aqueous chemistry related species data
    type(aq_spec_data_t), intent(in), target :: aq_spec_data
    !> Aqueous chemistry state
    type(aq_state_t), intent(inout) :: aq_state
    !> Aerosol data
    type(aero_data_t), intent(in) :: aero_data
    !> Particle state
    type(aero_particle_t), intent(inout) :: particle

    ! aerosol-phase water (L)
    real(kind=dp) :: water_conc

    ! Water does not condense or evaporate in the aqueous-phase mechanism
    ! so use the original PMC water concentration to calculate new
    ! aerosol-phase species concentrations

    ! Convert water concentration from m^3 (PMC) to L (aq. chemistry)
    water_conc = particle%vol(aero_data%i_water) * real(1d3, kind=dp)

    ! Convert aerosol-phase species from mol/L (aq. chemistry) to m^3 (PMC)
    where (aq_spec_data%pmc_aero_index.ne.0)
        particle%vol(aq_spec_data%pmc_aero_index) = &
            aq_state%mix_rat / aero_data%density(aq_spec_data%pmc_aero_index) &
            * aero_data%molec_weight(aq_spec_data%pmc_aero_index) * water_conc
    endwhere

  end subroutine aq_chem_aero_spec_to_PMC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aq_chem









