! Copyright (C) 2010 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_nucleate module.

!> Aerosol nucleation functions.
module pmc_nucleate

  use pmc_env_state
  use pmc_bin_grid
  use pmc_aero_state
  use pmc_aero_data
  use pmc_gas_data
  use pmc_gas_state
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do nucleation of the type given by the first argument.
  subroutine nucleate(nucleate_type, bin_grid, env_state, gas_data, &
       aero_data, aero_state, gas_state, del_t)

    !> Type of nucleation.
    character(len=*), intent(in) :: nucleate_type
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Time to perform nucleation for.
    real(kind=dp), intent(in) :: del_t

    if (nucleate_type == "sulf_acid") then
       call nucleate_sulf_acid(bin_grid, env_state, gas_data, aero_data, &
            aero_state, gas_state, del_t)
    elseif (nucleate_type == "none") then
       ! do nothing
    else
       call die_msg(983831728, &
            "unknown nucleation type: " // trim(nucleate_type))
    end if

  end subroutine nucleate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Nucleate sulfuric acid into aerosol particles, using a power-law
  !> dependence, for time \c del_t.
  !!
  !! The modeled emission rate is \f$ J = K H^2 \f$.
  !!
  !! The reference is:
  !!
  !! C. Kuang, P. H. McMurry, A. V. McCormick, and F. L. Eisele
  !! (2008), Dependence of nucleation rates on sulfuric acid vapor
  !! concentration in diverse atmospheric locations, J. Geophys. Res.,
  !! 113, D10209, doi:10.1029/2007JD009253.
  subroutine nucleate_sulf_acid(bin_grid, env_state, gas_data, aero_data, &
       aero_state, gas_state, del_t)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Time to perform nucleation for.
    real(kind=dp), intent(in) :: del_t

    real(kind=dp), parameter :: nucleate_coeff = 1d-18 ! K (m^3 s^{-1})
    real(kind=dp), parameter :: nucleate_diam = 1d-9   ! diameter of new particles (m)

    integer :: i_gas_h2so4, i_aero_so4, n_samp, i_samp, i_bin
    real(kind=dp) :: sulf_acid_conc, nucleate_rate, n_samp_avg
    real(kind=dp) :: total_so4_vol, vol, h2so4_removed_conc
    type(aero_particle_t) :: aero_particle

    ! look up the species numbers
    i_gas_h2so4 = gas_data_spec_by_name(gas_data, "H2SO4")
    call assert_msg(886839228, i_gas_h2so4 > 0, &
         "nucleate_sulf_acid requires H2SO4 as a gas species")
    i_aero_so4 = aero_data_spec_by_name(aero_data, "SO4")
    call assert_msg(551966998, i_aero_so4 > 0, &
         "nucleate_sulf_acid requires SO4 as an aerosol species")

    ! H2SO4 concentration in molecules / m^3
    sulf_acid_conc = env_state_ppb_to_conc(env_state, &
         gas_state%mix_rat(i_gas_h2so4))

    ! particle nucleation rate in (particles m^{-3} s^{-1})
    nucleate_rate = nucleate_coeff * sulf_acid_conc**2

    ! determine number of nucleated particles
    n_samp_avg = nucleate_rate * aero_state%comp_vol * del_t
    n_samp = rand_poisson(n_samp_avg)

    ! create the particles
    call aero_particle_allocate_size(aero_particle, aero_data%n_spec)
    total_so4_vol = 0d0
    do i_samp = 1,n_samp
       vol = diam2vol(nucleate_diam)
       total_so4_vol = total_so4_vol + vol

       aero_particle%vol(i_aero_so4) = vol
       call aero_particle_new_id(aero_particle)
       call aero_particle_set_create_time(aero_particle, &
            env_state%elapsed_time)
       i_bin = aero_particle_in_bin(aero_particle, bin_grid)
       call aero_state_add_particle(aero_state, i_bin, aero_particle)
    end do
    call aero_particle_deallocate(aero_particle)

    ! remove gases that formed new particles
    h2so4_removed_conc = &
         total_so4_vol * aero_data%density(i_aero_so4) &
         / aero_data%molec_weight(i_aero_so4) & ! moles of SO4
         * const%avagadro &                     ! molecules of SO4
         / aero_state%comp_vol                  ! molecules / m^3
    gas_state%mix_rat(i_gas_h2so4) = gas_state%mix_rat(i_gas_h2so4) &
         - env_state_conc_to_ppb(env_state, h2so4_removed_conc)
    if (gas_state%mix_rat(i_gas_h2so4) < 0d0) then
       gas_state%mix_rat(i_gas_h2so4) = 0d0
    end if

  end subroutine nucleate_sulf_acid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_nucleate
