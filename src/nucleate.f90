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

  !> Type code for unknown or invalid nucleation type.
  integer, parameter :: NUCLEATE_TYPE_INVALID   = 0
  !> Type code for no nucleation.
  integer, parameter :: NUCLEATE_TYPE_NONE      = 1
  !> Type code for H2SO4 to SO4 nucleation with quadratic rate.
  integer, parameter :: NUCLEATE_TYPE_SULF_ACID = 2
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do nucleation of the type given by the first argument.
  subroutine nucleate(nucleate_type, bin_grid, env_state, gas_data, &
       aero_data, aero_weight, aero_state, gas_state, del_t)

    !> Type of nucleation.
    integer, intent(in) :: nucleate_type
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Time to perform nucleation for.
    real(kind=dp), intent(in) :: del_t

    if (nucleate_type == NUCLEATE_TYPE_SULF_ACID) then
       call nucleate_sulf_acid(bin_grid, env_state, gas_data, aero_data, &
            aero_weight, aero_state, gas_state, del_t)
    elseif (nucleate_type == NUCLEATE_TYPE_NONE) then
       ! do nothing
    else
       call die_msg(983831728, &
            "unknown nucleation type: " // integer_to_string(nucleate_type))
    end if

  end subroutine nucleate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Nucleate sulfuric acid into aerosol particles, using a power-law
  !> dependence, for time \c del_t.
  !!
  !! The modeled emission rate is \f$ J = K H^2 \f$, where \f$H\f$ is
  !! the concentration of \f$ \rm H_2SO_4 \f$ and \f$K\f$ is a
  !! constant coefficient.
  !!
  !! The reference is:
  !!
  !! C. Kuang, P. H. McMurry, A. V. McCormick, and F. L. Eisele
  !! (2008), Dependence of nucleation rates on sulfuric acid vapor
  !! concentration in diverse atmospheric locations,
  !! <i>J. Geophys. Res.</i>, 113, D10209, doi:<a
  !! href="http://dx.doi.org/10.1029/2007JD009253">10.1029/2007JD009253</a>.
  subroutine nucleate_sulf_acid(bin_grid, env_state, gas_data, aero_data, &
       aero_weight, aero_state, gas_state, del_t)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
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
    real(kind=dp) :: nucleate_comp_vol
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

    ! computational volume at the size of nucleated particles (only
    ! valid for mono-disperse nucleation)
    nucleate_comp_vol = aero_state%comp_vol &
         / aero_weight_value(aero_weight, nucleate_diam / 2d0)

    ! determine number of nucleated particles
    n_samp_avg = nucleate_rate * nucleate_comp_vol * del_t
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
         / nucleate_comp_vol                    ! molecules / m^3
    gas_state%mix_rat(i_gas_h2so4) = gas_state%mix_rat(i_gas_h2so4) &
         - env_state_conc_to_ppb(env_state, h2so4_removed_conc)
    if (gas_state%mix_rat(i_gas_h2so4) < 0d0) then
       gas_state%mix_rat(i_gas_h2so4) = 0d0
    end if

  end subroutine nucleate_sulf_acid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_file_read_nucleate(file, nucleate_type)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aerosol weight.
    integer, intent(out) :: nucleate_type

    character(len=SPEC_LINE_MAX_VAR_LEN) :: nucleate_type_name

    !> \page input_format_nucleate Input File Format: Nucleation Parameterization
    !!
    !! The nucleation parameterization is specified by the parameters:
    !!   - \b nucleate (string): the type of nucleation
    !!     parameterization --- must be one of: "none" for no
    !!     nucleation; or "sulf_acid" for the nucleate_sulf_acid()
    !!     parameterization
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    call spec_file_read_string(file, 'nucleate', nucleate_type_name)
    if (nucleate_type_name == 'none') then
       nucleate_type = NUCLEATE_TYPE_NONE
    elseif (nucleate_type_name == 'sulf_acid') then
       nucleate_type = NUCLEATE_TYPE_SULF_ACID
    else
       call spec_file_die_msg(707263678, file, "unknown nucleate type: " &
            // trim(nucleate_type_name))
    end if

  end subroutine spec_file_read_nucleate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_nucleate
