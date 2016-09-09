! Copyright (C) 2010-2015 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_nucleate module.

!> Aerosol nucleation functions.
module pmc_nucleate

  use pmc_env_state
  use pmc_aero_state
  use pmc_aero_data
  use pmc_gas_data
  use pmc_gas_state

  !> Type code for unknown or invalid nucleation type.
  integer, parameter :: NUCLEATE_TYPE_INVALID   = 0
  !> Type code for H2SO4 to SO4 nucleation with quadratic rate.
  integer, parameter :: NUCLEATE_TYPE_SULF_ACID = 1

  !> Source name for nucleated particles.
  character(len=AERO_SOURCE_NAME_LEN), parameter :: NUCLEATE_SOURCE_NAME &
       = "nucleate"

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do nucleation of the type given by the first argument.
  subroutine nucleate(nucleate_type, nucleate_source, env_state, gas_data, &
       aero_data, aero_state, gas_state, del_t, allow_doubling, allow_halving)

    !> Type of nucleation.
    integer, intent(in) :: nucleate_type
    !> Nucleate source number.
    integer, intent(in) :: nucleate_source
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
    !> Whether to allow doubling of the population.
    logical, intent(in) :: allow_doubling
    !> Whether to allow halving of the population.
    logical, intent(in) :: allow_halving

    if (nucleate_type == NUCLEATE_TYPE_SULF_ACID) then
       call nucleate_sulf_acid(nucleate_source, env_state, gas_data, &
            aero_data, aero_state, gas_state, del_t, allow_doubling, &
            allow_halving)
    else
       call die_msg(983831728, &
            "unknown nucleation type: " &
            // trim(integer_to_string(nucleate_type)))
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
  subroutine nucleate_sulf_acid(nucleate_source, env_state, gas_data, &
       aero_data, aero_state, gas_state, del_t, allow_doubling, allow_halving)

    !> Nucleate source number.
    integer, intent(in) :: nucleate_source
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
    !> Whether to allow doubling of the population.
    logical, intent(in) :: allow_doubling
    !> Whether to allow halving of the population.
    logical, intent(in) :: allow_halving

    real(kind=dp), parameter :: nucleate_coeff = 1d-18 ! K (m^3 s^{-1})
    real(kind=dp), parameter :: nucleate_diam = 1d-9   ! diameter of new
                                                       ! particles (m)

    integer :: i_gas_h2so4, i_aero_so4, n_samp, i_samp, i_bin, i_group, n_group
    integer :: i_class
    real(kind=dp) :: sulf_acid_conc, nucleate_rate, n_samp_avg
    real(kind=dp) :: total_so4_vol, so4_vol, h2so4_removed_conc
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

    ! weight class to nucleate into
    i_class = aero_state_weight_class_for_source(aero_state, nucleate_source)

    ! add particles to each weight group
    total_so4_vol = 0d0
    do i_group = 1,aero_weight_array_n_group(aero_state%awa)
       ! adjust weight if necessary
       n_samp_avg = nucleate_rate * del_t / aero_weight_num_conc_at_radius( &
            aero_state%awa%weight(i_group, i_class), diam2rad(nucleate_diam))
       call aero_state_prepare_weight_for_add(aero_state, aero_data, &
            i_group, i_class, n_samp_avg, allow_doubling, allow_halving)

       ! determine number of nucleated particles
       n_samp_avg = nucleate_rate * del_t / aero_weight_num_conc_at_radius( &
            aero_state%awa%weight(i_group, i_class), diam2rad(nucleate_diam))
       n_samp = rand_poisson(n_samp_avg)

       ! create the particles
       do i_samp = 1,n_samp
          so4_vol = aero_data_diam2vol(aero_data, nucleate_diam)
          total_so4_vol = total_so4_vol + so4_vol

          call aero_particle_zero(aero_particle, aero_data)
          call aero_particle_set_create_time(aero_particle, &
               env_state%elapsed_time)
          aero_particle%vol(i_aero_so4) = so4_vol
          call aero_particle_new_id(aero_particle)
          call aero_particle_set_weight(aero_particle, i_group, i_class)
          call aero_state_add_particle(aero_state, aero_particle, aero_data)
       end do
    end do

    ! remove gases that formed new particles
    h2so4_removed_conc = &
         total_so4_vol &                        ! volume of SO4
         * aero_weight_array_num_conc_at_radius(aero_state%awa, i_class, &
         diam2rad(nucleate_diam)) &             ! volume conc of SO4
         * aero_data%density(i_aero_so4) &      ! mass conc of SO4
         / aero_data%molec_weight(i_aero_so4) & ! mole conc of SO4
         * const%avagadro                       ! molecule conc of SO4
    gas_state%mix_rat(i_gas_h2so4) = gas_state%mix_rat(i_gas_h2so4) &
         - env_state_conc_to_ppb(env_state, h2so4_removed_conc)
    if (gas_state%mix_rat(i_gas_h2so4) < 0d0) then
       gas_state%mix_rat(i_gas_h2so4) = 0d0
    end if

  end subroutine nucleate_sulf_acid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_file_read_nucleate_type(file, aero_data, nucleate_type, &
       nucleate_source)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aerosol data.
    type(aero_data_t), intent(inout) :: aero_data
    !> Nucleate type.
    integer, intent(out) :: nucleate_type
    !> Nucleate source number.
    integer, intent(out) :: nucleate_source

    character(len=SPEC_LINE_MAX_VAR_LEN) :: nucleate_type_name

    !> \page input_format_nucleate Input File Format: Nucleation Parameterization
    !!
    !! The nucleation parameterization is specified by the parameter:
    !!   - \b nucleate (string): the type of nucleation
    !!     parameterization --- must be one of: "none" for no
    !!     nucleation; or "sulf_acid" for the nucleate_sulf_acid()
    !!     parameterization
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    call spec_file_read_string(file, 'nucleate', nucleate_type_name)
    if (nucleate_type_name == 'sulf_acid') then
       nucleate_type = NUCLEATE_TYPE_SULF_ACID
       nucleate_source = aero_data_source_by_name(aero_data, &
            NUCLEATE_SOURCE_NAME)
    else
       call spec_file_die_msg(707263678, file, "unknown nucleate type: " &
            // trim(nucleate_type_name))
    end if

  end subroutine spec_file_read_nucleate_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_nucleate
