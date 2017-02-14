! Copyright (C) 2011-2012, 2016-2017 Jian Tian, Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_chamber module.

module pmc_chamber

  use pmc_aero_data
  use pmc_constants
  use pmc_env_state
  use pmc_spec_file

  !> Unit translational diffusion coefficient (m^2 s^{-1}).
  real(kind=dp), parameter :: CHAMBER_UNIT_DIFF_COEF = 1d0

  type chamber_t
     !> Chamber volume (m^3).
     real(kind=dp) :: volume = 0d0
     !> Diffusional deposition area (m^2).
     real(kind=dp) :: area_diffuse = 0d0
     !> Sedimentational deposition area (m^2).
     real(kind=dp) :: area_sedi = 0d0
     !> Prefactor in dissusive boundary layer thickness (m).
     real(kind=dp) :: prefactor_BL = 0d0
     !> Exponent in dissusive boundary layer thickness.
     real(kind=dp) :: exponent_BL = 0d0
  end type chamber_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Based on Eq. 23 in Naumann 2003 J. Aerosol. Sci.
  real(kind=dp) function chamber_diff_coef(vol, aero_data, temp, pressure)

    !> Particle volume (m^3).
    real(kind=dp), intent(in) :: vol
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: pressure

    real(kind=dp) :: R_eff, R_me_c, C

    R_eff = fractal_vol_to_effective_rad(aero_data%fractal, vol)
    R_me_c = fractal_vol_to_mobility_rad_in_continuum(aero_data%fractal, vol)
    C = fractal_slip_correct(R_eff, temp, pressure)

    chamber_diff_coef = (const%boltzmann * temp * C) &
         / (6d0 * const%pi * const%air_dyn_visc * R_me_c)

  end function chamber_diff_coef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate diffusional boundary layer thickness.
  !> Based on Eq. 40 in Naumann 2003 J. Aerosol. Sci.
  real(kind=dp) function chamber_diff_BL_thick(chamber, vol, aero_data, temp, &
       pressure)

    !> Chamber parameters.
    type(chamber_t) :: chamber
    !> Particle volume (m^3).
    real(kind=dp), intent(in) :: vol
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: pressure

    real(kind=dp) :: D

    D = chamber_diff_coef(vol, aero_data, temp, pressure)
    
    chamber_diff_BL_thick = chamber%prefactor_BL &
         * (D / CHAMBER_UNIT_DIFF_COEF)**chamber%exponent_BL

  end function chamber_diff_BL_thick

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the loss rate due to wall diffusion in chamber.
  !> Based on Eq. 37a in Naumann 2003 J. Aerosol. Sci.
  real(kind=dp) function chamber_loss_rate_wall(chamber, vol, aero_data, &
       env_state)

    !> Chamber parameters.
    type(chamber_t), intent(in) :: chamber
    !> Particle volume (m^3).
    real(kind=dp), intent(in) :: vol
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    real(kind=dp) :: D, delta

    D = chamber_diff_coef(vol, aero_data, env_state%temp, env_state%pressure)
    delta = chamber_diff_BL_thick(chamber, vol, aero_data, env_state%temp, &
         env_state%pressure)

    chamber_loss_rate_wall = (D * chamber%area_diffuse) &
         / (delta * chamber%volume)

  end function chamber_loss_rate_wall

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the loss rate due to sedimentation in chamber.
  !> Based on Eq. 37b in Naumann 2003 J. Aerosol. Sci.
  real(kind=dp) function chamber_loss_rate_sedi(chamber, vol, density, &
       aero_data, env_state)

    !> Chamber parameters.
    type(chamber_t), intent(in) :: chamber
    !> Particle volume (m^3).
    real(kind=dp), intent(in) :: vol
    !> Particle density (kg/m^3).
    real(kind=dp), intent(in) :: density
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    real(kind=dp) :: D

    D = chamber_diff_coef(vol, aero_data, env_state%temp, env_state%pressure)

    chamber_loss_rate_sedi &
         = (density * vol * const%std_grav * D * chamber%area_sedi) &
         / (const%boltzmann * env_state%temp * chamber%volume)

  end function chamber_loss_rate_sedi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read chamber specification from a spec file.
  subroutine spec_file_read_chamber(file, chamber)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Chamber data.
    type(chamber_t), intent(inout) :: chamber

    !> \page input_format_chamber Input File Format: Chamber
    !!
    !! The chamber model is specified by the parameters:
    !! - \b chamber_vol (real, unit m^3): the volume of the chamber
    !! - \b area_diffuse (real, unit m^2): the surface area in the chamber
    !!   available for wall diffusion deposition (the total surface area)
    !! - \b area_sedi (real, unit m^2): the surface area in the chamber
    !!   available for sedimentation deposition (the floor area)
    !! - \b prefactor_BL (real, unit m): the coefficient \f$k_{\rm D}\f$ in
    !!   the model \f$ \delta = k_{\rm D}(D/D_0)^a \f$ for boundary-layer
    !!   thickness \f$ \delta \f$
    !! - \b exponent_BL (real, dimensionless): the exponent \f$a\f$ in
    !!   the model \f$ \delta = k_{\rm D}(D/D_0)^a \f$ for boundary-layer
    !!   thickness \f$ \delta \f$
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_scenario --- the prescribed profiles of
    !!     other environment data

    call spec_file_read_real(file, 'chamber_vol', &
         chamber%volume)
    call spec_file_read_real(file, 'area_diffuse', &
         chamber%area_diffuse)
    call spec_file_read_real(file, 'area_sedi', &
         chamber%area_sedi)
    call spec_file_read_real(file, 'prefactor_BL', &
         chamber%prefactor_BL)
    call spec_file_read_real(file, 'exponent_BL', &
         chamber%exponent_BL)

  end subroutine spec_file_read_chamber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_chamber(val)

    !> Value to pack.
    type(chamber_t), intent(in) :: val

    integer :: total_size, i, n

    pmc_mpi_pack_size_chamber = &
         pmc_mpi_pack_size_real(val%volume) &
         + pmc_mpi_pack_size_real(val%area_diffuse) &
         + pmc_mpi_pack_size_real(val%area_sedi) &
         + pmc_mpi_pack_size_real(val%prefactor_BL) &
         + pmc_mpi_pack_size_real(val%exponent_BL)

  end function pmc_mpi_pack_size_chamber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_chamber(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(chamber_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_pack_real(buffer, position, val%volume)
    call pmc_mpi_pack_real(buffer, position, val%area_diffuse)
    call pmc_mpi_pack_real(buffer, position, val%area_sedi)
    call pmc_mpi_pack_real(buffer, position, val%prefactor_BL)
    call pmc_mpi_pack_real(buffer, position, val%exponent_BL)
    call assert(623603534, &
         position - prev_position <= pmc_mpi_pack_size_chamber(val))
#endif

  end subroutine pmc_mpi_pack_chamber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_chamber(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(chamber_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_unpack_real(buffer, position, val%volume)
    call pmc_mpi_unpack_real(buffer, position, val%area_diffuse)
    call pmc_mpi_unpack_real(buffer, position, val%area_sedi)
    call pmc_mpi_unpack_real(buffer, position, val%prefactor_BL)
    call pmc_mpi_unpack_real(buffer, position, val%exponent_BL)
    call assert(986998595, &
         position - prev_position <= pmc_mpi_pack_size_chamber(val))
#endif

  end subroutine pmc_mpi_unpack_chamber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_chamber
