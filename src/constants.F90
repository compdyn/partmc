! Copyright (C) 2005-2009, 2012, 2016 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_constants module.

!> Physical constants.
module pmc_constants

  !> Kind of a double precision real number.
  integer, parameter :: dp = kind(0.d0)
  !> Kind of a double precision complex number.
  integer, parameter :: dc = dp

  !> Kind of an integer
  integer, parameter :: i_kind = kind(1)

  !> Physical constants.
  !!
  !! These are all absolute constants. Environment parameters that in
  !! principle could change should go in scenario.f90.
  !!
  !! To access a constant in a subroutine you should <tt>use
  !! pmc_constants</tt> and then the constant value is accessed with
  !! \c const%%pi or similar. Note that the type is called \c const_t
  !! (with a trailing _t) but the single saved variable to access them
  !! is just called \c const.
  type const_t
     !> Pi.
     real(kind=dp) :: pi = 3.14159265358979323846d0
     !> Boltzmann constant (J K^{-1}).
     real(kind=dp) :: boltzmann = 1.3806505d-23
     !> Avogadro's number (mole^{-1}).
     real(kind=dp) :: avagadro = 6.02214179d23
     !> Universal gas constant (J mole^{-1} K^{-1}).
     real(kind=dp) :: univ_gas_const = 8.314472d0
     !> Acceleration due to gravity (m s^{-2}).
     real(kind=dp) :: std_grav = 9.80665d0
     !> Accomodation coefficient (have also used 0.045).
     real(kind=dp) :: accom_coeff = 1d0

     !> Equilibrium water vapor pressure at 273 K (Pa).
     real(kind=dp) :: water_eq_vap_press = 611d0
     !> Freezing point of water (K).
     real(kind=dp) :: water_freeze_temp = 273.15d0
     !> Surface energy of water (J m^{-2}).
     real(kind=dp) :: water_surf_eng = 0.073d0
     !> Latent heat of water (J kg^{-1}).
     real(kind=dp) :: water_latent_heat = 2.272d6
     !> Specific heat of air (J kg^{-1} K^{-1}).
     real(kind=dp) :: air_spec_heat = 1005d0
     !> Molecular weight of water (kg mole^{-1}).
     real(kind=dp) :: water_molec_weight = 18d-3
     !> Density of water (kg m^{-3}).
     real(kind=dp) :: water_density = 1d3

     !> Molecular weight of air (kg mole^{-1}).
     real(kind=dp) :: air_molec_weight = 2.89644d-2
     !> Atmospheric standard sea level pressure (Pa).
     real(kind=dp) :: air_std_press = 101325d0
     !> Dynamic viscosity of air (kg m^{-1} s^{-1}).
     real(kind=dp) :: air_dyn_visc = 1.78d-5
  end type const_t

  !> Fixed variable for accessing the constant's values.
  !!
  !! To access a constant in a subroutine you should <tt>use
  !! pmc_constants</tt> and then the constant value is accessed with
  !! \c const%%pi or similar. Note that the type is called \c const_t
  !! (with a trailing _t) but the single saved variable to access them
  !! is just called \c const.
  type(const_t), save :: const

end module pmc_constants
