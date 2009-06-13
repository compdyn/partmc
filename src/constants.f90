! Copyright (C) 2005-2008 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_constants module.

!> Physical constants.
module pmc_constants

  !> Physical constants.
  !!
  !! These are all absolute constants. Environment parameters that in
  !! principle could change should go in env_data.f90.
  !!
  !! To access a constant in a subroutine you should <tt>use
  !! pmc_constants</tt> and then the constant value is accessed with
  !! \c const%%pi or similar. Note that the type is called \c const_t
  !! (with a trailing _t) but the single saved variable to access them
  !! is just called \c const.
  type const_t
     !> Pi.
     real*8 :: pi = 3.14159265358979323846d0
     !> Boltzmann constant (J K^{-1}).
     real*8 :: boltzmann = 1.3806505d-23
     !> Avogadro's number (mole^{-1}).
     real*8 :: avagadro = 6.02214179d23
     !> Universal gas constant (J mole^{-1} K^{-1}).
     real*8 :: univ_gas_const = 8.314472d0
     !> Accomodation coefficient (have also used 0.045).
     real*8 :: accom_coeff = 1d0

     !> Equilibrium water vapor pressure at 273 K (Pa).
     real*8 :: water_eq_vap_press = 611d0
     !> Freezing point of water (K).
     real*8 :: water_freeze_temp = 273.15d0
     !> Surface energy of water (J m^{-2}).
     real*8 :: water_surf_eng = 0.073d0
     !> Latent heat of water (J kg^{-1}).
     real*8 :: water_latent_heat = 2.272d6
     !> Specific heat of water (J kg^{-1} K^{-1}).
     real*8 :: water_spec_heat = 1005d0
     !> Molecular weight of water (kg mole^{-1}).
     real*8 :: water_molec_weight = 18d-3
     !> Density of water (kg m^{-3}).
     real*8 :: water_density = 1d3

     !> Molecular weight of air (kg mole^{-1}).
     real*8 :: air_molec_weight = 2.89644d-2
     !> Atmospheric standard sea level pressure (Pa).
     real*8 :: air_std_press = 101325d0
     !> Dynamic viscosity of air (kg m^{-1} s^{-1}).
     real*8 :: air_dyn_visc = 1.78d-5
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
