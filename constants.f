! -*- mode: f90; -*-
! Copyright (C) 2005,2006 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
! -*- mode: f90;-*-

module mod_constants

  type consts
     real*8 :: pi = 3.14159265358979323846d0
     real*8 :: p00 = 611d0   ! equilibrium water vapor pressure at 273 K (Pa)
     real*8 :: T0 = 273.15d0 ! freezing point of water (K)
     real*8 :: sig = 0.073d0 ! surface energy (J m^{-2})
     real*8 :: R = 8.314d0   ! universal gas constant (J mole^{-1} K^{-1})
     real*8 :: L_v = 2.5d6   ! latent heat (J kg^{-1})
     real*8 :: alpha = 1d0   ! accomodation coefficient
                             ! (0.045 is sometimes used)
     real*8 :: cp = 1005d0   ! specific heat of water (J kg^{-1} K^{-1})
     real*8 :: M_a = 28d-3   ! molecular weight of air (kg mole^{-1})
     real*8 :: rho_a = 1.25d0 ! air density (kg m^{-3})
     real*8 :: atm = 101325d0 ! atmospheric pressure (Pa)
  end type consts

  type(consts), save :: const

end module mod_constants
