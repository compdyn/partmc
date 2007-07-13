! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Useful physical constants. These should all be absolute constants
! and environment parameters that in principle could change should go
! in environ.f.
!
! To access a constant in a subroutine you should "use mod_constants"
! and then the constant value is accessed with const%pi or
! similar. Note that the type is called consts (with a trailing s) but
! the single saved variable to access them is called const (without a
! trailing s).

module mod_constants

  type const_t
     real*8 :: pi = 3.14159265358979323846d0
     real*8 :: p00 = 611d0              ! equil. water vapor press at 273 K (Pa)
     real*8 :: T0 = 273.15d0            ! freezing point of water (K)
     real*8 :: sig = 0.073d0            ! surface energy of water (J m^{-2})
     real*8 :: R = 8.314472d0           ! univ. gas const (J mole^{-1} K^{-1})
     real*8 :: L_v = 2.272d6            ! latent heat of water (J kg^{-1})
     real*8 :: alpha = 1d0              ! accomodation coeff (0.045 also used)
     real*8 :: cp = 1005d0              ! spec. heat of water (J kg^{-1} K^{-1})
     real*8 :: M_a = 2.89644d-2         ! molecular weight of air (kg mole^{-1})
     real*8 :: atm = 101325d0           ! atm. standard sea level pressure (Pa)
     real*8 :: k_b = 1.3806505d-23      ! Boltzmann constant in J K^{-1}
     real*8 :: mu = 1.78d-5             ! dynamic visc. air (kg m^{-1} s^{-1})
     real*8 :: N_A = 6.02214179d23      ! Avogadro's number (mole^{-1})
  end type const_t

  type(const_t), save :: const

end module mod_constants
