
module mod_constants

  real*8, parameter :: p00 = 611d0   ! equilibrium water vapor pressure
                                     ! at 273 K (Pa)
  real*8, parameter :: T0 = 273.15d0 ! freezing point of water (K)
  real*8, parameter :: sig = 0.073d0 ! surface energy (J m^{-2})
  real*8, parameter :: R = 8.314d0   ! universal gas constant (J mole^{-1} K^{-1})
  real*8, parameter :: L_v = 2.5d6   ! latent heat (J kg^{-1})
  real*8, parameter :: alpha = 1d0   ! accomodation coefficient (the value 0.045 is also used sometimes)
  real*8, parameter :: cp = 1005d0   ! specific heat of water (J kg^{-1} K^{-1})
  real*8, parameter :: M_a = 28d-3   ! molecular weight of air (kg mole^{-1})
  real*8, parameter :: rho_a = 1.25d0 ! air density (kg m^{-3})
  real*8, parameter :: pi = 3.14159265358979323846d0

end module mod_constants
