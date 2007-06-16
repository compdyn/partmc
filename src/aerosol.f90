! Copyright (C) 2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

module mod_aerosol

  type bin_p
     real*8, dimension(:,:), pointer :: p ! particle volumes (m^3)
     ! dimension of p is (# particles in bin) x n_spec
  end type bin_p

  type aerosol
     integer, dimension(:), pointer :: MH  ! number of particles per bin
     type(bin_p), dimension(:), pointer :: VH ! particle volumes (m^3)
  end type aerosol

end module mod_aerosol
