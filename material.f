! -*- mode: f90;-*-

module mod_material

  type material
     integer n_spec
     integer i_water ! water species number
     real*8, dimension(:), pointer ::  rho ! densities (kg m^{-3})
     integer, dimension(:), pointer :: nu ! number of ions in the solute
     real*8, dimension(:), pointer :: eps ! solubilities (1)
     real*8, dimension(:), pointer :: M_w ! molecular weights (kg mole^{-1})
  end type material

end module mod_material
