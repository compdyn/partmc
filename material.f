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

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function average_solute_quantity(mat, quantity)

    ! returns the average of the non-water elements of quantity

    type(material), intent(in) :: mat    ! material properties
    real*8, dimension(:), intent(in) :: quantity       ! quantity to average

    real*8 total
    integer i

    total = 0d0
    do i = 1,mat%n_spec
       if (i .ne. mat%i_water) then
          total = total + quantity(i)
       end if
    end do
    average_solute_quantity = total / dble(mat%n_spec - 1)

  end function average_solute_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function average_solute_material_weight(mat)

    ! returns the average of the non-water molecular weights

    type(material), intent(in) :: mat    ! material properties

    average_solute_material_weight = average_solute_quantity(mat, mat%M_w)

  end function average_solute_material_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_material
