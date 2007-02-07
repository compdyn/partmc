! -*- mode: f90; -*-
! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Exact solution output.

module mod_exact
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine exact(n_bin, n_spec, bin_v, bin_g, bin_gs, &
       bin_n, dlnr, N_0, V_0, rho_p, soln, t_max, t_output, &
       env, mat)
    ! FIXME: N_0 and V_0 are really parameters for the initial value
    ! of the particle distribution. They should be replaced by a n_param,
    ! params() pair.
    
    use mod_bin
    use mod_array
    use mod_environ
    use mod_material
    
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    real*8, intent(in) :: bin_v(n_bin)  ! volume of bins
    real*8, intent(out) :: bin_g(n_bin)  ! volume in bins
    integer, intent(out) :: bin_n(n_bin) ! number in bins
    real*8, intent(out) :: bin_gs(n_bin,n_spec) ! number in bins by species
    real*8, intent(in) :: dlnr          ! bin scale factor
    real*8, intent(in) :: N_0           ! particle number concentration (#/m^3)
    real*8, intent(in) :: V_0           ! 
    real*8, intent(in) :: rho_p         ! particle density (kg/m^3)
    real*8, intent(in) :: t_max         ! total simulation time
    real*8, intent(in) :: t_output      ! interval to print info (seconds)
    type(environ), intent(inout) :: env  ! environment state
    type(material), intent(in) :: mat    ! material properties
    
    integer i_time, n_time,i
    real*8 time
    
    interface
       subroutine soln(n_bin, bin_v, &
            bin_g, bin_n, dlnr, &
            time, N_0, V_0, rho_p, env)

         use mod_environ

         integer, intent(in) :: n_bin             !  number of bins
         real*8, intent(in) :: bin_v(n_bin)       !  volume of particles in bins
         real*8, intent(out) :: bin_g(n_bin)       !  volume in bins
         integer, intent(out) :: bin_n(n_bin)      !  number in bins
         real*8, intent(in) :: dlnr               !  bin scale factor
         
         real*8, intent(in) :: time               !  cubin_rent time
         real*8, intent(in) :: N_0                !  particle number concentration (#/m^3)
         real*8, intent(in) :: V_0                ! 
         real*8, intent(in) :: rho_p              !  particle density (kg/m^3)
         type(environ), intent(in) :: env         ! environment state
       end subroutine soln
    end interface
    
    n_time = int(t_max / t_output)
    do i_time = 0,n_time
       time = dble(i_time) / dble(n_time) * dble(t_max)
       call soln(n_bin, bin_v, bin_g, bin_n, dlnr, &
            time, N_0, V_0, rho_p, env)
       
       do i=1,n_bin
          bin_gs(i,1) = bin_g(i)
       end do
       
       call print_info(time, n_spec, &
            n_bin, bin_v, bin_g, bin_gs,bin_n, dlnr, env, mat)
    end do
    
  end subroutine exact
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_exact
