! -*- mode: f90; -*-
! Copyright (C) 2005,2006 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Exact solution output.

module mod_mc_exact
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine mc_exact(n_bin, n_spec, bin_v, bin_g, bin_gs, &
       bin_n, dlnr, N_0, V_0, rho_p, soln, t_max, t_print, loop, &
       V_comp, env, mat)
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
    real*8, intent(in) :: t_print       ! interval to print info (seconds)
    integer, intent(in) :: loop         ! loop number of run
    real*8, intent(in) :: V_comp        ! computational volume
    type(environ), intent(inout) :: env  ! environment state
    type(material), intent(in) :: mat    ! material properties
    
    integer i_time, n_time,i
    real*8 time
    
    interface
       subroutine soln(n_bin, bin_v, &
            bin_g, bin_n, dlnr, &
            time, N_0, V_0, rho_p, V_comp, env)

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
         real*8, intent(in) :: V_comp             !  computational volume
         type(environ), intent(in) :: env         ! environment state
       end subroutine soln
    end interface
    
    n_time = int(t_max / t_print)
    do i_time = 0,n_time
       time = dble(i_time) / dble(n_time) * dble(t_max)
       call soln(n_bin, bin_v, bin_g, bin_n, dlnr, &
            time, N_0, V_0, rho_p, V_comp, env)
       
       do i=1,n_bin
          bin_gs(i,1) = bin_g(i)
       enddo
       
       call print_info(time, V_comp,n_spec, &
            n_bin, bin_v, bin_g, bin_gs,bin_n, dlnr, env, mat)
    enddo
    
  end subroutine mc_exact
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_mc_exact
