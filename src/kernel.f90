! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Generic kernel functions.

module mod_kernel

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine bin_kernel(n_bin, bin_v, kernel, env, k)
   
    use mod_environ
    
    ! Computes an array of kernel values for each bin pair. k(i,j) is
    ! the kernel value at the midpoint of bins i and j.
    
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    real*8, intent(out) :: k(n_bin,n_bin) ! kernel values
    type(environ), intent(in) :: env    ! environment state

    interface
       subroutine kernel(v1, v2, env, k)
         use mod_environ
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(environ), intent(in) :: env
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
    
    integer i, j
    
    do i = 1,n_bin
       do j = 1,n_bin
          call kernel(bin_v(i), bin_v(j), env, k(i,j))
       end do
    end do
    
  end subroutine bin_kernel
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine est_k_max_binned(n_bin, bin_v, kernel, env, k_max)

    ! Computes an array of maximum kernel values. Given particles v1
    ! in bin b1 and v2 in bin b2, it is approximately true that
    ! kernel(v1,v2) <= k_max(b1,b2).

    use mod_environ
    
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    type(environ), intent(in) :: env  ! environment state
    real*8, intent(out) :: k_max(n_bin,n_bin) ! maximum kernel values
    
    interface
       subroutine kernel(v1, v2, env, k)
         use mod_environ
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(environ), intent(in) :: env
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
    
    integer i, j
    
    do i = 1,n_bin
       do j = 1,n_bin
          call est_k_max_for_bin(n_bin, bin_v, kernel, i, j, &
               env, k_max(i,j))
       end do
    end do
    
  end subroutine est_k_max_binned
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine est_k_max_for_bin(n_bin, bin_v, kernel, b1, b2, env, k_max)

    ! Samples within bins b1 and b2 to find the maximum value of the
    ! kernel between particles from the two bins.
   
    use mod_environ
    use mod_bin
 
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    integer, intent(in) :: b1           ! first bin
    integer, intent(in) :: b2           ! second bin
    type(environ), intent(in) :: env    ! environment state    
    real*8, intent(out) :: k_max        ! maximum kernel values
    
    interface
       subroutine kernel(v1, v2, env, k)
         use mod_environ
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(environ), intent(in) :: env
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
    
    real*8 v1, v2, v1_high, v1_low, v2_high, v2_low, k
    integer i, j
    
    integer, parameter :: n_sample = 10  ! number of sample points per bin
    
    ! v1_low < bin_v(b1) < v1_high
    call bin_edge(n_bin, bin_v, b1, v1_low)
    call bin_edge(n_bin, bin_v, b1 + 1, v1_high)
    
    ! v2_low < bin_v(b2) < v2_high
    call bin_edge(n_bin, bin_v, b2, v2_low)
    call bin_edge(n_bin, bin_v, b2 + 1, v2_high)
    
    k_max = 0d0
    do i = 1,n_sample
       do j = 1,n_sample
          v1 = v1_high * dble(n_sample - i) / dble(n_sample - 1) + &
               v1_low * dble(i - 1) / dble(n_sample - 1)
          v2 = v2_high * dble(n_sample - j) / dble(n_sample - 1) + &
               v2_low * dble(j - 1) / dble(n_sample - 1)
          call kernel(v1, v2, env, k)
          if (k .gt. k_max) k_max = k
       end do
    end do
    
  end subroutine est_k_max_for_bin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_kernel
