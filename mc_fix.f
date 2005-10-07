C mc_fix.f
C
C Monte Carlo with fixed timestep

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mc_fix(MM, M, M_comp, V, V_comp, kernel, n_bin, vv,
     &     rr, dp, dlnr, t_max, del_t)

      integer MM        ! INPUT: physical dimension of V
      integer M         ! INPUT/OUTPUT: number of particles
      integer M_comp    ! INPUT/OUTPUT: logical dimension of V
      real*8 V(MM)      ! INPUT/OUTPUT: particle volumes
      real*8 V_comp     ! INPUT/OUTPUT: computational volume
      external kernel   ! INPUT: kernel procedure
      integer n_bin     ! INPUT: number of bins
      real*8 vv(n_bin)  ! INPUT: volume of bins
      real*8 rr(n_bin)  ! INPUT: radius of bins
      real*8 dp(n_bin)  ! INPUT: diameter of bins
      real*8 dlnr       ! INPUT: scale factor
      real*8 t_max      ! INPUT: total simulation time
      real*8 del_t      ! INPUT: timestep

      integer i_top, nt, n_samp, i_samp
      real*8 g(n_bin), n_ln(n_bin), k_max, time

      time = 0.
      M_comp = M

      call moments(MM, V, n_bin, M_comp, V_comp, vv, dlnr, g, n_ln)
      call print_info(n_bin, time, dp, g, n_ln)

      nt = t_max / del_t
      do i_top = 1,nt
         time = real(i_top) / real(nt) * t_max
         
         call coagmax(n_bin, rr, n_ln, dlnr, k_max)
         call compute_n_samp(M, k_max, V_comp, del_t, n_samp)
         do i_samp = 1,n_samp
            call maybe_coag_pair(V, MM, M, M_comp, V_comp,
     &           del_t, n_samp)
            if (M .lt. MM / 2) then
               call double(MM, M, M_comp, V, V_comp)
            endif
         enddo

         call moments(MM, V, n_bin, M_comp, V_comp, vv, dlnr, g, n_ln)
         call print_info(n_bin, time, dp, g, n_ln)
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine compute_n_samp(M, k_max, V_comp,
     &     del_T, n_samp)

      integer M            ! INPUT: number of particles
      real*8 k_max  ! INPUT: maximum kernel value
      real*8 V_comp        ! INPUT: computational volume
      real*8 del_T         ! INPUT: timestep
      integer n_samp       ! OUTPUT: number of samples to take

      real*8 p_max, r_samp
      parameter (p_max = 0.01)

      r_samp = - (k_max * 1/V_comp *del_T/log(1-p_max))
      n_samp = r_samp * M*(M-1)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     
      subroutine maybe_coag_pair(V, MM, M, M_comp, V_comp,
     &     del_T, n_samp)
      
      integer n_samp, M, MM, M_comp ! INPUT
      real*8 V(MM), V_comp, del_T   ! INPUT

      integer s1, s2
      real*8 expo, p, k

      call find_rand_pair(MM, V, M_comp, s1, s2) ! test particles s1, s2
      call kernel_sedi(V(s1), V(s2), k)
      expo = k * 1.0/V_comp * del_T * M*(M-1)/n_samp
      p = 1 - exp(-expo) ! probability of coagulation
      if (rand() .lt. p) call coagulate(MM, M, V, s1, s2)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
