C     Functions to deal with a particle array VH that is stored by
C     bin. VH(i_bin,i) is the i-th particle in the i_bin-th bin.

C     FIXME: MH and bin_n are pretty much identical. Probably best to
C     ignore it, for symmetry with non-hybrid code, and because in
C     superparticle code there is a difference.

      module mod_array_hybrid
      contains

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine array_to_hybrid(MM, M, V, n_spec, n_bin, bin_v, TDV, MH
     $     ,VH)

C     Convert a standard linear particle array V to a hybrid particle
C     array VH stored by bins.

      use mod_array

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT: logical dimension of V
      integer n_spec       ! INPUT: number of species
      real*8 V(MM,n_spec)  ! INPUT/OUTPUT: particle volumes
      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      integer TDV          ! INPUT: trailing dimension of VH
      integer MH(n_bin)    ! OUTPUT: number of particles per bin
      real*8 VH(n_bin,TDV,n_spec) ! OUTPUT: particle volumes in hybrid array

      integer i, j, k
      real*8 pv

      do k = 1,n_bin
         MH(k) = 0
      enddo

      write(6,*)'in array ',M,n_spec
      do i = 1,M
         call particle_vol(MM,n_spec,V,i,pv)     
         call particle_in_bin(pv, n_bin, bin_v, k)
         MH(k) = MH(k) + 1
         if (MH(k) .gt. TDV) then
            write(*,*)'ERROR: TDV too small'
            call exit(2)
         endif
         do j=1,n_spec
            VH(k, MH(k),j) = V(i,j)
         enddo
      enddo

      do k=1,n_bin
         write(6,*)'k,MH(k) ',k,MH(k)
      enddo

      write(6,*)'array ende '
      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine maybe_coag_pair_hybrid(M, n_bin, TDV, MH, VH, V_comp,
     $     n_spec, bin_v, bin_r, bin_g, bin_gs, bin_n, dlnr, b1, b2,
     $     del_t, k_max, kernel, did_coag, bin_change)

C     Choose a random pair for potential coagulation and test its
C     probability of coagulation. If it happens, do the coagulation and
C     update all structures. The probability of a coagulation will be
C     taken as (kernel / k_max).

      integer M            ! INPUT/OUTPUT: number of particles
      integer n_bin        ! INPUT: number of bins
      integer TDV          ! INPUT: trailing dimension of VH
      integer n_spec       ! INPUT: number of species
      integer MH(n_bin)    ! INPUT/OUTPUT: number of particles per bin
      real*8 VH(n_bin,TDV,n_spec) ! INPUT/OUTPUT: particle volumes
      real*8 V_comp        ! INPUT: computational volume

      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins
      real*8 bin_gs(n_bin,n_spec)  ! INPUT/OUTPUT: species mass in bins
      integer bin_n(n_bin) ! INPUT/OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor

      integer b1           ! INPUT: bin of first particle
      integer b2           ! INPUT: bin of second particle
      real*8 del_t         ! INPUT: timestep
      real*8 k_max         ! INPUT: k_max scale factor
C      external kernel      ! INPUT: kernel function
      logical did_coag     ! OUTPUT: whether a coagulation occured
      logical bin_change   ! OUTPUT: whether bin structure changed

      interface
         subroutine kernel(v1, v2, k)
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         real*8, intent(out) :: k
         end subroutine
      end interface

      integer s1, s2
      real*8 p, k, pv1, pv2

      bin_change = .false.
      did_coag = .false.

      if ((MH(b1) .le. 0) .or. (MH(b2) .le. 0)) then
         return
      endif

      call find_rand_pair_hybrid(n_bin, MH, b1, b2, s1, s2)
      call particle_vol_hybrid(n_bin,TDV,n_spec,VH,b1,s1,pv1)
      call particle_vol_hybrid(n_bin,TDV,n_spec,VH,b2,s2,pv2)
      call kernel(pv1, pv2, k)
      p = k / k_max

      if (dble(rand()) .lt. p) then
         call coagulate_hybrid(M, n_bin, TDV, MH, VH, V_comp, n_spec
     $        ,bin_v,bin_r,bin_g, bin_gs, bin_n, dlnr, b1, s1, b2, s2,
     $        bin_change)
         did_coag = .true.
      endif

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     
      subroutine find_rand_pair_hybrid(n_bin, MH, b1, b2, s1, s2)

C     Find a random pair of particles (b1, s1) and (b2, s2).
      
      integer n_bin     ! INPUT: number of bins
      integer MH(n_bin) ! INPUT: number particles per bin
      integer b1        ! INPUT: bin number of first particle
      integer b2        ! INPUT: bin number of second particle
      integer s1        ! OUTPUT: first random particle 1 <= s1 <= M(b1)
      integer s2        ! OUTPUT: second random particle 1 <= s2 <= M(b2)
                        !         (b1,s1) != (b2,s2)

 100  s1 = int(rand() * MH(b1)) + 1
      if ((s1 .lt. 1) .or. (s1 .gt. MH(b1))) goto 100
 101  s2 = int(rand() * MH(b2)) + 1
      if ((s2 .lt. 1) .or. (s2 .gt. MH(b2))) goto 101
      if ((b1 .eq. b2) .and. (s1 .eq. s2)) goto 101

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine coagulate_hybrid(M, n_bin, TDV, MH, VH, V_comp, n_spec
     $     ,bin_v,bin_r, bin_g, bin_gs,bin_n, dlnr, b1, s1, b2, s2,
     $     bin_change)

C     Join together particles (b1, s1) and (b2, s2), updating all
C     particle and bin structures to reflect the change. bin_change is
C     true if the used bin set changed due to the coagulation (i.e. an
C     empty bin filled or a filled bin became empty).

      use mod_array

      integer M            ! INPUT/OUTPUT: number of particles
      integer n_bin        ! INPUT: number of bins
      integer TDV          ! INPUT: trailing dimension of VH
      integer n_spec       ! INPUT: number of species
      integer MH(n_bin)    ! INPUT/OUTPUT: number of particles per bin
      real*8 VH(n_bin,TDV,n_spec) ! INPUT/OUTPUT: particle volumes
      real*8 V_comp        ! INPUT: computational volume

      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins
      real*8 bin_gs(n_bin,n_spec)  ! INPUT/OUTPUT: species mass in bins
      integer bin_n(n_bin) ! INPUT/OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor
      
      integer b1           ! INPUT: first particle (bin number)
      integer s1           ! INPUT: first particle (number in bin)
      integer b2           ! INPUT: second particle (bin number)
      integer s2           ! INPUT: second particle (number in bin)
      logical bin_change   ! OUTPUT: whether an empty bin filled,
                           !         or a filled bin became empty

      integer bn, i, j
      real*8 new_v(n_spec), pv1, pv2, new_v_tot

      bin_change = .false.
      new_v_tot = 0.d0

      call particle_vol_hybrid(n_bin,TDV,n_spec,VH,b1,s1,pv1)
      call particle_vol_hybrid(n_bin,TDV,n_spec,VH,b2,s2,pv2)

      ! remove s1 and s2 from bins
      bin_n(b1) = bin_n(b1) - 1
      bin_n(b2) = bin_n(b2) - 1
      bin_g(b1) = bin_g(b1) - pv1
      bin_g(b2) = bin_g(b2) - pv2
      do j=1,n_spec
         bin_gs(b1,j) = bin_gs(b1,j) - VH(b1,s1,j)
         bin_gs(b2,j) = bin_gs(b2,j) - VH(b2,s2,j)
      enddo
      if ((bin_n(b1) .lt. 0) .or. (bin_n(b2) .lt. 0)) then
         write(*,*)'ERROR: invalid bin_n'
         call exit(2)
      endif

      ! do coagulation in MH, VH arrays
      do i=1,n_spec
         new_v(i) = VH(b1,s1,i) + VH(b2,s2,i)   ! add particle volumes
      enddo

      do i=1,n_spec
         new_v_tot = new_v_tot + new_v(i)
      enddo

      call particle_in_bin(new_v_tot, n_bin, bin_v, bn)  ! find new bin

      do i=1,n_spec
         VH(b1, s1,i) = VH(b1, MH(b1),i) ! shift last particle into empty slot
      enddo

      MH(b1) = MH(b1) - 1          ! decrease length of array
      do i=1,n_spec
         VH(b2, s2,i) = VH(b2, MH(b2),i) ! same for second particle
      enddo
      MH(b2) = MH(b2) - 1
      if ((MH(b1) .lt. 0) .or. (MH(b2) .lt. 0)) then
         write(*,*)'ERROR: invalid MH'
         call exit(2)
      endif
      MH(bn) = MH(bn) + 1          ! increase the length of array
      do i=1,n_spec
         VH(bn, MH(bn),i) = new_v(i) ! add the new particle at the end
      enddo
      M = M - 1                    ! decrease the total number of particles

      ! add new particle to bins
      bin_n(bn) = bin_n(bn) + 1
      do i=1,n_spec
         bin_g(bn) = bin_g(bn) + VH(bn, MH(bn),i)
         bin_gs(bn,i) = bin_gs(bn,i)+VH(bn,MH(bn),i)
      enddo

      ! did we empty a bin?
      if ((bin_n(b1) .eq. 0) .or. (bin_n(b2) .eq. 0))
     &     bin_change = .true.
      ! did we start a new bin?
      if ((bin_n(bn) .eq. 1) .and. (bn .ne. b1) .and. (bn .ne. b2))
     &     bin_change = .true.

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine double_hybrid(M, n_bin, TDV, MH, VH, V_comp, n_spec
     $     ,bin_v,bin_r, bin_g, bin_gs, bin_n, dlnr)

C     Double number of particles in a hybrid array.

      integer M            ! INPUT/OUTPUT: number of particles
      integer n_bin        ! INPUT: number of bins
      integer TDV          ! INPUT: trailing dimension of VH
      integer n_spec       ! INPUT: number of species
      integer MH(n_bin)    ! INPUT/OUTPUT: number of particles per bin
      real*8 VH(n_bin,TDV,n_spec) ! INPUT/OUTPUT: particle volumes
      real*8 V_comp        ! INPUT/OUTPUT: computational volume

      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins
      real*8 bin_gs(n_bin,n_spec) ! INPUT/OUTPUT: species mass in bins
      integer bin_n(n_bin) ! INPUT/OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor

      integer i, j, k, i_spec

      ! double VH and associated structures
      do k = 1,n_bin
         ! only double if we have enough space to do so
         if (2 * MH(k) .gt. TDV) then
            write(*,*)'ERROR: double without enough space'
            call exit(2)
         endif
         do i = 1,MH(k)
            do i_spec = 1,n_spec
               VH(k, i + MH(k),i_spec) = VH(k, i, i_spec)
            enddo
         enddo
         MH(k) = 2 * MH(k)
      enddo
      M = 2 * M
      V_comp = 2d0 * V_comp

      ! double bin structures
      do i = 1,n_bin
         bin_g(i) = bin_g(i) * 2d0
         bin_n(i) = bin_n(i) * 2
         do j=1,n_spec
            bin_gs(i,j) = bin_gs(i,j) * 2d0
         enddo
      enddo

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine check_hybrid(M, n_bin, n_spec, TDV, MH, VH, bin_v,
     $     bin_r)

C     Check that V has all particles in the correct bins.

      use mod_array

      integer M            ! INPUT: number of particles
      integer n_bin        ! INPUT: number of bins
      integer TDV          ! INPUT: trailing dimension of VH
      integer n_spec       ! INPUT: number of species
      integer MH(n_bin)    ! INPUT: number of particles per bin
      real*8 VH(n_bin,TDV,n_spec) ! INPUT: particle volumes

      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins

      real*8 pv
      integer i, k, k_check, M_check
      logical error

      error = .false.
      M_check = 0
      do k = 1,n_bin
         M_check = M_check + MH(k)
         do i = 1,MH(k)
            call particle_vol_hybrid(n_bin,TDV,n_spec,VH,k,i,pv)
            call particle_in_bin(pv, n_bin, bin_v, k_check)
            if (k .ne. k_check) then
               write(*,'(a10,a10,a12,a10)') 'k', 'i', 'VH(k, i)',
     $              'k_check'
               write(*,'(i10,i10,e12.5,i10)') k, i, pv, k_check
               error = .true.
            endif
         enddo
      enddo
      if (M .ne. M_check) then
         write(*,'(a10,a10)') 'M', 'M_check'
         write(*,'(i10,i10)') M, M_check
         error = .true.
      endif
      if (error) then
         write(*,*)'ERROR: check_hybrid() failed'
         call exit(2)
      endif

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine particle_vol_hybrid(n_bin,TDV,n_spec,VH,k,i,pv)

      integer n_spec       ! INPUT: number of species
      integer n_bin        ! INPUT: number of bins
      integer TDV          ! INPUT: trailing dimension of VH      
      real*8 VH(n_bin,TDV,n_spec)  ! INPUT: particle volumes (m^3)
      integer i            ! INPUT: particle index
      integer k            ! INPUT: bin index
      real*8 pv            ! OUPUT: total volume of particle

      integer j

      pv = 0d0
cn this is maybe not quite right! 
      do j=1,n_spec
         pv = pv + VH(k,i,j)
      enddo
      return
      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      end module
