C     Functions to deal with a particle array V that is stored by
C     bin. VH(i_bin,i) is the i-th particle in the i_bin-th bin.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      call array_to_hybrid(MM, M, V, n_bin, bin_v, MH, VH)

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT: logical dimension of V
      real*8 V(MM)         ! INPUT/OUTPUT: particle volumes
      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      integer MH(n_bin)    ! OUTPUT: number of particles per bin
      real*8 VH(n_bin,MM)  ! OUTPUT: particle volumes in hybrid array

      integer i, k

      do k = 1,n_bin
         MH(k) = 0
      enddo
      do i = 1,M
         call particle_in_bin(V(i), n_bin, bin_v, k)
         MH(k) = MH(k) + 1
         VH(k, MH(k)) = V(i)
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine maybe_coag_pair(MM, M, n_bin, MH, VH, V_comp, bin_v,
     $     bin_r, bin_g, bin_n, dlnr, del_t, k_max, kernel, did_coag,
     $     bin_change)

      integer MM           ! INPUT: trailing dimension of V
      integer M            ! INPUT/OUTPUT: number of particles
      integer n_bin        ! INPUT: number of bins
      integer MH(n_bin)    ! INPUT/OUTPUT: number of particles per bin
      real*8 VH(n_bin,MM)  ! INPUT/OUTPUT: particle volumes
      real*8 V_comp        ! INPUT: computational volume

      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins
      integer bin_n(n_bin) ! INPUT/OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor
      
      real*8 del_t         ! INPUT: timestep
      real*8 k_max         ! INPUT: k_max scale factor
      external kernel      ! INPUT: kernel function
      logical did_coag     ! OUTPUT: whether a coagulation occured
      logical bin_change   ! OUTPUT: whether bin structure changed

      integer b1, s1, b2, s2
      real*8 p, k

      call find_rand_pair_hybrid(n_bin, MH, b1, s1, b2, s2)

      call kernel(VH(b1, s1), VH(b2, s2), k)
      p = k / k_max

      bin_change = .false.
      if (dble(rand()) .lt. p) then
         call coagulate_hybrid(MM, M, n_bin, MH, VH, V_comp, bin_v,
     $        bin_r,bin_g, bin_n, dlnr, b1, s1, b2, s2, bin_change)
         did_coag = .true.
      else
         did_coag = .false.
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     
      subroutine find_rand_pair_hybrid(n_bin, M, b1, b2, s1, s2)
      
      integer n_bin     ! INPUT: number of bins
      integer M(n_bin)  ! INPUT: number particles per bin
      integer b1        ! INPUT: bin number of first particle
      integer b2        ! INPUT: bin number of second particle
      integer s1        ! OUTPUT: first random particle 1 <= s1 <= M(b1)
      integer s2        ! OUTPUT: second random particle 1 <= s2 <= M(b2)
                        !         (b1,s1) != (b2,s2)

 100  s1 = int(rand() * M(b1)) + 1
      if ((s1 .lt. 1) .or. (s1 .gt. M(b1))) goto 100
 101  s2 = int(rand() * M(b2)) + 1
      if ((s2 .lt. 1) .or. (s2 .gt. M(b2))) goto 101
      if ((b1 .eq. b2) .and. (s1 .eq. s2)) goto 101

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine coagulate_hybrid(MM, n_bin, M, V, V_comp, bin_v, bin_r,
     $     bin_g, bin_n, dlnr, b1, s1, b2, s2, bin_change)

      integer MM           ! INPUT: physical dimension of V
      integer n_bin        ! INPUT: number of bins
      integer M(n_bin)     ! INPUT/OUTPUT: number particles per bin
      real*8 VH(MM)         ! INPUT/OUTPUT: particle volumes
      real*8 V_comp        ! INPUT: computational volume

      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins
      integer bin_n(n_bin) ! INPUT/OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor
      
      integer b1           ! INPUT: first particle (bin number)
      integer s1           ! INPUT: first particle (number in bin)
      integer b2           ! INPUT: second particle (bin number)
      integer s3           ! INPUT: second particle (number in bin)
      logical bin_change   ! OUTPUT: whether an empty bin filled,
                           !         or a filled bin became empty

      integer k1, k2, kn
      real*8 tmp

      bin_change = .false.

      ! remove s1 and s2 from bins
      call particle_in_bin(VH(s1), n_bin, bin_v, k1)
      call particle_in_bin(VH(s2), n_bin, bin_v, k2)
      bin_n(k1) = bin_n(k1) - 1
      bin_n(k2) = bin_n(k2) - 1
      bin_g(k1) = bin_g(k1) - VH(s1)
      bin_g(k2) = bin_g(k2) - VH(s2)

      if ((bin_n(k1) .lt. 0) .or. (bin_n(k2) .lt. 0)) then
         write(*,*)'ERROR: invalid bin_n'
         call exit(2)
      endif

      ! add particles
      VH(s1) = VH(s1) + VH(s2)    ! add particle 2 onto particle 1
      if (k2 .le. s_bin) then  ! s2 was a small particle
         VH(s2) = VH(MS)         ! move last small particle into empty slot
         VH(MS) = VH(M)          ! move last big particle into empty slot
         MS = MS - 1           ! decrease the number of small particles
         M = M - 1             ! decrease the number of total particles
      else                     ! s2 was a big particle
         VH(s2) = VH(M)          ! shift the last particle into empty slot
         M = M - 1             ! decrease the number of total particles
      endif

      ! add new particle to bins
      call particle_in_bin(VH(s1), n_bin, bin_v, kn)
      bin_n(kn) = bin_n(kn) + 1
      bin_g(kn) = bin_g(kn) + VH(s1)

      if ((bin_n(k1) .eq. 0) .or. (bin_n(k2) .eq. 0))
     &     bin_change = .true.
      if ((bin_n(kn) .eq. 1) .and. (kn .ne. k1) .and. (kn .ne. k2))
     &     bin_change = .true.

      ! fix up big/small ordering of s1
      if ((k1 .le. s_bin) .and. (kn .gt. s_bin)) then
         ! if s1 was small but is now big
         tmp = VH(s1)       ! swap s1 and the last small particle
         VH(s1) = VH(MS)
         VH(MS) = tmp
         MS = MS - 1       ! decrease the number of small particles
      endif
      ! we can't have the case where s1 was big but is now small

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine double_hybrid(MM, n_bin, M, V, V_comp, bin_v, bin_r,
     $     bin_g, bin_n, dlnr)

C     Double number of particles in a hybrid array.

      integer MM           ! INPUT: physical dimension of V
      integer n_bin        ! INPUT: number of bins
      integer M(n_bin)     ! INPUT/OUTPUT: number particles per bin
      real*8 VH(MM)         ! INPUT/OUTPUT: particle volumes
      real*8 V_comp        ! INPUT: computational volume

      integer n_bin        ! INPUT: number of bins
      integer s_bin        ! INPUT: split point in bin grid
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins
      integer bin_n(n_bin) ! INPUT/OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor

      integer i

      ! only double if we have enough space to do so
      if (M .gt. MM / 2) then
         write(*,*)'ERROR: double without enough space'
         call exit(2)
      endif
      
      ! double V and associated structures
      do i = (MS + 1),M        ! copy bigs to top of new array
         VH(i + M) = VH(i)
      enddo
      do i = M,(MS + 1),-1     ! copy bigs to second-top of new array
         VH(i + MS) = VH(i + M)
      enddo
      do i = 1,MS              ! copy smalls above smalls
         VH(i + MS) = VH(i)
      enddo
      M = 2 * M
      MS = 2 * MS
      V_comp = 2d0 * V_comp
      
      ! double bin structures
      do i = 1,n_bin
         bin_g(i) = bin_g(i) * 2d0
         bin_n(i) = bin_n(i) * 2
      enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine check_hybrid(MM, M, MS, V, n_bin, bin_v, bin_r)

C     Check that V has all particles in the correct bins.

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT/OUTPUT: logical dimension of V
      integer MS           ! INPUT/OUTPUT: split point in V
      real*8 VH(MM)         ! INPUT/OUTPUT: particle volumes

      integer n_bin        ! INPUT: number of bins
      integer s_bin        ! INPUT: split point in bin grid
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins

      integer i, k
      logical error

      error = .false.
      do i = 1,M
         call particle_in_bin(VH(i), n_bin, bin_v, k)
         if (((i .le. MS) .and. (k .gt. s_bin)) .or. ((i .gt. MS) .and.
     $        (k .le. s_bin))) then
            write(*,'(a10,a10,a12,a8,a8)') 'i', 'MS', 'VH(i)', 'k',
     $           's_bin'
            write(*,'(i10,i10,e12.5,i8,i8)') i, MS, VH(i), k, s_bin
            error = .true.
         endif
      enddo
      if (error) call exit(2)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
