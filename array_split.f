C     Functions to deal with a particle array V that is divided into two
C     groups, consisting of big particles and small particles. Big and
C     small particles are determined by a cutoff in the bin grid.
C
C     The bin grid 1,...,n_bin is divided into small bins (1,...,s_bin)
C     and big bins ((s_bin+1),...,n_bin). Similarly, the V array is
C     divided into small particles (1,...,MS) and big particles
C     ((MS+1),...,M).

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine est_k_max_split(n_bin, s_bin, bin_v, bin_n, kernel,
     $     k_max_small, k_max_big)

      integer n_bin         ! INPUT: number of bins
      integer s_bin         ! INPUT: split point in bin grid
      real*8 bin_v(n_bin)   ! INPUT: volume of particles in bins (m^3)
      integer bin_n(n_bin)  ! INPUT: number in each bin
      external kernel       ! INPUT: kernel function
      real*8 k_max_small    ! OUTPUT: maximum kernel for small-small
      real*8 k_max_big      ! OUTPUT: maximum kernel for big-big or big-small

      real*8 k
      integer i, j
      logical use_bin(n_bin)
! DEBUG
      integer min_use_bin, max_use_bin, max_i, max_j
! DEBUG

      ! use_bin starts as non-empty bins
      do i = 1,n_bin
         use_bin(i) = (bin_n(i) .gt. 0)
      enddo

      ! add all bins downstream of non-empty bins
      do i = 2,n_bin
         if (use_bin(i)) use_bin(i-1) = .true.
      enddo

      ! add all bins upstream of non-empty bins
      do i = (n_bin-1),1,-1
         if (use_bin(i)) use_bin(i+1) = .true.
      enddo

! DEBUG
      min_use_bin = n_bin
      max_use_bin = 1
      do i = 1,n_bin
         if (use_bin(i)) then
            if (i .lt. min_use_bin) min_use_bin = i
            if (i .gt. max_use_bin) max_use_bin = i
         endif
      enddo
c      write(*,*)'min_use_bin = ', min_use_bin
c      write(*,*)'max_use_bin = ', max_use_bin
! DEBUG

      ! compute k_max_small
      k_max_small = 0d0
      do i = 1,s_bin  ! FIXME: go to s_bin+1 to avoid underestimation
         if (use_bin(i)) then
            do j = 1,i
               if (use_bin(j)) then
                  call kernel(bin_v(i), bin_v(j), k)
                  if (k .gt. k_max_small) then
                     k_max_small = k
! DEBUG
                     max_i = i
                     max_j = j
! DEBUG
                  endif
               endif
            enddo
         endif
      enddo
! DEBUG
c      write(*,*)'k_max_small = ', k_max_small
c      write(*,*)'using i,j = ', max_i, max_j
! DEBUG

      ! compute k_max_big
      k_max_big = 0d0
      do i = (s_bin + 1),n_bin
         if (use_bin(i)) then
            do j = 1,i
               if (use_bin(j)) then
                  call kernel(bin_v(i), bin_v(j), k)
                  if (k .gt. k_max_big) then
                     k_max_big = k
! DEBUG
                     max_i = i
                     max_j = j
! DEBUG
                  endif
               endif
            enddo
         endif
      enddo
! DEBUG
c      write(*,*)'k_max_big = ', k_max_big
c      write(*,*)'using i,j = ', max_i, max_j
! DEBUG

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine maybe_coag_pair_small(MM, M, MS, V, V_comp, n_bin,
     $     s_bin, bin_v, bin_r, bin_g, bin_n, dlnr, del_t, n_samp,
     $     kernel, did_coag, bin_change)

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT/OUTPUT: logical dimension of V
      integer MS           ! INPUT/OUTPUT: split point in V
      real*8 V(MM)         ! INPUT/OUTPUT: particle volumes
      real*8 V_comp        ! INPUT: computational volume

      integer n_bin        ! INPUT: number of bins
      integer s_bin        ! INPUT: split point in bin grid
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins
      integer bin_n(n_bin) ! INPUT/OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor
      
      real*8 del_t         ! INPUT: timestep
      integer n_samp       ! INPUT: number of samples per timestep
      external kernel      ! INPUT: kernel function
      logical did_coag     ! OUTPUT: whether a coagulation occured
      logical bin_change   ! OUTPUT: whether bin structure changed

      integer s1, s2
      real*8 p, k
      real*8 M_small, n_possible ! use real*8 to avoid integer overflow

      call find_rand_pair(MS - 1, s1, s2) ! test particles s1, s2
      call kernel(V(s1), V(s2), k)

      M_small = MS
      n_possible = M_small * (M_small - 1d0) / 2d0

      p = k * 1d0/V_comp * del_t * n_possible / n_samp
! DEBUG
c      write(*,*)'s1, s2, k, p, p_scale = ', s1, s2, k, p, 1d0/V_comp
c     $     * del_t * n_possible / n_samp
! DEBUG

      bin_change = .false.
      if (dble(rand()) .lt. p) then
         call coagulate_split(MM, M, MS, V, V_comp,
     &        n_bin, s_bin, bin_v, bin_r, bin_g, bin_n, dlnr,
     &        s1, s2, bin_change)
         did_coag = .true.
      else
         did_coag = .false.
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine maybe_coag_pair_big(MM, M, MS, V, V_comp, n_bin,
     $     s_bin, bin_v, bin_r, bin_g, bin_n, dlnr, del_t, n_samp,
     $     kernel, did_coag, bin_change)

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT/OUTPUT: logical dimension of V
      integer MS           ! INPUT/OUTPUT: split point in V
      real*8 V(MM)         ! INPUT/OUTPUT: particle volumes
      real*8 V_comp        ! INPUT: computational volume

      integer n_bin        ! INPUT: number of bins
      integer s_bin        ! INPUT: split point in bin grid
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins
      integer bin_n(n_bin) ! INPUT/OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor
      
      real*8 del_t         ! INPUT: timestep
      integer n_samp       ! INPUT: number of samples per timestep
      external kernel      ! INPUT: kernel function
      logical did_coag     ! OUTPUT: whether a coagulation occured
      logical bin_change   ! OUTPUT: whether bin structure changed

      integer s1, s2
      real*8 p, k
      real*8 M_small, M_big, n_possible ! use real*8 to avoid integer overflow

      call find_rand_pair_big(M, MS, s1, s2) ! test particles s1, s2
      call kernel(V(s1), V(s2), k)

      M_small = MS
      M_big = M - MS
      n_possible = M_big * (M_big - 1d0) / 2d0 + M_big * M_small

      p = k * 1d0/V_comp * del_t * n_possible / n_samp

      bin_change = .false.
      if (dble(rand()) .lt. p) then
         call coagulate_split(MM, M, MS, V, V_comp,
     &        n_bin, s_bin, bin_v, bin_r, bin_g, bin_n, dlnr,
     &        s1, s2, bin_change)
         did_coag = .true.
      else
         did_coag = .false.
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     
      subroutine find_rand_pair_big(M, MS, s1, s2)
      
      integer M       ! INPUT: number of particles
      integer MS      ! INPUT: split point in V
      integer s1, s2  ! OUTPUT: s1 and s2 are not equal, random
                      !         particles with MS < s1 <= M
                      !                    and  1 < s2 <= M

 100  s1 = int(rand() * (M - MS)) + 1 + MS
      if ((s1 .lt. MS + 1) .or. (s1 .gt. M)) goto 100
 101  s2 = int(rand() * M) + 1
      if ((s2 .lt. 1) .or. (s2 .gt. M)) goto 101

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine coagulate_split(MM, M, MS, V, V_comp, n_bin, s_bin,
     $     bin_v, bin_r, bin_g, bin_n, dlnr, s1, s2, bin_change)

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT/OUTPUT: logical dimension of V
      integer MS           ! INPUT/OUTPUT: split point in V
      real*8 V(MM)         ! INPUT/OUTPUT: particle volumes  (m^3)
      real*8 V_comp        ! INPUT: computational volume   (m^3)

      integer n_bin        ! INPUT: number of bins
      integer s_bin        ! INPUT; split point in bin grid
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins (m^3)
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins (m)
      real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins 
      integer bin_n(n_bin) ! INPUT/OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor

      integer s1           ! INPUT: first particle to coagulate
      integer s2           ! INPUT: second particle to coagulate
      logical bin_change   ! OUTPUT: whether an empty bin filled,
                           !         or a filled bin became empty

      integer k1, k2, kn
      real*8 tmp

      bin_change = .false.

      ! remove s1 and s2 from bins
      call particle_in_bin(V(s1), n_bin, bin_v, k1)
      call particle_in_bin(V(s2), n_bin, bin_v, k2)
      bin_n(k1) = bin_n(k1) - 1
      bin_n(k2) = bin_n(k2) - 1
      bin_g(k1) = bin_g(k1) - V(s1)
      bin_g(k2) = bin_g(k2) - V(s2)

      if ((bin_n(k1) .lt. 0) .or. (bin_n(k2) .lt. 0)) then
         write(*,*)'ERROR: invalid bin_n'
         call exit(2)
      endif

      ! add particles
      V(s1) = V(s1) + V(s2)    ! add particle 2 onto particle 1
      if (k2 .le. s_bin) then  ! s2 was a small particle
         V(s2) = V(MS)         ! move last small particle into empty slot
         V(MS) = V(M)          ! move last big particle into empty slot
         MS = MS - 1           ! decrease the number of small particles
         M = M - 1             ! decrease the number of total particles
      else                     ! s2 was a big particle
         V(s2) = V(M)          ! shift the last particle into empty slot
         M = M - 1             ! decrease the number of total particles
      endif

      ! add new particle to bins
      call particle_in_bin(V(s1), n_bin, bin_v, kn)
      bin_n(kn) = bin_n(kn) + 1
      bin_g(kn) = bin_g(kn) + V(s1)

      if ((bin_n(k1) .eq. 0) .or. (bin_n(k2) .eq. 0))
     &     bin_change = .true.
      if ((bin_n(kn) .eq. 1) .and. (kn .ne. k1) .and. (kn .ne. k2))
     &     bin_change = .true.

      ! fix up big/small ordering of s1
      if ((k1 .le. s_bin) .and. (kn .gt. s_bin)) then
         ! if s1 was small but is now big
         tmp = V(s1)       ! swap s1 and the last small particle
         V(s1) = V(MS)
         V(MS) = tmp
         MS = MS - 1       ! decrease the number of small particles
      endif
      ! we can't have the case where s1 was big but is now small

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine double_split(MM, M, MS, V, V_comp, n_bin, s_bin, bin_v,
     $     bin_r, bin_g, bin_n, dlnr)

C     Double number of particles in a split array.

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT/OUTPUT: logical dimension of V
      integer MS           ! INPUT/OUTPUT: split point in V
      real*8 V(MM)         ! INPUT/OUTPUT: particle volumes
      real*8 V_comp        ! INPUT/OUTPUT: computational volume

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
         V(i + M) = V(i)
      enddo
      do i = M,(MS + 1),-1     ! copy bigs to second-top of new array
         V(i + MS) = V(i + M)
      enddo
      do i = 1,MS              ! copy smalls above smalls
         V(i + MS) = V(i)
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

      subroutine find_split_bin(n_bin, bin_v, r_split, s_bin)

C     Find the bin to split at, given an approximate radius r_split to
C     split at. The split boundary will be the smallest boundary greater
C     than r_split.

      integer n_bin        ! INPUT: number of bins
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins (m^3)
      real*8 r_split       ! INPUT: desired split radius (m)
      integer s_bin        ! OUTPUT: split point in bin grid

      real*8 v_split, v_mid

      real*8 pi
      parameter (pi = 3.14159265358979323846d0)

      v_split = 4d0/3d0 * pi * r_split**3

      if (v_split .lt. bin_v(1)) then
         write(*,*)'ERROR: fell off bottom in find_split_bin()'
         call exit(2)
      endif

C     for (s_bin = 1; s_bin < n_bin; s_bin++)
C        if (test(s_bin + 1))
C           break;
C     if (s_bin == n_bin)
C        error();
      s_bin = 1
 100  if (.not. (s_bin .lt. n_bin)) goto 101
      v_mid = (bin_v(s_bin) + bin_v(s_bin + 1)) / 2d0
      if (v_split .lt. v_mid) goto 101
      s_bin = s_bin + 1
      goto 100
 101  if (s_bin .eq. n_bin) then
         write(*,*)'ERROR: fell off top in find_split_bin()'
         call exit(2)
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine find_split_particle(MM, M, V, n_bin, s_bin, bin_v, MS)

C     Determine the point MS at which the V array is split. MS will be
C     the last small entry. V must be sortest smallest to largest on
C     entry.

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT/OUTPUT: logical dimension of V
      real*8 V(MM)         ! INPUT/OUTPUT: particle volumes

      integer n_bin        ! INPUT: number of bins
      integer s_bin        ! INPUT: split point in bin grid
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins

      integer MS           ! OUTPUT: split point in V

      integer k

      MS = 0
 200  call particle_in_bin(V(MS + 1), n_bin, bin_v, k)
      if (k .gt. s_bin) goto 201    ! MS + 1 is big
      MS = MS + 1
      if (MS .ge. M) goto 201       ! no big particles
      goto 200
 201  continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine check_split(MM, M, MS, V, n_bin, s_bin, bin_v, bin_r)

C     Check that V is split correctly.

      integer MM           ! INPUT: physical dimension of V
      integer M            ! INPUT/OUTPUT: logical dimension of V
      integer MS           ! INPUT/OUTPUT: split point in V
      real*8 V(MM)         ! INPUT/OUTPUT: particle volumes

      integer n_bin        ! INPUT: number of bins
      integer s_bin        ! INPUT: split point in bin grid
      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins

      integer i, k
      logical error

      error = .false.
      do i = 1,M
         call particle_in_bin(V(i), n_bin, bin_v, k)
         if (((i .le. MS) .and. (k .gt. s_bin)) .or. ((i .gt. MS) .and.
     $        (k .le. s_bin))) then
            write(*,'(a10,a10,a12,a8,a8)') 'i', 'MS', 'V(i)', 'k',
     $           's_bin'
            write(*,'(i10,i10,e12.5,i8,i8)') i, MS, V(i), k, s_bin
            error = .true.
         endif
      enddo
      if (error) call exit(2)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
