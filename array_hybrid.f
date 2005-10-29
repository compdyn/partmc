C     Functions to deal with a particle array VH that is stored by
C     bin. VH(i_bin,i) is the i-th particle in the i_bin-th bin.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine array_to_hybrid(MM, M, V, n_bin, bin_v, MH, VH)

C     Convert a standard linear particle array V to a hybrid particle
C     array VH stored by bins.

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

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine maybe_coag_pair_hybrid(MM, M, n_bin, MH, VH, V_comp,
     $     bin_v, bin_r, bin_g, bin_n, dlnr, del_t, k_max, kernel,
     $     did_coag, bin_change)

C     Choose a random pair for potential coagulation and test its
C     probability of coagulation. If it happens, do the coagulation and
C     update all structures. The probability of a coagulation will be
C     taken as (kernel / k_max).

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

      bin_change = .false.
      did_coag = .false.

      if ((MH(b1) .le. 0) .or. (MH(b2) .le. 0)) then
         return
      endif

      call find_rand_pair_hybrid(n_bin, MH, b1, s1, b2, s2)

      call kernel(VH(b1, s1), VH(b2, s2), k)
      p = k / k_max

      if (dble(rand()) .lt. p) then
         call coagulate_hybrid(MM, M, n_bin, MH, VH, V_comp, bin_v,
     $        bin_r,bin_g, bin_n, dlnr, b1, s1, b2, s2, bin_change)
         did_coag = .true.
      endif

      end

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

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine coagulate_hybrid(MM, M, n_bin, MH, VH, V_comp, bin_v,
     $     bin_r, bin_g, bin_n, dlnr, b1, s1, b2, s2, bin_change)

C     Join together particles (b1, s1) and (b2, s2), updating all
C     particle and bin structures to reflect the change. bin_change is
C     true if the used bin set changed due to the coagulation (i.e. an
C     empty bin filled or a filled bin became empty).

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
      
      integer b1           ! INPUT: first particle (bin number)
      integer s1           ! INPUT: first particle (number in bin)
      integer b2           ! INPUT: second particle (bin number)
      integer s2           ! INPUT: second particle (number in bin)
      logical bin_change   ! OUTPUT: whether an empty bin filled,
                           !         or a filled bin became empty

      integer bn
      real*8 new_v

      bin_change = .false.

      ! remove s1 and s2 from bins
      bin_n(b1) = bin_n(b1) - 1
      bin_n(b2) = bin_n(b2) - 1
      bin_g(b1) = bin_g(b1) - VH(b1, s1)
      bin_g(b2) = bin_g(b2) - VH(b2, s2)
      if ((bin_n(b1) .lt. 0) .or. (bin_n(b2) .lt. 0)) then
         write(*,*)'ERROR: invalid bin_n'
         call exit(2)
      endif

      ! do coagulation in MH, VH arrays
      new_v = VH(b1, s1) + VH(b2, s2)   ! add particle volumes
      call particle_in_bin(new_v, n_bin, bin_v, bn)  ! find new bin
      VH(b1, s1) = VH(b1, MH(b1))  ! shift last particle into empty slot
      MH(b1) = MH(b1) - 1          ! decrease length of array
      VH(b2, s2) = VH(b2, MH(b2))  ! same for second particle
      MH(b2) = MH(b2) - 1
      if ((MH(b1) .lt. 0) .or. (MH(b2) .lt. 0)) then
         write(*,*)'ERROR: invalid MH'
         call exit(2)
      endif
      MH(bn) = MH(bn) + 1          ! increase the length of array
      VH(bn, MH(bn)) = new_v       ! add the new particle at the end
      M = M - 1                    ! decrease the total number of particles

      ! add new particle to bins
      bin_n(bn) = bin_n(bn) + 1
      bin_g(bn) = bin_g(bn) + VH(bn, MH(bn))

      ! did we empty a bin?
      if ((bin_n(b1) .eq. 0) .or. (bin_n(b2) .eq. 0))
     &     bin_change = .true.
      ! did we start a new bin?
      if ((bin_n(bn) .eq. 1) .and. (bn .ne. b1) .and. (bn .ne. b2))
     &     bin_change = .true.

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine double_hybrid(MM, M, n_bin, MH, VH, V_comp, bin_v,
     $     bin_r, bin_g, bin_n, dlnr)

C     Double number of particles in a hybrid array.

      integer MM           ! INPUT: trailing dimension of V
      integer M            ! INPUT/OUTPUT: number of particles
      integer n_bin        ! INPUT: number of bins
      integer MH(n_bin)    ! INPUT/OUTPUT: number of particles per bin
      real*8 VH(n_bin,MM)  ! INPUT/OUTPUT: particle volumes
      real*8 V_comp        ! INPUT/OUTPUT: computational volume

      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins
      integer bin_n(n_bin) ! INPUT/OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor

      integer i, k

      ! only double if we have enough space to do so
      if (M .gt. MM / 2) then
         write(*,*)'ERROR: double without enough space'
         call exit(2)
      endif
      
      ! double VH and associated structures
      do k = 1,n_bin
         do i = 1,MH(k)
            VH(k, i + MH(k)) = VH(k, i)
         enddo
         MH(k) = 2 * MH(k)
      enddo
      M = 2 * M
      V_comp = 2d0 * V_comp

      ! double bin structures
      do i = 1,n_bin
         bin_g(i) = bin_g(i) * 2d0
         bin_n(i) = bin_n(i) * 2
      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine check_hybrid(MM, M, n_bin, MH, VH, bin_v, bin_r)

C     Check that V has all particles in the correct bins.

      integer MM           ! INPUT: trailing dimension of V
      integer M            ! INPUT: number of particles
      integer n_bin        ! INPUT: number of bins
      integer MH(n_bin)    ! INPUT: number of particles per bin
      real*8 VH(n_bin,MM)  ! INPUT: particle volumes

      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins

      integer i, k, k_check, M_check
      logical error

      error = .false.
      M_check = 0
      do k = 1,n_bin
         M_check = M_check + MH(k)
         do i = 1,MH(k)
            call particle_in_bin(VH(k, i), n_bin, bin_v, k_check)
            if (k .ne. k_check) then
               write(*,'(a10,a10,a12,a10)') 'k', 'i', 'VH(k, i)',
     $              'k_check'
               write(*,'(i10,i10,e12.5,i10)') k, i, VH(k, i), k_check
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

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
