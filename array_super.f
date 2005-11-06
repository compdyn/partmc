C     Functions to deal with a particle array VH that is stored by bin
C     superparticle size. The allowable superparticle factors are 1, a,
C     a^2, a^3, etc. This is the number of physical particles
C     represented by the superparticle. In the code the base a is the
C     variable sup_base.
C
C     MS(n_bin, n_fact): MS(b, f) is the number of computational
C     particles in bin b with factor exponent f.
C
C     VS(n_bin, n_fact, TDV): VS(b, f, i) is the volume of the physical
C     particles that make up the i-th superparticle in bin b with factor
C     exponent f. This superparticle represents a^{f - 1} physical
C     particles of this volume.
C
C     We use the terminology "size" to refer to the radius/volume/mass
C     of a physical particle, and "factor" to refer to the number of
C     physical particles represented by a single computational
C     superparticle.
C
C     single particles v0, v1 with rate k
C
C     super-N particle v0, single particle v1: really N events, so
C     expect rate kN. One super-N/single collision will produce N
C     collisions, so use rate k.
C
C     super-N particles v0, v1: really N^2 events, so expect rate
C     kN^2. One super-N/super-N collision will produce N collisions, so
C     use rate kN.
C
C     2 <--> 10      20 per sec  ==>  k = 0.05
C     1 <--> 5        5 per sec  ==>  k = 0.2
C
C     100 <--> 100       k = 0.001  ==>  10 collisions per sec
C     10x10 <--> 10x10   k = 0.01   ==>  1x10 collisions per sec
C
C     super-N particle v0 colliding with super-M particle v1 with N < M:
C     really NM events, so expect rate kNM. One super-N/super-M
C     collision will produce M physical collisions, so use rate kN.
C
C     100 <--> 1000      k = 0.001  ==>  100 collisions per sec
C     1x100 <--> 1x1000  k = 0.1    ==>  0.1x1000 collisions per sec
C
C     There are two possibilites for how to do the collisions:
C
C     1. Do everything per-bin, per-fact, which is just like hybrid, but
C     with more compartmentalization.
C
C     2. Do things just per-bin, then randomly collide either single
C     particles or superparticles, and account for this fact.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine init_to_super(n_bin, n_fact, TDV, sup_base, MS, VS,
     $     n_ini, min_fill)

C     Convert an initial number distribution into a superparticle
C     array. We try to ensure that we have at least min_fill
C     computational particles per bin, but within that restriction we
C     use the biggest superparticles we can. A simple greedy algorithm
C     is used, which is slightly suboptimal.

      integer n_bin               ! INPUT: number of bins
      integer n_fact              ! INPUT: number of allowable factors
      integer TDV                 ! INPUT: trailing dimension of VS
      integer sup_base            ! INPUT: factor base of a superparticle
      integer MS(n_bin,n_fact)    ! OUTPUT: number of superparticles
      real*8 VS(n_bin,n_fact,TDV) ! OUTPUT: volume of physical particles
      real*8 bin_v(n_bin)         ! INPUT: volume of particles in bins
      integer bin_n(n_bin)        ! INPUT: desired number of particles in bins
      integer min_fill            ! INPUT: desired minimum number of
                                  !        computational particles in each bin

      integer b, f, i

C DEBUG
      do b = 1,n_bin
         do f = 1,n_fact
            MS(b, f) = 99999999
            do i = 1,TDV
               VS(b, f, i) = 9e200
            enddo
         enddo
      enddo
C DEBUG

      ! make MS
      do b = 1,n_bin
         ! find the maximum factor we can use which will still result
         ! in min_fill computational particles of that maximum factor
         max_f = 0
         done = .false.
         do while (.not. done)
            max_f = max_f + 1
            if (max_f .ge. n_fact) done = .true.
            factor = sup_base**((max_f + 1) - 1)
            if (bin_n(b) / factor .lt. min_fill) done = .true.
         enddo
         ! work our way down the factor sizes, starting from the maximum
         ! factor, with a greedy algorithm that makes as many particles
         ! of each factor as possible. We finish with factor = 1, so we
         ! will always exactly have bin_n(b) physical particles represented.
         n_left = bin_n(b)
         do f = max_f,1,-1
            factor = sup_base**(f - 1)
            MS(b,f) = n_left / factor
            n_left = n_left - MS(b,f) * factor
         enddo
      enddo

      ! make VS
      do b = 1,n_bin
         do f = 1,n_fact
            do i = 1,MS(b,f)
               VS(b,f,i) = bin_v(b)
            enddo
         enddo
      enddo
         
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine maybe_coag_pair_super(M, n_bin, n_fact, TDV, sup_base,
     $     MS, VS, min_fill, V_comp, bin_v, bin_r, bin_g, bin_n, dlnr,
     $     b1, b2, f1, f2, del_t, k_max, kernel, did_coag, bin_change)

C     Choose a random pair for potential coagulation and test its
C     probability of coagulation. If it happens, do the coagulation and
C     update all structures. The probability of a coagulation will be
C     taken as (kernel / k_max).

      integer M                   ! INPUT/OUTPUT: number of particles
      integer n_bin               ! INPUT: number of bins
      integer n_fact              ! INPUT: number of allowable factors
      integer TDV                 ! INPUT: trailing dimension of VS
      integer sup_base            ! INPUT: factor base of a superparticle
      integer MS(n_bin,n_fact)    ! INPUT/OUTPUT: number of superparticles
      real*8 VS(n_bin,n_fact,TDV) ! INPUT/OUTPUT: volume of physical particles
      integer min_fill            ! INPUT: minimum comp. part. per bin
      real*8 V_comp               ! INPUT: computational volume
      real*8 bin_v(n_bin)         ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)         ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)         ! INPUT/OUTPUT: mass in bins
      integer bin_n(n_bin)        ! INPUT/OUTPUT: number in bins
      real*8 dlnr                 ! INPUT: bin scale factor

      integer b1                  ! INPUT: bin of first particle
      integer b2                  ! INPUT: bin of second particle
      integer f1                  ! INPUT: factor step of first particle
      integer f2                  ! INPUT: factor step of second particle
      real*8 del_t                ! INPUT: timestep
      real*8 k_max                ! INPUT: k_max scale factor
      external kernel             ! INPUT: kernel function
      logical did_coag            ! OUTPUT: whether a coagulation occured
      logical bin_change          ! OUTPUT: whether bin structure changed

      integer s1, s2
      real*8 p, k

      bin_change = .false.
      did_coag = .false.

      if ((MH(b1) .le. 0) .or. (MH(b2) .le. 0)) then
         return
      endif

      call find_rand_pair_hybrid(n_bin, MH, b1, b2, s1, s2)

      call kernel(VH(b1, s1), VH(b2, s2), k)
      p = k / k_max

      if (dble(rand()) .lt. p) then
         call coagulate_hybrid(M, n_bin, TDV, MH, VH, V_comp, bin_v,
     $        bin_r,bin_g, bin_n, dlnr, b1, s1, b2, s2, bin_change)
         did_coag = .true.
      endif

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     
      subroutine find_rand_pair_super(n_bin, n_fact, MS, b1, b2, f1, f2,
     $     s1, s2)

C     Find a random pair s1, s2 so that we have a random non-equal pair
C     of particles (b1, f1, s1) and (b2, f2, s2).

      integer n_bin            ! INPUT: number of bins
      integer n_fact           ! INPUT: number of allowable factors
      integer MS(n_bin,n_fact) ! INPUT/OUTPUT: number of superparticles
      integer b1               ! INPUT: bin of first particle
      integer b2               ! INPUT: bin of second particle
      integer f1               ! INPUT: factor step of first particle
      integer f2               ! INPUT: factor step of second particle
      integer s1               ! OUTPUT: first random particle
                               !         1 <= s1 <= M(b1)
      integer s2               ! OUTPUT: second random particle
                               !         1 <= s2 <= M(b2)
                               !         (b1,s1) != (b2,s2)

 100  s1 = int(rand() * MS(b1,f1)) + 1
      if ((s1 .lt. 1) .or. (s1 .gt. MS(b1,f1))) goto 100
 101  s2 = int(rand() * MS(b2,f2)) + 1
      if ((s2 .lt. 1) .or. (s2 .gt. MS(b2,f2))) goto 101
      if ((b1 .eq. b2) .and. (s1 .eq. s2)) goto 101

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine coagulate_super(M, n_bin, TDV, MH, VH, V_comp, bin_v,
     $     bin_r, bin_g, bin_n, dlnr, b1, s1, b2, s2, bin_change)

C     Join together particles (b1, f1, s1) and (b2, f2, s2), updating all
C     particle and bin structures to reflect the change. bin_change is
C     true if the used bin set changed due to the coagulation (i.e. an
C     empty bin filled or a filled bin became empty).

      integer M                   ! INPUT/OUTPUT: number of particles
      integer n_bin               ! INPUT: number of bins
      integer n_fact              ! INPUT: number of allowable factors
      integer TDV                 ! INPUT: trailing dimension of VS
      integer sup_base            ! INPUT: factor base of a superparticle
      integer MS(n_bin,n_fact)    ! INPUT/OUTPUT: number of superparticles
      real*8 VS(n_bin,n_fact,TDV) ! INPUT/OUTPUT: volume of physical particles
      integer min_fill            ! INPUT: minimum comp. part. per bin
      real*8 V_comp               ! INPUT: computational volume
      real*8 bin_v(n_bin)         ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)         ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)         ! INPUT/OUTPUT: mass in bins
      integer bin_n(n_bin)        ! INPUT/OUTPUT: number in bins
      real*8 dlnr                 ! INPUT: bin scale factor

      integer b1                  ! INPUT: bin of first particle
      integer b2                  ! INPUT: bin of second particle
      integer f1                  ! INPUT: factor step of first particle
      integer f2                  ! INPUT: factor step of second particle
      real*8 del_t                ! INPUT: timestep
      real*8 k_max                ! INPUT: k_max scale factor
      external kernel             ! INPUT: kernel function
      logical did_coag            ! OUTPUT: whether a coagulation occured
      logical bin_change          ! OUTPUT: whether bin structure changed

      integer M            ! INPUT/OUTPUT: number of particles
      integer n_bin        ! INPUT: number of bins
      integer TDV          ! INPUT: trailing dimension of VH
      integer MH(n_bin)    ! INPUT/OUTPUT: number of particles per bin
      real*8 VH(n_bin,TDV) ! INPUT/OUTPUT: particle volumes
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

      subroutine double_hybrid(M, n_bin, TDV, MH, VH, V_comp, bin_v,
     $     bin_r, bin_g, bin_n, dlnr)

C     Double number of particles in a hybrid array.

      integer M            ! INPUT/OUTPUT: number of particles
      integer n_bin        ! INPUT: number of bins
      integer TDV          ! INPUT: trailing dimension of VH
      integer MH(n_bin)    ! INPUT/OUTPUT: number of particles per bin
      real*8 VH(n_bin,TDV) ! INPUT/OUTPUT: particle volumes
      real*8 V_comp        ! INPUT/OUTPUT: computational volume

      real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins
      integer bin_n(n_bin) ! INPUT/OUTPUT: number in bins
      real*8 dlnr          ! INPUT: bin scale factor

      integer i, k

      ! double VH and associated structures
      do k = 1,n_bin
         ! only double if we have enough space to do so
         if (2 * MH(k) .gt. TDV) then
            write(*,*)'ERROR: double without enough space'
            call exit(2)
         endif
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

      subroutine check_hybrid(M, n_bin, TDV, MH, VH, bin_v, bin_r)

C     Check that V has all particles in the correct bins.

      integer M            ! INPUT: number of particles
      integer n_bin        ! INPUT: number of bins
      integer TDV          ! INPUT: trailing dimension of VH
      integer MH(n_bin)    ! INPUT: number of particles per bin
      real*8 VH(n_bin,TDV) ! INPUT: particle volumes

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
