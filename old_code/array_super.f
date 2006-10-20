C     Functions to deal with a particle array VH that is stored by bin
C     superparticle size. The allowable superparticle factors are 1, a,
C     a^2, a^3, etc. This is the number of physical particles
C     represented by the superparticle. In the code the base a is the
C     variable fac_base.
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
C     super-N particle v0 colliding with super-M particle v1 with N < M:
C     really NM events, so expect rate kNM. One super-N/super-M
C     collision will produce M physical collisions, so use rate kN.
C
C     n_phys_1 = n_comp_1 * factor_1
C     n_phys_2 = n_comp_2 * factor_2
C     n_phys_coll = k_phys * n_phys_1 * n_phys_2
C
C     n_comp_coll = k_comp * n_comp_1 * n_comp_2
C     n_phys_coll_per_comp_coll = max(factor_1, factor_2)
C     n_phys_coll = n_phys_coll_per_comp_coll * n_comp_coll
C
C     k_phys * n_phys_1 * n_phys_2
C              = n_phys_coll_per_comp_coll * n_comp_coll
C              = max(factor_1, factor_2) * k_comp * n_comp_1 * n_comp_2
C              = max(factor_1, factor_2)
C                        * k_comp * n_phys_1 / factor_1 * n_phys_2 / factor_2
C
C     k_phys = max(factor_1, factor_2) * k_comp / factor_1 / factor_2
C     k_comp = k_phys * factor_1 * factor_2 / max(factor_1, factor_2)
C            = k_phys * min(factor_1, factor_2)
C
C     2 <--> 10          k = 0.05   ==>  20 collisions per sec
C     1 <--> 5           k = 0.2    ==>  5 collisions per sec
C
C     100 <--> 100       k = 0.001  ==>  10 collisions per sec
C     10x10 <--> 10x10   k = 0.01   ==>  1x10 collisions per sec
C
C     100x1 <--> 1000x1  k = 0.001  ==>  0.001*100*1000 = 100 collisions per sec
C     1x100 <--> 1x1000  k = 0.1    ==>  0.1*1*1 = 0.1x1000 collisions per sec
C     10x10 <--> 10x100  k = 0.01   ==>  1x100 collisions per sec
C     10x10 <--> 100x10  k = 0.01   ==>  10x10 collisions per sec
C
C     Need to handle self-collision case properly.
C
C     There are two possibilites for how to do the collisions:
C
C     1. Do everything per-bin, per-fact, which is just like hybrid, but
C     with more compartmentalization.
C
C     2. Do things just per-bin, then randomly collide either single
C     particles or superparticles, and account for this fact.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine init_to_super(n_bin, n_fact, TDV, fac_base, MS, VS,
     $     bin_v, bin_n, min_fill)

C     Convert an initial number distribution into a superparticle
C     array. We try to ensure that we have at least min_fill
C     computational particles per bin, but within that restriction we
C     use the biggest superparticles we can. A simple greedy algorithm
C     is used, which is slightly suboptimal.

      integer n_bin               ! INPUT: number of bins
      integer n_fact              ! INPUT: number of allowable factors
      integer TDV                 ! INPUT: trailing dimension of VS
      integer fac_base            ! INPUT: factor base of a superparticle
      integer MS(n_bin,n_fact)    ! OUTPUT: number of superparticles
      real*8 VS(n_bin,n_fact,TDV) ! OUTPUT: volume of physical particles
      real*8 bin_v(n_bin)         ! INPUT: volume of particles in bins
      integer bin_n(n_bin)        ! INPUT: desired number of particles in bins
      integer min_fill            ! INPUT: desired minimum number of
                                  !        computational particles in each bin

      integer b, f, i, max_f, factor, n_left
      logical done

C DEBUG
      do b = 1,n_bin
         do f = 1,n_fact
            MS(b, f) = 99999999
            do i = 1,TDV
               VS(b, f, i) = 9d200
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
            factor = fac_base**((max_f + 1) - 1)
            if (bin_n(b) / factor .lt. min_fill) done = .true.
         enddo
         ! work our way down the factor sizes, starting from the maximum
         ! factor, with a greedy algorithm that makes as many particles
         ! of each factor as possible. We finish with factor = 1, so we
         ! will always exactly have bin_n(b) physical particles represented.
         do f = n_fact,(max_f + 1),-1
            MS(b,f) = 0
         enddo
         n_left = bin_n(b)
         do f = max_f,1,-1
            factor = fac_base**(f - 1)
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

      subroutine maybe_coag_pair_super(M, n_bin, n_fact, TDV, fac_base,
     $     MS, VS, min_fill, V_comp, bin_v, bin_r, bin_g, bin_n, dlnr,
     $     b1, b2, f1, f2, del_t, k_max, kernel, did_coag, bin_change)

C     Choose a random pair for potential coagulation and test its
C     probability of coagulation. If it happens, do the coagulation and
C     update all structures. The probability of a coagulation will be
C     taken as (factor * kernel / k_max).

      integer M                   ! INPUT/OUTPUT: number of particles
      integer n_bin               ! INPUT: number of bins
      integer n_fact              ! INPUT: number of allowable factors
      integer TDV                 ! INPUT: trailing dimension of VS
      integer fac_base            ! INPUT: factor base of a superparticle
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

      integer s1, s2, min_f, factor
      real*8 p, k

      bin_change = .false.
      did_coag = .false.

      if ((MS(b1, f1) .le. 0) .or. (MS(b2, f2) .le. 0)) then
         return
      endif

      call find_rand_pair_super(n_bin, n_fact, MS, b1, b2, f1, f2,
     $     s1, s2)

      call kernel(VS(b1, f1, s1), VS(b2, f2, s2), k)
      min_f = min(f1, f2)
      factor = fac_base**(min_f - 1)
      p = factor * k / k_max

      if (dble(rand()) .lt. p) then
         call coagulate_super(M, n_bin, n_fact, TDV, fac_base, MS, VS,
     $        min_fill, V_comp, bin_v, bin_r, bin_g, bin_n, dlnr, b1, f1
     $        , s1, b2, f2, s2, bin_change)
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
      integer MS(n_bin,n_fact) ! INPUT: number of superparticles
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

      subroutine coagulate_super(M, n_bin, n_fact, TDV, fac_base, MS, VS
     $     , min_fill, V_comp, bin_v, bin_r, bin_g, bin_n, dlnr, b1, f1,
     $     s1, b2, f2, s2, bin_change)

C     Join together particles (b1, f1, s1) and (b2, f2, s2), updating all
C     particle and bin structures to reflect the change. bin_change is
C     true if the used bin set changed due to the coagulation (i.e. an
C     empty bin filled or a filled bin became empty).

      integer M                   ! INPUT/OUTPUT: number of particles
      integer n_bin               ! INPUT: number of bins
      integer n_fact              ! INPUT: number of allowable factors
      integer TDV                 ! INPUT: trailing dimension of VS
      integer fac_base            ! INPUT: factor base of a superparticle
      integer MS(n_bin,n_fact)    ! INPUT/OUTPUT: number of superparticles
      real*8 VS(n_bin,n_fact,TDV) ! INPUT/OUTPUT: volume of physical particles
      integer min_fill            ! INPUT: minimum comp. part. per bin
      real*8 V_comp               ! INPUT: computational volume
      real*8 bin_v(n_bin)         ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)         ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)         ! INPUT/OUTPUT: mass in bins
      integer bin_n(n_bin)        ! INPUT/OUTPUT: number in bins
      real*8 dlnr                 ! INPUT: bin scale factor

      integer b1                  ! INPUT: first particle (bin number)
      integer f1                  ! INPUT: first particle (factor step)
      integer s1                  ! INPUT: first particle (particle number)
      integer b2                  ! INPUT: second particle (bin number)
      integer f2                  ! INPUT: second particle (factor step)
      integer s2                  ! INPUT: second particle (particle number)

      logical bin_change          ! OUTPUT: whether bin structure changed

      integer bn, fn, n1, n2, nn, f, f_red, n_comp, n_comp_add
      integer over1, over2, i
      real*8 new_v

      bin_change = .false.

      ! remove s1 and s2 from bins
      n1 = fac_base**(f1 - 1)
      n2 = fac_base**(f2 - 1)
      bin_n(b1) = bin_n(b1) - n1
      bin_n(b2) = bin_n(b2) - n2
      bin_g(b1) = bin_g(b1) - n1 * VS(b1, f1, s1)
      bin_g(b2) = bin_g(b2) - n2 * VS(b2, f2, s2)
      M = M - n1 - n2      ! decrease the total number of particles
      if ((bin_n(b1) .lt. 0) .or. (bin_n(b2) .lt. 0)) then
         write(*,*)'ERROR: invalid bin_n'
         call exit(2)
      endif

      ! find bin of new particle
      fn = min(f1, f2) ! new super-particle has factor of smaller of inputs
      over1 = fac_base**(f1 - fn) ! multiple of the new factor
      over2 = fac_base**(f2 - fn)
      ! add particle volumes
      new_v = over1 * VS(b1, f1, s1) + over2 * VS(b2, f2, s2)
      call particle_in_bin(new_v, n_bin, bin_v, bn)  ! find new bin
      
      ! remove particles from MS, VS arrays
      VS(b1, f1, s1) = VS(b1, f1, MS(b1, f1)) ! last particle into empty slot
      MS(b1, f1) = MS(b1, f1) - 1             ! decrease length of array
      VS(b2, f2, s2) = VS(b2, f2, MS(b2, f2)) ! same for second particle
      MS(b2, f2) = MS(b2, f2) - 1
      if ((MS(b1, f1) .lt. 0) .or. (MS(b2, f2) .lt. 0)) then
         write(*,*)'ERROR: invalid MS'
         call exit(2)
      endif

      ! did we empty a bin?
      if ((bin_n(b1) .eq. 0) .or. (bin_n(b2) .eq. 0))
     &     bin_change = .true.
      ! did we start a new bin?
      if ((bin_n(bn) .eq. 0) .and. (bn .ne. b1) .and. (bn .ne. b2))
     &     bin_change = .true.

      ! add new particle to MS, VS arrays
      ! if we don't have min_fill computational particles in the new bin
      ! then split the one we are adding
      f_red = 0 ! factor to split by
 200  n_comp = 0 ! number of computational particles
      do f = 1,n_fact
         n_comp = n_comp + MS(bn, f)
      enddo
      n_comp = n_comp + fac_base**f_red ! include the new ones we will add
      if ((n_comp .lt. min_fill) .and. (fn - f_red .gt. 1)) then
         f_red = f_red + 1
         goto 200 ! try again with one more level of splitting
      endif
      n_comp_add = fac_base**f_red  ! number of computational particles to add
      fn = fn - f_red               ! modified factor of new particles
      do i = 1,n_comp_add
         VS(bn, fn, MS(bn, fn) + i) = new_v  ! add new particles at the end
      enddo
      MS(bn, fn) = MS(bn, fn) + n_comp_add   ! increase the length of array

      ! add new particle to bins
      nn = n_comp_add * fac_base**(fn - 1)  ! number of new physical particles
      bin_n(bn) = bin_n(bn) + nn
      bin_g(bn) = bin_g(bn) + nn * new_v
      M = M + nn      ! increase the total number of particles

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine double_super(M, n_bin, n_fact, TDV, MS, VS, V_comp,
     $     bin_v, bin_r, bin_g, bin_n, dlnr)

C     Double number of particles in a super array.

      integer M                   ! INPUT/OUTPUT: number of particles
      integer n_bin               ! INPUT: number of bins
      integer n_fact              ! INPUT: number of allowable factors
      integer TDV                 ! INPUT: trailing dimension of VS
      integer MS(n_bin,n_fact)    ! INPUT/OUTPUT: number of superparticles
      real*8 VS(n_bin,n_fact,TDV) ! INPUT/OUTPUT: volume of physical particles
      real*8 V_comp               ! INPUT: computational volume
      real*8 bin_v(n_bin)         ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)         ! INPUT: radius of particles in bins
      real*8 bin_g(n_bin)         ! INPUT/OUTPUT: mass in bins
      integer bin_n(n_bin)        ! INPUT/OUTPUT: number in bins
      real*8 dlnr                 ! INPUT: bin scale factor

      integer i, f, b

      ! double VS and associated structures
      do b = 1,n_bin
         do f = 1,n_fact
            ! only double if we have enough space to do so
            if (2 * MS(b, f) .gt. TDV) then
               write(*,*)'ERROR: double without enough space'
               call exit(2)
            endif
            do i = 1,MS(b, f)
               VS(b, f, i + MS(b, f)) = VS(b, f, i)
            enddo
            MS(b, f) = 2 * MS(b, f)
         enddo
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

      subroutine check_super(M, n_bin, n_fact, TDV, fac_base, MS, VS,
     $     bin_v, bin_r)

C     Check that V has all particles in the correct bins.

      integer M                   ! INPUT: number of particles
      integer n_bin               ! INPUT: number of bins
      integer n_fact              ! INPUT: number of allowable factors
      integer TDV                 ! INPUT: trailing dimension of VS
      integer fac_base            ! INPUT: factor base of a superparticle
      integer MS(n_bin,n_fact)    ! INPUT: number of superparticles
      real*8 VS(n_bin,n_fact,TDV) ! INPUT: volume of physical particles

      real*8 bin_v(n_bin)         ! INPUT: volume of particles in bins
      real*8 bin_r(n_bin)         ! INPUT: radius of particles in bins

      integer i, f, b, b_check, M_check
      logical error

      error = .false.
      M_check = 0
      do b = 1,n_bin
         do f = 1,n_fact
            M_check = M_check + MS(b, f) * fac_base**(f - 1)
            do i = 1,MS(b, f)
               call particle_in_bin(VS(b, f, i), n_bin, bin_v, b_check)
               if (b .ne. b_check) then
                  write(*,'(a10,a10,a10,a12,a10)') 'b', 'f', 'i',
     $                 'VS(b, f, i)', 'b_check'
                  write(*,'(i10,i10,i10,e12.5,i10)') b, f, i, VS(b, f,
     $                 i), b_check
                  error = .true.
               endif
            enddo
         enddo
      enddo
      if (M .ne. M_check) then
         write(*,'(a10,a10)') 'M', 'M_check'
         write(*,'(i10,i10)') M, M_check
         error = .true.
      endif
      if (error) then
         write(*,*)'ERROR: check_super() failed'
         call exit(2)
      endif

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine est_k_max_super(n_bin, n_fact, fac_base, bin_v, kernel,
     $     k_max)

      integer n_bin               ! INPUT: number of bins
      integer n_fact              ! INPUT: number of allowable factors
      integer fac_base            ! INPUT: factor base of a superparticle
      real*8 bin_v(n_bin)         ! INPUT: volume of particles in bins (m^3)
      external kernel             ! INPUT: kernel function
      real*8 k_max(n_bin,n_fact,n_bin,n_fact) ! OUTPUT: maximum kernel values

      integer b1, b2, f1, f2, min_f, factor
      real*8 k
      
      do b1 = 1,n_bin
         do b2 = 1,n_bin
            call est_k_max_for_bin(n_bin, bin_v, kernel, b1, b2, k)
            do f1 = 1,n_fact
               do f2 = 1,n_fact
                  min_f = min(f1, f2)
                  factor = fac_base**(min_f - 1)
                  k_max(b1, f1, b2, f2) = factor * k
               enddo
            enddo
         enddo
      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
