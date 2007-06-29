! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Coagulation subroutines.

module mod_coagulate

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine maybe_coag_pair(M, n_bin, MH, VH, &
       n_spec, bin_v, bin_g, bin_gs, bin_n, dlnr, b1, b2, &
       del_t, k_max, kernel, env, did_coag, bin_change)
    
    ! Choose a random pair for potential coagulation and test its
    ! probability of coagulation. If it happens, do the coagulation and
    ! update all structures. The probability of a coagulation will be
    ! taken as (kernel / k_max).

    use mod_aero_data
    use mod_util
    use mod_environ
    use mod_aero_state

    integer, intent(inout) :: M         ! number of particles
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    integer, intent(inout) :: MH(n_bin) ! number of particles per bin
    type(bin_p), intent(inout) :: VH(n_bin) ! particle volumes
     
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins
    real*8, intent(inout) :: bin_g(n_bin) ! volume in bins
    real*8, intent(inout) :: bin_gs(n_bin,n_spec) ! species volume in bins
    integer, intent(inout) :: bin_n(n_bin) ! number in bins
    real*8, intent(in) :: dlnr          ! bin scale factor
    type(environ), intent(in) :: env    ! environment state
    
    integer, intent(in) :: b1           ! bin of first particle
    integer, intent(in) :: b2           ! bin of second particle
    real*8, intent(in) :: del_t         ! timestep
    real*8, intent(in) :: k_max         ! k_max scale factor
    logical, intent(out) :: did_coag    ! whether a coagulation occured
    logical, intent(out) :: bin_change  ! whether bin structure changed
    
    interface
       subroutine kernel(v1, v2, env, k)
         use mod_environ
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(environ), intent(in) :: env 
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
    
    integer s1, s2
    real*8 p, k, pv1, pv2
    
    bin_change = .false.
    did_coag = .false.
    
    if ((MH(b1) .le. 0) .or. (MH(b2) .le. 0)) then
       return
    end if
    
    call find_rand_pair(n_bin, MH, b1, b2, s1, s2)
    pv1 = particle_volume(VH(b1)%p(s1,:))
    pv2 = particle_volume(VH(b2)%p(s2,:))
    call kernel(pv1, pv2, env, k)
    p = k / k_max
    
    if (util_rand() .lt. p) then
       call coagulate(M, n_bin, MH, VH, n_spec &
            ,bin_v,bin_g, bin_gs, bin_n, dlnr, b1, s1, b2, s2, &
            env, bin_change)
       did_coag = .true.
    end if
    
  end subroutine maybe_coag_pair
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_rand_pair(n_bin, MH, b1, b2, s1, s2)
    
    ! Find a random pair of particles (b1, s1) and (b2, s2).
    
    use mod_util
    use mod_aero_state

    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: MH(n_bin)    ! number particles per bin
    integer, intent(in) :: b1           ! bin number of first particle
    integer, intent(in) :: b2           ! bin number of second particle
    integer, intent(out) :: s1          ! first rand particle 1 <= s1 <= M(b1)
    integer, intent(out) :: s2          ! second rand particle 1 <= s2 <= M(b2)
                                        !         (b1,s1) != (b2,s2)
    
    ! FIXME: rand() only returns a REAL*4, so we might not be able to
    ! generate all integers between 1 and M if M is too big.
100 s1 = int(util_rand() * dble(MH(b1))) + 1
    if ((s1 .lt. 1) .or. (s1 .gt. MH(b1))) goto 100
101 s2 = int(util_rand() * dble(MH(b2))) + 1
    if ((s2 .lt. 1) .or. (s2 .gt. MH(b2))) goto 101
    if ((b1 .eq. b2) .and. (s1 .eq. s2)) goto 101
    
  end subroutine find_rand_pair
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine coagulate(M, n_bin, MH, VH, n_spec &
       ,bin_v, bin_g, bin_gs,bin_n, dlnr, b1, s1, b2, s2, &
       env, bin_change)

    ! Join together particles (b1, s1) and (b2, s2), updating all
    ! particle and bin structures to reflect the change. bin_change is
    ! true if the used bin set changed due to the coagulation (i.e. an
    ! empty bin filled or a filled bin became empty).

    use mod_aero_data
    use mod_bin
    use mod_environ
    use mod_aero_state
    
    integer, intent(inout) :: M         ! number of particles
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    integer, intent(inout) :: MH(n_bin) ! number of particles per bin
    type(bin_p), intent(inout) :: VH(n_bin) ! particle volumes
     
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins
    real*8, intent(inout) :: bin_g(n_bin) ! volume in bins
    real*8, intent(inout) :: bin_gs(n_bin,n_spec) ! species volume in bins
    integer, intent(inout) :: bin_n(n_bin) ! number in bins
    real*8, intent(in) :: dlnr          ! bin scale factor
    
    integer, intent(in) :: b1           ! first particle (bin number)
    integer, intent(in) :: s1           ! first particle (number in bin)
    integer, intent(in) :: b2           ! second particle (bin number)
    integer, intent(in) :: s2           ! second particle (number in bin)
    type(environ), intent(in) :: env    ! environment state
    logical, intent(out) :: bin_change  ! whether an empty bin filled,
                                        ! or a filled bin became empty
    
    integer bn, i, j
    real*8 new_v(n_spec), pv1, pv2, new_v_tot

    bin_change = .false.
    
    pv1 = particle_volume(VH(b1)%p(s1,:))
    pv2 = particle_volume(VH(b2)%p(s2,:))

    ! remove s1 and s2 from bins
    bin_n(b1) = bin_n(b1) - 1
    bin_n(b2) = bin_n(b2) - 1
    bin_g(b1) = bin_g(b1) - pv1
    bin_g(b2) = bin_g(b2) - pv2
    bin_gs(b1,:) = bin_gs(b1,:) - VH(b1)%p(s1,:)
    bin_gs(b2,:) = bin_gs(b2,:) - VH(b2)%p(s2,:)
     if ((bin_n(b1) .lt. 0) .or. (bin_n(b2) .lt. 0)) then
       write(0,*)'ERROR: invalid bin_n'
       call exit(2)
    end if

    ! do coagulation in MH, VH arrays
    new_v(:) = VH(b1)%p(s1,:) + VH(b2)%p(s2,:)   ! add particle volumes
    new_v_tot = sum(new_v)

    call particle_in_bin(new_v_tot, n_bin, bin_v, bn)  ! find new bin

    ! handle a tricky corner case
    if ((b1 .eq. b2) .and. (s2 .eq. MH(b2))) then
       VH(b1)%p(s1,:) = VH(b1)%p(MH(b1)-1,:)
       MH(b1) = MH(b1) - 2
    else
       VH(b1)%p(s1,:) = VH(b1)%p(MH(b1),:) ! shift last particle into empty slot
       MH(b1) = MH(b1) - 1                 ! decrease length of array
       VH(b2)%p(s2,:) = VH(b2)%p(MH(b2),:) ! same for second particle
       MH(b2) = MH(b2) - 1
    end if
    if ((MH(b1) .lt. 0) .or. (MH(b2) .lt. 0)) then
       write(0,*)'ERROR: invalid MH'
       call exit(2)
    end if
    MH(bn) = MH(bn) + 1          ! increase the length of array
    call enlarge_bin_to(VH(bn), MH(bn))
    VH(bn)%p(MH(bn),:) = new_v   ! add the new particle at the end
    M = M - 1                    ! decrease the total number of particles
    
    ! add new particle to bins
    bin_n(bn) = bin_n(bn) + 1
    bin_g(bn) = bin_g(bn) + sum(VH(bn)%p(MH(bn),:))
    bin_gs(bn,:) = bin_gs(bn,:) + VH(bn)%p(MH(bn),:)

    ! did we empty a bin?
    if ((bin_n(b1) .eq. 0) .or. (bin_n(b2) .eq. 0)) &
         bin_change = .true.
    ! did we start a new bin?
    if ((bin_n(bn) .eq. 1) .and. (bn .ne. b1) .and. (bn .ne. b2)) &
         bin_change = .true.

    ! possibly repack memory
    call shrink_bin(MH(b1), VH(b1))
    call shrink_bin(MH(b2), VH(b2))

  end subroutine coagulate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine double(M, n_bin, MH, VH, n_spec &
       ,bin_v, bin_g, bin_gs, bin_n, dlnr, env)
    
    ! Doubles number of particles in an array.

    use mod_environ
    
    integer, intent(inout) :: M         ! number of particles
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    integer, intent(inout) :: MH(n_bin) ! number of particles per bin
    type(bin_p), intent(inout) :: VH(n_bin) ! particle volumes
     
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins
    real*8, intent(inout) :: bin_g(n_bin) ! volume in bins
    real*8, intent(inout) :: bin_gs(n_bin,n_spec) ! species volume in bins
    integer, intent(inout) :: bin_n(n_bin) ! number in bins
    type(environ), intent(inout) :: env ! environment state 
    real*8, intent(in) :: dlnr          ! bin scale factor
    
    integer i, k, i_spec
    
    ! double VH and associated structures
    do k = 1,n_bin
       call enlarge_bin_to(VH(k), 2 * MH(k))
       VH(k)%p((MH(k)+1):(2*MH(k)), :) = VH(k)%p(1:MH(k), :)
       MH(k) = 2 * MH(k)
    end do
    M = 2 * M
    env%V_comp = 2d0 * env%V_comp
    
    ! double bin structures
    bin_g = bin_g * 2d0
    bin_n = bin_n * 2
    bin_gs = bin_gs * 2d0
    
  end subroutine double
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_coagulate
