! -*- mode: f90; -*-
! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Functions to deal with a particle array VH that is stored by
! bin. VH(i_bin,i) is the i-th particle in the i_bin-th bin.
!
! FIXME: MH and bin_n are pretty much identical. Probably best to
! ignore it, for symmetry with non-hybrid code, and because in
! superparticle code there is a difference.

module mod_array_hybrid

  type bin_p
     real*8, dimension(:,:), pointer :: p ! particle volumes
     ! dimension of p is (# particles in bin) x n_spec
  end type bin_p

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_hybrid(n_spec, MH, VH)

    ! Initializes a hybrid array to have zero particles in each bin.

    integer, intent(in) :: n_spec       ! number of species
    integer, intent(out) :: MH(:)       ! number of particles per bin
    type(bin_p), intent(inout) :: VH(size(MH)) ! particle volumes
    
    integer :: n_bin
    integer i

    n_bin = size(VH)
    do i = 1,n_bin
       allocate(VH(i)%p(0, n_spec))
    end do
    MH = 0

  end subroutine init_hybrid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine enlarge_bin(bin)

    ! Enlarges the given bin (which must be allocated) by at least one
    ! element (currently doubles the length).

    type(bin_p), intent(inout) :: bin  ! bin data

    integer :: n_part, n_spec, new_n_part
    real*8, dimension(:,:), pointer :: new_p

    n_part = size(bin%p, 1)
    n_spec = size(bin%p, 2)
    new_n_part = max(n_part * 2, n_part + 1)
    allocate(new_p(new_n_part, n_spec))
    new_p(1:n_part,:) = bin%p
    deallocate(bin%p)
    bin%p => new_p
    
  end subroutine enlarge_bin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine shrink_bin(n_used, bin)

    ! Possibly shrinks the storage of the given bin, ensuring that it
    ! is at least of length n_used (currently halves the length or
    ! does nothing).

    integer, intent(in) :: n_used      ! number of used entries in bin
    type(bin_p), intent(inout) :: bin  ! bin data

    integer :: n_part, n_spec, new_n_part
    real*8, dimension(:,:), pointer :: new_p

    n_part = size(bin%p, 1)
    n_spec = size(bin%p, 2)
    new_n_part = n_part / 2
    if (n_used <= new_n_part) then
       allocate(new_p(new_n_part, n_spec))
       new_p(:,:) = bin%p(1:new_n_part,:)
       deallocate(bin%p)
       bin%p => new_p
    end if

  end subroutine shrink_bin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine array_to_hybrid(MM, M, V, n_spec, n_bin, bin_v, MH &
       ,VH)
    
    ! Convert a standard linear particle array V to a hybrid particle
    ! array VH stored by bins.

    use mod_material
    use mod_array
    use mod_bin
    
    integer, intent(in) :: MM           !  physical dimension of V
    integer, intent(in) :: M            !  logical dimension of V
    integer, intent(in) :: n_spec       !  number of species
    real*8, intent(inout) :: V(MM,n_spec)  !  particle volumes
    integer, intent(in) :: n_bin        !  number of bins
    real*8, intent(in) :: bin_v(n_bin)  !  volume of particles in bins
    integer, intent(out) :: MH(n_bin)    !  number of particles per bin
    type(bin_p), intent(out) :: VH(n_bin) !  particle volumes in hybrid array
    
    integer i, j, k
    real*8 pv
    
    do k = 1,n_bin
       MH(k) = 0
    enddo
    
    do i = 1,M
       pv = particle_volume(V(i,:))     
       call particle_in_bin(pv, n_bin, bin_v, k)
       MH(k) = MH(k) + 1
       if (MH(k) .gt. size(VH(k)%p,1)) then
          call enlarge_bin(VH(k))
       endif
       VH(k)%p(MH(k),:) = V(i,:)
    enddo
    
  end subroutine array_to_hybrid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine moments_hybrid(n_bin, n_spec, MH, VH, bin_v, &
       bin_g,bin_gs, bin_n, dlnr)
    
    ! Create the bin number and mass arrays from VH.

    use mod_material
    use mod_array
    
    integer, intent(in) :: n_bin        !  number of bins
    integer, intent(in) :: n_spec       !  number of species
    integer, intent(in) :: MH(n_bin)    !  number of particles per bin
    type(bin_p), intent(in) :: VH(n_bin) !  particle volumes
    real*8, intent(in) :: bin_v(n_bin)  !  volume of particles in bins
    real*8, intent(out) :: bin_g(n_bin)  !  volume in bins
    real*8, intent(out) :: bin_gs(n_bin,n_spec)  !  species volume in bins
    integer, intent(out) :: bin_n(n_bin) !  number in bins
    real*8, intent(in) :: dlnr          !  bin scale factor
    
    integer b, j, s
    
    bin_g = 0d0
    bin_gs = 0d0
    do b = 1,n_bin
       do j = 1,MH(b)
          bin_g(b) = bin_g(b) + particle_volume(VH(b)%p(j,:))
          bin_gs(b,:) = bin_gs(b,:) + VH(b)%p(j,:)
       enddo
    enddo
    bin_n = MH
   
  end subroutine moments_hybrid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine resort_array_hybrid(n_bin, n_spec, MH, VH, bin_v, &
        dlnr)
    
    ! Takes a VH array where the particle volumes might no longer be
    ! correct for the bins they are in and resorts it so that every
    ! particle is in the correct bin.

    use mod_material
    use mod_array
    use mod_bin
    
    integer, intent(in) :: n_bin ! number of bins
    integer, intent(in) :: n_spec ! number of species
    integer, intent(inout) :: MH(n_bin) ! number of particles per bin
    type(bin_p), intent(inout) :: VH(n_bin) ! particle volumes (m^3)
    real*8, intent(in) :: bin_v(n_bin) ! volume of particles in bins (m^3)
    real*8, intent(in) :: dlnr ! bin scale factor
    
    integer bin, j, new_bin, k
    real*8 pv
    
    ! FIXME: the approach here is inefficient because we might
    ! reprocess particles. For example, if we are doing bin 1 and we
    ! shift a particle up to bin 2, when we do bin 2 we will reprocess
    ! it. It seems to be more trouble than it's worth to worry about
    ! this yet, however.
    
    do bin = 1,n_bin
       j = 1
       do while (j .le. MH(bin))
          ! find the new volume and new bin
          pv = particle_volume(VH(bin)%p(j,:))
          call particle_in_bin(pv, n_bin, bin_v, new_bin)
          
          ! if the bin number has changed, move the particle
          if (bin .ne. new_bin) then
             ! move the particle to the new bin, leaving a hole
             MH(new_bin) = MH(new_bin) + 1
             if (MH(new_bin) .gt. size(VH(new_bin)%p,1)) then
                call enlarge_bin(VH(new_bin))
             end if
             VH(new_bin)%p(MH(new_bin),:) = VH(bin)%p(j,:)
             
             ! copy the last particle in the current bin into the hole
             ! if the hole isn't in fact the last particle
             if (j .lt. MH(bin)) then
                VH(bin)%p(j,:) = VH(bin)%p(MH(bin),:)
             end if
             MH(bin) = MH(bin) - 1
             if (MH(bin) .lt. 0) then
                write(*,*) 'ERROR: invalid MH in bin ', bin
                call exit(2)
             end if
             
             ! in this case, don't advance j, so that we will still
             ! process the particle we just moved into the hole
          else
             ! if we didn't move the particle, advance j to process
             ! the next particle
             j = j + 1
          end if
       end do
    end do

    ! now shrink the bin storage if necessary
    do bin = 1,n_bin
       call shrink_bin(MH(bin), VH(bin))
    end do
    
  end subroutine resort_array_hybrid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine maybe_coag_pair_hybrid(M, n_bin, MH, VH, V_comp, &
       n_spec, bin_v, bin_g, bin_gs, bin_n, dlnr, b1, b2, &
       del_t, k_max, kernel, env, did_coag, bin_change)
    
    ! Choose a random pair for potential coagulation and test its
    ! probability of coagulation. If it happens, do the coagulation and
    ! update all structures. The probability of a coagulation will be
    ! taken as (kernel / k_max).

    use mod_material
    use mod_util
    use mod_environ

    integer, intent(inout) :: M            !  number of particles
    integer, intent(in) :: n_bin        !  number of bins
    integer, intent(in) :: n_spec       !  number of species
    integer, intent(inout) :: MH(n_bin)    !  number of particles per bin
    type(bin_p), intent(inout) :: VH(n_bin) !  particle volumes
    real*8, intent(in) :: V_comp        !  computational volume
    
    real*8, intent(in) :: bin_v(n_bin)  !  volume of particles in bins
    real*8, intent(inout) :: bin_g(n_bin)  !  volume in bins
    real*8, intent(inout) :: bin_gs(n_bin,n_spec)  !  species volume in bins
    integer, intent(inout) :: bin_n(n_bin) !  number in bins
    real*8, intent(in) :: dlnr          !  bin scale factor
    type(environ), intent(in) :: env        ! environment state
    
    integer, intent(in) :: b1           !  bin of first particle
    integer, intent(in) :: b2           !  bin of second particle
    real*8, intent(in) :: del_t         !  timestep
    real*8, intent(in) :: k_max         !  k_max scale factor
    ! external, intent(in) :: kernel      !  kernel function
    logical, intent(out) :: did_coag     !  whether a coagulation occured
    logical, intent(out) :: bin_change   !  whether bin structure changed
    
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
    endif
    
    call find_rand_pair_hybrid(n_bin, MH, b1, b2, s1, s2)
    pv1 = particle_volume(VH(b1)%p(s1,:))
    pv2 = particle_volume(VH(b2)%p(s2,:))
    call kernel(pv1, pv2, env, k)
    p = k / k_max
    
    if (util_rand() .lt. p) then
       call coagulate_hybrid(M, n_bin, MH, VH, V_comp, n_spec &
            ,bin_v,bin_g, bin_gs, bin_n, dlnr, b1, s1, b2, s2, &
            bin_change)
       did_coag = .true.
    endif
    
  end subroutine maybe_coag_pair_hybrid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_rand_pair_hybrid(n_bin, MH, b1, b2, s1, s2)
    
    ! Find a random pair of particles (b1, s1) and (b2, s2).
    
    use mod_util

    integer, intent(in) :: n_bin     !  number of bins
    integer, intent(in) :: MH(n_bin) !  number particles per bin
    integer, intent(in) :: b1        !  bin number of first particle
    integer, intent(in) :: b2        !  bin number of second particle
    integer, intent(out) :: s1       !  first random particle 1 <= s1 <= M(b1)
    integer, intent(out) :: s2       !  second random particle 1 <= s2 <= M(b2)
                                     !         (b1,s1) != (b2,s2)
    
    ! FIXME: rand() only returns a REAL*4, so we might not be able to
    ! generate all integers between 1 and M if M is too big.
100 s1 = int(util_rand() * dble(MH(b1))) + 1
    if ((s1 .lt. 1) .or. (s1 .gt. MH(b1))) goto 100
101 s2 = int(util_rand() * dble(MH(b2))) + 1
    if ((s2 .lt. 1) .or. (s2 .gt. MH(b2))) goto 101
    if ((b1 .eq. b2) .and. (s1 .eq. s2)) goto 101
    
  end subroutine find_rand_pair_hybrid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine coagulate_hybrid(M, n_bin, MH, VH, V_comp, n_spec &
       ,bin_v, bin_g, bin_gs,bin_n, dlnr, b1, s1, b2, s2, &
       bin_change)
    
    ! Join together particles (b1, s1) and (b2, s2), updating all
    ! particle and bin structures to reflect the change. bin_change is
    ! true if the used bin set changed due to the coagulation (i.e. an
    ! empty bin filled or a filled bin became empty).

    use mod_material
    use mod_array
    use mod_bin
    
    integer, intent(inout) :: M            !  number of particles
    integer, intent(in) :: n_bin        !  number of bins
    integer, intent(in) :: n_spec       !  number of species
    integer, intent(inout) :: MH(n_bin)    !  number of particles per bin
    type(bin_p), intent(inout) :: VH(n_bin) !  particle volumes
    real*8, intent(in) :: V_comp        !  computational volume
    
    real*8, intent(in) :: bin_v(n_bin)  !  volume of particles in bins
    real*8, intent(inout) :: bin_g(n_bin)  !  volume in bins
    real*8, intent(inout) :: bin_gs(n_bin,n_spec)  !  species volume in bins
    integer, intent(inout) :: bin_n(n_bin) !  number in bins
    real*8, intent(in) :: dlnr          !  bin scale factor
    
    integer, intent(in) :: b1           !  first particle (bin number)
    integer, intent(in) :: s1           !  first particle (number in bin)
    integer, intent(in) :: b2           !  second particle (bin number)
    integer, intent(in) :: s2           !  second particle (number in bin)
    logical, intent(out) :: bin_change   !  whether an empty bin filled,
    !         or a filled bin became empty
    
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
       write(*,*)'ERROR: invalid bin_n'
       call exit(2)
    endif

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
       write(*,*)'ERROR: invalid MH'
       call exit(2)
    end if
    MH(bn) = MH(bn) + 1          ! increase the length of array
    if (MH(bn) > size(VH(bn)%p,1)) then
       call enlarge_bin(VH(bn))
    end if
    VH(bn)%p(MH(bn),:) = new_v  ! add the new particle at the end
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

  end subroutine coagulate_hybrid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine double_hybrid(M, n_bin, MH, VH, V_comp, n_spec &
       ,bin_v, bin_g, bin_gs, bin_n, dlnr)
    
    ! Double number of particles in a hybrid array.
    
    integer, intent(inout) :: M            !  number of particles
    integer, intent(in) :: n_bin        !  number of bins
    integer, intent(in) :: n_spec       !  number of species
    integer, intent(inout) :: MH(n_bin)    !  number of particles per bin
    type(bin_p), intent(inout) :: VH(n_bin) !  particle volumes
    real*8, intent(inout) :: V_comp        !  computational volume
    
    real*8, intent(in) :: bin_v(n_bin)  !  volume of particles in bins
    real*8, intent(inout) :: bin_g(n_bin)  !  volume in bins
    real*8, intent(inout) :: bin_gs(n_bin,n_spec) !  species volume in bins
    integer, intent(inout) :: bin_n(n_bin) !  number in bins
    real*8, intent(in) :: dlnr          !  bin scale factor
    
    integer i, k, i_spec
    
    ! double VH and associated structures
    do k = 1,n_bin
       do i = 1,MH(k)
          if (i + MH(k) > size(VH(k)%p,1)) then
             call enlarge_bin(VH(k))
          end if
          VH(k)%p(i + MH(k), :) = VH(k)%p(i, :)
       enddo
       MH(k) = 2 * MH(k)
    enddo
    M = 2 * M
    V_comp = 2d0 * V_comp
    
    ! double bin structures
    bin_g = bin_g * 2d0
    bin_n = bin_n * 2
    bin_gs = bin_gs * 2d0
    
  end subroutine double_hybrid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine check_hybrid(M, n_bin, n_spec, MH, VH, bin_v, &
        bin_g, bin_gs, bin_n, dlnr)
    
    ! Check that VH has all particles in the correct bins and that the
    ! bin numbers and masses are correct.

    use mod_material
    use mod_array
    use mod_util
    use mod_bin
    
    integer, intent(in) :: M            !  number of particles
    integer, intent(in) :: n_bin        !  number of bins
    integer, intent(in) :: n_spec       !  number of species
    integer, intent(in) :: MH(n_bin)    !  number of particles per bin
    type(bin_p), intent(in) :: VH(n_bin) !  particle volumes
    
    real*8, intent(in) :: bin_v(n_bin)       !  volume of particles in bins (m^3)
    real*8, intent(out) :: bin_g(n_bin)       !  volume in bins  
    real*8, intent(out) :: bin_gs(n_bin,n_spec) !  species volume in bins             
    integer, intent(out) :: bin_n(n_bin)      !  number in bins
    real*8, intent(in) :: dlnr               !  bin scale factor
    
    real*8 pv, check_bin_g, check_bin_gs(n_spec), vol_tol
    integer i, k, k_check, M_check, s
    logical error
    
    error = .false.
    
    ! check that all particles are in the correct bins
    do k = 1,n_bin
       do i = 1,MH(k)
          pv = particle_volume(VH(k)%p(i,:))
          call particle_in_bin(pv, n_bin, bin_v, k_check)
          if (k .ne. k_check) then
             write(*,'(a10,a10,a12,a10)') 'k', 'i', 'VH(k, i)', &
                  'k_check'
             write(*,'(i10,i10,e12.5,i10)') k, i, pv, k_check
             error = .true.
          endif
       enddo
    enddo
    
    ! check that the total number of particles is correct
    M_check = 0
    do k = 1,n_bin
       M_check = M_check + MH(k)
    enddo
    if (M .ne. M_check) then
       write(*,'(a10,a10)') 'M', 'M_check'
       write(*,'(i10,i10)') M, M_check
       error = .true.
    endif
    
    ! check the bin_n array
    do k = 1,n_bin
       if (MH(k) .ne. bin_n(k)) then
          write(*,'(a10,a10,a10)') 'k', 'MH(k)', 'bin_n(k)'
          write(*,'(i10,i10,i10)') k, MH(k), bin_n(k)
       endif
    enddo
    
    ! check the bin_g array
    do k = 1,n_bin
       check_bin_g = 0d0
       do i = 1,MH(k)
          pv = particle_volume(VH(k)%p(i,:))
          check_bin_g = check_bin_g + pv
       enddo
       vol_tol = bin_v(k) / 1d6 ! abs tolerance 1e6 less than single particle
       if (.not. almost_equal_abs(check_bin_g, bin_g(k), vol_tol)) then
          write(*,'(a10,a15,a15)') 'k', 'check_bin_g', 'bin_g(k)'
          write(*,'(i10,e15.5,e15.5)') k, check_bin_g, bin_g(k)
          error = .true.
       endif
    enddo
    
    ! check the bin_gs array
    do k = 1,n_bin
       check_bin_gs = sum(VH(k)%p(1:MH(k),:), 1)
       vol_tol = bin_v(k) / 1d3 ! abs tolerance 1e3 less than single particle
       do s = 1,n_spec
          if (.not. almost_equal_abs(check_bin_gs(s), bin_gs(k,s), &
                                     vol_tol)) then
             write(*,'(a10,a10,a20,a15)') 'k', 's', 'check_bin_gs(s)', &
                  'bin_gs(k,s)'
             write(*,'(i10,i10,e20.5,e15.5)') k, s, check_bin_gs(s), &
                  bin_gs(k,s)
             error = .true.
          endif
       enddo
    enddo
    
    if (error) then
       write(*,*) 'ERROR: check_hybrid() failed'
       call exit(2)
    endif
    
  end subroutine check_hybrid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_array_hybrid
