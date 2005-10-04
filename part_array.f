C part_array.f
C
C utility functions for handling V array of particle volumes

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     
      subroutine find_rand_pair(V, MM, M_comp, s1, s2)
      
      integer MM      ! INPUT: dimension of V
      integer M_comp  ! INPUT: maximum index of non-zero entries in V
      real*8 V(MM)    ! INPUT: array of particle volumes
      integer s1, s2  ! OUTPUT: s1 and s2 are not equal, random
                      !         particles with V(s1/s2) != 0

 100  s1 = int(rand() * M_comp) + 1
      s2 = int(rand() * M_comp) + 1
      if ((s1 .gt. M_comp) .or. (V(s1) .eq. 0) .or.
     &     (s2 .gt. M_comp) .or. (V(s2) .eq. 0) .or.
     &     (s1 .eq. s2)) then
         goto 100
      endif

      return
      end

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine compress(MM, V)

      integer MM   ! INPUT
      real*8 V(MM) ! INPUT/OUTPUT: on exit, all non-zero entries are
                   !               at the beginning, followed by all
                   !               zeros.

      integer i, i_w, i_v

      i_w = 1
      do i_v = 1,MM
         if (V(i_v) .ne. 0.) then
            V(i_w) = V(i_v)
            i_w = i_w + 1
         endif
      enddo

      do i = i_w + 1, MM
         V(i) = 0.
      enddo

      return
      end
