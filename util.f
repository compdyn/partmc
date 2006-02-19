C     Common utility functions.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine sum_int_1d(size1, array, sum)

      integer size1             ! INPUT: size of the array
      integer array(size1)      ! INPUT: array of numbers
      integer sum               ! OUTPUT: sum of array(i)
      
      integer i
      
      sum = 0
      do i = 1,size1
         sum = sum + array(i)
      enddo
      
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine sum_int_2d(size1, size2, array, sum)

      integer size1               ! INPUT: leading size of the array
      integer size2               ! INPUT: trailing size of the array
      integer array(size1,size2)  ! INPUT: array of numbers
      integer sum                 ! OUTPUT: sum of array(i,j)

      integer i, j

      sum = 0
      do i = 1,size1
         do j = 1,size2
            sum = sum + array(i, j)
         enddo
      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine max_int_2d(size1, size2, array, max)

      integer size1               ! INPUT: leading size of the array
      integer size2               ! INPUT: trailing size of the array
      integer array(size1,size2)  ! INPUT: array of numbers
      integer max                 ! OUTPUT: max of array(i,j)

      integer i, j

      max = 0
      do i = 1,size1
         do j = 1,size2
            if (array(i, j) .gt. max) then
               max = array(i, j)
            endif
         enddo
      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine check_event(time, interval, last_time, do_event)

      real*8 time       ! INPUT: cubin_rent time
      real*8 interval   ! INPUT: how often the event should be done
      real*8 last_time  ! INPUT/OUTPUT: when the event was last done
      logical do_event  ! OUTPUT: whether the event should be done

      real*8 interval_below

      if (time .eq. 0d0) then
         last_time = 0d0
         do_event = .true.
      else
         interval_below = aint(time / interval) * interval
         if (last_time .lt. interval_below) then
            last_time = time
            do_event = .true.
         else
            do_event = .false.
         endif
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
