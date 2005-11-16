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
