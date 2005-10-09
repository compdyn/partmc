C Process output data files.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program process_out

      integer n_bin_max, n_loop_max, n_time_max
      parameter (n_bin_max = 400)  ! maximum number of bins
      parameter (n_loop_max = 50)  ! maximum number of loops
      parameter (n_time_max = 100) ! maximum number of times
      
      integer f_in, f_out_num, f_out_mass
      integer f_out_num_avg, f_out_mass_avg
      parameter (f_in = 20)           ! input
      parameter (f_out_num = 21)      ! output number
      parameter (f_out_mass = 22)     ! output mass
      parameter (f_out_num_avg = 23)  ! output number average
      parameter (f_out_mass_avg = 24) ! output mass average
      
      integer n_bin, n_loop, n_time
      character name_in*50, name_out_num*50, name_out_mass*50
      character name_out_num_avg*50, name_out_mass_avg*50
      character dum*100, n_loop_str*10, n_time_str*10

      real*8 rr(n_bin_max), g(n_loop_max, n_time_max, n_bin_max)
      real*8 n(n_loop_max, n_time_max, n_bin_max)
      real*8 g_avg(n_time_max, n_bin_max), n_avg(n_time_max, n_bin_max)

      real*8 time
      integer i, i_loop, i_time, i_bin

      ! check there is exactly one commandline argument
      if (iargc() .ne. 1) then
         write(6,*)'Usage: process_out <filename.d>'
         call exit(2)
      endif

      ! get and check first commandline argument (must be "filename.d")
      call getarg(1, name_in)
      i = len_trim(name_in)
      if (i .gt. 40) then
         write(6,*)'ERROR: filename too long'
         call exit(2)
      endif
      if ((name_in(i:i) .ne. 'd') .or.
     &     (name_in((i-1):(i-1)) .ne. '.')) then
         write(6,*)'ERROR: Filename must end in .d'
         call exit(2)
      endif

      ! compute names of output files
      name_out_num = name_in
      name_out_mass = name_in
      name_out_num_avg = name_in
      name_out_mass_avg = name_in
      name_out_num((i-1):) = '_num.d'
      name_out_mass((i-1):) = '_mass.d'
      name_out_num_avg((i-1):) = '_num_avg.d'
      name_out_mass_avg((i-1):) = '_mass_avg.d'

      write(6,*)'name_in = ', name_in
      write(6,*)'name_out_num = ', name_out_num
      write(6,*)'name_out_mass = ', name_out_mass
      write(6,*)'name_out_num_avg = ', name_out_num_avg
      write(6,*)'name_out_mass_avg = ', name_out_mass_avg

      ! open files
      open(f_in, file=name_in)
      open(f_out_num, file=name_out_num)
      open(f_out_mass, file=name_out_mass)
      open(f_out_num_avg, file=name_out_num_avg)
      open(f_out_mass_avg, file=name_out_mass_avg)

      ! read and check dimensions
      read(f_in, '(a10,i10)'), dum, n_loop
      read(f_in, '(a10,i10)'), dum, n_bin
      read(f_in, '(a10,i10)'), dum, n_time

      if (n_loop .gt. n_loop_max) then
         write(6,*)'ERROR: n_loop too large'
         call exit(2)
      endif
      if (n_bin .gt. n_bin_max) then
         write(6,*)'ERROR: n_bin too large'
         call exit(2)
      endif
      if (n_time .gt. n_time_max) then
         write(6,*)'ERROR: n_time too large'
         call exit(2)
      endif

      write(6,*)'n_loop = ', n_loop
      write(6,*)'n_bin = ', n_bin
      write(6,*)'n_time = ', n_time

      write(n_loop_str, '(i10)'), n_loop
      write(n_time_str, '(i10)'), n_time

      ! read all data
      do i_loop = 1,n_loop
         do i_time = 1,n_time
            read(f_in, '(a10,e14.5)'), dum, time
            do i_bin = 1,n_bin
               read(f_in, '(i8,3e14.5)'), i, rr(i_bin),
     &              n(i_loop, i_time, i_bin), g(i_loop, i_time, i_bin)
            enddo
         enddo
      enddo

      ! compute averages
      do i_time = 1,n_time
         do i_bin = 1,n_bin
            g_avg(i_time, i_bin) = 0.
            n_avg(i_time, i_bin) = 0.
            do i_loop = 1,n_loop
               g_avg(i_time, i_bin) = g_avg(i_time, i_bin)
     &              + g(i_loop, i_time, i_bin)
               n_avg(i_time, i_bin) = n_avg(i_time, i_bin)
     &              + n(i_loop, i_time, i_bin)
            enddo
            g_avg(i_time, i_bin) = g_avg(i_time, i_bin) / n_loop
            n_avg(i_time, i_bin) = n_avg(i_time, i_bin) / n_loop
         enddo
      enddo

      ! output data
      do i_time = 1,n_time
         write(f_out_num, '(//,a10,i10)'), 'time', i_time - 1
         write(f_out_mass, '(//,a10,i10)'), 'time', i_time - 1
         do i_bin = 1,n_bin
            write(f_out_num, '(i8,e14.5,'//n_loop_str//'e14.5)'),
     &           i_bin, rr(i_bin),
     &           (n(i_loop, i_time, i_bin), i_loop = 1,n_loop)
            write(f_out_mass, '(i8,e14.5,'//n_loop_str//'e14.5)'),
     &           i_bin, rr(i_bin),
     &           (g(i_loop, i_time, i_bin), i_loop = 1,n_loop)
         enddo
      enddo

      do i_bin = 1,n_bin
         write(f_out_num_avg, '(i8,e14.5,'//n_time_str//'e14.5)'),
     &        i_bin, rr(i_bin),
     &        (n_avg(i_time, i_bin), i_time = 1,n_time)
         write(f_out_mass_avg, '(i8,e14.5,'//n_time_str//'e14.5)'),
     &        i_bin, rr(i_bin),
     &        (g_avg(i_time, i_bin), i_time = 1,n_time)
      enddo
      

c       integer n_bin,i_loop,tmax
c       parameter (n_bin=160,i_loop=1,tmax=11)
c
c       real*8 eta(n_bin),g(i_loop,tmax,n_bin),n(i_loop,tmax,n_bin)
c       real*8 g_bar(tmax,n_bin),n_bar(tmax,n_bin)
c       integer k, l, i, mm
c       character*20 dum
c
c       open(20,file='mc_fix.d')
c       open(30,file='mc_n.dtest2')
c       open(40,file='mc_m.dtest2')
c       open(50,file='mcr_n.dtest2')
c       open(60,file='mcr_m.dtest2')
c      
c
c       do k=1,tmax
c          do l=1,n_bin
c             n_bar(k,l) = 0.
c             g_bar(k,l) = 0.
c          enddo
c       enddo
c
c       do i=1,i_loop
c
c          read(20,'(a)')dum
c          write(6,'(a)')dum
c          do k=1,tmax
c             read(20,'(a)')dum
c             write(6,'(a)')dum
c                 do l=1,n_bin
c                    read(20,'(i4,5e14.5)')
c     *                      mm,eta(l),n(i,k,l),g(i,k,l)
c                    write(6,*)mm,eta(l),n(i,k,l),g(i,k,l)
c                  enddo
c           enddo 
c        enddo
c
c       do k=1,tmax
c          do l=1,n_bin
c             do i=1,i_loop
c                    g_bar(k,l) = g_bar(k,l)+g(i,k,l)
c                    n_bar(k,l) = n_bar(k,l)+n(i,k,l)
c             enddo
c          
c             g_bar(k,l) = g_bar(k,l)/i_loop
c             n_bar(k,l) = n_bar(k,l)/i_loop
c          enddo
c        enddo
c
c        do k=1,tmax
c           write(50,*)'time= ',k-1
c           write(60,*)'time= ',k-1
c       
c           do l=1,n_bin
c              write(50,'(i4,100e14.5)')l,eta(l),(n(i,k,l),i=1,i_loop),
c     &                 n_bar(k,l)
c              write(60,'(i4,100e14.5)')l,eta(l),(g(i,k,l),i=1,i_loop),
c     &                 g_bar(k,l)
c           enddo
c        enddo
c 
c       do l=1,n_bin
c            write(30,'(i4,100e14.5)')l,eta(l)*1.e+06
c     &           ,(n_bar(k,l),k=1,tmax)
c            write(40,'(i4,100e14.5)')l,eta(l)*1.e+06
c     &           ,(g_bar(k,l)*1000.,k=1,tmax)
c       enddo


        end
       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
