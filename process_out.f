C Process output data files.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program process_out

      integer n_bin_max, n_loop_max, n_time_max
      integer n_spec_max
      parameter (n_bin_max = 400)   ! maximum number of bins
      parameter (n_loop_max = 100)  ! maximum number of loops
      parameter (n_time_max = 200)  ! maximum number of times
      parameter (n_spec_max = 5)    ! maximum number of species
      
      integer f_in
      integer f_out_time, f_out_time_avg
      integer f_out_num, f_out_mass
      integer f_out_num_avg, f_out_mass_avg
      integer f_out_temp, f_out_temp_avg
      integer f_out_rh, f_out_rh_avg
      parameter (f_in = 20)           ! input
      parameter (f_out_num = 21)      ! output number
      parameter (f_out_mass = 22)     ! output mass
      parameter (f_out_temp = 23)     ! output temperature
      parameter (f_out_rh = 24)       ! output relative humidity
      parameter (f_out_time = 25)     ! output time
      parameter (f_out_num_avg = 26)  ! output number average
      parameter (f_out_mass_avg = 27) ! output mass average
      parameter (f_out_temp_avg = 28) ! output temperature average
      parameter (f_out_rh_avg = 29)   ! output relative humidity average
      parameter (f_out_time_avg = 30) ! output time average
      
      integer n_bin, n_loop, n_time, n_spec
      character name_in*50
      character name_out_time*50, name_out_time_avg*50
      character name_out_num*50, name_out_mass*50
      character name_out_temp*50, name_out_rh*50
      character name_out_num_avg*50, name_out_mass_avg*50
      character name_out_temp_avg*50, name_out_rh_avg*50
      character dum*100, n_loop_str*10, n_time_str*10

      real*8 time(n_loop_max, n_time_max), time_avg(n_time_max)
      real*8 bin_r(n_bin_max), bin_g(n_loop_max, n_time_max, n_bin_max)
      real*8 bin_gs(n_loop_max, n_time_max, n_bin_max, n_spec_max)
      real*8 n(n_loop_max, n_time_max, n_bin_max)
      real*8 g_avg(n_time_max, n_bin_max)
      real*8 gs_avg(n_time_max,n_bin_max,n_spec_max)
      real*8 n_avg(n_time_max, n_bin_max)
      real*8 temp(n_loop_max, n_time_max), temp_avg(n_time_max)
      real*8 rh(n_loop_max, n_time_max), rh_avg(n_time_max)

      integer i, i_loop, i_time, i_bin, i_spec

      ! check there is exactly one commandline argument
      if (iargc() .ne. 1) then
         write(6,*) 'Usage: process_out <filename.d>'
         call exit(2)
      endif

      ! get and check first commandline argument (must be "filename.d")
      call getarg(1, name_in)
      i = len_trim(name_in)
      if (i .gt. 40) then
         write(6,*) 'ERROR: filename too long'
         call exit(2)
      endif
      if ((name_in(i:i) .ne. 'd') .or.
     &     (name_in((i-1):(i-1)) .ne. '.')) then
         write(6,*) 'ERROR: Filename must end in .d'
         call exit(2)
      endif

      ! compute names of output files
      name_out_num = name_in
      name_out_mass = name_in
      name_out_temp = name_in
      name_out_rh = name_in
      name_out_time = name_in
      name_out_num_avg = name_in
      name_out_mass_avg = name_in
      name_out_temp_avg = name_in
      name_out_rh_avg = name_in
      name_out_time_avg = name_in
      name_out_num((i-1):) = '_num.d'
      name_out_mass((i-1):) = '_mass.d'
      name_out_temp((i-1):) = '_temp.d'
      name_out_rh((i-1):) = '_rh.d'
      name_out_time((i-1):) = '_time.d'
      name_out_num_avg((i-1):) = '_num_avg.d'
      name_out_mass_avg((i-1):) = '_mass_avg.d'
      name_out_temp_avg((i-1):) = '_temp_avg.d'
      name_out_rh_avg((i-1):) = '_rh_avg.d'
      name_out_time_avg((i-1):) = '_time_avg.d'

      write(6,*) 'name_in = ', name_in
      write(6,*) 'name_out_num = ', name_out_num
      write(6,*) 'name_out_mass = ', name_out_mass
      write(6,*) 'name_out_temp = ', name_out_temp
      write(6,*) 'name_out_rh = ', name_out_rh
      write(6,*) 'name_out_time = ', name_out_time
      write(6,*) 'name_out_num_avg = ', name_out_num_avg
      write(6,*) 'name_out_mass_avg = ', name_out_mass_avg
      write(6,*) 'name_out_temp_avg = ', name_out_temp_avg
      write(6,*) 'name_out_rh_avg = ', name_out_rh_avg
      write(6,*) 'name_out_time_avg = ', name_out_time_avg

      ! open files
      open(f_in, file=name_in)
      open(f_out_num, file=name_out_num)
      open(f_out_mass, file=name_out_mass)
      open(f_out_temp, file=name_out_temp)
      open(f_out_rh, file=name_out_rh)
      open(f_out_time, file=name_out_time)
      open(f_out_num_avg, file=name_out_num_avg)
      open(f_out_mass_avg, file=name_out_mass_avg)
      open(f_out_temp_avg, file=name_out_temp_avg)
      open(f_out_rh_avg, file=name_out_rh_avg)
      open(f_out_time_avg, file=name_out_time_avg)

      ! read and check dimensions
      read(f_in, '(a10,i10)') dum, n_loop
      read(f_in, '(a10,i10)') dum, n_bin
      read(f_in, '(a10,i10)') dum, n_time
      read(f_in, '(a10,i10)') dum, n_spec

      if (n_loop .gt. n_loop_max) then
         write(6,*) 'ERROR: n_loop too large'
         call exit(2)
      endif
      if (n_bin .gt. n_bin_max) then
         write(6,*) 'ERROR: n_bin too large'
         call exit(2)
      endif
      if (n_time .gt. n_time_max) then
         write(6,*) 'ERROR: n_time too large'
         call exit(2)
      endif

      write(6,*) 'n_loop = ', n_loop
      write(6,*) 'n_bin =  ', n_bin
      write(6,*) 'n_time = ', n_time
      write(6,*) 'n_spec = ', n_spec

      write(n_loop_str, '(i10)') n_loop
      write(n_time_str, '(i10)') n_time

      ! read all data
      do i_loop = 1,n_loop
         do i_time = 1,n_time
            read(f_in, '(a10,e20.10)') dum, time(i_loop, i_time)
            read(f_in, '(a10,e20.10)') dum, temp(i_loop, i_time)
            read(f_in, '(a10,e20.10)') dum, rh(i_loop, i_time)
            do i_bin = 1,n_bin
               read(f_in, '(i10,50e20.10)') i, bin_r(i_bin),
     &              n(i_loop, i_time, i_bin),
     &              bin_g(i_loop, i_time, i_bin),
     &              (bin_gs(i_loop,i_time, i_bin, i_spec),
     &              i_spec=1,n_spec)
            enddo
         enddo
      enddo

      ! compute simple loop averages
      do i_time = 1,n_time
         temp_avg(i_time) = 0d0
         rh_avg(i_time) = 0d0
         time_avg(i_time) = 0d0
         do i_loop = 1,n_loop
            temp_avg(i_time) = temp_avg(i_time) + temp(i_loop, i_time)
            rh_avg(i_time) = rh_avg(i_time) + rh(i_loop, i_time)
            time_avg(i_time) = time_avg(i_time) + time(i_loop, i_time)
         enddo
         temp_avg(i_time) = temp_avg(i_time) / dble(n_loop)
         rh_avg(i_time) = rh_avg(i_time) / dble(n_loop)
         time_avg(i_time) = time_avg(i_time) / dble(n_loop)
      enddo

      ! compute binned loop averages
      do i_time = 1,n_time
         do i_bin = 1,n_bin
            g_avg(i_time, i_bin) = 0d0
            n_avg(i_time, i_bin) = 0d0
            do i_spec = 1,n_spec
               gs_avg(i_time,i_bin, i_spec) = 0d0
            enddo
            do i_loop = 1,n_loop
               g_avg(i_time, i_bin) = g_avg(i_time, i_bin)
     &              + bin_g(i_loop, i_time, i_bin)
               n_avg(i_time, i_bin) = n_avg(i_time, i_bin)
     &              + n(i_loop, i_time, i_bin)
               do i_spec = 1,n_spec
                  gs_avg(i_time, i_bin, i_spec) = 
     &                 gs_avg(i_time, i_bin, i_spec)
     &                 +bin_gs(i_loop,i_time,i_bin,i_spec)
               enddo
            enddo
            g_avg(i_time, i_bin) = g_avg(i_time, i_bin) / dble(n_loop)
            n_avg(i_time, i_bin) = n_avg(i_time, i_bin) / dble(n_loop)
            do i_spec=1,n_spec
               gs_avg(i_time,i_bin,i_spec) = 
     &              gs_avg(i_time,i_bin,i_spec) / dble(n_loop)
            enddo
         enddo
      enddo

      ! output raw number and mass data
      do i_time = 1,n_time
         write(f_out_num, '(//,a10,i10)') 'time', i_time - 1
         write(f_out_mass, '(//,a10,i10)') 'time', i_time - 1
         do i_bin = 1,n_bin
            write(f_out_num, '(i10,e20.10,'//n_loop_str//'e20.10)')
     &           i_bin, bin_r(i_bin),
     &           (n(i_loop, i_time, i_bin), i_loop = 1,n_loop)
            write(f_out_mass, '(i10,e20.10,'//n_loop_str//'e20.10)')
     &           i_bin, bin_r(i_bin),
     &           (bin_g(i_loop, i_time, i_bin), i_loop = 1,n_loop)
         enddo
         do i_spec = 1,n_spec
            write(f_out_mass,*)
            write(f_out_mass,*)
            do i_bin=1,n_bin
               write(f_out_mass, '(i10,e20.10,'//n_loop_str//'e20.10)')
     &              i_bin, bin_r(i_bin),
     &              (bin_gs(i_loop, i_time, i_bin,i_spec)
     &              ,i_loop =1,n_loop)
            enddo
         enddo
      enddo

      ! output averaged number and mass data
      do i_bin = 1,n_bin
         write(f_out_num_avg, '(i10,e20.10,'//n_time_str//'e20.10)')
     &        i_bin, bin_r(i_bin),
     &        (n_avg(i_time, i_bin), i_time = 1,n_time)
         write(f_out_mass_avg, '(i10,e20.10,'//n_time_str//'e20.10)')
     &        i_bin, bin_r(i_bin),
     &        (g_avg(i_time, i_bin), i_time = 1,n_time)
      enddo
      do i_spec = 1,n_spec
         write(f_out_mass_avg,*)
         write(f_out_mass_avg,*)
         do i_bin=1,n_bin
            write(f_out_mass_avg, '(i10,e20.10,'//n_time_str//'e20.10)')
     &           i_bin, bin_r(i_bin),
     &           (gs_avg(i_time, i_bin,i_spec),i_time =1,n_time)
         enddo
      enddo      

      ! output time, temperature and relative humidity data
      do i_time = 1,n_time
         write(f_out_temp, '(e20.10,'//n_loop_str//'e20.10)')
     &        time_avg(i_time),
     &        (temp(i_loop, i_time), i_loop = 1,n_loop)
         write(f_out_rh, '(e20.10,'//n_loop_str//'e20.10)')
     &        time_avg(i_time),
     &        (rh(i_loop, i_time), i_loop = 1,n_loop)
         write(f_out_time, '(i10,'//n_loop_str//'e20.10)')
     &        i_time,
     &        (time(i_loop, i_time), i_loop = 1,n_loop)
         write(f_out_temp_avg, '(e20.10,e20.10)')
     &        time_avg(i_time), temp_avg(i_time)
         write(f_out_rh_avg, '(e20.10,e20.10)')
     &        time_avg(i_time), rh_avg(i_time)
         write(f_out_time_avg, '(i10,e20.10)')
     &        i_time, time_avg(i_time)
      enddo
      
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
