
program binomial_sample

  use pmc_rand
  use pmc_util

  integer :: n, n_samp, i_samp, k
  real(kind=dp) :: p
  real(kind=dp), allocatable :: count_pdf(:), pdf(:)
  integer, allocatable :: count(:)
  character(len=1000) :: tmp

  ! process commandline arguments
  if (command_argument_count() .ne. 3) then
     write(6,*) 'Usage: binomial_sample <n> <p> <n_samp>'
     stop 2
  endif
  call get_command_argument(1, tmp)
  n = string_to_integer(tmp)
  call get_command_argument(2, tmp)
  p = string_to_real(tmp)
  call get_command_argument(3, tmp)
  n_samp = string_to_integer(tmp)

  allocate(count(0:n), count_pdf(0:n), pdf(0:n))

  ! initialize the RNG with a random seed
  call pmc_srand(0, 0)
  
  ! compute exact PDF
  if (n_samp == 0) then
     do k = 0,n
        pdf(k) = p**k * (1d0 - p)**(n - k) &
             * gamma(real(n + 1, kind=dp)) &
             / (gamma(real(k + 1, kind=dp)) &
             * gamma(real(n - k + 1, kind=dp)))
     end do
  end if

  ! compute sampled PDF
  if (n_samp > 0) then
     count = 0
     do i_samp = 1,n_samp
        k = rand_binomial(n, p)
        count(k) = count(k) + 1
     end do
     count_pdf = dble(count) / dble(n_samp)
  end if

  ! write results
  do k = 0,n
     if (n_samp == 0) then
        write(*,'(i6,e20.10)') k, pdf(k)
     else
        write(*,'(i6,e20.10)') k, count_pdf(k)
     end if
  end do

  ! cleanup the RNG
  call pmc_rand_finalize()

end program binomial_sample
