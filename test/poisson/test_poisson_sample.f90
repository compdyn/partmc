
program poisson_sample

  use pmc_rand
  use pmc_util

  integer :: k_max, n_samp, i_samp, k
  real(kind=dp) :: lambda
  real(kind=dp), allocatable :: count_pdf(:), pdf(:)
  integer, allocatable :: count(:)
  character(len=1000) :: tmp

  ! process commandline arguments
  if (command_argument_count() .ne. 3) then
     write(6,*) 'Usage: poisson_sample <lambda> <k_max> <n_samp>'
     stop 2
  endif
  call get_command_argument(1, tmp)
  lambda = string_to_real(tmp)
  call get_command_argument(2, tmp)
  k_max = string_to_integer(tmp)
  call get_command_argument(3, tmp)
  n_samp = string_to_integer(tmp)

  allocate(count(0:k_max), count_pdf(0:k_max), pdf(0:k_max))

  ! initialize the RNG with a random seed
  call pmc_srand(0, 0)
  
  ! compute exact PDF
  if (n_samp == 0) then
     do k = 0,k_max
        pdf(k) = exp(-lambda) * lambda**k / gamma(real(k + 1, kind=dp))
     end do
  end if

  ! compute sampled PDF
  if (n_samp > 0) then
     count = 0
     do i_samp = 1,n_samp
        k = rand_poisson(lambda)
        if (k <= k_max) count(k) = count(k) + 1
     end do
     count_pdf = dble(count) / dble(n_samp)
  end if

  ! write results
  do k = 0,k_max
     if (n_samp == 0) then
        write(*,'(i6,e20.10)') k, pdf(k)
     else
        write(*,'(i6,e20.10)') k, count_pdf(k)
     end if
  end do

  ! cleanup the RNG
  call pmc_rand_finalize()

end program poisson_sample
