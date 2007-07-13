
program poisson_sample

  use mod_rand_poisson
  use mod_util

  integer :: k_max, n_samp, i_samp, k, i
  real*8 :: lambda
  real*8, allocatable :: count_pdf(:), pdf(:)
  integer, allocatable :: count(:)
  character(len=1000) :: tmp

  ! process commandline arguments
  if (iargc() .ne. 3) then
     write(6,*) 'Usage: process_summary <lambda> <k_max> <n_samp>'
     call exit(2)
  endif
  call getarg(1, tmp)
  lambda = string_to_real(tmp)
  call getarg(2, tmp)
  k_max = string_to_integer(tmp)
  call getarg(3, tmp)
  n_samp = string_to_integer(tmp)

  allocate(count(0:k_max), count_pdf(0:k_max), pdf(0:k_max))
  
  ! compute exact PDF
  do k = 0,k_max
     pdf(k) = exp(-lambda) * lambda**k
     do i = 1,k
        pdf(k) = pdf(k) / dble(i)
     end do
  end do

  ! compute sampled PDF
  count = 0
  do i_samp = 1,n_samp
     k = rand_poisson(lambda)
     if (k <= k_max) count(k) = count(k) + 1
  end do
  count_pdf = dble(count) / dble(n_samp)

  ! write results
  do k = 0,k_max
     write(*,'(i6,e20.10,e20.10)') k, pdf(k), count_pdf(k)
  end do

end program poisson_sample
