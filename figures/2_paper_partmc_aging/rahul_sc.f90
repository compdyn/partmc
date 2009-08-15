	program sc
      implicit none

      real Dd, Dd_micron, k
      real sigma_s, Temp, Mw, rhow, Ru, A
	real Sca, Scn
      integer iter

      write(6,*)'Enter Ddry in micron and avg kappa'
      read(5,*)Dd_micron, k

! find dry and wet diameters
	  Dd = Dd_micron * 1.e-4	! cm

! calculate critical supersaturation
        sigma_s = 72	  ! dynes/cm
        Temp    = 298.15  ! K
	  Mw      = 18.0	  ! g/mol
	  rhow    = 1.0 	  ! g/cc
	  Ru      = 8.3144e7! erg/K/mol
        A       = (4.0*sigma_s*Mw)/(Ru*Temp*rhow)

        call critical_supersat(Dd, k, A, Sca, Scn, iter)

	  write(6,1010)Dd_micron, k, Sca, Scn, iter
1010	  format('Ddry = ',f5.2,' um', 2x, 'kappa =', f6.3, 2x, 'Sc-analytical = ', f7.4, 2x, 'Sc-numerical = ', f7.4,2x, 'iter = ', i5 )



	end



!------------------------------------------------------------------------------------------
      subroutine critical_supersat(Dd, k, A, Sca, Scn, iter)
      implicit none
      real Dc, Dd, k, A, Sca, Scn
      real y1, y2, ymid, yc, fd1, fd2, fmid, relerr
      real fd
      integer iter

! analytical solution for k >= 0.2
      if(k .gt. 0.0)then
        Sca = 100.*(exp(sqrt((4.0*A**3.0)/(27.0*Dd**3.0 * k))) - 1.0)	! % critical supersaturation
      elseif(k .eq. 0.0)then
        Sca = 100.*(exp(A/Dd) - 1.0)
      endif

! numerical bisection solution (y = D^3 - Dd^3)
      y1 = 0.0
      y2 = 2.*Dd**3.
      iter = 0

10    fd1 = fd(y1, Dd, k, A)
      fd2 = fd(y2, Dd, k, A)
      iter = iter + 1
      if(fd1 .gt. 0.0 .and. fd2 .gt. 0.0)then
        y1 = y2
        y2 = 2.*y2
        goto 10
      endif
      
! the two ends are set, now begin bisecting
20    ymid = 0.5*(y1 + y2)
      fmid = fd(ymid, Dd, k, A)
      iter = iter + 1
      if(fmid .lt. 0.0)then
        y2 = ymid
      else
        y1 = ymid
      endif

      relerr = (y2 - y1)/y2
      if(relerr .gt. 1.e-5)goto 20
   
! bisect one last time to get the critical values
      yc = 0.5*(y1 + y2)
      Dc = (yc + Dd**3.0)**0.33333333333
      Scn = 100.*(yc/(yc + k*Dd**3.)*exp(A/Dc) - 1.0)		! % critical supersaturation

      return
      end subroutine critical_supersat



!-----------------------------------------------------------------------
      function fd(y, Dd, k, A)
      implicit none
      real Dd, D, y, k, A, dum
      real fd

      D   = (y + Dd**3.0)**0.33333333333
      dum = y + k*Dd**3.

      fd = -3.*exp(A/D)*D*D*y/dum**2. + 3.*exp(A/D)*D*D/dum - A*exp(A/D)*y/(D**2. * dum)

      return
      end function fd


