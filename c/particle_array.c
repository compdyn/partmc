/* particle_array.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

/* Table of constant values */

double c_b2 = 2.;
int c__9 = 9;
int c__1 = 1;
int c__2 = 2;

/*     Utility functions for handling V array of particle volumes. */

/*     There are two different representations of particle size */
/*     distributions used throughout this code: a sectional */
/*     representation and an explicit particle representation. */

/*     The sectional representation stores the number and mass of */
/*     particles in bins, which are logarithmicly spaced. The bins are */
/*     described by the bin_v(n_bin) and bin_r(n_bin) arrays, which store the */
/*     volume and radius of the centerpoint of each bin. The variable */
/*     dlnr ... FIXME */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<       subroutine make_grid(n_bin, scal, rho_p, bin_v, bin_r, dlnr, e, r) >*/
/* Subroutine */ int make_grid__(int n_bin, int *scal, double *
	rho_p__, double *bin_v, double *bin_r, double dlnr, 
	double *e, double *r__)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1;

    /* Builtin functions */
    double log(double), pow_dd(double *, double *), pow_di(
	    double *, int *), exp(double);

    /* Local variables */
    int i__;
    double ax, emin;

/*<       int n_bin        ! INPUT: number of bins >*/
/*<       int scal         ! INPUT: scale factor >*/
/*<       real*8 rho_p         ! INPUT: density >*/
/*<       real*8 bin_v(n_bin)  ! OUTPUT: volume of particles in bins >*/
/*<       real*8 bin_r(n_bin)  ! OUTPUT: radius of particles in bins >*/
/*<       real*8 dlnr          ! OUTPUT: scale factor >*/
/*<       int i >*/
/*<       real*8 ax, e(n_bin), r(n_bin), emin >*/
/*<       real*8 pi >*/
/*<       parameter (pi = 3.14159265358979323846d0) >*/
/*<       dlnr = dlog(2d0) / (3d0 * scal) >*/
    /* Parameter adjustments */
    --r__;
    --e;
    --bin_r;
    --bin_v;

    /* Function Body */
    dlnr = log(2.) / (*scal * 3.);
/*<       ax = 2d0**(1d0 / scal) >*/
    d__1 = 1. / *scal;
    ax = pow_dd(&c_b2, &d__1);
/*<       emin = 1d-15 >*/
    emin = 1e-15;
/* FIXME: rewrite in a sane way */
/*<       do i = 1,n_bin >*/
    i__1 = n_bin;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* mass (FIXME: units?) */
/*<          e(i) = emin * 0.5d0 * (ax + 1d0) * ax**(i - 1) >*/
	i__2 = i__ - 1;
	e[i__] = emin * .5 * (ax + 1.) * pow_di(&ax, &i__2);
/* radius (um) */
/* FIXME: following line assumes rho_p = 1000 */
/*<          r(i) = 1000d0 * dexp(dlog(3d0 * e(i) / (4d0 * pi)) / 3d0) >*/
	r__[i__] = exp(log(e[i__] * 3. / 12.566370614359172) / 3.) * 1e3;
/* volume (FIXME: units?) */
/*<          bin_v(i) = 1d-6 * e(i) / rho_p >*/
	bin_v[i__] = e[i__] * 1e-6 / rho_p;
/* radius (m) */
/*<          bin_r(i) = 1d-6 * r(i) >*/
	bin_r[i__] = r__[i__] * 1e-6;
/*<       enddo >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* make_grid__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<       subroutine compute_volumes(n_bin, MM, n_ini, bin_r, dlnr, V, M) >*/
/* Subroutine */ int compute_volumes__(int n_bin, int *mm, 
	double *n_ini__, double *bin_r, double dlnr, 
	double *v, int *m)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1;

    /* Local variables */
    int i__, k, sum_a__, sum_e__, delta_n__;

/*<       int n_bin        ! INPUT: number of bins >*/
/*<       int MM           ! INPUT: physical size of V >*/
/*<       real*8 n_ini(n_bin)  ! INPUT: initial number distribution >*/
/*<       real*8 bin_r(n_bin)  ! INPUT: diameter of particles in bins >*/
/*<       real*8 dlnr          ! INPUT: scale factor >*/
/*<       real*8 V(MM)         ! OUTPUT: particle volumes >*/
/*<       int M            ! OUTPUT: logical dimension of V >*/
/*<       int k, i, sum_e, sum_a, delta_n >*/
/*<       real*8 pi >*/
/*<       parameter (pi = 3.14159265358979323846d0) >*/
/*<       sum_e = 0 >*/
    /* Parameter adjustments */
    --bin_r;
    --n_ini__;
    --v;

    /* Function Body */
    sum_e__ = 0;
/*<       do k = 1,n_bin >*/
    i__1 = n_bin;
    for (k = 1; k <= i__1; ++k) {
/*<          delta_n = int(n_ini(k) * dlnr) >*/
	delta_n__ = (int) (n_ini__[k] * dlnr);
/*<          sum_a = sum_e + 1 >*/
	sum_a__ = sum_e__ + 1;
/*<          sum_e = sum_e + delta_n >*/
	sum_e__ += delta_n__;
/*<          do i = sum_a,sum_e >*/
	i__2 = sum_e__;
	for (i__ = sum_a__; i__ <= i__2; ++i__) {
/*<             V(i) = 4d0/3d0 * pi * bin_r(k)**3 >*/
/* Computing 3rd power */
	    d__1 = bin_r[k];
	    v[i__] = d__1 * (d__1 * d__1) * 4.1887902047863905;
/*<          enddo >*/
	}
/*<       enddo >*/
    }
/*<       M = sum_e >*/
    *m = sum_e__;
/*<       return >*/
    return 0;
/*<       end >*/
} /* compute_volumes__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<       subroutine find_rand_pair(M, s1, s2) >*/
/* Subroutine */ int find_rand_pair__(int *m, int *s1, int *s2)
{
    extern double rand_(void);

/*<       int M       ! INPUT: number of particles >*/
/*<       int s1, s2  ! OUTPUT: s1 and s2 are not equal, random >*/
/*         particles with (1 <= s1,s2 <= M) */
/*<  100  s1 = int(rand() * M) + 1 >*/
L100:
    *s1 = (int) (rand_() * *m) + 1;
/*<       if ((s1 .lt. 1) .or. (s1 .gt. M)) goto 100 >*/
    if (*s1 < 1 || *s1 > *m) {
	goto L100;
    }
/*<  101  s2 = int(rand() * M) + 1 >*/
L101:
    *s2 = (int) (rand_() * *m) + 1;
/*<       if ((s2 .lt. 1) .or. (s2 .gt. M)) goto 101 >*/
    if (*s2 < 1 || *s2 > *m) {
	goto L101;
    }
/*<       if (s1 .eq. s2) goto 101 >*/
    if (*s1 == *s2) {
	goto L101;
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* find_rand_pair__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<    >*/
/* Subroutine */ int find_rand_pair_acc_rej__(int *mm, int *m, 
	double *v, double *max_k__, S_fp kernel, int *s1, int 
	*s2)
{
    double k, p;
    extern double rand_(void);
    extern /* Subroutine */ int find_rand_pair__(int *, int *, 
	    int *);

/*<       int MM      ! INPUT: physical dimension of V >*/
/*<       int M       ! INPUT: logical dimension of V >*/
/*<       real*8 V(MM)    ! INPUT: array of particle volumes >*/
/*<       real*8 max_k    ! INPUT: maximum value of the kernel >*/
/*<       external kernel ! INPUT: kernel function >*/
/*<       int s1, s2  ! OUTPUT: s1 and s2 are not equal, random >*/
/*         particles with V(s1/s2) != 0 */
/*<       real*8 k, p >*/
/*<  200  continue >*/
    /* Parameter adjustments */
    --v;

    /* Function Body */
L200:
/*<       call find_rand_pair(M, s1, s2) ! test particles s1, s2 >*/
    find_rand_pair__(m, s1, s2);
/*<       call kernel(V(s1), V(s2), k) >*/
    (*kernel)(&v[*s1], &v[*s2], &k);
/*<       p = k / max_k     ! collision probability    >*/
    p = k / *max_k__;
/*<       if (dble(rand()) .gt. p ) goto 200 >*/
    if ((double) rand_() > p) {
	goto L200;
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* find_rand_pair_acc_rej__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<    >*/
/* Subroutine */ int coagulate_(int *mm, int *m, double *v, 
	double v_comp, int n_bin, double *bin_v, 
	double *bin_r, double *bin_g, int *bin_n, 
	double dlnr, int *s1, int *s2, logical *bin_change__)
{
    /* Builtin functions */
    int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    int k1, k2, kn;
    extern /* Subroutine */ int exit_(int *), particle_in_bin__(
	    double *, int *, double *, int *);

    /* Fortran I/O blocks */
    cilist io___13 = { 0, 6, 0, 0, 0 };


/*<       int MM           ! INPUT: physical dimension of V >*/
/*<       int M            ! INPUT/OUTPUT: logical dimension of V >*/
/*<       real*8 V(MM)         ! INPUT/OUTPUT: particle volumes >*/
/*<       real*8 V_comp        ! INPUT: computational volume >*/
/*<       int n_bin        ! INPUT: number of bins >*/
/*<       real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins >*/
/*<       real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins >*/
/*<       real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins >*/
/*<       int bin_n(n_bin) ! INPUT/OUTPUT: number in bins >*/
/*<       real*8 dlnr          ! INPUT: bin scale factor >*/
/*<       int s1           ! INPUT: first particle to coagulate >*/
/*<       int s2           ! INPUT: second particle to coagulate >*/
/*<       logical bin_change   ! OUTPUT: whether an empty bin filled, >*/
/*         or a filled bin became empty */
/*<       int k1, k2, kn >*/
/*<       bin_change = .false. >*/
    /* Parameter adjustments */
    --v;
    --bin_n;
    --bin_g;
    --bin_r;
    --bin_v;

    /* Function Body */
    *bin_change__ = FALSE_;
/* remove s1 and s2 from bins */
/*<       call particle_in_bin(V(s1), n_bin, bin_v, k1) >*/
    particle_in_bin__(&v[*s1], n_bin__, &bin_v[1], &k1);
/*<       call particle_in_bin(V(s2), n_bin, bin_v, k2) >*/
    particle_in_bin__(&v[*s2], n_bin__, &bin_v[1], &k2);
/*<       bin_n(k1) = bin_n(k1) - 1 >*/
    --bin_n[k1];
/*<       bin_n(k2) = bin_n(k2) - 1 >*/
    --bin_n[k2];
/* DEBUG */
/*      write(*,*)'V(s1),V(s2),k1,k2,n(k1),n(k2) = ', */
/*     &     V(s1), V(s2), k1, k2, bin_n(k1), bin_n(k2) */
/*      write(*,*)'r(k1),r(k2) = ', bin_r(k1), bin_r(k2) */
/*      if (bin_n(k1) .eq. 0) write(*,*)'bin became empty:', k1 */
/*      if (bin_n(k2) .eq. 0) write(*,*)'bin became empty:', k2 */
/* DEBUG */
/*<       bin_g(k1) = bin_g(k1) - V(s1) >*/
    bin_g[k1] -= v[*s1];
/*<       bin_g(k2) = bin_g(k2) - V(s2) >*/
    bin_g[k2] -= v[*s2];
/*<       if ((bin_n(k1) .lt. 0) .or. (bin_n(k2) .lt. 0)) then >*/
    if (bin_n[k1] < 0 || bin_n[k2] < 0) {
/*<          write(*,*)'ERROR: invalid bin_n' >*/
	s_wsle(&io___13);
	do_lio(&c__9, &c__1, "ERROR: invalid bin_n", (ftnlen)20);
	e_wsle();
/*<          call exit(2) >*/
	exit_(&c__2);
/*<       endif >*/
    }
/*<       V(s1) = V(s1) + V(s2) ! add particle 2 onto particle 1 >*/
    v[*s1] += v[*s2];
/*<       V(s2) = V(M)          ! shift the last particle into empty slot >*/
    v[*s2] = v[*m];
/*<       M = M - 1             ! shorten array >*/
    --(*m);
/* add new particle to bins */
/*<       call particle_in_bin(V(s1), n_bin, bin_v, kn) >*/
    particle_in_bin__(&v[*s1], n_bin__, &bin_v[1], &kn);
/* DEBUG */
/*      write(*,*)'V(s1),kn,n(kn) = ', V(s1), kn, bin_n(kn) */
/*      write(*,*)'r(kn) = ', bin_r(kn) */
/*      if (bin_n(kn) .le. 0) write(*,*)'bin became full:', kn */
/* DEBUG */
/*<       bin_n(kn) = bin_n(kn) + 1 >*/
    ++bin_n[kn];
/*<       bin_g(kn) = bin_g(kn) + V(s1) >*/
    bin_g[kn] += v[*s1];
/*<    >*/
    if (bin_n[k1] == 0 || bin_n[k2] == 0) {
	*bin_change__ = TRUE_;
    }
/*<    >*/
    if (bin_n[kn] == 1 && kn != k1 && kn != k2) {
	*bin_change__ = TRUE_;
    }
/* DEBUG */
/*      if ((bin_n(k1) .eq. 0) .or. (bin_n(k2) .eq. 0)) then */
/*         write(*,*)'k1,k2,n(k1),n(k2) = ', */
/*     &        k1, k2, bin_n(k1), bin_n(k2) */
/*         write(*,*)'kn,n(kn) = ', kn, bin_n(kn) */
/*       endif */
/*      if ((bin_n(kn) .eq. 1) .and. (kn .ne. k1) .and. (kn .ne. k2)) then */
/*         write(*,*)'k1,k2,n(k1),n(k2) = ', */
/*     &        k1, k2, bin_n(k1), bin_n(k2) */
/*         write(*,*)'kn,n(kn) = ', kn, bin_n(kn) */
/*       endif */
/* DEBUG */
/*      if (bin_change) write(*,*)'bin_change is true' */
/*<       return >*/
    return 0;
/*<       end >*/
} /* coagulate_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<    >*/
/* Subroutine */ int maybe_coag_pair__(int *mm, int *m, double *v,
	 double v_comp, int n_bin, double *bin_v, 
	double *bin_r, double *bin_g, int *bin_n, 
	double dlnr, double *del_t__, int *n_samp__, S_fp kernel,
	 logical *did_coag__, logical *bin_change__)
{
    /* Builtin functions */
    double exp(double);

    /* Local variables */
    double k, p;
    extern /* Subroutine */ int coagulate_(int *, int *, double *,
	     double *, int *, double *, double *, double *
	    , int *, double *, int *, int *, logical *);
    int s1, s2;
    extern double rand_(void);
    double expo;
    extern /* Subroutine */ int find_rand_pair__(int *, int *, 
	    int *);

/*<       int MM           ! INPUT: physical dimension of V >*/
/*<       int M            ! INPUT/OUTPUT: logical dimension of V >*/
/*<       real*8 V(MM)         ! INPUT/OUTPUT: particle volumes >*/
/*<       real*8 V_comp        ! INPUT: computational volume >*/
/*<       int n_bin        ! INPUT: number of bins >*/
/*<       real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins >*/
/*<       real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins >*/
/*<       real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins >*/
/*<       int bin_n(n_bin) ! INPUT/OUTPUT: number in bins >*/
/*<       real*8 dlnr          ! INPUT: bin scale factor >*/
/*<       real*8 del_t         ! INPUT: timestep >*/
/*<       int n_samp       ! INPUT: number of samples per timestep >*/
/*<       external kernel      ! INPUT: kernel function >*/
/*<       logical did_coag     ! OUTPUT: whether a coagulation occured >*/
/*<       logical bin_change   ! OUTPUT: whether bin structure changed >*/
/*<       int s1, s2 >*/
/*<       real*8 expo, p, k >*/
/*<       call find_rand_pair(M, s1, s2) ! test particles s1, s2 >*/
    /* Parameter adjustments */
    --v;
    --bin_n;
    --bin_g;
    --bin_r;
    --bin_v;

    /* Function Body */
    find_rand_pair__(m, &s1, &s2);
/*<       call kernel(V(s1), V(s2), k) >*/
    (*kernel)(&v[s1], &v[s2], &k);
/*<       expo = k * 1d0/V_comp * del_t * (M*(M-1)/2d0) / n_samp >*/
    expo = k * 1. / v_comp * *del_t__ * (*m * (*m - 1) / 2.) / *n_samp__;
/*<       p = 1d0 - exp(-expo) ! probability of coagulation >*/
    p = 1. - exp(-expo);
/*<       bin_change = .false. >*/
    *bin_change__ = FALSE_;
/*<       if (dble(rand()) .lt. p) then >*/
    if ((double) rand_() < p) {
/*<    >*/
	coagulate_(mm, m, &v[1], v_comp__, n_bin__, &bin_v[1], &bin_r[1], 
		&bin_g[1], &bin_n[1], dlnr, &s1, &s2, bin_change__);
/*<          did_coag = .true. >*/
	*did_coag__ = TRUE_;
/*<       else >*/
    } else {
/*<          did_coag = .false. >*/
	*did_coag__ = FALSE_;
/*<       endif >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* maybe_coag_pair__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<       subroutine kernel_avg(MM, M, V, kernel, n_samp, k_avg) >*/
/* Subroutine */ int kernel_avg__(int *mm, int *m, double *v, 
	S_fp kernel, int *n_samp__, double *k_avg__)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    int i__;
    double k;
    int s1, s2;
    double k_sum__;
    extern /* Subroutine */ int find_rand_pair__(int *, int *, 
	    int *);

/* FIXME: use binned data instead */
/*<       int MM      ! INPUT: physical dimension of V >*/
/*<       int M       ! INPUT: logical dimension of V >*/
/*<       real*8 V(MM)    ! INPUT: array of particle volumes >*/
/*<       external kernel ! INPUT: kernel function >*/
/*<       int n_samp  ! INPUT: number of samples to use (squared) >*/
/*<       real*8 k_avg    ! OUTPUT: estimated average of kernel values >*/
/*<       int i, s1, s2 >*/
/*<       real*8 k, k_sum >*/
/*<       k_sum = 0d0 >*/
    /* Parameter adjustments */
    --v;

    /* Function Body */
    k_sum__ = 0.;
/*<       do i = 1,(n_samp**2) >*/
/* Computing 2nd power */
    i__2 = *n_samp__;
    i__1 = i__2 * i__2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          call find_rand_pair(M, s1, s2) >*/
	find_rand_pair__(m, &s1, &s2);
/*<          call kernel(V(s1), V(s2), k) >*/
	(*kernel)(&v[s1], &v[s2], &k);
/*<          k_sum = k_sum + k >*/
	k_sum__ += k;
/*<       enddo >*/
    }
/*<       k_avg  = k_sum / (n_samp**2) >*/
/* Computing 2nd power */
    i__1 = *n_samp__;
    *k_avg__ = k_sum__ / (i__1 * i__1);
/*<       return >*/
    return 0;
/*<       end >*/
} /* kernel_avg__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<    >*/
/* Subroutine */ int double_(int *mm, int *m, double *v, 
	double v_comp, int n_bin, double *bin_v, 
	double *bin_r, double *bin_g, int *bin_n, 
	double dlnr)
{
    /* System generated locals */
    int i__1;

    /* Builtin functions */
    int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    int i__;
    extern /* Subroutine */ int exit_(int *);

    /* Fortran I/O blocks */
    cilist io___25 = { 0, 6, 0, 0, 0 };


/*<       int MM           ! INPUT: physical dimension of V >*/
/*<       int M            ! INPUT/OUTPUT: logical dimension of V >*/
/*<       real*8 V(MM)         ! INPUT/OUTPUT: particle volumes >*/
/*<       real*8 V_comp        ! INPUT/OUTPUT: computational volume >*/
/*<       int n_bin        ! INPUT: number of bins >*/
/*<       real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins >*/
/*<       real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins >*/
/*<       real*8 bin_g(n_bin)  ! INPUT/OUTPUT: mass in bins >*/
/*<       int bin_n(n_bin) ! INPUT/OUTPUT: number in bins >*/
/*<       real*8 dlnr          ! INPUT: bin scale factor >*/
/*<       int i >*/
/* only double if we have enough space to do so */
/*<       if (M .gt. MM / 2) then >*/
    /* Parameter adjustments */
    --v;
    --bin_n;
    --bin_g;
    --bin_r;
    --bin_v;

    /* Function Body */
    if (*m > *mm / 2) {
/*<          write(*,*)'ERROR: double without enough space' >*/
	s_wsle(&io___25);
	do_lio(&c__9, &c__1, "ERROR: double without enough space", (ftnlen)34)
		;
	e_wsle();
/*<          call exit(2) >*/
	exit_(&c__2);
/*<       endif >*/
    }
/* double V and associated structures */
/*<       do i = 1,M >*/
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          V(i + M) = V(i) >*/
	v[i__ + *m] = v[i__];
/*<       enddo >*/
    }
/*<       M = 2 * M >*/
    *m <<= 1;
/*<       V_comp = 2d0 * V_comp >*/
    v_comp *= 2.;
/* double bin structures */
/*<       do i = 1,n_bin >*/
    i__1 = n_bin;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          bin_g(i) = bin_g(i) * 2d0 >*/
	bin_g[i__] *= 2.;
/*<          bin_n(i) = bin_n(i) * 2 >*/
	bin_n[i__] <<= 1;
/*<       enddo >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* double_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<       subroutine est_k_max(n_bin, bin_v, bin_n, kernel, k_max, use_bin) >*/
/* Subroutine */ int est_k_max__(int n_bin, double *bin_v, 
	int *bin_n, S_fp kernel, double *k_max__, logical *
	use_bin__)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    int i__, j;
    double k;

/*<       int n_bin         ! INPUT: number of bins >*/
/*<       real*8 bin_v(n_bin)   ! INPUT: volume of particles in bins >*/
/*<       int bin_n(n_bin)  ! INPUT: number in each bin >*/
/*<       external kernel       ! INPUT: kernel function >*/
/*<       real*8 k_max          ! OUTPUT: maximum kernel value >*/
/*<       real*8 k >*/
/*<       int i, j >*/
/*<       logical use_bin(n_bin) >*/
/*<       real*8 pi >*/
/*<       parameter (pi = 3.14159265358979323846d0) >*/
/*      write(*,*)'RECALCULATING K_MAX' */
/* use_bin starts as non-empty bins */
/*<       do i = 1,n_bin >*/
    /* Parameter adjustments */
    --use_bin__;
    --bin_n;
    --bin_v;

    /* Function Body */
    i__1 = n_bin;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          use_bin(i) = (bin_n(i) .gt. 0) >*/
	use_bin__[i__] = bin_n[i__] > 0;
/*<       enddo >*/
    }
/* add all bins downstream of non-empty bins */
/*<       do i = 2,n_bin >*/
    i__1 = n_bin;
    for (i__ = 2; i__ <= i__1; ++i__) {
/*<          if (use_bin(i)) use_bin(i-1) = .true. >*/
	if (use_bin__[i__]) {
	    use_bin__[i__ - 1] = TRUE_;
	}
/*<       enddo >*/
    }
/* add all bins upstream of non-empty bins */
/*<       do i = (n_bin-1),1,-1 >*/
    for (i__ = n_bin - 1; i__ >= 1; --i__) {
/*<          if (use_bin(i)) use_bin(i+1) = .true. >*/
	if (use_bin__[i__]) {
	    use_bin__[i__ + 1] = TRUE_;
	}
/*<       enddo >*/
    }
/*<       k_max = 0d0 >*/
    *k_max__ = 0.;
/*<       do i = 1,n_bin >*/
    i__1 = n_bin;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          if (use_bin(i)) then >*/
	if (use_bin__[i__]) {
/*<             do j = 1,i >*/
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
/*<                if (use_bin(j)) then >*/
		if (use_bin__[j]) {
/*<                   call kernel(bin_v(i), bin_v(j), k) >*/
		    (*kernel)(&bin_v[i__], &bin_v[j], &k);
/*<                   if (k .gt. k_max) then >*/
		    if (k > *k_max__) {
/*<                      k_max = k >*/
			*k_max__ = k;
/*<                   endif >*/
		    }
/*<                endif >*/
		}
/*<             enddo >*/
	    }
/*<          endif >*/
	}
/*<       enddo >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* est_k_max__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<       subroutine particle_in_bin(v, n_bin, bin_v, k) >*/
/* Subroutine */ int particle_in_bin__(double *v, int n_bin, 
	double *bin_v, int *k)
{
/* FIXME: for log-spaced bins we can do this without search */
/*<       real*8 v             ! INPUT: volume of particle >*/
/*<       int n_bin        ! INPUT: number of bins >*/
/*<       real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins >*/
/*<       int k            ! OUTPUT: bin number containing particle >*/
/*<       k = 0 >*/
    /* Parameter adjustments */
    --bin_v;

    /* Function Body */
    *k = 0;
/*<  300  k = k + 1 >*/
L300:
    ++(*k);
/*      write(*,*)'k,bin_v(k) = ', k, bin_v(k) */
/*<    >*/
    if (*k < n_bin && *v > (bin_v[*k] + bin_v[*k + 1]) / 2.) {
	goto L300;
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* particle_in_bin__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<    >*/
/* Subroutine */ int moments_(int *mm, int *m, double *v, 
	double v_comp, int n_bin, double *bin_v, 
	double *bin_r, double *bin_g, int *bin_n, 
	double dlnr)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    int i__, k;
    extern /* Subroutine */ int particle_in_bin__(double *, int *, 
	    double *, int *);

/*<       int MM           ! INPUT: physical dimension of V >*/
/*<       int M            ! INPUT: logical dimension of V >*/
/*<       real*8 V(MM)         ! INPUT: particle volumes >*/
/*<       real*8 V_comp        ! INPUT: computational volume >*/
/*<       int n_bin        ! INPUT: number of bins >*/
/*<       real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins >*/
/*<       real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins >*/
/*<       real*8 bin_g(n_bin)  ! OUTPUT: mass in bins >*/
/*<       int bin_n(n_bin) ! OUTPUT: number in bins >*/
/*<       real*8 dlnr          ! INPUT: bin scale factor >*/
/*<       int i, k >*/
/*<       do k = 1,n_bin >*/
    /* Parameter adjustments */
    --v;
    --bin_n;
    --bin_g;
    --bin_r;
    --bin_v;

    /* Function Body */
    i__1 = n_bin;
    for (k = 1; k <= i__1; ++k) {
/*<          bin_g(k) = 0d0 >*/
	bin_g[k] = 0.;
/*<          bin_n(k) = 0 >*/
	bin_n[k] = 0;
/*<       enddo >*/
    }
/*<       do i = 1,M >*/
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          call particle_in_bin(V(i), n_bin, bin_v, k) >*/
	particle_in_bin__(&v[i__], n_bin__, &bin_v[1], &k);
/*<          bin_g(k) = bin_g(k) + V(i) >*/
	bin_g[k] += v[i__];
/*<          bin_n(k) = bin_n(k) + 1 >*/
	++bin_n[k];
/*<       enddo >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* moments_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<       subroutine check_event(time, interval, last_time, do_event) >*/
/* Subroutine */ int check_event__(double time, double *interval, 
	double *last_time__, logical *do_event__)
{
    /* System generated locals */
    double d__1;

    /* Builtin functions */
    double d_int(double *);

    /* Local variables */
    double interval_below__;

/*<       real*8 time       ! INPUT: cubin_rent time >*/
/*<       real*8 interval   ! INPUT: how often the event should be done >*/
/*<       real*8 last_time  ! INPUT/OUTPUT: when the event was last done >*/
/*<       logical do_event  ! OUTPUT: whether the event should be done >*/
/*<       real*8 interval_below >*/
/*<       if (time .eq. 0d0) then >*/
    if (time == 0.) {
/*<          last_time = 0d0 >*/
	*last_time__ = 0.;
/*<          do_event = .true. >*/
	*do_event__ = TRUE_;
/*<       else >*/
    } else {
/*<          interval_below = aint(time / interval) * interval >*/
	d__1 = time / *interval;
	interval_below__ = d_int(&d__1) * *interval;
/*<          if (last_time .lt. interval_below) then >*/
	if (*last_time__ < interval_below__) {
/*<             last_time = time >*/
	    *last_time__ = time;
/*<             do_event = .true. >*/
	    *do_event__ = TRUE_;
/*<          else >*/
	} else {
/*<             do_event = .false. >*/
	    *do_event__ = FALSE_;
/*<          endif >*/
	}
/*<       endif >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* check_event__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<       subroutine print_header(n_loop, n_bin, n_time) >*/
/* Subroutine */ int print_header__(int *n_loop__, int n_bin, 
	int *n_time__)
{
    /* Builtin functions */
    int s_wsfe(cilist *), do_fio(int *, char *, ftnlen), e_wsfe(void);

    /* Fortran I/O blocks */
    cilist io___33 = { 0, 30, 0, "(a10,i10)", 0 };
    cilist io___34 = { 0, 30, 0, "(a10,i10)", 0 };
    cilist io___35 = { 0, 30, 0, "(a10,i10)", 0 };


/*<       int n_loop  ! INPUT: number of loops >*/
/*<       int n_bin   ! INPUT: number of bins >*/
/*<       int n_time  ! INPUT: number of times >*/
/*<       write(30,'(a10,i10)') 'n_loop', n_loop >*/
    s_wsfe(&io___33);
    do_fio(&c__1, "n_loop", (ftnlen)6);
    do_fio(&c__1, (char *)&(*n_loop__), (ftnlen)sizeof(int));
    e_wsfe();
/*<       write(30,'(a10,i10)') 'n_bin', n_bin >*/
    s_wsfe(&io___34);
    do_fio(&c__1, "n_bin", (ftnlen)5);
    do_fio(&c__1, (char *)&(n_bin), (ftnlen)sizeof(int));
    e_wsfe();
/*<       write(30,'(a10,i10)') 'n_time', n_time >*/
    s_wsfe(&io___35);
    do_fio(&c__1, "n_time", (ftnlen)6);
    do_fio(&c__1, (char *)&(*n_time__), (ftnlen)sizeof(int));
    e_wsfe();
/*<       return >*/
    return 0;
/*<       end >*/
} /* print_header__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<    >*/
/* Subroutine */ int print_info__(double time, double v_comp, 
	int n_bin, double *bin_v, double *bin_r, 
	double *bin_g, int *bin_n, double dlnr)
{
    /* System generated locals */
    int i__1;
    double d__1, d__2;

    /* Builtin functions */
    int s_wsfe(cilist *), do_fio(int *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    int k;

    /* Fortran I/O blocks */
    cilist io___36 = { 0, 30, 0, "(a10,e14.5)", 0 };
    cilist io___38 = { 0, 30, 0, "(i8,3e14.5)", 0 };


/*<       real*8 time          ! INPUT: cubin_rent simulation time >*/
/*<       real*8 V_comp        ! INPUT: computational volume >*/
/*<       int n_bin        ! INPUT: number of bins >*/
/*<       real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins >*/
/*<       real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins >*/
/*<       real*8 bin_g(n_bin)  ! INPUT: mass in bins >*/
/*<       int bin_n(n_bin) ! INPUT: number in bins >*/
/*<       real*8 dlnr          ! INPUT: bin scale factor >*/
/*<       int k >*/
/*<       write(30,'(a10,e14.5)') 'time', time >*/
    /* Parameter adjustments */
    --bin_n;
    --bin_g;
    --bin_r;
    --bin_v;

    /* Function Body */
    s_wsfe(&io___36);
    do_fio(&c__1, "time", (ftnlen)4);
    do_fio(&c__1, (char *)&(time), (ftnlen)sizeof(double));
    e_wsfe();
/*<       do k = 1,n_bin >*/
    i__1 = n_bin;
    for (k = 1; k <= i__1; ++k) {
/*<    >*/
	s_wsfe(&io___38);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(int));
	do_fio(&c__1, (char *)&bin_r[k], (ftnlen)sizeof(double));
	d__1 = bin_n[k] / v_comp / dlnr;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(double));
	d__2 = bin_g[k] / v_comp / dlnr;
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(double));
	e_wsfe();
/*<       enddo >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* print_info__ */

