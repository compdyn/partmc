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

#include "f2c.h"

/* Table of constant values */

static doublereal c_b2 = 2.;
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__2 = 2;

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
/* Subroutine */ int make_grid__(integer *n_bin__, integer *scal, doublereal *
	rho_p__, doublereal *bin_v__, doublereal *bin_r__, doublereal *dlnr, 
	doublereal *e, doublereal *r__)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), pow_dd(doublereal *, doublereal *), pow_di(
	    doublereal *, integer *), exp(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal ax, emin;

/* INPUT: number of bins */
/* INPUT: scale factor */
/* INPUT: density */
/* OUTPUT: volume of particles in bins */
/* OUTPUT: radius of particles in bins */
/* OUTPUT: scale factor */
    /* Parameter adjustments */
    --r__;
    --e;
    --bin_r__;
    --bin_v__;

    /* Function Body */
    *dlnr = log(2.) / (*scal * 3.);
    d__1 = 1. / *scal;
    ax = pow_dd(&c_b2, &d__1);
    emin = 1e-15;
/* FIXME: rewrite in a sane way */
    i__1 = *n_bin__;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* mass (FIXME: units?) */
	i__2 = i__ - 1;
	e[i__] = emin * .5 * (ax + 1.) * pow_di(&ax, &i__2);
/* radius (um) */
/* FIXME: following line assumes rho_p = 1000 */
	r__[i__] = exp(log(e[i__] * 3. / 12.566370614359172) / 3.) * 1e3;
/* volume (FIXME: units?) */
	bin_v__[i__] = e[i__] * 1e-6 / *rho_p__;
/* radius (m) */
	bin_r__[i__] = r__[i__] * 1e-6;
    }
    return 0;
} /* make_grid__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int compute_volumes__(integer *n_bin__, integer *mm, 
	doublereal *n_ini__, doublereal *bin_r__, doublereal *dlnr, 
	doublereal *v, integer *m)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, k, sum_a__, sum_e__, delta_n__;

/* INPUT: number of bins */
/* INPUT: physical size of V */
/* INPUT: initial number distribution */
/* INPUT: diameter of particles in bins */
/* INPUT: scale factor */
/* OUTPUT: particle volumes */
/* OUTPUT: logical dimension of V */
    /* Parameter adjustments */
    --bin_r__;
    --n_ini__;
    --v;

    /* Function Body */
    sum_e__ = 0;
    i__1 = *n_bin__;
    for (k = 1; k <= i__1; ++k) {
	delta_n__ = (integer) (n_ini__[k] * *dlnr);
	sum_a__ = sum_e__ + 1;
	sum_e__ += delta_n__;
	i__2 = sum_e__;
	for (i__ = sum_a__; i__ <= i__2; ++i__) {
/* Computing 3rd power */
	    d__1 = bin_r__[k];
	    v[i__] = d__1 * (d__1 * d__1) * 4.1887902047863905;
	}
    }
    *m = sum_e__;
    return 0;
} /* compute_volumes__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int find_rand_pair__(integer *m, integer *s1, integer *s2)
{
    extern doublereal rand_(void);

/* INPUT: number of particles */
/*         particles with (1 <= s1,s2 <= M) */
/* OUTPUT: s1 and s2 are not equal, random */
L100:
    *s1 = (integer) (rand_() * *m) + 1;
    if (*s1 < 1 || *s1 > *m) {
	goto L100;
    }
L101:
    *s2 = (integer) (rand_() * *m) + 1;
    if (*s2 < 1 || *s2 > *m) {
	goto L101;
    }
    if (*s1 == *s2) {
	goto L101;
    }
    return 0;
} /* find_rand_pair__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int find_rand_pair_acc_rej__(integer *mm, integer *m, 
	doublereal *v, doublereal *max_k__, S_fp kernel, integer *s1, integer 
	*s2)
{
    static doublereal k, p;
    extern doublereal rand_(void);
    extern /* Subroutine */ int find_rand_pair__(integer *, integer *, 
	    integer *);

/* INPUT: physical dimension of V */
/* INPUT: logical dimension of V */
/* INPUT: array of particle volumes */
/* INPUT: maximum value of the kernel */
/* INPUT: kernel function */
/*         particles with V(s1/s2) != 0 */
/* OUTPUT: s1 and s2 are not equal, random */
    /* Parameter adjustments */
    --v;

    /* Function Body */
L200:
    find_rand_pair__(m, s1, s2);
/* test particles s1, s2 */
    (*kernel)(&v[*s1], &v[*s2], &k);
    p = k / *max_k__;
/* collision probability */
    if ((doublereal) rand_() > p) {
	goto L200;
    }
    return 0;
} /* find_rand_pair_acc_rej__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int coagulate_(integer *mm, integer *m, doublereal *v, 
	doublereal *v_comp__, integer *n_bin__, doublereal *bin_v__, 
	doublereal *bin_r__, doublereal *bin_g__, integer *bin_n__, 
	doublereal *dlnr, integer *s1, integer *s2, logical *bin_change__)
{
    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static integer k1, k2, kn;
    extern /* Subroutine */ int exit_(integer *), particle_in_bin__(
	    doublereal *, integer *, doublereal *, integer *);

    /* Fortran I/O blocks */
    static cilist io___13 = { 0, 6, 0, 0, 0 };


/* INPUT: physical dimension of V */
/* INPUT/OUTPUT: logical dimension of V */
/* INPUT/OUTPUT: particle volumes */
/* INPUT: computational volume */
/* INPUT: number of bins */
/* INPUT: volume of particles in bins */
/* INPUT: radius of particles in bins */
/* INPUT/OUTPUT: mass in bins */
/* INPUT/OUTPUT: number in bins */
/* INPUT: bin scale factor */
/* INPUT: first particle to coagulate */
/* INPUT: second particle to coagulate */
/*         or a filled bin became empty */
/* OUTPUT: whether an empty bin filled, */
    /* Parameter adjustments */
    --v;
    --bin_n__;
    --bin_g__;
    --bin_r__;
    --bin_v__;

    /* Function Body */
    *bin_change__ = FALSE_;
/* remove s1 and s2 from bins */
    particle_in_bin__(&v[*s1], n_bin__, &bin_v__[1], &k1);
    particle_in_bin__(&v[*s2], n_bin__, &bin_v__[1], &k2);
    --bin_n__[k1];
    --bin_n__[k2];
/* DEBUG */
/*      write(*,*)'V(s1),V(s2),k1,k2,n(k1),n(k2) = ', */
/*     &     V(s1), V(s2), k1, k2, bin_n(k1), bin_n(k2) */
/*      write(*,*)'r(k1),r(k2) = ', bin_r(k1), bin_r(k2) */
/*      if (bin_n(k1) .eq. 0) write(*,*)'bin became empty:', k1 */
/*      if (bin_n(k2) .eq. 0) write(*,*)'bin became empty:', k2 */
/* DEBUG */
    bin_g__[k1] -= v[*s1];
    bin_g__[k2] -= v[*s2];
    if (bin_n__[k1] < 0 || bin_n__[k2] < 0) {
	s_wsle(&io___13);
	do_lio(&c__9, &c__1, "ERROR: invalid bin_n", (ftnlen)20);
	e_wsle();
	exit_(&c__2);
    }
    v[*s1] += v[*s2];
/* add particle 2 onto particle 1 */
    v[*s2] = v[*m];
/* shift the last particle into empty slot */
    --(*m);
/* add new particle to bins */
/* shorten array */
    particle_in_bin__(&v[*s1], n_bin__, &bin_v__[1], &kn);
/* DEBUG */
/*      write(*,*)'V(s1),kn,n(kn) = ', V(s1), kn, bin_n(kn) */
/*      write(*,*)'r(kn) = ', bin_r(kn) */
/*      if (bin_n(kn) .le. 0) write(*,*)'bin became full:', kn */
/* DEBUG */
    ++bin_n__[kn];
    bin_g__[kn] += v[*s1];
    if (bin_n__[k1] == 0 || bin_n__[k2] == 0) {
	*bin_change__ = TRUE_;
    }
    if (bin_n__[kn] == 1 && kn != k1 && kn != k2) {
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
    return 0;
} /* coagulate_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int maybe_coag_pair__(integer *mm, integer *m, doublereal *v,
	 doublereal *v_comp__, integer *n_bin__, doublereal *bin_v__, 
	doublereal *bin_r__, doublereal *bin_g__, integer *bin_n__, 
	doublereal *dlnr, doublereal *del_t__, integer *n_samp__, S_fp kernel,
	 logical *did_coag__, logical *bin_change__)
{
    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal k, p;
    extern /* Subroutine */ int coagulate_(integer *, integer *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    , integer *, doublereal *, integer *, integer *, logical *);
    static integer s1, s2;
    extern doublereal rand_(void);
    static doublereal expo;
    extern /* Subroutine */ int find_rand_pair__(integer *, integer *, 
	    integer *);

/* INPUT: physical dimension of V */
/* INPUT/OUTPUT: logical dimension of V */
/* INPUT/OUTPUT: particle volumes */
/* INPUT: computational volume */
/* INPUT: number of bins */
/* INPUT: volume of particles in bins */
/* INPUT: radius of particles in bins */
/* INPUT/OUTPUT: mass in bins */
/* INPUT/OUTPUT: number in bins */
/* INPUT: bin scale factor */
/* INPUT: timestep */
/* INPUT: number of samples per timestep */
/* INPUT: kernel function */
/* OUTPUT: whether a coagulation occured */
/* OUTPUT: whether bin structure changed */
    /* Parameter adjustments */
    --v;
    --bin_n__;
    --bin_g__;
    --bin_r__;
    --bin_v__;

    /* Function Body */
    find_rand_pair__(m, &s1, &s2);
/* test particles s1, s2 */
    (*kernel)(&v[s1], &v[s2], &k);
    expo = k * 1. / *v_comp__ * *del_t__ * (*m * (*m - 1) / 2.) / *n_samp__;
    p = 1. - exp(-expo);
/* probability of coagulation */
    *bin_change__ = FALSE_;
    if ((doublereal) rand_() < p) {
	coagulate_(mm, m, &v[1], v_comp__, n_bin__, &bin_v__[1], &bin_r__[1], 
		&bin_g__[1], &bin_n__[1], dlnr, &s1, &s2, bin_change__);
	*did_coag__ = TRUE_;
    } else {
	*did_coag__ = FALSE_;
    }
    return 0;
} /* maybe_coag_pair__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int kernel_avg__(integer *mm, integer *m, doublereal *v, 
	S_fp kernel, integer *n_samp__, doublereal *k_avg__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__;
    static doublereal k;
    static integer s1, s2;
    static doublereal k_sum__;
    extern /* Subroutine */ int find_rand_pair__(integer *, integer *, 
	    integer *);

/* FIXME: use binned data instead */
/* INPUT: physical dimension of V */
/* INPUT: logical dimension of V */
/* INPUT: array of particle volumes */
/* INPUT: kernel function */
/* INPUT: number of samples to use (squared) */
/* OUTPUT: estimated average of kernel values */
    /* Parameter adjustments */
    --v;

    /* Function Body */
    k_sum__ = 0.;
/* Computing 2nd power */
    i__2 = *n_samp__;
    i__1 = i__2 * i__2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	find_rand_pair__(m, &s1, &s2);
	(*kernel)(&v[s1], &v[s2], &k);
	k_sum__ += k;
    }
/* Computing 2nd power */
    i__1 = *n_samp__;
    *k_avg__ = k_sum__ / (i__1 * i__1);
    return 0;
} /* kernel_avg__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int double_(integer *mm, integer *m, doublereal *v, 
	doublereal *v_comp__, integer *n_bin__, doublereal *bin_v__, 
	doublereal *bin_r__, doublereal *bin_g__, integer *bin_n__, 
	doublereal *dlnr)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int exit_(integer *);

    /* Fortran I/O blocks */
    static cilist io___25 = { 0, 6, 0, 0, 0 };


/* INPUT: physical dimension of V */
/* INPUT/OUTPUT: logical dimension of V */
/* INPUT/OUTPUT: particle volumes */
/* INPUT/OUTPUT: computational volume */
/* INPUT: number of bins */
/* INPUT: volume of particles in bins */
/* INPUT: radius of particles in bins */
/* INPUT/OUTPUT: mass in bins */
/* INPUT/OUTPUT: number in bins */
/* INPUT: bin scale factor */
/* only double if we have enough space to do so */
    /* Parameter adjustments */
    --v;
    --bin_n__;
    --bin_g__;
    --bin_r__;
    --bin_v__;

    /* Function Body */
    if (*m > *mm / 2) {
	s_wsle(&io___25);
	do_lio(&c__9, &c__1, "ERROR: double without enough space", (ftnlen)34)
		;
	e_wsle();
	exit_(&c__2);
    }
/* double V and associated structures */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[i__ + *m] = v[i__];
    }
    *m <<= 1;
    *v_comp__ *= 2.;
/* double bin structures */
    i__1 = *n_bin__;
    for (i__ = 1; i__ <= i__1; ++i__) {
	bin_g__[i__] *= 2.;
	bin_n__[i__] <<= 1;
    }
    return 0;
} /* double_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int est_k_max__(integer *n_bin__, doublereal *bin_v__, 
	integer *bin_n__, S_fp kernel, doublereal *k_max__, logical *
	use_bin__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal k;

/* INPUT: number of bins */
/* INPUT: volume of particles in bins */
/* INPUT: number in each bin */
/* INPUT: kernel function */
/* OUTPUT: maximum kernel value */
/*      write(*,*)'RECALCULATING K_MAX' */
/* use_bin starts as non-empty bins */
    /* Parameter adjustments */
    --use_bin__;
    --bin_n__;
    --bin_v__;

    /* Function Body */
    i__1 = *n_bin__;
    for (i__ = 1; i__ <= i__1; ++i__) {
	use_bin__[i__] = bin_n__[i__] > 0;
    }
/* add all bins downstream of non-empty bins */
    i__1 = *n_bin__;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (use_bin__[i__]) {
	    use_bin__[i__ - 1] = TRUE_;
	}
    }
/* add all bins upstream of non-empty bins */
    for (i__ = *n_bin__ - 1; i__ >= 1; --i__) {
	if (use_bin__[i__]) {
	    use_bin__[i__ + 1] = TRUE_;
	}
    }
    *k_max__ = 0.;
    i__1 = *n_bin__;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (use_bin__[i__]) {
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		if (use_bin__[j]) {
		    (*kernel)(&bin_v__[i__], &bin_v__[j], &k);
		    if (k > *k_max__) {
			*k_max__ = k;
		    }
		}
	    }
	}
    }
    return 0;
} /* est_k_max__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int particle_in_bin__(doublereal *v, integer *n_bin__, 
	doublereal *bin_v__, integer *k)
{
/* FIXME: for log-spaced bins we can do this without search */
/* INPUT: volume of particle */
/* INPUT: number of bins */
/* INPUT: volume of particles in bins */
/* OUTPUT: bin number containing particle */
    /* Parameter adjustments */
    --bin_v__;

    /* Function Body */
    *k = 0;
L300:
    ++(*k);
/*      write(*,*)'k,bin_v(k) = ', k, bin_v(k) */
    if (*k < *n_bin__ && *v > (bin_v__[*k] + bin_v__[*k + 1]) / 2.) {
	goto L300;
    }
    return 0;
} /* particle_in_bin__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int moments_(integer *mm, integer *m, doublereal *v, 
	doublereal *v_comp__, integer *n_bin__, doublereal *bin_v__, 
	doublereal *bin_r__, doublereal *bin_g__, integer *bin_n__, 
	doublereal *dlnr)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k;
    extern /* Subroutine */ int particle_in_bin__(doublereal *, integer *, 
	    doublereal *, integer *);

/* INPUT: physical dimension of V */
/* INPUT: logical dimension of V */
/* INPUT: particle volumes */
/* INPUT: computational volume */
/* INPUT: number of bins */
/* INPUT: volume of particles in bins */
/* INPUT: radius of particles in bins */
/* OUTPUT: mass in bins */
/* OUTPUT: number in bins */
/* INPUT: bin scale factor */
    /* Parameter adjustments */
    --v;
    --bin_n__;
    --bin_g__;
    --bin_r__;
    --bin_v__;

    /* Function Body */
    i__1 = *n_bin__;
    for (k = 1; k <= i__1; ++k) {
	bin_g__[k] = 0.;
	bin_n__[k] = 0;
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	particle_in_bin__(&v[i__], n_bin__, &bin_v__[1], &k);
	bin_g__[k] += v[i__];
	++bin_n__[k];
    }
    return 0;
} /* moments_ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int check_event__(doublereal *time, doublereal *interval, 
	doublereal *last_time__, logical *do_event__)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double d_int(doublereal *);

    /* Local variables */
    static doublereal interval_below__;

/* INPUT: cubin_rent time */
/* INPUT: how often the event should be done */
/* INPUT/OUTPUT: when the event was last done */
/* OUTPUT: whether the event should be done */
    if (*time == 0.) {
	*last_time__ = 0.;
	*do_event__ = TRUE_;
    } else {
	d__1 = *time / *interval;
	interval_below__ = d_int(&d__1) * *interval;
	if (*last_time__ < interval_below__) {
	    *last_time__ = *time;
	    *do_event__ = TRUE_;
	} else {
	    *do_event__ = FALSE_;
	}
    }
    return 0;
} /* check_event__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int print_header__(integer *n_loop__, integer *n_bin__, 
	integer *n_time__)
{
    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Fortran I/O blocks */
    static cilist io___33 = { 0, 30, 0, "(a10,i10)", 0 };
    static cilist io___34 = { 0, 30, 0, "(a10,i10)", 0 };
    static cilist io___35 = { 0, 30, 0, "(a10,i10)", 0 };


/* INPUT: number of loops */
/* INPUT: number of bins */
/* INPUT: number of times */
    s_wsfe(&io___33);
    do_fio(&c__1, "n_loop", (ftnlen)6);
    do_fio(&c__1, (char *)&(*n_loop__), (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___34);
    do_fio(&c__1, "n_bin", (ftnlen)5);
    do_fio(&c__1, (char *)&(*n_bin__), (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___35);
    do_fio(&c__1, "n_time", (ftnlen)6);
    do_fio(&c__1, (char *)&(*n_time__), (ftnlen)sizeof(integer));
    e_wsfe();
    return 0;
} /* print_header__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int print_info__(doublereal *time, doublereal *v_comp__, 
	integer *n_bin__, doublereal *bin_v__, doublereal *bin_r__, 
	doublereal *bin_g__, integer *bin_n__, doublereal *dlnr)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer k;

    /* Fortran I/O blocks */
    static cilist io___36 = { 0, 30, 0, "(a10,e14.5)", 0 };
    static cilist io___38 = { 0, 30, 0, "(i8,3e14.5)", 0 };


/* INPUT: cubin_rent simulation time */
/* INPUT: computational volume */
/* INPUT: number of bins */
/* INPUT: volume of particles in bins */
/* INPUT: radius of particles in bins */
/* INPUT: mass in bins */
/* INPUT: number in bins */
/* INPUT: bin scale factor */
    /* Parameter adjustments */
    --bin_n__;
    --bin_g__;
    --bin_r__;
    --bin_v__;

    /* Function Body */
    s_wsfe(&io___36);
    do_fio(&c__1, "time", (ftnlen)4);
    do_fio(&c__1, (char *)&(*time), (ftnlen)sizeof(doublereal));
    e_wsfe();
    i__1 = *n_bin__;
    for (k = 1; k <= i__1; ++k) {
	s_wsfe(&io___38);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&bin_r__[k], (ftnlen)sizeof(doublereal));
	d__1 = bin_n__[k] / *v_comp__ / *dlnr;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	d__2 = bin_g__[k] / *v_comp__ / *dlnr;
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    return 0;
} /* print_info__ */

