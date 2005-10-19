/* mc_adapt.f -- translated by f2c (version 20031025).
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

static integer c__1 = 1;

/* Monte Carlo with adaptive timestep. */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int mc_adapt__(integer *mm, integer *m, doublereal *v, 
	doublereal *v_comp__, integer *n_bin__, doublereal *bin_v__, 
	doublereal *bin_r__, doublereal *bin_g__, integer *bin_n__, 
	doublereal *dlnr, U_fp kernel, doublereal *t_max__, doublereal *
	t_print__, doublereal *p_max__, doublereal *r_samp_max__, doublereal *
	del_t_max__, integer *loop)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static logical did_coag__;
    static doublereal last_print_time__;
    extern /* Subroutine */ int cpu_time__(doublereal *);
    static logical do_print__;
    extern /* Subroutine */ int est_k_max__(integer *, doublereal *, integer *
	    , U_fp, doublereal *);
    static logical bin_change__;
    static doublereal t_per_samp__;
    extern /* Subroutine */ int print_info__(doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *), check_event__(doublereal *, doublereal *, 
	    doublereal *, logical *);
    static doublereal time;
    extern /* Subroutine */ int compute_n_samp_del_t__(integer *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, integer 
	    *, doublereal *);
    static doublereal del_t__, t_end__, k_max__;
    static integer n_coag__, i_samp__;
    extern /* Subroutine */ int double_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *);
    static integer n_samp__;
    static doublereal t_loop__;
    extern /* Subroutine */ int maybe_coag_pair__(integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    U_fp, logical *, logical *);
    static doublereal t_start__;
    extern /* Subroutine */ int moments_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___15 = { 0, 6, 0, "(a6,a6,a6,a6,a10,a9,a11,a9)", 0 };
    static cilist io___16 = { 0, 6, 0, "(i6,f6.1,f6.3,i6,e10.3,i9,e11.3,i9)", 
	    0 };


/* INPUT: physical dimension of V */
/* INPUT/OUTPUT: logical dimension of V */
/* INPUT/OUTPUT: particle volumes */
/* INPUT/OUTPUT: computational volume */
/* INPUT: number of bins */
/* INPUT: volume of particles in bins */
/* INPUT: radius of particles in bins */
/* OUTPUT: mass in bins */
/* OUTPUT: number in bins */
/* INPUT: bin scale factor */
/* INPUT: kernel function */
/* INPUT: final time (seconds) */
/* INPUT: interval to print info (seconds) */
/* INPUT: maximum coagulation probability */
/* INPUT: maximum sampling ratio per timestep */
/* INPUT: maximum timestep */
/* INPUT: loop number of run */
    /* Parameter adjustments */
    --v;
    --bin_n__;
    --bin_g__;
    --bin_r__;
    --bin_v__;

    /* Function Body */
    time = 0.;
    n_coag__ = 0;
    moments_(mm, m, &v[1], v_comp__, n_bin__, &bin_v__[1], &bin_r__[1], &
	    bin_g__[1], &bin_n__[1], dlnr);
    check_event__(&time, t_print__, &last_print_time__, &do_print__);
    if (do_print__) {
	print_info__(&time, v_comp__, n_bin__, &bin_v__[1], &bin_r__[1], &
		bin_g__[1], &bin_n__[1], dlnr);
    }
    est_k_max__(n_bin__, &bin_v__[1], &bin_n__[1], (U_fp)kernel, &k_max__);
    while(time < *t_max__) {
	cpu_time__(&t_start__);
	compute_n_samp_del_t__(m, v_comp__, &k_max__, p_max__, r_samp_max__, 
		del_t_max__, &n_samp__, &del_t__);
	i__1 = n_samp__;
	for (i_samp__ = 1; i_samp__ <= i__1; ++i_samp__) {
	    maybe_coag_pair__(mm, m, &v[1], v_comp__, n_bin__, &bin_v__[1], &
		    bin_r__[1], &bin_g__[1], &bin_n__[1], dlnr, &del_t__, &
		    n_samp__, (U_fp)kernel, &did_coag__, &bin_change__);
	    if (did_coag__) {
		++n_coag__;
	    }
	    if (bin_change__) {
		est_k_max__(n_bin__, &bin_v__[1], &bin_n__[1], (U_fp)kernel, &
			k_max__);
	    }
	    if (*m < *mm / 2) {
		double_(mm, m, &v[1], v_comp__, n_bin__, &bin_v__[1], &
			bin_r__[1], &bin_g__[1], &bin_n__[1], dlnr);
	    }
	}
	time += del_t__;
	check_event__(&time, t_print__, &last_print_time__, &do_print__);
	if (do_print__) {
	    print_info__(&time, v_comp__, n_bin__, &bin_v__[1], &bin_r__[1], &
		    bin_g__[1], &bin_n__[1], dlnr);
	}
	cpu_time__(&t_end__);
	t_loop__ = t_end__ - t_start__;
	t_per_samp__ = t_loop__ / n_samp__;
	s_wsfe(&io___15);
	do_fio(&c__1, "loop", (ftnlen)4);
	do_fio(&c__1, "time", (ftnlen)4);
	do_fio(&c__1, "del_t", (ftnlen)5);
	do_fio(&c__1, "M", (ftnlen)1);
	do_fio(&c__1, "k_max", (ftnlen)5);
	do_fio(&c__1, "n_samp", (ftnlen)6);
	do_fio(&c__1, "t_per_samp", (ftnlen)10);
	do_fio(&c__1, "n_coag", (ftnlen)6);
	e_wsfe();
	s_wsfe(&io___16);
	do_fio(&c__1, (char *)&(*loop), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&time, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&del_t__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*m), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&k_max__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&n_samp__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&t_per_samp__, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&n_coag__, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    return 0;
} /* mc_adapt__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int compute_n_samp_del_t__(integer *m, doublereal *v_comp__, 
	doublereal *k_max__, doublereal *p_max__, doublereal *r_samp_max__, 
	doublereal *del_t_max__, integer *n_samp__, doublereal *del_t__)
{
    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static doublereal c__, r_samp__;

/* INPUT: number of particles */
/* INPUT: computational volume */
/* INPUT: maximum kernel value */
/* INPUT: maximum coagulation probability */
/* INPUT: maximum sampling ratio per timestep */
/* INPUT: maximum timestep */
/* OUTPUT: number of samples per timestep */
/* OUTPUT: timestep */
    c__ = -(*k_max__ * 1. / *v_comp__ / log(1 - *p_max__));
    *del_t__ = *r_samp_max__ / c__;
    if (*del_t__ > *del_t_max__) {
	*del_t__ = *del_t_max__;
    }
    r_samp__ = *del_t__ * c__;
    *n_samp__ = (integer) (r_samp__ * *m * (*m - 1) / 2);
    return 0;
} /* compute_n_samp_del_t__ */

