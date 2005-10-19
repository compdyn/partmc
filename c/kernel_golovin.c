/* kernel_golovin.f -- translated by f2c (version 20031025).
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

static doublereal c_b3 = 1.;
static doublereal c_b4 = 0.;

/* Golovin coagulation kernel. */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int kernel_golovin__(doublereal *a, doublereal *b, 
	doublereal *k)
{
/* INPUT: volume of first particle */
/* INPUT: volume of second particle */
/* OUTPUT: coagulation kernel */
    *k = (*a + *b) * 1e3;
    return 0;
} /* kernel_golovin__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int soln_golovin_exp__(integer *n_bin__, doublereal *bin_v__,
	 doublereal *bin_r__, doublereal *bin_g__, integer *bin_n__, 
	doublereal *dlnr, doublereal *time, doublereal *n_0__, doublereal *
	v_0__, doublereal *rho_p__, doublereal *v_comp__)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal b;
    static integer k;
    static doublereal t, x, nn, tau, rat_v__, beta_1__;
    extern /* Subroutine */ int bessi1_(doublereal *, doublereal *), 
	    kernel_golovin__(doublereal *, doublereal *, doublereal *);

/* INPUT: number of bins */
/* INPUT: volume of particles in bins */
/* INPUT: radius of particles in bins */
/* OUTPUT: mass in bins */
/* OUTPUT: number in bins */
/* INPUT: bin scale factor */
/* INPUT: cubin_rent time */
/* INPUT: particle number concentration (#/m^3 */
/* INPUT: */
/* INPUT: particle density (kg/m^3) */
/* INPUT: computational volume */
    /* Parameter adjustments */
    --bin_n__;
    --bin_g__;
    --bin_r__;
    --bin_v__;

    /* Function Body */
    kernel_golovin__(&c_b3, &c_b4, &beta_1__);
    if (*time == 0.) {
	i__1 = *n_bin__;
	for (k = 1; k <= i__1; ++k) {
/* Computing 3rd power */
	    d__1 = bin_r__[k] * 2.;
	    bin_n__[k] = (integer) (d__1 * (d__1 * d__1) * 1.5707963267948966 
		    * *n_0__ / *v_0__ * exp(-(bin_v__[k] / *v_0__)));
	}
    } else {
	tau = *n_0__ * *v_0__ * beta_1__ * *time;
	t = 1 - exp(-tau);
	i__1 = *n_bin__;
	for (k = 1; k <= i__1; ++k) {
	    rat_v__ = bin_v__[k] / *v_0__;
	    x = rat_v__ * 2. * sqrt(t);
	    if (x < 500.) {
		bessi1_(&x, &b);
	    } else {
		b = 0.;
	    }
	    nn = *n_0__ / bin_v__[k] * (1. - t) / sqrt(t) * exp(-((t + 1.) * 
		    rat_v__)) * b;
/* Computing 3rd power */
	    d__1 = bin_r__[k] * 2.;
	    bin_n__[k] = (integer) (d__1 * (d__1 * d__1) * 1.5707963267948966 
		    * nn);
	}
    }
    i__1 = *n_bin__;
    for (k = 1; k <= i__1; ++k) {
/* Computing 3rd power */
	d__1 = bin_r__[k] * 2.;
	bin_g__[k] = *rho_p__ * .52359877559829882 * (d__1 * (d__1 * d__1)) * 
		bin_n__[k];
    }
    i__1 = *n_bin__;
    for (k = 1; k <= i__1; ++k) {
	bin_g__[k] = bin_g__[k] * *dlnr * *v_comp__;
	bin_n__[k] = (integer) (bin_n__[k] * *dlnr * *v_comp__);
    }
    return 0;
} /* soln_golovin_exp__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int bessi1_(doublereal *x, doublereal *r__)
{
    /* Initialized data */

    static doublereal p2 = .87890594;
    static doublereal p3 = .51498869;
    static doublereal p4 = .15084934;
    static doublereal p5 = .02658733;
    static doublereal p6 = .00301532;
    static doublereal p7 = 3.2411e-4;
    static doublereal q1 = .39894228;
    static doublereal q2 = -.03988024;
    static doublereal q3 = -.00362018;
    static doublereal q4 = .00163801;
    static doublereal q5 = -.01031555;
    static doublereal q6 = .02282967;
    static doublereal q7 = -.02895312;
    static doublereal q8 = .01787654;
    static doublereal q9 = -.00420059;
    static doublereal p1 = .5;

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal y, ax;

/* bessel function */
/* INPUT: function argument */
/* OUTPUT: function value */
    if (abs(*x) < 3.75) {
/* Computing 2nd power */
	d__1 = *x / 3.75;
	y = d__1 * d__1;
	*r__ = *x * (p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y 
		* p7))))));
    } else {
	ax = abs(*x);
	y = 3.75 / ax;
	*r__ = exp(ax) / sqrt(ax) * (q1 + y * (q2 + y * (q3 + y * (q4 + y * (
		q5 + y * (q6 + y * (q7 + y * (q8 + y * q9))))))));
	if (*x < 0.) {
	    *r__ = -(*r__);
	}
    }
    return 0;
} /* bessi1_ */

