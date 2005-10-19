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

static double c_b3 = 1.;
static double c_b4 = 0.;

/* Golovin coagulation kernel. */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<       subroutine kernel_golovin(a, b, k) >*/
/* Subroutine */ int kernel_golovin__(double *a, double *b, 
	double *k)
{
/*<       real*8 a  ! INPUT: volume of first particle >*/
/*<       real*8 b  ! INPUT: volume of second particle >*/
/*<       real*8 k  ! OUTPUT: coagulation kernel >*/
/*<       real*8 beta_1 >*/
/*<       parameter (beta_1 = 1000d0) >*/
/*<       k = beta_1 * (a + b) >*/
    *k = (*a + *b) * 1e3;
/*<       return >*/
    return 0;
/*<       end >*/
} /* kernel_golovin__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<    >*/
/* Subroutine */ int soln_golovin_exp__(int *n_bin__, double *bin_v__,
	 double *bin_r__, double *bin_g__, int *bin_n__, 
	double *dlnr, double *time, double *n_0__, double *
	v_0__, double *rho_p__, double *v_comp__)
{
    /* System generated locals */
    int i__1;
    double d__1;

    /* Builtin functions */
    double exp(double), sqrt(double);

    /* Local variables */
    static double b;
    static int k;
    static double t, x, nn, tau, rat_v__, beta_1__;
    extern /* Subroutine */ int bessi1_(double *, double *), 
	    kernel_golovin__(double *, double *, double *);

/*<       int n_bin        ! INPUT: number of bins >*/
/*<       real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins >*/
/*<       real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins >*/
/*<       real*8 bin_g(n_bin)  ! OUTPUT: mass in bins >*/
/*<       int bin_n(n_bin) ! OUTPUT: number in bins >*/
/*<       real*8 dlnr          ! INPUT: bin scale factor >*/
/*<       real*8 time          ! INPUT: cubin_rent time >*/
/*<       real*8 N_0           ! INPUT: particle number concentration (#/m^3 >*/
/*<       real*8 V_0           ! INPUT: >*/
/*<       real*8 rho_p         ! INPUT: particle density (kg/m^3) >*/
/*<       real*8 V_comp        ! INPUT: computational volume >*/
/*<       real*8 beta_1, tau, T, rat_v, nn, b, x >*/
/*<       int k >*/
/*<       real*8 pi >*/
/*<       parameter (pi = 3.14159265358979323846d0) >*/
/*<       call kernel_golovin(1d0, 0d0, beta_1) >*/
    /* Parameter adjustments */
    --bin_n__;
    --bin_g__;
    --bin_r__;
    --bin_v__;

    /* Function Body */
    kernel_golovin__(&c_b3, &c_b4, &beta_1__);
/*<       if (time .eq. 0d0) then >*/
    if (*time == 0.) {
/*<          do k = 1,n_bin >*/
	i__1 = *n_bin__;
	for (k = 1; k <= i__1; ++k) {
/*<    >*/
/* Computing 3rd power */
	    d__1 = bin_r__[k] * 2.;
	    bin_n__[k] = (int) (d__1 * (d__1 * d__1) * 1.5707963267948966 
		    * *n_0__ / *v_0__ * exp(-(bin_v__[k] / *v_0__)));
/*<          enddo >*/
	}
/*<       else >*/
    } else {
/*<          tau = N_0 * V_0 * beta_1 * time >*/
	tau = *n_0__ * *v_0__ * beta_1__ * *time;
/*<          T = 1 - exp(-tau) >*/
	t = 1 - exp(-tau);
/*<          do k = 1,n_bin >*/
	i__1 = *n_bin__;
	for (k = 1; k <= i__1; ++k) {
/*<             rat_v = bin_v(k) / V_0 >*/
	    rat_v__ = bin_v__[k] / *v_0__;
/*<             x = 2d0 * rat_v * sqrt(T) >*/
	    x = rat_v__ * 2. * sqrt(t);
/*<             if (x .lt. 500d0) then >*/
	    if (x < 500.) {
/*<                call bessi1(x, b) >*/
		bessi1_(&x, &b);
/*<             else >*/
	    } else {
/*<                b = 0d0 >*/
		b = 0.;
/*<             endif >*/
	    }
/*<    >*/
	    nn = *n_0__ / bin_v__[k] * (1. - t) / sqrt(t) * exp(-((t + 1.) * 
		    rat_v__)) * b;
/*<             bin_n(k) = pi/2d0 * (2d0*bin_r(k))**3 * nn >*/
/* Computing 3rd power */
	    d__1 = bin_r__[k] * 2.;
	    bin_n__[k] = (int) (d__1 * (d__1 * d__1) * 1.5707963267948966 
		    * nn);
/*<          enddo >*/
	}
/*<       endif >*/
    }
/*<       do k = 1,n_bin >*/
    i__1 = *n_bin__;
    for (k = 1; k <= i__1; ++k) {
/*<          bin_g(k) = pi/6d0 * rho_p * (2d0*bin_r(k))**3 * bin_n(k) >*/
/* Computing 3rd power */
	d__1 = bin_r__[k] * 2.;
	bin_g__[k] = *rho_p__ * .52359877559829882 * (d__1 * (d__1 * d__1)) * 
		bin_n__[k];
/*<       enddo >*/
    }
/*<       do k = 1,n_bin >*/
    i__1 = *n_bin__;
    for (k = 1; k <= i__1; ++k) {
/*<          bin_g(k) = bin_g(k) * dlnr * V_comp >*/
	bin_g__[k] = bin_g__[k] * *dlnr * *v_comp__;
/*<          bin_n(k) = bin_n(k) * dlnr * V_comp >*/
	bin_n__[k] = (int) (bin_n__[k] * *dlnr * *v_comp__);
/*<       enddo >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* soln_golovin_exp__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<       subroutine bessi1(x, r) >*/
/* Subroutine */ int bessi1_(double *x, double *r__)
{
    /* Initialized data */

    static double p2 = .87890594;
    static double p3 = .51498869;
    static double p4 = .15084934;
    static double p5 = .02658733;
    static double p6 = .00301532;
    static double p7 = 3.2411e-4;
    static double q1 = .39894228;
    static double q2 = -.03988024;
    static double q3 = -.00362018;
    static double q4 = .00163801;
    static double q5 = -.01031555;
    static double q6 = .02282967;
    static double q7 = -.02895312;
    static double q8 = .01787654;
    static double q9 = -.00420059;
    static double p1 = .5;

    /* System generated locals */
    double d__1;

    /* Builtin functions */
    double exp(double), sqrt(double);

    /* Local variables */
    static double y, ax;

/* bessel function */
/*<       real*8 x   ! INPUT: function argument >*/
/*<       real*8 r   ! OUTPUT: function value >*/
/*<       real*8 ax >*/
/*<       real*8 p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y >*/
/*<    >*/
/*<    >*/
/*<       if(abs(x) .lt. 3.75d0) then >*/
    if (abs(*x) < 3.75) {
/*<          y = (x / 3.75d0)**2 >*/
/* Computing 2nd power */
	d__1 = *x / 3.75;
	y = d__1 * d__1;
/*<          r = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))) >*/
	*r__ = *x * (p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y 
		* p7))))));
/*<       else >*/
    } else {
/*<          ax = abs(x) >*/
	ax = abs(*x);
/*<          y = 3.75d0 / ax >*/
	y = 3.75 / ax;
/*<    >*/
	*r__ = exp(ax) / sqrt(ax) * (q1 + y * (q2 + y * (q3 + y * (q4 + y * (
		q5 + y * (q6 + y * (q7 + y * (q8 + y * q9))))))));
/*<          if (x .lt. 0d0) r = -r >*/
	if (*x < 0.) {
	    *r__ = -(*r__);
	}
/*<       endif >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* bessi1_ */

