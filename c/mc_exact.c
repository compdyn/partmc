/* mc_exact.f -- translated by f2c (version 20031025).
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

/* Exact solution output. */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine */ int mc_exact__(integer *n_bin__, doublereal *bin_v__, 
	doublereal *bin_r__, doublereal *bin_g__, integer *bin_n__, 
	doublereal *dlnr, doublereal *n_0__, doublereal *v_0__, doublereal *
	rho_p__, S_fp soln, doublereal *t_max__, doublereal *t_print__, 
	integer *loop, doublereal *v_comp__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int print_info__(doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *);
    static doublereal time;
    static integer i_time__, n_time__;

/* FIXME: N_0 and V_0 are really parameters for the initial value */
/* of the particle distribution. They should be replaced by a n_par */
/* params() pair. */
/* INPUT: number of bins */
/* INPUT: volume of bins */
/* INPUT: radius of bins */
/* OUTPUT: mass in bins */
/* OUTPUT: number in bins */
/* INPUT: bin scale factor */
/* INPUT: particle number concentration (#/m^3 */
/* INPUT: */
/* INPUT: particle density (kg/m^3) */
/* INPUT: exact solution procedure */
/* INPUT: total simulation time */
/* INPUT: interval to print info (seconds) */
/* INPUT: loop number of run */
/* INPUT: computational volume */
    /* Parameter adjustments */
    --bin_n__;
    --bin_g__;
    --bin_r__;
    --bin_v__;

    /* Function Body */
    n_time__ = (integer) (*t_max__ / *t_print__);
    i__1 = n_time__;
    for (i_time__ = 0; i_time__ <= i__1; ++i_time__) {
	time = (doublereal) i_time__ / (doublereal) n_time__ * *t_max__;
	(*soln)(n_bin__, &bin_v__[1], &bin_r__[1], &bin_g__[1], &bin_n__[1], 
		dlnr, &time, n_0__, v_0__, rho_p__, v_comp__);
	print_info__(&time, v_comp__, n_bin__, &bin_v__[1], &bin_r__[1], &
		bin_g__[1], &bin_n__[1], dlnr);
    }
    return 0;
} /* mc_exact__ */

