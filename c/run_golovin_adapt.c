/* run_golovin_adapt.f -- translated by f2c (version 20031025).
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

static doublereal c_b2 = 10.;
static integer c__10 = 10;
static integer c__160 = 160;
static integer c__3 = 3;
static doublereal c_b8 = 1e3;
static integer c__10000 = 10000;
static doublereal c_b13 = 600.;
static doublereal c_b14 = 60.;
static doublereal c_b15 = .01;
static doublereal c_b16 = .005;
static doublereal c_b17 = 1.;

/* Simulation with Golovin kernel and adaptive timestepping. */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    olist o__1;

    /* Builtin functions */
    integer f_open(olist *), i_dnnt(doublereal *);
    double exp(doublereal);

    /* Local variables */
    extern /* Subroutine */ int mc_adapt__(integer *, integer *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    , integer *, doublereal *, U_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *), 
	    compute_volumes__(integer *, integer *, doublereal *, doublereal *
	    , doublereal *, doublereal *, integer *);
    static integer k, m;
    static doublereal v[10000];
    extern /* Subroutine */ int make_grid__(integer *, integer *, doublereal *
	    , doublereal *, doublereal *, doublereal *);
    static doublereal dlnr;
    extern /* Subroutine */ int print_header__(integer *, integer *, integer *
	    );
    static doublereal bin_g__[160];
    static integer bin_n__[160];
    static doublereal bin_r__[160], n_ini__[160], bin_v__[160];
    extern /* Subroutine */ int srand_(integer *);
    static integer i_loop__;
    static doublereal v_comp__;
    extern /* Subroutine */ int kernel_golovin__();

/* number of particles */
/* number of bins */
/* number of loops */
/* scale factor for bins */
/* total simulation time (seconds) */
/* particle density (kg/m^3) */
/* particle number concentration ( */
/* interval between printing (s) */
/* maximum coagulation probability */
/* maximum sampling ratio per time */
/* maximum timestep */
/* mean volume of initial distribu */
    o__1.oerr = 0;
    o__1.ounit = 30;
    o__1.ofnmlen = 19;
    o__1.ofnm = "out_golovin_adapt.d";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    i__1 = i_dnnt(&c_b2) + 1;
    print_header__(&c__10, &c__160, &i__1);
    srand_(&c__10);
    for (i_loop__ = 1; i_loop__ <= 10; ++i_loop__) {
	make_grid__(&c__160, &c__3, &c_b8, bin_v__, bin_r__, &dlnr);
/* define initial exponential distribution */
	for (k = 1; k <= 160; ++k) {
/* Computing 3rd power */
	    d__1 = bin_r__[k - 1] * 2.;
	    n_ini__[k - 1] = d__1 * (d__1 * d__1) * 1.5707963267948966 * 
		    10000 / 4.1886e-15 * exp(-(bin_v__[k - 1] / 4.1886e-15));
	}
	compute_volumes__(&c__160, &c__10000, n_ini__, bin_r__, &dlnr, v, &m);
	v_comp__ = m / 1e9;
	mc_adapt__(&c__10000, &m, v, &v_comp__, &c__160, bin_v__, bin_r__, 
		bin_g__, bin_n__, &dlnr, (U_fp)kernel_golovin__, &c_b13, &
		c_b14, &c_b15, &c_b16, &c_b17, &i_loop__);
    }
    return 0;
} /* MAIN__ */

/* Main program alias */ int montecarlo_ () { MAIN__ (); return 0; }
