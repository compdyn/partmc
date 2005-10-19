/* run_golovin_exact.f -- translated by f2c (version 20031025).
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
static integer c__1 = 1;
static integer c__160 = 160;
static integer c__3 = 3;
static doublereal c_b7 = 1e3;
static doublereal c_b9 = 1e9;
static doublereal c_b10 = 4.1886e-15;
static doublereal c_b12 = 600.;
static doublereal c_b13 = 60.;
static doublereal c_b14 = 1.;

/* Exact solution with Golovin kernel. */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    integer i__1;
    olist o__1;

    /* Builtin functions */
    integer f_open(olist *), i_dnnt(doublereal *);

    /* Local variables */
    extern /* Subroutine */ int mc_exact__(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, U_fp, doublereal *, 
	    doublereal *, integer *, doublereal *);
    extern /* Subroutine */ int soln_golovin_exp__();
    extern /* Subroutine */ int make_grid__(integer *, integer *, doublereal *
	    , doublereal *, doublereal *, doublereal *);
    static doublereal dlnr;
    extern /* Subroutine */ int print_header__(integer *, integer *, integer *
	    );
    static doublereal bin_g__[160], bin_n__[160], bin_r__[160], bin_v__[160];
    static integer i_loop__;

/* number of bins */
/* number of loops */
/* scale factor for bins */
/* total simulation time (seconds) */
/* particle density (kg/m^3) */
/* particle number concentration (#/ */
/* interval between printing (s) */
/* mean volume of initial distributi */
/* computational volume (dummy value */
    o__1.oerr = 0;
    o__1.ounit = 30;
    o__1.ofnmlen = 19;
    o__1.ofnm = "out_golovin_exact.d";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    i__1 = i_dnnt(&c_b2) + 1;
    print_header__(&c__1, &c__160, &i__1);
    for (i_loop__ = 1; i_loop__ <= 1; ++i_loop__) {
	make_grid__(&c__160, &c__3, &c_b7, bin_v__, bin_r__, &dlnr);
	mc_exact__(&c__160, bin_v__, bin_r__, bin_g__, bin_n__, &dlnr, &c_b9, 
		&c_b10, &c_b7, (U_fp)soln_golovin_exp__, &c_b12, &c_b13, &
		i_loop__, &c_b14);
    }
    return 0;
} /* MAIN__ */

/* Main program alias */ int montecarlo_ () { MAIN__ (); return 0; }
