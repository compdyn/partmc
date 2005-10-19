/* process_out.f -- translated by f2c (version 20031025).
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

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;

/* Process output data files. */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    address a__1[3];
    integer i__1, i__2, i__3, i__4[3];
    char ch__1[26];
    cilist ci__1;
    olist o__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_open(olist *), s_rsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_rsfe(void), s_wsfi(icilist *), e_wsfi(void), s_wsfe(
	    cilist *), e_wsfe(void);
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    extern integer len_trim__(char *, ftnlen);
    static char name_out_num_avg__[50];
    static integer i__;
    static doublereal n[2000000]	/* was [50][100][400] */;
    static char name_out_mass_avg__[50], n_time_str__[10], n_loop_str__[10], 
	    dum[100];
    static doublereal time;
    extern /* Subroutine */ int exit_(integer *);
    static char name_out_num__[50];
    static doublereal bin_g__[2000000]	/* was [50][100][400] */;
    static integer i_bin__;
    static doublereal g_avg__[40000]	/* was [100][400] */;
    static integer n_bin__;
    extern integer iargc_(void);
    static doublereal bin_r__[400], n_avg__[40000]	/* was [100][400] */;
    static char name_out_mass__[50];
    static integer i_time__;
    extern /* Subroutine */ int getarg_(integer *, char *, ftnlen);
    static integer n_time__, i_loop__, n_loop__;
    static char name_in__[50];

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };
    static cilist io___4 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___10 = { 0, 6, 0, 0, 0 };
    static cilist io___11 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };
    static cilist io___13 = { 0, 6, 0, 0, 0 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };
    static cilist io___15 = { 0, 20, 0, "(a10,i10)", 0 };
    static cilist io___18 = { 0, 20, 0, "(a10,i10)", 0 };
    static cilist io___20 = { 0, 20, 0, "(a10,i10)", 0 };
    static cilist io___22 = { 0, 6, 0, 0, 0 };
    static cilist io___23 = { 0, 6, 0, 0, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___25 = { 0, 6, 0, 0, 0 };
    static cilist io___26 = { 0, 6, 0, 0, 0 };
    static cilist io___27 = { 0, 6, 0, 0, 0 };
    static icilist io___29 = { 0, n_loop_str__, 0, "(i10)", 10, 1 };
    static icilist io___31 = { 0, n_time_str__, 0, "(i10)", 10, 1 };
    static cilist io___34 = { 0, 20, 0, "(a10,e14.5)", 0 };
    static cilist io___37 = { 0, 20, 0, "(i8,3e14.5)", 0 };
    static cilist io___43 = { 0, 21, 0, "(//,a10,i10)", 0 };
    static cilist io___44 = { 0, 22, 0, "(//,a10,i10)", 0 };


/* maximum number of bins */
/* maximum number of loops */
/* maximum number of times */
/* input */
/* output number */
/* output mass */
/* output number average */
/* output mass average */
/* check there is exactly one commandline argument */
    if (iargc_() != 1) {
	s_wsle(&io___1);
	do_lio(&c__9, &c__1, "Usage: process_out <filename.d>", (ftnlen)31);
	e_wsle();
	exit_(&c__2);
    }
/* get and check first commandline argument (must be "filename.d") */
    getarg_(&c__1, name_in__, (ftnlen)50);
    i__ = len_trim__(name_in__, (ftnlen)50);
    if (i__ > 40) {
	s_wsle(&io___4);
	do_lio(&c__9, &c__1, "ERROR: filename too long", (ftnlen)24);
	e_wsle();
	exit_(&c__2);
    }
    i__1 = i__ - 2;
    if (*(unsigned char *)&name_in__[i__ - 1] != 'd' || s_cmp(name_in__ + 
	    i__1, ".", i__ - 1 - i__1, (ftnlen)1) != 0) {
	s_wsle(&io___5);
	do_lio(&c__9, &c__1, "ERROR: Filename must end in .d", (ftnlen)30);
	e_wsle();
	exit_(&c__2);
    }
/* compute names of output files */
    s_copy(name_out_num__, name_in__, (ftnlen)50, (ftnlen)50);
    s_copy(name_out_mass__, name_in__, (ftnlen)50, (ftnlen)50);
    s_copy(name_out_num_avg__, name_in__, (ftnlen)50, (ftnlen)50);
    s_copy(name_out_mass_avg__, name_in__, (ftnlen)50, (ftnlen)50);
    i__1 = i__ - 2;
    s_copy(name_out_num__ + i__1, "_num.d", 50 - i__1, (ftnlen)6);
    i__1 = i__ - 2;
    s_copy(name_out_mass__ + i__1, "_mass.d", 50 - i__1, (ftnlen)7);
    i__1 = i__ - 2;
    s_copy(name_out_num_avg__ + i__1, "_num_avg.d", 50 - i__1, (ftnlen)10);
    i__1 = i__ - 2;
    s_copy(name_out_mass_avg__ + i__1, "_mass_avg.d", 50 - i__1, (ftnlen)11);
    s_wsle(&io___10);
    do_lio(&c__9, &c__1, "name_in = ", (ftnlen)10);
    do_lio(&c__9, &c__1, name_in__, (ftnlen)50);
    e_wsle();
    s_wsle(&io___11);
    do_lio(&c__9, &c__1, "name_out_num = ", (ftnlen)15);
    do_lio(&c__9, &c__1, name_out_num__, (ftnlen)50);
    e_wsle();
    s_wsle(&io___12);
    do_lio(&c__9, &c__1, "name_out_mass = ", (ftnlen)16);
    do_lio(&c__9, &c__1, name_out_mass__, (ftnlen)50);
    e_wsle();
    s_wsle(&io___13);
    do_lio(&c__9, &c__1, "name_out_num_avg = ", (ftnlen)19);
    do_lio(&c__9, &c__1, name_out_num_avg__, (ftnlen)50);
    e_wsle();
    s_wsle(&io___14);
    do_lio(&c__9, &c__1, "name_out_mass_avg = ", (ftnlen)20);
    do_lio(&c__9, &c__1, name_out_mass_avg__, (ftnlen)50);
    e_wsle();
/* open files */
    o__1.oerr = 0;
    o__1.ounit = 20;
    o__1.ofnmlen = 50;
    o__1.ofnm = name_in__;
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    o__1.oerr = 0;
    o__1.ounit = 21;
    o__1.ofnmlen = 50;
    o__1.ofnm = name_out_num__;
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    o__1.oerr = 0;
    o__1.ounit = 22;
    o__1.ofnmlen = 50;
    o__1.ofnm = name_out_mass__;
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    o__1.oerr = 0;
    o__1.ounit = 23;
    o__1.ofnmlen = 50;
    o__1.ofnm = name_out_num_avg__;
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    o__1.oerr = 0;
    o__1.ounit = 24;
    o__1.ofnmlen = 50;
    o__1.ofnm = name_out_mass_avg__;
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
/* read and check dimensions */
    s_rsfe(&io___15);
    do_fio(&c__1, dum, (ftnlen)100);
    do_fio(&c__1, (char *)&n_loop__, (ftnlen)sizeof(integer));
    e_rsfe();
    s_rsfe(&io___18);
    do_fio(&c__1, dum, (ftnlen)100);
    do_fio(&c__1, (char *)&n_bin__, (ftnlen)sizeof(integer));
    e_rsfe();
    s_rsfe(&io___20);
    do_fio(&c__1, dum, (ftnlen)100);
    do_fio(&c__1, (char *)&n_time__, (ftnlen)sizeof(integer));
    e_rsfe();
    if (n_loop__ > 50) {
	s_wsle(&io___22);
	do_lio(&c__9, &c__1, "ERROR: n_loop too large", (ftnlen)23);
	e_wsle();
	exit_(&c__2);
    }
    if (n_bin__ > 400) {
	s_wsle(&io___23);
	do_lio(&c__9, &c__1, "ERROR: n_bin too large", (ftnlen)22);
	e_wsle();
	exit_(&c__2);
    }
    if (n_time__ > 100) {
	s_wsle(&io___24);
	do_lio(&c__9, &c__1, "ERROR: n_time too large", (ftnlen)23);
	e_wsle();
	exit_(&c__2);
    }
    s_wsle(&io___25);
    do_lio(&c__9, &c__1, "n_loop = ", (ftnlen)9);
    do_lio(&c__3, &c__1, (char *)&n_loop__, (ftnlen)sizeof(integer));
    e_wsle();
    s_wsle(&io___26);
    do_lio(&c__9, &c__1, "n_bin = ", (ftnlen)8);
    do_lio(&c__3, &c__1, (char *)&n_bin__, (ftnlen)sizeof(integer));
    e_wsle();
    s_wsle(&io___27);
    do_lio(&c__9, &c__1, "n_time = ", (ftnlen)9);
    do_lio(&c__3, &c__1, (char *)&n_time__, (ftnlen)sizeof(integer));
    e_wsle();
    s_wsfi(&io___29);
    do_fio(&c__1, (char *)&n_loop__, (ftnlen)sizeof(integer));
    e_wsfi();
    s_wsfi(&io___31);
    do_fio(&c__1, (char *)&n_time__, (ftnlen)sizeof(integer));
    e_wsfi();
/* read all data */
    i__1 = n_loop__;
    for (i_loop__ = 1; i_loop__ <= i__1; ++i_loop__) {
	i__2 = n_time__;
	for (i_time__ = 1; i_time__ <= i__2; ++i_time__) {
	    s_rsfe(&io___34);
	    do_fio(&c__1, dum, (ftnlen)100);
	    do_fio(&c__1, (char *)&time, (ftnlen)sizeof(doublereal));
	    e_rsfe();
	    i__3 = n_bin__;
	    for (i_bin__ = 1; i_bin__ <= i__3; ++i_bin__) {
		s_rsfe(&io___37);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&bin_r__[i_bin__ - 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&n[i_loop__ + (i_time__ + i_bin__ * 100)
			 * 50 - 5051], (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&bin_g__[i_loop__ + (i_time__ + i_bin__ 
			* 100) * 50 - 5051], (ftnlen)sizeof(doublereal));
		e_rsfe();
	    }
	}
    }
/* compute averages */
    i__1 = n_time__;
    for (i_time__ = 1; i_time__ <= i__1; ++i_time__) {
	i__2 = n_bin__;
	for (i_bin__ = 1; i_bin__ <= i__2; ++i_bin__) {
	    g_avg__[i_time__ + i_bin__ * 100 - 101] = 0.;
	    n_avg__[i_time__ + i_bin__ * 100 - 101] = 0.;
	    i__3 = n_loop__;
	    for (i_loop__ = 1; i_loop__ <= i__3; ++i_loop__) {
		g_avg__[i_time__ + i_bin__ * 100 - 101] += bin_g__[i_loop__ + 
			(i_time__ + i_bin__ * 100) * 50 - 5051];
		n_avg__[i_time__ + i_bin__ * 100 - 101] += n[i_loop__ + (
			i_time__ + i_bin__ * 100) * 50 - 5051];
	    }
	    g_avg__[i_time__ + i_bin__ * 100 - 101] /= n_loop__;
	    n_avg__[i_time__ + i_bin__ * 100 - 101] /= n_loop__;
	}
    }
/* output data */
    i__1 = n_time__;
    for (i_time__ = 1; i_time__ <= i__1; ++i_time__) {
	s_wsfe(&io___43);
	do_fio(&c__1, "time", (ftnlen)4);
	i__2 = i_time__ - 1;
	do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___44);
	do_fio(&c__1, "time", (ftnlen)4);
	i__2 = i_time__ - 1;
	do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	e_wsfe();
	i__2 = n_bin__;
	for (i_bin__ = 1; i_bin__ <= i__2; ++i_bin__) {
	    ci__1.cierr = 0;
	    ci__1.ciunit = 21;
/* Writing concatenation */
	    i__4[0] = 10, a__1[0] = "(i8,e14.5,";
	    i__4[1] = 10, a__1[1] = n_loop_str__;
	    i__4[2] = 6, a__1[2] = "e14.5)";
	    ci__1.cifmt = (s_cat(ch__1, a__1, i__4, &c__3, (ftnlen)26), ch__1)
		    ;
	    s_wsfe(&ci__1);
	    do_fio(&c__1, (char *)&i_bin__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&bin_r__[i_bin__ - 1], (ftnlen)sizeof(
		    doublereal));
	    i__3 = n_loop__;
	    for (i_loop__ = 1; i_loop__ <= i__3; ++i_loop__) {
		do_fio(&c__1, (char *)&n[i_loop__ + (i_time__ + i_bin__ * 100)
			 * 50 - 5051], (ftnlen)sizeof(doublereal));
	    }
	    e_wsfe();
	    ci__1.cierr = 0;
	    ci__1.ciunit = 22;
/* Writing concatenation */
	    i__4[0] = 10, a__1[0] = "(i8,e14.5,";
	    i__4[1] = 10, a__1[1] = n_loop_str__;
	    i__4[2] = 6, a__1[2] = "e14.5)";
	    ci__1.cifmt = (s_cat(ch__1, a__1, i__4, &c__3, (ftnlen)26), ch__1)
		    ;
	    s_wsfe(&ci__1);
	    do_fio(&c__1, (char *)&i_bin__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&bin_r__[i_bin__ - 1], (ftnlen)sizeof(
		    doublereal));
	    i__3 = n_loop__;
	    for (i_loop__ = 1; i_loop__ <= i__3; ++i_loop__) {
		do_fio(&c__1, (char *)&bin_g__[i_loop__ + (i_time__ + i_bin__ 
			* 100) * 50 - 5051], (ftnlen)sizeof(doublereal));
	    }
	    e_wsfe();
	}
    }
    i__1 = n_bin__;
    for (i_bin__ = 1; i_bin__ <= i__1; ++i_bin__) {
	ci__1.cierr = 0;
	ci__1.ciunit = 23;
/* Writing concatenation */
	i__4[0] = 10, a__1[0] = "(i8,e14.5,";
	i__4[1] = 10, a__1[1] = n_time_str__;
	i__4[2] = 6, a__1[2] = "e14.5)";
	ci__1.cifmt = (s_cat(ch__1, a__1, i__4, &c__3, (ftnlen)26), ch__1);
	s_wsfe(&ci__1);
	do_fio(&c__1, (char *)&i_bin__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&bin_r__[i_bin__ - 1], (ftnlen)sizeof(
		doublereal));
	i__2 = n_time__;
	for (i_time__ = 1; i_time__ <= i__2; ++i_time__) {
	    do_fio(&c__1, (char *)&n_avg__[i_time__ + i_bin__ * 100 - 101], (
		    ftnlen)sizeof(doublereal));
	}
	e_wsfe();
	ci__1.cierr = 0;
	ci__1.ciunit = 24;
/* Writing concatenation */
	i__4[0] = 10, a__1[0] = "(i8,e14.5,";
	i__4[1] = 10, a__1[1] = n_time_str__;
	i__4[2] = 6, a__1[2] = "e14.5)";
	ci__1.cifmt = (s_cat(ch__1, a__1, i__4, &c__3, (ftnlen)26), ch__1);
	s_wsfe(&ci__1);
	do_fio(&c__1, (char *)&i_bin__, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&bin_r__[i_bin__ - 1], (ftnlen)sizeof(
		doublereal));
	i__2 = n_time__;
	for (i_time__ = 1; i_time__ <= i__2; ++i_time__) {
	    do_fio(&c__1, (char *)&g_avg__[i_time__ + i_bin__ * 100 - 101], (
		    ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }
    return 0;
} /* MAIN__ */

/* Main program alias */ int process_out__ () { MAIN__ (); return 0; }
