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

int c__9 = 9;
int c__1 = 1;
int c__2 = 2;
int c__3 = 3;

/* Process output data files. */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<       program process_out >*/
/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    address a__1[3];
    int i__1, i__2, i__3, i__4[3];
    char ch__1[26];
    cilist ci__1;
    olist o__1;

    /* Builtin functions */
    int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen), 
	    e_wsle(void), s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    int f_open(olist *), s_rsfe(cilist *), do_fio(int *, char *, 
	    ftnlen), e_rsfe(void), s_wsfi(icilist *), e_wsfi(void), s_wsfe(
	    cilist *), e_wsfe(void);
    /* Subroutine */ int s_cat(char *, char **, int *, int *, ftnlen);

    /* Local variables */
    extern int len_trim__(char *, ftnlen);
    char name_out_num_avg__[50];
    int i__;
    double n[2000000]	/* was [50][100][400] */;
    char name_out_mass_avg__[50], n_time_str__[10], n_loop_str__[10], 
	    dum[100];
    double time;
    extern /* Subroutine */ int exit_(int *);
    char name_out_num__[50];
    double bin_g[2000000]	/* was [50][100][400] */;
    int i_bin__;
    double g_avg__[40000]	/* was [100][400] */;
    int n_bin__;
    extern int iargc_(void);
    double bin_r[400], n_avg__[40000]	/* was [100][400] */;
    char name_out_mass__[50];
    int i_time__;
    extern /* Subroutine */ int getarg_(int *, char *, ftnlen);
    int n_time__, i_loop__, n_loop__;
    char name_in__[50];

    /* Fortran I/O blocks */
    cilist io___1 = { 0, 6, 0, 0, 0 };
    cilist io___4 = { 0, 6, 0, 0, 0 };
    cilist io___5 = { 0, 6, 0, 0, 0 };
    cilist io___10 = { 0, 6, 0, 0, 0 };
    cilist io___11 = { 0, 6, 0, 0, 0 };
    cilist io___12 = { 0, 6, 0, 0, 0 };
    cilist io___13 = { 0, 6, 0, 0, 0 };
    cilist io___14 = { 0, 6, 0, 0, 0 };
    cilist io___15 = { 0, 20, 0, "(a10,i10)", 0 };
    cilist io___18 = { 0, 20, 0, "(a10,i10)", 0 };
    cilist io___20 = { 0, 20, 0, "(a10,i10)", 0 };
    cilist io___22 = { 0, 6, 0, 0, 0 };
    cilist io___23 = { 0, 6, 0, 0, 0 };
    cilist io___24 = { 0, 6, 0, 0, 0 };
    cilist io___25 = { 0, 6, 0, 0, 0 };
    cilist io___26 = { 0, 6, 0, 0, 0 };
    cilist io___27 = { 0, 6, 0, 0, 0 };
    icilist io___29 = { 0, n_loop_str__, 0, "(i10)", 10, 1 };
    icilist io___31 = { 0, n_time_str__, 0, "(i10)", 10, 1 };
    cilist io___34 = { 0, 20, 0, "(a10,e14.5)", 0 };
    cilist io___37 = { 0, 20, 0, "(i8,3e14.5)", 0 };
    cilist io___43 = { 0, 21, 0, "(//,a10,i10)", 0 };
    cilist io___44 = { 0, 22, 0, "(//,a10,i10)", 0 };


/*<       int n_bin_max, n_loop_max, n_time_max >*/
/*<       parameter (n_bin_max = 400)  ! maximum number of bins >*/
/*<       parameter (n_loop_max = 50)  ! maximum number of loops >*/
/*<       parameter (n_time_max = 100) ! maximum number of times >*/
/*<       int f_in, f_out_num, f_out_mass >*/
/*<       int f_out_num_avg, f_out_mass_avg >*/
/*<       parameter (f_in = 20)           ! input >*/
/*<       parameter (f_out_num = 21)      ! output number >*/
/*<       parameter (f_out_mass = 22)     ! output mass >*/
/*<       parameter (f_out_num_avg = 23)  ! output number average >*/
/*<       parameter (f_out_mass_avg = 24) ! output mass average >*/
/*<       int n_bin, n_loop, n_time >*/
/*<       character name_in*50, name_out_num*50, name_out_mass*50 >*/
/*<       character name_out_num_avg*50, name_out_mass_avg*50 >*/
/*<       character dum*100, n_loop_str*10, n_time_str*10 >*/
/*<       real*8 bin_r(n_bin_max), bin_g(n_loop_max, n_time_max, n_bin_max) >*/
/*<       real*8 n(n_loop_max, n_time_max, n_bin_max) >*/
/*<       real*8 g_avg(n_time_max, n_bin_max) >*/
/*<       real*8 n_avg(n_time_max, n_bin_max) >*/
/*<       real*8 time >*/
/*<       int i, i_loop, i_time, i_bin >*/
/* check there is exactly one commandline argument */
/*<       if (iargc() .ne. 1) then >*/
    if (iargc_() != 1) {
/*<          write(6,*) 'Usage: process_out <filename.d>' >*/
	s_wsle(&io___1);
	do_lio(&c__9, &c__1, "Usage: process_out <filename.d>", (ftnlen)31);
	e_wsle();
/*<          call exit(2) >*/
	exit_(&c__2);
/*<       endif >*/
    }
/* get and check first commandline argument (must be "filename.d") */
/*<       call getarg(1, name_in) >*/
    getarg_(&c__1, name_in__, (ftnlen)50);
/*<       i = len_trim(name_in) >*/
    i__ = len_trim__(name_in__, (ftnlen)50);
/*<       if (i .gt. 40) then >*/
    if (i__ > 40) {
/*<          write(6,*) 'ERROR: filename too long' >*/
	s_wsle(&io___4);
	do_lio(&c__9, &c__1, "ERROR: filename too long", (ftnlen)24);
	e_wsle();
/*<          call exit(2) >*/
	exit_(&c__2);
/*<       endif >*/
    }
/*<    >*/
    i__1 = i__ - 2;
    if (*(unsigned char *)&name_in__[i__ - 1] != 'd' || s_cmp(name_in__ + 
	    i__1, ".", i__ - 1 - i__1, (ftnlen)1) != 0) {
/*<          write(6,*) 'ERROR: Filename must end in .d' >*/
	s_wsle(&io___5);
	do_lio(&c__9, &c__1, "ERROR: Filename must end in .d", (ftnlen)30);
	e_wsle();
/*<          call exit(2) >*/
	exit_(&c__2);
/*<       endif >*/
    }
/* compute names of output files */
/*<       name_out_num = name_in >*/
    s_copy(name_out_num__, name_in__, (ftnlen)50, (ftnlen)50);
/*<       name_out_mass = name_in >*/
    s_copy(name_out_mass__, name_in__, (ftnlen)50, (ftnlen)50);
/*<       name_out_num_avg = name_in >*/
    s_copy(name_out_num_avg__, name_in__, (ftnlen)50, (ftnlen)50);
/*<       name_out_mass_avg = name_in >*/
    s_copy(name_out_mass_avg__, name_in__, (ftnlen)50, (ftnlen)50);
/*<       name_out_num((i-1):) = '_num.d' >*/
    i__1 = i__ - 2;
    s_copy(name_out_num__ + i__1, "_num.d", 50 - i__1, (ftnlen)6);
/*<       name_out_mass((i-1):) = '_mass.d' >*/
    i__1 = i__ - 2;
    s_copy(name_out_mass__ + i__1, "_mass.d", 50 - i__1, (ftnlen)7);
/*<       name_out_num_avg((i-1):) = '_num_avg.d' >*/
    i__1 = i__ - 2;
    s_copy(name_out_num_avg__ + i__1, "_num_avg.d", 50 - i__1, (ftnlen)10);
/*<       name_out_mass_avg((i-1):) = '_mass_avg.d' >*/
    i__1 = i__ - 2;
    s_copy(name_out_mass_avg__ + i__1, "_mass_avg.d", 50 - i__1, (ftnlen)11);
/*<       write(6,*) 'name_in = ', name_in >*/
    s_wsle(&io___10);
    do_lio(&c__9, &c__1, "name_in = ", (ftnlen)10);
    do_lio(&c__9, &c__1, name_in__, (ftnlen)50);
    e_wsle();
/*<       write(6,*) 'name_out_num = ', name_out_num >*/
    s_wsle(&io___11);
    do_lio(&c__9, &c__1, "name_out_num = ", (ftnlen)15);
    do_lio(&c__9, &c__1, name_out_num__, (ftnlen)50);
    e_wsle();
/*<       write(6,*) 'name_out_mass = ', name_out_mass >*/
    s_wsle(&io___12);
    do_lio(&c__9, &c__1, "name_out_mass = ", (ftnlen)16);
    do_lio(&c__9, &c__1, name_out_mass__, (ftnlen)50);
    e_wsle();
/*<       write(6,*) 'name_out_num_avg = ', name_out_num_avg >*/
    s_wsle(&io___13);
    do_lio(&c__9, &c__1, "name_out_num_avg = ", (ftnlen)19);
    do_lio(&c__9, &c__1, name_out_num_avg__, (ftnlen)50);
    e_wsle();
/*<       write(6,*) 'name_out_mass_avg = ', name_out_mass_avg >*/
    s_wsle(&io___14);
    do_lio(&c__9, &c__1, "name_out_mass_avg = ", (ftnlen)20);
    do_lio(&c__9, &c__1, name_out_mass_avg__, (ftnlen)50);
    e_wsle();
/* open files */
/*<       open(f_in, file=name_in) >*/
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
/*<       open(f_out_num, file=name_out_num) >*/
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
/*<       open(f_out_mass, file=name_out_mass) >*/
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
/*<       open(f_out_num_avg, file=name_out_num_avg) >*/
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
/*<       open(f_out_mass_avg, file=name_out_mass_avg) >*/
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
/*<       read(f_in, '(a10,i10)') dum, n_loop >*/
    s_rsfe(&io___15);
    do_fio(&c__1, dum, (ftnlen)100);
    do_fio(&c__1, (char *)&n_loop__, (ftnlen)sizeof(int));
    e_rsfe();
/*<       read(f_in, '(a10,i10)') dum, n_bin >*/
    s_rsfe(&io___18);
    do_fio(&c__1, dum, (ftnlen)100);
    do_fio(&c__1, (char *)&n_bin__, (ftnlen)sizeof(int));
    e_rsfe();
/*<       read(f_in, '(a10,i10)') dum, n_time >*/
    s_rsfe(&io___20);
    do_fio(&c__1, dum, (ftnlen)100);
    do_fio(&c__1, (char *)&n_time__, (ftnlen)sizeof(int));
    e_rsfe();
/*<       if (n_loop .gt. n_loop_max) then >*/
    if (n_loop__ > 50) {
/*<          write(6,*) 'ERROR: n_loop too large' >*/
	s_wsle(&io___22);
	do_lio(&c__9, &c__1, "ERROR: n_loop too large", (ftnlen)23);
	e_wsle();
/*<          call exit(2) >*/
	exit_(&c__2);
/*<       endif >*/
    }
/*<       if (n_bin .gt. n_bin_max) then >*/
    if (n_bin__ > 400) {
/*<          write(6,*) 'ERROR: n_bin too large' >*/
	s_wsle(&io___23);
	do_lio(&c__9, &c__1, "ERROR: n_bin too large", (ftnlen)22);
	e_wsle();
/*<          call exit(2) >*/
	exit_(&c__2);
/*<       endif >*/
    }
/*<       if (n_time .gt. n_time_max) then >*/
    if (n_time__ > 100) {
/*<          write(6,*) 'ERROR: n_time too large' >*/
	s_wsle(&io___24);
	do_lio(&c__9, &c__1, "ERROR: n_time too large", (ftnlen)23);
	e_wsle();
/*<          call exit(2) >*/
	exit_(&c__2);
/*<       endif >*/
    }
/*<       write(6,*) 'n_loop = ', n_loop >*/
    s_wsle(&io___25);
    do_lio(&c__9, &c__1, "n_loop = ", (ftnlen)9);
    do_lio(&c__3, &c__1, (char *)&n_loop__, (ftnlen)sizeof(int));
    e_wsle();
/*<       write(6,*) 'n_bin = ', n_bin >*/
    s_wsle(&io___26);
    do_lio(&c__9, &c__1, "n_bin = ", (ftnlen)8);
    do_lio(&c__3, &c__1, (char *)&n_bin__, (ftnlen)sizeof(int));
    e_wsle();
/*<       write(6,*) 'n_time = ', n_time >*/
    s_wsle(&io___27);
    do_lio(&c__9, &c__1, "n_time = ", (ftnlen)9);
    do_lio(&c__3, &c__1, (char *)&n_time__, (ftnlen)sizeof(int));
    e_wsle();
/*<       write(n_loop_str, '(i10)') n_loop >*/
    s_wsfi(&io___29);
    do_fio(&c__1, (char *)&n_loop__, (ftnlen)sizeof(int));
    e_wsfi();
/*<       write(n_time_str, '(i10)') n_time >*/
    s_wsfi(&io___31);
    do_fio(&c__1, (char *)&n_time__, (ftnlen)sizeof(int));
    e_wsfi();
/* read all data */
/*<       do i_loop = 1,n_loop >*/
    i__1 = n_loop__;
    for (i_loop__ = 1; i_loop__ <= i__1; ++i_loop__) {
/*<          do i_time = 1,n_time >*/
	i__2 = n_time__;
	for (i_time__ = 1; i_time__ <= i__2; ++i_time__) {
/*<             read(f_in, '(a10,e14.5)') dum, time >*/
	    s_rsfe(&io___34);
	    do_fio(&c__1, dum, (ftnlen)100);
	    do_fio(&c__1, (char *)&time, (ftnlen)sizeof(double));
	    e_rsfe();
/*<             do i_bin = 1,n_bin >*/
	    i__3 = n_bin__;
	    for (i_bin__ = 1; i_bin__ <= i__3; ++i_bin__) {
/*<    >*/
		s_rsfe(&io___37);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&bin_r[i_bin__ - 1], (ftnlen)sizeof(
			double));
		do_fio(&c__1, (char *)&n[i_loop__ + (i_time__ + i_bin__ * 100)
			 * 50 - 5051], (ftnlen)sizeof(double));
		do_fio(&c__1, (char *)&bin_g[i_loop__ + (i_time__ + i_bin__ 
			* 100) * 50 - 5051], (ftnlen)sizeof(double));
		e_rsfe();
/*<             enddo >*/
	    }
/*<          enddo >*/
	}
/*<       enddo >*/
    }
/* compute averages */
/*<       do i_time = 1,n_time >*/
    i__1 = n_time__;
    for (i_time__ = 1; i_time__ <= i__1; ++i_time__) {
/*<          do i_bin = 1,n_bin >*/
	i__2 = n_bin__;
	for (i_bin__ = 1; i_bin__ <= i__2; ++i_bin__) {
/*<             g_avg(i_time, i_bin) = 0d0 >*/
	    g_avg__[i_time__ + i_bin__ * 100 - 101] = 0.;
/*<             n_avg(i_time, i_bin) = 0d0 >*/
	    n_avg__[i_time__ + i_bin__ * 100 - 101] = 0.;
/*<             do i_loop = 1,n_loop >*/
	    i__3 = n_loop__;
	    for (i_loop__ = 1; i_loop__ <= i__3; ++i_loop__) {
/*<    >*/
		g_avg__[i_time__ + i_bin__ * 100 - 101] += bin_g[i_loop__ + 
			(i_time__ + i_bin__ * 100) * 50 - 5051];
/*<    >*/
		n_avg__[i_time__ + i_bin__ * 100 - 101] += n[i_loop__ + (
			i_time__ + i_bin__ * 100) * 50 - 5051];
/*<             enddo >*/
	    }
/*<             g_avg(i_time, i_bin) = g_avg(i_time, i_bin) / n_loop >*/
	    g_avg__[i_time__ + i_bin__ * 100 - 101] /= n_loop__;
/*<             n_avg(i_time, i_bin) = n_avg(i_time, i_bin) / n_loop >*/
	    n_avg__[i_time__ + i_bin__ * 100 - 101] /= n_loop__;
/*<          enddo >*/
	}
/*<       enddo >*/
    }
/* output data */
/*<       do i_time = 1,n_time >*/
    i__1 = n_time__;
    for (i_time__ = 1; i_time__ <= i__1; ++i_time__) {
/*<          write(f_out_num, '(//,a10,i10)') 'time', i_time - 1 >*/
	s_wsfe(&io___43);
	do_fio(&c__1, "time", (ftnlen)4);
	i__2 = i_time__ - 1;
	do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(int));
	e_wsfe();
/*<          write(f_out_mass, '(//,a10,i10)') 'time', i_time - 1 >*/
	s_wsfe(&io___44);
	do_fio(&c__1, "time", (ftnlen)4);
	i__2 = i_time__ - 1;
	do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(int));
	e_wsfe();
/*<          do i_bin = 1,n_bin >*/
	i__2 = n_bin__;
	for (i_bin__ = 1; i_bin__ <= i__2; ++i_bin__) {
/*<    >*/
	    ci__1.cierr = 0;
	    ci__1.ciunit = 21;
/* Writing concatenation */
	    i__4[0] = 10, a__1[0] = "(i8,e14.5,";
	    i__4[1] = 10, a__1[1] = n_loop_str__;
	    i__4[2] = 6, a__1[2] = "e14.5)";
	    ci__1.cifmt = (s_cat(ch__1, a__1, i__4, &c__3, (ftnlen)26), ch__1)
		    ;
	    s_wsfe(&ci__1);
	    do_fio(&c__1, (char *)&i_bin__, (ftnlen)sizeof(int));
	    do_fio(&c__1, (char *)&bin_r[i_bin__ - 1], (ftnlen)sizeof(
		    double));
	    i__3 = n_loop__;
	    for (i_loop__ = 1; i_loop__ <= i__3; ++i_loop__) {
		do_fio(&c__1, (char *)&n[i_loop__ + (i_time__ + i_bin__ * 100)
			 * 50 - 5051], (ftnlen)sizeof(double));
	    }
	    e_wsfe();
/*<    >*/
	    ci__1.cierr = 0;
	    ci__1.ciunit = 22;
/* Writing concatenation */
	    i__4[0] = 10, a__1[0] = "(i8,e14.5,";
	    i__4[1] = 10, a__1[1] = n_loop_str__;
	    i__4[2] = 6, a__1[2] = "e14.5)";
	    ci__1.cifmt = (s_cat(ch__1, a__1, i__4, &c__3, (ftnlen)26), ch__1)
		    ;
	    s_wsfe(&ci__1);
	    do_fio(&c__1, (char *)&i_bin__, (ftnlen)sizeof(int));
	    do_fio(&c__1, (char *)&bin_r[i_bin__ - 1], (ftnlen)sizeof(
		    double));
	    i__3 = n_loop__;
	    for (i_loop__ = 1; i_loop__ <= i__3; ++i_loop__) {
		do_fio(&c__1, (char *)&bin_g[i_loop__ + (i_time__ + i_bin__ 
			* 100) * 50 - 5051], (ftnlen)sizeof(double));
	    }
	    e_wsfe();
/*<          enddo >*/
	}
/*<       enddo >*/
    }
/*<       do i_bin = 1,n_bin >*/
    i__1 = n_bin__;
    for (i_bin__ = 1; i_bin__ <= i__1; ++i_bin__) {
/*<    >*/
	ci__1.cierr = 0;
	ci__1.ciunit = 23;
/* Writing concatenation */
	i__4[0] = 10, a__1[0] = "(i8,e14.5,";
	i__4[1] = 10, a__1[1] = n_time_str__;
	i__4[2] = 6, a__1[2] = "e14.5)";
	ci__1.cifmt = (s_cat(ch__1, a__1, i__4, &c__3, (ftnlen)26), ch__1);
	s_wsfe(&ci__1);
	do_fio(&c__1, (char *)&i_bin__, (ftnlen)sizeof(int));
	do_fio(&c__1, (char *)&bin_r[i_bin__ - 1], (ftnlen)sizeof(
		double));
	i__2 = n_time__;
	for (i_time__ = 1; i_time__ <= i__2; ++i_time__) {
	    do_fio(&c__1, (char *)&n_avg__[i_time__ + i_bin__ * 100 - 101], (
		    ftnlen)sizeof(double));
	}
	e_wsfe();
/*<    >*/
	ci__1.cierr = 0;
	ci__1.ciunit = 24;
/* Writing concatenation */
	i__4[0] = 10, a__1[0] = "(i8,e14.5,";
	i__4[1] = 10, a__1[1] = n_time_str__;
	i__4[2] = 6, a__1[2] = "e14.5)";
	ci__1.cifmt = (s_cat(ch__1, a__1, i__4, &c__3, (ftnlen)26), ch__1);
	s_wsfe(&ci__1);
	do_fio(&c__1, (char *)&i_bin__, (ftnlen)sizeof(int));
	do_fio(&c__1, (char *)&bin_r[i_bin__ - 1], (ftnlen)sizeof(
		double));
	i__2 = n_time__;
	for (i_time__ = 1; i_time__ <= i__2; ++i_time__) {
	    do_fio(&c__1, (char *)&g_avg__[i_time__ + i_bin__ * 100 - 101], (
		    ftnlen)sizeof(double));
	}
	e_wsfe();
/*<       enddo >*/
    }
/*<       end >*/
    return 0;
} /* MAIN__ */

/* Main program alias */ int process_out__ () { MAIN__ (); return 0; }
