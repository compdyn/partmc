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

/* Table of constant values */

double c_b2 = 10.;
int c__1 = 1;
int c__160 = 160;
int c__3 = 3;
double c_b7 = 1e3;
double c_b9 = 1e9;
double c_b10 = 4.1886e-15;
double c_b12 = 600.;
double c_b13 = 60.;
double c_b14 = 1.;

/* Exact solution with Golovin kernel. */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<       program MonteCarlo >*/
/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    int i__1;
    olist o__1;

    /* Builtin functions */
    int f_open(olist *), i_dnnt(double *);

    /* Local variables */
    extern /* Subroutine */ int mc_exact__(int *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, U_fp, double *, 
	    double *, int *, double *);
    extern /* Subroutine */ int soln_golovin_exp__();
    extern /* Subroutine */ int make_grid__(int *, int *, double *
	    , double *, double *, double *);
    double dlnr;
    extern /* Subroutine */ int print_header__(int *, int *, int *
	    );
    double bin_g[160], bin_n[160], bin_r[160], bin_v[160];
    int i_loop__;

/*<       int n_bin, n_loop, scal >*/
/*<       real*8 t_max, rho_p, N_0, t_print, V_0, V_comp >*/
/*<       parameter (n_bin = 160)        ! number of bins >*/
/*<       parameter (n_loop = 1)         ! number of loops >*/
/*<       parameter (scal = 3)           ! scale factor for bins >*/
/*<       parameter (t_max = 600d0)      ! total simulation time (seconds) >*/
/*<       parameter (rho_p = 1000d0)     ! particle density (kg/m^3) >*/
/*<       parameter (N_0 = 1d9)          ! particle number concentration (#/ >*/
/*<       parameter (t_print = 60d0)     ! interval between printing (s) >*/
/*<       parameter (V_0 = 4.1886d-15)   ! mean volume of initial distributi >*/
/*<       parameter (V_comp = 1d0)       ! computational volume (dummy value >*/
/*<       int i_loop >*/
/*<       real*8 dlnr >*/
/*<       real*8 bin_v(n_bin), bin_r(n_bin), bin_g(n_bin), bin_n(n_bin) >*/
/*<       external soln_golovin_exp >*/
/*<       open(30,file='out_golovin_exact.d') >*/
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
/*<       call print_header(n_loop, n_bin, nint(t_max / t_print) + 1) >*/
    i__1 = i_dnnt(&c_b2) + 1;
    print_header__(&c__1, &c__160, &i__1);
/*<       do i_loop = 1,n_loop >*/
    for (i_loop__ = 1; i_loop__ <= 1; ++i_loop__) {
/*<          call make_grid(n_bin, scal, rho_p, bin_v, bin_r, dlnr) >*/
	make_grid__(&c__160, &c__3, &c_b7, bin_v, bin_r, &dlnr);
/*<    >*/
	mc_exact__(&c__160, bin_v, bin_r, bin_g, bin_n, &dlnr, &c_b9, 
		&c_b10, &c_b7, (U_fp)soln_golovin_exp__, &c_b12, &c_b13, &
		i_loop__, &c_b14);
/*<       enddo >*/
    }
/*<       end >*/
    return 0;
} /* MAIN__ */

/* Main program alias */ int montecarlo_ () { MAIN__ (); return 0; }
