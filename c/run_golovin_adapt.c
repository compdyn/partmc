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

/* Table of constant values */

double c_b2 = 10.;
int c__10 = 10;
int c__160 = 160;
int c__3 = 3;
double c_b8 = 1e3;
int c__10000 = 10000;
double c_b13 = 600.;
double c_b14 = 60.;
double c_b15 = .01;
double c_b16 = .005;
double c_b17 = 1.;

/* Simulation with Golovin kernel and adaptive timestepping. */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<       program MonteCarlo >*/
/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    int i__1;
    double d__1;
    olist o__1;

    /* Builtin functions */
    int f_open(olist *), i_dnnt(double *);
    double exp(double);

    /* Local variables */
    extern /* Subroutine */ int mc_adapt__(int *, int *, double *,
	     double *, int *, double *, double *, double *
	    , int *, double *, U_fp, double *, double *, 
	    double *, double *, double *, int *), 
	    compute_volumes__(int *, int *, double *, double *
	    , double *, double *, int *);
    int k, m;
    double v[10000];
    extern /* Subroutine */ int make_grid__(int *, int *, double *
	    , double *, double *, double *);
    double dlnr;
    extern /* Subroutine */ int print_header__(int *, int *, int *
	    );
    double bin_g[160];
    int bin_n[160];
    double bin_r[160], n_ini__[160], bin_v[160];
    extern /* Subroutine */ int srand_(int *);
    int i_loop__;
    double v_comp__;
    extern /* Subroutine */ int kernel_golovin__();

/*<       int MM, n_bin, n_loop, scal >*/
/*<       real*8 t_max, rho_p, N_0, t_print, p_max, V_0 >*/
/*<       real*8 r_samp_max, del_t_max >*/
/*<       parameter (MM = 10000)           ! number of particles >*/
/*<       parameter (n_bin = 160)          ! number of bins >*/
/*<       parameter (n_loop = 10)          ! number of loops >*/
/*<       parameter (scal = 3)             ! scale factor for bins >*/
/*<       parameter (t_max = 600d0)        ! total simulation time (seconds) >*/
/*<       parameter (rho_p = 1000d0)       ! particle density (kg/m^3) >*/
/*<       parameter (N_0 = 1d9)            ! particle number concentration ( >*/
/*<       parameter (t_print = 60d0)       ! interval between printing (s) >*/
/*<       parameter (p_max = 0.01d0)       ! maximum coagulation probability >*/
/*<       parameter (r_samp_max = 0.005d0) ! maximum sampling ratio per time >*/
/*<       parameter (del_t_max = 1d0)      ! maximum timestep >*/
/*<       parameter (V_0 = 4.1886d-15)     ! mean volume of initial distribu >*/
/*<       int M, i_loop, k >*/
/*<       real*8 V(MM), V_comp, dlnr >*/
/*<       real*8 n_ini(n_bin), bin_v(n_bin),  bin_r(n_bin) >*/
/*<       real*8 bin_g(n_bin) >*/
/*<       int bin_n(n_bin) >*/
/*<       real*8 pi >*/
/*<       parameter (pi = 3.14159265358979323846d0) >*/
/*<       external kernel_golovin >*/
/*<       open(30,file='out_golovin_adapt.d') >*/
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
/*<       call print_header(n_loop, n_bin, nint(t_max / t_print) + 1) >*/
    i__1 = i_dnnt(&c_b2) + 1;
    print_header__(&c__10, &c__160, &i__1);
/*<       call srand(10) >*/
    srand_(&c__10);
/*<       do i_loop = 1,n_loop >*/
    for (i_loop__ = 1; i_loop__ <= 10; ++i_loop__) {
/*<          call make_grid(n_bin, scal, rho_p, bin_v, bin_r, dlnr) >*/
	make_grid__(&c__160, &c__3, &c_b8, bin_v, bin_r, &dlnr);
/* define initial exponential distribution */
/*<          do k = 1,n_bin >*/
	for (k = 1; k <= 160; ++k) {
/*<    >*/
/* Computing 3rd power */
	    d__1 = bin_r[k - 1] * 2.;
	    n_ini__[k - 1] = d__1 * (d__1 * d__1) * 1.5707963267948966 * 
		    10000 / 4.1886e-15 * exp(-(bin_v[k - 1] / 4.1886e-15));
/*<          enddo >*/
	}
/*<          call compute_volumes(n_bin, MM, n_ini, bin_r, dlnr, V, M) >*/
	compute_volumes__(&c__160, &c__10000, n_ini__, bin_r, &dlnr, v, &m);
/*<          V_comp = M / N_0 >*/
	v_comp__ = m / 1e9;
/*<    >*/
	mc_adapt__(&c__10000, &m, v, &v_comp__, &c__160, bin_v, bin_r, 
		bin_g, bin_n, &dlnr, (U_fp)kernel_golovin__, &c_b13, &
		c_b14, &c_b15, &c_b16, &c_b17, &i_loop__);
/*<       enddo >*/
    }
/*<       end >*/
    return 0;
} /* MAIN__ */

/* Main program alias */ int montecarlo_ () { MAIN__ (); return 0; }
