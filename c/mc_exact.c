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

/* Exact solution output. */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<    >*/
/* Subroutine */ int mc_exact__(int n_bin, double *bin_v, 
	double *bin_r, double *bin_g, int *bin_n, 
	double dlnr, double n_0, double v_0, double *
	rho_p__, S_fp soln, double *t_max__, double *t_print__, 
	int *loop, double v_comp)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    extern /* Subroutine */ int print_info__(double *, double *, 
	    int *, double *, double *, double *, int *, 
	    double *);
    double time;
    int i_time__, n_time__;

/* FIXME: N_0 and V_0 are really parameters for the initial value */
/* of the particle distribution. They should be replaced by a n_par */
/* params() pair. */
/*<       int n_bin        ! INPUT: number of bins >*/
/*<       real*8 bin_v(n_bin)  ! INPUT: volume of bins >*/
/*<       real*8 bin_r(n_bin)  ! INPUT: radius of bins >*/
/*<       real*8 bin_g(n_bin)  ! OUTPUT: mass in bins >*/
/*<       int bin_n(n_bin) ! OUTPUT: number in bins >*/
/*<       real*8 dlnr          ! INPUT: bin scale factor >*/
/*<       real*8 N_0           ! INPUT: particle number concentration (#/m^3 >*/
/*<       real*8 V_0           ! INPUT: >*/
/*<       real*8 rho_p         ! INPUT: particle density (kg/m^3) >*/
/*<       external soln        ! INPUT: exact solution procedure >*/
/*<       real*8 t_max         ! INPUT: total simulation time >*/
/*<       real*8 t_print       ! INPUT: interval to print info (seconds) >*/
/*<       int loop         ! INPUT: loop number of run >*/
/*<       real*8 V_comp        ! INPUT: computational volume >*/
/*<       int i_time, n_time >*/
/*<       real*8 time >*/
/*<       n_time = int(t_max / t_print) >*/
    /* Parameter adjustments */
    --bin_n;
    --bin_g;
    --bin_r;
    --bin_v;

    /* Function Body */
    n_time__ = (int) (*t_max__ / *t_print__);
/*<       do i_time = 0,n_time >*/
    i__1 = n_time__;
    for (i_time__ = 0; i_time__ <= i__1; ++i_time__) {
/*<          time = dble(i_time) / dble(n_time) * dble(t_max) >*/
	time = (double) i_time__ / (double) n_time__ * *t_max__;
/*<    >*/
	(*soln)(n_bin__, &bin_v[1], &bin_r[1], &bin_g[1], &bin_n[1], 
		dlnr, &time, n_0__, v_0, rho_p__, v_comp__);
/*<    >*/
	print_info__(&time, v_comp__, n_bin__, &bin_v[1], &bin_r[1], &
		bin_g[1], &bin_n[1], dlnr);
/*<       enddo >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* mc_exact__ */

