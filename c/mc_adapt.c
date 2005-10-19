/* Monte Carlo with adaptive timestep. */

/***************************************************************************/
void mc_adapt__(int *mm, int *m, double *v, 
	double v_comp, int n_bin, double *bin_v, 
	double *bin_r, double *bin_g, int *bin_n, 
	double dlnr, U_fp kernel, double *t_max__, double *
	t_print__, double *p_max__, double *r_samp_max__, double *
	del_t_max__, int *loop)
{
    /* System generated locals */
    int i__1;

    /* Builtin functions */
    int s_wsfe(cilist *), do_fio(int *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    logical did_coag__;
    double last_print_time__;
    extern /* Subroutine */ int cpu_time__(double *);
    logical do_print__;
    extern /* Subroutine */ int est_k_max__(int *, double *, int *
	    , U_fp, double *);
    logical bin_change__;
    double t_per_samp__;
    extern /* Subroutine */ int print_info__(double *, double *, 
	    int *, double *, double *, double *, int *, 
	    double *), check_event__(double *, double *, 
	    double *, logical *);
    double time;
    extern /* Subroutine */ int compute_n_samp_del_t__(int *, double *
	    , double *, double *, double *, double *, int 
	    *, double *);
    double del_t__, t_end__, k_max__;
    int n_coag__, i_samp__;
    extern /* Subroutine */ int double_(int *, int *, double *, 
	    double *, int *, double *, double *, double *,
	     int *, double *);
    int n_samp__;
    double t_loop__;
    extern /* Subroutine */ int maybe_coag_pair__(int *, int *, 
	    double *, double *, int *, double *, double *,
	     double *, int *, double *, double *, int *, 
	    U_fp, logical *, logical *);
    double t_start__;
    extern /* Subroutine */ int moments_(int *, int *, double *, 
	    double *, int *, double *, double *, double *,
	     int *, double *);

    /* Fortran I/O blocks */
    cilist io___15 = { 0, 6, 0, "(a6,a6,a6,a6,a10,a9,a11,a9)", 0 };
    cilist io___16 = { 0, 6, 0, "(i6,f6.1,f6.3,i6,e10.3,i9,e11.3,i9)", 
	    0 };


/*<       int MM           ! INPUT: physical dimension of V >*/
/*<       int M            ! INPUT/OUTPUT: logical dimension of V >*/
/*<       real*8 V(MM)         ! INPUT/OUTPUT: particle volumes >*/
/*<       real*8 V_comp        ! INPUT/OUTPUT: computational volume >*/
/*<       int n_bin        ! INPUT: number of bins >*/
/*<       real*8 bin_v(n_bin)  ! INPUT: volume of particles in bins >*/
/*<       real*8 bin_r(n_bin)  ! INPUT: radius of particles in bins >*/
/*<       real*8 bin_g(n_bin)  ! OUTPUT: mass in bins >*/
/*<       int bin_n(n_bin) ! OUTPUT: number in bins >*/
/*<       real*8 dlnr          ! INPUT: bin scale factor >*/
/*<       external kernel      ! INPUT: kernel function >*/
/*<       real*8 t_max         ! INPUT: final time (seconds) >*/
/*<       real*8 t_print       ! INPUT: interval to print info (seconds) >*/
/*<       real*8 p_max         ! INPUT: maximum coagulation probability >*/
/*<       real*8 r_samp_max    ! INPUT: maximum sampling ratio per timestep >*/
/*<       real*8 del_t_max     ! INPUT: maximum timestep >*/
/*<       int loop         ! INPUT: loop number of run >*/
/*<       real*8 del_t, time, last_print_time, k_max >*/
/*<       int n_samp, i_samp, n_coag >*/
/*<       logical do_print, did_coag, bin_change >*/
/*<       real*8 t_start, t_end, t_loop, t_per_samp >*/
/*<       time = 0 >*/
    /* Parameter adjustments */
    --v;
    --bin_n;
    --bin_g;
    --bin_r;
    --bin_v;

    /* Function Body */
    time = 0.;
/*<       n_coag = 0 >*/
    n_coag__ = 0;
/*<    >*/
    moments_(mm, m, &v[1], v_comp__, n_bin__, &bin_v[1], &bin_r[1], &
	    bin_g[1], &bin_n[1], dlnr);
/*<       call check_event(time, t_print, last_print_time, do_print) >*/
    check_event__(&time, t_print__, &last_print_time__, &do_print__);
/*<    >*/
    if (do_print__) {
	print_info__(&time, v_comp__, n_bin__, &bin_v[1], &bin_r[1], &
		bin_g[1], &bin_n[1], dlnr);
    }
/*<       call est_k_max(n_bin, bin_v, bin_n, kernel, k_max) >*/
    est_k_max__(n_bin__, &bin_v[1], &bin_n[1], (U_fp)kernel, &k_max__);
/*<       do while (time < t_max) >*/
    while(time < *t_max__) {
/*<          call cpu_time(t_start) >*/
	cpu_time__(&t_start__);
/*<    >*/
	compute_n_samp_del_t__(m, v_comp__, &k_max__, p_max__, r_samp_max__, 
		del_t_max__, &n_samp__, &del_t__);
/*<          do i_samp = 1,n_samp >*/
	i__1 = n_samp__;
	for (i_samp__ = 1; i_samp__ <= i__1; ++i_samp__) {
/*<    >*/
	    maybe_coag_pair__(mm, m, &v[1], v_comp__, n_bin__, &bin_v[1], &
		    bin_r[1], &bin_g[1], &bin_n[1], dlnr, &del_t__, &
		    n_samp__, (U_fp)kernel, &did_coag__, &bin_change__);
/*<             if (did_coag) n_coag = n_coag + 1 >*/
	    if (did_coag__) {
		++n_coag__;
	    }
/*<    >*/
	    if (bin_change__) {
		est_k_max__(n_bin__, &bin_v[1], &bin_n[1], (U_fp)kernel, &
			k_max__);
	    }
/*<             if (M .lt. MM / 2) then >*/
	    if (*m < *mm / 2) {
/*<    >*/
		double_(mm, m, &v[1], v_comp__, n_bin__, &bin_v[1], &
			bin_r[1], &bin_g[1], &bin_n[1], dlnr);
/*<             endif >*/
	    }
/*<          enddo >*/
	}
/*<          time = time + del_t >*/
	time += del_t__;
/*<          call check_event(time, t_print, last_print_time, do_print) >*/
	check_event__(&time, t_print__, &last_print_time__, &do_print__);
/*<    >*/
	if (do_print__) {
	    print_info__(&time, v_comp__, n_bin__, &bin_v[1], &bin_r[1], &
		    bin_g[1], &bin_n[1], dlnr);
	}
/*<          call cpu_time(t_end) >*/
	cpu_time__(&t_end__);
/*<          t_loop = t_end - t_start >*/
	t_loop__ = t_end__ - t_start__;
/*<          t_per_samp = t_loop / n_samp >*/
	t_per_samp__ = t_loop__ / n_samp__;
/*<    >*/
	s_wsfe(&io___15);
	do_fio(&c__1, "loop", (ftnlen)4);
	do_fio(&c__1, "time", (ftnlen)4);
	do_fio(&c__1, "del_t", (ftnlen)5);
	do_fio(&c__1, "M", (ftnlen)1);
	do_fio(&c__1, "k_max", (ftnlen)5);
	do_fio(&c__1, "n_samp", (ftnlen)6);
	do_fio(&c__1, "t_per_samp", (ftnlen)10);
	do_fio(&c__1, "n_coag", (ftnlen)6);
	e_wsfe();
/*<    >*/
	s_wsfe(&io___16);
	do_fio(&c__1, (char *)&(*loop), (ftnlen)sizeof(int));
	do_fio(&c__1, (char *)&time, (ftnlen)sizeof(double));
	do_fio(&c__1, (char *)&del_t__, (ftnlen)sizeof(double));
	do_fio(&c__1, (char *)&(*m), (ftnlen)sizeof(int));
	do_fio(&c__1, (char *)&k_max__, (ftnlen)sizeof(double));
	do_fio(&c__1, (char *)&n_samp__, (ftnlen)sizeof(int));
	do_fio(&c__1, (char *)&t_per_samp__, (ftnlen)sizeof(double));
	do_fio(&c__1, (char *)&n_coag__, (ftnlen)sizeof(int));
	e_wsfe();
/*<       enddo >*/
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* mc_adapt__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<    >*/
/* Subroutine */ int compute_n_samp_del_t__(int *m, double v_comp, 
	double *k_max__, double *p_max__, double *r_samp_max__, 
	double *del_t_max__, int *n_samp__, double *del_t__)
{
    /* Builtin functions */
    double log(double);

    /* Local variables */
    double c__, r_samp__;

/*<       int M            ! INPUT: number of particles >*/
/*<       real*8 V_comp        ! INPUT: computational volume >*/
/*<       real*8 k_max         ! INPUT: maximum kernel value >*/
/*<       real*8 p_max         ! INPUT: maximum coagulation probability >*/
/*<       real*8 r_samp_max    ! INPUT: maximum sampling ratio per timestep >*/
/*<       real*8 del_t_max     ! INPUT: maximum timestep >*/
/*<       int n_samp       ! OUTPUT: number of samples per timestep >*/
/*<       real*8 del_t         ! OUTPUT: timestep >*/
/*<       real*8 r_samp, c >*/
/*<       c = - (k_max * 1d0/V_comp / log(1 - p_max)) >*/
    c__ = -(*k_max__ * 1. / v_comp / log(1 - *p_max__));
/*<       del_t = r_samp_max / c >*/
    *del_t__ = *r_samp_max__ / c__;
/*<       if (del_t .gt. del_t_max) del_t = del_t_max >*/
    if (*del_t__ > *del_t_max__) {
	*del_t__ = *del_t_max__;
    }
/*<       r_samp = del_t * c >*/
    r_samp__ = *del_t__ * c__;
/*<       n_samp = int(r_samp * M*(M-1)/2) >*/
    *n_samp__ = (int) (r_samp__ * *m * (*m - 1) / 2);
/*<       return >*/
    return 0;
/*<       end >*/
} /* compute_n_samp_del_t__ */

