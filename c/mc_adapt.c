/* Monte Carlo with adaptive timestep. */

/***************************************************************************/

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
void mc_adapt__(int MM, int M, double V, 
	double v_comp, int n_bin, double *bin_v, 
	double *bin_r, double *bin_g, int *bin_n, 
	double dlnr, U_fp kernel, double t_max, double *
	t_print, double p_max, double r_samp_max, double *
	del_t_max, int loop)
{
	int did_coag, do_print, bin_change;
	double last_print_time, t_per_samp, time,  del_t, t_end, k_max;
	int n_coag, i_samp, n_samp;
	double t_loop, t_start;

	/* Function Body */
	time = 0.;
/*<       n_coag = 0 >*/
	n_coag = 0;
/*<    >*/
	moments_(mm, m, &v[1], v_comp__, n_bin__, &bin_v[1], &bin_r[1], &
		 bin_g[1], &bin_n[1], dlnr);
/*<       call check_event(time, t_print, last_print_time, do_print) >*/
	check_event__(&time, t_print, &last_print_time, &do_print);
/*<    >*/
	if (do_print) {
		print_info__(&time, v_comp__, n_bin__, &bin_v[1], &bin_r[1], &
			     bin_g[1], &bin_n[1], dlnr);
	}
/*<       call est_k_max(n_bin, bin_v, bin_n, kernel, k_max) >*/
	est_k_max(n_bin__, &bin_v[1], &bin_n[1], (U_fp)kernel, &k_max);
/*<       do while (time < t_max) >*/
	while(time < t_max) {
/*<          call cpu_time(t_start) >*/
		cpu_time__(&t_start);
/*<    >*/
		compute_n_samp_del_t(m, v_comp__, &k_max, p_max__, r_samp_max__, 
				     del_t_max, &n_samp, &del_t);
/*<          do i_samp = 1,n_samp >*/
		i__1 = n_samp;
		for (i_samp = 1; i_samp <= i__1; ++i_samp) {
/*<    >*/
			maybe_coag_pair__(mm, m, &v[1], v_comp__, n_bin__, &bin_v[1], &
					  bin_r[1], &bin_g[1], &bin_n[1], dlnr, &del_t, &
					  n_samp, (U_fp)kernel, &did_coag, &bin_change);
/*<             if (did_coag) n_coag = n_coag + 1 >*/
			if (did_coag) {
				++n_coag;
			}
/*<    >*/
			if (bin_change) {
				est_k_max(n_bin__, &bin_v[1], &bin_n[1], (U_fp)kernel, &
					  k_max);
			}
/*<             if (M .lt. MM / 2) then >*/
			if (M < MM / 2) {
/*<    >*/
				double_(mm, m, &v[1], v_comp__, n_bin__, &bin_v[1], &
					bin_r[1], &bin_g[1], &bin_n[1], dlnr);
/*<             endif >*/
			}
/*<          enddo >*/
		}
/*<          time = time + del_t >*/
		time += del_t;
/*<          call check_event(time, t_print, last_print_time, do_print) >*/
		check_event__(&time, t_print, &last_print_time, &do_print);
/*<    >*/
		if (do_print) {
			print_info__(&time, v_comp__, n_bin__, &bin_v[1], &bin_r[1], &
				     bin_g[1], &bin_n[1], dlnr);
		}
/*<          call cpu_time(t_end) >*/
		cpu_time__(&t_end);
/*<          t_loop = t_end - t_start >*/
		t_loop = t_end - t_start;
/*<          t_per_samp = t_loop / n_samp >*/
		t_per_samp = t_loop / n_samp;
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
		do_fio(&c__1, (char *)&(loop), (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&time, (ftnlen)sizeof(double));
		do_fio(&c__1, (char *)&del_t, (ftnlen)sizeof(double));
		do_fio(&c__1, (char *)&(M), (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&k_max, (ftnlen)sizeof(double));
		do_fio(&c__1, (char *)&n_samp, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&t_per_samp, (ftnlen)sizeof(double));
		do_fio(&c__1, (char *)&n_coag, (ftnlen)sizeof(int));
		e_wsfe();
/*<       enddo >*/
	}
/*<       return >*/
	return 0;
/*<       end >*/
} /* mc_adapt__ */

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<    >*/
/* Subroutine */ int compute_n_samp_del_t(int M, double v_comp, 
	double *k_max, double p_max, double r_samp_max, 
	double del_t_max, int *n_samp, double *del_t)
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
    c__ = -(*k_max * 1. / v_comp / log(1 - p_max));
/*<       del_t = r_samp_max / c >*/
    *del_t = r_samp_max / c__;
/*<       if (del_t .gt. del_t_max) del_t = del_t_max >*/
    if (*del_t > del_t_max) {
	*del_t = del_t_max;
    }
/*<       r_samp = del_t * c >*/
    r_samp__ = *del_t * c__;
/*<       n_samp = int(r_samp * M*(M-1)/2) >*/
    *n_samp = (int) (r_samp__ * M * (M - 1) / 2);
/*<       return >*/
    return 0;
/*<       end >*/
} /* compute_n_samp_del_t */

