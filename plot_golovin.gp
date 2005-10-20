set logscale
plot [1e-7:1e-3] [1e2:1e10] "out_golovin_exact_num_avg.d" using 2:3 w l
replot "out_golovin_adapt_num_avg.d" using 2:3
replot "out_golovin_fix_num_avg.d" using 2:3
replot "out_golovin_var_num_avg.d" using 2:3
replot "out_golovin_exact_num_avg.d" using 2:8 w l
replot "out_golovin_adapt_num_avg.d" using 2:8
replot "out_golovin_fix_num_avg.d" using 2:8
replot "out_golovin_var_num_avg.d" using 2:8
replot "out_golovin_exact_num_avg.d" using 2:13 w l
replot "out_golovin_adapt_num_avg.d" using 2:13
replot "out_golovin_fix_num_avg.d" using 2:13
replot "out_golovin_var_num_avg.d" using 2:13
