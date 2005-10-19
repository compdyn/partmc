#ifndef KERNEL_GOLOVIN_H
#define KERNEL_GOLOVIN_H

/*!
  \file kernel_golovin.h
  \brief Golovin coagulation kernel and exact solutions.
*/

//! Golovin coagulation kernel.
/*!
\param v1       [Input] volume of first particle
\param v2       [Input] volume of second particle
\return         the coagulation kernel between \c v1 and \c v2
*/
double kernel_golovin(double v1, double v2);

//! Exact solution to the Golovin kernel with an exponential initial condition.
/*!
\param n_bin    [Input] number of bins
\param bin_v    [Input] volume of particles in bins
\param bin_r    [Input] radius of particles in bins
\param bin_g    [Output] mass in bins
\param bin_n    [Output] number in bins
\param dlnr     [Input] bin scale factor
\param time     [Input] cubin_rent time
\param n_0      [Input] particle number concentration (#/m^3)
\param v_0      [Input]
\param rho_p    [Input] particle density (kg/m^3)
\param v_comp   [Input] computational volume
*/
void soln_golovin_exp(int n_bin, double *bin_v,
		      double *bin_r, double *bin_g, int *bin_n, 
		      double dlnr, double time, double n_0,
		      double v_0, double rho_p, double v_comp);

//! Bessel function.
/*
\param x        [Input] function argument
\output         function value
*/
double bessi1(double x);

#endif
