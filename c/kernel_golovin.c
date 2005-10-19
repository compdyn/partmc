
#include <math.h>

#include "kernel_golovin.h"

/***************************************************************************/

double kernel_golovin(double v1, double v2)
{
	double beta_1 = 1000;

	return beta_1 * (v1 + v2);
}

/***************************************************************************/

void soln_golovin_exp(int n_bin, double *bin_v,
		      double *bin_r, double *bin_g, int *bin_n, 
		      double dlnr, double time, double n_0,
		      double v_0, double rho_p, double v_comp)
{
	int k;
	double b, t, x, nn, tau, rat_v, beta_1;

	beta_1 = kernel_golovin(1, 0);
	if (time == 0) {
		for (k = 0; k < n_bin; k++) {
			bin_n[k] = 4 * M_PI * pow(bin_r[k],3) * n_0 / v_0
				* exp(-(bin_v[k] / v_0));
		}
	} else {
		tau = n_0 * v_0 * beta_1 * time;
		t = 1 - exp(-tau);
		for (k = 0; k < n_bin; k++) {
			rat_v = bin_v[k] / v_0;
			x = 2 * rat_v * sqrt(t);
			if (x < 500) {
				b = bessi1(x);
			} else {
				b = 0;
			}
			nn = n_0 / bin_v[k] * (1 - t) / sqrt(t)
				* exp(-((t + 1) * rat_v)) * b;
			bin_n[k] = 4 * M_PI * pow(bin_r[k],3) * nn;
		}
	}
	for (k = 0; k < n_bin; k++) {
		bin_g[k] = 4/3 * M_PI * rho_p * pow(bin_r[k],3) * bin_n[k];
	}
	for (k = 0; k < n_bin; k++) {
		bin_g[k] = bin_g[k] * dlnr * v_comp;
		bin_n[k] = bin_n[k] * dlnr * v_comp;
	}
}

/***************************************************************************/

double bessi1(double x)
{
    double p1 = .5;
    double p2 = .87890594;
    double p3 = .51498869;
    double p4 = .15084934;
    double p5 = .02658733;
    double p6 = .00301532;
    double p7 = 3.2411e-4;

    double q1 = .39894228;
    double q2 = -.03988024;
    double q3 = -.00362018;
    double q4 = .00163801;
    double q5 = -.01031555;
    double q6 = .02282967;
    double q7 = -.02895312;
    double q8 = .01787654;
    double q9 = -.00420059;

    double r, y, ax;

    if (fabs(x) < 3.75) {
	    y = pow(x / 3.75, 2);
	    r = ((((((p7 * y + p6) * y + p5) * y + p4) * y
		   + p3) * y + p2) * y + p1) * x;
    } else {
	ax = fabs(x);
	y = 3.75 / ax;
	r = ((((((((q9 * y + q8) * y + q7) * y + q6) * y + q5) * y
		+ q4) * y + q3) * y + q2) * y + q1) * exp(ax) / sqrt(ax);
	if (x < 0)
	    r = -r;
    }
    return r;
}

/***************************************************************************/
