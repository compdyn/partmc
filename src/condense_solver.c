// Copyright (C) 2009 Matthew West
// Licensed under the GNU General Public License version 2 or (at your
// option) any later version. See the file COPYING for details.

// \file
// Helper file to call the SUNDIALS solver for condensation.

#include <stdlib.h>
#include <stdio.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_sptfqmr.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>

static int condense_vf(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int condense_jtimes(N_Vector v, N_Vector Jv, realtype t, N_Vector y,
			   N_Vector fy, void *user_data, N_Vector tmp);

//static void (*condense_vf_f)(int neq, double t, double *y_f, double *ydot_f);

//static void (*condense_jtimes_f)(int neq, double t, double *y_f, double *v_f, double *Jv_f);

static int condense_check_flag(void *flagvalue, char *funcname, int opt);

int condense_solver(int neq, double *x_f, double *abstol_f, double reltol_f,
		    double t_initial_f, double t_final_f)
//		    void (*condense_vf_f_p)(int neq, double t, double *y_f, double *ydot_f),
//		    void (*condense_jtimes_f_p)(int neq, double t, double *y_f, double *v_f, double *Jv_f))
{
	realtype reltol, t_initial, t_final, t, tout;
	N_Vector y, abstol;
	void *cvode_mem;
	int flag, i, pretype, maxl;
	realtype *y_data, *abstol_data;

	y = abstol = NULL;
	cvode_mem = NULL;

	//condense_vf_f = condense_vf_f_p;
	//condense_jtimes_f = condense_jtimes_f_p;

	printf("neq = %d\n", neq);

	y = N_VNew_Serial(neq);
	if (condense_check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

	abstol = N_VNew_Serial(neq); 
	if (condense_check_flag((void *)abstol, "N_VNew_Serial", 0)) return(1);
	
	y_data = NV_DATA_S(y);
	abstol_data = NV_DATA_S(abstol);
	for (i = 0; i < neq; i++) {
		y_data[i] = x_f[i];
		abstol_data[i] = abstol_f[i];
	}
	i = 0;
	printf("condense_solver: i %d, x %g, abstol %g\n", i, x_f[i], abstol_f[i]);
	i = neq - 1;
	printf("condense_solver: i %d, x %g, abstol %g\n", i, x_f[i], abstol_f[i]);
	
	reltol = reltol_f;
	t_initial = t_initial_f;
	t_final = t_final_f;
	printf("condense_solver: t_initial = %f, t_final = %f\n", t_initial, t_final);

	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (condense_check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

	flag = CVodeInit(cvode_mem, condense_vf, t_initial, y);
	if (condense_check_flag(&flag, "CVodeInit", 1)) return(1);

	flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
	if (condense_check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

	flag = CVodeSetMaxNumSteps(cvode_mem, 100000);
	if (condense_check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return(1);

	pretype = PREC_NONE;
	maxl = 0;
	flag = CVSptfqmr(cvode_mem, pretype, maxl); 
	if (condense_check_flag(&flag, "CVSptfqmr", 1)) return(1);

	flag = CVSpilsSetJacTimesVecFn(cvode_mem, condense_jtimes); 
	if (condense_check_flag(&flag, "CVSpilsSetJacTimesVecFn", 1)) return(1);

	t = t_initial;
	flag = CVode(cvode_mem, t_final, y, &t, CV_NORMAL);
	if (condense_check_flag(&flag, "CVode", 1)) return(1);

	for (i = 0; i < neq; i++) {
		x_f[i] = y_data[i];
	}

	N_VDestroy_Serial(y);
	N_VDestroy_Serial(abstol);
	CVodeFree(&cvode_mem);
	return(0);
}

static int condense_vf(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	realtype *y_data, *ydot_data;
	int i, neq;
	double *y_f, *ydot_f;

	neq = NV_LENGTH_S(y);
	y_data = NV_DATA_S(y);
	ydot_data = NV_DATA_S(ydot);

	y_f = malloc(neq * sizeof(double));
	ydot_f = malloc(neq * sizeof(double));

	for (i = 0; i < neq; i++) {
		y_f[i] = y_data[i];
	}
	condense_vf_f(neq, t, y_f, ydot_f);
	for (i = 0; i < neq; i++) {
		ydot_data[i] = ydot_f[i];
	}
	
	free(y_f);
	free(ydot_f);
	return(0);
}

static int condense_jtimes(N_Vector v, N_Vector Jv, realtype t, N_Vector y,
			   N_Vector fy, void *user_data, N_Vector tmp)
{
	realtype *y_data, *v_data, *Jv_data;
	int i, neq;
	double *y_f, *v_f, *Jv_f;

	neq = NV_LENGTH_S(y);
	y_data = NV_DATA_S(y);
	v_data = NV_DATA_S(v);
	Jv_data = NV_DATA_S(Jv);

	y_f = malloc(neq * sizeof(double));
	v_f = malloc(neq * sizeof(double));
	Jv_f = malloc(neq * sizeof(double));

	for (i = 0; i < neq; i++) {
		y_f[i] = y_data[i];
		v_f[i] = v_data[i];
	}
	condense_jtimes_f(neq, t, y_f, v_f, Jv_f);
	for (i = 0; i < neq; i++) {
		Jv_data[i] = Jv_f[i];
	}
	
	free(y_f);
	free(v_f);
	free(Jv_f);
	return(0);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int condense_check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
