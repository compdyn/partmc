// Copyright (C) 2009-2010 Matthew West
// Licensed under the GNU General Public License version 2 or (at your
// option) any later version. See the file COPYING for details.

// \file
// Helper file to call the SUNDIALS solver for condensation.

#include <stdlib.h>
#include <stdio.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <cvode/cvode_impl.h>

static int condense_vf(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int condense_check_flag(void *flagvalue, char *funcname, int opt);

/*******************************************************/
// solver block
static int condense_solver_Init(CVodeMem cv_mem);

static int condense_solver_Setup(CVodeMem cv_mem, int convfail, N_Vector ypred,
				 N_Vector fpred, booleantype *jcurPtr, 
				 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

static int condense_solver_Solve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
				 N_Vector ycur, N_Vector fcur);

static void condense_solver_Free(CVodeMem cv_mem);
/*******************************************************/

int condense_solver(int neq, double *x_f, double *abstol_f, double reltol_f,
		    double t_initial_f, double t_final_f)
{
	realtype reltol, t_initial, t_final, t, tout;
	N_Vector y, abstol;
	void *cvode_mem;
	CVodeMem cv_mem;
	int flag, i, pretype, maxl;
	realtype *y_data, *abstol_data;

	y = abstol = NULL;
	cvode_mem = NULL;

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
	
	reltol = reltol_f;
	t_initial = t_initial_f;
	t_final = t_final_f;

	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (condense_check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

	flag = CVodeInit(cvode_mem, condense_vf, t_initial, y);
	if (condense_check_flag(&flag, "CVodeInit", 1)) return(1);

	flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
	if (condense_check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

	flag = CVodeSetMaxNumSteps(cvode_mem, 100000);
	if (condense_check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return(1);

        /*******************************************************/
	// dense solver
	//flag = CVDense(cvode_mem, neq);
	//if (condense_check_flag(&flag, "CVDense", 1)) return(1);
        /*******************************************************/

        /*******************************************************/
	// iterative solver
	//pretype = PREC_LEFT;
	//maxl = 0;
	//flag = CVSptfqmr(cvode_mem, pretype, maxl); 
	//if (condense_check_flag(&flag, "CVSptfqmr", 1)) return(1);

	//flag = CVSpilsSetJacTimesVecFn(cvode_mem, condense_jtimes); 
	//if (condense_check_flag(&flag, "CVSpilsSetJacTimesVecFn", 1)) return(1);

	//flag = CVSpilsSetPreconditioner(cvode_mem, NULL, condense_prec);
	//if (condense_check_flag(&flag, "CVSpilsSetPreconditioner", 1)) return(1);
        /*******************************************************/

        /*******************************************************/
	// explicit solver
	cv_mem = (CVodeMem)cvode_mem;
	cv_mem->cv_linit = condense_solver_Init;
	cv_mem->cv_lsetup = condense_solver_Setup;
	cv_mem->cv_lsolve = condense_solver_Solve;
	cv_mem->cv_lfree = condense_solver_Free;
        /*******************************************************/

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

static int condense_solver_Init(CVodeMem cv_mem)
{
	return(0);
}

static int condense_solver_Setup(CVodeMem cv_mem, int convfail, N_Vector ypred,
				 N_Vector fpred, booleantype *jcurPtr, 
				 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
	return(0);
}

static int condense_solver_Solve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
				 N_Vector ycur, N_Vector fcur)
{
	realtype *b_data, *ycur_data, *fcur_data;
	int i, neq;
	double *b_f, *ycur_f, *fcur_f;
	double t, gamma;

	neq = NV_LENGTH_S(b);
	b_data = NV_DATA_S(b);
	ycur_data = NV_DATA_S(ycur);
	fcur_data = NV_DATA_S(fcur);

	b_f = malloc(neq * sizeof(double));
	ycur_f = malloc(neq * sizeof(double));
	fcur_f = malloc(neq * sizeof(double));

	t = cv_mem->cv_tn;
	gamma = cv_mem->cv_gamma;

	for (i = 0; i < neq; i++) {
		b_f[i] = b_data[i];
		ycur_f[i] = ycur_data[i];
		fcur_f[i] = fcur_data[i];
	}
	condense_jac_solve_f(neq, t, ycur_f, fcur_f, b_f, gamma);
	for (i = 0; i < neq; i++) {
		b_data[i] = b_f[i];
	}
	
	free(b_f);
	free(ycur_f);
	free(fcur_f);
	return(0);
}

static void condense_solver_Free(CVodeMem cv_mem)
{
}
