
#include <stdlib.h>
#include <stdio.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>

static int vf(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int jac(int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

void vf_f(int neq, double t_f, double *y_f, double *ydot_f);

static void PrintFinalStats(void *cvode_mem);

static int check_flag(void *flagvalue, char *funcname, int opt);

int do_run(int neq, double *x_f, double *abstol_f, double reltol_f,
	   double t_initial_f, double t_final_f)
{
	realtype reltol, t_initial, t_final, t, tout;
	N_Vector y, abstol;
	void *cvode_mem;
	int flag, i;
	realtype *y_data, *abstol_data;
	
	y = abstol = NULL;
	cvode_mem = NULL;

	printf("neq = %d\n", neq);
	y = N_VNew_Serial(neq);
	if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
	abstol = N_VNew_Serial(neq); 
	if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return(1);
	
	y_data = NV_DATA_S(y);
	abstol_data = NV_DATA_S(abstol);
	for (i = 0; i < neq; i++) {
		y_data[i] = x_f[i];
		abstol_data[i] = abstol_f[i];
		printf("i %d, x %f, abstol %f\n", i, x_f[i], abstol_f[i]);
	}
	
	reltol = reltol_f;
	t_initial = t_initial_f;
	t_final = t_final_f;
	printf("t_initial = %f, t_final = %f\n", t_initial, t_final);

	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
	flag = CVodeInit(cvode_mem, vf, t_initial, y);
	if (check_flag(&flag, "CVodeInit", 1)) return(1);
	flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);
	flag = CVDense(cvode_mem, neq);
	if (check_flag(&flag, "CVDense", 1)) return(1);
	flag = CVDlsSetDenseJacFn(cvode_mem, jac);
	if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);

	t = t_initial;
	flag = CVode(cvode_mem, t_final, y, &t, CV_NORMAL);
	if (check_flag(&flag, "CVode", 1)) return(1);

	for (i = 0; i < neq; i++) {
		x_f[i] = y_data[i];
	}

	PrintFinalStats(cvode_mem);
	N_VDestroy_Serial(y);
	N_VDestroy_Serial(abstol);
	CVodeFree(&cvode_mem);
	return(0);
}


/*
 * f routine. Compute function f(t,y) = - y.
 */

static int vf(realtype t, N_Vector y, N_Vector ydot, void *user_data)
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
	vf_f(neq, t, y_f, ydot_f);
	for (i = 0; i < neq; i++) {
		ydot_data[i] = ydot_f[i];
	}
	
	free(y_f);
	free(ydot_f);
	return(0);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy = -1. *
 */

static int jac(int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	realtype *y_data;
	int i, j, neq;
	double *y_f, *df_f;

	neq = NV_LENGTH_S(y);
	y_data = NV_DATA_S(y);

	y_f = malloc(neq * sizeof(double));
	df_f = malloc(neq * neq * sizeof(double));

	for (i = 0; i < neq; i++) {
		y_f[i] = y_data[i];
	}
	jac_f(neq, t, y_f, df_f);
	for (i = 0; i < neq; i++) {
		for (j = 0; j < neq; j++) {
			DENSE_ELEM(J,i,j) = df_f[i + j * neq];
		}
	}
	
	free(y_f);
	free(df_f);
	return(0);
}

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
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

static int check_flag(void *flagvalue, char *funcname, int opt)
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
