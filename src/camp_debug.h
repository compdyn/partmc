/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file with some debugging functions for use with camp_solver.c
 *
 */
#ifndef CAMP_DEBUG_H
#define CAMP_DEBUG_H

// file name prefix
int file_name_prefix = 1;

// Maximum size of output file names
#define MAX_FILE_NAME 256

// Number of points to advance state for output
#define N_OUTPUT_STATES 100

#ifdef PMC_DEBUG
#define PMC_DEBUG_SPEC_ 0
#define PMC_DEBUG_PRINT(x) \
  pmc_debug_print(sd->cvode_mem, x, false, 0, __LINE__, __func__)
#define PMC_DEBUG_PRINT_INT(x, y) \
  pmc_debug_print(sd->cvode_mem, x, false, y, __LINE__, __func__)
#define PMC_DEBUG_PRINT_FULL(x) \
  pmc_debug_print(sd->cvode_mem, x, true, 0, __LINE__, __func__)
#define PMC_DEBUG_JAC_STRUCT(J, x) pmc_debug_print_jac_struct((void *)sd, J, x)
#define PMC_DEBUG_JAC(J, x) pmc_debug_print_jac((void *)sd, J, x)
void pmc_debug_print(void *cvode_mem, const char *message, bool do_full,
                     const int int_val, const int line, const char *func) {
#ifdef PMC_USE_SUNDIALS
  CVodeMem cv_mem = (CVodeMem)cvode_mem;
  if (!(cv_mem->cv_debug_out)) return;
  printf(
      "\n[DEBUG] line %4d in %-20s(): %-25s %-4.0d t_n = %le h = %le q = %d "
      "hin = %le species %d(zn[0] = %le zn[1] = %le tempv = %le tempv1 = %le "
      "tempv2 = %le acor_init = %le last_yn = %le",
      line, func, message, int_val, cv_mem->cv_tn, cv_mem->cv_h, cv_mem->cv_q,
      cv_mem->cv_hin, PMC_DEBUG_SPEC_,
      NV_DATA_S(cv_mem->cv_zn[0])[PMC_DEBUG_SPEC_],
      NV_DATA_S(cv_mem->cv_zn[1])[PMC_DEBUG_SPEC_],
      NV_DATA_S(cv_mem->cv_tempv)[PMC_DEBUG_SPEC_],
      NV_DATA_S(cv_mem->cv_tempv1)[PMC_DEBUG_SPEC_],
      NV_DATA_S(cv_mem->cv_tempv2)[PMC_DEBUG_SPEC_],
      NV_DATA_S(cv_mem->cv_acor_init)[PMC_DEBUG_SPEC_],
      NV_DATA_S(cv_mem->cv_last_yn)[PMC_DEBUG_SPEC_]);
  if (do_full) {
    for (int i = 0; i < NV_LENGTH_S(cv_mem->cv_y); i++) {
      printf(
          "\n  zn[0][%3d] = % -le zn[1][%3d] = % -le tempv[%3d] = % -le "
          "tempv1[%3d] = % -le tempv2[%3d] = % -le acor_init[%3d] = % -le "
          "last_yn[%3d] = % -le",
          i, NV_DATA_S(cv_mem->cv_zn[0])[i], i, NV_DATA_S(cv_mem->cv_zn[1])[i],
          i, NV_DATA_S(cv_mem->cv_tempv)[i], i, NV_DATA_S(cv_mem->cv_tempv1)[i],
          i, NV_DATA_S(cv_mem->cv_tempv2)[i], i,
          NV_DATA_S(cv_mem->cv_acor_init)[i], i,
          NV_DATA_S(cv_mem->cv_last_yn)[i]);
    }
  }
#endif
}
void pmc_debug_print_jac_struct(void *solver_data, SUNMatrix J,
                                const char *message) {
#ifdef PMC_USE_SUNDIALS
  SolverData *sd = (SolverData *)solver_data;

  if (!(sd->debug_out)) return;
  int n_state_var = SM_COLUMNS_S(J);
  int i_elem = 0;
  int next_col = 0;
  printf("\n\n   Jacobian structure (↓ind →dep) - %s\n     ", message);
  for (int i_dep = 0; i_dep < n_state_var; i_dep++) printf("[%3d]", i_dep);
  for (int i_ind = 0; i_ind < n_state_var; i_ind++) {
    printf("\n[%3d]", i_ind);
    next_col = SM_INDEXPTRS_S(J)[i_ind + 1];
    for (int i_dep = 0; i_dep < n_state_var; i_dep++) {
      if (i_dep == SM_INDEXVALS_S(J)[i_elem] && i_elem < next_col) {
        printf(" %3d ", i_elem++);
      } else {
        printf("  -  ");
      }
    }
  }
#endif
}
void pmc_debug_print_jac(void *solver_data, SUNMatrix J, const char *message) {
#ifdef PMC_USE_SUNDIALS
  SolverData *sd = (SolverData *)solver_data;

  if (!(sd->debug_out)) return;
  int n_state_var = SM_COLUMNS_S(J);
  int i_elem = 0;
  int next_col = 0;
  printf("\n\n   Jacobian (↓ind →dep) - %s\n     ", message);
  for (int i_dep = 0; i_dep < n_state_var; i_dep++)
    printf("      [%3d]", i_dep);
  for (int i_ind = 0; i_ind < n_state_var; i_ind++) {
    printf("\n[%3d]   ", i_ind);
    next_col = SM_INDEXPTRS_S(J)[i_ind + 1];
    for (int i_dep = 0; i_dep < n_state_var; i_dep++) {
      if (i_dep == SM_INDEXVALS_S(J)[i_elem] && i_elem < next_col) {
        printf(" % -1.2le ", SM_DATA_S(J)[i_elem++]);
      } else {
        printf("     -     ");
      }
    }
  }
#endif
}

realtype pmc_jac_elem(SUNMatrix J, unsigned int j, unsigned int i) {
  for (int i_elem = SM_INDEXPTRS_S(J)[j]; i_elem < SM_INDEXPTRS_S(J)[j + 1];
       ++i_elem) {
    if (i == SM_INDEXVALS_S(J)[i_elem]) return SM_DATA_S(J)[i_elem];
  }
  return 0.0;
}
#else
#define PMC_DEBUG_PRINT(x)
#define PMC_DEBUG_PRINT_INT(x, y)
#define PMC_DEBUG_PRINT_FULL(x)
#define PMC_DEBUG_JAC_STRUCT(J, x)
#define PMC_DEBUG_JAC(J, x)
#endif

/** \brief Print some camp-chem data sizes
 *
 * \param md Pointer to the model data
 */
static void print_data_sizes(ModelData *md) {
  int *ptr = md->rxn_int_data;
  int n_rxn = ptr[0];

  printf("n_rxn: %d ", n_rxn);
  printf("n_state_var: %d", md->n_per_cell_state_var * md->n_cells);
  printf("n_dep_var: %d", md->n_per_cell_dep_var * md->n_cells);
}

/** \brief Print Jacobian matrix in format KLU SPARSE
 *
 * \param M Jacobian matrix
 */
static void print_jacobian(SUNMatrix M) {
  printf("\n NNZ JAC: %lld \n", SM_NNZ_S(M));
  printf("DATA | INDEXVALS:\n");
  for (int i = 0; i < SM_NNZ_S(M); i++) {
    printf("% -le \n", (SM_DATA_S(M))[i]);
    printf("%lld \n", (SM_INDEXVALS_S(M))[i]);
  }
  printf("PTRS:\n");
  for (int i = 0; i <= SM_NP_S(M); i++) {
    printf("%lld \n", (SM_INDEXPTRS_S(M))[i]);
  }
}

/** \brief Print derivative array
 *
 * \param deriv Derivative array
 */
static void print_derivative(N_Vector deriv) {
  // printf(" deriv length: %d\n", NV_LENGTH_S(deriv));
  for (int i = 0; i < NV_LENGTH_S(deriv); i++) {  // NV_LENGTH_S(deriv)
    printf(" deriv: % -le", NV_DATA_S(deriv)[i]);
    printf(" index: %d \n", i);
  }
}

/** \brief Evaluate the derivative and Jacobian near a given state
 *         for a specified species
 *
 * \param curr_time Current time
 * \param state State array
 * \param deriv Derivative array
 * \param solver_data Void pointer to solver data
 * \param f Pointer to derivative function
 * \param i_dep Dependent species index
 * \param i_ind Independent species index
 * \param d_rate_d_ind Change in rate for dependent species with change in
 *                         independent species
 * \param d_ind Increment to use in plot of rate_dep vs conc_ind
 */
void output_deriv_local_state(realtype curr_time, N_Vector state,
                              N_Vector deriv, void *solver_data,
                              int (*f)(realtype, N_Vector, N_Vector, void *),
                              int i_dep, int i_ind, double d_rate_d_ind,
                              double d_ind) {
  realtype *deriv_data = NV_DATA_S(deriv);
  realtype *state_data = NV_DATA_S(state);
  realtype rate_orig = deriv_data[i_dep];
  realtype ind_orig = state_data[i_ind];

  SolverData *sd = (SolverData *)solver_data;
  ModelData *md = &(sd->model_data);

  FILE *f_output;
  char file_name[MAX_FILE_NAME];
  sprintf(file_name, "local_%d_i%d_d%d", file_name_prefix++, i_ind, i_dep);
  f_output = fopen(file_name, "w");
  printf("\nOutputting deriv local state file: %s", file_name);

  // Output the loss of precision for all species
  ((SolverData *)solver_data)->output_precision = 1;
  if (f(curr_time, state, deriv, solver_data) != 0) {
    printf("\nERROR: Derivative failure\n\n");
  }
  ((SolverData *)solver_data)->output_precision = 0;

  // Output the current state
  fprintf(f_output, "#time %1.30le", curr_time);
  for (int i = 0; i < NV_LENGTH_S(state); ++i)
    fprintf(f_output, " [%3d] %1.30le", i, state_data[i]);
  fprintf(f_output, "\n");

  // Vary the model state and recalculate the derivative directly and using
  // the partial derivative provided
  for (int i = 0; i < N_OUTPUT_STATES; ++i) {
    state_data[i_ind] -= d_ind;
    if (state_data[i_ind] < 0.0) break;

    if (f(curr_time, state, deriv, solver_data) != 0) {
      printf("\nERROR: Derivative failure\n\n");
      break;
    }
    fprintf(f_output, "%d %1.30le %1.30le %1.30le\n", -i, state_data[i_ind],
            deriv_data[i_dep],
            (state_data[i_ind] - ind_orig) * d_rate_d_ind + rate_orig);
  }
  state_data[i_ind] = ind_orig;
  for (int i = 0; i < N_OUTPUT_STATES; ++i) {
    state_data[i_ind] += d_ind;
    if (f(curr_time, state, deriv, solver_data) != 0) {
      printf("\nERROR: Derivative failure\n\n");
      break;
    }
    fprintf(f_output, "%d %1.30le %1.30le %1.30le\n", i, state_data[i_ind],
            deriv_data[i_dep],
            (state_data[i_ind] - ind_orig) * d_rate_d_ind + rate_orig);
  }
  state_data[i_ind] = ind_orig;
  if (f(curr_time, state, deriv, solver_data) != 0) {
    printf("\nERROR: Derivative failure\n\n");
    EXIT_FAILURE;
  }

  printf("\nEnd output deriv local state file: %s", file_name);
  fclose(f_output);
}

#endif
