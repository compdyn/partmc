/* Copyright (C) 2019 Matthew Dawson
 * Licensed under the GNU General Public License version 1 or (at your
 * option) any later version. See the file COPYING for details.
 */
/** \file
 * \brief Tests for the Jacobian struct and related functions
 */
#include <stdio.h>
#include <stdlib.h>
#include "../test_common.h"
#include "../../src/Jacobian.h"

// Number of variables
#define NUM_VAR 5

// Number of Jacobian elements
#define NUM_ELEM 4

// Value for tests
#define REF_VAL (-9999.0)

int main(int argc, char * argv[]) {

  int errors = 0;

  Jacobian jac;
  errors+=ASSERT_MSG(jacobian_initialize_empty(&jac, NUM_VAR)==1, "281479700");

  // register elements
  jacobian_register_element(&jac, 3, 2);
  jacobian_register_element(&jac, 0, 0);
  jacobian_register_element(&jac, 2, 4);
  jacobian_register_element(&jac, 1, 2);

  // build the sparse matrix
  errors+=ASSERT_MSG(jacobian_build_matrix(&jac)==1, "403536047");

  errors+=ASSERT_MSG(jacobian_number_of_elements(jac)==4, "501315968");

  errors+=ASSERT_MSG(jacobian_column_pointer_value(jac,0)==0, "157643387");
  errors+=ASSERT_MSG(jacobian_column_pointer_value(jac,1)==1, "316270733");
  errors+=ASSERT_MSG(jacobian_column_pointer_value(jac,2)==1, "146113829");
  errors+=ASSERT_MSG(jacobian_column_pointer_value(jac,3)==3, "593481675");
  errors+=ASSERT_MSG(jacobian_column_pointer_value(jac,4)==3, "205857922");
  errors+=ASSERT_MSG(jacobian_column_pointer_value(jac,5)==4, "935701017");

  errors+=ASSERT_MSG(jacobian_row_index(jac,0)==0, "824835911");
  errors+=ASSERT_MSG(jacobian_row_index(jac,1)==1, "483068864");
  errors+=ASSERT_MSG(jacobian_row_index(jac,2)==3, "647961461");
  errors+=ASSERT_MSG(jacobian_row_index(jac,3)==2, "260337708");

  errors+=ASSERT_MSG(jacobian_get_element_id(jac, 3, 2)==2, "271491762");
  errors+=ASSERT_MSG(jacobian_get_element_id(jac, 0, 0)==0, "718859608");
  errors+=ASSERT_MSG(jacobian_get_element_id(jac, 2, 4)==3, "266227455");
  errors+=ASSERT_MSG(jacobian_get_element_id(jac, 1, 2)==1, "431120052");

  jacobian_reset(jac);
  double out_vals[NUM_ELEM+1];
  for(int i=0; i<NUM_ELEM+1; ++i)out_vals[i] = REF_VAL;
  jacobian_output(jac, out_vals);

  errors+=ASSERT_MSG(out_vals[0]==0.0, "181518712");
  errors+=ASSERT_MSG(out_vals[1]==0.0, "518473747");
  errors+=ASSERT_MSG(out_vals[2]==0.0, "965841593");
  errors+=ASSERT_MSG(out_vals[3]==0.0, "295742591");
  errors+=ASSERT_MSG(out_vals[4]==REF_VAL, "285213977");

  jacobian_add_value(jac, 0, 1, 70.0);
  jacobian_add_value(jac, 3, 1, 20.0);
  jacobian_add_value(jac, 2, 0, 40.0);
  jacobian_add_value(jac, 3, 0, 80.0);
  jacobian_add_value(jac, 3, 1, 10.0);
  jacobian_add_value(jac, 1, 0, 60.0);

  for(int i=0; i<NUM_ELEM+1; ++i)out_vals[i] = REF_VAL;
  jacobian_output(jac, out_vals);

  errors+=ASSERT_CLOSE_MSG(out_vals[0], -70.0, "819935518");
  errors+=ASSERT_CLOSE_MSG(out_vals[1],  60.0, "649778614");
  errors+=ASSERT_CLOSE_MSG(out_vals[2],  40.0, "197146461");
  errors+=ASSERT_CLOSE_MSG(out_vals[3],  50.0, "991997956");
  errors+=ASSERT_MSG(out_vals[4]==REF_VAL, "256890554");

  jacobian_free(&jac);

  // check Jacobian with a column with more than the buffer size of rows
  errors+=ASSERT_MSG(jacobian_initialize_empty(&jac, 20)==1, "895645677");
  for (int i=0; i<20; ++i) jacobian_register_element(&jac, i, 3);
  errors+=ASSERT_MSG(jacobian_build_matrix(&jac)==1, "354320568");
  errors+=ASSERT_MSG(jacobian_number_of_elements(jac)==20, "794970868");
  for (int i=0; i<=3; ++i) {
    errors+=ASSERT_MSG(jacobian_column_pointer_value(jac, i)==0, "160199958");
  }
  for (int i=4; i<20; ++i) {
    errors+=ASSERT_MSG(jacobian_column_pointer_value(jac, i)==20, "388647716");
  }
  for (int i=0; i<20; ++i) {
    errors+=ASSERT_MSG(jacobian_row_index(jac, i)==i, "374760329");
  }
  for (int i=0; i<20; ++i) {
    errors+=ASSERT_MSG(jacobian_get_element_id(jac, i, 3)==i, "317374065");
  }
  for (int i=0; i<20; ++i) {
    jacobian_add_value(jac, i, 0, 10.0+1.0*i);
    jacobian_add_value(jac, i, 1, 10.0+2.0*i);
  }
  double out_vals2[21];
  for(int i=0; i<21; ++i)out_vals2[i] = REF_VAL;
  jacobian_output(jac, out_vals2);
  for (int i=0; i<20; ++i) {
    errors+=ASSERT_CLOSE_MSG(out_vals2[i], -1.0*i, "125158989");
  }
  errors+=ASSERT_MSG(out_vals2[20]==REF_VAL, "174374468");

  jacobian_free(&jac);

  if (errors==0) {
    printf("\nPASS\n");
  } else {
    printf("\nFAIL\n");
  }

}
