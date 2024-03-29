#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "assert.h"
#include "cannon_matrix_mpi.h"
#include "cannon_matrix_omp.h"
#include "matrix.h"
#include "types.h"

char *get_arg_value(int argc, char *argv[], const char *arg_name) {
  for (int i = 1; i < argc; ++i) {
    if (strstr(argv[i], arg_name)) {
      char *value = strchr(argv[i], '=');
      ++value;
      char *ans = (char *)calloc(strlen(value) + 1, sizeof(char));
      strcpy(ans, value);
      return ans;
    }
  }
  return NULL;
}

int main(int argc, char *argv[]) {
  char *calc_mode = get_arg_value(argc, argv, "-mode");
  char *init_matrix_dim_str = get_arg_value(argc, argv, "-n");
  size_t init_matrix_dim;
  if (init_matrix_dim_str) {
    init_matrix_dim = strtol(init_matrix_dim_str, NULL, 10);
  }
  char *parallel_count_str = get_arg_value(argc, argv, "-parallel");
  int parallel_count;
  if (parallel_count_str) {
    parallel_count = strtol(parallel_count_str, NULL, 10);
  }
  size_t block_matrix_dim = (int)sqrt(parallel_count);
  size_t matrix_dim =
      ((init_matrix_dim + block_matrix_dim - 1) / block_matrix_dim) *
      block_matrix_dim;

  element_type **matrix1 = generate_square_matrix(matrix_dim);
  init_random_square_matrix(matrix1, init_matrix_dim);
  element_type **matrix1_copy = get_copy(matrix1, matrix_dim);

  element_type **matrix2 = generate_square_matrix(matrix_dim);
  init_random_square_matrix(matrix2, init_matrix_dim);
  element_type **matrix2_copy = get_copy(matrix2, matrix_dim);

  element_type **result = generate_square_matrix(matrix_dim);
  elapsed_time_type elapsed_time;

  if (!calc_mode) {
    printf("No calculation mode were specified, exiting the program\n");
    return 0;
  } else if (!strcmp(calc_mode, "serial")) {
    serial_matrix_multiplication(matrix1, matrix2, result, matrix_dim);
  } else if (!strcmp(calc_mode, "omp")) {
    cannon_matrix_multiplication_omp(matrix1, matrix2, result, matrix_dim,
                                     block_matrix_dim);
  } else if (!strcmp(calc_mode, "mpi")) {
    cannon_matrix_multiplication_mpi(matrix1, matrix2, result, matrix_dim,
                                     block_matrix_dim);
  } else {
    printf("Invalid passed mode argument \"%s\", exiting the program\n",
           calc_mode);
    return 1;
  }

  matrix_free(matrix1, matrix_dim);
  matrix_free(matrix2, matrix_dim);
  matrix_free(result, matrix_dim);

  if (parallel_count_str) {
    free(parallel_count_str);
  }
  if (init_matrix_dim_str) {
    free(init_matrix_dim_str);
  }
  if (calc_mode) {
    free(calc_mode);
  }
  return 0;
}
