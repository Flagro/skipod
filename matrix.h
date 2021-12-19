#pragma once

#include <stdlib.h>
#include <math.h>

#include "types.h"

element_type **generate_square_matrix(size_t matrix_dim) {
  element_type **result =
      (element_type **)calloc(matrix_dim, sizeof(element_type *));
  for (size_t i = 0; i < matrix_dim; ++i) {
    result[i] = (element_type *)calloc(matrix_dim, sizeof(element_type));
  }
  return result;
}

double get_random_double(double min, double max) {
  return min + (rand() / (RAND_MAX / (max - min)));
}

void init_random_square_matrix(element_type **matrix, size_t matrix_dim) {
  for (size_t i = 0; i < matrix_dim; ++i) {
    for (size_t j = 0; j < matrix_dim; ++j) {
      matrix[i][j] = get_random_double(-500.0, 500.0);
    }
  }
}

int matrix_equal(element_type **matrix1, element_type **matrix2,
               size_t matrix_dim) {
  static const double EPS = 1e-9;
  for (size_t i = 0; i < matrix_dim; ++i) {
    for (size_t j = 0; j < matrix_dim; ++j) {
      if (fabs(matrix1[i][j] - matrix2[i][j]) > EPS) {
          return 0;
      }
    }
  }
  return 1;
}

void matrix_free(element_type **matrix, size_t matrix_dim) {
  for (size_t i = 0; i < matrix_dim; ++i) {
    free(matrix[i]);
  }
  free(matrix);
}

time_type serial_matrix_multiplication(element_type **matrix1,
                                       element_type **matrix2,
                                       element_type **result,
                                       size_t matrix_dim) {}
