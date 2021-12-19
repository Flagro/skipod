#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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
      matrix[i][j] = rand() % 5; // get_random_double(-500.0, 500.0);
    }
  }
}

void print_matrix(element_type **matrix, size_t matrix_dim) {
  for (size_t i = 0; i < matrix_dim; ++i) {
    for (size_t j = 0; j < matrix_dim; ++j) {
      printf("%lf ", matrix[i][j]);
    }
    printf("\n");
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

void shift_single_row_left(element_type **matrix, size_t matrix_dim,
                           size_t row_id) {
  for (size_t i = 0; i + 1 < matrix_dim; ++i) {
    element_type tmp = matrix[row_id][i];
    matrix[row_id][i] = matrix[row_id][i + 1];
    matrix[row_id][i + 1] = tmp;
  }
}

void shift_single_column_up(element_type **matrix, size_t matrix_dim,
                            size_t column_id) {
  for (size_t i = 0; i + 1 < matrix_dim; ++i) {
    element_type tmp = matrix[i][column_id];
    matrix[i][column_id] = matrix[i + 1][column_id];
    matrix[i + 1][column_id] = tmp;
  }
}

void shift_rows_left(element_type **matrix, size_t matrix_dim) {
  for (size_t i = 0; i < matrix_dim; ++i) {
    shift_single_row_left(matrix, matrix_dim, i);
  }
}

void shift_columns_up(element_type **matrix, size_t matrix_dim) {
  for (size_t i = 0; i < matrix_dim; ++i) {
    shift_single_column_up(matrix, matrix_dim, i);
  }
}

void matrix_skewing_horizontal(element_type **matrix, size_t matrix_dim) {
  for (size_t i = 0; i < matrix_dim; ++i) {
    for (size_t j = 0; j < i; ++j) {
      shift_single_row_left(matrix, matrix_dim, i);
    }
  }
}

void matrix_skewing_vertical(element_type **matrix, size_t matrix_dim) {
  for (size_t i = 0; i < matrix_dim; ++i) {
    for (size_t j = 0; j < i; ++j) {
      shift_single_column_up(matrix, matrix_dim, i);
    }
  }
}

void multiply_matrix_blocks(element_type **matrix1, element_type **matrix2,
                            element_type **result, size_t block_dim,
                            size_t i_start, size_t j_start) {
  for (size_t i = 0; i < block_dim; ++i) {
    for (size_t k = 0; k < block_dim; ++k) {
      double r = matrix1[i_start + i][j_start + k];
      for (size_t j = 0; j < block_dim; ++j) {
        result[i_start + i][j_start + j] +=
            r * matrix2[i_start + k][j_start + j];
      }
    }
  }
}

time_type serial_matrix_multiplication(element_type **matrix1,
                                       element_type **matrix2,
                                       element_type **result,
                                       size_t matrix_dim) {
  clock_t begin = clock();
  for (size_t i = 0; i < matrix_dim; ++i) {
    for (size_t k = 0; k < matrix_dim; ++k) {
      double r = matrix1[i][k];
      for (size_t j = 0; j < matrix_dim; ++j) {
        result[i][j] += r * matrix2[k][j];
      }
    }
  }
  clock_t end = clock();
  return (float)(end - begin) / CLOCKS_PER_SEC;
  ;
}
