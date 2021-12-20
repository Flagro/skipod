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
  for (size_t i = 0; i < matrix_dim; ++i) {
    for (size_t j = 0; j < matrix_dim; ++j) {
      result[i][j] = 0;
    }
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

element_type **get_copy(element_type **matrix, size_t matrix_dim) {
  element_type **result = generate_square_matrix(matrix_dim);
  for (size_t i = 0; i < matrix_dim; ++i) {
    for (size_t j = 0; j < matrix_dim; ++j) {
      result[i][j] = matrix[i][j];
    }
  }
  return result;
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

void shift_block_matrix_row_left(element_type **matrix, element_type *buff,
                                 size_t matrix_dim, size_t block_size,
                                 size_t shift_row_id, size_t shift_block_cnt) {
  for (size_t i = shift_row_id * block_size;
       i < shift_row_id * block_size + block_size; ++i) {
    size_t shift_normal_size = shift_block_cnt * block_size;
    for (size_t j = 0; j < matrix_dim; ++j) {
      buff[j] = matrix[i][(j + shift_normal_size) % matrix_dim];
    }
    for (size_t j = 0; j < matrix_dim; ++j) {
      matrix[i][j] = buff[j];
    }
  }
}

void shift_block_matrix_column_up(element_type **matrix, element_type *buff,
                                  size_t matrix_dim, size_t block_size,
                                  size_t shift_column_id,
                                  size_t shift_block_cnt) {
  for (size_t i = shift_column_id * block_size;
       i < shift_column_id * block_size + block_size; ++i) {
    size_t shift_normal_size = shift_block_cnt * block_size;
    for (size_t j = 0; j < matrix_dim; ++j) {
      buff[j] = matrix[(j + shift_normal_size) % matrix_dim][i];
    }
    for (size_t j = 0; j < matrix_dim; ++j) {
      matrix[j][i] = buff[j];
    }
  }
}

void block_matrix_skewing_horizontal(element_type **matrix, element_type *buff,
                                     size_t matrix_dim, size_t block_size) {
  for (size_t i = 0; i < matrix_dim / block_size; ++i) {
    shift_block_matrix_row_left(matrix, buff, matrix_dim, block_size, i, i);
  }
}

void block_matrix_skewing_vertical(element_type **matrix, element_type *buff,
                                   size_t matrix_dim, size_t block_size) {
  for (size_t i = 0; i < matrix_dim / block_size; ++i) {
    shift_block_matrix_column_up(matrix, buff, matrix_dim, block_size, i, i);
  }
}

int verify_result(element_type **matrix1, element_type **matrix2,
                  element_type **result, size_t matrix_dim) {
  for (size_t i = 0; i < matrix_dim; ++i) {
    for (size_t j = 0; j < matrix_dim; ++j) {
      cur_result = 0;
      for (size_t k = 0; k < matrix_dim; ++k) {
        cur_result += matrix1[i][k] * matrix2[k][j];
      }
      if (cur_result != result[i][j]) {
        return 0;
      }
    }
  }
  return 1;
}

void serial_results_to_json(const char *calc_mode, int parallel_count,
                            size_t matrix_dim, elapsed_time_type elapsed_time) {
  return;
}

void serial_matrix_multiplication(element_type **matrix1,
                                  element_type **matrix2, element_type **result,
                                  size_t matrix_dim) {
  clock_t begin = clock();
  for (size_t i = 0; i < matrix_dim; ++i) {
    for (size_t j = 0; j < matrix_dim; ++j) {
      for (size_t k = 0; k < matrix_dim; ++k) {
        result[i][j] += matrix1[i][k] * matrix2[k][j];
      }
    }
  }
  clock_t end = clock();
  serial_results_to_json("serial", 1, matrix_dim, end - begin);
}
