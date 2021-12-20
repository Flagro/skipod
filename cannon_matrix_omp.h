#pragma once

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "omp.h"

#include "matrix.h"
#include "types.h"

void omp_results_to_json(const char *calc_mode, int parallel_count,
                     size_t matrix_dim, elapsed_time_type elapsed_time) {
  return;
}

void add_to_result(element_type **matrix1, element_type **matrix2,
                   element_type **result_accumulator, size_t matrix_dim,
                   size_t block_size, int block_matrix_dim) {
  omp_set_num_threads(block_matrix_dim * block_matrix_dim);
#pragma omp parallel shared(matrix1, matrix2, result_accumulator)
  {
    int thread_id = omp_get_thread_num();
    size_t cur_block_matrix_row = thread_id / block_matrix_dim;
    size_t cur_block_matrix_col = thread_id % block_matrix_dim;
    size_t i_start = cur_block_matrix_row * block_size;
    size_t j_start = cur_block_matrix_col * block_size;
    for (size_t i = 0; i < block_size; ++i) {
      for (size_t j = 0; j < block_size; ++j) {
        for (size_t k = 0; k < block_size; ++k) {
          result_accumulator[i_start + i][j_start + j] +=
              matrix1[i_start + i][j_start + k] *
              matrix2[i_start + k][j_start + j];
        }
      }
    }
  }
}

void cannon_matrix_multiplication_omp(element_type **matrix1,
                                                   element_type **matrix2,
                                                   element_type **result,
                                                   size_t matrix_dim,
                                                   int block_matrix_dim) {
  size_t block_size = matrix_dim / block_matrix_dim;
  element_type *buff = (element_type *)calloc(matrix_dim, sizeof(element_type));

  elapsed_time_type begin = omp_get_wtime();
  block_matrix_skewing_horizontal(matrix1, buff, matrix_dim, block_size);
  block_matrix_skewing_vertical(matrix2, buff, matrix_dim, block_size);
  for (int i = 0; i < block_matrix_dim; ++i) {
    add_to_result(matrix1, matrix2, result, matrix_dim, block_size,
                  block_matrix_dim);
    for (int j = 0; j < block_matrix_dim; ++j) {
      shift_block_matrix_row_left(matrix1, buff, matrix_dim, block_size, j, 1);
      shift_block_matrix_column_up(matrix2, buff, matrix_dim, block_size, j, 1);
    }
  }
  elapsed_time_type end = omp_get_wtime();
  free(buff);
  assert(verify_result(matrix1, matrix2, result, matrix_dim) &&
         "the omp answer is incorrect");
  omp_results_to_json("omp", block_matrix_dim * block_matrix_dim, matrix_dim, end - begin);
}
