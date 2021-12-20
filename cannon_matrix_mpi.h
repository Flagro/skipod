#pragma once

#include <assert.h>
#include <mpi.h>

#include "matrix.h"
#include "types.h"

void mpi_results_to_json(const char *calc_mode, int parallel_count,
                         size_t matrix_dim, elapsed_time_type elapsed_time) {
  return;
}

elapsed_time_type cannon_matrix_multiplication_mpi(element_type **matrix1,
                                                   element_type **matrix2,
                                                   element_type **result,
                                                   size_t matrix_dim,
                                                   int block_matrix_dim) {
  size_t block_size = matrix_dim / block_matrix_dim;

  int initialized;
  MPI_Initialized(&initialized);
  if (!initialized) {
    MPI_Init(NULL, NULL);
  }

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  element_type *vectored_matrix1 =
      calloc(matrix_dim * matrix_dim, sizeof(element_type));
  element_type *vectored_matrix2 =
      calloc(matrix_dim * matrix_dim, sizeof(element_type));
  if (!world_rank) {
    // making only 0th process generated matrix valid
    for (size_t i = 0; i < matrix_dim; ++i) {
      for (size_t j = 0; j < matrix_dim; ++j) {
        vectored_matrix1[i * matrix_dim + j] = matrix1[i][j];
        vectored_matrix2[i * matrix_dim + j] = matrix2[i][j];
      }
    }
    MPI_Bcast(&vectored_matrix1, matrix_dim * matrix_dim, MPI_DOUBLE, 0,
              MPI_COMM_WORLD);
    MPI_Bcast(&vectored_matrix2, matrix_dim * matrix_dim, MPI_DOUBLE, 0,
              MPI_COMM_WORLD);

    int dims[2];
    dims[0] = dims[1] = block_matrix_dim;
    int periods[2] = {1, 1};
    MPI_Comm comm_2d_cart;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm_2d_cart);

    int comm_2d_cart_rank;
    MPI_Comm_rank(comm_2d_cart, &comm_2d_cart_rank);

    int comm_2d_cart_coords[2];
    MPI_Cart_coords(comm_2d_cart, comm_2d_cart_rank, 2, comm_2d_cart_coords);

    size_t block_row = comm_2d_cart_coords[0];
    size_t block_col = comm_2d_cart_coords[1];

    element_type *vectored_matrix1_block =
        calloc(block_size * block_size, sizeof(element_type));
    element_type *vectored_matrix2_block =
        calloc(block_size * block_size, sizeof(element_type));
    element_type *vectored_result_block =
        calloc(block_size * block_size, sizeof(element_type));

    for (size_t i = 0; i < block_size; ++i) {
      for (size_t j = 0; j < block_size; ++j) {
        vectored_matrix1_block[i * block_size + j] =
            vectored_matrix1[(block_row * block_size + i) * block_size +
                             block_col * block_size + j];
        vectored_matrix2_block[i * block_size + j] =
            vectored_matrix2[(block_row * block_size + i) * block_size +
                             block_col * block_size + j];
        vectored_result_block[i * block_size + j] = 0;
      }
    }

    int rank_left, rank_right, rank_up, rank_down;
    MPI_Cart_shift(comm_2d_cart, 0, -1, &rank_right, &rank_left);
    MPI_Cart_shift(comm_2d_cart, 1, -1, &rank_down, &rank_up);

    int rank_source, rank_dest;
    MPI_Cart_shift(comm_2d_cart, 0, -comm_2d_cart_coords[0], &rank_source,
                   &rank_dest);
    MPI_Sendrecv_replace(matrix1, );

    elapsed_time_type begin = MPI_Wtime();
    for (int block_id = 0; block_id < block_matrix_dim; ++block_id) {
      size_t i_start = cur_block_matrix_row * block_size;
      size_t j_start = cur_block_matrix_col * block_size;
      for (size_t i = 0; i < block_size; ++i) {
        for (size_t j = 0; j < block_size; ++j) {
          for (size_t k = 0; k < block_size; ++k) {
            result[i_start + i][j_start + j] +=
                matrix1[i_start + i][j_start + k] *
                matrix2[i_start + k][j_start + j];
          }
        }
      }
      MPI_Status status;
      MPI_Sendrecv_replace(matrix1, nlocal * nlocal, MPI_DOUBLE, left, 1, right,
                           1, comm_2d_cart, &status);
      MPI_Sendrecv_replace(matrix2, nlocal * nlocal, MPI_DOUBLE, up, 1, down, 1,
                           comm_2d_cart, &status);
    }
    MPI_Barrier(comm_2d_cart);

    elapsed_time_type end = MPI_Wtime();

    MPI_Comm_free(comm_2d_cart);

    int finalized;
    MPI_Finalized(&finalized);
    if (!finalized) {
      MPI_Finalize();
    }
    return end - begin;
  }
