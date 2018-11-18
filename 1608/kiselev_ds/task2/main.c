#include <stdio.h>
#include "mpi.h"

#define EPSILON 10e-7

double** read_matrix(FILE* input, int x, int y) {
	double** matrix = (double**)malloc(sizeof(double*) * y);
	for (int i = 0; i < y; i++) {
		matrix[i] = (double*)malloc(sizeof(double) * x);
	}
	for (int i = 0; (i < y) && (!feof(input)); i++) {
		for (int j = 0; (j < x + 1) && (!feof(input)); j++) {
			fscanf(input, "%lf", &matrix[i][j]);
		}
	}
	return matrix;
}

void swap_rows(double** matrix, int row1, int row2) {
	double* temp = matrix[row1];
	matrix[row1] = matrix[row2];
	matrix[row2] = temp;
}

int main(int argc, char *argv[]) {
	int proc_num, proc_rank;
	int matrix_size;
	double** matrix;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	if (proc_rank == 0) {
		// Creating matrix from input file
		FILE* input_file = fopen(argv[1], "r");
		fscanf(input_file, "%d", &matrix_size);
		matrix = read_matrix(input_file, matrix_size, matrix_size + 1);

		for (int i = 0; i < matrix_size; i++) {
			// Swap rows if first element of sub-matrix is zero
			if (matrix[i][i] == 0) {
				for (int k = i + 1; i < matrix_size; k++) {
					if (fabs(matrix[k][i]) < EPSILON) {
						swap_rows(matrix, i, k);
						break;
					}
				}
			}

			int curr_matrix_size = matrix_size - i;
			// The rest_size is needed to handle the case when matrix size % proc_num != 0
			int devision_size = (curr_matrix_size - 1) / proc_num + 1;
			int rest_size = curr_matrix_size - (proc_num - 1) * devision_size;


		}
	}

	//free(matrix);
	MPI_Finalize();
	return 0;
}