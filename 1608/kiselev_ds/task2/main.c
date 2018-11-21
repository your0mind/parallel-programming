#include <stdio.h>
#include <math.h>
#include "mpi.h"

#define EPSILON 10e-7

double* read_matrix(FILE* input, int x, int y) {
	int vector_size = y * x;
	double* matrix = (double*)malloc(sizeof(double) * vector_size);
	for (int i = 0; (i < vector_size) && (!feof(input)); i++) {
		fscanf(input, "%lf", &matrix[i]);
	}
	return matrix;
}

void swap_vectors(double* v1, double* v2, int size) {
	for (int i = 0; i < size; i++) {
		double temp = v1[i];
		v1[i] = v2[i];
		v2[i] = temp;
	}
}

void multiply_vector2value(double* vector, int size, double value) {
	for (int i = 0; i < size; i++) {
		vector[i] *= value;
	}
}

void divide_vector2value(double* vector, int size, double value) {
	for (int i = 0; i < size; i++) {
		vector[i] /= value;
	}
}

void diff_vector2vector_with_coef(double* v1, double* v2, int size, double coef) {
	for (int i = 0; i < size; i++) {
		v1[i] -= v2[i] * coef;
	}
}

void vector_copy(double* vector, double* dist, int size) {

}

int main(int argc, char *argv[]) {
	int proc_num, proc_rank, proc_num_cur;
	MPI_Status status;
	int *blocklen, *displacement;
	int matrix_h, matrix_w;
	int division_size, rest_size;
	double* matrix;



	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	if (proc_rank == 0) {
		FILE* input_file = fopen(argv[1], "r");
		fscanf(input_file, "%d", &matrix_h);

		// Sending matrix_size to all processes
		MPI_Bcast(&matrix_h, 1, MPI_INT, 0, MPI_COMM_WORLD);
		matrix_w = matrix_h + 1;

		// Reading matrix from input file
		matrix = read_matrix(input_file, matrix_h, matrix_w);
		fclose(input_file);

		int diag_offset;
		MPI_Datatype* types;
		for (int i = 0, w_cur = matrix_w, h_cur = matrix_h; i < matrix_h; i++, w_cur--, h_cur--) {
			diag_offset = matrix_w * i + i;
			proc_num_cur = (proc_num <= h_cur - 1) ? proc_num : h_cur - 1;

			// Swap rows if first element of sub-matrix is zero
			if (fabs(matrix[diag_offset]) < EPSILON) {
				for (int k = i + 1; k < matrix_h; k++) {
					if (fabs(matrix[matrix_w * k + i]) >= EPSILON) {
						swap_vectors(&matrix[matrix_w * k], &matrix[matrix_w * i], matrix_w);
						break;
					}
				}
			}

			if (fabs(matrix[diag_offset] - 1) >= EPSILON) {
				double value = matrix[diag_offset];
				divide_vector2value(&matrix[diag_offset], w_cur, matrix[diag_offset]);
			}

			// The rest_size is needed to handle the case when matrix size % proc_num != 0
			division_size = (int)ceil((h_cur - 2.0) / proc_num_cur);
			rest_size = (h_cur - 1) - (proc_num_cur - 1) * division_size;

			types = (MPI_Datatype*)malloc(sizeof(MPI_Datatype) * (proc_num_cur - 1));

			blocklen = (int*)malloc(sizeof(int) * (division_size + 1));
			displacement = (int*)malloc(sizeof(int) * (division_size + 1));

			for (int j = 0; j < division_size + 1; j++) {
				blocklen[j] = w_cur;
			}
			displacement[0] = diag_offset;
			for (int proc = 1; proc < proc_num_cur; proc++) {
				displacement[1] = displacement[0] + matrix_w * (rest_size + 1 + (proc - 1) * division_size);
				for (int m = 2; m < division_size + 1; m++) {
					displacement[m] = displacement[m - 1] + matrix_w;
				}
				MPI_Type_indexed(division_size + 1, blocklen, displacement, MPI_DOUBLE, &types[proc - 1]);
				MPI_Type_commit(&types[proc - 1]);
				MPI_Send(matrix, 1, types[proc - 1], proc, 0, MPI_COMM_WORLD);
			}

			for (int m = 0; m < rest_size; m++) {
				double coef = matrix[diag_offset + (m + 1) * matrix_w];
				if (fabs(coef) >= EPSILON) {
					diff_vector2vector_with_coef(&matrix[diag_offset + (m + 1) * matrix_w],
						&matrix[diag_offset], w_cur, coef);
				}
			}

			for (int proc = 1; proc < proc_num_cur; proc++) {
				MPI_Recv(matrix, 1, types[proc - 1], proc, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}

			free(blocklen);
			free(displacement);
			free(types);
		}

		for (int i = 0; i < matrix_h * matrix_w; i++) {
			printf("%.1lf ", matrix[i]);
			if ((i + 1) % matrix_w == 0) {
				printf("\n");
			}
		}
		printf("\n");
	}
	else {
		// Getting the size of source matrix
		MPI_Bcast(&matrix_h, 1, MPI_INT, 0, MPI_COMM_WORLD);
		matrix_w = matrix_h + 1;
		matrix = (double*)malloc(sizeof(double) * matrix_h * matrix_w);

		MPI_Status status;
		MPI_Datatype type;
		int iterations = matrix_h - proc_rank - 1;

		int diag_offset;
		for (int i = 0, w_cur = matrix_w, h_cur = matrix_h; i < iterations; i++, w_cur--, h_cur--) {
			diag_offset = matrix_w * i + i;
			proc_num_cur = (proc_num <= h_cur - 1) ? proc_num : h_cur - 1;
			division_size = (int)ceil((h_cur - 2.0) / proc_num_cur);
			rest_size = (h_cur - 1) - (proc_num_cur - 1) * division_size;

			blocklen = (int*)malloc(sizeof(int) * (division_size + 1));
			displacement = (int*)malloc(sizeof(int) * (division_size + 1));

			blocklen[0] = blocklen[1] = w_cur;
			displacement[0] = diag_offset;
			displacement[1] = displacement[0] + matrix_w * (rest_size + 1 + (proc_rank - 1) * division_size);
			for (int m = 2; m < division_size + 1; m++) {
				displacement[m] = displacement[m - 1] + matrix_w;
				blocklen[m] = w_cur;
			}

			MPI_Type_indexed(division_size + 1, blocklen, displacement, MPI_DOUBLE, &type);
			MPI_Type_commit(&type);
			MPI_Recv(matrix, 1, type, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			for (int i = 0; i < division_size; i++) {
				double coef = matrix[displacement[i + 1]];
				if (fabs(coef) >= EPSILON) {
					diff_vector2vector_with_coef(&matrix[displacement[i + 1]], &matrix[displacement[0]], w_cur, coef);
				}
			}

			MPI_Send(matrix, 1, type, 0, 0, MPI_COMM_WORLD);

			free(blocklen);
			free(displacement);
		}
	}

	free(matrix);
	MPI_Finalize();
	return 0;
}
