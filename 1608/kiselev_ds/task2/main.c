#include <stdio.h>
#include <math.h>
#include "mpi.h"

#define EPSILON 10e-7

double* read_matrix(FILE* input, int x, int y) {
	int vector_size = y * x;
	double* matrix = (double*)malloc(sizeof(double) * vector_size);
	for (int i = 0; (i < vector_size) && (!feof(input)); i++)
		fscanf(input, "%lf", &matrix[i]);
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
	for (int i = 0; i < size; i++)
		vector[i] *= value;
}

void divide_vector2value(double* vector, int size, double value) {
	for (int i = 0; i < size; i++)
		vector[i] /= value;
}

void diff_vector2vector_with_coef(double* v1, double* v2, int size, double coef) {
	for (int i = 0; i < size; i++)
		v1[i] -= v2[i] * coef;
}

void print_matrix(double* matrix, int x, int y) {
	for (int i = 0; i < y; i++) {
		for (int j = 0; j < x; j++)
			printf("%.1lf\t", matrix[i * x + j]);
		printf("\n");
	}
}

int main(int argc, char *argv[]) {
	int proc_num, proc_rank, proc_num_used;
	int *blocklen, *displacement;
	int division_size, rest_size;
	int matrix_h, matrix_w;
	int submatr_h, submatr_w;
	double* matrix;
	MPI_Status status;

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
		submatr_h = matrix_h;
		submatr_w = matrix_w;
		MPI_Datatype* indexed_types;

		for (int i = 0; i < matrix_h; i++, submatr_h--, submatr_w--) {
			proc_num_used = (proc_num <= submatr_h - 1) ? proc_num : submatr_h - 1;
			diag_offset = matrix_w * i + i;
			
			// Swap rows if first element of sub-matrix is zero
			if (fabs(matrix[diag_offset]) < EPSILON) {
				for (int k = i + 1; k < matrix_h; k++)
					if (fabs(matrix[matrix_w * k + i]) >= EPSILON) {
						swap_vectors(&matrix[matrix_w * k], &matrix[matrix_w * i], matrix_w);
						break;
					}
			}

			if (fabs(matrix[diag_offset] - 1) >= EPSILON) {
				double value = matrix[diag_offset];
				divide_vector2value(&matrix[diag_offset], submatr_w, matrix[diag_offset]);
			}

			// The rest_size is needed to handle the case when matrix size % proc_num != 0
			division_size = (int)ceil((submatr_h - 2.0) / proc_num_used);
			rest_size = (submatr_h - 1) - (proc_num_used - 1) * division_size;

			indexed_types = (MPI_Datatype*)malloc(sizeof(MPI_Datatype) * (proc_num_used - 1));

			blocklen = (int*)malloc(sizeof(int) * (division_size + 1));
			displacement = (int*)malloc(sizeof(int) * (division_size + 1));

			for (int j = 0; j < division_size + 1; j++)
				blocklen[j] = submatr_w;

			displacement[0] = diag_offset;
			for (int proc = 1; proc < proc_num_used; proc++) {
				displacement[1] = displacement[0] + matrix_w * (rest_size + 1 + (proc - 1) * division_size);

				for (int m = 2; m < division_size + 1; m++)
					displacement[m] = displacement[m - 1] + matrix_w;

				MPI_Type_indexed(division_size + 1, blocklen, displacement, MPI_DOUBLE, &indexed_types[proc - 1]);
				MPI_Type_commit(&indexed_types[proc - 1]);
				MPI_Send(matrix, 1, indexed_types[proc - 1], proc, 0, MPI_COMM_WORLD);
			}

			for (int m = 0; m < rest_size; m++) {
				double coef = matrix[diag_offset + (m + 1) * matrix_w];
				if (fabs(coef) >= EPSILON) {
					diff_vector2vector_with_coef(&matrix[diag_offset + (m + 1) * matrix_w],
						&matrix[diag_offset], submatr_w, coef);
				}
			}

			for (int proc = 1; proc < proc_num_used; proc++)
				MPI_Recv(matrix, 1, indexed_types[proc - 1], proc, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			free(blocklen);
			free(displacement);
			free(indexed_types);
		}

		print_matrix(matrix, matrix_w, matrix_h);
		printf("\n");
	}
	else {
		// Getting the size of source matrix
		MPI_Bcast(&matrix_h, 1, MPI_INT, 0, MPI_COMM_WORLD);
		matrix_w = matrix_h + 1;
		matrix = (double*)malloc(sizeof(double) * matrix_h * matrix_w);

		// Calculate how many times the root process will
		// use this process during the forward or reverse course
		int actions = matrix_h - proc_rank - 1;

		int diag_offset;
		submatr_h = matrix_h;
		submatr_w = matrix_w;
		MPI_Datatype type;

		for (int i = 0; i < actions; i++, submatr_w--, submatr_h--) {
			proc_num_used = (proc_num <= submatr_h - 1) ? proc_num : submatr_h - 1;
			division_size = (int)ceil((submatr_h - 2.0) / proc_num_used);
			rest_size = (submatr_h - 1) - (proc_num_used - 1) * division_size;
			diag_offset = matrix_w * i + i;

			blocklen = (int*)malloc(sizeof(int) * (division_size + 1));
			displacement = (int*)malloc(sizeof(int) * (division_size + 1));

			blocklen[0] = blocklen[1] = submatr_w;
			displacement[0] = diag_offset;
			displacement[1] = displacement[0] + matrix_w * (rest_size + 1 + (proc_rank - 1) * division_size);
			for (int m = 2; m < division_size + 1; m++) {
				displacement[m] = displacement[m - 1] + matrix_w;
				blocklen[m] = submatr_w;
			}

			MPI_Type_indexed(division_size + 1, blocklen, displacement, MPI_DOUBLE, &type);
			MPI_Type_commit(&type);
			MPI_Recv(matrix, 1, type, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			for (int i = 0; i < division_size; i++) {
				double coef = matrix[displacement[i + 1]];
				if (fabs(coef) >= EPSILON) {
					diff_vector2vector_with_coef(&matrix[displacement[i + 1]], &matrix[displacement[0]], submatr_w, coef);
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
