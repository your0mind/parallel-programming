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

int main(int argc, char *argv[]) {
	int proc_num, proc_rank;
	MPI_Datatype submatrix_type;
	MPI_Status status;
	int blocklen[2], displacement[2];
	int matrix_h, matrix_w;
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

		for (int i = 0; i < 1; i++) {
			// Swap rows if first element of sub-matrix is zero
			if (fabs(matrix[matrix_w * i + i]) < EPSILON) {
				for (int k = 0; i < matrix_h; k++) {
					if (fabs(matrix[matrix_w * k + i]) >= EPSILON) {
						swap_vectors(matrix + i * matrix_w, matrix + k * matrix_w, matrix_w);
						break;
					}
				}
			}

			if (fabs(matrix[matrix_w * i + i] - 1) >= EPSILON) {
				divide_vector2value(matrix + i * matrix_w, matrix_w, matrix[matrix_w * i + i]);
			}
			
			int matrix_h_cur = matrix_h - i;
			int proc_num_cur = (proc_num <= matrix_h_cur - 1) ? proc_num : matrix_h_cur - 1;
			// The rest_size is needed to handle the case when matrix size % proc_num != 0
			int division_size = (int)ceil((matrix_h_cur - 2.0) / proc_num_cur);
			int rest_size = (matrix_h_cur - 1) - (proc_num - 1) * division_size;

			MPI_Datatype* types = (MPI_Datatype*)malloc(sizeof(MPI_Datatype) * (proc_num_cur - 1));
			for (int proc = 1; proc < proc_num_cur; proc++) {
				blocklen[0] = matrix_w;		
				blocklen[1] = matrix_w * division_size;
				displacement[0] = i * matrix_h;
				displacement[1] = displacement[0] + matrix_w * (rest_size + 1 + (proc - 1) * division_size);
				MPI_Type_indexed(2, blocklen, displacement, MPI_DOUBLE, types + proc - 1);
				MPI_Type_commit(types + proc - 1);
				MPI_Send(matrix, 1, types[proc - 1], 1, 0, MPI_COMM_WORLD);
			}

			free(types);
		}
	}
	else {
		// Getting the size of source matrix
		MPI_Bcast(&matrix_h, 1, MPI_INT, 0, MPI_COMM_WORLD);
		matrix_w = matrix_h + 1;

		MPI_Status status;
		MPI_Datatype type;
		matrix = (double*)malloc(sizeof(double) * matrix_h * matrix_w);

		blocklen[0] = 4;
		blocklen[1] = 4;
		displacement[0] = 0;
		displacement[1] = 8;
		MPI_Type_indexed(2, blocklen, displacement, MPI_DOUBLE, &type);
		MPI_Type_commit(&type);
		MPI_Recv(matrix, 1, type, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		for (int i = 0; i < matrix_h * matrix_w; i++) {
			printf("%lf ", matrix[i]);
		}
		
	}

	//free(matrix);
	MPI_Finalize();
	return 0;
}
