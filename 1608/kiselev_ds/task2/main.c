#include <stdio.h>
#include "mpi.h"

#define EPSILON 10e-7

double* read_matrix(FILE* input, int x, int y) {
	int vector_size = y * x;
	double* matrix = (double**)malloc(sizeof(double) * vector_size);
	for (int i = 0; (i < vector_size) && (!feof(input)); i++) {
		fscanf(input, "%lf", &matrix[i]);
	}
	return matrix;
}

void swap_rows(double* matrix, int y, int x, int row1, int row2) {
	for (int i = 0; i < x; i++) {
		double temp = matrix[y * (row1 - 1) + i];
		matrix[y * (row1 - 1) + i] = matrix[y * (row2 - 1) + i];
		matrix[y * (row2 - 1) + i] = temp;
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
	int blocklen[2] = { 4, 4 };
	int displacement[2] = { 0, 8 };
	int matrix_size;
	double* matrix;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	MPI_Type_indexed(2, blocklen, displacement, MPI_DOUBLE, &submatrix_type);
	MPI_Type_commit(&submatrix_type);

	if (proc_rank == 0) {
		FILE* input_file = fopen(argv[1], "r");
		fscanf(input_file, "%d", &matrix_size);

		// Sending matrix_size to all processes
		MPI_Bcast(&matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// Reading matrix from input file
		matrix = read_matrix(input_file, matrix_size, matrix_size + 1);
		fclose(input_file);

		MPI_Send(matrix, 1, submatrix_type, 1, 0, MPI_COMM_WORLD);

		//for (int i = 0; i < 0; i++) {
		//	// Swap rows if first element of sub-matrix is zero
		//	if (fabs(matrix[i][i]) < EPSILON) {
		//		for (int k = i + 1; i < matrix_size; k++) {
		//			if (fabs(matrix[k][i]) >= EPSILON) {
		//				swap_rows(matrix, i, k);
		//				break;
		//			}
		//		}
		//	}

		//	if (fabs(matrix[i][i] - 1) >= EPSILON) {
		//		divide_vector2value(matrix[i], matrix_size + 1, matrix[i][i]);
		//	}

		//	int curr_matrix_size = matrix_size - i;
		//	// The rest_size is needed to handle the case when matrix size % proc_num != 0
		//	int devision_size = (curr_matrix_size - 2) / proc_num + 1;
		//	int rest_size = (curr_matrix_size - 1) - (proc_num - 1) * devision_size;

		//}
	}
	else {
		// Getting the size of source matrix
		MPI_Bcast(&matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

		MPI_Status status;
		matrix = (double*)malloc(sizeof(double) * matrix_size * (matrix_size + 1));
		MPI_Recv(matrix, 1, submatrix_type, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		
		for ()
	}

	//free(matrix);
	MPI_Finalize();
	return 0;
}
