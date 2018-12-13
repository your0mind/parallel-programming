// Lab2: Gauss–Jordan elimination

#include <stdio.h>
#include <math.h>
#include "mpi.h"

#define EPSILON 10e-7

double* generate_simple_matrix(int height, int width) {
	double *matrix = (double*)malloc(sizeof(double) * height * width);
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
			if (y == x)
				matrix[y * width + x] = 2.0;
			else if (x == height)
				matrix[y * width + x] = height + 1.0;
			else
				matrix[y * width + x] = 1.0;
	return matrix;
}

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

void substract_vectors_with_coef(double* v1, double* v2, int size, double coef) {
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
	int proc_num, proc_rank;
	int *division_sizes, *rest_sizes, *proc_num_used;
	int *blocklen, *displacement;
	int matrix_h, matrix_w;
	int subm_h, subm_w;
	int elems_in_type;
	double *matrix;
	double start_time, end_time;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	if (proc_rank == 0) {
		if (argc == 3) {
			FILE* input_file = fopen(argv[1], "r");
			fscanf(input_file, "%d", &matrix_h);
			matrix_w = matrix_h + 1;

			// Reading matrix from input file
			matrix = read_matrix(input_file, matrix_h, matrix_w);
			fclose(input_file);
		}
		else {
			// Generatin simple matrix
			matrix_h = atoi(argv[1]);
			matrix_w = matrix_h + 1;
			matrix = generate_simple_matrix(matrix_h, matrix_w);
		}

		double *matrix_copy = (double*)malloc(sizeof(double) * matrix_h * matrix_w);
		memcpy(matrix_copy, matrix, matrix_h * matrix_w * sizeof(double));
	start:
		if (matrix == NULL)
			matrix = matrix_copy;

		start_time = MPI_Wtime();

		subm_h = matrix_h;
		subm_w = matrix_w;
		int diag_offset = 0;
		for (int i = 0; i < matrix_h; i++, subm_h--, subm_w--, diag_offset += matrix_w + 1) {
			if (fabs(matrix_copy[diag_offset]) < EPSILON)
				for (int j = i + 1; j < matrix_h; j++)
					if (fabs(matrix_copy[matrix_w * j + i]) >= EPSILON) {
						swap_vectors(&matrix_copy[matrix_w * j], &matrix_copy[matrix_w * i], matrix_w);
						break;
					}

			if (fabs(matrix_copy[diag_offset] - 1) >= EPSILON)
				divide_vector2value(&matrix_copy[diag_offset], subm_w, matrix_copy[diag_offset]);

			for (int j = 0; j < subm_h - 1; j++) {
				double coef = matrix_copy[diag_offset + (j + 1) * matrix_w];
				if (fabs(coef) >= EPSILON)
					substract_vectors_with_coef(&matrix_copy[diag_offset + (j + 1) * matrix_w],
						&matrix_copy[diag_offset], subm_w, coef);
			}
		}

		//printf("\nForward elimination's result:\n");
		//print_matrix(matrix_copy, matrix_w, matrix_h);
		//printf("\n\n");

		int pivot_elem_offset = matrix_h * matrix_w - 1;
		for (int i = 0; i < matrix_h - 1; i++, pivot_elem_offset -= matrix_w) {
			int start_of_col = matrix_w - 1;
			for (int k = 0; k < matrix_h - i - 1; k++, start_of_col += matrix_w) {
				matrix_copy[start_of_col] -= matrix_copy[pivot_elem_offset] * matrix_copy[start_of_col - (i + 1)];
			}
		}

		end_time = MPI_Wtime();
		printf("Sequential %f sec\n", end_time - start_time);
		FILE* output_file = fopen("output.txt", "w");
		for (int i = 0; i < matrix_h; i++)
			fprintf(output_file, "%lf\n", matrix_copy[(matrix_w - 1) + matrix_w * i]);
		fclose(output_file);

		//printf("\nBack substitution's result:\n");
		//print_matrix(matrix_copy, matrix_w, matrix_h);
		//printf("\n\n");

		// Parallel algorithm
		start_time = MPI_Wtime();

		// Sending matrix_size to all processes
		MPI_Bcast(&matrix_h, 1, MPI_INT, 0, MPI_COMM_WORLD);

		division_sizes = (int*)malloc(sizeof(int) * (matrix_h - 1));
		rest_sizes = (int*)malloc(sizeof(int) * (matrix_h - 1));
		proc_num_used = (int*)malloc(sizeof(int) * (matrix_h - 1));

		subm_h = matrix_h;
		subm_w = matrix_w;
		MPI_Datatype* indexed_types;

		// Forward elimination steps
		diag_offset = 0;
		for (int i = 0; i < matrix_h - 1; i++, subm_h--, subm_w--, diag_offset += matrix_w + 1) {
			diag_offset = matrix_w * i + i;
			// Calculate the number of processes that we will use on the current iteration
			proc_num_used[i] = (proc_num <= subm_h - 1) ? proc_num : subm_h - 1;
			
			// Swap rows if first element of sub-matrix is zero
			if (fabs(matrix[diag_offset]) < EPSILON)
				for (int j = i + 1; j < matrix_h; j++)
					if (fabs(matrix[matrix_w * j + i]) >= EPSILON) {
						swap_vectors(&matrix[matrix_w * j], &matrix[matrix_w * i], matrix_w);
						break;
					}

			// Divide the row into itself if first element != 1
			if (fabs(matrix[diag_offset] - 1) >= EPSILON)
				divide_vector2value(&matrix[diag_offset], subm_w, matrix[diag_offset]);

			// The rest size is needed to handle the case when matrix size % proc_num != 0
			division_sizes[i] = (subm_h - 1) / proc_num_used[i];
			rest_sizes[i] = (subm_h - 1) - (proc_num_used[i] - 1) * division_sizes[i];

			// We will store there the mpi indexed types that we create
			// in order to transfer the matrix parts to the processes
			indexed_types = (MPI_Datatype*)malloc(sizeof(MPI_Datatype) * (proc_num_used[i] - 1));

			// Creating the necessary blocklens and displacements for indexed types
			elems_in_type = division_sizes[i] + 1;
			blocklen = (int*)malloc(sizeof(int) * elems_in_type);
			displacement = (int*)malloc(sizeof(int) * elems_in_type);

			for (int j = 0; j < elems_in_type; j++)
				blocklen[j] = subm_w;

			displacement[0] = diag_offset;
			for (int proc = 1; proc < proc_num_used[i]; proc++) {
				displacement[1] = displacement[0] + matrix_w * (rest_sizes[i] + 1 + (proc - 1) * division_sizes[i]);
				for (int j = 2; j < elems_in_type; j++)
					displacement[j] = displacement[j - 1] + matrix_w;

				MPI_Type_indexed(division_sizes[i] + 1, blocklen, displacement, MPI_DOUBLE, &indexed_types[proc - 1]);
				MPI_Type_commit(&indexed_types[proc - 1]);
				// Send the matrix parts to the processes
				MPI_Send(matrix, 1, indexed_types[proc - 1], proc, 0, MPI_COMM_WORLD);
			}

			// Processing own part of matrix
			for (int j = 0; j < rest_sizes[i]; j++) {
				double coef = matrix[diag_offset + (j + 1) * matrix_w];
				if (fabs(coef) >= EPSILON)
					// Subtract from the current row the first row multiplied by first element
					substract_vectors_with_coef(&matrix[diag_offset + (j + 1) * matrix_w],
						&matrix[diag_offset], subm_w, coef);
			}

			// Getting updated matrix and release indexed types
			for (int proc = 1; proc < proc_num_used[i]; proc++) {
				MPI_Recv(matrix, 1, indexed_types[proc - 1], proc, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				MPI_Type_free(&indexed_types[proc - 1]);
			}

			free(blocklen);
			free(displacement);
			free(indexed_types);
		}

		matrix[diag_offset + 1] /= matrix[diag_offset];
		matrix[diag_offset] = 1.0;

		//printf("\nForward elimination's result:\n");
		//print_matrix(matrix, matrix_w, matrix_h);
		//printf("\n\n");

		// Back substitution
		pivot_elem_offset = matrix_h * matrix_w - 1;
		for (int i = 0; i < matrix_h - 1; i++, pivot_elem_offset -= matrix_w) {
			indexed_types = (MPI_Datatype*)malloc(sizeof(MPI_Datatype) * (proc_num_used[i] - 1));

			elems_in_type = division_sizes[i] * 2 + 1;
			blocklen = (int*)malloc(sizeof(int) * elems_in_type);
			displacement = (int*)malloc(sizeof(int) * elems_in_type);
			displacement[elems_in_type - 1] = pivot_elem_offset;

			for (int j = 0; j < elems_in_type; j++)
				blocklen[j] = 1;

			for (int proc = 1; proc < proc_num_used[i]; proc++) {
				for (int j = 0; j < elems_in_type - 1; j += 2) {
					displacement[j] = (matrix_w - i - 2 + division_sizes[i] * matrix_w * (proc - 1)) + ((j >> 1)) * matrix_w;
					displacement[j + 1] = displacement[j] + i + 1;
				}
				MPI_Type_indexed(elems_in_type, blocklen, displacement, MPI_DOUBLE, &indexed_types[proc - 1]);
				MPI_Type_commit(&indexed_types[proc - 1]);
				MPI_Send(matrix, 1, indexed_types[proc - 1], proc, 0, MPI_COMM_WORLD);
			}

			int start_of_col = pivot_elem_offset - matrix_w * rest_sizes[i];
			for (int k = 0; k < rest_sizes[i]; k++, start_of_col += matrix_w) {
				matrix[start_of_col] -= matrix[pivot_elem_offset] * matrix[start_of_col - (i + 1)];
				//matrix[start_of_col - (i + 1)] = 0.0;
			}

			for (int proc = 1; proc < proc_num_used[i]; proc++) {
				MPI_Recv(matrix, 1, indexed_types[proc - 1], proc, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				MPI_Type_free(&indexed_types[proc - 1]);
			}

			free(blocklen);
			free(displacement);
			free(indexed_types);
		}

		//printf("\nBack substitution's result:\n");
		//print_matrix(matrix, matrix_w, matrix_h);
		//printf("\n\n");

		if (argc == 3) {
			FILE* output_file = fopen(argv[2], "w");
			for (int i = 0; i < matrix_h; i++)
				fprintf(output_file, "%lf\n", matrix[(matrix_w - 1) + matrix_w * i]);
			fclose(output_file);
		}

		end_time = MPI_Wtime();
		printf("Parallel: %f sec\n", end_time - start_time);
		
		int successful_validation = 1;
		for (int last_elem = matrix_w - 1; last_elem < matrix_h * matrix_w; last_elem += matrix_w)
			if (fabs(matrix_copy[last_elem] - matrix[last_elem]) > EPSILON)
				successful_validation = 0;
		(successful_validation)
			? printf("SUCCESSFUL VALIDATION CHECK\n")
			: printf("UNSUCCESSFUL VALIDATION CHECK\n");
	}
	else {
		// Getting the size of source matrix
		MPI_Bcast(&matrix_h, 1, MPI_INT, 0, MPI_COMM_WORLD);
		matrix_w = matrix_h + 1;
		matrix = (double*)malloc(sizeof(double) * matrix_h * matrix_w);
		double* pivot_row = (double*)malloc(sizeof(double) * matrix_w);

		// Calculate how many times the root process will use this 
		// process during the forward elimination or back substitution
		int actions = matrix_h - proc_rank - 1;

		division_sizes = (int*)malloc(sizeof(int) * actions);
		rest_sizes = (int*)malloc(sizeof(int) * actions);
		proc_num_used = (int*)malloc(sizeof(int) * actions);

		int diag_offset;
		subm_h = matrix_h;
		subm_w = matrix_w;
		MPI_Datatype indexed_type;
		// Forward elimination steps
		for (int i = 0; i < actions; i++, subm_w--, subm_h--) {
			proc_num_used[i] = (proc_num <= subm_h - 1) ? proc_num : subm_h - 1;
			division_sizes[i] = (subm_h - 1) / proc_num_used[i];
			rest_sizes[i] = (subm_h - 1) - (proc_num_used[i] - 1) * division_sizes[i];
			diag_offset = matrix_w * i + i;

			elems_in_type = division_sizes[i] + 1;
			blocklen = (int*)malloc(sizeof(int) * elems_in_type);
			displacement = (int*)malloc(sizeof(int) * elems_in_type);

			blocklen[0] = blocklen[1] = subm_w;
			displacement[0] = diag_offset;
			displacement[1] = displacement[0] + matrix_w * ((rest_sizes[i] + 1) + (proc_rank - 1) * division_sizes[i]);
			for (int j = 2; j < elems_in_type; j++) {
				displacement[j] = displacement[j - 1] + matrix_w;
				blocklen[j] = subm_w;
			}

			MPI_Type_indexed(division_sizes[i] + 1, blocklen, displacement, MPI_DOUBLE, &indexed_type);
			MPI_Type_commit(&indexed_type);
			MPI_Recv(matrix, 1, indexed_type, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			for (int k = 0; k < division_sizes[i]; k++) {
				double coef = matrix[displacement[k + 1]];
				if (fabs(coef) >= EPSILON)
					substract_vectors_with_coef(&matrix[displacement[k + 1]], &matrix[displacement[0]], subm_w, coef);
			}

			MPI_Send(matrix, 1, indexed_type, 0, 0, MPI_COMM_WORLD);
			MPI_Type_free(&indexed_type);
			free(blocklen);
			free(displacement);
		}

		int pivot_elem_offset = matrix_h * matrix_w - 1;
		// Back substitution
		for (int i = 0; i < actions; i++, pivot_elem_offset -= matrix_w) {
			elems_in_type = division_sizes[i] * 2 + 1;
			blocklen = (int*)malloc(sizeof(int) * (division_sizes[i] * 2 + 1));
			displacement = (int*)malloc(sizeof(int) * (division_sizes[i] * 2 + 1));
			displacement[elems_in_type - 1] = pivot_elem_offset;

			blocklen[elems_in_type - 1] = 1;
			for (int j = 0; j < elems_in_type - 1; j += 2) {
				blocklen[j] = blocklen[j + 1] = 1;
				displacement[j] = (matrix_w - i - 2 + division_sizes[i] * matrix_w * (proc_rank - 1)) + ((j >> 1)) * matrix_w;
				displacement[j + 1] = displacement[j] + i + 1;
			}
			MPI_Type_indexed(elems_in_type, blocklen, displacement, MPI_DOUBLE, &indexed_type);
			MPI_Type_commit(&indexed_type);
			MPI_Recv(matrix, 1, indexed_type, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			for (int k = 0; k < elems_in_type - 1; k += 2) {
				matrix[displacement[k + 1]] -= matrix[displacement[k]] * matrix[pivot_elem_offset];
				//matrix[displacement[k]] = 0.0;
			}

			MPI_Send(matrix, 1, indexed_type, 0, 0, MPI_COMM_WORLD);
			MPI_Type_free(&indexed_type);
			free(blocklen);
			free(displacement);
		}
	}

	free(division_sizes);
	free(rest_sizes);
	free(proc_num_used);
	free(matrix);
	free(matrix_copy);
	MPI_Finalize();
	return 0;
}
