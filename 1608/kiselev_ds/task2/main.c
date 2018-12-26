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

void constr_row_and_subm_type(
	MPI_Datatype *source_type,
	int matrix_w,
	int row_offset,
	int subm_offset,
	int subm_h,
	int subm_w
) {
	// Creating the necessary blocklens and displacements for indexed types
	int elems_in_type = subm_h + 1;
	int *blocklen = (int*)malloc(sizeof(int) * elems_in_type);
	int *displacement = (int*)malloc(sizeof(int) * elems_in_type);

	for (int i = 0; i < elems_in_type; i++)
		blocklen[i] = subm_w;
	displacement[0] = row_offset;
	for (int j = 0; j < subm_h; j++)
		displacement[j + 1] = subm_offset + matrix_w * j;

	MPI_Type_indexed(elems_in_type, blocklen, displacement, MPI_DOUBLE, source_type);
	MPI_Type_commit(source_type);
	free(blocklen);
	free(displacement);
}

void constr_cols_and_elem_type(
	MPI_Datatype *source_type,
	int matrix_w,
	int pivot_elem_offset,
	int first_col_offset,
	int second_col_offset,
	int cols_h
) {
	int elems_in_type = cols_h * 2 + 1;
	int *blocklen = (int*)malloc(sizeof(int) * elems_in_type);
	int *displacement = (int*)malloc(sizeof(int) * elems_in_type);

	for (int j = 0; j < elems_in_type; j++)
		blocklen[j] = 1;
	for (int j = 0; j < elems_in_type - 1; j += 2) {
		int jump_size = matrix_w * (j >> 1);
		displacement[j] = first_col_offset + jump_size;
		displacement[j + 1] = second_col_offset + jump_size;
	}
	displacement[elems_in_type - 1] = pivot_elem_offset;

	MPI_Type_indexed(elems_in_type, blocklen, displacement, MPI_DOUBLE, source_type);
	MPI_Type_commit(source_type);
	free(blocklen);
	free(displacement);
}

int main(int argc, char *argv[]) {
	int proc_num, proc_rank;
	int *div_sizes, *rest_sizes, *proc_num_used;
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
		
		// Sequential algorithm
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
		double sequential_time = end_time - start_time;

		// Parallel algorithm
		start_time = MPI_Wtime();

		// Sending matrix_size to all processes
		MPI_Bcast(&matrix_h, 1, MPI_INT, 0, MPI_COMM_WORLD);

		div_sizes = (int*)malloc(sizeof(int) * (matrix_h - 1));
		rest_sizes = (int*)malloc(sizeof(int) * (matrix_h - 1));
		proc_num_used = (int*)malloc(sizeof(int) * (matrix_h - 1));

		subm_h = matrix_h;
		subm_w = matrix_w;
		MPI_Datatype* spec_types;

		// Forward elimination steps
		diag_offset = 0;
		for (int i = 0; i < matrix_h - 1; i++, subm_h--, subm_w--, diag_offset += matrix_w + 1) {
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
			div_sizes[i] = (subm_h - 1) / proc_num_used[i];
			rest_sizes[i] = (subm_h - 1) - (proc_num_used[i] - 1) * div_sizes[i];

			// We will store there the mpi indexed types that we create
			// in order to transfer the matrix parts to the processes
			spec_types = (MPI_Datatype*)malloc(sizeof(MPI_Datatype) * (proc_num_used[i] - 1));

			int subm_offset = diag_offset + (rest_sizes[i] + 1) * matrix_w;
			for (int proc = 1; proc < proc_num_used[i]; proc++) {
				constr_row_and_subm_type(&spec_types[proc - 1], matrix_w, diag_offset, subm_offset, div_sizes[i], subm_w);
				// Send the matrix parts to the processes
				MPI_Send(matrix, 1, spec_types[proc - 1], proc, 0, MPI_COMM_WORLD);
				subm_offset += div_sizes[i] * matrix_w;
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
				MPI_Recv(matrix, 1, spec_types[proc - 1], proc, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				MPI_Type_free(&spec_types[proc - 1]);
			}

			free(spec_types);
		}
		divide_vector2value(&matrix[diag_offset], 2, matrix[diag_offset]);

		//printf("\nForward elimination's result:\n");
		//print_matrix(matrix, matrix_w, matrix_h);
		//printf("\n\n");

		// Back substitution
		pivot_elem_offset = matrix_h * matrix_w - 1;
		for (int i = 0; i < matrix_h - 1; i++, pivot_elem_offset -= matrix_w) {
			spec_types = (MPI_Datatype*)malloc(sizeof(MPI_Datatype) * (proc_num_used[i] - 1));

			int first_col_offset = matrix_w - i - 2;
			for (int proc = 1; proc < proc_num_used[i]; proc++) {
				constr_cols_and_elem_type(
					&spec_types[proc - 1],
					matrix_w,
					pivot_elem_offset,
					first_col_offset,
					first_col_offset + i + 1,
					div_sizes[i]
				);
				MPI_Send(matrix, 1, spec_types[proc - 1], proc, 0, MPI_COMM_WORLD);
				first_col_offset += div_sizes[i] * matrix_w;
			}

			int start_of_col = pivot_elem_offset - matrix_w * rest_sizes[i];
			for (int k = 0; k < rest_sizes[i]; k++, start_of_col += matrix_w) {
				matrix[start_of_col] -= matrix[pivot_elem_offset] * matrix[start_of_col - (i + 1)];
				//matrix[start_of_col - (i + 1)] = 0.0;
			}

			for (int proc = 1; proc < proc_num_used[i]; proc++) {
				MPI_Recv(matrix, 1, spec_types[proc - 1], proc, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				MPI_Type_free(&spec_types[proc - 1]);
			}

			free(spec_types);
		}

		printf("\nBack substitution's result:\n");
		print_matrix(matrix, matrix_w, matrix_h);
		printf("\n\n");

		if (argc == 3) {
			FILE* output_file = fopen(argv[2], "w");
			for (int i = 0; i < matrix_h; i++)
				fprintf(output_file, "%lf\n", matrix[(matrix_w - 1) + matrix_w * i]);
			fclose(output_file);
		}

		end_time = MPI_Wtime();
		double parallel_time = end_time - start_time;
		
		printf("Parallel alg: %f sec\n", parallel_time);
		printf("Sequential alg: %f sec\n", sequential_time);
		int successful_validation = 1;
		for (int last_elem = matrix_w - 1; last_elem < matrix_h * matrix_w; last_elem += matrix_w)
			if (fabs(matrix[last_elem] - matrix_copy[last_elem]) > EPSILON)
				successful_validation = 0;
		(successful_validation)
			? printf("SUCCESSFUL VALIDATION CHECK\n")
			: printf("UNSUCCESSFUL VALIDATION CHECK\n");
		free(matrix_copy);
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

		div_sizes = (int*)malloc(sizeof(int) * actions);
		rest_sizes = (int*)malloc(sizeof(int) * actions);
		proc_num_used = (int*)malloc(sizeof(int) * actions);

		int diag_offset;
		subm_h = matrix_h;
		subm_w = matrix_w;
		MPI_Datatype spec_type;
		// Forward elimination steps
		for (int i = 0; i < actions; i++, subm_w--, subm_h--) {
			proc_num_used[i] = (proc_num <= subm_h - 1) ? proc_num : subm_h - 1;
			div_sizes[i] = (subm_h - 1) / proc_num_used[i];
			rest_sizes[i] = (subm_h - 1) - (proc_num_used[i] - 1) * div_sizes[i];
			diag_offset = matrix_w * i + i;

			int subm_offset = diag_offset + matrix_w * ((rest_sizes[i] + 1) + (proc_rank - 1) * div_sizes[i]);
			constr_row_and_subm_type(&spec_type, matrix_w, diag_offset, subm_offset, div_sizes[i], subm_w);
			MPI_Recv(matrix, 1, spec_type, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			for (int k = 0; k < div_sizes[i]; k++) {
				double coef = matrix[subm_offset + k * matrix_w];
				if (fabs(coef) >= EPSILON)
					substract_vectors_with_coef(&matrix[subm_offset + k * matrix_w], &matrix[diag_offset], subm_w, coef);
			}

			MPI_Send(matrix, 1, spec_type, 0, 0, MPI_COMM_WORLD);
			MPI_Type_free(&spec_type);
		}
		// Back substitution
		int pivot_elem_offset = matrix_h * matrix_w - 1;
		for (int i = 0; i < actions; i++, pivot_elem_offset -= matrix_w) {
			int first_col_offset = (matrix_w - i - 2) + (proc_rank - 1) * matrix_w;
			int second_col_offset = first_col_offset + i + 1;
			constr_cols_and_elem_type(
				&spec_type,
				matrix_w,
				pivot_elem_offset,
				first_col_offset,
				second_col_offset,
				div_sizes[i]
			);
			MPI_Recv(matrix, 1, spec_type, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			for (int k = 0; k < rest_sizes[i]; k++, second_col_offset += matrix_w) {
				matrix[second_col_offset] -= matrix[pivot_elem_offset] * matrix[second_col_offset - (i + 1)];
				//matrix[start_of_col - (i + 1)] = 0.0;
			}

			MPI_Send(matrix, 1, spec_type, 0, 0, MPI_COMM_WORLD);
			MPI_Type_free(&spec_type);
		}
	}

	free(div_sizes);
	free(rest_sizes);
	free(proc_num_used);
	free(matrix);
	MPI_Finalize();
	return 0;
}
