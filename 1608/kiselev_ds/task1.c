// Lab 1: Calculating average value from vector's elements

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include "mpi.h"

int* generate_vector(int size, int lower_value, int upper_value) {
	int *vector = (int*)malloc(sizeof(int) * size);
	for (int i = 0; i < size; i++)
		vector[i] = rand() % (upper_value + 1 - lower_value) + lower_value;
	return vector;
}

int get_vector_sum(int *vector, int size) {
	int sum = 0;
	for (int i = 0; i < size; i++)
		sum += vector[i];
	return sum;
}

int main(int argc, char *argv[]) {
	int proc_num, proc_rank;
	int is_error = 0;
	int vector_size;
	int *vector;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	if (proc_rank == 0)  {
		// Input processing
		if ((argc != 4) || (atoi(argv[1]) < proc_num)) {
			if (argc != 4) {
				printf("Usage:\n\tmpi-lab1.exe <vector size> <lower value> <upper value>\n\n");
			} 
			else {
				printf("Vector size shouldn't be less than the number of processes.");
			}

			is_error = 1;
			MPI_Bcast(&is_error, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Finalize();
			return 1;
		} 
		else {
			MPI_Bcast(&is_error, 1, MPI_INT, 0, MPI_COMM_WORLD);
		}

		srand(time(NULL));
		double start_time = MPI_Wtime();

		// Creating random vector with certain size (argv[1])
		// with int values in range [argv[2], argv[3]]
		vector_size = atoi(argv[1]);
		vector = generate_vector(vector_size, atoi(argv[2]), atoi(argv[3]));

		// The rest_size is needed to handle the case when vector_size % proc_num != 0
		int devision_size = (vector_size - 1) / proc_num + 1;
		int rest_size = vector_size - (proc_num - 1) * devision_size;

		// Sending devision_size to all processes
		MPI_Bcast(&devision_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// Sending parts of vector to other processes
		for (int i = 1; i < proc_num; i++)
			MPI_Send(vector + rest_size + devision_size * (i - 1),
				devision_size, MPI_INT, i, 0, MPI_COMM_WORLD);

		int sum = get_vector_sum(vector, vector_size);

		int final_sum;
		// Final sum of vector's elements
		MPI_Reduce(&sum, &final_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		printf("Average value: %f\n", final_sum / (double)vector_size);
		printf("Time: %f sec\n", MPI_Wtime() - start_time);
	}
	else {
		// Getting information about input processing
		MPI_Bcast(&is_error, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (is_error == 1) {
			MPI_Finalize();
			return 1;
		}

		// Getting the size of part of vector that will be processed
		MPI_Bcast(&vector_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		MPI_Status status;
		vector = (int*)malloc(sizeof(int) * vector_size);
		// Getting the part of the source vector
		MPI_Recv(vector, vector_size, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		int sum = get_vector_sum(vector, vector_size);
		MPI_Reduce(&sum, NULL, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	free(vector);
	MPI_Finalize();
	return 0;
}
