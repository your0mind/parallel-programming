// Task 2: Calculating average value from vector

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include "mpi.h"

int main(int argc, char *argv[]) {
	int proc_num, proc_rank;
	int vector_size;
	int *vector;

	srand(time(NULL));

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	if (proc_rank == 0)  {
		double start_time = MPI_Wtime();

		// Creating random vector with certain size (argv[1])
		// with int values in range [argv[2], argv[3]]
		vector_size = atoi(argv[1]);
		MPI_Bcast(&vector_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
		vector = (int*)malloc(sizeof(int) * vector_size);

		printf("Source vector: ");
		for (int i = 0; i < vector_size; i++) {
			vector[i] = rand() % (atoi(argv[3]) + 1 - atoi(argv[2])) + atoi(argv[2]);
			printf("%d ", vector[i]);
		}
		printf("\n");

		// Sending parts of vector to other processes
		// The rest_size is needed to handle the case when vector_size % proc_num != 0
		int devision_size = (vector_size - 1) / proc_num + 1;
		int rest_size = vector_size - (proc_num - 1) * devision_size;

		for (int i = 1; i < proc_num; i++)
			MPI_Send(vector + rest_size + devision_size * (i - 1),
				devision_size, MPI_INT, i, 0, MPI_COMM_WORLD);

		int sum = 0;
		for (int i = 0; i < rest_size; i++)
			sum += vector[i];

		int final_sum;
		// Final sum of vector's elements
		MPI_Reduce(&sum, &final_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		printf("Average value: %f\n", final_sum / (double)vector_size);
		printf("Time: %f sec\n", MPI_Wtime() - start_time);
	}
	else {
		// Getting size of source vector
		MPI_Bcast(&vector_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// Calculating the size of the part of the vector that will be processed
		int size = (vector_size - 1) / proc_num + 1;
		
		MPI_Status status;
		vector = (int*)malloc(sizeof(int) * size);
		// Getting the part of the source vector
		MPI_Recv(vector, size, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		int sum = 0;
		for (int i = 0; i < size; i++)
			sum += vector[i];

		MPI_Reduce(&sum, NULL, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	free(vector);
	MPI_Finalize();
	return 0;
}
