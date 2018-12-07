#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

const int NUMBER_OF_INT_DIGITS = sizeof(int) * 8;

int* generate_arr(int size, int lower_value, int upper_value) {
	int *arr = (int*)malloc(sizeof(int) * size);
	for (int i = 0; i < size; i++)
		arr[i] = rand() % (upper_value + 1 - lower_value) + lower_value;
	return arr;
}

int get_digit_place(int number, int digit_number) {
	return number % (int)pow(10.0, digit_number) / (int)pow(10.0, digit_number - 1);
}

void swap_int(int *elem1, int *elem2) {
	int tmp = *elem1;
	*elem1 = *elem2;
	*elem2 = tmp;
}

void qsort_by_digit_place(int *v, int left, int right, int digit_number) {
	if (left >= right)
		return;
	swap_int(&v[left], &v[(left + right) / 2]);
	int last = left;
	for (int i = left + 1; i <= right; i++)
		if (get_digit_place(v[i], digit_number) <= get_digit_place(v[left], digit_number))
			swap_int(&v[++last], &v[i]);
	swap_int(&v[left], &v[last]);
	qsort_by_digit_place(v, left, last - 1, digit_number);
	qsort_by_digit_place(v, last + 1, right, digit_number);
}

int main(int argc, char *argv[]) {
	int proc_num, proc_rank;
	int size;
	int *arr;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	if (proc_rank == 0) {
		size = atoi(argv[1]);
		srand(time(NULL));
		arr = generate_arr(size, atoi(argv[2]), atoi(argv[3]));
	}

	MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int div_size = size / proc_num;
	int rest_size = size - proc_num * div_size;
	int stages = 0;
	for (int num = proc_num; num > 1; num >>= 1, stages++)
		;

	for (int digit_num = 1; digit_num <= NUMBER_OF_INT_DIGITS; digit_num++) {
		if (proc_rank == 0) {
			for (int proc = 1; proc < proc_num; proc++)
				MPI_Send(&arr[(proc - 1)* div_size], div_size, MPI_INT, proc, 0, MPI_COMM_WORLD);
			qsort_by_digit_place(&arr[size - rest_size], 0, rest_size - 1, digit_num);
		}
		else {
			arr = (int*)malloc(sizeof(int) * div_size);
			MPI_Status status;
			MPI_Recv(arr, div_size, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			qsort_by_digit_place(arr, 0, div_size - 1, digit_num);
		}

		for (int step = 1; step < proc_num; step <<= 1) {
			if ((proc_rank - step) % (step * 2) == 0) {
				// Then send
			}
			else if (proc_rank % (step * 2)) {
				// Then recv
			}
		}
	}

	MPI_Finalize();
	getchar();
	return 0;
}