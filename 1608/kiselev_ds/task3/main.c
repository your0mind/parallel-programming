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

int* bond_arrs_by_digit_place(int *arr1, int size1, int *arr2, int size2, int digit_number) {
	int size = size1 + size2;
	int *arr = (int*)malloc(sizeof(int) * size);
	for (int i = 0, i1 = 0, i2 = 0; i < size; i++) {
		if (get_digit_place(arr1[i1], digit_number) <= get_digit_place(arr2[i1], digit_number))
			arr[i] = arr1[i1++];
		else
			arr[i] = arr2[i2++];
	}
	return arr;
}

int main(int argc, char *argv[]) {
	int proc_num, proc_rank;
	int size;
	int *src_arr, *arr, *mate_arr;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	if (proc_rank == 0) {
		size = atoi(argv[1]);
		srand(time(NULL));
		src_arr = generate_arr(size, atoi(argv[2]), atoi(argv[3]));
	}

	MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int div_size = size / proc_num;
	int rest_size = size - proc_num * div_size;

	for (int digit_num = 1; digit_num <= NUMBER_OF_INT_DIGITS; digit_num++) {
		int cur_size;
		if (proc_rank == 0) {
			cur_size = rest_size;
			for (int proc = 1; proc < proc_num; proc++)
				MPI_Send(&src_arr[(proc - 1) * div_size], div_size, MPI_INT, proc, 0, MPI_COMM_WORLD);
			arr = (int*)malloc(sizeof(int) * cur_size);
			memcpy(&src_arr[size - rest_size], arr, cur_size);
			free(src_arr);
		}
		else {
			cur_size = div_size;
			arr = (int*)malloc(sizeof(int) * cur_size);
			MPI_Recv(arr, cur_size, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}
		qsort_by_digit_place(arr, 0, cur_size - 1, digit_num);

		for (int step = 1; step < proc_num; step <<= 1) {
			int mate_size = div_size * step;
			if ((proc_rank - step) % (step * 2) == 0) {
				MPI_Send(arr, mate_size, MPI_INT, proc_rank - step, 0, MPI_COMM_WORLD);
			}
			else if (proc_rank % (step * 2) == 0) {
				mate_arr = (int*)mallpc(sizeof(int) * mate_size);
				MPI_Recv(mate_arr, mate_size, MPI_INT, proc_rank + step, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				int *temp = arr;
				arr = bond_arrs_by_digit_place(arr, cur_size, mate_arr, mate_size, digit_num);
				free(temp);
			}
			else
				free(arr);
		}
	}

	MPI_Finalize();
	getchar();
	return 0;
}