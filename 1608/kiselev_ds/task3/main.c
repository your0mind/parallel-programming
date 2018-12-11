// Lab3: bitwise sorting int numbers using merge sorting
// The code works only for the number of processes
// equal to the powers of two (2, 4, 8, ...)

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

void print_arr(int *arr, int size) {
	for (int i = 0; i < size; i++)
		printf("%d ", arr[i]);
	printf("\n");
}

int* generate_arr(int size, int lower_value, int upper_value) {
	int *arr = (int*)malloc(sizeof(int) * size);
	for (int i = 0; i < size; i++)
		arr[i] = rand() % (upper_value + 1 - lower_value) + lower_value;
	return arr;
}

int digit_place(int number, int digit_number) {
	return number % (int)pow(10.0, digit_number) /
		(int)pow(10.0, digit_number - 1);
}

int compare_by_digit_place(int num1, int num2, int digit_num) {
	int res, digit1, digit2;
	do {
		int digit1 = digit_place(num1, digit_num);
		int digit2 = digit_place(num2, digit_num);
		res = digit1 - digit2;
	}  while ((res == 0) && (--digit_num > 0));
	return res;
}

void swap_int(int *elem1, int *elem2) {
	int tmp = *elem1;
	*elem1 = *elem2;
	*elem2 = tmp;
}

int min_arr(int *arr, int size) {
	int min = INT_MAX;
	for (int i = 0; i < size; i++)
		if (arr[i] < min)
			min = arr[i];
	return min;
}

int max_arr(int *arr, int size) {
	int max = INT_MIN;
	for (int i = 0; i < size; i++)
		if (arr[i] > max)
			max = arr[i];
	return max;
}

int digit_count(int number) {
	int count = 0;
	for (count; number > 0; count++, number /= 10)
		;
	return count;
}

void qsort_by_digit_place(int *v, int left, int right, int digit_num) {
	if (left >= right)
		return;
	swap_int(&v[left], &v[(left + right) / 2]);
	int last = left;
	for (int i = left + 1; i <= right; i++)
		if (compare_by_digit_place(v[i], v[left], digit_num) < 0)
			swap_int(&v[++last], &v[i]);
	swap_int(&v[left], &v[last]);
	qsort_by_digit_place(v, left, last - 1, digit_num);
	qsort_by_digit_place(v, last + 1, right, digit_num);
}

int* merge_by_digit_place(
	int *arr1,
	int size1,
	int *arr2,
	int size2,
	int digit_num
) {
	int size = size1 + size2;
	int *arr = (int*)malloc(sizeof(int) * size);
	int i = 0, i1 = 0, i2 = 0;
	while ((i1 < size1) && (i2 < size2)) {
		if (compare_by_digit_place(arr1[i1], arr2[i2], digit_num) < 0)
			arr[i++] = arr1[i1++];
		else
			arr[i++] = arr2[i2++];
	}
	if (i1 < size1)
		memcpy(&arr[i], &arr1[i1], (size - i) * sizeof(int));
	else if (i2 < size2)
		memcpy(&arr[i], &arr2[i2], (size - i) * sizeof(int));
	return arr;
}

int main(int argc, char *argv[]) {
	int proc_num, proc_rank;
	int max_digit_count;
	int size;
	int *src_arr, *arr = NULL, *mate_arr;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	if (proc_rank == 0) {
		srand(time(NULL));
		size = atoi(argv[1]);
		arr = generate_arr(size, atoi(argv[2]), atoi(argv[3]));

		// Counting the maximum number of digits to be processed
		int in_min = digit_count(min_arr(arr, size));
		int in_max = digit_count(max_arr(arr, size));
		max_digit_count = (in_min > in_max) ? in_min : in_max;

		printf("Source array:\n");
		print_arr(arr, size);
	}

	MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&max_digit_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// The rest_size is needed to handle the
	// case when size % proc_num != 0
	int div_size = size / proc_num;
	int rest_size = size - (proc_num - 1) * div_size;

	for (int digit_num = 1; digit_num <= max_digit_count; digit_num++) {
		int cur_size;
		if (proc_rank == 0) {
			src_arr = arr;
			cur_size = rest_size;
			for (int proc = 1; proc < proc_num; proc++)
				MPI_Send(
					&src_arr[(proc - 1) * div_size + cur_size],
					div_size,
					MPI_INT,
					proc,
					0,
					MPI_COMM_WORLD
				);
			arr = (int*)malloc(sizeof(int) * cur_size);
			memcpy(arr, src_arr, cur_size * sizeof(int));
			free(src_arr);
		}
		else {
			cur_size = div_size;
			arr = (int*)malloc(sizeof(int) * cur_size);
			MPI_Recv(
				arr,
				cur_size,
				MPI_INT,
				0,
				MPI_ANY_TAG,
				MPI_COMM_WORLD,
				&status
			);
		}
		qsort_by_digit_place(arr, 0, cur_size - 1, digit_num);

		for (int step = 1; step < proc_num; step <<= 1) {
			int mate_size = div_size * step;
			if ((proc_rank - step) % (step * 2) == 0) {
				// Sending own sorted part of the array to a mate
				MPI_Send(
					arr,
					mate_size,
					MPI_INT,
					proc_rank - step,
					0,
					MPI_COMM_WORLD
				);
				free(arr);
				break;
			}
			else if (proc_rank % (step * 2) == 0) {
				mate_arr = (int*)malloc(sizeof(int) * mate_size);
				// Receiving sorted part and merge it with own sorted part
				MPI_Recv(
					mate_arr,
					mate_size,
					MPI_INT,
					proc_rank + step,
					MPI_ANY_TAG,
					MPI_COMM_WORLD,
					&status
				);
				int *temp = arr;
				arr = merge_by_digit_place(
					arr,
					cur_size,
					mate_arr,
					mate_size,
					digit_num
				);
				cur_size += mate_size;
				free(temp);
			}
		}
	}

	if (proc_rank == 0) {
		printf("Sorted array:\n");
		print_arr(arr, size);
		free(arr);
	}

	MPI_Finalize();
	return 0;
}
