#include <stdio.h>
#include <math.h>
#include "mpi.h"

const int NUMBER_OF_DIGITS = sizeof(int) * 8;

int get_digit_place(int number, int digit_number) {
	return number % (int)pow(10.0, digit_number) / (int)pow(10.0, digit_number - 1);
}

void swap_int(int *elem1, int *elem2) {
	int *tmp = elem1;
	elem1 = elem2;
	elem2 = tmp;
}

void qsort_by_digit_place(int *v, int left, int right, int digit_number) {
	if (left >= right)
		return;
	swap_int(&v[left], &v[(left + right) / 2]);
	int last = left;
	for (int i = left + 1; i <= right; i++)
		//if (get_digit_place(v[i], digit_number) < get_digit_place(v[left], digit_number))
		if (v[i] - v[left] < 0) {
			printf("%d %d\n", v[i], v[left]);
			printf("%d %d\n", v[last + 1], v[i]);
			swap_int(&v[++last], &v[i]);
			printf("%d %d\n", v[last], v[i]);
		}
	swap_int(&v[left], &v[last]);
	qsort_by_digit_place(v, left, last - 1, digit_number);
	qsort_by_digit_place(v, last + 1, right, digit_number);
}

int main(int argc, char *argv[]) {
	int proc_num, proc_rank;
	int* num = (int*)malloc(sizeof(int) * 10);
	for (int i = 0; i < 10; i++) num[i] = 10 - i;
	//qsort_by_digit_place(num, 0, 9, 1);
	swap_int(&num[0], &num[9]);
	for (int i = 0; i <10; i++) printf("%d ", num[i]);
	/*MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	MPI_Finalize();*/
	getchar();
	return 0;
}