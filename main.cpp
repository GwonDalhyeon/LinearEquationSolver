#include "operators.h"

void main()
{
	int num = 4;
	double** A = new double*[num];
	double* b = new double[num];
	double* x = new double[num];
	for (int i = 0; i < num; i++)
	{
		A[i] = new double[num];
	}
	A[0][0] = 4;
	A[0][1] = 3;
	A[0][2] = 2;
	A[0][3] = 1;

	A[1][0] = 3;
	A[1][1] = 3;

	A[1][2] = 2;
	A[1][3] = 1;
	
	A[2][0] = 2;
	A[2][1] = 2;
	A[2][2] = 2;
	A[2][3] = 1;
	
	A[3][0] = 1;
	A[3][1] = 1;
	A[3][2] = 1;
	A[3][3] = 1;

	b[0] = 16, b[1] = 17;
	b[2] = 15, b[3] = 9;

	x = PCG(num,A,b);

	for (int i = 0; i < num; i++)
	{
		delete A[i];
	}
	delete A, b, x;
}