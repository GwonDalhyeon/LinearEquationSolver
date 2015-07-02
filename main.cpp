#include "operators.h"

void main()
{
	int num = 4;
	double* A = new double[num*num];
	double* b = new double[num];
	double* x = new double[num];

	A[0*num + 0] = 4;
	A[0*num + 1] = 3;
	A[0*num + 2] = 2;
	A[0*num + 3] = 1;

	A[1*num + 0] = 3;
	A[1*num + 1] = 3;

	A[1*num + 2] = 2;
	A[1*num + 3] = 1;
	
	A[2*num + 0] = 2;
	A[2*num + 1] = 2;
	A[2*num + 2] = 2;
	A[2*num + 3] = 1;
	
	A[3*num + 0] = 1;
	A[3*num + 1] = 1;
	A[3*num + 2] = 1;
	A[3*num + 3] = 1;

	b[0] = 16, b[1] = 17;
	b[2] = 15, b[3] = 9;

	x = PCG(num,A,b);

	for (int i = 0; i < num; i++)
	{
		cout<<x[i]<<endl;
	}
	delete A, b, x;
}