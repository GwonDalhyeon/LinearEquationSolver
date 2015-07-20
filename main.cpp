#include "operators.h"

void main()
{
	int num = 4;
	double* A = new double[num*num];
	double* b = new double[num];
	double* x = new double[num];

	for (int i = 0; i < num; i++)
	{

	}

	b[0] = 16, b[1] = 17;
	b[2] = 15, b[3] = 9;
	csr sparseA = csr(num, num, A);
	x = CG(sparseA,b);
	//x = PCG(num,A,b);

	for (int i = 0; i < num; i++)
	{
		cout<<x[i]<<endl;
	}
	delete A, b, x;
}