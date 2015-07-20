#include <iostream>
#include <cmath>
#include <fstream>
//#include <math.h>

using namespace std;

class csr
{
public:
	double* val;
	double* col;
	int* indptr;
	int colNum;
	int rowNum;

	csr(int num1,int num2, double* matrix);

	csr();
	~csr();

private:

};
csr::csr()
{

}
csr::csr(int num1,int num2, double* matrix)
{
	rowNum = num1;
	colNum = num2;
	double* tempVal = new double[int(floor(sqrt(double(rowNum*colNum))))*10];
	double* tempCol = new double[int(floor(sqrt(double(rowNum*colNum))))*10];
	indptr = new int[rowNum+1];

	for (int i = 0; i < rowNum+1; i++)
	{
		indptr[i] = -1;

	}
	int tempIndex = 0;


	for (int i = 0; i < rowNum; i++)
	{
		for (int j = 0; j < colNum; j++)
		{
			if (matrix[i*colNum + j] != 0)
			{
				tempVal[tempIndex] = matrix[i*colNum + j];
				tempCol[tempIndex] = j;
				if (indptr[i]<0)
				{
					indptr[i] = tempIndex;
				}
				tempIndex = tempIndex + 1;
			}
		}
	}
	indptr[rowNum] = tempIndex;

	val = new double[tempIndex];
	col = new double[tempIndex];

	for (int i = 0; i < tempIndex; i++)
	{
		val[i] = tempVal[i];
		col[i] = tempCol[i];
	}

	delete[] tempVal, tempCol;
}

csr::~csr()
{
	delete[] val, col, indptr;
}