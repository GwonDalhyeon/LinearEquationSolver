#include "csr.h"

using namespace std;


double* CG(csr& A, double* b)
{
	int num = A.colNum;
	double tolerance = 1000*DBL_EPSILON;
	double* rOld = new double[num];
	double* p = new double[num];
	double* rNew = new double[num];
	double* x = new double[num];
	int j;

	for (int i = 0; i < num; i++)
	{
		rOld[i] = b[i];
		p[i] = b[i];
		x[i] = 0;
	}

	double alpha = 0;
	double beta = 0;
	double temp1 = 0, temp2 = 0;
	double temp = 0;
	double residual;
	double residualOld = 0;
	for (int i = 0; i < num; i++)
	{
		residualOld = residualOld + rOld[i]*rOld[i];
	}


	
	for (int k = 0; k < 2*num; k++)
	{
		temp1 = 0;
		temp2 = 0;
		for (int i = 0; i < num; i++)
		{
			for (int n = A.indptr[i]; n < A.indptr[i+1]; n++)
			{
				j = A.col[n];
				temp2 = temp2 + p[i]*A.val[n]*p[j];
			}
		}
		alpha = residualOld/temp2;
	
		for (int i = 0; i < num; i++)
		{
			x[i] = x[i] +alpha*p[i];
			temp = 0;
			for (int n = A.indptr[i]; n < A.indptr[i+1]; n++)
			{
				j = A.col[n];
				temp = temp + A.val[n]*p[j];
			}
			rNew[i] = rOld[i] - alpha*temp;
		}

		residual = 0;
		for (int i = 0; i < num; i++)
		{
			residual = residual + rNew[i]*rNew[i];
			//cout<<rNew[i]<<endl;
		}

		temp = sqrt(abs(residual));
		if (k%10==0)
		{
			cout<<".";
		}
		cout<<k<<" "<<temp<<endl;
		if (temp < tolerance)
		{
			delete rNew, rOld, p;
			cout<<"CG iterataion : "<<k<<endl;
			cout<< endl;
			return x;
		}

		beta = residual/residualOld;

		for (int i = 0; i < num; i++)
		{
			p[i] = rNew[i] + beta*p[i];
			rOld[i] = rNew[i];
		}
		residualOld = residual;
	
	}
	delete rNew, rOld, p;
	return x;
	
}

double* CG(int num, double** A, double* b)
{
	double tolerance = 10*DBL_EPSILON;
	double* rOld = new double[num];
	double* p = new double[num];
	double* rNew = new double[num];
	double* x = new double[num];

	for (int i = 0; i < num; i++)
	{
		rOld[i] = b[i];
		p[i] = b[i];
		x[i] = 0;
	}

	double alpha = 0;
	double beta = 0;
	double temp1 = 0, temp2 = 0;
	double temp = 0;
	double residual;

	for (int k = 0; k < 1000; k++)
	{
		temp1 = 0;
		temp2 = 0;
		for (int i = 0; i < num; i++)
		{
			temp1 = temp1 + rOld[i]*rOld[i];
			for (int j = 0; j < num; j++)
			{
				temp2 = temp2 + p[i]*A[i][j]*p[j];
			}
		}
		alpha = temp1/temp2;
	
		for (int i = 0; i < num; i++)
		{
			x[i] = x[i] +alpha*p[i];
			temp = 0;
			for (int j = 0; j < num; j++)
			{
				temp = temp + A[i][j]*p[j];
			}
			rNew[i] = rOld[i] - alpha*temp;
		}

		residual = 0;
		for (int i = 0; i < num; i++)
		{
			residual = residual + rNew[i]*rNew[i];
		}
		residual = sqrt(residual);

		if (residual < tolerance)
		{
			delete rNew, rOld, p;
			return x;
		}

		temp1 = 0;
		temp2 = 0;
		for (int i = 0; i < num; i++)
		{
			temp1 = temp1 + rNew[i]*rNew[i];
			temp2 = temp2 + rOld[i]*rOld[i];
		}
		beta = temp1/temp2;

		for (int i = 0; i < num; i++)
		{
			p[i] = rNew[i] + beta*p[i];
		}
	
		for (int i = 0; i < num; i++)
		{
			rOld[i] = rNew[i];
		}

	}
	delete rNew, rOld, p;
	return x;
	
}


double* CG(int num, double* A, double* b)
{
	double tolerance = 1000*DBL_EPSILON;
	double* rOld = new double[num];
	double* p = new double[num];
	double* rNew = new double[num];
	double* x = new double[num];

	for (int i = 0; i < num; i++)
	{
		rOld[i] = b[i];
		p[i] = b[i];
		x[i] = 0;
	}

	double alpha = 0;
	double beta = 0;
	double temp1 = 0, temp2 = 0;
	double temp = 0;
	double residual;
	double residualOld = 0;
	for (int i = 0; i < num; i++)
	{
		residualOld = residualOld + rOld[i]*rOld[i];
	}

	for (int k = 0; k < 2*num; k++)
	{
		temp1 = 0;
		temp2 = 0;
		for (int i = 0; i < num; i++)
		{
			for (int j = 0; j < num; j++)
			{
				temp2 = temp2 + p[i]*A[i*num + j]*p[j];
			}
		}
		alpha = residualOld/temp2;
	
		for (int i = 0; i < num; i++)
		{
			x[i] = x[i] +alpha*p[i];
			temp = 0;
			for (int j = 0; j < num; j++)
			{
				temp = temp + A[i*num + j]*p[j];
			}
			rNew[i] = rOld[i] - alpha*temp;
		}

		residual = 0;
		for (int i = 0; i < num; i++)
		{
			residual = residual + rNew[i]*rNew[i];
			//cout<<rNew[i]<<endl;
		}

		temp = sqrt(abs(residual));

		cout<<k<<" "<<temp<<endl;
		if (temp < tolerance)
		{
			delete rNew, rOld, p;
			cout<<"CG iterataion : "<<k<<endl;
			cout<< endl;
			return x;
		}

		beta = residual/residualOld;

		for (int i = 0; i < num; i++)
		{
			p[i] = rNew[i] + beta*p[i];
			rOld[i] = rNew[i];
		}
		residualOld = residual;
	
	}
	delete rNew, rOld, p;
	return x;
	
}
double** incompleteCholesky(int num, double** A)
{
	double** L;
	L = new double*[num];
	for (int i = 0; i < num; i++)
	{
		L[i] = new double[num];
		for (int j = 0; j < num; j++)
		{
			L[i][j] = 0;
		}
	}

	double temp;
	for (int i = 0; i < num; i++)
	{
		temp = 0;
		if (A[i][i]!=0)
		{
			for (int j = 0; j <= i-1; j++)
			{
				temp = temp + L[i][j]*L[i][j];
			}
			L[i][i] = sqrt(A[i][i] - temp);
		}
		
		for (int j = i+1; j < num; j++)
		{
			if (A[i][j]!=0 && L[i][i]!=0)
			{
				temp = 0;
				for (int k = 0; k <= i-1; k++)
				{
					temp = temp + L[i][k]*L[j][k];
				}
				L[j][i] = (A[i][j] - temp)/L[i][i];
			}
		}
	}

	return L;
}

double* incompleteCholesky(int num, double* A)
{
	double* L= new double[num*num];
	for (int i = 0; i < num; i++)
	{
		for (int j = 0; j < num; j++)
		{
			L[i*num + j] = 0;
		}
	}

	double temp;
	for (int i = 0; i < num; i++)
	{
		temp = 0;
		if (A[i*num + i]!=0)
		{
			for (int j = 0; j <= i-1; j++)
			{
				temp = temp + L[i*num + j]*L[i*num + j];
			}
			L[i*num + i] = sqrt(A[i*num + i] - temp);
		}
		
		for (int j = i+1; j < num; j++)
		{
			if (A[i*num + j]!=0 && L[i*num + i]!=0)
			{
				temp = 0;
				for (int k = 0; k <= i-1; k++)
				{
					temp = temp + L[i*num + k]*L[j*num + k];
				}
				L[j*num + i] = (A[i*num + j] - temp)/L[i*num + i];
			}
		}
	}

	return L;
}



double* PCG(int num, double** A, double* b)
{
	double tolerance = 10*DBL_EPSILON;
	double* rOld = new double[num];
	double* p = new double[num];
	double* rNew = new double[num];
	double* zOld = new double[num];
	double* zNew = new double[num];
	double* x = new double[num];

	double** L = incompleteCholesky(num,A);
	double** M = new double*[num];
	for (int i = 0; i < num; i++)
	{
		M[i] = new double[num];
		for (int j = 0; j < num; j++)
		{
			M[i][j]  = 0;
			for (int k = 0; k < num; k++)
			{
				M[i][j] = M[i][j] + L[i][k]*L[j][k];
			}
		}
	}

	double temp = 0;
	for (int i = 0; i < num; i++)
	{
		temp = 0;
		for (int j = 0; j < num; j++)
		{
			temp = temp + M[j][i]*rOld[j];
		}
		zOld[i] = temp;
		p[i] = zOld[i];
		rOld[i] = b[i];
		x[i] = 0;

	}

	double alpha = 0;
	double beta = 0;
	double temp1 = 0, temp2 = 0;
	double residual;

	for (int k = 0; k < 1000; k++)
	{
		temp1 = 0;
		temp2 = 0;
		for (int i = 0; i < num; i++)
		{
			temp1 = temp1 + rOld[i]*zOld[i];
			for (int j = 0; j < num; j++)
			{
				temp2 = temp2 + p[i]*A[i][j]*p[j];
			}
		}
		alpha = temp1/temp2;
	
		for (int i = 0; i < num; i++)
		{
			x[i] = x[i] +alpha*p[i];
			temp = 0;
			for (int j = 0; j < num; j++)
			{
				temp = temp + A[i][j]*p[j];
			}
			rNew[i] = rOld[i] - alpha*temp;
		}

		residual = 0;
		for (int i = 0; i < num; i++)
		{
			residual = residual + rNew[i]*rNew[i];
		}
		residual = sqrt(residual);

		if (residual < tolerance)
		{
			delete zNew, zOld, rNew, rOld, p;
			for (int i = 0; i < num; i++)
			{
				delete M[i], L[i];
			}
			delete M,L;

			return x;
		}

		for (int i = 0; i < num; i++)
		{
			temp = 0;
			for (int j = 0; j < num; j++)
			{
				temp = temp + M[j][i]*rNew[j];
			}
			zNew[i] = temp;
		}

		temp1 = 0;
		temp2 = 0;
		for (int i = 0; i < num; i++)
		{
			temp1 = temp1 + rNew[i]*zNew[i];
			temp2 = temp2 + rOld[i]*zOld[i];
		}
		beta = temp1/temp2;

		for (int i = 0; i < num; i++)
		{
			p[i] = zNew[i] + beta*p[i];
		}
	
		for (int i = 0; i < num; i++)
		{
			rOld[i] = rNew[i];
			zOld[i] = zNew[i];
		}

	}

	delete zNew, zOld, rNew, rOld, p;
	for (int i = 0; i < num; i++)
	{
		delete M[i], L[i];
	}
	delete M,L;

	return x;
	
}

double* PCG(int num, double* A, double* b)
{
	double tolerance = 10*DBL_EPSILON;
	double* rOld = new double[num];
	double* p = new double[num];
	double* rNew = new double[num];
	double* zOld = new double[num];
	double* zNew = new double[num];
	double* x = new double[num];

	double* L = incompleteCholesky(num,A);
	double* M = new double[num*num];
	for (int i = 0; i < num; i++)
	{
		for (int j = 0; j < num; j++)
		{
			M[i*num + j]  = 0;
			for (int k = 0; k < num; k++)
			{
				M[i*num + j] = M[i*num + j] + L[i*num + k]*L[j*num + k];
			}
		}
	}

	double temp = 0;
	for (int i = 0; i < num; i++)
	{
		temp = 0;
		for (int j = 0; j < num; j++)
		{
			temp = temp + M[j*num + i]*rOld[j];
		}
		zOld[i] = temp;
		p[i] = zOld[i];
		rOld[i] = b[i];
		x[i] = 0;

	}

	double alpha = 0;
	double beta = 0;
	double temp1 = 0, temp2 = 0;
	double residual;

	for (int k = 0; k < 1000; k++)
	{
		temp1 = 0;
		temp2 = 0;
		for (int i = 0; i < num; i++)
		{
			temp1 = temp1 + rOld[i]*zOld[i];
			for (int j = 0; j < num; j++)
			{
				temp2 = temp2 + p[i]*A[i*num + j]*p[j];
			}
		}
		alpha = temp1/temp2;
	
		for (int i = 0; i < num; i++)
		{
			x[i] = x[i] +alpha*p[i];
			temp = 0;
			for (int j = 0; j < num; j++)
			{
				temp = temp + A[i*num + j]*p[j];
			}
			rNew[i] = rOld[i] - alpha*temp;
		}

		residual = 0;
		for (int i = 0; i < num; i++)
		{
			residual = residual + rNew[i]*rNew[i];
		}
		residual = sqrt(residual);

		if (residual < tolerance)
		{
			delete zNew, zOld, rNew, rOld, p;
			delete M,L;

			return x;
		}

		for (int i = 0; i < num; i++)
		{
			temp = 0;
			for (int j = 0; j < num; j++)
			{
				temp = temp + M[j*num + i]*rNew[j];
			}
			zNew[i] = temp;
		}

		temp1 = 0;
		temp2 = 0;
		for (int i = 0; i < num; i++)
		{
			temp1 = temp1 + rNew[i]*zNew[i];
			temp2 = temp2 + rOld[i]*zOld[i];
		}
		beta = temp1/temp2;

		for (int i = 0; i < num; i++)
		{
			p[i] = zNew[i] + beta*p[i];
		}
	
		for (int i = 0; i < num; i++)
		{
			rOld[i] = rNew[i];
			zOld[i] = zNew[i];
		}

	}

	delete zNew, zOld, rNew, rOld, p;
	delete M,L;

	return x;
	
}

