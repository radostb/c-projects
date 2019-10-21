#include "LM.h"
#include <vector>
#include <functional>
#include <iostream>
#define N 2

double z(double v)
{
	return (-(1. / ((v - 1.)*(v - 1.) + 0.2)) - (1. / (2. * (v - 2.)*(v - 2.) + 0.15)) - (1. / (3. * (v - 3.)*(v - 3.) + 0.3)));
}


double f1(double* x)
{
	/*return x[0]*x[0];*/
	return z(x[0]);
}

double f2(double* x)
{
	return z(x[1]);
	/*return x[1] * x[1];*/
}

//double f(double* x)
//{
//	return x[0] * x[0] + x[1] * x[1];
//	//return x[0] * x[0] - x[0] * x[1] + x[1] * x[1] + 9 * x[0] - 6 * x[1] + 20;
//}


double f(double *x)
{
	/*return 0.5*(1. - x[0])*(1. - x[0]) + 0.5*(x[1] - x[0] * x[0])*(x[1] - x[0] * x[0]);*/
	return -z(x[0])*z(x[1]);
	//return x[0] * x[0] - x[0] * x[1] + x[1] * x[1] + 9 * x[0] - 6 * x[1] + 20;
	/*return 2 * x1*x1 + x1*x2 + x2*x2;*/
	/*double VHML_Result;
	double z1 = -(1. / ((x - 1.)*(x - 1.) + 0.2)) - (1. / (2.*(x - 2.)*(x - 2.) + 0.15)) - (1. / (3.*(x - 3.)*(x
	- 3.) + 0.3));
	double z2 = -(1. / ((y - 1.)*(y - 1.) + 0.2)) - (1. / (2.*(y - 2.)*(y - 2.) + 0.15)) - (1. / (3.*(y - 3.)*(y
	- 3.) + 0.3));
	VHML_Result = z1 + z2;
	return VHML_Result;*/
	/*return (x1 - 5)*(x1 - 5)*(x2 - 4)*(x2 - 4) + (x1 - 5)*(x1 - 5) + (x2 - 4)*(x2 - 4) + 1;*/
}


void jacobi(double* x, double** jac)
{
	 /*Íàïèñàíî äëÿ äâóìåðíîãî ñëó÷àÿ*/
	double h = 1e-8;
	double y[N];
	y[0] = x[0] + h;
	y[1] = x[1];
	jac[0][0] = (f1(y) - f1(x)) / h;
	jac[1][0] = (f2(y) - f2(x)) / h;
	y[0] = x[0];
	y[1] = x[1] + h;
	jac[0][1] = (f1(y) - f1(x)) / h;
	jac[1][1] = (f2(y) - f2(x)) / h;
	y[1] = x[1];
	/*jac[0][0] = -1;
	jac[0][1] = 0;
	jac[1][0] = -2 * x[0];
	jac[1][1] = 1;*/
	
}

void transponateMatrix(double** matrix, double** matrixT)
{
	// Íàïèñàíî äëÿ äâóìåðíîãî ñëó÷àÿ
	matrixT[0][0] = matrix[0][0];
	matrixT[0][1] = matrix[1][0];
	matrixT[1][0] = matrix[0][1];
	matrixT[1][1] = matrix[1][1];
}

void multiplyMatrix(double** A, double** B, double** C)
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			C[i][j] = 0;
			for (int k = 0; k < N; k++)
				C[i][j] += A[i][k] * B[k][j];
		}
}

void multiplyVector(double** matrix, double* vector, double* solution)
{
	for (int i = 0; i < N; i++)
	{
		solution[i] = 0;
		for (int j = 0; j < N; j++)
			solution[i] += matrix[i][j] * vector[j];
	}
}

void solve(double** A, double* b, double* x)
{
	//double max = 0, buf = 0;
	//int i = 0;
	//for (int k = 0; k < N; k++) //Поиск максимального элемента в первом столбце
	//{
	//	max = A[k][k];
	//	i = k;
	//	for (int m = k + 1; m < N; m++)
	//		if (abs(A[m][k]) > max)
	//		{
	//			i = m;
	//			max = abs(A[m][k]);
	//		}
	//
	//	if (i != k) // перестановка i-ой строки, содержащей главный элемент k-ой строки
	//	{
	//		for (int j = k; j < N; j++)
	//		{
	//			buf = A[k][j];
	//			A[k][j] = A[i][j];
	//			A[i][j] = buf;
	//		}
	//	}

	//	max = A[k][k];

	//	//преобразование k-ой строки (Вычисление масштабирующих множителей)

	//	for (int j = k; j < N; j++)
	//		A[k][j] = A[k][j] / max;


	//	for (int i = k + 1; i < N; i++)//преобразование строк с помощью k-ой строки
	//	{
	//		buf = A[i][k];
	//		for (int j = k; j < N; j++)
	//		{
	//			A[i][j] = A[i][j] - buf * A[k][j];
	//			std::cout << A[i][j] << " ";
	//		}
	//		b[i] = b[i] - buf * b[k];
	//		std::cout << b[i] << std::endl;
	//	}

	//	//for (int i = 0; i < N; i++)//преобразование строк с помощью k-ой строки
	//	//{
	//	//	for (int j = 0; j < N; j++)
	//	//	{
	//	//		std::cout << A[i][j] << " ";
	//	//	}
	//	//	std::cout << std::endl;
	//	//}
	//	//for (int i = 0; i < N; i++)
	//	//	std::cout << b[i] << " ";

	//}

	//for (int i = N - 1; i >= 0; i--) //Нахождение решений СЛАУ

	//{

	//	x[i] = b[i] / A[i][i];

	//	for (int j = i + 1; j < N; j++)
	//		x[i] = x[i] - A[i][j] * x[j] / A[i][i];
	//}

	int   i, j, k, m, rowx;
	double xfac, temp, temp1, amax;

	rowx = 0;
	for (k = 0; k < N - 1; ++k) {
		amax = (double)fabs(A[k][k]);
		m = k;
		for (i = k + 1; i < N; i++) {
			xfac = (double)fabs(A[i][k]);
			if (xfac > amax) { amax = xfac; m = i; }
		}
		if (m != k) {
			rowx = rowx + 1;
			temp1 = b[k];
			b[k] = b[m];
			b[m] = temp1;
			for (j = k; j < N; j++) {
				temp = A[k][j];
				A[k][j] = A[m][j];
				A[m][j] = temp;
			}
		}
		for (i = k + 1; i < N; ++i) {
			xfac = A[i][k] / A[k][k];

			for (j = k + 1; j < N; ++j) {
				A[i][j] = A[i][j] - xfac * A[k][j];
			}
			b[i] = b[i] - xfac * b[k];
		}

	}

	for (j = 0; j < N; ++j) {
		k = N - j - 1;
		x[k] = b[k];
		for (i = k + 1; i < N; ++i) {
			x[k] = x[k] - A[k][i] * x[i];
		}
		x[k] = x[k] / A[k][k];
	}
}
double determinant(double** matrix)
{
	// Íàïèñàíî äëÿ äâóìåðíîãî ñëó÷àÿ
	return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
}

void inverseMatrix(double** matrix, double** inverse)
{
	// Íàïèñàíî äëÿ äâóìåðíîãî ñëó÷àÿ
	double **temp = new double *[N];
	for (int i = 0; i < N; i++) {
		temp[i] = new double[N];
	}
	double det = determinant(matrix);
	int n = 2;
	temp[0][0] = pow(-1, n++)*matrix[1][1] / det;
	temp[1][0] = pow(-1, n++)*matrix[0][1] / det;
	temp[0][1] = pow(-1, n++)*matrix[1][0] / det;
	temp[1][1] = pow(-1, n++)*matrix[0][0] / det;
	transponateMatrix(temp, inverse);

	for (int i = 0; i < N; i++) {
		delete[]temp[i];
	}
	delete[]temp;
}

//void findJacobi(double** jacobi, std::vector<double> &x, std::vector<std::function<double(double*)>>& f, std::vector<double> &unit_vector, double step)
//{
//	//ñäåëàòü íîðìàëüíî
//	std::vector<double> temp_vector(N);
//	//std::vector<double>::iterator it;
//	temp_vector = x;
//	for (int i = 0; i < N; i++)
//		for (int j = 0; j < N; j++)
//		{
//			temp_vector[j] += step;
//			//unit_vector+j += step;
//			jacobi[i][j] = (f[i](temp_vector) - f[i](x))/step;
//			
//		}
//	
//}

int main()
{
	double mu = 0, k = 0;
	double eps1 = 1e-8;
	double eps2 = 1e-8;
	int kmax = 1000;
	double gainRatio = 0;
	double max = 0;
	int nu = 3;
	int found = 0;

	double *x = new double[N];
	double *xnew = new double[N];
	double *g = new double[N];
	double *r = new double[N];
	double *d = new double[N];
	double *l = new double[N];
	
	double L = 0, L0 = 0;
	double **jac = new double *[N];
	double **jacT = new double *[N];
	double **A = new double *[N];
	double **H = new double *[N];
	double **invA = new double *[N];

	for (int i = 0; i < N; i++)
	{
		jac[i] = new double[N];
		jacT[i] = new double[N];
		A[i] = new double[N];
		H[i] = new double[N];
		invA[i] = new double[N];
		g[i] = 0;
		x[i] = 0;
		r[i] = 0;
		l[i] = 0;
		d[i] = 0;
		xnew[i] = 0;
		for (int j = 0; j < N; j++)
		{
			A[i][j] = 0;
			H[i][j] = 0;
			invA[i][j] = 0;
			jac[i][j] = 0;
			jacT[i][j] = 0;
		}
	}
	
	x[0] = 1.7;
	x[1] = 1.7;
	
	jacobi(x, jac);
	transponateMatrix(jac, jacT);
	multiplyMatrix(jac, jacT, A);
	r[0] = f1(x);
	r[1] = f2(x);
	multiplyVector(jacT, r, g);
	for (int i = 0; i < N; i++)
		if (A[i][i] > max)
			max = A[i][i];
	mu = 0.005*max;//fudge factor
	if (g[0] * g[0] + g[1] * g[1] <= eps1)
		found = 1;
	while ((found == 0) && (k < kmax))
	{
		k++;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
				H[i][j] = A[i][j];
			H[i][i] = A[i][i] + mu * A[i][i];
		}
		solve(H, g, d);
		if ((d[0] * d[0] + d[1] * d[1]) < eps2*(x[0] * x[0] + x[1] * x[1] + eps2))
			found = 1;
		else
		{
			xnew[0] = x[0] + d[0];
			xnew[1] = x[1] + d[1];
			r[0] = f1(xnew) - f1(x);
			r[1] = f2(xnew) - f2(x);
			multiplyVector(jac, d, l);
			l[0] += r[0];
			l[1] += r[1];
			for (int i = 0; i < N; i++)
			{
				L += l[i] * l[i];
				L0 += r[i] * r[i];
			}
			L = 0.5*L;
			L0 = 0.5*L0;
			gainRatio = (f(x) - f(xnew)) / (L0 - L) ;
			L = 0;
			if (gainRatio > 0)
			{
				x[0] = xnew[0];
				x[1] = xnew[1];
				std::cout << x[0] << " " << x[1] << " " << f(x) << std::endl;
				jacobi(x, jac);
				transponateMatrix(jac, jacT);

				multiplyMatrix(jac, jacT, A);
				r[0] = f1(x);
				r[1] = f2(x);

				multiplyVector(jacT, r, g);
				if ((g[0] * g[0] + g[1] * g[1]) <= eps1)
					found = 1;
				if (0.33 > (1. - pow(2. * gainRatio - 1., 3.)))
					mu = mu * 0.33;
				else mu = mu * (1. - pow(2. * gainRatio - 1., 3.));
				nu = 3.;
			}
			else
			{
				mu = mu * nu;
				nu = 2. * nu;
			}
		}
	}
	std::cout << k << std::endl;
	/*double *x = new double[N];
	double *xnew = new double[N];
	double *g = new double[N];
	double *r = new double[N];
	double *d = new double[N];
	double *l = new double[N];

	double L = 0, L0 = 0;
	double **jac = new double *[N];
	double **jacT = new double *[N];
	double **A = new double *[N];
	double **invA = new double *[N];*/

	delete[] x;
	delete[] xnew;
	delete[] g;
	delete[] r;
	delete[] d;
	delete[] l;

	for (int i = 0; i < N; i++) {
		delete[]jac[i];
		delete[]jacT[i];
		delete[]A[i];
		delete[]invA[i];
	}
	delete[]jac;
	delete[]jacT;
	delete[]A;
	delete[]invA;

	system("pause");
}