#include <vector>
#include <iostream>
#include "NelderMead.h"
#define N 4

double f(const std::vector<double>&x)
{
	/*double result = 0;
	for (int i = 0; i<N; i++)
		result += x[i] * x[i] - cos(18.*x[i] * x[i]);
	result += N;
	return result;*/
	return 24 * x[0] + 27 * x[1] + 138 * x[2] + 78 * x[3];
}

//double f(const std::vector<double>&x)
//{
//	/*
//	Функция многих переменных: функция Gaussian quartic.
//	Тестовая функция вещественной оптимизации.
//	Входные параметры:
//	x - указатель на исходный массив;
//	VHML_N - размер массива x.
//	Возвращаемое значение:
//	Значение тестовой функции в точке x.
//	*/
//	double VHML_Result = 0;
//	for (int i = 0; i < N; i++) VHML_Result += (i + 1)*x[i] * x[i] * x[i] * x[i];
//	/*VHML_Result += HML_RandomNormal(0, 1);*/
//	return VHML_Result;
//}

//double f(const std::vector<double>&x)
//{
//	/*
//	28
//	Функция многих переменных: Гипер-эллипсоид.
//	Тестовая функция вещественной оптимизации.
//	Входные параметры:
//	x - указатель на исходный массив;
//	VHML_N - размер массива x.
//	Возвращаемое значение:
//	Значение тестовой функции в точке x.
//	*/
//	double VHML_Result = 0;
//	for (int i = 0; i < N; i++)
//		VHML_Result += (i + 1)*(i + 1)*x[i] * x[i];
//	return VHML_Result;
//}

void nelder_mead()
{
	int maxIter = 50;
	double alpha = 1., beta = 0.5, gamma = 2.;
	
	std::vector<double>b(2);
	std::vector<double>g(2);
	std::vector<double>w(2);
	std::vector<double>mid(4);
	std::vector<double>xr(4);
	std::vector<double>xe(4);
	std::vector<double>xc(4);
	std::vector<double>c(4);

	b[0] = 0.; b[1] = 1.; g[0] = 1.; g[0] = 0.; w[0] = 1.; w[1] = 1.;

	for (int i = 0; i < maxIter; i++)
	{
		if (f(w) < f(b))
			std::swap(w, b);
		if (f(w) < f(g))
			std::swap(w, g);
		if (f(g) < f(b))
			std::swap(b, g);
		mid[0] = (g[0] + b[0]) / 2;
		mid[1] = (g[1] + b[1]) / 2;


		//reflection
		xr[0] = mid[0] + alpha * (mid[0] - w[0]);
		xr[1] = mid[1] + alpha * (mid[1] - w[1]);
		if (f(xr) < f(g))
		{
			w[0] = xr[0];
			w[1] = xr[1];
		}
		else
		{
			if (f(xr) < f(w))
				w = xr;
			c[0] = (w[0] + mid[0]) / 2;
			c[1] = (w[1] + mid[1]) / 2;
			if (f(c) < f(w))
				w = c;
		}
		//expansion
		if (f(xr) < f(b))
		{
			xe[0] = mid[0] + gamma*(xr[0] - mid[0]);
			if (f(xe) < f(xr))
				w = xe;
			else
				w = xr;
		}
		//contraction
		if (f(xr) < f(g))
		{
			xc[0] = mid[0] + beta * (w[0] - mid[0]);
			if (f(xc) < f(w))
				w = xc;
		}
	}
	std::cout << b[0] << " " << b[1] << std::endl;
	
}
int main()
{
	nelder_mead();
	system("pause");
	return 0;
}