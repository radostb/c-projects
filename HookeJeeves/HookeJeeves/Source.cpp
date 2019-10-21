#include <vector>
#include <iostream>
#include "HookeJeeves.h"
#define N 2

//double f(const std::vector<double>&x)
//{
//	double result = 0;
//	for (int i = 0; i<N; i++)
//		result += x[i] * x[i] - cos(18.*x[i] * x[i]);
//	result += N;
//	return result;
//}

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

double f(const std::vector<double>&x)
{
	/*
	28
	Функция многих переменных: Гипер-эллипсоид.
	Тестовая функция вещественной оптимизации.
	Входные параметры:
	x - указатель на исходный массив;
	VHML_N - размер массива x.
	Возвращаемое значение:
	Значение тестовой функции в точке x.
	*/
	double VHML_Result = 0;
	for (int i = 0; i < N; i++)
		VHML_Result += (i + 1)*(i + 1)*x[i] * x[i];
	return VHML_Result;
}

double search(const std::vector<double>&x, const std::vector<double>&delta, int d1, int k)
{
	std::vector<double>y(N);
	y = x;
	y[k] = y[k] + delta[k] * d1;
	if (f(y) < f(x))
		return y[k];
	else
	{
		y[k] = y[k] - 2* delta[k] * d1;
		if (f(y) < f(x))
			return y[k];
		else
		{
			y[k] = y[k] + delta[k] * d1;
			return y[k];
		}
	}
}

void hooke_jeeves()
{
	
	std::vector<double>x(N);
	std::vector<double>x0(N);
	std::vector<double>y(N);
	std::vector<double>delta(N);
	

	double eps = 0.1, lambda = 1.5, d1 = 1., alpha = 4.;
	int k = 0, i = 0, n=10;

	x0[0] = 0.5; x0[1] = 1.;
	y[0] = x0[0]; y[1] = x0[1];
	delta[0] = 0.2; delta[1] = 0.4;
	bool flag = true;
	while(flag)
	{
		y[0] = search(y, delta, d1, 0);
		y[1] = search(y, delta, d1, 1);
		if (f(y) < f(x0))
		{
			x = y;
			y[0] = x[0] + lambda*(x[0] - x0[0]);
			y[1] = x[1] + lambda*(x[1] - x0[1]);
			x0 = x;
		}
		else
		{
			if ((delta[0] < eps) || (delta[1] < eps))
			{
				std::cout << y[0] << " " << y[1] << std::endl;
				flag = false;
			}
			else
			{
				if (delta[0] > eps)
				{
					delta[0] /= alpha;
				}
				if (delta[1] > eps)
				{
					delta[1] /= alpha;
				}
				y = x0;
				x = x0;
			}
		}
	}
	
}

int main()
{
	hooke_jeeves();
	system("pause");
	return(0);

}