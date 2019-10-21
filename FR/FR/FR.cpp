#include <iostream>
#include <cmath>
#include <vector>
//#include "FR.h"
#define N 2
double h = 1e-8;

double z(double v)
{
	return (-(1. / ((v - 1.)*(v - 1.) + 0.2)) - (1. / (2. * (v - 2.)*(v - 2.) + 0.15)) - (1. / (3. * (v - 3.)*(v - 3.) + 0.3)));
}

double f(double *x)
{
	return -z(x[0])*z(x[1]);
	/*return 2 * x[0]*x[0] + x[0]*x[1] + x[1]*x[1];*/
	/*return x[0] * x[0] * x[0] - x[0] * x[1] + x[1] * x[1] - 2 * x[0] + 3 * x[1] - 4;*/
	/*return (x[1] - x[0] * x[0])*(x[1] - x[0] * x[0]) + (1 - x[0])*(1 - x[0]);*/
	/*double VHML_Result;
	double z1 = -(1. / ((x - 1.)*(x - 1.) + 0.2)) - (1. / (2.*(x - 2.)*(x - 2.) + 0.15)) - (1. / (3.*(x - 3.)*(x
		- 3.) + 0.3));
	double z2 = -(1. / ((y - 1.)*(y - 1.) + 0.2)) - (1. / (2.*(y - 2.)*(y - 2.) + 0.15)) - (1. / (3.*(y - 3.)*(y
		- 3.) + 0.3));
	VHML_Result = z1 + z2;
	return VHML_Result;*/
	/*return (x1 - 5)*(x1 - 5)*(x2 - 4)*(x2 - 4) + (x1 - 5)*(x1 - 5) + (x2 - 4)*(x2 - 4) + 1;*/
}

double dfdx1(double *x)
{
	double y[2];
	y[0] = x[0] + h;
	y[1] = x[1];
	return (f(y) - f(x)) / h;
	//return (6.*(-3. + x[0]) / ((0.3 + 3.*(-3. + x[0])*(-3. + x[0]))*(0.3 + 3.*(-3. + x[0])*(-3. + x[0])))) + (4.*(-2. + x[0]) / ((0.15 + 2.*(-2. + x[0])*(-2. + x[0]))*(0.15 + 4.*(-2. + x[0])*(-2. + x[0])))) + (2.*(-1. + x[0]) / ((0.2 + (-1. + x[0])*(-1. + x[0]))*(0.2 + (-1. + x[0])*(-1. + x[0]))));
}

double dfdx2(double *x)
{
	double y[2];
	y[0] = x[0];
	y[1] = x[1] + h;
	return (f(y) - f(x)) / h;
	//return (6.*(-3. + x[1]) / ((0.3 + 3.*(-3. + x[1])*(-3. + x[1]))*(0.3 + 3.*(-3. + x[1])*(-3. + x[1])))) + (4.*(-2. + x[1]) / ((0.15 + 2.*(-2. + x[1])*(-2. + x[1]))*(0.15 + 4.*(-2. + x[1])*(-2. + x[1])))) + (2.*(-1. + x[1]) / ((0.2 + (-1. + x[1])*(-1. + x[1]))*(0.2 + (-1. + x[1])*(-1. + x[1]))));
}

double dfdt(double *x)
{
	double t[2];
	t[0] = x[0] + h;
	t[1] = x[1] + h;
	return (f(t) - f(x)) / h;
}

double dihotomia(double *x, double *d)
{
	double a = -10, b = 10;
	double eps = 0.02, delta = 0.01;
	double t[2]; t[0] = 0.0; t[1] = 0.0;
	double x1 = (a + b - delta) / 2;
	double x2 = (a + b + delta) / 2;
	double y1, y2;
	while ((b - a) > eps)
	{
		t[0] = x[0] + d[0] * x1;
		t[1] = x[1] + d[1] * x1;
		y1 = f(t);
		t[0] = x[0] + d[0] * x2;
		t[1] = x[1] + d[1] * x2;
		y2 = f(t);
		if (y1 > y2)
			a = x1;
		else
			b = x2;
		x1 = (a + b - delta) / 2;
		x2 = (a + b + delta) / 2;
		
	}
	return (a + b) / 2;
}

double scalarProduct(double* vector1, double* vector2)
{
	double scalar = 0.;
	for (int i = 0; i < N; i++)
		scalar += vector1[i] * vector2[i];
	return scalar;
}

double wolfe(double *x, double *p, double max_iter)
{
	double alpha = 0, beta = 1000, step = 0.1, c1 = 0.0001, c2 = 0.9;
	double leftf, rightf;
	int i = 0;

	double *valf = new double[N];
	double *grad = new double[N];
	double *gradf = new double[N];

	grad[0] = dfdx1(x);
	grad[1] = dfdx2(x);

	while (i <= max_iter)
	{

		valf[0] = x[0] + step * c1*scalarProduct(grad, p);
		valf[1] = x[1] + step * c1*scalarProduct(grad, p);
		rightf = f(valf);
		valf[0] = x[0] + step * p[0];
		valf[1] = x[1] + step * p[1];
		leftf = f(valf);
		gradf[0] = dfdx1(valf);
		gradf[1] = dfdx2(valf);

		if (leftf > rightf)
		{
			beta = step;
			step = .5*(alpha + beta);

		}
		else if (scalarProduct(gradf, p) <= c2 * scalarProduct(p, grad))
		{
			alpha = step;
			if (beta > 100)
				step = 2 * alpha;
			else
				step = 0.5*(alpha + beta);
		}
		else
			break;

		i += 1;

	}

	return step;
}
//double golden(double infinum, double supremum, double epsilon)
//{
//	double x1, x2;
//
//	while (fabs(supremum - infinum) > epsilon)
//	{
//		x1 = infinum + (1 - T)*(supremum - infinum);
//		x2 = infinum + T*(supremum - infinum);
//
//		if (func(x1) > func(x2))
//		{
//			infinum = x1;
//			x1 = x2;
//			x2 = infinum + T*(supremum - infinum);
//
//		}
//		else
//		{
//			supremum = x2;
//			x2 = x1;
//			x1 = infinum + (1 - T)*(supremum - infinum);
//
//		}
//	}
//	return (supremum + infinum) / 2;
//}
//double FindMin(double *s, double *p)
//{
//	const double eps = 1e-8;
//	const double tay = 1.618;
//	double a = 0;
//	double b = 1e5;
//	double x0, x1, xf1, xf2;
//	x0 = b - (b - a) / tay; // Ðàñ÷èòûâàåì òî÷êè äåëåíèÿ
//	x1 = a + (b - a) / tay;         //
//P:
//	double t1[2], t2[2];
//	t1[0] = s[0] + x0 * p[0];
//	t1[1] = s[1] + x0 * p[1];
//	t2[0] = s[0] + x1 * p[0];
//	t2[1] = s[1] + x1 * p[1];
//	xf1 = f(t1); // Ðàñ÷èòûâàåì â òî÷êàõ äåëåíèÿ çíà÷åíèå öåëåâîé ôóíêöèè
//	xf2 = f(t2); //
//	if (xf1 >= xf2)
//	{
//		a = x0;
//		x0 = x1;
//		t2[0] = s[0] + x1 * p[0];
//		t2[1] = s[1] + x1 * p[1];
//		xf1 = f(t2);
//		t2[0] = s[0] + x1 * p[0];
//		t2[1] = s[1] + x1 * p[1];
//		x1 = a + (b - a) / tay;
//		xf2 = f(t2);
//	}
//	else
//	{
//		b = x1;
//		x1 = x0;
//		xf2 = xf1;
//		x0 = b - (b - a) / tay;
//		t1[0] = s[0] + x0 * p[0];
//		t1[1] = s[1] + x0 * p[1];
//		xf1 = f(t1);
//	}
//	if (fabs(b - a) < eps)
//	{
//		return (a + b) / 2;
//	}
//	else
//		goto P;
//}

void fletcher_reeves()
{

	double x[2] = { 1.7, 1.7 };
	double min = 10., minx = 10., miny = 10.;
	double xnew[2] = { 0., 0. };

	double e1 = 1e-8;
	double e2 = 1e-8;
	double b = 0.;
	int found = 0, M = 1000, k = 0;

	double grad[2];
	double newgrad[2];

	grad[0] = dfdx1(x);
	grad[1] = dfdx2(x);

	if (sqrt(grad[0] * grad[0] + grad[1] * grad[1]) < e1)
		found = 1;
	
	double d[2] = { -grad[0], -grad[1] };
	double h[2] = { d[0], d[1] };
	double t = wolfe(x, d, 100);

	for (int i = 0; i < N; i++)
		xnew[i] = x[i] + t * h[i];

	while ((found == 0) && (k < M))
	{
		k++;
		newgrad[0] = dfdx1(xnew);
		newgrad[1] = dfdx2(xnew);

		if (k % (10) == 0)
			b = 0; //Обновление
		else
			b = scalarProduct(newgrad, newgrad) / scalarProduct(grad, grad);

		d[0] = -newgrad[0] + b * d[0];
		d[1] = -newgrad[1] + b * d[1];

		if (((xnew[0] - x[0])*(xnew[0] - x[0]) + (xnew[1] - x[1])*(xnew[1] - x[1]) < e2) && (abs(f(xnew) - f(x)) < e2))
			found = 1;
		for (int i = 0; i < N; i++)
		{
			x[i] = xnew[i];
			grad[i] = newgrad[i];
		}
		if (f(x) < min)
		{
			min = f(x);
			minx = x[0];
			miny = x[1];
		}
		t = wolfe(x, d, 100);
		for (int i = 0; i < N; i++)
			xnew[i] = x[i] + t * d[i];
	}
	std::cout << x[0] << " " << x[1] << " " << f(x) << std::endl;
}
int main()
{
	fletcher_reeves();

	system("pause");
	return 0;
}