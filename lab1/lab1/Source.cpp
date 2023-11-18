#include<iostream>
#include<cmath>
using namespace std;

void input(double** a, double* b, int n, double **a1, double **a2, double *b2)
{
	cout << "Еnter matrix: " << endl;
	for (int i = 0; i < n; ++i)
	{
		a[i] = new double[n];
		for (int j = 0; j < n; ++j)
		{
			cin >> a[i][j];
		}
	}
	for (int i = 0; i < n; ++i)
	{
		cin >> b[i];
	}
	for (int i = 0; i < n; ++i)
	{
		b2[i] = b[i];
	}
	for (int i = 0; i < n; ++i)
	{
		a1[i] = new double[n];
		for (int j = 0; j < n; ++j)
		{
			a1[i][j] = a[i][j];
		}
	}
	for (int i = 0; i < n; ++i)
	{
		a2[i] = new double[n];
		for (int j = 0; j < n; ++j)
		{
			a2[i][j] = a[i][j];
		}
	}
}

void output(double** a, double* b, int n)
{
	cout << "Matrix: " << endl;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			cout << a[i][j] << " ";
		}
		cout << b[i] << endl;
	}
}

void gauss(double** a, double* b, int n, double* x)
{
	double max;
	int ind;
	int k = 0;
	for (k = 0; k < n; ++k)
	{
		max = abs(a[k][k]);
		ind = k;

		for (int i = k + 1; i < n; ++i)
		{
			if (abs(a[i][k]) > max)
			{
				max = abs(a[i][k]);
				ind = i;
			}
		}

		for (int j = 0; j < n; ++j)
		{
			double t = a[k][j];
			a[k][j] = a[ind][j];
			a[ind][j] = t;
		}

		double t = b[k];
		b[k] = b[ind];
		b[ind] = t;

		for (int i = k; i < n; ++i)
		{
			double t = a[i][k];
			for (int j = 0; j < n; ++j)
			{
				a[i][j] = a[i][j] / t;
			}
			b[i] = b[i] / t;
			if (i == k)  continue;
			for (int j = 0; j < n; ++j)
			{
				a[i][j] = a[i][j] - a[k][j];
			}
			b[i] = b[i] - b[k];
		}
	}

	for (k = n - 1; k >= 0; --k)
	{
		x[k] = b[k];
		for (int i = 0; i < k; ++i)
			b[i] = b[i] - a[i][k] * x[k];
	}

	cout << "X: " << endl;
	for (int i = 0; i < n; ++i)
	{
		cout << x[i] << " ";
	}
	cout << endl;
}

void vector(double** a2, double* b2, int n, double* x)
{
	double* v;
	v = new double[n];
	double maxp = 0;
	for (int i = 0; i < n; ++i)
	{
		v[i] = b2[i];
		for (int j = 0; j < n; ++j)
		{
			v[i] -= a2[i][j] * x[j];
		}
		if (abs(v[i]) > maxp)
		{
			maxp = abs(v[i]);
		}
	}
	cout << "Vector norm: " << maxp << endl;

	cout << "Vector: " << endl;
	for (int i = 0; i < n; ++i)
	{
		cout << v[i] << " ";
	}
	cout << endl;
}

void relative_error(int n, double* x, double** a1, double *x1)
{

	double max;
	int ind;
	int k = 0;
	double *b1;
	b1 = new double[n];
	for (int i = 0; i < n; ++i)
	{
		b1[i] = 0;
		for (int j = 0; j < n; ++j)
		{
			b1[i] += a1[i][j] * x[j];
		}
	}
	for (k = 0; k < n; ++k)
	{
		max = abs(a1[k][k]);
		ind = k;

		for (int i = k + 1; i < n; ++i)
		{
			if (abs(a1[i][k]) > max)
			{
				max = abs(a1[i][k]);
				ind = i;
			}
		}

		for (int j = 0; j < n; ++j)
		{
			double t = a1[k][j];
			a1[k][j] = a1[ind][j];
			a1[ind][j] = t;
		}

		double t = b1[k];
		b1[k] = b1[ind];
		b1[ind] = t;

		for (int i = k; i < n; ++i)
		{
			double t = a1[i][k];
			for (int j = 0; j < n; ++j)
			{
				a1[i][j] = a1[i][j] / t;
			}
			b1[i] = b1[i] / t;
			if (i == k)  continue;
			for (int j = 0; j < n; ++j)
			{
				a1[i][j] = a1[i][j] - a1[k][j];
			}
			b1[i] = b1[i] - b1[k];
		}
	}
	for (k = n - 1; k >= 0; --k)
	{
		x1[k] = b1[k];
		for (int i = 0; i < k; ++i)
			b1[i] = b1[i] - a1[i][k] * b1[k];
	}

	cout << "X1: " << endl;
	for (int i = 0; i < n; ++i)
	{
		cout << x1[i] << " ";
	}
	cout << endl;

	double max1 = 0;
	double max2 = 0;
	for (int i = 0; i < n; ++i)
	{
		if (abs(x1[i] - x[i]) > max2)
		{
			max2 = abs(x1[i] - x[i]);
		}
		if (abs(x[i]) > max1)
		{
			max1 = x[i];
		}
	}
	cout << "Error: " << max2 / max1 << endl;
}

int main()
{
	int n = 3;
	double* x, ** a, * b, **a1, *x1, **a2, *b2;
	a = new double* [n];
	a1 = new double* [n];
	a2 = new double* [n];
	b = new double[n];
	b2 = new double[n];
	x = new double[n];
	x1 = new double[n];
	input(a, b, n, a1, a2, b2);
	cout << endl;
	output(a, b, n);
	cout << endl;
	gauss(a, b, n, x);
	cout << endl;
	vector(a2, b2, n, x);
	cout << endl;
	relative_error(n, x, a1, x1);
}
