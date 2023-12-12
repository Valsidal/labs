#include <iostream>
#include <math.h>
#include<iomanip>
using namespace std;

double f1(double x1, double x2);
double f2(double x1, double x2);
double f1d11(double x1, double x2, double dif);
double f1d12(double x1, double x2, double dif);
double f2d21(double x1, double x2, double dif);
double f2d22(double x1, double x2, double dif);
void resV(double* v1, double x1, double x2);
double* gauss(double** m, double* y, const int n, float accur);
void output(double** m, double* y, const int n);

void main()
{
    double x1 = 1, x2 = -1;
    int k = 1, N = 100;
    double accur = 1e-9, dif = 0.01;
    const int n = 2;
    double* v1, * x, * d, * sol, * v2;
    v1 = new double[n];
    v2 = new double[n];
    d = new double[n];
    x = new double[n];
    sol = new double[n];
    sol[0] = x1;
    sol[1] = x2;
    double** a = new double* [n];
    for (int i = 0; i < n; i++)
    { 
        a[i] = new double[n];
    }
    double** a2 = new double* [n];
    for (int i = 0; i < n; i++) 
    {
        a2[i] = new double[n];
    }
    do
    {
        resV(v1, x1, x2); 
        a[0][0] = f1d11(x1, x2, dif);
        a[1][0] = f1d12(x1, x2, dif);
        a[0][1] = f2d21(x1, x2, dif);
        a[1][1] = f2d22(x1, x2, dif);
        for (int i = 0; i < n; i++)
        { 
            for (int j = 0; j < n; j++)
            {
                a2[i][j] = a[i][j];
            }
        }
        for (int i = 0; i < n; i++)
        {
            v2[i] = v1[i];
        }
        output(a2, v2, n);
        d = gauss(a2, v1, n, accur); 
        for (int i = 0; i < n; i++)
            sol[i] += d[i]; 
        for (int i = 0; i < n; i++)
        {
            cout << "delta" << i + 1 << " = " << d[i] << " " << endl;
        }
        double max1 = 0;
        double max2 = 0;
        for (int i = 0; i < n; i++)  
        {
            if (abs(v1[i]) > max1)
            {
                max1 = abs(v1[i]);
            }
            if (abs(sol[i]) < 1)
            {
                if (abs(d[i]) > max2)
                    max2 = abs(d[i]);
            }
            if (abs(d[i] >= 1))
            {
                if (abs(d[i] / sol[i]) > max2)
                    max2 = abs(d[i]);
            }
        }

        d[0] = max1;
        d[1] = max2;
        cout << endl;
        x1 = sol[0];
        x2 = sol[1];
        k++;
        if (k >= N)
        {
            break;
        }

    } 
    while (d[0] > accur || d[1] > accur);

    x1 = sol[0];
    x2 = sol[1];
    cout << "____________________________________________________________" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << sol[i] << endl;
    }
    delete[] v1;
    delete[] v2;
    delete[] x;
    delete[] d;
    delete[] sol;
    for (int i = 0; i < n; i++)
    {
        delete[] a[i];
    }
    for (int i = 0; i < n; i++)
    {
        delete[] a2[i];
    }
    delete[] a;
    delete[] a2;
}

double f1(double x1, double x2)
{
    return (cos(0.4*x2+pow(x1,2)) + pow(x2, 2) + pow(x1, 2) - 1.6);
}
double f2(double x1, double x2)
{
    return (1.5*pow(x1,2)-(pow(x2,2)/0.36)-1);
}
double f1d11(double x1, double x2, double dif)
{
    cout << ((f1(x1 + dif, x2) - f1(x1, x2)) / (dif));
    return((f1(x1 + dif, x2) - f1(x1, x2)) / (dif));
}
double f1d12(double x1, double x2, double dif)
{
    cout << ((f1(x1, x2 + dif) - f1(x1, x2)) / (dif));
    return((f1(x1, x2 + dif) - f1(x1, x2)) / (dif));
}
double f2d21(double x1, double x2, double dif)
{
    cout << ((f2(x1 + dif, x2) - f2(x1, x2)) / (dif));
    return((f2(x1 + dif, x2) - f2(x1, x2)) / (dif));
}
double f2d22(double x1, double x2, double dif)
{
    cout << ((f2(x1, x2 + dif) - f2(x1, x2)) / (dif));
    return((f2(x1, x2 + dif) - f2(x1, x2)) / (dif));
}

void resV(double* vNev, double x1, double x2)
{
    vNev[0] = -f1(x1, x2);
    vNev[1] = -f2(x1, x2);
}
void output(double** m, double* y, const int n)
{
    cout << endl;
    cout << "matrix" << endl;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cout << m[i][j] << " ";
        }
        cout << y[i] << endl;
    }
}
double* gauss(double** m, double* y, const int n, float accur)
{
    double* x, max, roundH;
    int k, index;
    const double eps = 0.00001;  
    x = new double[n];
    k = 0;
    while (k < n)
    {
        max = abs(m[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++)
        {
            if (abs(m[i][k]) > max)
            {
                max = abs(m[i][k]);
                index = i;
            }
        }
        for (int j = 0; j < n; j++)
        {
            swap(m[k][j], m[index][j]);
        }
        swap(y[k], y[index]);
        for (int i = k; i < n; i++)
        {
            double temp = m[i][k];
            if (abs(temp) < eps) continue; 
            for (int j = 0; j < n; j++) 
            {
                m[i][j] = m[i][j] / temp;
            }
            y[i] = y[i] / temp;
            if (i == k)  continue;
            for (int j = 0; j < n; j++)
                m[i][j] = m[i][j] - m[k][j];
            y[i] = y[i] - y[k];
        }
        k++;
    }
    for (k = n - 1; k >= 0; k--)
    {
        x[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - m[i][k] * x[k];
    }
    return x;
}
