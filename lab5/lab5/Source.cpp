#include <iostream>
#include <math.H>
using namespace std;

double f(double x) 
{
    return (sqrt(1+2*pow(x, 3)));
}

double sInt(double a, double b, int n) 
{
    double h = (b - a) / n;
    double p1 = 0, p2 = 0;
    for (int i = 1; i < n; ++i) {
        if (i % 2 == 0)
            p1 += f(a + h * i);
        if (i % 2 == 1)
            p2 += f(a + (i + 1) * h);
    }
    return h / 3 * (f(a) + f(b) + 4 * p2 + 2 * p1);
}

double tInt(double a, double b, int n) 
{
    double h = (b - a) / n;
    double p = 0;
    for (int i = 1; i < n; ++i) 
    {
        p += f(a + h * i);
    }

    return h / 2 * (f(a) + 2 * p + f(a + h * n));
}

double intT(double a, double b, double eps, int n) 
{
    double Int1, Int2;
    Int1 = tInt(a, b, n);
    do 
    {
        Int2 = Int1;
        n = n * 2;
        Int1 = tInt(a, b, n);
    } 
    while (abs(Int1 - Int2) <= 3 * eps);
    return Int1;
}

double intS(double a, double b, double eps, int n) 
{
    double Int1, Int2;
    Int1 = sInt(a, b, n);
    do 
    {
        Int2 = Int1;
        n = n * 2;
        Int1 = sInt(a, b, n);
    } 
    while (abs(Int1 - Int2) >= 15 * eps);
    return Int1;
}

double f(double x);
double sInt(double a, double b, int n);
double tInt(double a, double b, int n);
double intT(double a, double b, double eps, int n);
double intS(double a, double b, double eps, int n);

void main()
{
    setlocale(LC_ALL, "Russian");
    double eps = 0.00001, a = 1.2, b = 2.471, int2, int1, e;
    int n = 2;
    int m;
    cout << "1. Интеграл по формуле трапеций" << endl;
    cout << "2. Интеграл по формуле Симпсона" << endl;
    cin >> m;
    cout << endl;
    if (m == 1)
    {
        int2 = intT(a, b, eps, n);
        cout << "Интеграл по формуле трапеций: " << int2;
        int1 = intT(a, b, eps, n * 2);
        e = (int1 - int2) * (pow(0.5, 2) - 1);
        cout << endl << "Погрешность: " << e;
    }
    if (m == 2)
    {
        int2 = intS(a, b, eps, n);
        cout << "Интеграл по формуле Симпсона: " << int2;
        int1 = intS(a, b, eps, n * 2);
        e = (int1 - int2) * (pow(0.5, 4) - 1);
        cout << endl << "Погрешность: " << e;
    }
    if (m != 2 && m != 1)
    {
        cout << "Неверный ввод";
    }
}
