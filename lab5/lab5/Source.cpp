#include <iostream>
#include <math.H>
using namespace std;

double f(double x) 
{
    return (sqrt(1+2*pow(x, 3)));
}

double f2(double x, double y) 
{
    return (pow(x, 2)/(1 + pow(y, 2)));
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

double ksInt(double N, double M, double a1, double b1, double c, double d) 
{
    double h, k;
    const int n = 50;
    double z[n][n], ax[n], ans;

    h = (b1 - a1) / (N - 1);
    k = (d - c) / (M - 1);

    for (int i = 0; i < N; ++i) 
    { 
        for (int j = 0; j < M; ++j) 
        {
            z[i][j] = f2(a1 + i * h, c + j * k);
        }
    }
    for (int i = 0; i < N; ++i) 
    {
        ax[i] = 0;
        for (int j = 0; j < M; ++j) 
        {
            if (j == 0 or j == M - 1)
                ax[i] += z[i][j];
            else if (j % 2 == 0)
                ax[i] += 2 * z[i][j]; 
            else
                ax[i] += 4 * z[i][j]; 
        }
        ax[i] *= (k / 3);
    }
    ans = 0;
    for (int i = 0; i < N; ++i) 
    {
        if (i == 0 or i == N - 1)
            ans += ax[i]; 
        else if (i % 2 == 0)  
            ans += 2 * ax[i];
        else
            ans += 4 * ax[i];
    }
    ans *= (h / 3);
    return ans;
}

double f(double x);
double f2(double x, double y);
double sInt(double a, double b, int n);
double tInt(double a, double b, int n);
double intT(double a, double b, double eps, int n);
double intS(double a, double b, double eps, int n);
double ksInt(double N, double M, double a1, double b1, double c, double d);

void main()
{
    setlocale(LC_ALL, "Russian");
    double eps = 0.00001, a = 1.2, b = 2.471, int2, int1, e;
    int n = 2;
    int m;
    cout << "1. Интеграл по формуле трапеций" << endl;
    cout << "2. Интеграл по формуле Симпсона" << endl;
    cout << "3. Интеграл по кубаторной формуле Симпсона" << endl;
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
    if (m == 3)
    {
        double h, k, a1, b1, c, d;
        int N, M;
        cout << "Введите N" << endl;
        cin >> N;
        cout << "Введите M" << endl;
        cin >> M;
        a1 = 0, b1 = 4, c = 1, d = 2;
        cout << "Интеграл по кубаторной формуле Симпсона: " << ksInt(N, M, a1, b1, c, d);
    }
    if (m != 2 && m != 1 && m != 3)
    {
        cout << "Неверный ввод";
    }
    
}
