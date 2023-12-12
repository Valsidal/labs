#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

void gauss(vector<vector<double>> A, vector<double> B, vector<double>& x)
{
    const int n = 2;
    for (int i = 0; i < n; i++)
    {
        B[i] *= -1;
    }
    for (int i = 0; i < n; i++)
    {
        int maxRow = i;
        for (int k = i + 1; k < n; k++)
        {
            if (abs(A[k][i]) > abs(A[maxRow][i]))
            {
                maxRow = k;
            }
        }
        swap(A[maxRow], A[i]);
        swap(B[maxRow], B[i]);

        for (int k = i + 1; k < n; k++)
        {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; j++)
            {
                A[k][j] -= factor * A[i][j];
            }
            B[k] -= factor * B[i];
        }
    }
    for (int i = 0; i < 2; i++)
    {
        x.push_back(B[i]);
    }
    for (int i = 1; i >= 0; i--)
    {
        for (int j = i + 1; j < 2; j++)
        {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
}

vector<double> f(vector<double> x)
{
    vector<double> F;
    F.push_back(cos(0.4 * x[1] + pow(x[0], 2)) + pow(x[1], 2) + pow(x[0], 2) - 1.6);
    F.push_back(1.5 * pow(x[0], 2) - (pow(x[1], 2) / 0.36) - 1);
    return F;
}

void jacobi(double x1, double x2, vector<vector<double>>& jac)
{
    jac[0][0] = -2 * x1 * sin(0.4 * x2 + pow(x1, 2)) + 2 * x1;
    jac[0][1] = -0.4 * sin(0.4 * x2 + pow(x1, 2)) + 2 * x2;
    jac[1][0] = 3 * x1;
    jac[1][1] = -2 * x2 / 0.36;
}

void anotherJacobi(vector<double> x, vector<vector<double>>& jac, double m)
{
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            vector<double> xPlus = x;
            xPlus[j] += xPlus[j] * m;
            jac[i][j] = (f(xPlus)[i] - f(x)[i]) / (x[j] * m);
        }
    }
}

void newton(double x1, double x2, double m = 0.0, double eps = 1e-9, int max_iter = 100)
{
    vector<double> x = { x1,x2 };
    cout << "Начальное приближение: " << x1 << ", " << x2 << endl;
    cout << "ε1 = " << eps << "; ε2 = " << eps << "; max_iter = " << max_iter << endl;
    cout << "итерация  1  2 " << endl;
    int k;
    for (k = 0; k < max_iter; k++)
    {
        vector<vector<double>> jac = { {0,0},{0,0} };
        vector<double> F = f(x);
        vector<double> dx;
        if (m > 0)
        {
            anotherJacobi(x, jac, m);
        }
        else
        {
            jacobi(x[0], x[1], jac);
        }
        gauss(jac, F, dx);
        double maxF;
        double maxGap = 0;
        if (abs(F[0]) > abs(F[1]))
        {
            maxF = abs(F[0]);
        }
        else
        {
            maxF = abs(F[1]);
        }
        for (int i = 0; i < 2; i++)
        {
            double gap;
            if (abs(x[i] + dx[i]) < 1)
            {
                gap = abs(dx[i]);
            }
            else
            {
                gap = abs((dx[i]) / x[i] + dx[i]);
            }
            if (maxGap < gap)
            {
                maxGap = gap;
            }
        }
        for (int i = 0; i < 2; i++)
        {
            x[i] += dx[i];
        }
        cout << k << " " << maxF << " " << maxGap << endl;
        if (maxF <= eps || maxGap <= eps)
        {
            cout << "Решение: " << x[0] << "; " << x[1] << endl;
            break;
        }
    }
    if (k == max_iter)
    {
        cout << "IER = 2";
    }
}

int main()
{
    setlocale(LC_ALL, "RUS");
    newton(1, -1);
    cout << endl;
    newton(-1, 1);
    cout << endl;
    cout << "M = 0.01" << endl;
    newton(1, -1, 0.01);
    cout << endl;
    newton(-1, 1, 0.01);
    cout << endl;
    cout << "M = 0.05" << endl;
    newton(1, -1, 0.05);
    cout << endl;
    newton(-1, 1, 0.05);
    cout << endl;
    cout << "M = 0.1" << endl;
    newton(1, -1, 0.1);
    cout << endl;
    newton(-1, 1, 0.1);
}
