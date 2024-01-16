#include <iostream>
#include <cmath>
// Cholesky decomposition
void cholesky_decomposition(float **A, int n)
{
    float **L = new float *[n];
    for (int i = 0; i < n; i++)
        L[i] = new float[n];

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            float sum = 0;
            if (j == i)
            {
                for (int k = 0; k < j; k++)
                    sum += pow(L[j][k], 2);
                L[j][j] = sqrt(A[j][j] - sum);
            }
            else
            {
                for (int k = 0; k < j; k++)
                    sum += (L[i][k] * L[j][k]);
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
    }

    for (int i = 0; i < n; i++)
        delete[] L[i];
    delete[] L;
}


// LU decomposition
void lu_decomposition(float **A, int n)
{
    float **L = new float *[n];
    float **U = new float *[n];
    for (int i = 0; i < n; i++)
    {
        L[i] = new float[n];
        U[i] = new float[n];
    }

    for (int i = 0; i < n; i++)
    {
        for (int k = i; k < n; k++)
        {
            float sum = 0;
            for (int j = 0; j < i; j++)
                sum += (L[i][j] * U[j][k]);
            U[i][k] = A[i][k] - sum;
        }

        for (int k = i; k < n; k++)
        {
            if (i == k)
                L[i][i] = 1;
            else
            {
                float sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (L[k][j] * U[j][i]);
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        delete[] L[i];
        delete[] U[i];
    }
    delete[] L;
    delete[] U;
}
