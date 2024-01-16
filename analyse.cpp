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