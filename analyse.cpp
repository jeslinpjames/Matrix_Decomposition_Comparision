#include <iostream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <iomanip>

using namespace std;
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


// QR decomposition
void qr_decomposition(float **A, int n)
{
    float **Q = new float *[n];
    float **R = new float *[n];
    for (int i = 0; i < n; i++)
    {
        Q[i] = new float[n];
        R[i] = new float[n];
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            Q[i][j] = 0;
            R[i][j] = 0;
        }
    }

    for (int i = 0; i < n; i++)
    {
        float norm = 0;
        for (int j = 0; j < n; j++)
            norm += A[j][i] * A[j][i];
        norm = sqrt(norm);
        R[i][i] = norm;

        for (int j = 0; j < n; j++)
            Q[j][i] = A[j][i] / norm;

        for (int j = i + 1; j < n; j++)
        {
            float dot = 0;
            for (int k = 0; k < n; k++)
                dot += Q[k][i] * A[k][j];
            R[i][j] = dot;

            for (int k = 0; k < n; k++)
                A[k][j] -= Q[k][i] * dot;
        }
    }

    for (int i = 0; i < n; i++)
    {
        delete[] Q[i];
        delete[] R[i];
    }
    delete[] Q;
    delete[] R;
}

void read_positive_definite_matrix(float **resultMatrix, int dimension, const string &filename)
{
    // Open the CSV file
    ifstream inputFile(filename);
    string line;

    // Read the file line by line
    for (int i = 0; i < dimension && getline(inputFile, line); ++i)
    {
        stringstream lineStream(line);
        string value;

        // Parse each value and store it in the matrix
        for (int j = 0; j < dimension && getline(lineStream, value, ','); ++j)
        {
            resultMatrix[i][j] = stof(value);
        }
    }

    // Close the file
    inputFile.close();
}


int main()
{
    const int matrixSize = 5; // Change this to the desired matrix dimension
    float **matrix = new float *[matrixSize];
    for (int i = 0; i < matrixSize; ++i)
    {
        matrix[i] = new float[matrixSize];
    }

    const string filename = "your_matrix_file.csv"; // Change this to your CSV file path

    // Read positive definite matrix from CSV
    read_positive_definite_matrix(matrix, matrixSize, filename);

    // Measure Cholesky decomposition time
    auto startCholesky = chrono::high_resolution_clock::now();
    cholesky_decomposition(matrix, matrixSize);
    auto endCholesky = chrono::high_resolution_clock::now();
    auto durationCholesky = chrono::duration_cast<chrono::microseconds>(endCholesky - startCholesky).count();

    // Measure QR decomposition time
    auto startQR = chrono::high_resolution_clock::now();
    qr_decomposition(matrix, matrixSize);
    auto endQR = chrono::high_resolution_clock::now();
    auto durationQR = chrono::duration_cast<chrono::microseconds>(endQR - startQR).count();

    // Measure LU decomposition time
    auto startLU = chrono::high_resolution_clock::now();
    lu_decomposition(matrix, matrixSize);
    auto endLU = chrono::high_resolution_clock::now();
    auto durationLU = chrono::duration_cast<chrono::microseconds>(endLU - startLU).count();

    // Output results
    cout << "Cholesky Decomposition Time: " << durationCholesky << " microseconds" << endl;
    cout << "QR Decomposition Time: " << durationQR << " microseconds" << endl;
    cout << "LU Decomposition Time: " << durationLU << " microseconds" << endl;

    // Calculate ratio of run times
    double ratioCholeskyQR = static_cast<double>(durationCholesky) / durationQR;
    double ratioCholeskyLU = static_cast<double>(durationCholesky) / durationLU;

    cout << "Ratio of Cholesky to QR to LU Decomposition Time: "
         << fixed << setprecision(2) << ratioCholeskyQR << " : " << ratioCholeskyLU << " : " << ratioCholeskyLU <<endl;

    for (int i = 0; i < matrixSize; ++i)
    {
        delete[] matrix[i];
    }
    delete[] matrix;

    return 0;
}