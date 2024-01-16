import numpy as np
import csv 

def generate_symmetric_positive_matrix(n):
    # Generate a random matrix
    A = np.random.randint(0, 10, size=(n, n))
    # Make the matrix symmetric by adding its transpose
    symmetric_matrix = A + A.T

    # Make the matrix positive definite
    positive_definite_matrix = np.dot(symmetric_matrix, symmetric_matrix.T)

    return positive_definite_matrix
 
n = int(input("Enter the dimension of matrix : "))
spd_matrix = generate_symmetric_positive_matrix(n)

with open('spd_matrix.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for row in spd_matrix:
        writer.writerow(row)