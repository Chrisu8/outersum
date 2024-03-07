import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from .utils import hann_window, one_block, anti_diagonal_sums  # Adjust based on your actual utility functions


# Data Reading and Processing
def read_data(file_path):
    data = []
    positions = set()

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            pos1 = int(parts[2])
            pos2 = int(parts[3])
            shrinkage_cov = float(parts[7])
            empirical_cov = float(parts[6])

            data.append((pos1, pos2, empirical_cov,shrinkage_cov))
            positions.update([pos1, pos2])

    return data, sorted(positions)

#Matrix Creation
def create_cov_matrix(data, position_map):
    size = len(position_map)
    matrix = np.zeros((size, size))

    for pos1, pos2, empirical_cov, shrinkage_cov in data:
        i = position_map[pos1]
        j = position_map[pos2]

        if i == j:
            matrix[i][j] = empirical_cov
        else:
            matrix[i][j] = shrinkage_cov
            matrix[j][i] = shrinkage_cov

    return matrix

def cov_cor_matrix(cov_matrix,):
    correlation_matrix = np.zeros(cov_matrix.shape)
    # Iterate over the rows
    for i in range(cov_matrix.shape[0]):
        # Iterate over the columns
        for j in range(cov_matrix.shape[1]):
            # Compute r^2_ij for off-diagonal elements
            if i != j:
                correlation_matrix[i][j] = (cov_matrix[i][j] ** 2) / (cov_matrix[i][i] * cov_matrix[j][j])
            else:
                # Set diagonal elements to 1
                correlation_matrix[i][j] = 1

    return correlation_matrix

def calculate_outer_sums(matrix):
    n = len(matrix)
    outer_sums = []

    for i in range(n):
        sum = 0
        for row in range(i):
            for col in range(i + 1, n):
                sum += matrix[row][col]
        outer_sums.append(sum)

    return outer_sums

def convolve_with_hann_window(data, window_size):
    window = hann_window(window_size)
    normalized_window = window / window.sum()
    convolved_data = np.convolve(normalized_window, data, mode='same')
    return convolved_data

def plot_data(positions, original_data, convolved_data):
    plt.figure(figsize=(12, 6))
    plt.plot(positions, original_data, label='Original Outer Sums', marker='o')
    plt.plot(positions, convolved_data, label='Convolved Data', marker='x')
    plt.title("Outer Sums of Chr2.39967768.40067768 and Convolved Data")
    plt.xlabel("Chr2_Position")
    plt.ylabel("outer sum")
    plt.legend()
    plt.grid(True)
    plt.show()