# the main function here is used for each party t get their outersum
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

def hann_window(N):
    # N is the filter width
    return 0.5 * (1 - np.cos(2 * np.pi * np.arange(N) / (N - 1)))

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

#  create a square block
def one_block(size):
    return np.ones((size, size), dtype=int)
def create_ld_matrix(matrix_size, block_sizes):
    ld_matrix = np.zeros((matrix_size, matrix_size), dtype=int)

    current_pos = 0
    for block_size in block_sizes:
        ld_matrix[current_pos:current_pos+block_size, current_pos:current_pos+block_size] = one_block(block_size)
        current_pos += block_size  # move the starting position for the next block

    return ld_matrix

def anti_diagonal_sums(matrix):
    n = len(matrix)
    anti_diag_sums = []

    # Process each anti-diagonal
    for k in range(1, 2 * n):
        vk = 0
        for i in range(1, k + 1):
            j = k - i + 1
            # Check the bounds
            if 1 <= i <= n and 1 <= j <= n:
                vk += matrix[i - 1][j - 1]  # Subtracting 1 as matrix indices start from 0
        anti_diag_sums.append(vk)

    return np.array(anti_diag_sums)

# Function to plot the anti-diagonal sums
def plot_anti_diagonal_sums(anti_diag_sums):
    # Create an array of indices, starting from 1
    indices = np.arange(1, len(anti_diag_sums) + 1)

    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(indices, anti_diag_sums, marker='o')

    # Add title and labels
    plt.title('Anti-diagonal Sums of the Matrix')
    plt.xlabel('position on chromosome')
    plt.ylabel('Sum r^2')

    # Show grid lines
    plt.grid(True)

    # Show the plot
    plt.show()



file_path = 'example_data/chr2.39967768.40067768'
data, positions = read_data(file_path)
position_map = {pos: idx for idx, pos in enumerate(positions)}
cov_matrix = create_cov_matrix(data, position_map)
matrix = cov_cor_matrix(cov_matrix)
outer_sums = calculate_outer_sums(matrix)

