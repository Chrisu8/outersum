import numpy as np
def hann_window(N):
    # N is the filter width
    return 0.5 * (1 - np.cos(2 * np.pi * np.arange(N) / (N - 1)))
def one_block(size):
    return np.ones((size, size), dtype=int)

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