import os
import unittest
from outersum.core import read_data, create_cov_matrix, cov_cor_matrix, calculate_outer_sums, convolve_with_hann_window
import numpy as np


class TestOuterSum(unittest.TestCase):

    def setUp(self):

        base_dir = os.path.dirname(os.path.abspath(__file__))
        self.sample_data_path = os.path.join(base_dir, 'sample_data')

    def test_read_data(self):
        data, positions = read_data(self.sample_data_path)
        expected_data_length = 10
        expected_positions_length = 10
        self.assertEqual(len(data), expected_data_length)
        self.assertEqual(len(positions), expected_positions_length)

    def test_create_cov_matrix(self):
        data, positions = read_data(self.sample_data_path)
        position_map = {pos: idx for idx, pos in enumerate(positions)}
        cov_matrix = create_cov_matrix(data, position_map)
        expected_matrix_shape = (len(positions), len(positions))
        self.assertEqual(cov_matrix.shape, expected_matrix_shape)


    def test_cov_cor_matrix(self):
        data, positions = read_data(self.sample_data_path)
        position_map = {pos: idx for idx, pos in enumerate(positions)}
        cov_matrix = create_cov_matrix(data, position_map)
        cor_matrix = cov_cor_matrix(cov_matrix)
        np.testing.assert_array_almost_equal(np.diag(cor_matrix), np.ones(cor_matrix.shape[0]))

    def test_calculate_outer_sums(self):
        matrix = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        outer_sums = calculate_outer_sums(matrix)
        expected_outer_sums = [0, 3, 0]
        self.assertEqual(outer_sums, expected_outer_sums)

    def test_convolve_with_hann_window(self):
        data = np.array([1, 2, 3, 4, 5])  # Example data
        window_size = 3
        convolved_data = convolve_with_hann_window(data, window_size)
        self.assertTrue(len(convolved_data), len(data))


if __name__ == '__main__':
    unittest.main()
