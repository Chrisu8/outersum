from main import*

file_path = 'example_data/chr2.39967768.40067768'
data, positions = read_data(file_path)
position_map = {pos: idx for idx, pos in enumerate(positions)}
cov_matrix = create_cov_matrix(data, position_map)
matrix = cov_cor_matrix(cov_matrix)
outer_sums = calculate_outer_sums(matrix)

window_size = 200  # You can adjust this as needed
convolved_sums = convolve_with_hann_window(outer_sums, window_size)

print("Original Outer Sums:", outer_sums)
print("Convolved Outer Sums:", convolved_sums)

plot_data(positions, outer_sums, convolved_sums)

minima_a = sig.argrelextrema(convolved_sums, np.less)[0]
minima_a_vals = [convolved_sums[i] for i in minima_a]
print("index",minima_a)
print("position",positions[int(minima_a)])
print("value",minima_a_vals)
