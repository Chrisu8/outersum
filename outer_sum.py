from main import*
# Define the size of LD matrix
matrix_size = 200

# Define the sizes of the blocks of 1s
block_sizes = [4, 5, 10, 15, 3, 8, 10, 20,35,10,27,6,28,19]

# Create the LD matrix with the defined block structure
ld_matrix = create_ld_matrix(matrix_size, block_sizes)

outer_sums = calculate_outer_sums(ld_matrix)
print("Outer sums for diagonal elements:", outer_sums)

print(hann_window(20))

# Generate Hann window
window_size = 20
hann = hann_window(window_size)

# Perform the convolution
convolved_signal = np.convolve(hann/hann.sum(),outer_sums, mode='same')


minima_a = sig.argrelextrema(convolved_signal, np.less)[0]
minima_a_vals = [convolved_signal[i] for i in minima_a]
print(minima_a)
print(minima_a_vals)

# Plot the result
plt.figure(figsize=(14, 6))
plt.plot(outer_sums, label='Anti-diagonal Sums')
plt.plot(hann, label='Hann Window')
plt.plot(convolved_signal, label='Convolved Signal')
plt.scatter(minima_a, minima_a_vals)
plt.legend()
plt.title('Convolution of Anti-diagonal Sums with Hann Window')
plt.xlabel('position on chromosome')
plt.ylabel('sum r^2')
plt.grid(True)
plt.show()
