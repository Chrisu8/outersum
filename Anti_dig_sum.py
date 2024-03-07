import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import convolve
import scipy.signal as sig

# Load your data
data = np.loadtxt('example_data/vector-EUR-chr2-39967768-40067768.txt', delimiter='\t')
first_column = data[:, 0]
second_column = data[:, 1]

# Create a Hann window of a desired size
window_size = 25  # Example size, adjust as needed
hann_window = np.hanning(window_size)

# Perform the convolution
convolved_data = convolve(second_column, hann_window, mode='same') / sum(hann_window)

minima_indices = sig.argrelextrema(convolved_data, np.less)[0]
minima_x = first_column[minima_indices]
minima_y = convolved_data[minima_indices]
# Create a scatter plot of original data and a plot of convolved data
plt.scatter(first_column, second_column, marker='o', s=5, label='Original Data')
plt.plot(first_column, convolved_data, color='red', label='Convolved with Hann Window')

# Add labels, a title, and a legend
plt.xlabel('Position')
plt.ylabel('Antisum')
plt.title('Chr2, 39967768-40067768 - Convolution with Hann Window')
plt.scatter(minima_x, minima_y, color='green', marker='x', label='Local Minima')
plt.legend()

# Show the plot
plt.show()


print("index",minima_x)
print("antisum",minima_y)



