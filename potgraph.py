from gridData import Grid
import numpy as np
import matplotlib.pyplot as plt

# Conversion factor from kT/e to Volts
CONVERSION_FACTOR = 0.0258

# Function to calculate the mean and standard deviation of electrostatic potentials in the XY-plane for each Z-slice
def calculate_potential_stats(file_path):
    grid = Grid(file_path)  # Load the .dx file
    mean_xy_per_z = np.mean(grid.grid, axis=(0, 1))  # Mean over XY-plane for each Z
    std_xy_per_z = np.std(grid.grid, axis=(0, 1))   # Standard deviation over XY-plane for each Z
    return mean_xy_per_z, std_xy_per_z

# Process both .dx files
mean_with_ca2, std_with_ca2 = calculate_potential_stats("ca.dx")  # System with Ca2+
mean_without_ca2, std_without_ca2 = calculate_potential_stats("noca.dx")  # System without Ca2+

# Convert values from kT/e to Volts
mean_with_ca2 *= CONVERSION_FACTOR
std_with_ca2 *= CONVERSION_FACTOR
mean_without_ca2 *= CONVERSION_FACTOR
std_without_ca2 *= CONVERSION_FACTOR

# Generate adjusted Z indices
num_slices = len(mean_with_ca2)  # Assuming both grids have the same number of Z-slices
z_values = np.linspace(-30, 30, num_slices)  # Adjusted Z values centered around 0

# Plot the results
plt.figure(figsize=(10, 12))

# System with Ca2+
plt.plot(mean_with_ca2, z_values, label='With Ca2+', color='red', linestyle='-')
plt.fill_betweenx(z_values, mean_with_ca2 - std_with_ca2, mean_with_ca2 + std_with_ca2, color='red', alpha=0.3)

# System without Ca2+
plt.plot(mean_without_ca2, z_values, label='Without Ca2+', color='blue', linestyle='-')
plt.fill_betweenx(z_values, mean_without_ca2 - std_without_ca2, mean_without_ca2 + std_without_ca2, color='blue', alpha=0.3)

# Configure the plot
plt.xlabel('Average electrostatic potential (Volts)')
plt.ylabel('Z (A°)')
plt.title('Electrostatic Surface Potential')
plt.legend()
plt.grid()

plt.show()
