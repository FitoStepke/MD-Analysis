from gridData import Grid
import numpy as np
from matplotlib import pyplot as plt

# Values for the iterations
t_ca = 3
max_ca = 3  # Assuming ca starts from 1
max_rep = 2

# Initialize lists for cumulative values and individual datasets
cumulative_v = None
datasets = []

for rep in range(max_rep + 1):
    for ca in range(1, max_ca + 1):
        # Construct the DX file name based on rep and ca
        dx_file = "MP_" + str(t_ca) + "_C" + str(ca) + "_R" + str(rep) + ".dx"

        # Read data from the current DX file
        g = Grid(dx_file)

        # Get the number of z slides
        mz = int(g.grid.shape[2]) - 1

        # Initialize the list for the current combination
        v = []

        # Iterate over z slides and sum values
        n = 0
        while n < mz:
            iz = g.grid[:, :, n]
            v.append(np.sum(iz))
            n += 1

        # If it's the first iteration, initialize the cumulative list
        if cumulative_v is None:
            cumulative_v = v
        else:
            # Sum the corresponding elements of the current and cumulative lists
            cumulative_v = [x + y for x, y in zip(cumulative_v, v)]

        # Store individual dataset for calculating standard deviation later
        datasets.append(v)

# Divide each value in cumulative_v by t_ca * 3
cumulative_v_normalized = [value / (t_ca * 2) for value in cumulative_v]

# Calculate standard deviation for each point along the z-axis
std_deviation = np.std(datasets, axis=0)

# Plot the normalized cumulative values
plt.figure(figsize=(10, 15))
plt.plot(cumulative_v_normalized, range(mz), color="C0", label='Cumulative Normalized')

# Plot the standard deviation using fill_between
plt.fill_betweenx(range(mz), np.subtract(cumulative_v_normalized, std_deviation), np.add(cumulative_v_normalized, std_deviation), color="red", alpha=0.5, label='Standard Deviation')

plt.xlim(0, 0.3)
plt.axhline(y=20, color='red', linestyle='dashed', linewidth=1)
plt.axhline(y=40, color='red', linestyle='dashed', linewidth=1) 
plt.title("3 $\mathregular{Ca^{2+}}$ model with Standard Deviation", fontsize=28)
plt.xlabel("Probability in Z axis", fontsize=24)
plt.ylabel("Z ($\mathregular{A^{0}}$)", fontsize=24)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend()
plt.show()
