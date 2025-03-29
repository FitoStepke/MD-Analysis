import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

def calculate_and_plot_distances(psf_file, dcd_file, output_file, smoothing_sigma=5):
    """
    Calculate distances over time between pairs of protein segment centers of mass
    and between each segment and the combined center of mass. Generate corresponding
    plots and calculate average and maximum distances.

    Parameters:
        psf_file (str): Path to the PSF file.
        dcd_file (str): Path to the DCD file.
        output_file (str): Path to save the plots and results.
        smoothing_sigma (float): Sigma value for Gaussian smoothing.
    """
    # Load the trajectory
    u = mda.Universe(psf_file, dcd_file)

    # Define protein selections
    protein_selections = {
        "AB": "protein and (segid PROA or segid PROB) and resid 42 43 47 48",
        "BC": "protein and (segid PROB or segid PROC) and resid 42 43 47 48",
        "CD": "protein and (segid PROC or segid PROD) and resid 42 43 47 48",
        "DE": "protein and (segid PROD or segid PROE) and resid 42 43 47 48",
        "EF": "protein and (segid PROE or segid PROF) and resid 42 43 47 48",
        "FA": "protein and (segid PROF or segid PROA) and resid 42 43 47 48"
    }

    # Parse selections into MDAnalysis format
    protein_groups = {key: u.select_atoms(sel) for key, sel in protein_selections.items()}

    # Define pairs to calculate distances
    pairs = [("AB", "BC"), ("BC", "CD"), ("CD", "DE"), ("DE", "EF"), ("EF", "FA"), ("FA", "AB")]

    # Initialize distance arrays
    distances = {pair: [] for pair in pairs}
    distances_to_combined = {key: [] for key in protein_selections.keys()}

    # Iterate over each frame in the trajectory
    for ts in u.trajectory:
        # Calculate distances between pairs
        for pair in pairs:
            group1, group2 = pair
            com1 = protein_groups[group1].center_of_mass()
            com2 = protein_groups[group2].center_of_mass()
            distance = np.linalg.norm(com1 - com2)
            distances[pair].append(distance)

        # Calculate the combined center of mass
        combined_com = np.mean([protein_groups[key].center_of_mass() for key in protein_selections.keys()], axis=0)

        # Calculate distances of each selection to the combined center
        for key in protein_selections.keys():
            group_com = protein_groups[key].center_of_mass()
            distance_to_combined = np.linalg.norm(group_com - combined_com)
            distances_to_combined[key].append(distance_to_combined)

    # Apply Gaussian smoothing to both sets of distances
    smoothed_distances = {pair: gaussian_filter1d(dist, sigma=smoothing_sigma) for pair, dist in distances.items()}
    smoothed_distances_to_combined = {key: gaussian_filter1d(dist, sigma=smoothing_sigma) for key, dist in distances_to_combined.items()}

    # Calculate average and maximum distances
    results = []
    for pair, dist in distances.items():
        avg_distance = np.mean(dist)
        max_distance = np.max(dist)
        results.append(f"{pair[0]}-{pair[1]}: Avg = {avg_distance:.2f} Å, Max = {max_distance:.2f} Å")
    for key, dist in distances_to_combined.items():
        avg_distance = np.mean(dist)
        max_distance = np.max(dist)
        results.append(f"{key} to Combined Center: Avg = {avg_distance:.2f} Å, Max = {max_distance:.2f} Å")

    # Save the results to a text file
    results_file = output_file.replace(".png", "_distances.txt")
    with open(results_file, "w") as f:
        f.write("\n".join(results))
    print(f"Results saved to {results_file}")

    # Plot distances between pairs
    plt.figure(figsize=(12, 8))
    for pair, dist in smoothed_distances.items():
        plt.plot(dist, label=f"{pair[0]}-{pair[1]}")
    plt.title("Smoothed Distances Between Protein Segments Over Time")
    plt.xlabel("Frame")
    plt.ylabel("Distance (Å)")
    plt.legend()
    plt.grid(True)
    plt.savefig(output_file.replace(".png", "_pair_distances.png"), dpi=300)
    plt.show()

    # Plot distances to combined center of mass
    plt.figure(figsize=(12, 8))
    for key, dist in smoothed_distances_to_combined.items():
        plt.plot(dist, label=f"{key} to Combined Center")
    plt.title("Distances to Combined Center of Mass Over Time")
    plt.xlabel("Frame")
    plt.ylabel("Distance (Å)")
    plt.legend()
    plt.grid(True)
    plt.savefig(output_file.replace(".png", "_to_combined_center.png"), dpi=300)
    plt.show()

    print(f"Plots saved: pair distances and distances to combined center.")

# Example usage
if __name__ == "__main__":
    # Define input files
    psf_file = "wr_6Ca.psf"
    dcd_file = "wr_500ns_Cx46_6Ca_rep1.dcd"
    output_file = "distances_analysis.png"

    # Calculate and plot distances with smoothing
    calculate_and_plot_distances(psf_file, dcd_file, output_file, smoothing_sigma=200)

