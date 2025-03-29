import MDAnalysis as mda
import numpy as np

def analyze_trajectory(psf_file, dcd_file, num_calcium, output_file, max_distance):
    """
    Analyze a trajectory to determine which protein selection is closest
    to calcium ions in each frame. Marks calcium as "OUT" if it is beyond
    a defined distance from all selections.
    
    Parameters:
        psf_file (str): Path to the PSF file.
        dcd_file (str): Path to the DCD file.
        num_calcium (int): Number of calcium ions to analyze (3 or 6).
        output_file (str): Path to the output file.
        max_distance (float): Maximum distance to consider a calcium "near" a selection.
    """
    # Load the trajectory
    u = mda.Universe(psf_file, dcd_file)

    # Define calcium selections
    cal_selections = [
        "resname CAL and resid 341",
        "resname CAL and resid 342",
        "resname CAL and resid 343",
        "resname CAL and resid 344",
        "resname CAL and resid 345",
        "resname CAL and resid 346"
    ]

    # Limit to the specified number of calcium ions
    cal_selections = cal_selections[:num_calcium]

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
    calcium_groups = [u.select_atoms(sel) for sel in cal_selections]

    # Open the output file
    with open(output_file, "w") as f_out:
        # Write header
        header = "Frame\t" + "\t".join([f"CAL{i+1}" for i in range(num_calcium)]) + "\n"
        f_out.write(header)

        # Iterate over each frame in the trajectory
        for ts in u.trajectory:
            results = [ts.frame]

            # Calculate distances for each calcium ion
            for cal_group in calcium_groups:
                min_distance = np.inf
                closest_selection = "OUT"  # Default to "OUT" if no selection is close

                # Find the closest protein selection
                for key, protein_group in protein_groups.items():
                    distance = np.linalg.norm(cal_group.center_of_mass() - protein_group.center_of_mass())
                    if distance < min_distance:
                        min_distance = distance
                        closest_selection = key

                # Check if the minimum distance is within the defined threshold
                if min_distance > max_distance:
                    closest_selection = "OUT"

                # Append the closest selection or "OUT" to the results
                results.append(closest_selection)

            # Write results to file
            f_out.write("\t".join(map(str, results)) + "\n")

    print(f"Analysis complete. Results saved to {output_file}")

# Example usage
if __name__ == "__main__":
    # Define input parameters
    psf_file = "wr_6Ca.psf"
    dcd_file = "wr_500ns_Cx46_6Ca_rep1.dcd"
    output_file = "calcium_sites_6R1.txt"
    
    # Prompt the user to specify the number of calcium ions and max distance
    num_calcium = 6
    max_distance = 5
    
    # Run the analysis
    analyze_trajectory(psf_file, dcd_file, num_calcium, output_file, max_distance)
