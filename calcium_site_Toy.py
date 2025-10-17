import MDAnalysis as mda
import numpy as np

def analyze_trajectory(psf_file, dcd_file, output_file, max_distance):
    """
    Analyze the trajectory to determine which protein group (site)
    is closest to the single calcium ion in each frame.
    Marks as "OUT" if the calcium is farther than the defined threshold.
    """
    # Load the trajectory
    u = mda.Universe(psf_file, dcd_file)

    # Calcium selection for a single ion
    cal_selection = "resname CAL"
    calcium_group = u.select_atoms(cal_selection)

    # New protein groups (ordered clockwise)
    protein_selections = {
      "AB": "(protein and (segid S5 or segid S6)) or (protein and (segid S4 or segid S11))",
        "BC": "(protein and (segid S4 or segid S11)) or (protein and (segid S3 or segid S9))",
        "CD": "(protein and (segid S3 or segid S9)) or (protein and (segid S1 or segid S10))",
        "DE": "(protein and (segid S1 or segid S10)) or (protein and (segid S8 or segid S12))",
        "EF": "(protein and (segid S8 or segid S12)) or (protein and (segid S2 or segid S13))",
        "FA": "(protein and (segid S2 or segid S13)) or (protein and (segid S5 or segid S6))",
    }
    protein_groups = {key: u.select_atoms(sel) for key, sel in protein_selections.items()}

    # Open output file for writing results
    with open(output_file, "w") as f_out:
        # Write header
        f_out.write("Frame\tGroup\n")
        # Iterate over frames in the trajectory
        for ts in u.trajectory:
            frame = ts.frame
            ca_com = calcium_group.center_of_mass()
            min_distance = np.inf
            closest_group = "OUT"
            # Calculate distance to each protein group
            for group_name, protein_group in protein_groups.items():
                distance = np.linalg.norm(ca_com - protein_group.center_of_mass())
                if distance < min_distance:
                    min_distance = distance
                    closest_group = group_name
            # Only assign if within threshold, else OUT
            if min_distance > max_distance:
                closest_group = "OUT"
            f_out.write(f"{frame}\t{closest_group}\n")
    print(f"Analysis complete. Results saved to {output_file}")

# Example usage
if __name__ == "__main__":
    psf_file = "wr_6Ca.psf"
    dcd_file = "wr_500ns_Cx46_6Ca_rep1.dcd"
    output_file = "calcium_site_simple.txt"
    max_distance = 5  # Change as needed
    analyze_trajectory(psf_file, dcd_file, output_file, max_distance)

