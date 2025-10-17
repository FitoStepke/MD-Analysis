import os

def ratio(output_file, timestep_fs, max_out_frames=100, min_valid_frames=100):
    """
    Analyze a single Ca2+ position assignment file to compute the ratio (events per microsecond)
    of site transitions, with tolerance for short OUT periods (max_out_frames).

    Parameters:
        timestep_fs (float): Trajectory timestep in femtoseconds (fs). Use 1000 if each frame is 1 ps.
        max_out_frames (int): Maximum consecutive OUT frames tolerated before site is considered abandoned.
        min_valid_frames (int): Minimum frames spent in a site for a transition to be counted as an event.
    """
    # Read analysis file
    with open(output_file, "r") as f:
        lines = f.readlines()
    group_col = 1  # Group assignment column

    # Assignment and event counters
    last_site = None
    site_start = 0
    site_time = 0
    out_count = 0
    event_count = 0
    valid_frames = 0

    # Trackers
    site_history = []

    # Iterate frames, skip header
    for i, line in enumerate(lines[1:]):
        values = line.strip().split("\t")
        group = values[group_col]

        # Start site if first
        if last_site is None and group != "OUT":
            last_site = group
            site_start = i
            site_time = 1
            out_count = 0
            continue

        if group == last_site:
            # Continues in the same site
            if group != "OUT":
                site_time += 1
                out_count = 0
                valid_frames += 1
            else:
                # Out but within tolerance: count towards site_time
                out_count += 1
                if out_count <= max_out_frames:
                    site_time += 1
                else:
                    # Too many OUT frames: site considered abandoned
                    if site_time >= min_valid_frames:
                        site_history.append((last_site, site_time))
                    last_site = "OUT"
                    site_start = i
                    site_time = 1
                    out_count = 1
        else:
            # Switching sites or OUT after a different site
            if group != "OUT":
                # Site change detected
                if site_time >= min_valid_frames and last_site != "OUT":
                    event_count += 1
                    site_history.append((last_site, site_time))
                last_site = group
                site_start = i
                site_time = 1
                out_count = 0
                valid_frames += 1
            else:
                # OUT after leaving site; start counting OUT sequence
                out_count = 1
                last_site = "OUT"
                site_start = i
                site_time = 1

    # Add last site if ended in valid position
    if last_site != "OUT" and site_time >= min_valid_frames:
        site_history.append((last_site, site_time))

    # Calculate simulation time in microseconds
    time_fs = valid_frames * timestep_fs
    time_us = time_fs / 1e6

    # Ratio (events per microsecond)
    ratio_value = event_count / time_us if time_us > 0 else 0

    # Write to output
    output_analysis_file = os.path.splitext(output_file)[0] + "_ratio.txt"
    with open(output_analysis_file, "w") as f_out:
        f_out.write("Ca2+ site transition event analysis (single ion):\n")
        f_out.write(f"Event count: {event_count}\n")
        f_out.write(f"Valid frames (not OUT): {valid_frames}\n")
        f_out.write(f"Ratio (events/µs): {ratio_value:.3f}\n")
        f_out.write("\nSite durations:\n")
        for site, duration in site_history:
            f_out.write(f"{site}: {duration} frames\n")

    print(f"Analysis complete. Ratio saved to {output_analysis_file}")

# Example usage:
if __name__ == "__main__":
    output_file = "calcium_site.txt"
    timestep_fs = 1000   # Each frame is 1 ps
    max_out_frames = 5000000 # Set to your tolerance for an abandoned site
    min_valid_frames = 500000 # Minimum frames for site occupancy to count transitions
    ratio(output_file, timestep_fs, max_out_frames, min_valid_frames)






















    
   
