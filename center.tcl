# Define input and output file names
set psf_file "wr_3Ca.psf"       		;# Input PSF file name
set dcd_file "wr_500ns_Cx46_3Ca_rep0.dcd"       ;# Input DCD file name
set output_dcd "centered.dcd"  			;# Output DCD file name for the centered trajectory

# Load the PSF file
mol new $psf_file type psf

# Load the DCD file
mol addfile $dcd_file type dcd waitfor all

# Place all subunits in the same periodic cell
package require pbctools
pbc unwrap -sel "protein or (resname CAL and not protein)" -all
pbc wrap -centersel "protein and backbone" -center com -compound residue -all

# Get the total number of frames in the trajectory
set nf [molinfo top get numframes]

# Process each frame
for {set i 0} {$i < $nf} {incr i} {
    # Select the protein backbone for centering
    set sel [atomselect top "protein and backbone" frame $i]
    
    # Select all atoms for the translation
    set all [atomselect top all frame $i]

    # Calculate the center of mass (COM) of the protein backbone
    set com [measure center $sel weight mass]

    # Move all atoms so that the protein backbone COM is at the origin
    $all moveby [vecinvert $com]
}

# Save the processed trajectory
animate write dcd $output_dcd

# Clean up
mol delete all
exit

