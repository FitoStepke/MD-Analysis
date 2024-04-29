pbc unwrap -sel "protein or (resname CAL and not protein)" -all
pbc wrap -centersel "protein and backbone" -center com -compound residue -all
set nf [molinfo top get numframes]
for {set i 0} {$i <= $nf} {incr i} {
set sel [atomselect top "protein and backbone" frame $i]
set all [atomselect top all frame $i] 
$all moveby [vecinvert [measure center $sel]]}


