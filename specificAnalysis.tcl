set dcdname [lindex $argv 0]

# Load trajectory
mol new pca_${dcdname}_cMD100ns_CA.psf type psf
mol addfile pca_${dcdname}_cMD100ns_CA.dcd type dcd first 0 last -1 step 1 waitfor all

# Distance definition
proc distanceM {seltext1 seltext2 f_r_out} {
	set sel1 [atomselect top "$seltext1"]
	set sel2 [atomselect top "$seltext2"] 
	set nf [molinfo top get numframes]
	set outfile [open $f_r_out w]
	for {set i 0} {$i < $nf} {incr i} {
		puts "frame $i of $nf"
		$sel1 frame $i
		$sel2 frame $i
		set com1 [measure center $sel1 weight mass]
		set com2 [measure center $sel2 weight mass]
		set simdata($i.r) [veclength [vecsub $com1 $com2]]
		puts $outfile "$i,$simdata($i.r)"
	}
	close $outfile
}

# Distance btw IgFLNa3-4
distanceM "resid 478 to 575" "resid 576 to 666" dist_IgFLNa3_4_${dcdname}.dat

# Distance btw IgFLNa3-5
distanceM "resid 478 to 575" "resid 667 to 766" dist_IgFLNa3_5_${dcdname}.dat

# Distance btw IgFLNa4-5
distanceM "resid 576 to 666" "resid 667 to 766" dist_IgFLNa4_5_${dcdname}.dat

# Distance btw IgFLNa3-4 (only interaction surface residues)
distanceM "resid 482 to 507" "resid 640 to 662" dist_IgFLNa3_4_onlyInteractionSurfaceResidues_${dcdname}.dat

# SASA overall
set sel [atomselect top "protein"]
set protein [atomselect top "protein"]
set n [molinfo top get numframes]
set output [open SASA_overall_${dcdname}.dat w]
# sasa calculation loop
for {set i 0} {$i < $n} {incr i} {
	molinfo top set frame $i
	set sasa [measure sasa 1.4 $protein -restrict $sel]
	puts "\t \t progress: $i/$n"
	puts $output "$sasa"
}
close $output

# SASA btw IgFLNa3 and IgFLNa4
set sel [atomselect top "resid 482 to 507 or resid 640 to 662"]
set protein [atomselect top "protein"]
set n [molinfo top get numframes]
set output [open SASA_btwIgFLNa3_4_${dcdname}.dat w]
# sasa calculation loop
for {set i 0} {$i < $n} {incr i} {
	molinfo top set frame $i
	set sasa [measure sasa 1.4 $protein -restrict $sel]
	puts "\t \t progress: $i/$n"
	puts $output "$sasa"
}
close $output

# Rg Definitions
proc gyr_radius {sel} {
  # make sure this is a proper selection and has atoms
  if {[$sel num] <= 0} {
    error "gyr_radius: must have at least one atom in selection"
  }
  # gyration is sqrt( sum((r(i) - r(center_of_mass))^2) / N)
  set com [center_of_mass $sel]
  set sum 0
  foreach coord [$sel get {x y z}] {
    set sum [vecadd $sum [veclength2 [vecsub $coord $com]]]
  }
  return [expr sqrt($sum / ([$sel num] + 0.0))]
}
proc center_of_mass {selection} {
        # some error checking
        if {[$selection num] <= 0} {
                error "center_of_mass: needs a selection with atoms"
        }
        # set the center of mass to 0
        set com [veczero]
        # set the total mass to 0
        set mass 0
        # [$selection get {x y z}] returns the coordinates {x y z} 
        # [$selection get {mass}] returns the masses
        # so the following says "for each pair of {coordinates} and masses,
	#  do the computation ..."
        foreach coord [$selection get {x y z}] m [$selection get mass] {
           # sum of the masses
           set mass [expr $mass + $m]
           # sum up the product of mass and coordinate
           set com [vecadd $com [vecscale $m $coord]]
        }
        # and scale by the inverse of the number of atoms
        if {$mass == 0} {
                error "center_of_mass: total mass is zero"
        }
        # The "1.0" can't be "1", since otherwise integer division is done
        return [vecscale [expr 1.0/$mass] $com]
}

# Rg for IgFLNa3
set outfile [open Rg_IgFLNa3_${dcdname}.dat w]
set nf [molinfo top get numframes] 
set i 0
while {$i < $nf} {
    set prot [atomselect top "resid 478 to 575" frame $i]
    set i [expr {$i + 1}]
    set rog [gyr_radius $prot]
    puts $outfile "$i	$rog"
} 
close $outfile

# Rg for IgFLNa4
set outfile [open Rg_IgFLNa4_${dcdname}.dat w]
set nf [molinfo top get numframes] 
set i 0
while {$i < $nf} {
    set prot [atomselect top "resid 576 to 666" frame $i]
    set i [expr {$i + 1}]
    set rog [gyr_radius $prot]
    puts $outfile "$i	$rog"
} 
close $outfile

# Rg for IgFLNa5
set outfile [open Rg_IgFLNa5_${dcdname}.dat w]
set nf [molinfo top get numframes] 
set i 0
while {$i < $nf} {
    set prot [atomselect top "resid 667 to 766" frame $i]
    set i [expr {$i + 1}]
    set rog [gyr_radius $prot]
    puts $outfile "$i	$rog"
} 
close $outfile

# Number of H-bonds btw domains 3-4
package require hbonds
hbonds -sel1 [atomselect top "resid 478 to 575"] -writefile yes -plot no -outfile hbondsBtwIgFLNa3_4_${dcdname}.dat

exit

