set dcdname [lindex $argv 0]

# Load trajectory
mol new pca_${dcdname}_cMD100ns_CA.psf type psf
mol addfile pca_${dcdname}_cMD100ns_CA.dcd type dcd first 0 last -1 step 1 waitfor all

# Salt bridges
package require saltbr
mkdir saltbridges_${dcdname}
saltbr -sel [atomselect top protein] -log saltbridges_${dcdname}.log -outdir ./saltbridges_${dcdname}

# Rg calculation
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
set outfile [open rog_${dcdname}.dat w]
set nf [molinfo top get numframes] 
set i 0
while {$i < $nf} {

    set prot [atomselect top "protein" frame $i]

    set i [expr {$i + 1}]
    set rog [gyr_radius $prot]

    puts $outfile "$i	$rog"

} 
close $outfile

# Align to first frame
set reference [atomselect top "backbone" frame 0]
		# the frame being compared
set compare [atomselect top "backbone"]
set num_steps [molinfo top get numframes]
set num [expr {$num_steps - 1}]
for {set frame 0} {$frame < $num} {incr frame} {
		# get the correct frame
		$compare frame $frame
		# compute the transformation
		set trans_mat [measure fit $compare $reference]
		# do the alignment
		$compare move $trans_mat
}
# Backbone RMSD
set outfile [open rmsd_${dcdname}.dat w]
for {set frame 0} {$frame < $num_steps} {incr frame} {
		# get the correct frame
   $compare frame $frame
		# compute the transformation
   set trans_mat [measure fit $compare $reference]
		# do the alignment
   $compare move $trans_mat
		# compute the RMSD
   set rmsd [measure rmsd $compare $reference]
		# print the RMSD
   puts $outfile "$frame $rmsd"
}
close $outfile

# RMSF calculation
set outfile [open rmsf_${dcdname}.dat w]
set sel [atomselect top "name CA"]
set rmsf [measure rmsf $sel first 0 last $num step 1]
for {set i 0} {$i < [$sel num]} {incr i} {
  puts $outfile "[expr {$i+1}] [lindex $rmsf $i]"
} 
close $outfile

# ddG
mol delete all
mol new pca_${dcdname}_cMD100ns_CA.psf type psf
mol addfile pca_${dcdname}_cMD100ns_CA.dcd type dcd first 0 last -1 step 1 waitfor all

mkdir ddG_${dcdname}
cd ddG_${dcdname}
set frameNumber [molinfo 0 get frame]
for {set i 0} {$i <= $frameNumber} {incr i 100} {
	[atomselect top "resname HSD or resname HSE"] set resname HIS
    [atomselect top all frame $i] writepdb $i.pdb
	foldx --command=Stability --pdb=$i.pdb
}
set i [expr {$i - 1}]
[atomselect top all frame $i] writepdb $i.pdb
foldx --command=Stability --pdb=$i.pdb

set outfile [open stability_${dcdname}.dat w]
set i 0
while {$i < $frameNumber} {
	puts $outfile "$i	[exec gawk {FNR==1 {print $2}} ${i}_0_ST.fxout]"
	set i [expr {$i + 100}]
}
set i [expr {$i - 1}]
puts $outfile "$i	[exec gawk {FNR==1 {print $2}} ${i}_0_ST.fxout]"
close $outfile
#mv stability_${dcdname}.dat ../.
file rename stability_${dcdname}.dat ../.
cd ..
#rm -r ./ddG_${dcdname}
file delete -force ./ddG_${dcdname}


exit

