set dcdname [lindex $argv 0]

# Load trajectory
mol new pca_${dcdname}_cMD100ns_CA.psf type psf
mol addfile pca_${dcdname}_cMD100ns_CA.dcd type dcd first 0 last -1 step 1 waitfor all

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

set reference [atomselect top "backbone and (resid 607 to 623)" frame 0]
		# the frame being compared
set compare [atomselect top "backbone and (resid 607 to 623)"]
# Backbone RMSD
set outfile [open rmsd_ligandBinding_${dcdname}.dat w]
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

exit

