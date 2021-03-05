set dcdname [lindex $argv 0]

mol new pca_${dcdname}_cMD100ns_CA.psf type psf
mol addfile pca_${dcdname}_cMD100ns_CA.dcd type dcd first 0 last -1 step 50 waitfor all

set frameNumber [molinfo 0 get frame]
set i 0
mkdir forStride

set sel [atomselect 0 "protein" frame $i]
$sel writepdb ./forStride/${dcdname}_prot_${i}ns.pdb

set i 9
while {$i <= $frameNumber} {
	set sel [atomselect 0 "protein" frame $i]
	$sel writepdb ./forStride/${dcdname}_prot_[expr {$i + 1}]ns.pdb
	set i [expr {$i + 10}]
}

exit