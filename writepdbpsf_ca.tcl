set dcdname [lindex $argv 0]

mol new ./$dcdname/2_min_equ/${dcdname}_ionized.psf type psf
mol addfile ./$dcdname/2_min_equ/${dcdname}_ionized.pdb type pdb

set sel [atomselect 0 "protein"]

$sel writepdb pca_${dcdname}_cMD100ns_CA.pdb
$sel writepsf pca_${dcdname}_cMD100ns_CA.psf

exit