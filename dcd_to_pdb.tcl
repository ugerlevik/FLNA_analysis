set dcdname [lindex $argv 0]

mol new pca_${dcdname}_cMD100ns_CA.psf type psf
mol addfile pca_${dcdname}_cMD100ns_CA.dcd type dcd first 0 last -1 step 5 waitfor all

animate write pdb 4M9P_${dcdname}_100ns_repeat1_SezermanLab.pdb beg 0 end -1 skip 1 sel [atomselect top protein]

exit