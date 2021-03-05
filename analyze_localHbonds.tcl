#set dcdname [lindex $argv 0]
set wtname "4m9p"
set mutname "L656F"

set case "resid 656"

# Load trajectory
mol new pca_${wtname}_cMD100ns_CA.psf type psf
mol addfile pca_${wtname}_cMD100ns_CA.dcd type dcd first 0 last -1 step 1 waitfor all

# Number of H-bonds
package require hbonds
hbonds -sel1 [atomselect top "same residue as (within 10 of ${case})"] -writefile yes -plot no -outfile hbonds10of${case}_wt.dat

# Load trajectory
mol new pca_${mutname}_cMD100ns_CA.psf type psf
mol addfile pca_${mutname}_cMD100ns_CA.dcd type dcd first 0 last -1 step 1 waitfor all

# Number of H-bonds
package require hbonds
hbonds -sel1 [atomselect top "same residue as (within 10 of ${case})"] -writefile yes -plot no -outfile hbonds10of${case}_${mutname}.dat

exit

