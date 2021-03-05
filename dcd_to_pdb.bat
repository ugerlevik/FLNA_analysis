@echo off
for /F "tokens=1" %%a in (dcdnames.txt) do (
echo %%a
vmd -dispdev text -e dcd_to_pdb.tcl -args %%a
)