@echo off
for /F "tokens=1" %%a in (dcdnames.txt) do (
echo %%a
vmd -dispdev text -e writepdb_each10ns.tcl -args %%a
)