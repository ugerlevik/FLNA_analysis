@echo off
for /F "tokens=1" %%a in (dcdnames.txt) do (
echo %%a
prody catdcd .\%%a\3_production_repeat1\%%a_ionized_production_repeat1_stride2.dcd --pdb .\%%a\2_min_equ\%%a_ionized.pdb --psf .\%%a\2_min_equ\%%a_ionized.psf -o pca_%%a_cMD100ns_CA.dcd -s "protein" --stride 2
)