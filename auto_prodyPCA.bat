@echo off
for /F "tokens=1" %%a in (dcdnames.txt) do (
echo %%a
pca.py 4m9p %%a
)