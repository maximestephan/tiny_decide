@echo off
setlocal

set "MINICONDA_DIR=%USERPROFILE%\miniconda3"
set "APP_DIR=%USERPROFILE%\tiny_decide"

winget install --id Git.Git -e --source winget
git clone https://github.com/maximestephan/tiny_decide.git "%APP_DIR%"

winget install --id Anaconda.Miniconda3 -e --source winget
"%MINICONDA_DIR%\Scripts\conda.exe" tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
"%MINICONDA_DIR%\Scripts\conda.exe" tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
"%MINICONDA_DIR%\Scripts\conda.exe" tos accept --override-channels --channel https://repo.anaconda.com/pkgs/msys2

conda env list | findstr /b /c:"tinydecide" >nul 2>&1
if errorlevel 1 (
    echo Environment "tinydecide" not found. Creating it now...
    conda create -y -n tinydecide python=3.11 rdkit flask pandas requests -c conda-forge
) else (
    echo Environment tinydecide already exists.
)

