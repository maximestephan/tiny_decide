@echo off
setlocal


winget install --id Git.Git -e --source winget
REM refresh path
for /f "tokens=2*" %%a in ('reg query HKLM\System\CurrentControlSet\Control\Session\Manager\Environment /v Path') do set "PATH=%%b;%PATH%"
git clone https://github.com/maximestephan/tiny_decide.git "%USERPROFILE%\tiny_decide"

winget install --id Anaconda.Miniconda3 -e --source winget
"%USERPROFILE%\miniconda3\Scripts\conda.exe"  tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
"%USERPROFILE%\miniconda3\Scripts\conda.exe" tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
"%USERPROFILE%\miniconda3\Scripts\conda.exe" tos accept --override-channels --channel https://repo.anaconda.com/pkgs/msys2
"%USERPROFILE%\miniconda3\Scripts\conda.exe" env list | findstr /b /c:"tinydecide" >nul 2>&1
if errorlevel 1 (
    echo Environment "tinydecide" not found. Creating it now...
    "%USERPROFILE%\miniconda3\Scripts\conda.exe" create -y -n tinydecide python=3.11 rdkit flask pandas requests -c conda-forge
) else (
    echo Environment tinydecide already exists.
)
pause

start "" "%USERPROFILE%\tiny_decide\run.bat"
exit