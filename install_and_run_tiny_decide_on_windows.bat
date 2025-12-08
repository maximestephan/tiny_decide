@echo off
setlocal

echo =========================================
echo tiny_decide setup and launcher
echo =========================================
echo.

REM -----------------------------------------
REM 0) Compute user-specific Miniconda path
REM -----------------------------------------
set "MINICONDA_DIR=%USERPROFILE%\miniconda3"

REM -----------------------------------------
REM 1) Check Git installation
REM -----------------------------------------
echo [1/5] Checking Git...

git --version >nul 2>&1
if %errorlevel% neq 0 (
    echo Git not found. Installing Git using winget...
    winget install --id Git.Git -e --source winget


    echo.
    echo Git installation completed.
    echo PLEASE CLOSE THIS WINDOW COMPLETELY AND RUN THE SCRIPT AGAIN.
    echo (This is necessary for PATH updates.)
    echo.
    pause
    goto :end
) else (
    echo Git is installed:
    git --version
    echo.
)

REM -----------------------------------------
REM 2) Check Miniconda installation (folder-based detection)
REM -----------------------------------------
echo [2/5] Checking Miniconda installation...

if not exist "%MINICONDA_DIR%\Scripts\activate.bat" (
    echo Miniconda not found at: %MINICONDA_DIR%
    echo Installing Miniconda via winget...

    winget install --id Anaconda.Miniconda3 -e --source winget

"%MINICONDA_DIR%\Scripts\conda.exe" tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
"%MINICONDA_DIR%\Scripts\conda.exe" tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
"%MINICONDA_DIR%\Scripts\conda.exe" tos accept --override-channels --channel https://repo.anaconda.com/pkgs/msys2

    echo.
    echo Miniconda installation completed.
    echo PLEASE CLOSE THIS WINDOW COMPLETELY AND RUN THE SCRIPT AGAIN.
    echo (This is necessary for PATH and base activation.)
    echo.
    pause
    goto :end

)

echo Miniconda found at: %MINICONDA_DIR%
echo.

REM -----------------------------------------
REM 3) Activate base conda environment
REM -----------------------------------------
echo [3/5] Activating base conda environment...
call "%MINICONDA_DIR%\Scripts\activate.bat"
if errorlevel 1 (
    echo ERROR: Failed to activate base conda environment.
    pause

)
echo Base conda activated.
echo.

REM -----------------------------------------
REM 4) Create / update env tinydecide
REM -----------------------------------------
echo [4/5] Checking conda environment tinydecide...

conda env list | findstr /b /c:"tinydecide" >nul 2>&1
if errorlevel 1 (
    echo Environment "tinydecide" not found. Creating it now...

    conda create -y -n tinydecide python=3.11 rdkit flask pandas requests -c conda-forge
) else (
    echo Environment tinydecide already exists.
)

echo Activating environment tinydecide...
call "%MINICONDA_DIR%\Scripts\activate.bat" tinydecide
if errorlevel 1 (
    echo ERROR: Failed to activate tinydecide.
    pause
)



REM -----------------------------------------
REM 5) Clone or update repository
REM -----------------------------------------
set "APP_DIR=%USERPROFILE%\tiny_decide"

echo [5/5] Preparing application folder: %APP_DIR%

if not exist "%APP_DIR%" (
    echo Cloning repository for the first time...
    git clone https://github.com/maximestephan/tiny_decide.git "%APP_DIR%"
) else (
    if not exist "%APP_DIR%\.git" (
        echo ERROR: Folder exists but is NOT a git repository.
        echo Please delete or rename "%APP_DIR%" and retry.
        pause
        goto :end

    )
    echo Repository exists. Pulling latest updates...
    cd /d "%APP_DIR%"
    git pull
)

cd /d "%APP_DIR%"
echo App directory ready.
echo.

REM -----------------------------------------
REM Start Flask and open browser
REM -----------------------------------------
echo Launching Flask server on http://127.0.0.1:8080 ...

start "tiny_decide Flask" cmd /k "call %MINICONDA_DIR%\Scripts\activate.bat tinydecide && cd /d %APP_DIR% && python -m flask --app main.py run --debug --port 8080"

echo Waiting before opening browser...
timeout /t 15 /nobreak >nul

start "" "http://127.0.0.1:8080"

echo.
echo INSTALLATION COMPLETE - Flask is running in a separate window.
echo.


:END
pause
endlocal
exit /b 0
