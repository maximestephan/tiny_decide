setlocal

REM Target directory
set TARGET_DIR=%USERPROFILE%\tiny_decide

REM Git repository URL (replace this)
set REPO_URL=https://github.com/maximestephan/tiny_decide.git

REM Check if target directory exists
if exist "%TARGET_DIR%" (
    REM Check if it's a git repository
    if exist "%TARGET_DIR%\.git" (
        echo Git repository found. Pulling latest changes...
        cd /d "%TARGET_DIR%"
        git pull
    ) else (
        echo Directory exists but is not a git repository.
        echo Aborting to avoid overwriting existing files.
    )
) else (
    echo Directory not found. Cloning repository...
    git clone "%REPO_URL%" "%TARGET_DIR%"
)


cd /d "%USERPROFILE%\tiny_decide"
call "%USERPROFILE%\miniconda3\Scripts\activate.bat" tinydecide
start "" /b python -m flask --app main.py run --port 8080
timeout /t 3 /nobreak > nul
start "" "http://127.0.0.1:8080"