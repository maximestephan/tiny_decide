@echo off
setlocal

cd /d "%USERPROFILE%\tiny_decide"
git pull
"%USERPROFILE%\miniconda3\Scripts\activate.bat" tinydecide
python -m flask --app main.py run  --port 8080
start "" "http://127.0.0.1:8080"