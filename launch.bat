@echo off
title Motif Scanner
cd /d "%~dp0"

echo.
echo   Starting Motif Scanner...
echo.

start "" "http://127.0.0.1:5000"
timeout /t 2 /nobreak >nul

python app.py
pause