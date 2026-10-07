@echo off
cd /d "%~dp0"
set "PATH=%~dp0..\.conda;%~dp0..\.conda\Scripts;%~dp0..\.conda\Library\bin;%PATH%"

echo Running Step 1 vs Step 2 Id-Vg comparison...
"%~dp0..\.conda\python.exe" compare_models.py

if errorlevel 1 (
    echo.
    echo Comparison failed. Check the error message above.
    pause
    exit /b 1
)

pause
