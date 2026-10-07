@echo off
cd /d "%~dp0"
set "PATH=%~dp0..\.conda;%~dp0..\.conda\Scripts;%~dp0..\.conda\Library\bin;%PATH%"

echo [1/2] Running DEVSIM C-V simulation...
"%~dp0..\.conda\python.exe" mosfet.py cv
if errorlevel 1 (
    echo.
    echo Simulation failed. Check the error message above.
    pause
    exit /b 1
)

echo.
echo [2/2] Opening C-V graph...
"%~dp0..\.conda\python.exe" plot_cv.py

pause
