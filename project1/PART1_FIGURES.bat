@echo off
chcp 65001 >nul
cd /d "%~dp0"
set "PATH=%~dp0..\.conda;%~dp0..\.conda\Scripts;%~dp0..\.conda\Library\bin;%PATH%"

echo ============================================================
echo  [3] 구조 그림 / Id-Vg 특성 곡선
echo ============================================================
echo  결과: results\fig_*.html (브라우저로 열림) 과 results\curves_A.csv
echo.
"%~dp0..\.conda\python.exe" -u -m part1.figures
if errorlevel 1 (
  echo.
  echo [FAILED] 위 오류 메시지를 확인하세요.
)
echo.
pause
