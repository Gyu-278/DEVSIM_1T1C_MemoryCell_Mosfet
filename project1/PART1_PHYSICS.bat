@echo off
chcp 65001 >nul
cd /d "%~dp0"
set "PATH=%~dp0..\.conda;%~dp0..\.conda\Scripts;%~dp0..\.conda\Library\bin;%PATH%"

echo ============================================================
echo  [5] 추가한 물리 전후 비교 (보고서 1번 항목)
echo ============================================================
echo  결과: results\physics_check_baseline.csv  ·  10분 정도
echo.
"%~dp0..\.conda\python.exe" -u -m part1.physics_check
if errorlevel 1 (
  echo.
  echo [FAILED] 위 오류 메시지를 확인하세요.
)
echo.
pause
