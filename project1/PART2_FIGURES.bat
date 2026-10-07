@echo off
chcp 65001 >nul
cd /d "%~dp0"
set "PATH=%~dp0..\.conda;%~dp0..\.conda\Scripts;%~dp0..\.conda\Library\bin;%PATH%"

echo ============================================================
echo  [P2-3] 셀 구조 그림 · ΔV_BL(t) 과도 곡선
echo ============================================================
echo  결과: results\fig_p2_*.html, p2_read_curves_G2.csv
echo.
"%~dp0..\.conda\python.exe" -u -m part2.figures_p2
if errorlevel 1 (
  echo.
  echo [FAILED] 위 오류 메시지를 확인하세요.
)
echo.
pause
