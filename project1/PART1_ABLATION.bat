@echo off
chcp 65001 >nul
cd /d "%~dp0"
set "PATH=%~dp0..\.conda;%~dp0..\.conda\Scripts;%~dp0..\.conda\Library\bin;%PATH%"

echo ============================================================
echo  [4] 파라미터별 ablation (보고서 3번 항목)
echo ============================================================
echo  결과: results\ablation_final.csv  ·  18개 조합, 1시간 정도 걸립니다
echo.
"%~dp0..\.conda\python.exe" -u -m part1.ablation
if errorlevel 1 (
  echo.
  echo [FAILED] 위 오류 메시지를 확인하세요.
)
echo.
pause
