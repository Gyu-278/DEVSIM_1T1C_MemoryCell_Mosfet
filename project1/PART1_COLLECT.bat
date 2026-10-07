@echo off
chcp 65001 >nul
cd /d "%~dp0"
set "PATH=%~dp0..\.conda;%~dp0..\.conda\Scripts;%~dp0..\.conda\Library\bin;%PATH%"

echo ============================================================
echo  [6] 탐색 기록표 모으기 (보고서 2번 항목)
echo ============================================================
echo  결과: results\search_record.csv  ·  이미 잰 결과만 모으므로 즉시 끝납니다
echo.
"%~dp0..\.conda\python.exe" -u -m part1.collect
if errorlevel 1 (
  echo.
  echo [FAILED] 위 오류 메시지를 확인하세요.
)
echo.
pause
