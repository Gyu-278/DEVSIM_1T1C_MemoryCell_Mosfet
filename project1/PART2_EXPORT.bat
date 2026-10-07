@echo off
chcp 65001 >nul
cd /d "%~dp0"
set "PATH=%~dp0..\.conda;%~dp0..\.conda\Scripts;%~dp0..\.conda\Library\bin;%PATH%"

echo ============================================================
echo  [P2-1] 제출 파일 생성 + 공식 검사기
echo ============================================================
echo  결과: part2_2021144048.devsim  ·  OK 가 뜨면 그 화면을 캡처하세요
echo.
"%~dp0..\.conda\python.exe" -u -m part2.export_devsim2
if errorlevel 1 (
  echo.
  echo [FAILED] 위 오류 메시지를 확인하세요.
)
echo.
pause
